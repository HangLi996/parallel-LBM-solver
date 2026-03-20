#include "Simulation.h"
#include "node.h"
#include <iostream>
#include <sstream>
#include <mpi.h>
#include <vector>
#include <map>
#include <cstdlib>      // free()
#include <unistd.h>     // usleep

#include "../LBM/parallel.h"
#include <omp.h>


namespace LBM {

    Simulation::Simulation(Initializer_Ptr initializer)
    :
        d_domain(initializer->domain()),
        d_stream_time(0.0),
        d_communicate_time(0.0),
        d_collision_time(0.0),
        d_poststream_time(0.0),
        d_step_count(0)
    {
        // Synchronize all processes before starting simulation
        MPI_Barrier(d_domain->comm);
    }

    Simulation::~Simulation()
    {
        // cleanup handled by vectors
    }

    void Simulation::step()
    {
        // Time the stream phase
        double t_start = MPI_Wtime();
        stream();
        d_stream_time += MPI_Wtime() - t_start;
        
        // Time the communication phase
        t_start = MPI_Wtime();
        communicate();
        d_communicate_time += MPI_Wtime() - t_start;
        
        // Time the post-stream processing phase
        t_start = MPI_Wtime();
        postStreamProcess();
        d_poststream_time += MPI_Wtime() - t_start;
        
        // Time the collision phase
        t_start = MPI_Wtime();
        collission();
        d_collision_time += MPI_Wtime() - t_start;
        
        d_step_count++;
    }

    void Simulation::stream()
    {
        size_t nDirections = d_domain->set->nDirections;
        size_t nNodes = d_domain->distribution_values.size() / nDirections;
        size_t null_neighbour_count = 0;

        #pragma omp parallel for reduction(+:null_neighbour_count) schedule(static)
        for (size_t idx = 0; idx < nNodes; ++idx)
        {
            size_t offset = idx * nDirections;
            
            for (size_t dir = 0; dir < nDirections; ++dir)
            {
                if (d_domain->distribution_neighbours[offset + dir] != nullptr)
                {
                    *d_domain->distribution_neighbours[offset + dir] = d_domain->distribution_values[offset + dir];
                }
                else
                {
                    null_neighbour_count++;
                }
            }
        }
    }

    void Simulation::collission()
    {
        double omega = d_domain->omega;
        size_t nDirections = d_domain->set->nDirections;
        size_t nNodes = d_domain->distribution_values.size() / nDirections;

        #pragma omp parallel for schedule(static)   
        for (size_t idx = 0; idx < nNodes; ++idx)
        {
            size_t offset = idx * nDirections;

            // switch to the newly streamed distribution values
            for (size_t dir = 0; dir < nDirections; ++dir)
            {
                if (d_domain->distribution_nextValues[offset + dir] >= 0.0) {
                    d_domain->distribution_values[offset + dir] = d_domain->distribution_nextValues[offset + dir];
                }
                // Reset nextValue for next iteration
                d_domain->distribution_nextValues[offset + dir] = -1.0;
            }

            // apply BGK approximation
            auto node_equilibrium = equilibrium(d_domain->set, d_domain->distribution_values, idx);
            for (size_t dir = 0; dir < nDirections; ++dir)
                d_domain->distribution_values[offset + dir] = d_domain->distribution_values[offset + dir] -
                    omega * (d_domain->distribution_values[offset + dir] - node_equilibrium[dir]);

            delete[] node_equilibrium;
        }
    }

    void Simulation::postStreamProcess()
    {
        for (size_t idx = 0; idx < d_domain->post_processors.size(); ++idx)
            d_domain->post_processors[idx]->process();
    }


    void Simulation::report(::Reporting::MatlabReporter reporter)
    {
        // reporter.reportOnTimeStep(d_domain->set, d_domain->nodes);
    }

    void Simulation::communicate()
    {
        int rank, size;
        MPI_Comm comm = d_domain->comm;
        MPI_Comm_rank(comm, &rank);
        MPI_Comm_size(comm, &size);

        std::vector<MPI_Request> send_requests;
        // Use d_domain->messengers directly
        auto& messengers = d_domain->messengers;
        send_requests.reserve(messengers.size() * 2);
        
        size_t nNodes = d_domain->distribution_values.size() / d_domain->set->nDirections;
        size_t nDimensions = d_domain->set->nDimensions;

        size_t skipped_sends = 0;
        size_t inferred_sends = 0;
        size_t total_sends = 0;
        
        for (size_t msg_idx = 0; msg_idx < messengers.size(); ++msg_idx)
        {
            auto& messenger = messengers[msg_idx];
            size_t node_idx = messenger.d_tag[0];
            size_t dir = messenger.d_tag[1];
            
            if (node_idx == SIZE_MAX || node_idx >= nNodes)
            {
                if (d_domain->messenger_to_sender.find(msg_idx) != d_domain->messenger_to_sender.end())
                {
                    auto sender_info = d_domain->messenger_to_sender[msg_idx];
                    size_t sender_node_idx = sender_info.first;
                    size_t sender_dir = sender_info.second;
                    
                    if (sender_node_idx < nNodes && 
                        sender_dir < d_domain->set->nDirections)
                    {
                        std::vector<int> sender_position;
                        size_t pos_offset = sender_node_idx * nDimensions;
                        // Use flat position vector
                        // Note: nDimensions check? position vector size is nNodes * nDimensions?
                        // Assuming position vector is populated.
                        if (pos_offset + nDimensions <= d_domain->position.size()) {
                            for (size_t dim = 0; dim < nDimensions; ++dim)
                                sender_position.push_back(static_cast<int>(d_domain->position[pos_offset + dim]));
                        }
                        
                        std::vector<int> receiver_position;
                        auto dir_vec = d_domain->set->direction(sender_dir);
                        
                        if (sender_position.size() <= d_domain->domain_size.size())
                        {
                            for (size_t dim = 0; dim < sender_position.size(); ++dim)
                            {
                                int domain_size_dim = static_cast<int>(d_domain->domain_size[dim]);
                                int coord = (sender_position[dim] + dir_vec[dim] + domain_size_dim) % domain_size_dim;
                                receiver_position.push_back(coord);
                            }
                            
                            size_t receiver_hash = receiver_position[0];
                            size_t multiplier = d_domain->domain_size[0];
                            for (size_t dim = 1; dim < receiver_position.size(); ++dim)
                            {
                                receiver_hash += receiver_position[dim] * multiplier;
                                multiplier *= d_domain->domain_size[dim];
                            }
                            
                            size_t tag_data[2] = {receiver_hash, dir};
                            MPI_Request req;
                            MPI_Isend(tag_data, 2, MPI_UNSIGNED_LONG, 
                                     static_cast<int>(messenger.d_p), 0, comm, &req);
                            send_requests.push_back(req);
                            
                            MPI_Isend(&messenger.d_src, 1, MPI_DOUBLE, 
                                     static_cast<int>(messenger.d_p), 1, comm, &req);
                            send_requests.push_back(req);
                            inferred_sends++;
                            continue;
                        }
                    }
                }
                
                skipped_sends++;
                continue;
            }
            
            if (dir >= d_domain->set->nDirections) {
                skipped_sends++;
                continue;
            }
            
            MPI_Request req;
            size_t tag_data[2] = {node_idx, dir};
            
            MPI_Isend(tag_data, 2, MPI_UNSIGNED_LONG, 
                     static_cast<int>(messenger.d_p), 0, comm, &req);
            send_requests.push_back(req);
            
            MPI_Isend(&messenger.d_src, 1, MPI_DOUBLE, 
                     static_cast<int>(messenger.d_p), 1, comm, &req);
            send_requests.push_back(req);
            total_sends++;
        }
        
        if (!send_requests.empty())
        {
            std::vector<MPI_Status> send_statuses(send_requests.size());
            MPI_Waitall(static_cast<int>(send_requests.size()), send_requests.data(), send_statuses.data());
        }
        
        MPI_Barrier(comm);
        
        std::vector<size_t> recv_tags_node;
        std::vector<size_t> recv_tags_dir;
        std::vector<double> recv_values;
        
        int received_pairs = 0;
        int probe_count = 0;
        int consecutive_empty_probes = 0;
        const int MAX_EMPTY_PROBES = 100;
        
         while (consecutive_empty_probes < MAX_EMPTY_PROBES)
         {
             MPI_Status probe_status;
             int flag = 0;
             MPI_Iprobe(MPI_ANY_SOURCE, 0, comm, &flag, &probe_status);
             probe_count++;
             
             if (flag)
             {
                 consecutive_empty_probes = 0;
                 int source = probe_status.MPI_SOURCE;
                 
                 size_t tag_data[2];
                 MPI_Recv(tag_data, 2, MPI_UNSIGNED_LONG, source, 0, comm, MPI_STATUS_IGNORE);
                 
                 MPI_Status value_probe_status;
                 int value_flag = 0;
                 MPI_Iprobe(source, 1, comm, &value_flag, &value_probe_status);
                 
                 if (!value_flag) {
                     continue;
                 }
                 
                 double value;
                 MPI_Recv(&value, 1, MPI_DOUBLE, source, 1, comm, MPI_STATUS_IGNORE);
                 
                 recv_tags_node.push_back(tag_data[0]);
                 recv_tags_dir.push_back(tag_data[1]);
                 recv_values.push_back(value);
                 
                 received_pairs++;
             }
             else
             {
                 consecutive_empty_probes++;
                 if (consecutive_empty_probes < 5) {
                 } else if (consecutive_empty_probes < 20) {
                     usleep(50);
                 } else {
                     usleep(200);
                 }
             }
         }
        
        MPI_Status final_probe_status;
        int final_flag = 0;
        MPI_Iprobe(MPI_ANY_SOURCE, 0, comm, &final_flag, &final_probe_status);
        
        int additional_received = 0;
        while (final_flag && additional_received < 1000)
        {
            int source = final_probe_status.MPI_SOURCE;
            size_t tag_data[2];
            MPI_Recv(tag_data, 2, MPI_UNSIGNED_LONG, source, 0, comm, MPI_STATUS_IGNORE);
            double value;
            MPI_Recv(&value, 1, MPI_DOUBLE, source, 1, comm, MPI_STATUS_IGNORE);
            recv_tags_node.push_back(tag_data[0]);
            recv_tags_dir.push_back(tag_data[1]);
            recv_values.push_back(value);
            received_pairs++;
            additional_received++;
            MPI_Iprobe(MPI_ANY_SOURCE, 0, comm, &final_flag, &final_probe_status);
        }

        size_t nDirections = d_domain->set->nDirections;
        size_t valid_updates = 0;
        size_t invalid_updates = 0;
        size_t hash_conversions = 0;
        
        for (size_t i = 0; i < recv_tags_node.size(); ++i)
        {
            size_t node_idx = recv_tags_node[i];
            size_t dir = recv_tags_dir[i];
            
            if (node_idx >= nNodes && !d_domain->map_to_index.empty())
            {
                if (d_domain->map_to_index.find(node_idx) != d_domain->map_to_index.end())
                {
                    node_idx = d_domain->map_to_index[node_idx];
                    hash_conversions++;
                }
            }
            
            if (node_idx < nNodes && dir < nDirections)
            {
                // SoA update
                d_domain->distribution_nextValues[node_idx * nDirections + dir] = recv_values[i];
                valid_updates++;
            }
            else
            {
                invalid_updates++;
            }
        }
        
        MPI_Barrier(comm);
    }

    void Simulation::report()
    {
        int total_p, s;
        MPI_Comm comm = d_domain->comm;
        MPI_Comm_size(comm, &total_p);
        MPI_Comm_rank(comm, &s);
        
        double current_density = 0.0;
        size_t nNodes = d_domain->distribution_values.size() / d_domain->set->nDirections;

        #pragma omp parallel for reduction(+:current_density) schedule(static)
        for (size_t idx = 0; idx < nNodes; ++idx)
        {
            current_density += density(d_domain->set, d_domain->distribution_values, idx);
        }

        double *densities = new double[total_p];
        MPI_Gather(&current_density, 1, MPI_DOUBLE, 
                   densities, 1, MPI_DOUBLE, 0, comm);

        if (s == 0)
        {
            double total_density = 0;
            for (int t = 0; t < total_p; t++)
                total_density += densities[t];
            std::cout << "Total density: " << total_density << '\n';
        }

        delete[] densities;
    }

}
