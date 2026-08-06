#include "Simulation.h"
#include "node.h"
#include <iostream>
#include <sstream>
#include <mpi.h>
#include <vector>
#include <map>
#include <cstdlib>      // free()

#include "../LBM/parallel.h"
#include <omp.h>


namespace LBM {

namespace {

struct OutboundMsg {
    int dest;
    size_t tag[2];
    double *value;
};

// Same send-set used at init (for recv counts) and every communicate() step.
std::vector<OutboundMsg> buildOutbound(Domain *domain)
{
    auto& messengers = domain->messengers;
    size_t nDirections = domain->set->nDirections;
    size_t nDimensions = domain->set->nDimensions;
    size_t nNodes = domain->distribution_values.size() / nDirections;

    std::vector<OutboundMsg> outbound;
    outbound.reserve(messengers.size());

    for (size_t msg_idx = 0; msg_idx < messengers.size(); ++msg_idx)
    {
        auto& messenger = messengers[msg_idx];
        size_t node_idx = messenger.d_tag[0];
        size_t dir = messenger.d_tag[1];
        int dest = static_cast<int>(messenger.d_p);

        if (node_idx == SIZE_MAX || node_idx >= nNodes)
        {
            auto it = domain->messenger_to_sender.find(msg_idx);
            if (it == domain->messenger_to_sender.end())
                continue;

            size_t sender_node_idx = it->second.first;
            size_t sender_dir = it->second.second;
            if (sender_node_idx >= nNodes || sender_dir >= nDirections)
                continue;

            size_t pos_offset = sender_node_idx * nDimensions;
            if (pos_offset + nDimensions > domain->position.size())
                continue;

            std::vector<int> sender_position;
            sender_position.reserve(nDimensions);
            for (size_t dim = 0; dim < nDimensions; ++dim)
                sender_position.push_back(static_cast<int>(domain->position[pos_offset + dim]));

            auto dir_vec = domain->set->direction(sender_dir);
            std::vector<int> receiver_position;
            receiver_position.reserve(sender_position.size());
            for (size_t dim = 0; dim < sender_position.size(); ++dim)
            {
                int domain_size_dim = static_cast<int>(domain->domain_size[dim]);
                int coord = (sender_position[dim] + dir_vec[dim] + domain_size_dim) % domain_size_dim;
                receiver_position.push_back(coord);
            }

            size_t receiver_hash = receiver_position[0];
            size_t multiplier = domain->domain_size[0];
            for (size_t dim = 1; dim < receiver_position.size(); ++dim)
            {
                receiver_hash += receiver_position[dim] * multiplier;
                multiplier *= domain->domain_size[dim];
            }

            outbound.push_back(OutboundMsg{dest, {receiver_hash, dir}, &messenger.d_src});
            continue;
        }

        if (dir >= nDirections)
            continue;

        outbound.push_back(OutboundMsg{dest, {node_idx, dir}, &messenger.d_src});
    }

    return outbound;
}

} // namespace

    Simulation::Simulation(Initializer_Ptr initializer)
    :
        d_domain(initializer->domain()),
        d_stream_time(0.0),
        d_communicate_time(0.0),
        d_collision_time(0.0),
        d_poststream_time(0.0),
        d_step_count(0)
    {
        cacheHaloRecvCounts();
        MPI_Barrier(d_domain->comm);
    }

    Simulation::~Simulation()
    {
        // cleanup handled by vectors
    }

    void Simulation::cacheHaloRecvCounts()
    {
        int size;
        MPI_Comm comm = d_domain->comm;
        MPI_Comm_size(comm, &size);

        auto outbound = buildOutbound(d_domain.get());
        std::vector<int> send_counts(static_cast<size_t>(size), 0);
        for (const auto &msg : outbound)
        {
            if (msg.dest >= 0 && msg.dest < size)
                send_counts[static_cast<size_t>(msg.dest)]++;
        }

        d_domain->halo_recv_counts.assign(static_cast<size_t>(size), 0);
        MPI_Alltoall(send_counts.data(), 1, MPI_INT,
                     d_domain->halo_recv_counts.data(), 1, MPI_INT, comm);
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
        int size;
        MPI_Comm comm = d_domain->comm;
        MPI_Comm_size(comm, &size);

        size_t nNodes = d_domain->distribution_values.size() / d_domain->set->nDirections;
        size_t nDirections = d_domain->set->nDirections;

        auto outbound = buildOutbound(d_domain.get());
        const auto &recv_counts = d_domain->halo_recv_counts;

        int n_recv = 0;
        for (int c : recv_counts)
            n_recv += c;

        std::vector<size_t> recv_tags(static_cast<size_t>(n_recv) * 2);
        std::vector<double> recv_values(static_cast<size_t>(n_recv));
        std::vector<MPI_Request> requests;
        requests.reserve(static_cast<size_t>(n_recv) * 2 + outbound.size() * 2);

        // Post all Irecvs first (deadlock-free with subsequent Isends).
        size_t recv_slot = 0;
        for (int src = 0; src < size; ++src)
        {
            int count = (static_cast<size_t>(src) < recv_counts.size())
                            ? recv_counts[static_cast<size_t>(src)] : 0;
            for (int i = 0; i < count; ++i)
            {
                MPI_Request req_tag, req_val;
                MPI_Irecv(recv_tags.data() + recv_slot * 2, 2, MPI_UNSIGNED_LONG,
                         src, 0, comm, &req_tag);
                MPI_Irecv(recv_values.data() + recv_slot, 1, MPI_DOUBLE,
                         src, 1, comm, &req_val);
                requests.push_back(req_tag);
                requests.push_back(req_val);
                recv_slot++;
            }
        }

        for (auto &msg : outbound)
        {
            MPI_Request req_tag, req_val;
            MPI_Isend(msg.tag, 2, MPI_UNSIGNED_LONG, msg.dest, 0, comm, &req_tag);
            MPI_Isend(msg.value, 1, MPI_DOUBLE, msg.dest, 1, comm, &req_val);
            requests.push_back(req_tag);
            requests.push_back(req_val);
        }

        if (!requests.empty())
            MPI_Waitall(static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE);

        for (int i = 0; i < n_recv; ++i)
        {
            size_t node_idx = recv_tags[static_cast<size_t>(i) * 2];
            size_t dir = recv_tags[static_cast<size_t>(i) * 2 + 1];

            if (node_idx >= nNodes && !d_domain->map_to_index.empty())
            {
                auto it = d_domain->map_to_index.find(node_idx);
                if (it != d_domain->map_to_index.end())
                    node_idx = it->second;
            }

            if (node_idx < nNodes && dir < nDirections)
                d_domain->distribution_nextValues[node_idx * nDirections + dir] = recv_values[static_cast<size_t>(i)];
        }
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
