#include "DomainInitializer.h"
#include <iostream>
#include <sstream>
#include <math.h>       /* ceil */
#include <mpi.h>
#include <mpi.h>
#include <cstdlib>      // posix_memalign
#include <new>          // std::bad_alloc
#include <iterator>

#include "../LBM/parallel.h"

namespace Domains {

    DomainInitializer::DomainInitializer(VelocitySet *set, std::vector<size_t> domainSize, size_t p, size_t totalProcessors,
                                         MPI_Comm cart_comm, int* cart_dims)
    :
        d_p(p),
        d_total_processors(totalProcessors),
        d_set(set),
        d_domain_size(domainSize),
        d_cart_comm(cart_comm),
        d_use_cart(cart_comm != MPI_COMM_NULL && cart_dims != nullptr)
    {
        if (d_use_cart) {
            d_cart_dims[0] = cart_dims[0];
            d_cart_dims[1] = cart_dims[1];
        } else {
            d_cart_dims[0] = 0;
            d_cart_dims[1] = 0;
        }
    }

    DomainInitializer::~DomainInitializer()
    {}

    std::unique_ptr<Domain> DomainInitializer::domain()
    {
        MPI_Comm comm = d_use_cart ? d_cart_comm : MPI_COMM_WORLD;
        // Synchronize all processes before domain initialization
        MPI_Barrier(comm);

        createNodes();

        std::unique_ptr<Domain> domain(new Domain);
        
        // Move SoA vectors to domain
        domain->distribution_values = std::move(d_distribution_values);
        domain->distribution_nextValues = std::move(d_distribution_nextValues);
        domain->distribution_neighbours = std::move(d_distribution_neighbours);
        domain->position = std::move(d_position);

        // Now that we've moved the nodes to our domain object, we can create post
        // processors which can point to these nodes
        createPostProcessors(domain.get());
        domain->post_processors = std::move(d_post_processors);

        domain->set    = d_set;
        domain->omega = omega();
        
        // Copy the hash-to-index mapping for use in communication
        domain->map_to_index = d_map_to_index;
        
        // Copy domain size for use in communication
        domain->domain_size = d_domain_size;
        
        // Store communicator for use in Simulation
        domain->comm = d_use_cart ? d_cart_comm : MPI_COMM_WORLD;

        // Setup messengers (for the parallelisation of the code)
        // Synchronize all processes before messenger setup
        MPI_Barrier(comm);
        
        // In MPI, we need to receive messages sent during sendLocationOfDistribution
        int rank, size;
        MPI_Comm_rank(comm, &rank);
        MPI_Comm_size(comm, &size);

        if (size > 1)
        {
            // Probe for all incoming messages (multi-process case)
            int flag = 1;
            int received_count = 0;
            int matched_count = 0;
            int unmatched_count = 0;
            
            while (flag)
            {
                MPI_Status probe_status;
                MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &probe_status);
                
                if (probe_status.MPI_TAG == 0)
                {
                    size_t tag_hash;
                    size_t localIdx;
                    
                    MPI_Recv(&tag_hash, 1, MPI_UNSIGNED_LONG, 
                            probe_status.MPI_SOURCE, 0, comm, &probe_status);
                    
                    MPI_Recv(&localIdx, 1, MPI_UNSIGNED_LONG, 
                            probe_status.MPI_SOURCE, 0, comm, &probe_status);
                    
                    received_count++;
                    
                    if (d_map_to_messenger.find(tag_hash) != d_map_to_messenger.end())
                    {
                        size_t messenger_idx = d_map_to_messenger[tag_hash];
                        
                        size_t multiplier = d_domain_size[0];
                        for (size_t dim = 1; dim < d_domain_size.size(); ++dim)
                            multiplier *= d_domain_size[dim];
                        
                        size_t dir = tag_hash / multiplier;
                        size_t sender_position_hash = tag_hash % multiplier;
                        
                        std::vector<int> sender_position;
                        size_t temp_hash = sender_position_hash;
                        for (size_t dim = 0; dim < d_domain_size.size(); ++dim)
                        {
                            sender_position.push_back(static_cast<int>(temp_hash % d_domain_size[dim]));
                            temp_hash /= d_domain_size[dim];
                        }
                        
                        std::vector<int> receiver_position;
                        auto dir_vec = d_set->direction(dir);
                        for (size_t dim = 0; dim < d_domain_size.size(); ++dim)
                        {
                            int coord = (static_cast<int>(sender_position[dim]) + dir_vec[dim] + static_cast<int>(d_domain_size[dim])) % static_cast<int>(d_domain_size[dim]);
                            receiver_position.push_back(coord);
                        }
                        
                        size_t receiver_hash = hashIdxOf(receiver_position);
                        if (d_map_to_index.find(receiver_hash) != d_map_to_index.end()) {
                            size_t our_local_idx = d_map_to_index[receiver_hash];
                            d_messengers[messenger_idx].d_tag[0] = our_local_idx;
                            matched_count++;
                        } else {
                            unmatched_count++;
                        }
                    }
                    else
                    {
                        unmatched_count++;
                    }
                }
                
                MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &flag, &probe_status);
            }
        }
        
        // Pointer fixup logic
        struct SenderPointer {
            size_t node_idx;
            size_t dir;
            size_t messenger_idx;
        };
        std::vector<SenderPointer> sender_pointers;
        size_t nDirections = d_set->nDirections;
        size_t nNodes = domain->distribution_values.size() / nDirections; // Calculated from size

        // Use d_messengers logic. 
        // Note: domain->distribution_neighbours already contains pointers. 
        // Some point to d_messengers. We need to identify them.

        for (size_t node_idx = 0; node_idx < nNodes; ++node_idx)
        {
            for (size_t dir = 0; dir < nDirections; ++dir)
            {
                double* neighbour_ptr = domain->distribution_neighbours[node_idx * nDirections + dir];
                if (neighbour_ptr == nullptr) continue;
                
                for (size_t msg_idx = 0; msg_idx < d_messengers.size(); ++msg_idx)
                {
                    if (neighbour_ptr == &d_messengers[msg_idx].d_src)
                    {
                        SenderPointer sp;
                        sp.node_idx = node_idx;
                        sp.dir = dir;
                        sp.messenger_idx = msg_idx;
                        sender_pointers.push_back(sp);
                        break;
                    }
                }
            }
        }

        // Move messengers from deque to vector
        domain->messengers.assign(std::make_move_iterator(d_messengers.begin()), 
                                 std::make_move_iterator(d_messengers.end()));

        for (const auto& sender : sender_pointers)
        {
            if (sender.messenger_idx < domain->messengers.size())
            {
                domain->messenger_to_sender[sender.messenger_idx] = std::make_pair(sender.node_idx, sender.dir);
            }
        }
        
        // Update sender pointers
        for (const auto& sender : sender_pointers)
        {
            if (sender.messenger_idx < domain->messengers.size())
            {
                 domain->distribution_neighbours[sender.node_idx * nDirections + sender.dir] = &domain->messengers[sender.messenger_idx].d_src;
            }
        }
        
        // Update receiver pointers
        for (size_t idx = 0; idx < domain->messengers.size(); ++idx)
        {
            size_t node_idx = domain->messengers[idx].d_tag[0];
            size_t dir = domain->messengers[idx].d_tag[1];
            
            if (node_idx != SIZE_MAX && node_idx < nNodes && dir < nDirections)
            {
                domain->distribution_neighbours[node_idx * nDirections + dir] = &domain->messengers[idx].d_src;
            }
        }
        
        // Final verification (simplified from original for brevity, but retaining core logic)
        for (size_t node_idx = 0; node_idx < nNodes; ++node_idx)
        {
            for (size_t dir = 0; dir < nDirections; ++dir)
            {
                double* neighbour_ptr = domain->distribution_neighbours[node_idx * nDirections + dir];
                if (neighbour_ptr == nullptr) continue;
                
                // Check if stale
                for (size_t msg_idx = 0; msg_idx < d_messengers.size(); ++msg_idx)
                {
                    if (neighbour_ptr == &d_messengers[msg_idx].d_src)
                    {
                        if (msg_idx < domain->messengers.size()) {
                            domain->distribution_neighbours[node_idx * nDirections + dir] = &domain->messengers[msg_idx].d_src;
                        }
                        break;
                    }
                }
            }
        }

        return domain;
    }

    double DomainInitializer::omega()
    {
        return 1.0;
    }

    void DomainInitializer::createNodes()
    {
        size_t potentialTotalNodes = 1;
        for (auto sizeOfDimension : d_domain_size)
            potentialTotalNodes *= sizeOfDimension;

        // Reserve memory? potentialTotalNodes is upper bound? No, it's exact size of domain box.
        // d_nodes.reserve(...) would be good.
        // For SoA vectors:
        // size_t estNodes = potentialTotalNodes / d_total_processors; // Rough estimate
        
        for (size_t idx = 0; idx < potentialTotalNodes; ++idx)
        {
            size_t currentIndex = idx;

            std::vector<int> position;
            position.push_back(currentIndex % d_domain_size[0]);
            for (size_t dim = 1; dim < d_domain_size.size(); ++dim)
            {
                currentIndex = (currentIndex - position[dim - 1]) / d_domain_size[dim - 1];
                position.push_back(currentIndex % d_domain_size[dim]);
            }

            if (isInDomain(position) && processorOfNode(position) == d_p)
            {
                d_map_to_index[idx] = d_position.size() / d_domain_size.size(); // Current node count
                initializeNodeAt(position);
            }
        }

        size_t nNodes = d_position.size() / d_domain_size.size();
        for (size_t idx = 0; idx < nNodes; ++idx)
            connectNodeToNeighbours(idx);

        // Connect messengers
        for (size_t idx = 0; idx < d_messengers.size(); ++idx)
        {
            size_t node_idx = d_messengers[idx].d_tag[0];
            size_t dir = d_messengers[idx].d_tag[1];
            if (node_idx != SIZE_MAX && node_idx < nNodes && dir < d_set->nDirections)
            {
                d_distribution_neighbours[node_idx * d_set->nDirections + dir] = &d_messengers[idx].d_src;
            }
        }
    }

    void DomainInitializer::initializeNodeAt(std::vector<int> &position)
    {
        size_t nDirections = d_set->nDirections;
        size_t nDimensions = d_set->nDimensions;

        // Set position (flat)
        for (size_t dim = 0; dim < nDimensions; ++dim)
            d_position.push_back(position[dim]);

        // Set distributions
        for (size_t dir = 0; dir < nDirections; ++dir)
        {
            d_distribution_values.push_back(d_set->weight(dir));
            d_distribution_nextValues.push_back(d_set->weight(dir));
            d_distribution_neighbours.push_back(nullptr);
        }
    }

    void DomainInitializer::connectNodeToNeighbours(size_t idx)
    {
        size_t nDirections = d_set->nDirections;
        size_t offset = idx * nDirections;
        
        // Reconstruct position from flat vector
        std::vector<int> node_position;
        size_t pos_offset = idx * d_domain_size.size();
        for(size_t dim=0; dim<d_domain_size.size(); ++dim)
            node_position.push_back(d_position[pos_offset + dim]);

        for (size_t dir = 0; dir < nDirections; ++dir)
        {
            std::vector<int> neighbour;
            for (size_t dim = 0; dim < d_domain_size.size(); ++dim)
            {
                neighbour.push_back((
                    node_position[dim] + d_set->direction(dir)[dim] + d_domain_size[dim]
                ) % d_domain_size[dim]);
            }

            sendLocationOfDistribution(idx, dir);

            size_t p = processorOfNode(neighbour);
            if (p == d_p)
            {
                size_t neighbour_hash = hashIdxOf(neighbour);
                if (d_map_to_index.find(neighbour_hash) != d_map_to_index.end())
                {
                    size_t neighbour_idx = d_map_to_index[neighbour_hash];
                    // Verify range?
                    
                    // Point to nextValue of neighbor
                    size_t neighbor_offset = neighbour_idx * nDirections;
                    d_distribution_neighbours[offset + dir] = &d_distribution_nextValues[neighbor_offset + dir];
                }
            }
            else
            {
                std::vector<int> sender_position = node_position;
                
                size_t messenger_hash = hashIdxOf(sender_position, dir);
                d_map_to_messenger[messenger_hash] = d_messengers.size();
                d_messengers.push_back(create_messenger(p, dir));
                
                if (d_map_to_index.find(hashIdxOf(neighbour)) != d_map_to_index.end()) {
                    d_messengers.back().d_tag[0] = d_map_to_index[hashIdxOf(neighbour)];
                } else {
                    d_messengers.back().d_tag[0] = SIZE_MAX;
                }
                
                d_distribution_neighbours[offset + dir] = &d_messengers.back().d_src;
            }
        }
    }

    void DomainInitializer::sendLocationOfDistribution(size_t node_idx, size_t dir)
    {
        // Reconstruct position
        std::vector<int> node_position;
        size_t pos_offset = node_idx * d_domain_size.size();
        for(size_t dim=0; dim<d_domain_size.size(); ++dim)
            node_position.push_back(d_position[pos_offset + dim]);

        std::vector<int> neighbour;
        for (size_t dim = 0; dim < d_domain_size.size(); ++dim)
        {
            neighbour.push_back((
                node_position[dim] - d_set->direction(dir)[dim] + d_domain_size[dim]
            ) % d_domain_size[dim]);
        }

        size_t p = processorOfNode(neighbour);
        if (p == d_p)
            return;

        auto tag = hashIdxOf(node_position, dir);
        size_t src = idxOf(node_position); // should be same as node_idx
        
        MPI_Comm comm = d_use_cart ? d_cart_comm : MPI_COMM_WORLD;
        MPI_Send(&tag, 1, MPI_UNSIGNED_LONG, static_cast<int>(p), 0, comm);
        MPI_Send(&src, 1, MPI_UNSIGNED_LONG, static_cast<int>(p), 0, comm);
    }

    size_t DomainInitializer::processorOfNode(std::vector<int> &position)
    {
        if (d_use_cart && d_domain_size.size() >= 2)
        {
            int px = d_cart_dims[0];
            int py = d_cart_dims[1];
            int block_x = static_cast<int>((static_cast<double>(position[0]) * px) / d_domain_size[0]);
            int block_y = static_cast<int>((static_cast<double>(position[1]) * py) / d_domain_size[1]);
            
            if (block_x >= px) block_x = px - 1;
            if (block_y >= py) block_y = py - 1;
            if (block_x < 0) block_x = 0;
            if (block_y < 0) block_y = 0;
            
            int coords[2] = {block_x, block_y};
            int rank;
            MPI_Cart_rank(d_cart_comm, coords, &rank);
            return static_cast<size_t>(rank);
        }
        
#ifdef DECOMP_HORIZONTAL
        double p = static_cast<double>(d_total_processors * position[1]) / d_domain_size[1];
        return static_cast<size_t>(floor(p));
#elif defined(DECOMP_2D)
        int px = static_cast<int>(sqrt(static_cast<double>(d_total_processors)));
        int py = px;
        if (px * py != static_cast<int>(d_total_processors))
        {
            double p = static_cast<double>(d_total_processors * position[0]) / d_domain_size[0];
            return static_cast<size_t>(floor(p));
        }
        int block_x = static_cast<int>(position[0] / (d_domain_size[0] / px));
        int block_y = static_cast<int>(position[1] / (d_domain_size[1] / py));
        if (block_x >= px) block_x = px - 1;
        if (block_y >= py) block_y = py - 1;
        return static_cast<size_t>(block_y * px + block_x);
#else
        double p = static_cast<double>(d_total_processors * position[0]) / d_domain_size[0];
        return static_cast<size_t>(floor(p));
#endif
    }

    void DomainInitializer::createPostProcessors(Domain* domain)
    {
    }

    bool DomainInitializer::isInDomain(std::vector<int> &position)
    {
        return true;
    }

    size_t DomainInitializer::hashIdxOf(std::vector<int> & position, size_t direction)
    {
        if (position.size() != d_domain_size.size())
            throw std::string("Position is not compatible with domain size");

        size_t hashIdx = position[0];
        size_t multiplier = d_domain_size[0];

        for (size_t dim = 1; dim < d_domain_size.size(); ++dim)
        {
            hashIdx += multiplier * position[dim];
            multiplier *= d_domain_size[dim];
        }
        hashIdx += direction * multiplier;

        return hashIdx;
    }

    size_t DomainInitializer::idxOf(std::vector<int> &position)
    {
        size_t hashIdx = hashIdxOf(position);
        return d_map_to_index[hashIdx];
    }
}
