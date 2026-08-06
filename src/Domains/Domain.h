#ifndef INCLUDED_DOMAINS_DOMAIN
#define INCLUDED_DOMAINS_DOMAIN

#include "../LBM/node.h"
#include "../LBM/Messenger.h"
#include "../VelocitySets/VelocitySet.h"
#include "../BoundaryConditions/BoundaryNode.h"
#include "../BoundaryConditions/PostProcessor.h"
#include <memory>
#include <unordered_map>
#include <mpi.h>

struct Domain {
    // SoA data structures
    std::vector<double> distribution_values;      // nNodes * nDirections
    std::vector<double> distribution_nextValues;  // nNodes * nDirections
    std::vector<double*> distribution_neighbours; // nNodes * nDirections
    std::vector<size_t> position;                 // nNodes * nDimensions (flat)

    // Keep this for now to match other parts of the code until refactored? 
    // Actually plan is to remove it, so I am removing it.
    // std::vector<Node> nodes; 
    
    // Change this to PostProcessor *
    std::vector<BoundaryNode> b_nodes;
    VelocitySet *set;
    double omega;
    std::vector<std::unique_ptr<BoundaryConditions::PostProcessor>> post_processors;

    // maybe add an ifdef statement ?
    std::vector<Messenger> messengers;
    
    // hashIdx -> nodeIdx mapping for converting received hash values to local node indices
    std::unordered_map<size_t, size_t> map_to_index;
    
    // messenger_idx -> (sender_node_idx, dir) mapping for inferring receiver nodes
    std::unordered_map<size_t, std::pair<size_t, size_t>> messenger_to_sender;
    
    // Domain size for calculating positions
    std::vector<size_t> domain_size;
    
    // MPI communicator (Cartesian topology if available, otherwise MPI_COMM_WORLD)
    MPI_Comm comm;

    // Halo exchange: per-source recv counts, filled once at init via Alltoall
    std::vector<int> halo_recv_counts;
};

#endif