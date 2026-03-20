#ifndef INCLUDED_DOMAINS_DOMAIN_INITIALIZER
#define INCLUDED_DOMAINS_DOMAIN_INITIALIZER

#include <vector>
#include <unordered_map>
#include <memory>
#include <deque>

#include "Domain.h"

#include "../VelocitySets/VelocitySet.h"
#include "../LBM/node.h"
#include "../BoundaryConditions/ZouHe.h"
#include "../BoundaryConditions/PostProcessor.h"
#include "../LBM/Messenger.h"

namespace Domains {
    class DomainInitializer {
        protected:
            size_t d_p;
            size_t d_total_processors;

            VelocitySet *d_set;
            std::vector<size_t> d_domain_size;
            
            // MPI Cartesian topology
            MPI_Comm d_cart_comm;
            int d_cart_dims[2];
            bool d_use_cart;

            // hashIdx -> nodeIdx
            std::unordered_map<size_t, size_t> d_map_to_index;
            std::unordered_map<size_t, size_t> d_map_to_messenger;

            // std::vector<Node> d_nodes;
            // SoA data structures
            std::vector<double> d_distribution_values;      // nNodes * nDirections
            std::vector<double> d_distribution_nextValues;  // nNodes * nDirections
            std::vector<double*> d_distribution_neighbours; // nNodes * nDirections
            std::vector<size_t> d_position;                 // nNodes * nDimensions (flat)

            std::vector<std::unique_ptr<BoundaryConditions::PostProcessor>> d_post_processors;

            // parallel:
            std::deque<Messenger> d_messengers;

        public:
            DomainInitializer(VelocitySet *set, std::vector<size_t> domainSize, size_t p = 0, size_t totalProcessors = 1, 
                             MPI_Comm cart_comm = MPI_COMM_NULL, int* cart_dims = nullptr);
            virtual ~DomainInitializer();
            std::unique_ptr<Domain> domain();

        protected:
            void createNodes();
            virtual void initializeNodeAt(std::vector<int> &position);
            // virtual void connectNodeToNeighbours(Node &node);

            virtual void connectNodeToNeighbours(size_t idx);
            virtual bool isInDomain(std::vector<int> &position);

            virtual void createPostProcessors(Domain* domain = nullptr);

            virtual double omega();

            // In order to use a multi dimensional array, we actually
            // convert each vector onto a unique integer
            size_t hashIdxOf(std::vector<int> &position, size_t direction = 0);
            size_t idxOf(std::vector<int> &position);

            // get the destrination of a distribution after streaming step
            // can be a pointer to a node's distribution, or a pointer to a messenger
            virtual void sendLocationOfDistribution(size_t node_idx, size_t direction);
            virtual size_t processorOfNode(std::vector<int> &position);

    };
}

#endif