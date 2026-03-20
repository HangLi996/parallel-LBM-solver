#include "BoxedDomain.h"
#include "../LBM/parallel.h"
#include <mpi.h>
#include <cmath>

namespace Domains {

    void BoxedDomain::connectNodeToNeighbours(size_t idx)
    {
        size_t nDirections = d_set->nDirections;
        size_t offset = idx * nDirections;

        // Reconstruct position
        std::vector<int> position;
        size_t pos_offset = idx * d_domain_size.size();
        for (size_t dim = 0; dim < d_domain_size.size(); ++dim)
            position.push_back(d_position[pos_offset + dim]);

        for (size_t dir = 0; dir < nDirections; ++dir)
        {
            std::vector<int> neighbour;
            for (size_t dim = 0; dim < d_domain_size.size(); ++dim)
                neighbour.push_back(
                    position[dim] + d_set->direction(dir)[dim]
                );

            if (isBounceBack(neighbour))
            {
                size_t op_dir =d_set->oppositeDirectionOf(dir);
                d_distribution_neighbours[offset + dir] = &d_distribution_nextValues[offset + op_dir];
            }
            // else (todo)
            //     d_nodes[idx].distributions[dir].neighbour = destination(neighbour, dir, idx);
            sendLocationOfDistribution(idx, dir);
        }
    }

    // get the node pointing to this distribution and if it is not in
    // the current processor, then send the source of this distribution
    // to that processor
    void BoxedDomain::sendLocationOfDistribution(size_t node_idx, size_t dir)
    {
        // Reconstruct position
        std::vector<int> position;
        size_t pos_offset = node_idx * d_domain_size.size();
        for (size_t dim = 0; dim < d_domain_size.size(); ++dim)
            position.push_back(d_position[pos_offset + dim]);

        std::vector<int> neighbour;
        for (size_t dim = 0; dim < d_domain_size.size(); ++dim)
        {
            // get the neighbour in this direction, using periodic boundary
            neighbour.push_back(position[dim] - d_set->direction(dir)[dim]);
        }

        size_t p = processorOfNode(neighbour);
        if (p == d_p || isBounceBack(neighbour))
            return;

        // we send the local index of the node to the messenger
        
        // tag should contain the position and direction
        // the tag tells us where the messenger is located
        auto tag = hashIdxOf(position, dir);
        size_t src = idxOf(position);
        
        // Use MPI instead of BSP
        MPI_Send(&tag, 1, MPI_UNSIGNED_LONG, static_cast<int>(p), 0, MPI_COMM_WORLD); // Blocking send tag hash
        MPI_Send(&src, 1, MPI_UNSIGNED_LONG, static_cast<int>(p), 0, MPI_COMM_WORLD); // Blocking send sender's local index
    }

    size_t BoxedDomain::processorOfNode(std::vector<int> &position)
    {
        // Use the same domain decomposition strategy as DomainInitializer
        // This should match the strategy used in DomainInitializer
#ifdef DECOMP_HORIZONTAL
        // Horizontal decomposition: divide along y-direction
        double p = static_cast<double>(d_total_processors * position[1]) / d_domain_size[1];
        return static_cast<size_t>(floor(p));
#elif defined(DECOMP_2D)
        // 2D block decomposition
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
        // Default: vertical decomposition (divide along x-direction)
        if (d_total_processors < 2)
            return 0;
        // Simple vertical split for compatibility
        if (position[0] < d_domain_size[0] / 2)
            return 0;
        else
            return 1;
#endif
    }

    bool BoxedDomain::isInDomain(std::vector<int> &position)
    {
        return (
            position[0] >= 0 && position[0] < d_domain_size[0] &&
            position[1] >= 0 && position[1] < d_domain_size[1]
        );
    }

    bool BoxedDomain::isBounceBack(std::vector<int> position)
    {
        return not isInDomain(position);
    }
}
