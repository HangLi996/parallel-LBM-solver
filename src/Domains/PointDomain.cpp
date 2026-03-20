#include "PointDomain.h"

namespace Domains {

    double PointDomain::omega()
    {
        return 1.0 / 0.6;
    }

    void PointDomain::initializeNodeAt(std::vector<int> &position)
    {
        // Call base class to push default values
        DomainInitializer::initializeNodeAt(position);
        
        // Modify the last pushed values if condition met
        if (
            (position[0] == 10 && position[1] == 10) ||
            (position[0] == 10 && position[1] == 30) ||
            (position[0] == 30 && position[1] == 30) ||
            (position[0] == 30 && position[1] == 10) ||
            (position[0] == 5 && position[1] == 5) ||
            (position[0] == 5 && position[1] == 25) ||
            (position[0] == 25 && position[1] == 25) ||
            (position[0] == 25 && position[1] == 5)
            )
        {
            size_t nDirections = d_set->nDirections;
            size_t nNodes = d_distribution_values.size() / nDirections;
            size_t idx = nNodes - 1; // The node we just added
            size_t offset = idx * nDirections;

            // Set distributions
            for (size_t dir = 0; dir < nDirections; ++dir)
            {
                d_distribution_values[offset + dir] *= 10;
                d_distribution_nextValues[offset + dir] *= 10;
            }
        }
    }

    void PointDomain::connectNodeToNeighbours(size_t idx)
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
            else
            {
                size_t neighbour_idx = idxOf(neighbour);
                // d_nodes[idx].distributions[dir].neighbour = &d_nodes[neighbour_idx].distributions[dir].nextValue;
                 size_t neighbor_offset = neighbour_idx * nDirections;
                d_distribution_neighbours[offset + dir] = &d_distribution_nextValues[neighbor_offset + dir];
            }
        }
    }

    bool PointDomain::isInDomain(std::vector<int> &position)
    {
        return (
            position[0] >= 0 && position[0] < d_domain_size[0] &&
            position[1] >= 0 && position[1] < d_domain_size[1]
        ) && true;
    }

    bool PointDomain::isBounceBack(std::vector<int> position)
    {
        return not isInDomain(position);
    }
}