#include "LidDrivenCavity.h"
#include "../BoundaryConditions/ZouHeVelocityBoundaryCondition.h"


using namespace BoundaryConditions;

namespace Domains {

    void LidDrivenCavityDomain::connectNodeToNeighbours(size_t idx)
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

            if (isZouHe(position) && pointsOutwards(neighbour))
            {
                d_distribution_neighbours[offset + dir] = nullptr;
            }
            else if (isBounceBack(neighbour))
            {
                size_t op_dir = d_set->oppositeDirectionOf(dir);
                // Point to own nextValue in opposite direction
                d_distribution_neighbours[offset + dir] = &d_distribution_nextValues[offset + op_dir];
            }
            else
            {
                // Periodic boundary
                for (size_t dim = 0; dim < d_domain_size.size(); ++dim)
                    neighbour[dim] = (neighbour[dim] + d_domain_size[dim]) % d_domain_size[dim];

                size_t neighbour_idx = idxOf(neighbour);
                // Point to neighbor's nextValue in same direction
                // d_nodes[neighbour_idx].distributions[dir].nextValue
                size_t neighbor_offset = neighbour_idx * nDirections;
                d_distribution_neighbours[offset + dir] = &d_distribution_nextValues[neighbor_offset + dir];
            }
        }
    }

    bool LidDrivenCavityDomain::isInDomain(std::vector<int> &position)
    {
        // square domain
        return (
            position[0] >= 0 && position[0] < d_domain_size[0] &&
            position[1] >= 0 && position[1] < d_domain_size[1]
        );
    }

    bool LidDrivenCavityDomain::isBounceBack(std::vector<int> position)
    {
        return not isInDomain(position) && not isZouHe(position);
    }

    bool LidDrivenCavityDomain::isZouHe(std::vector<int> position)
    {
        return position[1] == d_domain_size[1]; // top wall
    }

    bool LidDrivenCavityDomain::pointsOutwards(std::vector<int> position)
    {
        return position[1] < 0; // anything above our top wall is not periodic
    }

    void LidDrivenCavityDomain::createPostProcessors(Domain* domain)
    {
        // Get all Nodes that are on the moving wall
        std::vector<size_t> acts_on_indices;
        int y = d_domain_size[1] - 1;
        std::vector<int> position {0, y};
        for (size_t x = 0; x < d_domain_size[0]; ++x)
        {
            position[0] = x;
            acts_on_indices.push_back(idxOf(position));
        }
        // Create a new post processor
        std::vector<double> velocity = {0.05, 0};
        d_post_processors.push_back(
            std::unique_ptr<PostProcessor> (new ZouHeVelocityNorthBoundary(velocity, acts_on_indices, domain))
        );
    }
}