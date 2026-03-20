#include "ZouHeVelocityBoundaryCondition.h"
#include <iostream>
#include "../Domains/Domain.h"

namespace BoundaryConditions {
    ZouHeVelocityBoundary::ZouHeVelocityBoundary(std::vector<double> velocity, std::vector<size_t> acts_on_indices, Domain *domain)
    :
        d_acts_on_indices(acts_on_indices),
        d_velocity(velocity),
        d_domain(domain)
    {}

    ZouHeVelocityNorthBoundary::ZouHeVelocityNorthBoundary(std::vector<double> velocity, std::vector<size_t> acts_on_indices, Domain *domain)
    :
        ZouHeVelocityBoundary(velocity, acts_on_indices, domain)
    {}

    void ZouHeVelocityNorthBoundary::process()
    {
        double u_x = d_velocity[0];
        double u_y = d_velocity[1];
        double rho = 0;
        size_t nDirections = d_domain->set->nDirections;

        for (auto node_idx : d_acts_on_indices)
        {
            size_t offset = node_idx * nDirections;
            // Access nextValues
            // node->distributions[dir].nextValue becomes d_domain->distribution_nextValues[offset + dir]

            rho = (1 / (1 + u_y)) * (
                d_domain->distribution_nextValues[offset + 0] +
                d_domain->distribution_nextValues[offset + 1] +
                d_domain->distribution_nextValues[offset + 3] +
                2 * (
                    d_domain->distribution_nextValues[offset + 2] +
                    d_domain->distribution_nextValues[offset + 5] +
                    d_domain->distribution_nextValues[offset + 6]
                )
            );

            d_domain->distribution_nextValues[offset + 4] = d_domain->distribution_nextValues[offset + 2];

            d_domain->distribution_nextValues[offset + 7] = d_domain->distribution_nextValues[offset + 5] +
                0.5 *(d_domain->distribution_nextValues[offset + 1] - d_domain->distribution_nextValues[offset + 3] ) -
                0.5 * rho * u_x;

            d_domain->distribution_nextValues[offset + 8] = d_domain->distribution_nextValues[offset + 6] +
                0.5 * (d_domain->distribution_nextValues[offset + 3] - d_domain->distribution_nextValues[offset + 1]) +
                0.5 * rho * u_x;
        }
    }
}