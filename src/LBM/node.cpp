#include "node.h"

// Computes the total density on a node
double density(VelocitySet *set, const std::vector<double>& start_values, size_t node_idx)
{
    size_t nDirections = set->nDirections;
    size_t offset = node_idx * nDirections;
    double density = 0;
    for (size_t dir = 0; dir < nDirections; ++dir)
        density += start_values[offset + dir];

    return density;
}

// Computes the projected velocity in each dimension
double *velocity(VelocitySet *set, const std::vector<double>& start_values, size_t node_idx)
{
    size_t nDirections = set->nDirections;
    size_t nDimensions = set->nDimensions;
    size_t offset = node_idx * nDirections;
    
    double nodeDensity = density(set, start_values, node_idx);
    double *velocity = new double[nDimensions]();

    // compute the velocity in each dimension taking in account
    // the form of our velocity set
    for (size_t dir = 0; dir < nDirections; ++dir)
    {
        double distribution = start_values[offset + dir];
        for (size_t dim = 0; dim < nDimensions; ++dim)
            velocity[dim] += distribution * set->weight(dir) * set->direction(dir)[dim];
    }

    if (nodeDensity != 0) {
        for (size_t dim = 0; dim < nDimensions; ++dim)
            velocity[dim] /= nodeDensity;
    }

    return velocity;
}

double *equilibrium(VelocitySet *set, const std::vector<double>& start_values, size_t node_idx)
{
    size_t nDirections = set->nDirections;
    size_t nDimensions = set->nDimensions;
    double speedOfSoundSquared = set->speedOfSoundSquared();
    double node_density = density(set, start_values, node_idx);
    double *node_velocity = velocity(set, start_values, node_idx);

    // Pre calculate the speed of the node
    double speedSquared = 0;
    for (size_t dim = 0; dim < nDimensions; ++dim)
        speedSquared += node_velocity[dim] * node_velocity[dim];
    speedSquared /= (2 * speedOfSoundSquared);

    double *equilibrium = new double[nDirections];
    for (size_t dir = 0; dir < nDirections; ++dir)
    {
        double cu = 0;
        for (size_t dim = 0; dim < nDimensions; ++dim)
            cu = set->direction(dir)[dim] * node_velocity[dim];
        cu /= speedOfSoundSquared;

        // compressible
        equilibrium[dir] = node_density * set->weight(dir) * (
            1.0 +
            cu +
            0.5 * cu * cu -
            speedSquared
        );
    }

    delete[] node_velocity;
    return equilibrium;
}