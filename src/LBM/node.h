#ifndef INCLUDED_LBM_NODE
#define INCLUDED_LBM_NODE

#include <cstring>
#include <vector>
#include "../VelocitySets/VelocitySet.h"

// AoS structs - Deprecated/Removed for SoA refactor
// struct Distribution
// {
//     double value;
//     double nextValue;
//     double * neighbour; // stream_destination , should point to Distribution.value
// };
// 
// struct Node
// {
//     Distribution *distributions;
//     size_t *position;
// };

// Helper functions for SoA
// Computes density from distribution values
double density(VelocitySet *set, const std::vector<double>& start_values, size_t node_idx);

// Computes velocity from distribution values
// Returns a pointer to a new array (needs delete[])
double *velocity(VelocitySet *set, const std::vector<double>& start_values, size_t node_idx);

// Computes equilibrium distribution
// Returns a pointer to a new array (needs delete[])
double *equilibrium(VelocitySet *set, const std::vector<double>& start_values, size_t node_idx);

// Legacy compatibility (if needed, but trying to remove)
// double *equilibrium(VelocitySet *set, Node node);
// double *velocity(VelocitySet *set, Node node);
// double density(VelocitySet *set, Node node);

#endif