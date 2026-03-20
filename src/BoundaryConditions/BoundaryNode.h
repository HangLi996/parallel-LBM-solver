#ifndef INCLUDED_BOUNDARY_CONDITIONS_BOUNDARY_NODE
#define INCLUDED_BOUNDARY_CONDITIONS_BOUNDARY_NODE

// Temporary (until i have a better idea of what a boundarynode is)
#include "ZouHe.h"

struct BoundaryNode {
    // Node * node;
    size_t node_idx;
    Direction direction;
    union {
        double *velocity;
        double density;
        // Node *interpolateFrom;
        size_t interpolateFrom_idx;
    };
};

#endif