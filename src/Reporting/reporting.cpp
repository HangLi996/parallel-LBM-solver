#include <iostream>
#include <fstream>
#include <sstream>
#include "reporting.h"
#include "MatlabReporter.h"
#include "../LBM/node.h"

namespace Reporting {

    void reportOnInitialSetup(VelocitySet *set, size_t dx, size_t dy)
    {
        std::cout << "LBM D" << set->nDimensions << "Q" << set->nDirections << " simulation." << '\n'
            << "Grid size " << dx << "x" << dy
            << '\n' << '\n';
    }

    /*
    void reportOnInitialSetup(VelocitySet *set, Node *nodes, size_t dx, size_t dy)
    {
        reportOnInitialSetup(set, dx, dy);

        for (int x = 0; x < dx; ++x)
            for (int y = 0; y < dy; ++y)
            {
                int idx = x + y * dx;
                std::cout << "(" << x << ", " << y << ")" << ": ";
                reportOnNode(set, nodes[idx]);
            }

    }
    

    void reportOnNode(VelocitySet *set, Node &node)
    {
        double node_density = density(set, node);
        double *node_velocity = velocity(set, node);
        size_t nDimensions = set->nDimensions;

        std::cout << "Density: " << node_density << '\n';
        std::cout << "Velocity: ";
        for (size_t dir = 0; dir < nDimensions; ++dir)
            std::cout << node_velocity[dir] << ", ";
        std::cout << '\n';
        std::cout << "Position: ";
        for (size_t dir = 0; dir < nDimensions; ++dir)
            std::cout << node.position[dir] << ", ";
        
        std::cout << '\n' << '\n';
        reportOnDistributions(set, node);
        std::cout << "---------------------------------------" << '\n';

        delete[] node_velocity;
    }

    void reportOnDistributions(VelocitySet *set, Node &node)
    {
        size_t nDirections = set->nDirections;
        std::stringstream ss;
        ss << "Distributions for node at (" << node.position[0] << ", " << node.position[1] << "), (";
        for (size_t dir = 0; dir < nDirections; ++dir)
            ss << node.distributions[dir].value << ", ";
        ss << ")" << '\n';

        ss << "Next distributions for node at (" << node.position[0] << ", " << node.position[1] << "), (";
        for (size_t dir = 0; dir < nDirections; ++dir)
            ss << node.distributions[dir].nextValue << ", ";
        ss << ")" << '\n';

        ss << "Distributions neighbours for node at (" << node.position[0] << ", " << node.position[1] << "), (";
        for (size_t dir = 0; dir < nDirections; ++dir)
            ss << *node.distributions[dir].neighbour << ", ";
        ss << ")" << '\n';

        std::cout << ss.str();
    }

    void report(VelocitySet *set, Node *nodes, size_t dx, size_t dy)
    {
        report(set, nodes, dx * dy);
    }
    
    void report(VelocitySet *set, Node *nodes, size_t totalNodes)
    {
        double total_density = 0;
        
        std::cout << "---------------------------------------" << '\n';
        std::cout << "Reporting on simulation" << '\n';
        std::cout << "---------------------------------------" << '\n';

        // we don't know the rank of the current processor, so we can't
        // print it here.
        // int s;
        // MPI_Comm_rank(MPI_COMM_WORLD, &s);
        // std::cout << "Processor " << s << " reporting..." << '\n';

        
        #pragma omp parallel for reduction(+:total_density) schedule(static)
        for (size_t idx = 0; idx < totalNodes; ++idx)
        {
            total_density += density(set, nodes[idx]);
            double *node_velocity = velocity(set, nodes[idx]); //Might create race condition if used. Currently safe bcs its not used
            delete[] node_velocity;
        }

        std::cout << "Total density: " << total_density << '\n';
        std::cout << "---------------------------------------" << '\n';

    }

    void report(std::string outputFileName, VelocitySet *set, Node *nodes, size_t totalNodes)
    {
        std::ofstream output;
        output.open(outputFileName, std::ios::out);
        MatlabReporter reporter(output);
        reporter.reportOnTimeStep(set, nodes, totalNodes);
        output.close();
    }
    */
}