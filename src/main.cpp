#include "main.h"
#include <ctime>
#include <fstream>
#include <sstream>
#include <iostream>

#include "Domains/DomainInitializer.h"
#include "Domains/BoxedDomain.h"
#include "Domains/LidDrivenCavity.h"
#include "Domains/PointDomain.h"
#include "VelocitySets/d2q9.h"

#include "LBM/node.h"
#include "LBM/parallel.h"

#include <memory>
#include <mpi.h>

using namespace Domains;
using std::make_unique;

// Global variables
size_t dx = 80;
size_t dy = 80;
size_t ITERATIONS;
size_t REPORT_PER_ITERATION = 10;
size_t P;

int main(int argc, char **argv)
try
{
    // Initialize MPI environment
    MPI_Init(&argc, &argv);
    
    int rank, world_size;
    // Get current process rank in MPI_COMM_WORLD
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // Get total number of processes in MPI_COMM_WORLD
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    
    P = askForProcessors(argc, argv);
    
    // Read domain size from command line arguments
    // Format: ./bin/main <P> <dx> <dy> <iterations>
    if (argc > 2) {
        std::istringstream iss(argv[2]);
        if (!(iss >> dx)) {
            dx = 80;  // Default value
        }
    }
    if (argc > 3) {
        std::istringstream iss(argv[3]);
        if (!(iss >> dy)) {
            dy = 80;  // Default value
        }
    }
    
    // Read iterations (now at argv[4] in new format)
    ITERATIONS = askForIterations(argc, argv);

    // Check if we have enough processes
    if (P > static_cast<size_t>(world_size))
    {
        if (world_size == 0)
        {
            // Finalize MPI environment and exit
            MPI_Finalize();
            return 1;
        }
        P = world_size;
    }

    simulate();
    
    // Finalize MPI environment
    MPI_Finalize();
    return 0;
}
catch(int x)
{
    // Finalize MPI environment on error
    MPI_Finalize();
    return x;
}
catch(...)
{
    // Finalize MPI environment on exception
    MPI_Finalize();
    return 1;
}

void simulate()
{
    int p, s;
    // Get total number of processes in MPI_COMM_WORLD
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    // Get current process rank in MPI_COMM_WORLD
    MPI_Comm_rank(MPI_COMM_WORLD, &s);

    // Use the actual number of MPI processes
    size_t total_processors = static_cast<size_t>(p);
    size_t current_process = static_cast<size_t>(s);

    std::vector<size_t> domainSize{dx, dy};

    if (s == 0) {
        logSimulationData(domainSize);
    }

    // Create MPI Cartesian topology
    MPI_Comm cart_comm;
    int ndims = 2;  // 2D domain
    int dims[2];
    int periods[2] = {1, 1};  // Periodic boundaries in both dimensions
    int reorder = 1;  // Allow MPI to reorder processes for optimization
    
    // Determine optimal grid dimensions for Cartesian topology
    dims[0] = 0;
    dims[1] = 0;
    // Calculate optimal dimensions for Cartesian grid decomposition
    MPI_Dims_create(p, ndims, dims);
    
    // Create Cartesian communicator with periodic boundaries
    // Creates a new communicator with Cartesian topology from MPI_COMM_WORLD
    int result = MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &cart_comm);
    if (result != MPI_SUCCESS) {
        // Finalize MPI environment on error
        MPI_Finalize();
        return;
    }
    
    // If this process is not in the Cartesian topology, skip it
    if (cart_comm == MPI_COMM_NULL) {
        // Finalize MPI environment if excluded from topology
        MPI_Finalize();
        return;
    }
    
    // Get Cartesian coordinates for this process
    int coords[2];
    // Get Cartesian coordinates of the calling process in the Cartesian communicator
    MPI_Cart_coords(cart_comm, s, ndims, coords);

    // Synchronize all processes before initialization
    MPI_Barrier(cart_comm);
    
    // Get current wall-clock time for timing measurements
    // Returns elapsed wall-clock time since some time in the past
    double initialization_time = MPI_Wtime();

    D2Q9 set;
    LBM::Simulation sim(make_unique<DomainInitializer>(&set, domainSize, current_process, total_processors, cart_comm, dims));

    // Synchronize all processes after initialization
    MPI_Barrier(cart_comm);
    // Get current wall-clock time
    double current_time = MPI_Wtime();
    double init_time = current_time - initialization_time;
    
    if (s == 0)
    {
        // Output simulation configuration and performance metrics
        std::cout << "=== LBM Simulation Configuration ===" << std::endl;
        std::cout << "Processes: " << p << std::endl;
        std::cout << "Domain size: " << dx << "x" << dy << std::endl;
        std::cout << "Iterations: " << ITERATIONS << std::endl;
        std::cout << "Relaxation parameter (omega): " << sim.getOmega() << std::endl;
        std::cout << "Cartesian topology: " << dims[0] << "x" << dims[1] << std::endl;
        std::cout << "Initialization time: " << init_time << " seconds" << std::endl;
        
        std::ofstream out("logs/timings.log", std::ios::out | std::ios::app);
        // Initialization time
        out << "IT: " << init_time << " sec, ";
    }
    
    // Get current wall-clock time for computation timing
    double process_time = MPI_Wtime();

    run_simulation(sim);

    // Synchronize all processes after computation
    MPI_Barrier(cart_comm);
    // Get current wall-clock time
    current_time = MPI_Wtime();
    double comp_time = current_time - process_time;
    double total_time = current_time - initialization_time;
    
    if (s == 0)
    {
        // Output detailed computation performance metrics
        std::cout << "=== Computation Performance Breakdown ===" << std::endl;
        std::cout << "Total computation time: " << comp_time << " seconds" << std::endl;
        std::cout << "  - Stream time: " << sim.getStreamTime() << " seconds (" 
                  << (sim.getStreamTime() / comp_time * 100.0) << "%)" << std::endl;
        std::cout << "  - Communication time: " << sim.getCommunicateTime() << " seconds (" 
                  << (sim.getCommunicateTime() / comp_time * 100.0) << "%)" << std::endl;
        std::cout << "  - Post-stream time: " << sim.getPostStreamTime() << " seconds (" 
                  << (sim.getPostStreamTime() / comp_time * 100.0) << "%)" << std::endl;
        std::cout << "  - Collision time: " << sim.getCollisionTime() << " seconds (" 
                  << (sim.getCollisionTime() / comp_time * 100.0) << "%)" << std::endl;
        std::cout << "Total time: " << total_time << " seconds" << std::endl;
        std::cout << "==========================================" << std::endl;
        
        std::ofstream out("logs/timings.log", std::ios::out | std::ios::app);
        // computation time
        out << "CT: " << comp_time << " sec" << '\n';
    }
    
    // Free the Cartesian communicator after simulation is complete
    // Note: Domain stores a copy of the communicator, but MPI_Comm_free
    // only decrements the reference count, so it's safe to call here
    // Marks the communicator for deallocation and decrements its reference count
    MPI_Comm_free(&cart_comm);
}

/**
 * Run ITERATIONS step and periodically report the current state of the simulation
 * @param sim
 */
void run_simulation(LBM::Simulation &sim)
{
    int rank;
    // Get current process rank in MPI_COMM_WORLD
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    for (size_t iter = 0; iter < ITERATIONS; ++iter)
    {
        sim.step();

        if (iter % REPORT_PER_ITERATION == 0) {
            sim.report();
        }
    }
}

size_t askForIterations(int argc, char** argv)
{
    size_t iterations = 11;
    // New format: ./bin/main <P> <dx> <dy> <iterations>
    // So iterations is at argv[4]
    if (argc > 4) {
        std::istringstream iss(argv[4]);
        if ((iss >> iterations) )
            return iterations;
    }
    // Backward compatibility: if only 3 args, assume old format <P> <iterations>
    else if (argc > 3) {
        std::istringstream iss(argv[3]);
        size_t val;
        if ((iss >> val)) {
            // If it's a reasonable iteration count (< 10000), use it
            if (val < 10000) {
                return val;
            }
        }
    }

    return iterations;
}

size_t askForProcessors(int argc, char** argv)
{
    size_t P;
    if (argc > 1)
    {
        std::istringstream iss(argv[1]);
        if ((iss >> P) )
            return P;
    }
    else
    {
        // In MPI, we get the number of processes from MPI_Comm_size
        int world_size;
        // Get total number of processes in MPI_COMM_WORLD
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        return static_cast<size_t>(world_size);
    }

    // Check if requested processes exceed available
    int world_size;
    // Get total number of processes in MPI_COMM_WORLD
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    if (P > static_cast<size_t>(world_size))
    {
        std::cout << "Warning: Requested " << P << " processors, but only " 
                  << world_size << " available. Using " << world_size << " processors." << std::endl;
        return static_cast<size_t>(world_size);
    }
    return P;
}

void logSimulationData(std::vector<size_t> domainSize)
{
    std::ofstream out("logs/timings.log", std::ios::out | std::ios::app);
    int p;
    // Get total number of processes in MPI_COMM_WORLD
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    // Start by writing basic info to the file
    // out << "LBM simulation using " << p <<
    //     " processors to perform " << ITERATIONS << " iterations on the 'dummy' domain " <<
    //     ", with set: " << "D2Q9" << " and domain size: (";
    out << "LBM, p: " << p << ", it: " << ITERATIONS << ", ds (";
    showVector(domainSize, out);
    out << ", ";
}

void showVector(std::vector<size_t> vector, std::ofstream &out)
{
    for (size_t dim = 0; dim < (vector.size() - 1); ++dim)
        out << vector[dim] << ", ";
    out << vector[vector.size() - 1] << ")";
}

void createMatlabReport(LBM::Simulation &sim, size_t iter, std::vector<size_t> domainSize)
{
    // Create a output file
    std::ofstream out(createFileName(iter, "D2Q9", "PERIODIC", domainSize), std::ios::out | std::ios::app);
    Reporting::MatlabReporter reporter(out);
    sim.report(reporter);
}

std::string createFileName(size_t iteration, std::string setName, std::string domainName, std::vector<size_t> domainSize)
{
    // (iter, "D2Q9", "Lid_Driven_Cavity_5000_", domainSize);
    // logs/D2Q9_Lid_Driven_Cavity_5000_50_100
    // logs/D2Q9_Boxed_10x10_100.txt
    std::stringstream ss;
    ss << "logs/" << setName << "_" << domainName << "_";
    for (size_t idx = 0; idx < domainSize.size() - 1; ++idx)
        ss << domainSize[idx] << "x";
    ss << domainSize[domainSize.size() - 1];
    ss << "_" << ITERATIONS << "_" << iteration << ".txt";
    return ss.str();
}
