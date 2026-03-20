#!/bin/bash

# Configuration
# Load MPI module
module load openmpi/4.1.8-gcc13

ITERATIONS=50
OUTPUT_FILE="profiling/benchmark_results.csv"

# Create output file with header
echo "domain_size,computation_time" > $OUTPUT_FILE

echo "Running benchmark..."
echo "Iterations: $ITERATIONS"
echo "Processes: 1"

# Loop from 100 to 800 with step 100
for SIZE in {100..800..100}
do
    echo -n "Running domain size ${SIZE}x${SIZE}... "
    
    # Run simulation and capture output
    # Using 1 process as requested
    OUTPUT=$(mpirun -np 1 ./bin/main 1 $SIZE $SIZE $ITERATIONS)
    
    # Extract computation time using grep and awk
    # Looking for "Total computation time: X seconds"
    TIME=$(echo "$OUTPUT" | grep "Total computation time:" | awk '{print $4}')
    
    if [ -z "$TIME" ]; then
        echo "Error: Could not extract time for size $SIZE"
        continue
    fi
    
    echo "$TIME seconds"
    
    # Append to result file
    echo "$SIZE,$TIME" >> $OUTPUT_FILE
done

echo "Benchmark complete. Results saved to $OUTPUT_FILE"
