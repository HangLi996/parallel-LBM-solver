#!/bin/bash

# Load Likwid module
module load likwid/5.3.0-gcc13-perf
module load openmpi/4.1.8-gcc13

# Check if likwid-perfctr is available
if ! command -v likwid-perfctr &> /dev/null; then
    echo "Error: likwid-perfctr could not be found. Please ensure the likelihood module is loaded correctly."
    exit 1
fi

# Compile the code
echo "Compiling..."
make lbm

# Configuration
# Choose group: FLOPS_DP (Double Precision), FLOPS_AVX (AVX instructions), MEM (Memory bandwidth)
# Default to FLOPS_DP as MEM might require privileges not available on all nodes
GROUP=${1:-"FLOPS_DP"}

DOMAIN_SIZE=100
ITERATIONS=50
OUTPUT_FILE="likwid_${GROUP}.txt"

echo "Running Likwid profiling for group: $GROUP"
echo "Output will be saved to: $OUTPUT_FILE"

# Run with likwid-perfctr
# -C 0: Pin to core 0
# -g $GROUP: Performance group
# -m: Enable marker API (if code is instrumented, otherwise it profiles the whole execution)
# -f: Force overwrite of output file (implied by -o usually, but explicit good)

likwid-perfctr -C 0 -g "$GROUP" -f -o "$OUTPUT_FILE" ./bin/main 1 "$DOMAIN_SIZE" "$DOMAIN_SIZE" "$ITERATIONS"

echo "Likwid profiling complete."
echo "Results:"
grep "AVX DP" "$OUTPUT_FILE" || grep "Packed" "$OUTPUT_FILE" || head -n 20 "$OUTPUT_FILE"
