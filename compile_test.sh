#!/bin/bash

# Compilation test script
# Used to verify that the MPI version can compile successfully

echo "=== LBM MPI Version Compilation Test ==="
echo ""

# Check if mpic++ is available
if ! command -v mpic++ &> /dev/null
then
    echo "Error: mpic++ compiler not found"
    echo "Please install MPI development libraries (e.g., sudo apt-get install libopenmpi-dev)"
    exit 1
fi

echo "✓ Found mpic++ compiler"
echo ""

# Clean previous compilation
echo "Cleaning previous compilation..."
make clean 2>/dev/null || true
echo ""

# Test compilation of vertical decomposition version
echo "Testing compilation of vertical decomposition version (default)..."
if make lbm 2>&1 | tee /tmp/lbm_compile.log; then
    echo "✓ Vertical decomposition version compiled successfully"
else
    echo "✗ Vertical decomposition version compilation failed"
    echo "View compilation log: /tmp/lbm_compile.log"
    exit 1
fi
echo ""

# Test compilation of horizontal decomposition version
echo "Testing compilation of horizontal decomposition version..."
if make lbm_horizontal 2>&1 | tee /tmp/lbm_horizontal_compile.log; then
    echo "✓ Horizontal decomposition version compiled successfully"
else
    echo "✗ Horizontal decomposition version compilation failed"
    echo "View compilation log: /tmp/lbm_horizontal_compile.log"
    exit 1
fi
echo ""

# Test compilation of 2D block decomposition version
echo "Testing compilation of 2D block decomposition version..."
if make lbm_2d 2>&1 | tee /tmp/lbm_2d_compile.log; then
    echo "✓ 2D block decomposition version compiled successfully"
else
    echo "✗ 2D block decomposition version compilation failed"
    echo "View compilation log: /tmp/lbm_2d_compile.log"
    exit 1
fi
echo ""

echo "=== All versions compiled successfully! ==="
echo ""
echo "Executable locations:"
echo "  - Vertical decomposition: bin/main"
echo "  - Horizontal decomposition: bin/lbm_horizontal (if Makefile supports it)"
echo "  - 2D block decomposition: bin/lbm_2d (if Makefile supports it)"
echo ""
echo "Run example:"
echo "  mpirun -np 4 ./bin/main 4 80 80 100"
