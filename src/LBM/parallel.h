#ifndef INCLUDED_LBM_PARALLEL
#define INCLUDED_LBM_PARALLEL

#include <mpi.h>

// MPI compatibility layer - replaces BSP functions
// This file provides MPI equivalents for BSP functions used in the codebase

// MPI initialization (replaces bsp_init)
// Note: MPI_Init should be called in main() before using any MPI functions

// MPI process management (replaces bsp_begin/bsp_end)
// Note: These are handled directly with MPI_Comm_size and MPI_Comm_rank

// MPI communication functions
// bsp_send() -> MPI_Send() or MPI_Isend()
// bsp_sync() -> MPI_Barrier()
// bsp_get_tag()/bsp_move() -> MPI_Recv() or MPI_Irecv()

// Helper macros for compatibility (if needed)
#define MPI_COMM_WORLD_DEFAULT MPI_COMM_WORLD

#endif
