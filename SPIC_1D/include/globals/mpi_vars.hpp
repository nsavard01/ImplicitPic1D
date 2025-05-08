#pragma once
#include <mpi.h>

// global physical constants
namespace mpi_vars {
    extern int mpi_rank, mpi_size;
    extern MPI_Datatype mpi_size_t_type;
}