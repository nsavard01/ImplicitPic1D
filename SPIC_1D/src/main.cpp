#include "globals/constants.hpp"
#include "globals/mpi_vars.hpp"
#include <stdio.h>
#include <iostream>
#include <mpi.h>



int main(int argc, char** argv) {

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_vars::mpi_size);  // Get total processes
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_vars::mpi_rank);  // Get current rank
    std::cout << "mpi_rank " << mpi_vars::mpi_rank << std::endl;
    MPI_Finalize();
    return 0;
}