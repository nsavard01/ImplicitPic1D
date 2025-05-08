#include "globals/constants.hpp"
#include "globals/mpi_vars.hpp"
#include "domain/uniform_domain.hpp"
#include <stdio.h>
#include <iostream>
#include <mpi.h>



int main(int argc, char** argv) {

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_vars::mpi_size);  // Get total processes
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_vars::mpi_rank);  // Get current rank
    uniform_domain world(10, 1.0, 0, 0); // Create a uniform domain with 10 cells and length 1.0
    MPI_Finalize();
    return 0;
}