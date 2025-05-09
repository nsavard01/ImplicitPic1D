#include "globals/constants.hpp"
#include "globals/mpi_vars.hpp"
#include "domain/domain.hpp"
#include "domain/uniform_domain.hpp"
#include "domain/non_uniform_domain.hpp"
#include <stdio.h>
#include <iostream>
#include <mpi.h>



int main(int argc, char** argv) {

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_vars::mpi_size);  // Get total processes
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_vars::mpi_rank);  // Get current rank
    std::unique_ptr<domain> world = domain::create_from_file("../inputs/geometry.inp");
    if (mpi_vars::mpi_rank == 0) {
        world->print_out();
    }
    MPI_Finalize();
    return 0;
}