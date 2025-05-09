#include "globals/constants.hpp"
#include "globals/mpi_vars.hpp"
#include "domain/domain.hpp"
#include "domain/uniform_domain.hpp"
#include "domain/non_uniform_domain.hpp"
#include "rand_gen/pcg_rng.hpp"
#include "particles/charged_particle.hpp"
#include <stdio.h>
#include <iostream>
#include <mpi.h>
#include <omp.h>
#include <fstream>
#include <sstream>



int main(int argc, char** argv) {

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_vars::mpi_size);  // Get total processes
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_vars::mpi_rank);  // Get current rank

    // Determine the MPI datatype for size_t based on its size
    if (sizeof(size_t) == sizeof(unsigned int)) {
        mpi_vars::mpi_size_t_type = MPI_UNSIGNED;
    } else if (sizeof(size_t) == sizeof(unsigned long)) {
        mpi_vars::mpi_size_t_type = MPI_UNSIGNED_LONG;
    } else if (sizeof(size_t) == sizeof(unsigned long long)) {
        mpi_vars::mpi_size_t_type = MPI_UNSIGNED_LONG_LONG;
    } else {
        std::cout << "Unsupported size_t type!" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    int max_threads = omp_get_max_threads();

    // Get number of OpenMP threads from the file
    int number_omp_threads;
    if (mpi_vars::mpi_rank == 0) {
        std::ifstream file("../inputs/initial_setup.inp");
        if (!file) {
            std::cout << "Error: Unable to open file initial_setup.inp" << std::endl;
            exit(EXIT_FAILURE);
        }
        std::string line;
        std::getline(file, line);
        std::istringstream iss(line);
        iss >> number_omp_threads;
        file.close();
        if (number_omp_threads > max_threads) {
            std::cerr << "Error: Requested number of OpenMP threads (" << number_omp_threads
                    << ") exceeds the maximum available (" << max_threads << ")." << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&number_omp_threads, 1, MPI_INT, 0, MPI_COMM_WORLD);
    omp_set_num_threads(number_omp_threads); // Set the number of OpenMP threads
    for (int i = 0; i < mpi_vars::mpi_size; i++) {
        if (i == mpi_vars::mpi_rank) {
            std::cout << "MPI rank: " << mpi_vars::mpi_rank << " using " << omp_get_max_threads() << " OpenMP threads out of maximum of " << max_threads << std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    initialize_pcg(false); // Initialize the PCG RNG with a non-deterministic seed
    std::unique_ptr<domain> world = domain::create_from_file("../inputs/geometry.inp");
    world->print_out();
    MPI_Finalize();
    return 0;
}