#include <mpi.h>
#include <omp.h>
#include <iostream>

int main(int argc, char** argv) {
    // Initialize MPI
    MPI_Init(&argc, &argv);

    int world_size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Set the number of OpenMP threads per MPI rank
    const int num_threads = 4; // Change this to control OpenMP threads per MPI rank
    omp_set_num_threads(num_threads);

    // Start parallel region using OpenMP
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        int total_threads = omp_get_num_threads();

        // Each thread prints its info
        #pragma omp critical
        std::cout << "MPI Rank " << rank << " of " << world_size
                  << " | OpenMP Thread " << thread_id << " of " << total_threads
                  << std::endl;
    }

    // Finalize MPI
    MPI_Finalize();
    return 0;
}