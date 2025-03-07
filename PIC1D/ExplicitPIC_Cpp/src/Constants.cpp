#include "Constants.h"
#include <cmath>


namespace Constants {
    const double sin_third_rot = std::sin(2.0 * M_PI / 3.0);
    const double cos_third_rot = std::cos(2.0 * M_PI / 3.0);
    int mpi_rank, mpi_size;
    MPI_Datatype mpi_size_t_type;
}