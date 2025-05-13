
#include <vector>
#include <omp.h>
#include "ES_solvers/ES_solver_MC.hpp"
#include "solvers/poisson_solver_1D_tridiag.hpp"

ES_solver_MC::ES_solver_MC(const domain& world) {
    this->work_space.resize(omp_get_max_threads());
    for (int i = 0; i < omp_get_max_threads(); i++){
        this->work_space[i].resize(world.get_number_nodes(), 0.0);
    }
    this->phi.resize(world.get_number_nodes(), 0.0);
    this->rho.resize(world.get_number_nodes(), 0.0);
    this->poisson_solver = std::make_unique<poisson_solver_1D_tridiag>(world);
}



