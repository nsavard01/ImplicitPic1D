
#pragma once
#include <vector>
#include "solvers/poisson_solver_1D.hpp"
#include "ES_solvers/ES_solver.hpp"
#include "domain/domain.hpp"

class ES_solver_MC : public ES_solver {


public:
    ES_solver_MC(const domain& world);
    // void solve(double current_time, double del_t) override;
    // void interpolate_particles_to_grid(std::vector<charged_particle>& particle_list) override;
};

