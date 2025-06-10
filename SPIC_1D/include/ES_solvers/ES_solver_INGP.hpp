
#pragma once
#include <vector>
#include "solvers/poisson_solver_1D.hpp"
#include "non_linear_solvers/AA_solver.hpp"
#include "ES_solvers/ES_solver.hpp"
#include "domain/domain.hpp"

class ES_solver_INGP : public ES_solver {


public:
    std::vector<double> phi_past, J, rho_past; // future potential for implicit solver
    std::unique_ptr<non_linear_solver> implicit_solver; // pointer to the implicit non-linear solver
    ES_solver_INGP(const domain& world);
    void print_out() override;
    void make_EField(const domain& world) override;
    void push_particles(int thread_id, double del_t, std::vector<charged_particle>& particle_list, const domain& world) override;
    void integrate_time_step(int thread_id, double del_t, double current_time, const domain& world, std::vector<charged_particle>& particle_list) override;
    // void integrate_time(double current_time, double del_t) override;
    // void interpolate_particles_to_grid(std::vector<charged_particle>& particle_list) override;
};

