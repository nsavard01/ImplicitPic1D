#ifndef POTENTIAL_SOLVER_H
#define POTENTIAL_SOLVER_H

#include <vector>
#include "domain.h"
#include "particle.h"

class Potential_Solver {
public:
    std::vector<double> upper_tri, lower_tri, center_tri, phi, rho, EField, work_space;
    double rho_const, RF_rad_frequency, RF_half_amplitude;
    Potential_Solver(const std::string& filename, const Domain& world);
    void deposit_rho(std::vector<Particle> &particle_list, const Domain& world);
    void solve_potential_tridiag(const Domain& world, double time);
};

#endif // PARTICLE_H