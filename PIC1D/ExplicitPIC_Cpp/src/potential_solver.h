#ifndef POTENTIAL_SOLVER_H
#define POTENTIAL_SOLVER_H

#include <vector>
#include "domain.h"
#include "particle.h"

class Potential_Solver {
public:
    std::vector<double> upper_tri, lower_tri, center_tri, phi, rho, EField, work_space;
    double rho_const, RF_rad_frequency, RF_half_amplitude;
    Potential_Solver();
    void read_from_file(const std::string& filename, const Domain& world);
    void deposit_rho(std::vector<Particle> &particle_list, const Domain& world);
    void solve_potential_tridiag(const Domain& world, const double& time);
    void make_EField(const Domain& world);
    void initial_v_rewind(std::vector<Particle> &particle_list, const double& time_step);
    inline double get_EField_at_loc(const double& xi) const;
    void move_particles(std::vector<Particle> &particle_list, const Domain& world, const double& time_step);
};

#endif // PARTICLE_H