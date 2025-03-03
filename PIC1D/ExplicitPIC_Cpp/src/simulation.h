#ifndef SIMULATION_H
#define SIMULATION_H

#include "particle.h"
#include "potential_solver.h"
#include "null_collision.h"
#include "domain.h"
#include <string>


class Simulation {
public:
    std::string directory_name;
    size_t current_diag_step, diag_step_diff, current_step;
    double diag_time_division, current_time, diag_time, elapsed_time, charge_loss, energy_loss, momentum_total[3],
        tot_particle_time, tot_potential_time, tot_collision_time, collision_energy_loss, momentum_loss;
    double start_time_total, end_time_total;
    Simulation(const std::string dir_name);
    void initialize_data_files(Potential_Solver& solver, std::vector<Particle>& particle_list, std::vector<Target_Particle>& target_particle_list,
        std::vector<Null_Collision>& binary_collision_list, const Domain& world);
    void reset_diag(std::vector<Particle>& particle_list, std::vector<Null_Collision>& binary_collision_list);
    void diag_write(Potential_Solver& solver, std::vector<Particle>& particle_list, std::vector<Target_Particle>& target_particle_list,
    std::vector<Null_Collision>& binary_collision_list, const Domain& world);
    void run(Potential_Solver& solver, std::vector<Particle>& particle_list, std::vector<Target_Particle>& target_particle_list,
        std::vector<Null_Collision>& binary_collision_list, const Domain& world);
};

#endif 
