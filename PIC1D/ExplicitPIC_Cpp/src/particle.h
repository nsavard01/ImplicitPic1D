#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include <cstddef>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include "domain.h"
#include "pcg_rng.h"

class Particle {
public:
    double mass, charge, weight, q_over_m, q_times_wp;
    int accum_wall_loss[2];
    std::string name;
    double accum_energy_loss[2];
    size_t final_idx;
    std::vector<size_t> number_particles, number_collidable_particles;
    std::vector<std::vector<double>> phase_space, work_space, momentum_loss, energy_loss;
    std::vector<std::vector<int>> wall_loss;
    std::vector<double> densities;  // Flattened 3D array
    Particle(double mass_in, double charge_in, size_t number_in, size_t final_in, std::string name_in);
    void initialize_weight(double n_ave, double L_domain);
    void initialize_rand_uniform(double T_ave, const Domain& world);
    double get_KE_ave() const;
    double get_KE_total() const;
    void interpolate_particles();
    double get_momentum_total() const;
    void write_cell_temperature(const std::string& dir_name, int diag_num) const;
    void initialize_diagnostic_file(const std::string& dir_name) const; 
    // Accessor for 3D indexing
};
std::vector<Particle> read_particle_inputs(const std::string& filename, const Domain& world);

#endif // PARTICLE_H
