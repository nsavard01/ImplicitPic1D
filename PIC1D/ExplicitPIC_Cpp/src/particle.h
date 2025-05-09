#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include <cstddef>
#include <string>
#include "domain.h"
#include "pcg_rng.h"

class Particle {
public:
    double mass, charge, weight, q_over_m, q_times_wp;
    int accum_wall_loss[2];
    std::string name;
    double accum_energy_loss[2], accum_momentum_loss[2];
    size_t final_idx, total_number_particles;
    double total_sum_v_square, total_sum_v_x;
    std::vector<size_t> number_particles, number_collidable_particles;
    std::vector<std::vector<double>> phase_space, work_space, momentum_loss, energy_loss;
    std::vector<std::vector<int>> wall_loss;
    std::vector<double> density;  // Flattened 3D array
    Particle(double mass_in, double charge_in, size_t number_in, size_t final_in, std::string name_in);
    void initialize_weight(double n_ave, double L_domain);
    void initialize_rand_uniform(double T_ave, const Domain& world);
    double get_KE_ave() const;
    double get_KE_total() const;
    void interpolate_particles();
    double get_momentum_total() const;
    void write_cell_temperature(const std::string& dir_name, int diag_num) const;
    void initialize_diagnostic_file(const std::string& dir_name) const; 
    void diag_write(const std::string& dir_name, const double& time_diff, const double& current_time, bool average_bool = false) const;
    void load_density(bool reset_bool);
    void write_density(const std::string& dir_name, const Domain& world, size_t current_diag, bool average_bool);
    void gather_mpi();
    void write_phase_space(const std::string& dir_name, int diag_num) const;
    // Accessor for 3D indexing
};
std::vector<Particle> read_particle_inputs(const std::string& filename, const Domain& world);

#endif // PARTICLE_H
