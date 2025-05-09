
#pragma once

#include <vector>
#include <cstddef>
#include <string>
#include "domain/domain.hpp"

class charged_particle {
private:
    std::string name;
    size_t total_number_particles, number_cells;
    double mass, charge, weight, q_over_m, q_times_wp;
    double total_sum_v_square, total_sum_v_x;
    double accum_wall_energy_loss[2], accum_wall_momentum_loss[2];
    size_t accum_wall_loss[2];
    std::vector<size_t> number_particles_per_cell;
    std::vector<std::vector<double>> x, y, z, v_x, v_y, v_z;
    std::vector<std::vector<double>> cell_v_sqr, momentum_loss, energy_loss;
    std::vector<std::vector<size_t>> number_particles, number_collidable_particles, wall_loss, final_idx, cell_idx_array;
    std::vector<double> density;
public:
    charged_particle(double mass_in, double charge_in, size_t number_in, size_t final_in, std::string name_in, int number_nodes);
    void initialize_weight(double n_ave, double L_domain);
    void initialize_rand_uniform(double T_ave, const domain& world);
    double get_KE_ave() const;
    double get_KE_total() const;
    void interpolate_particles();
    double get_momentum_total() const;
    void write_cell_temperature(const std::string& dir_name, int diag_num) const;
    void initialize_diagnostic_file(const std::string& dir_name) const; 
    void diag_write(const std::string& dir_name, const double& time_diff, const double& current_time, bool average_bool = false) const;
    void load_density(bool reset_bool);
    void write_density(const std::string& dir_name, const domain& world, size_t current_diag, bool average_bool);
    void gather_mpi();
    void write_phase_space(const std::string& dir_name, int diag_num) const;
};

std::unique_ptr<std::vector<charged_particle>> read_charged_particle_inputs(const std::string& filename, const domain& world);


