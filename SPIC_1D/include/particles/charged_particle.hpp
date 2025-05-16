
#pragma once

#include <vector>
#include <cstddef>
#include <string>
#include "domain/domain.hpp"

class charged_particle {
    
public:
    std::string name;
    size_t total_number_particles;
    int number_space_coordinates, number_velocity_coordinates;
    double mass, charge, weight, q_over_m, q_times_wp, average_density, average_temperature;
    double total_sum_v_square;
    double accum_wall_energy_loss[2];
    size_t accum_wall_loss[2];
    std::vector<size_t> number_particles_per_cell;
    std::vector<std::vector<double>> weights;
    std::vector<std::vector<double>> xi, y, z, v_x, v_y, v_z;
    std::vector<std::vector<double>> energy_loss, accum_wall_momentum_loss;
    std::vector<std::vector<std::vector<double>>> momentum_loss;
    std::vector<std::vector<size_t>> number_particles, number_collidable_particles, wall_loss, final_idx;
    std::vector<double> density, total_sum_v, temperature;
    // Static variables for sorting, will be used per particle
    static std::vector<std::vector<double>> xi_sorted, y_sorted, z_sorted, v_x_sorted, v_y_sorted, v_z_sorted;
    static std::vector<std::vector<size_t>> sorted_number_particles_per_cell, cell_indices; //, particle_cell;
    charged_particle(double mass_in, double charge_in, size_t number_in, size_t final_in, std::string name_in, int number_nodes);
    void get_diagnostics();
    void sort_particles(int thread_id, int number_cells);
    void print_out() const;
    void initialize_number_coordinates(int space, int velocity);
    void initialize_weight(double n_ave, double L_domain);
    void initialize_rand_maxwellian(double T_ave, double v_drift);
    void initialize_rand_position_uniform(const domain& world);
    void ES_push_MC(int thread_id, double del_t, const std::vector<double>& E_field, 
        const double dx, const int left_boundary, const int right_boundary, int number_cells);
    void ES_push_EC_uniform(int thread_id, double del_t, const std::vector<double>& E_field, 
        const double dx, const int left_boundary, const int right_boundary, int number_cells);
    void ES_push_EC_non_uniform(int thread_id, double del_t, const std::vector<double>& E_field, 
        const std::vector<double>& dx_dxi, const std::vector<double>& grid, const int left_boundary, const int right_boundary, int number_cells);    
    void deposit_particles_linear(int thread_id, std::vector<double>& work_space);
    // double get_KE_ave() const;
    // double get_KE_total() const;
    // void interpolate_particles();
    // double get_momentum_total() const;
    // void write_cell_temperature(const std::string& dir_name, int diag_num) const;
    // void initialize_diagnostic_file(const std::string& dir_name) const; 
    // void diag_write(const std::string& dir_name, const double& time_diff, const double& current_time, bool average_bool = false) const;
    // void load_density(bool reset_bool);
    // void write_density(const std::string& dir_name, const domain& world, size_t current_diag, bool average_bool);
    // void gather_mpi();
    // void write_phase_space(const std::string& dir_name, int diag_num) const;
};

std::vector<charged_particle> read_charged_particle_inputs(const std::string& filename, const domain& world);


