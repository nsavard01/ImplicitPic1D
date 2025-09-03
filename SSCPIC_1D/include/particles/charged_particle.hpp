
#pragma once

#include <vector>
#include <cstddef>
#include <string>
#include "domain/domain.hpp"
#include "domain/non_uniform_domain.hpp"
#include "domain/uniform_domain.hpp"
#include "particles/target_particle.hpp"

class charged_particle {
    
public:
    std::string name;
    size_t total_number_particles;
    int target_particle_idx = -1;
    double mass, charge, q_over_m, q_times_wp, average_density, average_temperature, del_t_max;
    double v_sqr_min, v_sqr_max; // for bin diagnostics
    double total_sum_v_square[3];
    double total_sum_v[3];
    double accum_wall_energy_loss[2], accum_wall_freq_rel[2];
    size_t accum_wall_loss[2];
    bool null_collision_bool;
    std::vector<size_t> final_idx, number_particles;
    std::vector<std::vector<std::vector<double>>> particle_components; // release_frequency_weight, xi, y, z, v_x, v_y, v_z ; in 1D weight is # particles/m^2/s
    std::vector<std::vector<double>> energy_loss, wall_freq_rel;
    std::vector<std::vector<size_t>> wall_loss;
    std::vector<std::vector<double>> density_grid;
    std::vector<std::vector<std::vector<double>>> v_sqr_grid, v_drift_grid;
    std::vector<double> total_density_grid;
    std::vector<std::vector<double>> total_v_sqr_grid, total_v_drift_grid;

    
    charged_particle(double mass_in, double charge_in, size_t number_in, std::string name_in, int number_nodes);
    void gather_mpi();
    void print_out() const;
    void load_to_target(std::vector<target_particle>& target_particle_list, const domain& world);
    void read_initial_state(const std::string& dir_name, const domain& world);
    void initialize_diagnostic_files(const std::string& dir_name) const; 
    void write_diagnostics(const std::string& dir_name, int diag_number) const; 
    void reset_diagnostics(int thread_id);

};

std::vector<charged_particle> read_charged_particle_inputs(const std::string& filename, const domain& world);
void find_corresponding_targets(std::vector<charged_particle>& charged_particle_list, std::vector<target_particle>& target_particle_list);



