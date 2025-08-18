
#pragma once

#include <vector>
#include <cstddef>
#include <string>
#include "domain/domain.hpp"

class charged_particle {
    
public:
    std::string name;
    size_t total_number_particles;
    double mass, charge, q_over_m, q_times_wp, average_density, average_temperature;
    double v_sqr_min, v_sqr_max; // for bin diagnostics
    double total_sum_v_square[3];
    double total_sum_v[3];
    double accum_wall_energy_loss[2];
    size_t accum_wall_loss[2];
    std::vector<size_t> final_idx, number_particles;
    std::vector<std::vector<std::vector<double>>> particle_components; // release_frequency_weight x, y, z, v_x, v_y, v_z
    std::vector<std::vector<double>> energy_loss;
    std::vector<std::vector<size_t>> wall_loss;
    std::vector<std::vector<double>> density_grid, v_sqr_grid;

    
    charged_particle(double mass_in, double charge_in, size_t number_in, std::string name_in, int number_nodes);
    void gather_mpi();
    void print_out() const;  
    void initialize_diagnostic_files(const std::string& dir_name) const; 
    void write_diagnostics(const std::string& dir_name, int diag_number) const; 
    void reset_diagnostics(int thread_id);

};

std::vector<charged_particle> read_charged_particle_inputs(const std::string& filename, const domain& world);


