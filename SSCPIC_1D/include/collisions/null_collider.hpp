
#pragma once
#include <vector>
#include <cmath>
#include "particles/charged_particle.hpp"
#include "particles/target_particle.hpp"


class null_collider {
    // null collision primary particle with background gas
public:
    int primary_idx, number_targets, number_cells;
    std::vector<std::vector<std::vector<double>>> sigma_array; // sigma and energy array for each collision type [target, coll type, sigma], for coulomb collisions sigma is size cell_number with coulomb logarithm
    std::vector<std::vector<double>> energy_threshold;
    std::vector<std::vector<std::vector<double>>> total_incident_energy, total_energy_loss; // diagnostics for each collision type [cell, target, coll type]
    std::vector<std::vector<std::vector<std::vector<double>>>> total_incident_energy_thread, total_energy_loss_thread; // diagnostics for each collision type [thread, cell, target, coll type]
    std::vector<double> energy_array, reduced_mass, reduced_mass_ionization; 
    std::vector<double> max_time_step; // inverse nu_max for each cell
    std::vector<std::vector<std::vector<int>>> product_indices; // [target, numb coll, indices]
    std::vector<std::vector<int>> collision_type_per_target; // product indices for each array [target, coll type]
    std::vector<std::vector<std::vector<int>>> collision_id; // id for collision (keep track collision ordering for each cell) [cell, target, coll type]
    std::vector<int> target_idx, number_collisions_per_target; // collision type identifier 
    std::vector<std::vector<std::vector<size_t>>> total_amount_collisions; // total amount of collisions for each type [cell, target, coll type]
    std::vector<std::vector<std::vector<std::vector<size_t>>>> total_amount_collisions_thread; // total amount of collisions for each type [thread, cell, target, coll type]
    size_t total_amount_collidable_particles;
    null_collider(int primary_idx, int number_targets, int number_cells, const std::vector<int>& target_idx, const std::vector<int>& number_collisions_per_target, 
    const std::vector<std::vector<std::vector<double>>> &sigma_array, 
    const std::vector<double> &energy_array, const std::vector<std::vector<double>> &energy_threshold,
    const std::vector<std::vector<int>> &collision_type_per_target, const std::vector<std::vector<std::vector<int>>> &product_indices, 
    const std::vector<double>& reduced_mass, const std::vector<double>& reduced_mass_ionization);
    void set_initial_null_frequency(const std::vector<charged_particle>& particle_list, const std::vector<target_particle>& target_particle_list);
    void print_out(const std::vector<charged_particle>& particle_list, const std::vector<target_particle>& target_particle_list) const;
    // // void generate_null_collisions(const int thread_id, std::vector<charged_particle> &particle_list, const std::vector<target_particle> &target_particle_list, const double time_step);
    // void order_collisions();
    // inline void double_product_isotropic(const double &primary_mass, const double &target_mass, const double &del_E, 
    //     double (&incident_velocity)[3], double (&target_velocity)[3]);
    // inline void triple_product_isotropic(const double &primary_mass, const double &ion_mass, const double &target_mass, const double &del_E, 
    //     double (&incident_velocity)[3], double (&target_velocity)[3], double (&third_velocity)[3]);
    // void initialize_diagnostic_files(const std::string& dir_name, const std::vector<charged_particle>& particle_list, const std::vector<target_particle>& target_particle_list) const;
    // void write_diagnostics(const std::string& dir_name, const std::vector<charged_particle>& particle_list, const std::vector<target_particle>& target_particle_list) const;
    // void write_diagnostics_average(const std::string& dir_name, const std::vector<charged_particle>& particle_list, const std::vector<target_particle>& target_particle_list) const;
    // void diag_write(const std::string& dir_name, std::vector<Particle>& particle_list, std::vector<Target_Particle>& target_particle_list, const double& time_diff, bool average_bool = false);
    // void gather_mpi();
    // void reset_diagnostics(int thread_id);
};

std::vector<null_collider> read_null_collision_inputs(const std::string& directory, const std::vector<charged_particle> &particle_list, const std::vector<target_particle> &target_particle_list);

