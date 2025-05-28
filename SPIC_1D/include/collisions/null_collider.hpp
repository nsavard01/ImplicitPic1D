
#pragma once
#include <vector>
#include <cmath>
#include "particles/charged_particle.hpp"
#include "particles/target_particle.hpp"


class null_collider {
    // null collision primary particle with background gas
public:
    int number_collisions, primary_idx;
    std::vector<std::vector<double>> sigma_array, energy_array; // sigma and energy array for each collision type
    std::vector<double> energy_threshold, total_incident_energy, total_energy_loss; // diagnostics for each collision type
    double null_frequency;
    std::vector<std::vector<int>> products_indx; // product indices for each array
    std::vector<int> collision_type, target_idx; // collision type identifier
    std::vector<size_t> total_amount_collisions; // total amount of collisions for each type
    size_t total_amount_collidable_particles;
    null_collider(int number_collisions, int primary_idx, const std::vector<int>& target_idx, 
    const std::vector<std::vector<double>> &sigma_array, 
    const std::vector<std::vector<double>> &energy_array, const std::vector<double> &energy_threshold,
    const std::vector<int> &collision_type, const std::vector<std::vector<int>> &products_indx);
    void set_null_frequency(const std::vector<charged_particle>& particle_list, const std::vector<target_particle>& target_particle_list);
    // void generate_null_collisions(std::vector<Particle> &particle_list, std::vector<Target_Particle> &target_particle_list, double time_step);
    // inline void double_product_isotropic(const double &primary_mass, const double &target_mass, const double &del_E, 
    //     double (&incident_velocity)[3], double (&target_velocity)[3]);
    // inline void triple_product_isotropic(const double &primary_mass, const double &ion_mass, const double &target_mass, const double &del_E, 
    //     double (&incident_velocity)[3], double (&target_velocity)[3], double (&third_velocity)[3]);
    // void initialize_data_files(const std::string& dir_name, std::vector<Particle>& particle_list, std::vector<Target_Particle>& target_particle_list) const;
    // void diag_write(const std::string& dir_name, std::vector<Particle>& particle_list, std::vector<Target_Particle>& target_particle_list, const double& time_diff, bool average_bool = false);
    // void gather_mpi();
};

std::vector<null_collider> read_null_collision_inputs(const std::string& directory, const std::vector<charged_particle> &particle_list, const std::vector<target_particle> &target_particle_list);

