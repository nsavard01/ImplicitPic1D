#ifndef NULL_COLLISION_H
#define NULL_COLLISION_H

#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include "pcg_rng.h"
#include "particle.h"
#include "target_particle.h"


class Null_Collision {
    // collision primary particle with background gas
public:
    int number_collisions, length_arrays, primary_idx, target_idx;
    std::vector<std::vector<double>> sigma_array;
    std::vector<double> energy_array, energy_threshold, total_incident_energy, total_energy_loss;
    double sigma_v_max, min_energy, max_energy, mass_sum, reduced_mass, reduced_mass_triple;
    std::vector<std::vector<int>> products_indx;
    std::vector<int> collision_type;
    std::vector<size_t> total_amount_collisions;
    size_t total_amount_collidable_particles;
    Null_Collision(int number_collisions, int length_arrays, std::vector<int> reactant_idx, 
    const std::vector<std::vector<double>> &sigma_array, 
    const std::vector<double> &energy_array, const std::vector<double> &energy_threshold,
    const std::vector<int> &collision_type, const std::vector<std::vector<int>> &products_indx, double mass_inputs[3]);
    void generate_null_collisions(std::vector<Particle> &particle_list, std::vector<Target_Particle> &target_particle_list, double time_step);
    inline void double_product_isotropic(const double &primary_mass, const double &target_mass, const double &del_E, 
        double (&incident_velocity)[3], double (&target_velocity)[3]);
    inline void triple_product_isotropic(const double &primary_mass, const double &ion_mass, const double &target_mass, const double &del_E, 
        double (&incident_velocity)[3], double (&target_velocity)[3], double (&third_velocity)[3]);
};

std::vector<Null_Collision> read_null_collision_inputs(const std::string& filename, const std::vector<Particle> &particle_list, const std::vector<Target_Particle> &target_particle_list);

#endif
