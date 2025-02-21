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
    int number_collisions, length_arrays;
    std::vector<std::vector<double>> sigma_array;
    std::vector<double> energy_array, energy_threshold, total_incident_energy, total_energy_loss;
    double nu_max, min_energy, max_energy;
    std::vector<int> collision_type, reactant_indx, number_products;
    std::vector<std::vector<int>> products_indx;
    std::vector<size_t> total_amount_collisions;
};

void read_null_collision_inputs(const std::string& filename, const std::vector<Particle> &particle_list, const std::vector<Target_Particle> &target_particle_list);

#endif
