#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include <cstddef>
#include <iostream>
#include <string>

class Particle {
public:
    double mass, charge, weight, q_over_m, q_times_wp;
    int accum_wall_loss[2];
    std::string name;
    double accum_energy_loss[2];
    size_t final_idx;
    std::vector<size_t> number_particles;
    std::vector<double> phase_space, wall_loss, energy_loss, momentum_loss, densities, work_space;  // Flattened 3D array

    // Accessor for 3D indexing
    inline double& phase_space_at(size_t phase_idx, size_t part_idx, size_t thread_idx);
};
void read_particle_inputs(const std::string& filename);

#endif // PARTICLE_H
