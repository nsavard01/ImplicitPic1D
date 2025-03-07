#ifndef TARGET_PARTICLE_H
#define TARGET_PARTICLE_H

#include <cmath>
#include "Constants.h"
#include "global_inputs.h"
#include "pcg_rng.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <array>

class Target_Particle {
public:
    double mass, density, temperature, v_therm, accum_energy_change;
    std::string name;
    // Accessor for 3D indexing
    Target_Particle(double mass_in, double temp_in, double density_in, std::string name_in);
    inline void generate_maxwellian_velocity(double (&v_array)[3]) const {
        // Pass a reference to array and change in place
        double R_1, R_2, R_3, R_4;
        R_1 = pcg32_random_r();
        R_2 = pcg32_random_r();
        R_3 = pcg32_random_r();
        R_4 = pcg32_random_r();
        v_array[0] = this->v_therm * std::sqrt(-2.0 * std::log(R_1)) * std::cos(2.0 * M_PI * R_2);
        v_array[1] = this->v_therm * std::sqrt(-2.0 * std::log(R_1)) * std::sin(2.0 * M_PI * R_2);
        v_array[2] = this->v_therm * std::sqrt(-2.0 * std::log(R_3)) * std::cos(2.0 * M_PI * R_4);

    }
};
std::vector<Target_Particle> read_target_particle_inputs(const std::string& filename);

#endif // PARTICLE_H
