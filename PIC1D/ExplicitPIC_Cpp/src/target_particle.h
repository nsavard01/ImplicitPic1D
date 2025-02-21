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
    void generate_maxwellian_velocity(double (&v_array)[3]) const;
};
std::vector<Target_Particle> read_target_particle_inputs(const std::string& filename);

#endif // PARTICLE_H
