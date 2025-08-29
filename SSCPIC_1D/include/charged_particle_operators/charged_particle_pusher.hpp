
#pragma once
#include <vector>
#include <cmath>
#include "domain/domain.hpp"
#include "domain/uniform_domain.hpp"
#include "domain/non_uniform_domain.hpp"
#include "globals/mpi_vars.hpp"
#include "particles/charged_particle.hpp"
#include "particles/target_particle.hpp"
#include <cmath>
#include "globals/constants.hpp"
#include <omp.h>
#include <fstream>
#include <sstream>
#include <iomanip>

class charged_particle_pusher {

// general charged particle operator base
// input charged particle list, thread_id, and domain, and apply with run function using internal variables
public:

    
    charged_particle_pusher();
    void push_particle_trajectories_non_uniform(const int thread_id, std::vector<charged_particle>& particle_list, const std::vector<target_particle>& target_particle_list, const non_uniform_domain& world, const std::vector<double>& E_field);
    void push_particle_trajectories_uniform(const int thread_id, std::vector<charged_particle>& particle_list, const std::vector<target_particle>& target_particle_list, const uniform_domain& world, const std::vector<double>& E_field);
};

