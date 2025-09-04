
#pragma once
#include <vector>
#include <cmath>
#include "domain/domain.hpp"
#include "domain/uniform_domain.hpp"
#include "domain/non_uniform_domain.hpp"
#include "globals/mpi_vars.hpp"
#include "particles/charged_particle.hpp"
#include "particles/target_particle.hpp"
#include "collisions/null_collider.hpp"
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

    double max_cpu_time, min_cpu_time;
    std::vector<double> thread_time; // [thread]
    std::vector<std::vector<double>> integration_time;// [thread, charged_particle_number]
    std::vector<std::vector<size_t>> number_cell_crossings, number_time_steps;// [thread, charged_particle_number]
    std::vector<double> total_integration_time;// [charged_particle_number]
    std::vector<size_t> total_number_cell_crossings, total_number_time_steps;// [charged_particle_number]
    
    charged_particle_pusher();
    void set_diagnostic_vectors(const std::vector<charged_particle>& particle_list);
    void gather_mpi();
    void print_out(const std::vector<charged_particle>& particle_list) const;
    void push_particle_trajectories_non_uniform(const int thread_id, std::vector<charged_particle>& particle_list, std::vector<null_collider>& null_collider_list, const std::vector<target_particle>& target_particle_list, const non_uniform_domain& world, const std::vector<double>& E_field);
    void push_particle_trajectories_uniform(const int thread_id, std::vector<charged_particle>& particle_list, std::vector<null_collider>& null_collider_list, const std::vector<target_particle>& target_particle_list, const uniform_domain& world, const std::vector<double>& E_field);
};

