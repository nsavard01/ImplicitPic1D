
#pragma once

#include "globals/constants.hpp"
#include "globals/mpi_vars.hpp"
#include "domain/domain.hpp"
#include "domain/uniform_domain.hpp"
#include "domain/non_uniform_domain.hpp"
#include "rand_gen/pcg_rng.hpp"
#include "particles/charged_particle.hpp"
#include "particles/target_particle.hpp"
#include "solvers/poisson_solver_1D_tridiag.hpp"
#include "ES_solvers/ES_solver_MC.hpp"
#include "ES_solvers/ES_solver_EC.hpp"
#include "collisions/null_collider.hpp"
#include <stdio.h>
#include <iostream>
#include <mpi.h>
#include <omp.h>
#include <fstream>
#include <sstream>

class simulation {
    
public:
    // Everything needed to simulate particle in cell
    double del_t, simulation_time, averaging_time;
    double total_particle_momentum[3];
    double total_particle_KE[3];
    double total_field_energy;
    double inv_plasma_freq_fraction;
    bool restarted_simulation;
    std::string save_file_folder, save_file_path;
    int number_omp_threads, scheme_type, number_diagnostics;
    std::unique_ptr<domain> world;
    std::vector<charged_particle> charged_particle_list;
    std::vector<target_particle> target_particle_list;
    std::vector<null_collider> null_collider_list;
    std::unique_ptr<ES_solver> field_solver;
    simulation();
    void setup();
};




