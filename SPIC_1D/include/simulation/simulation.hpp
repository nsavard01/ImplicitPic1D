
#pragma once

#include "globals/constants.hpp"
#include "globals/mpi_vars.hpp"
#include "domain/domain.hpp"
#include "domain/uniform_domain.hpp"
#include "domain/non_uniform_domain.hpp"
#include "rand_gen/pcg_rng.hpp"
#include "particles/charged_particle.hpp"
#include "solvers/poisson_solver_1D_tridiag.hpp"
#include "ES_solvers/ES_solver_MC.hpp"
#include "ES_solvers/ES_solver_EC.hpp"
#include <stdio.h>
#include <iostream>
#include <mpi.h>
#include <omp.h>
#include <fstream>
#include <sstream>

class simulation {
    
public:
    double del_t;
    int number_omp_threads, scheme_type;
    std::unique_ptr<domain> world;
    std::vector<charged_particle> charged_particle_list;
    std::unique_ptr<ES_solver> field_solver;
    simulation();
};




