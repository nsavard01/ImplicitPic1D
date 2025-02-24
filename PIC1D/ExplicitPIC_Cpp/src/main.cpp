#include <iostream>
#include <stdio.h>
#include "pcg_rng.h"
#include <omp.h>
#include <iomanip>
#include <cmath>
#include "Constants.h"
#include "global_inputs.h"
#include "domain.h"
#include "particle.h"
#include "target_particle.h"
#include "null_collision.h"

int main() {
    int thread_id, i;
    int number_stats = 10;
    double sum = 0.0;
    double var = 0.0;
    double val;
    
    global_inputs::read_global_inputs("../InputData/InitialConditions.inp");
    initialize_pcg(false);
    Domain world("../InputData/Geometry.inp");
    std::vector<Particle> particle_list = read_particle_inputs("../InputData/ParticleTypes.inp", world);
    std::vector<Target_Particle> target_particle_list = read_target_particle_inputs("../InputData/ParticleTypes.inp");
    std::vector<Null_Collision> binary_collision_list = read_null_collision_inputs("../InputData/collision.inp", particle_list, target_particle_list);
    
    #pragma omp parallel private(i,val) reduction(+:sum,var)
    {
        for (i=0; i < 10000; i++){
            val = pcg32_random_r();
            sum += val;
            var += (val - 0.5) * (val - 0.5);
        }
        #pragma omp critical
        {
            std::cout <<"local sum thread " << omp_get_thread_num() << " is " << sum / (10000.0) << std::endl;
            std::cout <<"local std thread " << omp_get_thread_num() << " is " << std::sqrt(var/10000.0) << std::endl;
        }
    }
    std::cout << "Total average "<< sum / (10000.0) /omp_get_max_threads() << std::endl;
    std::cout << "Total std "<< std::sqrt(var / (10000.0) /omp_get_max_threads()) << std::endl;
    return 0;
}
