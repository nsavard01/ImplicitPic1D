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
#include "potential_solver.h"

int main() {
    int i;
    
    global_inputs::read_global_inputs("../InputData/InitialConditions.inp");
    initialize_pcg(false);
    Domain world("../InputData/Geometry.inp");
    std::vector<Particle> particle_list = read_particle_inputs("../InputData/ParticleTypes.inp", world);
    Potential_Solver solver("../InputData/Geometry.inp", world);
    std::vector<Target_Particle> target_particle_list = read_target_particle_inputs("../InputData/ParticleTypes.inp");
    std::vector<Null_Collision> binary_collision_list = read_null_collision_inputs("../InputData/collision.inp", particle_list, target_particle_list);
    solver.deposit_rho(particle_list, world);
    solver.solve_potential_tridiag(world, 0.0);
    solver.make_EField(world);
    solver.initial_v_rewind(particle_list, global_inputs::time_step);
    double P_before = 0.0;
    for (i=0;i<global_inputs::number_charged_particles;i++){
        P_before += particle_list[i].get_momentum_total();
        std::cout << "KE " << particle_list[i].get_KE_ave() << std::endl;
    }
    std::cout << "P_before " << P_before << std::endl;
    solver.move_particles(particle_list, world, global_inputs::time_step);
    double P_after = 0.0;
    for (i=0;i<global_inputs::number_charged_particles;i++){
        P_after += particle_list[i].get_momentum_total();
        std::cout << "KE " << particle_list[i].get_KE_ave() << std::endl;
    }
    std::cout << "P_after " << P_after << std::endl;

    // for (i=0;i<global_inputs::number_binary_collisions;i++){
    //     binary_collision_list[i].generate_null_collisions(particle_list, target_particle_list, global_inputs::time_step);
    // }
    
    // #pragma omp parallel private(i,val) reduction(+:sum,var)
    // {
    //     for (i=0; i < 10000; i++){
    //         val = pcg32_random_r();
    //         sum += val;
    //         var += (val - 0.5) * (val - 0.5);./main
    //     }
    //     #pragma omp critical
    //     {
    //         std::cout <<"local sum thread " << omp_get_thread_num() << " is " << sum / (10000.0) << std::endl;
    //         std::cout <<"local std thread " << omp_get_thread_num() << " is " << std::sqrt(var/10000.0) << std::endl;
    //     }
    // }
    // std::cout << "Total average "<< sum / (10000.0) /omp_get_max_threads() << std::endl;
    // std::cout << "Total std "<< std::sqrt(var / (10000.0) /omp_get_max_threads()) << std::endl;
    return 0;
}
