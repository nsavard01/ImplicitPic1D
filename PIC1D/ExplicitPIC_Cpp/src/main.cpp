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
    double KE_before = 0.0;
    double KE_after = 0.0;
    double coll_loss = 0.0;
    
    global_inputs::read_global_inputs("../InputData/InitialConditions.inp");
    initialize_pcg(false);
    Domain world("../InputData/Geometry.inp");
    std::vector<Particle> particle_list = read_particle_inputs("../InputData/ParticleTypes.inp", world);
    for (i=0;i<global_inputs::number_charged_particles;i++){
        KE_before += particle_list[i].get_KE_total();
        std::cout << "Number " << particle_list[i].name << " is: " << std::endl;
        for (int iter = 0; iter < omp_get_max_threads(); iter++) {
            std::cout << " " << particle_list[i].number_particles[iter];
        }
        std::cout << std::endl;
        std::cout << "Ave KE " << particle_list[i].get_KE_ave() << std::endl;
    }
    std::cout << "KE_before " << KE_before << std::endl;
    std::vector<Target_Particle> target_particle_list = read_target_particle_inputs("../InputData/ParticleTypes.inp");
    std::vector<Null_Collision> binary_collision_list = read_null_collision_inputs("../InputData/collision.inp", particle_list, target_particle_list);
    for (i=0;i<global_inputs::number_binary_collisions;i++){
        binary_collision_list[i].generate_null_collisions(particle_list, target_particle_list, global_inputs::time_step);
        std::cout << "Binary Collision number " << i << std::endl;
        for (int iter=0; iter < binary_collision_list[i].number_collisions; iter++) {
            std::cout << "Coll number " << iter << " with energy loss " << 0.5 * binary_collision_list[i].total_energy_loss[iter] * particle_list[binary_collision_list[i].primary_idx].weight << std::endl;
            coll_loss += 0.5 * binary_collision_list[i].total_energy_loss[iter] * particle_list[binary_collision_list[i].primary_idx].weight;
        }
    }
    for (i=0;i<global_inputs::number_charged_particles;i++){
        KE_after += particle_list[i].get_KE_total();
        std::cout << "Number " << particle_list[i].name << " is: " << std::endl;
        for (int iter = 0; iter < omp_get_max_threads(); iter++) {
            std::cout << " " << particle_list[i].number_particles[iter];
        }
        std::cout << std::endl;
        std::cout << "Ave KE " << particle_list[i].get_KE_ave() << std::endl;
    }
    std::cout << "Coll_loss " << (coll_loss + (KE_after - KE_before))/(KE_after-KE_before) << std::endl;
    // #pragma omp parallel private(i,val) reduction(+:sum,var)
    // {
    //     for (i=0; i < 10000; i++){
    //         val = pcg32_random_r();
    //         sum += val;
    //         var += (val - 0.5) * (val - 0.5);
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
