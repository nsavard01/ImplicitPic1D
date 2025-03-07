#include <iostream>
#include <stdio.h>
#include "pcg_rng.h"
#include <omp.h>
#include <mpi.h>
#include <iomanip>
#include <cmath>
#include "Constants.h"
#include "basic_tools.h"
#include "global_inputs.h"
#include "domain.h"
#include "particle.h"
#include "target_particle.h"
#include "null_collision.h"
#include "potential_solver.h"
#include "simulation.h"





int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    
    int i;
    
    MPI_Comm_size(MPI_COMM_WORLD, &Constants::mpi_size);  // Get total processes
    MPI_Comm_rank(MPI_COMM_WORLD, &Constants::mpi_rank);  // Get current rank
    
    Domain world;
    std::vector<Particle> particle_list;
    Potential_Solver solver;
    std::vector<Target_Particle> target_particle_list;
    std::vector<Null_Collision> binary_collision_list;
    global_inputs::read_global_inputs("../InputData/InitialConditions.inp");

    initialize_pcg(false);
    world.read_from_file("../InputData/Geometry.inp");
    particle_list = read_particle_inputs("../InputData/ParticleTypes.inp", world);
    solver.read_from_file("../InputData/Geometry.inp", world);
    target_particle_list = read_target_particle_inputs("../InputData/ParticleTypes.inp");
    binary_collision_list = read_null_collision_inputs("../InputData/collision.inp", particle_list, target_particle_list);
    
    solver.deposit_rho(particle_list, world);
    solver.solve_potential_tridiag(world, 0.0);
    solver.make_EField(world);
    solver.initial_v_rewind(particle_list, global_inputs::time_step);
    

    Simulation simulator(global_inputs::save_folder + global_inputs::save_filename);

    simulator.initialize_data_files(solver, particle_list, target_particle_list, binary_collision_list, world);
    simulator.reset_diag(particle_list, binary_collision_list);
    simulator.run(solver, particle_list, target_particle_list, binary_collision_list, world);
    simulator.averaging(solver, particle_list, target_particle_list, binary_collision_list, world);

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

    


    MPI_Finalize();
    return 0;
}
