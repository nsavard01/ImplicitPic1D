#ifndef GLOBAL_INPUTS_H
#define GLOBAL_INPUTS_H
#include <string>
#include <iostream>

namespace global_inputs {
    extern int number_omp_threads,number_nodes, number_cells, number_charged_particles, number_diagnostic_steps, number_target_particles, number_primary_colliders;
    extern double temp_electrons, temp_ions, initial_density;
    extern double time_step, averaging_time, simulation_time, start_simulation_time;
    extern std::string save_folder, save_filename;

    void read_global_inputs(const std::string& filename);
}


#endif // PCG_RNG_H