#include "global_inputs.h"
#include <fstream>
#include <cmath>
#include <sstream>
#include <omp.h>
#include "Constants.h"


namespace global_inputs {
    int number_omp_threads, number_nodes, number_cells, number_charged_particles, number_diagnostic_steps, number_target_particles, number_binary_collisions;
    double temp_electrons, temp_ions, initial_density;
    double time_step, averaging_time, simulation_time, start_simulation_time;
    std::string save_folder, save_filename;

    void read_global_inputs(const std::string& filename){
        std::ifstream file(filename);
        if (!file) {
            std::cerr << "Error: Unable to open file " << filename << std::endl;
            return;
        }
        double temp_real;
        std::string line;
        if (Constants::mpi_rank==0) {
            std::cout << "Reading global inputs: "  << std::endl;
            std::cout << "-------------------------- "  << std::endl;
        }
        std::getline(file, line);
        std::istringstream iss(line);
        iss >> number_omp_threads;
        iss.clear();

        std::getline(file, line);
        iss.str(line);
        iss >> number_cells;
        iss.clear();
        number_nodes = number_cells + 1;

        std::getline(file, line);
        iss.str(line);
        iss >> simulation_time >> start_simulation_time;
        iss.clear();

        std::getline(file, line);
        iss.str(line);
        iss >> initial_density;
        iss.clear();

        std::getline(file, line);
        iss.str(line);
        iss >> temp_electrons;
        iss.clear();

        std::getline(file, line);
        iss.str(line);
        iss >> temp_ions;
        iss.clear();

        std::getline(file, line);
        iss.str(line);
        iss >> number_diagnostic_steps;
        iss.clear();

        std::getline(file, line);
        iss.str(line);
        iss >> temp_real >> time_step;
        iss.clear();

        temp_real = temp_real / std::sqrt(initial_density * Constants::elementary_charge * Constants::elementary_charge / Constants::electron_mass / Constants::epsilon_0);
        time_step = std::min(temp_real, time_step);

        std::getline(file, line);
        iss.str(line);
        iss >> averaging_time;
        iss.clear();

        std::getline(file, line);
        iss.str(line);
        iss >> save_folder;
        iss.clear();

        std::getline(file, line);
        iss.str(line);
        iss >> save_filename;
        iss.clear();

        
        
        file.close();
        omp_set_num_threads(number_omp_threads);
        if (Constants::mpi_rank == 0) {
            std::cout << "Number of mpi ranks: " << Constants::mpi_size << std::endl;
            std::cout << "Number openmp threads: " << omp_get_max_threads() << std::endl;
            std::cout << "Number of nodes: " << global_inputs::number_nodes << std::endl;
            std::cout << "Simulation time: " << simulation_time << std::endl;
            std::cout << "Initial density: " << initial_density << std::endl;
            std::cout << "T_e: " << temp_electrons << std::endl;
            std::cout << "T_i: " << temp_ions << std::endl;
            std::cout << "Number diagnostic steps: " << number_diagnostic_steps << std::endl;
            std::cout << "Time step: " << time_step << std::endl;
            std::cout << "Averaging time: " << averaging_time << std::endl;
            std::cout << "Relative save folder: " << save_folder << std::endl;
            std::cout << "Save filename: " << save_filename << std::endl;
            std::cout << "-----------------------"  << std::endl;
        }

    }


}
