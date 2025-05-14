
#include <vector>
#include <omp.h>
#include <cmath>
#include "ES_solvers/ES_solver.hpp"
#include "domain/domain.hpp"
#include "domain/uniform_domain.hpp"
#include "domain/non_uniform_domain.hpp"
#include "ES_solvers/ES_solver_EC.hpp"
#include "ES_solvers/ES_solver_MC.hpp"
#include "globals/mpi_vars.hpp"
#include <iostream>
#include <fstream>
#include <sstream>


void ES_solver::set_phi(double left_voltage, double right_voltage, double RF_frequency, int left_boundary, int right_boundary) {
    this->left_voltage = left_voltage;
    this->right_voltage = right_voltage;
    this->RF_rad_frequency = RF_frequency * 2.0 * M_PI; // Convert to radians
    this->RF_indx = -1; // Initialize RF index to -1
    if (left_boundary != 1 && left_boundary != 4) {
        this->left_voltage = 0.0; // Set left boundary voltage
    } 
    if (right_boundary != 1 && right_boundary != 4) {
        this->right_voltage = 0.0; // Set left boundary voltage to zero
    } 
    if (left_boundary == 4) {
        this->RF_half_amplitude = this->left_voltage; // Set RF half amplitude for left boundary
        this->left_voltage = 0.0;
        this->RF_indx = 0; // Set RF index for left boundary
    } else if (right_boundary == 4) {
        this->RF_half_amplitude = this->right_voltage; // Set RF half amplitude for right boundary
        this->right_voltage = 0.0;
        this->RF_indx = this->phi.size()-1; // Set RF index for right boundary
    }
    this->phi[0] = this->left_voltage; // Set left boundary voltage in phi vector
    this->phi[this->phi.size()-1] = this->right_voltage; // Set right boundary voltage in phi vector
}

void ES_solver::deposit_charge_density(std::vector<charged_particle>& particle_list) {
    // Loop over all particles and deposit charge density
    
    int total_thread_count = omp_get_max_threads();
    int total_rho_size = this->rho.size();
    int num_particles = 1; //particle_list.size();
    #pragma omp parallel
    {   
        int thread_id = omp_get_thread_num();
        for (int i = 0; i < num_particles; ++i) {
            charged_particle& particle = particle_list[i];
            particle.clear_work_space(thread_id); // Clear work space for each thread
            particle.deposit_particles_linear(thread_id); // Deposit charge density for each particle
            double q_time_wp = particle.get_q_times_wp(); // Get charge density for each particle
            std::vector<double>& work_space = particle.get_work_space(thread_id);
            for (int j = 0; j < total_rho_size; ++j) {
                work_space[j] *= q_time_wp; // Accumulate charge density from all particles
            }
        }
    }
    // Accumulate charge density from all threads
    std::fill(this->rho.begin(), this->rho.end(), 0.0);
    for (int i = 0; i < num_particles; ++i) {
        charged_particle& particle = particle_list[i];
        for (int j = 0; j < total_thread_count; ++j) {
            std::vector<double>& work_space = particle.get_work_space(j);
            for (int k = 0; k < total_rho_size; ++k) {
                this->rho[k] += work_space[k]; // Accumulate charge density from all threads
            }
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, this->rho.data(), total_rho_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); // Synchronize charge density across all processes
    if (mpi_vars::mpi_rank == 0) {
        for (int i = 0; i < total_rho_size; ++i) {
            std::cout << "rho[" << i << "] = " << this->rho[i] << std::endl; // Print charge density for debugging
        }
    }
    
}

std::unique_ptr<ES_solver> read_voltage_inputs(const std::string& filename, int scheme_type, const domain& world){
    double left_voltage, right_voltage, RF_frequency;
    if (mpi_vars::mpi_rank == 0) {
        std::cout << " "  << std::endl;
        std::cout << "Reading electric potential inputs: "  << std::endl;
        std::cout << "-------------------------- "  << std::endl;
        std::string line;
        std::ifstream file(filename);
        if (!file) {
            std::cerr << "Error: Unable to open file " << filename << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        
        
        std::getline(file, line);
        std::getline(file, line);
        std::getline(file, line);
        std::getline(file, line);
        std::getline(file, line);
        std::istringstream iss(line);
        iss >> left_voltage >> right_voltage;
        iss.clear();
        std::getline(file, line);
        iss.str(line);
        iss >> RF_frequency;
        iss.clear();
        file.close();
    }
    MPI_Bcast(&left_voltage, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&right_voltage, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&RF_frequency, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    std::unique_ptr<ES_solver> es_solver;
    if (scheme_type == 0) {
        es_solver = std::make_unique<ES_solver_MC>(world);
    } else if (scheme_type == 1) {
        es_solver = std::make_unique<ES_solver_EC>(world);
    } else {
        throw std::invalid_argument("Invalid scheme type for ES solver.");
    }
    es_solver->set_phi(left_voltage, right_voltage, RF_frequency, world.get_left_boundary_condition(), world.get_right_boundary_condition());
    es_solver->print_out();
    return es_solver;

}



