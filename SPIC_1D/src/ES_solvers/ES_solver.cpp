
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
#include "globals/constants.hpp"
#include <iostream>
#include <fstream>
#include <sstream>


void ES_solver::set_phi(double left_voltage, double right_voltage, double RF_frequency, int left_boundary, int right_boundary) {
    this->left_voltage = left_voltage;
    this->right_voltage = right_voltage;
    this->RF_rad_frequency = 0.0; // Convert to radians
    if (left_boundary != 1 && left_boundary != 4) {
        this->left_voltage = 0.0; // Set left boundary voltage
    } 
    if (right_boundary != 1 && right_boundary != 4) {
        this->right_voltage = 0.0; // Set left boundary voltage to zero
    } 
    if (left_boundary == 4) {
        this->RF_half_amplitude = this->left_voltage; // Set RF half amplitude for left boundary
        this->left_voltage = 0.0;
    } else if (right_boundary == 4) {
        this->RF_half_amplitude = this->right_voltage; // Set RF half amplitude for right boundary
        this->right_voltage = 0.0;
    } 
    this->phi[0] = this->left_voltage; // Set left boundary voltage in phi vector
    this->phi[this->phi.size()-1] = this->right_voltage; // Set right boundary voltage in phi vector
}

void ES_solver::deposit_charge_density(std::vector<charged_particle>& particle_list) {
    // Loop over all particles and deposit charge density
    
    int total_thread_count = omp_get_max_threads();
    int total_rho_size = this->rho.size();
    int num_particles = particle_list.size();
    #pragma omp parallel
    {   
        int thread_id = omp_get_thread_num();
        // use xi_sorted as workspace, since N_p >> number_nodes
        std::vector<double>& part_work_space = charged_particle::xi_sorted[thread_id];
        // local work_space to accumulate over each particle
        std::vector<double>& local_work_space = this->work_space[thread_id];
        std::fill(local_work_space.begin(), local_work_space.end(), 0.0);
        for (int i = 0; i < num_particles; ++i) {
            charged_particle& particle = particle_list[i];
            // Set work space to 0
            std::fill(part_work_space.begin(), part_work_space.begin() + total_rho_size, 0.0);
            particle.deposit_particles_linear(thread_id, part_work_space); // Deposit charge density for each particle
            double q_time_wp = particle.q_times_wp; // Get charge density for each particle
            for (int j = 0; j < total_rho_size; ++j) {
                local_work_space[j] += part_work_space[j] * q_time_wp; // Accumulate charge density from all particles
            }
        }
    }
    // Accumulate charge density from all threads
    std::fill(this->rho.begin(), this->rho.end(), 0.0);
    for (int j = 0; j < total_thread_count; ++j) {
        std::vector<double>& work_space = this->work_space[j];
        for (int k = 0; k < total_rho_size; ++k) {
            this->rho[k] += work_space[k]; // Accumulate charge density from all threads
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, this->rho.data(), total_rho_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); // Synchronize charge density across all processes
    
}

void ES_solver::solve_potential(double current_time, const domain& world) {
    // generate the right-hand side of the Poisson equation
    double inv_epsilon_0 = 1.0 / constants::epsilon_0; // Inverse of permittivity
    int left_boundary = world.left_boundary_condition; // Get left boundary condition
    int right_boundary = world.right_boundary_condition; // Get right boundary condition
    int number_unknowns = this->poisson_solver->number_unknowns; // Number of unknowns in the system
    if (left_boundary == 2) {
        this->phi[0] = -this->rho[0] * inv_epsilon_0; // Change boundary phi
    } else if (left_boundary == 4) {
        this->phi[0] = this->RF_half_amplitude * std::sin(this->RF_rad_frequency * current_time); // Change boundary phi
    } 

    if (right_boundary == 2) {
        this->phi[number_unknowns-1] = -this->rho[number_unknowns-1] * inv_epsilon_0; // Change boundary phi
    } else if (left_boundary == 4) {
        this->phi[number_unknowns-1] = this->RF_half_amplitude * std::sin(this->RF_rad_frequency * current_time); // Change boundary phi
    } 

    for (int i = 1; i < number_unknowns-1; ++i) {
        this->phi[i] = -this->rho[i] * inv_epsilon_0; // Set right-hand side of the Poisson equation
    }

    this->poisson_solver->solve(this->phi, this->phi); // replace phi with solution

    



}

void ES_solver::make_EField(const domain& world) {
    // Calculate the electric field from the potential
    int number_nodes = world.number_nodes; // Number of cells in the domain
    double inv_dx = 1.0/world.min_dx; // Cell size
    int left_boundary = world.left_boundary_condition; // Get left boundary condition
    int right_boundary = world.right_boundary_condition; // Get right boundary condition
    for (int i = 1; i < number_nodes-1; ++i) {
        this->E_field[i] = 0.5 * (this->phi[i-1] - this->phi[i+1]) * inv_dx; // Electric field calculation
    }
    if (left_boundary == 1 || left_boundary == 4) {
        // First order at boundary consistent with rho = 0
        this->E_field[0] = (this->phi[0] - this->phi[1])*inv_dx; // Electric field at left boundary
    } else if (left_boundary == 3){
        this->E_field[0] = 0.5 * (this->phi[number_nodes-2] - this->phi[1]) * inv_dx; 
        this->E_field[number_nodes-1] = this->E_field[0]; 
    }
    if (right_boundary == 1 || right_boundary == 4) {
        this->E_field[number_nodes-1] = (this->phi[number_nodes-2] - this->phi[number_nodes-1]) * inv_dx; 
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



