
#include <vector>
#include <omp.h>
#include <cmath>
#include "ES_solvers/ES_solver.hpp"
#include "solvers/poisson_solver_1D_tridiag.hpp"
#include "domain/domain.hpp"
#include "domain/uniform_domain.hpp"
#include "domain/non_uniform_domain.hpp"
#include "globals/mpi_vars.hpp"
#include "globals/constants.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include "globals/write_functions.hpp"
#include <iomanip>



ES_solver::ES_solver(double left_voltage, double right_voltage, const domain& world){
    this->phi.resize(world.number_nodes, 0.0);
    this->rho.resize(world.number_nodes, 0.0);
    this->E_field.resize(world.number_cells, 0.0);
    this->left_voltage = left_voltage;
    this->right_voltage = right_voltage;
    this->poisson_solver = std::make_unique<poisson_solver_1D_tridiag>(world);
    if (world.left_boundary_condition != 1 && world.left_boundary_condition != 4) {
        this->left_voltage = 0.0; // Set left boundary voltage
    } 
    if (world.right_boundary_condition != 1 && world.right_boundary_condition != 4) {
        this->right_voltage = 0.0; // Set left boundary voltage to zero
    } 
    this->phi[0] = this->left_voltage; // Set left boundary voltage in phi vector
    this->phi[this->phi.size()-1] = this->right_voltage; // Set right boundary voltage in phi vector
}


void ES_solver::print_out() {
    if (mpi_vars::mpi_rank == 0) {
        std::cout << "ES_solver: " << std::endl;
        std::cout << "-------------------------- " << std::endl;
        std::cout << "Number of phi nodes: " << this->phi.size() << std::endl;
        std::cout << "Number of field nodes: " << this->E_field.size() << std::endl;
        std::cout << "Left voltage: " << this->left_voltage << std::endl;
        std::cout << "Right voltage: " << this->right_voltage << std::endl;
        std::cout << "-------------------------- " << std::endl;
    }
}


// void ES_solver::solve_potential(double current_time, const domain& world) {
//     // generate the right-hand side of the Poisson equation
//     double inv_epsilon_0 = 1.0 / constants::epsilon_0; // Inverse of permittivity
//     int left_boundary = world.left_boundary_condition; // Get left boundary condition
//     int right_boundary = world.right_boundary_condition; // Get right boundary condition
//     int number_unknowns = this->poisson_solver->number_unknowns; // Number of unknowns in the system
//     if (left_boundary == 2) {
//         this->phi[0] = -this->rho[0] * inv_epsilon_0; // Change boundary phi
//     } else if (left_boundary == 4) {
//         this->phi[0] = this->RF_half_amplitude * std::sin(this->RF_rad_frequency * current_time); // Change boundary phi
//     } 

//     if (right_boundary == 2) {
//         this->phi[number_unknowns-1] = -this->rho[number_unknowns-1] * inv_epsilon_0; // Change boundary phi
//     } else if (right_boundary == 4) {
//         this->phi[number_unknowns-1] = this->RF_half_amplitude * std::sin(this->RF_rad_frequency * current_time); // Change boundary phi
//     } 

//     for (int i = 1; i < number_unknowns-1; ++i) {
//         this->phi[i] = -this->rho[i] * inv_epsilon_0; // Set right-hand side of the Poisson equation
//     }

//     this->poisson_solver->solve(this->phi, this->phi); // replace phi with solution
    
// }

// void ES_solver::make_EField(const domain& world) {
//     // Calculate the electric field from the potential
//     int number_nodes = world.number_nodes; // Number of cells in the domain
//     int number_cells = world.number_cells;
//     double inv_dx = 1.0/world.min_dx; // Cell size
//     int left_boundary = world.left_boundary_condition; // Get left boundary condition
//     int right_boundary = world.right_boundary_condition; // Get right boundary condition
//     for (int i = 1; i < number_nodes-1; ++i) {
//         this->E_field[i] = 0.5 * (this->phi[i-1] - this->phi[i+1]) * inv_dx; // Electric field calculation
//     }
//     if (left_boundary == 1 || left_boundary == 4) {
//         // First order at boundary consistent with rho = 0
//         this->E_field[0] = (this->phi[0] - this->phi[1])*inv_dx; // Electric field at left boundary
//     } else if (left_boundary == 3){
//         this->E_field[0] = 0.5 * (this->phi[number_cells-1] - this->phi[1]) * inv_dx; 
//         this->E_field[number_cells] = this->E_field[0]; 
//     }
//     if (right_boundary == 1 || right_boundary == 4) {
//         this->E_field[number_cells] = (this->phi[number_cells-1] - this->phi[number_cells]) * inv_dx; 
//     } 
     
// }



std::unique_ptr<ES_solver> read_voltage_inputs(const std::string& filename, const domain& world){
    double left_voltage, right_voltage;
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
        file.close();
    }
    MPI_Bcast(&left_voltage, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&right_voltage, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    std::unique_ptr<ES_solver> es_solver;
    es_solver = std::make_unique<ES_solver>(left_voltage, right_voltage, world);
    es_solver->print_out();
    
    return es_solver;

}



