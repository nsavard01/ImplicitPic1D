#include "potential_solver.h"
#include "global_inputs.h"
#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <fstream>
#include <sstream>
#include <omp.h>
#include "Constants.h"

Potential_Solver::Potential_Solver(const std::string& filename, const Domain& world){
    this->lower_tri.resize(global_inputs::number_nodes-1, 0.0);
    this->upper_tri.resize(global_inputs::number_nodes-1, 0.0);
    this->center_tri.resize(global_inputs::number_nodes, 0.0);
    this->phi.resize(global_inputs::number_nodes, 0.0);
    this->EField.resize(global_inputs::number_cells, 0.0);
    this->rho.resize(global_inputs::number_nodes, 0.0);
    this->work_space.resize(global_inputs::number_nodes, 0.0);
    this->RF_half_amplitude = 0.0;
    this->rho_const = 0.0;
    double left_voltage, right_voltage, RF_frequency;
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        return;
    }
    
    std::string line;
    std::cout << " "  << std::endl;
    std::cout << "Reading potential inputs: "  << std::endl;
    std::cout << "-------------------------- "  << std::endl;

    std::getline(file, line);  
    std::getline(file, line);
    std::getline(file, line);
    std::istringstream iss(line);
    iss >> left_voltage >> right_voltage;
    iss.clear();
    std::getline(file, line);
    iss.str(line);
    iss >> RF_frequency;
    file.close();

    if (world.boundary_conditions[0] == 1) {this->phi[0] = left_voltage;}
    if (world.boundary_conditions[0] == 4) {this->RF_half_amplitude = left_voltage;}
    if (world.boundary_conditions[global_inputs::number_cells] == 1) {this->phi[global_inputs::number_cells] = right_voltage;}
    if (world.boundary_conditions[global_inputs::number_cells] == 4) {
        if (this->RF_half_amplitude != 0) {
            std::cout << "Half amplitude voltage for RF already set, have two RF boundaries!" << std::endl;
            exit(EXIT_FAILURE);
        } else {
            this->RF_half_amplitude = right_voltage;
        }
    }
    if (world.boundary_conditions[0] == 3) {
        this->phi[0] = left_voltage;
        this->phi[global_inputs::number_cells] = left_voltage;
    }
    this->RF_rad_frequency = 2.0 * M_PI * RF_frequency; 
    std::cout << "Left voltage " << this->phi[0] << std::endl;
    std::cout << "Right voltage " << this->phi[global_inputs::number_cells] << std::endl;
    std::cout << "Rf frequency " << this->RF_rad_frequency/2.0/M_PI << std::endl;
    std::cout << "RF half amplitude " << this->RF_half_amplitude << std::endl;
    std::cout << "-------------------------- "  << std::endl;
    std::cout << " "  << std::endl;

    for (int i = 0; i < global_inputs::number_nodes; i++){
        switch (world.boundary_conditions[i]) {
            case 0:
                if (i < global_inputs::number_cells){
                    this->upper_tri[i] = 1.0/world.del_X;
                } 
                if (i > 0) {
                    this->lower_tri[i-1] = 1.0/world.del_X;
                }
                this->center_tri[i] = -2.0/world.del_X;
                break;
            case 1:
                this->center_tri[i] = 1.0;
                break;
            case 2:
                if (i == 0) {
                    this->upper_tri[i] = 1.0/world.del_X;
                    this->center_tri[i] = -1.0/world.del_X;
                } else if (i == global_inputs::number_cells) {
                    this->lower_tri[i-1] = 1.0/world.del_X;
                    this->center_tri[i] = -1.0/world.del_X; 
                } else {
                    std::cout << "Neumann boundary not on left or right index!" << std::endl;
                    exit(EXIT_FAILURE);
                }
                break;
            case 3:
                this->center_tri[i] = 1.0;
                break;
            case 4:
                this->center_tri[i] = 1.0;
                break;
        }
    }
}

void Potential_Solver::deposit_rho(std::vector<Particle> &particle_list, const Domain& world) {
    #pragma omp parallel
    {   
        for (int j = 0; j<global_inputs::number_charged_particles;j++){
            particle_list[j].interpolate_particles();
        }  

        #pragma omp barrier

        int thread_id = omp_get_thread_num();
        int left_thread_indx = world.thread_node_indx[thread_id][0];
        int right_thread_indx = world.thread_node_indx[thread_id][1];
        
        for (int j = 0; j<global_inputs::number_charged_particles;j++){
            for (int i_thread = 0; i_thread < global_inputs::number_omp_threads; i_thread++) {
                for (int i = left_thread_indx; i <= right_thread_indx; i++) {
                    this->rho[i] += particle_list[j].work_space[i_thread][i] * particle_list[j].q_times_wp;
                }
            }
        }
    }
}

void Potential_Solver::solve_potential_tridiag(const Domain& world, double time) {
    double temp;
    if (world.boundary_conditions[0] == 2) {
        this->phi[0] = -this->rho[0]/Constants::epsilon_0;
    } else if (world.boundary_conditions[0] == 4) {
        this->phi[0] = this->RF_half_amplitude * std::sin(this->RF_rad_frequency * time);
    }
    if (world.boundary_conditions[global_inputs::number_cells] == 2) {
        this->phi[global_inputs::number_cells] = -this->rho[global_inputs::number_cells]/Constants::epsilon_0;
    } else if (world.boundary_conditions[global_inputs::number_cells] == 4) {
        this->phi[global_inputs::number_cells] = this->RF_half_amplitude * std::sin(this->RF_rad_frequency * time);
    }
    for (int i = 1;i<global_inputs::number_cells;i++){
        this->phi[i] = -this->rho[i]/Constants::epsilon_0;
    }
    this->work_space[0] = this->upper_tri[0]/this->center_tri[0];
    this->phi[0] = this->phi[0]/this->center_tri[0];
    for (int i=1;i<global_inputs::number_cells;i++){
        temp = this->center_tri[i] - this->work_space[i-1] * this->lower_tri[i-1];
        this->work_space[i] = this->upper_tri[i]/temp;
        this->phi[i] = (this->phi[i] - this->phi[i-1]*this->lower_tri[i-1])/temp;
    }
    temp = this->center_tri[global_inputs::number_cells] - this->work_space[global_inputs::number_cells-1] * this->lower_tri[global_inputs::number_cells-1];
    this->phi[global_inputs::number_cells] = (this->phi[global_inputs::number_cells] - this->phi[global_inputs::number_cells-1]*this->lower_tri[global_inputs::number_cells-1])/temp;
    for (int i=global_inputs::number_cells-1;i>=0;i--){
        this->phi[i] -= this->work_space[i] * this->phi[i+1];
    }

}