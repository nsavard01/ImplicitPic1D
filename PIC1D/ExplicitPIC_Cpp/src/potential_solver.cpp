#include "potential_solver.h"
#include "global_inputs.h"
#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <fstream>
#include <sstream>
#include <omp.h>
#include <mpi.h>
#include "Constants.h"

Potential_Solver::Potential_Solver(){
    
}

void Potential_Solver::read_from_file(const std::string& filename, const Domain& world){
    this->lower_tri.resize(global_inputs::number_nodes-1, 0.0);
    this->upper_tri.resize(global_inputs::number_nodes-1, 0.0);
    this->center_tri.resize(global_inputs::number_nodes, 0.0);
    this->phi.resize(global_inputs::number_nodes, 0.0);
    this->EField.resize(global_inputs::number_nodes, 0.0);
    this->rho.resize(global_inputs::number_nodes, 0.0);
    this->work_space.resize(global_inputs::number_nodes, 0.0);
    this->RF_half_amplitude = 0.0;
    this->rho_const = 0.0;
    double left_voltage, right_voltage, RF_frequency;
    std::string line;
    if (Constants::mpi_rank==0){
        std::cout << " "  << std::endl;
        std::cout << "Reading potential inputs: "  << std::endl;
        std::cout << "-------------------------- "  << std::endl;
    }

    for (int rank_num=0;rank_num < Constants::mpi_size; rank_num++) {
        if (rank_num == Constants::mpi_rank) {
            std::ifstream file(filename);
            if (!file) {
                std::cerr << "Error: Unable to open file " << filename << std::endl;
                return;
            }
            
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
        }
    }

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
    if (Constants::mpi_rank==0) { 
        std::cout << "Left voltage " << this->phi[0] << std::endl;
        std::cout << "Right voltage " << this->phi[global_inputs::number_cells] << std::endl;
        std::cout << "Rf frequency " << this->RF_rad_frequency/2.0/M_PI << std::endl;
        std::cout << "RF half amplitude " << this->RF_half_amplitude << std::endl;
        std::cout << "-------------------------- "  << std::endl;
        std::cout << " "  << std::endl;
    }

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
        
        for (int i = left_thread_indx; i <= right_thread_indx; i++) {
            this->rho[i] = 0.0;
        }
        for (int j = 0; j<global_inputs::number_charged_particles;j++){
            for (int i_thread = 0; i_thread < global_inputs::number_omp_threads; i_thread++) {
                for (int i = left_thread_indx; i <= right_thread_indx; i++) {
                    this->rho[i] += particle_list[j].work_space[i_thread][i] * particle_list[j].q_times_wp;
                }
            }
        }
    }
    // Now reduce the rho vector  across all ranks
    MPI_Allreduce(MPI_IN_PLACE, this->rho.data(), global_inputs::number_nodes, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

void Potential_Solver::solve_potential_tridiag(const Domain& world, const double& time) {
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

void Potential_Solver::make_EField(const Domain& world) {
    for (int i = 1; i<global_inputs::number_cells;i++){
        this->EField[i] = 0.5 * (this->phi[i-1] - this->phi[i+1])/world.del_X;
    }
    if (world.boundary_conditions[0]==1 || world.boundary_conditions[0]==4){
        this->EField[0] = (this->phi[0] - this->phi[1])/world.del_X;
    } else if (world.boundary_conditions[0]==3) {
        this->EField[0] = 0.5 * (this->phi[global_inputs::number_cells-1] - this->phi[1])/world.del_X;
        this->EField[global_inputs::number_cells] = this->EField[0];
    }

    if (world.boundary_conditions[global_inputs::number_cells]==1 || world.boundary_conditions[global_inputs::number_cells]==4){
        this->EField[global_inputs::number_cells] = (this->phi[global_inputs::number_cells-1] - this->phi[global_inputs::number_cells])/world.del_X;
    } 

}

inline double Potential_Solver::get_EField_at_loc(const double& xi) const{
    int l_left = int(xi);
    double d = xi - l_left;
    return this->EField[l_left] * (1.0 - d) + this->EField[l_left+1] * d;
}

void Potential_Solver::initial_v_rewind(std::vector<Particle> &particle_list, const double& time_step) {
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        for (int part_type = 0;part_type < global_inputs::number_charged_particles;part_type++){
            Particle& part = particle_list[part_type];
            size_t number_particles = part.number_particles[thread_id];
            std::vector<double>& phase_space = part.phase_space[thread_id];
            size_t part_indx;
            for (size_t part_num = 0; part_num < number_particles; part_num++){
                part_indx = part_num*4;
                phase_space[part_indx+1] -= 0.5 * part.q_over_m * this->get_EField_at_loc(phase_space[part_indx]) * time_step;
            }
        }
    }

}

void Potential_Solver::move_particles(std::vector<Particle> &particle_list, const Domain& world, const double& time_step) {
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        for (int part_type = 0;part_type < global_inputs::number_charged_particles;part_type++){
            Particle& part = particle_list[part_type];
            size_t last_idx = part.number_particles[thread_id]*4;
            part.wall_loss[thread_id][0] = 0; part.wall_loss[thread_id][1] = 0;
            part.energy_loss[thread_id][0] = 0.0; part.energy_loss[thread_id][1] = 0.0;
            std::vector<double>& phase_space = part.phase_space[thread_id];
            size_t space_delete = 0, new_idx;
            double v_x, xi, v_y, v_z;
            for (size_t part_indx= 0; part_indx < last_idx; part_indx += 4){
                xi = phase_space[part_indx];
                v_x = phase_space[part_indx+1];
                v_x += part.q_over_m * this->get_EField_at_loc(xi) * time_step;
                xi += v_x * time_step / world.del_X;

                if (xi <= 0) {
                    switch (world.boundary_conditions[0]){
                        case 1:
                        case 4:
                            v_y = phase_space[part_indx+2];
                            v_z = phase_space[part_indx+3];
                            part.energy_loss[thread_id][0] += v_x*v_x + v_y*v_y + v_z*v_z;
                            part.wall_loss[thread_id][0]++;
                            part.momentum_loss[thread_id][0] += v_x;
                            space_delete = space_delete + 4;
                            break;
                        case 2:
                            xi = - xi;
                            v_x = - v_x;
                            break;
                        case 3:
                            xi = global_inputs::number_cells + xi;
                            break;
                    }
                } else if (xi >= global_inputs::number_cells) {
                    switch (world.boundary_conditions[global_inputs::number_cells]){
                        case 1:
                        case 4:
                            v_y = phase_space[part_indx+2];
                            v_z = phase_space[part_indx+3];
                            part.energy_loss[thread_id][1] += v_x*v_x + v_y*v_y + v_z*v_z;
                            part.wall_loss[thread_id][1]++;
                            part.momentum_loss[thread_id][1] += v_x;
                            space_delete = space_delete + 4;
                            break;
                        case 2:
                            xi = 2.0 * global_inputs::number_cells - xi;
                            v_x = - v_x;
                            break;
                        case 3:
                            xi = xi - global_inputs::number_cells;
                            break;
                    }

                }
                if (xi > 0 && xi < global_inputs::number_cells) {
                    new_idx = part_indx-space_delete;
                    phase_space[new_idx] = xi;
                    phase_space[new_idx+1] = v_x;
                    phase_space[new_idx+2] = phase_space[part_indx+2];
                    phase_space[new_idx+3] = phase_space[part_indx+3];
                }
            }
            part.number_particles[thread_id] = (last_idx - space_delete)/4;
        }
    }
    for (int part_num = 0; part_num < global_inputs::number_charged_particles;part_num++){
        for (int i_thread = 0; i_thread < global_inputs::number_omp_threads; i_thread++){
            particle_list[part_num].number_collidable_particles[i_thread] = particle_list[part_num].number_particles[i_thread];
            particle_list[part_num].accum_wall_loss[0] += particle_list[part_num].wall_loss[i_thread][0];
            particle_list[part_num].accum_wall_loss[1] += particle_list[part_num].wall_loss[i_thread][1];
            particle_list[part_num].accum_energy_loss[0] += particle_list[part_num].energy_loss[i_thread][0];
            particle_list[part_num].accum_energy_loss[1] += particle_list[part_num].energy_loss[i_thread][1];
        }
    }

}