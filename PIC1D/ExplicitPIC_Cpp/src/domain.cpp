#include "domain.h"
#include "Constants.h"
#include <cmath>
#include "global_inputs.h"
#include <fstream>
#include <sstream>
#include <mpi.h>

Domain::Domain(){

}

void Domain::read_from_file(const std::string& filename)
    {
    std::string line;
    int left, right;
    if (Constants::mpi_rank == 0) {
        std::cout << " "  << std::endl;
        std::cout << "Reading domain inputs: "  << std::endl;
        std::cout << "-------------------------- "  << std::endl;
    }

    for (int rank_num=0;rank_num<Constants::mpi_size; rank_num++) {
        if (rank_num == Constants::mpi_rank) {
            std::ifstream file(filename);
            if (!file) {
                std::cerr << "Error: Unable to open file " << filename << std::endl;
                return;
            }
            
        
            std::getline(file, line);
            std::istringstream iss(line);
            iss >> this->L_domain;
            iss.clear();
            std::getline(file, line);
            iss.str(line);
            iss >> left >> right;
            file.close();
        }
    }

    int i;
    grid.resize(global_inputs::number_nodes);
    boundary_conditions.resize(global_inputs::number_nodes);
    grid[0] = 0.0;
    grid[global_inputs::number_nodes-1] = this->L_domain;
    boundary_conditions[0] = left;
    boundary_conditions[global_inputs::number_nodes-1] = right;

    
    // Initialize grid spacing
    this->del_X = this->L_domain / static_cast<double>(global_inputs::number_nodes-1);
    for (i = 1; i < global_inputs::number_nodes-1; i++){
        boundary_conditions[i] = 0;
        grid[i] = grid[i-1] + del_X;
    }

    if (global_inputs::number_omp_threads < global_inputs::number_nodes){
        this->num_thread_node_indx = global_inputs::number_omp_threads;
    } else {
        this->num_thread_node_indx = global_inputs::number_nodes;
    }
    this->thread_node_indx.resize(this->num_thread_node_indx);
    int spacing_thread = global_inputs::number_nodes/this->num_thread_node_indx - 1;
    int mod_thread = global_inputs::number_nodes % this->num_thread_node_indx;
    int k = 0;
    for (int i =0; i < this->num_thread_node_indx;i++){
        this->thread_node_indx[i].resize(2);
        this->thread_node_indx[i][0] = k;
        if (i <= mod_thread-1) {
            k = k + spacing_thread + 1;
        } else {
            k = k + spacing_thread;
        }
        this->thread_node_indx[i][1] = k;
        k++;
    }

    if (Constants::mpi_rank == 0) {
        std::cout << "Left boundary type: " << boundary_conditions[0] << std::endl;
        std::cout << "Right boundary type: " << boundary_conditions[global_inputs::number_nodes-1] << std::endl;
        std::cout << "Grid length: " << L_domain << std::endl;
        std::cout << "delta x is: " << del_X << std::endl;
        std::cout << "-----------------------"  << std::endl;
    }
}

double Domain::get_x_from_xi(double xi) const{

    return xi * this->del_X;
}

double Domain::get_xi_from_x(double x) const{

    return x/this->del_X;
}


