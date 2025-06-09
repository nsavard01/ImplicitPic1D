#include "non_linear_solvers/non_linear_solver.hpp"
#include "globals/mpi_vars.hpp"
#include <iostream>
#include <fstream>
#include <sstream> 


void read_non_linear_solver_inputs(const std::string& filename, std::vector<int>& int_params, std::vector<double>& double_params){
    if (mpi_vars::mpi_rank == 0) {
        std::cout << " "  << std::endl;
        std::cout << "Reading non-linear inputs: "  << std::endl;
        std::cout << "-------------------------- "  << std::endl;
        std::string line;
        std::ifstream file(filename);
        if (!file) {
            std::cerr << "Error: Unable to open file " << filename << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        
        
        std::getline(file, line);
        std::istringstream iss(line);
        iss >> int_params[0];
        iss.clear();
        std::getline(file, line);
        iss.str(line);
        iss >> int_params[1];
        iss.clear();
        std::getline(file, line);
        iss.str(line);
        iss >> int_params[2];
        iss.clear();
        std::getline(file, line);
        iss.str(line);
        iss >> double_params[0];
        iss.clear();
        std::getline(file, line);
        iss.str(line);
        iss >> double_params[1];
        iss.clear();
        std::getline(file, line);
        iss.str(line);
        iss >> double_params[2];
        iss.clear();
        file.close();
        std::cout << " "  << std::endl;
    }
}


