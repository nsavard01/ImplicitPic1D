#include "domain.h"
#include <cmath>
#include "global_inputs.h"
#include <fstream>
#include <sstream>

// Constructor definition
Domain::Domain(const std::string& filename)
    {
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        return;
    }
    
    std::string line;
    int line_count = 0;
    int left, right;
    std::cout << "Reading domain inputs: "  << std::endl;
    std::cout << "-------------------------- "  << std::endl;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        if (line.find("END") != std::string::npos) break; // Stop at END marker
        
        if (line_count == 0) {
            iss >> number_nodes;
        } else if (line_count == 1) {
            iss >> L_domain;
        } else if (line_count == 2) {
            
            iss >> left >> right;
        }
        line_count++;
    }
        
    int i;
    grid.resize(number_nodes);
    boundary_conditions.resize(number_nodes);
    grid[0] = 0.0;
    grid[number_nodes-1] = L_domain;
    boundary_conditions[0] = left;
    boundary_conditions[number_nodes-1] = right;

    
    // Initialize grid spacing
    del_X = L_domain / static_cast<double>(number_nodes-1);
    for (i = 1; i < number_nodes-1; i++){
        boundary_conditions[i] = 0;
        grid[i] = grid[i-1] + del_X;
    }

    std::cout << "Number of nodes: " << number_nodes << std::endl;
    std::cout << "Left boundary type: " << boundary_conditions[0] << std::endl;
    std::cout << "Right boundary type: " << boundary_conditions[number_nodes-1] << std::endl;
    std::cout << "Grid length: " << L_domain << std::endl;
    std::cout << "delta x is: " << del_X << std::endl;
    std::cout << "-----------------------"  << std::endl;
}

double Domain::get_x_from_xi(double xi) {

    return xi * del_X;
}

double Domain::get_xi_from_x(double x) {

    return x/del_X;
}

// Method definition
void Domain::readWorld() {
    std::cout << "Reading world data..." << std::endl;
    // Placeholder for file reading or other initialization
}
