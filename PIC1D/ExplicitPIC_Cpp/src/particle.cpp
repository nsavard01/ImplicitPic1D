
#include <vector>
#include "particle.h"
#include "global_inputs.h"
#include "domain.h"
#include <fstream>
#include <sstream>


// Accessor for 3D indexing
inline double& Particle::phase_space_at(size_t phase_idx, size_t part_idx, size_t thread_idx) {
    return phase_space[phase_idx + part_idx * 4 + thread_idx * 4 * final_idx];
}

void read_particle_inputs(const std::string& filename){
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        return;
    }

    std::string line;
    int line_count = 0;
    while (std::getline(file, line)) {
        std::cout << line << std::endl;
        // std::istringstream iss(line);
        // if (line.find("END") != std::string::npos) break; // Stop at END marker
        
        // if (line_count == 0) {
        //     iss >> number_nodes;
        // } else if (line_count == 1) {
        //     iss >> L_domain;
        // } else if (line_count == 2) {
            
        //     iss >> left >> right;
        // }
        line_count++;
    }


    

}


