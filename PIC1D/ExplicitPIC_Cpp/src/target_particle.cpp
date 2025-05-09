
#include "target_particle.h"


Target_Particle::Target_Particle(double mass_in, double temp_in, double density_in, std::string name_in)
    {
    this->name = name_in;
    this->mass = mass_in;
    this->density = density_in;
    this->temperature = temp_in;
    this->v_therm = std::sqrt(this->temperature * Constants::k_boltz / this->mass);
    
 
}






std::vector<Target_Particle> read_target_particle_inputs(const std::string& filename){
    if (Constants::mpi_rank==0) {
        std::cout << " " << std::endl;
        std::cout << "Reading Target Particles " << std::endl;
        std::cout << " -------------------------------- " << std::endl;
    }
    std::vector<Target_Particle> target_particle_list;
    int count_number_particles = 0;

    for (int rank_num = 0;rank_num < Constants::mpi_size; rank_num++){
        if (rank_num == Constants::mpi_rank) {
            std::ifstream file(filename);
            if (!file) {
                std::cerr << "Error: Unable to open file " << filename << std::endl;
            }
            
            std::string line;
            // First check which particles exist and how many there are
            while (std::getline(file, line)) {
                if (line.find("NEUTRALS") != std::string::npos){
                    std::getline(file, line);
                    std::getline(file, line);
                    std::getline(file, line);
                    while (line.find("-------") == std::string::npos) {
                        std::getline(file, line);
                        count_number_particles++;
                    }
                }
            }
            file.clear();
            file.seekg(0, std::ios::beg);
            global_inputs::number_target_particles = count_number_particles;
            target_particle_list.reserve(global_inputs::number_target_particles);
            while (std::getline(file, line)) {
                if (line.find("NEUTRALS") != std::string::npos){
                    std::getline(file, line);
                    std::getline(file, line);
                    std::getline(file, line);
                    while (line.find("-------") == std::string::npos) {
                        std::istringstream iss(line);
                        std::string name;
                        double mass_in, temp_in, density_in;
                        iss >> name >> mass_in >> temp_in >> density_in;
                        mass_in = mass_in * Constants::mass_amu;
                        target_particle_list.emplace_back(mass_in, temp_in, density_in, name);
                        std::getline(file, line);
                    }
                }
            }
            
            file.close();
        }
    }

    int i;
    if (Constants::mpi_rank==0) {
        std::cout << "TNumber of target particles: " << global_inputs::number_target_particles << std::endl;
        for (i = 0; i < global_inputs::number_target_particles; i++){
            std::cout << "----------------- " << std::endl;
            std::cout << "Target particle #: " << i << std::endl;
            std::cout << "Name: " << target_particle_list[i].name << std::endl;
            std::cout << "Mass: " << target_particle_list[i].mass << std::endl;
            std::cout << "Temperature: " << target_particle_list[i].temperature << std::endl;
            std::cout << "Density: " << target_particle_list[i].density << std::endl;
            std::cout << "----------------- " << std::endl;
        }
    }


    return target_particle_list;

}


