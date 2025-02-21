
#include "target_particle.h"


Target_Particle::Target_Particle(double mass_in, double temp_in, double density_in, std::string name_in)
    {
    this->name = name_in;
    this->mass = mass_in;
    this->density = density_in;
    this->temperature = temp_in;
    this->v_therm = std::sqrt(this->temperature * Constants::k_boltz / this->mass);
    
 
}

void Target_Particle::generate_maxwellian_velocity(double (&v_array)[3]) const {
    // Pass a reference to array and change in place
    double R_1, R_2, R_3, R_4;
    R_1 = pcg32_random_r();
    R_2 = pcg32_random_r();
    R_3 = pcg32_random_r();
    R_4 = pcg32_random_r();
    v_array[0] = this->v_therm * std::sqrt(-2.0 * std::log(R_1)) * std::cos(2.0 * M_PI * R_2);
    v_array[1] = this->v_therm * std::sqrt(-2.0 * std::log(R_1)) * std::sin(2.0 * M_PI * R_2);
    v_array[2] = this->v_therm * std::sqrt(-2.0 * std::log(R_3)) * std::cos(2.0 * M_PI * R_4);

}




std::vector<Target_Particle> read_target_particle_inputs(const std::string& filename){
    std::cout << " " << std::endl;
    std::cout << "Reading Target Particles " << std::endl;
    std::cout << " -------------------------------- " << std::endl;
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
    }
    std::vector<Target_Particle> target_particle_list;
    int count_number_particles = 0;
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

    int i;
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


    return target_particle_list;

}


