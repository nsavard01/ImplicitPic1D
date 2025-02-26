
#include <vector>
#include "particle.h"
#include "global_inputs.h"
#include "Constants.h"
#include <omp.h>
#include <fstream>
#include <sstream>
#include <math.h>
#include <numeric>
#include <algorithm>


Particle::Particle(double mass_in, double charge_in, uint32_t number_in, size_t final_in, std::string name_in)
    {
    int number_threads = omp_get_max_threads();
    int i;
    name = name_in;
    mass = mass_in;
    charge = charge_in;
    weight = 0.0;
    q_over_m = charge/mass;
    q_times_wp = 0.0;
    final_idx = final_in;
    accum_wall_loss[0] = accum_wall_loss[1] = 0;
    accum_energy_loss[0] = accum_energy_loss[0] = 0.0;
    number_particles.resize(number_threads);
    number_collidable_particles.resize(number_threads);
    size_t expanded_size = static_cast<size_t>(4) * final_idx;
    phase_space.resize(number_threads);
    work_space.resize(number_threads);
    momentum_loss.resize(number_threads);
    energy_loss.resize(number_threads);
    wall_loss.resize(number_threads);
    for (i = 0; i < number_threads; i++){
        number_particles[i] = number_in;
        number_collidable_particles[i] = number_in;
        phase_space[i].resize(expanded_size, 0.0);
        work_space[i].resize(global_inputs::number_nodes, 0.0);
        momentum_loss[i].resize(2, 0);
        energy_loss[i].resize(2, 0);
        wall_loss[i].resize(2, 0);
    }
    
 
}



void Particle::initialize_weight(double n_ave, double L_domain) {
    size_t total_part;
    int i;
    int number_threads = omp_get_max_threads();
    total_part = 0;
    for (i = 0; i < number_threads; i++){
       total_part += this->number_particles[i];
    }
    weight = n_ave * L_domain / static_cast<double>(total_part);
    q_times_wp = charge * weight;
}

void Particle::initialize_rand_uniform(double T_ave, const Domain& world) {
    double v_therm = std::sqrt(T_ave * Constants::elementary_charge / this->mass);
    #pragma omp parallel
    {   
        int thread_id = omp_get_thread_num();
        size_t part_num;
        double x_pos, R_1, R_2, R_3, R_4;
        for (part_num = 0; part_num < this->number_particles[thread_id]; part_num++){
            x_pos = pcg32_random_r() * world.L_domain;
            this->phase_space[thread_id][part_num * 4] = world.get_xi_from_x(x_pos);
            R_1 = pcg32_random_r();
            R_2 = pcg32_random_r();
            R_3 = pcg32_random_r();
            R_4 = pcg32_random_r();
            this->phase_space[thread_id][part_num*4+1] = v_therm * std::sqrt(-2.0 * std::log(R_1)) * std::cos(2.0 * M_PI * R_2);
            this->phase_space[thread_id][part_num*4+2] = v_therm * std::sqrt(-2.0 * std::log(R_1)) * std::sin(2.0 * M_PI * R_2);
            this->phase_space[thread_id][part_num*4+3] = v_therm * std::sqrt(-2.0 * std::log(R_3)) * std::cos(2.0 * M_PI * R_4);

        }
        
    }
}

void Particle::interpolate_particles(){
    
    int thread_id = omp_get_thread_num();
    size_t part_num;
    size_t number_particles = this->number_particles[thread_id];
    double d, xi;
    size_t xi_left, xi_right;
    std::fill(this->work_space[thread_id].begin(), this->work_space[thread_id].end(), 0.0);
    for (part_num = 0; part_num < number_particles; part_num++){
        xi = this->phase_space[thread_id][part_num * 4];
        xi_left = static_cast<size_t>(xi);
        xi_right = xi_left+1;
        d = xi - xi_left;
        this->work_space[thread_id][xi_left] += (1.0 - d);
        this->work_space[thread_id][xi_right] += d;
    }   
}

double Particle::get_KE_ave() const{
    // in eV
    double sum = 0.0;
    double num_part_total = 0.0;
    #pragma omp parallel reduction(+:sum, num_part_total)
    {   
        int thread_id = omp_get_thread_num();
        double v_x, v_y, v_z;
        size_t part_num;
        for (part_num = 0; part_num < this->number_particles[thread_id]; part_num++){
            v_x = this->phase_space[thread_id][part_num * 4 + 1];
            v_y = this->phase_space[thread_id][part_num * 4 + 2];
            v_z = this->phase_space[thread_id][part_num * 4 + 3];
            sum += v_x * v_x + v_y*v_y + v_z * v_z;
        }
        num_part_total += this->number_particles[thread_id];
    }
    
    return sum * this->mass * 0.5 / num_part_total / Constants::elementary_charge;
}

double Particle::get_KE_total() const{
    // in eV
    double sum = 0.0;
    #pragma omp parallel reduction(+:sum)
    {   
        int thread_id = omp_get_thread_num();
        double v_x, v_y, v_z;
        size_t part_num;
        for (part_num = 0; part_num < this->number_particles[thread_id]; part_num++){
            v_x = this->phase_space[thread_id][part_num * 4 + 1];
            v_y = this->phase_space[thread_id][part_num * 4 + 2];
            v_z = this->phase_space[thread_id][part_num * 4 + 3];
            sum += v_x * v_x + v_y*v_y + v_z * v_z;
        }
    }
    
    return sum * this->mass * 0.5 * this->weight;
}



std::vector<Particle> read_particle_inputs(const std::string& filename, const Domain& world){
    std::cout << " "<< std::endl;
    std::cout << "Reading charged particle inputs "<< std::endl;
    std::cout << "---------------------------------------- "<< std::endl;
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
    }
    std::vector<Particle> particle_list;
    int count_number_particles = 0;
    std::string line;
    // First check which particles exist and how many there are
    while (std::getline(file, line)) {
        if (line.find("ELECTRONS") != std::string::npos){
            std::getline(file, line);
            std::getline(file, line);
            std::getline(file, line);
            std::istringstream iss(line);
            std::string name;
            uint32_t num_part_thread;
            size_t factor;
            iss >> name >> num_part_thread >> factor;
            std::cout << name << " " << num_part_thread << " " << factor << std::endl;
            std::getline(file, line);
            count_number_particles++;
        }

        if (line.find("IONS") != std::string::npos){
            std::getline(file, line);
            std::getline(file, line);
            std::getline(file, line);
            while (line.find("-------") == std::string::npos) {
                std::istringstream iss(line);
                std::string name;
                uint32_t num_part_thread;
                double mass_in, charge_in;
                size_t factor;
                iss >> name >> mass_in >> charge_in >> num_part_thread >> factor;
                std::cout << name << " " << " " << mass_in << " " << charge_in << " " << num_part_thread << " " << factor << std::endl;
                std::getline(file, line);
                count_number_particles++;
            }
        }
    }
    
    global_inputs::number_charged_particles = count_number_particles;
    std::cout << global_inputs::number_charged_particles << std::endl;
    particle_list.reserve(global_inputs::number_charged_particles);
    // get electrons
    file.clear();
    file.seekg(0, std::ios::beg);
    while (std::getline(file, line)) {
        if (line.find("ELECTRONS") != std::string::npos){
            std::getline(file, line);
            std::getline(file, line);
            std::getline(file, line);
            std::istringstream iss(line);
            std::string name;
            uint32_t num_part_thread;
            size_t factor;
            iss >> name >> num_part_thread >> factor;
            factor = factor * static_cast<size_t>(num_part_thread);
            particle_list.emplace_back(Constants::electron_mass, -Constants::elementary_charge, num_part_thread, factor, name);
            std::getline(file, line);
        }
    }

    file.clear();
    file.seekg(0, std::ios::beg);
    while (std::getline(file, line)) {
        if (line.find("IONS") != std::string::npos){
            std::getline(file, line);
            std::getline(file, line);
            std::getline(file, line);
            while (line.find("-------") == std::string::npos) {
                std::istringstream iss(line);
                std::string name;
                uint32_t num_part_thread;
                double mass_in, charge_in;
                size_t factor;
                iss >> name >> mass_in >> charge_in >> num_part_thread >> factor;
                factor = factor * static_cast<size_t>(num_part_thread);
                mass_in = mass_in * Constants::mass_amu - charge_in * Constants::electron_mass;
                charge_in = charge_in * Constants::elementary_charge;
                particle_list.emplace_back(mass_in, charge_in, num_part_thread, factor, name);
                std::getline(file, line);
            }
        }
    }
    file.close();

    int i;
    double T_ave;
    for (i = 0; i < global_inputs::number_charged_particles; i++){
        particle_list[i].initialize_weight(global_inputs::initial_density, world.L_domain);
        if (particle_list[i].mass == Constants::electron_mass) {
            T_ave = global_inputs::temp_electrons;
        }
        else {
            T_ave = global_inputs::temp_ions;
        }
        particle_list[i].initialize_rand_uniform(T_ave, world);
        std::cout << "Particle #: " << i << std::endl;
        std::cout << "Name: " << particle_list[i].name << std::endl;
        std::cout << "Mass: " << particle_list[i].mass << std::endl;
        std::cout << "Charge: " << particle_list[i].charge << std::endl;
        std::cout << "Weight: " << particle_list[i].weight << std::endl;
        size_t tot_sum = std::accumulate(particle_list[i].number_particles.begin(),
            particle_list[i].number_particles.end(), 0);
        std::cout << "total number of particles: " << tot_sum << std::endl;
        std::cout << "Allocated total amount of particles per thread " << particle_list[i].final_idx << std::endl;
        std::cout << "Averge KE: " << particle_list[i].get_KE_ave() << std::endl;
        std::cout << "----------------- " << std::endl;
    }

    return particle_list;
    

}


