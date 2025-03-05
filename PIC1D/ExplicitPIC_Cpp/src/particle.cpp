
#include <vector>
#include "particle.h"
#include "global_inputs.h"
#include "Constants.h"
#include "basic_tools.h"
#include <omp.h>
#include <fstream>
#include <sstream>
#include <math.h>
#include <mpi.h>
#include <numeric>
#include <algorithm>
#include <iomanip>


Particle::Particle(double mass_in, double charge_in, size_t number_in, size_t final_in, std::string name_in)
    {
    int number_threads = omp_get_max_threads();
    int i;
    this->name = name_in;
    this->mass = mass_in;
    this->charge = charge_in;
    this->weight = 0.0;
    this->q_over_m = charge_in/mass_in;
    this->q_times_wp = 0.0;
    this->final_idx = final_in;
    this->accum_wall_loss[0] = this->accum_wall_loss[1] = 0;
    this->accum_energy_loss[0] = this->accum_energy_loss[0] = 0.0;
    this->accum_momentum_loss[0] = this->accum_momentum_loss[0] = 0.0;
    this->number_particles.resize(number_threads);
    this->number_collidable_particles.resize(number_threads);
    size_t expanded_size = static_cast<size_t>(4) * final_idx;
    this->phase_space.resize(number_threads);
    this->work_space.resize(number_threads);
    this->momentum_loss.resize(number_threads);
    this->energy_loss.resize(number_threads);
    this->wall_loss.resize(number_threads);
    this->density.resize(global_inputs::number_nodes, 0.0);
    for (i = 0; i < number_threads; i++){
        this->number_particles[i] = number_in;
        this->number_collidable_particles[i] = number_in;
        this->phase_space[i].resize(expanded_size, 0.0);
        this->work_space[i].resize(global_inputs::number_nodes, 0.0);
        this->momentum_loss[i].resize(2, 0);
        this->energy_loss[i].resize(2, 0);
        this->wall_loss[i].resize(2, 0);
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
    MPI_Allreduce(&total_part, &this->total_number_particles, 1, Constants::mpi_size_t_type, MPI_SUM, MPI_COMM_WORLD);
    this->weight = n_ave * L_domain / static_cast<double>(this->total_number_particles);
    this->q_times_wp = this->charge * this->weight;
}

void Particle::initialize_rand_uniform(double T_ave, const Domain& world) {
    double v_therm = std::sqrt(T_ave * Constants::elementary_charge / this->mass);
    double sum_v_sq = 0.0, sum_v = 0;
    #pragma omp parallel reduction(+:sum_v_sq, sum_v)
    {   
        int thread_id = omp_get_thread_num();
        size_t part_num;
        double x_pos, R_1, R_2, R_3, R_4, v_x, v_y, v_z;
        for (part_num = 0; part_num < this->number_particles[thread_id]; part_num++){
            x_pos = pcg32_random_r() * world.L_domain;
            this->phase_space[thread_id][part_num * 4] = world.get_xi_from_x(x_pos);
            R_1 = pcg32_random_r();
            R_2 = pcg32_random_r();
            R_3 = pcg32_random_r();
            R_4 = pcg32_random_r();
            v_x = v_therm * std::sqrt(-2.0 * std::log(R_1)) * std::cos(2.0 * M_PI * R_2);
            v_y = v_therm * std::sqrt(-2.0 * std::log(R_1)) * std::sin(2.0 * M_PI * R_2);
            v_z = v_therm * std::sqrt(-2.0 * std::log(R_3)) * std::cos(2.0 * M_PI * R_4);
            this->phase_space[thread_id][part_num*4+1] = v_x;
            this->phase_space[thread_id][part_num*4+2] = v_y;
            this->phase_space[thread_id][part_num*4+3] = v_z;
            sum_v += v_x;
            sum_v_sq += v_x * v_x + v_y*v_y + v_z * v_z;

        }
        
    }
    MPI_Allreduce(&sum_v, &this->total_sum_v_x, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&sum_v_sq, &this->total_sum_v_square, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

void Particle::interpolate_particles(){
    
    int thread_id = omp_get_thread_num();
    size_t last_idx = this->number_particles[thread_id]*4;
    double d, xi;
    size_t xi_left, xi_right;
    std::vector<double>& work_space = this->work_space[thread_id];
    std::vector<double>& phase_space = this->phase_space[thread_id];
    std::fill(work_space.begin(), work_space.end(), 0.0);
    for (size_t part_indx = 0; part_indx < last_idx; part_indx += 4){
        xi = phase_space[part_indx];
        xi_left = int(xi);
        xi_right = xi_left+1;
        d = xi - xi_left;
        work_space[xi_left] += (1.0 - d);
        work_space[xi_right] += d;
    } 
}

void Particle::gather_mpi(){
    double sum_v_sq = 0.0, sum_v = 0;
    #pragma omp parallel reduction(+:sum_v_sq, sum_v)
    {
        int thread_id = omp_get_thread_num();
        double v_x, v_y, v_z;
        size_t final_part = this->number_particles[thread_id]*4;
        size_t part_num;
        for (part_num = 0; part_num < final_part; part_num += 4){
            v_x = this->phase_space[thread_id][part_num + 1];
            v_y = this->phase_space[thread_id][part_num + 2];
            v_z = this->phase_space[thread_id][part_num + 3];
            sum_v += v_x;
            sum_v_sq += v_x * v_x + v_y*v_y + v_z * v_z;
        }
    }
    size_t num = std::accumulate(this->number_particles.begin(),this->number_particles.end(), 0);
    MPI_Allreduce(&sum_v, &this->total_sum_v_x, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&sum_v_sq, &this->total_sum_v_square, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&num, &this->total_number_particles, 1, Constants::mpi_size_t_type, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, this->accum_energy_loss, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, this->accum_momentum_loss, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, this->accum_wall_loss, 2, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
}

void Particle::load_density(bool reset_bool){
    if (reset_bool) {
        #pragma omp parallel for
        for (size_t i = 0; i < global_inputs::number_nodes;i++) {
            this->density[i] = 0.0;
        }
    }
    #pragma omp parallel
    {
        this->interpolate_particles();
        #pragma omp barrier
        int thread_id = omp_get_thread_num();
        #pragma omp critical
        {   
            for (int i = 0; i<global_inputs::number_nodes; i++) {
                this->density[i] += this->work_space[thread_id][i];
            }
        }
    }

}

void Particle::write_phase_space(const std::string& dir_name, int diag_num) const{
    for (int rank_num = 0; rank_num < Constants::mpi_size; rank_num++){
        if (Constants::mpi_rank == rank_num) {
            for (int i_thread = 0; i_thread<global_inputs::number_omp_threads;i_thread++){
                size_t size_vec = this->number_particles[i_thread]*4;
                std::vector<double> temp_phase_space(this->phase_space[i_thread].begin(), this->phase_space[i_thread].begin() + size_vec);
                bool append = (Constants::mpi_rank != 0);
                
                std::string filename = dir_name + "/PhaseSpace/phaseSpace_" + this->name + "_thread" + std::to_string(i_thread+1) + ".dat";
                write_vector_to_binary_file(temp_phase_space, filename, 0, 0x00, append);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

void Particle::write_density(const std::string& dir_name, const Domain& world, size_t current_diag, bool average_bool){
    if (world.boundary_conditions[0] == 3) {
        this->density[0] = this->density[0] + this->density[global_inputs::number_cells];
        this->density[global_inputs::number_cells] = this->density[0];
    } else {
        this->density[0] = 2.0*this->density[0];
        this->density[global_inputs::number_cells] = 2.0 * this->density[global_inputs::number_cells];
    }
    for (int i = 0;i<global_inputs::number_nodes; i++) {
        this->density[i] = this->density[i] * this->weight/ world.del_X;
    }
    MPI_Allreduce(MPI_IN_PLACE, this->density.data(), global_inputs::number_nodes, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    std::string filename;
    if (average_bool) {
        filename = dir_name + "/Density/density_" + this->name + "_Average.dat";
        for (int i = 0; i<global_inputs::number_nodes;i++){
            this->density[i] /= current_diag; //current diag is averaging number if averaging
        }
    } else {
        filename = dir_name + "/Density/density_" + this->name + "_" + std::to_string(current_diag) + ".dat";
    }
    if (Constants::mpi_rank == 0) {write_vector_to_binary_file(this->density, filename, 4);}
}


double Particle::get_KE_ave() const{
    // in eV
    return this->total_sum_v_square * this->mass * 0.5 / this->total_number_particles / Constants::elementary_charge;
}

double Particle::get_KE_total() const{
    return this->total_sum_v_square * this->mass * 0.5 * this->weight;
}

double Particle::get_momentum_total() const{
    return this->total_sum_v_x * this->mass;
}

void Particle::write_cell_temperature(const std::string& dir_name, int diag_num) const {
    std::vector<double> temp(global_inputs::number_cells, 0.0);
    std::vector<size_t> counter(global_inputs::number_cells, 0);
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        std::vector<double> temp_thread(global_inputs::number_cells, 0.0);
        std::vector<size_t> counter_thread(global_inputs::number_cells, 0);
        size_t last_idx = this->number_particles[thread_id]*4;
        double v_x, v_y, v_z;
        size_t xi_cell;
        for (size_t part_idx = 0; part_idx < last_idx; part_idx += 4){
            xi_cell = int(this->phase_space[thread_id][part_idx]);
            v_x = this->phase_space[thread_id][part_idx + 1];
            v_y = this->phase_space[thread_id][part_idx + 2];
            v_z = this->phase_space[thread_id][part_idx + 3];
            temp_thread[xi_cell] += v_x*v_x + v_y*v_y + v_z*v_z;
            counter_thread[xi_cell]++;
        }
        #pragma omp critical
        {
            for (int i = 0; i<global_inputs::number_cells;i++){
                temp[i] += temp_thread[i];
                counter[i] += counter_thread[i];
            }
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, temp.data(), global_inputs::number_cells, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, counter.data(), global_inputs::number_cells, Constants::mpi_size_t_type, MPI_SUM, MPI_COMM_WORLD);
    for (int i = 0; i<global_inputs::number_cells;i++){
        if (counter[i] > 0) {
            temp[i] = temp[i] * this->mass / counter[i] / 3.0/Constants::elementary_charge;
        } else {
            temp[i] = 0.0;
        }
    }
    std::string filename = dir_name + "/Temperature/Temp_" + this->name + "_" + std::to_string(diag_num) + ".dat";
    if (Constants::mpi_rank == 0) {write_vector_to_binary_file<double>(temp, filename, 4);}

}

void Particle::initialize_diagnostic_file(const std::string& dir_name) const {
    std::ofstream file(dir_name + "/ParticleDiagnostic_" + this->name + ".dat");
    if (!file) {
        std::cerr << "Error opening file for DateTime \n";
        return;
    }

    file << "Time (s), leftCurrentLoss(A/m^2), rightCurrentLoss(A/m^2), leftPowerLoss(W/m^2), rightPowerLoss(W/m^2), N_p, Temp \n";

    file.close();

}

void Particle::diag_write(const std::string& dir_name, const double& time_diff, const double& current_time, bool average_bool) const {
    double KE_ave = this->get_KE_ave() * 2.0/3.0;
    if (Constants::mpi_rank == 0) {
        if (!average_bool) {
            std::ofstream file(dir_name + "/ParticleDiagnostic_" + this->name + ".dat", std::ios::app);
            if (!file) {
                std::cerr << "Error opening file\n";
                return;
            }
            
            file << std::scientific << std::setprecision(8);
            file << current_time << "\t"
                << this->accum_wall_loss[0] * this->q_times_wp/time_diff << "\t"
                << this->accum_wall_loss[1] * this->q_times_wp/time_diff << "\t"
                << this->accum_energy_loss[0] * this->mass * this->weight * 0.5 / time_diff << "\t"
                << this->accum_energy_loss[1] * this->mass * this->weight * 0.5 / time_diff << "\t"
                << this->total_number_particles << "\t"
                << KE_ave
                <<"\n";

            file.close();
            std::cout << "Number of " << this->name << " is " << total_number_particles << std::endl;
            std::cout << "  " << std::endl;
        } else {
            std::ofstream file(dir_name + "/ParticleAveDiagnostic_" + this->name + ".dat");
            if (!file) {
                std::cerr << "Error opening file\n";
                return;
            }
            file << "Left curr (A/m^2), right curr (A/m^2), left power (W/m^2), right power (W/m^2) \n";
            file << std::scientific << std::setprecision(8);
            file << this->accum_wall_loss[0] * this->q_times_wp/time_diff << "\t"
                << this->accum_wall_loss[1] * this->q_times_wp/time_diff << "\t"
                << this->accum_energy_loss[0] * this->mass * this->weight * 0.5 / time_diff << "\t"
                << this->accum_energy_loss[1] * this->mass * this->weight * 0.5 / time_diff 
                <<"\n";

            file.close();
        }
    }

}



std::vector<Particle> read_particle_inputs(const std::string& filename, const Domain& world){
    if (Constants::mpi_rank == 0) {
        std::cout << " "<< std::endl;
        std::cout << "Reading charged particle inputs "<< std::endl;
        std::cout << "---------------------------------------- "<< std::endl;
    }
    std::vector<Particle> particle_list;
    int count_number_particles = 0;
    std::string line;
    for (int rank_num = 0;rank_num < Constants::mpi_size;rank_num++) {
        if (rank_num == Constants::mpi_rank) {
            std::ifstream file(filename);
            if (!file) {
                std::cerr << "Error: Unable to open file " << filename << std::endl;
            }
            
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
                        std::getline(file, line);
                        count_number_particles++;
                    }
                }
            }
            
            global_inputs::number_charged_particles = count_number_particles;
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
                    size_t num_part_thread;
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
                        size_t num_part_thread;
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
        }
    }

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
        if (Constants::mpi_rank == 0) {
            std::cout << "Particle #: " << i << std::endl;
            std::cout << "Name: " << particle_list[i].name << std::endl;
            std::cout << "Mass: " << particle_list[i].mass << std::endl;
            std::cout << "Charge: " << particle_list[i].charge << std::endl;
            std::cout << "Weight: " << particle_list[i].weight << std::endl;
        }
        double ave_KE = particle_list[i].get_KE_ave();
        if (Constants::mpi_rank == 0) {
            std::cout << "total number of particles: " << particle_list[i].total_number_particles << std::endl;
            std::cout << "Allocated total amount of particles per thread " << particle_list[i].final_idx << std::endl;
            std::cout << "Averge KE: " << ave_KE << std::endl;
            std::cout << "----------------- " << std::endl;
        }
    }

    return particle_list;
    

}


