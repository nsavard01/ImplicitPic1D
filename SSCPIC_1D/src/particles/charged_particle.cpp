
#include <vector>
#include "particles/charged_particle.hpp"
#include "globals/mpi_vars.hpp"
#include "globals/constants.hpp"
#include "globals/write_functions.hpp"
#include <omp.h>
#include <fstream>
#include <sstream>
#include <math.h>
#include <mpi.h>
#include <numeric>
#include <algorithm>
#include <iomanip>
#include <dirent.h>
#include "rand_gen/maxwell_generator.hpp"


charged_particle::charged_particle(double mass_in, double charge_in, size_t number_in, std::string name_in, int number_nodes){
    this->mass = mass_in;
    this->charge = charge_in;
    this->name = name_in;
    this->q_over_m = charge_in/mass_in;
    this->final_idx.resize(omp_get_max_threads(), number_in);
    this->number_particles.resize(omp_get_max_threads(), 0);
    this->density_grid.resize(omp_get_max_threads());
    this->v_sqr_grid.resize(omp_get_max_threads());
    this->v_drift_grid.resize(omp_get_max_threads());
    this->wall_loss.resize(omp_get_max_threads());
    this->energy_loss.resize(omp_get_max_threads());
    this->wall_freq_rel.resize(omp_get_max_threads());
    this->particle_components.resize(omp_get_max_threads());
    for (int i_thread = 0; i_thread < omp_get_max_threads(); i_thread++){
        this->density_grid[i_thread].resize(number_nodes, 0.0);
        this->v_sqr_grid[i_thread].resize(number_nodes);
        this->v_drift_grid[i_thread].resize(number_nodes);
        for (int node = 0; node < number_nodes; node++) {
            this->v_sqr_grid[i_thread][node].resize(3, 0.0);
            this->v_drift_grid[i_thread][node].resize(3, 0.0);
        }
        this->particle_components[i_thread].resize(number_in);
        this->wall_loss[i_thread].resize(2, 0);
        this->energy_loss[i_thread].resize(2,0.0);
        this->wall_freq_rel[i_thread].resize(2, 0.0);
        for (int part_idx = 0; part_idx < number_in; part_idx++) {
            this->particle_components[i_thread][part_idx].resize(7, 0);
        }
    }
    this->total_v_drift_grid.resize(number_nodes);
    this->total_density_grid.resize(number_nodes, 0.0);
    this->total_v_sqr_grid.resize(number_nodes);
    for (int node = 0; node < number_nodes; node++) {
        this->total_v_sqr_grid[node].resize(3, 0.0);
        this->total_v_drift_grid[node].resize(3, 0.0);
    }
    this->null_collision_bool = false;
}

void charged_particle::print_out() const {
    if (mpi_vars::mpi_rank == 0) {
        std::cout << "Particle name: " << this->name << std::endl;
        std::cout << "Mass (kg): " << this->mass << std::endl;
        std::cout << "Charge (C): " << this->charge << std::endl;
        std::cout << "q/m: " << this->q_over_m << std::endl;
        std::cout << "Maximum number of particles per thread: " << this->final_idx[0] << std::endl;
        std::cout << "number initial particles per thread: " << this->number_particles[0] << std::endl;
        std::cout << " " << std::endl;
    }
}

void charged_particle::read_initial_state(const std::string& dir_name, const domain& world) {
    int total_thread_count = omp_get_max_threads() * mpi_vars::mpi_size;
    for (int rank_num = 0; rank_num < mpi_vars::mpi_size; rank_num++){
        if (mpi_vars::mpi_rank == rank_num) {
            
                std::string filename = dir_name + this->name + ".inp";
                std::string line;
                std::ifstream file(filename);
                if (!file) {
                    std::cerr << "Error: Unable to open file " << filename << std::endl;
                    exit(EXIT_FAILURE);
                }
                while (std::getline(file, line)) {
                        if (line.find("----") != std::string::npos) {
                            std::getline(file, line);
                            std::istringstream iss(line);
                            std::string inj_name;
                            iss >> inj_name;
                            std::getline(file, line);
                            if (inj_name == "Wall" || inj_name == "wall"){
                                // wall injection
                                std::getline(file, line);
                                iss.str(line);
                                int node;
                                iss >> node;
                                if (node != 0 && node != world.number_nodes) {
                                    std::cout << "ERROR: Node for particle injection not on boundary!" << std::endl;
                                    MPI_Abort(MPI_COMM_WORLD, 1);
                                } else if (node == 0) {
                                    if ((world.left_boundary_condition != 1) && (world.left_boundary_condition != 4)){
                                        std::cout << "WARNING: Leftmost node for particle injection not on metallic boundary!" << std::endl;
                                    }
                                } else if (node == world.number_nodes) {
                                    if ((world.right_boundary_condition != 1) && (world.right_boundary_condition != 4)){
                                        std::cout << "WARNING: Rightmost node for particle injection not on metallic boundary!" << std::endl;
                                    }
                                }
                                iss.clear();
                                std::getline(file, line);
                                iss.str(line);
                                double J;
                                iss >> J;
                                J = std::abs(J)/double(total_thread_count);
                                iss.clear();
                                std::getline(file, line);
                                iss.str(line);
                                double v_x;
                                iss >> v_x;
                                if (v_x > 0 && node == world.number_nodes) {
                                    std::cout << "ERROR: V_x > 0 put at rightmost node for particle injection!" << std::endl;
                                    MPI_Abort(MPI_COMM_WORLD, 1);
                                } else if (v_x < 0 && node == 0) {
                                    std::cout << "ERROR: V_x < 0 put at leftmost node for particle injection!" << std::endl;
                                    MPI_Abort(MPI_COMM_WORLD, 1);
                                }
                                // if v_x == 0, then make tiny number so direction is known
                                if (v_x == 0.0 && node == world.number_nodes){
                                    v_x = - constants::machine_eps;
                                } else if (v_x == 0.0 && node == 0){
                                    v_x = constants::machine_eps;
                                }
                                // To avoid issues with other particle operations near boundary, put slightly within domain
                                double new_position;
                                if (v_x > 0) {
                                    new_position = std::nextafter(double(node), node+1);
                                } else {
                                    new_position = std::nextafter(double(node), node-1);
                                }
                                iss.clear();
                                std::getline(file, line);
                                iss.str(line);
                                double v_y;
                                iss >> v_y;
                                iss.clear();
                                std::getline(file, line);
                                iss.str(line);
                                double v_z;
                                iss >> v_z;
                                iss.clear();
                                std::getline(file, line);
                                iss.str(line);
                                double temperature;
                                iss >> temperature;
                                double v_therm_local = std::sqrt(temperature * std::abs(this->q_over_m));
                                iss.clear();
                                std::getline(file, line);
                                iss.str(line);
                                size_t amount_released;
                                iss >> amount_released;
                                iss.clear();
                                std::getline(file, line);


                                // Add to particle list
                                double part_flux = J / std::abs(this->charge) / double(amount_released);
                                for (int i_thread = 0; i_thread < omp_get_max_threads(); i_thread++) {
                                    this->number_particles[i_thread] = amount_released;
                                    double v_x_temp, v_y_temp, v_z_temp;
                                    for (size_t part_idx = 0; part_idx < amount_released; part_idx++) {
                                        std::vector<double>& part_component = this->particle_components[i_thread][part_idx];
                                        part_component[0] = part_flux;
                                        part_component[1] = new_position;
                                        maxwellian_3D_flux(v_x_temp, v_y_temp, v_z_temp, v_therm_local, v_x);
                                        part_component[4] = v_x_temp;
                                        part_component[5] = v_y_temp + v_y;
                                        part_component[6] = v_z_temp + v_z;

                                    }
                                }
                            } else {
                                break;
                            }
                            iss.clear();
                        }  
                }
                file.close();
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}



void charged_particle::write_diagnostics(const std::string& dir_name, int diag_number) const {
    // if (mpi_vars::mpi_rank == 0) {
        
    //     write_vector_to_binary_file(this->temperature, this->temperature.size(), dir_name + "/charged_particles/" + this->name + "/temperature/cell_temp_" + std::to_string(diag_number) + ".dat", 0);

    //     std::ofstream file(dir_name + "/charged_particles/" + this->name + "/momentum_diagnostics.dat", std::ios::app);
    //     if (!file) {
    //         std::cerr << "Error opening file for momentum particle \n";
    //         MPI_Abort(MPI_COMM_WORLD, 1);
    //     }
    //     file << std::scientific << std::setprecision(8);
    //     file << this->total_sum_v[0] << "\t"
    //         << this->total_sum_v[1]  << "\t"
    //         << this->total_sum_v[2] << "\t"
    //         << this->accum_wall_momentum_loss[0][0] << "\t"
    //         << this->accum_wall_momentum_loss[0][1] << "\t"
    //         << this->accum_wall_momentum_loss[0][2] << "\t"
    //         << this->accum_wall_momentum_loss[1][0] << "\t"
    //         << this->accum_wall_momentum_loss[1][1] << "\t"
    //         << this->accum_wall_momentum_loss[1][2]
    //         <<"\n";

    //     file.close();

    //     file.open(dir_name + "/charged_particles/" + this->name + "/energy_diagnostics.dat", std::ios::app);
    //     if (!file) {
    //         std::cerr << "Error opening file for energy particle \n";
    //         MPI_Abort(MPI_COMM_WORLD, 1);
    //     }

    //     file << std::scientific << std::setprecision(8);
    //     double sum_v_sq = this->total_sum_v_square[0] + this->total_sum_v_square[1] + this->total_sum_v_square[2]; 
    //     file << this->total_sum_v_square[0] << "\t"
    //         << this->total_sum_v_square[1] << "\t"
    //         << this->total_sum_v_square[2] << "\t"
    //         << sum_v_sq << "\t"
    //         << this->accum_wall_energy_loss[0] << "\t"
    //         << this->accum_wall_energy_loss[1]
    //         <<"\n";

    //     file.close();

    //     file.open(dir_name + "/charged_particles/" + this->name + "/number_diagnostics.dat", std::ios::app);
    //     if (!file) {
    //         std::cerr << "Error opening file for number particle \n";
    //         MPI_Abort(MPI_COMM_WORLD, 1);
    //     }

    //     file << this->total_number_particles << "\t"
    //         << this->accum_wall_loss[0] << "\t"
    //         << this->accum_wall_loss[1]
    //         <<"\n";

    //     file.close();
    // }
    // this->write_phase_space(dir_name);
}



void charged_particle::reset_diagnostics(int thread_id) {
    this->wall_loss[thread_id][0] = this->wall_loss[thread_id][1] = 0;
    this->energy_loss[thread_id][0] = this->energy_loss[thread_id][1] = 0;
    this->wall_freq_rel[thread_id][0] = this->wall_freq_rel[thread_id][1] = 0;
    const int grid_size = this->density_grid[thread_id].size();
    for (int i = 0; i < grid_size; i++) {
        this->density_grid[thread_id][i] = 0.0;
        for (int j = 0; j < 3; j++) {
            this->v_sqr_grid[thread_id][i][j] = 0.0;
            this->v_drift_grid[thread_id][i][j] = 0.0;
        }
    }
}



void charged_particle::gather_mpi(){
    int size_grid = this->total_density_grid.size();
    this->total_number_particles = 0;
    this->accum_wall_loss[0] = this->accum_wall_loss[1] = 0;
    this->accum_wall_energy_loss[0] = this->accum_wall_energy_loss[1] = 0;
    this->accum_wall_freq_rel[0] = this->accum_wall_freq_rel[1] = 0;
    std::fill(this->total_density_grid.begin(), this->total_density_grid.end(), 0.0);
    for (int node = 0; node < size_grid; node++) {
        std::fill(this->total_v_sqr_grid[node].begin(), this->total_v_sqr_grid[node].end(), 0.0);
        std::fill(this->total_v_drift_grid[node].begin(), this->total_v_drift_grid[node].end(), 0.0);
    }
    for (int i = 0; i < omp_get_max_threads(); i++){
        this->total_number_particles += this->number_particles[i];
        this->accum_wall_loss[0] += this->wall_loss[i][0];
        this->accum_wall_loss[1] += this->wall_loss[i][1];
        this->accum_wall_energy_loss[0] += this->energy_loss[i][0];
        this->accum_wall_energy_loss[1] += this->energy_loss[i][1];
        this->accum_wall_freq_rel[0] += this->wall_freq_rel[i][0];
        this->accum_wall_freq_rel[1] += this->wall_freq_rel[i][1];
        for (int j = 0; j < size_grid; j++) {
            this->total_density_grid[j] += this->density_grid[i][j];
            for (int u = 0; u < 3; u++) {
                this->total_v_sqr_grid[j][u] += this->v_sqr_grid[i][j][u];
                this->total_v_drift_grid[j][u] += this->v_drift_grid[i][j][u];
            }
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, &this->total_number_particles, 1, mpi_vars::mpi_size_t_type, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, this->accum_wall_energy_loss, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, this->accum_wall_loss, 2, mpi_vars::mpi_size_t_type, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, this->accum_wall_freq_rel, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, this->total_density_grid.data(), size_grid, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    for (int node = 0; node < size_grid; node++) {
        MPI_Allreduce(MPI_IN_PLACE, this->total_v_sqr_grid[node].data(), 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, this->total_v_drift_grid[node].data(), 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
    
    for (int j = 0; j < size_grid; j++) {
        for (int u = 0; u < 3; u++) {
            this->total_v_sqr_grid[j][u] /= this->total_density_grid[j];
            this->total_v_drift_grid[j][u] /= this->total_density_grid[j];
        }
    }
    
}

void charged_particle::load_to_target(std::vector<target_particle>& target_particle_list, const domain& world) {
    // Load particle properties to target if exist
    target_particle& local_target = target_particle_list[this->target_particle_idx];
    if (this->name == local_target.name) {
        // replace properties with last iteration from charged particle
        int size_grid = this->total_density_grid.size();
        std::vector<double>& target_density = local_target.density;
        std::vector<std::vector<double>>& target_v_therm_sqr = local_target.v_therm_sqr;
        std::vector<std::vector<double>>& target_v_drift = local_target.v_drift;
        if (world.domain_type == 0) {

        } else {
            const std::vector<double>& dx_dxi = world.dx_dxi;
            double dx;
            dx = 0.5 * dx_dxi[0];
            target_density[0] = this->total_density_grid[0] / dx;
            // if (mpi_vars::mpi_rank == 0) {
            //     std::cout << "i: " << 0 << " J: " << target_density[0] << std::endl;
            // }
            for (int j = 1; j< size_grid-1; j++) {
                dx = 0.5 * (dx_dxi[j-1] + dx_dxi[j]);
                target_density[j] = this->total_density_grid[j] / dx;
                // if (mpi_vars::mpi_rank == 0) {
                //     std::cout << "j: " << j << " J: " << target_density[j] << std::endl;
                // }
            }
            dx = 0.5 * dx_dxi[size_grid-2];
            target_density[size_grid-1] = this->total_density_grid[size_grid-1] / dx;
            // if (mpi_vars::mpi_rank == 0) {
            //     std::cout << "j: " << size_grid-1 << " J: " << target_density[size_grid-1] << std::endl;
            // }
        }
        // v_therm based on v_max^2 =  <v^2> - <v>^2 (maxwellian is variance)
        // then v_therm = sqrt(v_max^2/3) 
        for (int node = 0; node < size_grid; node++) {
            for (int u = 0; u < 3; u++){
                double v_temp = this->total_v_drift_grid[node][u];
                target_v_therm_sqr[node][u] = (this->total_v_sqr_grid[node][u] - v_temp * v_temp) / 3.0;
                target_v_drift[node][u] = v_temp;
            }
        }
    }
}


std::vector<charged_particle> read_charged_particle_inputs(const std::string& directory_path, const domain& world){
    
    std::vector<charged_particle> particle_list;
    std::vector<double> mass_in, charge_in;
    std::vector<size_t> num_part_thread;
    std::vector<std::string> particle_names;
    std::vector<int> index_order, number_space_coordinates, number_velocity_coordinates;
    int count_number_particles = 0;
    if (mpi_vars::mpi_rank == 0) {
        std::cout << " "<< std::endl;
        std::cout << "Reading charged particle inputs "<< std::endl;
        std::cout << "---------------------------------------- "<< std::endl;
    }
    for (int rank_num = 0; rank_num < mpi_vars::mpi_size; rank_num++){
        if (mpi_vars::mpi_rank == rank_num) {
            // Open the directory
            DIR* dir = opendir(directory_path.c_str());
            if (!dir) {
                perror("opendir");
                exit(EXIT_FAILURE);
            }

            struct dirent* entry;
            while ((entry = readdir(dir)) != nullptr) {
                if (entry->d_type == DT_REG) {  // regular file
                    count_number_particles++;
                    std::string filename = directory_path + entry->d_name;
                    std::string line;
                    std::ifstream file(filename);
                    if (!file) {
                        std::cerr << "Error: Unable to open file " << filename << std::endl;
                        exit(EXIT_FAILURE);
                    }
                    std::getline(file, line);
                    std::istringstream iss(line);
                    std::string name;
                    iss >> name;
                    particle_names.push_back(name);
                    iss.clear();
                    std::getline(file, line);
                    iss.str(line);
                    double charge;
                    iss >> charge;
                    charge_in.push_back(charge * constants::elementary_charge);
                    iss.clear();
                    std::getline(file, line);
                    iss.str(line);
                    double mass;
                    iss >> mass;
                    mass = mass * constants::mass_amu;
                    if (std::abs(constants::electron_mass - mass ) /constants::electron_mass < 1e-3 ) {
                        mass = constants::electron_mass;
                    } else {
                        mass = mass - charge * constants::electron_mass; // assume put in neutral mass, so subtract electron mass for momentum/energy conservation in collisions
                    }
                    mass_in.push_back(mass);
                    iss.clear();
                    std::getline(file, line);
                    iss.str(line);
                    size_t number_part;
                    iss >> number_part;
                    num_part_thread.push_back(number_part);
                    file.close();
                }
            }   
        
            closedir(dir);

            // Initialize index_order with indices [0, 1, 2, ..., count_number_particles - 1]
            index_order.resize(count_number_particles);
            std::iota(index_order.begin(), index_order.end(), 0);

            
            // Sort indices based on charge-to-mass ratio (q/m)
            std::sort(index_order.begin(), index_order.end(), [&](int i, int j) {
                return (charge_in[i] / mass_in[i]) < (charge_in[j] / mass_in[j]);
            });
        
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    

    //Create the particles in the sorted order 
    for (int num_part = 0; num_part < count_number_particles; num_part++){
        int i = index_order[num_part];
        std::string name = particle_names[i];
        double mass = mass_in[i];
        double charge = charge_in[i];
        size_t number_in = num_part_thread[i];
        charged_particle temp_particle(mass, charge, number_in, name, world.number_nodes);
        temp_particle.read_initial_state(directory_path, world);
        particle_list.push_back(temp_particle);
    }

    for (int i = 0; i < particle_list.size(); i++) {
        particle_list[i].print_out();
    }

    return particle_list;
    

}


void find_corresponding_targets(std::vector<charged_particle>& charged_particle_list, std::vector<target_particle>& target_particle_list) {
    for (int c_idx = 0; c_idx < charged_particle_list.size(); c_idx++) {
        bool found = false;
        for (int t_idx = 0; t_idx < target_particle_list.size(); t_idx++) {
            if (target_particle_list[t_idx].name == charged_particle_list[c_idx].name) {
                found = true;
                target_particle_list[t_idx].charged_particle_idx = c_idx;
                charged_particle_list[c_idx].target_particle_idx = t_idx;
            }
        }
        if (mpi_vars::mpi_rank == 0) {
            if (!found) {
                std::cout << "Charged particle " << charged_particle_list[c_idx].name << " doesn't have a corresponding target for deposition!" << std::endl;
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }
    }
    
}


