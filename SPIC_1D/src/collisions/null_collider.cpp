#include "collisions/null_collider.hpp"
#include "globals/constants.hpp"
#include "globals/mpi_vars.hpp"
#include <regex>
#include <cstdlib>
#include <omp.h>
#include <algorithm>
#include <numeric>
#include <mpi.h>
#include <iostream>
#include <iomanip>
#include <dirent.h>
#include <fstream>
#include <sstream>




null_collider::null_collider(int number_collisions, int primary_idx, const std::vector<int>& target_idx, 
    const std::vector<std::vector<double>> &sigma_array, 
    const std::vector<std::vector<double>> &energy_array, const std::vector<double> &energy_threshold,
    const std::vector<int> &collision_type, const std::vector<std::vector<int>> &products_indx)
    {
    this->number_collisions = number_collisions;
    this->primary_idx = primary_idx;
    this->target_idx = target_idx;
    this->sigma_array = sigma_array;
    this->energy_array = energy_array;
    this->energy_threshold = energy_threshold;
    this->collision_type = collision_type;
    this->products_indx = products_indx;
    // this->reduced_mass = mass_inputs[0];
    // this->mass_sum = mass_inputs[1];
    // this->reduced_mass_triple = mass_inputs[2];
    if (number_collisions > 0) {
        this->total_incident_energy.resize(this->number_collisions, 0.0);
        this->total_energy_loss.resize(this->number_collisions, 0.0);
        this->total_amount_collisions.resize(this->number_collisions, 0);
    }
    
    // int l, u;
    // this->sigma_v_max = 0.0;
    // double sum_sigma, v_r;
    // for (l=0;l<this->length_arrays;l++){
    //     sum_sigma = 0.0;
    //     v_r = std::sqrt(2.0 * this->energy_array[l] * Constants::elementary_charge / this->reduced_mass);
    //     for (u=0;u<this->number_collisions; u++){
    //         sum_sigma += this->sigma_array[u][l];
    //     }
    //     this->sigma_v_max = std::max(this->sigma_v_max, sum_sigma * v_r);
    // }

    // if (Constants::mpi_rank == 0) {
    //     std::cout << "------------------------------ " << std::endl;
    //     std::cout << "Primary particle is idx is " << this->primary_idx << " " << this->target_idx << std::endl;
    //     std::cout << "Amount of collisions is " << this->number_collisions << std::endl;
    //     std::cout << "Reduced mass: " << this->reduced_mass << std::endl;
    //     std::cout << "Sum mass: " << this->mass_sum << std::endl;
    //     std::cout << "Reduced mass triple product: " << this->reduced_mass_triple << std::endl;
    //     std::cout << "Sigma_v_max: " << this->sigma_v_max << std::endl;
    //     std::cout << "----- " << std::endl;
    //     std::cout << std::endl;
    //     for (l=0; l<this->number_collisions;l++){
    //         std::cout << "Collision # " << l << std::endl;
    //         std::cout << "Collision type integer is " << this->collision_type[l] << std::endl;
    //         std::cout << "Threshold energy (eV) " << this->energy_threshold[l] << std::endl;
    //         std::cout << "Products indices are: ";
    //         for (u=0;u<this->products_indx[l].size();u++){
    //             std::cout << this->products_indx[l][u] << " ";
    //         }
    //         std::cout << std::endl;
    //         // std::cout << std::endl;
    //         // for (u=0; u < this->length_arrays;u++){
    //         //     std::cout << this->energy_array[u] << " " << this->sigma_array[l][u] << std::endl;  
    //         // }
    //         std::cout << "----- " << std::endl;
    //         std::cout << std::endl;
    //     }
    //     std::cout << "------------------------------ " << std::endl;
    // }
    
 
}

void set_null_frequency(const std::vector<charged_particle>& particle_list, const std::vector<target_particle>& target_particle_list) {
    // double primary_mass = particle_list[this->primary_idx].mass;
    // double freq_sum = 0.0;
    // std::vector<double> reduced_mass(this->number_collisions);
    // std::vector<double> energy_current(this->number_collisions, 0.0);
    // std::vector<int> idx_current(this->number_collisions, 0);
    // std::vector<double> target_density(this->number_collisions);
    // for (int coll_idx = 0; coll_idx < this->number_collisions; coll_idx++) {
    //     double target_mass = target_particle_list[this->target_idx[coll_idx]].mass;
    //     target_density[coll_idx] = target_particle_list[this->target_idx[coll_idx]].ave_density;
    //     reduced_mass[coll_idx] = (primary_mass * target_mass)/ (primary_mass + target_mass);
    // }

    // bool not_complete = true;
    // double v_r;
    // while (not_complete) {
    //     freq_sum = 0.0;

    //     for (int coll_idx = 0; coll_idx < this->number_collisions; coll_idx++){
    //         int idx_coll = idx_current[coll_idx];
    //         double E = this->energy_array[idx_coll]
    //         v_r = std::sqrt(2.0 * this->energy_array[l] * Constants::elementary_charge / this->reduced_mass);
    //         sum_sigma += this->sigma_array[u][l];
    //     }
    //     this->sigma_v_max = std::max(this->sigma_v_max, sum_sigma * v_r);
    // }
}

// void Null_Collision::initialize_data_files(const std::string& dir_name, std::vector<Particle>& particle_list, std::vector<Target_Particle>& target_particle_list) const {

//     if (Constants::mpi_rank == 0) {
//         std::string binary_folder = dir_name + "/BinaryCollisions/" + particle_list[this->primary_idx].name + "_on_" + target_particle_list[this->target_idx].name;
//         if (!createDirectory(binary_folder)) {
//             std::cerr << "Save directory not successfully created!" << std::endl;
//             MPI_Abort(MPI_COMM_WORLD, 1);
//         }
//         std::ofstream file(binary_folder + "/CollisionProperties.dat");
//         if (!file) {
//             std::cerr << "Error opening file \n";
//             return;
//         }

//         file << "Coll #, collType, E_thres (eV), maxSigma (m^2), EatMaxSigma (eV) \n";
//         file << std::scientific << std::setprecision(8);
        
//         for (int i=0; i< this->number_collisions;i++){
//             size_t max_indx;
//             double max_val = 0.0;
//             for (int j= 0; j< this->length_arrays;j++){
//                 if (this->sigma_array[i][j] > max_val) {
//                     max_indx = j;
//                     max_val = this->sigma_array[i][j];
//                 }
//             }
//             file << i << "\t"
//                 << this->collision_type[i] << "\t"
//                 << this->energy_threshold[i] << "\t"
//                 << this->sigma_array[i][max_indx] << "\t"
//                 << this->energy_array[max_indx]
//                 <<"\n";
//         }
        
//         file.close();
//         for (int i=0; i< this->number_collisions;i++){
//             file.open(binary_folder + "/CollisionDiag_" + std::to_string(i+1) + ".dat");
//             file << "CollRatio, AveEnergyLoss (eV), AveIncidentEnergy (eV), P_loss(W/m^2), aveCollFreq (Hz/m^2) \n";
//             file.close();
//         }
        
//     }
// }


// void Null_Collision::gather_mpi() {
//     MPI_Allreduce(MPI_IN_PLACE, this->total_incident_energy.data(), this->number_collisions, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//     MPI_Allreduce(MPI_IN_PLACE, this->total_energy_loss.data(), this->number_collisions, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//     MPI_Allreduce(MPI_IN_PLACE, this->total_amount_collisions.data(), this->number_collisions, Constants::mpi_size_t_type, MPI_SUM, MPI_COMM_WORLD);
//     MPI_Allreduce(MPI_IN_PLACE, &this->total_amount_collidable_particles, 1, Constants::mpi_size_t_type, MPI_SUM, MPI_COMM_WORLD);
// }

// void Null_Collision::diag_write(const std::string& dir_name, std::vector<Particle>& particle_list, std::vector<Target_Particle>& target_particle_list, const double& time_diff, bool average_bool) {
//     if (Constants::mpi_rank == 0) {
//         std::string binary_folder = dir_name + "/BinaryCollisions/" + particle_list[this->primary_idx].name + "_on_" + target_particle_list[this->target_idx].name;
//         if (!average_bool) {   
//             for (int i=0; i< this->number_collisions;i++){
//                 std::ofstream file(binary_folder + "/CollisionDiag_" + std::to_string(i+1) + ".dat", std::ios::app);
//                 file << std::scientific << std::setprecision(8);
//                 file << double(this->total_amount_collisions[i])/double(this->total_amount_collidable_particles) << "\t"
//                     << this->total_energy_loss[i] * 0.5 / Constants::elementary_charge/double(this->total_amount_collisions[i]) << "\t"
//                     << this->total_incident_energy[i] * particle_list[this->primary_idx].mass * 0.5 / Constants::elementary_charge/double(this->total_amount_collisions[i]) << "\t"
//                     << this->total_energy_loss[i] * 0.5 * particle_list[this->primary_idx].weight / time_diff << "\t"
//                     << double(this->total_amount_collisions[i]) * particle_list[this->primary_idx].weight / time_diff
//                     << "\n";
//                 file.close();
//             }
//         } else {
//             std::ofstream file(binary_folder + "/AveCollisionDiag.dat");
//             file << "Coll #, CollRatio, AveEnergyLoss (eV), AveIncidentEnergy (eV), P_loss(W/m^2), aveCollFreq (Hz/m^2) \n";
//             file << std::scientific << std::setprecision(8);
//             for (int i=0; i< this->number_collisions;i++){
//                 file << i+1 << "\t"
//                     << double(this->total_amount_collisions[i])/double(this->total_amount_collidable_particles) << "\t"
//                     << this->total_energy_loss[i] * 0.5 / Constants::elementary_charge/double(this->total_amount_collisions[i]) << "\t"
//                     << this->total_incident_energy[i] * particle_list[this->primary_idx].mass * 0.5 / Constants::elementary_charge/double(this->total_amount_collisions[i]) << "\t"
//                     << this->total_energy_loss[i] * 0.5 * particle_list[this->primary_idx].weight / time_diff << "\t"
//                     << double(this->total_amount_collisions[i]) * particle_list[this->primary_idx].weight / time_diff
//                     << "\n";
//             }
//             file.close();
//         }
        
//     }
// }



// void Null_Collision::generate_null_collisions(std::vector<Particle> &particle_list, std::vector<Target_Particle> &target_particle_list, double time_step){
    
    
//     Particle& primary_particle = particle_list[primary_idx];
//     Target_Particle& target_particle = target_particle_list[target_idx];
//     double primary_mass = primary_particle.mass;
//     double target_mass = target_particle.mass;
    
//     this->total_amount_collidable_particles += std::accumulate(primary_particle.number_collidable_particles.begin(),
//         primary_particle.number_collidable_particles.end(), 0);

//     double P_null = 1.0 - std::exp(-this->sigma_v_max * target_particle.density * time_step);
//     if (P_null > 0.05){
//         std::cout << "P_null greater than 5% " << std::endl;
//         std::exit(EXIT_FAILURE);
//     }

//     #pragma omp parallel
//     {
//         int thread_id = omp_get_thread_num();
//         int iter, coll_idx;
//         std::vector<size_t> total_collisions(this->number_collisions,0);
//         std::vector<double> energy_loss(this->number_collisions, 0.0), tot_incident_energy(this->number_collisions, 0.0);
//         size_t number_total_particles = primary_particle.number_collidable_particles[thread_id];
//         double number_selected_real = P_null * static_cast<double>(number_total_particles);
//         size_t number_selected = static_cast<size_t>(number_selected_real);
//         double Rand = pcg32_random_r();
//         if (Rand < (number_selected_real - number_selected)) {
//             number_selected++;
//         }
        
//         size_t particle_indx, indx_low, indx_high, indx_middle, end_indx;
//         double particle_location, incident_velocity[3], target_velocity[3], velocity_CM[3], speed_CM, energy_CM, interp_d,
//             sigma_v, incident_energy, del_E;
//         for (size_t part_num=0;part_num<number_selected;part_num++){
//             particle_indx = static_cast<size_t>(number_total_particles * pcg32_random_r()) * 4;
//             particle_location = primary_particle.phase_space[thread_id][particle_indx];
//             incident_velocity[0] = primary_particle.phase_space[thread_id][particle_indx+1];
//             incident_velocity[1] = primary_particle.phase_space[thread_id][particle_indx+2];
//             incident_velocity[2] = primary_particle.phase_space[thread_id][particle_indx+3];
//             target_particle.generate_maxwellian_velocity(target_velocity);
//             speed_CM = 0.0;
//             for (iter=0;iter<3;iter++){
//                 velocity_CM[iter] = incident_velocity[iter] - target_velocity[iter];
//                 speed_CM += velocity_CM[iter] * velocity_CM[iter];
//             }
//             energy_CM = speed_CM * 0.5 * this->reduced_mass / Constants::elementary_charge;
//             speed_CM = std::sqrt(speed_CM);
//             indx_low = 0;
//             indx_high = this->length_arrays-1;
//             if (energy_CM <= this->min_energy){
//                 interp_d = 0.0;
//             } else if (energy_CM >= this->max_energy) {
//                 indx_low = indx_high-1;
//                 interp_d = 1.0;
//             } else {
//                 while (indx_low != indx_high-1) {
//                     indx_middle = (indx_low + indx_high)/2;
//                     if (this->energy_array[indx_middle] < energy_CM) {
//                         indx_low = indx_middle;
//                     } else if (this->energy_array[indx_middle] > energy_CM) {
//                         indx_high = indx_middle;
//                     } else {
//                         indx_low = indx_middle;
//                         indx_high = indx_low + 1;
//                     }
//                 }
//                 interp_d = (energy_CM - this->energy_array[indx_low])/(this->energy_array[indx_high] - this->energy_array[indx_low]);
//             }
//             Rand = pcg32_random_r();
//             sigma_v = 0.0;
//             for (coll_idx=0;coll_idx<this->number_collisions;coll_idx++){
//                 if (energy_CM > this->energy_threshold[coll_idx]) {
//                     sigma_v += (this->sigma_array[coll_idx][indx_low] * (1.0 - interp_d) + this->sigma_array[coll_idx][indx_high] * interp_d) * speed_CM;
//                     if (Rand <= sigma_v/this->sigma_v_max) {
//                         total_collisions[coll_idx]++;
//                         incident_energy = incident_velocity[0]*incident_velocity[0] + incident_velocity[1]*incident_velocity[1] + incident_velocity[2]*incident_velocity[2];
//                         tot_incident_energy[coll_idx] += incident_energy;
//                         del_E = (energy_CM - this->energy_threshold[coll_idx]) * Constants::elementary_charge;
//                         switch (this->collision_type[coll_idx]) {
//                             case 1:
//                                 this->double_product_isotropic(primary_mass, target_mass, del_E, incident_velocity, target_velocity);
//                                 break;
//                             case 2:
//                                 int secondary_product_idx = this->products_indx[coll_idx][1];
//                                 int third_product_idx = this->products_indx[coll_idx][2];
//                                 Particle& secondary_particle = particle_list[secondary_product_idx];
//                                 Particle& third_particle = particle_list[third_product_idx];
//                                 this->triple_product_isotropic(primary_mass, secondary_particle.mass, target_mass, del_E, 
//                                     incident_velocity, target_velocity, velocity_CM);
//                                 third_particle.number_particles[thread_id]++;
//                                 secondary_particle.number_particles[thread_id]++;

//                                 // electron set
//                                 size_t third_idx = third_particle.number_particles[thread_id]-1;
//                                 size_t secondary_idx = secondary_particle.number_particles[thread_id]-1;
//                                 third_particle.phase_space[thread_id][third_idx*4] = particle_location;
//                                 for (iter=0;iter<3;iter++){
//                                     third_particle.phase_space[thread_id][third_idx*4+iter+1] = velocity_CM[iter];
//                                 }

//                                 // ion set
//                                 secondary_particle.phase_space[thread_id][secondary_idx*4] = particle_location;
//                                 for (iter=0;iter<3;iter++){
//                                     secondary_particle.phase_space[thread_id][secondary_idx*4+iter+1] = target_velocity[iter];
//                                 }
//                                 energy_loss[coll_idx] += (- third_particle.mass * (velocity_CM[0]*velocity_CM[0] +
//                                     velocity_CM[1]*velocity_CM[1] + velocity_CM[2]*velocity_CM[2]) - secondary_particle.mass * (target_velocity[0]*target_velocity[0] + 
//                                     target_velocity[1]*target_velocity[1] + target_velocity[2]*target_velocity[2]));
//                                 break;
//                             case 3:
//                                 this->double_product_isotropic(primary_mass, target_mass, del_E, incident_velocity, target_velocity);
//                                 break;
//                             case 4:
//                                 incident_velocity[0] = target_velocity[0];
//                                 incident_velocity[1] = target_velocity[1];
//                                 incident_velocity[2] = target_velocity[2];
//                                 break;
//                         }
                        
//                         energy_loss[coll_idx] += primary_mass * (incident_energy - incident_velocity[0]*incident_velocity[0] - 
//                             incident_velocity[1]*incident_velocity[1] - incident_velocity[2]*incident_velocity[2]);
//                         break;
//                     }
//                 }
//             }
//             // switch particles location
//             end_indx = (number_total_particles-1)*4;
//             for (iter=0;iter<4;iter++)  {
//                 primary_particle.phase_space[thread_id][particle_indx+iter] = primary_particle.phase_space[thread_id][end_indx+iter];
//             }
//             primary_particle.phase_space[thread_id][end_indx] = particle_location;
//             for (iter=0;iter<3;iter++)  {
//                 primary_particle.phase_space[thread_id][end_indx+iter+1] = incident_velocity[iter];
//             }
//             number_total_particles--;
//         }
//         primary_particle.number_collidable_particles[thread_id] = number_total_particles;
//         #pragma omp critical
//         {
//            for (int iter=0;iter<this->number_collisions;iter++){
//                 this->total_incident_energy[iter] += tot_incident_energy[iter];
//                 this->total_energy_loss[iter] += energy_loss[iter];
//                 this->total_amount_collisions[iter] += total_collisions[iter];
//             } 
//         }
//     }
    
  

// }



// inline void Null_Collision::double_product_isotropic(const double &primary_mass, const double &target_mass, const double &del_E, double (&incident_velocity)[3], double (&target_velocity)[3]) {
    
//     double e_vector[3];
//     double cos_theta, speed_per_particle, phi, sin_theta, cos_phi, sin_phi, P_beginning, V_cm;
//     int i;

//     speed_per_particle = std::sqrt(2.0 * del_E * this->reduced_mass / primary_mass/primary_mass);

//     cos_theta = 1.0 - 2.0 * pcg32_random_r();
//     phi = pcg32_random_r() * 2.0 * M_PI;
//     sin_theta = std::sqrt(1.0 - cos_theta*cos_theta);
//     cos_phi = std::cos(phi);
//     sin_phi = std::sin(phi);
    

//     e_vector[0] = cos_phi * sin_theta;
//     e_vector[1] = sin_phi * sin_theta;
//     e_vector[2] = cos_theta;

//     for (i=0;i<3;i++){
//         P_beginning = primary_mass * incident_velocity[i] + target_mass * target_velocity[i];
//         V_cm = P_beginning/this->mass_sum;
//         incident_velocity[i] = e_vector[i] * speed_per_particle + V_cm;
//         target_velocity[i] = (P_beginning - primary_mass * incident_velocity[i])/target_mass;
//     }

// }

// inline void Null_Collision::triple_product_isotropic(const double &primary_mass, const double &ion_mass, const double &target_mass, const double &del_E, 
//     double (&incident_velocity)[3], double (&target_velocity)[3], double (&third_velocity)[3]) {
    
//     double e_vector[3], y_vector[3], u_vector[3];
//     double cos_theta, speed_per_particle, phi, sin_theta, cos_phi, sin_phi, cos_theta_new, sin_theta_new, speed_electron, P_beginning, V_cm;
//     int i;


//     speed_per_particle = std::sqrt(2.0 * del_E * this->reduced_mass_triple / primary_mass/primary_mass);

//     cos_theta = 1.0 - 2.0 * pcg32_random_r();
//     phi = pcg32_random_r() * 2.0 * M_PI;
//     sin_theta = std::sqrt(1.0 - cos_theta*cos_theta);
//     cos_phi = std::cos(phi);
//     sin_phi = std::sin(phi);
//     e_vector[0] = cos_phi * sin_theta;
//     e_vector[1] = sin_phi * sin_theta;
//     e_vector[2] = cos_theta;

//     cos_theta_new = Constants::cos_third_rot * cos_theta - Constants::sin_third_rot * sin_theta;
//     sin_theta_new = cos_theta * Constants::sin_third_rot + Constants::cos_third_rot * sin_theta;
//     y_vector[0] = cos_phi * sin_theta_new;
//     y_vector[1] = sin_phi * sin_theta_new;
//     y_vector[2] = cos_theta_new;

//     phi = pcg32_random_r() * 2.0 * M_PI;
//     cos_phi = std::cos(phi);
//     sin_phi = std::sin(phi);

//     u_vector[0] = e_vector[1] * y_vector[2] - e_vector[2] * y_vector[1];
//     u_vector[1] = e_vector[2] * y_vector[0] - e_vector[0] * y_vector[2];
//     u_vector[2] = e_vector[0] * y_vector[1] - e_vector[1] * y_vector[0];

//     speed_electron = std::sqrt(2.0 * del_E * this->reduced_mass_triple/ Constants::electron_mass/Constants::electron_mass);

//     for (i=0;i<3;i++){
//         u_vector[i] = -e_vector[i] * Constants::cos_third_rot * (cos_phi - 1.0) + y_vector[i] * cos_phi + u_vector[i] * sin_phi;
//         P_beginning = primary_mass * incident_velocity[i] + target_mass * target_velocity[i];
//         V_cm = P_beginning/this->mass_sum;
//         incident_velocity[i] = e_vector[i] * speed_per_particle + V_cm;
//         third_velocity[i] = u_vector[i] * speed_electron + V_cm;
//         target_velocity[i] = (P_beginning - primary_mass * incident_velocity[i] - Constants::electron_mass * third_velocity[i])/ion_mass;
//     }


// }

std::vector<null_collider> read_null_collision_inputs(const std::string& directory_path, const std::vector<charged_particle> &particle_list, const std::vector<target_particle> &target_particle_list){
    if (mpi_vars::mpi_rank==0){
        std::cout << " "<< std::endl;
        std::cout << "Reading null collision inputs "<< std::endl;
        std::cout << "---------------------------------------- "<< std::endl;
    }

    
    int number_charged_particles = particle_list.size();
    int number_target_particles = target_particle_list.size();
    std::vector<std::vector<int>> target_indices(number_charged_particles);
    std::vector<std::vector<std::vector<int>>> product_indices(number_charged_particles);
    std::vector<int> number_collisions_per_primary(number_charged_particles, 0);
    std::vector<std::vector<int>> collision_type_per_primary(number_charged_particles);
    std::vector<std::vector<double>> E_threshold_per_primary_collision(number_charged_particles);
    std::vector<std::vector<std::vector<double>>> energy_arrays(number_charged_particles);
    std::vector<std::vector<std::vector<double>>> sigma_arrays(number_charged_particles);
    double E_threshold, temp_var, E_scaling, sigma_scaling;

    for (int rank_num = 0;rank_num < mpi_vars::mpi_size; rank_num++) {
        if (rank_num == mpi_vars::mpi_rank){
            // Open the directory
            DIR* dir = opendir(directory_path.c_str());
            if (!dir) {
                perror("opendir");
                exit(EXIT_FAILURE);
            }

            struct dirent* entry;
            if ((entry = readdir(dir)) == nullptr && mpi_vars::mpi_rank == 0) {
                std::cout << "No binary collisions included!" << std::endl;
            }
            while ((entry = readdir(dir)) != nullptr) {
                if (entry->d_type == DT_REG) {  // regular file
                    std::string filename = directory_path + entry->d_name;
                    std::ifstream file(filename);
                    if (!file) {
                        std::cerr << "Error: Unable to open file " << filename << std::endl;
                    }
                    std::string line;
                    std::getline(file, line);
                    std::istringstream iss(line);
                    while (line.find("END") == std::string::npos) {
                        if (line.find("REACTION") != std::string::npos) {
                            std::string reaction_string;
                            std::getline(file, line);
                            iss.clear();
                            iss.str(line);
                            iss >> reaction_string;
                            if (mpi_vars::mpi_rank==0) {std::cout << reaction_string << std::endl;}
                            size_t arrow_pos = reaction_string.find("->");
                            std::string reactant_string = reaction_string.substr(0,arrow_pos);
                            std::string product_string = reaction_string.substr(arrow_pos+2);
                            std::regex pattern(R"(\[([^\]]+)\])");
                            std::smatch match;
                            int primary_idx, target_idx; 
                            int i;
                            // Get indices of the primary and target particle
                            while (std::regex_search(reactant_string, match, pattern)){
                                
                                for (i=0;i<number_charged_particles;i++){
                                    if (particle_list[i].name == match[1]){
                                        primary_idx = i;     
                                        break;
                                    }
                                }
                                for (i=0;i<number_target_particles;i++){
                                    if (target_particle_list[i].name == match[1]){
                                        target_idx = i;
                                        break;
                                    }
                                }
                                reactant_string = match.suffix().str();
                            }
                        
                            
                            target_indices[primary_idx].push_back(target_idx);
                            number_collisions_per_primary[primary_idx]++;

                            std::vector<int> local_product_indx;
                            target_idx = -1;
                            while (std::regex_search(product_string, match, pattern)){
                                for (i=0;i<number_charged_particles;i++){
                                    if (particle_list[i].name == match[1]){
                                        local_product_indx.push_back(i);
                                    }
                                }
                                for (i=0;i<number_target_particles;i++){
                                    if (target_particle_list[i].name == match[1]){
                                        target_idx = i;
                                    }
                                }
                                product_string = match.suffix().str();
                            }
                            if (target_idx > - 1) {
                                local_product_indx.push_back(target_idx); // put target index as last
                            }
                            
                            std::getline(file, line);
                            iss.clear();
                            iss.str(line);
                            iss >> E_threshold >> temp_var;
                            std::getline(file, line);
                            iss.clear();
                            iss.str(line);
                            iss >> E_scaling >> sigma_scaling;
                            std::string coll_string;
                            std::getline(file, line);
                            iss.clear();
                            iss.str(line);
                            iss >> coll_string;
                            int collision_type = 0;
                            if (coll_string == "ELASTIC") {
                                collision_type = 1;
                            } else if (coll_string == "EXCITATION"){
                                collision_type = 3;
                            } else if (coll_string == "IONIZATION"){
                                collision_type = 2;
                                // make sure primary particle is first, then followed by by-product
                                for (i=0;i<3;i++){
                                    if (primary_idx == local_product_indx[i]){
                                        // switch primary to front
                                        local_product_indx[3] = local_product_indx[0];
                                        local_product_indx[0] = primary_idx;
                                        local_product_indx[i] = local_product_indx[3];
                                        break;
                                    }
                                }
                                for (i=1;i<3;i++){
                                    if (particle_list[local_product_indx[i]].mass == constants::electron_mass){
                                        // switch electron to back
                                        local_product_indx[3] = local_product_indx[2];
                                        local_product_indx[2] = local_product_indx[i];
                                        local_product_indx[i] = local_product_indx[3];
                                        break;
                                    }
                                }
                                local_product_indx.pop_back();
                            } else if (coll_string == "CHARGEEXCHANGE"){
                                collision_type = 4;
                            }
                            collision_type_per_primary[primary_idx].push_back(collision_type);
                            E_threshold_per_primary_collision[primary_idx].push_back(E_threshold);
                            product_indices[primary_idx].push_back(local_product_indx);
                            
                            
                            std::getline(file, line);
                            while (line.find("---------") == std::string::npos) {
                                std::getline(file, line);
                            }
                            double energy_point, sigma_point;
                            std::vector<double> local_energy_array;
                            std::vector<double> local_sigma_array;
                            std::getline(file, line);
                            while (line.find("---------") == std::string::npos) {
                                iss.clear();
                                iss.str(line);
                                iss >> energy_point >> sigma_point;
                                local_energy_array.push_back(energy_point * E_scaling);
                                local_sigma_array.push_back(sigma_point * sigma_scaling);
                                std::getline(file, line);
                            }
                            // find which particle indices they match to
                            energy_arrays[primary_idx].push_back(local_energy_array);
                            sigma_arrays[primary_idx].push_back(local_sigma_array);
                
                        }

                        // go to next line
                        std::getline(file, line);
                    }
                    file.close();
                    
                }
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    


    std::vector<null_collider> binary_collision_list;
    binary_collision_list.reserve(number_charged_particles);
    for (int k=0; k<number_charged_particles; k++){
        binary_collision_list.emplace_back(number_collisions_per_primary[k], k, target_indices[k], 
            sigma_arrays[k], energy_arrays[k], E_threshold_per_primary_collision[k],
            collision_type_per_primary[k], product_indices[k]);
    }
    
    return binary_collision_list;

}