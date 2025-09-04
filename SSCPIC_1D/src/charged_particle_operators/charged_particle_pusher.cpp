
#include "charged_particle_operators/charged_particle_pusher.hpp"
#include "rand_gen/pcg_rng.hpp"
#include <iomanip>
#include <dirent.h>
#include <fstream>
#include <sstream>
#include <algorithm>

charged_particle_pusher::charged_particle_pusher() {
    
}

void charged_particle_pusher::set_diagnostic_vectors(const std::vector<charged_particle>& particle_list) {
    int number_threads = omp_get_max_threads();
    int size_particles = particle_list.size();
    this->thread_time.resize(number_threads, 0.0); // [thread]
    this->integration_time.resize(number_threads);// [thread, charged_particle_number]
    this->number_cell_crossings.resize(number_threads);
    this->number_time_steps.resize(number_threads);
    for (int i_thread = 0; i_thread < number_threads; i_thread++) {
        this->integration_time[i_thread].resize(size_particles, 0);
        this->number_cell_crossings[i_thread].resize(size_particles, 0);
        this->number_time_steps[i_thread].resize(size_particles, 0);
    }

    this->total_integration_time.resize(size_particles);// [charged_particle_number]
    this->total_number_cell_crossings.resize(size_particles);
    this->total_number_time_steps.resize(size_particles);
}

void charged_particle_pusher::gather_mpi() {
    const int number_charged_particles = this->total_integration_time.size();
    std::fill(this->total_integration_time.begin(), this->total_integration_time.end(), 0);
    std::fill(this->total_number_cell_crossings.begin(), this->total_number_cell_crossings.end(), 0);
    std::fill(this->total_number_time_steps.begin(), this->total_number_time_steps.end(), 0);
    for (int i_thread = 0; i_thread < omp_get_max_threads(); i_thread++) {
        for (int part_num = 0; part_num < number_charged_particles; part_num++) {
            this->total_integration_time[part_num] += this->integration_time[i_thread][part_num];
            this->total_number_cell_crossings[part_num] += this->number_cell_crossings[i_thread][part_num];
            this->total_number_time_steps[part_num] += this->number_time_steps[i_thread][part_num];
        }
    }
    double local_min_time = *std::min_element(this->thread_time.begin(), this->thread_time.end());
    double local_max_time = *std::max_element(this->thread_time.begin(), this->thread_time.end());
    MPI_Allreduce(MPI_IN_PLACE, this->total_integration_time.data(), number_charged_particles, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, this->total_number_cell_crossings.data(), number_charged_particles, mpi_vars::mpi_size_t_type, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, this->total_number_time_steps.data(), number_charged_particles, mpi_vars::mpi_size_t_type, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_min_time, &this->min_cpu_time, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&local_max_time, &this->max_cpu_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
}

void charged_particle_pusher::print_out(const std::vector<charged_particle>& particle_list) const {
    if (mpi_vars::mpi_rank ==0 ) {
        const int number_charged_particles = particle_list.size();
        std::cout << "Maximum cpu time: " << this->max_cpu_time << " seconds" << std::endl;
        std::cout << "Mminimum cpu time: " << this->min_cpu_time << " seconds" << std::endl;
        for (int part_num = 0; part_num < number_charged_particles; part_num++) {
            std::cout << std::endl;
            std::cout << "Particle #: " << particle_list[part_num].name << std::endl;
            size_t particle_number = particle_list[part_num].total_number_particles;
            std::cout << "Total integration time per particle: " << this->total_integration_time[part_num]/double(particle_number) << std::endl;
            std::cout << "Total number cell crossings per particle: " << double(this->total_number_cell_crossings[part_num])/double(particle_number) << std::endl;
            std::cout << "Total number time steps per particle: " << double(this->total_number_time_steps[part_num])/double(particle_number) << std::endl;
        }
    }
}


void charged_particle_pusher::push_particle_trajectories_non_uniform(const int thread_id, std::vector<charged_particle>& particle_list, std::vector<null_collider>& null_collider_list, const std::vector<target_particle>& target_particle_list, const non_uniform_domain& world, const std::vector<double>& E_field){
    int number_charged_particles = particle_list.size();
    // make vector of 
    bool push_continue = true;
    std::vector<double> particle_component(7);
    std::vector<size_t> last_particle_idx_vector(number_charged_particles), prev_last_particle_idx_vector(number_charged_particles);
    std::vector<size_t> start_particle_idx_vector(number_charged_particles, 0);
    // initialize amount of particles
    for (int part_num = 0; part_num < number_charged_particles; part_num++) {
        last_particle_idx_vector[part_num] = particle_list[part_num].number_particles[thread_id];
        prev_last_particle_idx_vector[part_num] = last_particle_idx_vector[part_num];
        this->integration_time[thread_id][part_num] = 0;
        this->number_cell_crossings[thread_id][part_num] = 0;
        this->number_time_steps[thread_id][part_num] = 0;
        particle_list[part_num].reset_diagnostics(thread_id);
    }
    double start_time_total = MPI_Wtime();
    while (push_continue) {
        for (int part_num = 0; part_num < number_charged_particles; part_num++) {
            charged_particle& particle_local = particle_list[part_num];
            size_t& last_particle_idx = last_particle_idx_vector[part_num];
            std::vector<std::vector<double>>& thread_particle_components = particle_local.particle_components[thread_id];
            const double q_over_m = particle_local.q_over_m;
            const int number_cells = world.number_cells;
            const int left_boundary = world.left_boundary_condition;
            const int right_boundary = world.right_boundary_condition;
            const std::vector<double>& dx_dxi = world.dx_dxi;
            const double simp_coeff_end = 1.0 / 6.0;
            const double simp_coeff_center = 2.0 / 3.0;
            double time_passed = 0;
            size_t num_cell_crossings = 0;
            size_t num_time_steps = 0;
            std::vector<double>& thread_density = particle_local.density_grid[thread_id];
            std::vector<std::vector<double>>& thread_v_sqr = particle_local.v_sqr_grid[thread_id];
            std::vector<std::vector<double>>& thread_v_drift = particle_local.v_drift_grid[thread_id];
            null_collider& particle_collider = null_collider_list[part_num];
            const std::vector<double> max_time_step = particle_collider.max_time_step;
            bool particle_collider_bool = (particle_collider.number_targets > 0);
            const size_t start_idx = start_particle_idx_vector[part_num];
            for (size_t part_idx = start_idx; part_idx < last_particle_idx; part_idx++){
                particle_component = thread_particle_components[part_idx];
                double& freq_rel = particle_component[0];
                double& xi_i = particle_component[1];
                double& v_x_i = particle_component[4];
                double& v_y = particle_component[5];
                double& v_z = particle_component[6];
                double v_x_f, xi_f, v_x_half;
                int xi_boundary;
                int cell_num = int(xi_i);
                double E_field_local = E_field[cell_num];
                double dx = dx_dxi[cell_num];
                double accel = q_over_m * E_field_local;
                double v_i_sqr = v_x_i*v_x_i;
                double v_f_sqr;
                double del_tau;
                double del_tau_boundary;
                double del_tau_collision = 1.0; // maximum allowable time
                double del_tau_collision_cell;
                if (particle_collider_bool) {
                    del_tau_collision_cell = max_time_step[cell_num];
                }
                double interp_xi, interp_num, interp_v_sqr, interp_v_x;
                int v_sign = (v_x_i > 0) - (v_x_i < 0);
                bool boundary_bool;
                bool wall_bool = false;
                while (! wall_bool) {
                    // get time to wall based on starting condition
                    num_time_steps += 1;
                    xi_boundary = cell_num + (v_sign + 1)/2;
                    xi_f = double(xi_boundary);
                    v_f_sqr = 2.0 * accel * (xi_f - xi_i) * dx + v_i_sqr;
                    if (v_f_sqr > 0) {
                        // minimum del_tau to boundary in direction of propogation
                        v_x_f = v_sign * std::sqrt(v_f_sqr);
                        if (std::abs(v_x_f - v_x_i) > std::abs(v_x_f + v_x_i)) {
                            del_tau_boundary = (v_x_f - v_x_i)/accel;
                        } else {
                            del_tau_boundary = 2.0 * (xi_f - xi_i)/ (v_x_f + v_x_i) * dx;
                        }
                    } else {
                        // minimum del_tau oppposite boundary
                        xi_boundary = 2*cell_num + 1 - xi_boundary;
                        xi_f = double(xi_boundary);
                        v_f_sqr = 2.0 * accel * (xi_f - xi_i) * dx + v_i_sqr;
                        v_x_f = - v_sign * std::sqrt(v_f_sqr);
                        del_tau_boundary = (v_x_f - v_x_i) / accel;
                    }

                    
                    if (particle_collider_bool) {
                        del_tau_collision = - std::log(pcg32_random_r()) * del_tau_collision_cell;
                    }
                    // If goes to boundary take del_tau boundary
                    boundary_bool = (del_tau_boundary < del_tau_collision);
                    if (boundary_bool) {
                        del_tau = del_tau_boundary;
                    } else {
                        del_tau = del_tau_collision;
                    }

                    // get properties at half time in trajectory
                    v_x_half = 0.5 * (v_x_i + v_x_f);
                    // xi_half = xi_i + 0.25 * (v_x_i + v_x_half) * del_tau / dx;
                    

                    // interpolate half way point to densities and v_sqr
                    interp_num = freq_rel * del_tau;
                    interp_xi = 0.5 * (xi_f + xi_i) - cell_num;
                    interp_v_sqr = interp_num * (simp_coeff_end * (v_i_sqr + v_f_sqr) + simp_coeff_center * v_x_half * v_x_half);
                    interp_v_x = interp_num * v_x_half;
                    thread_density[cell_num] += interp_num * (1.0 - interp_xi);
                    thread_density[cell_num+1] += interp_num * interp_xi;
                    
                    thread_v_sqr[cell_num][0] += interp_v_sqr * (1.0 - interp_xi);
                    thread_v_sqr[cell_num+1][0] += interp_v_sqr * (interp_xi);

                    thread_v_drift[cell_num][0] += interp_v_x * (1.0 - interp_xi);
                    thread_v_drift[cell_num+1][0] += interp_v_x * (interp_xi);


                    thread_v_drift[cell_num][1] += interp_num * v_y * (1.0 - interp_xi);
                    thread_v_drift[cell_num+1][1] += interp_num * v_y *  (interp_xi);

                    thread_v_drift[cell_num][2] += interp_num * v_z * (1.0 - interp_xi);
                    thread_v_drift[cell_num+1][2] += interp_num * v_z *  (interp_xi);

                    
                    thread_v_sqr[cell_num][1] += thread_v_drift[cell_num][1] * v_y;
                    thread_v_sqr[cell_num+1][1] += thread_v_drift[cell_num+1][1] * v_y;

                    thread_v_sqr[cell_num][2] += thread_v_drift[cell_num][2] * v_z;
                    thread_v_sqr[cell_num+1][2] += thread_v_drift[cell_num+1][2] * v_z;

                    



                    // reset for next trajectory
                    v_sign = (v_x_f > 0) - (v_x_f < 0);
                    if (boundary_bool) {
                        num_cell_crossings += 1;
                        cell_num = cell_num + v_sign;
                        if (particle_collider_bool) {
                            del_tau_collision_cell  = max_time_step[cell_num];
                        }
                        if (xi_boundary == 0) {
                            switch (left_boundary){
                                case 1:
                                case 4:
                                    cell_num = 0;
                                    wall_bool = true;
                                    particle_local.energy_loss[thread_id][0] += v_x_f*v_x_f + v_y*v_y + v_z*v_z;
                                    particle_local.wall_loss[thread_id][0]++;
                                    particle_local.wall_freq_rel[thread_id][0] += freq_rel;
                                    break;
                                case 2:
                                    cell_num = 0;
                                    v_x_f = - v_x_f;
                                    v_sign = -v_sign;
                                    break;
                                case 3:
                                    xi_f = double(number_cells);
                                    cell_num = number_cells-1;
                                    break;
                            }
                        } else if (xi_boundary == number_cells) {
                            switch (right_boundary){
                                case 1:
                                case 4:
                                    cell_num = number_cells-1;
                                    wall_bool = true;
                                    particle_local.energy_loss[thread_id][1] += v_x_f*v_x_f + v_y*v_y + v_z*v_z;
                                    particle_local.wall_loss[thread_id][1]++;
                                    particle_local.wall_freq_rel[thread_id][1] += freq_rel;
                                    break;
                                case 2:
                                    cell_num = number_cells-1;
                                    v_x_f = -v_x_f;
                                    v_sign = -v_sign;
                                    break;
                                case 3:
                                    xi_f = 0.0;
                                    cell_num = 0;
                                    break;
                            }        
                        }
                        E_field_local = E_field[cell_num];
                        dx = dx_dxi[cell_num];
                        accel = q_over_m * E_field_local;
                    } else {
                        // if doesn't reach boundary then test collision
                        if (particle_collider_bool) {
                            particle_collider.generate_null_collision(thread_id, cell_num, particle_component, particle_list, target_particle_list);
                        }
                    }
                    v_x_i = v_x_f;
                    xi_i = xi_f;
                    time_passed += del_tau;
                    v_i_sqr = v_x_i*v_x_i;
                }
            }
            // get final diagnostics
            this->integration_time[thread_id][part_num] += time_passed;
            this->number_cell_crossings[thread_id][part_num] += num_cell_crossings;
            this->number_time_steps[thread_id][part_num] += num_time_steps;
        }
        push_continue = false;
        for (int part_num = 0; part_num < number_charged_particles; part_num++) {
            if (last_particle_idx_vector[part_num] != prev_last_particle_idx_vector[part_num]) {
                // ionization has led to more particles, need to reset and re-loop
                start_particle_idx_vector[part_num] = last_particle_idx_vector[part_num];
                prev_last_particle_idx_vector[part_num] = last_particle_idx_vector[part_num];
                push_continue = true;
            }
        }
    }
    this->thread_time[thread_id] = MPI_Wtime() - start_time_total;
    
}

void charged_particle_pusher::push_particle_trajectories_uniform(const int thread_id, std::vector<charged_particle>& particle_list, std::vector<null_collider>& null_collider_list, const std::vector<target_particle>& target_particle_list, const uniform_domain& world, const std::vector<double>& E_field){
    // size_t last_particle_idx = this->number_particles[thread_id];
    // std::vector<std::vector<double>>& thread_particle_components = this->particle_components[thread_id];
    // const double q_over_m = this->q_over_m;
    // const double charge = this->charge;
    // for (size_t part_idx = 0; part_idx < last_particle_idx; part_idx++){
    //     double freq_rel = thread_particle_components[part_idx][0];
    //     double xi = thread_particle_components[part_idx][1];
    //     double v_x = thread_particle_components[part_idx][4];
    //     double v_y = thread_particle_components[part_idx][5];
    //     double v_z = thread_particle_components[part_idx][6];
    //     int cell_num = int(xi);
    //     double E_field_local = E_field[cell_num];

    //     bool wall_bool = false;
    //     while (! wall_bool) {
    //     }
    // }
}
