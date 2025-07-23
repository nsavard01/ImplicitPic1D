
#include "charged_particle_operators/charged_particle_wall_injector.hpp"
#include "rand_gen/maxwell_generator.hpp"


charged_particle_wall_injector::charged_particle_wall_injector(std::vector<int>& particle_indx, std::vector<double>& current_density, std::vector<std::vector<double>>& v_3D, 
    std::vector<double>& v_therm, std::vector<double>& particle_location) {
    this->particle_indx = particle_indx;
    this->current_density = current_density;
    this->v_therm = v_therm; 
    this->v_3D = v_3D;
    int total_thread_count = omp_get_max_threads() * mpi_vars::mpi_size;
    this->particle_location = particle_location;
    for (int i = 0; i < this->current_density.size(); i++) {
        // Divide current density by number of threads, so have equal amount of particle introduced in each thread
        this->current_density[i] = this->current_density[i]/double(total_thread_count);
    }
}

void charged_particle_wall_injector::print_out() {
    if (mpi_vars::mpi_rank == 0) {
        std::cout  << std::endl;
        std::cout << " ------------ " << std::endl;
        std::cout << "Particle Injection Operator " << std::endl;
        for (int i = 0; i < this->particle_indx.size(); i++) {
            std::cout << "Particle index " << this->particle_indx[i] << std::endl;
            std::cout << "With current density: " << this->current_density[i] << " A/m^2 per thread" << std::endl;
            std::cout << "Particle location: " << this->particle_location[i] << std::endl;
            std::cout << "thermal velocity: " << this->v_therm[i] << std::endl;
            std::cout << "v_x: " << v_3D[i][0] <<  ", v_y: " << v_3D[i][1] << " , v_z: " << v_3D[i][2] << std::endl;
            std::cout << std::endl;
        }
        std::cout << " ------------ " << std::endl;
        std::cout  << std::endl;
    }
}


void charged_particle_wall_injector::run(const int thread_id, const double current_time, const double del_t, std::vector<charged_particle>& particle_list, const domain& world) {
    for (int i = 0; i < this->particle_indx.size(); i++){
        int part_indx = this->particle_indx[i];
        charged_particle& part = particle_list[part_indx];
        const double particle_position = this->particle_location[i];
        const double v_x_drift = this->v_3D[i][0];
        const double v_y_drift = this->v_3D[i][1];
        const double v_z_drift = this->v_3D[i][2];
        const double v_therm_local = this->v_therm[i] * ((v_x_drift>0) - (v_x_drift<0));
        double number_particles_inject = std::abs(del_t * this->current_density[i]/part.q_times_wp); 
        size_t number_selected = size_t(number_particles_inject); 
        if (pcg32_random_r() < (number_particles_inject - number_selected)) {
            number_selected++;
        }
        size_t number_particles = part.number_particles[thread_id][0];
        // references to particle arrays
        std::vector<double>& xi_local = part.xi[thread_id];
        std::vector<double>& v_x_local = part.v_x[thread_id];
        bool use_vy = (part.number_velocity_coordinates > 1);
        bool use_vz = (part.number_velocity_coordinates > 2);
        double v_x_temp, v_y_temp, v_z_temp;
        for (size_t part_indx = 0; part_indx < number_selected; part_indx++){
            maxwellian_3D_flux(v_x_temp, v_y_temp, v_z_temp, v_therm_local, v_x_drift);
            xi_local[number_particles] = particle_position;
            v_x_local[number_particles] = v_x_temp;
            if (use_vy) {
                part.v_y[thread_id][number_particles] = v_y_temp + v_y_drift;
            }
            if (use_vz) {
                part.v_z[thread_id][number_particles] = v_z_temp + v_z_drift;
            }
            number_particles++;
        }
        part.number_particles[thread_id][0] = number_particles;
        
    }
    

}


