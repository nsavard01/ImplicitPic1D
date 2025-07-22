
#include "charged_particle_operators/charged_particle_wall_injector.hpp"


charged_particle_wall_injector::charged_particle_wall_injector(std::vector<int>& particle_indx, std::vector<double>& current_density, std::vector<double>& v_x, std::vector<int>& wall_node_location) {
    this->particle_indx = particle_indx;
    this->current_density = current_density;
    this->v_x = v_x;
    this->wall_node_location = wall_node_location;
    for (int i = 0; i < this->v_x.size(); i++) {
        this->direction.push_back((v_x[i] > 0.0) - (v_x[i] < 0.0));
    }
}

void charged_particle_wall_injector::print_out() {
    if (mpi_vars::mpi_rank == 0) {
        std::cout  << std::endl;
        std::cout << " ------------ " << std::endl;
        std::cout << "Particle Injection Operator " << std::endl;
        for (int i = 0; i < this->particle_indx.size(); i++) {
            std::cout << "Particle index " << this->particle_indx[i] << std::endl;
            std::cout << "With current density: " << this->current_density[i] << " A/m^2, v_x: " << this->v_x[i] << std::endl;
            std::cout << "On boundary number: " << this->wall_node_location[i] << std::endl;
            std::cout << std::endl;
        }
        std::cout << " ------------ " << std::endl;
        std::cout  << std::endl;
    }
}


void charged_particle_wall_injector::run(const int thread_id, const double current_time, const double del_t, std::vector<charged_particle>& particle_list, const domain& world) {

}


