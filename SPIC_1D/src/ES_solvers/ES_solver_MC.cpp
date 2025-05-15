
#include <vector>
#include <omp.h>
#include "ES_solvers/ES_solver_MC.hpp"
#include "solvers/poisson_solver_1D_tridiag.hpp"

ES_solver_MC::ES_solver_MC(const domain& world) {
    this->phi.resize(world.get_number_nodes(), 0.0);
    this->rho.resize(world.get_number_nodes(), 0.0);
    this->E_field.resize(world.get_number_nodes(), 0.0);
    this->poisson_solver = std::make_unique<poisson_solver_1D_tridiag>(world);
}

void ES_solver_MC::print_out() {
    if (mpi_vars::mpi_rank == 0) {
        std::cout << "ES_solver_MC: " << std::endl;
        std::cout << "-------------------------- " << std::endl;
        std::cout << "Number of phi nodes: " << this->phi.size() << std::endl;
        std::cout << "Number of field nodes: " << this->E_field.size() << std::endl;
        std::cout << "Left voltage: " << this->left_voltage << std::endl;
        std::cout << "Right voltage: " << this->right_voltage << std::endl;
        std::cout << "RF frequency: " << this->RF_rad_frequency / (2.0 * M_PI) << std::endl;
        if (this->RF_half_amplitude != 0.0) {
            std::cout << "RF half amplitude " << this->RF_half_amplitude << std::endl;
        } else {
            std::cout << "No RF set." << std::endl;
        }
        std::cout << "-------------------------- " << std::endl;
    }
}

// void ES_solver_MC::push_particles(std::vector<charged_particle>& particle_list) {
//     // Loop over all particles and push them to the grid
//     int total_thread_count = omp_get_max_threads();
//     int num_particles = particle_list.size();
//     #pragma omp parallel
//     {   
//         int thread_id = omp_get_thread_num();
//         for (int i = 0; i < num_particles; ++i) {
//             charged_particle& particle = particle_list[i];
//             particle.push_particles_linear(thread_id); // Push particles to the grid
//         }
//     }
// }



