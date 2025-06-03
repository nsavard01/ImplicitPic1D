
#include <vector>
#include <omp.h>
#include "ES_solvers/ES_solver_MC.hpp"
#include "solvers/poisson_solver_1D_tridiag.hpp"

ES_solver_MC::ES_solver_MC(const domain& world) {
    this->phi.resize(world.number_nodes, 0.0);
    this->rho.resize(world.number_nodes, 0.0);
    this->E_field.resize(world.number_nodes, 0.0);
    int number_threads = omp_get_max_threads();
    this->work_space.resize(number_threads);
    for (int i = 0; i < number_threads; i++) {
        this->work_space[i].resize(world.number_nodes, 0.0);
    }
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

void ES_solver_MC::push_particles(int thread_id, double del_t, std::vector<charged_particle>& particle_list, const domain& world) {
    // Loop over all particles and push them to the grid
    int total_thread_count = omp_get_max_threads();
    int num_particles = particle_list.size();
    double inv_dx = 1.0 / world.min_dx; // Cell size
    int left_boundary = world.left_boundary_condition; // Get left boundary condition
    int right_boundary = world.right_boundary_condition; // Get right boundary condition
    int number_cells = world.number_cells; // Number of cells in the domain
    
    for (int i = 0; i < num_particles; ++i) {
        charged_particle& particle = particle_list[i];
        particle.ES_push_MC(thread_id, del_t, this->E_field, inv_dx, left_boundary, right_boundary, number_cells); // Push particles to the grid
    }

    #pragma omp barrier
    #pragma omp master
    {
        for (int part_num = 0; part_num < num_particles; part_num++){
            charged_particle& particle = particle_list[part_num];
            for (int i_thread = 0; i_thread < total_thread_count; i_thread++){
                particle.accum_wall_loss[0] += particle.wall_loss[i_thread][0];
                particle.accum_wall_loss[1] += particle.wall_loss[i_thread][1];
                particle.accum_wall_energy_loss[0] += particle.energy_loss[i_thread][0];
                particle.accum_wall_energy_loss[1] += particle.energy_loss[i_thread][1];
                for (int coord = 0; coord < particle.number_velocity_coordinates; coord++) {
                    particle.accum_wall_momentum_loss[0][coord] += particle.momentum_loss[i_thread][0][coord];
                    particle.accum_wall_momentum_loss[1][coord] += particle.momentum_loss[i_thread][1][coord];
                }
            }
        }
    }
    #pragma omp barrier
}



