
#include <vector>
#include <omp.h>
#include "ES_solvers/ES_solver_INGP.hpp"
#include "non_linear_solvers/non_linear_solver.hpp"
#include "solvers/poisson_solver_1D_tridiag.hpp"
#include "globals/mpi_vars.hpp"
#include "globals/constants.hpp"

ES_solver_INGP::ES_solver_INGP(const domain& world) {
    this->phi.resize(world.number_nodes, 0.0);
    this->phi_past.resize(world.number_nodes, 0.0);
    this->rho.resize(world.number_nodes, 0.0);
    this->J.resize(world.number_cells, 0.0);
    this->E_field.resize(world.number_cells, 0.0);
    int number_threads = omp_get_max_threads();
    this->work_space.resize(number_threads);
    for (int i = 0; i < number_threads; i++) {
        this->work_space[i].resize(world.number_nodes, 0.0);
    }
    this->poisson_solver = std::make_unique<poisson_solver_1D_tridiag>(world);
    std::vector<double> double_params(3);
    std::vector<int> int_params(3);
    read_non_linear_solver_inputs("../inputs/implicit_solver.inp", int_params, double_params); // Read non-linear solver inputs
    MPI_Bcast(double_params.data(), 3, MPI_DOUBLE, 0, MPI_COMM_WORLD); // Broadcast double parameters
    MPI_Bcast(int_params.data(), 3, MPI_INT, 0, MPI_COMM_WORLD); // Broadcast integer parameters
    this->implicit_solver = std::make_unique<AA_solver>(double_params[0], double_params[1], double_params[3], int_params[2], int_params[1], world.number_nodes);
    this->implicit_solver->print_out(); // Print implicit solver parameters
}

void ES_solver_INGP::print_out() {
    if (mpi_vars::mpi_rank == 0) {
        std::cout << "ES_solver_INGP: " << std::endl;
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


void ES_solver_INGP::make_EField(const domain& world) {
    // Calculate the electric field from the potential
    int number_cells = world.number_cells; // Number of cells in the domain
    if (world.domain_type == 0) {
        double inv_dx = 1.0/world.min_dx; // Cell size for uniform domain
        inv_dx = 1.0/world.min_dx; // Cell size for uniform domain
        for (int i = 0; i < number_cells; ++i) {
            this->E_field[i] = 0.5 * (this->phi[i] + this->phi_past[i] - this->phi[i+1] - this->phi_past[i+1]) * inv_dx; // Electric field calculation
        }
    } else if (world.domain_type == 1) {
        const std::vector<double>& dx = world.dx_dxi; // Cell size for non-uniform domain
        for (int i = 0; i < number_cells; ++i) {
            this->E_field[i] = 0.5 * (this->phi[i] + this->phi_past[i] - this->phi[i+1] - this->phi_past[i+1])/ dx[i]; // Electric field calculation
        }
    }

}




void ES_solver_INGP::push_particles(const int thread_id, double del_t, std::vector<charged_particle>& particle_list, const domain& world){
    // Loop over all particles and push them to the grid
    int num_particles = particle_list.size();
    int left_boundary = world.left_boundary_condition; // Get left boundary condition
    int right_boundary = world.right_boundary_condition; // Get right boundary condition
    int number_cells = world.number_cells; // Number of cells in the domain
    int domain_type = world.domain_type;
    
    std::vector<double>& local_work_space = this->work_space[thread_id];
    std::vector<double>& part_work_space = charged_particle::xi_sorted[thread_id];
    std::fill(local_work_space.begin(), local_work_space.end(), 0.0); // Reset work space for this thread
    if (domain_type == 0) {
        double inv_dx = 1.0 / world.min_dx; // Cell size
        for (int i = 0; i < num_particles; ++i) {
            charged_particle& particle = particle_list[i];
            std::fill(part_work_space.begin(), part_work_space.begin() + world.number_nodes, 0.0); // Reset work space for this particle
            particle.ES_push_deposit_INGP_uniform(thread_id, del_t, this->E_field, part_work_space, 
                inv_dx, left_boundary, right_boundary, number_cells); // Push particles to the grid
            for (int i = 0; i < world.number_nodes; i++) {
                local_work_space[i] += part_work_space[i] * particle.q_times_wp;
            }
        }
    } else {
        for (int i = 0; i < num_particles; ++i) {
            charged_particle& particle = particle_list[i];
            std::fill(part_work_space.begin(), part_work_space.begin() + world.number_nodes, 0.0); // Reset work space for this particle
            particle.ES_push_deposit_INGP_non_uniform(thread_id, del_t, this->E_field, part_work_space, 
                world.dx_dxi, left_boundary, right_boundary, number_cells); // Push particles to the grid
            for (int i = 0; i < world.number_nodes; i++) {
                local_work_space[i] += part_work_space[i] * particle.q_times_wp;
            }
        }
    }
}

void ES_solver_INGP::integrate_time_step(const int thread_id, double del_t, double current_time, const domain& world, std::vector<charged_particle>& particle_list) {

    auto integral_function = [&](std::vector<double>& res_output) {
        int total_thread_count = omp_get_max_threads();
        double part_timer_start;
        #pragma omp master
        {   
            this->make_EField(world);
            part_timer_start = MPI_Wtime();
        }
        this->push_particles(thread_id, del_t, particle_list, world);
        #pragma omp barrier
        #pragma omp for
        for (int i = 0; i < world.number_nodes; i++) {
            double sum = 0.0;
            for (int i_thread = 0; i_thread < total_thread_count; i_thread++) {
                sum += this->work_space[i_thread][i];
            }
            this->rho[i] = sum; // Set charge density for each cell
        }
        #pragma omp barrier
        #pragma omp master
        {
            MPI_Allreduce(MPI_IN_PLACE, this->rho.data(), world.number_nodes, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); // Synchronize charge density across all processes
        }
        #pragma omp barrier
        #pragma omp master
        {
            double end_time = MPI_Wtime();
            this->particle_timer += (end_time - part_timer_start);
            double start_time = MPI_Wtime();
            this->solve_potential(current_time + del_t, world);
            end_time = MPI_Wtime();
            this->potential_timer += (end_time - start_time);
        }
    };
    #pragma omp master
    {
        this->phi_past = this->phi; // Copy current potential to future potential
        this->particle_timer = 0.0;
        this->potential_timer = 0.0;
        if (mpi_vars::mpi_rank == 0) {
            std::cout << "ES_solver_INGP: Integrating time step at time " << current_time << " with del_t = " << del_t << std::endl;
        }
    }
    this->implicit_solver->solve(this->phi, integral_function); // Solve the non-linear system
    #pragma omp master
    {
        this->make_EField(world); // Calculate the electric field from the potential
    }
    #pragma omp barrier
    int num_particles = particle_list.size();
    int left_boundary = world.left_boundary_condition; // Get left boundary condition
    int right_boundary = world.right_boundary_condition; // Get right boundary condition
    int number_cells = world.number_cells; // Number of cells in the domain
    int domain_type = world.domain_type;
    std::vector<int> number_sub_steps(num_particles, 0); // Initialize number of sub-steps for each particle
    // Final push of particles
    if (domain_type == 0) {
        double inv_dx = 1.0 / world.min_dx; // Cell size
        for (int i = 0; i < num_particles; ++i) {
            charged_particle& particle = particle_list[i];
            particle.ES_push_INGP_uniform(thread_id, del_t, this->E_field, number_sub_steps[i], 
                inv_dx, left_boundary, right_boundary, number_cells); // Push particles to the grid
        }
    } else {
        for (int i = 0; i < num_particles; ++i) {
            charged_particle& particle = particle_list[i];
            particle.ES_push_INGP_non_uniform(thread_id, del_t, this->E_field, number_sub_steps[i], 
                world.dx_dxi, left_boundary, right_boundary, number_cells); // Push particles to the grid
        }
    }
    MPI_Abort(MPI_COMM_WORLD, 1); // Ensure all processes have completed the push operation
    
    

}





