
#include <vector>
#include <omp.h>
#include "ES_solvers/ES_solver_EC.hpp"
#include "solvers/poisson_solver_1D_tridiag.hpp"
#include "globals/mpi_vars.hpp"
#include "globals/constants.hpp"

ES_solver_EC::ES_solver_EC(const domain& world) {
    this->phi.resize(world.get_number_nodes(), 0.0);
    this->rho.resize(world.get_number_nodes(), 0.0);
    this->E_field.resize(world.get_number_cells(), 0.0);
    this->poisson_solver = std::make_unique<poisson_solver_1D_tridiag>(world);
}

void ES_solver_EC::print_out() {
    if (mpi_vars::mpi_rank == 0) {
        std::cout << "ES_solver_EC: " << std::endl;
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


void ES_solver_EC::make_EField(const domain& world) {
    // Calculate the electric field from the potential
    int number_cells = world.number_cells; // Number of cells in the domain
    if (world.domain_type == 0) {
        double inv_dx = 1.0/world.min_dx; // Cell size for uniform domain
        inv_dx = 1.0/world.min_dx; // Cell size for uniform domain
        for (int i = 0; i < number_cells; ++i) {
            this->E_field[i] = (this->phi[i] - this->phi[i+1]) * inv_dx; // Electric field calculation
        }
    } else if (world.domain_type == 1) {
        const std::vector<double>& dx = world.dx_dxi; // Cell size for non-uniform domain
        for (int i = 0; i < number_cells; ++i) {
            this->E_field[i] = (this->phi[i] - this->phi[i+1])/ dx[i]; // Electric field calculation
        }
    }

    if (mpi_vars::mpi_rank == 0) {
        for (int i = 0; i < number_cells; ++i) {
            std::cout << "EField[" << i << "] = " << this->E_field[i] << "should be " << -1e14 * constants::elementary_charge * (world.cell_centers[i] - 0.5 * world.length_domain)  / constants::epsilon_0 << std::endl; // Print charge density for debugging
        }
    } 
}



