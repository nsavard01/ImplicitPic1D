#include "solvers/poisson_solver_1D.hpp"
#include "solvers/poisson_solver_1D_tridiag.hpp"
#include "globals/mpi_vars.hpp"

poisson_solver_1D_tridiag::poisson_solver_1D_tridiag(const domain& world) {
    this->number_unknowns = world.get_number_nodes();
    this->diagonal.resize(this->number_unknowns, 0.0);
    this->upper.resize(this->number_unknowns-1, 0.0);
    this->lower.resize(this->number_unknowns-1, 0.0);
    this->source_term.resize(this->number_unknowns, 0.0);
    this->solution.resize(this->number_unknowns, 0.0);
    if (typeid(world) == typeid(uniform_domain)) {
        double dx = world.get_min_dx();
        switch (world.get_left_boundary_condition()) {
            case 1:
            case 4:
            case 3:
                this->diagonal[0] = 1.0;
                this->upper[0] = 0.0;
                break;
            case 2:
                this->diagonal[0] = -1.0 / dx;
                this->upper[0] = 1.0 / dx;
                break;
            default:
                throw std::invalid_argument("Invalid left boundary condition.");
        }
    
        switch (world.get_right_boundary_condition()) {
            case 1:
            case 4:
            case 3:
                this->diagonal[this->number_unknowns-1] = 1.0;
                this->lower[this->number_unknowns-2] = 0.0;
                break;
            case 2:
                this->diagonal[this->number_unknowns-1] = -1.0 / dx;
                this->lower[this->number_unknowns-2] = 1.0 / dx;
                break;
            default:
                throw std::invalid_argument("Invalid right boundary condition.");
        }

        for (int i = 1; i < this->number_unknowns - 1; ++i) {
            this->diagonal[i] = -2.0 / dx;
            this->upper[i] = 1.0 / dx;
            this->lower[i-1] = 1.0 / dx;
        }
    } else if (typeid(world) == typeid(non_uniform_domain)) {
        const std::vector<double>& dx_dxi = world.get_dx_dxi();
        switch (world.get_left_boundary_condition()) {
            case 1:
            case 4:
            case 3:
                this->diagonal[0] = 1.0;
                this->upper[0] = 0.0;
                break;
            case 2:
                this->diagonal[0] = -1.0 / dx_dxi[0];
                this->upper[0] = 1.0 / dx_dxi[0];
                break;
            default:
                throw std::invalid_argument("Invalid left boundary condition.");
        }
    
        switch (world.get_right_boundary_condition()) {
            case 1:
            case 4:
            case 3:
                this->diagonal[this->number_unknowns-1] = 1.0;
                this->lower[this->number_unknowns-2] = 0.0;
                break;
            case 2:
                this->diagonal[this->number_unknowns-1] = -1.0 / dx_dxi[this->number_unknowns-2];
                this->lower[this->number_unknowns-2] = 1.0 / dx_dxi[this->number_unknowns-2];
                break;
            default:
                throw std::invalid_argument("Invalid right boundary condition.");
        }

        for (int i = 1; i < this->number_unknowns - 1; ++i) {
            this->diagonal[i] = -(1.0 / dx_dxi[i-1] + 1.0 / dx_dxi[i]);
            this->upper[i] = 1.0 / dx_dxi[i];
            this->lower[i-1] = 1.0 / dx_dxi[i-1];
        }
    }

    // if (mpi_vars::mpi_rank == 0) {
    //     for (int i = 0; i < this->number_unknowns; ++i) {
    //         std::cout << "Diagonal[" << i << "] = " << this->diagonal[i] << std::endl;
    //         if (i < this->number_unknowns - 1) {
    //             std::cout << "Upper[" << i << "] = " << this->upper[i] << std::endl;
    //             std::cout << "Lower[" << i << "] = " << this->lower[i] << std::endl;
    //         }
    //     }
    // }
    
    
}

void poisson_solver_1D_tridiag::solve() {
    // Solve the Poisson equation using the Thomas algorithm, with solution vector as source then turned into solution
    // use source_term as working vector
    double m;
    this->source_term[0] = this->upper[0] / this->diagonal[0];
    this->solution[0] = this->solution[0] / this->diagonal[0];
    for (int i = 1; i < this->number_unknowns-1; i++) {
        m = this->diagonal[i] - this->lower[i-1] * this->source_term[i-1];
        this->source_term[i] = this->upper[i] / m;
        this->solution[i] = (this->solution[i] - this->lower[i-1] * this->solution[i-1]) / m;
    }

    m = this->diagonal[this->number_unknowns-1] - this->lower[this->number_unknowns-2] * this->source_term[this->number_unknowns-2];
    this->solution[this->number_unknowns-1] = (this->solution[this->number_unknowns-1] - this->lower[this->number_unknowns-2] * this->solution[this->number_unknowns-2]) / m; 
    for (int i = this->number_unknowns-2; i >= 0; i--) {
        this->solution[i] = this->solution[i] - this->source_term[i] * this->solution[i+1];
    }
}


