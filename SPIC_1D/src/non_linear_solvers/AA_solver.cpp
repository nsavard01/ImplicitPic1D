#include "non_linear_solvers/AA_solver.hpp"
#include "globals/mpi_vars.hpp"



AA_solver::AA_solver(double beta, double eps_r, double eps_a, int max_iterations, int m_anderson, int number_unknowns){
    this->beta = beta;
    this->eps_r = eps_r;
    this->eps_a = eps_a;
    this->max_iterations = max_iterations;
    this->m_anderson = m_anderson;
    this->number_unknowns = number_unknowns;
    this->norm_residual.resize(this->m_anderson+1, 0.0);
    this->min_matrix.resize(this->number_unknowns * this->m_anderson);
    this->residual_k.resize(this->m_anderson + 1);
    this->x_k.resize(this->m_anderson + 1);
    for (int i = 0; i < this->m_anderson + 1; i++) {
        this->residual_k[i].resize(this->number_unknowns, 0.0);
        this->x_k[i].resize(this->number_unknowns, 0.0);
    }
}

void AA_solver::print_out() const {
    if (mpi_vars::mpi_rank == 0) {
        std::cout << "AA_solver: " << std::endl;
        std::cout << "-------------------------- " << std::endl;
        std::cout << "Beta: " << this->beta << std::endl;
        std::cout << "Epsilon_r: " << this->eps_r << std::endl;
        std::cout << "Epsilon_a: " << this->eps_a << std::endl;
        std::cout << "Max iterations: " << this->max_iterations << std::endl;
        std::cout << "M Anderson: " << this->m_anderson << std::endl;
        std::cout << "Number of unknowns: " << this->number_unknowns << std::endl;
        std::cout << "-------------------------- " << std::endl;
    }
}

