#pragma once
#include <vector>
#include <cmath>
#include <mkl.h>
#include <stdexcept>
#include <iostream>

std::vector<double> solveNormalEquationMKL(const std::vector<double>& A, const std::vector<double>& b, int m, int n) {
    // if (A.size() != size_t(m * n) || b.size() != size_t(m)) {
    //     throw std::invalid_argument("Matrix/vector size mismatch.");
    // }

    // Compute Aᵗ * A (size n×n)
    std::vector<double> AtA(n * n, 0.0);
    cblas_dsyrk(CblasColMajor, CblasUpper, CblasTrans, 
                n, m, 
                1.0, A.data(), m,
                0.0, AtA.data(), n);

    // Compute Aᵗ * b (size n)
    std::vector<double> Atb(n, 0.0);
    cblas_dgemv(CblasColMajor, CblasTrans,
                m, n,
                1.0, A.data(), m,
                b.data(), 1,
                0.0, Atb.data(), 1);

    // Solve (Aᵗ A) x = Aᵗ b using LAPACK's dposv (since AtA is symmetric positive-definite)
    lapack_int info = LAPACKE_dposv(LAPACK_COL_MAJOR, 'U', n, 1,
                                    AtA.data(), n,
                                    Atb.data(), n);

    if (info != 0) {
        throw std::runtime_error("LAPACKE_dposv failed with error code " + std::to_string(info));
    }

    // Atb now contains the solution vector x
    return Atb;
}


class AA_solver { // Anderson accelerated 
public:

    double beta, eps_r, eps_a;
    int max_iterations, m_anderson;
    int number_iterations, number_unknowns;
    std::vector<double> norm_residual, min_matrix;
    std::vector<std::vector<double>> residual_k, x_k;
    AA_solver(double beta, double eps_r, double eps_a, int max_iterations, int m_anderson, int number_unknowns);

    template <typename Func> // fixed point
    void solve(std::vector<double>& x_result, // x_result first with initial guess, then pass actual result
               Func&& fixed_point_function) { //Fixed point function returns next value x_k F(x_k, x_k+1), always uses x_result
        this->x_k[0] = x_result;
        fixed_point_function(x_result);
        this->x_k[1] = x_result;
        double sum_residual = 0.0;
        for (int i = 0; i < this->number_unknowns; i++){
            this->residual_k[0][i] = this->x_k[1][i] - this->x_k[0][i];
            sum_residual += this->residual_k[0][i] * this->residual_k[0][i];
        }
        this->norm_residual[0] = std::sqrt(sum_residual);
        double eps_tol = this->eps_r * this->norm_residual[0] + this->eps_a * std::sqrt(double(this->number_unknowns));
        int index, m_k;
        for (int iter = 1; iter< this->max_iterations; iter++) {
            if (iter < this->m_anderson) {
                m_k = iter;
            } else {
                m_k = this->m_anderson;
            }
            index = iter % (this->m_anderson+1);
            fixed_point_function(x_result); // fixed  
            sum_residual = 0.0;
            for (int i = 0; i < this->number_unknowns; i++){
                this->residual_k[index][i] = x_result[i] - this->x_k[index][i];
                sum_residual += this->residual_k[index][i] * this->residual_k[index][i];
            }
            this->norm_residual[index] = std::sqrt(sum_residual);
            if (this->norm_residual[index] < eps_tol) {
                this->number_iterations = iter + 1;
                break;
            }
            for (int j = 0; j < m_k; j++) {
                int past_indx = (index - m_k + j) % (this->m_anderson + 1);
                size_t start_indx = j * this->number_unknowns; // flattened index start
                for (int i = 0; i < this->number_unknowns; i++) {
                    this->min_matrix[start_indx+i] = this->residual_k[index][i] - this->residual_k[past_indx][i];
                }      
            }

            // solve minimization problem
            std::vector<double> alpha = solveNormalEquationMKL(min_matrix, residual_k[index], this->number_unknowns, m_k);
            
            int next_idx = (index+1) % (this->m_anderson + 1);
            int past_indx = (index - m_k) % (this->m_anderson + 1);
            double alpha_last = alpha[0];
            // initially set next x_k to first component
            for (int i = 0; i < this->number_unknowns; i++) {
                this->x_k[next_idx][i] = alpha[0] * (this->beta * residual_k[past_indx][i] + x_k[past_indx][i]);
            }   
            for (int j = 1; j < m_k; j++) {
                // Sum alphas
                alpha_last += alpha[j];
                past_indx = (index - m_k + j) % (this->m_anderson + 1);
                for (int i = 0; i < this->number_unknowns; i++) {
                    this->x_k[next_idx][i] += alpha[j] * (this->beta * residual_k[past_indx][i] + x_k[past_indx][i]);
                }   
            }
            alpha_last = 1.0 - alpha_last; // close coefficients so add to 1
            for (int i = 0; i < this->number_unknowns; i++) {
                this->x_k[next_idx][i] += alpha_last * (this->beta * residual_k[index][i] + x_k[index][i]); // add current component
                x_result[i] = this->x_k[next_idx][i];
            }

            

        }

    }
};


