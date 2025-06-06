#include "non_linear_solver/AA_solver.hpp"



AA_solver::AA_solver(double beta, double eps_r, double eps_a, int max_iterations, int m_anderson, int number_unknowns){
    this->beta = beta;
    this->eps_r = eps_r;
    this->eps_a = eps_a;
    this->max_iterations = max_iterations;
    this->m_anderson = m_anderson;
    this->number_unknowns = number_unknowns;
}

// std::vector<double> solveNormalEquationMKL(const std::vector<double>& A, const std::vector<double>& b, int m, int n) {
//     if (A.size() != size_t(m * n) || b.size() != size_t(m)) {
//         throw std::invalid_argument("Matrix/vector size mismatch.");
//     }

//     // Compute Aᵗ * A (size n×n)
//     std::vector<double> AtA(n * n, 0.0);
//     cblas_dsyrk(CblasColMajor, CblasUpper, CblasTrans, 
//                 n, m, 
//                 1.0, A.data(), m,
//                 0.0, AtA.data(), n);

//     // Compute Aᵗ * b (size n)
//     std::vector<double> Atb(n, 0.0);
//     cblas_dgemv(CblasColMajor, CblasTrans,
//                 m, n,
//                 1.0, A.data(), m,
//                 b.data(), 1,
//                 0.0, Atb.data(), 1);

//     // Solve (Aᵗ A) x = Aᵗ b using LAPACK's dposv (since AtA is symmetric positive-definite)
//     lapack_int info = LAPACKE_dposv(LAPACK_COL_MAJOR, 'U', n, 1,
//                                     AtA.data(), n,
//                                     Atb.data(), n);

//     if (info != 0) {
//         throw std::runtime_error("LAPACKE_dposv failed with error code " + std::to_string(info));
//     }

//     // Atb now contains the solution vector x
//     return Atb;
// }

