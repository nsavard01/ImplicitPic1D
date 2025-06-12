#pragma once
#include <vector>
#include <cmath>
#include <mkl.h>
#include <stdexcept>
#include <iostream>
#include <omp.h>
#include "globals/mpi_vars.hpp"
#include <iostream>

inline std::vector<double> solveNormalEquationMKL(const std::vector<double>& A, const std::vector<double>& b, int m, int n) {
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

inline std::vector<double> solveNormalEquationManual(
    const std::vector<double>& A,
    const std::vector<double>& b,
    int m, int n)
{
    // if (A.size() != size_t(m * n) || b.size() != size_t(m)) {
    //     throw std::invalid_argument("Matrix/vector size mismatch.");
    // }

    // Step 1: Compute AtA = Aᵗ * A (n×n)
    std::vector<double> AtA(n * n, 0.0);
    for (int i = 0; i < n; ++i) {
        for (int j = i; j < n; ++j) {
            double sum = 0.0;
            for (int k = 0; k < m; ++k) {
                sum += A[k + i * m] * A[k + j * m];  // column-major access
            }
            AtA[i + j * n] = sum;
            if (i != j) AtA[j + i * n] = sum;  // Symmetric fill
        }
    }

    // Step 2: Compute Atb = Aᵗ * b (n)
    std::vector<double> Atb(n, 0.0);
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int k = 0; k < m; ++k) {
            sum += A[k + i * m] * b[k];
        }
        Atb[i] = sum;
    }

    // Step 3: Cholesky decomposition: AtA = L * Lᵗ
    std::vector<double> L(n * n, 0.0);  // Lower triangular matrix
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            double sum = AtA[i + j * n];
            for (int k = 0; k < j; ++k)
                sum -= L[i + k * n] * L[j + k * n];

            if (i == j) {
                if (sum <= 0.0) throw std::runtime_error("Matrix is not positive definite.");
                L[i + j * n] = std::sqrt(sum);
            } else {
                L[i + j * n] = sum / L[j + j * n];
            }
        }
    }

    // Step 4: Solve L y = Atb (forward substitution)
    std::vector<double> y(n, 0.0);
    for (int i = 0; i < n; ++i) {
        double sum = Atb[i];
        for (int j = 0; j < i; ++j) {
            sum -= L[i + j * n] * y[j];
        }
        y[i] = sum / L[i + i * n];
    }

    // Step 5: Solve Lᵗ x = y (back substitution)
    std::vector<double> x(n, 0.0);
    for (int i = n - 1; i >= 0; --i) {
        double sum = y[i];
        for (int j = i + 1; j < n; ++j) {
            sum -= L[j + i * n] * x[j];
        }
        x[i] = sum / L[i + i * n];
    }

    return x;
}

inline std::vector<double> solveNormalEquationGauss(
    const std::vector<double>& A,
    const std::vector<double>& b,
    int m, int n)
{
    // if (A.size() != size_t(m * n) || b.size() != size_t(m)) {
    //     throw std::invalid_argument("Matrix/vector size mismatch.");
    // }

    // Step 1: Compute AtA = Aᵗ * A (n×n) and Atb = Aᵗ * b (n)
    std::vector<double> AtA(n * n, 0.0);
    std::vector<double> Atb(n, 0.0);

    for (int i = 0; i < n; ++i) {
        for (int j = i; j < n; ++j) {
            double sum = 0.0;
            for (int k = 0; k < m; ++k) {
                sum += A[k + i * m] * A[k + j * m];
            }
            if (!std::isfinite(sum)) {
                std::cout << "error" << std::endl;
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            AtA[i + j * n] = sum;
            AtA[j + i * n] = sum;  // symmetric
        }

        double sum_b = 0.0;
        for (int k = 0; k < m; ++k) {
            sum_b += A[k + i * m] * b[k];
        }
        if (!std::isfinite(sum_b)) {
            std::cout << "error" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        Atb[i] = sum_b;
    }
    // Step 2: Gaussian elimination on AtA and Atb
    for (int k = 0; k < n - 1; ++k) {
        double pivot = AtA[k + k * n];
        if (!std::isfinite(pivot)) {
            std::cout << "error" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        if (pivot == 0.0) {
            throw std::runtime_error("Zero pivot encountered in Gaussian elimination.");
        }

        for (int i = k + 1; i < n; ++i) {
            double factor = AtA[i + k * n] / pivot;
            if (!std::isfinite(factor)) {
                std::cout << "error" << std::endl;
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            AtA[i + k * n] = factor;  // store L

            for (int j = k + 1; j < n; ++j) {
                if (!std::isfinite(factor * AtA[k + j * n])) {
                    std::cout << "error" << std::endl;
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
                AtA[i + j * n] -= factor * AtA[k + j * n];
            }
        }
    }
    
    // Step 3: Forward substitution (solve L * y = Atb)
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < i; ++j) {
            if (!std::isfinite(AtA[i + j * n] * Atb[j])) {
                std::cout << "error" << std::endl;
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            Atb[i] -= AtA[i + j * n] * Atb[j];
        }
    }
    // Step 4: Backward substitution (solve U * x = y)
    for (int i = n - 1; i >= 0; --i) {
        for (int j = i + 1; j < n; ++j) {
            if (!std::isfinite(AtA[i + j * n] * Atb[j])) {
                std::cout << "error" << std::endl;
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            Atb[i] -= AtA[i + j * n] * Atb[j];
        }

        double diag = AtA[i + i * n];
        if (diag == 0.0) {
            throw std::runtime_error("Zero diagonal encountered in back substitution.");
        }
        Atb[i] /= diag;
    }

    // Atb now contains the solution vector x
    return Atb;
}


