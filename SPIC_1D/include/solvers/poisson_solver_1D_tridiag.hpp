
#pragma once
#include "solvers/poisson_solver_1D.hpp"
#include "domain/uniform_domain.hpp"
#include "domain/non_uniform_domain.hpp"
#include <string>
#include <memory>

class poisson_solver_1D_tridiag : public poisson_solver_1D {
protected:
    std::vector<double> diagonal, upper, lower; // diagonal and upper matrix elements
public:
    poisson_solver_1D_tridiag(const domain& world);
    void solve() override;
};

