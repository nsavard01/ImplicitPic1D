
#pragma once
#include <vector>
#include "solvers/poisson_solver_1D.hpp"
#include "particles/charged_particle.hpp"

class ES_solver {

protected:
    std::vector<std::vector<double>> work_space; // work space for openmp calculations
    std::vector<double> phi, rho; // diagonal and upper matrix elements
    double RF_half_amplitude, RF_rad_frequency; // RF amplitude and frequency
    std::unique_ptr<poisson_solver_1D> poisson_solver; // pointer to the Poisson solver
public:

    
    virtual ~ES_solver() = default;
    
    std::vector<double>& get_phi() {
        return this->phi;
    };

    std::vector<double>& get_rho() {
        return this->rho;
    };

    
};

