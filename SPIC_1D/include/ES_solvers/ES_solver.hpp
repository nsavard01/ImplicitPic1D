
#pragma once
#include <vector>
#include <cmath>
#include "solvers/poisson_solver_1D.hpp"
#include "domain/domain.hpp"
#include "domain/uniform_domain.hpp"
#include "domain/non_uniform_domain.hpp"
#include "globals/mpi_vars.hpp"
#include "particles/charged_particle.hpp"
#include <cmath>

class ES_solver {

protected:
    std::vector<double> phi, rho, E_field; // diagonal and upper matrix elements
    double RF_half_amplitude, RF_rad_frequency, left_voltage, right_voltage; // RF amplitude and frequency
    int RF_indx = -1;
    std::unique_ptr<poisson_solver_1D> poisson_solver; // pointer to the Poisson solver
public:

    
    virtual ~ES_solver() = default;
    
    std::vector<double>& get_phi() {
        return this->phi;
    };

    std::vector<double>& get_rho() {
        return this->rho;
    };

    void set_phi(double left_voltage, double right_voltage, double RF_frequency, int left_boundary, int right_boundary);

    virtual void print_out() = 0;
    virtual void deposit_charge_density(std::vector<charged_particle>& particle_list);

};

std::unique_ptr<ES_solver> read_voltage_inputs(const std::string& filename, int scheme_type, const domain& world);

