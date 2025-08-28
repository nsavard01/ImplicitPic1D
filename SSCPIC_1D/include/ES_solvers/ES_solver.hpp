
#pragma once
#include <vector>
#include <cmath>
#include "solvers/poisson_solver_1D.hpp"
#include "domain/domain.hpp"
#include "domain/uniform_domain.hpp"
#include "domain/non_uniform_domain.hpp"
#include "globals/mpi_vars.hpp"
#include "particles/target_particle.hpp"
#include <cmath>

class ES_solver {

public:
    std::vector<double> phi, rho, E_field; // diagonal and upper matrix elements
    double gauss_error = 0.0;
    int field_interpolation;
    double left_voltage, right_voltage; // RF amplitude and frequency
    std::unique_ptr<poisson_solver_1D> poisson_solver; // pointer to the Poisson solver

    // diagnostics
    double total_field_energy;
    double potential_timer;
    
    ES_solver(double left_voltage, double right_voltage, const domain& world);
    
    std::vector<double>& get_phi() {
        return this->phi;
    };

    std::vector<double>& get_rho() {
        return this->rho;
    };


    void print_out();
    // void write_particle_densities(const std::string file_path, const std::string filename, std::vector<charged_particle>& particle_list, const domain& world) const; // since density determinined by potential solver type
    void deposit_charge_density(const domain& world, std::vector<target_particle>& particle_list);
    // void deposit_density(std::vector<charged_particle>& particle_list, int thread_id);
    void solve_potential(const domain& world);
    // void solve_field_energy(const domain& world);
    void make_EField(const domain& world);
    // // general integration through time step which solves for fields after time step given initial fields
    // // general enough that it can include non-linear processes as well (implicit)
    // void get_diagnostics(const domain& world, std::vector<charged_particle>& particle_list);
    // void integrate_time_step(const int thread_id, double del_t, double current_time, const domain& world, std::vector<charged_particle>& particle_list) = 0; 
    // void push_particles(const int thread_id, double del_t, std::vector<charged_particle>& particle_list, const domain& world) = 0;

};

std::unique_ptr<ES_solver> read_voltage_inputs(const std::string& filename, const domain& world);

