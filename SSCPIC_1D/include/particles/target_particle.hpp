#pragma once

#include <vector>
#include <cstddef>
#include <string>
#include "domain/domain.hpp"

class target_particle {
public:
    double charge, mass, average_density, average_temperature, accum_energy_change, diffusion_coeff;
    std::vector<double> density;
    int charged_particle_idx = -1;
    std::vector<std::vector<double>> v_therm_sqr, v_drift;
    std::string name;
    // Accessor for 3D indexing
    target_particle(double charge_in, double mass_in, double temp_in, double density_in, double v_drift_in, size_t number_cells, std::string name_in);
    double get_ave_v_sqr_cell(int cell) const;
    void print_out() const;
    void initialize_diagnostic_files(const std::string& dir_name) const;
    void write_diagnostics(const std::string& dir_name, int diag_number) const;
};
std::vector<target_particle> read_target_particle_inputs(const std::string& filename, const domain& world);

