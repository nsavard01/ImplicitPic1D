#ifndef CONSTANTS_H
#define CONSTANTS_H

namespace Constants {
    constexpr double speed_of_light = 299792458.0; // m/s
    constexpr double electron_mass = 9.1093837015e-31; // kg
    constexpr double elementary_charge = 1.602176634e-19; // C
    constexpr double epsilon_0 = 8.8541878128e-12 ; // Vacuum permittivity (F/m)
    constexpr double mu_0 = 1.25663706212e-6; // Vacuum permeability (N/A^2)
    constexpr double mass_amu = 1.66053906660e-27; // kg
    constexpr double k_boltz = 1.380649e-23; // m^2 kg s^-2 K^-1
    extern const double cos_third_rot, sin_third_rot;
    extern int mpi_rank, mpi_size;
}

#endif // CONSTANTS_H