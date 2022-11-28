module constants
    use iso_fortran_env, only: real64
    implicit none

    
    real(real64), parameter :: c = 299792458.0 ! m/s
    real(real64), parameter :: eps_0 = 8.8541878128e-12 !F/m
    real(real64), parameter :: m_p = 1.67262192369e-27 ! kg
    real(real64), parameter :: m_e = 9.1093837015e-31 ! kg
    real(real64), parameter :: k_B = 1.380649e-23 ! m^2 kg s^-2 K^-1
    real(real64), parameter :: e = 1.602176634e-19 ! C
    real(real64), parameter :: mu_0 = 1.25663706212e-6 ! m kg s^-2 A^-2
    real(real64), parameter :: pi = 3.141592653589793


end module constants