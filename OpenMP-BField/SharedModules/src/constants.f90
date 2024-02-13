module constants
    use iso_fortran_env, only: real64, int32
    implicit none

    ! physical constants and global variables used among all modules
    real(real64), parameter :: c = 299792458.0d0 ! m/s
    real(real64), parameter :: eps_0 = 8.8541878128d-12 !F/m
    real(real64), parameter :: m_p = 1.67262192369d-27 ! kg
    real(real64), parameter :: m_amu = 1.66053906660d-27 ! kg
    real(real64), parameter :: m_e = 9.1093837015d-31 ! kg
    real(real64), parameter :: k_B = 1.380649d-23 ! m^2 kg s^-2 K^-1
    real(real64), parameter :: e = 1.602176634d-19 ! C
    real(real64), parameter :: mu_0 = 1.25663706212d-6 ! m kg s^-2 A^-2
    real(real64), parameter :: pi = 4.0d0*atan(1.0d0) ! pi from atan
    ! Essential parameters set that is important for entire simulation
    integer(int32) :: maxIter, numDiagnosticSteps, schemeNum, numThread
    real(real64) :: eps_r, fractionFreq, n_ave, T_e, T_i
    real(real64) :: del_t, simulationTime, averagingTime
    character(len=:), allocatable :: directoryName ! Name of save directory

end module constants