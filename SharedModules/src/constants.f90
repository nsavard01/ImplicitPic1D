module constants
    use iso_fortran_env, only: real64, int32
    implicit none

    ! physical constants and global variables used among all modules
    real(real64), parameter :: c = 299792458.0d0 ! m/s
    real(real64), parameter :: eps_0 = 8.8541878128d-12 !F/m
    real(real64), parameter :: m_p = 1.67262192369d-27 ! kg
    real(real64), parameter :: m_e = 9.1093837015d-31 ! kg
    real(real64), parameter :: k_B = 1.380649d-23 ! m^2 kg s^-2 K^-1
    real(real64), parameter :: e = 1.602176634d-19 ! C
    real(real64), parameter :: mu_0 = 1.25663706212d-6 ! m kg s^-2 A^-2
    real(real64), parameter :: pi = 4.0d0*atan(1.0d0) ! pi from atan
    integer(int32) :: NumberXNodes, numberChargedParticles, maxIter, numDiagnosticSteps
    integer(int32) :: irand = 987888 !Global variable that will need to be changed in-out
    real(real64) :: eps_r, fractionFreq, n_ave, T_e, T_i, simulationTime, averagingTime
    character(len=:), allocatable :: directoryName

contains

    subroutine readInitialConditions(InitFileName)
        character(len=*), intent(in) :: InitFilename
        integer(int32) :: io
        character(len=100) :: tempName
        print *, ""
        print *, "--------------------"
        print *, "Reading initial inputs:"
        open(10,file='../../SharedModules/InputData/'//InitFilename, IOSTAT=io)
        read(10, *, IOSTAT = io) simulationTime
        read(10, *, IOSTAT = io) n_ave
        read(10, *, IOSTAT = io) T_e
        read(10, *, IOSTAT = io) T_i
        read(10, *, IOSTAT = io) numDiagnosticSteps
        read(10, *, IOSTAT = io) fractionFreq
        read(10, *, IOSTAT = io) averagingTime
        read(10, *, IOSTAT = io) 
        read(10, *, IOSTAT = io) 
        read(10, *, IOSTAT = io) 
        read(10, *, IOSTAT = io) tempName
        close(10)
        directoryName = trim(tempName)
        if (len(directoryName) < 2) then
            stop "Directory name length less than 2 characters!"
        end if
        print *, "Save data folder: ", directoryName
        print *, "Average initial particle density:", n_ave
        print *, "Initial T_e:", T_e
        print *, "Initial T_i:", T_i
        print *, "Number of diagnostic steps is:", numDiagnosticSteps
        print *, "Fraction of 1/w_p for time step:", fractionFreq
        print *, "Final averaging time is:", averagingTime
        print *, "--------------------"
        print *, ""
    end subroutine readInitialConditions

end module constants