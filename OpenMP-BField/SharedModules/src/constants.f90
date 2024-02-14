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
    ! Essential parameters set that is important for entire simulation state
    integer(int32), protected :: numDiagnosticSteps, numThread
    real(real64), protected :: fractionFreq, n_ave, T_e, T_i
    real(real64), protected :: del_t, simulationTime, averagingTime
    character(len=:), allocatable, protected :: directoryName ! Name of save directory

contains

    subroutine readInitialInputs(InitFilename, saveFolderName)
        use omp_lib
        implicit none
        character(len=*), intent(in) :: InitFilename, saveFolderName
        integer(int32) :: io
        character(len=100) :: tempName
        real(real64) :: plasmaFreqTemp
        print *, ""
        print *, "Reading initial inputs:"
        print *, "------------------"
        open(10,file='../InputData/'//InitFilename, IOSTAT=io)
        read(10, *, IOSTAT = io) numThread
        read(10, *, IOSTAT = io) simulationTime
        read(10, *, IOSTAT = io) n_ave
        read(10, *, IOSTAT = io) T_e
        read(10, *, IOSTAT = io) T_i
        read(10, *, IOSTAT = io) numDiagnosticSteps
        read(10, *, IOSTAT = io) fractionFreq, del_t
        read(10, *, IOSTAT = io) averagingTime
        read(10, *, IOSTAT = io) tempName
        close(10)
        directoryName = trim(tempName)
        if (len(directoryName) < 2) then
            stop "Directory name length less than 2 characters!"
        end if
        directoryName = saveFolderName//directoryName
        if (numThread > omp_get_num_procs()) then
            print *, "Number of threads set is larger than the maximum number of threads which is", omp_get_num_procs()
            stop
        end if
        call omp_set_num_threads(numThread)
        plasmaFreqTemp = SQRT(n_ave * (e**2) / m_e / eps_0)
        del_t = MIN(fractionFreq * 1.0d0 / plasmaFreqTemp, del_t)
        print *, "Save data folder: ", directoryName
        print *, "Number of threads is:", numThread
        print *, "Average initial electron density:", n_ave
        print *, "Initial electron temperature:", T_e
        print *, "Initial ion temperature:", T_i
        print *, "Number of diagnostic steps is:", numDiagnosticSteps
        print *, "Fraction of 1/w_p for time step:", fractionFreq
        print *, 'del_t will be:', del_t
        print *, "Simulation time is:", simulationTime
        print *, "Final averaging time is:", averagingTime
        print *, "------------------"
        print *, ""
    end subroutine readInitialInputs

end module constants