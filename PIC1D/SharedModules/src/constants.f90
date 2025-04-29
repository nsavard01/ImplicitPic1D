module constants
    use iso_fortran_env, only: real64, int32
    use ifport, only: makedirqq
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
    integer(int32), protected :: numDiagnosticSteps , numThread
    logical, protected :: restartBool = .false.
    integer(int32), protected :: oldDiagStep, ionStepMult
    real(real64), protected :: fractionFreq, n_ave, T_e, T_i
    real(real64), protected :: del_t, simulationTime, averagingTime
    character(len=:), allocatable, protected :: directoryName, restartDirectory ! Name of save directory
    real(real64), protected :: startSimulationTime

contains

    subroutine readInitialInputs(InitFilename)
        use omp_lib
        implicit none
        character(len=*), intent(in) :: InitFilename
        integer(int32) :: io, k, u
        character(len=100) :: tempName, restartName, otherTemp, saveFolderName
        real(real64) :: plasmaFreqTemp, tempReal
        logical :: fileExists
        ! Read input parameters
        print *, ""
        print *, "Reading initial inputs:"
        print *, "------------------"
        open(10,file='../InputData/'//InitFilename, IOSTAT=io)
        read(10, *, IOSTAT = io) numThread
        read(10, *, IOSTAT = io) simulationTime, startSimulationTime
        read(10, *, IOSTAT = io) n_ave
        read(10, *, IOSTAT = io) T_e
        read(10, *, IOSTAT = io) T_i
        read(10, *, IOSTAT = io) numDiagnosticSteps
        read(10, *, IOSTAT = io) fractionFreq, del_t, ionStepMult
        read(10, *, IOSTAT = io) averagingTime
        read(10, '(A)', IOSTAT = io) saveFolderName
        read(10, *, IOSTAT = io) tempName
        read(10, *, IOSTAT = io) restartName, otherTemp
        close(10)
        oldDiagStep = 1
        do k = 1, len(saveFolderName)
            if (saveFolderName(k:k) == ' ') then
                exit
            end if
        end do
        restartBool = trim(restartName) == 'yes' .or. trim(restartName) == 'Yes' .or. trim(restartName) == 'YES'
        if (restartBool) then
            ! Read from restart input file
            restartDirectory = trim(otherTemp)
            restartDirectory = saveFolderName(1:k-1)//restartDirectory
            print *, 'Restart directory is:', restartDirectory
            print *, ""
            print *, "Reading restart inputs:"
            print *, "------------------"
            open(10,file=restartDirectory//'/InputData/'//InitFilename, IOSTAT=io)
            read(10, *, IOSTAT = io) numThread
            read(10, *, IOSTAT = io) tempReal, startSimulationTime
            read(10, *, IOSTAT = io) n_ave
            read(10, *, IOSTAT = io) T_e
            read(10, *, IOSTAT = io) T_i
            read(10, *, IOSTAT = io) u
            read(10, *, IOSTAT = io) fractionFreq, del_t, ionStepMult
            read(10, *, IOSTAT = io) tempReal
            read(10, '(A)', IOSTAT = io) saveFolderName
            read(10, *, IOSTAT = io) tempName
            close(10)
            INQUIRE( file=restartDirectory//"/"//"GlobalDiagnosticData.dat", EXIST=fileExists) 
            if (fileExists) then
                !Set new diagnostic step number
                open(10,file=restartDirectory//"/"//"GlobalDiagnosticData.dat", IOSTAT=io)
                io = 0
                read(10, *, IOSTAT = io)
                k = 0
                do while (io == 0)
                    k = k + 1
                    read(10, *, IOSTAT = io) tempReal
                end do
                close(10)
                oldDiagStep = k
                startSimulationTime = tempReal 
            else
                print *, 'Restart global Diagnostic file does not exist!'
                stop
            end if
            directoryName = restartDirectory
        else
            if (len(trim(tempName)) < 2) then
                stop "Directory name length less than 2 characters!"
            end if
            directoryName = saveFolderName(1:k-1)//trim(tempName)
        end if
        if (numThread > omp_get_num_procs()) then
            print *, "Number of threads set is larger than the maximum number of threads which is", omp_get_num_procs()
            stop
        end if
        call omp_set_num_threads(numThread)
        plasmaFreqTemp = SQRT(n_ave * (e**2) / m_e / eps_0)
        del_t = MIN(fractionFreq * 1.0d0 / plasmaFreqTemp, del_t)
        fractionFreq = del_t * plasmaFreqTemp
        simulationTime = simulationTime + startSimulationTime
        if (ionStepMult < 1 .or. ionStepMult > 1000) ionStepMult = 1
        print *, "Save data folder: ", directoryName
        print *, 'Restart Bool:', restartBool
        print *, "Number of threads is:", omp_get_max_threads()
        print *, 'Total allowable threads:', omp_get_num_procs()
        print *, "Average initial electron density:", n_ave
        print *, "Initial electron temperature:", T_e
        print *, "Initial ion temperature:", T_i
        print *, "Number of diagnostic steps is:", numDiagnosticSteps
        print *, "Fraction of 1/w_p for time step:", fractionFreq
        print *, 'del_t will be:', del_t
        print *, 'Simulation start time is:', startSimulationTime
        print *, "Simulation time is:", simulationTime
        print *, "Final averaging time is:", averagingTime
        print *, "------------------"
        print *, ""
    end subroutine readInitialInputs

    subroutine generateSaveDirectory(dirName, collNumber)
        ! Generate save directory with appropriate folders
        character(*), intent(in) :: dirName
        integer(int32), intent(in) :: collNumber
        logical :: bool
        integer(int32) :: io
        character(len=10) :: buf
        bool = makedirqq(dirName)
        if (.not. bool) then
            ! User input for rewriting same directory name
            print *, "Save directory ", dirName, " already exists. Are you sure you want to continue(yes/no)?"
            read *, buf
            if (buf(1:3) /= 'yes') then
                stop "You have decided to create a new directory for the save files I suppose"
            end if
            call execute_command_line("rm -r "//dirName//"/*", EXITSTAT = io)
            if (io /= 0) then
                stop 'Issue removing old data'
            end if
        end if
        bool = makedirqq(dirName//'/Density')
        if (.not. bool) then
            stop "Save directory not successfully created!"
        end if
        bool = makedirqq(dirName//'/PhaseSpace')
        if (.not. bool) then
            stop "Save directory not successfully created!"
        end if
        bool = makedirqq(dirName//'/Phi')
        if (.not. bool) then
            stop "Save directory not successfully created!"
        end if
        bool = makedirqq(dirName//'/Temperature')
        if (.not. bool) then
            stop "Save directory not successfully created!"
        end if
        if (collNumber > 0) then
            bool = makedirqq(dirName//'/BinaryCollisions')
            if (.not. bool) then
                stop "Save directory not successfully created!"
            end if
        end if
        ! Copy input files
        call execute_command_line("cp -Tr ../InputData "//dirName//"/InputData", EXITSTAT = io)
        if (io /= 0) then
            stop 'Issue copying input data deck'
        end if
    end subroutine generateSaveDirectory

    subroutine changeDelT(dt)
        ! Subroutine to use if del_t needs to be changed
        real(real64), intent(in) :: dt
        del_t = dt
    end subroutine changeDelT

end module constants