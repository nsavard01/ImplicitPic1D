module mod_readInputs
    use iso_fortran_env, only: int32, real64
    use constants
    use mod_BasicFunctions
    use mod_domain
    use mod_particle
    use mod_potentialSolver
    use mod_Scheme
    use mod_particleInjection
    use omp_lib
    implicit none


contains

    ! subroutine readRestart(InitFilename, world, solver, particleList)
    !     type(Domain), intent(in out) :: world
    !     type(potentialSolver), intent(in out) :: solver
    !     type(Particle), allocatable, intent(out) :: particleList(:)
    !     character(len=*), intent(in) :: InitFilename
    !     integer(int32), allocatable :: workInt(:)
    !     real(real64), allocatable :: workReal(:)
    !     character(len=:), allocatable :: filename
    !     character(len=5) :: char_i
    !     character(len=10) :: name
    !     character(len=100) :: buf
    !     real(real64) :: mass, charge, w_p, temp
    !     integer(int32) :: io, i, finalIdx, N_p, j
    !     logical :: bool
    !     print *, "So you've decided to create a simulation based on the ending of another simulation."
    !     print *, "Please type in the name of the file to extract data from."
    !     print *, "If you're having second thoughts, put in a file that doesn't exist, and I'll quit on you."
    !     read *, buf
    !     filename = trim(buf)
    !     INQUIRE(DIRECTORY = '../'//filename, EXIST = bool)
    !     if (bool) then
    !         open(10,file='../'//filename//'/InitialConditions.dat')
    !         read(10, *)
    !         read(10, *) NumberXNodes, simulationTime, del_t, FractionFreq, world%delX, n_ave, T_e, T_i, numDiagnosticSteps, numberChargedParticles
    !         close(10)
    !         allocate(workInt(NumberXNodes), workReal(NumberXNodes))
    !         open(10, file = '../'//filename//'/domainBoundaryConditions.dat', form = 'unformatted')
    !         read(10) workInt
    !         close(10)
    !         world = Domain(workInt(1), workInt(NumberXNodes))
    !         open(10, file = '../'//filename//'/domainGrid.dat', form = 'unformatted')
    !         read(10) workReal
    !         close(10)
    !         world%grid = workReal
    !         world%delX = world%grid(2) - world%grid(1)

    !         write(char_i, '(I3)') numDiagnosticSteps
    !         open(10, file = '../'//filename//'/Phi/phi_'//trim(adjustl(char_i))//'.dat', form = 'unformatted')
    !         read(10) workReal
    !         close(10)
    !         solver = potentialSolver(world, workReal(1), workReal(NumberXNodes))
    !         call solver%construct_diagMatrix(world)
    !         solver%phi = workReal
    !         allocate(particleList(numberChargedParticles))
    !         print *, 'size particleList is:', numberChargedParticles
    !         open(10,file='../'//filename//'/ParticleProperties.dat')
    !         read(10, *)
    !         do i = 1, numberChargedParticles
    !             read(10, *) name, mass, charge, w_p, finalIdx
    !             open(13, file='../'//filename//'/ParticleDiagnostic_'//trim(name)//'.dat')
    !             do j = 1, numDiagnosticSteps
    !                 read(13, *)
    !             end do
    !             read(13, *) temp, temp, temp, temp, temp, N_p
    !             particleList(i) = Particle(mass, charge, w_p, N_p, finalIdx, trim(name))
    !             open(15, file='../'//filename//'/PhaseSpace/phaseSpace_'//trim(name)//'_'//trim(adjustl(char_i))//'.dat', form = 'unformatted')
    !             read(15) particleList(i)%phaseSpace(:, 1:particleList(i)%N_p)
    !             close(15)
    !         end do
    !         close(13)
    !         close(10)
    !         open(10,file='../InputData/'//InitFilename, IOSTAT=io)
    !         read(10, *, IOSTAT = io) simulationTime
    !         read(10, *, IOSTAT = io) 
    !         read(10, *, IOSTAT = io) 
    !         read(10, *, IOSTAT = io) 
    !         read(10, *, IOSTAT = io) numDiagnosticSteps
    !         read(10, *, IOSTAT = io) fractionFreq
    !         read(10, *, IOSTAT = io) averagingTime
    !         read(10, *, IOSTAT = io)
    !         read(10, *, IOSTAT = io) 
    !         read(10, *, IOSTAT = io) 
    !         read(10, *, IOSTAT = io) buf
    !         close(10)
    !         directoryName = trim(buf)
    !         if (len(directoryName) < 2) then
    !             stop "Directory name length less than 2 characters, be more descriptive!"
    !         end if
    !         print *, ""
    !         print *, "--------------------"
    !         print *, "Save data folder: ", directoryName
    !         print *, "Average initial electron density:", n_ave
    !         print *, "Initial electron temperature:", T_e
    !         print *, "Initial ion temperature:", T_i
    !         print *, "Number of diagnostic steps is:", numDiagnosticSteps
    !         print *, "Fraction of 1/w_p for time step:", fractionFreq
    !         print *, "Final averaging time is:", averagingTime
    !         print *, "------------------"
    !         print *, ""
    !         print *, "------------------"
    !         print *, "Number of nodes:", NumberXNodes
    !         print *, "Grid length:", world%grid(NumberXNodes) - world%grid(1)
    !         print *, "Left boundary type:", world%boundaryConditions(1)
    !         print *, "Right boundary type:", world%boundaryConditions(NumberXNodes)
    !         print *, "------------------"
    !         print *, ""
    !         do j=1, numberChargedParticles
    !             print *, 'Particle name ', particleList(j) % name
    !             print *, 'Amount of macroparticles is:', particleList(j) % N_p
    !             print *, "Particle mass is:", particleList(j)%mass
    !             print *, "Particle charge is:", particleList(j)%q
    !             print *, "Particle mean KE is:", particleList(j)%getKEAve()
    !             print *, ""
    !         end do
    !     else
    !         stop "File you're trying to restart from doesn't exist."
    !     end if
    ! end subroutine readRestart

    subroutine readInjectionInputs(InjFilename, addLostPartBool, refluxPartBool, injectionBool, injectionFlux, w_p)
        logical, intent(in out) :: addLostPartBool, refluxPartBool, injectionBool
        real(real64), intent(in out) :: injectionFlux
        real(real64), intent(in) :: w_p
        character(len=*), intent(in) :: InjFilename
        integer(int32) :: tempInt, io
        real(real64) :: numFluxPart
        print *, ""
        print *, "Reading initial inputs for particle injection:"
        print *, "------------------"
        open(10,file='../../SharedModules/InputData/'//InjFilename, IOSTAT=io)
        read(10, *, IOSTAT = io) tempInt
        addLostPartBool = (tempInt == 1)
        read(10, *, IOSTAT = io) tempInt
        refluxPartBool = (tempInt == 1)
        read(10, *, IOSTAT = io) tempInt, injectionFlux
        injectionBool = (tempInt == 1)
        close(10)
        print *, "Particle lost is reinjected:", addLostPartBool
        print *, "Particle refluxing activated on neumann boundary:", refluxPartBool
        print *, "Particle injection on neumann boundary", injectionBool
        if (injectionBool) then
            print *, 'Particle injection flux:', injectionFlux
            numFluxPart = injectionFlux * del_t / w_p/real(numThread) ! particles injected per thread
            numFluxParticlesLow = floor(numFluxPart)
            numFluxParticlesHigh = numFluxParticlesLow + 1
            print *, 'Low end of flux particles:', numFluxParticlesLow
            print *, 'High end of flux particles:', numFluxParticlesHigh
            injectionR = numFluxPart - real(numFluxParticlesLow)
            print *, 'Number for selection of flux particles is:', injectionR
        end if
        print *, "------------------"
        print *, ""
    end subroutine readInjectionInputs

    subroutine readInitialInputs(InitFilename, simulationTime, n_ave, T_e, T_i, numDiagnosticSteps, fractionFreq, averagingTime, numThread, irand)
        real(real64), intent(in out) :: fractionFreq, n_ave, simulationTime, averagingTime, T_e, T_i
        integer(int32), intent(in out) :: numDiagnosticSteps, numThread
        integer(int32), allocatable, intent(out) :: irand(:)
        character(len=*), intent(in) :: InitFilename
        integer(int32) :: io, i
        character(len=100) :: tempName
        print *, ""
        print *, "Reading initial inputs:"
        print *, "------------------"
        open(10,file='../../SharedModules/InputData/'//InitFilename, IOSTAT=io)
        read(10, *, IOSTAT = io) numThread
        read(10, *, IOSTAT = io) simulationTime
        read(10, *, IOSTAT = io) n_ave
        read(10, *, IOSTAT = io) T_e
        read(10, *, IOSTAT = io) T_i
        read(10, *, IOSTAT = io) numDiagnosticSteps
        read(10, *, IOSTAT = io) fractionFreq
        read(10, *, IOSTAT = io) averagingTime
        read(10, *, IOSTAT = io) tempName
        close(10)
        directoryName = trim(tempName)
        if (len(directoryName) < 2) then
            stop "Directory name length less than 2 characters!"
        end if
        directoryName = '../../../../ImplicitData-BField/'//directoryName
        if (numThread > omp_get_num_procs()) then
            print *, "Number of threads set is larger than the maximum number of threads which is", omp_get_num_procs()
            stop
        end if
        call omp_set_num_threads(numThread)
        allocate(irand(numThread))
        do i = 1, numThread
            irand(i) = 123456*i*11
        end do
        del_t = fractionFreq * 1.0d0 / getPlasmaFreq(n_ave)
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

    subroutine readGeometry(world, solver, GeomFilename)
        type(Domain), intent(in out) :: world
        type(potentialSolver), intent(in out) :: solver
        character(len=*), intent(in) :: GeomFilename
        integer(int32) :: io, leftBoundary, rightBoundary, gridType
        real(real64) :: leftVoltage, rightVoltage, debyeLength, L_domain, BFieldMag, angle

        print *, ""
        print *, "Reading domain inputs:"
        print *, "------------------"
        open(10,file='../../SharedModules/InputData/'//GeomFilename)
        read(10, *, IOSTAT = io) NumberXNodes
        read(10, *, IOSTAT = io) L_domain
        read(10, *, IOSTAT = io) gridType
        read(10, *, IOSTAT = io) leftBoundary, rightBoundary
        read(10, *, IOSTAT = io) leftVoltage, rightVoltage
        read(10, *, IOSTAT = io) BFieldMag, angle
        close(10)
        debyeLength = getDebyeLength(T_e, n_ave)
        if ((leftBoundary == 3) .or. (rightBoundary == 3)) then
            leftBoundary = 3
            rightBoundary = 3
            leftVoltage = rightVoltage
        end if
        world = Domain(leftBoundary, rightBoundary)
        call world % constructGrid(debyeLength, L_domain, gridType)
        print *, "Number of nodes:", NumberXNodes
        print *, "Grid length:", world%grid(NumberXNodes) - world%grid(1)
        print *, "Left boundary type:", leftBoundary
        print *, "Right boundary type:", rightBoundary
        print *, 'Grid type is:', gridType
        solver = potentialSolver(world, leftVoltage, rightVoltage, BFieldMag, angle)
        print *, "BField vector:", solver%BField
        print *, "------------------"
        print *, ""
        call solver%construct_diagMatrix(world)

    end subroutine readGeometry

    function readParticleInputs(filename, numberChargedParticles, irand, T_e, T_i, numThread, world) result(particleList)
        type(Particle), allocatable :: particleList(:)
        type(Domain) :: world
        character(len=*), intent(in) :: filename
        integer(int32), intent(in) :: numThread
        integer(int32), intent(in out) :: numberChargedParticles, irand(numThread)
        real(real64), intent(in) :: T_e, T_i
        integer(int32) :: j, numSpecies = 0, numParticles(100), particleIdxFactor(100)
        character(len=15) :: name
        character(len=8) :: particleNames(100)
        real(real64) :: mass(100), charge(100), Ti(100)

        print *, "Reading particle inputs:"
        open(10,file='../../SharedModules/InputData/'//filename, action = 'read')
        do j=1, 10000
            read(10,*,END=101,ERR=100) name

            if( name(1:9).eq.'ELECTRONS') then
                read(10,*,END=101,ERR=100) name
                read(10,*,END=101,ERR=100) name
                read(10,'(A4)',END=101,ERR=100, ADVANCE = 'NO') name(1:2)
                numSpecies = numSpecies + 1
                read(10,*,END=101,ERR=100) numParticles(numSpecies), particleIdxFactor(numSpecies)
                Ti(numSpecies) = T_e
                mass(numSpecies) = m_e
                charge(numSpecies) = -1.0
                particleNames(numSpecies) = 'e'
                read(10,*,END=101,ERR=100) name
                read(10,*,END=101,ERR=100) name
            endif


            if(name(1:4).eq.'IONS' .or. name(1:4).eq.'Ions' .or. name(1:4).eq.'ions' ) then
                do while(name(1:4).ne.'----')
                    read(10,*,END=101,ERR=100) name
                end do
200             read(10,'(A6)',END=101,ERR=100, ADVANCE = 'NO') name
                if (name(1:4).eq.'----') then
                    close(10)
                else
                    numSpecies = numSpecies + 1
                    read(10,*,END=101,ERR=100) mass(numSpecies),charge(numSpecies), numParticles(numSpecies), particleIdxFactor(numSpecies)
                    Ti(numSpecies) = T_i
                    mass(numSpecies) = mass(numSpecies) * m_p
                    particleNames(numSpecies) = trim(name)
                    goto 200
                end if
            endif
            ! Take care of extra text I guess        

            if (name(1:7) == 'ENDFILE') then
                close(10)
            end if

        end do
100     continue
101     continue
        numberChargedParticles = numSpecies
        allocate(particleList(numberChargedParticles))

        do j=1, numberChargedParticles
            particleList(j) = Particle(mass(j), e * charge(j), 1.0d0, numParticles(j), numParticles(j) * particleIdxFactor(j), trim(particleNames(j)), numThread)
            particleList(j) % w_p = n_ave * (world%grid(NumberXNodes) - world%grid(1)) / SUM(particleList(j) % N_p)
            call particleList(j) % generate3DMaxwellian(Ti(j), irand)
            if (j == 1) then
                call initialize_randUniform(particleList(j), world, irand)
            else
                call initialize_QuasiNeutral(particleList(j), particleList(1))
            end if
            print *, 'Initializing ', particleList(j) % name
            print *, 'Amount of macroparticles is:', SUM(particleList(j) % N_p)
            print *, "Particle mass is:", particleList(j)%mass
            print *, "Particle charge is:", particleList(j)%q
            print *, "Particle weight is:", particleList(j)%w_p
            print *, "Particle mean KE is:", particleList(j)%getKEAve(), ", should be", Ti(j) * 1.5
        end do
        
        print *, "---------------"
        print *, ""



    end function readParticleInputs

end module mod_readInputs