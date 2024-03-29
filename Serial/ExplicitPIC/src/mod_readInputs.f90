module mod_readInputs
    use iso_fortran_env, only: int32, real64
    use constants
    use mod_BasicFunctions
    use mod_domain
    use mod_particle
    use mod_potentialSolver
    use mod_collisions
    use mod_simulation
    implicit none


contains

    subroutine readRestart(InitFilename, world, solver, particleList)
        type(Domain), intent(in out) :: world
        type(potentialSolver), intent(in out) :: solver
        type(Particle), allocatable, intent(out) :: particleList(:)
        character(len=*), intent(in) :: InitFilename
        integer(int32), allocatable :: workInt(:)
        real(real64), allocatable :: workReal(:)
        character(len=:), allocatable :: filename
        character(len=5) :: char_i
        character(len=10) :: name
        character(len=100) :: buf
        real(real64) :: mass, charge, w_p, temp
        integer(int32) :: io, i, finalIdx, N_p, j
        logical :: bool
        print *, "So you've decided to create a simulation based on the ending of another simulation."
        print *, "Please type in the name of the file to extract data from."
        print *, "If you're having second thoughts, put in a file that doesn't exist, and I'll quit on you."
        read *, buf
        filename = trim(buf)
        INQUIRE(DIRECTORY = '../'//filename, EXIST = bool)
        if (bool) then
            open(10,file='../'//filename//'/InitialConditions.dat')
            read(10, *)
            read(10, *) NumberXNodes, simulationTime, del_t, FractionFreq, world%delX, n_ave, T_e, T_i, numDiagnosticSteps, numberChargedParticles
            close(10)
            allocate(workInt(NumberXNodes), workReal(NumberXNodes))
            open(10, file = '../'//filename//'/domainBoundaryConditions.dat', form = 'unformatted')
            read(10) workInt
            close(10)
            world = Domain(workInt(1), workInt(NumberXNodes))
            open(10, file = '../'//filename//'/domainGrid.dat', form = 'unformatted')
            read(10) workReal
            close(10)
            world%grid = workReal
            world%delX = world%grid(2) - world%grid(1)

            write(char_i, '(I3)') numDiagnosticSteps
            open(10, file = '../'//filename//'/Phi/phi_'//trim(adjustl(char_i))//'.dat', form = 'unformatted')
            read(10) workReal
            close(10)
            solver = potentialSolver(world, workReal(1), workReal(NumberXNodes))
            call solver%construct_diagMatrix(world)
            solver%phi = workReal
            allocate(particleList(numberChargedParticles))
            print *, 'size particleList is:', numberChargedParticles
            open(10,file='../'//filename//'/ParticleProperties.dat')
            read(10, *)
            do i = 1, numberChargedParticles
                read(10, *) name, mass, charge, w_p, finalIdx
                open(13, file='../'//filename//'/ParticleDiagnostic_'//trim(name)//'.dat')
                do j = 1, numDiagnosticSteps
                    read(13, *)
                end do
                read(13, *) temp, temp, temp, temp, temp, N_p
                particleList(i) = Particle(mass, charge, w_p, N_p, finalIdx, trim(name))
                open(15, file='../'//filename//'/PhaseSpace/phaseSpace_'//trim(name)//'_'//trim(adjustl(char_i))//'.dat', form = 'unformatted')
                read(15) particleList(i)%phaseSpace(:, 1:particleList(i)%N_p)
                close(15)
            end do
            close(13)
            close(10)
            open(10,file='../InputData/'//InitFilename, IOSTAT=io)
            read(10, *, IOSTAT = io) simulationTime
            read(10, *, IOSTAT = io) 
            read(10, *, IOSTAT = io) 
            read(10, *, IOSTAT = io) 
            read(10, *, IOSTAT = io) numDiagnosticSteps
            read(10, *, IOSTAT = io) fractionFreq
            read(10, *, IOSTAT = io) averagingTime
            read(10, *, IOSTAT = io)
            read(10, *, IOSTAT = io) 
            read(10, *, IOSTAT = io) 
            read(10, *, IOSTAT = io) buf
            close(10)
            directoryName = trim(buf)
            if (len(directoryName) < 2) then
                stop "Directory name length less than 2 characters, be more descriptive!"
            end if
            print *, ""
            print *, "--------------------"
            print *, "Save data folder: ", directoryName
            print *, "Average initial electron density:", n_ave
            print *, "Initial electron temperature:", T_e
            print *, "Initial ion temperature:", T_i
            print *, "Number of diagnostic steps is:", numDiagnosticSteps
            print *, "Fraction of 1/w_p for time step:", fractionFreq
            print *, "Final averaging time is:", averagingTime
            print *, "------------------"
            print *, ""
            print *, "------------------"
            print *, "Number of nodes:", NumberXNodes
            print *, "Grid length:", world%grid(NumberXNodes) - world%grid(1)
            print *, "Left boundary type:", world%boundaryConditions(1)
            print *, "Right boundary type:", world%boundaryConditions(NumberXNodes)
            print *, "------------------"
            print *, ""
            do j=1, numberChargedParticles
                print *, 'Particle name ', particleList(j) % name
                print *, 'Amount of macroparticles is:', particleList(j) % N_p
                print *, "Particle mass is:", particleList(j)%mass
                print *, "Particle charge is:", particleList(j)%q
                print *, "Particle mean KE is:", particleList(j)%getKEAve()
                print *, ""
            end do
        else
            stop "File you're trying to restart from doesn't exist."
        end if
    end subroutine readRestart

    subroutine readInitialInputs(InitFilename, simulationTime, n_ave, T_e, T_i, numDiagnosticSteps, fractionFreq, averagingTime)
        real(real64), intent(in out) :: fractionFreq, n_ave, simulationTime, averagingTime, T_e, T_i
        integer(int32), intent(in out) :: numDiagnosticSteps
        character(len=*), intent(in) :: InitFilename
        integer(int32) :: io
        character(len=100) :: tempName
        print *, ""
        print *, "Reading initial inputs:"
        print *, "------------------"
        open(10,file='../InputData/'//InitFilename, IOSTAT=io)
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
        print *, "Average initial electron density:", n_ave
        print *, "Initial electron temperature:", T_e
        print *, "Initial ion temperature:", T_i
        print *, "Number of diagnostic steps is:", numDiagnosticSteps
        print *, "Fraction of 1/w_p for time step:", fractionFreq
        print *, "Final averaging time is:", averagingTime
        print *, "------------------"
        print *, ""
    end subroutine readInitialInputs

    subroutine readGeometry(world, solver, GeomFilename)
        type(Domain), intent(in out) :: world
        type(potentialSolver), intent(in out) :: solver
        character(len=*), intent(in) :: GeomFilename
        integer(int32) :: io, leftBoundary, rightBoundary
        real(real64) :: leftVoltage, rightVoltage, debyeLength, L_domain

        print *, ""
        print *, "Reading domain inputs:"
        print *, "------------------"
        open(10,file='../InputData/'//GeomFilename)
        read(10, *, IOSTAT = io) NumberXNodes
        read(10, *, IOSTAT = io) L_domain
        read(10, *, IOSTAT = io) leftBoundary, rightBoundary
        read(10, *, IOSTAT = io) leftVoltage, rightVoltage
        close(10)
        debyeLength = getDebyeLength(T_e, n_ave)
        if (L_domain / (NumberXNodes-1) > debyeLength) then
            print *, "Insufficient amount of nodes to resolve initial debyeLength"
            print *, "Changing amount of nodes to have 0.75 * debyeLength"
            NumberXNodes = NINT(L_domain/debyeLength/0.75d0) + 1
        end if
        print *, "Number of nodes:", NumberXNodes
        print *, "Grid length:", L_domain
        print *, "Left boundary type:", leftBoundary
        print *, "Right boundary type:", rightBoundary
        print *, "------------------"
        print *, ""
        if ((leftBoundary == 3) .or. (rightBoundary == 3)) then
            leftBoundary = 3
            rightBoundary = 3
            leftVoltage = rightVoltage
        end if
        world = Domain(leftBoundary, rightBoundary)
        call world % constructGrid(L_domain)
        solver = potentialSolver(world, leftVoltage, rightVoltage)
        call solver%construct_diagMatrix(world)

    end subroutine readGeometry

    function readParticleInputs(filename, numberChargedParticles, irand, T_e, T_i) result(particleList)
        type(Particle), allocatable :: particleList(:)
        character(len=*), intent(in) :: filename
        integer(int32), intent(in out) :: numberChargedParticles, irand
        real(real64), intent(in) :: T_e, T_i
        integer(int32) :: j, numSpecies = 0, numParticles(100), particleIdxFactor(100)
        character(len=15) :: name
        character(len=8) :: particleNames(100)
        real(real64) :: mass(100), charge(100), Ti(100)

        print *, "Reading particle inputs:"
        open(10,file='../InputData/'//filename, action = 'read')
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
            particleList(j) = Particle(mass(j), e * charge(j), 1.0d0, numParticles(j) * (NumberXNodes-1), numParticles(j) * particleIdxFactor(j) * (NumberXNodes - 1), trim(particleNames(j)))
            call particleList(j) % generate3DMaxwellian(Ti(j), irand)
            print *, 'Initializing ', particleList(j) % name
            print *, 'Amount of macroparticles is:', particleList(j) % N_p
            print *, "Particle mass is:", particleList(j)%mass
            print *, "Particle charge is:", particleList(j)%q
            print *, "Particle mean KE is:", particleList(j)%getKEAve(), ", should be", Ti(j) * 1.5
        end do
        
        print *, "---------------"
        print *, ""



    end function readParticleInputs

end module mod_readInputs