module mod_simulation
    ! module which actually uses others to run a simulation (single or multiple time steps)
    ! will contain diagnostics as well
    use iso_fortran_env, only: int32, real64, int64
    use constants
    use mod_BasicFunctions
    use mod_particle
    use mod_domain
    use mod_potentialSolver
    use mod_collisions
    use ifport, only: makedirqq
    implicit none

    integer(int32) :: numTimeSteps
    real(real64) :: del_t, simulationTime, averagingTime
    character(len=:), allocatable :: directoryName

contains

    ! ------------------------- Reading Input data --------------------------------

!     subroutine readInputs(NumberXNodes, numDiagnosticSteps, averagingTime, fractionFreq, n_ave, world, solver, simulationTime, Power, heatSkipSteps, nu_h, T_e, T_i, GeomFilename, InitFilename)
!         ! Set initial conditions and global constants based on read input from txt file, create world and solver from these inputs
!         integer(int32), intent(in out) :: NumberXNodes, numDiagnosticSteps, heatSkipSteps
!         real(real64), intent(in out) :: fractionFreq, n_ave, simulationTime, Power, nu_h, averagingTime
!         real(real64), intent(in out) :: T_e, T_i
!         character(len=*), intent(in) :: GeomFilename, InitFilename
!         integer(int32) :: io, leftBoundary, rightBoundary
!         real(real64) :: leftVoltage, rightVoltage, L_domain, debyeLength
!         character(len=100) :: tempName
!         type(Domain) :: world
!         type(potentialSolver) :: solver

!         print *, "Reading initial inputs:"
!         open(10,file='../InputData/'//InitFilename, IOSTAT=io)
!         read(10, *, IOSTAT = io) simulationTime
!         read(10, *, IOSTAT = io) n_ave
!         read(10, *, IOSTAT = io) T_e
!         read(10, *, IOSTAT = io) T_i
!         read(10, *, IOSTAT = io) numDiagnosticSteps
!         read(10, *, IOSTAT = io) fractionFreq
!         read(10, *, IOSTAT = io) averagingTime
!         read(10, *, IOSTAT = io) Power
!         read(10, *, IOSTAT = io) heatSkipSteps
!         read(10, *, IOSTAT = io) nu_h
!         read(10, *, IOSTAT = io) tempName
!         close(10)
!         directoryName = trim(tempName)
!         if (len(directoryName) < 2) then
!             stop "Directory name length less than 2 characters!"
!         end if
!         print *, "Save data folder: ", directoryName
!         print *, "Average initial electron density:", n_ave
!         print *, "Initial electron temperature:", T_e
!         print *, "Initial ion temperature:", T_i
!         print *, "Number of diagnostic steps is:", numDiagnosticSteps
!         print *, "Fraction of 1/w_p for time step:", fractionFreq
!         print *, "Final averaging time is:", averagingTime
!         print *, "Power input (W/m^2):", Power
!         print *, "Steps to skip for heating:", heatSkipSteps
!         print *, "Heating frequency (Hz):", nu_h
!         print *, "------------------"
!         print *, ""
!         print *, "Reading domain inputs:"
!         open(10,file='../InputData/'//GeomFilename)
!         read(10, *, IOSTAT = io) NumberXNodes
!         read(10, *, IOSTAT = io) L_domain
!         read(10, *, IOSTAT = io) leftBoundary, rightBoundary
!         read(10, *, IOSTAT = io) leftVoltage, rightVoltage
!         close(10)
!         debyeLength = getDebyeLength(T_e, n_ave)
!         if (L_domain / (NumberXNodes-1) > debyeLength) then
!             print *, "Insufficient amount of nodes to resolve initial debyeLength"
!             print *, "Changing amount of nodes to have 0.75 * debyeLength"
!             NumberXNodes = NINT(L_domain/debyeLength/0.75d0) + 1
!         end if
!         print *, "Number of nodes:", NumberXNodes
!         print *, "Grid length:", L_domain
!         print *, "Left boundary type:", leftBoundary
!         print *, "Right boundary type:", rightBoundary
!         print *, "------------------"
!         print *, ""
!         ! if one boundary is periodic, other must also be
!         if ((leftBoundary == 3) .or. (rightBoundary == 3)) then
!             leftBoundary = 3
!             rightBoundary = 3
!             leftVoltage = rightVoltage
!         end if
!         world = Domain(leftBoundary, rightBoundary)
!         call world % constructGrid(L_domain)
!         solver = potentialSolver(world, leftVoltage, rightVoltage)
!         call solver%construct_diagMatrix(world)
        
!     end subroutine readInputs

!     function readParticleInputs(filename, numberChargedParticles, irand, T_e, T_i) result(particleList)
!         type(Particle), allocatable :: particleList(:)
!         character(len=*), intent(in) :: filename
!         integer(int32), intent(in out) :: numberChargedParticles, irand
!         real(real64), intent(in) :: T_e, T_i
!         integer(int32) :: j, numSpecies = 0, numParticles(100), particleIdxFactor(100)
!         character(len=15) :: name
!         character(len=8) :: particleNames(100)
!         real(real64) :: mass(100), charge(100), Ti(100)

!         print *, "Reading particle inputs:"
!         open(10,file='../InputData/'//filename, action = 'read')
!         do j=1, 10000
!             read(10,*,END=101,ERR=100) name

!             if( name(1:9).eq.'ELECTRONS') then
!                 read(10,*,END=101,ERR=100) name
!                 read(10,*,END=101,ERR=100) name
!                 read(10,'(A4)',END=101,ERR=100, ADVANCE = 'NO') name(1:2)
!                 numSpecies = numSpecies + 1
!                 read(10,*,END=101,ERR=100) numParticles(numSpecies), particleIdxFactor(numSpecies)
!                 Ti(numSpecies) = T_e
!                 mass(numSpecies) = m_e
!                 charge(numSpecies) = -1.0
!                 particleNames(numSpecies) = 'e'
!                 read(10,*,END=101,ERR=100) name
!                 read(10,*,END=101,ERR=100) name
!             endif


!             if(name(1:4).eq.'IONS' .or. name(1:4).eq.'Ions' .or. name(1:4).eq.'ions' ) then
!                 do while(name(1:4).ne.'----')
!                     read(10,*,END=101,ERR=100) name
!                 end do
! 200             read(10,'(A6)',END=101,ERR=100, ADVANCE = 'NO') name
!                 if (name(1:4).eq.'----') then
!                     close(10)
!                 else
!                     numSpecies = numSpecies + 1
!                     read(10,*,END=101,ERR=100) mass(numSpecies),charge(numSpecies), numParticles(numSpecies), particleIdxFactor(numSpecies)
!                     Ti(numSpecies) = T_i
!                     mass(numSpecies) = mass(numSpecies) * m_p
!                     particleNames(numSpecies) = trim(name)
!                     goto 200
!                 end if
!             endif
!             ! Take care of extra text I guess        

!             if (name(1:7) == 'ENDFILE') then
!                 close(10)
!             end if

!         end do
! 100     continue
! 101     continue
!         numberChargedParticles = numSpecies
!         allocate(particleList(numberChargedParticles))
!         do j=1, numberChargedParticles
!             particleList(j) = Particle(mass(j), e * charge(j), 1.0d0, numParticles(j) * (NumberXNodes-1), numParticles(j) * particleIdxFactor(j) * (NumberXNodes - 1), trim(particleNames(j)))
!             call particleList(j) % generate3DMaxwellian(Ti(j), irand)
!             print *, 'Initializing ', particleList(j) % name
!             print *, 'Amount of macroparticles is:', particleList(j) % N_p
!             print *, "Particle mass is:", particleList(j)%mass
!             print *, "Particle charge is:", particleList(j)%q
!             print *, "Particle mean KE is:", particleList(j)%getKEAve(), ", should be", Ti(j) * 1.5
!         end do
        
!         print *, "---------------"
!         print *, ""



!     end function readParticleInputs

    subroutine addMaxwellianLostParticles(particleList, T_e, T_i, irand, world)
        ! Add power to all particles in domain
        type(Particle), intent(in out) :: particleList(2)
        type(Domain), intent(in) :: world
        integer(int32), intent(in out) :: irand
        real(real64), intent(in) :: T_e, T_i
        integer(int32) :: i,j
        do i=1, particleList(1)%delIdx
            call getMaxwellianSample(particleList(1)%phaseSpace(2:4, particleList(1)%N_p + i), particleList(1)%mass, T_e, irand)
            call getMaxwellianFluxSample(particleList(2)%phaseSpace(2:4, particleList(2)%N_p + i), particleList(2)%mass, T_i, irand)
            particleList(1)%phaseSpace(1, particleList(1)%N_p + i) = ran2(irand) * (NumberXNodes - 1) + 1.0d0
            particleList(2)%phaseSpace(1, particleList(2)%N_p + i) = particleList(1)%phaseSpace(1, particleList(1)%N_p + i)
        end do
        particleList(2)%N_p = particleList(2)%N_p + particleList(1)%delIdx
        particleList(1)%N_p = particleList(1)%N_p + particleList(1)%delIdx
        do j = 1, 2
            do i = 1, particleList(j)%refIdx
                if (j == 1) then
                    call getMaxwellianFluxSample(particleList(j)%phaseSpace(2:4, particleList(j)%refRecordIdx(i)), particleList(j)%mass, T_e, irand)
                else
                    call getMaxwellianFluxSample(particleList(j)%phaseSpace(2:4, particleList(j)%refRecordIdx(i)), particleList(j)%mass, T_i, irand)
                end if
                if (world%boundaryConditions(1) == 2) then
                    particleList(j)%phaseSpace(2, particleList(j)%refRecordIdx(i)) = ABS(particleList(j)%phaseSpace(2, particleList(j)%refRecordIdx(i)))
                else
                    particleList(j)%phaseSpace(2, particleList(j)%refRecordIdx(i)) = -ABS(particleList(j)%phaseSpace(2, particleList(j)%refRecordIdx(i)))
                end if
            end do
        end do

    end subroutine addMaxwellianLostParticles
    
    subroutine writePhi(phi, CurrentDiagStep, boolAverage, dirName) 
        ! For diagnostics, deposit single particle density
        ! Re-use rho array since it doesn't get used after first Poisson
        real(real64), intent(in) :: phi(:)
        integer(int32), intent(in) :: CurrentDiagStep
        character(*), intent(in) :: dirName
        character(len=5) :: char_i
        logical, intent(in) :: boolAverage
        write(char_i, '(I3)') CurrentDiagStep
        if (boolAverage) then
            open(41,file='../'//dirName//'/Phi/phi_Average.dat', form='UNFORMATTED')
        else
            open(41,file='../'//dirName//'/Phi/phi_'//trim(adjustl(char_i))//".dat", form='UNFORMATTED')
        end if
        write(41) phi
        close(41)
        
    end subroutine writePhi


    subroutine initialize_randUniform(part, irand)
        ! place particles randomly in each dx_dl based on portion of volume it take up
        type(Particle), intent(in out) :: part
        integer(int32), intent(in out) :: irand
        call getRandom(part%phaseSpace(1,1:part%N_p), irand)
        part%phaseSpace(1,1:part%N_p) = part%phaseSpace(1,1:part%N_p) * (NumberXNodes-1) + 1
        
    end subroutine initialize_randUniform

    ! --------------------------- Diagnostics ------------------------------------




    subroutine loadParticleDensity(densities, particleList, world)
        type(Particle), intent(in) :: particleList(:)
        type(Domain), intent(in) :: world
        real(real64), intent(in out) :: densities(:,:)
        integer(int32) :: i,j, l_left, l_right
        real(real64) :: d
        do i=1, numberChargedParticles
            do j = 1, particleList(i)%N_p
                l_left = INT(particleList(i)%phaseSpace(1,j))
                l_right = l_left + 1
                d = particleList(i)%phaseSpace(1,j) - l_left
                if (l_left > 1) then
                    densities(l_left, i) = densities(l_left, i) + particleList(i)%w_p * (1.0d0-d)
                else
                    SELECT CASE (world%boundaryConditions(l_left))
                    CASE(1)
                        densities(l_left, i) = densities(l_left, i) + particleList(i)%w_p * (1.0d0-d)
                    CASE(2)
                        densities(l_left, i) = densities(l_left, i) + 2.0d0 * particleList(i)%w_p * (1.0d0-d)
                    CASE(3)
                        densities(l_left, i) = densities(l_left, i) + particleList(i)%w_p * (1.0d0-d)
                        densities(NumberXNodes, i) = densities(NumberXNodes, i) + particleList(i)%w_p * (1.0d0-d)
                    END SELECT
                end if
                if (l_right < NumberXNodes) then
                    densities(l_right, i) = densities(l_right, i) + particleList(i)%w_p * d
                else
                    SELECT CASE (world%boundaryConditions(l_right))
                    CASE(1)
                        densities(l_right, i) = densities(l_right, i) + particleList(i)%w_p * d
                    CASE(2)
                        densities(l_right, i) = densities(l_right, i) + 2.0d0 * particleList(i)%w_p * d
                    CASE(3)
                        densities(l_right, i) = densities(l_right, i) + particleList(i)%w_p * d
                        densities(1, i) = densities(1, i) + particleList(i)%w_p * (1.0d0-d)
                    END SELECT
                end if
            end do
        end do

    end subroutine loadParticleDensity

    subroutine WriteParticleDensity(densities, particleList, world, CurrentDiagStep, boolAverage, dirName) 
        ! For diagnostics, deposit single particle density
        ! Re-use rho array since it doesn't get used after first Poisson
        real(real64), intent(in out) :: densities(:,:)
        type(Particle), intent(in) :: particleList(:)
        type(Domain), intent(in) :: world
        character(*), intent(in) :: dirName
        integer(int32), intent(in) :: CurrentDiagStep
        integer(int32) :: i
        logical, intent(in) :: boolAverage
        character(len=5) :: char_i
        
        do i=1, numberChargedParticles
            densities(:,i) = densities(:,i)/world%delX
            write(char_i, '(I3)'), CurrentDiagStep
            if (boolAverage) then
                open(41,file='../'//dirName//'/Density/density_'//particleList(i)%name//"_Average.dat", form='UNFORMATTED')
            else
                open(41,file='../'//dirName//'/Density/density_'//particleList(i)%name//"_"//trim(adjustl(char_i))//".dat", form='UNFORMATTED')
            end if
            write(41) densities(:,i)
            close(41)
        end do
        
    end subroutine WriteParticleDensity


    subroutine generateSaveDirectory(dirName)
        character(*), intent(in) :: dirName
        logical :: bool
        character(len=10) :: buf
        bool = makedirqq('../'//dirName)
        if (.not. bool) then
            print *, "Save directory ", dirName, " already exists. Are you sure you want to continue(yes/no)?"
            read *, buf
            if (buf(1:3) /= 'yes') then
                stop "You have decided to create a new directory for the save files I suppose"
            end if
        else
            bool = makedirqq('../'//dirName//'/Density')
            if (.not. bool) then
                stop "Save directory not successfully created!"
            end if
            bool = makedirqq('../'//dirName//'/ElectronTemperature')
            if (.not. bool) then
                stop "Save directory not successfully created!"
            end if
            bool = makedirqq('../'//dirName//'/PhaseSpace')
            if (.not. bool) then
                stop "Save directory not successfully created!"
            end if
            bool = makedirqq('../'//dirName//'/Phi')
            if (.not. bool) then
                stop "Save directory not successfully created!"
            end if
        end if
    end subroutine generateSaveDirectory

    subroutine solveSimulation(solver, particleList, world, del_t, irand, simulationTime, heatSkipSteps)
        ! Perform certain amount of timesteps, with diagnostics taken at first and last time step
        ! Impliment averaging for final select amount of timeSteps, this will be last data dump
        type(Particle), intent(in out) :: particleList(:)
        type(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t, simulationTime
        integer(int32), intent(in) :: heatSkipSteps
        integer(int32), intent(in out) :: irand
        integer(int32) :: i, j, CurrentDiagStep, diagStepDiff, unitPart1
        real(real64) :: currentTime, densities(NumberXNodes, numberChargedParticles), diagTimeDivision, diagTime, elapsedTime, chargeTotal, energyLoss, elapsed_time
        integer(int64) :: startTime, endTime, timingRate, collisionTime, potentialTime, startTotal, endTotal
        CurrentDiagStep = 1
        unitPart1 = 100
        call generateSaveDirectory(directoryName)
        !Wrtie Initial conditions
        open(15,file='../'//directoryName//'/InitialConditions.dat')
        write(15,'("Number Grid Nodes, Final Expected Time(s), Delta t(s), FractionFreq, delX, n_ave, T_e, T_i, numDiag, NumChargedPart")')
        write(15,"((I3.3, 1x), 7(es16.8,1x), 2(I3.3, 1x))") NumberXNodes, simulationTime, del_t, FractionFreq, world%delX, n_ave, T_e, T_i, numDiagnosticSteps, numberChargedParticles
        close(15)
        call system_clock(count_rate = timingRate)
        ! Write Particle properties
        open(9,file='../'//directoryName//'/ParticleProperties.dat')
        write(9,'("Particle Symbol, Particle Mass (kg), Particle Charge (C), Particle Weight (N/m^2), maxIdx")')
        do j=1, numberChargedParticles
            write(9,"((A, 1x), 3(es16.8,1x), (I6, 1x))") particleList(j)%name//'       ', particleList(j)%mass, particleList(j)%q, particleList(j)%w_p, particleList(j)%finalIdx
        end do
        close(9)
        do i = 1, numberChargedParticles
            particleList(i)%energyLoss = 0.0d0
            particleList(i)%wallLoss = 0
            open(unitPart1+i,file='../'//directoryName//'/ParticleDiagnostic_'//particleList(i)%name//'.dat')
            write(unitPart1+i,'("Time (s), leftCurrentLoss(A/m^2), rightCurrentLoss(A/m^2), leftPowerLoss(W/m^2), rightPowerLoss(W/m^2), N_p")')
        end do
        collisionTime = 0
        potentialTime = 0
        i = 0
        elapsedTime = 0.0d0
        diagStepDiff = 1
        diagTimeDivision = simulationTime/real(numDiagnosticSteps)
        diagTime = diagTimeDivision
        101 format(20(1x,es16.8))
        open(22,file='../'//directoryName//'/GlobalDiagnosticData.dat')
        write(22,'("Time (s), Collision Loss (W/m^2), ParticleCurrentLoss (A/m^2), ParticlePowerLoss(W/m^2)")')
        
        !Save initial particle/field data, along with domain
        inelasticEnergyLoss = 0.0d0
        do i = 1, numberChargedParticles
            particleList(i)%energyLoss = 0.0d0
            particleList(i)%wallLoss = 0
        end do
        call solver%solvePotential(particleList, world)
        !call solver%initialVRewind(particleList, del_t)
        densities = 0.0d0
        call loadParticleDensity(densities, particleList, world)
        call writeParticleDensity(densities, particleList, world, 0, .false., directoryName) 
        call writePhi(solver%phi, 0, .false., directoryName)
        call particleList(1)%writeLocalTemperature(0, directoryName)
        call world%writeDomain(directoryName)
        do j=1, numberChargedParticles
            call particleList(j)%writePhaseSpace(0, directoryName)
        end do
        currentTime = 0.0d0
        call system_clock(startTotal)
        do while(currentTime < simulationTime)
            if (currentTime < diagTime) then
                call system_clock(startTime)
                call solver%moveParticles(particleList, world, del_t)
                call solver%solvePotential(particleList, world)
                call system_clock(endTime)
                potentialTime = potentialTime + (endTime - startTime)
                !call ionizationCollisionIsotropic(particleList(1), particleList(2), 1.0d20, 1.0d-20, del_t, 15.8d0, 0.0d0, irand)
                call system_clock(startTime)
                call addMaxwellianLostParticles(particleList, T_e, T_i, irand, world)
                call system_clock(endTime)
                collisionTime = collisionTime + (endTime - startTime)
            else  
                ! Data dump with diagnostics
                print *, "Simulation is", currentTime/simulationTime * 100.0, "percent done"
                call system_clock(startTime)
                call solver%moveParticles(particleList, world, del_t)
                call solver%solvePotential(particleList, world)
                call system_clock(endTime)
                potentialTime = potentialTime + (endTime - startTime)
                !call ionizationCollisionIsotropic(particleList(1), particleList(2), 1.0d20, 1.0d-20, del_t, 15.8d0, 0.0d0, irand)
                call system_clock(startTime)
                call addMaxwellianLostParticles(particleList, T_e, T_i, irand, world)
                call system_clock(endTime)
                collisionTime = collisionTime + (endTime - startTime)
                densities = 0.0d0
                call loadParticleDensity(densities, particleList, world)
                call writeParticleDensity(densities, particleList, world, CurrentDiagStep, .false., directoryName) 
                call writePhi(solver%phi, CurrentDiagStep, .false., directoryName)
                call particleList(1)%writeLocalTemperature(CurrentDiagStep, directoryName)
                chargeTotal = 0.0d0
                energyLoss = 0.0d0
                do j=1, numberChargedParticles
                    call particleList(j)%writePhaseSpace(CurrentDiagStep, directoryName)
                    chargeTotal = chargeTotal + SUM(particleList(j)%wallLoss) * particleList(j)%q * particleList(j)%w_p
                    energyLoss = energyLoss + SUM(particleList(j)%energyLoss)
                    write(unitPart1+j,"(5(es16.8,1x), (I6,1x))") currentTime, &
                        particleList(j)%wallLoss(1) * particleList(j)%q * particleList(j)%w_p/del_t/diagStepDiff, particleList(j)%wallLoss(2) * particleList(j)%q * particleList(j)%w_p/del_t/diagStepDiff, &
                        particleList(j)%energyLoss(1)/del_t/diagStepDiff, particleList(j)%energyLoss(2)/del_t/diagStepDiff, particleList(j)%N_p
                    particleList(j)%energyLoss = 0.0d0
                    particleList(j)%wallLoss = 0
                end do
                write(22,"(4(es16.8,1x))") currentTime, inelasticEnergyLoss/del_t/diagStepDiff, &
                chargeTotal/del_t/diagStepDiff, energyLoss/del_t/diagStepDiff
                CurrentDiagStep = CurrentDiagStep + 1
                print *, "Number of electrons is:", particleList(1)%N_p
                inelasticEnergyLoss = 0.0d0
                diagStepDiff = 0
                diagTime = diagTime + diagTimeDivision
            end if
            ! if (MODULO(i+1, heatSkipSteps) == 0) then
            !     call addUniformPowerMaxwellian(particleList(1), Power, nu_h, irand, heatSkipSteps*del_t)
            ! end if
            currentTime = currentTime + del_t
            i = i + 1
            diagStepDiff = diagStepDiff + 1
        end do
        call system_clock(startTime)
        call solver%moveParticles(particleList, world, del_t)
        call solver%solvePotential(particleList, world)
        call system_clock(endTime)
        potentialTime = potentialTime + (endTime - startTime)
        !call ionizationCollisionIsotropic(particleList(1), particleList(2), 1.0d20, 1.0d-20, del_t, 15.8d0, 0.0d0, irand)
        call system_clock(startTime)
        call addMaxwellianLostParticles(particleList, T_e, T_i, irand, world)
        call system_clock(endTime)
        collisionTime = collisionTime + (endTime-startTime)
        densities = 0.0d0
        call loadParticleDensity(densities, particleList, world)
        call writeParticleDensity(densities, particleList, world, CurrentDiagStep, .false., directoryName) 
        call writePhi(solver%phi, CurrentDiagStep, .false., directoryName)
        call particleList(1)%writeLocalTemperature(CurrentDiagStep, directoryName)
        chargeTotal = 0.0d0
        energyLoss = 0.0d0
        do j=1, numberChargedParticles
            call particleList(j)%writePhaseSpace(CurrentDiagStep, directoryName)
            chargeTotal = chargeTotal + SUM(particleList(j)%wallLoss) * particleList(j)%q * particleList(j)%w_p
            energyLoss = energyLoss + SUM(particleList(j)%energyLoss)
            write(unitPart1+j,"(5(es16.8,1x), (I6, 1x))") currentTime, &
                particleList(j)%wallLoss(1) * particleList(j)%q * particleList(j)%w_p/del_t/diagStepDiff, particleList(j)%wallLoss(2) * particleList(j)%q * particleList(j)%w_p/del_t/diagStepDiff, &
                particleList(j)%energyLoss(1)/del_t/diagStepDiff, particleList(j)%energyLoss(2)/del_t/diagStepDiff, particleList(j)%N_p
        end do
        write(22,"(4(es16.8,1x))") currentTime, inelasticEnergyLoss/del_t/diagStepDiff, &
        chargeTotal/del_t/diagStepDiff, energyLoss/del_t/diagStepDiff
        close(22)
        call system_clock(endTotal)
        elapsed_time = real((endTotal - startTotal), kind = real64) / real(timingRate, kind = real64)
        ! Write Final Data
        open(9,file='../'//directoryName//'/SimulationFinalData.dat')
        write(9,'("Elapsed Times(s), Potential Time (s), Collision Time (s), Total Steps")')
        write(9,"(3(es16.8,1x), 1(I6, 1x))") elapsed_time, real(potentialTime, kind = real64) / real(timingRate, kind = real64), real(collisionTime, kind = real64) / real(timingRate, kind = real64), i+1
        close(9)
        print *, "Elapsed time for simulation is:", elapsed_time, "seconds"

    end subroutine solveSimulation


    subroutine solveSimulationFinalAverage(solver, particleList, world, del_t, irand, averagingTime, binNumber)
        ! Perform certain amount of timesteps, with diagnostics taken at first and last time step
        ! Impliment averaging for final select amount of timeSteps, this will be last data dump
        type(Particle), intent(in out) :: particleList(:)
        type(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t, averagingTime
        integer(int32), intent(in) :: binNumber
        integer(int32), intent(in out) :: irand
        integer(int32) :: i, stepsAverage, windowNum, windowDivision, j, intPartV, VHist(2*binNumber)
        real(real64) :: phi_average(NumberXNodes), densities(NumberXNodes, numberChargedParticles), currentTime, chargeTotal, energyLoss, meanLoss, stdLoss
        real(real64) :: E_max, VMax
        real(real64), allocatable :: wallLoss(:)


        
        
        !Save initial particle/field data, along with domain
        do i = 1, numberChargedParticles
            particleList(i)%energyLoss = 0.0d0
            particleList(i)%wallLoss = 0
        end do
        inelasticEnergyLoss = 0.0d0
        phi_average = 0.0d0
        densities = 0.0d0
        currentTime = 0.0d0
        i = 0
        windowDivision = INT(200.0d0 / fractionFreq)
        allocate(wallLoss(2 * windowDivision))
        windowNum = 0
        do while(currentTime < averagingTime)
            call solver%moveParticles(particleList, world, del_t)
            call solver%solvePotential(particleList, world)
            !call ionizationCollisionIsotropic(particleList(1), particleList(2), 1.0d20, 1.0d-20, del_t, 15.8d0, 0.0d0, irand)
            call addMaxwellianLostParticles(particleList, T_e, T_i, irand, world)
            call loadParticleDensity(densities, particleList, world)
            phi_average = phi_average + solver%phi
            ! if (MODULO(i, heatSkipSteps) == 0) then
            !     call addUniformPowerMaxwellian(particleList(1), Power, nu_h, irand, heatSksipSteps*del_t)
            ! end if
            i = i + 1
            currentTime = currentTime + del_t
            windowNum = windowNum + 1
            wallLoss(windowNum) = 0.0d0
            do j = 1, numberChargedParticles
                wallLoss(windowNum) = wallLoss(windowNum) + SUM(particleList(j)%energyLoss)
            end do
            wallLoss(windowNum) = wallLoss(windowNum)/currentTime
            if (windowNum > windowDivision) then
                meanLoss = SUM(wallLoss(1:windowNum))/real(windowNum)
                stdLoss = SQRT(SUM( (wallLoss(1:windowNum) - meanLoss)**2 )/real(windowNum))
                if (stdLoss/meanLoss < 1d-5) exit
                windowNum = 0
            end if
        end do
        print *, "Averaging finished over", currentTime, 'simulation time (s)'
        stepsAverage = i
        densities = densities/i
        call writeParticleDensity(densities, particleList, world, 0, .true., directoryName) 
        phi_average = phi_average/stepsAverage
        call writePhi(phi_average, 0, .true., directoryName)
        chargeTotal = 0.0d0
        energyLoss = 0.0d0
        do i=1, numberChargedParticles
            chargeTotal = chargeTotal + SUM(particleList(i)%wallLoss) * particleList(i)%q * particleList(i)%w_p
            energyLoss = energyLoss + SUM(particleList(i)%energyLoss)
        end do
        open(22,file='../'//directoryName//'/GlobalDiagnosticDataAveraged.dat')
        write(22,'("Number Steps, Collision Loss (W/m^2), ParticleCurrentLoss (A/m^2), ParticlePowerLoss(W/m^2)")')
        write(22,"((I6, 1x), 3(es16.8,1x))") stepsAverage, inelasticEnergyLoss/currentTime, chargeTotal/currentTime, energyLoss/currentTime
        close(22)
        print *, "Electron average wall loss:", SUM(particleList(1)%wallLoss)* particleList(1)%w_p * particleList(1)%q/currentTime
        print *, "Ion average wall loss:", SUM(particleList(2)%wallLoss)* particleList(2)%w_p * particleList(2)%q/currentTime
        print *, "Performing average for EEDF over 50/omega_p"
        E_max = 3.0d0 * (MAXVAL(phi_average) - minval(phi_average))
        VMax = SQRT(2.0d0 * E_max *e/ m_e)
        windowDivision = INT(50.0d0/fractionFreq)
        VHist = 0.0d0
        do i = 1, windowDivision
            call solver%moveParticles(particleList, world, del_t)
            call solver%solvePotential(particleList, world)
            call addMaxwellianLostParticles(particleList, T_e, T_i, irand, world)
            do j = 1, particleList(1)%N_p
                intPartV = INT(particleList(1)%phaseSpace(2, j) * (binNumber) / VMax + binNumber + 1)
                VHist(intPartV) = VHist(intPartV) + 1
            end do
        end do
        open(10,file='../'//directoryName//'/ElectronTemperature/EVDF_average.dat', form='UNFORMATTED')
        write(10) real(VHist, kind = real64), VMax
        close(10)

    end subroutine solveSimulationFinalAverage




end module mod_simulation