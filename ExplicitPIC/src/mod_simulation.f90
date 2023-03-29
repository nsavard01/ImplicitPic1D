module mod_simulation
    ! module which actually uses others to run a simulation (single or multiple time steps)
    ! will contain diagnostics as well
    use iso_fortran_env, only: int32, real64
    use constants
    use mod_BasicFunctions
    use mod_particle
    use mod_domain
    use mod_potentialSolver
    use mod_collisions

    integer(int32) :: numTimeSteps
    real(real64) :: del_t, simulationTime

contains

    ! ------------------------- Reading Input data --------------------------------

    subroutine readInputs(electrons, NumberXNodes, numDiagnosticSteps, stepsAverage, fractionFreq, n_ave, world, solver, simulationTime)
        ! Set initial conditions and global constants based on read input from txt file, create world and solver from these inputs
        integer(int32), intent(in out) :: NumberXNodes, numDiagnosticSteps, stepsAverage
        real(real64), intent(in out) :: fractionFreq, n_ave, simulationTime
        integer(int32) :: io, leftBoundary, rightBoundary
        real(real64) :: leftVoltage, rightVoltage, L_domain, debyeL
        type(Particle) :: electrons
        type(Domain) :: world
        type(potentialSolver) :: solver

        print *, "Reading solver inputs:"
        open(10,file='../InputData/InitialConditions.inp', IOSTAT=io)
        read(10, *, IOSTAT = io) simulationTime
        read(10, *, IOSTAT = io) n_ave
        read(10, *, IOSTAT = io) numDiagnosticSteps
        read(10, *, IOSTAT = io) fractionFreq
        read(10, *, IOSTAT = io) stepsAverage
        close(10)
        print *, "Average initial particle density:", n_ave
        print *, "Number of diagnostic steps is:", numDiagnosticSteps
        print *, "Fraction of 1/w_p for time step:", fractionFreq
        print *, "Number of final steps for averaging:", stepsAverage
        print *, "------------------"
        debyeL = getDebyeLength(electrons%getKEAve()/1.5d0, n_ave)
        print *, ""
        print *, "Reading domain inputs:"
        open(10,file='../InputData/Geometry.inp')
        read(10, *, IOSTAT = io) NumberXNodes
        read(10, *, IOSTAT = io) L_domain
        read(10, *, IOSTAT = io) leftBoundary, rightBoundary
        read(10, *, IOSTAT = io) leftVoltage, rightVoltage
        close(10)
        if (L_domain/debyeL > NumberXNodes) then
            print *, "STARTING DEBYE LENGTH NOT RESOLVED, REDUCTION BEING DONE OF 0.8 DEBYE LENGTH"
            NumberXNodes = CEILING(L_domain/debyeL/0.8d0)
        end if
        print *, "Number of nodes:", NumberXNodes
        print *, "Grid length:", L_domain
        print *, "Left boundary type:", leftBoundary
        print *, "Right boundary type:", rightBoundary
        print *, "------------------"
        print *, ""

        ! if one boundary is periodic, other must also be
        if ((leftBoundary == 3) .or. (rightBoundary == 3)) then
            leftBoundary = 3
            rightBoundary = 3
            leftVoltage = rightVoltage
        end if
        world = Domain(leftBoundary, rightBoundary)
        call world % constructGrid(L_domain)
        solver = potentialSolver(world, leftVoltage, rightVoltage)
        
    end subroutine readInputs

    function readParticleInputs(filename, numberChargedParticles, irand) result(particleList)
        type(Particle), allocatable :: particleList(:)
        character(len=*), intent(in) :: filename
        integer(int32), intent(in out) :: numberChargedParticles, irand
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
                read(10,'(A4)',END=101,ERR=100, ADVANCE = 'NO') name(1:4)
                numSpecies = numSpecies + 1
                read(10,*,END=101,ERR=100) Ti(numSpecies), numParticles(numSpecies), particleIdxFactor(numSpecies)
                mass(numSpecies) = m_e
                charge(numSpecies) = -1.0
                particleNames(numSpecies) = '[e]'
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
                    read(10,*,END=101,ERR=100) mass(numSpecies),charge(numSpecies), Ti(numSpecies), numParticles(numSpecies), particleIdxFactor(numSpecies)
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
            particleList(j) = Particle(mass(j), e * charge(j), 1.0d0, numParticles(j), numParticles(j) * particleIdxFactor(j), trim(particleNames(j)))
            call particleList(j) % generate3DMaxwellian(Ti(j), irand)
            print *, 'Initializing ', particleList(j) % name
            print *, "Particle mass is:", particleList(j)%mass
            print *, "Particle charge is:", particleList(j)%q
            print *, "Particle mean KE is:", particleList(j)%getKEAve(), ", should be", Ti(j) * 1.5
        end do
        
        print *, "---------------"
        print *, ""



    end function readParticleInputs
    
    subroutine writePhi(phi, CurrentDiagStep, boolAverage) 
        ! For diagnostics, deposit single particle density
        ! Re-use rho array since it doesn't get used after first Poisson
        real(real64), intent(in) :: phi(:)
        integer(int32), intent(in) :: CurrentDiagStep
        character(len=5) :: char_i
        logical, intent(in) :: boolAverage
        write(char_i, '(I3)') CurrentDiagStep
        if (boolAverage) then
            open(41,file='../Data/Phi/phi_Average.dat', form='UNFORMATTED')
        else
            open(41,file='../Data/Phi/phi_'//trim(adjustl(char_i))//".dat", form='UNFORMATTED')
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




    subroutine loadParticleDensity(densities, particleList)
        type(Particle), intent(in) :: particleList(:)
        real(real64), intent(in out) :: densities(:,:)
        integer(int32) :: i,j, l_left
        real(real64) :: d
        do i=1, numberChargedParticles
            do j = 1, particleList(i)%N_p
                l_left = INT(particleList(i)%phaseSpace(1,j))
                d = MOD(particleList(i)%phaseSpace(1,j), 1.0d0)
                densities(l_left, i) = densities(l_left, i) + particleList(i)%w_p * (1.0d0-d)
                densities(l_left + 1, i) = densities(l_left + 1, i) + particleList(i)%w_p * d
            end do
        end do

    end subroutine loadParticleDensity

    subroutine WriteParticleDensity(densities, particleList, world, CurrentDiagStep, boolAverage) 
        ! For diagnostics, deposit single particle density
        ! Re-use rho array since it doesn't get used after first Poisson
        real(real64), intent(in out) :: densities(:,:)
        type(Particle), intent(in) :: particleList(:)
        type(Domain), intent(in) :: world
        integer(int32), intent(in) :: CurrentDiagStep
        integer(int32) :: i
        logical, intent(in) :: boolAverage
        character(len=5) :: char_i
        
        do i=1, numberChargedParticles
            densities(:,i) = densities(:,i)/world%delX
            if (world%boundaryConditions(1) == 3) then
                densities(1,i) = densities(1,i) + densities(NumberXNodes, i)
                densities(NumberXNodes, i) = densities(1, i)
            end if
            write(char_i, '(I3)'), CurrentDiagStep
            if (boolAverage) then
                open(41,file='../Data/Density/density_'//particleList(i)%name//"_Average.dat", form='UNFORMATTED')
            else
                open(41,file='../Data/Density/density_'//particleList(i)%name//"_"//trim(adjustl(char_i))//".dat", form='UNFORMATTED')
            end if
            write(41) densities(:,i)
            close(41)
        end do
        
    end subroutine WriteParticleDensity
    
    ! -------------------------- Simulation ------------------------------------------

    subroutine solveSimulationOnlyPotential(solver, particleList, world, del_t, simulationTime)
        ! Perform certain amount of timesteps, with diagnostics taken at first and last time step
        ! Impliment averaging for final select amount of timeSteps, this will be last data dump
        type(Particle), intent(in out) :: particleList(:)
        type(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: simulationTime, del_t
        integer(int32) :: j, CurrentDiagStep
        real(real64) :: currentTime, densities(NumberXNodes, numberChargedParticles), diagTimeDivision, diagTime
        CurrentDiagStep = 1
        currentTime = 0.0d0
        !Wrtie Initial conditions
        open(15,file='../Data/InitialConditions.dat')
        write(15,'("Number Grid Nodes, Final Expected Time(s), Delta t(s), FractionFreq, delX")')
        write(15,"((I3.3, 1x), 4(es16.8,1x))") NumberXNodes, simulationTime, del_t, FractionFreq, world%delX
        close(15)

        ! Write Particle properties
        open(9,file='../Data/ParticleProperties.dat')
        write(9,'("Particle Symbol, Particle Mass (kg), Particle Charge (C), Particle Weight (N/m^2)")')
        do j=1, numberChargedParticles
            write(9,"((A, 1x), 3(es16.8,1x))") particleList(j)%name, particleList(j)%mass, particleList(j)%q, particleList(j)%w_p
        end do
        close(9)

        diagTimeDivision = simulationTime/real(numDiagnosticSteps)
        diagTime = diagTimeDivision
        101 format(20(1x,es16.8))
        open(22,file='../Data/GlobalDiagnosticData.dat')
        write(22,'("Time (s), Collision Loss (W/m^2), ParticleCurrentLoss (A/m^2), ParticlePowerLoss(W/m^2)")')
        
        !Save initial particle/field data, along with domain
        call solver%solvePotential(particleList, world)
        call solver%initialVRewind(particleList, del_t)
        densities = 0.0d0
        call loadParticleDensity(densities, particleList)
        call writeParticleDensity(densities, particleList, world, 0, .false.) 
        call writePhi(solver%phi, 0, .false.)
        call particleList(1)%writeLocalTemperature(0)
        call world%writeDomain()
        do j=1, numberChargedParticles
            call particleList(j)%writePhaseSpace(0)
        end do

        do while(currentTime < simulationTime)
            if (currentTime < diagTime) then
                call solver%moveParticles(particleList, world, del_t)
                call solver%solvePotential(particleList, world)
            else  
                ! Data dump with diagnostics
                print *, "Simulation is", currentTime/simulationTime * 100.0, "percent done"
                inelasticEnergyLoss = 0.0d0
                solver%particleEnergyLoss = 0.0d0
                solver%particleChargeLoss = 0.0d0
                call solver%moveParticles(particleList, world, del_t)
                call solver%solvePotential(particleList, world)
                densities = 0.0d0
                call loadParticleDensity(densities, particleList)
                call writeParticleDensity(densities, particleList, world, CurrentDiagStep, .false.) 
                call writePhi(solver%phi, CurrentDiagStep, .false.)
                call particleList(1)%writeLocalTemperature(CurrentDiagStep)
                do j=1, numberChargedParticles
                    call particleList(j)%writePhaseSpace(CurrentDiagStep)
                end do
                write(22,"(4(es16.8,1x))") currentTime, inelasticEnergyLoss/del_t, &
                SUM(solver%particleChargeLoss)/del_t, solver%particleEnergyLoss/del_t
                CurrentDiagStep = CurrentDiagStep + 1
                diagTime = diagTime + diagTimeDivision
            end if
            currentTime = currentTime + del_t
        end do
        inelasticEnergyLoss = 0.0d0
        solver%particleEnergyLoss = 0.0d0
        solver%particleChargeLoss = 0.0d0
        call solver%moveParticles(particleList, world, del_t)
        call solver%solvePotential(particleList, world)
        densities = 0.0d0
        call loadParticleDensity(densities, particleList)
        call writeParticleDensity(densities, particleList, world, CurrentDiagStep, .false.) 
        call writePhi(solver%phi, CurrentDiagStep, .false.)
        call particleList(1)%writeLocalTemperature(CurrentDiagStep)
        do j=1, numberChargedParticles
            call particleList(j)%writePhaseSpace(CurrentDiagStep)
        end do
        write(22,"(4(es16.8,1x))") currentTime, inelasticEnergyLoss/del_t, &
                SUM(solver%particleChargeLoss)/del_t, solver%particleEnergyLoss/del_t
        close(22)

    end subroutine solveSimulationOnlyPotential

    subroutine solveSimulation(solver, particleList, world, del_t, irand, simulationTime, heatSkipSteps)
        ! Perform certain amount of timesteps, with diagnostics taken at first and last time step
        ! Impliment averaging for final select amount of timeSteps, this will be last data dump
        type(Particle), intent(in out) :: particleList(:)
        type(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t, simulationTime
        integer(int32), intent(in) :: heatSkipSteps
        integer(int32), intent(in out) :: irand
        integer(int32) :: i, j, CurrentDiagStep
        real(real64) :: currentTime, densities(NumberXNodes, numberChargedParticles), diagTimeDivision, diagTime
        CurrentDiagStep = 1

        !Wrtie Initial conditions
        open(15,file='../Data/InitialConditions.dat')
        write(15,'("Number Grid Nodes, Final Expected Time(s), Delta t(s), FractionFreq, delX")')
        write(15,"((I3.3, 1x), 4(es16.8,1x))") NumberXNodes, simulationTime, del_t, FractionFreq, world%delX
        close(15)

        ! Write Particle properties
        open(9,file='../Data/ParticleProperties.dat')
        write(9,'("Particle Symbol, Particle Mass (kg), Particle Charge (C), Particle Weight (N/m^2)")')
        do j=1, numberChargedParticles
            write(9,"((A, 1x), 3(es16.8,1x))") particleList(j)%name, particleList(j)%mass, particleList(j)%q, particleList(j)%w_p
        end do
        close(9)
        i = 0
        diagTimeDivision = simulationTime/real(numDiagnosticSteps)
        diagTime = diagTimeDivision
        101 format(20(1x,es16.8))
        open(22,file='../Data/GlobalDiagnosticData.dat')
        write(22,'("Time (s), Collision Loss (W/m^2), ParticleCurrentLoss (A/m^2), ParticlePowerLoss(W/m^2)")')
        
        !Save initial particle/field data, along with domain
        call solver%solvePotential(particleList, world)
        call solver%initialVRewind(particleList, del_t)
        densities = 0.0d0
        call loadParticleDensity(densities, particleList)
        call writeParticleDensity(densities, particleList, world, 0, .false.) 
        call writePhi(solver%phi, 0, .false.)
        call particleList(1)%writeLocalTemperature(0)
        call world%writeDomain()
        do j=1, numberChargedParticles
            call particleList(j)%writePhaseSpace(0)
        end do
        currentTime = 0.0d0
        do while(currentTime < simulationTime)
            if (currentTime < diagTime) then
                call solver%moveParticles(particleList, world, del_t)
                call solver%solvePotential(particleList, world)
                call ionizationCollisionIsotropic(particleList(1), particleList(2), 1.0d20, 1.0d-20, del_t, 15.8d0, 0.0d0, irand)
            else  
                ! Data dump with diagnostics
                print *, "Simulation is", currentTime/simulationTime * 100.0, "percent done"
                inelasticEnergyLoss = 0.0d0
                solver%particleEnergyLoss = 0.0d0
                solver%particleChargeLoss = 0.0d0
                call solver%moveParticles(particleList, world, del_t)
                call solver%solvePotential(particleList, world)
                call ionizationCollisionIsotropic(particleList(1), particleList(2), 1.0d20, 1.0d-20, del_t, 15.8d0, 0.0d0, irand)
                densities = 0.0d0
                call loadParticleDensity(densities, particleList)
                call writeParticleDensity(densities, particleList, world, CurrentDiagStep, .false.) 
                call writePhi(solver%phi, CurrentDiagStep, .false.)
                call particleList(1)%writeLocalTemperature(CurrentDiagStep)
                do j=1, numberChargedParticles
                    call particleList(j)%writePhaseSpace(CurrentDiagStep)
                end do
                write(22,"(4(es16.8,1x))") currentTime, inelasticEnergyLoss/del_t, &
                SUM(solver%particleChargeLoss)/del_t, solver%particleEnergyLoss/del_t
                CurrentDiagStep = CurrentDiagStep + 1
                print *, "Number of electrons is:", particleList(1)%N_p
                diagTime = diagTime + diagTimeDivision
            end if
            if (MODULO(i+1, heatSkipSteps) == 0) then
                call addUniformPowerMaxwellian(particleList(1), Power, nu_h, irand, heatSkipSteps*del_t)
            end if
            currentTime = currentTime + del_t
            i = i + 1
        end do
        inelasticEnergyLoss = 0.0d0
        solver%particleEnergyLoss = 0.0d0
        solver%particleChargeLoss = 0.0d0
        call solver%moveParticles(particleList, world, del_t)
        call solver%solvePotential(particleList, world)
        call ionizationCollisionIsotropic(particleList(1), particleList(2), 1.0d20, 1.0d-20, del_t, 15.8d0, 0.0d0, irand)
        densities = 0.0d0
        call loadParticleDensity(densities, particleList)
        call writeParticleDensity(densities, particleList, world, CurrentDiagStep, .false.) 
        call writePhi(solver%phi, CurrentDiagStep, .false.)
        call particleList(1)%writeLocalTemperature(CurrentDiagStep)
        do j=1, numberChargedParticles
            call particleList(j)%writePhaseSpace(CurrentDiagStep)
        end do
        write(22,"(4(es16.8,1x))") currentTime, inelasticEnergyLoss/del_t, &
        SUM(solver%particleChargeLoss)/del_t, solver%particleEnergyLoss/del_t

    end subroutine solveSimulation


    subroutine solveSimulationFinalAverage(solver, particleList, world, del_t, irand, stepsAverage, heatSkipSteps)
        ! Perform certain amount of timesteps, with diagnostics taken at first and last time step
        ! Impliment averaging for final select amount of timeSteps, this will be last data dump
        type(Particle), intent(in out) :: particleList(:)
        type(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t
        integer(int32), intent(in) :: stepsAverage, heatSkipSteps
        integer(int32), intent(in out) :: irand
        integer(int32) :: i
        real(real64) :: phi_average(NumberXNodes), densities(NumberXNodes, numberChargedParticles)


        open(22,file='../Data/GlobalDiagnosticDataAveraged.dat')
        write(22,'("Number Steps, Collision Loss (W/m^2), ParticleCurrentLoss (A/m^2), ParticlePowerLoss(W/m^2)")')
        
        !Save initial particle/field data, along with domain
        solver%particleEnergyLoss = 0.0d0
        solver%particleChargeLoss = 0.0d0
        inelasticEnergyLoss = 0.0d0
        phi_average = 0.0d0
        densities = 0.0d0
        do i =1, stepsAverage
            call solver%moveParticles(particleList, world, del_t)
            call solver%solvePotential(particleList, world)
            call ionizationCollisionIsotropic(particleList(1), particleList(2), 1.0d20, 1.0d-20, del_t, 15.8d0, 0.0d0, irand)
            call loadParticleDensity(densities, particleList)
            phi_average = phi_average + solver%phi
            if (MODULO(i, heatSkipSteps) == 0) then
                call addUniformPowerMaxwellian(particleList(1), Power, nu_h, irand, heatSkipSteps*del_t)
            end if
        end do
        densities = densities/stepsAverage
        call writeParticleDensity(densities, particleList, world, 0, .true.) 
        call writePhi(phi_average/stepsAverage, 0, .true.)
        write(22,"((I6, 1x), 3(es16.8,1x))") stepsAverage, inelasticEnergyLoss/del_t/stepsAverage, SUM(solver%particleChargeLoss)/del_t/stepsAverage, solver%particleEnergyLoss/del_t/stepsAverage
        close(22)
        print *, "Electron average wall loss:", SUM(solver%particleChargeLoss(:, 1))/del_t/stepsAverage
        print *, "Ion average wall loss:", SUM(solver%particleChargeLoss(:, 2))/del_t/stepsAverage

    end subroutine solveSimulationFinalAverage




end module mod_simulation