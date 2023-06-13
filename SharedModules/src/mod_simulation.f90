module mod_simulation
    ! module which actually uses others to run a simulation (single or multiple time steps)
    ! will contain diagnostics as well
    use iso_fortran_env, only: int32, real64, int64
    use constants
    use mod_BasicFunctions
    use mod_particle
    use mod_domain
    use mod_potentialSolver
    use mod_particleMover
    use mod_collisions
    use mod_Scheme
    use mod_nonLinSolvers
    implicit none

    integer(int32) :: numTimeSteps
    real(real64) :: del_t, simulationTime, averagingTime


contains

    ! ------------------------- Reading Input data --------------------------------

    subroutine readInputs(NumberXNodes, numDiagnosticSteps, averagingTime, fractionFreq, n_ave, world, solver, simulationTime, Power, heatSkipSteps, nu_h, T_e, GeomFilename, InitFilename)
        ! Set initial conditions and global constants based on read input from txt file, create world and solver from these inputs
        integer(int32), intent(in out) :: NumberXNodes, numDiagnosticSteps, heatSkipSteps
        real(real64), intent(in out) :: fractionFreq, n_ave, simulationTime, Power, nu_h, averagingTime
        real(real64), intent(in) :: T_e
        character(len=*), intent(in) :: GeomFilename, InitFilename
        integer(int32) :: io, leftBoundary, rightBoundary, gridType!, schemeType
        real(real64) :: leftVoltage, rightVoltage, L_domain, debyeLength
        type(Domain) :: world
        type(potentialSolver) :: solver

        print *, "Reading initial inputs:"
        open(10,file='../../SharedModules/InputData/'//InitFilename, IOSTAT=io)
        read(10, *, IOSTAT = io) simulationTime
        read(10, *, IOSTAT = io) n_ave
        read(10, *, IOSTAT = io) numDiagnosticSteps
        read(10, *, IOSTAT = io) fractionFreq
        read(10, *, IOSTAT = io) averagingTime
        read(10, *, IOSTAT = io) Power
        read(10, *, IOSTAT = io) heatSkipSteps
        read(10, *, IOSTAT = io) nu_h
        close(10)
        print *, "Average initial particle density:", n_ave
        print *, "Number of diagnostic steps is:", numDiagnosticSteps
        print *, "Fraction of 1/w_p for time step:", fractionFreq
        print *, "Final averaging time is:", averagingTime
        print *, "Power input (W/m^2):", Power
        print *, "Steps to skip for heating:", heatSkipSteps
        print *, "Heating frequency (Hz):", nu_h
        print *, "------------------"
        print *, ""
        print *, "Reading domain inputs:"
        open(10,file='../../SharedModules/InputData/'//GeomFilename)
        read(10, *, IOSTAT = io) NumberXNodes
        read(10, *, IOSTAT = io) L_domain
        read(10, *, IOSTAT = io) gridType
        read(10, *, IOSTAT = io) leftBoundary, rightBoundary
        read(10, *, IOSTAT = io) leftVoltage, rightVoltage
        close(10)
        if (schemeNum == 0) then
            print *, "Scheme is NGP"
        else
            print *, "Scheme is CIC"
        end if
        print *, "Number of nodes:", NumberXNodes
        print *, "Grid length:", L_domain
        print *, "Grid type is:", gridType
        print *, "Left boundary type:", leftBoundary
        print *, "Right boundary type:", rightBoundary
        print *, "------------------"
        print *, ""
        debyeLength = getDebyeLength(T_e, n_ave)
        ! if one boundary is periodic, other must also be
        if ((leftBoundary == 3) .or. (rightBoundary == 3)) then
            leftBoundary = 3
            rightBoundary = 3
            leftVoltage = rightVoltage
        end if
        world = Domain(leftBoundary, rightBoundary)
        call world % constructGrid(debyeLength, L_domain, gridType)
        solver = potentialSolver(world, leftVoltage, rightVoltage)
        
    end subroutine readInputs

    function readParticleInputs(filename, numberChargedParticles, irand, T_e) result(particleList)
        type(Particle), allocatable :: particleList(:)
        character(len=*), intent(in) :: filename
        integer(int32), intent(in out) :: numberChargedParticles, irand
        real(real64), intent(in out) :: T_e
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
                read(10,'(A4)',END=101,ERR=100, ADVANCE = 'NO') name(1:4)
                numSpecies = numSpecies + 1
                read(10,*,END=101,ERR=100) Ti(numSpecies), numParticles(numSpecies), particleIdxFactor(numSpecies)
                T_e = Ti(numSpecies)
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
    
    ! -------------------------- Simulation ------------------------------------------

    subroutine solveSingleTimeStepDiagnostic(solver, particleList, world, del_t, maxIter, eps_r)
        ! Single time step solver with Divergence of ampere, followed by adding of power, followed by collisions
        type(Particle), intent(in out) :: particleList(:)
        type(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t, eps_r
        integer(int32), intent(in) :: maxIter
        integer(int32) :: j!, k
        real(real64) :: KE_i, KE_f, PE_i, PE_f!, rho_f(NumberXNodes)

        ! Get charge/energy conservation error
        solver%particleEnergyLoss = 0.0d0
        PE_i = solver%getTotalPE(world, .false.)
        KE_i = 0.0d0
        do j=1, numberChargedParticles
            KE_i = KE_i + particleList(j)%getTotalKE()
        end do
        !call solver%depositRho(particleList, world) 
        call solvePotential(solver, particleList, world, del_t, maxIter, eps_r)
        KE_f = solver%particleEnergyLoss
        do j=1, numberChargedParticles
            KE_f = KE_f + particleList(j)%getTotalKE()
        end do
        PE_f = solver%getTotalPE(world, .false.)
        solver%energyError = ABS((KE_i + PE_i - KE_f - PE_f)/(KE_i + PE_i))
        ! call depositRhoDiag(rho_f, particleList, world)
        ! solver%chargeError = 0.0d0
        ! j = 0
        ! if (world%boundaryConditions(1) == 3) then
        !     j = j + 1
        !     solver%chargeError = solver%chargeError + (1 + (solver%J(1) - solver%J(NumberXNodes-1)) *del_t/ world%dx_dl(1)/(rho_f(1) - solver%rho(1)))**2
        ! end if
        ! do k = 1, NumberXNodes -2
        !     if ((rho_f(k+1) - solver%rho(k+1)) /= 0) then
        !         solver%chargeError = solver%chargeError + (1 + (solver%J(k + 1) - solver%J(k)) *del_t/ world%dx_dl(k+1)/(rho_f(k+1) - solver%rho(k+1)))**2
        !         j = j + 1
        !     end if
        ! end do
        ! solver%chargeError = SQRT(solver%chargeError/j)

    end subroutine solveSingleTimeStepDiagnostic

    subroutine solveSimulationOnlyPotential(solver, particleList, world, del_t, maxIter, eps_r, simulationTime)
        ! Perform certain amount of timesteps, with diagnostics taken at first and last time step
        ! Impliment averaging for final select amount of timeSteps, this will be last data dump
        type(Particle), intent(in out) :: particleList(:)
        type(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t, eps_r, simulationTime
        integer(int32), intent(in) :: maxIter
        integer(int32) :: i, j, CurrentDiagStep, diagStepDiff, startTime, endTime, timingRate
        real(real64) :: currentTime, densities(NumberXNodes, numberChargedParticles), diagTimeDivision, diagTime, Etotal, elapsed_time
        real(real64) :: particleEnergyLossTemp, E_initial, particleEnergyLossTotal
        integer(int64) :: totalTime
        CurrentDiagStep = 1
        !Wrtie Initial conditions
        open(15,file='../Data/InitialConditions.dat')
        write(15,'("Number Grid Nodes, Final Expected Time(s), Delta t(s), FractionFreq, Power(W/m^2), heatSteps, nu_h")')
        write(15,"((I3.3, 1x), 4(es16.8,1x), (I3.3, 1x), (es16.8,1x))") NumberXNodes, simulationTime, del_t, FractionFreq, Power, heatSkipSteps, nu_h
        close(15)

        call system_clock(count_rate = timingRate)
        ! Write Particle properties
        open(9,file='../Data/ParticleProperties.dat')
        write(9,'("Particle Symbol, Particle Mass (kg), Particle Charge (C), Particle Weight (N/m^2)")')
        do j=1, numberChargedParticles
            write(9,"((A, 1x), 3(es16.8,1x))") particleList(j)%name, particleList(j)%mass, particleList(j)%q, particleList(j)%w_p
        end do
        close(9)
        totalTime = 0
        i = 0

        E_initial = solver%getTotalPE(world, .false.)
        do j=1, numberChargedParticles
            call particleList(j)%writePhaseSpace(CurrentDiagStep)
            E_initial = E_initial + particleList(j)%getTotalKE()
        end do
        solver%particleEnergyLoss = 0.0d0
        solver%particleChargeLoss = 0.0d0
        particleEnergyLossTotal = 0.0d0
        inelasticEnergyLoss = 0.0d0
        diagStepDiff = 1
        diagTimeDivision = simulationTime/real(numDiagnosticSteps)
        diagTime = diagTimeDivision
        101 format(20(1x,es16.8))
        open(22,file='../Data/GlobalDiagnosticData.dat')
        write(22,'("Time (s), Collision Loss (W/m^2), ParticleCurrentLoss (A/m^2), ParticlePowerLoss(W/m^2), EnergyTotal (J/m^2), chargeError (a.u), energyError(a.u), Picard Iteration Number, diagSteps")')
        
        !Save initial particle/field data, along with domain
        densities = 0.0d0
        call loadParticleDensity(densities, particleList, world)
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
                call system_clock(startTime)
                call solvePotential(solver, particleList, world, del_t, maxIter, eps_r)
                call system_clock(endTime)
                totalTime = totalTime + (endTime - startTime)
            else  
                ! Data dump with diagnostics
                print *, "Simulation is", currentTime/simulationTime * 100.0, "percent done"
                particleEnergyLossTemp = solver%particleEnergyLoss
                call system_clock(startTime)
                call solveSingleTimeStepDiagnostic(solver, particleList, world, del_t, maxIter, eps_r)
                particleEnergyLossTemp = particleEnergyLossTemp + solver%particleEnergyLoss
                particleEnergyLossTotal = particleEnergyLossTotal + particleEnergyLossTemp
                call system_clock(endTime)
                totalTime = totalTime + (endTime - startTime)
                densities = 0.0d0
                call loadParticleDensity(densities, particleList, world)
                call writeParticleDensity(densities, particleList, world, CurrentDiagStep, .false.) 
                call writePhi(solver%phi, CurrentDiagStep, .false.)
                solver%rho = 0.0d0
                do j=1, numberChargedParticles
                    solver%rho = solver%rho + densities(:, j) * particleList(j)%q
                end do

                ! Get error gauss' law
                call solver%construct_diagMatrix(world)
                solver%chargeError = solver%getError_tridiag_Poisson(world)
                solver%chargeError = solver%chargeError / SQRT(SUM(solver%rho**2))
                call solver%construct_diagMatrix_Ampere(world)

                ! Stop program if catch abnormally large error
                if (solver%energyError > eps_r) then
                    print *, "-------------------------WARNING------------------------"
                    print *, "Energy error is:", solver%energyError
                    print *, "Total energy not conserved over time step in sub-step procedure!"
                end if
                
                if (solver%chargeError > 1.0d-4) then
                    print *, "-------------------------WARNING------------------------"
                    print *, "Charge error is:", solver%chargeError
                    stop "Total charge not conserved over time step in sub-step procedure!"
                end if
                
                
                call particleList(1)%writeLocalTemperature(CurrentDiagStep)
                Etotal = solver%getTotalPE(world, .false.)
                do j=1, numberChargedParticles
                    call particleList(j)%writePhaseSpace(CurrentDiagStep)
                    Etotal = Etotal + particleList(j)%getTotalKE()
                end do
                write(22,"(7(es16.8,1x), 2(I4, 1x))") currentTime, inelasticEnergyLoss/del_t/diagStepDiff, &
                SUM(solver%particleChargeLoss)/del_t/diagStepDiff, particleEnergyLossTemp/del_t/diagStepDiff, Etotal, &
                solver%chargeError, solver%energyError, iterNumPicard, diagStepDiff
                CurrentDiagStep = CurrentDiagStep + 1
                solver%particleEnergyLoss = 0.0d0
                solver%particleChargeLoss = 0.0d0
                inelasticEnergyLoss = 0.0d0
                print *, "Number of electrons is:", particleList(1)%N_p
                diagStepDiff = 0
                diagTime = diagTime + diagTimeDivision
            end if
            currentTime = currentTime + del_t
            i = i + 1
            diagStepDiff = diagStepDiff + 1
        end do
        call system_clock(startTime)
        particleEnergyLossTemp = solver%particleEnergyLoss
        call solveSingleTimeStepDiagnostic(solver, particleList, world, del_t, maxIter, eps_r)
        particleEnergyLossTemp = particleEnergyLossTemp + solver%particleEnergyLoss
        particleEnergyLossTotal = particleEnergyLossTotal + particleEnergyLossTemp
        call ionizationCollisionIsotropic(particleList(1), particleList(2), 1.0d20, 1.0d-20, del_t, 15.8d0, 0.0d0, irand)
        call system_clock(endTime)
        totalTime = totalTime + (endTime - startTime)
        densities = 0.0d0
        call loadParticleDensity(densities, particleList, world)
        call writeParticleDensity(densities, particleList, world, CurrentDiagStep, .false.) 
        call writePhi(solver%phi, CurrentDiagStep, .false.)
        solver%rho = 0.0d0
        do j=1, numberChargedParticles
            solver%rho = solver%rho + densities(:, j) * particleList(j)%q
        end do

        ! Get error gauss' law
        call solver%construct_diagMatrix(world)
        solver%chargeError = solver%getError_tridiag_Poisson(world)
        solver%chargeError = solver%chargeError / SQRT(SUM(solver%rho**2))
        call solver%construct_diagMatrix_Ampere(world)

        ! Stop program if catch abnormally large error
        if (solver%energyError > eps_r) then
            print *, "-------------------------WARNING------------------------"
            print *, "Energy error is:", solver%energyError
            stop "Total energy not conserved over time step in sub-step procedure!"
        end if
      
        if (solver%chargeError > 1.0d-4) then
            print *, "-------------------------WARNING------------------------"
            print *, "Charge error is:", solver%chargeError
            stop "Total charge not conserved over time step in sub-step procedure!"
        end if

        call particleList(1)%writeLocalTemperature(CurrentDiagStep)
        Etotal = solver%getTotalPE(world, .false.)
        do j=1, numberChargedParticles
            call particleList(j)%writePhaseSpace(CurrentDiagStep)
            Etotal = Etotal + particleList(j)%getTotalKE()
        end do

        ! See what final charge error is
        print *, "Final error in energy is:", ABS((E_initial - particleEnergyLossTotal - Etotal)/E_initial)

        write(22,"(7(es16.8,1x), 2(I4, 1x))") currentTime, inelasticEnergyLoss/del_t/diagStepDiff, &
        SUM(solver%particleChargeLoss)/del_t/diagStepDiff, particleEnergyLossTemp/del_t/diagStepDiff, &
        Etotal, solver%chargeError, solver%energyError, iterNumPicard, diagStepDiff
        elapsed_time = real(totalTime, kind = real64) / real(timingRate, kind = real64)
        print *, "Elapsed time for simulation is:", elapsed_time, "seconds"
        print *, "Percentage of steps adaptive is:", 100.0d0 * real(amountTimeSplits)/real(i + 1)

        ! Write Particle properties
        open(9,file='../Data/SimulationFinalData.dat')
        write(9,'("Simulation time (s), Total Steps, Number Adaptive Steps")')
        write(9,"(1(es16.8,1x), 2(I6, 1x))") elapsed_time, i+1, amountTimeSplits
        close(9)

    end subroutine solveSimulationOnlyPotential

    subroutine solveSimulation(solver, particleList, world, del_t, maxIter, eps_r, irand, simulationTime, heatSkipSteps)
        ! Perform certain amount of timesteps, with diagnostics taken at first and last time step
        ! Impliment averaging for final select amount of timeSteps, this will be last data dump
        type(Particle), intent(in out) :: particleList(:)
        type(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t, eps_r, simulationTime
        integer(int32), intent(in) :: maxIter, heatSkipSteps
        integer(int32), intent(in out) :: irand
        integer(int32) :: i, j, CurrentDiagStep, diagStepDiff, startTime, endTime, timingRate
        real(real64) :: currentTime, densities(NumberXNodes, numberChargedParticles), diagTimeDivision, diagTime, Etotal, elapsed_time
        real(real64) :: particleEnergyLossTemp
        integer(int64) :: totalTime
        CurrentDiagStep = 1
        !Wrtie Initial conditions
        open(15,file='../Data/InitialConditions.dat')
        write(15,'("Number Grid Nodes, Final Expected Time(s), Delta t(s), FractionFreq, Power(W/m^2), heatSteps, nu_h")')
        write(15,"((I3.3, 1x), 4(es16.8,1x), (I3.3, 1x), (es16.8,1x))") NumberXNodes, simulationTime, del_t, FractionFreq, Power, heatSkipSteps, nu_h
        close(15)

        call system_clock(count_rate = timingRate)
        ! Write Particle properties
        open(9,file='../Data/ParticleProperties.dat')
        write(9,'("Particle Symbol, Particle Mass (kg), Particle Charge (C), Particle Weight (N/m^2)")')
        do j=1, numberChargedParticles
            write(9,"((A, 1x), 3(es16.8,1x))") particleList(j)%name, particleList(j)%mass, particleList(j)%q, particleList(j)%w_p
        end do
        close(9)
        totalTime = 0
        i = 0
        solver%particleEnergyLoss = 0.0d0
        solver%particleChargeLoss = 0.0d0
        inelasticEnergyLoss = 0.0d0
        diagStepDiff = 1
        diagTimeDivision = simulationTime/real(numDiagnosticSteps)
        diagTime = diagTimeDivision
        101 format(20(1x,es16.8))
        open(22,file='../Data/GlobalDiagnosticData.dat')
        write(22,'("Time (s), Collision Loss (W/m^2), ParticleCurrentLoss (A/m^2), ParticlePowerLoss(W/m^2), EnergyTotal (J/m^2), chargeError (a.u), energyError(a.u), Picard Iteration Number, diagSteps")')
        
        !Save initial particle/field data, along with domain
        densities = 0.0d0
        call loadParticleDensity(densities, particleList, world)
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
                call system_clock(startTime)
                call solvePotential(solver, particleList, world, del_t, maxIter, eps_r)
                !call ionizationCollisionIsotropic(particleList(1), particleList(2), 1.0d20, 1.0d-20, del_t, 15.8d0, 0.0d0, irand)
                !call addMaxwellianLostParticles(globalParticleList, T_e, 0.1d0, irand, delIdx, idxReFlux, reFluxMaxIdx, 0.03d0, globalWorld)
                call system_clock(endTime)
                totalTime = totalTime + (endTime - startTime)
            else  
                ! Data dump with diagnostics
                print *, "Simulation is", currentTime/simulationTime * 100.0, "percent done"
                particleEnergyLossTemp = solver%particleEnergyLoss
                call system_clock(startTime)
                call solveSingleTimeStepDiagnostic(solver, particleList, world, del_t, maxIter, eps_r)
                particleEnergyLossTemp = particleEnergyLossTemp + solver%particleEnergyLoss
            
                !call ionizationCollisionIsotropic(particleList(1), particleList(2), 1.0d20, 1.0d-20, del_t, 15.8d0, 0.0d0, irand)
                !call addMaxwellianLostParticles(globalParticleList, T_e, 0.1d0, irand, delIdx, idxReFlux, reFluxMaxIdx, 0.03d0, globalWorld)
                call system_clock(endTime)
                totalTime = totalTime + (endTime - startTime)
                densities = 0.0d0
                call loadParticleDensity(densities, particleList, world)
                call writeParticleDensity(densities, particleList, world, CurrentDiagStep, .false.) 
                call writePhi(solver%phi, CurrentDiagStep, .false.)
                solver%rho = 0.0d0
                do j=1, numberChargedParticles
                    solver%rho = solver%rho + densities(:, j) * particleList(j)%q
                end do

                ! Get error gauss' law
                call solver%construct_diagMatrix(world)
                solver%chargeError = solver%getError_tridiag_Poisson(world)
                solver%chargeError = solver%chargeError / SQRT(SUM(solver%rho**2))
                call solver%construct_diagMatrix_Ampere(world)

                ! Stop program if catch abnormally large error
                if (solver%energyError > eps_r) then
                    print *, "-------------------------WARNING------------------------"
                    print *, "Energy error is:", solver%energyError
                    print *, "Total energy not conserved over time step in sub-step procedure!"
                end if
                
                if (solver%chargeError > 1.0d-4) then
                    print *, "-------------------------WARNING------------------------"
                    print *, "Charge error is:", solver%chargeError
                    stop "Total charge not conserved over time step in sub-step procedure!"
                end if
                
                
                call particleList(1)%writeLocalTemperature(CurrentDiagStep)
                Etotal = solver%getTotalPE(world, .false.)
                do j=1, numberChargedParticles
                    call particleList(j)%writePhaseSpace(CurrentDiagStep)
                    Etotal = Etotal + particleList(j)%getTotalKE()
                end do
                write(22,"(7(es16.8,1x), 2(I4, 1x))") currentTime, inelasticEnergyLoss/del_t/diagStepDiff, &
                SUM(solver%particleChargeLoss)/del_t/diagStepDiff, particleEnergyLossTemp/del_t/diagStepDiff, Etotal, &
                solver%chargeError, solver%energyError, iterNumPicard, diagStepDiff
                CurrentDiagStep = CurrentDiagStep + 1
                solver%particleEnergyLoss = 0.0d0
                solver%particleChargeLoss = 0.0d0
                inelasticEnergyLoss = 0.0d0
                print *, "Number of electrons is:", particleList(1)%N_p
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
        particleEnergyLossTemp = solver%particleEnergyLoss
        call solveSingleTimeStepDiagnostic(solver, particleList, world, del_t, maxIter, eps_r)
        particleEnergyLossTemp = particleEnergyLossTemp + solver%particleEnergyLoss
        
        !call ionizationCollisionIsotropic(particleList(1), particleList(2), 1.0d20, 1.0d-20, del_t, 15.8d0, 0.0d0, irand)
        !call addMaxwellianLostParticles(globalParticleList, T_e, 0.1d0, irand, delIdx, idxReFlux, reFluxMaxIdx, 0.03d0, globalWorld)
        call system_clock(endTime)
        totalTime = totalTime + (endTime - startTime)
        densities = 0.0d0
        call loadParticleDensity(densities, particleList, world)
        call writeParticleDensity(densities, particleList, world, CurrentDiagStep, .false.) 
        call writePhi(solver%phi, CurrentDiagStep, .false.)
        solver%rho = 0.0d0
        do j=1, numberChargedParticles
            solver%rho = solver%rho + densities(:, j) * particleList(j)%q
        end do

        ! Get error gauss' law
        call solver%construct_diagMatrix(world)
        solver%chargeError = solver%getError_tridiag_Poisson(world)
        solver%chargeError = solver%chargeError / SQRT(SUM(solver%rho**2))
        call solver%construct_diagMatrix_Ampere(world)

        ! Stop program if catch abnormally large error
        if (solver%energyError > eps_r) then
            print *, "-------------------------WARNING------------------------"
            print *, "Energy error is:", solver%energyError
            print*, "Total energy not conserved over time step in sub-step procedure!"
        end if
      
        if (solver%chargeError > 1.0d-4) then
            print *, "-------------------------WARNING------------------------"
            print *, "Charge error is:", solver%chargeError
            stop "Total charge not conserved over time step in sub-step procedure!"
        end if

        call particleList(1)%writeLocalTemperature(CurrentDiagStep)
        Etotal = solver%getTotalPE(world, .false.)
        do j=1, numberChargedParticles
            call particleList(j)%writePhaseSpace(CurrentDiagStep)
            Etotal = Etotal + particleList(j)%getTotalKE()
        end do
        write(22,"(7(es16.8,1x), 2(I4, 1x))") currentTime, inelasticEnergyLoss/del_t/diagStepDiff, &
        SUM(solver%particleChargeLoss)/del_t/diagStepDiff, particleEnergyLossTemp/del_t/diagStepDiff, &
        Etotal, solver%chargeError, solver%energyError, iterNumPicard, diagStepDiff
        elapsed_time = real(totalTime, kind = real64) / real(timingRate, kind = real64)
        print *, "Elapsed time for simulation is:", elapsed_time, "seconds"
        print *, "Percentage of steps adaptive is:", 100.0d0 * real(amountTimeSplits)/real(i + 1)

        ! Write Particle properties
        open(9,file='../Data/SimulationFinalData.dat')
        write(9,'("Simulation time (s), Total Steps, Number Adaptive Steps")')
        write(9,"(1(es16.8,1x), 2(I6, 1x))") elapsed_time, i+1, amountTimeSplits
        close(9)

    end subroutine solveSimulation


    subroutine solveSimulationFinalAverage(solver, particleList, world, del_t, maxIter, eps_r, irand, averagingTime, heatSkipSteps)
        ! Perform certain amount of timesteps, with diagnostics taken at first and last time step
        ! Impliment averaging for final select amount of timeSteps, this will be last data dump
        type(Particle), intent(in out) :: particleList(:)
        type(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t, eps_r, averagingTime
        integer(int32), intent(in) :: maxIter, heatSkipSteps
        integer(int32), intent(in out) :: irand
        integer(int32) :: i
        real(real64) :: phi_average(NumberXNodes), densities(NumberXNodes, numberChargedParticles), currentTime
        
        !Save initial particle/field data, along with domain
        solver%particleEnergyLoss = 0.0d0
        solver%particleChargeLoss = 0.0d0
        inelasticEnergyLoss = 0.0d0
        phi_average = 0.0d0
        densities = 0.0d0
        i = 0
        currentTime = 0.0d0
        do while(currentTime < averagingTime)
            call solvePotential(solver, particleList, world, del_t, maxIter, eps_r)
            !call ionizationCollisionIsotropic(particleList(1), particleList(2), 1.0d20, 1.0d-20, del_t, 15.8d0, 0.0d0, irand)
            !call addMaxwellianLostParticles(globalParticleList, T_e, 0.1d0, irand, delIdx, idxReFlux, reFluxMaxIdx, 0.03d0, globalWorld)
            call loadParticleDensity(densities, particleList, world)
            phi_average = phi_average + solver%phi
            ! if (MODULO(i+1, heatSkipSteps) == 0) then
            !     call addUniformPowerMaxwellian(particleList(1), Power, nu_h, irand, heatSkipSteps*del_t)
            ! end if
            currentTime = currentTime + del_t
            i = i + 1
        end do
        densities = densities/i
        call writeParticleDensity(densities, particleList, world, 0, .true.) 
        call writePhi(phi_average/i, 0, .true.)
        open(22,file='../Data/GlobalDiagnosticDataAveraged.dat')
        write(22,'("Steps Averaged, Averaging Time, Collision Loss (W/m^2), ParticleCurrentLoss (A/m^2), ParticlePowerLoss(W/m^2)")')
        write(22,"((I6, 1x), 4(es16.8,1x))") i, currentTime, inelasticEnergyLoss/currentTime, SUM(solver%particleChargeLoss)/currentTime, solver%particleEnergyLoss/currentTime
        close(22)
        print *, "Electron average wall loss:", SUM(solver%particleChargeLoss(:, 1))/currentTime
        print *, "Ion average wall loss:", SUM(solver%particleChargeLoss(:, 2))/currentTime






    end subroutine solveSimulationFinalAverage




end module mod_simulation