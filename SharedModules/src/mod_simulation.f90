module mod_simulation
    ! module which actually uses others to run a simulation (single or multiple time steps)
    ! will contain diagnostics as well
    use iso_fortran_env, only: int32, real64, int64
    use constants
    use mod_BasicFunctions
    use mod_particle
    use mod_domain
    use mod_collisions
    use mod_potentialSolver
    use mod_particleMover
    use mod_nonLinSolvers
    use ifport, only: makedirqq
    implicit none

    integer(int32) :: numTimeSteps
    real(real64) :: del_t, energyError, chargeError, gaussError
    


contains


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

    subroutine solveSimulationTest(solver, particleList, world, del_t, maxIter, eps_r, irand, simulationTime, heatSkipSteps)
        ! Perform certain amount of timesteps, with diagnostics taken at first and last time step
        ! Impliment averaging for final select amount of timeSteps, this will be last data dump
        type(Particle), intent(in out) :: particleList(:)
        type(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t, eps_r, simulationTime
        integer(int32), intent(in) :: maxIter, heatSkipSteps
        integer(int32), intent(in out) :: irand
        integer(int32) :: i, j, CurrentDiagStep, k
        integer(int64) :: startTime, endTime, startTotal, endTotal, timingRate
        real(real64) :: currentTime, densities(NumberXNodes, numberChargedParticles), diagTimeDivision, diagTime, Etotal, chargeTotal, elapsed_time, pastDiagTime
        real(real64) :: particleEnergyLossTemp, currDel_t, remainDel_t
        real(real64) :: KE_i, KE_f, PE_i, PE_f
        integer(int64) :: potentialTime, collisionTime, unitPart1
        real(real64) :: rho_i(NumberXNodes)
        CurrentDiagStep = 1
        unitPart1 = 100
        !Wrtie Initial conditions


        call system_clock(count_rate = timingRate)
        ! Write Particle properties
        collisionTime = 0
        potentialTime = 0
        do i = 1, numberChargedParticles
            particleList(i)%energyLoss = 0.0d0
            particleList(i)%wallLoss = 0
        end do
        inelasticEnergyLoss = 0.0d0
        diagTimeDivision = simulationTime/real(numDiagnosticSteps)
        diagTime = diagTimeDivision
       
        !Save initial particle/field data, along with domain
        densities = 0.0d0
        do j=1, numberChargedParticles
            call particleList(j)%loadParticleDensity(densities(:,j), world)
        end do
        currentTime = 0.0d0
        pastDiagTime = 0.0d0
        currDel_t = del_t
        remainDel_t = del_t
        i = 0
        call system_clock(startTotal)
        do while(currentTime < simulationTime) 
            print *, ""
            print *, 'step i =', i
                ! Data dump with diagnostics
            PE_i = solver%getTotalPE(world, .false.)
            KE_i = 0.0d0
            rho_i = 0.0d0
            do j=1, numberChargedParticles
                call particleList(j)%depositRho(rho_i, world)
                KE_i = KE_i + particleList(j)%getTotalKE()
            end do
            call system_clock(startTime)
            call solvePotential(solver, particleList, world, del_t, remainDel_t, currDel_t, maxIter, eps_r)
            call system_clock(endTime)
            KE_f = 0.0d0
            solver%rho = 0.0d0
            do j=1, numberChargedParticles
                call particleList(j)%depositRho(solver%rho, world)
                KE_f = KE_f + particleList(j)%getTotalKE() + SUM(particleList(j)%energyLoss)
            end do
            PE_f = solver%getTotalPE(world, .false.)
            energyError = ABS((KE_i + PE_i - KE_f - PE_f)/(KE_i + PE_i))
            potentialTime = potentialTime + (endTime - startTime)
            !charge conservation directly
            print *, "number electrons after pot solver is:", particleList(1)%N_p
            print *, "number ions after pot solver is:", particleList(2)%N_p
            
            !call ionizationCollisionIsotropic(particleList(1), particleList(2), 1.0d20, 1.0d-20, currDel_t, 15.8d0, 0.0d0, irand)
            call system_clock(startTime)
            call addMaxwellianLostParticles(particleList, T_e, T_i, irand, world)
            call system_clock(endTime)

            collisionTime = collisionTime + (endTime - startTime)
            solver%rho = 0.0d0
            densities = 0.0d0
            do j=1, numberChargedParticles
                call particleList(j)%mergeAndSplit(irand)
                call particleList(j)%loadParticleDensity(densities(:,j), world)
                solver%rho = solver%rho + densities(:, j) * particleList(j)%q
            end do
            chargeError = 0.0d0
            k = 0
            do j = 1, NumberXNodes
                SELECT CASE (world%boundaryConditions(j))
                CASE(0)
                    chargeError = chargeError + (1.0 + currDel_t * (solver%J(j) - solver%J(j-1))/(solver%rho(j) - rho_i(j)))**2
                    k = k + 1
                CASE(2)
                    if (j == 1) then
                        chargeError = chargeError + (1.0 + currDel_t * solver%J(1)/(solver%rho(j) - rho_i(j)))**2
                    else
                        chargeError = chargeError + (1.0 - currDel_t * solver%J(NumberXNodes-1)/(solver%rho(j) - rho_i(j)))**2
                    end if
                    k = k + 1
                CASE(3)
                    chargeError = chargeError + (1.0 + currDel_t * (solver%J(1) - solver%J(NumberXNodes-1))/(solver%rho(j) - rho_i(j)))**2
                    k = k + 1
                END SELECT
            end do
            chargeError = SQRT(chargeError/ k)

            ! Get error gauss' law
            gaussError = solver%getError_tridiag_Poisson(world)
            !chargeError = chargeError / SQRT(SUM(solver%rho**2))
            
            ! Stop program if catch abnormally large error
            print *, 'number of iterations is:', iterNumPicard
            print *, "Energy error is:", energyError
            print *, "Guass error is:", gaussError
            print *, "Charge error is:", chargeError
            print *, 'charge error full is:'
            print *, (1.0 + currDel_t * (solver%J(2:) - solver%J(1:NumberXNodes-2))/(solver%rho(2:NumberXNodes-1) - rho_i(2:NumberXNodes-1)))
            print *, "Number of electrons is:", particleList(1)%N_p
            print *, 'Total number of electrons is:', SUM(particleList(1)%N_p)
            print *, "Number of electrons lost to walls was:", particleList(1)%delIdx
            print *, "Number of electrons refluxed was:", particleList(1)%refIdx
            Etotal = solver%getTotalPE(world, .false.)
            chargeTotal = 0.0d0
            particleEnergyLossTemp = 0.0d0
            do j=1, numberChargedParticles
                chargeTotal = chargeTotal + SUM(particleList(j)%wallLoss) * particleList(j)%q
                particleEnergyLossTemp = particleEnergyLossTemp + SUM(particleList(j)%energyLoss)
                Etotal = Etotal + particleList(j)%getTotalKE()
                particleList(j)%energyLoss = 0.0d0
                particleList(j)%wallLoss = 0.0d0
            end do
            inelasticEnergyLoss = 0.0d0
            pastDiagTime = currentTime + currDel_t
            diagTime = diagTime + diagTimeDivision
            ! if (MODULO(i+1, heatSkipSteps) == 0) then
            !     call addUniformPowerMaxwellian(particleList(1), Power, nu_h, irand, heatSkipSteps*del_t)
            ! end if
            currentTime = currentTime + currDel_t
            i = i + 1
            if (i==400) stop
        end do
    
        

    end subroutine solveSimulationTest

    subroutine solveSimulation(solver, particleList, world, del_t, maxIter, eps_r, irand, simulationTime, heatSkipSteps)
        ! Perform certain amount of timesteps, with diagnostics taken at first and last time step
        ! Impliment averaging for final select amount of timeSteps, this will be last data dump
        type(Particle), intent(in out) :: particleList(:)
        type(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t, eps_r, simulationTime
        integer(int32), intent(in) :: maxIter, heatSkipSteps
        integer(int32), intent(in out) :: irand
        integer(int32) :: i, j, CurrentDiagStep, k
        integer(int64) :: startTime, endTime, startTotal, endTotal, timingRate
        real(real64) :: currentTime, densities(NumberXNodes, numberChargedParticles), diagTimeDivision, diagTime, Etotal, chargeTotal, elapsed_time, pastDiagTime
        real(real64) :: particleEnergyLossTemp, currDel_t, remainDel_t
        real(real64) :: KE_i, KE_f, PE_i, PE_f
        integer(int64) :: potentialTime, collisionTime, unitPart1
        real(real64) :: rho_i(NumberXNodes)
        CurrentDiagStep = 1
        unitPart1 = 100
        call generateSaveDirectory(directoryName)
        !Wrtie Initial conditions
        open(15,file='../'//directoryName//'/InitialConditions.dat')
        write(15,'("Scheme, Number Grid Nodes, T_e, T_i, n_ave, Final Expected Time(s), Delta t(s), FractionFreq, Power(W/m^2), heatSteps, nu_h, numDiag")')
        write(15,"(2(I3.3, 1x), 7(es16.8,1x), (I3.3, 1x), (es16.8,1x), (I3.3, 1x))") schemeNum, NumberXNodes, T_e, T_i, n_ave, simulationTime, del_t, FractionFreq, Power, heatSkipSteps, nu_h, numDiagnosticSteps
        close(15)

        open(15,file='../'//directoryName//'/SolverState.dat')
        write(15,'("Solver Type, eps_r, m_Anderson, Beta_k, maximum iterations")')
        write(15,"(1(I3.3, 1x), 1(es16.8,1x), 1(I3.3, 1x), 1(es16.8,1x), 1(I3.3, 1x))") solverType, eps_r, m_Anderson, Beta_k, maxIter
        close(15)

        call system_clock(count_rate = timingRate)
        ! Write Particle properties
        open(9,file='../'//directoryName//'/ParticleProperties.dat')
        write(9,'("Particle Symbol, Particle Mass (kg), Particle Charge (C)")')
        do j=1, numberChargedParticles
            write(9,"((A, 1x), 2(es16.8,1x))") particleList(j)%name, particleList(j)%mass, particleList(j)%q
        end do
        close(9)
        collisionTime = 0
        potentialTime = 0
        do i = 1, numberChargedParticles
            particleList(i)%energyLoss = 0.0d0
            particleList(i)%wallLoss = 0
            open(unitPart1+i,file='../'//directoryName//'/ParticleDiagnostic_'//particleList(i)%name//'.dat')
            write(unitPart1+i,'("Time (s), leftCurrentLoss(A/m^2), rightCurrentLoss(A/m^2), leftPowerLoss(W/m^2), rightPowerLoss(W/m^2)")')
        end do
        inelasticEnergyLoss = 0.0d0
        diagTimeDivision = simulationTime/real(numDiagnosticSteps)
        diagTime = diagTimeDivision
        open(22,file='../'//directoryName//'/GlobalDiagnosticData.dat')
        write(22,'("Time (s), Collision Loss (W/m^2), ParticleCurrentLoss (A/m^2), ParticlePowerLoss(W/m^2), EnergyTotal (J/m^2), gaussError (a.u), chargeError (a.u), energyError(a.u), Picard Iteration Number")')
        
        !Save initial particle/field data, along with domain
        densities = 0.0d0
        do j=1, numberChargedParticles
            call particleList(j)%loadParticleDensity(densities(:,j), world)
            call particleList(j)%writeParticleDensity(densities(:,j), world, 0, .false., directoryName) 
            call particleList(j)%writePhaseSpace(0, directoryName)
        end do
        call writePhi(solver%phi, 0, .false., directoryName)
        call particleList(1)%writeLocalTemperature(0, directoryName)
        call world%writeDomain(directoryName)
        currentTime = 0.0d0
        pastDiagTime = 0.0d0
        currDel_t = del_t
        remainDel_t = del_t
        call system_clock(startTotal)
        do while(currentTime < simulationTime)
            if (currentTime < diagTime) then
                call system_clock(startTime)
                call solvePotential(solver, particleList, world, del_t, remainDel_t, currDel_t, maxIter, eps_r)
                call system_clock(endTime)
                potentialTime = potentialTime + (endTime - startTime)
                !call ionizationCollisionIsotropic(particleList(1), particleList(2), 1.0d20, 1.0d-20, currDel_t, 15.8d0, 0.0d0, irand)
                call system_clock(startTime)
                call addMaxwellianLostParticles(particleList, T_e, T_i, irand, world)
                call system_clock(endTime)
                collisionTime = collisionTime + (endTime - startTime)
            else  
                ! Data dump with diagnostics
                print *, "Simulation is", currentTime/simulationTime * 100.0, "percent done"
                particleEnergyLossTemp = 0.0d0
                PE_i = solver%getTotalPE(world, .false.)
                KE_i = 0.0d0
                rho_i = 0.0d0
                do j=1, numberChargedParticles
                    call particleList(j)%depositRho(rho_i, world)
                    particleEnergyLossTemp = particleEnergyLossTemp + SUM(particleList(j)%energyLoss)
                    KE_i = KE_i + particleList(j)%getTotalKE()
                end do
                call system_clock(startTime)
                call solvePotential(solver, particleList, world, del_t, remainDel_t, currDel_t, maxIter, eps_r)
                call system_clock(endTime)
                KE_f = -particleEnergyLossTemp
                solver%rho = 0.0d0
                do j=1, numberChargedParticles
                    call particleList(j)%depositRho(solver%rho, world)
                    KE_f = KE_f + particleList(j)%getTotalKE() + SUM(particleList(j)%energyLoss)
                end do
                PE_f = solver%getTotalPE(world, .false.)
                energyError = ABS((KE_i + PE_i - KE_f - PE_f)/(KE_i + PE_i))
                potentialTime = potentialTime + (endTime - startTime)
                !charge conservation directly
                chargeError = 0.0d0
                k = 0
                do j = 1, NumberXNodes
                    SELECT CASE (world%boundaryConditions(j))
                    CASE(0)
                        chargeError = chargeError + (1.0 + currDel_t * (solver%J(j) - solver%J(j-1))/(solver%rho(j) - rho_i(j))/world%nodeVol(j))**2
                        k = k + 1
                    CASE(2)
                        if (j == 1) then
                            chargeError = chargeError + (1.0 + currDel_t * 2.0d0 * solver%J(1)/(solver%rho(j) - rho_i(j))/world%nodeVol(j))**2
                        else
                            chargeError = chargeError + (1.0 - currDel_t * 2.0d0 * solver%J(NumberXNodes-1)/(solver%rho(j) - rho_i(j))/world%nodeVol(j))**2
                        end if
                        k = k + 1
                    CASE(3)
                        chargeError = chargeError + (1.0 + currDel_t * (solver%J(1) - solver%J(NumberXNodes-1))/(solver%rho(j) - rho_i(j))/world%nodeVol(j))**2
                        k = k + 1
                    END SELECT
                end do
                chargeError = SQRT(chargeError/ k)
                
                
                !call ionizationCollisionIsotropic(particleList(1), particleList(2), 1.0d20, 1.0d-20, currDel_t, 15.8d0, 0.0d0, irand)
                call system_clock(startTime)
                call addMaxwellianLostParticles(particleList, T_e, T_i, irand, world)
                call system_clock(endTime)
                collisionTime = collisionTime + (endTime - startTime)
                call writePhi(solver%phi, CurrentDiagStep, .false., directoryName)
                solver%rho = 0.0d0
                densities = 0.0d0
                do j=1, numberChargedParticles
                    call particleList(j)%loadParticleDensity(densities(:,j), world)
                    call particleList(j)%writeParticleDensity(densities(:,j), world, CurrentDiagStep, .false., directoryName) 
                    solver%rho = solver%rho + densities(:, j) * particleList(j)%q
                end do

                ! Get error gauss' law
                gaussError = solver%getError_tridiag_Poisson(world)
                !chargeError = chargeError / SQRT(SUM(solver%rho**2))
                
                ! Stop program if catch abnormally large error
                if (energyError > eps_r) then
                    print *, "-------------------------WARNING------------------------"
                    print *, "Energy error is:", energyError
                    print *, "Total energy not conserved over time step in sub-step procedure!"
                end if
                
                if (gaussError > 1.0d-3) then
                    print *, "-------------------------WARNING------------------------"
                    print *, "Charge error is:", gaussError
                    ! print *, "Re-solve Gauss law"
                    ! call solver%solve_tridiag_Poisson(world)
                end if
                
                call particleList(1)%writeLocalTemperature(CurrentDiagStep, directoryName)
                Etotal = solver%getTotalPE(world, .false.)
                chargeTotal = 0.0d0
                particleEnergyLossTemp = 0.0d0
                do j=1, numberChargedParticles
                    call particleList(j)%writePhaseSpace(CurrentDiagStep, directoryName)
                    chargeTotal = chargeTotal + SUM(particleList(j)%wallLoss) * particleList(j)%q
                    particleEnergyLossTemp = particleEnergyLossTemp + SUM(particleList(j)%energyLoss)
                    Etotal = Etotal + particleList(j)%getTotalKE()
                    write(unitPart1+j,"(5(es16.8,1x))") currentTime + currDel_t, &
                        particleList(j)%wallLoss(1) * particleList(j)%q * particleList(j)%w_p/(currentTime + currDel_t - pastDiagTime), particleList(j)%wallLoss(2) * particleList(j)%q * particleList(j)%w_p/(currentTime + currDel_t - pastDiagTime), &
                        particleList(j)%energyLoss(1)/(currentTime + currDel_t - pastDiagTime), particleList(j)%energyLoss(2)/(currentTime + currDel_t - pastDiagTime)
                    particleList(j)%energyLoss = 0.0d0
                    particleList(j)%wallLoss = 0.0d0
                end do
                write(22,"(8(es16.8,1x), 1(I4, 1x))") currentTime + currDel_t, inelasticEnergyLoss/(currentTime + currDel_t - pastDiagTime), &
                chargeTotal/(currentTime + currDel_t - pastDiagTime), particleEnergyLossTemp/(currentTime + currDel_t - pastDiagTime), Etotal, &
                gaussError, chargeError, energyError, iterNumPicard
                CurrentDiagStep = CurrentDiagStep + 1
                inelasticEnergyLoss = 0.0d0
                print *, "Number of electrons is:", particleList(1)%N_p
                print *, "Number of electrons lost to walls was:", particleList(1)%delIdx
                print *, "Number of electrons refluxed was:", particleList(1)%refIdx
                pastDiagTime = currentTime + currDel_t
                diagTime = diagTime + diagTimeDivision
            end if
            ! if (MODULO(i+1, heatSkipSteps) == 0) then
            !     call addUniformPowerMaxwellian(particleList(1), Power, nu_h, irand, heatSkipSteps*del_t)
            ! end if
            currentTime = currentTime + currDel_t
            i = i + 1
        end do
        
        particleEnergyLossTemp = 0.0d0
        PE_i = solver%getTotalPE(world, .false.)
        KE_i = 0.0d0
        rho_i = 0.0d0
        do j=1, numberChargedParticles
            call particleList(j)%depositRho(rho_i, world)
            particleEnergyLossTemp = particleEnergyLossTemp + SUM(particleList(j)%energyLoss)
            KE_i = KE_i + particleList(j)%getTotalKE()
        end do
        call system_clock(startTime)
        call solvePotential(solver, particleList, world, del_t, remainDel_t, currDel_t, maxIter, eps_r)
        call system_clock(endTime)
        KE_f = -particleEnergyLossTemp
        solver%rho = 0.0d0
        do j=1, numberChargedParticles
            call particleList(j)%depositRho(solver%rho, world)
            KE_f = KE_f + particleList(j)%getTotalKE() + SUM(particleList(j)%energyLoss)
        end do
        PE_f = solver%getTotalPE(world, .false.)
        energyError = ABS((KE_i + PE_i - KE_f - PE_f)/(KE_i + PE_i))
        potentialTime = potentialTime + (endTime - startTime)
        !charge conservation directly
        chargeError = 0.0d0
        k = 0
        do j = 1, NumberXNodes
            SELECT CASE (world%boundaryConditions(j))
            CASE(0)
                chargeError = chargeError + (1.0 + currDel_t * (solver%J(j) - solver%J(j-1))/(solver%rho(j) - rho_i(j))/world%nodeVol(j))**2
                k = k + 1
            CASE(2)
                if (j == 1) then
                    chargeError = chargeError + (1.0 + currDel_t * 2.0d0 * solver%J(1)/(solver%rho(j) - rho_i(j))/world%nodeVol(j))**2
                else
                    chargeError = chargeError + (1.0 - currDel_t * 2.0d0 * solver%J(NumberXNodes-1)/(solver%rho(j) - rho_i(j))/world%nodeVol(j))**2
                end if
                k = k + 1
            CASE(3)
                chargeError = chargeError + (1.0 + currDel_t * (solver%J(1) - solver%J(NumberXNodes-1))/(solver%rho(j) - rho_i(j))/world%nodeVol(j))**2
                k = k + 1
            END SELECT
        end do
        chargeError = SQRT(chargeError/ k)

        !call ionizationCollisionIsotropic(particleList(1), particleList(2), 1.0d20, 1.0d-20, currDel_t, 15.8d0, 0.0d0, irand)
        call system_clock(startTime)
        call addMaxwellianLostParticles(particleList, T_e, T_i, irand, world)
        call system_clock(endTime)
        collisionTime = collisionTime + (endTime-startTime)
        call writePhi(solver%phi, CurrentDiagStep, .false., directoryName)
        solver%rho = 0.0d0
        densities = 0.0d0
        do j=1, numberChargedParticles
            call particleList(j)%loadParticleDensity(densities(:,j), world)
            call particleList(j)%writeParticleDensity(densities(:,j), world, CurrentDiagStep, .false., directoryName) 
            solver%rho = solver%rho + densities(:, j) * particleList(j)%q
        end do

        ! Get error gauss' law
        gaussError = solver%getError_tridiag_Poisson(world)
        !chargeError = chargeError / SQRT(SUM(solver%rho**2))

        ! Stop program if catch abnormally large error
        if (energyError > eps_r) then
            print *, "-------------------------WARNING------------------------"
            print *, "Energy error is:", energyError
            print*, "Total energy not conserved over time step in sub-step procedure!"
        end if
      
        if (gaussError > 1.0d-3) then
            print *, "-------------------------WARNING------------------------"
            print *, "Charge error is:", gaussError
            ! print *, "Re-solve Gauss law"
            ! call solver%solve_tridiag_Poisson(world)
        end if

        call particleList(1)%writeLocalTemperature(CurrentDiagStep, directoryName)
        Etotal = solver%getTotalPE(world, .false.)
        chargeTotal = 0.0d0
        particleEnergyLossTemp = 0.0d0
        do j=1, numberChargedParticles
            call particleList(j)%writePhaseSpace(CurrentDiagStep, directoryName)
            particleEnergyLossTemp = particleEnergyLossTemp + SUM(particleList(j)%energyLoss)
            chargeTotal = chargeTotal + SUM(particleList(j)%wallLoss) * particleList(j)%q
            Etotal = Etotal + particleList(j)%getTotalKE()
            write(unitPart1+j,"(5(es16.8,1x))") currentTime + currDel_t, &
                particleList(j)%wallLoss(1) * particleList(j)%q * particleList(j)%w_p/(currentTime + currDel_t - pastDiagTime), particleList(j)%wallLoss(2) * particleList(j)%q * particleList(j)%w_p/(currentTime + currDel_t - pastDiagTime), &
                particleList(j)%energyLoss(1)/(currentTime + currDel_t - pastDiagTime), particleList(j)%energyLoss(2)/(currentTime + currDel_t - pastDiagTime)
            close(unitPart1+j)
        end do
        write(22,"(8(es16.8,1x), 1(I4, 1x))") currentTime + currDel_t, inelasticEnergyLoss/(currentTime + currDel_t - pastDiagTime), &
        chargeTotal/(currentTime + currDel_t - pastDiagTime), particleEnergyLossTemp/(currentTime + currDel_t - pastDiagTime), &
        Etotal, gaussError, chargeError, energyError, iterNumPicard
        close(22)
        call system_clock(endTotal)
        elapsed_time = real((endTotal - startTotal), kind = real64) / real(timingRate, kind = real64)
        print *, "Elapsed time for simulation is:", elapsed_time, "seconds"
        print *, "Percentage of steps adaptive is:", 100.0d0 * real(amountTimeSplits)/real(i + 1)

        ! Write Particle properties
        open(9,file='../'//directoryName//'/SimulationFinalData.dat')
        write(9,'("Elapsed Times(s), Potential Time (s), Collision Time (s), Total Steps, Number Adaptive Steps")')
        write(9,"(3(es16.8,1x), 2(I6, 1x))") elapsed_time, real(potentialTime, kind = real64) / real(timingRate, kind = real64), real(collisionTime, kind = real64) / real(timingRate, kind = real64), i+1, amountTimeSplits
        close(9)

    end subroutine solveSimulation


    subroutine solveSimulationFinalAverage(solver, particleList, world, del_t, maxIter, eps_r, irand, averagingTime, binNumber)
        ! Perform certain amount of timesteps, with diagnostics taken at first and last time step
        ! Impliment averaging for final select amount of timeSteps, this will be last data dump
        type(Particle), intent(in out) :: particleList(:)
        type(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t, eps_r, averagingTime
        integer(int32), intent(in) :: maxIter, binNumber
        integer(int32), intent(in out) :: irand
        integer(int32) :: i, j, windowNum, VHist(2*binNumber), intPartV
        real(real64) :: phi_average(NumberXNodes), densities(NumberXNodes, numberChargedParticles), currentTime, currDel_t, remainDel_t
        real(real64) :: chargeLossTotal, ELossTotal, lastCheckTime, checkTimeDivision, meanLoss, stdLoss
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
        i = 0
        currentTime = 0.0d0
        currDel_t = del_t
        remainDel_t = del_t
        lastCheckTime = 0.0d0
        checkTimeDivision = 200.0d0 * del_t/fractionFreq
        windowNum = 0
        allocate(wallLoss(2 * INT(checkTimeDivision/del_t)))
        do while(currentTime < averagingTime)
            call solvePotential(solver, particleList, world, del_t, remainDel_t, currDel_t, maxIter, eps_r)
            !call ionizationCollisionIsotropic(particleList(1), particleList(2), 1.0d20, 1.0d-20, currDel_t, 15.8d0, 0.0d0, irand)
            call addMaxwellianLostParticles(particleList, T_e, T_i, irand, world)
            phi_average = phi_average + solver%phi
            ! if (MODULO(i+1, heatSkipSteps) == 0) then
            !     call addUniformPowerMaxwellian(particleList(1), Power, nu_h, irand, heatSkipSteps*del_t)
            ! end if
            currentTime = currentTime + currDel_t
            i = i + 1
            windowNum = windowNum + 1
            wallLoss(windowNum) = 0.0d0
            do j = 1, numberChargedParticles
                call particleList(j)%loadParticleDensity(densities(:,j), world)
                wallLoss(windowNum) = wallLoss(windowNum) + SUM(particleList(j)%energyLoss)
            end do
            wallLoss(windowNum) = wallLoss(windowNum)/currentTime
            if ((currentTime - lastCheckTime) > checkTimeDivision) then
                meanLoss = SUM(wallLoss(1:windowNum))/real(windowNum)
                stdLoss = SQRT(SUM( (wallLoss(1:windowNum) - meanLoss)**2 )/real(windowNum))
                if (stdLoss/meanLoss < 1d-5) exit
                windowNum = 0
                lastCheckTime = currentTime
            end if
        end do
        deallocate(wallLoss)
        print *, "Averaging finished over", currentTime, 'simulation time (s)'
        densities = densities/i
        call writePhi(phi_average/i, 0, .true., directoryName)
        chargeLossTotal = 0.0d0
        ELossTotal = 0.0d0
        do j = 1, numberChargedParticles
            call particleList(j)%writeParticleDensity(densities(:,j), world, 0, .true., directoryName) 
            chargeLossTotal = chargeLossTotal + SUM(particleList(j)%wallLoss) * particleList(j)%q
            ELossTotal = ELossTotal + SUM(particleList(j)%energyLoss)
        end do
        solver%phi = phi_average/i
        solver%rho = 0.0d0
        do j=1, numberChargedParticles
            solver%rho = solver%rho + densities(:, j) * particleList(j)%q
        end do
        call solver%construct_diagMatrix(world)
        gaussError = solver%getError_tridiag_Poisson(world)
        print *, 'gaussError average is:', gaussError
        open(22,file='../'//directoryName//'/GlobalDiagnosticDataAveraged.dat')
        write(22,'("Steps Averaged, Averaging Time, Collision Loss (W/m^2), ParticleCurrentLoss (A/m^2), ParticlePowerLoss(W/m^2), gaussError")')
        write(22,"((I6, 1x), 5(es16.8,1x))") i, currentTime, inelasticEnergyLoss/currentTime, chargeLossTotal/currentTime, ELossTotal/currentTime, gaussError
        close(22)
        print *, "Electron average wall loss:", SUM(particleList(1)%wallLoss)* particleList(1)%q * particleList(1)%w_p/currentTime
        print *, "Ion average wall loss:", SUM(particleList(2)%wallLoss)* particleList(2)%q * particleList(2)%w_p/currentTime
        print *, "Performing average for EEDF over 50/omega_p"
        E_max = 3.0d0 * (MAXVAL(phi_average/i) - minval(phi_average/i))
        VMax = SQRT(2.0d0 * E_max *e/ m_e)
        print *, "V_max is:", VMax
        checkTimeDivision = 50.0d0 * del_t/fractionFreq
        currentTime = 0.0d0
        VHist = 0.0d0
        do while(currentTime < checkTimeDivision)
            call solvePotential(solver, particleList, world, del_t, remainDel_t, currDel_t, maxIter, eps_r)
            !call ionizationCollisionIsotropic(particleList(1), particleList(2), 1.0d20, 1.0d-20, currDel_t, 15.8d0, 0.0d0, irand)
            call addMaxwellianLostParticles(particleList, T_e, T_i, irand, world)
            do j = 1, NumberXNodes-1
                do i=1, particleList(1)%N_p(j)
                    intPartV = INT(particleList(1)%phaseSpace(2, i, j) * (binNumber) / VMax + binNumber + 1)
                    VHist(intPartV) = VHist(intPartV) + 1
                end do
            end do
            currentTime = currentTime + currDel_t
        end do
        open(10,file='../'//directoryName//'/ElectronTemperature/EVDF_average.dat', form='UNFORMATTED')
        write(10) real(VHist, kind = real64), VMax
        close(10)





    end subroutine solveSimulationFinalAverage




end module mod_simulation