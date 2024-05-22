module mod_simulation
    ! module which actually uses others to run a simulation (single or multiple time steps)
    ! will contain diagnostics as well
    use iso_fortran_env, only: int32, real64, int64
    use constants
    use mod_BasicFunctions
    use mod_particle
    use mod_targetParticle
    use mod_NullCollision
    use mod_domain
    use mod_potentialSolver
    use mod_particleMover
    use mod_particleInjection
    use mod_Scheme
    use mod_nonLinSolvers
    use ifport, only: makedirqq
    implicit none

    integer(int32) :: numTimeSteps, heatSkipSteps
    real(real64), private :: energyError, chargeError, gaussError, Power, nu_h, inelasticEnergyLoss
    real(real64) :: currentTime


contains

    ! ------------------------- Reading Input data --------------------------------
    
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
            open(41,file=dirName//'/Phi/phi_Average.dat', form='UNFORMATTED')
        else
            open(41,file=dirName//'/Phi/phi_'//trim(adjustl(char_i))//".dat", form='UNFORMATTED')
        end if
        write(41) phi
        close(41)
        
    end subroutine writePhi
    

    subroutine generateSaveDirectory(dirName)
        character(*), intent(in) :: dirName
        logical :: bool
        integer(int32) :: io
        character(len=10) :: buf
        bool = makedirqq(dirName)
        if (.not. bool) then
            print *, "Save directory ", dirName, " already exists. Are you sure you want to continue(yes/no)?"
            read *, buf
            if (buf(1:3) /= 'yes') then
                stop "You have decided to create a new directory for the save files I suppose"
            end if
        else
            bool = makedirqq(dirName//'/Density')
            if (.not. bool) then
                stop "Save directory not successfully created!"
            end if
            bool = makedirqq(dirName//'/ElectronTemperature')
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
            bool = makedirqq(dirName//'/BinaryCollisions')
            if (.not. bool) then
                stop "Save directory not successfully created!"
            end if
        end if
        call execute_command_line("cp -Tr ../InputData "//dirName//"/InputData", EXITSTAT = io)
        if (io /= 0) then
            stop 'Issue copying input data deck'
        end if
    end subroutine generateSaveDirectory

    subroutine solveSimulation(solver, particleList, targetParticleList, nullCollisionList, world, del_t, irand, simulationTime)
        ! Perform certain amount of timesteps, with diagnostics taken at first and last time step
        ! Impliment averaging for final select amount of timeSteps, this will be last data dump
        type(Particle), intent(in out) :: particleList(numberChargedParticles)
        type(potentialSolver), intent(in out) :: solver
        type(targetParticle), intent(in) :: targetParticleList(numberNeutralParticles)
        type(nullCollision), intent(in out) :: nullCollisionList(numberBinaryCollisions)
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t, simulationTime
        integer(int32), intent(in out) :: irand(numThread)
        integer(int32) :: i, j, CurrentDiagStep, smoothInt
        integer(int64) :: startTimer, endTimer, startTotal, endTotal, timingRate
        real(real64) :: diagTimeDivision, diagTime, Etotal, chargeTotal, elapsed_time, pastDiagTime, energyLoss
        real(real64) :: currDel_t, remainDel_t
        real(real64) :: KE_i, KE_f, PE_i, PE_f, momentum_total(3)
        integer(int64) :: potentialTime, collisionTime, unitPart1
        character(len=8) :: date_char
        character(len=10) :: time_char
        real(real64) :: rho_i(NumberXNodes)
        allocate(energyAddColl(numThread))
        CurrentDiagStep = 1
        unitPart1 = 100
        call generateSaveDirectory(directoryName)
        !Wrtie Initial conditions
        smoothInt = 0
        if (world%gridSmoothBool) smoothInt = 1
        call date_and_time(DATE=date_char, TIME = time_char)
        open(15,file=directoryName//'/DateTime.dat')
        write(15,'("UTC Date, UTC Time")')
        write(15,"(2(A,1x))") date_char, time_char 
        close(15)
        open(15,file=directoryName//'/InitialConditions.dat')
        write(15,'("Scheme, Number Grid Nodes, T_e, T_i, n_ave, Final Expected Time(s), Delta t(s), FractionFreq, Power(W/m^2), smoothInt, nu_h, numDiag, numThread, RF_rad_frequency, RF_half_amplitude")')
        write(15,"(2(I6, 1x), 7(es16.8,1x), (I6, 1x), (es16.8,1x), 2(I6, 1x), 2(es16.8,1x))") schemeNum, NumberXNodes, T_e, T_i, n_ave, simulationTime, del_t, FractionFreq, Power, smoothInt, nu_h, numDiagnosticSteps, numThread, solver%RF_rad_frequency, solver%RF_half_amplitude
        close(15)

        call writeSolverState(directoryName)
        call writeParticleInjectionInputs(directoryName)
        call system_clock(count_rate = timingRate)
        ! Write Particle properties
        open(9,file=directoryName//'/ParticleProperties.dat')
        write(9,'("Particle Symbol, Particle Mass (kg), Particle Charge (C), Particle Weight (N/m^2), maxIdx")')
        do j=1, numberChargedParticles
            write(9,"((A, 1x), 3(es16.8,1x), (I10, 1x))") particleList(j)%name, particleList(j)%mass, particleList(j)%q, particleList(j)%w_p, particleList(j)%finalIdx
        end do
        close(9)
        ! Write Neutral particle properties
        open(9,file=directoryName//'/TargetProperties.dat')
        write(9,'("Target Symbol, Target Mass (kg), Target Density (m^-3), Target Temp (K)")')
        do j=1, numberNeutralParticles
            write(9,"((A, 1x), 3(es16.8,1x))") targetParticleList(j)%name, targetParticleList(j)%mass, targetParticleList(j)%density, targetParticleList(j)%temperature
        end do
        close(9)
        do i = 1, numberBinaryCollisions
            call nullCollisionList(i)%writeCollisionProperties(particleList, targetParticleList, directoryName)
        end do
        collisionTime = 0
        potentialTime = 0
        momentum_total = 0.0d0
        do i = 1, numberChargedParticles
            particleList(i)%accumEnergyLoss = 0.0d0
            particleList(i)%accumWallLoss = 0
            particleList(i)%momentumLoss = 0.0d0
            momentum_total = momentum_total + particleList(i)%getTotalMomentum()
            open(unitPart1+i,file=directoryName//'/ParticleDiagnostic_'//particleList(i)%name//'.dat')
            write(unitPart1+i,'("Time (s), leftCurrentLoss(A/m^2), rightCurrentLoss(A/m^2), leftPowerLoss(W/m^2), rightPowerLoss(W/m^2), Number Particles, Temperature (eV), Ave. Num. SubSteps., Ave. Num. Func. Evals.")')
        end do
        inelasticEnergyLoss = 0.0d0
        energyAddColl = 0.0d0
        diagTimeDivision = (simulationTime - startSimulationTime)/real(numDiagnosticSteps)
        diagTime = diagTimeDivision + startSimulationTime
        open(22,file=directoryName//'/GlobalDiagnosticData.dat')
        write(22,'("Time (s), Collision Loss (W/m^2), ParticleCurrentLoss (A/m^2), ParticlePowerLoss(W/m^2), EnergyTotal (J/m^2), MomentumTotal (kg/s/m), gaussError (a.u), chargeError (a.u), energyError(a.u), Picard Iteration Number")')
        
        !Save initial particle/field data, along with domain
        call loadParticleDensity(particleList, world, .true.)
        call writeParticleDensity(particleList, world, 0, .false., directoryName) 
        call writePhi(solver%phi, 0, .false., directoryName)
        call world%writeDomain(directoryName)
        do j=1, numberChargedParticles
            call particleList(j)%writeLocalTemperature(0, directoryName, MIN(NumberXNodes, NumberXHalfNodes))
            ! call particleList(j)%writePhaseSpace(0, directoryName)
        end do
        currentTime = startSimulationTime
        pastDiagTime = startSimulationTime
        currDel_t = del_t
        remainDel_t = del_t
        i = 0
        call system_clock(startTotal)
        do while(currentTime < simulationTime)
            if (currentTime < diagTime) then
                call system_clock(startTimer)
                call solvePotential(solver, particleList, world, del_t, remainDel_t, currDel_t, currentTime)
                call system_clock(endTimer)
                potentialTime = potentialTime + (endTimer - startTimer)
                !call ionizationCollisionIsotropic(particleList(1), particleList(2), 1.0d20, 1.0d-20, currDel_t, 15.8d0, 0.0d0, irand)
                call system_clock(startTimer)
                if (heatingBool) call maxwellianHeating(particleList(1), irand, fractionFreq, T_e, currDel_t, del_t)
                if (addLostPartBool) call addMaxwellianLostParticles(particleList, T_e, T_i, irand, world)
                if (refluxPartBool) call refluxParticles(particleList, T_e, T_i, irand, world)
                if (injectionBool) call injectAtBoundary(particleList, T_e, T_i, irand, world, currDel_t, solver%BFieldAngle)
                if (uniformInjectionBool) call injectUniformFlux(particleList, T_e, T_i, irand, world)
                do j = 1, numberBinaryCollisions
                    call nullCollisionList(j)%generateCollision(particleList, targetParticleList, numberChargedParticles, numberBinaryCollisions, irand, currDel_t)
                end do
                call system_clock(endTimer)
                collisionTime = collisionTime + (endTimer - startTimer)
                
            else  
                ! Data dump with diagnostics
                print *, "Simulation is", (currentTime-startSimulationTime)/(simulationTime - startSimulationTime) * 100.0, "percent done"
                KE_i = 0.0d0
                do j=1, numberChargedParticles
                    KE_i = KE_i + particleList(j)%getTotalKE()
                end do
                call solver%depositRho(particleList, world)
                rho_i = solver%rho
                call system_clock(startTimer)
                call solvePotential(solver, particleList, world, del_t, remainDel_t, currDel_t, currentTime)
                call system_clock(endTimer)
                PE_i = solver%getTotalPE(world, .false.)
                KE_f = 0.0d0
                do j=1, numberChargedParticles
                    KE_f = KE_f + particleList(j)%getTotalKE() + SUM(particleList(j)%energyLoss) * particleList(j)%mass * particleList(j)%w_p * 0.5d0
                end do
                PE_f = solver%getTotalPE(world, .true.) - solver%getEnergyFromBoundary(world, currDel_t) 
                energyError = ABS((KE_i + PE_i - KE_f - PE_f)/(KE_i + PE_i))
                potentialTime = potentialTime + (endTimer - startTimer)
                
                !call ionizationCollisionIsotropic(particleList(1), particleList(2), 1.0d20, 1.0d-20, currDel_t, 15.8d0, 0.0d0, irand)
                call system_clock(startTimer)
                if (heatingBool) call maxwellianHeating(particleList(1), irand, fractionFreq, T_e, currDel_t, del_t)
                if (addLostPartBool) call addMaxwellianLostParticles(particleList, T_e, T_i, irand, world)
                if (refluxPartBool) call refluxParticles(particleList, T_e, T_i, irand, world)
                if (injectionBool) call injectAtBoundary(particleList, T_e, T_i, irand, world, currDel_t, solver%BFieldAngle)
                if (uniformInjectionBool) call injectUniformFlux(particleList, T_e, T_i, irand, world)
                do j = 1, numberBinaryCollisions
                    call nullCollisionList(j)%generateCollision(particleList, targetParticleList, numberChargedParticles, numberBinaryCollisions, irand, currDel_t)
                end do
                call system_clock(endTimer)
                collisionTime = collisionTime + (endTimer - startTimer)
                call loadParticleDensity(particleList, world, .true.)
                call writeParticleDensity(particleList, world, CurrentDiagStep, .false., directoryName) 
                call writePhi(solver%phi_f, CurrentDiagStep, .false., directoryName)
                call solver%depositRho(particleList, world)

                !charge conservation directly
                chargeError = solver%getChargeContinuityError(rho_i, world, currDel_t)

                ! Get error gauss' law
                gaussError = solver%getError_tridiag_Poisson(world)
                !chargeError = chargeError / SQRT(SUM(solver%rho**2))
                
                ! Stop program if catch abnormally large error
                if (energyError > 1.0d-5) then
                    print *, "-------------------------WARNING------------------------"
                    print *, "Energy error is:", energyError
                end if
                
                if (gaussError > 1.0d-5) then
                    print *, "-------------------------WARNING------------------------"
                    print *, "Gauss error is:", gaussError
                    ! print *, "Re-solve Gauss law"
                    ! call solver%solve_tridiag_Poisson(world)
                end if
                
                Etotal = solver%getTotalPE(world, .true.)
                chargeTotal = 0.0d0
                energyLoss = 0.0d0
                momentum_total = 0.0d0
                do j=1, numberChargedParticles
                    call particleList(j)%writeLocalTemperature(CurrentDiagStep, directoryName, MIN(NumberXNodes, NumberXHalfNodes))
                    call particleList(j)%writePhaseSpace(directoryName)
                    momentum_total = momentum_total + particleList(j)%getTotalMomentum()
                    momentum_total(1) = momentum_total(1) + SUM(particleList(j)%momentumLoss(:,:)) * particleList(j)%mass * particleList(j)%w_p
                    chargeTotal = chargeTotal + SUM(particleList(j)%accumWallLoss) * particleList(j)%q * particleList(j)%w_p
                    Etotal = Etotal + particleList(j)%getTotalKE()
                    energyLoss = energyLoss + SUM(particleList(j)%accumEnergyLoss) * particleList(j)%w_p * particleList(j)%mass * 0.5d0
                    write(unitPart1+j,"(5(es16.8,1x), (I10, 1x), 3(es16.8,1x))") currentTime + currDel_t, &
                        particleList(j)%accumWallLoss(1) * particleList(j)%q * particleList(j)%w_p/(currentTime + currDel_t - pastDiagTime), particleList(j)%accumWallLoss(2) * particleList(j)%q * particleList(j)%w_p/(currentTime + currDel_t - pastDiagTime), &
                        particleList(j)%accumEnergyLoss(1)* particleList(j)%mass * particleList(j)%w_p * 0.5d0/(currentTime + currDel_t - pastDiagTime), &
                        particleList(j)%accumEnergyLoss(2) * particleList(j)%mass * particleList(j)%w_p * 0.5d0/(currentTime + currDel_t - pastDiagTime), SUM(particleList(j)%N_p), &
                        particleList(j)%getKEAve()*2.0d0/3.0d0, particleList(j)%numSubStepsAve, particleList(j)%numFuncEvalAve
                    particleList(j)%accumEnergyLoss = 0.0d0
                    particleList(j)%accumWallLoss = 0.0d0
                end do
                do j = 1, numberBinaryCollisions
                    inelasticEnergyLoss = inelasticEnergyLoss + nullCollisionList(j)%totalEnergyLoss * e * particleList(nullCollisionList(j)%reactantsIndx(1))%w_p
                end do
                write(22,"(9(es16.8,1x), 1(I4, 1x))") currentTime + currDel_t, inelasticEnergyLoss/(currentTime + currDel_t - pastDiagTime), &
                chargeTotal/(currentTime + currDel_t - pastDiagTime), energyLoss/(currentTime + currDel_t - pastDiagTime), Etotal, momentum_total(1), &
                gaussError, chargeError, energyError, iterNumPicard
                do j = 1, numberBinaryCollisions
                    nullCollisionList(j)%totalEnergyLoss = 0
                    nullCollisionList(j)%totalAmountCollisions = 0
                end do
                print *, 'Total wall power loss:', energyLoss/(currentTime + currDel_t - pastDiagTime)
                print *, 'Total collision energy loss:', inelasticEnergyLoss/(currentTime + currDel_t - pastDiagTime)
                CurrentDiagStep = CurrentDiagStep + 1
                inelasticEnergyLoss = 0.0d0
                energyAddColl = 0.0d0
                print *, "Number of electrons is:", SUM(particleList(1)%N_p)
                print *, "Number of ions is:", SUM(particleList(2)%N_p)
                print *, "Percentage of steps adaptive is:", 100.0d0 * real(amountTimeSplits)/real(i + 1)
                pastDiagTime = currentTime + currDel_t
                diagTime = diagTime + diagTimeDivision
            end if
            ! if (MODULO(i+1, heatSkipSteps) == 0) then
            !     call addUniformPowerMaxwellian(particleList(1), Power, nu_h, irand, heatSkipSteps*del_t)
            ! end if
            currentTime = currentTime + currDel_t
            i = i + 1
        end do
        
        KE_i = 0.0d0
        do j=1, numberChargedParticles
            KE_i = KE_i + particleList(j)%getTotalKE()
        end do
        call solver%depositRho(particleList, world)
        rho_i = solver%rho
        call system_clock(startTimer)
        call solvePotential(solver, particleList, world, del_t, remainDel_t, currDel_t, currentTime)
        call system_clock(endTimer)
        PE_i = solver%getTotalPE(world, .false.)
        KE_f = 0.0d0
        do j=1, numberChargedParticles
            KE_f = KE_f + particleList(j)%getTotalKE() + SUM(particleList(j)%energyLoss) * particleList(j)%mass * particleList(j)%w_p * 0.5d0
        end do
        PE_f = solver%getTotalPE(world, .true.) - solver%getEnergyFromBoundary(world, currDel_t) 
        energyError = ABS((KE_i + PE_i - KE_f - PE_f)/(KE_i + PE_i))
        potentialTime = potentialTime + (endTimer - startTimer)

        !call ionizationCollisionIsotropic(particleList(1), particleList(2), 1.0d20, 1.0d-20, currDel_t, 15.8d0, 0.0d0, irand)
        call system_clock(startTimer)
        if (heatingBool) call maxwellianHeating(particleList(1), irand, fractionFreq, T_e, currDel_t, del_t)
        if (addLostPartBool) call addMaxwellianLostParticles(particleList, T_e, T_i, irand, world)
        if (refluxPartBool) call refluxParticles(particleList, T_e, T_i, irand, world)
        if (injectionBool) call injectAtBoundary(particleList, T_e, T_i, irand, world, currDel_t, solver%BFieldAngle)
        if (uniformInjectionBool) call injectUniformFlux(particleList, T_e, T_i, irand, world)
        do j = 1, numberBinaryCollisions
            call nullCollisionList(j)%generateCollision(particleList, targetParticleList, numberChargedParticles, numberBinaryCollisions, irand, currDel_t)
        end do
        call system_clock(endTimer)
        collisionTime = collisionTime + (endTimer-startTimer)
        call loadParticleDensity(particleList, world, .true.)
        call writeParticleDensity(particleList, world, CurrentDiagStep, .false., directoryName) 
        call writePhi(solver%phi_f, CurrentDiagStep, .false., directoryName)
        call solver%depositRho(particleList, world)

        chargeError = solver%getChargeContinuityError(rho_i, world, currDel_t)

        ! Get error gauss' law
        gaussError = solver%getError_tridiag_Poisson(world)
        !chargeError = chargeError / SQRT(SUM(solver%rho**2))

        ! Stop program if catch abnormally large error
        if (energyError > 1.0d-5) then
            print *, "-------------------------WARNING------------------------"
            print *, "Energy error is:", energyError
            print*, "Total energy not conserved over time step in sub-step procedure!"
        end if
      
        if (gaussError > 1.0d-5) then
            print *, "-------------------------WARNING------------------------"
            print *, "Gauss error is:", gaussError
            ! print *, "Re-solve Gauss law"
            ! call solver%solve_tridiag_Poisson(world)
        end if

        Etotal = solver%getTotalPE(world, .false.)
        chargeTotal = 0.0d0
        energyLoss = 0.0d0
        momentum_total = 0.0d0
        do j=1, numberChargedParticles
            call particleList(j)%writeLocalTemperature(CurrentDiagStep, directoryName, MIN(NumberXNodes, NumberXHalfNodes))
            call particleList(j)%writePhaseSpace(directoryName)
            momentum_total = momentum_total + particleList(j)%getTotalMomentum()
            momentum_total(1) = momentum_total(1) + SUM(particleList(j)%momentumLoss(:,:)) * particleList(j)%mass * particleList(j)%w_p
            chargeTotal = chargeTotal + SUM(particleList(j)%accumWallLoss) * particleList(j)%q_times_wp
            Etotal = Etotal + particleList(j)%getTotalKE()
            energyLoss = energyLoss + SUM(particleList(j)%accumEnergyLoss) * particleList(j)%w_p * particleList(j)%mass * 0.5d0
            write(unitPart1+j,"(5(es16.8,1x), (I10, 1x), 3(es16.8,1x))") currentTime + currDel_t, &
                particleList(j)%accumWallLoss(1) * particleList(j)%q * particleList(j)%w_p/(currentTime + currDel_t - pastDiagTime), particleList(j)%accumWallLoss(2) * particleList(j)%q * particleList(j)%w_p/(currentTime + currDel_t - pastDiagTime), &
                particleList(j)%accumEnergyLoss(1)* particleList(j)%mass * particleList(j)%w_p * 0.5d0/(currentTime + currDel_t - pastDiagTime), &
                particleList(j)%accumEnergyLoss(2)* particleList(j)%mass * particleList(j)%w_p * 0.5d0/(currentTime + currDel_t - pastDiagTime), SUM(particleList(j)%N_p), &
                particleList(j)%getKEAve()*2.0d0/3.0d0, particleList(j)%numSubStepsAve, particleList(j)%numFuncEvalAve
            close(unitPart1+j)
        end do
        do j = 1, numberBinaryCollisions
            inelasticEnergyLoss = inelasticEnergyLoss + nullCollisionList(j)%totalEnergyLoss * e * particleList(nullCollisionList(j)%reactantsIndx(1))%w_p
        end do
        write(22,"(9(es16.8,1x), 1(I4, 1x))") currentTime + currDel_t, inelasticEnergyLoss/(currentTime + currDel_t - pastDiagTime), &
        chargeTotal/(currentTime + currDel_t - pastDiagTime), energyLoss/(currentTime + currDel_t - pastDiagTime), &
        Etotal, momentum_total(1), gaussError, chargeError, energyError, iterNumPicard
        close(22)
        call system_clock(endTotal)
        currentTime = currentTime + currDel_t
        elapsed_time = real((endTotal - startTotal), kind = real64) / real(timingRate, kind = real64)
        print *, "Elapsed time for simulation is:", elapsed_time, "seconds"
        print *, "Percentage of steps adaptive is:", 100.0d0 * real(amountTimeSplits)/real(i + 1)
        print *, 'potentialTime is:', real(potentialTime, kind = real64) / real(timingRate, kind = real64)
        print *, 'collision time is:', real(collisionTime, kind=real64)/real(timingRate, kind = real64)

        ! Write Particle properties
        open(9,file=directoryName//'/SimulationFinalData.dat')
        write(9,'("Elapsed Times(s), Potential Time (s), Collision Time (s), Total Steps, Number Adaptive Steps")')
        write(9,"(3(es16.8,1x), 2(I10, 1x))") elapsed_time, real(potentialTime, kind = real64) / real(timingRate, kind = real64), real(collisionTime, kind = real64) / real(timingRate, kind = real64), i+1, amountTimeSplits
        close(9)

    end subroutine solveSimulation


    subroutine solveSimulationFinalAverage(solver, particleList, targetParticleList, nullCollisionList, world, del_t, irand, averagingTime, binNumber)
        ! Perform certain amount of timesteps, with diagnostics taken at first and last time step
        ! Impliment averaging for final select amount of timeSteps, this will be last data dump
        type(Particle), intent(in out) :: particleList(numberChargedParticles)
        type(potentialSolver), intent(in out) :: solver
        type(targetParticle), intent(in) :: targetParticleList(numberNeutralParticles)
        type(nullCollision), intent(in out) :: nullCollisionList(numberBinaryCollisions)
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t, averagingTime
        integer(int32), intent(in) :: binNumber
        integer(int32), intent(in out) :: irand(numThread)
        integer(int32) :: i, j, windowNum, VHist(2*binNumber), intPartV, k, iThread
        real(real64) :: startTime, phi_average(NumberXNodes), currDel_t, remainDel_t
        real(real64) :: chargeLossTotal, ELossTotal, lastCheckTime, checkTimeDivision, meanLoss, stdLoss, RF_ave
        real(real64) :: E_max, VMax
        real(real64), allocatable :: wallLoss(:)
        
        !Save initial particle/field data, along with domain
        do i = 1, numberChargedParticles
            particleList(i)%accumEnergyLoss = 0.0d0
            particleList(i)%accumWallLoss = 0
            particleList(i)%densities = 0
        end do
        inelasticEnergyLoss = 0.0d0
        do j = 1, numberBinaryCollisions
            nullCollisionList(j)%totalEnergyLoss = 0
            nullCollisionList(j)%totalAmountCollisions = 0
        end do
        phi_average = 0.0d0
        i = 0
        energyAddColl = 0.0d0
        startTime = currentTime
        currDel_t = del_t
        remainDel_t = del_t
        lastCheckTime = startTime
        checkTimeDivision = 200.0d0 * del_t/fractionFreq
        windowNum = 0
        RF_ave = 0
        allocate(wallLoss(2 * INT(checkTimeDivision/del_t)))
        do while((currentTime - startTime) < averagingTime)
            call solvePotential(solver, particleList, world, del_t, remainDel_t, currDel_t, currentTime)
            !call ionizationCollisionIsotropic(particleList(1), particleList(2), 1.0d20, 1.0d-20, currDel_t, 15.8d0, 0.0d0, irand)
            if (heatingBool) call maxwellianHeating(particleList(1), irand, fractionFreq, T_e, currDel_t, del_t)
            if (addLostPartBool) call addMaxwellianLostParticles(particleList, T_e, T_i, irand, world)
            if (refluxPartBool) call refluxParticles(particleList, T_e, T_i, irand, world)
            if (injectionBool) call injectAtBoundary(particleList, T_e, T_i, irand, world, currDel_t, solver%BFieldAngle)
            if (uniformInjectionBool) call injectUniformFlux(particleList, T_e, T_i, irand, world)
            do j = 1, numberBinaryCollisions
                call nullCollisionList(j)%generateCollision(particleList, targetParticleList, numberChargedParticles, numberBinaryCollisions, irand, currDel_t)
            end do
            call loadParticleDensity(particleList, world, .false.)
            call solver%aveRFVoltage(.true., phi_average, RF_ave, i, world)
            ! if (MODULO(i+1, heatSkipSteps) == 0) then
            !     call addUniformPowerMaxwellian(particleList(1), Power, nu_h, irand, heatSkipSteps*del_t)
            ! end if
            currentTime = currentTime + currDel_t
            i = i + 1
            windowNum = windowNum + 1
            wallLoss(windowNum) = 0.0d0
            do j = 1, numberChargedParticles
                wallLoss(windowNum) = wallLoss(windowNum) + SUM(particleList(j)%energyLoss)
            end do
            wallLoss(windowNum) = wallLoss(windowNum)/(currentTime - startTime)
            if ((currentTime - lastCheckTime) > checkTimeDivision) then
                meanLoss = SUM(wallLoss(1:windowNum))/real(windowNum)
                stdLoss = SQRT(SUM( (wallLoss(1:windowNum) - meanLoss)**2 )/real(windowNum))
                if (stdLoss/meanLoss < 1d-6) exit
                windowNum = 0
                lastCheckTime = currentTime
            end if
        end do
        deallocate(wallLoss)
        print *, "Averaging finished over", (currentTime - startTime), 'simulation time (s)'
        do j=1, numberChargedParticles
            !$OMP parallel private(iThread)
            iThread = omp_get_thread_num() + 1
            particleList(j)%densities(:, iThread) = particleList(j)%densities(:, iThread) /real(i)
            !$OMP end parallel
        end do
        call solver%aveRFVoltage(.false., phi_average, RF_ave, i, world)
        call writePhi(phi_average, 0, .true., directoryName)
        chargeLossTotal = 0.0d0
        ELossTotal = 0.0d0
        do j = 1, numberChargedParticles
            chargeLossTotal = chargeLossTotal + SUM(particleList(j)%accumWallLoss) * particleList(j)%q_times_wp
            ELossTotal = ELossTotal + SUM(particleList(j)%accumEnergyLoss) * particleList(j)%mass * particleList(j)%w_p * 0.5d0
        end do
        do j = 1, numberBinaryCollisions
            inelasticEnergyLoss = inelasticEnergyLoss + nullCollisionList(j)%totalEnergyLoss * e * particleList(nullCollisionList(j)%reactantsIndx(1))%w_p
        end do
        call writeParticleDensity(particleList, world, 0, .true., directoryName) 
        solver%rho = 0.0d0
        do j=1, numberChargedParticles
            solver%rho = solver%rho + particleList(j)%densities(:,1) * particleList(j)%q
        end do
        gaussError = solver%getError_tridiag_Poisson(world)
        print *, 'gaussError average is:', gaussError
        open(22,file=directoryName//'/GlobalDiagnosticDataAveraged.dat')
        write(22,'("Steps Averaged, Averaging Time, Collision Loss (W/m^2), ParticleCurrentLoss (A/m^2), ParticlePowerLoss(W/m^2), gaussError")')
        write(22,"((I10, 1x), 5(es16.8,1x))") i, (currentTime - startTime), inelasticEnergyLoss/(currentTime-startTime), chargeLossTotal/(currentTime-startTime), ELossTotal/(currentTime-startTime), gaussError
        close(22)
        print *, 'Power loss to walls is:', ELossTotal/(currentTime - startTime)
        print *, 'Power gain in plasma is:', SUM(energyAddColl)/(currentTime - startTime)
        print *, 'Electron total power loss in W/m^2:', SUM(particleList(1)%accumEnergyLoss) * particleList(1)%mass * particleList(1)%w_p * 0.5d0 / (currentTime - startTime)
        print *, 'Ion total power loss in W/m^2:', SUM(particleList(2)%accumEnergyLoss) * particleList(2)%mass * particleList(2)%w_p * 0.5d0 / (currentTime - startTime)
        print *, "Electron average wall loss flux:", SUM(particleList(1)%accumWallLoss)* particleList(1)%w_p/(currentTime - startTime)
        print *, "Ion average wall loss flux:", SUM(particleList(2)%accumWallLoss)* particleList(2)%w_p/(currentTime - startTime)
        print *, "Performing average for EEDF over 50/omega_p"
        E_max = 3.0d0 * (MAXVAL(phi_average) - minval(phi_average))
        print *, 'E_max is:', E_max
        VMax = SQRT(2.0d0 * E_max *e/ m_e)
        print *, "V_max is:", VMax
        checkTimeDivision = 50.0d0 * del_t/fractionFreq
        startTime = currentTime
        VHist = 0.0d0
        do while((currentTime - startTime) < checkTimeDivision)
            call solvePotential(solver, particleList, world, del_t, remainDel_t, currDel_t, currentTime)
            !call ionizationCollisionIsotropic(particleList(1), particleList(2), 1.0d20, 1.0d-20, currDel_t, 15.8d0, 0.0d0, irand)
            if (heatingBool) call maxwellianHeating(particleList(1), irand, fractionFreq, T_e, currDel_t, del_t)
            if (addLostPartBool) call addMaxwellianLostParticles(particleList, T_e, T_i, irand, world)
            if (refluxPartBool) call refluxParticles(particleList, T_e, T_i, irand, world)
            if (injectionBool) call injectAtBoundary(particleList, T_e, T_i, irand, world, currDel_t, solver%BFieldAngle)
            if (uniformInjectionBool) call injectUniformFlux(particleList, T_e, T_i, irand, world)
            do j = 1, numberBinaryCollisions
                call nullCollisionList(j)%generateCollision(particleList, targetParticleList, numberChargedParticles, numberBinaryCollisions, irand, currDel_t)
            end do
            do k = 1, numThread
                do i=1, particleList(1)%N_p(k)
                    intPartV = INT(particleList(1)%phaseSpace(2, i, k) * (binNumber) / VMax + binNumber + 1)
                    if (intPartV > 0 .and. intPartV < binNumber*2+1) VHist(intPartV) = VHist(intPartV) + 1
                end do
            end do
            currentTime = currentTime + currDel_t
        end do
        open(10,file=directoryName//'/ElectronTemperature/EVDF_average.dat', form='UNFORMATTED')
        write(10) real(VHist, kind = real64), VMax
        close(10)





    end subroutine solveSimulationFinalAverage




end module mod_simulation