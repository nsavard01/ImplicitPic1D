module mod_simulation
    ! module which runs actual simulation over several time steps
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
    use iso_c_binding
    use ifport, only: makedirqq
    implicit none

    integer(int32) :: numTimeSteps, heatSkipSteps
    real(real64), private :: energyError, chargeError, gaussError, Power, nu_h, inelasticEnergyLoss
    real(real64) :: currentTime


contains

    ! ------------------------- Reading Input data --------------------------------
    
    subroutine writePhi(phi, CurrentDiagStep, boolAverage, dirName) 
        ! Save phi at single time step or averaged over several time steps
        real(real64), intent(in) :: phi(NumberXNodes)
        integer(int32), intent(in) :: CurrentDiagStep
        character(*), intent(in) :: dirName
        character(len=5) :: char_i
        logical, intent(in) :: boolAverage
        write(char_i, '(I4)') CurrentDiagStep
        if (boolAverage) then
            open(41,file=dirName//'/Phi/phi_Average.dat', form='UNFORMATTED')
        else
            open(41,file=dirName//'/Phi/phi_'//trim(adjustl(char_i))//".dat", form='UNFORMATTED')
        end if
        write(41) phi
        close(41)
        
    end subroutine writePhi
    

    subroutine solveSimulation(solver, particleList, targetParticleList, nullCollisionList, world, del_t, irand, simulationTime)
        ! Perform simulation
        type(Particle), intent(in out) :: particleList(numberChargedParticles)
        type(potentialSolver), intent(in out) :: solver
        type(targetParticle), intent(in) :: targetParticleList(numberNeutralParticles)
        type(nullCollision), intent(in out) :: nullCollisionList(numberBinaryCollisions)
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t, simulationTime
        integer(c_int64_t), intent(in out) :: irand(numThread)
        integer(int32) :: i, j, CurrentDiagStep, smoothInt
        integer(int64) :: startTimer, endTimer, startTotal, endTotal, timingRate
        real(real64) :: diagTimeDivision, diagTime, Etotal, chargeTotal, elapsed_time, pastDiagTime, energyLoss, totCollTime, totPotTime, totSolverTime, totMoverTime
        real(real64) :: currDel_t, remainDel_t
        real(real64) :: KE_i, KE_f, PE_i, PE_f, momentum_total(3)
        integer(int64) :: potentialTime, collisionTime, unitPart1
        character(len=8) :: date_char
        character(len=10) :: time_char
        real(real64) :: rho_i(NumberXNodes)
        allocate(energyAddColl(numThread))
        CurrentDiagStep = oldDiagStep
        unitPart1 = 100
        if (.not. restartBool) call generateSaveDirectory(directoryName, numberBinaryCollisions)
        !Wrtie Initial conditions
        smoothInt = 0
        if (world%gridSmoothBool) smoothInt = 1
        ! Date-time
        call date_and_time(DATE=date_char, TIME = time_char)
        open(15,file=directoryName//'/DateTime.dat')
        write(15,'("UTC Date, UTC Time")')
        write(15,"(2(A,1x))") date_char, time_char 
        close(15)
        ! Save Initial conditions
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

        ! Initialize variables to 0
        collisionTime = 0
        potentialTime = 0
        momentum_total = 0.0d0
        do i = 1, numberChargedParticles
            particleList(i)%accumEnergyLoss = 0.0d0
            particleList(i)%accumWallLoss = 0
            particleList(i)%momentumLoss = 0.0d0
            momentum_total = momentum_total + particleList(i)%getTotalMomentum()
            call particleList(i)%writeLocalTemperature(0, directoryName, MIN(NumberXNodes, NumberXHalfNodes))
            open(unitPart1+i,file=directoryName//'/ParticleDiagnostic_'//particleList(i)%name//'.dat', access = 'APPEND')
            if (.not. restartBool) write(unitPart1+i,'("Time (s), leftCurrentLoss(A/m^2), rightCurrentLoss(A/m^2), leftPowerLoss(W/m^2), rightPowerLoss(W/m^2), Number Particles, Temperature (eV), Ave. Num. SubSteps., Ave. Num. Func. Evals.")')
        end do
        do i = 1, numberBinaryCollisions
            if (.not. restartBool) call nullCollisionList(i)%writeCollisionProperties(particleList, targetParticleList, directoryName)
        end do

        if (restartBool) then
            ! Continue timing from previous simulation
            open(303,file=directoryName//'/SimulationTimeData.dat')
            read(303, *, IOSTAT = j)
            do while (j == 0)
                read(303, *, IOSTAT = j) elapsed_time, totSolverTime, totPotTime, totMoverTime, totCollTime, i
            end do
            close(303)
            open(303,file=directoryName//'/SimulationTimeData.dat', access = 'APPEND')
        else
            elapsed_time = 0
            totPotTime = 0
            totSolverTime = 0
            totMoverTime = 0
            totCollTime = 0
            i = 0
            open(303,file=directoryName//'/SimulationTimeData.dat')  
            write(303,'("Elapsed Times(s), Total Solver Time (s), Potential Time (s), Mover Time (s), Collision Time (s), Total Steps")')
        end if

        inelasticEnergyLoss = 0.0d0
        energyAddColl = 0.0d0
        ! Diagnostic times
        diagTimeDivision = (simulationTime - startSimulationTime)/real(numDiagnosticSteps)
        diagTime = diagTimeDivision + startSimulationTime
        ! Initialize file storing global data
        open(202,file=directoryName//'/GlobalDiagnosticData.dat', access = 'APPEND')
        if (.not. restartBool) write(202,'("Time (s), Collision Loss (W/m^2), ParticleCurrentLoss (A/m^2), ParticlePowerLoss(W/m^2), EnergyTotal (J/m^2), MomentumTotal (kg/s/m), gaussError (a.u), chargeError (a.u), energyError(a.u), Picard Iteration Number")')
        
        !Save initial particle/field data, along with domain
        call loadParticleDensity(particleList, world, .true.)
        call writeParticleDensity(particleList, world, 0, .false., directoryName) 
        call writePhi(solver%phi, 0, .false., directoryName)
        call world%writeDomain(directoryName)
        currentTime = startSimulationTime
        pastDiagTime = startSimulationTime
        currDel_t = del_t
        remainDel_t = del_t
        call system_clock(startTotal)
        do while(currentTime < simulationTime)
            if (currentTime + remainDel_t < diagTime) then
                ! Non-diagnostic operations
                call system_clock(startTimer)
                call solvePotential(solver, particleList, world, del_t, remainDel_t, currDel_t, currentTime)
                call system_clock(endTimer)
                potentialTime = potentialTime + (endTimer - startTimer)
                call system_clock(startTimer)
                if (heatingBool) call maxwellianHeating(particleList(1), irand, fractionFreq, T_e, currDel_t, del_t)
                if (addLostPartBool) call addMaxwellianLostParticles(particleList, T_e, T_i, irand, world)
                if (refluxPartBool) call refluxParticles(particleList, T_e, T_i, irand, world)
                if (injectionBool) call injectAtBoundary(particleList, T_e, T_i, irand, world, currDel_t)
                if (uniformInjectionBool) call injectUniformFlux(particleList, T_e, T_i, irand, world)
                do j = 1, numberBinaryCollisions
                    call nullCollisionList(j)%generateCollision(particleList, targetParticleList, numberChargedParticles, numberBinaryCollisions, irand, currDel_t)
                end do
                call system_clock(endTimer)
                collisionTime = collisionTime + (endTimer - startTimer)
                
            else  
                ! Data dump with diagnostics
                print *, "Simulation is", (currentTime-startSimulationTime)/(simulationTime - startSimulationTime) * 100.0, "percent done"
                ! Get kinetic energy total and charge density before solving potential
                KE_i = 0.0d0
                do j=1, numberChargedParticles
                    KE_i = KE_i + particleList(j)%getTotalKE()
                end do
                call solver%depositRho(particleList, world)
                rho_i = solver%rho
                call system_clock(startTimer)
                call solvePotential(solver, particleList, world, del_t, remainDel_t, currDel_t, currentTime)
                call system_clock(endTimer)
                ! Calculate error in energy
                PE_i = solver%getTotalPE(world, .false.)
                KE_f = 0.0d0
                do j=1, numberChargedParticles
                    KE_f = KE_f + particleList(j)%getTotalKE() + SUM(particleList(j)%energyLoss) * particleList(j)%mass * particleList(j)%w_p * 0.5d0
                end do
                PE_f = solver%getTotalPE(world, .true.) - solver%getEnergyFromBoundary(world, currDel_t) 
                energyError = ABS((KE_i + PE_i - KE_f - PE_f)/(KE_i + PE_i))
                potentialTime = potentialTime + (endTimer - startTimer)
                
                
                call system_clock(startTimer)
                if (heatingBool) call maxwellianHeating(particleList(1), irand, fractionFreq, T_e, currDel_t, del_t)
                if (addLostPartBool) call addMaxwellianLostParticles(particleList, T_e, T_i, irand, world)
                if (refluxPartBool) call refluxParticles(particleList, T_e, T_i, irand, world)
                if (injectionBool) call injectAtBoundary(particleList, T_e, T_i, irand, world, currDel_t)
                if (uniformInjectionBool) call injectUniformFlux(particleList, T_e, T_i, irand, world)
                do j = 1, numberBinaryCollisions
                    call nullCollisionList(j)%generateCollision(particleList, targetParticleList, numberChargedParticles, numberBinaryCollisions, irand, currDel_t)
                end do
                call system_clock(endTimer)
                collisionTime = collisionTime + (endTimer - startTimer)
                ! Save current diagnostics
                call loadParticleDensity(particleList, world, .true.)
                call writeParticleDensity(particleList, world, CurrentDiagStep, .false., directoryName) 
                call writePhi(solver%phi_f, CurrentDiagStep, .false., directoryName)
                call solver%depositRho(particleList, world)

                !charge conservation directly
                chargeError = solver%getChargeContinuityError(rho_i, world, currDel_t)

                ! Get error gauss' law
                gaussError = solver%getError_tridiag_Poisson(world)
                
                ! Print if errors are above 1e-5
                if (energyError > 1.0d-5) then
                    print *, "-------------------------WARNING------------------------"
                    print *, "Energy error is:", energyError
                end if
                
                if (gaussError > 1.0d-5) then
                    print *, "-------------------------WARNING------------------------"
                    print *, "Gauss error is:", gaussError
                end if
                
                Etotal = solver%getTotalPE(world, .true.)
                chargeTotal = 0.0d0
                energyLoss = 0.0d0
                momentum_total = 0.0d0
                do j=1, numberChargedParticles
                    ! Save particle-dependent diagnostics
                    call particleList(j)%writeLocalTemperature(CurrentDiagStep, directoryName, MIN(NumberXNodes, NumberXHalfNodes))
                    ! Save phase space once per diagnostic, overwrite
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
                    ! Save diagnostics for each binary collision
                    inelasticEnergyLoss = inelasticEnergyLoss + 0.5d0 * SUM(nullCollisionList(j)%totalEnergyLoss) * particleList(nullCollisionList(j)%reactantsIndx(1))%w_p
                    call nullCollisionList(j)%writeCollisionDiag(particleList, targetParticleList, directoryName, (currentTime + currDel_t - pastDiagTime))
                end do
                ! Save global diagnostics
                write(202,"(9(es16.8,1x), 1(I4, 1x))") currentTime + currDel_t, inelasticEnergyLoss/(currentTime + currDel_t - pastDiagTime), &
                chargeTotal/(currentTime + currDel_t - pastDiagTime), energyLoss/(currentTime + currDel_t - pastDiagTime), Etotal, momentum_total(1), &
                gaussError, chargeError, energyError, iterNumPicard
                call system_clock(endTotal)
                ! Save accumulated time diagnostics
                write(303,"(5(es16.8,1x), 1(I8, 1x))") elapsed_time + real((endTotal - startTotal), kind = real64) / real(timingRate, kind = real64), &
                totSolverTime + real(potentialTime, kind = real64) / real(timingRate, kind = real64), &
                totPotTime + real(potTimer, kind = real64) / real(timingRate, kind = real64), &
                totMoverTime + real(moverTimer, kind = real64) / real(timingRate, kind = real64), &
                totCollTime + real(collisionTime, kind = real64) / real(timingRate, kind = real64), i+1
                ! Reset accumulated diagnostic variables for collisions
                do j = 1, numberBinaryCollisions
                    nullCollisionList(j)%totalNumberCollidableParticles = 0
                    nullCollisionList(j)%totalEnergyLoss = 0
                    nullCollisionList(j)%totalAmountCollisions = 0
                    nullCollisionList(j)%totalIncidentEnergy = 0
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
            currentTime = currentTime + currDel_t
            i = i + 1
        end do
        
        close(202)
        close(303)
        call system_clock(endTotal)
        ! Calculate and save final timings
        elapsed_time = elapsed_time + real((endTotal - startTotal), kind = real64) / real(timingRate, kind = real64)
        totSolverTime = totSolverTime + real(potentialTime, kind = real64) / real(timingRate, kind = real64)
        totPotTime = totPotTime + real(potTimer, kind = real64) / real(timingRate, kind = real64)
        totMoverTime = totMoverTime + real(moverTimer, kind = real64) / real(timingRate, kind = real64)
        totCollTime = totCollTime + real(collisionTime, kind = real64) / real(timingRate, kind = real64)
        currentTime = currentTime + currDel_t
        print *, "Elapsed time for simulation is:", elapsed_time, "seconds"
        print *, "Percentage of steps adaptive is:", 100.0d0 * real(amountTimeSplits)/real(i + 1)
        print *, 'total solver time is:', totSolverTime
        print *, 'total collision time is:', totCollTime

        ! Write Particle properties
        open(9,file=directoryName//'/SimulationFinalData.dat')
        write(9,'("Elapsed Times(s), Solver Time (s), Collision Time (s), Total Steps, Number Adaptive Steps")')
        write(9,"(3(es16.8,1x), 2(I10, 1x))") elapsed_time, totSolverTime, totCollTime, i+1, amountTimeSplits
        close(9)

    end subroutine solveSimulation


    subroutine solveSimulationFinalAverage(solver, particleList, targetParticleList, nullCollisionList, world, del_t, irand, averagingTime, binNumber)
        ! Averaging over a certain amount of simulated time
        type(Particle), intent(in out) :: particleList(numberChargedParticles)
        type(potentialSolver), intent(in out) :: solver
        type(targetParticle), intent(in) :: targetParticleList(numberNeutralParticles)
        type(nullCollision), intent(in out) :: nullCollisionList(numberBinaryCollisions)
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t, averagingTime
        integer(int32), intent(in) :: binNumber
        integer(c_int64_t), intent(in out) :: irand(numThread)
        integer(int32) :: i, j, intPartV, k, iThread
        real(real64) :: startTime, phi_average(NumberXNodes), currDel_t, remainDel_t
        real(real64) :: chargeLossTotal, ELossTotal, lastCheckTime, checkTimeDivision, meanLoss, stdLoss, RF_ave, partV, partE
        real(real64) :: VMax(numberChargedParticles), EMax(numberChargedParticles), Emin(numberChargedParticles), &
        E_grid(binNumber, numberChargedParticles), diffE(numberChargedParticles), Emin_log(numberChargedParticles), VHist(2*binNumber, numberChargedParticles), &
        EHist(binNumber, numberChargedParticles)
        
        ! Initialize accumulation data
        do i = 1, numberChargedParticles
            particleList(i)%accumEnergyLoss = 0.0d0
            particleList(i)%accumWallLoss = 0
            particleList(i)%densities = 0
        end do
        inelasticEnergyLoss = 0.0d0
        do j = 1, numberBinaryCollisions
            nullCollisionList(j)%totalEnergyLoss = 0
            nullCollisionList(j)%totalAmountCollisions = 0
            nullCollisionList(j)%totalIncidentEnergy = 0
            nullCollisionList(j)%totalNumberCollidableParticles = 0
        end do
        phi_average = 0.0d0
        i = 0
        energyAddColl = 0.0d0
        startTime = currentTime
        currDel_t = del_t
        remainDel_t = del_t
        lastCheckTime = startTime
        checkTimeDivision = 200.0d0 * del_t/fractionFreq
        RF_ave = 0

        do while((currentTime - startTime) < averagingTime)
            call solvePotential(solver, particleList, world, del_t, remainDel_t, currDel_t, currentTime)
            
            if (heatingBool) call maxwellianHeating(particleList(1), irand, fractionFreq, T_e, currDel_t, del_t)
            if (addLostPartBool) call addMaxwellianLostParticles(particleList, T_e, T_i, irand, world)
            if (refluxPartBool) call refluxParticles(particleList, T_e, T_i, irand, world)
            if (injectionBool) call injectAtBoundary(particleList, T_e, T_i, irand, world, currDel_t)
            if (uniformInjectionBool) call injectUniformFlux(particleList, T_e, T_i, irand, world)
            do j = 1, numberBinaryCollisions
                call nullCollisionList(j)%generateCollision(particleList, targetParticleList, numberChargedParticles, numberBinaryCollisions, irand, currDel_t)
            end do
            if (remainDel_t == del_t) then
                call loadParticleDensity(particleList, world, .false.)
                call solver%aveRFVoltage(.true., phi_average, RF_ave, i, world)
                i = i + 1
            end if
            currentTime = currentTime + currDel_t
            
        end do
        print *, "Averaging finished over", (currentTime - startTime), 'simulation time (s)'
        ! Average densities
        do j=1, numberChargedParticles
            !$OMP parallel private(iThread)
            iThread = omp_get_thread_num() + 1
            particleList(j)%densities(:, iThread) = particleList(j)%densities(:, iThread) /real(i, kind = 8)
            !$OMP end parallel
            open(22,file=directoryName//'/ParticleAveDiagnostic_'//particleList(j)%name//'.dat')
            write(22, '("Left curr (A/m^2), right curr (A/m^2), left power (W/m^2), right power (W/m^2)")')
            write(22,"(4(es16.8, 1x))") particleList(j)%accumWallLoss(1) * particleList(j)%q_times_wp /(currentTime-startTime), &
                particleList(j)%accumWallLoss(2) * particleList(j)%q_times_wp /(currentTime-startTime), &
                particleList(j)%accumEnergyLoss(1)* particleList(j)%mass * particleList(j)%w_p * 0.5d0/(currentTime-startTime), &
                particleList(j)%accumEnergyLoss(2)* particleList(j)%mass * particleList(j)%w_p * 0.5d0/(currentTime-startTime)
            close(22)
        end do
        ! Ave voltage (including RF)
        call solver%aveRFVoltage(.false., phi_average, RF_ave, i, world)
        call writePhi(phi_average, 0, .true., directoryName)
        chargeLossTotal = 0.0d0
        ELossTotal = 0.0d0
        do j = 1, numberChargedParticles
            chargeLossTotal = chargeLossTotal + SUM(particleList(j)%accumWallLoss) * particleList(j)%q_times_wp
            ELossTotal = ELossTotal + SUM(particleList(j)%accumEnergyLoss) * particleList(j)%mass * particleList(j)%w_p * 0.5d0
        end do
        ! Averaged binary Collision data dump
        do j = 1, numberBinaryCollisions
            inelasticEnergyLoss = inelasticEnergyLoss + SUM(nullCollisionList(j)%totalEnergyLoss) * 0.5d0 * particleList(nullCollisionList(j)%reactantsIndx(1))%w_p
            open(202,file=directoryName//'/BinaryCollisions/'//particleList(nullCollisionList(j)%reactantsIndx(1))%name//"_on_"//targetParticleList(nullCollisionList(j)%reactantsIndx(2))%name//"/AveCollisionDiag.dat")
            write(202, '("Coll #, CollRatio, AveEnergyLoss (eV), AveIncidentEnergy (eV), P_loss(W/m^2), aveCollFreq (Hz/m^2)")')
            do k = 1, nullCollisionList(j)%numberCollisions
                write(202,"((I3, 1x), 5(es16.8,1x))") k, real(nullCollisionList(j)%totalAmountCollisions(k))/real(nullCollisionList(j)%totalNumberCollidableParticles), &
                nullCollisionList(j)%totalEnergyLoss(k) * 0.5d0 / e/ real(nullCollisionList(j)%totalAmountCollisions(k)), &
                nullCollisionList(j)%totalIncidentEnergy(k) * particleList(nullCollisionList(j)%reactantsIndx(1))%mass * 0.5d0 / e/ real(nullCollisionList(j)%totalAmountCollisions(k)), &
                nullCollisionList(j)%totalEnergyLoss(k) * 0.5d0 * particleList(nullCollisionList(j)%reactantsIndx(1))%w_p / (currentTime-startTime), &
                real(nullCollisionList(j)%totalAmountCollisions(k)) * particleList(nullCollisionList(j)%reactantsIndx(1))%w_p / (currentTime-startTime)
            end do
            close(202)
        end do
        call writeParticleDensity(particleList, world, 0, .true., directoryName) 
        solver%rho = 0.0d0
        do j=1, numberChargedParticles
            solver%rho = solver%rho + particleList(j)%densities(:,1) * particleList(j)%q
        end do
        ! averaged error in gauss from averaged phi and rho
        gaussError = solver%getError_tridiag_Poisson(world)
        print *, 'gaussError average is:', gaussError
        open(202,file=directoryName//'/GlobalDiagnosticDataAveraged.dat')
        write(202,'("Steps Averaged, Averaging Time, Collision Loss (W/m^2), ParticleCurrentLoss (A/m^2), ParticlePowerLoss(W/m^2), gaussError")')
        write(202,"((I10, 1x), 5(es16.8,1x))") i, (currentTime - startTime), inelasticEnergyLoss/(currentTime-startTime), chargeLossTotal/(currentTime-startTime), ELossTotal/(currentTime-startTime), gaussError
        close(202)
        print *, 'Power loss to walls is:', ELossTotal/(currentTime - startTime)
        print *, 'Power gain in plasma is:', SUM(energyAddColl)/(currentTime - startTime)
        print *, 'Electron total power loss in W/m^2:', SUM(particleList(1)%accumEnergyLoss) * particleList(1)%mass * particleList(1)%w_p * 0.5d0 / (currentTime - startTime)
        print *, 'Ion total power loss in W/m^2:', SUM(particleList(2)%accumEnergyLoss) * particleList(2)%mass * particleList(2)%w_p * 0.5d0 / (currentTime - startTime)
        print *, "Electron average wall loss flux:", SUM(particleList(1)%accumWallLoss)* particleList(1)%w_p/(currentTime - startTime)
        print *, "Ion average wall loss flux:", SUM(particleList(2)%accumWallLoss)* particleList(2)%w_p/(currentTime - startTime)
        print *, "Performing average for EEDF over 50/omega_p"
        VMax = 0
        EMax = 0
        do j = 1, numberChargedParticles
            Emin(j) = 0.1d0 * e * 2.0d0 / particleList(j)%mass
        end do
        ! Create histogram array for velocity and energy of each particle
        !$OMP parallel private(iThread, j, i, partV) reduction(max:VMax, EMax) reduction(min:Emin)
        iThread = omp_get_thread_num() + 1
        do j = 1, numberChargedParticles
            do i = 1, particleList(j)%N_p(iThread)
                partV = ABS(particleList(j)%phaseSpace(2,i, iThread))
                VMax(j) = MAX(VMax(j), partV)
                partV = SUM(particleList(j)%phaseSpace(2:4,i, iThread)**2)
                EMax(j) = MAX(EMax(j), partV)
                Emin(j) = MIN(Emin(j), partV)
            end do
        end do
        !$OMP end parallel
        ! Use log mapping for energy array
        Emin_log = LOG(Emin)
        do j = 1, numberChargedParticles
            diffE(j) = (LOG(EMax(j)) - LOG(Emin(j)))/real(binNumber - 1)
            partV = LOG(Emin(j))
            do i = 1, binNumber
                E_grid(i, j) = EXP(partV)
                partV = partV + diffE(j)
            end do
        end do
        ! Either over 50/omega_pe or over an RF period for averaging distribution function
        if (solver%RF_bool) then
            checkTimeDivision = 2.0d0 * pi / solver%RF_rad_frequency
        else if (EField_heating_bool) then
            checkTimeDivision = 1.0d0 / Efield_heating_freq
        else
            checkTimeDivision = 50.0d0 * del_t/fractionFreq
        end if
        startTime = currentTime
        VHist = 0.0d0
        EHist = 0.0d0
        do while((currentTime - startTime) < checkTimeDivision)
            call solvePotential(solver, particleList, world, del_t, remainDel_t, currDel_t, currentTime)
            !call ionizationCollisionIsotropic(particleList(1), particleList(2), 1.0d20, 1.0d-20, currDel_t, 15.8d0, 0.0d0, irand)
            if (heatingBool) call maxwellianHeating(particleList(1), irand, fractionFreq, T_e, currDel_t, del_t)
            if (addLostPartBool) call addMaxwellianLostParticles(particleList, T_e, T_i, irand, world)
            if (refluxPartBool) call refluxParticles(particleList, T_e, T_i, irand, world)
            if (injectionBool) call injectAtBoundary(particleList, T_e, T_i, irand, world, currDel_t)
            if (uniformInjectionBool) call injectUniformFlux(particleList, T_e, T_i, irand, world)
            do j = 1, numberBinaryCollisions
                call nullCollisionList(j)%generateCollision(particleList, targetParticleList, numberChargedParticles, numberBinaryCollisions, irand, currDel_t)
            end do
            if (remainDel_t == del_t) then
                !$OMP parallel private(iThread, j, k, partV, intPartV, partE) reduction(+:VHist, EHist)
                iThread = omp_get_thread_num() + 1
                do j = 1, numberChargedParticles
                    do k = 1, particleList(j)%N_p(iThread)
                        partV = particleList(j)%phaseSpace(2, k, iThread) * (2.0d0 * binNumber - 1) / 2.0d0 / VMax(j) + binNumber + 0.5d0
                        if (partV > 1 .and. partV < binNumber*2) then
                            intPartV = INT(partV)
                            VHist(intPartV, j) = VHist(intPartV, j) + (1.0d0 - (partV - intPartV))
                            VHist(intPartV+1, j) = VHist(intPartV + 1, j) + partV - intPartV
                        end if
                        partE = SUM(particleList(j)%phaseSpace(2:4, k, iThread)**2)
                        if (partE > Emin(j) .and. partE < EMax(j)) then
                            partV = (LOG(partE) - Emin_log(j))/diffE(j) + 1.0d0
                            intPartV = INT(partV)
                            partV = partV - intPartV
                            EHist(intPartV, j) = EHist(intPartV, j) + (1.0d0 - partV)
                            EHist(intPartV+1, j) = EHist(intPartV + 1, j) + partV
                        end if
                    end do
                end do
                !$OMP end parallel
            end if
            currentTime = currentTime + currDel_t
        end do
        ! Save particle distribution data
        do j = 1, numberChargedParticles
            open(10,file=directoryName//'/Temperature/Temp_'//particleList(j)%name//'_average.dat', form='UNFORMATTED')
            write(10) real(VHist(:, j), kind = real64), VMax(j)
            close(10)
            open(10,file=directoryName//'/Temperature/TempEnergy_'//particleList(j)%name//'_average.dat', form='UNFORMATTED')
            write(10) real(EHist(:, j), kind = real64), E_grid(:,j) * 0.5d0 * particleList(j)%mass / e
            close(10)
        end do





    end subroutine solveSimulationFinalAverage




end module mod_simulation