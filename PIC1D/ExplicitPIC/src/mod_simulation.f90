module mod_simulation
    ! Module for running simulation over several time steps
    use iso_fortran_env, only: int32, real64, int64
    use constants
    use mod_BasicFunctions
    use mod_particle
    use mod_domain
    use mod_potentialSolver
    use mod_particleInjection
    use mod_targetParticle
    use mod_NullCollision
    use mod_particle_operations
    use iso_c_binding
    use omp_lib
    implicit none

    real(real64) :: inelasticEnergyLoss, currentTime

contains


    ! --------------------------- Diagnostics ------------------------------------

    subroutine writePhi(phi, CurrentDiagStep, boolAverage, dirName) 
        ! Write phi array 
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


    subroutine WriteParticleDensity(particleList, world, CurrentDiagStep, boolAverage, dirName) 
        ! Write particle densities, option for averaged
        type(Particle), intent(in) :: particleList(numberChargedParticles)
        type(Domain), intent(in) :: world
        character(*), intent(in) :: dirName
        integer(int32), intent(in) :: CurrentDiagStep
        logical, intent(in) :: boolAverage
        real(real64) :: densities(NumberXNodes)
        integer(int32) :: i, iThread, leftThreadIndx, rightThreadIndx
        character(len=5) :: char_i
        
        do i=1, numberChargedParticles
            densities = 0.0d0
            !$OMP parallel private(iThread, leftThreadIndx, rightThreadIndx)
            iThread = omp_get_thread_num() + 1
            leftThreadIndx = world%threadNodeIndx(1,iThread)
            rightThreadIndx = world%threadNodeIndx(2,iThread)
            densities(leftThreadIndx:rightThreadIndx) = SUM(particleList(i)%densities(leftThreadIndx:rightThreadIndx, :), DIM=2) * particleList(i)%w_p
            !$OMP end parallel
            densities = densities/world%delX
            if (world%boundaryConditions(1) == 3) then
                densities(1) = (densities(1) + densities(NumberXNodes))
                densities(NumberXNodes) = densities(1)
            else
                densities(1) = 2.0d0 * densities(1)
                densities(NumberXNodes) = 2.0d0 * densities(NumberXNodes)
            end if
            write(char_i, '(I3)'), CurrentDiagStep
            if (boolAverage) then
                open(41,file=dirName//'/Density/density_'//particleList(i)%name//"_Average.dat", form='UNFORMATTED')
            else
                open(41,file=dirName//'/Density/density_'//particleList(i)%name//"_"//trim(adjustl(char_i))//".dat", form='UNFORMATTED')
            end if
            write(41) densities
            close(41)
        end do
        
    end subroutine WriteParticleDensity


    subroutine solveSimulation(solver, particleList, targetParticleList, nullCollisionList, world, del_t, irand, simulationTime)
        ! Run simulation over several time steps with included diagnostic data dumps
        type(Particle), intent(in out) :: particleList(numberChargedParticles)
        type(targetParticle), intent(in) :: targetParticleList(numberNeutralParticles)
        type(nullCollision), intent(in out) :: nullCollisionList(numberBinaryCollisions)
        type(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t, simulationTime
        integer(c_int64_t), intent(in out) :: irand(numThread)
        integer(int32) :: i, j, CurrentDiagStep, diagStepDiff, unitPart1
        real(real64) :: diagTimeDivision, diagTime, elapsedTime, chargeTotal, energyLoss, elapsed_time, momentum_total(3), totMoverTime, totPotTime, totCollTime
        integer(int64) :: startTime, endTime, timingRate, collisionTime, potentialTime, moverTime, startTotal, endTotal, sub_step
        character(len=8) :: date_char
        character(len=10) :: time_char
        allocate(energyAddColl(numThread))
        CurrentDiagStep = oldDiagStep
        unitPart1 = 100
        ! Create save directory
        if (.not. restartBool) call generateSaveDirectory(directoryName, numberBinaryCollisions)
        ! Write initial conditions
        call date_and_time(DATE=date_char, TIME = time_char)
        open(15,file=directoryName//'/DateTime.dat')
        write(15,'("UTC Date, UTC Time")')
        write(15,"(2(A,1x))") date_char, time_char 
        close(15)
        open(15,file=directoryName//'/InitialConditions.dat')
        write(15,'("Number Grid Nodes, Final Expected Time(s), Delta t(s), FractionFreq, delX, &
            n_ave, T_e, T_i, numDiag, NumChargedPart, numThread, RF_frequency, RF_half_amplitude, ion_step")')
        write(15,"((I6, 1x), 7(es16.8,1x), 3(I6, 1x), 2(es16.8,1x), 2(I3, 1x))") NumberXNodes, simulationTime, del_t, FractionFreq, world%delX, n_ave, T_e, T_i, numDiagnosticSteps, &
            numberChargedParticles, numThread, solver%RF_rad_frequency, solver%RF_half_amplitude, ionStepMult, 0
        close(15)
        call system_clock(count_rate = timingRate)
        ! Write Particle properties
        open(9,file=directoryName//'/ParticleProperties.dat')
        write(9,'("Particle Symbol, Particle Mass (kg), Particle Charge (C), Particle Weight (N/m^2), maxIdx")')
        do j=1, numberChargedParticles
            write(9,"((A, 1x), 3(es16.8,1x), (I10, 1x))") particleList(j)%name//'       ', particleList(j)%mass, particleList(j)%q, particleList(j)%w_p, particleList(j)%finalIdx
        end do
        close(9)
        open(9,file=directoryName//'/TargetParticleProperties.dat')
        write(9,'("Particle Symbol, Particle Mass (kg), Particle Temp (K), Particle Density (m^-3)")')
        do j=1, numberNeutralParticles
            write(9,"((A, 1x), 3(es16.8,1x))") targetParticleList(j)%name//'       ', targetParticleList(j)%mass, targetParticleList(j)%temperature, targetParticleList(j)%density
        end do
        close(9)
        ! Initialize particle diagnostics
        do i = 1, numberChargedParticles
            particleList(i)%accumEnergyLoss = 0.0d0
            particleList(i)%accumWallLoss = 0
            particleList(i)%momentumLoss = 0.0d0
            call particleList(i)%writeLocalTemperature(0, directoryName, NumberXHalfNodes)
            open(unitPart1+i,file=directoryName//'/ParticleDiagnostic_'//particleList(i)%name//'.dat', access = 'APPEND')
            if (.not. restartBool) write(unitPart1+i,'("Time (s), leftCurrentLoss(A/m^2), rightCurrentLoss(A/m^2), leftPowerLoss(W/m^2), rightPowerLoss(W/m^2), N_p, Temp")')
        end do
        do i = 1, numberBinaryCollisions
            if (.not. restartBool) call nullCollisionList(i)%writeCollisionProperties(particleList, targetParticleList, directoryName)
        end do

        if (restartBool) then
            ! If restart, continue timing from previous simulation
            open(303,file=directoryName//'/SimulationTimeData.dat')
            read(303, *, IOSTAT = j)
            do while (j == 0)
                read(303, *, IOSTAT = j) elapsed_time, totPotTime, totMoverTime, totCollTime, i
            end do
            close(303)
            open(303,file=directoryName//'/SimulationTimeData.dat', access = 'APPEND')
        else
            elapsed_time = 0
            totPotTime = 0
            totMoverTime = 0
            totCollTime = 0
            i = 0
            open(303,file=directoryName//'/SimulationTimeData.dat')  
            write(303,'("Elapsed Times(s), Potential Time (s), Mover Time (s), Collision Time (s), Total Steps")')
        end if
        collisionTime = 0
        potentialTime = 0
        moverTime = 0
        currentTime = startSimulationTime
        elapsedTime = 0.0d0
        diagStepDiff = 1
        diagTimeDivision = (simulationTime-currentTime)/real(numDiagnosticSteps)
        diagTime = diagTimeDivision + currentTime
        open(202,file=directoryName//'/GlobalDiagnosticData.dat', access = 'APPEND')
        if (.not. restartBool) write(202,'("Time (s), Collision Loss (W/m^2), ParticleCurrentLoss (A/m^2), ParticlePowerLoss(W/m^2), TotalMomentum(kg/m/s), TotalEnergy(J/m^2)")')

        
        
        ! Save initial grid properties
        energyAddColl = 0.0d0
        inelasticEnergyLoss = 0.0d0
        call loadParticleDensity(particleList, world, .true.)
        call writeParticleDensity(particleList, world, 0, .false., directoryName)  
        call writePhi(solver%phi, 0, .false., directoryName)
        call world%writeDomain(directoryName)
        
        call system_clock(startTotal)
        do while(currentTime < simulationTime)
 
            ! move ions forward to ionStepMult * del_t
            call system_clock(startTime)
            if (numberChargedParticles > 1) call moveParticles(particleList(2:numberChargedParticles), solver%EField, world, ionStepMult * del_t)
            call system_clock(endTime)
            moverTime = moverTime + (endTime - startTime)
           
            ! Now go extra (ionStepMult -1)/2 steps electrons
            
            do sub_step = 1, ionStepMult
                ! move electrons
                call system_clock(startTime)
                call moveParticles(particleList(1:1), solver%EField, world, del_t)
                call system_clock(endTime)
                moverTime = moverTime + (endTime - startTime)
                currentTime = currentTime + del_t

                do j = 2, numberChargedParticles
                    ! For ionStepMult-1/2 steps, use ion rho from initial step, then accumulate through collisions
                    ! For center step, will then after accumulate total ion rho from ionStepMult*del_t by resetting startIdx
                    ! Then remaining two will again accumulate, with base rho from ionStepMult*del_t
                    if (sub_step /= (ionStepMult-1)/2 + 1) then
                        particleList(j)%startIdx = particleList(j)%N_p + 1
                    else
                        particleList(j)%startIdx = 1
                    end if
                end do
                
                ! undergo collisions for each time step
                call system_clock(startTime)
                if (particleList(1)%mass == m_e .and. numberChargedParticles > 1) then
                    if (heatingBool) call maxwellianHeating(particleList(1), irand, fractionFreq, T_e, del_t, del_t)
                    if (EField_heating_bool) call Efield_heating(particleList(1), currentTime - del_t, currentTime, world)
                    if (addLostPartBool) call addMaxwellianLostParticles(particleList, T_e, T_i, irand, world)
                    if (refluxPartBool) call refluxParticles(particleList, T_e, T_i, irand, world)
                    if (injectionBool) call injectAtBoundary(particleList, T_e, T_i, irand, world, del_t)
                    if (uniformInjectionBool) call injectUniformFlux(particleList, T_e, T_i, irand, world)
                end if
                do j = 1, numberBinaryCollisions
                    if (nullCollisionList(j)%reactantsIndx(1) == 1) then
                        call nullCollisionList(j)%generateCollision(particleList, targetParticleList, numberChargedParticles, numberBinaryCollisions, irand, del_t)  
                    else if (sub_step == ionStepMult) then
                        call nullCollisionList(j)%generateCollision(particleList, targetParticleList, numberChargedParticles, numberBinaryCollisions, irand, ionStepMult * del_t)
                    end if
                end do
                call system_clock(endTime)
                collisionTime = collisionTime + (endTime - startTime)
        
                ! Re-evaluate rho with added particles
                call system_clock(startTime)
                call depositRho(solver, particleList, world)
                call system_clock(endTime)
                moverTime = moverTime + (endTime - startTime)

                ! re-evaluate potential
                call system_clock(startTime)
                call solver%solve_tridiag_Poisson(world, currentTime)
                call solver%makeEField(world)
                call system_clock(endTime)
                potentialTime = potentialTime + (endTime - startTime)
            end do
        

            if (currentTime >= diagTime) then
                ! Operations with diagnostics
                print *, "Simulation is", currentTime/simulationTime * 100.0, "percent done"

                ! Save grid quantities
                call loadParticleDensity(particleList, world, .true.)
                call writeParticleDensity(particleList, world, CurrentDiagStep, .false., directoryName)  
                call writePhi(solver%phi, CurrentDiagStep, .false., directoryName)
                chargeTotal = 0.0d0
                energyLoss = 0.0d0
                momentum_total = 0.0d0
                do j=1, numberChargedParticles
                    ! Save particle diagnostics
                    call particleList(j)%writeLocalTemperature(CurrentDiagStep, directoryName, NumberXHalfNodes)
                    call particleList(j)%writePhaseSpace(directoryName)
                    chargeTotal = chargeTotal + SUM(particleList(j)%accumWallLoss) * particleList(j)%q * particleList(j)%w_p
                    energyLoss = energyLoss + SUM(particleList(j)%accumEnergyLoss) * particleList(j)%mass * particleList(j)%w_p * 0.5d0
                    momentum_total = momentum_total + particleList(j)%getTotalMomentum()
                    momentum_total(1) = momentum_total(1) + SUM(particleList(j)%momentumLoss(:,:)) * particleList(j)%mass * particleList(j)%w_p
                    write(unitPart1+j,"(5(es16.8,1x), (I10,1x), (es16.8,1x))") currentTime, &
                        particleList(j)%accumWallLoss(1) * particleList(j)%q * particleList(j)%w_p/del_t/diagStepDiff, particleList(j)%accumWallLoss(2) * particleList(j)%q * particleList(j)%w_p/del_t/diagStepDiff, &
                        particleList(j)%accumEnergyLoss(1)* particleList(j)%mass * particleList(j)%w_p * 0.5d0/del_t/diagStepDiff, particleList(j)%accumEnergyLoss(2)* particleList(j)%mass * particleList(j)%w_p * 0.5d0/del_t/diagStepDiff, &
                        SUM(particleList(j)%N_p), particleList(j)%getKEAve()*2.0d0/3.0d0
                    particleList(j)%accumEnergyLoss = 0.0d0
                    particleList(j)%accumWallLoss = 0
                end do
                do j = 1, numberBinaryCollisions
                    ! Save collision diagnostics
                    inelasticEnergyLoss = inelasticEnergyLoss + 0.5d0 * SUM(nullCollisionList(j)%totalEnergyLoss) * particleList(nullCollisionList(j)%reactantsIndx(1))%w_p
                    call nullCollisionList(j)%writeCollisionDiag(particleList, targetParticleList, directoryName, del_t * diagStepDiff)
                end do
                ! write global diagnostics
                write(202,"(6(es16.8,1x))") currentTime, inelasticEnergyLoss/del_t/diagStepDiff, &
                chargeTotal/del_t/diagStepDiff, energyLoss/del_t/diagStepDiff, momentum_total(1), solver%getTotalE(particleList, world)
                call system_clock(endTotal)
                ! Write timing data
                write(303,"(4(es16.8,1x), 1(I8, 1x))") elapsed_time + real((endTotal - startTotal), kind = real64) / real(timingRate, kind = real64), &
                totPotTime + real(potentialTime, kind = real64) / real(timingRate, kind = real64), &
                totMoverTime + real(moverTime, kind = real64) / real(timingRate, kind = real64), &
                totCollTime + real(collisionTime, kind = real64) / real(timingRate, kind = real64), i+1

                do j = 1, numberBinaryCollisions
                    ! reset accumulation variables
                    nullCollisionList(j)%totalNumberCollidableParticles = 0
                    nullCollisionList(j)%totalEnergyLoss = 0
                    nullCollisionList(j)%totalAmountCollisions = 0
                    nullCollisionList(j)%totalIncidentEnergy = 0
                end do
                CurrentDiagStep = CurrentDiagStep + 1
                print *, "Number of electrons is:", SUM(particleList(1)%N_p)
                print *, 'Number of ions is:', SUM(particleList(2)%N_p)
                energyAddColl = 0.0d0
                inelasticEnergyLoss = 0.0d0
                diagStepDiff = 0
                diagTime = diagTime + diagTimeDivision
            end if
           
            i = i + ionStepMult
            diagStepDiff = diagStepDiff + ionStepMult
        end do
        
        close(202)
        close(303)
        call system_clock(endTotal)
        
        
        ! Write Final Data
        elapsed_time = elapsed_time + real((endTotal - startTotal), kind = real64) / real(timingRate, kind = real64)
        totPotTime = totPotTime + real(potentialTime, kind = real64) / real(timingRate, kind = real64)
        totMoverTime = totMoverTime + real(moverTime, kind = real64) / real(timingRate, kind = real64)
        totCollTime = totCollTime + real(collisionTime, kind = real64) / real(timingRate, kind = real64)
        open(9,file=directoryName//'/SimulationFinalData.dat')
        write(9,'("Elapsed Times(s), Potential Time (s), Mover Time (s), Collision Time (s), Total Steps")')
        write(9,"(4(es16.8,1x), 1(I8, 1x))") elapsed_time, totPotTime, totMoverTime, totCollTime, i
        close(9)
        print *, "Elapsed time for simulation is:", elapsed_time, "seconds"
        print *, 'moverTime is:', totMoverTime
        print *, 'potentialTime is:', totPotTime
        print *, 'collision time is:', totCollTime

    end subroutine solveSimulation


    subroutine solveSimulationFinalAverage(solver, particleList, targetParticleList, nullCollisionList, world, del_t, irand, averagingTime, binNumber)
        ! Perform averaging over certain time period for diagnostics
        type(Particle), intent(in out) :: particleList(numberChargedParticles)
        type(targetParticle), intent(in) :: targetParticleList(numberNeutralParticles)
        type(nullCollision), intent(in out) :: nullCollisionList(numberBinaryCollisions)
        type(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t, averagingTime
        integer(int32), intent(in) :: binNumber
        integer(c_int64_t), intent(in out) :: irand(numThread)
        integer(int32) :: i, stepsAverage, windowDivision, j, k, iThread, intPartV, sub_step
        real(real64) :: startTime, phi_average(NumberXNodes), chargeTotal, energyLoss, meanLoss, stdLoss, VHist(2*binNumber, numberChargedParticles), EHist(binNumber, numberChargedParticles), partV, partE
        real(real64) :: VMax(numberChargedParticles), EMax(numberChargedParticles), Emin(numberChargedParticles), &
        E_grid(binNumber, numberChargedParticles), diffE(numberChargedParticles), Emin_log(numberChargedParticles)
        real(real64) :: Efield_RF_energy_total, J_particles_heat_total, v_ave_heat_tot


        
        
        !Save initial particle/field data, along with domain
        do i = 1, numberChargedParticles
            particleList(i)%accumEnergyLoss = 0.0d0
            particleList(i)%accumWallLoss = 0
            particleList(i)%densities = 0
        end do
        inelasticEnergyLoss = 0.0d0
        energyAddColl = 0.0d0
        phi_average = 0.0d0
        Efield_RF_energy_total = 0.0d0
        J_particles_heat_total = 0.0d0
        v_ave_heat_tot = 0.0d0
        startTime = currentTime
        i = 0
        do j = 1, numberBinaryCollisions
            nullCollisionList(j)%totalEnergyLoss = 0
            nullCollisionList(j)%totalAmountCollisions = 0
            nullCollisionList(j)%totalIncidentEnergy = 0
            nullCollisionList(j)%totalNumberCollidableParticles = 0
        end do
        windowDivision = INT(200.0d0 / fractionFreq)
        do while(currentTime-startTime < averagingTime)
            ! move ions forward to ionStepMult * del_t
            if (numberChargedParticles > 1) call moveParticles(particleList(2:numberChargedParticles), solver%EField, world, ionStepMult * del_t)
            ! Now go extra (ionStepMult -1)/2 steps electrons
            
            do sub_step = 1, ionStepMult
                ! move electrons
                call moveParticles(particleList(1:1), solver%EField, world, del_t)
                currentTime = currentTime + del_t
                do j = 2, numberChargedParticles
                    ! For ionStepMult-1/2 steps, use ion rho from initial step, then accumulate through collisions
                    ! For center step, will then after accumulate total ion rho from ionStepMult*del_t by resetting startIdx
                    ! Then remaining two will again accumulate, with base rho from ionStepMult*del_t
                    if (sub_step /= (ionStepMult-1)/2 + 1) then
                        particleList(j)%startIdx = particleList(j)%N_p + 1
                    else
                        particleList(j)%startIdx = 1
                    end if
                end do
                
                ! undergo collisions for each time step
                if (particleList(1)%mass == m_e .and. numberChargedParticles > 1) then
                    if (heatingBool) call maxwellianHeating(particleList(1), irand, fractionFreq, T_e, del_t, del_t)
                    if (EField_heating_bool) call Efield_heating(particleList(1), currentTime - del_t, currentTime, world)
                    if (addLostPartBool) call addMaxwellianLostParticles(particleList, T_e, T_i, irand, world)
                    if (refluxPartBool) call refluxParticles(particleList, T_e, T_i, irand, world)
                    if (injectionBool) call injectAtBoundary(particleList, T_e, T_i, irand, world, del_t)
                    if (uniformInjectionBool) call injectUniformFlux(particleList, T_e, T_i, irand, world)
                end if
                do j = 1, numberBinaryCollisions
                    if (nullCollisionList(j)%reactantsIndx(1) == 1) then
                        call nullCollisionList(j)%generateCollision(particleList, targetParticleList, numberChargedParticles, numberBinaryCollisions, irand, del_t)  
                    else if (sub_step == ionStepMult) then
                        call nullCollisionList(j)%generateCollision(particleList, targetParticleList, numberChargedParticles, numberBinaryCollisions, irand, ionStepMult * del_t)
                    end if
                end do
        
                ! Re-evaluate rho with added particles
                call depositRho(solver, particleList, world)
                call solver%solve_tridiag_Poisson(world, currentTime)
                call solver%makeEField(world)

                call loadParticleDensity(particleList, world, .false.)
                phi_average = phi_average + solver%phi
                i = i + 1
            end do
        end do
        print *, "Averaging finished over", (currentTime - startTime), 'simulation time (s)'
        stepsAverage = i 
        phi_average = phi_average/stepsAverage
        J_particles_heat_total = J_particles_heat_total/stepsAverage
        v_ave_heat_tot = v_ave_heat_tot/stepsAverage
        call writePhi(phi_average, 0, .true., directoryName)
        chargeTotal = 0.0d0
        energyLoss = 0.0d0
        do i=1, numberChargedParticles
            chargeTotal = chargeTotal + SUM(particleList(i)%accumWallLoss) * particleList(i)%q * particleList(i)%w_p
            energyLoss = energyLoss + SUM(particleList(i)%accumEnergyLoss) * particleList(i)%mass * particleList(i)%w_p * 0.5d0
            !$OMP parallel private(iThread)
            iThread = omp_get_thread_num() + 1
            particleList(i)%densities(:, iThread) = particleList(i)%densities(:, iThread) /real(stepsAverage)
            !$OMP end parallel
            open(22,file=directoryName//'/ParticleAveDiagnostic_'//particleList(i)%name//'.dat')
            write(22, '("Left curr (A/m^2), right curr (A/m^2), left power (W/m^2), right power (W/m^2)")')
            write(22,"(4(es16.8, 1x))") particleList(i)%accumWallLoss(1) * particleList(i)%q_times_wp /(currentTime-startTime), &
                particleList(i)%accumWallLoss(2) * particleList(i)%q_times_wp /(currentTime-startTime), &
                particleList(i)%accumEnergyLoss(1)* particleList(i)%mass * particleList(i)%w_p * 0.5d0/(currentTime-startTime), &
                particleList(i)%accumEnergyLoss(2)* particleList(i)%mass * particleList(i)%w_p * 0.5d0/(currentTime-startTime)
            close(22)
        end do
        call writeParticleDensity(particleList, world, 0, .true., directoryName) 
        do j = 1, numberBinaryCollisions
            ! Save binary collision statistics
            inelasticEnergyLoss = inelasticEnergyLoss + SUM(nullCollisionList(j)%totalEnergyLoss) * 0.5d0 * particleList(nullCollisionList(j)%reactantsIndx(1))%w_p
            open(22,file=directoryName//'/BinaryCollisions/'//particleList(nullCollisionList(j)%reactantsIndx(1))%name//"_on_"//targetParticleList(nullCollisionList(j)%reactantsIndx(2))%name//"/AveCollisionDiag.dat")
            write(22, '("Coll #, CollRatio, AveEnergyLoss (eV), AveIncidentEnergy (eV), P_loss(W/m^2), aveCollFreq (Hz/m^2)")')
            do i = 1, nullCollisionList(j)%numberCollisions
                write(22,"((I3, 1x), 5(es16.8,1x))") i, real(nullCollisionList(j)%totalAmountCollisions(i))/real(nullCollisionList(j)%totalNumberCollidableParticles), &
                nullCollisionList(j)%totalEnergyLoss(i) * 0.5d0 / e/ real(nullCollisionList(j)%totalAmountCollisions(i)), &
                nullCollisionList(j)%totalIncidentEnergy(i) * particleList(nullCollisionList(j)%reactantsIndx(1))%mass * 0.5d0 / e/ real(nullCollisionList(j)%totalAmountCollisions(i)), &
                nullCollisionList(j)%totalEnergyLoss(i) * 0.5d0 * particleList(nullCollisionList(j)%reactantsIndx(1))%w_p / (currentTime-startTime), &
                real(nullCollisionList(j)%totalAmountCollisions(i)) * particleList(nullCollisionList(j)%reactantsIndx(1))%w_p / (currentTime-startTime)
            end do
            close(22)
        end do
        open(22,file=directoryName//'/GlobalDiagnosticDataAveraged.dat')
        write(22,'("Number Steps, Collision Loss (W/m^2), ParticleCurrentLoss (A/m^2), ParticlePowerLoss(W/m^2)")')
        write(22,"((I10, 1x), 3(es16.8,1x))") stepsAverage, inelasticEnergyLoss/(currentTime-startTime), chargeTotal/(currentTime-startTime), energyLoss/(currentTime-startTime)
        close(22)
        if (world%boundaryConditions(1) == 4) then
            print *, 'Final RF boundary values:', solver%phi(1)
        else if (world%boundaryConditions(NumberXNodes) == 4) then
            print *, 'Final RF boundary values:', solver%phi(NumberXNodes)
        end if
        if (EField_heating_bool) then
            print*, 'Power average for Efield heating is:', Efield_RF_energy_total / (currentTime - startTime), 'W'
            print*, 'Power per m^2 for particles is:', Efield_RF_energy_total * particleList(1)%w_p / (currentTime - startTime)
            print*, 'J average:', SQRT(J_particles_heat_total)
            print*, 'v average:', SQRT(v_ave_heat_tot)
        end if
        print *, 'Power loss to wall is:', energyLoss/(currentTime-startTime)
        print *, 'Power loss to inelastic collisions:', inelasticEnergyLoss/(currentTime-startTime)
        print *, "Electron average wall loss (A/m^2):", SUM(particleList(1)%accumWallLoss)* particleList(1)%w_p * particleList(1)%q/(currentTime-startTime)
        print *, "Ion average wall loss (A/m^2):", SUM(particleList(2)%accumWallLoss)* particleList(2)%w_p * particleList(2)%q/(currentTime-startTime)
        print *, "Performing average for EEDF over 50/omega_p"
        VMax = 0
        EMax = 0
        do j = 1, numberChargedParticles
            Emin(j) = 0.1d0 * e * 2.0d0 / particleList(j)%mass
        end do
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
        Emin_log = LOG(Emin)
        do j = 1, numberChargedParticles
            diffE(j) = (LOG(EMax(j)) - LOG(Emin(j)))/real(binNumber - 1)
            partV = LOG(Emin(j))
            do i = 1, binNumber
                E_grid(i, j) = EXP(partV)
                partV = partV + diffE(j)
            end do
        end do
      
        if (world%boundaryConditions(1) == 4 .or. world%boundaryConditions(NumberXNodes) == 4) then
            ! If RF boundary then average over RF cycle
            windowDivision = INT(2.0d0 * pi / solver%RF_rad_frequency/del_t/ionStepMult)
        else if (Efield_heating_bool) then
            windowDivision = INT(1.0d0 / Efield_heating_freq/del_t/ionStepMult)
        else
            windowDivision = INT(50.0d0/fractionFreq/ionStepMult)
        end if
        VHist = 0.0d0
        EHist = 0.0d0
        do i = 1, windowDivision
             ! move ions forward to ionStepMult * del_t
            if (numberChargedParticles > 1) call moveParticles(particleList(2:numberChargedParticles), solver%EField, world, ionStepMult * del_t)
            ! Now go extra (ionStepMult -1)/2 steps electrons
            
            do sub_step = 1, ionStepMult
                ! move electrons
                call moveParticles(particleList(1:1), solver%EField, world, del_t)
                currentTime = currentTime + del_t
                do j = 2, numberChargedParticles
                    ! For ionStepMult-1/2 steps, use ion rho from initial step, then accumulate through collisions
                    ! For center step, will then after accumulate total ion rho from ionStepMult*del_t by resetting startIdx
                    ! Then remaining two will again accumulate, with base rho from ionStepMult*del_t
                    if (sub_step /= (ionStepMult-1)/2 + 1) then
                        particleList(j)%startIdx = particleList(j)%N_p + 1
                    else
                        particleList(j)%startIdx = 1
                    end if
                end do
                
                ! undergo collisions for each time step
                if (particleList(1)%mass == m_e .and. numberChargedParticles > 1) then
                    if (heatingBool) call maxwellianHeating(particleList(1), irand, fractionFreq, T_e, del_t, del_t)
                    if (EField_heating_bool) call Efield_heating(particleList(1), currentTime - del_t, currentTime, world)
                    if (addLostPartBool) call addMaxwellianLostParticles(particleList, T_e, T_i, irand, world)
                    if (refluxPartBool) call refluxParticles(particleList, T_e, T_i, irand, world)
                    if (injectionBool) call injectAtBoundary(particleList, T_e, T_i, irand, world, del_t)
                    if (uniformInjectionBool) call injectUniformFlux(particleList, T_e, T_i, irand, world)
                end if
                do j = 1, numberBinaryCollisions
                    if (nullCollisionList(j)%reactantsIndx(1) == 1) then
                        call nullCollisionList(j)%generateCollision(particleList, targetParticleList, numberChargedParticles, numberBinaryCollisions, irand, del_t)  
                    else if (sub_step == ionStepMult) then
                        call nullCollisionList(j)%generateCollision(particleList, targetParticleList, numberChargedParticles, numberBinaryCollisions, irand, ionStepMult * del_t)
                    end if
                end do
        
                ! Re-evaluate rho with added particles
                call depositRho(solver, particleList, world)
                call solver%solve_tridiag_Poisson(world, currentTime)
                call solver%makeEField(world)

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
            end do
            
        end do
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