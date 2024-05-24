module mod_simulation
    ! module which actually uses others to run a simulation (single or multiple time steps)
    ! will contain diagnostics as well
    use iso_fortran_env, only: int32, real64, int64
    use constants
    use mod_BasicFunctions
    use mod_particle
    use mod_domain
    use mod_potentialSolver
    use mod_particleInjection
    use mod_targetParticle
    use mod_NullCollision
    !use mod_collisions
    use omp_lib
    use ifport, only: makedirqq
    implicit none

    real(real64) :: inelasticEnergyLoss, currentTime

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

    ! --------------------------- Diagnostics ------------------------------------




    subroutine loadParticleDensity(particleList, world, reset)
        type(Particle), intent(in out) :: particleList(numberChargedParticles)
        type(Domain), intent(in) :: world
        logical, intent(in) :: reset
        integer(int32) :: i,j, l_left, l_right, iThread
        real(real64) :: d
        !$OMP parallel private(iThread, j, i, l_left, l_right, d)
        iThread = omp_get_thread_num() + 1
        do i=1, numberChargedParticles
            if (reset) then
                particleList(i)%densities(:,iThread) = 0.0d0
            end if
            do j = 1, particleList(i)%N_p(iThread)
                l_left = INT(particleList(i)%phaseSpace(1,j, iThread))
                l_right = l_left + 1
                d = particleList(i)%phaseSpace(1,j, iThread) - l_left
                particleList(i)%densities(l_left, iThread) = particleList(i)%densities(l_left, iThread) + (1.0d0-d)
                particleList(i)%densities(l_right, iThread) = particleList(i)%densities(l_right, iThread) + d
            end do
        end do
        !$OMP end parallel

    end subroutine loadParticleDensity

    subroutine WriteParticleDensity(particleList, world, CurrentDiagStep, boolAverage, dirName) 
        ! For diagnostics, deposit single particle density
        ! Re-use rho array since it doesn't get used after first Poisson
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
            densities(1) = 2.0d0 * densities(1)
            densities(NumberXNodes) = 2.0d0 * densities(NumberXNodes)
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
            call execute_command_line("rm -r "//dirName//"/*", EXITSTAT = io)
            if (io /= 0) then
                stop 'Issue removing old data'
            end if
        end if
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
        if (numberBinaryCollisions > 0) then
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
        type(targetParticle), intent(in) :: targetParticleList(numberNeutralParticles)
        type(nullCollision), intent(in out) :: nullCollisionList(numberBinaryCollisions)
        type(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t, simulationTime
        integer(int32), intent(in out) :: irand(numThread)
        integer(int32) :: i, j, CurrentDiagStep, diagStepDiff, unitPart1
        real(real64) :: diagTimeDivision, diagTime, elapsedTime, chargeTotal, energyLoss, elapsed_time, momentum_total(3), totMoverTime, totPotTime, totCollTime
        integer(int64) :: startTime, endTime, timingRate, collisionTime, potentialTime, moverTime, startTotal, endTotal
        character(len=8) :: date_char
        character(len=10) :: time_char
        allocate(energyAddColl(numThread))
        CurrentDiagStep = oldDiagStep
        unitPart1 = 100
        if (.not. restartBool) call generateSaveDirectory(directoryName)
        !Wrtie Initial conditions
        call date_and_time(DATE=date_char, TIME = time_char)
        open(15,file=directoryName//'/DateTime.dat')
        write(15,'("UTC Date, UTC Time")')
        write(15,"(2(A,1x))") date_char, time_char 
        close(15)
        open(15,file=directoryName//'/InitialConditions.dat')
        write(15,'("Number Grid Nodes, Final Expected Time(s), Delta t(s), FractionFreq, delX, n_ave, T_e, T_i, numDiag, NumChargedPart, numThread, RF_frequency, RF_half_amplitude")')
        write(15,"((I6, 1x), 7(es16.8,1x), 3(I6, 1x), 2(es16.8,1x))") NumberXNodes, simulationTime, del_t, FractionFreq, world%delX, n_ave, T_e, T_i, numDiagnosticSteps, numberChargedParticles, numThread, solver%RF_rad_frequency, solver%RF_half_amplitude
        close(15)
        call system_clock(count_rate = timingRate)
        ! Write Particle properties
        open(9,file=directoryName//'/ParticleProperties.dat')
        write(9,'("Particle Symbol, Particle Mass (kg), Particle Charge (C), Particle Weight (N/m^2), maxIdx")')
        do j=1, numberChargedParticles
            write(9,"((A, 1x), 3(es16.8,1x), (I10, 1x))") particleList(j)%name//'       ', particleList(j)%mass, particleList(j)%q, particleList(j)%w_p, particleList(j)%finalIdx
        end do
        close(9)
        do i = 1, numberChargedParticles
            particleList(i)%accumEnergyLoss = 0.0d0
            particleList(i)%accumWallLoss = 0
            open(unitPart1+i,file=directoryName//'/ParticleDiagnostic_'//particleList(i)%name//'.dat', access = 'APPEND')
            if (.not. restartBool) write(unitPart1+i,'("Time (s), leftCurrentLoss(A/m^2), rightCurrentLoss(A/m^2), leftPowerLoss(W/m^2), rightPowerLoss(W/m^2), N_p, Temp")')
        end do
        do i = 1, numberBinaryCollisions
            if (.not. restartBool) call nullCollisionList(i)%writeCollisionProperties(particleList, targetParticleList, directoryName)
        end do

        if (restartBool) then
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

        
        
        !Save initial particle/field data, along with domain
        energyAddColl = 0.0d0
        inelasticEnergyLoss = 0.0d0
        !call solver%initialVRewind(particleList, del_t)
        call loadParticleDensity(particleList, world, .true.)
        call writeParticleDensity(particleList, world, 0, .false., directoryName)  
        call writePhi(solver%phi, 0, .false., directoryName)
        call world%writeDomain(directoryName)
        do j=1, numberChargedParticles
            particleList(j)%momentumLoss = 0.0d0
            call particleList(j)%writeLocalTemperature(0, directoryName, NumberXHalfNodes)
            !call particleList(j)%writePhaseSpace(0, directoryName)
        end do
        call system_clock(startTotal)
        do while(currentTime < simulationTime)
            currentTime = currentTime + del_t
            if (currentTime < diagTime) then
                call system_clock(startTime)
                call solver%moveParticles(particleList, world, del_t)
                call solver%depositRho(particleList, world)
                call system_clock(endTime)
                moverTime = moverTime + (endTime - startTime)
                call system_clock(startTime)
                call solver%solve_tridiag_Poisson(world, currentTime)
                call solver%makeEField(world)
                call system_clock(endTime)
                potentialTime = potentialTime + (endTime - startTime)
                !call ionizationCollisionIsotropic(particleList(1), particleList(2), 1.0d20, 1.0d-20, del_t, 15.8d0, 0.0d0, irand)
                call system_clock(startTime)
                if (heatingBool) call maxwellianHeating(particleList(1), irand, fractionFreq, T_e, del_t, del_t)
                if (addLostPartBool) call addMaxwellianLostParticles(particleList, T_e, T_i, irand, world)
                if (refluxPartBool) call refluxParticles(particleList, T_e, T_i, irand, world)
                if (injectionBool) call injectAtBoundary(particleList, T_e, T_i, irand, world, del_t)
                if (uniformInjectionBool) call injectUniformFlux(particleList, T_e, T_i, irand, world)
                do j = 1, numberBinaryCollisions
                    call nullCollisionList(j)%generateCollision(particleList, targetParticleList, numberChargedParticles, numberBinaryCollisions, irand, del_t)
                end do
                call system_clock(endTime)
                collisionTime = collisionTime + (endTime - startTime)
            else  
                ! Data dump with diagnostics
                print *, "Simulation is", currentTime/simulationTime * 100.0, "percent done"
                call system_clock(startTime)
                call solver%moveParticles(particleList, world, del_t)
                call solver%depositRho(particleList, world)
                call system_clock(endTime)
                moverTime = moverTime + (endTime - startTime)
                call system_clock(startTime)
                call solver%solve_tridiag_Poisson(world, currentTime)
                call solver%makeEField(world)
                call system_clock(endTime)
                potentialTime = potentialTime + (endTime - startTime)
                !call ionizationCollisionIsotropic(particleList(1), particleList(2), 1.0d20, 1.0d-20, del_t, 15.8d0, 0.0d0, irand)
                call system_clock(startTime)
                if (heatingBool) call maxwellianHeating(particleList(1), irand, fractionFreq, T_e, del_t, del_t)
                if (addLostPartBool) call addMaxwellianLostParticles(particleList, T_e, T_i, irand, world)
                if (refluxPartBool) call refluxParticles(particleList, T_e, T_i, irand, world)
                if (injectionBool) call injectAtBoundary(particleList, T_e, T_i, irand, world, del_t)
                if (uniformInjectionBool) call injectUniformFlux(particleList, T_e, T_i, irand, world)
                do j = 1, numberBinaryCollisions
                    call nullCollisionList(j)%generateCollision(particleList, targetParticleList, numberChargedParticles, numberBinaryCollisions, irand, del_t)
                end do
                call system_clock(endTime)
                collisionTime = collisionTime + (endTime - startTime)
                call loadParticleDensity(particleList, world, .true.)
                call writeParticleDensity(particleList, world, CurrentDiagStep, .false., directoryName)  
                call writePhi(solver%phi, CurrentDiagStep, .false., directoryName)
                chargeTotal = 0.0d0
                energyLoss = 0.0d0
                momentum_total = 0.0d0
                do j=1, numberChargedParticles
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
                    inelasticEnergyLoss = inelasticEnergyLoss + 0.5d0 * SUM(nullCollisionList(j)%totalEnergyLoss) * particleList(nullCollisionList(j)%reactantsIndx(1))%w_p
                    call nullCollisionList(j)%writeCollisionDiag(particleList, targetParticleList, directoryName, del_t * diagStepDiff)
                end do
                write(202,"(6(es16.8,1x))") currentTime, inelasticEnergyLoss/del_t/diagStepDiff, &
                chargeTotal/del_t/diagStepDiff, energyLoss/del_t/diagStepDiff, momentum_total(1), solver%getTotalE(particleList, world)
                call system_clock(endTotal)
                write(303,"(4(es16.8,1x), 1(I8, 1x))") elapsed_time + real((endTotal - startTotal), kind = real64) / real(timingRate, kind = real64), &
                totPotTime + real(potentialTime, kind = real64) / real(timingRate, kind = real64), &
                totMoverTime + real(moverTime, kind = real64) / real(timingRate, kind = real64), &
                totCollTime + real(collisionTime, kind = real64) / real(timingRate, kind = real64), i+1

                do j = 1, numberBinaryCollisions
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
            ! if (MODULO(i+1, heatSkipSteps) == 0) then
            !     call addUniformPowerMaxwellian(particleList(1), Power, nu_h, irand, heatSkipSteps*del_t)
            ! end if
            i = i + 1
            diagStepDiff = diagStepDiff + 1
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
        write(9,"(4(es16.8,1x), 1(I8, 1x))") elapsed_time, totPotTime, totMoverTime, totCollTime, i+1
        close(9)
        print *, "Elapsed time for simulation is:", elapsed_time, "seconds"
        print *, 'moverTime is:', totMoverTime
        print *, 'potentialTime is:', totPotTime
        print *, 'collision time is:', totCollTime

    end subroutine solveSimulation


    subroutine solveSimulationFinalAverage(solver, particleList, targetParticleList, nullCollisionList, world, del_t, irand, averagingTime, binNumber)
        ! Perform certain amount of timesteps, with diagnostics taken at first and last time step
        ! Impliment averaging for final select amount of timeSteps, this will be last data dump
        type(Particle), intent(in out) :: particleList(numberChargedParticles)
        type(targetParticle), intent(in) :: targetParticleList(numberNeutralParticles)
        type(nullCollision), intent(in out) :: nullCollisionList(numberBinaryCollisions)
        type(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t, averagingTime
        integer(int32), intent(in) :: binNumber
        integer(int32), intent(in out) :: irand(numThread)
        integer(int32) :: i, stepsAverage, windowNum, windowDivision, j, intPartV, VHist(2*binNumber), k, iThread
        real(real64) :: startTime, phi_average(NumberXNodes), chargeTotal, energyLoss, meanLoss, stdLoss
        real(real64) :: E_max, VMax
        real(real64), allocatable :: wallLoss(:)


        
        
        !Save initial particle/field data, along with domain
        do i = 1, numberChargedParticles
            particleList(i)%accumEnergyLoss = 0.0d0
            particleList(i)%accumWallLoss = 0
            particleList(i)%densities = 0
        end do
        inelasticEnergyLoss = 0.0d0
        energyAddColl = 0.0d0
        phi_average = 0.0d0
        startTime = currentTime
        i = 0
        do j = 1, numberBinaryCollisions
            nullCollisionList(j)%totalEnergyLoss = 0
            nullCollisionList(j)%totalAmountCollisions = 0
            nullCollisionList(j)%totalIncidentEnergy = 0
            nullCollisionList(j)%totalNumberCollidableParticles = 0
        end do
        windowDivision = INT(200.0d0 / fractionFreq)
        allocate(wallLoss(2 * windowDivision))
        windowNum = 0
        do while(currentTime-startTime < averagingTime)
            currentTime = currentTime + del_t
            call solver%moveParticles(particleList, world, del_t)
            call solver%solvePotential(particleList, world, currentTime)
            !call ionizationCollisionIsotropic(particleList(1), particleList(2), 1.0d20, 1.0d-20, del_t, 15.8d0, 0.0d0, irand)
            if (heatingBool) call maxwellianHeating(particleList(1), irand, fractionFreq, T_e, del_t, del_t)
            if (addLostPartBool) call addMaxwellianLostParticles(particleList, T_e, T_i, irand, world)
            if (refluxPartBool) call refluxParticles(particleList, T_e, T_i, irand, world)
            if (injectionBool) call injectAtBoundary(particleList, T_e, T_i, irand, world, del_t)
            if (uniformInjectionBool) call injectUniformFlux(particleList, T_e, T_i, irand, world)
            do j = 1, numberBinaryCollisions
                call nullCollisionList(j)%generateCollision(particleList, targetParticleList, numberChargedParticles, numberBinaryCollisions, irand, del_t)
            end do
            call loadParticleDensity(particleList, world, .false.)
            phi_average = phi_average + solver%phi
            ! if (MODULO(i, heatSkipSteps) == 0) then
            !     call addUniformPowerMaxwellian(particleList(1), Power, nu_h, irand, heatSksipSteps*del_t)
            ! end if
            i = i + 1
            windowNum = windowNum + 1
            wallLoss(windowNum) = 0.0d0
            do j = 1, numberChargedParticles
                wallLoss(windowNum) = wallLoss(windowNum) + SUM(particleList(j)%accumEnergyLoss)
            end do
            wallLoss(windowNum) = wallLoss(windowNum)/(currentTime - startTime)
            if (windowNum > windowDivision) then
                meanLoss = SUM(wallLoss(1:windowNum))/real(windowNum)
                stdLoss = SQRT(SUM( (wallLoss(1:windowNum) - meanLoss)**2 )/real(windowNum))
                if (stdLoss/meanLoss < 1d-6) exit
                windowNum = 0
            end if
        end do
        print *, "Averaging finished over", (currentTime - startTime), 'simulation time (s)'
        stepsAverage = i 
        phi_average = phi_average/stepsAverage
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
        end do
        call writeParticleDensity(particleList, world, 0, .true., directoryName) 
        do j = 1, numberBinaryCollisions
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
        
        print *, 'Power loss to wall is:', energyLoss/(currentTime-startTime)
        print *, 'Power loss to inelastic collisions:', inelasticEnergyLoss/(currentTime-startTime)
        print *, "Electron average wall loss (A/m^2):", SUM(particleList(1)%accumWallLoss)* particleList(1)%w_p * particleList(1)%q/(currentTime-startTime)
        print *, "Ion average wall loss (A/m^2):", SUM(particleList(2)%accumWallLoss)* particleList(2)%w_p * particleList(2)%q/(currentTime-startTime)
        print *, "Performing average for EEDF over 50/omega_p"
        E_max = 3.0d0 * (MAXVAL(phi_average) - minval(phi_average))
        VMax = SQRT(2.0d0 * E_max *e/ m_e)
        windowDivision = INT(50.0d0/fractionFreq)
        VHist = 0.0d0
        do i = 1, windowDivision
            currentTime = currentTime + del_t
            call solver%moveParticles(particleList, world, del_t)
            call solver%solvePotential(particleList, world, currentTime)
            if (heatingBool) call maxwellianHeating(particleList(1), irand, fractionFreq, T_e, del_t, del_t)
            if (addLostPartBool) call addMaxwellianLostParticles(particleList, T_e, T_i, irand, world)
            if (refluxPartBool) call refluxParticles(particleList, T_e, T_i, irand, world)
            if (injectionBool) call injectAtBoundary(particleList, T_e, T_i, irand, world, del_t)
            do j = 1, numberBinaryCollisions
                call nullCollisionList(j)%generateCollision(particleList, targetParticleList, numberChargedParticles, numberBinaryCollisions, irand, del_t)
            end do
            do k = 1, numThread
                do j = 1, particleList(1)%N_p(k)
                    intPartV = INT(particleList(1)%phaseSpace(2, j, k) * (binNumber) / VMax + binNumber + 1)
                    if (intPartV > 0 .and. intPartV < binNumber*2+1) VHist(intPartV) = VHist(intPartV) + 1
                end do
            end do
        end do
        open(10,file=directoryName//'/ElectronTemperature/EVDF_average.dat', form='UNFORMATTED')
        write(10) real(VHist, kind = real64), VMax
        close(10)

    end subroutine solveSimulationFinalAverage




end module mod_simulation