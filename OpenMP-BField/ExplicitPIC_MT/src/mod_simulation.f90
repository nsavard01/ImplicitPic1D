module mod_simulation
    ! module which actually uses others to run a simulation (single or multiple time steps)
    ! will contain diagnostics as well
    use iso_fortran_env, only: int32, real64, int64
    use constants
    use mod_BasicFunctions
    use mod_particle
    use mod_domain
    use mod_potentialSolver
    use mod_mt19937
    use mod_targetParticle
    use mod_NullCollision
    !use mod_collisions
    use omp_lib
    use ifport, only: makedirqq
    implicit none

    real(real64) :: inelasticEnergyLoss, currentTime
    real(real64), allocatable :: energyAddColl(:)

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




    subroutine loadParticleDensity(densities, particleList, world)
        type(Particle), intent(in) :: particleList(numberChargedParticles)
        type(Domain), intent(in) :: world
        real(real64), intent(in out) :: densities(NumberXNodes,numberChargedParticles)
        integer(int32) :: i,j, l_left, l_right, iThread
        real(real64) :: d, tempDensity(NumberXNodes, numThread)
        do i=1, numberChargedParticles
            tempDensity = 0.0d0
            !$OMP parallel private(iThread, j, l_left, l_right, d)
            iThread = omp_get_thread_num() + 1
            do j = 1, particleList(i)%N_p(iThread)
                l_left = INT(particleList(i)%phaseSpace(1,j, iThread))
                l_right = l_left + 1
                d = particleList(i)%phaseSpace(1,j, iThread) - l_left
                tempDensity(l_left, iThread) = tempDensity(l_left, iThread) + (1.0d0-d)
                tempDensity(l_right, iThread) = tempDensity(l_right, iThread) + d
            end do
            !$OMP end parallel
            densities(:, i) = densities(:, i) + SUM(tempDensity, DIM = 2) * particleList(i)%w_p
        end do

    end subroutine loadParticleDensity

    subroutine WriteParticleDensity(densities, particleList, world, CurrentDiagStep, boolAverage, dirName) 
        ! For diagnostics, deposit single particle density
        ! Re-use rho array since it doesn't get used after first Poisson
        real(real64), intent(in out) :: densities(NumberXNodes,numberChargedParticles)
        type(Particle), intent(in) :: particleList(numberChargedParticles)
        type(Domain), intent(in) :: world
        character(*), intent(in) :: dirName
        integer(int32), intent(in) :: CurrentDiagStep
        integer(int32) :: i
        logical, intent(in) :: boolAverage
        character(len=5) :: char_i
        
        do i=1, numberChargedParticles
            densities(:, i) = densities(:,i)/world%delX
            densities(1, i) = 2.0d0 * densities(1, i)
            densities(NumberXNodes, i) = 2.0d0 * densities(NumberXNodes, i)
            write(char_i, '(I3)'), CurrentDiagStep
            if (boolAverage) then
                open(41,file=dirName//'/Density/density_'//particleList(i)%name//"_Average.dat", form='UNFORMATTED')
            else
                open(41,file=dirName//'/Density/density_'//particleList(i)%name//"_"//trim(adjustl(char_i))//".dat", form='UNFORMATTED')
            end if
            write(41) densities(:, i)
            close(41)
        end do
        
    end subroutine WriteParticleDensity


    subroutine generateSaveDirectory(dirName)
        character(*), intent(in) :: dirName
        logical :: bool
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
            bool = makedirqq(dirName//'/CrossSections')
            if (.not. bool) then
                stop "Save directory not successfully created!"
            end if
        end if
    end subroutine generateSaveDirectory

    subroutine solveSimulation(solver, particleList, targetParticleList, nullCollisionList, world, del_t, randGen, simulationTime)
        ! Perform certain amount of timesteps, with diagnostics taken at first and last time step
        ! Impliment averaging for final select amount of timeSteps, this will be last data dump
        type(Particle), intent(in out) :: particleList(numberChargedParticles)
        type(targetParticle), intent(in) :: targetParticleList(numberNeutralParticles)
        type(nullCollision), intent(in out) :: nullCollisionList(numberBinaryCollisions)
        type(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t, simulationTime
        type(mt19937), intent(in out) :: randGen(numThread)
        integer(int32) :: i, j, CurrentDiagStep, diagStepDiff, unitPart1
        real(real64) :: densities(NumberXNodes, numberChargedParticles), diagTimeDivision, diagTime, elapsedTime, chargeTotal, energyLoss, elapsed_time
        integer(int64) :: startTime, endTime, timingRate, collisionTime, potentialTime, moverTime, startTotal, endTotal
        allocate(energyAddColl(numThread))
        CurrentDiagStep = 1
        unitPart1 = 100
        call generateSaveDirectory(directoryName)
        !Wrtie Initial conditions
        open(15,file=directoryName//'/InitialConditions.dat')
        write(15,'("Number Grid Nodes, Final Expected Time(s), Delta t(s), FractionFreq, delX, n_ave, T_e, T_i, numDiag, NumChargedPart, numThread, RF_frequency, RF_half_amplitude")')
        write(15,"((I6, 1x), 7(es16.8,1x), 3(I6, 1x), 2(es16.8,1x))") NumberXNodes, simulationTime, del_t, FractionFreq, world%delX, n_ave, T_e, T_i, numDiagnosticSteps, numberChargedParticles, numThread, solver%RF_rad_frequency, solver%RF_half_amplitude
        close(15)
        call system_clock(count_rate = timingRate)
        ! Write Particle properties
        open(9,file=directoryName//'/ParticleProperties.dat')
        write(9,'("Particle Symbol, Particle Mass (kg), Particle Charge (C), Particle Weight (N/m^2), maxIdx")')
        do j=1, numberChargedParticles
            write(9,"((A, 1x), 3(es16.8,1x), (I6, 1x))") particleList(j)%name//'       ', particleList(j)%mass, particleList(j)%q, particleList(j)%w_p, particleList(j)%finalIdx
        end do
        close(9)
        do i = 1, numberChargedParticles
            particleList(i)%accumEnergyLoss = 0.0d0
            particleList(i)%accumWallLoss = 0
            open(unitPart1+i,file=directoryName//'/ParticleDiagnostic_'//particleList(i)%name//'.dat')
            write(unitPart1+i,'("Time (s), leftCurrentLoss(A/m^2), rightCurrentLoss(A/m^2), leftPowerLoss(W/m^2), rightPowerLoss(W/m^2), N_p, Temp")')
        end do
        do i = 1, numberBinaryCollisions
            call nullCollisionList(i)%writeCollisionCrossSection(particleList(nullCollisionList(i)%reactantsIndx(1)), directoryName)
        end do
        collisionTime = 0
        potentialTime = 0
        moverTime = 0
        i = 0
        elapsedTime = 0.0d0
        diagStepDiff = 1
        diagTimeDivision = simulationTime/real(numDiagnosticSteps)
        diagTime = diagTimeDivision
        101 format(20(1x,es16.8))
        open(22,file=directoryName//'/GlobalDiagnosticData.dat')
        write(22,'("Time (s), Collision Loss (W/m^2), ParticleCurrentLoss (A/m^2), ParticlePowerLoss(W/m^2)")')
        
        !Save initial particle/field data, along with domain
        energyAddColl = 0.0d0
        inelasticEnergyLoss = 0.0d0
        !call solver%initialVRewind(particleList, del_t)
        densities = 0.0d0
        call loadParticleDensity(densities, particleList, world)
        call writeParticleDensity(densities, particleList, world, 0, .false., directoryName) 
        call writePhi(solver%phi, 0, .false., directoryName)
        call world%writeDomain(directoryName)
        do j=1, numberChargedParticles
            call particleList(j)%writeLocalTemperature(0, directoryName)
            !call particleList(j)%writePhaseSpace(0, directoryName)
        end do
        currentTime = 0.0d0
        call system_clock(startTotal)
        do while(currentTime < simulationTime)
            currentTime = currentTime + del_t
            if (currentTime < diagTime) then
                call system_clock(startTime)
                call solver%moveParticles(particleList, world, del_t)
                call system_clock(endTime)
                moverTime = moverTime + (endTime - startTime)
                call system_clock(startTime)
                call solver%solvePotential(particleList, world, currentTime)
                call system_clock(endTime)
                potentialTime = potentialTime + (endTime - startTime)
                !call ionizationCollisionIsotropic(particleList(1), particleList(2), 1.0d20, 1.0d-20, del_t, 15.8d0, 0.0d0, irand)
                call system_clock(startTime)
                do j = 1, numberBinaryCollisions
                    call nullCollisionList(j)%generateCollision(particleList, targetParticleList, numberChargedParticles, numberBinaryCollisions, randGen, del_t)
                end do
                call system_clock(endTime)
                collisionTime = collisionTime + (endTime - startTime)
            else  
                ! Data dump with diagnostics
                print *, "Simulation is", currentTime/simulationTime * 100.0, "percent done"
                call system_clock(startTime)
                call solver%moveParticles(particleList, world, del_t)
                call system_clock(endTime)
                moverTime = moverTime + (endTime - startTime)
                call system_clock(startTime)
                call solver%solvePotential(particleList, world, currentTime)
                call system_clock(endTime)
                potentialTime = potentialTime + (endTime - startTime)
                !call ionizationCollisionIsotropic(particleList(1), particleList(2), 1.0d20, 1.0d-20, del_t, 15.8d0, 0.0d0, irand)
                call system_clock(startTime)
                do j = 1, numberBinaryCollisions
                    call nullCollisionList(j)%generateCollision(particleList, targetParticleList, numberChargedParticles, numberBinaryCollisions, randGen, del_t)
                end do
                call system_clock(endTime)
                collisionTime = collisionTime + (endTime - startTime)
                densities = 0.0d0
                call loadParticleDensity(densities, particleList, world)
                call writeParticleDensity(densities, particleList, world, CurrentDiagStep, .false., directoryName) 
                call writePhi(solver%phi, CurrentDiagStep, .false., directoryName)
                chargeTotal = 0.0d0
                energyLoss = 0.0d0
                do j=1, numberChargedParticles
                    call particleList(j)%writeLocalTemperature(CurrentDiagStep, directoryName)
                    !call particleList(j)%writePhaseSpace(CurrentDiagStep, directoryName)
                    chargeTotal = chargeTotal + SUM(particleList(j)%accumWallLoss) * particleList(j)%q * particleList(j)%w_p
                    energyLoss = energyLoss + SUM(particleList(j)%accumEnergyLoss) * particleList(j)%mass * particleList(j)%w_p * 0.5d0
                    write(unitPart1+j,"(5(es16.8,1x), (I6,1x), (es16.8,1x))") currentTime, &
                        particleList(j)%accumWallLoss(1) * particleList(j)%q * particleList(j)%w_p/del_t/diagStepDiff, particleList(j)%accumWallLoss(2) * particleList(j)%q * particleList(j)%w_p/del_t/diagStepDiff, &
                        particleList(j)%accumEnergyLoss(1)* particleList(j)%mass * particleList(j)%w_p * 0.5d0/del_t/diagStepDiff, particleList(j)%accumEnergyLoss(2)* particleList(j)%mass * particleList(j)%w_p * 0.5d0/del_t/diagStepDiff, &
                        SUM(particleList(j)%N_p), particleList(j)%getKEAve()*2.0d0/3.0d0
                    particleList(j)%accumEnergyLoss = 0.0d0
                    particleList(j)%accumWallLoss = 0
                end do
                do j = 1, numberBinaryCollisions
                    inelasticEnergyLoss = inelasticEnergyLoss + nullCollisionList(j)%totalEnergyLoss * e * particleList(nullCollisionList(j)%reactantsIndx(1))%w_p
                end do
                write(22,"(4(es16.8,1x))") currentTime, inelasticEnergyLoss/del_t/diagStepDiff, &
                chargeTotal/del_t/diagStepDiff, energyLoss/del_t/diagStepDiff
                do j = 1, numberBinaryCollisions
                    nullCollisionList(j)%totalEnergyLoss = 0
                    nullCollisionList(j)%totalAmountCollisions = 0
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
        
        close(22)
        call system_clock(endTotal)
        elapsed_time = real((endTotal - startTotal), kind = real64) / real(timingRate, kind = real64)
        ! Write Final Data
        open(9,file=directoryName//'/SimulationFinalData.dat')
        write(9,'("Elapsed Times(s), Potential Time (s), Collision Time (s), Total Steps")')
        write(9,"(3(es16.8,1x), 1(I8, 1x))") elapsed_time, real(potentialTime + moverTime, kind = real64) / real(timingRate, kind = real64), real(collisionTime, kind = real64) / real(timingRate, kind = real64), i+1
        close(9)
        print *, "Elapsed time for simulation is:", elapsed_time, "seconds"
        print *, 'moverTime is:', real(moverTime, kind = real64)/real(timingRate, kind = real64)
        print *, 'potentialTime is:', real(potentialTime, kind = real64)/real(timingRate, kind = real64)
        print *, 'collision time is:', real(collisionTime, kind=real64)/real(timingRate, kind = real64)

    end subroutine solveSimulation


    subroutine solveSimulationFinalAverage(solver, particleList, targetParticleList, nullCollisionList, world, del_t, randGen, averagingTime, binNumber)
        ! Perform certain amount of timesteps, with diagnostics taken at first and last time step
        ! Impliment averaging for final select amount of timeSteps, this will be last data dump
        type(Particle), intent(in out) :: particleList(numberChargedParticles)
        type(targetParticle), intent(in) :: targetParticleList(numberNeutralParticles)
        type(nullCollision), intent(in out) :: nullCollisionList(numberBinaryCollisions)
        type(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t, averagingTime
        integer(int32), intent(in) :: binNumber
        type(mt19937), intent(in out) :: randGen(numThread)
        integer(int32) :: i, stepsAverage, windowNum, windowDivision, j, intPartV, VHist(2*binNumber), k
        real(real64) :: startTime, phi_average(NumberXNodes), densities(NumberXNodes, numberChargedParticles), chargeTotal, energyLoss, meanLoss, stdLoss
        real(real64) :: E_max, VMax
        real(real64), allocatable :: wallLoss(:)


        
        
        !Save initial particle/field data, along with domain
        do i = 1, numberChargedParticles
            particleList(i)%accumEnergyLoss = 0.0d0
            particleList(i)%accumWallLoss = 0
        end do
        inelasticEnergyLoss = 0.0d0
        energyAddColl = 0.0d0
        phi_average = 0.0d0
        densities = 0.0d0
        startTime = currentTime
        i = 0
        do j = 1, numberBinaryCollisions
            nullCollisionList(j)%totalEnergyLoss = 0
            nullCollisionList(j)%totalAmountCollisions = 0
        end do
        windowDivision = INT(200.0d0 / fractionFreq)
        allocate(wallLoss(2 * windowDivision))
        windowNum = 0
        do while(currentTime-startTime < averagingTime)
            currentTime = currentTime + del_t
            call solver%moveParticles(particleList, world, del_t)
            call solver%solvePotential(particleList, world, currentTime)
            !call ionizationCollisionIsotropic(particleList(1), particleList(2), 1.0d20, 1.0d-20, del_t, 15.8d0, 0.0d0, irand)
            do j = 1, numberBinaryCollisions
                call nullCollisionList(j)%generateCollision(particleList, targetParticleList, numberChargedParticles, numberBinaryCollisions, randGen, del_t)
            end do
            call loadParticleDensity(densities, particleList, world)
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
        densities = densities/i
        call writeParticleDensity(densities, particleList, world, 0, .true., directoryName) 
        phi_average = phi_average/stepsAverage
        call writePhi(phi_average, 0, .true., directoryName)
        chargeTotal = 0.0d0
        energyLoss = 0.0d0
        do i=1, numberChargedParticles
            chargeTotal = chargeTotal + SUM(particleList(i)%accumWallLoss) * particleList(i)%q * particleList(i)%w_p
            energyLoss = energyLoss + SUM(particleList(i)%accumEnergyLoss) * particleList(i)%mass * particleList(i)%w_p * 0.5d0
        end do
        do j = 1, numberBinaryCollisions
            inelasticEnergyLoss = inelasticEnergyLoss + nullCollisionList(j)%totalEnergyLoss * e * particleList(nullCollisionList(j)%reactantsIndx(1))%w_p
        end do
        open(22,file=directoryName//'/GlobalDiagnosticDataAveraged.dat')
        write(22,'("Number Steps, Collision Loss (W/m^2), ParticleCurrentLoss (A/m^2), ParticlePowerLoss(W/m^2)")')
        write(22,"((I8, 1x), 3(es16.8,1x))") stepsAverage, inelasticEnergyLoss/(currentTime-startTime), chargeTotal/(currentTime-startTime), energyLoss/(currentTime-startTime)
        close(22)
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
            do j = 1, numberBinaryCollisions
                call nullCollisionList(j)%generateCollision(particleList, targetParticleList, numberChargedParticles, numberBinaryCollisions, randGen, del_t)
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