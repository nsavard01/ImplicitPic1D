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
    real(real64) :: del_t

contains

    ! ------------------------- Reading Input data --------------------------------

    subroutine readInputs(NumberXNodes, maxIter, numDiagnosticSteps, stepsAverage, eps_r, fractionFreq, n_ave, T_e, T_i, L_domain, del_l, world, solver)
        ! Set initial conditions and global constants based on read input from txt file, create world and solver from these inputs
        integer(int32), intent(in out) :: NumberXNodes, maxIter, numDiagnosticSteps, stepsAverage
        real(real64), intent(in out) :: eps_r, fractionFreq, n_ave, T_e, T_i, L_domain, del_l
        integer(int32) :: io, leftBoundary, rightBoundary
        real(real64) :: leftVoltage, rightVoltage
        type(Domain) :: world
        type(potentialSolver) :: solver
        open(10,file='../InputData/InitialConditions.inp', IOSTAT=io)
        read(10, *, IOSTAT = io) T_e
        read(10, *, IOSTAT = io) T_i
        read(10, *, IOSTAT = io) eps_r
        read(10, *, IOSTAT = io) n_ave
        read(10, *, IOSTAT = io) numDiagnosticSteps
        read(10, *, IOSTAT = io) maxIter
        read(10, *, IOSTAT = io) fractionFreq
        read(10, *, IOSTAT = io) stepsAverage
        close(10)

        open(10,file='../InputData/Geometry.inp')
        read(10, *, IOSTAT = io) NumberXNodes
        read(10, *, IOSTAT = io) L_domain
        read(10, *, IOSTAT = io) del_l
        read(10, *, IOSTAT = io) leftBoundary, rightBoundary
        read(10, *, IOSTAT = io) leftVoltage, rightVoltage
        print *, "Left voltage is:", leftVoltage
        print *, "Right voltage is:", rightVoltage
        close(10)

        ! if one boundary is periodic, other must also be
        if ((leftBoundary == -3) .or. (rightBoundary == -3)) then
            leftBoundary = -3
            rightBoundary = -3
            leftVoltage = rightVoltage
        end if
        world = Domain()
        call world % constructSineGrid(del_l, L_domain)
        solver = potentialSolver(world, leftBoundary, rightBoundary, leftVoltage, rightVoltage)
        
    end subroutine readInputs

    function readParticleInputs(numberChargedParticles) result(particleList)
        type(Particle), allocatable :: particleList(:)
        integer(int32), intent(in out) :: numberChargedParticles
        integer(int32) :: j, numSpecies = 0, numParticles(100), particleIdxFactor(100)
        character(len=15) :: name
        character(len=8) :: particleNames(100)
        real(real64) :: mass(100), charge(100), Ti(100)

        print *, "Reading particle inputs"
        open(10,file='../InputData/BoundExample.dat', action = 'read')

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
                goto 20
            endif
            ! Take care of extra text I guess
20          do while(name(1:4).ne.'----')
                read(10,*,END=101,ERR=100) name
                print *, name
            end do
            
200         read(10,'(A6)',END=101,ERR=100, ADVANCE = 'NO') name
            if (name(1:4).eq.'----') then
                close(10)
            else
                numSpecies = numSpecies + 1
                read(10,*,END=101,ERR=100) mass(numSpecies),charge(numSpecies), Ti(numSpecies), numParticles(numSpecies), particleIdxFactor(numSpecies)
                mass(numSpecies) = mass(numSpecies) * m_p
                particleNames(numSpecies) = trim(name)
                goto 200
            end if



        end do
100     continue
101     continue
        numberChargedParticles = numSpecies
        allocate(particleList(numberChargedParticles))
        do j=1, numberChargedParticles
            particleList(j) = Particle(mass(j), e * charge(j), 1.0d0, numParticles(j), numParticles(j) * particleIdxFactor(j), trim(particleNames(j)))
            print *, mass(j), e * charge(j), numParticles(j), numParticles(j) * particleIdxFactor(j), trim(particleNames(j))
        end do
        



    end function readParticleInputs

    ! --------------------------- Diagnostics ------------------------------------

    subroutine depositRhoDiag(rho, particleList, world) 
        real(real64), intent(in out) :: rho(:)
        type(Particle), intent(in) :: particleList(:)
        type(Domain), intent(in) :: world
        integer(int32) :: i, j, l_left
        real(real64) :: d
        rho = 0.0d0
        do i=1, numberChargedParticles
            do j = 1, particleList(i)%N_p
                l_left = INT(particleList(i)%phaseSpace(1, j))
                d = MOD(particleList(i)%phaseSpace(1, j), 1.0d0)
                rho(l_left) = rho(l_left) + particleList(i)%q * particleList(i)%w_p * (1.0d0-d)
                rho(l_left + 1) = rho(l_left + 1) + particleList(i)%q * particleList(i)%w_p * d
            end do
        end do
        rho = rho / world%nodeVol
    end subroutine depositRhoDiag

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

    subroutine WriteParticleDensity(densities, particleList, world, CurrentDiagStep) 
        ! For diagnostics, deposit single particle density
        ! Re-use rho array since it doesn't get used after first Poisson
        real(real64), intent(in) :: densities(:,:)
        type(Particle), intent(in) :: particleList(:)
        type(Domain), intent(in) :: world
        integer(int32), intent(in) :: CurrentDiagStep
        integer(int32) :: i
        character(len=5) :: char_i
        
        do i=1, numberChargedParticles
            write(char_i, '(I3)'), CurrentDiagStep
            open(41,file='../Data/Density/density_'//particleList(i)%name//"_"//trim(adjustl(char_i))//".dat", form='UNFORMATTED')
            write(41) densities(:,i)/world%nodeVol
            close(41)
        end do
        
    end subroutine WriteParticleDensity

    subroutine writePhi(phi, CurrentDiagStep) 
        ! For diagnostics, deposit single particle density
        ! Re-use rho array since it doesn't get used after first Poisson
        real(real64), intent(in) :: phi(:)
        integer(int32), intent(in) :: CurrentDiagStep
        character(len=5) :: char_i
        write(char_i, '(I3)') CurrentDiagStep
        open(41,file='../Data/Phi/phi_'//trim(adjustl(char_i))//".dat", form='UNFORMATTED')
        write(41) phi
        close(41)
        
    end subroutine writePhi

    
    ! -------------------------- Simulation ------------------------------------------

    subroutine solveSingleTimeStep(solver, particleList, world, del_t, maxIter, eps_r, irand, boolDiagnostic)
        ! Single time step solver with Divergence of ampere, followed by adding of power, followed by collisions
        type(Particle), intent(in out) :: particleList(:)
        type(potentialSolver), intent(in out) :: solver
        logical, intent(in) :: boolDiagnostic
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t, eps_r
        integer(int32), intent(in) :: maxIter
        integer(int32), intent(in out) :: irand
        integer(int32) :: j, k
        real(real64) :: KE_i, KE_f, PE_i, PE_f, rho_f(NumberXNodes)

        if (boolDiagnostic) then

            ! Get charge/energy conservation error
            solver%particleEnergyLoss = 0.0d0
            solver%particleChargeLoss = 0.0d0
            PE_i = solver%getTotalPE(world, .false.)
            KE_i = 0.0d0
            do j=1, numberChargedParticles
                KE_i = KE_i + particleList(j)%getTotalKE()
            end do
            print *, "KE_i is:", KE_i
            call solver%depositRho(particleList, world) 
            call solver%adaptiveSolveDivAmperePicard(particleList, world, del_t, maxIter, eps_r)
            KE_f = solver%particleEnergyLoss
            do j=1, numberChargedParticles
                KE_f = KE_f + particleList(j)%getTotalKE()
            end do
            print *, "EnergyLoss is:", solver%particleEnergyLoss
            print *, "KE_f is:", KE_f
            PE_f = solver%getTotalPE(world, .false.)
            print *, "PE_i is:", PE_i
            print *, "PE_f is:", PE_f
            stop
            solver%energyError = ABS((KE_i + PE_i - KE_f - PE_f)/(KE_i + PE_i))
            call depositRhoDiag(rho_f, particleList, world)
            solver%chargeError = 0.0d0
            j = 0
            do k = 1, NumberXNodes -2
                if ((solver%J(k + 1) - solver%J(k) /= 0) .and. (solver%rho(k+1) /= 0)) then
                    solver%chargeError = solver%chargeError + (1 + (solver%J(k + 1) - solver%J(k)) *del_t/ world%nodeVol(k+1)/(rho_f(k+1) - solver%rho(k+1)))**2
                    j = j + 1
                end if
            end do
            solver%chargeError = SQRT(solver%chargeError/j)

        else
            call solver%adaptiveSolveDivAmperePicard(particleList, world, del_t, maxIter, eps_r)
        end if

        call addUniformPowerMaxwellian(particleList(1), Power, nu_h, irand, del_t)
        call ionizationCollisionIsotropic(particleList(1), particleList(2), 1.0d20, 1.0d-20, del_t, 15.8d0, irand)

    end subroutine solveSingleTimeStep


    subroutine solveSimulation(solver, particleList, world, del_t, maxIter, eps_r, irand, numTimeSteps, stepsAverage)
        ! Perform certain amount of timesteps, with diagnostics taken at first and last time step
        ! Impliment averaging for final select amount of timeSteps, this will be last data dump
        type(Particle), intent(in out) :: particleList(:)
        type(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t, eps_r
        integer(int32), intent(in) :: maxIter, numTimeSteps, stepsAverage
        integer(int32), intent(in out) :: irand
        integer(int32) :: numSkipSteps, i, j, CurrentDiagStep
        real(real64) :: currentTime, phi_average(NumberXNodes), densities(NumberXNodes, numberChargedParticles)
        CurrentDiagStep = 1

        open(15,file='../Data/InitialConditions.dat',access='APPEND')
        write(15,'("Number Grid Nodes, Final Expected Time(s), Delta t(s)")')
        write(15,"((I3.3, 1x), 2(es16.8,1x))") NumberXNodes, numTimeSteps*del_t, del_t
        close(15)
        close(22)
        numSkipSteps = numTimeSteps/(numDiagnosticSteps)
        101 format(20(1x,es16.8))
        open(22,file='../Data/GlobalDiagnosticData.dat',access='APPEND')
        write(22,'("Time (s), Collision Loss (W/m^2), ParticleCurrentLoss (A/m^2), ParticlePowerLoss(W/m^2), chargeError (a.u), energyError(a.u)")')
        
        !Save initial particle/field data, along with domain
        densities = 0.0d0
        call loadParticleDensity(densities, particleList)
        call writeParticleDensity(densities, particleList, world, 0) 
        call writePhi(solver%phi, 0)
        call particleList(1)%writeLocalTemperature(0)
        call world%writeDomain()
        do j=1, numberChargedParticles
            call particleList(j)%writePhaseSpace(0)
        end do

        do i = 1, numTimeSteps-stepsAverage
            if (MOD((i-1), numSkipSteps) /= 0) then
                call solveSingleTimeStep(solver, particleList, world, del_t, maxIter, eps_r, irand, .false.)
            else  
                ! Data dump with diagnostics
                print *, "Simulation is", real(i)/numTimeSteps * 100.0, "percent done"
                currentTime = (i) * del_t
                inelasticEnergyLoss = 0.0d0
                call solveSingleTimeStep(solver, particleList, world, del_t, maxIter, eps_r, irand, .true.)
                
                ! Stop program if catch abnormally large error
                if (solver%energyError > 1e-8) then
                    print *, "-------------------------WARNING------------------------"
                    print *, "Energy error is:", solver%energyError
                    stop "Total energy not conserved over time step in sub-step procedure!"
                end if
                
                if (solver%chargeError > 1e-6) then
                    print *, "-------------------------WARNING------------------------"
                    print *, "Charge error is:", solver%chargeError
                    stop "Total charge not conserved over time step in sub-step procedure!"
                end if

                densities = 0.0d0
                call loadParticleDensity(densities, particleList)
                call writeParticleDensity(densities, particleList, world, CurrentDiagStep) 
                call writePhi(solver%phi, CurrentDiagStep)
                call particleList(1)%writeLocalTemperature(CurrentDiagStep)
                do j=1, numberChargedParticles
                    call particleList(j)%writePhaseSpace(CurrentDiagStep)
                end do
                write(22,101) currentTime, inelasticEnergyLoss*e/del_t, SUM(solver%particleChargeLoss)/del_t, solver%particleEnergyLoss/del_t, solver%chargeError, solver%energyError
                CurrentDiagStep = CurrentDiagStep + 1
            end if
        end do
        currentTime = (numTimeSteps) * del_t
        print *, "Starting averaging in last", stepsAverage, "steps"
        solver%particleEnergyLoss = 0.0d0
        solver%particleChargeLoss = 0.0d0
        inelasticEnergyLoss = 0.0d0
        phi_average = 0.0d0
        densities = 0.0d0
        do i =1, stepsAverage-1
            call solveSingleTimeStep(solver, particleList, world, del_t, maxIter, eps_r, irand, .false.)
            call loadParticleDensity(densities, particleList)
            phi_average = phi_average + solver%phi
        end do
        call solveSingleTimeStep(solver, particleList, world, del_t, maxIter, eps_r, irand, .true.)
        call loadParticleDensity(densities, particleList)
        phi_average = phi_average + solver%phi
        call writeParticleDensity(densities/stepsAverage, particleList, world, CurrentDiagStep) 
        call writePhi(phi_average/stepsAverage, CurrentDiagStep)
        call particleList(1)%writeLocalTemperature(CurrentDiagStep)
        do j=1, numberChargedParticles
            call particleList(j)%writePhaseSpace(CurrentDiagStep)
        end do
        write(22,101) currentTime, inelasticEnergyLoss*e/del_t/stepsAverage, SUM(solver%particleChargeLoss)/del_t/stepsAverage, solver%particleEnergyLoss/del_t/stepsAverage, solver%chargeError, solver%energyError
        close(22)
        print *, "Electron average wall loss:", solver%particleChargeLoss(1)/del_t/stepsAverage
        print *, "Ion average wall loss:", solver%particleChargeLoss(2)/del_t/stepsAverage






    end subroutine solveSimulation







end module mod_simulation