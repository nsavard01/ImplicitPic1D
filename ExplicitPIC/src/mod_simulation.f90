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

    subroutine readInputs(NumberXNodes, maxIter, numDiagnosticSteps, stepsAverage, eps_r, fractionFreq, n_ave, world, solver)
        ! Set initial conditions and global constants based on read input from txt file, create world and solver from these inputs
        integer(int32), intent(in out) :: NumberXNodes, maxIter, numDiagnosticSteps, stepsAverage
        real(real64), intent(in out) :: eps_r, fractionFreq, n_ave
        integer(int32) :: io, leftBoundary, rightBoundary, gridType
        real(real64) :: leftVoltage, rightVoltage, L_domain, del_l
        type(Domain) :: world
        type(potentialSolver) :: solver

        print *, "Reading solver inputs:"
        open(10,file='../InputData/InitialConditions.inp', IOSTAT=io)
        read(10, *, IOSTAT = io) eps_r
        read(10, *, IOSTAT = io) n_ave
        read(10, *, IOSTAT = io) numDiagnosticSteps
        read(10, *, IOSTAT = io) maxIter
        read(10, *, IOSTAT = io) fractionFreq
        read(10, *, IOSTAT = io) stepsAverage
        close(10)
        print *, "Relative error:", eps_r
        print *, "Average initial particle density:", n_ave
        print *, "Number of diagnostic steps is:", numDiagnosticSteps
        print *, "Maximum iteration number:", maxIter
        print *, "Fraction of 1/w_p for time step:", fractionFreq
        print *, "Number of final steps for averaging:", stepsAverage
        print *, "------------------"
        print *, ""
        print *, "Reading domain inputs:"
        open(10,file='../InputData/Geometry.inp')
        read(10, *, IOSTAT = io) NumberXNodes
        read(10, *, IOSTAT = io) L_domain
        read(10, *, IOSTAT = io) gridType
        read(10, *, IOSTAT = io) del_l
        read(10, *, IOSTAT = io) leftBoundary, rightBoundary
        read(10, *, IOSTAT = io) leftVoltage, rightVoltage
        close(10)
        print *, "Number of nodes:", NumberXNodes
        print *, "Grid length:", L_domain
        print *, "Type of grid:", gridType
        print *, "Minimum logical space unit if sinusoidal grid:", del_l
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
        call world % constructGrid(del_l, L_domain, gridType)
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

    ! --------------------------- Diagnostics ------------------------------------

    subroutine depositRhoDiag(rho, particleList, world) 
        real(real64), intent(in out) :: rho(:)
        type(Particle), intent(in) :: particleList(:)
        type(Domain), intent(in) :: world
        integer(int32) :: i, j, l_center
        real(real64) :: d
        rho = 0.0d0
        do i=1, numberChargedParticles
            do j = 1, particleList(i)%N_p
                l_center = NINT(particleList(i)%phaseSpace(1, j))
                d = particleList(i)%phaseSpace(1, j) - l_center
                if (world%boundaryConditions(l_center) == 0) then
                    ! Inside domain
                    rho(l_center) = rho(l_center) + particleList(i)%q * particleList(i)%w_p * (0.75 - d**2)
                    rho(l_center + 1) = rho(l_center + 1) + particleList(i)%q * particleList(i)%w_p * 0.5d0 * (0.5d0 + d)**2
                    rho(l_center - 1) = rho(l_center - 1) + particleList(i)%q * particleList(i)%w_p * 0.5d0 * (0.5d0 - d)**2
                else if (world%boundaryConditions(l_center) == 1) then
                    !Dirichlet
                    rho(l_center) = rho(l_center) + particleList(i)%q * particleList(i)%w_p * (1.0d0-ABS(d))
                    rho(l_center + INT(SIGN(1.0, d))) = rho(l_center + INT(SIGN(1.0, d))) + particleList(i)%q * particleList(i)%w_p * ABS(d)
                
                else if (world%boundaryConditions(l_center) == 3) then
                    ! Periodic
                    rho(l_center) = rho(l_center) + particleList(i)%q * particleList(i)%w_p * (0.75 - d**2)
                    ! towards domain
                    rho(l_center+INT(SIGN(1.0, d))) = rho(l_center+INT(SIGN(1.0, d))) + particleList(i)%q * particleList(i)%w_p * 0.5d0 * (0.5d0 + ABS(d))**2
                    ! across periodic boundary
                    rho(MODULO(l_center-2*INT(SIGN(1.0, d)),NumberXNodes)) = rho(MODULO(l_center-2*INT(SIGN(1.0, d)),NumberXNodes)) + particleList(i)%q * particleList(i)%w_p * 0.5d0 * (0.5d0 - ABS(d))**2
                end if
            end do
        end do
        if (world%boundaryConditions(1) == 3) then
            rho(1) = rho(1) + rho(NumberXNodes)
            rho(NumberXNodes) = rho(1)
        end if
        rho = rho / world%dx_dl
    end subroutine depositRhoDiag

    subroutine loadParticleDensity(densities, particleList, world)
        type(Particle), intent(in) :: particleList(:)
        type(Domain), intent(in) :: world
        real(real64), intent(in out) :: densities(:,:)
        integer(int32) :: i,j, l_center
        real(real64) :: d
        do i=1, numberChargedParticles
            do j = 1, particleList(i)%N_p
                l_center = NINT(particleList(i)%phaseSpace(1, j))
                d = particleList(i)%phaseSpace(1, j) - l_center
                if (world%boundaryConditions(l_center) == 0) then
                    ! Inside domain
                    densities(l_center, i) = densities(l_center, i) + particleList(i)%w_p * (0.75 - d**2)
                    densities(l_center + 1, i) = densities(l_center + 1, i) + particleList(i)%w_p * 0.5d0 * (0.5d0 + d)**2
                    densities(l_center - 1, i) = densities(l_center - 1, i) + particleList(i)%w_p * 0.5d0 * (0.5d0 - d)**2
                else if (world%boundaryConditions(l_center) == 1) then
                    !Dirichlet
                    densities(l_center, i) = densities(l_center, i) + particleList(i)%w_p * (1.0d0-ABS(d))
                    densities(l_center + INT(SIGN(1.0, d)), i) = densities(l_center + INT(SIGN(1.0, d)), i) + particleList(i)%w_p * ABS(d)
                
                else if (world%boundaryConditions(l_center) == 3) then
                    ! Periodic
                    densities(l_center, i) = densities(l_center, i) + particleList(i)%w_p * (0.75 - d**2)
                    ! towards domain
                    densities(l_center+INT(SIGN(1.0, d)), i) = densities(l_center+INT(SIGN(1.0, d)), i) + particleList(i)%w_p * 0.5d0 * (0.5d0 + ABS(d))**2
                    ! across periodic boundary
                    densities(MODULO(l_center-2*INT(SIGN(1.0, d)),NumberXNodes), i) = densities(MODULO(l_center-2*INT(SIGN(1.0, d)),NumberXNodes), i) + particleList(i)%w_p * 0.5d0 * (0.5d0 - ABS(d))**2
                end if
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
            if (world%boundaryConditions(1) == 3) then
                densities(1,i) = densities(1,i) + densities(NumberXNodes, i)
                densities(NumberXNodes, i) = densities(1, i)
            end if
            densities(:,i) = densities(:,i)/world%dx_dl
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
        integer(int32) :: j, k
        real(real64) :: KE_i, KE_f, PE_i, PE_f, rho_f(NumberXNodes)

        ! Get charge/energy conservation error
        solver%particleEnergyLoss = 0.0d0
        PE_i = solver%getTotalPE(world, .false.)
        KE_i = 0.0d0
        do j=1, numberChargedParticles
            KE_i = KE_i + particleList(j)%getTotalKE()
        end do
        call solver%depositRho(particleList, world) 
        call solver%adaptiveSolveDivAmperePicard(particleList, world, del_t, maxIter, eps_r)
        KE_f = solver%particleEnergyLoss
        do j=1, numberChargedParticles
            KE_f = KE_f + particleList(j)%getTotalKE()
        end do
        PE_f = solver%getTotalPE(world, .false.)
        solver%energyError = ABS((KE_i + PE_i - KE_f - PE_f)/(KE_i + PE_i))
        call depositRhoDiag(rho_f, particleList, world)
        solver%chargeError = 0.0d0
        j = 0
        if (world%boundaryConditions(1) == 3) then
            j = j + 1
            solver%chargeError = solver%chargeError + (1 + (solver%J(1) - solver%J(NumberXNodes-1)) *del_t/ world%dx_dl(1)/(rho_f(1) - solver%rho(1)))**2
        end if
        do k = 1, NumberXNodes -2
            if ((solver%rho(k+1) /= 0)) then
                solver%chargeError = solver%chargeError + (1 + (solver%J(k + 1) - solver%J(k)) *del_t/ world%dx_dl(k+1)/(rho_f(k+1) - solver%rho(k+1)))**2
                j = j + 1
            end if
        end do
        solver%chargeError = SQRT(solver%chargeError/j)

    end subroutine solveSingleTimeStepDiagnostic

    subroutine solveSimulationOnlyPotential(solver, particleList, world, del_t, maxIter, eps_r, numTimeSteps)
        ! Perform certain amount of timesteps, with diagnostics taken at first and last time step
        ! Impliment averaging for final select amount of timeSteps, this will be last data dump
        type(Particle), intent(in out) :: particleList(:)
        type(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t, eps_r
        integer(int32), intent(in) :: maxIter, numTimeSteps
        integer(int32) :: numSkipSteps, i, j, CurrentDiagStep
        real(real64) :: currentTime, densities(NumberXNodes, numberChargedParticles)
        CurrentDiagStep = 1

        !Wrtie Initial conditions
        open(15,file='../Data/InitialConditions.dat')
        write(15,'("Number Grid Nodes, Final Expected Time(s), Delta t(s), FractionFreq")')
        write(15,"((I3.3, 1x), 3(es16.8,1x))") NumberXNodes, numTimeSteps*del_t, del_t, FractionFreq
        close(15)

        ! Write Particle properties
        open(9,file='../Data/ParticleProperties.dat')
        write(9,'("Particle Symbol, Particle Mass (kg), Particle Charge (C), Particle Weight (N/m^2)")')
        do j=1, numberChargedParticles
            write(9,"((A, 1x), 3(es16.8,1x))") particleList(j)%name, particleList(j)%mass, particleList(j)%q, particleList(j)%w_p
        end do
        close(9)

        numSkipSteps = numTimeSteps/(numDiagnosticSteps)
        101 format(20(1x,es16.8))
        open(22,file='../Data/GlobalDiagnosticData.dat')
        write(22,'("Time (s), Collision Loss (W/m^2), ParticleCurrentLoss (A/m^2), ParticlePowerLoss(W/m^2), chargeError (a.u), energyError(a.u), Picard Iteration Number")')
        
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

        do i = 1, numTimeSteps-1
            if (MOD((i-1), numSkipSteps) /= 0) then
                call solver%adaptiveSolveDivAmperePicard(particleList, world, del_t, maxIter, eps_r)
            else  
                ! Data dump with diagnostics
                print *, "Simulation is", real(i)/numTimeSteps * 100.0, "percent done"
                currentTime = (i) * del_t
                inelasticEnergyLoss = 0.0d0
                solver%particleEnergyLoss = 0.0d0
                call solveSingleTimeStepDiagnostic(solver, particleList, world, del_t, maxIter, eps_r)
                
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
                call loadParticleDensity(densities, particleList, world)
                call writeParticleDensity(densities, particleList, world, CurrentDiagStep, .false.) 
                call writePhi(solver%phi, CurrentDiagStep, .false.)
                call particleList(1)%writeLocalTemperature(CurrentDiagStep)
                do j=1, numberChargedParticles
                    call particleList(j)%writePhaseSpace(CurrentDiagStep)
                end do
                write(22,"(6(es16.8,1x), (I4, 1x))") currentTime, inelasticEnergyLoss/del_t, &
                SUM(solver%particleChargeLoss)/del_t, solver%particleEnergyLoss/del_t, solver%chargeError, solver%energyError, solver%iterNumPicard
                CurrentDiagStep = CurrentDiagStep + 1
            end if
        end do
        currentTime = (i) * del_t
        inelasticEnergyLoss = 0.0d0
        solver%particleEnergyLoss = 0.0d0
        call solveSingleTimeStepDiagnostic(solver, particleList, world, del_t, maxIter, eps_r)
        
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
        call loadParticleDensity(densities, particleList, world)
        call writeParticleDensity(densities, particleList, world, CurrentDiagStep, .false.) 
        call writePhi(solver%phi, CurrentDiagStep, .false.)
        call particleList(1)%writeLocalTemperature(CurrentDiagStep)
        do j=1, numberChargedParticles
            call particleList(j)%writePhaseSpace(CurrentDiagStep)
        end do
        write(22,"(6(es16.8,1x), (I4, 1x))") currentTime, inelasticEnergyLoss/del_t, &
                SUM(solver%particleChargeLoss)/del_t, solver%particleEnergyLoss/del_t, solver%chargeError, solver%energyError, solver%iterNumPicard
        close(22)

    end subroutine solveSimulationOnlyPotential

    subroutine solveSimulation(solver, particleList, world, del_t, maxIter, eps_r, irand, numTimeSteps, heatSkipSteps)
        ! Perform certain amount of timesteps, with diagnostics taken at first and last time step
        ! Impliment averaging for final select amount of timeSteps, this will be last data dump
        type(Particle), intent(in out) :: particleList(:)
        type(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t, eps_r
        integer(int32), intent(in) :: maxIter, numTimeSteps, heatSkipSteps
        integer(int32), intent(in out) :: irand
        integer(int32) :: numSkipSteps, i, j, CurrentDiagStep
        real(real64) :: currentTime, densities(NumberXNodes, numberChargedParticles)
        CurrentDiagStep = 1

        !Wrtie Initial conditions
        open(15,file='../Data/InitialConditions.dat')
        write(15,'("Number Grid Nodes, Final Expected Time(s), Delta t(s), FractionFreq")')
        write(15,"((I3.3, 1x), 3(es16.8,1x))") NumberXNodes, numTimeSteps*del_t, del_t, FractionFreq
        close(15)

        ! Write Particle properties
        open(9,file='../Data/ParticleProperties.dat')
        write(9,'("Particle Symbol, Particle Mass (kg), Particle Charge (C), Particle Weight (N/m^2)")')
        do j=1, numberChargedParticles
            write(9,"((A, 1x), 3(es16.8,1x))") particleList(j)%name, particleList(j)%mass, particleList(j)%q, particleList(j)%w_p
        end do
        close(9)

        numSkipSteps = numTimeSteps/(numDiagnosticSteps)
        101 format(20(1x,es16.8))
        open(22,file='../Data/GlobalDiagnosticData.dat')
        write(22,'("Time (s), Collision Loss (W/m^2), ParticleCurrentLoss (A/m^2), ParticlePowerLoss(W/m^2), chargeError (a.u), energyError(a.u), Picard Iteration Number")')
        
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

        do i = 1, numTimeSteps-1
            if (MOD((i-1), numSkipSteps) /= 0) then
                call solver%adaptiveSolveDivAmperePicard(particleList, world, del_t, maxIter, eps_r)
                call ionizationCollisionIsotropic(particleList(1), particleList(2), 1.0d20, 1.0d-20, del_t, 15.8d0, 0.025d0, irand)
            else  
                ! Data dump with Diagnostics
                print *, "Simulation is", real(i)/numTimeSteps * 100.0, "percent done"
                currentTime = (i) * del_t
                inelasticEnergyLoss = 0.0d0
                solver%particleEnergyLoss = 0.0d0
                call solveSingleTimeStepDiagnostic(solver, particleList, world, del_t, maxIter, eps_r)
                
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
                call ionizationCollisionIsotropic(particleList(1), particleList(2), 1.0d20, 1.0d-20, del_t, 15.8d0, 0.025d0, irand)
                densities = 0.0d0
                call loadParticleDensity(densities, particleList, world)
                call writeParticleDensity(densities, particleList, world, CurrentDiagStep, .false.) 
                call writePhi(solver%phi, CurrentDiagStep, .false.)
                call particleList(1)%writeLocalTemperature(CurrentDiagStep)
                do j=1, numberChargedParticles
                    call particleList(j)%writePhaseSpace(CurrentDiagStep)
                end do
                write(22,"(6(es16.8,1x), (I4, 1x))") currentTime, inelasticEnergyLoss/del_t, &
                SUM(solver%particleChargeLoss)/del_t, solver%particleEnergyLoss/del_t, solver%chargeError, solver%energyError, solver%iterNumPicard
                CurrentDiagStep = CurrentDiagStep + 1
            end if
            if (MOD(i, heatSkipSteps) == 0) then
                call addUniformPowerMaxwellian(particleList(1), Power, nu_h, irand, del_t, heatSkipSteps)
            end if
        end do

        currentTime = (i) * del_t
        inelasticEnergyLoss = 0.0d0
        solver%particleEnergyLoss = 0.0d0
        call solveSingleTimeStepDiagnostic(solver, particleList, world, del_t, maxIter, eps_r)
        
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
        call ionizationCollisionIsotropic(particleList(1), particleList(2), 1.0d20, 1.0d-20, del_t, 15.8d0, 0.0d0, irand)
        densities = 0.0d0
        call loadParticleDensity(densities, particleList, world)
        call writeParticleDensity(densities, particleList, world, CurrentDiagStep, .false.) 
        call writePhi(solver%phi, CurrentDiagStep, .false.)
        call particleList(1)%writeLocalTemperature(CurrentDiagStep)
        do j=1, numberChargedParticles
            call particleList(j)%writePhaseSpace(CurrentDiagStep)
        end do
        write(22,"(6(es16.8,1x), (I4, 1x))") currentTime, inelasticEnergyLoss/del_t, &
        SUM(solver%particleChargeLoss)/del_t, solver%particleEnergyLoss/del_t, solver%chargeError, solver%energyError, solver%iterNumPicard
        CurrentDiagStep = CurrentDiagStep + 1

    end subroutine solveSimulation


    subroutine solveSimulationFinalAverage(solver, particleList, world, del_t, maxIter, eps_r, irand, stepsAverage, heatSkipSteps)
        ! Perform certain amount of timesteps, with diagnostics taken at first and last time step
        ! Impliment averaging for final select amount of timeSteps, this will be last data dump
        type(Particle), intent(in out) :: particleList(:)
        type(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t, eps_r
        integer(int32), intent(in) :: maxIter, stepsAverage, heatSkipSteps
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
            call solver%adaptiveSolveDivAmperePicard(particleList, world, del_t, maxIter, eps_r)
            call ionizationCollisionIsotropic(particleList(1), particleList(2), 1.0d20, 1.0d-20, del_t, 15.8d0, 0.0d0, irand)
            call loadParticleDensity(densities, particleList, world)
            phi_average = phi_average + solver%phi
            if (MOD(i, heatSkipSteps) == 0) then
                call addUniformPowerMaxwellian(particleList(1), Power, nu_h, irand, del_t, heatSkipSteps)
            end if
        end do
        densities = densities/stepsAverage
        call writeParticleDensity(densities, particleList, world, 0, .true.) 
        call writePhi(phi_average/stepsAverage, 0, .true.)
        write(22,"((I4, 1x), 3(es16.8,1x))") stepsAverage, inelasticEnergyLoss/del_t/stepsAverage, SUM(solver%particleChargeLoss)/del_t/stepsAverage, solver%particleEnergyLoss/del_t/stepsAverage
        close(22)
        print *, "Electron average wall loss on left wall is:", solver%particleChargeLoss(1, 1)/del_t/stepsAverage
        print *, "Electron average wall loss on right wall is:", solver%particleChargeLoss(2, 1)/del_t/stepsAverage
        print *, "Ion average wall loss on left wall is:", solver%particleChargeLoss(1, 2)/del_t/stepsAverage
        print *, "Ion average wall loss on right wall is:", solver%particleChargeLoss(2, 2)/del_t/stepsAverage

    end subroutine solveSimulationFinalAverage




end module mod_simulation