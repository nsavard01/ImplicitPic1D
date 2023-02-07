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

    subroutine readInputs(n_x, maxIter, numDiagnosticSteps, stepsAverage, eps_r, fractionFreq, n_ave, T_e, T_i, L_domain, del_l, numParticles, particleIdxFactor)
        ! Set initial conditions and global constants based on read input from txt file
        integer(int32), intent(in out) :: n_x, maxIter, numDiagnosticSteps, stepsAverage, numParticles, particleIdxFactor
        real(real64), intent(in out) :: eps_r, fractionFreq, n_ave, T_e, T_i, L_domain, del_l
        integer(int32) :: io
        open(10,file='../InputData/InitialConditions.inp', IOSTAT=io)
        read(10, *, IOSTAT = io) T_e
        read(10, *, IOSTAT = io) T_i
        read(10, *, IOSTAT = io) numParticles
        read(10, *, IOSTAT = io) eps_r
        read(10, *, IOSTAT = io) particleIdxFactor
        read(10, *, IOSTAT = io) n_ave
        read(10, *, IOSTAT = io) numDiagnosticSteps
        read(10, *, IOSTAT = io) maxIter
        read(10, *, IOSTAT = io) fractionFreq
        read(10, *, IOSTAT = io) stepsAverage
        close(10)

        open(10,file='../InputData/Geometry.inp')
        read(10, *, IOSTAT = io) n_x
        read(10, *, IOSTAT = io) L_domain
        read(10, *, IOSTAT = io) del_l
        close(10)
        
    end subroutine readInputs

    function readParticleInputs(numberChargedParticles) result(particleList)
        type(Particle), allocatable :: particleList(:)
        integer(int32), intent(in out) :: numberChargedParticles
        integer(int32) :: j, numSpecies = 0
        character(len=6) :: name, particleNames(100)
        real(real64) :: mass(100), charge(100), Ti(100)
        numSpecies = 1 ! option in future to not have electrons by default
        mass(1) = m_e
        charge(1) = -1.0d0
        Ti(1) = T_e
        particleNames(1) = '[e]'
        print *, "Reading particle inputs"
        open(10,file='../InputData/BoundExample.dat')

        do j=1, 10000
            read(10,*,END=101,ERR=100) name

            if( name(1:4).eq.'IONS' .or. name(1:4).eq.'Ions' .or. name(1:4).eq.'ions' ) then
                goto 20
            endif
            ! Take care of extra text I guess
20          do while(name.ne.'------')
                read(10,*,END=101,ERR=100) name
            end do
            
200         read(10,'(A6)',END=101,ERR=100, ADVANCE = 'NO') name
            if (name.eq.'------') then
                close(10)
            else
                numSpecies = numSpecies + 1
                read(10,*,END=101,ERR=100) mass(numSpecies),charge(numSpecies), Ti(numSpecies)
                mass(numSpecies) = mass(numSpecies) * m_p
                particleNames(numSpecies) = trim(name)
                goto 200
            end if

        end do
100     continue
101     continue
        numberChargedParticles = numSpecies
        print *, numberChargedParticles
        allocate(particleList(numberChargedParticles))
        do j=1, numberChargedParticles
            particleList(j) = Particle(mass(j), e * charge(j), 1.0d0, numParticles, numParticles * particleIdxFactor, trim(particleNames(j)))
        end do
        



    end function readParticleInputs

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
        type(potSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t, eps_r
        logical, intent(in) :: boolDiagnostic
        integer(int32), intent(in) :: maxIter
        integer(int32), intent(in out) :: irand
        !real(real64) :: P_before(3), P_after(3), E_before, E_after

        call solver%adaptiveSolveDivAmperePicard(particleList, world, del_t, maxIter, eps_r, boolDiagnostic)

        !call addUniformPowerMaxwellianNicolas(particleList(1), Power, 0.05d0, irand, del_t)
        call addUniformPowerMaxwellian(particleList(1), Power, nu_h, irand, del_t)
        ! P_before = m_e * SUM(particleList(1)%phaseSpace(2:4,:), DIM = 2) + particleList(2)%mass *SUM(particleList(2)%phaseSpace(2:4,:), DIM = 2)
        ! E_before = m_e * SUM(particleList(1)%phaseSpace(2:4,:)**2) * 0.5d0 / e + particleList(2)%mass * SUM(particleList(2)%phaseSpace(2:4, :)**2) * 0.5d0/e
        call ionizationCollisionIsotropic(particleList(1), particleList(2), 1.0d20, 1.0d-20, del_t, 15.8d0, irand)
        ! P_after = m_e * SUM(particleList(1)%phaseSpace(2:4,:), DIM = 2) + particleList(2)%mass *SUM(particleList(2)%phaseSpace(2:4,:), DIM = 2)
        ! E_after = m_e * SUM(particleList(1)%phaseSpace(2:4,:)**2) * 0.5d0 / e + particleList(2)%mass * SUM(particleList(2)%phaseSpace(2:4, :)**2) * 0.5d0/e + inelasticEnergyLoss
        ! print *, "P_before is:", P_before
        ! print *, "P_after is:", P_after
        ! print *, "E_before is:", E_before
        ! print *, "E_after is:", E_after
        ! print *, "Inelastic Energy loss is:", inelasticEnergyLoss
        ! stop "Just to check collision"


    end subroutine solveSingleTimeStep


    subroutine solveSimulation(solver, particleList, world, del_t, maxIter, eps_r, irand, numTimeSteps, stepsAverage)
        ! Perform certain amount of timesteps, with diagnostics taken at first and last time step
        ! Impliment averaging for final select amount of timeSteps, this will be last data dump
        type(Particle), intent(in out) :: particleList(:)
        type(potSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t, eps_r
        integer(int32), intent(in) :: maxIter, numTimeSteps, stepsAverage
        integer(int32), intent(in out) :: irand
        integer(int32) :: numSkipSteps, i, j,CurrentDiagStep
        real(real64) :: currentTime, phi_average(n_x), densities(n_x, numberChargedParticles)
        CurrentDiagStep = 1

        open(15,file='../Data/InitialConditions.dat',access='APPEND')
        write(15,'("Number Grid Nodes, Final Expected Time(s), Delta t(s)")')
        write(15,"((I3.3, 1x), 2(es16.8,1x))") n_x, numTimeSteps*del_t, del_t
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
        do j=1, size(particleList)
            call particleList(j)%writePhaseSpace(0)
        end do

        do i = 1, numTimeSteps-stepsAverage
            if (MOD((i-1), numSkipSteps) /= 0) then
                call solveSingleTimeStep(solver, particleList, world, del_t, maxIter, eps_r, irand, .false.)
            else  
                print *, "Simulation is", real(i)/numTimeSteps * 100.0, "percent done"
                currentTime = (i) * del_t
                solver%particleEnergyLoss = 0.0d0
                solver%particleChargeLoss = 0.0d0
                inelasticEnergyLoss = 0.0d0
                call solveSingleTimeStep(solver, particleList, world, del_t, maxIter, eps_r, irand, .true.)
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
        do i =1, stepsAverage
            call solveSingleTimeStep(solver, particleList, world, del_t, maxIter, eps_r, irand, .true.)
            call loadParticleDensity(densities, particleList)
            phi_average = phi_average + solver%phi
        end do
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