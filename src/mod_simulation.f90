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

    integer(int32) :: maxIter = 50, numDiagnosticSteps = 8, numTimeSteps = 24
    real(real64) :: eps_r = 1e-8, del_t, fractionFreq = 0.5d0

contains

    ! -------------------------------- Diagnostic ------------------------------------
    subroutine WriteParticleDensity(solver, particleList, world, CurrentDiagStep) 
        ! For diagnostics, deposit single particle density
        ! Re-use rho array since it doesn't get used after first Poisson
        type(potSolver), intent(in out) :: solver
        type(Particle), intent(in) :: particleList(:)
        type(Domain), intent(in) :: world
        integer(int32), intent(in) :: CurrentDiagStep
        integer(int32) :: i, j, l_left
        character(len=5) :: char_i
        real(real64) :: d
        
        do i=1, size(particleList)
            solver % rho = 0.0d0
            do j = 1, particleList(i)%N_p
                l_left = INT(particleList(i)%l_p(j))
                d = MOD(particleList(i)%l_p(j), 1.0d0)
                solver % rho(l_left) = solver % rho(l_left) + particleList(i)%w_p * (1.0d0-d)
                solver % rho(l_left + 1) = solver % rho(l_left + 1) + particleList(i)%w_p * d
            end do
            solver % rho = solver % rho / world%nodeVol
            write(char_i, '(I3)'), CurrentDiagStep
            open(41,file='../Data/density_'//particleList(i)%name//"_"//trim(adjustl(char_i))//".dat", form='UNFORMATTED')
            write(41) solver%rho
            close(41)
        end do
        
    end subroutine WriteParticleDensity




    ! -------------------------- Simulation ------------------------------------------


    subroutine solveInitialPotential(particleList, solver, world)
        ! Solve for initial potential
        type(Particle), intent(in) :: particleList(:)
        type(Domain), intent(in) :: world
        type(potSolver), intent(in out) :: solver
        call solver%depositRho(particleList, world)
        call solver%solve_tridiag_Poisson()
        ! Assume only use potential solver once, then need to generate matrix for Div-Ampere
        call solver%construct_diagMatrix_Ampere(world)

    end subroutine solveInitialPotential

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

        call solver%solveDivAmperePicard(particleList, world, del_t, maxIter, eps_r, boolDiagnostic)

        call addUniformPowerMaxwellianNicolas(particleList(1), Power, 0.05d0, irand, del_t)
        ! P_before = m_e * SUM(particleList(1)%v_p, DIM = 1) + particleList(2)%mass *SUM(particleList(2)%v_p, DIM = 1)
        ! E_before = m_e * SUM(particleList(1)%v_p**2) * 0.5d0 / e + particleList(2)%mass * SUM(particleList(2)%v_p**2) * 0.5d0/e
        call ionizationCollisionIsotropic(particleList(1), particleList(2), 1.0d20, 1.0d-20, del_t, 15.8d0, irand, 0.0d0)
        ! P_after = m_e * SUM(particleList(1)%v_p, DIM = 1) + particleList(2)%mass * SUM(particleList(2)%v_p, DIM = 1)
        ! E_after = m_e * SUM(particleList(1)%v_p**2) * 0.5d0 / e + particleList(2)%mass * SUM(particleList(2)%v_p**2) * 0.5d0/e + inelasticEnergyLoss
        ! print *, "P_before is:", P_before
        ! print *, "P_after is:", P_after
        ! print *, "E_before is:", E_before
        ! print *, "E_after is:", E_after


    end subroutine solveSingleTimeStep


    subroutine solveSimulation(solver, particleList, world, del_t, maxIter, eps_r, irand, numTimeSteps)
        ! Perform certain amount of timesteps, with diagnostics taken at first and last time step
        type(Particle), intent(in out) :: particleList(:)
        type(potSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t, eps_r
        integer(int32), intent(in) :: maxIter, numTimeSteps
        integer(int32), intent(in out) :: irand
        integer(int32) :: numSkipSteps, i, CurrentDiagStep = 1
        real(real64) :: currentTime
        open(15,file='../Data/InitialConditions.dat',access='APPEND')
        write(15,'("Number Grid Nodes, Final Expected Time(s)")')
        write(15,"((I3.3, 1x), (es16.8,1x))") n_x, numTimeSteps*del_t
        close(15)
        close(22)
        numSkipSteps = numTimeSteps/(numDiagnosticSteps)
        101 format(20(1x,es16.8))
        open(22,file='../Data/ScalarDiagnosticData.dat',access='APPEND')
        write(22,'("Time (s), InelasticEnergyLoss (eV), ParticleCurrentLoss (A/m^2), ParticlePowerLoss(W/m^2)")')

        !Save initial particle/field data, along with domain
        call WriteParticleDensity(solver, particleList, world, 0) 
        call world%writeDomain()

        do i = 1, numTimeSteps-1
            print *, "Current Step num is:", i
            if (MOD((i-1), numSkipSteps) /= 0) then
                call solveSingleTimeStep(solver, particleList, world, del_t, maxIter, eps_r, irand, .true.)
            else
                currentTime = (i) * del_t
                call solveSingleTimeStep(solver, particleList, world, del_t, maxIter, eps_r, irand, .true.)
                call WriteParticleDensity(solver, particleList, world, CurrentDiagStep) 
                write(22,101) currentTime, inelasticEnergyLoss, solver%particleCurrentLoss, solver%particlePowerLoss
                CurrentDiagStep = CurrentDiagStep + 1
            end if
        end do
        currentTime = (numTimeSteps) * del_t
        call solveSingleTimeStep(solver, particleList, world, del_t, maxIter, eps_r, irand, .true.)
        call WriteParticleDensity(solver, particleList, world, CurrentDiagStep) 
        write(22,101) currentTime, inelasticEnergyLoss, solver%particleCurrentLoss, solver%particlePowerLoss
        close(22)






    end subroutine solveSimulation







end module mod_simulation