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

    integer(int32) :: maxIter = 50, numDiagnosticSteps = 4, numTimeSteps = 10
    real(real64) :: eps_r = 1e-8, del_t, fractionFreq = 0.5d0
    real(real64), allocatable :: particleDensity(:,:,:), electricPotential(:, :), simulationTime(:)

contains

    ! -------------------------------- Diagnostic ------------------------------------
    subroutine depositSingleParticleDensity(particleDensity, particleList, world) 
        ! For diagnostics, deposit single particle density
        ! Re-use rho array since it doesn't get used after first Poisson
        type(Particle), intent(in) :: particleList(:)
        type(Domain), intent(in) :: world
        real(real64), intent(in out) :: particleDensity(:,:) !Input Density array without time dimension
        integer(int32) :: i, j, l_left
        real(real64) :: d
        particleDensity = 0.0d0
        do i=1, size(particleList)
            do j = 1, particleList(i)%N_p
                l_left = INT(particleList(i)%l_p(j))
                d = MOD(particleList(i)%l_p(j), 1.0d0)
                particleDensity(i, l_left) = particleDensity(i, l_left) + particleList(i)%w_p * (1.0d0-d)
                particleDensity(i, l_left + 1) = particleDensity(i, l_left + 1) + particleList(i)%w_p * d
            end do
            particleDensity(i, :) = particleDensity(i, :)/world%nodeVol
        end do
    end subroutine depositSingleParticleDensity




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


    ! subroutine solveSimulation(solver, particleList, world, del_t, maxIter, eps_r, irand, numTimeSteps)
    !     type(Particle), intent(in out) :: particleList(:)
    !     type(potSolver), intent(in out) :: solver
    !     type(Domain), intent(in) :: world
    !     real(real64), intent(in) :: del_t, eps_r
    !     logical, intent(in) :: boolDiagnostic
    !     integer(int32), intent(in) :: maxIter, numTimeSteps
    !     integer(int32), intent(in out) :: irand
    !     integer(int32) :: numSkipSteps, currentStepNumber = 1, i
    !     allocate(numDiagnosticSteps, particleDensity(size(particleList), n_x), electricPotential(numDiagnosticSteps, n_x), simulationTime(numDiagnosticSteps))
    !     numSkipSteps = numTimeSteps/(numDiagnosticSteps-1)
    !     do i = 1, numTimeSteps
    !         if (MOD((i-1), numSkipSteps) /= 0) then
    !             call solveSingleTimeStep(solver, particleList, world, del_t, maxIter, eps_r, irand, .false.)
    !         else
    !             call solver%depositSingleParticleDensity(particleList(1), world)

    !         end if
    !     end do



    ! end subroutine solveSimulation







end module mod_simulation