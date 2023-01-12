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

    integer(int32) :: maxIter = 50, numRunSteps = 100, numTimeSteps = 0
    real(real64) :: eps_r = 1e-8, del_t, fractionFreq = 0.5d0
    real(real64), allocatable :: electronDensity(:,:), electricPotential(:)

contains
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

    subroutine solveSingleTimeStep(solver, particleList, world, del_t, maxIter, eps_r, irand)
        ! Single time step solver with Divergence of ampere, followed by adding of power, followed by collisions
        type(Particle), intent(in out) :: particleList(:)
        type(potSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t, eps_r
        integer(int32), intent(in) :: maxIter
        integer(int32), intent(in out) :: irand

        call solver%solveDivAmperePicard(particleList, world, del_t, maxIter, eps_r, .true.)

        call addUniformPowerMaxwellianNicolas(particleList(1), Power, 0.05d0, irand, del_t)
        call ionizationCollisionIsotropic(particleList(1), particleList(2), 1.0d20, 1.0d-20, del_t, 15.8d0, irand, 300.0d0 * k_B/e)


    end subroutine solveSingleTimeStep





end module mod_simulation