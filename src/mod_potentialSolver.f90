module mod_potentialSolver
    use iso_fortran_env, only: int32, real64
    use constants
    use mod_particle
    use mod_domain
    implicit none

    ! This module will be for the potential solver
    ! This uses the interaction between particles (particle module) and grid (domain module) in order to solve for the potential
    ! Will contain particle to mesh gathers (J, rho) and potential which comes from them (phi)
    ! Will also contain particle mover, since needed to store to J, and cannot be separate
    private
    public :: potSolver

    type :: potSolver
        real(real64), allocatable :: phi(:), J(:), rho(:)
        type(Domain) :: world
        type(Particle), allocatable :: particleList(:)

    end type

    interface potSolver
        module procedure :: potSolver_constructor
    end interface potSolver

contains

    type(potSolver) function potSolver_constructor(particleList, world) result(self)
        ! Construct domain object, initialize grid, dx_dl, and nodeVol.
        type(Domain), intent(in) :: world
        type(Particle), intent(in) :: particleList(:)
        allocate(self % J(world%n_x-1), self % rho(world%n_x), self % phi(world%n_x))
        self % world = world
        self % particleList = particleList
    end function potSolver_constructor


end module mod_potentialSolver