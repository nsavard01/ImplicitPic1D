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
    ! Assume dirichlet boundaries at ends for now, so matrix solver can go down by two dimensions
    private
    public :: potSolver

    type :: potSolver
        real(real64), allocatable :: phi(:), J(:), rho(:)
        !real(real64) :: phi_left, phi_right
        real(real64), allocatable :: a_tri(:), b_tri(:), c_tri(:) !for thomas algorithm potential solver

    contains
        procedure, public, pass(self) :: depositRho

    end type

    interface potSolver
        module procedure :: potSolver_constructor
    end interface potSolver

contains

    type(potSolver) function potSolver_constructor(num_nodes) result(self)
        ! Construct domain object, initialize grid, dx_dl, and nodeVol.
        integer(int32), intent(in) :: num_nodes
        allocate(self % J(num_nodes-1), self % rho(num_nodes), self % phi(num_nodes), self%a_tri(num_nodes-3), &
        self%b_tri(num_nodes-2), self%c_tri(num_nodes-3))
        self % a_tri = 0
        self % b_tri = 0
        self % c_tri = 0
        self % rho = 0
        self % J = 0
        self % phi = 0

    end function potSolver_constructor

    subroutine depositRho(self, particleList, world) 
        class(potSolver), intent(in out) :: self
        type(Particle), intent(in) :: particleList(:)
        type(Domain), intent(in) :: world
        integer(int32) :: i, j, l_left
        real(real64) :: d
        self % rho = 0
        do i=1, size(particleList)
            do j = 1, particleList(i)%N_p
                l_left = INT(particleList(i)%l_p(j))
                d = MOD(particleList(i)%l_p(j), 1.0)
                self % rho(l_left) = self % rho(l_left) + particleList(i)%q * particleList(i)%w_p * (1.0-d)
                self % rho(l_left + 1) = self % rho(l_left + 1) + particleList(i)%q * particleList(i)%w_p * d
            end do
        end do
        self % rho = self % rho / world%nodeVol
    end subroutine depositRho





end module mod_potentialSolver