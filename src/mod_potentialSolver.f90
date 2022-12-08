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
        real(real64), allocatable :: phi(:), J(:), rho(:), phi_f(:) !phi_f is final phi, will likely need to store two arrays for phi, can't be avoided
        real(real64) :: phi_left, phi_right
        real(real64) :: coeff_left, coeff_right ! these are coefficients (from world dimensions) needed with phi_left and phi_right in rhs of matrix equation
        real(real64), allocatable :: a_tri(:), b_tri(:), c_tri(:) !for thomas algorithm potential solver, a_tri is lower diagonal, b_tri middle, c_tri upper

    contains
        procedure, public, pass(self) :: depositRho
        procedure, public, pass(self) :: solve_tridiag_Poisson
        procedure, public, pass(self) :: getEField
        procedure, private, pass(self) :: construct_diagMatrix
    end type

    interface potSolver
        module procedure :: potSolver_constructor
    end interface potSolver

contains

    type(potSolver) function potSolver_constructor(num_nodes, world) result(self)
        ! Construct domain object, initialize grid, dx_dl, and nodeVol.
        integer(int32), intent(in) :: num_nodes
        type(Domain), intent(in) :: world
        allocate(self % J(num_nodes-1), self % rho(num_nodes), self % phi(num_nodes), self % phi_f(num_nodes), self%a_tri(num_nodes-3), &
        self%b_tri(num_nodes-2), self%c_tri(num_nodes-3))
        call construct_diagMatrix(self, world)
        self % rho = 0
        self % J = 0
        self % phi = 0
        self % phi_left = 0
        self % phi_right = 0
        self % coeff_left = 2/(world%dx_dl(1) + world%dx_dl(2))/world%dx_dl(1)
        self % coeff_right = 2/(world%dx_dl(size(world%dx_dl)-1) + world%dx_dl(size(world%dx_dl)))/world%dx_dl(size(world%dx_dl))

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

    subroutine construct_diagMatrix(self, world)
        ! construct diagonal components for thomas algorithm
        class(potSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        self % a_tri = 2/(world%dx_dl(2:world%n_x-2) + world%dx_dl(3:)) / world%dx_dl(2:world%n_x-2)
        self % c_tri = 2/(world%dx_dl(1:world%n_x-3) + world%dx_dl(2:world%n_x-2))/world%dx_dl(2:world%n_x-2)
        self % b_tri = -2/(world%dx_dl(1:world%n_x-2) + world%dx_dl(2:))/world%dx_dl(1:world%n_x-2) - 2/(world%dx_dl(1:world%n_x-2) + world%dx_dl(2:))/world%dx_dl(2:world%n_x-1)

    end subroutine construct_diagMatrix

    subroutine solve_tridiag_Poisson(self)
        ! Tridiagonal (Thomas algorithm) solver for initial Poisson
        class(potSolver), intent(in out) :: self
        integer(int32) :: i, n !n size dependent on how many points are boundary conditions and thus moved to rhs equation
        real(real64) :: m, d(size(self%phi) - 2), cp(size(self%phi) - 2),dp(size(self%phi) - 2)
        n = size(self%phi) - 2

        d = -self%rho / eps_0
        d(1) = d(1) - 2 * self%phi_left * self%coeff_left
        d(n) = d(n) - 2 * self%phi_right * self%coeff_right
    ! initialize c-prime and d-prime
        cp(1) = self%c_tri(1)/self%b_tri(1)
        dp(1) = d(1)/self%b_tri(1)
    ! solve for vectors c-prime and d-prime
        do i = 2,n
            m = self%b_tri(i)-cp(i-1)*self%a_tri(i-1)
            cp(i) = self%c_tri(i)/m
            dp(i) = (d(i)-dp(i-1)*self%a_tri(i-1))/m
        end do
    ! initialize x
        self%phi(2:n+1) = dp(n)
    ! solve for x from the vectors c-prime and d-prime
        do i = n-1, 1, -1
            self%phi(i+1) = dp(i)-cp(i)*self%phi(i+2)
        end do
        self%phi_f = self%phi

    end subroutine solve_tridiag_Poisson

    pure function getEField(self, l_p, world) result(EField)
        !return EField of particle at logical position l_p (this is half distance), per particle since particle mover loop will be per particle
        class(potSolver), intent(in) :: self
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: l_p
        real(real64) :: EField
        EField = (self%phi_f(INT(l_p)) + self%phi(INT(l_p)) - self%phi(INT(l_p)+1) - self%phi_f(INT(l_p) + 1)) / world%dx_dl(INT(l_p))/2
    end function getEField





end module mod_potentialSolver