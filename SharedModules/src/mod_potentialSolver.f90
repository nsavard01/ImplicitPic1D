module mod_potentialSolver
    use iso_fortran_env, only: int32, real64, output_unit
    use constants
    use mod_BasicFunctions
    use mod_particle
    use mod_domain
    implicit none

    ! This module will be for the potential solver
    ! This uses the interaction between particles (particle module) and grid (domain module) in order to solve for the potential
    ! Will contain particle to mesh gathers (J, rho) and potential which comes from them (phi)
    ! Will also contain particle mover, since needed to store to J, and cannot be separate
    ! Assume dirichlet boundaries at ends for now, so matrix solver can go down by two dimensions
    private
    public :: potentialSolver

    type :: potentialSolver
        real(real64), allocatable :: phi(:), J(:), rho(:), phi_f(:), particleChargeLoss(:,:) !phi_f is final phi, will likely need to store two arrays for phi, can't be avoided
        real(real64) :: energyError, chargeError, particleEnergyLoss
        real(real64), allocatable :: a_tri(:), b_tri(:), c_tri(:) !for thomas algorithm potential solver, a_tri is lower diagonal, b_tri middle, c_tri upper


    contains
        procedure, public, pass(self) :: depositRho
        procedure, public, pass(self) :: solve_tridiag_Poisson
        procedure, public, pass(self) :: solve_tridiag_Ampere
        procedure, public, pass(self) :: getTotalPE
        procedure, public, pass(self) :: getError_tridiag_Ampere
        procedure, public, pass(self) :: getError_tridiag_Poisson
        procedure, public, pass(self) :: construct_diagMatrix_Ampere
        procedure, public, pass(self) :: solveInitialPotential
        procedure, public, pass(self) :: construct_diagMatrix
    end type

    interface potentialSolver
        module procedure :: potentialSolver_constructor
    end interface potentialSolver

contains

    type(potentialSolver) function potentialSolver_constructor(world, leftVoltage, rightVoltage) result(self)
        ! Construct domain object, initialize grid, dx_dl, and nodeVol.
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: leftVoltage, rightVoltage
        allocate(self % J(NumberXNodes-1), self % rho(NumberXNodes), self % phi(NumberXNodes), self % phi_f(NumberXNodes), self%a_tri(NumberXNodes-1), &
        self%b_tri(NumberXNodes), self%c_tri(NumberXNodes-1), self%particleChargeLoss(2, numberChargedParticles))
        call construct_diagMatrix(self, world)
        self % rho = 0
        self % J = 0
        self % phi = 0
        ! self%coeff_left = 0.0d0
        ! self%coeff_left = 0.0d0
        self%particleEnergyLoss = 0.0d0
        self%particleChargeLoss = 0.0d0
        self%energyError = 0.0d0
        self%chargeError = 0.0d0
        SELECT CASE (world%boundaryConditions(1))
        CASE(1)
            self%phi(1) = leftVoltage
        CASE(2)
            self%phi(1) = 0.0d0
        CASE(3)
            self%phi(1) = leftVoltage
            self%phi(NumberXNodes) = leftVoltage
        END SELECT

        SELECT CASE (world%boundaryConditions(NumberXNodes))
        CASE(1,3)
            self%phi(NumberXNodes) = rightVoltage
        CASE(2)
            self%phi(NumberXNodes) = 0.0d0
        END SELECT
        self%phi_f = self%phi  

    end function potentialSolver_constructor

    subroutine depositRho(self, particleList, world) 
        class(potentialSolver), intent(in out) :: self
        type(Particle), intent(in) :: particleList(:)
        type(Domain), intent(in) :: world
        integer(int32) :: i, j, l_left
        real(real64) :: d
        self % rho = 0.0d0
        do i=1, numberChargedParticles
            do j = 1, particleList(i)%N_p
                l_left = INT(particleList(i)%phaseSpace(1, j))
                d = MOD(particleList(i)%phaseSpace(1, j), 1.0d0)
                self % rho(l_left) = self % rho(l_left) + particleList(i)%q * particleList(i)%w_p * (1.0d0-d)
                self % rho(l_left + 1) = self % rho(l_left + 1) + particleList(i)%q * particleList(i)%w_p * d
            end do
        end do
        self % rho = self % rho / world%nodeVol
        if (world%boundaryConditions(1) == 3) then
            self%rho(1) = self%rho(1) + self%rho(NumberXNodes)
            self%rho(NumberXNodes) = self%rho(1)
        end if
    end subroutine depositRho

    subroutine construct_diagMatrix(self, world)
        ! construct diagonal components for thomas algorithm
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        integer(int32) :: i
        do i = 1, NumberXNodes
            SELECT CASE (world%boundaryConditions(i))
            CASE(0)
                if (i < NumberXNodes) then
                    self % c_tri(i) = 1.0d0/world%nodeVol(i)/world%dx_dl(i)
                end if
                if (i > 1) then
                    self%a_tri(i-1) = 1.0d0/world%nodeVol(i)/ world%dx_dl(i-1)
                end if
                self%b_tri(i) = - (1.0d0/world%dx_dl(i-1)  + 1.0d0/world%dx_dl(i)) /world%nodeVol(i)
            CASE(1)
                self%b_tri(i) = 1.0d0
            CASE(2)
                if (i == 1) then
                    self % c_tri(i) = 2.0d0/world%nodeVol(i)/world%dx_dl(i)
                    !self%a_tri(i - leftNodeIdx) = 2.0d0/(world%dx_dl(i-1) + world%dx_dl(i)) / world%dx_dl(i-1)
                    self%b_tri(i) = - (2.0d0/world%dx_dl(i)) /world%nodeVol(i)
                else if (i == NumberXNodes) then
                    self % a_tri(i-1) = 2.0d0/world%nodeVol(i)/ world%dx_dl(i-1)
                    !self % c_tri(i - leftNodeIdx) = 2.0d0/(world%dx_dl(i-2) + world%dx_dl(i-1))/world%dx_dl(i-1)
                    self%b_tri(i) = - (2.0d0/world%dx_dl(i-1)) /world%nodeVol(i)
                else
                    print *, "Neumann boundary not on left or right most index!"
                    stop
                end if
            CASE(3)
                self%b_tri(i) = 1.0d0
            CASE default
                print *, "Error when constructing poisson matrix, inner nodes not plasma or neumann!"
            END SELECT
        end do

    end subroutine construct_diagMatrix

    subroutine construct_diagMatrix_Ampere(self, world)
        ! construct diagonal components for thomas algorithm, for Ampere (after initial Poisson)
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        self % a_tri = -1.0d0/world%dx_dl(2:NumberXNodes-2)
        self % c_tri = -1.0d0/world%dx_dl(2:NumberXNodes-2)
        self % b_tri = (world%dx_dl(1:NumberXNodes-2) + world%dx_dl(2:NumberXNodes-1))/ (world%dx_dl(1:NumberXNodes-2) * world%dx_dl(2:NumberXNodes-1))

    end subroutine construct_diagMatrix_Ampere

    subroutine solve_tridiag_Poisson(self, world)
        ! Tridiagonal (Thomas algorithm) solver for initial Poisson
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        integer(int32) :: i
        real(real64) :: m, d(NumberXNodes), cp(NumberXNodes-1),dp(NumberXNodes)

        do i=1, NumberXNodes
            SELECT CASE (world%boundaryConditions(i))
            CASE(0,2)
                d(i) = -self%rho(i) / eps_0
            CASE(1,3)
                d(i) = self%phi(i)
            END SELECT
        end do
    ! initialize c-prime and d-prime
        cp(1) = self%c_tri(1)/self%b_tri(1)
        dp(1) = d(1)/self%b_tri(1)
    ! solve for vectors c-prime and d-prime
        do i = 2,NumberXNodes-1
            m = self%b_tri(i)-cp(i-1)*self%a_tri(i-1)
            cp(i) = self%c_tri(i)/m
            dp(i) = (d(i)-dp(i-1)*self%a_tri(i-1))/m
        end do
        m = self%b_tri(NumberXNodes)-cp(NumberXNodes-1)*self%a_tri(NumberXNodes-1)
        dp(NumberXNodes) = (d(NumberXNodes)-dp(NumberXNodes-1)*self%a_tri(NumberXNodes-1))/m
    ! initialize x
        self%phi = dp
    ! solve for x from the vectors c-prime and d-prime
        do i = NumberXNodes-1, 1, -1
            self%phi(i) = dp(i)-cp(i)*self%phi(i+1)
        end do
        self%phi_f = self%phi

    end subroutine solve_tridiag_Poisson

    function getError_tridiag_Poisson(self) result(res)
        class(potentialSolver), intent(in) :: self
        integer(int32) :: i
        real(real64) :: Ax(NumberXNodes-2), d(NumberXNodes-2), res
        Ax(1) = self%b_tri(1)*self%phi_f(2) + self%c_tri(1) * self%phi_f(3)
        do i=2, NumberXNodes-3
            Ax(i) = self%b_tri(i)*self%phi_f(i+1) + self%c_tri(i) * self%phi_f(i+2) + self%a_tri(i-1) * self%phi_f(i)
        end do
        Ax(NumberXNodes-2) = self%b_tri(NumberXNodes-2)*self%phi_f(NumberXNodes-1) + self%a_tri(NumberXNodes-3) * self%phi_f(NumberXNodes-2)
        d = -self%rho(2:NumberXNodes-1)
        !res = Ax*eps_0 - d
        res = SQRT(SUM((Ax*eps_0 - d)**2))

    end function getError_tridiag_Poisson

    subroutine solve_tridiag_Ampere(self, world, del_t)
        ! Tridiagonal (Thomas algorithm) solver for Ampere
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t
        integer(int32) :: i, n !n size dependent on how many points are inside (not boundary), so how many in rhs equation
        real(real64) :: m, d(NumberXNodes - 2), cp(NumberXNodes - 3),dp(NumberXNodes- 2)
        n = NumberXNodes - 2

        d = (-self%J(2:) + self%J(1:n)) * del_t / eps_0 + arrayDiff(self%phi(1:n+1), n+1)/world%dx_dl(1:n) - arrayDiff(self%phi(2:), n+1)/world%dx_dl(2:)
    ! initialize c-prime and d-prime
        cp(1) = self%c_tri(1)/self%b_tri(1)
        dp(1) = d(1)/self%b_tri(1)
    ! solve for vectors c-prime and d-prime
        do i = 2,n-1
            m = self%b_tri(i)-cp(i-1)*self%a_tri(i-1)
            cp(i) = self%c_tri(i)/m
            dp(i) = (d(i)-dp(i-1)*self%a_tri(i-1))/m
        end do
        dp(n) = (d(n)-dp(n-1)*self%a_tri(n-1))/(self%b_tri(n)-cp(n-1)*self%a_tri(n-1))
    ! initialize phi_f
        self%phi_f(2:n+1) = dp(n)
    ! solve for x from the vectors c-prime and d-prime
        do i = n-1, 1, -1
            self%phi_f(i+1) = dp(i)-cp(i)*self%phi_f(i+2)
        end do

    end subroutine solve_tridiag_Ampere

    function getError_tridiag_Ampere(self, world, del_t) result(res)
        class(potentialSolver), intent(in) :: self
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t
        integer(int32) :: i
        real(real64) :: Ax(NumberXNodes-2), d(NumberXNodes-2), res(NumberXNodes-2)
        Ax(1) = self%b_tri(1)*self%phi_f(2) + self%c_tri(1) * self%phi_f(3)
        do i=2, NumberXNodes-3
            Ax(i) = self%b_tri(i)*self%phi_f(i+1) + self%c_tri(i) * self%phi_f(i+2) + self%a_tri(i-1) * self%phi_f(i)
        end do
        Ax(NumberXNodes-2) = self%b_tri(NumberXNodes-2)*self%phi_f(NumberXNodes-1) + self%a_tri(NumberXNodes-3) * self%phi_f(NumberXNodes-2)
        d = (-self%J(2:) + self%J(1:NumberXNodes-2)) * del_t / eps_0 + arrayDiff(self%phi(1:NumberXNodes-1), NumberXNodes-1)/world%dx_dl(1:NumberXNodes-2) - arrayDiff(self%phi(2:), NumberXNodes-1)/world%dx_dl(2:)
        !res = SQRT(SUM(((Ax- d)/self%minEField)**2)/(NumberXNodes-2))
        res = Ax- d

    end function getError_tridiag_Ampere

    function getTotalPE(self, world, future) result(res)
        ! Get energy in electric fields, future true, then derive from phi_f, otherwise phi
        ! In 1D is J/m^2
        class(potentialSolver), intent(in) :: self
        type(Domain), intent(in) :: world
        logical :: future
        real(real64) :: res
        if (future) then
            res = 0.5 * eps_0 * SUM(arrayDiff(self%phi_f, NumberXNodes)**2 / world%dx_dl)
        else
            res = 0.5 * eps_0 * SUM(arrayDiff(self%phi, NumberXNodes)**2 / world%dx_dl)
        end if
    end function getTotalPE

    ! ---------------- Initial Poisson Solver -------------------------------------------------

    subroutine solveInitialPotential(self, particleList, world)
        ! Solve for initial potential
        class(potentialSolver), intent(in out) :: self
        type(Particle), intent(in) :: particleList(:)
        type(Domain), intent(in) :: world
        call self%depositRho(particleList, world)
        call self%solve_tridiag_Poisson(world)
        ! Assume only use potential solver once, then need to generate matrix for Div-Ampere
        call self%construct_diagMatrix_Ampere(world)

    end subroutine solveInitialPotential



end module mod_potentialSolver