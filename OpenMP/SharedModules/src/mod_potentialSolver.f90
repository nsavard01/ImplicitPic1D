module mod_potentialSolver
    use iso_fortran_env, only: int32, real64, output_unit
    use constants
    use mod_BasicFunctions
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
        real(real64), allocatable :: phi(:), J(:, :), rho(:), phi_f(:) !phi_f is final phi, will likely need to store two arrays for phi, can't be avoided
        real(real64) :: rho_const
        real(real64), allocatable :: a_tri(:), b_tri(:), c_tri(:) !for thomas algorithm potential solver, a_tri is lower diagonal, b_tri middle, c_tri upper


    contains
        procedure, public, pass(self) :: solve_tridiag_Poisson
        procedure, public, pass(self) :: solve_tridiag_Ampere
        procedure, public, pass(self) :: getTotalPE
        procedure, public, pass(self) :: getError_tridiag_Ampere
        procedure, public, pass(self) :: solve_CG_Ampere
        procedure, public, pass(self) :: getError_tridiag_Poisson
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
        allocate(self % J(NumberXNodes-1, numThread), self%rho(NumberXNodes), self % phi(NumberXNodes), self % phi_f(NumberXNodes), self%a_tri(NumberXNodes-1), &
        self%b_tri(NumberXNodes), self%c_tri(NumberXNodes-1))
        call construct_diagMatrix(self, world)
        self % rho = 0.0d0
        self % J = 0.0d0
        self % phi = 0.0d0
        self % rho_const = 0.0d0
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

    subroutine construct_diagMatrix(self, world)
        ! construct diagonal components for thomas algorithm
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        integer(int32) :: i
        do i = 1, NumberXNodes
            SELECT CASE (world%boundaryConditions(i))
            CASE(0)
                if (i < NumberXNodes) then
                    self % c_tri(i) = 1.0d0/world%dx_dl(i)
                end if
                if (i > 1) then
                    self%a_tri(i-1) = 1.0d0/ world%dx_dl(i-1)
                end if
                self%b_tri(i) = - (1.0d0/world%dx_dl(i-1)  + 1.0d0/world%dx_dl(i))
            CASE(1)
                self%b_tri(i) = 1.0d0
            CASE(2)
                if (i == 1) then
                    self % c_tri(i) = 1.0d0/world%dx_dl(i)
                    !self%a_tri(i - leftNodeIdx) = 2.0d0/(world%dx_dl(i-1) + world%dx_dl(i)) / world%dx_dl(i-1)
                    self%b_tri(i) = - (1.0d0/world%dx_dl(i))
                else if (i == NumberXNodes) then
                    self % a_tri(i-1) = 1.0d0/ world%dx_dl(i-1)
                    !self % c_tri(i - leftNodeIdx) = 2.0d0/(world%dx_dl(i-2) + world%dx_dl(i-1))/world%dx_dl(i-1)
                    self%b_tri(i) = - (1.0d0/world%dx_dl(i-1))
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

    subroutine solve_tridiag_Poisson(self, world)
        ! Tridiagonal (Thomas algorithm) solver for initial Poisson
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        integer(int32) :: i
        real(real64) :: m, d(NumberXNodes), cp(NumberXNodes-1),dp(NumberXNodes)

        do i=1, NumberXNodes
            SELECT CASE (world%boundaryConditions(i))
            CASE(0,2)
                d(i) = (-self%rho(i) - self%rho_const) / eps_0
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

    function getError_tridiag_Poisson(self, world) result(res)
        class(potentialSolver), intent(in) :: self
        type(Domain), intent(in) :: world
        integer(int32) :: i
        real(real64) :: Ax(NumberXNodes), res
        Ax = triMul(NumberXNodes, self%a_tri, self%c_tri, self%b_tri, self%phi)
        res = 0.0d0
        do i = 1, NumberXNodes
            SELECT CASE (world%boundaryConditions(i))
            CASE(0,2)
                res = res + (Ax(i)*eps_0/(-self%rho(i) - self%rho_const + 1d-15) - 1.0d0)**2
                !d(i) = -self%rho(i) - self%rho_const
            CASE(1,3)
                res = res + ((Ax(i) + 1d-15)/(self.phi_f(i) + 1d-15) - 1.0d0)**2
                !d(i) = self%phi_f(i)*eps_0
            END SELECT
        end do
        !res = Ax*eps_0 - 
        res = SQRT(res/NumberXNodes)

    end function getError_tridiag_Poisson

    subroutine solve_tridiag_Ampere(self, world, del_t)
        ! Tridiagonal (Thomas algorithm) solver for Ampere
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t
        integer(int32) :: i !n size dependent on how many points are boundary conditions and thus moved to rhs equation
        real(real64) :: m, d(NumberXNodes), cp(NumberXNodes-1),dp(NumberXNodes)

        do i =1, NumberXNodes
            SELECT CASE (world%boundaryConditions(i))
            CASE(0)
                d(i) = (SUM(self%J(i, :)) - SUM(self%J(i-1, :))) * del_t / eps_0 - (self%phi(i) - self%phi(i-1))/world%dx_dl(i-1) + (self%phi(i+1) - self%phi(i))/world%dx_dl(i)
            CASE(1,3)
                d(i) = self%phi_f(i)
            CASE(2)
                if (i == 1) then
                    d(i) = (del_t * SUM(self%J(i, :))/eps_0 - (self%phi(i) - self%phi(i+1))/world%dx_dl(i))
                else if (i == NumberXNodes) then
                    d(i) = (-del_t * SUM(self%J(i-1, :))/eps_0 - (self%phi(i) - self%phi(i-1))/world%dx_dl(i-1))
                end if
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
    ! initialize phi_f
        !d = self%phi_f
        self%phi_f(NumberXNodes) = dp(NumberXNodes)
    ! solve for x from the vectors c-prime and d-prime
        do i = NumberXNodes-1, 1, -1
            self%phi_f(i) = dp(i)-cp(i)*self%phi_f(i+1)
        end do
    end subroutine solve_tridiag_Ampere

    subroutine solve_CG_Ampere(self, world, del_t)
        ! Tridiagonal (Thomas algorithm) solver for initial Poisson
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t
        integer(int32) :: i
        real(real64) :: b(NumberXNodes), RPast(NumberXNodes), RFuture(NumberXNodes), D(NumberXNodes), beta, alpha, resPast, resFuture, Ax(NumberXNodes)
        logical :: converge
        converge = .false.
        do i =1, NumberXNodes
            SELECT CASE (world%boundaryConditions(i))
            CASE(0)
                b(i) = (SUM(self%J(i, :)) - SUM(self%J(i-1, :))) * del_t / eps_0 - (self%phi(i) - self%phi(i-1))/world%dx_dl(i-1) + (self%phi(i+1) - self%phi(i))/world%dx_dl(i)
            CASE(1,3)
                b(i) = self%phi_f(i)
            CASE(2)
                if (i == 1) then
                    b(i) = (del_t * SUM(self%J(i, :))/eps_0 - (self%phi(i) - self%phi(i+1))/world%dx_dl(i))
                else if (i == NumberXNodes) then
                    b(i) = (-del_t * SUM(self%J(i-1, :))/eps_0 - (self%phi(i) - self%phi(i-1))/world%dx_dl(i-1))
                end if
            END SELECT
        end do
        RPast = b - triMul(NumberXNodes, self%a_tri, self%c_tri, self%b_tri, self%phi_f)
        resPast = SQRT(SUM(RPast**2))
        D = RPast
        do i = 1, 1000
            alpha = resPast**2 / SUM(D * triMul(NumberXNodes, self%a_tri, self%c_tri, self%b_tri, D))
            self%phi_f = self%phi_f + alpha * D
            Ax = triMul(NumberXNodes, self%a_tri, self%c_tri, self%b_tri, self%phi_f)
            RFuture = b - Ax
            resFuture = SQRT(SUM(RFuture**2))
            if (SUM(ABS(RFuture/(b + 1.d-15)))/NumberXNodes < eps_r * 1.d-2) then
                converge = .true.
                exit
            end if
            beta = (resFuture**2)/(resPast**2)
            D = RFuture + beta * D
            RPast = RFuture
            resPast = resFuture
        end do
        if (.not. converge) then
            stop 'Max iterations reached CG solver!'
        end if
    end subroutine solve_CG_Ampere

    function getError_tridiag_Ampere(self, world, del_t) result(res)
        class(potentialSolver), intent(in) :: self
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t
        integer(int32) :: i
        real(real64) :: Ax(NumberXNodes), d(NumberXNodes), res(NumberXNodes)
        Ax = triMul(NumberXNodes, self%a_tri, self%c_tri, self%b_tri, self%phi_f)
        do i =1, NumberXNodes
            SELECT CASE (world%boundaryConditions(i))
            CASE(0)
                d(i) = (SUM(self%J(i, :)) - SUM(self%J(i-1, :))) * del_t / eps_0 - (self%phi(i) - self%phi(i-1))/world%dx_dl(i-1) + (self%phi(i+1) - self%phi(i))/world%dx_dl(i)
            CASE(1,3)
                d(i) = self%phi_f(i)
            CASE(2)
                if (i == 1) then
                    d(i) = -(del_t / eps_0) * (-SUM(self%J(i,:)) + (self%phi(i) - self%phi(i+1))/world%dx_dl(i))
                else if (i == NumberXNodes) then
                    d(i) = -(del_t / eps_0) * (-SUM(self%J(i-1, :)) + (self%phi(i-1) - self%phi(i))/world%dx_dl(i-1))
                end if
            END SELECT
        end do
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


end module mod_potentialSolver