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
        real(real64), allocatable :: phi(:), J(:, :), rho(:), phi_f(:), EField(:) !phi_f is final phi, will likely need to store two arrays for phi, can't be avoided
        real(real64) :: rho_const, BFieldMag, BField(3), BFieldAngle
        real(real64), allocatable :: a_tri(:), b_tri(:), c_tri(:), dirichletVals(:), sourceTermVals(:) !for thomas algorithm potential solver, a_tri is lower diagonal, b_tri middle, c_tri upper
        integer(int32), allocatable :: dirichletIndx(:)
        integer(int32) :: numDirichletNodes
        logical :: BFieldBool


    contains
        procedure, public, pass(self) :: solve_tridiag_Poisson
        procedure, public, pass(self) :: solve_tridiag_Ampere
        procedure, public, pass(self) :: getTotalPE
        procedure, public, pass(self) :: evaluateEFieldHalfTime
        procedure, public, pass(self) :: evaluateEFieldCurrTime
        procedure, public, pass(self) :: getError_tridiag_Ampere
        procedure, public, pass(self) :: solve_CG_Ampere
        procedure, public, pass(self) :: getError_tridiag_Poisson
        procedure, public, pass(self) :: construct_diagMatrix
    end type

    interface potentialSolver
        module procedure :: potentialSolver_constructor
    end interface potentialSolver

contains

    type(potentialSolver) function potentialSolver_constructor(world, leftVoltage, rightVoltage, BFieldMag, angle) result(self)
        ! Construct domain object, initialize grid, dx_dl, and nodeVol.
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: leftVoltage, rightVoltage
        real(real64), intent(in) :: BFieldMag, angle
        real(real64) :: angle_rad
        integer(int32) :: i
        allocate(self % J(NumberXNodes+1, numThread), self%rho(NumberXNodes), self % phi(NumberXNodes), self % phi_f(NumberXNodes), self%a_tri(NumberXNodes-1), &
        self%b_tri(NumberXNodes), self%c_tri(NumberXNodes-1), self%EField(NumberXNodes+1))
        self % rho = 0.0d0
        self % J = 0.0d0
        self % phi = 0.0d0
        self%EField = 0.0d0
        self % rho_const = 0.0d0
        self%phi_f = self%phi 
        self%BFieldMag = BFieldMag
        self%BFieldBool = (BFieldMag /= 0.0d0)
        angle_rad = angle * pi / 180.0d0
        self%BFieldAngle = angle_rad
        self%BField(1) = BFieldMag * COS(angle_rad)
        self%BField(2) = BFieldMag * SIN(angle_rad)
        self%BField(3) = 0.0d0 
        self%EField = 0.0d0
        self%numDirichletNodes = 0
        if (ABS(world%boundaryConditions(1)) == 1) then
            self%numDirichletNodes = self%numDirichletNodes + 1
        end if
        if (ABS(world%boundaryConditions(NumberXNodes+1)) == 1) then
            self%numDirichletNodes = self%numDirichletNodes + 1
        end if

        allocate(self%dirichletIndx(self%numDirichletNodes), self%dirichletVals(self%numDirichletNodes), self%sourceTermVals(self%numDirichletNodes))
        i = 0
        if (world%boundaryConditions(1) == 1) then
            i = i + 1
            self%dirichletIndx(i) = 1
            self%dirichletVals(i) = leftVoltage
        end if
        if (world%boundaryConditions(NumberXNodes+1) == 1) then
            i = i + 1
            self%dirichletIndx(i) = NumberXNodes
            self%dirichletVals(i) = rightVoltage
        end if

    end function potentialSolver_constructor

    subroutine construct_diagMatrix(self, world)
        ! construct diagonal components for thomas algorithm
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        integer(int32) :: i
        ! Construct primary matrix coefficients
        do i = 1, NumberXNodes
            if (world%boundaryConditions(i) == 0 .and. world%boundaryConditions(i+1) == 0) then
                self%b_tri(i) = -2.0d0 * (1.0d0 / (world%dx_dl(i-1) + world%dx_dl(i)) + 1.0d0/(world%dx_dl(i+1) + world%dx_dl(i)))
                self%a_tri(i-1) = 2.0d0 / (world%dx_dl(i-1) + world%dx_dl(i))
                self%c_tri(i) = 2.0d0 / (world%dx_dl(i) + world%dx_dl(i+1))
            else if (world%boundaryConditions(i+1) == 1) then
                self%b_tri(i) = -2.0d0 * (1.0d0 / (world%dx_dl(i-1) + world%dx_dl(i)) + 1.0d0/world%dx_dl(i))
                self%a_tri(i-1) = 2.0d0 / (world%dx_dl(i-1) + world%dx_dl(i))
            else if (world%boundaryConditions(i) == 1) then
                self%b_tri(i) = -2.0d0 * (1.0d0 / (world%dx_dl(i+1) + world%dx_dl(i)) + 1.0d0/world%dx_dl(i))
                self%c_tri(i) = 2.0d0 / (world%dx_dl(i+1) + world%dx_dl(i))
            else if (world%boundaryConditions(i+1) == 2) then
                self%b_tri(i) = -2.0d0 / (world%dx_dl(i-1) + world%dx_dl(i))
                self%a_tri(i-1) = 2.0d0 / (world%dx_dl(i-1) + world%dx_dl(i))
            else if (world%boundaryConditions(i) == 2) then
                self%b_tri(i) = -2.0d0 / (world%dx_dl(i+1) + world%dx_dl(i))
                self%c_tri(i) = 2.0d0 / (world%dx_dl(i+1) + world%dx_dl(i))
            end if
        end do

        ! Construct Source term additions due to boundary conditions
        if (self%numDirichletNodes > 0) then
            do i = 1, self%numDirichletNodes
                self%sourceTermVals(i) = - self%dirichletVals(i) * 2.0d0 / world%dx_dl(self%dirichletIndx(i))
            end do
        end if


    end subroutine construct_diagMatrix

    subroutine solve_tridiag_Poisson(self, world)
        ! Tridiagonal (Thomas algorithm) solver for initial Poisson
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        integer(int32) :: i
        real(real64) :: m, d(NumberXNodes), cp(NumberXNodes-1),dp(NumberXNodes)

        d = (-self%rho - self%rho_const) / eps_0
        do i=1, self%numDirichletNodes
            d(self%dirichletIndx(i)) = d(self%dirichletIndx(i)) + self%sourceTermVals(i)
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
        !self%EField = (self%phi(1:NumberXNodes-1) - self%phi(2:NumberXNodes)) / world%dx_dl
    end subroutine solve_tridiag_Poisson

    function getError_tridiag_Poisson(self, world) result(res)
        class(potentialSolver), intent(in) :: self
        type(Domain), intent(in) :: world
        integer(int32) :: i
        real(real64) :: Ax(NumberXNodes), d(NumberXNodes), res
        Ax = triMul(NumberXNodes, self%a_tri, self%c_tri, self%b_tri, self%phi)
        d = 1.0d0 - Ax/((-self%rho - self%rho_const) / eps_0 + 1.d-15)
        
        !res = Ax*eps_0 - 
        res = SQRT(SUM(d**2)/NumberXNodes)

    end function getError_tridiag_Poisson

    subroutine solve_tridiag_Ampere(self, world, del_t)
        ! Tridiagonal (Thomas algorithm) solver for Ampere
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t
        integer(int32) :: i !n size dependent on how many points are boundary conditions and thus moved to rhs equation
        real(real64) :: m, d(NumberXNodes), cp(NumberXNodes-1),dp(NumberXNodes)

        call self%evaluateEFieldCurrTime(world)
        do i =1, NumberXNodes
            SELECT CASE (world%boundaryConditions(i+1) - world%boundaryConditions(i))
            CASE(0)
                d(i) = (SUM(self%J(i+1, :)) - SUM(self%J(i, :))) * del_t / eps_0 - self%EField(i+1) + self%EField(i)
            CASE(-1)
                d(i) = (SUM(self%J(i+1, :)) - 2.0d0 * SUM(self%J(i, :))) * del_t / eps_0 - self%EField(i+1) + self%EField(i)
            CASE(1)
                d(i) = (2.0d0 * SUM(self%J(i+1, :)) - SUM(self%J(i, :))) * del_t / eps_0 - self%EField(i+1) + self%EField(i)
            CASE(-2)
                d(i) = SUM(self%J(i+1, :)) * del_t / eps_0 - self%EField(i+1)
            CASE(2)
                d(i) = -SUM(self%J(i, :)) * del_t / eps_0 + self%EField(i)
            END SELECT
        end do

        do i=1, self%numDirichletNodes
            d(self%dirichletIndx(i)) = d(self%dirichletIndx(i)) + self%sourceTermVals(i)
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
        !self%EField = 0.5d0 * (self%phi(1:NumberXNodes-1) + self%phi_f(1:NumberXNodes-1) - self%phi(2:NumberXNodes) - self%phi_f(2:NumberXNodes)) / world%dx_dl
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

    subroutine evaluateEFieldHalfTime(self, world)
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        integer(int32) :: i
        ! 'logical' field, dphi_dl
        self%EField(2:NumberXNodes) = 0.5d0 * (self%phi(1:NumberXNodes-1) + self%phi_f(1:NumberXNodes-1) - self%phi(2:NumberXNodes) - self%phi_f(2:NumberXNodes))

        do i = 1, self%numDirichletNodes
            if (self%dirichletIndx(i) == 1) then
                self%EField(1) = (2.0d0 * self%dirichletVals(i) - self%phi(1) - self%phi_f(1))
            else
                self%EField(NumberXNodes+1) = (self%phi(NumberXNodes) + self%phi_f(NumberXNodes) - 2.0d0 * self%dirichletVals(i))
            end if
        end do
    end subroutine evaluateEFieldHalfTime

    subroutine evaluateEFieldCurrTime(self, world)
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        integer(int32) :: i
        self%EField(2:NumberXNodes) = (self%phi(1:NumberXNodes-1) - self%phi(2:NumberXNodes)) &
            /world%centerDiff

        do i = 1, self%numDirichletNodes
            if (self%dirichletIndx(i) == 1) then
                self%EField(1) = 2.0d0 * (self%dirichletVals(i) - self%phi(1))/world%dx_dl(1)
            else
                self%EField(NumberXNodes+1) = 2.0d0 * (self%phi(NumberXNodes) - self%dirichletVals(i))/world%dx_dl(NumberXNodes)
            end if
        end do
    end subroutine evaluateEFieldCurrTime

    function getTotalPE(self, world, future) result(res)
        ! Get energy in electric fields, future true, then derive from phi_f, otherwise phi
        ! In 1D is J/m^2
        class(potentialSolver), intent(in) :: self
        type(Domain), intent(in) :: world
        logical :: future
        integer(int32) :: i
        real(real64) :: res
        if (future) then
            res = 0.5 * eps_0 * SUM(arrayDiff(self%phi_f, NumberXNodes)**2 / world%centerDiff)
            do i = 1, self%numDirichletNodes
                res = res + eps_0 * (self%dirichletVals(i) - self%phi_f(i))**2 / world%dx_dl(i)
            end do
        else
            res = 0.5 * eps_0 * SUM(arrayDiff(self%phi, NumberXNodes)**2 / world%centerDiff)
            do i = 1, self%numDirichletNodes
                res = res + eps_0 * (self%dirichletVals(i) - self%phi(i))**2 / world%dx_dl(i)
            end do
        end if
    end function getTotalPE


end module mod_potentialSolver