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
    public :: potentialSolver, readSolver

    type :: potentialSolver
        real(real64), allocatable :: phi(:), J(:), rho(:), phi_f(:), EField(:) !phi_f is final phi, will likely need to store two arrays for phi, can't be avoided
        real(real64) :: rho_const, BFieldMag, BField(3), BFieldAngle, RF_rad_frequency, RF_half_amplitude, RF_voltage
        real(real64), allocatable :: a_tri(:), b_tri(:), c_tri(:), dirichletVals(:), sourceTermVals(:) !for thomas algorithm potential solver, a_tri is lower diagonal, b_tri middle, c_tri upper
        integer(int32), allocatable :: dirichletIndx(:)
        logical, allocatable :: dirichletIsRFBool(:)
        integer(int32) :: numDirichletNodes
        logical :: BFieldBool, RF_bool


    contains
        procedure, public, pass(self) :: solve_tridiag_Poisson
        procedure, public, pass(self) :: solve_tridiag_Ampere
        procedure, public, pass(self) :: getChargeContinuityError
        procedure, public, pass(self) :: getTotalPE
        procedure, public, pass(self) :: setRFVoltage
        procedure, public, pass(self) :: resetVoltage
        procedure, public, pass(self) :: aveRFVoltage
        procedure, public, pass(self) :: getEnergyFromBoundary
        procedure, public, pass(self) :: evaluateEFieldHalfTime
        procedure, public, pass(self) :: evaluateEFieldCurrTime
        procedure, public, pass(self) :: evaluateEFieldFutureTime
        procedure, public, pass(self) :: getError_tridiag_Ampere
        ! procedure, public, pass(self) :: solve_CG_Ampere
        procedure, public, pass(self) :: getError_tridiag_Poisson
        procedure, public, pass(self) :: construct_diagMatrix
    end type

    interface potentialSolver
        module procedure :: potentialSolver_constructor
    end interface potentialSolver

contains

    type(potentialSolver) function potentialSolver_constructor(world, leftVoltage, rightVoltage, BFieldMag, angle, RF_frequency) result(self)
        ! Construct domain object, initialize grid, dx_dl, and nodeVol.
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: leftVoltage, rightVoltage
        real(real64), intent(in) :: BFieldMag, angle, RF_frequency
        real(real64) :: angle_rad
        integer(int32) :: i
        allocate(self % J(NumberXHalfNodes), self%rho(NumberXNodes), self % phi(NumberXNodes), self % phi_f(NumberXNodes), self%a_tri(NumberXNodes-1), &
        self%b_tri(NumberXNodes), self%c_tri(NumberXNodes-1), self%EField(NumberXHalfNodes))
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
        self%RF_half_amplitude = 0.0d0
        if (world%boundaryConditions(1)== 1 .or. world%boundaryConditions(1)== 4) then
            self%numDirichletNodes = self%numDirichletNodes + 1
        end if
        if (world%boundaryConditions(NumberXHalfNodes)== 1 .or. world%boundaryConditions(NumberXHalfNodes)== 4) then
            self%numDirichletNodes = self%numDirichletNodes + 1
        end if
        

        allocate(self%dirichletIndx(self%numDirichletNodes), self%dirichletVals(self%numDirichletNodes), self%sourceTermVals(self%numDirichletNodes), self%dirichletIsRFBool(self%numDirichletNodes))
        i = 0
        if (world%boundaryConditions(1) == 1) then
            i = i + 1
            self%dirichletIndx(i) = 1
            self%dirichletVals(i) = leftVoltage
            self%dirichletIsRFBool(i) = .false.
        else if (world%boundaryConditions(1) == 4) then
            i = i + 1
            self%dirichletIndx(i) = 1
            self%RF_half_amplitude = leftVoltage
            self%dirichletVals(i) = 0.0d0
            self%dirichletIsRFBool(i) = .true.
            self%RF_voltage = 0.0d0
        end if
        if (world%boundaryConditions(NumberXHalfNodes) == 1) then
            i = i + 1
            self%dirichletIndx(i) = NumberXNodes
            self%dirichletVals(i) = rightVoltage
            self%dirichletIsRFBool(i) = .false.
        else if (world%boundaryConditions(NumberXHalfNodes) == 4) then
            if (self%RF_half_amplitude /= 0) then
                print *, 'Half amplitude voltage for RF already set, have two RF boundaries!'
                stop
            else
                i = i + 1
                self%dirichletIndx(i) = NumberXNodes
                self%RF_half_amplitude = rightVoltage
                self%dirichletVals(i) = 0.0d0
                self%dirichletIsRFBool(i) = .true.
                self%RF_voltage = 0.0d0
            end if
        end if
        self%RF_rad_frequency = 2.0d0 * pi * RF_frequency
        self%RF_bool = self%RF_rad_frequency > 0 .and. self%RF_half_amplitude > 0 .and. (world%boundaryConditions(1) == 4 .or. world%boundaryConditions(NumberXHalfNodes) == 4)
        call self%construct_diagMatrix(world)
    end function potentialSolver_constructor

    subroutine construct_diagMatrix(self, world)
        ! construct diagonal components for thomas algorithm
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        integer(int32) :: i
        ! Construct primary matrix coefficients
        self%a_tri = 0.0
        self%b_tri = 0.0
        self%c_tri = 0.0
        do i = 1, NumberXNodes
            if (world%boundaryConditions(i) == 0 .and. world%boundaryConditions(i+1) == 0) then
                self%b_tri(i) = -2.0d0 * (1.0d0 / (world%dx_dl(i-1) + world%dx_dl(i)) + 1.0d0/(world%dx_dl(i+1) + world%dx_dl(i)))
                self%a_tri(i-1) = 2.0d0 / (world%dx_dl(i-1) + world%dx_dl(i))
                self%c_tri(i) = 2.0d0 / (world%dx_dl(i) + world%dx_dl(i+1))
            else if (world%boundaryConditions(i+1) == 1 .or. world%boundaryConditions(i+1) == 4 .or. world%boundaryConditions(i+1) == 3) then
                self%b_tri(i) = -2.0d0 * (1.0d0 / (world%dx_dl(i-1) + world%dx_dl(i)) + 1.0d0/world%dx_dl(i))
                self%a_tri(i-1) = 2.0d0 / (world%dx_dl(i-1) + world%dx_dl(i))
            else if (world%boundaryConditions(i) == 1 .or. world%boundaryConditions(i) == 4 .or. world%boundaryConditions(i) == 3) then
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

    subroutine solve_tridiag_Poisson(self, world, timeCurrent)
        ! Tridiagonal (Thomas algorithm) solver for initial Poisson
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: timeCurrent
        integer(int32) :: i
        real(real64) :: m, d(NumberXNodes), cp(NumberXNodes-1),dp(NumberXNodes)

        d = (-self%rho - self%rho_const) / eps_0
        do i=1, self%numDirichletNodes
            if (self%dirichletIsRFBool(i)) then
                self%dirichletVals(i) = self%RF_half_amplitude * SIN(self%RF_rad_frequency * timeCurrent)
                self%sourceTermVals(i) = -self%dirichletVals(i) * 2.0d0 / world%dx_dl(self%dirichletIndx(i))
            end if
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
        Ax = triMul(NumberXNodes, self%a_tri, self%c_tri, self%b_tri, self%phi_f)
        d = (-self%rho - self%rho_const) / eps_0
        do i=1, self%numDirichletNodes
            d(self%dirichletIndx(i)) = d(self%dirichletIndx(i)) + self%sourceTermVals(i)
        end do
        d = 1.0d0 - Ax/(d + 1.d-15)
        
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
                d(i) = (self%J(i+1) - self%J(i)) * del_t / eps_0 - self%EField(i+1) + self%EField(i)
            CASE(-1,-4)
                d(i) = (self%J(i+1) - 2.0d0 * self%J(i)) * del_t / eps_0 - self%EField(i+1) + self%EField(i)
            CASE(1,4)
                d(i) = (2.0d0 * self%J(i+1) - self%J(i)) * del_t / eps_0 - self%EField(i+1) + self%EField(i)
            CASE(-2)
                d(i) = self%J(i+1) * del_t / eps_0 - self%EField(i+1)
            CASE(2)
                d(i) = -self%J(i) * del_t / eps_0 + self%EField(i)
            CASE(-3)
                d(i) = (self%J(2) - self%J(1) - self%J(NumberXHalfNodes)) * del_t / eps_0 - self%EField(2) + self%EField(1)
            CASE(3)
                d(i) = (self%J(1) + self%J(NumberXHalfNodes) - self%J(NumberXNodes)) * del_t / eps_0 - self%EField(NumberXHalfNodes) + self%EField(NumberXNodes)
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

    ! subroutine solve_CG_Ampere(self, world, del_t)
    !     ! Tridiagonal (Thomas algorithm) solver for initial Poisson
    !     class(potentialSolver), intent(in out) :: self
    !     type(Domain), intent(in) :: world
    !     real(real64), intent(in) :: del_t
    !     integer(int32) :: i
    !     real(real64) :: b(NumberXNodes), RPast(NumberXNodes), RFuture(NumberXNodes), D(NumberXNodes), beta, alpha, resPast, resFuture, Ax(NumberXNodes)
    !     logical :: converge
    !     converge = .false.
    !     do i =1, NumberXNodes
    !         SELECT CASE (world%boundaryConditions(i))
    !         CASE(0)
    !             b(i) = (SUM(self%J(i, :)) - SUM(self%J(i-1, :))) * del_t / eps_0 - (self%phi(i) - self%phi(i-1))/world%dx_dl(i-1) + (self%phi(i+1) - self%phi(i))/world%dx_dl(i)
    !         CASE(1,3)
    !             b(i) = self%phi_f(i)
    !         CASE(2)
    !             if (i == 1) then
    !                 b(i) = (del_t * SUM(self%J(i, :))/eps_0 - (self%phi(i) - self%phi(i+1))/world%dx_dl(i))
    !             else if (i == NumberXNodes) then
    !                 b(i) = (-del_t * SUM(self%J(i-1, :))/eps_0 - (self%phi(i) - self%phi(i-1))/world%dx_dl(i-1))
    !             end if
    !         END SELECT
    !     end do
    !     RPast = b - triMul(NumberXNodes, self%a_tri, self%c_tri, self%b_tri, self%phi_f)
    !     resPast = SQRT(SUM(RPast**2))
    !     D = RPast
    !     do i = 1, 1000
    !         alpha = resPast**2 / SUM(D * triMul(NumberXNodes, self%a_tri, self%c_tri, self%b_tri, D))
    !         self%phi_f = self%phi_f + alpha * D
    !         Ax = triMul(NumberXNodes, self%a_tri, self%c_tri, self%b_tri, self%phi_f)
    !         RFuture = b - Ax
    !         resFuture = SQRT(SUM(RFuture**2))
    !         if (SUM(ABS(RFuture/(b + 1.d-15)))/NumberXNodes < eps_r * 1.d-2) then
    !             converge = .true.
    !             exit
    !         end if
    !         beta = (resFuture**2)/(resPast**2)
    !         D = RFuture + beta * D
    !         RPast = RFuture
    !         resPast = resFuture
    !     end do
    !     if (.not. converge) then
    !         stop 'Max iterations reached CG solver!'
    !     end if
    ! end subroutine solve_CG_Ampere

    function getError_tridiag_Ampere(self, world, del_t) result(res)
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t
        integer(int32) :: i
        real(real64) :: Ax(NumberXNodes), d(NumberXNodes), res(NumberXNodes)
        Ax = triMul(NumberXNodes, self%a_tri, self%c_tri, self%b_tri, self%phi_f)
        call self%evaluateEFieldCurrTime(world)
        do i =1, NumberXNodes
            SELECT CASE (world%boundaryConditions(i+1) - world%boundaryConditions(i))
            CASE(0)
                d(i) = (self%J(i+1) - self%J(i)) * del_t / eps_0 - self%EField(i+1) + self%EField(i)
            CASE(-1,-4)
                d(i) = (self%J(i+1) - 2.0d0 * self%J(i)) * del_t / eps_0 - self%EField(i+1) + self%EField(i)
            CASE(1,4)
                d(i) = (2.0d0 * self%J(i+1) - self%J(i)) * del_t / eps_0 - self%EField(i+1) + self%EField(i)
            CASE(-2)
                d(i) = self%J(i+1) * del_t / eps_0 - self%EField(i+1)
            CASE(2)
                d(i) = -self%J(i) * del_t / eps_0 + self%EField(i)
            CASE(-3)
                d(i) = (self%J(2) - self%J(1) - self%J(NumberXHalfNodes)) * del_t / eps_0 - self%EField(2) + self%EField(1)
            CASE(3)
                d(i) = (self%J(1) + self%J(NumberXHalfNodes) - self%J(NumberXNodes)) * del_t / eps_0 - self%EField(NumberXHalfNodes) + self%EField(NumberXNodes)
            END SELECT
        end do

        do i=1, self%numDirichletNodes
            d(self%dirichletIndx(i)) = d(self%dirichletIndx(i)) + self%sourceTermVals(i)
        end do
        res = Ax- d

    end function getError_tridiag_Ampere

    function getChargeContinuityError(self, rho_i, world, del_t) result(chargeError)
        class(potentialSolver), intent(in) :: self
        real(real64), intent(in) :: del_t, rho_i(NumberXNodes)
        type(Domain), intent(in) :: world
        integer(int32) :: i, k
        real(real64) :: chargeError, del_Rho
        chargeError = 0.0d0
        k = 0
        do i = 1, NumberXNodes
            del_Rho = self%rho(i) - rho_i(i)
            if (del_Rho /= 0) then
                SELECT CASE (world%boundaryConditions(i+1) - world%boundaryConditions(i))
                CASE(0)
                    chargeError = chargeError + (1.0d0 + del_t * (self%J(i+1) - self%J(i))/(del_Rho))**2
                CASE(-1,-4)
                    chargeError = chargeError + (1.0d0 + del_t * (self%J(i+1) - 2.0d0 * self%J(i))/(del_Rho))**2
                CASE(1,4)
                    chargeError = chargeError + (1.0d0 + del_t * (2.0d0 * self%J(i+1) - self%J(i))/del_Rho)**2
                CASE(2)
                    chargeError = chargeError + (1.0d0 + del_t * (-self%J(i))/del_Rho)**2 
                CASE(-2)
                    chargeError = chargeError + (1.0d0 + del_t * (self%J(i+1))/del_Rho)**2
                CASE(3)
                    chargeError = chargeError + (1.0d0 + del_t * (self%J(NumberXHalfNodes) + self%J(1) - self%J(NumberXNodes))/(del_Rho))**2
                CASE(-3)
                    chargeError = chargeError + (1.0d0 + del_t * (self%J(2) - self%J(1) - self%J(NumberXHalfNodes))/(del_Rho))**2
                END SELECT
                k = k + 1
            end if
        end do
        chargeError = SQRT(chargeError/k)
    end function getChargeContinuityError

    subroutine evaluateEFieldHalfTime(self, boundaryConditionsLeftEdge)
        class(potentialSolver), intent(in out) :: self
        integer(int32), intent(in) :: boundaryConditionsLeftEdge
        integer(int32) :: i
        ! 'logical' field, dphi_dl, used in particle mover
        self%EField(2:NumberXNodes) = 0.5d0 * (self%phi(1:NumberXNodes-1) + self%phi_f(1:NumberXNodes-1) - self%phi(2:NumberXNodes) - self%phi_f(2:NumberXNodes))

        do i = 1, self%numDirichletNodes
            if (self%dirichletIndx(i) == 1) then
                if (self%dirichletIsRFBool(i)) then
                    self%EField(1) = (self%dirichletVals(i) + self%RF_voltage - self%phi(1) - self%phi_f(1))
                else
                    self%EField(1) = (2.0d0 * self%dirichletVals(i) - self%phi(1) - self%phi_f(1))
                end if
            else
                if (self%dirichletIsRFBool(i)) then
                    self%EField(NumberXHalfNodes) = self%phi(NumberXNodes) + self%phi_f(NumberXNodes) - self%dirichletVals(i) - self%RF_voltage 
                else
                    self%EField(NumberXHalfNodes) = self%phi(NumberXNodes) + self%phi_f(NumberXNodes) - 2.0d0 * self%dirichletVals(i)
                end if
            end if
        end do

        if (boundaryConditionsLeftEdge == 3) then
            self%EField(1) = 0.5d0 * (self%phi(NumberXNodes) + self%phi_f(NumberXNodes) - self%phi(1) - self%phi_f(1))
            self%EField(NumberXHalfNodes) = self%EField(1)
        end if

    end subroutine evaluateEFieldHalfTime

    subroutine evaluateEFieldCurrTime(self, world)
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        integer(int32) :: i
        ! used as source term in divergence ampere solver
        self%EField(2:NumberXNodes) = (self%phi(1:NumberXNodes-1) - self%phi(2:NumberXNodes)) &
            /world%centerDiff

        do i = 1, self%numDirichletNodes
            if (self%dirichletIndx(i) == 1) then
                self%EField(1) = 2.0d0 * (self%dirichletVals(i) - self%phi(1))/world%dx_dl(1)
            else
                self%EField(NumberXHalfNodes) = 2.0d0 * (self%phi(NumberXNodes) - self%dirichletVals(i))/world%dx_dl(NumberXNodes)
            end if
        end do

        if (world%boundaryConditions(1) == 3) then
            self%EField(1) = (self%phi(NumberXNodes) - self%phi(1)) / (0.5d0 * (world%dx_dl(1) + world%dx_dl(NumberXNodes)))
            self%EField(NumberXHalfNodes) = self%EField(1)
        end if
    end subroutine evaluateEFieldCurrTime

    subroutine evaluateEFieldFutureTime(self, world)
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        integer(int32) :: i
        ! used as source term in divergence ampere solver
        self%EField(2:NumberXNodes) = (self%phi_f(1:NumberXNodes-1) - self%phi_f(2:NumberXNodes)) &
            /world%centerDiff

        do i = 1, self%numDirichletNodes
            if (self%dirichletIndx(i) == 1) then
                if (self%dirichletIsRFBool(i)) then
                    self%EField(1) = 2.0d0 * (self%RF_voltage - self%phi_f(1))/world%dx_dl(1)
                else
                    self%EField(1) = 2.0d0 * (self%dirichletVals(i) - self%phi_f(1))/world%dx_dl(1)
                end if
            else
                if (self%dirichletIsRFBool(i)) then
                    self%EField(NumberXHalfNodes) = 2.0d0 * (self%phi_f(NumberXNodes) - self%RF_voltage)/world%dx_dl(NumberXNodes) 
                else
                    self%EField(NumberXHalfNodes) = 2.0d0 * (self%phi_f(NumberXNodes) - self%dirichletVals(i))/world%dx_dl(NumberXNodes)
                end if
            end if
        end do
        if (world%boundaryConditions(1) == 3) then
            self%EField(1) = (self%phi_f(NumberXNodes) - self%phi_f(1)) / (0.5d0 * (world%dx_dl(1) + world%dx_dl(NumberXNodes)))
            self%EField(NumberXHalfNodes) = self%EField(1)
        end if
    end subroutine evaluateEFieldFutureTime

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
                if (self%dirichletIsRFBool(i)) then
                    res = res + eps_0 * (self%RF_voltage - self%phi_f(self%dirichletIndx(i)))**2 / world%dx_dl(self%dirichletIndx(i))
                else
                    res = res + eps_0 * (self%dirichletVals(i) - self%phi_f(self%dirichletIndx(i)))**2 / world%dx_dl(self%dirichletIndx(i))
                end if
            end do
            if (world%boundaryConditions(1) == 3) then
                res = res + eps_0 * (self%phi_f(NumberXNodes) - self%phi_f(1))**2 / ((world%dx_dl(1) + world%dx_dl(NumberXNodes)))
            end if
        else
            res = 0.5 * eps_0 * SUM(arrayDiff(self%phi, NumberXNodes)**2 / world%centerDiff)
            do i = 1, self%numDirichletNodes
                res = res + eps_0 * (self%dirichletVals(i) - self%phi(self%dirichletIndx(i)))**2 / world%dx_dl(self%dirichletIndx(i))
            end do
            if (world%boundaryConditions(1) == 3) then
                res = res + eps_0 * (self%phi(NumberXNodes) - self%phi(1))**2 / ((world%dx_dl(1) + world%dx_dl(NumberXNodes)))
            end if
        end if
    end function getTotalPE

    function getEnergyFromBoundary(self, world, del_t) result(res)
        ! Get energy input into domain from dirichlet boundary
        ! In 1D is J/m^2
        class(potentialSolver), intent(in) :: self
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t
        real(real64) :: res, leftPhiAve, rightPhiAve
        res = 0.0d0
        if (self%numDirichletNodes > 1) then
            res = del_t * self%J(2)
            res = res + eps_0 * ((self%phi_f(1) - self%phi_f(2)) -(self%phi(1) - self%phi(2))) /world%centerDiff(1)
            if (self%dirichletIsRFBool(1)) then
                leftPhiAve = 0.5d0 * (self%RF_voltage + self%dirichletVals(1))
            else
                leftPhiAve = self%dirichletVals(1)
            end if

            if (self%dirichletIsRFBool(2)) then
                rightPhiAve = 0.5d0 * (self%RF_voltage + self%dirichletVals(2))
            else
                rightPhiAve = self%dirichletVals(2)
            end if
            res = res * (leftPhiAve - rightPhiAve)
        end if
    end function getEnergyFromBoundary

    subroutine setRFVoltage(self, world, timeFuture)
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: timeFuture 
        integer(int32) :: i

        do i = 1, self%numDirichletNodes
            if (self%dirichletIsRFBool(i)) then
                self%RF_voltage = self%RF_half_amplitude * SIN(self%RF_rad_frequency * (timeFuture))
                self%sourceTermVals(i) = -self%RF_voltage * 2.0d0 / world%dx_dl(self%dirichletIndx(i))
            end if
        end do

    end subroutine setRFVoltage

    subroutine resetVoltage(self)
        class(potentialSolver), intent(in out) :: self
        integer(int32) :: i
        self%phi = self%phi_f
        do i = 1, self%numDirichletNodes
            if (self%dirichletIsRFBool(i)) then
                self%dirichletVals(i) = self%RF_voltage
            end if
        end do
    end subroutine resetVoltage

    subroutine aveRFVoltage(self, accumBool, phi_average, RF_ave, numSteps, world)
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        logical, intent(in) :: accumBool
        real(real64), intent(in out) :: phi_average(NumberXNodes)
        real(real64), intent(in out) :: RF_ave
        integer(int32), intent(in) :: numSteps
        integer(int32) :: j

        if (accumBool) then
            phi_average = phi_average + self%phi_f
            ! average RF voltage
            do j = 1, self%numDirichletNodes
                if (self%dirichletIsRFBool(j)) then
                    RF_ave = RF_ave + self%RF_voltage
                end if
            end do
        else
            phi_average = phi_average/numSteps
            self%phi_f = phi_average
            do j = 1, self%numDirichletNodes
                if (self%dirichletIsRFBool(j)) then
                    RF_ave = RF_ave/numSteps
                    self%sourceTermVals(j) = -RF_ave * 2.0d0 / world%dx_dl(self%dirichletIndx(j))
                end if
            end do
        end if

    end subroutine aveRFVoltage

    ! ------------------ read solver input ---------------------------

    subroutine readSolver(GeomFilename, solver, world)
        type(Domain), intent(in) :: world
        type(potentialSolver), intent(in out) :: solver
        character(len=*), intent(in) :: GeomFilename
        integer(int32) :: io, i
        real(real64) :: leftVoltage, rightVoltage, BFieldMag, angle, RF_frequency, realTemp(20)

        print *, ""
        print *, "Reading solver inputs:"
        print *, "------------------"
        open(10,file='../InputData/'//GeomFilename)
        read(10, *, IOSTAT = io)
        read(10, *, IOSTAT = io)
        read(10, *, IOSTAT = io)
        read(10, *, IOSTAT = io)
        read(10, *, IOSTAT = io) leftVoltage, rightVoltage
        read(10, *, IOSTAT = io) RF_frequency
        close(10)
        if (restartBool) then
            leftVoltage = 0.0d0
            rightVoltage = 0.0d0
            open(10,file=restartDirectory//"/"//"InitialConditions.dat", IOSTAT=io)
            read(10, *, IOSTAT = io)
            read(10, *, IOSTAT = io) realTemp(1:15)
            close(10) 
            RF_frequency = realTemp(14)/2.0d0/pi
            if (world%boundaryConditions(1) == 4) then
                leftVoltage = realTemp(15)
            else if (world%boundaryConditions(NumberXHalfNodes) == 4) then
                rightVoltage = realTemp(15)
            end if   
        end if
        BFieldMag = 0.0d0
        angle = 0.0d0
        solver = potentialSolver(world, leftVoltage, rightVoltage, BFieldMag, angle, RF_frequency)
        if (solver%numDirichletNodes > 0) then
            do i = 1, solver%numDirichletNodes
                print *, 'Dirichlet potential near boundary', solver%dirichletIndx(i), 'is', solver%dirichletVals(i)
            end do
        end if
        print *, 'RF frequency:', solver%RF_rad_frequency/2.0d0/pi
        print *, 'RF_half_amplitude:', solver%RF_half_amplitude
        print *, 'RF_rad_frequency:', solver%RF_rad_frequency
        print *, "BField vector:", solver%BField
        print *, "------------------"
        print *, ""

    end subroutine readSolver

end module mod_potentialSolver