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
        real(real64), allocatable :: phi(:), J(:), rho(:), EField(:), EField_f(:), EField_half(:) !phi_f is final phi, will likely need to store two arrays for phi, can't be avoided
        real(real64) :: rho_const, BFieldMag, BField(3), BFieldAngle, RF_rad_frequency, RF_half_amplitude, aveJ
        real(real64), allocatable :: a_tri(:), b_tri(:), c_tri(:) !for thomas algorithm potential solver, a_tri is lower diagonal, b_tri middle, c_tri upper
        logical :: BFieldBool, RF_bool, globJBool


    contains
        procedure, public, pass(self) :: solve_tridiag_Poisson
        procedure, public, pass(self) :: solve_tridiag_Ampere
        procedure, public, pass(self) :: getChargeContinuityError
        procedure, public, pass(self) :: getTotalPE
        procedure, public, pass(self) :: setRFVoltage
        procedure, public, pass(self) :: aveRFVoltage
        procedure, public, pass(self) :: getEnergyFromBoundary
        procedure, public, pass(self) :: getError_tridiag_Ampere
        ! procedure, public, pass(self) :: solve_CG_Ampere
        procedure, public, pass(self) :: resetVoltage
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
        allocate(self % J(NumberXHalfNodes), self%rho(NumberXNodes), self % phi(NumberXNodes), self % EField_f(NumberXHalfNodes), self%EField_half(NumberXHalfNodes), self%a_tri(NumberXHalfNodes), &
        self%b_tri(NumberXNodes), self%c_tri(NumberXHalfNodes), self%EField(NumberXHalfNodes))
        call construct_diagMatrix(self, world)
        self % rho = 0.0d0
        self % J = 0.0d0
        self % phi = 0.0d0
        self % rho_const = 0.0d0
        self%RF_half_amplitude = 0.0d0
        if (world%boundaryConditions(1) == 1) self%phi(1) = leftVoltage
        if (world%boundaryConditions(1) == 4) self%RF_half_amplitude = leftVoltage
        if (world%boundaryConditions(NumberXNodes) == 1) self%phi(NumberXNodes) = rightVoltage
        if (world%boundaryConditions(NumberXNodes) == 4) then
            if (self%RF_half_amplitude /= 0) then
                print *, 'Half amplitude voltage for RF already set, have two RF boundaries!'
                stop
            else
                self%RF_half_amplitude = rightVoltage
            end if
        end if
        if (world%boundaryConditions(1) == 3) then
            self%phi(1) = leftVoltage
            self%phi(NumberXNodes) = leftVoltage
        end if
        self%RF_rad_frequency = 2.0d0 * pi * RF_frequency
        self%RF_bool = self%RF_rad_frequency > 0 .and. self%RF_half_amplitude > 0 .and. (world%boundaryConditions(1) == 4 .or. world%boundaryConditions(NumberXNodes) == 4)
        self%BFieldMag = BFieldMag
        self%BFieldBool = (BFieldMag /= 0.0d0)
        angle_rad = angle * pi / 180.0d0
        self%BFieldAngle = angle_rad
        self%BField(1) = BFieldMag * COS(angle_rad)
        self%BField(2) = BFieldMag * SIN(angle_rad)
        self%BField(3) = 0.0d0 
        self%EField = 0.0d0
        self%EField_f = 0.0d0
        self%aveJ = 0.0
        self%EField_half = 0.0d0
        self%globJBool = .false.
    end function potentialSolver_constructor

    subroutine construct_diagMatrix(self, world)
        ! construct diagonal components for thomas algorithm
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        integer(int32) :: i
        self%a_tri = 0.0d0
        self%b_tri = 0.0d0
        self%c_tri = 0.0d0
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
            CASE(1,4)
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

    subroutine solve_tridiag_Poisson(self, world, timeCurrent)
        ! Tridiagonal (Thomas algorithm) solver for initial Poisson
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: timeCurrent
        integer(int32) :: i
        real(real64) :: m, d(NumberXNodes), cp(NumberXHalfNodes),dp(NumberXNodes)

        do i=1, NumberXNodes
            SELECT CASE (world%boundaryConditions(i))
            CASE(0,2)
                d(i) = (-self%rho(i) - self%rho_const) / eps_0
            CASE(1,3)
                d(i) = self%phi(i)
            CASE(4)
                d(i) = self%RF_half_amplitude * SIN(self%RF_rad_frequency * timeCurrent)
            END SELECT
        end do
    ! initialize c-prime and d-prime
        cp(1) = self%c_tri(1)/self%b_tri(1)
        dp(1) = d(1)/self%b_tri(1)
    ! solve for vectors c-prime and d-prime
        do i = 2,NumberXHalfNodes
            m = self%b_tri(i)-cp(i-1)*self%a_tri(i-1)
            cp(i) = self%c_tri(i)/m
            dp(i) = (d(i)-dp(i-1)*self%a_tri(i-1))/m
        end do
        m = self%b_tri(NumberXNodes)-cp(NumberXHalfNodes)*self%a_tri(NumberXHalfNodes)
        dp(NumberXNodes) = (d(NumberXNodes)-dp(NumberXHalfNodes)*self%a_tri(NumberXHalfNodes))/m
    ! initialize x
        self%phi = dp
    ! solve for x from the vectors c-prime and d-prime
        do i = NumberXHalfNodes, 1, -1
            self%phi(i) = dp(i)-cp(i)*self%phi(i+1)
        end do
        self%EField = (self%phi(1:NumberXHalfNodes) - self%phi(2:NumberXNodes)) / world%dx_dl
        self%EField_f = self%EField
    end subroutine solve_tridiag_Poisson

    function getError_tridiag_Poisson(self, world) result(res)
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        integer(int32) :: i
        real(real64) :: Ax(NumberXNodes), res
        res = (eps_0 * (self%EField_f(1) - self%EField_f(NumberXHalfNodes))/(self%rho(1) + self%rho(NumberXNodes)) - 1.0d0)**2
        do i = 2, NumberXHalfNodes
            res = res + (eps_0 * (self%EField_f(i) - self%EField_f(i-1))/(self%rho(i)) - 1.0d0)**2
        end do
        res = SQRT(res/NumberXHalfNodes)

    end function getError_tridiag_Poisson

    subroutine solve_tridiag_Ampere(self, world, del_t)
        ! Tridiagonal (Thomas algorithm) solver for Ampere
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t
        if (.not. self%globJBool) self%aveJ = (SUM(self%J * world%dx_dl)/world%L_domain)
        self%EField_f = (self%aveJ - self%J) * del_t / eps_0 + self%EField
        !self%EField = 0.5d0 * (self%phi(1:NumberXHalfNodes) + self%phi_f(1:NumberXHalfNodes) - self%phi(2:NumberXNodes) - self%phi_f(2:NumberXNodes)) / world%dx_dl
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
        class(potentialSolver), intent(in) :: self
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t
        integer(int32) :: i
        real(real64) :: res(NumberXHalfNodes)
        res = eps_0 * (self%EField_f - self%EField)/del_t + self%J - SUM(self%J * world%dx_dl)/world%L_domain

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
                SELECT CASE (world%boundaryConditions(i))
                CASE(0)
                    chargeError = chargeError + (1.0d0 + del_t * (self%J(i) - self%J(i-1))/del_Rho)**2
                    k = k + 1
                CASE(1)
                    continue
                CASE(2)
                    if (i == 1) then
                        chargeError = chargeError + (1.0d0 + del_t * self%J(1)/del_Rho)**2
                    else
                        chargeError = chargeError + (1.0d0 - del_t * self%J(NumberXHalfNodes)/del_Rho)**2
                    end if
                    k = k + 1
                CASE(3)
                    if (i == 1) then
                        del_Rho = del_Rho + self%rho(NumberXNodes) - rho_i(NumberXNodes)
                        chargeError = chargeError + (1.0d0 + del_t * (self%J(1) - self%J(NumberXHalfNodes))/del_Rho)**2
                    end if
                END SELECT
            end if
        end do
        chargeError = SQRT(chargeError/k)
    end function getChargeContinuityError

    function getTotalPE(self, world, future) result(res)
        ! Get energy in electric fields, future true, then derive from phi_f, otherwise phi
        ! In 1D is J/m^2
        class(potentialSolver), intent(in) :: self
        type(Domain), intent(in) :: world
        logical :: future
        real(real64) :: res
        if (future) then
            res = 0.5 * eps_0 * SUM(self%EField_f**2 * world%dx_dl)
        else
            res = 0.5 * eps_0 * SUM(self%EField**2 * world%dx_dl)
        end if
    end function getTotalPE

    function getEnergyFromBoundary(self, world, del_t) result(res)
        ! Get energy input into domain from dirichlet boundary
        ! In 1D is J/m^2
        class(potentialSolver), intent(in) :: self
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t
        real(real64) :: res
        res = 0.0d0
    end function getEnergyFromBoundary

    subroutine setRFVoltage(self, world, timeFuture)
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: timeFuture 
    end subroutine setRFVoltage

    subroutine resetVoltage(self)
        class(potentialSolver), intent(in out) :: self
        ! reset current phi value to be future phi value in last time step
        self%EField = self%EField_f
    end subroutine resetVoltage

    subroutine aveRFVoltage(self, accumBool, phi_average, RF_ave, numSteps, world)
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        logical, intent(in) :: accumBool
        real(real64), intent(in out) :: phi_average(NumberXNodes)
        real(real64), intent(in out) :: RF_ave
        integer(int32), intent(in) :: numSteps

    end subroutine aveRFVoltage

    ! ---------------------- read solver input ----------------------------

    subroutine readSolver(GeomFilename, solver, world)
        type(Domain), intent(in) :: world
        type(potentialSolver), intent(in out) :: solver
        character(len=*), intent(in) :: GeomFilename
        integer(int32) :: io, intTemp
        real(real64) :: leftVoltage, rightVoltage, BFieldMag, angle, RF_frequency, realTemp(20)
        real(real64), allocatable :: allocReal(:)

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
            allocate(allocReal(NumberXNodes))
            open(10,file=restartDirectory//"/Phi/phi_0.dat", form='UNFORMATTED', IOSTAT=io)
            read(10, IOSTAT = io) allocReal
            close(10)
            leftVoltage = allocReal(1)
            rightVoltage = allocReal(NumberXNodes)
            deallocate(allocReal)
            open(10,file=restartDirectory//"/"//"InitialConditions.dat", IOSTAT=io)
            read(10, *, IOSTAT = io)
            read(10, *, IOSTAT = io) realTemp(1:15)
            close(10) 
            RF_frequency = realTemp(14)/2.0d0/pi
            if (world%boundaryConditions(1) == 4) then
                leftVoltage = realTemp(15)
            else if (world%boundaryConditions(NumberXNodes) == 4) then
                rightVoltage = realTemp(15)
            end if
            
        end if
        BFieldMag = 0.0d0
        angle = 0.0d0 
        solver = potentialSolver(world, leftVoltage, rightVoltage, BFieldMag, angle, RF_frequency)
        print *, 'Left voltage:', solver%phi(1)
        print *, 'Right voltage:', solver%phi(NumberXNodes)
        print *, 'RF frequency:', solver%RF_rad_frequency/2.0d0/pi
        print *, 'RF_half_amplitude:', solver%RF_half_amplitude
        print *, 'RF_rad_frequency:', solver%RF_rad_frequency
        print *, "BField vector:", solver%BField
        print *, "------------------"
        print *, ""

    end subroutine readSolver


end module mod_potentialSolver