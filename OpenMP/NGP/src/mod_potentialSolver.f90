module mod_potentialSolver
    use iso_fortran_env, only: int32, real64, output_unit
    use constants
    use mod_BasicFunctions
    use mod_domain
    use mod_Scheme
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
        real(real64) :: rho_const, BFieldMag, BField(3), BFieldAngle, RF_rad_frequency, RF_half_amplitude
        real(real64), allocatable :: a_tri(:), b_tri(:), c_tri(:) !for thomas algorithm potential solver, a_tri is lower diagonal, b_tri middle, c_tri upper
        logical :: BFieldBool, RF_bool


    contains
        procedure, public, pass(self) :: solve_tridiag_Poisson
        procedure, public, pass(self) :: solve_tridiag_Ampere
        procedure, public, pass(self) :: getChargeContinuityError
        procedure, public, pass(self) :: getTotalPE
        procedure, public, pass(self) :: setRFVoltage
        procedure, public, pass(self) :: aveRFVoltage
        procedure, public, pass(self) :: getEnergyFromBoundary
        procedure, public, pass(self) :: getError_tridiag_Ampere
        procedure, public, pass(self) :: depositRho
        procedure, public, pass(self) :: makeConstSourceTerm
        procedure, public, pass(self) :: makeHalfTimeEField
        procedure, public, pass(self) :: solveInitialPotential
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
        allocate(self % J(NumberXHalfNodes), self%rho(NumberXNodes), self % phi(NumberXNodes), self % phi_f(NumberXNodes), self%a_tri(NumberXHalfNodes), &
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
        self%phi_f = self%phi 
        self%BFieldMag = BFieldMag
        self%BFieldBool = (BFieldMag /= 0.0d0)
        angle_rad = angle * pi / 180.0d0
        self%BFieldAngle = angle_rad
        self%BField(1) = BFieldMag * COS(angle_rad)
        self%BField(2) = BFieldMag * SIN(angle_rad)
        self%BField(3) = 0.0d0 
        self%EField = 0.0d0
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
                    self % c_tri(i) = 2.0d0/world%dx_dl(i)
                    !self%a_tri(i - leftNodeIdx) = 2.0d0/(world%dx_dl(i-1) + world%dx_dl(i)) / world%dx_dl(i-1)
                    self%b_tri(i) = - (2.0d0/world%dx_dl(i))
                else if (i == NumberXNodes) then
                    self % a_tri(i-1) = 2.0d0/ world%dx_dl(i-1)
                    !self % c_tri(i - leftNodeIdx) = 2.0d0/(world%dx_dl(i-2) + world%dx_dl(i-1))/world%dx_dl(i-1)
                    self%b_tri(i) = - (2.0d0/world%dx_dl(i-1))
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

    subroutine depositRho(self, particleList, world) 
        class(potentialSolver), intent(in out) :: self
        type(Particle), intent(in out) :: particleList(numberChargedParticles)
        type(Domain), intent(in) :: world
        integer(int32) :: i, iThread, rightThreadIndx, leftThreadIndx
        !$OMP parallel private(iThread, i, rightThreadIndx, leftThreadIndx)
        iThread = omp_get_thread_num() + 1 
        leftThreadIndx = world%threadNodeIndx(1,iThread)
        rightThreadIndx = world%threadNodeIndx(2,iThread)
        do i = 1, numberChargedParticles
            particleList(i)%densities(:,iThread) = 0.0d0
            call interpolateParticleToNodes(particleList(i), world, iThread) 
        end do
        !$OMP barrier
        ! Concatenate density array among threads
        particleList(1)%densities(leftThreadIndx:rightThreadIndx, 1) = SUM(particleList(1)%densities(leftThreadIndx:rightThreadIndx, :), DIM=2) * particleList(1)%q_times_wp
        do i = 2, numberChargedParticles
            particleList(1)%densities(leftThreadIndx:rightThreadIndx, 1) = particleList(1)%densities(leftThreadIndx:rightThreadIndx, 1) &
                + SUM(particleList(i)%densities(leftThreadIndx:rightThreadIndx, :), DIM=2) * particleList(i)%q_times_wp
        end do
        !$OMP end parallel
        SELECT CASE (world%boundaryConditions(1))
        CASE(1,4)
            particleList(1)%densities(1, 1) = 0.0d0
        CASE(2)
            particleList(1)%densities(1,1) = 2.0d0 * particleList(1)%densities(1,1)
        CASE(3)
            particleList(1)%densities(1,1) = particleList(1)%densities(1,1) + particleList(1)%densities(NumberXNodes, 1)
        END SELECT
        SELECT CASE (world%boundaryConditions(NumberXNodes))
        CASE(1,4)
            particleList(1)%densities(NumberXNodes, 1) = 0.0d0
        CASE(2)
            particleList(1)%densities(NumberXNodes,1) = 2.0d0 * particleList(1)%densities(NumberXNodes,1)
        CASE(3)
            particleList(1)%densities(NumberXNodes,1) = particleList(1)%densities(1,1)
        END SELECT
        if (world%gridSmoothBool) then
            do i = 1, NumberXNodes
                SELECT CASE(world%boundaryConditions(i))
                CASE(0)
                    self%rho(i) = 0.25d0 * (particleList(1)%densities(i-1, 1) + 2.0d0 * particleList(1)%densities(i, 1) + particleList(1)%densities(i+1, 1))
                CASE(1,4)
                    self%rho(i) = 0.0d0
                CASE(2)
                    if (i ==1) then
                        self%rho(i) = 0.25d0 * (2.0d0 * particleList(1)%densities(1,1) + 2.0d0 * particleList(1)%densities(2,1))
                    else    
                        self%rho(i) = 0.25d0 * (2.0d0 * particleList(1)%densities(NumberXNodes,1) + 2.0d0 * particleList(1)%densities(NumberXHalfNodes,1))
                    end if
                CASE(3)
                    self%rho(i) = 0.25d0 * (particleList(1)%densities(NumberXHalfNodes, 1) + 2.0d0 * particleList(1)%densities(1, 1) + particleList(1)%densities(2, 1))
                END SELECT
            end do
        else
            self%rho = particleList(1)%densities(:, 1)
        end if
    end subroutine depositRho

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
        self%phi_f = self%phi
        !self%EField = (self%phi(1:NumberXHalfNodes) - self%phi(2:NumberXNodes)) / world%dx_dl
    end subroutine solve_tridiag_Poisson

    subroutine solveInitialPotential(self, particleList, world, timeCurrent)
        ! Solve for initial potential
        class(potentialSolver), intent(in out) :: self
        type(Particle), intent(in out) :: particleList(numberChargedParticles)
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: timeCurrent
        call self%depositRho(particleList, world)
        call self%solve_tridiag_Poisson(world, timeCurrent)
        ! Assume only use potential solver once, then need to generate matrix for Div-Ampere

    end subroutine solveInitialPotential

    function getError_tridiag_Poisson(self, world) result(res)
        class(potentialSolver), intent(in) :: self
        type(Domain), intent(in) :: world
        integer(int32) :: i
        real(real64) :: Ax(NumberXNodes), res, tiny
        Ax = triMul(NumberXNodes, self%a_tri, self%c_tri, self%b_tri, self%phi_f)
        res = 0.0d0
        do i = 1, NumberXNodes
            SELECT CASE (world%boundaryConditions(i))
            CASE(0,2)
                if (self%rho(i) + self%rho_const /= 0) then
                    res = res + (Ax(i)*eps_0/(-self%rho(i) - self%rho_const) - 1.0d0)**2
                else
                    res = res + (Ax(i)*eps_0)**2
                end if
            CASE(1,3,4)
                if (self%phi_f(i) /= 0) then
                    res = res + ((Ax(i))/(self.phi_f(i)) - 1.0d0)**2
                else
                    res = res + (Ax(i))**2
                end if
            END SELECT
        end do
        res = SQRT(res/NumberXNodes)

    end function getError_tridiag_Poisson

    subroutine makeConstSourceTerm(self, world)
        ! Store constant source term for ampere into rho array
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        integer(int32) :: i
        do i =1, NumberXNodes
            SELECT CASE (world%boundaryConditions(i))
            CASE(0)
                self%rho(i) = - (self%phi(i) - self%phi(i-1))/world%dx_dl(i-1) + (self%phi(i+1) - self%phi(i))/world%dx_dl(i)
            CASE(1,3,4)
                self%rho(i) = self%phi_f(i)
            CASE(2)
                if (i == 1) then
                    self%rho(i) = -2.0d0 * (self%phi(i) - self%phi(i+1))/world%dx_dl(i)
                else if (i == NumberXNodes) then
                    self%rho(i) = -2.0d0 * (self%phi(i) - self%phi(i-1))/world%dx_dl(i-1)
                end if
            END SELECT
        end do
    end subroutine makeConstSourceTerm

    subroutine solve_tridiag_Ampere(self, world, del_t)
        ! Tridiagonal (Thomas algorithm) solver for Ampere
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t
        integer(int32) :: i !n size dependent on how many points are boundary conditions and thus moved to rhs equation
        real(real64) :: m, d(NumberXNodes), cp(NumberXHalfNodes),dp(NumberXNodes)
        do i =1, NumberXNodes
            SELECT CASE (world%boundaryConditions(i))
            CASE(0)
                d(i) = (self%J(i) - self%J(i-1)) * del_t / eps_0
            CASE(1,3,4)
                d(i) = 0.0d0
            CASE(2)
                if (i == 1) then
                    d(i) = 2.0d0 * (del_t * self%J(i)/eps_0)
                else if (i == NumberXNodes) then
                    d(i) = 2.0d0 * (-del_t * self%J(i-1)/eps_0)
                end if
            END SELECT
            d(i) = d(i) + self%rho(i)
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
    ! initialize phi_f
        !d = self%phi_f
        self%phi_f(NumberXNodes) = dp(NumberXNodes)

    ! solve for x from the vectors c-prime and d-prime
        do i = NumberXHalfNodes, 1, -1
            self%phi_f(i) = dp(i)-cp(i)*self%phi_f(i+1)
        end do
        !self%EField = 0.5d0 * (self%phi(1:NumberXHalfNodes) + self%phi_f(1:NumberXHalfNodes) - self%phi(2:NumberXNodes) - self%phi_f(2:NumberXNodes)) / world%dx_dl
    end subroutine solve_tridiag_Ampere


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
                d(i) = (self%J(i) - self%J(i-1)) * del_t / eps_0 - (self%phi(i) - self%phi(i-1))/world%dx_dl(i-1) + (self%phi(i+1) - self%phi(i))/world%dx_dl(i)
            CASE(1,3,4)
                d(i) = self%phi_f(i)
            CASE(2)
                if (i == 1) then
                    d(i) = 2.0d0 * (del_t * self%J(i)/eps_0 - (self%phi(i) - self%phi(i+1))/world%dx_dl(i))
                else if (i == NumberXNodes) then
                    d(i) = 2.0d0 * (-del_t * self%J(i-1)/eps_0 - (self%phi(i) - self%phi(i-1))/world%dx_dl(i-1))
                end if
            END SELECT
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
                SELECT CASE (world%boundaryConditions(i))
                CASE(0)
                    chargeError = chargeError + (1.0d0 + del_t * (self%J(i) - self%J(i-1))/del_Rho)**2
                    k = k + 1
                CASE(1,4)
                    continue
                CASE(2)
                    if (i == 1) then
                        chargeError = chargeError + (1.0d0 + 2.0d0 * del_t * self%J(1)/del_Rho)**2
                    else
                        chargeError = chargeError + (1.0d0 - 2.0d0 * del_t * self%J(NumberXHalfNodes)/del_Rho)**2
                    end if
                    k = k + 1
                CASE(3)
                    if (i == 1) then
                        chargeError = chargeError + (1.0d0 + del_t * (self%J(1) - self%J(NumberXHalfNodes))/del_Rho)**2
                    end if
                END SELECT
            end if
        end do
        chargeError = SQRT(chargeError/k)
    end function getChargeContinuityError

    subroutine makeHalfTimeEField(self, workArray, world)
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        real(real64), intent(in out) :: workArray(NumberXHalfNodes)
        ! Place temporary field into particle workspace
        workArray = 0.5d0 * (self%phi_f(1:NumberXHalfNodes) + self%phi(1:NumberXHalfNodes) - &
            self%phi_f(2:NumberXNodes) - self%phi(2:NumberXNodes))
        ! Now apply smoothing if on
        if (world%gridSmoothBool) then
            call world%smoothField(workArray, self%EField)
            self%EField = self%EField/world%dx_dl
        else
            self%EField = workArray/world%dx_dl
        end if
    end subroutine makeHalfTimeEField

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

    function getEnergyFromBoundary(self, world, del_t) result(res)
        ! Get energy input into domain from dirichlet boundary
        ! In 1D is J/m^2
        class(potentialSolver), intent(in) :: self
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t
        real(real64) :: res
        res = del_t * self%J(1)
        res = res + eps_0 * ((self%phi_f(1) - self%phi_f(2)) -(self%phi(1) - self%phi(2))) /world%dx_dl(1)
        res = res * (self%phi_f(1) + self%phi(1) - self%phi_f(NumberXNodes) - self%phi(NumberXNodes)) * 0.5d0
    end function getEnergyFromBoundary

    subroutine setRFVoltage(self, world, timeFuture)
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: timeFuture 
        if (world%boundaryConditions(1) == 4) then
            self%phi_f(1) = self%RF_half_amplitude * SIN(self%RF_rad_frequency * (timeFuture))
        else
            self%phi_f(NumberXNodes) = self%RF_half_amplitude * SIN(self%RF_rad_frequency * (timeFuture))
        end if
    end subroutine setRFVoltage

    subroutine resetVoltage(self)
        class(potentialSolver), intent(in out) :: self
        ! reset current phi value to be future phi value in last time step
        self%phi = self%phi_f
    end subroutine resetVoltage

    subroutine aveRFVoltage(self, accumBool, phi_average, RF_ave, numSteps, world)
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        logical, intent(in) :: accumBool
        real(real64), intent(in out) :: phi_average(NumberXNodes)
        real(real64), intent(in out) :: RF_ave
        integer(int32), intent(in) :: numSteps

        if (accumBool) then
            phi_average = phi_average + self%phi_f  
            if (world%boundaryConditions(1) == 4) then
                RF_ave = RF_ave + self%phi_f(1)
            else if (world%boundaryConditions(NumberXNodes) == 4) then
                RF_ave = RF_ave + self%phi_f(NumberXNodes)
            end if
        else
            phi_average = phi_average/numSteps
            self%phi_f = phi_average
            if (self%RF_bool) then
                RF_ave = RF_ave/numSteps
            end if
        end if

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