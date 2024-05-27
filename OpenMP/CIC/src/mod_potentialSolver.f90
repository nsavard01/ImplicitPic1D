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
        real(real64) :: rho_const, RF_rad_frequency, RF_half_amplitude, RF_voltage, boundPhi(2), sourceTermVals(2)
        real(real64), allocatable :: a_tri(:), b_tri(:), c_tri(:) !for thomas algorithm potential solver, a_tri is lower diagonal, b_tri middle, c_tri upper
        integer(int32) :: RF_indx
        logical :: RF_bool


    contains
        procedure, public, pass(self) :: solve_tridiag_Poisson
        procedure, public, pass(self) :: solve_tridiag_Ampere
        procedure, public, pass(self) :: getChargeContinuityError
        procedure, public, pass(self) :: getTotalPE
        procedure, public, pass(self) :: getJdotE
        procedure, public, pass(self) :: setRFVoltage
        procedure, public, pass(self) :: resetVoltage
        procedure, public, pass(self) :: aveRFVoltage
        procedure, public, pass(self) :: getEnergyFromBoundary
        procedure, public, pass(self) :: getError_tridiag_Ampere
        procedure, public, pass(self) :: depositRho
        procedure, public, pass(self) :: makeConstSourceTerm
        procedure, public, pass(self) :: makeHalfTimeEField
        procedure, public, pass(self) :: solveInitialPotential
        procedure, public, pass(self) :: getError_tridiag_Poisson
        procedure, public, pass(self) :: construct_diagMatrix
    end type

    interface potentialSolver
        module procedure :: potentialSolver_constructor
    end interface potentialSolver

contains

    type(potentialSolver) function potentialSolver_constructor(world, leftVoltage, rightVoltage, RF_frequency) result(self)
        ! Construct domain object, initialize grid, dx_dl, and nodeVol.
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: leftVoltage, rightVoltage
        real(real64), intent(in) :: RF_frequency
        allocate(self % J(NumberXHalfNodes), self%rho(NumberXNodes), self % phi(NumberXNodes), self % phi_f(NumberXNodes), self%a_tri(NumberXNodes-1), &
        self%b_tri(NumberXNodes), self%c_tri(NumberXNodes-1), self%EField(NumberXHalfNodes))
        self % rho = 0.0d0
        self % J = 0.0d0
        self % phi = 0.0d0
        self%EField = 0.0d0
        self % rho_const = 0.0d0
        self%phi_f = self%phi 
        self%EField = 0.0d0
        self%RF_half_amplitude = 0.0d0
        self%boundPhi = 0.0d0 
        self%sourceTermVals = 0.0d0

        if (world%boundaryConditions(1)== 1) then
            self%boundPhi(1) = leftVoltage
            self%sourceTermVals(1) = -self%boundPhi(1) * 2.0d0 / world%dx_dl(1)
        else if (world%boundaryConditions(1) == 4) then
            self%RF_half_amplitude = leftVoltage
            self%RF_indx = 1
        end if

        if (world%boundaryConditions(NumberXHalfNodes) == 1) then
            self%boundPhi(2) = rightVoltage
            self%sourceTermVals(2) = -self%boundPhi(2) * 2.0d0 / world%dx_dl(NumberXNodes)
        else if (world%boundaryConditions(NumberXHalfNodes) == 4) then
            self%RF_half_amplitude = rightVoltage
            self%RF_indx = 2
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
            SELECT CASE (world%boundaryConditions(i+1) - world%boundaryConditions(i))
            CASE(0)
                self%b_tri(i) = -2.0d0 * (1.0d0 / (world%dx_dl(i-1) + world%dx_dl(i)) + 1.0d0/(world%dx_dl(i+1) + world%dx_dl(i)))
                self%a_tri(i-1) = 2.0d0 / (world%dx_dl(i-1) + world%dx_dl(i))
                self%c_tri(i) = 2.0d0 / (world%dx_dl(i) + world%dx_dl(i+1))
            CASE(-1,-4, -3)
                self%b_tri(i) = -2.0d0 * (1.0d0 / (world%dx_dl(i+1) + world%dx_dl(i)) + 1.0d0/world%dx_dl(i))
                self%c_tri(i) = 2.0d0 / (world%dx_dl(i+1) + world%dx_dl(i))
                self%sourceTermVals(1) = - self%boundPhi(1) * 2.0d0 / world%dx_dl(1)
            CASE(1,4, 3)
                self%b_tri(i) = -2.0d0 * (1.0d0 / (world%dx_dl(i-1) + world%dx_dl(i)) + 1.0d0/world%dx_dl(i))
                self%a_tri(i-1) = 2.0d0 / (world%dx_dl(i-1) + world%dx_dl(i))
                self%sourceTermVals(2) = - self%boundPhi(2) * 2.0d0 / world%dx_dl(NumberXNodes)
            CASE(-2)
                self%b_tri(i) = -2.0d0 / (world%dx_dl(i+1) + world%dx_dl(i))
                self%c_tri(i) = 2.0d0 / (world%dx_dl(i+1) + world%dx_dl(i))
            CASE(2)
                self%b_tri(i) = -2.0d0 / (world%dx_dl(i-1) + world%dx_dl(i))
                self%a_tri(i-1) = 2.0d0 / (world%dx_dl(i-1) + world%dx_dl(i))
            END SELECT
        end do


    end subroutine construct_diagMatrix

    subroutine depositRho(self, particleList, world) 
        class(potentialSolver), intent(in out) :: self
        type(Particle), intent(in out) :: particleList(numberChargedParticles)
        type(Domain), intent(in) :: world
        integer(int32) :: i, iThread, leftThreadIndx, rightThreadIndx
        !$OMP parallel private(iThread, i, leftThreadIndx, rightThreadIndx)
        iThread = omp_get_thread_num() + 1 
        do i = 1, numberChargedParticles
            particleList(i)%densities(:,iThread) = 0.0d0
            call interpolateParticleToNodes(particleList(i), world, iThread)
        end do
        !$OMP barrier
        leftThreadIndx = world%threadNodeIndx(1,iThread)
        rightThreadIndx = world%threadNodeIndx(2,iThread)
        particleList(1)%densities(leftThreadIndx:rightThreadIndx, 1) = SUM(particleList(1)%densities(leftThreadIndx:rightThreadIndx, :), DIM=2) * particleList(1)%q_times_wp
        do i = 2, numberChargedParticles
            particleList(1)%densities(leftThreadIndx:rightThreadIndx, 1) = particleList(1)%densities(leftThreadIndx:rightThreadIndx, 1) &
                + SUM(particleList(i)%densities(leftThreadIndx:rightThreadIndx, :), DIM=2) * particleList(i)%q_times_wp
        end do
        !$OMP end parallel
        if (world%gridSmoothBool) then
            do i = 1, NumberXNodes
                SELECT CASE(world%boundaryConditions(i+1) - world%boundaryConditions(i))
                CASE(0)
                    self%rho(i) = 0.25d0 * (particleList(1)%densities(i-1, 1) + 2.0d0 * particleList(1)%densities(i, 1) + particleList(1)%densities(i+1, 1))
                CASE(-1,-4)
                    !Dirichlet to left
                    self%rho(i) = 0.25d0 * (particleList(1)%densities(1,1) + particleList(1)%densities(2,1))
                CASE(1,4)
                    ! Dirichlet to right
                    self%rho(i) = 0.25d0 * (particleList(1)%densities(NumberXNodes,1) + particleList(1)%densities(NumberXNodes-1,1))
                CASE(-2)
                    !Neumann to left
                    self%rho(i) = 0.25d0 * (3.0d0 * particleList(1)%densities(1, 1) + particleList(1)%densities(2, 1))
                CASE(2)
                    !Neumann to right
                    self%rho(i) = 0.25d0 * (3.0d0 * particleList(1)%densities(NumberXNodes, 1) + particleList(1)%densities(NumberXNodes-1, 1))
                CASE(-3)
                    !periodic to left
                    self%rho(i) = 0.25d0 * (particleList(1)%densities(NumberXNodes, 1) + 2.0d0 * particleList(1)%densities(1, 1) + particleList(1)%densities(2, 1))
                CASE(3)
                    !periodic to right
                    self%rho(i) = 0.25d0 * (particleList(1)%densities(1, 1) + 2.0d0 * particleList(1)%densities(NumberXNodes, 1) + particleList(1)%densities(NumberXNodes-1, 1))
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
        real(real64) :: m, d(NumberXNodes), cp(NumberXNodes-1),dp(NumberXNodes)
        
        if (self%RF_bool) then
            if (self%RF_indx == 1) then
                self%boundPhi(1) = self%RF_half_amplitude * SIN(self%RF_rad_frequency * timeCurrent)
                self%sourceTermVals(1) = -self%boundPhi(1) * 2.0d0 / world%dx_dl(1)
            else
                self%boundPhi(2) = self%RF_half_amplitude * SIN(self%RF_rad_frequency * timeCurrent)
                self%sourceTermVals(2) = -self%boundPhi(2) * 2.0d0 / world%dx_dl(NumberXNodes)
            end if
        end if
        d = (-self%rho - self%rho_const) / eps_0
        d(1) = d(1) + self%sourceTermVals(1)
        d(NumberXNodes) = d(NumberXNodes) + self%sourceTermVals(2)

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
        real(real64) :: Ax(NumberXNodes), d(NumberXNodes), res
        Ax = triMul(NumberXNodes, self%a_tri, self%c_tri, self%b_tri, self%phi_f)
        d = (-self%rho - self%rho_const) / eps_0
        d(1) = d(1) + self%sourceTermVals(1)
        d(NumberXNodes) = d(NumberXNodes) + self%sourceTermVals(2)
        d = 1.0d0 - Ax/d
        
        !res = Ax*eps_0 - 
        res = SQRT(SUM(d**2)/NumberXNodes)

    end function getError_tridiag_Poisson

    subroutine makeConstSourceTerm(self, world)
        ! Store constant source term for ampere into rho array
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        integer(int32) :: i

        do i =1, NumberXNodes
            SELECT CASE (world%boundaryConditions(i+1) - world%boundaryConditions(i))
            CASE(0)
                self%rho(i) = -(self%phi(i) - self%phi(i+1))/(world%grid(i+1) - world%grid(i)) + (self%phi(i-1) - self%phi(i))/(world%grid(i) - world%grid(i-1))
            CASE(-1,-4)     
                self%rho(i) = -(self%phi(i) - self%phi(i+1))/(world%grid(i+1) - world%grid(i)) + 2.0d0 * (self%boundPhi(1) - self%phi(1))/world%dx_dl(1) + self%sourceTermVals(1)
            CASE(1,4)
                self%rho(i) = -2.0d0 * (self%phi(NumberXNodes) - self%boundPhi(2))/world%dx_dl(NumberXNodes) + (self%phi(i-1) - self%phi(i))/(world%grid(i) - world%grid(i-1)) + self%sourceTermVals(2)
            CASE(-2)
                self%rho(i) = - (self%phi(i) - self%phi(i+1))/(world%grid(i+1) - world%grid(i))
            CASE(2)
                self%rho(i) = (self%phi(i-1) - self%phi(i))/(world%grid(i) - world%grid(i-1))
            CASE(-3)
                self%rho(i) = -(self%phi(i) - self%phi(i+1))/(world%grid(i+1) - world%grid(i)) + 2.0d0 * (self%phi(NumberXNodes) - self%phi(1))/(world%dx_dl(1) + world%dx_dl(NumberXNodes))
            CASE(3)
                self%rho(i) = -2.0d0 * (self%phi(NumberXNodes) - self%phi(1))/(world%dx_dl(1) + world%dx_dl(NumberXNodes)) + (self%phi(i-1) - self%phi(i))/(world%grid(i) - world%grid(i-1))
            END SELECT
        end do
    end subroutine makeConstSourceTerm

    subroutine solve_tridiag_Ampere(self, world, del_t)
        ! Tridiagonal (Thomas algorithm) solver for Ampere
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t
        integer(int32) :: i !n size dependent on how many points are boundary conditions and thus moved to rhs equation
        real(real64) :: d(NumberXNodes),dp(NumberXNodes), cp(NumberXNodes-1), m
        
        d = (self%J(2:NumberXHalfNodes) - self%J(1:NumberXNodes)) * del_t/eps_0 + self%rho
        
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

    function getError_tridiag_Ampere(self, world, del_t) result(res)
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t
        integer(int32) :: i
        real(real64) :: Ax(NumberXNodes), d(NumberXNodes), res(NumberXNodes)
        Ax = triMul(NumberXNodes, self%a_tri, self%c_tri, self%b_tri, self%phi_f)
        d = (self%J(2:NumberXHalfNodes) - self%J(1:NumberXNodes)) * del_t/eps_0 + self%rho
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
                chargeError = chargeError + (1.0d0 + del_t * (self%J(i+1) - self%J(i))/(del_Rho))**2
                k = k + 1
            end if
        end do
        chargeError = SQRT(chargeError/k)
    end function getChargeContinuityError

    subroutine makeHalfTimeEField(self, workArray, world)
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        real(real64), intent(in out) :: workArray(NumberXHalfNodes)
        SELECT CASE (world%boundaryConditions(1))
        CASE(1)
            workArray(1) = (2.0d0 * self%boundPhi(1) - self%phi(1) - self%phi_f(1))  
        CASE(2)
            workArray(1) = 0.0d0
        CASE(3)
            workArray(1) = 0.5d0 * (self%phi(NumberXNodes) + self%phi_f(NumberXNodes) - self%phi(1) - self%phi_f(1))
        CASE(4)
            workArray(1) = (self%boundPhi(1) + self%RF_voltage - self%phi(1) - self%phi_f(1))
        END SELECT
        SELECT CASE (world%boundaryConditions(NumberXHalfNodes))
        CASE(1)
            workArray(NumberXHalfNodes) = self%phi(NumberXNodes) + self%phi_f(NumberXNodes) - 2.0d0 * self%boundPhi(2) 
        CASE(2)
            workArray(NumberXHalfNodes) = 0.0d0
        CASE(3)
            workArray(NumberXHalfNodes) = workArray(1)
        CASE(4)
            workArray(NumberXHalfNodes) = self%phi(NumberXNodes) + self%phi_f(NumberXNodes) - self%boundPhi(2) - self%RF_voltage
        END SELECT
        ! Place temporary field into particle workspace
        workArray(2:NumberXNodes) = 0.5d0 * (self%phi_f(1:NumberXNodes-1) + self%phi(1:NumberXNodes-1) - &
            self%phi_f(2:NumberXNodes) - self%phi(2:NumberXNodes))
        ! Now apply smoothing if on
        if (world%gridSmoothBool) then
            call world%smoothField(workArray, self%EField)
        else
            self%EField = workArray
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
            res = 0.5 * eps_0 * SUM(arrayDiff(self%phi_f, NumberXNodes)**2 / arrayDiff(world%grid, NumberXNodes))
            SELECT CASE (world%boundaryConditions(1))
            CASE(1)
                res = res + eps_0 * (self%boundPhi(1) - self%phi_f(1))**2 / world%dx_dl(1)
            CASE(2)
                continue
            CASE(3)
                res = res + eps_0 * (self%phi_f(NumberXNodes) - self%phi_f(1))**2 / (world%dx_dl(1) + world%dx_dl(NumberXNodes))
            CASE(4)
                res = res + eps_0 * (self%RF_voltage - self%phi_f(1))**2 / world%dx_dl(1)
            END SELECT
            SELECT CASE (world%boundaryConditions(numberXHalfNodes))
            CASE(1)
                res = res + eps_0 * (self%boundPhi(2) - self%phi_f(NumberXNodes))**2 / world%dx_dl(NumberXNodes)
            CASE(2,3)
                continue
            CASE(4)
                res = res + eps_0 * (self%RF_voltage - self%phi_f(NumberXNodes))**2 / world%dx_dl(NumberXNodes)
            END SELECT
        else
            res = 0.5 * eps_0 * SUM(arrayDiff(self%phi, NumberXNodes)**2 / arrayDiff(world%grid, NumberXNodes))
            SELECT CASE (world%boundaryConditions(1))
            CASE(1,4)
                res = res + eps_0 * (self%boundPhi(1) - self%phi(1))**2 / world%dx_dl(1)
            CASE(2)
                continue
            CASE(3)
                res = res + eps_0 * (self%phi(NumberXNodes) - self%phi(1))**2 / ((world%dx_dl(1) + world%dx_dl(NumberXNodes)))
            END SELECT
            SELECT CASE (world%boundaryConditions(numberXHalfNodes))
            CASE(1,4)
                res = res + eps_0 * (self%boundPhi(2) - self%phi(NumberXNodes))**2 / world%dx_dl(NumberXNodes)
            CASE(2,3)
                continue
            END SELECT
        end if
    end function getTotalPE

    function getJdotE(self, world, del_t) result(res)
        ! Get energy in electric fields, future true, then derive from phi_f, otherwise phi
        ! In 1D is J/m^2
        class(potentialSolver), intent(in) :: self
        type(Domain), intent(in) :: world
        real(real64), intent(in) ::del_t
        real(real64) :: res
        
        res = SUM(self%J(2:NumberXNodes) * self%EField(2:NumberXNodes))
        res = res + 0.5d0 * self%J(1) * self%EField(1) + self%J(NumberXHalfNodes) * self%EField(NumberXHalfNodes) * 0.5d0
        res = res * del_t
    end function getJdotE

    function getEnergyFromBoundary(self, world, del_t) result(res)
        ! Get energy input into domain from dirichlet boundary
        ! In 1D is J/m^2
        class(potentialSolver), intent(in) :: self
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t
        real(real64) :: res, leftPhiAve, rightPhiAve
      
        
        SELECT CASE (world%boundaryConditions(1))
        CASE(1)
            leftPhiAve = self%boundPhi(1)
        CASE(2,3)
            leftPhiAve = 0.0
        CASE(4)
            leftPhiAve = 0.5d0 * (self%RF_voltage + self%boundPhi(1))
        END SELECT
        SELECT CASE (world%boundaryConditions(NumberXHalfNodes))
        CASE(1)
            rightPhiAve = self%boundPhi(2)
        CASE(2,3)
            rightPhiAve = 0.0
        CASE(4)
            rightPhiAve = 0.5d0 * (self%RF_voltage + self%boundPhi(2))
        END SELECT
        res = del_t * self%J(2) + eps_0 * ((self%phi_f(1) - self%phi_f(2)) -(self%phi(1) - self%phi(2))) /(world%grid(2) - world%grid(1))
        res = res * (leftPhiAve - rightPhiAve)
        
    end function getEnergyFromBoundary

    subroutine setRFVoltage(self, world, timeFuture)
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: timeFuture 
        
        self%RF_voltage = self%RF_half_amplitude * SIN(self%RF_rad_frequency * (timeFuture))
        if (self%RF_indx == 1) then 
            self%sourceTermVals(1) = -self%RF_voltage * 2.0d0 / world%dx_dl(1)
        else
            self%sourceTermVals(2) = -self%RF_voltage * 2.0d0 / world%dx_dl(NumberXNodes)
        end if

    end subroutine setRFVoltage

    subroutine resetVoltage(self)
        class(potentialSolver), intent(in out) :: self
        self%phi = self%phi_f
        if (self%RF_bool) then
            self%boundPhi(self%RF_indx) = self%RF_voltage
        end if
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
            ! average RF voltage
            if (self%RF_bool) then
                RF_ave = RF_ave + self%RF_voltage
            end if
        else
            phi_average = phi_average/numSteps
            self%phi_f = phi_average
            if (self%RF_bool) then
                RF_ave = RF_ave/numSteps
                if (world%boundaryConditions(1) == 4) then 
                    self%sourceTermVals(1) = -RF_ave * 2.0d0 / world%dx_dl(1)
                else if (world%boundaryConditions(NumberXHalfNodes) == 4) then
                    self%sourceTermVals(2) = -RF_ave * 2.0d0 / world%dx_dl(NumberXNodes)
                end if
            end if
        end if

    end subroutine aveRFVoltage

    ! ------------------ read solver input ---------------------------

    subroutine readSolver(GeomFilename, solver, world)
        type(Domain), intent(in) :: world
        type(potentialSolver), intent(in out) :: solver
        character(len=*), intent(in) :: GeomFilename
        integer(int32) :: io
        real(real64) :: leftVoltage, rightVoltage, BFieldMag, angle, RF_frequency, realTemp(20)

        print *, ""
        print *, "Reading solver inputs:"
        print *, "------------------"
        if (.not. restartBool) then
            open(10,file='../InputData/'//GeomFilename)
        else
            open(10,file=restartDirectory//'/InputData/'//GeomFilename)
        end if
        read(10, *, IOSTAT = io)
        read(10, *, IOSTAT = io)
        read(10, *, IOSTAT = io)
        read(10, *, IOSTAT = io)
        read(10, *, IOSTAT = io)
        read(10, *, IOSTAT = io) leftVoltage, rightVoltage
        read(10, *, IOSTAT = io) RF_frequency
        close(10)
        solver = potentialSolver(world, leftVoltage, rightVoltage, RF_frequency)
        print *, 'Left Voltage val:', solver%boundPhi(1)
        print *, 'Right Voltage val:', solver%boundPhi(2)
        print *, 'RF_bool:', solver%RF_bool
        print *, 'RF frequency:', solver%RF_rad_frequency/2.0d0/pi
        print *, 'RF_half_amplitude:', solver%RF_half_amplitude
        print *, 'RF_rad_frequency:', solver%RF_rad_frequency
        print *, "------------------"
        print *, ""

    end subroutine readSolver

end module mod_potentialSolver