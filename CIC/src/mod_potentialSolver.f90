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
    public :: potentialSolver

    type :: potentialSolver
        real(real64), allocatable :: phi(:), J(:), rho(:), phi_f(:), particleChargeLoss(:,:) !phi_f is final phi, will likely need to store two arrays for phi, can't be avoided
        real(real64) :: energyError, chargeError, particleEnergyLoss
        real(real64) :: coeff_left, coeff_right ! these are coefficients (from world dimensions) needed with phi_left and phi_right in rhs of matrix equation
        real(real64), allocatable :: a_tri(:), b_tri(:), c_tri(:) !for thomas algorithm potential solver, a_tri is lower diagonal, b_tri middle, c_tri upper


    contains
        procedure, public, pass(self) :: depositRho
        procedure, public, pass(self) :: solve_tridiag_Poisson
        procedure, public, pass(self) :: solve_tridiag_Ampere
        procedure, public, pass(self) :: getEField
        procedure, public, pass(self) :: getEFieldLogical
        procedure, public, pass(self) :: getEFieldPeriodic
        procedure, public, pass(self) :: getEFieldDirichlet
        procedure, public, pass(self) :: getTotalPE
        procedure, public, pass(self) :: depositJ
        procedure, public, pass(self) :: getError_tridiag_Ampere
        procedure, public, pass(self) :: construct_diagMatrix_Ampere
        ! procedure, public, pass(self) :: solveDivAmperePicard
        ! procedure, public, pass(self) :: adaptiveSolveDivAmperePicard
        procedure, public, pass(self) :: solveInitialPotential
        procedure, public, pass(self) :: getDelTauInitialSubStep
        procedure, public, pass(self) :: getDelTauSubStep
        procedure, public, pass(self) :: getDelTauInitialSubStepPeriodic
        procedure, public, pass(self) :: getDelTauSubStepPeriodic
        procedure, public, pass(self) :: getDelTauInitialSubStepDirichlet
        procedure, public, pass(self) :: getDelTauSubStepDirichlet
        procedure, public, pass(self) :: picardIterParticles
        procedure, public, pass(self) :: analyticalParticleMover
        procedure, public, pass(self) :: analyticalParticleMoverPeriodic
        procedure, public, pass(self) :: depositJSubStep
        procedure, public, pass(self) :: depositJSubStepPeriodic
        procedure, public, pass(self) :: depositJSubStepDirichlet
        procedure, public, pass(self) :: moveParticles
        procedure, private, pass(self) :: construct_diagMatrix
        ! procedure, public, pass(self) :: adaptiveSolveDivAmpereAnderson
        ! procedure, public, pass(self) :: solveDivAmpereAnderson
    end type

    interface potentialSolver
        module procedure :: potentialSolver_constructor
    end interface potentialSolver
   
contains

    type(potentialSolver) function potentialSolver_constructor(world, leftVoltage, rightVoltage) result(self)
        ! Construct domain object, initialize grid, dx_dl, and dx_dl.
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: leftVoltage, rightVoltage
        allocate(self % J(NumberXNodes-1), self % rho(NumberXNodes), self % phi(NumberXNodes), self % phi_f(NumberXNodes), self%a_tri(NumberXNodes-3), &
        self%b_tri(NumberXNodes-2), self%c_tri(NumberXNodes-3), self%particleChargeLoss(2, numberChargedParticles))
        call construct_diagMatrix(self, world)
        self % rho = 0
        self % J = 0
        self % phi = 0
        self % coeff_left = 2/(world%dx_dl(1) + world%dx_dl(2))/world%dx_dl(2)
        self % coeff_right = 2/(world%dx_dl(NumberXNodes-1) + world%dx_dl(NumberXNodes))/world%dx_dl(NumberXNodes-1)
        self%particleEnergyLoss = 0.0d0
        self%particleChargeLoss = 0.0d0
        self%energyError = 0.0d0
        self%chargeError = 0.0d0
        if (world%boundaryConditions(1) == 1) self%phi(1) = leftVoltage
        if (world%boundaryConditions(NumberXNodes) == 1) self%phi(NumberXNodes) = rightVoltage
        if (world%boundaryConditions(1) == 3) then
            self%phi(1) = leftVoltage
            self%phi(NumberXNodes) = leftVoltage
        end if
        self%phi_f = self%phi  


    end function potentialSolver_constructor


    subroutine depositRho(self, particleList, world) 
        class(potentialSolver), intent(in out) :: self
        type(Particle), intent(in) :: particleList(:)
        type(Domain), intent(in) :: world
        integer(int32) :: i, j, l_center
        real(real64) :: d
        self % rho = 0.0d0
        do i=1, numberChargedParticles
            do j = 1, particleList(i)%N_p
                l_center = NINT(particleList(i)%phaseSpace(1, j))
                d = particleList(i)%phaseSpace(1, j) - l_center
                if (world%boundaryConditions(l_center) == 0) then
                    ! Inside domain
                    self % rho(l_center) = self % rho(l_center) + particleList(i)%q * particleList(i)%w_p * (0.75 - d**2)
                    self % rho(l_center + 1) = self % rho(l_center + 1) + particleList(i)%q * particleList(i)%w_p * 0.5d0 * (0.5d0 + d)**2
                    self % rho(l_center - 1) = self % rho(l_center - 1) + particleList(i)%q * particleList(i)%w_p * 0.5d0 * (0.5d0 - d)**2
                else if (world%boundaryConditions(l_center) == 1) then
                    !Dirichlet
                    self % rho(l_center) = self % rho(l_center) + particleList(i)%q * particleList(i)%w_p * (1.0d0-ABS(d))
                    self % rho(l_center + INT(SIGN(1.0, d))) = self % rho(l_center + INT(SIGN(1.0, d))) + particleList(i)%q * particleList(i)%w_p * ABS(d)
                
                else if (world%boundaryConditions(l_center) == 3) then
                    ! Periodic
                    self % rho(l_center) = self % rho(l_center) + particleList(i)%q * particleList(i)%w_p * (0.75 - d**2)
                    ! towards domain
                    self % rho(l_center+INT(SIGN(1.0, d))) = self % rho(l_center+INT(SIGN(1.0, d))) + particleList(i)%q * particleList(i)%w_p * 0.5d0 * (0.5d0 + ABS(d))**2
                    ! across periodic boundary
                    self % rho(MODULO(l_center-2*INT(SIGN(1.0, d)),NumberXNodes)) = self % rho(MODULO(l_center-2*INT(SIGN(1.0, d)),NumberXNodes)) + particleList(i)%q * particleList(i)%w_p * 0.5d0 * (0.5d0 - ABS(d))**2
                end if
            end do
        end do
        if (world%boundaryConditions(1) == 3) then
            self%rho(1) = self%rho(1) + self%rho(NumberXNodes)
            self%rho(NumberXNodes) = self%rho(1)
        end if
        self % rho = self % rho / world%dx_dl
    end subroutine depositRho

    subroutine construct_diagMatrix(self, world)
        ! construct diagonal components for thomas algorithm
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        self % a_tri = 2.0d0/(world%dx_dl(2:NumberXNodes-2) + world%dx_dl(3:NumberXNodes-1)) / world%dx_dl(3:NumberXNodes-1)
        self % c_tri = 2.0d0/(world%dx_dl(3:NumberXNodes-1) + world%dx_dl(2:NumberXNodes-2))/world%dx_dl(2:NumberXNodes-2)
        self % b_tri = -2.0d0/(world%dx_dl(1:NumberXNodes-2) + world%dx_dl(2:NumberXNodes-1))/world%dx_dl(2:NumberXNodes-1) - 2.0d0/(world%dx_dl(2:NumberXNodes-1) + world%dx_dl(3:))/world%dx_dl(2:NumberXNodes-1)

    end subroutine construct_diagMatrix

    subroutine construct_diagMatrix_Ampere(self, world)
        ! construct diagonal components for thomas algorithm, for Ampere (after initial Poisson)
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        self % a_tri = -2.0d0/(world%dx_dl(2:NumberXNodes-2) + world%dx_dl(3:NumberXNodes-1))
        self % c_tri = -2.0d0/(world%dx_dl(2:NumberXNodes-2) + world%dx_dl(3:NumberXNodes-1))
        self % b_tri = 2.0d0/(world%dx_dl(2:NumberXNodes-1) + world%dx_dl(1:NumberXNodes-2)) + 2.0d0/(world%dx_dl(2:NumberXNodes-1) + world%dx_dl(3:NumberXNodes))
        self % coeff_left = 2.0d0/(world%dx_dl(1) + world%dx_dl(2))
        self % coeff_right = 2.0d0/(world%dx_dl(NumberXNodes) + world%dx_dl(NumberXNodes-1))

    end subroutine construct_diagMatrix_Ampere

    subroutine solve_tridiag_Poisson(self)
        ! Tridiagonal (Thomas algorithm) solver for initial Poisson
        class(potentialSolver), intent(in out) :: self
        integer(int32) :: i, n !n size dependent on how many points are boundary conditions and thus moved to rhs equation
        real(real64) :: m, d(size(self%phi) - 2), cp(size(self%phi) - 3),dp(size(self%phi) - 2)
        n = NumberXNodes - 2

        d = -self%rho(2:NumberXNodes-1) / eps_0
        d(1) = d(1) - 2 * self%phi(1) * self%coeff_left
        d(n) = d(n) - 2 * self%phi(NumberXNodes) * self%coeff_right
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
    ! initialize x
        self%phi(2:n+1) = dp(n)
    ! solve for x from the vectors c-prime and d-prime
        do i = n-1, 1, -1
            self%phi(i+1) = dp(i)-cp(i)*self%phi(i+2)
        end do
        self%phi_f = self%phi

    end subroutine solve_tridiag_Poisson

    pure function singleRho(l_p, w_p, q, world) result(rho)
        ! for diagnostic in substep routine
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: l_p, w_p, q
        real(real64) :: rho(NumberXNodes), d
        integer(int32) :: l_center
        rho = 0.0d0
        l_center = NINT(l_p)
        d = l_p - l_center
        if (world%boundaryConditions(l_center) == 0) then
            ! Inside domain
            rho(l_center) = rho(l_center) + q * w_p * (0.75 - d**2)
            rho(l_center + 1) = rho(l_center + 1) + q * w_p * 0.5d0 * (0.5d0 + d)**2
            rho(l_center - 1) = rho(l_center - 1) + q * w_p * 0.5d0 * (0.5d0 - d)**2
        else if (world%boundaryConditions(l_center) == 1) then
            !Dirichlet
            rho(l_center) = rho(l_center) + q * w_p * (1.0d0-ABS(d))
            rho(l_center + INT(SIGN(1.0, d))) = rho(l_center + INT(SIGN(1.0, d))) + q * w_p * ABS(d)
        
        else if (world%boundaryConditions(l_center) == 3) then
            ! Periodic
            rho(l_center) = rho(l_center) + q * w_p * (0.75 - d**2)
            ! towards domain
            rho(l_center+INT(SIGN(1.0, d))) = rho(l_center+INT(SIGN(1.0, d))) + q * w_p * 0.5d0 * (0.5d0 + d)**2
            ! across periodic boundary
            rho(MOD(l_center-2*INT(SIGN(1.0, d)),NumberXNodes)) = rho(MOD(l_center-2*INT(SIGN(1.0, d)),NumberXNodes)) + q * w_p * 0.5d0 * (0.5d0 - d)**2
        end if
        rho = rho / world%dx_dl
    end function singleRho

    subroutine singleJ(J, l_half, v_half, del_tau, del_t, w_p, q, world)
        ! for diagnostic in substep routine
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: l_half, v_half, w_p, q, del_tau, del_t
        real(real64), intent(in out) :: J(NumberXNodes-1)
        real(real64) :: d
        if (world%boundaryConditions(NINT(l_half)) == 0) then
            d = (l_half) - REAL(NINT(l_half), kind = real64) + 0.5d0
            J(NINT(l_half)-1) = J(NINT(l_half)-1) + (1.0d0 - d) * w_p * q * v_half*del_tau/world%dx_dl(NINT(l_half))/del_t
            J(NINT(l_half)) = J(NINT(l_half)) + d * w_p * q * v_half*del_tau/world%dx_dl(NINT(l_half))/del_t
        else if (world%boundaryConditions(NINT(l_half)) == 1) then
            J(INT(l_half)) = J(INT(l_half)) + w_p * q * v_half*del_tau/world%dx_dl(NINT(l_half))/del_t
        end if  
    end subroutine singleJ

    subroutine solve_tridiag_Ampere(self, world, del_t)
        ! Tridiagonal (Thomas algorithm) solver for Ampere
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t
        integer(int32) :: i, n !n size dependent on how many points are inside (not boundary), so how many in rhs equation
        real(real64) :: m, d(NumberXNodes - 2), cp(NumberXNodes - 3),dp(NumberXNodes- 2)
        n = NumberXNodes - 2

        d = (-self%J(2:) + self%J(1:n)) * del_t / eps_0 &
        + arrayDiff(self%phi(1:n+1), n+1)*2.0d0/(world%dx_dl(1:n) + world%dx_dl(2:n+1)) &
        - arrayDiff(self%phi(2:), n+1)*2.0d0/(world%dx_dl(3:n+2) + world%dx_dl(2:n+1))
        d(1) = d(1) + self%phi(1) * self%coeff_left
        d(n) = d(n) + self%phi(NumberXNodes) * self%coeff_right
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
        real(real64) :: Ax(NumberXNodes-2), d(NumberXNodes-2), res
        Ax(1) = self%b_tri(1)*self%phi_f(2) + self%c_tri(1) * self%phi_f(3)
        do i=2, NumberXNodes-3
            Ax(i) = self%b_tri(i)*self%phi_f(i+1) + self%c_tri(i) * self%phi_f(i+2) + self%a_tri(i-1) * self%phi_f(i)
        end do
        Ax(NumberXNodes-2) = self%b_tri(NumberXNodes-2)*self%phi_f(NumberXNodes-1) + self%a_tri(NumberXNodes-3) * self%phi_f(NumberXNodes-2)
        d = (-self%J(2:) + self%J(1:NumberXNodes-2)) * del_t / eps_0 &
        + arrayDiff(self%phi(1:NumberXNodes-1), NumberXNodes-1)*2.0d0/(world%dx_dl(1:NumberXNodes-2) + world%dx_dl(2:NumberXNodes-1)) &
        - arrayDiff(self%phi(2:), NumberXNodes-1)*2.0d0/(world%dx_dl(3:NumberXNodes) + world%dx_dl(2:NumberXNodes-1))
        d(1) = d(1) + self%phi(1) * self%coeff_left
        d(NumberXNodes-2) = d(NumberXNodes-2) + self%phi(NumberXNodes) * self%coeff_right
        !res = SQRT(SUM(((Ax- d)/self%minEField)**2)/(NumberXNodes-2))
        res = SQRT(SUM((Ax- d)**2))

    end function getError_tridiag_Ampere

    function getTotalPE(self, world, future) result(res)
        ! Get energy in electric fields, future true, then derive from phi_f, otherwise phi
        ! In 1D is J/m^2
        class(potentialSolver), intent(in) :: self
        type(Domain), intent(in) :: world
        logical :: future
        real(real64) :: res
        if (future) then
            res = 0.5 * eps_0 * SUM(arrayDiff(self%phi_f, NumberXNodes)**2 * 2.0d0 / (world%dx_dl(1:NumberXNodes-1) + world%dx_dl(2:NumberXNodes)))
        else
            res = 0.5 * eps_0 * SUM(arrayDiff(self%phi, NumberXNodes)**2 * 2.0d0 / (world%dx_dl(1:NumberXNodes-1) + world%dx_dl(2:NumberXNodes)))
        end if
    end function getTotalPE


    


    !----------------- Particle mover/sub-stepping procedures -------------------------------------------------------

    pure function getEField(self, l_p, l_cell, world) result(EField)
        !return EField of particle at logical position l_p (this is half distance), per particle since particle mover loop will be per particle
        class(potentialSolver), intent(in) :: self
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: l_p
        integer(int32), intent(in) :: l_cell
        real(real64) :: EField, d
        d = l_p - real(l_cell, kind = real64) + 0.5d0
        EField = ((self%phi_f(l_cell) + self%phi(l_cell) - self%phi(l_cell+1) - self%phi_f(l_cell + 1)) * d/2.0d0 +  &
        (self%phi_f(l_cell-1) + self%phi(l_cell-1) - self%phi(l_cell) - self%phi_f(l_cell)) * (1.0d0 - d)/2.0d0)/world%dx_dl(l_cell)
    end function getEField

    pure function getEFieldPeriodic(self, l_p, l_cell, world) result(EField)
        !return EField of particle at logical position l_p (this is half distance), per particle since particle mover loop will be per particle
        class(potentialSolver), intent(in) :: self
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: l_p
        integer(int32), intent(in) :: l_cell
        real(real64) :: EField, d
        d = l_p - real(l_cell, kind = real64) + 0.5d0
        EField = ((self%phi_f(1) + self%phi(1) - self%phi(2) - self%phi_f(2)) * d/2.0d0 +  &
        (self%phi_f(NumberXNodes-1) + self%phi(NumberXNodes-1) - self%phi(1) - self%phi_f(1)) * (1.0d0 - d)/2.0d0)/world%dx_dl(l_cell)
    end function getEFieldPeriodic

    pure function getEFieldLogical(self, l_p) result(EField)
        !return EField of particle at logical position l_p (this is half distance), per particle since particle mover loop will be per particle
        class(potentialSolver), intent(in) :: self
        real(real64), intent(in) :: l_p
        integer(int32) :: l_cell
        real(real64) :: EField, d
        l_cell = NINT(l_p)
        d = l_p - real(l_cell, kind = real64) + 0.5d0
        EField = (self%phi_f(l_cell) + self%phi(l_cell) - self%phi(l_cell+1) - self%phi_f(l_cell + 1)) * d/2.0d0 +  &
        (self%phi_f(l_cell-1) + self%phi(l_cell-1) - self%phi(l_cell) - self%phi_f(l_cell)) * (1.0d0 - d)/2.0d0
    end function getEFieldLogical

    pure function getEFieldDirichlet(self, l_p, l_cell, world) result(EField)
        !return EField of particle at logical position l_p (this is half distance), per particle since particle mover loop will be per particle
        class(potentialSolver), intent(in) :: self
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: l_p
        integer(int32), intent(in) :: l_cell
        real(real64) :: EField, d
        d = l_p - real(l_cell, kind = real64)
        EField = SIGN(1.0, d) * (self%phi_f(l_cell) + self%phi(l_cell) - self%phi(l_cell+ INT(SIGN(1.0, d))) - self%phi_f(l_cell + INT(SIGN(1.0, d))))/world%dx_dl(l_cell)/2.0d0
    end function getEFieldDirichlet

    subroutine getDelTauInitialSubStep(self, world, part, l_sub, l_f, v_sub, v_f, del_tau, del_t, l_alongV, l_awayV, l_cell, a, c)
        class(potentialSolver), intent(in) :: self
        type(Domain), intent(in) :: world
        type(Particle), intent(in) :: part
        real(real64), intent(in out) :: l_sub, l_f, v_sub, v_f, del_tau, l_alongV, l_awayV, a, c
        real(real64), intent(in) :: del_t
        integer(int32), intent(in) :: l_cell
        real(real64) :: del_tau_tmp
        del_tau = del_t
        a = (part%q / part%mass / 2.0d0) * self%getEField((l_sub + l_alongV)/2.0d0, l_cell, world)
        ! Particle first between nodes, so solve quadratic for that particle depending on conditions
        if (v_sub/=0.0d0) then ! make first case, since pretty much always likely to be the case (could just not have, assume always field, never have to check)
            c = (l_sub - l_alongV) * world%dx_dl(l_cell)
            if (a*v_sub > 0) then
                ! velocity and acceleration in same direction
                del_tau = 2.0d0 * ABS(c)/(ABS(v_sub) + SQRT(v_sub**2 - 4.0d0*a*c))
                l_f = l_alongV
                if (del_tau <= 0.0d0) then
                    stop "Have issue with del_tau for v,a in same direction"
                end if
            else if (v_sub**2 - 4.0d0*a*c > 0.0d0) then
                ! v and a opposite direction, but particle can still reach boundary along v
                del_tau = 2.0d0 * ABS(c)/(ABS(v_sub) + SQRT(v_sub**2 - 4.0d0*a*c))
                l_f = l_alongV
                if (del_tau <= 0.0d0) then
                    print *, "l_sub is:", l_sub
                    print *, "v_sub is:", v_sub
                    print *, "a is:", a
                    print *, "In initial sub-step routine"
                    stop "Have issue with del_tau for v,a in opposite direction, but still reach boundary along v"
                end if
            end if
            if (del_tau >= del_t) then
                ! boundary opposite direction of v
                a = (part%q / part%mass / 2.0d0) * self%getEField((l_sub + l_awayV)/2.0d0, l_cell, world)
                if (a*v_sub < 0) then
                    c = (l_sub - l_awayV) * world%dx_dl(l_cell)
                    del_tau_tmp = 2.0d0 * ABS(c)/(SQRT(v_sub**2 - 4.0d0*a*c) - ABS(v_sub))
                    if (del_tau_tmp < del_tau) then
                        ! If del_tau isn't reduced, then want to keep saved l_f since might make picard iteration more stable as initial condition
                        del_tau = del_tau_tmp
                        l_f = l_awayV
                    end if
                    if (del_tau <= 0) then
                        stop "Have issue with del_tau for v,a in opposite direction, boundary opposite v"
                    end if
                end if  
            end if
        else
            ! only in direction of field: USE l_alongV AS BOUNDARY ALONG DIRECTION OF a SINCE VELOCITY = 0!!!
            if (world%boundaryConditions(l_cell) == 0) then
                l_alongV = real(l_cell, kind = real64) + SIGN(0.5d0, a)
            else if (world%boundaryConditions(l_cell) == 1) then
                !Dirichlet
                l_alongV = NINT(2.0d0*l_sub + SIGN(0.5d0, a))/2.0d0
            else if (world%boundaryConditions(l_cell) == 1) then
                l_alongV = real(l_cell, kind = real64) + SIGN(0.5d0, a)
            end if
            c = (l_sub - l_alongV) * world%dx_dl(l_cell)
            del_tau = SQRT(-c/a)
            l_f = l_alongV
            if (del_tau <= 0.0d0) then
                stop "Have issue with del_tau for v = 0"
            end if
        end if

    end subroutine getDelTauInitialSubstep

    subroutine getDelTauSubStep(self, world, part, l_sub, l_f, v_sub, v_f, del_tau, timePassed, del_t, l_alongV, l_cell, a, c)
        class(potentialSolver), intent(in) :: self
        type(Domain), intent(in) :: world
        type(Particle), intent(in) :: part
        real(real64), intent(in out) :: l_sub, l_f, v_sub, v_f, del_tau, l_alongV, a, c
        real(real64), intent(in) :: del_t, timePassed
        integer(int32), intent(in) :: l_cell
        real(real64) :: del_tau_tmp
        del_tau = del_t - timePassed
        ! get index cell where field and dx_dl is evaluated
        a = (part%q / part%mass / 2.0d0) * self%getEField((l_sub + l_alongV)/2.0d0, l_cell, world)
        ! Particle first between nodes, so solve quadratic for that particle depending on conditions
        if ((a/=0.0d0)) then ! make first case, since pretty much always likely to be the case (could just not have, assume always field exists, never have to check)
            c = (l_sub - l_alongV) * world%dx_dl(l_cell)
            l_f = l_alongV
            if (a*v_sub > 0) then
                ! velocity and acceleration in same direction
                del_tau = 2.0d0 * ABS(c)/(ABS(v_sub) + SQRT(v_sub**2 - 4.0d0*a*c))
                l_f = l_alongV
                if (del_tau <= 0.0d0) then
                    stop "Have issue with del_tau for v,a in same direction"
                end if
            else if (v_sub**2 - 4.0d0*a*c > 0.0d0) then
                ! v and a opposite direction, but particle can still reach boundary along v
                del_tau = 2.0d0 * ABS(c)/(ABS(v_sub) + SQRT(v_sub**2 - 4.0d0*a*c))
                l_f = l_alongV
                if (del_tau <= 0.0d0) then
                    print *, "In regular substep routine"
                    stop "Have issue with del_tau for v,a in opposite direction, but still reach boundary along v"
                end if
            end if
            if (del_tau >= del_t-timePassed) then
                ! For particles returning to boundary, make guess for maximum distance could go before v_f = 0 in half sub-step
                if (ABS((del_t - timePassed) *v_sub/4.0d0/world%dx_dl(l_cell)) < 1.0) then
                    a = (part%q / part%mass / 2.0d0) * self%getEField(l_sub + (del_t - timePassed) *v_sub/8.0d0/world%dx_dl(l_cell), l_cell, world)
                    if (a*v_sub < 0) then
                        del_tau_tmp = ABS(v_sub/a)
                        if (del_tau_tmp < del_tau) then
                            del_tau = del_tau_tmp
                            l_f = l_sub
                        end if
                        if (del_tau <= 0) then
                            stop "Have issue with del_tau for v,a in opposite direction, boundary opposite v"
                        end if
                    end if
                end if
                ! try for minimum distance particle can go
                if (del_tau >= del_t - timePassed) then
                    a = (part%q / part%mass / 2.0d0) * self%getEField(l_sub + SIGN(1.0d-8, v_sub), l_cell, world)
                    if (a*v_sub < 0) then
                        del_tau_tmp = ABS(v_sub/a)
                        if (del_tau_tmp < del_tau) then
                            del_tau = del_tau_tmp
                            l_f = l_sub
                        end if
                        if (del_tau <= 0) then
                            stop "Have issue with del_tau for v,a in opposite direction, boundary opposite v, second try close to boundary"
                        end if
                    end if
                end if
            end if
        else
            !Free particle drift
            del_tau = (l_alongV - l_sub) * world%dx_dl(l_cell)/v_sub
            v_f = (l_alongV - l_sub) * world%dx_dl(l_cell) / del_tau
            l_f = l_alongV
            if (del_tau <= 0.0d0) then
                stop "Have issue with del_tau for a = 0"
            end if
        end if

    end subroutine getDelTauSubstep

    subroutine getDelTauInitialSubStepPeriodic(self, world, part, l_sub, l_f, v_sub, v_f, del_tau, del_t, l_alongV, l_awayV, l_cell, a, c)
        class(potentialSolver), intent(in) :: self
        type(Domain), intent(in) :: world
        type(Particle), intent(in) :: part
        real(real64), intent(in out) :: l_sub, l_f, v_sub, v_f, del_tau, l_alongV, l_awayV, a, c
        real(real64), intent(in) :: del_t
        integer(int32), intent(in) :: l_cell
        real(real64) :: del_tau_tmp
        del_tau = del_t
        a = (part%q / part%mass / 2.0d0) * self%getEFieldPeriodic((l_sub + l_alongV)/2.0d0, l_cell, world)
        ! Particle first between nodes, so solve quadratic for that particle depending on conditions
        if ((v_sub/=0.0d0)) then ! make first case, since pretty much always likely to be the case (could just not have, assume always field, never have to check)
            c = (l_sub - l_alongV) * world%dx_dl(l_cell)
            if (a*v_sub > 0) then
                ! velocity and acceleration in same direction
                del_tau = 2.0d0 * ABS(c)/(ABS(v_sub) + SQRT(v_sub**2 - 4.0d0*a*c))
                l_f = l_alongV
                if (del_tau <= 0.0d0) then
                    stop "Have issue with del_tau for v,a in same direction"
                end if
            else if (v_sub**2 - 4.0d0*a*c > 0.0d0) then
                ! v and a opposite direction, but particle can still reach boundary along v
                del_tau = 2.0d0 * ABS(c)/(ABS(v_sub) + SQRT(v_sub**2 - 4.0d0*a*c))
                l_f = l_alongV
                if (del_tau <= 0.0d0) then
                    print *, "l_sub is:", l_sub
                    print *, "v_sub is:", v_sub
                    print *, "a is:", a
                    print *, "In initial sub-step routine"
                    stop "Have issue with del_tau for v,a in opposite direction, but still reach boundary along v"
                end if
            end if
            if (del_tau >= del_t) then
                ! boundary opposite direction of v
                a = (part%q / part%mass / 2.0d0) * self%getEFieldPeriodic((l_sub + l_awayV)/2.0d0, l_cell, world)
                if (a*v_sub < 0) then
                    c = (l_sub - l_awayV) * world%dx_dl(l_cell)
                    del_tau_tmp = 2.0d0 * ABS(c)/(SQRT(v_sub**2 - 4.0d0*a*c) - ABS(v_sub))
                    if (del_tau_tmp < del_tau) then
                        ! If del_tau isn't reduced, then want to keep saved l_f since might make picard iteration more stable as initial condition
                        del_tau = del_tau_tmp
                        l_f = l_awayV
                    end if
                    if (del_tau <= 0) then
                        stop "Have issue with del_tau for v,a in opposite direction, boundary opposite v"
                    end if
                end if  
            end if

        else
            ! only in direction of field: USE l_alongV AS BOUNDARY ALONG DIRECTION OF a SINCE VELOCITY = 0!!!
            if (world%boundaryConditions(l_cell) == 0) then
                l_alongV = real(l_cell, kind = real64) + SIGN(0.5d0, a)
            else if (world%boundaryConditions(l_cell) == 1) then
                !Dirichlet
                l_alongV = NINT(2.0d0*l_sub + SIGN(0.5d0, a))/2.0d0
            else if (world%boundaryConditions(l_cell) == 3) then
                l_alongV = real(l_cell, kind = real64) + SIGN(0.5d0, a)
            end if
            c = (l_sub - l_alongV) * world%dx_dl(l_cell)
            del_tau = SQRT(-c/a)
            l_f = l_alongV
            if (del_tau <= 0.0d0) then
                stop "Have issue with del_tau for v = 0"
            end if
        end if

    end subroutine getDelTauInitialSubstepPeriodic

    subroutine getDelTauSubStepPeriodic(self, world, part, l_sub, l_f, v_sub, v_f, del_tau, timePassed, del_t, l_alongV, l_cell, a, c)
        class(potentialSolver), intent(in) :: self
        type(Domain), intent(in) :: world
        type(Particle), intent(in) :: part
        real(real64), intent(in out) :: l_sub, l_f, v_sub, v_f, del_tau, l_alongV, a, c
        real(real64), intent(in) :: del_t, timePassed
        integer(int32), intent(in) :: l_cell
        real(real64) :: del_tau_tmp
        del_tau = del_t - timePassed
        ! get index cell where field and dx_dl is evaluated
        a = (part%q / part%mass / 2.0d0) * self%getEFieldPeriodic((l_sub + l_alongV)/2.0d0, l_cell, world)
        ! Particle first between nodes, so solve quadratic for that particle depending on conditions
        if ((a/=0.0d0)) then ! make first case, since pretty much always likely to be the case (could just not have, assume always field exists, never have to check)
            c = (l_sub - l_alongV) * world%dx_dl(l_cell)
            l_f = l_alongV
            if (a*v_sub > 0) then
                ! velocity and acceleration in same direction
                del_tau = 2.0d0 * ABS(c)/(ABS(v_sub) + SQRT(v_sub**2 - 4.0d0*a*c))
                l_f = l_alongV
                if (del_tau <= 0.0d0) then
                    stop "Have issue with del_tau for v,a in same direction"
                end if
            else if (v_sub**2 - 4.0d0*a*c > 0.0d0) then
                ! v and a opposite direction, but particle can still reach boundary along v
                del_tau = 2.0d0 * ABS(c)/(ABS(v_sub) + SQRT(v_sub**2 - 4.0d0*a*c))
                l_f = l_alongV
                if (del_tau <= 0.0d0) then
                    print *, "In regular substep routine"
                    stop "Have issue with del_tau for v,a in opposite direction, but still reach boundary along v"
                end if
            end if
            if (del_tau >= del_t-timePassed) then
                ! For particles returning to boundary, make guess for maximum distance could go before v_f = 0 in half sub-step
                if (ABS((del_t - timePassed) *v_sub/4.0d0/world%dx_dl(l_cell)) < 1.0) then
                    a = (part%q / part%mass / 2.0d0) * self%getEFieldPeriodic(l_sub + (del_t - timePassed) *v_sub/8.0d0/world%dx_dl(l_cell), l_cell, world)
                    if (a*v_sub < 0) then
                        del_tau_tmp = ABS(v_sub/a)
                        if (del_tau_tmp < del_tau) then
                            del_tau = del_tau_tmp
                            l_f = l_sub
                        end if
                        if (del_tau <= 0) then
                            stop "Have issue with del_tau for v,a in opposite direction, boundary opposite v"
                        end if
                    end if
                end if
                ! try for minimum distance particle can go
                if (del_tau >= del_t - timePassed) then
                    a = (part%q / part%mass / 2.0d0) * self%getEFieldPeriodic(l_sub + SIGN(1.0d-10, v_sub), l_cell, world)
                    if (a*v_sub < 0) then
                        del_tau_tmp = ABS(v_sub/a)
                        if (del_tau_tmp < del_tau) then
                            del_tau = del_tau_tmp
                            l_f = l_sub
                        end if
                        if (del_tau <= 0) then
                            stop "Have issue with del_tau for v,a in opposite direction, boundary opposite v, second try close to boundary"
                        end if
                    end if
                end if
            end if
        else
            !Free particle drift
            del_tau = (l_alongV - l_sub) * world%dx_dl(l_cell)/v_sub
            v_f = (l_alongV - l_sub) * world%dx_dl(l_cell) / del_tau
            l_f = l_alongV
            if (del_tau <= 0.0d0) then
                stop "Have issue with del_tau for a = 0"
            end if
        end if

    end subroutine getDelTauSubstepPeriodic

    subroutine getDelTauInitialSubStepDirichlet(self, world, part, l_sub, l_f, v_sub, v_f, del_tau, l_awayV, l_alongV, l_cell, a, c)
        class(potentialSolver), intent(in) :: self
        type(Domain), intent(in) :: world
        type(Particle), intent(in) :: part
        real(real64), intent(in out) :: l_sub, l_f, v_sub, v_f, del_tau, l_awayV, l_alongV, a, c
        integer(int32), intent(in) :: l_cell
        a = (part%q / part%mass / 2.0d0) * self%getEFieldDirichlet(l_sub, l_cell, world)
        ! Particle first between nodes, so solve quadratic for that particle depending on conditions
        if ((v_sub/=0.0d0)) then ! make first case, since pretty much always likely to be the case (could just not have, assume always field, never have to check)
            c = (l_sub - l_alongV) * world%dx_dl(l_cell)
            if (a*v_sub > 0) then
                ! velocity and acceleration in same direction
                del_tau = 2.0d0 * ABS(c)/(ABS(v_sub) + SQRT(v_sub**2 - 4.0d0*a*c))
                l_f = l_alongV
                if (del_tau <= 0.0d0) then
                    stop "Have issue with del_tau for v,a in same direction"
                end if
            else if (v_sub**2 - 4.0d0*a*c > 0.0d0) then
                ! v and a opposite direction, but particle can still reach boundary along v
                del_tau = 2.0d0 * ABS(c)/(ABS(v_sub) + SQRT(v_sub**2 - 4.0d0*a*c))
                l_f = l_alongV
                if (del_tau <= 0.0d0) then
                    stop "Have issue with del_tau for v,a in opposite direction, but still reach boundary along v"
                end if
            else
                ! v and a opposite direction, boundary opposite direction of v
                c = (l_sub - l_awayV) * world%dx_dl(l_cell)
                del_tau = 2.0d0 * ABS(c)/(SQRT(v_sub**2 - 4.0d0*a*c) - ABS(v_sub))
                l_f = l_awayV
                if (del_tau <= 0) then
                    stop "Have issue with del_tau for v,a in opposite direction, boundary opposite v"
                end if
            end if
        else
            ! only in direction of field: USE l_alongV AS BOUNDARY ALONG DIRECTION OF a SINCE VELOCITY = 0!!!
            if (a > 0) then
                l_alongV = real(INT(l_sub) + 1, kind = real64)
            else
                l_alongV = real(INT(l_sub), kind = real64)
            end if
            c = (l_sub - l_alongV) * world%dx_dl(l_cell)
            del_tau = SQRT(-c/a)
            l_f = l_alongV
            if (del_tau <= 0.0d0) then
                stop "Have issue with del_tau for v = 0"
            end if
        end if
    end subroutine getDelTauInitialSubStepDirichlet

    subroutine getDelTauSubStepDirichlet(self, world, part, l_sub, l_f, v_sub, v_f, del_tau, l_alongV, l_cell, a, c)
        class(potentialSolver), intent(in) :: self
        type(Domain), intent(in) :: world
        type(Particle), intent(in) :: part
        real(real64), intent(in out) :: l_sub, l_f, v_sub, v_f, del_tau, l_alongV, a, c
        integer(int32), intent(in) :: l_cell
        a = (part%q / part%mass / 2.0d0) * self%getEFieldDirichlet(l_sub, l_cell, world)
        ! Particle first between nodes, so solve quadratic for that particle depending on conditions
        if ((a/=0.0d0)) then ! make first case, since pretty much always likely to be the case (could just not have, assume always field exists, never have to check)
            c = (l_sub - l_alongV) * world%dx_dl(l_cell)
            if (a*v_sub > 0.0d0) then
                ! velocity and acceleration in same direction
                del_tau = 2.0d0 * ABS(c)/(ABS(v_sub) + SQRT(v_sub**2 - 4.0d0*a*c))
                l_f = l_alongV
                if (del_tau <= 0.0d0) then
                    stop "Have issue with del_tau for v,a in same direction"
                end if
            else if (v_sub**2 - 4.0d0*a*c > 0.0d0) then
                ! v and a opposite direction, but particle can still reach boundary along v
                del_tau = 2.0d0 * ABS(c)/(ABS(v_sub) + SQRT(v_sub**2 - 4.0d0*a*c))
                l_f = l_alongV
                if (del_tau <= 0.0d0) then
                    stop "Have issue with del_tau for v,a in opposite direction, but still reach boundary along v"
                end if
            else
                ! v and a opposite direction, reverses back to initial position
                del_tau = ABS(v_sub)/ABS(a)
                l_f = l_sub
                if (del_tau <= 0.0d0) then
                    stop "Have issue with del_tau for v,a in opposite direction, boundary opposite v"
                end if
            end if
        else
            !Free particle drift
            del_tau = (l_alongV - l_sub) * world%dx_dl(l_cell)/v_sub
            v_f = (l_alongV - l_sub) * world%dx_dl(l_cell) / del_tau
            l_f = l_alongV
            if (del_tau <= 0.0d0) then
                stop "Have issue with del_tau for a = 0"
            end if
        end if
    end subroutine getDelTauSubStepDirichlet

    function getAlpha(l_sub, l_f, world) result(res)
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: l_sub, l_f
        real(real64) :: res, x_f, x_i
        if (l_sub /= l_f) then
            x_i = world%grid(NINT(l_sub)) + (l_sub - NINT(l_sub)) * world%dx_dl(NINT(l_sub))
            x_f = world%grid(NINT(l_f)) + (l_f - NINT(l_f)) * world%dx_dl(NINT(l_f))
            res = (l_f - l_sub)/(x_f - x_i)
        else
            res = 1.0d0/world%dx_dl(NINT(l_sub))
        end if

    end function getAlpha

    subroutine analyticalParticleMover(self, world, q, mass, l_sub, l_f, v_sub, v_f, l_cell, del_t, timePassed)
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: l_sub, v_sub, del_t, timePassed, q, mass
        real(real64), intent(in out) :: l_f, v_f
        integer(int32), intent(in) :: l_cell
        real(real64) :: del_tau, del_tau_sqr, real_l_cell, dx
        dx = world%dx_dl(l_cell)
        del_tau = del_t - timePassed
        del_tau_sqr = del_tau**2
        real_l_cell = real(l_cell, kind = real64)
        l_f = (-4.0d0*del_tau**2 *real_l_cell*self%phi(l_cell)*q + 2.0d0*del_tau**2 *real_l_cell *self%phi(l_cell+1)*q + & 
        2.0d0* del_tau_sqr *real_l_cell*self%phi(l_cell-1)*q - 4.0d0*del_tau**2 *real_l_cell*self%phi_f(l_cell)*q + &
        2.0d0*del_tau_sqr *real_l_cell*self%phi_f(l_cell+1)*q + 2.0d0*del_tau_sqr *real_l_cell*self%phi_f(l_cell-1)*q + &
        2.0d0*del_tau_sqr  *l_sub*self%phi(l_cell)*q - del_tau_sqr *l_sub*self%phi(l_cell+1)*q - &
        del_tau_sqr *l_sub*self%phi(l_cell-1)*q + 2.0d0*del_tau_sqr *l_sub*self%phi_f(l_cell)*q - &
        del_tau_sqr *l_sub*self%phi_f(l_cell+1)*q - del_tau_sqr *l_sub*self%phi_f(l_cell-1)*q - &
        del_tau_sqr *self%phi(l_cell+1)*q + del_tau_sqr *self%phi(l_cell-1)*q - del_tau_sqr *self%phi_f(l_cell+1)*q + &
        del_tau_sqr *self%phi_f(l_cell-1)*q + 8.0d0*(del_t - timePassed)*mass*v_sub*dx + &
        8.0d0*l_sub*mass*dx**2 )/(-2.0d0*del_tau_sqr *self%phi(l_cell)*q + del_tau_sqr *self%phi(l_cell+1)*q + &
        del_tau_sqr *self%phi(l_cell-1)*q - 2.0d0*del_tau_sqr *self%phi_f(l_cell)*q + del_tau_sqr *self%phi_f(l_cell+1)*q + &
        del_tau_sqr *self%phi_f(l_cell-1)*q + 8.0d0*mass*dx**2)
        v_f = 2.0d0 * (l_f - l_sub) * dx / del_tau - v_sub
        if (NINT(l_f) /= l_cell) then
            print *, "Have final l_f outside initial cell"
            print *, "particle charge is:", q
            print *, "l_sub is:", l_sub
            print *, "l_f is:", l_f
            print *, "v_sub is:", v_sub
            print *, "v_f is:", v_f
        end if

    end subroutine analyticalParticleMover

    subroutine analyticalParticleMoverPeriodic(self, world, q, mass, l_sub, l_f, v_sub, v_f, l_cell, del_t, timePassed)
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: l_sub, v_sub, del_t, timePassed, q, mass
        real(real64), intent(in out) :: l_f, v_f
        integer(int32), intent(in) :: l_cell
        real(real64) :: del_tau, del_tau_sqr, real_l_cell, dx
        dx = world%dx_dl(l_cell)
        del_tau = del_t - timePassed
        del_tau_sqr = del_tau**2
        real_l_cell = real(l_cell, kind = real64)
        l_f = (-4.0d0*del_tau**2 *real_l_cell*self%phi(1)*q + 2.0d0*del_tau**2 *real_l_cell *self%phi(2)*q + & 
        2.0d0* del_tau_sqr *real_l_cell*self%phi(NumberXNodes-1)*q - 4.0d0*del_tau**2 *real_l_cell*self%phi_f(1)*q + &
        2.0d0*del_tau_sqr *real_l_cell*self%phi_f(2)*q + 2.0d0*del_tau_sqr *real_l_cell*self%phi_f(NumberXNodes-1)*q + &
        2.0d0*del_tau_sqr  *l_sub*self%phi(1)*q - del_tau_sqr *l_sub*self%phi(2)*q - &
        del_tau_sqr *l_sub*self%phi(NumberXNodes-1)*q + 2.0d0*del_tau_sqr *l_sub*self%phi_f(1)*q - &
        del_tau_sqr *l_sub*self%phi_f(2)*q - del_tau_sqr *l_sub*self%phi_f(NumberXNodes-1)*q - &
        del_tau_sqr *self%phi(2)*q + del_tau_sqr *self%phi(NumberXNodes-1)*q - del_tau_sqr *self%phi_f(2)*q + &
        del_tau_sqr *self%phi_f(NumberXNodes-1)*q + 8.0d0*(del_t - timePassed)*mass*v_sub*dx + &
        8.0d0*l_sub*mass*dx**2 )/(-2.0d0*del_tau_sqr *self%phi(1)*q + del_tau_sqr *self%phi(2)*q + &
        del_tau_sqr *self%phi(NumberXNodes-1)*q - 2.0d0*del_tau_sqr *self%phi_f(1)*q + del_tau_sqr *self%phi_f(2)*q + &
        del_tau_sqr *self%phi_f(NumberXNodes-1)*q + 8.0d0*mass*dx**2)
        v_f = 2.0d0 * (l_f - l_sub) * dx / del_tau - v_sub
        if (NINT(l_f) /= l_cell) then
            print *, "Have final l_f outside initial cell"
            print *, "particle charge is:", q
            stop 
        end if
    end subroutine analyticalParticleMoverPeriodic

    subroutine picardIterParticles(self, world, q, mass, l_sub, l_f, v_sub, v_f, l_cell, del_t, timePassed)
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        real(real64), intent(in out) :: l_sub, l_f, v_sub, v_f
        real(real64), intent(in) :: del_t, timePassed, q, mass
        integer(int32), intent(in) :: l_cell
        integer(int32) :: i
        real(real64) :: l_f_previous
        l_f_previous = l_f
        v_f = v_sub + (q/mass) * self%getEField((l_sub + l_f)/2.0d0, l_cell, world) * (del_t - timePassed)
        l_f = (v_f + v_sub) * (del_t - timePassed) /world%dx_dl(l_cell) / 2.0d0 + l_sub
        do i = 1, 50
            l_f_previous = l_f
            v_f = v_sub + (q/mass) * self%getEField((l_sub + l_f)/2.0d0, l_cell, world) * (del_t - timePassed)
            l_f = (v_f + v_sub) * (del_t - timePassed) /world%dx_dl(l_cell) / 2.0d0 + l_sub
            if (NINT(l_f) /= l_cell) then
                l_f = l_cell + SIGN(0.5d0, l_f - l_cell)
            end if
            if (ABS(l_f - l_f_previous) < eps_r) exit
        end do
        if (NINT(l_f) /= l_cell) then
            print *, "Have final l_f outside initial cell"
            print *, "particle charge is:", q
            print *, "a direction is:", (q/mass/2.0d0) * self%getEField((l_sub + l_f)/2.0d0, l_cell, world)
            print *, "del_t remaining is:", del_t - timePassed
            print *, "l_sub is:", l_sub
            print *, "l_f is:", l_f
            print *, "v_sub is:", v_sub
            print *, "v_f is:", v_f
            v_f = v_sub + (q/mass) * self%getEField((l_sub + l_f)/2.0d0, l_cell, world) * (del_t - timePassed)
            l_f = (v_f + v_sub) * (del_t - timePassed) / world%dx_dl(l_cell) / 2.0d0 + l_sub
            print *, "next l_f is:", l_f
            print *, "next v_f is:", v_f
            v_f = v_sub + (q/mass) * self%getEField((l_sub + l_f)/2.0d0, l_cell, world) * (del_t - timePassed)
            l_f = (v_f + v_sub) * (del_t - timePassed) / world%dx_dl(l_cell) / 2.0d0 + l_sub
            print *, "next l_f is:", l_f
            print *, "next v_f is:", v_f
            call self%analyticalParticleMover(world, q, mass, l_sub, l_f, v_sub, v_f, l_cell, del_t, timePassed)
            print *, "l_f analytical is:", l_f
            stop 
        end if
        if (i == 51) then
            stop "Picard not converged"
        end if
    end subroutine picardIterParticles

    subroutine depositJSubStep(self, world, q, w_p, l_sub, l_f, v_f, v_sub, l_cell, del_tau, del_t)
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: q, w_p, v_f, v_sub, del_tau, del_t, l_sub, l_f
        integer(int32), intent(in) :: l_cell
        real(real64) :: d
        d = (l_sub + l_f)/2.0d0 - REAL(l_cell, kind = real64) + 0.5d0
        self%J(l_cell-1) = self%J(l_cell-1) + (1.0d0 - d) * w_p * q * (v_f + v_sub)*del_tau/2.0d0/world%dx_dl(l_cell)/del_t
        self%J(l_cell) = self%J(l_cell) + d * w_p * q * (v_f + v_sub)*del_tau/2.0d0/world%dx_dl(l_cell)/del_t
    end subroutine depositJSubStep

    subroutine depositJSubStepPeriodic(self, world, q, w_p, l_sub, l_f, v_f, v_sub, l_cell, del_tau, del_t)
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: q, w_p, v_f, v_sub, del_tau, del_t, l_sub, l_f
        integer(int32), intent(in) :: l_cell
        real(real64) :: d
        d = (l_sub + l_f)/2.0d0 - REAL(l_cell, kind = real64) + 0.5d0
        self%J(NumberXNodes-1) = self%J(NumberXNodes-1) + (1.0d0 - d) * w_p * q * (v_f + v_sub)*del_tau/2.0d0/world%dx_dl(l_cell)/del_t
        self%J(1) = self%J(1) + d * w_p * q * (v_f + v_sub)*del_tau/2.0d0/world%dx_dl(l_cell)/del_t

    end subroutine depositJSubStepPeriodic

    subroutine depositJSubStepDirichlet(self, world, q, w_p, l_sub, l_cell, v_f, v_sub, del_tau, del_t)
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: q, w_p, v_f, v_sub, del_tau, del_t, l_sub
        integer(int32), intent(in) :: l_cell
        self%J(INT(l_sub)) = self%J(INT(l_sub)) + w_p * q * (v_f + v_sub)*del_tau/2.0d0/world%dx_dl(l_cell)/del_t
    end subroutine depositJSubStepDirichlet


    subroutine depositJ(self, particleList, world, del_t)
        ! particle substepping procedure which deposits J, also used for general moving of particles
        ! boolDepositJ = true : deposit J
        ! boolDepositJ = false: only move particles (done after non-linear iteration), also delete particles for wall collisions (saves extra loop in collisions)
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        type(Particle), intent(in out) :: particleList(:)
        real(real64), intent(in) :: del_t
        !a and c correspond to quadratic equations | l_alongV is nearest integer boundary along velocity component, away is opposite
        real(real64) :: l_f, l_sub, v_sub, v_f, timePassed, del_tau, l_alongV, l_awayV, a, c
        integer(int32) :: subStepNum, j, i, delIdx, l_cell
        self%J = 0.0d0
        loopSpecies: do j = 1, numberChargedParticles
            delIdx = 0
            loopParticles: do i = 1, particleList(j)%N_p
                v_sub = particleList(j)%phaseSpace(2,i)
                l_sub = particleList(j)%phaseSpace(1,i)
                v_f = v_sub
                l_f = l_sub
                timePassed = 0.0d0
                subStepNum = 0

                ! First substep
                l_cell = NINT(l_sub)
                firstStep: SELECT CASE (world%boundaryConditions(l_cell))
                CASE(0)
                    ! Within Domain
                    l_alongV = real(l_cell, kind = real64) + SIGN(0.5d0, v_sub)
                    l_awayV = real(l_cell, kind = real64) - SIGN(0.5d0, v_sub)
                    call self%getDelTauInitialSubStep(world, particleList(j), l_sub, l_f, v_sub, v_f, del_tau, del_t, l_alongV, l_awayV, l_cell, a, c)
                    if (del_tau >= del_t) then
                        call self%analyticalParticleMover(world, particleList(j)%q, particleList(j)%mass, l_sub, l_f, v_sub, v_f, l_cell, del_t, timePassed)
                        !call self%picardIterParticles(world, particleList(j)%q, particleList(j)%mass, l_sub, l_f, v_sub, v_f, l_cell, del_t, timePassed)
                        call self%depositJSubStep(world, particleList(j)%q, particleList(j)%w_p, l_sub, l_f, v_f, v_sub, l_cell, del_t, del_t) 
                        if (NINT(l_f) /= NINT(l_sub)) then
                            print *, "After initial domain substep, l_f is not in correct cell"
                            stop
                        end if
                        timePassed = del_t  
                    else
                        v_f = 2.0d0 * (l_f - l_sub) * world%dx_dl(l_cell) / del_tau - v_sub
                        call self%depositJSubStep(world, particleList(j)%q, particleList(j)%w_p, l_sub, l_f, v_f, v_sub, l_cell, del_tau, del_t)
                        timePassed = timePassed + del_tau
                        if ((l_f /= l_alongV) .and. (l_f /= l_awayV)) then
                            print *, "l_sub is:", l_sub
                            print *, "l_f is:", l_f
                            stop "l_f is not correct boundary after initial domain substep"
                        end if
                        ! now final position/velocity becomes next starting position/velocity
                        l_sub = l_f
                        v_sub = v_f
                    end if
                CASE(1)
                    !Dirichlet
                    l_alongV = NINT(2.0d0*l_sub + SIGN(0.5d0, v_sub))/2.0d0
                    l_awayV = NINT(2.0d0*l_sub - SIGN(0.5d0, v_sub))/2.0d0
                    call self%getDelTauInitialSubStepDirichlet(world, particleList(j), l_sub, l_f, v_sub, v_f, del_tau, l_awayV, l_alongV, l_cell, a, c)
                    if (del_tau >= del_t-timePassed) then
                        ! Add directly to J with no substep
                        l_f = v_sub * del_t / world%dx_dl(l_cell) + (a/ world%dx_dl(l_cell)) * del_t**2 + l_sub
                        v_f = 2.0d0 * (l_f - l_sub) * world%dx_dl(l_cell) / del_t - v_sub
                        call self%depositJSubStepDirichlet(world, particleList(j)%q, particleList(j)%w_p, l_sub, l_cell, v_f, v_sub, del_t, del_t)
                        if (NINT(l_f) /= NINT(l_sub)) then
                            print *, "After initial Dirichlet substep, l_f is not in correct cell"
                            stop
                        end if
                        timePassed = del_t  
                    else
                        v_f = 2.0d0 * (l_f - l_sub) * world%dx_dl(l_cell) / del_tau - v_sub
                        call self%depositJSubStepDirichlet(world, particleList(j)%q, particleList(j)%w_p, l_sub, l_cell, v_f, v_sub, del_tau, del_t)
                        timePassed = del_tau
                        if (l_f == l_cell) then ! if particle is now on node, must be boundary, exit
                            timePassed = del_t
                        end if
                        if ((l_f /= l_alongV) .and. (l_f /= l_awayV)) then
                            print *, "l_sub is:", l_sub
                            print *, "l_f is:", l_f
                            stop "l_f is not correct boundary after initial dirichlet domain substep"
                        end if
                        ! now final position/velocity becomes next starting position/velocity
                        l_sub = l_f
                        v_sub = v_f
                    end if
                CASE(3)
                    l_alongV = real(l_cell, kind = real64) + SIGN(0.5d0, v_sub)
                    l_awayV = real(l_cell, kind = real64) - SIGN(0.5d0, v_sub)
                    call self%getDelTauInitialSubStepPeriodic(world, particleList(j), l_sub, l_f, v_sub, v_f, del_tau, del_t, l_alongV, l_awayV, l_cell, a, c)
                    if (del_tau >= del_t-timePassed) then
                        call self%analyticalParticleMoverPeriodic(world, particleList(j)%q, particleList(j)%mass, l_sub, l_f, v_sub, v_f, l_cell, del_t, timePassed)    
                        call self%depositJSubStepPeriodic(world, particleList(j)%q, particleList(j)%w_p, l_sub, l_f, v_f, v_sub, l_cell, del_t, del_t) 
                        if (NINT(l_f) /= NINT(l_sub)) then
                            print *, "After initial periodic substep, l_f is not in correct cell"
                            stop
                        end if
                        if (l_f < 1) then
                            l_f = NumberXNodes + (l_f - l_cell)
                        else if (l_f > NumberXNodes) then
                            l_f = l_f - l_cell + 1
                        end if
                        timePassed = del_t  
                    else
                        v_f = 2.0d0 * (l_f - l_sub) * world%dx_dl(l_cell) / del_tau - v_sub
                        call self%depositJSubStepPeriodic(world, particleList(j)%q, particleList(j)%w_p, l_sub, l_f, v_f, v_sub, l_cell, del_tau, del_t)
                        if ((l_f /= l_alongV) .and. (l_f /= l_awayV)) then
                            print *, "l_sub is:", l_sub
                            print *, "l_f is:", l_f
                            stop "l_f is not correct boundary after initial periodic substep"
                        end if
                        if (l_f < 1) then
                            l_f = NumberXNodes + (l_f - l_cell)
                        else if (l_f > NumberXNodes) then
                            l_f = l_f - l_cell + 1
                        end if
                        timePassed = del_tau
                        ! now final position/velocity becomes next starting position/velocity
                        l_sub = l_f
                        v_sub = v_f
                    end if
                CASE default
                    print *, "The boundary condition case doesn't exist!"
                    print *, "Happened within initial substep"
                    stop
                END SELECT firstStep
                do while((timePassed < del_t))
                    
                    l_cell = NINT(l_sub + SIGN(0.1d0, v_sub))
                    subStep: SELECT CASE (world%boundaryConditions(l_cell))
                    CASE(0)
                        l_alongV = l_sub + SIGN(1.0d0, v_sub)
                        call self%getDelTauSubStep(world, particleList(j), l_sub, l_f, v_sub, v_f, del_tau, timePassed, del_t, l_alongV, l_cell, a, c)
                        if (del_tau >= del_t-timePassed) then
                            ! Add directly to J with no substep
                            call self%analyticalParticleMover(world, particleList(j)%q, particleList(j)%mass, l_sub, l_f, v_sub, v_f, l_cell, del_t, timePassed)
                            !call self%picardIterParticles(world, particleList(j)%q, particleList(j)%mass, l_sub, l_f, v_sub, v_f, l_cell, del_t, timePassed)
                            call self%depositJSubStep(world, particleList(j)%q, particleList(j)%w_p, l_sub, l_f, v_f, v_sub, l_cell, del_t - timePassed, del_t)
                            if (NINT(l_f) /= l_cell) then
                                print *, "After ongoing substep, l_f is not in correct cell"
                            end if
                            timePassed = del_t  
                        else
                            v_f = 2.0d0 * (l_f - l_sub) * world%dx_dl(l_cell) / del_tau - v_sub
                            call self%depositJSubStep(world, particleList(j)%q, particleList(j)%w_p, l_sub, l_f, v_f, v_sub, l_cell, del_tau, del_t) 
                            timePassed = timePassed + del_tau
                            if ((MOD(l_f, 0.5d0) /= 0.0d0) .or. (ABS(l_f - l_sub) > 1.0d0)) then
                                print *, l_f
                                stop "l_f is not half integer after subStep or it is too far away"
                            end if
                            ! now final position/velocity becomes next starting position/velocity
                            l_sub = l_f
                            v_sub = v_f
                        end if
                    CASE(1)
                        !Near Dirichlet node
                        l_alongV = l_sub + SIGN(0.5d0, v_sub)
                        call self%getDelTauSubStepDirichlet(world, particleList(j), l_sub, l_f, v_sub, v_f, del_tau, l_alongV, l_cell, a, c)
                        if (del_tau >= del_t-timePassed) then
                            ! Add directly to J with no substep
                            l_f = v_sub * (del_t-timePassed) / world%dx_dl(l_cell) + (a/ world%dx_dl(l_cell)) * (del_t - timePassed)**2 + l_sub
                            v_f = 2.0d0 * (l_f - l_sub) * world%dx_dl(l_cell) / (del_t - timePassed) - v_sub
                            call self%depositJSubStepDirichlet(world, particleList(j)%q, particleList(j)%w_p, l_sub, l_cell, v_f, v_sub, del_t-timePassed, del_t)
                            if (NINT(l_f) /= l_cell) then
                                print *, "a is:", a
                                print *, world%dx_dl
                                print *, 'v_sub is:', v_sub
                                print *, "remaining time:", del_t - timePassed
                                print *, "l_f is:", l_f
                                print *, "l_sub is:", l_sub
                                print *, "l_cell is:", l_cell
                                print *, "After ongoing Dirichlet last substep, l_f is not in correct cell"
                                stop
                            end if
                            timePassed = del_t  
                        else
                            v_f = 2.0d0 * (l_f - l_sub) * world%dx_dl(l_cell) / del_tau - v_sub
                            call self%depositJSubStepDirichlet(world, particleList(j)%q, particleList(j)%w_p, l_sub, l_cell, v_f, v_sub, del_tau, del_t) 
                            if (l_f == l_cell) exit ! if particle is now on node, must be boundary, exit
                            if ((MOD(l_f, 0.5d0) /= 0.0d0) .or. (ABS(l_f - l_sub) > 0.5d0)) then
                                print *, "l_sub is:", l_sub
                                print *, "l_f is:", l_f
                                print *, "a is:", a 
                                print *, "v_sub is:", v_sub
                                print *, "v_f is:", v_f
                                stop "l_f is not correct boundary after ongoing dirichlet substep"
                            end if
                            ! now final position/velocity becomes next starting position/velocity
                            l_sub = l_f
                            v_sub = v_f
                        end if
                    CASE(3)
                        l_alongV = l_sub + SIGN(1.0d0, v_sub)
                        call self%getDelTauSubStepPeriodic(world, particleList(j), l_sub, l_f, v_sub, v_f, del_tau, timePassed, del_t, l_alongV, l_cell, a, c)
                        if (del_tau >= del_t-timePassed) then
                            ! Add directly to J with no substep
                            call self%analyticalParticleMoverPeriodic(world, particleList(j)%q, particleList(j)%mass, l_sub, l_f, v_sub, v_f, l_cell, del_t, timePassed)
                            call self%depositJSubStepPeriodic(world, particleList(j)%q, particleList(j)%w_p, l_sub, l_f, v_f, v_sub, l_cell, del_t - timePassed, del_t)
                            if (NINT(l_f) /= l_cell) then
                                print *, "After final periodic substep, l_f is not in correct cell"
                                stop
                            end if
                            if (l_f < 1) then
                                l_f = NumberXNodes + (l_f - l_cell)
                            else if (l_f > NumberXNodes) then
                                l_f = l_f - l_cell + 1
                            end if
                            timePassed = del_t  
                        else
                            v_f = 2.0d0 * (l_f - l_sub) * world%dx_dl(l_cell) / del_tau - v_sub
                            call self%depositJSubStepPeriodic(world, particleList(j)%q, particleList(j)%w_p, l_sub, l_f, v_f, v_sub, l_cell, del_tau, del_t) 
                            if ((MOD(l_f, 0.5d0) /= 0.0d0) .or. (ABS(l_f - l_sub) > 1.0d0)) then
                                print *, l_f
                                stop "l_f is not half integer after ogoing periodic subStep or it is too far away"
                            end if
                            if (l_f < 1) then
                                l_f = NumberXNodes + (l_f - l_cell)
                            else if (l_f > NumberXNodes) then
                                l_f = l_f - l_cell + 1
                            end if
                            timePassed = timePassed + del_tau
                            ! now final position/velocity becomes next starting position/velocity
                            l_sub = l_f
                            v_sub = v_f
                        end if
                    CASE default
                        print *, "The boundary condition case doesn't exist!"
                        print *, "Happened within ongoing substep"
                        stop
                    END SELECT subStep
                   
                    subStepNum = subStepNum + 1
                end do
                if ((l_f < 1) .or. (l_f > NumberXNodes)) then
                    stop "Have particles travelling outside domain!"
                end if
            end do loopParticles
            
        end do loopSpecies
        
    end subroutine depositJ


    ! -------------------------------------------- Particle mover without boolean checks for depositing J -----------------------------------------------------------

    subroutine moveParticles(self, particleList, world, del_t)
        ! particle mover to avoid the boolean checks which mostly don't happen when depositing J
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        type(Particle), intent(in out) :: particleList(:)
        real(real64), intent(in) :: del_t
        !a and c correspond to quadratic equations | l_alongV is nearest integer boundary along velocity component, away is opposite
        real(real64) :: l_f, l_sub, v_sub, v_f, timePassed, del_tau, l_alongV, l_awayV, a, c
        integer(int32) :: subStepNum, j, i, delIdx, l_cell
        logical :: delParticle
        delParticle = .false.
        loopSpecies: do j = 1, numberChargedParticles
            delIdx = 0
            loopParticles: do i = 1, particleList(j)%N_p
                v_sub = particleList(j)%phaseSpace(2,i)
                l_sub = particleList(j)%phaseSpace(1,i)
                v_f = v_sub
                l_f = l_sub
                timePassed = 0.0d0
                subStepNum = 0
                ! First substep
                l_cell = NINT(l_sub)
                firstStep: SELECT CASE (world%boundaryConditions(l_cell))
                CASE(0)
                    ! Within Domain
                    l_alongV = real(l_cell, kind = real64) + SIGN(0.5d0, v_sub)
                    l_awayV = real(l_cell, kind = real64) - SIGN(0.5d0, v_sub)
                    call self%getDelTauInitialSubStep(world, particleList(j), l_sub, l_f, v_sub, v_f, del_tau, del_t, l_alongV, l_awayV, l_cell, a, c)
                    if (del_tau >= del_t) then
                        call self%analyticalParticleMover(world, particleList(j)%q, particleList(j)%mass, l_sub, l_f, v_sub, v_f, l_cell, del_t, timePassed)
                        !call self%picardIterParticles(world, particleList(j)%q, particleList(j)%mass, l_sub, l_f, v_sub, v_f, l_cell, del_t, timePassed)
                        if (NINT(l_f) /= NINT(l_sub)) then
                            print *, "After initial domain substep, l_f is not in correct cell"
                            stop
                        end if
                        timePassed = del_t 
                        particleList(j)%phaseSpace(1, i-delIdx) = l_f
                        particleList(j)%phaseSpace(2,i-delIdx) = v_f
                        particleList(j)%phaseSpace(3:4, i-delIdx) = particleList(j)%phaseSpace(3:4, i) 
                    else
                        v_f = 2.0d0 * (l_f - l_sub) * world%dx_dl(l_cell) / del_tau - v_sub
                        timePassed = timePassed + del_tau
                        if ((l_f /= l_alongV) .and. (l_f /= l_awayV)) then
                            print *, "l_sub is:", l_sub
                            print *, "l_f is:", l_f
                            stop "l_f is not correct boundary after initial domain substep"
                        end if
                        ! now final position/velocity becomes next starting position/velocity
                        l_sub = l_f
                        v_sub = v_f
                    end if
                CASE(1)
                    !Dirichlet
                    l_alongV = NINT(2.0d0*l_sub + SIGN(0.5d0, v_sub))/2.0d0
                    l_awayV = NINT(2.0d0*l_sub - SIGN(0.5d0, v_sub))/2.0d0
                    call self%getDelTauInitialSubStepDirichlet(world, particleList(j), l_sub, l_f, v_sub, v_f, del_tau, l_awayV, l_alongV, l_cell, a, c)
                    if (del_tau >= del_t-timePassed) then
                        ! Add directly to J with no substep
                        l_f = v_sub * del_t / world%dx_dl(l_cell) + (a/ world%dx_dl(l_cell)) * del_t**2 + l_sub
                        v_f = 2.0d0 * (l_f - l_sub) * world%dx_dl(l_cell) / del_t - v_sub
                        if (NINT(l_f) /= NINT(l_sub)) then
                            print *, "After initial Dirichlet substep, l_f is not in correct cell"
                            stop
                        end if
                        timePassed = del_t 
                        particleList(j)%phaseSpace(1, i-delIdx) = l_f
                        particleList(j)%phaseSpace(2,i-delIdx) = v_f
                        particleList(j)%phaseSpace(3:4, i-delIdx) = particleList(j)%phaseSpace(3:4, i) 
                    else
                        v_f = 2.0d0 * (l_f - l_sub) * world%dx_dl(l_cell) / del_tau - v_sub
                        if (l_f == l_cell) then ! if particle is now on node, must be boundary, exit
                            delIdx = delIdx + 1
                            timePassed = del_t
                            self%particleEnergyLoss = self%particleEnergyLoss + particleList(j)%w_p * (v_f**2 + SUM(particleList(j)%phaseSpace(3:4, i)**2)) * particleList(j)%mass * 0.5d0 !J/m^2 in 1D
                            if (l_f == 1) then
                                self%particleChargeLoss(1, j) = self%particleChargeLoss(1, j) + particleList(j)%q * particleList(j)%w_p !C/m^2 in 1D
                            else
                                self%particleChargeLoss(2, j) = self%particleChargeLoss(2, j) + particleList(j)%q * particleList(j)%w_p !C/m^2 in 1D
                            end if
                        else
                            timePassed = del_tau
                        end if
                        if ((l_f /= l_alongV) .and. (l_f /= l_awayV)) then
                            print *, "l_sub is:", l_sub
                            print *, "l_f is:", l_f
                            stop "l_f is not correct boundary after initial dirichlet domain substep"
                        end if
                        ! now final position/velocity becomes next starting position/velocity
                        l_sub = l_f
                        v_sub = v_f
                    end if
                CASE(3)
                    l_alongV = real(l_cell, kind = real64) + SIGN(0.5d0, v_sub)
                    l_awayV = real(l_cell, kind = real64) - SIGN(0.5d0, v_sub)
                    call self%getDelTauInitialSubStepPeriodic(world, particleList(j), l_sub, l_f, v_sub, v_f, del_tau, del_t, l_alongV, l_awayV, l_cell, a, c)
                    if (del_tau >= del_t-timePassed) then
                        call self%analyticalParticleMoverPeriodic(world, particleList(j)%q, particleList(j)%mass, l_sub, l_f, v_sub, v_f, l_cell, del_t, timePassed)    
                        if (NINT(l_f) /= NINT(l_sub)) then
                            print *, "After initial periodic substep, l_f is not in correct cell"
                            stop
                        end if
                        if (l_f < 1) then
                            l_f = NumberXNodes + (l_f - l_cell)
                        else if (l_f > NumberXNodes) then
                            l_f = l_f - l_cell + 1
                        end if
                        timePassed = del_t
                        particleList(j)%phaseSpace(1, i-delIdx) = l_f
                        particleList(j)%phaseSpace(2,i-delIdx) = v_f
                        particleList(j)%phaseSpace(3:4, i-delIdx) = particleList(j)%phaseSpace(3:4, i)  
                    else
                        v_f = 2.0d0 * (l_f - l_sub) * world%dx_dl(l_cell) / del_tau - v_sub
                        if ((l_f /= l_alongV) .and. (l_f /= l_awayV)) then
                            print *, "l_sub is:", l_sub
                            print *, "l_f is:", l_f
                            stop "l_f is not correct boundary after initial periodic substep"
                        end if
                        if (l_f < 1) then
                            l_f = NumberXNodes + (l_f - l_cell)
                        else if (l_f > NumberXNodes) then
                            l_f = l_f - l_cell + 1
                        end if
                        timePassed = del_tau
                        ! now final position/velocity becomes next starting position/velocity
                        l_sub = l_f
                        v_sub = v_f
                    end if
                CASE default
                    print *, "The boundary condition case doesn't exist!"
                    print *, "Happened within initial substep"
                    stop
                END SELECT firstStep
                do while((timePassed < del_t))
                    
                    l_cell = NINT(l_sub + SIGN(0.1d0, v_sub))
                    subStep: SELECT CASE (world%boundaryConditions(l_cell))
                    CASE(0)
                        l_alongV = l_sub + SIGN(1.0d0, v_sub)
                        call self%getDelTauSubStep(world, particleList(j), l_sub, l_f, v_sub, v_f, del_tau, timePassed, del_t, l_alongV, l_cell, a, c)
                        if (del_tau >= del_t-timePassed) then
                            ! Add directly to J with no substep
                            call self%analyticalParticleMover(world, particleList(j)%q, particleList(j)%mass, l_sub, l_f, v_sub, v_f, l_cell, del_t, timePassed)
                            !call self%picardIterParticles(world, particleList(j)%q, particleList(j)%mass, l_sub, l_f, v_sub, v_f, l_cell, del_t, timePassed)
                            if (NINT(l_f) /= l_cell) then
                                print *, "After ongoing substep, l_f is not in correct cell"
                            end if
                            timePassed = del_t 
                            particleList(j)%phaseSpace(1, i-delIdx) = l_f
                            particleList(j)%phaseSpace(2,i-delIdx) = v_f
                            particleList(j)%phaseSpace(3:4, i-delIdx) = particleList(j)%phaseSpace(3:4, i) 
                        else
                            v_f = 2.0d0 * (l_f - l_sub) * world%dx_dl(l_cell) / del_tau - v_sub
                            timePassed = timePassed + del_tau
                            if ((MOD(l_f, 0.5d0) /= 0.0d0) .or. (ABS(l_f - l_sub) > 1.0d0)) then
                                print *, l_f
                                stop "l_f is not half integer after subStep or it is too far away"
                            end if
                            ! now final position/velocity becomes next starting position/velocity
                            l_sub = l_f
                            v_sub = v_f
                        end if
                    CASE(1)
                        !Near Dirichlet node
                        l_alongV = l_sub + SIGN(0.5d0, v_sub)
                        call self%getDelTauSubStepDirichlet(world, particleList(j), l_sub, l_f, v_sub, v_f, del_tau, l_alongV, l_cell, a, c)
                        if (del_tau >= del_t-timePassed) then
                            ! Add directly to J with no substep
                            l_f = v_sub * (del_t-timePassed) / world%dx_dl(l_cell) + (a/ world%dx_dl(l_cell)) * (del_t - timePassed)**2 + l_sub
                            v_f = 2.0d0 * (l_f - l_sub) * world%dx_dl(l_cell) / (del_t - timePassed) - v_sub
                            if (NINT(l_f) /= l_cell) then
                                print *, "After ongoing Dirichlet last substep, l_f is not in correct cell"
                                stop
                            end if
                            timePassed = del_t 
                            particleList(j)%phaseSpace(1, i-delIdx) = l_f
                            particleList(j)%phaseSpace(2,i-delIdx) = v_f
                            particleList(j)%phaseSpace(3:4, i-delIdx) = particleList(j)%phaseSpace(3:4, i) 
                        else
                            v_f = 2.0d0 * (l_f - l_sub) * world%dx_dl(l_cell) / del_tau - v_sub
                            if (l_f == l_cell) then 
                                delIdx = delIdx + 1
                                self%particleEnergyLoss = self%particleEnergyLoss + particleList(j)%w_p * (v_f**2 + SUM(particleList(j)%phaseSpace(3:4, i)**2)) * particleList(j)%mass * 0.5d0 !J/m^2 in 1D
                                if (l_f == 1) then
                                    self%particleChargeLoss(1, j) = self%particleChargeLoss(1, j) + particleList(j)%q * particleList(j)%w_p !C/m^2 in 1D
                                else
                                    self%particleChargeLoss(2, j) = self%particleChargeLoss(2, j) + particleList(j)%q * particleList(j)%w_p !C/m^2 in 1D
                                end if
                                exit ! if particle is now on node, must be boundary, exit
                            end if
                            if ((MOD(l_f, 0.5d0) /= 0.0d0) .or. (ABS(l_f - l_sub) > 0.5d0)) then
                                print *, "l_sub is:", l_sub
                                print *, "l_f is:", l_f
                                print *, "a is:", a 
                                print *, "v_sub is:", v_sub
                                print *, "v_f is:", v_f
                                stop "l_f is not correct boundary after ongoing dirichlet substep"
                            end if
                            ! now final position/velocity becomes next starting position/velocity
                            l_sub = l_f
                            v_sub = v_f
                        end if
                    CASE(3)
                        l_alongV = l_sub + SIGN(1.0d0, v_sub)
                        call self%getDelTauSubStepPeriodic(world, particleList(j), l_sub, l_f, v_sub, v_f, del_tau, timePassed, del_t, l_alongV, l_cell, a, c)
                        if (del_tau >= del_t-timePassed) then
                            ! Add directly to J with no substep
                            call self%analyticalParticleMoverPeriodic(world, particleList(j)%q, particleList(j)%mass, l_sub, l_f, v_sub, v_f, l_cell, del_t, timePassed)
                            if (NINT(l_f) /= l_cell) then
                                print *, "After final periodic substep, l_f is not in correct cell"
                                stop
                            end if
                            if (l_f < 1) then
                                l_f = NumberXNodes + (l_f - l_cell)
                            else if (l_f > NumberXNodes) then
                                l_f = l_f - l_cell + 1
                            end if
                            timePassed = del_t
                            particleList(j)%phaseSpace(1, i-delIdx) = l_f
                            particleList(j)%phaseSpace(2,i-delIdx) = v_f
                            particleList(j)%phaseSpace(3:4, i-delIdx) = particleList(j)%phaseSpace(3:4, i)  
                        else
                            v_f = 2.0d0 * (l_f - l_sub) * world%dx_dl(l_cell) / del_tau - v_sub
                            if ((MOD(l_f, 0.5d0) /= 0.0d0) .or. (ABS(l_f - l_sub) > 1.0d0)) then
                                print *, l_f
                                stop "l_f is not half integer after ogoing periodic subStep or it is too far away"
                            end if
                            if (l_f < 1) then
                                l_f = NumberXNodes + (l_f - l_cell)
                            else if (l_f > NumberXNodes) then
                                l_f = l_f - l_cell + 1
                            end if
                            timePassed = timePassed + del_tau
                            ! now final position/velocity becomes next starting position/velocity
                            l_sub = l_f
                            v_sub = v_f
                        end if
                    CASE default
                        print *, "The boundary condition case doesn't exist!"
                        print *, "Happened within ongoing substep"
                        stop
                    END SELECT subStep
                   
                    subStepNum = subStepNum + 1
                end do
                
                if ((l_f < 1) .or. (l_f > NumberXNodes)) then
                    stop "Have particles travelling outside domain!"
                end if
                ! When not depositing, then updating particles, overwrite deleted indices

            end do loopParticles
            particleList(j)%N_p = particleList(j)%N_p - delIdx
        end do loopSpecies
    end subroutine moveParticles

    ! ---------------- Initial Poisson Solver -------------------------------------------------

    subroutine solveInitialPotential(self, particleList, world)
        ! Solve for initial potential
        class(potentialSolver), intent(in out) :: self
        type(Particle), intent(in) :: particleList(:)
        type(Domain), intent(in) :: world
        call self%depositRho(particleList, world)
        call self%solve_tridiag_Poisson()
        ! Assume only use potential solver once, then need to generate matrix for Div-Ampere
        call self%construct_diagMatrix_Ampere(world)

    end subroutine solveInitialPotential




end module mod_potentialSolver