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
        real(real64) :: energyError, chargeError, particleEnergyLoss, Beta_k
        integer(int32) :: iterNumPicard, iterNumParticle, iterNumAdaptiveSteps, m_Anderson, amountTimeSplits
        real(real64) :: coeff_left, coeff_right ! these are coefficients (from world dimensions) needed with phi_left and phi_right in rhs of matrix equation
        real(real64), allocatable :: a_tri(:), b_tri(:), c_tri(:) !for thomas algorithm potential solver, a_tri is lower diagonal, b_tri middle, c_tri upper


    contains
        procedure, public, pass(self) :: depositRho
        procedure, public, pass(self) :: solve_tridiag_Poisson
        procedure, public, pass(self) :: solve_tridiag_Ampere
        procedure, public, pass(self) :: getEField
        procedure, public, pass(self) :: getTotalPE
        procedure, public, pass(self) :: depositJ
        procedure, public, pass(self) :: getError_tridiag_Ampere
        procedure, public, pass(self) :: construct_diagMatrix_Ampere
        procedure, public, pass(self) :: solveDivAmperePicard
        procedure, public, pass(self) :: solveDivAmpereAnderson
        procedure, public, pass(self) :: adaptiveSolveDivAmperePicard
        procedure, public, pass(self) :: adaptiveSolveDivAmpereAnderson
        procedure, public, pass(self) :: solveInitialPotential
        procedure, public, pass(self) :: moveParticles
        procedure, private, pass(self) :: construct_diagMatrix
        procedure, private, pass(self) :: particleSubStepInitialTau
        procedure, private, pass(self) :: particleSubStepTau
    end type

    interface potentialSolver
        module procedure :: potentialSolver_constructor
    end interface potentialSolver

contains

    type(potentialSolver) function potentialSolver_constructor(world, leftVoltage, rightVoltage, m_Anderson, Beta_k) result(self)
        ! Construct domain object, initialize grid, dx_dl, and nodeVol.
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: leftVoltage, rightVoltage, Beta_k
        integer(int32), intent(in) :: m_Anderson
        allocate(self % J(NumberXNodes-1), self % rho(NumberXNodes), self % phi(NumberXNodes), self % phi_f(NumberXNodes), self%a_tri(NumberXNodes-3), &
        self%b_tri(NumberXNodes-2), self%c_tri(NumberXNodes-3), self%particleChargeLoss(2, numberChargedParticles))
        call construct_diagMatrix(self, world)
        self % rho = 0
        self % J = 0
        self % phi = 0
        self % amountTimeSplits = 0
        self % Beta_k = Beta_k
        self % m_Anderson = m_Anderson
        self % coeff_left = 2/(world%dx_dl(1) + world%dx_dl(2))/world%dx_dl(1)
        self % coeff_right = 2/(world%dx_dl(size(world%dx_dl)-1) + world%dx_dl(size(world%dx_dl)))/world%dx_dl(size(world%dx_dl))
        self%iterNumPicard = 0
        self%iterNumAdaptiveSteps = 0
        self%iterNumParticle = 0
        self%particleEnergyLoss = 0.0d0
        self%particleChargeLoss = 0.0d0
        self%energyError = 0.0d0
        self%chargeError = 0.0d0
        if (world%boundaryConditions(1) == 1) self%phi(1) = leftVoltage
        if (world%boundaryConditions(NumberXNodes) == 1) self%phi(NumberXNodes) = rightVoltage
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
        self % a_tri = 2.0d0/(world%dx_dl(2:NumberXNodes-2) + world%dx_dl(3:)) / world%dx_dl(2:NumberXNodes-2)
        self % c_tri = 2.0d0/(world%dx_dl(1:NumberXNodes-3) + world%dx_dl(2:NumberXNodes-2))/world%dx_dl(2:NumberXNodes-2)
        self % b_tri = -2.0d0/(world%dx_dl(1:NumberXNodes-2) + world%dx_dl(2:))/world%dx_dl(1:NumberXNodes-2) - 2.0d0/(world%dx_dl(1:NumberXNodes-2) + world%dx_dl(2:))/world%dx_dl(2:NumberXNodes-1)

    end subroutine construct_diagMatrix

    subroutine construct_diagMatrix_Ampere(self, world)
        ! construct diagonal components for thomas algorithm, for Ampere (after initial Poisson)
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        self % a_tri = -1.0d0/world%dx_dl(2:NumberXNodes-2)
        self % c_tri = -1.0d0/world%dx_dl(2:NumberXNodes-2)
        self % b_tri = (world%dx_dl(1:NumberXNodes-2) + world%dx_dl(2:NumberXNodes-1))/ (world%dx_dl(1:NumberXNodes-2) * world%dx_dl(2:NumberXNodes-1))
        self % coeff_left = 1.0d0/world%dx_dl(1)
        self % coeff_right = 1.0d0/world%dx_dl(NumberXNodes-1)

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

    subroutine solve_tridiag_Ampere(self, world, del_t)
        ! Tridiagonal (Thomas algorithm) solver for Ampere
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t
        integer(int32) :: i, n !n size dependent on how many points are inside (not boundary), so how many in rhs equation
        real(real64) :: m, d(NumberXNodes - 2), cp(NumberXNodes - 3),dp(NumberXNodes- 2)
        n = NumberXNodes - 2

        d = (-self%J(2:) + self%J(1:n)) * del_t / eps_0 + arrayDiff(self%phi(1:n+1), n+1)/world%dx_dl(1:n) - arrayDiff(self%phi(2:), n+1)/world%dx_dl(2:)
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
        d = (-self%J(2:) + self%J(1:NumberXNodes-2)) * del_t / eps_0 + arrayDiff(self%phi(1:NumberXNodes-1), NumberXNodes-1)/world%dx_dl(1:NumberXNodes-2) - arrayDiff(self%phi(2:), NumberXNodes-1)/world%dx_dl(2:)
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
            res = 0.5 * eps_0 * SUM(arrayDiff(self%phi_f, NumberXNodes)**2 / world%dx_dl)
        else
            res = 0.5 * eps_0 * SUM(arrayDiff(self%phi, NumberXNodes)**2 / world%dx_dl)
        end if
    end function getTotalPE


    


    !----------------- Particle mover/sub-stepping procedures -------------------------------------------------------

    pure function getEField(self, l_cell, world) result(EField)
        !return EField of particle at logical position l_p (this is half distance), per particle since particle mover loop will be per particle
        class(potentialSolver), intent(in) :: self
        type(Domain), intent(in) :: world
        integer(int32), intent(in) :: l_cell
        real(real64) :: EField
        EField = (self%phi_f(l_cell) + self%phi(l_cell) - self%phi(l_cell+1) - self%phi_f(l_cell + 1)) / world%dx_dl(l_cell)/2
    end function getEField

    subroutine getl_BoundaryInitial(l_sub, v_sub, l_alongV, l_awayV)
        ! get point in l-space on boundary which is away or towards boundary based on velocity direction, when particle between nodes
        real(real64), intent(in out) :: l_alongV, l_awayV
        real(real64), intent(in) :: l_sub, v_sub
        if (v_sub > 0.0) then
            l_alongV = real(INT(l_sub) + 1, kind = real64)
            l_awayV = real(INT(l_sub), kind = real64)
        else
            l_alongV = real(INT(l_sub), kind = real64)
            l_awayV = real(INT(l_sub) + 1, kind = real64)
        end if

    end subroutine getl_BoundaryInitial


    subroutine particleSubStepInitialTau(self, world, part, l_sub, l_f, v_sub, v_f, del_tau, l_alongV, l_awayV, l_cell, a)
        ! Do initial substep, where particles start between nodes
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        type(Particle), intent(in) :: part
        integer(int32), intent(in) :: l_cell
        real(real64), intent(in out) :: l_f, v_f, del_tau, l_alongV, l_awayV, a
        real(real64), intent(in) :: l_sub, v_sub
        real(real64) :: c !rho_i(NumberXNodes), rho_f(NumberXNodes), gradJ(NumberXNodes-2), test(NumberXNodes-2), testConserv, 
        !integer(int32) :: k
        
        ! Particle first between nodes, so solve quadratic for that particle depending on conditions
        if ((v_sub/=0.0d0)) then ! make first case, since pretty much always likely to be the case (could just not have, assume always field, never have to check)
            c = (l_sub - l_alongV) * world%dx_dl(l_cell)
            if (a*v_sub > 0) then
                ! velocity and acceleration in same direction
                del_tau = 2.0d0 * ABS(c)/(ABS(v_sub) + SQRT(v_sub**2 - 4.0d0*a*c))
                l_f = l_alongV
                if (del_tau <= 0.0d0) then
                    print *, "del_tau is:", del_tau
                    print *, "a is:", a
                    print *, "c is:", c
                    print *, "l_sub is:", l_sub
                    print *, "v_sub is:", v_sub
                    stop "Have issue with del_tau for v,a in same direction"
                end if
            else if (v_sub**2 - 4.0d0*a*c > 0.0d0) then
                ! v and a opposite direction, but particle can still reach boundary along v
                del_tau = 2.0d0 * ABS(c)/(ABS(v_sub) + SQRT(v_sub**2 - 4.0d0*a*c))
                l_f = l_alongV
                if (del_tau <= 0.0d0) then
                    print *, "del_tau is:", del_tau
                    print *, "a is:", a
                    print *, "c is:", c
                    print *, "l_sub is:", l_sub
                    print *, "v_sub is:", v_sub
                    stop "Have issue with del_tau for v,a in opposite direction, but still reach boundary along v, initial tau"
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
                l_alongV = real(l_cell + 1, kind = real64)
            else
                l_alongV = real(l_cell, kind = real64)
            end if
            c = (l_sub - l_alongV) * world%dx_dl(l_cell)
            del_tau = SQRT(-c/a)
            l_f = l_alongV
            if (del_tau <= 0.0d0) then
                stop "Have issue with del_tau for v = 0"
            end if
        end if

    end subroutine particleSubStepInitialTau

    subroutine particleSubStepTau(self, world, part, l_sub, l_f, v_sub, v_f, del_tau, l_alongV, l_cell, a)
        ! Substeps, where particles start at nodes
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        type(Particle), intent(in) :: part
        integer(int32), intent(in) :: l_cell
        real(real64), intent(in out) :: l_f, v_f, del_tau, l_alongV, a
        real(real64), intent(in) :: l_sub, v_sub
        real(real64) :: c
        !integer(int32) :: k
        ! get index cell where field and dx_dl is evaluated
        ! Particle first between nodes, so solve quadratic for that particle depending on conditions
        c = (l_sub - l_alongV) * world%dx_dl(l_cell)
        if (a*v_sub > 0.0d0) then
            ! velocity and acceleration in same direction
            del_tau = 2.0d0 * ABS(c)/(ABS(v_sub) + SQRT(v_sub**2 - 4.0d0*a*c))
            l_f = l_alongV
            if (del_tau <= 0.0d0) then
                stop "Have issue with del_tau for v,a in same direction, ongoing substep"
            end if
        else if (v_sub**2 - 4.0d0*a*c > 0.0d0) then
            ! v and a opposite direction, but particle can still reach boundary along v
            del_tau = 2.0d0 * ABS(c)/(ABS(v_sub) + SQRT(v_sub**2 - 4.0d0*a*c))
            l_f = l_alongV
            if (del_tau <= 0.0d0) then
                print *, "del_tau is:", del_tau
                print *, "a is:", a
                print *, "c is:", c
                print *, "l_sub is:", l_sub
                print *, "v_sub is:", v_sub
                stop "Have issue with del_tau for v,a in opposite direction, but still reach boundary along v, ongoing substep"
            end if
        else
            ! v and a opposite direction, reverses back to initial position
            del_tau = ABS(v_sub)/ABS(a)
            l_f = l_sub
            if (del_tau <= 0.0d0) then
                stop "Have issue with del_tau for v,a in opposite direction, boundary opposite v"
            end if
        end if

    end subroutine particleSubStepTau

    subroutine depositJ(self, particleList, world, del_t)
        ! particle substepping procedure which deposits J, also used for general moving of particles
        ! boolDepositJ = true : deposit J
        ! boolDepositJ = false: only move particles (done after non-linear iteration), also delete particles for wall collisions (saves extra loop in collisions)
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        type(Particle), intent(in out) :: particleList(:)
        real(real64), intent(in) :: del_t
        !a and c correspond to quadratic equations | l_alongV is nearest integer boundary along velocity component, away is opposite
        real(real64) :: l_f, l_sub, v_sub, v_f, timePassed, del_tau, l_alongV, l_awayV, a
        integer(int32) :: subStepNum, j, i, delIdx, l_cell
        self%J = 0.0d0
        loopSpecies: do j = 1, numberChargedParticles
            delIdx = 0
            loopParticles: do i = 1, particleList(j)%N_p
                v_sub = particleList(j)%phaseSpace(2,i)
                l_sub = particleList(j)%phaseSpace(1,i)
                timePassed = 0.0d0
                subStepNum = 0
                l_cell = INT(l_sub)
                l_alongV = INT(l_sub) + 0.5d0 + SIGN(0.5d0, v_sub)
                l_awayV = INT(l_sub) + 0.5d0 - SIGN(0.5d0, v_sub)
                a = (particleList(j)%q / particleList(j)%mass / 2.0d0) * self%getEField(l_cell, world)
                call self%particleSubStepInitialTau(world, particleList(j), l_sub, l_f, v_sub, v_f, del_tau, l_alongV, l_awayV, l_cell, a)
                if (del_tau >= del_t) then
                    ! Add directly to J with no substep
                    l_f = v_sub * del_t / world%dx_dl(l_cell) + (a/ world%dx_dl(l_cell)) * del_t**2 + l_sub
                    v_f = 2.0d0 * (l_f - l_sub) * world%dx_dl(l_cell) / del_t - v_sub
                    timePassed = del_t
                    if (INT(l_f) /= l_cell) then
                        stop "l_f has crossed boundary when condition says is shouldn't have any substeps"
                    end if
                    self%J(l_cell) = self%J(l_cell) + particleList(j)%w_p * particleList(j)%q * (v_f + v_sub)/2.0d0/world%dx_dl(l_cell)
                else
                    v_f = 2.0d0 * (l_f - l_sub) * world%dx_dl(l_cell) / del_tau - v_sub
                    timePassed = timePassed + del_tau
                    self%J(l_cell) = self%J(l_cell) + particleList(j)%w_p * particleList(j)%q * (v_f + v_sub)*del_tau/2.0d0/world%dx_dl(l_cell)/del_t
                    if (MOD(l_f, 1.0d0) /= 0.0d0) then
                        print *, l_f
                        stop "l_f is not integer after subStep"
                    end if
                    firstStep: SELECT CASE (world%boundaryConditions(INT(l_f)))
                    CASE(0)
                        continue
                    CASE(1)
                       timePassed = del_t
                    CASE(3)
                        l_f = ABS(l_f - real(NumberXNodes, kind = real64) - 1.0d0)
                    CASE default
                        print *, "Case does not exist in first substep, depositJ"
                        stop
                    END SELECT firstStep
                    ! now final position/velocity becomes next starting position/velocity
                    l_sub = l_f
                    v_sub = v_f
                end if
                do while((timePassed < del_t))
                    l_alongV = l_sub + SIGN(1.0d0, v_sub)
                    l_cell = INT(l_sub + SIGN(0.5d0, v_sub))
                    a = (particleList(j)%q / particleList(j)%mass / 2.0d0) * self%getEField(l_cell, world)
                    call self%particleSubStepTau(world, particleList(j), l_sub, l_f, v_sub, v_f, del_tau, l_alongV, l_cell, a)
                    if (del_tau >= del_t-timePassed) then
                        ! Add directly to J with no substep
                        l_f = v_sub * (del_t-timePassed) / world%dx_dl(l_cell) + (a/ world%dx_dl(l_cell)) * (del_t-timePassed)**2 + l_sub
                        v_f = 2.0d0 * (l_f - l_sub) * world%dx_dl(l_cell) / (del_t - timePassed) - v_sub
                        if (ABS(l_f - l_sub) >= 1) then
                            print *, "l_sub is:", l_sub
                            print *, "v_sub is:", v_sub
                            print *, "a is:", a
                            print *, "relative error is:", a*(del_t-timePassed)/v_sub
                            print *, "c is:", (l_sub - l_alongV) * world%dx_dl(l_cell)
                            print *, "l_f is:", l_f
                            print *, "v_f is:", v_f
                            print *, "del_tau is:", del_tau
                            print *, "remaining is:", del_t - timePassed
                            print *, "del_tau with other inverse formula:", 2.0d0 * ABS((l_sub - l_alongV) * world%dx_dl(l_cell))/(SQRT(v_sub**2 - 4.0d0*a*(l_sub - l_alongV) * world%dx_dl(l_cell)) + ABS(v_sub))
                            stop "l_f has crossed boundary when condition says is shouldn't have any substeps"
                        end if
                        self%J(l_cell) = self%J(l_cell) + particleList(j)%w_p * particleList(j)%q * (v_f + v_sub)*(del_t - timePassed)/2.0d0/world%dx_dl(l_cell)/del_t
                        timePassed = del_t
                        
                    else
                        v_f = 2.0d0 * (l_f - l_sub) * world%dx_dl(l_cell) / del_tau - v_sub
                        timePassed = timePassed + del_tau
                        self%J(l_cell) = self%J(l_cell) + particleList(j)%w_p * particleList(j)%q * (v_f + v_sub)*del_tau/2.0d0/world%dx_dl(l_cell)/del_t
                        if (MOD(l_f, 1.0d0) /= 0.0d0) then
                            print *, l_f
                            stop "l_f is not integer after subStep"
                        end if
                        subStep: SELECT CASE (world%boundaryConditions(INT(l_f)))
                        CASE(0)
                            continue
                        CASE(1)
                            exit
                        CASE(3)
                            l_f = ABS(l_f - real(NumberXNodes, kind = real64) - 1.0d0)
                        CASE default
                            print *, "Case does not exist in ongoing substep, depositJ"
                            stop
                        END SELECT subStep
                        ! now final position/velocity becomes next starting position/velocity
                        l_sub = l_f
                        v_sub = v_f
                    end if
                    subStepNum = subStepNum + 1
                end do
                if ((l_f < 1) .or. (l_f > NumberXNodes)) then
                    stop "Have particles travelling oremainDel_tutside domain!"
                end if
            end do loopParticles
            
        end do loopSpecies
        
    end subroutine depositJ


    ! -------------------------------------------- Particle mover without boolean checks for depositing J ------------------------------------------------------------
   
    subroutine moveParticles(self, particleList, world, del_t)
        ! particle mover to avoid the boolean checks which mostly don't happen when depositing J
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        type(Particle), intent(in out) :: particleList(:)
        real(real64), intent(in) :: del_t
        logical :: delParticle
        !a and c correspond to quadratic equations | l_alongV is nearest integer boundary along velocity component, away is opposite
        real(real64) :: l_f, l_sub, v_sub, v_f, timePassed, del_tau, l_alongV, l_awayV, a
        integer(int32) :: subStepNum, j, i, delIdx, l_cell
        delParticle = .false.
        loopSpecies: do j = 1, numberChargedParticles
            delIdx = 0
            loopParticles: do i = 1, particleList(j)%N_p
                v_sub = particleList(j)%phaseSpace(2,i)
                l_sub = particleList(j)%phaseSpace(1,i)
                timePassed = 0
                subStepNum = 0
                l_cell = INT(l_sub)
                l_alongV = INT(l_sub) + 0.5d0 + SIGN(0.5d0, v_sub)
                l_awayV = INT(l_sub) + 0.5d0 - SIGN(0.5d0, v_sub)
                a = (particleList(j)%q / particleList(j)%mass / 2.0d0) * self%getEField(l_cell, world)
                call self%particleSubStepInitialTau(world, particleList(j), l_sub, l_f, v_sub, v_f, del_tau, l_alongV, l_awayV, l_cell, a)
                if (del_tau >= del_t) then
                    ! Add directly to J with no substep
                    l_f = v_sub * del_t / world%dx_dl(l_cell) + (a/ world%dx_dl(l_cell)) * del_t**2 + l_sub
                    v_f = 2.0d0 * (l_f - l_sub) * world%dx_dl(l_cell) / del_t - v_sub
                    particleList(j)%phaseSpace(1, i-delIdx) = l_f
                    particleList(j)%phaseSpace(2,i-delIdx) = v_f
                    particleList(j)%phaseSpace(3:4, i-delIdx) = particleList(j)%phaseSpace(3:4, i)
                    timePassed = del_t
                    if (INT(l_f) /= l_cell) then
                        stop "l_f has crossed boundary when condition says is shouldn't have any substeps"
                    end if
                    
                else
                    v_f = 2.0d0 * (l_f - l_sub) * world%dx_dl(l_cell) / del_tau - v_sub
                    timePassed = timePassed + del_tau
                    if (MOD(l_f, 1.0d0) /= 0.0d0) then
                        print *, l_f
                        stop "l_f is not integer after subStep"
                    end if
                    firstStep: SELECT CASE (world%boundaryConditions(INT(l_f)))
                    CASE(0)
                        continue
                    CASE(1)
                        timePassed = del_t
                        delIdx = delIdx + 1
                        self%particleEnergyLoss = self%particleEnergyLoss + particleList(j)%w_p * (v_f**2 + SUM(particleList(j)%phaseSpace(3:4, i)**2)) * particleList(j)%mass * 0.5d0 !J/m^2 in 1D
                        if (l_f == 1) then
                            self%particleChargeLoss(1, j) = self%particleChargeLoss(1, j) + particleList(j)%q * particleList(j)%w_p !C/m^2 in 1D
                        else
                            self%particleChargeLoss(2, j) = self%particleChargeLoss(2, j) + particleList(j)%q * particleList(j)%w_p !C/m^2 in 1D
                        end if
                    CASE(3)
                        l_f = ABS(l_f - real(NumberXNodes, kind = real64) - 1.0d0)
                    CASE default
                        print *, "Case does not exist in first substep, depositJ"
                        stop
                    END SELECT firstStep
                    ! now final position/velocity becomes next starting position/velocity
                    l_sub = l_f
                    v_sub = v_f
                end if
                do while((timePassed < del_t))
                    l_alongV = l_sub + SIGN(1.0d0, v_sub)
                    l_cell = INT(l_sub + SIGN(0.5d0, v_sub))
                    a = (particleList(j)%q / particleList(j)%mass / 2.0d0) * self%getEField(l_cell, world)
                    call self%particleSubStepTau(world, particleList(j), l_sub, l_f, v_sub, v_f, del_tau, l_alongV, l_cell, a)
                    if (del_tau >= del_t-timePassed) then
                        ! Add directly to J with no substep
                        l_f = v_sub * (del_t-timePassed) / world%dx_dl(l_cell) + (a/ world%dx_dl(l_cell)) * (del_t-timePassed)**2 + l_sub
                        v_f = 2.0d0 * (l_f - l_sub) * world%dx_dl(l_cell) / (del_t - timePassed) - v_sub
                        if (ABS(l_f - l_sub) >= 1) then
                            stop "l_f has crossed boundary when condition says is shouldn't have any substeps"
                        end if
                        particleList(j)%phaseSpace(1, i-delIdx) = l_f
                        particleList(j)%phaseSpace(2,i-delIdx) = v_f
                        particleList(j)%phaseSpace(3:4, i-delIdx) = particleList(j)%phaseSpace(3:4, i)
                        timePassed = del_t
                        
                    else
                        v_f = 2.0d0 * (l_f - l_sub) * world%dx_dl(l_cell) / del_tau - v_sub
                        timePassed = timePassed + del_tau
                        if (MOD(l_f, 1.0d0) /= 0.0d0) then
                            print *, l_f
                            stop "l_f is not integer after subStep"
                        end if
                        subStep: SELECT CASE (world%boundaryConditions(INT(l_f)))
                        CASE(0)
                            continue
                        CASE(1)
                            delIdx = delIdx + 1
                            self%particleEnergyLoss = self%particleEnergyLoss + particleList(j)%w_p * (v_f**2 + SUM(particleList(j)%phaseSpace(3:4, i)**2)) * particleList(j)%mass * 0.5d0 !J/m^2 in 1D
                            if (l_f == 1) then
                                self%particleChargeLoss(1, j) = self%particleChargeLoss(1, j) + particleList(j)%q * particleList(j)%w_p !C/m^2 in 1D
                            else
                                self%particleChargeLoss(2, j) = self%particleChargeLoss(2, j) + particleList(j)%q * particleList(j)%w_p !C/m^2 in 1D
                            end if
                            exit
                        CASE(3)
                            l_f = ABS(l_f - real(NumberXNodes, kind = real64) - 1.0d0)
                        CASE default
                            print *, "Case does not exist in ongoing substep, depositJ"
                            stop
                        END SELECT subStep
                        ! now final position/velocity becomes next starting position/velocity
                        l_sub = l_f
                        v_sub = v_f
                    end if
                    subStepNum = subStepNum + 1
                end do
                if ((l_f < 1) .or. (l_f > NumberXNodes)) then
                    stop "Have particles travelling oremainDel_tutside domain!"
                end if
                
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

    !--------------------------- non-linear solver for time step using divergence of ampere ------------------------

    subroutine solveDivAmperePicard(self, particleList, world, del_t, maxIter, eps_r)
        ! Solve for divergence of ampere using picard iterations
        class(potentialSolver), intent(in out) :: self
        type(Particle), intent(in out) :: particleList(:)
        type(Domain), intent(in) :: world
        integer(int32), intent(in) :: maxIter
        real(real64), intent(in) :: del_t, eps_r
        real(real64) :: errorCurrent, errorInitial
        integer(int32) :: i
        call self%depositJ(particleList, world, del_t)
        errorInitial = self%getError_tridiag_Ampere(world, del_t)
        do i = 1, maxIter
            call self%solve_tridiag_Ampere(world, del_t)
            call self%depositJ(particleList, world, del_t)
            errorCurrent = self%getError_tridiag_Ampere(world, del_t)
            if (i > 2) then
                if (errorCurrent < eps_r*errorInitial) then
                    call self%moveParticles(particleList, world, del_t)
                    self%phi = self%phi_f
                    exit
                end if
            end if
        end do
        self%iterNumPicard = i-1
        

    end subroutine solveDivAmperePicard

    subroutine solveDivAmpereAnderson(self, particleList, world, del_t, maxIter, eps_r)
        ! Solve for divergence of ampere using picard iterations
        class(potentialSolver), intent(in out) :: self
        type(Particle), intent(in out) :: particleList(:)
        type(Domain), intent(in) :: world
        integer(int32), intent(in) :: maxIter
        real(real64), intent(in) :: del_t, eps_r
        real(real64) :: initialR, sumPastResiduals, initialNorm
        real(real64) :: Residual_k(NumberXNodes-2, self%m_Anderson+1), phi_k(NumberXNodes-2, self%m_Anderson+1), fitMat(NumberXNodes-2, self%m_Anderson), alpha(NumberXNodes-2)
        integer(int32) :: lwork, work(NumberXNodes -2 + self%m_Anderson)
        integer(int32) :: i, j, index, m_k, info
        lwork=(NumberXNodes -2)+ self%m_Anderson

        phi_k(:,1) = self%phi(2:NumberXNodes-1)
        initialNorm = SQRT(SUM(self%phi(2:NumberXNodes-1)**2))
        call self%depositJ(particleList, world, del_t)
        call self%solve_tridiag_Ampere(world, del_t)
        phi_k(:,2) = self%phi_f(2:NumberXNodes-1)
        initialR = SQRT(SUM((self%phi_f(2:NumberXNodes-1) - phi_k(:,1))**2))
        Residual_k(:,1) = phi_k(:,2) - phi_k(:,1)
        do i = 1, maxIter
            index = MODULO(i, self%m_Anderson+1) + 1
            m_k = MIN(i, self%m_Anderson)
            call self%depositJ(particleList, world, del_t)
            if (SQRT(SUM((self%phi_f(2:NumberXNodes-1) - phi_k(:,MODULO(i-1, self%m_Anderson+1) + 1))**2)) < eps_r*(initialR + initialNorm)) then
                call self%moveParticles(particleList, world, del_t)
                self%phi = self%phi_f
                exit
            end if
            call self%solve_tridiag_Ampere(world, del_t)
            Residual_k(:, index) = self%phi_f(2:NumberXNodes-1) - phi_k(:,index)
            if (i > self%m_Anderson) then
                if (self%m_anderson > 1) then
                    sumPastResiduals = SUM(Residual_k(:, MODULO(i-1, self%m_Anderson+1) + 1)**2) &
                    + SUM(Residual_k(:, MODULO(i-2, self%m_Anderson+1) + 1)**2)
                    sumPastResiduals = sumPastResiduals/2.0d0
                else
                    sumPastResiduals = SUM(Residual_k(:, MODULO(i-1, self%m_Anderson+1) + 1)**2)
                end if
                if (SUM(Residual_k(:, index)**2) > sumPastResiduals) then
                    self%iterNumPicard = maxIter
                    goto 75
                end if
            end if
            do j = 0, m_k-1
                fitMat(:,j+1) = Residual_k(:, MODULO(i - m_k + j, self%m_Anderson+1) + 1) - Residual_k(:, index)
            end do
            alpha = -Residual_k(:, index)
            call dgels('N', NumberXNodes-2, m_k, 1, fitMat(:, 1:m_k), NumberXNodes-2, alpha, NumberXNodes-2, work, lwork, info)
            if (info /= 0) then
                print *, "Issue with minimization procedure dgels in Anderson Acceleration!"
                stop
            end if
            alpha(m_k+1) = 1.0d0 - SUM(alpha(1:m_k)) 
            phi_k(:, MODULO(i+1, self%m_Anderson+1) + 1) = alpha(1) * (self%Beta_k*Residual_k(:, MODULO(i-m_k, self%m_Anderson+1) + 1) + phi_k(:, MODULO(i-m_k,self%m_Anderson+1) + 1))
            do j=1, m_k
                phi_k(:, MODULO(i+1, self%m_Anderson+1) + 1) = phi_k(:, MODULO(i+1, self%m_Anderson+1) + 1) + alpha(j + 1) * (self%Beta_k*Residual_k(:, MODULO(i-m_k + j, self%m_Anderson+1) + 1) + phi_k(:, MODULO(i-m_k + j, self%m_Anderson+1) + 1))
            end do
            self%phi_f(2:NumberXNodes-1) = phi_k(:, MODULO(i+1, self%m_Anderson+1) + 1)
    
        end do
        self%iterNumPicard = i-1
        75 continue
        

    end subroutine solveDivAmpereAnderson

    subroutine adaptiveSolveDivAmperePicard(self, particleList, world, del_t, maxIter, eps_r)
        ! Solve for divergence of ampere's law with picard
        ! cut del_t in two if non-convergence after maxIter, repeat until convergence
        class(potentialSolver), intent(in out) :: self
        type(Particle), intent(in out) :: particleList(:)
        type(Domain), intent(in) :: world
        integer(int32), intent(in) :: maxIter
        real(real64), intent(in) :: del_t, eps_r
        real(real64) :: remainDel_t, currDel_t, adaptiveJ(NumberXNodes-1)
        call self%solveDivAmperePicard(particleList, world, del_t, maxIter, eps_r)
        remainDel_t = del_t  
        if (self%iterNumPicard == maxIter) then
            do while (self%iterNumPicard == maxIter)
                print *, "Reducing time step adaptively"
                self%iterNumAdaptiveSteps = 0
                currDel_t = remainDel_t
                adaptiveJ = 0.0d0
                do while (self%iterNumPicard == maxIter)
                    currDel_t = currDel_t/2.0d0
                    self%iterNumAdaptiveSteps = self%iterNumAdaptiveSteps + 1
                    if (self%iterNumAdaptiveSteps > 4) then
                        stop "ALREADY REDUCED TIME STEP MORE THAN 3 TIMES, REDUCE INITIAL TIME STEP!!!"
                    end if
                    call self%solveDivAmperePicard(particleList, world, currDel_t, maxIter, eps_r)   
                end do
                adaptiveJ = adaptiveJ + self%J * currDel_t/del_t
                remainDel_t = remainDel_t - currDel_t 
                call self%solveDivAmperePicard(particleList, world, remainDel_t, maxIter, eps_r)
            end do
            adaptiveJ = adaptiveJ + self%J * remainDel_t/del_t
            self%J = adaptiveJ
        end if
    end subroutine adaptiveSolveDivAmperePicard

    subroutine adaptiveSolveDivAmpereAnderson(self, particleList, world, del_t, maxIter, eps_r)
        ! Solve for divergence of ampere's law with picard
        ! cut del_t in two if non-convergence after maxIter, repeat until convergence
        class(potentialSolver), intent(in out) :: self
        type(Particle), intent(in out) :: particleList(:)
        type(Domain), intent(in) :: world
        integer(int32), intent(in) :: maxIter
        real(real64), intent(in) :: del_t, eps_r
        real(real64) :: remainDel_t, currDel_t, adaptiveJ(NumberXNodes-1)
        call self%solveDivAmpereAnderson(particleList, world, del_t, maxIter, eps_r)
        remainDel_t = del_t  
        if (self%iterNumPicard == maxIter) then
            do while (self%iterNumPicard == maxIter)
                self%iterNumAdaptiveSteps = 0
                currDel_t = remainDel_t
                adaptiveJ = 0.0d0
                do while (self%iterNumPicard == maxIter)
                    self%amountTimeSplits = self%amountTimeSplits + 1
                    currDel_t = currDel_t/2.0d0
                    self%iterNumAdaptiveSteps = self%iterNumAdaptiveSteps + 1
                    if (self%iterNumAdaptiveSteps > 4) then
                        stop "ALREADY REDUCED TIME STEP MORE THAN 3 TIMES, REDUCE INITIAL TIME STEP!!!"
                    end if
                    call self%solveDivAmpereAnderson(particleList, world, currDel_t, maxIter, eps_r)   
                end do
                adaptiveJ = adaptiveJ + self%J * currDel_t/del_t
                remainDel_t = remainDel_t - currDel_t 
                call self%solveDivAmpereAnderson(particleList, world, remainDel_t, maxIter, eps_r)
            end do
            adaptiveJ = adaptiveJ + self%J * remainDel_t/del_t
            self%J = adaptiveJ
        end if
    end subroutine adaptiveSolveDivAmpereAnderson



end module mod_potentialSolver