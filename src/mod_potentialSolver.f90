module mod_potentialSolver
    use iso_fortran_env, only: int32, real64
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
    public :: potSolver

    type :: potSolver
        real(real64), allocatable :: phi(:), J(:), rho(:), phi_f(:) !phi_f is final phi, will likely need to store two arrays for phi, can't be avoided
        real(real64) :: phi_left, phi_right, energyError, chargeError, particleEnergyLoss, particleChargeLoss
        integer(int32) :: iterNumPicard, iterNumParticle, iterNumAdaptiveSteps
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
        procedure, public, pass(self) :: adaptiveSolveDivAmperePicard
        procedure, public, pass(self) :: solveInitialPotential
        procedure, public, pass(self) :: moveParticles
        procedure, private, pass(self) :: construct_diagMatrix
    end type

    interface potSolver
        module procedure :: potSolver_constructor
    end interface potSolver

contains

    type(potSolver) function potSolver_constructor(world) result(self)
        ! Construct domain object, initialize grid, dx_dl, and nodeVol.
        type(Domain), intent(in) :: world
        allocate(self % J(n_x-1), self % rho(n_x), self % phi(n_x), self % phi_f(n_x), self%a_tri(n_x-3), &
        self%b_tri(n_x-2), self%c_tri(n_x-3))
        call construct_diagMatrix(self, world)
        self % rho = 0
        self % J = 0
        self % phi = 0
        self % phi_f = 0
        self % phi_left = 0.0d0
        self % phi_right = 0.0d0
        self % coeff_left = 2/(world%dx_dl(1) + world%dx_dl(2))/world%dx_dl(1)
        self % coeff_right = 2/(world%dx_dl(size(world%dx_dl)-1) + world%dx_dl(size(world%dx_dl)))/world%dx_dl(size(world%dx_dl))
        self%iterNumPicard = 0
        self%iterNumAdaptiveSteps = 0
        self%iterNumParticle = 0
        self%particleEnergyLoss = 0.0d0
        self%particleChargeLoss = 0.0d0
        self%energyError = 0.0d0
        self%chargeError = 0.0d0

    end function potSolver_constructor

    subroutine depositRho(self, particleList, world) 
        class(potSolver), intent(in out) :: self
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
    end subroutine depositRho

    subroutine construct_diagMatrix(self, world)
        ! construct diagonal components for thomas algorithm
        class(potSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        self % a_tri = 2.0d0/(world%dx_dl(2:n_x-2) + world%dx_dl(3:)) / world%dx_dl(2:n_x-2)
        self % c_tri = 2.0d0/(world%dx_dl(1:n_x-3) + world%dx_dl(2:n_x-2))/world%dx_dl(2:n_x-2)
        self % b_tri = -2.0d0/(world%dx_dl(1:n_x-2) + world%dx_dl(2:))/world%dx_dl(1:n_x-2) - 2.0d0/(world%dx_dl(1:n_x-2) + world%dx_dl(2:))/world%dx_dl(2:n_x-1)

    end subroutine construct_diagMatrix

    subroutine construct_diagMatrix_Ampere(self, world)
        ! construct diagonal components for thomas algorithm, for Ampere (after initial Poisson)
        class(potSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        self % a_tri = -1.0d0/world%dx_dl(2:n_x-2)
        self % c_tri = -1.0d0/world%dx_dl(2:n_x-2)
        self % b_tri = (world%dx_dl(1:n_x-2) + world%dx_dl(2:n_x-1))/ (world%dx_dl(1:n_x-2) * world%dx_dl(2:n_x-1))
        self % coeff_left = 1.0d0/world%dx_dl(1)
        self % coeff_right = 1.0d0/world%dx_dl(n_x-1)

    end subroutine construct_diagMatrix_Ampere

    subroutine solve_tridiag_Poisson(self)
        ! Tridiagonal (Thomas algorithm) solver for initial Poisson
        class(potSolver), intent(in out) :: self
        integer(int32) :: i, n !n size dependent on how many points are boundary conditions and thus moved to rhs equation
        real(real64) :: m, d(size(self%phi) - 2), cp(size(self%phi) - 3),dp(size(self%phi) - 2)
        n = n_x - 2

        d = -self%rho(2:n_x-1) / eps_0
        d(1) = d(1) - 2 * self%phi_left * self%coeff_left
        d(n) = d(n) - 2 * self%phi_right * self%coeff_right
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
        class(potSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t
        integer(int32) :: i, n !n size dependent on how many points are inside (not boundary), so how many in rhs equation
        real(real64) :: m, d(n_x - 2), cp(n_x - 3),dp(n_x- 2)
        n = n_x - 2

        d = (-self%J(2:) + self%J(1:n)) * del_t / eps_0 + arrayDiff(self%phi(1:n+1))/world%dx_dl(1:n) - arrayDiff(self%phi(2:))/world%dx_dl(2:)
        d(1) = d(1) + self%phi_left * self%coeff_left
        d(n) = d(n) + self%phi_right * self%coeff_right
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
        class(potSolver), intent(in) :: self
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t
        integer(int32) :: i
        real(real64) :: Ax(n_x-2), d(n_x-2), res
        Ax(1) = self%b_tri(1)*self%phi_f(2) + self%c_tri(1) * self%phi_f(3)
        do i=2, n_x-3
            Ax(i) = self%b_tri(i)*self%phi_f(i+1) + self%c_tri(i) * self%phi_f(i+2) + self%a_tri(i-1) * self%phi_f(i)
        end do
        Ax(n_x-2) = self%b_tri(n_x-2)*self%phi_f(n_x-1) + self%a_tri(n_x-3) * self%phi_f(n_x-2)
        d = (-self%J(2:) + self%J(1:n_x-2)) * del_t / eps_0 + arrayDiff(self%phi(1:n_x-1))/world%dx_dl(1:n_x-2) - arrayDiff(self%phi(2:))/world%dx_dl(2:)
        d(1) = d(1) + self%phi_left * self%coeff_left
        d(n_x-2) = d(n_x-2) + self%phi_right * self%coeff_right
        !res = SQRT(SUM(((Ax- d)/self%minEField)**2)/(n_x-2))
        res = SQRT(SUM((Ax- d)**2))

    end function getError_tridiag_Ampere

    function getTotalPE(self, world, future) result(res)
        ! Get energy in electric fields, future true, then derive from phi_f, otherwise phi
        ! In 1D is J/m^2
        class(potSolver), intent(in) :: self
        type(Domain), intent(in) :: world
        logical :: future
        real(real64) :: res
        if (future) then
            res = 0.5 * eps_0 * SUM(arrayDiff(self%phi_f)**2 / world%dx_dl)
        else
            res = 0.5 * eps_0 * SUM(arrayDiff(self%phi)**2 / world%dx_dl)
        end if
    end function getTotalPE


    


    !----------------- Particle mover/sub-stepping procedures -------------------------------------------------------

    pure function getEField(self, l_p, world) result(EField)
        !return EField of particle at logical position l_p (this is half distance), per particle since particle mover loop will be per particle
        class(potSolver), intent(in) :: self
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: l_p
        real(real64) :: EField
        if (l_p < 1 .or. l_p >= n_x) then
            EField = 0.0
        else
            EField = (self%phi_f(INT(l_p)) + self%phi(INT(l_p)) - self%phi(INT(l_p)+1) - self%phi_f(INT(l_p) + 1)) / world%dx_dl(INT(l_p))/2
        end if
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


    subroutine particleSubStepInitial(solver, world, part, l_sub, l_f, v_sub, v_f, timePassed, del_tau, del_t, l_alongV, l_awayV, boolDepositJ)
        ! Do initial substep, where particles start between nodes
        type(potSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        type(Particle), intent(in) :: part
        real(real64), intent(in out) :: l_sub, l_f, v_sub, v_f, timePassed, del_tau, l_alongV, l_awayV
        real(real64), intent(in) :: del_t
        logical, intent(in) :: boolDepositJ
        real(real64) :: a, c !rho_i(n_x), rho_f(n_x), gradJ(n_x-2), test(n_x-2), testConserv, 
        !integer(int32) :: k
        
        a = (part%q / part%mass / 2) * solver%getEField(l_sub, world)
        ! Particle first between nodes, so solve quadratic for that particle depending on conditions
        if ((a/=0) .and. (v_sub/=0)) then ! make first case, since pretty much always likely to be the case (could just not have, assume always field, never have to check)
            c = (l_sub - l_alongV) * world%dx_dl(INT(l_sub))
            if (a*v_sub > 0) then
                ! velocity and acceleration in same direction
                del_tau = (-ABS(v_sub) + SQRT(v_sub**2 - 4*a*c))/2/ABS(a)
                l_f = l_alongV
                if (del_tau <= 0) then
                    stop "Have issue with del_tau for v,a in same direction"
                end if
            else if (v_sub**2 - 4*a*c > 0) then
                ! v and a opposite direction, but particle can still reach boundary along v
                del_tau = (ABS(v_sub) - SQRT(v_sub**2 - 4*a*c))/2/ABS(a)
                l_f = l_alongV
                if (del_tau <= 0) then
                    stop "Have issue with del_tau for v,a in opposite direction, but still reach boundary along v"
                end if
            else
                ! v and a opposite direction, boundary opposite direction of v
                c = (l_sub - l_awayV) * world%dx_dl(INT(l_sub))
                del_tau = (ABS(v_sub) + SQRT(v_sub**2 - 4*a*c))/2/ABS(a)
                l_f = l_awayV
                if (del_tau <= 0) then
                    stop "Have issue with del_tau for v,a in opposite direction, boundary opposite v"
                end if
            end if
        else if ((a == 0.0) .and. (v_sub /= 0.0)) then
            !Free particle drift
            del_tau = (l_alongV - l_sub) * world%dx_dl(INT(l_sub))/v_sub
            v_f = (l_alongV - l_sub) * world%dx_dl(INT(l_sub)) / del_tau
            l_f = l_alongV
            if (del_tau <= 0) then
                stop "Have issue with del_tau for a = 0"
            end if

        else
            ! only in direction of field: USE l_alongV AS BOUNDARY ALONG DIRECTION OF a SINCE VELOCITY = 0!!!
            if (a > 0) then
                l_alongV = real(INT(l_sub) + 1, kind = real64)
            else
                l_alongV = real(INT(l_sub), kind = real64)
            end if
            c = (l_sub - l_alongV) * world%dx_dl(INT(l_sub))
            del_tau = SQRT(-c/a)
            l_f = l_alongV
            if (del_tau <= 0) then
                stop "Have issue with del_tau for v = 0"
            end if
        end if

        if (del_tau == del_t) then
            stop "Somehow del_tau didn't change in one of the conditions"
        end if

        if (del_tau >= del_t) then
            ! Add directly to J with no substep
            l_f = v_sub * del_t / world%dx_dl(INT(l_sub)) + (a/ world%dx_dl(INT(l_sub))) * del_t**2 + l_sub
            v_f = 2 * (l_f - l_sub) * world%dx_dl(INT(l_sub)) / del_t - v_sub
            timePassed = del_t
            ! if ((ABS((v_f - 2.0d0*a*del_t - v_sub)/v_sub) > 1e-3)) then
            !     print *, "WARNING: Kinematic equation doesn't match for no substep, subStepNum = 0"
            !     print *, "Particle is:", part%name
            !     print *, "v_sub is:", v_subdelParticle = .false.
            !     print *, "l_sub is:", l_sub
            !     print *, "l_f is:", l_f
            !     print *, "v_f is:", v_f
            !     print *, "error is:", ABS((v_f - 2.0d0*a*del_t - v_sub)/v_sub)
            if (INT(l_f) /= INT(l_sub)) then
                stop "l_f has crossed boundary when condition says is shouldn't have any substeps"
            end if
            if (boolDepositJ) solver%J(INT(l_sub)) = solver%J(INT(l_sub)) + part%w_p * part%q * (v_f + v_sub)/2/world%dx_dl(INT(l_sub))
            
        else
            v_f = 2 * (l_f - l_sub) * world%dx_dl(INT(l_sub)) / del_tau - v_sub
            timePassed = timePassed + del_tau
            if (boolDepositJ) solver%J(INT(l_sub)) = solver%J(INT(l_sub)) + part%w_p * part%q * (v_f + v_sub)*del_tau/2/world%dx_dl(INT(l_sub))/del_t
            if (MOD(l_f, 1.0) /= 0.0) then
                print *, l_f
                stop "l_f is not integer after subStep"
            end if
            ! now final position/velocity becomes next starting position/velocity
            l_sub = l_f
            v_sub = v_f
    
        end if

    end subroutine particleSubStepInitial

    subroutine particleSubStep(solver, world, part, l_sub, l_f, v_sub, v_f, timePassed, del_tau, del_t, l_alongV, l_cell, boolDepositJ)
        ! Substeps, where particles start at nodes
        type(potSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        type(Particle), intent(in) :: part
        real(real64), intent(in out) :: l_sub, l_f, v_sub, v_f, timePassed, del_tau, l_alongV, l_cell
        real(real64), intent(in) :: del_t
        logical, intent(in) :: boolDepositJ
        real(real64) :: a, c !rho_i(n_x), rho_f(n_x), gradJ(n_x-2), test(n_x-2), testConserv, 
        !integer(int32) :: k
        ! get index cell where field and dx_dl is evaluated
        if (v_sub > 0) then
            l_cell = l_sub
        else
            l_cell = l_sub - 1.0d0
        end if
        a = (part%q / part%mass / 2) * solver%getEField(l_cell, world)
        ! Particle first between nodes, so solve quadratic for that particle depending on conditions
        if ((a/=0)) then ! make first case, since pretty much always likely to be the case (could just not have, assume always field exists, never have to check)
            c = (l_sub - l_alongV) * world%dx_dl(INT(l_cell))
            if (a*v_sub > 0) then
                ! velocity and acceleration in same direction
                del_tau = (-ABS(v_sub) + SQRT(v_sub**2 - 4*a*c))/2/ABS(a)
                l_f = l_alongV
                if (del_tau <= 0) then
                    stop "Have issue with del_tau for v,a in same direction"
                end if
            else if (v_sub**2 - 4*a*c > 0) then
                ! v and a opposite direction, but particle can still reach boundary along v
                del_tau = (ABS(v_sub) - SQRT(v_sub**2 - 4*a*c))/2/ABS(a)
                l_f = l_alongV
                if (del_tau <= 0) then
                    stop "Have issue with del_tau for v,a in opposite direction, but still reach boundary along v"
                end if
            else
                ! v and a opposite direction, reverses back to initial position
                del_tau = ABS(v_sub)/ABS(a)
                l_f = l_sub
                if (del_tau <= 0) then
                    stop "Have issue with del_tau for v,a in opposite direction, boundary opposite v"
                end if
            end if
        else
            !Free particle drift
            del_tau = (l_alongV - l_sub) * world%dx_dl(INT(l_cell))/v_sub
            v_f = (l_alongV - l_sub) * world%dx_dl(INT(l_cell)) / del_tau
            l_f = l_alongV
            if (del_tau <= 0) then
                stop "Have issue with del_tau for a = 0"
            end if
        end if


        if (del_tau >= del_t-timePassed) then
            ! Add directly to J with no substep
            l_f = v_sub * (del_t - timePassed) / world%dx_dl(INT(l_cell)) + (a/ world%dx_dl(INT(l_cell))) * (del_t - timePassed)**2 + l_sub
            v_f = 2 * (l_f - l_sub) * world%dx_dl(INT(l_cell)) / (del_t - timePassed) - v_sub
            ! if (ABS((v_f - 2*a*(del_t - timePassed) - v_sub)/v_sub) > 1e-3) then
            !     print *, "WARNING: Kinematic equation doesn't match for no substep, subStepNum > 0"
            !     print *, "Particle is:", part%name
            !     print *, "v_sub is:", v_sub
            !     print *, "l_sub is:", l_sub
            !     print *, "l_f is:", l_f
            !     print *, "v_f is:", v_f
            !     print *, "error is:", ABS((v_f - 2.0d0*a*del_t - v_sub)/v_sub)
            if (ABS(l_f - l_sub) >= 1) then
                stop "l_f has crossed boundary when condition says is shouldn't have any substeps"
            end if
            if (boolDepositJ) solver%J(INT(l_cell)) = solver%J(INT(l_cell)) + part%w_p * part%q * (v_f + v_sub)*(del_t - timePassed)/2/world%dx_dl(INT(l_cell))/del_t
            timePassed = del_t
            
        else
            v_f = 2 * (l_f - l_sub) * world%dx_dl(INT(l_cell)) / del_tau - v_sub
            timePassed = timePassed + del_tau
            if (boolDepositJ) solver%J(INT(l_cell)) = solver%J(INT(l_cell)) + part%w_p * part%q * (v_f + v_sub)*del_tau/2/world%dx_dl(INT(l_cell))/del_t
            if (MOD(l_f, 1.0) /= 0.0) then
                print *, l_f
                stop "l_f is not integer after subStep"
            end if
            ! now final position/velocity becomes next starting position/velocity
            l_sub = l_f
            v_sub = v_f
    
        end if

    end subroutine particleSubStep

    subroutine depositJ(self, particleList, world, del_t, boolDepositJ, boolDiagnostic)
        ! particle substepping procedure which deposits J, also used for general moving of particles
        ! boolDepositJ = true : deposit J
        ! boolDepositJ = false: only move particles (done after non-linear iteration), also delete particles for wall collisions (saves extra loop in collisions)
        class(potSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        type(Particle), intent(in out) :: particleList(:)
        real(real64), intent(in) :: del_t
        logical, intent(in) :: boolDiagnostic, boolDepositJ
        logical :: boolDepositJTemp, delParticle
        !a and c correspond to quadratic equations | l_alongV is nearest integer boundary along velocity component, away is opposite
        real(real64) :: l_f, l_sub, v_sub, v_f, timePassed, del_tau, l_alongV, l_awayV, l_cell
        integer(int32) :: subStepNum, j, i, delIdx
        real(real64), allocatable :: rho_f(:), KE_i, KE_f, PE_i, PE_f
        boolDepositJTemp = boolDepositJ
        delParticle = .false.
        if (boolDiagnostic) then
            allocate(rho_f(n_x))
            call self%depositRho(particleList, world)
            rho_f = 0
            KE_i = 0
            KE_f = 0
            PE_i = self%getTotalPE(world, .false.)
            PE_f = 0
            boolDepositJTemp = .true.
        end if
        if (boolDepositJTemp) self%J = 0
        loopSpecies: do j = 1, numberChargedParticles
            delIdx = 0
            loopParticles: do i = 1, particleList(j)%N_p
                v_sub = particleList(j)%phaseSpace(2,i)
                l_sub = particleList(j)%phaseSpace(1,i)
                timePassed = 0
                subStepNum = 0
                if (boolDiagnostic) KE_i = KE_i + (v_sub**2) * particleList(j)%mass * 0.5d0 * particleList(j)%w_p
                do while((timePassed < del_t))
                    if (subStepNum == 0) then
                        ! Initial sub-step
                        call getl_BoundaryInitial(l_sub, v_sub, l_alongV, l_awayV)
                        call particleSubStepInitial(self, world, particleList(j), l_sub, l_f, v_sub, v_f, timePassed, del_tau, del_t, l_alongV, l_awayV, boolDepositJTemp)
                    else
                        ! Further sub-steps, particles start on grid nodes
                        l_alongV = l_sub + SIGN(1.0, v_sub)
                        call particleSubStep(self, world, particleList(j), l_sub, l_f, v_sub, v_f, timePassed, del_tau, del_t, l_alongV, l_cell, boolDepositJTemp)
                    end if
                    subStepNum = subStepNum + 1
                    if ((l_sub == 1.0) .or. (l_sub == n_x)) then
                        !check if particle has hit boundary, in which case delete
                        !will eventually need to replace with subroutine which takes in boundary inputs from world
                        delParticle = .true.
                        exit
                    end if
                end do
                if (boolDiagnostic) call singleRhoPass(rho_f, l_f, particleList(j)%w_p, particleList(j)%q, world%nodeVol) 
                if (boolDiagnostic) KE_f = KE_f + (v_f**2) * particleList(j)%mass * 0.5d0 * particleList(j)%w_p
                if ((l_f < 1) .or. (l_f > n_x)) then
                    stop "Have particles travelling oremainDel_tutside domain!"
                end if
                if (.not. boolDepositJ) then
                    ! When not depositing, then updating particles, overwrite deleted indices
                    if (delParticle) then
                        delIdx = delIdx + 1
                        delParticle = .false.
                        if (boolDiagnostic) then
                            self%particleEnergyLoss = self%particleEnergyLoss + particleList(j)%w_p * SUM(particleList(j)%phaseSpace(2:4, i)**2) * particleList(j)%mass * 0.5d0 !J/m^2 in 1D
                            self%particleChargeLoss = self%particleChargeLoss + particleList(j)%q * particleList(j)%w_p !C/m^2 in 1D
                        end if
                    else
                        particleList(j)%phaseSpace(1, i-delIdx) = l_f
                        particleList(j)%phaseSpace(2,i-delIdx) = v_f
                    end if
                end if
            end do loopParticles
            
            if (.not. boolDepositJ) then
                particleList(j)%phaseSpace(1, particleList(j)%N_p+1-delIdx:particleList(j)%N_p+1) = 0.0d0
                particleList(j)%phaseSpace(2:4,particleList(j)%N_p+1-delIdx:particleList(j)%N_p+1) = 0.0d0
                particleList(j)%N_p = particleList(j)%N_p - delIdx
            end if
        end do loopSpecies
        if (boolDiagnostic) PE_f = self%getTotalPE(world, .true.)
        ! Check final charge conservation
        if (boolDiagnostic) then
            self%chargeError = 0.0d0
            j = 0
            do i = 1, size(self%J) -1
                if ((self%J(i + 1) - self%J(i) /= 0) .and. (self%rho(i+1) /= 0)) then
                    self%chargeError = self%chargeError + (1 + (self%J(i + 1) - self%J(i)) *del_t/ world%nodeVol(i+1)/(rho_f(i+1) - self%rho(i+1)))**2
                    j = j + 1
                end if
            end do

            self%chargeError = SQRT(self%chargeError/j)
            if (self%chargeError > 1e-6) then
                print *, "-------------------------WARNING------------------------"
                print *, "Charge error is:", self%chargeError
                print *, "Charge error array is:"
                print *, ABS(1 + (self%J(2:n_x-1) - self%J(1:n_x-2)) *del_t/ world%nodeVol(2:n_x-1)/(rho_f(2:n_x-1) - self%rho(2:n_x-1)))
                stop "Total charge not conserved over time step in sub-step procedure!"
            end if

            self%energyError = ABS((KE_i + PE_i - KE_f - PE_f)/(KE_i + PE_i))
            if (self%energyError > 1e-8) then
                print *, "-------------------------WARNING------------------------"
                print *, "Energy error is:", self%energyError
                stop "Total energy not conserved over time step in sub-step procedure!"
            end if
            deallocate(rho_f)
        end if
    end subroutine depositJ

    subroutine moveParticles(self, particleList, world, del_t, boolDiagnostic)
        ! particle substepping procedure which deposits J, also used for general moving of particles
        ! boolDepositJ = true : deposit J
        ! boolDepositJ = false: only move particles (done after non-linear iteration), also delete particles for wall collisions (saves extra loop in collisions)
        class(potSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        type(Particle), intent(in out) :: particleList(:)
        real(real64), intent(in) :: del_t
        logical, intent(in) :: boolDiagnostic
        logical :: delParticle
        !a and c correspond to quadratic equations | l_alongV is nearest integer boundary along velocity component, away is opposite
        real(real64) :: l_f, l_sub, v_sub, v_f, timePassed, del_tau, l_alongV, l_awayV, l_cell, KE_i, KE_f, PE_i, PE_f
        integer(int32) :: subStepNum, j, i, delIdx
        real(real64), allocatable :: rho_f(:)
        delParticle = .false.
        if (boolDiagnostic) then
            allocate(rho_f(n_x))
            call self%depositRho(particleList, world)
            rho_f = 0
            KE_i = 0
            KE_f = 0
            PE_i = self%getTotalPE(world, .false.)
            PE_f = 0
            self%J = 0
        end if
        loopSpecies: do j = 1, numberChargedParticles
            delIdx = 0
            loopParticles: do i = 1, particleList(j)%N_p
                v_sub = particleList(j)%phaseSpace(2,i)
                l_sub = particleList(j)%phaseSpace(1,i)
                timePassed = 0
                subStepNum = 0
                if (boolDiagnostic) KE_i = KE_i + (v_sub**2) * particleList(j)%mass * 0.5d0 * particleList(j)%w_p
                do while((timePassed < del_t))
                    if (subStepNum == 0) then
                        ! Initial sub-step
                        call getl_BoundaryInitial(l_sub, v_sub, l_alongV, l_awayV)
                        call particleSubStepInitial(self, world, particleList(j), l_sub, l_f, v_sub, v_f, timePassed, del_tau, del_t, l_alongV, l_awayV, boolDiagnostic)
                    else
                        ! Further sub-steps, particles start on grid nodes
                        l_alongV = l_sub + SIGN(1.0, v_sub)
                        call particleSubStep(self, world, particleList(j), l_sub, l_f, v_sub, v_f, timePassed, del_tau, del_t, l_alongV, l_cell, boolDiagnostic)
                    end if
                    subStepNum = subStepNum + 1
                    if ((l_sub == 1.0) .or. (l_sub == n_x)) then
                        !check if particle has hit boundary, in which case delete
                        !will eventually need to replace with subroutine which takes in boundary inputs from world
                        delParticle = .true.
                        exit
                    end if
                end do
                if (boolDiagnostic) then 
                    call singleRhoPass(rho_f, l_f, particleList(j)%w_p, particleList(j)%q, world%nodeVol) 
                    KE_f = KE_f + (v_f**2) * particleList(j)%mass * 0.5d0 * particleList(j)%w_p
                end if 
                if ((l_f < 1) .or. (l_f > n_x)) then
                    stop "Have particles travelling oremainDel_tutside domain!"
                end if
                
                ! When not depositing, then updating particles, overwrite deleted indices
                if (delParticle) then
                    delIdx = delIdx + 1
                    delParticle = .false.
                    if (boolDiagnostic) then
                        self%particleEnergyLoss = self%particleEnergyLoss + particleList(j)%w_p * SUM(particleList(j)%phaseSpace(2:4, i)**2) * particleList(j)%mass * 0.5d0 !J/m^2 in 1D
                        self%particleChargeLoss = self%particleChargeLoss + particleList(j)%q * particleList(j)%w_p !C/m^2 in 1D
                    end if
                else
                    particleList(j)%phaseSpace(1, i-delIdx) = l_f
                    particleList(j)%phaseSpace(2,i-delIdx) = v_f
                end if
                
            end do loopParticles
            

            particleList(j)%phaseSpace(1, particleList(j)%N_p+1-delIdx:particleList(j)%N_p+1) = 0.0d0
            particleList(j)%phaseSpace(2:4,particleList(j)%N_p+1-delIdx:particleList(j)%N_p+1) = 0.0d0
            particleList(j)%N_p = particleList(j)%N_p - delIdx
            
        end do loopSpecies
        ! Check final charge conservation
        if (boolDiagnostic) then
            PE_f = self%getTotalPE(world, .true.)
            self%chargeError = 0.0d0
            j = 0
            do i = 1, size(self%J) -1
                if ((self%J(i + 1) - self%J(i) /= 0) .and. (self%rho(i+1) /= 0)) then
                    self%chargeError = self%chargeError + (1 + (self%J(i + 1) - self%J(i)) *del_t/ world%nodeVol(i+1)/(rho_f(i+1) - self%rho(i+1)))**2
                    j = j + 1
                end if
            end do

            self%chargeError = SQRT(self%chargeError/j)
            if (self%chargeError > 1e-6) then
                print *, "-------------------------WARNING------------------------"
                print *, "Charge error is:", self%chargeError
                print *, "Charge error array is:"
                print *, ABS(1 + (self%J(2:n_x-1) - self%J(1:n_x-2)) *del_t/ world%nodeVol(2:n_x-1)/(rho_f(2:n_x-1) - self%rho(2:n_x-1)))
                stop "Total charge not conserved over time step in sub-step procedure!"
            end if

            self%energyError = ABS((KE_i + PE_i - KE_f - PE_f)/(KE_i + PE_i))
            if (self%energyError > 1e-8) then
                print *, "-------------------------WARNING------------------------"
                print *, "Energy error is:", self%energyError
                stop "Total energy not conserved over time step in sub-step procedure!"
            end if
            deallocate(rho_f)
        end if
    end subroutine moveParticles

    ! ---------------- Initial Poisson Solver -------------------------------------------------

    subroutine solveInitialPotential(self, particleList, world)
        ! Solve for initial potential
        class(potSolver), intent(in out) :: self
        type(Particle), intent(in) :: particleList(:)
        type(Domain), intent(in) :: world
        call self%depositRho(particleList, world)
        call self%solve_tridiag_Poisson()
        ! Assume only use potential solver once, then need to generate matrix for Div-Ampere
        call self%construct_diagMatrix_Ampere(world)

    end subroutine solveInitialPotential

    !--------------------------- non-linear solver for time step using divergence of ampere ------------------------

    subroutine solveDivAmperePicard(self, particleList, world, del_t, maxIter, eps_a, boolDiagnostic)
        ! Solve for divergence of ampere using picard iterations
        class(potSolver), intent(in out) :: self
        type(Particle), intent(in out) :: particleList(:)
        type(Domain), intent(in) :: world
        integer(int32), intent(in) :: maxIter
        real(real64), intent(in) :: del_t, eps_a
        logical, intent(in) :: boolDiagnostic
        real(real64) :: errorCurrent, errorInitial
        integer(int32) :: i
        call self%depositJ(particleList, world, del_t, .true., .false.)
        errorInitial = self%getError_tridiag_Ampere(world, del_t)
        do i = 1, maxIter
            call self%solve_tridiag_Ampere(world, del_t)
            call self%depositJ(particleList, world, del_t, .true., .false.)
            errorCurrent = self%getError_tridiag_Ampere(world, del_t)
            if (i > 2) then
                if (errorCurrent < eps_a*errorInitial) then
                    call self%moveParticles(particleList, world, del_t, boolDiagnostic)
                    self%phi = self%phi_f
                    exit
                end if
            end if
        end do
        self%iterNumPicard = i-1
        

    end subroutine solveDivAmperePicard

    subroutine adaptiveSolveDivAmperePicard(self, particleList, world, del_t, maxIter, eps_a, boolDiagnostic)
        ! Solve for divergence of ampere's law with picard
        ! cut del_t in two if non-convergence after maxIter, repeat until convergence
        class(potSolver), intent(in out) :: self
        type(Particle), intent(in out) :: particleList(:)
        type(Domain), intent(in) :: world
        integer(int32), intent(in) :: maxIter
        real(real64), intent(in) :: del_t, eps_a
        logical, intent(in) :: boolDiagnostic
        real(real64) :: errorInitial, remainDel_t, currDel_t
        call self%solveDivAmperePicard(particleList, world, del_t, maxIter, eps_a, boolDiagnostic)
        remainDel_t = del_t  
        do while (self%iterNumPicard == maxIter)
            self%iterNumAdaptiveSteps = 0
            currDel_t = remainDel_t
            do while (self%iterNumPicard == maxIter)
                print *, "Reached maximum picard iteration count, reducing step size."
                currDel_t = currDel_t/2.0d0
                self%iterNumAdaptiveSteps = self%iterNumAdaptiveSteps + 1
                if (self%iterNumAdaptiveSteps > 2) then
                    stop "ALREADY REDUCED TIME STEP MORE THAN 2 TIMES, REDUCE INITIAL TIME STEP!!!"
                end if
                call self%solveDivAmperePicard(particleList, world, currDel_t, maxIter, eps_a, boolDiagnostic)   
            end do
            remainDel_t = remainDel_t - currDel_t   
            call self%solveDivAmperePicard(particleList, world, remainDel_t, maxIter, eps_a, boolDiagnostic)
        end do
       
        

    end subroutine adaptiveSolveDivAmperePicard



end module mod_potentialSolver