module mod_potentialSolver
    use iso_fortran_env, only: int32, real64
    use constants
    use mod_particle
    use mod_BasicFunctions
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
        real(real64) :: phi_left, phi_right
        real(real64) :: coeff_left, coeff_right ! these are coefficients (from world dimensions) needed with phi_left and phi_right in rhs of matrix equation
        real(real64), allocatable :: a_tri(:), b_tri(:), c_tri(:) !for thomas algorithm potential solver, a_tri is lower diagonal, b_tri middle, c_tri upper

    contains
        procedure, public, pass(self) :: depositRho
        procedure, public, pass(self) :: solve_tridiag_Poisson
        procedure, public, pass(self) :: getEField
        procedure, public, pass(self) :: depositJ
        procedure, private, pass(self) :: construct_diagMatrix
    end type

    interface potSolver
        module procedure :: potSolver_constructor
    end interface potSolver

contains

    type(potSolver) function potSolver_constructor(num_nodes, world) result(self)
        ! Construct domain object, initialize grid, dx_dl, and nodeVol.
        integer(int32), intent(in) :: num_nodes
        type(Domain), intent(in) :: world
        allocate(self % J(num_nodes-1), self % rho(num_nodes), self % phi(num_nodes), self % phi_f(num_nodes), self%a_tri(num_nodes-3), &
        self%b_tri(num_nodes-2), self%c_tri(num_nodes-3))
        call construct_diagMatrix(self, world)
        self % rho = 0
        self % J = 0
        self % phi = 0
        self % phi_left = 0
        self % phi_right = 0
        self % coeff_left = 2/(world%dx_dl(1) + world%dx_dl(2))/world%dx_dl(1)
        self % coeff_right = 2/(world%dx_dl(size(world%dx_dl)-1) + world%dx_dl(size(world%dx_dl)))/world%dx_dl(size(world%dx_dl))

    end function potSolver_constructor

    subroutine depositRho(self, particleList, world) 
        class(potSolver), intent(in out) :: self
        type(Particle), intent(in) :: particleList(:)
        type(Domain), intent(in) :: world
        integer(int32) :: i, j, l_left
        real(real64) :: d
        self % rho = 0.0d0
        do i=1, size(particleList)
            do j = 1, particleList(i)%N_p
                l_left = INT(particleList(i)%l_p(j))
                d = MOD(particleList(i)%l_p(j), 1.0d0)
                self % rho(l_left) = self % rho(l_left) + particleList(i)%q * particleList(i)%w_p * (1.0d0-d)
                self % rho(l_left + 1) = self % rho(l_left + 1) + particleList(i)%q * particleList(i)%w_p * d
            end do
        end do
        self % rho = self % rho / world%nodeVol
    end subroutine depositRho

    pure function singleRho(n, l_p, w_p, q, nodeVol) result(rho)
        ! for diagnostic in substep routine
        integer(int32), intent(in) :: n
        real(real64), intent(in) :: l_p, w_p, q, nodeVol(:)
        real(real64) :: rho(n), d
        integer(int32) :: l_left
        rho = 0.0d0
        l_left = INT(l_p)
        d = MOD(l_p, 1.0d0)
        rho(l_left) = q * w_p * (1.0d0-d) / nodeVol(l_left)
        rho(l_left + 1) =  q * w_p * d / nodeVol(l_left+1)
    end function singleRho

    pure function singleGradJ(dx_dl, l_half, v_half, w_p, q, nodeVol) result(gradJ)
        ! for diagnostic in substep routine
        real(real64), intent(in) :: dx_dl(:), l_half, v_half, w_p, q, nodeVol(:)
        real(real64) :: gradJ(size(dx_dl)-1), J(size(dx_dl))
        integer(int32) :: i
        J = 0.0d0
        gradJ = 0.0d0
        J(INT(l_half)) = J(INT(l_half)) + w_p * q * (v_half)/dx_dl(INT(l_half))
        do i = 1, size(J) -1
            gradJ(i) = (J(i + 1) - J(i)) / nodeVol(i+1)
        end do
    end function singleGradJ

    subroutine construct_diagMatrix(self, world)
        ! construct diagonal components for thomas algorithm
        class(potSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        self % a_tri = 2/(world%dx_dl(2:world%n_x-2) + world%dx_dl(3:)) / world%dx_dl(2:world%n_x-2)
        self % c_tri = 2/(world%dx_dl(1:world%n_x-3) + world%dx_dl(2:world%n_x-2))/world%dx_dl(2:world%n_x-2)
        self % b_tri = -2/(world%dx_dl(1:world%n_x-2) + world%dx_dl(2:))/world%dx_dl(1:world%n_x-2) - 2/(world%dx_dl(1:world%n_x-2) + world%dx_dl(2:))/world%dx_dl(2:world%n_x-1)

    end subroutine construct_diagMatrix

    subroutine solve_tridiag_Poisson(self)
        ! Tridiagonal (Thomas algorithm) solver for initial Poisson
        class(potSolver), intent(in out) :: self
        integer(int32) :: i, n !n size dependent on how many points are boundary conditions and thus moved to rhs equation
        real(real64) :: m, d(size(self%phi) - 2), cp(size(self%phi) - 2),dp(size(self%phi) - 2)
        n = size(self%phi) - 2

        d = -self%rho / eps_0
        d(1) = d(1) - 2 * self%phi_left * self%coeff_left
        d(n) = d(n) - 2 * self%phi_right * self%coeff_right
    ! initialize c-prime and d-prime
        cp(1) = self%c_tri(1)/self%b_tri(1)
        dp(1) = d(1)/self%b_tri(1)
    ! solve for vectors c-prime and d-prime
        do i = 2,n
            m = self%b_tri(i)-cp(i-1)*self%a_tri(i-1)
            cp(i) = self%c_tri(i)/m
            dp(i) = (d(i)-dp(i-1)*self%a_tri(i-1))/m
        end do
    ! initialize x
        self%phi(2:n+1) = dp(n)
    ! solve for x from the vectors c-prime and d-prime
        do i = n-1, 1, -1
            self%phi(i+1) = dp(i)-cp(i)*self%phi(i+2)
        end do
        self%phi_f = self%phi

    end subroutine solve_tridiag_Poisson

    !----------------- Particle mover/sub-stepping procedures -------------------------------------------------------

    pure function getEField(self, l_p, world) result(EField)
        !return EField of particle at logical position l_p (this is half distance), per particle since particle mover loop will be per particle
        class(potSolver), intent(in) :: self
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: l_p
        real(real64) :: EField
        if (l_p < 1 .or. l_p >= world%n_x) then
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


    subroutine particleSubStepInitial(solver, world, part, l_sub, l_f, v_sub, v_f, timePassed, del_tau, del_t, l_alongV, l_awayV)
        ! Do initial substep, where particles start between nodes
        type(potSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        type(Particle), intent(in) :: part
        real(real64), intent(in out) :: l_sub, l_f, v_sub, v_f, timePassed, del_tau, l_alongV, l_awayV
        real(real64), intent(in) :: del_t
        real(real64) :: rho_i(world%n_x), rho_f(world%n_x), gradJ(world%n_x-2), test(world%n_x-2), testConserv, a, c
        integer(int32) :: k
        
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
            timePassed = timePassed + del_t
            if (ABS((v_f - 2*a*del_t - v_sub)/v_sub) > 1e-8) then
                stop "Kinematic equation doesn't match for no substep"
            else if (INT(l_f) /= INT(l_sub)) then
                stop "l_f has crossed boundary when condition says is shouldn't have any substeps"
            end if
            solver%J(INT(l_sub)) = solver%J(INT(l_sub)) + part%w_p * part%q * (v_f + v_sub)/2/world%dx_dl(INT(l_sub))
            rho_i = singleRho(world%n_x, l_sub, part%w_p, part%q, world%nodeVol)
            rho_f = singleRho(world%n_x, l_f, part%w_p, part%q, world%nodeVol)
            gradJ = singleGradJ(world%dx_dl, l_sub, (v_f + v_sub)/2, part%w_p, part%q, world%nodeVol)
            testConserv = 0.0d0
            test = (rho_f(2:size(rho_f)-1) - rho_i(2:size(rho_f)-1)) + gradJ*del_t
            do k = 1, world%n_x - 2
                if (test(k) /= 0.0) then
                    testConserv = testConserv + (test(k)/(rho_f(k+1) - rho_i(k+1)))**2
                end if
            end do
            if (SQRT(testConserv) > 1e-8) then
                stop "Failed non-substep charge conservation, subStepNum = 0"
            end if
            
        else
            v_f = 2 * (l_f - l_sub) * world%dx_dl(INT(l_sub)) / del_tau - v_sub
            timePassed = timePassed + del_tau
            solver%J(INT(l_sub)) = solver%J(INT(l_sub)) + part%w_p * part%q * (v_f + v_sub)*del_tau/2/world%dx_dl(INT(l_sub))/del_t
            rho_i = singleRho(world%n_x, l_sub, part%w_p, part%q, world%nodeVol)
            rho_f = singleRho(world%n_x, l_f, part%w_p, part%q, world%nodeVol)
            gradJ = singleGradJ(world%dx_dl, l_sub, (v_f + v_sub)*del_tau/2/del_t, part%w_p, part%q, world%nodeVol)
            testConserv = 0.0d0
            test = (rho_f(2:size(rho_f)-1) - rho_i(2:size(rho_f)-1)) + gradJ*del_t
            do k = 1, world%n_x - 2
                if (test(k) /= 0.0) then
                    testConserv = testConserv + (test(k)/(rho_f(k+1) - rho_i(k+1)))**2
                end if
            end do
            if (SQRT(testConserv) > 1e-8) then
                stop "Failed substep charge conservation, subStepNum = 0"
            end if
            if (MOD(l_f, 1.0) /= 0.0) then
                print *, l_f
                stop "l_f is not integer after subStep"
            end if
            ! now final position/velocity becomes next starting position/velocity
            l_sub = l_f
            v_sub = v_f
    
        end if

    end subroutine particleSubStepInitial

    subroutine depositJ(self, particleList, world, del_t)
        ! particle substepping procedure which deposits J
        class(potSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        type(Particle), intent(in) :: particleList(:)
        real(real64), intent(in) :: del_t
        !a and c correspond to quadratic equations | l_alongV is nearest integer boundary along velocity component, away is opposite
        real(real64) :: l_f, l_sub, v_sub, v_f, timePassed, del_tau, l_alongV, l_awayV
        integer(int32) :: subStepNum, j, i
        self%J = 0
        loopSpecies: do j = 1, size(particleList)
            loopParticles: do i = 1, particleList(j)%N_p
                v_sub = particleList(j)%v_p(i, 1)
                l_sub = particleList(j)%l_p(i)
                timePassed = 0
                subStepNum = 0
                del_tau = 0
                print *, i
                do while((timePassed < del_t) .and. (l_sub /= 1.0) .and. (l_sub /= world%n_x) .and. (subStepNum < 2))
                    if (subStepNum == 0) then
                        ! Initial sub-step
                        print *, "l_sub is:", l_sub
                        print *, "v_sub is:", v_sub
                        call getl_BoundaryInitial(l_sub, v_sub, l_alongV, l_awayV)
                        call particleSubStepInitial(self, world, particleList(j), l_sub, l_f, v_sub, v_f, timePassed, del_tau, del_t, l_alongV, l_awayV)
                        print *, "l_f is now:", l_f
                        print *, "v_f is now::", v_f
                    else
                        ! Further sub-steps, particles start on grid nodes
                        l_alongV = l_sub + SIGN(1.0, v_sub)
                        print *, "l_sub is:", l_sub
                        print *, "l_alongV is:", l_alongV
                    end if
                    subStepNum = subStepNum + 1
                end do

            end do loopParticles
        end do loopSpecies

    end subroutine depositJ





end module mod_potentialSolver