module mod_particleMover
    use iso_fortran_env, only: int32, real64, output_unit
    use constants
    use mod_BasicFunctions
    use mod_particle
    use mod_domain
    use mod_potentialSolver
    implicit none

contains
    
    function getEField(solver, l_p, l_cell) result(EField)
        !return EField of particle at logical position l_p (this is half distance), per particle since particle mover loop will be per particle
        type(potentialSolver), intent(in) :: solver
        real(real64), intent(in) :: l_p
        integer(int32), intent(in) :: l_cell
        real(real64) :: EField, d
        d = l_p - real(l_cell, kind = real64) + 0.5d0
        EField = ((solver%phi_f(l_cell) + solver%phi(l_cell) - solver%phi(l_cell+1) - solver%phi_f(l_cell + 1)) * d +  &
        (solver%phi_f(l_cell-1) + solver%phi(l_cell-1) - solver%phi(l_cell) - solver%phi_f(l_cell)) * (1.0d0 - d)) * 0.5d0
    end function getEField

    function getEFieldNeumann(solver, l_p, l_cell) result(EField)
        !return EField of particle at logical position l_p (this is half distance), per particle since particle mover loop will be per particle
        type(potentialSolver), intent(in) :: solver
        real(real64), intent(in) :: l_p
        integer(int32), intent(in) :: l_cell
        real(real64) :: EField, d
        integer(int32) :: int_l_p
        int_l_p = INT(l_p)
        if (int_l_p == NumberXNodes) int_l_p = int_l_p - 1
        d = l_p - real(l_cell, kind = real64)
        EField = (solver%phi_f(int_l_p) + solver%phi(int_l_p) - solver%phi(int_l_p + 1) - solver%phi_f(int_l_p + 1)) * ABS(d)
    end function getEFieldNeumann

    pure function getEFieldPeriodic(solver, l_p, l_cell) result(EField)
        !return EField of particle at logical position l_p (this is half distance), per particle since particle mover loop will be per particle
        type(potentialSolver), intent(in) :: solver
        real(real64), intent(in) :: l_p
        integer(int32), intent(in) :: l_cell
        real(real64) :: EField, d
        d = l_p - real(l_cell, kind = real64) + 0.5d0
        EField = 0.5d0*((solver%phi_f(1) + solver%phi(1) - solver%phi(2) - solver%phi_f(2)) * d +  &
        (solver%phi_f(NumberXNodes-1) + solver%phi(NumberXNodes-1) - solver%phi(1) - solver%phi_f(1)) * (1.0d0 - d))
    end function getEFieldPeriodic

    pure function getEFieldDirichlet(solver, int_l_sub) result(EField)
        !return EField of particle at logical position l_p (this is half distance), per particle since particle mover loop will be per particle
        type(potentialSolver), intent(in) :: solver
        integer(int32), intent(in) :: int_l_sub
        real(real64) :: EField
        EField = (solver%phi_f(int_l_sub) + solver%phi(int_l_sub) - solver%phi(int_l_sub+ 1) - solver%phi_f(int_l_sub + 1)) * 0.5d0
    end function getEFieldDirichlet

    function getDx(l_sub, l_f, world, int_l_sub) result(res)
        ! Keep in mind, int_l_sub is the half-node integer for field/current density
        ! int_l_sub is evaluated at beginning of substep
        ! int_l_f may change during calculations based on calculations of different l_alongV, l_awayV, picard, etc, so will have change done here
        ! int_l_f only needs to be changed from INT(l_f) when l_f == NumberXNodes
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: l_sub, l_f
        integer(int32), intent(in) :: int_l_sub
        integer(int32) :: int_l_f
        real(real64) :: res, x_f, x_i
        int_l_f = INT(l_f)
        if (int_l_f == NumberXNodes) int_l_f = int_l_f - 1
        if (int_l_f == int_l_sub) then
            res = world%dx_dl(int_l_sub) * (l_f - l_sub)
        else
            x_i = world%grid(int_l_sub) + (l_sub - real(int_l_sub, kind = real64)) * world%dx_dl(int_l_sub)
            x_f = world%grid(int_l_f) + (l_f - real(int_l_f, kind = real64)) * world%dx_dl(int_l_f)
            res = (x_f - x_i)
        end if
    end function getDx

    subroutine getDelTauAlongVelocity(delX, a, l_alongV, l_sub, l_f, v_sub, l_cell, del_tau)
        real(real64), intent(in) :: l_sub, v_sub, l_alongV, a, delX
        real(real64), intent(in out) :: del_tau, l_f
        integer(int32), intent(in) :: l_cell
        if (a*v_sub > 0) then
            ! velocity and acceleration in same direction
            del_tau = 2.0d0 * ABS(delX)/(ABS(v_sub) + SQRT(v_sub**2 + 4.0d0*a*delX))
            l_f = l_alongV
            if (del_tau <= 0.0d0) then
                stop "Have issue with del_tau for v,a in same direction"
            end if
        else if (v_sub**2 + 4.0d0*a*delX > 0.0d0) then
            ! v and a opposite direction, but particle can still reach boundary along v
            del_tau = 2.0d0 * ABS(delX)/(ABS(v_sub) + SQRT(v_sub**2 + 4.0d0*a*delX))
            l_f = l_alongV
            if (del_tau <= 0.0d0) then
                print *, "l_sub is:", l_sub
                print *, "v_sub is:", v_sub
                print *, "a is:", a
                print *, "In initial sub-step routine"
                stop "Have issue with del_tau for v,a in opposite direction, but still reach boundary along v"
            end if
        end if
    end subroutine getDelTauAlongVelocity

    subroutine getDelTauInitialSubStep(solver, world, part, l_sub, l_f, v_sub, del_tau, del_t, l_alongV, l_awayV, int_l_sub, l_cell, a, c, fieldFunc)
    type(potentialSolver), intent(in) :: solver
    type(Domain), intent(in) :: world
    type(Particle), intent(in) :: part
    real(real64), intent(in out) :: l_sub, l_f, v_sub, del_tau, l_alongV, l_awayV, a, c
    real(real64), intent(in) :: del_t
    integer(int32), intent(in) :: l_cell, int_l_sub
    interface
        function fieldFunc(solver, l_p, l_cell) result(res)
            use iso_fortran_env, only: int32, real64
            use mod_potentialSolver
            type(potentialSolver), intent(in) :: solver
            real(real64), intent(in) :: l_p
            integer(int32), intent(in) :: l_cell
            real(real64) :: res
        end function fieldFunc
    end interface
    !real(real64), external :: fieldFunc
    real(real64) :: del_tau_tmp, delX
    integer(int32) :: int_l_alongV, int_l_awayV
    del_tau = del_t
    delX = getDx(l_sub, l_alongV, world, int_l_sub)
    a = (part%q / part%mass / 2.0d0) * fieldFunc(solver,(l_sub + l_alongV)*0.5d0, l_cell) * (l_alongV - l_sub) / delX
    ! Particle first between nodes, so solve quadratic for that particle depending on conditions
    if (v_sub/=0.0d0) then ! make first case, since pretty much always likely to be the case (could just not have, assume always field, never have to check)
        call getDelTauAlongVelocity(delX, a, l_alongV, l_sub, l_f, v_sub, l_cell, del_tau)
        if (del_tau >= del_t) then
            ! boundary opposite direction of v
            delX = getDx(l_sub, l_awayV, world, int_l_sub)
            a = (part%q / part%mass / 2.0d0) * fieldFunc(solver,(l_sub + l_awayV)*0.5d0, l_cell) * (l_awayV - l_sub) / delX
            if (a*v_sub < 0) then
                del_tau_tmp = 2.0d0 * ABS(delX)/(SQRT(v_sub**2 + 4.0d0*a*delX) - ABS(v_sub))
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
        l_alongV = real(l_cell, kind = real64) + SIGN(0.5d0, a)
        delX = getDX(l_sub, l_alongV, world, int_l_sub)
        del_tau = SQRT(delX/a)
        l_f = l_alongV
        if (del_tau <= 0.0d0) then
            stop "Have issue with del_tau for v = 0"
        end if
    end if

    end subroutine getDelTauInitialSubstep


    subroutine getDelTauSubStep(solver, world, part, l_sub, l_f, v_sub, del_tau, timePassed, del_t, l_alongV, int_l_sub, l_cell, a, c, fieldFunc)
        type(potentialSolver), intent(in) :: solver
        type(Domain), intent(in) :: world
        type(Particle), intent(in) :: part
        real(real64), intent(in out) :: l_sub, l_f, v_sub, del_tau, l_alongV, a, c
        real(real64), intent(in) :: del_t, timePassed
        integer(int32), intent(in) :: l_cell, int_l_sub
        real(real64) :: del_tau_tmp, delX
        integer(int32) :: int_l_alongV
        interface
            function fieldFunc(solver, l_p, l_cell) result(res)
                use iso_fortran_env, only: int32, real64
                use mod_potentialSolver
                type(potentialSolver), intent(in) :: solver
                real(real64), intent(in) :: l_p
                integer(int32), intent(in) :: l_cell
                real(real64) :: res
            end function fieldFunc
        end interface
        !real(real64), external :: fieldFunc
        del_tau = del_t - timePassed
        delX = getDx(l_sub, l_alongV, world, int_l_sub)
        ! get index cell where field and dx_dl is evaluated
        a = 0.5d0*(part%q / part%mass) * fieldFunc(solver,(l_sub + l_alongV)*0.5d0, l_cell) * (l_alongV - l_sub) / delX
        ! Particle first between nodes, so solve quadratic for that particle depending on conditions
        if ((a/=0.0d0)) then ! make first case, since pretty much always likely to be the case (could just not have, assume always field exists, never have to check)
            call getDelTauAlongVelocity(delX, a, l_alongV, l_sub, l_f, v_sub, l_cell, del_tau)
            if (del_tau >= del_t-timePassed) then
                ! For particles returning to boundary, make guess for maximum distance could go before v_f = 0 in half sub-step
                if (0.25d0*ABS((del_t - timePassed) *v_sub/delX) < 1.0) then
                    delX = 0.25d0*v_sub * (del_t - timePassed)
                    if (abs(delX) > 0.5 * world%dx_dl(int_l_sub)) then
                        int_l_alongV = INT(l_alongV)
                        if (int_l_alongV == NumberXNodes) int_l_alongV = int_l_alongV - 1
                        l_f = real(l_cell, kind = real64) + (delX - SIGN(0.5d0, v_sub) * world%dx_dl(int_l_sub))/world%dx_dl(int_l_alongV)
                        a = 0.5d0*(part%q / part%mass) * fieldFunc(solver,(l_sub + l_f)*0.5d0, l_cell) * (l_f - l_sub)/delX
                    else
                        l_f = l_sub + delX/world%dx_dl(int_l_sub)
                        if (INT(l_f) /= int_l_sub) then
                            print *, "In ongoing substep, when finding l_f for maximum distance for particle turn-around"
                            print *, "have issue with INT(l_f) == INT(l_sub)"
                            print *, "l_sub is:", l_sub
                            print *, "l_f is:", l_f
                            stop
                        end if
                        a = 0.5d0*(part%q / part%mass) * fieldFunc(solver,(l_sub + l_f) * 0.5d0, l_cell)/world%dx_dl(int_l_sub)
                    end if
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
                    a = 0.5d0*(part%q / part%mass) * fieldFunc(solver,l_sub, l_cell)/world%dx_dl(int_l_sub)
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
            del_tau = getDx(l_sub, l_alongV, world, int_l_sub)/v_sub
            l_f = l_alongV
            if (del_tau <= 0.0d0) then
                stop "Have issue with del_tau for a = 0"
            end if
        end if

    end subroutine getDelTauSubstep

    subroutine getDelTauInitialSubStepDirichlet(solver, world, part, l_sub, l_f, v_sub, del_tau, l_awayV, l_alongV, int_l_sub, a, c)
        type(potentialSolver), intent(in) :: solver
        type(Domain), intent(in) :: world
        type(Particle), intent(in) :: part
        real(real64), intent(in out) :: l_sub, l_f, v_sub, del_tau, l_awayV, l_alongV, a, c
        integer(int32), intent(in) :: int_l_sub
        a = (part%q / part%mass / 2.0d0) * getEFieldDirichlet(solver,int_l_sub)/world%dx_dl(int_l_sub)
        ! Particle first between nodes, so solve quadratic for that particle depending on conditions
        if ((v_sub/=0.0d0)) then ! make first case, since pretty much always likely to be the case (could just not have, assume always field, never have to check)
            c = (l_sub - l_alongV) * world%dx_dl(int_l_sub)
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
                c = (l_sub - l_awayV) * world%dx_dl(int_l_sub)
                del_tau = 2.0d0 * ABS(c)/(SQRT(v_sub**2 - 4.0d0*a*c) - ABS(v_sub))
                l_f = l_awayV
                if (del_tau <= 0) then
                    stop "Have issue with del_tau for v,a in opposite direction, boundary opposite v"
                end if
            end if
        else
            ! only in direction of field: USE l_alongV AS BOUNDARY ALONG DIRECTION OF a SINCE VELOCITY = 0!!!
            if (a > 0) then
                l_alongV = real(int_l_sub + 1, kind = real64)
            else
                l_alongV = real(int_l_sub, kind = real64)
            end if
            c = (l_sub - l_alongV) * world%dx_dl(int_l_sub)
            del_tau = SQRT(-c/a)
            l_f = l_alongV
            if (del_tau <= 0.0d0) then
                stop "Have issue with del_tau for v = 0"
            end if
        end if
    end subroutine getDelTauInitialSubStepDirichlet

    subroutine getDelTauSubStepDirichlet(solver, world, part, l_sub, l_f, v_sub, del_tau, l_alongV, int_l_sub, a, c)
        type(potentialSolver), intent(in) :: solver
        type(Domain), intent(in) :: world
        type(Particle), intent(in) :: part
        real(real64), intent(in out) :: l_sub, l_f, v_sub, del_tau, l_alongV, a, c
        integer(int32), intent(in) :: int_l_sub
        a = 0.5d0*(part%q / part%mass) * getEFieldDirichlet(solver,int_l_sub)/world%dx_dl(int_l_sub)
        ! Particle first between nodes, so solve quadratic for that particle depending on conditions
        if ((a/=0.0d0)) then ! make first case, since pretty much always likely to be the case (could just not have, assume always field exists, never have to check)
            c = (l_sub - l_alongV) * world%dx_dl(int_l_sub)
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
            del_tau = (l_alongV - l_sub) * world%dx_dl(int_l_sub)/v_sub
            l_f = l_alongV
            if (del_tau <= 0.0d0) then
                stop "Have issue with del_tau for a = 0"
            end if
        end if
    end subroutine getDelTauSubStepDirichlet

    subroutine picardIterParticles(solver, world, q, mass, l_sub, l_f, v_sub, v_f, int_l_sub, l_cell, del_t, timePassed, fieldFunc)
        type(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        real(real64), intent(in out) :: l_sub, l_f, v_sub, v_f
        real(real64), intent(in) :: del_t, timePassed, q, mass
        integer(int32), intent(in) :: l_cell, int_l_sub
        integer(int32) :: i
        real(real64) :: l_initialR, alpha, l_f_previous
        interface
            function fieldFunc(solver, l_p, l_cell) result(res)
                use iso_fortran_env, only: int32, real64
                use mod_potentialSolver
                type(potentialSolver), intent(in) :: solver
                real(real64), intent(in) :: l_p
                integer(int32), intent(in) :: l_cell
                real(real64) :: res
            end function fieldFunc
        end interface
        l_f_previous = l_sub
        alpha = 1.0d0/world%dx_dl(int_l_sub)
        v_f = v_sub + (q/mass) * fieldFunc(solver,(l_sub + l_f)*0.5d0, l_cell) * (del_t - timePassed) * alpha
        l_f = (v_f + v_sub) * (del_t - timePassed) * 0.5d0 * alpha + l_sub
        l_initialR = ABS(l_f - l_sub)
        do i = 1, 50
            l_f_previous = l_f
            alpha = (l_f - l_sub) / getDx(l_sub, l_f, world, int_l_sub)
            v_f = v_sub + (q/mass) * fieldFunc(solver,(l_sub + l_f)*0.5d0, l_cell) * (del_t - timePassed) * alpha
            l_f = (v_f + v_sub) * (del_t - timePassed) * 0.5d0 * alpha + l_sub
            if (ABS(l_f - l_f_previous) < 1.0d-8*(l_initialR + 1.0d0)) exit
        end do
        if (NINT(l_f) /= l_cell) then
            print *, "Have final l_f outside initial cell"
            print *, "particle charge is:", q
            print *, "a direction is:", (q/mass/2.0d0) * getEField(solver,(l_sub + l_f)/2.0d0, l_cell) * alpha
            print *, "del_t remaining is:", del_t - timePassed
            print *, "l_sub is:", l_sub
            print *, "l_f is:", l_f
            print *, "v_sub is:", v_sub
            print *, "v_f is:", v_f
            print *, "l_f analytical is:", l_f
            stop 
        end if
        if (i == 51) then
            stop "Picard not converged"
        end if
    end subroutine picardIterParticles

    subroutine depositJSubStep(solver, world, q, w_p, l_sub, l_f, v_f, v_sub, int_l_sub, l_cell, del_tau, del_t)
        type(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: q, w_p, v_f, v_sub, del_tau, del_t, l_sub, l_f
        integer(int32), intent(in) :: l_cell, int_l_sub
        real(real64) :: d
        d = (l_sub + l_f)/2.0d0 - REAL(l_cell, kind = real64) + 0.5d0
        solver%J(l_cell-1) = solver%J(l_cell-1) + (1.0d0 - d) * w_p * q * 0.5d0 * (v_f + v_sub)*del_tau * (l_f - l_sub)/getDx(l_sub, l_f, world, int_l_sub)/del_t
        solver%J(l_cell) = solver%J(l_cell) + d * w_p * q * 0.5d0 * (v_f + v_sub)*del_tau * (l_f - l_sub)/getDx(l_sub, l_f, world, int_l_sub)/del_t
    end subroutine depositJSubStep

    subroutine depositJSubStepNeumann(solver, world, q, w_p, l_sub, l_f, v_f, v_sub, int_l_sub, l_cell, del_tau, del_t)
        type(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: q, w_p, v_f, v_sub, del_tau, del_t, l_sub, l_f
        integer(int32), intent(in) :: l_cell, int_l_sub
        real(real64) :: l_half
        l_half = 0.5d0*(l_sub + l_f)
        solver%J(int_l_sub) = solver%J(int_l_sub) + ABS(l_half - real(l_cell, kind = real64)) * w_p * q * (v_f + v_sub)*del_tau/world%dx_dl(int_l_sub)/del_t
    end subroutine depositJSubStepNeumann

    subroutine depositJSubStepPeriodic(solver, world, q, w_p, l_sub, l_f, v_f, v_sub, int_l_sub, l_cell, del_tau, del_t)
        type(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: q, w_p, v_f, v_sub, del_tau, del_t, l_sub, l_f
        integer(int32), intent(in) :: l_cell, int_l_sub
        real(real64) :: d
        d = (l_sub + l_f)*0.5d0 - REAL(l_cell, kind = real64) + 0.5d0
        solver%J(NumberXNodes-1) = solver%J(NumberXNodes-1) + (1.0d0 - d) * w_p * q * (v_f + v_sub)*0.5d0*del_tau/world%dx_dl(int_l_sub)/del_t
        solver%J(1) = solver%J(1) + d * w_p * q * (v_f + v_sub)*0.5d0*del_tau/world%dx_dl(int_l_sub)/del_t

    end subroutine depositJSubStepPeriodic

    subroutine depositJSubStepDirichlet(solver, world, q, w_p, int_l_sub, v_f, v_sub, del_tau, del_t)
        type(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: q, w_p, v_f, v_sub, del_tau, del_t
        integer(int32), intent(in) :: int_l_sub
        solver%J(int_l_sub) = solver%J(int_l_sub) + w_p * q * (v_f + v_sub)*0.5d0*del_tau/world%dx_dl(int_l_sub)/del_t
    end subroutine depositJSubStepDirichlet

    subroutine depositRhoSingle(rho, l_sub, q, w_p, world) 
        real(real64), intent(in out) :: rho(NumberXNodes)
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: l_sub, q, w_p
        integer(int32) :: l_center, l_left, l_right
        real(real64) :: d
        l_center = NINT(l_sub)
        d = l_sub - l_center
        l_left = l_center -1
        l_right = l_center + 1
        SELECT CASE (world%boundaryConditions(l_center))
        CASE(0)
            ! Inside domain
            rho(l_center) = rho(l_center) + q * w_p * (0.75 - d**2)
            if (l_right /= NumberXNodes) then
                rho(l_right) = rho(l_right) + q * w_p * 0.5d0 * (0.5d0 + d)**2
            else
                SELECT CASE (world%boundaryConditions(l_right))
                CASE(1)
                    rho(l_right) = rho(l_right) + q * w_p * 0.5d0 * (0.5d0 + d)**2
                CASE(2)
                    rho(l_right) = rho(l_right) + q * w_p * (0.5d0 + d)**2
                CASE(3)
                    rho(l_right) = rho(l_right) + q * w_p * 0.5d0 * (0.5d0 + d)**2
                    rho(1) = rho(1) + q * w_p * 0.5d0 * (0.5d0 + d)**2
                END SELECT
            end if
            if (l_left /= 1) then
                rho(l_left) = rho(l_left) + q * w_p * 0.5d0 * (0.5d0 - d)**2
            else
                SELECT CASE (world%boundaryConditions(l_left))
                CASE(1)
                    rho(l_left) = rho(l_left) + q * w_p * 0.5d0 * (0.5d0 - d)**2
                CASE(2)
                    rho(l_left) = rho(l_left) + q * w_p * (0.5d0 - d)**2
                CASE(3)
                    rho(l_left) = rho(l_left) + q * w_p * 0.5d0 * (0.5d0 - d)**2
                    rho(NumberXNodes) = rho(NumberXNodes) + q * w_p * 0.5d0 * (0.5d0 - d)**2
                END SELECT
            end if
        CASE(1)
            !Dirichlet
            rho(l_center) = rho(l_center) + q * w_p * (1.0d0-ABS(d))
            rho(l_center + INT(SIGN(1.0, d))) = rho(l_center + INT(SIGN(1.0, d))) + q * w_p * ABS(d)
        CASE(2)
            !Neumann symmetric
            rho(l_center) = rho(l_center) + 2.0d0 * q * w_p * (0.75 - d**2)
            if ((INT(l_sub) /= NumberXNodes)) then
                rho(l_center + INT(SIGN(1.0, d))) = rho(l_center + INT(SIGN(1.0, d))) + q * w_p * (0.25d0 + d**2)
            else
                rho(NumberXNodes-1) = rho(NumberXNodes-1) + q * w_p * (0.25d0)
            end if
        CASE(3)
            ! Periodic
            rho(l_center) = rho(l_center) + q * w_p * (0.75 - d**2)
            ! towards domain
            rho(l_center+INT(SIGN(1.0, d))) = rho(l_center+INT(SIGN(1.0, d))) + q * w_p * 0.5d0 * (0.5d0 + ABS(d))**2
            ! across periodic boundary
            rho(MODULO(l_center-2*INT(SIGN(1.0, d)),NumberXNodes)) = rho(MODULO(l_center-2*INT(SIGN(1.0, d)),NumberXNodes)) + q * w_p * 0.5d0 * (0.5d0 - ABS(d))**2
        END SELECT
    end subroutine depositRhoSingle



    subroutine depositJ(solver, particleList, world, del_t)
    ! particle substepping procedure which deposits J, also used for general moving of particles
    ! boolDepositJ = true : deposit J
    ! boolDepositJ = false: only move particles (done after non-linear iteration), also delete particles for wall collisions (saves extra loop in collisions)
    type(potentialSolver), intent(in out) :: solver
    type(Domain), intent(in) :: world
    type(Particle), intent(in out) :: particleList(:)
    real(real64), intent(in) :: del_t
    !a and c correspond to quadratic equations | l_alongV is nearest integer boundary along velocity component, away is opposite
    real(real64) :: l_f, l_sub, v_sub, v_f, timePassed, del_tau, l_alongV, l_awayV, a, c
    integer(int32) :: subStepNum, j, i, l_cell, int_l_sub
    solver%J = 0.0d0
    loopSpecies: do j = 1, numberChargedParticles
        loopParticles: do i = 1, particleList(j)%N_p
            v_sub = particleList(j)%phaseSpace(2,i)
            l_sub = particleList(j)%phaseSpace(1,i)
            v_f = v_sub
            l_f = l_sub
            timePassed = 0.0d0
            subStepNum = 0
            ! First substep
            l_cell = NINT(l_sub)
            int_l_sub = INT(l_sub)
            firstStep: SELECT CASE (world%boundaryConditions(l_cell))
            CASE(0)
                ! Within Domain
                l_alongV = real(l_cell, kind = real64) + SIGN(0.5d0, v_sub)
                l_awayV = real(l_cell, kind = real64) - SIGN(0.5d0, v_sub)
                call getDelTauInitialSubStep(solver,world, particleList(j), l_sub, l_f, v_sub, del_tau, del_t, l_alongV, l_awayV, int_l_sub, l_cell, a, c, getEField)
                if (del_tau >= del_t) then
                    call picardIterParticles(solver, world, particleList(j)%q, particleList(j)%mass, l_sub, l_f, v_sub, v_f, int_l_sub, l_cell, del_t, timePassed, getEField)
                    !call solver%picardIterParticles(world, particleList(j)%q, particleList(j)%mass, l_sub, l_f, v_sub, v_f, l_cell, del_t, timePassed)
                    call depositJSubStep(solver,world, particleList(j)%q, particleList(j)%w_p, l_sub, l_f, v_f, v_sub, int_l_sub, l_cell, del_t, del_t) 
                    if (NINT(l_f) /= l_cell) then
                        print *, "After initial domain substep, l_f is not in correct cell"
                        stop
                    end if
                    timePassed = del_t
                else
                    v_f = 2.0d0 * getDx(l_sub, l_f, world, int_l_sub) / del_tau - v_sub
                    call depositJSubStep(solver,world, particleList(j)%q, particleList(j)%w_p, l_sub, l_f, v_f, v_sub, int_l_sub, l_cell, del_tau, del_t)
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
                l_alongV = NINT(2.0d0*l_sub + SIGN(0.5d0, v_sub))*0.5d0
                l_awayV = NINT(2.0d0*l_sub - SIGN(0.5d0, v_sub))*0.5d0
                call getDelTauInitialSubStepDirichlet(solver,world, particleList(j), l_sub, l_f, v_sub, del_tau, l_awayV, l_alongV, int_l_sub, a, c)
                if (del_tau >= del_t-timePassed) then
                    ! Add directly to J with no substep
                    l_f = v_sub * del_t / world%dx_dl(int_l_sub) + (a/ world%dx_dl(int_l_sub)) * del_t**2 + l_sub
                    v_f = 2.0d0 * (l_f - l_sub) * world%dx_dl(int_l_sub) / del_t - v_sub
                    call depositJSubStepDirichlet(solver,world, particleList(j)%q, particleList(j)%w_p, int_l_sub, v_f, v_sub, del_t, del_t)
                    if (NINT(l_f) /= l_cell) then
                        print *, "After initial Dirichlet substep, l_f is not in correct cell"
                        stop
                    end if
                    timePassed = del_t  
                else
                    v_f = 2.0d0 * (l_f - l_sub) * world%dx_dl(int_l_sub) / del_tau - v_sub
                    call depositJSubStepDirichlet(solver,world, particleList(j)%q, particleList(j)%w_p, int_l_sub, v_f, v_sub, del_tau, del_t)
                    if (l_f == l_cell) then ! if particle is now on node, must be boundary, exit
                        timePassed = del_t
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
            CASE(2)
                !Neumann-symmetric
                l_alongV = NINT(2.0d0*l_sub + SIGN(0.5d0, v_sub))*0.5d0
                l_awayV = NINT(2.0d0*l_sub - SIGN(0.5d0, v_sub))*0.5d0
                call getDelTauInitialSubStep(solver,world, particleList(j), l_sub, l_f, v_sub, del_tau, del_t, l_alongV, l_awayV, int_l_sub, l_cell, a, c, getEFieldNeumann)
                if (del_tau >= del_t) then
                    call picardIterParticles(solver, world, particleList(j)%q, particleList(j)%mass, l_sub, l_f, v_sub, v_f, int_l_sub, l_cell, del_t, timePassed, getEFieldNeumann)
                    call depositJSubStepNeumann(solver, world, particleList(j)%q, particleList(j)%w_p, l_sub, l_f, v_f, v_sub, int_l_sub, l_cell, del_t, del_t)  
                    if (NINT(l_f) /= l_cell) then
                        print *, "After initial domain substep, l_f is not in correct cell"
                        stop
                    end if
                    timePassed = del_t  
                else
                    v_f = 2.0d0 * (l_f - l_sub) * world%dx_dl(int_l_sub) / del_tau - v_sub
                    call depositJSubStepNeumann(solver,world, particleList(j)%q, particleList(j)%w_p, l_sub, l_f, v_f, v_sub, int_l_sub, l_cell, del_tau, del_t)
                    if (INT(l_f) == l_cell) v_f = -v_f
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
            CASE(3)
                l_alongV = NINT(2.0d0*l_sub + SIGN(0.5d0, v_sub))*0.5d0
                l_awayV = NINT(2.0d0*l_sub - SIGN(0.5d0, v_sub))*0.5d0
                call getDelTauInitialSubStep(solver,world, particleList(j), l_sub, l_f, v_sub, del_tau, del_t, l_alongV, l_awayV, int_l_sub, l_cell, a, c, getEFieldPeriodic)
                if (del_tau >= del_t-timePassed) then
                    call picardIterParticles(solver, world, particleList(j)%q, particleList(j)%mass, l_sub, l_f, v_sub, v_f, int_l_sub, l_cell, del_t, timePassed, getEFieldPeriodic)   
                    call depositJSubStepPeriodic(solver,world, particleList(j)%q, particleList(j)%w_p, l_sub, l_f, v_f, v_sub, int_l_sub, l_cell, del_t, del_t) 
                    if (NINT(l_f) /= l_cell) then
                        print *, "After initial periodic substep, l_f is not in correct cell"
                        stop
                    end if
                    timePassed = del_t  
                else
                    v_f = 2.0d0 * (l_f - l_sub) * world%dx_dl(int_l_sub) / del_tau - v_sub
                    call depositJSubStepPeriodic(solver,world, particleList(j)%q, particleList(j)%w_p, l_sub, l_f, v_f, v_sub, int_l_sub, l_cell, del_tau, del_t)
                    if ((l_f /= l_alongV) .and. (l_f /= l_awayV)) then
                        print *, "l_sub is:", l_sub
                        print *, "l_f is:", l_f
                        stop "l_f is not correct boundary after initial periodic substep"
                    end if
                    if (INT(l_f) == l_cell) l_f = ABS(l_f - NumberXNodes) + 1
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
                int_l_sub = INT(l_sub)
                if (int_l_sub == NumberXNodes) int_l_sub = int_l_sub - 1
                subStep: SELECT CASE (world%boundaryConditions(l_cell))
                CASE(0)
                    l_alongV = l_sub + SIGN(1.0d0, v_sub)
                    call getDelTauSubStep(solver,world, particleList(j), l_sub, l_f, v_sub, del_tau, timePassed, del_t, l_alongV, int_l_sub, l_cell, a, c, getEField)
                    if (del_tau >= del_t-timePassed) then
                        ! Add directly to J with no substep
                        call picardIterParticles(solver, world, particleList(j)%q, particleList(j)%mass, l_sub, l_f, v_sub, v_f, int_l_sub, l_cell, del_t, timePassed, getEField)
                        call depositJSubStep(solver,world, particleList(j)%q, particleList(j)%w_p, l_sub, l_f, v_f, v_sub, int_l_sub, l_cell, del_t - timePassed, del_t)
                        if (NINT(l_f) /= l_cell) then
                            print *, "After ongoing substep, l_f is not in correct cell"
                        end if
                        timePassed = del_t  
                    else 
                        if (l_f /= l_sub) then
                            v_f = 2.0d0 * getDx(l_sub, l_f, world, int_l_sub) / del_tau - v_sub
                            call depositJSubStep(solver,world, particleList(j)%q, particleList(j)%w_p, l_sub, l_f, v_f, v_sub, int_l_sub, l_cell, del_tau, del_t) 
                        else
                            v_f = -v_sub
                        end if
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
                    call getDelTauSubStepDirichlet(solver,world, particleList(j), l_sub, l_f, v_sub, del_tau, l_alongV, int_l_sub, a, c)
                    if (del_tau >= del_t-timePassed) then
                        ! Add directly to J with no substep
                        l_f = v_sub * (del_t-timePassed) / world%dx_dl(int_l_sub) + (a/ world%dx_dl(int_l_sub)) * (del_t - timePassed)**2 + l_sub
                        v_f = 2.0d0 * (l_f - l_sub) * world%dx_dl(int_l_sub) / (del_t - timePassed) - v_sub
                        call depositJSubStepDirichlet(solver,world, particleList(j)%q, particleList(j)%w_p, int_l_sub, v_f, v_sub, del_t-timePassed, del_t)
                        if (NINT(l_f) /= l_cell) then
                            print *, "a is:", a
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
                        v_f = 2.0d0 * (l_f - l_sub) * world%dx_dl(int_l_sub) / del_tau - v_sub
                        call depositJSubStepDirichlet(solver,world, particleList(j)%q, particleList(j)%w_p, int_l_sub, v_f, v_sub, del_tau, del_t)
                        if (l_f == l_cell) then
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
                CASE(2)
                    !Neumann Symmetric
                    l_alongV = l_sub + SIGN(0.5d0, v_sub)
                    call getDelTauSubStep(solver,world, particleList(j), l_sub, l_f, v_sub, del_tau, timePassed, del_t, l_alongV, int_l_sub, l_cell, a, c, getEFieldNeumann)
                    if (del_tau >= del_t-timePassed) then
                        ! Add directly to J with no substep
                        call picardIterParticles(solver, world, particleList(j)%q, particleList(j)%mass, l_sub, l_f, v_sub, v_f, int_l_sub, l_cell, del_t, timePassed, getEFieldNeumann)
                        call depositJSubStepNeumann(solver,world, particleList(j)%q, particleList(j)%w_p, l_sub, l_f, v_f, v_sub, int_l_sub, l_cell, del_t - timePassed, del_t)
                        if (NINT(l_f) /= l_cell) then
                            print *, "After ongoing substep, l_f is not in correct cell"
                        end if
                        timePassed = del_t  
                    else
                        if (l_f /= l_sub) then
                            v_f = 2.0d0 * (l_f - l_sub) * world%dx_dl(int_l_sub) / del_tau - v_sub
                            call depositJSubStepNeumann(solver,world, particleList(j)%q, particleList(j)%w_p, l_sub, l_f, v_f, v_sub, int_l_sub, l_cell, del_tau, del_t)
                        else
                            v_f = -v_sub
                        end if
                        if (INT(l_f) == l_cell) v_f = -v_f
                        timePassed = timePassed + del_tau
                        if ((MOD(l_f, 0.5d0) /= 0.0d0) .or. (ABS(l_f - l_sub) > 1.0d0)) then
                            print *, l_f
                            stop "l_f is not half integer after subStep or it is too far away"
                        end if
                        ! now final position/velocity becomes next starting position/velocity
                        l_sub = l_f
                        v_sub = v_f
                    end if
                CASE(3)
                    l_alongV = l_sub + SIGN(0.5d0, v_sub)
                    call getDelTauSubStep(solver,world, particleList(j), l_sub, l_f, v_sub, del_tau, timePassed, del_t, l_alongV, int_l_sub, l_cell, a, c, getEFieldPeriodic)
                    if (del_tau >= del_t-timePassed) then
                        ! Add directly to J with no substep
                        call picardIterParticles(solver, world, particleList(j)%q, particleList(j)%mass, l_sub, l_f, v_sub, v_f, int_l_sub, l_cell, del_t, timePassed, getEFieldPeriodic)
                        call depositJSubStepPeriodic(solver,world, particleList(j)%q, particleList(j)%w_p, l_sub, l_f, v_f, v_sub, int_l_sub, l_cell, del_t - timePassed, del_t)
                        if (NINT(l_f) /= l_cell) then
                            print *, "After final periodic substep, l_f is not in correct cell"
                            stop
                        end if
                        timePassed = del_t  
                    else
                        if (l_f /= l_sub) then
                            v_f = 2.0d0 * (l_f - l_sub) * world%dx_dl(int_l_sub) / del_tau - v_sub
                            call depositJSubStepPeriodic(solver,world, particleList(j)%q, particleList(j)%w_p, l_sub, l_f, v_f, v_sub, int_l_sub, l_cell, del_tau, del_t) 
                        else
                            v_f = -v_sub
                        end if
                        if ((MOD(l_f, 0.5d0) /= 0.0d0) .or. (ABS(l_f - l_sub) > 1.0d0)) then
                            print *, l_f
                            stop "l_f is not half integer after ogoing periodic subStep or it is too far away"
                        end if
                        if (INT(l_f) == l_cell) l_f = ABS(l_f - NumberXNodes) + 1
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

    subroutine moveParticles(solver, particleList, world, del_t)
    ! particle mover to avoid the boolean checks which mostly don't happen when depositing J
    type(potentialSolver), intent(in out) :: solver
    type(Domain), intent(in) :: world
    type(Particle), intent(in out) :: particleList(:)
    real(real64), intent(in) :: del_t
    !a and c correspond to quadratic equations | l_alongV is nearest integer boundary along velocity component, away is opposite
    real(real64) :: l_f, l_sub, v_sub, v_f, timePassed, del_tau, l_alongV, l_awayV, a, c
    integer(int32) :: subStepNum, j, i, delIdx, l_cell, int_l_sub
    loopSpecies: do j = 1, numberChargedParticles
        delIdx = 0
        loopParticles: do i = 1, particleList(j)%N_p
            v_sub = particleList(j)%phaseSpace(2,i)
            l_sub = particleList(j)%phaseSpace(1,i)
            v_f = v_sub
            l_f = l_sub
            timePassed = 0.0d0
            subStepNum = 0

            l_cell = NINT(l_sub)
            int_l_sub = INT(l_sub)
            firstStep: SELECT CASE (world%boundaryConditions(l_cell))
            CASE(0)
                ! Within Domain
                l_alongV = real(l_cell, kind = real64) + SIGN(0.5d0, v_sub)
                l_awayV = real(l_cell, kind = real64) - SIGN(0.5d0, v_sub)
                call getDelTauInitialSubStep(solver,world, particleList(j), l_sub, l_f, v_sub, del_tau, del_t, l_alongV, l_awayV, int_l_sub, l_cell, a, c, getEField)
                if (del_tau >= del_t) then
                    call picardIterParticles(solver, world, particleList(j)%q, particleList(j)%mass, l_sub, l_f, v_sub, v_f, int_l_sub, l_cell, del_t, timePassed, getEField)
                    if (NINT(l_f) /= l_cell) then
                        print *, "After initial domain substep, l_f is not in correct cell"
                        stop
                    end if
                    timePassed = del_t  
                    particleList(j)%phaseSpace(1, i-delIdx) = l_f
                    particleList(j)%phaseSpace(2,i-delIdx) = v_f
                    particleList(j)%phaseSpace(3:4, i-delIdx) = particleList(j)%phaseSpace(3:4, i) 
                else
                    v_f = 2.0d0 * getDx(l_sub, l_f, world, int_l_sub) / del_tau - v_sub
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
                l_alongV = NINT(2.0d0*l_sub + SIGN(0.5d0, v_sub))*0.5d0
                l_awayV = NINT(2.0d0*l_sub - SIGN(0.5d0, v_sub))*0.5d0
                call getDelTauInitialSubStepDirichlet(solver,world, particleList(j), l_sub, l_f, v_sub, del_tau, l_awayV, l_alongV, int_l_sub, a, c)
                if (del_tau >= del_t-timePassed) then
                    ! Add directly to J with no substep
                    l_f = v_sub * del_t / world%dx_dl(int_l_sub) + (a/ world%dx_dl(int_l_sub)) * del_t**2 + l_sub
                    v_f = 2.0d0 * (l_f - l_sub) * world%dx_dl(int_l_sub) / del_t - v_sub
                    if (NINT(l_f) /= l_cell) then
                        print *, "After initial Dirichlet substep, l_f is not in correct cell"
                        stop
                    end if
                    timePassed = del_t 
                    particleList(j)%phaseSpace(1, i-delIdx) = l_f
                    particleList(j)%phaseSpace(2,i-delIdx) = v_f
                    particleList(j)%phaseSpace(3:4, i-delIdx) = particleList(j)%phaseSpace(3:4, i)  
                else
                    v_f = 2.0d0 * (l_f - l_sub) * world%dx_dl(int_l_sub) / del_tau - v_sub
                    if (l_f == l_cell) then ! if particle is now on node, must be boundary, exit
                        delIdx = delIdx + 1
                        timePassed = del_t
                        solver%particleEnergyLoss = solver%particleEnergyLoss + particleList(j)%w_p * (v_f**2 + SUM(particleList(j)%phaseSpace(3:4, i)**2)) * particleList(j)%mass * 0.5d0 !J/m^2 in 1D
                        if (l_cell == 1) then
                            solver%particleChargeLoss(1, j) = solver%particleChargeLoss(1, j) + particleList(j)%q * particleList(j)%w_p !C/m^2 in 1D
                        else
                            solver%particleChargeLoss(2, j) = solver%particleChargeLoss(2, j) + particleList(j)%q * particleList(j)%w_p !C/m^2 in 1D
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
            CASE(2)
                !Neumann-symmetric
                l_alongV = NINT(2.0d0*l_sub + SIGN(0.5d0, v_sub))*0.5d0
                l_awayV = NINT(2.0d0*l_sub - SIGN(0.5d0, v_sub))*0.5d0
                call getDelTauInitialSubStep(solver,world, particleList(j), l_sub, l_f, v_sub, del_tau, del_t, l_alongV, l_awayV, int_l_sub, l_cell, a, c, getEFieldNeumann)
                if (del_tau >= del_t) then
                    call picardIterParticles(solver, world, particleList(j)%q, particleList(j)%mass, l_sub, l_f, v_sub, v_f, int_l_sub, l_cell, del_t, timePassed, getEFieldNeumann)
                    if (NINT(l_f) /= l_cell) then
                        print *, "After initial domain substep, l_f is not in correct cell"
                        stop
                    end if
                    timePassed = del_t  
                    particleList(j)%phaseSpace(1, i-delIdx) = l_f
                    particleList(j)%phaseSpace(2,i-delIdx) = v_f
                    particleList(j)%phaseSpace(3:4, i-delIdx) = particleList(j)%phaseSpace(3:4, i) 
                else
                    v_f = 2.0d0 * (l_f - l_sub) * world%dx_dl(int_l_sub) / del_tau - v_sub
                    if (INT(l_f) == l_cell) v_f = -v_f
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
            CASE(3)
                l_alongV = NINT(2.0d0*l_sub + SIGN(0.5d0, v_sub))*0.5d0
                l_awayV = NINT(2.0d0*l_sub - SIGN(0.5d0, v_sub))*0.5d0
                call getDelTauInitialSubStep(solver,world, particleList(j), l_sub, l_f, v_sub, del_tau, del_t, l_alongV, l_awayV, int_l_sub, l_cell, a, c, getEFieldPeriodic)
                if (del_tau >= del_t-timePassed) then
                    call picardIterParticles(solver, world, particleList(j)%q, particleList(j)%mass, l_sub, l_f, v_sub, v_f, int_l_sub, l_cell, del_t, timePassed, getEFieldPeriodic)   
                    if (NINT(l_f) /= l_cell) then
                        print *, "After initial periodic substep, l_f is not in correct cell"
                        stop
                    end if
                    timePassed = del_t  
                    particleList(j)%phaseSpace(1, i-delIdx) = l_f
                    particleList(j)%phaseSpace(2,i-delIdx) = v_f
                    particleList(j)%phaseSpace(3:4, i-delIdx) = particleList(j)%phaseSpace(3:4, i)  
                else
                    v_f = 2.0d0 * (l_f - l_sub) * world%dx_dl(int_l_sub) / del_tau - v_sub
                    if ((l_f /= l_alongV) .and. (l_f /= l_awayV)) then
                        print *, "l_sub is:", l_sub
                        print *, "l_f is:", l_f
                        stop "l_f is not correct boundary after initial periodic substep"
                    end if
                    if (INT(l_f) == l_cell) l_f = ABS(l_f - NumberXNodes) + 1
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
                int_l_sub = INT(l_sub)
                if (int_l_sub == NumberXNodes) int_l_sub = int_l_sub - 1
                subStep: SELECT CASE (world%boundaryConditions(l_cell))
                CASE(0)
                    l_alongV = l_sub + SIGN(1.0d0, v_sub)
                    call getDelTauSubStep(solver,world, particleList(j), l_sub, l_f, v_sub, del_tau, timePassed, del_t, l_alongV, int_l_sub, l_cell, a, c, getEField)
                    if (del_tau >= del_t-timePassed) then
                        ! Add directly to J with no substep
                        call picardIterParticles(solver, world, particleList(j)%q, particleList(j)%mass, l_sub, l_f, v_sub, v_f, int_l_sub, l_cell, del_t, timePassed, getEField)
                        if (NINT(l_f) /= l_cell) then
                            print *, "After ongoing substep, l_f is not in correct cell"
                        end if
                        timePassed = del_t 
                        particleList(j)%phaseSpace(1, i-delIdx) = l_f
                        particleList(j)%phaseSpace(2,i-delIdx) = v_f
                        particleList(j)%phaseSpace(3:4, i-delIdx) = particleList(j)%phaseSpace(3:4, i)  
                    else
                        if (l_f /= l_sub) then
                            v_f = 2.0d0 * getDx(l_sub, l_f, world, int_l_sub) / del_tau - v_sub
                        else
                            v_f = -v_sub
                        end if
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
                    call getDelTauSubStepDirichlet(solver,world, particleList(j), l_sub, l_f, v_sub, del_tau, l_alongV, int_l_sub, a, c)
                    if (del_tau >= del_t-timePassed) then
                        ! Add directly to J with no substep
                        l_f = v_sub * (del_t-timePassed) / world%dx_dl(int_l_sub) + (a/ world%dx_dl(int_l_sub)) * (del_t - timePassed)**2 + l_sub
                        v_f = 2.0d0 * (l_f - l_sub) * world%dx_dl(int_l_sub) / (del_t - timePassed) - v_sub
                        if (NINT(l_f) /= l_cell) then
                            print *, "a is:", a
                            print *, 'v_sub is:', v_sub
                            print *, "remaining time:", del_t - timePassed
                            print *, "l_f is:", l_f
                            print *, "l_sub is:", l_sub
                            print *, "l_cell is:", l_cell
                            print *, "After ongoing Dirichlet last substep, l_f is not in correct cell"
                            stop
                        end if
                        timePassed = del_t 
                        particleList(j)%phaseSpace(1, i-delIdx) = l_f
                        particleList(j)%phaseSpace(2,i-delIdx) = v_f
                        particleList(j)%phaseSpace(3:4, i-delIdx) = particleList(j)%phaseSpace(3:4, i)   
                    else
                        v_f = 2.0d0 * (l_f - l_sub) * world%dx_dl(int_l_sub) / del_tau - v_sub
                        if (l_f == l_cell) then 
                            delIdx = delIdx + 1
                            solver%particleEnergyLoss = solver%particleEnergyLoss + particleList(j)%w_p * (v_f**2 + SUM(particleList(j)%phaseSpace(3:4, i)**2)) * particleList(j)%mass * 0.5d0 !J/m^2 in 1D
                            if (l_cell == 1) then
                                solver%particleChargeLoss(1, j) = solver%particleChargeLoss(1, j) + particleList(j)%q * particleList(j)%w_p !C/m^2 in 1D
                            else
                                solver%particleChargeLoss(2, j) = solver%particleChargeLoss(2, j) + particleList(j)%q * particleList(j)%w_p !C/m^2 in 1D
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
                CASE(2)
                    !Neumann Symmetric
                    l_alongV = l_sub + SIGN(0.5d0, v_sub)
                    call getDelTauSubStep(solver,world, particleList(j), l_sub, l_f, v_sub, del_tau, timePassed, del_t, l_alongV, int_l_sub, l_cell, a, c, getEFieldNeumann)
                    if (del_tau >= del_t-timePassed) then
                        ! Add directly to J with no substep
                        call picardIterParticles(solver, world, particleList(j)%q, particleList(j)%mass, l_sub, l_f, v_sub, v_f, int_l_sub, l_cell, del_t, timePassed, getEFieldNeumann)
                        if (NINT(l_f) /= l_cell) then
                            print *, "After ongoing substep, l_f is not in correct cell"
                        end if
                        timePassed = del_t  
                        particleList(j)%phaseSpace(1, i-delIdx) = l_f
                        particleList(j)%phaseSpace(2,i-delIdx) = v_f
                        particleList(j)%phaseSpace(3:4, i-delIdx) = particleList(j)%phaseSpace(3:4, i) 
                    else
                        if (l_f /= l_sub) then
                            v_f = 2.0d0 * (l_f - l_sub) * world%dx_dl(int_l_sub) / del_tau - v_sub
                        else
                            v_f = -v_sub
                        end if
                        if (INT(l_f) == l_cell) v_f = -v_f
                        timePassed = timePassed + del_tau
                        if ((MOD(l_f, 0.5d0) /= 0.0d0) .or. (ABS(l_f - l_sub) > 1.0d0)) then
                            print *, l_f
                            stop "l_f is not half integer after subStep or it is too far away"
                        end if
                        ! now final position/velocity becomes next starting position/velocity
                        l_sub = l_f
                        v_sub = v_f
                    end if
                CASE(3)
                    l_alongV = l_sub + SIGN(0.5d0, v_sub)
                    call getDelTauSubStep(solver,world, particleList(j), l_sub, l_f, v_sub, del_tau, timePassed, del_t, l_alongV, int_l_sub, l_cell, a, c, getEFieldPeriodic)
                    if (del_tau >= del_t-timePassed) then
                        ! Add directly to J with no substep
                        call picardIterParticles(solver, world, particleList(j)%q, particleList(j)%mass, l_sub, l_f, v_sub, v_f, int_l_sub, l_cell, del_t, timePassed, getEFieldPeriodic)
                        if (NINT(l_f) /= l_cell) then
                            print *, "After final periodic substep, l_f is not in correct cell"
                            stop
                        end if
                        timePassed = del_t
                        particleList(j)%phaseSpace(1, i-delIdx) = l_f
                        particleList(j)%phaseSpace(2,i-delIdx) = v_f
                        particleList(j)%phaseSpace(3:4, i-delIdx) = particleList(j)%phaseSpace(3:4, i)   
                    else
                        v_f = 2.0d0 * (l_f - l_sub) * world%dx_dl(int_l_sub) / del_tau - v_sub
                        if ((MOD(l_f, 0.5d0) /= 0.0d0) .or. (ABS(l_f - l_sub) > 1.0d0)) then
                            print *, l_f
                            stop "l_f is not half integer after ogoing periodic subStep or it is too far away"
                        end if
                        if (INT(l_f) == l_cell) l_f = ABS(l_f - NumberXNodes) + 1
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
        particleList(j)%N_p = particleList(j)%N_p - delIdx
    end do loopSpecies
    end subroutine moveParticles

end module mod_particleMover