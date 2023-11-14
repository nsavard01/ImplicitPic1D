module mod_particleMover
    use iso_fortran_env, only: int32, real64, output_unit
    use constants
    use mod_BasicFunctions
    use mod_particle
    use mod_domain
    use mod_potentialSolver
    use omp_lib
    implicit none


contains
    
    pure function getEField(solver, l_p, l_cell, world) result(EField)
        !return EField of particle at logical position l_p (this is half distance), per particle since particle mover loop will be per particle
        type(potentialSolver), intent(in) :: solver
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: l_p
        integer(int32), intent(in) :: l_cell
        real(real64) :: EField, d
        d = l_p - real(l_cell, kind = real64) + 0.5d0
        EField = (solver%EField(l_cell) * d + solver%EField(l_cell-1) * (1.0d0 - d)) / world%nodeVol(l_cell)
    end function getEField

    function getEFieldNeumann(solver, l_p, l_cell, world) result(EField)
        !return EField of particle at logical position l_p (this is half distance), per particle since particle mover loop will be per particle
        type(potentialSolver), intent(in) :: solver
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: l_p
        integer(int32), intent(in) :: l_cell
        real(real64) :: EField, d
        integer(int32) :: l_inner
        d = 2.0d0 * ABS(l_p - real(l_cell, kind = real64))
        l_inner = INT(l_p)
        EField = (solver%EField(l_inner) * d)/world%nodeVol(l_cell)
    end function getEFieldNeumann

    pure function getEFieldPeriodic(solver, l_p, l_cell, world) result(EField)
        !return EField of particle at logical position l_p (this is half distance), per particle since particle mover loop will be per particle
        type(potentialSolver), intent(in) :: solver
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: l_p
        integer(int32), intent(in) :: l_cell
        real(real64) :: EField, d
        d = l_p - real(l_cell, kind = real64) + 0.5d0
        EField = (solver%EField(1) * d +  solver%EField(NumberXNodes-1) * (1.0d0 - d))/world%nodeVol(l_cell)
    end function getEFieldPeriodic

    pure function getEFieldDirichlet(solver, l_p, l_cell, world) result(EField)
        !return EField of particle at logical position l_p (this is half distance), per particle since particle mover loop will be per particle
        type(potentialSolver), intent(in) :: solver
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: l_p
        integer(int32), intent(in) :: l_cell
        real(real64) :: EField
        integer(int32) :: l_inner
        l_inner = INT(l_p)
        EField = solver%EField(l_inner)/world%nodeVol(l_cell)
    end function getEFieldDirichlet

    subroutine getDelTauAlongVelocity(nodeVol, a, l_alongV, l_sub, l_f, v_sub, l_cell, del_tau)
        real(real64), intent(in) :: l_sub, v_sub, l_alongV, a, nodeVol(NumberXNodes)
        real(real64), intent(in out) :: del_tau, l_f
        integer(int32), intent(in) :: l_cell
        real(real64) :: c
        c = (l_sub - l_alongV) * nodeVol(l_cell)
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
    end subroutine getDelTauAlongVelocity

    subroutine getDelTauInitialSubStep(solver, world, part, l_sub, l_f, v_sub, del_tau, del_t, l_alongV, l_awayV, l_cell, a, c, fieldFunc)
    type(potentialSolver), intent(in) :: solver
    type(Domain), intent(in) :: world
    type(Particle), intent(in) :: part
    real(real64), intent(in out) :: l_sub, l_f, v_sub, del_tau, l_alongV, l_awayV, a, c
    real(real64), intent(in) :: del_t
    integer(int32), intent(in) :: l_cell
    interface
        function fieldFunc(solver, l_p, l_cell, world) result(res)
            use iso_fortran_env, only: int32, real64
            use mod_domain
            use mod_potentialSolver
            type(potentialSolver), intent(in) :: solver
            type(Domain), intent(in) :: world
            real(real64), intent(in) :: l_p
            integer(int32), intent(in) :: l_cell
            real(real64) :: res
        end function fieldFunc
    end interface
    !real(real64), external :: fieldFunc
    real(real64) :: del_tau_tmp
    del_tau = del_t
    a = 0.5d0 * (part%q / part%mass) * fieldFunc(solver,(l_sub + l_alongV)*0.5d0, l_cell, world)
    ! Particle first between nodes, so solve quadratic for that particle depending on conditions
    if (v_sub/=0.0d0) then ! make first case, since pretty much always likely to be the case (could just not have, assume always field, never have to check)
        call getDelTauAlongVelocity(world%nodeVol, a, l_alongV, l_sub, l_f, v_sub, l_cell, del_tau)
        if (del_tau >= del_t) then
            ! boundary opposite direction of v
            a = 0.5d0 * (part%q / part%mass) * fieldFunc(solver,(l_sub + l_awayV)*0.5d0, l_cell, world)
            if (a*v_sub < 0) then
                c = (l_sub - l_awayV) * world%nodeVol(l_cell)
                del_tau_tmp = 2.0d0 * ABS(c)/(SQRT(v_sub**2 - 4.0d0*a*c) - ABS(v_sub))
                if (del_tau_tmp < del_tau) then
                    ! If del_tau isn't reduced, then want to keep saved l_f since might make picard iteration more stable as initial condition
                    del_tau = del_tau_tmp
                    l_f = l_awayV
                end if
                if (del_tau <= 0) then
                    print *, "l_sub is:", l_sub
                    print *, 'l_f is:', l_f
                    print *, 'v_sub is:', v_sub
                    stop "Have issue with initial substep del_tau for v,a in opposite direction, boundary opposite v"
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
        c = (l_sub - l_alongV) * world%nodeVol(l_cell)
        del_tau = SQRT(-c/a)
        l_f = l_alongV
        if (del_tau <= 0.0d0) then
            stop "Have issue with del_tau for v = 0"
        end if
    end if

    end subroutine getDelTauInitialSubstep


    subroutine getDelTauSubStep(solver, world, part, l_sub, l_f, v_sub, del_tau, timePassed, del_t, l_alongV, l_cell, a, c, fieldFunc)
        type(potentialSolver), intent(in) :: solver
        type(Domain), intent(in) :: world
        type(Particle), intent(in) :: part
        real(real64), intent(in out) :: l_sub, l_f, v_sub, del_tau, l_alongV, a, c
        real(real64), intent(in) :: del_t, timePassed
        integer(int32), intent(in) :: l_cell
        real(real64) :: del_tau_tmp
        interface
            function fieldFunc(solver, l_p, l_cell, world) result(res)
                use iso_fortran_env, only: int32, real64
                use mod_domain
                use mod_potentialSolver
                type(potentialSolver), intent(in) :: solver
                type(Domain), intent(in) :: world
                real(real64), intent(in) :: l_p
                integer(int32), intent(in) :: l_cell
                real(real64) :: res
            end function fieldFunc
        end interface
        !real(real64), external :: fieldFunc
        del_tau = del_t - timePassed
        ! get index cell where field and dx_dl is evaluated
        a = 0.5d0 * (part%q / part%mass) * fieldFunc(solver,(l_sub + l_alongV)*0.5d0, l_cell, world)
        ! Particle first between nodes, so solve quadratic for that particle depending on conditions
        if ((a/=0.0d0)) then ! make first case, since pretty much always likely to be the case (could just not have, assume always field exists, never have to check)
            call getDelTauAlongVelocity(world%nodeVol, a, l_alongV, l_sub, l_f, v_sub, l_cell, del_tau)
            if (del_tau >= del_t-timePassed) then
                ! For particles returning to boundary, make guess for maximum distance could go before v_f = 0 in half sub-step
                if (ABS((del_t - timePassed) *v_sub/4.0d0/world%nodeVol(l_cell)) < 1.0) then
                    a = 0.5d0 * (part%q / part%mass) * fieldFunc(solver,l_sub + (del_t - timePassed) *v_sub/8.0d0/world%nodeVol(l_cell), l_cell, world)
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
                    a = (part%q / part%mass / 2.0d0) * fieldFunc(solver,l_sub + SIGN(1.0d-8, v_sub), l_cell, world)
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
            del_tau = (l_alongV - l_sub) * world%nodeVol(l_cell)/v_sub
            l_f = l_alongV
            if (del_tau <= 0.0d0) then
                stop "Have issue with del_tau for a = 0"
            end if
        end if

    end subroutine getDelTauSubstep

    subroutine getDelTauInitialSubStepDirichlet(solver, world, part, l_sub, l_f, v_sub, del_tau, l_awayV, l_alongV, l_cell, a, c)
        type(potentialSolver), intent(in) :: solver
        type(Domain), intent(in) :: world
        type(Particle), intent(in) :: part
        real(real64), intent(in out) :: l_sub, l_f, v_sub, del_tau, l_awayV, l_alongV, a, c
        integer(int32), intent(in) :: l_cell
        a = 0.5d0 * (part%q / part%mass) * getEFieldDirichlet(solver,l_sub, l_cell, world)
        ! Particle first between nodes, so solve quadratic for that particle depending on conditions
        if ((v_sub/=0.0d0)) then ! make first case, since pretty much always likely to be the case (could just not have, assume always field, never have to check)
            c = (l_sub - l_alongV) * world%nodeVol(l_cell)
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
                c = (l_sub - l_awayV) * world%nodeVol(l_cell)
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
            c = (l_sub - l_alongV) * world%nodeVol(l_cell)
            del_tau = SQRT(-c/a)
            l_f = l_alongV
            if (del_tau <= 0.0d0) then
                stop "Have issue with del_tau for v = 0"
            end if
        end if
    end subroutine getDelTauInitialSubStepDirichlet

    subroutine getDelTauSubStepDirichlet(solver, world, part, l_sub, l_f, v_sub, del_tau, l_alongV, l_cell, a, c)
        type(potentialSolver), intent(in) :: solver
        type(Domain), intent(in) :: world
        type(Particle), intent(in) :: part
        real(real64), intent(in out) :: l_sub, l_f, v_sub, del_tau, l_alongV, a, c
        integer(int32), intent(in) :: l_cell
        a = 0.5d0 * (part%q / part%mass) * getEFieldDirichlet(solver,l_sub, l_cell, world)
        ! Particle first between nodes, so solve quadratic for that particle depending on conditions
        if ((a/=0.0d0)) then ! make first case, since pretty much always likely to be the case (could just not have, assume always field exists, never have to check)
            c = (l_sub - l_alongV) * world%nodeVol(l_cell)
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
            del_tau = (l_alongV - l_sub) * world%nodeVol(l_cell)/v_sub
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
        x_i = world%grid(NINT(l_sub)) + (l_sub - NINT(l_sub)) * world%nodeVol(NINT(l_sub))
        x_f = world%grid(NINT(l_f)) + (l_f - NINT(l_f)) * world%nodeVol(NINT(l_f))
        res = (l_f - l_sub)/(x_f - x_i)
    else
        res = 1.0d0/world%nodeVol(NINT(l_sub))
    end if

    end function getAlpha

    subroutine analyticalParticleMover(solver, world, q, mass, l_sub, l_f, v_sub, v_f, l_cell, del_t, timePassed)
        type(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: l_sub, v_sub, del_t, timePassed, q, mass
        real(real64), intent(in out) :: l_f, v_f
        integer(int32), intent(in) :: l_cell
        real(real64) :: del_tau, del_tau_sqr, real_l_cell, dx
        dx = world%nodeVol(l_cell)
        del_tau = del_t - timePassed
        del_tau_sqr = del_tau**2
        real_l_cell = real(l_cell, kind = real64)
        l_f = (-2.0d0 * del_tau_sqr * real_l_cell*q*solver%EField(l_cell) + 2.0d0 * del_tau_sqr * real_l_cell*q*solver%EField(l_cell-1) + &
            del_tau_sqr * l_sub*q*solver%EField(l_cell) - del_tau_sqr * l_sub*q*solver%EField(l_cell-1) + &
            del_tau_sqr * q*solver%EField(l_cell) + del_tau_sqr * q*solver%EField(l_cell-1) + 4.0*del_tau*dx*mass*v_sub + &
            4.0 * dx**2 * l_sub*mass)/(-del_tau_sqr * q*solver%EField(l_cell) + &
            del_tau_sqr * q*solver%EField(l_cell-1) + 4.0 * dx**2 * mass)
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

    subroutine analyticalParticleMoverPeriodic(solver, world, q, mass, l_sub, l_f, v_sub, v_f, l_cell, del_t, timePassed)
        type(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: l_sub, v_sub, del_t, timePassed, q, mass
        real(real64), intent(in out) :: l_f, v_f
        integer(int32), intent(in) :: l_cell
        real(real64) :: del_tau, del_tau_sqr, real_l_cell, dx
        dx = world%nodeVol(l_cell)
        del_tau = del_t - timePassed
        del_tau_sqr = del_tau**2
        real_l_cell = real(l_cell, kind = real64)
        l_f = (-2.0d0 * del_tau_sqr * real_l_cell*q*solver%EField(1) + 2.0d0 * del_tau_sqr * real_l_cell*q*solver%EField(NumberXNodes-1) + &
            del_tau_sqr * l_sub*q*solver%EField(1) - del_tau_sqr * l_sub*q*solver%EField(NumberXNodes-1) + &
            del_tau_sqr * q*solver%EField(1) + del_tau_sqr * q*solver%EField(NumberXNodes-1) + 4.0*del_tau*dx*mass*v_sub + &
            4.0 * dx**2 * l_sub*mass)/(-del_tau_sqr * q*solver%EField(1) + &
            del_tau_sqr * q*solver%EField(NumberXNodes-1) + 4.0 * dx**2 * mass)
        v_f = 2.0d0 * (l_f - l_sub) * dx / del_tau - v_sub
        if (NINT(l_f) /= l_cell) then
            print *, "Have final l_f outside initial cell"
            print *, "particle charge is:", q
            stop 
        end if
    end subroutine analyticalParticleMoverPeriodic

    subroutine analyticalParticleMoverNeumann(solver, world, q, mass, l_sub, l_f, v_sub, v_f, l_cell, del_t, timePassed)
        type(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: l_sub, v_sub, del_t, timePassed, q, mass
        real(real64), intent(in out) :: l_f, v_f
        integer(int32), intent(in) :: l_cell
        real(real64) :: del_tau, del_tau_sqr, real_l_cell, dx
        integer(int32) :: l_inner
        dx = world%nodeVol(l_cell)
        del_tau = del_t - timePassed
        del_tau_sqr = del_tau**2
        real_l_cell = real(l_cell, kind = real64)
        if (l_cell == NumberXNodes) then
            l_inner = NumberXNodes-1
            l_f = (2.0d0*del_tau_sqr*real_l_cell*q*solver%EField(l_inner) - del_tau_sqr*l_sub*q*solver%EField(l_inner) + &
                2.0d0*del_tau*mass*v_sub*dx + 2.0d0*l_sub*mass* dx**2)/&
                (del_tau_sqr*q*solver%EField(l_inner) + 2.0d0*mass * dx**2)
        else
            l_inner = 1
            l_f = (2.0d0*del_tau_sqr*real_l_cell*q*solver%EField(l_inner) - del_tau_sqr*l_sub*q*solver%EField(l_inner) - &
                2.0d0*del_tau*mass*v_sub*dx - 2.0d0*l_sub*mass* dx**2)/&
                (del_tau_sqr*q*solver%EField(l_inner) - 2.0d0*mass * dx**2)
        end if
        v_f = 2.0d0 * (l_f - l_sub) * dx / del_tau - v_sub
        if (NINT(l_f) /= l_cell) then
            print *, "Have final l_f outside initial cell analytical mover neumann"
            print *, "particle charge is:", q
            print *, "l_sub is:", l_sub
            print *, "l_f is:", l_f
            print *, "v_sub is:", v_sub
            print *, "v_f is:", v_f
        end if
    end subroutine analyticalParticleMoverNeumann

    subroutine picardIterParticles(solver, world, q, mass, l_sub, l_f, v_sub, v_f, l_cell, del_t, timePassed)
    type(potentialSolver), intent(in out) :: solver
    type(Domain), intent(in) :: world
    real(real64), intent(in out) :: l_sub, l_f, v_sub, v_f
    real(real64), intent(in) :: del_t, timePassed, q, mass
    integer(int32), intent(in) :: l_cell
    integer(int32) :: i
    real(real64) :: l_f_previous
    l_f_previous = l_f
    v_f = v_sub + (q/mass) * getEField(solver,(l_sub + l_f)/2.0d0, l_cell, world) * (del_t - timePassed)
    l_f = (v_f + v_sub) * (del_t - timePassed) /world%nodeVol(l_cell) / 2.0d0 + l_sub
    do i = 1, 50
        l_f_previous = l_f
        v_f = v_sub + (q/mass) * getEField(solver,(l_sub + l_f)/2.0d0, l_cell, world) * (del_t - timePassed)
        l_f = (v_f + v_sub) * (del_t - timePassed) /world%nodeVol(l_cell) / 2.0d0 + l_sub
        if (NINT(l_f) /= l_cell) then
            l_f = l_cell + SIGN(0.5d0, l_f - l_cell)
        end if
        if (ABS(l_f - l_f_previous) < eps_r) exit
    end do
    if (NINT(l_f) /= l_cell) then
        print *, "Have final l_f outside initial cell"
        print *, "particle charge is:", q
        print *, "a direction is:", (q/mass/2.0d0) * getEField(solver,(l_sub + l_f)/2.0d0, l_cell, world)
        print *, "del_t remaining is:", del_t - timePassed
        print *, "l_sub is:", l_sub
        print *, "l_f is:", l_f
        print *, "v_sub is:", v_sub
        print *, "v_f is:", v_f
        v_f = v_sub + (q/mass) * getEField(solver,(l_sub + l_f)/2.0d0, l_cell, world) * (del_t - timePassed)
        l_f = (v_f + v_sub) * (del_t - timePassed) / world%nodeVol(l_cell) / 2.0d0 + l_sub
        print *, "next l_f is:", l_f
        print *, "next v_f is:", v_f
        v_f = v_sub + (q/mass) * getEField(solver,(l_sub + l_f)/2.0d0, l_cell, world) * (del_t - timePassed)
        l_f = (v_f + v_sub) * (del_t - timePassed) / world%nodeVol(l_cell) / 2.0d0 + l_sub
        print *, "next l_f is:", l_f
        print *, "next v_f is:", v_f
        call analyticalParticleMover(solver,world, q, mass, l_sub, l_f, v_sub, v_f, l_cell, del_t, timePassed)
        print *, "l_f analytical is:", l_f
        stop 
    end if
    if (i == 51) then
        stop "Picard not converged"
    end if
    end subroutine picardIterParticles

    subroutine depositJSubStep(solver, world, q, w_p, l_sub, l_f, v_f, v_sub, l_cell, del_tau, del_t, iThread)
        type(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: q, w_p, v_f, v_sub, del_tau, del_t, l_sub, l_f
        integer(int32), intent(in) :: l_cell, iThread
        real(real64) :: d
        d = (l_sub + l_f)/2.0d0 - REAL(l_cell, kind = real64) + 0.5d0
        solver%J(l_cell-1, iThread) = solver%J(l_cell-1, iThread) + 0.5d0 * (1.0d0 - d) * w_p * q * (v_f + v_sub)*del_tau/world%nodeVol(l_cell)/del_t
        solver%J(l_cell, iThread) = solver%J(l_cell, iThread) + 0.5d0 * d * w_p * q * (v_f + v_sub)*del_tau/world%nodeVol(l_cell)/del_t
    end subroutine depositJSubStep

    subroutine depositJSubStepNeumann(solver, world, q, w_p, l_sub, l_f, v_f, v_sub, l_cell, del_tau, del_t, iThread)
        type(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: q, w_p, v_f, v_sub, del_tau, del_t, l_sub, l_f
        integer(int32), intent(in) :: l_cell, iThread
        real(real64) :: l_half
        l_half = 0.5d0*(l_sub + l_f)
        solver%J(INT(l_half), iThread) = solver%J(INT(l_half), iThread) + ABS(l_half - real(l_cell, kind = real64)) * w_p * q * (v_f + v_sub)*del_tau/world%nodeVol(l_cell)/del_t
        ! solver%J(INT(l_half), iThread) = solver%J(INT(l_half), iThread) + 0.5d0 * (ABS(l_half - real(l_cell, kind = real64)) + 0.5d0) * w_p * q * (v_f + v_sub)*del_tau/world%nodeVol(l_cell)/del_t
    end subroutine depositJSubStepNeumann

    subroutine depositJSubStepPeriodic(solver, world, q, w_p, l_sub, l_f, v_f, v_sub, l_cell, del_tau, del_t, iThread)
        type(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: q, w_p, v_f, v_sub, del_tau, del_t, l_sub, l_f
        integer(int32), intent(in) :: l_cell, iThread
        real(real64) :: d
        d = (l_sub + l_f)/2.0d0 - REAL(l_cell, kind = real64) + 0.5d0
        solver%J(NumberXNodes-1, iThread) = solver%J(NumberXNodes-1, iThread) + 0.5d0 * (1.0d0 - d) * w_p * q * (v_f + v_sub)*del_tau/world%nodeVol(l_cell)/del_t
        solver%J(1, iThread) = solver%J(1, iThread) + 0.5d0 * d * w_p * q * (v_f + v_sub)*del_tau/world%nodeVol(l_cell)/del_t
    end subroutine depositJSubStepPeriodic

    subroutine depositJSubStepDirichlet(solver, world, q, w_p, l_sub, l_cell, v_f, v_sub, del_tau, del_t, iThread)
        type(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: q, w_p, v_f, v_sub, del_tau, del_t, l_sub
        integer(int32), intent(in) :: l_cell, iThread
        solver%J(INT(l_sub), iThread) = solver%J(INT(l_sub), iThread) + 0.5d0 * w_p * q * (v_f + v_sub)*del_tau/world%nodeVol(l_cell)/del_t
    end subroutine depositJSubStepDirichlet


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
    integer(int32) :: subStepNum, j, i, l_cell, iThread
    solver%J = 0.0d0
    solver%EField = 0.5d0 * (solver%phi(1:NumberXNodes-1) + solver%phi_f(1:NumberXNodes-1) - solver%phi(2:NumberXNodes) - solver%phi_f(2:NumberXNodes))
    loopSpecies: do j = 1, numberChargedParticles
        !$OMP parallel private(iThread, i, l_f, l_sub, v_sub, v_f, timePassed, del_tau, l_alongV, l_awayV, a, c,&
                    subStepNum, l_cell)
        iThread = omp_get_thread_num() + 1 
        loopParticles: do i = 1, particleList(j)%N_p(iThread)
            v_sub = particleList(j)%phaseSpace(2,i, iThread)
            l_sub = particleList(j)%phaseSpace(1,i, iThread)
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
                call getDelTauInitialSubStep(solver,world, particleList(j), l_sub, l_f, v_sub, del_tau, del_t, l_alongV, l_awayV, l_cell, a, c, getEField)
                if (del_tau >= del_t) then
                    call analyticalParticleMover(solver,world, particleList(j)%q, particleList(j)%mass, l_sub, l_f, v_sub, v_f, l_cell, del_t, timePassed)
                    !call solver%picardIterParticles(world, particleList(j)%q, particleList(j)%mass, l_sub, l_f, v_sub, v_f, l_cell, del_t, timePassed)
                    call depositJSubStep(solver,world, particleList(j)%q, particleList(j)%w_p, l_sub, l_f, v_f, v_sub, l_cell, del_t, del_t, iThread) 
                    if (NINT(l_f) /= NINT(l_sub)) then
                        print *, "After initial domain substep in regular domain, l_f is not in correct cell"
                        stop
                    end if
                    timePassed = del_t  
                else
                    v_f = 2.0d0 * (l_f - l_sub) * world%nodeVol(l_cell) / del_tau - v_sub
                    call depositJSubStep(solver,world, particleList(j)%q, particleList(j)%w_p, l_sub, l_f, v_f, v_sub, l_cell, del_tau, del_t, iThread)
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
                call getDelTauInitialSubStepDirichlet(solver,world, particleList(j), l_sub, l_f, v_sub, del_tau, l_awayV, l_alongV, l_cell, a, c)
                if (del_tau >= del_t-timePassed) then
                    ! Add directly to J with no substep
                    l_f = v_sub * del_t / world%nodeVol(l_cell) + (a/ world%nodeVol(l_cell)) * del_t**2 + l_sub
                    v_f = 2.0d0 * (l_f - l_sub) * world%nodeVol(l_cell) / del_t - v_sub
                    call depositJSubStepDirichlet(solver,world, particleList(j)%q, particleList(j)%w_p, l_sub, l_cell, v_f, v_sub, del_t, del_t, iThread)
                    if (NINT(l_f) /= NINT(l_sub)) then
                        print *, "After initial Dirichlet substep, l_f is not in correct cell"
                        stop
                    end if
                    timePassed = del_t  
                else
                    v_f = 2.0d0 * (l_f - l_sub) * world%nodeVol(l_cell) / del_tau - v_sub
                    call depositJSubStepDirichlet(solver,world, particleList(j)%q, particleList(j)%w_p, l_sub, l_cell, v_f, v_sub, del_tau, del_t, iThread)
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
            CASE(2)
                !Neumann-symmetric
                l_alongV = NINT(2.0d0*l_sub + SIGN(0.5d0, v_sub))/2.0d0
                l_awayV = NINT(2.0d0*l_sub - SIGN(0.5d0, v_sub))/2.0d0
                call getDelTauInitialSubStep(solver,world, particleList(j), l_sub, l_f, v_sub, del_tau, del_t, l_alongV, l_awayV, l_cell, a, c, getEFieldNeumann)
                if (del_tau >= del_t) then
                    ! print *, "-------------------------"
                    ! print *, "Entering Neumann final deposit:"
                    ! print *, 'l_sub is:', l_sub
                    ! print *, "v_sub is:", v_sub
                    call analyticalParticleMoverNeumann(solver,world, particleList(j)%q, particleList(j)%mass, l_sub, l_f, v_sub, v_f, l_cell, del_t, timePassed)
                    ! print *, "l_f is:", l_f
                    ! print *, "v_f is:", v_f
                    ! stop
                    !call solver%picardIterParticles(world, particleList(j)%q, particleList(j)%mass, l_sub, l_f, v_sub, v_f, l_cell, del_t, timePassed)
                    call depositJSubStepNeumann(solver,world, particleList(j)%q, particleList(j)%w_p, l_sub, l_f, v_f, v_sub, l_cell, del_t, del_t, iThread) 
                    if (NINT(l_f) /= NINT(l_sub)) then
                        print *, "After initial domain substep for neumann-symmetric, l_f is not in correct cell"
                        stop
                    end if
                    timePassed = del_t  
                else
                    v_f = 2.0d0 * (l_f - l_sub) * world%nodeVol(l_cell) / del_tau - v_sub
                    call depositJSubStepNeumann(solver,world, particleList(j)%q, particleList(j)%w_p, l_sub, l_f, v_f, v_sub, l_cell, del_tau, del_t, iThread)
                    if (l_f == l_cell) v_f = -v_f
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
                l_alongV = real(l_cell, kind = real64) + SIGN(0.5d0, v_sub)
                l_awayV = real(l_cell, kind = real64) - SIGN(0.5d0, v_sub)
                call getDelTauInitialSubStep(solver,world, particleList(j), l_sub, l_f, v_sub, del_tau, del_t, l_alongV, l_awayV, l_cell, a, c, getEFieldPeriodic)
                if (del_tau >= del_t-timePassed) then
                    call analyticalParticleMoverPeriodic(solver,world, particleList(j)%q, particleList(j)%mass, l_sub, l_f, v_sub, v_f, l_cell, del_t, timePassed)    
                    call depositJSubStepPeriodic(solver,world, particleList(j)%q, particleList(j)%w_p, l_sub, l_f, v_f, v_sub, l_cell, del_t, del_t, iThread) 
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
                    v_f = 2.0d0 * (l_f - l_sub) * world%nodeVol(l_cell) / del_tau - v_sub
                    call depositJSubStepPeriodic(solver,world, particleList(j)%q, particleList(j)%w_p, l_sub, l_f, v_f, v_sub, l_cell, del_tau, del_t, iThread)
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
                    call getDelTauSubStep(solver,world, particleList(j), l_sub, l_f, v_sub, del_tau, timePassed, del_t, l_alongV, l_cell, a, c, getEField)
                    if (del_tau >= del_t-timePassed) then
                        ! Add directly to J with no substep
                        call analyticalParticleMover(solver,world, particleList(j)%q, particleList(j)%mass, l_sub, l_f, v_sub, v_f, l_cell, del_t, timePassed)
                        !call solver%picardIterParticles(world, particleList(j)%q, particleList(j)%mass, l_sub, l_f, v_sub, v_f, l_cell, del_t, timePassed)
                        call depositJSubStep(solver,world, particleList(j)%q, particleList(j)%w_p, l_sub, l_f, v_f, v_sub, l_cell, del_t - timePassed, del_t, iThread)
                        if (NINT(l_f) /= l_cell) then
                            print *, "After ongoing substep, l_f is not in correct cell"
                        end if
                        timePassed = del_t  
                    else
                        v_f = 2.0d0 * (l_f - l_sub) * world%nodeVol(l_cell) / del_tau - v_sub
                        call depositJSubStep(solver,world, particleList(j)%q, particleList(j)%w_p, l_sub, l_f, v_f, v_sub, l_cell, del_tau, del_t, iThread) 
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
                    call getDelTauSubStepDirichlet(solver,world, particleList(j), l_sub, l_f, v_sub, del_tau, l_alongV, l_cell, a, c)
                    if (del_tau >= del_t-timePassed) then
                        ! Add directly to J with no substep
                        l_f = v_sub * (del_t-timePassed) / world%nodeVol(l_cell) + (a/ world%nodeVol(l_cell)) * (del_t - timePassed)**2 + l_sub
                        v_f = 2.0d0 * (l_f - l_sub) * world%nodeVol(l_cell) / (del_t - timePassed) - v_sub
                        call depositJSubStepDirichlet(solver,world, particleList(j)%q, particleList(j)%w_p, l_sub, l_cell, v_f, v_sub, del_t-timePassed, del_t, iThread)
                        if (NINT(l_f) /= l_cell) then
                            print *, "a is:", a
                            print *, world%nodeVol
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
                        v_f = 2.0d0 * (l_f - l_sub) * world%nodeVol(l_cell) / del_tau - v_sub
                        call depositJSubStepDirichlet(solver,world, particleList(j)%q, particleList(j)%w_p, l_sub, l_cell, v_f, v_sub, del_tau, del_t, iThread) 
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
                CASE(2)
                    !Neumann Symmetric
                    l_alongV = l_sub + SIGN(0.5d0, v_sub)
                    call getDelTauSubStep(solver,world, particleList(j), l_sub, l_f, v_sub, del_tau, timePassed, del_t, l_alongV, l_cell, a, c, getEFieldNeumann)
                    if (del_tau >= del_t-timePassed) then
                        ! Add directly to J with no substep
                        call analyticalParticleMoverNeumann(solver,world, particleList(j)%q, particleList(j)%mass, l_sub, l_f, v_sub, v_f, l_cell, del_t, timePassed)
                        !call solver%picardIterParticles(world, particleList(j)%q, particleList(j)%mass, l_sub, l_f, v_sub, v_f, l_cell, del_t, timePassed)
                        call depositJSubStepNeumann(solver,world, particleList(j)%q, particleList(j)%w_p, l_sub, l_f, v_f, v_sub, l_cell, del_t - timePassed, del_t, iThread)
                        if (NINT(l_f) /= l_cell) then
                            print *, "After ongoing substep, l_f is not in correct cell"
                        end if
                        timePassed = del_t  
                    else
                        v_f = 2.0d0 * (l_f - l_sub) * world%nodeVol(l_cell) / del_tau - v_sub
                        call depositJSubStepNeumann(solver,world, particleList(j)%q, particleList(j)%w_p, l_sub, l_f, v_f, v_sub, l_cell, del_tau, del_t, iThread) 
                        if (l_f == l_cell) v_f = -v_f
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
                    l_alongV = l_sub + SIGN(1.0d0, v_sub)
                    call getDelTauSubStep(solver,world, particleList(j), l_sub, l_f, v_sub, del_tau, timePassed, del_t, l_alongV, l_cell, a, c, getEFieldPeriodic)
                    if (del_tau >= del_t-timePassed) then
                        ! Add directly to J with no substep
                        call analyticalParticleMoverPeriodic(solver,world, particleList(j)%q, particleList(j)%mass, l_sub, l_f, v_sub, v_f, l_cell, del_t, timePassed)
                        call depositJSubStepPeriodic(solver,world, particleList(j)%q, particleList(j)%w_p, l_sub, l_f, v_f, v_sub, l_cell, del_t - timePassed, del_t, iThread)
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
                        v_f = 2.0d0 * (l_f - l_sub) * world%nodeVol(l_cell) / del_tau - v_sub
                        call depositJSubStepPeriodic(solver,world, particleList(j)%q, particleList(j)%w_p, l_sub, l_f, v_f, v_sub, l_cell, del_tau, del_t, iThread) 
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
        !$OMP end parallel
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
    integer(int32) :: subStepNum, j, i, l_cell, iThread, delIdx
    solver%EField = 0.5d0 * (solver%phi(1:NumberXNodes-1) + solver%phi_f(1:NumberXNodes-1) - solver%phi(2:NumberXNodes) - solver%phi_f(2:NumberXNodes))
    loopSpecies: do j = 1, numberChargedParticles
        !$OMP parallel private(iThread, i, l_f, l_sub, v_sub, v_f, timePassed, del_tau, l_alongV, l_awayV, a, c,&
                    subStepNum, delIdx, l_cell)
        iThread = omp_get_thread_num() + 1 
        delIdx = 0
        particleList(j)%refIdx(iThread) = 0
        particleList(j)%energyLoss(:, iThread) = 0.0d0
        particleList(j)%wallLoss(:, iThread) = 0.0d0
        loopParticles: do i = 1, particleList(j)%N_p(iThread)
            v_sub = particleList(j)%phaseSpace(2,i, iThread)
            l_sub = particleList(j)%phaseSpace(1,i, iThread)
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
                call getDelTauInitialSubStep(solver,world, particleList(j), l_sub, l_f, v_sub, del_tau, del_t, l_alongV, l_awayV, l_cell, a, c, getEField)
                if (del_tau >= del_t) then
                    call analyticalParticleMover(solver,world, particleList(j)%q, particleList(j)%mass, l_sub, l_f, v_sub, v_f, l_cell, del_t, timePassed)
                    !call solver%picardIterParticles(world, particleList(j)%q, particleList(j)%mass, l_sub, l_f, v_sub, v_f, l_cell, del_t, timePassed)
                    if (NINT(l_f) /= NINT(l_sub)) then
                        print *, "After initial domain substep, l_f is not in correct cell"
                        stop
                    end if
                    timePassed = del_t 
                    particleList(j)%phaseSpace(1, i-delIdx, iThread) = l_f
                    particleList(j)%phaseSpace(2,i-delIdx, iThread) = v_f
                    particleList(j)%phaseSpace(3:4, i-delIdx, iThread) = particleList(j)%phaseSpace(3:4, i, iThread) 
                else
                    v_f = 2.0d0 * (l_f - l_sub) * world%nodeVol(l_cell) / del_tau - v_sub
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
                call getDelTauInitialSubStepDirichlet(solver,world, particleList(j), l_sub, l_f, v_sub, del_tau, l_awayV, l_alongV, l_cell, a, c)
                if (del_tau >= del_t-timePassed) then
                    ! Add directly to J with no substep
                    l_f = v_sub * del_t / world%nodeVol(l_cell) + (a/ world%nodeVol(l_cell)) * del_t**2 + l_sub
                    v_f = 2.0d0 * (l_f - l_sub) * world%nodeVol(l_cell) / del_t - v_sub
                    if (NINT(l_f) /= NINT(l_sub)) then
                        print *, "After initial Dirichlet substep, l_f is not in correct cell"
                        stop
                    end if
                    timePassed = del_t 
                    particleList(j)%phaseSpace(1, i-delIdx, iThread) = l_f
                    particleList(j)%phaseSpace(2,i-delIdx, iThread) = v_f
                    particleList(j)%phaseSpace(3:4, i-delIdx, iThread) = particleList(j)%phaseSpace(3:4, i, iThread) 
                else
                    v_f = 2.0d0 * (l_f - l_sub) * world%nodeVol(l_cell) / del_tau - v_sub
                    if (l_f == l_cell) then ! if particle is now on node, must be boundary, exit
                        delIdx = delIdx + 1
                        timePassed = del_t
                        if (l_f == 1) then
                            particleList(j)%energyLoss(1, iThread) = particleList(j)%energyLoss(1, iThread) + (v_f**2 + SUM(particleList(j)%phaseSpace(3:4, i, iThread)**2))!J/m^2 in 1D
                            particleList(j)%wallLoss(1, iThread) = particleList(j)%wallLoss(1, iThread) + 1 !C/m^2 in 1D
                        else if (l_f == NumberXNodes) then
                            particleList(j)%energyLoss(2, iThread) = particleList(j)%energyLoss(2, iThread) + (v_f**2 + SUM(particleList(j)%phaseSpace(3:4, i, iThread)**2)) !J/m^2 in 1D
                            particleList(j)%wallLoss(2, iThread) = particleList(j)%wallLoss(2, iThread) + 1 !C/m^2 in 1D
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
                l_alongV = NINT(2.0d0*l_sub + SIGN(0.5d0, v_sub))/2.0d0
                l_awayV = NINT(2.0d0*l_sub - SIGN(0.5d0, v_sub))/2.0d0
                call getDelTauInitialSubStep(solver,world, particleList(j), l_sub, l_f, v_sub, del_tau, del_t, l_alongV, l_awayV, l_cell, a, c, getEFieldNeumann)
                if (del_tau >= del_t) then
                    call analyticalParticleMoverNeumann(solver,world, particleList(j)%q, particleList(j)%mass, l_sub, l_f, v_sub, v_f, l_cell, del_t, timePassed)
                    !call solver%picardIterParticles(world, particleList(j)%q, particleList(j)%mass, l_sub, l_f, v_sub, v_f, l_cell, del_t, timePassed)
                    if (NINT(l_f) /= NINT(l_sub)) then
                        print *, "After initial domain substep, l_f is not in correct cell"
                        stop
                    end if
                    timePassed = del_t
                    particleList(j)%phaseSpace(1, i-delIdx, iThread) = l_f
                    particleList(j)%phaseSpace(2,i-delIdx, iThread) = v_f
                    particleList(j)%phaseSpace(3:4, i-delIdx, iThread) = particleList(j)%phaseSpace(3:4, i, iThread)  
                else
                    v_f = 2.0d0 * (l_f - l_sub) * world%nodeVol(l_cell) / del_tau - v_sub
                    if (l_f == l_cell) then
                        v_f = -v_f
                        particleList(j)%refIdx(iThread) = particleList(j)%refIdx(iThread) + 1
                        particleList(j)%refRecordIdx(particleList(j)%refIdx(iThread), iThread) = i - delIdx
                    end if
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
                l_alongV = real(l_cell, kind = real64) + SIGN(0.5d0, v_sub)
                l_awayV = real(l_cell, kind = real64) - SIGN(0.5d0, v_sub)
                call getDelTauInitialSubStep(solver,world, particleList(j), l_sub, l_f, v_sub, del_tau, del_t, l_alongV, l_awayV, l_cell, a, c, getEFieldPeriodic)
                if (del_tau >= del_t-timePassed) then
                    call analyticalParticleMoverPeriodic(solver,world, particleList(j)%q, particleList(j)%mass, l_sub, l_f, v_sub, v_f, l_cell, del_t, timePassed)    
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
                    particleList(j)%phaseSpace(1, i-delIdx, iThread) = l_f
                    particleList(j)%phaseSpace(2,i-delIdx, iThread) = v_f
                    particleList(j)%phaseSpace(3:4, i-delIdx, iThread) = particleList(j)%phaseSpace(3:4, i, iThread)  
                else
                    v_f = 2.0d0 * (l_f - l_sub) * world%nodeVol(l_cell) / del_tau - v_sub
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
                    call getDelTauSubStep(solver,world, particleList(j), l_sub, l_f, v_sub, del_tau, timePassed, del_t, l_alongV, l_cell, a, c, getEField)
                    if (del_tau >= del_t-timePassed) then
                        ! Add directly to J with no substep
                        call analyticalParticleMover(solver,world, particleList(j)%q, particleList(j)%mass, l_sub, l_f, v_sub, v_f, l_cell, del_t, timePassed)
                        !call solver%picardIterParticles(world, particleList(j)%q, particleList(j)%mass, l_sub, l_f, v_sub, v_f, l_cell, del_t, timePassed)
                        if (NINT(l_f) /= l_cell) then
                            print *, "After ongoing substep, l_f is not in correct cell"
                        end if
                        timePassed = del_t 
                        particleList(j)%phaseSpace(1, i-delIdx, iThread) = l_f
                        particleList(j)%phaseSpace(2,i-delIdx, iThread) = v_f
                        particleList(j)%phaseSpace(3:4, i-delIdx, iThread) = particleList(j)%phaseSpace(3:4, i, iThread) 
                    else
                        v_f = 2.0d0 * (l_f - l_sub) * world%nodeVol(l_cell) / del_tau - v_sub
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
                    call getDelTauSubStepDirichlet(solver,world, particleList(j), l_sub, l_f, v_sub, del_tau, l_alongV, l_cell, a, c)
                    if (del_tau >= del_t-timePassed) then
                        ! Add directly to J with no substep
                        l_f = v_sub * (del_t-timePassed) / world%nodeVol(l_cell) + (a/ world%nodeVol(l_cell)) * (del_t - timePassed)**2 + l_sub
                        v_f = 2.0d0 * (l_f - l_sub) * world%nodeVol(l_cell) / (del_t - timePassed) - v_sub
                        if (NINT(l_f) /= l_cell) then
                            print *, "After ongoing Dirichlet last substep, l_f is not in correct cell"
                            stop
                        end if
                        timePassed = del_t 
                        particleList(j)%phaseSpace(1, i-delIdx, iThread) = l_f
                        particleList(j)%phaseSpace(2,i-delIdx, iThread) = v_f
                        particleList(j)%phaseSpace(3:4, i-delIdx, iThread) = particleList(j)%phaseSpace(3:4, i, iThread) 
                    else
                        v_f = 2.0d0 * (l_f - l_sub) * world%nodeVol(l_cell) / del_tau - v_sub
                        if (l_f == l_cell) then 
                            delIdx = delIdx + 1
                            if (l_f == 1) then
                                particleList(j)%energyLoss(1, iThread) = particleList(j)%energyLoss(1, iThread) + (v_f**2 + SUM(particleList(j)%phaseSpace(3:4, i, iThread)**2))!J/m^2 in 1D
                                particleList(j)%wallLoss(1, iThread) = particleList(j)%wallLoss(1, iThread) + 1 !C/m^2 in 1D
                            else if (l_f == NumberXNodes) then
                                particleList(j)%energyLoss(2, iThread) = particleList(j)%energyLoss(2, iThread) + (v_f**2 + SUM(particleList(j)%phaseSpace(3:4, i, iThread)**2)) !J/m^2 in 1D
                                particleList(j)%wallLoss(2, iThread) = particleList(j)%wallLoss(2, iThread) + 1 !C/m^2 in 1D
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
                    l_alongV = l_sub + SIGN(0.5d0, v_sub)
                    call getDelTauSubStep(solver,world, particleList(j), l_sub, l_f, v_sub, del_tau, timePassed, del_t, l_alongV, l_cell, a, c, getEFieldNeumann)
                    if (del_tau >= del_t-timePassed) then
                        ! Add directly to J with no substep
                        call analyticalParticleMoverNeumann(solver,world, particleList(j)%q, particleList(j)%mass, l_sub, l_f, v_sub, v_f, l_cell, del_t, timePassed)
                        !call solver%picardIterParticles(world, particleList(j)%q, particleList(j)%mass, l_sub, l_f, v_sub, v_f, l_cell, del_t, timePassed)
                        if (NINT(l_f) /= l_cell) then
                            print *, "After ongoing substep, l_f is not in correct cell"
                        end if
                        timePassed = del_t 
                        particleList(j)%phaseSpace(1, i-delIdx, iThread) = l_f
                        particleList(j)%phaseSpace(2,i-delIdx, iThread) = v_f
                        particleList(j)%phaseSpace(3:4, i-delIdx, iThread) = particleList(j)%phaseSpace(3:4, i, iThread) 
                    else
                        v_f = 2.0d0 * (l_f - l_sub) * world%nodeVol(l_cell) / del_tau - v_sub
                        if (l_f == l_cell) then
                            v_f = -v_f
                            particleList(j)%refIdx(iThread) = particleList(j)%refIdx(iThread) + 1
                            particleList(j)%refRecordIdx(particleList(j)%refIdx(iThread), iThread) = i - delIdx
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
                CASE(3)
                    l_alongV = l_sub + SIGN(1.0d0, v_sub)
                    call getDelTauSubStep(solver,world, particleList(j), l_sub, l_f, v_sub, del_tau, timePassed, del_t, l_alongV, l_cell, a, c, getEFieldPeriodic)
                    if (del_tau >= del_t-timePassed) then
                        ! Add directly to J with no substep
                        call analyticalParticleMoverPeriodic(solver,world, particleList(j)%q, particleList(j)%mass, l_sub, l_f, v_sub, v_f, l_cell, del_t, timePassed)
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
                        particleList(j)%phaseSpace(1, i-delIdx, iThread) = l_f
                        particleList(j)%phaseSpace(2,i-delIdx, iThread) = v_f
                        particleList(j)%phaseSpace(3:4, i-delIdx, iThread) = particleList(j)%phaseSpace(3:4, i, iThread)  
                    else
                        v_f = 2.0d0 * (l_f - l_sub) * world%nodeVol(l_cell) / del_tau - v_sub
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
        particleList(j)%N_p(iThread) = particleList(j)%N_p(iThread) - delIdx
        particleList(j)%delIdx(iThread) = delIdx
        !$OMP end parallel
        particleList(j)%accumEnergyLoss = particleList(j)%accumEnergyLoss + SUM(particleList(j)%energyLoss, DIM = 2)
        particleList(j)%accumWallLoss = particleList(j)%accumWallLoss + SUM(particleList(j)%wallLoss, DIM = 2)
    end do loopSpecies
    end subroutine moveParticles

end module mod_particleMover