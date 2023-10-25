module mod_particleMover
    use iso_fortran_env, only: int32, real64, output_unit
    use constants
    use mod_BasicFunctions
    use mod_particle
    use mod_domain
    use mod_potentialSolver
    use mod_particleInjection
    use omp_lib
    implicit none
    ! Procedures for moving particles and depositing J 
    integer(int32), parameter :: m_Anderson_Particle = 2

contains

    subroutine subStepSolverAA(l_sub, v_sub, l_f, v_f, v_half, del_tau, timePassed, del_t, E_x, q_over_m, l_cell, BField, B_mag, dx_dl, f_tol, AtBoundaryBool)
        real(real64), intent(in) :: l_sub, v_sub(3), del_t, E_x, q_over_m, BField(3), B_mag, dx_dl, f_tol
        integer(int32), intent(in) :: l_cell
        real(real64), intent(in out) :: l_f, v_f(3), timePassed, del_tau, v_half(3)
        logical, intent(in out) :: AtBoundaryBool
        logical :: convergeBool = .false.
        real(real64) :: Res_k(m_Anderson_Particle+1), del_tau_k(m_Anderson_Particle+1), fitMat(m_Anderson_Particle), minCoeff(m_Anderson_Particle+1), del_tau_max
        integer(int32) :: m_k, index, k, l_boundary, u
        real(real64) :: coeffAccel, v_prime(3), curr_Res
        ! print *, ''
        ! print *, 'ENTER SUBSTEP SOLVER'
        ! print *, "--------------"
        ! print *, 'initial v_sub:', v_sub
        ! print *, 'initial l_sub:', l_sub
        ! print *, 'AtBoundaryBool:', AtBoundaryBool
        ! print *, 'l_cell:', l_cell
        ! print *, 'E_x:', E_x
        v_half = v_sub
        v_prime = v_sub
        del_tau_max = del_t - timePassed
        del_tau_k(1) = del_tau_max
        k = 0
        convergeBool = .false.
        do while (.not. convergeBool)
            index = MODULO(k, m_Anderson_Particle+1) + 1
            m_k = MIN(k, m_Anderson_Particle)
            coeffAccel = 0.5d0 * del_tau_k(index) * q_over_m
            v_prime(1) = v_sub(1) + coeffAccel * E_x
            v_half = v_prime + coeffAccel * (crossProduct(v_prime, BField) + coeffAccel* SUM(v_prime * BField) * BField)
            v_half = v_half / (1.0d0 + (coeffAccel*B_mag)**2)
            if (v_half(1) > 0) then
                l_boundary = l_cell + 1
            else
                l_boundary = l_cell
            end if
            del_tau = MIN(del_tau_max, (real(l_boundary) - l_sub) * dx_dl / v_half(1))
            Res_k(index) = del_tau - del_tau_k(index)
            print *, 'k is:', k, 'with del_tau_k:', del_tau_k
            curr_Res = ABS(Res_k(index))
            if (curr_Res < f_tol) then
                if (del_tau < del_tau_max) then
                    AtBoundaryBool = .true.
                    l_f = real(l_boundary)
                else
                    AtBoundaryBool = .false.
                    l_f = l_sub + v_half(1) * del_tau / dx_dl
                end if
                v_f = 2.0d0 * v_half - v_sub
                timePassed = timePassed + del_tau
                convergeBool = .true.
                exit
            end if
            if (k > m_Anderson_Particle) then
                if (curr_Res >= SUM(abs(Res_k))/real(m_Anderson_Particle+1)) then
                    del_tau_max = 0.5d0 * del_tau_max
                    del_tau_k(1) = del_tau_max
                    k = -1
                end if
            end if
            if (k > 0) then
                do u = 0, m_k-1
                    fitMat(u+1) = Res_k(MODULO(k - m_k + u, m_Anderson_Particle+1) + 1) - Res_k(index)
                end do
                minCoeff(1:m_k) = (fitMat(1:m_k) * -Res_k(index))/SUM(fitMat(1:m_k)**2)
                minCoeff(m_k+1) = 1.0d0 - SUM(minCoeff(1:m_k))
                del_tau_k(MODULO(k+1, m_Anderson_Particle+1) + 1) = minCoeff(1) * (Res_k(MODULO(k-m_k, m_Anderson_Particle+1) + 1) + del_tau_k(MODULO(k-m_k,m_Anderson_Particle+1) + 1))
                do u=1, m_k
                    del_tau_k(MODULO(k+1, m_Anderson_Particle+1) + 1) = del_tau_k(MODULO(k+1, m_Anderson_Particle+1) + 1) + minCoeff(u + 1) * (Res_k(MODULO(k-m_k + u, m_Anderson_Particle+1) + 1) + del_tau_k(MODULO(k-m_k + u, m_Anderson_Particle+1) + 1))
                end do
            else
                del_tau_k(2) = del_tau
            end if
            k = k + 1
        end do
        ! print *, 'v_f is:', v_f
        ! print *, 'l_f is:', l_f
        ! print *, 'k is:', k
        ! print *, ''

    end subroutine subStepSolverAA

    subroutine depositJ(solver, particleList, world, del_t)
        ! particle substepping procedure which deposits J, also used for general moving of particles
        ! boolDepositJ = true : deposit J
        ! boolDepositJ = false: only move particles (done after non-linear iteration), also delete particles for wall collisions (saves extra loop in collisions)
        class(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        type(Particle), intent(in out) :: particleList(:)
        real(real64), intent(in) :: del_t
        !a and c correspond to quadratic equations | l_alongV is nearest integer boundary along velocity component, away is opposite
        real(real64) :: l_f, l_sub, v_sub(3), v_f(3), timePassed, del_tau, q_over_m, f_tol, v_half(3), dx_dl, E_x, q_times_wp
        integer(int32) :: j, i, l_cell, iThread
        logical :: AtBoundaryBool
        logical :: convergeBool = .false.
        real(real64) :: del_tau_max, del_tau_prev
        integer(int32) :: k, l_boundary
        real(real64) :: coeffAccel, v_prime(3), curr_Res, prev_Res
        solver%J = 0.0d0
        solver%EField = 0.5d0 * (solver%phi(1:NumberXNodes-1) + solver%phi_f(1:NumberXNodes-1) - solver%phi(2:NumberXNodes) - solver%phi_f(2:NumberXNodes)) / world%dx_dl
        f_tol = del_t * 1.d-8
        loopSpecies: do j = 1, numberChargedParticles
            q_over_m = particleList(j)%q/particleList(j)%mass
            q_times_wp = particleList(j)%q * particleList(j)%w_p
            !$OMP parallel private(iThread, i, l_f, l_sub, v_sub, v_f, v_half, timePassed, del_tau, l_cell, AtBoundaryBool, dx_dl, E_x, &
            convergeBool, del_tau_max, del_tau_prev, k, l_boundary, coeffAccel, v_prime, curr_Res, prev_Res)
            iThread = omp_get_thread_num() + 1 
            loopParticles: do i = 1, particleList(j)%N_p(iThread)
                v_sub = particleList(j)%phaseSpace(2:4,i,iThread)
                l_sub = particleList(j)%phaseSpace(1,i,iThread)
                timePassed = 0.0d0
                do while((timePassed < del_t))
                    l_cell = INT(l_sub)
                    E_x = solver%EField(l_cell)
                    dx_dl = world%dx_dl(l_cell)
                    v_half = v_sub
                    v_prime = v_sub
                    del_tau_max = del_t - timePassed
                    del_tau_prev = del_tau_max
                    prev_Res = del_tau_max
                    k = 0
                    convergeBool = .false.
                    do while (.not. convergeBool)
                        coeffAccel = 0.5d0 * del_tau_prev * q_over_m
                        v_prime(1) = v_sub(1) + coeffAccel * E_x
                        v_half = v_prime + coeffAccel * (crossProduct(v_prime, solver%BField) + coeffAccel* SUM(v_prime * solver%BField) * solver%BField)
                        v_half = v_half / (1.0d0 + (coeffAccel*solver%BFieldMag)**2)
                        if (v_half(1) > 0) then
                            l_boundary = l_cell + 1
                        else
                            l_boundary = l_cell
                        end if
                        del_tau = MIN(del_tau_max, (real(l_boundary) - l_sub) * dx_dl / v_half(1))
                        curr_Res = ABS(del_tau - del_tau_prev)
                        if (curr_Res < f_tol) then
                            v_f = 2.0d0 * v_half - v_sub
                            AtBoundaryBool = (del_tau < del_tau_max)
                            if(.not. AtBoundaryBool) then
                                l_f = l_sub + v_half(1) * del_tau / dx_dl
                                if (MOD(l_f, 1.0d0) == 0.0d0) then
                                    AtBoundaryBool = .true.
                                    l_boundary = NINT(l_f)
                                end if
                                if (INT(l_f) /= l_cell) then
                                    print *, 'l_f not in cell after going del_tau_max!'
                                    print *, 'del_tau is:', del_tau
                                    print *, 'l_sub:', l_sub
                                    print *, 'l_f:', l_f
                                    print *, 'v_sub:', v_sub
                                    print *, 'v_half:', v_half
                                    print *, 'v_f:',v_f
                                    print *, 'del_tau_max:', del_tau_max
                                    stop
                                end if
                            end if
                            if (del_tau <= 0.0d0) then
                                print *, 'AtBoundaryBool:', AtBoundaryBool
                                print *, 'del_tau is 0 or less'
                                print *, 'del_tau is:', del_tau
                                print *, 'l_sub:', l_sub
                                print *, 'l_f:', l_f
                                print *, 'v_sub:', v_sub
                                print *, 'v_half:', v_half
                                print *, 'v_f:',v_f
                                print *, 'del_tau_max:', del_tau_max
                                stop
                            end if
                            timePassed = timePassed + del_tau
                            convergeBool = .true.
                            exit
                        end if

                        if (curr_Res >= prev_Res) then
                            del_tau_max = 0.5d0 * del_tau_max
                            del_tau = del_tau_max
                            k = -1
                        end if
                        del_tau_prev = del_tau
                        prev_Res = curr_Res
                        k = k + 1
                    end do
                    solver%J(l_cell, iThread) = solver%J(l_cell, iThread) + q_times_wp * (v_half(1))*del_tau/world%dx_dl(l_cell)/del_t
                    
                    if (AtBoundaryBool) then
                        SELECT CASE (world%boundaryConditions(l_boundary))
                        CASE(0)
                            l_f = real(l_boundary) + SIGN(1.d-12, v_f(1))
                        CASE(1,4)
                            exit
                        CASE(2)
                            if (l_boundary == NumberXNodes) then
                                v_f(1) = -ABS(v_f(1))
                            else
                                v_f(1) = ABS(v_f(1))
                            end if
                            v_f(2:3) = -v_f(2:3)
                            l_f = real(l_boundary) + SIGN(1.d-12, v_f(1))
                        CASE(3)
                            l_f = REAL(ABS(l_boundary - real(NumberXNodes, kind = real64) - 1)) + SIGN(1.d-12, v_f(1))
                        CASE default
                            print *, "l_sub is:", l_sub
                            print *, 'l_f is:', l_f
                            print *, world%boundaryConditions(INT(l_f))
                            print *, "Case does not exist in ongoing substep, depositJ"
                            stop
                        END SELECT
                    end if
                    ! now final position/velocity becomes next starting position/velocity
                    l_sub = l_f
                    v_sub = v_f
                    if ((l_f < 1) .or. (l_f > NumberXNodes)) then
                        print *, 'l_f is:', l_f
                        stop "Have particles travelling outside the domain in depositJ!"
                    end if
                end do
            end do loopParticles
            !$OMP end parallel
        end do loopSpecies
        
    end subroutine depositJ

    subroutine depositJOther(solver, particleList, world, del_t)
        ! particle substepping procedure which deposits J, also used for general moving of particles
        ! boolDepositJ = true : deposit J
        ! boolDepositJ = false: only move particles (done after non-linear iteration), also delete particles for wall collisions (saves extra loop in collisions)
        class(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        type(Particle), intent(in out) :: particleList(:)
        real(real64), intent(in) :: del_t
        !a and c correspond to quadratic equations | l_alongV is nearest integer boundary along velocity component, away is opposite
        real(real64) :: l_f, l_sub, v_sub(3), v_f(3), timePassed, del_tau, q_over_m, f_tol, v_half(3), dx_dl, E_x
        integer(int32) :: j, i, l_cell, iThread
        logical :: AtBoundaryBool
        logical :: convergeBool = .false.
        real(real64) :: Res_k(m_Anderson_Particle+1), del_tau_k(m_Anderson_Particle+1), fitMat(m_Anderson_Particle), minCoeff(m_Anderson_Particle+1), del_tau_max
        integer(int32) :: m_k, index, k, l_boundary, u
        real(real64) :: coeffAccel, v_prime(3), curr_Res
        solver%J = 0.0d0
        solver%EField = 0.5d0 * (solver%phi(1:NumberXNodes-1) + solver%phi_f(1:NumberXNodes-1) - solver%phi(2:NumberXNodes) - solver%phi_f(2:NumberXNodes)) / world%dx_dl
        f_tol = del_t * 1.d-8
        loopSpecies: do j = 1, numberChargedParticles
            q_over_m = particleList(j)%q/particleList(j)%mass
            !$OMP parallel private(iThread, i, l_f, l_sub, v_sub, v_f, v_half, timePassed, del_tau, l_cell, AtBoundaryBool, dx_dl, E_x, &
            convergeBool, Res_k, del_tau_k, fitMat, minCoeff, del_tau_max, m_k, index, k, l_boundary, u, coeffAccel, v_prime, curr_Res)
            iThread = omp_get_thread_num() + 1 
            loopParticles: do i = 1, particleList(j)%N_p(iThread)
                v_sub = particleList(j)%phaseSpace(2:4,i,iThread)
                l_sub = particleList(j)%phaseSpace(1,i,iThread)
                timePassed = 0.0d0
                AtBoundaryBool = (MODULO(l_sub, 1.0d0) == 0)
                if (AtBoundaryBool) then
                    print *, 'Initial particle location at boundary!'
                    print *, 'l_sub:', l_sub
                    print *, 'v_sub:', v_sub
                    stop
                end if
                do while((timePassed < del_t))
                    if (.not. AtBoundaryBool) then
                        l_cell = INT(l_sub)
                    else
                        l_cell = INT(l_sub + SIGN(0.5d0, v_sub(1)))
                        l_boundary = INT(l_sub + SIGN(1.0d0, v_sub(1)))
                    end if
                    E_x = solver%EField(l_cell)
                    dx_dl = world%dx_dl(l_cell)
                    v_half = v_sub
                    v_prime = v_sub
                    del_tau_max = del_t - timePassed
                    del_tau_k(1) = del_tau_max
                    k = 0
                    convergeBool = .false.
                    do while (.not. convergeBool)
                        index = MODULO(k, m_Anderson_Particle+1) + 1
                        m_k = MIN(k, m_Anderson_Particle)
                        coeffAccel = 0.5d0 * del_tau_k(index) * q_over_m
                        v_prime(1) = v_sub(1) + coeffAccel * E_x
                        v_half = v_prime + coeffAccel * (crossProduct(v_prime, solver%BField) + coeffAccel* SUM(v_prime * solver%BField) * solver%BField)
                        v_half = v_half / (1.0d0 + (coeffAccel*solver%BFieldMag)**2)
                        if (.not. AtBoundaryBool) then
                            if (v_half(1) > 0) then
                                l_boundary = l_cell + 1
                            else
                                l_boundary = l_cell
                            end if
                            del_tau = MIN(del_tau_max, (real(l_boundary) - l_sub) * dx_dl / v_half(1))
                            if (del_tau == 0.0d0) then
                                print *, 'del_tau is 0 for not at boundary'
                                print *, 'del_tau in iteration is:', del_tau
                                print *, 'del_tau_max:', del_tau_max
                                print *, 'l_cell:', l_cell
                                print *, 'l_sub:', l_sub
                                print *, 'l_sub - int(l_sub):', l_sub - l_cell
                                print *, 'l_boundary:', l_boundary
                                print *, 'v_sub:', v_sub
                                print *, 'v_half:', v_half
                                print *, 'AtBoundaryBool:', AtBoundaryBool
                                stop
                            end if
                        else
                            if (INT(SIGN(1.0d0, v_half(1))) == INT(SIGN(1.0d0, v_sub(1)))) then
                                del_tau = MIN(del_tau_max, (real(l_boundary) - l_sub) * dx_dl / v_half(1))
                            else
                                del_tau_max = del_tau_max * 0.5d0
                                del_tau = del_tau_max
                            end if
                            if (del_tau == 0.0d0) then
                                print *, 'del_tau is 0 for at boundary'
                                print *, 'del_tau in iteration is:', del_tau
                                print *, 'del_tau_max:', del_tau_max
                                print *, 'l_sub:', l_sub
                                print *, 'v_sub:', v_sub
                                print *, 'v_half:', v_half
                                print *, 'AtBoundaryBool:', AtBoundaryBool
                                stop
                            end if
                        end if
                        Res_k(index) = del_tau - del_tau_k(index)
                        curr_Res = ABS(Res_k(index))
                        if (curr_Res < f_tol) then
                            v_f = 2.0d0 * v_half - v_sub
                            l_f = l_sub + v_half(1) * del_tau / dx_dl
                            ! if (del_tau <= 0.0d0) then
                            !     print *, 'AtBoundaryBool:', AtBoundaryBool
                            !     print *, 'del_tau is 0 or less'
                            !     print *, 'del_tau is:', del_tau
                            !     print *, 'l_sub:', l_sub
                            !     print *, 'l_f:', l_f
                            !     print *, 'v_sub:', v_sub
                            !     print *, 'v_half:', v_half
                            !     print *, 'v_f:',v_f
                            !     print *, 'del_tau_max:', del_tau_max
                            !     stop
                            ! end if
                            ! AtBoundaryBool = (del_tau < del_tau_max)
                            ! if (AtBoundaryBool) then
                            !     l_f = real(l_boundary)
                            !     if (l_f > NumberXNodes .or. l_f < 1) then
                            !         print *, 'l_f not in domain!'
                            !         print *, 'l_f:', l_f
                            !         print *, 'l_sub:', l_sub
                            !         print *, 'v_sub:', v_sub
                            !         print *, 'v_half:', v_half
                            !         print *, 'del_tau:', del_tau
                            !         print *, 'del_tau_max:', del_tau_max
                            !         stop
                            !     end if
                            ! else
                            !     l_f = l_sub + v_half(1) * del_tau / dx_dl
                            !     if (MOD(l_f, 1.0d0) == 0.0d0) then
                            !         print *, 'l_f should not be at boundary, but is at boundary!'
                            !         print *, 'l_f:', l_f
                            !         print *, 'l_sub:', l_sub
                            !         print *, 'v_sub:', v_sub
                            !         print *, 'v_half:', v_half
                            !         print *, 'v_half(1) * del_tau / dx_dl', v_half(1) * del_tau / dx_dl
                            !         print *, 'del_tau:', del_tau
                            !         print *, 'del_tau_max:', del_tau_max
                            !         stop
                            !     end if
                            ! end if
                            timePassed = timePassed + del_tau
                            AtBoundaryBool = (ABS(l_f - real(NINT(l_f))) < 1.d-12)
                            if ((.not. AtBoundaryBool) .and. (del_tau < del_tau_max)) then
                                print *, 'particle not at boundary, yet del_tau < del_tau_max'
                                print *, 'l_boundary:', l_boundary
                                print *, 'del_tau is:', del_tau
                                print *, 'l_sub:', l_sub
                                print *, 'l_f:', l_f
                                print *, '(ABS(l_f - real(l_boundary))', ABS(l_f - real(l_boundary))
                                print *, 'v_sub:', v_sub
                                print *, 'v_half:', v_half
                                print *, 'v_f:',v_f
                                print *, 'del_tau_max:', del_tau_max
                                stop
                            end if
                            convergeBool = .true.
                            exit
                        end if
                        if (k > m_Anderson_Particle) then
                            if (curr_Res >= SUM(abs(Res_k))/real(m_Anderson_Particle+1)) then
                                del_tau_max = 0.5d0 * del_tau_max
                                del_tau_k(1) = del_tau_max
                                k = -1
                            end if
                        end if
                        if (k > 0) then
                            do u = 0, m_k-1
                                fitMat(u+1) = Res_k(MODULO(k - m_k + u, m_Anderson_Particle+1) + 1) - Res_k(index)
                            end do
                            minCoeff(1:m_k) = (fitMat(1:m_k) * -Res_k(index))/SUM(fitMat(1:m_k)**2)
                            minCoeff(m_k+1) = 1.0d0 - SUM(minCoeff(1:m_k))
                            del_tau_k(MODULO(k+1, m_Anderson_Particle+1) + 1) = minCoeff(1) * (Res_k(MODULO(k-m_k, m_Anderson_Particle+1) + 1) + del_tau_k(MODULO(k-m_k,m_Anderson_Particle+1) + 1))
                            do u=1, m_k
                                del_tau_k(MODULO(k+1, m_Anderson_Particle+1) + 1) = del_tau_k(MODULO(k+1, m_Anderson_Particle+1) + 1) + minCoeff(u + 1) * (Res_k(MODULO(k-m_k + u, m_Anderson_Particle+1) + 1) + del_tau_k(MODULO(k-m_k + u, m_Anderson_Particle+1) + 1))
                            end do
                        else
                            del_tau_k(2) = del_tau
                        end if
                        k = k + 1
                    end do
                    solver%J(l_cell, iThread) = solver%J(l_cell, iThread) + particleList(j)%w_p * particleList(j)%q * (v_half(1))*del_tau/world%dx_dl(l_cell)/del_t
                    
                    if (AtBoundaryBool) then
                        l_f = real(NINT(l_f))
                        SELECT CASE (world%boundaryConditions(INT(l_f)))
                        CASE(0)
                            continue
                        CASE(1,4)
                            exit
                        CASE(2)
                            if (INT(l_f) == NumberXNodes) then
                                v_f(1) = -ABS(v_f(1))
                            else
                                v_f(1) = ABS(v_f(1))
                            end if
                            v_f(2:3) = -v_f(2:3)
                        CASE(3)
                            l_f = ABS(l_f - real(NumberXNodes, kind = real64) - 1.0d0)
                        CASE default
                            print *, "l_sub is:", l_sub
                            print *, 'l_f is:', l_f
                            print *, world%boundaryConditions(INT(l_f))
                            print *, "Case does not exist in ongoing substep, depositJ"
                            stop
                        END SELECT
                    end if
                    ! now final position/velocity becomes next starting position/velocity
                    if ((l_f == NumberXNodes) .and. (v_f(1) > 0.0d0)) then
                        print *, 'l_f is at edge and v_f > 0'
                        print *, 'AtBoundaryBool:', AtBoundaryBool
                        print *, "l_sub is:", l_sub
                        print *, 'v_sub is:', v_sub
                        print *, 'l_f is:', l_f
                        print *, 'v_f is:', v_f
                        print *, 'del_tau is:', del_tau
                        print *, 'del_tau_max is:', del_tau_max
                        stop
                    end if
                    l_sub = l_f
                    v_sub = v_f
                    if ((l_f < 1) .or. (l_f > NumberXNodes)) then
                        print *, 'l_f is:', l_f
                        stop "Have particles travelling outside the domain in depositJ!"
                    end if
                end do
            end do loopParticles
            !$OMP end parallel
        end do loopSpecies
        
    end subroutine depositJOther


    ! -------------------------------------------- Particle mover without boolean checks for depositing J ------------------------------------------------------------

    subroutine moveParticles(solver, particleList, world, del_t)
        ! particle mover to avoid the boolean checks which mostly don't happen when depositing J
        class(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        type(Particle), intent(in out) :: particleList(:)
        real(real64), intent(in) :: del_t
        !a and c correspond to quadratic equations | l_alongV is nearest integer boundary along velocity component, away is opposite
        real(real64) :: l_f, l_sub, v_sub(3), v_f(3), timePassed, del_tau, q_over_m, f_tol, v_half(3), dx_dl, E_x
        integer(int32) :: j, i, l_cell, iThread, delIdx
        logical :: AtBoundaryBool
        logical :: convergeBool = .false.
        real(real64) :: del_tau_max, del_tau_prev
        integer(int32) :: k, l_boundary
        real(real64) :: coeffAccel, v_prime(3), curr_Res, prev_Res
        solver%EField = 0.5d0 * (solver%phi(1:NumberXNodes-1) + solver%phi_f(1:NumberXNodes-1) - solver%phi(2:NumberXNodes) - solver%phi_f(2:NumberXNodes)) / world%dx_dl
        f_tol = del_t * 1.d-8
        loopSpecies: do j = 1, numberChargedParticles
            q_over_m = particleList(j)%q/particleList(j)%mass
            !$OMP parallel private(iThread, i, l_f, l_sub, v_sub, v_f, v_half, timePassed, del_tau, l_cell, AtBoundaryBool, delIdx, dx_dl, E_x, &
            convergeBool, del_tau_max, del_tau_prev, k, l_boundary, coeffAccel, v_prime, curr_Res, prev_Res)
            iThread = omp_get_thread_num() + 1 
            delIdx = 0
            particleList(j)%refIdx(iThread) = 0
            particleList(j)%energyLoss(:, iThread) = 0.0d0
            particleList(j)%wallLoss(:, iThread) = 0.0d0
            loopParticles: do i = 1, particleList(j)%N_p(iThread)
                v_sub = particleList(j)%phaseSpace(2:4,i,iThread)
                l_sub = particleList(j)%phaseSpace(1,i,iThread)
                timePassed = 0.0d0
                do while((timePassed < del_t))
                    l_cell = INT(l_sub)
                    E_x = solver%EField(l_cell)
                    dx_dl = world%dx_dl(l_cell)
                    v_half = v_sub
                    v_prime = v_sub
                    del_tau_max = del_t - timePassed
                    del_tau_prev = del_tau_max
                    prev_Res = del_tau_max
                    k = 0
                    convergeBool = .false.
                    do while (.not. convergeBool)
                        coeffAccel = 0.5d0 * del_tau_prev * q_over_m
                        v_prime(1) = v_sub(1) + coeffAccel * E_x
                        v_half = v_prime + coeffAccel * (crossProduct(v_prime, solver%BField) + coeffAccel* SUM(v_prime * solver%BField) * solver%BField)
                        v_half = v_half / (1.0d0 + (coeffAccel*solver%BFieldMag)**2)
                        if (v_half(1) > 0) then
                            l_boundary = l_cell + 1
                        else
                            l_boundary = l_cell
                        end if
                        del_tau = MIN(del_tau_max, (real(l_boundary) - l_sub) * dx_dl / v_half(1))
                        curr_Res = ABS(del_tau - del_tau_prev)
                        if (curr_Res < f_tol) then
                            v_f = 2.0d0 * v_half - v_sub
                            AtBoundaryBool = (del_tau < del_tau_max)
                            if(.not. AtBoundaryBool) then
                                l_f = l_sub + v_half(1) * del_tau / dx_dl
                                if (MOD(l_f, 1.0d0) == 0.0d0) then
                                    AtBoundaryBool = .true.
                                    l_boundary = NINT(l_f)
                                end if
                            end if
                            timePassed = timePassed + del_tau
                            convergeBool = .true.
                            exit
                        end if

                        if (curr_Res >= prev_Res) then
                            del_tau_max = 0.5d0 * del_tau_max
                            del_tau = del_tau_max
                            k = -1
                        end if
                        del_tau_prev = del_tau
                        prev_Res = curr_Res
                        k = k + 1
                    end do
                    ! End AA-Picard
                    
                    if (AtBoundaryBool) then
                        SELECT CASE (world%boundaryConditions(l_boundary))
                        CASE(0)
                            l_f = real(l_boundary) + SIGN(1.d-12, v_f(1))
                        CASE(1)
                            delIdx = delIdx + 1
                            if (l_boundary == 1) then
                                particleList(j)%energyLoss(1, iThread) = particleList(j)%energyLoss(1, iThread) + (SUM(v_f**2))!J/m^2 in 1D
                                particleList(j)%wallLoss(1, iThread) = particleList(j)%wallLoss(1, iThread) + 1 !C/m^2 in 1D
                            else if (l_boundary == NumberXNodes) then
                                particleList(j)%energyLoss(2, iThread) = particleList(j)%energyLoss(2, iThread) + (SUM(v_f**2)) !J/m^2 in 1D
                                particleList(j)%wallLoss(2, iThread) = particleList(j)%wallLoss(2, iThread) + 1 !C/m^2 in 1D
                            end if
                            exit
                        CASE(2)
                            if (l_boundary == NumberXNodes) then
                                v_f(1) = -ABS(v_f(1))
                            else
                                v_f(1) = ABS(v_f(1))
                            end if
                            v_f(2:3) = -v_f(2:3)
                            l_f = real(l_boundary) + SIGN(1.d-12, v_f(1))
                            particleList(j)%refIdx(iThread) = particleList(j)%refIdx(iThread) + 1
                            particleList(j)%refRecordIdx(particleList(j)%refIdx(iThread), iThread) = i - delIdx
                        CASE(3)
                            l_f = REAL(ABS(l_boundary - real(NumberXNodes, kind = real64) - 1)) + SIGN(1.d-12, v_f(1))
                        CASE(4)
                            if (refluxPartBool .or. injectionBool) then
                                particleList(j)%phaseSpace(1, i-delIdx, iThread) = l_f
                                particleList(j)%phaseSpace(2:4, i-delIdx, iThread) = v_f
                                particleList(j)%refIdx(iThread) = particleList(j)%refIdx(iThread) + 1
                                particleList(j)%refRecordIdx(particleList(j)%refIdx(iThread), iThread) = i - delIdx
                            else
                                delIdx = delIdx + 1
                                if (l_f == 1) then
                                    particleList(j)%energyLoss(1, iThread) = particleList(j)%energyLoss(1, iThread) + (SUM(v_f**2))!J/m^2 in 1D
                                    particleList(j)%wallLoss(1, iThread) = particleList(j)%wallLoss(1, iThread) + 1 !C/m^2 in 1D
                                else if (l_f == NumberXNodes) then
                                    particleList(j)%energyLoss(2, iThread) = particleList(j)%energyLoss(2, iThread) + (SUM(v_f**2)) !J/m^2 in 1D
                                    particleList(j)%wallLoss(2, iThread) = particleList(j)%wallLoss(2, iThread) + 1 !C/m^2 in 1D
                                end if
                                exit
                            end if
                        CASE default
                            print *, "l_sub is:", l_sub
                            print *, 'l_f is:', l_f
                            print *, world%boundaryConditions(INT(l_f))
                            print *, "Case does not exist in ongoing substep, depositJ"
                            stop
                        END SELECT
                    else
                        particleList(j)%phaseSpace(1, i-delIdx, iThread) = l_f
                        particleList(j)%phaseSpace(2:4,i-delIdx, iThread) = v_f
                    end if
                    if ((l_f < 1) .or. (l_f > NumberXNodes)) then
                        print *, "Have particles travelling outside domain in moveparticles!"
                        print *, 'del_tau is:', del_tau
                        print *, 'l_sub:', l_sub
                        print *, 'l_f:', l_f
                        print *, 'v_sub:', v_sub
                        print *, 'v_half:', v_half
                        print *, 'v_f:',v_f
                        print *, 'del_tau_max:', del_tau_max
                        stop
                    end if
                    ! now final position/velocity becomes next starting position/velocity
                    l_sub = l_f
                    v_sub = v_f
                end do
                
            end do loopParticles
            particleList(j)%N_p(iThread) = particleList(j)%N_p(iThread) - delIdx
            particleList(j)%delIdx(iThread) = delIdx
            !$OMP end parallel
            particleList(j)%accumEnergyLoss = particleList(j)%accumEnergyLoss + SUM(particleList(j)%energyLoss, DIM = 2)
            particleList(j)%accumWallLoss = particleList(j)%accumWallLoss + SUM(particleList(j)%wallLoss, DIM = 2)
        end do loopSpecies
    end subroutine moveParticles

    subroutine moveParticlesOther(solver, particleList, world, del_t)
        ! particle mover to avoid the boolean checks which mostly don't happen when depositing J
        class(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        type(Particle), intent(in out) :: particleList(:)
        real(real64), intent(in) :: del_t
        !a and c correspond to quadratic equations | l_alongV is nearest integer boundary along velocity component, away is opposite
        real(real64) :: l_f, l_sub, v_sub(3), v_f(3), timePassed, del_tau, q_over_m, f_tol, v_half(3), dx_dl, E_x
        integer(int32) :: j, i, l_cell, iThread, delIdx
        logical :: AtBoundaryBool
        logical :: convergeBool = .false.
        real(real64) :: Res_k(m_Anderson_Particle+1), del_tau_k(m_Anderson_Particle+1), fitMat(m_Anderson_Particle), minCoeff(m_Anderson_Particle+1), del_tau_max
        integer(int32) :: m_k, index, k, l_boundary, u
        real(real64) :: coeffAccel, v_prime(3), curr_Res
        solver%EField = 0.5d0 * (solver%phi(1:NumberXNodes-1) + solver%phi_f(1:NumberXNodes-1) - solver%phi(2:NumberXNodes) - solver%phi_f(2:NumberXNodes)) / world%dx_dl
        f_tol = del_t * 1.d-8
        loopSpecies: do j = 1, numberChargedParticles
            q_over_m = particleList(j)%q/particleList(j)%mass
            !$OMP parallel private(iThread, i, l_f, l_sub, v_sub, v_f, v_half, timePassed, del_tau, l_cell, AtBoundaryBool, delIdx, dx_dl, E_x, &
            convergeBool, Res_k, del_tau_k, fitMat, minCoeff, del_tau_max, m_k, index, k, l_boundary, u, coeffAccel, v_prime, curr_Res)
            iThread = omp_get_thread_num() + 1 
            delIdx = 0
            particleList(j)%refIdx(iThread) = 0
            particleList(j)%energyLoss(:, iThread) = 0.0d0
            particleList(j)%wallLoss(:, iThread) = 0.0d0
            loopParticles: do i = 1, particleList(j)%N_p(iThread)
                v_sub = particleList(j)%phaseSpace(2:4,i,iThread)
                l_sub = particleList(j)%phaseSpace(1,i,iThread)
                timePassed = 0.0d0
                AtBoundaryBool = (MODULO(l_sub, 1.0d0) == 0)
                do while((timePassed < del_t))
                    if (.not. AtBoundaryBool) then
                        l_cell = INT(l_sub)
                    else
                        l_cell = INT(l_sub + SIGN(0.5d0, v_sub(1)))
                        l_boundary = INT(l_sub + SIGN(1.0d0, v_sub(1)))
                    end if

                    ! Anderson-Accelerated Picard
                    E_x = solver%EField(l_cell)
                    dx_dl = world%dx_dl(l_cell)
                    v_half = v_sub
                    v_prime = v_sub
                    del_tau_max = del_t - timePassed
                    del_tau_k(1) = del_tau_max
                    k = 0
                    convergeBool = .false.
                    do while (.not. convergeBool)
                        index = MODULO(k, m_Anderson_Particle+1) + 1
                        m_k = MIN(k, m_Anderson_Particle)
                        coeffAccel = 0.5d0 * del_tau_k(index) * q_over_m
                        v_prime(1) = v_sub(1) + coeffAccel * E_x
                        v_half = v_prime + coeffAccel * (crossProduct(v_prime, solver%BField) + coeffAccel* SUM(v_prime * solver%BField) * solver%BField)
                        v_half = v_half / (1.0d0 + (coeffAccel*solver%BFieldMag)**2)
                        if (.not. AtBoundaryBool) then
                            if (v_half(1) > 0) then
                                l_boundary = l_cell + 1
                            else
                                l_boundary = l_cell
                            end if
                            del_tau = MIN(del_tau_max, (real(l_boundary) - l_sub) * dx_dl / v_half(1))
                        else
                            if (INT(SIGN(1.0d0, v_half(1))) == INT(SIGN(1.0d0, v_sub(1)))) then
                                del_tau = MIN(del_tau_max, (real(l_boundary) - l_sub) * dx_dl / v_half(1))
                            else
                                del_tau_max = del_tau_max * 0.5d0
                                del_tau = del_tau_max
                            end if
                        end if
                        Res_k(index) = del_tau - del_tau_k(index)
                        curr_Res = ABS(Res_k(index))
                        if (curr_Res < f_tol) then
                            v_f = 2.0d0 * v_half - v_sub
                            l_f = l_sub + v_half(1) * del_tau / dx_dl
                            ! AtBoundaryBool = (del_tau < del_tau_max)
                            ! if (AtBoundaryBool) then
                            !     l_f = real(l_boundary)
                            ! else
                            !     l_f = l_sub + v_half(1) * del_tau / dx_dl
                            ! end if
                            timePassed = timePassed + del_tau
                            AtBoundaryBool = (ABS(l_f - real(NINT(l_f))) < 1.d-12)
                            convergeBool = .true.
                            exit
                        end if
                        if (k > m_Anderson_Particle) then
                            if (curr_Res >= SUM(abs(Res_k))/real(m_Anderson_Particle+1)) then
                                del_tau_max = 0.5d0 * del_tau_max
                                del_tau_k(1) = del_tau_max
                                k = -1
                            end if
                        end if
                        if (k > 0) then
                            do u = 0, m_k-1
                                fitMat(u+1) = Res_k(MODULO(k - m_k + u, m_Anderson_Particle+1) + 1) - Res_k(index)
                            end do
                            minCoeff(1:m_k) = (fitMat(1:m_k) * -Res_k(index))/SUM(fitMat(1:m_k)**2)
                            minCoeff(m_k+1) = 1.0d0 - SUM(minCoeff(1:m_k))
                            del_tau_k(MODULO(k+1, m_Anderson_Particle+1) + 1) = minCoeff(1) * (Res_k(MODULO(k-m_k, m_Anderson_Particle+1) + 1) + del_tau_k(MODULO(k-m_k,m_Anderson_Particle+1) + 1))
                            do u=1, m_k
                                del_tau_k(MODULO(k+1, m_Anderson_Particle+1) + 1) = del_tau_k(MODULO(k+1, m_Anderson_Particle+1) + 1) + minCoeff(u + 1) * (Res_k(MODULO(k-m_k + u, m_Anderson_Particle+1) + 1) + del_tau_k(MODULO(k-m_k + u, m_Anderson_Particle+1) + 1))
                            end do
                        else
                            del_tau_k(2) = del_tau
                        end if
                        k = k + 1
                    end do
                    ! End AA-Picard
                    
                    if (AtBoundaryBool) then
                        l_f = real(NINT(l_f))
                        SELECT CASE (world%boundaryConditions(INT(l_f)))
                        CASE(0)
                            continue
                        CASE(1)
                            delIdx = delIdx + 1
                            if (l_f == 1) then
                                particleList(j)%energyLoss(1, iThread) = particleList(j)%energyLoss(1, iThread) + (SUM(v_f**2))!J/m^2 in 1D
                                particleList(j)%wallLoss(1, iThread) = particleList(j)%wallLoss(1, iThread) + 1 !C/m^2 in 1D
                            else if (l_f == NumberXNodes) then
                                particleList(j)%energyLoss(2, iThread) = particleList(j)%energyLoss(2, iThread) + (SUM(v_f**2)) !J/m^2 in 1D
                                particleList(j)%wallLoss(2, iThread) = particleList(j)%wallLoss(2, iThread) + 1 !C/m^2 in 1D
                            end if
                            exit
                        CASE(2)
                            if (INT(l_f) == NumberXNodes) then
                                v_f(1) = -ABS(v_f(1))
                            else
                                v_f(1) = ABS(v_f(1))
                            end if
                            v_f(2:3) = -v_f(2:3)
                            particleList(j)%refIdx(iThread) = particleList(j)%refIdx(iThread) + 1
                            particleList(j)%refRecordIdx(particleList(j)%refIdx(iThread), iThread) = i - delIdx
                        CASE(3)
                            l_f = ABS(l_f - real(NumberXNodes, kind = real64) - 1.0d0)
                        CASE(4)
                            if (refluxPartBool .or. injectionBool) then
                                particleList(j)%phaseSpace(1, i-delIdx, iThread) = l_f
                                particleList(j)%phaseSpace(2:4, i-delIdx, iThread) = v_f
                                particleList(j)%refIdx(iThread) = particleList(j)%refIdx(iThread) + 1
                                particleList(j)%refRecordIdx(particleList(j)%refIdx(iThread), iThread) = i - delIdx
                            else
                                delIdx = delIdx + 1
                                if (l_f == 1) then
                                    particleList(j)%energyLoss(1, iThread) = particleList(j)%energyLoss(1, iThread) + (SUM(v_f**2))!J/m^2 in 1D
                                    particleList(j)%wallLoss(1, iThread) = particleList(j)%wallLoss(1, iThread) + 1 !C/m^2 in 1D
                                else if (l_f == NumberXNodes) then
                                    particleList(j)%energyLoss(2, iThread) = particleList(j)%energyLoss(2, iThread) + (SUM(v_f**2)) !J/m^2 in 1D
                                    particleList(j)%wallLoss(2, iThread) = particleList(j)%wallLoss(2, iThread) + 1 !C/m^2 in 1D
                                end if
                                exit
                            end if
                        CASE default
                            print *, "l_sub is:", l_sub
                            print *, 'l_f is:', l_f
                            print *, world%boundaryConditions(INT(l_f))
                            print *, "Case does not exist in ongoing substep, depositJ"
                            stop
                        END SELECT
                    else
                        particleList(j)%phaseSpace(1, i-delIdx, iThread) = l_f
                        particleList(j)%phaseSpace(2:4,i-delIdx, iThread) = v_f
                    end if
                    ! now final position/velocity becomes next starting position/velocity
                    l_sub = l_f
                    v_sub = v_f
                end do
                if ((l_f < 1) .or. (l_f > NumberXNodes)) then
                    stop "Have particles travelling outside domain!"
                end if
                
            end do loopParticles
            particleList(j)%N_p(iThread) = particleList(j)%N_p(iThread) - delIdx
            particleList(j)%delIdx(iThread) = delIdx
            !$OMP end parallel
            particleList(j)%accumEnergyLoss = particleList(j)%accumEnergyLoss + SUM(particleList(j)%energyLoss, DIM = 2)
            particleList(j)%accumWallLoss = particleList(j)%accumWallLoss + SUM(particleList(j)%wallLoss, DIM = 2)
        end do loopSpecies
    end subroutine moveParticlesOther


end module mod_particleMover