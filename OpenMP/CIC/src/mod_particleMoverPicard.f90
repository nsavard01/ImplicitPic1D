module mod_particleMover
    use iso_fortran_env, only: int32, real64, output_unit
    use constants
    use mod_BasicFunctions
    use mod_particle
    use mod_domain
    use mod_potentialSolver
    use omp_lib
    implicit none
    integer(int32), parameter, private :: m_Anderson_Particle = 2
    integer(int32), parameter, private :: maxPartIter = 50
    real(real64), allocatable, private :: particleTimeStep(:,:), particleAccelGrad(:, :), particleAccelLeft(:,:)


contains

    subroutine particleSubStepPicard(l_sub, v_sub, l_f, v_f, d_half, del_tau, accel_grad, accel_left, dx_dl, f_tol, FutureAtBoundaryBool, l_boundary, numIter)
        real(real64), intent(in) :: dx_dl, l_sub, v_sub, f_tol, accel_left, accel_grad
        real(real64), intent(in out) :: l_f, del_tau, v_f, d_half
        logical, intent(in out) :: FutureAtBoundaryBool
        integer(int32), intent(in out) :: l_boundary, numIter
        integer(int32) :: i
        real(real64) :: accel, del_tau_temp, l_f_prev, v_half
        ! get index cell where field and dx_dl is evaluated
        l_f_prev = l_sub
        do i = 1, maxPartIter
            d_half = (l_sub + l_f_prev) * 0.5d0
            accel = (accel_left - accel_grad * d_half)
            v_half = v_sub + 0.5d0 * accel * del_tau/dx_dl
            l_f = l_sub + v_half * del_tau / dx_dl
            if (v_half < 0) then
                FutureAtBoundaryBool = (l_f < 0)
                if (FutureAtBoundaryBool) l_f = 0.0d0
            else
                FutureAtBoundaryBool = (l_f > 1)
                if (FutureAtBoundaryBool) l_f = 1.0d0
            end if
            if (ABS(l_f - l_f_prev) < 1.d-10) exit
            l_f_prev = l_f
        end do
        numIter = i
        if (.not. FutureAtBoundaryBool) then
            v_f = 2.0d0 * v_half - v_sub
        else
            if (l_f /= l_sub) then
                v_f = SIGN(1.0d0, v_half) * SQRT(2.0d0 * accel * (l_f - l_sub) + v_sub**2)
                del_tau = 2.0d0 * (l_f - l_sub) * dx_dl / (v_f + v_sub)
                numIter = numIter + 1
            else
                del_tau = ABS(v_sub)/ABS(accel)
                v_f = -v_sub
            end if
            l_boundary = INT(l_f)
        end if

    end subroutine particleSubStepPicard

    subroutine allocateParticleMoverData()
        allocate(particleTimeStep(NumberXNodes,numberChargedParticles), particleAccelGrad(NumberXNodes, numberChargedParticles), &
        particleAccelLeft(NumberXNodes,numberChargedParticles))

    end subroutine allocateParticleMoverData


    subroutine initializeParticleMoverData(EField, dx_dl, particleList)
        ! For picard solver
        real(real64), intent(in) :: EField(NumberXHalfNodes), dx_dl(NumberXNodes)
        type(Particle), intent(in) :: particleList(numberChargedParticles)
        integer(int32) :: i, j
        do j = 1, numberChargedParticles
            do i = 1, NumberXNodes
                particleAccelGrad(i, j) = particleList(j)%q_over_m * (EField(i) - EField(i+1))
                particleAccelLeft(i, j) = particleList(j)%q_over_m * (EField(i))
                particleTimeStep(i, j) = 0.1d0 / (SQRT(ABS(particleAccelGrad(i, j)))/dx_dl(i))
            end do
        end do

    end subroutine initializeParticleMoverData

    subroutine particleSubStepNoBField(l_sub, v_sub, del_tau, l_f, v_f, accel_left, accel_grad, omega, l_accel, dx_dl, AtBoundaryBool, FutureAtBoundaryBool, l_boundary, v_sign)
        ! Do initial substep, where particles start between nodes
        integer(int32), intent(in) :: v_sign
        real(real64), intent(in) :: dx_dl, l_sub, v_sub, accel_left, accel_grad, omega, l_accel
        real(real64), intent(in out) :: del_tau, l_f, v_f
        logical, intent(in) :: AtBoundaryBool
        logical, intent(in out) :: FutureAtBoundaryBool
        integer(int32), intent(in out) :: l_boundary
        real(real64) :: diff_PE, v_i_sqr, accel_grad_root, u_i, l_sub_sqr, C_1, C_2, del_tau_temp, phase, cos_phase, sin_phase
        ! real(real64) :: l_f_test, v_f_test, del_t_test, l_f_prev, l_i_test, v_i_test, l_half_test
        logical :: peakBool
        ! integer(int32) :: i, k
        ! print *, 'del_tau is:', del_tau
        ! Particle first between nodes, so solve quadratic for that particle depending on conditions
        ! accel_left = E_left * q_over_m
        ! accel_right = E_right * q_over_m
        ! accel_grad = (accel_left - accel_right)
        l_sub_sqr = l_sub**2
        v_i_sqr = v_sub**2
        ! l_accel = -accel_left/accel_grad
        accel_grad_root = omega*dx_dl
        u_i = l_sub + l_accel
        if (accel_grad > 0) then
            ! On potential well
            if (v_sign > 0) then
                l_boundary = 1
            else
                l_boundary = 0
            end if
            l_f = real(l_boundary)
            diff_PE = 2.0d0 * (accel_left * (l_f - l_sub) - 0.5d0 * accel_grad * (l_f**2 - l_sub_sqr))
            FutureAtBoundaryBool = (diff_PE > -v_i_sqr)
            if (FutureAtBoundaryBool) then
                v_f = SQRT(diff_PE + v_i_sqr) * v_sign
                C_2 = u_i**2 * accel_grad + v_i_sqr
                C_1 = u_i * (l_f + l_accel) * accel_grad + v_sub * v_f
                del_tau_temp = ACOS(C_1/C_2)/omega
                FutureAtBoundaryBool = (del_tau_temp < del_tau)
            else
                ! Particle does u-turn
                l_boundary = 1 - l_boundary
                l_f = real(l_boundary)
                if (.not. AtBoundaryBool) then
                    diff_PE = 2.0d0 * (accel_left * (l_f - l_sub) - 0.5d0 * accel_grad * (l_f**2 - l_sub_sqr))
                    FutureAtBoundaryBool = (diff_PE > -v_i_sqr)
                    if (FutureAtBoundaryBool) then
                        v_f = -SQRT(diff_PE + v_i_sqr) * v_sign
                        phase = v_sub/accel_grad_root/u_i
                        if (phase < 0) then
                            phase = 2.0d0 * (ATAN(phase) + pi)
                        else
                            phase = 2 * ATAN(phase)
                        end if
                        C_2 = u_i**2 * accel_grad + v_i_sqr
                        C_1 = u_i * (l_f + l_accel) * accel_grad - v_sub * v_f
                        phase = phase + ACOS(C_1/C_2)
                        del_tau_temp = phase/omega
                        FutureAtBoundaryBool = (del_tau_temp < del_tau)
                    end if
                else
                    v_f = -v_sub
                    phase = v_sub/accel_grad_root/u_i
                    if (phase < 0) then
                        phase = 2.0d0 * (ATAN(phase) + pi)
                    else
                        phase = 2 * ATAN(phase)
                    end if
                    del_tau_temp = phase/omega
                    FutureAtBoundaryBool = (del_tau_temp < del_tau)
                end if
            end if
            if (.not. FutureAtBoundaryBool) then
                phase = omega*del_tau
                cos_phase = COS(phase)
                sin_phase = SIN(phase)
                l_f = -l_accel + u_i * cos_phase + (v_sub/accel_grad_root) * sin_phase
                v_f = -accel_grad_root * u_i * sin_phase + v_sub * cos_phase
            end if
        else
            ! On potential hill
            if (v_sign > 0) then
                l_boundary = 1
                peakBool = (-l_accel > l_sub .and. -l_accel < 1)
            else
                l_boundary = 0
                peakBool = (-l_accel < l_sub .and. -l_accel > 0)
            end if
            l_f = real(l_boundary)
            if (peakBool) then
                diff_PE = 2.0d0 * (accel_left * (-l_accel - l_sub) - 0.5d0 * accel_grad * (l_accel**2 - l_sub_sqr))
            else
                diff_PE = 2.0d0 * (accel_left * (l_f - l_sub) - 0.5d0 * accel_grad * (l_f**2 - l_sub_sqr))
            end if
            FutureAtBoundaryBool = (diff_PE > -v_i_sqr)
            if (FutureAtBoundaryBool) then
                if (peakBool) diff_PE = 2.0d0 * (accel_left * (l_f - l_sub) - 0.5d0 * accel_grad * (l_f**2 - l_sub_sqr))
                v_f = SQRT(diff_PE + v_i_sqr) * v_sign
                C_1 = accel_grad_root * (v_sub * (l_f + l_accel) - u_i * v_f)
                C_2 = v_i_sqr + u_i**2 * accel_grad
                del_tau_temp = ASINH(C_1/C_2)/omega
                FutureAtBoundaryBool = (del_tau_temp < del_tau)
            else
                l_boundary = 1 - l_boundary
                l_f = real(l_boundary)
                if (.not. AtBoundaryBool) then
                    diff_PE = 2.0d0 * (accel_left * (l_f - l_sub) - 0.5d0 * accel_grad * (l_f**2 - l_sub_sqr))
                    v_f = -SQRT(diff_PE + v_i_sqr) * v_sign
                    C_1 = accel_grad_root * (v_sub * (l_f + l_accel) - u_i * v_f)
                else
                    v_f = -v_sub
                    C_1 = 2.0d0 * accel_grad_root * (v_sub * u_i)
                end if
                C_2 = v_i_sqr + u_i**2 * accel_grad
                del_tau_temp = ASINH(C_1/C_2)/omega
                FutureAtBoundaryBool = (del_tau_temp < del_tau)
            end if
            if (.not. FutureAtBoundaryBool) then
                phase = omega*del_tau
                cos_phase = COSH(phase)
                sin_phase = SINH(phase)
                l_f = -l_accel + u_i * cos_phase + (v_sub/accel_grad_root) * sin_phase
                v_f = accel_grad_root * u_i * sin_phase + v_sub * cos_phase
            end if
        end if
        if (FutureAtBoundaryBool) del_tau = del_tau_temp
        ! del_t_test = del_tau/1000
        ! l_i_test = l_sub
        ! v_i_test = v_sub
        ! do i = 1, 1000
        !     l_f_prev = l_i_test
        !     do k = 1, 50
        !         l_half_test = 0.5d0 * (l_i_test + l_f_prev)
        !         v_f_test = del_t_test * (accel_left * (1.0d0 - l_half_test) + accel_right * (l_half_test))/dx_dl + v_i_test
        !         l_f_test = 0.5d0 * (v_f_test + v_i_test) * del_t_test/dx_dl + l_i_test
        !         if (ABS(l_f_test - l_f_prev) < 1.d-8) exit
        !         l_f_prev = l_f_test
        !     end do
        !     if (k > 49) then
        !         print *, 'Picard taking too long!'
        !     end if
        !     v_i_test = v_f_test
        !     l_i_test = l_f_test
        ! end do
        ! if (ABS(l_f_test - l_f) > 1.d-6) then
        !     print *, 'Issue with solver'
        !     print *, 'v_sub:', v_sub
        !     print *, 'l_sub:', l_sub
        !     print *, 'accel_left:', accel_left
        !     print *, 'accel_right:', accel_right
        !     print *, 'AtBoundaryBool:', AtBoundaryBool
        !     print *, 'FutureAtBoundaryBool:', FutureAtBoundaryBool
        !     print *, 'del_tau:', del_tau
        !     print *, 'l_f:', l_f
        !     print *, 'l_f_test:', l_f_test
        !     print *, 'v_f:', v_f
        !     print *, 'v_f_test:', v_f_test
        !     print *, 'l_boundary:', l_boundary
        !     diff_PE = 2.0d0 * (accel_left * (l_f_test - l_sub) - 0.5d0 * accel_grad * (l_f_test**2 - l_sub_sqr))
        !     print *, 'Error in picard is:', (diff_PE + v_i_sqr - v_f_test**2)/(diff_PE + v_i_sqr)
        !     stop
        ! end if


        
    end subroutine particleSubStepNoBField

    subroutine getDelTauSubStepNoBField(l_sub, v_sub, del_tau, d_half, q_over_m, l_cell, E_left, E_right, E_x, dx_dl, AtBoundaryBool, FutureAtBoundaryBool, l_boundary)
        integer(int32), intent(in) :: l_cell
        real(real64), intent(in) :: q_over_m, dx_dl, E_left, E_right, l_sub, v_sub
        real(real64), intent(in out) :: del_tau, d_half, E_x
        logical, intent(in) :: AtBoundaryBool
        logical, intent(in out) :: FutureAtBoundaryBool
        integer(int32), intent(in out) :: l_boundary
        integer(int32) :: l_alongV, l_awayV
        real(real64) :: del_tau_temp, a, b, c
        logical :: goingForwardBool
        ! get index cell where field and dx_dl is evaluated
        goingForwardBool = (v_sub > 0)
        if (goingForwardBool) then
            l_alongV = l_cell + 1
            l_awayV = l_cell
        else
            l_alongV = l_cell
            l_awayV = l_cell + 1
        end if
       
        d_half = (l_sub + real(l_alongV)) * 0.5d0 - real(l_cell)
        E_x = (E_right * d_half + E_left * (1.0d0 - d_half))/dx_dl
        a = 0.5d0 * q_over_m * E_x
        ! Particle first between nodes, so solve quadratic for that particle depending on conditions
        c = (l_sub - real(l_alongV)) * dx_dl
        b = v_sub**2 - 4.0d0*a*c
        FutureAtBoundaryBool = (b >= 0)
        if (FutureAtBoundaryBool) then
            del_tau_temp = 2.0d0 * ABS(c)/(ABS(v_sub) + SQRT(b))
            FutureAtBoundaryBool = (del_tau_temp < del_tau)
            if (FutureAtBoundaryBool) then
                del_tau = del_tau_temp
                l_boundary = l_alongV
            end if
        end if
        if (.not. FutureAtBoundaryBool) then
            ! v and a opposite direction, boundary opposite direction of v
            d_half = (l_sub + real(l_awayV)) * 0.5d0 - real(l_cell)
            E_x = (E_right * d_half + E_left * (1.0d0 - d_half))/dx_dl
            a = 0.5d0 * q_over_m * E_x
            FutureAtBoundaryBool = ((a > 0) .neqv. goingForwardBool)
            if (FutureAtBoundaryBool) then
                if (.not. AtBoundaryBool) then
                    c = (l_sub - real(l_awayV)) * dx_dl
                    del_tau_temp = 2.0d0 * ABS(c)/(SQRT(v_sub**2 - 4.0d0*a*c) - ABS(v_sub))
                else
                    del_tau_temp = ABS(v_sub)/ABS(a)
                end if
                FutureAtBoundaryBool = (del_tau_temp < del_tau)
                if (FutureAtBoundaryBool) then
                    del_tau = del_tau_temp
                    l_boundary = l_awayV
                end if
            end if
        end if
        

    end subroutine getDelTauSubStepNoBField


    subroutine getDelTauSubStepNoBFieldNoUTurn(l_sub, v_sub, del_tau, timePassed, d_half, q_over_m, l_cell, E_left, E_right, E_x, dx_dl, AtBoundaryBool, FutureAtBoundaryBool, l_boundary)
        integer(int32), intent(in) :: l_cell
        real(real64), intent(in) :: q_over_m, dx_dl, E_left, E_right, l_sub
        real(real64), intent(in out) :: del_tau, d_half, E_x, timePassed, v_sub
        logical, intent(in) :: AtBoundaryBool
        logical, intent(in out) :: FutureAtBoundaryBool
        integer(int32), intent(in out) :: l_boundary
        real(real64) :: del_tau_temp, a, b, c, real_l_boundary
        logical :: goingForwardBool

        ! get index cell where field and dx_dl is evaluated
        goingForwardBool = (v_sub > 0)    
        if (goingForwardBool) then
            l_boundary = l_cell + 1
        else
            l_boundary = l_cell
        end if
        real_l_boundary = real(l_boundary)
        d_half = (l_sub + real_l_boundary) * 0.5d0 - real(l_cell)
        E_x = (E_right * d_half + E_left * (1.0d0 - d_half))/dx_dl
        a = 0.5d0 * q_over_m * E_x
        ! Particle first between nodes, so solve quadratic for that particle depending on conditions
        c = (l_sub - real_l_boundary) * dx_dl
        b = v_sub**2 - 4.0d0*a*c
        FutureAtBoundaryBool = (b >= 0)
        if (FutureAtBoundaryBool) then
            del_tau_temp = 2.0d0 * ABS(c)/(ABS(v_sub) + SQRT(b))
            FutureAtBoundaryBool = (del_tau_temp < del_tau)
            if (FutureAtBoundaryBool) del_tau = del_tau_temp
        end if
        if (.not. FutureAtBoundaryBool) then
            ! Check goes back around
            d_half = l_sub - real(l_cell)
            E_x = (E_right * d_half + E_left * (1.0d0 - d_half))/dx_dl
            a = 0.5d0 * q_over_m * E_x
            if ((a > 0) .neqv. goingForwardBool) then
                del_tau_temp = ABS(v_sub)/ABS(a)
                if (del_tau_temp < del_tau) then
                    if (.not. AtBoundaryBool) then
                        ! Now velocity and acceleration is same direction at point, so no more turn around
                        ! Try finding substep toward new boundary
                        timePassed = timePassed + del_tau_temp
                        del_tau = del_tau - del_tau_temp
                        v_sub = -v_sub
                        goingForwardBool = .not. goingForwardBool
                        if (goingForwardBool) then
                            l_boundary = l_cell + 1
                        else
                            l_boundary = l_cell
                        end if
                        real_l_boundary = real(l_boundary)
                        d_half = (l_sub + real_l_boundary) * 0.5d0 - real(l_cell)
                        E_x = (E_right * d_half + E_left * (1.0d0 - d_half))/dx_dl
                        a = 0.5d0 * q_over_m * E_x
                        ! Particle first between nodes, so solve quadratic for that particle depending on conditions
                        c = (l_sub - real_l_boundary) * dx_dl
                        b = v_sub**2 - 4.0d0*a*c
                        FutureAtBoundaryBool = (b >= 0)
                        if (FutureAtBoundaryBool) then
                            del_tau_temp = 2.0d0 * ABS(c)/(ABS(v_sub) + SQRT(b))
                            FutureAtBoundaryBool = (del_tau_temp < del_tau)
                            if (FutureAtBoundaryBool) del_tau = del_tau_temp
                        end if
                    else
                        FutureAtBoundaryBool = .true.
                        l_boundary = INT(l_sub)
                        del_tau = del_tau_temp
                    end if
                end if
            end if
        end if
        

    end subroutine getDelTauSubStepNoBFieldNoUTurn

    ! subroutine analyticalParticleMoverNoBField(l_sub, v_sub, l_f, del_tau, d_half, q_over_m, l_cell, E_left, E_right, dx_dl)
    !     real(real64), intent(in) :: l_sub, v_sub, del_tau, q_over_m, E_left, E_right, dx_dl
    !     real(real64), intent(in out) :: l_f, d_half
    !     integer(int32), intent(in) :: l_cell
    !     real(real64) :: del_tau_sqr, real_l_cell, dx_sqr
    !     del_tau_sqr = del_tau**2
    !     real_l_cell = real(l_cell, kind = real64)
    !     dx_sqr = dx_dl**2
    !     l_f = (2.0d0*E_left*del_tau_sqr*real_l_cell*q_over_m - E_left*del_tau_sqr*l_sub*q_over_m + 2.0d0*E_left*del_tau_sqr*q_over_m - &
    !         2.0d0*E_right*del_tau_sqr*real_l_cell*q_over_m + E_right*del_tau_sqr*l_sub*q_over_m + 4.0d0*del_tau*dx_dl*v_sub + 4.0d0*dx_sqr*l_sub)&
    !         /(E_left*del_tau_sqr*q_over_m - E_right*del_tau_sqr*q_over_m + 4.0d0*dx_sqr)
    !     d_half = (l_f + l_sub)*0.5d0 - real(l_cell)

    
    ! end subroutine analyticalParticleMoverNoBField

    !-------------------------------------------- Picard particle movers -----------------------------------------------

    subroutine depositJ(solver, particleList, world, del_t)
        ! particle substepping procedure which deposits J, also used for general moving of particles
        ! boolDepositJ = true : deposit J
        ! boolDepositJ = false: only move particles (done after non-linear iteration), also delete particles for wall collisions (saves extra loop in collisions)
        class(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        type(Particle), intent(in out) :: particleList(:)
        real(real64), intent(in) :: del_t
        !a and c correspond to quadratic equations | l_alongV is nearest integer boundary along velocity component, away is opposite
        real(real64) :: l_f, l_sub, v_sub, v_f, timePassed, del_tau, f_tol, dx_dl, J_part, d_half
        integer(int32) :: j, i, l_cell, iThread, l_boundary, numIter
        logical :: FutureAtBoundaryBool, AtBoundaryBool, timeNotConvergedBool

        call solver%evaluateEFieldHalfTime(world%boundaryConditions(1))
        call initializeParticleMoverData(solver%EField, world%dx_dl, particleList)
        f_tol = del_t * 1.d-10
        solver%J = 0
        !$OMP parallel private(iThread, i, j, l_f, l_sub, v_sub, v_f, timePassed, del_tau, l_cell, FutureAtBoundaryBool, AtBoundaryBool, &
                dx_dl, l_boundary, numIter, timeNotConvergedBool, J_part, d_half)
        iThread = omp_get_thread_num() + 1 
        loopSpecies: do j = 1, numberChargedParticles 
            particleList(j)%workSpace(1:numberXHalfNodes, iThread) = 0.0d0
            loopParticles: do i = 1, particleList(j)%N_p(iThread)
                v_sub = particleList(j)%phaseSpace(2,i,iThread)
                l_sub = particleList(j)%phaseSpace(1,i,iThread)
                timePassed = 0.0d0
                AtBoundaryBool = MOD(l_sub, 1.0d0) == 0.0d0
                timeNotConvergedBool = .true.
                del_tau = del_t
                do while(timeNotConvergedBool)
                    if (.not. AtBoundaryBool) then
                        l_cell = INT(l_sub)
                    else
                        l_cell = INT(l_sub) + (INT(SIGN(1.0d0, v_sub)) - 1)/2
                    end if
                    
                    dx_dl = world%dx_dl(l_cell)
                    l_sub = l_sub - real(l_cell)
                    del_tau = MIN(del_tau, particleTimeStep(l_cell,j))
                    call particleSubStepPicard(l_sub, v_sub, l_f, v_f, d_half, del_tau, particleAccelGrad(l_cell, j), particleAccelLeft(l_cell, j), dx_dl, f_tol, FutureAtBoundaryBool, l_boundary, numIter)
                    
                    J_part = l_f - l_sub
                    particleList(j)%workSpace(l_cell, iThread) = particleList(j)%workSpace(l_cell, iThread) + (1.0d0-d_half) * J_part
                    particleList(j)%workSpace(l_cell + 1, iThread) = particleList(j)%workSpace(l_cell + 1, iThread) + J_part * d_half
                    
                    ! solver%J(l_cell, iThread) = solver%J(l_cell, iThread) + J_part * (1.0d0 - d_half)
                    ! solver%J(l_cell+1, iThread) = solver%J(l_cell+1, iThread) + J_part * (d_half)
                    l_f = real(l_cell) + l_f
                    if (FutureAtBoundaryBool) then
                        l_boundary = l_cell + l_boundary
                        SELECT CASE (world%boundaryConditions(l_boundary))
                        CASE(0)
                            continue
                        CASE(1,4)
                            exit
                        CASE(2)
                            if (l_boundary == NumberXHalfNodes) then
                                if (v_f > 0) then
                                    v_f = -v_f
                                end if
                            else
                                if (v_f < 0) then
                                    v_f = -v_f
                                end if
                            end if
                        CASE(3)
                            l_f = REAL(ABS(l_boundary - real(NumberXHalfNodes, kind = real64) - 1))
                        ! CASE default
                        !     print *, "l_sub is:", l_sub
                        !     print *, 'l_f is:', l_f
                        !     print *, world%boundaryConditions(INT(l_f))
                        !     print *, "Case does not exist in ongoing substep, depositJ"
                        !     stop
                        END SELECT
                    end if
                    ! if ((l_f < l_cell) .or. (l_f > l_cell+1)) then
                    !     print *, 'l_f is:', l_f
                    !     print *, 'l_sub:', l_sub
                    !     print *, 'del_tau:', del_tau
                    !     print *, 'v_sub:', v_sub
                    !     print *, 'l_cell:', l_cell
                    !     print *, 'AtBoundaryBool:', AtBoundaryBool
                    !     print *, 'FutureAtBoundaryBool:', FutureAtBoundaryBool
                    !     stop "Have particles travelling outside local cell!"
                    ! end if
                    timePassed = timePassed + del_tau
                    del_tau = del_t - timePassed
                    timeNotConvergedBool = (del_tau > f_tol)
                    ! now final position/velocity becomes next starting position/velocity
                    l_sub = l_f
                    v_sub = v_f
                    AtBoundaryBool = FutureAtBoundaryBool
                end do
            end do loopParticles
        end do loopSpecies
        !$OMP barrier
        do j = 1, numberChargedParticles
            call world%addThreadedDomainArray(solver%J, particleList(j)%workSpace, NumberXHalfNodes, NumberXHalfNodes, iThread, particleList(j)%q_times_wp)
        end do
        !$OMP end parallel
        solver%J = solver%J/del_t
    end subroutine depositJ


    subroutine moveParticles(solver, particleList, world, del_t)
        ! particle mover to avoid the boolean checks which mostly don't happen when depositing J
        class(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        type(Particle), intent(in out) :: particleList(:)
        real(real64), intent(in) :: del_t
        !a and c correspond to quadratic equations | l_alongV is nearest integer boundary along velocity component, away is opposite
        real(real64) :: l_f, l_sub, v_sub, v_f, timePassed, del_tau, f_tol, d_half, dx_dl
        integer(int32) :: j, i, l_cell, iThread, l_boundary, numSubStepAve(numberChargedParticles), numIter, funcEvalCounter(numberChargedParticles), delIdx, refIdx, N_p
        logical :: FutureAtBoundaryBool, AtBoundaryBool, refluxedBool, timeNotConvergedBool
        call solver%evaluateEFieldHalfTime(world%boundaryConditions(1))
        call initializeParticleMoverData(solver%EField, world%dx_dl, particleList)
        f_tol = del_t * 1.d-10
        numSubStepAve = 0
        funcEvalCounter = 0
        !$OMP parallel private(iThread, i, j,l_f, l_sub, v_sub, v_f, d_half, timePassed, del_tau, l_cell, FutureAtBoundaryBool, &
                AtBoundaryBool, dx_dl, l_boundary, numIter, delIdx, refIdx, refluxedBool, timeNotConvergedBool, N_p) reduction(+:numSubStepAve, funcEvalCounter)
        iThread = omp_get_thread_num() + 1 
        loopSpecies: do j = 1, numberChargedParticles
            delIdx = 0
            refIdx = 0
            particleList(j)%wallLoss(:, iThread) = 0
            particleList(j)%energyLoss(:, iThread) = 0.0d0
            N_p = particleList(j)%N_p(iThread)
            loopParticles: do i = 1, N_p
                v_sub = particleList(j)%phaseSpace(2,i,iThread)
                l_sub = particleList(j)%phaseSpace(1,i,iThread)
                timePassed = 0.0d0
                AtBoundaryBool = MOD(l_sub, 1.0d0) == 0.0d0
                refluxedBool = .false.
                timeNotConvergedBool = .true.
                del_tau = del_t
                do while(timeNotConvergedBool)
                    numSubStepAve(j) = numSubStepAve(j) + 1
                    if (.not. AtBoundaryBool) then
                        l_cell = INT(l_sub)
                    else
                        l_cell = INT(l_sub) + (INT(SIGN(1.0d0, v_sub)) - 1)/2
                    end if
                    dx_dl = world%dx_dl(l_cell)
                    l_sub = l_sub - real(l_cell)
                    ! AA particle mover
                    ! call getDelTauSubStepNoBField(l_sub, v_sub, del_tau, d_half, particleList(j)%q_over_m, l_cell, solver%EField(l_cell), solver%EField(l_cell+1), E_x, dx_dl, AtBoundaryBool, FutureAtBoundaryBool, l_boundary)
                    ! if (FutureAtBoundaryBool) then
                    !     l_f = real(l_boundary)
                    ! else
                    !     call analyticalParticleMoverNoBField(l_sub, v_sub, l_f, del_tau, d_half, particleList(j)%q_over_m, l_cell, solver%EField(l_cell), solver%EField(l_cell+1), dx_dl)
                    ! end if
                    ! v_half = (l_f - l_sub) * dx_dl/del_tau
                    ! ! del_tau = MIN(del_tau, 0.1 * dx_dl / SQRT(ABS(q_over_m * (solver%EField(l_cell+1) - solver%EField(l_cell)))))
                    ! ! call getDelTauSubStepPicard(l_sub, v_sub, l_f, d_half, v_half, del_tau, solver%EField(l_cell), solver%EField(l_cell+1), E_x, q_over_m, l_cell, dx_dl, FutureAtBoundaryBool, l_boundary, AtBoundaryBool)
                    ! v_f = 2.0d0 * v_half - v_sub
                    del_tau = MIN(del_tau, particleTimeStep(l_cell, j))
                    call particleSubStepPicard(l_sub, v_sub, l_f, v_f, d_half, del_tau, particleAccelGrad(l_cell, j), particleAccelLeft(l_cell, j), dx_dl, f_tol, FutureAtBoundaryBool, l_boundary, numIter)
                    funcEvalCounter(j) = funcEvalCounter(j) + numIter
                    l_f = real(l_cell) + l_f
                    if (FutureAtBoundaryBool) then
                        l_boundary = l_cell + l_boundary
                        SELECT CASE (world%boundaryConditions(l_boundary))
                        CASE(0)
                            continue
                        CASE(1,4)
                            delIdx = delIdx + 1
                            if (l_boundary == 1) then
                                particleList(j)%energyLoss(1,iThread) = particleList(j)%energyLoss(1,iThread) + v_f**2 + (SUM(particleList(j)%phaseSpace(3:4,i,iThread)**2))!J/m^2 in 1D
                                particleList(j)%wallLoss(1,iThread) = particleList(j)%wallLoss(1,iThread) + 1 !C/m^2 in 1D
                            else if (l_boundary == NumberXHalfNodes) then
                                particleList(j)%energyLoss(2,iThread) = particleList(j)%energyLoss(2,iThread) + v_f**2 + (SUM(particleList(j)%phaseSpace(3:4,i,iThread)**2)) !J/m^2 in 1D
                                particleList(j)%wallLoss(2,iThread) = particleList(j)%wallLoss(2,iThread) + 1 !C/m^2 in 1D
                            end if
                            exit
                        CASE(2)
                            if (l_boundary == NumberXHalfNodes) then
                                if (v_f > 0) then
                                    v_f = -v_f
                                end if
                            else
                                if (v_f < 0) then
                                    v_f = -v_f
                                end if
                            end if
                            if (.not. refluxedBool) then
                                refIdx = refIdx + 1
                                particleList(j)%refRecordIdx(refIdx, iThread) = i - delIdx
                                refluxedBool = .true.
                            end if
                        CASE(3)
                            l_f = REAL(ABS(l_boundary - real(NumberXHalfNodes, kind = real64) - 1))
                        ! CASE default
                        !     print *, "l_sub is:", l_sub
                        !     print *, 'l_f is:', l_f
                        !     print *, world%boundaryConditions(INT(l_f))
                        !     print *, "Case does not exist in ongoing substep, depositJ"
                        !     stop
                        END SELECT
                    end if
                    ! if ((l_f < l_cell) .or. (l_f > l_cell+1)) then
                    !     print *, 'l_f is:', l_f
                    !     print *, 'l_sub:', l_sub
                    !     print *, 'del_tau:', del_tau
                    !     print *, 'v_sub:', v_sub
                    !     print *, 'l_cell:', l_cell
                    !     print *, 'AtBoundaryBool:', AtBoundaryBool
                    !     print *, 'FutureAtBoundaryBool:', FutureAtBoundaryBool
                    !     stop "Have particles travelling outside local cell!"
                    ! end if
                    timePassed = timePassed + del_tau
                    del_tau = del_t - timePassed
                    timeNotConvergedBool = (del_tau > f_tol) ! Substraction of bool may not have del_t == timePassed, so just use f_tol
                    ! now final position/velocity becomes next starting position/velocity
                    l_sub = l_f
                    v_sub = v_f
                    AtBoundaryBool = FutureAtBoundaryBool  
                end do
            
                if (.not. timeNotConvergedBool) then
                    particleList(j)%phaseSpace(1, i-delIdx, iThread) = l_f
                    particleList(j)%phaseSpace(2,i-delIdx, iThread) = v_f
                    particleList(j)%phaseSpace(3:4,i-delIdx, iThread) = particleList(j)%phaseSpace(3:4,i,iThread)
                end if
            end do loopParticles
            particleList(j)%N_p(iThread) = particleList(j)%N_p(iThread) - delIdx
            particleList(j)%delIdx(iThread) = delIdx
            particleList(j)%refIdx(iThread) = refIdx
        end do loopSpecies
        !$OMP end parallel
        do j = 1, numberChargedParticles
            particleList(j)%numToCollide = particleList(j)%N_p
            particleList(j)%numSubStepsAve = real(numSubStepAve(j)) / real(SUM(particleList(j)%N_p) + SUM(particleList(j)%delIdx))
            particleList(j)%numFuncEvalAve = real(funcEvalCounter(j)) / real(SUM(particleList(j)%N_p) + SUM(particleList(j)%delIdx))
            particleList(j)%accumEnergyLoss = particleList(j)%accumEnergyLoss + SUM(particleList(j)%energyLoss, DIM = 2)
            particleList(j)%accumWallLoss = particleList(j)%accumWallLoss + SUM(particleList(j)%wallLoss, DIM=2)
        end do
    end subroutine moveParticles

    

end module mod_particleMover