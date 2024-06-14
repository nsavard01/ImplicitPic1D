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
    real(real64), allocatable, private :: particleOmega(:,:), particleAccelGrad(:, :), particleLAccel(:,:), particleAccelLeft(:,:)


contains



    subroutine allocateParticleMoverData()
        allocate(particleOmega(NumberXNodes,numberChargedParticles), particleAccelGrad(NumberXNodes, numberChargedParticles), &
            particleLAccel(NumberXNodes,numberChargedParticles), particleAccelLeft(NumberXNodes,numberChargedParticles))

    end subroutine allocateParticleMoverData

    subroutine initializeParticleMoverData(EField, dx_dl, particleList)
        ! For analytical solver
        real(real64), intent(in) :: EField(NumberXHalfNodes), dx_dl(NumberXNodes)
        type(Particle), intent(in) :: particleList(numberChargedParticles)
        integer(int32) :: i, j
        do j = 1, numberChargedParticles
            do i = 1, NumberXNodes
                particleAccelGrad(i, j) = particleList(j)%q_over_m * (EField(i) - EField(i+1))
                particleAccelLeft(i, j) = particleList(j)%q_over_m * (EField(i))
                particleOmega(i, j) = SQRT(ABS(particleAccelGrad(i, j)))/dx_dl(i)
                particleLAccel(i,j) = -particleAccelLeft(i,j)/particleAccelGrad(i,j)
            end do
        end do

    end subroutine initializeParticleMoverData


    subroutine particleSubStepWell(l_sub, v_sub, del_tau, l_f, v_f, accel_left, accel_grad, omega, l_accel, dx_dl, AtBoundaryBool, FutureAtBoundaryBool, l_boundary, v_sign)
        integer(int32), intent(in) :: v_sign
        real(real64), intent(in) :: dx_dl, l_sub, v_sub, accel_left, accel_grad, omega, l_accel
        real(real64), intent(in out) :: del_tau, l_f, v_f
        logical, intent(in) :: AtBoundaryBool
        logical, intent(in out) :: FutureAtBoundaryBool
        integer(int32), intent(in out) :: l_boundary
        real(real64) :: diff_PE, v_i_sqr, accel_grad_root, u_i, l_sub_sqr, C_1, C_2, phase, cos_phase, sin_phase
        logical :: inCellBool, equalVSignBool, highPhaseBool
        l_sub_sqr = l_sub**2
        v_i_sqr = v_sub**2
        accel_grad_root = omega*dx_dl
        u_i = l_sub + l_accel
        phase = omega*del_tau
        cos_phase = COS(phase)
        sin_phase = SIN(phase)
        l_f = -l_accel + u_i * cos_phase + (v_sub/accel_grad_root) * sin_phase
        v_f = -accel_grad_root * u_i * sin_phase + v_sub * cos_phase
        inCellBool = (INT(l_f+1) == 1)
        equalVSignBool = (v_sign == INT(SIGN(1.0d0, v_f)))
        highPhaseBool = (phase > 2.0d0 * pi)
        FutureAtBoundaryBool = (.not. inCellBool) .or. (INT(SIGN(1.0d0, l_f - l_sub)) /= v_sign) .or. (.not. equalVSignBool) .or. highPhaseBool
        if (FutureAtBoundaryBool) then
            if (v_sign > 0) then
                l_boundary = 1
            else
                l_boundary = 0
            end if
            diff_PE = 2.0d0 * (accel_left * (real(l_boundary) - l_sub) - 0.5d0 * accel_grad * (real(l_boundary)**2 - l_sub_sqr))
            FutureAtBoundaryBool = (diff_PE > -v_i_sqr)
            if (FutureAtBoundaryBool) then
                l_f = real(l_boundary)
                v_f = SQRT(diff_PE + v_i_sqr) * v_sign
                C_2 = u_i**2 * accel_grad + v_i_sqr
                C_1 = u_i * (l_f + l_accel) * accel_grad + v_sub * v_f
                del_tau = ACOS(C_1/C_2)/omega
            else if (highPhaseBool .or. equalVSignBool .or. (.not. inCellBool)) then
                l_boundary = 1-l_boundary
                if (.not. AtBoundaryBool) then
                    diff_PE = 2.0d0 * (accel_left * (real(l_boundary) - l_sub) - 0.5d0 * accel_grad * (real(l_boundary)**2 - l_sub_sqr))
                    FutureAtBoundaryBool = (diff_PE > -v_i_sqr)
                    if (FutureAtBoundaryBool) then
                        l_f = real(l_boundary)
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
                        del_tau = phase/omega
                    end if
                else
                    FutureAtBoundaryBool = .true.
                    l_f = real(l_boundary)
                    v_f = -v_sub
                    phase = v_sub/accel_grad_root/u_i
                    if (phase < 0) then
                        phase = 2.0d0 * (ATAN(phase) + pi)
                    else
                        phase = 2 * ATAN(phase)
                    end if
                    del_tau = phase/omega
                end if
            end if
        end if
    end subroutine particleSubStepWell

    subroutine particleSubStepHill(l_sub, v_sub, del_tau, l_f, v_f, accel_left, accel_grad, omega, l_accel, dx_dl, AtBoundaryBool, FutureAtBoundaryBool, l_boundary, v_sign)
        integer(int32), intent(in) :: v_sign
        real(real64), intent(in) :: dx_dl, l_sub, v_sub, accel_left, accel_grad, omega, l_accel
        real(real64), intent(in out) :: del_tau, l_f, v_f
        logical, intent(in) :: AtBoundaryBool
        logical, intent(in out) :: FutureAtBoundaryBool
        integer(int32), intent(in out) :: l_boundary
        real(real64) :: diff_PE, v_i_sqr, accel_grad_root, u_i, l_sub_sqr, C_1, C_2, phase, cos_phase, sin_phase
        logical :: inCellBool, equalVSignBool

        l_sub_sqr = l_sub**2
        v_i_sqr = v_sub**2
        accel_grad_root = omega*dx_dl
        u_i = l_sub + l_accel
        phase = omega*del_tau
        cos_phase = COSH(phase)
        sin_phase = SINH(phase)
        l_f = -l_accel + u_i * cos_phase + (v_sub/accel_grad_root) * sin_phase
        v_f = accel_grad_root * u_i * sin_phase + v_sub * cos_phase
        inCellBool = (INT(l_f+1) == 1)
        equalVSignBool = (v_sign == INT(SIGN(1.0d0, v_f)))
        FutureAtBoundaryBool = (.not. inCellBool) .or. (.not. equalVSignBool)
        if (FutureAtBoundaryBool) then
            if (v_sign > 0) then
                l_boundary = 1
                FutureAtBoundaryBool = l_f > 1
            else
                l_boundary = 0
                FutureAtBoundaryBool = l_f < 0
            end if
            if (FutureAtBoundaryBool) then
                diff_PE = 2.0d0 * (accel_left * (real(l_boundary) - l_sub) - 0.5d0 * accel_grad * (real(l_boundary)**2 - l_sub_sqr))
            else if ((.not. equalVSignBool) .and. (INT(-l_accel+1) /= 1)) then
                diff_PE = 2.0d0 * (accel_left * (real(l_boundary) - l_sub) - 0.5d0 * accel_grad * (real(l_boundary)**2 - l_sub_sqr))
                FutureAtBoundaryBool = (diff_PE > -v_i_sqr)
            end if
            if (FutureAtBoundaryBool) then
                l_f = real(l_boundary)
                v_f = SQRT(diff_PE + v_i_sqr) * v_sign
                C_1 = accel_grad_root * (v_sub * (l_f + l_accel) - u_i * v_f)
                C_2 = v_i_sqr + u_i**2 * accel_grad
                del_tau = ASINH(C_1/C_2)/omega
            else if  (.not. inCellBool) then
                FutureAtBoundaryBool = .true.
                l_boundary = 1-l_boundary
                l_f = real(l_boundary)
                if (.not. AtBoundaryBool) then
                    diff_PE = 2.0d0 * (accel_left * (l_f - l_sub) - 0.5d0 * accel_grad * (l_f**2 - l_sub_sqr))
                    v_f = -SQRT(diff_PE + v_i_sqr) * v_sign
                    C_1 = accel_grad_root * (v_sub * (l_f + l_accel) - u_i * v_f)
                else
                    v_f = -v_sub
                    C_1 = 2.0d0 * accel_grad_root * v_sub * u_i
                end if
                C_2 = v_i_sqr + u_i**2 * accel_grad
                del_tau = ASINH(C_1/C_2)/omega
            end if
        end if
end subroutine particleSubStepHill

    subroutine particleSubStepNoBField(l_sub, v_sub, del_tau, l_f, v_f, accel_left, accel_grad, omega, l_accel, dx_dl, AtBoundaryBool, FutureAtBoundaryBool, l_boundary, v_sign)
        ! Do initial substep, where particles start between nodes
        integer(int32), intent(in) :: v_sign
        real(real64), intent(in) :: dx_dl, l_sub, v_sub, accel_left, accel_grad, omega, l_accel
        real(real64), intent(in out) :: del_tau, l_f, v_f
        logical, intent(in) :: AtBoundaryBool
        logical, intent(in out) :: FutureAtBoundaryBool
        integer(int32), intent(in out) :: l_boundary
        real(real64) :: diff_PE, v_i_sqr, accel_grad_root, u_i, l_sub_sqr, C_1, C_2, phase, cos_phase, sin_phase
        logical :: inCellBool, equalVSignBool, highPhaseBool
        ! real(real64) :: v_f_test, l_f_test, l_f_prev, dt, v_i_test, l_i_test, l_half, a
        ! integer(int32) :: i, j
        l_sub_sqr = l_sub**2
        v_i_sqr = v_sub**2
        accel_grad_root = omega*dx_dl
        u_i = l_sub + l_accel
        phase = omega*del_tau
        if (accel_grad > 0) then
            ! On potential well
            cos_phase = COS(phase)
            sin_phase = SIN(phase)
            l_f = -l_accel + u_i * cos_phase + (v_sub/accel_grad_root) * sin_phase
            v_f = -accel_grad_root * u_i * sin_phase + v_sub * cos_phase
            inCellBool = (INT(l_f+1) == 1)
            equalVSignBool = (v_sign == INT(SIGN(1.0d0, v_f)))
            highPhaseBool = (phase > 2.0d0 * pi)
            FutureAtBoundaryBool = (.not. inCellBool) .or. (INT(SIGN(1.0d0, l_f - l_sub)) /= v_sign) .or. (.not. equalVSignBool) .or. highPhaseBool
            if (FutureAtBoundaryBool) then
                if (v_sign > 0) then
                    l_boundary = 1
                else
                    l_boundary = 0
                end if
                diff_PE = 2.0d0 * (accel_left * (real(l_boundary) - l_sub) - 0.5d0 * accel_grad * (real(l_boundary)**2 - l_sub_sqr))
                FutureAtBoundaryBool = (diff_PE > -v_i_sqr)
                if (FutureAtBoundaryBool) then
                    l_f = real(l_boundary)
                    v_f = SQRT(diff_PE + v_i_sqr) * v_sign
                    C_2 = u_i**2 * accel_grad + v_i_sqr
                    C_1 = u_i * (l_f + l_accel) * accel_grad + v_sub * v_f
                    del_tau = ACOS(C_1/C_2)/omega
                else if (highPhaseBool .or. equalVSignBool .or. (.not. inCellBool)) then
                    l_boundary = 1-l_boundary
                    if (.not. AtBoundaryBool) then
                        diff_PE = 2.0d0 * (accel_left * (real(l_boundary) - l_sub) - 0.5d0 * accel_grad * (real(l_boundary)**2 - l_sub_sqr))
                        FutureAtBoundaryBool = (diff_PE > -v_i_sqr)
                        if (FutureAtBoundaryBool) then
                            l_f = real(l_boundary)
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
                            del_tau = phase/omega
                        end if
                    else
                        FutureAtBoundaryBool = .true.
                        l_f = real(l_boundary)
                        v_f = -v_sub
                        phase = v_sub/accel_grad_root/u_i
                        if (phase < 0) then
                            phase = 2.0d0 * (ATAN(phase) + pi)
                        else
                            phase = 2 * ATAN(phase)
                        end if
                        del_tau = phase/omega
                    end if
                end if
            end if
        else
            ! On potential hill
            cos_phase = COSH(phase)
            sin_phase = SINH(phase)
            l_f = -l_accel + u_i * cos_phase + (v_sub/accel_grad_root) * sin_phase
            v_f = accel_grad_root * u_i * sin_phase + v_sub * cos_phase
            inCellBool = (INT(l_f+1) == 1)
            equalVSignBool = (v_sign == INT(SIGN(1.0d0, v_f)))
            FutureAtBoundaryBool = (.not. inCellBool) .or. (.not. equalVSignBool)
            if (FutureAtBoundaryBool) then
                if (v_sign > 0) then
                    l_boundary = 1
                    FutureAtBoundaryBool = l_f > 1
                else
                    l_boundary = 0
                    FutureAtBoundaryBool = l_f < 0
                end if
                if (FutureAtBoundaryBool) then
                    diff_PE = 2.0d0 * (accel_left * (real(l_boundary) - l_sub) - 0.5d0 * accel_grad * (real(l_boundary)**2 - l_sub_sqr))
                else if ((.not. equalVSignBool) .and. (INT(-l_accel+1) /= 1)) then
                    diff_PE = 2.0d0 * (accel_left * (real(l_boundary) - l_sub) - 0.5d0 * accel_grad * (real(l_boundary)**2 - l_sub_sqr))
                    FutureAtBoundaryBool = (diff_PE > -v_i_sqr)
                end if
                if (FutureAtBoundaryBool) then
                    l_f = real(l_boundary)
                    v_f = SQRT(diff_PE + v_i_sqr) * v_sign
                    C_1 = accel_grad_root * (v_sub * (l_f + l_accel) - u_i * v_f)
                    C_2 = v_i_sqr + u_i**2 * accel_grad
                    del_tau = ASINH(C_1/C_2)/omega
                else if  (.not. inCellBool) then
                    FutureAtBoundaryBool = .true.
                    l_boundary = 1-l_boundary
                    l_f = real(l_boundary)
                    if (.not. AtBoundaryBool) then
                        diff_PE = 2.0d0 * (accel_left * (l_f - l_sub) - 0.5d0 * accel_grad * (l_f**2 - l_sub_sqr))
                        v_f = -SQRT(diff_PE + v_i_sqr) * v_sign
                        C_1 = accel_grad_root * (v_sub * (l_f + l_accel) - u_i * v_f)
                    else
                        v_f = -v_sub
                        C_1 = 2.0d0 * accel_grad_root * v_sub * u_i
                    end if
                    C_2 = v_i_sqr + u_i**2 * accel_grad
                    del_tau = ASINH(C_1/C_2)/omega
                end if
            end if
        end if

        ! dt = del_tau/real(1000)
        ! l_i_test = l_sub
        ! v_i_test = v_sub
        ! do i=1, 1000
        !     l_f_prev = l_i_test
        !     do j =1,50
        !         l_half = 0.5d0 * (l_i_test + l_f_prev)
        !         a = (accel_left * (1.0d0 - l_half) + (accel_left - accel_grad) * l_half)/dx_dl
        !         v_f_test = dt * a + v_i_test
        !         l_f_test = 0.5d0 * (v_i_test + v_f_test) * dt / dx_dl + l_i_test
        !         if (ABS(l_f_test - l_f_prev) < 1d-8) exit
        !         l_f_prev = l_f_test
        !     end do
        !     v_i_test = v_f_test
        !     l_i_test = l_f_test
        ! end do

        ! if (ABS(l_f_test - l_f) > 1e-3) then
        !     print *, 'l_f error too large'
        !     print *, 'error is:', ABS(l_f_test - l_f)
        !     print *, 'l_sub:', l_sub
        !     print *, 'v_sub:', v_sub
        !     print *, 'accel_grad:', accel_grad
        !     print *, 'del_tau:', del_tau
        !     print *, 'l_f:', l_f
        !     print *, 'v_f:', v_f
        !     print *, 'l_f_test:', l_f_test
        !     print *, 'v_f_test:', v_f_test
        !     print *, 'dt:', dt
        !     stop
        ! else if (ABS((v_f_test- v_f)/v_f) > 1e-3) then
        !     print *, 'v_f error too large'
        !     print *, 'error is:', ABS((v_f_test- v_f)/v_f) 
        !     print *, 'l_sub:', l_sub
        !     print *, 'v_sub:', v_sub
        !     print *, 'accel_grad:', accel_grad
        !     print *, 'del_tau:', del_tau
        !     print *, 'l_f:', l_f
        !     print *, 'v_f:', v_f
        !     print *, 'l_f_test:', l_f_test
        !     print *, 'v_f_test:', v_f_test
        !     print *, 'dt:', dt
        !     stop
        !     stop
        ! end if
    end subroutine particleSubStepNoBField

  


    !-------------------------------------- Analytical particle movers ----------------------------------------------    

    subroutine depositJ(solver, particleList, world, del_t)
        ! particle substepping procedure which deposits J, also used for general moving of particles
        ! boolDepositJ = true : deposit J
        ! boolDepositJ = false: only move particles (done after non-linear iteration), also delete particles for wall collisions (saves extra loop in collisions)
        class(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        type(Particle), intent(in out) :: particleList(:)
        real(real64), intent(in) :: del_t
        !a and c correspond to quadratic equations | l_alongV is nearest integer boundary along velocity component, away is opposite
        real(real64) :: l_f, l_sub, v_sub, v_f, timePassed, del_tau, f_tol, dx_dl, J_part
        integer(int32) :: j, i, l_cell, iThread, l_boundary, numIter, v_sign, leftThreadIndx, rightThreadIndx
        logical :: FutureAtBoundaryBool, AtBoundaryBool, timeNotConvergedBool

        call solver%makeHalfTimeEField(particleList(1)%workSpace(1:NumberXHalfNodes,1), world)
        call initializeParticleMoverData(solver%EField, world%dx_dl, particleList)
        f_tol = del_t * 1.d-10
        !$OMP parallel private(iThread, i, j, l_f, l_sub, v_sub, v_f, timePassed, del_tau, l_cell, FutureAtBoundaryBool, AtBoundaryBool, &
                dx_dl, l_boundary, numIter, timeNotConvergedBool, J_part, v_sign, leftThreadIndx, rightThreadIndx)
        iThread = omp_get_thread_num() + 1 
        loopSpecies: do j = 1, numberChargedParticles 
            particleList(j)%workSpace(1:numberXHalfNodes, iThread) = 0.0d0
            loopParticles: do i = 1, particleList(j)%N_p(iThread)
                v_sub = particleList(j)%phaseSpace(2,i,iThread)
                l_sub = particleList(j)%phaseSpace(1,i,iThread)
                timePassed = 0.0d0
                AtBoundaryBool = MOD(l_sub, 1.0d0) == 0.0d0
                if (.not. AtBoundaryBool) then
                    l_cell = INT(l_sub)
                else
                    l_cell = INT(l_sub) + (INT(SIGN(1.0d0, v_sub)) - 1)/2
                end if
                timeNotConvergedBool = .true.
                del_tau = del_t
                do while(timeNotConvergedBool)
                    v_sign = INT(SIGN(1.0d0, v_sub))

                    dx_dl = world%dx_dl(l_cell)
                    l_sub = l_sub - real(l_cell)
                    
                    ! call particleSubStepNoBField(l_sub, v_sub, del_tau, l_f, v_f, particleAccelLeft(l_cell, j), particleAccelGrad(l_cell, j), &
                    !     particleOmega(l_cell, j), particleLAccel(l_cell, j), dx_dl, AtBoundaryBool, FutureAtBoundaryBool, l_boundary, v_sign)
                    if (particleAccelGrad(l_cell, j) > 0) then
                        call particleSubStepWell(l_sub, v_sub, del_tau, l_f, v_f, particleAccelLeft(l_cell, j), particleAccelGrad(l_cell, j), particleOmega(l_cell, j), particleLAccel(l_cell, j), &
                        dx_dl, AtBoundaryBool, FutureAtBoundaryBool, l_boundary, v_sign)
                    else
                        call particleSubStepHill(l_sub, v_sub, del_tau, l_f, v_f, particleAccelLeft(l_cell, j), particleAccelGrad(l_cell, j), particleOmega(l_cell, j), particleLAccel(l_cell, j), &
                            dx_dl, AtBoundaryBool, FutureAtBoundaryBool, l_boundary, v_sign)
                    end if
                    J_part = 0.5d0 * (l_f**2 - l_sub**2)
                    particleList(j)%workSpace(l_cell, iThread) = particleList(j)%workSpace(l_cell, iThread) + l_f - l_sub - J_part
                    particleList(j)%workSpace(l_cell + 1, iThread) = particleList(j)%workSpace(l_cell + 1, iThread) + J_part
                    ! solver%J(l_cell, iThread) = solver%J(l_cell, iThread) + J_part * (1.0d0 - d_half)
                    ! solver%J(l_cell+1, iThread) = solver%J(l_cell+1, iThread) + J_part * (d_half)
                    l_f = real(l_cell) + l_f
                    if (FutureAtBoundaryBool) then
                        l_boundary = l_cell + l_boundary
                        SELECT CASE (world%boundaryConditions(l_boundary))
                        CASE(0)
                            l_cell = l_cell + INT(SIGN(1.0d0, v_f))
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
                            l_cell = ABS(l_cell - NumberXHalfNodes)
                        ! CASE default
                        !     print *, "l_sub is:", l_sub
                        !     print *, 'l_f is:', l_f
                        !     print *, world%boundaryConditions(INT(l_f))
                        !     print *, "Case does not exist in ongoing substep, depositJ"
                        !     stop
                        END SELECT
                    end if
                    
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
        leftThreadIndx = world%threadHalfNodeIndx(1,iThread)
        rightThreadIndx = world%threadHalfNodeIndx(2,iThread)
        particleList(1)%workSpace(leftThreadIndx:rightThreadIndx, 1) = SUM(particleList(1)%workSpace(leftThreadIndx:rightThreadIndx, :), DIM=2) * particleList(1)%q_times_wp
        do j = 2, numberChargedParticles
            particleList(1)%workSpace(leftThreadIndx:rightThreadIndx, 1) = particleList(1)%workSpace(leftThreadIndx:rightThreadIndx, 1) &
                + SUM(particleList(j)%workSpace(leftThreadIndx:rightThreadIndx, :), DIM=2) * particleList(j)%q_times_wp
        end do
        !$OMP end parallel
        SELECT CASE (world%boundaryConditions(1))
        CASE(1,4)
            particleList(1)%workSpace(1,1) = 2.0d0 * particleList(1)%workSpace(1,1)
        CASE(2)
            particleList(1)%workSpace(1,1) = 0.0d0
        CASE(3)
            particleList(1)%workSpace(1,1) = particleList(1)%workSpace(1,1) + particleList(1)%workSpace(NumberXHalfNodes,1)
        END SELECT
        SELECT CASE (world%boundaryConditions(NumberXHalfNodes))
        CASE(1,4)
            particleList(1)%workSpace(NumberXHalfNodes,1) = 2.0d0 * particleList(1)%workSpace(NumberXHalfNodes,1)
        CASE(2)
            particleList(1)%workSpace(NumberXHalfNodes,1) = 0.0d0
        CASE(3)
            particleList(1)%workSpace(NumberXHalfNodes,1) = particleList(1)%workSpace(1,1)
        END SELECT
        if (world%gridSmoothBool) then
            call world%smoothField(particleList(1)%workSpace(:,1), solver%J)
            solver%J = solver%J/del_t
        else
            solver%J = particleList(1)%workSpace(:,1)/del_t
        end if
    end subroutine depositJ


    subroutine moveParticles(solver, particleList, world, del_t)
        ! particle mover to avoid the boolean checks which mostly don't happen when depositing J
        class(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        type(Particle), intent(in out) :: particleList(:)
        real(real64), intent(in) :: del_t
        !a and c correspond to quadratic equations | l_alongV is nearest integer boundary along velocity component, away is opposite
        real(real64) :: l_f, l_sub, v_sub, v_f, timePassed, del_tau, f_tol, v_half, dx_dl
        integer(int32) :: j, i, l_cell, iThread, l_boundary, numSubStepAve(numberChargedParticles), numIter, funcEvalCounter(numberChargedParticles), delIdx, refIdx, N_p, v_sign
        logical :: FutureAtBoundaryBool, AtBoundaryBool, refluxedBool, timeNotConvergedBool

        call solver%makeHalfTimeEField(particleList(1)%workSpace(1:NumberXHalfNodes,1), world)
        call initializeParticleMoverData(solver%EField, world%dx_dl, particleList)
        f_tol = del_t * 1.d-10
        numSubStepAve = 0
        funcEvalCounter = 0
        !$OMP parallel private(iThread, i, j,l_f, l_sub, v_sub, v_f, v_half, timePassed, del_tau, l_cell, FutureAtBoundaryBool, &
                AtBoundaryBool, dx_dl, l_boundary, numIter, delIdx, refIdx, refluxedBool, timeNotConvergedBool, N_p, v_sign) reduction(+:numSubStepAve, funcEvalCounter)
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
                if (.not. AtBoundaryBool) then
                    l_cell = INT(l_sub)
                else
                    l_cell = INT(l_sub) + (INT(SIGN(1.0d0, v_sub)) - 1)/2
                end if
                refluxedBool = .false.
                timeNotConvergedBool = .true.
                del_tau = del_t
                do while(timeNotConvergedBool)
                    v_sign = INT(SIGN(1.0d0, v_sub))
                    numSubStepAve(j) = numSubStepAve(j) + 1
                    dx_dl = world%dx_dl(l_cell)
                    l_sub = l_sub - real(l_cell)

                    if (particleAccelGrad(l_cell, j) > 0) then
                        call particleSubStepWell(l_sub, v_sub, del_tau, l_f, v_f, particleAccelLeft(l_cell, j), particleAccelGrad(l_cell, j), particleOmega(l_cell, j), particleLAccel(l_cell, j), &
                        dx_dl, AtBoundaryBool, FutureAtBoundaryBool, l_boundary, v_sign)
                    else
                        call particleSubStepHill(l_sub, v_sub, del_tau, l_f, v_f, particleAccelLeft(l_cell, j), particleAccelGrad(l_cell, j), particleOmega(l_cell, j), particleLAccel(l_cell, j), &
                            dx_dl, AtBoundaryBool, FutureAtBoundaryBool, l_boundary, v_sign)
                    end if
                  
                    ! call particleSubStepNoBField(l_sub, v_sub, del_tau, l_f, v_f, particleAccelLeft(l_cell, j), particleAccelGrad(l_cell, j), &
                    !     particleOmega(l_cell, j), particleLAccel(l_cell, j), dx_dl, AtBoundaryBool, FutureAtBoundaryBool, l_boundary, v_sign)
                    l_f = real(l_cell) + l_f
                    if (FutureAtBoundaryBool) then
                        l_boundary = l_cell + l_boundary
                        SELECT CASE (world%boundaryConditions(l_boundary))
                        CASE(0)
                            l_cell = l_cell + INT(SIGN(1.0d0, v_f))
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
                            l_cell = ABS(l_cell - NumberXHalfNodes)
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