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
    real(real64), allocatable, private :: particleTimeStep(:,:)


contains

    ! subroutine particleSubStepPicard(l_sub, v_sub, l_f, v_f, d_half, del_tau, E_left, E_right, dx_dl, q_over_m, f_tol, FutureAtBoundaryBool, l_boundary, numIter)
    !     real(real64), intent(in) :: dx_dl, l_sub, v_sub, f_tol, E_left, E_right, q_over_m
    !     real(real64), intent(in out) :: l_f, del_tau, v_f, d_half
    !     logical, intent(in out) :: FutureAtBoundaryBool
    !     integer(int32), intent(in out) :: l_boundary, numIter
    !     integer(int32) :: i
    !     real(real64) :: accel, del_tau_temp, l_f_prev, v_half
    !     ! get index cell where field and dx_dl is evaluated
    !     l_f_prev = l_sub
    !     del_tau = MIN(del_tau, 0.1d0 * dx_dl / SQRT(ABS(q_over_m * (E_left - E_right))))
    !     do i = 1, maxPartIter
    !         d_half = (l_sub + l_f_prev) * 0.5d0
    !         accel = q_over_m * (E_left * (1.0d0 - d_half) + E_right * d_half)
    !         v_half = v_sub + 0.5d0 * accel * del_tau/dx_dl
    !         l_f = l_sub + v_half * del_tau / dx_dl
    !         if (v_half < 0) then
    !             FutureAtBoundaryBool = (l_f < 0)
    !             if (FutureAtBoundaryBool) l_f = 0.0d0
    !         else
    !             FutureAtBoundaryBool = (l_f > 1)
    !             if (FutureAtBoundaryBool) l_f = 1.0d0
    !         end if
    !         if (ABS(l_f - l_f_prev) < 1.d-10) exit
    !         l_f_prev = l_f
    !     end do
    !     numIter = i
    !     if (.not. FutureAtBoundaryBool) then
    !         v_f = 2.0d0 * v_half - v_sub
    !     else
    !         if (l_f /= l_sub) then
    !             v_f = SIGN(1.0d0, v_half) * SQRT(2.0d0 * accel * (l_f - l_sub) + v_sub**2)
    !             del_tau = 2.0d0 * (l_f - l_sub) * dx_dl / (v_f + v_sub)
    !             numIter = numIter + 1
    !         else
    !             del_tau = -2.0d0 * v_sub * dx_dl/accel
    !             v_f = -v_sub
    !         end if
    !         l_boundary = INT(l_f)
    !     end if
    ! end subroutine particleSubStepPicard
    subroutine particleSubStepPicard(l_sub, v_sub, l_f, v_f, d_half, del_tau, E_left, E_right, dx_dl, q_over_m, f_tol, FutureAtBoundaryBool, l_boundary, numIter)
        real(real64), intent(in) :: dx_dl, l_sub, v_sub, f_tol, E_left, E_right, q_over_m
        real(real64), intent(in out) :: l_f, del_tau, v_f, d_half
        logical, intent(in out) :: FutureAtBoundaryBool
        integer(int32), intent(in out) :: l_boundary, numIter
        integer(int32) :: i
        real(real64) :: accel, del_tau_temp, l_f_prev, temp_var, del_tau_prev, v_f_test, l_turnAround
        logical :: vForwardBool, reverseBool
        ! get index cell where field and dx_dl is evaluated
        l_f_prev = l_sub
        del_tau_temp = del_tau
        vForwardBool = (v_sub>0)
        l_boundary = (INT(SIGN(1.0d0, v_sub)) + 1)/2
        l_turnAround = l_sub + 0.5d0 * v_sub * del_tau/dx_dl
        if (INT(l_turnAround + 1)/=1) then
            del_tau = (real(l_boundary) - l_sub) * dx_dl / v_sub
        end if
        do i = 1, maxPartIter
            d_half = (l_sub + l_f_prev) * 0.5d0
            accel = q_over_m * (E_left * (1.0d0 - d_half) + E_right * d_half)
            v_f = v_sub + accel * del_tau/dx_dl
            reverseBool = (v_f > 0) .ne. vForwardBool
            if (reverseBool) v_f = 0.0d0
            l_f = l_sub + 0.5d0 * (v_f + v_sub) * del_tau / dx_dl
            FutureAtBoundaryBool = ((l_f - real(l_boundary) > 0) .eq. vForwardBool)
            if (FutureAtBoundaryBool) l_f = real(l_boundary)
            if (ABS(l_f - l_f_prev) < 1.d-10) exit
            l_f_prev = l_f
        end do
        numIter = i
        if (.not. reverseBool .and. FutureAtBoundaryBool) then
            v_f = SIGN(1.0d0, v_sub) * SQRT(2.0d0 * accel * (l_f - l_sub) + v_sub**2)
            del_tau = 2.0d0 * (l_f - l_sub) * dx_dl / (v_f + v_sub)
            numIter = numIter + 1
        else if (reverseBool) then
            del_tau_prev = 0.0d0
            do i = 1, maxPartIter
                l_f = l_sub + 0.5d0 * v_sub * del_tau_prev / dx_dl
                d_half = (l_sub + l_f) * 0.5d0
                accel = q_over_m * (E_left * (1.0d0 - d_half) + E_right * d_half)
                del_tau = MAX(-v_sub * dx_dl/accel, 0.0)
                if (ABS(del_tau - del_tau_prev) < f_tol) exit
                del_tau_prev = del_tau
            end do
            v_f = -1.d-15 * SIGN(1.0d0, v_sub)
            if (i >= maxPartIter) then
                print *, 'issue converge reverseBool'
                stop
            end if
            numIter = numIter + i
        end if
        v_f_test = v_sub + accel * del_tau/dx_dl
        l_f_prev = l_sub + 0.5d0 * (v_f_test + v_sub) * del_tau / dx_dl
        if ((l_f - l_sub > 0) .ne. vForwardBool) then
            print *, 'final direction not correct'
            stop
        end if
        if (ABS(v_f_test - v_f) > SQRT(e*T_e/m_e) * 1.d-8) then
            print *, 'issue with v_f'
            stop
        end if
        if (reverseBool .and. (.not. FutureAtBoundaryBool) .and. ABS(v_f) > 1.d-8) then
            print *, 'issue with v_f'
            print *, 'v_f:', v_f
            stop
        end if
        if (ABS(l_f_prev - l_f) > 1.d-8) then
            print *, 'issue with l_f'
            stop
        end if
        if (del_tau <= f_tol .or. del_tau > del_tau_temp) then
            print *, 'del_tau off'
            stop
        end if
        if (del_tau <= f_tol) then
            print *, 'del_tau is super small'
            stop
        end if
        ! if (((v_f > 0) .ne. vForwardBool) .and. ABS(v_f) > 1.d-14) then
        !     print *, 'reversed velocity'
        !     stop
        ! end if
        ! if (l_f < 0 .or. l_f > 1) then
        !     print *, 'l_f outside cell'
        !     stop
        ! end if
        
    end subroutine particleSubStepPicard

    subroutine allocateParticleMoverData()

    end subroutine allocateParticleMoverData


    

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
        integer(int32) :: j, i, l_cell, iThread, l_boundary, numIter, leftThreadIndx, rightThreadIndx
        logical :: FutureAtBoundaryBool

        f_tol = del_t * 1.d-10
        call solver%makeHalfTimeEField(particleList(1)%workSpace(1:NumberXHalfNodes,1), world)
        !$OMP parallel private(iThread, i, j, l_f, l_sub, v_sub, v_f, timePassed, del_tau, l_cell, FutureAtBoundaryBool, &
                dx_dl, l_boundary, numIter, J_part, d_half, leftThreadIndx, rightThreadIndx)
        iThread = omp_get_thread_num() + 1 
        loopSpecies: do j = 1, numberChargedParticles 
            particleList(j)%workSpace(1:numberXHalfNodes, iThread) = 0.0d0
            loopParticles: do i = 1, particleList(j)%N_p(iThread)
                v_sub = particleList(j)%phaseSpace(2,i,iThread)
                l_sub = particleList(j)%phaseSpace(1,i,iThread)
                timePassed = 0.0d0
                l_cell = INT(l_sub)
                del_tau = del_t
                do while(del_tau > f_tol)
                    
                    dx_dl = world%dx_dl(l_cell)
                    l_sub = l_sub - real(l_cell)
                    call particleSubStepPicard(l_sub, v_sub, l_f, v_f, d_half, del_tau, solver%EField(l_cell), solver%EField(l_cell+1), dx_dl, &
                        particleList(j)%q_over_m, f_tol, FutureAtBoundaryBool, l_boundary, numIter)
                    
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
                            l_f = REAL(ABS(l_boundary - NumberXHalfNodes - 1))
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
                    ! now final position/velocity becomes next starting position/velocity
                    l_sub = l_f
                    v_sub = v_f
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
        real(real64) :: l_f, l_sub, v_sub, v_f, timePassed, del_tau, f_tol, d_half, dx_dl
        integer(int32) :: j, i, l_cell, iThread, l_boundary, numSubStepAve(numberChargedParticles), numIter, funcEvalCounter(numberChargedParticles), delIdx, refIdx, N_p
        logical :: FutureAtBoundaryBool, refluxedBool
        f_tol = del_t * 1.d-10
        numSubStepAve = 0
        funcEvalCounter = 0
        call solver%makeHalfTimeEField(particleList(1)%workSpace(1:NumberXHalfNodes,1), world)
        !$OMP parallel private(iThread, i, j,l_f, l_sub, v_sub, v_f, d_half, timePassed, del_tau, l_cell, FutureAtBoundaryBool, &
            dx_dl, l_boundary, numIter, delIdx, refIdx, refluxedBool, N_p) reduction(+:numSubStepAve, funcEvalCounter)
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
                l_cell = INT(l_sub)
                refluxedBool = .false.
                del_tau = del_t
                do while(del_tau > f_tol)
                    numSubStepAve(j) = numSubStepAve(j) + 1
                    
                    dx_dl = world%dx_dl(l_cell)
                    l_sub = l_sub - real(l_cell)
                    call particleSubStepPicard(l_sub, v_sub, l_f, v_f, d_half, del_tau, solver%EField(l_cell), solver%EField(l_cell+1), dx_dl, &
                        particleList(j)%q_over_m, f_tol, FutureAtBoundaryBool, l_boundary, numIter)
                    funcEvalCounter(j) = funcEvalCounter(j) + numIter
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
                                particleList(j)%momentumLoss(1,iThread) = particleList(j)%momentumLoss(1,iThread) + v_f
                            else if (l_boundary == NumberXHalfNodes) then
                                particleList(j)%energyLoss(2,iThread) = particleList(j)%energyLoss(2,iThread) + v_f**2 + (SUM(particleList(j)%phaseSpace(3:4,i,iThread)**2)) !J/m^2 in 1D
                                particleList(j)%wallLoss(2,iThread) = particleList(j)%wallLoss(2,iThread) + 1 !C/m^2 in 1D
                                particleList(j)%momentumLoss(2,iThread) = particleList(j)%momentumLoss(2,iThread) + v_f
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
                            l_f = REAL(ABS(l_boundary - NumberXHalfNodes - 1))
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
                    ! now final position/velocity becomes next starting position/velocity
                    l_sub = l_f
                    v_sub = v_f
                end do
                if (.not. FutureAtBoundaryBool) then
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