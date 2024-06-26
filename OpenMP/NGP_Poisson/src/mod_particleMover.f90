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
    integer(int32), parameter, private :: m_Anderson_Particle = 2
    integer(int32), parameter, private :: maxPartIter = 50

contains

    subroutine allocateParticleMoverData()

    end subroutine allocateParticleMoverData
    
    subroutine particleSubStepNoBField(l_sub, v_sub, l_f, v_f, del_tau, E_x, q_over_m, dx_dl, FutureAtBoundaryBool, l_boundary)
        ! Do initial substep, where particles start between nodes
        real(real64), intent(in) :: q_over_m, dx_dl, l_sub, v_sub, E_x
        real(real64), intent(in out) :: l_f, del_tau, v_f
        logical, intent(in out) :: FutureAtBoundaryBool
        integer(int32), intent(in out) :: l_boundary
        real(real64) :: diff_PE, a, v_i_sqr
        !integer(int32) :: l_alongV, l_awayV
        logical :: inCellBool, equalVSignBool
        ! print *, 'del_tau is:', del_tau
        ! Particle first between nodes, so solve quadratic for that particle depending on conditions

        a = q_over_m * E_x
        v_f = v_sub + a * del_tau
        l_f = l_sub + 0.5d0 * (v_sub + v_f) * del_tau / dx_dl
        inCellBool = (INT(l_f+1.0d0) == 1)
        equalVSignBool = ((v_sub > 0) .eq. (v_f > 0))
        FutureAtBoundaryBool = (.not. inCellBool) .or. (.not. equalVSignBool)
        if (FutureAtBoundaryBool) then
            if (v_sub > 0) then
                l_boundary = 1
            else
                l_boundary = 0
            end if
            diff_PE = 2.0d0 * a * (real(l_boundary) - l_sub) * dx_dl
            v_i_sqr = v_sub**2
            FutureAtBoundaryBool = (diff_PE > -v_i_sqr)
            if (FutureAtBoundaryBool) then
                l_f = real(l_boundary)
                v_f = SIGN(1.0d0, v_sub) * SQRT(diff_PE + v_i_sqr)
                del_tau = (v_f - v_sub)/a
            else if (.not. inCellBool) then
                FutureAtBoundaryBool = .true.
                l_boundary = 1-l_boundary
                l_f = real(l_boundary)
                diff_PE = 2.0d0 * a * (l_f - l_sub) * dx_dl
                v_f = -SIGN(1.0d0, v_sub) * SQRT(diff_PE + v_i_sqr)
                del_tau = (v_f - v_sub)/a
            end if
        end if
       
    end subroutine particleSubStepNoBField



    subroutine depositJ(solver, particleList, world, del_t)
        ! particle substepping procedure which deposits J, also used for general moving of particles
        ! boolDepositJ = true : deposit J
        ! boolDepositJ = false: only move particles (done after non-linear iteration), also delete particles for wall collisions (saves extra loop in collisions)
        class(potentialSolver), intent(in out), target :: solver
        type(Domain), intent(in) :: world
        type(Particle), intent(in out) :: particleList(:)
        real(real64), intent(in) :: del_t
        !a and c correspond to quadratic equations | l_alongV is nearest integer boundary along velocity component, away is opposite
        real(real64) :: l_f, l_sub, v_sub, v_f, timePassed, del_tau, f_tol, dx_dl, E_x, d
        integer(int32) :: j, i, l_cell, iThread, l_boundary, numIter, startTime, endTime
        logical :: FutureAtBoundaryBool
        f_tol = del_t * 1.d-10
        call solver%makeHalfTimeEField(particleList(1)%workSpace(1:NumberXHalfNodes,1), world)
        !$OMP parallel private(iThread, i, j, l_f, l_sub, v_sub, v_f, timePassed, del_tau, l_cell, FutureAtBoundaryBool, dx_dl, E_x, l_boundary, &
                numIter, d)
        iThread = omp_get_thread_num() + 1 
        loopSpecies: do j = 1, numberChargedParticles
        
            particleList(j)%workSpace(1:NumberXNodes,iThread) = 0.0d0
            loopParticles: do i = 1, particleList(j)%N_p(iThread)
                v_sub = particleList(j)%phaseSpace(2,i,iThread)
                l_sub = particleList(j)%phaseSpace(1,i,iThread)
                timePassed = 0.0d0
                l_cell = INT(l_sub)
                del_tau = del_t
                do while(del_tau > f_tol)
                    
                    E_x = solver%EField(l_cell)
                    dx_dl = world%dx_dl(l_cell)
                    l_sub = l_sub - real(l_cell)
                    call particleSubStepNoBField(l_sub, v_sub, l_f, v_f, del_tau, E_x, particleList(j)%q_over_m, dx_dl, FutureAtBoundaryBool, l_boundary)
                    
                    !J_temp(l_cell) = J_temp(l_cell) + (l_f - l_sub)
                    !particleList(j)%workSpace(l_cell, iThread) = particleList(j)%workSpace(l_cell, iThread) + l_f - l_sub
                    l_f = l_f + real(l_cell)
                    if (FutureAtBoundaryBool) then
                        l_boundary = l_boundary + l_cell
                        SELECT CASE (world%boundaryConditions(l_boundary))
                        CASE(0)
                            l_cell = l_cell + INT(SIGN(1.0d0, v_f))
                        CASE(1,4)
                            exit
                        CASE(2)
                            if (l_boundary == NumberXNodes) then
                                if (v_f > 0) then
                                    v_f = -v_f
                                end if
                            else
                                if (v_f < 0) then
                                    v_f = -v_f
                                end if
                            end if
                        CASE(3)
                            l_f = real(ABS(l_boundary - NumberXNodes- 1))
                            l_cell = ABS(l_cell - NumberXNodes)
                        ! CASE default
                        !     print *, "l_sub is:", l_sub
                        !     print *, 'l_f is:', l_f
                        !     print *, world%boundaryConditions(INT(l_f))
                        !     print *, "Case does not exist in ongoing substep, depositJ"
                        !     stop
                        END SELECT
                    end if
                    ! if ((l_f < 1) .or. (l_f > NumberXNodes)) then
                    !     print *, 'l_f is:', l_f
                    !     stop "Have particles travelling outside the domain in depositJ!"
                    ! end if
                    timePassed = timePassed + del_tau
                    del_tau = del_t - timePassed
                    ! now final position/velocity becomes next starting position/velocity
                    l_sub = l_f
                    v_sub = v_f
                    
                end do
                if (.not. FutureAtBoundaryBool) then
                    d = l_f - real(l_cell)
                    particleList(j)%workSpace(l_cell, iThread) = particleList(j)%workSpace(l_cell, iThread) + (1.0d0-d)
                    particleList(j)%workSpace(l_cell+1, iThread) = particleList(j)%workSpace(l_cell+1, iThread) + d
                end if
            end do loopParticles
            !J_temp = J_temp + particleList(j)%J_particle(:, iThread)* particleList(j)%q_times_wp
        end do loopSpecies
        !$OMP barrier
        particleList(1)%workSpace(world%threadNodeIndx(1,iThread):world%threadNodeIndx(2,iThread), 1) = SUM(particleList(1)%workSpace(world%threadNodeIndx(1,iThread):world%threadNodeIndx(2,iThread), :), DIM=2) * particleList(1)%q_times_wp
        do j = 2, numberChargedParticles
            particleList(1)%workSpace(world%threadNodeIndx(1,iThread):world%threadNodeIndx(2,iThread), 1) = particleList(1)%workSpace(world%threadNodeIndx(1,iThread):world%threadNodeIndx(2,iThread), 1) &
                + SUM(particleList(j)%workSpace(world%threadNodeIndx(1,iThread):world%threadNodeIndx(2,iThread), :), DIM=2) * particleList(j)%q_times_wp
        end do
        !$OMP end parallel
        SELECT CASE (world%boundaryConditions(1))
        CASE(1,4)
            particleList(1)%workSpace(1, 1) = 0.0d0
        CASE(2)
            particleList(1)%workSpace(1,1) = 2.0d0 * particleList(1)%workSpace(1,1)
        CASE(3)
            particleList(1)%workSpace(1,1) = particleList(1)%workSpace(1,1) + particleList(1)%workSpace(NumberXNodes, 1)
        END SELECT
        SELECT CASE (world%boundaryConditions(NumberXNodes))
        CASE(1,4)
            particleList(1)%workSpace(NumberXNodes, 1) = 0.0d0
        CASE(2)
            particleList(1)%workSpace(NumberXNodes,1) = 2.0d0 * particleList(1)%workSpace(NumberXNodes,1)
        CASE(3)
            particleList(1)%workSpace(NumberXNodes,1) = particleList(1)%workSpace(1,1)
        END SELECT
        if (world%gridSmoothBool) then
            do i = 1, NumberXNodes
                SELECT CASE(world%boundaryConditions(i))
                CASE(0)
                    solver%rho_f(i) = 0.25d0 * (particleList(1)%workSpace(i-1, 1) + 2.0d0 * particleList(1)%workSpace(i, 1) + particleList(1)%workSpace(i+1, 1))
                CASE(1,4)
                    solver%rho_f(i) = 0.0d0
                CASE(2)
                    if (i ==1) then
                        solver%rho_f(i) = 0.25d0 * (2.0d0 * particleList(1)%workSpace(1,1) + 2.0d0 * particleList(1)%workSpace(2,1))
                    else    
                        solver%rho_f(i) = 0.25d0 * (2.0d0 * particleList(1)%workSpace(NumberXNodes,1) + 2.0d0 * particleList(1)%workSpace(NumberXHalfNodes,1))
                    end if
                CASE(3)
                    solver%rho_f(i) = 0.25d0 * (particleList(1)%workSpace(NumberXHalfNodes, 1) + 2.0d0 * particleList(1)%workSpace(1, 1) + particleList(1)%workSpace(2, 1))
                END SELECT
            end do
        else
            solver%rho_f = particleList(1)%workSpace(:, 1)
        end if
    end subroutine depositJ


    ! -------------------------------------------- Particle mover without boolean checks for depositing J ------------------------------------------------------------

    subroutine moveParticles(solver, particleList, world, del_t)
        ! particle mover to avoid the boolean checks which mostly don't happen when depositing J
        class(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        type(Particle), intent(in out) :: particleList(:)
        real(real64), intent(in) :: del_t
        !a and c correspond to quadratic equations | l_alongV is nearest integer boundary along velocity component, away is opposite
        real(real64) :: l_f, l_sub, v_sub, v_f, timePassed, del_tau, f_tol, dx_dl, E_x
        integer(int32) :: j, i, l_cell, iThread, delIdx, l_boundary, numIter, numSubStepAve(numberChargedParticles), funcEvalCounter(numberChargedParticles), refIdx, N_p
        logical :: refluxedBool, FutureAtBoundaryBool
        f_tol = del_t * 1.d-10
        numSubStepAve = 0
        funcEvalCounter = 0
        call solver%makeHalfTimeEField(particleList(1)%workSpace(1:NumberXHalfNodes,1), world)
        !$OMP parallel private(iThread, i, j, l_f, l_sub, v_sub, v_f, timePassed, del_tau, l_cell, delIdx, dx_dl, E_x, l_boundary, numIter, &
                refIdx, refluxedBool, FutureAtBoundaryBool, N_p) reduction(+:numSubStepAve, funcEvalCounter)
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
                refluxedBool = .false.
                l_cell = INT(l_sub)
                del_tau = del_t
                do while(del_tau > f_tol)
                    numSubStepAve(j) = numSubStepAve(j) + 1
                    
                    E_x = solver%EField(l_cell)
                    dx_dl = world%dx_dl(l_cell)
                    l_sub = l_sub - real(l_cell)
                    
                    call particleSubStepNoBField(l_sub, v_sub, l_f, v_f, del_tau, E_x, particleList(j)%q_over_m, dx_dl, FutureAtBoundaryBool, l_boundary)
                    
                    l_f = l_f + real(l_cell)
                    if (FutureAtBoundaryBool) then
                        l_boundary = l_boundary + l_cell
                        SELECT CASE (world%boundaryConditions(l_boundary))
                        CASE(0)
                            l_cell = l_cell + INT(SIGN(1.0d0, v_f))
                        CASE(1,4)
                            delIdx = delIdx + 1
                            if (l_boundary == 1) then
                                particleList(j)%energyLoss(1, iThread) = particleList(j)%energyLoss(1, iThread) + v_f**2 + (SUM(particleList(j)%phaseSpace(3:4,i,iThread)**2))
                                particleList(j)%wallLoss(1, iThread) = particleList(j)%wallLoss(1, iThread) + 1 !C/m^2 in 1D
                                particleList(j)%momentumLoss(1, iThread) = particleList(j)%momentumLoss(1, iThread) + v_f
                            else if (l_boundary == NumberXNodes) then
                                particleList(j)%energyLoss(2, iThread) = particleList(j)%energyLoss(2, iThread) + v_f**2 + (SUM(particleList(j)%phaseSpace(3:4,i,iThread)**2)) !J/m^2 in 1D
                                particleList(j)%wallLoss(2, iThread) = particleList(j)%wallLoss(2, iThread) + 1 !C/m^2 in 1D
                                particleList(j)%momentumLoss(2, iThread) = particleList(j)%momentumLoss(2, iThread) + v_f
                            end if
                            exit
                        CASE(2)
                            if (l_boundary == NumberXNodes) then
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
                            l_f = real(ABS(l_boundary - NumberXNodes- 1))
                            l_cell = ABS(l_cell - NumberXNodes)
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
            particleList(j)%N_p(iThread) = N_p - delIdx
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