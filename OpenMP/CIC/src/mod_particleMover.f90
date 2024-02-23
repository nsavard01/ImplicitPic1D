module mod_particleMover
    use iso_fortran_env, only: int32, real64, output_unit
    use constants
    use mod_BasicFunctions
    use mod_particle
    use mod_domain
    use mod_potentialSolver
    use omp_lib
    implicit none
    integer(int32), parameter :: m_Anderson_Particle = 2
    integer(int32), parameter :: maxPartIter = 50


contains

    subroutine getDelTauSubStepNoBField(l_sub, v_sub, del_tau, d_half, q_over_m, l_cell, E_left, E_right, E_x, dx_dl, AtBoundaryBool, FutureAtBoundaryBool, l_boundary)
        integer(int32), intent(in) :: l_cell
        real(real64), intent(in) :: q_over_m, dx_dl, E_left, E_right, l_sub, v_sub
        real(real64), intent(in out) :: del_tau, d_half, E_x
        logical, intent(in) :: AtBoundaryBool
        logical, intent(in out) :: FutureAtBoundaryBool
        integer(int32), intent(in out) :: l_boundary
        integer(int32) :: l_alongV, l_awayV
        real(real64) :: del_tau_temp, a, b, c, l_f_alongV, l_f_awayV
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

    subroutine analyticalParticleMoverNoBField(l_sub, v_sub, l_f, del_tau, d_half, q_over_m, l_cell, E_left, E_right, dx_dl)
        real(real64), intent(in) :: l_sub, v_sub, del_tau, q_over_m, E_left, E_right, dx_dl
        real(real64), intent(in out) :: l_f, d_half
        integer(int32), intent(in) :: l_cell
        real(real64) :: del_tau_sqr, real_l_cell, dx_sqr
        integer(int32) :: k
        del_tau_sqr = del_tau**2
        real_l_cell = real(l_cell, kind = real64)
        dx_sqr = dx_dl**2
        l_f = (2.0d0*E_left*del_tau_sqr*real_l_cell*q_over_m - E_left*del_tau_sqr*l_sub*q_over_m + 2.0d0*E_left*del_tau_sqr*q_over_m - &
            2.0d0*E_right*del_tau_sqr*real_l_cell*q_over_m + E_right*del_tau_sqr*l_sub*q_over_m + 4.0d0*del_tau*dx_dl*v_sub + 4.0d0*dx_sqr*l_sub)&
            /(E_left*del_tau_sqr*q_over_m - E_right*del_tau_sqr*q_over_m + 4.0d0*dx_sqr)
        d_half = (l_f + l_sub)*0.5d0 - real(l_cell)

    
    end subroutine analyticalParticleMoverNoBField

    

    subroutine depositJ(solver, particleList, world, del_t)
        ! particle substepping procedure which deposits J, also used for general moving of particles
        ! boolDepositJ = true : deposit J
        ! boolDepositJ = false: only move particles (done after non-linear iteration), also delete particles for wall collisions (saves extra loop in collisions)
        class(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        type(Particle), intent(in out) :: particleList(:)
        real(real64), intent(in) :: del_t
        !a and c correspond to quadratic equations | l_alongV is nearest integer boundary along velocity component, away is opposite
        real(real64) :: l_f, l_sub, v_sub, v_f, timePassed, del_tau, q_over_m, f_tol, v_half, dx_dl, q_times_wp, d_half, E_x, J_part, J_temp(NumberXNodes+1)
        integer(int32) :: j, i, l_cell, iThread, l_boundary, numIter
        logical :: FutureAtBoundaryBool, AtBoundaryBool, timeNotConvergedBool

        call solver%evaluateEFieldHalfTime()
        f_tol = del_t * 1.d-10
        solver%J = 0
        loopSpecies: do j = 1, numberChargedParticles
            q_over_m = particleList(j)%q/particleList(j)%mass
            q_times_wp = particleList(j)%q * particleList(j)%w_p
            J_temp = 0
            !$OMP parallel private(iThread, i, l_f, l_sub, v_sub, v_f, v_half, timePassed, del_tau, E_x, l_cell, FutureAtBoundaryBool, AtBoundaryBool, &
                dx_dl, l_boundary, d_half, numIter, J_part, timeNotConvergedBool) reduction(+:J_temp)
            iThread = omp_get_thread_num() + 1 
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
                    ! Start AA

                    call getDelTauSubStepNoBField(l_sub, v_sub, del_tau, d_half, q_over_m, l_cell, solver%EField(l_cell), solver%EField(l_cell+1), E_x, dx_dl, AtBoundaryBool, FutureAtBoundaryBool, l_boundary)
                    if (FutureAtBoundaryBool) then
                        l_f = real(l_boundary)
                    else
                        call analyticalParticleMoverNoBField(l_sub, v_sub, l_f, del_tau, d_half, q_over_m, l_cell, solver%EField(l_cell), solver%EField(l_cell+1), dx_dl)
                    end if
                    v_half = (l_f - l_sub) * dx_dl/del_tau
                    v_f = 2.0d0 * v_half - v_sub
                    
        
                    J_part = q_times_wp * (v_half)*del_tau/dx_dl/del_t
                    J_temp(l_cell) = J_temp(l_cell) + J_part * (1.0d0 - d_half)
                    J_temp(l_cell + 1) = J_temp(l_cell + 1) + J_part * (d_half)
                    ! solver%J(l_cell, iThread) = solver%J(l_cell, iThread) + J_part * (1.0d0 - d_half)
                    ! solver%J(l_cell+1, iThread) = solver%J(l_cell+1, iThread) + J_part * (d_half)

                    if (FutureAtBoundaryBool) then
                        SELECT CASE (world%boundaryConditions(l_boundary))
                        CASE(0)
                            l_f = real(l_boundary)
                        CASE(1,4)
                            exit
                        CASE(2)
                            if (l_boundary == NumberXNodes+1) then
                                if (v_f > 0) then
                                    v_f = -v_f
                                end if
                            else
                                if (v_f < 0) then
                                    v_f = -v_f
                                end if
                            end if
                            l_f = real(l_boundary) 
                        CASE(3)
                            l_f = REAL(ABS(l_boundary - real(NumberXNodes+1, kind = real64) - 1))
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
            !$OMP end parallel
            solver%J = solver%J + J_temp
        end do loopSpecies
       
    end subroutine depositJ


    ! -------------------------------------------- Particle mover without boolean checks for depositing J ------------------------------------------------------------

    subroutine moveParticles(solver, particleList, world, del_t)
        ! particle mover to avoid the boolean checks which mostly don't happen when depositing J
        class(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        type(Particle), intent(in out) :: particleList(:)
        real(real64), intent(in) :: del_t
        !a and c correspond to quadratic equations | l_alongV is nearest integer boundary along velocity component, away is opposite
        real(real64) :: l_f, l_sub, v_sub, v_f, timePassed, del_tau, q_over_m, f_tol, d_half, v_half, dx_dl, E_x, energyLoss(2)
        integer(int32) :: j, i, l_cell, iThread, l_boundary, numSubStepAve, numIter, funcEvalCounter, delIdx, refIdx, wallLoss(2)
        logical :: FutureAtBoundaryBool, AtBoundaryBool, refluxedBool, timeNotConvergedBool
        call solver%evaluateEFieldHalfTime()
        f_tol = del_t * 1.d-10
        loopSpecies: do j = 1, numberChargedParticles
            q_over_m = particleList(j)%q/particleList(j)%mass
            numSubStepAve = 0
            funcEvalCounter = 0
            energyLoss = 0.0d0
            wallLoss = 0
            !$OMP parallel private(iThread, i, l_f, l_sub, v_sub, v_f, v_half, timePassed, del_tau, E_x, l_cell, FutureAtBoundaryBool, &
                AtBoundaryBool, dx_dl, l_boundary, d_half, numIter, delIdx, refIdx, refluxedBool, timeNotConvergedBool) reduction(+:numSubStepAve, funcEvalCounter, energyLoss, wallLoss)
            iThread = omp_get_thread_num() + 1 
            delIdx = 0
            refIdx = 0
            loopParticles: do i = 1, particleList(j)%N_p(iThread)
                v_sub = particleList(j)%phaseSpace(2,i,iThread)
                l_sub = particleList(j)%phaseSpace(1,i,iThread)
                timePassed = 0.0d0
                AtBoundaryBool = MOD(l_sub, 1.0d0) == 0.0d0
                refluxedBool = .false.
                timeNotConvergedBool = .true.
                del_tau = del_t
                do while(timeNotConvergedBool)
                    numSubStepAve = numSubStepAve + 1
                    if (.not. AtBoundaryBool) then
                        l_cell = INT(l_sub)
                    else
                        l_cell = INT(l_sub) + (INT(SIGN(1.0d0, v_sub)) - 1)/2
                    end if
                    dx_dl = world%dx_dl(l_cell)
            
                    ! AA particle mover
                    call getDelTauSubStepNoBField(l_sub, v_sub, del_tau, d_half, q_over_m, l_cell, solver%EField(l_cell), solver%EField(l_cell+1), E_x, dx_dl, AtBoundaryBool, FutureAtBoundaryBool, l_boundary)
                    if (FutureAtBoundaryBool) then
                        l_f = real(l_boundary)
                    else
                        call analyticalParticleMoverNoBField(l_sub, v_sub, l_f, del_tau, d_half, q_over_m, l_cell, solver%EField(l_cell), solver%EField(l_cell+1), dx_dl)
                    end if
                    v_half = (l_f - l_sub) * dx_dl/del_tau
                    v_f = 2.0d0 * v_half - v_sub
                
        
                    if (FutureAtBoundaryBool) then
                        SELECT CASE (world%boundaryConditions(l_boundary))
                        CASE(0)
                            l_f = real(l_boundary)
                        CASE(1,4)
                            delIdx = delIdx + 1
                            if (l_boundary == 1) then
                                energyLoss(1) = energyLoss(1) + v_f**2 + (SUM(particleList(j)%phaseSpace(3:4,i,iThread)**2))!J/m^2 in 1D
                                wallLoss(1) = wallLoss(1) + 1 !C/m^2 in 1D
                            else if (l_boundary == NumberXNodes+1) then
                                energyLoss(2) = energyLoss(2) + v_f**2 + (SUM(particleList(j)%phaseSpace(3:4,i,iThread)**2)) !J/m^2 in 1D
                                wallLoss(2) = wallLoss(2) + 1 !C/m^2 in 1D
                            end if
                            exit
                        CASE(2)
                            if (l_boundary == NumberXNodes+1) then
                                if (v_f > 0) then
                                    v_f = -v_f
                                end if
                            else
                                if (v_f < 0) then
                                    v_f = -v_f
                                end if
                            end if
                            l_f = real(l_boundary) 
                            if (.not. refluxedBool) then
                                refIdx = refIdx + 1
                                particleList(j)%refRecordIdx(refIdx, iThread) = i - delIdx
                                refluxedBool = .true.
                            end if
                        CASE(3)
                            l_f = REAL(ABS(l_boundary - real(NumberXNodes+1, kind = real64) - 1))
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
            !$OMP end parallel
            particleList(j)%numSubStepsAve = real(numSubStepAve) / real(SUM(particleList(j)%N_p) + SUM(particleList(j)%delIdx))
            particleList(j)%numFuncEvalAve = real(funcEvalCounter) / real(SUM(particleList(j)%N_p) + SUM(particleList(j)%delIdx))
            particleList(j)%accumEnergyLoss = particleList(j)%accumEnergyLoss + energyLoss
            particleList(j)%accumWallLoss = particleList(j)%accumWallLoss + wallLoss
            particleList(j)%energyLoss = energyLoss
            particleList(j)%wallLoss = wallLoss
        end do loopSpecies
    end subroutine moveParticles

end module mod_particleMover