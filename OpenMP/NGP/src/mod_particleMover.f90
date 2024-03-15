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

    subroutine particleSubStepPicard(l_sub, v_sub, l_f, v_half, del_tau, E_x, q_over_m, l_cell, dx_dl, FutureAtBoundaryBool, l_boundary)
        integer(int32), intent(in) :: l_cell
        real(real64), intent(in) :: q_over_m, dx_dl, l_sub, v_sub, E_x
        real(real64), intent(in out) :: l_f, v_half, del_tau
        logical, intent(in out) :: FutureAtBoundaryBool
        integer(int32), intent(in out) :: l_boundary
        real(real64) :: c, b, a, del_tau_temp
        
        a = 0.5d0 * q_over_m * E_x
        v_half = v_sub + a * del_tau
        l_f = l_sub + v_half * del_tau / dx_dl
        FutureAtBoundaryBool = (INT(l_f) /= l_cell)
        if (FutureAtBoundaryBool) then
            if (v_half > 0) then
                l_boundary = l_cell + 1
            else
                l_boundary = l_cell
            end if
            l_f = real(l_boundary)
            c = (l_sub - l_f) * dx_dl
            b = v_sub**2 - 4.0d0 * a * c
            if (v_sub * v_half > 0) then
                del_tau_temp= 2.0d0 * ABS(c)/(ABS(v_sub) + SQRT(b))
            else
                if (l_sub /= l_f) then
                    del_tau_temp = 2.0d0 * ABS(c)/(SQRT(b) - ABS(v_sub))
                else
                ! v and a opposite direction, reverses back to initial position
                    del_tau_temp = ABS(v_sub)/ABS(a)
                end if
            end if
            del_tau = del_tau_temp
            v_half = v_sub + a * del_tau
        end if



    end subroutine particleSubStepPicard
    
    subroutine particleSubStepNoBField(l_sub, v_sub, l_f, v_half, del_tau, E_x, q_over_m, l_cell, dx_dl, FutureAtBoundaryBool, l_boundary)
        ! Do initial substep, where particles start between nodes
        integer(int32), intent(in) :: l_cell
        real(real64), intent(in) :: q_over_m, dx_dl, l_sub, v_sub, E_x
        real(real64), intent(in out) :: l_f, v_half, del_tau
        logical, intent(in out) :: FutureAtBoundaryBool
        integer(int32), intent(in out) :: l_boundary
        real(real64) :: c,a, b, del_tau_temp
        integer(int32) :: l_alongV, l_awayV
        ! print *, 'del_tau is:', del_tau
        ! Particle first between nodes, so solve quadratic for that particle depending on conditions
        a = 0.5d0 * q_over_m * E_x
        if (v_sub > 0) then
            l_alongV = l_cell + 1
            l_awayV = l_cell
        else
            l_alongV = l_cell
            l_awayV = l_cell + 1
        end if
        c = (l_sub - real(l_alongV)) * dx_dl
        b = v_sub**2 - 4.0d0*a*c
        if (b > 0) then
            ! v and a opposite direction, but particle can still reach boundary along v
            del_tau_temp = 2.0d0 * ABS(c)/(ABS(v_sub) + SQRT(b))
            FutureAtBoundaryBool = (del_tau_temp < del_tau)
            if (FutureAtBoundaryBool) then
                del_tau = del_tau_temp
                l_boundary = l_alongV
            end if
        else
            ! v and a opposite direction, boundary opposite direction of v
            if (l_sub /= real(l_awayV)) then
                c = (l_sub - real(l_awayV)) * dx_dl
                del_tau_temp = 2.0d0 * ABS(c)/(SQRT(v_sub**2 - 4.0d0*a*c) - ABS(v_sub))
                FutureAtBoundaryBool = (del_tau_temp < del_tau) 
            else
            ! v and a opposite direction, reverses back to initial position
                del_tau_temp = ABS(v_sub)/ABS(a)
                FutureAtBoundaryBool = (del_tau_temp < del_tau)
            end if
            if (FutureAtBoundaryBool) then
                del_tau = del_tau_temp
                l_boundary = l_awayV
            end if
        end if
        ! if (del_tau < 0.0d0) then
        !     print *, 'del_tau < 0'
        !     stop
        ! end if
        v_half = v_sub + a * del_tau
        if (.not. FutureAtBoundaryBool) then
            l_f = l_sub + v_half * del_tau / dx_dl
        else
            l_f = real(l_boundary)
        end if 
        ! if (l_f < real(l_cell) .or. l_f > real(l_cell + 1)) then
        !     print *, 'l_f out of bounds'
        !     stop
        ! end if
       
    end subroutine particleSubStepNoBField



    subroutine depositJ(solver, particleList, world, del_t)
        ! particle substepping procedure which deposits J, also used for general moving of particles
        ! boolDepositJ = true : deposit J
        ! boolDepositJ = false: only move particles (done after non-linear iteration), also delete particles for wall collisions (saves extra loop in collisions)
        class(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        type(Particle), intent(in out) :: particleList(:)
        real(real64), intent(in) :: del_t
        !a and c correspond to quadratic equations | l_alongV is nearest integer boundary along velocity component, away is opposite
        real(real64) :: l_f, l_sub, v_sub, v_f, timePassed, del_tau, q_over_m, f_tol, v_half, dx_dl, E_x, q_times_wp, J_temp(NumberXNodes-1)
        integer(int32) :: j, i, l_cell, iThread, l_boundary, numIter
        logical :: AtBoundaryBool, FutureAtBoundaryBool, timeNotConvergedBool
        solver%EField = 0.5d0 * (solver%phi(1:NumberXNodes-1) + solver%phi_f(1:NumberXNodes-1) - solver%phi(2:NumberXNodes) - solver%phi_f(2:NumberXNodes)) / world%dx_dl
        f_tol = del_t * 1.d-10
        solver%J = 0
        loopSpecies: do j = 1, numberChargedParticles
            q_over_m = particleList(j)%q/particleList(j)%mass
            q_times_wp = particleList(j)%q * particleList(j)%w_p
            J_temp = 0
            !$OMP parallel private(iThread, i, l_f, l_sub, v_sub, v_f, v_half, timePassed, del_tau, l_cell, AtBoundaryBool, FutureAtBoundaryBool, dx_dl, E_x, l_boundary, &
                numIter, timeNotConvergedBool) reduction(+:J_temp)
            iThread = omp_get_thread_num() + 1 
            loopParticles: do i = 1, particleList(j)%N_p(iThread)
                v_sub = particleList(j)%phaseSpace(2,i,iThread)
                l_sub = particleList(j)%phaseSpace(1,i,iThread)
                timePassed = 0.0d0
                AtBoundaryBool = MOD(l_sub, 1.0d0) == 0.0d0
                timeNotConvergedBool = .true.
                do while(timeNotConvergedBool)
                    if (.not. AtBoundaryBool) then
                        l_cell = INT(l_sub)
                    else
                        l_cell = INT(l_sub) + (INT(SIGN(1.0d0, v_sub)) - 1)/2
                    end if
                   
                    E_x = solver%EField(l_cell)
                    dx_dl = world%dx_dl(l_cell)
                    del_tau = del_t - timePassed

                    call particleSubStepNoBField(l_sub, v_sub, l_f, v_half, del_tau, E_x, q_over_m, l_cell, dx_dl, FutureAtBoundaryBool, l_boundary)
    
                    v_f = 2.0d0 * v_half - v_sub
                    
                    J_temp(l_cell) = J_temp(l_cell) + (l_f - l_sub)
                    !solver%J(l_cell, iThread) = solver%J(l_cell, iThread) + q_times_wp * (v_half(1))*del_tau/world%dx_dl(l_cell)/del_t
                    
                    if (FutureAtBoundaryBool) then
                        SELECT CASE (world%boundaryConditions(l_boundary))
                        CASE(0)
                            l_f = real(l_boundary)
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
                            l_f = real(l_boundary) 
                        CASE(3)
                            l_f = REAL(ABS(l_boundary - real(NumberXNodes, kind = real64) - 1))
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
                    timeNotConvergedBool = (del_t - timePassed > f_tol)
                    ! now final position/velocity becomes next starting position/velocity
                    l_sub = l_f
                    v_sub = v_f
                    AtBoundaryBool = FutureAtBoundaryBool
                    
                end do
            end do loopParticles
            !solver%J(:, iThread) = solver%J(:, iThread) + J_temp
            !$OMP end parallel
            solver%J = solver%J + J_temp * q_times_wp
        end do loopSpecies
        solver%J = solver%J/del_t
    end subroutine depositJ


    ! -------------------------------------------- Particle mover without boolean checks for depositing J ------------------------------------------------------------

    subroutine moveParticles(solver, particleList, world, del_t)
        ! particle mover to avoid the boolean checks which mostly don't happen when depositing J
        class(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        type(Particle), intent(in out) :: particleList(:)
        real(real64), intent(in) :: del_t
        !a and c correspond to quadratic equations | l_alongV is nearest integer boundary along velocity component, away is opposite
        real(real64) :: l_f, l_sub, v_sub, v_f, timePassed, del_tau, q_over_m, f_tol, v_half, dx_dl, E_x, energyLoss(2)
        integer(int32) :: j, i, l_cell, iThread, delIdx, l_boundary, numIter, numSubStepAve, funcEvalCounter, refIdx, wallLoss(2)
        logical :: AtBoundaryBool, refluxedBool, FutureAtBoundaryBool, timeNotConvergedBool
        solver%EField = 0.5d0 * (solver%phi(1:NumberXNodes-1) + solver%phi_f(1:NumberXNodes-1) - solver%phi(2:NumberXNodes) - solver%phi_f(2:NumberXNodes)) / world%dx_dl
        f_tol = del_t * 1.d-10
        loopSpecies: do j = 1, numberChargedParticles
            q_over_m = particleList(j)%q/particleList(j)%mass
            numSubStepAve = 0
            funcEvalCounter = 0
            energyLoss = 0.0d0
            wallLoss = 0
            !$OMP parallel private(iThread, i, l_f, l_sub, v_sub, v_f, v_half, timePassed, del_tau, l_cell, AtBoundaryBool, delIdx, dx_dl, E_x, l_boundary, numIter, &
                refIdx, refluxedBool, FutureAtBoundaryBool, timeNotConvergedBool) reduction(+:numSubStepAve, funcEvalCounter, energyLoss, wallLoss)
            iThread = omp_get_thread_num() + 1 
            delIdx = 0
            refIdx = 0
            loopParticles: do i = 1, particleList(j)%N_p(iThread)
                v_sub = particleList(j)%phaseSpace(2,i,iThread)
                l_sub = particleList(j)%phaseSpace(1,i,iThread)
                timePassed = 0.0d0
                refluxedBool = .false.
                AtBoundaryBool = MOD(l_sub, 1.0d0) == 0.0d0
                timeNotConvergedBool = .true.
                do while(timeNotConvergedBool)
                    numSubStepAve = numSubStepAve + 1
                    if (.not. AtBoundaryBool) then
                        l_cell = INT(l_sub)
                    else
                        l_cell = INT(l_sub) + (INT(SIGN(1.0d0, v_sub)) - 1)/2
                    end if
                
                    E_x = solver%EField(l_cell)
                    dx_dl = world%dx_dl(l_cell)
                    del_tau = del_t - timePassed
                    
                    call particleSubStepNoBField(l_sub, v_sub, l_f, v_half, del_tau, E_x, q_over_m, l_cell, dx_dl, FutureAtBoundaryBool, l_boundary)
                    
                    
                    v_f = 2.0d0 * v_half - v_sub
                    
                    if (FutureAtBoundaryBool) then
                        SELECT CASE (world%boundaryConditions(l_boundary))
                        CASE(0)
                            l_f = real(l_boundary)
                        CASE(1,4)
                            delIdx = delIdx + 1
                            if (l_boundary == 1) then
                                energyLoss(1) = energyLoss(1) + v_f**2 + (SUM(particleList(j)%phaseSpace(3:4,i,iThread)**2))
                                wallLoss(1) = wallLoss(1) + 1 !C/m^2 in 1D
                                ! particleList(j)%energyLoss(1, iThread) = particleList(j)%energyLoss(1, iThread) + (SUM(v_f**2))!J/m^2 in 1D
                                ! particleList(j)%wallLoss(1, iThread) = particleList(j)%wallLoss(1, iThread) + 1 !C/m^2 in 1D
                            else if (l_boundary == NumberXNodes) then
                                energyLoss(2) = energyLoss(2) + v_f**2 + (SUM(particleList(j)%phaseSpace(3:4,i,iThread)**2)) !J/m^2 in 1D
                                wallLoss(2) = wallLoss(2) + 1 !C/m^2 in 1D
                                ! particleList(j)%energyLoss(2, iThread) = particleList(j)%energyLoss(2, iThread) + (SUM(v_f**2)) !J/m^2 in 1D
                                ! particleList(j)%wallLoss(2, iThread) = particleList(j)%wallLoss(2, iThread) + 1 !C/m^2 in 1D
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
                            l_f = real(l_boundary) 
                            if (.not. refluxedBool) then
                                refIdx = refIdx + 1
                                particleList(j)%refRecordIdx(refIdx, iThread) = i - delIdx
                                refluxedBool = .true.
                            end if
                        CASE(3)
                            l_f = REAL(ABS(l_boundary - real(NumberXNodes, kind = real64) - 1))
                        ! CASE default
                        !     print *, "l_sub is:", l_sub
                        !     print *, 'l_f is:', l_f
                        !     print *, world%boundaryConditions(INT(l_f))
                        !     print *, "Case does not exist in ongoing substep, depositJ"
                        !     stop
                        END SELECT
                    end if
                    
                    timePassed = timePassed + del_tau
                    timeNotConvergedBool = (del_t - timePassed > f_tol)
                    
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