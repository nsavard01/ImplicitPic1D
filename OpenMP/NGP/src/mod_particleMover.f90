module mod_particleMover
    use iso_fortran_env, only: int32, real64, output_unit
    use constants
    use mod_BasicFunctions
    use mod_particle
    use mod_domain
    use mod_potentialSolver
    use omp_lib
    implicit none
    ! Procedures for moving particles and depositing J 

contains

    ! subroutine getl_BoundaryInitial(l_sub, v_sub, l_alongV, l_awayV)
    !     ! get point in l-space on boundary which is away or towards boundary based on velocity direction, when particle between nodes
    !     real(real64), intent(in out) :: l_alongV, l_awayV
    !     real(real64), intent(in) :: l_sub, v_sub
    !     if (v_sub > 0.0) then
    !         l_alongV = real(INT(l_sub) + 1, kind = real64)
    !         l_awayV = real(INT(l_sub), kind = real64)
    !     else
    !         l_alongV = real(INT(l_sub), kind = real64)
    !         l_awayV = real(INT(l_sub) + 1, kind = real64)
    !     end if
    ! end subroutine getl_BoundaryInitial

    subroutine particleSubStepInitialTau(world, l_sub, l_f, v_sub, del_tau, l_alongV, l_awayV, l_cell, a)
        ! Do initial substep, where particles start between nodes
        type(Domain), intent(in) :: world
        integer(int32), intent(in) :: l_cell
        real(real64), intent(in out) :: l_f, del_tau, l_alongV, l_awayV, a
        real(real64), intent(in) :: l_sub, v_sub
        real(real64) :: c !rho_i(NumberXNodes), rho_f(NumberXNodes), gradJ(NumberXNodes-2), test(NumberXNodes-2), testConserv, 
        !integer(int32) :: k
        
        ! Particle first between nodes, so solve quadratic for that particle depending on conditions
        if ((v_sub/=0.0d0)) then ! make first case, since pretty much always likely to be the case (could just not have, assume always field, never have to check)
            c = (l_sub - l_alongV) * world%dx_dl(l_cell)
            if (a*v_sub > 0) then
                ! velocity and acceleration in same direction
                del_tau = 2.0d0 * ABS(c)/(ABS(v_sub) + SQRT(v_sub**2 - 4.0d0*a*c))
                l_f = l_alongV
                if (del_tau <= 0.0d0) then
                    print *, "del_tau is:", del_tau
                    print *, "a is:", a
                    print *, "c is:", c
                    print *, "l_sub is:", l_sub
                    print *, "v_sub is:", v_sub
                    stop "Have issue with del_tau for v,a in same direction"
                end if
            else if (v_sub**2 - 4.0d0*a*c > 0.0d0) then
                ! v and a opposite direction, but particle can still reach boundary along v
                del_tau = 2.0d0 * ABS(c)/(ABS(v_sub) + SQRT(v_sub**2 - 4.0d0*a*c))
                l_f = l_alongV
                if (del_tau <= 0.0d0) then
                    print *, "del_tau is:", del_tau
                    print *, "a is:", a
                    print *, "c is:", c
                    print *, "l_sub is:", l_sub
                    print *, "v_sub is:", v_sub
                    stop "Have issue with del_tau for v,a in opposite direction, but still reach boundary along v, initial tau"
                end if
            else
                ! v and a opposite direction, boundary opposite direction of v
                c = (l_sub - l_awayV) * world%dx_dl(l_cell)
                del_tau = 2.0d0 * ABS(c)/(SQRT(v_sub**2 - 4.0d0*a*c) - ABS(v_sub))
                l_f = l_awayV
                if (del_tau <= 0) then
                    stop "Have issue with del_tau for v,a in opposite direction, boundary opposite v"
                end if
            end if
        else
            ! only in direction of field: USE l_alongV AS BOUNDARY ALONG DIRECTION OF a SINCE VELOCITY = 0!!!
            if (a > 0) then
                l_alongV = real(l_cell + 1, kind = real64)
            else
                l_alongV = real(l_cell, kind = real64)
            end if
            c = (l_sub - l_alongV) * world%dx_dl(l_cell)
            del_tau = SQRT(-c/a)
            l_f = l_alongV
            if (del_tau <= 0.0d0) then
                stop "Have issue with del_tau for v = 0"
            end if
        end if

    end subroutine particleSubStepInitialTau

    subroutine particleSubStepTau(world, l_sub, l_f, v_sub, del_tau, l_alongV, l_cell, a)
        ! Substeps, where particles start at nodes
        type(Domain), intent(in) :: world
        integer(int32), intent(in) :: l_cell
        real(real64), intent(in out) :: l_f, del_tau, l_alongV, a
        real(real64), intent(in) :: l_sub, v_sub
        real(real64) :: c
        !integer(int32) :: k
        ! get index cell where field and dx_dl is evaluated
        ! Particle first between nodes, so solve quadratic for that particle depending on conditions
        c = (l_sub - l_alongV) * world%dx_dl(l_cell)
        if (a*v_sub > 0.0d0) then
            ! velocity and acceleration in same direction
            del_tau = 2.0d0 * ABS(c)/(ABS(v_sub) + SQRT(v_sub**2 - 4.0d0*a*c))
            l_f = l_alongV
            if (del_tau <= 0.0d0) then
                stop "Have issue with del_tau for v,a in same direction, ongoing substep"
            end if
        else if (v_sub**2 - 4.0d0*a*c > 0.0d0) then
            ! v and a opposite direction, but particle can still reach boundary along v
            del_tau = 2.0d0 * ABS(c)/(ABS(v_sub) + SQRT(v_sub**2 - 4.0d0*a*c))
            l_f = l_alongV
            if (del_tau <= 0.0d0) then
                print *, "del_tau is:", del_tau
                print *, "a is:", a
                print *, "c is:", c
                print *, "l_sub is:", l_sub
                print *, "v_sub is:", v_sub
                stop "Have issue with del_tau for v,a in opposite direction, but still reach boundary along v, ongoing substep"
            end if
        else
            ! v and a opposite direction, reverses back to initial position
            del_tau = ABS(v_sub)/ABS(a)
            l_f = l_sub
            if (del_tau <= 0.0d0) then
                stop "Have issue with del_tau for v,a in opposite direction, boundary opposite v"
            end if
        end if

    end subroutine particleSubStepTau

    subroutine depositJ(solver, particleList, world, del_t)
        ! particle substepping procedure which deposits J, also used for general moving of particles
        ! boolDepositJ = true : deposit J
        ! boolDepositJ = false: only move particles (done after non-linear iteration), also delete particles for wall collisions (saves extra loop in collisions)
        class(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        type(Particle), intent(in out) :: particleList(:)
        real(real64), intent(in) :: del_t
        !a and c correspond to quadratic equations | l_alongV is nearest integer boundary along velocity component, away is opposite
        real(real64) :: l_f, l_sub, v_sub, v_f, timePassed, del_tau, l_alongV, l_awayV, a
        integer(int32) :: subStepNum, j, i, l_cell, iThread
        solver%J = 0.0d0
        loopSpecies: do j = 1, numberChargedParticles
            !$OMP parallel private(iThread, i, l_f, l_sub, v_sub, v_f, timePassed, del_tau, l_alongV, l_awayV, a, &
                    subStepNum, l_cell)
            iThread = omp_get_thread_num() + 1 
            loopParticles: do i = 1, particleList(j)%N_p(iThread)
                v_sub = particleList(j)%phaseSpace(2,i,iThread)
                l_sub = particleList(j)%phaseSpace(1,i,iThread)
                timePassed = 0.0d0
                subStepNum = 0
                l_cell = INT(l_sub)
                l_alongV = INT(l_sub) + 0.5d0 + SIGN(0.5d0, v_sub)
                l_awayV = INT(l_sub) + 0.5d0 - SIGN(0.5d0, v_sub)
                a = 0.5d0 * (particleList(j)%q / particleList(j)%mass) * solver%EField(l_cell)
                call particleSubStepInitialTau(world, l_sub, l_f, v_sub, del_tau, l_alongV, l_awayV, l_cell, a)
                if (del_tau >= del_t) then
                    ! Add directly to J with no substep
                    l_f = v_sub * del_t / world%dx_dl(l_cell) + (a/ world%dx_dl(l_cell)) * del_t**2 + l_sub
                    v_f = 2.0d0 * (l_f - l_sub) * world%dx_dl(l_cell) / del_t - v_sub
                    timePassed = del_t
                    if (INT(l_f) /= l_cell) then
                        stop "l_f has crossed boundary when condition says is shouldn't have any substeps"
                    end if
                    solver%J(l_cell, iThread) = solver%J(l_cell, iThread) + particleList(j)%w_p * particleList(j)%q * (v_f + v_sub)/2.0d0/world%dx_dl(l_cell)
                else
                    v_f = 2.0d0 * (l_f - l_sub) * world%dx_dl(l_cell) / del_tau - v_sub
                    timePassed = timePassed + del_tau
                    solver%J(l_cell, iThread) = solver%J(l_cell, iThread) + particleList(j)%w_p * particleList(j)%q * (v_f + v_sub)*del_tau/2.0d0/world%dx_dl(l_cell)/del_t
                    if (MOD(l_f, 1.0d0) /= 0.0d0) then
                        print *, l_f
                        stop "l_f is not integer after subStep"
                    end if
                    firstStep: SELECT CASE (world%boundaryConditions(INT(l_f)))
                    CASE(0)
                        continue
                    CASE(1)
                        timePassed = del_t
                    CASE(2)
                        v_f = -v_f
                    CASE(3)
                        l_f = ABS(l_f - real(NumberXNodes, kind = real64) - 1.0d0)
                    CASE default
                        print *, "Case does not exist in first substep, depositJ"
                        stop
                    END SELECT firstStep
                    ! now final position/velocity becomes next starting position/velocity
                    l_sub = l_f
                    v_sub = v_f
                end if
                do while((timePassed < del_t))
                    l_alongV = l_sub + SIGN(1.0d0, v_sub)
                    l_cell = INT(l_sub + SIGN(0.5d0, v_sub))
                    a = (particleList(j)%q / particleList(j)%mass / 2.0d0) * solver%EField(l_cell)
                    call particleSubStepTau(world, l_sub, l_f, v_sub, del_tau, l_alongV, l_cell, a)
                    if (del_tau >= del_t-timePassed) then
                        ! Add directly to J with no substep
                        l_f = v_sub * (del_t-timePassed) / world%dx_dl(l_cell) + (a/ world%dx_dl(l_cell)) * (del_t-timePassed)**2 + l_sub
                        v_f = 2.0d0 * (l_f - l_sub) * world%dx_dl(l_cell) / (del_t - timePassed) - v_sub
                        if (ABS(l_f - l_sub) >= 1) then
                            print *, "l_sub is:", l_sub
                            print *, "v_sub is:", v_sub
                            print *, "a is:", a
                            print *, "relative error is:", a*(del_t-timePassed)/v_sub
                            print *, "c is:", (l_sub - l_alongV) * world%dx_dl(l_cell)
                            print *, "l_f is:", l_f
                            print *, "v_f is:", v_f
                            print *, "del_tau is:", del_tau
                            print *, "remaining is:", del_t - timePassed
                            print *, "del_tau with other inverse formula:", 2.0d0 * ABS((l_sub - l_alongV) * world%dx_dl(l_cell))/(SQRT(v_sub**2 - 4.0d0*a*(l_sub - l_alongV) * world%dx_dl(l_cell)) + ABS(v_sub))
                            stop "l_f has crossed boundary when condition says is shouldn't have any substeps"
                        end if
                        solver%J(l_cell, iThread) = solver%J(l_cell, iThread) + particleList(j)%w_p * particleList(j)%q * (v_f + v_sub)*(del_t - timePassed)/2.0d0/world%dx_dl(l_cell)/del_t
                        timePassed = del_t
                        
                    else
                        v_f = 2.0d0 * (l_f - l_sub) * world%dx_dl(l_cell) / del_tau - v_sub
                        timePassed = timePassed + del_tau
                        solver%J(l_cell, iThread) = solver%J(l_cell, iThread) + particleList(j)%w_p * particleList(j)%q * (v_f + v_sub)*del_tau/2.0d0/world%dx_dl(l_cell)/del_t
                        if (MOD(l_f, 1.0d0) /= 0.0d0) then
                            print *, l_f
                            stop "l_f is not integer after subStep"
                        end if
                        subStep: SELECT CASE (world%boundaryConditions(INT(l_f)))
                        CASE(0)
                            continue
                        CASE(1)
                            exit
                        CASE(2)
                            v_f = -v_f
                        CASE(3)
                            l_f = ABS(l_f - real(NumberXNodes, kind = real64) - 1.0d0)
                        CASE default
                            print *, "l_sub is:", l_sub
                            print *, 'l_f is:', l_f
                            print *, world%boundaryConditions(INT(l_f))
                            print *, "Case does not exist in ongoing substep, depositJ"
                            stop
                        END SELECT subStep
                        ! now final position/velocity becomes next starting position/velocity
                        l_sub = l_f
                        v_sub = v_f
                    end if
                    subStepNum = subStepNum + 1
                end do
                if ((l_f < 1) .or. (l_f > NumberXNodes)) then
                    stop "Have particles travelling oremainDel_tutside domain!"
                end if
            end do loopParticles
            !$OMP end parallel
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
        real(real64) :: l_f, l_sub, v_sub, v_f, timePassed, del_tau, l_alongV, l_awayV, a
        integer(int32) :: subStepNum, j, i, l_cell, delIdx, iThread
        loopSpecies: do j = 1, numberChargedParticles
            !$OMP parallel private(iThread, i, l_f, l_sub, v_sub, v_f, timePassed, del_tau, l_alongV, l_awayV, a, &
                    subStepNum, delIdx, l_cell)
            iThread = omp_get_thread_num() + 1 
            delIdx = 0
            particleList(j)%refIdx(iThread) = 0
            particleList(j)%energyLoss(:, iThread) = 0.0d0
            particleList(j)%wallLoss(:, iThread) = 0.0d0
            loopParticles: do i = 1, particleList(j)%N_p(iThread)
                v_sub = particleList(j)%phaseSpace(2,i,iThread)
                l_sub = particleList(j)%phaseSpace(1,i,iThread)
                timePassed = 0
                subStepNum = 0
                l_cell = INT(l_sub)
                l_alongV = INT(l_sub) + 0.5d0 + SIGN(0.5d0, v_sub)
                l_awayV = INT(l_sub) + 0.5d0 - SIGN(0.5d0, v_sub)
                a = 0.5d0 * (particleList(j)%q / particleList(j)%mass) * solver%EField(l_cell)
                call particleSubStepInitialTau(world, l_sub, l_f, v_sub, del_tau, l_alongV, l_awayV, l_cell, a)
                if (del_tau >= del_t) then
                    ! Add directly to J with no substep
                    l_f = v_sub * del_t / world%dx_dl(l_cell) + (a/ world%dx_dl(l_cell)) * del_t**2 + l_sub
                    v_f = 2.0d0 * (l_f - l_sub) * world%dx_dl(l_cell) / del_t - v_sub
                    particleList(j)%phaseSpace(1, i-delIdx, iThread) = l_f
                    particleList(j)%phaseSpace(2,i-delIdx, iThread) = v_f
                    particleList(j)%phaseSpace(3:4, i-delIdx, iThread) = particleList(j)%phaseSpace(3:4, i, iThread)
                    timePassed = del_t
                    if (INT(l_f) /= l_cell) then
                        stop "l_f has crossed boundary when condition says is shouldn't have any substeps"
                    end if
                    
                else
                    v_f = 2.0d0 * (l_f - l_sub) * world%dx_dl(l_cell) / del_tau - v_sub
                    timePassed = timePassed + del_tau
                    if (MOD(l_f, 1.0d0) /= 0.0d0) then
                        print *, l_f
                        stop "l_f is not integer after subStep"
                    end if
                    firstStep: SELECT CASE (world%boundaryConditions(INT(l_f)))
                    CASE(0)
                        continue
                    CASE(1)
                        timePassed = del_t
                        delIdx = delIdx + 1
                        if (l_f == 1) then
                            particleList(j)%energyLoss(1, iThread) = particleList(j)%energyLoss(1, iThread) + (v_f**2 + SUM(particleList(j)%phaseSpace(3:4, i, iThread)**2))!J/m^2 in 1D
                            particleList(j)%wallLoss(1, iThread) = particleList(j)%wallLoss(1, iThread) + 1 !C/m^2 in 1D
                        else if (l_f == NumberXNodes) then
                            particleList(j)%energyLoss(2, iThread) = particleList(j)%energyLoss(2, iThread) + (v_f**2 + SUM(particleList(j)%phaseSpace(3:4, i, iThread)**2)) !J/m^2 in 1D
                            particleList(j)%wallLoss(2, iThread) = particleList(j)%wallLoss(2, iThread) + 1 !C/m^2 in 1D
                        end if
                    CASE(2)
                        v_f = -v_f
                        particleList(j)%refIdx(iThread) = particleList(j)%refIdx(iThread) + 1
                        particleList(j)%refRecordIdx(particleList(j)%refIdx(iThread), iThread) = i - delIdx
                    CASE(3)
                        l_f = ABS(l_f - real(NumberXNodes, kind = real64) - 1.0d0)
                    CASE default
                        print *, "Case does not exist in first substep, depositJ"
                        stop
                    END SELECT firstStep
                    ! now final position/velocity becomes next starting position/velocity
                    l_sub = l_f
                    v_sub = v_f
                end if
                do while((timePassed < del_t))
                    l_alongV = l_sub + SIGN(1.0d0, v_sub)
                    l_cell = INT(l_sub + SIGN(0.5d0, v_sub))
                    a = (particleList(j)%q / particleList(j)%mass / 2.0d0) * solver%EField(l_cell)
                    call particleSubStepTau(world, l_sub, l_f, v_sub, del_tau, l_alongV, l_cell, a)
                    if (del_tau >= del_t-timePassed) then
                        ! Add directly to J with no substep
                        l_f = v_sub * (del_t-timePassed) / world%dx_dl(l_cell) + (a/ world%dx_dl(l_cell)) * (del_t-timePassed)**2 + l_sub
                        v_f = 2.0d0 * (l_f - l_sub) * world%dx_dl(l_cell) / (del_t - timePassed) - v_sub
                        if (ABS(l_f - l_sub) >= 1) then
                            stop "l_f has crossed boundary when condition says is shouldn't have any substeps"
                        end if
                        particleList(j)%phaseSpace(1, i-delIdx, iThread) = l_f
                        particleList(j)%phaseSpace(2,i-delIdx, iThread) = v_f
                        particleList(j)%phaseSpace(3:4, i-delIdx, iThread) = particleList(j)%phaseSpace(3:4, i, iThread)
                        timePassed = del_t
                        
                    else
                        v_f = 2.0d0 * (l_f - l_sub) * world%dx_dl(l_cell) / del_tau - v_sub
                        timePassed = timePassed + del_tau
                        if (MOD(l_f, 1.0d0) /= 0.0d0) then
                            print *, l_f
                            stop "l_f is not integer after subStep"
                        end if
                        subStep: SELECT CASE (world%boundaryConditions(INT(l_f)))
                        CASE(0)
                            continue
                        CASE(1)
                            delIdx = delIdx + 1
                            if (l_f == 1) then
                                particleList(j)%energyLoss(1, iThread) = particleList(j)%energyLoss(1, iThread) + (v_f**2 + SUM(particleList(j)%phaseSpace(3:4, i, iThread)**2)) !J/m^2 in 1D
                                particleList(j)%wallLoss(1, iThread) = particleList(j)%wallLoss(1, iThread) + 1 !C/m^2 in 1D
                            else if (l_f == NumberXNodes) then
                                particleList(j)%energyLoss(2, iThread) = particleList(j)%energyLoss(2, iThread) + (v_f**2 + SUM(particleList(j)%phaseSpace(3:4, i, iThread)**2))!J/m^2 in 1D
                                particleList(j)%wallLoss(2, iThread) = particleList(j)%wallLoss(2, iThread) + 1 !C/m^2 in 1D
                            end if
                            exit
                        CASE(2)
                            v_f = -v_f
                            particleList(j)%refIdx(iThread) = particleList(j)%refIdx(iThread) + 1
                            particleList(j)%refRecordIdx(particleList(j)%refIdx(iThread), iThread) = i - delIdx
                        CASE(3)
                            l_f = ABS(l_f - real(NumberXNodes, kind = real64) - 1.0d0)
                        CASE default
                            print *, "Case does not exist in ongoing substep, depositJ"
                            stop
                        END SELECT subStep
                        ! now final position/velocity becomes next starting position/velocity
                        l_sub = l_f
                        v_sub = v_f
                    end if
                    subStepNum = subStepNum + 1
                end do
                if ((l_f < 1) .or. (l_f > NumberXNodes)) then
                    stop "Have particles travelling oremainDel_tutside domain!"
                end if
                
            end do loopParticles
            particleList(j)%N_p(iThread) = particleList(j)%N_p(iThread) - delIdx
            particleList(j)%delIdx(iThread) = delIdx
            !$OMP end parallel
            particleList(j)%accumEnergyLoss = particleList(j)%accumEnergyLoss + SUM(particleList(j)%energyLoss, DIM = 2)
            particleList(j)%accumWallLoss = particleList(j)%accumWallLoss + SUM(particleList(j)%wallLoss, DIM = 2)
        end do loopSpecies
    end subroutine moveParticles




end module mod_particleMover