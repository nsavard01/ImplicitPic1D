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
    integer(int32), parameter, private :: maxPartIter = 50

contains

    subroutine allocateParticleMoverData()

    end subroutine allocateParticleMoverData
    
    subroutine particleSubStepNoBField(l_sub, v_sub, l_f, v_f, del_tau, E_x, q_over_m, dx_dl, l_cell, FutureAtBoundaryBool, l_boundary)
        ! Procedure for a single substep
        real(real64), intent(in) :: q_over_m, dx_dl, l_sub, v_sub, E_x
        real(real64), intent(in out) :: l_f, del_tau, v_f
        logical, intent(in out) :: FutureAtBoundaryBool
        integer(int32), intent(in out) :: l_boundary
        integer(int32), intent(in) :: l_cell
        real(real64) :: diff_PE, a, v_i_sqr
        logical :: inCellBool, equalVSignBool

        ! Solve for particle v_f, l_f over remaining del_tau
        a = q_over_m * E_x
        v_f = v_sub + a * del_tau
        l_f = l_sub + 0.5d0 * (v_sub + v_f) * del_tau / dx_dl
        ! Check if particle is outside cell or flipped direction
        ! This is typically faster for general conditions than checking for substep time to each cell boundary
        ! Faster because many particles will not need these checks
        inCellBool = (INT(l_f) == l_cell)
        equalVSignBool = ((v_sub > 0) .eq. (v_f > 0))
        FutureAtBoundaryBool = (.not. inCellBool) .or. (.not. equalVSignBool)
        if (FutureAtBoundaryBool) then
            ! Get location of boundary it's heading towards initially
            if (v_sub > 0) then
                l_boundary = l_cell + 1
            else
                l_boundary = l_cell
            end if
            diff_PE = 2.0d0 * a * (real(l_boundary) - l_sub) * dx_dl
            v_i_sqr = v_sub**2
            ! Will the the particle be able to reach this boundary?
            FutureAtBoundaryBool = (diff_PE > -v_i_sqr)
            if (FutureAtBoundaryBool) then
                ! calculate v_f and del_tau at boundary along v_sub
                l_f = real(l_boundary)
                v_f = SIGN(1.0d0, v_sub) * SQRT(diff_PE + v_i_sqr)
                del_tau = (v_f - v_sub)/a
            else if (.not. inCellBool) then
                ! calculate v_f and del_tau at boundary in opposite direction
                FutureAtBoundaryBool = .true.
                l_boundary = 2*l_cell + 1 - l_boundary
                l_f = real(l_boundary)
                diff_PE = 2.0d0 * a * (l_f - l_sub) * dx_dl
                v_f = -SIGN(1.0d0, v_sub) * SQRT(diff_PE + v_i_sqr)
                del_tau = (v_f - v_sub)/a
            end if
        end if
       
    end subroutine particleSubStepNoBField



    subroutine depositJ(solver, particleList, world, del_t)
        ! Solve for calculation of J on domain as particle moves along trajectory in current iteration of electric fields
        class(potentialSolver), intent(in out), target :: solver
        type(Domain), intent(in) :: world
        type(Particle), intent(in out) :: particleList(:)
        real(real64), intent(in) :: del_t
        real(real64) :: l_f, l_sub, v_sub, v_f, timePassed, del_tau, f_tol, dx_dl, E_x
        integer(int32) :: j, i, l_cell, iThread, l_boundary, numIter, startTime, endTime
        logical :: FutureAtBoundaryBool
        f_tol = del_t * 1.d-10
        ! Solve for current iteration of E^n+1/2
        call solver%makeHalfTimeEField(particleList(1)%workSpace(1:NumberXHalfNodes,1), world)
        !$OMP parallel private(iThread, i, j, l_f, l_sub, v_sub, v_f, timePassed, del_tau, l_cell, FutureAtBoundaryBool, dx_dl, E_x, l_boundary, &
                numIter)
        iThread = omp_get_thread_num() + 1 
        loopSpecies: do j = 1, numberChargedParticles
            ! Reset workspace for specific thread for adding J to domain
            particleList(j)%workSpace(1:NumberXHalfNodes,iThread) = 0.0d0
            ! loop over particles
            loopParticles: do i = 1, particleList(j)%N_p(iThread)
                ! initialize variables
                v_sub = particleList(j)%phaseSpace(2,i,iThread)
                l_sub = particleList(j)%phaseSpace(1,i,iThread)
                timePassed = 0.0d0
                ! Very small probability particle trajectories end on boundaries, so look along cell to determine
                l_cell = INT(l_sub + SIGN(1.d-12, v_sub))
                del_tau = del_t
                do while(del_tau > f_tol)
                    
                    E_x = solver%EField(l_cell)
                    dx_dl = world%dx_dl(l_cell)

                    call particleSubStepNoBField(l_sub, v_sub, l_f, v_f, del_tau, E_x, particleList(j)%q_over_m, dx_dl, l_cell, FutureAtBoundaryBool, l_boundary)
                    
                    ! Add to J
                    particleList(j)%workSpace(l_cell, iThread) = particleList(j)%workSpace(l_cell, iThread) + l_f - l_sub
                    
                    if (FutureAtBoundaryBool) then
                        ! Determine new cell based on boundary particle is on
                        SELECT CASE (world%boundaryConditions(l_boundary))
                        CASE(0)
                            l_cell = l_cell + INT(SIGN(1.0d0, v_f))
                        CASE(1,4)
                            ! Dirichlet absorbing
                            exit
                        CASE(2)
                            ! Neumann-symmetric
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
                            ! Periodic
                            l_f = real(ABS(l_boundary - NumberXNodes- 1))
                            l_cell = ABS(l_cell - NumberXNodes)
                        END SELECT
                    end if
                    ! Update parameters of particle trajectory for next substep
                    timePassed = timePassed + del_tau
                    del_tau = del_t - timePassed
                    l_sub = l_f
                    v_sub = v_f
                    
                end do
            end do loopParticles
        end do loopSpecies
        !$OMP barrier
        ! Concatenate workspace array to get J array
        particleList(1)%workSpace(world%threadHalfNodeIndx(1,iThread):world%threadHalfNodeIndx(2,iThread), 1) = SUM(particleList(1)%workSpace(world%threadHalfNodeIndx(1,iThread):world%threadHalfNodeIndx(2,iThread), :), DIM=2) * particleList(1)%q_times_wp
        do j = 2, numberChargedParticles
            particleList(1)%workSpace(world%threadHalfNodeIndx(1,iThread):world%threadHalfNodeIndx(2,iThread), 1) = particleList(1)%workSpace(world%threadHalfNodeIndx(1,iThread):world%threadHalfNodeIndx(2,iThread), 1) &
                + SUM(particleList(j)%workSpace(world%threadHalfNodeIndx(1,iThread):world%threadHalfNodeIndx(2,iThread), :), DIM=2) * particleList(j)%q_times_wp
        end do
        !$OMP end parallel
        ! apply smoothing
        if (world%gridSmoothBool) then
            call world%smoothField(particleList(1)%workSpace(1:NumberXHalfNodes,1), solver%J)
            solver%J = solver%J/del_t
        else
            solver%J = particleList(1)%workSpace(1:NumberXHalfNodes,1)/del_t
        end if
    end subroutine depositJ

    subroutine moveParticles(solver, particleList, world, del_t)
        ! Solve for phasespace of particles after convergence of potential solver
        class(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        type(Particle), intent(in out) :: particleList(:)
        real(real64), intent(in) :: del_t
        real(real64) :: l_f, l_sub, v_sub, v_f, timePassed, del_tau, f_tol, dx_dl, E_x
        integer(int32) :: j, i, l_cell, iThread, delIdx, l_boundary, numIter, numSubStepAve(numberChargedParticles), funcEvalCounter(numberChargedParticles), refIdx, N_p
        logical :: refluxedBool, FutureAtBoundaryBool, convergeTimeBool
        f_tol = del_t * 1.d-10
        ! Initialize counters for function evaluations and substep number
        numSubStepAve = 0
        funcEvalCounter = 0
        ! E^1/2 
        call solver%makeHalfTimeEField(particleList(1)%workSpace(1:NumberXHalfNodes,1), world)
        ! Make counters private, otherwise have false sharing which causes significant slowdown! In general, keep variables private whenever possible
        !$OMP parallel private(iThread, i, j, l_f, l_sub, v_sub, v_f, timePassed, del_tau, l_cell, delIdx, dx_dl, E_x, l_boundary, numIter, &
                refIdx, refluxedBool, FutureAtBoundaryBool, N_p, convergeTimeBool) reduction(+:numSubStepAve, funcEvalCounter)
        iThread = omp_get_thread_num() + 1 
        loopSpecies: do j = 1, numberChargedParticles
            ! initialize deletion index at dirichlet boundaries and reflux if applicable to neumann boundaries
            delIdx = 0
            refIdx = 0
            ! initialize loss at walls 
            particleList(j)%wallLoss(:, iThread) = 0
            particleList(j)%energyLoss(:, iThread) = 0.0d0
            N_p = particleList(j)%N_p(iThread)
            loopParticles: do i = 1, N_p
                v_sub = particleList(j)%phaseSpace(2,i,iThread)
                l_sub = particleList(j)%phaseSpace(1,i,iThread)
                timePassed = 0.0d0
                refluxedBool = .false.
                l_cell = INT(l_sub + SIGN(1.d-12, v_sub))
                del_tau = del_t
                convergeTimeBool = .true.
                do while(convergeTimeBool)
                    ! accumulate statistics for num sub steps
                    numSubStepAve(j) = numSubStepAve(j) + 1
                    
                    E_x = solver%EField(l_cell)
                    dx_dl = world%dx_dl(l_cell)
                    
                    call particleSubStepNoBField(l_sub, v_sub, l_f, v_f, del_tau, E_x, particleList(j)%q_over_m, dx_dl, l_cell, FutureAtBoundaryBool, l_boundary)
                    
                    if (FutureAtBoundaryBool) then
                        ! Boundaries similar to deposit J
                        SELECT CASE (world%boundaryConditions(l_boundary))
                        CASE(0)
                            l_cell = l_cell + INT(SIGN(1.0d0, v_f))
                        CASE(1,4)
                            ! delete index up by one, so saved phasespace overwrites current particle number
                            delIdx = delIdx + 1
                            ! accumulate number, energy (v^2), and momentum loss (v) for each particles
                            ! These multiplied by other factors after concatenation
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
                            ! If reflux, save particle index 
                            if (.not. refluxedBool) then
                                refIdx = refIdx + 1
                                particleList(j)%refRecordIdx(refIdx, iThread) = i - delIdx
                                refluxedBool = .true.
                            end if
                        CASE(3)
                            l_f = real(ABS(l_boundary - NumberXNodes- 1))
                            l_cell = ABS(l_cell - NumberXNodes)
                        END SELECT
                    end if
                    
                    timePassed = timePassed + del_tau
                    del_tau = del_t - timePassed
                    ! Save convergeTimeBool to determine if particle has been elminated or not
                    convergeTimeBool = del_tau > f_tol
                    l_sub = l_f
                    v_sub = v_f
                    
                end do
                ! Save final particle properties
                if (.not. convergeTimeBool) then
                    particleList(j)%phaseSpace(1, i-delIdx, iThread) = l_f
                    particleList(j)%phaseSpace(2,i-delIdx, iThread) = v_f
                    particleList(j)%phaseSpace(3:4,i-delIdx, iThread) = particleList(j)%phaseSpace(3:4,i,iThread)
                end if
            end do loopParticles
            ! Per thread number of particles, along with number deleted and number refluxed
            particleList(j)%N_p(iThread) = N_p - delIdx
            particleList(j)%delIdx(iThread) = delIdx
            particleList(j)%refIdx(iThread) = refIdx
        end do loopSpecies
        !$OMP end parallel
        do j = 1, numberChargedParticles
            ! Accumulate statistics
            particleList(j)%numToCollide = particleList(j)%N_p
            particleList(j)%numSubStepsAve = real(numSubStepAve(j)) / real(SUM(particleList(j)%N_p) + SUM(particleList(j)%delIdx))
            particleList(j)%numFuncEvalAve = real(funcEvalCounter(j)) / real(SUM(particleList(j)%N_p) + SUM(particleList(j)%delIdx))
            particleList(j)%accumEnergyLoss = particleList(j)%accumEnergyLoss + SUM(particleList(j)%energyLoss, DIM = 2)
            particleList(j)%accumWallLoss = particleList(j)%accumWallLoss + SUM(particleList(j)%wallLoss, DIM=2)
        end do
    end subroutine moveParticles

end module mod_particleMover