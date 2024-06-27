module mod_particleMover
    use iso_fortran_env, only: int32, real64, output_unit
    use constants
    use mod_BasicFunctions
    use mod_particle
    use mod_domain
    use mod_potentialSolver
    use omp_lib
    implicit none
    integer(int32), parameter, private :: maxPartIter = 50
    real(real64), allocatable, private :: particleTimeStep(:,:)


contains

    subroutine particleSubStepPicard(l_sub, v_sub, l_f, v_f, d_half, del_tau, E_left, E_right, dx_dl, q_over_m, l_cell, FutureAtBoundaryBool, l_boundary, numIter)
        ! Substeps using picard iterations
        real(real64), intent(in) :: dx_dl, l_sub, v_sub, E_left, E_right, q_over_m
        real(real64), intent(in out) :: l_f, del_tau, v_f, d_half
        logical, intent(in out) :: FutureAtBoundaryBool
        integer(int32), intent(in out) :: l_boundary, numIter
        integer(int32), intent(in) :: l_cell
        integer(int32) :: i
        real(real64) :: accel, l_f_prev, v_half
        ! Make sure l_sub is entered normalized to cell values (l_sub - l_cell)
        l_f_prev = l_sub
        do i = 1, maxPartIter
            ! Use previous estimate of l_f to find new l_f
            d_half = (l_sub + l_f_prev) * 0.5d0
            accel = q_over_m * (E_left * (1.0d0 - d_half) + E_right * d_half)
            v_half = v_sub + 0.5d0 * accel * del_tau/dx_dl
            l_f = l_sub + v_half * del_tau / dx_dl
            FutureAtBoundaryBool = ((l_f <= 0) .or. (l_f >= 1))
            if (FutureAtBoundaryBool) then
                ! If particle is outside (or on boundary of) cell, then set particle position to the cell
                if (v_half < 0) then
                    l_f = 0.0d0
                else
                    l_f = 1.0d0
                end if
            end if
            ! Convergence criteria is 1.d-10 within cell, should be quite strict
            if (ABS(l_f - l_f_prev) < 1.d-10) exit
            l_f_prev = l_f
        end do
        numIter = i
        if (.not. FutureAtBoundaryBool) then
            v_f = 2.0d0 * v_half - v_sub
        else
            if (l_f /= l_sub) then
                v_f = SIGN(1.0d0, v_half) * SQRT(2.0d0 * accel * (l_f - l_sub) + v_sub**2)
                del_tau = (l_f - l_sub) * dx_dl / v_half
                numIter = numIter + 1
            else
                ! In the case that the particle goes back to starting position on boundary
                del_tau = -2.0d0 * v_sub * dx_dl/accel
                v_f = -v_sub
            end if
            l_boundary = INT(l_f) + l_cell
        end if
        ! Need to add cell number to l_f
        l_f = l_f + real(l_cell)
    end subroutine particleSubStepPicard

    subroutine allocateParticleMoverData()
        ! Generate grid data for each cell of minimum time step for each particle
        allocate(particleTimeStep(NumberXNodes,numberChargedParticles))

    end subroutine allocateParticleMoverData


    

    !-------------------------------------------- Picard particle movers -----------------------------------------------

    subroutine depositJ(solver, particleList, world, del_t)
        ! Solve for calculation of J on domain as particle moves along trajectory in current iteration of electric fields
        class(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        type(Particle), intent(in out) :: particleList(:)
        real(real64), intent(in) :: del_t
        !a and c correspond to quadratic equations | l_alongV is nearest integer boundary along velocity component, away is opposite
        real(real64) :: l_f, l_sub, v_sub, v_f, timePassed, del_tau, f_tol, dx_dl, J_part, d_half
        integer(int32) :: j, i, l_cell, iThread, l_boundary, numIter, leftThreadIndx, rightThreadIndx
        logical :: FutureAtBoundaryBool

        f_tol = del_t * 1.d-10
        ! Solve for current iteration of E^n+1/2
        call solver%makeHalfTimeEField(particleList(1)%workSpace(1:NumberXHalfNodes,1), world)
        ! Generate minimum substep for each cell for each particle based on gradient of acceleration along the cell
        do j = 1, numberChargedParticles
            particleTimeStep(:, j) = 0.1d0 * world%dx_dl / SQRT(ABS(particleList(j)%q_over_m * (solver%EField(1:NumberXNodes) - solver%EField(2:NumberXHalfNodes))))
        end do
        !$OMP parallel private(iThread, i, j, l_f, l_sub, v_sub, v_f, timePassed, del_tau, l_cell, FutureAtBoundaryBool, &
                dx_dl, l_boundary, numIter, J_part, d_half, leftThreadIndx, rightThreadIndx)
        iThread = omp_get_thread_num() + 1 
        loopSpecies: do j = 1, numberChargedParticles 
        ! Reset workspace for specific thread for adding J to domain
            particleList(j)%workSpace(1:numberXHalfNodes, iThread) = 0.0d0
            loopParticles: do i = 1, particleList(j)%N_p(iThread)
            ! initialize variables
                v_sub = particleList(j)%phaseSpace(2,i,iThread)
                l_sub = particleList(j)%phaseSpace(1,i,iThread)
                timePassed = 0.0d0
                ! Very small probability particle trajectories end on boundaries, so look along cell to determine
                l_cell = INT(l_sub + SIGN(1.d-12, v_sub))
                del_tau = del_t
                do while(del_tau > f_tol)
                    
                    dx_dl = world%dx_dl(l_cell)
                    ! Calculate maximum substep 
                    del_tau = MIN(del_tau, particleTimeStep(l_cell,j))
                    call particleSubStepPicard(l_sub-real(l_cell), v_sub, l_f, v_f, d_half, del_tau, solver%EField(l_cell), solver%EField(l_cell+1), dx_dl, &
                        particleList(j)%q_over_m, l_cell, FutureAtBoundaryBool, l_boundary, numIter)
                    
                    ! Add to J
                    J_part = l_f - l_sub
                    particleList(j)%workSpace(l_cell, iThread) = particleList(j)%workSpace(l_cell, iThread) + (1.0d0-d_half) * J_part
                    particleList(j)%workSpace(l_cell + 1, iThread) = particleList(j)%workSpace(l_cell + 1, iThread) + J_part * d_half
                    
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
                        END SELECT
                    end if
                    timePassed = timePassed + del_tau
                    del_tau = del_t - timePassed
                    ! Update parameters of particle trajectory for next substep
                    l_sub = l_f
                    v_sub = v_f
                end do
            end do loopParticles
        end do loopSpecies
        !$OMP barrier
        ! Concatenate workspace array to get J array
        leftThreadIndx = world%threadHalfNodeIndx(1,iThread)
        rightThreadIndx = world%threadHalfNodeIndx(2,iThread)
        particleList(1)%workSpace(leftThreadIndx:rightThreadIndx, 1) = SUM(particleList(1)%workSpace(leftThreadIndx:rightThreadIndx, :), DIM=2) * particleList(1)%q_times_wp
        do j = 2, numberChargedParticles
            particleList(1)%workSpace(leftThreadIndx:rightThreadIndx, 1) = particleList(1)%workSpace(leftThreadIndx:rightThreadIndx, 1) &
                + SUM(particleList(j)%workSpace(leftThreadIndx:rightThreadIndx, :), DIM=2) * particleList(j)%q_times_wp
        end do
        !$OMP end parallel
        ! Change J at boundaries based on boundary conditions
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
        ! apply smoothing
        if (world%gridSmoothBool) then
            call world%smoothField(particleList(1)%workSpace(:,1), solver%J)
            solver%J = solver%J/del_t
        else
            solver%J = particleList(1)%workSpace(:,1)/del_t
        end if
    end subroutine depositJ


    subroutine moveParticles(solver, particleList, world, del_t)
        ! Solve for phasespace of particles after convergence of potential solver
        class(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        type(Particle), intent(in out) :: particleList(:)
        real(real64), intent(in) :: del_t
        real(real64) :: l_f, l_sub, v_sub, v_f, timePassed, del_tau, f_tol, d_half, dx_dl
        integer(int32) :: j, i, l_cell, iThread, l_boundary, numSubStepAve(numberChargedParticles), numIter, funcEvalCounter(numberChargedParticles), delIdx, refIdx, N_p
        logical :: FutureAtBoundaryBool, refluxedBool, convergeTimeBool
        f_tol = del_t * 1.d-10
        ! Initialize counters for function evaluations and substep number
        numSubStepAve = 0
        funcEvalCounter = 0
        call solver%makeHalfTimeEField(particleList(1)%workSpace(1:NumberXHalfNodes,1), world)
        do j = 1, numberChargedParticles
            particleTimeStep(:, j) = 0.1d0 * world%dx_dl / SQRT(ABS(particleList(j)%q_over_m * (solver%EField(1:NumberXNodes) - solver%EField(2:NumberXHalfNodes))))
        end do
        ! Make counters private, otherwise have false sharing which causes significant slowdown! In general, keep variables private whenever possible
        !$OMP parallel private(iThread, i, j,l_f, l_sub, v_sub, v_f, d_half, timePassed, del_tau, l_cell, FutureAtBoundaryBool, &
            dx_dl, l_boundary, numIter, delIdx, refIdx, refluxedBool, N_p, convergeTimeBool) reduction(+:numSubStepAve, funcEvalCounter)
        iThread = omp_get_thread_num() + 1 
        loopSpecies: do j = 1, numberChargedParticles
            !initialize deletion index at dirichlet boundaries and reflux if applicable to neumann boundaries
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
                l_cell = INT(l_sub + SIGN(1.d-12, v_sub))
                refluxedBool = .false.
                del_tau = del_t
                convergeTimeBool = .true.
                do while(convergeTimeBool)
                    ! accumulate statistics for num sub steps
                    numSubStepAve(j) = numSubStepAve(j) + 1
                    
                    dx_dl = world%dx_dl(l_cell)
                    del_tau = MIN(del_tau, particleTimeStep(l_cell, j))
                    call particleSubStepPicard(l_sub - real(l_cell), v_sub, l_f, v_f, d_half, del_tau, solver%EField(l_cell), solver%EField(l_cell+1), dx_dl, &
                        particleList(j)%q_over_m, l_cell, FutureAtBoundaryBool, l_boundary, numIter)
                    ! Accumulate functional evaluation counter for each picard iteration
                    funcEvalCounter(j) = funcEvalCounter(j) + numIter
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
                            ! If reflux, save particle index 
                            if (.not. refluxedBool) then
                                refIdx = refIdx + 1
                                particleList(j)%refRecordIdx(refIdx, iThread) = i - delIdx
                                refluxedBool = .true.
                            end if
                        CASE(3)
                            l_f = REAL(ABS(l_boundary - NumberXHalfNodes - 1))
                            l_cell = ABS(l_cell - NumberXHalfNodes)
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
            particleList(j)%N_p(iThread) = particleList(j)%N_p(iThread) - delIdx
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