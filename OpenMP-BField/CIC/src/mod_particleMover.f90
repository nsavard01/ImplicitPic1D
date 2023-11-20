module mod_particleMover
    use iso_fortran_env, only: int32, real64, output_unit
    use constants
    use mod_BasicFunctions
    use mod_particle
    use mod_domain
    use mod_potentialSolver
    use omp_lib
    implicit none


contains
    
    pure function getEField(solver, l_p, l_cell, world) result(EField)
        !return EField of particle at logical position l_p (this is half distance), per particle since particle mover loop will be per particle
        type(potentialSolver), intent(in) :: solver
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: l_p
        integer(int32), intent(in) :: l_cell
        real(real64) :: EField, d
        d = l_p - real(l_cell, kind = real64)
        EField = (solver%EField(l_cell+1) * d + solver%EField(l_cell) * (1.0d0 - d))
    end function getEField

    subroutine subStepFreeDrift(l_sub, v_sub, l_f, d_half, v_f, v_half, del_tau, timePassed, E_x, q_over_m, l_cell, BField, B_mag, dx_dl, f_tol, AtBoundaryBool, &
        del_t, l_boundary)
        ! picard iteration particle mover, like Chen in 2015 paper
        integer(int32), intent(in) :: l_cell
        real(real64), intent(in) :: E_x, q_over_m, BField(3), B_mag, dx_dl, f_tol, del_t
        real(real64), intent(in out) :: l_sub, v_sub(3), l_f, v_f(3), v_half(3), del_tau, timePassed, d_half
        logical, intent(in out) :: AtBoundaryBool
        integer(int32), intent(in out) :: l_boundary
        real(real64) :: del_tau_max
        
        del_tau_max = del_t - timePassed
        if (v_sub(1) > 0) then
            l_boundary = l_cell + 1
        else
            l_boundary = l_cell
        end if
        del_tau = MIN(del_tau_max, (real(l_boundary) - l_sub) * dx_dl / v_sub(1))
        v_half = v_sub
        v_f = v_sub
        l_f = l_sub + v_half(1) * del_tau / dx_dl
        d_half = 0.5d0 * (l_f + l_sub) - real(l_cell)
        if (d_half > 1 .or. d_half < 0) then
            print *, 'd_half not good!'
            stop
        end if
        AtBoundaryBool = (del_tau < del_tau_max)
        if(.not. AtBoundaryBool) then
            if (MOD(l_f, 1.0d0) == 0.0d0) then
                AtBoundaryBool = .true.
                l_boundary = NINT(l_f)
            end if
        end if
        if (ABS(l_f - (real(l_cell) + 0.5d0)) > 0.5d0 + 1.d-12) then
            print *, 'Final l_f in drift not inside'
            print *, 'val sub is:', ABS(l_f - (real(l_cell) + 0.5d0))
            print *, 'l_f is:', l_f
            print *, 'l_sub is:', l_sub
            print *, 'l_cell is:', l_cell
            stop
        end if 
    end subroutine subStepFreeDrift
    
    subroutine GetRho(rho, l, world) 
        real(real64), intent(in out) :: rho(NumberXNodes)
        real(real64), intent(in) :: l
        type(Domain), intent(in) :: world
        integer(int32) :: i, j, l_center, l_left, l_right
        real(real64) :: d
        rho = 0.0d0
        l_center = INT(l)
        d = l - real(l_center)
        rho(l_center) = rho(l_center) + (-d**2 + d + 0.5d0)
        l_right = l_center + 1
        l_left = l_center - 1
        if (world%boundaryConditions(l_center) == 0 .and. world%boundaryConditions(l_right) == 0) then
            ! No Boundary either side
            rho(l_left) = rho(l_left) + 0.5d0 * (1.0d0 - d)**2
            rho(l_right) = rho(l_right) + 0.5d0 * d**2
        else if (world%boundaryConditions(l_right) == 1) then
            ! Dirichlet to right
            rho(l_center) = rho(l_center) - 0.5d0 * d**2
            rho(l_left) = rho(l_left) + 0.5d0 * (1.0d0 - d)**2
        else if (world%boundaryConditions(l_center) == 1) then
            !Dirichlet to left
            rho(l_center) = rho(l_center) - 0.5d0 * (1.0d0 - d)**2
            rho(l_right) = rho(l_right) + 0.5d0 * d**2
        else if (world%boundaryConditions(l_right) == 2) then
            !Neumann to right
            rho(l_center) = rho(l_center) + 0.5d0 * d**2
            rho(l_left) = rho(l_left) + 0.5d0 * (1.0d0 - d)**2
        else if (world%boundaryConditions(l_center) == 2) then
            !Neumann to left
            rho(l_center) = rho(l_center) + 0.5d0 * (1.0d0 - d)**2
            rho(l_right) = rho(l_right) + 0.5d0 * d**2
        end if
    end subroutine GetRho

    subroutine depositJ(solver, particleList, world, del_t)
        ! particle substepping procedure which deposits J, also used for general moving of particles
        ! boolDepositJ = true : deposit J
        ! boolDepositJ = false: only move particles (done after non-linear iteration), also delete particles for wall collisions (saves extra loop in collisions)
        class(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        type(Particle), intent(in out) :: particleList(:)
        real(real64), intent(in) :: del_t
        !a and c correspond to quadratic equations | l_alongV is nearest integer boundary along velocity component, away is opposite
        real(real64) :: l_f, l_sub, v_sub(3), v_f(3), timePassed, del_tau, q_over_m, f_tol, v_half(3), dx_dl, E_x, q_times_wp, J_part, d_half
        integer(int32) :: j, i, l_cell, iThread, l_boundary
        logical :: AtBoundaryBool
        solver%J = 0.0d0
        !Solve EField
        f_tol = del_t * 1.d-8
        loopSpecies: do j = 1, numberChargedParticles
            q_over_m = particleList(j)%q/particleList(j)%mass
            q_times_wp = particleList(j)%q * particleList(j)%w_p
            !$OMP parallel private(iThread, i, l_f, l_sub, v_sub, v_f, v_half, timePassed, del_tau, l_cell, AtBoundaryBool, dx_dl, E_x, l_boundary, J_part, d_half)
            iThread = omp_get_thread_num() + 1 
            loopParticles: do i = 1, particleList(j)%N_p(iThread)
                v_sub = particleList(j)%phaseSpace(2:4,i,iThread)
                l_sub = particleList(j)%phaseSpace(1,i,iThread)
                timePassed = 0.0d0
                do while((timePassed < del_t))
                    l_cell = INT(l_sub)
                    E_x = solver%EField(l_cell)
                    dx_dl = world%dx_dl(l_cell)

                    ! Start AA
                    call subStepFreeDrift(l_sub, v_sub, l_f, d_half, v_f, v_half, del_tau, timePassed, E_x, q_over_m, l_cell, solver%BField, solver%BFieldMag, dx_dl, f_tol, AtBoundaryBool, &
                        del_t, l_boundary)
                    
                    J_part = q_times_wp * (v_half(1))*del_tau/world%dx_dl(l_cell)/del_t
                    solver%J(l_cell, iThread) = solver%J(l_cell, iThread) + J_part * (1.0d0 - d_half)
                    solver%J(l_cell+1, iThread) = solver%J(l_cell+1, iThread) + J_part * (d_half)
                    if (AtBoundaryBool) then
                        SELECT CASE (world%boundaryConditions(l_boundary))
                        CASE(0)
                            l_f = real(l_boundary) + SIGN(1.d-12, v_f(1))
                        CASE(1,4)
                            exit
                        CASE(2)
                            if (l_boundary == NumberXNodes+1) then
                                v_f(1) = -ABS(v_f(1))
                                l_f = real(l_boundary) - 1.d-12
                            else
                                v_f(1) = ABS(v_f(1))
                                l_f = real(l_boundary) + 1.d-12
                            end if
                            v_f(2:3) = -v_f(2:3)  
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
                    timePassed = timePassed + del_tau
                    ! now final position/velocity becomes next starting position/velocity
                    l_sub = l_f
                    v_sub = v_f
                    if ((l_f < 1) .or. (l_f > NumberXNodes+1)) then
                        print *, 'l_f is:', l_f
                        stop "Have particles travelling outside the domain in depositJ!"
                    end if
                end do
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
        real(real64) :: l_f, l_sub, v_sub(3), v_f(3), timePassed, del_tau, q_over_m, f_tol, d_half, v_half(3), dx_dl, E_x
        integer(int32) :: j, i, l_cell, iThread, delIdx, l_boundary
        logical :: AtBoundaryBool
        ! Solve EField
        f_tol = del_t * 1.d-8
        loopSpecies: do j = 1, numberChargedParticles
            q_over_m = particleList(j)%q/particleList(j)%mass
            !$OMP parallel private(iThread, i, l_f, l_sub, d_half, v_sub, v_f, v_half, timePassed, del_tau, l_cell, AtBoundaryBool, delIdx, dx_dl, E_x, l_boundary)
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
                    
                    ! AA particle mover
                    call subStepFreeDrift(l_sub, v_sub, l_f, d_half, v_f, v_half, del_tau, timePassed, E_x, q_over_m, l_cell, solver%BField, solver%BFieldMag, dx_dl, f_tol, AtBoundaryBool, &
                        del_t, l_boundary)
                    
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
                            if (l_boundary == NumberXNodes+1) then
                                v_f(1) = -ABS(v_f(1))
                                l_f = real(l_boundary) - 1.d-12
                            else
                                v_f(1) = ABS(v_f(1))
                                l_f = real(l_boundary) + 1.d-12
                            end if
                            v_f(2:3) = -v_f(2:3)
                            ! particleList(j)%refIdx(iThread) = particleList(j)%refIdx(iThread) + 1
                            ! particleList(j)%refRecordIdx(particleList(j)%refIdx(iThread), iThread) = i - delIdx
                        CASE(3)
                            l_f = REAL(ABS(l_boundary - real(NumberXNodes, kind = real64) - 1)) + SIGN(1.d-12, v_f(1))
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
                    if ((l_f < 1) .or. (l_f > NumberXNodes+1)) then
                        print *, "Have particles travelling outside domain in moveparticles!"
                        print *, 'del_tau is:', del_tau
                        print *, 'l_sub:', l_sub
                        print *, 'l_f:', l_f
                        print *, 'v_sub:', v_sub
                        print *, 'v_half:', v_half
                        print *, 'v_f:',v_f
                        stop
                    end if
                    timePassed = timePassed + del_tau
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

end module mod_particleMover