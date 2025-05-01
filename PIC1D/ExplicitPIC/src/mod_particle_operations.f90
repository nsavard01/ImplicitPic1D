module mod_particle_operations
    ! Module for running simulation over several time steps
    use iso_fortran_env, only: int32, real64, int64
    use constants
    use mod_BasicFunctions
    use mod_particle
    use mod_domain
    use mod_potentialSolver
    use omp_lib
    implicit none

contains


    subroutine loadParticleDensity(particleList, world, reset)
        ! load particle densities
        type(Particle), intent(in out) :: particleList(:)
        type(Domain), intent(in) :: world
        logical, intent(in) :: reset
        integer(int32) :: i,j, l_left, l_right, iThread, startIdx
        real(real64) :: d
        !$OMP parallel private(iThread, j, i, l_left, l_right, d)
        iThread = omp_get_thread_num() + 1
        do i=1, size(particleList)
            if (reset) then
                particleList(i)%densities(:,iThread) = 0.0d0
            end if
            do j = 1, particleList(i)%N_p(iThread)
                l_left = INT(particleList(i)%phaseSpace(1,j, iThread))
                l_right = l_left + 1
                d = particleList(i)%phaseSpace(1,j, iThread) - l_left
                particleList(i)%densities(l_left, iThread) = particleList(i)%densities(l_left, iThread) + (1.0d0-d)
                particleList(i)%densities(l_right, iThread) = particleList(i)%densities(l_right, iThread) + d
            end do
        end do
        !$OMP end parallel

    end subroutine loadParticleDensity

    subroutine move_particle_thread(part, EField, world, del_t, iThread)
        ! particle mover using the fields
        type(Particle), intent(in out) :: part
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: EField(NumberXNodes)
        real(real64), intent(in) :: del_t
        integer(int32), intent(in) :: iThread
        integer(int32) :: i, delIdx, refIdx, N_p, xi_left
        real(real64) :: v_prime, q_over_m, partLoc, d, EField_part

        q_over_m = part%q_over_m
        delIdx = 0
        refIdx = 0
        part%wallLoss(:, iThread) = 0
        part%energyLoss(:, iThread) = 0.0d0
        N_p = part%N_p(iThread)
        loopParticles: do i = 1, N_p
            ! First velocity change
            partLoc = part%phaseSpace(1, i, iThread)
            xi_left = INT(partLoc)
            d = partLoc - xi_left
            EField_part = EField(xi_left) * (1.0d0 - d) + EField(xi_left + 1) * d
            v_prime = part%phaseSpace(2, i, iThread) + q_over_m * EField_part * del_t
            ! Get new position
            partLoc = partLoc + v_prime * del_t/world%delX

            ! Check if outside boundary
            if (partLoc <= 1) then
                SELECT CASE (world%boundaryConditions(1))
                CASE(1,4)
                    part%energyLoss(1, iThread) = part%energyLoss(1, iThread) + v_prime**2 + SUM(part%phaseSpace(3:4, i, iThread)**2)
                    part%wallLoss(1, iThread) = part%wallLoss(1, iThread) + 1 !C/m^2 in 1D
                    part%momentumLoss(1,iThread) = part%momentumLoss(1,iThread) + v_prime
                    delIdx = delIdx + 1
                CASE(2)
                    partLoc = 2.0d0 - partLoc
                    v_prime = -v_prime
                    refIdx = refIdx + 1
                    part%refRecordIdx(refIdx, iThread) = i - delIdx
                CASE(3)
                    partLoc = MODULO(partLoc - 2.0d0, real(NumberXNodes, kind = real64)) + 1
                END SELECT
            else if ((partLoc >= NumberXNodes)) then
                SELECT CASE (world%boundaryConditions(NumberXNodes))
                CASE(1,4)
                    part%energyLoss(2, iThread) = part%energyLoss(2, iThread) + v_prime**2 + SUM(part%phaseSpace(3:4, i, iThread)**2)
                    part%wallLoss(2, iThread) = part%wallLoss(2, iThread) + 1 !C/m^2 in 1D
                    part%momentumLoss(2,iThread) = part%momentumLoss(2,iThread) + v_prime
                    delIdx = delIdx + 1
                CASE(2)
                    partLoc = 2.0d0 * NumberXNodes - partLoc
                    v_prime = -v_prime
                    refIdx = refIdx + 1
                    part%refRecordIdx(refIdx, iThread) = i - delIdx
                CASE(3)
                    partLoc = MODULO(partLoc, real(NumberXNodes, kind = real64)) + 1
                END SELECT
            end if
            if (partLoc > 1 .and. partLoc < NumberXNodes) then
                part%phaseSpace(1, i-delIdx, iThread) = partLoc
                part%phaseSpace(2, i-delIdx, iThread) = v_prime
                part%phaseSpace(3:4, i-delIdx, iThread) = part%phaseSpace(3:4, i, iThread)
            end if
        end do loopParticles
        part%N_p(iThread) = N_p - delIdx
        part%delIdx(iThread) = delIdx
        part%refIdx(iThread) = refIdx
    end subroutine move_particle_thread

    subroutine moveParticles(particleList, EField, world, del_t)
        type(Particle), intent(in out) :: particleList(:)
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: EField(NumberXNodes)
        real(real64), intent(in) :: del_t
        integer(int32) :: iThread, j, size_list
        size_list = size(particleList)
        !$OMP parallel private(iThread, j)
        iThread = omp_get_thread_num() + 1
        do j = 1, size_list
            call move_particle_thread(particleList(j), EField, world, del_t, iThread)
        end do
        !$OMP end parallel
        ! Update particle accumulation stats
        do j = 1, size_list
            particleList(j)%numToCollide = particleList(j)%N_p
            particleList(j)%accumEnergyLoss = particleList(j)%accumEnergyLoss + SUM(particleList(j)%energyLoss, DIM=2)
            particleList(j)%accumWallLoss = particleList(j)%accumWallLoss + SUM(particleList(j)%wallLoss, DIM=2)
        end do


    end subroutine moveParticles


    subroutine interpolate_particle_thread(part, iThread)
        type(Particle), intent(in out) :: part
        integer(int32), intent(in) :: iThread
        integer(int32) :: j, l_left, l_right, startIdx, endIdx
        real(real64) :: d 
        startIdx = part%startIdx(iThread)
        if (startIdx == 1) part%workSpace(:,iThread) = 0.0d0 ! reset if going through all particles
        endIdx = part%N_p(iThread)
        do j = startIdx, endIdx
            l_left = INT(part%phaseSpace(1, j, iThread))
            l_right = l_left + 1
            d = part%phaseSpace(1, j, iThread) - real(l_left)
            part%workSpace(l_left, iThread) = part%workSpace(l_left, iThread) + (1.0d0-d)
            part%workSpace(l_right, iThread) = part%workSpace(l_right, iThread) + d
        end do 
    end subroutine interpolate_particle_thread

    subroutine depositRho(solver, particleList, world) 
        ! calculate rho from particle locations
        class(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        type(Particle), intent(in out) :: particleList(numberChargedParticles)
        integer(int32) :: i, leftThreadIndx, rightThreadIndx, iThread
        solver % rho = solver%rho_const
        !$OMP parallel private(i, leftThreadIndx, rightThreadIndx, iThread)
        iThread = omp_get_thread_num() + 1
        do i=1, numberChargedParticles
            call interpolate_particle_thread(particleList(i), iThread)
        end do
        !$OMP barrier
        leftThreadIndx = world%threadNodeIndx(1,iThread)
        rightThreadIndx = world%threadNodeIndx(2,iThread)
        solver%rho(leftThreadIndx:rightThreadIndx) = SUM(particleList(1)%workSpace(leftThreadIndx:rightThreadIndx, :), DIM=2) * particleList(1)%q_times_wp
        do i = 2, numberChargedParticles
            solver%rho(leftThreadIndx:rightThreadIndx) = solver%rho(leftThreadIndx:rightThreadIndx) &
                + SUM(particleList(i)%workSpace(leftThreadIndx:rightThreadIndx, :), DIM=2) * particleList(i)%q_times_wp
        end do
        !$OMP end parallel
    end subroutine depositRho

    subroutine initialVRewind(solver, particleList, del_t, ionStepMult)
        ! Rewind particle velocities by half step using E-Field
        type(Particle), intent(in out) :: particleList(:)
        type(potentialSolver), intent(in) :: solver
        real(real64), intent(in) :: del_t
        integer(int32), intent(in) :: ionStepMult
        real(real64) :: del_t_temp
        real(real64) :: q_m_ratio
        integer(int32) :: j, i, iThread
        loopSpecies: do j = 1, numberChargedParticles
            q_m_ratio = particleList(j)%q/particleList(j)%mass
            if (j == 1) then
                del_t_temp = del_t
            else
                del_t_temp = ionStepMult * del_t
            end if
            !$OMP parallel private(iThread, i)
            iThread = omp_get_thread_num() + 1
            loopParticles: do i = 1, particleList(j)%N_p(iThread)
                particleList(j)%phaseSpace(2, i, iThread) = particleList(j)%phaseSpace(2, i, iThread) - 0.5d0 * (q_m_ratio) * solver%getEField(particleList(j)%phaseSpace(1, i, iThread)) * del_t_temp
            end do loopParticles
            !$OMP end parallel
        end do loopSpecies
    end subroutine initialVRewind




end module mod_particle_operations