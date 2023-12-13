module mod_particleInjection
    use iso_fortran_env, only: int32, real64
    use constants
    use mod_BasicFunctions
    use mod_domain
    use mod_particle
    use mod_Scheme
    use omp_lib
    implicit none
    real(real64), allocatable :: energyAddColl(:)
    logical :: addLostPartBool, refluxPartBool, injectionBool
    integer(int32) :: numFluxParticlesHigh, numFluxParticlesLow
    real(real64) :: injectionFlux, injectionR

contains

    subroutine addMaxwellianLostParticles(particleList, T_e, T_i, irand, world)
        ! Add power to all particles in domain
        type(Particle), intent(in out) :: particleList(2)
        type(Domain), intent(in) :: world
        integer(int32), intent(in out) :: irand(numThread)
        real(real64), intent(in) :: T_e, T_i
        integer(int32) :: i, iThread
        real(real64) :: x_random
        !$OMP parallel private(iThread, i, x_random)
        iThread = omp_get_thread_num() + 1
        do i=1, particleList(1)%delIdx(iThread)
            call getMaxwellianSample(particleList(1)%phaseSpace(2:4, particleList(1)%N_p(iThread) + i, iThread), particleList(1)%mass, T_e, irand(iThread))
            call getMaxwellianFluxSample(particleList(2)%phaseSpace(2:4, particleList(2)%N_p(iThread) + i, iThread), particleList(2)%mass, T_i, irand(iThread))
            energyAddColl(iThread) = energyAddColl(iThread) + (SUM(particleList(1)%phaseSpace(2:4, particleList(1)%N_p(iThread) + i, iThread)**2) * particleList(1)%mass * particleList(1)%w_p &
            + SUM(particleList(2)%phaseSpace(2:4, particleList(2)%N_p(iThread) + i, iThread)**2) * particleList(2)%mass * particleList(2)%w_p) * 0.5d0
            if (schemeNum == 0) then
                x_random = world%grid(1) + (world%grid(NumberXNodes) - world%grid(NumberXNodes)) * ran2(irand(iThread))
            else
                x_random = world%grid(1) - 0.5d0*world%dx_dl(1) + (world%grid(NumberXNodes) + 0.5d0*world%dx_dl(NumberXNodes) - world%grid(1) + 0.5d0 * world%dx_dl(1)) * ran2(irand(iThread))
            end if
            particleList(1)%phaseSpace(1, particleList(1)%N_p(iThread) + i, iThread) = world%getLFromX(x_random)
            particleList(2)%phaseSpace(1, particleList(2)%N_p(iThread) + i, iThread) = particleList(1)%phaseSpace(1, particleList(1)%N_p(iThread) + i, iThread)
        end do
        particleList(2)%N_p(iThread) = particleList(2)%N_p(iThread) + particleList(1)%delIdx(iThread)
        particleList(1)%N_p(iThread) = particleList(1)%N_p(iThread) + particleList(1)%delIdx(iThread)
        !$OMP end parallel
    end subroutine addMaxwellianLostParticles

    subroutine refluxParticles(particleList, T_e, T_i, irand, world)
        type(Particle), intent(in out) :: particleList(2)
        type(Domain), intent(in) :: world
        integer(int32), intent(in out) :: irand(numThread)
        real(real64), intent(in) :: T_e, T_i
        integer(int32) :: i,j, iThread
        do j = 1, 2
            !$OMP parallel private(iThread, i)
            iThread = omp_get_thread_num() + 1
            do i = 1, particleList(j)%refIdx(iThread)
                energyAddColl(iThread) = energyAddColl(iThread) - (SUM(particleList(j)%phaseSpace(2:4, particleList(j)%refRecordIdx(i, iThread), iThread)**2) * particleList(j)%mass * particleList(j)%w_p) * 0.5d0
                if (j == 1) then
                    call getMaxwellianFluxSample(particleList(j)%phaseSpace(2:4, particleList(j)%refRecordIdx(i, iThread), iThread), particleList(j)%mass, T_e, irand(iThread))
                else
                    call getMaxwellianFluxSample(particleList(j)%phaseSpace(2:4, particleList(j)%refRecordIdx(i, iThread), iThread), particleList(j)%mass, T_i, irand(iThread))
                end if
                if (world%boundaryConditions(1) == 2) then
                    particleList(j)%phaseSpace(2, particleList(j)%refRecordIdx(i, iThread), iThread) = ABS(particleList(j)%phaseSpace(2, particleList(j)%refRecordIdx(i, iThread), iThread))
                else
                    particleList(j)%phaseSpace(2, particleList(j)%refRecordIdx(i, iThread), iThread) = -ABS(particleList(j)%phaseSpace(2, particleList(j)%refRecordIdx(i, iThread), iThread))
                end if
                energyAddColl(iThread) = energyAddColl(iThread) + (SUM(particleList(j)%phaseSpace(2:4, particleList(j)%refRecordIdx(i, iThread), iThread)**2) * particleList(j)%mass * particleList(j)%w_p) * 0.5d0
            end do
            !$OMP end parallel
        end do
    end subroutine refluxParticles

    subroutine injectAtBoundary(particleList, T_e, T_i, irand, world, del_t, BFieldAngle)
        type(Particle), intent(in out) :: particleList(2)
        type(Domain), intent(in) :: world
        integer(int32), intent(in out) :: irand(numThread)
        real(real64), intent(in) :: T_e, T_i, del_t, BFieldAngle
        integer(int32) :: i,j, iThread, numInjected
        real(real64) :: ionSoundSpeed, shift_rand, cosAngle, sinAngle
        cosAngle = COS(BFieldAngle)
        sinAngle = SIN(BFieldAngle)
        ionSoundSpeed = SQRT(T_e * e / particleList(2)%mass)
        do j = 1, 2
            !$OMP parallel private(iThread, i)
            iThread = omp_get_thread_num() + 1
            do i = 1, particleList(j)%refIdx(iThread)
                !energyAddColl(iThread) = energyAddColl(iThread) - (SUM(particleList(j)%phaseSpace(2:4, particleList(j)%refRecordIdx(i, iThread), iThread)**2) * particleList(j)%mass * particleList(j)%w_p) * 0.5d0
                if (j == 1) then
                    call getMaxwellianFluxSample(particleList(j)%phaseSpace(2:4, particleList(j)%refRecordIdx(i, iThread), iThread), particleList(j)%mass, T_e, irand(iThread))
                else
                    particleList(j)%phaseSpace(2, particleList(j)%refRecordIdx(i, iThread), iThread) = ionSoundSpeed * cosAngle
                    particleList(j)%phaseSpace(3, particleList(j)%refRecordIdx(i, iThread), iThread) = -ionSoundSpeed * sinAngle
                    particleList(j)%phaseSpace(4, particleList(j)%refRecordIdx(i, iThread), iThread) = 0.0d0
                end if
                if (world%boundaryConditions(1) == 2) then
                    particleList(j)%phaseSpace(2, particleList(j)%refRecordIdx(i, iThread), iThread) = ABS(particleList(j)%phaseSpace(2, particleList(j)%refRecordIdx(i, iThread), iThread))
                else
                    particleList(j)%phaseSpace(2, particleList(j)%refRecordIdx(i, iThread), iThread) = -ABS(particleList(j)%phaseSpace(2, particleList(j)%refRecordIdx(i, iThread), iThread))
                end if
                !energyAddColl(iThread) = energyAddColl(iThread) + (SUM(particleList(j)%phaseSpace(2:4, particleList(j)%refRecordIdx(i, iThread), iThread)**2) * particleList(j)%mass * particleList(j)%w_p) * 0.5d0
            end do
            !$OMP end parallel
        end do
        ! Injection of particles

        !$OMP parallel private(iThread, i, numInjected, shift_rand)
        iThread = omp_get_thread_num() + 1
        if (ran2(irand(iThread)) > injectionR) then
            numInjected = numFluxParticlesLow
        else
            numInjected = numFluxParticlesHigh
        end if
        do i = 1, numInjected
            shift_rand = del_t * ran2(irand(iThread)) * ionSoundSpeed * cosAngle
            call getMaxwellianFluxSample(particleList(1)%phaseSpace(2:4, particleList(1)%N_p(iThread) + i, iThread), particleList(1)%mass, T_e, irand(iThread))
            particleList(2)%phaseSpace(2, particleList(2)%N_p(iThread) + i, iThread) = ionSoundSpeed * cosAngle
            particleList(2)%phaseSpace(3, particleList(2)%N_p(iThread) + i, iThread) = -ionSoundSpeed * sinAngle
            particleList(2)%phaseSpace(4, particleList(2)%N_p(iThread) + i, iThread) = 0.0d0
            if (world%boundaryConditions(1) == 2 .or. world%boundaryConditions(1) == 4) then
                particleList(1)%phaseSpace(2, particleList(1)%N_p(iThread) + i, iThread) = ABS(particleList(1)%phaseSpace(2, particleList(1)%N_p(iThread) + i, iThread))
                particleList(2)%phaseSpace(2, particleList(2)%N_p(iThread) + i, iThread) = ABS(particleList(2)%phaseSpace(2, particleList(2)%N_p(iThread) + i, iThread))
                if (schemeNum == 0) then
                    particleList(1)%phaseSpace(1, particleList(1)%N_p(iThread) + i, iThread) = world%getLFromX(world%grid(1) + shift_rand)
                else
                    particleList(1)%phaseSpace(1, particleList(1)%N_p(iThread) + i, iThread) = world%getLFromX(world%grid(1) - 0.5d0 * world%dx_dl(1) + shift_rand)
                end if
                particleList(2)%phaseSpace(1, particleList(2)%N_p(iThread) + i, iThread) = particleList(1)%phaseSpace(1, particleList(1)%N_p(iThread) + i, iThread)
            else
                particleList(1)%phaseSpace(2, particleList(1)%N_p(iThread) + i, iThread) = -ABS(particleList(1)%phaseSpace(2, particleList(1)%N_p(iThread) + i, iThread))
                particleList(2)%phaseSpace(2, particleList(2)%N_p(iThread) + i, iThread) = -ABS(particleList(2)%phaseSpace(2, particleList(2)%N_p(iThread) + i, iThread))
                if (schemeNum == 0) then
                    particleList(1)%phaseSpace(1, particleList(1)%N_p(iThread) + i, iThread) = world%getLFromX(world%grid(NumberXNodes) - shift_rand)
                else
                    particleList(1)%phaseSpace(1, particleList(1)%N_p(iThread) + i, iThread) = world%getLFromX(world%grid(NumberXNodes) + 0.5d0 * world%dx_dl(NumberXNodes) - shift_rand)
                end if
                particleList(2)%phaseSpace(1, particleList(2)%N_p(iThread) + i, iThread) = particleList(1)%phaseSpace(1, particleList(1)%N_p(iThread) + i, iThread)
            end if
        end do
        particleList(1)%N_p(iThread) = particleList(1)%N_p(iThread) + numInjected
        particleList(2)%N_p(iThread) = particleList(2)%N_p(iThread) + numInjected
        !$OMP end parallel
    end subroutine injectAtBoundary

    subroutine injectUniformFlux(particleList, T_e, T_i, irand, world)
        type(Particle), intent(in out) :: particleList(2)
        type(Domain), intent(in) :: world
        integer(int32), intent(in out) :: irand(numThread)
        real(real64), intent(in) :: T_e, T_i
        integer(int32) :: i,j, iThread, numInjected
        real(real64) :: ionSoundSpeed, x_random
        ! Reflux particles
        ionSoundSpeed = SQRT(T_e * e / particleList(2)%mass)
        do j = 1, 2
            !$OMP parallel private(iThread, i)
            iThread = omp_get_thread_num() + 1
            do i = 1, particleList(j)%refIdx(iThread)
                !energyAddColl(iThread) = energyAddColl(iThread) - (SUM(particleList(j)%phaseSpace(2:4, particleList(j)%refRecordIdx(i, iThread), iThread)**2) * particleList(j)%mass * particleList(j)%w_p) * 0.5d0
                if (j == 1) then
                    call getMaxwellianFluxSample(particleList(j)%phaseSpace(2:4, particleList(j)%refRecordIdx(i, iThread), iThread), particleList(j)%mass, T_e, irand(iThread))
                else
                    particleList(j)%phaseSpace(2, particleList(j)%refRecordIdx(i, iThread), iThread) = ionSoundSpeed
                    particleList(j)%phaseSpace(3, particleList(j)%refRecordIdx(i, iThread), iThread) = 0.0d0
                    particleList(j)%phaseSpace(4, particleList(j)%refRecordIdx(i, iThread), iThread) = 0.0d0
                end if
                if (world%boundaryConditions(1) == 2 .or. world%boundaryConditions(1) == 4) then
                    particleList(j)%phaseSpace(2, particleList(j)%refRecordIdx(i, iThread), iThread) = ABS(particleList(j)%phaseSpace(2, particleList(j)%refRecordIdx(i, iThread), iThread))
                else
                    particleList(j)%phaseSpace(2, particleList(j)%refRecordIdx(i, iThread), iThread) = -ABS(particleList(j)%phaseSpace(2, particleList(j)%refRecordIdx(i, iThread), iThread))
                end if
                !energyAddColl(iThread) = energyAddColl(iThread) + (SUM(particleList(j)%phaseSpace(2:4, particleList(j)%refRecordIdx(i, iThread), iThread)**2) * particleList(j)%mass * particleList(j)%w_p) * 0.5d0
            end do
            !$OMP end parallel
        end do
        ! Injection of particles
        
        !$OMP parallel private(iThread, i, numInjected, x_random)
        iThread = omp_get_thread_num() + 1
        if (ran2(irand(iThread)) > injectionR) then
            numInjected = numFluxParticlesLow
        else
            numInjected = numFluxParticlesHigh
        end if
        do i = 1, numInjected
            call getMaxwellianSample(particleList(1)%phaseSpace(2:4, particleList(1)%N_p(iThread) + i, iThread), particleList(1)%mass, T_e, irand(iThread))
            if (world%boundaryConditions(NumberXNodes) == 2) then
                particleList(2)%phaseSpace(2, particleList(2)%N_p(iThread) + i, iThread) = -ionSoundSpeed
            else
                particleList(2)%phaseSpace(2, particleList(2)%N_p(iThread) + i, iThread) = ionSoundSpeed
            end if
            particleList(2)%phaseSpace(3, particleList(2)%N_p(iThread) + i, iThread) = 0.0d0
            particleList(2)%phaseSpace(4, particleList(2)%N_p(iThread) + i, iThread) = 0.0d0
            x_random = world%grid(NumberXNodes) * ran2(irand(iThread))
            particleList(1)%phaseSpace(1, particleList(1)%N_p(iThread) + i, iThread) = world%getLFromX(x_random)
            particleList(2)%phaseSpace(1, particleList(2)%N_p(iThread) + i, iThread) = particleList(1)%phaseSpace(1, particleList(1)%N_p(iThread) + i, iThread)
            energyAddColl(iThread) = energyAddColl(iThread) + (SUM(particleList(1)%phaseSpace(2:4, particleList(1)%N_p(iThread) + i, iThread)**2) * particleList(1)%mass * particleList(1)%w_p &
            + SUM(particleList(2)%phaseSpace(2:4, particleList(2)%N_p(iThread) + i, iThread)**2) * particleList(2)%mass * particleList(2)%w_p) * 0.5d0
        end do
        particleList(1)%N_p(iThread) = particleList(1)%N_p(iThread) + numInjected
        particleList(2)%N_p(iThread) = particleList(2)%N_p(iThread) + numInjected
        !$OMP end parallel
    end subroutine injectUniformFlux


end module mod_particleInjection