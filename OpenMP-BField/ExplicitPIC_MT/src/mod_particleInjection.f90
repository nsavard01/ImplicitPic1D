module mod_particleInjection
    use iso_fortran_env, only: int32, real64
    use constants
    use mod_BasicFunctions
    use mod_domain
    use mod_particle
    use omp_lib
    use mod_mt19937
    implicit none
    logical :: addLostPartBool, refluxPartBool, injectionBool, uniformInjectionBool, heatingBool
    integer(int32) :: numFluxParticlesHigh, numFluxParticlesLow
    real(real64) :: injectionFlux, injectionR, FractionFreqHeating

contains

    subroutine addMaxwellianLostParticles(particleList, T_e, T_i, randGen, world)
        ! Add power to all particles in domain
        type(Particle), intent(in out) :: particleList(2)
        type(Domain), intent(in) :: world
        type(mt19937), intent(in out) :: randGen(numThread)
        real(real64), intent(in) :: T_e, T_i
        integer(int32) :: i, iThread
        !$OMP parallel private(iThread, i)
        iThread = omp_get_thread_num() + 1
        do i=1, particleList(1)%delIdx(iThread)
            call getMaxwellianSample(particleList(1)%phaseSpace(2:4, particleList(1)%N_p(iThread) + i, iThread), particleList(1)%mass, T_e, irand(iThread))
            call getMaxwellianFluxSample(particleList(2)%phaseSpace(2:4, particleList(2)%N_p(iThread) + i, iThread), particleList(2)%mass, T_i, irand(iThread))
            particleList(1)%phaseSpace(1, particleList(1)%N_p(iThread) + i, iThread) = ran2(irand(iThread)) * real(NumberXNodes - 1) + 1.0d0
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
        ! Reflux particles
        ionSoundSpeed = SQRT(T_e * e / particleList(2)%mass)
        cosAngle = COS(BFieldAngle)
        sinAngle = SIN(BFieldAngle)
        do j = 1, 2
            !$OMP parallel private(iThread, i)
            iThread = omp_get_thread_num() + 1
            do i = 1, particleList(j)%refIdx(iThread)
                !energyAddColl(iThread) = energyAddColl(iThread) - (SUM(particleList(j)%phaseSpace(2:4, particleList(j)%refRecordIdx(i, iThread), iThread)**2) * particleList(j)%mass * particleList(j)%w_p) * 0.5d0
                if (j == 1) then
                    call getMaxwellianFluxSample(particleList(j)%phaseSpace(2:4, particleList(j)%refRecordIdx(i, iThread), iThread), particleList(j)%mass, T_e, irand(iThread))
                    ! tempVec(1) = cosAngle * particleList(j)%phaseSpace(2, particleList(j)%refRecordIdx(i, iThread), iThread) - sinAngle * particleList(j)%phaseSpace(3, particleList(j)%refRecordIdx(i, iThread), iThread)
                    ! tempVec(2) = sinAngle * particleList(j)%phaseSpace(2, particleList(j)%refRecordIdx(i, iThread), iThread) + cosAngle * particleList(j)%phaseSpace(3, particleList(j)%refRecordIdx(i, iThread), iThread)
                    ! particleList(j)%phaseSpace(2:3, particleList(j)%refRecordIdx(i, iThread), iThread) = tempVec
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
            shift_rand = del_t * ran2(irand(iThread)) * ionSoundSpeed * cosAngle/ world%delX
            call getMaxwellianFluxSample(particleList(1)%phaseSpace(2:4, particleList(1)%N_p(iThread) + i, iThread), particleList(1)%mass, T_e, irand(iThread))
            ! tempVec(1) = cosAngle * particleList(1)%phaseSpace(2, particleList(1)%N_p(iThread) + i, iThread) - sinAngle * particleList(1)%phaseSpace(3, particleList(1)%N_p(iThread) + i, iThread)
            ! tempVec(2) = sinAngle * particleList(1)%phaseSpace(2, particleList(1)%N_p(iThread) + i, iThread) + cosAngle * particleList(1)%phaseSpace(3, particleList(1)%N_p(iThread) + i, iThread)
            ! particleList(1)%phaseSpace(2:3, particleList(1)%N_p(iThread) + i, iThread) = tempVec
            particleList(2)%phaseSpace(2, particleList(2)%N_p(iThread) + i, iThread) = ionSoundSpeed * cosAngle
            particleList(2)%phaseSpace(3, particleList(2)%N_p(iThread) + i, iThread) = -ionSoundSpeed * sinAngle
            particleList(2)%phaseSpace(4, particleList(2)%N_p(iThread) + i, iThread) = 0.0d0
            if (world%boundaryConditions(1) == 2) then
                particleList(1)%phaseSpace(2, particleList(1)%N_p(iThread) + i, iThread) = ABS(particleList(1)%phaseSpace(2, particleList(1)%N_p(iThread) + i, iThread))
                particleList(2)%phaseSpace(2, particleList(2)%N_p(iThread) + i, iThread) = ABS(particleList(2)%phaseSpace(2, particleList(2)%N_p(iThread) + i, iThread))
                particleList(1)%phaseSpace(1, particleList(1)%N_p(iThread) + i, iThread) = 1.0d0 + shift_rand
                particleList(2)%phaseSpace(1, particleList(2)%N_p(iThread) + i, iThread) = particleList(1)%phaseSpace(1, particleList(1)%N_p(iThread) + i, iThread)
            else
                particleList(1)%phaseSpace(2, particleList(1)%N_p(iThread) + i, iThread) = -ABS(particleList(1)%phaseSpace(2, particleList(1)%N_p(iThread) + i, iThread))
                particleList(2)%phaseSpace(2, particleList(2)%N_p(iThread) + i, iThread) = -ABS(particleList(2)%phaseSpace(2, particleList(2)%N_p(iThread) + i, iThread))
                particleList(1)%phaseSpace(1, particleList(1)%N_p(iThread) + i, iThread) = real(NumberXNodes) - shift_rand
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
        integer(int32) :: i, iThread, numInjected
        ! Reflux particles
        
        !$OMP parallel private(iThread, i, numInjected)
        iThread = omp_get_thread_num() + 1
        if (ran2(irand(iThread)) > injectionR) then
            numInjected = numFluxParticlesLow
        else
            numInjected = numFluxParticlesHigh
        end if
        do i = 1, numInjected
            call getMaxwellianSample(particleList(1)%phaseSpace(2:4, particleList(1)%N_p(iThread) + i, iThread), particleList(1)%mass, 20.0d0 * T_e, irand(iThread))
            call getMaxwellianSample(particleList(2)%phaseSpace(2:4, particleList(2)%N_p(iThread) + i, iThread), particleList(2)%mass, T_i, irand(iThread))
            particleList(1)%phaseSpace(1, particleList(1)%N_p(iThread) + i, iThread) = ran2(irand(iThread)) * real(NumberXNodes - 1) + 1.0d0
            particleList(2)%phaseSpace(1, particleList(2)%N_p(iThread) + i, iThread) = particleList(1)%phaseSpace(1, particleList(1)%N_p(iThread) + i, iThread)
        end do
        particleList(1)%N_p(iThread) = particleList(1)%N_p(iThread) + numInjected
        particleList(2)%N_p(iThread) = particleList(2)%N_p(iThread) + numInjected
        !$OMP end parallel
    end subroutine injectUniformFlux


    subroutine maxwellianHeating(species, FractionFreqHeating, irand, FracFreq, T)
        ! Add power to all particles in domain
        type(Particle), intent(in out) :: species
        integer(int32), intent(in out) :: irand(numThread)
        real(real64), intent(in) :: FractionFreqHeating, FracFreq, T
        real(real64) :: R, compValue
        integer(int32) :: i, iThread
        compValue = FracFreq/FractionFreqHeating
        !$OMP parallel private(iThread, i, R)
        iThread = omp_get_thread_num() + 1
        do i = 1, species%N_p(iThread)
            R = ran2(irand(iThread))
            if (R < compValue) then
                call getMaxwellianSample(species%phaseSpace(2:4, i, iThread), species%mass, T, irand(iThread))
            end if
        end do
        !$OMP end parallel
    end subroutine maxwellianHeating


end module mod_particleInjection