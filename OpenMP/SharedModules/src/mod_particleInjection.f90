module mod_particleInjection
    use iso_fortran_env, only: int32, real64
    use constants
    use mod_BasicFunctions
    use mod_domain
    use mod_particle
    use omp_lib
    implicit none
    real(real64), allocatable, public :: energyAddColl(:)
    logical, public, protected :: addLostPartBool, refluxPartBool, injectionBool, uniformInjectionBool, heatingBool
    integer(int32), private :: numFluxParticlesHigh, numFluxParticlesLow
    real(real64), private :: injectionFlux, injectionR, FractionFreqHeating

contains

    subroutine addMaxwellianLostParticles(particleList, T_e, T_i, irand, world)
        ! Add power to all particles in domain
        type(Particle), intent(in out) :: particleList(2)
        type(Domain), intent(in) :: world
        integer(int32), intent(in out) :: irand(numThread)
        real(real64), intent(in) :: T_e, T_i
        integer(int32) :: i, iThread, irand_thread
        real(real64) :: x_random
        !$OMP parallel private(iThread, i, x_random, irand_thread)
        iThread = omp_get_thread_num() + 1
        irand_thread = irand(iThread)
        do i=1, particleList(1)%delIdx(iThread)
            call getMaxwellianSample(particleList(1)%phaseSpace(2:4, particleList(1)%N_p(iThread) + i, iThread), particleList(1)%mass, T_e, irand_thread)
            call getMaxwellianFluxSample(particleList(2)%phaseSpace(2:4, particleList(2)%N_p(iThread) + i, iThread), particleList(2)%mass, T_i, irand_thread)
            energyAddColl(iThread) = energyAddColl(iThread) + (SUM(particleList(1)%phaseSpace(2:4, particleList(1)%N_p(iThread) + i, iThread)**2) * particleList(1)%mass * particleList(1)%w_p &
            + SUM(particleList(2)%phaseSpace(2:4, particleList(2)%N_p(iThread) + i, iThread)**2) * particleList(2)%mass * particleList(2)%w_p) * 0.5d0
            x_random = ran2(irand_thread) * world%L_domain + world%startX
            particleList(1)%phaseSpace(1, particleList(1)%N_p(iThread) + i, iThread) = world%getLFromX(x_random)
            particleList(2)%phaseSpace(1, particleList(2)%N_p(iThread) + i, iThread) = particleList(1)%phaseSpace(1, particleList(1)%N_p(iThread) + i, iThread)
        end do
        particleList(2)%N_p(iThread) = particleList(2)%N_p(iThread) + particleList(1)%delIdx(iThread)
        particleList(1)%N_p(iThread) = particleList(1)%N_p(iThread) + particleList(1)%delIdx(iThread)
        irand(iThread) = irand_thread
        !$OMP end parallel
    end subroutine addMaxwellianLostParticles

    subroutine refluxParticles(particleList, T_e, T_i, irand, world)
        type(Particle), intent(in out) :: particleList(2)
        type(Domain), intent(in) :: world
        integer(int32), intent(in out) :: irand(numThread)
        real(real64), intent(in) :: T_e, T_i
        integer(int32) :: i,j, iThread, irand_thread
        do j = 1, 2
            !$OMP parallel private(iThread, i, irand_thread)
            iThread = omp_get_thread_num() + 1
            irand_thread = irand(iThread)
            do i = 1, particleList(j)%refIdx(iThread)
                energyAddColl(iThread) = energyAddColl(iThread) - (SUM(particleList(j)%phaseSpace(2:4, particleList(j)%refRecordIdx(i, iThread), iThread)**2) * particleList(j)%mass * particleList(j)%w_p) * 0.5d0
                if (j == 1) then
                    call getMaxwellianFluxSample(particleList(j)%phaseSpace(2:4, particleList(j)%refRecordIdx(i, iThread), iThread), particleList(j)%mass, T_e, irand_thread)
                else
                    call getMaxwellianFluxSample(particleList(j)%phaseSpace(2:4, particleList(j)%refRecordIdx(i, iThread), iThread), particleList(j)%mass, T_i, irand_thread)
                end if
                if (world%boundaryConditions(1) == 2) then
                    particleList(j)%phaseSpace(2, particleList(j)%refRecordIdx(i, iThread), iThread) = ABS(particleList(j)%phaseSpace(2, particleList(j)%refRecordIdx(i, iThread), iThread))
                else
                    particleList(j)%phaseSpace(2, particleList(j)%refRecordIdx(i, iThread), iThread) = -ABS(particleList(j)%phaseSpace(2, particleList(j)%refRecordIdx(i, iThread), iThread))
                end if
                energyAddColl(iThread) = energyAddColl(iThread) + (SUM(particleList(j)%phaseSpace(2:4, particleList(j)%refRecordIdx(i, iThread), iThread)**2) * particleList(j)%mass * particleList(j)%w_p) * 0.5d0
            end do
            irand(iThread) = irand_thread
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
            if (world%boundaryConditions(1) == 2) then
                particleList(1)%phaseSpace(2, particleList(1)%N_p(iThread) + i, iThread) = ABS(particleList(1)%phaseSpace(2, particleList(1)%N_p(iThread) + i, iThread))
                particleList(2)%phaseSpace(2, particleList(2)%N_p(iThread) + i, iThread) = ABS(particleList(2)%phaseSpace(2, particleList(2)%N_p(iThread) + i, iThread))
                particleList(1)%phaseSpace(1, particleList(1)%N_p(iThread) + i, iThread) = world%getLFromX(world%startX + shift_rand)
                particleList(2)%phaseSpace(1, particleList(2)%N_p(iThread) + i, iThread) = particleList(1)%phaseSpace(1, particleList(1)%N_p(iThread) + i, iThread)
            else
                particleList(1)%phaseSpace(2, particleList(1)%N_p(iThread) + i, iThread) = -ABS(particleList(1)%phaseSpace(2, particleList(1)%N_p(iThread) + i, iThread))
                particleList(2)%phaseSpace(2, particleList(2)%N_p(iThread) + i, iThread) = -ABS(particleList(2)%phaseSpace(2, particleList(2)%N_p(iThread) + i, iThread))
                particleList(1)%phaseSpace(1, particleList(1)%N_p(iThread) + i, iThread) = world%getLFromX(world%endX - shift_rand)
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
        real(real64) :: x_random
        ! Reflux particles
        
        !$OMP parallel private(iThread, i, numInjected, x_random)
        iThread = omp_get_thread_num() + 1
        if (ran2(irand(iThread)) > injectionR) then
            numInjected = numFluxParticlesLow
        else
            numInjected = numFluxParticlesHigh
        end if
        do i = 1, numInjected
            call getMaxwellianSample(particleList(1)%phaseSpace(2:4, particleList(1)%N_p(iThread) + i, iThread), particleList(1)%mass, 20.0d0 * T_e, irand(iThread))
            call getMaxwellianSample(particleList(2)%phaseSpace(2:4, particleList(2)%N_p(iThread) + i, iThread), particleList(2)%mass, T_i, irand(iThread))
            x_random = ran2(irand(iThread)) * world%L_domain + world%startX
            particleList(1)%phaseSpace(1, particleList(1)%N_p(iThread) + i, iThread) = world%getLFromX(x_random)
            particleList(2)%phaseSpace(1, particleList(2)%N_p(iThread) + i, iThread) = particleList(1)%phaseSpace(1, particleList(1)%N_p(iThread) + i, iThread)
        end do
        particleList(1)%N_p(iThread) = particleList(1)%N_p(iThread) + numInjected
        particleList(2)%N_p(iThread) = particleList(2)%N_p(iThread) + numInjected
        !$OMP end parallel
    end subroutine injectUniformFlux

    subroutine maxwellianHeating(species, irand, FracFreq, T, currDel_t, del_t)
        ! Add power to all particles in domain
        type(Particle), intent(in out) :: species
        integer(int32), intent(in out) :: irand(numThread)
        real(real64), intent(in) :: FracFreq, T, currDel_t, del_t
        real(real64) :: R, compValue
        integer(int32) :: i, iThread
        compValue = FracFreq * (currDel_t/del_t)/FractionFreqHeating
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


    ! -------------------- read inputs artificial collisions/particle injection ------------------------------------

    subroutine readInjectionInputs(InjFilename, w_p, angleRad)
        real(real64), intent(in) :: w_p, angleRad
        character(len=*), intent(in) :: InjFilename
        integer(int32) :: tempInt, io
        real(real64) :: numFluxPart
        print *, ""
        print *, "Reading initial inputs for particle injection:"
        print *, "------------------"
        if (.not. restartBool) then
            open(10,file='../InputData/'//InjFilename, IOSTAT=io)
            read(10, *, IOSTAT = io) tempInt
            addLostPartBool = (tempInt == 1)
            read(10, *, IOSTAT = io) tempInt
            refluxPartBool = (tempInt == 1)
            read(10, *, IOSTAT = io) tempInt, injectionFlux
            injectionBool = (tempInt == 1)
            if (.not. injectionBool) then
                read(10, *, IOSTAT = io) tempInt, injectionFlux
                uniformInjectionBool = (tempInt == 1)
            end if
            read(10, *, IOSTAT = io) tempInt, FractionFreqHeating
            heatingBool = (tempInt == 1)
            close(10)
        else
            open(10,file=restartDirectory//"/"//"ParticleInjectionInput.dat", IOSTAT=io)
            read(10, *, IOSTAT = io)
            read(10, *, IOSTAT = io) addLostPartBool, refluxPartBool, injectionBool, uniformInjectionBool, heatingBool, injectionFlux, FractionFreqHeating
            close(10)
        end if
        print *, "Particle lost is reinjected:", addLostPartBool
        print *, "Particle refluxing activated on neumann boundary:", refluxPartBool
        print *, "Particle injection on neumann boundary", injectionBool
        print *, "Particle injection unformly with maxwellian", uniformInjectionBool
        print *, "Electron maxwellian heating", heatingBool
        if (injectionBool) then
            injectionFlux = injectionFlux * COS(angleRad)
            print *, 'Particle injection flux:', injectionFlux
            numFluxPart = injectionFlux * del_t / w_p/real(numThread) ! particles injected per thread
            numFluxParticlesLow = floor(numFluxPart)
            numFluxParticlesHigh = numFluxParticlesLow + 1
            print *, 'Low end of flux particles:', numFluxParticlesLow
            print *, 'High end of flux particles:', numFluxParticlesHigh
            injectionR = numFluxPart - real(numFluxParticlesLow)
            print *, 'Number for selection of flux particles is:', injectionR
        else if (uniformInjectionBool) then
            print *, 'Particle injection flux:', injectionFlux
            numFluxPart = injectionFlux * del_t / w_p/real(numThread) ! particles injected per thread
            numFluxParticlesLow = floor(numFluxPart)
            numFluxParticlesHigh = numFluxParticlesLow + 1
            print *, 'Low end of flux particles:', numFluxParticlesLow
            print *, 'High end of flux particles:', numFluxParticlesHigh
            injectionR = numFluxPart - real(numFluxParticlesLow)
            print *, 'Number for selection of flux particles is:', injectionR
        else if (heatingBool) then
            print *, 'FractionFreqHeating for maxwellian heating is:', FractionFreqHeating
        end if
        print *, "------------------"
        print *, ""
    end subroutine readInjectionInputs

    subroutine writeParticleInjectionInputs(dirName)
        character(*), intent(in) :: dirName

        open(15,file=dirName//'/ParticleInjectionInput.dat')
        write(15,'("AddLostPartBool, refluxPartBool, injectionBool, uniformInjectionBool, heatingBool, injectionFlux, FractionFreqHeating")')
        write(15,"(5(L, 1x), 2(es16.8,1x))") addLostPartBool, refluxPartBool, injectionBool, uniformInjectionBool, heatingBool, injectionFlux, FractionFreqHeating
        close(15)
    end subroutine writeParticleInjectionInputs

end module mod_particleInjection