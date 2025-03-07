module mod_particleInjection
    use iso_fortran_env, only: int32, real64
    use constants
    use mod_BasicFunctions
    use mod_domain
    use mod_particle
    use iso_c_binding
    use omp_lib
    implicit none
    ! Module for artificial particle injection at boundaries or artificial collisions
    real(real64), allocatable, public :: energyAddColl(:)
    logical, public, protected :: addLostPartBool, refluxPartBool, injectionBool, uniformInjectionBool, heatingBool, EField_heating_bool
    integer(int32), private :: numFluxParticlesHigh, numFluxParticlesLow
    real(real64), protected :: injectionFlux, injectionR, FractionFreqHeating, Efield_heating_freq, J_total_RF, power_left_bound_x, power_right_bound_x, Efield_RF_past, Efield_RF_energy, J_particles_heat, v_ave_heat
    integer(int32) :: Number_RF_heat_electrons
    real(real64) :: v_sum_RF_electrons, del_v_RF_electrons
    real(real64) :: current_RF_time = 0.0

contains

    subroutine addMaxwellianLostParticles(particleList, T_e, T_i, irand, world)
        ! Add power to all particles in domain
        type(Particle), intent(in out) :: particleList(2)
        type(Domain), intent(in) :: world
        integer(c_int64_t), intent(in out) :: irand(numThread)
        real(real64), intent(in) :: T_e, T_i
        integer(int32) :: i, iThread
        integer(c_int64_t) :: irand_thread
        real(real64) :: x_random
        !$OMP parallel private(iThread, i, x_random, irand_thread)
        iThread = omp_get_thread_num() + 1
        irand_thread = irand(iThread)
        do i=1, particleList(1)%delIdx(iThread)
            call getMaxwellianSample(particleList(1)%phaseSpace(2:4, particleList(1)%N_p(iThread) + i, iThread), particleList(1)%mass, T_e, irand_thread)
            call getMaxwellianFluxSample(particleList(2)%phaseSpace(2:4, particleList(2)%N_p(iThread) + i, iThread), particleList(2)%mass, T_i, irand_thread)
            energyAddColl(iThread) = energyAddColl(iThread) + (SUM(particleList(1)%phaseSpace(2:4, particleList(1)%N_p(iThread) + i, iThread)**2) * particleList(1)%mass * particleList(1)%w_p &
            + SUM(particleList(2)%phaseSpace(2:4, particleList(2)%N_p(iThread) + i, iThread)**2) * particleList(2)%mass * particleList(2)%w_p) * 0.5d0
            x_random = pcg32_random_r(irand_thread) * world%L_domain + world%startX
            particleList(1)%phaseSpace(1, particleList(1)%N_p(iThread) + i, iThread) = world%getLFromX(x_random)
            particleList(2)%phaseSpace(1, particleList(2)%N_p(iThread) + i, iThread) = particleList(1)%phaseSpace(1, particleList(1)%N_p(iThread) + i, iThread)
        end do
        particleList(2)%N_p(iThread) = particleList(2)%N_p(iThread) + particleList(1)%delIdx(iThread)
        particleList(1)%N_p(iThread) = particleList(1)%N_p(iThread) + particleList(1)%delIdx(iThread)
        irand(iThread) = irand_thread
        !$OMP end parallel
    end subroutine addMaxwellianLostParticles

    subroutine refluxParticles(particleList, T_e, T_i, irand, world)
        !reflux particles at the boundary
        type(Particle), intent(in out) :: particleList(2)
        type(Domain), intent(in) :: world
        integer(c_int64_t), intent(in out) :: irand(numThread)
        real(real64), intent(in) :: T_e, T_i
        integer(int32) :: i,j, iThread
        integer(c_int64_t) :: irand_thread
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

    subroutine injectAtBoundary(particleList, T_e, T_i, irand, world, del_t)
        ! inject particles at boundary at ion sound speed
        type(Particle), intent(in out) :: particleList(2)
        type(Domain), intent(in) :: world
        integer(c_int64_t), intent(in out) :: irand(numThread)
        real(real64), intent(in) :: T_e, T_i, del_t
        integer(int32) :: i,j, iThread, numInjected
        real(real64) :: ionSoundSpeed, shift_rand
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
        if (pcg32_random_r(irand(iThread)) > injectionR) then
            numInjected = numFluxParticlesLow
        else
            numInjected = numFluxParticlesHigh
        end if
        do i = 1, numInjected
            shift_rand = del_t * pcg32_random_r(irand(iThread)) * ionSoundSpeed
            call getMaxwellianFluxSample(particleList(1)%phaseSpace(2:4, particleList(1)%N_p(iThread) + i, iThread), particleList(1)%mass, T_e, irand(iThread))
            particleList(2)%phaseSpace(2, particleList(2)%N_p(iThread) + i, iThread) = ionSoundSpeed
            particleList(2)%phaseSpace(3, particleList(2)%N_p(iThread) + i, iThread) = 0.0d0
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
        ! volume injection at set rate
        type(Particle), intent(in out) :: particleList(2)
        type(Domain), intent(in) :: world
        integer(c_int64_t), intent(in out) :: irand(numThread)
        real(real64), intent(in) :: T_e, T_i
        integer(int32) :: i, iThread, numInjected
        real(real64) :: x_random
        
        !$OMP parallel private(iThread, i, numInjected, x_random)
        iThread = omp_get_thread_num() + 1
        if (pcg32_random_r(irand(iThread)) > injectionR) then
            numInjected = numFluxParticlesLow
        else
            numInjected = numFluxParticlesHigh
        end if
        do i = 1, numInjected
            call getMaxwellianSample(particleList(1)%phaseSpace(2:4, particleList(1)%N_p(iThread) + i, iThread), particleList(1)%mass, 20.0d0 * T_e, irand(iThread))
            call getMaxwellianSample(particleList(2)%phaseSpace(2:4, particleList(2)%N_p(iThread) + i, iThread), particleList(2)%mass, T_i, irand(iThread))
            x_random = pcg32_random_r(irand(iThread)) * world%L_domain + world%startX
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
        integer(c_int64_t), intent(in out) :: irand(numThread)
        real(real64), intent(in) :: FracFreq, T, currDel_t, del_t
        real(real64) :: R, compValue
        integer(int32) :: i, iThread
        compValue = FracFreq * (currDel_t/del_t)/FractionFreqHeating
        !$OMP parallel private(iThread, i, R)
        iThread = omp_get_thread_num() + 1
        do i = 1, species%N_p(iThread)
            R = pcg32_random_r(irand(iThread))
            if (R < compValue) then
                call getMaxwellianSample(species%phaseSpace(2:4, i, iThread), species%mass, T, irand(iThread))
            end if
        end do
        !$OMP end parallel
    end subroutine maxwellianHeating

    subroutine EField_heating(electron, old_time, current_time, world)
        ! Add velocity along y direction to particles from given oscillating field (sinusoidal)
        type(Particle), intent(in out) :: electron
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: current_time, old_time
        real(real64) :: del_v, coeff_time, xi, Efield_future, v_sum, del_t, xi_left, xi_right
        integer(int32) :: i, iThread, number_part
        
        coeff_time = 2.0d0 * pi * Efield_heating_freq

        ! try leap-frog to make it easier
        xi_left = world%getLfromX(power_left_bound_x)
        xi_right = world%getLfromX(power_right_bound_x)
        del_t = current_time - old_time
        del_v = electron%q_over_m * Efield_RF_past * del_t
        v_sum = 0
        number_part = 0
        !$OMP parallel private(iThread, i) reduction(+:v_sum, number_part)
        iThread = omp_get_thread_num() + 1
        do i = 1, electron%N_p(iThread)
            xi = electron%phaseSpace(1, i, iThread)
            if (xi > xi_left .and. xi < xi_right) then
                electron%phaseSpace(3, i, iThread) = electron%phaseSpace(3, i, iThread) + del_v
                v_sum = v_sum + electron%phaseSpace(3, i, iThread)
                number_part = number_part + 1
            end if
        end do
        !$OMP end parallel

        
        Efield_future = J_total_RF * (cos(coeff_time * old_time) - cos(coeff_time * current_time)) / coeff_time / eps_0 - electron%q_times_wp * v_sum * del_t / (power_right_bound_x - power_left_bound_x) /eps_0 + Efield_RF_past
       
      
        J_particles_heat = electron%q_times_wp * v_sum / (power_right_bound_x - power_left_bound_x)
        
        v_ave_heat = v_sum / real(number_part, kind = 8)
        Efield_RF_energy = (v_sum * electron%q * Efield_RF_past * del_t + real(number_part, kind = 8) * 0.5d0 * electron%q**2 * Efield_Rf_past**2 * del_t**2 / electron%mass) ! Joules
   
        Efield_RF_past = Efield_future
        
    end subroutine EField_heating

    subroutine EField_heating_implicit(electron, del_t, world)
        ! Add velocity along y direction to particles from given oscillating field (sinusoidal)
        type(Particle), intent(in out) :: electron
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t
        real(real64) :: coeff_time, Efield_future, EField_half, other_coeff, new_time
        
        coeff_time = 2.0d0 * pi * Efield_heating_freq

        ! try leap-frog to make it easier
        ! xi_left = world%getLfromX(power_left_bound_x)
        ! xi_right = world%getLfromX(power_right_bound_x)
        new_time = current_RF_time + del_t
        ! v_sum = 0
        ! number_part = 0
        ! !$OMP parallel private(iThread, i) reduction(+:v_sum, number_part)
        ! iThread = omp_get_thread_num() + 1
        ! do i = 1, electron%N_p(iThread)
        !     xi = electron%phaseSpace(1, i, iThread)
        !     if (xi > xi_left .and. xi < xi_right) then
        !         v_sum = v_sum + electron%phaseSpace(3, i, iThread)
        !         number_part = number_part + 1
        !     end if
        ! end do
        ! !$OMP end parallel

        other_coeff = 0.25d0 * electron%q * electron%q_times_wp * real(Number_RF_heat_electrons, kind = 8) * del_t**2 / (power_right_bound_x - power_left_bound_x) / electron%mass / eps_0
        Efield_future = J_total_RF * (cos(coeff_time * current_RF_time) - cos(coeff_time * new_time)) / coeff_time / eps_0 - &
            electron%q_times_wp * v_sum_RF_electrons * del_t / (power_right_bound_x - power_left_bound_x) /eps_0 + &
                (1.0d0 - other_coeff) * Efield_RF_past
        Efield_future = Efield_future / (1.0d0 + other_coeff)
                
        EField_half = 0.5d0 * (Efield_future + Efield_RF_past)
        del_v_RF_electrons = electron%q_over_m * Efield_half * del_t
        
        ! !$OMP parallel private(iThread, i)
        ! iThread = omp_get_thread_num() + 1
        ! do i = 1, electron%N_p(iThread)
        !     xi = electron%phaseSpace(1, i, iThread)
        !     if (xi > xi_left .and. xi < xi_right) then
        !         electron%phaseSpace(3, i, iThread) = electron%phaseSpace(3, i, iThread) + del_v
        !     end if
        ! end do
        ! !$OMP end parallel

        J_particles_heat = electron%q_times_wp * (v_sum_RF_electrons + del_v_RF_electrons) / (power_right_bound_x - power_left_bound_x)
        
        v_ave_heat = (v_sum_RF_electrons + del_v_RF_electrons) / real(Number_RF_heat_electrons, kind = 8)
        Efield_RF_energy = 0.5d0 * electron%mass * (2 * del_v_RF_electrons * v_sum_RF_electrons + real(Number_RF_heat_electrons, kind = 8) * del_v_RF_electrons**2) ! Joules
   
        Efield_RF_past = Efield_future
        current_RF_time = new_time
        
    end subroutine EField_heating_implicit


    ! -------------------- read inputs artificial collisions/particle injection ------------------------------------

    subroutine readInjectionInputs(InjFilename, w_p)
        ! Read input for artificial particle collisions
        real(real64), intent(in) :: w_p
        character(len=*), intent(in) :: InjFilename
        integer(int32) :: tempInt, io
        real(real64) :: numFluxPart
        print *, ""
        print *, "Reading initial inputs for particle injection:"
        print *, "------------------"
        if (.not. restartBool) then
            open(10,file='../InputData/'//InjFilename, IOSTAT=io)
        else
            open(10,file=restartDirectory//"/InputData/"//InjFilename, IOSTAT=io)
        end if
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
        read(10, *, IOSTAT = io) Efield_heating_freq, J_total_RF, power_left_bound_x, power_right_bound_x
        if (EField_heating_freq /= 0) then
            EField_heating_bool = .true.
        else
            EField_heating_bool = .false.
        end if
        close(10)
        print *, "Particle lost is reinjected:", addLostPartBool
        print *, "Particle refluxing activated on neumann boundary:", refluxPartBool
        print *, "Particle injection on neumann boundary", injectionBool
        print *, "Particle injection unformly with maxwellian", uniformInjectionBool
        print *, "Electron maxwellian heating", heatingBool
        print *, 'Electron E-field heating', EField_heating_bool
        if (injectionBool) then
            injectionFlux = injectionFlux
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
        if (EField_heating_bool) then
            print *, 'EField frequency:', EField_heating_freq, 'Hz'
            print *, 'J_total current density:', J_total_RF, 'A/m^2'
            print *, 'Heating region boundary between', power_left_bound_x, 'and', power_right_bound_x, ' meters'
            Efield_RF_past = 0.0d0
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