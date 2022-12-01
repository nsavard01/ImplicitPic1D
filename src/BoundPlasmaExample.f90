program BoundPlasmaExample
    use constants
    use mod_BasicFunctions
    use iso_fortran_env, only: int32, real64, output_unit
    use mod_domain
    use mod_particle
    implicit none

    integer(int32) :: particleIdxFactor = 2, i, irand = 12345, tclock1, tclock2, clock_rate
    integer(int32), parameter :: num_grid_nodes = 32, numParticles = 10000, maxIter = 50
    real(real64), parameter :: L_domain = 0.1, del_l = 0.005
    real(real64) :: w_p = 1.0, n_ave = 5e14, T_e = 5.0, T_i = 0.025, T, R_test(100000), t1, t2, elapsed_cpu_time, elapsed_time
    type(Domain) :: world
    type(Particle) :: particleList(2)

    R_test = 0
    ! create the world the particles live in
    world = Domain(num_grid_nodes)
    call world % constructSineGrid(del_l, L_domain)

    ! initialize the particles in this world
    particleList(1) = Particle(m_e, -e, w_p, numParticles, numParticles * particleIdxFactor, "electron")
    particleList(2) = Particle(m_p, e, w_p, numParticles, numParticles * particleIdxFactor, "proton")
    do i = 1, 2
        print *, particleList(i) % name
        call particleList(i) % initialize_randUniform(L_domain, world%dx_dl, world%n_x)
        call particleList(i) % initialize_n_ave(n_ave, L_domain)
        if (i == 1) then
            T = T_e
        else
            T = T_i
        end if
        call particleList(i) % generate3DMaxwellian(T)
        print *, particleList(i) % w_p * particleList(i) % N_p / L_domain
    end do

    print *, "Debye length is:", getDebyeLength(T_e, n_ave)
    print *, "Plasma frequency is:", getPlasmaFreq(n_ave)
    print *, "Mean temperature of electron is:", particleList(1)%getTemperature()
    print *, "Mean temperature of proton is:", particleList(2)%getTemperature()

    call cpu_time(t1)
    call system_clock(tclock1)
    do i = 1, 100
        call random_number(R_test)
    end do
    call system_clock(tclock2, clock_rate)
    call cpu_time(t2)
    elapsed_cpu_time = t2 - t1
    elapsed_time = float(tclock2 - tclock1) / float(clock_rate)

    print *, "Fortran random_number took ", elapsed_cpu_time, "CPU seconds"
    print *, "Fortran random_number took ", elapsed_time, "wall seconds"
    print *, "First few numbers in array are:"
    print *, R_test(1:10)
    print *, R_test(size(R_test)-10:size(R_test))
    print *, "Average is:", SUM(R_test)/size(R_test)

    call cpu_time(t1)
    call system_clock(tclock1)
    do i = 1, 100
        call getRandom(R_test, irand)
    end do
    call system_clock(tclock2, clock_rate)
    call cpu_time(t2)
    elapsed_cpu_time = t2 - t1
    elapsed_time = float(tclock2 - tclock1) / float(clock_rate)
    print *, "Gwenael random_number took ", elapsed_cpu_time, "CPU seconds"
    print *, "Gwenael random_number took ", elapsed_time, "wall seconds"
    print *, "First few numbers in array are:"
    print *, R_test(1:10)
    print *, R_test(size(R_test)-10:size(R_test))
    print *, "Average is:", SUM(R_test)/size(R_test)

end program BoundPlasmaExample