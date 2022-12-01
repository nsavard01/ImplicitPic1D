program BoundPlasmaExample
    use constants
    use mod_BasicFunctions
    use iso_fortran_env, only: int32, real64, output_unit
    use mod_domain
    use mod_particle
    use mod_TestClass
    implicit none

    integer(int32) :: particleIdxFactor = 2, i, irand = 9872364, tclock1, tclock2, clock_rate
    integer(int32), parameter :: num_grid_nodes = 32, numParticles = 100000, maxIter = 50
    real(real64), parameter :: L_domain = 0.1, del_l = 0.005
    real(real64) :: w_p = 1.0, n_ave = 5e14, T_e = 5.0, T_i = 0.025, T, R_test(100000), t1, t2, &
    elapsed_cpu_time, elapsed_time
    real(real64), allocatable :: v_test(:,:)
    type(Domain) :: world
    type(Particle) :: particleList(2)
    type(testClass) :: testing

    allocate(v_test(numParticles * particleIdxFactor, 3))
    v_test = 0
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
        call particleList(i) % generate3DMaxwellian(T, irand)
        print *, particleList(i) % w_p * particleList(i) % N_p / L_domain
    end do

    testing = testClass(particleList)
    print *, "Debye length is:", getDebyeLength(T_e, n_ave)
    print *, "Plasma frequency is:", getPlasmaFreq(n_ave)
    print *, "Mean temperature of electron is:", testing%particleList(1)%getTemperature(), "should be", 7.5
    print *, "Mean temperature of proton is:", testing%particleList(2)%getTemperature()


    ! Test timing for random number generator
    print *, "------------------- Test Random Number Generation -------------------"
    call cpu_time(t1)
    call system_clock(tclock1)
    do i = 1, 1000
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
    do i = 1, 1000
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

    ! Test maxwellian generation inside a particle class or outside
    print *, "--------------- Test objects vs functional speed -------------"

    call cpu_time(t1)
    call system_clock(tclock1)
    call testing % testFunction(1, T_e, irand)
    call system_clock(tclock2, clock_rate)
    call cpu_time(t2)
    elapsed_cpu_time = t2 - t1
    elapsed_time = float(tclock2 - tclock1) / float(clock_rate)

    print *, "OOP maxwell took ", elapsed_cpu_time, "CPU seconds"
    print *, "OOP maxwell took ", elapsed_time, "wall seconds"
    print *, "Mean temperature is:", testing%particleList(1)%getTemperature(), "should be ", 7.5
    
    call cpu_time(t1)
    call system_clock(tclock1)
    do i = 1, 1000
        call test3DMaxwellian(v_test, T_e, irand, numParticles, m_e)
    end do
    call system_clock(tclock2, clock_rate)
    call cpu_time(t2)
    elapsed_cpu_time = t2 - t1
    elapsed_time = float(tclock2 - tclock1) / float(clock_rate)

    print *, "Functional maxwell took ", elapsed_cpu_time, "CPU seconds"
    print *, "Functional maxwell took ", elapsed_time, "wall seconds"
    print *, "Mean temperature is:", testTemperature(v_test, m_e, numParticles), "should be ", 7.5

    


contains

    subroutine test3DMaxwellian(v, T, irand, v_size, mass)
        ! random velocity generator for the particle for temperature T (eV)
        ! Use box-muller method for random guassian variable, same as gwenael but doesn't have factor 2? Maybe factored into v_th
        real(real64), intent(in out) :: v(:,:)
        real(real64), intent(in) :: T, mass
        integer(int32), intent(in) :: v_size
        integer(int32), intent(in out) :: irand
        real(real64) :: U1(v_size), U2(v_size), U3(v_size), U4(v_size)
        call getRandom(U1, irand)
        call getRandom(U2, irand)
        call getRandom(U3, irand)
        call getRandom(U4, irand)
        v(1:v_size, 1) = SQRT(T*e/ mass) * SQRT(-2 * LOG(U1)) * COS(2 * pi * U2)
        v(1:v_size, 2) = SQRT(T*e/ mass) * SQRT(-2 * LOG(U1)) * SIN(2 * pi * U2)
        v(1:v_size, 3) = SQRT(T*e/ mass) * SQRT(-2 * LOG(U3)) * SIN(2 * pi * U4)
    end subroutine test3DMaxwellian

    pure function testTemperature(v, mass, N_p) result(res)
        ! calculate average kinetic energy (temperature) in eV
        real(real64), intent(in) :: v(:,:), mass
        integer(int32), intent(in) :: N_p
        real(real64) :: res
        res = SUM(v(1:N_p, :)**2) * mass * 0.5 / e / N_p

    end function testTemperature

end program BoundPlasmaExample