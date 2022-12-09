program BoundPlasmaExample
    use constants
    use mod_BasicFunctions
    use iso_fortran_env, only: int32, real64, output_unit
    use mod_domain
    use mod_particle
    use mod_potentialSolver
    implicit none

    integer(int32) :: particleIdxFactor = 2, i, irand = 9872364, tclock1, tclock2, clock_rate
    integer(int32), parameter :: num_grid_nodes = 64, numParticles = 200000, maxIter = 50
    real(real64), parameter :: L_domain = 0.1d0, del_l = 0.005d0
    real(real64) :: w_p = 1.0d0, n_ave = 5d14, T_e = 5.0d0, T_i = 0.025d0, T, del_t, fractionFreq = 0.5d0, elapsed_time
    type(Domain) :: world
    type(Particle) :: particleList(2)
    type(potSolver) :: solver
    

    ! create the world the particles live in
    world = Domain(num_grid_nodes)
    call world % constructSineGrid(del_l, L_domain)
    ! initialize the particles in this world, at some point will be read from input file or something
    particleList(1) = Particle(m_e, -e, w_p, numParticles, numParticles * particleIdxFactor, "electron")
    particleList(2) = Particle(m_p, e, w_p, numParticles, numParticles * particleIdxFactor, "proton")
    do i = 1, 2
        print *, 'Initializing ', particleList(i) % name
        call particleList(i) % initialize_randUniform(L_domain, world%dx_dl, world%n_x)
        call particleList(i) % initialize_n_ave(n_ave, L_domain)
        if (i == 1) then
            T = T_e
        else
            T = T_i
        end if
        call particleList(i) % generate3DMaxwellian(T, irand)
    end do

    
    print *, "Debye length is:", getDebyeLength(T_e, n_ave)
    print *, "Plasma frequency is:", getPlasmaFreq(n_ave)
    del_t = fractionFreq/getPlasmaFreq(n_ave)
    print *, "Time step (sec) is:", del_t
    print *, "Mean temperature of electron is:", particleList(1)%getTemperature(), "should be", T_e * 1.5
    print *, "Mean temperature of proton is:", particleList(2)%getTemperature(), "should be", T_i * 1.5
    solver = potSolver(num_grid_nodes, world)
    call solver%depositRho(particleList(1:1), world)
    call solver%solve_tridiag_Poisson()

    call system_clock(tclock1)
    call solver%depositJ(particleList(1:1), world, del_t)
    call system_clock(tclock2, clock_rate)
    elapsed_time = float(tclock2 - tclock1) / float(clock_rate)
    print *, "Elapsed time is:", elapsed_time, "seconds"

    

    





end program BoundPlasmaExample