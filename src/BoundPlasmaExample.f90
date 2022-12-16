program BoundPlasmaExample
    use iso_fortran_env, only: int32, real64, output_unit
    use constants
    use mod_BasicFunctions
    use mod_domain
    use mod_particle
    use mod_potentialSolver
    implicit none

    integer(int32) :: particleIdxFactor = 2, i, irand = 9872364, numParticles = 10000, num_grid_nodes = 64, maxIter = 50, tclock1, tclock2, clock_rate
    real(real64) :: w_p = 1.0d0, n_ave = 5d14, T_e = 5.0d0, T_i = 0.025d0, T, del_t, fractionFreq = 0.5d0, L_domain = 0.1d0, del_l = 0.005d0, eps_a = 1e-8, elapsed_time
    type(Domain) :: world
    type(Particle) :: particleList(2)
    type(potSolver) :: solver
    
    n_x = num_grid_nodes
    ! create the world the particles live in
    world = Domain()
    call world % constructSineGrid(del_l, L_domain)
    !initialize the particles in this world, at some point will be read from input file or something
    particleList(1) = Particle(m_e, -e, w_p, numParticles, numParticles * particleIdxFactor, "electron")
    particleList(2) = Particle(m_p, e, w_p, numParticles, numParticles * particleIdxFactor, "proton")
    do i = 1, size(particleList)
        print *, 'Initializing ', particleList(i) % name
        call particleList(i) % initialize_randUniform(L_domain, world%dx_dl)
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
    solver = potSolver(world)
    call solver%depositRho(particleList, world)
    call solver%solve_tridiag_Poisson(world)
    call solver%construct_diagMatrix_Ampere(world)

    call system_clock(tclock1)
    call solver%solveDivAmperePicard(particleList, world, del_t, maxIter, eps_a)
    call system_clock(tclock2, clock_rate)
    elapsed_time = float(tclock2 - tclock1) / float(clock_rate)
    print *, "Elapsed time is:", elapsed_time, "seconds"

    

    



    
    

    

    





end program BoundPlasmaExample