program BoundPlasmaExample
    use iso_fortran_env, only: int32, real64, output_unit
    use constants
    use mod_BasicFunctions
    use mod_domain
    use mod_particle
    use mod_potentialSolver
    use mod_collisions
    implicit none

    integer(int32) :: particleIdxFactor = 2, i, numParticles = 20000, num_grid_nodes = 64, maxIter = 50, tclock1, tclock2, clock_rate
    real(real64) :: w_p = 1.0d0, n_ave = 5d14, T_e = 5.0d0, T_i = 0.025d0, T, del_t, fractionFreq = 0.5d0, L_domain = 0.1d0, del_l = 0.005d0, eps_a = 1e-8, elapsed_time
    real(real64) :: KE_e_i
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
    print *, "Mean KE of electron is:", particleList(1)%getKEAve(), "should be", T_e * 1.5
    print *, "Mean KE of proton is:", particleList(2)%getKEAve(), "should be", T_i * 1.5
    solver = potSolver(world)
    call solver%depositRho(particleList, world)
    call solver%solve_tridiag_Poisson(world)
    call solver%construct_diagMatrix_Ampere(world)

    call system_clock(tclock1)
    call solver%solveDivAmperePicard(particleList, world, del_t, maxIter, eps_a, .true.)
    call system_clock(tclock2, clock_rate)
    elapsed_time = float(tclock2 - tclock1) / float(clock_rate)
    print *, "Elapsed time is:", elapsed_time, "seconds"
    print *, "Energy error is:", solver%energyError
    print *, "Charge error is:", solver%chargeError

    KE_e_i = particleList(1)%getTotalKE()
    print *, "KE sum of electron before adding power is:", KE_e_i
    call addUniformPowerMaxwellianNicolas(particleList(1), Power, 0.05d0, irand, del_t)

    print *, "KE sum of electron after adding power is:", particleList(1)%getTotalKE()
    print *, "Energy loss is:", Power*del_t
    print *, "Loss of electron energy is:", particleList(1)%getTotalKE() - KE_e_i

    !call ionizationCollisionIsotropic(particleList(1), particleList(2), 1.0d20, 1.0d-20, del_t, 15.8d0, irand, 300.0d0 * k_B/e)

    



    

    



    
    

    

    





end program BoundPlasmaExample