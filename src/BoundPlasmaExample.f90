program BoundPlasmaExample
    use iso_fortran_env, only: int32, real64, output_unit
    use constants
    use mod_BasicFunctions
    use mod_domain
    use mod_particle
    use mod_potentialSolver
    use mod_collisions
    use mod_simulation
    implicit none

    integer(int32) :: particleIdxFactor = 2, i, numParticles = 10000, num_grid_nodes = 64
    real(real64) :: w_p = 1.0d0, n_ave = 5.0d14, T_e = 5.0d0, T_i = 0.1d0, T, L_domain = 0.1d0, del_l = 0.005d0
    type(Domain) :: world
    type(Particle) :: particleList(2)
    type(potSolver) :: solver
    
    n_x = num_grid_nodes
    ! create the world the particles live in
    world = Domain()
    call world % constructSineGrid(del_l, L_domain)

    !initialize the particles in this world, at some point will be read from input file or something
    particleList(1) = Particle(m_e, -e, w_p, numParticles, numParticles * particleIdxFactor, "e")
    particleList(2) = Particle(m_p, e, w_p, numParticles, numParticles * particleIdxFactor, "H+")
    do i = 1, size(particleList)
        print *, 'Initializing ', particleList(i) % name
        call particleList(i) % initialize_randUniform(L_domain, world%dx_dl, irand)
        call particleList(i) % initialize_n_ave(n_ave, L_domain)
        if (i == 1) then
            T = T_e
        else
            T = T_i
        end if
        call particleList(i) % generate3DMaxwellian(T, irand)
    end do

    print *, "w_p is:", particleList(1)%w_p
    print *, "Debye length is:", getDebyeLength(T_e, n_ave)
    print *, "Plasma frequency is:", getPlasmaFreq(n_ave)
    del_t = fractionFreq/getPlasmaFreq(n_ave)   
    print *, "Time step (sec) is:", del_t
    print *, "Mean initial KE of electron is:", particleList(1)%getKEAve(), "should be", T_e * 1.5
    print *, "Mean initial KE of proton is:", particleList(2)%getKEAve(), "should be", T_i * 1.5
    ! ! Generate solver object, and then solve for initial rho/potential
    solver = potSolver(world)
    call solver%solveInitialPotential(particleList, world)
    
    numTimeSteps = NINT(22.0d-6 / del_t)
    !call solver%solveDivAmperePicard(particleList, world, del_t, maxIter, eps_r, .false.)
    call solveSimulation(solver, particleList, world, del_t, maxIter, eps_r, irand, numTimeSteps)
    

    

    
end program BoundPlasmaExample