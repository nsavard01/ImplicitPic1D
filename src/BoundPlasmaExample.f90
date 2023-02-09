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

    integer(int32) :: i, tclock1, tclock2, clock_rate
    real(real64) :: T, elapsed_time
    type(Domain) :: world
    type(Particle), allocatable :: particleList(:)
    type(potentialSolver) :: solver

    !initialize the particles in this world, at some point will be read from input file or something
    particleList = readParticleInputs(numberChargedParticles) 
    ! Initialize constants with inputs
    ! create the world the particles live in
    call readInputs(NumberXNodes, maxIter, numDiagnosticSteps, stepsAverage, eps_r, fractionFreq, n_ave, T_e, T_i, L_domain, del_l, world, solver)
    
    do i = 1, numberChargedParticles
        print *, 'Initializing ', particleList(i) % name
        print *, "Particle mass is:", particleList(i)%mass
        print *, "Particle charge is:", particleList(i)%q
        call particleList(i) % initialize_randUniform(L_domain, world%dx_dl, irand)
        call particleList(i) % initialize_n_ave(n_ave, L_domain)
        if (i == 1) then
            T = T_e
        else
            T = T_i
        end if
        call particleList(i) % generate3DMaxwellian(T, irand)
    end do
    print *, "w_p is:", particleList(2)%w_p
    print *, "Debye length is:", getDebyeLength(T_e, n_ave)
    print *, "Plasma frequency is:", getPlasmaFreq(n_ave)
    del_t = fractionFreq/getPlasmaFreq(n_ave)   
    print *, "Time step (sec) is:", del_t
    print *, "Mean initial KE of electron is:", particleList(1)%getKEAve(), "should be", T_e * 1.5
    print *, "Mean initial KE of proton is:", particleList(2)%getKEAve(), "should be", T_i * 1.5
    ! ! Generate solver object, and then solve for initial rho/potential
    
    call solver%solveInitialPotential(particleList, world)
    numTimeSteps = NINT(22.0d-6 / del_t)
    call system_clock(tclock1)
    call solveSimulation(solver, particleList, world, del_t, maxIter, eps_r, irand, numTimeSteps, stepsAverage)
    call system_clock(tclock2, clock_rate)
    elapsed_time = float(tclock2 - tclock1) / float(clock_rate)
    print *, "Elapsed time for simulation is:", elapsed_time/60.0d0, "minutes"
    

    

    
end program BoundPlasmaExample