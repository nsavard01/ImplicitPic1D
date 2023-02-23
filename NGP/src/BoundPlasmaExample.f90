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
    real(real64) :: elapsed_time, KE_i, KE_f
    type(Domain) :: world
    type(Particle), allocatable :: particleList(:)
    type(potentialSolver) :: solver

    !initialize the particles in this world, at some point will be read from input file or something
    particleList = readParticleInputs('BoundExample.dat',numberChargedParticles, irand) 
    ! Initialize constants with inputs
    ! create the world the particles live in
    call readInputs(NumberXNodes, maxIter, numDiagnosticSteps, stepsAverage, eps_r, fractionFreq, n_ave, world, solver, simulationTime)
    do i = 1, numberChargedParticles
        call particleList(i) % initialize_randUniform(world%grid(NumberXNodes) - world%grid(1), world%dx_dl, irand)
        call particleList(i) % initialize_n_ave(n_ave, world%grid(NumberXNodes) - world%grid(1))
    end do

    print *, "Calulated values:"
    print *, "Number of particles is:", particleList(1)%N_p
    print *, "w_p is:", particleList(1)%w_p
    print *, "Debye length is:", getDebyeLength(particleList(1)%getKEAve()*2.0d0/3.0d0, n_ave)
    print *, "Plasma frequency is:", getPlasmaFreq(n_ave)
    print *, "Average density is ", particleList(1)%N_p * particleList(1)%w_p / (world%grid(NumberXNodes) - world%grid(1)), "should be", n_ave
    del_t = fractionFreq/getPlasmaFreq(n_ave)   
    print *, "Time step (sec) is:", del_t
    print *, "----------------"
    print *, ""
    ! ! Generate solver object, and then solve for initial rho/potential
    
    call solver%solveInitialPotential(particleList, world)
    numTimeSteps = NINT(simulationTime / del_t)
    
    stop
    call system_clock(tclock1)
    call solveSimulationOnlyPotential(solver, particleList, world, del_t, maxIter, eps_r, numTimeSteps)
    call system_clock(tclock2, clock_rate)
    elapsed_time = float(tclock2 - tclock1) / float(clock_rate)
    print *, "Elapsed time for simulation is:", elapsed_time/60.0d0, "minutes"
    print *, "Averaging over", stepsAverage, "time steps"
    !call solveSimulationFinalAverage(solver, particleList, world, del_t, maxIter, eps_r, irand, stepsAverage, heatSkipSteps)
    

    

    
end program BoundPlasmaExample