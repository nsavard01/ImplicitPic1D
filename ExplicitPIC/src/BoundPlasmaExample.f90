program BoundPlasmaExample
    use iso_fortran_env, only: int32, real64, output_unit
    use constants
    use mod_BasicFunctions
    use mod_domain
    use mod_particle
    use mod_collisions
    use mod_potentialSolver
    use mod_simulation
    implicit none

    integer(int32) :: i
    real(real64) :: KE_i, KE_f
    type(Domain) :: world
    type(Particle), allocatable :: particleList(:)
    type(potentialSolver) :: solver
    
    particleList = readParticleInputs('BoundExample.dat',numberChargedParticles, irand) 
    ! Initialize constants with inputs
    ! create the world the particles live in
    call readInputs(particleList(1), NumberXNodes, numDiagnosticSteps, averagingTime, fractionFreq, n_ave, world, solver, simulationTime, heatSkipSteps, Power, nu_h)
    do i = 1, numberChargedParticles
        call initialize_randUniform(particleList(i), irand)
        call particleList(i) % initialize_n_ave(n_ave, world%grid(NumberXNodes) - world%grid(1))
    end do

    print *, "Calulated values:"
    print *, "Number of particles is:", particleList(1)%N_p
    print *, "w_p is:", particleList(1)%w_p
    print *, "Debye length is:", getDebyeLength(particleList(1)%getKEAve()*2.0d0/3.0d0, n_ave)
    print *, "delX is:", world%delX
    print *, "Plasma frequency is:", getPlasmaFreq(n_ave)
    print *, "Average density is ", particleList(1)%N_p * particleList(1)%w_p / (world%grid(NumberXNodes) - world%grid(1)), "should be", n_ave
    del_t = fractionFreq/getPlasmaFreq(n_ave)   
    print *, "Time step (sec) is:", del_t
    print *, "----------------"
    print *, ""
    
    ! Generate solver object, and then solve for initial rho/potential
    call solveSimulation(solver, particleList, world, del_t, irand, simulationTime, heatSkipSteps)
    print *, "Averaging over", averagingTime, "seconds"
    call solveSimulationFinalAverage(solver, particleList, world, del_t, irand, averagingTime, heatSkipSteps)


    

    

    
end program BoundPlasmaExample