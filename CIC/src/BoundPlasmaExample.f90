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
    real(real64) :: elapsed_time, KE_i, KE_f, PE_i, PE_f
    type(Domain) :: world
    type(Particle), allocatable :: particleList(:)
    type(potentialSolver) :: solver
    
    !initialize the particles in this world, at some point will be read from input file or something
    particleList = readParticleInputs('BoundExample.dat',numberChargedParticles, irand) 
    ! Initialize constants with inputs
    ! create the world the particles live in
    call readInputs(NumberXNodes, maxIter, numDiagnosticSteps, stepsAverage, eps_r, fractionFreq, n_ave, world, solver)
    do i = 1, numberChargedParticles
        call particleList(i) % initialize_randUniform(world, irand)
        call particleList(i) % initialize_n_ave(n_ave, world%grid(NumberXNodes) - world%grid(1))
    end do

    print *, "Calulated values:"
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
    inelasticEnergyLoss = 0.0d0
    PE_i = SUM(particleList(1)%phaseSpace(2:4, 1:particleList(1)%N_p)) * m_e + SUM(particleList(2)%phaseSpace(2:4, 1:particleList(2)%N_p)) * particleList(2)%mass
    KE_i = particleList(1)%getTotalKE() + particleList(2)%getTotalKE()
    call ionizationCollisionIsotropic(particleList(1), particleList(2), 1.0d20, 1.0d-20, del_t*4.0d0, 15.8d0, 0.0d0, irand)
    KE_f = particleList(1)%getTotalKE() + particleList(2)%getTotalKE()
    PE_f = SUM(particleList(1)%phaseSpace(2:4, 1:particleList(1)%N_p)) * m_e + SUM(particleList(2)%phaseSpace(2:4, 1:particleList(2)%N_p)) * particleList(2)%mass
    print *, "Difference in energy is:"
    print *, ABS((inelasticEnergyLoss - (KE_i - KE_f)) / (KE_i - KE_f))
    print *, "Difference in momentum is:"
    print *, ABS((PE_i - PE_f)/PE_i)
    numTimeSteps = NINT(22.0d-6 / del_t)
    stop
    call system_clock(tclock1)
    call solveSimulation(solver, particleList, world, del_t, maxIter, eps_r, irand, numTimeSteps, 4)
    call system_clock(tclock2, clock_rate)
    elapsed_time = float(tclock2 - tclock1) / float(clock_rate)
    print *, "Elapsed time for simulation is:", elapsed_time/60.0d0, "minutes"
    

    

    
end program BoundPlasmaExample