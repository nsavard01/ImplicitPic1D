program BoundPlasmaExample
    use iso_fortran_env, only: int32, real64, output_unit, int64
    use constants
    use mod_BasicFunctions
    use mod_domain
    use mod_particle
    use mod_potentialSolver
    use mod_collisions
    use mod_nonLinSolvers
    use mod_Scheme
    use mod_simulation
    implicit none

    integer(int32) :: i
    
    globalParticleList = readParticleInputs('BoundExample.dat',numberChargedParticles, irand, T_e) 
    ! Initialize constants with inputs
    ! create the world the particles live in
    call readInputs(NumberXNodes, numDiagnosticSteps, averagingTime, fractionFreq, n_ave, globalWorld, globalSolver, simulationTime, Power, heatSkipSteps, nu_h, T_e)
    do i = 1, numberChargedParticles
        call initialize_randUniform(globalParticleList(i), globalWorld, irand)
        call globalParticleList(i) % initialize_n_ave(n_ave, globalWorld%grid(NumberXNodes) - globalWorld%grid(1))
    end do
    call initializeSolver(eps_r, solverType, m_Anderson, Beta_k, maxIter)
    print *, "Calulated values:"
    print *, "Number of particles is:", globalParticleList(1)%N_p
    print *, "w_p is:", globalParticleList(1)%w_p
    print *, "Debye length is:", getDebyeLength(globalParticleList(1)%getKEAve()*2.0d0/3.0d0, n_ave)
    print *, "Plasma frequency is:", getPlasmaFreq(n_ave)
    print *, "Average density is ", globalParticleList(1)%N_p * globalParticleList(1)%w_p / (globalWorld%grid(NumberXNodes) - globalWorld%grid(1)), "should be", n_ave
    del_t = fractionFreq/getPlasmaFreq(n_ave)   
    print *, "Time step (sec) is:", del_t
    print *, "----------------"
    print *, ""
    print *, ""

    ! ! Generate solver object, and then solve for initial rho/potential
    

    call globalSolver%solveInitialPotential(globalParticleList, globalWorld)
    ! call system_clock(count_rate = clock_rate)
    ! call system_clock(time1)
    ! call cpu_time(cpu1)
    ! do i=1, 100
    !     call solvePotential(globalSolver, globalParticleList, globalWorld, del_t, maxIter, eps_r)
    ! end do
    ! call cpu_time(cpu2)
    ! call system_clock(time2)
    ! print *, "Total time is:", real(time2-time1, kind = real64) / real(clock_rate, kind = real64)
    ! print *, "Total CPU time is:", cpu2-cpu1
    ! print *, "Number of time splits is:", amountTimeSplits
    ! stop
    call solveSimulation(globalSolver, globalParticleList, globalWorld, del_t, maxIter, eps_r, irand, simulationTime, heatSkipSteps)
    print *, "Averaging over", averagingTime, "seconds"
    call solveSimulationFinalAverage(globalSolver, globalParticleList, globalWorld, del_t, maxIter, eps_r, irand, averagingTime, heatSkipSteps)
    

    
end program BoundPlasmaExample