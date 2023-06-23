program BoundPlasmaExample
    use iso_fortran_env, only: int32, real64, output_unit, int64
    use constants
    use mod_BasicFunctions
    use mod_domain
    use mod_particle
    use mod_potentialSolver
    use mod_particleMover
    use mod_collisions
    use mod_nonLinSolvers
    use mod_Scheme
    use mod_simulation
    implicit none

    integer(int32) :: i
    call initializeScheme(schemeNum)
    ! Initialize constants with inputs
    ! create the world the particles live in
    call readInputs(NumberXNodes, numDiagnosticSteps, averagingTime, fractionFreq, n_ave, globalWorld, globalSolver, simulationTime, Power, heatSkipSteps, nu_h, T_e, T_i, 'Geometry.inp', 'InitialConditions.inp')
    globalParticleList = readParticleInputs('BoundExample.inp',numberChargedParticles, irand, T_e, T_i) 
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
    ! Generate solver object, and then solve for initial rho/potential
    call solveInitialPotential(globalSolver, globalParticleList, globalWorld)
    ! call solveSingleTimeStepDiagnostic(globalSolver, globalParticleList, globalWorld, del_t, maxIter, eps_r)
    ! print *, "Energy error is:", globalSolver%energyError
    ! print *, "Took", iterNumPicard, "iterations"
    ! print *, 'electron number:', globalParticleList(1)%N_P
    ! print *, 'ion number:', globalParticleList(2)%N_p
    ! call addMaxwellianLostParticles(globalParticleList, T_e, 0.1d0, irand, delIdx, idxReFlux, reFluxMaxIdx, 0.03d0, globalWorld)
    ! print *, 'electron number:', globalParticleList(1)%N_P
    ! print *, 'ion number:', globalParticleList(2)%N_p
    ! ! Get error gauss' law
    ! call depositRho(globalSolver%rho, globalParticleList, globalWorld)
    ! call globalSolver%construct_diagMatrix(globalWorld)
    ! globalSolver%chargeError = globalSolver%getError_tridiag_Poisson(globalWorld)
    ! globalSolver%chargeError = globalSolver%chargeError / SQRT(SUM(globalSolver%rho**2))
    ! call globalSolver%construct_diagMatrix_Ampere(globalWorld)
    ! print *, "Charge error is:", globalSolver%chargeError
    ! stop 
    call solveSimulation(globalSolver, globalParticleList, globalWorld, del_t, maxIter, eps_r, irand, simulationTime, heatSkipSteps)

    print *, "Averaging up to", averagingTime, "simulation seconds"
    call solveSimulationFinalAverage(globalSolver, globalParticleList, globalWorld, del_t, maxIter, eps_r, irand, averagingTime, heatSkipSteps)

    

    
end program BoundPlasmaExample