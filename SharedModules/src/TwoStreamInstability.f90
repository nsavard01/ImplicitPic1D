program TwoStreamInstability
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
    globalParticleList = readParticleInputs('ParticlesTwoStream.inp',numberChargedParticles, irand, T_e) 
    ! Initialize constants with inputs
    ! create the world the particles live in
    call readInputs(NumberXNodes, numDiagnosticSteps, averagingTime, fractionFreq, n_ave, globalWorld, globalSolver, simulationTime, Power, heatSkipSteps, nu_h, T_e, 'TwoStreamGeometry.inp', 'TwoStreamInitialConditions.inp')
    do i = 1, numberChargedParticles
        call initialize_randUniform(globalParticleList(i), globalWorld, irand)
        call globalParticleList(i) % initialize_n_ave(n_ave, globalWorld%grid(NumberXNodes) - globalWorld%grid(1))
    end do
    do i = 1, globalParticleList(1)%N_p
        globalParticleList(1)%phaseSpace(1, i + globalParticleList(1)%N_p) = globalParticleList(1)%phaseSpace(1, i)
        globalParticleList(1)%phaseSpace(2, i) = -2000.0d0
        globalParticleList(1)%phaseSpace(2, i + globalParticleList(1)%N_p) = 2000.0d0
    end do
    globalParticleList(1)%N_p = 2 * globalParticleList(1)%N_p
    call initializeSolver(eps_r, solverType, m_Anderson, Beta_k, maxIter)
    print *, "Calulated values:"
    print *, "Number of particles is:", globalParticleList(1)%N_p
    print *, "w_p is:", globalParticleList(1)%w_p
    print *, "Debye length is:", getDebyeLength(T_e, 2.0d0 * n_ave)
    print *, "Plasma frequency is:", getPlasmaFreq(2.0d0 * n_ave)
    print *, "Average density is ", globalParticleList(1)%N_p * globalParticleList(1)%w_p / (globalWorld%grid(NumberXNodes) - globalWorld%grid(1)), "should be", 2.0d0 * n_ave
    del_t = fractionFreq/getPlasmaFreq(2.0d0 * n_ave)   
    print *, "Time step (sec) is:", del_t
    print *, "----------------"
    ! Generate solver object, and then solve for initial rho/potential
    globalSolver%rho_const = 2.0d0 * n_ave * e
    call solveInitialPotential(globalSolver, globalParticleList, globalWorld)
    ! call solveSingleTimeStepDiagnostic(globalSolver, globalParticleList, globalWorld, del_t, maxIter, eps_r)
    ! print *, "Energy error is:", globalSolver%energyError
    ! ! Get error gauss' law
    ! call depositRho(globalSolver%rho, globalParticleList, globalWorld)
    ! call globalSolver%construct_diagMatrix(globalWorld)
    ! globalSolver%chargeError = globalSolver%getError_tridiag_Poisson(globalWorld)
    ! globalSolver%chargeError = globalSolver%chargeError / SQRT(SUM(globalSolver%rho**2))
    ! call globalSolver%construct_diagMatrix_Ampere(globalWorld)
    ! print *, "Charge error is:", globalSolver%chargeError
    ! stop 
    if (Power /= 0.0d0) then
        call solveSimulation(globalSolver, globalParticleList, globalWorld, del_t, maxIter, eps_r, irand, simulationTime, heatSkipSteps)
    else
        call solveSimulationOnlyPotential(globalSolver, globalParticleList, globalWorld, del_t, maxIter, eps_r, simulationTime)
    end if
    if (averagingTime /= 0.0d0) then
        print *, "Averaging over", averagingTime, "seconds"
        call solveSimulationFinalAverage(globalSolver, globalParticleList, globalWorld, del_t, maxIter, eps_r, irand, averagingTime, heatSkipSteps)
    end if
            


end program TwoStreamInstability