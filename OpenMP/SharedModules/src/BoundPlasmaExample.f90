program BoundPlasmaExample
    use iso_fortran_env, only: int32, real64, output_unit, int64
    use constants
    use mod_BasicFunctions
    use mod_particle
    use mod_domain
    use mod_potentialSolver
    use mod_readInputs
    use mod_Scheme
    ! use mod_particleMover
    ! use mod_collisions
    ! use mod_nonLinSolvers
    ! use mod_simulation
    implicit none
    type(Domain) :: globalWorld
    type(Particle), allocatable :: globalParticleList(:)
    type(potentialSolver) :: globalSolver
    integer(int32) :: i, j
    real(real64) :: remainDel_t, currDel_t, E_i, E_f
    call initializeScheme(schemeNum)
    call readInitialInputs('InitialConditions.inp', simulationTime, n_ave, T_e, T_i, numDiagnosticSteps, fractionFreq, averagingTime, numThread, irand)
    call readGeometry(globalWorld, globalSolver, 'Geometry.inp')
    globalParticleList = readParticleInputs('BoundExample.inp', numberChargedParticles, irand, T_e, T_i, numThread, globalWorld)
    !call initializeSolver(eps_r, solverType, m_Anderson, Beta_k, maxIter)
    call depositRho(globalSolver%rho, globalParticleList, globalWorld)
    call globalSolver%solve_tridiag_Poisson(globalWorld)

    
    ! call solveSimulation(globalSolver, globalParticleList, globalWorld, del_t, maxIter, eps_r, irand, simulationTime)
    ! call solveSimulationFinalAverage(globalSolver, globalParticleList, globalWorld, del_t, maxIter, eps_r, irand, averagingTime, 100)
end program BoundPlasmaExample