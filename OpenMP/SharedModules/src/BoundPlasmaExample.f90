program BoundPlasmaExample
    use iso_fortran_env, only: int32, real64, output_unit, int64
    use constants
    use mod_BasicFunctions
    use mod_particle
    use mod_domain
    use mod_Scheme
    use mod_potentialSolver
    use mod_readInputs
    use mod_particleMover
    ! use mod_collisions
    use mod_nonLinSolvers
    ! use mod_simulation
    implicit none

    integer(int32) :: i, j
    real(real64) :: remainDel_t, currDel_t
    call initializeScheme(schemeNum)
    call readInitialInputs('InitialConditions.inp', simulationTime, n_ave, T_e, T_i, numDiagnosticSteps, fractionFreq, averagingTime, numThread, irand)
    call readGeometry(globalWorld, globalSolver, 'Geometry.inp')
    globalParticleList = readParticleInputs('BoundExample.inp', numberChargedParticles, irand, T_e, T_i, numThread, globalWorld)
    call initializeSolver(eps_r, solverType, m_Anderson, Beta_k, maxIter)
    call solveInitialPotential(globalSolver, globalParticleList, globalWorld)
    remainDel_t = del_t
    currDel_t = del_t
    call solvePotential(globalSolver, globalParticleList, globalWorld, del_t, remainDel_t, currDel_t, maxIter, eps_r)
    print *, 'Number of iterations was:', iterNumPicard
end program BoundPlasmaExample