program BoundPlasmaExample
    use iso_fortran_env, only: int32, real64, output_unit, int64
    use constants
    use mod_BasicFunctions
    use mod_particle
    use mod_domain
    use mod_Scheme
    use mod_potentialSolver
    use mod_readInputs
    ! use mod_particleMover
    ! use mod_collisions
    ! use mod_nonLinSolvers
    ! use mod_simulation
    implicit none

    integer(int32) :: i, j
    type(Domain) :: globalWorld
    type(Particle), allocatable :: globalParticleList(:)
    type(potentialSolver) :: globalSolver
    call initializeScheme(schemeNum)
    call readInitialInputs('InitialConditions.inp', simulationTime, n_ave, T_e, T_i, numDiagnosticSteps, fractionFreq, averagingTime, numThread, irand)
    call readGeometry(globalWorld, globalSolver, 'Geometry.inp')
    globalParticleList = readParticleInputs('BoundExample.inp', numberChargedParticles, irand, T_e, T_i, numThread, globalWorld)
    call depositRho(globalSolver%rho, globalParticleList, globalWorld)
    call globalSolver%solve_tridiag_Poisson(globalWorld)
    print *, globalSolver%phi
    
    

 
end program BoundPlasmaExample