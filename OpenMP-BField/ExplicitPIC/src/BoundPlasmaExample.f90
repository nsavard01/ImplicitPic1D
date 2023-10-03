program BoundPlasmaExample
    use iso_fortran_env, only: int32, real64, output_unit
    use constants
    use omp_lib 
    use mod_BasicFunctions
    use mod_domain
    use mod_particle
    use mod_potentialSolver
    ! use mod_collisions
    use mod_readInputs
    use mod_simulation
    implicit none

    integer(int32) :: i,j!, tclock1, tclock2, clock_rate
    type(Domain) :: world
    type(Particle), allocatable :: particleList(:)
    type(potentialSolver) :: solver
    call readInitialInputs('InitialConditions.inp', simulationTime, n_ave, T_e, T_i, numDiagnosticSteps, fractionFreq, averagingTime, numThread, irand)
    call readGeometry(world, solver, 'Geometry.inp')
    particleList = readParticleInputs('BoundExample.inp', numberChargedParticles, irand, T_e, T_i, numThread, world)
    call solver%depositRho(particleList)
    call solver%solve_tridiag_Poisson(world)
    ! Assume only use potential solver once, then need to generate matrix for Div-Ampere
    call solver%makeEField(world)
    call solver%initialVRewind(particleList, del_t)
    call solveSimulation(solver, particleList, world, del_t, irand, simulationTime)
    if (averagingTime /= 0.0d0) then
        print *, "Averaging over", averagingTime, "seconds"
        call solveSimulationFinalAverage(solver, particleList, world, del_t, irand, averagingTime, 100)
    end if


    

    

    
end program BoundPlasmaExample