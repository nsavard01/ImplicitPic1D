program BoundPlasmaExample
    use iso_fortran_env, only: int32, real64, output_unit
    use constants
    use omp_lib 
    use mod_BasicFunctions
    use mod_domain
    use mod_particle
    use mod_particleInjection
    use mod_potentialSolver
    use mod_NullCollision
    use mod_simulation
    implicit none

    integer(int32) :: i,j, iThread!, tclock1, tclock2, clock_rate
    type(Domain) :: world
    type(Particle), allocatable :: particleList(:)
    type(targetParticle), allocatable :: targetParticleList(:)
    type(nullCollision), allocatable :: nullCollisionList(:)
    type(potentialSolver) :: solver
    real(real64), allocatable :: temp(:)
    integer(int32), allocatable :: numTemp(:)   
   
    call readInitialInputs('InitialConditions.inp')
    call initializeRandomGenerators(numThread, stateRan0, stateRanNew, .false.)
    call readWorld('Geometry.inp', world, T_e, n_ave)
    call readSolver('Geometry.inp', solver, world)
    call readChargedParticleInputs('ParticleTypes.inp', stateRan0, T_e, T_i, numThread, world, particleList)
    call readNeutralParticleInputs('ParticleTypes.inp', targetParticleList)
    call readNullCollisionInputs('collision.inp', nullCollisionList, particleList, targetParticleList)
    ! do i = 1, numberChargedParticles
    !     particleList(i)%N_p = 0
    ! end do
    call readInjectionInputs('ParticleInjection.inp', particleList(1)%w_p)
    ! if (injectionBool) call injectAtBoundary(particleList, T_e, T_i, irand, world, del_t)
    call solver%depositRho(particleList, world)
    call solver%solve_tridiag_Poisson(world, startSimulationTime)

    ! Assume only use potential solver once, then need to generate matrix for Div-Ampere
    call solver%makeEField(world)
    if (.not. restartBool) call solver%initialVRewind(particleList, del_t)
    
    call solveSimulation(solver, particleList, targetParticleList, nullCollisionList, world, del_t, stateRan0, simulationTime)
    if (averagingTime /= 0.0d0) then
        print *, "Averaging over", averagingTime, "seconds"
        call solveSimulationFinalAverage(solver, particleList, targetParticleList, nullCollisionList, world, del_t, stateRan0, averagingTime, 100)
    end if


    

    

    
end program BoundPlasmaExample