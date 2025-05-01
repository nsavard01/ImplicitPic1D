program BoundPlasmaExample
    use iso_fortran_env, only: int32, real64, output_unit
    use constants
    use omp_lib 
    use mod_BasicFunctions
    use mod_domain
    use mod_particle
    use mod_particleInjection
    use mod_potentialSolver
    use mod_particle_operations
    use mod_NullCollision
    use mod_simulation
    implicit none

    integer(int32) :: i,j, iThread
    type(Domain) :: world
    type(Particle), allocatable :: particleList(:)
    type(targetParticle), allocatable :: targetParticleList(:)
    type(nullCollision), allocatable :: nullCollisionList(:)
    type(potentialSolver) :: solver
    real(real64), allocatable :: temp(:)
    integer(int32), allocatable :: numTemp(:)   
    
    ! Read initial inputs
    call readInitialInputs('InitialConditions.inp')
    call initializeRandomGenerators(numThread, stateRan0, state_PCG, .false.)
    call readWorld('Geometry.inp', world, T_e, n_ave)
    call readSolver('Geometry.inp', solver, world)
    call readChargedParticleInputs('ParticleTypes.inp', state_PCG, T_e, T_i, numThread, world, particleList)
    call readNeutralParticleInputs('ParticleTypes.inp', targetParticleList)
    call readNullCollisionInputs('collision.inp', nullCollisionList, particleList, targetParticleList)
    call readInjectionInputs('ParticleInjection.inp', particleList(1)%w_p)
    ! Solve initial potential
    call depositRho(solver, particleList, world) 
    call solver%solve_tridiag_Poisson(world, startSimulationTime)

    ! Make field, and if new simulation, rewind particle velocities a half step
    call solver%makeEField(world)
    if (.not. restartBool) call initialVRewind(solver, particleList, del_t, ionStepMult)
    
    ! Solve Simulation
   
    call solveSimulation(solver, particleList, targetParticleList, nullCollisionList, world, del_t, state_PCG, simulationTime)
    if (averagingTime /= 0.0d0) then
        print *, "Averaging over", averagingTime, "seconds"
        call solveSimulationFinalAverage(solver, particleList, targetParticleList, nullCollisionList, world, del_t, state_PCG, averagingTime, 100)
    end if


    

    

    
end program BoundPlasmaExample