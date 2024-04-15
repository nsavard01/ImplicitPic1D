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
    real(real64) :: l_temp, v_temp
    
   
    call readInitialInputs('InitialConditions.inp')
    call initializeRandomGenerators(numThread, stateRan0, stateRanNew, .false.)
    call readWorld('Geometry.inp', world, T_e, n_ave)
    call readSolver('Geometry.inp', solver, world)
    call readChargedParticleInputs('BoundExample.inp', stateRan0, T_e, T_i, numThread, world, particleList)
    call readNeutralParticleInputs('BoundExample.inp', targetParticleList)
    call readNullCollisionInputs('collision.inp', nullCollisionList, particleList, targetParticleList)
    ! do i = 1, numberChargedParticles
    !     particleList(i)%N_p = 0
    ! end do
    call readInjectionInputs('ParticleInjection.inp', particleList(1)%w_p, solver%BFieldAngle)
    ! if (injectionBool) call injectAtBoundary(particleList, T_e, T_i, irand, world, del_t)
    call solver%depositRho(particleList, world)
    call solver%solve_tridiag_Poisson(world, 0.0d0)
    ! print *, 'mass ratio:', particleList(2)%mass / particleList(1)%mass
    ! print *, SUM(particleList(2)%workSpace, DIM = 2)
    ! print *, ""
    ! allocate(temp(NumberXNodes-1), numTemp(NumberXNodes-1))
    ! temp = 0.0
    ! numTemp = 0
    ! do i = 1, numThread
    !     do j = 1, particleList(2)%N_p(i)
    !         l_temp = particleList(2)%phaseSpace(1, j, i)
    !         v_temp = particleList(2)%phaseSpace(2, j, i)
    !         temp(INT(l_temp)) = temp(INT(l_temp)) + v_temp
    !         numTemp(INT(l_temp)) = numTemp(INT(l_temp)) + 1
    !     end do
    ! end do

    ! print *, temp/numTemp
    ! stop
    ! do j = 1, numberChargedParticles
    !     print *, "--------------"
    !     print *, 'particle is:', particleList(j)%name
    !     print *, "--------------"
    !     do iThread = 1, numThread
    !         if (particleList(j)%N_p(iThread) > 0) then
    !             print *, 'l is:', particleList(j)%phaseSpace(1, particleList(j)%N_p(iThread), iThread)
    !             print *, 'v_x is:', particleList(j)%phaseSpace(2, particleList(j)%N_p(iThread), iThread)
    !         end if
    !     end do
    ! end do
    ! stop
    ! Assume only use potential solver once, then need to generate matrix for Div-Ampere
    call solver%makeEField(world)
    call solver%initialVRewind(particleList, del_t)
    
    call solveSimulation(solver, particleList, targetParticleList, nullCollisionList, world, del_t, stateRan0, simulationTime)
    if (averagingTime /= 0.0d0) then
        print *, "Averaging over", averagingTime, "seconds"
        call solveSimulationFinalAverage(solver, particleList, targetParticleList, nullCollisionList, world, del_t, stateRan0, averagingTime, 100)
    end if


    

    

    
end program BoundPlasmaExample