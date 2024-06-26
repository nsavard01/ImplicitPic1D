program BoundPlasmaExample
    use iso_fortran_env, only: int32, real64, output_unit, int64
    use constants
    use mod_BasicFunctions
    use mod_particle
    use mod_targetParticle
    use mod_NullCollision
    use mod_domain
    use mod_potentialSolver
    use mod_Scheme
    use mod_particleInjection
    use mod_particleMover
    ! use mod_collisions
    use mod_nonLinSolvers
    use mod_simulation
    use omp_lib
    implicit none
    
    integer(int32) :: i, j, iThread, startTime, endTime
    real(real64) :: remainDel_t, currDel_t, PE_i, KE_i, PE_f, KE_f, Momentum_i(3), Momentum_f(3)
    ! real(real64), allocatable :: rho_i(:), J_total(:)
    type(targetParticle), allocatable :: targetParticleList(:)
    type(nullCollision), allocatable :: nullCollisionList(:)

    ! Read all inputs and create objects
    call initializeScheme()
    call readInitialInputs('InitialConditions.inp')
    call initializeRandomGenerators(numThread, stateRan0, stateRanNew, .false.)
    call readWorld('Geometry.inp', globalWorld, T_e, n_ave)
    call readSolver('Geometry.inp', globalSolver, globalWorld, startSimulationTime)
    call readChargedParticleInputs('ParticleTypes.inp', stateRan0, T_e, T_i, numThread, globalWorld, globalParticleList)
    call readNeutralParticleInputs('ParticleTypes.inp', targetParticleList)
    call readNullCollisionInputs('collision.inp', nullCollisionList, globalParticleList, targetParticleList)
    call allocateParticleMoverData()
    call readInjectionInputs('ParticleInjection.inp', globalParticleList(1)%w_p)
    call initializeSolver()
    ! Solve initial potential
    call globalSolver%solveInitialPotential(globalParticleList, globalWorld, startSimulationTime)
    ! Used for testing, comment out
    ! allocate(rho_i(NumberXNodes), J_total(NumberXHalfNodes))
    ! rho_i = globalSolver%rho
    ! KE_i = 0.0d0
    ! do j=1, numberChargedParticles
    !     KE_i = KE_i + globalParticleList(j)%getTotalKE()
    ! end do
    ! remainDel_t = del_t
    ! currDel_t = del_t
    ! Momentum_i = 0
    ! do j = 1, numberChargedParticles
    !     Momentum_i = Momentum_i + globalParticleList(j)%getTotalMomentum()
    ! end do
    
    ! call system_clock(startTime)
    ! call solvePotential(globalSolver, globalParticleList, globalWorld, del_t, remainDel_t, currDel_t, startSimulationTime)
    ! call system_clock(endTime)
   
    ! PE_i = globalSolver%getTotalPE(globalWorld, .false.)
    ! KE_f = 0.0d0
    ! do i = 1, numberChargedParticles
    !     print *, 'Particle:', globalParticleList(i)%name
    !     print *, 'Remaining number of particles:', SUM(globalParticleList(i)%N_p)
    !     print *, 'Ave. num substeps:', globalParticleList(i)%numSubStepsAve
    !     print *, 'Ave. num func evals:', globalParticleList(i)%numFuncEvalAve
    !     print *, 'Kept KE:', globalParticleList(i)%getTotalKE()
    !     print *, 'Lost KE', SUM(globalParticleList(i)%energyLoss) * globalParticleList(i)%mass * globalParticleList(i)%w_p * 0.5d0
    !     KE_f = KE_f + globalParticleList(i)%getTotalKE() + SUM(globalParticleList(i)%energyLoss) * globalParticleList(i)%mass * globalParticleList(i)%w_p * 0.5d0
    ! end do
    ! PE_f = globalSolver%getTotalPE(globalWorld, .true.) - globalSolver%getEnergyFromBoundary(globalWorld, currDel_t) 
    ! print *, 'PE_f:', PE_f
    ! Momentum_f = 0
    ! do j = 1, numberChargedParticles
    !     Momentum_f = Momentum_f + globalParticleList(j)%getTotalMomentum()
    ! end do
    ! print *, 'Time step:', currDel_t
    ! print *, 'Energy error is:', ABS((PE_i + KE_i - PE_f - KE_f)/(PE_i + KE_i))
    ! print *, 'Total KE loss:', KE_f - KE_i
    ! print *, 'Total PE loss:', PE_f - PE_i
    ! print *, 'JE:', globalSolver%getJdotE(globalWorld, currDel_t)
    ! print *, 'Energy in:', globalSolver%getEnergyFromBoundary(globalWorld, currDel_t)
    ! print *, 'took', iterNumPicard, 'iterations'
    ! print *, 'took amount of integer time:', endTime - startTime
    ! print *, 'Momentum percent diff:', 100 * ABS((Momentum_i(1) - Momentum_f(1))/Momentum_i(1))
    ! call globalSolver%depositRho(globalParticleList, globalWorld)
    ! print *, 'gauss error first is:', globalSolver%getError_tridiag_Poisson(globalWorld)
    ! do j = 1, numberBinaryCollisions
    !     KE_i = globalParticleList(j)%getTotalKE()
    !     call nullCollisionList(j)%generateCollision(globalParticleList, targetParticleList, numberChargedParticles, numberBinaryCollisions, stateRan0, del_t)
    !     KE_f = globalParticleList(j)%getTotalKE()
    !     print *, 'energy loss for collision', j
    !     print *, nullCollisionList(j)%totalEnergyLoss * e * globalParticleList(j)%w_p
    !     print *, 'energy loss from particles:', KE_f - KE_i
    ! end do
    ! call globalSolver%depositRho(globalParticleList, globalWorld)
    ! print *, 'gauss error is:', globalSolver%getError_tridiag_Poisson(globalWorld)
    ! print *, 'charge error is:', globalSolver%getChargeContinuityError(rho_i, globalWorld, currDel_t)
    ! stop

    ! Run simulations
    call solveSimulation(globalSolver, globalParticleList, targetParticleList, nullCollisionList, globalWorld, del_t, stateRan0, simulationTime)
    call solveSimulationFinalAverage(globalSolver, globalParticleList, targetParticleList, nullCollisionList, globalWorld, del_t, stateRan0, averagingTime, 100)
end program BoundPlasmaExample