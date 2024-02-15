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
    
    integer(int32) :: i, j, iThread
    real(real64) :: remainDel_t, currDel_t, PE_i, KE_i, PE_f, KE_f
    real(real64), allocatable :: rho_i(:)
    character(:), allocatable :: saveFolderName !name of the particle
    type(targetParticle), allocatable :: targetParticleList(:)
    type(nullCollision), allocatable :: nullCollisionList(:)

    saveFolderName = '../../../../ImplicitData-BField/'
    call initializeScheme()
    call readInitialInputs('InitialConditions.inp', saveFolderName)
    call initializeRandomGenerators(numThread, stateRan0, stateRanNew)
    call readWorld('Geometry.inp', globalWorld, T_e, n_ave)
    call readSolver('Geometry.inp', globalSolver, globalWorld)
    call readChargedParticleInputs('BoundExample.inp', stateRan0, T_e, T_i, numThread, globalWorld, globalParticleList)
    call readNeutralParticleInputs('BoundExample.inp', targetParticleList)
    call readNullCollisionInputs('collision.inp', nullCollisionList, globalParticleList, targetParticleList)
    
    ! do i = 1, numberChargedParticles
    !     globalParticleList(i)%N_p = 0
    ! end do
    call readInjectionInputs('ParticleInjection.inp', globalParticleList(1)%w_p, globalSolver%BFieldAngle)
    call initializeSolver()
    ! if (injectionBool) call injectAtBoundary(globalParticleList, T_e, T_i, irand, globalWorld, del_t, globalSolver%BFieldAngle)
    call solveInitialPotential(globalSolver, globalParticleList, globalWorld, 0.0d0)
    ! print *, 'electron particle number:', SUM(globalParticleList(1)%N_p)
    ! print *, 'ion particle number:', SUM(globalParticleList(2)%N_p)
    ! E_i = globalSolver%getTotalPE(globalWorld, .false.)
    ! do j=1, numberChargedParticles
    !     E_i = E_i + globalParticleList(j)%getTotalKE()
    ! end do
    ! call moveParticles(globalSolver, globalParticleList, globalWorld, del_t)
    ! print *, 'electron particle number:', SUM(globalParticleList(1)%N_p)
    ! print *, 'ion particle number:', SUM(globalParticleList(2)%N_p)
    ! E_f = globalSolver%getTotalPE(globalWorld, .false.)
    ! do j=1, numberChargedParticles
    !     E_f = E_f + globalParticleList(j)%getTotalKE() + SUM(globalParticleList(j)%energyLoss) * globalParticleList(j)%mass * globalParticleList(j)%w_p * 0.5d0
    ! end do
    ! print *, ABS((E_i - E_f)/(E_i))
    ! stop
  
    ! allocate(rho_i(NumberXNodes))
    ! rho_i = globalSolver%rho
    ! KE_i = 0.0d0
    ! do j=1, numberChargedParticles
    !     KE_i = KE_i + globalParticleList(j)%getTotalKE()
    ! end do
    ! remainDel_t = del_t
    ! currDel_t = del_t
    ! call solvePotential(globalSolver, globalParticleList, globalWorld, del_t, remainDel_t, currDel_t, maxIter, eps_r, 0.0d0)
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
    
    ! print *, 'Time step:', currDel_t
    ! print *, 'Energy error is:', ABS((PE_i + KE_i - PE_f - KE_f)/(PE_i + KE_i))
    ! print *, 'took', iterNumPicard, 'iterations'

    ! do j = 1, numberBinaryCollisions
    !     KE_i = globalParticleList(j)%getTotalKE()
    !     call nullCollisionList(j)%generateCollision(globalParticleList, targetParticleList, numberChargedParticles, numberBinaryCollisions, stateRan0, del_t)
    !     KE_f = globalParticleList(j)%getTotalKE()
    !     print *, 'energy loss for collision', j
    !     print *, nullCollisionList(j)%totalEnergyLoss * e * globalParticleList(j)%w_p
    !     print *, 'energy loss from particles:', KE_f - KE_i
    ! end do
    ! call depositRho(globalSolver%rho, globalParticleList, globalWorld)
    ! print *, 'gauss error is:', globalSolver%getError_tridiag_Poisson(globalWorld)
    ! print *, 'charge error is:', getChargeContinuityError(rho_i, globalSolver%rho, globalSolver%J, globalWorld, currDel_t)
    ! stop
    call solveSimulation(globalSolver, globalParticleList, targetParticleList, nullCollisionList, globalWorld, del_t, maxIter, eps_r, stateRan0, simulationTime)
    call solveSimulationFinalAverage(globalSolver, globalParticleList, targetParticleList, nullCollisionList, globalWorld, del_t, maxIter, eps_r, stateRan0, averagingTime, 100)
end program BoundPlasmaExample