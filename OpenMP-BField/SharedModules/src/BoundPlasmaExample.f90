program BoundPlasmaExample
    use iso_fortran_env, only: int32, real64, output_unit, int64
    use constants
    use mod_BasicFunctions
    use mod_particle
    use mod_domain
    use mod_potentialSolver
    use mod_readInputs
    use mod_Scheme
    use mod_particleInjection
    use mod_particleMover
    ! use mod_collisions
    use mod_nonLinSolvers
    use mod_simulation
    use omp_lib
    implicit none
    
    integer(int32) :: i, j, iThread
    real(real64) :: remainDel_t, currDel_t, E_i, E_f, EJ
    call initializeScheme(schemeNum)
    call readInitialInputs('InitialConditions.inp', simulationTime, n_ave, T_e, T_i, numDiagnosticSteps, fractionFreq, averagingTime, numThread, irand)
    call readGeometry(globalWorld, globalSolver, 'Geometry.inp')
    globalParticleList = readParticleInputs('BoundExample.inp', numberChargedParticles, irand, T_e, T_i, numThread, globalWorld)
    ! do i = 1, numberChargedParticles
    !     globalParticleList(i)%N_p = 0
    ! end do
    call readInjectionInputs('ParticleInjection.inp', addLostPartBool, refluxPartBool, injectionBool, injectionFlux, globalParticleList(1)%w_p)
    call initializeSolver(eps_r, solverType, m_Anderson, Beta_k, maxIter)
    ! if (injectionBool) call injectAtBoundary(globalParticleList, T_e, T_i, irand, globalWorld, del_t)
    call solveInitialPotential(globalSolver, globalParticleList, globalWorld)
    ! print *, 'electron particle number:', SUM(globalParticleList(1)%N_p)
    ! print *, 'ion particle number:', SUM(globalParticleList(2)%N_p)
    ! E_i = globalSolver%getTotalPE(globalWorld, .false.)
    ! do j=1, numberChargedParticles
    !     E_i = E_i + globalParticleList(j)%getTotalKE()
    ! end do
    ! call depositJ(globalSolver, globalParticleList, globalWorld, del_t)
    ! print *, 'electron particle number:', SUM(globalParticleList(1)%N_p)
    ! print *, 'ion particle number:', SUM(globalParticleList(2)%N_p)
    ! E_f = globalSolver%getTotalPE(globalWorld, .false.)
    ! do j=1, numberChargedParticles
    !     E_f = E_f + globalParticleList(j)%getTotalKE() + SUM(globalParticleList(j)%energyLoss) * globalParticleList(j)%mass * globalParticleList(j)%w_p * 0.5d0
    ! end do
    ! print *, ABS((E_i - E_f)/(E_i))
    ! stop
    E_i = 0.0d0 !globalSolver%getTotalPE(globalWorld, .false.)
    do j=1, numberChargedParticles
        E_i = E_i + globalParticleList(j)%getTotalKE()
    end do
    remainDel_t = del_t
    currDel_t = del_t
    call solvePotential(globalSolver, globalParticleList, globalWorld, del_t, remainDel_t, currDel_t, maxIter, eps_r)
    E_f = 0.0d0 !globalSolver%getTotalPE(globalWorld, .false.)
    do j=1, numberChargedParticles
        E_f = E_f + globalParticleList(j)%getTotalKE() + SUM(globalParticleList(j)%energyLoss) * globalParticleList(j)%mass * globalParticleList(j)%w_p * 0.5d0
    end do
    EJ = SUM(globalSolver%EField * SUM(globalSolver%J, DIM = 2))
    print *, EJ, (E_f - E_i) * globalParticleList(1)%w_p
    print *, ABS((E_i - E_f)/(E_i))
    print *, 'took', iterNumPicard, 'iterations'
    call depositRho(globalSolver%rho, globalParticleList, globalWorld)
    print *, 'gauss error is:', globalSolver%getError_tridiag_Poisson(globalWorld)
    stop
    call solveSimulation(globalSolver, globalParticleList, globalWorld, del_t, maxIter, eps_r, irand, simulationTime)
    call solveSimulationFinalAverage(globalSolver, globalParticleList, globalWorld, del_t, maxIter, eps_r, irand, averagingTime, 100)
end program BoundPlasmaExample