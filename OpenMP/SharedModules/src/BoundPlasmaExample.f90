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
    real(real64) :: remainDel_t, currDel_t, E_i, E_f
    call initializeScheme(schemeNum)
    call readInitialInputs('InitialConditions.inp', simulationTime, n_ave, T_e, T_i, numDiagnosticSteps, fractionFreq, averagingTime, numThread, irand)
    call readGeometry(globalWorld, globalSolver, 'Geometry.inp')
    globalParticleList = readParticleInputs('BoundExample.inp', numberChargedParticles, irand, T_e, T_i, numThread, globalWorld)
    do i = 1, numberChargedParticles
        globalParticleList(i)%N_p = 0
    end do
    call readInjectionInputs('ParticleInjection.inp', addLostPartBool, refluxPartBool, injectionBool, injectionFlux, globalParticleList(1)%w_p)
    call initializeSolver(eps_r, solverType, m_Anderson, Beta_k, maxIter)
    if (injectionBool) call injectAtBoundary(globalParticleList, T_e, T_i, irand, globalWorld, del_t)
    call solveInitialPotential(globalSolver, globalParticleList, globalWorld)
    ! E_i = globalSolver%getTotalPE(globalWorld, .false.)
    ! do j=1, numberChargedParticles
    !     E_i = E_i + globalParticleList(j)%getTotalKE()
    ! end do
    ! remainDel_t = del_t
    ! currDel_t = del_t
    ! call solvePotential(globalSolver, globalParticleList, globalWorld, del_t, remainDel_t, currDel_t, maxIter, eps_r)
    ! E_f = globalSolver%getTotalPE(globalWorld, .false.)
    ! do j=1, numberChargedParticles
    !     E_f = E_f + globalParticleList(j)%getTotalKE() + SUM(globalParticleList(j)%energyLoss) * globalParticleList(j)%mass * globalParticleList(j)%w_p * 0.5d0
    ! end do
    ! print *, ABS((E_i - E_f)/(E_i))
    ! print *, 'took', iterNumPicard, 'iterations'
    ! call depositRho(globalSolver%rho, globalParticleList, globalWorld)
    ! print *, 'gauss error is:', globalSolver%getError_tridiag_Poisson(globalWorld)
    ! if (injectionBool) call injectAtBoundary(globalParticleList, T_e, T_i, irand, globalWorld, del_t)
    ! do j = 1, numberChargedParticles
    !     print *, "--------------"
    !     print *, 'particle is:', globalParticleList(j)%name
    !     print *, "--------------"
    !     do iThread = 1, numThread
    !         if (globalParticleList(j)%N_p(iThread) > 0) then
    !             print *, 'l is:', globalParticleList(j)%phaseSpace(1, globalParticleList(j)%N_p(iThread), iThread)
    !             print *, 'v_x is:', globalParticleList(j)%phaseSpace(2, globalParticleList(j)%N_p(iThread), iThread)
    !         end if
    !     end do
    ! end do
    ! stop
    call solveSimulation(globalSolver, globalParticleList, globalWorld, del_t, maxIter, eps_r, irand, simulationTime)
    call solveSimulationFinalAverage(globalSolver, globalParticleList, globalWorld, del_t, maxIter, eps_r, irand, averagingTime, 100)
end program BoundPlasmaExample