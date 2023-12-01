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
    ! ! use mod_collisions
    use mod_nonLinSolvers
    ! use mod_simulation
    use omp_lib
    implicit none
    
    integer(int32) :: i, j, iThread
    real(real64) :: remainDel_t, currDel_t, KE_i, KE_f, PE_i, PE_f, EJ, delJ_phi
    real(real64), allocatable :: rho_i(:), Field_past(:)
    call initializeScheme(schemeNum)
    call readInitialInputs('InitialConditions.inp', simulationTime, n_ave, T_e, T_i, numDiagnosticSteps, fractionFreq, averagingTime, numThread, irand)
    call readGeometry(globalWorld, globalSolver, 'Geometry.inp')
    globalParticleList = readParticleInputs('BoundExample.inp', numberChargedParticles, irand, T_e, T_i, numThread, globalWorld)
    call initializeSolver(eps_r, solverType, m_Anderson, Beta_k, maxIter)
    call depositRho(globalSolver%rho, globalParticleList, globalWorld)
    call globalSolver%construct_diagMatrix(globalWorld)
    call globalSolver%solve_tridiag_Poisson(globalWorld)
    print *, 'past phi:'
    print *, globalSolver%phi
    print *, ''
    call globalSolver%evaluateEFieldCurrTime(globalWorld)
    print *, globalSolver%EField
    print *, ''
    allocate(rho_i(NumberXNodes))
    rho_i = globalSolver%rho
    print *, 'past rho is:'
    print *, rho_i
    print *, ''
    
    ! call depositJ(globalSolver, globalParticleList, globalWorld, del_t)
    ! call moveParticles(globalSolver, globalParticleList, globalWorld, del_t)
    ! print *, 'J is:'
    ! print *, SUM(globalSolver%J, dim = 2)
    ! call depositRho(globalSolver%rho, globalParticleList, globalWorld)
    ! print *, 'Future rho is:'
    ! print *, globalSolver%rho
    ! print *, 'Number electrons left:', SUM(globalParticleList(1)%N_p)
    ! print *, 'Number ions left:', SUM(globalParticleList(2)%N_p)
    ! print *, ''
    ! do i = 1, NumberXNodes
    !     if (i == 1) then
    !         print *, 1.0d0 + del_t * (SUM(globalSolver%J(i+1, :)) - 2.0d0 * SUM(globalSolver%J(i,:)))/(globalSolver%rho(i) - rho_i(i))
    !     else if (i == NumberXNodes) then
    !         print *, 1.0d0 + del_t * ( - SUM(globalSolver%J(i,:)))/(globalSolver%rho(i) - rho_i(i))
    !     else
    !         print *, 1.0d0 + del_t * (SUM(globalSolver%J(i+1, :)) - SUM(globalSolver%J(i,:)))/(globalSolver%rho(i) - rho_i(i))
    !     end if
    ! end do

    ! do i = 1, numberChargedParticles
    !     globalParticleList(i)%N_p = 0
    ! end do
    ! call readInjectionInputs('ParticleInjection.inp', addLostPartBool, refluxPartBool, injectionBool, injectionFlux, globalParticleList(1)%w_p, globalSolver%BFieldAngle)
    ! call initializeSolver(eps_r, solverType, m_Anderson, Beta_k, maxIter)
    ! ! if (injectionBool) call injectAtBoundary(globalParticleList, T_e, T_i, irand, globalWorld, del_t, globalSolver%BFieldAngle)
    ! call solveInitialPotential(globalSolver, globalParticleList, globalWorld)
    
    ! ! print *, 'electron particle number:', SUM(globalParticleList(1)%N_p)
    ! ! print *, 'ion particle number:', SUM(globalParticleList(2)%N_p)
    ! ! E_i = globalSolver%getTotalPE(globalWorld, .false.)
    ! ! do j=1, numberChargedParticles
    ! !     E_i = E_i + globalParticleList(j)%getTotalKE()
    ! ! end do
    ! ! call moveParticles(globalSolver, globalParticleList, globalWorld, del_t)
    ! ! print *, 'electron particle number:', SUM(globalParticleList(1)%N_p)
    ! ! print *, 'ion particle number:', SUM(globalParticleList(2)%N_p)
    ! ! E_f = globalSolver%getTotalPE(globalWorld, .false.)
    ! ! do j=1, numberChargedParticles
    ! !     E_f = E_f + globalParticleList(j)%getTotalKE() + SUM(globalParticleList(j)%energyLoss) * globalParticleList(j)%mass * globalParticleList(j)%w_p * 0.5d0
    ! ! end do
    ! ! print *, ABS((E_i - E_f)/(E_i))
    ! ! stop
    allocate(Field_past(NumberXNodes+1))
    PE_i = globalSolver%getTotalPE(globalWorld, .false.)
    KE_i = 0.0d0
    do j=1, numberChargedParticles
        KE_i = KE_i + globalParticleList(j)%getTotalKE()
    end do
    remainDel_t = del_t
    currDel_t = del_t
    call globalSolver%evaluateEFieldCurrTime(globalWorld)
    Field_past = globalSolver%EField
    call solvePotential(globalSolver, globalParticleList, globalWorld, del_t, remainDel_t, currDel_t, maxIter, eps_r)
    print *, ''
    print *, 'new phi:'
    print *, globalSolver%phi
    PE_f = globalSolver%getTotalPE(globalWorld, .false.)
    KE_f = 0.0d0
    do j=1, numberChargedParticles
        KE_f = KE_f + globalParticleList(j)%getTotalKE() + SUM(globalParticleList(j)%energyLoss) * globalParticleList(j)%mass * globalParticleList(j)%w_p * 0.5d0
    end do
    print *, ''
    print *, 'J is:'
    print *, SUM(globalSolver%J, DIM=2)
    print *, ''
    print *, 'Field half logical is:'
    print *, globalSolver%EField
    print *, 'KE_i:', KE_i
    print *, 'PE_i:', PE_i
    print *, 'KE_f:', KE_f
    print *, 'PE_f:', PE_f
    print *, 'Energy conservation:'
    print *, ABS((KE_i + PE_i - KE_f - PE_f)/(KE_i + PE_i))
    EJ = del_t * SUM(SUM(globalSolver%J, DIM = 2) * globalSolver%EField)
    print *, 'E*J:', SUM(globalSolver%J, DIM = 2) * globalSolver%EField
    print *, 'Total E*J is:', EJ
    print *, 'KE_f - KE_i:', KE_f - KE_i
    print *, 'PE_f - PE_i:', PE_f - PE_i
    print *, 'offset energy is:', PE_f + KE_f - KE_i - PE_i
    call globalSolver%evaluateEFieldCurrTime(globalWorld)
    print *, 'JV is:', (2.0d0 * SUM(globalSolver%J(1,:)) + (eps_0/del_t) * (globalSolver%EField(1) - Field_past(1))) * del_t * (globalSolver%dirichletVals(1) - globalSolver%dirichletVals(2))
    print *, 'took', iterNumPicard, 'iterations'
    call depositRho(globalSolver%rho, globalParticleList, globalWorld)
    do i = 1, NumberXNodes
        if (i == 1) then
            print *, 1.0d0 + del_t * (SUM(globalSolver%J(i+1, :)) - 2.0d0 * SUM(globalSolver%J(i,:)))/(globalSolver%rho(i) - rho_i(i))
        else if (i == NumberXNodes) then
            print *, 1.0d0 + del_t * ( - SUM(globalSolver%J(i,:)))/(globalSolver%rho(i) - rho_i(i))
        else
            print *, 1.0d0 + del_t * (SUM(globalSolver%J(i+1, :)) - SUM(globalSolver%J(i,:)))/(globalSolver%rho(i) - rho_i(i))
        end if
    end do
    print *, 'gauss error is:', globalSolver%getError_tridiag_Poisson(globalWorld)
    stop
    ! call solveSimulation(globalSolver, globalParticleList, globalWorld, del_t, maxIter, eps_r, irand, simulationTime)
    ! call solveSimulationFinalAverage(globalSolver, globalParticleList, globalWorld, del_t, maxIter, eps_r, irand, averagingTime, 100)
end program BoundPlasmaExample