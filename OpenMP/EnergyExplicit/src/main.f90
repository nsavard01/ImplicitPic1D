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
    ! use mod_nonLinSolvers
    use mod_simulation
    use omp_lib
    implicit none
    
    integer(int32) :: i, j, iThread, startTime, endTime, loop_var
    real(real64) :: remainDel_t, currDel_t, PE_i, KE_i, PE_f, KE_f, Momentum_i(3), Momentum_f(3), E_i, E_f, PE2_i, PE2_f, E2_i, E2_f
    real(real64), allocatable :: rho_i(:), J_total(:), phi_calc(:), E_calc(:)
    type(targetParticle), allocatable :: targetParticleList(:)
    type(nullCollision), allocatable :: nullCollisionList(:)


    call readInitialInputs('InitialConditions.inp')
    call initializeRandomGenerators(numThread, stateRan0, stateRanNew, .false.)
    call readWorld('Geometry.inp', globalWorld, T_e, n_ave)
    call readSolver('Geometry.inp', globalSolver, globalWorld, startSimulationTime)
    call readChargedParticleInputs('ParticleTypes.inp', stateRan0, T_e, T_i, numThread, globalWorld, globalParticleList)
    call readNeutralParticleInputs('ParticleTypes.inp', targetParticleList)
    call readNullCollisionInputs('collision.inp', nullCollisionList, globalParticleList, targetParticleList)
    ! do i = 1, numberChargedParticles
    !     globalParticleList(i)%N_p = 0
    ! end do
    call readInjectionInputs('ParticleInjection.inp', globalParticleList(1)%w_p)
    ! if (injectionBool) call injectAtBoundary(globalParticleList, T_e, T_i, irand, globalWorld, del_t, globalSolver%BFieldAngle)

    call globalSolver%initialVRewind(globalParticleList, del_t)
    call globalSolver%solvePotential(globalParticleList, globalWorld, startSimulationTime)
    
    ! allocate(phi_calc(NumberXNodes))
    ! do loop_var = 1, NumberXNodes
    !     phi_calc(loop_var) = -0.5 * E * n_ave * globalWorld%grid(loop_var) * (globalWorld%grid(NumberXNodes) - globalWorld%grid(loop_var)) / eps_0
    ! end do
    ! print*,"PHI: ", globalSolver%phi
    ! print*,"PHI_CALC: ", phi_calc
    ! print*,"DIFFRENCE: ", globalSolver%phi - phi_calc
    ! print*,"RHO: ", globalSolver%rho/globalWorld%dx_dl(1)
    ! print*,"dx_dl: ", globalWorld%dx_dl


    ! allocate(E_calc(NumberXHalfNodes))
    ! do loop_var = 1, NumberXHalfNodes
    !     E_calc(loop_var) = -e * n_ave * (globalWorld%grid(loop_var) - (globalWorld%grid(NumberXNodes))/2) / eps_0
    ! end do
    ! print*,"E: ", globalSolver%EField
    ! print*,"PHI_CALC: ", E_calc
    ! print*,"DIFFRENCE: ", (globalSolver%EField - E_calc)/globalSolver%EField
    ! stop

    ! print*,"A_TRI: ", globalSolver%a_tri
    ! print*,"B_TRI: ", globalSolver%b_tri
    ! print*,"C_TRI: ", globalSolver%c_tri
    !stop





    ! print *, 'electron particle number:', SUM(globalParticleList(1)%N_p)
    ! print *, 'ion particle number:', SUM(globalParticleList(2)%N_p)
    ! PE_i = globalSolver%getTotalPE(globalWorld)
    ! PE2_i = globalSolver%getTotalPE2(globalWorld)
    ! KE_i = globalSolver%getTotalKE(del_t, globalParticleList,globalWorld)
    ! E_i = PE_i + KE_i
    ! E2_i = PE2_i + KE_i

    ! print*, globalWorld%grid(1)
    ! call moveParticles(globalSolver, globalParticleList, globalWorld, del_t)
    ! call globalSolver%solvePotential(globalParticleList, globalWorld, startSimulationTime + del_t)
    ! print *, 'electron particle number:', SUM(globalParticleList(1)%N_p)
    ! print *, 'ion particle number:', SUM(globalParticleList(2)%N_p)
    ! stop
    ! PE_f = globalSolver%getTotalPE(globalWorld)
    ! PE2_f = globalSolver%getTotalPE2(globalWorld)
    ! KE_f = globalSolver%getTotalKE(del_t, globalParticleList, globalWorld)
    ! do j=1, numberChargedParticles
    !     KE_f = KE_f + SUM(globalParticleList(j)%energyLoss) * globalParticleList(j)%mass * globalParticleList(j)%w_p * 0.5d0
    ! end do
    ! E_f = PE_f + KE_f
    ! E2_f = PE2_f + KE_f
    ! print *, 'PE_i, PE_f: ', PE_i, PE_f
    ! print *, 'PE2_i, PE2_f: ', PE2_i, PE2_f
    ! print *, 'KE_i, KE_f: ', KE_i, KE_f
    ! print *, 'E_i, E2_i: ', E_i, E2_i
    ! print *, 'E_f, E2_f: ', E_f, E2_f
    ! print *, 'E error: ', E_f - E_i
    ! print *, 'E2 error: ', E2_f - E2_i
    ! print *, 'E error adjusted: ', ABS((E_i - E_f)/(E_i))
    ! print *, 'E2 error adjusted: ', ABS((E2_i - E2_f)/(E2_i))
    ! stop
    
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
    ! call system_clock(endTime)
    
    ! ! do i = 1, NumberXHalfNodes
    ! !     J_total(i) = eps_0 * (((globalSolver%phi_f(i) - globalSolver%phi_f(i+1)) -(globalSolver%phi(i) - globalSolver%phi(i+1))) /globalWorld%dx_dl(i)) / currDel_t + globalSolver%J(i)
    ! ! end do
    ! ! print *, J_total
    ! ! print *, ''
    ! ! print *, 'S is:', (SUM(globalSolver%J * globalWorld%dx_dl) + (eps_0/currDel_t) * &
    ! !     (globalSolver%phi_f(NumberXNodes) - globalSolver%phi_f(1) - globalSolver%phi(NumberXNodes) + globalSolver%phi(1)))/globalWorld%L_domain
   
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
    call solveSimulation(globalSolver, globalParticleList, targetParticleList, nullCollisionList, globalWorld, del_t, stateRan0, simulationTime)
    call solveSimulationFinalAverage(globalSolver, globalParticleList, targetParticleList, nullCollisionList, globalWorld, del_t, stateRan0, averagingTime, 100)
end program BoundPlasmaExample