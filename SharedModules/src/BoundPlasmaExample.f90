program BoundPlasmaExample
    use iso_fortran_env, only: int32, real64, output_unit, int64
    use constants
    use mod_BasicFunctions
    use mod_domain
    use mod_particle
    use mod_collisions
    use mod_potentialSolver
    use mod_particleMover
    use mod_nonLinSolvers
    ! use mod_Scheme
    use mod_simulation
    implicit none

    integer(int32) :: i
    ! type(Domain) :: globalWorld
    ! type(potentialSolver) :: globalSolver
    ! type(Particle), allocatable :: globalParticleList(:)
    real(real64) :: currDel_t, remainDel_t, E_i, E_f
    call readInitialConditions('InitialConditions.inp')
    globalWorld = Domain('Geometry.inp')
    globalParticleList =  readParticleInputs('BoundExample.inp',numberChargedParticles, irand)
    do i = 1, numberChargedParticles
        call globalParticleList(i)%initialize_randUniform(globalWorld, irand)
    end do
    globalSolver = potentialSolver('Geometry.inp', globalWorld)
    call initializeSolver(eps_r, solverType, m_Anderson, Beta_k, maxIter)
    print *, ""
    print *, "Calulated values:"
    print *, "Number of particles is:", SUM(globalParticleList(1)%N_p)
    print *, "Debye length is:", getDebyeLength(globalParticleList(1)%getKEAve()*2.0d0/3.0d0, n_ave)
    print *, "Plasma frequency is:", getPlasmaFreq(n_ave)
    del_t = fractionFreq/getPlasmaFreq(n_ave)   
    print *, "Time step (sec) is:", del_t
    print *, "----------------"
    print *, ""
    call solveInitialPotential(globalSolver, globalParticleList, globalWorld)
    ! currDel_t = del_t
    ! remainDel_t = del_t
    ! E_i = globalSolver%getTotalPE(globalWorld, .false.)
    ! do i=1, numberChargedParticles
    !     E_i = E_i + globalParticleList(i)%getTotalKE()
    ! end do
    ! call solvePotential(globalSolver, globalParticleList, globalWorld, del_t, remainDel_t, currDel_t, maxIter, eps_r)
    ! E_f = globalSolver%getTotalPE(globalWorld, .false.)
    ! do i=1, numberChargedParticles
    !     E_f = E_f + globalParticleList(i)%getTotalKE() + SUM(globalParticleList(i)%energyLoss)
    ! end do
    ! print *, "solve in:", iterNumPicard, "iterations"
    ! print *, "Energy error is:", ABS((E_f - E_i)/E_i)
    ! do i = 1, numberChargedParticles
    !     print *, "For particles", globalParticleList(i)%name
    !     print *, 'Number of refluxed particles:', globalParticleList(i)%refIdx
    !     print *, 'Number of wall lost particles:', globalParticleList(i)%delIdx
    !     print *, 'Number of particles still in:', globalParticleList(i)%N_p
    !     print *, 'Total weight loss:', SUM(globalParticleList(i)%wallLoss)
    ! end do
    ! call addMaxwellianLostParticles(globalParticleList, T_e, T_i, irand, globalWorld)
    ! globalSolver%rho = 0.0d0
    ! do i =1, numberChargedParticles
    !     call globalParticleList(i)%depositRho(globalSolver%rho, globalWorld)
    ! end do
    ! print *, "gauss error is:", globalSolver%getError_tridiag_Poisson(globalWorld)
    
    call solveSimulationTest(globalSolver, globalParticleList, globalWorld, del_t, maxIter, eps_r, irand, simulationTime, heatSkipSteps)

    
end program BoundPlasmaExample