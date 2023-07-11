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
    real(real64) :: del_t, currDel_t, remainDel_t, E_i, E_f, l_f_test, l_del
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
    currDel_t = del_t
    remainDel_t = del_t
    E_i = globalSolver%getTotalPE(globalWorld, .false.)
    do i=1, numberChargedParticles
        E_i = E_i + globalParticleList(i)%getTotalKE()
    end do
    ! print *, 'electron total weights before:', globalParticleList(1)%getSumWeights()
    call solvePotential(globalSolver, globalParticleList, globalWorld, del_t, remainDel_t, currDel_t, maxIter, eps_r)
    ! print *, 'electron total weights after:', globalParticleList(1)%getSumWeights() + SUM(globalParticleList(1)%wallLoss)
    ! globalParticleList(1)%N_p(1) = globalParticleList(1)%N_p(1) + 1
    ! globalParticleList(1)%phaseSpace(:, globalParticleList(1)%N_p(1), 1) = globalParticleList(1)%phaseSpace(:, 50, 1)
    ! globalParticleList(1)%w_p(50, 1) = globalParticleList(1)%w_p(50, 1)/2.0d0
    ! globalParticleList(1)%w_p(globalParticleList(1)%N_p(1), 1) = globalParticleList(1)%w_p(50, 1)
    ! l_del = MIN(ABS(globalParticleList(1)%phaseSpace(1, 50, 1) - 1), ABS(globalParticleList(1)%phaseSpace(1, 50, 1) - 2))
    ! l_del = l_del - 0.01
    ! print *, 'OG position:', globalParticleList(1)%phaseSpace(1, 50, 1), globalParticleList(1)%phaseSpace(1, globalParticleList(1)%N_p(1), 1)
    ! globalParticleList(1)%phaseSpace(1, 50, 1) = globalParticleList(1)%phaseSpace(1, 50, 1) + l_del
    ! globalParticleList(1)%phaseSpace(1, globalParticleList(1)%N_p(1), 1) = globalParticleList(1)%phaseSpace(1, globalParticleList(1)%N_p(1), 1) - l_del
    ! print *, 'final positions:', globalParticleList(1)%phaseSpace(1, 50, 1), globalParticleList(1)%phaseSpace(1, globalParticleList(1)%N_p(1), 1)
    E_f = globalSolver%getTotalPE(globalWorld, .false.)
    do i=1, numberChargedParticles
        E_f = E_f + globalParticleList(i)%getTotalKE() + SUM(globalParticleList(i)%energyLoss)
    end do
    print *, "solve in:", iterNumPicard, "iterations"
    print *, "Energy error is:", ABS((E_f - E_i)/E_i)
    do i = 1, numberChargedParticles
        print *, "For particles", globalParticleList(i)%name
        print *, 'Number of refluxed particles:', globalParticleList(i)%refIdx
        print *, 'Number of wall lost particles:', globalParticleList(i)%delIdx
        print *, 'Number of particles still in:', globalParticleList(i)%N_p
        print *, 'Total weight loss:', SUM(globalParticleList(i)%wallLoss)
    end do
    call addMaxwellianLostParticles(globalParticleList, T_e, T_i, irand, globalWorld)
    globalSolver%rho = 0.0d0
    do i =1, numberChargedParticles
        call globalParticleList(i)%depositRho(globalSolver%rho, globalWorld)
    end do
    print *, "gauss error is:", globalSolver%getError_tridiag_Poisson(globalWorld)
    
    !real(real64) :: remainDel_t, currDel_t
    ! call initializeScheme(schemeNum)
    ! ! Initialize constants with inputs
    ! ! create the world the particles live in
    ! call readInputs(NumberXNodes, numDiagnosticSteps, averagingTime, fractionFreq, n_ave, globalWorld, globalSolver, simulationTime, Power, heatSkipSteps, nu_h, T_e, T_i, 'Geometry.inp', 'InitialConditions.inp')
    ! globalParticleList = readParticleInputs('BoundExample.inp',numberChargedParticles, irand, T_e, T_i) 
    ! do i = 1, numberChargedParticles
    !     call initialize_randUniform(globalParticleList(i), globalWorld, irand)
    !     call globalParticleList(i) % initialize_n_ave(n_ave, globalWorld%grid(NumberXNodes) - globalWorld%grid(1))
    ! end do
    ! call initializeSolver(eps_r, solverType, m_Anderson, Beta_k, maxIter)
    ! print *, "Calulated values:"
    ! print *, "Number of particles is:", globalParticleList(1)%N_p
    ! print *, "w_p is:", globalParticleList(1)%w_p
    ! print *, "Debye length is:", getDebyeLength(globalParticleList(1)%getKEAve()*2.0d0/3.0d0, n_ave)
    ! print *, "Plasma frequency is:", getPlasmaFreq(n_ave)
    ! print *, "Average density is ", globalParticleList(1)%N_p * globalParticleList(1)%w_p / (globalWorld%grid(NumberXNodes) - globalWorld%grid(1)), "should be", n_ave
    ! del_t = fractionFreq/getPlasmaFreq(n_ave)   
    ! print *, "Time step (sec) is:", del_t
    ! print *, "----------------"
    ! ! Generate solver object, and then solve for initial rho/potential
    ! call solveInitialPotential(globalSolver, globalParticleList, globalWorld)
    ! ! currDel_t = del_t
    ! ! remainDel_t = del_t
    ! ! call solvePotential(globalSolver, globalParticleList, globalWorld, del_t, remainDel_t, currDel_t, maxIter, eps_r)
    ! ! print *, "Took", iterNumPicard, "iterations"
    ! ! print *, 'electron number:', globalParticleList(1)%N_P
    ! ! print *, 'ion number:', globalParticleList(2)%N_p
    ! ! call addMaxwellianLostParticles(globalParticleList, T_e, T_i, irand, globalWorld)
    ! ! print *, 'electron number:', globalParticleList(1)%N_P
    ! ! print *, 'ion number:', globalParticleList(2)%N_p
    ! ! ! Get error gauss' law
    ! ! call depositRho(globalSolver%rho, globalParticleList, globalWorld)
    ! ! call globalSolver%construct_diagMatrix(globalWorld)
    ! ! chargeError = globalSolver%getError_tridiag_Poisson(globalWorld)
    ! ! !chargeError = chargeError / SQRT(SUM(globalSolver%rho**2))
    ! ! call globalSolver%construct_diagMatrix_Ampere(globalWorld)
    ! ! print *, "Charge error is:", chargeError
    ! ! stop 
    ! call solveSimulation(globalSolver, globalParticleList, globalWorld, del_t, maxIter, eps_r, irand, simulationTime, heatSkipSteps)

    ! print *, "Averaging up to", averagingTime, "simulation seconds"
    ! call solveSimulationFinalAverage(globalSolver, globalParticleList, globalWorld, del_t, maxIter, eps_r, irand, averagingTime, 100)

    
end program BoundPlasmaExample