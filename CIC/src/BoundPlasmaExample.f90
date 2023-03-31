program BoundPlasmaExample
    use iso_fortran_env, only: int32, real64, output_unit
    use constants
    use mod_BasicFunctions
    use mod_domain
    use mod_particle
    use mod_potentialSolver
    use mod_collisions
    use mod_simulation
    implicit none

    integer(int32) :: i
    real(real64) :: E_i, E_f
    type(Domain) :: world
    type(Particle), allocatable :: particleList(:)
    type(potentialSolver) :: solver
    
    particleList = readParticleInputs('BoundExample.dat',numberChargedParticles, irand) 
    ! Initialize constants with inputs
    ! create the world the particles live in
    call readInputs(NumberXNodes, maxIter, numDiagnosticSteps, stepsAverage, eps_r, fractionFreq, n_ave, world, solver, simulationTime, Power, heatSkipSteps, nu_h)
    do i = 1, numberChargedParticles
        call initialize_randUniform(particleList(i), world, irand)
        call particleList(i) % initialize_n_ave(n_ave, world%grid(NumberXNodes) - world%grid(1))
    end do

    print *, "Calulated values:"
    print *, "Number of particles is:", particleList(1)%N_p
    print *, "w_p is:", particleList(1)%w_p
    print *, "Debye length is:", getDebyeLength(particleList(1)%getKEAve()*2.0d0/3.0d0, n_ave)
    print *, "Plasma frequency is:", getPlasmaFreq(n_ave)
    print *, "Average density is ", particleList(1)%N_p * particleList(1)%w_p / (world%grid(NumberXNodes) - world%grid(1)), "should be", n_ave
    del_t = fractionFreq/getPlasmaFreq(n_ave)   
    print *, "Time step (sec) is:", del_t
    print *, "----------------"
    print *, ""
    ! ! Generate solver object, and then solve for initial rho/potential
    
    call solver%solveInitialPotential(particleList, world)
    ! call solveSingleTimeStepDiagnostic(solver, particleList, world, del_t, maxIter, eps_r)
    ! print *, "Total number of iterations is:", solver%iterNumPicard
    ! print *, "Anderson number is:", solver%m_Anderson
    ! print *, "Relaxation parameter is:", solver%Beta_k
    ! print *, "Energy error is:", solver%energyError
    ! print *, "Charge error is:", solver%chargeError
    ! if (solver%energyError > eps_r) then
    !     print *, "-------------------------WARNING------------------------"
    !     print *, "Energy error is:", solver%energyError
    !     stop "Total energy not conserved over time step in sub-step procedure!"
    ! end if
    ! if (solver%chargeError > eps_r) then
    !     print *, "-------------------------WARNING------------------------"
    !     print *, "Charge error is:", solver%chargeError
    !     stop "Total charge not conserved over time step in sub-step procedure!"
    ! end if
    ! stop
    
    E_i = particleList(1)%getTotalKE() + particleList(2)%getTotalKE() + solver%getTotalPE(world, .false.)
    call solver%solveDivAmpereAnderson(particleList, world, del_t, maxIter, eps_r)
    E_f = particleList(1)%getTotalKE() + particleList(2)%getTotalKE() + solver%getTotalPE(world, .true.) + solver%particleEnergyLoss
    print *, "took iteration number:", solver%iterNumPicard
    print *, "Percent difference is:", ABS((E_f - E_i)/E_i) * 100.0d0
    stop
    call solveSimulation(solver, particleList, world, del_t, maxIter, eps_r, irand, simulationTime, heatSkipSteps)
    print *, "Averaging over", stepsAverage, "time steps"
    call solveSimulationFinalAverage(solver, particleList, world, del_t, maxIter, eps_r, irand, stepsAverage, heatSkipSteps)
    

    

    
end program BoundPlasmaExample