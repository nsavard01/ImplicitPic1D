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
    real(real64) :: E_i, E_f, initialNorm, ftol, stptol
    
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
    ! Set Nitsol parameters
    iplvl = 4 ! maximum print
    iterm = 0
    kdmax = 60 ! maximum krylov subspace dimension
    input = 0
    input(1) = 50 ! maximum iterations
    input(2) = 0 !ijacv
    input(4) = kdmax
    input(5) = 0 !ipre
    input(10) = 2
    etamax = 0.8d0
    choice2_exp = 1.5d0
    choice2_coef = 0.9d0
    ftol = 1.0d-8
    stptol = 1.0d-8

    call solver%depositRho(particleList, world)
    call solver%solve_tridiag_Poisson()
    print *, solver%phi_f
    solver%phi_f = 0.0d0
    !call solver%solveInitialPotential(particleList, world)
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
    allocate(fcur(NumberXNodes-2), rwork((NumberXNodes-2)*(kdmax+5)+kdmax*(kdmax+3)))
    call funcNitsol(NumberXNodes-2, solver%phi_f(2:NumberXNodes-1), fcur, rpar, ipar, itrmf)
    initialNorm = SUM(SQRT(fcur**2))
    call nitsol(NumberXNodes-2, solver%phi_f, funcNitsol, jacNitsol, ftol*initialNorm, stptol,input, info, rwork, rpar, ipar, iterm, ddot, dnrm2)
    print *, solver%phi_f
    write(6,*) 
    write(6,880) iterm
    write(6,900) info(1)
    write(6,910) info(2)
    write(6,920) info(3)
    write(6,930) info(4)
    write(6,940) info(5)
    write(6,950) info(6)
    880  format(' Termination flag iterm:       ', i9)
    890  format(' Final f-norm:                 ', t36, 1pe9.3)
    900  format(' No. function evaluations:     ', i9)
    910  format(' No. J*v evaluations:          ', i9) 
    920  format(' No. P(inverse)*v evaluations: ', i9)
    930  format(' No. linear iterations:        ', i9)
    940  format(' No. nonlinear iterations:     ', i9)
    950  format(' No. backtracks:               ', i9)
    stop
    
    call solveSimulation(solver, particleList, world, del_t, maxIter, eps_r, irand, simulationTime, heatSkipSteps)
    print *, "Averaging over", stepsAverage, "time steps"
    call solveSimulationFinalAverage(solver, particleList, world, del_t, maxIter, eps_r, irand, stepsAverage, heatSkipSteps)
    


    
end program BoundPlasmaExample