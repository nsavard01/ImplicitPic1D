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

    integer(int32) :: i, tclock1, tclock2, clock_rate, j, k
    real(real64) :: elapsed_time
    type(Domain) :: world
    type(Particle), allocatable :: particleList(:)
    type(potentialSolver) :: solver
    real(real64), allocatable :: rho_f(:), rho_i(:)
    
    !initialize the particles in this world, at some point will be read from input file or something
    particleList = readParticleInputs('BoundExample.dat',numberChargedParticles, irand) 
    ! Initialize constants with inputs
    ! create the world the particles live in
    call readInputs(NumberXNodes, maxIter, numDiagnosticSteps, stepsAverage, eps_r, fractionFreq, n_ave, world, solver)
    do i = 1, numberChargedParticles
        call particleList(i) % initialize_randUniform(world, irand)
        call particleList(i) % initialize_n_ave(n_ave, world%grid(NumberXNodes) - world%grid(1))
    end do

    print *, "Calulated values:"
    print *, "w_p is:", particleList(1)%w_p
    print *, "Debye length is:", getDebyeLength(particleList(1)%getKEAve()*2.0d0/3.0d0, n_ave)
    print *, "Plasma frequency is:", getPlasmaFreq(n_ave)
    print *, "Average density is ", particleList(1)%N_p * particleList(1)%w_p / (world%grid(NumberXNodes) - world%grid(1)), "should be", n_ave
    del_t = fractionFreq/getPlasmaFreq(n_ave)   
    print *, "Time step (sec) is:", del_t
    print *, "----------------"
    print *, ""
    ! ! Generate solver object, and then solve for initial rho/potential
    allocate(rho_f(NumberXNodes), rho_i(NumberXNodes))
    call solver%solveInitialPotential(particleList, world)
    call depositRhoDiag(rho_i, particleList, world)
    call solver%depositJ(particleList, world, del_t, maxIter, eps_r)
    call depositRhoDiag(rho_f, particleList, world)
    solver%chargeError = 0.0d0
    j = 0
    do k = 1, NumberXNodes -2
        if ((solver%J(k + 1) - solver%J(k) /= 0) .and. (rho_f(k+1) - rho_i(k+1))/= 0) then
            solver%chargeError = solver%chargeError + (1 + (solver%J(k + 1) - solver%J(k)) *del_t/ world%dx_dl(k+1)/(rho_f(k+1) - rho_i(k+1)))**2
            j = j + 1
        end if
    end do
    print *, 1 + (solver%J(2:) - solver%J(1:NumberXNodes-2)) *del_t/ world%dx_dl(2:NumberXNodes-1)/(rho_f(2:NumberXNodes-2) - rho_i(2:NumberXNodes-2))
    solver%chargeError = SQRT(solver%chargeError/real(j, kind = real64))
    print *, "Charge error is:", solver%chargeError
    stop
    numTimeSteps = NINT(22.0d-6 / del_t)
    call solveSingleTimeStepDiagnostic(solver, particleList, world, del_t, maxIter, eps_r)
    call system_clock(tclock1)
    call solveSimulation(solver, particleList, world, del_t, maxIter, eps_r, irand, numTimeSteps)
    call system_clock(tclock2, clock_rate)
    elapsed_time = float(tclock2 - tclock1) / float(clock_rate)
    print *, "Elapsed time for simulation is:", elapsed_time/60.0d0, "minutes"
    

    

    
end program BoundPlasmaExample