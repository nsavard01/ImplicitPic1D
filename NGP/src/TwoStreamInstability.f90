program TwoStreamInstability
    use iso_fortran_env, only: int32, real64, output_unit
    use constants
    use mod_BasicFunctions
    use mod_domain
    use mod_particle
    use mod_potentialSolver
    use mod_collisions
    use mod_simulation
    implicit none

    integer(int32) :: i!, tclock1, tclock2, clock_rate
    real(real64) :: v_init = 2000.0d0!elapsed_time
    type(Domain) :: world
    type(Particle), allocatable :: particleList(:)
    type(potentialSolver) :: solver

    !initialize the particles in this world, at some point will be read from input file or something
    particleList = readParticleInputs('TwoStream.dat',numberChargedParticles, irand) 
    ! Initialize constants with inputs
    ! create the world the particles live in
    call readInputs(NumberXNodes, maxIter, numDiagnosticSteps, stepsAverage, eps_r, fractionFreq, n_ave, world, solver)
    do i = 1, numberChargedParticles
        call particleList(i) % initialize_randUniform(world%grid(NumberXNodes) - world%grid(1), world%dx_dl, irand)
        call particleList(i) % initialize_n_ave(n_ave, world%grid(NumberXNodes) - world%grid(1))
    end do

    print *, "Calulated values:"
    print *, "w_p is:", particleList(1)%w_p
    print *, "Debye length is:", getDebyeLength(particleList(1)%getTotalKE()*2.0d0/3.0d0, n_ave)
    print *, "Plasma frequency is:", getPlasmaFreq(n_ave)
    print *, "Average density is ", particleList(1)%N_p * particleList(1)%w_p / (world%grid(NumberXNodes) - world%grid(1)), "should be", n_ave
    del_t = fractionFreq/getPlasmaFreq(n_ave)   
    print *, "Time step (sec) is:", del_t
    print *, "----------------"
    print *, ""
    
    do i=1, particleList(1)%N_p
        if (MOD(i,2) /= 0) then
            particleList(1)%phaseSpace(2,i) = v_init
        else
            particleList(1)%phaseSpace(2,i) = -v_init
        end if
    end do
    call solver%depositRho(particleList, world)
  
    solver%rho = solver%rho + n_ave*e
    call solver%solve_tridiag_Poisson()
    ! Assume only use potential solver once, then need to generate matrix for Div-Ampere
    call solver%construct_diagMatrix_Ampere(world)
    numTimeSteps = NINT(10.0d0 * 2.0d0 * pi /getPlasmaFreq(n_ave)/del_t)
    call solveSimulationOnlyPotential(solver, particleList, world, del_t, maxIter, eps_r, numTimeSteps, stepsAverage)
            


end program TwoStreamInstability