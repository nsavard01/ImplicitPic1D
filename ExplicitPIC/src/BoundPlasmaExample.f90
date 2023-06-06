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

    integer(int32) :: i!, tclock1, tclock2, clock_rate
    type(Domain) :: world
    type(Particle), allocatable :: particleList(:)
    type(potentialSolver) :: solver

    !initialize the particles in this world, at some point will be read from input file or something
    particleList = readParticleInputs('BoundExample.inp',numberChargedParticles, irand, T_e) 
    ! Initialize constants with inputs
    ! create the world the particles live in
    call readInputs(NumberXNodes, numDiagnosticSteps, averagingTime, fractionFreq, n_ave, world, solver, simulationTime, Power, heatSkipSteps, nu_h, T_e, 'Geometry.inp', 'InitialConditions.inp')
    do i = 1, numberChargedParticles
        call getRandom(particleList(i)%phaseSpace(1, 1:particleList(i)%N_p), irand)
        particleList(i)%phaseSpace(1, 1:particleList(i)%N_p) = particleList(i)%phaseSpace(1, 1:particleList(i)%N_p) * real(NumberXNodes-1) + 1.0d0
        call particleList(i) % initialize_n_ave(n_ave, world%grid(NumberXNodes) - world%grid(1))
    end do
    print *, "Calulated values:"
    print *, "w_p is:", particleList(1)%w_p
    print *, "Debye length is:", getDebyeLength(particleList(1)%getTotalKE()*2.0d0/3.0d0, 2.0d0 * n_ave)
    print *, 'delX grid is:', world%delX
    print *, "Plasma frequency is:", getPlasmaFreq(2.0d0 *  n_ave)
    print *, "Average density is ", particleList(1)%N_p * particleList(1)%w_p / (world%grid(NumberXNodes) - world%grid(1)), "should be", 2.0d0 * n_ave
    del_t = fractionFreq/getPlasmaFreq(2.0d0 * n_ave)  
    print *, "rho_const is:", solver%rho_const 
    print *, "Time step (sec) is:", del_t
    print *, "----------------"
    print *, ""
    call solver%solvePotential(particleList, world)
    stop

    ! call solver%construct_diagMatrix_Ampere(world)
    call solveSimulationOnlyPotential(solver, particleList, world, del_t, simulationTime)


    

    

    
end program BoundPlasmaExample