program BoundPlasmaExample
    use iso_fortran_env, only: int32, real64, output_unit, int64
    use constants
    use mod_BasicFunctions
    use mod_domain
    use mod_particle
    use mod_potentialSolver
    implicit none

    type(Domain) :: world
    type(potentialSolver) :: solver
    integer(int32) :: gridType, i
    boolCIC = .true.
    gridType = 1
    NumberXNodes = 32
    world = Domain(1, 1)
    call world%constructGrid(7.0d-4, 0.1d0, gridType)
    print *, ""
    solver = potentialSolver(world, 0.0d0, 0.0d0)
    print *, "a_tri:"
    print *, solver%a_tri
    print *, "b_tri:"
    print *, solver%b_tri
    print *, "c_tri:"
    print *, solver%c_tri
    solver%rho = 1.0d12 * e
    call solver%solve_tridiag_Poisson(world)
    print *, "Discrete phi:"
    print *, solver%phi
    print *, "error is:"
    print *, solver%getError_tridiag_Poisson(world)
    print *, "Analytical phi:"
    print *, (solver%rho/2.0/eps_0) * world%grid * (world%grid(NumberXNodes) - world%grid)
    
    

    
end program BoundPlasmaExample