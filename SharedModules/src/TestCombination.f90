program BoundPlasmaExample
    use iso_fortran_env, only: int32, real64, output_unit, int64
    use constants
    use mod_BasicFunctions
    use mod_domain
    implicit none

    type(Domain) :: world
    integer(int32) :: gridType
    boolCIC = .true.
    gridType = 2
    NumberXNodes = 32
   
    world = Domain(1, 1)
    call world%constructGrid(7.0d-4, 0.1d0, gridType)
    print *, "grid is:"
    print *, world%grid
    print *, "dx_dl is:"
    print *, world%dx_dl
    print *, "sum of dx_dl is:"
    print *, SUM(world%dx_dl)
    print *, "nodeVol is:"
    print *, world%nodeVol
    

    
end program BoundPlasmaExample