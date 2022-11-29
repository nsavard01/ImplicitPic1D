program BoundPlasmaExample
    use constants
    use iso_fortran_env, only: int32, real64, output_unit
    use mod_domain
    implicit none

    !integer(int32) :: i
    integer(int32), parameter :: num_grid_nodes = 32, numParticles = 10000, maxIter = 50
    real(real64), parameter :: L_domain = 0.1, n_i = 5e14, T_e = 5, del_l = 0.005
    type(Domain) :: world
    world = Domain(num_grid_nodes)
    call world % constructUniformGrid(L_domain)
    print *, world%grid
    print *, world%dx_dl
    print *, world%nodeVol


    

end program BoundPlasmaExample