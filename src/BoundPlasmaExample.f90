program BoundPlasmaExample
    use constants
    use iso_fortran_env, only: int32, real64, output_unit
    implicit none

    integer(int32) :: i
    integer(int32), parameter :: num_grid_nodes = 32, numParticles = 10000, maxIter = 50
    real(real64), parameter :: L_domain = 0.1, n_i = 5e14, T_e = 5
    integer(int32) :: l_grid(num_grid_nodes)
    write(output_unit, '(e10.5)') n_i

    l_grid = (/(i, i=1, num_grid_nodes)/)

    do i = 1, num_grid_nodes
        print *, l_grid(i)
    end do

    

end program BoundPlasmaExample