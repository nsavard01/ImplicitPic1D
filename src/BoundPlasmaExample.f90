program BoundPlasmaExample
    use constants
    use iso_fortran_env, only: int32, real64
    implicit none

    integer(int32) :: num_grid_nodes, numParticles, maxIter
    num_grid_nodes = 32
    numParticles = 10000
    maxIter = 50

    print *, num_grid_nodes
    print *, numParticles
    

end program BoundPlasmaExample