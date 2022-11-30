program BoundPlasmaExample
    use constants
    use mod_BasicFunctions
    use iso_fortran_env, only: int32, real64, output_unit
    use mod_domain
    use mod_particle
    implicit none

    integer(int32) :: particleIdxFactor = 2, i
    integer(int32), parameter :: num_grid_nodes = 32, numParticles = 10000, maxIter = 50
    real(real64), parameter :: L_domain = 0.1, del_l = 0.005
    real(real64) :: w_p = 1.0 !, n_ave = 5e14
    !real(real64) :: n_i = 5.0e14, T_e = 5.0, T_i = 0.025
    type(Domain) :: world
    type(Particle) :: particleList(2)
    ! character(len=50) :: name1, name2

    ! name1 = "electron"
    ! name2 = "proton"

    ! create the world the particles live in
    world = Domain(num_grid_nodes)
    call world % constructSineGrid(del_l, L_domain)

    ! initialize the particles in this world
    particleList(1) = Particle(m_e, -e, w_p, numParticles, numParticles * particleIdxFactor, "electron")
    !call electron % initialize_randUniform(L_domain, world%dx_dl, world%n_x)
    particleList(2) = Particle(m_p, e, w_p, numParticles, numParticles * particleIdxFactor, "proton")
    do i = 1, 2
        print *, particleList(i) % name
        call particleList(i) % initialize_randUniform(L_domain, world%dx_dl, world%n_x)
    end do

    print *, particleList(2)%l_p(particleList(2)%N_p-10:particleList(2)%N_p+1)

    


    

end program BoundPlasmaExample