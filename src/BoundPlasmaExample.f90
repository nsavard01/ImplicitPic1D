program BoundPlasmaExample
    use constants
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
    type(Particle) :: electron
    type(Particle) :: ion
    type(Particle) :: particleList(2)

    ! create the world the particles live in
    world = Domain(num_grid_nodes)
    call world % constructSineGrid(del_l, L_domain)

    ! initialize the particles in this world
    electron = Particle(mass = m_e, q = -e, w_p = w_p, N_p = numParticles, finalIdx = numParticles * particleIdxFactor)
    call electron % initialize_randUniform(L_domain, world%dx_dl, world%n_x)
    ion = Particle(mass = m_p, q = e, w_p = w_p, N_p = numParticles, finalIdx = numParticles * particleIdxFactor)
    particleList = [electron, ion]
    do i = 1, 2
        print *, particleList(i) % mass
    end do

    print *, electron%l_p(electron%N_p-10:electron%N_p+1)

    


    

end program BoundPlasmaExample