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
    real(real64) :: w_p = 1.0, n_ave = 5e14, T_e = 5.0
    type(Domain) :: world
    type(Particle) :: particleList(2)

    ! create the world the particles live in
    world = Domain(num_grid_nodes)
    call world % constructSineGrid(del_l, L_domain)

    ! initialize the particles in this world
    particleList(1) = Particle(m_e, -e, w_p, numParticles, numParticles * particleIdxFactor, "electron")
    particleList(2) = Particle(m_p, e, w_p, numParticles, numParticles * particleIdxFactor, "proton")
    do i = 1, 2
        print *, particleList(i) % name
        call particleList(i) % initialize_randUniform(L_domain, world%dx_dl, world%n_x)
        call particleList(i) % initialize_n_ave(n_ave, L_domain)
        print *, particleList(i) % w_p * particleList(i) % N_p / L_domain
    end do

    print *, "Debye length is:", getDebyeLength(T_e, n_ave)
    print *, "Plasma frequency is:", getPlasmaFreq(n_ave)
    


    

end program BoundPlasmaExample