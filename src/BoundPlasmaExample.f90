program BoundPlasmaExample
    use constants
    use iso_fortran_env, only: int32, real64, output_unit
    use mod_domain
    use mod_particle
    implicit none

    integer(int32) :: particleIdxFactor = 2
    integer(int32), parameter :: num_grid_nodes = 32, numParticles = 10000, maxIter = 50
    real(real64), parameter :: L_domain = 0.1, del_l = 0.005
    real(real64) :: w_p = 1.0
    !real(real64) :: n_i = 5.0e14, T_e = 5.0, T_i = 0.025
    type(Domain) :: world
    type(Particle) :: electron

    world = Domain(num_grid_nodes)
    call world % constructSineGrid(del_l, L_domain)

    electron = Particle(mass = m_p, q = e, w_p = w_p, N_p = numParticles, finalIdx = numParticles * particleIdxFactor)
    print *, electron%q
    print *, electron%mass
    print *, electron%w_p
    print *, electron%N_p
    print *, electron%finalIdx
    print *, size(electron%x)


    

end program BoundPlasmaExample