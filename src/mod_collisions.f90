module mod_collisions
    use iso_fortran_env, only: int32, real64
    use constants
    use mod_BasicFunctions
    use mod_particle
    use mod_domain
    implicit none

    real(real64) :: inelasticEnergyLoss = 0.0d0
    real(real64) :: Power = 100 !W/m^2 in 1D
    real(real64) :: nu_h = 1e8 !Hz
    ! Type of collisions, will likely need arrays which you can loop through which has source, target particle, choose from particle list
    ! Will likely have construct method which takes in collision data and turns into array

contains

    ! ----------------------- For 1D example, keep simple for now --------------------------------------
    subroutine ionizationCollisionIsotropic(electron, ion, n_g, sigma, del_t, E_thres, irand, T_gas)
        ! Temporary subroutine for artificial ionization, using MCC, energy units in eV
        ! divide energy among electrons, isotropic scatter
        type(Particle), intent(in out) :: electron, ion
        real(real64), intent(in) :: n_g, sigma, del_t, E_thres, T_gas
        integer(int32), intent(in out) :: irand
        real(real64) :: electronEnergy, electronSpeed, Pcoll, R, energyPerElectron, speedPerElectron, U, theta, e_vector(3), V_neutral(3), V_cm(3)
        integer(int32) :: i, addIdx
        addIdx = 0
        do i=1, electron%N_p
            electronEnergy = SUM(electron%v_p(i, :)**2) * 0.5d0 * electron%mass / e ! in eV
            if (electronEnergy > E_thres) then
                R = ran2(irand)
                electronSpeed = SQRT(SUM(electron%v_p(i, :)**2))
                Pcoll = 1.0d0 - EXP(-n_g * sigma * del_t * electronSpeed)
                if (Pcoll > R) then
                    addIdx = addIdx + 1
                    energyPerElectron = (electronEnergy - E_thres)/2.0d0
                    speedPerElectron = SQRT(2.0d0 * e * energyPerElectron/electron%mass)
                    call getMaxwellianSample(V_neutral, ion%mass, T_gas, irand)
                    V_cm = (m_e *electron%v_p(i, :) + ion%mass * V_neutral) / (ion%mass)

                    U = 1.0d0 - 2.0d0*ran2(irand)
                    theta = ran2(irand) * 2 * pi
                    e_vector(1) = SQRT(1-U**2) * COS(theta) ! unit normal in x direction
                    e_vector(2) = SQRT(1 - U**2) * SIN(theta)
                    e_vector(3) = U
                    electron%v_p(i, :) = e_vector * speedPerElectron

                    U = 1.0d0 - 2.0d0*ran2(irand)
                    theta = ran2(irand) * 2 * pi
                    e_vector(1) = SQRT(1-U**2) * COS(theta) ! unit normal in x direction
                    e_vector(2) = SQRT(1 - U**2) * SIN(theta)
                    e_vector(3) = U
                    electron%v_p(electron%N_p + addIdx, :) = e_vector * speedPerElectron
                    electron%l_p(electron%N_p + addIdx) = electron%l_p(i)

                    ! Calculate new ion velocity
                    ion%v_p(ion%N_p + addIdx, :) = -(m_e/ion%mass) * (electron%v_p(i, :) + electron%v_p(electron%N_p + addIdx, :)) + V_cm
                    ion%l_p(ion%N_p + addIdx) = electron%l_p(i)
                    
                    

                end if
            end if


        end do

        electron%N_p = electron%N_p + addIdx
        ion%N_p = ion%N_p + addIdx
        inelasticEnergyLoss = inelasticEnergyLoss + E_thres * addIdx
    end subroutine ionizationCollisionIsotropic

    subroutine getMaxwellianSample(v, M, T_gas, irand)
        ! Use to sample a single neutral background particle at some temperature T_gas (eV)
        real(real64), intent(in out) :: v(3)
        real(real64), intent(in) :: T_gas, M
        integer(int32), intent(in out) :: irand
        real(real64) :: U(4)
        call getRandom(U, irand)
        v(1) = SQRT(T_gas*e/ M) * SQRT(-2 * LOG(U(1))) * COS(2 * pi * U(2))
        v(2) = SQRT(T_gas*e/ M) * SQRT(-2 * LOG(U(1))) * SIN(2 * pi * U(2))
        v(3) = SQRT(T_gas*e/ M) * SQRT(-2 * LOG(U(3))) * SIN(2 * pi * U(4))

    end subroutine getMaxwellianSample


    !-------------------------- add Power in form of re-maxwellianizing collision ---------------------------------------------
    subroutine addUniformPowerMaxwellian(species, Power, nu_h, irand, del_t)
        ! Add power to all particles in domain
        type(Particle), intent(in out) :: species
        integer(int32), intent(in out) :: irand
        real(real64), intent(in) :: Power, nu_h, del_t
        real(real64) :: T_h, idxChange(NINT(nu_h * del_t * species%N_p)), v_replace(3)
        integer(int32) :: j, i
        T_h = (species%getKEAve() + Power/e/species%N_p/species%w_p/nu_h) * 2.0d0 / 3.0d0
        call getRandom(idxChange, irand)
        idxChange = INT(idxChange*(species%N_p-1) + 1)
        do i=1, size(idxChange)
            j = idxChange(i)
            call getMaxwellianSample(v_replace, species%mass, T_h, irand)
            species%v_p(j, :) = v_replace
        end do
    end subroutine addUniformPowerMaxwellian

    subroutine addUniformPowerMaxwellianNicolas(species, Power, percent, irand, del_t)
        ! Gwenael's method does not lead to conservation, very shady
        ! Implement uniform power with conservation properties in mind (still statistical) to get improvement
        ! Percent is percentage of particles to sample, isotropic scatter
        type(Particle), intent(in out) :: species
        integer(int32), intent(in out) :: irand
        real(real64), intent(in) :: Power, percent, del_t
        integer(int32) :: j, i
        real(real64) :: KE_new, idxChange(NINT(percent * species%N_p)), v_replace(3), E_distribute, U, theta
        call getRandom(idxChange, irand)
        idxChange = INT(idxChange*(species%N_p-1) + 1)
        E_distribute = Power*del_t/e/species%w_p/size(idxChange)
        do i=1, size(idxChange)
            j = idxChange(i)
            KE_new = (SUM(species%v_p(j,:)**2) * species%mass * 0.5/e + E_distribute)
            U = 1.0d0 - 2.0d0*ran2(irand)
            theta = ran2(irand) * 2 * pi
            v_replace(1) = SQRT(1-U**2) * COS(theta) ! unit normal in x direction
            v_replace(2) = SQRT(1 - U**2) * SIN(theta)
            v_replace(3) = U
            species%v_p(j, :) = v_replace * SQRT(2.0d0 * e * KE_new/species%mass)
        end do
    end subroutine addUniformPowerMaxwellianNicolas



end module mod_collisions