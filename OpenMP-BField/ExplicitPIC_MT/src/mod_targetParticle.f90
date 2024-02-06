module mod_targetParticle

    use iso_fortran_env, only: int32, real64, output_unit
    use constants
    use mod_BasicFunctions
    use mod_mt19937
    implicit none

    private
    public :: targetParticle

    ! target particle contains basic properties of basic neutral particle for target null collision
    type :: targetParticle
        character(:), allocatable :: name !name of the particle
        real(real64) :: mass, density, temperature, v_therm
    contains
        procedure, public, pass(self) :: generate3DMaxwellianVelocity
    end type targetParticle


    interface targetParticle
        module procedure :: targetParticle_constructor
    end interface targetParticle

contains

    type(targetParticle) function targetParticle_constructor(particleName, mass, density, temperature) result(self)
        ! Construct particle object, sizeIncrease is fraction larger stored array compared to initial amount of particles
        ! In future, use hash function for possible k = 1 .. Nx, m amount of boundaries, p = prime number  m < p < N_x. h(k) = (k%p)%m
        real(real64), intent(in) :: mass, density, temperature
        character(*), intent(in) :: particleName
        self % name = particleName
        self % mass = mass
        self % density = density ! #N/m^3
        self % temperature = temperature !temperature
        self % v_therm = SQRT(self%temperature*k_B/ self%mass)
    end function targetParticle_constructor


    ! ------------------------ Generating and initializing velocity with Maxwellian --------------------------


    function generate3DMaxwellianVelocity(self, randGen) result(res)
        ! random velocity generator for the particle for temperature T (eV)
        ! Use box-muller method for random guassian variable, same as gwenael but doesn't have factor 2? Maybe factored into v_th
        class(targetParticle), intent(in) :: self
        type(mt19937), intent(in out) :: randGen
        real(real64) :: U1, U2, U3, U4, res(3)
        U1 = randGen%genrand64_real1()
        U2 = randGen%genrand64_real1()
        U3 = randGen%genrand64_real1()
        U4 = randGen%genrand64_real1()
        res(1) = self%v_therm * SQRT(-2 * LOG(U1)) * COS(2 * pi * U2)
        res(2) = self%v_therm * SQRT(-2 * LOG(U1)) * SIN(2 * pi * U2)
        res(3) = self%v_therm * SQRT(-2 * LOG(U3)) * SIN(2 * pi * U4)
    end function generate3DMaxwellianVelocity


end module mod_targetParticle