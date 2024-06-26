module mod_targetParticle

    use iso_fortran_env, only: int32, real64, output_unit
    use constants
    implicit none

    private
    public :: targetParticle

    ! target particle contains basic properties of basic neutral particle for target null collision
    type :: targetParticle
        character(:), allocatable :: name !name of the particle
        real(real64) :: mass, q, density, temperature
    contains
        procedure, public, pass(self) :: generate3DMaxwellianVelocity
    end type targetParticle


    interface targetParticle
        module procedure :: targetParticle_constructor
    end interface targetParticle

contains

    type(targetParticle) function targetParticle_constructor(particleName, mass, q, density, temperature) result(self)
        ! Construct particle object, sizeIncrease is fraction larger stored array compared to initial amount of particles
        ! In future, use hash function for possible k = 1 .. Nx, m amount of boundaries, p = prime number  m < p < N_x. h(k) = (k%p)%m
        real(real64), intent(in) :: mass, q, density, temperature
        character(*), intent(in) :: particleName
        self % name = particleName
        self % mass = mass
        self % q = q
        self % density = density ! #N/m^3
        self % temperature = temperature !temperature
    end function targetParticle_constructor


    ! ------------------------ Generating and initializing velocity with Maxwellian --------------------------


    function generate3DMaxwellianVelocity(self, irand) result(res)
        ! random velocity generator for the particle for temperature T (eV)
        ! Use box-muller method for random guassian variable, same as gwenael but doesn't have factor 2? Maybe factored into v_th
        class(Particle), intent(in out) :: self
        integer(int32), intent(in out) :: irand
        real(real64) :: U1, U2, U3, U4, res(3), v_therm
        U1 = ran2(irand)
        U2 = ran2(irand)
        U3 = ran2(irand)
        U4 = ran2(irand)
        v_term = SQRT(self.temperature*k_B/ self%mass)
        res(1) = v_therm * SQRT(-2 * LOG(U1)) * COS(2 * pi * U2)
        res(2) = v_therm * SQRT(-2 * LOG(U1)) * SIN(2 * pi * U2)
        res(3) = v_therm * SQRT(-2 * LOG(U3)) * SIN(2 * pi * U4)
    end function generate3DMaxwellianVelocity


end module mod_targetParticle