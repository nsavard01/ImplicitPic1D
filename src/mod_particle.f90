module mod_particle

    use iso_fortran_env, only: int32, real64, output_unit
    use constants
    implicit none

    private
    public :: Particle

    ! Particle contains particle properties and stored values in phase space
    type :: Particle
        integer(int32) :: N_p, finalIdx !N_p is the current last index of particle, final idx is the last allowable in memory. Index starts at 1
        real(real64), allocatable :: x(:,:) !positions of particles in m
        real(real64), allocatable :: v(:, :) !velocities of particles in m/s
        real(real64) :: mass, q, w_p ! mass (kg), charge(C), and weight (N/m^2 in 1D) of particles. Assume constant weight for moment

    ! contains
    !     procedure, private, pass(self) :: derive_DxDl_NodeVol
    !     procedure, public, pass(self) :: constructSineGrid
    !     procedure, public, pass(self) :: constructUniformGrid
    end type Particle


    interface Particle
        module procedure :: particle_constructor
    end interface Particle

contains

    type(Particle) function particle_constructor(mass, q, w_p, N_p, finalIdx) result(self)
        ! Construct particle object, sizeIncrease is fraction larger stored array compared to initial amount of particles
        real(real64), intent(in) :: mass, q, w_p
        integer(int32), intent(in) :: N_p, finalIdx
        self % mass = mass
        self % q = q
        self % w_p = w_p
        self % N_p = N_p
        self % finalIdx = finalIdx
        allocate(self%x(self%finalIdx, 3), self%v(self%finalIdx, 3))
        self%x = 0
        self%v = 0

    end function particle_constructor

end module mod_particle