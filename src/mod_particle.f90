module mod_particle

    use iso_fortran_env, only: int32, real64, output_unit
    use constants
    implicit none

    private
    public :: Particle

    ! Particle contains particle properties and stored values in phase space
    type :: Particle
        character(:), allocatable :: name !name of the particle
        integer(int32) :: N_p, finalIdx !N_p is the current last index of particle, final idx is the last allowable in memory. Index starts at 1
        real(real64), allocatable :: l_p(:) !positions of particles in logical space
        real(real64), allocatable :: v_p(:, :) !velocities of particles in m/s
        real(real64) :: mass, q, w_p ! mass (kg), charge(C), and weight (N/m^2 in 1D) of particles. Assume constant weight for moment

    contains
        procedure, public, pass(self) :: initialize_n_ave
        procedure, public, pass(self) :: initialize_randUniform
        procedure, public, pass(self) :: generate3DMaxwellian
        procedure, public, pass(self) :: getTemperature
        procedure, public, pass(self) :: getVrms
    end type Particle


    interface Particle
        module procedure :: particle_constructor
    end interface Particle

contains

    type(Particle) function particle_constructor(mass, q, w_p, N_p, finalIdx, particleName) result(self)
        ! Construct particle object, sizeIncrease is fraction larger stored array compared to initial amount of particles
        real(real64), intent(in) :: mass, q, w_p
        integer(int32), intent(in) :: N_p, finalIdx
        character(*), intent(in) :: particleName
        self % name = particleName
        self % mass = mass
        self % q = q
        self % w_p = w_p
        self % N_p = N_p
        self % finalIdx = finalIdx
        allocate(self%l_p(self%finalIdx), self%v_p(self%finalIdx, 3))
        self%l_p = 0
        self%v_p = 0

    end function particle_constructor

    pure subroutine initialize_n_ave(self, n_ave, L_domain)
        ! initialize w_p based on initial average n (N/m^3) over domain
        class(Particle), intent(in out) :: self
        real(real64), intent(in) :: n_ave, L_domain
        self % w_p = n_ave * L_domain / self % N_p
    end subroutine initialize_n_ave

    subroutine initialize_randUniform(self, L_domain, dx_dl, n_x)
        ! place particles randomly in each dx_dl based on portion of volume it take up
        class(Particle), intent(in out) :: self
        real(real64), intent(in) :: dx_dl(:), L_domain
        integer(int32) :: i, numInCell, idxLower, numPerCell(n_x-1)
        real(real64) :: sumDxDl
        integer(int32), intent(in):: n_x
        idxLower = 1
        sumDxDl = 0
        do i=1, n_x-1
            ! Use int to make sure always have a bit left over, other will fill up before getting to end
            numInCell = INT(self%N_p * dx_dl(i)/L_domain)
            if (idxLower + numInCell > self % N_P + 1) then
                stop "You are putting too many particles for the uniform particle case"
            end if

            
            call random_number(self%l_p(idxLower:idxLower + numInCell-1))
            self%l_p(idxLower:idxLower + numInCell - 1) = self%l_p(idxLower:idxLower + numInCell - 1) + i
            idxLower = idxLower + numInCell
            numPerCell(i) = numInCell
            
        end do
        if (idxLower < self%N_p + 1) then
            call random_number(self%l_p(idxLower:self%N_p))
            self%l_p(idxLower:self%N_p) = self%l_p(idxLower:self%N_p) * (n_x - 1) + 1
        end if
        
    end subroutine initialize_randUniform

    subroutine generate3DMaxwellian(self, T)
        ! random velocity generator for the particle for temperature T (eV)
        ! Use box-muller method for random guassian variable, same as gwenael but doesn't have factor 2? Maybe factored into v_th
        class(Particle), intent(in out) :: self
        real(real64), intent(in) :: T
        real(real64) :: U1(self%N_p), U2(self%N_p), U3(self%N_p), U4(self%N_p)
        call random_number(U1)
        call random_number(U2)
        call random_number(U3)
        call random_number(U4)
        self%v_p(1:self%N_p, 1) = SQRT(T*e/ self%mass) * SQRT(-2 * LOG(U1)) * COS(2 * pi * U2)
        self%v_p(1:self%N_p, 2) = SQRT(T*e/ self%mass) * SQRT(-2 * LOG(U1)) * SIN(2 * pi * U2)
        self%v_p(1:self%N_p, 3) = SQRT(T*e/ self%mass) * SQRT(-2 * LOG(U3)) * SIN(2 * pi * U4)
    end subroutine generate3DMaxwellian

    pure function getTemperature(self) result(res)
        ! calculate average kinetic energy (temperature) in eV
        class(Particle), intent(in) :: self
        real(real64) :: res
        res = SUM(self%v_p(1:self%N_p, :)**2) * self % mass * 0.5 / e / self%N_p

    end function getTemperature

    pure function getVrms(self) result(res)
        ! get rms velocity for checking
        class(Particle), intent(in) :: self
        real(real64) :: res
        res = SQRT(SUM(self%v_p(1:self%N_p, :)**2)/ self%N_p)
    end function getVrms

end module mod_particle