module mod_particle

    use iso_fortran_env, only: int32, real64, output_unit
    use constants
    implicit none

    private
    public :: Particle

    ! Particle contains particle properties and stored values in phase space
    type :: Particle
        integer(int32) :: N_p, finalIdx !N_p is the current last index of particle, final idx is the last allowable in memory. Index starts at 1
        real(real64), allocatable :: l_p(:) !positions of particles in logical space
        real(real64), allocatable :: v_p(:, :) !velocities of particles in m/s
        real(real64) :: mass, q, w_p ! mass (kg), charge(C), and weight (N/m^2 in 1D) of particles. Assume constant weight for moment

    contains
        procedure, public, pass(self) :: initialize_n_ave
        procedure, public, pass(self) :: initialize_randUniform
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
            numInCell = NINT(self%N_p * dx_dl(i)/L_domain)
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
        
        print *, numPerCell
        print *, "sum of numPerCell is:", sum(numPerCell)
        
    end subroutine initialize_randUniform

end module mod_particle