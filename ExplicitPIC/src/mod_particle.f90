module mod_particle

    use iso_fortran_env, only: int32, real64, output_unit
    use constants
    use mod_BasicFunctions
    implicit none

    private
    public :: Particle

    ! Particle contains particle properties and stored values in phase space df
    type :: Particle
        character(:), allocatable :: name !name of the particle
        integer(int32) :: N_p, finalIdx !N_p is the current last index of particle, final idx is the last allowable in memory. Index starts at 1
        real(real64), allocatable :: phaseSpace(:,:) !particle phase space, represents [l_x, v_x, v_y, v_z] in first index
        real(real64) :: mass, q, w_p ! mass (kg), charge(C), and weight (N/m^2 in 1D) of particles. Assume constant weight for moment
        real(real64) :: delIdx, wallLoss(2), energyLoss(2), refIdx !keep track particle losses at boundaries
        real(real64), allocatable :: refRecordIdx(:) 

    contains
        procedure, public, pass(self) :: initialize_n_ave
        procedure, public, pass(self) :: generate3DMaxwellian
        procedure, public, pass(self) :: getKEAve
        procedure, public, pass(self) :: getTotalKE
        procedure, public, pass(self) :: getTotalKE1D
        procedure, public, pass(self) :: getTotalMomentum
        procedure, public, pass(self) :: getVrms
        procedure, public, pass(self) :: writePhaseSpace
        procedure, public, pass(self) :: writeLocalTemperature
    end type Particle


    interface Particle
        module procedure :: particle_constructor
    end interface Particle

contains

    type(Particle) function particle_constructor(mass, q, w_p, N_p, finalIdx, particleName) result(self)
        ! Construct particle object, sizeIncrease is fraction larger stored array compared to initial amount of particles
        ! In future, use hash function for possible k = 1 .. Nx, m amount of boundaries, p = prime number  m < p < N_x. h(k) = (k%p)%m
        real(real64), intent(in) :: mass, q, w_p
        integer(int32), intent(in) :: N_p, finalIdx
        character(*), intent(in) :: particleName
        self % name = particleName
        self % mass = mass
        self % q = q
        self % w_p = w_p
        self % N_p = N_p
        self % finalIdx = finalIdx
        allocate(self%phaseSpace(4,finalIdx), self%refRecordIdx(INT(N_p/10)))
    end function particle_constructor

    pure subroutine initialize_n_ave(self, n_ave, L_domain)
        ! initialize w_p based on initial average n (N/m^3) over domain
        class(Particle), intent(in out) :: self
        real(real64), intent(in) :: n_ave, L_domain
        self % w_p = n_ave * L_domain / self % N_p
    end subroutine initialize_n_ave

    ! ------------------------ Generating and initializing velocity with Maxwellian --------------------------


    subroutine generate3DMaxwellian(self, T, irand)
        ! random velocity generator for the particle for temperature T (eV)
        ! Use box-muller method for random guassian variable, same as gwenael but doesn't have factor 2? Maybe factored into v_th
        class(Particle), intent(in out) :: self
        real(real64), intent(in) :: T
        integer(int32), intent(in out) :: irand
        real(real64) :: U1(self%N_p), U2(self%N_p), U3(self%N_p), U4(self%N_p)
        call getRandom(U1, irand)
        call getRandom(U2, irand)
        call getRandom(U3, irand)
        call getRandom(U4, irand)
        self%phaseSpace(2, 1:self%N_p) = SQRT(T*e/ self%mass) * SQRT(-2 * LOG(U1)) * COS(2 * pi * U2)
        self%phaseSpace(3, 1:self%N_p) = SQRT(T*e/ self%mass) * SQRT(-2 * LOG(U1)) * SIN(2 * pi * U2)
        self%phaseSpace(4, 1:self%N_p) = SQRT(T*e/ self%mass) * SQRT(-2 * LOG(U3)) * SIN(2 * pi * U4)
    end subroutine generate3DMaxwellian

    pure function getKEAve(self) result(res)
        ! calculate average kinetic energy (temperature) in eV
        class(Particle), intent(in) :: self
        real(real64) :: res
        res = SUM(self%phaseSpace(2:4, 1:self%N_p)**2) * self % mass * 0.5d0 / e / self%N_p

    end function getKEAve

    pure function getTotalMomentum(self) result(res)
        class(Particle), intent(in) :: self
        real(real64) :: res(3)
        res = self%w_p * self%mass * SUM(self%phaseSpace(2:4, 1:self%N_p), DIM = 2)
    end function getTotalMomentum

    pure function getVrms(self) result(res)
        ! get rms velocity for checking
        class(Particle), intent(in) :: self
        real(real64) :: res
        res = SQRT(SUM(self%phaseSpace(2:4, 1:self%N_p)**2)/ self%N_p)
    end function getVrms

    pure function getTotalKE(self) result(res)
        ! calculate total KE in Joules/m^2
        class(Particle), intent(in) :: self
        real(real64) :: res
        res = SUM(self%phaseSpace(2:4, 1:self%N_p)**2) * self % mass * 0.5d0 * self%w_p

    end function getTotalKE

    pure function getTotalKE1D(self) result(res)
        ! calculate total KE in Joules/m^2 for domain dimension
        class(Particle), intent(in) :: self
        real(real64) :: res
        res = SUM(self%phaseSpace(2, 1:self%N_p)**2) * self % mass * 0.5d0 * self%w_p

    end function getTotalKE1D

    ! ------------------------ delete/add particles --------------------------------------------
    

    ! --------------------------- Writing Particle Data to File -----------------------------------

    subroutine writePhaseSpace(self, CurrentDiagStep)
        ! Writes particle phase space into binary file
        class(Particle), intent(in) :: self
        integer(int32), intent(in) :: CurrentDiagStep
        character(len=5) :: char_i
        write(char_i, '(I3)'), CurrentDiagStep
        open(10,file='../Data/PhaseSpace/phaseSpace_'//self%name//"_"//trim(adjustl(char_i))//".dat", form='UNFORMATTED')
        write(10) self%phaseSpace(:, 1:self%N_p)
        close(10)
    end subroutine writePhaseSpace

    subroutine writeLocalTemperature(self, CurrentDiagStep)
        ! Write particle temperature averaged over local grid
        class(Particle), intent(in) :: self
        integer(int32), intent(in) :: CurrentDiagStep
        character(len=5) :: char_i
        integer(int32) :: j, index, counter(NumberXNodes-1)
        real(real64) :: temp(NumberXNodes-1)
        temp = 0.0d0
        counter = 0
    
        do j = 1, self%N_p
            index = INT(self%phaseSpace(1, j))
            temp(index) = temp(index) + SUM(self%phaseSpace(2:4, j)**2) * 0.5d0 * self%mass/e
            counter(index) = counter(index) + 1
        end do
        do j = 1, NumberXNodes-1
            if (counter(j) > 0) then
                temp(j) = temp(j)*2.0d0/counter(j)/3.0d0
            end if
        end do
        write(char_i, '(I3)'), CurrentDiagStep
        open(10,file='../Data/ElectronTemperature/eTemp_'//trim(adjustl(char_i))//".dat", form='UNFORMATTED')
        write(10) temp
        close(10)
        
    end subroutine writeLocalTemperature

end module mod_particle