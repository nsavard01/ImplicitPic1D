module mod_particle

    use iso_fortran_env, only: int32, real64, output_unit
    use constants
    use mod_BasicFunctions
    use omp_lib
    implicit none

    private
    public :: Particle

    ! Particle contains particle properties and stored values in phase space df
    type :: Particle
        character(:), allocatable :: name !name of the particle
        integer(int32), allocatable :: N_p(:), refIdx(:), refRecordIdx(:, :), delIdx(:) !N_p is the current last index of particle, final idx is the last allowable in memory. Index starts at 1
        integer(int32) :: finalIdx, accumWallLoss(2)
        real(real64), allocatable :: phaseSpace(:,:, :) !particle phase space, represents [l_x, v_x, v_y, v_z] in first index
        real(real64) :: mass, q, w_p ! mass (kg), charge(C), and weight (N/m^2 in 1D) of particles. Assume constant weight for moment
        real(real64), allocatable :: wallLoss(:, :), energyLoss(:, :)
        real(real64), allocatable :: densities(:, :)
        real(real64) :: numFuncEvalAve, numSubStepsAve
        real(real64) :: accumEnergyLoss(2)

    contains
        procedure, public, pass(self) :: initialize_n_ave
        procedure, public, pass(self) :: initializeRandUniform
        procedure, public, pass(self) :: generate3DMaxwellian
        procedure, public, pass(self) :: getKEAve
        procedure, public, pass(self) :: getTotalKE
        procedure, public, pass(self) :: getTotalMomentum
        procedure, public, pass(self) :: writePhaseSpace
        procedure, public, pass(self) :: writeLocalTemperature
    end type Particle


    interface Particle
        module procedure :: particle_constructor
    end interface Particle

contains

    type(Particle) function particle_constructor(mass, q, w_p, N_p, finalIdx, particleName, numThread) result(self)
        ! Construct particle object, sizeIncrease is fraction larger stored array compared to initial amount of particles
        ! In future, use hash function for possible k = 1 .. Nx, m amount of boundaries, p = prime number  m < p < N_x. h(k) = (k%p)%m
        real(real64), intent(in) :: mass, q, w_p
        integer(int32), intent(in) :: N_p, finalIdx, numThread
        character(*), intent(in) :: particleName
        self % name = particleName
        self % mass = mass
        self % q = q
        self % w_p = w_p
        self % finalIdx = finalIdx
        self%accumWallLoss = 0.0d0
        self%accumEnergyLoss = 0.0d0
        self%numFuncEvalAve = 0.0d0
        self%numSubStepsAve = 0.0d0
        allocate(self%phaseSpace(4,finalIdx, numThread), self%refRecordIdx(INT(self%finalIdx/10), numThread), self%N_p(numThread), &
            self%delIdx(numThread), self%wallLoss(2, numThread), self%energyLoss(2, numThread), self%refIdx(numThread), self%densities(NumberXNodes, numThread))
        self%refIdx = 0
        self%delIdx = 0
        self%N_p = N_p
        self%energyLoss = 0.0d0
        self%wallLoss = 0
    end function particle_constructor

    pure subroutine initialize_n_ave(self, n_ave, L_domain)
        ! initialize w_p based on initial average n (N/m^3) over domain
        class(Particle), intent(in out) :: self
        real(real64), intent(in) :: n_ave, L_domain
        self % w_p = n_ave * L_domain / SUM(self % N_p)
    end subroutine initialize_n_ave

    subroutine initializeRandUniform(self, irand)
        class(Particle), intent(in out) :: self
        integer(int32), intent(in out) :: irand(2,numThread)
        integer(int32) :: iThread, i
        !$OMP parallel private(iThread, i)
        iThread = omp_get_thread_num() + 1
        do i = 1, self%N_p(iThread)
            self%phaseSpace(1, i, iThread) = randNew(irand(:, iThread)) * real(NumberXNodes-1) + 1.0d0
        end do
        !$OMP end parallel
    end subroutine initializeRandUniform

    ! ------------------------ Generating and initializing velocity with Maxwellian --------------------------


    subroutine generate3DMaxwellian(self, T, irand)
        ! random velocity generator for the particle for temperature T (eV)
        ! Use box-muller method for random guassian variable, same as gwenael but doesn't have factor 2? Maybe factored into v_th
        class(Particle), intent(in out) :: self
        real(real64), intent(in) :: T
        integer(int32), intent(in out) :: irand(2,numThread)
        integer(int32) :: iThread, i
        integer(int32) :: irand_thread(2)
        real(real64) :: U1, U2, U3, U4, v_therm
        v_therm = SQRT(T*e/self%mass)
        !$OMP PARALLEL PRIVATE(iThread, i, U1, U2, U3, U4, irand_thread)
        iThread = omp_get_thread_num() + 1
        irand_thread = irand(:,iThread)
        do i = 1, self%N_p(iThread)
            U1 = randNew(irand_thread)
            U2 = randNew(irand_thread)
            U3 = randNew(irand_thread)
            U4 = randNew(irand_thread)
            self%phaseSpace(2, i, iThread) = v_therm * SQRT(-2.0d0 * LOG(U1)) * COS(2.0d0 * pi * U2)
            self%phaseSpace(3, i, iThread) = v_therm * SQRT(-2.0d0 * LOG(U1)) * SIN(2.0d0 * pi * U2)
            self%phaseSpace(4, i, iThread) = v_therm * SQRT(-2.0d0 * LOG(U3)) * SIN(2.0d0 * pi * U4)
        end do
        irand(:,iThread) = irand_thread
        !$OMP END PARALLEL
    end subroutine generate3DMaxwellian

    function getKEAve(self) result(res)
        ! calculate average kinetic energy (temperature) in eV
        class(Particle), intent(in) :: self
        integer(int32) :: iThread
        real(real64) :: res
        res = 0.0d0
        !$OMP parallel private(iThread) reduction(+:res)
        iThread = omp_get_thread_num() + 1
        res = res + SUM(self%phaseSpace(2:4, 1:self%N_p(iThread), iThread)**2) 
        !$OMP end parallel
        res = res * self % mass * 0.5d0 / e / SUM(self%N_p)
    end function getKEAve

    function getTotalMomentum(self) result(res)
        class(Particle), intent(in) :: self
        real(real64) :: res(3), temp(3, numThread)
        integer(int32) :: iThread
        !$OMP parallel private(iThread)
        iThread = omp_get_thread_num() + 1
        temp(:, iThread) = SUM(self%phaseSpace(2:4, 1:self%N_p(iThread), iThread), DIM = 2)
        !$OMP end parallel
        res = SUM(temp, DIM = 2) * self%w_p * self%mass
    end function getTotalMomentum

    function getTotalKE(self) result(res)
        ! calculate total KE in Joules/m^2
        class(Particle), intent(in) :: self
        real(real64) :: res
        integer(int32) :: iThread
        res = 0.0d0
        !$OMP parallel private(iThread) reduction(+:res)
        iThread = omp_get_thread_num() + 1
        res = res + SUM(self%phaseSpace(2:4, 1:self%N_p(iThread), iThread)**2) 
        !$OMP end parallel
        res = res * self % mass * 0.5d0 * self%w_p

    end function getTotalKE
    

    ! --------------------------- Writing Particle Data to File -----------------------------------

    subroutine writePhaseSpace(self, CurrentDiagStep, dirName)
        ! Writes particle phase space into binary file
        class(Particle), intent(in) :: self
        integer(int32), intent(in) :: CurrentDiagStep
        character(*), intent(in) :: dirName
        character(len=5) :: char_i
        integer(int32) :: i
        write(char_i, '(I3)'), CurrentDiagStep
        open(10,file=dirName//'/PhaseSpace/phaseSpace_'//self%name//"_"//trim(adjustl(char_i))//".dat", form='UNFORMATTED', access = 'STREAM', status = 'REPLACE')
        do i = 1, numThread
            write(10) self%phaseSpace(:, 1:self%N_p(i), i)
        end do
        close(10)
    end subroutine writePhaseSpace

    subroutine writeLocalTemperature(self, CurrentDiagStep, dirName)
        ! Write particle temperature averaged over local grid
        class(Particle), intent(in) :: self
        integer(int32), intent(in) :: CurrentDiagStep
        character(*), intent(in) :: dirName
        character(len=5) :: char_i
        integer(int32) :: j, index, counter(NumberXNodes-1, numThread), iThread
        real(real64) :: temp(NumberXNodes-1, numThread), EHist(NumberXNodes-1)
        temp = 0.0d0
        counter = 0
        !$OMP parallel private(iThread, j, index) 
        iThread = omp_get_thread_num() + 1
        do j = 1, self%N_p(iThread)
            index = INT(self%phaseSpace(1, j, iThread))
            if (index < NumberXNodes) then
                temp(index, iThread) = temp(index, iThread) + SUM(self%phaseSpace(2:4, j, iThread)**2)
                counter(index, iThread) = counter(index, iThread) + 1
            end if
        end do
        !$OMP end parallel
        do j = 1, NumberXNodes-1
            if (SUM(counter(j, :)) > 0) then
                EHist(j) = SUM(temp(j,:))*self%mass/SUM(counter(j, :))/3.0d0/e
            else
                EHist(j) = 0.0d0
            end if
        end do
        write(char_i, '(I3)'), CurrentDiagStep
        open(10,file=dirName//'/Temperature/Temp_'//self%name//"_"//trim(adjustl(char_i))//".dat", form='UNFORMATTED')
        write(10) EHist
        close(10)
        
    end subroutine writeLocalTemperature

end module mod_particle