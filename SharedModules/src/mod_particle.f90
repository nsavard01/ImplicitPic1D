module mod_particle

    use iso_fortran_env, only: int32, real64, output_unit
    use constants
    use mod_BasicFunctions
    use mod_domain
    implicit none

    private
    public :: Particle

    ! Particle contains particle properties and stored values in phase space df
    type :: Particle
        character(:), allocatable :: name !name of the particle
        integer(int32) :: finalIdx, partPerCell, delIdx, refIdx !N_p is the current last index of particle, final idx is the last allowable in memory. Index starts at 1
        real(real64), allocatable :: phaseSpace(:,:,:) !particle phase space, represents [l_x, v_x, v_y, v_z] in first index, particle number in second index, cell number in final index
        real(real64), allocatable :: w_p(:,:)
        integer(int32), allocatable :: N_p(:)
        real(real64) :: mass, q ! mass (kg), charge(C), and weight (N/m^2 in 1D) of particles. Assume constant weight for moment
        real(real64) :: wallLoss(2), energyLoss(2) !keep track particle losses at boundaries
        real(real64), allocatable :: refPhaseSpace(:, :), refw_p(:) !saved particle location and weight 

    contains
        procedure, public, pass(self) :: initialize_randUniform
        procedure, public, pass(self) :: generate3DMaxwellian
        procedure, public, pass(self) :: getKEAve
        procedure, public, pass(self) :: getLocalTemperature
        procedure, public, pass(self) :: getTotalKE
        procedure, public, pass(self) :: getTotalKE1D
        procedure, public, pass(self) :: getSumWeights
        procedure, public, pass(self) :: getTotalMomentum
        procedure, public, pass(self) :: getVrms
        procedure, public, pass(self) :: writePhaseSpace
        procedure, public, pass(self) :: loadParticleDensity
        procedure, public, pass(self) :: depositRho
        procedure, public, pass(self) :: writeLocalTemperature
    end type Particle


    interface Particle
        module procedure :: particle_constructor
    end interface Particle

contains

    type(Particle) function particle_constructor(mass, q, partPerCell, partIdxFactor, particleName) result(self)
        ! Construct particle object, sizeIncrease is fraction larger stored array compared to initial amount of particles
        ! In future, use hash function for possible k = 1 .. Nx, m amount of boundaries, p = prime number  m < p < N_x. h(k) = (k%p)%m
        real(real64), intent(in) :: mass, q
        integer(int32), intent(in) :: partPerCell, partIdxFactor
        character(*), intent(in) :: particleName
        self % name = particleName
        self % mass = mass
        self % q = q
        self % partPerCell = partPerCell
        self % finalIdx = partPerCell * partIdxFactor
        self%energyLoss = 0.0d0
        self%wallLoss = 0.0d0
        self%delIdx = 0
        self%refIdx = 0
        allocate(self%phaseSpace(4,self%finalIdx, NumberXNodes-1), self%refPhaseSpace(4, partPerCell), self%refw_p(partPerCell),&
            self%w_p(self%finalIdx, NumberXNodes-1), self%N_p(NumberXNodes-1))
        self%N_p = self%partPerCell
    end function particle_constructor

    subroutine initialize_randUniform(self, world, irand)
        ! place particles randomly in each dx_dl based on portion of volume it take up
        class(Particle), intent(in out) :: self
        type(Domain), intent(in) :: world
        integer(int32), intent(in out) :: irand
        integer(int32) :: i, j
        real(real64) :: w_0, L_domain
        L_domain = world%grid(NumberXNodes) - world%grid(1)
        w_0 = n_ave / real(self%partPerCell, kind = 8)
        do j = 1, NumberXNodes-1
            do i = 1, self%partPerCell
                self%phaseSpace(1, i, j) = j + ran2(irand)
            end do
            self%w_p(1:self%partPerCell, j) = w_0 * world%dx_dl(j)
        end do
        self%N_p = self%partPerCell
        
    end subroutine initialize_randUniform


    ! ------------------------ Generating and initializing velocity with Maxwellian --------------------------


    subroutine generate3DMaxwellian(self, T, irand)
        ! random velocity generator for the particle for temperature T (eV)
        ! Use box-muller method for random guassian variable, same as gwenael but doesn't have factor 2? Maybe factored into v_th
        class(Particle), intent(in out) :: self
        real(real64), intent(in) :: T
        integer(int32), intent(in out) :: irand
        real(real64) :: U(4)
        integer(int32) :: i, j
        do i = 1, NumberXNodes-1
            do j = 1, self%N_p(i)
                call getRandom(U, irand)
                self%phaseSpace(2, j, i) = SQRT(T*e/ self%mass) * SQRT(-2 * LOG(U(1))) * COS(2 * pi * U(2))
                self%phaseSpace(3, j, i) = SQRT(T*e/ self%mass) * SQRT(-2 * LOG(U(1))) * SIN(2 * pi * U(2))
                self%phaseSpace(4, j, i) = SQRT(T*e/ self%mass) * SQRT(-2 * LOG(U(3))) * SIN(2 * pi * U(4))
            end do
        end do
    end subroutine generate3DMaxwellian

    pure function getKEAve(self) result(res)
        ! calculate average kinetic energy (temperature) in eV
        class(Particle), intent(in) :: self
        real(real64) :: res
        integer(int32) :: i
        res = 0.0d0
        do i = 1, NumberXNodes-1
            res = res + SUM(self%phaseSpace(2:4, 1:self%N_p(i), i)**2)
        end do
        res = res * self%mass * 0.5d0 / e / SUM(self%N_p)

    end function getKEAve

    function getSumWeights(self) result(res)
        ! calculate average kinetic energy (temperature) in eV
        class(Particle), intent(in) :: self
        real(real64) :: res
        integer(int32) :: i
        res = 0.0d0
        do i = 1, NumberXNodes-1
            res = res + SUM(self%w_p(1:self%N_p(i), i))
        end do
        res = res

    end function getSumWeights

    pure function getTotalMomentum(self) result(res)
        class(Particle), intent(in) :: self
        real(real64) :: res(3)
        integer(int32) :: i, j
        res = 0.0d0
        do i = 1, NumberXNodes-1
            do j = 1, self%N_p(i)
                res = res + self%w_p(j, i) * self%phaseSpace(2:4, j, i)
            end do
        end do
        res = res * self%mass
    end function getTotalMomentum

    pure function getVrms(self) result(res)
        ! get rms velocity for checking
        class(Particle), intent(in) :: self
        real(real64) :: res
        integer(int32) :: i
        res = 0.0d0
        do i = 1, NumberXNodes-1
            res = res + SUM(self%phaseSpace(2:4, 1:self%N_p(i), i)**2)
        end do
        res = SQRT(res / SUM(self%N_p))
    end function getVrms

    pure function getTotalKE(self) result(res)
        ! calculate total KE in Joules/m^2
        class(Particle), intent(in) :: self
        real(real64) :: res
        integer(int32) :: i, j
        res = 0.0d0
        do i = 1, NumberXNodes-1
            do j = 1, self%N_p(i)
                res = res + SUM(self%phaseSpace(2:4, j, i)**2) * self%w_p(j, i)
            end do
        end do
        do i = 1, self%refIdx
            res = res + SUM(self%refPhaseSpace(2:4, i)**2) * self%refw_p(i)
        end do
        res = res * self%mass * 0.5d0
    end function getTotalKE

    pure function getTotalKE1D(self) result(res)
        ! calculate total KE in Joules/m^2 for domain dimension
        class(Particle), intent(in) :: self
        real(real64) :: res
        integer(int32) :: i, j
        res = 0.0d0
        do i = 1, NumberXNodes-1
            do j = 1, self%N_p(i)
                res = res + self%phaseSpace(2, j, i)**2 * self%w_p(j, i)
            end do
        end do
        res = res * self%mass * 0.5d0

    end function getTotalKE1D

    function getLocalTemperature(self) result(res)
        class(Particle), intent(in) :: self
        real(real64) :: res(NumberXNodes-1)
        integer(int32) :: i
        do i = 1, NumberXNodes-1
            res(i) = SUM(self%phaseSpace(2:4, 1:self%N_p(i), i)**2) * self % mass * 0.5d0 / e / self%N_p(i)
        end do
    end function getLocalTemperature

    ! ------------------------ Gather to mesh  --------------------------------------------
    subroutine loadParticleDensity(self, densities, world)
        class(Particle), intent(in) :: self
        type(Domain), intent(in) :: world
        real(real64), intent(in out) :: densities(NumberXNodes)
        integer(int32) :: i, j
        real(real64) :: d, w_p
        do i = 1, NumberXNodes-1
            do j = 1, self%N_p(i)
                w_p = self%w_p(j, i)
                d = MOD(self%phaseSpace(1,j, i), 1.0d0)
                densities(i) = densities(i) + w_p * (1.0d0-d)
                densities(i+1) = densities(i + 1) + w_p * d
            end do
        end do
        if (world%boundaryConditions(1) > 4) then
            print *, 'Just excuse to use world object since needed in CIC scheme'
        end if

    end subroutine loadParticleDensity

    subroutine depositRho(self, rho, world) 
        class(Particle), intent(in) :: self
        real(real64), intent(in out) :: rho(NumberXNodes)
        type(Domain), intent(in) :: world
        real(real64) :: densities(NumberXNodes)
        densities = 0.0d0
        call self%loadParticleDensity(densities, world)
        if (world%boundaryConditions(1) == 3) then
            rho(1) = rho(1) + rho(NumberXNodes)
            rho(NumberXNodes) = rho(1)
        end if
        rho = rho + densities * self%q
    end subroutine depositRho

    ! --------------------------- Writing Particle Data to File -----------------------------------

    subroutine writePhaseSpace(self, CurrentDiagStep, dirName)
        ! Writes particle phase space into binary file
        class(Particle), intent(in) :: self
        integer(int32), intent(in) :: CurrentDiagStep
        character(*), intent(in) :: dirName
        character(len=5) :: char_i
        integer(int32) :: i
        write(char_i, '(I3)'), CurrentDiagStep
        open(10,file='../'//dirName//'/PhaseSpace/phaseSpace_'//self%name//"_"//trim(adjustl(char_i))//".dat", form='UNFORMATTED')
        do i = 1, NumberXNodes-1
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
        real(real64) :: temp(NumberXNodes-1)
        temp = self%getLocalTemperature()
        write(char_i, '(I3)'), CurrentDiagStep
        open(10,file='../'//dirName//'/ElectronTemperature/eTemp_'//trim(adjustl(char_i))//".dat", form='UNFORMATTED')
        write(10) temp
        close(10)
        
    end subroutine writeLocalTemperature

    subroutine writeParticleDensity(self, densities, world, CurrentDiagStep, boolAverage, dirName) 
        ! For diagnostics, deposit single particle density
        ! Re-use rho array since it doesn't get used after first Poisson
        real(real64), intent(in out) :: densities(NumberXNodes)
        class(Particle), intent(in) :: self
        type(Domain), intent(in) :: world
        integer(int32), intent(in) :: CurrentDiagStep
        character(*), intent(in) :: dirName
        logical, intent(in) :: boolAverage
        character(len=5) :: char_i
        
        densities = densities/world%nodeVol
        if (world%boundaryConditions(1) == 3) then
            densities(1) = densities(1) + densities(NumberXNodes)
            densities(NumberXNodes) = densities(1)
        end if

        write(char_i, '(I3)'), CurrentDiagStep
        if (boolAverage) then
            open(41,file='../'//dirName//'/Density/density_'//self%name//"_Average.dat", form='UNFORMATTED')
        else
            open(41,file='../'//dirName//'/Density/density_'//self%name//"_"//trim(adjustl(char_i))//".dat", form='UNFORMATTED')
        end if
        write(41) densities
        close(41)
        
    end subroutine writeParticleDensity

    

end module mod_particle