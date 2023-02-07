module mod_domain
    use iso_fortran_env, only: int32, real64
    use constants
    use mod_particle
    implicit none

    private
    public :: Domain

    ! domain contains arrays and values related to physical, logical dimensions of the spatial grid
    type :: Domain
        real(real64), allocatable :: grid(:) !array representing values at grid in m
        real(real64), allocatable :: dx_dl(:) !ratio of grid differences from physical to logical, assume logical separated by 1
        real(real64), allocatable :: nodeVol(:) !node vol, or difference between half grid in physical space
        !real(real64), allocatable :: phi(:), rho(:), J(:) ! array of phi, rho, and J values in space

    contains
        procedure, private, pass(self) :: derive_DxDl_NodeVol
        procedure, public, pass(self) :: constructSineGrid
        procedure, public, pass(self) :: constructUniformGrid
        ! procedure, public, pass(self) :: depositRho
        ! procedure, public, pass(self) :: writeRho
        procedure, public, pass(self) :: writeDomain
    end type Domain


    interface Domain
        module procedure :: domain_constructor
    end interface Domain

contains

    ! Initialization procedures
    type(Domain) function domain_constructor() result(self)
        ! Construct domain object, initialize grid, dx_dl, and nodeVol.
        integer(int32) :: i
        allocate(self % grid(NumberXNodes), self % dx_dl(NumberXNodes-1), self % nodeVol(NumberXNodes))
        self % grid = (/(i, i=1, NumberXNodes)/)
        self % dx_dl = 1
        self % nodeVol = 1
    end function domain_constructor

    pure subroutine derive_DxDl_NodeVol(self)
        class(Domain), intent(in out) :: self
        integer(int32) :: i
        do i = 1, NumberXNodes-1
            self%dx_dl(i) = self%grid(i+1) - self%grid(i)
        end do

        self%nodeVol(1) = self%dx_dl(1)/2
        self%nodeVol(NumberXNodes) = self%dx_dl(NumberXNodes-1)/2

        do i = 2, NumberXNodes-1
            self%nodeVol(i) = (self%dx_dl(i-1) + self%dx_dl(i))/2
        end do
    end subroutine derive_DxDl_NodeVol

    pure subroutine constructSineGrid(self, del_l, L_domain)
        class(Domain), intent(in out) :: self
        real(real64), intent(in) :: del_l, L_domain
        integer(int32) :: i
        self%grid(1) = 0
        self%grid(NumberXNodes) = L_domain
        do concurrent (i = 2:NumberXNodes-1)
            self % grid(i) = L_domain * ((i-1) - (NumberXNodes - 1)*(1 - del_l)/pi/2 &
            * SIN(2 * pi * (i-1) / (NumberXNodes - 1))) / (NumberXNodes - 1)
        end do
        call derive_DxDl_NodeVol(self)
    end subroutine constructSineGrid

    pure subroutine constructUniformGrid(self, L_domain)
        class(Domain), intent(in out) :: self
        real(real64), intent(in) :: L_domain
        integer(int32) :: i
        self%grid(1) = 0d0
        self%grid(NumberXNodes) = L_domain
        do i = 2, NumberXNodes-1
            self % grid(i) =  (i-1) * L_domain / (NumberXNodes - 1)
        end do
        call derive_DxDl_NodeVol(self)
    end subroutine constructUniformGrid

    ! subroutine depositRho(self, particleList) 
    !     class(Domain), intent(in out) :: self
    !     type(Particle), intent(in) :: particleList(:)
    !     integer(int32) :: i, j, l_left
    !     real(real64) :: d
    !     do i=1, size(particleList)
    !         do j = 1, particleList(i)%N_p
    !             l_left = INT(particleList(i)%l_p(j))
    !             d = MOD(particleList(i)%l_p(j), 1.0)
    !             self % rho(l_left) = self % rho(l_left) + particleList(i)%q * particleList(i)%w_p * (1.0-d)
    !             self % rho(l_left + 1) = self % rho(l_left + 1) + particleList(i)%q * particleList(i)%w_p * d
    !         end do
    !     end do
    !     self % rho = self % rho / self % nodeVol
    ! end subroutine depositRho

    ! ! Write data from rho

    ! subroutine writeRho(self)
    !     ! Writes rho into a binary file.
    !     class(Domain), intent(in out) :: self
    !     integer(int32) :: fileunit, record_length
    !     character(100) :: filename
    !     filename = 'record_Rho.dat'
    !     record_length = 2 * size(self%rho)
    !     open(newunit=fileunit, file=filename, access='direct', recl= record_length)
    !     write(unit=fileunit, rec=1) self%rho
    !     close(fileunit)
    !   end subroutine writeRho

      subroutine writeDomain(self)
        ! Writes domain data into binary file under Data
        class(Domain), intent(in) :: self
        open(41,file="../Data/domainGrid.dat", form='UNFORMATTED')
        write(41) self%grid
        close(41)
      end subroutine writeDomain





end module mod_domain