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
        integer(int32), allocatable :: boundaryConditions(:) ! Boundary condition flags for fields and particles
        ! (>0 dirichlet, -2 Neumann, -3 periodic, <=-4 dielectric), 0 is default in-body condition 

    contains
        procedure, private, pass(self) :: derive_DxDl_NodeVol
        procedure, public, pass(self) :: constructSineGrid
        procedure, public, pass(self) :: constructUniformGrid
        procedure, public, pass(self) :: constructGrid
        procedure, public, pass(self) :: writeDomain
    end type Domain


    interface Domain
        module procedure :: domain_constructor
    end interface Domain

contains

    ! Initialization procedures
    type(Domain) function domain_constructor(leftBoundary, rightBoundary) result(self)
        ! Construct domain object, initialize grid, dx_dl, and nodeVol.
        integer(int32), intent(in) :: leftBoundary, rightBoundary
        integer(int32) :: i
        allocate(self % grid(NumberXNodes), self % dx_dl(NumberXNodes-1), self % nodeVol(NumberXNodes), self%boundaryConditions(NumberXNodes))
        self % grid = (/(i, i=1, NumberXNodes)/)
        self % dx_dl = 1.0d0
        self % nodeVol = 1.0d0
        self % boundaryConditions = 0
        self%boundaryConditions(1) = leftBoundary
        self%boundaryConditions(NumberXNodes) = rightBoundary
    end function domain_constructor

    pure subroutine derive_DxDl_NodeVol(self)
        class(Domain), intent(in out) :: self
        integer(int32) :: i
        do i = 1, NumberXNodes-1
            self%dx_dl(i) = self%grid(i+1) - self%grid(i)
        end do

        if (self%boundaryConditions(1) == 1) then
            self%nodeVol(1) = self%dx_dl(1)
        else if (self%boundaryConditions(1) == 3) then
            self%nodeVol(1) = (self%dx_dl(1) + self%dx_dl(NumberXNodes-1))/2
        end if

        if (self%boundaryConditions(NumberXNodes) == 1) then
            self%nodeVol(NumberXNodes) = self%dx_dl(NumberXNodes-1)
        else if (self%boundaryConditions(NumberXNodes) == 3) then
            self%nodeVol(NumberXNodes) = (self%dx_dl(1) + self%dx_dl(NumberXNodes-1))/2
        end if

        do i = 2, NumberXNodes-1
            self%nodeVol(i) = (self%dx_dl(i-1) + self%dx_dl(i))/2
        end do
    end subroutine derive_DxDl_NodeVol

    subroutine constructGrid(self, del_l, L_domain, gridType)
        class(Domain), intent(in out) :: self
        real(real64), intent(in) :: del_l, L_domain
        integer(int32), intent(in) :: gridType
        if (gridType == 0) then
            call self%constructUniformGrid(L_domain)
        else
            call self%constructSineGrid(del_l, L_domain)
        end if
    end subroutine constructGrid

    subroutine constructSineGrid(self, del_l, L_domain)
        class(Domain), intent(in out) :: self
        real(real64), intent(in) :: del_l, L_domain
        integer(int32) :: i
        self%grid(1) = 0.0d0
        self%grid(NumberXNodes) = L_domain
        do concurrent (i = 2:NumberXNodes-1)
            self % grid(i) = L_domain * ((i-1) - (NumberXNodes - 1)*(1 - del_l)/pi/2 &
            * SIN(2 * pi * (i-1) / (NumberXNodes - 1))) / (NumberXNodes - 1)
        end do
        call self%derive_DxDl_NodeVol()
    end subroutine constructSineGrid

    subroutine constructUniformGrid(self, L_domain)
        class(Domain), intent(in out) :: self
        real(real64), intent(in) :: L_domain
        integer(int32) :: i
        self%grid(1) = 0.0d0
        self%grid(NumberXNodes) = L_domain
        do i = 2, NumberXNodes-1
            self % grid(i) =  (i-1) * L_domain / (NumberXNodes - 1)
        end do
        call self%derive_DxDl_NodeVol()
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