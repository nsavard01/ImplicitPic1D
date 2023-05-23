module mod_domain
    use iso_fortran_env, only: int32, real64
    use constants
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
        procedure, public, pass(self) :: constructHalfSineGrid
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


    subroutine constructGrid(self, del_x, L_domain, gridType)
        class(Domain), intent(in out) :: self
        real(real64), intent(in) :: del_x, L_domain
        integer(int32), intent(in) :: gridType
        SELECT CASE (gridType)
        CASE(0)
            call self%constructUniformGrid(L_domain)
        CASE(1)
            call self%constructSineGrid(del_x, L_domain)
        CASE(2)
            call self%constructHalfSineGrid(del_x, L_domain)
        CASE default
            print *, "Gridtype", gridType, "doesn't exist!"
            stop
        END SELECT
    end subroutine constructGrid

    subroutine constructSineGrid(self, del_x, L_domain)
        class(Domain), intent(in out) :: self
        real(real64), intent(in) :: del_x, L_domain
        integer(int32) :: i
        real(real64) :: gridField(NumberXNodes-1)
        if (del_x/L_domain >= 1.0d0/(real(NumberXNodes) - 1.0d0)) then
            print *, "The debyeLength is really large, less nodes needed!"
            stop
        end if
        self%grid(1) = 0.0d0
        self%grid(NumberXNodes) = L_domain
        if (schemeNum == 1) then
            self%grid(2) = self%grid(1) + del_x
            self%grid(NumberXNodes-1) = self%grid(NumberXNodes) - del_x
            do i = 1,NumberXNodes-1
                gridField(i) = (L_domain - del_x) * (real(i-1)/(real(NumberXNodes) - 2.0d0) - (1.0d0/(real(NumberXNodes) - 2.0d0) - del_x/(L_domain - del_x)) &
                * SIN(2.0d0 * pi * real(i-1)/(real(NumberXNodes)-2.0d0)) / SIN(2.0d0 * pi / (real(NumberXNodes) - 2.0d0) )) + del_x/2.0d0
            end do
            self%grid(3:NumberXNodes-2) = (gridField(2:NumberXNodes-3) + gridField(3:NumberXNodes-2))/2.0d0
            self%nodeVol(1) = del_x
            self%nodeVol(2) = del_x
            self%nodeVol(NumberXNodes) = del_x
            self%nodeVol(NumberXNodes-1) = del_x
            do i = 3,NumberXNodes-2
                self%nodeVol(i) = gridField(i) - gridField(i-1)
            end do
            do i = 1, NumberXNodes-1
                self%dx_dl(i) = (self%nodeVol(i+1) + self%nodeVol(i))/2.0d0
            end do
        else
            do i = 2,NumberXNodes-1
                self % grid(i) = L_domain * ((real(i)-1.0d0)/(real(NumberXNodes) - 1.0d0) - (1.0d0/(real(NumberXNodes) - 1.0d0) - del_x/L_domain) &
                * SIN(2 * pi * (i-1) / (NumberXNodes - 1)) / SIN(2 * pi / (NumberXNodes - 1)) )
            end do
            do i = 1, NumberXNodes-1
                self%dx_dl(i) = self%grid(i+1) - self%grid(i)
            end do
            do i = 2, NumberXNodes-1
                self%nodeVol(i) = (self%dx_dl(i-1) + self%dx_dl(i))/2.0d0
            end do
            self%nodeVol(1) = self%dx_dl(1)
            self%nodeVol(NumberXNodes) = self%dx_dl(NumberXNodes-1)
        end if

    end subroutine constructSineGrid

    subroutine constructHalfSineGrid(self, del_x, L_domain)
        class(Domain), intent(in out) :: self
        real(real64), intent(in) :: del_x, L_domain
        integer(int32) :: i
        real(real64) :: gridField(NumberXNodes-1)
        if (del_x/L_domain >= 1.0d0/(real(NumberXNodes) - 1.0d0)) then
            print *, "The debyeLength is really large, less nodes needed!"
            stop
        end if
        self%grid(1) = 0.0d0
        self%grid(NumberXNodes) = L_domain/2.0d0
        if (schemeNum == 1) then
            self%grid(2) = self%grid(1) + del_x
            do i = 1,NumberXNodes-1
                gridField(i) = (L_domain - del_x) * (real(i-1)/(2.0d0 * real(NumberXNodes) - 3.0d0) - (1.0d0/(2.0d0 * real(NumberXNodes) - 3.0d0) - del_x/(L_domain - del_x)) &
                * SIN(2.0d0 * pi * real(i-1)/(2.0d0 * real(NumberXNodes) - 3.0d0)) / SIN(2.0d0 * pi / (2.0d0 * real(NumberXNodes) - 3.0d0))) + del_x/2.0d0
            end do
            self%grid(3:NumberXNodes-1) = (gridField(2:NumberXNodes-2) + gridField(3:NumberXNodes-1))/2.0d0
            self%nodeVol(1) = del_x
            self%nodeVol(NumberXNodes) = 2.0d0 * (self%grid(NumberXNodes) - gridField(NumberXNodes-1))
            do i = 2,NumberXNodes-1
                self%nodeVol(i) = gridField(i) - gridField(i-1)
            end do
            do i = 1, NumberXNodes-1
                self%dx_dl(i) = self%grid(i+1) - self%grid(i)
            end do
        else
            do i = 2,NumberXNodes-1
                self % grid(i) = L_domain * ((real(i)-1.0d0)/(real(NumberXNodes) - 1.0d0)/2.0d0 - (1.0d0/(real(NumberXNodes) - 1.0d0)/2.0d0 - del_x/L_domain) &
                * SIN(pi * (real(i)-1)/(real(NumberXNodes) - 1)) / SIN(pi / (real(NumberXNodes) - 1.0d0)) )
            end do
            do i = 1, NumberXNodes-1
                self%dx_dl(i) = self%grid(i+1) - self%grid(i)
            end do
            do i = 2, NumberXNodes-1
                self%nodeVol(i) = (self%dx_dl(i-1) + self%dx_dl(i))/2
            end do
            self%nodeVol(1) = self%dx_dl(1)
            self%nodeVol(NumberXNodes) = self%dx_dl(NumberXNodes-1)
        end if

        if (self%boundaryConditions(1) == 3 .or. self%boundaryConditions(NumberXNodes) == 3) then
            print *, "Mesh is not periodic, cannot have periodic boundary!"
            stop
        end if
    end subroutine constructHalfSineGrid

    subroutine constructUniformGrid(self, L_domain)
        class(Domain), intent(in out) :: self
        real(real64), intent(in) :: L_domain
        integer(int32) :: i
        self%grid(1) = 0.0d0
        self%grid(NumberXNodes) = L_domain
        do i = 2, NumberXNodes-1
            self % grid(i) =  (i-1) * L_domain / (NumberXNodes - 1)
        end do
        do i = 1, NumberXNodes-1
            self%dx_dl(i) = self%grid(i+1) - self%grid(i)
        end do
        do i = 2, NumberXNodes-1
            self%nodeVol(i) = (self%dx_dl(i-1) + self%dx_dl(i))/2
        end do
        self%nodeVol(1) = self%dx_dl(1)
        self%nodeVol(NumberXNodes) = self%dx_dl(NumberXNodes-1)
    end subroutine constructUniformGrid


    subroutine writeDomain(self)
    ! Writes domain data into binary file under Data
    class(Domain), intent(in) :: self
    open(41,file="../Data/domainGrid.dat", form='UNFORMATTED')
    write(41) self%grid
    close(41)
    open(41,file="../Data/domainDxDl.dat", form='UNFORMATTED')
    write(41) self%dx_dl
    close(41)
    open(41,file="../Data/domainNodeVol.dat", form='UNFORMATTED')
    write(41) self%nodeVol
    close(41)
    end subroutine writeDomain





end module mod_domain