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
        procedure, public, pass(self) :: constructExpHalfGrid
        procedure, public, pass(self) :: getLFromX
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
        CASE(3)
            call self%constructExpHalfGrid(del_x, L_domain)
        CASE default
            print *, "Gridtype", gridType, "doesn't exist!"
            stop
        END SELECT
    end subroutine constructGrid

    subroutine constructSineGrid(self, del_x, L_domain)
        class(Domain), intent(in out) :: self
        real(real64), intent(in) :: del_x, L_domain
        integer(int32) :: i
        if (del_x/L_domain >= 1.0d0/(real(NumberXNodes) - 1.0d0)) then
            print *, "The debyeLength is really large, less nodes needed!"
            print *, "debyeLength is:", del_x
            print *, "L_domain is:", L_domain
            stop
        end if
        self%grid(1) = 0.0d0
        self%grid(NumberXNodes) = L_domain
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
        self%nodeVol(1) = 0.5d0 * self%dx_dl(1)
        self%nodeVol(NumberXNodes) = 0.5d0 * self%dx_dl(NumberXNodes-1)

    end subroutine constructSineGrid

    subroutine constructHalfSineGrid(self, del_x, L_domain)
        class(Domain), intent(in out) :: self
        real(real64), intent(in) :: del_x, L_domain
        integer(int32) :: i
        if (del_x/L_domain >= 1.0d0/(real(NumberXNodes) - 1.0d0)) then
            print *, "The debyeLength is really large, less nodes needed!"
            stop
        end if
        self%grid(1) = 0.0d0
        self%grid(NumberXNodes) = L_domain/2.0d0
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
        self%nodeVol(1) = 0.5d0 * self%dx_dl(1)
        self%nodeVol(NumberXNodes) = 0.5d0 * self%dx_dl(NumberXNodes-1)

        if (self%boundaryConditions(1) == 3 .or. self%boundaryConditions(NumberXNodes) == 3) then
            print *, "Mesh is not periodic, cannot have periodic boundary!"
            stop
        end if
    end subroutine constructHalfSineGrid

    subroutine constructExpHalfGrid(self, del_x, L_domain)
        class(Domain), intent(in out) :: self
        real(real64), intent(in) :: del_x, L_domain
        integer(int32) :: i
        real(real64) :: k, A, F
        self%grid(1) = 0.0d0
        self%grid(NumberXNodes) = L_domain
        F = 100.0d0 * del_x
        if (F/L_domain > 0.5d0) then
            stop "Debye length too small for exponential increase "
        end if
          
        F = 100.0d0 * del_x
        if (F/L_domain > 0.5d0) then

        end if
        k = 2.0d0 * LOG(L_domain/F - 1.0d0)/real(NumberXNodes-1)
        A = L_domain / (EXP(k * real(NumberXNodes-1, kind = real64)) -1.0d0)
        do i = 2,NumberXNodes-1
            self % grid(i) = A * (EXP(k * real(i-1, kind = real64)) - 1.0d0)
        end do
        do i = 1, NumberXNodes-1
            self%dx_dl(i) = self%grid(i+1) - self%grid(i)
        end do
        do i = 2, NumberXNodes-1
            self%nodeVol(i) = (self%dx_dl(i-1) + self%dx_dl(i))/2
        end do
        self%nodeVol(1) = 0.5d0 * self%dx_dl(1)
        self%nodeVol(NumberXNodes) = 0.5d0 * self%dx_dl(NumberXNodes-1)

        if (self%boundaryConditions(1) == 3 .or. self%boundaryConditions(NumberXNodes) == 3) then
            print *, "Mesh is not periodic, cannot have periodic boundary!"
            stop
        end if
    end subroutine constructExpHalfGrid

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
        self%nodeVol(1) = 0.5d0 * self%dx_dl(1)
        self%nodeVol(NumberXNodes) = 0.5d0 * self%dx_dl(NumberXNodes-1)
    end subroutine constructUniformGrid

    function getLFromX(self, x) result(l)
        class(Domain), intent(in) :: self
        real(real64), intent(in) :: x
        integer(int32) :: idxLower, idxHigher, idxMiddle
        real(real64) :: l
        idxLower = 1
        idxHigher = NumberXNodes
        if ((x<self%grid(1)) .or. (x > self%grid(NumberXNodes))) then
            print *, 'x value outside of grid range in getLFromX!'
            stop
        end if
        do while (idxLower /= idxHigher-1)
            idxMiddle = (idxLower + idxHigher)/2
            if (self%grid(idxMiddle) <= x) then
                idxLower = idxMiddle
            else
                idxHigher = idxMiddle
            end if
        end do
        l = idxLower + (x - self%grid(idxLower))/self%dx_dl(idxLower)
        
    end function getLFromX


    subroutine writeDomain(self, dirName)
        ! Writes domain data into binary file under Data
        class(Domain), intent(in) :: self
        character(*), intent(in) :: dirName
        open(41,file=dirName//'/domainGrid.dat', form='UNFORMATTED')
        write(41) self%grid
        close(41)
        open(41,file=dirName//'/domainDxDl.dat', form='UNFORMATTED')
        write(41) self%dx_dl
        close(41)
        open(41,file=dirName//'/domainNodeVol.dat', form='UNFORMATTED')
        write(41) self%nodeVol
        close(41)
    end subroutine writeDomain





end module mod_domain