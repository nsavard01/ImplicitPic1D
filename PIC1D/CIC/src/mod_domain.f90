module mod_domain
    use iso_fortran_env, only: int32, real64
    use mod_BasicFunctions
    use constants
    implicit none

    private
    public :: Domain, readWorld
    integer(int32), public, protected :: NumberXNodes = 10
    integer(int32), public, protected :: NumberXHalfNodes = 11

    ! domain contains arrays and values related to physical, logical dimensions of the spatial grid
    type :: Domain
        real(real64), allocatable :: grid(:) !array representing values at grid in m
        real(real64), allocatable :: dx_dl(:), centerDiff(:) !ratio of grid differences from physical to logical, assume logical separated by 1
        integer(int32), allocatable :: boundaryConditions(:), threadNodeIndx(:,:), threadHalfNodeIndx(:,:) ! Boundary condition flags for fields and particles
        real(real64) :: L_domain, startX, endX
        integer(int32) :: numThreadNodeIndx, numThreadHalfNodeIndx
        logical :: gridSmoothBool
        ! (>0 dirichlet, -2 Neumann, -3 periodic, <=-4 dielectric), 0 is default in-body condition 

    contains
        procedure, public, pass(self) :: constructHalfSineGrid
        procedure, public, pass(self) :: constructSineGrid
        procedure, public, pass(self) :: constructInvSineGrid
        procedure, public, pass(self) :: constructUniformGrid
        procedure, public, pass(self) :: constructGrid
        procedure, public, pass(self) :: constructExpHalfGrid
        procedure, public, pass(self) :: smoothField
        procedure, public, pass(self) :: getLFromX
        procedure, public, pass(self) :: getXFromL
        procedure, public, pass(self) :: writeDomain
    end type Domain


    interface Domain
        module procedure :: domain_constructor
    end interface Domain

contains
    ! Grid points at phi-calculation points, dx_dl is the distance between grid points in NGP, or distance between 
    ! Initialization procedures
    type(Domain) function domain_constructor(leftBoundary, rightBoundary) result(self)
        ! Construct domain object, initialize grid, dx_dl, and nodeVol.
        integer(int32), intent(in) :: leftBoundary, rightBoundary
        integer(int32) :: i, k, spacingThread, modThread
        allocate(self % grid(NumberXNodes), self % dx_dl(NumberXNodes), self%centerDiff(NumberXNodes-1), self%boundaryConditions(NumberXHalfNodes))
        self % grid = (/(i, i=1, NumberXNodes)/)
        self % dx_dl = 1.0d0
        self%gridSmoothBool = .false.
        self % boundaryConditions = 0
        self%boundaryConditions(1) = leftBoundary
        self%boundaryConditions(NumberXHalfNodes) = rightBoundary
        if (leftBoundary == 3 .or. rightBoundary == 3) then
            self%boundaryConditions(1) = 3
            self%boundaryConditions(NumberXHalfNodes) = 3
        end if

        if (numThread < NumberXNodes) then
            self%numThreadNodeIndx = numThread
        else
            self%numThreadNodeIndx = NumberXNodes
        end if

        allocate(self%threadNodeIndx(2, self%numThreadNodeIndx))
        spacingThread = NumberXNodes/self%numThreadNodeIndx - 1
        modThread = MOD(NumberXNodes, self%numThreadNodeIndx)
        k = 1
        do i = 1, self%numThreadNodeIndx
            self%threadNodeIndx(1, i) = k
            if (i <= modThread) then
                k = k + spacingThread + 1
            else
                k = k + spacingThread
            end if
            self%threadNodeIndx(2,i) = k
            k = k + 1
        end do
        
        if (numThread < NumberXHalfNodes) then
            self%numThreadHalfNodeIndx = numThread
        else
            self%numThreadHalfNodeIndx = NumberXHalfNodes
        end if
        allocate(self%threadHalfNodeIndx(2, self%numThreadHalfNodeIndx))
        spacingThread = NumberXHalfNodes/self%numThreadHalfNodeIndx - 1
        modThread = MOD(NumberXHalfNodes, self%numThreadHalfNodeIndx)
        k = 1
        do i = 1, self%numThreadHalfNodeIndx
            self%threadHalfNodeIndx(1, i) = k
            if (i <= modThread) then
                k = k + spacingThread + 1
            else
                k = k + spacingThread
            end if
            self%threadHalfNodeIndx(2,i) = k
            k = k + 1
        end do
        
    end function domain_constructor


    subroutine constructGrid(self, del_x, L_domain, gridType)
        class(Domain), intent(in out) :: self
        real(real64), intent(in) :: del_x, L_domain
        integer(int32), intent(in) :: gridType
        self%L_domain = L_domain
        SELECT CASE (gridType)
        CASE(0)
            call self%constructUniformGrid(L_domain)
        CASE(1)
            call self%constructSineGrid(del_x, L_domain)
        CASE(2)
            call self%constructHalfSineGrid(del_x, L_domain)
        CASE(3)
            call self%constructExpHalfGrid(del_x, L_domain)
        CASE(5)
            call self%constructInvSineGrid(del_x, L_domain)
        CASE default
            print *, "Gridtype", gridType, "doesn't exist!"
            stop
        END SELECT
        self%startX = self%grid(1) - 0.5d0 * self%dx_dl(1)
        self%endX = self%grid(NumberXNodes) + 0.5d0 * self%dx_dl(NumberXNodes)
        self%L_domain = self%endX - self%startX
    end subroutine constructGrid

    subroutine constructSineGrid(self, del_x, L_domain)
        class(Domain), intent(in out) :: self
        real(real64), intent(in) :: del_x, L_domain
        integer(int32) :: i, NumberNodes
        real(real64) :: gridField(NumberXNodes+1)
        NumberNodes = NumberXNodes + 1
        if (del_x/L_domain >= 1.0d0/(real(NumberNodes) - 1.0d0)) then
            print *, "The debyeLength is really large, less nodes needed!"
            print *, "debyeLength is:", del_x
            print *, "L_domain is:", L_domain
            stop
        end if
        gridField(1) = 0.0d0
        gridField(NumberNodes) = L_domain
        do i = 2,NumberNodes-1
            gridField(i) = L_domain * ((real(i)-1.0d0)/(real(NumberNodes) - 1.0d0) - (1.0d0/(real(NumberNodes) - 1.0d0) - del_x/L_domain) &
            * SIN(2 * pi * (i-1) / (NumberNodes - 1)) / SIN(2 * pi / (NumberNodes - 1)) )
        end do
        do i = 1, NumberNodes-1
            self%dx_dl(i) = gridField(i+1) - gridField(i)
        end do
        self%grid = 0.5d0 * (gridField(1:NumberNodes-1) + gridField(2:))
        self%centerDiff = self%grid(2:NumberXNodes) - self%grid(1:NumberXNodes-1)

    end subroutine constructSineGrid

    subroutine constructInvSineGrid(self, del_x, L_domain)
        class(Domain), intent(in out) :: self
        real(real64), intent(in) :: del_x, L_domain
        integer(int32) :: i
        real(real64) :: gridField(NumberXHalfNodes), phase
        if (del_x/L_domain >= 1.0d0/(real(NumberXHalfNodes) - 1.0d0)) then
            print *, "The debyeLength is really large, less nodes needed!"
            print *, "debyeLength is:", del_x
            print *, "L_domain is:", L_domain
            stop
        end if
        gridField(1) = 0.0d0
        gridField(NumberXHalfNodes) = L_domain
        do i = 2,NumberXNodes
            phase = real(i-1)/real(NumberXNodes)
            gridField(i) = L_domain * (phase + (1.0d0 - real(NumberXHalfNodes) * del_x/L_domain)* SIN(2.0d0 * pi * phase)/ 2.0d0 / pi )
        end do
        do i = 1, NumberXNodes
            self%dx_dl(i) = gridField(i+1) - gridField(i)
        end do
        self%grid = 0.5d0 * (gridField(1:NumberXNodes) + gridField(2:))
        self%centerDiff = self%grid(2:NumberXNodes) - self%grid(1:NumberXNodes-1)

    end subroutine constructInvSineGrid

    subroutine constructHalfSineGrid(self, del_x, L_domain)
        class(Domain), intent(in out) :: self
        real(real64), intent(in) :: del_x, L_domain
        integer(int32) :: i, NumberNodes
        real(real64) :: gridField(NumberXNodes+1)
        NumberNodes = NumberXNodes + 1
        if (del_x/L_domain >= 1.0d0/(real(NumberNodes) - 1.0d0)) then
            print *, "The debyeLength is really large, less nodes needed!"
            stop
        end if
        gridField(1) = 0.0d0
        gridField(NumberNodes) = L_domain/2.0d0
        do i = 2,NumberNodes-1
            gridField(i) = L_domain * ((real(i)-1.0d0)/(real(NumberNodes) - 1.0d0)/2.0d0 - (1.0d0/(real(NumberNodes) - 1.0d0)/2.0d0 - del_x/L_domain) &
            * SIN(pi * (real(i)-1)/(real(NumberNodes) - 1)) / SIN(pi / (real(NumberNodes) - 1.0d0)) )
        end do
        do i = 1, NumberNodes-1
            self%dx_dl(i) = gridField(i+1) - gridField(i)
        end do
        self%grid = 0.5d0 * (gridField(1:NumberNodes-1) + gridField(2:))
        self%centerDiff = self%grid(2:NumberXNodes) - self%grid(1:NumberXNodes-1)

    end subroutine constructHalfSineGrid

    subroutine constructExpHalfGrid(self, del_x, L_domain)
        class(Domain), intent(in out) :: self
        real(real64), intent(in) :: del_x, L_domain
        integer(int32) :: i, NumberNodes
        real(real64) :: k, A, F, gridField(NumberXNodes+1)
        NumberNodes = NumberXNodes + 1
        gridField(1) = 0.0d0
        gridField(NumberNodes) = L_domain
        F = 100.0d0 * del_x
        if (F/L_domain > 0.5d0) then
            stop "Debye length too small for exponential increase "
        end if
          
        F = 100.0d0 * del_x
        if (F/L_domain > 0.5d0) then

        end if
        k = 2.0d0 * LOG(L_domain/F - 1.0d0)/real(NumberNodes-1)
        A = L_domain / (EXP(k * real(NumberNodes-1, kind = real64)) -1.0d0)
        do i = 2,NumberNodes-1
            gridField(i) = A * (EXP(k * real(i-1, kind = real64)) - 1.0d0)
        end do
        do i = 1, NumberNodes-1
            self%dx_dl(i) = self%grid(i+1) - self%grid(i)
        end do
        self%grid = 0.5d0 * (gridField(1:NumberNodes-1) + gridField(2:))
        self%centerDiff = self%grid(2:NumberXNodes) - self%grid(1:NumberXNodes-1)
        if (self%boundaryConditions(1) == 3 .or. self%boundaryConditions(NumberXNodes) == 3) then
            print *, "Mesh is not periodic, cannot have periodic boundary!"
            stop
        end if
    end subroutine constructExpHalfGrid

    subroutine constructUniformGrid(self, L_domain)
        class(Domain), intent(in out) :: self
        real(real64), intent(in) :: L_domain
        integer(int32) :: i, NumberNodes
        real(real64) :: gridField(NumberXNodes+1)
        NumberNodes = NumberXNodes + 1
        gridField(1) = 0.0d0
        gridField(NumberNodes) = L_domain
        do i = 2, NumberNodes-1
            gridField(i) =  (i-1) * L_domain / (NumberNodes - 1)
        end do
        do i = 1, NumberNodes-1
            self%dx_dl(i) = gridField(i+1) - gridField(i)
        end do
        self%grid = 0.5d0 * (gridField(1:NumberNodes-1) + gridField(2:))
        self%centerDiff = self%grid(2:NumberXNodes) - self%grid(1:NumberXNodes-1)
    end subroutine constructUniformGrid

    function getLFromX(self, x) result(l)
        class(Domain), intent(in) :: self
        real(real64), intent(in) :: x
        integer(int32) :: idxLower, idxHigher, idxMiddle
        real(real64) :: l
        idxLower = 1
        idxHigher = NumberXNodes
        if ((x<self%grid(1) - 0.5d0 * self%dx_dl(1)) .or. (x > self%grid(NumberXNodes) + 0.5d0 * self%dx_dl(NumberXNodes))) then
            print *, 'x value outside of grid range in getLFromX!'
            print *, 'x is:', x
            stop
        end if
        do while (idxLower /= idxHigher-1)
            idxMiddle = (idxLower + idxHigher)/2
            if (self%grid(idxMiddle) - 0.5d0 * self%dx_dl(idxMiddle) <= x) then
                idxLower = idxMiddle
            else
                idxHigher = idxMiddle
            end if
        end do
        if (self%grid(idxLower) + 0.5d0 * self%dx_dl(idxLower) <= x) then
            idxLower = idxHigher
        end if
        l = idxLower + 0.5d0 + (x - self%grid(idxLower))/self%dx_dl(idxLower)
        
    end function getLFromX

    function getXFromL(self, l) result(x)
        class(Domain), intent(in) :: self
        real(real64), intent(in) :: l
        real(real64) :: x
        x = self%grid(INT(l)) + self%dx_dl(INT(l)) * (l - INT(l) - 0.5d0)
        
    end function getXFromL

    subroutine smoothField(self, rawField, newField)
        class(Domain), intent(in) :: self
        real(real64), intent(in) :: rawField(NumberXHalfNodes)
        real(real64), intent(in out) :: newField(NumberXHalfNodes)
        integer(int32) :: i, boundVal
        ! smoothing of fields
        SELECT CASE (self%boundaryConditions(1))
        CASE(1,4)
            newField(1) = 0.25d0 * (2.0d0 * rawField(1) + 2.0d0 * rawField(2))
        CASE(2)
            newField(1) = 0.0d0
        CASE(3)
            newField(1) = 0.25d0 * (rawField(NumberXHalfNodes-1) + 2.0d0 * rawField(1) + rawField(2))
        END SELECT
        SELECT CASE (self%boundaryConditions(NumberXHalfNodes))
        CASE(1,4)
            newField(NumberXHalfNodes) = 0.25d0 * (2.0d0 * rawField(NumberXHalfNodes) + 2.0d0 * rawField(NumberXHalfNodes-1))
        CASE(2)
            newField(NumberXHalfNodes) = 0.0d0
        CASE(3)
            newField(NumberXHalfNodes) = newField(1)
        END SELECT
        newField(2:NumberXNodes) = 0.25d0 * (rawField(1:NumberXNodes-1) + 2.0d0 * rawField(2:NumberXNodes) + rawField(3:NumberXHalfNodes))

    end subroutine smoothField

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
        open(41,file=dirName//"/domainBoundaryConditions.dat", form='UNFORMATTED')
        write(41) self%boundaryConditions
        close(41)
    end subroutine writeDomain

    ! ------------------------ read world inputs ---------------------------------

    subroutine readWorld(GeomFilename, world, T_e, n_ave)
        type(Domain), intent(in out) :: world
        character(len=*), intent(in) :: GeomFilename
        real(real64), intent(in) :: T_e, n_ave
        integer(int32) :: io, leftBoundary, rightBoundary, gridType, tempInt, i, smoothInt
        real(real64) :: debyeLength, L_domain
        integer(int32), allocatable :: boundArray(:)
        print *, ""
        print *, "Reading domain inputs:"
        print *, "------------------"
        if (.not. restartBool) then
            open(10,file='../InputData/'//GeomFilename)
        else
            open(10,file=restartDirectory//'/InputData/'//GeomFilename)
        end if
        read(10, *, IOSTAT = io) NumberXNodes
        read(10, *, IOSTAT = io) L_domain
        read(10, *, IOSTAT = io) gridType, debyeLength
        read(10, *, IOSTAT = io) leftBoundary, rightBoundary
        read(10, *, IOSTAT = io) smoothInt
        read(10, *, IOSTAT = io)
        read(10, *, IOSTAT = io)
        close(10)
        NumberXHalfNodes = NumberXNodes + 1
        debyeLength = MAX(0.1d0 * getDebyeLength(T_e, n_ave), debyeLength)
        if ((leftBoundary == 3) .or. (rightBoundary == 3)) then
            leftBoundary = 3
            rightBoundary = 3
        end if
        world = Domain(leftBoundary, rightBoundary)
        call world % constructGrid(debyeLength, L_domain, gridType)
        if (smoothInt /= 0) world%gridSmoothBool = .true.
        print *, "Number of nodes:", NumberXNodes
        print *, 'Number of half nodes:', NumberXHalfNodes
        print *, "Grid length:", world%L_domain
        print *, 'gridType:', gridType
        print *, "Left boundary type:", world%boundaryConditions(1)
        print *, "Right boundary type:", world%boundaryConditions(NumberXNodes+1)
        print *, 'smallest delX:', MINVAL(world%dx_dl)
        print *, 'Binomial smoothing:', world%gridSmoothBool
        print *, "------------------"
        print *, ""

    end subroutine readWorld



end module mod_domain