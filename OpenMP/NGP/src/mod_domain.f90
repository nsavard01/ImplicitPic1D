module mod_domain
    use iso_fortran_env, only: int32, real64
    use mod_BasicFunctions
    use constants
    implicit none

    private
    public :: Domain, readWorld
    integer(int32), public, protected :: NumberXNodes = 10
    integer(int32), public, protected :: NumberXHalfNodes = 9

    ! domain contains arrays and values related to physical, logical dimensions of the spatial grid
    type :: Domain
        real(real64), allocatable :: grid(:) !array representing values at grid in m
        real(real64), allocatable :: dx_dl(:) !ratio of grid differences from physical to logical, assume logical separated by 1
        real(real64), allocatable :: nodeVol(:) !node vol, or difference between half grid in physical space
        integer(int32), allocatable :: boundaryConditions(:), threadNodeIndx(:,:), threadHalfNodeIndx(:,:) ! Boundary condition flags for fields and particles
        real(real64) :: L_domain, startX, endX
        integer(int32) :: numThreadNodeIndx, numThreadHalfNodeIndx
        ! (>0 dirichlet, -2 Neumann, -3 periodic, <=-4 dielectric), 0 is default in-body condition 

    contains
        procedure, public, pass(self) :: constructHalfSineGrid
        procedure, public, pass(self) :: constructSineGrid
        procedure, public, pass(self) :: constructUniformGrid
        procedure, public, pass(self) :: constructGrid
        procedure, public, pass(self) :: constructExpHalfGrid
        procedure, public, pass(self) :: constructHalfEvenHalfSinusoid
        procedure, public, pass(self) :: makeArraysFromGrid
        procedure, public, pass(self) :: addThreadedDomainArray
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
        integer(int32) :: i, k, spacingThread, modThread
        allocate(self % grid(NumberXNodes), self % dx_dl(NumberXHalfNodes), self % nodeVol(NumberXNodes), self%boundaryConditions(NumberXNodes))
        self % grid = (/(i, i=1, NumberXNodes)/)
        self % dx_dl = 1.0d0
        self % nodeVol = 1.0d0
        self % boundaryConditions = 0
        self%boundaryConditions(1) = leftBoundary
        self%boundaryConditions(NumberXNodes) = rightBoundary
        if (leftBoundary == 3 .or. rightBoundary == 3) then
            self%boundaryConditions(1) = 3
            self%boundaryConditions(NumberXNodes) = 3
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
        SELECT CASE (gridType)
        CASE(0)
            call self%constructUniformGrid(L_domain)
        CASE(1)
            call self%constructSineGrid(del_x, L_domain)
        CASE(2)
            call self%constructHalfSineGrid(del_x, L_domain)
        CASE(3)
            call self%constructExpHalfGrid(del_x, L_domain)
        CASE(4)
            call self%constructHalfEvenHalfSinusoid(del_x, L_domain)
        CASE default
            print *, "Gridtype", gridType, "doesn't exist!"
            stop
        END SELECT
        call self%makeArraysFromGrid()
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
        do i = 2,NumberXHalfNodes
            self % grid(i) = L_domain * ((real(i)-1.0d0)/(real(NumberXNodes) - 1.0d0) - (1.0d0/(real(NumberXNodes) - 1.0d0) - del_x/L_domain) &
            * SIN(2 * pi * (i-1) / (NumberXNodes - 1)) / SIN(2 * pi / (NumberXNodes - 1)) )
        end do

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
        do i = 2,NumberXHalfNodes
            self % grid(i) = L_domain * ((real(i)-1.0d0)/(real(NumberXNodes) - 1.0d0)/2.0d0 - (1.0d0/(real(NumberXNodes) - 1.0d0)/2.0d0 - del_x/L_domain) &
            * SIN(pi * (real(i)-1)/(real(NumberXNodes) - 1)) / SIN(pi / (real(NumberXNodes) - 1.0d0)) )
        end do
        

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
        k = 2.0d0 * LOG(L_domain/F - 1.0d0)/real(NumberXHalfNodes)
        A = L_domain / (EXP(k * real(NumberXHalfNodes, kind = real64)) -1.0d0)
        do i = 2,NumberXHalfNodes
            self % grid(i) = A * (EXP(k * real(i-1, kind = real64)) - 1.0d0)
        end do

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
        do i = 2, NumberXHalfNodes
            self % grid(i) =  (i-1) * L_domain / (NumberXHalfNodes)
        end do
    end subroutine constructUniformGrid

    subroutine constructHalfEvenHalfSinusoid(self, del_x, L_domain)
        class(Domain), intent(in out) :: self
        real(real64), intent(in) :: del_x, L_domain
        integer(int32) :: numEvenCells, numSinCells, i
        real(real64) :: leftX, rightX, evenDelX, middleDelX
        numEvenCells = (NumberXNodes-1)/4 
        numSinCells = NumberXNodes - 1 - 2 * numEvenCells
        evenDelX = del_x/numEvenCells
        self%grid(1) = 0.0d0
        self%grid(NumberXNodes) = L_domain
        do i = 1, numEvenCells
            self%grid(i+1) = self%grid(i) + evenDelX
            self%grid(NumberXNodes-i) = self%grid(NumberXNodes-i+1) - evenDelX
        end do
        middleDelX = self%grid(NumberXNodes-numEvenCells) - self%grid(numEvenCells+1)
        do i = 2,numSinCells
            self % grid(i+numEvenCells) = self%grid(numEvenCells+1) + middleDelX * ((real(i)-1.0d0)/(real(numSinCells)) - (1.0d0/(real(numSinCells)) - evenDelX/middleDelX) &
            * SIN(2 * pi * (i-1) / (numSinCells)) / SIN(2 * pi / (numSinCells)) )
        end do
        


    end subroutine constructHalfEvenHalfSinusoid

    subroutine makeArraysFromGrid(self)
        class(Domain), intent(in out) :: self
        integer(int32) :: i
        do i = 1, NumberXHalfNodes
            self%dx_dl(i) = self%grid(i+1) - self%grid(i)
        end do
        do i = 2, NumberXHalfNodes
            self%nodeVol(i) = (self%dx_dl(i-1) + self%dx_dl(i))/2.0d0
        end do
        self%nodeVol(1) = 0.5d0 * self%dx_dl(1)
        self%nodeVol(NumberXNodes) = 0.5d0 * self%dx_dl(NumberXHalfNodes)
        self%startX = self%grid(1)
        self%endX = self%grid(NumberXNodes)
        self%L_domain = self%endX - self%startX
    end subroutine makeArraysFromGrid

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

    subroutine addThreadedDomainArray(self, array_add, x, N_x, iThread, const)
        ! Take array on grid nodes of half nodes x with second dimension thread count and add to array_add of same domain dimension using Openmp
        class(Domain), intent(in) :: self
        real(real64), intent(in out) :: array_add(N_x)
        real(real64), intent(in) :: x(NumberXNodes, numThread), const
        integer(int32), intent(in) :: iThread, N_x
        if (N_x == NumberXNodes .and. iThread <= self%numThreadNodeIndx) then
            array_add(self%threadNodeIndx(1,iThread):self%threadNodeIndx(2,iThread)) = array_add(self%threadNodeIndx(1,iThread):self%threadNodeIndx(2,iThread)) &
            + SUM(x(self%threadNodeIndx(1,iThread):self%threadNodeIndx(2,iThread), :), DIM=2) * const
        else if (N_x == NumberXHalfNodes .and. iThread <= self%numThreadHalfNodeIndx) then
            array_add(self%threadHalfNodeIndx(1,iThread):self%threadHalfNodeIndx(2,iThread)) = array_add(self%threadHalfNodeIndx(1,iThread):self%threadHalfNodeIndx(2,iThread)) &
            + SUM(x(self%threadHalfNodeIndx(1,iThread):self%threadHalfNodeIndx(2,iThread), :), DIM=2) * const
        end if
    end subroutine addThreadedDomainArray


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
        open(41,file=dirName//"/domainBoundaryConditions.dat", form='UNFORMATTED')
        write(41) self%boundaryConditions
        close(41)
    end subroutine writeDomain

    ! ----------------------- read inputs ------------------------------------

    subroutine readWorld(GeomFilename, world, T_e, n_ave)
        type(Domain), intent(in out) :: world
        character(len=*), intent(in) :: GeomFilename
        real(real64), intent(in) :: T_e, n_ave
        integer(int32) :: io, leftBoundary, rightBoundary, gridType, intArray(20), i
        real(real64) :: debyeLength, L_domain
        integer(int32), allocatable :: boundArray(:)
        print *, ""
        print *, "Reading domain inputs:"
        print *, "------------------"
        open(10,file='../InputData/'//GeomFilename)
        read(10, *, IOSTAT = io) NumberXNodes
        read(10, *, IOSTAT = io) L_domain
        read(10, *, IOSTAT = io) gridType, debyeLength
        read(10, *, IOSTAT = io) leftBoundary, rightBoundary
        read(10, *, IOSTAT = io)
        read(10, *, IOSTAT = io)
        read(10, *, IOSTAT = io)
        close(10)
        if (restartBool) then
            open(10,file=restartDirectory//"/"//"InitialConditions.dat", IOSTAT=io)
            read(10, *, IOSTAT = io)
            read(10, *, IOSTAT = io) intArray(1), NumberXNodes
            close(10) 
            allocate(boundArray(NumberXNodes))
            open(10,file=restartDirectory//"/"//"domainBoundaryConditions.dat", form='UNFORMATTED', IOSTAT=io)
            read(10, IOSTAT = io) boundArray
            close(10)
            leftBoundary = boundArray(1)
            rightBoundary = boundArray(NumberXNodes)
            deallocate(boundArray)
        end if
        NumberXHalfNodes = NumberXNodes-1
        debyeLength = MAX(getDebyeLength(T_e, n_ave), debyeLength)
        if ((leftBoundary == 3) .or. (rightBoundary == 3)) then
            leftBoundary = 3
            rightBoundary = 3
        end if
        world = Domain(leftBoundary, rightBoundary)
        if (.not. restartBool) then
            call world % constructGrid(debyeLength, L_domain, gridType)
        else
            open(10,file=restartDirectory//"/"//"domainGrid.dat", form='UNFORMATTED', IOSTAT=io)
            read(10, IOSTAT = io) world%grid
            close(10)
            call world%makeArraysFromGrid()
        end if
        print *, "Number of nodes:", NumberXNodes
        print *, "Number of half nodes:", NumberXHalfNodes
        print *, "Grid length:", world%L_domain
        print *, 'gridType:', gridType
        print *, "Left boundary type:", world%boundaryConditions(1)
        print *, "Right boundary type:", world%boundaryConditions(NumberXNodes)
        print *, 'smallest delX:', MINVAL(world%dx_dl)
        print *, "------------------"
        print *, ""

    end subroutine readWorld



end module mod_domain