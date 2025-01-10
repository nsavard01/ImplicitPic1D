module mod_domain
    use iso_fortran_env, only: int32, real64
    use mod_BasicFunctions
    use constants
    implicit none

    private
    public :: Domain, readWorld
    integer(int32), public, protected :: NumberXNodes = 10
    integer(int32), public, protected :: NumberXHalfNodes = 9

    ! Class which stores properties of domain 
    type :: Domain
        real(real64), allocatable :: grid(:) !array representing values at grid in m
        real(real64), allocatable :: dx_dl(:) ! l for logical, cell sizes
        real(real64), allocatable :: nodeVol(:) ! node volume, size of 
        integer(int32), allocatable :: boundaryConditions(:), threadNodeIndx(:,:), threadHalfNodeIndx(:,:) ! Boundary conditions, node indices divided into threads for OpenMP
        real(real64) :: L_domain, startX, endX
        real(real64), allocatable :: curv_params(:)
        integer(int32) :: numThreadNodeIndx, numThreadHalfNodeIndx, num_even_edge_cells, grid_type
        logical :: gridSmoothBool

    contains
        procedure, public, pass(self) :: get_jacobian
        procedure, public, pass(self) :: smoothField
        procedure, public, pass(self) :: getLFromX
        procedure, public, pass(self) :: getXFromL
        procedure, public, pass(self) :: writeDomain
    end type Domain


    interface Domain
        module procedure :: domain_constructor
    end interface Domain

contains

    type(Domain) function domain_constructor(leftBoundary, rightBoundary, min_del_x, grid_type, num_even_edge_cells, L_domain) result(self)
        ! Construct domain object
        integer(int32), intent(in) :: leftBoundary, rightBoundary, num_even_edge_cells, grid_type
        real(real64), intent(in) :: min_del_x, L_domain
        integer(int32) :: i, k, spacingThread, modThread
        allocate(self % grid(NumberXNodes), self % dx_dl(NumberXHalfNodes), self % nodeVol(NumberXNodes), self%boundaryConditions(NumberXNodes))
        self % grid = (/(i, i=1, NumberXNodes)/)
        self%gridSmoothBool = .false.
        self % dx_dl = 1.0d0
        self % nodeVol = 1.0d0
        self % boundaryConditions = 0
        self%L_domain = L_domain
        self%num_even_edge_cells = num_even_edge_cells
        self%grid_type = grid_type
        self%min_del_x = min_del_x
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
        self%startX = 0.0d0
        self%endX = self%L_domain
        self%grid(1) = self%startX
        self%grid(NumberXNodes) = self%L_domain
        if (self%grid_type == 4) then
            self%grid(1 + self%num_even_edge_cells) = self%min_del_x
            self%grid(NumberXNodes - self%num_even_edge_cells) = self%L_domain - self%min_del_x
            self%min_del_x = self%min_del_x / self%num_even_edge_cells
        end if
        
        do i = 2, NumberXHalfNodes
            self%grid(i) = self%getXFromL(real(i, kind = 8))
        end do
        
        do i = 1, NumberXHalfNodes
            self%dx_dl(i) = self%get_jacobian(real(i, kind = 8) + 0.5d0)
        end do
        do i = 1, NumberXNodes
            self%nodeVol(i) = self%get_jacobian(real(i, kind = 8))
        end do
        

    end function domain_constructor


    function get_jacobian(self, xi) result(res)
        class(Domain), intent(in) :: self
        real(real64), intent(in) :: xi
        real(real64) :: res, N_cells, length_scale

        select case (self%grid_type)
        case(0)
            res = self%min_del_x
        case(1)
            res = self%L_domain * (1.0d0/real(NumberXHalfNodes) - (2.0d0 * pi / real(NumberXHalfNodes)) * (1.0d0/real(NumberXHalfNodes) - self%min_del_x/self%L_domain) &
            * COS(2 * pi * (xi-1.0d0) / real(NumberXHalfNodes)) / SIN(2 * pi / real(NumberXHalfNodes)) )
        case(4)
            if (xi > self%num_even_edge_cells + 1 .and. xi < NumberXNodes - self%num_even_edge_cells) then
                N_cells = real(NumberXHalfNodes - 2 * self%num_even_edge_cells, kind = 8)
                length_scale = (self%grid(NumberXNodes - self%num_even_edge_cells) - self%grid(1 + self%num_even_edge_cells)) / N_cells

                res = length_scale * ( 1.0d0 - cos(2.0d0 * pi * (xi - 1.0d0 - self%num_even_edge_cells) / N_cells) * (1.0d0 - self%min_del_x/length_scale) )
            else
                res = self%min_del_x
            end if

        end select
    end function get_jacobian


    function getXFromL(self, l) result(x)
        ! get physical coordinate from logical
        class(Domain), intent(in) :: self
        real(real64), intent(in) :: l
        real(real64) :: x, N_cells, length_scale

        select case (self%grid_type)
        case(0)
            x = self%min_del_x * (l - 1.0d0) + self%startX
        case(1)
            x = self%L_domain * ((l-1.0d0)/real(NumberXHalfNodes) - (1.0d0/real(NumberXHalfNodes) - self%min_del_x/self%L_domain) &
            * SIN(2 * pi * (l-1.0d0) / real(NumberXHalfNodes)) / SIN(2 * pi / real(NumberXHalfNodes)) )
        case(4)
            ! Assume have already defined grid points where two grids divide
            if (l > self%num_even_edge_cells + 1 .and. l < NumberXNodes - self%num_even_edge_cells) then
                N_cells = real(NumberXHalfNodes - 2 * self%num_even_edge_cells, kind = 8)
                length_scale = (self%grid(NumberXNodes - self%num_even_edge_cells) - self%grid(1 + self%num_even_edge_cells)) / N_cells

                x = self%grid(self%num_even_edge_cells+1) + length_scale * ( l-1.0d0 - self%num_even_edge_cells - 0.5d0 * N_cells * sin(2.0d0 * pi * (l-1.0d0 - self%num_even_edge_cells)/N_cells) * &
                    (1.0d0 - self%min_del_x / length_scale)/pi)
            else if (l <= self%num_even_edge_cells + 1) then
                x = self%min_del_x * (l - 1.0d0) + self%startX
            else
                x = self%endX - (real(NumberXNodes, kind = 8) - l) * self%min_del_x
            end if
        end select

    end function getXFromL

    function getLFromX(self, x) result(l)
        ! get logical from physical space, use binary search
        ! then 
        class(Domain), intent(in) :: self
        real(real64), intent(in) :: x
        integer(int32) :: idxLower, idxHigher, idxMiddle
        real(real64) :: l, l_prev
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
        l_prev = idxLower + (x - self%grid(idxLower))/self%dx_dl(idxLower)
        l = l_prev + (x - self%getXFromL(l_prev))/self%get_jacobian(l_prev)
        !picard iteration (assume close enough initial condition to get to solution)
        do while (abs(l - l_prev) > 1.d-8)
            l_prev = l
            l = l_prev + (x - self%getXFromL(l_prev))/self%get_jacobian(l_prev)
        end do

        
    end function getLFromX

    subroutine smoothField(self, rawField, newField)
        ! Smooth fields on grid with quadratic smoothing
        class(Domain), intent(in) :: self
        real(real64), intent(in) :: rawField(NumberXHalfNodes)
        real(real64), intent(in out) :: newField(NumberXHalfNodes)
        integer(int32) :: i, boundVal
        do i = 1, NumberXHalfNodes
            boundVal = self%boundaryConditions(i+1) - self%boundaryConditions(i)
            SELECT CASE(boundVal)
            CASE(0)
                newField(i) = 0.25d0 * (rawField(i-1) + 2.0d0 * rawField(i) + rawField(i+1))
            CASE(-1,-4)
                newField(i) = 0.25d0 * (3.0d0 * rawField(1) + rawField(2))
            CASE(1,4)
                newField(i) = 0.25d0 * (3.0d0 * rawField(NumberXHalfNodes) + rawField(NumberXHalfNodes-1))
            CASE(-2)
                newField(i) = 0.25d0 * (rawField(1) + rawField(2))
            CASE(2)
                newField(i) = 0.25d0 * (rawField(NumberXHalfNodes) + rawField(NumberXHalfNodes-1))
            CASE(-3)
                newField(i) = 0.25d0 * (rawField(NumberXHalfNodes) + 2.0d0 * rawField(1) + rawField(2))
            CASE(3)
                newField(i) = 0.25d0 * (rawField(NumberXHalfNodes-1) + 2.0d0 * rawField(NumberXHalfNodes) + rawField(1))
            END SELECT
        end do

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
        open(41,file=dirName//'/domainNodeVol.dat', form='UNFORMATTED')
        write(41) self%nodeVol
        close(41)
        open(41,file=dirName//"/domainBoundaryConditions.dat", form='UNFORMATTED')
        write(41) self%boundaryConditions
        close(41)
    end subroutine writeDomain

    ! ----------------------- read inputs ------------------------------------

    subroutine readWorld(GeomFilename, world, T_e, n_ave)
        ! read world input
        type(Domain), intent(in out) :: world
        character(len=*), intent(in) :: GeomFilename
        real(real64), intent(in) :: T_e, n_ave
        integer(int32) :: io, leftBoundary, rightBoundary, gridType, intArray(20), i, smoothInt, extra_int
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
        read(10, *, IOSTAT = io) gridType, debyeLength, extra_int
        read(10, *, IOSTAT = io) leftBoundary, rightBoundary
        read(10, *, IOSTAT = io) smoothInt
        read(10, *, IOSTAT = io)
        read(10, *, IOSTAT = io)
        close(10)
        NumberXHalfNodes = NumberXNodes-1
        debyeLength = MAX(0.1d0 * getDebyeLength(T_e, n_ave), debyeLength)
        if ((leftBoundary == 3) .or. (rightBoundary == 3)) then
            leftBoundary = 3
            rightBoundary = 3
        end if
        world = Domain(leftBoundary, rightBoundary, debyeLength, gridType, extra_int, L_domain)
        if (smoothInt /= 0) world%gridSmoothBool = .true.
        print *, "Number of nodes:", NumberXNodes
        print *, "Number of half nodes:", NumberXHalfNodes
        print *, "Grid length:", world%L_domain
        print *, 'gridType:', gridType
        print *, "Left boundary type:", world%boundaryConditions(1)
        print *, "Right boundary type:", world%boundaryConditions(NumberXNodes)
        print *, 'smallest delX:', MINVAL(world%dx_dl)
        print *, 'Binomial smoothing:', world%gridSmoothBool
        print *, "------------------"
        print *, ""

    end subroutine readWorld



end module mod_domain