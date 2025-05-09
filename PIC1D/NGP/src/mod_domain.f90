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
        integer(int32) :: numThreadNodeIndx, numThreadHalfNodeIndx
        logical :: gridSmoothBool

    contains
        procedure, public, pass(self) :: constructHalfSineGrid
        procedure, public, pass(self) :: constructSineGrid
        procedure, public, pass(self) :: constructInvSineGrid
        procedure, public, pass(self) :: constructUniformGrid
        procedure, public, pass(self) :: constructGrid
        procedure, public, pass(self) :: constructExpHalfGrid
        procedure, public, pass(self) :: construct_uniformCenter_InvSineGrid
        procedure, public, pass(self) :: constructHalfEvenHalfSinusoid
        procedure, public, pass(self) :: construct_even_grid_doubling_region
        procedure, public, pass(self) :: construct_tanh_grid
        procedure, public, pass(self) :: makeArraysFromGrid
        procedure, public, pass(self) :: smoothField
        procedure, public, pass(self) :: getLFromX
        procedure, public, pass(self) :: getXFromL
        procedure, public, pass(self) :: writeDomain
    end type Domain


    interface Domain
        module procedure :: domain_constructor
    end interface Domain

contains

    type(Domain) function domain_constructor(leftBoundary, rightBoundary) result(self)
        ! Construct domain object
        integer(int32), intent(in) :: leftBoundary, rightBoundary
        integer(int32) :: i, k, spacingThread, modThread
        allocate(self % grid(NumberXNodes), self % dx_dl(NumberXHalfNodes), self % nodeVol(NumberXNodes), self%boundaryConditions(NumberXNodes))
        self % grid = (/(i, i=1, NumberXNodes)/)
        self%gridSmoothBool = .false.
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


    subroutine constructGrid(self, del_x, L_domain, gridType, numEvenCells, extra_real)
        ! Construct Grid
        class(Domain), intent(in out) :: self
        real(real64), intent(in) :: del_x, L_domain, extra_real
        integer(int32), intent(in) :: gridType, numEvenCells
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
            call self%constructHalfEvenHalfSinusoid(del_x, L_domain, numEvenCells)
        CASE(5)
            call self%constructInvSineGrid(del_x, L_domain)
        CASE(6)
            call self%construct_even_grid_doubling_region(del_x, L_domain, numEvenCells, extra_real)
        CASE(7)
            call self%construct_uniformCenter_InvSineGrid(del_x, L_domain, numEvenCells)
        CASE(8)
            call self%construct_tanh_grid(del_x, L_domain, numEvenCells, extra_real)
        CASE default
            print *, "Gridtype", gridType, "doesn't exist!"
            stop
        END SELECT
        call self%makeArraysFromGrid()
    end subroutine constructGrid

    subroutine constructSineGrid(self, del_x, L_domain)
        ! Sine grid with grid compressed at edges, large in center
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

    subroutine constructInvSineGrid(self, del_x, L_domain)
        ! Sin grid compressed in center, larger at edges
        class(Domain), intent(in out) :: self
        real(real64), intent(in) :: del_x, L_domain
        integer(int32) :: i
        real(real64) :: phase
        if (del_x/L_domain >= 1.0d0/(real(NumberXNodes) - 1.0d0)) then
            print *, "The debyeLength is really large, less nodes needed!"
            print *, "debyeLength is:", del_x
            print *, "L_domain is:", L_domain
            stop
        end if
        self%grid(1) = 0.0d0
        self%grid(NumberXNodes) = L_domain
        do i = 2,NumberXHalfNodes
            phase = real(i-1)/real(NumberXHalfNodes)
            self % grid(i) = L_domain * (phase + (1.0d0 - real(NumberXNodes) * del_x/L_domain)* SIN(2.0d0 * pi * phase)/ 2.0d0 / pi )
        end do
    end subroutine constructInvSineGrid

    subroutine construct_uniformCenter_InvSineGrid(self, del_x, L_domain, numEvenCells)
        ! Sin grid compressed in center, larger at edges
        ! Add center section which is uniform for 
        class(Domain), intent(in out) :: self
        real(real64), intent(in) :: del_x, L_domain
        integer(int32), intent(in) :: numEvenCells
        integer(int32) :: i, num_edge_cells
        real(real64) :: phase, uniform_size, edge_size, uniform_left_location, uniform_right_location, num_edge_real, del_x_min, top, bottom
        if (del_x/L_domain >= 1.0d0/(real(NumberXNodes) - 1.0d0)) then
            print *, "The debyeLength is really large, less nodes needed!"
            print *, "debyeLength is:", del_x
            print *, "L_domain is:", L_domain
            stop
        end if

        if (MOD(NumberXHalfNodes,2) /= MOD(numEvenCells,2)) then
            print *, 'Number of center cells must match (in terms of odd/even) number of total cells'
            stop
        end if

        self%grid(1) = 0.0d0
        self%grid(NumberXNodes) = L_domain
        if (MOD(NumberXHalfNodes,2) == 0) then
            ! even amount of cells
            self%grid(NumberXHalfNodes/2 + 1) = 0.5d0 * L_domain
            do i = 1, numEvenCells/2
                self%grid(NumberXHalfNodes/2 + 1 + i) = self%grid(NumberXHalfNodes/2 + 1) + del_x * i
                self%grid(NumberXHalfNodes/2 + 1 - i) = self%grid(NumberXHalfNodes/2 + 1) - del_x * i
            end do
            edge_size = 2.0d0 * self%grid(NumberXHalfNodes/2 + 1 - numEvenCells/2)
        else
            ! odd number of cells
            self%grid(NumberXNodes/2) = 0.5d0 * (L_domain - del_x)
            self%grid(NumberXNodes/2 + 1) = 0.5d0 * (L_domain + del_x)
            do i = 1, (numEvenCells-1)/2
                self%grid(NumberXNodes/2 + 1 + i) = self%grid(NumberXNodes/2 + 1) + del_x * i
                self%grid(NumberXNodes/2 - i) = self%grid(NumberXNodes/2) - del_x * i
            end do
            edge_size = 2.0d0 * self%grid(NumberXNodes/2 - (numEvenCells-1)/2)
        end if
        num_edge_cells = (NumberXHalfNodes - numEvenCells)
        num_edge_real = real(num_edge_cells, kind = 8)
        top = 2.0d0 * pi * (del_x/edge_size - 1.0d0 / num_edge_real)
        bottom = sin(pi + 2.0d0 * pi / num_edge_real)
        del_x_min = (1.0d0 - top/bottom) * edge_size / num_edge_real

        do i = 2,num_edge_cells/2
            phase = real(i-1, kind = 8)/num_edge_real
            self % grid(i) = edge_size* (phase + (1.0d0 - num_edge_real * del_x_min/edge_size)* SIN(2.0d0 * pi * phase)/ 2.0d0 / pi )
            self% grid(NumberXNodes + 1 - i) = self%grid(NumberXNodes) - self%grid(i)
        end do
    end subroutine construct_uniformCenter_InvSineGrid

    subroutine constructHalfSineGrid(self, del_x, L_domain)
        ! First sine grid but half of the grid
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

    subroutine construct_tanh_grid(self, del_x, L_domain, numEvenCells, delta_xi_1)
        ! Thanks to Denis Eremin for providing this mapping
        class(Domain), intent(in out) :: self
        real(real64), intent(in) :: del_x, L_domain, delta_xi_1
        integer(int32), intent(in) :: numEvenCells
        real(real64) :: delta_xi_2, xi_pivot, edge_del_x, xi
        real(real64) :: top, bottom, norm
        integer(int32) :: i, mid_xi

        if (MOD(NumberXHalfNodes,2) == 0) then
            mid_xi = NumberXHalfNodes/2 + 1
        else
            mid_xi = (NumberXNodes)/2
        end if
        edge_del_x = del_x / real(numEvenCells, kind = 8)
        xi_pivot = real(numEvenCells, kind = 8) / real(NumberXHalfNodes, kind = 8)
        bottom = 1.0d0 - real(NumberXHalfNodes, kind = 8) * edge_del_x / L_domain
        top = 1.0d0 + 2.0d0 * delta_xi_1 * log( cosh((0.5d0 - xi_pivot) / delta_xi_1) / cosh( xi_pivot / delta_xi_1 ) )
        delta_xi_2 = edge_del_x * real(NumberXHalfNodes, kind = 8) * top / L_domain / bottom
        
        self%grid(1) = 0.0d0
        self%grid(NumberXNodes) = L_domain
        norm = 1.0d0 + delta_xi_2 + 2.0d0 * delta_xi_1 * log( cosh((0.5d0 - xi_pivot) / delta_xi_1) / cosh( xi_pivot / delta_xi_1 ) )
        do i = 2, mid_xi
            xi = real(i-1, kind = 8)
            top = (1.0d0 + delta_xi_2) * xi / real(NumberXHalfNodes, kind = 8)
            bottom = delta_xi_1 * log( cosh((xi/real(NumberXHalfNodes, kind = 8) - xi_pivot) / delta_xi_1) / cosh( xi_pivot / delta_xi_1 ) )
            self%grid(i) = L_domain * (top + bottom) / norm
            self%grid(NumberXNodes + 1 - i) = L_domain - self%grid(i)
        end do


    end subroutine construct_tanh_grid

    subroutine constructExpHalfGrid(self, del_x, L_domain)
        ! Exponential increase of grid
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
        ! uniform grid
        class(Domain), intent(in out) :: self
        real(real64), intent(in) :: L_domain
        integer(int32) :: i
        self%grid(1) = 0.0d0
        self%grid(NumberXNodes) = L_domain
        do i = 2, NumberXHalfNodes
            self % grid(i) =  (i-1) * L_domain / (NumberXHalfNodes)
        end do
    end subroutine constructUniformGrid

    subroutine constructHalfEvenHalfSinusoid(self, del_x, L_domain, numEvenCells)
        ! Half of grid is even, half increases like sinusoid
        class(Domain), intent(in out) :: self
        real(real64), intent(in) :: del_x, L_domain
        integer(int32), intent(in) :: numEvenCells
        integer(int32) :: numSinCells, i
        real(real64) :: leftX, rightX, evenDelX, middleDelX, num_sin_cells_real
        numSinCells = NumberXNodes - 1 - 2 * numEvenCells
        evenDelX = del_x/real(numEvenCells, kind = 8)
        self%grid(1) = 0.0d0
        self%grid(NumberXNodes) = L_domain
        do i = 1, numEvenCells
            self%grid(i+1) = self%grid(i) + evenDelX
            self%grid(NumberXNodes-i) = self%grid(NumberXNodes-i+1) - evenDelX
        end do
        middleDelX = self%grid(NumberXNodes-numEvenCells) - self%grid(numEvenCells+1)
        num_sin_cells_real = real(numSinCells, kind = 8)
        do i = 2,numSinCells
            self % grid(i+numEvenCells) = self%grid(numEvenCells+1) + middleDelX * ((real(i, kind = 8)-1.0d0)/(num_sin_cells_real) - (1.0d0/(num_sin_cells_real) - evenDelX/middleDelX) &
            * SIN(2.0d0 * pi * (real(i, kind = 8)-1.0d0) / num_sin_cells_real) / SIN(2.0d0 * pi / num_sin_cells_real) )
        end do

    end subroutine constructHalfEvenHalfSinusoid

    subroutine construct_even_grid_doubling_region(self, min_del_x, L_domain, numEvenEdgeCells, max_factor)
        ! Half of grid cells used for even cells on either side up to del_x
        ! Other half used to center, have increasing cell size until cell are at maximum doubled
        class(Domain), intent(in out) :: self
        real(real64), intent(in) :: min_del_x, L_domain, max_factor
        integer(int32), intent(in) :: numEvenEdgeCells
        integer(int32) :: num_max_cells, i, factor_power
        real(real64) :: middle_L, edge_L, num_max_cell_real, max_del_x, doubled_total_L, temp_real
        edge_L = min_del_x * real(numEvenEdgeCells, kind = 8)
        middle_L = L_domain - 2.0d0 * edge_L
        factor_power = NINT(log(max_factor)/log(2.0d0))
        max_del_x = max_factor * min_del_x

        ! Just take from file, this is for benchmark and no need to take more time to get algorithm to create the map
        open(10,file='../../SharedModules/Analysis/ec-pic_ccp_data_for_nicolas/ccp_30cm_30mTorr_1000V_ec_refined32_ion_density.dat')
        do i = 1, NumberXNodes
            read(10, *) temp_real
            self%grid(i) = temp_real
        end do
        close(10)

        if (abs(self%grid(NumberXNodes) - L_domain) / L_domain > 1.d-6) then
            print *, 'For doubling grid for comparison to Tasman paper, incorrect grid end!'
            stop
        end if
       

        
        ! numMiddleCells = NumberXNodes-1 - 2 * numEvenEdgeCells
        ! edge_del_x = del_x/numEvenEdgeCells
        ! print *, edge_del_x
        ! middle_L = L_domain - 2.0d0 * del_x
        ! print *, middle_L
        ! stop
        ! numberIncreaseCells = 0
        ! middle_del_x = middle_L/numMiddleCells
        ! self%grid = 0.0d0
        ! do while (middle_del_x > (2**numberIncreaseCells) * edge_del_x) 
        !     numberIncreaseCells = numberIncreaseCells + 1
        !     center_even_L = middle_L - 2.0d0 * (2**(numberIncreaseCells+1)-2) * edge_del_x
        !     num_center_even_cells = numMiddleCells - 2 * numberIncreaseCells
        !     middle_del_x = center_even_L/num_center_even_cells
        ! end do
        ! numberIncreaseCells = numberIncreaseCells -1
        ! center_even_L = middle_L - 2.0d0 * (2**(numberIncreaseCells+1)-2) * edge_del_x
        ! num_center_even_cells = numMiddleCells - 2 * numberIncreaseCells
        ! middle_del_x = center_even_L/num_center_even_cells
        
        ! self%grid(1) = 0.0d0
        ! self%grid(NumberXNodes) = L_domain
        ! do i = 1, numEvenEdgeCells
        !     self%grid(i+1) = self%grid(i) + edge_del_x
        !     self%grid(NumberXNodes-i) = self%grid(NumberXNodes-i+1) - edge_del_x
        ! end do

        ! do i = 1, numberIncreaseCells
        !     self%grid(numEvenEdgeCells + i + 1) = self%grid(numEvenEdgeCells + i) + (2**i) * edge_del_x
        !     self%grid(NumberXNodes-numEvenEdgeCells - i) = self%grid(NumberXNodes-numEvenEdgeCells - i + 1) - (2**i) * edge_del_x
        ! end do

        ! do i = numEvenEdgeCells + numberIncreaseCells + 2, NumberXNodes - numEvenEdgeCells - numberIncreaseCells - 1
        !     self%grid(i) = self%grid(i-1) + middle_del_x
        ! end do

    end subroutine construct_even_grid_doubling_region

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

    function getXFromL(self, l) result(x)
        ! get physical coordinate from logical
        class(Domain), intent(in) :: self
        real(real64), intent(in) :: l
        real(real64) :: x
        x = self%grid(INT(l)) + self%dx_dl(INT(l)) * (l - INT(l))

    end function getXFromL

    function getLFromX(self, x) result(l)
        ! get logical from physical space, use binary search
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
        real(real64) :: debyeLength, L_domain, extra_real
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
        read(10, *, IOSTAT = io) gridType, debyeLength, extra_int, extra_real
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
        world = Domain(leftBoundary, rightBoundary)
        call world % constructGrid(debyeLength, L_domain, gridType, extra_int, extra_real)
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