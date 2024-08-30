module mod_domain
    use iso_fortran_env, only: int32, real64
    use constants
    use mod_BasicFunctions
    implicit none

    private
    public :: Domain, readWorld
    integer(int32), public, protected :: NumberXNodes = 10 ! Number of nodes, where phi and E defined
    integer(int32), public, protected :: NumberXHalfNodes = 10 ! Number of cells

    ! domain contains arrays and values related to physical, logical dimensions of the spatial grid
    type :: Domain
        real(real64), allocatable :: grid(:) !Physical grid where phi is evaluated
        real(real64) :: delX, L_domain, startX, endX ! Physical end points, cell size, domain size
        integer(int32), allocatable :: boundaryConditions(:), threadNodeIndx(:,:) ! Boundary condition flags for fields and particles
        integer(int32) :: numThreadNodeIndx
        logical :: gridSmoothBool
        ! (>0 dirichlet, -2 Neumann, -3 periodic, <=-4 dielectric), 0 is default in-body condition 

    contains
        procedure, public, pass(self) :: constructHalfSineGrid
        procedure, public, pass(self) :: constructSineGrid
        procedure, public, pass(self) :: constructInvSineGrid
        procedure, public, pass(self) :: constructUniformGrid
        procedure, public, pass(self) :: constructExpHalfGrid
        procedure, public, pass(self) :: constructHalfEvenHalfSinusoid

        procedure, public, pass(self) :: smoothDensities
        procedure, public, pass(self) :: smoothField
        procedure, public, pass(self) :: smoothRho
        procedure, public, pass(self) :: constructGrid
        procedure, public, pass(self) :: writeDomain
        procedure, public, pass(self) :: getLFromX
        procedure, public, pass(self) :: getXFromL
        procedure, public, pass(self) :: addThreadedDomainArray
    end type Domain


    interface Domain
        module procedure :: domain_constructor
    end interface Domain

contains

    type(Domain) function domain_constructor(leftBoundary, rightBoundary) result(self)
        ! Construct domain object, initialize grid, dx_dl
        integer(int32), intent(in) :: leftBoundary, rightBoundary
        integer(int32) :: i, k, spacingThread, modThread
        allocate(self % grid(NumberXNodes), self%boundaryConditions(NumberXNodes))
        self % grid = (/(i, i=0, NumberXNodes-1)/)
        self % delX = 0.0d0
        self % L_domain = 0.0d0
        self % boundaryConditions = 0
        self%boundaryConditions(1) = leftBoundary
        self%boundaryConditions(NumberXNodes) = rightBoundary
        self%gridSmoothBool = .true.
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
    end function domain_constructor

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

    subroutine constructInvSineGrid(self, del_x, L_domain)
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



    subroutine smoothDensities(self, densities)
        ! Smooth Rho on grid with quadratic smoothing
        class(Domain), intent(in) :: self
        real(real64), intent(in out) :: densities(NumberXHalfNodes)
        real(real64) :: work(NumberXNodes)

        work(2:NumberXHalfNodes) = 0.25d0 * (densities(1:NumberXHalfNodes-1) + 2.0d0 * densities(2:NumberXHalfNodes) + densities(3:NumberXNodes))
        SELECT CASE(self%boundaryConditions(1))
            CASE(1,2,4)
                work(1) = 0.25d0 * (2.0d0 * densities(1) + densities(2))
                work(2) = work(2) + 0.25d0 * densities(1)
            CASE(3)
                work(1) = 0.25d0 * (densities(NumberXHalfNodes) + 2.0d0 * densities(1) + densities(2))
        END SELECT

        SELECT CASE(self%boundaryConditions(NumberXNodes))
            CASE(1,2,4)
                work(NumberXNodes) = 0.25d0 * (2.0d0 * densities(NumberXNodes) + densities(NumberXHalfNodes))
                work(NumberXHalfNodes) = work(NumberXHalfNodes) + 0.25d0 * densities(NumberXNodes)
            CASE(3)
                work(NumberXNodes) = 0.25d0 * (densities(NumberXHalfNodes) + 2.0d0 * densities(1) + densities(2))
        END SELECT
        densities = work
    end subroutine smoothDensities    

    subroutine smoothRho(self, Rho)
        ! Smooth Rho on grid with quadratic smoothing
        class(Domain), intent(in) :: self
        real(real64), intent(in out) :: Rho(NumberXHalfNodes)
        real(real64) :: work(NumberXNodes)   
    
        work(2:NumberXHalfNodes) = 0.25d0 * (Rho(1:NumberXHalfNodes-1) + 2.0d0 * Rho(2:NumberXHalfNodes) + Rho(3:NumberXHalfNodes+1))
        
        SELECT CASE(self%boundaryConditions(1))
            CASE(1,4)
                work(1) = 0.0d0
            CASE(2)
                work(1) = 0.25d0 * (2.0d0 * Rho(1) + 2.0d0 * Rho(2))
            CASE(3)
                work(1) = 0.25d0 * (Rho(NumberXHalfNodes) + 2.0d0 * Rho(1) + Rho(2))
        END SELECT

        SELECT CASE(self%boundaryConditions(NumberXNodes))
            CASE(1,4)
                work(NumberXNodes) = 0.0d0
            CASE(2)
                work(NumberXNodes) = 0.25d0 * (2.0d0 * Rho(NumberXNodes) + 2.0d0 * Rho(NumberXHalfNodes))
            CASE(3)
                work(NumberXNodes) = 0.25d0 * (Rho(NumberXHalfNodes) + 2.0d0 * Rho(NumberXNodes) + Rho(2))
        END SELECT
    end subroutine smoothRho


    subroutine smoothField(self, rawField, newField)
        ! Smooth fields on grid with quadratic smoothing
        class(Domain), intent(in) :: self
        real(real64), intent(in) :: rawField(NumberXHalfNodes)
        real(real64), intent(in out) :: newField(NumberXHalfNodes)

        newField(2:NumberXHalfNodes-1) = 0.25d0 * (rawField(1:NumberXHalfNodes-2) + 2.0d0 * rawField(2:NumberXHalfNodes-1) + rawField(3:NumberXHalfNodes))
        
        SELECT CASE(self%boundaryConditions(1))
            CASE(1,4)
                newField(1) = 0.25d0 * (3.0d0 * rawField(1) + rawField(2))
            CASE(2)
                newField(1) = 0.25d0 * (rawField(1) + rawField(2))
            CASE(3)
                newField(1) = 0.25d0 * (rawField(NumberXHalfNodes) + 2.0d0 * rawField(1) + rawField(2))
        END SELECT

        SELECT CASE(self%boundaryConditions(NumberXNodes))
            CASE(1,4)
                newField(NumberXHalfNodes) = 0.25d0 * (3.0d0 * rawField(NumberXHalfNodes) + rawField(NumberXHalfNodes-1))
            CASE(2)
                newField(NumberXHalfNodes) = 0.25d0 * (rawField(NumberXHalfNodes) + rawField(NumberXHalfNodes-1))
            CASE(3)
                newField(NumberXHalfNodes) = 0.25d0 * (rawField(NumberXHalfNodes-1) + 2.0d0 * rawField(NumberXHalfNodes) + rawField(1))
        END SELECT

    end subroutine smoothField

    subroutine constructGrid(self, L_domain)
        ! Make grid
        class(Domain), intent(in out) :: self
        real(real64), intent(in) :: L_domain
        self%L_domain = L_domain
        self%delX = self%L_domain/(NumberXNodes - 1)
        self%grid = self%grid * self%delX
        self%startX = self%grid(1)
        self%endX = self%grid(NumberXNodes)
    end subroutine constructGrid

    function getLFromX(self, x) result(l)
        ! get computational location from physical location
        class(Domain), intent(in) :: self
        real(real64), intent(in) :: x
        real(real64) :: l
        l = x/self%delX + 1.0d0
        
    end function getLFromX

    function getXFromL(self, l) result(x)
        ! get physical location from computational
        class(Domain), intent(in) :: self
        real(real64), intent(in) :: l
        real(real64) :: x
        x = (l-1) * self%delX + self%startX
        
    end function getXFromL

    subroutine addThreadedDomainArray(self, array_add, x, iThread, const)
        ! Take array on grid nodes of half nodes x with second dimension thread count and add to array_add of same domain dimension using Openmp
        class(Domain), intent(in) :: self
        real(real64), intent(in out) :: array_add(NumberXNodes)
        real(real64), intent(in) :: x(NumberXNodes, numThread), const
        integer(int32), intent(in) :: iThread
        
        array_add(self%threadNodeIndx(1,iThread):self%threadNodeIndx(2,iThread)) = array_add(self%threadNodeIndx(1,iThread):self%threadNodeIndx(2,iThread)) &
            + SUM(x(self%threadNodeIndx(1,iThread):self%threadNodeIndx(2,iThread), :), DIM=2) * const

    end subroutine addThreadedDomainArray

    subroutine writeDomain(self, dirName)
        ! Writes domain data into binary file under Data
        class(Domain), intent(in) :: self
        character(*), intent(in) :: dirName
        open(41,file=dirName//"/domainGrid.dat", form='UNFORMATTED')
        write(41) self%grid
        close(41)
        open(41,file=dirName//"/domainBoundaryConditions.dat", form='UNFORMATTED')
        write(41) self%boundaryConditions
        close(41)
    end subroutine writeDomain


    ! -------------- read in world components -----------------------------------


    subroutine readWorld(GeomFilename, world, T_e, n_ave)
        ! Read inputs for domain
        type(Domain), intent(in out) :: world
        character(len=*), intent(in) :: GeomFilename
        real(real64), intent(in) :: T_e, n_ave
        integer(int32) :: io, leftBoundary, rightBoundary
        real(real64) :: debyeLength, L_domain, fracDebye

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
        read(10, *, IOSTAT = io) leftBoundary, rightBoundary
        read(10, *, IOSTAT = io)
        read(10, *, IOSTAT = io) fracDebye
        read(10, *, IOSTAT = io) 
        close(10)
        debyeLength = getDebyeLength(T_e, n_ave)
        if (L_domain / (NumberXNodes-1) > fracDebye * debyeLength) then
            print *, "Insufficient amount of nodes to resolve initial debyeLength"
            NumberXNodes = NINT(L_domain/debyeLength/fracDebye) + 1
        end if
        NumberXHalfNodes = NumberXNodes-1
        if ((leftBoundary == 3) .or. (rightBoundary == 3)) then
            leftBoundary = 3
            rightBoundary = 3
        end if
        world = Domain(leftBoundary, rightBoundary)
        call world % constructGrid(L_domain) 
        print *, "Number of nodes:", NumberXNodes
        print *, "Left boundary type:", world%boundaryConditions(1)
        print *, "Right boundary type:", world%boundaryConditions(NumberXNodes)
        print *, 'Fraction of debye length:', world%delX/getDebyeLength(T_e, n_ave)
        print *, "Grid length:", world%L_domain
        print *, 'DelX is:', world%delX
        print *, "------------------"
        print *, ""

    end subroutine readWorld


end module mod_domain