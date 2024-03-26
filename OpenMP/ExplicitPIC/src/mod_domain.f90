module mod_domain
    use iso_fortran_env, only: int32, real64
    use constants
    use mod_BasicFunctions
    implicit none

    private
    public :: Domain, readWorld
    integer(int32), public, protected :: NumberXNodes = 10
    integer(int32), public, protected :: NumberXHalfNodes = 10

    ! domain contains arrays and values related to physical, logical dimensions of the spatial grid
    type :: Domain
        real(real64), allocatable :: grid(:) !Physical grid where phi is evaluated
        real(real64) :: delX, L_domain, startX, endX !Physical grid where E-Fields are
        integer(int32), allocatable :: boundaryConditions(:), threadNodeIndx(:,:) ! Boundary condition flags for fields and particles
        integer(int32) :: numThreadNodeIndx
        ! (>0 dirichlet, -2 Neumann, -3 periodic, <=-4 dielectric), 0 is default in-body condition 

    contains
        procedure, public, pass(self) :: constructGrid
        procedure, public, pass(self) :: writeDomain
        procedure, public, pass(self) :: getLFromX
        procedure, public, pass(self) :: addThreadedDomainArray
    end type Domain


    interface Domain
        module procedure :: domain_constructor
    end interface Domain

contains

    ! Initialization procedures
    type(Domain) function domain_constructor(leftBoundary, rightBoundary) result(self)
        ! Construct domain object, initialize grid, dx_dl, and dx_dl.
        integer(int32), intent(in) :: leftBoundary, rightBoundary
        integer(int32) :: i, k, spacingThread, modThread
        allocate(self % grid(NumberXNodes), self%boundaryConditions(NumberXNodes))
        self % grid = (/(i, i=0, NumberXNodes-1)/)
        self % delX = 0.0d0
        self % L_domain = 0.0d0
        self % boundaryConditions = 0
        self%boundaryConditions(1) = leftBoundary
        self%boundaryConditions(NumberXNodes) = rightBoundary
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


    subroutine constructGrid(self, L_domain)
        class(Domain), intent(in out) :: self
        real(real64), intent(in) :: L_domain
        self%L_domain = L_domain
        self%delX = self%L_domain/(NumberXNodes - 1)
        self%grid = self%grid * self%delX
        self%startX = self%grid(1)
        self%endX = self%grid(NumberXNodes)
    end subroutine constructGrid

    function getLFromX(self, x) result(l)
        class(Domain), intent(in) :: self
        real(real64), intent(in) :: x
        real(real64) :: l
        l = x/self%delX + 1.0d0
        
    end function getLFromX

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
        type(Domain), intent(in out) :: world
        character(len=*), intent(in) :: GeomFilename
        real(real64), intent(in) :: T_e, n_ave
        integer(int32) :: io, leftBoundary, rightBoundary
        real(real64) :: debyeLength, L_domain, fracDebye

        print *, ""
        print *, "Reading domain inputs:"
        print *, "------------------"
        open(10,file='../InputData/'//GeomFilename)
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