module mod_domain
    use iso_fortran_env, only: int32, real64
    use constants
    implicit none

    private
    public :: Domain

    ! domain contains arrays and values related to physical, logical dimensions of the spatial grid
    type :: Domain
        real(real64), allocatable :: grid(:) !Physical grid where phi is evaluated
        real(real64) :: delX !Physical grid where E-Fields are
        integer(int32), allocatable :: boundaryConditions(:) ! Boundary condition flags for fields and particles
        ! (>0 dirichlet, -2 Neumann, -3 periodic, <=-4 dielectric), 0 is default in-body condition 

    contains
        procedure, public, pass(self) :: constructGrid
        procedure, public, pass(self) :: writeDomain
    end type Domain


    interface Domain
        module procedure :: domain_constructor
    end interface Domain

contains

    ! Initialization procedures
    type(Domain) function domain_constructor(leftBoundary, rightBoundary) result(self)
        ! Construct domain object, initialize grid, dx_dl, and dx_dl.
        integer(int32), intent(in) :: leftBoundary, rightBoundary
        integer(int32) :: i
        allocate(self % grid(NumberXNodes), self%boundaryConditions(NumberXNodes+1))
        self % grid = (/(i, i=0, NumberXNodes-1)/)
        self % delX = 0.0d0
        self % boundaryConditions = 0
        self%boundaryConditions(1) = leftBoundary
        self%boundaryConditions(NumberXNodes) = rightBoundary
    end function domain_constructor


    subroutine constructGrid(self, L_domain)
        class(Domain), intent(in out) :: self
        real(real64), intent(in) :: L_domain
        self%delX = L_domain/(NumberXNodes - 1)
        self%grid = self%grid * self%delX
    end subroutine constructGrid

   

    subroutine writeDomain(self)
        ! Writes domain data into binary file under Data
        class(Domain), intent(in) :: self
        open(41,file="../Data/domainGrid.dat", form='UNFORMATTED')
        write(41) self%grid
        close(41)
    end subroutine writeDomain





end module mod_domain