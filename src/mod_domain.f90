module mod_domain
    use iso_fortran_env, only: int32, real64
    use constants
    use mod_particle
    implicit none

    private
    public :: Domain

    ! domain contains arrays and values related to physical, logical dimensions of the spatial grid
    type :: Domain
        integer(int32) :: n_x !number grid nodes
        real(real64), allocatable :: grid(:) !array representing values at grid in m
        real(real64), allocatable :: dx_dl(:) !ratio of grid differences from physical to logical, assume logical separated by 1
        real(real64), allocatable :: nodeVol(:) !node vol, or difference between half grid in physical space
        real(real64), allocatable :: phi(:), rho(:), J(:) ! array of phi, rho, and J values in space

    contains
        procedure, private, pass(self) :: derive_DxDl_NodeVol
        procedure, public, pass(self) :: constructSineGrid
        procedure, public, pass(self) :: constructUniformGrid
        procedure, public, pass(self) :: depositRho
    end type Domain


    interface Domain
        module procedure :: domain_constructor
    end interface Domain

contains

    ! Initialization procedures
    type(Domain) function domain_constructor(num) result(self)
        ! Construct domain object, initialize grid, dx_dl, and nodeVol.
        integer(int32), intent(in) :: num
        integer(int32) :: i
        self % n_x = num
        allocate(self % grid(num), self % dx_dl(num-1), self % nodeVol(num) &
        , self % phi(num), self % rho(num), self % J(num-1))
        self % grid = (/(i, i=1, num)/)
        self % dx_dl = 1
        self % nodeVol = 1
        self % phi = 0
        self % rho = 0
        self % J = 0
    end function domain_constructor

    pure subroutine derive_DxDl_NodeVol(self)
        class(Domain), intent(in out) :: self
        integer(int32) :: i
        do i = 1, self%n_x-1
            self%dx_dl(i) = self%grid(i+1) - self%grid(i)
        end do

        self%nodeVol(1) = self%dx_dl(1)/2
        self%nodeVol(self%n_x) = self%dx_dl(self%n_x-1)/2

        do i = 2, self%n_x-1
            self%nodeVol(i) = (self%dx_dl(i-1) + self%dx_dl(i))/2
        end do
    end subroutine derive_DxDl_NodeVol

    pure subroutine constructSineGrid(self, del_l, L_domain)
        class(Domain), intent(in out) :: self
        real(real64), intent(in) :: del_l, L_domain
        integer(int32) :: i
        self%grid(1) = 0
        self%grid(self%n_x) = L_domain
        do concurrent (i = 2:self % n_x-1)
            self % grid(i) = L_domain * ((i-1) - (self%n_x - 1)*(1 - del_l)/pi/2 &
            * SIN(2 * pi * (i-1) / (self%n_x - 1))) / (self%n_x - 1)
        end do
        call derive_DxDl_NodeVol(self)
    end subroutine constructSineGrid

    pure subroutine constructUniformGrid(self, L_domain)
        class(Domain), intent(in out) :: self
        real(real64), intent(in) :: L_domain
        integer(int32) :: i
        self%grid(1) = 0
        self%grid(self%n_x) = L_domain
        do concurrent (i = 2:self % n_x-1)
            self % grid(i) =  (i-1) * L_domain / (self%n_x - 1)
        end do
        call derive_DxDl_NodeVol(self)
    end subroutine constructUniformGrid

    subroutine depositRho(self, particleList) 
        class(Domain), intent(in out) :: self
        type(Particle), intent(in) :: particleList(1)
        integer(int32) :: i, j, l_left
        real(real64) :: d
        do i=1, size(particleList)
            do j = 1, particleList(i)%N_p
                l_left = INT(particleList(i)%l_p(j))
                d = MOD(particleList(i)%l_p(j), 1.0)
                self % rho(l_left) = self % rho(l_left) + particleList(i)%q * particleList(i)%w_p * (1.0-d)
                self % rho(l_left + 1) = self % rho(l_left + 1) + particleList(i)%q * particleList(i)%w_p * d
            end do
        end do
        self % rho = self % rho / self % nodeVol
    end subroutine depositRho





end module mod_domain