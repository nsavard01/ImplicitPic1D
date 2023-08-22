module mod_multiGridStorage
    use iso_fortran_env, only: int32, real64
    implicit none
    ! temporary class for system A*x = b
    private
    public :: multiGridStorage

    type :: multiGridStorage
        integer(int32) :: N_x
        real(real64), allocatable :: phi_C(:), R_C(:), res(:) !phi_f is final phi, will likely need to store two arrays for phi, can't be avoided
        real(real64), allocatable :: a_tri(:), b_tri(:), c_tri(:) !for thomas algorithm potential solver, a_tri is lower diagonal, b_tri middle, c_tri upper


    contains
        procedure, public, pass(self) :: constructPoissonEven
    end type

    interface multiGridStorage
        module procedure :: multiGridStorage_constructor
    end interface multiGridStorage
contains

    type(multiGridStorage) function multiGridStorage_constructor(N_x) result(self)
        ! Construct domain object, initialize grid, dx_dl, and dx_dl.
        integer(int32), intent(in) :: N_x
        allocate(self % phi_C(N_x), self%R_C(N_x), self%a_tri(N_x - 1), &
        self%b_tri(N_x), self%c_tri(N_x-1), self%res(N_x))
        self % a_tri = 0.0d0
        self % c_tri = 0.0d0
        self % b_tri = 0.0d0
        self % res = 0.0d0
        self % phi_C = 0.0d0
        self % R_C = 0.0d0
        self % N_x = N_x

    end function multiGridStorage_constructor

    subroutine constructPoissonEven(self, dx)
        class(multiGridStorage), intent(in out) :: self
        real(real64), intent(in) :: dx
        integer(int32) :: i
        self%c_tri(1) = 0.0d0
        self%b_tri(1) = 1.0d0
        do i = 2, self%N_x-1
            self%c_tri(i) = 1.0d0/dx**2
            self%a_tri(i-1) = 1.0d0/dx**2
            self%b_tri(i) = - 2.0d0/dx**2
        end do
        self%a_tri(self%N_x-1) = 2.0d0/dx**2
        self%b_tri(self%N_x) = -2.0d0/dx**2
    end subroutine constructPoissonEven

end module mod_multiGridStorage