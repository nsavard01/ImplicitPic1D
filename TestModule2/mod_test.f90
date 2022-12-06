module mod_test
    use iso_fortran_env, only: int32, real64
    implicit none

    public :: type1
    type :: type1
        real(real64), allocatable :: x(:, :)

    end type type1

    interface type1
        module procedure :: type1_constructor
    end interface type1

contains

    type(type1) function type1_constructor(m, n) result(self)
        ! Construct particle object, sizeIncrease is fraction larger stored array compared to initial amount of particles
        integer(int32), intent(in) :: m,n
        allocate(self%x(m, n))
        self%x = 0
    end function type1_constructor


end module mod_test