module mod_test
    use iso_fortran_env, only: int32, real64
    implicit none

    type :: type1
        real(real64), allocatable :: x(:, :, :)
        real(real64) :: m

    contains
        procedure, public, pass(self) :: typeFunc

    end type type1

    interface type1
        module procedure :: type1_constructor
    end interface type1

contains

    type(type1) function type1_constructor(m, numberPerThread, nproc) result(self)
        ! Construct particle object, sizeIncrease is fraction larger stored array compared to initial amount of particles
        real(real64), intent(in) :: m
        integer(int32), intent(in) :: numberPerThread, nproc
        allocate(self%x(6, numberPerThread, nproc))
        self%x = 0
        self%m = m
        
    end function type1_constructor

    subroutine typeFunc(self)
        class(type1), intent(in out) :: self
        self%x = self%x + self%m

    end subroutine typeFunc


end module mod_test