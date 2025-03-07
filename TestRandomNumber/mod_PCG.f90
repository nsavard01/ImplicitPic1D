module mod_PCG
    use, intrinsic :: iso_fortran_env
    use iso_c_binding
    implicit none

    type,public :: PCG_type

    !! main class for random number generator

    private

    integer(c_int64_t) :: state

  contains

    private

    
    procedure, public :: initialize

    procedure, public :: get_rand

  end type PCG_type


    interface 
        real(c_double) function pcg32_random_r(state) bind(c)
        use iso_c_binding
        integer(c_int64_t) :: state
        end function
    end interface

contains

    subroutine initialize(self,seed)
    !! Initialize with seed

    class(PCG_type),intent(inout) :: self
    integer(c_int64_t), intent(in) :: seed
    self%state = seed

  end subroutine initialize

    function get_rand(self) result(res)
        class(PCG_type), intent(inout) :: self
        real(real64) :: res
        res = pcg32_random_r(self%state)
    end function get_rand

end module mod_PCG