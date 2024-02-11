module mod_PCG
    use, intrinsic :: iso_fortran_env
    implicit none
    
    public :: pcg_init, getPCGRand


    integer(int64), parameter :: seed = 5573589319906701683 ! Seed-dependent initial state
    integer(int64), parameter :: multiplier = 6364136223846793005
    integer(int64), parameter :: increment = 1442695040888963407
    integer(int64), parameter :: max_int32 = 2147483648
    real(real64), parameter :: norm_PCG = 4294967296.0d0

contains

    function getPCGRand(state) result(result)
        integer(int64), intent(inout) :: state
        integer(int64) :: x
        integer(int32) :: count, x_shifted, random_number
        real(real64) :: result

        x = state
        count = int(shiftR(x,59), kind = int32)

        state = x * multiplier + increment
        x = ieor(x, shiftR(x, 18))
        x_shifted = int(shiftR(x, 27), kind = int32)
        random_number = ior(shiftR(x_shifted, count), shiftL(x_shifted, iand(-count, 31)))
        result = real(random_number, kind = real64) + max_int32
        result = result / norm_PCG
    end function getPCGRand

    subroutine pcg_init(state)
        integer(int64), intent(in out) :: state

        state = state + increment
        state = state * multiplier + increment
    end subroutine pcg_init

end module mod_PCG