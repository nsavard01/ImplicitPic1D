module mod_test
    use iso_fortran_env, only: int32, real64

    real(real64), allocatable :: h(:,:)

contains

    subroutine initializeH(h, n_x, n_y)
        real(real64), allocatable, intent(in out) :: h(:,:)
        integer(int32), intent(in) :: n_x, n_y
        allocate(h(n_x, n_y))
        h = 0
    end subroutine initializeH

    subroutine add10(h)
        real(real64), intent(in out) :: h(:,:)
        h = h + 10
    end subroutine

end module mod_test