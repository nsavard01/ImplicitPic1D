program main
    use iso_fortran_env, only: int32, real64
    use omp_lib
    implicit none

    integer(int32) :: nproc, i, i_omp
    real(real64), allocatable :: x(:, :)
    real(real64) :: m = 10.0

    nproc = omp_get_max_threads()
    allocate(x(4, nproc))
    do i = 1, nproc
        x(:, i) = i
    end do

    !$OMP PARALLEL PRIVATE(i_omp)
    i_omp = omp_get_thread_num() + 1
    call testFunc(x, m, i_omp)
    !$OMP END PARALLEL
    ! do i_omp = 1, nproc
    !     call testFunc(x, m, i_omp)
    ! end do

    

contains

    subroutine testFunc(x, m, i_omp)
        real(real64), intent(in out) :: x(:, :)
        real(real64), intent(in) :: m
        integer(int32), intent(in) :: i_omp
        x(:, i_omp) = x(:, i_omp) + m

    end subroutine testFunc

end program main