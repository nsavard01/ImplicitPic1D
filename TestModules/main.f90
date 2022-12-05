program main
    use iso_fortran_env, only: int32, real64
    use omp_lib
    use mod_test
    implicit none

    integer(int32) :: nproc, i, i_omp, tclock1, tclock2, clock_rate, n = 100
    real(real64), allocatable :: x(:, :)
    real(real64) :: m = 10.0, elapsed_time
    type(type1), allocatable :: object(:)

    nproc = omp_get_max_threads()
    allocate(x(n, nproc))
    allocate(object(nproc))
    do i = 1, nproc
        x(:, i) = i
        object(i) = type1(m, n)
        object(i)%x = i
    end do

    call system_clock(tclock1)
    do i =1, 10000
        !$OMP PARALLEL PRIVATE(i_omp)
        i_omp = omp_get_thread_num() + 1
        call testFunc(x, m, i_omp)
        !$OMP END PARALLEL
    end do
    call system_clock(tclock2, clock_rate)
    elapsed_time = float(tclock2 - tclock1) / float(clock_rate)
    print *, "Elapsed time for functional is:", elapsed_time, "seconds"
    print *, x(1:10, 1)
    print *, x(1:10, nproc)

    call system_clock(tclock1)
    do i =1, 10000
        !$OMP PARALLEL PRIVATE(i_omp)
        i_omp = omp_get_thread_num() + 1
        call object(i_omp)%typeFunc()
        !$OMP END PARALLEL
    end do
    call system_clock(tclock2, clock_rate)
    elapsed_time = float(tclock2 - tclock1) / float(clock_rate)
    print *, "Elapsed time for OOP is:", elapsed_time, "seconds"
    print *, object(1)%x(1:10)
    print *, object(nproc)%x(1:10)

contains

    subroutine testFunc(x, m, i_omp)
        real(real64), intent(in out) :: x(:, :)
        real(real64), intent(in) :: m
        integer(int32), intent(in) :: i_omp
        x(:, i_omp) = x(:, i_omp) + m

    end subroutine testFunc

end program main