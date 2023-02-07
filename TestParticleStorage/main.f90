program main
    use iso_fortran_env, only: int32, real64
    use omp_lib
    use mod_test
    implicit none

    integer(int32) :: nproc, i, j, i_omp, tclock1, tclock2, clock_rate, numPartPerThread
    real(real64), allocatable :: x(:, :, :, :), R(:, :)
    real(real64) :: m(50), a, elapsed_time
    type(type1) :: object(50)
    numPartPerThread = 10000
    nproc = omp_get_max_threads() 
    print *, "Total amount of threads is:", nproc
    allocate(x(6, numPartPerThread, nproc, 50), R(6, numPartPerThread))
    a = 1458.98
    do j = 1, 50
        m(j) = j
        object(j) = type1(real(j, real64), numPartPerThread, nproc)
        do i = 1, nproc
            call random_number(R)
            x(:, :, i, j) = R
            object(j)%x(:, :, i) = R
        end do
    end do

    print *, m(25), object(25)%m
    print *, x(1:3, 4000, 40, 25)
    print *, object(25)%x(1:3, 4000, 40)

    call system_clock(tclock1)
    do i =1, 10
        !$OMP PARALLEL PRIVATE(i_omp)
        i_omp = omp_get_thread_num() + 1
        do j = 1, 50
            call testFuncObject(object(j), i_omp, numPartPerThread, a)
        end do
        !$OMP END PARALLEL
    end do
    call system_clock(tclock2, clock_rate)
    elapsed_time = float(tclock2 - tclock1) / float(clock_rate)
    print *, "Elapsed time for object is:", elapsed_time, "seconds"
    print *, object(25)%x(1:3, 4000, 40)

    call system_clock(tclock1)
    do i =1, 10
        !$OMP PARALLEL PRIVATE(i_omp)
        i_omp = omp_get_thread_num() + 1
        do j = 1, 50
            call testFunc(x(:,:,:,j), m(j), i_omp, numPartPerThread, a)
        end do
        !$OMP END PARALLEL
    end do
    call system_clock(tclock2, clock_rate)
    elapsed_time = float(tclock2 - tclock1) / float(clock_rate)
    print *, "Elapsed time for functional is:", elapsed_time, "seconds"
    print *, x(1:3, 4000, 40, 25)

    

contains

    subroutine testFunc(x, m, i_omp, numPartPerThread, a)
        real(real64), intent(in out) :: x(:, :, :)
        real(real64), intent(in) :: m, a
        integer(int32), intent(in) :: i_omp, numPartPerThread
        integer(int32) :: i
        do i=1, numPartPerThread
            x(1:3, i, i_omp) = x(1:3, i, i_omp) + x(4:6, i, i_omp)*a * m
        end do

    end subroutine testFunc

    subroutine testFuncObject(object, i_omp, numPartPerThread, a)
        type(type1), intent(in out) :: object
        real(real64), intent(in) :: a
        integer(int32), intent(in) :: i_omp, numPartPerThread
        integer(int32) :: i
        do i=1, numPartPerThread
            object%x(1:3, i, i_omp) = object%x(1:3, i, i_omp) + object%x(4:6, i, i_omp)*a * object%m
        end do

    end subroutine testFuncObject

end program main