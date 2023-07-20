program main
    use iso_fortran_env, only: int32, real64
    use mod_test
    implicit none

    integer(int32) :: nproc = 100, i, tclock1, tclock2, clock_rate, n = 5000
    real(real64), allocatable :: x(:, :)
    real(real64) :: m = 10.0, elapsed_time
    type(type1), allocatable :: object(:)

    allocate(x(n, nproc))
    allocate(object(nproc))
    do i = 1, nproc
        x(:, i) = i
        object(i) = type1(m, n)
        object(i)%x = i
    end do

    open(41,file='dataObject.dat', form='UNFORMATTED')
    call system_clock(tclock1)
    do i=1, nproc
        write(41) object(i)%x
    end do
    call system_clock(tclock2, clock_rate)
    close(41)
    elapsed_time = float(tclock2 - tclock1) / float(clock_rate)
    print *, "Elapsed time for OOP is:", elapsed_time, "seconds"
    print *, object(1)%x(1:10)
    print *, object(nproc)%x(1:10)


    open(41,file='dataNormal.dat',form='UNFORMATTED')
    call system_clock(tclock1)
    write(41) x
    call system_clock(tclock2, clock_rate)
    close(41)
    elapsed_time = float(tclock2 - tclock1) / float(clock_rate)
    print *, "Elapsed time for functional is:", elapsed_time, "seconds"
    print *, x(1:10, 1)
    print *, x(1:10, nproc)

    
contains

    subroutine testFunc(x, m, i_omp)
        real(real64), intent(in out) :: x(:, :)
        real(real64), intent(in) :: m
        integer(int32), intent(in) :: i_omp
        x(:, i_omp) = x(:, i_omp) + m

    end subroutine testFunc

end program main