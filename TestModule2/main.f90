program main
    use iso_fortran_env, only: int32, real64
    use omp_lib
    use mod_test
    implicit none

    integer(int32) :: nproc, i, i_omp, sizek = 1000000,tclock1, tclock2, clock_rate
    real(real64), allocatable :: a(:,:), b(:,:)
    real(real64) :: elapsed_time
    type(type1) :: type_a, type_b
    nproc = omp_get_max_threads()
    type_a = type1(sizek, nproc)
    type_b = type1(sizek-1, nproc)
    
    print *, "Have", nproc, "processors to use."
    allocate(a(sizek,nproc), b(sizek-1, nproc))
    !allocate(object(nproc))
    do i_omp = 1, nproc
        do i = 1, sizek
            a(i, i_omp) = i
            type_a%x(i, i_omp) = i
        end do
    end do
    b = 0.0

    call system_clock(tclock1)
    do i = 1, 1000
        !$OMP PARALLEL PRIVATE(i_omp)
        i_omp = omp_get_thread_num() + 1
        call testFunc(sizek, nproc, type_a%x, type_b%x, i_omp)
        !$OMP END PARALLEL
    end do
    call system_clock(tclock2, clock_rate)
    elapsed_time = float(tclock2 - tclock1) / float(clock_rate)
    print *, "Elapsed time for partial OOP is:", elapsed_time, "seconds"

    call system_clock(tclock1)
    do i = 1, 1000
        !$OMP PARALLEL PRIVATE(i_omp)
        i_omp = omp_get_thread_num() + 1
        call testFuncOOP(sizek,type_a, type_b, i_omp)
        !$OMP END PARALLEL
    end do
    call system_clock(tclock2, clock_rate)
    elapsed_time = float(tclock2 - tclock1) / float(clock_rate)
    print *, "Elapsed time for full OOP is:", elapsed_time, "seconds"

    call system_clock(tclock1)
    do i = 1, 1000
        !$OMP PARALLEL PRIVATE(i_omp)
        i_omp = omp_get_thread_num() + 1
        call testFunc(sizek, nproc, a, b, i_omp)
        !$OMP END PARALLEL
    end do
    call system_clock(tclock2, clock_rate)
    elapsed_time = float(tclock2 - tclock1) / float(clock_rate)
    print *, "Elapsed time for functional is:", elapsed_time, "seconds"

    

    


    

contains

    subroutine testFunc(n1, n2, a, b, i_omp)
        integer(int32), intent(in) :: i_omp, n1, n2
        real(real64), intent(in out) :: b(n1-1, n2)
        real(real64), intent(in) :: a(n1,n2)
        integer(int32) :: i
        b(:, i_omp) = 0.0

        do i = 2, n1
            b(i-1, i_omp) = a(i-1, i_omp) * a(i, i_omp)
        end do

    end subroutine testFunc

    subroutine testFuncOOP(n1, type_a, type_b, i_omp)
        integer(int32), intent(in) :: i_omp, n1
        type(type1), intent(in out) :: type_b
        type(type1), intent(in) :: type_a
        integer(int32) :: i
        type_b%x(:, i_omp) = 0.0

        do i = 2, n1
            type_b%x(i-1, i_omp) = type_a%x(i-1, i_omp) * type_a%x(i, i_omp)
        end do

    end subroutine testFuncOOP

end program main