

program pardiso_unsym
    use iso_fortran_env, only: int32, real64
    implicit none
    integer(int32) :: mtype, pt(64), solver, iparm(64), mnum = 1, maxfct = 1, nrhs = 1, msglvl = 0, matDimension = 64, j, phase, idum, error, tclock1, tclock2, clock_rate,i
    real(real64), allocatable :: a(:), b(:), x(:), x_other(:), a_tri(:), b_tri(:), c_tri(:)
    integer(int32), allocatable :: ja(:), ia(:)
    real(real64) :: elapsed_time, cpu1, cpu2

    allocate(a((matDimension - 2)*3 + 4), ja((matDimension - 2)*3 + 4), ia(matDimension + 1), b(matDimension), x(matDimension), a_tri(matDimension-1), b_tri(matDimension), c_tri(matDimension-1), x_other(matDimension))
    x = 0
    b = 0
    a_tri = 1
    b_tri = -2
    c_tri = 1
    a(1:2) = (/-2, 1/)
    ja(1:2) = (/1,2/)
    a(size(a)-1:size(a)) = (/1,-2/)
    ja(size(a)-1:size(a)) = (/matDimension-1, matDimension/)
    ia(1:2) = (/1, 3/)
    ia(matDimension + 1) = (matDimension - 2)*3 + 5
    print *, size(a_tri)
    print *, size(b_tri)
    print *, size(c_tri)
    do j = 2, matDimension-1
        a(3 * (j-1): 3*(j-1) + 2) = (/1, -2, 1/)
        ja(3 * (j-1): 3*(j-1) + 2) = (/j-1, j, j+1/)
        ia(j + 1) = j*3
    end do
    do j = 1, matDimension
        b(j) = j-1
    end do
    mtype = 1
    solver = 0
    pt=0 ! important !
    phase = 11
    call pardisoinit(pt,mtype,iparm)

    call pardiso(pt, maxfct, mnum, mtype, phase, matDimension, a, ia, ja, idum, nrhs, iparm, msglvl, b, b, error)
    print *, error

    phase = 22
    call pardiso(pt, maxfct, mnum, mtype, phase, matDimension, a, ia, ja, idum, nrhs, iparm, msglvl, b, b, error)
    print *, error

    phase = 33
    !call cpu_time(cpu1)
    call system_clock(tclock1)
    do i = 1, 100000
        call pardiso(pt, maxfct, mnum, mtype, phase, matDimension, a, ia, ja, idum, nrhs, iparm, msglvl, b, x, error)
    end do
    call system_clock(tclock2, clock_rate)
    !call cpu_time(cpu2)
    elapsed_time = float(tclock2 - tclock1) / float(clock_rate)
    print *, "Elapsed wall time is:", elapsed_time, "seconds"
    print *, "CPU time is:", cpu2-cpu1, "seconds"
    print *, x
    
    ! call system_clock(tclock1)
    ! do i = 1, 100000
    !     call solve_tridiag(a_tri, b_tri, c_tri, b, x_other, matDimension)
    ! end do
    ! call system_clock(tclock2, clock_rate)
    ! elapsed_time = float(tclock2 - tclock1) / float(clock_rate)
    ! print *, "Elapsed time is:", elapsed_time, "seconds"
    ! print *, x_other



contains

    subroutine solve_tridiag(a,b,c,d,x,n)
        implicit none
    !	 a - sub-diagonal (means it is the diagonal below the main diagonal)
    !	 b - the main diagonal
    !	 c - sup-diagonal (means it is the diagonal above the main diagonal)
    !	 d - right part
    !	 x - the answer
    !	 n - number of equations

        !integer,parameter :: r8 = kind(1.d0)

        integer(int32),intent(in) :: n
        real(real64), intent(in) :: a(n-1),b(n),c(n-1),d(n)
        real(real64), intent(in out) :: x(n)
        real(real64) :: cp(n),dp(n)
        real(real64) :: m
        integer(int32) :: i

    ! initialize c-prime and d-prime
        cp(1) = c(1)/b(1)
        dp(1) = d(1)/b(1)
    ! solve for vectors c-prime and d-prime
        do i = 2,n
            m = b(i)-cp(i-1)*a(i-1)
            cp(i) = c(i)/m
            dp(i) = (d(i)-dp(i-1)*a(i-1))/m
        end do
    ! initialize x
        x(n) = dp(n)
    ! solve for x from the vectors c-prime and d-prime
        do i = n-1, 1, -1
            x(i) = dp(i)-cp(i)*x(i+1)
        end do

    end subroutine solve_tridiag

end program pardiso_unsym

