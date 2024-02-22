

program main
    use iso_fortran_env, only: int32, real64
    implicit none
    interface
        subroutine dgttrf(n, dl, d, du, du2, ipiv, info)
            integer, intent(in) :: n
            integer, intent(out) :: ipiv(n), info
            real(kind=8), intent(in out) :: dl(n-1), d(n), du(n-1), du2(n-2)

        end subroutine dgttrf

        subroutine dgttrs(trans, n, nrhs, dl, d, du, du2, ipiv, b, ldb, info)
            character*1, intent(in) :: trans
            integer, intent(in) :: n, nrhs, ldb, ipiv(n)
            integer, intent(in out) :: info
            real(kind=8), intent(in) :: dl(n-1), d(n), du(n-1), du2(n-2)
            real(kind=8), intent(in out) :: b(n, nrhs)

        end subroutine dgttrs

    end interface
    integer(int32) :: matDimension = 10000
    integer(int32) :: j, phase, idum, error, tclock1, tclock2, clock_rate,i, info, mtype, solver, pt(64), iparm(64), maxfct, mnum, msglvl
    real(real64), allocatable :: b(:), b_other(:), a_tri(:), b_tri(:), c_tri(:), a(:), x_other(:)
    integer(int32), allocatable :: ja(:), ia(:)
    real(real64) :: elapsed_time

    allocate(a((matDimension - 2)*3 + 4), ja((matDimension - 2)*3 + 4), ia(matDimension + 1), b(matDimension), b_other(matDimension), a_tri(matDimension-1), b_tri(matDimension), c_tri(matDimension-1), x_other(matDimension))
    
    a_tri = 1
    b_tri = -2
    c_tri = 1
    call system_clock(count_rate = clock_rate)
    do j = 1, matDimension
        b(j) = j-1
        b_other(j) = j-1
    end do
    call mkl_set_num_threads(32)

    ! call dgttrf(matDimension, a_tri, b_tri, c_tri, du2, ipiv, info)

    ! call system_clock(tclock1)
    ! call dgttrs('N', matDimension, 1, a_tri, b_tri, c_tri, du2, ipiv, b_other, matDimension, info)
    ! call system_clock(tclock2)
    ! elapsed_time = real(tclock2 - tclock1) / real(clock_rate)
    ! print *, "Elapsed wall time for dgttrs is:", elapsed_time, "seconds"

    a(1:2) = (/-2, 1/)
    ja(1:2) = (/1,2/)
    a(size(a)-1:size(a)) = (/1,-2/)
    ja(size(a)-1:size(a)) = (/matDimension-1, matDimension/)
    ia(1:2) = (/1, 3/)
    ia(matDimension + 1) = (matDimension - 2)*3 + 5
    do j = 2, matDimension-1
        a(3 * (j-1): 3*(j-1) + 2) = (/1, -2, 1/)
        ja(3 * (j-1): 3*(j-1) + 2) = (/j-1, j, j+1/)
        ia(j + 1) = j*3
    end do
    mtype = 1
    solver = 0
    pt=0 ! important !
    iparm = 0
    phase = 11
    maxfct = 1
    mnum = 1
    msglvl = 0
    call pardisoinit(pt,mtype,iparm)
    call pardiso(pt, maxfct, mnum, mtype, phase, matDimension, a, ia, ja, idum, 1, iparm, msglvl, b_other, x_other, error)
  
   
    phase = 22
    call pardiso(pt, maxfct, mnum, mtype, phase, matDimension, a, ia, ja, idum, 1, iparm, msglvl, b_other, x_other, error)
   
    !call cpu_time(cpu1)
    call system_clock(tclock1)
    call solve_tridiag(matDimension, a_tri, c_tri, b_tri, b)
    call system_clock(tclock2)
    !call cpu_time(cpu2)
    print *, tclock2-tclock1
    elapsed_time = real(tclock2 - tclock1) / real(clock_rate)
    print *, "Elapsed wall time for triDiag is:", elapsed_time, "seconds"
   
    phase = 33
    !call cpu_time(cpu1)
    call system_clock(tclock1)
    call pardiso(pt, maxfct, mnum, mtype, phase, matDimension, a, ia, ja, idum, 1, iparm, msglvl, b_other, x_other, error)
    call system_clock(tclock2)
    !call cpu_time(cpu2)
    elapsed_time = real(tclock2 - tclock1) / real(clock_rate)
    print *, (tclock2 - tclock1)
    print *, "Elapsed wall time for pardiso is:", elapsed_time, "seconds"
    

    
    


contains

    subroutine solve_tridiag(n, diagLower, diagUpper, diag, b)
        ! General tridiagonal solver replace solution to b
        integer(int32), intent(in) :: n
        real(real64), intent(in out) :: b(n)
        real(real64), intent(in) :: diagLower(n-1), diagUpper(n-1), diag(n)
        integer(int32) :: i !n size dependent on how many points are inside (not boundary), so how many in rhs equation
        real(real64) :: m, cp(n-1),dp(n)

    ! initialize c-prime and d-prime
        cp(1) = diagUpper(1)/diag(1)
        dp(1) = b(1)/diag(1)
    ! solve for vectors c-prime and d-prime
        do i = 2,n-1
            m = diag(i)-cp(i-1)*diagLower(i-1)
            cp(i) = diagUpper(i)/m
            dp(i) = (b(i)-dp(i-1)*diagLower(i-1))/m
        end do
        dp(n) = (b(n)-dp(n-1)*diagLower(n-1))/(diag(n)-cp(n-1)*diagLower(n-1))
        b(n) = dp(n)
        do i = n-1, 1, -1
            b(i) = dp(i)-cp(i)*b(i+1)
        end do
    end subroutine solve_tridiag

    function triMul(n, diagLower, diagUpper, diag, x) result(res)
        integer(int32), intent(in) :: n
        real(real64), intent(in) :: diag(n), diagUpper(n-1), diagLower(n-1)
        real(real64), intent(in) :: x(n)
        integer(int32) :: i
        real(real64) :: res(n)
        res(1) = x(1) * diag(1) + x(2) * diagUpper(1)
        do i = 2, n-1
            res(i) = x(i) * diag(i) + x(i-1) * diagLower(i-1) + x(i+1) * diagUpper(i)
        end do
        res(n) = x(n) * diag(n) + x(n-1) * diagLower(n-1)
    end function triMul

end program main