

program pardiso_unsym
    use iso_fortran_env, only: int32, real64
    implicit none
    integer(int32) :: mtype, pt(64), solver, iparm(64), mnum = 1, maxfct = 1, nrhs = 1, msglvl = 0, matDimension = 8, j, phase, idum, error
    real(real64), allocatable :: a(:), b(:), x(:)
    integer(int32), allocatable :: ja(:), ia(:)

    allocate(a((matDimension - 2)*3 + 4), ja((matDimension - 2)*3 + 4), ia(matDimension + 1), b(matDimension), x(matDimension))
    x = 0
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
    do j = 1, matDimension
        b(j) = j-1
    end do
    print *, a
    print *, ja
    print *, ia
    print *, b
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
    call pardiso(pt, maxfct, mnum, mtype, phase, matDimension, a, ia, ja, idum, nrhs, iparm, msglvl, b, x, error)
    print *, error

    print *, x

end program pardiso_unsym

