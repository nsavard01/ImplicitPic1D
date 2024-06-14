subroutine nitjv(n, xcur, fcur, f, jacv, rpar, ipar, ijacv, &
    ifdord, itask, nfe, njve, nrpre, v, z, &
    rwork1, rwork2, dnorm, itrmjv)

    implicit none

    integer, intent(in) :: n, ijacv, ifdord, itask
    integer, intent(inout) :: nfe, njve, nrpre, itrmjv
    integer, intent(inout) :: ipar(*)
    double precision, intent(inout) :: xcur(n), fcur(n), rpar(*)
    double precision, intent(inout) :: v(n), z(n), rwork1(n), rwork2(n)
    external f, jacv, dnorm

    ! Other declarations

    ! ...

    integer :: ijob

    ! ------------------------------------------------------------------------
    ! Explanation:
    ! ...
    ! ------------------------------------------------------------------------

    ! Common block (if necessary)
    ! include 'nitprint.h'

    ! ------------------------------------------------------------------------
    ! If z = J*v is desired (itask = 0), then copy v into rwork1; if 
    ! z = J*P(inverse)*v or z = P(inverse)*v is desired (itask = 1,2), 
    ! then compute P(inverse)*v in rwork1. 
    ! ------------------------------------------------------------------------
    if (itask == 0) then
        call dcopy(n, v, 1, rwork1, 1)
    else
        ijob = 1
        call jacv(n, xcur, fcur, ijob, v, rwork1, rpar, ipar, itrmjv)
        nrpre = nrpre + 1
        if (itrmjv /= 0) go to 900
    end if

    ! ------------------------------------------------------------------------
    ! If only z = P(inverse)*v is desired (itask = 2), then copy rwork1 into 
    ! z and exit.
    ! ------------------------------------------------------------------------
    if (itask == 2) then
        call dcopy(n, rwork1, 1, z, 1)
        go to 900
    end if

    ! ------------------------------------------------------------------------
    ! If z = J*v or z = J*P(inverse)*v is desired (itask = 0, 1), then 
    ! compute J*rwork1 in z by either analytic evaluation (ijacv = 1) or 
    ! finite-differences (ijacv = 0, -1). 
    ! ------------------------------------------------------------------------
    if (ijacv == 1) then
        ijob = 0
        call jacv(n, xcur, fcur, ijob, rwork1, z, rpar, ipar, itrmjv)
    else
        call nitfd(n, xcur, fcur, f, rpar, ipar, ijacv, ifdord, &
            nfe, rwork1, z, rwork2, dnorm, itrmjv)
    end if
    njve = njve + 1

    ! ------------------------------------------------------------------------
    ! All returns made here.
    ! ------------------------------------------------------------------------
 900 continue
    return
end subroutine nitjv

subroutine nitfd(n, xcur, fcur, f, rpar, ipar, ijacv, ifdord, &
    nfe, v, z, rwork, dnorm, itrmjv)

    implicit none

    integer, intent(in) :: n, ijacv, ifdord
    integer, intent(inout) :: nfe, itrmjv
    integer, intent(inout) :: ipar(*)
    double precision, intent(inout) :: xcur(n), fcur(n), rpar(*)
    double precision, intent(inout) :: v(n), z(n), rwork(n)
    external f, dnorm

    ! Other declarations

    ! ...

    integer :: i
    double precision :: eps, epsmach, temp

    ! ------------------------------------------------------------------------
    ! Explanation:
    ! ...
    ! ------------------------------------------------------------------------

    ! Common block (if necessary)
    ! include 'nitprint.h'

    ! ------------------------------------------------------------------------
    ! Set epsmach (machine epsilon) on first call. 
    ! ------------------------------------------------------------------------
    if (ncall == 0) epsmach = 2.0d0 * dlamch('e')
    ncall = 1

    ! ------------------------------------------------------------------------
    ! Compute z = J*v by finite-differences: First, set eps = ||v||for later 
    ! use in computing the difference step; then evaluate the difference 
    ! formula according to ijacv and ifdord. 
    ! ------------------------------------------------------------------------
    eps = dnorm(n, v, 1)
    if (eps == 0.d0) then
        itrmjv = 1
        go to 900
    end if

    ! ------------------------------------------------------------------------
    ! Here ijacv = 0 or ifdord = 1 => first-order forward difference. 
    ! ------------------------------------------------------------------------
    if (ijacv == 0 .or. ifdord == 1) then
        eps = sqrt((1.d0 + dnorm(n, xcur, 1)) * epsmach) / eps
        do i = 1, n
            v(i) = xcur(i) + eps * v(i)
        end do
        call f(n, v, z, rpar, ipar, itrmjv)
        if (itrmjv /= 0) then
            itrmjv = 1
            go to 900
        end if
        nfe = nfe + 1
        do i = 1, n
            z(i) = (z(i) - fcur(i)) / eps
        end do
        itrmjv = 0
        go to 900
    end if

    ! ------------------------------------------------------------------------
    ! Here ijacv = -1 and ifdord = 2 => second-order central difference. 
    ! ------------------------------------------------------------------------
    if (ifdord == 2) then
        eps = (((1.d0 + dnorm(n, xcur, 1)) * epsmach) ** (1.d0 / 3.d0)) / eps
        do i = 1, n
            rwork(i) = xcur(i) + eps * v(i)
        end do
        call f(n, rwork, z, rpar, ipar, itrmjv)
        nfe = nfe + 1
        if (itrmjv /= 0) then
            itrmjv = 1
            go to 900
        end if
        do i = 1, n
            rwork(i) = xcur(i) - eps * v(i)
        end do
        call f(n, rwork, v, rpar, ipar, itrmjv)
        nfe = nfe + 1
        if (itrmjv /= 0) then
            itrmjv = 1
            go to 900
        end if
        temp = 2.d0 * eps
        do i = 1, n
            z(i) = (z(i) - v(i)) / temp
        end do
        itrmjv = 0
        go to 900
    end if

    ! ------------------------------------------------------------------------
    ! Here ijacv = -1 and ifdord = 4 => fourth-order difference. 
    ! ------------------------------------------------------------------------
    if (ifdord == 4) then
        eps = (((1.d0 + dnorm(n, xcur, 1)) * epsmach) ** (1.d0 / 5.d0)) / eps
        do i = 1, n
            rwork(i) = xcur(i) + eps * v(i)
        end do
        call f(n, rwork, z, rpar, ipar, itrmjv)
        nfe = nfe + 1
        if (itrmjv /= 0) then
            itrmjv = 1
            go to 900
        end if
        temp = -eps
        call daxpy(n, temp, v, 1, xcur, 1)
        call f(n, xcur, rwork, rpar, ipar, itrmjv)
        nfe = nfe + 1
        if (itrmjv /= 0) then
            itrmjv = 1
            go to 900
        end if
        do i = 1, n
            z(i) = rwork(i) - z(i)
        end do
        temp = eps / 2.d0
        call daxpy(n, temp, v, 1, xcur, 1)
        call f(n, xcur, rwork, rpar, ipar, itrmjv)
        nfe = nfe + 1
        if (itrmjv /= 0) then
            itrmjv = 1
            go to 900
        end if
        temp = -8.d0
        call daxpy(n, temp, rwork, 1, z, 1)
        call daxpy(n, eps, v, 1, xcur, 1)
        call f(n, xcur, rwork, rpar, ipar, itrmjv)
        nfe = nfe + 1
        if (itrmjv /= 0) then
            itrmjv = 1
            go to 900
        end if
        temp = 8.d0
        call daxpy(n, temp, rwork, 1, z, 1)
        temp = 1.d0 / (6.d0 * eps)
        call dscal(n, temp, z, 1)
        temp = -eps / 2.d0
        call daxpy(n, temp, v, 1, xcur, 1)
        itrmjv = 0
    end if

    ! ------------------------------------------------------------------------
    ! All returns made here.
    ! ------------------------------------------------------------------------
 900 continue
    return
end subroutine nitfd