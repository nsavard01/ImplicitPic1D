subroutine nitstb(n, xcur, fcur, fcnrm, step, eta, f, jacv, rpar, &
    ipar, ijacv, irpre, iksmax, ifdord, nfe, njve, nrpre, nli, &
    r, rtil, p, phat, v, t, rwork1, rwork2, rsnrm, dinpr, dnorm, itrmks)

    implicit none

    integer, intent(in) :: n, ijacv, irpre, iksmax, ifdord, nfe, njve, nrpre, nli, itrmks
    double precision, intent(inout) :: xcur(n), fcur(n), fcnrm, step(n), eta, rpar(*)
    double precision, intent(inout) :: r(n), rtil(n), p(n), phat(n), v(n), t(n)
    double precision, intent(out) :: rsnrm
    external :: f, jacv, dinpr, dnorm

    ! Other declarations

    ! ...

    integer :: i, istb, itask, itrmjv
    double precision :: abstol, alpha, beta, omega, rho, rhomns, tau, temp
    double precision, external :: dlamch
    double precision, parameter :: sfmin = 0.0d0

    ! ------------------------------------------------------------------------
    ! Explanation:
    ! ...
    ! ------------------------------------------------------------------------

    ! If finite-differences are used to evaluate J*v products (ijacv= 0), then 
    ! ijacv is set to -1 within this subroutine to signal to nitjv that the 
    ! order of the finite-difference formula is to be determined by ifdord. 
    ! The original value ijacv= 0 is restored on return. 
    ! ------------------------------------------------------------------------
    if (ijacv == 0) ijacv = -1

    ! ------------------------------------------------------------------------
    ! Set the stopping tolerance, initialize the step, etc. 
    ! ------------------------------------------------------------------------
    rsnrm = fcnrm
    abstol = eta * rsnrm
    step = 0.0d0
    istb = 0

    ! ------------------------------------------------------------------------ 
    ! For printing:
    ! ------------------------------------------------------------------------
    if (iplvl >= 3) then
        write(ipunit, *)
        write(ipunit, 800) eta
    endif
800 format('nitstb:  eta =', 1pd10.3)
    if (iplvl >= 4) then
        write(ipunit, 810)
        write(ipunit, *)
        write(ipunit, 820) istb, rsnrm
    endif
810 format('nitstb:  BiCGSTAB iteration no. (parts a and b)',
     $ ' linear residual norm, ')
820 format(5x, i4, 5x, 1pd10.3)

    ! ------------------------------------------------------------------------
    ! Top of the iteration loop. 
    ! ------------------------------------------------------------------------
100 continue
    istb = istb + 1
    nli = nli + 1

    ! ------------------------------------------------------------------------
    ! Perform the first "half-iteration". 
    ! ------------------------------------------------------------------------
    rho = dinpr(n, rtil, 1, r, 1)
    if (istb == 1) then
        p = r
    else
        if (abs(rhomns) < sfmin * abs(rho)) then
            itrmks = 4
            goto 900
        else
            beta = (rho / rhomns) * (alpha / omega)
            p = -omega * v + beta * p + r
        endif
    endif

    if (irpre == 0) then
        phat = p
    else
        itask = 2
        call nitjv(n, xcur, fcur, f, jacv, rpar, ipar, ijacv, &
            ifdord, itask, nfe, njve, nrpre, p, phat, &
            rwork1, rwork2, dnorm, itrmjv)
        if (itrmjv > 0) then
            itrmks = 2
            goto 900
        endif
    endif

    itask = 0
    call nitjv(n, xcur, fcur, f, jacv, rpar, ipar, ijacv, &
        ifdord, itask, nfe, njve, nrpre, phat, v, &
        rwork1, rwork2, dnorm, itrmjv)

    if (itrmjv > 0) then
        itrmks = 1
        goto 900
    endif

    tau = dinpr(n, rtil, 1, v, 1)

    if (abs(tau) < sfmin * abs(rho)) then
        itrmks = 4
        goto 900
    else
        alpha = rho / tau
    endif

    r = r - alpha * v
    step = step + alpha * phat
    rsnrm = dnorm(n, r, 1)

    ! ------------------------------------------------------------------------ 
    ! For printing:
    ! ------------------------------------------------------------------------
    if (iplvl >= 4) then
        write(ipunit, 830) istb, rsnrm
830 format(5x, i4, '.a', 3x, 1pd10.3)
    endif

    ! ------------------------------------------------------------------------
    ! Test for termination. 
    ! ------------------------------------------------------------------------
    if (rsnrm <= abstol) then
        itrmks = 0
        goto 900
    endif

    ! ------------------------------------------------------------------------
    ! Perform the second "half-iteration". 
    ! ------------------------------------------------------------------------
    if (irpre == 0) then
        phat = r
    else
        itask = 2
        call nitjv(n, xcur, fcur, f, jacv, rpar, ipar, ijacv, &
            ifdord, itask, nfe, njve, nrpre, r, phat, &
            rwork1, rwork2, dnorm, itrmjv)
        if (itrmjv > 0) then
            itrmks = 2
            goto 900
        endif
    endif

    itask = 0
    call nitjv(n, xcur, fcur, f, jacv, rpar, ipar, ijacv, &
        ifdord, itask, nfe, njve, nrpre, phat, t, &
        rwork1, rwork2, dnorm, itrmjv)

    if (itrmjv > 0) then
        itrmks = 1
        goto 900
    endif

    tau = dnorm(n, t, 1)
    tau = tau * tau
    temp = dinpr(n, t, 1, r, 1)

    if (tau <= sfmin * abs(temp)) then
        itrmks = 4
        goto 900
    else
        omega = temp / tau
    endif

    if (abs(omega) < sfmin * abs(alpha)) then
        itrmks = 4
        goto 900
    endif

    r = r - omega * t
    step = step + omega * phat
    rsnrm = dnorm(n, r, 1)

    ! ------------------------------------------------------------------------ 
    ! For printing:
    ! ------------------------------------------------------------------------
    if (iplvl >= 4) then
        write(ipunit, 840) istb, rsnrm
840 format(5x, i4, '.b', 3x, 1pd10.3)
    endif

    ! ------------------------------------------------------------------------
    ! Test for termination. 
    ! ------------------------------------------------------------------------
    if (rsnrm <= abstol) then
        itrmks = 0
        goto 900
    endif

    if (istb >= iksmax) then
        itrmks = 3
        goto 900
    endif

    ! ------------------------------------------------------------------------
    ! If continuing, update and return to the top of the iteration loop. 
    ! ------------------------------------------------------------------------
    rhomns = rho
    goto 100

    ! ------------------------------------------------------------------------
    ! Bottom of the iteration loop. 
    ! ------------------------------------------------------------------------
    ! ------------------------------------------------------------------------
    ! All returns made here.
    ! ------------------------------------------------------------------------
900 continue

    ! ------------------------------------------------------------------------ 
    ! For printing:
    ! ------------------------------------------------------------------------
    if (iplvl >= 3) then
        write(ipunit, *)
        if (itrmks /= 1 .and. itrmks /= 2) then
            write(ipunit, 850) itrmks, rsnrm
850         format('nitstb:  itrmks =', i2, '   final lin. res. norm =', 1pd10.3)
        else
            write(ipunit, 860) itrmks
860         format('nitstb: itrmks:', i4)
        endif
    endif

    ! ------------------------------------------------------------------------
    ! If ijacv = -1, then restore it to the original value ijacv = 0. 
    ! ------------------------------------------------------------------------
    if (ijacv == -1) ijacv = 0

    return
end subroutine nitstb
