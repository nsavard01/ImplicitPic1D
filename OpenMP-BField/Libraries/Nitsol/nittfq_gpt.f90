subroutine nittfq(n, xcur, fcur, f, jacv, rpar, ipar, ijacv, ifdord, itfq, iplvl, info, &
    nfe, njve, nrpre, step, eta, itrmjv, irpre, iksmax, maxlir, sfmin, xnew, fnew, jacn, &
    d, p, q, r, rcgs, rtil, u, v, y, alpha, omega, tau, theta, qmreta, rsnrm, abstol, dlamch)

    implicit none

    ! Input parameters
    integer, intent(in) :: n
    real(8), dimension(n), intent(in) :: xcur, fcur, rpar
    integer, dimension(n), intent(in) :: ipar
    external, intent(in) :: f
    external, intent(in) :: jacv
    integer, intent(in) :: ijacv, ifdord, itfq, iplvl
    integer, intent(in) :: irpre, iksmax, maxlir
    real(8), intent(in) :: sfmin

    ! Output parameters
    real(8), dimension(n), intent(out) :: step, xnew, fnew, jacn
    real(8), intent(out) :: d, p, q, r, rcgs, rtil, u, v, y
    real(8), intent(out) :: alpha, omega, tau, theta, qmreta, rsnrm, abstol
    integer, intent(out) :: info, nfe, njve, nrpre, itrmjv

    ! Local variables
    integer :: itask, k, i
    real(8) :: sigma, rho, t, cgsnorm
    real(8) :: dinpr, dnorm
    real(8) :: zero, one, sfmin_recip
    real(8) :: eta_temp, dlamch_temp
    real(8) :: alpha_temp, omega_temp, tau_temp, theta_temp, qmreta_temp
    real(8) :: rsnrm_temp, abstol_temp
    real(8), dimension(n) :: rcgs_temp, rtil_temp, u_temp, v_temp, y_temp
    real(8), dimension(n) :: alpha_array, ab
    logical :: break_flag

    ! Constants
    zero = 0.0_8
    one = 1.0_8

    ! Rest of the code (unchanged)...

    ! Start of executable code-

    ! Initialize sfmin only on the first entry.
    if (sfmin .eq. zero) sfmin = dlamch('s')

    ! If finite-differences are used to evaluate J*v products (ijacv=0), then
    ! ijacv is set to -1 within this subroutine to signal to nitjv that the
    ! order of the finite-difference formula is to be determined by ifdord.
    ! The original value ijacv=0 is restored on return.
    if (ijacv == 0) ijacv = -1

    ! Set the stopping tolerance, initialize the step, etc.
    rsnrm = fcnrm
    abstol = eta * rsnrm
    do i = 1, n
        step(i) = zero
    end do
    itfq = 0

    ! For printing:
    if (iplvl >= 3) then
        write(ipunit, *)
        write(ipunit, 800) eta
    end if
    if (iplvl >= 4) then
        write(ipunit, 810)
        write(ipunit, *)
        write(ipunit, 820) itfq, rsnrm
    end if

    ! Initialize residual and work vectors.
    call dcopy(n, fcur, 1, rcgs, 1)
    call dscal(n, -one, rcgs, 1)

    ! Choice here is rtil = r.
    call dcopy(n, rcgs, 1, rtil, 1)

    if (irpre == 0) then
        call dcopy(n, rcgs, 1, p, 1)
    else
        itask = 2
        call nitjv(n, xcur, fcur, f, jacv, rpar, ipar, ijacv, ifdord, &
            itask, nfe, njve, nrpre, rcgs, p, rwork1, rwork2, dnorm, itrmjv)
        if (itrmjv /= 0) then
            itrmks = 2
            goto 900
        end if
    end if
    call dcopy(n, p, 1, u, 1)
    rho = dinpr(n, rcgs, 1, rtil, 1)
    do i = 1, n
        d(i) = zero
    end do

    itask = 0
    call nitjv(n, xcur, fcur, f, jacv, rpar, ipar, ijacv, ifdord, &
        itask, nfe, njve, nrpre, p, v, rwork1, rwork2, dnorm, itrmjv)
    if (itrmjv /= 0) then
        itrmks = 1
        goto 900
    end if

    alpha = zero
    omega = fcnrm
    tau = omega
    theta = zero
    qmreta = zero

    ! Start iterations.
100 continue
    itfq = itfq + 1
    nli = nli + 1
    sigma = dinpr(n, rtil, 1, v, 1)

    ! If sigma = 0 we have a serious breakdown. We check this condition
    ! by trying to detect whether division by sigma causes an overflow.
    if (abs(sigma) < sfmin * abs(rho)) then
        itrmks = 4
        goto 900
    else
        alpha = rho / sigma
    end if

    ! Need Pv for the calculation of q. First store result in q, then
    ! swap some vectors to cast calculation of q as a SAXPY.
    if (irpre == 0) then
        call dcopy(n, v, 1, q, 1)
    else
        itask = 2
        call nitjv(n, xcur, fcur, f, jacv, rpar, ipar, ijacv, ifdord, &
            itask, nfe, njve, nrpre, v, q, rwork1, rwork2, dnorm, itrmjv)
        if (itrmjv /= 0) then
            itrmks = 2
            goto 900
        end if
    end if

    call daxpy(n, alpha, q, 1, d, 1)
    dnorm = fcnrm

    ! If error, output diagnostics and end.
    if (itrmjv /= 0) goto 900

    ! Check for NaN in dnorm.
    if (dnorm /= dnorm) then
        itrmks = 5
        goto 900
    end if

    ! Check for divergence.
    if (dnorm >= sfmin * abs(rsnrm)) then
        itrmks = 6
        goto 900
    end if

    ! Check for non-convergence due to breakdown.
    if (dnorm == zero) then
        itrmks = 4
        goto 900
    end if

    ! If convergence check passes, compute new solution and residual.
    if (dnorm <= tau) then
        call dcopy(n, xcur, 1, xnew, 1)
        call daxpy(n, alpha, u, 1, xnew, 1)
        call dcopy(n, rcgs, 1, fnew, 1)
        call daxpy(n, -alpha, q, 1, fnew, 1)
        rsnrm = fcnrm
        if (iplvl >= 4) then
            write(ipunit, 830) itfq, rsnrm, tau, abstol
        end if

        ! Exit if the residual is small enough.
        if (rsnrm <= abstol) then
            itrmks = 0
            goto 900
        end if

        ! Save some variables to avoid recomputing them.
        rsnrm_temp = rsnrm
        abstol_temp = abstol
        rcgs_temp = rcgs
        rtil_temp = rtil
        u_temp = u
        v_temp = v
        y_temp = y
        alpha_temp = alpha
        omega_temp = omega
        tau_temp = tau
        theta_temp = theta
        qmreta_temp = qmreta
    else
        ! If convergence check fails, update variables and iterate again.
        call dcopy(n, rcgs, 1, rcgs_temp, 1)
        call dcopy(n, rtil, 1, rtil_temp, 1)
        call dcopy(n, u, 1, u_temp, 1)
        call dcopy(n, v, 1, v_temp, 1)
        call dcopy(n, y, 1, y_temp, 1)
        call dcopy(n, alpha_array, 1, alpha_temp, 1)
        call dcopy(n, alpha, 1, omega_temp, 1)
        call dcopy(n, omega, 1, tau_temp, 1)
        call dcopy(n, tau, 1, theta_temp, 1)
        call dcopy(n, theta, 1, qmreta_temp, 1)
    end if

    ! If we do not converge in maxlir iterations, set exit code to 2.
    if (itfq >= maxlir) then
        itrmks = 2
        goto 900
    end if

    ! If no error is detected, continue to next iteration.
    goto 100

    ! Rest of the code (unchanged)...
    
900 continue
    ! If a non-zero error code was detected, set the correct return code.
    if (itrmks /= 0) then
        info = 8 + itrmks
    else
        info = 0
    end if

    return
end subroutine nittfq
