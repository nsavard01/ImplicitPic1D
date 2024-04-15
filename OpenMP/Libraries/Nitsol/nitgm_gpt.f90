module nitgm_mod
    implicit none
  
    private
  
    public :: nitgm
  
    contains
  
    subroutine nitgm(n, xcur, fcur, fcnrm, step, eta, f, jacv, rpar, &
                     ipar, ijacv, irpre, iksmax, iresup, ifdord, nfe, njve, &
                     nrpre, nli, kdmax, kdmaxp1, vv, rr, svbig, svsml, w, rwork, &
                     rsnrm, dinpr, dnorm, itrmks)
      integer, intent(in) :: n, ipar(*)
      integer, intent(in) :: ijacv, irpre, iksmax, iresup, ifdord, nfe, njve
      integer, intent(inout) :: nrpre, nli, kdmax, kdmaxp1, itrmks
      real(8), intent(inout) :: xcur(n), fcur(n), fcnrm, step(n), eta
      real(8), intent(out) :: rsnrm
      real(8), intent(in) :: rpar(*)
      real(8), intent(inout) :: vv(n, kdmaxp1), rr(kdmax, kdmax), svbig(kdmax), svsml(kdmax), w(kdmax)
      real(8), intent(inout) :: rwork(n)
      external :: f, jacv, dinpr, dnorm
  
      real(8) :: abstol, big, cs, cndmax, epsmach, rsnrm0, sestpr, small, sn, temp
      integer :: i, igm, ijob, itask, itrmjv, kd, kdp1, ncall
  
      real(8) :: dlamch
      external dlamch
  
      ! Initialize local variables
  
      ! Initialize
      if (ncall .eq. 0) then
        epsmach = 2.0d0 * dlamch('e')
        ncall = 1
        cndmax = 1.d0 / (100.d0 * epsmach)
      end if
  
      do i = 1, n
        step(i) = 0.d0
      end do
  
      igm = 0
  
      ! Set the stopping tolerance, etc.
      rsnrm0 = fcnrm
      abstol = eta * rsnrm0
  
      ! Place the normalized initial residual in the first column of the vv array.
      call dcopy(n, fcur, 1, vv(1, 1), 1)
      temp = -1.d0 / fcnrm
      call dscal(n, temp, vv(1, 1), 1)
  
      ! Top of the outer GMRES loop.
    outer_loop: do
  
        kd = 0
        rsnrm = 1.d0
  
        ! Top of the inner GMRES loop.
      inner_loop: do
  
          kd = kd + 1
          kdp1 = kd + 1
          nli = nli + 1
  
          ! Evaluate J*(kd-th Krylov subspace basis vector) in vv(.,kdp1).
          if (irpre .eq. 0) then
            itask = 0
          else
            itask = 1
          end if
          call nitjv(n, xcur, fcur, f, jacv, rpar, ipar, ijacv, ifdord, itask, nfe, njve, nrpre, &
                     vv(1, kd), vv(1, kdp1), rwork, rwork, dnorm, itrmjv)
          if (itrmjv .gt. 0) then
            if (itrmjv .eq. 1) itrmks = 1
            if (itrmjv .eq. 2) itrmks = 2
            exit outer_loop
          end if
  
          ! Do modified Gram-Schmidt.
          do i = 2, kd
            rr(i - 1, kd) = dinpr(n, vv(1, i), 1, vv(1, kdp1), 1)
            call daxpy(n, -rr(i - 1, kd), vv(1, i), 1, vv(1, kdp1), 1)
          end do
          rr(kd, kd) = dnorm(n, vv(1, kdp1), 1)
  
          ! Update the estimates of the largest and smallest singular values.
          if (kd == 1) then
            big = rr(1, 1)
            small = big
            svbig(1) = 1.d0
            svsml(1) = 1.d0
          else
            ijob = 1
            call dlaic1(ijob, kd - 1, svbig, big, rr(1, kd), rr(kd, kd), sestpr, sn, cs)
            big = sestpr
            call dscal(kd - 1, sn, svbig, 1)
            svbig(kd) = cs
            ijob = 2
            call dlaic1(ijob, kd - 1, svsml, small, rr(1, kd), rr(kd, kd), sestpr, sn, cs)
            small = sestpr
            call dscal(kd - 1, sn, svsml, 1)
            svsml(kd) = cs
          end if
  
          ! Terminate if the estimated condition number is too great.
          if (big .ge. small * cndmax) then
            if (kd == 1) then
              itrmks = 5
              exit outer_loop
            else
              kdp1 = kd
              kd = kd - 1
              call daxpy(n, w(kd), vv(1, kdp1), 1, vv(1, 1), 1)
              exit inner_loop
            end if
          end if
  
          ! Normalize vv(.,kdp1).
          temp = 1.d0 / rr(kd, kd)
          call dscal(n, temp, vv(1, kdp1), 1)
  
          ! Update w and the residual norm.
          w(kd) = dinpr(n, vv(1, 1), 1, vv(1, kdp1), 1)
          temp = max(min(w(kd) / rsnrm, 1.0D0), -1.0d0)
          rsnrm = rsnrm * dsin(dacos(temp))
  
          ! For printing: Print iteration information.
          if (iplvl .ge. 4) then
            write(ipunit, *) igm, rsnrm * rsnrm0, big / small
          end if
  
          ! Test for termination of the inner loop.
          if (rsnrm .le. abstol) then
            itrmks = 0
            exit outer_loop
          end if
  
          ! If not terminating the inner loop, update the residual vector and go to the top of the inner loop.
          if (iresup .eq. 0) then
            call daxpy(n, -w(kd), vv(1, kdp1), 1, fcur, 1)
          else
            call nitjv(n, xcur, fcur, f, jacv, rpar, ipar, ijacv, ifdord, 1, nfe, njve, nrpre, &
                       fcur, vv(1, kdp1), rwork, rwork, dnorm, itrmjv)
            if (itrmjv .gt. 0) then
              if (itrmjv .eq. 1) itrmks = 1
              if (itrmjv .eq. 2) itrmks = 2
              exit outer_loop
            end if
          end if
  
        end do ! Inner loop
  
        ! Bottom of inner loop.
  
        ! For printing: Print information after the inner loop.
        if (iplvl .ge. 4) then
          write(ipunit, *) ' GMRES: Outer = ', igm, ', Inner = ', nli, ', Residual Norm = ', rsnrm * rsnrm0
        end if
  
        ! Compute the solution.
        call dgels('N', n, kd, 1, vv, n, w, n, rwork, itask)
  
        ! Use svbig for storage of the original components of w.
        call dcopy(kd, w, 1, svbig(1), 1)
  
        ! Overwrite w with the solution of the upper triangular system.
        call dtrsm('L', 'U', 'N', 'N', kd, 1, 1.d0, rr, kd, w, kd)
  
        ! Now form the linear combination to accumulate the correction in the work vector.
        do i = 1, kd
          call daxpy(n, w(i), vv(1, i + 1), 1, fcur, 1)
        end do
  
        ! If iresup .eq. 0, then update the residual vector by linear combination.
        if (iresup .eq. 0) then
          call dgemv('N', n, kd, 1.d0, vv, n, w, 1, 1.d0, fcur, 1)
        end if
  
        ! If right preconditioning is used, overwrite correction <-- P(inverse)*correction.
        if (irpre .eq. 1) then
          itask = 2
          call nitjv(n, xcur, fcur, f, jacv, rpar, ipar, ijacv, ifdord, itask, nfe, njve, nrpre, &
                     fcur, w, rwork, rwork, dnorm, itrmjv)
          if (itrmjv .gt. 0) then
            if (itrmjv .eq. 1) itrmks = 1
            if (itrmjv .eq. 2) itrmks = 2
            exit outer_loop
          end if
        end if
  
        ! Update the step.
        call daxpy(n, 1.d0, w, 1, step, 1)
  
        ! If iresup .eq. 1, then update the residual vector by direct evaluation.
        if (iresup .eq. 1) then
          call nitjv(n, xcur, fcur, f, jacv, rpar, ipar, ijacv, ifdord, 0, nfe, njve, nrpre, &
                     fcur, fcur, rwork, rwork, dnorm, itrmjv)
          if (itrmjv .gt. 0) then
            if (itrmjv .eq. 1) itrmks = 1
            if (itrmjv .eq. 2) itrmks = 2
            exit outer_loop
          end if
        end if
  
        ! Test for termination.
        if (rsnrm .le. abstol) then
          itrmks = 0
          exit outer_loop
        end if
  
        ! If not terminating, then normalize the initial residual, etc., and return to the top of the outer loop.
        call dcopy(n, fcur, 1, vv(1, 1), 1)
        temp = -1.d0 / fcnrm
        call dscal(n, temp, vv(1, 1), 1)
        igm = igm + 1
  
      end do ! Outer loop
  
      ! All returns made here.
  
    end subroutine nitgm
  
    subroutine nitjv(n, x, f, fjac, fjacv, rpar, ipar, ijacv, ifdord, itask, nfe, njve, nrpre, &
                     vec, prod, work, twork, dnorm, itrmjv)
      ! Implementation of nitjv is required.
    end subroutine nitjv
  
  end module nitgm_mod
  