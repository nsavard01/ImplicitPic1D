subroutine nitdrv(n, xcur, fcur, xpls, fpls, &
    step, f, jacv, rpar, ipar, ftol, stptol, nnimax, &
    ijacv, ikrysl, kdmax, irpre, iksmax, iresup, ifdord, &
    ibtmax, ieta, iterm, nfe, njve, nrpre, nli, nni, nbt, &
    rwork, dinpr, dnorm)

 implicit none

 integer, intent(in) :: n, nnimax, ijacv, ikrysl, kdmax, irpre, iksmax, iresup, ifdord, ibtmax, ieta, iterm
 double precision, intent(inout) :: xcur(n), fcur(n), xpls(n), fpls(n), step(n)
 double precision, intent(inout) :: rpar(*), rwork(*)
 integer, intent(inout) :: ipar(*)
 double precision, intent(in) :: ftol, stptol, dinpr, dnorm
 external f, jacv

 ! Remaining declarations
 double precision :: alpha, epsmach, eta, etamin, fcnrm, flmnrm, fpnrm, gamma, oftjs, oftlm, redfac, rsnrm
 double precision :: stpnrm, temp
 integer :: ibt, itrmbt, itrmf, itrmks, lrr, lsvbig, lsvsml, lvv, lw, lr, ld, lrcgs, lrtil, lp, lphat, lq, lu
 integer :: lv, lt, lrwork, ly, kdmaxp1, ncall

 double precision, external :: dlamch
 external nitbd

 ! Optional common blocks
 integer :: infe, injve, inrpre, inli

 ! Data initialization
 data ncall / 0 /
 save ncall, alpha, gamma, epsmach

 ! Initialize
 if ( ncall .eq. 0 ) epsmach = 2.0d0 * dlamch( 'e' )
 ncall = 1
 if (ieta .eq. 0) alpha = choice1_exp
 if (ieta .eq. 2) alpha = choice2_exp
 if (ieta .eq. 2) gamma = choice2_coef
 if (ieta .eq. 3) then 
   eta = etafixed
 else
   eta = 0.5d0
 endif
 nfe = 0
 njve = 0
 nrpre = 0
 nli = 0
 nni = 0
 nbt = 0

 ! Nonlinear iteration loop
 do
   ! Evaluate f at initial x and initialize eta
   call f(n, xcur, fcur, rpar, ipar, itrmf)
   if ( itrmf /= 0 ) then
     iterm = 2
     exit
   endif
   nfe = nfe + 1
   fcnrm = dnorm(n, fcur, 1)

   ! For printing
   if (iplvl >= 1) then 
     write(ipunit,*)
     write(ipunit,*)
     write(ipunit,*) 'nitdrv:  Beginning nonlinear iterations.'
   endif

   ! Test for stopping
   if (fcnrm <= ftol) then 
     iterm = 0
     exit
   endif
   if (nni > 0 .and. stpnrm <= stptol .and. itrmks == 0) then 
     iterm = 0
     exit
   endif
   if (nni >= nnimax) then 
     iterm = 1
     exit
   endif

   ! Compute the (trial) inexact Newton step with the Krylov solver
   ! Update data in nitinfo to mark the start of a new inexact Newton step
   newstep = 0
   instep = nni
   fcurnrm = fcnrm

   ! If ikrysl = 0, apply GMRES, using fpls as a work array
   if (ikrysl == 0) then 
     kdmaxp1 = kdmax + 1
     lvv = 1
     lrr = lvv + n * kdmaxp1
     lsvbig = lrr + kdmax * kdmax
     lsvsml = lsvbig + kdmax
     lw = lsvsml + kdmax
     call nitgm(n, xcur, fcur, fcnrm, step, eta, f, jacv, nnimax, ijacv, kdmax, rpar, ipar, lrwork(lvv), lrwork(lrr), lrwork(lsvbig), lrwork(lsvsml), lrwork(lw), irpre, iresup, ifdord, ibtmax, ieta, iterm, stptol, ljcr, nfe, nli, nni, nbt)
   else
     ! Otherwise, apply Newton-GMRES, using fpls as a work array
     lv = 1
     lt = lv + n
     lr = lt + n * n
     ld = lr + n * n
     lrcgs = ld + n
     lrtil = lrcgs + n
     lp = lrtil + n * n
     lphat = lp + n * n
     lq = lphat + n * n
     lu = lq + n * n
     lrwork(ld + n) = 0.0d0
     lrwork(ld + n + 1) = 1.0d0
     call nitngm(n, xcur, fcur, fcnrm, step, eta, f, jacv, nnimax, ijacv, kdmax, rpar, ipar, lrwork(lv), lrwork(lt), lrwork(lr), lrwork(ld), lrwork(lrcgs), lrwork(lrtil), lrwork(lp), lrwork(lphat), lrwork(lq), lrwork(lu), irpre, iresup, ifdord, ibtmax, ieta, iterm, stptol, lsvbig, lsvsml, lrwork(lw), nfe, njve, nrpre, nli, nni, nbt, iksmax)
   endif

   ! If the Krylov solver fails or the Newton step is too small
   if (iterm > 1 .or. newstep == 0) then 
     iterm = 1
     exit
   endif
 end do

 ! End nonlinear iteration loop
 return

end subroutine nitdrv
