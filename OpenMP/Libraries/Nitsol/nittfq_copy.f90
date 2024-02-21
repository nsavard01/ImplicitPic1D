subroutine nittfq (n, xcur, fcur, fcnrm, step, eta, f, jacv, rpar, &
    ipar, ijacv, irpre, iksmax, ifdord, nfe, njve, nrpre, nli, r, &
    rcgs, rtil, d, p, q, u, v, y, rwork1, rwork2, rsnrm, dinpr, &  
        dnorm, itrmks )

implicit none

integer ifdord, ijacv, iksmax, irpre, itrmks, n, nfe, njve, nrpre, nli

integer ipar(*)

double precision eta, fcnrm, rsnrm

double precision d(n), fcur(n), p(n), q(n), rcgs(n), r(n), &
    rpar(*), rtil(n), rwork1(n), rwork2(n), step(n), & 
    u(n), v(n), xcur(n), y(n)

double precision dinpr, dnorm

external f, jacv, dinpr, dnorm

! ------------------------------------------------------------------------
!
! This is nittfq v0.1, the TFQMR routine for determining (trial) inexact 
! Newton steps. The original reference is R. W. Freund, "A Transpose-Free
! Quasi-Minimal Residual Algorithm for Non-Hermitian Linear Systems", 
! SIAM J. Sci. Comput., 14 (1993), pp. 470-482.  The implementation here
! is based on the right preconditioned algorithm that appears in 
! J. N. Shadid and R. S. Tuminaro, "A Comparison of Preconditioned
! Nonsymmetric Krylov Methods on a Large-Scale MIMD Machine", SIAM J.
! Sci. Comput., 15 (1994), pp. 440-459.
!
! ------------------------------------------------------------------------
!
! Explanation:
!
!  n       = dimension of the problem.
!
!  xcur    = vector of length n, current approximate solution. 
!
!  fcur    = vector of length n, value of f at xcur.
!
!  fcnrm   = norm of fcur. 
!
!  step    = vector of length n, (trial) step. 
!
!  eta     = relative residual reduction factor. 
!
!  f      = name of user-supplied subroutine for evaluating the function 
!           the zero of which is sought; this routine has the form
!
!                 subroutine f(n, xcur, fcur, rpar, ipar, itrmf)
!
!           where xcur is the array containing the current x value, fcur 
!           is f(xcur) on output, rpar and ipar are, respectively, real 
!           and integer parameter/work arrays for use by the subroutine,
!           and itrmf is an integer termination flag.  The meaning of
!           itrmf is as follows:
!             0 => normal termination; desired function value calculated.
!             1 => failure to produce f(xcur).
! 
!  jacv   = name of user-supplied subroutine for evaluating J*v or 
!           P(inverse)*v, where J is the Jacobian of f and P is a 
!           right preconditioning operator. If neither analytic J*v 
!           evaluations nor right preconditioning is used, this can 
!           be a dummy subroutine; if right preconditioning is used but 
!           not analytic J*v evaluations, this need only evaluate 
!           P(inverse)*v. The form is 
!
!           subroutine jacv(n, xcur, fcur, ijob, v, z, rpar, ipar, itrmjv)
!
!           where xcur and fcur are vectors of length n containing the 
!           current x and f values, ijob is an integer flag indicating 
!           which product is desired, v is a vector of length n to be 
!           multiplied, z is a vector of length n containing the desired 
!           product on output, rpar and ipar are, respectively, real 
!           and integer parameter/work arrays for use by the subroutine, 
!           and itrmjv is an integer termination 
!           flag. The meaning of ijob is as follows: 
!             0 => z = J*v
!             1 => z = P(inverse)*v 
!           The meaning of itrmjv is as follows:
!             0 => normal termination; desired product evaluated. 
!             1 => failure to produce J*v.
!             2 => failure to produce P(inverse)*v. 
!           This subroutine is called only from nitjv, and is always 
!           called with v .ne. z. 
!
!  rpar    = real parameter/work array passed to the f and jacv routines. 
!
!  ipar    = integer parameter/work array passed to the f and jacv routines. 
!
!  ijacv   = flag for determining method of J*v evaluation.
!              0 => finite-difference evaluation (default). 
!              1 => analytic evaluation. 
!
!  irpre   = flag for right preconditioning. 
!              0 => no right preconditioning
!              1 => right preconditioning
!
!  iksmax  = maximum allowable number of TFQMR iterations. 
!
!  ifdord  = order of the finite-difference formula used in BiCGSTAB 
!            when J*v products are evaluated using finite-differences. 
!            When ijacv = 0 on input to nitsol, ifdord is set to 1, 2, or 
!            4 in nitsol; otherwise, it is irrelevant. When ijacv = 0 on 
!            input to this subroutine, ifdord determines the order of the 
!            finite-difference formula used at each BiCGSTAB iteration 
!            (default 1). In this case, ijacv is set to -1 below to 
!            signal to nitjv that the order of the finite-difference 
!            formula is to be determined by ifdord. The original value 
!            ijacv = 0 is restored on return. 
!            
!  nfe     = number of function evaluations.
!
!  njve    = number of J*v evaluations. 
!
!  nrpre   = number of P(inverse)*v evaluations.
!
!  nli     = number of linear iterations.
!
!  r       = residual vector (for the QMR process)
!
!  rcgs    = residual vector (of the underlying CGS process)
!
!  rtil    = 'shadow' residual vector used in bi-orthogonalization
!
!  d       = vector used in TFQMR
!
!  p       = vector used in TFQMR
!
!  q       = vector used in TFQMR
!
!  u       = vector used in TFQMR
!
!  v       = vector used in TFQMR
!
!  y       = vector used in TFQMR
!
!  rwork1  = work vector, passed on to nitjv
!
!  rwork2  = work vector, passed on to nitjv
!
!  rsnrm   = TFQMR residual norm on return. 
!
!  dinpr   = inner-product routine, either user-supplied or blas ddot. 
!
!  dnorm   = norm routine, either user-supplied or blas dnrm2. 
!
!  itrmks  = termination flag; values have the following meanings: 
!              0 => normal termination: acceptable step found. 
!              1 => J*v failure in nitjv. 
!              2 => P(inverse)*v failure in nitjv. 
!              3 => acceptable step not found in iksmax TFQMR iterations. 
!              4 => TFQMR breakdown.
!              5 => floating point error (the underlying CGS iteration
!                   has probably blown up)
!
!             Note: On return, nitsol terminates if itrmks is 1 or 2. 
!             If itrmks is 3 or 4, nitsol may terminate or continue. 
!             In this event, the step returned is a meaningful inexact 
!             Newton step only if the residual norm has been reduced. 
!             A decision on termination/continuation is made in nitdrv 
!             according to whether there is sufficient residual norm 
!             reduction, even though the desired inexact Newton condition 
!             may not hold.  
!
! -------------------------------------------------------------------------

!  Subroutines required by this and all called routines:

!    user supplied:  nitjv

!    nitsol routines: none

!    BLAS routines- dcopy, daxpy, dscal, dswap

!    LAPACK routines - dlamch

!    user supplied or BLAS:  dinpr, dnorm

!    explanation: In nitsol, dinpr and dnorm are set to either the BLAS 
!    ddot and dnrm2 routines or the user-supplied routines. 

! This subroutine called by: nitdrv

! Subroutines called by this subroutine: daxpy, dcopy, dscal, dswap, dinpr,
!    dlamch, dnorm, nitjv, dlamch

! Common block: 

include 'nitprint.h'

! If diagnostic information is desired, include this common block in the 
! main program and set iplvl and ipunit according to the following: 
!
!     iplvl = 0 => no printout
!           = 1 => iteration numbers and F-norms
!           = 2 => ... + some stats, step norms, and linear model norms
!           = 3 => ... + some Krylov solver and backtrack information
!           = 4 => ... + more Krylov solver and backtrack information
!
!     ipunit = printout unit number.

!  Parameters-

double precision zero,       one
parameter      ( zero=0.0d0, one=1.0d0 )

!  Local variables-

integer i
integer itask
integer itrmjv
integer itfq
integer k

double precision alpha
double precision abstol
double precision beta
double precision c
double precision cgsnorm
double precision qmreta
double precision omega
double precision rho
double precision rho_old
double precision sigma
double precision t
double precision tau
double precision theta

character*2 ab(0:1)
data ab /'.a','.b'/

double precision sfmin
data sfmin /zero/

!  External subroutines-

double precision dlamch

external daxpy
external dcopy
external dlamch
external dscal
external dswap
external nitjv

!  Intrinsics-

intrinsic abs
intrinsic dble
intrinsic sqrt

!  Start of executable code-

!  Initialize sfmin only on first entry.

if ( sfmin .eq. zero ) sfmin = dlamch( 's' )

! If finite-differences are used to evaluate J*v products (ijacv= 0), then 
! ijacv is set to -1 within this subroutine to signal to nitjv that the 
! order of the finite-difference formula is to be determined by ifdord. 
! The original value ijacv= 0 is restored on return. 

if (ijacv .eq. 0) ijacv = -1 

! Set the stopping tolerance, initialize the step, etc. 

rsnrm = fcnrm
abstol = eta*rsnrm
do 10 i = 1, n
step(i) = zero
10   continue
itfq = 0

! For printing:

if ( iplvl .ge. 3 ) then 
write(ipunit,*) 
write(ipunit,800) eta 
endif
if ( iplvl .ge. 4 ) then 
write(ipunit,810) 
write(ipunit,*) 
write(ipunit,820) itfq, rsnrm 
endif

!  Initialize residual and work vectors.

call dcopy( n, fcur, 1, rcgs, 1 )
call dscal( n, -one, rcgs, 1 )

!  Choice here is rtil = r.

call dcopy( n, rcgs, 1, rtil, 1 )

if ( irpre .eq. 0 ) then
call dcopy( n, rcgs, 1, p, 1 )
else
itask = 2
call nitjv( n, xcur, fcur, f, jacv, rpar, ipar, ijacv, ifdord,
&               itask, nfe, njve, nrpre, rcgs, p, rwork1, rwork2,
&                                                   dnorm, itrmjv )    
if ( itrmjv .ne. 0 ) then
itrmks = 2
goto 900
endif
endif
call dcopy( n, p, 1, u, 1 )
rho = dinpr( n, rcgs, 1, rtil, 1 )
do 20 i = 1, n
d(i) = zero
20   continue

itask = 0
call nitjv( n, xcur, fcur, f, jacv, rpar, ipar, ijacv, ifdord,
&            itask, nfe, njve, nrpre, p, v, rwork1, rwork2,
&                                                   dnorm, itrmjv )
if ( itrmjv .ne. 0 ) then
itrmks = 1
goto 900
end if

alpha = zero
omega = fcnrm
tau = omega
theta = zero
qmreta = zero

!  Start iterations.

100  continue
itfq = itfq + 1
nli = nli + 1
sigma = dinpr( n, rtil, 1, v, 1 )

!  If sigma = 0 we have a serious breakdown.  We check this condition
!  by trying to detect whether division by sigma causes an overflow.

if ( abs(sigma) .lt. sfmin*abs(rho) ) then
itrmks = 4
goto 900
else
alpha = rho/sigma
endif

!  Need Pv for calculation of q.  First store result in q, then
!  swap some vectors to cast calculation of q as a SAXPY.

if ( irpre .eq. 0 ) then
call dcopy( n, v, 1, q, 1 )
else
itask = 2
call nitjv( n, xcur, fcur, f, jacv, rpar, ipar, ijacv, ifdord,
&               itask, nfe, njve, nrpre, v, q, rwork1, rwork2,
&                                                   dnorm, itrmjv )    
if ( itrmjv .ne. 0 ) then
itrmks = 2
goto 900
endif
endif
call dcopy( n, u, 1, y, 1 )
call dswap( n, q, 1, u, 1 )
call daxpy( n, -alpha, u, 1, q, 1 )
call dcopy( n, y, 1, u, 1 )

!  Update residual.

do 30 i = 1, n
y(i) = u(i) + q(i)
30   continue
itask = 0
call nitjv( n, xcur, fcur, f, jacv, rpar, ipar, ijacv, ifdord,
&            itask, nfe, njve, nrpre, y, v, rwork1, rwork2,
&                                                   dnorm, itrmjv )
if ( itrmjv .ne. 0 ) then
itrmks = 1
goto 900
end if
call daxpy( n, -alpha, v, 1, rcgs, 1 )
cgsnorm = dnorm( n, rcgs, 1 )

!  Check for cgsnorm = NaN.

if ( cgsnorm .ne. cgsnorm ) then
itrmks = 5
goto 900
endif

!  QMR section.

do 60 k = 0, 1

!  Use weighting strategy from (5.11) of Freund reference.

t = qmreta*theta**2
t = t/alpha

if ( k .eq. 0 ) then
do 40 i = 1, n
d(i) = u(i) + t*d(i)
40         continue
omega = sqrt(omega*cgsnorm)
else if ( k .eq. 1 ) then
do 50 i = 1, n
d(i) = q(i) + t*d(i)
50         continue
omega = cgsnorm
endif

theta = omega/tau
c = one/sqrt(one + theta**2)
tau = tau*theta*c
qmreta = alpha*c**2

! For printing:

if ( iplvl .ge. 4 .and. tau .gt. abstol ) then 
write(ipunit,830) itfq, ab(k), tau, '     (estimated)'
endif

call daxpy( n, qmreta, d, 1, step, 1 )

!  Convergence check.  Do a cheap test to save on Jacobi-vector products.
!  In case residual history is requested by iplvl, we must calculate 
!  the QMR residual from scratch.  Note termination is always determined
!  by the smoothed residual, calculated from scratch.

if ( tau .le. abstol ) then

itask = 0
call nitjv( n, xcur, fcur, f, jacv, rpar, ipar, ijacv,
&            ifdord, itask, nfe, njve, nrpre, step, r, rwork1,
&                                           rwork2, dnorm, itrmjv )

!  This calculation of the QMR residual is off by a factor
!  of -1, but we don't care until we return from this routine.

call daxpy( n, one, fcur, 1, r, 1 )
rsnrm = dnorm( n, r, 1 )

! For printing:

if ( iplvl .ge. 4 ) then 
write(ipunit,830) itfq, ab(k), rsnrm, '  (from scratch)'
endif

!  Check for rsnrm = NaN.

if ( rsnrm .ne. rsnrm ) then
itrmks = 5
goto 900
endif

!  If rsnrm is small enough, exit.

if ( rsnrm .lt. abstol ) then
itrmks = 0
goto 900
endif
endif

60   continue

rho_old = rho
rho = dinpr( n, rtil, 1, rcgs, 1 )

!  If rho_old = 0 we have a serious breakdown.  We check this condition
!  by trying to detect whether division by rho_old causes an overflow.

if ( abs(rho_old) .lt. sfmin*abs(rho) ) then
itrmks = 4
goto 900
else
beta = rho/rho_old
endif

if ( irpre .eq. 0 ) then
call dcopy( n, rcgs, 1, v, 1 )
else
itask = 2
call nitjv( n, xcur, fcur, f, jacv, rpar, ipar, ijacv, ifdord,
&               itask, nfe, njve, nrpre, rcgs, v, rwork1, rwork2,
&                                                   dnorm, itrmjv )    
if ( itrmjv .ne. 0 ) then
itrmks = 2
goto 900
endif
endif
call daxpy( n, beta, q, 1, v, 1 )
call dcopy( n, v, 1, u, 1 )
call daxpy( n, beta, p, 1, q, 1 )
call daxpy( n, beta, q, 1, v, 1 )
call dcopy( n, v, 1, p, 1 )

itask = 0
call nitjv( n, xcur, fcur, f, jacv, rpar, ipar, ijacv, ifdord,
&            itask, nfe, njve, nrpre, p, v, rwork1, rwork2,
&                                                   dnorm, itrmjv )
if ( itrmjv .ne. 0 ) then
itrmks = 1
goto 900
end if
if ( itfq .ge. iksmax ) then
itrmks = 3
goto 900
end if

!  Do again

goto 100

!  All returns made here.

900 continue

!  If residual hasn't been updated, force
!  computation of residual from scratch.

if ( rsnrm .eq. fcnrm ) then

itask = 0
call nitjv( n, xcur, fcur, f, jacv, rpar, ipar, ijacv,
&         ifdord, itask, nfe, njve, nrpre, step, r, rwork1,
&                                        rwork2, dnorm, itrmjv )

call daxpy( n, one, fcur, 1, r, 1 )
rsnrm = dnorm( n, r, 1 )
if ( rsnrm .le. abstol ) itrmks = 0

end if

!  Correct residual before returning.

call dscal( n, -one, r, 1 )

! If ijacv = -1, then restore it to the original value ijacv = 0. 

if (ijacv .eq. -1) ijacv = 0 

! For printing:

if ( iplvl .ge. 3 ) then 
write(ipunit,*) 
if (itrmks .ne. 1 .and. itrmks .ne. 2) then 
write(ipunit,840) itrmks, rsnrm 
else
write(ipunit,850) itrmks
endif
endif

return

800  format('nittfq:  eta =', 1pd10.3)
810  format('nittfq:  TFQMR iteration no. (parts a and b),',
&                                        ' linear residual norm')
820  format(5x,i4,5x,1pd10.3)
830  format(5x,i4,a2,3x,1pd10.3,a16)
840  format('nittfq:  itrmks =', i2, '   final lin. res. norm =', 
&                                                          1pd10.3)
850  format('nittfq: itrmks:', i4) 

!  End of nittfq.

end subroutine nittfq
