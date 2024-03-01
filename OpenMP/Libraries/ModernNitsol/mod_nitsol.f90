module mod_nitsol
    use iso_fortran_env, only: int32, real64
    implicit none

    interface

        function ddot(n, x, sx, y, sy)
            use iso_fortran_env, only: int32, real64
            implicit none
            integer(int32), intent(in)                :: n
            real(real64), dimension(*), intent(in) :: x
            integer(int32), intent(in)                :: sx
            real(real64), dimension(*), intent(in) :: y
            integer(int32), intent(in)                :: sy
            real(real64)                           :: ddot
          end function ddot
      
        function dnrm2(n, x, sx)
            use iso_fortran_env, only: int32, real64
            implicit none
            integer(int32), intent(in)                :: n
            real(real64), dimension(*), intent(in) :: x
            integer(int32), intent(in)                :: sx
            real(real64)                           :: dnrm2
        end function dnrm2

        function dlamch(cmach)
            character*1, intent(in) :: cmach
            double precision :: dlamch
        end function

        subroutine daxpy(n, a, x, incx, y, incy)
            integer, intent(in) :: n, incx, incy
            double precision, intent(in) :: a
            double precision, intent(in) :: x(*)
            double precision, intent(in out) :: y(*)
        end subroutine daxpy

        subroutine dscal(n, a, x, incx)
            integer, intent(in) :: n, incx
            double precision, intent(in) :: a
            double precision, intent(in out) :: x(*)
        end subroutine dscal

        subroutine dcopy(n, x, incx, y, incy)
            integer, intent(in) :: n, incx,incy
            double precision, intent(in) :: x(*)
            double precision, intent(in out) :: y(*)
        end subroutine dcopy

        subroutine dswap(n, x, incx, y, incy)
            integer, intent(in) :: n, incx, incy
            double precision, intent(in out) :: x(:), y(:)
        end subroutine dswap

        subroutine dlaic1( job, j, x, sest, w, gamma, sestpr, s, c )
            integer, intent(in) :: job, j
            double precision, intent(in) :: x(j), w(j), gamma, sest
            double precision, intent(in out) :: sestpr, s, c
        end subroutine dlaic1

        subroutine funcEval(n, xcur, fcur, rpar, ipar, itrmf)
            ! Use solver and whatnot as global inputs, I'm sure as hell not combining all data into one rpar and then creating new routines on those!
            integer, intent(in) :: n
            integer, intent(in out) :: itrmf, ipar(*)
            double precision, intent(in) :: xcur(n)
            double precision, intent(in out) :: rpar(*), fcur(n)
        end subroutine funcEval
    
        subroutine jacv(n, xcur, fcur, ijob, v, z, rpar, ipar, itrmjv) 
            ! If analytical jacobian matrix-vector product or preconditioner needed
            integer, intent(in) :: ijob, n
            integer, intent(in out) :: itrmjv, ipar(*)
            double precision, intent(in) :: fcur(n), v(n), xcur(n)
            double precision, intent(in out) :: rpar(*), z(n)
        end subroutine jacv

    end interface

    private

    ! default parameters
    double precision, parameter :: one=1.0d0, two=2.0d0, rtfiv=2.23606797749978981, tenth=0.10d0, half=0.50d0, fournines=one-1.0d-4 
    double precision, parameter :: DFLT_CHOICE1_EXP=(one+rtfiv)*half, DFLT_CHOICE2_EXP=two, DFLT_CHOICE2_COEF=one, DFLT_ETA_CUTOFF=tenth, DFLT_ETA_MAX=fournines, &
    DFLT_THMIN=tenth, DFLT_THMAX=half, DFLT_ETA_FIXED=tenth
    integer, parameter :: DFLT_PRLVL=0, STDOUT=6

    ! shared variables
    integer :: iplvl, ipunit
    double precision :: choice1_exp, choice2_exp, choice2_coef, eta_cutoff, etamax, thmin, thmax, etafixed, epsmach, cndmax

    ! public variables
    integer, public :: instep, newstep, krystat
    double precision, public :: avrate, fcurnrm

contains

subroutine initializeNitsol()
    epsmach = 2.0d0*dlamch( 'e' )
    cndmax = 1.0/100.0d0/epsmach
end subroutine initializeNitsol

subroutine nitfd(n, xcur, fcur, rpar, ipar, ijacv, ifdord, nfe, v, z, rwork, itrmjv) 

     integer, intent(in) :: n, ijacv
     integer, intent(in out) :: ipar(*),ifdord, nfe, itrmjv
     double precision, intent(in) :: fcur(n)
     double precision, intent(in out) :: xcur(n), rpar(*), v(n), z(n), rwork(n) 


    ! c ------------------------------------------------------------------------
    ! c
    ! c This is nitfd v0.3, the routine for finite-difference evaluation of 
    ! c products z = J*v, where J is the Jacobian of f. 
    ! c
    ! c ------------------------------------------------------------------------
    ! c 
    ! c Explanation: 
    ! c
    ! c  n       = dimension of the problem.
    ! c
    ! c  xcur    = vector of length n, initial guess on input and final 
    ! c            approximate solution on output. 
    ! c
    ! c  fcur    = vector of length n, value of f at xcur. 
    ! c
    ! c  f       = name of user-supplied subroutine for evaluating the function 
    ! c            the zero of which is sought; this routine has the form 
    ! c
    ! c                  subroutine funcEval(n, xcur, fcur, rpar, ipar, itrmf)
    ! c
    ! c           where xcur is the array containing the current x value, fcur 
    ! c           is funcEval(xcur) on output, rpar and ipar are, respectively, real 
    ! c           and integer parameter/work arrays for use by the subroutine,
    ! c           and itrmf is an integer termination flag.  The meaning of
    ! c           itrmf is as follows:
    ! c             0 => normal termination; desired function value calculated.
    ! c             1 => failure to produce funcEval(xcur).
    ! c 
    ! c  rpar    = real parameter/work array passed to the f routine. 
    ! c
    ! c  ipar    = integer parameter/work array passed to the f routine. 
    ! c
    ! c  ijacv   = flag for determining method of J*v evaluation. In this 
    ! c            subroutine, this should be 0 or -1 on input, as follows: 
    ! c             -1 => finite-difference evaluation of order ifdord. 
    ! c              0 => first-order finite-difference evaluation. 
    ! c
    ! c  ifdord  = order of the finite-difference formula used when ijacv = -1.  
    ! c            When ijacv = 0 on input to nitsol, ifdord is set to 1, 2, or 
    ! c            4 in nitsol; subsequently, ijacv may be temporarily reset to 
    ! c            -1 in the Krylov solver to force a finite-difference evaluation 
    ! c            of order ifdord. If ijacv = 1 on input to nitsol, then ifdord 
    ! c            is irrelevant. When ijacv = 0 on input to nitsol, the precise 
    ! c            meaning of ifdord is as follows: 
    ! c
    ! c            If GMRES is used, then ifdord matters only if iresup = 1, in  
    ! c            which case it determines the order of the finite-difference 
    ! c            formula used in evaluating the initial residual at each GMRES 
    ! c            restart (default 2). NOTE: This only affects initial residuals 
    ! c            at restarts; first-order differences are always used within 
    ! c            each GMRES cycle. Using higher-order differences at restarts 
    ! c            only should give the same accuracy as if higher-order 
    ! c            differences were used throughout; see K. Turner and H. F. 
    ! c            Walker, "Efficient high accuracy solutions with GMRES(m)," 
    ! c            SIAM J. Sci. Stat. Comput., 13 (1992), pp. 815--825. 
    ! c               
    ! c            If BiCGSTAB or TFQMR is used, then ifdord determines the 
    ! c            order of the finite-difference formula used at each 
    ! c            iteration (default 1). 
    ! c
    ! c  nfe     = number of function evaluations.
    ! c
    ! c  v       = vector to be multiplied in the product J*v. 
    ! c
    ! c  z       = desired product J*v. 
    ! c
    ! c  rwork   = vector of length n, work array. 
    ! c
    ! c  dnrm2   = norm routine, either user-supplied or blas dnrm2. 
    ! c
    ! c  itrmjv  = termination flag; values have the following meanings: 
    ! c              0 => normal termination; desired product evaluated. 
    ! c              1 => failure to produce J*v.
    ! c
    ! c ------------------------------------------------------------------------
    ! c
    ! c Remark: In the selection of the difference step eps in the perturbation 
    ! c x + eps*v for approximating J*v by a finite difference, we assume the 
    ! c following: 
    ! c     1. f and a few derivatives (up to five when ifdord = 4) all 
    ! c        have about the same scale. 
    ! c     2. The relative error in f-evaluations is about epsmach (machine 
    ! c        epsilon). 
    ! c     3. The computed value of the sum of two vectors y and z has error 
    ! c        bounded by epsmach*(||y|| + ||z||). 
    ! c
    ! c The choice of eps is 
    ! c
    ! c             eps = {[(1 + ||x||)*epsmach]**(1/(ifdord+1)} /||v||, 
    ! c
    ! c which approximately minimizes a bound on the relative error in the 
    ! c difference approximation. 
    ! c
    ! c ------------------------------------------------------------------------
    ! c
    ! c Subroutines required by this and all called routines: 
    ! c
    ! c    user supplied: f 
    ! c
    ! c    nitsol routines: none 
    ! c
    ! c    blas routine: daxpy, dscal 
    ! c
    ! c    lapack routine:  dlamch

    ! c    user supplied or blas: dnrm2 
    ! c
    ! c    explanation: In nitsol, dnrm2 is set to either the blas 
    ! c    dnrm2 routine or the user-supplied usrnrm routine. 
    ! c
    ! c This subroutine called by: nitjv
    ! c
    ! c Subroutines called by this subroutine: daxpy, dlamch, dscal, dnrm2, f 
    ! c
    ! c Common block: 
    ! c
    !     
    ! c
    ! c If diagnostic information is desired, include this common block in the 
    ! c main program and set iplvl and ipunit according to the following: 
    ! c
    ! c     iplvl = 0 => no printout
    ! c           = 1 => iteration numbers and F-norms
    ! c           = 2 => ... + some stats, step norms, and linear model norms
    ! c           = 3 => ... + some Krylov solver and backtrack information
    ! c           = 4 => ... + more Krylov solver and backtrack information
    ! c
    ! c     ipunit = printout unit number.
    ! c
    ! c ------------------------------------------------------------------------
    ! c 
     double precision ::  eps, temp
     integer :: i, itrmf


    ! c ------------------------------------------------------------------------
    ! c Compute z = J*v by finite-differences: First, set eps = ||v||for later 
    ! c use in computing the difference step; then evaluate the difference 
    ! c formula according to ijacv and ifdord. 
    ! c ------------------------------------------------------------------------
     eps = dnrm2(n, v, 1)
     if (eps .eq. 0.d0) then 
        itrmjv = 1
        go to 900
     endif
    ! c ------------------------------------------------------------------------
    ! c Here ijacv = 0 or ifdord = 1 => first-order forward difference. 
    ! c ------------------------------------------------------------------------
     if (ijacv .eq. 0 .or. ifdord .eq. 1) then 
        eps = dsqrt((1.d0 + dnrm2(n,xcur,1))*epsmach)/eps
        ! could replace with just vectorized code?
        do 100 i = 1, n
           v(i) = xcur(i) + eps*v(i)
100     continue
        call funcEval(n, v, z, rpar, ipar, itrmf)
        if (itrmf .ne. 0) then
           itrmjv = 1
           goto 900
        endif 
        nfe = nfe + 1
        ! replace with vectorized code?
        do 110 i = 1, n
           z(i) = (z(i) - fcur(i))/eps
110     continue
        itrmjv = 0 
        go to 900
     endif
! c ------------------------------------------------------------------------
! c Here ijacv = -1 and ifdord = 2 => second-order central difference. 
! c ------------------------------------------------------------------------
     if (ifdord .eq. 2) then 
        eps = (((1.d0 + dnrm2(n,xcur,1))*epsmach)**(1.d0/3.d0))/eps
        do 200 i = 1, n
           rwork(i) = xcur(i) + eps*v(i)
200     continue
        call funcEval(n, rwork, z, rpar, ipar, itrmf) 
        nfe = nfe + 1
        if (itrmf .ne. 0) then
           itrmjv = 1
           goto 900
        endif 
        do 210 i = 1, n
           rwork(i) = xcur(i) - eps*v(i)
210     continue
        call funcEval(n, rwork, v, rpar, ipar, itrmf) 
        nfe = nfe + 1
        if (itrmf .ne. 0) then
           itrmjv = 1
           goto 900
        endif 
        temp = 2.d0*eps
        do 220 i = 1, n
           z(i) = (z(i) - v(i))/temp
220     continue
        itrmjv = 0 
        go to 900
     endif

     if (ifdord .eq. 4) then 
        eps = (((1.d0 + dnrm2(n,xcur,1))*epsmach)**(1.d0/5.d0))/eps
        do 300 i = 1, n
           rwork(i) = xcur(i) + eps*v(i)
300     continue
        call funcEval(n, rwork, z, rpar, ipar, itrmf) 
        nfe = nfe + 1
        if (itrmf .ne. 0) then
           itrmjv = 1
           goto 900
        endif 
        temp = -eps 
        call daxpy(n, temp, v, 1, xcur, 1)
        call funcEval(n, xcur, rwork, rpar, ipar, itrmf) 
        nfe = nfe + 1
        if (itrmf .ne. 0) then
           itrmjv = 1
           goto 900
        endif 
        do 310 i = 1, n
           z(i) = rwork(i) - z(i)
310     continue
        temp = eps/2.d0
        call daxpy(n, temp, v, 1, xcur, 1)
        call funcEval(n, xcur, rwork, rpar, ipar, itrmf) 
        nfe = nfe + 1
        if (itrmf .ne. 0) then
           itrmjv = 1
           goto 900
        endif 
        temp = -8.d0
        call daxpy(n, temp, rwork, 1, z, 1)
        call daxpy(n, eps, v, 1, xcur, 1)
        call funcEval(n, xcur, rwork, rpar, ipar, itrmf) 
        nfe = nfe + 1
        if (itrmf .ne. 0) then
           itrmjv = 1
           goto 900
        endif 
        temp = 8.d0
        call daxpy(n, temp, rwork, 1, z, 1)
        temp = 1.d0/(6.d0*eps)
        call dscal(n, temp, z, 1)
        temp = -eps/2.d0
        call daxpy(n, temp, v, 1, xcur, 1)
        itrmjv = 0 
     endif

900  continue
end subroutine nitfd

subroutine nitjv(n, xcur, fcur, rpar, ipar, ijacv, ifdord, itask, nfe, njve, nrpre, v, z, rwork1, rwork2, itrmjv) 

     integer, intent(in) :: n, ijacv, itask
     integer, intent(in out) :: ifdord, ipar(*), nfe, njve, nrpre, itrmjv
     double precision, intent(in) :: fcur(n)
     double precision, intent(in out) :: xcur(n), rpar(*), v(n), z(n), rwork1(n), rwork2(n)

    ! c ------------------------------------------------------------------------
    ! c
    ! c This is nitjv v0.3, the routine for controlling evaluation of products 
    ! c J*v or J*P(inverse)*v or P(inverse)*v, where J is the Jacobian of f 
    ! c and P is a right preconditioning operator. 
    ! c
    ! c ------------------------------------------------------------------------
    ! c 
    ! c Explanation: 
    ! c
    ! c  n       = dimension of the problem.
    ! c
    ! c  xcur    = vector of length n, initial guess on input and final 
    ! c            approximate solution on output. 
    ! c
    ! c  fcur    = vector of length n, value of f at xcur. 
    ! c
    ! c  f       = name of user-supplied subroutine for evaluating the function 
    ! c            the zero of which is sought; this routine has the form 
    ! c
    ! c                 subroutine f(n, xcur, fcur, rpar, ipar, itrmf)
    ! c
    ! c            where xcur is the array containing the current x value, fcur 
    ! c            is f(xcur) on output, rpar and ipar are, respectively, real 
    ! c            and integer parameter/work arrays for use by the subroutine,
    ! c            and itrmf is an integer termination flag.  The meaning of
    ! c            itrmf is as follows:
    ! c              0 => normal termination; desired function value calculated.
    ! c              1 => failure to produce f(xcur).
    ! c
    ! c  jacv    = name of user-supplied subroutine for evaluating J*v or 
    ! c            P(inverse)*v, where J is the Jacobian of f and P is a 
    ! c            right preconditioning operator. If neither analytic J*v 
    ! c            evaluations nor right preconditioning is used, this can 
    ! c            be a dummy subroutine; if right preconditioning is used but 
    ! c            not analytic J*v evaluations, this need only evaluate 
    ! c            P(inverse)*v. The form is 
    ! c
    ! c            subroutine jacv(n, xcur, fcur, ijob, v, z, rpar, ipar, itrmjv)
    ! c
    ! c            where xcur and fcur are vectors of length n containing the 
    ! c            current x and f values, ijob is an integer flag indicating 
    ! c            which product is desired, v is a vector of length n to be 
    ! c            multiplied, z is a vector of length n containing the desired 
    ! c            product on output, rpar and ipar are, respectively, real 
    ! c            and integer parameter/work arrays for use by the subroutine, 
    ! c            and itrmjv is an integer termination 
    ! c            flag. The meaning of ijob is as follows: 
    ! c              0 => z = J*v
    ! c              1 => z = P(inverse)*v 
    ! c            The meaning of itrmjv is as follows:
    ! c              0 => normal termination; desired product evaluated. 
    ! c              1 => failure to produce J*v.
    ! c              2 => failure to produce P(inverse)*v. 
    ! c            This subroutine is called only from nitjv, and is always 
    ! c            called with v .ne. z. 
    ! c
    ! c  rpar    = real parameter/work array passed to the f and jacv routines. 
    ! c
    ! c  ipar    = integer parameter/work array passed to the f and jacv routines. 
    ! c
    ! c  ijacv   = flag for determining method of J*v evaluation.
    ! c             -1 => finite-difference evaluation of order ifdord. 
    ! c              0 => first-order finite-difference evaluation. 
    ! c              1 => analytic evaluation. 
    ! c
    ! c  ifdord  = order of the finite-difference formula (sometimes) used when 
    ! c            J*v products are evaluated using finite-differences. When 
    ! c            ijacv = 0 on input to nitsol, ifdord is set to 1, 2, or 4 in 
    ! c            nitsol; subsequently, ijacv may be temporarily reset to -1 in 
    ! c            the Krylov solver to force a finite-difference evaluation of 
    ! c            order ifdord. If ijacv = 1 on input to nitsol, then ifdord is 
    ! c            irrelevant. When ijacv = 0 on input to nitsol, the precise 
    ! c            meaning of ifdord is as follows: 
    ! c
    ! c            If GMRES is used, then ifdord matters only if iresup = 1, in  
    ! c            which case it determines the order of the finite-difference 
    ! c            formula used in evaluating the initial residual at each GMRES 
    ! c            restart (default 2). NOTE: This only affects initial residuals 
    ! c            at restarts; first-order differences are always used within 
    ! c            each GMRES cycle. Using higher-order differences at restarts 
    ! c            only should give the same accuracy as if higher-order 
    ! c            differences were used throughout; see K. Turner and H. F. 
    ! c            Walker, "Efficient high accuracy solutions with GMRES(m)," 
    ! c            SIAM J. Sci. Stat. Comput., 13 (1992), pp. 815--825. 
    ! c               
    ! c            If BiCGSTAB or TFQMR is used, then ifdord determines the 
    ! c            order of the finite-difference formula used at each 
    ! c            iteration (default 1). 
    ! c
    ! c  itask   = flag for determining which product is produced.
    ! c              0 => z = J*v
    ! c              1 => z = J*P(inverse)*v 
    ! c              2 => z = P(inverse)*v 
    ! c
    ! c  nfe     = number of function evaluations.
    ! c
    ! c  njve    = number of J*v evaluations. 
    ! c
    ! c  nrpre   = number of P(inverse)*v evaluations.
    ! c
    ! c  v       = vector to be multiplied. 
    ! c
    ! c  z       = desired product. 
    ! c
    ! c  rwork1  = vector of length n, work array. 
    ! c
    ! c  rwork2  = vector of length n, work array. Note: rwork2 is only referenced 
    ! c            when a J-product is evaluated with a finite-difference of 
    ! c            order 2 or 4, i.e., when itask = 0 or 1, ijacv = 0, and 
    ! c            ifdord = 2 or 4. 
    ! c
    ! c  dnorm   = norm routine, either user-supplied or blas dnrm2. 
    ! c
    ! c  itrmjv  = termination flag; values have the following meanings: 
    ! c              0 => normal termination; desired product evaluated. 
    ! c              1 => failure to produce J*v.
    ! c              2 => failure to produce P(inverse)*v. 
    ! c
    ! c ------------------------------------------------------------------------
    ! c
    ! c Subroutines required by this and all called routines: 
    ! c
    ! c    user supplied: f, jacv 
    ! c
    ! c    nitsol routines: nitfd 
    ! c
    ! c    blas routine: daxpy, dcopy, dscal 
    ! c
    ! c    user supplied or blas: dnorm 
    ! c
    ! c    explanation: In nitsol, dnorm is set to either the blas 
    ! c    dnrm2 routine or the user-supplied usrnrm routine. 
    ! c
    ! c This subroutine called by: nitgm, nitstb, nittfq
    ! c
    ! c Subroutines called by this subroutine: dcopy, jacv, nitfd 
    ! c
    ! c Common block: 
    ! c
    !      include 'nitprint.h'
    ! c
    ! c If diagnostic information is desired, include this common block in the 
    ! c main program and set iplvl and ipunit according to the following: 
    ! c
    ! c     iplvl = 0 => no printout
    ! c           = 1 => iteration numbers and F-norms
    ! c           = 2 => ... + some stats, step norms, and linear model norms
    ! c           = 3 => ... + some Krylov solver and backtrack information
    ! c           = 4 => ... + more Krylov solver and backtrack information
    ! c
    ! c     ipunit = printout unit number.
    ! c
    ! c ------------------------------------------------------------------------

     integer :: ijob 
    ! c
    ! c ------------------------------------------------------------------------ 
    ! c
    ! c ------------------------------------------------------------------------
    ! c If z = J*v is desired (itask = 0), then copy v into rwork1; if 
    ! c z = J*P(inverse)*v or z = P(inverse)*v is desired (itask = 1,2), 
    ! c then compute P(inverse)*v in rwork1. 
    ! c ------------------------------------------------------------------------
     if (itask .eq. 0) then 
        call dcopy(n, v, 1, rwork1, 1)
     else
        ijob = 1
        call jacv(n, xcur, fcur, ijob, v, rwork1, rpar, ipar, itrmjv)
        nrpre = nrpre + 1
        if (itrmjv .ne. 0) go to 900
     endif
    ! c ------------------------------------------------------------------------
    ! c If only z = P(inverse)*v is desired (itask = 2), then copy rwork1 into 
    ! c z and exit.
    ! c ------------------------------------------------------------------------
     if (itask .eq. 2) then 
        call dcopy(n, rwork1, 1, z, 1)
        go to 900
     endif
    ! c ------------------------------------------------------------------------
    ! c If z = J*v or z = J*P(inverse)*v is desired (itask = 0, 1), then 
    ! c compute J*rwork1 in z by either analytic evaluation (ijacv = 1) or 
    ! c finite-differences (ijacv = 0, -1). 
    ! c ------------------------------------------------------------------------
     if (ijacv .eq. 1) then 
        ijob = 0
        call jacv(n, xcur, fcur, ijob, rwork1, z, rpar, ipar, itrmjv)
     else	
        call nitfd(n, xcur, fcur, rpar, ipar, ijacv, ifdord, nfe, rwork1, z, rwork2, itrmjv)
     endif
     njve = njve + 1
    ! c ------------------------------------------------------------------------
    ! c All returns made here.
    ! c ------------------------------------------------------------------------
900  continue

end subroutine nitjv


subroutine nitgm (n, xcur, fcur, fcnrm, step, eta, rpar, &
        ipar, ijacv, irpre, iksmax, iresup, ifdord, nfe, njve, & 
        nrpre, nli, kdmax, kdmaxp1, vv, rr, svbig, svsml, w, rwork, &
        rsnrm, itrmks)

     integer :: n, ipar(*), ijacv, irpre, iksmax, iresup, ifdord, nfe, njve, nrpre, nli, kdmax, kdmaxp1, itrmks
     double precision, intent(in out) :: xcur(n), fcur(n), fcnrm, step(n), eta, rpar(*), &
        vv(n,kdmaxp1), rr(kdmax,kdmax), svbig(kdmax), svsml(kdmax), &
        w(kdmax), rwork(n), rsnrm

! c ------------------------------------------------------------------------
! c
! c This is nitgm v0.3, the GMRES routine for determining (trial) inexact 
! c Newton steps. This implementation is the "simpler" Gram-Schmidt GMRES 
! c implementation from L. Zhou and H. F. Walker, "A simpler GMRES," 
! c J. Numerical Lin. Alg. Appl., 1 (1994), pp. 571-581. 
! c
! c ------------------------------------------------------------------------
! c 	
! c Explanation: 
! c
! c  n       = dimension of the problem.
! c
! c  xcur    = vector of length n, current approximate solution. 
! c
! c  fcur    = vector of length n, value of f at xcur. 
! c
! c  fcnrm   = norm of fcur. 
! c
! c  step    = vector of length n, (trial) step. 
! c
! c  eta     = relative residual reduction factor. 
! c
! c  f       = name of user-supplied subroutine for evaluating the function 
! c            the zero of which is sought; this routine has the form 
! c
! c                  subroutine f(n, xcur, fcur, rpar, ipar, itrmf)
! c
! c           where xcur is the array containing the current x value, fcur 
! c           is f(xcur) on output, rpar and ipar are, respectively, real 
! c           and integer parameter/work arrays for use by the subroutine,
! c           and itrmf is an integer termination flag.  The meaning of
! c           itrmf is as follows:
! c             0 => normal termination; desired function value calculated.
! c             1 => failure to produce f(xcur).
! c 
! c  jacv    = name of user-supplied subroutine for evaluating J*v or 
! c            P(inverse)*v, where J is the Jacobian of f and P is a 
! c            right preconditioning operator. If neither analytic J*v 
! c            evaluations nor right preconditioning is used, this can 
! c            be a dummy subroutine; if right preconditioning is used but 
! c            not analytic J*v evaluations, this need only evaluate 
! c            P(inverse)*v. The form is 
! c
! c            subroutine jacv(n, xcur, fcur, ijob, v, z, rpar, ipar, itrmjv)
! c
! c            where xcur and fcur are vectors of length n containing the 
! c            current x and f values, ijob is an integer flag indicating 
! c            which product is desired, v is a vector of length n to be 
! c            multiplied, z is a vector of length n containing the desired 
! c            product on output, rpar and ipar are, respectively, real 
! c            and integer parameter/work arrays for use by the subroutine, 
! c            and itrmjv is an integer termination 
! c            flag. The meaning of ijob is as follows: 
! c              0 => z = J*v
! c              1 => z = P(inverse)*v 
! c            The meaning of itrmjv is as follows:
! c              0 => normal termination; desired product evaluated. 
! c              1 => failure to produce J*v.
! c              2 => failure to produce P(inverse)*v. 
! c            This subroutine is called only from nitjv, and is always 
! c            called with v .ne. z. 
! c
! c  rpar    = real parameter/work array passed to the f and jacv routines. 
! c
! c  ipar    = integer parameter/work array passed to the f and jacv routines. 
! c
! c  ijacv   = flag for determining method of J*v evaluation.
! c              0 => finite-difference evaluation (default). 
! c              1 => analytic evaluation. 
! c
! c  irpre   = flag for right preconditioning. 
! c              0 => no right preconditioning
! c              1 => right preconditioning
! c
! c  iksmax  = maximum allowable number of GMRES iterations. 
! c
! c  iresup  = residual update flag; on GMRES restarts, the residual is 
! c            updated as follows: 
! c               0 => linear combination (default) 
! c               1 => direct evaluation
! c            The first is cheap (one n-vector saxpy) but may lose 
! c            accuracy with extreme residual reduction; the second 
! c            retains accuracy better but costs one J*v product 
! c
! c  ifdord  = order of the finite-difference formula (sometimes) used on 
! c            GMRES restarts when J*v products are evaluated using finite-
! c            differences. When ijacv = 0 on input to nitsol, ifdord is set 
! c            to 1, 2, or 4 in nitsol; otherwise, it is irrelevant. When 
! c            ijacv = 0 on input to this subroutine, the precise meaning is 
! c            as follows: 
! c
! c            With GMRES, ifdord matters only if iresup = 1, in which case 
! c            it determines the order of the finite-difference formula used 
! c            in evaluating the initial residual at each GMRES restart 
! c            (default 2). If iresup = 1 and ijacv = 0 on input to this 
! c            subroutine, then ijacv is temporarily reset to -1 at each 
! c            restart below to force a finite-difference evaluation of order 
! c            ifdord. NOTE: This only affects initial residuals at restarts; 
! c            first-order differences are always used within each GMRES 
! c            cycle. Using higher-order differences at restarts only should 
! c            give the same accuracy as if higher-order differences were 
! c            used throughout; see K. Turner and H. F. Walker, "Efficient 
! c            high accuracy solutions with GMRES(m)," SIAM J. Sci. Stat. 
! c            Comput., 13 (1992), pp. 815--825. 
! c               
! c  nfe     = number of function evaluations.
! c
! c  njve    = number of J*v evaluations. 
! c
! c  nrpre   = number of P(inverse)*v evaluations.
! c
! c  nli     = number of linear iterations.
! c
! c  kdmax   = maximum Krylov subspace dimension; default 10. 
! c
! c  kdmaxp1 = kdmax + 1. 
! c
! c  vv      = n x (kdmax+1) matrix for storage of Krylov basis in GMRES;
! c            on return, the residual vector is contained in the first 
! c            column.
! c
! c  rr      = kdmax x kdmax matrix for storage of triangular matrix in GMRES. 
! c
! c  svbig   = vector of length kdmax for storage of estimate of singular 
! c            vector of rr with largest singular value. 
! c
! c  svsml   = vector of length kdmax for storage of estimate of singular 
! c            vector of rr with smallest singular value. 
! c
! c  w       = vector of length kdmax, contains right-hand side of 
! c            triangular system and least-squares residual norm in GMRES. 
! c
! c  rwork   = vector of length n, work array. 
! c
! c  rsnrm   = GMRES residual norm on return. 
! c
! c  dinpr   = inner-product routine, either user-supplied or blas ddot. 
! c
! c  dnorm   = norm routine, either user-supplied or blas dnrm2. 
! c
! c  itrmks  = termination flag; values have the following meanings: 
! c              0 => normal termination: acceptable step found. 
! c              1 => J*v failure in nitjv. 
! c              2 => P(inverse)*v failure in nitjv. 
! c              3 => acceptable step not found in iksmax GMRES iterations. 
! c              4 => insufficient residual norm reduction over a cycle 
! c                   of kdmax steps (stagnation) before an acceptable step 
! c                   has been found. 
! c              5 => dangerous ill-conditioning detected before an acceptable 
! c                   step has been found. 
! c
! c             Note: On return, nitsol terminates if itrmks is 1 or 2. If  
! c             itrmks is 3, 4, or 5, nitsol may terminate or continue. In 
! c             this event, a meaningful inexact Newton step is returned, 
! c             even though the desired inexact Newton condition may not 
! c             hold, and a decision on termination/continuation is made 
! c             in nitdrv according to whether there is sufficient residual 
! c             norm reduction.
! c
! c -------------------------------------------------------------------------
! c
! c Subroutines required by this and all called routines: 
! c
! c    user supplied: f, jacv 
! c
! c    nitsol routines: nitjv, nitfd
! c
! c    lapack routine: dlaic1, dlamch
! c
! c    blas routines: daxpy, dcopy, dscal
! c
! c    user supplied or blas: ddot, dnorm 
! c
! c    explanation: In nitsol, ddot and dnorm are set to either the blas 
! c    ddot and dnrm2 routines or the user-supplied usrnpr and usrnrm 
! c    routines. 
! c
! c This subroutine called by: nitdrv
! c
! c Subroutines called by this subroutine: daxpy, dcopy, dscal, dlaic1, 
! c    dlamch, ddot, dnorm, nitjv
! c
! c Common block: 
! c
!      include 'nitprint.h'
! c
! c If diagnostic information is desired, include this common block in the 
! c main program and set iplvl and ipunit according to the following: 
! c
! c     iplvl = 0 => no printout
! c           = 1 => iteration numbers and F-norms
! c           = 2 => ... + some stats, step norms, and linear model norms
! c           = 3 => ... + some Krylov solver and backtrack information
! c           = 4 => ... + more Krylov solver and backtrack information
! c
! c     ipunit = printout unit number.
! c
! c ------------------------------------------------------------------------
! c 
     double precision :: abstol, big, cs, rsnrm0, sestpr, small, sn, temp 
     integer :: i, igm, ijob, itask, itrmjv, kd, kdp1



! c ------------------------------------------------------------------------
! c Initialize. 
! c ------------------------------------------------------------------------
     do 20 i = 1, n
        step(i) = 0.d0
20   continue
     igm = 0
! c ------------------------------------------------------------------------ 
! c For printing:
     if (iplvl .ge. 3) then 
        write(ipunit,*) 
        write(ipunit,800) eta 
     endif
800  format('nitgm:  eta =', 1pd10.3)
     if (iplvl .ge. 4) then 
        write(ipunit,810) 
        write(ipunit,*) 
        write(ipunit,820) igm, fcnrm 
     endif
810  format('nitgm:  GMRES iteration no., linear residual norm, ','condition no. estimate')
820  format(5x,i4,2(5x,1pd10.3))
! c ------------------------------------------------------------------------
! c
! c ------------------------------------------------------------------------
! c Set the stopping tolerance, etc. 
! c ------------------------------------------------------------------------
     rsnrm0 = fcnrm
     abstol = eta*rsnrm0
! c ------------------------------------------------------------------------
! c Place the normalized initial residual in the first column of the vv array.
! c ------------------------------------------------------------------------
     call dcopy(n, fcur, 1, vv(1,1), 1)
     temp = -1.d0/fcnrm
     call dscal(n, temp, vv(1,1), 1)
! c ------------------------------------------------------------------------
! c Top of the outer GMRES loop. 
! c ------------------------------------------------------------------------
100  continue
     kd = 0
     rsnrm = 1.d0
! c ------------------------------------------------------------------------
! c Top of the inner GMRES loop.
! c ------------------------------------------------------------------------
200  continue
     kd = kd + 1
     kdp1 = kd + 1
     nli = nli + 1
     igm = igm + 1
! c ------------------------------------------------------------------------
! c Evaluate J*(kd-th Krylov subspace basis vector) in vv(.,kdp1). 
! c Note: rwork can be used for both work arrays in this call because 
! c the second is not referenced within nitjv. 
! c ------------------------------------------------------------------------
     if (irpre .eq. 0) then 
        itask = 0
     else
        itask = 1
     endif
     call nitjv(n, xcur, fcur, rpar, ipar, ijacv, ifdord, itask, nfe, njve, nrpre, vv(1,kd), vv(1,kdp1), rwork, rwork, itrmjv)
     if (itrmjv .gt. 0) then 
        if (itrmjv .eq. 1) itrmks = 1
        if (itrmjv .eq. 2) itrmks = 2
        go to 900
     endif
! c ------------------------------------------------------------------------
! c Do modified Gram-Schmidt. 
! c ------------------------------------------------------------------------
     do 210 i = 2, kd
        rr(i-1,kd) = ddot(n, vv(1,i), 1, vv(1,kdp1), 1)
        call daxpy(n, -rr(i-1,kd), vv(1,i), 1, vv(1,kdp1), 1)
210  continue
     rr(kd,kd) = dnrm2(n, vv(1,kdp1), 1)
! c ------------------------------------------------------------------------
! c Update the estimates of the largest and smallest singular values. 
! c ------------------------------------------------------------------------
     if (kd .eq. 1) then
        big = rr(1,1)
        small = big
        svbig(1) = 1.d0
        svsml(1) = 1.d0
     else
        ijob = 1
        call dlaic1(ijob, kd-1, svbig, big, rr(1,kd), rr(kd,kd),sestpr, sn, cs)
        big = sestpr
        call dscal(kd-1, sn, svbig, 1)
        svbig(kd) = cs
        ijob = 2
        call dlaic1(ijob, kd-1, svsml, small, rr(1,kd), rr(kd,kd),sestpr, sn, cs)
        small = sestpr
        call dscal(kd-1, sn, svsml, 1)
        svsml(kd) = cs
     endif
! c ------------------------------------------------------------------------
! c Terminate if the estimated condition number is too great. 
! c ------------------------------------------------------------------------
     if (big .ge. small*cndmax) then 
        if (kd .eq. 1) then 
           itrmks = 5
           go to 900
        else 
           kdp1 = kd
           kd = kd - 1
           call daxpy(n, w(kd), vv(1,kdp1), 1, vv(1,1), 1)
           go to 300
        endif
     endif
! c ------------------------------------------------------------------------
! c Normalize vv(.,kdp1). 
! c ------------------------------------------------------------------------
     temp = 1.d0/rr(kd,kd) 
     call dscal(n, temp, vv(1,kdp1), 1)
! c ------------------------------------------------------------------------
! c Update w and the residual norm by rsnrm <- rsnrm*dsin(dacos(w(kd)/rsnrm). 
! c ------------------------------------------------------------------------
     w(kd) = ddot(n, vv(1,1), 1, vv(1,kdp1), 1)
     temp = max(min(w(kd)/rsnrm,1.0D0),-1.0d0) 
     rsnrm = rsnrm*dsin(dacos(temp))
! c ------------------------------------------------------------------------ 
! c For printing:
     if (iplvl .ge. 4) then 
        write(ipunit,820) igm, rsnrm*rsnrm0, big/small
! c         if (kd .eq. kdmax) write(ipunit,*) 
     endif
! c ------------------------------------------------------------------------
! c
! c ------------------------------------------------------------------------
! c Test for termination of the inner loop.
! c ------------------------------------------------------------------------
     if ( (rsnrm0*rsnrm .le. abstol) .or. (kd .eq. kdmax) .or.  (igm .ge. iksmax) )  go to 300
! c ------------------------------------------------------------------------
! c If not terminating the inner loop, update the residual vector 
! c and go to the top of the inner loop. 
! c ------------------------------------------------------------------------
     call daxpy(n, -w(kd), vv(1,kdp1), 1, vv(1,1), 1)
     go to 200
! c ------------------------------------------------------------------------
! c Bottom of inner loop.
! c ------------------------------------------------------------------------
300  continue
! c ------------------------------------------------------------------------ 
! c For printing:
     if (iplvl .ge. 4) then 
        write(ipunit,*) 
     endif
! c ------------------------------------------------------------------------
! c
! c ------------------------------------------------------------------------
! c Compute the solution: 
! c ------------------------------------------------------------------------
! c
! c ------------------------------------------------------------------------
! c Use svbig for storage of the original components of w. 
! c ------------------------------------------------------------------------
     call dcopy(kd, w, 1, svbig, 1)
! c ------------------------------------------------------------------------
! c Overwrite w with the solution of the upper triangular system.
! c ------------------------------------------------------------------------
     do 310 i = kd, 1, -1
        w(i) = w(i)/rr(i,i)
        if (i .gt. 1) call daxpy(i-1, -w(i), rr(1,i), 1, w, 1)
310  continue
! c ------------------------------------------------------------------------
! c Now form the linear combination to accumulate the correction in 
! c the work vector.
! c ------------------------------------------------------------------------
     call dcopy(n, vv(1,1), 1, rwork, 1)
     call dscal(n, w(1), rwork, 1)
     if (kd .gt. 1) then 
        call daxpy(kd-1, w(1), svbig, 1, w(2), 1)
        do 320 i = 2, kd
           call daxpy(n, w(i), vv(1,i), 1, rwork, 1)
320     continue
     endif
! c ------------------------------------------------------------------------
! c If iresup .eq. 0, then update the residual vector by linear 
! c combination. This frees vv(.,kdp1) for use as a work array. 
! c ------------------------------------------------------------------------
     if (iresup .eq. 0) then 
        call daxpy(n, -svbig(kd), vv(1,kdp1), 1, vv(1,1), 1)
     endif
! c ------------------------------------------------------------------------
! c If right preconditioning is used, overwrite 
! c correction <-- P(inverse)*correction, using vv(.,kdp1) as a work array. 
! c Note: vv(.,kdp1) can be used for both work arrays in this call because 
! c the second is not referenced within nitjv. 
! c ------------------------------------------------------------------------
     if (irpre .gt. 0) then 
        itask = 2
        call nitjv(n, xcur, fcur, rpar, ipar, ijacv, ifdord, itask, nfe, njve, nrpre, rwork, rwork, vv(1,kdp1), vv(1,kdp1), itrmjv)
        if (itrmjv .gt. 0) then 
           itrmks = 2
           go to 900
        endif
     endif
! c ------------------------------------------------------------------------
! c Update the step. This frees rwork for use as a work array.
! c ------------------------------------------------------------------------
     call daxpy(n, rsnrm0, rwork, 1, step, 1)
! c ------------------------------------------------------------------------
! c If iresup .eq. 1, then update the residual vector by direct evaluation, 
! c using rwork and vv(.,kdp1) as work arrays. Note: Two distinct work  
! c arrays are needed in this call because both are referenced within nitjv 
! c if the J*step product is evaluated with a finite-difference of order 
! c two or higher. If finite-differences are used (ijacv= 0), then ijacv 
! c is temporarily set to -1 to signal to nitjv that the order of the 
! c finite-difference formula is to be determined by ifdord. 
! c ------------------------------------------------------------------------
     if (iresup .eq. 1) then 
        itask = 0
        if (ijacv .eq. 0) ijacv = -1 
        call nitjv(n, xcur, fcur, rpar, ipar, ijacv, ifdord, itask, nfe, njve, nrpre, step, vv(1,1), rwork, vv(1,kdp1), itrmjv)
        if (ijacv .eq. -1) ijacv = 0
        if (itrmjv .gt. 0) then 
           itrmks = 1
           go to 900
        endif
        do 330 i = 1, n
           vv(i,1) = -fcur(i) - vv(i,1)
330     continue
     endif
! c ------------------------------------------------------------------------
! c Test for termination.  
! c ------------------------------------------------------------------------
     if (rsnrm0*rsnrm .le. abstol) then 
        itrmks = 0
        go to 900
     endif
     if (igm .ge. iksmax) then 
        itrmks = 3
        go to 900
     endif
     if (big .ge. small*cndmax) then 
        itrmks = 5
        go to 900
     endif
     temp = dfloat(kd)*dlog(abstol/(rsnrm0*rsnrm))/dlog(rsnrm/(1.d0 + 10.d0*epsmach))
     if (temp .ge. 1000.d0*dfloat(iksmax - igm)) then 
        itrmks = 4
        go to 900
     endif
! c ------------------------------------------------------------------------
! c If not terminating, then normalize the initial residual, etc., and 
! c return to the top of the outer loop. 
! c ------------------------------------------------------------------------
     if (iresup .eq. 0) then 
        rsnrm0 = rsnrm0*rsnrm
        temp = 1.d0/rsnrm
     else
        rsnrm0 = dnrm2(n, vv(1,1), 1)
        temp = 1.d0/rsnrm0
     endif
     call dscal(n, temp, vv(1,1), 1)
     go to 100
! c ------------------------------------------------------------------------
! c All returns made here.
! c ------------------------------------------------------------------------
900  continue
     if (itrmks .ne. 1 .and. itrmks .ne. 2) then 
        if (iresup .eq. 0) then 
           call dscal(n, rsnrm0, vv(1,1), 1)
           rsnrm = rsnrm0*rsnrm
        else
           rsnrm = dnrm2(n, vv(1,1), 1)
        endif
     endif
! c ------------------------------------------------------------------------ 
! c For printing:
     if (iplvl .ge. 3) then 
        if (itrmks .ne. 1 .and. itrmks .ne. 2) then 
           write(ipunit,830) itrmks, rsnrm 
830        format('nitgm:  itrmks =', i2, '    final lin. res. norm =', 1pd10.3)
        else
           write(ipunit,840) itrmks
840        format('nitgm: itrmks:', i4) 
        endif
     endif

end subroutine nitgm

end module mod_nitsol