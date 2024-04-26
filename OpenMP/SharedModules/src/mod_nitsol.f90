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

    end interface

   abstract interface

      subroutine funcInterface(n, xcur, fcur, rpar, ipar, itrmf)
         ! Use solver and whatnot as global inputs, I'm sure as hell not combining all data into one rpar and then creating new routines on those!
         integer, intent(in) :: n
         integer, intent(in out) :: itrmf, ipar(*)
         double precision, intent(in) :: xcur(n)
         double precision, intent(in out) :: rpar(*), fcur(n)
      end subroutine funcInterface

      subroutine jacvInterface(n, xcur, fcur, ijob, v, z, rpar, ipar, itrmjv) 
            ! If analytical jacobian matrix-vector product or preconditioner needed
            integer, intent(in) :: ijob, n
            integer, intent(in out) :: itrmjv, ipar(*)
            double precision, intent(in) :: fcur(n), v(n), xcur(n)
            double precision, intent(in out) :: rpar(*), z(n)
      end subroutine jacvInterface

   end interface

    private
    public :: nitsol, initializeNitsol

    ! default parameters
    double precision, parameter :: one=1.0d0, two=2.0d0, rtfiv=2.23606797749978981, tenth=0.10d0, half=0.50d0, fournines=one-1.0d-4, zero = 0.0d0 
    double precision, parameter :: DFLT_CHOICE1_EXP=(one+rtfiv)*half, DFLT_CHOICE2_EXP=two, DFLT_CHOICE2_COEF=one, DFLT_ETA_CUTOFF=tenth, DFLT_ETA_MAX=fournines, &
    DFLT_THMIN=tenth, DFLT_THMAX=half, DFLT_ETA_FIXED=tenth
    integer, parameter :: DFLT_PRLVL=0, STDOUT=6

    ! shared variables
    integer :: iplvl, ipunit
    double precision :: choice1_exp, choice2_exp, choice2_coef, eta_cutoff, etamax, thmin, thmax, etafixed, epsmach, cndmax, sfmin
    double precision, allocatable :: rworkNitsol(:)

    logical :: nitsolInitBool = .false.

    ! Set parameters used in driver routine
    integer :: nnimax_glob, ijacv_glob, ifdord_glob, ikrysl_glob, kdmax_glob, irpre_glob, iksmax_glob, iresup_glob, ibtmax_glob, ieta_glob


    ! public variables
    integer, public :: instep, newstep, krystat
    double precision, public :: avrate, fcurnrm

contains

subroutine initializeNitsol(maxIter, m_Anderson, n_size)
   integer, intent(in) :: maxIter, m_Anderson, n_size
   integer :: input(10), sizeWorkArray
    ! initialize variables upon first call
    epsmach = 2.0d0*dlamch( 'e' )
    sfmin = dlamch('s')
    cndmax = 1.0/100.0d0/epsmach
    iplvl = DFLT_CHOICE1_EXP
    ipunit = STDOUT
    choice1_exp = DFLT_CHOICE1_EXP
    choice2_exp = DFLT_CHOICE2_EXP
    choice2_coef = DFLT_CHOICE2_COEF
    eta_cutoff = DFLT_ETA_CUTOFF
    etamax = DFLT_ETA_MAX
    thmin = DFLT_THMIN
    thmax = DFLT_THMAX
    etafixed = DFLT_ETA_FIXED

   iplvl = 1 ! print level
   input = 0
   input(1) = maxIter ! maximum outer newton iterations
   input(2) = 0 !ijacv
   input(3) = 0 ! krylov solver (0 == GMRES)
   input(4) = m_Anderson ! maximum krylov subspace dimension
   input(5) = 0 !ipre
   input(9) = m_Anderson !number backtracks
   input(6) = m_Anderson*2 ! maximum iterations to kyrlov routine
   input(10) = 2 ! eta with gamma and alpha
   etamax = 0.8d0 ! eta max
   choice2_exp = 1.5d0 ! alpha
   choice2_coef = 0.9d0 ! gamma 

   ! c ------------------------------------------------------------------------ 
   ! c
   ! c Begin executable code. 
   ! c 
   ! c ------------------------------------------------------------------------
   ! c Check inputs and initialize parameters. 
   ! c ------------------------------------------------------------------------
   if (ipunit .gt. 6) open( unit=ipunit, status='unknown' )
   if (ipunit .eq. 0) ipunit = 6
   if (input(1) .eq. 0) then
      nnimax_glob = 200
   else
      if (input(1) .gt. 0) then 
         nnimax_glob = input(1)
      else
         stop "negative max iteration for nitsol"
      endif
   endif
   if (input(2) .eq. 0 .or. input(2) .eq. 1) then 
      ijacv_glob = input(2) 
   else
      stop "ijacv should be 0 or 1 in nitsol"
   endif
   if (input(3) .ge. 0 .and. input(3) .le. 2) then 
      ikrysl_glob = input(3)
   else 
      stop "input(3) in nitsol not correct"
   endif
   if (ikrysl_glob .eq. 0) then 
      if (input(4) .eq. 0) then 
         kdmax_glob = 20
      else
         if (input(4) .gt. 0) then 
            kdmax_glob = input(4) 
         else
            stop "input(4) in nitsol not correct"
         endif
      endif
   endif
   if (input(5) .eq. 0 .or. input(5) .eq. 1) then 
      irpre_glob = input(5) 
   else
      stop "input(5) in nitsol not correct"
   endif
   if (input(6) .eq. 0) then 
      iksmax_glob = 1000
   else
      if (input(6) .gt. 0) then 
         iksmax_glob = input(6)
      else
         stop "input(6) in nitsol incorrect"
      endif
   endif
   if (ikrysl_glob .eq. 0) then 
      if (input(7) .eq. 0 .or. input(7) .eq. 1) then 
         iresup_glob = input(7)
      else
         stop "input(7) in nitsol incorrect"
      endif
   endif
   if (ijacv_glob .eq. 0) then 
      if (input(8) .eq. 0) then 
         if (ikrysl_glob .eq. 0) then 
            ifdord_glob = 2
         else
            ifdord_glob = 1
         endif
      else
         if (input(8) .eq. 1 .or. input(8) .eq. 2 .or. input(8) .eq. 4) then 
            ifdord_glob = input(8)
         else
            stop "input(8) in nitsol not correct"
         endif
      endif
   endif
   if (input(9) .eq. 0) then 
      ibtmax_glob = 10
   else
      if (input(9) .gt. 0 .or. input(9) .eq. -1) then 
         ibtmax_glob = input(9)
      else
         stop "input(9) in nitsol not correct"
      endif
   endif
   if (input(10) .ge. 0 .and. input(10) .le. 3) then 
      ieta_glob = input(10)
   else
      stop "input(10) in nitsol not correct"
   endif
! c ------------------------------------------------------------------------
! c  Check possible invalid value for printout level.  In
! c  case the value is invalid the default is restored.
! c ------------------------------------------------------------------------
   if ( iplvl .lt. 0 .or. iplvl .gt. 4 ) iplvl = DFLT_PRLVL
! c ------------------------------------------------------------------------
! c  Check possible invalid values for various parameters.  In
! c  case the values are invalid the defaults are restored.
! c ------------------------------------------------------------------------
   if ( choice1_exp .le. 1.0d0 .or. choice1_exp .gt. 2.0d0 ) choice1_exp = DFLT_CHOICE1_EXP
   if ( choice2_exp .le. 1.0d0 .or. choice2_exp .gt. 2.0d0 ) choice2_exp = DFLT_CHOICE2_EXP
   if ( choice2_coef .lt. 0.0d0 .or. choice2_coef .gt. 1.0d0 ) choice2_coef = DFLT_CHOICE2_COEF
   if ( etamax .le. 0.0d0 ) etamax = DFLT_ETA_MAX
   if ( thmin .lt. 0.0d0 .or. thmin .gt. thmax ) thmin = DFLT_THMIN
   if ( thmax .gt. 1.0d0 .or. thmax .lt. thmin ) thmax = DFLT_THMAX
   if (ikrysl_glob == 0) then
      sizeWorkArray = n_size*(kdmax_glob+5)+kdmax_glob*(kdmax_glob+3)
   else if (ikrysl_glob == 1) then
      sizeWorkArray = n_size*11
   else if (ikrysl_glob == 2) then
      sizeWorkArray = n_size * 14
   end if
   allocate(rworkNitsol(sizeWorkArray))
   nitsolInitBool = .true.
end subroutine initializeNitsol

subroutine nitfd(n, xcur, fcur, f, rpar, ipar, nfe, v, z, rwork, itrmjv) 

     integer, intent(in) :: n
     integer, intent(in out) :: ipar(*), nfe, itrmjv
     double precision, intent(in) :: fcur(n)
     double precision, intent(in out) :: xcur(n), rpar(*), v(n), z(n), rwork(n) 
     procedure(funcInterface) :: f


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
     eps = SQRT(SUM(v**2))
     if (eps .eq. 0.d0) then 
        itrmjv = 1
        go to 900
     endif
    ! c ------------------------------------------------------------------------
    ! c Here ijacv = 0 or ifdord = 1 => first-order forward difference. 
    ! c ------------------------------------------------------------------------
     if (ijacv_glob .eq. 0 .or. ifdord_glob .eq. 1) then 
        eps = dsqrt((1.d0 + SQRT(SUM(xcur**2)))*epsmach)/eps
        ! could replace with just vectorized code?
        do 100 i = 1, n
           v(i) = xcur(i) + eps*v(i)
    100     continue
        call f(n, v, z, rpar, ipar, itrmf)
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
     if (ifdord_glob .eq. 2) then 
        eps = (((1.d0 + SQRT(SUM(xcur**2)))*epsmach)**(1.d0/3.d0))/eps
        do 200 i = 1, n
           rwork(i) = xcur(i) + eps*v(i)
    200     continue
        call f(n, rwork, z, rpar, ipar, itrmf) 
        nfe = nfe + 1
        if (itrmf .ne. 0) then
           itrmjv = 1
           goto 900
        endif 
        do 210 i = 1, n
           rwork(i) = xcur(i) - eps*v(i)
    210     continue
        call f(n, rwork, v, rpar, ipar, itrmf) 
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

     if (ifdord_glob .eq. 4) then 
        eps = (((1.d0 + SQRT(SUM(xcur**2)))*epsmach)**(1.d0/5.d0))/eps
        do 300 i = 1, n
           rwork(i) = xcur(i) + eps*v(i)
    300     continue
        call f(n, rwork, z, rpar, ipar, itrmf) 
        nfe = nfe + 1
        if (itrmf .ne. 0) then
           itrmjv = 1
           goto 900
        endif 
        temp = -eps 
        xcur = xcur + temp*v!call daxpy(n, temp, v, 1, xcur, 1)
        call f(n, xcur, rwork, rpar, ipar, itrmf) 
        nfe = nfe + 1
        if (itrmf .ne. 0) then
           itrmjv = 1
           goto 900
        endif 
        do 310 i = 1, n
           z(i) = rwork(i) - z(i)
    310     continue
        temp = eps/2.d0
        xcur = xcur + temp*v!call daxpy(n, temp, v, 1, xcur, 1)
        call f(n, xcur, rwork, rpar, ipar, itrmf) 
        nfe = nfe + 1
        if (itrmf .ne. 0) then
           itrmjv = 1
           goto 900
        endif 
        temp = -8.d0
        z = z + temp * rwork(1:n) !call daxpy(n, temp, rwork, 1, z, 1)
        xcur = xcur + eps * v!call daxpy(n, eps, v, 1, xcur, 1)
        call f(n, xcur, rwork, rpar, ipar, itrmf) 
        nfe = nfe + 1
        if (itrmf .ne. 0) then
           itrmjv = 1
           goto 900
        endif 
        temp = 8.d0
        z = z + temp * rwork(1:n)!call daxpy(n, temp, rwork, 1, z, 1)
        temp = 1.d0/(6.d0*eps)
        z = temp * z!call dscal(n, temp, z, 1)
        temp = -eps/2.d0
        xcur = xcur + temp*v!call daxpy(n, temp, v, 1, xcur, 1)
        itrmjv = 0 
     endif

    900  continue
end subroutine nitfd

subroutine nitjv(n, xcur, fcur, f, jacv, rpar, ipar, itask, nfe, njve, nrpre, v, z, rwork1, rwork2, itrmjv) 

     integer, intent(in) :: n, itask
     integer, intent(in out) :: ipar(*), nfe, njve, nrpre, itrmjv
     double precision, intent(in) :: fcur(n)
     double precision, intent(in out) :: xcur(n), rpar(*), v(n), z(n), rwork1(n), rwork2(n)
     procedure(funcInterface) :: f
     procedure(jacvInterface) :: jacv

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
        rwork1(1:n) = v!call dcopy(n, v, 1, rwork1, 1)
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
        z = rwork1(1:n)!call dcopy(n, rwork1, 1, z, 1)
        go to 900
     endif
    ! c ------------------------------------------------------------------------
    ! c If z = J*v or z = J*P(inverse)*v is desired (itask = 0, 1), then 
    ! c compute J*rwork1 in z by either analytic evaluation (ijacv = 1) or 
    ! c finite-differences (ijacv = 0, -1). 
    ! c ------------------------------------------------------------------------
     if (ijacv_glob .eq. 1) then 
        ijob = 0
        call jacv(n, xcur, fcur, ijob, rwork1, z, rpar, ipar, itrmjv)
     else	
        call nitfd(n, xcur, fcur, f, rpar, ipar, nfe, rwork1, z, rwork2, itrmjv)
     endif
     njve = njve + 1
    ! c ------------------------------------------------------------------------
    ! c All returns made here.
    ! c ------------------------------------------------------------------------
    900  continue

end subroutine nitjv


subroutine nitgm(n, xcur, fcur, fcnrm, step, eta, f, jacv, rpar, &
        ipar, nfe, njve, nrpre, nli, kdmaxp1, vv, rr, svbig, svsml, w, rwork, &
        rsnrm, itrmks)

   integer, intent(in) :: n, kdmaxp1
   integer, intent(in out) :: ipar(*), nfe, njve, nrpre, nli, itrmks
   double precision, intent(in) :: eta, fcnrm
   double precision, intent(in out) :: xcur(n), fcur(n), step(n), rpar(*), &
      vv(n,kdmaxp1), rr(kdmax_glob,kdmax_glob), svbig(kdmax_glob), svsml(kdmax_glob), &
      w(kdmax_glob), rwork(n), rsnrm
   procedure(funcInterface) :: f
   procedure(jacvInterface) :: jacv

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
     step(1:n) = 0
   !   do 20 i = 1, n
   !      step(i) = 0.d0
   !  20   continue
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
     vv(1:n, 1) = fcur!call dcopy(n, fcur, 1, vv(1,1), 1)
     temp = -1.d0/fcnrm
     vv(1:n, 1) = temp * vv(1:n, 1)!call dscal(n, temp, vv(1,1), 1)
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
     if (irpre_glob .eq. 0) then 
        itask = 0
     else
        itask = 1
     endif
     call nitjv(n, xcur, fcur, f, jacv, rpar, ipar, itask, nfe, njve, nrpre, vv(1,kd), vv(1,kdp1), rwork, rwork, itrmjv)
     if (itrmjv .gt. 0) then 
        if (itrmjv .eq. 1) itrmks = 1
        if (itrmjv .eq. 2) itrmks = 2
        go to 900
     endif
    ! c ------------------------------------------------------------------------
    ! c Do modified Gram-Schmidt. 
    ! c ------------------------------------------------------------------------
     do 210 i = 2, kd
        rr(i-1,kd) = SUM(vv(1:n,i) * vv(1:n, kdp1))!ddot(n, vv(1,i), 1, vv(1,kdp1), 1)
        vv(1:n, kdp1) = vv(1:n, kdp1) -rr(i-1, kd) * vv(1:n, i)!call daxpy(n, -rr(i-1,kd), vv(1,i), 1, vv(1,kdp1), 1)
    210  continue
     rr(kd,kd) = SQRT(SUM(vv(1:n, kdp1)**2)) !dnrm2(n, vv(1,kdp1), 1)
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
        svbig(1:kd-1) = svbig(1:kd-1) * sn !call dscal(kd-1, sn, svbig, 1)
        svbig(kd) = cs
        ijob = 2
        call dlaic1(ijob, kd-1, svsml, small, rr(1,kd), rr(kd,kd),sestpr, sn, cs)
        small = sestpr
        svsml(1:kd-1) = svsml(1:kd-1) * sn !call dscal(kd-1, sn, svsml, 1)
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
           vv(1:n, 1) = vv(1:n, 1) + w(kd) * vv(1:n, kdp1) !call daxpy(n, w(kd), vv(1,kdp1), 1, vv(1,1), 1)
           go to 300
        endif
     endif
    ! c ------------------------------------------------------------------------
    ! c Normalize vv(.,kdp1). 
    ! c ------------------------------------------------------------------------
     temp = 1.d0/rr(kd,kd) 
     vv(1:n, kdp1) = vv(1:n, kdp1) * temp !call dscal(n, temp, vv(1,kdp1), 1)
    ! c ------------------------------------------------------------------------
    ! c Update w and the residual norm by rsnrm <- rsnrm*dsin(dacos(w(kd)/rsnrm). 
    ! c ------------------------------------------------------------------------
     w(kd) = SUM(vv(1:n, 1) * vv(1:n, kdp1))!ddot(n, vv(1,1), 1, vv(1,kdp1), 1)
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
     if ( (rsnrm0*rsnrm .le. abstol) .or. (kd .eq. kdmax_glob) .or.  (igm .ge. iksmax_glob) )  go to 300
    ! c ------------------------------------------------------------------------
    ! c If not terminating the inner loop, update the residual vector 
    ! c and go to the top of the inner loop. 
    ! c ------------------------------------------------------------------------
     vv(1:n, 1) = vv(1:n, 1) - w(kd) * vv(1:n, kdp1)!call daxpy(n, -w(kd), vv(1,kdp1), 1, vv(1,1), 1)
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
     svbig(1:kd) = w(1:kd) !call dcopy(kd, w, 1, svbig, 1)
    ! c ------------------------------------------------------------------------
    ! c Overwrite w with the solution of the upper triangular system.
    ! c ------------------------------------------------------------------------
     do 310 i = kd, 1, -1
        w(i) = w(i)/rr(i,i)
        if (i .gt. 1) w(1:i-1) = w(1:i-1) - w(i) * rr(1:i-1, i) !call daxpy(i-1, -w(i), rr(1,i), 1, w, 1)
    310  continue
    ! c ------------------------------------------------------------------------
    ! c Now form the linear combination to accumulate the correction in 
    ! c the work vector.
    ! c ------------------------------------------------------------------------
     rwork(1:n) = vv(1:n, 1) !call dcopy(n, vv(1,1), 1, rwork, 1)
     rwork(1:n) = rwork(1:n) * w(1) !call dscal(n, w(1), rwork, 1)
     if (kd .gt. 1) then 
        w(2:kd) =  w(2:kd) + w(1) * svbig(1:kd-1)!call daxpy(kd-1, w(1), svbig, 1, w(2), 1)
        do 320 i = 2, kd
           rwork(1:n) = rwork(1:n) + w(i) * vv(1:n, i)!call daxpy(n, w(i), vv(1,i), 1, rwork, 1)
    320     continue
     endif
    ! c ------------------------------------------------------------------------
    ! c If iresup .eq. 0, then update the residual vector by linear 
    ! c combination. This frees vv(.,kdp1) for use as a work array. 
    ! c ------------------------------------------------------------------------
     if (iresup_glob .eq. 0) then 
        vv(1:n, 1) = vv(1:n, 1) -svbig(kd) * vv(1:n, kdp1) !call daxpy(n, -svbig(kd), vv(1,kdp1), 1, vv(1,1), 1)
     endif
    ! c ------------------------------------------------------------------------
    ! c If right preconditioning is used, overwrite 
    ! c correction <-- P(inverse)*correction, using vv(.,kdp1) as a work array. 
    ! c Note: vv(.,kdp1) can be used for both work arrays in this call because 
    ! c the second is not referenced within nitjv. 
    ! c ------------------------------------------------------------------------
     if (irpre_glob .gt. 0) then 
        itask = 2
        call nitjv(n, xcur, fcur,f, jacv, rpar, ipar, itask, nfe, njve, nrpre, rwork, rwork, vv(1,kdp1), vv(1,kdp1), itrmjv)
        if (itrmjv .gt. 0) then 
           itrmks = 2
           go to 900
        endif
     endif
    ! c ------------------------------------------------------------------------
    ! c Update the step. This frees rwork for use as a work array.
    ! c ------------------------------------------------------------------------
     step = step + rsnrm0 * rwork(1:n)!call daxpy(n, rsnrm0, rwork, 1, step, 1)
    ! c ------------------------------------------------------------------------
    ! c If iresup .eq. 1, then update the residual vector by direct evaluation, 
    ! c using rwork and vv(.,kdp1) as work arrays. Note: Two distinct work  
    ! c arrays are needed in this call because both are referenced within nitjv 
    ! c if the J*step product is evaluated with a finite-difference of order 
    ! c two or higher. If finite-differences are used (ijacv= 0), then ijacv 
    ! c is temporarily set to -1 to signal to nitjv that the order of the 
    ! c finite-difference formula is to be determined by ifdord. 
    ! c ------------------------------------------------------------------------
     if (iresup_glob .eq. 1) then 
        itask = 0
        if (ijacv_glob .eq. 0) ijacv_glob = -1 
        call nitjv(n, xcur, fcur, f, jacv, rpar, ipar, itask, nfe, njve, nrpre, step, vv(1,1), rwork, vv(1,kdp1), itrmjv)
        if (ijacv_glob .eq. -1) ijacv_glob = 0
        if (itrmjv .gt. 0) then 
           itrmks = 1
           go to 900
        endif
        vv(1:n, 1) = -fcur - vv(1:n, 1)
   !      do 330 i = 1, n
   !         vv(i,1) = -fcur(i) - vv(i,1)
   !  330     continue
     endif
    ! c ------------------------------------------------------------------------
    ! c Test for termination.  
    ! c ------------------------------------------------------------------------
     if (rsnrm0*rsnrm .le. abstol) then 
        itrmks = 0
        go to 900
     endif
     if (igm .ge. iksmax_glob) then 
        itrmks = 3
        go to 900
     endif
     if (big .ge. small*cndmax) then 
        itrmks = 5
        go to 900
     endif
     temp = dfloat(kd)*dlog(abstol/(rsnrm0*rsnrm))/dlog(rsnrm/(1.d0 + 10.d0*epsmach))
     if (temp .ge. 1000.d0*dfloat(iksmax_glob - igm)) then 
        itrmks = 4
        go to 900
     endif
    ! c ------------------------------------------------------------------------
    ! c If not terminating, then normalize the initial residual, etc., and 
    ! c return to the top of the outer loop. 
    ! c ------------------------------------------------------------------------
     if (iresup_glob .eq. 0) then 
        rsnrm0 = rsnrm0*rsnrm
        temp = 1.d0/rsnrm
     else
        rsnrm0 = SQRT(SUM(vv(1:n,1)**2)) !dnrm2(n, vv(1,1), 1)
        temp = 1.d0/rsnrm0
     endif
     vv(1:n, 1) = vv(1:n, 1) * temp !call dscal(n, temp, vv(1,1), 1)
     go to 100
    ! c ------------------------------------------------------------------------
    ! c All returns made here.
    ! c ------------------------------------------------------------------------
    900  continue
     if (itrmks .ne. 1 .and. itrmks .ne. 2) then 
        if (iresup_glob .eq. 0) then 
           vv(1:n, 1) = vv(1:n, 1) * rsnrm0 !call dscal(n, rsnrm0, vv(1,1), 1)
           rsnrm = rsnrm0*rsnrm
        else
           rsnrm = SQRT(SUM(vv(1:n,1)**2)) !dnrm2(n, vv(1,1), 1)
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


subroutine nittfq (n, xcur, fcur, fcnrm, step, eta, f, jacv, rpar, ipar, &
        nfe, njve, nrpre, nli, r,rcgs, rtil, d, p, q, u, v, y, rwork1, rwork2, &
        rsnrm, itrmks )

   integer :: itrmks, n, nfe, njve,nrpre, nli

   integer :: ipar(*)

   double precision :: eta, fcnrm, rsnrm

   double precision :: d(n), fcur(n), p(n), q(n), rcgs(n), r(n),rpar(*), rtil(n), rwork1(n), rwork2(n), step(n),u(n), v(n), xcur(n), y(n)
   procedure(funcInterface) :: f
   procedure(jacvInterface) :: jacv

    ! c ------------------------------------------------------------------------
    ! c
    ! c This is nittfq v0.1, the TFQMR routine for determining (trial) inexact 
    ! c Newton steps. The original reference is R. W. Freund, "A Transpose-Free
    ! c Quasi-Minimal Residual Algorithm for Non-Hermitian Linear Systems", 
    ! c SIAM J. Sci. Comput., 14 (1993), pp. 470-482.  The implementation here
    ! c is based on the right preconditioned algorithm that appears in 
    ! c J. N. Shadid and R. S. Tuminaro, "A Comparison of Preconditioned
    ! c Nonsymmetric Krylov Methods on a Large-Scale MIMD Machine", SIAM J.
    ! c Sci. Comput., 15 (1994), pp. 440-459.
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
    ! c  f      = name of user-supplied subroutine for evaluating the function 
    ! c           the zero of which is sought; this routine has the form
    ! c
    ! c                 subroutine f(n, xcur, fcur, rpar, ipar, itrmf)
    ! c
    ! c           where xcur is the array containing the current x value, fcur 
    ! c           is f(xcur) on output, rpar and ipar are, respectively, real 
    ! c           and integer parameter/work arrays for use by the subroutine,
    ! c           and itrmf is an integer termination flag.  The meaning of
    ! c           itrmf is as follows:
    ! c             0 => normal termination; desired function value calculated.
    ! c             1 => failure to produce f(xcur).
    ! c 
    ! c  jacv   = name of user-supplied subroutine for evaluating J*v or 
    ! c           P(inverse)*v, where J is the Jacobian of f and P is a 
    ! c           right preconditioning operator. If neither analytic J*v 
    ! c           evaluations nor right preconditioning is used, this can 
    ! c           be a dummy subroutine; if right preconditioning is used but 
    ! c           not analytic J*v evaluations, this need only evaluate 
    ! c           P(inverse)*v. The form is 
    ! c
    ! c           subroutine jacv(n, xcur, fcur, ijob, v, z, rpar, ipar, itrmjv)
    ! c
    ! c           where xcur and fcur are vectors of length n containing the 
    ! c           current x and f values, ijob is an integer flag indicating 
    ! c           which product is desired, v is a vector of length n to be 
    ! c           multiplied, z is a vector of length n containing the desired 
    ! c           product on output, rpar and ipar are, respectively, real 
    ! c           and integer parameter/work arrays for use by the subroutine, 
    ! c           and itrmjv is an integer termination 
    ! c           flag. The meaning of ijob is as follows: 
    ! c             0 => z = J*v
    ! c             1 => z = P(inverse)*v 
    ! c           The meaning of itrmjv is as follows:
    ! c             0 => normal termination; desired product evaluated. 
    ! c             1 => failure to produce J*v.
    ! c             2 => failure to produce P(inverse)*v. 
    ! c           This subroutine is called only from nitjv, and is always 
    ! c           called with v .ne. z. 
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
    ! c  iksmax  = maximum allowable number of TFQMR iterations. 
    ! c
    ! c  ifdord  = order of the finite-difference formula used in BiCGSTAB 
    ! c            when J*v products are evaluated using finite-differences. 
    ! c            When ijacv = 0 on input to nitsol, ifdord is set to 1, 2, or 
    ! c            4 in nitsol; otherwise, it is irrelevant. When ijacv = 0 on 
    ! c            input to this subroutine, ifdord determines the order of the 
    ! c            finite-difference formula used at each BiCGSTAB iteration 
    ! c            (default 1). In this case, ijacv is set to -1 below to 
    ! c            signal to nitjv that the order of the finite-difference 
    ! c            formula is to be determined by ifdord. The original value 
    ! c            ijacv = 0 is restored on return. 
    ! c            
    ! c  nfe     = number of function evaluations.
    ! c
    ! c  njve    = number of J*v evaluations. 
    ! c
    ! c  nrpre   = number of P(inverse)*v evaluations.
    ! c
    ! c  nli     = number of linear iterations.
    ! c
    ! c  r       = residual vector (for the QMR process)
    ! c
    ! c  rcgs    = residual vector (of the underlying CGS process)
    ! c
    ! c  rtil    = 'shadow' residual vector used in bi-orthogonalization
    ! c
    ! c  d       = vector used in TFQMR
    ! c
    ! c  p       = vector used in TFQMR
    ! c
    ! c  q       = vector used in TFQMR
    ! c
    ! c  u       = vector used in TFQMR
    ! c
    ! c  v       = vector used in TFQMR
    ! c
    ! c  y       = vector used in TFQMR
    ! c
    ! c  rwork1  = work vector, passed on to nitjv
    ! c
    ! c  rwork2  = work vector, passed on to nitjv
    ! c
    ! c  rsnrm   = TFQMR residual norm on return. 
    ! c
    ! c  dinpr   = inner-product routine, either user-supplied or blas ddot. 
    ! c
    ! c  dnorm   = norm routine, either user-supplied or blas dnrm2. 
    ! c
    ! c  itrmks  = termination flag; values have the following meanings: 
    ! c              0 => normal termination: acceptable step found. 
    ! c              1 => J*v failure in nitjv. 
    ! c              2 => P(inverse)*v failure in nitjv. 
    ! c              3 => acceptable step not found in iksmax TFQMR iterations. 
    ! c              4 => TFQMR breakdown.
    ! c              5 => floating point error (the underlying CGS iteration
    ! c                   has probably blown up)
    ! c
    ! c             Note: On return, nitsol terminates if itrmks is 1 or 2. 
    ! c             If itrmks is 3 or 4, nitsol may terminate or continue. 
    ! c             In this event, the step returned is a meaningful inexact 
    ! c             Newton step only if the residual norm has been reduced. 
    ! c             A decision on termination/continuation is made in nitdrv 
    ! c             according to whether there is sufficient residual norm 
    ! c             reduction, even though the desired inexact Newton condition 
    ! c             may not hold.  
    ! c
    ! c -------------------------------------------------------------------------

    ! c  Subroutines required by this and all called routines:

    ! c    user supplied:  nitjv

    ! c    nitsol routines: none

    ! c    BLAS routines- dcopy, daxpy, dscal, dswap

    ! c    LAPACK routines - dlamch

    ! c    user supplied or BLAS:  dinpr, dnorm

    ! c    explanation: In nitsol, dinpr and dnorm are set to either the BLAS 
    ! c    ddot and dnrm2 routines or the user-supplied routines. 

    ! c This subroutine called by: nitdrv

    ! c Subroutines called by this subroutine: daxpy, dcopy, dscal, dswap, dinpr,
    ! c    dlamch, dnorm, nitjv, dlamch

    ! c Common block: 

    !      include 'nitprint.h'

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


     integer :: i, itask, itrmjv, itfq, k

     double precision :: alpha, abstol, beta, c, cgsnorm, qmreta, omega, rho, rho_old, sigma, t, tau, theta

     character*2 :: ab(0:1)
     data ab /'.a','.b'/


    ! c If finite-differences are used to evaluate J*v products (ijacv= 0), then 
    ! c ijacv is set to -1 within this subroutine to signal to nitjv that the 
    ! c order of the finite-difference formula is to be determined by ifdord. 
    ! c The original value ijacv= 0 is restored on return. 

     if (ijacv_glob .eq. 0) ijacv_glob = -1 

    ! c Set the stopping tolerance, initialize the step, etc. 
     rsnrm = fcnrm
     abstol = eta*rsnrm
     step(1:n) = 0.0d0
   !   do 10 i = 1, n
   !      step(i) = zero
   !  10   continue
     itfq = 0

    ! c For printing:

     if ( iplvl .ge. 3 ) then 
        write(ipunit,*) 
        write(ipunit,800) eta 
     endif
     if ( iplvl .ge. 4 ) then 
        write(ipunit,810) 
        write(ipunit,*) 
        write(ipunit,820) itfq, rsnrm 
     endif

    ! c  Initialize residual and work vectors.

     rcgs = fcur !call dcopy( n, fcur, 1, rcgs, 1 )
     rcgs = -rcgs !call dscal( n, -one, rcgs, 1 )

    ! c  Choice here is rtil = r.

     rtil = rcgs!call dcopy( n, rcgs, 1, rtil, 1 )

     if ( irpre_glob .eq. 0 ) then
        p = rcgs !call dcopy( n, rcgs, 1, p, 1 )
     else
        itask = 2
        call nitjv( n, xcur, fcur, f, jacv, rpar, ipar,itask, nfe, njve, nrpre, rcgs, p, rwork1, rwork2, itrmjv )    
        if ( itrmjv .ne. 0 ) then
           itrmks = 2
           goto 900
        endif
     endif
     u = p!call dcopy( n, p, 1, u, 1 )
     rho = SUM(rcgs * rtil)!ddot( n, rcgs, 1, rtil, 1 )
     d(1:n) = 0.0d0
   !   do 20 i = 1, n
   !      d(i) = zero
   !  20   continue

     itask = 0
     call nitjv( n, xcur, fcur, f, jacv, rpar, ipar,itask, nfe, njve, nrpre, p, v, rwork1, rwork2, itrmjv )
     if ( itrmjv .ne. 0 ) then
        itrmks = 1
        goto 900
     end if

     alpha = zero
     omega = fcnrm
     tau = omega
     theta = zero
     qmreta = zero

    ! c  Start iterations.

    100  continue
     itfq = itfq + 1
     nli = nli + 1
     sigma = SUM(rtil * v)!ddot( n, rtil, 1, v, 1 )

    ! c  If sigma = 0 we have a serious breakdown.  We check this condition
    ! c  by trying to detect whether division by sigma causes an overflow.

     if ( abs(sigma) .lt. sfmin*abs(rho) ) then
        itrmks = 4
        goto 900
     else
        alpha = rho/sigma
     endif
    
    ! c  Need Pv for calculation of q.  First store result in q, then
    ! c  swap some vectors to cast calculation of q as a SAXPY.

     if ( irpre_glob .eq. 0 ) then
        q = v!call dcopy( n, v, 1, q, 1 )
     else
        itask = 2
        call nitjv( n, xcur, fcur, f, jacv, rpar, ipar,itask, nfe, njve, nrpre, v, q, rwork1, rwork2, itrmjv )    
        if ( itrmjv .ne. 0 ) then
           itrmks = 2
           goto 900
        endif
     endif
     y = u!call dcopy( n, u, 1, y, 1 )
     !call dswap( n, q, 1, u, 1 )
   !   print *, u
   !   print *, q
     do i = 1, n
         cgsnorm = u(i)
         u(i) = q(i)
         q(i) = cgsnorm
     end do
     q = q - alpha * u!call daxpy( n, -alpha, u, 1, q, 1 )
     u = y !call dcopy( n, y, 1, u, 1 )

    ! c  Update residual.
     y = u + q
   !   do 30 i = 1, n
   !      y(i) = u(i) + q(i)
   !  30   continue
     itask = 0
     call nitjv( n, xcur, fcur, f, jacv, rpar, ipar,itask, nfe, njve, nrpre, y, v, rwork1, rwork2, itrmjv )
     if ( itrmjv .ne. 0 ) then
        itrmks = 1
        goto 900
     end if
     rcgs = rcgs - alpha * v !call daxpy( n, -alpha, v, 1, rcgs, 1 )
     cgsnorm = SQRT(SUM(rcgs**2))!dnrm2( n, rcgs, 1 )

    ! c  Check for cgsnorm = NaN.

     if ( cgsnorm .ne. cgsnorm ) then
        itrmks = 5
        goto 900
     endif

    ! c  QMR section.

     do 60 k = 0, 1

    ! c  Use weighting strategy from (5.11) of Freund reference.

        t = qmreta*theta**2
        t = t/alpha

        if ( k .eq. 0 ) then
            d = u + t*d
   !         do 40 i = 1, n
   !            d(i) = u(i) + t*d(i)
   !  40         continue
           omega = sqrt(omega*cgsnorm)
        else if ( k .eq. 1 ) then
            d = q + t * d
   !         do 50 i = 1, n
   !            d(i) = q(i) + t*d(i)
   !  50         continue
           omega = cgsnorm
        endif

        theta = omega/tau
        c = one/sqrt(one + theta**2)
        tau = tau*theta*c
        qmreta = alpha*c**2

    ! c For printing:

        if ( iplvl .ge. 4 .and. tau .gt. abstol ) then 
           write(ipunit,830) itfq, ab(k), tau, '     (estimated)'
        endif

      step = step + qmreta * d!call daxpy( n, qmreta, d, 1, step, 1 )

    ! c  Convergence check.  Do a cheap test to save on Jacobi-vector products.
    ! c  In case residual history is requested by iplvl, we must calculate 
    ! c  the QMR residual from scratch.  Note termination is always determined
    ! c  by the smoothed residual, calculated from scratch.

        if ( tau .le. abstol ) then

           itask = 0
           call nitjv( n, xcur, fcur, f, jacv, rpar, ipar, itask, nfe, njve, nrpre, step, r, rwork1,rwork2, itrmjv )

    ! c  This calculation of the QMR residual is off by a factor
    ! c  of -1, but we dont care until we return from this routine.

           r = r + fcur !call daxpy( n, one, fcur, 1, r, 1 )
           rsnrm = SQRT(SUM(r**2))!dnrm2( n, r, 1 )
    ! c For printing:

           if ( iplvl .ge. 4 ) then 
              write(ipunit,830) itfq, ab(k), rsnrm, '  (from scratch)'
           endif

    ! c  Check for rsnrm = NaN.

           if ( rsnrm .ne. rsnrm ) then
              itrmks = 5
              goto 900
           endif

    ! c  If rsnrm is small enough, exit.

           if ( rsnrm .lt. abstol ) then
              itrmks = 0
              goto 900
           endif
        endif

    60   continue

     rho_old = rho
     rho = SUM(rtil * rcgs) !ddot( n, rtil, 1, rcgs, 1 )

    ! c  If rho_old = 0 we have a serious breakdown.  We check this condition
    ! c  by trying to detect whether division by rho_old causes an overflow.

     if ( abs(rho_old) .lt. sfmin*abs(rho) ) then
        itrmks = 4
        goto 900
     else
        beta = rho/rho_old
     endif

     if ( irpre_glob .eq. 0 ) then
        v = rcgs !call dcopy( n, rcgs, 1, v, 1 )
     else
        itask = 2
        call nitjv( n, xcur, fcur, f, jacv, rpar, ipar,itask, nfe, njve, nrpre, rcgs, v, rwork1, rwork2, itrmjv )    
        if ( itrmjv .ne. 0 ) then
           itrmks = 2
           goto 900
        endif
     endif
     v = v + beta * q !call daxpy( n, beta, q, 1, v, 1 )
     u = v!call dcopy( n, v, 1, u, 1 )
     q = q + beta * p !call daxpy( n, beta, p, 1, q, 1 )
     v = v + beta * q !call daxpy( n, beta, q, 1, v, 1 )
     p = v !call dcopy( n, v, 1, p, 1 )

     itask = 0
     call nitjv( n, xcur, fcur, f, jacv, rpar, ipar,itask, nfe, njve, nrpre, p, v, rwork1, rwork2, itrmjv )
     if ( itrmjv .ne. 0 ) then
        itrmks = 1
        goto 900
     end if
     if ( itfq .ge. iksmax_glob ) then
        itrmks = 3
        goto 900
     end if

    ! c  Do again

     goto 100

    ! c  All returns made here.

    900 continue

    ! c  If residual hasnt been updated, force
    ! c  computation of residual from scratch.

     if ( rsnrm .eq. fcnrm ) then

        itask = 0
        call nitjv( n, xcur, fcur, f, jacv, rpar, ipar, itask, nfe, njve, nrpre, step, r, rwork1,rwork2, itrmjv )

        r = r + fcur !call daxpy( n, one, fcur, 1, r, 1 )
        rsnrm = SQRT(SUM(r**2)) !dnrm2( n, r, 1 )
        if ( rsnrm .le. abstol ) itrmks = 0

     end if

    ! c  Correct residual before returning.

     r = -r !call dscal( n, -one, r, 1 )

    ! c If ijacv = -1, then restore it to the original value ijacv = 0. 

     if (ijacv_glob .eq. -1) ijacv_glob = 0 

    ! c For printing:

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
    810  format('nittfq:  TFQMR iteration no. (parts a and b),',' linear residual norm')
    820  format(5x,i4,5x,1pd10.3)
    830  format(5x,i4,a2,3x,1pd10.3,a16)
    840  format('nittfq:  itrmks =', i2, '   final lin. res. norm =', 1pd10.3)
    850  format('nittfq: itrmks:', i4) 


end subroutine nittfq

subroutine nitstb (n, xcur, fcur, fcnrm, step, eta, f, jacv, rpar, ipar, &
    nfe, njve, nrpre, nli, r, rtil, p, phat, v, &
    t, rwork1, rwork2, rsnrm, itrmks)


   integer :: n, ipar(*), nfe, njve, nrpre, nli, itrmks
   double precision :: xcur(n), fcur(n), fcnrm, step(n), eta, rpar(*), &
   r(n), rtil(n), p(n), phat(n), v(n), t(n), rwork1(n), rwork2(n), rsnrm
   procedure(funcInterface) :: f
   procedure(jacvInterface) :: jacv

    ! c ------------------------------------------------------------------------
    ! c
    ! c This is nitstb v0.3, the BiCGSTAB routine for determining (trial) inexact 
    ! c Newton steps. The original reference is H. van der Vorst, "Bi-CGSTAB: 
    ! c A fast and smoothly converging variant of Bi-CG for the soluton of 
    ! c nonsymmetric linear systems," SIAM J. Sci. Statist. Comput., 13 (1992), 
    ! c pp. 631--644. 
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
    ! c  f      = name of user-supplied subroutine for evaluating the function 
    ! c           the zero of which is sought; this routine has the form 
    ! c
    ! c                 subroutine f(n, xcur, fcur, rpar, ipar, itrmf)
    ! c
    ! c           where xcur is the array containing the current x value, fcur 
    ! c           is f(xcur) on output, rpar and ipar are, respectively, real 
    ! c           and integer parameter/work arrays for use by the subroutine,
    ! c           and itrmf is an integer termination flag.  The meaning of
    ! c           itrmf is as follows:
    ! c             0 => normal termination; desired function value calculated.
    ! c             1 => failure to produce f(xcur).
    ! c 
    ! c  jacv   = name of user-supplied subroutine for evaluating J*v or 
    ! c           P(inverse)*v, where J is the Jacobian of f and P is a 
    ! c           right preconditioning operator. If neither analytic J*v 
    ! c           evaluations nor right preconditioning is used, this can 
    ! c           be a dummy subroutine; if right preconditioning is used but 
    ! c           not analytic J*v evaluations, this need only evaluate 
    ! c           P(inverse)*v. The form is 
    ! c
    ! c           subroutine jacv(n, xcur, fcur, ijob, v, z, rpar, ipar, itrmjv)
    ! c
    ! c           where xcur and fcur are vectors of length n containing the 
    ! c           current x and f values, ijob is an integer flag indicating 
    ! c           which product is desired, v is a vector of length n to be 
    ! c           multiplied, z is a vector of length n containing the desired 
    ! c           product on output, rpar and ipar are, respectively, real 
    ! c           and integer parameter/work arrays for use by the subroutine, 
    ! c           and itrmjv is an integer termination 
    ! c           flag. The meaning of ijob is as follows: 
    ! c             0 => z = J*v
    ! c             1 => z = P(inverse)*v 
    ! c           The meaning of itrmjv is as follows:
    ! c             0 => normal termination; desired product evaluated. 
    ! c             1 => failure to produce J*v.
    ! c             2 => failure to produce P(inverse)*v. 
    ! c           This subroutine is called only from nitjv, and is always 
    ! c           called with v .ne. z. 
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
    ! c  iksmax  = maximum allowable number of BiCGSTAB iterations. 
    ! c
    ! c  ifdord  = order of the finite-difference formula used in BiCGSTAB 
    ! c            when J*v products are evaluated using finite-differences. 
    ! c            When ijacv = 0 on input to nitsol, ifdord is set to 1, 2, or 
    ! c            4 in nitsol; otherwise, it is irrelevant. When ijacv = 0 on 
    ! c            input to this subroutine, ifdord determines the order of the 
    ! c            finite-difference formula used at each BiCGSTAB iteration 
    ! c            (default 1). In this case, ijacv is set to -1 below to 
    ! c            signal to nitjv that the order of the finite-difference 
    ! c            formula is to be determined by ifdord. The original value 
    ! c            ijacv = 0 is restored on return. 
    ! c            
    ! c  nfe     = number of function evaluations.
    ! c
    ! c  njve    = number of J*v evaluations. 
    ! c
    ! c  nrpre   = number of P(inverse)*v evaluations.
    ! c
    ! c  nli     = number of linear iterations.
    ! c
    ! c  r       = residual vector 
    ! c
    ! c  rtil    = "r-tilde" vector used in BiCGSTAB
    ! c
    ! c  p       = vector used in BiCGSTAB
    ! c
    ! c  phat    = vector used in BiCGSTAB
    ! c
    ! c  v       = vector used in BiCGSTAB
    ! c
    ! c  t       = vector used in BiCGSTAB
    ! c
    ! c  rwork1  = work vector, passed on to nitjv
    ! c
    ! c  rwork2  = work vector, passed on to nitjv
    ! c
    ! c  rsnrm   = BiCGSTAB residual norm on return. 
    ! c
    ! c  dinpr   = inner-product routine, either user-supplied or blas ddot. 
    ! c
    ! c  dnorm   = norm routine, either user-supplied or blas dnrm2. 
    ! c
    ! c  itrmks  = termination flag; values have the following meanings: 
    ! c              0 => normal termination: acceptable step found. 
    ! c              1 => J*v failure in nitjv. 
    ! c              2 => P(inverse)*v failure in nitjv. 
    ! c              3 => acceptable step not found in iksmax BiCGSTAB iterations. 
    ! c              4 => BiCGSTAB breakdown. 
    ! c
    ! c             Note: On return, nitsol terminates if itrmks is 1 or 2. 
    ! c             If itrmks is 3 or 4, nitsol may terminate or continue. 
    ! c             In this event, the step returned is a meaningful inexact 
    ! c             Newton step only if the residual norm has been reduced. 
    ! c             A decision on termination/continuation is made in nitdrv 
    ! c             according to whether there is sufficient residual norm 
    ! c             reduction, even though the desired inexact Newton condition 
    ! c             may not hold.  
    ! c
    ! c -------------------------------------------------------------------------
    ! c
    ! c Subroutines required by this and all called routines: 
    ! c
    ! c    user supplied: f, jacv 
    ! c
    ! c    nitsol routines: nitjv, nitfd
    ! c
    ! c    blas routines: daxpy, dcopy, dscal
    ! c
    ! c    lapack routine:  dlamch
    ! c
    ! c    user supplied or blas: dinpr, dnorm 
    ! c
    ! c    explanation: In nitsol, dinpr and dnorm are set to either the blas 
    ! c    ddot and dnrm2 routines or the user-supplied usrnpr and usrnrm 
    ! c    routines. 
    ! c
    ! c This subroutine called by: nitdrv
    ! c
    ! c Subroutines called by this subroutine: daxpy, dcopy, dscal, dinpr, dlamch,
    ! c    dnorm, nitjv
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
     double precision :: abstol, alpha, beta, omega, rho, rhomns, tau, temp
     integer :: istb, itask, itrmjv

    ! c ------------------------------------------------------------------------
    ! c If finite-differences are used to evaluate J*v products (ijacv= 0), then 
    ! c ijacv is set to -1 within this subroutine to signal to nitjv that the 
    ! c order of the finite-difference formula is to be determined by ifdord. 
    ! c The original value ijacv= 0 is restored on return. 
    ! c ------------------------------------------------------------------------
     if (ijacv_glob .eq. 0) ijacv_glob = -1 
    ! c ------------------------------------------------------------------------
    ! c Set the stopping tolerance, initialize the step, etc. 
    ! c ------------------------------------------------------------------------
     rsnrm = fcnrm
     abstol = eta*rsnrm
     step(1:n) = 0.0d0
   !   do 10 i = 1, n
   !      step(i) = 0.d0
   !  10   continue
     istb = 0
    ! c ------------------------------------------------------------------------ 
    ! c For printing:
     if (iplvl .ge. 3) then 
        write(ipunit,*) 
        write(ipunit,800) eta 
     endif
    800  format('nitstb:  eta =', 1pd10.3)
     if (iplvl .ge. 4) then 
        write(ipunit,810) 
        write(ipunit,*) 
        write(ipunit,820) istb, rsnrm 
     endif
    810  format('nitstb:  BiCGSTAB iteration no. (parts a and b)',' linear residual norm, ')
    820  format(5x,i4,5x,1pd10.3)
    ! c ------------------------------------------------------------------------
    ! c
    ! c ------------------------------------------------------------------------
    ! c Set up r and rtil. 
    ! c ------------------------------------------------------------------------
     r = fcur!call dcopy(n,fcur,1,r,1)
     temp = -1.d0
     r = temp * r!call dscal(n,temp,r,1)
     rtil = r!call dcopy(n,r,1,rtil,1)
    ! c ------------------------------------------------------------------------
    ! c Top of the iteration loop. 
    ! c ------------------------------------------------------------------------
    100  continue
     istb = istb + 1
     nli = nli + 1
    ! c ------------------------------------------------------------------------
    ! c Perform the first "half-iteration". 
    ! c ------------------------------------------------------------------------
     rho = SUM(rtil * r)!ddot(n,rtil,1,r,1)
     if (istb .eq. 1) then 
        p = r!call dcopy(n,r,1,p,1)
     else
        if ( abs(rhomns) .lt. sfmin*abs(rho) ) then
           itrmks = 4
           goto 900
        else
           beta = (rho/rhomns)*(alpha/omega)
           p = p - omega*v!call daxpy(n,-omega,v,1,p,1)
           p = beta * p!call dscal(n,beta,p,1)
           p = p + r!call daxpy(n,1.d0,r,1,p,1)
        endif
     endif
     if (irpre_glob .eq. 0) then 
        phat = p!call dcopy(n,p,1,phat,1)
     else
        itask = 2
        call nitjv(n, xcur, fcur, f, jacv,rpar, ipar, itask, nfe, njve, nrpre, p, phat, rwork1, rwork2, itrmjv)
        if (itrmjv .gt. 0) then 
           itrmks = 2
           go to 900
        endif
     endif
     itask = 0
     call nitjv(n, xcur, fcur, f, jacv,rpar, ipar, itask, nfe, njve, nrpre, phat, v, rwork1, rwork2, itrmjv)
     if (itrmjv .gt. 0) then 
        itrmks = 1
        go to 900
     endif
     tau = SUM(rtil * v)!ddot(n,rtil,1,v,1)
     if ( abs(tau) .lt. sfmin*abs(rho) ) then
        itrmks = 4
        goto 900
     else
        alpha = rho/tau
     endif
     r = r - alpha*v!call daxpy(n,-alpha,v,1,r,1)
     step = step + alpha * phat!call daxpy(n,alpha,phat,1,step,1)
     rsnrm = SQRT(SUM(r**2))!dnrm2(n,r,1)
    ! c ------------------------------------------------------------------------ 
    ! c For printing:
     if (iplvl .ge. 4) then 
        write(ipunit,830) istb, rsnrm 
    830     format(5x,i4,'.a',3x,1pd10.3)
     endif
    ! c ------------------------------------------------------------------------
    ! c
    ! c ------------------------------------------------------------------------
    ! c Test for termination. 
    ! c ------------------------------------------------------------------------
     if (rsnrm .le. abstol) then 
        itrmks = 0
        go to 900
     endif
    ! c ------------------------------------------------------------------------
    ! c Perform the second "half-iteration". 
    ! c ------------------------------------------------------------------------
     if (irpre_glob .eq. 0) then 
        phat = r!call dcopy(n,r,1,phat,1)
     else
        itask = 2
        call nitjv(n, xcur, fcur, f, jacv,rpar, ipar, itask, nfe, njve, nrpre, r, phat, rwork1, rwork2, itrmjv)
        if (itrmjv .gt. 0) then 
           itrmks = 2
           go to 900
        endif
     endif
     itask = 0
     call nitjv(n, xcur, fcur, f, jacv,rpar, ipar, itask, nfe, njve, nrpre, phat, t, rwork1, rwork2, itrmjv)
     if (itrmjv .gt. 0) then 
        itrmks = 1
        go to 900
     endif
     tau = SQRT(SUM(t**2))!dnrm2(n,t,1)
     tau = tau*tau
     temp = SUM(t*r) !ddot(n,t,1,r,1)
     if ( tau .le. sfmin*abs(temp) ) then
        itrmks = 4
        goto 900
     else
        omega = temp/tau
     endif
     if ( abs(omega) .lt. sfmin*abs(alpha) ) then 
        itrmks = 4
        go to 900
     endif
     r = r -omega*t!call daxpy(n,-omega,t,1,r,1)
     step = step + omega * phat!call daxpy(n,omega,phat,1,step,1)
     rsnrm = SQRT(SUM(r**2)) !dnrm2(n,r,1)
    ! c ------------------------------------------------------------------------ 
    ! c For printing:
     if (iplvl .ge. 4) then 
        write(ipunit,840) istb, rsnrm 
    840     format(5x,i4,'.b',3x,1pd10.3)
     endif
    ! c ------------------------------------------------------------------------
    ! c
    ! c ------------------------------------------------------------------------
    ! c Test for termination. 
    ! c ------------------------------------------------------------------------
     if (rsnrm .le. abstol) then 
        itrmks = 0
        go to 900
     endif
     if (istb .ge. iksmax_glob) then 
        itrmks = 3
        go to 900
     endif
    ! c ------------------------------------------------------------------------
    ! c If continuing, update and return to the top of the iteration loop. 
    ! c ------------------------------------------------------------------------
     rhomns = rho
     go to 100
    ! c ------------------------------------------------------------------------
    ! c Bottom of the iteration loop. 
    ! c ------------------------------------------------------------------------
    ! c
    ! c ------------------------------------------------------------------------
    ! c All returns made here.
    ! c ------------------------------------------------------------------------
    900  continue
    ! c ------------------------------------------------------------------------ 
    ! c For printing:
     if (iplvl .ge. 3) then 
        write(ipunit,*) 
        if (itrmks .ne. 1 .and. itrmks .ne. 2) then 
           write(ipunit,850) itrmks, rsnrm 
    850        format('nitstb:  itrmks =', i2, '   final lin. res. norm =', 1pd10.3)
        else
           write(ipunit,860) itrmks
    860        format('nitstb: itrmks:', i4) 
        endif
     endif
    ! c ------------------------------------------------------------------------
    ! c
    ! c ------------------------------------------------------------------------
    ! c If ijacv = -1, then restore it to the original value ijacv = 0. 
    ! c ------------------------------------------------------------------------
     if (ijacv_glob .eq. -1) ijacv_glob = 0 
end subroutine nitstb

subroutine nitbt(n, xcur, fcnrm, step, eta, xpls, fpls, fpnrm, oftjs, redfac, nfe, ibt, f, rpar, ipar, itrmbt)

   integer, intent(in) :: n
   integer, intent(in out) :: ipar(*), itrmbt, ibt, nfe
   double precision, intent(in) :: xcur(n), fcnrm
   double precision, intent(in out) :: eta, xpls(n), fpls(n), fpnrm, oftjs, redfac, rpar(*), step(n)
   procedure(funcInterface) :: f

    ! c ------------------------------------------------------------------------
    ! c
    ! c This is nitbt v0.3, the backtracking routine for the (inexact) Newton 
    ! c iterations. 
    ! c
    ! c ------------------------------------------------------------------------
    ! c 
    ! c Explanation: 
    ! c
    ! c  n       = dimension of the problem
    ! c
    ! c  xcur    = vector of length n, current approximate solution (input). 
    ! c
    ! c  fcnrm   = norm of f(xcur) 
    ! c
    ! c  step    = vector of length n, initial (trial) step on input, final 
    ! c            acceptable step on output. 
    ! c
    ! c  eta     = initial inexact Newton level on input, final inexact Newton 
    ! c            level on output. 
    ! c
    ! c  xpls    = vector of length n, next approximate solution on output. 
    ! c
    ! c  fpls    = vector of length n, value of f at xpls. 
    ! c
    ! c  fpnrm   = norm of f(xpls) 
    ! c
    ! c  oftjs   = original value of f(transpose)*Js. 
    ! c
    ! c  redfac   = scalar factor by which the original step is reduced on output. 
    ! c
    ! c  nfe     = number of function evaluations.
    ! c
    ! c  ibt     = number of backtracks on this call. 
    ! c
    ! c  ibtmax   = maximum allowable number of backtracks (step reductions) 
    ! c             per call to nitbt (default 10). 
    ! c
    ! c             USAGE NOTE: If ibtmax = -1, then backtracking is turned 
    ! c             off. In this case, the only function of this subroutine 
    ! c             is to update xpls, fpls, and fpnrm. 
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
    ! c  rpar    = real parameter/work array passed to the f and jacv routines. 
    ! c
    ! c  ipar    = integer parameter/work array passed to the f and jacv routines. 
    ! c
    ! c  dnorm   = norm routine, either user-supplied or blas dnrm2. 
    ! c
    ! c  itrmbt  = termination flag; values have the following meanings: 
    ! c              0 => normal termination: acceptable step found. 
    ! c              1 => acceptable step not found in ibtmax reductions.
    ! c              2 => error in evaluation of f.
    ! c 
    ! c ------------------------------------------------------------------------
    ! c
    ! c Subroutines required by this and all called routines: 
    ! c
    ! c    user supplied: f
    ! c
    ! c    nitsol routines: none 
    ! c
    ! c    blas routine: dscal
    ! c
    ! c    user supplied or blas: dnorm 
    ! c
    ! c    explanation: In nitsol, dnorm is set to either the blas 
    ! c    dnrm2 routine or the user-supplied usrnrm routine. 
    ! c
    ! c This subroutine called by: nitdrv
    ! c
    ! c Subroutines called by this subroutine: dnorm, dscal, f
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
     double precision :: t, theta
     integer :: i, itrmf

    ! c ------------------------------------------------------------------------ 
    ! c
    ! c ------------------------------------------------------------------------
    ! c Initialize.
    ! c ------------------------------------------------------------------------
     t = 1.d-4
     ibt = 0
     redfac = 1.d0
    ! c ------------------------------------------------------------------------
    ! c Backtracking loop. 
    ! c ------------------------------------------------------------------------
    100  continue
     do 110 i = 1, n
        xpls(i) = xcur(i) + step(i) 
    110  continue
     call f(n, xpls, fpls, rpar, ipar, itrmf) 
     if (itrmf .ne. 0) then
        itrmbt = 2
        go to 900
     endif
     nfe = nfe + 1
     fpnrm = SQRT(SUM(fpls**2))!dnrm2(n, fpls, 1) 
    ! c ------------------------------------------------------------------------
    ! c If t-condition is met or backtracking is turned off, return. 
    ! c ------------------------------------------------------------------------
     if (fpnrm .le. (1.d0 - t*(1.d0-eta))*fcnrm .or. ibtmax_glob .eq. -1) then 
        itrmbt = 0
        go to 900 
     endif
    ! c ------------------------------------------------------------------------
    ! c Otherwise, ... 
    ! c ------------------------------------------------------------------------
     ibt = ibt + 1
     if (ibt .gt. ibtmax_glob) then 
        itrmbt = 1
        go to 900
     endif
    ! c ------------------------------------------------------------------------
    ! c ... choose theta ...
    ! c ------------------------------------------------------------------------
     theta = -(oftjs*redfac)/(fpnrm**2 - fcnrm**2 - 2.d0*oftjs*redfac)
     if(theta .lt. thmin) theta = thmin
     if(theta .gt. thmax) theta = thmax
    ! c ------------------------------------------------------------------------
    ! c ... then reduce the step, increase eta, update redfac ... 
    ! c ------------------------------------------------------------------------
     step = theta * step!call dscal(n, theta, step, 1)
     eta = 1.d0 - theta*(1.d0 - eta)
     redfac = theta*redfac
    ! c ------------------------------------------------------------------------
    ! c ... and return to the top of the loop. 
    ! c ------------------------------------------------------------------------
    ! c ------------------------------------------------------------------------ 
    ! c For printing:
     if (iplvl .ge. 4) then 
        if (ibt .eq. 1) then 
           write(ipunit,*)
           write(ipunit,800)
           write(ipunit,*)
        endif
        write(ipunit,810) ibt, fpnrm, theta
     endif
    800  format('nitbt:  Step reduction no., trial F norm, current ','reduction factor')
    810  format(5x,i4,2(5x,1pd10.3))
    ! c ------------------------------------------------------------------------
     go to 100
    ! c ------------------------------------------------------------------------
    ! c All returns made here.
    ! c ------------------------------------------------------------------------
    900  continue
    ! c ------------------------------------------------------------------------ 
    ! c For printing:
     if (iplvl .ge. 3 .and. ibtmax_glob .ne. -1) then 
        write(ipunit,*) 
        if (ibt .eq. 0) then 
           write(ipunit,820) 
        else
           write(ipunit,830) ibt, redfac
        endif
     endif
    820  format( 'nitbt:  no. of step reductions. = 0')
    830  format( 'nitbt:  no. of step reductions. =', i2, '    total reduction factor =', 1pd10.3)

end subroutine nitbt

subroutine nitdrv(n, xcur, fcur, xpls, fpls, step, f, jacv, rpar, ipar, abstol, reltol, stptol, &
   iterm, nfe, njve, nrpre, nli, nni, nbt, rwork)

   integer, intent(in) :: n
   integer, intent(in out) :: ipar(*), nfe, njve, nrpre, nli, nni, nbt, iterm
   double precision, intent(in) :: abstol, stptol, reltol
   double precision, intent(in out) :: xcur(n), fcur(n), xpls(n), fpls(n), step(n), rpar(*), rwork(*)
   procedure(funcInterface) :: f
   procedure(jacvInterface) :: jacv

    ! c ------------------------------------------------------------------------
    ! c
    ! c This is nitdrv v0.3, the driver routine for the Newton iterative 
    ! c method.  
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
    ! c  xpls    = vector of length n, next (trial) approximate solution. 
    ! c            Also used as a work array where indicated.  
    ! c
    ! c  fpls    = vector of length n, value of f at xpls. Also used as a 
    ! c            work array where indicated.  
    ! c
    ! c  step    = vector of length n, (trial) step. 
    ! c
    ! c  f       = name of user-supplied subroutine for evaluating the function 
    ! c            the zero of which is sought; this routine has the form
    ! c
    ! c                 subroutine f(n, xcur, fcur, rpar, ipar, itrmf)
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
    ! c            and itrmjv is an integer termination flag. The meaning of 
    ! c            ijob is as follows: 
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
    ! c  ftol    = stopping tolerance on the f-norm.
    ! c
    ! c  stptol  = stopping tolerance on the steplength.
    ! c
    ! c  nnimax  = maximum allowable number of nonlinear iterations (default 200). 
    ! c
    ! c  ijacv   = flag for determining the method of J*v evaluation: 
    ! c              0 => finite-difference evaluation (default) 
    ! c              1 => analytic evaluation 
    ! c
    ! c  ikrysl  = flag for determining the Krylov solver: 
    ! c              0 => GMRES (default) 
    ! c              1 => BiCGSTAB
    ! c              2 => TFQMR
    ! c
    ! c            For brief descriptions of the solvers plus references, 
    ! c            see the subroutines nitgm, nitstb, and nittfq. 
    ! c
    ! c  kdmax   = maximum Krylov subspace dimension when GMRES is used 
    ! c            (default 20).
    ! c
    ! c  irpre   = flag for right preconditioning: 
    ! c              0 => no right preconditioning
    ! c              1 => right preconditioning
    ! c
    ! c  iksmax  = maximum allowable number of iterations per call to the Krylov 
    ! c            solver (default 1000). 
    ! c
    ! c  iresup  = residual update flag when GMRES is used; on GMRES restarts, 
    ! c            the residual is updated as follows: 
    ! c              0 => linear combination (default) 
    ! c              1 => direct evaluation
    ! c            The first is cheap (one n-vector saxpy) but may lose 
    ! c            accuracy with extreme residual reduction; the second 
    ! c            retains accuracy better but costs one J*v product per 
    ! c            restart. 
    ! c
    ! c  ifdord  = order of the finite-difference formula (sometimes) used when 
    ! c            ijacv = 0. When ijacv = 0, this is set to 1, 2, or 4 in nitsol; 
    ! c            otherwise it is irrelevant. With ijacv = 0, the precise meaning 
    ! c            is as follows: 
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
    ! c  ibtmax  = maximum allowable number of backtracks (step 
    ! c            reductions) per call to nitbt (default 10). 
    ! c
    ! c            USAGE NOTE: Backtracking can be turned off by setting 
    ! c            ibtmax = -1. Other negative values of ibtmax are not 
    ! c            valid. 
    ! c
    ! c
    ! c  ieta    = flag determining the forcing term eta as follows: 
    ! c               0 => (||fcur|| - ||fprev+Jprev*sprev||)/||fprev||
    ! c                    (default)
    ! c               1 => (||fcur||/||fprev||)**2
    ! c               2 => gamma*(||fcur||/||fprev||)**alpha  
    ! c                    for user-supplied gamma in (0,1] and alpha in (1,2] 
    ! c               3 => user-supplied eta in [0,1). 
    ! c                 3 => fixed (constant) eta in (0,1), either 0.1 (default) 
    ! c		       or specified by the user (see USAGE NOTE below) 
    ! c            Here, fcur = current f, fprev = previous f, etc. The Krylov 
    ! c            iterations are terminated when an iterate s satisfies 
    ! c            an inexact Newton condition ||F + J*s|| .le. eta*||F||.
    ! c
    ! c            USAGE NOTE: If ieta = 2, then alpha and gamma must be set 
    ! c            in common block nitparam.h as described below. If 
    ! c	     ieta = 3, then the desired constant eta may be similarly 
    ! c	     set in nitparam.h if a value other than the default of 
    ! c	     0.1 is desired. 
    ! c               
    ! c            The first three expressions above are from S. C. Eisenstat 
    ! c            and H. F. Walker, "Choosing the forcing terms in an inexact 
    ! c            Newton method", SIAM J. Scientific Computing, 17 (1996), 
    ! c            pp. 16--32. (They may be modified according to certain 
    ! c            safeguards below.) The first gives convergence that is 
    ! c            q-superlinear and of r-order (1+sqrt(5))/2. The second gives 
    ! c            convergence that is r-quadratic and of q-order p for every p 
    ! c            in [1,2). The third gives convergence that is of q-order alpha 
    ! c            when gamma < 1 and, when gamma = 1, of r-order alpha and 
    ! c            q-order p for every p in [1,alpha). The fourth gives q-linear 
    ! c            convergence with asymptotic rate constant alpha in a certain 
    ! c            norm; see R. S. Dembo, S. C. Eisenstat, and T. Steihaug, 
    ! c            "Inexact Newton methods", SIAM J. Numer. Anal., 18 (1982), 
    ! c            pp. 400-408. 
    ! c
    ! c            Of these four choices, the 1st is usually satisfactory, 
    ! c            the 2nd or 3rd is sometimes preferred, and the 4th may be 
    ! c            useful in some situations, e.g., it may be desirable to 
    ! c            choose a fairly large fixed eta in (0,1), such as eta = .1, 
    ! c            when numerical inaccuracy prevents the Krylov solver 
    ! c            from obtaining much residual reduction. 
    ! c               
    ! c  iterm   = termination flag; values have the following meanings: 
    ! c              0 => normal termination: ||F||.le.ftol or ||step||.le.stptol.
    ! c              1 => nnimax nonlinear iterations reached without success. 
    ! c              2 => failure to evaluate F.
    ! c              3 => in nitjv, J*v failure. 
    ! c              4 => in nitjv, P(inverse)*v failure. 
    ! c              5 => in nitdrv, insufficient initial model norm reduction 
    ! c                   for adequate progress. NOTE: This can occur for several 
    ! c                   reasons; examine itrmks on return from the Krylov 
    ! c                   solver for further information. (This will be printed out 
    ! c                   if iplvl .ge. 3, see the discussion of optional 
    ! c                   common blocks below). 
    ! c              6 => in nitbt, failure to reach an acceptable step through 
    ! c                   backtracking. 
    ! c 
    ! c  nfe     = number of function evaluations.
    ! c
    ! c  njve    = number of J*v evaluations. 
    ! c
    ! c  nrpre   = number of P(inverse)*v evaluations.
    ! c
    ! c  nli     = number of linear iterations.
    ! c
    ! c  nni     = number of nonlinear iterations.
    ! c
    ! c  nbt     = number of backtracks. 
    ! c
    ! c  rwork   = real work vector for use by the Krylov solver. It is passed 
    ! c            in as the tail of the rwork vector in nitsol. On input to 
    ! c            nitsol, it should have length as follows: 
    ! c
    ! c             solver    rwork length
    ! c             GMRES     n*(kdmax+5)+kdmax*(kdmax+3), where kdmax is the 
    ! c                       maximum Krylov subspace dimension, either the 
    ! c                       default value of 20 or another value specified 
    ! c                       by the user. 
    ! c             BiCGSTAB  11*n
    ! c             TFQMR     14*n
    ! c
    ! c  dinpr   = inner-product routine, either user-supplied or BLAS ddot. 
    ! c
    ! c  dnorm   = norm routine, either user-supplied or BLAS dnrm2. 
    ! c
    ! c ------------------------------------------------------------------------
    ! c
    ! c Optional common blocks: 
    ! c
    ! c These can be used to control printing of diagnostic information by nitsol, 
    ! c to pass information about the nonlinear iterations to jacv or other user 
    ! c subroutines, or to control the default behavior of the nonlinear iterations. 
    ! c
    ! c For controlling printing of diagnostic information: 
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
    ! c     ipunit = printout unit number, e.g., ipunit = 6 => standard output. 
    ! c
    ! c For passing information about the nonlinear iterations to user-supplied 
    ! c subroutines: 
    ! c
    !      include 'nitinfo.h'
    ! c
    ! c If information on the current state of the nonlinear iteration is
    ! c desired in a user-supplied subroutine (for example, deciding 
    ! c whether to update a preconditioner), include this common block
    ! c in the subroutine. The variables are as follows: 
    ! c
    ! c     instep - inexact Newton step number. 
    ! c
    ! c    newstep - set to 0 at the beginning of an inexact Newton step.
    ! c              This may be checked in a user-supplied jacv to decide
    ! c              whether to update the preconditioner.  If you test on
    ! c              newstep .eq. 0 to determine whether to take some 
    ! c              special action at the beginning of a nonlinear iteration, 
    ! c              you must also set newstep to some nonzero value to
    ! c              subsequently avoid taking that action unnecessarily. 
    ! c
    ! c    krystat - status of the Krylov iteration; same as itrmks (see 
    ! c              the nitsol documentation). 
    ! c
    ! c    avrate  - average rate of convergence of the Krylov solver during
    ! c              the previous inexact Newton step.  This may be checked
    ! c              in a user-supplied jacv to decide when to update the
    ! c              preconditioner.
    ! c
    ! c    fcurnrm - ||f(xcur)||. 
    ! c
    ! c
    ! c  For controlling the default behavior of the nonlinear iterations:
    ! c
    !      include 'nitparam.h'

    ! c nitparam contains some parameters that control the nonlinear
    ! c iterations.  In some cases, the default values reflect prevailing  
    ! c practice; in other cases, they are chosen to produce good 
    ! c average-case behavior.  To change the default values, include this 
    ! c common block in the main program and set the desired variables 
    ! c according to the following:
    ! c
    ! c    choice1_exp -  parameter used in the update of the forcing term 
    ! c                   eta when ieta = 0 (default).  This is the exponent
    ! c                   for determining the etamin safeguard.  The default
    ! c                   value is choice1_exp = (1+sqrt(5))/2.  A larger
    ! c                   value will allow eta to decrease more rapidly,
    ! c                   while a smaller value will result in a larger 
    ! c                   value for the safeguard. 
    ! c
    ! c    choice2_exp  - parameter used in the update of the forcing term 
    ! c                   eta when ieta = 2.  This is the exponent alpha 
    ! c		    in the expression gamma*(||fcur||/||fprev||)**alpha; 
    ! c		    it is also used to determine the etamin safeguard.  
    ! c		    The default value is 2.0. Valid values are in the 
    ! c		    range (1.0, 2.0].
    ! c
    ! c    choice2_coef - parameter used in the update of the forcing term eta 
    ! c                   when ieta = 2.  This is the coefficient gamma used 
    ! c		    in the expression gamma*(||fcur||/||fprev||)**alpha;
    ! c                   it is also used to determine the etamin safeguard.
    ! c                   The default value is 1.0. Valid values are in the 
    ! c		    range (0.0, 1.0]. 
    ! c
    ! c    eta_cutoff   - parameter used to determine when to disable 
    ! c                   safeguarding the update of the forcing term.  It
    ! c                   only has meaning when ieta .ne. 3.  The default
    ! c                   value is 0.1.  A value of 0.0 will enable 
    ! c		    safeguarding always; a value of 1.0 will disable 
    ! c		    safeguarding always. 
    ! c
    ! c    etamax       - parameter used to provide an upper bound on the 
    ! c		    forcing terms when input(10) .ne. 3. This is 
    ! c		    necessary to ensure convergence of the inexact Newton 
    ! c		    iterates and is imposed whenever eta would otherwise 
    ! c		    be too large. (An overly large eta can result from 
    ! c		    the updating formulas when input(10) .ne. 3 or from 
    ! c                   safeguarding when the previous forcing term has been 
    ! c		    excessively increased during backtracking.) The 
    ! c		    default value of etamax is 1.0 - 1.e-4.  When 
    ! c		    backtracking occurs several times during a nonlinear 
    ! c		    solve the forcing term can remain near etamax for several
    ! c                   nonlinear steps and cause the nonlinear iterations
    ! c                   to nearly stagnate.  In such cases a smaller value of 
    ! c                   etamax may prevent this.  Valid values are in the 
    ! c                   range (0.0, 1.0).
    ! c
    ! c    etafixed     - this is the user-supplied fixed eta when ieta = 3.
    ! c                   The  default value is etafixed = 0.1.  Valid values
    ! c                   are in the range (0.0,1.0).
    ! c
    ! c    thmin        - when backtracking occurs, this is the smallest
    ! c                   reduction factor that will be applied to the current
    ! c                   step in a single backtracking reduction.  The default
    ! c                   value is 0.1.  Valid  values are in the range
    ! c                   [0.0, thmax].
    ! c
    ! c    thmax        - when backtracking occurs, this is the largest
    ! c                   reduction factor that will be applied to the current
    ! c                   step in a single backtracking reduction.  The default
    ! c                   value is 0.5.  Valid values are in the range
    ! c                   [thmin, 1.0).
    ! c
    ! c  The values in this common block are not checked here.  We assume
    ! c  that if you call nitdrv directly you know what you are doing.
    ! c
    ! c ------------------------------------------------------------------------
    ! c
    ! c Subroutines required by this and all called routines:
    ! c
    ! c    user supplied: f, jacv
    ! c
    ! c    nitsol routines: nitbd.f, nitbt, nitgm, nitjv, nitstb, nittfq, 
    ! c
    ! c    lapack routines: dlaic1, dlamch
    ! c
    ! c    blas routines: daxpy, dcopy, dscal, dswap 
    ! c
    ! c    user supplied or BLAS (see below): dinpr, dnorm 
    ! c
    ! c    Explanation: The call to nitsol specifies dinpr and dnorm as 
    ! c    either user-supplied inner-product and norm routines or the 
    ! c    BLAS ddot and dnrm2 routines. 
    ! c
    ! c This subroutine called by: nitsol
    ! c
    ! c Subroutines called by this subroutine: dcopy, dinpr, dlamch, dnorm, f,
    ! c    nitbt, nitgm, nitstb, nittfq 
    ! c
    ! c ------------------------------------------------------------------------
    ! c
    ! c Remaining declarations: 
    ! c
    ! c NOTE: In nitinfo.h, instep, newstep, krystat are declared integer, 
    ! c and avrate and fcurnrm are declared double precision. 
    ! c 
     double precision eta, etamin, fcnrm, flmnrm, fpnrm, oftjs, oftlm, redfac, rsnrm, stpnrm, temp, alpha, gamma, ftol
     integer ibt, itrmbt, itrmf, itrmks, lrr, lsvbig, lsvsml, lvv, lw, lr, ld, lrcgs, lrtil, lp, lphat, lq, lu, lv, lt, lrwork, ly, kdmaxp1



     integer :: infe, injve, inrpre, inli

    ! c ------------------------------------------------------------------------
    ! c Initialize.
    ! c ------------------------------------------------------------------------
     if (ieta_glob .eq. 0) alpha = choice1_exp
     if (ieta_glob .eq. 2) alpha = choice2_exp
     if (ieta_glob .eq. 2) gamma = choice2_coef
     if (ieta_glob .eq. 3) then 
        eta = etafixed
     else
        eta = .5d0
     endif
     nfe = 0
     njve = 0
     nrpre = 0
     nli = 0
     nni = 0
     nbt = 0
     avrate = 1.0d0
    ! c ------------------------------------------------------------------------
    ! c Evaluate f at initial x and initialize eta. 
    ! c ------------------------------------------------------------------------
     call f(n, xcur, fcur, rpar, ipar, itrmf)
     if ( itrmf .ne. 0 ) then
        iterm = 2
        go to 900
     endif
     nfe = nfe + 1
     fcnrm = SQRT(SUM(fcur**2))!dnrm2(n, fcur, 1) 
     ftol = abstol + reltol*fcnrm
     
    ! c ------------------------------------------------------------------------ 
    ! c For printing:
     if (iplvl .ge. 1) then 
        write(ipunit,*) 
        write(ipunit,*) 
        write(ipunit,*) 'nitdrv:  Beginning nonlinear iterations.'
     endif
    ! c ------------------------------------------------------------------------
    ! c
    ! c ------------------------------------------------------------------------
    ! c Nonlinear iteration loop.
    ! c      tstrt = etime(dummy)
    ! c ------------------------------------------------------------------------
    100  continue
    ! c ------------------------------------------------------------------------ 
    ! c For printing:
    ! c      if (iplvl .eq. 0) then
    ! c         if ( nni .eq. 0 ) then
    ! c            write(ipunit,799) nni, njve, 0.0d0, fcnrm
    ! c         else 
    ! c            write(ipunit,799) nni, njve, etime(dummy)-tstrt, fcnrm
    ! c         endif
    ! c      endif
    ! c 799  format( 2(2x,i4),2(2x,1pe12.3) )
     if (iplvl .ge. 1) then 
        write(ipunit,*)
        if (iplvl .eq. 1) write(ipunit,800) nni, fcnrm
        if (iplvl .ge. 2) write(ipunit,810) nni, fcnrm 
     endif
    800  format('    It. no.', i4, '      F norm =', 1pd12.3)
    810  format('--- It. no.', i4, '      F norm =', 1pd12.3,' ---------------------------')
     if (iplvl .ge. 2) then 
        write(ipunit,*) 
        write(ipunit,820) nfe, njve, nrpre, nli
        infe = nfe
        injve = njve
        inrpre = nrpre 
        inli = nli
     endif
    820  format('  Initial totals:  nfe =', i4, '    njve =', i4, '    nrpre =', i4, '    nli =', i4)
    ! c ------------------------------------------------------------------------
    ! c
    ! c------------------------------------------------------------------------
    ! c Test for stopping. 
    ! c------------------------------------------------------------------------
     if (fcnrm .le. ftol) then 
        iterm = 0
        go to 900
     endif
     if (nni .gt. 0 .and. stpnrm .le. stptol .and. itrmks .eq. 0) then 
        iterm = 0
        go to 900
     endif
     if (nni .ge. nnimax_glob) then 
        iterm = 1
        go to 900
     endif
    ! c ------------------------------------------------------------------------
    ! c Compute the (trial) inexact Newton step with the Krylov solver. 
    ! c ------------------------------------------------------------------------
    ! c Update data in nitinfo to mark the start of a new inexact Newton step.
    ! c ------------------------------------------------------------------------
     newstep = 0
     instep = nni
     fcurnrm = fcnrm
     SELECT CASE (ikrysl_glob)
     CASE(0)
      ! c ------------------------------------------------------------------------
      ! c If ikrysl = 0, apply GMRES, using fpls as a work array. 
      ! c ------------------------------------------------------------------------
      kdmaxp1 = kdmax_glob + 1 
      lvv = 1
      lrr = lvv + n*kdmaxp1
      lsvbig = lrr + kdmax_glob*kdmax_glob 
      lsvsml = lsvbig + kdmax_glob
      lw = lsvsml + kdmax_glob
      call nitgm(n, xcur, fcur, fcnrm, step, eta, f, jacv,rpar, ipar, nfe, njve, &
      nrpre, nli, kdmaxp1, rwork(lvv), rwork(lrr), rwork(lsvbig), rwork(lsvsml), rwork(lw), fpls, rsnrm, itrmks)
     CASE(1)
      ! c ------------------------------------------------------------------------
      ! c If ikrysl = 1, apply BiCGSTAB, using fpls as a work array. 
      ! c ------------------------------------------------------------------------
        lr = 1
        lrtil = lr + n
        lp = lrtil + n
        lphat = lp + n
        lv = lphat + n
        lt = lv + n
        lrwork = lt + n
        call nitstb (n, xcur, fcur, fcnrm, step, eta, f, jacv,rpar, &
        ipar, nfe, njve, &
        nrpre, nli, rwork(lr), rwork(lrtil), rwork(lp), &
        rwork(lphat), rwork(lv), rwork(lt), rwork(lrwork), fpls, rsnrm, itrmks)
     CASE(2)
      ! c ------------------------------------------------------------------------
      ! c If ikrysl = 2, apply TFQMR 
      ! c ------------------------------------------------------------------------
         lr = 1
         lrcgs = lr + n
         lrtil = lrcgs + n
         ld = lrtil + n
         lp = ld + n
         lq = lp + n
         lu = lq + n
         lv = lu + n
         ly = lv + n
         lrwork = ly + n
         call nittfq (n, xcur, fcur, fcnrm, step, eta, f, jacv,rpar, &
         ipar, nfe, njve, nrpre, nli, &
         rwork(lr), rwork(lrcgs), rwork(lrtil), rwork(ld), rwork(lp), rwork(lq), &
         rwork(lu), rwork(lv), rwork(ly), rwork(lrwork), fpls, rsnrm, itrmks)
     END SELECT
     
    ! c ------------------------------------------------------------------------
    ! c  Set values in nitinfo that reflect state of iterative solver.
    ! c ------------------------------------------------------------------------
     krystat = itrmks
     avrate = (rsnrm/fcnrm)**(1.0d0/dble(nli))
     ! USE avrate as additional constraint for solutions converging too slowly or stalling
   !   if (avrate > 0.98) then
   !       iterm = 5
   !       go to 900
   !   end if
     
    ! c ------------------------------------------------------------------------
    ! c Check itrmks and decide whether to terminate or continue: 
    ! c        0 => continue, inexact Newton condition successfully met 
    ! c   1 or 2 => terminate unconditionally, J*v or P(inverse)*v failure) 
    ! c   .ge. 3 => terminate if the model norm has increased or if reduction at 
    ! c             the current rate would at best require more than 1000 time 
    ! c             the maximum remaining number of nonlinear iterations. 
    ! c ------------------------------------------------------------------------
     if (itrmks .eq. 1 .or. itrmks .eq. 2) then
        iterm = itrmks + 2
        go to 900
     endif
     if (itrmks .ge. 3) then 
        if (rsnrm/fcnrm .ge. 1.d0) then 
           iterm = 5
           go to 900
        else
           temp = dlog(ftol/fcnrm)/dlog(rsnrm/((1.d0 + 10.d0*epsmach)*fcnrm))
           if (temp .gt. 1000.d0*dfloat(nnimax_glob - nni)) then 
              iterm = 5
              go to 900
           endif
        endif
     endif
    ! c ------------------------------------------------------------------------
    ! c Compute the original value of f(transpose)*Js for backtracking; the 
    ! c original value of f(transpose)*(linear model) is also computed for 
    ! c later use. NOTE: The first n components of rwork contain the residual  
    ! c vector for the Newton equation, which is -(linear model). 
    ! c ------------------------------------------------------------------------
     oftlm = - SUM(fcur*rwork(1:n))!ddot(n, fcur, 1, rwork, 1)
     oftjs = oftlm - fcnrm**2
    ! c ------------------------------------------------------------------------
    ! c Determine an acceptable step via backtracking. 
    ! c ------------------------------------------------------------------------
     call nitbt(n, xcur, fcnrm, step, eta, xpls, fpls, fpnrm, oftjs, redfac, nfe, ibt, f, rpar, ipar, itrmbt)
     if (itrmbt .eq. 1) then
        iterm = 6
        go to 900
     else if (itrmbt .eq. 2) then
        iterm = 2
        go to 900
     endif
     nbt = nbt + ibt
    ! c ------------------------------------------------------------------------
    ! c Set eta for next iteration. 
    ! c ------------------------------------------------------------------------
     if (ieta_glob .eq. 0) then 
        etamin = eta**alpha
        temp = 1.d0 - redfac
        flmnrm = dsqrt((temp*fcnrm)**2 + 2.d0*redfac*temp*oftlm + (redfac*rsnrm)**2)
        eta = dabs(fpnrm - flmnrm)/fcnrm 
     endif
     if (ieta_glob .eq. 1) then 
        etamin = eta**2
        eta = (fpnrm/fcnrm)**2
     endif
     if (ieta_glob .eq. 2) then 
        etamin = gamma*eta**alpha
        eta = gamma*(fpnrm/fcnrm)**alpha
     endif
     if (ieta_glob .ne. 3) then 
        if (etamin .le. eta_cutoff) etamin = 0.d0
    ! chfw         if (etamin .le. 1.d-1) etamin = 1.d-1
        if (eta .lt. etamin) eta = etamin 
        if (eta .gt. etamax) eta = etamax 
        if (eta*fpnrm .le. 2.d0*ftol) eta = (.8d0*ftol)/fpnrm
     endif
     if (ieta_glob .eq. 3) eta = etafixed
    ! c ------------------------------------------------------------------------
    ! c Update xcur, fcur, fcnrm, stpnrm, nni for next iteration.
    ! c ------------------------------------------------------------------------
     xcur = xpls!call dcopy(n, xpls, 1, xcur, 1)
     fcur = fpls!call dcopy(n, fpls, 1, fcur, 1)
     fcnrm = fpnrm
     stpnrm = SQRT(SUM(step**2))!dnrm2(n, step, 1)
     nni = nni + 1
    ! c ------------------------------------------------------------------------ 
    ! c For printing:
     if (iplvl .ge. 2) then 
        infe = nfe - infe 
        injve = njve - injve 
        inrpre = nrpre  - inrpre 
        inli = nli - inli 
        if (ieta_glob .gt. 0) then 
           temp = 1.d0 - redfac
           flmnrm = dsqrt((temp*fcnrm)**2 + 2.d0*redfac*temp*oftlm + (redfac*rsnrm)**2)
        endif
        write(ipunit,*) 
        write(ipunit,830) infe, injve, inrpre, inli, stpnrm, flmnrm 
     endif
    830  format('  At this step:   nfe =', i4, '    njve =', i4, '    nrpre =', i4, '    nli =', i4, /, 17x, &
    ' step norm =', 1pd10.3, 5x, 'final lin. model norm =', 1pd10.3)
    ! c ------------------------------------------------------------------------
    ! c
    ! c ------------------------------------------------------------------------
    ! c Return to top of loop for next iteration.
    ! c ------------------------------------------------------------------------
     go to 100
    ! c ------------------------------------------------------------------------
    ! c All returns made here.
    ! c ------------------------------------------------------------------------
    900  continue
    ! c ------------------------------------------------------------------------ 
    ! c For printing:
     if (iplvl .ge. 1) then 
        write(ipunit,*) 
        write(ipunit,*) 
        write(ipunit,*) 'nitdrv:  Terminating nonlinear iterations.'
     endif
end subroutine nitdrv

subroutine nitsol(n, x, f, jacv, abstol, reltol, stptol, info, rpar, ipar, iterm)


   integer, intent(in) :: n
   integer, intent(in out) :: info(6), ipar(*), iterm
   double precision, intent(in) :: abstol, stptol, reltol
   double precision, intent(in out) :: x(n), rpar(*)
   procedure(funcInterface) :: f
   procedure(jacvInterface) :: jacv

    ! c
    ! c ------------------------------------------------------------------------
    ! c
    ! c This is nitsol v0.3. The basic algorithm is an inexact Newton method with 
    ! c a backtracking globalization; the model is Algorithm INB of S. C. Eisenstat 
    ! c and H. F. Walker, "Globally convergent inexact Newton methods", SIAM J. 
    ! c Optimization, 4 (1994), pp. 393--422. Initial inexact Newton steps are 
    ! c obtained by approximately solving the Newton equation with a transpose-free
    ! c Krylov subspace method; the current choices are GMRES, BiCGSTAB, and 
    ! c TFQMR. Jacobian-vector products are evaluated either by finite-difference 
    ! c approximation or a user-supplied analytic-evaluation subroutine. An option 
    ! c is provided for user-supplied right preconditioning. Left preconditioning 
    ! c is not explicitly included as an option, but the user may provide this 
    ! c in the subroutines for evaluating the function and Jacobian-vector 
    ! c products. Various algorithmic options can be selected through the input 
    ! c vector. Optional common blocks are also available for printing diagnostic 
    ! c information, passing information about the nonlinear iterations to user 
    ! c subroutines, and controlling the behavior of the nonlinear iterations. 
    ! c Summary statistics are provided by the info vector on output. 
    ! c
    ! c This is the interface subroutine, which calls the driver subroutine nitdrv. 
    ! c
    ! c ------------------------------------------------------------------------
    ! c
    ! c Explanation: 
    ! c
    ! c  n      = dimension of the problem. 
    ! c
    ! c  x      = vector of length n, initial guess on input and final approximate 
    ! c           solution on output. 
    ! c
    ! c  f      = name of user-supplied subroutine for evaluating the function 
    ! c           the zero of which is sought; this routine has the form 
    ! c
    ! c                 subroutine f(n, xcur, fcur, rpar, ipar, itrmf)
    ! c
    ! c           where xcur is the array containing the current x value, fcur 
    ! c           is f(xcur) on output, rpar and ipar are, respectively, real 
    ! c           and integer parameter/work arrays for use by the subroutine,
    ! c           and itrmf is an integer termination flag.  The meaning of
    ! c           itrmf is as follows:
    ! c             0 => normal termination; desired function value calculated.
    ! c             1 => failure to produce f(xcur).
    ! c 
    ! c  jacv   = name of user-supplied subroutine for optionally evaluating J*v 
    ! c           or P(inverse)*v, where J is the Jacobian of f and P is a 
    ! c           right preconditioning operator. If neither analytic J*v 
    ! c           evaluations nor right preconditioning is used, this can 
    ! c           be a dummy subroutine; if right preconditioning is used but 
    ! c           not analytic J*v evaluations, this need only evaluate 
    ! c           P(inverse)*v. The form is 
    ! c
    ! c           subroutine jacv(n, xcur, fcur, ijob, v, z, rpar, ipar, itrmjv)
    ! c
    ! c           where xcur and fcur are vectors of length n containing the 
    ! c           current x and f values, ijob is an integer flag indicating 
    ! c           which product is desired, v is a vector of length n to be 
    ! c           multiplied, z is a vector of length n containing the desired 
    ! c           product on output, rpar and ipar are, respectively, real 
    ! c           and integer parameter/work arrays for use by the subroutine, 
    ! c           and itrmjv is an integer termination flag. The meaning of 
    ! c           ijob is as follows: 
    ! c             0 => z = J*v
    ! c             1 => z = P(inverse)*v 
    ! c           The meaning of itrmjv is as follows:
    ! c             0 => normal termination; desired product evaluated. 
    ! c             1 => failure to produce J*v.
    ! c             2 => failure to produce P(inverse)*v. 
    ! c           This subroutine is called only from nitjv, and is always 
    ! c           called with v .ne. z. 
    ! c
    ! c  ftol   = stopping tolerance on the f-norm.
    ! c
    ! c  stptol = stopping tolerance on the steplength.
    ! c
    ! c  input  = integer vector of length 10 containing various user-specified 
    ! c           inputs; see below. 
    ! c
    ! c  info   = integer vector of length 6 containing various outputs; 
    ! c           see below. 
    ! c
    ! c  rwork  = real work vector with length as follows: 
    ! c
    ! c             solver    rwork length
    ! c
    ! c             GMRES     n*(kdmax+5)+kdmax*(kdmax+3), where kdmax is the 
    ! c                       maximum Krylov subspace dimension, either the 
    ! c                       default value of 20 or another value specified 
    ! c                       by the user (see input(4) below). 
    ! c
    ! c             BiCGSTAB  11*n
    ! c
    ! c             TFQMR     14*n
    ! c
    ! c           The final f-value is contained in the first n components of 
    ! c           rwork on return. 
    ! c
    ! c  rpar   = real parameter/work array passed to the f and jacv routines. 
    ! c
    ! c  ipar   = integer parameter/work array passed to the f and jacv routines. 
    ! c               
    ! c  iterm   = termination flag; values have the following meanings: 
    ! c             -k => illegal value in input(k). 
    ! c              0 => normal termination: ||F||.le.ftol or ||step||.le.stptol.
    ! c              1 => nnimax nonlinear iterations reached without success. 
    ! c              2 => failure to evaluate F.
    ! c              3 => in nitjv, J*v failure. 
    ! c              4 => in nitjv, P(inverse)*v failure. 
    ! c              5 => in nitdrv, insufficient initial model norm reduction 
    ! c                   for adequate progress. NOTE: This can occur for several 
    ! c                   reasons; examine itrmks on return from the Krylov 
    ! c                   solver for further information. (This will be printed out 
    ! c                   if iplvl .ge. 3; see the discussion of optional 
    ! c                   common blocks below). 
    ! c              6 => in nitbt, failure to reach an acceptable step through 
    ! c                   backtracking. 
    ! c
    ! c  dinpr  = name of user-supplied function for calculating vector inner
    ! c           products.  This function must have the form
    ! c
    ! c              xdoty = dinpr( n, x, sx, y, sy )
    ! c
    ! c           where n is the length of the vectors, x and y are the
    ! c           starting locations of the vectors, and sx (sy) is the stride
    ! c           in memory between consecutive elements of x (y).  This is the
    ! c           same signature as the BLAS routine ddot; if the Euclidean
    ! c           inner product is desired the user can link to a local BLAS
    ! c           library and provide the name ddot to nitsol.  dinpr must be
    ! c           declared as an external function that returns a double
    ! c           precision value in the calling program.
    ! c
    ! c  dnorm  = name of user-supplied function for calculating vector norms.
    ! c           This function must have the form
    ! c
    ! c              xnorm = dnorm( n, x, sx )
    ! c
    ! c           where n is the length of the vector, x is the starting
    ! c           location of the vector, and sx is the stride in memory
    ! c           between consecutive elements of x.  This is the same
    ! c           signature as the BLAS routine dnrm2; if the Euclidean
    ! c           norm is desired the user can link to a local BLAS library
    ! c           and provide the name dnrm2 to nitsol.  dnorm must be
    ! c           declared as an external function that returns a double
    ! c           precision value in the calling program.
    ! c
    ! c ------------------------------------------------------------------------
    ! c 
    ! c Further explanation of input: 
    ! c
    ! c This array allows the user to specify various options. It should be 
    ! c declared an integer vector of length 11 in the calling program. To 
    ! c specify an option, set the appropriate input component to the desired 
    ! c value according to the specifications below. 
    ! c
    ! c USAGE NOTE: Setting a particular input component to zero gives the 
    ! c default option for that component in all cases. 
    ! c
    ! c The first five input components are things that every user might wish 
    ! c to modify; the remainder will usually be of interest only to more 
    ! c experienced users. 
    ! c
    ! c Optional every-user input:
    ! c
    ! c    input(1) = nnimax = maximum number of nonlinear iterations (default 200).
    ! c 
    ! c    input(2) = ijacv = flag for determining the method of J*v evaluation:
    ! c                 0 => finite-difference evaluation (default) 
    ! c                 1 => analytic evaluation
    ! c
    ! c    input(3) = ikrysl = flag for determining the Krylov solver: 
    ! c                 0 => GMRES (default)
    ! c                 1 => BiCGSTAB
    ! c                 2 => TFQMR
    ! c
    ! c               For brief descriptions of the solvers plus references, 
    ! c               see the subroutines nitgm, nitstb, and nittfq. 
    ! c
    ! c    input(4) = kdmax = maximum Krylov subspace dimension when GMRES is used 
    ! c               (default 20). 
    ! c
    ! c    input(5) = irpre = flag for right preconditioning: 
    ! c                 0 => no right preconditioning
    ! c                 1 => right preconditioning
    ! c
    ! c Optional experienced user input:
    ! c
    ! c    input(6) = iksmax = maximum allowable number of iterations per call 
    ! c               to the Krylov solver routine (default 1000). 
    ! c
    ! c    input(7) = iresup = residual update flag when GMRES is used; on 
    ! c               restarts, the residual is updated as follows: 
    ! c                 0 => linear combination (default) 
    ! c                 1 => direct evaluation
    ! c               The first is cheap (one n-vector saxpy) but may lose 
    ! c               accuracy with extreme residual reduction; the second 
    ! c               retains accuracy better but costs one J*v product per 
    ! c               restart. 
    ! c
    ! c    input(8) = ifdord = order of the finite-difference formula (sometimes) 
    ! c               used when input(2) = ijacv = 0. When input(2) = ijacv = 0, 
    ! c               this must be 0, 1, 2, or 4 on input; otherwise, it is 
    ! c               irrelevant. With input(2) = ijacv = 0, the precise 
    ! c               meaning is as follows: 
    ! c
    ! c               If GMRES is used, then ifdord matters only if input(7) = 
    ! c               iresup = 1, in which case it determines the order of 
    ! c               the finite-difference formula used in evaluating the 
    ! c               initial residual at each GMRES restart (default 2); if 
    ! c               ifdord = 0 on input, then it is set to 2 below. NOTE: This 
    ! c               only affects initial residuals at restarts; first-order 
    ! c               differences are always used within each GMRES cycle. Using 
    ! c               higher-order differences at restarts only should give 
    ! c               the same accuracy as if higher-order differences were 
    ! c               used throughout; see K. Turner and H. F. Walker, "Efficient 
    ! c               high accuracy solutions with GMRES(m)," SIAM J. Sci. 
    ! c               Stat. Comput., 13 (1992), pp. 815--825. 
    ! c               
    ! c               If BiCGSTAB or TFQMR is used, then ifdord determines the 
    ! c               order of the finite-difference formula used at each 
    ! c               iteration (default 1); if ifdord = 0 on input, then it 
    ! c               is set to 1 below. 
    ! c
    ! c    input(9) = ibtmax = maximum allowable number of backtracks (step 
    ! c               reductions) per call to nitbt (default 10). 
    ! c
    ! c               USAGE NOTE: Backtracking can be turned off by setting 
    ! c		ibtmax = -1. Other negative values of ibtmax are not 
    ! c               valid. 
    ! c
    ! c    input(10) = ieta = flag determining the forcing term eta as follows: 
    ! c                 0 => abs( ||fcur|| - ||fprev+Jprev*sprev|| )/||fprev||
    ! c                      (default) 
    ! c                 1 => (||fcur||/||fprev||)**2 
    ! c                 2 => gamma*(||fcur||/||fprev||)**alpha 
    ! c                      for user-supplied gamma in (0,1] and alpha in (1,2] 
    ! c                 3 => fixed (constant) eta in (0,1), either 0.1 (default) 
    ! c		       or specified by the user (see USAGE NOTE below) 
    ! c               Here, fcur = current f, fprev = previous f, etc. The Krylov 
    ! c               iterations are terminated when an iterate s satisfies 
    ! c               an inexact Newton condition ||F + J*s|| .le. eta*||F||.
    ! c
    ! c               USAGE NOTE: If input(10) = ieta = 2, then alpha and gamma 
    ! c               must be set in common block nitparam.h as described below. 
    ! c		If input(10) = ieta = 3, then the desired constant eta may 
    ! c		be similarly set in nitparam.h if a value other than the 
    ! c		default of 0.1 is desired. 
    ! c               
    ! c               The first three expressions above are from S. C. Eisenstat 
    ! c               and H. F. Walker, "Choosing the forcing terms in an inexact 
    ! c               Newton method", SIAM J. Scientific Computing, 17 (1996), 
    ! c               pp. 16--32. (They may be modified according to certain 
    ! c               safeguards in subroutine nitdrv.) The first gives convergence 
    ! c               that is q-superlinear and of r-order (1+sqrt(5))/2. The 
    ! c               second gives convergence that is r-quadratic and of q-order 
    ! c               p for every p in [1,2). The third gives convergence that is 
    ! c               of q-order alpha when gamma < 1 and, when gamma = 1, of 
    ! c               r-order alpha and q-order p for every p in [1,alpha). The 
    ! c               fourth gives q-linear convergence with asymptotic rate 
    ! c               constant eta in a certain norm; see R. S. Dembo, S. C. 
    ! c		Eisenstat, and T. Steihaug, "Inexact Newton methods", 
    ! c               SIAM J. Numer. Anal., 18 (1982), pp. 400-408. 
    ! c
    ! c               Of these four choices, the 1st is usually satisfactory, 
    ! c               the 2nd or 3rd is sometimes preferred, and the 4th may be 
    ! c               useful in some situations, e.g., it may be desirable to 
    ! c               choose a fairly large fixed eta in (0,1), such as eta = .1, 
    ! c               when numerical inaccuracy prevents the Krylov solver 
    ! c               from obtaining much residual reduction. 
    ! c
    ! c ------------------------------------------------------------------------
    ! c
    ! c Further explanation of info: On output, the components of info are 
    ! c as follows: 
    ! c
    ! c     info(1)   = nfe (number of function evaluations)
    ! c     info(2)   = njve (number of J*v evaluations)
    ! c     info(3)   = nrpre (number of P(inverse)*v evaluations)
    ! c     info(4)   = nli (number of linear iterations)
    ! c     info(5)   = nni (number of nonlinear iterations)
    ! c     info(6)   = nbt (number of backtracks)
    ! c
    ! c ------------------------------------------------------------------------
    ! c
    ! c Optional common blocks: 
    ! c
    ! c These can be used to control printing of diagnostic information by nitsol, 
    ! c to pass information about the nonlinear iterations to jacv or other user 
    ! c subroutines, or to control the default behavior of the nonlinear iterations. 
    ! c
    ! c For controlling printing of diagnostic information: 
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
    ! c     ipunit = printout unit number, e.g., ipunit = 6 => standard output. 
    ! c              NOTE: If ipunit = 0 on input, then it is set to 6 below.
    ! c
    ! c For passing information about the nonlinear iterations to user-supplied 
    ! c subroutines: 
    ! c
    !      include 'nitinfo.h'
    ! c
    ! c If information on the current state of the nonlinear iteration is
    ! c desired in a user-supplied subroutine (for example, deciding 
    ! c whether to update a preconditioner), include this common block
    ! c in the subroutine. The variables are as follows: 
    ! c
    ! c     instep - inexact Newton step number. 
    ! c
    ! c    newstep - set to 0 at the beginning of an inexact Newton step.
    ! c              This may be checked in a user-supplied jacv to decide
    ! c              whether to update the preconditioner.  If you test on
    ! c              newstep .eq. 0 to determine whether to take some 
    ! c              special action at the beginning of a nonlinear iteration, 
    ! c              you must also set newstep to some nonzero value to
    ! c              subsequently avoid taking that action unnecessarily. 
    ! c
    ! c    krystat - status of the Krylov iteration; same as itrmks (see 
    ! c              the nitsol documentation). 
    ! c
    ! c    avrate  - average rate of convergence of the Krylov solver during
    ! c              the previous inexact Newton step.  This may be checked
    ! c              in a user-supplied jacv to decide when to update the
    ! c              preconditioner.
    ! c
    ! c    fcurnrm - ||f(xcur)||. 
    ! c
    ! c        eta - forcing term. 
    ! c
    ! c
    ! c  For controlling the default behavior of the nonlinear iterations:
    ! c
    !      include 'nitparam.h'

    ! c nitparam contains some parameters that control the nonlinear
    ! c iterations.  In some cases, the default values reflect prevailing  
    ! c practice; in other cases, they are chosen to produce good 
    ! c average-case behavior.  To change the default values, include this 
    ! c common block in the main program and set the desired variables 
    ! c according to the following:
    ! c
    ! c    choice1_exp -  parameter used in the update of the forcing term 
    ! c                   eta when ieta = 0 (default).  This is the exponent
    ! c                   for determining the etamin safeguard.  The default
    ! c                   value is choice1_exp = (1+sqrt(5))/2.  A larger
    ! c                   value will allow eta to decrease more rapidly,
    ! c                   while a smaller value will result in a larger 
    ! c                   value for the safeguard. 
    ! c
    ! c    choice2_exp  - parameter used in the update of the forcing term 
    ! c                   eta when ieta = 2.  This is the exponent alpha 
    ! c		    in the expression gamma*(||fcur||/||fprev||)**alpha; 
    ! c		    it is also used to determine the etamin safeguard.  
    ! c		    The default value is 2.0. Valid values are in the 
    ! c		    range (1.0, 2.0].
    ! c
    ! c    choice2_coef - parameter used in the update of the forcing term eta 
    ! c                   when ieta = 2.  This is the coefficient gamma used 
    ! c		    in the expression gamma*(||fcur||/||fprev||)**alpha;
    ! c                   it is also used to determine the etamin safeguard.
    ! c                   The default value is 1.0. Valid values are in the 
    ! c		    range (0.0, 1.0]. 
    ! c
    ! c    eta_cutoff   - parameter used to determine when to disable 
    ! c                   safeguarding the update of the forcing term.  It
    ! c                   only has meaning when ieta .ne. 3.  The default
    ! c                   value is 0.1.  A value of 0.0 will enable 
    ! c		    safeguarding always; a value of 1.0 will disable 
    ! c		    safeguarding always. 
    ! c
    ! c    etamax       - parameter used to provide an upper bound on the 
    ! c		    forcing terms when input(10) .ne. 3. This is 
    ! c		    necessary to ensure convergence of the inexact Newton 
    ! c		    iterates and is imposed whenever eta would otherwise 
    ! c		    be too large. (An overly large eta can result from 
    ! c		    the updating formulas when input(10) .ne. 3 or from 
    ! c                   safeguarding when the previous forcing term has been 
    ! c		    excessively increased during backtracking.) The 
    ! c		    default value of etamax is 1.0 - 1.e-4.  When 
    ! c		    backtracking occurs several times during a nonlinear 
    ! c		    solve the forcing term can remain near etamax for several
    ! c                   nonlinear steps and cause the nonlinear iterations
    ! c                   to nearly stagnate.  In such cases a smaller value of 
    ! c                   etamax may prevent this.  Valid values are in the 
    ! c                   range (0.0, 1.0).
    ! c
    ! c    etafixed     - this is the user-supplied fixed eta when ieta = 3.
    ! c                   The  default value is etafixed = 0.1.  Valid values
    ! c                   are in the range (0.0,1.0).
    ! c
    ! c    thmin        - when backtracking occurs, this is the smallest
    ! c                   reduction factor that will be applied to the current
    ! c                   step in a single backtracking reduction.  The default
    ! c                   value is 0.1.  Valid  values are in the range
    ! c                   [0.0, thmax].
    ! c
    ! c    thmax        - when backtracking occurs, this is the largest
    ! c                   reduction factor that will be applied to the current
    ! c                   step in a single backtracking reduction.  The default
    ! c                   value is 0.5.  Valid values are in the range
    ! c                   [thmin, 1.0).
    ! c
    ! c  The values in this common block are checked once here in nitsol 
    ! c  before the main solution driver is called.  If any parameter has
    ! c  an invalid value, it is silently reset to the default value.
    ! c
    ! c ------------------------------------------------------------------------
    ! c
    ! c Subroutines required by this and all called routines: 
    ! c
    ! c     user supplied: f, jacv, dinpr, dnorm
    ! c
    ! c     nitsol routines: nitbd.f, nitbt, nitdrv, nitgm, nitjv, nitstb, nittfq
    ! c
    ! c     lapack routines: dlaic1, dlamch.f
    ! c
    ! c     blas routines: daxpy, dcopy, dscal, dswap 
    ! c
    ! c This subroutine called by: calling program 
    ! c
    ! c Subroutines called by this subroutine: nitdrv 
    ! c
    ! c ------------------------------------------------------------------------
    ! c
    ! c Remaining declarations:
    ! c

     integer :: lfcur, lfpls, lrwork, lstep, lxpls, nbt, nfe, njve, nli, nni, nrpre 
     if (.not. nitsolInitBool) then
         print *, "Nitsol has not been initialized!"
         stop
      end if

    
    ! c ------------------------------------------------------------------------
    ! c  Initialize some pointers into the rwork array.
    ! c ------------------------------------------------------------------------
     lfcur = 1
     lxpls = lfcur + n
     lfpls = lxpls + n
     lstep = lfpls + n
     lrwork = lstep + n
    ! c ------------------------------------------------------------------------
    ! c Call nitdrv. 
    ! c ------------------------------------------------------------------------
     call nitdrv(n, x, rworkNitsol(lfcur), rworkNitsol(lxpls), rworkNitsol(lfpls), rworkNitsol(lstep), f, jacv, rpar, ipar, abstol, reltol, stptol, &
      iterm, nfe, njve, nrpre, nli, nni, nbt, rworkNitsol(lrwork))
    ! c ------------------------------------------------------------------------
    ! c Set output for return. 
    ! c ------------------------------------------------------------------------
     info(1) = nfe 
     info(2) = njve 
     info(3) = nrpre 
     info(4) = nli 
     info(5) = nni 
     info(6) = nbt 
    ! c ------------------------------------------------------------------------
    ! c All returns made here.
    ! c ------------------------------------------------------------------------
    900  continue
     if ( ipunit .gt. 6 ) close( unit=ipunit )
end subroutine nitsol


end module mod_nitsol