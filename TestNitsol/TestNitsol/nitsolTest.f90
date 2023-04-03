program nitsolTest
    use iso_fortran_env, only: int32, real64
    implicit none

    integer iplvl, ipunit
    common /nitprint/ iplvl, ipunit
    double precision choice1_exp, choice2_exp, choice2_coef
    double precision eta_cutoff, etamax
    double precision thmin, thmax, etafixed

    common /nitparam/ choice1_exp, choice2_exp, choice2_coef, eta_cutoff, etamax, thmin, thmax, etafixed
    real(real64), parameter :: c = 299792458.0d0 ! m/s
    real(real64), parameter :: eps_0 = 8.8541878128d-12 !F/m
    real(real64), parameter :: m_p = 1.67262192369d-27 ! kg
    real(real64), parameter :: m_e = 9.1093837015d-31 ! kg
    real(real64), parameter :: k_B = 1.380649d-23 ! m^2 kg s^-2 K^-1
    real(real64), parameter :: e = 1.602176634d-19 ! C
    real(real64), parameter :: mu_0 = 1.25663706212d-6 ! m kg s^-2 A^-2
    real(real64), parameter :: pi = 4.0d0*atan(1.0d0) ! pi from atan
    integer(int32) :: num_grid_nodes = 64, i,j, ipar(3), itrmf, input(10), info(6), kdmax, iterm, infoLpack
    real(real64) :: phi(64), grid(64), rho(64), L_domain = 0.1d0, b(62), delx, phi_Pic(62), A_pic(62,62), fcur(62), diag(62), diagLower(61), diagUpper(61), rpar(62), ftol, initialNormR, stptol
    real(real64), allocatable :: rwork(:)
    double precision, external :: ddot
    double precision, external :: dnrm2
    iplvl = 4 ! maximum print
    iterm = 0
    kdmax = 60 ! maximum krylov subspace dimension
    input = 0
    input(1) = 50 ! maximum iterations
    input(2) = 0
    input(4) = kdmax
    input(5) = 1
    input(10) = 2
    etamax = 0.8d0
    choice2_exp = 1.5d0
    choice2_coef = 0.9d0
    allocate(rwork(62*(kdmax+5)+kdmax*(kdmax+3)))
    ftol = 1.0d-6
    stptol = 1.0d-8
    rho = e * 1.0d12
    delx = L_domain/63.0d0
    phi_Pic = 0.0d0
    b = -rho(2:63) * delx**2 /eps_0
    A_pic = 0.0d0
    diagLower = 1.0d0
    diagUpper = 1.0d0
    diag = -2.0d0
    rpar = solve_tridiag(62, diagLower, diagUpper, diag, b)
    do i=1, 62
        do j=1, 62
            if (ABS(j-i) == 1) then
                A_pic(i,j) = 1.0d0
            end if
            if (i == j) then
                A_pic(i,j) = -2.0d0
            end if
        end do
    end do
    fcur = b
    call func(62, phi_Pic, fcur, rpar, ipar, itrmf)
    initialNormR = dnrm2(62, fcur, 1)
    print *, "Initial norm is:", initialNormR
    call nitsol(62, phi_Pic, func, jacv, ftol*initialNormR, stptol,input, info, rwork, rpar, ipar, iterm, ddot, dnrm2)
    print *, ABS((phi_Pic - rpar)/rpar) * 100
    write(6,*) 
    write(6,880) iterm
    write(6,900) info(1)
    write(6,910) info(2)
    write(6,920) info(3)
    write(6,930) info(4)
    write(6,940) info(5)
    write(6,950) info(6)
    880  format(' Termination flag iterm:       ', i9)
    890  format(' Final f-norm:                 ', t36, 1pe9.3)
    900  format(' No. function evaluations:     ', i9)
    910  format(' No. J*v evaluations:          ', i9) 
    920  format(' No. P(inverse)*v evaluations: ', i9)
    930  format(' No. linear iterations:        ', i9)
    940  format(' No. nonlinear iterations:     ', i9)
    950  format(' No. backtracks:               ', i9)
contains

    function solve_tridiag(n, diagLower, diagUpper, diag, b) result(x)
        ! Tridiagonal (Thomas algorithm) solver for Ampere
        integer(int32), intent(in) :: n
        real(real64), intent(in) :: diagLower(n-1), diagUpper(n-1), diag(n), b(n)
        integer(int32) :: i !n size dependent on how many points are inside (not boundary), so how many in rhs equation
        real(real64) :: m, cp(n-1),dp(n), x(n)

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
        x(n) = dp(n)
        do i = n-1, 1, -1
            x(i) = dp(i)-cp(i)*x(i+1)
        end do
    end function solve_tridiag

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

    subroutine func(n, xcur, fcur, rpar, ipar, itrmf)
        integer(int32), intent(in) :: n
        integer(int32), intent(in out) :: itrmf, ipar(*)
        real(real64), intent(in) :: xcur(n)
        real(real64), intent(in) :: rpar(n)
        real(real64), intent(in out) :: fcur(n)
        fcur = triMul(n, diagLower, diagUpper, diag, xcur) - b
        itrmf = 0

    end subroutine

    subroutine jacv(n, xcur, fcur, ijob, v, z, rpar, ipar, itrmjv) 
        integer, intent(in) :: ijob
        integer, intent(in out) :: itrmjv
        integer, intent(in) :: n

        integer, intent(in out) :: ipar(*)

        real(real64), intent(in) :: fcur(n)
        real(real64), intent(in) :: rpar(n)
        real(real64), intent(in) :: v(n)
        real(real64), intent(in) ::  xcur(n)
        real(real64), intent(in out) :: z(n)
        if (ijob == 0) then
            z = triMul(n, diagLower, diagUpper, diag, v)
        else if (ijob == 1) then
            z = solve_tridiag(62, diagLower, diagUpper, diag, v)
        end if
        itrmjv = 0
    end subroutine jacv

end program nitsolTest