program generateModel

    use iso_fortran_env, only: int32, real64
    use constants
    use mod_BasicFunctions
    implicit none
    real(real64) :: L, tau, fracDebye, lambDebye, dt, dah, M_m, n_0, B, A, rpar, InitialNorm, currNorm
    integer(int32) :: maxNodes, i, j, inputJFNK(10), ipar, itrmf, info(6), iterm,k, u
    real(real64), allocatable :: phi(:), t(:), s(:), I_G(:), eta(:), temp(:), I_int(:), truncFact(:), fcur(:), rworkSolver(:)

    ! Common blocks for nitsol
    integer iplvl, ipunit
    common /nitprint/ iplvl, ipunit
    double precision choice1_exp, choice2_exp, choice2_coef
    double precision eta_cutoff, etamax
    double precision thmin, thmax, etafixed
    double precision, external :: ddot
    double precision, external :: dnrm2

    common /nitparam/ choice1_exp, choice2_exp, choice2_coef, eta_cutoff, etamax, thmin, thmax, etafixed

    integer instep, newstep, krystat
    double precision avrate, fcurnrm
    common /nitinfo/ avrate, fcurnrm, instep, newstep, krystat
    iplvl = 0
    inputJFNK = 0
    inputJFNK(1) = 100 ! maximum iterations
    inputJFNK(2) = 0 !ijacv
    inputJFNK(3) = 0 ! krylov solver
    inputJFNK(4) = 20 ! maximum krylov subspace dimension
    inputJFNK(5) = 0 !ipre
    inputJFNK(9) = -1
    inputJFNK(6) = 20
    inputJFNK(10) = 2 ! eta with gamma and alpha
    etamax = 0.8d0 ! eta max
    choice2_exp = 1.5d0 ! alpha
    choice2_coef = 0.9d0 ! gamma
    eps_r = 1.d-8
    maxIter = 10000


    fracDebye = 0.05d0
    T_e = 5.0d0
    T_i = 1.0d0
    n_ave = 2.5d16
    L = 0.05d0
    tau = T_e/T_i
    M_m = m_p/m_e
    lambDebye = getDebyeLength(T_e, n_ave)
    maxNodes = int(L/lambDebye)
    NumberXNodes = maxNodes
    dt = 1.0d0/real(NumberXNodes-1)
    dah = LOG(fracDebye * lambDebye/L)/LOG(1.0d0/real(NumberXNodes-1))
    allocate(phi(NumberXNodes), t(NumberXNodes), s(NumberXNodes), I_G(NumberXNodes), eta(NumberXNodes), temp(NumberXNodes), truncFact(NumberXNodes), fcur(NumberXNodes), &
    rworkSolver((NumberXNodes)*(inputJFNK(4)+5)+inputJFNK(4)*(inputJFNK(4)+3)))
    do i = 1, NumberXNodes
       t(i) = (i-1) * dt 
    end do
    s = 1.0d0 - (1.0d0 - t)**dah
    eta = -1.0d0
    eta(NumberXNodes) = 0.0d0
    B = 0.5d0 * SQRT(M_m * tau)

    do u = 1, 4
        print *, 'u is:', u
        print *, 'number of nodes is:', NumberXNodes
        temp = EXP(-eta) * dah * (1.0d0-t)**(dah-1)
        n_0 = n_ave/EXP(eta(1))/trapZ(temp, dt, NumberXNodes)
        A = eps_0 * T_e / n_0 / e / L**2
        do i = 1, NumberXNodes
            do j = 1, NumberXNodes
                if (eta(i) > eta(j)) then
                    temp(j) = (1.0d0 - t(j))**(dah-1.0d0) * EXP(-eta(j) * tau) * ERFC(SQRT((eta(i) - eta(j))*tau))
                else
                    temp(j) = (1.0d0 - t(j))**(dah-1.0d0) * EXP(-eta(j) * tau)
                end if
            end do
            I_G(i) = trapZ(temp, dt, NumberXNodes)
        end do
        I_G = I_G * dah
        truncFact = (1.0d0 - 0.5d0 * ERFC(SQRT(-eta)))/(1.0d0 - 0.5d0 * ERFC(SQRT(-eta(1))))
        call funcNitsol(NumberXNodes, eta, temp, rpar, ipar, itrmf) 
        InitialNorm = dnrm2(NumberXNodes, temp, 1)
        call nitsol(NumberXNodes, eta, funcNitsol, jacNitsol, eps_r*initialNorm, eps_r,inputJFNK, info, rworkSolver, rpar, ipar, iterm, ddot, dnrm2)
        do k = 1, maxIter
            temp = EXP(-eta) * dah * (1.0d0-t)**(dah-1)
            n_0 = n_ave/EXP(eta(1))/trapZ(temp, dt, NumberXNodes)
            A = eps_0 * T_e / n_0 / e / L**2
            do i = 1, NumberXNodes
                do j = 1, NumberXNodes
                    if (eta(i) > eta(j)) then
                        temp(j) = (1.0d0 - t(j))**(dah-1.0d0) * EXP(-eta(j) * tau) * ERFC(SQRT((eta(i) - eta(j))*tau))
                    else
                        temp(j) = (1.0d0 - t(j))**(dah-1.0d0) * EXP(-eta(j) * tau)
                    end if
                end do
                I_G(i) = trapZ(temp, dt, NumberXNodes)
            end do
            I_G = I_G * dah
            truncFact = (1.0d0 - 0.5d0 * ERFC(SQRT(-eta)))/(1.0d0 - 0.5d0 * ERFC(SQRT(-eta(1))))
            temp = eta
            call nitsol(NumberXNodes, eta, funcNitsol, jacNitsol, eps_r*InitialNorm, eps_r,inputJFNK, info, rworkSolver, rpar, ipar, iterm, ddot, dnrm2)
            currNorm = SQRT(SUM((eta - temp)**2) / real(NumberXNodes))
            if (currNorm < eps_r) then
                print *, 'Iteration has converged in', k, 'iterations'
                exit
            end if
        end do
        if (k == maxIter) then
            print *, '----------'
            print *, 'Did not converge!'
        end if
        deallocate(phi, t, temp)
        NumberXNodes = 2*NumberXNodes
        allocate(phi(NumberXNodes), t(NumberXNodes), temp(NumberXNodes))
        dt = 1.0d0/real(NumberXNodes-1)
        dah = LOG(fracDebye * lambDebye/L)/LOG(1.0d0/real(NumberXNodes-1))
        do i = 1, NumberXNodes
            t(i) = (i-1) * dt 
        end do
        temp = 1.0d0 - (1.0d0 - t)**dah
        j = 2
        do i = 1, NumberXNodes
            if (temp(i) == s(j)) then
                phi(i) = eta(j)
            else if (temp(i) > s(j)) then
                j = j + 1
                rpar = (temp(i) - s(j-1))/(s(j) - s(j-1))
                phi(i) = eta(j-1) * (1.0d0 - rpar) + rpar * eta(j)
            else
                rpar = (temp(i) - s(j-1))/(s(j) - s(j-1))
                phi(i) = eta(j-1) * (1.0d0 - rpar) + rpar * eta(j)
            end if
        end do
        deallocate(s, I_G, eta, truncFact, fcur, rworkSolver)
        allocate(s(NumberXNodes), I_G(NumberXNodes), eta(NumberXNodes), truncFact(NumberXNodes), fcur(NumberXNodes), rworkSolver((NumberXNodes)*(inputJFNK(4)+5)+inputJFNK(4)*(inputJFNK(4)+3)))
        s = temp
        eta = phi
    end do
    print *, eta(1:10) * -T_e
    
 
contains

    subroutine funcNitsol(n, xcur, fcur, rpar, ipar, itrmf)
        ! Use solver and whatnot as global inputs, I'm sure as hell not combining all data into one rpar and then creating new routines on those!
        integer(int32), intent(in) :: n
        integer(int32), intent(in out) :: itrmf, ipar
        real(real64), intent(in) :: xcur(n)
        real(real64), intent(in out) :: rpar
        real(real64), intent(in out) :: fcur(n)
        real(real64) :: laplace, firstDer, bracket
        integer(int32) :: i
        fcur(1) = EXP(-xcur(1)) - B* EXP(xcur(1)*tau) * I_G(1)
        do i = 2, NumberXNodes-1
            laplace = xcur(i-1) - 2.0d0 * xcur(i) + xcur(i+1)
            firstDer = xcur(i+1) - xcur(i-1)
            bracket = laplace * (1.0d0 - t(i))**(2.0d0 - 2.0d0 * dah) / dah**2 / dt**2 + firstDer * (1.0d0 - t(i))**(1.0d0 - 2.0d0 * dah) * (dah - 1.0d0) / dah**2 / dt
            fcur(i) = A * bracket * EXP(-xcur(1)) + EXP(-xcur(i))*truncFact(i) - B * EXP(xcur(i) * tau) * I_G(i)
        end do
        fcur(NumberXNodes) = xcur(NumberXNodes) - 0.0d0
        itrmf = 0

    end subroutine funcNitsol

    subroutine jacNitsol(n, xcur, fcur, ijob, v, z, rpar, ipar, itrmjv) 
        ! If analytical jacobian matrix-vector product or preconditioner needed
        integer, intent(in) :: ijob
        integer, intent(in out) :: itrmjv
        integer, intent(in) :: n

        integer, intent(in out) :: ipar

        real(real64), intent(in) :: fcur(n)
        real(real64), intent(in out) :: rpar
        real(real64), intent(in) :: v(n)
        real(real64), intent(in) ::  xcur(n)
        real(real64), intent(in out) :: z(n)
        if (ijob == 0) then
            print *, "in analytical jacobian if block"
        else if (ijob == 1) then
            print *, 'in right preconditioning block'
        end if
        itrmjv = 0
    end subroutine jacNitsol

end program generateModel