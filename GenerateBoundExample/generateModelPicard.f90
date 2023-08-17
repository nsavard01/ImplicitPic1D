program generateModel

    use iso_fortran_env, only: int32, real64
    use constants
    use mod_BasicFunctions
    implicit none
    real(real64) :: L, tau, fracDebye, lambDebye, dt, dah, M_m, n_0, B, A, rpar, InitialNorm, currNorm
    integer(int32) :: maxNodes, i, j, inputJFNK(10), ipar, itrmf, info(6), iterm,k
    real(real64), allocatable :: phi(:), t(:), s(:), a_tri(:), I_G(:), eta(:), temp(:), I_int(:), truncFact(:), b_tri(:), c_tri(:)

    eps_r = 1.d-8
    maxIter = 1


    fracDebye = 0.1d0
    T_e = 5.0d0
    T_i = 1.0d0
    n_ave = 5.d14
    L = 0.05d0
    tau = T_e/T_i
    M_m = m_p/m_e
    lambDebye = getDebyeLength(T_e, n_ave)
    maxNodes = int(L/lambDebye)
    NumberXNodes = maxNodes
    dt = 1.0d0/real(NumberXNodes-1)
    dah = LOG(fracDebye * lambDebye/L)/LOG(1.0d0/real(NumberXNodes-1))
    allocate(phi(NumberXNodes), t(NumberXNodes), s(NumberXNodes), I_G(NumberXNodes), eta(NumberXNodes), temp(NumberXNodes), truncFact(NumberXNodes), a_tri(NumberXNodes-1), &
    b_tri(NumberXNodes), c_tri(NumberXNodes-1))
    do i = 1, NumberXNodes
       t(i) = (i-1) * dt 
    end do
    s = 1.0d0 - (1.0d0 - t)**dah
    b_tri(1) = 1.0d0
    c_tri(1) = 0.0d0
    do i = 2, NumberXNodes-1
        c_tri(i) = ((1.0d0 - t(i))**(2.0d0 - 2.0d0 * dah) / dah**2 / dt**2 + (1.0d0 - t(i))**(1.0d0 - 2.0d0 * dah) * (dah - 1.0d0) / dah**2 / dt)
        a_tri(i-1) = ((1.0d0 - t(i))**(2.0d0 - 2.0d0 * dah) / dah**2 / dt**2 - (1.0d0 - t(i))**(1.0d0 - 2.0d0 * dah) * (dah - 1.0d0) / dah**2 / dt)
        b_tri(i) = -2.0d0 * (1.0d0 - t(i))**(2.0d0 - 2.0d0 * dah) / dah**2 / dt**2 
    end do
    b_tri(NumberXNodes) = 1.0d0
    a_tri(NumberXNodes-1) = 0.0d0
    eta = -1.0d0
    eta(NumberXNodes) = 0.0d0

    B = 0.5d0 * SQRT(M_m * tau)

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
    temp(1) = -LOG(I_G(1) * B)/(tau + 1.0d0)
    temp(2:NumberXNodes-1) = EXP(eta(1)) * (-EXP(-eta(2:NumberXNodes-1)) * truncFact(2:NumberXNodes-1) + B * EXP(eta(2:NumberXNodes-1) * tau) * I_G(2:NumberXNodes-1)) / A
    temp(NumberXNodes) = 0.0d0
    phi = eta
    call solve_tridiag(NumberXNodes, a_tri, c_tri, b_tri, temp, eta)
    InitialNorm = SQRT(SUM((eta - phi)**2)) 
    do k = 1, maxIter
        temp = EXP(-eta) * dah * (1.0d0-t)**(dah-1)
        print *, temp
        stop
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
        temp(1) = -LOG(I_G(1) * B)/(tau + 1.0d0)
        temp(2:NumberXNodes-1) = EXP(eta(1)) * (-EXP(-eta(2:NumberXNodes-1)) * truncFact(2:NumberXNodes-1) + B * EXP(eta(2:NumberXNodes-1) * tau) * I_G(2:NumberXNodes-1)) / A
        temp(NumberXNodes) = 0.0d0
        print *, temp
        phi = eta
        call solve_tridiag(NumberXNodes, a_tri, c_tri, b_tri, temp, eta)
        currNorm = SQRT(SUM((eta - phi)**2)) 
        print *, currNorm
        if (currNorm < eps_r * InitialNorm) then
            print *, 'Iteration has converged in', k, 'iterations'
            exit
        end if
    end do
    if (k == maxIter) then
        print *, '----------'
        print *, 'Did not converge!'
    end if
 
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