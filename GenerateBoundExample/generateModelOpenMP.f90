program generateModel

    use iso_fortran_env, only: int32, real64
    use constants
    use mod_BasicFunctions
    use omp_lib
    implicit none
    real(real64) :: L, tau, fracDebye, lambDebye, dt, dah, M_m, n_0, B, A, A_coeff, InitialNorm, currNorm, d
    integer(int32) :: maxNodes, i, j, inputJFNK(10), ipar, itrmf, info(6), iterm,k, u, numThreads, threadDiv, iThread
    integer(int32), allocatable :: threadIndx(:)
    real(real64), allocatable :: phi(:), t(:), s(:), I_G(:), eta(:), temp(:), truncFact(:), fcur(:), rworkSolver(:), rpar(:, :)

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

    numThreads = 32
    call omp_set_num_threads(numThreads)
    fracDebye = 0.05d0
    T_e = 5.0d0
    T_i = 1.0d0
    n_ave = 5.0d14
    L = 0.05d0
    tau = T_e/T_i
    M_m = m_p/m_e
    lambDebye = getDebyeLength(T_e, n_ave)
    maxNodes = int(L/lambDebye)
    threadDiv = INT(maxNodes/numThreads)
    NumberXNodes = (threadDiv+1) * numThreads
    dt = 1.0d0/real(NumberXNodes-1)
    dah = LOG(fracDebye * lambDebye/L)/LOG(1.0d0/real(NumberXNodes-1))
    allocate(phi(NumberXNodes), t(NumberXNodes), s(NumberXNodes), I_G(NumberXNodes), eta(NumberXNodes), temp(NumberXNodes), truncFact(NumberXNodes), fcur(NumberXNodes), &
    rworkSolver((NumberXNodes)*(inputJFNK(4)+5)+inputJFNK(4)*(inputJFNK(4)+3)), threadIndx(numThreads), rpar(NumberXNodes, 4))
    do i = 1, NumberXNodes
       t(i) = (i-1) * dt 
    end do
    s = 1.0d0 - (1.0d0 - t)**dah
    rpar(:, 3) = t
    eta = -1.0d0
    eta(NumberXNodes) = 0.0d0
    threadIndx = (/(i, i = 1, NumberXNodes - threadDiv, threadDiv + 1)/)
    B = 0.5d0 * SQRT(M_m * tau)
    A_coeff = eps_0 * T_e / e / L**2

    do u = 1, 2
        print *, 'u is:', u
        print *, 'number of nodes is:', NumberXNodes
        ! temp = EXP(-eta) * dah * (1.0d0-t)**(dah-1)
        ! n_0 = n_ave/EXP(eta(1))/trapZ(temp, dt, NumberXNodes)
        call getA(A, eta, temp)
        call getI_G(I_G, eta)
        call getTruncFact(truncFact, eta)
        rpar(:, 1) = I_G
        rpar(:, 2) = truncFact
        rpar(1, 4) = A
        rpar(2, 4) = B
        call funcNitsol(NumberXNodes, eta, temp, rpar, ipar, itrmf)
        InitialNorm = dnrm2(NumberXNodes, temp, 1)
        call nitsol(NumberXNodes, eta, funcNitsol, jacNitsol, eps_r*initialNorm, eps_r,inputJFNK, info, rworkSolver, rpar, ipar, iterm, ddot, dnrm2)
        print *, eta
        stop
        do k = 1, maxIter
            call getA(A, eta, temp)
            call getI_G(I_G, eta)
            call getTruncFact(truncFact, eta)
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
        threadDiv = 2*threadDiv + 1
        threadIndx = (/(i, i = 1, NumberXNodes - threadDiv, threadDiv + 1)/)
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
                d = (temp(i) - s(j-1))/(s(j) - s(j-1))
                phi(i) = eta(j-1) * (1.0d0 - d) + d * eta(j)
            else
                d = (temp(i) - s(j-1))/(s(j) - s(j-1))
                phi(i) = eta(j-1) * (1.0d0 - d) + d * eta(j)
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
        real(real64), intent(in out) :: rpar(n, 4)
        real(real64), intent(in out) :: fcur(n)
        real(real64) :: laplace, firstDer, bracket
        integer(int32) :: i, beg, end
        fcur(1) = EXP(-xcur(1)) - rpar(2,4)* EXP(xcur(1)*tau) * rpar(1, 1)
        !$OMP parallel private(iThread, i, beg, end, firstDer ,laplace, bracket)
        iThread = omp_get_thread_num() + 1
        if (iThread == 1) then
            beg = 2
            end = 1 + threadDiv
        else if (iThread == numThreads) then
            beg = threadIndx(iThread)
            end = threadIndx(iThread) + threadDiv - 1
        else
            beg = threadIndx(iThread)
            end = threadIndx(iThread) + threadDiv
        end if
        do i = beg, end
            laplace = xcur(i-1) - 2.0d0 * xcur(i) + xcur(i+1)
            firstDer = xcur(i+1) - xcur(i-1)
            bracket = laplace * (1.0d0 - rpar(i, 3))**(2.0d0 - 2.0d0 * dah) / dah**2 / dt**2 + firstDer * (1.0d0 - rpar(i,3))**(1.0d0 - 2.0d0 * dah) * (dah - 1.0d0) / dah**2 / dt
            fcur(i) = rpar(1,4) * bracket * EXP(-xcur(1)) + EXP(-xcur(i))*rpar(i,2) - rpar(2,4) * EXP(xcur(i) * tau) * rpar(i,1)
        end do
        !$OMP end parallel
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

    subroutine getA(A, eta, temp)
        real(real64), intent(in out) :: A
        real(real64), intent(in) :: eta(NumberXNodes)
        real(real64), intent(in out) :: temp(NumberXNodes)
        temp = EXP(-eta) * dah * (1.0d0-t)**(dah-1)
        !$OMP parallel private(iThread)
            iThread = omp_get_thread_num() + 1
            temp(threadIndx(iThread):threadIndx(iThread) + threadDiv) = EXP(-eta(threadIndx(iThread):threadIndx(iThread) + threadDiv)) * dah * (1.0d0-t(threadIndx(iThread):threadIndx(iThread) + threadDiv))**(dah-1)
        !$OMP end parallel
        n_0 = n_ave/EXP(eta(1))/trapZ(temp, dt, NumberXNodes)
        A = A_coeff/n_0

    end subroutine getA

    subroutine getI_G(I_G, eta)
        real(real64), intent(in out) :: I_G(NumberXNodes)
        real(real64), intent(in) :: eta(NumberXNodes)
        real(real64) :: array(NumberXNodes)
        integer(int32) :: i, j
        !$OMP parallel private(iThread, i, j, array)
        iThread = omp_get_thread_num() + 1
        do i = threadIndx(iThread), threadIndx(iThread) + threadDiv
            do j = 1, NumberXNodes
                if (eta(i) > eta(j)) then
                    array(j) = (1.0d0 - t(j))**(dah-1.0d0) * EXP(-eta(j) * tau) * ERFC(SQRT((eta(i) - eta(j))*tau))
                else
                    array(j) = (1.0d0 - t(j))**(dah-1.0d0) * EXP(-eta(j) * tau)
                end if
            end do
            I_G(i) = trapZ(array, dt, NumberXNodes) * dah
        end do
        !$OMP end parallel
    end subroutine getI_G

    subroutine getTruncFact(truncFact, eta)
        real(real64), intent(in out) :: truncFact(NumberXNodes)
        real(real64), intent(in) :: eta(NumberXNodes)
        !$OMP parallel private(iThread)
        iThread = omp_get_thread_num() + 1
        truncFact(threadIndx(iThread):threadIndx(iThread) + threadDiv) = (1.0d0 - 0.5d0 * ERFC(SQRT(-eta(threadIndx(iThread):threadIndx(iThread) + threadDiv))))/(1.0d0 - 0.5d0 * ERFC(SQRT(-eta(1))))
        !$OMP end parallel
    end subroutine getTruncFact

    function trapZ(y, dx, n) result(res)
        integer(int32), intent(in) :: n
        real(real64), intent(in) :: dx, y(n)
        real(real64) :: res
        res = 0.5d0 * SUM(y(2:n) + y(1:n-1)) * dx
    end function trapZ

end program generateModel