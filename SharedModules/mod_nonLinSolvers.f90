module mod_nonLinSolvers

    use iso_fortran_env, only: int32, real64
    use constants
    use mod_BasicFunctions
    use mod_particle
    use mod_domain
    use mod_potentialSolver
    implicit none

    ! Initialize objects needed
    type(Domain) :: world
    type(Particle), allocatable :: particleList(:)
    type(potentialSolver) :: solver
    integer(int32) :: iterNumPicard, iterNumParticle, iterNumAdaptiveSteps, m_Anderson, amountTimeSplits
    real(real64) :: Beta_k

    
    ! Common blocks for nitsol
    integer iplvl, ipunit
    common /nitprint/ iplvl, ipunit
    double precision choice1_exp, choice2_exp, choice2_coef
    double precision eta_cutoff, etamax
    double precision thmin, thmax, etafixed

    common /nitparam/ choice1_exp, choice2_exp, choice2_coef, eta_cutoff, etamax, thmin, thmax, etafixed

contains

    ! ---------------- Basic Picard ----------------------------
    subroutine solveDivAmperePicard(solver, particleList, world, del_t, maxIter, eps_r)
        ! Solve for divergence of ampere using picard iterations
        class(potentialSolver), intent(in out) :: solver
        type(Particle), intent(in out) :: particleList(:)
        type(Domain), intent(in) :: world
        integer(int32), intent(in) :: maxIter
        real(real64), intent(in) :: del_t, eps_r
        real(real64) :: errorCurrent, errorInitial
        integer(int32) :: i
        call solver%depositJ(particleList, world, del_t)
        errorInitial = solver%getError_tridiag_Ampere(world, del_t)
        do i = 1, maxIter
            call solver%solve_tridiag_Ampere(world, del_t)
            call solver%depositJ(particleList, world, del_t)
            errorCurrent = solver%getError_tridiag_Ampere(world, del_t)
            if (i > 2) then
                if (errorCurrent < eps_r*errorInitial) then
                    call solver%moveParticles(particleList, world, del_t)
                    solver%phi = solver%phi_f
                    exit
                end if
            end if
        end do
        iterNumPicard = i-1

    end subroutine solveDivAmperePicard

    subroutine adaptiveSolveDivAmperePicard(solver, particleList, world, del_t, maxIter, eps_r)
        ! Solve for divergence of ampere's law with picard
        ! cut del_t in two if non-convergence after maxIter, repeat until convergence
        class(potentialSolver), intent(in out) :: solver
        type(Particle), intent(in out) :: particleList(:)
        type(Domain), intent(in) :: world
        integer(int32), intent(in) :: maxIter
        real(real64), intent(in) :: del_t, eps_r
        real(real64) :: remainDel_t, currDel_t, adaptiveJ(NumberXNodes-1)
        call solveDivAmperePicard(solver, particleList, world, del_t, maxIter, eps_r)
        remainDel_t = del_t  
        if (iterNumPicard == maxIter) then
            amountTimeSplits = amountTimeSplits + 1
            do while (iterNumPicard == maxIter)
                print *, "Reducing time step adaptively"
                iterNumAdaptiveSteps = 0
                currDel_t = remainDel_t
                adaptiveJ = 0.0d0
                do while (iterNumPicard == maxIter)
                    currDel_t = currDel_t/2.0d0
                    iterNumAdaptiveSteps = iterNumAdaptiveSteps + 1
                    if (iterNumAdaptiveSteps > 4) then
                        stop "ALREADY REDUCED TIME STEP MORE THAN 3 TIMES, REDUCE INITIAL TIME STEP!!!"
                    end if
                    call solveDivAmperePicard(solver, particleList, world, currDel_t, maxIter, eps_r)   
                end do
                adaptiveJ = adaptiveJ + solver%J * currDel_t/del_t
                remainDel_t = remainDel_t - currDel_t 
                call solveDivAmperePicard(solver, particleList, world, remainDel_t, maxIter, eps_r)
            end do
            adaptiveJ = adaptiveJ + solver%J * remainDel_t/del_t
            solver%J = adaptiveJ
        end if
        

    end subroutine adaptiveSolveDivAmperePicard

    ! ----------- Picard with Anderson Acceleration -------------------------------

    subroutine solveDivAmpereAnderson(solver, particleList, world, del_t, maxIter, eps_r)
        ! Solve for divergence of ampere using picard iterations
        class(potentialSolver), intent(in out) :: solver
        type(Particle), intent(in out) :: particleList(:)
        type(Domain), intent(in) :: world
        integer(int32), intent(in) :: maxIter
        real(real64), intent(in) :: del_t, eps_r
        real(real64) :: initialR, sumPastResiduals, initialNorm
        real(real64) :: Residual_k(NumberXNodes-2, m_Anderson+1), phi_k(NumberXNodes-2, m_Anderson+1), fitMat(NumberXNodes-2, m_Anderson)
        integer(int32) :: lwork
        real(real64), allocatable :: work(:), alpha(:)
        integer(int32) :: i, j, index, m_k, info, ldb
        ldb = MAX(m_Anderson, (NumberXNodes -2))
        lwork= MIN((NumberXNodes -2),m_Anderson) + ldb
        allocate(alpha(ldb), work(lwork))
        phi_k(:,1) = solver%phi(2:NumberXNodes-1)
        call solver%depositJ(particleList, world, del_t)
        initialNorm = SQRT(SUM(solver%phi(2:NumberXNodes-1)**2))
        call solver%solve_tridiag_Ampere(world, del_t)
        phi_k(:,2) = solver%phi_f(2:NumberXNodes-1)
        initialR = SQRT(SUM((solver%phi_f(2:NumberXNodes-1) - phi_k(:,1))**2))
        Residual_k(:,1) = phi_k(:,2) - phi_k(:,1)
        
        do i = 1, maxIter
            index = MODULO(i, m_Anderson+1) + 1
            m_k = MIN(i, m_Anderson)
            ldb = MAX(m_k, (NumberXNodes -2))
            call solver%depositJ(particleList, world, del_t)
            if (SQRT(SUM((solver%phi_f(2:NumberXNodes-1) - phi_k(:,MODULO(i-1, m_Anderson+1) + 1))**2)) < eps_r*(initialR + initialNorm)) then
                call solver%moveParticles(particleList, world, del_t)
                solver%phi = solver%phi_f
                exit
            end if
            call solver%solve_tridiag_Ampere(world, del_t)
            Residual_k(:, index) = solver%phi_f(2:NumberXNodes-1) - phi_k(:,index)
            if (i > m_Anderson) then
                if (m_Anderson > 1) then
                    sumPastResiduals = SUM(Residual_k(:, MODULO(i-1, m_Anderson+1) + 1)**2) &
                    + SUM(Residual_k(:, MODULO(i-2, m_Anderson+1) + 1)**2)
                    sumPastResiduals = sumPastResiduals/2.0d0
                else
                    sumPastResiduals = SUM(Residual_k(:, MODULO(i-1, m_Anderson+1) + 1)**2)
                end if
                if (SUM(Residual_k(:, index)**2) > sumPastResiduals) then
                    iterNumPicard = maxIter
                    goto 75
                end if
            end if
            do j = 0, m_k-1
                fitMat(:,j+1) = Residual_k(:, MODULO(i - m_k + j, m_Anderson+1) + 1) - Residual_k(:, index)
            end do
            alpha(1:NumberXNodes-2) = -Residual_k(1:NumberXNodes-2, index)
            call dgels('N', NumberXNodes-2, m_k, 1, fitMat(:, 1:m_k), NumberXNodes-2, alpha(1:ldb), ldb, work, lwork, info)
            if (info /= 0) then
                print *, "Issue with minimization procedure dgels in Anderson Acceleration!"
                stop
            end if
            alpha(m_k+1) = 1.0d0 - SUM(alpha(1:m_k)) 
            phi_k(:, MODULO(i+1, m_Anderson+1) + 1) = alpha(1) * (Beta_k*Residual_k(:, MODULO(i-m_k, m_Anderson+1) + 1) + phi_k(:, MODULO(i-m_k,m_Anderson+1) + 1))
            do j=1, m_k
                phi_k(:, MODULO(i+1, m_Anderson+1) + 1) = phi_k(:, MODULO(i+1, m_Anderson+1) + 1) + alpha(j + 1) * (Beta_k*Residual_k(:, MODULO(i-m_k + j, m_Anderson+1) + 1) + phi_k(:, MODULO(i-m_k + j, m_Anderson+1) + 1))
            end do
            solver%phi_f(2:NumberXNodes-1) = phi_k(:, MODULO(i+1, m_Anderson+1) + 1)
    
        end do
        iterNumPicard = i-1
        75 continue
        deallocate(alpha, work)
        

    end subroutine solveDivAmpereAnderson

    subroutine adaptiveSolveDivAmpereAnderson(solver, particleList, world, del_t, maxIter, eps_r)
        ! Solve for divergence of ampere's law with picard
        ! cut del_t in two if non-convergence after maxIter, repeat until convergence
        class(potentialSolver), intent(in out) :: solver
        type(Particle), intent(in out) :: particleList(:)
        type(Domain), intent(in) :: world
        integer(int32), intent(in) :: maxIter
        real(real64), intent(in) :: del_t, eps_r
        real(real64) :: remainDel_t, currDel_t, adaptiveJ(NumberXNodes-1)
        call solveDivAmpereAnderson(solver, particleList, world, del_t, maxIter, eps_r)
        remainDel_t = del_t  
        if (iterNumPicard == maxIter) then
            amountTimeSplits = amountTimeSplits + 1
            do while (iterNumPicard == maxIter)
                iterNumAdaptiveSteps = 0
                currDel_t = remainDel_t
                adaptiveJ = 0.0d0
                do while (iterNumPicard == maxIter)
                    currDel_t = currDel_t/2.0d0
                    iterNumAdaptiveSteps = iterNumAdaptiveSteps + 1
                    if (iterNumAdaptiveSteps > 4) then
                        stop "ALREADY REDUCED TIME STEP MORE THAN 3 TIMES, REDUCE INITIAL TIME STEP!!!"
                    end if
                    call solveDivAmpereAnderson(solver, particleList, world, currDel_t, maxIter, eps_r)   
                end do
                adaptiveJ = adaptiveJ + solver%J * currDel_t/del_t
                remainDel_t = remainDel_t - currDel_t 
                call solveDivAmpereAnderson(solver, particleList, world, remainDel_t, maxIter, eps_r)
            end do
            adaptiveJ = adaptiveJ + solver%J * remainDel_t/del_t
            solver%J = adaptiveJ
        end if
    end subroutine adaptiveSolveDivAmpereAnderson

    ! -------------------- JFNK functions ---------------------------------------
    subroutine funcNitsol(n, xcur, fcur, rpar, ipar, itrmf)
        ! Use solver and whatnot as global inputs, I'm sure as hell not combining all data into one rpar and then creating new routines on those!
        integer(int32), intent(in) :: n
        integer(int32), intent(in out) :: itrmf, ipar(*)
        real(real64), intent(in) :: xcur(n)
        real(real64), intent(in) :: rpar
        real(real64), intent(in out) :: fcur(n)
        real(real64) :: d(n)
        solver%phi_f(2:NumberXNodes-1) = xcur
        call solver%depositJ(particleList, world, rpar)
        d = (-solver%J(2:) + solver%J(1:NumberXNodes-2)) * rpar / eps_0 &
        + arrayDiff(solver%phi(1:NumberXNodes-1), NumberXNodes-1)*2.0d0/(world%dx_dl(1:NumberXNodes-2) + world%dx_dl(2:NumberXNodes-1)) &
        - arrayDiff(solver%phi(2:), NumberXNodes-1)*2.0d0/(world%dx_dl(3:NumberXNodes) + world%dx_dl(2:NumberXNodes-1))
        d(1) = d(1) + solver%phi(1) * solver%coeff_left
        d(n) = d(n) + solver%phi(NumberXNodes) * solver%coeff_right
        d = triMul(n, solver%a_tri, solver%c_tri, solver%b_tri, xcur) - d
        call solve_tridiag(n, solver%a_tri, solver%c_tri, solver%b_tri, d, fcur)
        itrmf = 0

    end subroutine funcNitsol

    subroutine jacNitsol(n, xcur, fcur, ijob, v, z, rpar, ipar, itrmjv) 
        ! If analytical jacobian matrix-vector product or preconditioner needed
        integer, intent(in) :: ijob
        integer, intent(in out) :: itrmjv
        integer, intent(in) :: n

        integer, intent(in out) :: ipar(*)

        real(real64), intent(in) :: fcur(n)
        real(real64), intent(in) :: rpar
        real(real64), intent(in) :: v(n)
        real(real64), intent(in) ::  xcur(n)
        real(real64), intent(in out) :: z(n)
        if (ijob == 0) then
            print *, "in analytical jacobian if block"
        else if (ijob == 1) then
            call solve_tridiag(n, solver%a_tri, solver%c_tri, solver%b_tri, v, z)
        end if
        itrmjv = 0
    end subroutine jacNitsol

    subroutine solveJFNK(del_t, maxIter, eps_r)
        integer(int32), intent(in) :: maxIter
        real(real64), intent(in) :: del_t, eps_r
        real(real64) :: initialNorm, rpar(2), KE_i, KE_f, PE_i, PE_f
        real(real64), allocatable :: rwork(:), fcur(:), xcur(:)
        integer(int32) :: input(10), info(6), kdmax, iterm, ipar(2), itrmf, j
        double precision, external :: ddot
        double precision, external :: dnrm2

        ! Set Nitsol parameters
        iplvl = 1 ! maximum print
        iterm = 0
        kdmax = 20 ! maximum krylov subspace dimension
        input = 0
        input(1) = 100 ! maximum iterations
        input(2) = 0 !ijacv
        input(4) = kdmax ! maximum krylov subspace dimension
        input(5) = 0 !ipre
        input(9) = -1
        input(10) = 2 ! eta with gamma and alpha
        etamax = 0.8d0 ! eta max
        choice2_exp = 1.5d0 ! alpha
        choice2_coef = 0.9d0 ! gamma
        allocate(fcur(NumberXNodes-2), rwork((NumberXNodes-2)*(kdmax+5)+kdmax*(kdmax+3)), xcur(NumberXNodes-2))
        xcur = solver%phi_f(2:NumberXNodes-1)
        call funcNitsol(NumberXNodes-2, xcur, fcur, del_t, ipar, itrmf)
        initialNorm = dnrm2(NumberXNodes-2, fcur, 1)
        print *, "initialNorm is:", initialNorm
        solver%particleEnergyLoss = 0.0d0
        PE_i = solver%getTotalPE(world, .false.)
        KE_i = 0.0d0
        do j=1, numberChargedParticles
            KE_i = KE_i + particleList(j)%getTotalKE()
        end do
        call nitsol(NumberXNodes-2, xcur, funcNitsol, jacNitsol, eps_r*initialNorm, eps_r,input, info, rwork, del_t, ipar, iterm, ddot, dnrm2)
        call solver%moveParticles(particleList, world, del_t)
        solver%phi = solver%phi_f
        KE_f = solver%particleEnergyLoss
        do j=1, numberChargedParticles
            KE_f = KE_f + particleList(j)%getTotalKE()
        end do
        PE_f = solver%getTotalPE(world, .false.)
        solver%energyError = ABS((KE_i + PE_i - KE_f - PE_f)/(KE_i + PE_i))
        print *, "energy error is:", solver%energyError
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
        deallocate(fcur, rwork, xcur)
    end subroutine solveJFNK


end module mod_nonLinSolvers