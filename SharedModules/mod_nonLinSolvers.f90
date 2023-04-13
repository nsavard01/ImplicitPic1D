module mod_nonLinSolvers

    use iso_fortran_env, only: int32, real64
    use constants
    use mod_BasicFunctions
    use mod_particle
    use mod_domain
    use mod_potentialSolver
    implicit none

    ! Initialize objects needed
    type(Domain) :: globalWorld
    type(Particle), allocatable :: globalParticleList(:)
    type(potentialSolver) :: globalSolver
    integer(int32) :: iterNumPicard, iterNumParticle, iterNumAdaptiveSteps, m_Anderson, amountTimeSplits, solverType
    real(real64) :: Beta_k

    !allocatable arrays for JFNK or Anderson
    integer(int32), allocatable :: inputJFNK(:)

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

contains

    subroutine initializeSolver(eps_r, solverType, m_Anderson, Beta_k, maxIter)
        real(real64), intent(in out) :: eps_r, Beta_k
        integer(int32), intent(in out) :: m_Anderson, solverType, maxIter
        integer(int32) :: io
        print *, "Reading non-linear solver inputs:"
        open(10,file='../InputData/SolverState.inp', IOSTAT=io)
        read(10, *, IOSTAT = io) eps_r
        read(10, *, IOSTAT = io) solverType
        read(10, *, IOSTAT = io) m_Anderson
        read(10, *, IOSTAT = io) Beta_k
        read(10, *, IOSTAT = io) maxIter
        close(10)
        print *, "Relative error:", eps_r
        select case (solverType)
        case(0)
            print *, "Solver type is: Anderson Acceleration"
        case(1)
            print *, "Solver type is: JFNK"
        case default
            print *, "Solver type incorrect!"
            stop
        end select
        print *, "Maximum iteration number:", maxIter
        print *, "Anderson number m is:", m_Anderson
        print *, "Relaxation parameter is:", Beta_k
        print *, ""
        SELECT CASE (solverType)
        CASE(1)
            allocate(inputJFNK(10))
            iplvl = 0 ! print level
            inputJFNK = 0
            inputJFNK(1) = maxIter ! maximum iterations
            inputJFNK(2) = 0 !ijacv
            inputJFNK(3) = 0 ! krylov solver
            inputJFNK(4) = m_Anderson ! maximum krylov subspace dimension
            inputJFNK(5) = 0 !ipre
            inputJFNK(9) = -1
            inputJFNK(6) = m_Anderson
            inputJFNK(10) = 2 ! eta with gamma and alpha
            etamax = 0.8d0 ! eta max
            choice2_exp = 1.5d0 ! alpha
            choice2_coef = 0.9d0 ! gamma
        END SELECT

    end subroutine initializeSolver

    ! ---------------- Basic Picard ----------------------------
    subroutine solveDivAmperePicard(solver, particleList, world, del_t, maxIter, eps_r)
        ! Solve for divergence of ampere using picard iterations
        type(potentialSolver), intent(in out) :: solver
        type(Particle), intent(in out) :: particleList(:)
        type(Domain), intent(in) :: world
        integer(int32), intent(in) :: maxIter
        real(real64), intent(in) :: del_t, eps_r
        real(real64) :: errorCurrent, errorInitial
        integer(int32) :: i
        call solver%depositJ(particleList, world, del_t)
        errorInitial = SQRT(SUM(solver%getError_tridiag_Ampere(world, del_t)**2))
        do i = 1, maxIter
            call solver%solve_tridiag_Ampere(world, del_t)
            call solver%depositJ(particleList, world, del_t)
            errorCurrent = SQRT(SUM(solver%getError_tridiag_Ampere(world, del_t)**2))
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
        type(potentialSolver), intent(in out) :: solver
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
        type(potentialSolver), intent(in out) :: solver
        type(Particle), intent(in out) :: particleList(:)
        type(Domain), intent(in) :: world
        integer(int32), intent(in) :: maxIter
        real(real64), intent(in) :: del_t, eps_r
        real(real64) :: initialR, sumPastResiduals, initialNorm
        real(real64) :: Residual_k(NumberXNodes-2, m_Anderson+1), phi_k(NumberXNodes-2, m_Anderson+1), fitMat(NumberXNodes-2, m_Anderson), normResidual(m_Anderson+1)
        integer(int32) :: lwork!, ipar(2), itrmf
        real(real64) :: work(MAX(m_Anderson, (NumberXNodes -2)) + MIN((NumberXNodes -2),m_Anderson)), alpha(MAX(m_Anderson, (NumberXNodes -2)))!, fcurSolver(NumberXNodes-2)
        integer(int32) :: i, j, index, m_k, info, ldb
        ldb = MAX(m_Anderson, (NumberXNodes -2))
        lwork= MIN((NumberXNodes -2),m_Anderson) + ldb
        phi_k(:,1) = solver%phi(2:NumberXNodes-1)
        call solver%depositJ(particleList, world, del_t)
        initialNorm = SQRT(SUM(solver%phi(2:NumberXNodes-1)**2))
        call solver%solve_tridiag_Ampere(world, del_t)
        phi_k(:,2) = solver%phi_f(2:NumberXNodes-1)
        initialR = SQRT(SUM((solver%phi_f(2:NumberXNodes-1) - phi_k(:,1))**2))
        normResidual(1) = initialR
        Residual_k(:,1) = phi_k(:,2) - phi_k(:,1)
        !print *, "Initial norm is:", initialR
        do i = 1, maxIter
            index = MODULO(i, m_Anderson+1) + 1
            m_k = MIN(i, m_Anderson)
            ldb = MAX(m_k, (NumberXNodes -2))
            call solver%depositJ(particleList, world, del_t)
            call solver%solve_tridiag_Ampere(world, del_t)
            Residual_k(:, index) = solver%phi_f(2:NumberXNodes-1) - phi_k(:,index)
            normResidual(index) = SQRT(SUM(Residual_k(:, index)**2))
            if (normResidual(index) < eps_r*(initialR + initialNorm)) then
                call solver%moveParticles(particleList, world, del_t)
                solver%phi = solver%phi_f
                exit
            end if
            if (i > m_Anderson) then
                if (m_Anderson > 1) then
                    sumPastResiduals = normResidual(MODULO(i-1, m_Anderson+1) + 1) &
                    + normResidual(MODULO(i-2, m_Anderson+1) + 1)
                    sumPastResiduals = sumPastResiduals/2.0d0
                else
                    sumPastResiduals = normResidual(MODULO(i-1, m_Anderson+1) + 1)
                end if
                if (normResidual(index) > sumPastResiduals) then
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
        

    end subroutine solveDivAmpereAnderson

    subroutine adaptiveSolveDivAmpereAnderson(solver, particleList, world, del_t, maxIter, eps_r)
        ! Solve for divergence of ampere's law with picard
        ! cut del_t in two if non-convergence after maxIter, repeat until convergence
        type(potentialSolver), intent(in out) :: solver
        type(Particle), intent(in out) :: particleList(:)
        type(Domain), intent(in) :: world
        integer(int32), intent(in) :: maxIter
        real(real64), intent(in) :: del_t, eps_r
        real(real64) :: remainDel_t, currDel_t!, adaptiveJ(NumberXNodes-1)
        call solveDivAmpereAnderson(solver, particleList, world, del_t, maxIter, eps_r)
        remainDel_t = del_t  
        if (iterNumPicard == maxIter) then
            amountTimeSplits = amountTimeSplits + 1
            do while (iterNumPicard == maxIter)
                !print *, "Reducing time step"
                iterNumAdaptiveSteps = 0
                currDel_t = remainDel_t
                !adaptiveJ = 0.0d0
                do while (iterNumPicard == maxIter)
                    currDel_t = currDel_t/2.0d0
                    iterNumAdaptiveSteps = iterNumAdaptiveSteps + 1
                    if (iterNumAdaptiveSteps > 4) then
                        stop "ALREADY REDUCED TIME STEP MORE THAN 3 TIMES, REDUCE INITIAL TIME STEP!!!"
                    end if
                    call solveDivAmpereAnderson(solver, particleList, world, currDel_t, maxIter, eps_r)   
                end do
                !adaptiveJ = adaptiveJ + solver%J * currDel_t/del_t
                remainDel_t = remainDel_t - currDel_t 
                call solveDivAmpereAnderson(solver, particleList, world, remainDel_t, maxIter, eps_r)
            end do
            !adaptiveJ = adaptiveJ + solver%J * remainDel_t/del_t
            !solver%J = adaptiveJ
            
        end if
    end subroutine adaptiveSolveDivAmpereAnderson

    subroutine solvePotential(solver, particleList, world, del_t, maxIter, eps_r)
        type(potentialSolver), intent(in out) :: solver
        type(Particle), intent(in out) :: particleList(:)
        type(Domain), intent(in) :: world
        integer(int32), intent(in) :: maxIter
        real(real64), intent(in) :: del_t, eps_r
        SELECT CASE (solverType)
        CASE(0)
            call adaptiveSolveDivAmpereAnderson(solver, particleList, world, del_t, maxIter, eps_r)
        CASE(1)
            call adaptiveSolveDivAmpereJFNK(del_t, maxIter, eps_r)
        CASE default
            print *, "Solver type doesn't exit!"
            stop
        END SELECT
    end subroutine solvePotential

    ! -------------------- JFNK functions ---------------------------------------
    subroutine funcNitsol(n, xcur, fcur, rpar, ipar, itrmf)
        ! Use solver and whatnot as global inputs, I'm sure as hell not combining all data into one rpar and then creating new routines on those!
        integer(int32), intent(in) :: n
        integer(int32), intent(in out) :: itrmf, ipar(*)
        real(real64), intent(in) :: xcur(n)
        real(real64), intent(in) :: rpar
        real(real64), intent(in out) :: fcur(n)
        real(real64) :: d(n)
        globalSolver%phi_f(2:NumberXNodes-1) = xcur
        call globalSolver%depositJ(globalParticleList, globalWorld, rpar)
        d = globalSolver%getError_tridiag_Ampere(globalWorld, rpar)
        call solve_tridiag(n, globalSolver%a_tri, globalSolver%c_tri, globalSolver%b_tri, d, fcur)
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
            call solve_tridiag(n, globalSolver%a_tri, globalSolver%c_tri, globalSolver%b_tri, v, z)
        end if
        itrmjv = 0
    end subroutine jacNitsol

    subroutine solveJFNK(del_t, maxIter, eps_r)
        integer(int32), intent(in) :: maxIter
        real(real64), intent(in) :: del_t, eps_r
        real(real64) :: initialNorm
        integer(int32) :: info(6), iterm, ipar(2), itrmf
        real(real64) :: fcurSolver(NumberXNodes-2), rworkSolver((NumberXNodes-2)*(m_Anderson+5)+m_Anderson*(m_Anderson+3)), xcurSolver(NumberXNodes-2)

        ! Set Nitsol parameters
        iterm = 0
        xcurSolver = globalSolver%phi(2:NumberXNodes-1)
        call funcNitsol(NumberXNodes-2, xcurSolver, fcurSolver, del_t, ipar, itrmf)
        initialNorm = normFunc(NumberXNodes-2, fcurSolver, 1)
        !print *, "initial norm is:", initialNorm
        call nitsol(NumberXNodes-2, xcurSolver, funcNitsol, jacNitsol, eps_r*initialNorm, eps_r,inputJFNK, info, rworkSolver, del_t, ipar, iterm, ddot, dnrm2)
        SELECT CASE (iterm)
        CASE(0)
            iterNumPicard = info(4)
            call globalSolver%moveParticles(globalParticleList, globalWorld, del_t)
            globalSolver%phi = globalSolver%phi_f
        CASE(1)
            iterNumPicard = maxIter
        CASE default
            print *, "Nitsol error with iterm == ", iterm
            stop
        END SELECT
    end subroutine solveJFNK

    subroutine adaptiveSolveDivAmpereJFNK(del_t, maxIter, eps_r)
        ! Solve for divergence of ampere's law with picard
        ! cut del_t in two if non-convergence after maxIter, repeat until convergence
        integer(int32), intent(in) :: maxIter
        real(real64), intent(in) :: del_t, eps_r
        real(real64) :: remainDel_t, currDel_t!, adaptiveJ(NumberXNodes-1)
        call solveJFNK(del_t, maxIter, eps_r)
        remainDel_t = del_t  
        if (iterNumPicard == maxIter) then
            amountTimeSplits = amountTimeSplits + 1
            do while (iterNumPicard == maxIter)
                !print *, "Reducing time step"
                iterNumAdaptiveSteps = 0
                currDel_t = remainDel_t
                !adaptiveJ = 0.0d0
                do while (iterNumPicard == maxIter)
                    currDel_t = currDel_t/2.0d0
                    iterNumAdaptiveSteps = iterNumAdaptiveSteps + 1
                    if (iterNumAdaptiveSteps > 4) then
                        stop "ALREADY REDUCED TIME STEP MORE THAN 3 TIMES, REDUCE INITIAL TIME STEP!!!"
                    end if
                    call solveJFNK(currDel_t, maxIter, eps_r)
                end do
                !adaptiveJ = adaptiveJ + globalSolver%J * currDel_t/del_t
                remainDel_t = remainDel_t - currDel_t 
                call solveJFNK(remainDel_t, maxIter, eps_r)
            end do
            ! adaptiveJ = adaptiveJ + globalSolver%J * remainDel_t/del_t
            ! globalSolver%J = adaptiveJ
        end if
    end subroutine adaptiveSolveDivAmpereJFNK


end module mod_nonLinSolvers