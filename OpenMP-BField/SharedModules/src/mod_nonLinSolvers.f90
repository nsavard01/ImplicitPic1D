module mod_nonLinSolvers

    use iso_fortran_env, only: int32, real64
    use constants
    use mod_BasicFunctions
    use mod_particle
    use mod_domain
    use mod_potentialSolver
    use mod_particleMover
    use mod_Scheme
    implicit none
    
    interface
        subroutine dgels(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
            Character*1, intent(in) :: trans
            integer, intent(in) :: m, n, nrhs, lda, ldb, lwork
            integer, intent(out) :: info
            real(kind=8), intent(in out) :: a(lda, *), b(ldb, *), work(lwork)

        end subroutine dgels

        subroutine nitsol(n, x, f, jacv, ftol, stptol, input, info, rwork, rpar, ipar, iterm, dinpr, dnorm)
            integer, intent(in)                   :: n
            real(kind=8), dimension(n), intent(inout) :: x
            external :: f, jacv
            real(kind=8), intent(in)                  :: ftol
            real(kind=8), intent(in)                  :: stptol
            integer, dimension(10), intent(in)    :: input
            integer, dimension(6), intent(out)    :: info
            real(kind=8), dimension(*), intent(inout) :: rwork
            real(kind=8), intent(in) :: rpar(:)
            integer, dimension(*), intent(in out)  :: ipar
            integer, intent(out)                  :: iterm
            
            real(kind=8), external :: dinpr, dnorm

        end subroutine nitsol

        function ddot(n, x, sx, y, sy)
            implicit none
            integer, intent(in)                :: n
            real(kind=8), dimension(*), intent(in) :: x
            integer, intent(in)                :: sx
            real(kind=8), dimension(*), intent(in) :: y
            integer, intent(in)                :: sy
            real(kind=8)                           :: ddot
          end function ddot
      
          function dnrm2(n, x, sx)
            implicit none
            integer, intent(in)                :: n
            real(kind=8), dimension(*), intent(in) :: x
            integer, intent(in)                :: sx
            real(kind = 8)                           :: dnrm2
          end function dnrm2

    end interface

    ! Initialize objects needed
    type(Domain) :: globalWorld
    type(Particle), allocatable :: globalParticleList(:)
    type(potentialSolver) :: globalSolver
    integer(int32) :: iterNumPicard, iterNumParticle, iterNumAdaptiveSteps, m_Anderson, amountTimeSplits, solverType
    real(real64) :: Beta_k

    !allocatable arrays for JFNK or Anderson
    integer(int32), private, allocatable :: inputJFNK(:)
    real(real64), private, allocatable :: Residual_k(:, :), phi_k(:, :), fitMat(:, :), work(:), alpha(:)

    ! Common blocks for nitsol
    integer iplvl, ipunit
    common /nitprint/ iplvl, ipunit
    double precision choice1_exp, choice2_exp, choice2_coef
    double precision eta_cutoff, etamax
    double precision thmin, thmax, etafixed

    common /nitparam/ choice1_exp, choice2_exp, choice2_coef, eta_cutoff, etamax, thmin, thmax, etafixed

    integer instep, newstep, krystat
    double precision avrate, fcurnrm
    common /nitinfo/ avrate, fcurnrm, instep, newstep, krystat

contains

    ! ---------------- Initial Poisson Solver -------------------------------------------------

    subroutine solveInitialPotential(solver, particleList, world, timeCurrent)
        ! Solve for initial potential
        class(potentialSolver), intent(in out) :: solver
        type(Particle), intent(in out) :: particleList(:)
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: timeCurrent
        call depositRho(solver%rho, particleList, world)
        call solver%solve_tridiag_Poisson(world, timeCurrent)
        ! Assume only use potential solver once, then need to generate matrix for Div-Ampere

    end subroutine solveInitialPotential

    ! Non-linear solver stuff -------------------------------------------------------------

    subroutine initializeSolver(eps_r, solverType, m_Anderson, Beta_k, maxIter)
        real(real64), intent(in out) :: eps_r, Beta_k
        integer(int32), intent(in out) :: m_Anderson, solverType, maxIter
        integer(int32) :: io
        print *, "Reading non-linear solver inputs:"
        open(10,file='../../SharedModules/InputData/SolverState.inp', IOSTAT=io)
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
        CASE(0)
            allocate(Residual_k(NumberXNodes, m_Anderson+1), phi_k(NumberXNodes, m_Anderson+1), fitMat(NumberXNodes, m_Anderson), &
            work(MAX(m_Anderson, (NumberXNodes)) + MIN((NumberXNodes),m_Anderson)), alpha(MAX(m_Anderson, (NumberXNodes))))
        CASE(1)
            allocate(inputJFNK(10))
            iplvl = 4 ! print level
            inputJFNK = 0
            inputJFNK(1) = maxIter ! maximum iterations
            inputJFNK(2) = 0 !ijacv
            inputJFNK(3) = 0 ! krylov solver
            inputJFNK(4) = m_Anderson ! maximum krylov subspace dimension
            inputJFNK(5) = 0 !ipre
            inputJFNK(9) = -1 !number backtracks
            inputJFNK(6) = m_Anderson*10
            inputJFNK(10) = 2 ! eta with gamma and alpha
            etamax = 0.8d0 ! eta max
            choice2_exp = 1.5d0 ! alpha
            choice2_coef = 0.9d0 ! gamma
        END SELECT

    end subroutine initializeSolver

    

    ! ----------- Picard with Anderson Acceleration -------------------------------

    subroutine solveDivAmpereAnderson(solver, particleList, world, del_t, maxIter, eps_r, timeCurrent)
        ! Solve for divergence of ampere using picard iterations
        type(potentialSolver), intent(in out) :: solver
        type(Particle), intent(in out) :: particleList(:)
        type(Domain), intent(in) :: world
        integer(int32), intent(in) :: maxIter
        real(real64), intent(in) :: del_t, eps_r, timeCurrent
        real(real64) :: initialR, sumPastResiduals, initialNorm, timeFuture
        real(real64) :: normResidual(m_Anderson+1)
        integer(int32) :: lwork!, ipar(2), itrmf
        integer(int32) :: i, j, index, m_k, info, ldb
        timeFuture = timeCurrent + del_t
        ldb = MAX(m_Anderson, (NumberXNodes))
        lwork= MIN((NumberXNodes),m_Anderson) + ldb
        phi_k(:,1) = solver%phi
        call depositJ(solver, particleList, world, del_t)
        initialNorm = SQRT(SUM(solver%phi**2))
        call solver%solve_tridiag_Ampere(world, del_t, timeFuture)
        phi_k(:,2) = solver%phi_f
        initialR = SQRT(real(NumberXNodes))!SQRT(SUM((solver%phi_f - phi_k(:,1))**2))
        normResidual(1) = initialR
        Residual_k(:,1) = phi_k(:,2) - phi_k(:,1)
        !print *, "Initial norm is:", initialR
        do i = 1, maxIter
            index = MODULO(i, m_Anderson+1) + 1
            m_k = MIN(i, m_Anderson)
            ldb = MAX(m_k, (NumberXNodes))
            call depositJ(solver,particleList, world, del_t)
            call solver%solve_tridiag_Ampere(world, del_t, timeFuture)
            Residual_k(:, index) = solver%phi_f - phi_k(:,index)
            normResidual(index) = SQRT(SUM(Residual_k(:, index)**2))
            !print *, normResidual(index)
            if (normResidual(index) < eps_r*(initialR)) then
                call moveParticles(solver,particleList, world, del_t)
                solver%phi = solver%phi_f
                iterNumPicard = i
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
                    exit
                end if
            end if
            do j = 0, m_k-1
                fitMat(:,j+1) = Residual_k(:, MODULO(i - m_k + j, m_Anderson+1) + 1) - Residual_k(:, index)
            end do
            alpha(1:NumberXNodes) = -Residual_k(:, index)
            call dgels('N', NumberXNodes, m_k, 1, fitMat(:, 1:m_k), NumberXNodes, alpha(1:ldb), ldb, work, lwork, info)
            if (info /= 0) then
                print *, "Issue with minimization procedure dgels in Anderson Acceleration!"
                stop
            end if
            alpha(m_k+1) = 1.0d0 - SUM(alpha(1:m_k)) 
            phi_k(:, MODULO(i+1, m_Anderson+1) + 1) = alpha(1) * (Beta_k*Residual_k(:, MODULO(i-m_k, m_Anderson+1) + 1) + phi_k(:, MODULO(i-m_k,m_Anderson+1) + 1))
            do j=1, m_k
                phi_k(:, MODULO(i+1, m_Anderson+1) + 1) = phi_k(:, MODULO(i+1, m_Anderson+1) + 1) + alpha(j + 1) * (Beta_k*Residual_k(:, MODULO(i-m_k + j, m_Anderson+1) + 1) + phi_k(:, MODULO(i-m_k + j, m_Anderson+1) + 1))
            end do
            solver%phi_f = phi_k(:, MODULO(i+1, m_Anderson+1) + 1)
        end do
    end subroutine solveDivAmpereAnderson

    subroutine adaptiveSolveDivAmpereAnderson(solver, particleList, world, del_t, remainDel_t, currDel_t, maxIter, eps_r, timeCurrent)
        ! Solve for divergence of ampere's law with picard
        ! cut del_t in two if non-convergence after maxIter, repeat until convergence
        type(potentialSolver), intent(in out) :: solver
        type(Particle), intent(in out) :: particleList(:)
        type(Domain), intent(in) :: world
        integer(int32), intent(in) :: maxIter
        real(real64), intent(in) :: del_t, eps_r, timeCurrent
        real(real64), intent(in out) :: remainDel_t, currDel_t
        currDel_t = remainDel_t
        call solveDivAmpereAnderson(solver, particleList, world, remainDel_t, maxIter, eps_r, timeCurrent) 
        if (iterNumPicard < maxIter) then
            remainDel_t = del_t  
        else
            amountTimeSplits = amountTimeSplits + 1
            iterNumAdaptiveSteps = 0
            !adaptiveJ = 0.0d0
            do while (iterNumPicard == maxIter)
                currDel_t = currDel_t/2.0d0
                iterNumAdaptiveSteps = iterNumAdaptiveSteps + 1
                if (iterNumAdaptiveSteps > 4) then
                    stop "ALREADY REDUCED TIME STEP MORE THAN 3 TIMES, REDUCE INITIAL TIME STEP!!!"
                end if
                call solveDivAmpereAnderson(solver, particleList, world, currDel_t, maxIter, eps_r, timeCurrent)   
            end do 
            remainDel_t = remainDel_t - currDel_t 
        end if
    end subroutine adaptiveSolveDivAmpereAnderson

    subroutine solvePotential(solver, particleList, world, del_t, remainDel_t, currDel_t, maxIter, eps_r, timeCurrent)
        type(potentialSolver), intent(in out) :: solver
        type(Particle), intent(in out) :: particleList(:)
        type(Domain), intent(in) :: world
        integer(int32), intent(in) :: maxIter
        real(real64), intent(in) :: del_t, eps_r, timeCurrent
        real(real64), intent(in out) :: remainDel_t, currDel_t
        SELECT CASE (solverType)
        CASE(0)
            call adaptiveSolveDivAmpereAnderson(solver, particleList, world, del_t, remainDel_t, currDel_t, maxIter, eps_r, timeCurrent)
        CASE(1)
            call adaptiveSolveDivAmpereJFNK(del_t, remainDel_t, currDel_t, maxIter, eps_r, timeCurrent)
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
        real(real64), intent(in) :: rpar(2)
        real(real64), intent(in out) :: fcur(n)
        !real(real64) :: d(n)
        globalSolver%phi_f = xcur
        call depositJ(globalSolver, globalParticleList, globalWorld, rpar(1))
        ! d = globalSolver%getError_tridiag_Ampere(globalWorld, rpar)
        ! call solve_tridiag(n, globalSolver%a_tri, globalSolver%c_tri, globalSolver%b_tri, d, fcur)
        call globalSolver%solve_tridiag_Ampere(globalWorld, rpar(1), rpar(2))
        fcur = xcur - globalSolver%phi_f
        itrmf = 0

    end subroutine funcNitsol

    subroutine jacNitsol(n, xcur, fcur, ijob, v, z, rpar, ipar, itrmjv) 
        ! If analytical jacobian matrix-vector product or preconditioner needed
        integer, intent(in) :: ijob
        integer, intent(in out) :: itrmjv
        integer, intent(in) :: n

        integer, intent(in out) :: ipar(*)

        real(real64), intent(in) :: fcur(n)
        real(real64), intent(in) :: rpar(2)
        real(real64), intent(in) :: v(n)
        real(real64), intent(in) ::  xcur(n)
        real(real64), intent(in out) :: z(n)
        if (ijob == 0) then
            print *, "in analytical jacobian if block"
        else if (ijob == 1) then
            call solve_tridiag(n, globalSolver%a_tri, globalSolver%c_tri, globalSolver%b_tri, v, z)
            print *, 'In solve tridiag?'
        end if
        itrmjv = 0
    end subroutine jacNitsol

    subroutine solveJFNK(del_t, maxIter, eps_r, timeCurrent)
        integer(int32), intent(in) :: maxIter
        real(real64), intent(in) :: del_t, eps_r, timeCurrent
        real(real64) :: initialNorm
        integer(int32) :: info(6), iterm, ipar(2), itrmf
        real(real64) :: fcurSolver(NumberXNodes), rworkSolver((NumberXNodes)*(m_Anderson+5)+m_Anderson*(m_Anderson+3)), xcurSolver(NumberXNodes), rpar(2)

        ! Set Nitsol parameters
        iterm = 0
        rpar(1) = del_t
        rpar(2) = timeCurrent
        xcurSolver = globalSolver%phi
        call funcNitsol(NumberXNodes, xcurSolver, fcurSolver, rpar, ipar, itrmf)
        initialNorm = dnrm2(NumberXNodes, fcurSolver, 1)
        !print *, "initial norm is:", initialNorm
        call nitsol(NumberXNodes, xcurSolver, funcNitsol, jacNitsol, eps_r * SQRT(real(NumberXNodes)), 1.d-20,inputJFNK, info, rworkSolver, rpar, ipar, iterm, ddot, dnrm2)
        SELECT CASE (iterm)
        CASE(0)
            iterNumPicard = info(4)
            call moveParticles(globalSolver, globalParticleList, globalWorld, del_t)
            globalSolver%phi = globalSolver%phi_f
        CASE(1)
            iterNumPicard = maxIter
        CASE default
            print *, "Nitsol error with iterm == ", iterm
            stop
        END SELECT
        print *, 'iterNumPicard is:', iterNumPicard
    end subroutine solveJFNK

    subroutine adaptiveSolveDivAmpereJFNK(del_t, remainDel_t, currDel_t, maxIter, eps_r, timeCurrent)
        ! Solve for divergence of ampere's law with picard
        ! cut del_t in two if non-convergence after maxIter, repeat until convergence
        integer(int32), intent(in) :: maxIter
        real(real64), intent(in) :: del_t, eps_r, timeCurrent
        real(real64), intent(in out) :: remainDel_t, currDel_t
        currDel_t = remainDel_t
        call solveJFNK(remainDel_t, maxIter, eps_r, timeCurrent)
        if (iterNumPicard < maxIter) then
            remainDel_t = del_t
        else
            amountTimeSplits = amountTimeSplits + 1
            iterNumAdaptiveSteps = 0
            do while (iterNumPicard == maxIter)
                currDel_t = currDel_t/2.0d0
                iterNumAdaptiveSteps = iterNumAdaptiveSteps + 1
                if (iterNumAdaptiveSteps > 4) then
                    stop "ALREADY REDUCED TIME STEP MORE THAN 3 TIMES, REDUCE INITIAL TIME STEP!!!"
                end if
                call solveJFNK(currDel_t, maxIter, eps_r, timeCurrent)
            end do
            remainDel_t = remainDel_t - currDel_t 
        end if
    end subroutine adaptiveSolveDivAmpereJFNK


end module mod_nonLinSolvers