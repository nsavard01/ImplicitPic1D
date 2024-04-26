module mod_nonLinSolvers

    use iso_fortran_env, only: int32, real64
    use constants
    use mod_BasicFunctions
    use mod_particle
    use mod_domain
    use mod_potentialSolver
    use mod_particleMover
    use mod_Scheme
    use mod_nitsol
    implicit none
    
    interface
        subroutine dgels(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
            Character*1, intent(in) :: trans
            integer, intent(in) :: m, n, nrhs, lda, ldb, lwork
            integer, intent(out) :: info
            real(kind=8), intent(in out) :: a(lda, *), b(ldb, *), work(lwork)

        end subroutine dgels

    end interface

    ! Initialize objects needed
    type(Domain) :: globalWorld
    type(Particle), allocatable :: globalParticleList(:)
    type(potentialSolver) :: globalSolver
    integer(int32), protected :: iterNumPicard, iterNumParticle, iterNumAdaptiveSteps, amountTimeSplits
    integer(int32), private :: maxIter, solverType, m_Anderson
    real(real64), private :: Beta_k, eps_r, eps_a

    !allocatable arrays for Anderson
    real(real64), private, allocatable :: Residual_k(:, :), EField_k(:, :), fitMat(:, :)

contains

    ! ---------------- Initial Poisson Solver -------------------------------------------------

    subroutine solveInitialPotential(solver, particleList, world, timeCurrent)
        ! Solve for initial potential
        class(potentialSolver), intent(in out) :: solver
        type(Particle), intent(in out) :: particleList(numberChargedParticles)
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: timeCurrent
        call depositRho(solver%rho, particleList, world)
        call solver%solve_tridiag_Poisson(world, timeCurrent)
        ! Assume only use potential solver once, then need to generate matrix for Div-Ampere

    end subroutine solveInitialPotential

    ! Non-linear solver stuff -------------------------------------------------------------

    subroutine initializeSolver()
        integer(int32) :: io
        print *, "Reading non-linear solver inputs:"
        open(10,file='../InputData/SolverState.inp', IOSTAT=io)
        read(10, *, IOSTAT = io) eps_a
        read(10, *, IOSTAT = io) eps_r
        read(10, *, IOSTAT = io) solverType
        read(10, *, IOSTAT = io) m_Anderson
        read(10, *, IOSTAT = io) Beta_k
        read(10, *, IOSTAT = io) maxIter
        close(10)
        print *, 'Absolute error:', eps_a
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
            allocate(Residual_k(NumberXHalfNodes, m_Anderson+1), EField_k(NumberXHalfNodes, m_Anderson+1), fitMat(NumberXHalfNodes, m_Anderson) )
        CASE(1)
            call initializeNitsol(maxIter, m_Anderson, NumberXHalfNodes)
        END SELECT

    end subroutine initializeSolver

    

    ! ----------- Picard with Anderson Acceleration -------------------------------

    subroutine solveDivAmpereAnderson(solver, particleList, world, del_t)
        ! Solve for divergence of ampere using picard iterations
        type(potentialSolver), intent(in out) :: solver
        type(Particle), intent(in out) :: particleList(numberChargedParticles)
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t
        real(real64) :: eps_tol, sumPastResiduals
        real(real64) :: normResidual(m_Anderson+1), alpha(m_Anderson+1)
        integer(int32) :: i, j, index, m_k
        print *, 'Entering DivAmpere:'
        print *, ''
        
        EField_k(:,1) = solver%EField_f
        call depositJ(solver, particleList, world, del_t)
        call solver%solve_tridiag_Ampere(world, del_t)
        EField_k(:,2) = solver%EField_f
        Residual_k(:,1) = EField_k(:,2) - EField_k(:,1)
        normResidual(1) = SQRT(SUM(Residual_k(:,1)**2))
        eps_tol = eps_r * normResidual(1) + eps_a * SQRT(real(NumberXHalfNodes))
        !print *, "Initial norm is:", initialR
        do i = 1, maxIter
            index = MODULO(i, m_Anderson+1) + 1
            m_k = MIN(i, m_Anderson)
            call depositJ(solver,particleList, world, del_t)
            call solver%solve_tridiag_Ampere(world, del_t)
            Residual_k(:, index) = solver%EField_f - EField_k(:,index)
            normResidual(index) = SQRT(SUM(Residual_k(:, index)**2))
            print *, 'NormResidual:', normResidual(index)
            if (normResidual(index) < eps_tol) then
                call moveParticles(solver,particleList, world, del_t)
                iterNumPicard = i
                exit
            end if
            ! if (i > m_Anderson) then
            !     if (m_Anderson > 1) then
            !         sumPastResiduals = normResidual(MODULO(i-1, m_Anderson+1) + 1) &
            !         + normResidual(MODULO(i-2, m_Anderson+1) + 1)
            !         sumPastResiduals = sumPastResiduals/2.0d0
            !     else
            !         sumPastResiduals = normResidual(MODULO(i-1, m_Anderson+1) + 1)
            !     end if
            !     if (normResidual(index) > sumPastResiduals) then
            !         iterNumPicard = maxIter
            !         exit
            !     end if
            ! end if
            
            do j = 0, m_k-1
                fitMat(:,j+1) =  Residual_k(:, index) - Residual_k(:, MODULO(i - m_k + j, m_Anderson+1) + 1)
            end do
            !call dgels('N', NumberXNodes, m_k, 1, fitMat(:, 1:m_k), NumberXNodes, alpha(1:ldb), ldb, work, lwork, info)
            alpha(1:m_k) = solveNormalEquation(fitMat(:, 1:m_k), Residual_k(:, index), NumberXHalfNodes, m_k)
            alpha(m_k+1) = 1.0d0 - SUM(alpha(1:m_k)) 
            EField_k(:, MODULO(i+1, m_Anderson+1) + 1) = alpha(1) * (Beta_k*Residual_k(:, MODULO(i-m_k, m_Anderson+1) + 1) + EField_k(:, MODULO(i-m_k,m_Anderson+1) + 1))
            do j=1, m_k
                EField_k(:, MODULO(i+1, m_Anderson+1) + 1) = EField_k(:, MODULO(i+1, m_Anderson+1) + 1) + alpha(j + 1) * (Beta_k*Residual_k(:, MODULO(i-m_k + j, m_Anderson+1) + 1) + EField_k(:, MODULO(i-m_k + j, m_Anderson+1) + 1))
            end do
            solver%EField_f = EField_k(:, MODULO(i+1, m_Anderson+1) + 1)
        end do
       
    end subroutine solveDivAmpereAnderson

    subroutine adaptiveSolveDivAmpereAnderson(solver, particleList, world, del_t, remainDel_t, currDel_t, timeCurrent)
        ! Solve for divergence of ampere's law with picard
        ! cut del_t in two if non-convergence after maxIter, repeat until convergence
        type(potentialSolver), intent(in out) :: solver
        type(Particle), intent(in out) :: particleList(numberChargedParticles)
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t, timeCurrent
        real(real64), intent(in out) :: remainDel_t, currDel_t
        currDel_t = remainDel_t
        if (solver%RF_bool) then
            ! if RF, change value of future phi values at RF boundary
            call solver%setRFVoltage(world, timeCurrent + remainDel_t)
        end if
        call solveDivAmpereAnderson(solver, particleList, world, remainDel_t) 
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
                if (solver%RF_bool) then
                    ! if RF, change value of future phi values at RF boundary
                    call solver%setRFVoltage(world, timeCurrent + currDel_t)
                end if 
                call solveDivAmpereAnderson(solver, particleList, world, currDel_t)  
            end do 
            remainDel_t = remainDel_t - currDel_t 
        end if
    end subroutine adaptiveSolveDivAmpereAnderson

    subroutine solvePotential(solver, particleList, world, del_t, remainDel_t, currDel_t, timeCurrent)
        type(potentialSolver), intent(in out) :: solver
        type(Particle), intent(in out) :: particleList(numberChargedParticles)
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t, timeCurrent
        real(real64), intent(in out) :: remainDel_t, currDel_t
        ! make future phi now current phi
        call solver%resetVoltage()
        SELECT CASE (solverType)
        CASE(0)
            call adaptiveSolveDivAmpereAnderson(solver, particleList, world, del_t, remainDel_t, currDel_t, timeCurrent)
        CASE(1)
            call adaptiveSolveDivAmpereJFNK(del_t, remainDel_t, currDel_t, timeCurrent)
        CASE default
            print *, "Solver type doesn't exit!"
            stop
        END SELECT
    end subroutine solvePotential

    ! -------------------- JFNK functions ---------------------------------------
    subroutine funcNitsol(n, xcur, fcur, rpar, ipar, itrmf)
        ! Use solver and whatnot as global inputs, I'm sure as hell not combining all data into one rpar and then creating new routines on those!
        integer, intent(in) :: n
        integer, intent(in out) :: itrmf, ipar(*)
        double precision, intent(in) :: xcur(n)
        double precision, intent(in out) :: rpar(*), fcur(n)
        !real(real64) :: d(n)
        globalSolver%EField_f = xcur
        call depositJ(globalSolver, globalParticleList, globalWorld, rpar(1))
        call globalSolver%solve_tridiag_Ampere(globalWorld, rpar(1))
        fcur = xcur - globalSolver%EField_f
        ! fcur = globalSolver%getError_tridiag_Ampere(globalWorld, rpar(1))
        itrmf = 0

    end subroutine funcNitsol

    subroutine jacNitsol(n, xcur, fcur, ijob, v, z, rpar, ipar, itrmjv) 
        ! If analytical jacobian matrix-vector product or preconditioner needed
        integer, intent(in) :: ijob, n
        integer, intent(in out) :: itrmjv, ipar(*)
        double precision, intent(in) :: fcur(n), v(n), xcur(n)
        double precision, intent(in out) :: rpar(*), z(n)
        if (ijob == 0) then
            continue
        else if (ijob == 1) then
            call solve_tridiag(n, globalSolver%a_tri, globalSolver%c_tri, globalSolver%b_tri, v, z)
        end if
        itrmjv = 0
    end subroutine jacNitsol

    subroutine solveJFNK(del_t)
        real(real64), intent(in) :: del_t
        integer(int32) :: info(6), iterm, ipar(2), itrmf
        real(real64) :: fcurSolver(NumberXHalfNodes), xcurSolver(NumberXHalfNodes), rpar(1)

        ! Set Nitsol parameters
        iterm = 0
        xcurSolver = globalSolver%EField_f
        rpar(1) = del_t
        !dnrm2(NumberXNodes, fcurSolver, 1)
        !print *, "initial norm is:", initialNorm
        !eps_r * SQRT(real(NumberXNodes))
        call nitsol(NumberXHalfNodes, xcurSolver, funcNitsol, jacNitsol, eps_a, eps_r, 1.d-12, info, rpar, ipar, iterm)
        SELECT CASE (iterm)
        CASE(0)
            iterNumPicard = info(4)
            call moveParticles(globalSolver, globalParticleList, globalWorld, del_t)
        CASE(1,5,6)
            iterNumPicard = maxIter
        CASE default
            print *, "Nitsol error with iterm == ", iterm
            stop
        END SELECT
    end subroutine solveJFNK

    subroutine adaptiveSolveDivAmpereJFNK(del_t, remainDel_t, currDel_t, timeCurrent)
        ! Solve for divergence of ampere's law with picard
        ! cut del_t in two if non-convergence after maxIter, repeat until convergence
        real(real64), intent(in) :: del_t, timeCurrent
        real(real64), intent(in out) :: remainDel_t, currDel_t
        currDel_t = remainDel_t
        if (globalSolver%RF_bool) then
            ! if RF, change value of future phi values at RF boundary
            call globalSolver%setRFVoltage(globalWorld, timeCurrent + remainDel_t)
        end if
        call solveJFNK(remainDel_t)
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
                if (globalSolver%RF_bool) then
                    ! if RF, change value of future phi values at RF boundary
                    call globalSolver%setRFVoltage(globalWorld, timeCurrent + currDel_t)
                end if
                call solveJFNK(currDel_t)
            end do
            remainDel_t = remainDel_t - currDel_t 
        end if
    end subroutine adaptiveSolveDivAmpereJFNK

    subroutine writeSolverState(dirName)
        character(*), intent(in) :: dirName
        open(15,file=dirName//'/SolverState.dat')
        write(15,'("Solver Type, eps_a, eps_r, m_Anderson, Beta_k, maximum iterations")')
        write(15,"(1(I3.3, 1x), 2(es16.8,1x), 1(I3.3, 1x), 1(es16.8,1x), 1(I3.3, 1x))") solverType, eps_a, eps_r, m_Anderson, Beta_k, maxIter
        close(15)

    end subroutine writeSolverState

end module mod_nonLinSolvers