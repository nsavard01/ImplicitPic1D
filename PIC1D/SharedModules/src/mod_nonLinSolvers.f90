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

    ! Non-linear solvers for implicit PIC
    
    interface
        ! MKL subroutine to solve under or overdetermined linear system
        ! can be used 
        subroutine dgels(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
            Character*1, intent(in) :: trans
            integer, intent(in) :: m, n, nrhs, lda, ldb, lwork
            integer, intent(out) :: info
            real(kind=8), intent(in out) :: a(lda, *), b(ldb, *), work(lwork)

        end subroutine dgels

    end interface

    ! Initialize classes here since nitsol made to work with arrays of numbers instead of objects
    ! For this reason, just make global so nitsol can be used. Suffer potential consequences later
    ! Generally use AA anyway
    type(Domain) :: globalWorld
    type(Particle), allocatable :: globalParticleList(:)
    type(potentialSolver) :: globalSolver
    ! Initialize variables used for solver parameters and diagnostics
    integer(int32), protected :: iterNumPicard, iterNumParticle, iterNumAdaptiveSteps, amountTimeSplits
    integer(int64), protected :: potTimer, moverTimer
    integer(int32), private :: maxIter, solverType, m_Anderson
    real(real64), private :: Beta_k, eps_r, eps_a

    !allocatable arrays for Anderson Acceleration subroutines, these will be reused so have them initialized into heap arrays
    real(real64), private, allocatable :: Residual_k(:, :), phi_k(:, :), fitMat(:, :)

contains

    ! Non-linear solver stuff -------------------------------------------------------------

    subroutine initializeSolver()
        ! Initialize either AA or JFNK (nitsol) solver parameters
        integer(int32) :: io
        print *, "Reading non-linear solver inputs:"
        if (.not. restartBool) then
            open(10,file='../InputData/SolverState.inp', IOSTAT=io)
        else
            open(10,file=restartDirectory//'/InputData/SolverState.inp', IOSTAT=io)
        end if
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
            allocate(Residual_k(NumberXNodes, m_Anderson+1), phi_k(NumberXNodes, m_Anderson+1), fitMat(NumberXNodes, m_Anderson) )
        CASE(1)
            call initializeNitsol(maxIter, m_Anderson, NumberXNodes)
        END SELECT
        potTimer = 0
        moverTimer = 0
    end subroutine initializeSolver

    

    ! ----------- Picard with Anderson Acceleration -------------------------------

    subroutine solveDivAmpereAnderson(solver, particleList, world, del_t)
        ! Solve for divergence of ampere using picard iterations with anderson acceleration
        ! Allows for higher del_t than straight Picard
        type(potentialSolver), intent(in out) :: solver
        type(Particle), intent(in out) :: particleList(numberChargedParticles)
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t
        real(real64) :: eps_tol, sumPastResiduals
        real(real64) :: normResidual(m_Anderson+1), alpha(m_Anderson+1)
        integer(int32) :: i, j, index, m_k
        integer(int64) :: startTime, endTime
        
        ! put constaent source term into solver%rho so don't need to recalculate everytime
        call solver%makeConstSourceTerm(world)
        phi_k(:,1) = solver%phi_f
        call system_clock(startTime)
        call depositJ(solver, particleList, world, del_t)
        call system_clock(endTime)
        moverTimer = moverTimer + endTime - startTime
        call system_clock(startTime)
        call solver%solve_tridiag_Ampere(world, del_t)
        call system_clock(endTime)
        potTimer = potTimer + endTime - startTime
        ! Store solutions 
        phi_k(:,2) = solver%phi_f
        Residual_k(:,1) = phi_k(:,2) - phi_k(:,1)
        normResidual(1) = SQRT(SUM(Residual_k(:,1)**2))
        ! Find total convergence criterion eps_r * R_init + eps_a * sqrt(N_x)
        ! sqrt(N_x) so eps_a is properly normalized as we don't divied by sqrt(N_x) when evaluating RMS values of residual
        eps_tol = eps_r * normResidual(1) + eps_a * SQRT(real(NumberXNodes))
        do i = 1, maxIter
            ! store last m iterations with index cycling
            index = MODULO(i, m_Anderson+1) + 1
            m_k = MIN(i, m_Anderson)
            call system_clock(startTime)
            call depositJ(solver,particleList, world, del_t)
            call system_clock(endTime)
            moverTimer = moverTimer + endTime - startTime
            call system_clock(startTime)
            call solver%solve_tridiag_Ampere(world, del_t)
            call system_clock(endTime)
            potTimer = potTimer + endTime - startTime
            Residual_k(:, index) = solver%phi_f - phi_k(:,index)
            normResidual(index) = SQRT(SUM(Residual_k(:, index)**2))
            ! Check for convergence, done using (phi_k+1 - phi_k)**2 sum
            if (normResidual(index) < eps_tol) then
                ! Move particles
                call system_clock(startTime)
                call moveParticles(solver,particleList, world, del_t)
                call system_clock(endTime)
                moverTimer = moverTimer + endTime - startTime
                iterNumPicard = i+1
                exit
            end if
            if (i > m_Anderson) then
                !Stagnant if average previous 3 residuals > current residual
                if (m_Anderson > 3) then
                    sumPastResiduals = normResidual(MODULO(i-1, m_Anderson+1) + 1) &
                    + normResidual(MODULO(i-2, m_Anderson+1) + 1) + normResidual(MODULO(i-3, m_Anderson+1) + 1)
                    sumPastResiduals = sumPastResiduals/3.0d0
                else
                    sumPastResiduals = normResidual(MODULO(i-1, m_Anderson+1) + 1)
                end if
                ! If stagnating residual then exit so time step can be cut in half
                if (normResidual(index) > sumPastResiduals) then
                    iterNumPicard = maxIter
                    exit
                end if
            end if
            
            ! Form matrix for minimization procedure
            do j = 0, m_k-1
                fitMat(:,j+1) =  Residual_k(:, index) - Residual_k(:, MODULO(i - m_k + j, m_Anderson+1) + 1)
            end do
            !call dgels('N', NumberXNodes, m_k, 1, fitMat(:, 1:m_k), NumberXNodes, alpha(1:ldb), ldb, work, lwork, info)
            ! For small m_anderson quicker to solve normal equation rather than use dgels
            alpha(1:m_k) = solveNormalEquation(fitMat(:, 1:m_k), Residual_k(:, index), NumberXNodes, m_k)
            alpha(m_k+1) = 1.0d0 - SUM(alpha(1:m_k)) 
            ! update next phi_k
            phi_k(:, MODULO(i+1, m_Anderson+1) + 1) = alpha(1) * (Beta_k*Residual_k(:, MODULO(i-m_k, m_Anderson+1) + 1) + phi_k(:, MODULO(i-m_k,m_Anderson+1) + 1))
            do j=1, m_k
                phi_k(:, MODULO(i+1, m_Anderson+1) + 1) = phi_k(:, MODULO(i+1, m_Anderson+1) + 1) + alpha(j + 1) * (Beta_k*Residual_k(:, MODULO(i-m_k + j, m_Anderson+1) + 1) + phi_k(:, MODULO(i-m_k + j, m_Anderson+1) + 1))
            end do
            solver%phi_f = phi_k(:, MODULO(i+1, m_Anderson+1) + 1)
        end do
       
    end subroutine solveDivAmpereAnderson

    subroutine adaptiveSolveDivAmpereAnderson(solver, particleList, world, del_t, remainDel_t, currDel_t, timeCurrent)
        ! Solve for divergence of ampere's law with picard
        ! cut del_t in two if non-convergence after maxIter or is stagnating convergence, repeat until convergence
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
            ! Time step is split
            amountTimeSplits = amountTimeSplits + 1
            iterNumAdaptiveSteps = 0
            do while (iterNumPicard == maxIter)
                currDel_t = currDel_t/2.0d0
                iterNumAdaptiveSteps = iterNumAdaptiveSteps + 1
                if (iterNumAdaptiveSteps > 4) then
                    ! If you have to reduce by 4x, likely something wrong
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
        ! General potential solver
        type(potentialSolver), intent(in out) :: solver
        type(Particle), intent(in out) :: particleList(numberChargedParticles)
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t, timeCurrent
        real(real64), intent(in out) :: remainDel_t, currDel_t
        ! make future phi now current phi
        ! allows use of phi_f and phi for diagnostics before solving next time step
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
        ! Function evaluation for Nitsol mod
        ! Instead of reorganizing all data into rpar, use globals set in this module for simplicity
        integer, intent(in) :: n
        integer, intent(in out) :: itrmf, ipar(*)
        double precision, intent(in) :: xcur(n)
        double precision, intent(in out) :: rpar(*), fcur(n)
        integer(int64) :: endTimer, startTimer
        globalSolver%phi_f = xcur
        call system_clock(startTimer)
        call depositJ(globalSolver, globalParticleList, globalWorld, rpar(1))
        call system_clock(endTimer)
        moverTimer = moverTimer + endTimer - startTimer
        call system_clock(startTimer)
        call globalSolver%solve_tridiag_Ampere(globalWorld, rpar(1))
        call system_clock(endTimer)
        potTimer = potTimer + endTimer - startTimer
        fcur = xcur - globalSolver%phi_f
        itrmf = 0
        
    end subroutine funcNitsol

    subroutine jacNitsol(n, xcur, fcur, ijob, v, z, rpar, ipar, itrmjv) 
        ! If analytical jacobian matrix-vector product or right-preconditioner needed
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
        ! Solve divergence of ampere with Nitsol-mod
        real(real64), intent(in) :: del_t
        integer(int32) :: info(6), iterm, ipar(2), itrmf
        real(real64) :: fcurSolver(NumberXNodes), xcurSolver(NumberXNodes), rpar(1)
        integer(int64) :: startTimer, endTimer
        iterm = 0
        call globalSolver%makeConstSourceTerm(globalWorld)
        xcurSolver = globalSolver%phi_f
        rpar(1) = del_t
        call nitsol(NumberXNodes, xcurSolver, funcNitsol, jacNitsol, eps_a, eps_r, 1.d-12, info, rpar, ipar, iterm)
        SELECT CASE (iterm)
        CASE(0)
            ! Move final particles
            iterNumPicard = info(1)
            call system_clock(startTimer)
            call moveParticles(globalSolver, globalParticleList, globalWorld, del_t)
            call system_clock(endTimer)
            moverTimer = moverTimer + endTimer - startTimer
        CASE(1,5,6)
            ! Issue with convergence, set to maxIter for adaptive time step
            iterNumPicard = maxIter
        CASE default
            print *, "Nitsol error with iterm == ", iterm
            stop
        END SELECT
    end subroutine solveJFNK

    subroutine adaptiveSolveDivAmpereJFNK(del_t, remainDel_t, currDel_t, timeCurrent)
        ! Same function as AA applied to JFNK
        ! JFNK does not have stagnation check so far
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
        ! Write solver state 
        character(*), intent(in) :: dirName
        open(15,file=dirName//'/SolverState.dat')
        write(15,'("Solver Type, eps_a, eps_r, m_Anderson, Beta_k, maximum iterations")')
        write(15,"(1(I3.3, 1x), 2(es16.8,1x), 1(I3.3, 1x), 1(es16.8,1x), 1(I3.3, 1x))") solverType, eps_a, eps_r, m_Anderson, Beta_k, maxIter
        close(15)

    end subroutine writeSolverState

end module mod_nonLinSolvers