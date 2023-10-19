module mod_particleMover
    use iso_fortran_env, only: int32, real64, output_unit
    use constants
    use mod_BasicFunctions
    use mod_particle
    use mod_domain
    use mod_potentialSolver
    use mod_particleInjection
    use omp_lib
    implicit none
    ! Procedures for moving particles and depositing J 
    integer(int32), parameter :: m_Anderson_Particle = 2

contains

    subroutine subStepSolverAA(l_sub, v_sub, l_f, v_f, v_half, del_tau, timePassed, del_t, E_x, q_over_m, l_cell, BField, B_mag, dx_dl, f_tol, AtBoundaryBool)
        real(real64), intent(in) :: l_sub, v_sub(3), del_t, E_x, q_over_m, BField(3), B_mag, dx_dl, f_tol
        integer(int32), intent(in) :: l_cell
        real(real64), intent(in out) :: l_f, v_f(3), timePassed, del_tau, v_half(3)
        logical, intent(in out) :: AtBoundaryBool
        logical :: convergeBool = .false.
        real(real64) :: Res_k(m_Anderson_Particle+1), del_tau_k(m_Anderson_Particle+1), fitMat(m_Anderson_Particle), minCoeff(m_Anderson_Particle+1), del_tau_max
        integer(int32) :: m_k, index, k, l_boundary, u
        real(real64) :: coeffAccel, v_prime(3), curr_Res
        ! print *, ''
        ! print *, 'ENTER SUBSTEP SOLVER'
        ! print *, "--------------"
        ! print *, 'initial v_sub:', v_sub
        ! print *, 'initial l_sub:', l_sub
        ! print *, 'AtBoundaryBool:', AtBoundaryBool
        ! print *, 'l_cell:', l_cell
        ! print *, 'E_x:', E_x
        v_half = v_sub
        v_prime = v_sub
        del_tau_max = del_t - timePassed
        del_tau_k(1) = del_tau_max
        k = 0
        convergeBool = .false.
        do while (.not. convergeBool)
            index = MODULO(k, m_Anderson_Particle+1) + 1
            m_k = MIN(k, m_Anderson_Particle)
            coeffAccel = 0.5d0 * del_tau_k(index) * q_over_m
            v_prime(1) = v_sub(1) + coeffAccel * E_x
            v_half = v_prime + coeffAccel * (crossProduct(v_prime, BField) + coeffAccel* SUM(v_prime * BField) * BField)
            v_half = v_half / (1.0d0 + (coeffAccel*B_mag)**2)
            if (v_half(1) > 0) then
                l_boundary = l_cell + 1
            else
                l_boundary = l_cell
            end if
            del_tau = MIN(del_tau_max, (real(l_boundary) - l_sub) * dx_dl / v_half(1))
            Res_k(index) = del_tau - del_tau_k(index)
            print *, 'k is:', k, 'with del_tau_k:', del_tau_k
            curr_Res = ABS(Res_k(index))
            if (curr_Res < f_tol) then
                if (del_tau < del_tau_max) then
                    AtBoundaryBool = .true.
                    l_f = real(l_boundary)
                else
                    AtBoundaryBool = .false.
                    l_f = l_sub + v_half(1) * del_tau / dx_dl
                end if
                v_f = 2.0d0 * v_half - v_sub
                timePassed = timePassed + del_tau
                convergeBool = .true.
                exit
            end if
            if (k > m_Anderson_Particle) then
                if (curr_Res >= SUM(abs(Res_k))/real(m_Anderson_Particle+1)) then
                    del_tau_max = 0.5d0 * del_tau_max
                    del_tau_k(1) = del_tau_max
                    k = -1
                end if
            end if
            if (k > 0) then
                do u = 0, m_k-1
                    fitMat(u+1) = Res_k(MODULO(k - m_k + u, m_Anderson_Particle+1) + 1) - Res_k(index)
                end do
                minCoeff(1:m_k) = (fitMat(1:m_k) * -Res_k(index))/SUM(fitMat(1:m_k)**2)
                minCoeff(m_k+1) = 1.0d0 - SUM(minCoeff(1:m_k))
                del_tau_k(MODULO(k+1, m_Anderson_Particle+1) + 1) = minCoeff(1) * (Res_k(MODULO(k-m_k, m_Anderson_Particle+1) + 1) + del_tau_k(MODULO(k-m_k,m_Anderson_Particle+1) + 1))
                do u=1, m_k
                    del_tau_k(MODULO(k+1, m_Anderson_Particle+1) + 1) = del_tau_k(MODULO(k+1, m_Anderson_Particle+1) + 1) + minCoeff(u + 1) * (Res_k(MODULO(k-m_k + u, m_Anderson_Particle+1) + 1) + del_tau_k(MODULO(k-m_k + u, m_Anderson_Particle+1) + 1))
                end do
            else
                del_tau_k(2) = del_tau
            end if
            k = k + 1
        end do
        ! print *, 'v_f is:', v_f
        ! print *, 'l_f is:', l_f
        ! print *, 'k is:', k
        ! print *, ''

    end subroutine subStepSolverAA


    subroutine depositJ(solver, particleList, world, del_t)
        ! particle substepping procedure which deposits J, also used for general moving of particles
        ! boolDepositJ = true : deposit J
        ! boolDepositJ = false: only move particles (done after non-linear iteration), also delete particles for wall collisions (saves extra loop in collisions)
        class(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        type(Particle), intent(in out) :: particleList(:)
        real(real64), intent(in) :: del_t
        !a and c correspond to quadratic equations | l_alongV is nearest integer boundary along velocity component, away is opposite
        real(real64) :: l_f, l_sub, v_sub(3), v_f(3), timePassed, del_tau, q_over_m, f_tol, v_half(3), dx_dl, E_x
        integer(int32) :: j, i, l_cell, iThread
        logical :: AtBoundaryBool
        logical :: convergeBool = .false.
        real(real64) :: Res_k(m_Anderson_Particle+1), del_tau_k(m_Anderson_Particle+1), fitMat(m_Anderson_Particle), minCoeff(m_Anderson_Particle+1), del_tau_max
        integer(int32) :: m_k, index, k, l_boundary, u
        real(real64) :: coeffAccel, v_prime(3), curr_Res
        solver%J = 0.0d0
        solver%EField = 0.5d0 * (solver%phi(1:NumberXNodes-1) + solver%phi_f(1:NumberXNodes-1) - solver%phi(2:NumberXNodes) - solver%phi_f(2:NumberXNodes)) / world%dx_dl
        f_tol = del_t * 1.d-8
        loopSpecies: do j = 1, numberChargedParticles
            q_over_m = particleList(j)%q/particleList(j)%mass
            !$OMP parallel private(iThread, i, l_f, l_sub, v_sub, v_f, v_half, timePassed, del_tau, l_cell, AtBoundaryBool, dx_dl, E_x, &
            convergeBool, Res_k, del_tau_k, fitMat, minCoeff, del_tau_max, m_k, index, k, l_boundary, u, coeffAccel, v_prime, curr_Res)
            iThread = omp_get_thread_num() + 1 
            loopParticles: do i = 1, particleList(j)%N_p(iThread)
                v_sub = particleList(j)%phaseSpace(2:4,i,iThread)
                l_sub = particleList(j)%phaseSpace(1,i,iThread)
                timePassed = 0.0d0
                AtBoundaryBool = (MOD(l_sub, 1.0d0) == 0)
                do while((timePassed < del_t))
                    if (.not. AtBoundaryBool) then
                        l_cell = INT(l_sub)
                    else
                        l_cell = INT(l_sub + SIGN(0.5d0, v_sub(1)))
                    end if
                    E_x = solver%EField(l_cell)
                    dx_dl = world%dx_dl(l_cell)
                    v_half = v_sub
                    v_prime = v_sub
                    del_tau_max = del_t - timePassed
                    del_tau_k(1) = del_tau_max
                    k = 0
                    convergeBool = .false.
                    do while (.not. convergeBool)
                        index = MODULO(k, m_Anderson_Particle+1) + 1
                        m_k = MIN(k, m_Anderson_Particle)
                        coeffAccel = 0.5d0 * del_tau_k(index) * q_over_m
                        v_prime(1) = v_sub(1) + coeffAccel * E_x
                        v_half = v_prime + coeffAccel * (crossProduct(v_prime, solver%BField) + coeffAccel* SUM(v_prime * solver%BField) * solver%BField)
                        v_half = v_half / (1.0d0 + (coeffAccel*solver%BFieldMag)**2)
                        if (v_half(1) > 0) then
                            l_boundary = l_cell + 1
                        else
                            l_boundary = l_cell
                        end if
                        del_tau = MIN(del_tau_max, (real(l_boundary) - l_sub) * dx_dl / v_half(1))
                        Res_k(index) = del_tau - del_tau_k(index)
                        curr_Res = ABS(Res_k(index))
                        if (curr_Res < f_tol) then
                            if (del_tau <= del_tau_max) then
                                AtBoundaryBool = .true.
                                l_f = real(l_boundary)
                            else
                                AtBoundaryBool = .false.
                                l_f = l_sub + v_half(1) * del_tau / dx_dl
                            end if
                            v_f = 2.0d0 * v_half - v_sub
                            timePassed = timePassed + del_tau
                            convergeBool = .true.
                            exit
                        end if
                        if (k > m_Anderson_Particle) then
                            if (curr_Res >= SUM(abs(Res_k))/real(m_Anderson_Particle+1)) then
                                del_tau_max = 0.5d0 * del_tau_max
                                del_tau_k(1) = del_tau_max
                                k = -1
                            end if
                        end if
                        if (k > 0) then
                            do u = 0, m_k-1
                                fitMat(u+1) = Res_k(MODULO(k - m_k + u, m_Anderson_Particle+1) + 1) - Res_k(index)
                            end do
                            minCoeff(1:m_k) = (fitMat(1:m_k) * -Res_k(index))/SUM(fitMat(1:m_k)**2)
                            minCoeff(m_k+1) = 1.0d0 - SUM(minCoeff(1:m_k))
                            del_tau_k(MODULO(k+1, m_Anderson_Particle+1) + 1) = minCoeff(1) * (Res_k(MODULO(k-m_k, m_Anderson_Particle+1) + 1) + del_tau_k(MODULO(k-m_k,m_Anderson_Particle+1) + 1))
                            do u=1, m_k
                                del_tau_k(MODULO(k+1, m_Anderson_Particle+1) + 1) = del_tau_k(MODULO(k+1, m_Anderson_Particle+1) + 1) + minCoeff(u + 1) * (Res_k(MODULO(k-m_k + u, m_Anderson_Particle+1) + 1) + del_tau_k(MODULO(k-m_k + u, m_Anderson_Particle+1) + 1))
                            end do
                        else
                            del_tau_k(2) = del_tau
                        end if
                        k = k + 1
                    end do
                    solver%J(l_cell, iThread) = solver%J(l_cell, iThread) + particleList(j)%w_p * particleList(j)%q * (v_half(1))*del_tau/world%dx_dl(l_cell)/del_t
                    
                    if (AtBoundaryBool) then
                        SELECT CASE (world%boundaryConditions(INT(l_f)))
                        CASE(0)
                            continue
                        CASE(1,4)
                            exit
                        CASE(2)
                            if (INT(l_f) == NumberXNodes) then
                                v_f(1) = -ABS(v_f(1))
                            else
                                v_f(1) = ABS(v_f(1))
                            end if
                            if (v_f(1) > 0.0d0) then
                                print *, 'Issue with v_f after neumann'
                                stop
                            end if
                        CASE(3)
                            l_f = ABS(l_f - real(NumberXNodes, kind = real64) - 1.0d0)
                        CASE default
                            print *, "l_sub is:", l_sub
                            print *, 'l_f is:', l_f
                            print *, world%boundaryConditions(INT(l_f))
                            print *, "Case does not exist in ongoing substep, depositJ"
                            stop
                        END SELECT
                    end if
                    ! now final position/velocity becomes next starting position/velocity
                    l_sub = l_f
                    v_sub = v_f
                    if ((l_f < 1) .or. (l_f > NumberXNodes)) then
                        stop "Have particles travelling outside the domain!"
                    end if
                end do
            end do loopParticles
            !$OMP end parallel
        end do loopSpecies
        
    end subroutine depositJ


    ! -------------------------------------------- Particle mover without boolean checks for depositing J ------------------------------------------------------------

    subroutine moveParticlesTest(solver, particleList, world, del_t)
        ! particle mover to avoid the boolean checks which mostly don't happen when depositing J
        class(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        type(Particle), intent(in out) :: particleList(:)
        real(real64), intent(in) :: del_t
        !a and c correspond to quadratic equations | l_alongV is nearest integer boundary along velocity component, away is opposite
        real(real64) :: l_f, l_sub, v_sub(3), v_f(3), timePassed, del_tau, q_over_m, f_tol, v_half(3), dx_dl, E_x
        integer(int32) :: j, i, l_cell, iThread, delIdx
        logical :: AtBoundaryBool
        solver%EField = 0.5d0 * (solver%phi(1:NumberXNodes-1) + solver%phi_f(1:NumberXNodes-1) - solver%phi(2:NumberXNodes) - solver%phi_f(2:NumberXNodes)) / world%dx_dl
        f_tol = del_t * 1.d-8
        loopSpecies: do j = 1, numberChargedParticles
            q_over_m = particleList(j)%q/particleList(j)%mass
            !$OMP parallel private(iThread, i, l_f, l_sub, v_sub, v_f, v_half, timePassed, del_tau, l_cell, AtBoundaryBool, delIdx, dx_dl, E_x)
            iThread = omp_get_thread_num() + 1 
            delIdx = 0
            particleList(j)%refIdx(iThread) = 0
            particleList(j)%energyLoss(:, iThread) = 0.0d0
            particleList(j)%wallLoss(:, iThread) = 0.0d0
            loopParticles: do i = 1, particleList(j)%N_p(iThread)
                v_sub = particleList(j)%phaseSpace(2:4,i,iThread)
                l_sub = particleList(j)%phaseSpace(1,i,iThread)
                timePassed = 0.0d0
                AtBoundaryBool = (MOD(l_sub, 1.0d0) == 0)
                do while((timePassed < del_t))
                    if (.not. AtBoundaryBool) then
                        l_cell = INT(l_sub)
                    else
                        l_cell = INT(l_sub + SIGN(0.5d0, v_sub(1)))
                    end if
                    E_x = solver%EField(l_cell)
                    dx_dl = world%dx_dl(l_cell)
                    call subStepSolverAA(l_sub, v_sub, l_f, v_f, v_half, del_tau, timePassed, del_t, E_x, q_over_m, &
                        l_cell, solver%BField, solver%BFieldMag, dx_dl, f_tol, AtBoundaryBool)
                    if (AtBoundaryBool) then
                        SELECT CASE (world%boundaryConditions(INT(l_f)))
                        CASE(0)
                            continue
                        CASE(1)
                            delIdx = delIdx + 1
                            if (l_f == 1) then
                                particleList(j)%energyLoss(1, iThread) = particleList(j)%energyLoss(1, iThread) + (SUM(v_f**2))!J/m^2 in 1D
                                particleList(j)%wallLoss(1, iThread) = particleList(j)%wallLoss(1, iThread) + 1 !C/m^2 in 1D
                            else if (l_f == NumberXNodes) then
                                particleList(j)%energyLoss(2, iThread) = particleList(j)%energyLoss(2, iThread) + (SUM(v_f**2)) !J/m^2 in 1D
                                particleList(j)%wallLoss(2, iThread) = particleList(j)%wallLoss(2, iThread) + 1 !C/m^2 in 1D
                            end if
                            exit
                        CASE(2)
                            if (INT(l_f) == NumberXNodes) then
                                v_f(1) = -ABS(v_f(1))
                            else
                                v_f(1) = ABS(v_f(1))
                            end if
                            particleList(j)%refIdx(iThread) = particleList(j)%refIdx(iThread) + 1
                            particleList(j)%refRecordIdx(particleList(j)%refIdx(iThread), iThread) = i - delIdx
                        CASE(3)
                            l_f = ABS(l_f - real(NumberXNodes, kind = real64) - 1.0d0)
                        CASE(4)
                            if (refluxPartBool .or. injectionBool) then
                                particleList(j)%phaseSpace(1, i-delIdx, iThread) = l_f
                                particleList(j)%phaseSpace(2:4, i-delIdx, iThread) = v_f
                                particleList(j)%refIdx(iThread) = particleList(j)%refIdx(iThread) + 1
                                particleList(j)%refRecordIdx(particleList(j)%refIdx(iThread), iThread) = i - delIdx
                            else
                                delIdx = delIdx + 1
                                if (l_f == 1) then
                                    particleList(j)%energyLoss(1, iThread) = particleList(j)%energyLoss(1, iThread) + (SUM(v_f**2))!J/m^2 in 1D
                                    particleList(j)%wallLoss(1, iThread) = particleList(j)%wallLoss(1, iThread) + 1 !C/m^2 in 1D
                                else if (l_f == NumberXNodes) then
                                    particleList(j)%energyLoss(2, iThread) = particleList(j)%energyLoss(2, iThread) + (SUM(v_f**2)) !J/m^2 in 1D
                                    particleList(j)%wallLoss(2, iThread) = particleList(j)%wallLoss(2, iThread) + 1 !C/m^2 in 1D
                                end if
                                exit
                            end if
                        CASE default
                            print *, "l_sub is:", l_sub
                            print *, 'l_f is:', l_f
                            print *, world%boundaryConditions(INT(l_f))
                            print *, "Case does not exist in ongoing substep, depositJ"
                            stop
                        END SELECT
                    else
                        particleList(j)%phaseSpace(1, i-delIdx, iThread) = l_f
                        particleList(j)%phaseSpace(2:4,i-delIdx, iThread) = v_f
                    end if
                    ! now final position/velocity becomes next starting position/velocity
                    l_sub = l_f
                    v_sub = v_f
                end do
                stop
                if ((l_f < 1) .or. (l_f > NumberXNodes)) then
                    stop "Have particles travelling outside domain!"
                end if
                
            end do loopParticles
            particleList(j)%N_p(iThread) = particleList(j)%N_p(iThread) - delIdx
            particleList(j)%delIdx(iThread) = delIdx
            !$OMP end parallel
            particleList(j)%accumEnergyLoss = particleList(j)%accumEnergyLoss + SUM(particleList(j)%energyLoss, DIM = 2)
            particleList(j)%accumWallLoss = particleList(j)%accumWallLoss + SUM(particleList(j)%wallLoss, DIM = 2)
        end do loopSpecies
    end subroutine moveParticlesTest

    subroutine moveParticles(solver, particleList, world, del_t)
        ! particle mover to avoid the boolean checks which mostly don't happen when depositing J
        class(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        type(Particle), intent(in out) :: particleList(:)
        real(real64), intent(in) :: del_t
        !a and c correspond to quadratic equations | l_alongV is nearest integer boundary along velocity component, away is opposite
        real(real64) :: l_f, l_sub, v_sub(3), v_f(3), timePassed, del_tau, q_over_m, f_tol, v_half(3), dx_dl, E_x
        integer(int32) :: j, i, l_cell, iThread, delIdx
        logical :: AtBoundaryBool
        logical :: convergeBool = .false.
        real(real64) :: Res_k(m_Anderson_Particle+1), del_tau_k(m_Anderson_Particle+1), fitMat(m_Anderson_Particle), minCoeff(m_Anderson_Particle+1), del_tau_max
        integer(int32) :: m_k, index, k, l_boundary, u
        real(real64) :: coeffAccel, v_prime(3), curr_Res
        solver%EField = 0.5d0 * (solver%phi(1:NumberXNodes-1) + solver%phi_f(1:NumberXNodes-1) - solver%phi(2:NumberXNodes) - solver%phi_f(2:NumberXNodes)) / world%dx_dl
        f_tol = del_t * 1.d-8
        loopSpecies: do j = 1, numberChargedParticles
            q_over_m = particleList(j)%q/particleList(j)%mass
            !$OMP parallel private(iThread, i, l_f, l_sub, v_sub, v_f, v_half, timePassed, del_tau, l_cell, AtBoundaryBool, delIdx, dx_dl, E_x, &
            convergeBool, Res_k, del_tau_k, fitMat, minCoeff, del_tau_max, m_k, index, k, l_boundary, u, coeffAccel, v_prime, curr_Res)
            iThread = omp_get_thread_num() + 1 
            delIdx = 0
            particleList(j)%refIdx(iThread) = 0
            particleList(j)%energyLoss(:, iThread) = 0.0d0
            particleList(j)%wallLoss(:, iThread) = 0.0d0
            loopParticles: do i = 1, particleList(j)%N_p(iThread)
                v_sub = particleList(j)%phaseSpace(2:4,i,iThread)
                l_sub = particleList(j)%phaseSpace(1,i,iThread)
                timePassed = 0.0d0
                AtBoundaryBool = (MOD(l_sub, 1.0d0) == 0)
                do while((timePassed < del_t))
                    if (.not. AtBoundaryBool) then
                        l_cell = INT(l_sub)
                    else
                        l_cell = INT(l_sub + SIGN(0.5d0, v_sub(1)))
                    end if

                    ! Anderson-Accelerated Picard
                    E_x = solver%EField(l_cell)
                    dx_dl = world%dx_dl(l_cell)
                    v_half = v_sub
                    v_prime = v_sub
                    del_tau_max = del_t - timePassed
                    del_tau_k(1) = del_tau_max
                    k = 0
                    convergeBool = .false.
                    do while (.not. convergeBool)
                        index = MODULO(k, m_Anderson_Particle+1) + 1
                        m_k = MIN(k, m_Anderson_Particle)
                        coeffAccel = 0.5d0 * del_tau_k(index) * q_over_m
                        v_prime(1) = v_sub(1) + coeffAccel * E_x
                        v_half = v_prime + coeffAccel * (crossProduct(v_prime, solver%BField) + coeffAccel* SUM(v_prime * solver%BField) * solver%BField)
                        v_half = v_half / (1.0d0 + (coeffAccel*solver%BFieldMag)**2)
                        if (v_half(1) > 0) then
                            l_boundary = l_cell + 1
                        else
                            l_boundary = l_cell
                        end if
                        del_tau = MIN(del_tau_max, (real(l_boundary) - l_sub) * dx_dl / v_half(1))
                        Res_k(index) = del_tau - del_tau_k(index)
                        curr_Res = ABS(Res_k(index))
                        if (curr_Res < f_tol) then
                            if (del_tau <= del_tau_max) then
                                AtBoundaryBool = .true.
                                l_f = real(l_boundary)
                            else
                                AtBoundaryBool = .false.
                                l_f = l_sub + v_half(1) * del_tau / dx_dl
                            end if
                            v_f = 2.0d0 * v_half - v_sub
                            timePassed = timePassed + del_tau
                            convergeBool = .true.
                            exit
                        end if
                        if (k > m_Anderson_Particle) then
                            if (curr_Res >= SUM(abs(Res_k))/real(m_Anderson_Particle+1)) then
                                del_tau_max = 0.5d0 * del_tau_max
                                del_tau_k(1) = del_tau_max
                                k = -1
                            end if
                        end if
                        if (k > 0) then
                            do u = 0, m_k-1
                                fitMat(u+1) = Res_k(MODULO(k - m_k + u, m_Anderson_Particle+1) + 1) - Res_k(index)
                            end do
                            minCoeff(1:m_k) = (fitMat(1:m_k) * -Res_k(index))/SUM(fitMat(1:m_k)**2)
                            minCoeff(m_k+1) = 1.0d0 - SUM(minCoeff(1:m_k))
                            del_tau_k(MODULO(k+1, m_Anderson_Particle+1) + 1) = minCoeff(1) * (Res_k(MODULO(k-m_k, m_Anderson_Particle+1) + 1) + del_tau_k(MODULO(k-m_k,m_Anderson_Particle+1) + 1))
                            do u=1, m_k
                                del_tau_k(MODULO(k+1, m_Anderson_Particle+1) + 1) = del_tau_k(MODULO(k+1, m_Anderson_Particle+1) + 1) + minCoeff(u + 1) * (Res_k(MODULO(k-m_k + u, m_Anderson_Particle+1) + 1) + del_tau_k(MODULO(k-m_k + u, m_Anderson_Particle+1) + 1))
                            end do
                        else
                            del_tau_k(2) = del_tau
                        end if
                        k = k + 1
                    end do
                    ! End AA-Picard

                    if (AtBoundaryBool) then
                        SELECT CASE (world%boundaryConditions(INT(l_f)))
                        CASE(0)
                            continue
                        CASE(1)
                            delIdx = delIdx + 1
                            if (l_f == 1) then
                                particleList(j)%energyLoss(1, iThread) = particleList(j)%energyLoss(1, iThread) + (SUM(v_f**2))!J/m^2 in 1D
                                particleList(j)%wallLoss(1, iThread) = particleList(j)%wallLoss(1, iThread) + 1 !C/m^2 in 1D
                            else if (l_f == NumberXNodes) then
                                particleList(j)%energyLoss(2, iThread) = particleList(j)%energyLoss(2, iThread) + (SUM(v_f**2)) !J/m^2 in 1D
                                particleList(j)%wallLoss(2, iThread) = particleList(j)%wallLoss(2, iThread) + 1 !C/m^2 in 1D
                            end if
                            exit
                        CASE(2)
                            if (INT(l_f) == NumberXNodes) then
                                v_f(1) = -ABS(v_f(1))
                            else
                                v_f(1) = ABS(v_f(1))
                            end if
                            particleList(j)%refIdx(iThread) = particleList(j)%refIdx(iThread) + 1
                            particleList(j)%refRecordIdx(particleList(j)%refIdx(iThread), iThread) = i - delIdx
                        CASE(3)
                            l_f = ABS(l_f - real(NumberXNodes, kind = real64) - 1.0d0)
                        CASE(4)
                            if (refluxPartBool .or. injectionBool) then
                                particleList(j)%phaseSpace(1, i-delIdx, iThread) = l_f
                                particleList(j)%phaseSpace(2:4, i-delIdx, iThread) = v_f
                                particleList(j)%refIdx(iThread) = particleList(j)%refIdx(iThread) + 1
                                particleList(j)%refRecordIdx(particleList(j)%refIdx(iThread), iThread) = i - delIdx
                            else
                                delIdx = delIdx + 1
                                if (l_f == 1) then
                                    particleList(j)%energyLoss(1, iThread) = particleList(j)%energyLoss(1, iThread) + (SUM(v_f**2))!J/m^2 in 1D
                                    particleList(j)%wallLoss(1, iThread) = particleList(j)%wallLoss(1, iThread) + 1 !C/m^2 in 1D
                                else if (l_f == NumberXNodes) then
                                    particleList(j)%energyLoss(2, iThread) = particleList(j)%energyLoss(2, iThread) + (SUM(v_f**2)) !J/m^2 in 1D
                                    particleList(j)%wallLoss(2, iThread) = particleList(j)%wallLoss(2, iThread) + 1 !C/m^2 in 1D
                                end if
                                exit
                            end if
                        CASE default
                            print *, "l_sub is:", l_sub
                            print *, 'l_f is:', l_f
                            print *, world%boundaryConditions(INT(l_f))
                            print *, "Case does not exist in ongoing substep, depositJ"
                            stop
                        END SELECT
                    else
                        particleList(j)%phaseSpace(1, i-delIdx, iThread) = l_f
                        particleList(j)%phaseSpace(2:4,i-delIdx, iThread) = v_f
                    end if
                    ! now final position/velocity becomes next starting position/velocity
                    l_sub = l_f
                    v_sub = v_f
                end do
                if ((l_f < 1) .or. (l_f > NumberXNodes)) then
                    stop "Have particles travelling outside domain!"
                end if
                
            end do loopParticles
            particleList(j)%N_p(iThread) = particleList(j)%N_p(iThread) - delIdx
            particleList(j)%delIdx(iThread) = delIdx
            !$OMP end parallel
            particleList(j)%accumEnergyLoss = particleList(j)%accumEnergyLoss + SUM(particleList(j)%energyLoss, DIM = 2)
            particleList(j)%accumWallLoss = particleList(j)%accumWallLoss + SUM(particleList(j)%wallLoss, DIM = 2)
        end do loopSpecies
    end subroutine moveParticles


end module mod_particleMover