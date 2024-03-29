module mod_particleMover
    use iso_fortran_env, only: int32, real64, output_unit
    use constants
    use mod_BasicFunctions
    use mod_particle
    use mod_domain
    use mod_potentialSolver
    use omp_lib
    implicit none
    integer(int32), parameter :: m_Anderson_Particle = 2


contains

    subroutine subStepSolverPI(l_sub, v_sub, l_f, v_f, v_half, del_tau, timePassed, d_half, q_over_m, l_cell, E_left, E_right, BField, B_mag, dx_dl, f_tol, AtBoundaryBool, &
        del_t, l_boundary)
        ! picard iteration particle mover, like Chen in 2015 paper
        integer(int32), intent(in) :: l_cell
        real(real64), intent(in) :: q_over_m, BField(3), B_mag, dx_dl, f_tol, del_t, E_left, E_right
        real(real64), intent(in out) :: l_sub, v_sub(3), l_f, v_f(3), v_half(3), del_tau, timePassed, d_half
        logical, intent(in out) :: AtBoundaryBool
        integer(int32), intent(in out) :: l_boundary
        real(real64) :: v_prime(3), coeffAccel, curr_Res, l_f_prev, E_x, del_tau_prev
        integer(int32) :: k
        logical :: convergeBool
        
        v_prime = v_sub
        v_half = v_sub
        ! del_tau_prev being smaller (like fraction gyro-period) makes things way slower, more substeps
        del_tau_prev = del_t - timePassed
        l_f_prev = l_sub
        convergeBool = .false.
        k = 0
        do while (.not. convergeBool)
            ! First try full time, see where it ends up
            coeffAccel = 0.5d0 * del_tau_prev * q_over_m
            d_half = (l_sub + l_f_prev) * 0.5d0 - real(l_cell)
            E_x = E_right * d_half + E_left * (1.0d0 - d_half)
            v_prime(1) = v_sub(1) + coeffAccel * E_x
            v_half = v_prime + coeffAccel * (crossProduct(v_prime, BField) + coeffAccel* SUM(v_prime * BField) * BField)
            v_half = v_half / (1.0d0 + (coeffAccel*B_mag)**2)
            if (v_half(1) > 0) then
                l_boundary = l_cell + 1
            else
                l_boundary = l_cell
            end if
            del_tau = (real(l_boundary, kind = real64) - l_sub) * dx_dl / v_half(1)
            del_tau = MIN(del_tau, del_tau_prev)
            ! print *, 'del_tau in iteration is:', del_tau
            l_f = l_sub + v_half(1) * del_tau / dx_dl
            ! print *, 'l_f in iteration is:', l_f
            
            if (ABS(l_f_prev - l_f) < 1.d-8 .and. ABS(del_tau_prev - del_tau) < f_tol) then
                v_f = 2.0d0 * v_half - v_sub
                AtBoundaryBool = (ABS(l_f - real(l_boundary)) < 1.d-10)
                if (.not. AtBoundaryBool .and. INT(l_f) /= l_cell) then
                    print *, "not at boundary but l_f not in cell"
                    stop
                end if
                if (del_tau <= 0.0d0) then
                    print *, 'AtBoundaryBool:', AtBoundaryBool
                    print *, 'del_tau is 0 or less'
                    print *, 'del_tau is:', del_tau
                    print *, 'l_sub:', l_sub
                    print *, 'l_f:', l_f
                    print *, 'v_sub:', v_sub
                    print *, 'v_half:', v_half
                    print *, 'v_f:',v_f
                    stop
                end if
                exit
            end if
            ! if (k > 30) then
            !     print *, 'PI iteration is:', k
            !     print *, 'l_sub is:', l_sub
            !     print *, 'v_sub is:', v_sub
            !     print *, 'v_half is:', v_half
            !     print *, 'l_f is:', l_f
            !     print *, 'l_f prev is:', l_f_prev
            !     print *, 'del_tau_prev is:', del_tau_prev
            !     print *, 'del_tau is:', del_tau
            !     l_f_prev = l_f
            !     del_tau_prev = del_tau
            !     coeffAccel = 0.5d0 * del_tau_prev * q_over_m
            !     d_half = (l_sub + l_f_prev) * 0.5d0 - real(l_cell)
            !     E_x = E_right * d_half + E_left * (1.0d0 - d_half)
            !     v_prime(1) = v_sub(1) + coeffAccel * E_x
            !     v_half = v_prime + coeffAccel * (crossProduct(v_prime, BField) + coeffAccel* SUM(v_prime * BField) * BField)
            !     v_half = v_half / (1.0d0 + (coeffAccel*B_mag)**2)
            !     if (v_half(1) > 0) then
            !         l_boundary = l_cell + 1
            !     else
            !         l_boundary = l_cell
            !     end if
            !     del_tau = (real(l_boundary, kind = real64) - l_sub) * dx_dl / v_half(1)
            !     del_tau = MIN(del_tau, del_tau_prev)
            !     print *, 'next del_tau is:', del_tau
            !     stop
            ! end if
            del_tau_prev = del_tau
            l_f_prev = l_f
            k = k + 1
        end do

        ! if (del_tau/del_t > 1.d-4) then
        !     if (ABS((SUM(v_f**2) - SUM(v_sub**2) - 2.0d0 * q_over_m * E_x * dx_dl * (l_f - l_sub))/(SUM(v_f**2) - SUM(v_sub**2))) > 1.d-5) then
        !         print *, 'energy issue for large enough del_tau'
        !         print *, 'l_sub is:', l_sub
        !         print *, 'v_sub is:', v_sub
        !         print *, 'del_tau_max is:', 1.0d0 - timePassed/del_t
        !         print *, 'l_f is:', l_f
        !         print *, 'v_half is:', v_half
        !         print *, 'del_tau is:', del_tau/del_t
        !         print *, 'Energy is:', ABS((SUM(v_f**2) - SUM(v_sub**2) - 2.0d0 * q_over_m * E_x * dx_dl * (l_f - l_sub))/(SUM(v_f**2) - SUM(v_sub**2))) 
        !         stop
        !     end if
        ! end if
    end subroutine subStepSolverPI

    subroutine subStepSolverAAGetPosition(l_sub, v_sub, l_f, v_half, timePassed, d_half, q_over_m, l_cell, E_left, E_right, E_x, BField, B_mag, dx_dl, AtBoundaryBool, &
        del_t, l_boundary)
        ! Anderson Acceleration particle mover
        integer(int32), intent(in) :: l_cell
        real(real64), intent(in) :: l_sub, v_sub(3), timePassed, q_over_m, BField(3), B_mag, dx_dl, del_t, E_left, E_right
        real(real64), intent(in out) :: l_f, v_half(3), d_half, E_x
        logical, intent(in out) :: AtBoundaryBool
        integer(int32), intent(in out) :: l_boundary
        real(real64) :: v_prime(3), coeffAccel, Res_k(m_Anderson_Particle+1), param_k(m_Anderson_Particle+1), fitMat(m_Anderson_Particle), del_tau_max
        real(real64) :: minCoeff(m_Anderson_Particle+1)
        integer(int32) :: k, m_k, u, index, maxIterations
        logical :: convergeBool

        maxIterations = 50
        v_prime = v_sub
        v_half = v_sub
        del_tau_max = del_t - timePassed
        ! print *, 'l_sub:', l_sub
        ! print *, 'v_sub:', v_sub
        ! print *, 'del_tau_max:', del_tau_max
        param_k(1) = l_sub
        k = 0
        convergeBool = .false.
        ! Keep del_tau const, solve for l_f
        do k = 0, maxIterations
            index = MODULO(k, m_Anderson_Particle+1) + 1
            m_k = MIN(k, m_Anderson_Particle)
            ! print *, ''
            ! print *, 'Iteration k:', k
            ! print *, 'l_f:', param_k(index)
            coeffAccel = 0.5d0 * del_tau_max * q_over_m
            d_half = (l_sub + param_k(index)) * 0.5d0 - real(l_cell)
            E_x = E_right * d_half + E_left * (1.0d0 - d_half)
            v_prime(1) = v_sub(1) + coeffAccel * E_x
            v_half = v_prime + coeffAccel * (crossProduct(v_prime, BField) + coeffAccel* SUM(v_prime * BField) * BField)
            v_half = v_half / (1.0d0 + (coeffAccel*B_mag)**2)
            l_f = l_sub + del_tau_max * v_half(1) / dx_dl
            if (v_half(1) > 0) then
                l_boundary = l_cell + 1
                AtBoundaryBool = (l_f >= real(l_boundary))
                if (AtBoundaryBool) l_f = real(l_boundary)
            else
                l_boundary = l_cell
                AtBoundaryBool = (l_f <= real(l_boundary))
                if (AtBoundaryBool) l_f = real(l_boundary)
            end if
            Res_k(index) = l_f - param_k(index)
            if (ABS(Res_k(index)) < 1.d-8) then
                exit      
            end if

            ! if (k > m_Anderson_Particle) then
            !     if (curr_Res >= SUM(abs(Res_k))/real(m_Anderson_Particle+1)) then
            !         del_tau_max = 0.5d0 * del_tau_max
            !         del_tau_k(1) = del_tau_max
            !         k = -1
            !     end if
            ! end if
            ! if (k > 5) then
            !     del_tau_max = 0.5d0 * del_tau_max
            !     deltau_lf_k(1, 1) = del_tau_max
            !     deltau_lf_k(2,1) = l_sub
            !     k = -1
            ! end if
            if (k > 0) then
                do u = 0, m_k-1
                    fitMat(u+1) = Res_k(MODULO(k - m_k + u, m_Anderson_Particle+1) + 1) - Res_k(index)
                end do
                minCoeff(1:m_k) = (fitMat(1:m_k) * -Res_k(index))/SUM(fitMat(1:m_k)**2)
                minCoeff(m_k+1) = 1.0d0 - SUM(minCoeff(1:m_k))
                param_k(MODULO(k+1, m_Anderson_Particle+1) + 1) = minCoeff(1) * (Res_k(MODULO(k-m_k, m_Anderson_Particle+1) + 1) + param_k(MODULO(k-m_k,m_Anderson_Particle+1) + 1))
                do u=1, m_k
                    param_k(MODULO(k+1, m_Anderson_Particle+1) + 1) = param_k(MODULO(k+1, m_Anderson_Particle+1) + 1) + minCoeff(u + 1) * (Res_k(MODULO(k-m_k + u, m_Anderson_Particle+1) + 1) + param_k(MODULO(k-m_k + u, m_Anderson_Particle+1) + 1))
                end do
            else
                param_k(2) = l_f
            end if
        end do

        if (k == maxIterations+1) then
            print *, 'Get Position reached maximum iterations!'
            print *, 'l_sub is:', l_sub
            print *, 'v_sub is:', v_sub
            print *, 'v_half is:', v_half
            print *, 'l_f is:', l_f
            stop
        end if

    end subroutine subStepSolverAAGetPosition

    subroutine subStepSolverAAGetTimeBoundary(l_sub, v_sub, l_f, v_half, del_tau, timePassed, q_over_m, E_x, BField, B_mag, dx_dl, f_tol, del_t)
        ! Anderson Acceleration particle mover, get time at boundary not in place
        real(real64), intent(in) :: q_over_m, BField(3), B_mag, dx_dl, f_tol, del_t, E_x
        real(real64), intent(in out) :: l_sub, v_sub(3), l_f, v_half(3), del_tau, timePassed
        real(real64) :: v_prime(3), coeffAccel, Res_k(m_Anderson_Particle+1), param_k(m_Anderson_Particle+1), fitMat(m_Anderson_Particle), del_tau_max
        real(real64) :: minCoeff(m_Anderson_Particle+1)
        integer(int32) :: k, m_k, u, index
        logical :: convergeBool

        v_prime = v_sub
        v_half = v_sub
        k = 0
        param_k(1) = del_tau
        ! l_f at boundary, solve for del_tau
        convergeBool = .false.
        do while (.not. convergeBool)
            index = MODULO(k, m_Anderson_Particle+1) + 1
            m_k = MIN(k, m_Anderson_Particle)
            coeffAccel = 0.5d0 * param_k(index) * q_over_m
            v_prime(1) = v_sub(1) + coeffAccel * E_x
            v_half = v_prime + coeffAccel * (crossProduct(v_prime, BField) + coeffAccel* SUM(v_prime * BField) * BField)
            v_half = v_half / (1.0d0 + (coeffAccel*B_mag)**2)
            del_tau = (l_f - l_sub) * dx_dl / v_half(1)
            Res_k(index) = del_tau - param_k(index)
            if (ABS(Res_k(index)) < f_tol) then
                ! if (INT(l_f) /= l_cell) then
                !     print *, 'l_f not in cell after going del_tau_max!'
                !     print *, 'del_tau is:', del_tau
                !     print *, 'l_sub:', l_sub
                !     print *, 'l_f:', l_f
                !     print *, 'v_sub:', v_sub
                !     print *, 'v_half:', v_half
                !     print *, 'v_f:',v_f
                !     print *, 'del_tau_max:', del_tau_max
                !     stop
                ! end if
                if (del_tau <= f_tol) then
                    print *, 'del_tau is less than f_tol in getTime'
                    print *, 'iteration:', k
                    print *, 'f_tol is:', f_tol
                    print *, 'del_tau is:', del_tau
                    print *, 'param_k(index):', param_k(index)
                    print *, 'l_sub:', l_sub
                    print *, 'l_f:', l_f
                    print *, 'v_sub:', v_sub
                    print *, 'v_half OG:', v_half
                    print *, 'del_tau OG:', (l_f - l_sub) * dx_dl / v_half(1)
                    stop
                end if
                convergeBool = .true.
                exit
            end if

            ! if (k > m_Anderson_Particle) then
            !     if (curr_Res >= SUM(abs(Res_k))/real(m_Anderson_Particle+1)) then
            !         del_tau_max = 0.5d0 * del_tau_max
            !         del_tau_k(1) = del_tau_max
            !         k = -1
            !     end if
            ! end if
            ! if (k > 5) then
            !     del_tau_max = 0.5d0 * del_tau_max
            !     deltau_lf_k(1, 1) = del_tau_max
            !     deltau_lf_k(2,1) = l_sub
            !     k = -1
            ! end if
            if (k > 50) then
                print *, 'AA solver taking forever in del_tau'
                print *, 'l_sub is:', l_sub
                print *, 'v_sub is:', v_sub
                print *, 'del_tau is:', del_tau
                print *, 'del_tau faction is:', del_tau/del_t
                print *, 'v_half is:', v_half
                print *, 'l_f is:', l_f
                print *, 'f_tol is:', f_tol
                print *, 'param_k:', param_k
                stop
            end if
            if (k > 0) then
                do u = 0, m_k-1
                    fitMat(u+1) = Res_k(MODULO(k - m_k + u, m_Anderson_Particle+1) + 1) - Res_k(index)
                end do
                minCoeff(1:m_k) = (fitMat(1:m_k) * -Res_k(index))/SUM(fitMat(1:m_k)**2)
                minCoeff(m_k+1) = 1.0d0 - SUM(minCoeff(1:m_k))
                param_k(MODULO(k+1, m_Anderson_Particle+1) + 1) = minCoeff(1) * (Res_k(MODULO(k-m_k, m_Anderson_Particle+1) + 1) + param_k(MODULO(k-m_k,m_Anderson_Particle+1) + 1))
                do u=1, m_k
                    param_k(MODULO(k+1, m_Anderson_Particle+1) + 1) = param_k(MODULO(k+1, m_Anderson_Particle+1) + 1) + minCoeff(u + 1) * (Res_k(MODULO(k-m_k + u, m_Anderson_Particle+1) + 1) + param_k(MODULO(k-m_k + u, m_Anderson_Particle+1) + 1))
                end do
            else
                param_k(2) = del_tau
            end if
            k = k + 1
        end do
        ! if (ABS((SUM(v_f**2) - SUM(v_sub**2) - 2.0d0 * q_over_m * E_x * dx_dl * (l_f - l_sub))/(SUM(v_f**2) - SUM(v_sub**2))) > 1.d-5) then
        !     print *, 'energy issue!'
        !     print *, 'Energy is:', ABS((SUM(v_f**2) - SUM(v_sub**2) - 2.0d0 * q_over_m * E_x * dx_dl * (l_f - l_sub))/(SUM(v_f**2) - SUM(v_sub**2))) 
        !     stop
        ! end if

    end subroutine subStepSolverAAGetTimeBoundary

    subroutine subStepSolverAATwoStepAtBoundary(l_sub, v_sub, l_f, v_f, v_half, del_tau, timePassed, d_half, q_over_m, l_cell, E_left, E_right, BField, B_mag, dx_dl, f_tol, AtBoundaryBool, &
        del_t, l_boundary)
        ! Anderson Acceleration particle mover
        integer(int32), intent(in) :: l_cell
        real(real64), intent(in) :: q_over_m, BField(3), B_mag, dx_dl, f_tol, del_t, E_left, E_right
        real(real64), intent(in out) :: l_sub, v_sub(3), l_f, v_f(3), v_half(3), del_tau, timePassed, d_half
        logical, intent(in out) :: AtBoundaryBool
        integer(int32), intent(in out) :: l_boundary
        real(real64) :: v_prime(3), coeffAccel, Res_k(m_Anderson_Particle+1), param_k(m_Anderson_Particle+1), fitMat(m_Anderson_Particle), del_tau_max, E_x
        real(real64) :: minCoeff(m_Anderson_Particle+1)
        integer(int32) :: k, m_k, u, index
        logical :: convergeBool

        v_prime = v_sub
        v_half = v_sub
        del_tau_max = del_t - timePassed
        ! print *, 'l_sub:', l_sub
        ! print *, 'v_sub:', v_sub
        ! print *, 'del_tau_max:', del_tau_max
        param_k(1) = l_sub
        k = 0
        convergeBool = .false.
        ! Keep del_tau const, solve for l_f
        do while (.not. convergeBool)
            index = MODULO(k, m_Anderson_Particle+1) + 1
            m_k = MIN(k, m_Anderson_Particle)
            ! print *, ''
            ! print *, 'Iteration k:', k
            ! print *, 'l_f:', param_k(index)
            coeffAccel = 0.5d0 * del_tau_max * q_over_m
            d_half = (l_sub + param_k(index)) * 0.5d0 - real(l_cell)
            E_x = E_right * d_half + E_left * (1.0d0 - d_half)
            v_prime(1) = v_sub(1) + coeffAccel * E_x
            v_half = v_prime + coeffAccel * (crossProduct(v_prime, BField) + coeffAccel* SUM(v_prime * BField) * BField)
            v_half = v_half / (1.0d0 + (coeffAccel*B_mag)**2)
            l_f = l_sub + del_tau_max * v_half(1) / dx_dl
            if (v_half(1) > 0) then
                l_boundary = l_cell + 1
                AtBoundaryBool = (l_f >= real(l_boundary))
                if (AtBoundaryBool) l_f = real(l_boundary)
            else
                l_boundary = l_cell
                AtBoundaryBool = (l_f <= real(l_boundary))
                if (AtBoundaryBool) l_f = real(l_boundary)
            end if
            Res_k(index) = l_f - param_k(index)
            if (ABS(Res_k(index)) < 1.d-8) then
                convergeBool = .not. AtBoundaryBool
                if (convergeBool .and. INT(l_f) /= l_cell) then
                    print *, 'l_f should be in same cell since not at boundary!'
                    stop
                end if
                exit      
            end if

            ! if (k > m_Anderson_Particle) then
            !     if (curr_Res >= SUM(abs(Res_k))/real(m_Anderson_Particle+1)) then
            !         del_tau_max = 0.5d0 * del_tau_max
            !         del_tau_k(1) = del_tau_max
            !         k = -1
            !     end if
            ! end if
            ! if (k > 5) then
            !     del_tau_max = 0.5d0 * del_tau_max
            !     deltau_lf_k(1, 1) = del_tau_max
            !     deltau_lf_k(2,1) = l_sub
            !     k = -1
            ! end if
            if (k > 50) then
                print *, 'AA solver taking forever in l_f'
                print *, 'l_sub is:', l_sub
                print *, 'v_sub is:', v_sub
                print *, 'del_tau is:', del_tau
                print *, 'del_tau faction is:', del_tau/del_t
                print *, 'v_half is:', v_half
                print *, 'l_f is:', l_f
                stop
            end if
            if (k > 0) then
                do u = 0, m_k-1
                    fitMat(u+1) = Res_k(MODULO(k - m_k + u, m_Anderson_Particle+1) + 1) - Res_k(index)
                end do
                minCoeff(1:m_k) = (fitMat(1:m_k) * -Res_k(index))/SUM(fitMat(1:m_k)**2)
                minCoeff(m_k+1) = 1.0d0 - SUM(minCoeff(1:m_k))
                param_k(MODULO(k+1, m_Anderson_Particle+1) + 1) = minCoeff(1) * (Res_k(MODULO(k-m_k, m_Anderson_Particle+1) + 1) + param_k(MODULO(k-m_k,m_Anderson_Particle+1) + 1))
                do u=1, m_k
                    param_k(MODULO(k+1, m_Anderson_Particle+1) + 1) = param_k(MODULO(k+1, m_Anderson_Particle+1) + 1) + minCoeff(u + 1) * (Res_k(MODULO(k-m_k + u, m_Anderson_Particle+1) + 1) + param_k(MODULO(k-m_k + u, m_Anderson_Particle+1) + 1))
                end do
            else
                param_k(2) = l_f
            end if
            k = k + 1
        end do
        ! print *, 'Went through first stage'
        ! print *, 'd_half is:', d_half
        ! print *, 'E_x is:', E_x
        k = 0
        param_k(1) = del_tau_max
        ! l_f at boundary, solve for del_tau
        do while (.not. convergeBool)
            index = MODULO(k, m_Anderson_Particle+1) + 1
            m_k = MIN(k, m_Anderson_Particle)
            coeffAccel = 0.5d0 * param_k(index) * q_over_m
            v_prime(1) = v_sub(1) + coeffAccel * E_x
            v_half = v_prime + coeffAccel * (crossProduct(v_prime, BField) + coeffAccel* SUM(v_prime * BField) * BField)
            v_half = v_half / (1.0d0 + (coeffAccel*B_mag)**2)
            del_tau = MAX((l_f - l_sub) * dx_dl / v_half(1), f_tol)
            Res_k(index) = del_tau - param_k(index)
            if (ABS(Res_k(index)) < f_tol) then
                ! if (INT(l_f) /= l_cell) then
                !     print *, 'l_f not in cell after going del_tau_max!'
                !     print *, 'del_tau is:', del_tau
                !     print *, 'l_sub:', l_sub
                !     print *, 'l_f:', l_f
                !     print *, 'v_sub:', v_sub
                !     print *, 'v_half:', v_half
                !     print *, 'v_f:',v_f
                !     print *, 'del_tau_max:', del_tau_max
                !     stop
                ! end if
                if (del_tau <= f_tol) then
                    print *, 'AtBoundaryBool:', AtBoundaryBool
                    print *, 'del_tau is 0 or less in two step'
                    print *, 'iteration:', k
                    print *, 'f_tol is:', f_tol
                    print *, 'del_tau is:', del_tau
                    print *, 'param_k(index):', param_k(index)
                    print *, 'l_sub:', l_sub
                    print *, 'l_f:', l_f
                    print *, 'v_sub:', v_sub
                    print *, 'v_half OG:', v_half
                    print *, 'del_tau OG:', (l_f - l_sub) * dx_dl / v_half(1)
                    print *, 'del_tau_max:', del_tau_max
                    stop
                end if
                convergeBool = .true.
                exit
            end if

            ! if (k > m_Anderson_Particle) then
            !     if (curr_Res >= SUM(abs(Res_k))/real(m_Anderson_Particle+1)) then
            !         del_tau_max = 0.5d0 * del_tau_max
            !         del_tau_k(1) = del_tau_max
            !         k = -1
            !     end if
            ! end if
            ! if (k > 5) then
            !     del_tau_max = 0.5d0 * del_tau_max
            !     deltau_lf_k(1, 1) = del_tau_max
            !     deltau_lf_k(2,1) = l_sub
            !     k = -1
            ! end if
            if (k > 50) then
                print *, 'AA solver taking forever in del_tau'
                print *, 'l_sub is:', l_sub
                print *, 'v_sub is:', v_sub
                print *, 'del_tau is:', del_tau
                print *, 'del_tau_max is:', del_tau_max
                print *, 'del_tau faction is:', del_tau/del_t
                print *, 'v_half is:', v_half
                print *, 'l_f is:', l_f
                print *, 'f_tol is:', f_tol
                print *, 'param_k:', param_k
                stop
            end if
            if (k > 0) then
                do u = 0, m_k-1
                    fitMat(u+1) = Res_k(MODULO(k - m_k + u, m_Anderson_Particle+1) + 1) - Res_k(index)
                end do
                minCoeff(1:m_k) = (fitMat(1:m_k) * -Res_k(index))/SUM(fitMat(1:m_k)**2)
                minCoeff(m_k+1) = 1.0d0 - SUM(minCoeff(1:m_k))
                param_k(MODULO(k+1, m_Anderson_Particle+1) + 1) = minCoeff(1) * (Res_k(MODULO(k-m_k, m_Anderson_Particle+1) + 1) + param_k(MODULO(k-m_k,m_Anderson_Particle+1) + 1))
                do u=1, m_k
                    param_k(MODULO(k+1, m_Anderson_Particle+1) + 1) = param_k(MODULO(k+1, m_Anderson_Particle+1) + 1) + minCoeff(u + 1) * (Res_k(MODULO(k-m_k + u, m_Anderson_Particle+1) + 1) + param_k(MODULO(k-m_k + u, m_Anderson_Particle+1) + 1))
                end do
            else
                param_k(2) = del_tau
            end if
            k = k + 1
        end do
        v_f = 2.0d0 * v_half - v_sub
        ! if (ABS((SUM(v_f**2) - SUM(v_sub**2) - 2.0d0 * q_over_m * E_x * dx_dl * (l_f - l_sub))/(SUM(v_f**2) - SUM(v_sub**2))) > 1.d-5) then
        !     print *, 'energy issue!'
        !     print *, 'Energy is:', ABS((SUM(v_f**2) - SUM(v_sub**2) - 2.0d0 * q_over_m * E_x * dx_dl * (l_f - l_sub))/(SUM(v_f**2) - SUM(v_sub**2))) 
        !     stop
        ! end if

    end subroutine subStepSolverAATwoStepAtBoundary

    
    subroutine GetRho(rho, l, world) 
        real(real64), intent(in out) :: rho(NumberXNodes)
        real(real64), intent(in) :: l
        type(Domain), intent(in) :: world
        integer(int32) :: i, j, l_center, l_left, l_right
        real(real64) :: d
        rho = 0.0d0
        l_center = INT(l)
        d = l - real(l_center)
        rho(l_center) = rho(l_center) + (-d**2 + d + 0.5d0)
        l_right = l_center + 1
        l_left = l_center - 1
        if (world%boundaryConditions(l_center) == 0 .and. world%boundaryConditions(l_right) == 0) then
            ! No Boundary either side
            rho(l_left) = rho(l_left) + 0.5d0 * (1.0d0 - d)**2
            rho(l_right) = rho(l_right) + 0.5d0 * d**2
        else if (world%boundaryConditions(l_right) == 1) then
            ! Dirichlet to right
            rho(l_center) = rho(l_center) - 0.5d0 * d**2
            rho(l_left) = rho(l_left) + 0.5d0 * (1.0d0 - d)**2
        else if (world%boundaryConditions(l_center) == 1) then
            !Dirichlet to left
            rho(l_center) = rho(l_center) - 0.5d0 * (1.0d0 - d)**2
            rho(l_right) = rho(l_right) + 0.5d0 * d**2
        else if (world%boundaryConditions(l_right) == 2) then
            !Neumann to right
            rho(l_center) = rho(l_center) + 0.5d0 * d**2
            rho(l_left) = rho(l_left) + 0.5d0 * (1.0d0 - d)**2
        else if (world%boundaryConditions(l_center) == 2) then
            !Neumann to left
            rho(l_center) = rho(l_center) + 0.5d0 * (1.0d0 - d)**2
            rho(l_right) = rho(l_right) + 0.5d0 * d**2
        end if
    end subroutine GetRho

    subroutine depositJ(solver, particleList, world, del_t)
        ! particle substepping procedure which deposits J, also used for general moving of particles
        ! boolDepositJ = true : deposit J
        ! boolDepositJ = false: only move particles (done after non-linear iteration), also delete particles for wall collisions (saves extra loop in collisions)
        class(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        type(Particle), intent(in out) :: particleList(:)
        real(real64), intent(in) :: del_t
        !a and c correspond to quadratic equations | l_alongV is nearest integer boundary along velocity component, away is opposite
        real(real64) :: l_f, l_sub, v_sub(3), v_f(3), timePassed, del_tau, q_over_m, f_tol, v_half(3), dx_dl, q_times_wp, J_part, d_half, E_x
        integer(int32) :: j, i, l_cell, iThread, l_boundary, numSubStepAve(numThread), numSubStep
        logical :: AtBoundaryBool, PrevAtBoundaryBool
        solver%J = 0.0d0
        call solver%evaluateEFieldHalfTime(world)
        f_tol = del_t * 1.d-8
        loopSpecies: do j = 1, numberChargedParticles
            q_over_m = particleList(j)%q/particleList(j)%mass
            q_times_wp = particleList(j)%q * particleList(j)%w_p
            numSubStepAve = 0
            !$OMP parallel private(iThread, i, l_f, l_sub, v_sub, v_f, v_half, timePassed, del_tau, l_cell, AtBoundaryBool, PrevAtBoundaryBool, dx_dl, l_boundary, J_part, d_half, E_x, numSubStep)
            iThread = omp_get_thread_num() + 1 
            loopParticles: do i = 1, particleList(j)%N_p(iThread)
                v_sub = particleList(j)%phaseSpace(2:4,i,iThread)
                l_sub = particleList(j)%phaseSpace(1,i,iThread)
                PrevAtBoundaryBool = .false.
                timePassed = 0.0d0
                numSubStep = 0
                l_cell = INT(l_sub)
                dx_dl = world%dx_dl(l_cell)
                call subStepSolverAAGetPosition(l_sub, v_sub, l_f, v_half, timePassed, d_half, q_over_m, l_cell, solver%EField(l_cell), solver%EField(l_cell + 1), E_x, solver%BField, &
                        solver%BFieldMag, dx_dl, AtBoundaryBool, del_t, l_boundary)
                if (.not. AtBoundaryBool) then
                    timePassed = del_t
                    v_f = 2.0d0 * v_half - v_sub
                else
                    del_tau = del_t
                    call subStepSolverAAGetTimeBoundary(l_sub, v_sub, l_f, v_half, del_tau, timePassed, q_over_m, E_x, solver%BField, solver%BFieldMag, dx_dl, f_tol, del_t)
                end if

                
                ! do while((timePassed < del_t))
                !     numSubStep = numSubStep + 1
                !     if (.not. PrevAtBoundaryBool) then
                !         l_cell = INT(l_sub)
                !     else
                !         l_cell = NINT(l_sub) + (INT(SIGN(1.0d0, v_sub)) - 1)/2
                !     end if
                !     dx_dl = world%dx_dl(l_cell)
                !     ! print *, ""
                !     ! print *, 'l_sub:', l_sub
                !     ! print *, 'v_sub:', v_sub
                !     ! Start AA
                !     call subStepSolverPI(l_sub, v_sub, l_f, v_f, v_half, del_tau, timePassed, d_half, q_over_m, l_cell, solver%EField(l_cell), solver%EField(l_cell + 1), solver%BField, &
                !         solver%BFieldMag, dx_dl, f_tol, AtBoundaryBool, del_t, l_boundary)
                !     ! print *, 'v_half is:', v_half
                !     ! print *, 'del_tau is:', del_tau
                !     ! print *, 'l_f is:', l_f
                !     ! print *, 'v_f is:', v_f
                !     J_part = q_times_wp * (v_half(1))*del_tau/dx_dl/del_t
                !     solver%J(l_cell, iThread) = solver%J(l_cell, iThread) + J_part * (1.0d0 - d_half)
                !     solver%J(l_cell+1, iThread) = solver%J(l_cell+1, iThread) + J_part * (d_half)
                !     if (AtBoundaryBool) then
                !         SELECT CASE (world%boundaryConditions(l_boundary))
                !         CASE(0)
                !             l_f = real(l_boundary) + SIGN(1.d-12, v_f(1))
                !         CASE(1,4)
                !             exit
                !         CASE(2)
                !             if (l_boundary == NumberXNodes+1) then
                !                 v_f(1) = -ABS(v_f(1))
                !                 l_f = real(l_boundary) - 1.d-12
                !             else
                !                 v_f(1) = ABS(v_f(1))
                !                 l_f = real(l_boundary) + 1.d-12
                !             end if
                !             v_f(2:3) = -v_f(2:3)  
                !         CASE(3)
                !             l_f = REAL(ABS(l_boundary - real(NumberXNodes, kind = real64) - 1)) + SIGN(1.d-12, v_f(1))
                !         CASE default
                !             print *, "l_sub is:", l_sub
                !             print *, 'l_f is:', l_f
                !             print *, world%boundaryConditions(INT(l_f))
                !             print *, "Case does not exist in ongoing substep, depositJ"
                !             stop
                !         END SELECT
                !     end if
                !     timePassed = timePassed + del_tau
                !     ! now final position/velocity becomes next starting position/velocity
                !     l_sub = l_f
                !     v_sub = v_f
                !     if ((l_f < 1) .or. (l_f > NumberXNodes+1)) then
                !         print *, 'l_f is:', l_f
                !         stop "Have particles travelling outside the domain in depositJ!"
                !     end if
                ! end do
                numSubStepAve(iThread) = numSubStepAve(iThread) + numSubStep
            end do loopParticles
            !$OMP end parallel
            print *, 'Average substeps for', particleList(j)%name, 'is:'
            print *, real(SUM(numSubStepAve))/real(SUM(particleList(j)%N_p))
        end do loopSpecies
        
    end subroutine depositJ


    ! -------------------------------------------- Particle mover without boolean checks for depositing J ------------------------------------------------------------

    subroutine moveParticles(solver, particleList, world, del_t)
        ! particle mover to avoid the boolean checks which mostly don't happen when depositing J
        class(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        type(Particle), intent(in out) :: particleList(:)
        real(real64), intent(in) :: del_t
        !a and c correspond to quadratic equations | l_alongV is nearest integer boundary along velocity component, away is opposite
        real(real64) :: l_f, l_sub, v_sub(3), v_f(3), timePassed, del_tau, q_over_m, f_tol, d_half, v_half(3), dx_dl
        integer(int32) :: j, i, l_cell, iThread, delIdx, l_boundary
        logical :: AtBoundaryBool
        call solver%evaluateEFieldHalfTime(world)
        f_tol = del_t * 1.d-8
        loopSpecies: do j = 1, numberChargedParticles
            q_over_m = particleList(j)%q/particleList(j)%mass
            !$OMP parallel private(iThread, i, l_f, l_sub, d_half, v_sub, v_f, v_half, timePassed, del_tau, l_cell, AtBoundaryBool, delIdx, dx_dl, l_boundary)
            iThread = omp_get_thread_num() + 1 
            delIdx = 0
            particleList(j)%refIdx(iThread) = 0
            particleList(j)%energyLoss(:, iThread) = 0.0d0
            particleList(j)%wallLoss(:, iThread) = 0.0d0
            loopParticles: do i = 1, particleList(j)%N_p(iThread)
                v_sub = particleList(j)%phaseSpace(2:4,i,iThread)
                l_sub = particleList(j)%phaseSpace(1,i,iThread)
                timePassed = 0.0d0
                do while((timePassed < del_t))
                    
                    l_cell = INT(l_sub)
                    dx_dl = world%dx_dl(l_cell)
                    
                    ! AA particle mover
                    call subStepSolverPI(l_sub, v_sub, l_f, v_f, v_half, del_tau, timePassed, d_half, q_over_m, l_cell, solver%EField(l_cell), solver%EField(l_cell + 1), solver%BField, &
                        solver%BFieldMag, dx_dl, f_tol, AtBoundaryBool, del_t, l_boundary)
                    
                    if (AtBoundaryBool) then
                        SELECT CASE (world%boundaryConditions(l_boundary))
                        CASE(0)
                            l_f = real(l_boundary) + SIGN(1.d-12, v_f(1))
                        CASE(1)
                            delIdx = delIdx + 1
                            if (l_boundary == 1) then
                                particleList(j)%energyLoss(1, iThread) = particleList(j)%energyLoss(1, iThread) + (SUM(v_f**2))!J/m^2 in 1D
                                particleList(j)%wallLoss(1, iThread) = particleList(j)%wallLoss(1, iThread) + 1 !C/m^2 in 1D
                            else if (l_boundary == NumberXNodes) then
                                particleList(j)%energyLoss(2, iThread) = particleList(j)%energyLoss(2, iThread) + (SUM(v_f**2)) !J/m^2 in 1D
                                particleList(j)%wallLoss(2, iThread) = particleList(j)%wallLoss(2, iThread) + 1 !C/m^2 in 1D
                            end if
                            exit
                        CASE(2)
                            if (l_boundary == NumberXNodes+1) then
                                v_f(1) = -ABS(v_f(1))
                                l_f = real(l_boundary) - 1.d-12
                            else
                                v_f(1) = ABS(v_f(1))
                                l_f = real(l_boundary) + 1.d-12
                            end if
                            v_f(2:3) = -v_f(2:3)
                            ! particleList(j)%refIdx(iThread) = particleList(j)%refIdx(iThread) + 1
                            ! particleList(j)%refRecordIdx(particleList(j)%refIdx(iThread), iThread) = i - delIdx
                        CASE(3)
                            l_f = REAL(ABS(l_boundary - real(NumberXNodes, kind = real64) - 1)) + SIGN(1.d-12, v_f(1))
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
                    if ((l_f < 1) .or. (l_f > NumberXNodes+1)) then
                        print *, "Have particles travelling outside domain in moveparticles!"
                        print *, 'del_tau is:', del_tau
                        print *, 'l_sub:', l_sub
                        print *, 'l_f:', l_f
                        print *, 'v_sub:', v_sub
                        print *, 'v_half:', v_half
                        print *, 'v_f:',v_f
                        stop
                    end if
                    timePassed = timePassed + del_tau
                    ! now final position/velocity becomes next starting position/velocity
                    l_sub = l_f
                    v_sub = v_f
                end do
                
            end do loopParticles
            particleList(j)%N_p(iThread) = particleList(j)%N_p(iThread) - delIdx
            particleList(j)%delIdx(iThread) = delIdx
            !$OMP end parallel
            particleList(j)%accumEnergyLoss = particleList(j)%accumEnergyLoss + SUM(particleList(j)%energyLoss, DIM = 2)
            particleList(j)%accumWallLoss = particleList(j)%accumWallLoss + SUM(particleList(j)%wallLoss, DIM = 2)
        end do loopSpecies
    end subroutine moveParticles

end module mod_particleMover