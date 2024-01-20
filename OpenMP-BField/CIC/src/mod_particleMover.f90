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
    integer(int32), parameter :: maxPartIter = 50


contains

    subroutine getDelTauSubStepNoBField(l_sub, v_sub, del_tau, d_half, q_over_m, l_cell, E_left, E_right, E_x, dx_dl, AtBoundaryBool, FutureAtBoundaryBool, l_boundary)
        integer(int32), intent(in) :: l_cell
        real(real64), intent(in) :: q_over_m, dx_dl, E_left, E_right, l_sub, v_sub(3)
        real(real64), intent(in out) :: del_tau, d_half, E_x
        logical, intent(in) :: AtBoundaryBool
        logical, intent(in out) :: FutureAtBoundaryBool
        integer(int32), intent(in out) :: l_boundary
        integer(int32) :: l_alongV, l_awayV
        real(real64) :: del_tau_temp, a, b, c, l_f_alongV, l_f_awayV
        logical :: goingForwardBool
        ! get index cell where field and dx_dl is evaluated
        goingForwardBool = (v_sub(1) > 0)
        if (goingForwardBool) then
            l_alongV = l_cell + 1
            l_awayV = l_cell
        else
            l_alongV = l_cell
            l_awayV = l_cell + 1
        end if
       
        d_half = (l_sub + real(l_alongV)) * 0.5d0 - real(l_cell)
        E_x = (E_right * d_half + E_left * (1.0d0 - d_half))/dx_dl
        a = 0.5d0 * q_over_m * E_x
        ! Particle first between nodes, so solve quadratic for that particle depending on conditions
        c = (l_sub - real(l_alongV)) * dx_dl
        b = v_sub(1)**2 - 4.0d0*a*c
        FutureAtBoundaryBool = (b >= 0)
        if (FutureAtBoundaryBool) then
            del_tau_temp = 2.0d0 * ABS(c)/(ABS(v_sub(1)) + SQRT(b))
            FutureAtBoundaryBool = (del_tau_temp < del_tau)
            if (FutureAtBoundaryBool) then
                del_tau = del_tau_temp
                l_boundary = l_alongV
            end if
        end if
        if (.not. FutureAtBoundaryBool) then
            ! v and a opposite direction, boundary opposite direction of v
            d_half = (l_sub + real(l_awayV)) * 0.5d0 - real(l_cell)
            E_x = (E_right * d_half + E_left * (1.0d0 - d_half))/dx_dl
            a = 0.5d0 * q_over_m * E_x
            FutureAtBoundaryBool = XOR_op(a > 0, goingForwardBool)
            if (FutureAtBoundaryBool) then
                if (.not. AtBoundaryBool) then
                    c = (l_sub - real(l_awayV)) * dx_dl
                    del_tau_temp = 2.0d0 * ABS(c)/(SQRT(v_sub(1)**2 - 4.0d0*a*c) - ABS(v_sub(1)))
                else
                    del_tau_temp = ABS(v_sub(1))/ABS(a)
                end if
                FutureAtBoundaryBool = (del_tau_temp < del_tau)
                if (FutureAtBoundaryBool) then
                    del_tau = del_tau_temp
                    l_boundary = l_awayV
                end if
            end if
        end if
        

    end subroutine getDelTauSubStepNoBField

    subroutine analyticalParticleMoverNoBField(l_sub, v_sub, l_f, del_tau, d_half, q_over_m, l_cell, E_left, E_right, dx_dl)
        real(real64), intent(in) :: l_sub, v_sub(3), del_tau, q_over_m, E_left, E_right, dx_dl
        real(real64), intent(in out) :: l_f, d_half
        integer(int32), intent(in) :: l_cell
        real(real64) :: del_tau_sqr, real_l_cell, dx_sqr
        integer(int32) :: k
        del_tau_sqr = del_tau**2
        real_l_cell = real(l_cell, kind = real64)
        dx_sqr = dx_dl**2
        l_f = (2.0d0*E_left*del_tau_sqr*real_l_cell*q_over_m - E_left*del_tau_sqr*l_sub*q_over_m + 2.0d0*E_left*del_tau_sqr*q_over_m - &
            2.0d0*E_right*del_tau_sqr*real_l_cell*q_over_m + E_right*del_tau_sqr*l_sub*q_over_m + 4.0d0*del_tau*dx_dl*v_sub(1) + 4.0d0*dx_sqr*l_sub)&
            /(E_left*del_tau_sqr*q_over_m - E_right*del_tau_sqr*q_over_m + 4.0d0*dx_sqr)
        d_half = (l_f + l_sub)*0.5d0 - real(l_cell)

    
    end subroutine analyticalParticleMoverNoBField

    subroutine subStepSolverPI(l_sub, v_sub, l_f, v_f, v_half, del_tau, d_half, q_over_m, l_cell, E_left, E_right, BField, B_mag, dx_dl, f_tol, FutureAtBoundaryBool, l_boundary, numIter)
        ! picard iteration particle mover, like Chen in 2015 paper
        integer(int32), intent(in) :: l_cell
        real(real64), intent(in) :: q_over_m, BField(3), B_mag, dx_dl, f_tol, E_left, E_right
        real(real64), intent(in out) :: l_sub, v_sub(3), l_f, v_f(3), v_half(3), del_tau, d_half
        logical, intent(in out) :: FutureAtBoundaryBool
        integer(int32), intent(in out) :: l_boundary, numIter
        real(real64) :: v_prime(3), coeffAccel, curr_Res, l_f_prev, E_x, del_tau_prev
        integer(int32) :: k
        
        v_prime = v_sub
        del_tau_prev = del_tau
        l_f_prev = l_sub
        do k = 1, 100
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
            
            if (ABS(l_f_prev - l_f) < 1.d-10 .and. ABS(del_tau_prev - del_tau) < f_tol) then
                FutureAtBoundaryBool = (ABS(l_f - real(l_boundary)) < 1.d-10)
                if (FutureAtBoundaryBool) l_f = real(l_boundary)
                numIter = k
                if (.not. FutureAtBoundaryBool .and. INT(l_f) /= l_cell) then
                    print *, 'l_sub:', l_sub
                    print *, 'v_sub:', v_sub
                    print *, 'l_f:', l_f
                    print *, 'del_tau:', del_tau
                    print *, "not at boundary but l_f not in cell"
                    stop
                end if
                if (del_tau <= 0.0d0) then
                    print *, 'FutureAtBoundaryBool:', FutureAtBoundaryBool
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
            del_tau_prev = del_tau
            l_f_prev = l_f
        end do

        if (k > 30) then
            print *, 'more than 30 iterations for general PI solver'
            stop
        end if
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

    subroutine subStepSolverGetPosition(l_sub, v_sub, l_f, v_half, del_tau, d_half, q_over_m, l_cell, E_left, E_right, E_x, BField, B_mag, dx_dl, FutureAtBoundaryBool, l_boundary, numIter)
        ! Anderson Acceleration particle mover
        integer(int32), intent(in) :: l_cell
        real(real64), intent(in) :: q_over_m, BField(3), B_mag, dx_dl, E_left, E_right, l_sub, v_sub(3)
        real(real64), intent(in out) :: l_f, v_half(3), del_tau, d_half, E_x
        logical, intent(in out) :: FutureAtBoundaryBool
        integer(int32), intent(in out) :: l_boundary, numIter
        real(real64) :: v_prime(3), coeffAccel, Res_k(m_Anderson_Particle+1), param_k(m_Anderson_Particle+1), fitMat(m_Anderson_Particle), minCoeff(m_Anderson_Particle+1)
        integer(int32) :: k, index, m_k, u

        v_prime = v_sub
        v_half = v_sub
        ! print *, 'l_sub:', l_sub
        ! print *, 'v_sub:', v_sub
        ! print *, 'del_tau_max:', del_tau_max
        param_k(1) = l_sub
        coeffAccel = 0.5d0 * del_tau * q_over_m
        d_half = (l_sub + param_k(1)) * 0.5d0 - real(l_cell)
        E_x = (E_right * d_half + E_left * (1.0d0 - d_half))/dx_dl
        v_prime(1) = v_sub(1) + coeffAccel * E_x
        v_half = v_prime + coeffAccel * (crossProduct(v_prime, BField) + coeffAccel* SUM(v_prime * BField) * BField)
        v_half = v_half / (1.0d0 + (coeffAccel*B_mag)**2)
        l_f = l_sub + del_tau * v_half(1) / dx_dl
        if (v_half(1) > 0) then
            l_boundary = l_cell + 1
            FutureAtBoundaryBool = (l_f >= real(l_boundary))
            if (FutureAtBoundaryBool) l_f = real(l_boundary)
        else
            l_boundary = l_cell
            FutureAtBoundaryBool = (l_f <= real(l_boundary))
            if (FutureAtBoundaryBool) l_f = real(l_boundary)
        end if
        param_k(2) = l_f
        Res_k(1) = param_k(2) - param_k(1)

        do k = 1, maxPartIter
            index = MODULO(k, m_Anderson_Particle+1) + 1
            m_k = MIN(k, m_Anderson_Particle)
            coeffAccel = 0.5d0 * del_tau * q_over_m
            d_half = (l_sub + param_k(index)) * 0.5d0 - real(l_cell)
            E_x = (E_right * d_half + E_left * (1.0d0 - d_half))/dx_dl
            v_prime(1) = v_sub(1) + coeffAccel * E_x
            v_half = v_prime + coeffAccel * (crossProduct(v_prime, BField) + coeffAccel* SUM(v_prime * BField) * BField)
            v_half = v_half / (1.0d0 + (coeffAccel*B_mag)**2)
            l_f = l_sub + del_tau * v_half(1) / dx_dl
            if (v_half(1) > 0) then
                l_boundary = l_cell + 1
                FutureAtBoundaryBool = (l_f >= real(l_boundary))
                if (FutureAtBoundaryBool) l_f = real(l_boundary)
            else
                l_boundary = l_cell
                FutureAtBoundaryBool = (l_f <= real(l_boundary))
                if (FutureAtBoundaryBool) l_f = real(l_boundary)
            end if
            Res_k(index) = l_f - param_k(index)

            if (ABS(Res_k(index)) < 1.d-10) then
                numIter = k+1
                if (.not. FutureAtBoundaryBool .and. INT(l_f) /= l_cell) then
                    print *, 'l_f should be in same cell since not at boundary!'
                    stop
                end if
                exit      
            end if
            do u = 0, m_k-1
                fitMat(u+1) = Res_k(MODULO(k - m_k + u, m_Anderson_Particle+1) + 1) - Res_k(index)
            end do
            minCoeff(1:m_k) = (fitMat(1:m_k) * -Res_k(index))/SUM(fitMat(1:m_k)**2)
            minCoeff(m_k+1) = 1.0d0 - SUM(minCoeff(1:m_k))
            param_k(MODULO(k+1, m_Anderson_Particle+1) + 1) = minCoeff(1) * (Res_k(MODULO(k-m_k, m_Anderson_Particle+1) + 1) + param_k(MODULO(k-m_k,m_Anderson_Particle+1) + 1))
            do u=1, m_k
                param_k(MODULO(k+1, m_Anderson_Particle+1) + 1) = param_k(MODULO(k+1, m_Anderson_Particle+1) + 1) + minCoeff(u + 1) * (Res_k(MODULO(k-m_k + u, m_Anderson_Particle+1) + 1) + param_k(MODULO(k-m_k + u, m_Anderson_Particle+1) + 1))
            end do
        end do
        if (k > maxPartIter) numIter = maxPartIter


    end subroutine subStepSolverGetPosition

    subroutine subStepSolverGetPositionBrent(l_sub, v_sub, l_f, v_half, del_tau, d_half, q_over_m, l_cell, E_left, E_right, E_x, BField, B_mag, dx_dl, FutureAtBoundaryBool, l_boundary, numIter)
        ! Anderson Acceleration particle mover
        integer(int32), intent(in) :: l_cell
        real(real64), intent(in) :: q_over_m, BField(3), B_mag, dx_dl, E_left, E_right, l_sub, v_sub(3)
        real(real64), intent(in out) :: l_f, v_half(3), del_tau, d_half, E_x
        logical, intent(in out) :: FutureAtBoundaryBool
        integer(int32), intent(in out) :: l_boundary, numIter
        real(real64) :: v_prime(3), coeffAccel, Res_a, Res_b, Res_c, Res_s, l_f_a, l_f_b, l_f_c, l_f_d, l_f_s, tempVar
        integer(int32) :: k
        logical :: usedBiSection

        v_prime = v_sub
        coeffAccel = 0.5d0 * del_tau * q_over_m
        l_f_a = real(l_cell)
        l_f_b = real(l_cell + 1)

        
        d_half = (l_sub + l_f_a) * 0.5d0 - real(l_cell)
        E_x = (E_right * d_half + E_left * (1.0d0 - d_half))/dx_dl
        v_prime(1) = v_sub(1) + coeffAccel * E_x
        v_half = v_prime + coeffAccel * (crossProduct(v_prime, BField) + coeffAccel* SUM(v_prime * BField) * BField)
        v_half = v_half / (1.0d0 + (coeffAccel*B_mag)**2)
        l_f = l_sub + del_tau * v_half(1) / dx_dl
        if (l_f < l_f_a) then
            l_boundary = l_cell
            numIter = 1
            FutureAtBoundaryBool = .true.
            l_f = l_f_a
        else
            Res_a = l_f - l_f_a

            d_half = (l_sub + l_f_b) * 0.5d0 - real(l_cell)
            E_x = (E_right * d_half + E_left * (1.0d0 - d_half))/dx_dl
            v_prime(1) = v_sub(1) + coeffAccel * E_x
            v_half = v_prime + coeffAccel * (crossProduct(v_prime, BField) + coeffAccel* SUM(v_prime * BField) * BField)
            v_half = v_half / (1.0d0 + (coeffAccel*B_mag)**2)
            l_f = l_sub + del_tau * v_half(1) / dx_dl
            if (l_f > l_f_b) then
                l_boundary = l_cell+1
                numIter = 1
                FutureAtBoundaryBool = .true.
                l_f = l_f_b
            else
                FutureAtBoundaryBool = .false.
                Res_b = l_f - l_f_b
            end if
        end if

        if (.not. FutureAtBoundaryBool) then
            if (Res_a * Res_b >= 0) then
                print *, 'first res are not opposite!'
                stop
            end if

            if (abs(Res_a) < abs(Res_b) ) then
                ! if next Res higher, switch!
                tempVar = l_f_a
                l_f_a = l_f_b
                l_f_b = tempVar
                tempVar = Res_a
                Res_a = Res_b
                Res_b = tempVar
            end if

            l_f_c = l_f_a
            Res_c = Res_a
            usedBiSection = .true.

            do k = 1, maxPartIter
                if (Res_a /= Res_c .and. Res_b /= Res_c) then
                    l_f_s = l_f_a * Res_b * Res_c/((Res_a - Res_b) * (Res_a - Res_c)) + l_f_b * Res_a * Res_c/((Res_b - Res_a) * (Res_b - Res_c)) &
                        + l_f_c * Res_a * Res_b/((Res_c - Res_a) * (Res_c - Res_b))
                else
                    l_f_s = l_f_b - Res_b * (l_f_b - l_f_a)/(Res_b - Res_a)
                end if

                if ((l_f_s - 0.25d0 * (3.0d0 * l_f_a + l_f_b)) * (l_f_s - l_f_b) >= 0 &
                    .or.  (usedBiSection .and. ABS(l_f_s - l_f_b) >= 0.5d0 * ABS(l_f_b - l_f_c)) &
                    .or.  (.not. usedBiSection .and. ABS(l_f_s - l_f_b) >= 0.5d0 * ABS(l_f_c - l_f_d)) &
                    .or.  (usedBiSection .and. ABS(l_f_b - l_f_c) < 1.d-10) &
                    .or. (.not. usedBiSection .and. ABS(l_f_c - l_f_d) < 1.d-10)) then

                    l_f_s = 0.5d0 * (l_f_a + l_f_b)
                    usedBiSection = .true.
                else
                    usedBiSection = .false.
                end if
                if (ABS(l_f_s - l_f_b) < 1.d-10) then
                    l_f = l_f_s
                    numIter = k+1
                    exit
                end if
                d_half = (l_sub + l_f_s) * 0.5d0 - real(l_cell)
                E_x = (E_right * d_half + E_left * (1.0d0 - d_half))/dx_dl
                v_prime(1) = v_sub(1) + coeffAccel * E_x
                v_half = v_prime + coeffAccel * (crossProduct(v_prime, BField) + coeffAccel* SUM(v_prime * BField) * BField)
                v_half = v_half / (1.0d0 + (coeffAccel*B_mag)**2)
                l_f = l_sub + del_tau * v_half(1) / dx_dl
                Res_s = l_f - l_f_s

                l_f_d = l_f_c
                l_f_c = l_f_b
                Res_c = Res_b

                if (Res_a * Res_s < 0) then
                    l_f_b = l_f_s
                    Res_b = Res_s
                else
                    l_f_a = l_f_s
                    Res_a = Res_s
                end if

                if (abs(Res_a) < abs(Res_b) ) then
                    ! if next Res higher, switch!
                    tempVar = l_f_a
                    l_f_a = l_f_b
                    l_f_b = tempVar
                    tempVar = Res_a
                    Res_a = Res_b
                    Res_b = tempVar
                end if       
                
            end do

            if (k > maxPartIter) then
                print *, 'issue solving time with Brent solver in get position'
                print *, 'l_sub:', l_sub
                print *, 'del_tau_max:', del_tau
                print *, 'l_boundary:', l_boundary
                print *, 'v_sub:', v_sub
                stop
            end if
        end if

    end subroutine subStepSolverGetPositionBrent

    subroutine subStepSolverGetTimeBoundaryBrent(l_sub, v_sub, l_f, v_half, del_tau, q_over_m, E_x, BField, B_mag, dx_dl, f_tol, l_boundary, numIter)
        ! Anderson Acceleration particle mover
        real(real64), intent(in) :: q_over_m, l_sub, v_sub(3), BField(3), B_mag, dx_dl, E_x, f_tol
        real(real64), intent(in out) :: l_f, v_half(3), del_tau
        integer(int32), intent(in) :: l_boundary
        integer(int32), intent(in out) :: numIter
        real(real64) :: v_prime(3), coeffAccel, Res_a, Res_b, Res_c, Res_s, del_tau_a, del_tau_b, del_tau_c, del_tau_d, del_tau_s, tempVar
        integer(int32) :: k
        logical :: usedBiSection

        v_prime = v_sub

        del_tau_a = 0.0d0
        del_tau_b = del_tau

        ! coeffAccel = 0.5d0 * del_tau_a * q_over_m
        ! v_prime(1) = v_sub(1) + coeffAccel * E_x
        ! v_half = v_prime + coeffAccel * (crossProduct(v_prime, BField) + coeffAccel* SUM(v_prime * BField) * BField)
        ! v_half = v_half / (1.0d0 + (coeffAccel*B_mag)**2)
        ! l_f = l_sub + del_tau_a * v_half(1) / dx_dl
        Res_a = l_sub - real(l_boundary)
            
        coeffAccel = 0.5d0 * del_tau_b * q_over_m
        v_prime(1) = v_sub(1) + coeffAccel * E_x
        v_half = v_prime + coeffAccel * (crossProduct(v_prime, BField) + coeffAccel* SUM(v_prime * BField) * BField)
        v_half = v_half / (1.0d0 + (coeffAccel*B_mag)**2)
        l_f = l_sub + del_tau_b * v_half(1) / dx_dl
        Res_b = l_f - real(l_boundary)

        if (Res_a * Res_b >= 0) then
            print *, 'first res are not opposite!'
            stop
        end if

        if (abs(Res_a) < abs(Res_b) ) then
            ! if next Res higher, switch!
            tempVar = del_tau_a
            del_tau_a = del_tau_b
            del_tau_b = tempVar
            tempVar = Res_a
            Res_a = Res_b
            Res_b = tempVar
        end if

        del_tau_c = del_tau_a
        Res_c = Res_a
        usedBiSection = .true.

        do k = 1, maxPartIter
            if (Res_a /= Res_c .and. Res_b /= Res_c) then
                del_tau_s = del_tau_a * Res_b * Res_c/((Res_a - Res_b) * (Res_a - Res_c)) + del_tau_b * Res_a * Res_c/((Res_b - Res_a) * (Res_b - Res_c)) &
                    + del_tau_c * Res_a * Res_b/((Res_c - Res_a) * (Res_c - Res_b))
            else
                del_tau_s = del_tau_b - Res_b * (del_tau_b - del_tau_a)/(Res_b - Res_a)
            end if

            if ((del_tau_s - 0.25d0 * (3.0d0 * del_tau_a + del_tau_b)) * (del_tau_s - del_tau_b) >= 0 &
                .or.  (usedBiSection .and. ABS(del_tau_s - del_tau_b) >= 0.5d0 * ABS(del_tau_b - del_tau_c)) &
                .or.  (.not. usedBiSection .and. ABS(del_tau_s - del_tau_b) >= 0.5d0 * ABS(del_tau_c - del_tau_d)) &
                .or.  (usedBiSection .and. ABS(del_tau_b - del_tau_c) < f_tol) &
                .or. (.not. usedBiSection .and. ABS(del_tau_c - del_tau_d) < f_tol)) then

                del_tau_s = 0.5d0 * (del_tau_a + del_tau_b)
                usedBiSection = .true.
            else
                usedBiSection = .false.
            end if
            if (ABS(del_tau_s - del_tau_b) < f_tol) then
                del_tau = del_tau_s
                l_f = real(l_boundary)
                numIter = k
                exit
            end if
            coeffAccel = 0.5d0 * del_tau_s * q_over_m
            v_prime(1) = v_sub(1) + coeffAccel * E_x
            v_half = v_prime + coeffAccel * (crossProduct(v_prime, BField) + coeffAccel* SUM(v_prime * BField) * BField)
            v_half = v_half / (1.0d0 + (coeffAccel*B_mag)**2)
            l_f = l_sub + del_tau_s * v_half(1) / dx_dl
            Res_s = l_f - real(l_boundary)

            del_tau_d = del_tau_c
            del_tau_c = del_tau_b
            Res_c = Res_b

            if (Res_a * Res_s < 0) then
                del_tau_b = del_tau_s
                Res_b = Res_s
            else
                del_tau_a = del_tau_s
                Res_a = Res_s
            end if

            if (abs(Res_a) < abs(Res_b) ) then
                ! if next Res higher, switch!
                tempVar = del_tau_a
                del_tau_a = del_tau_b
                del_tau_b = tempVar
                tempVar = Res_a
                Res_a = Res_b
                Res_b = tempVar
            end if       
            
        end do

        if (k > maxPartIter) then
            print *, 'issue solving time with Brent solver'
            print *, 'l_sub:', l_sub
            print *, 'del_tau_max:', del_tau
            print *, 'l_boundary:', l_boundary
            print *, 'v_sub:', v_sub
            stop
        end if
        

    end subroutine subStepSolverGetTimeBoundaryBrent

    subroutine subStepSolverGetTimeSelfBoundaryBrent(l_sub, v_sub, v_half, del_tau, q_over_m, E_x, BField, B_mag, dx_dl, f_tol, l_boundary, numIter)
        ! When particle loops back on itself on boundary
        real(real64), intent(in) :: q_over_m, l_sub, v_sub(3), BField(3), B_mag, dx_dl, E_x, f_tol
        real(real64), intent(in out) :: v_half(3), del_tau
        integer(int32), intent(in) :: l_boundary
        integer(int32), intent(in out) :: numIter
        real(real64) :: v_prime(3), coeffAccel, Res_a, Res_b, Res_c, Res_s, del_tau_a, del_tau_b, del_tau_c, del_tau_d, del_tau_s, tempVar
        integer(int32) :: k
        logical :: usedBiSection

        v_prime = v_sub

        del_tau_a = 0.0d0
        del_tau_b = del_tau
       
        ! coeffAccel = 0.5d0 * del_tau_a * q_over_m
        ! v_prime(1) = v_sub(1) + coeffAccel * E_x
        ! v_half = v_prime + coeffAccel * (crossProduct(v_prime, BField) + coeffAccel* SUM(v_prime * BField) * BField)
        ! v_half = v_half / (1.0d0 + (coeffAccel*B_mag)**2)
        Res_a = v_sub(1)
          
        coeffAccel = 0.5d0 * del_tau_b * q_over_m
        v_prime(1) = v_sub(1) + coeffAccel * E_x
        v_half = v_prime + coeffAccel * (crossProduct(v_prime, BField) + coeffAccel* SUM(v_prime * BField) * BField)
        v_half = v_half / (1.0d0 + (coeffAccel*B_mag)**2)
        Res_b = v_half(1)

        if (Res_a * Res_b >= 0) then
            print *, 'first res are not opposite!'
            print *, 'l_sub:', l_sub
            print *, 'v_sub:', v_sub
            print *, 'Res_a:', Res_a
            print *, 'Res_b:', Res_b
            stop
        end if

        if (abs(Res_a) < abs(Res_b) ) then
            ! if next Res higher, switch!
            tempVar = del_tau_a
            del_tau_a = del_tau_b
            del_tau_b = tempVar
            tempVar = Res_a
            Res_a = Res_b
            Res_b = tempVar
        end if

        del_tau_c = del_tau_a
        Res_c = Res_a
        usedBiSection = .true.

        do k = 1, maxPartIter
            if (Res_a /= Res_c .and. Res_b /= Res_c) then
                del_tau_s = del_tau_a * Res_b * Res_c/((Res_a - Res_b) * (Res_a - Res_c)) + del_tau_b * Res_a * Res_c/((Res_b - Res_a) * (Res_b - Res_c)) &
                    + del_tau_c * Res_a * Res_b/((Res_c - Res_a) * (Res_c - Res_b))
            else
                del_tau_s = del_tau_b - Res_b * (del_tau_b - del_tau_a)/(Res_b - Res_a)
            end if

            if ((del_tau_s - 0.25d0 * (3.0d0 * del_tau_a + del_tau_b)) * (del_tau_s - del_tau_b) >= 0 &
                .or.  (usedBiSection .and. ABS(del_tau_s - del_tau_b) >= 0.5d0 * ABS(del_tau_b - del_tau_c)) &
                .or.  (.not. usedBiSection .and. ABS(del_tau_s - del_tau_b) >= 0.5d0 * ABS(del_tau_c - del_tau_d)) &
                .or.  (usedBiSection .and. ABS(del_tau_b - del_tau_c) < f_tol) &
                .or. (.not. usedBiSection .and. ABS(del_tau_c - del_tau_d) < f_tol)) then

                del_tau_s = 0.5d0 * (del_tau_a + del_tau_b)
                usedBiSection = .true.
            else
                usedBiSection = .false.
            end if

            if (ABS(del_tau_s - del_tau_b) < f_tol) then
                del_tau = del_tau_s
                numIter = k
                exit
            end if

            coeffAccel = 0.5d0 * del_tau_s * q_over_m
            v_prime(1) = v_sub(1) + coeffAccel * E_x
            v_half = v_prime + coeffAccel * (crossProduct(v_prime, BField) + coeffAccel* SUM(v_prime * BField) * BField)
            v_half = v_half / (1.0d0 + (coeffAccel*B_mag)**2)
            Res_s = v_half(1)
            ! if (ABS(del_tau_s-del_tau_b) < f_tol) then
            !     del_tau = del_tau_s
            !     print *, 'del_tau:', del_tau
            !     print *, 'iteration:', k
            !     nonLinCounter = nonLinCounter + k
            !     exit
            ! end if

            del_tau_d = del_tau_c
            del_tau_c = del_tau_b
            Res_c = Res_b

            if (Res_a * Res_s < 0) then
                del_tau_b = del_tau_s
                Res_b = Res_s
            else
                del_tau_a = del_tau_s
                Res_a = Res_s
            end if

            if (abs(Res_a) < abs(Res_b) ) then
                ! if next Res higher, switch!
                tempVar = del_tau_a
                del_tau_a = del_tau_b
                del_tau_b = tempVar
                tempVar = Res_a
                Res_a = Res_b
                Res_b = tempVar
            end if       
            
        end do

        if (k > maxPartIter) then
            print *, 'issue solving time with Brent solver'
            print *, 'l_sub:', l_sub
            print *, 'del_tau_max:', del_tau
            print *, 'l_boundary:', l_boundary
            print *, 'v_sub:', v_sub
            stop
        end if
        

    end subroutine subStepSolverGetTimeSelfBoundaryBrent
    
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
        else if (world%boundaryConditions(l_right) == 4) then
            !Neumann absorbing to right
            rho(l_left) = rho(l_left) + 0.5d0 * (1.0d0 - d)**2
        else if (world%boundaryConditions(l_center) == 4) then
            !Neumann absorbing to left
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
        real(real64) :: l_f, l_sub, v_sub(3), v_f(3), timePassed, del_tau, q_over_m, f_tol, v_half(3), dx_dl, q_times_wp, d_half, E_x, J_part
        integer(int32) :: j, i, l_cell, iThread, l_boundary, numIter
        logical :: FutureAtBoundaryBool, AtBoundaryBool, timeNotConvergedBool
        !$OMP parallel private(iThread)
        iThread = omp_get_thread_num() + 1 
        solver%J(:, iThread) = 0.0d0
        !$OMP end parallel
        call solver%evaluateEFieldHalfTime(world)
        f_tol = del_t * 1.d-10
        loopSpecies: do j = 1, numberChargedParticles
            q_over_m = particleList(j)%q/particleList(j)%mass
            q_times_wp = particleList(j)%q * particleList(j)%w_p
            !$OMP parallel private(iThread, i, l_f, l_sub, v_sub, v_f, v_half, timePassed, del_tau, E_x, l_cell, FutureAtBoundaryBool, AtBoundaryBool, &
                dx_dl, l_boundary, d_half, numIter, J_part, timeNotConvergedBool)
            iThread = omp_get_thread_num() + 1 
            loopParticles: do i = 1, particleList(j)%N_p(iThread)
                v_sub = particleList(j)%phaseSpace(2:4,i,iThread)
                l_sub = particleList(j)%phaseSpace(1,i,iThread)
                timePassed = 0.0d0
                AtBoundaryBool = MOD(l_sub, 1.0d0) == 0.0d0
                timeNotConvergedBool = .true.
                del_tau = del_t
                do while(timeNotConvergedBool)
                    if (.not. AtBoundaryBool) then
                        l_cell = INT(l_sub)
                    else
                        l_cell = INT(l_sub) + (INT(SIGN(1.0d0, v_sub(1))) - 1)/2
                    end if
                    dx_dl = world%dx_dl(l_cell)
                    ! Start AA

                    if (.not. solver%BFieldBool) then
                        call getDelTauSubStepNoBField(l_sub, v_sub, del_tau, d_half, q_over_m, l_cell, solver%EField(l_cell), solver%EField(l_cell+1), E_x, dx_dl, AtBoundaryBool, FutureAtBoundaryBool, l_boundary)
                        if (FutureAtBoundaryBool) then
                            l_f = real(l_boundary)
                        else
                            call analyticalParticleMoverNoBField(l_sub, v_sub, l_f, del_tau, d_half, q_over_m, l_cell, solver%EField(l_cell), solver%EField(l_cell+1), dx_dl)
                        end if
                        v_half(1) = (l_f - l_sub) * dx_dl/del_tau
                        v_f(1) = 2.0d0 * v_half(1) - v_sub(1)
                        v_f(2:3) = v_sub(2:3)
                    else
                        call subStepSolverGetPosition(l_sub, v_sub, l_f, v_half, del_tau, d_half, q_over_m, l_cell, solver%EField(l_cell), solver%EField(l_cell+1), E_x, solver%BField, &
                            solver%BFieldMag, dx_dl, FutureAtBoundaryBool, l_boundary, numIter)

                        if (numIter == maxPartIter) call subStepSolverGetPositionBrent(l_sub, v_sub, l_f, v_half, del_tau, d_half, q_over_m, l_cell, solver%EField(l_cell), solver%EField(l_cell+1), E_x, solver%BField, &
                            solver%BFieldMag, dx_dl, FutureAtBoundaryBool, l_boundary, numIter)

                        if (FutureAtBoundaryBool) then
                            if (.not. AtBoundaryBool .or. l_boundary /= INT(l_sub)) then
                                call subStepSolverGetTimeBoundaryBrent(l_sub, v_sub, l_f, v_half, del_tau, q_over_m, E_x, solver%BField, solver%BFieldMag, dx_dl, f_tol, l_boundary, numIter)
                            else
                                call subStepSolverGetTimeSelfBoundaryBrent(l_sub, v_sub, v_half, del_tau, q_over_m, E_x, solver%BField, solver%BFieldMag, dx_dl, f_tol, l_boundary, numIter)
                            end if
                        end if
                        v_f = 2.0d0 * v_half - v_sub
                    end if
        
                    J_part = q_times_wp * (v_half(1))*del_tau/dx_dl/del_t
                    solver%J(l_cell, iThread) = solver%J(l_cell, iThread) + J_part * (1.0d0 - d_half)
                    solver%J(l_cell+1, iThread) = solver%J(l_cell+1, iThread) + J_part * (d_half)

                    if (FutureAtBoundaryBool) then
                        SELECT CASE (world%boundaryConditions(l_boundary))
                        CASE(0)
                            l_f = real(l_boundary)
                        CASE(1,4)
                            exit
                        CASE(2)
                            if (l_boundary == NumberXNodes+1) then
                                if (v_f(1) > 0) then
                                    v_f = -v_f
                                end if
                            else
                                if (v_f(1) < 0) then
                                    v_f = -v_f
                                end if
                            end if
                            l_f = real(l_boundary) 
                        CASE(3)
                            l_f = REAL(ABS(l_boundary - real(NumberXNodes+1, kind = real64) - 1))
                        CASE default
                            print *, "l_sub is:", l_sub
                            print *, 'l_f is:', l_f
                            print *, world%boundaryConditions(INT(l_f))
                            print *, "Case does not exist in ongoing substep, depositJ"
                            stop
                        END SELECT
                    end if
                    if ((l_f < l_cell) .or. (l_f > l_cell+1)) then
                        print *, 'l_f is:', l_f
                        print *, 'l_sub:', l_sub
                        print *, 'del_tau:', del_tau
                        print *, 'v_sub:', v_sub
                        print *, 'l_cell:', l_cell
                        print *, 'AtBoundaryBool:', AtBoundaryBool
                        print *, 'FutureAtBoundaryBool:', FutureAtBoundaryBool
                        stop "Have particles travelling outside local cell!"
                    end if
                    timePassed = timePassed + del_tau
                    del_tau = del_t - timePassed
                    timeNotConvergedBool = (del_tau > f_tol)
                    ! now final position/velocity becomes next starting position/velocity
                    l_sub = l_f
                    v_sub = v_f
                    AtBoundaryBool = FutureAtBoundaryBool
                    
                end do
            end do loopParticles
            !$OMP end parallel
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
        real(real64) :: l_f, l_sub, v_sub(3), v_f(3), timePassed, del_tau, q_over_m, f_tol, d_half, v_half(3), dx_dl, E_x
        integer(int32) :: j, i, l_cell, iThread, l_boundary, numSubStepAve, numIter, funcEvalCounter, delIdx, refIdx
        logical :: FutureAtBoundaryBool, AtBoundaryBool, refluxedBool, timeNotConvergedBool
        call solver%evaluateEFieldHalfTime(world)
        f_tol = del_t * 1.d-10
        loopSpecies: do j = 1, numberChargedParticles
            q_over_m = particleList(j)%q/particleList(j)%mass
            numSubStepAve = 0
            funcEvalCounter = 0
            !$OMP parallel private(iThread, i, l_f, l_sub, v_sub, v_f, v_half, timePassed, del_tau, E_x, l_cell, FutureAtBoundaryBool, &
                AtBoundaryBool, dx_dl, l_boundary, d_half, numIter, delIdx, refIdx, refluxedBool, timeNotConvergedBool) reduction(+:numSubStepAve, funcEvalCounter)
            iThread = omp_get_thread_num() + 1 
            delIdx = 0
            refIdx = 0
            particleList(j)%energyLoss(:, iThread) = 0.0d0
            particleList(j)%wallLoss(:, iThread) = 0.0d0
            loopParticles: do i = 1, particleList(j)%N_p(iThread)
                v_sub = particleList(j)%phaseSpace(2:4,i,iThread)
                l_sub = particleList(j)%phaseSpace(1,i,iThread)
                timePassed = 0.0d0
                AtBoundaryBool = MOD(l_sub, 1.0d0) == 0.0d0
                refluxedBool = .false.
                timeNotConvergedBool = .true.
                del_tau = del_t
                do while(timeNotConvergedBool)
                    numSubStepAve = numSubStepAve + 1
                    if (.not. AtBoundaryBool) then
                        l_cell = INT(l_sub)
                    else
                        l_cell = INT(l_sub) + (INT(SIGN(1.0d0, v_sub(1))) - 1)/2
                    end if
                    dx_dl = world%dx_dl(l_cell)
            
                    ! AA particle mover
                    if (.not. solver%BFieldBool) then
                        call getDelTauSubStepNoBField(l_sub, v_sub, del_tau, d_half, q_over_m, l_cell, solver%EField(l_cell), solver%EField(l_cell+1), E_x, dx_dl, AtBoundaryBool, FutureAtBoundaryBool, l_boundary)
                        if (FutureAtBoundaryBool) then
                            l_f = real(l_boundary)
                        else
                            call analyticalParticleMoverNoBField(l_sub, v_sub, l_f, del_tau, d_half, q_over_m, l_cell, solver%EField(l_cell), solver%EField(l_cell+1), dx_dl)
                        end if
                        v_half(1) = (l_f - l_sub) * dx_dl/del_tau
                        v_f(1) = 2.0d0 * v_half(1) - v_sub(1)
                        v_f(2:3) = v_sub(2:3)
                        funcEvalCounter = funcEvalCounter + 1
                    else
                        call subStepSolverGetPosition(l_sub, v_sub, l_f, v_half, del_tau, d_half, q_over_m, l_cell, solver%EField(l_cell), solver%EField(l_cell+1), E_x, solver%BField, &
                            solver%BFieldMag, dx_dl, FutureAtBoundaryBool, l_boundary, numIter)

                        if (numIter /= maxPartIter) then
                            funcEvalCounter = funcEvalCounter + numIter
                        else
                            call subStepSolverGetPositionBrent(l_sub, v_sub, l_f, v_half, del_tau, d_half, q_over_m, l_cell, solver%EField(l_cell), solver%EField(l_cell+1), E_x, solver%BField, &
                            solver%BFieldMag, dx_dl, FutureAtBoundaryBool, l_boundary, numIter)
                            funcEvalCounter = funcEvalCounter + numIter + maxPartIter
                        end if
                        if (FutureAtBoundaryBool) then
                            if (.not. AtBoundaryBool .or. l_boundary /= INT(l_sub)) then
                                call subStepSolverGetTimeBoundaryBrent(l_sub, v_sub, l_f, v_half, del_tau, q_over_m, E_x, solver%BField, solver%BFieldMag, dx_dl, f_tol, l_boundary, numIter)
                            else
                                call subStepSolverGetTimeSelfBoundaryBrent(l_sub, v_sub, v_half, del_tau, q_over_m, E_x, solver%BField, solver%BFieldMag, dx_dl, f_tol, l_boundary, numIter)
                            end if
                            funcEvalCounter = funcEvalCounter + numIter
                        end if
                        v_f = 2.0d0 * v_half - v_sub
                    end if
        
                    if (FutureAtBoundaryBool) then
                        SELECT CASE (world%boundaryConditions(l_boundary))
                        CASE(0)
                            l_f = real(l_boundary)
                        CASE(1,4)
                            delIdx = delIdx + 1
                            if (l_boundary == 1) then
                                particleList(j)%energyLoss(1, iThread) = particleList(j)%energyLoss(1, iThread) + (SUM(v_f**2))!J/m^2 in 1D
                                particleList(j)%wallLoss(1, iThread) = particleList(j)%wallLoss(1, iThread) + 1 !C/m^2 in 1D
                            else if (l_boundary == NumberXNodes+1) then
                                particleList(j)%energyLoss(2, iThread) = particleList(j)%energyLoss(2, iThread) + (SUM(v_f**2)) !J/m^2 in 1D
                                particleList(j)%wallLoss(2, iThread) = particleList(j)%wallLoss(2, iThread) + 1 !C/m^2 in 1D
                            end if
                            exit
                        CASE(2)
                            if (l_boundary == NumberXNodes+1) then
                                if (v_f(1) > 0) then
                                    v_f = -v_f
                                end if
                            else
                                if (v_f(1) < 0) then
                                    v_f = -v_f
                                end if
                            end if
                            l_f = real(l_boundary) 
                            if (.not. refluxedBool) then
                                refIdx = refIdx + 1
                                particleList(j)%refRecordIdx(refIdx, iThread) = i - delIdx
                                refluxedBool = .true.
                            end if
                        CASE(3)
                            l_f = REAL(ABS(l_boundary - real(NumberXNodes+1, kind = real64) - 1))
                        CASE default
                            print *, "l_sub is:", l_sub
                            print *, 'l_f is:', l_f
                            print *, world%boundaryConditions(INT(l_f))
                            print *, "Case does not exist in ongoing substep, depositJ"
                            stop
                        END SELECT
                    end if
                    if ((l_f < l_cell) .or. (l_f > l_cell+1)) then
                        print *, 'l_f is:', l_f
                        print *, 'l_sub:', l_sub
                        print *, 'del_tau:', del_tau
                        print *, 'v_sub:', v_sub
                        print *, 'l_cell:', l_cell
                        print *, 'AtBoundaryBool:', AtBoundaryBool
                        print *, 'FutureAtBoundaryBool:', FutureAtBoundaryBool
                        stop "Have particles travelling outside local cell!"
                    end if
                    timePassed = timePassed + del_tau
                    del_tau = del_t - timePassed
                    timeNotConvergedBool = (del_tau > f_tol) ! Substraction of bool may not have del_t == timePassed, so just use f_tol
                    ! now final position/velocity becomes next starting position/velocity
                    l_sub = l_f
                    v_sub = v_f
                    AtBoundaryBool = FutureAtBoundaryBool
                    
                end do
            
                if (.not. timeNotConvergedBool) then
                    particleList(j)%phaseSpace(1, i-delIdx, iThread) = l_f
                    particleList(j)%phaseSpace(2:4,i-delIdx, iThread) = v_f
                end if
            end do loopParticles
            particleList(j)%N_p(iThread) = particleList(j)%N_p(iThread) - delIdx
            particleList(j)%delIdx(iThread) = delIdx
            particleList(j)%refIdx(iThread) = refIdx
            !$OMP end parallel
            particleList(j)%numSubStepsAve = real(numSubStepAve) / real(SUM(particleList(j)%N_p) + SUM(particleList(j)%delIdx))
            particleList(j)%numFuncEvalAve = real(funcEvalCounter) / real(SUM(particleList(j)%N_p) + SUM(particleList(j)%delIdx))
            particleList(j)%accumEnergyLoss = particleList(j)%accumEnergyLoss + SUM(particleList(j)%energyLoss, DIM = 2)
            particleList(j)%accumWallLoss = particleList(j)%accumWallLoss + SUM(particleList(j)%wallLoss, DIM = 2)
        end do loopSpecies
    end subroutine moveParticles

end module mod_particleMover