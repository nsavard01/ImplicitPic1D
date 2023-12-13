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
        real(real64) :: l_f, l_sub, v_sub(3), v_f(3), timePassed, del_tau, q_over_m, f_tol, v_half(3), dx_dl, q_times_wp, J_part, d_half, E_x, J_temp(NumberXNodes+1)
        integer(int32) :: j, i, l_cell, iThread, l_boundary, numIter
        logical :: FutureAtBoundaryBool, AtBoundaryBool
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
                dx_dl, l_boundary, J_part, d_half, numIter, J_temp)
            iThread = omp_get_thread_num() + 1 
            J_temp = 0.0d0
            loopParticles: do i = 1, particleList(j)%N_p(iThread)
                v_sub = particleList(j)%phaseSpace(2:4,i,iThread)
                l_sub = particleList(j)%phaseSpace(1,i,iThread)
                timePassed = 0.0d0
                AtBoundaryBool = .false.
                do while((timePassed < del_t))
                    if (.not. AtBoundaryBool) then
                        l_cell = INT(l_sub)
                    else
                        l_cell = l_boundary + (INT(SIGN(1.0d0, v_sub(1))) - 1)/2
                    end if
                    dx_dl = world%dx_dl(l_cell)
                    del_tau = del_t - timePassed
                    ! Start AA
                    ! call subStepSolverPI(l_sub, v_sub, l_f, v_f, v_half, del_tau, d_half, q_over_m, l_cell, solver%EField(l_cell), solver%EField(l_cell+1), solver%BField, &
                    !     solver%BFieldMag, dx_dl, f_tol, FutureAtBoundaryBool, l_boundary, numIter)
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
                    J_part = q_times_wp * (v_half(1))*del_tau/dx_dl/del_t
                    J_temp(l_cell) = J_temp(l_cell) + J_part * (1.0d0 - d_half)
                    J_temp(l_cell+1) = J_temp(l_cell+1) + J_part * (d_half)
                    if (FutureAtBoundaryBool) then
                        SELECT CASE (world%boundaryConditions(l_boundary))
                        CASE(0)
                            l_f = real(l_boundary)
                        CASE(1)
                            exit
                        CASE(2)
                            if (l_boundary == NumberXNodes+1) then
                                v_f(1) = -ABS(v_f(1))  
                            else
                                v_f(1) = ABS(v_f(1))
                            end if
                            l_f = real(l_boundary)
                            v_f(2:3) = -v_f(2:3)  
                        CASE(3)
                            l_f = REAL(ABS(l_boundary - real(NumberXNodes, kind = real64) - 1))
                        CASE(4)
                            J_temp(l_boundary) = J_temp(l_boundary) + q_times_wp * v_f(1) * (del_t - timePassed-del_tau) * dx_dl / del_t
                            exit
                        CASE default
                            print *, "l_sub is:", l_sub
                            print *, 'l_f is:', l_f
                            print *, world%boundaryConditions(INT(l_f))
                            print *, "Case does not exist in ongoing substep, depositJ"
                            stop
                        END SELECT
                    end if
                    timePassed = timePassed + del_tau
                    ! now final position/velocity becomes next starting position/velocity
                    l_sub = l_f
                    v_sub = v_f
                    AtBoundaryBool = FutureAtBoundaryBool
                    if ((l_f < 1) .or. (l_f > NumberXNodes+1)) then
                        print *, 'l_f is:', l_f
                        stop "Have particles travelling outside the domain in depositJ!"
                    end if
                end do
            end do loopParticles
            solver%J(:, iThread) = solver%J(:, iThread) + J_temp
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
        logical :: FutureAtBoundaryBool, AtBoundaryBool, refluxedBool
        call solver%evaluateEFieldHalfTime(world)
        f_tol = del_t * 1.d-10
        loopSpecies: do j = 1, numberChargedParticles
            q_over_m = particleList(j)%q/particleList(j)%mass
            numSubStepAve = 0
            funcEvalCounter = 0
            !$OMP parallel private(iThread, i, l_f, l_sub, v_sub, v_f, v_half, timePassed, del_tau, E_x, l_cell, FutureAtBoundaryBool, &
                AtBoundaryBool, dx_dl, l_boundary, d_half, numIter, delIdx, refIdx, refluxedBool) reduction(+:numSubStepAve, funcEvalCounter)
            iThread = omp_get_thread_num() + 1 
            delIdx = 0
            refIdx = 0
            particleList(j)%energyLoss(:, iThread) = 0.0d0
            particleList(j)%wallLoss(:, iThread) = 0.0d0
            loopParticles: do i = 1, particleList(j)%N_p(iThread)
                v_sub = particleList(j)%phaseSpace(2:4,i,iThread)
                l_sub = particleList(j)%phaseSpace(1,i,iThread)
                timePassed = 0.0d0
                AtBoundaryBool = .false.
                refluxedBool = .false.
                do while((timePassed < del_t))
                    numSubStepAve = numSubStepAve + 1
                    if (.not. AtBoundaryBool) then
                        l_cell = INT(l_sub)
                    else
                        l_cell = l_boundary + (INT(SIGN(1.0d0, v_sub(1))) - 1)/2
                    end if
                    dx_dl = world%dx_dl(l_cell)
                    del_tau = del_t - timePassed
                    
                    ! AA particle mover

                    ! call subStepSolverPI(l_sub, v_sub, l_f, v_f, v_half, del_tau, d_half, q_over_m, l_cell, solver%EField(l_cell), solver%EField(l_cell+1), solver%BField, &
                    !     solver%BFieldMag, dx_dl, f_tol, FutureAtBoundaryBool, l_boundary, numIter)
                    call subStepSolverGetPosition(l_sub, v_sub, l_f, v_half, del_tau, d_half, q_over_m, l_cell, solver%EField(l_cell), solver%EField(l_cell+1), E_x, solver%BField, &
                        solver%BFieldMag, dx_dl, FutureAtBoundaryBool, l_boundary, numIter)

                    if (numIter /= maxPartIter) then
                        funcEvalCounter = funcEvalCounter + numIter
                    else
                        call subStepSolverGetPosition(l_sub, v_sub, l_f, v_half, del_tau, d_half, q_over_m, l_cell, solver%EField(l_cell), solver%EField(l_cell+1), E_x, solver%BField, &
                        solver%BFieldMag, dx_dl, FutureAtBoundaryBool, l_boundary, numIter)
                        funcEvalCounter = funcEvalCounter + numIter + maxPartIter
                    end if
                    if (FutureAtBoundaryBool) then
                        if (.not. AtBoundaryBool .or. l_boundary /= INT(l_sub)) then
                            call subStepSolverGetTimeBoundaryBrent(l_sub, v_sub, l_f, v_half, del_tau, q_over_m, E_x, solver%BField, solver%BFieldMag, dx_dl, f_tol, l_boundary, numIter)
                        else
                            call subStepSolverGetTimeSelfBoundaryBrent(l_sub, v_sub, v_half, del_tau, q_over_m, E_x, solver%BField, solver%BFieldMag, dx_dl, f_tol, l_boundary, numIter)
                        end if
                    end if
                    funcEvalCounter = funcEvalCounter + numIter
                    v_f = 2.0d0 * v_half - v_sub
                    if (FutureAtBoundaryBool) then
                        SELECT CASE (world%boundaryConditions(l_boundary))
                        CASE(0)
                            l_f = real(l_boundary)
                        CASE(1)
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
                                v_f(1) = -ABS(v_f(1))  
                            else
                                v_f(1) = ABS(v_f(1))
                            end if
                            l_f = real(l_boundary)
                            v_f(2:3) = -v_f(2:3)  
                            if (.not. refluxedBool) then
                                refIdx = refIdx + 1
                                particleList(j)%refRecordIdx(refIdx, iThread) = i - delIdx
                                refluxedBool = .true.
                            end if
                        CASE(3)
                            l_f = REAL(ABS(l_boundary - real(NumberXNodes, kind = real64) - 1))
                        CASE(4)
                            delIdx = delIdx + 1
                            if (l_boundary == 1) then
                                particleList(j)%energyLoss(1, iThread) = particleList(j)%energyLoss(1, iThread) + (SUM(v_f**2))!J/m^2 in 1D
                                particleList(j)%wallLoss(1, iThread) = particleList(j)%wallLoss(1, iThread) + 1 !C/m^2 in 1D
                            else if (l_boundary == NumberXNodes+1) then
                                particleList(j)%energyLoss(2, iThread) = particleList(j)%energyLoss(2, iThread) + (SUM(v_f**2)) !J/m^2 in 1D
                                particleList(j)%wallLoss(2, iThread) = particleList(j)%wallLoss(2, iThread) + 1 !C/m^2 in 1D
                            end if
                            exit
                        CASE default
                            print *, "l_sub is:", l_sub
                            print *, 'l_f is:', l_f
                            print *, world%boundaryConditions(INT(l_f))
                            print *, "Case does not exist in ongoing substep, depositJ"
                            stop
                        END SELECT
                    end if
                    
                    timePassed = timePassed + del_tau
                    if (timePassed >= del_t) then
                        particleList(j)%phaseSpace(1, i-delIdx, iThread) = l_f
                        particleList(j)%phaseSpace(2:4,i-delIdx, iThread) = v_f
                        exit
                    end if
                    ! now final position/velocity becomes next starting position/velocity
                    l_sub = l_f
                    v_sub = v_f
                    AtBoundaryBool = FutureAtBoundaryBool
                    if ((l_f < 1) .or. (l_f > NumberXNodes+1)) then
                        print *, 'l_f is:', l_f
                        stop "Have particles travelling outside the domain in depositJ!"
                    end if
                end do
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