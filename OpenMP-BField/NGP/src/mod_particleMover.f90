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
    integer(int32), parameter :: maxPartIter = 50

contains

    subroutine subStepSolverAA(l_sub, v_sub, l_f, v_f, v_half, del_tau, timePassed, E_x, q_over_m, l_cell, BField, B_mag, dx_dl, f_tol, AtBoundaryBool, &
            del_t, l_boundary, numIter)
        ! Anderson Acceleration particle mover
        integer(int32), intent(in) :: l_cell
        real(real64), intent(in) :: E_x, q_over_m, BField(3), B_mag, dx_dl, f_tol, del_t
        real(real64), intent(in out) :: l_sub, v_sub(3), l_f, v_f(3), v_half(3), del_tau, timePassed
        logical, intent(in out) :: AtBoundaryBool
        integer(int32), intent(in out) :: l_boundary, numIter
        real(real64) :: coeffAccel, curr_Res, Res_k(m_Anderson_Particle+1), del_tau_k(m_Anderson_Particle+1), fitMat(m_Anderson_Particle), minCoeff(m_Anderson_Particle+1), del_tau_max
        integer(int32) :: k, m_k, u, index
        logical :: convergeBool
        
        del_tau_max = del_t - timePassed
        del_tau_k(1) = del_tau_max
        k = 0
        convergeBool = .false.
        do while (.not. convergeBool)
            index = MODULO(k, m_Anderson_Particle+1) + 1
            m_k = MIN(k, m_Anderson_Particle)
            coeffAccel = 0.5d0 * del_tau_k(index) * q_over_m
            v_half(1) = v_sub(1) + coeffAccel * E_x
            v_half(2:3) = v_sub(2:3)
            v_half = v_half + coeffAccel * (crossProduct(v_half, BField) + coeffAccel* SUM(v_half * BField) * BField)
            v_half = v_half / (1.0d0 + (coeffAccel*B_mag)**2)
            if (v_half(1) > 0) then
                l_boundary = l_cell + 1
            else
                l_boundary = l_cell
            end if
            del_tau = MIN(del_tau_max, (real(l_boundary) - l_sub) * dx_dl / v_half(1))
            Res_k(index) = del_tau - del_tau_k(index)
            curr_Res = ABS(Res_k(index))
            if (curr_Res < f_tol) then
                v_f = 2.0d0 * v_half - v_sub
                numIter = k+1
                AtBoundaryBool = (del_tau < del_tau_max)
                if(.not. AtBoundaryBool) then
                    l_f = l_sub + v_half(1) * del_tau / dx_dl
                    if (MOD(l_f, 1.0d0) == 0.0d0) then
                        AtBoundaryBool = .true.
                        l_boundary = NINT(l_f)
                    end if
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
                    print *, 'del_tau_max:', del_tau_max
                    stop
                end if
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

    end subroutine subStepSolverAA

    subroutine particleSubStepNoBField(l_sub, v_sub, l_f, v_half, del_tau, E_x, q_over_m, l_cell, dx_dl, FutureAtBoundaryBool, l_boundary)
        ! Do initial substep, where particles start between nodes
        integer(int32), intent(in) :: l_cell
        real(real64), intent(in) :: q_over_m, dx_dl, l_sub, v_sub(3), E_x
        real(real64), intent(in out) :: l_f, v_half(3), del_tau
        logical, intent(in out) :: FutureAtBoundaryBool
        integer(int32), intent(in out) :: l_boundary
        real(real64) :: c,a, b, del_tau_temp
        integer(int32) :: l_alongV, l_awayV
        
        ! Particle first between nodes, so solve quadratic for that particle depending on conditions
        a = 0.5d0 * q_over_m * E_x
        if (v_sub(1) > 0) then
            l_alongV = l_cell + 1
            l_awayV = l_cell
        else
            l_alongV = l_cell
            l_awayV = l_cell + 1
        end if
        c = (l_sub - real(l_alongV)) * dx_dl
        b = v_sub(1)**2 - 4.0d0*a*c
        if (b > 0) then
            ! v and a opposite direction, but particle can still reach boundary along v
            del_tau_temp = 2.0d0 * ABS(c)/(ABS(v_sub(1)) + SQRT(b))
            FutureAtBoundaryBool = (del_tau_temp < del_tau)
            if (FutureAtBoundaryBool) then
                del_tau = del_tau_temp
                l_boundary = l_alongV
            end if
        else
            ! v and a opposite direction, boundary opposite direction of v
            if (l_sub /= real(l_awayV)) then
                c = (l_sub - real(l_awayV)) * dx_dl
                del_tau_temp = 2.0d0 * ABS(c)/(SQRT(v_sub(1)**2 - 4.0d0*a*c) - ABS(v_sub(1)))
                FutureAtBoundaryBool = (del_tau_temp < del_tau) 
            else
            ! v and a opposite direction, reverses back to initial position
                del_tau_temp = ABS(v_sub(1))/ABS(a)
                FutureAtBoundaryBool = (del_tau_temp < del_tau)
            end if
            if (FutureAtBoundaryBool) then
                del_tau = del_tau_temp
                l_boundary = l_awayV
            end if
        end if
        if (del_tau < 0.0d0) then
            print *, 'del_tau < 0'
            stop
        end if
        v_half(1) = v_sub(1) + a * del_tau
        v_half(2:3) = v_sub(2:3)
        if (.not. FutureAtBoundaryBool) then
            l_f = l_sub + v_half(1) * del_tau / dx_dl
        else
            l_f = real(l_boundary)
        end if 
        if (l_f < real(l_cell) .or. l_f > real(l_cell + 1)) then
            print *, 'l_f out of bounds'
            stop
        end if
       
    end subroutine particleSubStepNoBField


    subroutine particleSubStepNoBFieldVStop(l_sub, v_sub, l_f, v_half, del_tau, E_x, q_over_m, l_cell, dx_dl, FutureAtBoundaryBool, l_boundary, f_tol)
        ! Do initial substep, where particles start between nodes
        integer(int32), intent(in) :: l_cell
        real(real64), intent(in) :: q_over_m, dx_dl, l_sub, v_sub(3), E_x, f_tol
        real(real64), intent(in out) :: l_f, v_half(3), del_tau
        logical, intent(in out) :: FutureAtBoundaryBool
        integer(int32), intent(in out) :: l_boundary
        real(real64) :: c,a, b, del_tau_temp
        integer(int32) :: l_alongV, l_awayV
        
        ! Particle first between nodes, so solve quadratic for that particle depending on conditions
        a = 0.5d0 * q_over_m * E_x
        if (v_sub(1) > 0) then
            l_alongV = l_cell + 1
            l_awayV = l_cell
        else
            l_alongV = l_cell
            l_awayV = l_cell + 1
        end if
        c = (l_sub - real(l_alongV)) * dx_dl
        b = v_sub(1)**2 - 4.0d0*a*c
        if (b > 0) then
            ! v and a opposite direction, but particle can still reach boundary along v
            del_tau_temp = 2.0d0 * ABS(c)/(ABS(v_sub(1)) + SQRT(b))
            FutureAtBoundaryBool = (del_tau_temp < del_tau)
            if (FutureAtBoundaryBool) then
                del_tau = del_tau_temp
                l_boundary = l_alongV
            end if
        else
            FutureAtBoundaryBool = .false.
            del_tau = MIN(- 0.5d0 * v_sub(1)/a + f_tol, del_tau)
        end if
        if (del_tau < 0.0d0) then
            print *, 'del_tau < 0'
            stop
        end if
        v_half(1) = v_sub(1) + a * del_tau
        v_half(2:3) = v_sub(2:3)
        if (.not. FutureAtBoundaryBool) then
            l_f = l_sub + v_half(1) * del_tau / dx_dl
        else
            l_f = real(l_boundary)
        end if 
        if (l_f < real(l_cell) .or. l_f > real(l_cell + 1)) then
            print *, 'l_f out of bounds'
            stop
        end if
        if (v_sub(1) > 0) then
            if (2.0d0 * v_half(1) - v_sub(1) < -1.0d0) then
                print *, 'issue with particle when v_sub positive'
            end if
        else
            if (2.0d0 * v_half(1) - v_sub(1) > 1.0d0) then
                print *, 'issue with particle when v_sub negative'
            end if
        end if
       
    end subroutine particleSubStepNoBFieldVStop



    subroutine subStepSolverGetPosition(l_sub, v_sub, l_f, v_half, del_tau, E_x, q_over_m, l_cell, BField, B_mag, dx_dl, FutureAtBoundaryBool, l_boundary)
        ! Anderson Acceleration particle mover
        integer(int32), intent(in) :: l_cell
        real(real64), intent(in) :: q_over_m, BField(3), B_mag, dx_dl, l_sub, v_sub(3), del_tau, E_x
        real(real64), intent(in out) :: l_f, v_half(3)
        logical, intent(in out) :: FutureAtBoundaryBool
        integer(int32), intent(in out) :: l_boundary
        real(real64) :: coeffAccel, Res_0, Res_1, del_tau_0, del_tau_1, del_tau_2, l_f_other
        integer(int32) :: l_boundary_test, k

        coeffAccel = 0.5d0 * del_tau * q_over_m
        v_half(1) = v_sub(1) + coeffAccel * E_x
        v_half(2:3) = v_sub(2:3)
        v_half = v_half + coeffAccel * (crossProduct(v_half, BField) + coeffAccel* SUM(v_half * BField) * BField)
        v_half = v_half / (1.0d0 + (coeffAccel*B_mag)**2)
        l_f = l_sub + v_half(1) * del_tau / dx_dl
        FutureAtBoundaryBool = (ABS(l_f - l_cell - 0.5d0) >= 0.5d0)
        if (FutureAtBoundaryBool) then
            if (l_f > l_cell + 1) then
                l_boundary = l_cell + 1
                l_f = real(l_boundary)
            else if (l_f < l_cell) then
                l_boundary = l_cell
                l_f = real(l_boundary)
            else
                FutureAtBoundaryBool = .false.
                if (l_f == l_cell + 1) then
                    l_f = real(l_cell) + 1.0d0 - 1.0d-12
                else if (l_f == l_cell) then
                    l_f = real(l_cell) + 1.0d-12
                else
                    print *, 'issue re-assigning particle outside boundary!'
                end if    
            end if
        end if
        
    end subroutine subStepSolverGetPosition

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
            print *, 'first res are not opposite in time Boundary!'
            print *, 'l_sub:', l_sub
            print *, 'l_f:', l_f
            print *, 'del_tau:', del_tau
            print *, 'v_sub:', v_sub
            print *, 'l_boundary:', l_boundary
            print *, 'v_half:', v_half
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
            print *, 'first res are not opposite self boundary!'
            print *, 'l_sub:', l_sub
            print *, 'del_tau:', del_tau
            print *, 'v_sub:', v_sub
            print *, 'l_boundary:', l_boundary
            print *, 'v_half:', v_half
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

    subroutine subStepSolverTotal(l_sub, v_sub, l_f, v_half, del_tau, q_over_m, E_x, BField, B_mag, dx_dl, f_tol, FutureAtBoundaryBool, l_cell, l_boundary, numIter, AtBoundaryBool)
        ! Anderson Acceleration particle mover
        real(real64), intent(in) :: q_over_m, l_sub, v_sub(3), BField(3), B_mag, dx_dl, E_x, f_tol
        real(real64), intent(in out) :: l_f, v_half(3), del_tau
        integer(int32), intent(in out) :: l_boundary, numIter
        integer(int32), intent(in) :: l_cell
        logical, intent(in) :: AtBoundaryBool
        logical, intent(in out) :: FutureAtBoundaryBool
        real(real64) :: coeffAccel, Res_a, Res_b, Res_c, Res_s, del_tau_a, del_tau_b, del_tau_c, del_tau_d, del_tau_s, tempVar
        integer(int32) :: k, u, l_boundary_up
        logical :: usedBiSection, dirVForwardBool, dirAccelForwardBool, convergeBool

        dirVForwardBool = (v_sub(1) > 0)
        tempVar = q_over_m * E_x 
        dirAccelForwardBool = (tempVar > 0)
        l_boundary_up = l_cell + 1

        if (dirVForwardBool) then
            l_boundary = l_boundary_up
        else
            l_boundary = l_cell
        end if

        del_tau_a = 0.0d0
        Res_a = l_sub - real(l_boundary)
        if (XOR_op(dirVForwardBool, dirAccelForwardBool)) then
            ! Accel and v_i opposite, take delta_tau from time to get v_f = 0
            del_tau_b = MIN(-v_sub(1) / tempVar, del_tau)
        else
            ! can take largest time step
            del_tau_b = del_tau
        end if
        ! print *, ''
        ! print *, 'start new loop'
        ! print *, 'del_tau:', del_tau
        ! print *, 'l_sub:', l_sub
        if (B_mag > 0) then
            del_tau_b = MIN(0.1d0/(ABS(q_over_m)*B_mag), del_tau_b)
        end if
        convergeBool = .false.
        ! Go forward with secant method until have res flip over boundary
        do k = 1, maxPartIter
            coeffAccel = 0.5d0 * del_tau_b * q_over_m
            v_half(1) = v_sub(1) + coeffAccel * E_x
            v_half(2:3) = v_sub(2:3)
            v_half = v_half + coeffAccel * (crossProduct(v_half, BField) + coeffAccel* SUM(v_half * BField) * BField)
            v_half = v_half / (1.0d0 + (coeffAccel*B_mag)**2)
            l_f = l_sub + v_half(1) * del_tau_b / dx_dl
            Res_b = l_f - real(l_boundary)
            ! print *, ''
            ! print *, 'del_tau_a:', del_tau_a
            ! print *, 'del_tau_b:', del_tau_b
            ! print *, 'l_f:', l_f
            ! print *, 'l_boundary:', l_boundary
            ! print *, 'Res_a:', Res_a
            ! print *, 'Res_b:', Res_b
            if (ABS(Res_b) < 1.0d-12) then
                FutureAtBoundaryBool = .true.
                numIter = k
                l_f = real(l_boundary)
                del_tau =  del_tau_b
                convergeBool = .true.
                exit
                ! print *, 'CONVERGED!'
            end if
            if (Res_a * Res_b < 0) then
                ! Sign flipped, cross boundary!
                ! print *, 'BOUNDARY CROSS!'
                FutureAtBoundaryBool = .true.
                numIter = k
                exit
            else
                if (ABS(Res_b) > ABS(Res_a)) then
                    ! going towards other boundary!
                        Res_a = Res_a + real(l_boundary)
                        if (l_boundary == l_boundary_up) then
                            l_boundary = l_cell
                        else
                            l_boundary = l_boundary_up
                        end if
                        Res_a = Res_a - real(l_boundary)
                        Res_b = l_f - real(l_boundary)
                        ! print *, 'Boundary flipped to:', l_boundary
                        ! print *, 'Res_a:', Res_a
                        ! print *, 'Res_b:', Res_b
                        if (Res_a * Res_b <= 0) then
                            ! Sign flipped w/respect to new boundary, cross boundary!
                            ! print *, 'BOUNDARY CROSS AFTER FLIP!'
                            FutureAtBoundaryBool = .true.
                            numIter = k
                            exit
                        end if
                end if
                del_tau_c = del_tau_b - Res_b * (del_tau_b - del_tau_a)/(Res_b - Res_a)
                if (del_tau_c >= del_tau) then
                    ! Solve for final del_tau
                    ! print *, 'FINAL DEL_TAU!'
                    coeffAccel = 0.5d0 * del_tau * q_over_m
                    v_half(1) = v_sub(1) + coeffAccel * E_x
                    v_half(2:3) = v_sub(2:3)
                    v_half = v_half + coeffAccel * (crossProduct(v_half, BField) + coeffAccel* SUM(v_half * BField) * BField)
                    v_half = v_half / (1.0d0 + (coeffAccel*B_mag)**2)
                    l_f = l_sub + v_half(1) * del_tau / dx_dl
                    ! print *, 'l_f:', l_f
                    if ((l_f < l_boundary_up) .and. (l_f > l_cell)) then
                        FutureAtBoundaryBool = .false.
                        convergeBool = .true.
                        ! print *, 'CONVERGED!'
                    else if (l_f >= l_boundary_up) then
                        ! print *, 'RIGHT BOUNDARY!'
                        Res_c = l_f - real(l_boundary_up)
                        convergeBool = (ABS(Res_c) < 1.0d-12)
                        if (.not. convergeBool) then
                            del_tau_a = del_tau_b
                            del_tau_b = del_tau
                            if (l_boundary == l_boundary_up) then
                                Res_a = Res_b
                            else
                                Res_a = Res_a + real(l_boundary)
                                l_boundary = l_boundary_up
                                Res_a = Res_a - real(l_boundary)
                            end if
                            Res_b = Res_c
                            FutureAtBoundaryBool = .true.
                        else
                            l_f = real(l_boundary_up) - 1.0d-12
                            FutureAtBoundaryBool = .false.
                        end if
                    else if (l_f <= l_cell) then
                        ! print *, 'LEFT BOUNDARY!'
                        Res_c = l_f - real(l_cell)
                        convergeBool = (ABS(Res_c) < 1.0d-12)
                        if (.not. convergeBool) then
                            del_tau_a = del_tau_b
                            del_tau_b = del_tau
                            if (l_boundary == l_cell) then
                                Res_a = Res_b
                            else
                                Res_a = Res_a + real(l_boundary)
                                l_boundary = l_cell
                                Res_a = Res_a - real(l_boundary)
                            end if
                            Res_b = Res_c
                            FutureAtBoundaryBool = .true.
                        else
                            l_f = real(l_cell) + 1.0d-12
                            FutureAtBoundaryBool = .false.
                        end if
                    else
                        print *, 'l_f not in possibilities!'
                        stop
                    end if
                    numIter = k + 1
                    exit
                end if
            end if
    
            
            del_tau_a = del_tau_b
            del_tau_b = del_tau_c
            Res_a = Res_b
        end do
        if (k > maxPartIter) then
            print *, 'maximum iterations for initial particle secant method!'
            print *, 'l_sub:', l_sub
            print *, 'v_sub:', v_sub
            print *, 'l_cell:', l_cell
            print *, 'del_tau:', del_tau
            stop
        end if
        
        if (.not. convergeBool) then
            if (.not. AtBoundaryBool .or. l_boundary /= INT(l_sub)) then
                ! Brent method for particle at boundary which is not at same location at initial position
                if (Res_a * Res_b >= 0) then
                    print *, 'first res are not opposite in time Boundary!'
                    print *, 'l_sub:', l_sub
                    print *, 'l_f:', l_f
                    print *, 'del_tau:', del_tau
                    print *, 'del_tau_a:', del_tau_a
                    print *, 'del_tau_b:', del_tau_b
                    print *, 'v_sub:', v_sub
                    print *, 'l_boundary:', l_boundary
                    print *, 'v_half:', v_half
                    print *, 'Res_a:', Res_a
                    print *, 'Res_b:', Res_b
                    stop
                end if

                if (ABS(Res_a) < ABS(Res_b) ) then
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
                        numIter = numIter + k -1
                        exit
                    end if
                    coeffAccel = 0.5d0 * del_tau_s * q_over_m
                    v_half(1) = v_sub(1) + coeffAccel * E_x
                    v_half(2:3) = v_sub(2:3)
                    v_half = v_half + coeffAccel * (crossProduct(v_half, BField) + coeffAccel* SUM(v_half * BField) * BField)
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
            else
                ! Brent method for particle coming back through initial boundary position
                
                del_tau_a = 0.0d0
                Res_a = v_sub(1)
                Res_b = v_half(1)

                if (Res_a * Res_b >= 0) then
                    print *, 'first res are not opposite self boundary!'
                    print *, 'l_sub:', l_sub
                    print *, 'del_tau:', del_tau
                    print *, 'v_sub:', v_sub
                    print *, 'l_boundary:', l_boundary
                    print *, 'v_half:', v_half
                    print *, 'Res_a:', Res_a
                    print *, 'Res_b:', Res_b
                    stop
                end if

                if (ABS(Res_a) < ABS(Res_b) ) then
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
                        numIter = numIter + k
                        exit
                    end if

                    coeffAccel = 0.5d0 * del_tau_s * q_over_m
                    v_half(1) = v_sub(1) + coeffAccel * E_x
                    v_half(2:3) = v_sub(2:3)
                    v_half = v_half + coeffAccel * (crossProduct(v_half, BField) + coeffAccel* SUM(v_half * BField) * BField)
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
            end if
        end if
        

    end subroutine subStepSolverTotal

    subroutine depositJ(solver, particleList, world, del_t)
        ! particle substepping procedure which deposits J, also used for general moving of particles
        ! boolDepositJ = true : deposit J
        ! boolDepositJ = false: only move particles (done after non-linear iteration), also delete particles for wall collisions (saves extra loop in collisions)
        class(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        type(Particle), intent(in out) :: particleList(:)
        real(real64), intent(in) :: del_t
        !a and c correspond to quadratic equations | l_alongV is nearest integer boundary along velocity component, away is opposite
        real(real64) :: l_f, l_sub, v_sub(3), v_f(3), timePassed, del_tau, q_over_m, f_tol, v_half(3), dx_dl, E_x, q_times_wp, l_f_other, del_tau_other, v_prime(3)
        integer(int32) :: j, i, l_cell, iThread, l_boundary, numIter, l_boundary_other
        logical :: AtBoundaryBool, FutureAtBoundaryBool, FutureAtBoundaryBool_other, timeNotConvergedBool
        !$OMP parallel private(iThread)
        iThread = omp_get_thread_num() + 1 
        solver%J(:, iThread) = 0.0d0
        !$OMP end parallel
        solver%EField = 0.5d0 * (solver%phi(1:NumberXNodes-1) + solver%phi_f(1:NumberXNodes-1) - solver%phi(2:NumberXNodes) - solver%phi_f(2:NumberXNodes)) / world%dx_dl
        f_tol = del_t * 1.d-10
        loopSpecies: do j = 1, numberChargedParticles
            q_over_m = particleList(j)%q/particleList(j)%mass
            q_times_wp = particleList(j)%q * particleList(j)%w_p
            !$OMP parallel private(iThread, i, l_f, l_sub, v_sub, v_f, v_half, timePassed, del_tau, l_cell, AtBoundaryBool, FutureAtBoundaryBool, dx_dl, E_x, l_boundary, numIter, timeNotConvergedBool)
            iThread = omp_get_thread_num() + 1 
            loopParticles: do i = 1, particleList(j)%N_p(iThread)
                v_sub = particleList(j)%phaseSpace(2:4,i,iThread)
                l_sub = particleList(j)%phaseSpace(1,i,iThread)
                timePassed = 0.0d0
                AtBoundaryBool = .false.
                timeNotConvergedBool = .true.
                do while(timeNotConvergedBool)
                    if (.not. AtBoundaryBool) then
                        l_cell = INT(l_sub)
                    else
                        l_cell = l_boundary + (INT(SIGN(1.0d0, v_sub(1))) - 1)/2
                    end if
                   
                    E_x = solver%EField(l_cell)
                    dx_dl = world%dx_dl(l_cell)
                    del_tau = del_t - timePassed
                    del_tau_other = del_t - timePassed
                    if (del_tau < f_tol) then
                        print *, 'del_tau less than f_tol!'
                        stop
                    end if

                    if (.not. solver%BFieldBool) then
                        call particleSubStepNoBField(l_sub, v_sub, l_f_other, v_prime, del_tau_other, E_x, q_over_m, l_cell, dx_dl, FutureAtBoundaryBool_other, l_boundary_other)
                        call subStepSolverTotal(l_sub, v_sub, l_f, v_half, del_tau, q_over_m, E_x, solver%BField, solver%BFieldMag, dx_dl, f_tol, &
                            FutureAtBoundaryBool, l_cell, l_boundary, numIter, AtBoundaryBool)
                    else
                        ! Start AA
                        call subStepSolverTotal(l_sub, v_sub, l_f, v_half, del_tau, q_over_m, E_x, solver%BField, solver%BFieldMag, dx_dl, f_tol, &
                            FutureAtBoundaryBool, l_cell, l_boundary, numIter, AtBoundaryBool)
                        ! call subStepSolverGetPosition(l_sub, v_sub, l_f, v_half, del_tau, E_x, q_over_m, l_cell, solver%BField, solver%BFieldMag, dx_dl, FutureAtBoundaryBool, l_boundary)
                        
                        ! if (FutureAtBoundaryBool) then
                        !     if (.not. AtBoundaryBool .or. l_boundary /= INT(l_sub)) then
                        !         call subStepSolverGetTimeBoundaryBrent(l_sub, v_sub, l_f, v_half, del_tau, q_over_m, E_x, solver%BField, solver%BFieldMag, dx_dl, f_tol, l_boundary, numIter)
                        !     else
                        !         call subStepSolverGetTimeSelfBoundaryBrent(l_sub, v_sub, v_half, del_tau, q_over_m, E_x, solver%BField, solver%BFieldMag, dx_dl, f_tol, l_boundary, numIter)
                        !     end if
                        ! end if
                    end if
                    
                    if (l_f - l_cell > 1 .or. l_f-l_cell < 0) then
                        print *, 'l_f outside cell'
                        print *, 'l_cell:', l_cell
                        print *, 'l_f:', l_f
                        stop
                    end if
                    if (ABS(l_f - l_f_other) > 1.0d-8) then
                        print *, 'discrepancy:'
                        print *, 'l_sub:', l_sub
                        print *, 'l_f:', l_f
                        print *, 'l_f_other:', l_f_other
                        stop
                    end if
                    if (ABS(del_tau - del_tau_other) > 10*f_tol) then
                        print *, 'del_tau is off'
                        print *, 'del_t - timePassed', del_t - timePassed
                        print *, 'del_tau - del_tau_other', del_tau - del_tau_other
                        print *, 'del_tau:', del_tau
                        print *, 'del_tau_other:', del_tau_other
                        print *, 'l_sub:', l_sub
                        print *, 'l_f:', l_f
                        print *, 'v_sub:', v_sub
                        print *, 'v_half:', v_half
                        print *, 'v_prime:', v_prime
                        print *, 'l_boundary:', l_boundary
                        stop
                    end if
                    if (FutureAtBoundaryBool_other /= FutureAtBoundaryBool) then
                        print *, 'Future at boundary boolean not same!'
                        stop
                    end if
                    
                    v_f = 2.0d0 * v_half - v_sub
                    v_prime = 2.0d0 * v_prime - v_sub
                    if (ABS((v_prime(1) - v_f(1))/v_f(1)) > 1.d-6) then
                        print *, 'mismatch in v_f(1)'
                        print *, 'l_sub:', l_sub
                        print *, 'l_f:', l_f
                        print *, 'v_sub:', v_sub
                        print *, 'v_f:', v_f
                        print *, 'v_prime:', v_prime
                        print *, 'per diff:', ABS((v_prime(1) - v_f(1))/v_f(1))
                        print *, 'l_boundary:', l_boundary
                        print *, 'del_tau:', del_tau
                        stop
                    end if
                    if (ABS((v_prime(2) - v_f(2))/v_f(2)) > 1.d-6) then
                        print *, 'l_sub:', l_sub
                        print *, 'l_f:', l_f
                        print *, 'v_sub:', v_sub
                        print *, 'l_boundary:', l_boundary
                        print *, 'del_tau:', del_tau
                        print *, 'mismatch in v_f(2)'
                        stop
                    end if

                    solver%J(l_cell, iThread) = solver%J(l_cell, iThread) + q_times_wp * (v_half(1))*del_tau/world%dx_dl(l_cell)/del_t
                    
                    if (FutureAtBoundaryBool) then
                        if (l_boundary_other /= l_boundary) then
                            print *, 'mismatch in l_boundary'
                            print *, 'l_boundary:', l_boundary
                            print *, 'l_boundary_other:', l_boundary_other
                            print *, 'l_sub:', l_sub
                            print *, 'l_f:', l_f
                            print *, 'l_f_other:', l_f
                            print *, 'del_tau:', del_tau
                            print *, 'del_tau_other:', del_tau_other
                            stop
                        end if
                        SELECT CASE (world%boundaryConditions(l_boundary))
                        CASE(0)
                            l_f = real(l_boundary)
                        CASE(1)
                            exit
                        CASE(2)
                            if (l_boundary == NumberXNodes) then
                                v_f(1) = -ABS(v_f(1))  
                            else
                                v_f(1) = ABS(v_f(1))
                            end if
                            l_f = real(l_boundary) 
                        CASE(3)
                            l_f = REAL(ABS(l_boundary - real(NumberXNodes, kind = real64) - 1))
                        CASE default
                            print *, "l_sub is:", l_sub
                            print *, 'l_f is:', l_f
                            print *, world%boundaryConditions(INT(l_f))
                            print *, "Case does not exist in ongoing substep, depositJ"
                            stop
                        END SELECT
                    end if
                    if ((l_f < 1) .or. (l_f > NumberXNodes)) then
                        print *, 'l_f is:', l_f
                        stop "Have particles travelling outside the domain in depositJ!"
                    end if
                    timePassed = timePassed + del_tau
                    timeNotConvergedBool = (del_t - timePassed > f_tol)
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
        real(real64) :: l_f, l_sub, v_sub(3), v_f(3), timePassed, del_tau, q_over_m, f_tol, v_half(3), dx_dl, E_x
        integer(int32) :: j, i, l_cell, iThread, delIdx, l_boundary, numIter, numSubStepAve, funcEvalCounter, refIdx
        logical :: AtBoundaryBool, refluxedBool, FutureAtBoundaryBool, timeNotConvergedBool
        solver%EField = 0.5d0 * (solver%phi(1:NumberXNodes-1) + solver%phi_f(1:NumberXNodes-1) - solver%phi(2:NumberXNodes) - solver%phi_f(2:NumberXNodes)) / world%dx_dl
        f_tol = del_t * 1.d-10
        loopSpecies: do j = 1, numberChargedParticles
            q_over_m = particleList(j)%q/particleList(j)%mass
            numSubStepAve = 0
            funcEvalCounter = 0
            !$OMP parallel private(iThread, i, l_f, l_sub, v_sub, v_f, v_half, timePassed, del_tau, l_cell, AtBoundaryBool, delIdx, dx_dl, E_x, l_boundary, numIter, &
                refIdx, refluxedBool, FutureAtBoundaryBool, timeNotConvergedBool) reduction(+:numSubStepAve, funcEvalCounter)
            iThread = omp_get_thread_num() + 1 
            delIdx = 0
            refIdx = 0
            particleList(j)%energyLoss(:, iThread) = 0.0d0
            particleList(j)%wallLoss(:, iThread) = 0.0d0
            loopParticles: do i = 1, particleList(j)%N_p(iThread)
                v_sub = particleList(j)%phaseSpace(2:4,i,iThread)
                l_sub = particleList(j)%phaseSpace(1,i,iThread)
                timePassed = 0.0d0
                refluxedBool = .false.
                AtBoundaryBool = .false.
                timeNotConvergedBool = .true.
                do while(timeNotConvergedBool)
                    numSubStepAve = numSubStepAve + 1
                    if (.not. AtBoundaryBool) then
                        l_cell = INT(l_sub)
                    else
                        l_cell = l_boundary + (INT(SIGN(1.0d0, v_sub(1))) - 1)/2
                    end if
                
                    E_x = solver%EField(l_cell)
                    dx_dl = world%dx_dl(l_cell)
                    del_tau = del_t - timePassed

                    if (.not. solver%BFieldBool) then
                        !call particleSubStepNoBField(l_sub, v_sub, l_f, v_half, del_tau, E_x, q_over_m, l_cell, dx_dl, FutureAtBoundaryBool, l_boundary)
                        call subStepSolverTotal(l_sub, v_sub, l_f, v_half, del_tau, q_over_m, E_x, solver%BField, solver%BFieldMag, dx_dl, f_tol, &
                            FutureAtBoundaryBool, l_cell, l_boundary, numIter, AtBoundaryBool)
                    else
                    
                        ! AA particle mover
                        call subStepSolverTotal(l_sub, v_sub, l_f, v_half, del_tau, q_over_m, E_x, solver%BField, solver%BFieldMag, dx_dl, f_tol, &
                            FutureAtBoundaryBool, l_cell, l_boundary, numIter, AtBoundaryBool)
                        ! call subStepSolverGetPosition(l_sub, v_sub, l_f, v_half, del_tau, E_x, q_over_m, l_cell, solver%BField, solver%BFieldMag, dx_dl, FutureAtBoundaryBool, l_boundary)
                        ! funcEvalCounter = funcEvalCounter + 1
                        ! if (FutureAtBoundaryBool) then
                        !     if (.not. AtBoundaryBool .or. l_boundary /= INT(l_sub)) then
                        !         call subStepSolverGetTimeBoundaryBrent(l_sub, v_sub, l_f, v_half, del_tau, q_over_m, E_x, solver%BField, solver%BFieldMag, dx_dl, f_tol, l_boundary, numIter)
                        !     else
                        !         call subStepSolverGetTimeSelfBoundaryBrent(l_sub, v_sub, v_half, del_tau, q_over_m, E_x, solver%BField, solver%BFieldMag, dx_dl, f_tol, l_boundary, numIter)
                        !     end if
                        !     funcEvalCounter = funcEvalCounter + numIter
                        ! end if
                    end if
                    
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
                            else if (l_boundary == NumberXNodes) then
                                particleList(j)%energyLoss(2, iThread) = particleList(j)%energyLoss(2, iThread) + (SUM(v_f**2)) !J/m^2 in 1D
                                particleList(j)%wallLoss(2, iThread) = particleList(j)%wallLoss(2, iThread) + 1 !C/m^2 in 1D
                            end if
                            exit
                        CASE(2)
                            if (l_boundary == NumberXNodes) then
                                v_f(1) = -ABS(v_f(1))  
                            else
                                v_f(1) = ABS(v_f(1))
                            end if
                            l_f = real(l_boundary) 
                            if (.not. refluxedBool) then
                                refIdx = refIdx + 1
                                particleList(j)%refRecordIdx(refIdx, iThread) = i - delIdx
                                refluxedBool = .true.
                            end if
                        CASE(3)
                            l_f = REAL(ABS(l_boundary - real(NumberXNodes, kind = real64) - 1))
                        CASE default
                            print *, "l_sub is:", l_sub
                            print *, 'l_f is:', l_f
                            print *, world%boundaryConditions(INT(l_f))
                            print *, "Case does not exist in ongoing substep, depositJ"
                            stop
                        END SELECT
                    end if
                    
                    timePassed = timePassed + del_tau
                    timeNotConvergedBool = (del_t - timePassed > f_tol)
                    
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