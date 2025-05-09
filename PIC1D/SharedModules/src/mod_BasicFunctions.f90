module mod_BasicFunctions

    ! Module containing basic functions and subroutines on arrays, etc sdf
    use iso_fortran_env, only: int32, real64, int64, output_unit
    use constants
    use omp_lib
    use iso_c_binding
    implicit none

    integer(int32), allocatable :: stateRan0(:) !Global variable that will need to be changed in-out for threaded random number generator
    integer(c_int64_t), allocatable :: state_PCG(:) !Random number save state for larger period numbers from PCG

    interface 
        ! Interface for PCG generator run using C
        real(c_double) function pcg32_random_r(state) bind(c)
        use iso_c_binding
        integer(c_int64_t) :: state
        end function
    end interface

contains

    !-------------------------- Basic plasma properties we may want to calculate----------------------------------------
    pure function getDebyeLength(T_e, n_e) result(res)
        real(real64), intent(in) :: T_e, n_e
        real(real64) :: res
        res = SQRT(eps_0 * T_e / n_e / e)
    end function getDebyeLength

    pure function getPlasmaFreq(n_e) result(res)
        real(real64), intent(in) :: n_e
        real(real64) :: res
        res = SQRT(n_e * (e**2) / m_e / eps_0)
    end function getPlasmaFreq

    pure function getMaxwellDistVx(v_x, T, m) result(res)
        ! maxwell distribution along one dimension in velocity
        ! T in eV, v_x and m in SI
        real(real64), intent(in) :: v_x(:), T, m
        real(real64) :: res(size(v_x))
        res = SQRT(m/2/pi/e/T) * EXP(- m * (v_x**2) /2 /e / T)
    end function getMaxwellDistVx

    pure function getMaxwellDistE(E_p, T) result(res)
        ! maxwell distribution along Energy (EEDF, 1/eV)
        ! T, E_p in eV and m in SI
        real(real64), intent(in) :: E_p(:), T
        real(real64) :: res(size(E_p))
        res = 2 * SQRT(E_p / pi) * (1/T)**(1.5) * EXP(- E_p/ T)
    end function getMaxwellDistE


    ! --------------- uniform Random Generators -----------------------------
    function ran2(irand) result(res)
        ! Function from from "Numerical recipes in fortran 77", basic rng, small period for 2^9 about
        ! Called ran0 actually in this book but already used ran2 so keep it
        integer(int32),parameter :: ia_ran2=16807,im_ran2=2147483647,iq_ran2=127773,ir_ran2=2836
        real(real64),parameter :: am_ran2=1.0d0/im_ran2
        integer(int32), intent(in out) :: irand
        real(real64):: res
        integer(int32) :: k
        
      
        k=irand/iq_ran2
        irand=ia_ran2*(irand-k*iq_ran2)-ir_ran2*k
        if (irand < 0) irand=irand+im_ran2
        res=am_ran2*irand
      end function ran2


    subroutine getRandom(x, irand)
        !use Gwenael's function to generate array of random number
        real(real64), intent(in out) :: x(:)
        integer(int32) :: i
        integer(c_int64_t), intent(in out) :: irand
        do i = 1, size(x)
            x(i) = pcg32_random_r(irand)
        end do
    end subroutine getRandom

    subroutine initializeRandomGenerators(numThread, stateRan0, state_PCG, preDeterminedBool)
        ! Initialize random number states for each thread
        integer(int32), allocatable, intent(out) :: stateRan0(:)
        integer(c_int64_t), allocatable, intent(out) :: state_PCG(:)
        integer(int32), intent(in) :: numThread
        logical, intent(in) :: preDeterminedBool
        real(real64) :: rando
        integer(int32) :: i
        
        allocate(stateRan0(numThread), state_PCG(numThread))
        if (.not. preDeterminedBool) call random_seed() ! Having seeding be randomly generated
        do i = 1, numThread
            call random_number(rando)
            stateRan0(i) = INT(rando * (huge(i)-1)) + 1 !12345 * i + 11 * (i-1) !
            call random_number(rando)
            state_PCG(i) = INT((rando-0.5d0) * 2 * (huge(state_PCG(i)-1)), kind = c_int64_t)
        end do
    end subroutine initializeRandomGenerators

    ! ------------------------------------------------

    subroutine get3DMaxwellianVelocity(V, mass, T, irand)
        ! T in eV
        real(real64), intent(in out) :: V(3)
        real(real64), intent(in) :: mass, t
        integer(c_int64_t), intent(in out) :: irand
        real(real64) :: U(4), v_therm
        integer(int32) :: i
        v_therm = SQRT(T*e/ mass)
        do i = 1, 4
            U(i) = pcg32_random_r(irand)
        end do
        V(1) = v_therm * SQRT(-2 * LOG(U(1))) * COS(2 * pi * U(2))
        V(2) = v_therm * SQRT(-2 * LOG(U(1))) * SIN(2 * pi * U(2))
        V(3) = v_therm * SQRT(-2 * LOG(U(3))) * SIN(2 * pi * U(4))
    end subroutine

    subroutine getMaxwellianSample(v, M, T_gas, irand)
        ! Use to sample a single neutral background particle at some temperature T_gas (eV)
        real(real64), intent(in out) :: v(3)
        real(real64), intent(in) :: T_gas, M
        integer(c_int64_t), intent(in out) :: irand
        real(real64) :: U(4), v_therm
        integer(int32) :: i
        v_therm = SQRT(T_gas*e/ M)
        do i = 1, 4
            U(i) = pcg32_random_r(irand)
        end do
        v(1) = v_therm * SQRT(-2.0d0 * LOG(U(1))) * COS(2.0d0 * pi * U(2))
        v(2) = v_therm * SQRT(-2.0d0 * LOG(U(1))) * SIN(2.0d0 * pi * U(2))
        v(3) = v_therm * SQRT(-2.0d0 * LOG(U(3))) * SIN(2.0d0 * pi * U(4))

    end subroutine getMaxwellianSample

    subroutine getMaxwellianFluxSample(v, M, T_gas, irand)
        ! Use to sample a single neutral background particle at some temperature T_gas (eV)
        real(real64), intent(in out) :: v(3)
        real(real64), intent(in) :: T_gas, M
        integer(c_int64_t), intent(in out) :: irand
        real(real64) :: U(4), v_therm
        integer(int32) :: i
        v_therm = SQRT(T_gas*e/ M)
        do i = 1, 4
            U(i) = pcg32_random_r(irand)
        end do
        v(1) = v_therm * SQRT(-2 * LOG(U(1))) * SIGN(1.0d0, U(2)-0.5d0)
        v(2) = v_therm * SQRT(-2 * LOG(U(3))) * COS(2 * pi * U(4))
        v(3) = v_therm * SQRT(-2 * LOG(U(3))) * SIN(2 * pi * U(4))

    end subroutine getMaxwellianFluxSample


    ! -------------------------General solver using arrays --------------------------------------

    subroutine solve_tridiag(n, diagLower, diagUpper, diag, b, x)
        ! General tridiagonal solver
        integer(int32), intent(in) :: n
        real(real64), intent(in out) :: x(n)
        real(real64), intent(in) :: diagLower(n-1), diagUpper(n-1), diag(n), b(n)
        integer(int32) :: i !n size dependent on how many points are inside (not boundary), so how many in rhs equation
        real(real64) :: m, cp(n-1),dp(n)

    ! initialize c-prime and d-prime
        cp(1) = diagUpper(1)/diag(1)
        dp(1) = b(1)/diag(1)
    ! solve for vectors c-prime and d-prime
        do i = 2,n-1
            m = diag(i)-cp(i-1)*diagLower(i-1)
            cp(i) = diagUpper(i)/m
            dp(i) = (b(i)-dp(i-1)*diagLower(i-1))/m
        end do
        dp(n) = (b(n)-dp(n-1)*diagLower(n-1))/(diag(n)-cp(n-1)*diagLower(n-1))
        x(n) = dp(n)
        do i = n-1, 1, -1
            x(i) = dp(i)-cp(i)*x(i+1)
        end do
    end subroutine solve_tridiag

    function triMul(n, diagLower, diagUpper, diag, x) result(res)
        ! for size nxn tridiagonal matrix with lower, upper and center diagonal matrix multiplied with x(n)
        integer(int32), intent(in) :: n
        real(real64), intent(in) :: diag(n), diagUpper(n-1), diagLower(n-1)
        real(real64), intent(in) :: x(n)
        integer(int32) :: i
        real(real64) :: res(n)
        res(1) = x(1) * diag(1) + x(2) * diagUpper(1)
        do i = 2, n-1
            res(i) = x(i) * diag(i) + x(i-1) * diagLower(i-1) + x(i+1) * diagUpper(i)
        end do
        res(n) = x(n) * diag(n) + x(n-1) * diagLower(n-1)
    end function triMul

    subroutine solveGaussElimination(A, b, n)
        ! solve system of equations (x = A^-1 b) using gaussian elimination
        ! return solution in b matrix, changes A matrix as well
        integer(int32), intent(in) :: n
        real(real64), intent(in out) :: A(n,n)
        real(real64), intent(in out) :: b(n)
        integer(int32) :: i, k, j
        real(real64) :: k1, k2
        do k = 1, n-1
            k1 = A(k,k)
            do i = k+1, n
                A(i,k) = A(i,k)/k1
            end do
            
            do j = k+1, n
                k2 = A(k,j)
                do i = k+1, n
                    A(i,j) = A(i,j) - A(i,k) * k2
                end do
            end do
        end do

        do i = 1, n
            do j = 1, i-1
                b(i) = b(i) - a(i,j) * b(j)
            end do
        end do

        do i = n, 1, -1
            do j = i+1, n
                b(i) = b(i) - a(i, j) * b(j)
            end do
    
            b(i) = b(i) / a(i, i)
        end do


    end subroutine solveGaussElimination

    function solveNormalEquation(A, b, m, n) result(b_new)
        ! normal equations for minimized A*x - b for matrix A(m,n), b(m)
        integer(int32), intent(in) :: m, n
        real(real64), intent(in) :: A(m,n), b(m)
        real(real64) :: A_new(n, n), b_new(n)
        b_new = MATMUL(TRANSPOSE(A), b)
        A_new = MATMUL(TRANSPOSE(A), A)
        call solveGaussElimination(A_new, b_new, n)

    end function solveNormalEquation

    ! double precision function normFunc(n, x, incx) result(res)
    !     integer(int32), intent(in) :: n, incx
    !     real(real64), intent(in) :: x(n)
    !     res = SQRT(SUM(x**2))
    ! end function normFunc

    ! double precision function innerProduct(n, x, incx, y, incy) result(res)
    !     integer(int32), intent(in) :: n, incx, incy
    !     real(real64), intent(in) :: x(n), y(n)
    !     res = SUM(x * y)
    ! end function innerProduct


    !------------------------ Array Functions -------------------------------------

    pure function arrayDiff(x, n) result(res)
        integer(int32), intent(in) :: n
        real(real64), intent(in) :: x(n)
        real(real64) :: res(n-1)
        res = x(2:) - x(1:n-1)
    end function arrayDiff

    pure function getArrayMean1D(x) result(res)
        real(real64), intent(in) :: x(:)
        real(real64) :: res
        res = sum(x)/size(x)

    end function getArrayMean1D

    pure function crossProduct(x, y) result(res)
        ! 3D vector cross product
        real(real64), intent(in) :: x(3), y(3)
        real(real64) :: res(3)
        res(1) = x(2) * y(3) - x(3) * y(2)
        res(2) = x(3) * y(1) - x(1) * y(3)
        res(3) = x(1) * y(2) - x(2) * y(1)
    end function crossProduct

    subroutine rotate2D(v, theta)
        ! rotate 2D vector by theta, theta going towards y axis
        real(real64), intent(in out) :: v(2)
        real(real64), intent(in) :: theta
        real(real64) :: res(2), cosTheta, sinTheta
        cosTheta = COS(theta)
        sinTheta = SIN(theta)
        res(1) = cosTheta * v(1) - sinTheta * v(2)
        res(2) = sinTheta * v(1) + cosTheta * v(2)
        v = res
    end subroutine rotate2D

    SUBROUTINE indexSortArray(n,arr,indx)
        !     ==============================================================
        !     VERSION:          /
        !     LAST MOD:      Years 1986-92
        !     MOD AUTHOR:    Numerical Recipes Software
        !     COMMENTS:         /
        !     --------------------------------------------------------------
          INTEGER:: n,indx(n),M,NSTACK
          REAL(KIND=8):: arr(n)
          PARAMETER (M=7,NSTACK=50)
          INTEGER:: i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
          REAL(KIND=8) a
        
          do  j=1,n
             indx(j)=j
          enddo
          jstack=0
          l=1
          ir=n
        1 if(ir-l.lt.M)then
             do j=l+1,ir
                indxt=indx(j)
                a=arr(indxt)
                do i=j-1,1,-1
                   if(arr(indx(i)).le.a)goto 2
                   indx(i+1)=indx(i)
                enddo
                i=0
        2       indx(i+1)=indxt
             enddo
             if(jstack.eq.0)return
             ir=istack(jstack)
             l=istack(jstack-1)
             jstack=jstack-2
          else
             k=(l+ir)/2
             itemp=indx(k)
             indx(k)=indx(l+1)
             indx(l+1)=itemp
             if(arr(indx(l+1)).gt.arr(indx(ir)))then
                itemp=indx(l+1)
                indx(l+1)=indx(ir)
                indx(ir)=itemp
             endif
             if(arr(indx(l)).gt.arr(indx(ir)))then
                itemp=indx(l)
                indx(l)=indx(ir)
                indx(ir)=itemp
             endif
             if(arr(indx(l+1)).gt.arr(indx(l)))then
                itemp=indx(l+1)
                indx(l+1)=indx(l)
                indx(l)=itemp
             endif
             i=l+1
             j=ir
             indxt=indx(l)
             a=arr(indxt)
        3    continue
             i=i+1
             if(arr(indx(i)).lt.a)goto 3
        4    continue
             j=j-1
             if(arr(indx(j)).gt.a)goto 4
             if(j.lt.i)goto 5
             itemp=indx(i)
             indx(i)=indx(j)
             indx(j)=itemp
             goto 3
        5    indx(l)=indx(j)
             indx(j)=indxt
             jstack=jstack+2
             if(jstack.gt.NSTACK) then 
                print*, 'NSTACK too small in indexx'
                stop
             endif
             if(ir-i+1.ge.j-l)then
                istack(jstack)=ir
                istack(jstack-1)=i
                ir=j-1
             else
                istack(jstack)=j-1
                istack(jstack-1)=l
                l=i
             endif
          endif
          goto 1
    END SUBROUTINE indexSortArray



end module mod_BasicFunctions