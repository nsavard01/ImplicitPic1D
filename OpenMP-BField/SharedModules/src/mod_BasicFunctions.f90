module mod_BasicFunctions

    ! Module containing basic functions and subroutines on arrays, etc sdf
    use iso_fortran_env, only: int32, real64, output_unit
    use constants
    implicit none


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

    function ran2(irand)
        ! Function from Gwenael for making random number
        integer,parameter :: ia=16807,im=2147483647,iq=127773,ir=2836
        real(real64),parameter :: am=1.0d0/im
        real(real64):: ran2
        integer(int32) :: k
        integer(int32), intent(in out) :: irand
      
        k=irand/iq
        irand=ia*(irand-k*iq)-ir*k
        if (irand < 0) irand=irand+im
        ran2=am*irand
        return
      end function ran2

    subroutine getRandom(x, irand)
        !use Gwenael's function to generate array of random number
        real(real64), intent(in out) :: x(:)
        integer(int32) :: i
        integer(int32), intent(in out) :: irand
        do i = 1, size(x)
            x(i) = ran2(irand)
        end do
    end subroutine getRandom

    subroutine get3DMaxwellianVelocity(V, mass, T, irand)
        ! T in eV
        real(real64), intent(in out) :: V(3)
        real(real64), intent(in) :: mass, t
        integer(int32), intent(in out) :: irand
        real(real64) :: U(4)
        call getRandom(U, irand)
        V(1) = SQRT(T*e/ mass) * SQRT(-2 * LOG(U(1))) * COS(2 * pi * U(2))
        V(2) = SQRT(T*e/ mass) * SQRT(-2 * LOG(U(1))) * SIN(2 * pi * U(2))
        V(3) = SQRT(T*e/ mass) * SQRT(-2 * LOG(U(3))) * SIN(2 * pi * U(4))
    end subroutine

    subroutine getMaxwellianSample(v, M, T_gas, irand)
        ! Use to sample a single neutral background particle at some temperature T_gas (eV)
        real(real64), intent(in out) :: v(3)
        real(real64), intent(in) :: T_gas, M
        integer(int32), intent(in out) :: irand
        real(real64) :: U(4)
        call getRandom(U, irand)
        v(1) = SQRT(T_gas*e/ M) * SQRT(-2 * LOG(U(1))) * COS(2 * pi * U(2))
        v(2) = SQRT(T_gas*e/ M) * SQRT(-2 * LOG(U(1))) * SIN(2 * pi * U(2))
        v(3) = SQRT(T_gas*e/ M) * SQRT(-2 * LOG(U(3))) * SIN(2 * pi * U(4))

    end subroutine getMaxwellianSample

    subroutine getMaxwellianFluxSample(v, M, T_gas, irand)
        ! Use to sample a single neutral background particle at some temperature T_gas (eV)
        real(real64), intent(in out) :: v(3)
        real(real64), intent(in) :: T_gas, M
        integer(int32), intent(in out) :: irand
        real(real64) :: U(4)
        call getRandom(U, irand)
        v(1) = SQRT(T_gas*e/ M) * SQRT(-2 * LOG(U(1))) * SIGN(1.0d0, U(2)-0.5d0)
        v(2) = SQRT(T_gas*e/ M) * SQRT(-2 * LOG(U(3))) * COS(2 * pi * U(4))
        v(3) = SQRT(T_gas*e/ M) * SQRT(-2 * LOG(U(3))) * SIN(2 * pi * U(4))

    end subroutine getMaxwellianFluxSample


    ! General solver using arrays 

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



end module mod_BasicFunctions