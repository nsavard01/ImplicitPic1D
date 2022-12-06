module mod_BasicFunctions

    ! Module containing basic functions and subroutines on arrays, etc
    use iso_fortran_env, only: int32, real64, output_unit
    use constants
    implicit none


contains

    ! Basic plasma properties we may want to calculate
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

    pure real function getArrayMean1D(x) result(res)
        real(real64), intent(in) :: x(:)
        res = sum(x)/size(x)

    end function getArrayMean1D



end module mod_BasicFunctions