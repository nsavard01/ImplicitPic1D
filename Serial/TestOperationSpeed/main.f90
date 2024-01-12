program main
    use iso_fortran_env, only: int32, real64
    implicit none

    integer(int32) :: n = 1000000000
    integer(int32) :: i, tclock1, tclock2, clock_rate, irand = 1234567, tempInt, otherInt
    integer(int32), allocatable :: l_boundary(:), l_array(:)
    logical, allocatable :: bool(:)
    real(real64), allocatable :: x(:), v(:), a(:)
    real(real64) :: elapsed_time

    allocate(x(n), v(n), l_boundary(n), l_array(n), bool(n), a(n))
    ! initialize x array
    do i = 1, n
        x(i) = ran2(irand) + irand
        v(i) = ran2(irand) - 0.5d0
    end do

    call system_clock(count_rate = clock_rate)
    
    
    print *, ''
    

    call system_clock(tclock1)
    do i = 1, n
        bool(i) = ABS(x(i) - 0.5d0) < 0.5d0
    end do
    call system_clock(tclock2)
    elapsed_time = real(tclock2 - tclock1, kind = real64) / real(clock_rate, kind = real64)
    print *, "Elapsed time for comparison ABS(x-0.5) < 0.5d0", elapsed_time, "seconds"

    call system_clock(tclock1)
    do i = 1, n
        tempInt = INT(x(i))
        l_array(i) = tempInt + (INT(SIGN(1.0d0, v(i))) - 1)/2 + 1
    end do
    call system_clock(tclock2)
    elapsed_time = real(tclock2 - tclock1, kind = real64) / real(clock_rate, kind = real64)
    print *, "Elapsed time for l_boundary on arithmetic is", elapsed_time, "seconds"


    call system_clock(tclock1)
    do i = 1, n
        tempInt = INT(x(i))
        if (v(i) > 0.0) then
            l_boundary(i) = tempInt + 1
        else
            l_boundary(i) = tempInt
        end if
    end do
    call system_clock(tclock2)
    elapsed_time = real(tclock2 - tclock1, kind = real64) / real(clock_rate, kind = real64)
    print *, "Elapsed time for l_boundary on condition is:", elapsed_time, "seconds"

    call system_clock(tclock1)
    do i = 1, n
        bool(i) = (x(i) < 0.5d0)
    end do
    call system_clock(tclock2)
    elapsed_time = real(tclock2 - tclock1, kind = real64) / real(clock_rate, kind = real64)
    print *, "Elapsed time for comparison x< 0.5", elapsed_time, "seconds"

    call system_clock(tclock1)
    do i = 1, n
        bool(i) = MOD(x(i), 1.0d0) == 0.0d0
    end do
    call system_clock(tclock2)
    elapsed_time = real(tclock2 - tclock1, kind = real64) / real(clock_rate, kind = real64)
    print *, "Elapsed time for comparison mod(x, 1.0d0) == 0.0d0", elapsed_time, "seconds"

    call system_clock(tclock1)
    do i = 1, n
        a(i) = x(i) + v(i) * 3.2d0
    end do
    call system_clock(tclock2)
    elapsed_time = real(tclock2 - tclock1, kind = real64) / real(clock_rate, kind = real64)
    print *, "Elapsed time for multiplication", elapsed_time, "seconds"

    call system_clock(tclock1)
    do i = 1, n
        a(i) = x(i) + v(i)
    end do
    call system_clock(tclock2)
    elapsed_time = real(tclock2 - tclock1, kind = real64) / real(clock_rate, kind = real64)
    print *, "Elapsed time for addition", elapsed_time, "seconds"
    

    
contains

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

end program main