program mtTest


    use iso_fortran_env
    
    implicit none
    

    real(real64), allocatable :: x(:), y(:), u(:)
    real(real64) :: a, count
    integer(int32) :: i, startTime, endTime, timingRate
    integer, parameter :: n = 10**8
    
    call system_clock(count_rate = timingRate) 
    allocate(x(n), y(n), u(n))
    call random_number(x)
    call random_number(y)
    a = 50.4d0
    
    


    
    call system_clock(startTime)
    do i = 1, n
        count = COS(x(i))
        u(i) = a * count + 0.5d0 * y(i)
    end do
    call system_clock(endTime)
    print *, 'Time cosine is:', real(endTime - startTime)/real(timingRate)

    call system_clock(startTime)
    do i = 1, n
        count = x(i)**2
        u(i) = a * count + 0.5d0 * y(i)
    end do
    call system_clock(endTime)
    print *, 'Time is square is:', real(endTime - startTime)/real(timingRate)

    
    
    

    
 
    


end program mtTest
