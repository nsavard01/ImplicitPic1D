program mtTest

    use mt19937_64
    use omp_lib
    use iso_fortran_env, only: output_unit, wp => real64, i4 => int32, i8 => int64
    
    implicit none
    

    type(mt19937), allocatable :: randGen(:)
    real(wp) :: r
    integer(i4) :: i, j, startTime, endTime, timingRate, iThread
    integer, parameter :: n = 10**6, numThread = 16
    integer(i4), allocatable :: irand(:)
    real(wp), allocatable :: x(:,:)
    call omp_set_num_threads(numThread)
    allocate(irand(numThread), randGen(numThread), x(n, numThread))
    call system_clock(count_rate = timingRate)
    do i = 1, numThread
        irand(i) = 12345*11 + i * 7
        call randGen(i)%initialize(12345*11 + i * 7)
    end do
    

    ! call random%initialize(42)


    ! do i = 1, 10
    !     r = random%genrand64_real1()
    !     write(output_unit, '(E30.16)') r
    ! end do
    
    ! call random%initialize(1776_i8)

    ! do i = 1, 10
    !     r = random%genrand64_real1()
    !     write(output_unit, '(E30.16)') r
    ! end do

    ! call random%initialize([1_i8,22_i8,333_i8])

    ! r = random%genrand64_real1()
    ! write(output_unit, '(E30.16)') r

    ! r = random%genrand64_real2()
    ! write(output_unit, '(E30.16)') r

    ! r = random%genrand64_real3()
    ! write(output_unit, '(E30.16)') r

    ! ! randomness tests:

    call system_clock(startTime)
    !$OMP parallel private(iThread, j,i) 
    iThread = omp_get_thread_num() + 1
    do j = 1, 10
        do i = 1, n
            x(i, iThread) = ran2(irand(iThread))
        end do
    end do
    !$OMP end parallel
    call system_clock(endTime)
    call print_results('ran2')
    print *, 'Time:', real(endTime - startTime)/real(timingRate)
    
    call system_clock(startTime)
    !$OMP parallel private(iThread, j,i) 
    iThread = omp_get_thread_num() + 1
    do j = 1, 10
        do i = 1, n
            x(i, iThread) = randGen(iThread)%genrand64_real1()
        end do
    end do
    !$OMP end parallel
    call system_clock(endTime)
    call print_results('genrand64_real1')
    print *, 'Time:', real(endTime - startTime)/real(timingRate)

    

    

    

    ! do i = 1, n
    !     x(i) = random%genrand64_real2()
    ! end do
    ! call print_results('genrand64_real2')

    ! do i = 1, n
    !     x(i) = random%genrand64_real3()
    ! end do
    ! call print_results('genrand64_real3')

    contains

        subroutine print_results(method)
            character(len=*),intent(in) :: method
            write(*,'(/A)') method//': test randomness '
            write(*,*) 'theory:', 0.5_wp,   1/12.0_wp
            write(*,*) 'actual:', sum(x(:,1))/n, sum((x(:,1)-0.5_wp)**2)/n
        end subroutine print_results

        function ran2(irand)
            ! Function from Gwenael for making random number
            integer(i4),parameter :: ia=16807,im=2147483647,iq=127773,ir=2836
            real(wp),parameter :: am=1.0d0/im
            real(wp):: ran2
            integer(i4) :: k
            integer(i4), intent(in out) :: irand
          
            k=irand/iq
            irand=ia*(irand-k*iq)-ir*k
            if (irand < 0) irand=irand+im
            ran2=am*irand
            return
        end function ran2

end program mtTest
