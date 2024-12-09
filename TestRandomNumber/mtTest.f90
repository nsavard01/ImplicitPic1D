program mtTest

    use mt19937_64
    use omp_lib
    use mod_Random
    use mod_PCG
    use iso_fortran_env, only: output_unit, wp => real64, i4 => int32, i8 => int64
    
    implicit none
    

    type(mt19937), allocatable :: randGen(:)
    type(randType), allocatable :: randOther(:)
    type(PCG_type), allocatable :: rand_PCG(:)
    real(wp) :: r, temp
    integer(i4) :: i, j, startTime, endTime, timingRate, iThread
    integer, parameter :: n = 10**7, numThread = 6, numBins = 200, outer_n = 100
    integer(int32) :: hist(numBins), thread_irand
    real(real64) :: var, mean
    real(c_double) :: test_c_double
    integer(c_int64_t) :: test_c_int
    real(real64), pointer :: point_x(:)
    integer(c_int64_t), allocatable, target :: state_PCG(:)
    integer(i4), allocatable, target :: irand(:)
    integer(int64), allocatable :: iStatePCG(:) 
    real(wp), allocatable, target :: x(:,:)
    

    call omp_set_num_threads(numThread)
    allocate(irand(numThread), randGen(numThread), x(n, numThread), randOther(numThread), iStatePCG(numThread), state_PCG(numThread), rand_PCG(numThread))
    call system_clock(count_rate = timingRate)
    call random_seed()
    do i = 1, numThread
        call random_number(r)
        irand(i) = INT(r * (huge(irand(i))) + 1)
        call random_number(r)
        iStatePCG(i) = INT((r-0.5d0) * (huge(iStatePCG(i))), kind = int64)
        call random_number(r)
        call randGen(i)%initialize(INT(r * (huge(irand(i))) + 1))
        call random_number(r)
        call randOther(i)%initialize(INT(r * (huge(irand(i))) + 1))
        call random_number(r)
        state_PCG(i) = INT((r-0.5d0) * (huge(state_PCG(i))), kind = c_int64_t)
        call random_number(r)
        call rand_PCG(i)%initialize(INT((r-0.5d0) * (huge(state_PCG(i))), kind = c_int64_t))
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

    var = 0
    mean = 0
    call system_clock(startTime)
    !$OMP parallel private(iThread, j,i) reduction(+:var, mean) 
    iThread = omp_get_thread_num() + 1
    do j = 1, outer_n
        do i = 1, n
            x(i, iThread) = randGen(iThread)%genrand64_real1()
        end do
        var = var + SUM((x(:, iThread) - 0.5d0)**2)
        mean = mean + SUM(x(:,iThread))
    end do
    !$OMP end parallel
    call system_clock(endTime)
    print *, 'results mt19973'
    print *, 'mean:', mean/(real(numThread) * n * outer_n)
    print *, 'var:', var/(real(numThread) * n * outer_n)
    print *, 'Time:', real(endTime - startTime)/real(timingRate)


    var = 0
    mean = 0
    print *, 'starting PCG:'
    print *, state_PCG
    call system_clock(startTime)
    !$OMP parallel private(iThread, j,i, test_c_int) reduction(+:var,mean)
    iThread = omp_get_thread_num() + 1
    test_c_int = state_PCG(iThread)
    do j = 1, outer_n
        do i = 1, n
            x(i, iThread) = pcg32_random_r(test_c_int)
        end do
        var = var + SUM((x(:,iThread) - 0.5d0)**2)
        mean = mean + SUM(x(:,iThread))
    end do
    state_PCG(iThread) = test_c_int
    !$OMP end parallel
    call system_clock(endTime)
    print *, 'results PCG'
    print *, 'mean:', mean/(real(numThread) * n * outer_n)
    print *, 'var:', var/(real(numThread) * n * outer_n)
    print *, 'Time:', real(endTime - startTime)/real(timingRate)
    print *, 'end PCG:'
    print *, state_PCG

    var = 0
    mean = 0
    print *, 'starting irand:'
    print *, irand
    call system_clock(startTime)
    !$OMP parallel private(iThread, j,i, thread_irand) reduction(+:var,mean)
    iThread = omp_get_thread_num() + 1
    thread_irand = irand(iThread)
    do j = 1, outer_n
        do i = 1, n
            x(i, iThread) = ran2(thread_irand)
        end do
        var = var + SUM((x(:,iThread) - 0.5d0)**2)
        mean = mean + SUM(x(:,iThread))
    end do
    irand(iThread) = thread_irand
    !$OMP end parallel
    call system_clock(endTime)
    print *, 'results ran0'
    print *, 'mean:', mean/(real(numThread) * n * outer_n)
    print *, 'var:', var/(real(numThread) * n * outer_n)
    print *, 'Time:', real(endTime - startTime)/real(timingRate)
    print *, 'end irand:'
    print *, irand

    var = 0
    mean = 0
    print *, 'starting irand type:'
    print *, irand
    call system_clock(startTime)
    !$OMP parallel private(iThread, j,i) reduction(+:var,mean)
    iThread = omp_get_thread_num() + 1
    do j = 1, outer_n
        do i = 1, n
            x(i, iThread) = randOther(iThread)%getRand()
        end do
        var = var + SUM((x(:,iThread) - 0.5d0)**2)
        mean = mean + SUM(x(:,iThread))
    end do
    !$OMP end parallel
    call system_clock(endTime)
    print *, 'results ran0 type'
    print *, 'mean:', mean/(real(numThread) * n * outer_n)
    print *, 'var:', var/(real(numThread) * n * outer_n)
    print *, 'Time:', real(endTime - startTime)/real(timingRate)




    





    
    

    
    
    stop
    ! var = 0
    ! mean = 0
    ! call system_clock(startTime)
    ! !$OMP parallel private(iThread, j,i) reduction(+:var, mean) 
    ! iThread = omp_get_thread_num() + 1
    ! do j = 1, outer_n
    !     do i = 1, n
    !         x(i, iThread) = getPCGRand(iStatePCG(iThread))
    !     end do
    !     var = var + SUM((x(:, iThread) - 0.5d0)**2)
    !     mean = mean + SUM(x(:,iThread))
    ! end do
    ! !$OMP end parallel
    ! call system_clock(endTime)
    ! print *, 'results PCG'
    ! print *, 'mean:', mean/(real(numThread) * n * outer_n)
    ! print *, 'var:', var/(real(numThread) * n * outer_n)
    ! print *, 'Time:', real(endTime - startTime)/real(timingRate)

    ! var = 0
    ! mean = 0
    ! call system_clock(startTime)
    ! !$OMP parallel private(iThread, j,i) reduction(+:var,mean)
    ! iThread = omp_get_thread_num() + 1
    ! do j = 1, outer_n
    !     do i = 1, n
    !         x(i, iThread) = randOther(iThread)%getRand()
    !     end do
    !     var = var + SUM((x(:, iThread) - 0.5d0)**2)
    !     mean = mean + SUM(x(:,iThread))
    ! end do
    ! !$OMP end parallel
    ! call system_clock(endTime)
    ! print *, 'results ranOther'
    ! print *, 'mean:', mean/(real(numThread) * n * outer_n)
    ! print *, 'var:', var/(real(numThread) * n * outer_n)
    ! print *, 'Time:', real(endTime - startTime)/real(timingRate)

    

    

    

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
