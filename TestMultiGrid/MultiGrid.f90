program MultiGrid

    use iso_fortran_env, only: int32, real64
    use constants
    use mod_BasicFunctions
    use mod_multiGridStorage
    implicit none
    integer(int32) :: i, N_C, totalSweepCount, numGrid_C, N_C_min, j, N_F, numSweeps, startTime, endTime, timingRate
    real(real64) :: L, dx, n_e, dx_C, initialRes, res_C
    real(real64), allocatable :: phi(:), a_tri(:), b_tri(:), c_tri(:), phi_test(:), x(:), rho(:), b(:), R_F(:)
    type(multiGridStorage), allocatable :: mGrid(:)

    call system_clock(count_rate = timingRate)
    NumberXNodes = 8193
    totalSweepCount = 0
    numGrid_C = 6
    numSweeps = 2
    N_C_min = 3
    call getMultiGridParam(numGrid_C, NumberXNodes, N_C_min)
    print *, 'Number of nodes', NumberXNodes
    print *, 'Number of coarse grids:', numGrid_C
    print *, 'Minimum coarse grid size:', N_C_min
    n_e = 5.d14
    L = 0.05d0
    eps_r = 1.d-6
    allocate(phi(NumberXNodes), b_tri(NumberXNodes), a_tri(NumberXNodes-1), c_tri(NumberXNodes-1), phi_test(NumberXNodes), x(NumberXNodes), rho(NumberXNodes), b(NumberXNodes), R_F(NumberXNodes), mGrid(numGrid_C))
    rho = e * n_e
    rho(1) = 0.0d0
    b = -rho/eps_0
    dx = L/real(NumberXNodes-1)
    do i = 1, NumberXNodes
        x(i) = dx * (i-1)
    end do
    c_tri(1) = 0.0d0
    b_tri(1) = 1.0d0
    do i = 2, NumberXNodes-1
        c_tri(i) = 1.0d0/dx**2
        a_tri(i-1) = 1.0d0/dx**2
        b_tri(i) = - 2.0d0/dx**2
    end do
    a_tri(NumberXNodes-1) = 2.0d0/dx**2
    b_tri(NumberXNodes) = -2.0d0/dx**2

    ! a_tri(NumberXNodes-1) = 1.0d0/dx**2
    ! b_tri(NumberXNodes) = -1.0d0/dx**2
    do i = 1, numGrid_C
        N_C = (NumberXNodes + (2**i - 1))/(2**i)
        mGrid(i) = multiGridStorage(N_C)
        print *, 'Creating multigrid with size', N_C
        dx_C = dx*(2**i)
        call mGrid(i)%constructPoissonEven(dx_C)
    end do
    
    phi_test = e*n_e * x * (2.0*L - x)/2.0d0 / eps_0
    totalSweepCount = 0
    phi = 0.0d0
    
    R_F = b - triMul(NumberXNodes, a_tri, c_tri, b_tri, phi)
    initialRes = SQRT(SUM(R_F**2))
    call system_clock(startTime)
    do i = 1, 1000000
        call sweep_gaussSiedel(phi, NumberXNodes, a_tri, b_tri, c_tri, b, numSweeps)
        R_F = b - triMul(NumberXNodes, a_tri, c_tri, b_tri, phi)
        res_C = SQRT(SUM((R_F**2)))
        if (res_C < eps_r * initialRes) then
            call system_clock(endTime)
            print *, 'converge in', i, 'iterations'
            print *, 'took', totalSweepCount, 'sweeps'
            print *, 'Took', real(endTime - startTime)/real(timingRate), 'seconds'
            exit
        end if
        ! Non-recursive 
        ! call restriction(mGrid(1)%R_C, R_F, mGrid(1)%N_x, NumberXNodes)
        ! mGrid(1)%phi_C = 0.0d0
        ! if (numGrid_C == 1) then
        !     ! call sweep_gaussSiedel(mGrid(1)%phi_C, mGrid(1)%N_x, mGrid(1)%a_tri, mGrid(1)%b_tri, &
        !     !     mGrid(1)%c_tri, mGrid(1)%R_C, numSweeps)
        !     call solve_tridiag(mGrid(1)%N_x, mGrid(1)%a_tri, mGrid(1)%c_tri, mGrid(1)%b_tri, mGrid(1)%R_C, mGrid(1)%phi_C)
        !     !call solve_gaussSiedel(mGrid(1)%phi_C, mGrid(1)%N_x, mGrid(1)%a_tri, mGrid(1)%b_tri, mGrid(1)%c_tri, mGrid(1)%R_C, eps_r, 10000000)
        ! else
        !     call V_Cycle_2(mGrid, numGrid_C, numSweeps)
        ! end if
        ! call prolongation(phi, mGrid(1)%phi_C, mGrid(1)%N_x, NumberXNodes)

        !Transfer fine-grid residual to course grid
        call restriction(mGrid(1)%R_C, R_F, mGrid(1)%N_x, NumberXNodes)
        mGrid(1)%phi_C = 0.0d0
        if (numGrid_C == 1) then
            ! call sweep_gaussSiedel(mGrid(1)%phi_C, mGrid(1)%N_x, mGrid(1)%a_tri, mGrid(1)%b_tri, &
            !     mGrid(1)%c_tri, mGrid(1)%R_C, 10)
            !call solve_tridiag(mGrid(1)%N_x, mGrid(1)%a_tri, mGrid(1)%c_tri, mGrid(1)%b_tri, mGrid(1)%R_C, mGrid(1)%phi_C)
            call solve_gaussSiedel(mGrid(1)%phi_C, mGrid(1)%N_x, mGrid(1)%a_tri, mGrid(1)%b_tri, mGrid(1)%c_tri, mGrid(1)%R_C, eps_r, 10000000)
        else
            call V_Cycle(mGrid, 1, numGrid_C, numSweeps)
        end if
        call prolongation(phi, mGrid(1)%phi_C, mGrid(1)%N_x, NumberXNodes)
        call sweep_gaussSiedel(phi, NumberXNodes, a_tri, b_tri, c_tri, b, numSweeps)
    end do
    

    
 
contains

subroutine getMultiGridParam(numGrid_C, NumberXNodes, N_C_min)
    integer(int32), intent(in out) :: numGrid_C, NumberXNodes, N_C_min
    integer(int32) :: temp
    temp = (NumberXNodes + (2** numGrid_C - 1))/(2**numGrid_C)
    if (temp < N_C_min) then
        stop "Number of coarse nodes is less than set minimum!"
    else
        if (MOD(temp, 2) /= 1) temp = temp + 1
        N_C_min = temp
    end if
    NumberXNodes = (2**numGrid_C) * N_C_min - (2**numGrid_C - 1)
end subroutine getMultiGridParam


subroutine prolongation(phi_f, phi_C, N_C, N_F)
    integer(int32), intent(in) :: N_C, N_F
    real(real64), intent(in out) :: phi_f(N_F)
    real(real64), intent(in) :: phi_C(N_C)
    integer(int32) :: i
    phi_f(1) = phi_f(1) + phi_C(1)
    do i = 2, N_C
        phi_f(2*i-2) = phi_f(2*i-2) + 0.5d0 * (phi_C(i-1) + phi_C(i))
        phi_f(2*i-1) = phi_f(2*i-1) + phi_C(i)
    end do
end subroutine prolongation

subroutine restriction(R_C, R_F, N_C, N_F)
    integer(int32), intent(in) :: N_C, N_F
    real(real64), intent(in out) :: R_C(N_C)
    real(real64), intent(in) :: R_F(N_F)
    integer(int32) :: i
    R_C(1) = R_F(1)
    do i = 2, N_C-1
        R_C(i) = 0.25d0 * (R_F(2*i-2) + 2.0d0 * R_F(2*i-1) + R_F(2*i))
    end do
    R_C(N_C) = 0.5d0*(R_F(N_F) + R_F(N_F-1))
end subroutine

recursive subroutine V_Cycle(multGrid, gridNum, numGridTot, numSweeps)
    type(multiGridStorage), intent(in out) :: multGrid(:)
    integer(int32), intent(in) :: gridNum, numGridTot, numSweeps
    ! print *, 'Starting sweep on', gridNum, 'grid'
    call sweep_gaussSiedel(multGrid(gridNum)%phi_C, multGrid(gridNum)%N_x, multGrid(gridNum)%a_tri, multGrid(gridNum)%b_tri, &
        multGrid(gridNum)%c_tri, multGrid(gridNum)%R_C, numSweeps)
    mGrid(gridNum)%res = mGrid(gridNum)%R_C - triMul(multGrid(gridNum)%N_x, multGrid(gridNum)%a_tri, multGrid(gridNum)%c_tri, &
        multGrid(gridNum)%b_tri, mGrid(gridNum)%phi_C)
    call restriction(mGrid(gridNum+1)%R_C, mGrid(gridNum)%res, multGrid(gridNum+1)%N_x, multGrid(gridNum)%N_x)
    multGrid(gridNum + 1)%phi_C = 0.0d0
    if (gridNum + 1 == numGridTot) then
        ! call sweep_gaussSiedel(multGrid(gridNum+1)%phi_C, multGrid(gridNum+1)%N_x, multGrid(gridNum+1)%a_tri, multGrid(gridNum+1)%b_tri, &
        !     multGrid(gridNum+1)%c_tri, multGrid(gridNum+1)%R_C, 3*numSweeps)
        !call solve_tridiag(mGrid(gridNum+1)%N_x, mGrid(gridNum+1)%a_tri, mGrid(gridNum+1)%c_tri, mGrid(gridNum+1)%b_tri, mGrid(gridNum+1)%R_C, mGrid(gridNum+1)%phi_C)
        call solve_gaussSiedel(mGrid(gridNum+1)%phi_C, mGrid(gridNum+1)%N_x, mGrid(gridNum+1)%a_tri, mGrid(gridNum+1)%b_tri, mGrid(gridNum+1)%c_tri, mGrid(gridNum+1)%R_C, eps_r, 10000000)
    else
        call V_Cycle(multGrid, gridNum + 1, numGridTot, numSweeps)
    end if
    call prolongation(multGrid(gridNum)%phi_C, multGrid(gridNum+1)%phi_C, multGrid(gridNum+1)%N_x, multGrid(gridNum)%N_x)
    call sweep_gaussSiedel(multGrid(gridNum)%phi_C, multGrid(gridNum)%N_x, multGrid(gridNum)%a_tri, multGrid(gridNum)%b_tri, &
        multGrid(gridNum)%c_tri, multGrid(gridNum)%R_C, numSweeps)
end subroutine V_Cycle

subroutine V_Cycle_2(multGrid, numGridTot, numSweeps)
    type(multiGridStorage), intent(in out) :: multGrid(:)
    integer(int32), intent(in) :: numGridTot, numSweeps
    integer(int32) :: gridNum
    ! print *, 'Starting sweep on', gridNum, 'grid'
    do gridNum = 1, numGridTot-1
        call sweep_gaussSiedel(multGrid(gridNum)%phi_C, multGrid(gridNum)%N_x, multGrid(gridNum)%a_tri, multGrid(gridNum)%b_tri, &
            multGrid(gridNum)%c_tri, multGrid(gridNum)%R_C, numSweeps)
        mGrid(gridNum)%res = mGrid(gridNum)%R_C - triMul(multGrid(gridNum)%N_x, multGrid(gridNum)%a_tri, multGrid(gridNum)%c_tri, &
            multGrid(gridNum)%b_tri, mGrid(gridNum)%phi_C)
        call restriction(mGrid(gridNum+1)%R_C, mGrid(gridNum)%res, multGrid(gridNum+1)%N_x, multGrid(gridNum)%N_x)
        multGrid(gridNum + 1)%phi_C = 0.0d0
    end do

    ! Last actual solve

    ! call sweep_gaussSiedel(multGrid(gridNum+1)%phi_C, multGrid(gridNum+1)%N_x, multGrid(gridNum+1)%a_tri, multGrid(gridNum+1)%b_tri, &
    !     multGrid(gridNum+1)%c_tri, multGrid(gridNum+1)%R_C, numSweeps)
    call solve_tridiag(mGrid(numGridTot)%N_x, mGrid(numGridTot)%a_tri, mGrid(numGridTot)%c_tri, mGrid(numGridTot)%b_tri, mGrid(numGridTot)%R_C, mGrid(numGridTot)%phi_C)
    !call solve_gaussSiedel(mGrid(gridNum+1)%phi_C, mGrid(gridNum+1)%N_x, mGrid(gridNum+1)%a_tri, mGrid(gridNum+1)%b_tri, mGrid(gridNum+1)%c_tri, mGrid(gridNum+1)%R_C, eps_r, 10000000)
    do gridNum = numGridTot-1, 1, -1
        call prolongation(multGrid(gridNum)%phi_C, multGrid(gridNum+1)%phi_C, multGrid(gridNum+1)%N_x, multGrid(gridNum)%N_x)
        call sweep_gaussSiedel(multGrid(gridNum)%phi_C, multGrid(gridNum)%N_x, multGrid(gridNum)%a_tri, multGrid(gridNum)%b_tri, &
            multGrid(gridNum)%c_tri, multGrid(gridNum)%R_C, numSweeps)
    end do
end subroutine V_Cycle_2

subroutine solve_gaussSiedel(x, n, a_tri, b_tri, c_tri, b, e_tol, maxIter)
    ! Tridiagonal (Thomas algorithm) solver for initial Poisson
    real(real64), intent(in out) :: x(n)
    integer(int32), intent(in) :: n, maxIter
    real(real64), intent(in) :: a_tri(n-1), b_tri(n), c_tri(n-1), b(n), e_tol
    real(real64) :: sigma, Ax(n), res
    integer(int32) :: i, j
    logical :: converge
    converge = .false.
    do i = 1, maxIter
        sigma = c_tri(1)*x(2)
        x(1) = x(1) + 1.4d0 * ((b(1) - sigma)/b_tri(1) - x(1))
        do j = 2, n-1
            sigma = c_tri(j)*x(j+1) + a_tri(j-1) * x(j-1)    
            x(j) = x(j) + 1.4d0 * ((b(j) - sigma)/b_tri(j) - x(j))
        end do
        sigma = a_tri(n-1) * x(n-1)
        x(n) = x(n) + 1.4d0 * ((b(n) - sigma)/b_tri(n) - x(n))
        Ax = (triMul(n, a_tri, c_tri, b_tri, x) + 1.d-15)/(b + 1.d-15) - 1.0d0
        res = SQRT(SUM(Ax**2)/n)
        if (res < e_tol) then
            converge = .true.
            totalSweepCount = totalSweepCount + i*n
            exit
        end if
    end do
    if (.not. converge) then
        print *, 'gauss Siedel has not converged!'
        stop
    end if
end subroutine solve_gaussSiedel

subroutine sweep_gaussSiedel(x, n, a_tri, b_tri, c_tri, b, sweepNumber)
    ! Tridiagonal (Thomas algorithm) solver for initial Poisson
    real(real64), intent(in out) :: x(n)
    integer(int32), intent(in) :: n, sweepNumber
    real(real64), intent(in) :: a_tri(n-1), b_tri(n), c_tri(n-1), b(n)
    real(real64) :: sigma
    integer(int32) :: i, j

    do i = 1, sweepNumber
        sigma = c_tri(1)*x(2)
        x(1) = (b(1) - sigma)/b_tri(1)
        do j = 2, n-1
            sigma = c_tri(j)*x(j+1) + a_tri(j-1) * x(j-1)    
            x(j) = (b(j) - sigma)/b_tri(j)
        end do
        sigma = a_tri(n-1) * x(n-1)
        x(n) = (b(n) - sigma)/b_tri(n)
    end do
    totalSweepCount = totalSweepCount + sweepNumber*n
end subroutine sweep_gaussSiedel

subroutine sweep_dampJacobi(x, n, a_tri, b_tri, c_tri, b, sweepNumber)
    ! Tridiagonal (Thomas algorithm) solver for initial Poisson
    real(real64), intent(in out) :: x(n)
    integer(int32), intent(in) :: n, sweepNumber
    real(real64), intent(in) :: a_tri(n-1), b_tri(n), c_tri(n-1), b(n)
    real(real64) :: sigma, temp(n)
    integer(int32) :: i, j

    do i = 1, sweepNumber
        sigma = c_tri(1)*x(2)
        temp(1) = x(1) + 0.66d0 * ((b(1) - sigma)/b_tri(1) - x(1))
        do j = 2, n-1
            sigma = c_tri(j)*x(j+1) + a_tri(j-1) * x(j-1)    
            temp(j) = x(j) + 0.66d0 * ((b(j) - sigma)/b_tri(j) - x(j))
        end do
        sigma = a_tri(n-1) * x(n-1)
        temp(n) = x(n) + 0.66d0 * ((b(n) - sigma)/b_tri(n) - x(n))
        x = temp
    end do
    totalSweepCount = totalSweepCount + sweepNumber*n
end subroutine sweep_dampJacobi

end program MultiGrid