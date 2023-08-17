program MultiGrid

    use iso_fortran_env, only: int32, real64
    use constants
    use mod_BasicFunctions
    implicit none
    integer(int32) :: i, N_course
    real(real64) :: L, dx, n_e, dx_course, initialRes, res_C
    real(real64), allocatable :: phi(:), a_tri(:), b_tri(:), c_tri(:), phi_test(:), x(:), rho(:), x_course(:), phi_course(:), b(:), R(:), a_tri_C(:), b_tri_C(:), c_tri_C(:), R_C(:)

    NumberXNodes = 61
    n_e = 5.d14
    L = 0.05d0
    eps_r = 1.d-8
    allocate(phi(NumberXNodes), b_tri(NumberXNodes), a_tri(NumberXNodes-1), c_tri(NumberXNodes-1), phi_test(NumberXNodes), x(NumberXNodes), rho(NumberXNodes), b(NumberXNodes), R(NumberXNodes))
    rho = e * n_e
    rho(NumberXNodes) = rho(NumberXNodes)*0.5d0
    rho(1) = 0.0d0
    b = -rho/eps_0
    dx = L/real(NumberXNodes-1)
    do i = 1, NumberXNodes
        x(i) = dx * (i-1)
    end do
    print *, 'dx is:', dx
    c_tri(1) = 0.0d0
    b_tri(1) = 1.0d0
    do i = 2, NumberXNodes-1
        c_tri(i) = 1.0d0/dx**2
        a_tri(i-1) = 1.0d0/dx**2
        b_tri(i) = - 2.0d0/dx**2
    end do
    a_tri(NumberXNodes-1) = 1.0d0/dx**2
    b_tri(NumberXNodes) = -1.0d0/dx**2

    N_course = (NumberXNodes + 1)/2
    allocate(x_course(N_course), phi_course(N_course), a_tri_C(N_course-1), b_tri_C(N_course), c_tri_C(N_course-1), R_C(N_course))
    phi_course = 0.0d0
    dx_course = 2 * dx
    c_tri_C(1) = 0.0d0
    b_tri_C(1) = 1.0d0
    do i = 2, N_course-1
        c_tri_C(i) = 1.0d0/dx_course**2
        a_tri_C(i-1) = 1.0d0/dx_course**2
        b_tri_C(i) = - 2.0d0/dx_course**2
    end do
    a_tri_C(N_course-1) = 1.0d0/dx_course**2
    b_tri_C(N_course) = -1.0d0/dx_course**2

    phi_test = e*n_e * x * (2.0d0*L - x)/2.0d0 / eps_0
    phi = 0.0d0
    R = b - triMul(NumberXNodes, a_tri, c_tri, b_tri, phi) 
    initialRes = SQRT(SUM(R**2))
    res_C = initialRes
    do i = 1, 100
        call solve_gaussSiedel(phi, NumberXNodes, a_tri, b_tri, c_tri, b, res_C*1.d-2, 10000)
        R = b - triMul(NumberXNodes, a_tri, c_tri, b_tri, phi) 
        res_C = SQRT(SUM(R**2))
        print *, res_C
        if (res_C < eps_r * initialRes) then
            print *, 'converge in', i, 'iterations'
            print *, phi - phi_test
            exit
        end if
        
        !Transfer fine-grid residual to course grid

        call restriction(R_C, R, N_course)
        res_C = SQRT(SUM((triMul(N_course, a_tri_C, c_tri_C, b_tri_C, phi_course) - R_C)**2))
        call solve_gaussSiedel(phi_course, N_course, a_tri_C, b_tri_C, c_tri_C, R_C, res_C*1.d-2, 10000)
        call prolongation(phi, phi_course, N_course)
    end do
    

    
 
contains

subroutine prolongation(phi_f, phi_C, N_C)
    integer(int32), intent(in) :: N_C
    real(real64), intent(in out) :: phi_f(2*N_C - 1)
    real(real64), intent(in) :: phi_C(N_C)
    integer(int32) :: i
    phi_f(1) = phi_f(1) + phi_C(1)
    phi_f(2*N_C - 1) = phi_f(2*N_C - 1) + phi_C(N_C)
    do i = 2, N_C-1
        phi_f(2*i-2) = phi_f(2*i-2) + 0.5d0 * (phi_C(i-1) + phi_C(i))
        phi_f(2*i-1) = phi_f(2*i-1) + phi_C(i)
        phi_f(2*i) = phi_f(2*i) + 0.5d0 * (phi_C(i) + phi_C(i+1))
    end do

end subroutine prolongation

subroutine restriction(R_C, R_F, N_C)
    integer(int32), intent(in) :: N_C
    real(real64), intent(in out) :: R_C(N_C)
    real(real64), intent(in) :: R_F(2*N_C - 1)
    integer(int32) :: i
    R_C(1) = R_F(1)
    do i = 2, N_C-1
        R_C(i) = 0.25d0 * (R_F(2*i-2) + 2.0d0 * R_F(2*i-1) + R_F(2*i))
    end do
    R_C(N_C) = R_F(2*N_C - 1)
end subroutine

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
        x(1) = (b(1) - sigma)/b_tri(1)
        do j = 2, n-1
            sigma = c_tri(j)*x(j+1) + a_tri(j-1) * x(j-1)    
            x(j) = (b(j) - sigma)/b_tri(j)
        end do
        sigma = a_tri(n-1) * x(n-1)
        x(n) = (b(n) - sigma)/b_tri(n)
        Ax = triMul(n, a_tri, c_tri, b_tri, x) - b
        res = SQRT(SUM(Ax**2))
        if (res < e_tol) then
            converge = .true.
            exit
        end if
    end do
    if (.not. converge) then
        stop 'Max iterations reached Gauss-siedel solver!'
    end if
end subroutine solve_gaussSiedel

end program MultiGrid