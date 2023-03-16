program main
    use iso_fortran_env, only: int32, real64
    implicit none

    real(real64), parameter :: c = 299792458.0d0 ! m/s
    real(real64), parameter :: eps_0 = 8.8541878128d-12 !F/m
    real(real64), parameter :: m_p = 1.67262192369d-27 ! kg
    real(real64), parameter :: m_e = 9.1093837015d-31 ! kg
    real(real64), parameter :: k_B = 1.380649d-23 ! m^2 kg s^-2 K^-1
    real(real64), parameter :: e = 1.602176634d-19 ! C
    real(real64), parameter :: mu_0 = 1.25663706212d-6 ! m kg s^-2 A^-2
    real(real64), parameter :: pi = 4.0d0*atan(1.0d0) ! pi from atan
    integer(int32) :: num_grid_nodes = 64, ipiv, info, i, j, m = 3, index, m_k
    real(real64) :: phi(64), grid(64), rho(64), L_domain = 0.1d0, d(62), du(61), dl(61), b(62), delx, phi_Pic(62), A_pic(62,62), initialR, res
    real(real64) :: F_k(62, 4), u_k(62, 4), mat(62, 3), alpha(62)
    integer(int32) :: lwork=65, work(65)
    rho = e * 1.0d12
    d = 0.0d0
    dl = 0.5d0
    du = 0.5d0
    delx = L_domain/63.0d0
    phi_Pic = 0.0d0
    b = -rho(2:63) * delx**2 /eps_0/2.0d0
    A_pic = 0.0d0
    do i=1, 62
        do j=1, 62
            if (ABS(j-i) == 1) then
                A_pic(i,j) = 0.5d0
            end if
        end do
    end do
    initialR = SQRT(SUM((MATMUL(A_pic, phi_Pic) - b - phi_Pic)**2))
    print *, initialR
    do i = 1, 50000
        phi_Pic = MATMUL(A_pic, phi_Pic) - b
        res = SQRT(SUM((MATMUL(A_pic, phi_Pic) - b - phi_Pic)**2))
        if (res < initialR * 1.0d-8) then
            print *, "Picard converges after", i, "iterations"
            exit
        end if
    end do
    
    u_k(:,1) = 0.0d0
    u_k(:,2) = func(u_k(:,1), b, A_pic)
    F_k(:,1) = u_k(:,2) - u_k(:,1)
    initialR = SQRT(SUM(F_k(:,1)**2))
    do i = 1, 10000
        index = MODULO(i, m+1) + 1
        m_k = MIN(i, m)
        F_k(:, index) = func(u_k(:,index), b, A_pic) - u_k(:,index)
        if (SQRT(SUM(F_k(:,index)**2)) < initialR * 1.0d-8) then
            print *, 'Anderson acceleration converged at', i, 'iteration'
            print *, 'phi is:'
            print *, u_k(:,index)
            exit
        end if
        do j = 0, m_k-1
            mat(:,j+1) = F_k(:, MODULO(i - m_k + j, m+1) + 1) - F_k(:, index)
        end do
        alpha = -F_k(:, index)
        call dgels('N', 62, m_k, 1, mat(:, 1:m_k), 62, alpha, 62, work, lwork, info)
        alpha(m_k+1) = 1.0d0 - SUM(alpha(1:m_k))
        u_k(:, MODULO(i+1, m+1) + 1) = alpha(1) * (F_k(:, MODULO(i-m_k, m+1) + 1) + u_k(:, MODULO(i-m_k, m+1) + 1))
        do j=1, m_k
            u_k(:, MODULO(i+1, m+1) + 1) = u_k(:, MODULO(i+1, m+1) + 1) + alpha(j + 1) * (F_k(:, MODULO(i-m_k + j, m+1) + 1) + u_k(:, MODULO(i-m_k + j, m+1) + 1))
        end do

    end do
    

    
    


    

contains
    function func(phi, b, A_pic) result(res)
        real(real64), intent(in) :: phi(62), b(62), A_pic(62,62)
        real(real64) :: res(62)
        res = MATMUL(A_pic, phi) - b
    end function func

end program main