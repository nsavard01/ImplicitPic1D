module mod_potentialSolver
    use iso_fortran_env, only: int32, real64, output_unit
    use constants
    use mod_BasicFunctions
    use mod_domain
    use mod_Scheme
    implicit none

    ! This module will be for the potential solver
    ! This uses the interaction between particles (particle module) and grid (domain module) in order to solve for the potential
    ! Will contain particle to mesh gathers (J, rho) and potential which comes from them (phi)
    ! Will also contain particle mover, since needed to store to J, and cannot be separate
    ! Assume dirichlet boundaries at ends for now, so matrix solver can go down by two dimensions
    private
    public :: potentialSolver, readSolver

    type :: potentialSolver
        real(real64), allocatable :: phi(:), J(:), rho(:), EField(:) !phi_f is final phi, will likely need to store two arrays for phi, can't be avoided
        real(real64) :: rho_const, RF_rad_frequency, RF_half_amplitude
        real(real64), allocatable :: a_tri(:), b_tri(:), c_tri(:) !for thomas algorithm potential solver, a_tri is lower diagonal, b_tri middle, c_tri upper
        logical :: RF_bool


    contains
        procedure, public, pass(self) :: solve_tridiag_Poisson
        procedure, public, pass(self) :: getTotalPE
        procedure, public, pass(self) :: depositRho
        !procedure, public, pass(self) :: makeConstSourceTerm
        procedure, public, pass(self) :: solvePotential
        procedure, public, pass(self) :: getError_tridiag_Poisson
        procedure, public, pass(self) :: construct_diagMatrix
        procedure, public, pass(self) :: getChargeContinuityError
        procedure, public, pass(self) :: initialVRewind
        procedure, public, pass(self) :: getEField
        procedure, public, pass(self) :: getJdotE
        procedure, public, pass(self) :: getTotalKE
        procedure, public, pass(self) :: getTotalPE2
        procedure, public, pass(self) :: getTotalE
    end type

    interface potentialSolver
        module procedure :: potentialSolver_constructor
    end interface potentialSolver

contains

    
    type(potentialSolver) function potentialSolver_constructor(world, leftVoltage, rightVoltage, RF_frequency, timeCurrent) result(self)
        ! Construct domain object, initialize grid, dx_dl, and nodeVol.
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: leftVoltage, rightVoltage
        real(real64), intent(in) :: RF_frequency, timeCurrent
        allocate(self % J(NumberXHalfNodes), self%rho(NumberXNodes), self % phi(NumberXNodes), self%a_tri(NumberXHalfNodes), &
        self%b_tri(NumberXNodes), self%c_tri(NumberXHalfNodes), self%EField(NumberXHalfNodes))
        call construct_diagMatrix(self, world)
        self % rho = 0.0d0
        self % J = 0.0d0
        self % phi = 0.0d0
        self % rho_const = 0.0d0
        self%RF_half_amplitude = 0.0d0
        self%RF_rad_frequency = 2.0d0 * pi * RF_frequency
        if (world%boundaryConditions(1) == 1) self%phi(1) = leftVoltage
        if (world%boundaryConditions(1) == 4) then
            self%RF_half_amplitude = leftVoltage
            self%phi(1) = self%RF_half_amplitude * SIN(timeCurrent * self%RF_rad_frequency)
        end if
        if (world%boundaryConditions(NumberXNodes) == 1) self%phi(NumberXNodes) = rightVoltage
        if (world%boundaryConditions(NumberXNodes) == 4) then
            if (self%RF_half_amplitude /= 0) then
                print *, 'Half amplitude voltage for RF already set, have two RF boundaries!'
                stop
            else
                self%RF_half_amplitude = rightVoltage
                self%phi(NumberXNodes) = self%RF_half_amplitude * SIN(timeCurrent * self%RF_rad_frequency)
            end if
        end if
        if (world%boundaryConditions(1) == 3) then
            self%phi(1) = leftVoltage
            self%phi(NumberXNodes) = leftVoltage
        end if
        
        self%RF_bool = self%RF_rad_frequency > 0 .and. self%RF_half_amplitude > 0 .and. (world%boundaryConditions(1) == 4 .or. world%boundaryConditions(NumberXNodes) == 4)
        self%EField = 0.0d0
    end function potentialSolver_constructor

    function getTotalE(self, del_t, particleList, world) result(res)
        type(Domain), intent(in) :: world
        class(potentialSolver), intent(in out) :: self
        type(Particle), intent(in out) :: particleList(:)
        real(real64), intent(in) :: del_t
        real(real64) :: res, PE, KE
        res = 0.0d0
        PE = self%getTotalPE2(world)
        KE = self%getTotalKE(del_t, particleList, world)
        res = PE + KE
    end function getTotalE

    function getTotalKE(self, del_t, particleList, world) result(res)
        ! calculate total KE in Joules/m^2
        type(Domain), intent(in) :: world
        class(potentialSolver), intent(in out) :: self
        type(Particle), intent(in out) :: particleList(:)
        real(real64), intent(in) :: del_t
        real(real64) :: res, v_i, v_f, E_x
        integer(int32) :: iThread, i, j, l_cell
        res = 0.0d0
        !$OMP parallel private(i, j, v_i, v_f, E_x, l_cell, iThread) reduction(+:res)
        iThread = omp_get_thread_num() + 1
        loopSpecies: do j = 1, numberChargedParticles
            loopParticles: do i=1, particleList(j)%N_p(iThread)
                ! cell rounded down
                l_cell = INT(particleList(j)%phaseSpace(1,i,ithread))
                ! e field
                E_x = self%EField(l_cell)
                ! initial and final velocities
                v_i = particleList(j)%phaseSpace(2,i,iThread)
                v_f = v_i + particleList(j)%q_over_m * E_x * del_t
                ! average velocity ^2
                res = res + ((v_i + v_f)/2)**2 * particleList(j)%mass * 0.5d0 * particleList(j)%w_p
            end do loopParticles
        end do loopSpecies
        !$OMP end parallel
    end function getTotalKE

    function getJdotE(self, world, del_t) result(res)
        ! Get energy in electric fields, future true, then derive from phi_f, otherwise phi
        ! In 1D is J/m^2
        class(potentialSolver), intent(in) :: self
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t
        real(real64) :: res
        res = SUM(self%EField * world%dx_dl * self%J) * del_t
    end function getJdotE

    pure function getEField(self, l_p) result(res)
        !return EField of particle at logical position l_p (this is half distance), per particle since particle mover loop will be per particle
        class(potentialSolver), intent(in) :: self
        real(real64), intent(in) :: l_p
        integer(int32) :: l_left
        real(real64) :: res
        l_left = INT(l_p)
        res = self%EField(l_left)
    end function getEField


    subroutine initialVRewind(self, particleList, del_t)
        class(potentialSolver), intent(in out) :: self
        type(Particle), intent(in out) :: particleList(:)
        real(real64), intent(in) :: del_t
        real(real64) :: q_m_ratio
        integer(int32) :: j, i, iThread
        loopSpecies: do j = 1, numberChargedParticles
            q_m_ratio = particleList(j)%q/particleList(j)%mass
            !$OMP parallel private(iThread, i)
            iThread = omp_get_thread_num() + 1
            loopParticles: do i = 1, particleList(j)%N_p(iThread)
                particleList(j)%phaseSpace(2, i, iThread) = particleList(j)%phaseSpace(2, i, iThread) - 0.5d0 * (q_m_ratio) * self%getEField(particleList(j)%phaseSpace(1, i, iThread)) * del_t
            end do loopParticles
            !$OMP end parallel
        end do loopSpecies
    end subroutine initialVRewind

    function getChargeContinuityError(self, rho_i, world, del_t) result(chargeError)
        class(potentialSolver), intent(in) :: self
        real(real64), intent(in) :: del_t, rho_i(NumberXNodes)
        type(Domain), intent(in) :: world
        integer(int32) :: i, k
        real(real64) :: chargeError, del_Rho
        chargeError = 0.0d0
        k = 0
        do i = 1, NumberXNodes
            del_Rho = self%rho(i) - rho_i(i)
            if (del_Rho /= 0) then
                SELECT CASE (world%boundaryConditions(i))
                CASE(0)
                    chargeError = chargeError + (1.0d0 + del_t * (self%J(i) - self%J(i-1))/del_Rho)**2
                    k = k + 1
                CASE(1,4)
                    continue
                CASE(2)
                    if (i == 1) then
                        chargeError = chargeError + (1.0d0 + 2.0d0 * del_t * self%J(1)/del_Rho)**2
                    else
                        chargeError = chargeError + (1.0d0 - 2.0d0 * del_t * self%J(NumberXHalfNodes)/del_Rho)**2
                    end if
                    k = k + 1
                CASE(3)
                    if (i == 1) then
                        chargeError = chargeError + (1.0d0 + del_t * (self%J(1) - self%J(NumberXHalfNodes))/del_Rho)**2
                    end if
                END SELECT
            end if
        end do
        chargeError = SQRT(chargeError/k)
    end function getChargeContinuityError

    subroutine construct_diagMatrix(self, world)
        ! construct diagonal components for thomas algorithm
        ! Same matrix used for gauss' law and divergence of ampere
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        integer(int32) :: i
        self%a_tri = 0.0d0
        self%b_tri = 0.0d0
        self%c_tri = 0.0d0
        do i = 1, NumberXNodes
            SELECT CASE (world%boundaryConditions(i))
            CASE(0)
                if (i < NumberXNodes) then
                    self % c_tri(i) = 1.0d0/world%dx_dl(i)
                end if
                if (i > 1) then
                    self%a_tri(i-1) = 1.0d0/ world%dx_dl(i-1)
                end if
                self%b_tri(i) = - (1.0d0/world%dx_dl(i-1)  + 1.0d0/world%dx_dl(i))
            CASE(1,4)
                self%b_tri(i) = 1.0d0
            CASE(2)
                if (i == 1) then
                    self % c_tri(i) = 2.0d0/world%dx_dl(i)
                    self%b_tri(i) = - (2.0d0/world%dx_dl(i))
                else if (i == NumberXNodes) then
                    self % a_tri(i-1) = 2.0d0/ world%dx_dl(i-1)
                    self%b_tri(i) = - (2.0d0/world%dx_dl(i-1))
                else
                    print *, "Neumann boundary not on left or right most index!"
                    stop
                end if
            CASE(3)
                self%b_tri(i) = 1.0d0
            CASE default
                print *, "Error when constructing poisson matrix, inner nodes not plasma or neumann!"
            END SELECT
        end do

    end subroutine construct_diagMatrix

    subroutine depositRho(self, particleList, world) 
        class(potentialSolver), intent(in out) :: self
        type(Particle), intent(in out) :: particleList(numberChargedParticles)
        type(Domain), intent(in) :: world
        integer(int32) :: i, iThread, rightThreadIndx, leftThreadIndx
        !$OMP parallel private(iThread, i, rightThreadIndx, leftThreadIndx)
        iThread = omp_get_thread_num() + 1 
        leftThreadIndx = world%threadNodeIndx(1,iThread)
        rightThreadIndx = world%threadNodeIndx(2,iThread)
        do i = 1, numberChargedParticles
            particleList(i)%workspace(:,iThread) = 0.0d0
            call interpolateParticleToNodes(particleList(i), world, iThread) 
        end do
        !$OMP barrier
        ! Concatenate density array among threads
        particleList(1)%workspace(leftThreadIndx:rightThreadIndx, 1) = SUM(particleList(1)%workspace(leftThreadIndx:rightThreadIndx, :), DIM=2) * particleList(1)%q_times_wp
        do i = 2, numberChargedParticles
            particleList(1)%workspace(leftThreadIndx:rightThreadIndx, 1) = particleList(1)%workspace(leftThreadIndx:rightThreadIndx, 1) &
                + SUM(particleList(i)%workspace(leftThreadIndx:rightThreadIndx, :), DIM=2) * particleList(i)%q_times_wp
        end do
        !$OMP end parallel
        SELECT CASE (world%boundaryConditions(1))
        CASE(1,4)
            particleList(1)%workspace(1, 1) = 0.0d0
        CASE(2)
            particleList(1)%workspace(1,1) = 2.0d0 * particleList(1)%workspace(1,1)
        CASE(3)
            particleList(1)%workspace(1,1) = particleList(1)%workspace(1,1) + particleList(1)%workspace(NumberXNodes, 1)
        END SELECT
        SELECT CASE (world%boundaryConditions(NumberXNodes))
        CASE(1,4)
            particleList(1)%workspace(NumberXNodes, 1) = 0.0d0
        CASE(2)
            particleList(1)%workspace(NumberXNodes,1) = 2.0d0 * particleList(1)%workspace(NumberXNodes,1)
        CASE(3)
            particleList(1)%workspace(NumberXNodes,1) = particleList(1)%workspace(1,1)
        END SELECT
        if (world%gridSmoothBool) then
            do i = 1, NumberXNodes
                SELECT CASE(world%boundaryConditions(i))
                CASE(0)
                    self%rho(i) = 0.25d0 * (particleList(1)%workspace(i-1, 1) + 2.0d0 * particleList(1)%workspace(i, 1) + particleList(1)%workspace(i+1, 1))
                CASE(1,4)
                    self%rho(i) = 0.0d0
                CASE(2)
                    if (i ==1) then
                        self%rho(i) = 0.25d0 * (2.0d0 * particleList(1)%workspace(1,1) + 2.0d0 * particleList(1)%workspace(2,1))
                    else    
                        self%rho(i) = 0.25d0 * (2.0d0 * particleList(1)%workspace(NumberXNodes,1) + 2.0d0 * particleList(1)%workspace(NumberXHalfNodes,1))
                    end if
                CASE(3)
                    self%rho(i) = 0.25d0 * (particleList(1)%workspace(NumberXHalfNodes, 1) + 2.0d0 * particleList(1)%workspace(1, 1) + particleList(1)%workspace(2, 1))
                END SELECT
            end do
        else
            self%rho = particleList(1)%workspace(:, 1)
        end if
    end subroutine depositRho

    subroutine solve_tridiag_Poisson(self, world, timeCurrent)
        ! Tridiagonal (Thomas algorithm) solver for initial Poisson
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: timeCurrent
        integer(int32) :: i
        real(real64) :: m, d(NumberXNodes), cp(NumberXHalfNodes),dp(NumberXNodes)

        do i=1, NumberXNodes
            SELECT CASE (world%boundaryConditions(i))
            CASE(0,2)
                d(i) = (-self%rho(i) - self%rho_const) / eps_0
            CASE(1,3)
                d(i) = self%phi(i)
            CASE(4)
                d(i) = self%RF_half_amplitude * SIN(self%RF_rad_frequency * timeCurrent)
            END SELECT
        end do
    ! initialize c-prime and d-prime
        cp(1) = self%c_tri(1)/self%b_tri(1)
        dp(1) = d(1)/self%b_tri(1)
    ! solve for vectors c-prime and d-prime
        do i = 2,NumberXHalfNodes
            m = self%b_tri(i)-cp(i-1)*self%a_tri(i-1)
            cp(i) = self%c_tri(i)/m
            dp(i) = (d(i)-dp(i-1)*self%a_tri(i-1))/m
        end do
        m = self%b_tri(NumberXNodes)-cp(NumberXHalfNodes)*self%a_tri(NumberXHalfNodes)
        dp(NumberXNodes) = (d(NumberXNodes)-dp(NumberXHalfNodes)*self%a_tri(NumberXHalfNodes))/m
    ! initialize x
        self%phi = dp
    ! solve for x from the vectors c-prime and d-prime
        do i = NumberXHalfNodes, 1, -1
            self%phi(i) = dp(i)-cp(i)*self%phi(i+1)
        end do
    self%EField = (self%phi(1:NumberXHalfNodes) - self%phi(2:NumberXNodes)) / world%dx_dl
    !self%EField(1) = (self%phi(1) - self%phi(2))/world%dx_dl(1)
    !self%EField(2:NumberXHalfNodes) = (self%phi(1:NumberXNodes-2) - self%phi(3:NumberXNodes))/world%dx_dl/2.0
    end subroutine solve_tridiag_Poisson

    subroutine solvePotential(self, particleList, world, timeCurrent)
        ! Solve for initial potential
        class(potentialSolver), intent(in out) :: self
        type(Particle), intent(in out) :: particleList(numberChargedParticles)
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: timeCurrent
        call self%depositRho(particleList, world)
        call self%solve_tridiag_Poisson(world, timeCurrent)
        
    end subroutine solvePotential

    function getError_tridiag_Poisson(self, world) result(res)
        class(potentialSolver), intent(in) :: self
        type(Domain), intent(in) :: world
        integer(int32) :: i
        real(real64) :: Ax(NumberXNodes), res, tiny
        Ax = triMul(NumberXNodes, self%a_tri, self%c_tri, self%b_tri, self%phi)
        res = 0.0d0
        do i = 1, NumberXNodes
            SELECT CASE (world%boundaryConditions(i))
            CASE(0,2)
                if (self%rho(i) + self%rho_const /= 0) then
                    res = res + (Ax(i)*eps_0/(-self%rho(i) - self%rho_const) - 1.0d0)**2
                else
                    res = res + (Ax(i)*eps_0)**2
                end if
            CASE(1,3,4)
                if (self%phi(i) /= 0) then
                    res = res + ((Ax(i))/(self.phi(i)) - 1.0d0)**2
                else
                    res = res + (Ax(i))**2
                end if
            END SELECT
        end do
        res = SQRT(res/NumberXNodes)

    end function getError_tridiag_Poisson

    ! subroutine makeConstSourceTerm(self, world)
    !     ! Store constant source term for ampere into rho array
    !     class(potentialSolver), intent(in out) :: self
    !     type(Domain), intent(in) :: world
    !     integer(int32) :: i
    !     do i =1, NumberXNodes
    !         SELECT CASE (world%boundaryConditions(i))
    !         CASE(0)
    !             self%rho(i) = - (self%phi(i) - self%phi(i-1))/world%dx_dl(i-1) + (self%phi(i+1) - self%phi(i))/world%dx_dl(i)
    !         CASE(1,3,4)
    !             self%rho(i) = self%phi_f(i)
    !         CASE(2)
    !             if (i == 1) then
    !                 self%rho(i) = -2.0d0 * (self%phi(i) - self%phi(i+1))/world%dx_dl(i)
    !             else if (i == NumberXNodes) then
    !                 self%rho(i) = -2.0d0 * (self%phi(i) - self%phi(i-1))/world%dx_dl(i-1)
    !             end if
    !         END SELECT
    !     end do
    ! end subroutine makeConstSourceTerm

    function getTotalPE(self, world) result(res)
        ! Get energy in electric fields, future true, then derive from phi_f, otherwise phi
        ! In 1D is J/m^2
        class(potentialSolver), intent(in) :: self
        type(Domain), intent(in) :: world
        real(real64) :: res
        res = 0.5 * eps_0 * SUM(arrayDiff(self%phi, NumberXNodes)**2 / world%dx_dl)
    end function getTotalPE

    function getTotalPE2(self, world) result(res)
        ! Get energy in electric fields, future true, then derive from phi_f, otherwise phi
        ! In 1D is J/m^2
        class(potentialSolver), intent(in) :: self
        type(Domain), intent(in) :: world
        real(real64) :: res
        res = 0.25 * (self%rho(1) * self%phi(1) * world%dx_dl(1) + self%rho(NumberXNodes) * self%phi(NumberXNodes) * world%dx_dl(NumberXHalfNodes))
        res = res + 0.25 * SUM(self%rho(2:NumberXHalfNodes)* self%phi(2:NumberXHalfNodes) * (world%dx_dl(1:NumberXHalfNodes-1) + world%dx_dl(2:NumberXHalfNodes)))
    end function getTotalPE2

    ! ---------------------- read solver input ----------------------------

    subroutine readSolver(GeomFilename, solver, world, timeCurrent)
        type(Domain), intent(in) :: world
        type(potentialSolver), intent(in out) :: solver
        character(len=*), intent(in) :: GeomFilename
        real(real64) :: timeCurrent
        integer(int32) :: io, intTemp
        real(real64) :: leftVoltage, rightVoltage, RF_frequency, realTemp(20)
        real(real64), allocatable :: allocReal(:)

        print *, ""
        print *, "Reading solver inputs:"
        print *, "------------------"
        if (.not. restartBool) then
            open(10,file='../InputData/'//GeomFilename)
        else
            open(10,file=restartDirectory//'/InputData/'//GeomFilename)
        end if
        read(10, *, IOSTAT = io)
        read(10, *, IOSTAT = io)
        read(10, *, IOSTAT = io)
        read(10, *, IOSTAT = io)
        read(10, *, IOSTAT = io)
        read(10, *, IOSTAT = io) leftVoltage, rightVoltage
        read(10, *, IOSTAT = io) RF_frequency
        close(10)
        solver = potentialSolver(world, leftVoltage, rightVoltage, RF_frequency, timeCurrent)
        print *, 'Left voltage:', solver%phi(1)
        print *, 'Right voltage:', solver%phi(NumberXNodes)
        print *, 'RF frequency:', solver%RF_rad_frequency/2.0d0/pi
        print *, 'RF_half_amplitude:', solver%RF_half_amplitude
        print *, 'RF_rad_frequency:', solver%RF_rad_frequency
        print *, "------------------"
        print *, ""

    end subroutine readSolver


end module mod_potentialSolver