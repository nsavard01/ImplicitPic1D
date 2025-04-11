module mod_potentialSolver
    use iso_fortran_env, only: int32, real64, output_unit
    use constants
    use mod_BasicFunctions
    use mod_particle
    use mod_domain
    use omp_lib
    implicit none

    ! Class for scalar and field quantities on the domain
    public :: potentialSolver, readSolver

    type :: potentialSolver
        real(real64), allocatable :: phi(:), rho(:), EField(:) ! Store phi, rho, and EField
        real(real64) :: rho_const, RF_rad_frequency, RF_half_amplitude
        real(real64), allocatable :: a_tri(:), b_tri(:), c_tri(:) !for thomas algorithm potential solver, a_tri is lower diagonal, b_tri middle, c_tri upper


    contains
        procedure, public, pass(self) :: depositRho
        procedure, public, pass(self) :: solve_tridiag_Poisson
        procedure, public, pass(self) :: makeEField
        procedure, public, pass(self) :: solvePotential
        procedure, public, pass(self) :: construct_diagMatrix
        procedure, public, pass(self) :: initialVRewind
        procedure, public, pass(self) :: getTotalPE
        procedure, public, pass(self) :: getTotalE
        procedure, public, pass(self) :: moveParticles
    end type

    interface potentialSolver
        module procedure :: potentialSolver_constructor
    end interface potentialSolver
   
contains

    type(potentialSolver) function potentialSolver_constructor(world, leftVoltage, rightVoltage, RF_frequency) result(self)
        ! Construct object
        real(real64), intent(in) :: leftVoltage, rightVoltage
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: RF_frequency
        real(real64) :: angle_rad
        allocate(self % rho(NumberXNodes), self % phi(NumberXNodes), self%EField(NumberXNodes-1), self%a_tri(NumberXNodes-1), &
        self%b_tri(NumberXNodes), self%c_tri(NumberXNodes-1))
        self % a_tri = 0.0d0
        self % c_tri = 0.0d0
        self % b_tri = 0.0d0
        self % rho = 0.0d0
        self % rho_const = 0.0d0
        self % phi = 0.0d0
        self % EField = 0.0d0
        self%RF_half_amplitude = 0.0d0
        if (world%boundaryConditions(1) == 1) self%phi(1) = leftVoltage
        if (world%boundaryConditions(1) == 4) self%RF_half_amplitude = leftVoltage
        if (world%boundaryConditions(NumberXNodes) == 1) self%phi(NumberXNodes) = rightVoltage
        if (world%boundaryConditions(NumberXNodes) == 4) then
            if (self%RF_half_amplitude /= 0) then
                print *, 'Half amplitude voltage for RF already set, have two RF boundaries!'
                stop
            else
                self%RF_half_amplitude = rightVoltage
            end if
        end if
        if (world%boundaryConditions(1) == 3) then
            self%phi(1) = leftVoltage
            self%phi(NumberXNodes) = leftVoltage
        end if
        self%RF_rad_frequency = 2.0d0 * pi * RF_frequency

    end function potentialSolver_constructor

    subroutine depositRho(self, particleList, world) 
        ! calculate rho from particle locations
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        type(Particle), intent(in out) :: particleList(:)
        integer(int32) :: i, j, l_left, l_right, iThread, leftThreadIndx, rightThreadIndx
        real(real64) :: d
        self % rho = self%rho_const
        !$OMP parallel private(iThread, j, i, l_left, l_right, d, leftThreadIndx, rightThreadIndx)
        iThread = omp_get_thread_num() + 1
        do i=1, numberChargedParticles
            particleList(i)%workSpace(:, iThread) = 0.0d0   
            do j = 1, particleList(i)%N_p(iThread)
                l_left = INT(particleList(i)%phaseSpace(1, j, iThread))
                l_right = l_left + 1
                d = particleList(i)%phaseSpace(1, j, iThread) - real(l_left)
                particleList(i)%workSpace(l_left, iThread) = particleList(i)%workSpace(l_left, iThread) + (1.0d0-d)
                particleList(i)%workSpace(l_right, iThread) = particleList(i)%workSpace(l_right, iThread) + d
            end do        
        end do
        !$OMP barrier
        leftThreadIndx = world%threadNodeIndx(1,iThread)
        rightThreadIndx = world%threadNodeIndx(2,iThread)
        self%rho(leftThreadIndx:rightThreadIndx) = SUM(particleList(1)%workSpace(leftThreadIndx:rightThreadIndx, :), DIM=2) * particleList(1)%q_times_wp
        do i = 2, numberChargedParticles
            self%rho(leftThreadIndx:rightThreadIndx) = self%rho(leftThreadIndx:rightThreadIndx) &
                + SUM(particleList(i)%workSpace(leftThreadIndx:rightThreadIndx, :), DIM=2) * particleList(i)%q_times_wp
        end do
        !$OMP end parallel
    end subroutine depositRho

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

    subroutine solve_tridiag_Poisson(self, world, timeCurrent)
        ! Solve for phi using rho
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: timeCurrent
        integer(int32) :: i
        real(real64) :: m, d(NumberXNodes), cp(NumberXNodes-1),dp(NumberXNodes)

        do i=1, NumberXNodes
            SELECT CASE (world%boundaryConditions(i))
            CASE(0,2)
                d(i) = -self%rho(i) / eps_0
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
        do i = 2,NumberXNodes-1
            m = self%b_tri(i)-cp(i-1)*self%a_tri(i-1)
            cp(i) = self%c_tri(i)/m
            dp(i) = (d(i)-dp(i-1)*self%a_tri(i-1))/m
        end do
        m = self%b_tri(NumberXNodes)-cp(NumberXNodes-1)*self%a_tri(NumberXNodes-1)
        dp(NumberXNodes) = (d(NumberXNodes)-dp(NumberXNodes-1)*self%a_tri(NumberXNodes-1))/m
    ! initialize x
        self%phi = dp
    ! solve for x from the vectors c-prime and d-prime
        do i = NumberXNodes-1, 1, -1
            self%phi(i) = dp(i)-cp(i)*self%phi(i+1)
        end do

    end subroutine solve_tridiag_Poisson


    function getTotalPE(self, world) result(res)
        ! Get field energy in J/m^2
        class(potentialSolver), intent(in) :: self
        type(Domain), intent(in) :: world
        real(real64) :: res
        res = 0.5 * eps_0 * SUM(arrayDiff(self%phi, NumberXNodes)**2 / world%dx_dl)
        
    end function getTotalPE


    


    !----------------- Particle mover/sub-stepping procedures -------------------------------------------------------

    subroutine makeEField(self, world)
        ! Generate EField at nodes
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        self%EField = (self%phi(1:NumberXNodes-1) - self%phi(2:NumberXNodes))/world%dx_dl
    end subroutine makeEField

    ! -------------------------------------------- Particle mover without boolean checks for depositing J -----------------------------------------------------------

    subroutine initialVRewind(self, particleList, del_t)
        ! Rewind particle velocities by half step using E-Field
        class(potentialSolver), intent(in out) :: self
        type(Particle), intent(in out) :: particleList(:)
        real(real64), intent(in) :: del_t
        real(real64) :: q_m_ratio, E_Field
        integer(int32) :: j, i, iThread
        loopSpecies: do j = 1, numberChargedParticles
            q_m_ratio = particleList(j)%q/particleList(j)%mass
            !$OMP parallel private(iThread, i, E_Field)
            iThread = omp_get_thread_num() + 1
            loopParticles: do i = 1, particleList(j)%N_p(iThread)
                E_Field = self%EField(int(particleList(j)%phaseSpace(1, i, iThread)))
                particleList(j)%phaseSpace(2, i, iThread) = particleList(j)%phaseSpace(2, i, iThread) - 0.5d0 * (q_m_ratio) * E_Field * del_t
            end do loopParticles
            !$OMP end parallel
        end do loopSpecies
    end subroutine initialVRewind

    function getTotalE(self, particleList, world) result(res)
        ! Get total kinetic + potential energy using average velocity (v^k-1/2 + v^k+1/2) / 2
        class(potentialSolver), intent(in) :: self
        type(Domain), intent(in) :: world
        type(Particle), intent(in) :: particleList(:)
        real(real64) :: res, temp(numThread), PE, v, q_m_ratio, E_Field
        integer(int32) :: i, j, iThread
        PE = self%getTotalPE(world)
        temp = 0.0d0
        res = 0.0d0
        do j = 1, numberChargedParticles
            !$OMP parallel private(iThread, i, v, E_Field)
            iThread = omp_get_thread_num() + 1
            loopParticles: do i = 1, particleList(j)%N_p(iThread)
                E_Field = self%EField(int(particleList(j)%phaseSpace(1, i, iThread)))
                v = particleList(j)%phaseSpace(2, i, iThread) + 0.5d0 * particleList(j)%q_over_m * E_Field * del_t
                temp(iThread) = temp(iThread) + (v**2 + SUM(particleList(j)%phaseSpace(3:4, i, iThread)**2)) * particleList(j)%w_p * particleList(j)%mass * 0.5d0
            end do loopParticles
            !$OMP end parallel
        end do
        res = SUM(temp) + PE
    end function getTotalE

    subroutine moveParticles(self, particleList, world, del_t)
        ! particle mover using the fields
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        type(Particle), intent(in out) :: particleList(:)
        real(real64), intent(in) :: del_t
        integer(int32) :: j, i, delIdx, refIdx, iThread, N_p, l_cell, direction
        real(real64) :: v_prime, q_over_m, x_f, x_i, l_i
        !$OMP parallel private(iThread, i, j,delIdx, v_prime, x_f, l_cell, l_i, x_i, refIdx, N_p, direction)
        iThread = omp_get_thread_num() + 1
        loopSpecies: do j = 1, numberChargedParticles
            delIdx = 0
            refIdx = 0
            particleList(j)%wallLoss(:, iThread) = 0
            particleList(j)%energyLoss(:, iThread) = 0.0d0
            N_p = particleList(j)%N_p(iThread)
            loopParticles: do i = 1, N_p
                l_i = particleList(j)%phaseSpace(1, i, iThread)
                l_cell = int(l_i)
                ! First velocity change
                v_prime = particleList(j)%phaseSpace(2, i, iThread) + particleList(j)%q_over_m * self%EField(l_cell) * del_t 
                ! Get new position
                x_i = world%grid(l_cell) + world%dx_dl(l_cell) * (l_i - real(l_cell, kind = 8))
                x_f = x_i + v_prime * del_t
                ! Check if outside boundary
                if (x_f <= world%startX) then
                    SELECT CASE (world%boundaryConditions(1))
                    CASE(1,4)
                        particleList(j)%energyLoss(1, iThread) = particleList(j)%energyLoss(1, iThread) + v_prime**2 + SUM(particleList(j)%phaseSpace(3:4, i, iThread)**2)
                        particleList(j)%wallLoss(1, iThread) = particleList(j)%wallLoss(1, iThread) + 1 !C/m^2 in 1D
                        particleList(j)%momentumLoss(1,iThread) = particleList(j)%momentumLoss(1,iThread) + v_prime
                        delIdx = delIdx + 1
                    CASE(2)
                        x_f = 2.0d0 * world%startX - x_f
                        v_prime = -v_prime
                        refIdx = refIdx + 1
                        particleList(j)%refRecordIdx(refIdx, iThread) = i - delIdx
                    CASE(3)
                        x_f = MODULO(x_f - world%startX, world%L_domain) + world%startX
                    END SELECT
                else if ((x_f >= world%endX)) then
                    SELECT CASE (world%boundaryConditions(NumberXNodes))
                    CASE(1,4)
                        particleList(j)%energyLoss(2, iThread) = particleList(j)%energyLoss(2, iThread) + v_prime**2 + SUM(particleList(j)%phaseSpace(3:4, i, iThread)**2)
                        particleList(j)%wallLoss(2, iThread) = particleList(j)%wallLoss(2, iThread) + 1 !C/m^2 in 1D
                        particleList(j)%momentumLoss(2,iThread) = particleList(j)%momentumLoss(2,iThread) + v_prime
                        delIdx = delIdx + 1
                    CASE(2)
                        x_f = 2.0d0 * world%endX - x_f
                        v_prime = -v_prime
                        refIdx = refIdx + 1
                        particleList(j)%refRecordIdx(refIdx, iThread) = i - delIdx
                    CASE(3)
                        x_f = MODULO(x_f - world%startX, world%L_domain) + world%startX
                    END SELECT
                end if
                
                if (x_f > world%startX .and. x_f < world%endX) then  
                    direction = int(sign(1.0d0, x_f - x_i))
                    do while (x_f < world%grid(l_cell) .or. x_f >= world%grid(l_cell+1))
                        l_cell = l_cell + direction
                    end do
                    particleList(j)%phaseSpace(1, i-delIdx, iThread) = real(l_cell, kind = 8) + (x_f - world%grid(l_cell))/world%dx_dl(l_cell)
                    particleList(j)%phaseSpace(2, i-delIdx, iThread) = v_prime
                    particleList(j)%phaseSpace(3:4, i-delIdx, iThread) = particleList(j)%phaseSpace(3:4, i, iThread)
                end if
            end do loopParticles
            particleList(j)%N_p(iThread) = N_p - delIdx
            particleList(j)%delIdx(iThread) = delIdx
            particleList(j)%refIdx(iThread) = refIdx
        end do loopSpecies
        !$OMP end parallel
        ! Update particle accumulation stats
        do j = 1, numberChargedParticles
            particleList(j)%numToCollide = particleList(j)%N_p
            particleList(j)%accumEnergyLoss = particleList(j)%accumEnergyLoss + SUM(particleList(j)%energyLoss, DIM=2)
            particleList(j)%accumWallLoss = particleList(j)%accumWallLoss + SUM(particleList(j)%wallLoss, DIM=2)
        end do
    end subroutine moveParticles

    ! ---------------- Initial Poisson Solver -------------------------------------------------

    subroutine solvePotential(self, particleList, world, timeCurrent)
        ! Solve for initial potential
        class(potentialSolver), intent(in out) :: self
        type(Particle), intent(in out) :: particleList(:)
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: timeCurrent
        call self%depositRho(particleList, world)
        call self%solve_tridiag_Poisson(world, timeCurrent)
        call self%makeEField(world)
    end subroutine solvePotential


    ! ---------------- read Solver parameters --------------------------------------

    subroutine readSolver(GeomFilename,solver,world)
        ! Read inputs for solver class
        type(Domain), intent(in) :: world
        type(potentialSolver), intent(in out) :: solver
        character(len=*), intent(in) :: GeomFilename
        integer(int32) :: io
        real(real64) :: leftVoltage, rightVoltage, BFieldMag, angle, RF_frequency

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
        solver = potentialSolver(world, leftVoltage, rightVoltage, RF_frequency)
        print *, "Left voltage:", solver%phi(1)
        print *, "Right  voltage:", solver%phi(NumberXNodes)
        print *, 'RF frequency:', solver%RF_rad_frequency/2.0d0/pi
        print *, 'RF half amplitude:', solver%RF_half_amplitude
        print *, "------------------"
        print *, ""
        call solver%construct_diagMatrix(world)

    end subroutine readSolver


end module mod_potentialSolver