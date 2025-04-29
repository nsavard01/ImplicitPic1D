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
        procedure, public, pass(self) :: solve_tridiag_Poisson
        procedure, public, pass(self) :: getEField
        procedure, public, pass(self) :: makeEField
        procedure, public, pass(self) :: solvePotential
        procedure, public, pass(self) :: construct_diagMatrix
        procedure, public, pass(self) :: getTotalPE
        procedure, public, pass(self) :: getTotalE
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
        allocate(self % rho(NumberXNodes), self % phi(NumberXNodes), self%EField(NumberXNodes), self%a_tri(NumberXNodes-1), &
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

    subroutine construct_diagMatrix(self, world)
        ! construct tridiagonal components of matrix for thomas algorithm
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        integer(int32) :: i
        do i = 1, NumberXNodes
            SELECT CASE (world%boundaryConditions(i))
            CASE(0)
                if (i < NumberXNodes) then
                    self % c_tri(i) = 1.0d0/(world%delX)
                end if
                if (i > 1) then
                    self%a_tri(i-1) = 1.0d0/(world%delX)
                end if
                self%b_tri(i) = - 2.0d0 / (world%delX)
            CASE(1,4)
                self%b_tri(i) = 1.0d0
            CASE(2)
                if (i == 1) then
                    self % c_tri(i) = 1.0d0/(world%delX)
                    self%b_tri(i) = - 1.0d0 / (world%delX)
                else if (i == NumberXNodes) then
                    self % a_tri(i-1) = 1.0d0/(world%delX)
                    self%b_tri(i) = - 1.0d0 / (world%delX)
                else
                    print *, "Neumann boundary reflecting not on left or right most index!"
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
        res = eps_0 * world%delX * (0.5d0 * (self%EField(1)**2 + self%EField(NumberXNodes)**2) + SUM(self%EField(2:NumberXNodes-1)**2))
        
    end function getTotalPE


    


    !----------------- Particle mover/sub-stepping procedures -------------------------------------------------------

    subroutine makeEField(self, world)
        ! Generate EField at nodes
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        self%EField(2:NumberXNodes-1) = (self%phi(1:NumberXNodes-2) - self%phi(3:NumberXNodes))/2.0d0/world%delX
        SELECT CASE (world%boundaryConditions(1))
        CASE(1,4)
            self%EField(1) = (self%phi(1) - self%phi(2))/world%delX !Use first order so in line with implicit, 0.5d0 * (3.0d0 * self%phi(1) - 4.0d0 * self%phi(2) + self%phi(3)) / world%delX
        CASE(2)
            self%EField(1) = 0.0d0
        CASE(3)
            self%EField(1) = (self%phi(NumberXNodes-1) - self%phi(2))/2.0d0/world%delX
        CASE default
            print *, "No case makeEField"
            stop
        END SELECT

        SELECT CASE (world%boundaryConditions(NumberXNodes))
        CASE(1,4)
            self%EField(NumberXNodes) = (self%phi(NumberXNodes-1) - self%phi(NumberXNodes))/world%delX !Use first order so in line with implicit, 0.5d0 * (-3.0d0 * self%phi(NumberXNodes) + 4.0d0 * self%phi(NumberXNodes-1) - self%phi(NumberXNodes-2))/ world%delX
        CASE(2)
            self%EField(NumberXNodes) = 0.0d0
        CASE(3)
            self%EField(NumberXNodes) = (self%phi(NumberXNodes-1) - self%phi(2))/2.0d0/world%delX
        CASE default
            print *, "No case makeEField"
            stop
        END SELECT
    end subroutine makeEField

    pure function getEField(self, l_p) result(res)
        !return EField of particle at computational position l_p
        class(potentialSolver), intent(in) :: self
        real(real64), intent(in) :: l_p
        integer(int32) :: l_left
        real(real64) :: res, d
        l_left = INT(l_p)
        d = l_p - l_left
        res = self%EField(l_left) * (1.0d0 - d) + self%EField(l_left + 1) * d
    end function getEField

    ! -------------------------------------------- Particle mover without boolean checks for depositing J -----------------------------------------------------------

    function getTotalE(self, particleList, world) result(res)
        ! Get total kinetic + potential energy using average velocity (v^k-1/2 + v^k+1/2) / 2
        class(potentialSolver), intent(in) :: self
        type(Domain), intent(in) :: world
        type(Particle), intent(in) :: particleList(:)
        real(real64) :: res, temp(numThread), PE, v, q_m_ratio
        integer(int32) :: i, j, iThread
        PE = self%getTotalPE(world)
        temp = 0.0d0
        res = 0.0d0
        do j = 1, numberChargedParticles
            !$OMP parallel private(iThread, i, v)
            iThread = omp_get_thread_num() + 1
            loopParticles: do i = 1, particleList(j)%N_p(iThread)
                v = particleList(j)%phaseSpace(2, i, iThread) + 0.5d0 * particleList(j)%q_over_m * self%getEField(particleList(j)%phaseSpace(1, i, iThread)) * del_t
                temp(iThread) = temp(iThread) + (v**2 + SUM(particleList(j)%phaseSpace(3:4, i, iThread)**2)) * particleList(j)%w_p * particleList(j)%mass * 0.5d0
            end do loopParticles
            !$OMP end parallel
        end do
        res = SUM(temp) + PE
    end function getTotalE


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
        read(10, *, IOSTAT = io) leftVoltage, rightVoltage
        read(10, *, IOSTAT = io) 
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