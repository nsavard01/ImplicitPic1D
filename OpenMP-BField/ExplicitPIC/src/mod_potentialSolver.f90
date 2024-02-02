module mod_potentialSolver
    use iso_fortran_env, only: int32, real64, output_unit
    use constants
    use mod_BasicFunctions
    use mod_particle
    use mod_domain
    use omp_lib
    implicit none

    ! This module will be for the potential solver
    ! This uses the interaction between particles (particle module) and grid (domain module) in order to solve for the potential
    ! Will contain particle to mesh gathers (J, rho) and potential which comes from them (phi)
    ! Will also contain particle mover, since needed to store to J, and cannot be separate
    ! Assume dirichlet boundaries at ends for now, so matrix solver can go down by two dimensions
    public :: potentialSolver

    type :: potentialSolver
        real(real64), allocatable :: phi(:), rho(:, :), EField(:) !phi_f is final phi, will likely need to store two arrays for phi, can't be avoided
        real(real64) :: energyError, rho_const, siedelIter, siedelEps, BFieldMag, BField(3), BFieldAngle, RF_rad_frequency, RF_half_amplitude
        real(real64), allocatable :: a_tri(:), b_tri(:), c_tri(:) !for thomas algorithm potential solver, a_tri is lower diagonal, b_tri middle, c_tri upper
        logical :: BFieldBool


    contains
        procedure, public, pass(self) :: depositRho
        procedure, public, pass(self) :: solve_tridiag_Poisson
        procedure, public, pass(self) :: solve_gaussSiedel_Poisson
        procedure, public, pass(self) :: solve_CG_Poisson
        procedure, public, pass(self) :: getEField
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

    type(potentialSolver) function potentialSolver_constructor(world, leftVoltage, rightVoltage, BFieldMag, angle, RF_frequency) result(self)
        ! Construct domain object, initialize grid, dx_dl, and dx_dl.
        real(real64), intent(in) :: leftVoltage, rightVoltage
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: BFieldMag, angle, RF_frequency
        real(real64) :: angle_rad
        allocate(self % rho(NumberXNodes, numThread), self % phi(NumberXNodes), self%EField(NumberXNodes), self%a_tri(NumberXNodes-1), &
        self%b_tri(NumberXNodes), self%c_tri(NumberXNodes-1))
        self % a_tri = 0.0d0
        self % c_tri = 0.0d0
        self % b_tri = 0.0d0
        self % rho = 0.0d0
        self % rho_const = 0.0d0
        self % phi = 0.0d0
        self % EField = 0.0d0
        self%energyError = 0.0d0
        self%siedelIter = 100000
        self%BFieldMag = BFieldMag
        self%BFieldBool = (BFieldMag /= 0.0d0)
        angle_rad = angle * pi / 180.0d0
        self%BFieldAngle = angle_rad
        self%BField(1) = BFieldMag * COS(angle_rad)
        self%BField(2) = BFieldMag * SIN(angle_rad)
        self%BField(3) = 0.0d0
        self%siedelEps = 1d-6
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

    subroutine depositRho(self, particleList) 
        class(potentialSolver), intent(in out) :: self
        type(Particle), intent(in) :: particleList(:)
        integer(int32) :: i, j, l_left, l_right, iThread
        real(real64) :: d
        self % rho = self%rho_const
        do i=1, numberChargedParticles
            !$OMP parallel private(iThread, j, l_left, l_right, d)
            iThread = omp_get_thread_num() + 1
            do j = 1, particleList(i)%N_p(iThread)
                l_left = INT(particleList(i)%phaseSpace(1, j, iThread))
                l_right = l_left + 1
                d = particleList(i)%phaseSpace(1, j, iThread) - real(l_left)
                self % rho(l_left, iThread) = self % rho(l_left, iThread) + (1.0d0-d) * particleList(i)%w_p * particleList(i)%q
                self % rho(l_right, iThread) = self % rho(l_right, iThread) + d * particleList(i)%w_p * particleList(i)%q
            end do
            !$OMP end parallel
        end do
    end subroutine depositRho

    subroutine construct_diagMatrix(self, world)
        ! construct diagonal components for thomas algorithm
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
                    !self%a_tri(i - leftNodeIdx) = 2.0d0/(world%dx_dl(i-1) + world%dx_dl(i)) / world%dx_dl(i-1)
                    self%b_tri(i) = - 1.0d0 / (world%delX)
                else if (i == NumberXNodes) then
                    self % a_tri(i-1) = 1.0d0/(world%delX)
                    !self % c_tri(i - leftNodeIdx) = 2.0d0/(world%dx_dl(i-2) + world%dx_dl(i-1))/world%dx_dl(i-1)
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
        ! Tridiagonal (Thomas algorithm) solver for initial Poisson
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: timeCurrent
        integer(int32) :: i
        real(real64) :: m, d(NumberXNodes), cp(NumberXNodes-1),dp(NumberXNodes)

        do i=1, NumberXNodes
            SELECT CASE (world%boundaryConditions(i))
            CASE(0,2)
                d(i) = (-SUM(self%rho(i, :))) / eps_0
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

    subroutine solve_gaussSiedel_Poisson(self, world)
        ! Tridiagonal (Thomas algorithm) solver for initial Poisson
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        integer(int32) :: i, j
        real(real64) :: b(NumberXNodes), sigma(NumberXNodes), Ax(NumberXNodes), res

        do i=1, NumberXNodes
            SELECT CASE (world%boundaryConditions(i))
            CASE(0,2)
                b(i) = (-SUM(self%rho(i, :))) / eps_0
            CASE(1,3)
                b(i) = self%phi(i)
            END SELECT
        end do
        do i = 1, self%siedelIter
            sigma(1) = self%c_tri(1)*self%phi(2)
            self%phi(1) = self%phi(1) + 1.4d0*((b(1) - sigma(1))/self%b_tri(1) - self%phi(1))
            do j = 2, NumberXNodes-1
                sigma(j) = self%c_tri(j)*self%phi(j+1) + self%a_tri(j-1) * self%phi(j-1)    
                self%phi(j) = self%phi(j) + 1.4d0*((b(j) - sigma(j))/self%b_tri(j) - self%phi(j))
            end do
            sigma(NumberXNodes) = self%a_tri(NumberXNodes-1) * self%phi(NumberXNodes-1)
            self%phi(NumberXNodes) = self%phi(NumberXNodes) + 1.4d0 * ((b(NumberXNodes) - sigma(NumberXNodes))/self%b_tri(NumberXNodes) - self%phi(NumberXNodes))
            if (MOD(i, 100) == 0) then
                Ax = triMul(NumberXNodes, self%a_tri, self%c_tri, self%b_tri, self%phi)
                res = SUM(((Ax + 1d-15)/(b + 1d-15) - 1.0d0)**2)
                res = SQRT(res/NumberXNodes)
                if (res < self%siedelEps) then
                    exit
                end if
            end if
        end do
        if (i == self%siedelIter+1) then
            stop 'Max iterations reached Gauss-siedel solver!'
        end if
    end subroutine solve_gaussSiedel_Poisson

    subroutine solve_CG_Poisson(self, world)
        ! Tridiagonal (Thomas algorithm) solver for initial Poisson
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        integer(int32) :: i
        real(real64) :: b(NumberXNodes), RPast(NumberXNodes), RFuture(NumberXNodes), D(NumberXNodes), beta, alpha, resPast, resFuture, Ax(NumberXNodes)
        logical :: converge
        converge = .false.
        do i=1, NumberXNodes
            SELECT CASE (world%boundaryConditions(i))
            CASE(0,2)
                b(i) = (-SUM(self%rho(i, :))) / eps_0
            CASE(1,3)
                b(i) = self%phi(i)
            END SELECT
        end do
        RPast = b - triMul(NumberXNodes, self%a_tri, self%c_tri, self%b_tri, self%phi)
        resPast = SQRT(SUM(RPast**2))
        D = RPast
        do i = 1, 1000
            alpha = resPast**2 / SUM(D * triMul(NumberXNodes, self%a_tri, self%c_tri, self%b_tri, D))
            self%phi = self%phi + alpha * D
            Ax = triMul(NumberXNodes, self%a_tri, self%c_tri, self%b_tri, self%phi)
            RFuture = b - Ax
            resFuture = SQRT(SUM(RFuture**2))
            if (SUM(ABS(RFuture/(b + 1.d-15)))/NumberXNodes < 1.d-8) then
                converge = .true.
                exit
            end if
            beta = (resFuture**2)/(resPast**2)
            D = RFuture + beta * D
            RPast = RFuture
            resPast = resFuture
        end do
        if (.not. converge) then
            stop 'Max iterations reached CG solver!'
        end if
    end subroutine solve_CG_Poisson


    function getTotalPE(self, world) result(res)
        ! Get energy in electric fields, future true, then derive from phi_f, otherwise phi
        ! In 1D is J/m^2
        class(potentialSolver), intent(in) :: self
        type(Domain), intent(in) :: world
        real(real64) :: res
        res = 0.5 * eps_0 * SUM(arrayDiff(self%phi, NumberXNodes)**2) / world%delX
        
    end function getTotalPE


    


    !----------------- Particle mover/sub-stepping procedures -------------------------------------------------------

    subroutine makeEField(self, world)
        !return EField of particle at logical position l_p (this is half distance), per particle since particle mover loop will be per particle
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        self%EField(2:NumberXNodes-1) = (self%phi(1:NumberXNodes-2) - self%phi(3:NumberXNodes))/2.0d0/world%delX
        SELECT CASE (world%boundaryConditions(1))
        CASE(1,4)
            self%EField(1) = 0.5d0 * (3.0d0 * self%phi(1) - 4.0d0 * self%phi(2) + self%phi(3)) / world%delX
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
            self%EField(NumberXNodes) = 0.5d0 * (-3.0d0 * self%phi(NumberXNodes) + 4.0d0 * self%phi(NumberXNodes-1) - self%phi(NumberXNodes-2))/ world%delX
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
        !return EField of particle at logical position l_p (this is half distance), per particle since particle mover loop will be per particle
        class(potentialSolver), intent(in) :: self
        real(real64), intent(in) :: l_p
        integer(int32) :: l_left
        real(real64) :: res, d
        l_left = INT(l_p)
        d = l_p - l_left
        res = self%EField(l_left) * (1.0d0 - d) + self%EField(l_left + 1) * d
    end function getEField

    ! -------------------------------------------- Particle mover without boolean checks for depositing J -----------------------------------------------------------

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

    function getTotalE(self, particleList, world) result(res)
        class(potentialSolver), intent(in) :: self
        type(Domain), intent(in) :: world
        type(Particle), intent(in) :: particleList(:)
        real(real64) :: res, temp(numThread), PE, v, q_m_ratio
        integer(int32) :: i, j, iThread
        PE = self%getTotalPE(world)
        temp = 0.0d0
        res = 0.0d0
        do j = 1, numberChargedParticles
            q_m_ratio = particleList(j)%q/particleList(j)%mass
            !$OMP parallel private(iThread, i, v)
            iThread = omp_get_thread_num() + 1
            loopParticles: do i = 1, particleList(j)%N_p(iThread)
                v = particleList(j)%phaseSpace(2, i, iThread) + 0.5d0 * (q_m_ratio) * self%getEField(particleList(j)%phaseSpace(1, i, iThread)) * del_t
                temp(iThread) = temp(iThread) + (v**2 + SUM(particleList(j)%phaseSpace(3:4, i, iThread)**2)) * particleList(j)%w_p * particleList(j)%mass * 0.5d0
            end do loopParticles
            !$OMP end parallel
        end do
        res = SUM(temp) + PE
    end function getTotalE

    subroutine moveParticles(self, particleList, world, del_t)
        ! particle mover to avoid the boolean checks which mostly don't happen when depositing J
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        type(Particle), intent(in out) :: particleList(:)
        real(real64), intent(in) :: del_t
        !a and c correspond to quadratic equations | l_alongV is nearest integer boundary along velocity component, away is opposite
        integer(int32) :: j, i, delIdx, iThread
        real(real64) :: v_minus(3), v_plus(3), coeffField, tBoris(3), v_prime(3), sBoris(3)
        loopSpecies: do j = 1, numberChargedParticles
            !$OMP parallel private(iThread, i, delIdx, v_minus, v_plus, coeffField, tBoris, v_prime, sBoris)
            iThread = omp_get_thread_num() + 1
            delIdx = 0
            particleList(j)%refIdx(iThread) = 0
            particleList(j)%energyLoss(:, iThread) = 0.0d0
            particleList(j)%wallLoss(:, iThread) = 0.0d0
            loopParticles: do i = 1, particleList(j)%N_p(iThread)
                ! First half step acceleration
                if (self%BFieldBool) then
                    coeffField = 0.5d0 * (particleList(j)%q/particleList(j)%mass) * self%getEField(particleList(j)%phaseSpace(1, i, iThread)) * del_t
                    ! First half step acceleration
                    v_minus(1) =  particleList(j)%phaseSpace(2, i, iThread) + coeffField
                    v_minus(2:3) = particleList(j)%phaseSpace(3:4, i, iThread)

                    !v_prime, first half rotation
                    tBoris = 0.5d0 * (particleList(j)%q/particleList(j)%mass) * self%BField * del_t
                    v_prime = v_minus + crossProduct(v_minus, tBoris)

                    !v_plus, second half rotation
                    sBoris = 2.0d0 * tBoris / (1.0d0 + SUM(tBoris**2))
                    v_plus = v_minus + crossProduct(v_prime, sBoris)

                    particleList(j)%phaseSpace(2, i-delIdx, iThread) = v_plus(1) + coeffField
                    particleList(j)%phaseSpace(3:4, i-delIdx, iThread) = v_plus(2:3)

                else    
                    particleList(j)%phaseSpace(2, i-delIdx, iThread) = particleList(j)%phaseSpace(2, i, iThread) + (particleList(j)%q/particleList(j)%mass) * self%getEField(particleList(j)%phaseSpace(1, i, iThread)) * del_t
                    particleList(j)%phaseSpace(3:4, i-delIdx, iThread) = particleList(j)%phaseSpace(3:4, i, iThread)
                end if
                ! Get new position
                particleList(j)%phaseSpace(1, i-delIdx, iThread) = particleList(j)%phaseSpace(1, i, iThread) + particleList(j)%phaseSpace(2, i-delIdx, iThread) * del_t/world%delX
                if (particleList(j)%phaseSpace(1, i-delIdx, iThread) <= 1) then
                    SELECT CASE (world%boundaryConditions(1))
                    CASE(1,4)
                        particleList(j)%energyLoss(1, iThread) = particleList(j)%energyLoss(1, iThread) + SUM(particleList(j)%phaseSpace(2:4, i-delIdx, iThread)**2)
                        particleList(j)%wallLoss(1, iThread) = particleList(j)%wallLoss(1, iThread) + 1 !C/m^2 in 1D
                        delIdx = delIdx + 1
                    CASE(2)
                        particleList(j)%phaseSpace(1, i-delIdx, iThread) = 2.0d0 - particleList(j)%phaseSpace(1, i-delIdx, iThread)
                        particleList(j)%phaseSpace(2:4, i-delIdx, iThread) = -particleList(j)%phaseSpace(2:4, i-delIdx, iThread)
                        particleList(j)%refIdx(iThread) = particleList(j)%refIdx(iThread) + 1
                        particleList(j)%refRecordIdx(particleList(j)%refIdx(iThread), iThread) = i - delIdx
                    CASE(3)
                        particleList(j)%phaseSpace(1, i-delIdx, iThread) = MODULO(particleList(j)%phaseSpace(1, i-delIdx, iThread) - 2.0d0, real(NumberXNodes, kind = real64)) + 1
                    CASE default
                        print *, 'no case, moveParticles'
                        stop
                    END SELECT
                else if ((particleList(j)%phaseSpace(1, i-delIdx, iThread) >= NumberXNodes)) then
                    SELECT CASE (world%boundaryConditions(NumberXNodes))
                    CASE(1,4)
                        particleList(j)%energyLoss(2, iThread) = particleList(j)%energyLoss(2, iThread) + SUM(particleList(j)%phaseSpace(2:4, i-delIdx, iThread)**2)
                        particleList(j)%wallLoss(2, iThread) = particleList(j)%wallLoss(2, iThread) + 1 !C/m^2 in 1D
                        delIdx = delIdx + 1
                    CASE(2)
                        particleList(j)%phaseSpace(1, i-delIdx, iThread) = 2.0d0 * NumberXNodes - particleList(j)%phaseSpace(1, i-delIdx, iThread)
                        particleList(j)%phaseSpace(2:4, i-delIdx, iThread) = -particleList(j)%phaseSpace(2:4, i-delIdx, iThread)
                        particleList(j)%refIdx(iThread) = particleList(j)%refIdx(iThread) + 1
                        particleList(j)%refRecordIdx(particleList(j)%refIdx(iThread), iThread) = i - delIdx
                    CASE(3)
                        particleList(j)%phaseSpace(1, i-delIdx, iThread) = MODULO(particleList(j)%phaseSpace(1, i-delIdx, iThread), real(NumberXNodes, kind = real64)) + 1
                    CASE default
                        print *, 'no case, moveParticles'
                        stop
                    END SELECT
                end if
            end do loopParticles
            particleList(j)%N_p(iThread) = particleList(j)%N_p(iThread) - delIdx
            particleList(j)%delIdx(iThread) = delIdx
            !$OMP end parallel
            particleList(j)%accumEnergyLoss = particleList(j)%accumEnergyLoss + SUM(particleList(j)%energyLoss, DIM = 2)
            particleList(j)%accumWallLoss = particleList(j)%accumWallLoss + SUM(particleList(j)%wallLoss, DIM = 2)
        end do loopSpecies
    end subroutine moveParticles

    ! ---------------- Initial Poisson Solver -------------------------------------------------

    subroutine solvePotential(self, particleList, world, timeCurrent)
        ! Solve for initial potential
        class(potentialSolver), intent(in out) :: self
        type(Particle), intent(in) :: particleList(:)
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: timeCurrent
        call self%depositRho(particleList)
        call self%solve_tridiag_Poisson(world, timeCurrent)
        ! Assume only use potential solver once, then need to generate matrix for Div-Ampere
        call self%makeEField(world)
    end subroutine solvePotential


end module mod_potentialSolver