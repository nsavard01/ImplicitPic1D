module mod_potentialSolver
    use iso_fortran_env, only: int32, real64, output_unit
    use constants
    use mod_BasicFunctions
    use mod_particle
    use mod_domain
    implicit none

    ! This module will be for the potential solver
    ! This uses the interaction between particles (particle module) and grid (domain module) in order to solve for the potential
    ! Will contain particle to mesh gathers (J, rho) and potential which comes from them (phi)
    ! Will also contain particle mover, since needed to store to J, and cannot be separate
    ! Assume dirichlet boundaries at ends for now, so matrix solver can go down by two dimensions
    integer(int32) :: idxReFlux(1000, 2), reFluxMaxIdx(2), delIdx(2)
    public :: potentialSolver

    type :: potentialSolver
        real(real64), allocatable :: phi(:), rho(:), EField(:), particleChargeLoss(:,:) !phi_f is final phi, will likely need to store two arrays for phi, can't be avoided
        real(real64) :: energyError, particleEnergyLoss, rho_const
        real(real64), allocatable :: a_tri(:), b_tri(:), c_tri(:) !for thomas algorithm potential solver, a_tri is lower diagonal, b_tri middle, c_tri upper


    contains
        procedure, public, pass(self) :: depositRho
        procedure, public, pass(self) :: solve_tridiag_Poisson
        procedure, public, pass(self) :: getEField
        procedure, public, pass(self) :: makeEField
        procedure, public, pass(self) :: solvePotential
        procedure, public, pass(self) :: construct_diagMatrix
        procedure, public, pass(self) :: initialVRewind
        procedure, public, pass(self) :: getTotalPE
        procedure, public, pass(self) :: moveParticles
    end type

    interface potentialSolver
        module procedure :: potentialSolver_constructor
    end interface potentialSolver
   
contains

    type(potentialSolver) function potentialSolver_constructor(world, leftVoltage, rightVoltage) result(self)
        ! Construct domain object, initialize grid, dx_dl, and dx_dl.
        real(real64), intent(in) :: leftVoltage, rightVoltage
        type(Domain), intent(in) :: world
        allocate(self % rho(NumberXNodes), self % phi(NumberXNodes), self%EField(NumberXNodes), self%a_tri(NumberXNodes-1), &
        self%b_tri(NumberXNodes), self%c_tri(NumberXNodes-1), self%particleChargeLoss(2, numberChargedParticles))
        self % a_tri = 0.0d0
        self % c_tri = 0.0d0
        self % b_tri = 0.0d0
        self % rho = 0.0d0
        self % rho_const = 0.0d0
        self % phi = 0.0d0
        self % EField = 0.0d0
        self%energyError = 0.0d0
        self%particleEnergyLoss = 0.0d0
        if (world%boundaryConditions(1) == 1) self%phi(1) = leftVoltage
        if (world%boundaryConditions(NumberXNodes) == 1) self%phi(NumberXNodes) = rightVoltage
        if (world%boundaryConditions(1) == 3) then
            self%phi(1) = leftVoltage
            self%phi(NumberXNodes) = leftVoltage
        end if

    end function potentialSolver_constructor

    subroutine depositRho(self, particleList, world) 
        class(potentialSolver), intent(in out) :: self
        type(Particle), intent(in) :: particleList(:)
        type(Domain), intent(in) :: world
        integer(int32) :: i, j, l_left, l_right
        real(real64) :: d
        self % rho = 0.0d0
        do i=1, numberChargedParticles
            do j = 1, particleList(i)%N_p
                l_left = INT(particleList(i)%phaseSpace(1, j))
                l_right = l_left + 1
                d = particleList(i)%phaseSpace(1, j) - l_left
                SELECT CASE (world%boundaryConditions(l_left))
                CASE(0)
                    self % rho(l_left) = self % rho(l_left) + particleList(i)%q * particleList(i)%w_p * (1.0d0-d)
                CASE(1)
                    self % rho(l_left) = self % rho(l_left) + particleList(i)%q * particleList(i)%w_p * (1.0d0-d)
                CASE(2)
                    self % rho(l_left) = self % rho(l_left) + 2.0d0 * particleList(i)%q * particleList(i)%w_p * (1.0d0-d)
                CASE(3)
                    self % rho(l_left) = self % rho(l_left) + particleList(i)%q * particleList(i)%w_p * (1.0d0-d)
                    self % rho(ABS(l_left - NumberXNodes)+1) = self % rho(ABS(l_left - NumberXNodes)+1) + particleList(i)%q * particleList(i)%w_p * (1.0d0-d)
                CASE(4)
                    self % rho(l_left) = self % rho(l_left) + 2.0d0 * particleList(i)%q * particleList(i)%w_p * (1.0d0-d)
                CASE default
                    print *, "case doesn't exist in deposit rho"
                    print *, 'particle position:', particleList(i)%phaseSpace(1, j)
                    stop
                END SELECT

                SELECT CASE (world%boundaryConditions(l_right))
                CASE(0)
                    self % rho(l_right) = self % rho(l_right) + particleList(i)%q * particleList(i)%w_p * d
                CASE(1)
                    self % rho(l_right) = self % rho(l_right) + particleList(i)%q * particleList(i)%w_p * d
                CASE(2)
                    self % rho(l_right) = self % rho(l_right) + 2.0d0 * particleList(i)%q * particleList(i)%w_p * d
                CASE(3)
                    self % rho(l_right) = self % rho(l_right) + particleList(i)%q * particleList(i)%w_p * d
                    self % rho(ABS(l_right - NumberXNodes)+1) = self % rho(ABS(l_right-NumberXNodes)+1) + particleList(i)%q * particleList(i)%w_p * d
                CASE(4)
                    self % rho(l_right) = self % rho(l_right) + 2.0d0 * particleList(i)%q * particleList(i)%w_p * d
                CASE default
                    print *, "case doesn't exist in deposit rho"
                    print *, 'particle position:', particleList(i)%phaseSpace(1, j)
                    stop
                END SELECT
            end do
        end do
        self % rho = self % rho / world%delX
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
                    self % c_tri(i) = 1.0d0/(world%delX**2)
                end if
                if (i > 1) then
                    self%a_tri(i-1) = 1.0d0/(world%delX**2)
                end if
                self%b_tri(i) = - 2.0d0 / (world%delX**2)
            CASE(1)
                self%b_tri(i) = 1.0d0
            CASE(2)
                if (i == 1) then
                    self % c_tri(i) = 2.0d0/(world%delX**2)
                    !self%a_tri(i - leftNodeIdx) = 2.0d0/(world%dx_dl(i-1) + world%dx_dl(i)) / world%dx_dl(i-1)
                    self%b_tri(i) = - 2.0d0 / (world%delX**2)
                else if (i == NumberXNodes) then
                    self % a_tri(i-1) = 2.0d0/(world%delX**2)
                    !self % c_tri(i - leftNodeIdx) = 2.0d0/(world%dx_dl(i-2) + world%dx_dl(i-1))/world%dx_dl(i-1)
                    self%b_tri(i) = - 2.0d0 / (world%delX**2)
                else
                    print *, "Neumann boundary reflecting not on left or right most index!"
                    stop
                end if
            CASE(3)
                self%b_tri(i) = 1.0d0
            CASE(4)
                if (i == 1) then
                    self % c_tri(i) = -1.0d0
                    !self%a_tri(i - leftNodeIdx) = 2.0d0/(world%dx_dl(i-1) + world%dx_dl(i)) / world%dx_dl(i-1)
                    self%b_tri(i) = 1.0d0
                else if (i == NumberXNodes) then
                    self % a_tri(i-1) = -1.0d0
                    !self % c_tri(i - leftNodeIdx) = 2.0d0/(world%dx_dl(i-2) + world%dx_dl(i-1))/world%dx_dl(i-1)
                    self%b_tri(i) = 1.0d0
                else
                    print *, "Neumann boundary absorbing not on left or right most index!"
                    stop
                end if
            CASE default
                print *, "Error when constructing poisson matrix, inner nodes not plasma or neumann!"
            END SELECT
        end do

    end subroutine construct_diagMatrix

    subroutine solve_tridiag_Poisson(self, world)
        ! Tridiagonal (Thomas algorithm) solver for initial Poisson
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        integer(int32) :: i
        real(real64) :: m, d(NumberXNodes), cp(NumberXNodes-1),dp(NumberXNodes)

        do i=1, NumberXNodes
            SELECT CASE (world%boundaryConditions(i))
            CASE(0,2)
                d(i) = (-self%rho(i) - self%rho_const) / eps_0
            CASE(1,3)
                d(i) = self%phi(i)
            CASE(4)
                d(i) = 0.0d0
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
        CASE(1)
            self%EField(1) = 0.5d0 * (3.0d0 * self%phi(1) - 4.0d0 * self%phi(2) + self%phi(3)) / world%delX
        CASE(2)
            self%EField(1) = 0.0d0
        CASE(3)
            self%EField(1) = (self%phi(NumberXNodes-1) - self%phi(2))/2.0d0/world%delX
        CASE(4)
            self%EField(1) = 0.0d0
        CASE default
            print *, "No case makeEField"
            stop
        END SELECT

        SELECT CASE (world%boundaryConditions(NumberXNodes))
        CASE(1)
            self%EField(NumberXNodes) = 0.5d0 * (-3.0d0 * self%phi(NumberXNodes) + 4.0d0 * self%phi(NumberXNodes-1) - self%phi(NumberXNodes-2))/ world%delX
        CASE(2)
            self%EField(NumberXNodes) = 0.0d0
        CASE(3)
            self%EField(NumberXNodes) = (self%phi(NumberXNodes-1) - self%phi(2))/2.0d0/world%delX
        CASE(4)
            self%EField(NumberXNodes) = 0.0d0
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
        integer(int32) :: j, i
        loopSpecies: do j = 1, numberChargedParticles
            loopParticles: do i = 1, particleList(j)%N_p
                particleList(j)%phaseSpace(2, i) = particleList(j)%phaseSpace(2, i) - 0.5 * (particleList(j)%q/particleList(j)%mass) * self%getEField(particleList(j)%phaseSpace(1, i)) * del_t
            end do loopParticles
        end do loopSpecies
    end subroutine initialVRewind

    subroutine moveParticles(self, particleList, world, del_t)
        ! particle mover to avoid the boolean checks which mostly don't happen when depositing J
        class(potentialSolver), intent(in out) :: self
        type(Domain), intent(in) :: world
        type(Particle), intent(in out) :: particleList(:)
        real(real64), intent(in) :: del_t
        !a and c correspond to quadratic equations | l_alongV is nearest integer boundary along velocity component, away is opposite
        integer(int32) :: j, i
        delIdx = 0
        reFluxMaxIdx = 0
        loopSpecies: do j = 1, numberChargedParticles
            loopParticles: do i = 1, particleList(j)%N_p
                particleList(j)%phaseSpace(2, i-delIdx(j)) = particleList(j)%phaseSpace(2, i) + (particleList(j)%q/particleList(j)%mass) * self%getEField(particleList(j)%phaseSpace(1, i)) * del_t
                particleList(j)%phaseSpace(1, i-delIdx(j)) = particleList(j)%phaseSpace(1, i) + particleList(j)%phaseSpace(2, i-delIdx(j)) * del_t/world%delX
                particleList(j)%phaseSpace(3:4, i-delIdx(j)) = particleList(j)%phaseSpace(3:4, i)
                if (particleList(j)%phaseSpace(1, i-delIdx(j)) <= 1) then
                    SELECT CASE (world%boundaryConditions(1))
                    CASE(1)
                        self%particleChargeLoss(1, j) = self%particleChargeLoss(1, j) + particleList(j)%q * particleList(j)%w_p !C/m^2 in 1D
                        self%particleEnergyLoss = self%particleEnergyLoss + particleList(j)%w_p * SUM(particleList(j)%phaseSpace(2:4, i-delIdx(j))**2) * particleList(j)%mass * 0.5d0
                        delIdx(j) = delIdx(j) + 1
                    CASE(2)
                        reFluxMaxIdx(j) = reFluxMaxIdx(j) + 1
                        particleList(j)%phaseSpace(1, i-delIdx(j)) = 2.0d0 - particleList(j)%phaseSpace(1, i-delIdx(j))
                        particleList(j)%phaseSpace(2, i-delIdx(j)) = -particleList(j)%phaseSpace(2, i-delIdx(j))
                        idxReFlux(reFluxMaxIdx(j), j) = i - delIdx(j)
                    CASE(3)
                        particleList(j)%phaseSpace(1, i-delIdx(j)) = MODULO(particleList(j)%phaseSpace(1, i-delIdx(j)) - 2.0d0, real(NumberXNodes, kind = real64)) + 1
                    CASE default
                        print *, 'no case, moveParticles'
                        stop
                    END SELECT
                else if ((particleList(j)%phaseSpace(1, i-delIdx(j)) >= NumberXNodes)) then
                    SELECT CASE (world%boundaryConditions(NumberXNodes))
                    CASE(1)
                        self%particleChargeLoss(2, j) = self%particleChargeLoss(2, j) + particleList(j)%q * particleList(j)%w_p !C/m^2 in 1D
                        self%particleEnergyLoss = self%particleEnergyLoss + particleList(j)%w_p * SUM(particleList(j)%phaseSpace(2:4, i-delIdx(j))**2) * particleList(j)%mass * 0.5d0
                        delIdx(j) = delIdx(j) + 1
                    CASE(2)
                        reFluxMaxIdx(j) = reFluxMaxIdx(j) + 1
                        particleList(j)%phaseSpace(1, i-delIdx(j)) = 2.0d0 * NumberXNodes - particleList(j)%phaseSpace(1, i-delIdx(j))
                        particleList(j)%phaseSpace(2, i-delIdx(j)) = -particleList(j)%phaseSpace(2, i-delIdx(j))
                        idxReFlux(reFluxMaxIdx(j), j) = i - delIdx(j)
                    CASE(3)
                        particleList(j)%phaseSpace(1, i-delIdx(j)) = MODULO(particleList(j)%phaseSpace(1, i-delIdx(j)), real(NumberXNodes, kind = real64)) + 1
                    CASE default
                        print *, 'no case, moveParticles'
                        stop
                    END SELECT
                end if
            end do loopParticles
            particleList(j)%N_p = particleList(j)%N_p - delIdx(j)
        end do loopSpecies
    end subroutine moveParticles

    ! ---------------- Initial Poisson Solver -------------------------------------------------

    subroutine solvePotential(self, particleList, world)
        ! Solve for initial potential
        class(potentialSolver), intent(in out) :: self
        type(Particle), intent(in) :: particleList(:)
        type(Domain), intent(in) :: world
        call self%depositRho(particleList, world)
        call self%solve_tridiag_Poisson(world)
        ! Assume only use potential solver once, then need to generate matrix for Div-Ampere
        call self%makeEField(world)

    end subroutine solvePotential


end module mod_potentialSolver