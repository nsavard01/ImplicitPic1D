module mod_NullCollision

    use iso_fortran_env, only: int32, real64, output_unit
    use constants
    use mod_BasicFunctions
    use mod_particle
    use mod_targetParticle
    use omp_lib
    implicit none

    private
    public :: nullCollision
    

    ! Particle contains particle properties and stored values in phase space df
    type :: nullCollision
        integer(int32) :: numberCollisions, numberReactants, lengthArrays
        real(real64), allocatable :: energyArray(:), sigmaVArray(:, :), energyThreshold(:)
        real(real64) :: sigmaVMax, reducedMass, minEnergy, maxEnergy
        integer(int32), allocatable :: collisionType(:), reactantsIndx(:), numberProducts(:), productsIndx(:,:)
    contains
        procedure, public, pass(self) :: generateCollision

    end type nullCollision


    interface nullCollision 
        procedure :: nullCollision_constructor
    end interface nullCollision

contains

    type(nullCollision) function nullCollision_constructor(numberReactants, numberCollisions, lengthArrays, red_mass, energyArray, sigmaVArray, energyThreshold, collisionType, reactantsIndx, numberProducts, productsIndx) result(self)
        ! Construct particle object, sizeIncrease is fraction larger stored array compared to initial amount of particles
        ! In future, use hash function for possible k = 1 .. Nx, m amount of boundaries, p = prime number  m < p < N_x. h(k) = (k%p)%m
        integer(int32), intent(in) :: numberCollisions, numberReactants, lengthArrays
        real(real64), intent(in) :: energyArray(lengthArrays), sigmaVArray(lengthArrays, numberCollisions), energyThreshold(numberCollisions), red_mass
        integer(int32), intent(in) :: collisionType(numberCollisions), reactantsIndx(numberReactants), numberProducts(numberCollisions), productsIndx(3,numberCollisions)
        integer(int32) :: i, j, indxOrder(numberCollisions), k
        real(real64) :: sigma_v_max(numberCollisions)
        self%numberCollisions = numberCollisions
        self%numberReactants = numberReactants
        self%lengthArrays = lengthArrays
        self%reducedMass = red_mass

        ! order collisions by maximum sigma_v in arrays
        do i = 1, numberCollisions
            sigma_v_max(i) = 0
            do j = 1, lengthArrays
                sigma_v_max(i) = MAX(sigma_v_max(i), sigmaVArray(j, i))
            end do
        end do
        call indexSortArray(self%numberCollisions, sigma_v_max, indxOrder)
        ! swap indx Order
        do i = 1, numberCollisions/2
            k = indxOrder(i)
            indxOrder(i) = indxOrder(numberCollisions-i + 1)
            indxOrder(numberCollisions-i + 1) = k
        end do
        allocate(self%energyArray(lengthArrays), self%sigmaVArray(lengthArrays, numberCollisions), self%energyThreshold(numberCollisions), self%collisionType(numberCollisions), &
            self%reactantsIndx(numberReactants), self%numberProducts(numberCollisions), self%productsIndx(3, numberCollisions))
        self%energyArray = energyArray
        self%reactantsIndx = reactantsIndx
        do i = 1, numberCollisions
            ! Go top down in order, assign arrays. This is just to put more 'likely' colliison first when looking at collision probabilities
            j = indxOrder(i)
            self%sigmaVArray(:, i) = sigmaVArray(:, j)
            self%energyThreshold(i) = energyThreshold(j)
            self%collisionType(i) = collisionType(j)
            self%numberProducts(i) = numberProducts(j)
            self%productsIndx(:,i) = productsIndx(:,j)
        end do
        
        self%sigmaVMax = MAXVAL(SUM(self%sigmaVArray, DIM=2))
        self%minEnergy = MINVAL(self%energyArray)
        self%maxEnergy = MAXVAL(self%energyArray)
    
    end function nullCollision_constructor

    subroutine generateCollision(self, particleList, targetParticleList, numberChargedParticles, numberBinaryCollisions, irand, del_t)
        class(nullCollision), intent(in) :: self
        integer(int32), intent(in) :: numberChargedParticles, numberBinaryCollisions
        type(Particle), intent(in out) :: particleList(numberChargedParticles)
        type(targetParticle), intent(in) :: targetParticleList(numberBinaryCollisions)
        integer(int32), intent(in out) :: irand(numThread)
        real(real64), intent(in) :: del_t
        logical :: collisionBool
        real(real64) :: P_null, numberSelectedReal, Rand, particleEnergy, targetVelocity(3), incidentVelocity(3), velocity_CM(3), energyCM, d_value, sigma_v, sigma_v_low
        integer(int32) :: numberSelected, iThread, i, particleIndx, numberTotalParticles, indxHigh, indxLow, indxMiddle, collIdx

        P_null = 1.0d0 - EXP(-self%sigmaVMax * targetParticleList(self%reactantsIndx(2))%density * del_t)
        if (P_null > 0.5d0) then
            print *, 'P_null greater than 50%'
            stop
        end if
        
        !$OMP parallel private(iThread, i,numberSelected, numberSelectedReal, Rand, particleIndx, numberTotalParticles, targetVelocity, incidentVelocity, &
            velocity_CM, energyCM, d_value, indxHigh, indxLow, indxMiddle, collIdx, sigma_v, sigma_v_low, collisionBool)
        iThread = omp_get_thread_num() + 1
        numberTotalParticles = particleList(self%reactantsIndx(1))%N_p(iThread)
        numberSelectedReal = P_null * real(numberTotalParticles)
        numberSelected = INT(numberSelectedReal)
        Rand = ran2(irand(iThread))
        if (Rand < (numberSelectedReal - numberSelected)) numberSelected = numberSelected + 1   
        do i = 1, numberSelected
            Rand = ran2(irand(iThread))
            particleIndx = INT((numberTotalParticles) * Rand) + 1
            incidentVelocity = particleList(self%reactantsIndx(1))%phaseSpace(2:4, particleIndx, iThread)
            targetVelocity = targetParticleList(self%reactantsIndx(2))%generate3DMaxwellianVelocity(irand(iThread))
            velocity_CM = incidentVelocity - targetVelocity
            energyCM = SUM(velocity_CM**2) * 0.5d0 * self%reducedMass / e !in eV
            indxLow = 1
            indxHigh = self%lengthArrays
            if (energyCM < self%minEnergy) then
                ! take minimum sigma_v
                d_value = 0.0d0
            else if (energyCM > self%maxEnergy) then
                ! take maximum sigma_v
                indxLow = indxHigh -1
                d_value = 1.0d0
            else
                ! binary search
                do while (indxLow /= indxHigh-1)
                    indxMiddle = (indxLow + indxHigh)/2
                    if (self%energyArray(indxMiddle) < energyCM) then
                        indxLow = indxMiddle
                    else if (self%energyArray(indxMiddle) > energyCM) then
                        indxHigh = indxMiddle
                    else
                        indxLow = indxMiddle
                        indxHigh = indxLow + 1
                    end if
                end do
                d_value = (energyCM - self%energyArray(indxLow))/(self%energyArray(indxHigh) - self%energyArray(indxLow))
            end if
            Rand = ran2(irand(iThread))
            sigma_v_low = 0.0d0
            do collIdx = 1, self%numberCollisions
                sigma_v = self%sigmaVArray(indxLow, collIdx) * (1.0d0 - d_value) + self%sigmaVArray(indxHigh, collIdx) * d_value + sigma_v_low
                collisionBool = (Rand <= sigma_v/self%sigmaVMax)
                if (collisionBool) then
                    SELECT CASE (self%collisionType(collIdx))
                    CASE(1)
                        print *, 'elastic'
                    CASE(2)
                        print *, 'ionization'
                    CASE(3)
                        print *, 'excitation'
                    CASE(4)
                        print *, 'charge exchange'
                    CASE(5)
                        print *, 'dissociation'
                    END SELECT
                    exit
                end if
                sigma_v_low = sigma_v
            end do
            if (.not. collisionBool) then
                print *, 'null collision'
            end if
        end do
        !$OMP end parallel


    end subroutine generateCollision

    ! subroutine ionizationCollisionIsotropic(electron, ion, n_g, sigma, del_t, E_thres, T, irand)
    !     ! Temporary subroutine for artificial ionization, using MCC, energy units in eV
    !     ! T in eV
    !     ! divide energy among electrons, isotropic scatter
    !     ! For moment assume null collision with gas temperature = 0
    !     type(Particle), intent(in out) :: electron, ion
    !     real(real64), intent(in) :: n_g, sigma, del_t, E_thres, T
    !     integer(int32), intent(in out) :: irand
    !     real(real64) :: electronSpeed, Pcoll, R, speedPerParticle, phi, theta, e_vector(3), V_cm(3), delE, V_neutral(3), mu, mass_neutral, sumMassInverse
    !     integer(int32) :: i, addIdx
    !     addIdx = 0
    !     mass_neutral = ion%mass + m_e
    !     sumMassInverse = (2.0d0/m_e + 1.0d0/ion%mass)
    !     mu = (m_e * mass_neutral)/ (m_e + mass_neutral)
    !     do i=1, electron%N_p
    !         call get3DMaxwellianVelocity(V_neutral, mass_neutral, T, irand)
    !         delE = 0.5d0 * mu * SUM((electron%phaseSpace(2:4, i) - V_neutral)**2) - E_thres*e
    !         if (delE > 0) then
    !             R = ran2(irand)
    !             electronSpeed = SQRT(SUM((electron%phaseSpace(2:4, i) - V_neutral)**2))
    !             Pcoll = 1.0d0 - EXP(-n_g * sigma * del_t * electronSpeed)
    !             if (Pcoll > R) then
    !                 addIdx = addIdx + 1
    !                 V_cm = (m_e *electron%phaseSpace(2:4, i) + V_neutral * mass_neutral) / (mass_neutral + m_e)

    !                 ! first add to electron
    !                 theta = ACOS(1.0d0 - 2.0d0*ran2(irand))
    !                 phi = ran2(irand) * 2 * pi
    !                 e_vector(1) = COS(theta)
    !                 e_vector(2) = SIN(phi) * SIN(theta)
    !                 e_vector(3) = COS(phi) * SIN(theta)
    !                 speedPerParticle = SQRT(2.0d0 * delE / sumMassInverse / m_e**2)
    !                 electron%phaseSpace(2:4, i) = e_vector * speedPerParticle + V_cm

    !                 ! second electron
    !                 e_vector(1) = COS(theta + 2.0d0 * pi / 3.0d0)
    !                 e_vector(2) = SIN(phi) * SIN(theta + 2.0d0 * pi / 3.0d0)
    !                 e_vector(3) = COS(phi) * SIN(theta + 2.0d0 * pi / 3.0d0)
    !                 electron%phaseSpace(2:4, electron%N_p + addIdx) = e_vector * speedPerParticle + V_cm
    !                 electron%phaseSpace(1,electron%N_p + addIdx) = electron%phaseSpace(1,i)

    !                 ! ion
    !                 e_vector(1) = COS(theta + 4.0d0 * pi / 3.0d0)
    !                 e_vector(2) = SIN(phi) * SIN(theta + 4.0d0 * pi / 3.0d0)
    !                 e_vector(3) = COS(phi) * SIN(theta + 4.0d0 * pi / 3.0d0)
    !                 speedPerParticle = SQRT(2.0d0 * delE / sumMassInverse / ion%mass**2)
    !                 ion%phaseSpace(2:4, ion%N_p + addIdx) = e_vector * speedPerParticle + V_cm
    !                 ion%phaseSpace(1,ion%N_p + addIdx) = electron%phaseSpace(1,i)
    !             end if
    !         end if
    !     end do
    !     electron%N_p = electron%N_p + addIdx
    !     ion%N_p = ion%N_p + addIdx
    !     inelasticEnergyLoss = inelasticEnergyLoss + e*E_thres * addIdx * electron%w_p
    ! end subroutine ionizationCollisionIsotropic


end module mod_NullCollision