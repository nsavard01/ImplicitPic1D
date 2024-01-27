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
        real(real64) :: P_null, numberSelectedReal, Rand, particleEnergy, targetVelocity(3), incidentVelocity(3), velocity_CM(3), energyCM, d_value, sigma_v, sigma_v_low, speedCM
        integer(int32) :: numberSelected, iThread, i, particleIndx, numberTotalParticles, indxHigh, indxLow, indxMiddle, collIdx, addAmount

        P_null = 1.0d0 - EXP(-self%sigmaVMax * targetParticleList(self%reactantsIndx(2))%density * del_t)
        if (P_null > 0.5d0) then
            print *, 'P_null greater than 50%'
            stop
        end if
        
        !$OMP parallel private(iThread, i,numberSelected, numberSelectedReal, Rand, particleIndx, numberTotalParticles, targetVelocity, incidentVelocity, &
            velocity_CM, energyCM, d_value, indxHigh, indxLow, indxMiddle, collIdx, sigma_v, sigma_v_low, collisionBool, speedCM, addAmount)
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
            speedCM = SQRT(SUM(velocity_CM**2))
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
                if (energyCM > self%energyThreshold(collIdx)) then
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
                end if
            end do
            if (.not. collisionBool) then
                print *, 'null collision'
            end if
        end do
        !$OMP end parallel
        print *, ''


    end subroutine generateCollision

    subroutine ionizationCollisionIsotropic(self, electron, ion, targetParticle, del_t, irand, energyCM, E_thres, speedCM, incidentVelocity, targetVelocity, velocityCM)
        ! Temporary subroutine for artificial ionization, using MCC, energy units in eV
        ! T in eV
        ! divide energy among electrons, isotropic scatter
        ! For moment assume null collision with gas temperature = 0
        class(nullCollision), intent(in) :: self
        type(Particle), intent(in out) :: electron, ion
        type(targetParticle), intent(in) :: targetParticle
        real(real64), intent(in) :: del_t, energyCM, E_thres, speedCM
        real(real64), intent(in out) :: incidentVelocity(3), targetVelocity(3), velocityCM(3)
        integer(int32), intent(in) :: partIndx, iThread
        integer(int32), intent(in out) :: irand, addAmount
        real(real64) :: speedPerParticle, phi, theta, e_vector(3), V_cm(3), delE, V_neutral(3), sumMassInverse, cos_theta, sin_theta, cos_phi, sin_phi, secTheta, thirdTheta
        integer(int32) :: i, addIdx
        
        sumMassInverse = (2.0d0/m_e + 1.0d0/ion%mass)
        delE = energyCM - E_thres
        
        V_cm = (m_e * incidentVelocity + targetVelocity * targetParticle%mass) / (targetParticle%mass + m_e)

        ! first add to electron
        theta = ACOS(1.0d0 - 2.0d0*ran2(irand))
        phi = ran2(irand) * 2 * pi
        cos_theta = COS(theta)
        sin_theta = SIN(theta)
        cos_phi = COS(phi)
        sin_phi = SIN(phi)
        e_vector(1) = cos_theta
        e_vector(2) = sin_phi * sin_theta
        e_vector(3) = cos_phi * sin_theta
        speedPerParticle = SQRT(2.0d0 * delE / sumMassInverse / m_e**2)
        incidentVelocity = e_vector * speedPerParticle + V_cm

        ! second electron
        secTheta = theta + 2.0d0 * pi/3.0d0
        e_vector(1) = COS(secTheta)
        e_vector(2) = sin_phi * SIN(secTheta)
        e_vector(3) = cos_phi * SIN(secTheta)
        velocityCM = e_vector * speedPerParticle + V_cm

        ! ion
        thirdTheta = theta + 4.0d0 * pi / 3.0d0
        e_vector(1) = COS(thirdTheta)
        e_vector(2) = sin_phi * SIN(thirdTheta)
        e_vector(3) = cos_phi * SIN(thirdTheta)
        speedPerParticle = SQRT(2.0d0 * delE / sumMassInverse / ion%mass**2)
        targetVelocity = e_vector * speedPerParticle + V_cm
           
    end subroutine ionizationCollisionIsotropic


end module mod_NullCollision