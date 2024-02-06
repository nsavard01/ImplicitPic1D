module mod_NullCollision

    use iso_fortran_env, only: int32, real64, output_unit
    use constants
    use mod_BasicFunctions
    use mod_particle
    use mod_targetParticle
    use omp_lib
    use mod_mt19937
    implicit none

    private
    public :: nullCollision
    

    ! Particle contains particle properties and stored values in phase space df
    type :: nullCollision
        integer(int32) :: numberCollisions, numberReactants, lengthArrays, totalAmountCollisions
        real(real64), allocatable :: energyArray(:), sigmaArray(:, :), energyThreshold(:)
        real(real64) :: sigmaVMax, reducedMass, minEnergy, maxEnergy, sumMass, reducedMassIonization, totalEnergyLoss
        integer(int32), allocatable :: collisionType(:), reactantsIndx(:), numberProducts(:), productsIndx(:,:)
    contains
        procedure, public, pass(self) :: generateCollision
        procedure, public, pass(self) :: elasticExcitCollisionIsotropic
        procedure, public, pass(self) :: ionizationCollisionIsotropic
        procedure, public, pass(self) :: ionizationCollisionIsotropicNanbul
        procedure, public, pass(self) :: writeCollisionCrossSection

    end type nullCollision


    interface nullCollision 
        procedure :: nullCollision_constructor
    end interface nullCollision

contains

    type(nullCollision) function nullCollision_constructor(numberReactants, numberCollisions, lengthArrays, red_mass, sumMass, redMassTripleProducts, energyArray, sigmaArray, energyThreshold, collisionType, reactantsIndx, numberProducts, productsIndx) result(self)
        ! Construct particle object, sizeIncrease is fraction larger stored array compared to initial amount of particles
        ! In future, use hash function for possible k = 1 .. Nx, m amount of boundaries, p = prime number  m < p < N_x. h(k) = (k%p)%m
        integer(int32), intent(in) :: numberCollisions, numberReactants, lengthArrays
        real(real64), intent(in) :: energyArray(lengthArrays), sigmaArray(lengthArrays, numberCollisions), energyThreshold(numberCollisions), red_mass, sumMass, redMassTripleProducts
        integer(int32), intent(in) :: collisionType(numberCollisions), reactantsIndx(numberReactants), numberProducts(numberCollisions), productsIndx(3,numberCollisions)
        integer(int32) :: i, j, indxOrder(numberCollisions), k
        real(real64) :: sigma_v_max(numberCollisions), v_r, sum_sigma
        self%numberCollisions = numberCollisions
        self%numberReactants = numberReactants
        self%lengthArrays = lengthArrays
        self%sumMass = sumMass
        self%reducedMass = red_mass
        self%reducedMassIonization = redMassTripleProducts
        self%totalEnergyLoss = 0.0d0
        self%totalAmountCollisions = 0

        ! order collisions by maximum sigma_v in each collision
        do i = 1, numberCollisions
            sigma_v_max(i) = 0
            do j = 1, lengthArrays
                v_r = SQRT(2.0d0 * energyArray(j) * e / self%reducedMass)
                sigma_v_max(i) = MAX(sigma_v_max(i), sigmaArray(j, i) * v_r)
            end do
        end do
        call indexSortArray(self%numberCollisions, sigma_v_max, indxOrder)
        ! swap indx Order
        do i = 1, numberCollisions/2
            k = indxOrder(i)
            indxOrder(i) = indxOrder(numberCollisions-i + 1)
            indxOrder(numberCollisions-i + 1) = k
        end do

        
        allocate(self%energyArray(lengthArrays), self%sigmaArray(lengthArrays, numberCollisions), self%energyThreshold(numberCollisions), self%collisionType(numberCollisions), &
            self%reactantsIndx(numberReactants), self%numberProducts(numberCollisions), self%productsIndx(3, numberCollisions))
        self%energyArray = energyArray
        self%reactantsIndx = reactantsIndx
        do i = 1, numberCollisions
            ! Go top down in order, assign arrays. This is just to put more 'likely' colliison first when looking at collision probabilities
            j = indxOrder(i)
            self%sigmaArray(:, i) = sigmaArray(:, j)
            self%energyThreshold(i) = energyThreshold(j)
            self%collisionType(i) = collisionType(j)
            self%numberProducts(i) = numberProducts(j)
            self%productsIndx(:,i) = productsIndx(:,j)
        end do

        ! get sigma_v_max total along energy array
        self%sigmaVMax = 0.0d0
        do i = 1, lengthArrays
            sum_sigma = 0
            v_r = SQRT(2.0d0 * energyArray(i) * e / self%reducedMass)
            do j = 1, numberCollisions
                sum_sigma = sum_sigma + sigmaArray(i, j)
            end do
            self%sigmaVMax = MAX(self%sigmaVMax, sum_sigma * v_r)
        end do
        
        self%minEnergy = MINVAL(self%energyArray)
        self%maxEnergy = MAXVAL(self%energyArray)
    
    end function nullCollision_constructor

    subroutine generateCollision(self, particleList, targetParticleList, numberChargedParticles, numberBinaryCollisions, randGen, del_t)
        class(nullCollision), intent(in out) :: self
        integer(int32), intent(in) :: numberChargedParticles, numberBinaryCollisions
        type(Particle), intent(in out) :: particleList(numberChargedParticles)
        type(targetParticle), intent(in) :: targetParticleList(numberBinaryCollisions)
        type(mt19937), intent(in out) :: randGen(numThread)
        real(real64), intent(in) :: del_t
        logical :: collisionBool
        real(real64) :: P_null, numberSelectedReal, Rand, targetVelocity(3), incidentVelocity(3), velocity_CM(3), energyCM, d_value, sigma_v, sigma_v_low, speedCM, particleLocation, energyLoss
        integer(int32) :: numberSelected, iThread, i, particleIndx, numberTotalParticles, indxHigh, indxLow, indxMiddle, collIdx, addIonizationIndx, totalCollisions

        P_null = self%sigmaVMax * targetParticleList(self%reactantsIndx(2))%density * del_t !1.0d0 - EXP(-self%sigmaVMax * targetParticleList(self%reactantsIndx(2))%density * del_t)
    
        if (P_null > 0.5d0) then
            print *, 'P_null greater than 50%'
            stop
        end if

        energyLoss = 0
        totalCollisions = 0
        !$OMP parallel private(iThread, i,numberSelected, numberSelectedReal, Rand, particleIndx, numberTotalParticles, targetVelocity, incidentVelocity, &
            velocity_CM, energyCM, d_value, indxHigh, indxLow, indxMiddle, collIdx, sigma_v, sigma_v_low, collisionBool, speedCM, addIonizationIndx, particleLocation) reduction(+:energyLoss,totalCollisions)
        iThread = omp_get_thread_num() + 1
        numberTotalParticles = particleList(self%reactantsIndx(1))%N_p(iThread)
        numberSelectedReal = P_null * real(numberTotalParticles)
        numberSelected = INT(numberSelectedReal)
        Rand = randGen(iThread)%genrand64_real1()
        if (Rand < (numberSelectedReal - numberSelected)) numberSelected = numberSelected + 1   
        addIonizationIndx = 0
        do i = 1, numberSelected
            Rand = randGen(iThread)%genrand64_real1()
            particleIndx = INT((numberTotalParticles) * Rand) + 1
            incidentVelocity = particleList(self%reactantsIndx(1))%phaseSpace(2:4, particleIndx, iThread)
            targetVelocity = targetParticleList(self%reactantsIndx(2))%generate3DMaxwellianVelocity(randGen(iThread))
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
            Rand = randGen(iThread)%genrand64_real1()
            sigma_v_low = 0.0d0
            do collIdx = 1, self%numberCollisions
                if (energyCM > self%energyThreshold(collIdx)) then
                    sigma_v = (self%sigmaArray(indxLow, collIdx) * (1.0d0 - d_value) + self%sigmaArray(indxHigh, collIdx) * d_value) * speedCM + sigma_v_low
                    collisionBool = (Rand <= sigma_v/self%sigmaVMax)
                    if (collisionBool) then
                        totalCollisions = totalCollisions + 1
                        SELECT CASE (self%collisionType(collIdx))
                        CASE(1)
                            velocity_CM = incidentVelocity
                            call self%elasticExcitCollisionIsotropic(particleList(self%reactantsIndx(1)), targetParticleList(self%reactantsIndx(2)), &
                                randGen(iThread), energyCM, self%energyThreshold(collIdx), incidentVelocity, targetVelocity)
                            particleList(self%reactantsIndx(1))%phaseSpace(2:4, particleIndx, iThread) = incidentVelocity
                            energyLoss = energyLoss + particleList(self%reactantsIndx(1))%mass * 0.5d0 * (SUM(velocity_CM**2) - SUM(incidentVelocity**2)) / e
                        CASE(2)
                            call self%ionizationCollisionIsotropic(particleList(self%reactantsIndx(1)), particleList(self%productsIndx(2, collIdx)), targetParticleList(self%reactantsIndx(2)), &
                                randGen(iThread), energyCM, self%energyThreshold(collIdx), incidentVelocity, targetVelocity, velocity_CM)
                            particleLocation = particleList(self%reactantsIndx(1))%phaseSpace(1, particleIndx, iThread)
                            particleList(self%reactantsIndx(1))%N_p(iThread) = particleList(self%reactantsIndx(1))%N_p(iThread) + 1
                            particleList(self%productsIndx(2, collIdx))%N_p(iThread) = particleList(self%productsIndx(2, collIdx))%N_p(iThread) + 1
                            ! set primary particle velocity
                            particleList(self%reactantsIndx(1))%phaseSpace(2:4, particleIndx, iThread) = incidentVelocity
                            ! set secondary particle velocity and position
                            particleList(self%reactantsIndx(1))%phaseSpace(2:4, particleList(self%reactantsIndx(1))%N_p(iThread), iThread) = velocity_CM
                            particleList(self%reactantsIndx(1))%phaseSpace(1, particleList(self%reactantsIndx(1))%N_p(iThread), iThread) = particleLocation
                            ! set ion velocity and position
                            particleList(self%productsIndx(2, collIdx))%phaseSpace(2:4, particleList(self%productsIndx(2, collIdx))%N_p(iThread), iThread) = targetVelocity
                            particleList(self%productsIndx(2, collIdx))%phaseSpace(1, particleList(self%productsIndx(2, collIdx))%N_p(iThread), iThread) = particleLocation
                            energyLoss = energyLoss + self%energyThreshold(collIdx)
                        CASE(3)
                            velocity_CM = incidentVelocity
                            call self%elasticExcitCollisionIsotropic(particleList(self%reactantsIndx(1)), targetParticleList(self%reactantsIndx(2)), &
                                randGen(iThread), energyCM, self%energyThreshold(collIdx), incidentVelocity, targetVelocity)
                            particleList(self%reactantsIndx(1))%phaseSpace(2:4, particleIndx, iThread) = incidentVelocity
                            energyLoss = energyLoss + particleList(self%reactantsIndx(1))%mass * 0.5d0 * (SUM(velocity_CM**2) - SUM(incidentVelocity**2)) / e
                            !energyLoss = energyLoss + self%energyThreshold(collIdx)
                        CASE(4)
                            particleList(self%reactantsIndx(1))%phaseSpace(2:4, particleIndx, iThread) = targetVelocity
                            energyLoss = energyLoss + particleList(self%reactantsIndx(1))%mass * 0.5d0 * (SUM(incidentVelocity**2) - SUM(targetVelocity**2)) / e
                        CASE(5)
                            print *, 'dissociation'
                        END SELECT
                        exit
                    end if
                    sigma_v_low = sigma_v
                end if
            end do
            ! Swap with last particle so no repeats
            ! incidentVelocity = particleList(self%reactantsIndx(1))%phaseSpace(2:4, particleIndx, iThread)
            ! particleLocation = particleList(self%reactantsIndx(1))%phaseSpace(1, particleIndx, iThread)
            ! particleList(self%reactantsIndx(1))%phaseSpace(2:4, particleIndx, iThread) = particleList(self%reactantsIndx(1))%phaseSpace(2:4, numberTotalParticles, iThread)
            ! particleList(self%reactantsIndx(1))%phaseSpace(1, particleIndx, iThread) = particleList(self%reactantsIndx(1))%phaseSpace(1, numberTotalParticles, iThread)
            ! particleList(self%reactantsIndx(1))%phaseSpace(2:4, numberTotalParticles, iThread) = incidentVelocity
            ! particleList(self%reactantsIndx(1))%phaseSpace(1, numberTotalParticles, iThread) = particleLocation
            ! numberTotalParticles = numberTotalParticles - 1
        end do
        !$OMP end parallel
        self%totalEnergyLoss = self%totalEnergyLoss + energyLoss
        self%totalAmountCollisions = self%totalAmountCollisions + totalCollisions

        


    end subroutine generateCollision

    subroutine ionizationCollisionIsotropic(self, electron, ion, targetPart, randGen, energyCM, E_thres, incidentVelocity, targetVelocity, velocityCM)
        ! Ionization subroutine with momentum/energy conservation
        ! Replace incidentVelocity, velocityCM, and targetVelocity with primary electron, secondary electron, and ion velocity
        class(nullCollision), intent(in) :: self
        type(Particle), intent(in out) :: electron, ion
        type(targetParticle), intent(in) :: targetPart
        real(real64), intent(in) :: energyCM, E_thres
        real(real64), intent(in out) :: incidentVelocity(3), targetVelocity(3), velocityCM(3)
        type(mt19937), intent(in out) :: randGen
        real(real64) :: speedPerParticle, phi, e_vector(3), u_vector(3), V_cm(3), delE, secTheta, cos_theta, sin_theta, cos_phi, sin_phi, P_beginning(3)!, E_beginning, E_end, P_end(3)
        !integer(int32) :: i

        delE = (energyCM - E_thres)*e
        !E_beginning = m_e * 0.5d0 * SUM(incidentVelocity**2) + 0.5d0 * targetPart%mass * SUM(targetVelocity**2)
        P_beginning = m_e * incidentVelocity + targetPart%mass * targetVelocity
        
        V_cm = (P_beginning) / (self%sumMass)
    
        ! first add to electron
        cos_theta = 1.0d0 - 2.0d0*randGen%genrand64_real1()
        phi = randGen%genrand64_real1() * 2.0d0 * pi
        sin_theta = SQRT(1.0d0 - cos_theta**2)
        cos_phi = COS(phi)
        sin_phi = SIN(phi)
        e_vector(1) = cos_theta
        e_vector(2) = sin_phi * sin_theta
        e_vector(3) = cos_phi * sin_theta
        speedPerParticle = SQRT(2.0d0 * delE * self%reducedMassIonization / m_e**2)
        incidentVelocity = e_vector * speedPerParticle + V_cm

        ! second electron
        cos_theta = COS(2.0d0 * pi/3.0d0)
        phi = randGen%genrand64_real1() * 2.0d0 * pi

        call scatterVector(u_vector, e_vector, cos_theta, phi)

        ! Keep in case scatter Vector doesn't work properly
        ! cos_theta = 1.0d0 - 2.0d0*ran2(irand)
        ! phi = ran2(irand) * 2.0d0 * pi
        ! sin_theta = SQRT(1.0d0 - cos_theta**2)
        ! cos_phi = COS(phi)
        ! sin_phi = SIN(phi)
        ! y_vector(1) = cos_theta
        ! y_vector(2) = sin_phi * sin_theta
        ! y_vector(3) = cos_phi * sin_theta

        ! x_vector = crossProduct(e_vector, y_vector)
        ! x_vector = crossProduct(x_vector, e_vector)/SQRT(SUM(x_vector**2))
       
        ! cos_theta = COS(2.0d0 * pi/3.0d0)
        ! sin_theta = SIN(2.0d0 * pi / 3.0d0)
        ! u_vector = e_vector * cos_theta + x_vector * sin_theta
        ! print *, 'u vector total:', SUM(u_vector**2)
        ! print *, 'dot product u_vector and e_vector:', SUM(u_vector * e_vector)
        ! print *, 'should be:', cos_theta
        ! print *, 'u_vector is:', u_vector


        ! secTheta = ACOS(cos_theta) + 2.0d0 * pi/3.0d0
        ! u_vector(1) = COS(secTheta)
        ! u_vector(2) = sin_phi * SIN(secTheta)
        ! u_vector(3) = cos_phi * SIN(secTheta)
        velocityCM = u_vector * speedPerParticle + V_cm
    
        ! ion
        targetVelocity = (P_beginning - m_e * (incidentVelocity + velocityCM))/ion%mass
        ! E_end = m_e * 0.5d0 * SUM(incidentVelocity**2) + m_e * 0.5d0 * SUM(velocityCM**2) + 0.5d0 * ion%mass * SUM(targetVelocity**2)
        ! P_end = m_e * incidentVelocity + m_e * velocityCM + ion%mass * targetVelocity
        
        
        ! if (ABS((E_beginning - E_end - E_thres*e)/E_beginning) > 1.d-8) then
        !     print *, 'issue energy conservation:'
        !     print *, 'E_beginning:', E_beginning
        !     print *, 'E_end:', E_end + E_thres*e
        !     stop
        ! end if
        ! do i = 1, 3
        !     if (ABS((P_beginning(i) - P_end(i))/P_beginning(i)) > 1.d-8) then
        !         print *, 'issue momentum conservation:'
        !         print *, 'momentum before:', P_beginning
        !         print *, 'momentum after:', P_end
        !         stop
        !     end if
        ! end do

        
    end subroutine ionizationCollisionIsotropic

    subroutine ionizationCollisionIsotropicNanbul(self, electron, ion, targetPart, randGen, energyCM, E_thres, incidentVelocity, targetVelocity, velocityCM)
        ! Ionization subroutine with only approximate energy/momentum conservation
        ! Replace incidentVelocity, velocityCM, and targetVelocity with primary electron, secondary electron, and ion velocity
        class(nullCollision), intent(in) :: self
        type(Particle), intent(in out) :: electron, ion
        type(targetParticle), intent(in) :: targetPart
        real(real64), intent(in) :: energyCM, E_thres
        real(real64), intent(in out) :: incidentVelocity(3), targetVelocity(3), velocityCM(3)
        type(mt19937), intent(in out) :: randGen
        real(real64) :: speedPerParticle, phi, e_vector(3), V_cm(3), delE, cos_theta, sin_theta, cos_phi, sin_phi, P_beginning(3)!, E_beginning, E_end, P_end(3)
        !integer(int32) :: i

        delE = (energyCM - E_thres)*e
        !E_beginning = m_e * 0.5d0 * SUM(incidentVelocity**2) + 0.5d0 * targetPart%mass * SUM(targetVelocity**2)
        P_beginning = m_e * incidentVelocity + targetPart%mass * targetVelocity
        
        V_cm = (P_beginning) / (self%sumMass)
    
        ! first add to electron
        cos_theta = 1.0d0 - 2.0d0*randGen%genrand64_real1()
        phi = randGen%genrand64_real1() * 2.0d0 * pi
        sin_theta = SQRT(1.0d0 - cos_theta**2)
        cos_phi = COS(phi)
        sin_phi = SIN(phi)
        e_vector(1) = cos_theta
        e_vector(2) = sin_phi * sin_theta
        e_vector(3) = cos_phi * sin_theta
        speedPerParticle = SQRT(delE/m_e)
        incidentVelocity = e_vector * speedPerParticle + V_cm

        ! second electron
        cos_theta = 1.0d0 - 2.0d0*randGen%genrand64_real1()
        phi = randGen%genrand64_real1() * 2.0d0 * pi
        sin_theta = SQRT(1.0d0 - cos_theta**2)
        cos_phi = COS(phi)
        sin_phi = SIN(phi)
        e_vector(1) = cos_theta
        e_vector(2) = sin_phi * sin_theta
        e_vector(3) = cos_phi * sin_theta
        velocityCM = e_vector * speedPerParticle + V_cm
    
        ! ion
        targetVelocity = (P_beginning - m_e * (incidentVelocity + velocityCM))/ion%mass
        
        ! E_end = m_e * 0.5d0 * SUM(incidentVelocity**2) + m_e * 0.5d0 * SUM(velocityCM**2) + 0.5d0 * ion%mass * SUM(targetVelocity**2)
        ! P_end = m_e * incidentVelocity + m_e * velocityCM + ion%mass * targetVelocity
        
        
        ! if (ABS((E_beginning - E_end - E_thres*e)/E_beginning) > 1.d-8) then
        !     print *, 'issue energy conservation:'
        !     print *, 'E_beginning:', E_beginning
        !     print *, 'E_end:', E_end + E_thres*e
        !     stop
        ! end if
        ! do i = 1, 3
        !     if (ABS((P_beginning(i) - P_end(i))/P_beginning(i)) > 1.d-8) then
        !         print *, 'issue momentum conservation:'
        !         print *, 'momentum before:', P_beginning
        !         print *, 'momentum after:', P_end
        !         stop
        !     end if
        ! end do

        
    end subroutine ionizationCollisionIsotropicNanbul

    subroutine elasticExcitCollisionIsotropic(self, primary, targetPart, randGen, energyCM, E_thres, incidentVelocity, targetVelocity)
        ! elastic collision routine with momentum/energy conservation
        class(nullCollision), intent(in) :: self
        type(Particle), intent(in out) :: primary
        type(targetParticle), intent(in) :: targetPart
        real(real64), intent(in) :: energyCM, E_thres
        real(real64), intent(in out) :: incidentVelocity(3), targetVelocity(3)
        type(mt19937), intent(in out) :: randGen
        real(real64) :: speedPerParticle, phi, e_vector(3), V_cm(3), delE, cos_theta, sin_theta, cos_phi, sin_phi, secTheta, P_beginning(3)!, E_beginning, E_end, P_end(3)
        !integer(int32) :: i

        delE = (energyCM - E_thres)*e
        !E_beginning = primary%mass * 0.5d0 * SUM(incidentVelocity**2) + 0.5d0 * targetPart%mass * SUM(targetVelocity**2)
        P_beginning = primary%mass * incidentVelocity + targetPart%mass * targetVelocity
        
        V_cm = (P_beginning) / self%sumMass
    
        ! first add to primary
        cos_theta = 1.0d0 - 2.0d0*randGen%genrand64_real1()
        phi = randGen%genrand64_real1() * 2.0d0 * pi
        sin_theta = SQRT(1.0d0 - cos_theta**2)
        cos_phi = COS(phi)
        sin_phi = SIN(phi)
        e_vector(1) = cos_theta
        e_vector(2) = sin_phi * sin_theta
        e_vector(3) = cos_phi * sin_theta
        speedPerParticle = SQRT(2.0d0 * delE * self%reducedMass / primary%mass**2)
        incidentVelocity = e_vector * speedPerParticle + V_cm

        ! second primary
    
        targetVelocity = (P_beginning - primary%mass * incidentVelocity)/targetPart%mass

        ! E_end = primary%mass * 0.5d0 * SUM(incidentVelocity**2) + targetPart%mass * 0.5d0 * SUM(targetVelocity**2)
        ! P_end = primary%mass * incidentVelocity + targetPart%mass * targetVelocity

        ! if (ABS((E_beginning - E_end - E_thres*e)/E_beginning) > 1.d-8) then
        !     print *, 'issue energy conservation:'
        !     print *, 'E_beginning:', E_beginning
        !     print *, 'E_end:', E_end + E_thres*e
        !     stop
        ! end if
        
        ! do i = 1, 3
        !     if (ABS((P_beginning(i) - P_end(i))/P_beginning(i)) > 1.d-8) then
        !         print *, 'issue momentum conservation:'
        !         print *, 'momentum before:', P_beginning
        !         print *, 'momentum after:', P_end
        !         stop
        !     end if
        ! end do

    
    end subroutine elasticExcitCollisionIsotropic

    subroutine writeCollisionCrossSection(self, part, dirName)
        class(nullCollision), intent(in) :: self
        type(Particle), intent(in) :: part
        character(*), intent(in) :: dirName
        integer(int32) :: i
        open(10,file=dirName//'/CrossSections/IncidentPart_'//part%name//"_energy.dat", form='UNFORMATTED', access = 'STREAM')
        write(10) self%energyArray
        close(10)

        open(10,file=dirName//'/CrossSections/IncidentPart_'//part%name//"_sigma.dat", form='UNFORMATTED', access = 'STREAM')
        write(10) self%sigmaArray
        close(10)

    end subroutine writeCollisionCrossSection

    subroutine scatterVector(u,e,costheta,phi)
        !     ===================================================================
        !     VERSION:         0.1
        !     LAST MOD:      April/19
        !     MOD AUTHOR:    G. Hagelaar
        !     COMMENTS:     
        !     NOTES:   Turn vector vx,vy,vz over a scattering angle theta and azimuthal angle phi
        !     -------------------------------------------------------------------
        real(real64), intent(in) :: costheta, phi
        real(real64), intent(in) :: e(3)
        real(real64), intent(in out) :: u(3)
        real(real64) :: sintheta,cosphi,sinphi, vv
        
        sintheta=SQRT(1d0-costheta**2)
        sinphi=dsin(phi)
        cosphi=dcos(phi)
        IF (dabs(e(2))>dabs(e(3))) THEN
            vv=dsqrt(e(1)**2+e(2)**2)
            u(1)=e(1)*costheta+(e(2)*sinphi+e(1)*e(3)*cosphi)/vv*sintheta
            u(2)=e(2)*costheta+(-e(1)*sinphi+e(2)*e(3)*cosphi)/vv*sintheta
            u(3)=e(3)*costheta-vv*cosphi*sintheta
        ELSE
            vv=dsqrt(e(1)**2+e(3)**2)
            u(1)=e(1)*costheta+(e(3)*sinphi-e(2)*e(1)*cosphi)/vv*sintheta
            u(2)=e(2)*costheta+vv*cosphi*sintheta
            u(3)=e(3)*costheta-(e(1)*sinphi+e(2)*e(3)*cosphi)/vv*sintheta
        ENDIF
    end subroutine scatterVector


end module mod_NullCollision