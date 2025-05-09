module mod_NullCollision

    use iso_fortran_env, only: int32, real64, output_unit
    use constants
    use mod_BasicFunctions
    use mod_particle
    use mod_targetParticle
    use omp_lib
    use iso_c_binding
    use ifport, only: makedirqq
    implicit none

    private
    public :: nullCollision, readNullCollisionInputs
    integer(int32), public, protected :: numberBinaryCollisions = 0
    real(real64), parameter :: sin_third_rot = sin(2.0d0 * pi / 3.0d0), cos_third_rot = cos(2.0d0 * pi / 3.0d0)
    

    ! Stores information for each binary collision of charged particle on target particle
    ! Uses null collision method to randomly select only a portion of particles to test for collision
    type :: nullCollision
        integer(int32) :: numberCollisions, numberReactants, lengthArrays ! total number collisions, reactant number, length of energy/cross section arrays
        real(real64), allocatable :: energyArray(:), sigmaArray(:, :), energyThreshold(:), totalIncidentEnergy(:), totalEnergyLoss(:)
        real(real64) :: sigmaVMax, reducedMass, minEnergy, maxEnergy, sumMass, reducedMassIonization 
        integer(int32), allocatable :: collisionType(:), reactantsIndx(:), numberProducts(:), productsIndx(:,:)
        integer(int64), allocatable :: totalAmountCollisions(:)
        integer(int64) :: totalNumberCollidableParticles
    contains
        procedure, public, pass(self) :: generateCollision
        procedure, public, pass(self) :: elasticExcitCollisionIsotropic
        procedure, public, pass(self) :: ionizationCollisionIsotropic
        procedure, public, pass(self) :: writeCollisionProperties
        procedure, public, pass(self) :: writeCollisionCrossSection
        procedure, public, pass(self) :: writeCollisionDiag


    end type nullCollision


    interface nullCollision 
        procedure :: nullCollision_constructor
    end interface nullCollision

contains

    type(nullCollision) function nullCollision_constructor(numberReactants, numberCollisions, lengthArrays, red_mass, sumMass, redMassTripleProducts, energyArray, sigmaArray, energyThreshold, collisionType, reactantsIndx, numberProducts, productsIndx) result(self)
        ! Construct null Collision object
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
            self%reactantsIndx(numberReactants), self%numberProducts(numberCollisions), self%productsIndx(3, numberCollisions), self%totalAmountCollisions(numberCollisions), &
            self%totalIncidentEnergy(numberCollisions), self%totalEnergyLoss(numberCollisions))
        self%energyArray = energyArray
        self%reactantsIndx = reactantsIndx
        do i = 1, numberCollisions
            ! Go top down in order, assign arrays. This is just to put more 'likely' collision first when looking at collision probabilities
            j = indxOrder(i)
            ! array of sigma for each collision
            self%sigmaArray(:, i) = sigmaArray(:, j)
            ! lower energy threshold
            self%energyThreshold(i) = energyThreshold(j)
            ! collision type
            self%collisionType(i) = collisionType(j)
            ! number of products
            self%numberProducts(i) = numberProducts(j)
            ! identifying indx for each product (in particle list)
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
        self%totalAmountCollisions = 0.0d0
        self%totalNumberCollidableParticles = 0
        self%totalIncidentEnergy = 0.0d0
        self%totalEnergyLoss = 0.0d0
        self%minEnergy = MINVAL(self%energyArray)
        self%maxEnergy = MAXVAL(self%energyArray)
    
    end function nullCollision_constructor

    subroutine generateCollision(self, particleList, targetParticleList, numberChargedParticles, numberBinaryCollisions, irand, del_t)
        ! Binary collision operator on particle with background gas
        class(nullCollision), intent(in out) :: self
        integer(int32), intent(in) :: numberChargedParticles, numberBinaryCollisions
        type(Particle), intent(in out) :: particleList(numberChargedParticles)
        type(targetParticle), intent(in) :: targetParticleList(numberBinaryCollisions)
        integer(c_int64_t), intent(in out) :: irand(numThread)
        real(real64), intent(in) :: del_t
        logical :: collisionBool
        real(real64) :: P_null, numberSelectedReal, Rand, targetVelocity(3), incidentVelocity(3), velocity_CM(3), energyCM, d_value, sigma_v, &
        sigma_v_low, speedCM, particleLocation, energyLoss(self%numberCollisions), primary_mass, target_mass, IncidentEnergy, totIncidentEnergy(self%numberCollisions)
        integer(int32) :: numberSelected, iThread, i, particleIndx, numberTotalParticles, indxHigh, indxLow, indxMiddle, collIdx, addIonizationIndx
        integer(int64) :: totalCollisions(self%numberCollisions)
        integer(c_int64_t) :: irand_thread

        ! Calculate null collision probability to determine percentage of particles to react
        P_null = 1.0d0 - EXP(-self%sigmaVMax * targetParticleList(self%reactantsIndx(2))%density * del_t)
    
        if (P_null > 0.1d0) then
            print *, 'P_null greater than 10% at ', P_null*100, '%'
            stop
        end if
        ! Initialize diagnostics
        energyLoss = 0
        totalCollisions = 0
        totIncidentEnergy = 0

        ! collision of primary particle on target particle

        primary_mass = particleList(self%reactantsIndx(1))%mass
        target_mass = targetParticleList(self%reactantsIndx(2))%mass
        self%totalNumberCollidableParticles = self%totalNumberCollidableParticles + SUM(particleList(self%reactantsIndx(1))%numToCollide) ! accumulate number of particles collided
        !$OMP parallel private(iThread, i,numberSelected, numberSelectedReal, Rand, particleIndx, numberTotalParticles, targetVelocity, incidentVelocity, &
        !$OMP&    velocity_CM, energyCM, d_value, indxHigh, indxLow, indxMiddle, collIdx, sigma_v, sigma_v_low, collisionBool, speedCM, &
        !$OMP&    addIonizationIndx, particleLocation, irand_thread, IncidentEnergy) reduction(+:energyLoss,totalCollisions, totIncidentEnergy)
        iThread = omp_get_thread_num() + 1
        irand_thread = irand(iThread) ! set to local variable to prevent false sharing
        ! Select amount of particles to sample
        numberTotalParticles = particleList(self%reactantsIndx(1))%numToCollide(iThread)
        numberSelectedReal = P_null * real(numberTotalParticles)
        numberSelected = INT(numberSelectedReal)
        ! Account for integer selection
        Rand = pcg32_random_r(irand_thread)
        if (Rand < (numberSelectedReal - numberSelected)) numberSelected = numberSelected + 1   
        addIonizationIndx = 0
        do i = 1, numberSelected
            ! Randomly sample particle and save properties
            Rand = pcg32_random_r(irand_thread)
            particleIndx = INT((numberTotalParticles) * Rand) + 1
            particleLocation = particleList(self%reactantsIndx(1))%phaseSpace(1, particleIndx, iThread)
            incidentVelocity = particleList(self%reactantsIndx(1))%phaseSpace(2:4, particleIndx, iThread)

            ! Generate velocity for target particle
            targetVelocity = targetParticleList(self%reactantsIndx(2))%generate3DMaxwellianVelocity(irand_thread)

            ! Get relative velocity and center of mass values
            velocity_CM = incidentVelocity - targetVelocity
            speedCM = SQRT(SUM(velocity_CM**2))
            energyCM = SUM(velocity_CM**2) * 0.5d0 * self%reducedMass / e !in eV
            indxLow = 1
            indxHigh = self%lengthArrays
            ! get index and linear interpolation value (d_value) for sigma array
            if (energyCM < self%minEnergy) then
                ! take minimum sigma
                d_value = 0.0d0
            else if (energyCM > self%maxEnergy) then
                ! take maximum sigma
                indxLow = indxHigh -1
                d_value = 1.0d0
            else
                ! binary search to find location of CM energy in energy array
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
            Rand = pcg32_random_r(irand_thread)
            sigma_v_low = 0.0d0
            do collIdx = 1, self%numberCollisions
                ! Cycle through each collision and test if in probabilty bracket
                if (energyCM > self%energyThreshold(collIdx)) then
                    ! linear interpolation of sigma between lower and higher index
                    sigma_v = (self%sigmaArray(indxLow, collIdx) * (1.0d0 - d_value) + self%sigmaArray(indxHigh, collIdx) * d_value) * speedCM  + sigma_v_low ! accumulate sigma_v
                    collisionBool = (Rand <= sigma_v/self%sigmaVMax) 
                    if (collisionBool) then
                        ! Collide particle
                        totalCollisions(collIdx) = totalCollisions(collIdx) + 1
                        IncidentEnergy = SUM(incidentVelocity**2)
                        totIncidentEnergy(collIdx) = totIncidentEnergy(collIdx) + IncidentEnergy
                        ! Collision type
                        SELECT CASE (self%collisionType(collIdx))
                        CASE(1)
                            ! elastic collision
                            call self%elasticExcitCollisionIsotropic(primary_mass, target_mass, &
                                irand_thread, energyCM, self%energyThreshold(collIdx), incidentVelocity, targetVelocity)
                        CASE(2)
                            ! ionization
                            call self%ionizationCollisionIsotropic(primary_mass, particleList(self%productsIndx(2, collIdx))%mass, target_mass, &
                                irand_thread, energyCM, self%energyThreshold(collIdx), incidentVelocity, targetVelocity, velocity_CM)
                            !Update total particle count
                            particleList(self%reactantsIndx(1))%N_p(iThread) = particleList(self%reactantsIndx(1))%N_p(iThread) + 1
                            particleList(self%productsIndx(2, collIdx))%N_p(iThread) = particleList(self%productsIndx(2, collIdx))%N_p(iThread) + 1
                            ! set new secondary particle (should be electron) velocity and position
                            particleList(self%reactantsIndx(1))%phaseSpace(2:4, particleList(self%reactantsIndx(1))%N_p(iThread), iThread) = velocity_CM
                            particleList(self%reactantsIndx(1))%phaseSpace(1, particleList(self%reactantsIndx(1))%N_p(iThread), iThread) = particleLocation
                            ! set new ion velocity and position
                            particleList(self%productsIndx(2, collIdx))%phaseSpace(2:4, particleList(self%productsIndx(2, collIdx))%N_p(iThread), iThread) = targetVelocity
                            particleList(self%productsIndx(2, collIdx))%phaseSpace(1, particleList(self%productsIndx(2, collIdx))%N_p(iThread), iThread) = particleLocation
                            energyLoss(collIdx) = energyLoss(collIdx) - SUM(targetVelocity**2) * target_mass - SUM(velocity_CM**2) * primary_mass
                        CASE(3)
                            ! excitation
                            call self%elasticExcitCollisionIsotropic(primary_mass, target_mass, &
                                irand_thread, energyCM, self%energyThreshold(collIdx), incidentVelocity, targetVelocity)
                        CASE(4)
                            ! charge exchange
                            incidentVelocity = targetVelocity
                        CASE(5)
                            continue
                        END SELECT
                        ! collect energy loss statistics
                        energyLoss(collIdx) = energyLoss(collIdx) + primary_mass * (IncidentEnergy - SUM(incidentVelocity**2))
                        exit
                    end if
                    sigma_v_low = sigma_v
                end if
            end do
            ! Swap with last particle so no repeats, incident velocity is new velocity for primary particle chosen, then lower amount of particles to choose from
            particleList(self%reactantsIndx(1))%phaseSpace(2:4, particleIndx, iThread) = particleList(self%reactantsIndx(1))%phaseSpace(2:4, numberTotalParticles, iThread)
            particleList(self%reactantsIndx(1))%phaseSpace(1, particleIndx, iThread) = particleList(self%reactantsIndx(1))%phaseSpace(1, numberTotalParticles, iThread)
            particleList(self%reactantsIndx(1))%phaseSpace(2:4, numberTotalParticles, iThread) = incidentVelocity
            particleList(self%reactantsIndx(1))%phaseSpace(1, numberTotalParticles, iThread) = particleLocation
            numberTotalParticles = numberTotalParticles - 1
        end do
        irand(iThread) = irand_thread
        particleList(self%reactantsIndx(1))%numToCollide(iThread) = numberTotalParticles
        !$OMP end parallel
        ! accumulate statistics
        self%totalAmountCollisions = self%totalAmountCollisions + totalCollisions
        self%totalEnergyLoss = self%totalEnergyLoss + energyLoss
        self%totalIncidentEnergy = self%totalIncidentEnergy + totIncidentEnergy
        

    end subroutine generateCollision


    subroutine ionizationCollisionIsotropic(self, primary_mass, ion_mass, target_mass, irand, energyCM, E_thres, incidentVelocity, targetVelocity, velocityCM)
        ! Ionization subroutine with momentum/energy conservation
        ! Replace incidentVelocity, velocityCM, and targetVelocity with primary electron, secondary electron, and ion velocity
        ! theta azimuthal angle, phi is inclination
        class(nullCollision), intent(in) :: self
        real(real64), intent(in) :: energyCM, E_thres, primary_mass, ion_mass, target_mass
        real(real64), intent(in out) :: incidentVelocity(3), targetVelocity(3), velocityCM(3)
        integer(c_int64_t), intent(in out) :: irand
        real(real64) :: speedPerParticle, phi, e_vector(3), u_vector(3), y_vector(3), &
        V_cm(3), delE, sin_theta_new, cos_theta_new, cos_theta, sin_theta, cos_phi, sin_phi, P_beginning(3)!, E_beginning, E_end, P_end(3)

        !get center of mass quantities
        delE = (energyCM - E_thres)*e
        P_beginning = primary_mass * incidentVelocity + target_mass * targetVelocity
        V_cm = (P_beginning) / (self%sumMass)
    
        ! Scatter primary particle isotropically in center of mass frame
        cos_theta = 1.0d0 - 2.0d0*pcg32_random_r(irand)
        phi = pcg32_random_r(irand) * 2.0d0 * pi
        sin_theta = SQRT(1.0d0 - cos_theta**2)
        cos_phi = COS(phi)
        sin_phi = SIN(phi)
        e_vector(1) = cos_phi * sin_theta ! x
        e_vector(2) = sin_phi * sin_theta ! y
        e_vector(3) = cos_theta ! z
        speedPerParticle = SQRT(2.0d0 * delE * self%reducedMassIonization / primary_mass**2)
        incidentVelocity = e_vector * speedPerParticle + V_cm

        ! get rotation vector rotating by third around theta
        cos_theta_new = cos_third_rot * cos_theta - sin_third_rot * sin_theta
        sin_theta_new = cos_theta * sin_third_rot + cos_third_rot * sin_theta
        y_vector(1) = cos_phi * sin_theta_new ! x
        y_vector(2) = sin_phi * sin_theta_new ! y
        y_vector(3) = cos_theta_new ! z
        
        phi = pcg32_random_r(irand) * 2.0d0 * pi
        cos_phi = cos(phi)
        sin_phi = sin(phi)

        ! Rotate angle randomly about the primary particle vector
        ! Source is "Rodrigues' rotation formula" on wikipedia
        u_vector = -e_vector * cos_third_rot * (cos_phi - 1.0d0) + y_vector* cos_phi + crossProduct(e_vector, y_vector) * sin_phi

        
        ! velocity of produced electron
        speedPerParticle = SQRT(2.0d0 * delE * self%reducedMassIonization / m_e**2)
        velocityCM = u_vector * speedPerParticle + V_cm
    
        ! ion velocity
        targetVelocity = (P_beginning - primary_mass * incidentVelocity - m_e * velocityCM)/ion_mass

        
    end subroutine ionizationCollisionIsotropic

    ! subroutine ionizationCollisionIsotropicNanbul(self, primary_mass, ion_mass, target_mass, irand, energyCM, E_thres, incidentVelocity, targetVelocity, velocityCM)
    !     ! Ionization subroutine with only approximate energy/momentum conservation
    !     ! Replace incidentVelocity, velocityCM, and targetVelocity with primary electron, secondary electron, and ion velocity
    !     class(nullCollision), intent(in) :: self
    !     real(real64), intent(in) :: energyCM, E_thres, primary_mass, ion_mass, target_mass
    !     real(real64), intent(in out) :: incidentVelocity(3), targetVelocity(3), velocityCM(3)
    !     integer(int64), intent(in out) :: irand
    !     real(real64) :: speedPerParticle, phi, e_vector(3), V_cm(3), delE, cos_theta, sin_theta, cos_phi, sin_phi, P_beginning(3)!, E_beginning, E_end, P_end(3)
    !     !integer(int32) :: i

    !     delE = (energyCM - E_thres)*e
    !     !E_beginning = m_e * 0.5d0 * SUM(incidentVelocity**2) + 0.5d0 * targetPart%mass * SUM(targetVelocity**2)
    !     P_beginning = primary_mass * incidentVelocity + target_mass * targetVelocity
        
    !     V_cm = (P_beginning) / (self%sumMass)
    
    !     ! first add to primary
    !     cos_theta = 1.0d0 - 2.0d0*randPCG(irand)
    !     phi = randPCG(irand) * 2.0d0 * pi
    !     sin_theta = SQRT(1.0d0 - cos_theta**2)
    !     cos_phi = COS(phi)
    !     sin_phi = SIN(phi)
    !     e_vector(1) = cos_theta
    !     e_vector(2) = sin_phi * sin_theta
    !     e_vector(3) = cos_phi * sin_theta
    !     speedPerParticle = SQRT(delE/primary_mass)
    !     incidentVelocity = e_vector * speedPerParticle + V_cm

    !     ! second electron
    !     cos_theta = 1.0d0 - 2.0d0*randPCG(irand)
    !     phi = randPCG(irand) * 2.0d0 * pi
    !     sin_theta = SQRT(1.0d0 - cos_theta**2)
    !     cos_phi = COS(phi)
    !     sin_phi = SIN(phi)
    !     e_vector(1) = cos_theta
    !     e_vector(2) = sin_phi * sin_theta
    !     e_vector(3) = cos_phi * sin_theta
    !     speedPerParticle = SQRT(delE/m_e)
    !     velocityCM = e_vector * speedPerParticle + V_cm
    
    !     ! ion
    !     targetVelocity = (P_beginning - primary_mass * incidentVelocity - m_e*velocityCM)/ion_mass
        
    !     ! E_end = m_e * 0.5d0 * SUM(incidentVelocity**2) + m_e * 0.5d0 * SUM(velocityCM**2) + 0.5d0 * ion%mass * SUM(targetVelocity**2)
    !     ! P_end = m_e * incidentVelocity + m_e * velocityCM + ion%mass * targetVelocity
        
        
    !     ! if (ABS((E_beginning - E_end - E_thres*e)/E_beginning) > 1.d-8) then
    !     !     print *, 'issue energy conservation:'
    !     !     print *, 'E_beginning:', E_beginning
    !     !     print *, 'E_end:', E_end + E_thres*e
    !     !     stop
    !     ! end if
    !     ! do i = 1, 3
    !     !     if (ABS((P_beginning(i) - P_end(i))/P_beginning(i)) > 1.d-8) then
    !     !         print *, 'issue momentum conservation:'
    !     !         print *, 'momentum before:', P_beginning
    !     !         print *, 'momentum after:', P_end
    !     !         stop
    !     !     end if
    !     ! end do

        
    ! end subroutine ionizationCollisionIsotropicNanbul

    subroutine elasticExcitCollisionIsotropic(self, primary_mass, target_mass, irand, energyCM, E_thres, incidentVelocity, targetVelocity)
        ! elastic collision routine with momentum/energy conservation
        class(nullCollision), intent(in) :: self
        real(real64), intent(in) :: energyCM, E_thres, primary_mass, target_mass
        real(real64), intent(in out) :: incidentVelocity(3), targetVelocity(3)
        integer(c_int64_t), intent(in out) :: irand
        real(real64) :: speedPerParticle, phi, e_vector(3), V_cm(3), delE, cos_theta, sin_theta, cos_phi, sin_phi, secTheta, P_beginning(3)!, E_beginning, E_end, P_end(3)

        delE = (energyCM - E_thres)*e
        P_beginning = primary_mass * incidentVelocity + target_mass * targetVelocity
        
        V_cm = (P_beginning) / self%sumMass
    
        ! scatter primary particle in CM frame
        cos_theta = 1.0d0 - 2.0d0*pcg32_random_r(irand)
        phi = pcg32_random_r(irand) * 2.0d0 * pi
        sin_theta = SQRT(1.0d0 - cos_theta**2)
        cos_phi = COS(phi)
        sin_phi = SIN(phi)
        
        
        e_vector(1) = cos_phi * sin_theta ! x
        e_vector(2) = sin_phi * sin_theta ! y
        e_vector(3) = cos_theta ! z
        speedPerParticle = SQRT(2.0d0 * delE * self%reducedMass / primary_mass**2)
        incidentVelocity = e_vector * speedPerParticle + V_cm

        ! scatter secondary particle
    
        targetVelocity = (P_beginning - primary_mass * incidentVelocity)/target_mass

    
    end subroutine elasticExcitCollisionIsotropic

    subroutine writeCollisionCrossSection(self, part, dirName)
        ! write cross sections generated
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

    subroutine writeCollisionProperties(self, particleList, targetParticleList, dirName)
        ! Generate folders in save directory for binary collisions
        class(nullCollision), intent(in) :: self
        type(Particle), intent(in out) :: particleList(numberChargedParticles)
        type(targetParticle), intent(in) :: targetParticleList(numberBinaryCollisions)
        character(*), intent(in) :: dirName
        logical :: bool
        integer(int32) :: i, maxIndx(1)
        character(len=5) :: char_i

        bool = makedirqq(dirName//'/BinaryCollisions/'//particleList(self%reactantsIndx(1))%name//"_on_"//targetParticleList(self%reactantsIndx(2))%name)
        if (.not. bool) then
            print *, 'Not generating folder', dirName//'/BinaryCollisions/'//particleList(self%reactantsIndx(1))%name//"_on_"//targetParticleList(self%reactantsIndx(2))%name
            stop "Save directory binary collision not successfully created!"
        end if
        open(10,file=dirName//'/BinaryCollisions/'//particleList(self%reactantsIndx(1))%name//"_on_"//targetParticleList(self%reactantsIndx(2))%name//"/CollisionProperties.dat")
        write(10, '("Coll #, collType, E_thres (eV), maxSigma (m^2), EatMaxSigma (eV)")')
        do i = 1, self%numberCollisions
            maxIndx = MAXLOC(self%sigmaArray(:, i))
            write(10,"(2(I3, 1x), 3(es16.8,1x))") i, self%collisionType(i), self%energyThreshold(i), self%sigmaArray(maxIndx(1), i), self%energyArray(maxIndx(1))
        end do
        close(10)

        ! For each collision create diagnostic file
        do i = 1, self%numberCollisions
            write(char_i, '(I3)'), i
            open(10,file=dirName//'/BinaryCollisions/'//particleList(self%reactantsIndx(1))%name//"_on_"//targetParticleList(self%reactantsIndx(2))%name//"/CollisionDiag_"//trim(adjustl(char_i))//".dat")
            write(10, '("CollRatio, AveEnergyLoss (eV), AveIncidentEnergy (eV), P_loss(W/m^2), aveCollFreq (Hz/m^2)")')
            close(10)
        end do

    end subroutine writeCollisionProperties

    subroutine writeCollisionDiag(self, particleList, targetParticleList, dirName, timeInterval)
        class(nullCollision), intent(in) :: self
        type(Particle), intent(in out) :: particleList(numberChargedParticles)
        type(targetParticle), intent(in) :: targetParticleList(numberBinaryCollisions)
        character(*), intent(in) :: dirName
        real(real64), intent(in) :: timeInterval
        integer(int32) :: i, maxIndx(1)
        character(len=5) :: char_i

        ! For each collision add to diag file
        do i = 1, self%numberCollisions
            write(char_i, '(I3)'), i
            open(10,file=dirName//'/BinaryCollisions/'//particleList(self%reactantsIndx(1))%name//"_on_"//targetParticleList(self%reactantsIndx(2))%name//"/CollisionDiag_"//trim(adjustl(char_i))//".dat", &
                STATUS = 'OLD', ACCESS = 'APPEND')
            write(10,"(5(es16.8,1x))") real(self%totalAmountCollisions(i))/real(self%totalNumberCollidableParticles), self%totalEnergyLoss(i) * 0.5d0 / e/ real(self%totalAmountCollisions(i)), &
            self%totalIncidentEnergy(i) * particleList(self%reactantsIndx(1))%mass * 0.5d0 / e/ real(self%totalAmountCollisions(i)), &
            self%totalEnergyLoss(i) * 0.5d0 * particleList(self%reactantsIndx(1))%w_p / timeInterval, &
            real(self%totalAmountCollisions(i)) * particleList(self%reactantsIndx(1))%w_p / timeInterval
            close(10)
        end do

    end subroutine writeCollisionDiag

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

    subroutine rotate_vector(vector_final, rotation_vector, rotation_axis,cos_theta_angle,phi)
        ! Remake scatter algorithm for clarity
        ! Input is e which is initial vector to scatter by angle theta
        ! theta input is already in cos and sin since might be predetermined input
        ! Basically make vector with angle theta to e simply by rotating along theta
        ! Have initial 'base' vector u with vector which to rotate by angle phi
        ! Source is "Rodrigues' rotation formula" on wikipedia
        real(real64), intent(in) :: cos_theta_angle, phi
        real(real64), intent(in) :: rotation_vector(3), rotation_axis(3)
        real(real64), intent(in out) :: vector_final(3)
        real(real64) :: cos_phi, sin_phi
        
        sin_phi= sin(phi)
        cos_phi= cos(phi)

        vector_final = -rotation_axis * cos_theta_angle * (cos_phi - 1.0d0) + rotation_vector* cos_phi + crossProduct(rotation_axis, rotation_vector) * sin_phi
    end subroutine rotate_vector


    ! --------------------------- Read in collisions --------------------------------------------

    subroutine readNullCollisionInputs(filename, nullCollisionList, particleList, targetParticleList)
        ! Reaction is any unique set of particle on target. For each reaction, you may have several different collision types
        character(len=*), intent(in) :: filename
        type(nullCollision), allocatable, intent(out) :: nullCollisionList(:)
        type(Particle), intent(in) :: particleList(numberChargedParticles)
        type(targetParticle), intent(in) :: targetParticleList(numberNeutralParticles)
        integer(int32), parameter :: maxNumberPoints = 1500, maxNumberColls = 65
        logical :: fileExists, foundParticleBool, reactantBool, incidentParticleBool, oldSetBool
        character(len=100) :: string
        character(:), allocatable :: collFileName
        real(real64) :: tempVar1, tempVar2, E_scaling, sigma_scaling, E_threshold(maxNumberColls), E_temp(maxNumberPoints, maxNumberColls), sigma_temp(maxNumberPoints, maxNumberColls), red_mass, sumMass, redMassTripleProducts
        integer(int32) :: i, j, k, numberCollisions, reactionIndx(maxNumberColls), lowerIndx, higherIndx, numberReactants(maxNumberColls), numberProducts(maxNumberColls), collisionReactantsIndx(2, maxNumberColls), &
            collisionProductsIndx(3, maxNumberColls), numberSigmaPoints(maxNumberColls), length_concatenate, collisionType(maxNumberColls), h, numberCollisionsPerReaction
        real(real64), allocatable :: E_concatenate_temp(:), E_concatenate(:), sigma_concatenate(:, :), E_threshold_array(:)
        integer(int32), allocatable :: indxSort(:), collisionTypeArray(:), numberProductsArray(:), productsIndxArray(:,:)
        

        print *, "Reading collision inputs:"
        print *, '---------------------------'
        numberBinaryCollisions = 0
        numberCollisions = 0
        if (.not. restartBool) then
            open(10,file='../InputData/'//filename, action = 'read')
        else
            open(10,file=restartDirectory//'/InputData/'//filename, action = 'read')
        end if
        read(10,*) string
        do while (.not. (string(1:3) == 'END' .or. string(1:3) == 'end'))
            collFileName = trim(string)
            inquire(file = '../../CollisionData/'//collFileName, exist = fileExists)
            ! If collision data file exists
            if (fileExists) then
                print *, 'Taking reactions from ', collFileName
                open(20,file='../../CollisionData/'//collFileName, action = 'read')
                read(20,*) string
                do while (.not. (string(1:3) == 'END' .or. string(1:3) == 'end'))
                    if( string(1:4).eq.'REAC' .or. string(1:4).eq.'Reac' .or. string(1:4).eq.'reac' ) then
                        ! add reaction
                        numberCollisions = numberCollisions + 1
                        numberReactants(numberCollisions) = 0
                        numberProducts(numberCollisions) = 0
                        read(20,*) string
                        if (string(1:1) .ne. '[') then
                            print *, 'Do not have chemical reaction after REACTION label'
                            stop
                        end if
                        lowerIndx = 0
                        higherIndx = 0
                        reactantBool = .true.
                        ! Get reactants and byproducts
                        do i = 1, 100
                            if (string(i:i) == '[') lowerIndx = i
                            if (string(i:i) == ']') higherIndx = i
                            if (string(i:i) == '>') then
                                reactantBool = .false.
                            end if
                            
                            ! Found particle name in reaction
                            if (higherIndx > lowerIndx) then
                                foundParticleBool = .false.
                                lowerIndx = lowerIndx + 1
                                higherIndx = higherIndx-1
                                ! Find particle in charged or neutral particles
                                do j = 1, numberChargedParticles
                                    if (particleList(j)%name == string(lowerIndx:higherIndx)) then
                                        foundParticleBool = .true.
                                        incidentParticleBool = .true.
                                        exit
                                    end if
                                end do
                                do k = 1, numberNeutralParticles
                                    if (targetParticleList(k)%name == string(lowerIndx:higherIndx)) then
                                        foundParticleBool = .true.
                                        incidentParticleBool = .false.
                                        exit
                                    end if
                                end do
                                ! for each found particle add to product or reactant count
                                if (foundParticleBool) then
                                    if (reactantBool) then
                                        numberReactants(numberCollisions) = numberReactants(numberCollisions) + 1
                                        if (incidentParticleBool) then
                                            collisionReactantsIndx(numberReactants(numberCollisions), numberCollisions) = j
                                        else
                                            collisionReactantsIndx(numberReactants(numberCollisions), numberCollisions) = k
                                        end if
                                    else
                                        numberProducts(numberCollisions) = numberProducts(numberCollisions) + 1
                                        if (incidentParticleBool) then
                                            collisionProductsIndx(numberProducts(numberCollisions), numberCollisions) = j
                                        else
                                            collisionProductsIndx(numberProducts(numberCollisions), numberCollisions) = k
                                        end if
                                    end if
                                else
                                    print *, 'Could not find particle for collision:', string(lowerIndx:higherIndx)
                                    stop
                                end if
                                lowerIndx = 0
                                higherIndx = 0
                            end if
                        end do
                       
                        ! If reactant pair belongs to set already brought up, add that to number of times set pops up, otherwise add another set for reactionIndx
                        oldSetBool = .false.
                        do i = 1, numberBinaryCollisions
                            oldSetBool = (collisionReactantsIndx(1, numberCollisions) == collisionReactantsIndx(1, i) .and. collisionReactantsIndx(2, numberCollisions) == collisionReactantsIndx(2, i))
                            if (oldSetBool) then
                                reactionIndx(numberCollisions) = i
                                exit
                            end if
                        end do
                        if (.not. oldSetBool) then
                            numberBinaryCollisions = numberBinaryCollisions + 1
                            collisionReactantsIndx(:, numberBinaryCollisions) = collisionReactantsIndx(:, numberCollisions)
                            reactionIndx(numberCollisions) = numberBinaryCollisions
                        end if

                        ! read reaction
                        read(20,*) E_threshold(numberCollisions), tempVar2
                        read(20,*) E_scaling, sigma_scaling
                        read(20,*) string
                        if( trim(string).eq.'ELASTIC' .or. trim(string).eq.'elastic' ) &
                            collisionType(numberCollisions)=1
                        if( trim(string).eq.'IONIZATION' .or. trim(string).eq.'ionization' ) &
                            collisionType(numberCollisions)=2
                        if( trim(string).eq.'EXCITATION' .or. trim(string).eq.'excitation' ) &
                            collisionType(numberCollisions)=3
                        if( trim(string).eq.'CHARGEEXCHANGE' .or. trim(string).eq.'chargeexchange' ) &
                            collisionType(numberCollisions)=4
                        if( trim(string).eq.'DISSOCIATION' .or. trim(string).eq.'dissociation' ) &
                            collisionType(numberCollisions)=5
                      
                        do while(string(1:4).ne.'----')
                            read(20,*) string
                        end do
                        numberSigmaPoints(numberCollisions) = 0
                        do i = 1, maxNumberPoints
                            read(20, '(4A)', ADVANCE = 'NO') string(1:4)
                            if (string(1:4) == '----') then
                                read(20, *) string
                                exit
                            else
                                backspace(20)
                                numberSigmaPoints(numberCollisions) = numberSigmaPoints(numberCollisions) + 1
                            end if
                            read(20, *) tempVar1, tempVar2
                            E_temp(numberSigmaPoints(numberCollisions), numberCollisions) = tempVar1 * E_scaling
                            sigma_temp(numberSigmaPoints(numberCollisions), numberCollisions) = tempVar2 * sigma_scaling
                        end do
                      
                    end if
                    read(20,*) string
                end do
                close(20)
            end if
            read(10,*) string
        end do
        close(10)

        if (numberBinaryCollisions > 0) then
            ! Binary collision for each primary on target set
            allocate(nullCollisionList(numberBinaryCollisions))
            do j = 1, numberBinaryCollisions
                ! get properties for each binary collision
                red_mass = particleList(collisionReactantsIndx(1, j))%mass * targetParticleList(collisionReactantsIndx(2, j))%mass / (particleList(collisionReactantsIndx(1, j))%mass + targetParticleList(collisionReactantsIndx(2, j))%mass)
                sumMass = (particleList(collisionReactantsIndx(1, j))%mass + targetParticleList(collisionReactantsIndx(2, j))%mass)
                redMassTripleProducts = 0
                numberCollisionsPerReaction = 0
                length_concatenate = 0
                ! accumulate total number of array points for energy/sigma
                do i = 1, numberCollisions
                    if (reactionIndx(i) == j) then
                        length_concatenate = length_concatenate + numberSigmaPoints(i)
                        numberCollisionsPerReaction = numberCollisionsPerReaction + 1
                    end if
                end do
                allocate(E_concatenate_temp(length_concatenate), indxSort(length_concatenate))
                k = 0
                do i = 1, numberCollisions
                    if (reactionIndx(i)==j) then
                        E_concatenate_temp(k+1:k+numberSigmaPoints(i)) = E_temp(1:numberSigmaPoints(i), i)
                        k = k + numberSigmaPoints(i)
                    end if
                end do

                ! Sort new concatenated array
                call indexSortArray(length_concatenate,E_concatenate_temp,indxSort)
                ! count amount of repeats to k
                k = 0
                do i = 1, length_concatenate-1
                    if (E_concatenate_temp(indxSort(i)) == E_concatenate_temp(indxSort(i+1))) then
                        k = k + 1
                    end if
                end do
                ! allocate proper size E_array
                allocate(E_concatenate(length_concatenate-k))
                ! place non-repeats in new array
                k = 0
                tempVar1 = 1.d10
                do i = 1, length_concatenate
                    if (E_concatenate_temp(indxSort(i)) /= tempVar1) then
                        k = k + 1
                        E_concatenate(k) = E_concatenate_temp(indxSort(i))
                        tempVar1 = E_concatenate(k)
                    end if
                end do
                length_concatenate = k
                ! Now interpolate each collision in set to new E_array
                allocate(sigma_concatenate(length_concatenate, numberCollisionsPerReaction), collisionTypeArray(numberCollisionsPerReaction), E_threshold_array(numberCollisionsPerReaction), &
                    numberProductsArray(numberCollisionsPerReaction), productsIndxArray(3, numberCollisionsPerReaction))
                h = 0
                do i = 1, numberCollisions
                    if (reactionIndx(i)==j) then
                        h = h + 1
                        collisionTypeArray(h) = collisionType(i)
                        E_threshold_array(h) = E_threshold(i)
                        numberProductsArray(h) = numberProducts(i)
                        productsIndxArray(:, h) = collisionProductsIndx(:, i)
                        if (collisionType(i) == 2) redMassTripleProducts = 1.0d0 / (1.0d0/particleList(collisionProductsIndx(1,i))%mass + 1.0d0/particleList(collisionProductsIndx(2,i))%mass + 1.0d0/particleList(collisionProductsIndx(3,i))%mass)
                        lowerIndx = 1
                        do k = 1, length_concatenate
                            if (E_threshold(i) > E_concatenate(k)) then
                                ! If Energy below threshold, then sigma is 0
                                sigma_concatenate(k, h) = 0.0d0
                            else if (E_concatenate(k) < E_temp(numberSigmaPoints(i), i)) then
                                ! Go up index in OG energy array until get to energy greater than current energy on concatenated array
                                do while (E_temp(lowerIndx, i) <= E_concatenate(k))
                                    lowerIndx = lowerIndx + 1
                                end do
                                lowerIndx = lowerIndx - 1
                                if (E_temp(lowerIndx, i) == E_concatenate(k)) then
                                    sigma_concatenate(k, h) = sigma_temp(lowerIndx, i)
                                else
                                    ! interpolate
                                    tempVar1 = E_concatenate(k) - E_temp(lowerIndx, i)
                                    tempVar2 = tempVar1/ (E_temp(lowerIndx+1,i) - E_temp(lowerIndx, i))
                                    sigma_concatenate(k, h) = (sigma_temp(lowerIndx, i) * (1 - tempVar2) + sigma_temp(lowerIndx+1, i) * (tempVar2))
                                end if
                            else
                                ! Energy outside range, use constant extrapolation to outer array
                                sigma_concatenate(k, h) = sigma_temp(numberSigmaPoints(i), i)
                            end if
                        end do
                    end if   
                end do
                ! create actual object
                nullCollisionList(j) = nullCollision(2, h, length_concatenate, red_mass, sumMass, redMassTripleProducts, E_concatenate, sigma_concatenate, E_threshold_array, collisionTypeArray, &
                    collisionReactantsIndx(:,j), numberProductsArray, productsIndxArray)
                deallocate(sigma_concatenate, collisionTypeArray, E_threshold_array, numberProductsArray, productsIndxArray)
                deallocate(E_concatenate)
                deallocate(E_concatenate_temp, indxSort)
                print *, ''
                print *, 'Null Collision generated'
                print *, 'Reactants are:', particleList(nullCollisionList(j)%reactantsIndx(1))%name, ' and ', targetParticleList(nullCollisionList(j)%reactantsIndx(2))%name
                print *, 'Amount of collisions:', nullCollisionList(j)%numberCollisions
                print *, 'Reduced mass:', nullCollisionList(j)%reducedMass
                print *, 'length of arrays:', nullCollisionList(j)%lengthArrays
                print *, 'Max sigma_v:', nullCollisionList(j)%sigmaVMax
                print *, 'reducedMass:', nullCollisionList(j)%reducedMass
                print *, 'reducedMassIonization:', nullCollisionList(j)%reducedMassIonization
                do i = 1, nullCollisionList(j)%numberCollisions
                    print *, ''
                    print *, 'For collision #:', i
                    print *, 'Energy threshold:', nullCollisionList(j)%energyThreshold(i)
                    print *, 'Collision type:', nullCollisionList(j)%collisionType(i)
                    print *, 'number Products are:', nullCollisionList(j)%numberProducts(i)
                    print *, 'products are:'
                    if (nullCollisionList(j)%collisionType(i) == 2) then
                        print *, particleList(nullCollisionList(j)%productsIndx(1, i))%name, ' and ', particleList(nullCollisionList(j)%productsIndx(2,i))%name, 'and ', particleList(nullCollisionList(j)%productsIndx(3,i))%name
                    else
                        print *, particleList(nullCollisionList(j)%productsIndx(1,i))%name, ' and ', targetParticleList(nullCollisionList(j)%productsIndx(2,i))%name
                    end if 
                end do
                print *, ''
            end do
        else
            print *, 'No binary collisions read!'
            allocate(nullCollisionList(0))
        end if

        print *, '---------------------------'

    end subroutine readNullCollisionInputs


end module mod_NullCollision