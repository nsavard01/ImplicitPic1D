module mod_readInputs
    use iso_fortran_env, only: int32, real64
    use constants
    use mod_BasicFunctions
    use mod_particle
    use mod_targetParticle
    use mod_NullCollision
    use omp_lib
    implicit none



contains

    subroutine readParticleInputs(filename, numberChargedParticles, irand, T_e, T_i, numThread, particleList, targetParticleList)
        type(Particle), allocatable, intent(out) :: particleList(:)
        type(targetParticle), allocatable, intent(out) :: targetParticleList(:)
        character(len=*), intent(in) :: filename
        integer(int32), intent(in) :: numThread
        integer(int32), intent(in out) :: numberChargedParticles, irand(numThread)
        real(real64), intent(in) :: T_e, T_i
        integer(int32) :: j, numSpecies = 0, numNeutral = 0, numParticles(100), particleIdxFactor(100)
        character(len=15) :: name
        character(len=8) :: particleNames(100), neutralName
        real(real64) :: mass(100), charge(100), Ti(100), mass_neutral, Temp_neutral, density_neutral

        print *, "Reading particle inputs:"
        open(10,file=filename, action = 'read')
        do j=1, 10000
            read(10,*) name

            if( name(1:9).eq.'ELECTRONS') then
                read(10,*) name
                read(10,*) name
                read(10,'(A4)', ADVANCE = 'NO') name(1:2)
                numSpecies = numSpecies + 1
                read(10,*) numParticles(numSpecies), particleIdxFactor(numSpecies)
                Ti(numSpecies) = T_e
                mass(numSpecies) = m_e
                charge(numSpecies) = -1.0
                particleNames(numSpecies) = 'e'
                read(10,*) name
                read(10,*) name
            endif


            if(name(1:4).eq.'IONS' .or. name(1:4).eq.'Ions' .or. name(1:4).eq.'ions' ) then
                do while(name(1:4).ne.'----')
                    read(10,*) name
                end do
                read(10,'(A6)', ADVANCE = 'NO') name
                do while (name(1:4).ne.'----')
                    numSpecies = numSpecies + 1
                    read(10,*) mass(numSpecies),charge(numSpecies), numParticles(numSpecies), particleIdxFactor(numSpecies)
                    Ti(numSpecies) = T_i
                    mass(numSpecies) = mass(numSpecies) * m_amu - charge(numSpecies) * m_e
                    particleNames(numSpecies) = trim(name)
                    read(10,'(A6)', ADVANCE = 'NO') name
                end do
            endif

            if( name(1:4).eq.'NEUT' .or. name(1:4).eq.'Neut' .or. name(1:4).eq.'neut' ) then
                do while(name(1:4).ne.'----')
                    read(10,*) name
                end do
                read(10,'(A6)', ADVANCE = 'NO') name
                do while (name(1:4).ne.'----')
                    numNeutral = numNeutral + 1
                    read(10,*) mass_neutral, Temp_neutral, density_neutral
                    mass_neutral = mass_neutral * m_amu
                    neutralName = trim(name)
                    read(10,'(A6)', ADVANCE = 'NO') name
                end do
            endif       

            if (name(1:7) == 'ENDFILE') then
                close(10)
                exit
            end if

        end do

        numberChargedParticles = numSpecies
        if (numberChargedParticles > 0) then
            allocate(particleList(numberChargedParticles))
            do j=1, numberChargedParticles
                particleList(j) = Particle(mass(j), e * charge(j), 1.0d0, numParticles(j), numParticles(j) * particleIdxFactor(j), trim(particleNames(j)), numThread)
                call particleList(j) % generate3DMaxwellian(Ti(j), irand)
                call particleList(j)% initializeRandUniform(irand)
                print *, 'Initializing ', particleList(j) % name
                print *, 'Amount of macroparticles is:', SUM(particleList(j) % N_p)
                print *, "Particle mass is:", particleList(j)%mass
                print *, "Particle charge is:", particleList(j)%q
                print *, "Particle weight is:", particleList(j)%w_p
                print *, "Particle mean KE is:", particleList(j)%getKEAve(), ", should be", Ti(j) * 1.5
            end do
        end if

        numberNeutralParticles = numNeutral
        if (numberNeutralParticles > 0) then
            allocate(targetParticleList(numberNeutralParticles))
            do j = 1, numberNeutralParticles
                targetParticleList(j) = targetParticle(neutralName, mass_neutral, density_neutral, Temp_neutral)
                print *, 'Initializing target particle:', targetParticleList(j)%name
                print *, 'Particle mass is:', targetParticleList(j)%mass
                print *, 'Particle temperature(K) is:', targetParticleList(j)%temperature
                print *, 'Particle density is:', targetParticleList(j)%density
            end do
        end if

        
        print *, "---------------"
        print *, ""



    end subroutine readParticleInputs


    subroutine readNullCollisionInputs(filename, nullCollisionList, particleList, targetParticleList, numberBinaryCollisions)
        ! Reaction is any unique set of particle on target. For each reaction, you may have several different collision types
        character(len=*), intent(in) :: filename
        type(nullCollision), allocatable, intent(out) :: nullCollisionList(:)
        type(Particle), intent(in) :: particleList(numberChargedParticles)
        type(targetParticle), intent(in) :: targetParticleList(numberNeutralParticles)
        integer(int32), intent(in out) :: numberBinaryCollisions
        integer(int32), parameter :: maxNumberPoints = 1500, maxNumberColls = 65
        logical :: fileExists, foundParticleBool, reactantBool, incidentParticleBool, oldSetBool
        character(len=100) :: string
        character(:), allocatable :: collFileName
        real(real64) :: tempVar1, tempVar2, E_scaling, sigma_scaling, E_threshold(maxNumberColls), E_temp(maxNumberPoints, maxNumberColls), sigma_temp(maxNumberPoints, maxNumberColls), v_r, red_mass, sumMass, redMassTripleProducts
        integer(int32) :: i, j, k, numberCollisions, reactionIndx(maxNumberColls), lowerIndx, higherIndx, numberReactants(maxNumberColls), numberProducts(maxNumberColls), collisionReactantsIndx(2, maxNumberColls), &
            collisionProductsIndx(3, maxNumberColls), numberSigmaPoints(maxNumberColls), length_concatenate, collisionType(maxNumberColls), h, numberCollisionsPerReaction
        real(real64), allocatable :: E_concatenate_temp(:), E_concatenate(:), sigma_v_concatenate(:, :), E_threshold_array(:)
        integer(int32), allocatable :: indxSort(:), collisionTypeArray(:), numberProductsArray(:), productsIndxArray(:,:)
        

        print *, "Reading collision inputs:"
        print *, '---------------------------'
        numberBinaryCollisions = 0
        numberCollisions = 0
        open(10,file=filename, action = 'read')
        read(10,*) string
        do while (.not. (string(1:3) == 'END' .or. string(1:3) == 'end'))
            collFileName = trim(string)
            inquire(file = '../CollisionData/'//collFileName, exist = fileExists)
            if (fileExists) then
                print *, 'Taking reactions from ', collFileName
                open(20,file='../CollisionData/'//collFileName, action = 'read')
                read(20,*) string
                do while (.not. (string(1:3) == 'END' .or. string(1:3) == 'end'))
                    if( string(1:4).eq.'REAC' .or. string(1:4).eq.'Reac' .or. string(1:4).eq.'reac' ) then
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
                
                            if (higherIndx > lowerIndx) then
                                foundParticleBool = .false.
                                lowerIndx = lowerIndx + 1
                                higherIndx = higherIndx-1
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
            else
                print *, 'Collision file:', collFileName, 'does not exist!'
                stop
            end if
            read(10,*) string
        end do
        close(10)

        if (numberBinaryCollisions > 0) then
            allocate(nullCollisionList(numberBinaryCollisions))
            do j = 1, numberBinaryCollisions
                red_mass = particleList(collisionReactantsIndx(1, j))%mass * targetParticleList(collisionReactantsIndx(2, j))%mass / (particleList(collisionReactantsIndx(1, j))%mass + targetParticleList(collisionReactantsIndx(2, j))%mass)
                sumMass = (particleList(collisionReactantsIndx(1, j))%mass + targetParticleList(collisionReactantsIndx(2, j))%mass)
                redMassTripleProducts = 0
                numberCollisionsPerReaction = 0
                length_concatenate = 0
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
                allocate(sigma_v_concatenate(length_concatenate, numberCollisionsPerReaction), collisionTypeArray(numberCollisionsPerReaction), E_threshold_array(numberCollisionsPerReaction), &
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
                            v_r = SQRT(2.0d0 * E_concatenate(k) * e / red_mass)
                            if (E_threshold(i) > E_concatenate(k)) then
                                ! If Energy below threshold, then sigma is 0
                                sigma_v_concatenate(k, h) = 0.0d0
                            else if (E_concatenate(k) < E_temp(numberSigmaPoints(i), i)) then
                                ! Go up index in OG energy array until get to energy greater than current energy on concatenated array
                                do while (E_temp(lowerIndx, i) <= E_concatenate(k))
                                    lowerIndx = lowerIndx + 1
                                end do
                                lowerIndx = lowerIndx - 1
                                if (E_temp(lowerIndx, i) == E_concatenate(k)) then
                                    sigma_v_concatenate(k, h) = sigma_temp(lowerIndx, i) * v_r
                                else
                                    ! interpolate
                                    tempVar1 = E_concatenate(k) - E_temp(lowerIndx, i)
                                    tempVar2 = tempVar1/ (E_temp(lowerIndx+1,i) - E_temp(lowerIndx, i))
                                    sigma_v_concatenate(k, h) = (sigma_temp(lowerIndx, i) * (1 - tempVar2) + sigma_temp(lowerIndx+1, i) * (tempVar2)) * v_r
                                end if
                            else
                                ! Energy outside range, use constant extrapolation to outer array
                                sigma_v_concatenate(k, h) = sigma_temp(numberSigmaPoints(i), i) * v_r
                            end if
                        end do
                    end if   
                end do
                nullCollisionList(j) = nullCollision(2, h, length_concatenate, red_mass, sumMass, redMassTripleProducts, E_concatenate, sigma_v_concatenate, E_threshold_array, collisionTypeArray, &
                    collisionReactantsIndx(:,j), numberProductsArray, productsIndxArray)
                deallocate(sigma_v_concatenate, collisionTypeArray, E_threshold_array, numberProductsArray, productsIndxArray)
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
                    print *, 'For collision #:', i
                    print *, 'Energy threshold:', nullCollisionList(j)%energyThreshold(i)
                    print *, 'Collision types:', nullCollisionList(j)%collisionType(i)
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
        end if

        print *, '---------------------------'

    end subroutine readNullCollisionInputs

end module mod_readInputs