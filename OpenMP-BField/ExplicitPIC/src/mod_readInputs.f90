module mod_readInputs
    use iso_fortran_env
    use constants
    use mod_BasicFunctions
    use mod_domain
    use mod_particle
    use mod_targetParticle
    use mod_NullCollision
    use mod_potentialSolver
    use mod_particleInjection
    use omp_lib
    implicit none



contains

    ! subroutine readRestart(InitFilename, world, solver, particleList)
    !     type(Domain), intent(in out) :: world
    !     type(potentialSolver), intent(in out) :: solver
    !     type(Particle), allocatable, intent(out) :: particleList(:)
    !     character(len=*), intent(in) :: InitFilename
    !     integer(int32), allocatable :: workInt(:)
    !     real(real64), allocatable :: workReal(:)
    !     character(len=:), allocatable :: filename
    !     character(len=5) :: char_i
    !     character(len=10) :: name
    !     character(len=100) :: buf
    !     real(real64) :: mass, charge, w_p, temp
    !     integer(int32) :: io, i, finalIdx, N_p, j
    !     logical :: bool
    !     print *, "So you've decided to create a simulation based on the ending of another simulation."
    !     print *, "Please type in the name of the file to extract data from."
    !     print *, "If you're having second thoughts, put in a file that doesn't exist, and I'll quit on you."
    !     read *, buf
    !     filename = trim(buf)
    !     INQUIRE(DIRECTORY = '../'//filename, EXIST = bool)
    !     if (bool) then
    !         open(10,file='../'//filename//'/InitialConditions.dat')
    !         read(10, *)
    !         read(10, *) NumberXNodes, simulationTime, del_t, FractionFreq, world%delX, n_ave, T_e, T_i, numDiagnosticSteps, numberChargedParticles
    !         close(10)
    !         allocate(workInt(NumberXNodes), workReal(NumberXNodes))
    !         open(10, file = '../'//filename//'/domainBoundaryConditions.dat', form = 'unformatted')
    !         read(10) workInt
    !         close(10)
    !         world = Domain(workInt(1), workInt(NumberXNodes))
    !         open(10, file = '../'//filename//'/domainGrid.dat', form = 'unformatted')
    !         read(10) workReal
    !         close(10)
    !         world%grid = workReal
    !         world%delX = world%grid(2) - world%grid(1)

    !         write(char_i, '(I3)') numDiagnosticSteps
    !         open(10, file = '../'//filename//'/Phi/phi_'//trim(adjustl(char_i))//'.dat', form = 'unformatted')
    !         read(10) workReal
    !         close(10)
    !         solver = potentialSolver(world, workReal(1), workReal(NumberXNodes))
    !         call solver%construct_diagMatrix(world)
    !         solver%phi = workReal
    !         allocate(particleList(numberChargedParticles))
    !         print *, 'size particleList is:', numberChargedParticles
    !         open(10,file='../'//filename//'/ParticleProperties.dat')
    !         read(10, *)
    !         do i = 1, numberChargedParticles
    !             read(10, *) name, mass, charge, w_p, finalIdx
    !             open(13, file='../'//filename//'/ParticleDiagnostic_'//trim(name)//'.dat')
    !             do j = 1, numDiagnosticSteps
    !                 read(13, *)
    !             end do
    !             read(13, *) temp, temp, temp, temp, temp, N_p
    !             particleList(i) = Particle(mass, charge, w_p, N_p, finalIdx, trim(name))
    !             open(15, file='../'//filename//'/PhaseSpace/phaseSpace_'//trim(name)//'_'//trim(adjustl(char_i))//'.dat', form = 'unformatted')
    !             read(15) particleList(i)%phaseSpace(:, 1:particleList(i)%N_p)
    !             close(15)
    !         end do
    !         close(13)
    !         close(10)
    !         open(10,file='../InputData/'//InitFilename, IOSTAT=io)
    !         read(10, *, IOSTAT = io) simulationTime
    !         read(10, *, IOSTAT = io) 
    !         read(10, *, IOSTAT = io) 
    !         read(10, *, IOSTAT = io) 
    !         read(10, *, IOSTAT = io) numDiagnosticSteps
    !         read(10, *, IOSTAT = io) fractionFreq
    !         read(10, *, IOSTAT = io) averagingTime
    !         read(10, *, IOSTAT = io)
    !         read(10, *, IOSTAT = io) 
    !         read(10, *, IOSTAT = io) 
    !         read(10, *, IOSTAT = io) buf
    !         close(10)
    !         directoryName = trim(buf)
    !         if (len(directoryName) < 2) then
    !             stop "Directory name length less than 2 characters, be more descriptive!"
    !         end if
    !         print *, ""
    !         print *, "--------------------"
    !         print *, "Save data folder: ", directoryName
    !         print *, "Average initial electron density:", n_ave
    !         print *, "Initial electron temperature:", T_e
    !         print *, "Initial ion temperature:", T_i
    !         print *, "Number of diagnostic steps is:", numDiagnosticSteps
    !         print *, "Fraction of 1/w_p for time step:", fractionFreq
    !         print *, "Final averaging time is:", averagingTime
    !         print *, "------------------"
    !         print *, ""
    !         print *, "------------------"
    !         print *, "Number of nodes:", NumberXNodes
    !         print *, "Grid length:", world%grid(NumberXNodes) - world%grid(1)
    !         print *, "Left boundary type:", world%boundaryConditions(1)
    !         print *, "Right boundary type:", world%boundaryConditions(NumberXNodes)
    !         print *, "------------------"
    !         print *, ""
    !         do j=1, numberChargedParticles
    !             print *, 'Particle name ', particleList(j) % name
    !             print *, 'Amount of macroparticles is:', particleList(j) % N_p
    !             print *, "Particle mass is:", particleList(j)%mass
    !             print *, "Particle charge is:", particleList(j)%q
    !             print *, "Particle mean KE is:", particleList(j)%getKEAve()
    !             print *, ""
    !         end do
    !     else
    !         stop "File you're trying to restart from doesn't exist."
    !     end if
    ! end subroutine readRestart

    subroutine readInjectionInputs(InjFilename, addLostPartBool, refluxPartBool, injectionBool, uniformInjectionBool, heatingBool, injectionFlux, w_p, angleRad, FractionFreqHeating)
        logical, intent(in out) :: addLostPartBool, refluxPartBool, injectionBool, uniformInjectionBool, heatingBool
        real(real64), intent(in out) :: injectionFlux, FractionFreqHeating
        real(real64), intent(in) :: w_p, angleRad
        character(len=*), intent(in) :: InjFilename
        integer(int32) :: tempInt, io
        real(real64) :: numFluxPart
        print *, ""
        print *, "Reading initial inputs for particle injection:"
        print *, "------------------"
        open(10,file='../InputData/'//InjFilename, IOSTAT=io)
        read(10, *, IOSTAT = io) tempInt
        addLostPartBool = (tempInt == 1)
        read(10, *, IOSTAT = io) tempInt
        refluxPartBool = (tempInt == 1)
        read(10, *, IOSTAT = io) tempInt, injectionFlux
        injectionBool = (tempInt == 1)
        if (.not. injectionBool) then
            read(10, *, IOSTAT = io) tempInt, injectionFlux
            uniformInjectionBool = (tempInt == 1)
        end if
        read(10, *, IOSTAT = io) tempInt, FractionFreqHeating
        heatingBool = (tempInt == 1)
        close(10)
        print *, "Particle lost is reinjected:", addLostPartBool
        print *, "Particle refluxing activated on neumann boundary:", refluxPartBool
        print *, "Particle injection on neumann boundary", injectionBool
        print *, "Particle injection unformly with maxwellian", uniformInjectionBool
        print *, "Electron maxwellian heating", heatingBool
        if (injectionBool) then
            injectionFlux = injectionFlux * COS(angleRad)
            print *, 'Particle injection flux:', injectionFlux
            numFluxPart = injectionFlux * del_t / w_p/real(numThread) ! particles injected per thread
            numFluxParticlesLow = floor(numFluxPart)
            numFluxParticlesHigh = numFluxParticlesLow + 1
            print *, 'Low end of flux particles:', numFluxParticlesLow
            print *, 'High end of flux particles:', numFluxParticlesHigh
            injectionR = numFluxPart - real(numFluxParticlesLow)
            print *, 'Number for selection of flux particles is:', injectionR
        else if (uniformInjectionBool) then
            print *, 'Particle injection flux:', injectionFlux
            numFluxPart = injectionFlux * del_t / w_p/real(numThread) ! particles injected per thread
            numFluxParticlesLow = floor(numFluxPart)
            numFluxParticlesHigh = numFluxParticlesLow + 1
            print *, 'Low end of flux particles:', numFluxParticlesLow
            print *, 'High end of flux particles:', numFluxParticlesHigh
            injectionR = numFluxPart - real(numFluxParticlesLow)
            print *, 'Number for selection of flux particles is:', injectionR
        else if (heatingBool) then
            print *, 'FractionFreqHeating for maxwellian heating is:', FractionFreqHeating
        end if
        print *, "------------------"
        print *, ""
    end subroutine readInjectionInputs

    subroutine readInitialInputs(InitFilename, simulationTime, n_ave, T_e, T_i, numDiagnosticSteps, fractionFreq, averagingTime, numThread, irand, stateRanNew)
        real(real64), intent(in out) :: fractionFreq, n_ave, simulationTime, averagingTime, T_e, T_i
        integer(int32), intent(in out) :: numDiagnosticSteps, numThread
        integer(int32), allocatable, intent(out) :: irand(:)
        integer(int32), allocatable, intent(out) :: stateRanNew(:,:)
        character(len=*), intent(in) :: InitFilename
        integer(int32) :: io, i
        character(len=100) :: tempName
        real(real64) :: rando
        print *, ""
        print *, "Reading initial inputs:"
        print *, "------------------"
        open(10,file='../InputData/'//InitFilename, IOSTAT=io)
        read(10, *, IOSTAT = io) numThread
        read(10, *, IOSTAT = io) simulationTime
        read(10, *, IOSTAT = io) n_ave
        read(10, *, IOSTAT = io) T_e
        read(10, *, IOSTAT = io) T_i
        read(10, *, IOSTAT = io) numDiagnosticSteps
        read(10, *, IOSTAT = io) fractionFreq, del_t
        read(10, *, IOSTAT = io) averagingTime
        read(10, *, IOSTAT = io) 
        read(10, *, IOSTAT = io) 
        read(10, *, IOSTAT = io) 
        read(10, *, IOSTAT = io) tempName
        close(10)
        directoryName = trim(tempName)
        if (len(directoryName) < 2) then
            stop "Directory name length less than 2 characters!"
        end if
        directoryName = '../../../../ExplicitData-BField/'//directoryName
        if (numThread > omp_get_num_procs()) then
            print *, "Number of threads set is larger than the maximum number of threads which is", omp_get_num_procs()
            stop
        end if
        call omp_set_num_threads(numThread)
        allocate(irand(numThread), stateRanNew(2,numThread))
        do i = 1, numThread
            call random_number(rando)
            irand(i) = INT(rando * (huge(i)))
            call random_number(rando)
            stateRanNew(1, i) = INT(2.0d0 * (rando-0.5d0) * (huge(stateRanNew(1, i))))
            stateRanNew(2, i) = INT(rando * (huge(stateRanNew(2, i))))
        end do
        del_t = MIN(fractionFreq * 1.0d0 / getPlasmaFreq(n_ave), del_t)
        print *, "Save data folder: ", directoryName
        print *, "Number of threads is:", numThread
        print *, "Average initial electron density:", n_ave
        print *, "Initial electron temperature:", T_e
        print *, "Initial ion temperature:", T_i
        print *, "Number of diagnostic steps is:", numDiagnosticSteps
        print *, "Fraction of 1/w_p for time step:", fractionFreq
        print *, 'del_t will be:', del_t
        print *, "Simulation time is:", simulationTime
        print *, "Final averaging time is:", averagingTime
        print *, "------------------"
        print *, ""
    end subroutine readInitialInputs

    subroutine readGeometry(world, solver, GeomFilename)
        type(Domain), intent(in out) :: world
        type(potentialSolver), intent(in out) :: solver
        character(len=*), intent(in) :: GeomFilename
        integer(int32) :: io, leftBoundary, rightBoundary
        real(real64) :: leftVoltage, rightVoltage, debyeLength, L_domain, fracDebye, BFieldMag, angle, RF_frequency

        print *, ""
        print *, "Reading domain inputs:"
        print *, "------------------"
        open(10,file='../InputData/'//GeomFilename)
        read(10, *, IOSTAT = io) NumberXNodes
        read(10, *, IOSTAT = io) L_domain
        read(10, *, IOSTAT = io) leftBoundary, rightBoundary
        read(10, *, IOSTAT = io) leftVoltage, rightVoltage
        read(10, *, IOSTAT = io) fracDebye
        read(10, *, IOSTAT = io) BFieldMag, angle
        read(10, *, IOSTAT = io) RF_frequency
        close(10)
        debyeLength = getDebyeLength(T_e, n_ave)
        if (L_domain / (NumberXNodes-1) > debyeLength) then
            print *, "Insufficient amount of nodes to resolve initial debyeLength"
            NumberXNodes = NINT(L_domain/debyeLength/fracDebye) + 1
        end if
        print *, "Number of nodes:", NumberXNodes
        print *, "Grid length:", L_domain
        print *, "Left boundary type:", leftBoundary
        print *, "Right boundary type:", rightBoundary
        print *, 'Fraction of debye length:', fracDebye
        print *, 'BField magnitude:', BFieldMag
        print *, 'BField angle:', angle
        print *, 'RF frequency:', RF_frequency
        if ((leftBoundary == 3) .or. (rightBoundary == 3)) then
            leftBoundary = 3
            rightBoundary = 3
            leftVoltage = rightVoltage
        end if
        world = Domain(leftBoundary, rightBoundary)
        call world % constructGrid(L_domain)
        solver = potentialSolver(world, leftVoltage, rightVoltage, BFieldMag, angle, RF_frequency)
        print *, 'delX is:', world%delX
        print *, "BField vector:", solver%BField
        print *, "------------------"
        print *, ""
        call solver%construct_diagMatrix(world)

    end subroutine readGeometry

    subroutine readParticleInputs(filename, numberChargedParticles, irand, T_e, T_i, numThread, world, particleList, targetParticleList)
        type(Particle), allocatable, intent(out) :: particleList(:)
        type(targetParticle), allocatable, intent(out) :: targetParticleList(:)
        type(Domain) :: world
        character(len=*), intent(in) :: filename
        integer(int32), intent(in) :: numThread
        integer(int32), intent(in out) :: numberChargedParticles
        integer(int32), intent(in out) :: irand(2,numThread)
        real(real64), intent(in) :: T_e, T_i
        integer(int32) :: j, numSpecies = 0, numNeutral = 0, numParticles(100), particleIdxFactor(100)
        character(len=15) :: name
        character(len=8) :: particleNames(100), neutralName
        real(real64) :: mass(100), charge(100), Ti(100), mass_neutral, Temp_neutral, density_neutral

        print *, "Reading particle inputs:"
        open(10,file='../InputData/'//filename, action = 'read')
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
                call particleList(j) % initialize_n_ave(n_ave, (world%grid(NumberXNodes) - world%grid(1)))
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
        real(real64) :: tempVar1, tempVar2, E_scaling, sigma_scaling, E_threshold(maxNumberColls), E_temp(maxNumberPoints, maxNumberColls), sigma_temp(maxNumberPoints, maxNumberColls), red_mass, sumMass, redMassTripleProducts
        integer(int32) :: i, j, k, numberCollisions, reactionIndx(maxNumberColls), lowerIndx, higherIndx, numberReactants(maxNumberColls), numberProducts(maxNumberColls), collisionReactantsIndx(2, maxNumberColls), &
            collisionProductsIndx(3, maxNumberColls), numberSigmaPoints(maxNumberColls), length_concatenate, collisionType(maxNumberColls), h, numberCollisionsPerReaction
        real(real64), allocatable :: E_concatenate_temp(:), E_concatenate(:), sigma_concatenate(:, :), E_threshold_array(:)
        integer(int32), allocatable :: indxSort(:), collisionTypeArray(:), numberProductsArray(:), productsIndxArray(:,:)
        

        print *, "Reading collision inputs:"
        print *, '---------------------------'
        numberBinaryCollisions = 0
        numberCollisions = 0
        open(10,file='../InputData/'//filename, action = 'read')
        read(10,*) string
        do while (.not. (string(1:3) == 'END' .or. string(1:3) == 'end'))
            collFileName = trim(string)
            inquire(file = '../../CollisionData/'//collFileName, exist = fileExists)
            if (fileExists) then
                print *, 'Taking reactions from ', collFileName
                open(20,file='../../CollisionData/'//collFileName, action = 'read')
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
        end if

        print *, '---------------------------'

    end subroutine readNullCollisionInputs

end module mod_readInputs