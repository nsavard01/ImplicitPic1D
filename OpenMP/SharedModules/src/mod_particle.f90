module mod_particle

    use iso_fortran_env, only: int32, real64, output_unit
    use constants
    use mod_BasicFunctions
    use mod_domain
    use omp_lib
    implicit none

    private
    public :: Particle, readChargedParticleInputs
    integer(int32), public, protected :: numberChargedParticles = 0

    ! Particle contains particle properties and stored values in phase space df
    type :: Particle
        character(:), allocatable :: name !name of the particle
        integer(int32), allocatable :: N_p(:), refIdx(:), refRecordIdx(:, :), delIdx(:), numToCollide(:) !N_p is the current last index of particle, final idx is the last allowable in memory. Index starts at 1
        integer(int32) :: finalIdx, accumWallLoss(2)
        real(real64), allocatable :: phaseSpace(:,:, :) !particle phase space, represents [l_x, v_x, v_y, v_z] in first index
        real(real64) :: mass, q, w_p, q_over_m, q_times_wp ! mass (kg), charge(C), and weight (N/m^2 in 1D) of particles. Assume constant weight for moment
        real(real64), allocatable :: wallLoss(:,:), energyLoss(:,:), momentumLoss(:,:)
        real(real64), allocatable :: densities(:, :), workSpace(:,:)
        real(real64) :: numFuncEvalAve, numSubStepsAve
        real(real64) :: accumEnergyLoss(2)

    contains
        procedure, public, pass(self) :: initialize_n_ave
        procedure, public, pass(self) :: initializeRandUniform
        procedure, public, pass(self) :: initializeRandCosine
        procedure, public, pass(self) :: initializeRandSine
        procedure, public, pass(self) :: generate3DMaxwellian
        procedure, public, pass(self) :: getKEAve
        procedure, public, pass(self) :: getTotalKE
        procedure, public, pass(self) :: getTotalMomentum
        procedure, public, pass(self) :: writePhaseSpace
        procedure, public, pass(self) :: writeLocalTemperature
    end type Particle


    interface Particle
        module procedure :: particle_constructor
    end interface Particle

contains

    type(Particle) function particle_constructor(mass, q, w_p, N_p, finalIdx, particleName, numThread) result(self)
        ! Construct particle object, sizeIncrease is fraction larger stored array compared to initial amount of particles
        ! In future, use hash function for possible k = 1 .. Nx, m amount of boundaries, p = prime number  m < p < N_x. h(k) = (k%p)%m
        real(real64), intent(in) :: mass, q, w_p
        integer(int32), intent(in) :: N_p, finalIdx, numThread
        character(*), intent(in) :: particleName
        integer(int32) :: maxNumNodes
        self % name = particleName
        self % mass = mass
        self % q = q
        self % w_p = w_p
        self%q_over_m = q/mass
        self%q_times_wp = q * w_p
        self % finalIdx = finalIdx
        self%accumWallLoss = 0.0d0
        self%accumEnergyLoss = 0.0d0
        self%numFuncEvalAve = 0.0d0
        self%numSubStepsAve = 0.0d0
        maxNumNodes = MAX(NumberXNodes, NumberXHalfNodes)
        allocate(self%phaseSpace(4,finalIdx, numThread), self%refRecordIdx(INT(self%finalIdx/10), numThread), self%N_p(numThread), &
            self%delIdx(numThread), self%wallLoss(2, numThread), self%energyLoss(2, numThread), self%refIdx(numThread), &
            self%densities(NumberXNodes, numThread), self%workSpace(maxNumNodes, numThread), self%numToCollide(numThread), self%momentumLoss(2, numThread))
        self%refIdx = 0
        self%delIdx = 0
        self%N_p = N_p
        self%energyLoss = 0.0d0
        self%wallLoss = 0
        self%momentumLoss = 0.0d0
    end function particle_constructor

    pure subroutine initialize_n_ave(self, n_ave, L_domain)
        ! initialize w_p based on initial average n (N/m^3) over domain
        class(Particle), intent(in out) :: self
        real(real64), intent(in) :: n_ave, L_domain
        self % w_p = n_ave * L_domain / SUM(self % N_p)
        self%q_times_wp = self%q * self%w_p
    end subroutine initialize_n_ave

    subroutine initializeRandUniform(self, world, irand)
        class(Particle), intent(in out) :: self
        type(Domain), intent(in) :: world
        integer(int32), intent(in out) :: irand(numThread)
        integer(int32) :: iThread, i
        real(real64) :: x_pos
        !$OMP parallel private(iThread, i, x_pos)
        iThread = omp_get_thread_num() + 1
        do i = 1, self%N_p(iThread)
            x_pos = ran2(irand(iThread)) * world%L_domain + world%startX
            self%phaseSpace(1,i, iThread) = world%getLFromX(x_pos)
        end do
        !$OMP end parallel
    end subroutine initializeRandUniform

    subroutine initializeRandCosine(self, world, irand, alpha)
        class(Particle), intent(in out) :: self
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: alpha
        integer(int32), intent(in out) :: irand(numThread)
        integer(int32) :: iThread, i, randNum
        real(real64) :: Rand1, Rand2, x_pos, val
        !$OMP parallel private(iThread, i, randNum, Rand1, Rand2, x_pos, val)
        iThread = omp_get_thread_num() + 1
        randNum = irand(iThread)
        do i = 1, self%N_p(iThread)
            Rand1 = ran2(randNum)
            Rand2 = ran2(randNum)
            val = (1.0d0 + alpha * COS(2 * pi * Rand2))/(1.0d0 + alpha)
            do while (Rand1 > val)
                Rand1 = ran2(randNum)
                Rand2 = ran2(randNum)
                val = (1.0d0 + alpha * COS(2 * pi * Rand2))/(1.0d0 + alpha)
            end do
            x_pos = Rand2 * world%L_domain + world%startX
            self%phaseSpace(1,i, iThread) = world%getLFromX(x_pos)
        end do
        irand(iThread) = randNum
        !$OMP end parallel
    end subroutine initializeRandCosine

    subroutine initializeRandSine(self, world, irand, alpha)
        class(Particle), intent(in out) :: self
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: alpha
        integer(int32), intent(in out) :: irand(numThread)
        integer(int32) :: iThread, i, randNum
        real(real64) :: Rand1, Rand2, x_pos, val
        !$OMP parallel private(iThread, i, randNum, Rand1, Rand2, x_pos, val)
        iThread = omp_get_thread_num() + 1
        randNum = irand(iThread)
        do i = 1, self%N_p(iThread)
            Rand1 = ran2(randNum)
            Rand2 = ran2(randNum)
            val = (1.0d0 + alpha * SIN(2 * pi * Rand2))/(1.0d0 + alpha)
            do while (Rand1 > val)
                Rand1 = ran2(randNum)
                Rand2 = ran2(randNum)
                val = (1.0d0 + alpha * SIN(2 * pi * Rand2))/(1.0d0 + alpha)
            end do
            x_pos = Rand2 * world%L_domain + world%startX
            self%phaseSpace(1,i, iThread) = world%getLFromX(x_pos)
        end do
        irand(iThread) = randNum
        !$OMP end parallel
    end subroutine initializeRandSine

    ! ------------------------ Generating and initializing velocity with Maxwellian --------------------------


    subroutine generate3DMaxwellian(self, T, world, irand, alpha, distType, v_drift)
        ! random velocity generator for the particle for temperature T (eV)
        ! Use box-muller method for random guassian variable, same as gwenael but doesn't have factor 2? Maybe factored into v_th
        class(Particle), intent(in out) :: self
        class(Domain), intent(in) :: world
        real(real64), intent(in) :: T, alpha, v_drift
        integer(int32), intent(in) :: distType
        integer(int32), intent(in out) :: irand(numThread)
        integer(int32) :: iThread, i
        integer(int32) :: irand_thread
        real(real64) :: U1, U2, U3, U4, v_therm, v_extra, x_pos
        v_therm = SQRT(T*e/self%mass)
        !$OMP PARALLEL PRIVATE(iThread, i, U1, U2, U3, U4, irand_thread, v_extra, x_pos)
        iThread = omp_get_thread_num() + 1
        irand_thread = irand(iThread)
        do i = 1, self%N_p(iThread)
            U1 = ran2(irand_thread)
            U2 = ran2(irand_thread)
            U3 = ran2(irand_thread)
            U4 = ran2(irand_thread)
            SELECT CASE (distType)
            CASE(0)
                v_extra = v_drift
            CASE(1)
                x_pos = world%getXFromL(self%phaseSpace(1, i, iThread))/world%L_domain
                v_extra = v_drift * (1.0d0 + alpha * COS(2 * pi * x_pos))
            CASE(2)
                x_pos = world%getXFromL(self%phaseSpace(1, i, iThread))/world%L_domain
                v_extra = v_drift * (1.0d0 - alpha * SIN(2 * pi * x_pos)) 
            END SELECT
            self%phaseSpace(2, i, iThread) = v_therm * SQRT(-2.0d0 * LOG(U1)) * COS(2.0d0 * pi * U2) + v_extra
            self%phaseSpace(3, i, iThread) = v_therm * SQRT(-2.0d0 * LOG(U1)) * SIN(2.0d0 * pi * U2)
            self%phaseSpace(4, i, iThread) = v_therm * SQRT(-2.0d0 * LOG(U3)) * SIN(2.0d0 * pi * U4)
            
        end do
        irand(iThread) = irand_thread
        !$OMP END PARALLEL
    end subroutine generate3DMaxwellian

    function getKEAve(self) result(res)
        ! calculate average kinetic energy (temperature) in eV
        class(Particle), intent(in) :: self
        integer(int32) :: iThread
        real(real64) :: res
        res = 0.0d0
        !$OMP parallel private(iThread) reduction(+:res)
        iThread = omp_get_thread_num() + 1
        res = res + SUM(self%phaseSpace(2:4, 1:self%N_p(iThread), iThread)**2) 
        !$OMP end parallel
        res = res * self % mass * 0.5d0 / e / SUM(self%N_p)
    end function getKEAve

    function getTotalMomentum(self) result(res)
        class(Particle), intent(in) :: self
        real(real64) :: res(3), temp(3)
        integer(int32) :: iThread
        temp = 0
        !$OMP parallel private(iThread) reduction(+:temp)
        iThread = omp_get_thread_num() + 1
        temp = temp +  SUM(self%phaseSpace(2:4, 1:self%N_p(iThread), iThread), DIM = 2)
        !$OMP end parallel
        res = temp * self%w_p * self%mass
    end function getTotalMomentum

    function getTotalKE(self) result(res)
        ! calculate total KE in Joules/m^2
        class(Particle), intent(in) :: self
        real(real64) :: res
        integer(int32) :: iThread
        res = 0.0d0
        !$OMP parallel private(iThread) reduction(+:res)
        iThread = omp_get_thread_num() + 1
        res = res + SUM(self%phaseSpace(2:4, 1:self%N_p(iThread), iThread)**2) 
        !$OMP end parallel
        res = res * self % mass * 0.5d0 * self%w_p

    end function getTotalKE
    

    ! --------------------------- Writing Particle Data to File -----------------------------------

    subroutine writePhaseSpace(self, dirName)
        ! Writes particle phase space into binary file
        class(Particle), intent(in) :: self
        character(*), intent(in) :: dirName
        character(len=5) :: char_i
        integer(int32) :: i
        !write(char_i, '(I3)'), CurrentDiagStep
        open(10,file=dirName//'/PhaseSpace/phaseSpace_'//self%name//".dat", form='UNFORMATTED', access = 'STREAM', status = 'REPLACE')
        do i = 1, numThread
            write(10) self%phaseSpace(:, 1:self%N_p(i), i)
        end do
        close(10)
    end subroutine writePhaseSpace

    subroutine writeLocalTemperature(self, CurrentDiagStep, dirName, numCells)
        ! Write particle temperature averaged over local grid
        class(Particle), intent(in) :: self
        integer(int32), intent(in) :: CurrentDiagStep, numCells
        character(*), intent(in) :: dirName
        character(len=5) :: char_i
        integer(int32) :: j, index, counter(numCells, numThread), iThread
        real(real64) :: temp(numCells, numThread), EHist(numCells)
        temp = 0.0d0
        counter = 0
        !$OMP parallel private(iThread, j, index) 
        iThread = omp_get_thread_num() + 1
        do j = 1, self%N_p(iThread)
            index = INT(self%phaseSpace(1, j, iThread))
            if (index < NumberXNodes) then
                temp(index, iThread) = temp(index, iThread) + SUM(self%phaseSpace(2:4, j, iThread)**2)
                counter(index, iThread) = counter(index, iThread) + 1
            end if
        end do
        !$OMP end parallel
        do j = 1, numCells
            if (SUM(counter(j, :)) > 0) then
                EHist(j) = SUM(temp(j,:))*self%mass/SUM(counter(j, :))/3.0d0/e
            else
                EHist(j) = 0.0d0
            end if
        end do
        write(char_i, '(I3)'), CurrentDiagStep
        open(10,file=dirName//'/Temperature/Temp_'//self%name//"_"//trim(adjustl(char_i))//".dat", form='UNFORMATTED')
        write(10) EHist
        close(10)
        
    end subroutine writeLocalTemperature


    ! ------------------- read and initialize particles using world type ---------------------------------

    subroutine readChargedParticleInputs(filename, irand, T_e, T_i, numThread, world, particleList)
        type(Particle), allocatable, intent(out) :: particleList(:)
        type(Domain) :: world
        character(len=*), intent(in) :: filename
        integer(int32), intent(in) :: numThread
        integer(int32), intent(in out) :: irand(numThread)
        real(real64), intent(in) :: T_e, T_i
        integer(int32) :: j, io, numSpecies = 0, numParticles(100), particleIdxFactor(100), i, tempInt, distType(100)
        character(len=15) :: name
        character(len=8) :: particleNames(100)
        real(real64) :: mass(100), charge(100), Ti(100), tempReal, alpha(100), v_drift(100)
        logical :: boolVal

        print *, "Reading particle inputs:"
        if (.not. restartBool) then
            open(10,file='../InputData/'//filename, action = 'read')
            do j=1, 10000
                read(10,*) name

                if( name(1:9).eq.'ELECTRONS') then
                    read(10,*) name
                    read(10,*) name
                    read(10,'(A4)', ADVANCE = 'NO') name(1:2)
                    numSpecies = numSpecies + 1
                    read(10,*) numParticles(numSpecies), particleIdxFactor(numSpecies), alpha(numSpecies), distType(numSpecies), v_drift(numSpecies)
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
                        read(10,*) mass(numSpecies),charge(numSpecies), numParticles(numSpecies), particleIdxFactor(numSpecies), alpha(numSpecies), distType(numSpecies), v_drift(numSpecies)
                        Ti(numSpecies) = T_i
                        mass(numSpecies) = mass(numSpecies) * m_amu - charge(numSpecies) * m_e
                        particleNames(numSpecies) = trim(name)
                        read(10,'(A6)', ADVANCE = 'NO') name
                    end do
                endif      

                if (name(1:7) == 'ENDFILE') then
                    close(10)
                    exit
                end if

            end do
        else
            open(10,file=restartDirectory//"/"//"ParticleProperties.dat", IOSTAT=io)
            read(10, *, IOSTAT = io)
            read(10, '(A)', Advance = 'NO', IOSTAT = io) name
            do while (io == 0)    
                do j = 1, len(name)
                    if (name(j:j) == ' ') then
                        exit
                    end if
                end do
                backspace(10)
                numSpecies = numSpecies + 1
                read(10, *, IOSTAT = io) name(1:j-1), mass(numSpecies), charge(numSpecies), Ti(numSpecies), particleIdxFactor(numSpecies)
                particleNames(numSpecies) = name(1:j-1)
                open(20,file=restartDirectory//"/"//"ParticleDiagnostic_"//trim(particleNames(numSpecies))//".dat", IOSTAT=io)
                read(20, *, IOSTAT = io)
                read(20, '(A)', Advance = 'NO', IOSTAT = io) name
                do while (io == 0)
                    read(20, *, IOSTAT = io)
                    read(20, '(A)', Advance = 'NO', IOSTAT = io) name
                end do
                backspace(20)
                backspace(20)
                read(20, *, IOSTAT = io) tempReal, tempReal, tempReal, tempReal, tempReal, numParticles(numSpecies)
                close(20)
                read(10, '(A)', Advance = 'NO', IOSTAT = io) name
            end do
            close(10) 

        end if

        numberChargedParticles = numSpecies
        print *, 'Amount charged particles:', numberChargedParticles
        if (numberChargedParticles > 0) then
            allocate(particleList(numberChargedParticles))
            do j=1, numberChargedParticles
                if (.not. restartBool) then
                    particleList(j) = Particle(mass(j), e * charge(j), 1.0d0, numParticles(j), numParticles(j) * particleIdxFactor(j), trim(particleNames(j)), numThread)
                    call particleList(j) % initialize_n_ave(n_ave, world%L_domain)
                    SELECT CASE(distType(j))
                    CASE(0)
                        call particleList(j)% initializeRandUniform(world, irand)
                    CASE(1)
                        call particleList(j)% initializeRandCosine(world, irand, alpha(j))
                    CASE(2)
                        call particleList(j)% initializeRandSine(world, irand, alpha(j))
                    CASE default
                        print *, 'Distribution type should be between 0 and 2!'
                        stop
                    END SELECT
                    if (j == 1) then
                        call particleList(j) % generate3DMaxwellian(Ti(j), world, irand, alpha(j), distType(j), v_drift(j))
                    else
                        call particleList(j) % generate3DMaxwellian(Ti(j), world, irand, 1.236d0, distType(j), v_drift(j))
                    end if
                else
                    particleList(j) = Particle(mass(j), charge(j), Ti(j), numParticles(j), particleIdxFactor(j), trim(particleNames(j)), numThread)
                    numParticles(j) = numParticles(j)/numThread
                    tempInt = MOD(particleList(j)%N_p(1), numThread)
                    INQUIRE(file=restartDirectory//"/PhaseSpace/phaseSpace_"//particleList(j)%name//".dat", exist = boolVal)
                    if (.not. boolVal) then
                        print *, 'Phase Space information for particle', particleList(j)%name, 'does not exist!'
                        stop
                    end if
                    open(10,file=restartDirectory//"/PhaseSpace/phaseSpace_"//particleList(j)%name//".dat", form = 'UNFORMATTED', access = 'stream', status = 'old', IOSTAT=io)
                    do i = 1, numThread
                        if (i <= tempInt) then
                            particleList(j)%N_p(i) = numParticles(j) + 1
                        else
                            particleList(j)%N_p(i) = numParticles(j)
                        end if
                        read(10,  IOSTAT = io) particleList(j)%phaseSpace(:, 1:particleList(j)%N_p(i), i)
                    end do
                    close(10)
                end if
                print *, 'Initializing ', particleList(j) % name
                print *, 'Amount of macroparticles is:', SUM(particleList(j) % N_p)
                print *, "Particle mass is:", particleList(j)%mass
                print *, "Particle charge is:", particleList(j)%q
                print *, "Particle weight is:", particleList(j)%w_p
                print *, "Particle mean KE is:", particleList(j)%getKEAve(), ", should be", Ti(j) * 1.5
                print *, 'Distribution type:', distType(j)
                print *, 'Drift velocity:', v_drift(j)
            end do
        end if

        
        print *, "---------------"
        print *, ""


    end subroutine readChargedParticleInputs

end module mod_particle