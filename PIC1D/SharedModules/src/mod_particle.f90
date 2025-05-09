module mod_particle

    use iso_fortran_env, only: int32, real64, output_unit
    use constants
    use mod_BasicFunctions
    use mod_domain
    use iso_c_binding
    use omp_lib
    implicit none

    private
    public :: Particle, readChargedParticleInputs
    integer(int32), public, protected :: numberChargedParticles = 0 ! number of charged particles total

    ! Particle contains particle properties and stored values in phase space for each charged particle
    type :: Particle
        character(:), allocatable :: name !name of the particle
        integer(int32), allocatable :: N_p(:), refIdx(:), refRecordIdx(:, :), delIdx(:), numToCollide(:), startIdx(:) !N_p is the current last index of particle, refIdx and delIdx are number to reflux and delete
        ! numToCollide saves number that can collided in nullCollision
        integer(int32) :: finalIdx, accumWallLoss(2) ! finalIdx is maximum particles per thread
        real(real64), allocatable :: phaseSpace(:,:, :) !particle phase space, represents [l_x, v_x, v_y, v_z] in first index
        real(real64) :: mass, q, w_p, q_over_m, q_times_wp ! mass (kg), charge(C), and weight (N/m^2 in 1D) of particles. Assume constant weight for moment
        real(real64), allocatable :: wallLoss(:,:), energyLoss(:,:), momentumLoss(:,:) ! Accumulate number of particles, v^2, and v of particles at wall
        real(real64), allocatable :: densities(:, :), workSpace(:,:) ! Per thread densities and general workspace over domain
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
        ! Construct particle object
        real(real64), intent(in) :: mass, q, w_p
        integer(int32), intent(in) :: N_p, finalIdx, numThread
        character(*), intent(in) :: particleName
        integer(int32) :: maxNumNodes
        integer(int32) :: iThread
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
        maxNumNodes = MAX(NumberXNodes, NumberXHalfNodes) ! Take maximum of nodes on domain
        allocate(self%phaseSpace(4,finalIdx, numThread), self%refRecordIdx(INT(self%finalIdx/10), numThread), self%N_p(numThread), &
            self%delIdx(numThread), self%wallLoss(2, numThread), self%energyLoss(2, numThread), self%refIdx(numThread), &
            self%densities(NumberXNodes, numThread), self%workSpace(maxNumNodes, numThread), self%numToCollide(numThread), self%momentumLoss(2, numThread), self%startIdx(numThread))
        
        !First touch initiation
        !$OMP parallel private(iThread)
        iThread = omp_get_thread_num() + 1
        self%phaseSpace(:,:, iThread) = 0
        self%refRecordIdx(:,iThread) = 0
        self%N_p(iThread) = N_p
        self%energyLoss(:,iThread) = 0
        self%wallLoss(:,iThread) = 0
        self%momentumLoss(:,iThread) = 0
        self%startIdx(iThread) = 1
        self%delIdx(iThread) = 1
        self%workSpace(:,iThread) = 0
        self%numToCollide(iThread) = 0
        self%densities(:, iThread) = 0
        !$OMP end parallel
    end function particle_constructor

    pure subroutine initialize_n_ave(self, n_ave, L_domain)
        ! initialize w_p based on initial average n (N/m^3) over domain
        class(Particle), intent(in out) :: self
        real(real64), intent(in) :: n_ave, L_domain
        self % w_p = n_ave * L_domain / SUM(self % N_p)
        self%q_times_wp = self%q * self%w_p
    end subroutine initialize_n_ave

    subroutine initializeRandUniform(self, world, irand)
        ! distribute particles randomly over the domain
        class(Particle), intent(in out) :: self
        type(Domain), intent(in) :: world
        integer(c_int64_t), intent(in out) :: irand(numThread)
        integer(int32) :: iThread, i
        integer(c_int64_t) :: irand_thread
        real(real64) :: x_pos
        !$OMP parallel private(iThread, i, x_pos, irand_thread)
        iThread = omp_get_thread_num() + 1
        irand_thread = irand(iThread)
        do i = 1, self%N_p(iThread)
            ! get random position within the domain
            x_pos = pcg32_random_r(irand_thread) * world%L_domain + world%startX
            ! convert to logical space
            self%phaseSpace(1,i, iThread) = world%getLFromX(x_pos)
        end do
        irand(iThread) = irand_thread
        !$OMP end parallel
    end subroutine initializeRandUniform

    subroutine initializeRandCosine(self, world, irand, alpha)
        ! Use acceptance-rejection method to distribute particles
        ! f(x) = 1 + alpha * Cos(2 pi x/L)
        class(Particle), intent(in out) :: self
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: alpha
        integer(c_int64_t), intent(in out) :: irand(numThread)
        integer(int32) :: iThread, i
        integer(c_int64_t) :: irand_thread
        real(real64) :: Rand1, Rand2, x_pos, val
        !$OMP parallel private(iThread, i, irand_thread, Rand1, Rand2, x_pos, val)
        iThread = omp_get_thread_num() + 1
        irand_thread = irand(iThread)
        do i = 1, self%N_p(iThread)
            Rand1 = pcg32_random_r(irand_thread)
            Rand2 = pcg32_random_r(irand_thread)
            val = (1.0d0 + alpha * COS(2 * pi * Rand2))/(1.0d0 + alpha)
            do while (Rand1 > val)
                Rand1 = pcg32_random_r(irand_thread)
                Rand2 = pcg32_random_r(irand_thread)
                val = (1.0d0 + alpha * COS(2 * pi * Rand2))/(1.0d0 + alpha)
            end do
            x_pos = Rand2 * world%L_domain + world%startX
            self%phaseSpace(1,i, iThread) = world%getLFromX(x_pos)
        end do
        irand(iThread) = irand_thread
        !$OMP end parallel
    end subroutine initializeRandCosine

    subroutine initializeRandSine(self, world, irand, alpha)
        ! Use acceptance-rejection method to distribute particles
        ! f(x) = 1 + alpha * Sin(2 pi x/L)
        class(Particle), intent(in out) :: self
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: alpha
        integer(c_int64_t), intent(in out) :: irand(numThread)
        integer(int32) :: iThread, i
        integer(c_int64_t) :: irand_thread
        real(real64) :: Rand1, Rand2, x_pos, val
        !$OMP parallel private(iThread, i, irand_thread, Rand1, Rand2, x_pos, val)
        iThread = omp_get_thread_num() + 1
        irand_thread = irand(iThread)
        do i = 1, self%N_p(iThread)
            Rand1 = pcg32_random_r(irand_thread)
            Rand2 = pcg32_random_r(irand_thread)
            val = (1.0d0 + alpha * SIN(2 * pi * Rand2))/(1.0d0 + alpha)
            do while (Rand1 > val)
                Rand1 = pcg32_random_r(irand_thread)
                Rand2 = pcg32_random_r(irand_thread)
                val = (1.0d0 + alpha * SIN(2 * pi * Rand2))/(1.0d0 + alpha)
            end do
            x_pos = Rand2 * world%L_domain + world%startX
            self%phaseSpace(1,i, iThread) = world%getLFromX(x_pos)
        end do
        irand(iThread) = irand_thread
        !$OMP end parallel
    end subroutine initializeRandSine

    ! ------------------------ Generating and initializing velocity with Maxwellian --------------------------


    subroutine generate3DMaxwellian(self, T, world, irand, alpha, distType, v_drift)
        ! random velocity generator for the particle for temperature T (eV)
        ! Use box-muller method for random guassian variable
        ! apply to all particles
        class(Particle), intent(in out) :: self
        class(Domain), intent(in) :: world
        real(real64), intent(in) :: T, alpha, v_drift
        integer(int32), intent(in) :: distType
        integer(c_int64_t), intent(in out) :: irand(numThread)
        integer(int32) :: iThread, i
        integer(c_int64_t) :: irand_thread
        real(real64) :: U1, U2, U3, U4, v_therm, v_extra, x_pos
        v_therm = SQRT(T*e/self%mass)
        !$OMP PARALLEL PRIVATE(iThread, i, U1, U2, U3, U4, irand_thread, v_extra, x_pos)
        iThread = omp_get_thread_num() + 1
        irand_thread = irand(iThread)
        do i = 1, self%N_p(iThread)
            U1 = pcg32_random_r(irand_thread)
            U2 = pcg32_random_r(irand_thread)
            U3 = pcg32_random_r(irand_thread)
            U4 = pcg32_random_r(irand_thread)
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
        ! Get total momentum in domain
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
        integer(int32) :: iThread
        !write(char_i, '(I3)'), CurrentDiagStep
        !$OMP parallel private(iThread, char_i)
        iThread = omp_get_thread_num() + 1
        write(char_i, '(I3)'), iThread
        open(iThread,file=dirName//'/PhaseSpace/phaseSpace_'//self%name//"_thread"//trim(adjustl(char_i))//".dat", form='UNFORMATTED', access = 'STREAM', status = 'REPLACE')
        write(iThread) self%phaseSpace(:, 1:self%N_p(iThread), iThread)
        close(ithread)
        !$OMP end parallel
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
        ! Read input file for particles
        type(Particle), allocatable, intent(out) :: particleList(:)
        type(Domain) :: world
        character(len=*), intent(in) :: filename
        integer(int32), intent(in) :: numThread
        integer(c_int64_t), intent(in out) :: irand(numThread)
        real(real64), intent(in) :: T_e, T_i
        integer(int32) :: j, io, numSpecies = 0, numParticles(100), particleIdxFactor(100), i, tempInt, distType(100), iThread
        character(len=15) :: name
        character(len=8) :: particleNames(100), char_i
        real(real64) :: mass(100), charge(100), Ti(100), tempReal, alpha(100), v_drift(100)
        logical :: boolVal

        print *, "Reading particle inputs:"
        if (.not. restartBool) then
            open(10,file='../InputData/'//filename, action = 'read')
        else
            open(10,file=restartDirectory//'/InputData/'//filename, action = 'read')
        end if
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
                exit
            end if

        end do
        close(10)
        if (mass(1) /= m_e) call changeIonStep(1) ! if no electron, then you have del_t over only ions
        numberChargedParticles = numSpecies
        print *, 'Amount charged particles:', numberChargedParticles
        if (numberChargedParticles > 0) then
            ! Initialize and generate particles
            allocate(particleList(numberChargedParticles))
            do j=1, numberChargedParticles
                particleList(j) = Particle(mass(j), e * charge(j), 1.0d0, numParticles(j), numParticles(j) * particleIdxFactor(j), trim(particleNames(j)), numThread)
                call particleList(j) % initialize_n_ave(n_ave, world%L_domain)
                if (.not. restartBool) then
                    if (j==2 .and. numberChargedParticles == 2 .and. charge(2) == -charge(1) .and. numParticles(1) == numParticles(2)) then
                        ! If only ions and electrons (electrons come first) then set ion positions same as electrons for neutral start
                        print *, 'Neutral charge start!'
                        particleList(2)%phaseSpace(1, :, :) = particleList(1)%phaseSpace(1,:,:)
                    else
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
                    end if
                    if (j == 1) then
                        call particleList(j) % generate3DMaxwellian(Ti(j), world, irand, alpha(j), distType(j), v_drift(j))
                    else
                        call particleList(j) % generate3DMaxwellian(Ti(j), world, irand, 1.236d0, distType(j), v_drift(j)) ! Used for IASW case, will likely change to general later
                    end if
                else
                    !$OMP parallel private(iThread, boolVal, char_i, io, i)
                    iThread = omp_get_thread_num() + 1
                    write(char_i, '(I3)'), iThread
                    INQUIRE(file=restartDirectory//'/PhaseSpace/phaseSpace_'//particleList(j)%name//"_thread"//trim(adjustl(char_i))//".dat", exist = boolVal)
                    if (.not. boolVal) then
                        print *, restartDirectory//'/PhaseSpace/phaseSpace_'//particleList(j)%name//"_thread"//trim(adjustl(char_i))//".dat", 'Does not exist'
                        stop
                    end if
                    open(iThread,file=restartDirectory//"/PhaseSpace/phaseSpace_"//particleList(j)%name//"_thread"//trim(adjustl(char_i))//".dat", form = 'UNFORMATTED', access = 'stream', status = 'old', IOSTAT=io)
                    i = 0
                    read(iThread,  IOSTAT = io) particleList(j)%phaseSpace(:, i+1, iThread)
                    do while (io == 0)
                        i = i + 1
                        read(iThread,  IOSTAT = io) particleList(j)%phaseSpace(:, i+1, iThread)
                    end do
                    particleList(j)%N_p(iThread) = i
                    close(iThread)
                    !$OMP end parallel
                end if
                print *, 'Initializing ', particleList(j) % name
                print *, 'Amount of macroparticles is:', SUM(particleList(j) % N_p)
                print *, "Particle mass is:", particleList(j)%mass
                print *, "Particle charge is:", particleList(j)%q
                print *, "Particle weight is:", particleList(j)%w_p
                print *, "Particle mean KE is:", particleList(j)%getKEAve()
                print *, 'Distribution type:', distType(j)
                print *, 'Drift velocity:', v_drift(j)
            end do
        end if

        
        print *, "---------------"
        print *, ""


    end subroutine readChargedParticleInputs

end module mod_particle