program BoundPlasmaExample
    use iso_fortran_env, only: int32, real64, output_unit, int64
    use constants
    use mod_BasicFunctions
    use mod_domain
    use mod_particle
    ! use mod_potentialSolver
    ! use mod_particleMover
    ! use mod_collisions
    ! use mod_nonLinSolvers
    ! use mod_Scheme
    ! use mod_simulation
    implicit none

    integer(int32) :: i
    real(real64), allocatable :: densities(:)
    type(Domain) :: globalWorld
    type(Particle), allocatable :: globalParticleList(:)
    call readInitialConditions('InitialConditions.inp')
    globalWorld = Domain('Geometry.inp')
    globalParticleList =  readParticleInputs('BoundExample.inp',numberChargedParticles, irand)
    do i = 1, numberChargedParticles
        call globalParticleList(i)%initialize_randUniform(globalWorld, irand)
    end do
    allocate(densities(NumberXNodes))
    densities = 0.0d0
    call globalParticleList(1)%loadParticleDensity(densities, globalWorld)
    print *, densities/globalWorld%nodeVol
    !real(real64) :: remainDel_t, currDel_t
    ! call initializeScheme(schemeNum)
    ! ! Initialize constants with inputs
    ! ! create the world the particles live in
    ! call readInputs(NumberXNodes, numDiagnosticSteps, averagingTime, fractionFreq, n_ave, globalWorld, globalSolver, simulationTime, Power, heatSkipSteps, nu_h, T_e, T_i, 'Geometry.inp', 'InitialConditions.inp')
    ! globalParticleList = readParticleInputs('BoundExample.inp',numberChargedParticles, irand, T_e, T_i) 
    ! do i = 1, numberChargedParticles
    !     call initialize_randUniform(globalParticleList(i), globalWorld, irand)
    !     call globalParticleList(i) % initialize_n_ave(n_ave, globalWorld%grid(NumberXNodes) - globalWorld%grid(1))
    ! end do
    ! call initializeSolver(eps_r, solverType, m_Anderson, Beta_k, maxIter)
    ! print *, "Calulated values:"
    ! print *, "Number of particles is:", globalParticleList(1)%N_p
    ! print *, "w_p is:", globalParticleList(1)%w_p
    ! print *, "Debye length is:", getDebyeLength(globalParticleList(1)%getKEAve()*2.0d0/3.0d0, n_ave)
    ! print *, "Plasma frequency is:", getPlasmaFreq(n_ave)
    ! print *, "Average density is ", globalParticleList(1)%N_p * globalParticleList(1)%w_p / (globalWorld%grid(NumberXNodes) - globalWorld%grid(1)), "should be", n_ave
    ! del_t = fractionFreq/getPlasmaFreq(n_ave)   
    ! print *, "Time step (sec) is:", del_t
    ! print *, "----------------"
    ! ! Generate solver object, and then solve for initial rho/potential
    ! call solveInitialPotential(globalSolver, globalParticleList, globalWorld)
    ! ! currDel_t = del_t
    ! ! remainDel_t = del_t
    ! ! call solvePotential(globalSolver, globalParticleList, globalWorld, del_t, remainDel_t, currDel_t, maxIter, eps_r)
    ! ! print *, "Took", iterNumPicard, "iterations"
    ! ! print *, 'electron number:', globalParticleList(1)%N_P
    ! ! print *, 'ion number:', globalParticleList(2)%N_p
    ! ! call addMaxwellianLostParticles(globalParticleList, T_e, T_i, irand, globalWorld)
    ! ! print *, 'electron number:', globalParticleList(1)%N_P
    ! ! print *, 'ion number:', globalParticleList(2)%N_p
    ! ! ! Get error gauss' law
    ! ! call depositRho(globalSolver%rho, globalParticleList, globalWorld)
    ! ! call globalSolver%construct_diagMatrix(globalWorld)
    ! ! chargeError = globalSolver%getError_tridiag_Poisson(globalWorld)
    ! ! !chargeError = chargeError / SQRT(SUM(globalSolver%rho**2))
    ! ! call globalSolver%construct_diagMatrix_Ampere(globalWorld)
    ! ! print *, "Charge error is:", chargeError
    ! ! stop 
    ! call solveSimulation(globalSolver, globalParticleList, globalWorld, del_t, maxIter, eps_r, irand, simulationTime, heatSkipSteps)

    ! print *, "Averaging up to", averagingTime, "simulation seconds"
    ! call solveSimulationFinalAverage(globalSolver, globalParticleList, globalWorld, del_t, maxIter, eps_r, irand, averagingTime, 100)

contains

    function readParticleInputs(filename, numberChargedParticles, irand) result(particleList)
        type(Particle), allocatable :: particleList(:)
        character(len=*), intent(in) :: filename
        integer(int32), intent(in out) :: numberChargedParticles, irand
        integer(int32) :: j, numSpecies = 0, numParticles(100), particleIdxFactor(100)
        character(len=15) :: name
        character(len=8) :: particleNames(100)
        real(real64) :: mass(100), charge(100), Ti(100)

        print *, "Reading particle inputs:"
        open(10,file='../../SharedModules/InputData/'//filename, action = 'read')

        do j=1, 10000
            read(10,*,END=101,ERR=100) name

            if( name(1:9).eq.'ELECTRONS') then
                read(10,*,END=101,ERR=100) name
                read(10,*,END=101,ERR=100) name
                read(10,'(A2)',END=101,ERR=100, ADVANCE = 'NO') name(1:2)
                numSpecies = numSpecies + 1
                read(10,*,END=101,ERR=100) numParticles(numSpecies), particleIdxFactor(numSpecies)
                Ti(numSpecies) = T_e
                mass(numSpecies) = m_e
                charge(numSpecies) = -1.0
                particleNames(numSpecies) = 'e'
                read(10,*,END=101,ERR=100) name
                read(10,*,END=101,ERR=100) name
            endif


            if(name(1:4).eq.'IONS' .or. name(1:4).eq.'Ions' .or. name(1:4).eq.'ions' ) then
                do while(name(1:4).ne.'----')
                    read(10,*,END=101,ERR=100) name
                end do
    200             read(10,'(A6)',END=101,ERR=100, ADVANCE = 'NO') name
                if (name(1:4).eq.'----') then
                    close(10)
                else
                    numSpecies = numSpecies + 1
                    read(10,*,END=101,ERR=100) mass(numSpecies),charge(numSpecies), numParticles(numSpecies), particleIdxFactor(numSpecies)
                    Ti(numSpecies) = T_i
                    mass(numSpecies) = mass(numSpecies) * m_p
                    particleNames(numSpecies) = trim(name)
                    goto 200
                end if
            endif
            ! Take care of extra text I guess        

            if (name(1:7) == 'ENDFILE') then
                close(10)
            end if

        end do
    100     continue
    101     continue
        numberChargedParticles = numSpecies
        allocate(particleList(numberChargedParticles))
        do j=1, numberChargedParticles
            particleList(j) = Particle(mass(j), e * charge(j), numParticles(j), numParticles(j) * particleIdxFactor(j), trim(particleNames(j)))
            call particleList(j) % generate3DMaxwellian(Ti(j), irand)
            print *, 'Initializing ', particleList(j) % name
            print *, 'Number of particles is:', SUM(particleList(j)%N_p)
            print *, "Particle mass is:", particleList(j)%mass
            print *, "Particle charge is:", particleList(j)%q
            print *, "Particle mean KE is:", particleList(j)%getKEAve(), ", should be", Ti(j) * 1.5
        end do
        
        print *, "---------------"
        print *, ""



    end function readParticleInputs

    
end program BoundPlasmaExample