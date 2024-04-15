program MCCExample
    use iso_fortran_env, only: int32, int64, real64, output_unit
    use constants
    use omp_lib 
    use mod_BasicFunctions
    use mod_particle
    use mod_targetParticle
    use mod_NullCollision
    use mod_readInputs
    implicit none

    integer(int32), parameter :: sizeArray = 800, numSteps = 100000, numStepsAveraged = 100000
    real(real64), parameter :: E_max = 80.0d0, E_over_N = 25.0d0, bolsigCollFreq = 0.7733d-13
    integer(int32) :: i,j, iThread, energyIndx
    integer(int64) :: E_hist(sizeArray)!, tclock1, tclock2, clock_rate, 
    real(real64) :: EField, velocity_old, accelAdd, delE, E_array(sizeArray), partEnergy
    type(Particle), allocatable :: particleList(:)
    type(targetParticle), allocatable :: targetParticleList(:)
    type(nullCollision), allocatable :: nullCollisionList(:)
    !integer(int32), allocatable :: E_hist(:, :)
    numThread = 32
    call omp_set_num_threads(numThread)
    allocate(irand(numThread))
    do i = 1, numThread
        irand(i) = 123456*i*11
    end do
   
    call readParticleInputs('BoundExample.inp', numberChargedParticles, irand, 2.25d0 * 2.0d0 / 3.0d0, 0.1d0, numThread, particleList, targetParticleList)
    call readNullCollisionInputs('collision.inp', nullCollisionList, particleList, targetParticleList, numberBinaryCollisions)

    del_t = 1.d-2 / (bolsigCollFreq * targetParticleList(1)%density)
    print *, 'del_t is:', del_t
    EField = E_over_N * 1.d-21 * targetParticleList(1)%density
    print *, 'EField is:', EField
 
    
    accelAdd =  particleList(1)%q * EField * del_t/particleList(1)%mass
    print *, 'P_null is:', 1.0d0 - EXP(-nullCollisionList(1)%sigmaVMax * targetParticleList(1)%density * del_t)
    print *, 'Ave energy before accel:', particleList(1)%getKEAve()
    do i = 1, numSteps
        !$OMP parallel private(iThread, j, velocity_old)
        iThread = omp_get_thread_num() + 1
        do j = 1, particleList(1)%N_p(iThread)
            velocity_old = particleList(1)%phaseSpace(2, j, iThread)
            particleList(1)%phaseSpace(2, j, iThread) = velocity_old + accelAdd
        end do
        !$OMP end parallel
    
        call nullCollisionList(1)%generateCollision(particleList, targetParticleList, numberChargedParticles, numberBinaryCollisions, irand, del_t)
        if (MOD(i, 10000) == 0) then
            print *, real(i) * 100.0d0 /numSteps, ' % done'
        end if

    end do

    print *, 'Ave energy after accel:', particleList(1)%getKEAve()

    delE = E_max/real(sizeArray)
    do i = 1, sizeArray
        E_array(i) = real(i) * delE
    end do
    E_hist = 0
    
    nullCollisionList(1)%totalAmountCollisions = 0
    do i = 1, numSteps
        !$OMP parallel private(iThread, j, velocity_old)
        iThread = omp_get_thread_num() + 1
        do j = 1, particleList(1)%N_p(iThread)
            velocity_old = particleList(1)%phaseSpace(2, j, iThread)
            particleList(1)%phaseSpace(2, j, iThread) = velocity_old + accelAdd
        end do
        !$OMP end parallel
    
        call nullCollisionList(1)%generateCollision(particleList, targetParticleList, numberChargedParticles, numberBinaryCollisions, irand, del_t)
        !$OMP parallel private(iThread, j, partEnergy, energyIndx) reduction(+:E_hist)
        iThread = omp_get_thread_num() + 1
        do j = 1, particleList(1)%N_p(iThread)
            partEnergy = particleList(1)%mass * 0.5d0 * SUM(particleList(1)%phaseSpace(2:4, j, iThread)**2) / e
            energyIndx = INT(partEnergy/delE)
            E_hist(energyIndx) = E_hist(energyIndx) + 1
        end do
        !$OMP end parallel

        if (MOD(i, 10000) == 0) then
            print *, real(i) * 100.0d0/numSteps, ' % done'
        end if
        
    end do

    print *, 'collision frequency is:', real(nullCollisionList(1)%totalAmountCollisions)/real(numSteps)/del_t/SUM(particleList(1)%N_p)
    print *, 'should be about:', bolsigCollFreq * targetParticleList(1)%density

    open(10,file='Hist_25Td_He.dat')
    do i = 1, sizeArray
        write(10, "((es16.8,1x), (I14, 1x))") E_array(i), E_hist(i)
    end do
    close(10)
    
    E_hist = 0
    !$OMP parallel private(iThread, j, partEnergy, energyIndx) reduction(+:E_hist)
    iThread = omp_get_thread_num() + 1
    do j = 1, particleList(1)%N_p(iThread)
        partEnergy = particleList(1)%mass * 0.5d0 * SUM(particleList(1)%phaseSpace(2:4, j, iThread)**2) / e
        energyIndx = INT(partEnergy/delE)
        E_hist(energyIndx) = E_hist(energyIndx) + 1
    end do
    !$OMP end parallel

    
end program MCCExample