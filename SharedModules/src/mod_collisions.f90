module mod_collisions
    use iso_fortran_env, only: int32, real64
    use constants
    use mod_BasicFunctions
    use mod_domain
    use mod_particle
    implicit none

    real(real64) :: inelasticEnergyLoss = 0.0d0
    real(real64) :: Power = 100.d0 !W/m^2 in 1D
    real(real64) :: nu_h = 1.0d8, fraction = 0.1d0 !Hz
    integer(int32) :: heatSkipSteps = 1
    ! Type of collisions, will likely need arrays which you can loop through which has source, target particle, choose from particle list
    ! Will likely have construct method which takes in collision data and turns into array

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
            particleList(j) = Particle(mass(j), e * charge(j), numParticles(j), particleIdxFactor(j), trim(particleNames(j)))
            call particleList(j) % generate3DMaxwellian(Ti(j), irand)
            print *, 'Initializing ', particleList(j) % name
            print *, 'Number of particles is:', SUM(particleList(j)%N_p)
            print *, "Particle mass is:", particleList(j)%mass
            print *, "Particle charge is:", particleList(j)%q
            print *, "Particle mean KE is:", particleList(j)%getKEAve(), ", should be", Ti(j) * 1.5d0
        end do
        
        print *, "---------------"
        print *, ""



    end function readParticleInputs

    subroutine addMaxwellianLostParticles(particleList, T_e, T_i, irand, world)
        ! Add power to all particles in domain
        type(Particle), intent(in out) :: particleList(2)
        type(Domain), intent(in) :: world
        integer(int32), intent(in out) :: irand
        real(real64), intent(in) :: T_e, T_i
        integer(int32) :: i,j, k, numToMerge, numNonMerged, numToSplit, numToSplitCell, idxRan, l_int
        real(real64) :: l_random, l_merged, w_merged, Temp(2), l_del, sumElectronLoss, w_coll, l_array(NumberXNodes-1), w_array(NumberXNodes-1)
        Temp(1) = T_e
        Temp(2) = T_i
        do j = 1, 2
            do i = 1, particleList(j)%refIdx
                l_int = INT(particleList(j)%refPhaseSpace(1, i))
                particleList(j)%N_p(l_int) = particleList(j)%N_p(l_int) + 1
                particleList(j)%phaseSpace(1, particleList(j)%N_p(l_int), l_int) = particleList(j)%refPhaseSpace(1, i)
                particleList(j)%w_p(particleList(j)%N_p(l_int), l_int) = particleList(j)%refw_p(i)
                call getMaxwellianFluxSample(particleList(j)%phaseSpace(2:4, particleList(j)%N_p(l_int), l_int), particleList(j)%mass, Temp(j), irand)
                particleList(j)%phaseSpace(2:4, particleList(j)%N_p(l_int), l_int) = -ABS(particleList(j)%phaseSpace(2:4, particleList(j)%N_p(l_int), l_int))
            end do
            particleList(j)%refIdx = 0
        end do
        sumElectronLoss = SUM(particleList(1)%wallLoss)
        if (sumElectronLoss > 0) then
            do i = 1, NumberXNodes-1
                l_random = ran2(irand) + real(i)
                do j = 1, 2
                    particleList(j)%N_p(i) = particleList(j)%N_p(i) + 1
                    particleList(j)%phaseSpace(1, particleList(j)%N_p(i), i) = l_random
                    particleList(j)%w_p(particleList(j)%N_p(i), i) = sumElectronLoss * world%dx_dl(i)/ (world%grid(NumberXNodes))
                end do
                call getMaxwellianSample(particleList(1)%phaseSpace(2:4, particleList(1)%N_p(i), i), particleList(1)%mass, T_e, irand)
                call getMaxwellianFluxSample(particleList(2)%phaseSpace(2:4, particleList(2)%N_p(i), i), particleList(2)%mass, T_i, irand)
            end do
        end if
        ! w_coll = sumElectronLoss / particleList(1)%delIdx
        ! do i=1, particleList(1)%delIdx
        !     l_random = world%getLFromX(world%grid(NumberXNodes) * ran2(irand))
        !     l_int = INT(l_random)
        !     do j = 1, 2
        !         particleList(j)%N_p(l_int) = particleList(j)%N_p(l_int) + 1
        !         particleList(j)%phaseSpace(1, particleList(j)%N_p(l_int), l_int) = l_random
        !         particleList(j)%w_p(particleList(j)%N_p(l_int), l_int) = w_coll
        !     end do
        !     call getMaxwellianSample(particleList(1)%phaseSpace(2:4, particleList(1)%N_p(l_int), l_int), particleList(1)%mass, T_e, irand)
        !     call getMaxwellianFluxSample(particleList(2)%phaseSpace(2:4, particleList(2)%N_p(l_int), l_int), particleList(2)%mass, T_i, irand)
        ! end do
    end subroutine addMaxwellianLostParticles

    ! subroutine addMaxwellianLostParticles(particleList, T_e, T_i, irand, world)
    !     ! Add power to all particles in domain
    !     type(Particle), intent(in out) :: particleList(2)
    !     type(Domain), intent(in) :: world
    !     integer(int32), intent(in out) :: irand
    !     real(real64), intent(in) :: T_e, T_i
    !     integer(int32) :: i,j, k, numToMerge, numNonMerged, numToSplit, numToSplitCell, idxRan, l_int
    !     real(real64) :: l_random, l_merged, w_merged, Temp(2), l_del, sumElectronLoss, w_coll, l_array(NumberXNodes-1), w_array(NumberXNodes-1)
    !     Temp(1) = T_e
    !     Temp(2) = T_i
    !     do j = 1, 2
    !         if (particleList(j)%refIdx > 0) then
    !             numToMerge = MIN(particleList(j)%refIdx - particleList(j)%partPerCell + particleList(j)%N_p(NumberXNodes-1) + 1, particleList(j)%refIdx)
    !             if (numToMerge > 1) then
    !                 numNonMerged = particleList(j)%refIdx - numToMerge
    !                 numToSplit = numToMerge - 1
    !                 do k = 1, NumberXNodes-1
    !                     if (particleList(j)%N_p(k) > 0) then
    !                         if (particleList(j)%partPerCell > particleList(j)%N_p(k)) then
    !                             numToSplitCell = MIN(numToSplit, particleList(j)%partPerCell - particleList(j)%N_p(k))
    !                             do i = 1, numToSplitCell
    !                                 idxRan = INT(ran2(irand) * particleList(j)%N_p(k) + 1)
    !                                 particleList(j)%phaseSpace(2:4,particleList(j)%N_p(k) + i, k) = particleList(j)%phaseSpace(2:4,idxRan, k)
    !                                 particleList(j)%w_p(idxRan, k) = 0.5d0 * particleList(j)%w_p(idxRan, k)
    !                                 particleList(j)%w_p(particleList(j)%N_p(k) + i, k) = particleList(j)%w_p(idxRan, k)
    !                                 l_del = MIN(ABS(particleList(j)%phaseSpace(1,idxRan, k)-k), ABS(particleList(j)%phaseSpace(1,idxRan, k)-k - 1))
    !                                 l_del = ran2(irand) * (l_del - 1d-8)
    !                                 particleList(j)%phaseSpace(1,particleList(j)%N_p(k) + i, k) = particleList(j)%phaseSpace(1,idxRan, k) + l_del
    !                                 particleList(j)%phaseSpace(1,idxRan, k) = particleList(j)%phaseSpace(1,idxRan, k) - l_del
    !                             end do
    !                             particleList(j)%N_p(k) = particleList(j)%N_p(k) + numToSplitCell
    !                             numToSplit = numToSplit - numToSplitCell
    !                             if (numToSplit == 0) exit
    !                         end if
    !                     end if
    !                 end do
    !                 w_merged = SUM(particleList(j)%refw_p(numNonMerged+1:particleList(j)%refIdx))
    !                 l_merged = SUM((particleList(j)%refPhaseSpace(1, numNonMerged+1:particleList(j)%refIdx) - real(NumberXNodes-1)) * particleList(j)%refw_p(numNonMerged+1:particleList(j)%refIdx) )
    !                 l_merged = l_merged / w_merged + real(NumberXNodes-1)
    !                 particleList(j)%N_p(NumberXNodes-1) = particleList(j)%N_p(NumberXNodes-1) + 1
    !                 particleList(j)%phaseSpace(1, particleList(j)%N_p(NumberXNodes-1), NumberXNodes-1) = l_merged
    !                 particleList(j)%w_p(particleList(j)%N_p(NumberXNodes-1), NumberXNodes-1) = w_merged
    !                 call getMaxwellianFluxSample(particleList(j)%phaseSpace(2:4, particleList(j)%N_p(NumberXNodes-1), NumberXNodes-1), particleList(j)%mass, Temp(j), irand)
    !                 particleList(j)%phaseSpace(2:4, particleList(j)%N_p(NumberXNodes-1), NumberXNodes-1) = -ABS(particleList(j)%phaseSpace(2:4, particleList(j)%N_p(NumberXNodes-1), NumberXNodes-1))
    !             else
    !                 numNonMerged = particleList(j)%refIdx
    !             end if
    !             do i = 1, numNonMerged
    !                 particleList(j)%N_p(NumberXNodes-1) = particleList(j)%N_p(NumberXNodes-1) + 1
    !                 particleList(j)%w_p(particleList(j)%N_p(NumberXNodes-1), NumberXNodes-1) = particleList(j)%refw_p(i)
    !                 particleList(j)%phaseSpace(1,particleList(j)%N_p(NumberXNodes-1), NumberXNodes-1) = particleList(j)%refPhaseSpace(1, i)
    !                 call getMaxwellianFluxSample(particleList(j)%phaseSpace(2:4, particleList(j)%N_p(NumberXNodes-1), NumberXNodes-1), particleList(j)%mass, Temp(j), irand)
    !                 particleList(j)%phaseSpace(2:4, particleList(j)%N_p(NumberXNodes-1), NumberXNodes-1) = -ABS(particleList(j)%phaseSpace(2:4, particleList(j)%N_p(NumberXNodes-1), NumberXNodes-1))
    !             end do
    !         end if
    !         particleList(j)%refIdx = 0
    !     end do
    !     sumElectronLoss = SUM(particleList(1)%wallLoss)
    !     w_coll = sumElectronLoss/particleList(1)%delIdx
    !     l_array = 0.0
    !     w_array = 0.0
    !     do i=1, particleList(1)%delIdx
    !         l_random = world%getLFromX(world%grid(NumberXNodes) * ran2(irand))
    !         l_int = INT(l_random)
    !         particleList(1)%N_p(l_int) = particleList(1)%N_p(l_int) + 1
    !         particleList(1)%phaseSpace(1, particleList(1)%N_p(l_int), l_int) = l_random
    !         particleList(1)%w_p(particleList(1)%N_p(l_int), l_int) = w_coll
    !         call getMaxwellianSample(particleList(1)%phaseSpace(2:4, particleList(1)%N_p(l_int), l_int), particleList(1)%mass, T_e, irand)
    !         l_array(l_int) = l_array(l_int) + (l_random - l_int) * w_coll
    !         w_array(l_int) = w_array(l_int) + w_coll
    !         !call getMaxwellianFluxSample(particleList(2)%phaseSpace(2:4, particleList(2)%N_p(l_int), l_int), particleList(2)%mass, T_i, irand)
    !     end do
    !     do i = 1, NumberXNodes-1
    !         if (w_array(i) > 0.0) then
    !             particleList(2)%N_p(i) = particleList(2)%N_p(i) + 1
    !             particleList(2)%phaseSpace(1, particleList(2)%N_p(i), i) = l_array(i)/w_array(i) + i
    !             particleList(2)%w_p(particleList(2)%N_p(i), i) = w_array(i)
    !             call getMaxwellianFluxSample(particleList(2)%phaseSpace(2:4, particleList(2)%N_p(i), i), particleList(2)%mass, T_i, irand)
    !         end if
    !     end do
    ! end subroutine addMaxwellianLostParticles



end module mod_collisions