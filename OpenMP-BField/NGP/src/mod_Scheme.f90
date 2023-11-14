module mod_Scheme
    ! module which has NGP specific procedures, typically having to do with the different domain grid structure
    use iso_fortran_env, only: int32, real64
    use constants
    use mod_BasicFunctions
    use mod_particle
    use mod_domain
    use omp_lib
    implicit none

contains

    subroutine initializeScheme(schemeNum)
        integer(int32), intent(in out) :: schemeNum
        schemeNum = 0
    end subroutine initializeScheme

    subroutine initialize_randUniform(part, world, irand)
        ! place particles randomly in each dx_dl based on portion of volume it take up
        type(Particle), intent(in out) :: part
        type(Domain), intent(in) :: world
        integer(int32), intent(in out) :: irand(numThread)
        integer(int32) :: i, iThread
        real(real64) :: L_domain
        L_domain = world%grid(NumberXNodes) - world%grid(1)
        !$OMP parallel private(iThread, i)
        iThread = omp_get_thread_num() + 1
        do i=1, part%N_p(iThread)
            part%phaseSpace(1,i, iThread) = ran2(irand(iThread)) * L_domain + world%grid(1)
            part%phaseSpace(1,i, iThread) = world%getLFromX(part%phaseSpace(1,i, iThread))
        end do
        !$OMP end parallel    
    end subroutine initialize_randUniform

    subroutine initialize_QuasiNeutral(ions, electrons)
        type(Particle), intent(in) :: electrons
        type(Particle), intent(in out) :: ions
        integer(int32) :: i, iThread
        !$OMP parallel private(iThread, i)
        iThread = omp_get_thread_num() + 1
        do i=1, ions%N_p(iThread)
            ions%phaseSpace(1,i, iThread) = electrons%phaseSpace(1,i, iThread)
        end do
        !$OMP end parallel  
    end subroutine initialize_QuasiNeutral

    function getLFromX(x, world) result(l)
        real(real64), intent(in) :: x
        type(Domain), intent(in) :: world
        integer(int32) :: idxLower, idxHigher, idxMiddle
        real(real64) :: l
        idxLower = 1
        idxHigher = NumberXNodes
        if ((x<world%grid(1)) .or. (x > world%grid(NumberXNodes))) then
            print *, 'x value outside of grid range in getLFromX!'
            stop
        end if
        do while (idxLower /= idxHigher-1)
            idxMiddle = (idxLower + idxHigher)/2
            if (world%grid(idxMiddle) <= x) then
                idxLower = idxMiddle
            else
                idxHigher = idxMiddle
            end if
        end do
        l = idxLower + (x - world%grid(idxLower))/world%dx_dl(idxLower)
        
    end function getLFromX

    ! --------------------------- Diagnostics ------------------------------------

    subroutine depositRho(rho, particleList, world) 
        type(Domain), intent(in) :: world
        type(Particle), intent(in) :: particleList(:)
        real(real64), intent(in out) :: rho(NumberXNodes)
        integer(int32) :: i, j, l_left, l_right, iThread
        real(real64) :: d, rhoTemp(NumberXNodes, numThread)
        rho = 0.0d0
        do i=1, numberChargedParticles
            rhoTemp = 0.0d0
            !$OMP parallel private(iThread, j, l_left, l_right, d)
            iThread = omp_get_thread_num() + 1
            do j = 1, particleList(i)%N_p(iThread)
                l_left = INT(particleList(i)%phaseSpace(1, j, iThread))
                l_right = l_left + 1
                d = particleList(i)%phaseSpace(1, j, iThread) - real(l_left)
                rhoTemp(l_left, iThread) = rhoTemp(l_left, iThread) + (1.0d0-d)
                if (l_left < NumberXNodes) rhoTemp(l_right, iThread) = rhoTemp(l_right, iThread) + d
            end do
            !$OMP end parallel
            rho = rho + SUM(rhoTemp, DIM = 2) * particleList(i)%w_p * particleList(i)%q
        end do
    end subroutine depositRho

    subroutine loadParticleDensity(densities, particleList, world)
        type(Particle), intent(in) :: particleList(numberChargedParticles)
        type(Domain), intent(in) :: world
        real(real64), intent(in out) :: densities(NumberXNodes,numberChargedParticles)
        integer(int32) :: i,j, l_left, l_right, iThread
        real(real64) :: d, tempDensity(NumberXNodes, numThread)
        do i=1, numberChargedParticles
            tempDensity = 0.0d0
            !$OMP parallel private(iThread, j, l_left, l_right, d)
            iThread = omp_get_thread_num() + 1
            do j = 1, particleList(i)%N_p(iThread)
                l_left = INT(particleList(i)%phaseSpace(1,j, iThread))
                l_right = l_left + 1
                d = particleList(i)%phaseSpace(1,j, iThread) - l_left
                tempDensity(l_left, iThread) = tempDensity(l_left, iThread) + (1.0d0-d)
                if (l_left < NumberXNodes) tempDensity(l_right, iThread) = tempDensity(l_right, iThread) + d
            end do
            !$OMP end parallel
            densities(:, i) = densities(:, i) + SUM(tempDensity, DIM = 2) * particleList(i)%w_p
        end do
    end subroutine loadParticleDensity

    subroutine WriteParticleDensity(densities, particleList, world, CurrentDiagStep, boolAverage, dirName) 
        ! For diagnostics, deposit single particle density
        ! Re-use rho array since it doesn't get used after first Poisson
        real(real64), intent(in out) :: densities(:,:)
        type(Particle), intent(in) :: particleList(:)
        type(Domain), intent(in) :: world
        integer(int32), intent(in) :: CurrentDiagStep
        character(*), intent(in) :: dirName
        integer(int32) :: i
        logical, intent(in) :: boolAverage
        character(len=5) :: char_i
        
        do i=1, numberChargedParticles
            densities(:, i) = densities(:,i)/world%nodeVol
            write(char_i, '(I3)'), CurrentDiagStep
            if (boolAverage) then
                open(41,file=dirName//'/Density/density_'//particleList(i)%name//"_Average.dat", form='UNFORMATTED')
            else
                open(41,file=dirName//'/Density/density_'//particleList(i)%name//"_"//trim(adjustl(char_i))//".dat", form='UNFORMATTED')
            end if
            write(41) densities(:, i)
            close(41)
        end do
        
    end subroutine WriteParticleDensity
    
end module mod_Scheme