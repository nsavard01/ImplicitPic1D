module mod_Scheme
    ! module which has NGP specific procedures, typically having to do with the different domain grid structure
    use iso_fortran_env, only: int32, real64
    use constants
    use mod_BasicFunctions
    use mod_particle
    use mod_domain
    use omp_lib
    implicit none
    integer(int32), protected ::schemeNum

contains

    subroutine initializeScheme()
        schemeNum = 0
        print *, "--Scheme---"
        print *, "NGP constant grid size between phi/rho nodes"
        print *, '----------'
    end subroutine initializeScheme

    ! subroutine initialize_randUniform(part, world, irand)
    !     ! place particles randomly in each dx_dl based on portion of volume it take up
    !     type(Particle), intent(in out) :: part
    !     type(Domain), intent(in) :: world
    !     integer(int32), intent(in out) :: irand(numThread)
    !     integer(int32) :: i, iThread
    !     real(real64) :: L_domain
    !     L_domain = world%grid(NumberXNodes) - world%grid(1)
    !     !$OMP parallel private(iThread, i)
    !     iThread = omp_get_thread_num() + 1
    !     do i=1, part%N_p(iThread)
    !         part%phaseSpace(1,i, iThread) = ran2(irand(iThread)) * L_domain + world%grid(1)
    !         part%phaseSpace(1,i, iThread) = world%getLFromX(part%phaseSpace(1,i, iThread))
    !     end do
    !     !$OMP end parallel    
    ! end subroutine initialize_randUniform

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

    subroutine interpolateParticleToNodes(part, world, iThread)
        ! Interpolate particles to logical grid for a single thread
        type(Particle), intent(in out) :: part
        type(Domain), intent(in) :: world
        integer(int32), intent(in) :: iThread
        integer(int32) :: j, l_left, l_right
        real(real64) :: d

        do j = 1, part%N_p(iThread)
            l_left = INT(part%phaseSpace(1, j, iThread))
            l_right = l_left + 1
            d = part%phaseSpace(1, j, iThread) - real(l_left)
            part%densities(l_left, iThread) = part%densities(l_left, iThread) + (1.0d0-d)
            part%densities(l_right, iThread) = part%densities(l_right, iThread) + d
        end do
    end subroutine interpolateParticleToNodes


    subroutine loadParticleDensity(particleList, world, reset)
        type(Particle), intent(in out) :: particleList(numberChargedParticles)
        type(Domain), intent(in) :: world
        logical, intent(in) :: reset
        integer(int32) :: i, iThread
        !$OMP parallel private(iThread, i)
        do i = 1, numberChargedParticles  
            iThread = omp_get_thread_num() + 1 
            if (reset) then
                particleList(i)%densities(:,iThread) = 0.0d0
            end if
            call interpolateParticleToNodes(particleList(i), world, iThread)
        end do
        !$OMP end parallel
    end subroutine

    subroutine depositRho(rho, particleList, world) 
        real(real64), intent(in out) :: rho(NumberXNodes)
        type(Particle), intent(in out) :: particleList(numberChargedParticles)
        type(Domain), intent(in) :: world
        integer(int32) :: i, iThread
        rho = 0.0d0
        !$OMP parallel private(iThread, i)
        do i = 1, numberChargedParticles
            iThread = omp_get_thread_num() + 1 
            particleList(i)%densities(:,iThread) = 0.0d0
            call interpolateParticleToNodes(particleList(i), world, iThread) 
        end do
        !$OMP barrier
        do i = 1, numberChargedParticles
            call world%addThreadedDomainArray(rho, particleList(i)%densities, NumberXNodes, iThread, particleList(i)%q_times_wp)
        end do
        !$OMP end parallel  
    end subroutine depositRho


    subroutine WriteParticleDensity(particleList, world, CurrentDiagStep, boolAverage, dirName) 
        ! For diagnostics, deposit single particle density
        ! Re-use rho array since it doesn't get used after first Poisson
        type(Particle), intent(in) :: particleList(numberChargedParticles)
        type(Domain), intent(in) :: world
        integer(int32), intent(in) :: CurrentDiagStep
        character(*), intent(in) :: dirName
        logical, intent(in) :: boolAverage
        real(real64) :: densities(NumberXNodes)
        integer(int32) :: i, iThread
        character(len=5) :: char_i
        
        do i=1, numberChargedParticles
            densities = 0.0d0
            !$OMP parallel private(iThread)
            iThread = omp_get_thread_num() + 1
            call world%addThreadedDomainArray(densities, particleList(i)%densities, NumberXNodes, iThread, particleList(i)%w_p)
            !$OMP end parallel
            densities = densities/world%nodeVol
            write(char_i, '(I4)'), CurrentDiagStep
            if (boolAverage) then
                open(41,file=dirName//'/Density/density_'//particleList(i)%name//"_Average.dat", form='UNFORMATTED')
            else
                open(41,file=dirName//'/Density/density_'//particleList(i)%name//"_"//trim(adjustl(char_i))//".dat", form='UNFORMATTED')
            end if
            write(41) densities
            close(41)
        end do
        
    end subroutine WriteParticleDensity
    
end module mod_Scheme