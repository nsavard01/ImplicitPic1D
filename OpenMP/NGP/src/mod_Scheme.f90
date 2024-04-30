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
        particleList(1)%densities(world%threadNodeIndx(1,iThread):world%threadNodeIndx(2,iThread), 1) = SUM(particleList(1)%densities(world%threadNodeIndx(1,iThread):world%threadNodeIndx(2,iThread), :), DIM=2) * particleList(1)%q_times_wp
        do i = 2, numberChargedParticles
            particleList(1)%densities(world%threadNodeIndx(1,iThread):world%threadNodeIndx(2,iThread), 1) = particleList(1)%densities(world%threadNodeIndx(1,iThread):world%threadNodeIndx(2,iThread), 1) &
                + SUM(particleList(i)%densities(world%threadNodeIndx(1,iThread):world%threadNodeIndx(2,iThread), :), DIM=2) * particleList(i)%q_times_wp
        end do
        !$OMP barrier
        !$OMP single
        SELECT CASE (world%boundaryConditions(1))
        CASE(1,4)
            particleList(1)%densities(1, 1) = 0.0d0
        CASE(2)
            particleList(1)%densities(1,1) = 2.0d0 * particleList(1)%densities(1,1)
        CASE(3)
            particleList(1)%densities(1,1) = particleList(1)%densities(1,1) + particleList(1)%densities(NumberXNodes, 1)
        END SELECT
        SELECT CASE (world%boundaryConditions(NumberXNodes))
        CASE(1,4)
            particleList(1)%densities(NumberXNodes, 1) = 0.0d0
        CASE(2)
            particleList(1)%densities(NumberXNodes,1) = 2.0d0 * particleList(1)%densities(NumberXNodes,1)
        CASE(3)
            particleList(1)%densities(NumberXNodes,1) = particleList(1)%densities(1,1)
        END SELECT
        !$OMP end single
        if (world%gridSmoothBool) then
            do i = world%threadNodeIndx(1, iThread), world%threadNodeIndx(2, iThread)
                SELECT CASE(world%boundaryConditions(i))
                CASE(0)
                    rho(i) = 0.25d0 * (particleList(1)%densities(i-1, 1) + 2.0d0 * particleList(1)%densities(i, 1) + particleList(1)%densities(i+1, 1))
                CASE(1,4)
                    rho(i) = 0.0d0
                CASE(2)
                    if (i ==1) then
                        rho(i) = 0.25d0 * (2.0d0 * particleList(1)%densities(1,1) + 2.0d0 * particleList(1)%densities(2,1))
                    else    
                        rho(i) = 0.25d0 * (2.0d0 * particleList(1)%densities(NumberXNodes,1) + 2.0d0 * particleList(1)%densities(NumberXHalfNodes,1))
                    end if
                CASE(3)
                    rho(i) = 0.25d0 * (particleList(1)%densities(NumberXHalfNodes, 1) + 2.0d0 * particleList(1)%densities(1, 1) + particleList(1)%densities(2, 1))
                END SELECT
            end do
        else
            rho(world%threadNodeIndx(1, iThread):world%threadNodeIndx(2, iThread)) = particleList(1)%densities(world%threadNodeIndx(1, iThread):world%threadNodeIndx(2, iThread), 1)
        end if
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
            !$OMP parallel private(iThread)
            iThread = omp_get_thread_num() + 1
            densities(world%threadNodeIndx(1,iThread):world%threadNodeIndx(2,iThread)) = SUM(particleList(i)%densities(world%threadNodeIndx(1,iThread):world%threadNodeIndx(2,iThread), :), DIM=2) * particleList(i)%w_p
            !$OMP end parallel
            if (world%boundaryConditions(1) == 3) then
                densities(1) = (densities(1) + densities(NumberXNodes)) * 0.5d0
                densities(NumberXNodes) = densities(1)
            end if
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