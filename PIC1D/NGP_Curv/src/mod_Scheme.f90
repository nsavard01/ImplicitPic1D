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
        ! Initialize scheme number
        schemeNum = 0
        print *, "--Scheme---"
        print *, "NGP constant grid size between phi/rho nodes"
        print *, '----------'
    end subroutine initializeScheme


    subroutine initialize_QuasiNeutral(ions, electrons)
        ! Initialze plasma where ions and electrons start at same location
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
        ! load particles to density array for each thread
        type(Particle), intent(in out) :: particleList(numberChargedParticles)
        type(Domain), intent(in) :: world
        logical, intent(in) :: reset
        integer(int32) :: i, iThread
        !$OMP parallel private(iThread, i)
        iThread = omp_get_thread_num() + 1 
        do i = 1, numberChargedParticles  
            if (reset) then
                particleList(i)%densities(:,iThread) = 0.0d0
            end if
            call interpolateParticleToNodes(particleList(i), world, iThread)
        end do
        !$OMP end parallel
    end subroutine


    subroutine WriteParticleDensity(particleList, world, CurrentDiagStep, boolAverage, dirName) 
        ! For diagnostics, deposit single particle density
        ! Can choose to have averaged as well
        type(Particle), intent(in out) :: particleList(numberChargedParticles)
        type(Domain), intent(in) :: world
        integer(int32), intent(in) :: CurrentDiagStep
        character(*), intent(in) :: dirName
        logical, intent(in) :: boolAverage
        real(real64) :: densities(NumberXNodes)
        integer(int32) :: i, iThread, j, leftThreadIndx, rightThreadIndx
        character(len=5) :: char_i
        
        do i=1, numberChargedParticles
            ! accumulate particle densities on different threads
            !$OMP parallel private(iThread, j, leftThreadIndx, rightThreadIndx)
            iThread = omp_get_thread_num() + 1
            leftThreadIndx = world%threadNodeIndx(1,iThread)
            rightThreadIndx = world%threadNodeIndx(2,iThread)
            particleList(i)%densities(leftThreadIndx:rightThreadIndx, 1) = SUM(particleList(i)%densities(leftThreadIndx:rightThreadIndx, :), DIM=2) * particleList(i)%w_p
            densities(leftThreadIndx:rightThreadIndx) = 0.0d0
            !$OMP end parallel
            if (world%boundaryConditions(1) == 3) then
                particleList(i)%densities(1,1) = (particleList(i)%densities(1,1) + particleList(i)%densities(NumberXNodes,1))
                particleList(i)%densities(NumberXNodes,1) = particleList(i)%densities(1,1)
            end if
            if (world%gridSmoothBool) then
                do j = 1, NumberXNodes
                    SELECT CASE(world%boundaryConditions(j))
                    CASE(0)
                        densities(j) = densities(j) + 0.25d0 * (particleList(i)%densities(j-1, 1) + 2.0d0 * particleList(i)%densities(j, 1) + particleList(i)%densities(j+1, 1))
                    CASE(1,4,2)
                        if (j == 1) then
                            densities(j) = densities(j) + 0.25d0 * (2.0d0 * particleList(i)%densities(1,1) + particleList(i)%densities(2,1))
                            densities(j+1) = densities(j+1) + 0.25d0 * particleList(i)%densities(1,1)
                        else
                            densities(j) = densities(j) + 0.25d0 * (2.0d0 * particleList(i)%densities(NumberXNodes,1) + particleList(i)%densities(NumberXHalfNodes,1))
                            densities(j-1) = densities(j-1) + 0.25d0 * particleList(i)%densities(NumberXNodes,1)
                        end if
                    CASE(3)
                        densities(j) = densities(j) + 0.25d0 * (particleList(i)%densities(NumberXHalfNodes, 1) + 2.0d0 * particleList(i)%densities(1, 1) + particleList(i)%densities(2, 1))
                    END SELECT
                end do
            else
                densities = particleList(i)%densities(:,1)
            end if
            SELECT CASE (world%boundaryConditions(1))
            CASE(1,2,4)
                densities(1) = 2.0d0 * densities(1)/world%nodeVol(1)
            CASE(3)
                densities(1) = 2.0d0 * densities(1)/(world%nodeVol(1) + world%nodeVol(NumberXNodes))
            END SELECT 
            do j = 2, NumberXHalfNodes
                densities(j) = densities(j) / (world%nodeVol(j))
            end do
            SELECT CASE (world%boundaryConditions(NumberXNodes))
            CASE(1,2,4)
                densities(NumberXNodes) = 2.0d0 * densities(NumberXNodes)/world%nodeVol(NumberXNodes)
            CASE(3)
                densities(NumberXNodes) = densities(1)
            END SELECT

            ! Write density
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