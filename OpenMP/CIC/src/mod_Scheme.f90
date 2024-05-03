module mod_Scheme
    use iso_fortran_env, only: int32, real64
    use constants
    use mod_BasicFunctions
    use mod_particle
    use mod_domain
    use omp_lib
    implicit none
    integer(int32), protected ::schemeNum
    ! Scheme module for CIC
contains

    subroutine initializeScheme()
        schemeNum = 1
        print *, "--Scheme---"
        print *, "CIC constant grid size between Field/J nodes"
        print *, '----------'
    end subroutine initializeScheme


    subroutine interpolateParticleToNodes(part, world, iThread)
        ! Interpolate particles to logical grid for a single thread
        type(Particle), intent(in out) :: part
        type(Domain), intent(in) :: world
        integer(int32), intent(in) :: iThread
        integer(int32) :: j, l_center, l_left, l_right
        real(real64) :: d

        do j = 1, part%N_p(iThread)
            l_center = INT(part%phaseSpace(1, j, iThread))
            d = part%phaseSpace(1, j, iThread) - real(l_center)
            part%densities(l_center, iThread) = part%densities(l_center, iThread) + (-d**2 + d + 0.5d0)
            l_right = l_center + 1
            l_left = l_center - 1
            SELECT CASE (world%boundaryConditions(l_right) - world%boundaryConditions(l_center))
            CASE(0)
                ! No Boundary either side
                part%densities(l_left, iThread) = part%densities(l_left, iThread) + 0.5d0 * (1.0d0 - d)**2
                part%densities(l_right, iThread) = part%densities(l_right, iThread) + 0.5d0 * d**2
            CASE(-1,-4)
                !Dirichlet to left
                part%densities(l_center, iThread) = part%densities(l_center, iThread) - 0.5d0 * (1.0d0 - d)**2
                part%densities(l_right, iThread) = part%densities(l_right, iThread) + 0.5d0 * d**2
            CASE(1,4)
                ! Dirichlet to right
                part%densities(l_center, iThread) = part%densities(l_center, iThread) - 0.5d0 * d**2
                part%densities(l_left, iThread) = part%densities(l_left, iThread) + 0.5d0 * (1.0d0 - d)**2
            CASE(-2)
                !Neumann to left
                part%densities(l_center, iThread) = part%densities(l_center, iThread) + 0.5d0 * (1.0d0 - d)**2
                part%densities(l_right, iThread) = part%densities(l_right, iThread) + 0.5d0 * d**2
            CASE(2)
                !Neumann to right
                part%densities(l_center, iThread) = part%densities(l_center, iThread) + 0.5d0 * d**2
                part%densities(l_left, iThread) = part%densities(l_left, iThread) + 0.5d0 * (1.0d0 - d)**2
            CASE(-3)
                !periodic to left
                part%densities(NumberXNodes, iThread) = part%densities(NumberXNodes, iThread) + 0.5d0 * (1.0d0 - d)**2
                part%densities(l_right, iThread) = part%densities(l_right, iThread) + 0.5d0 * d**2
            CASE(3)
                !periodic to right
                part%densities(l_left, iThread) = part%densities(l_left, iThread) + 0.5d0 * (1.0d0 - d)**2
                part%densities(1, iThread) = part%densities(1, iThread) + 0.5d0 * d**2
            END SELECT
        end do
    end subroutine interpolateParticleToNodes

    subroutine loadParticleDensity(particleList, world, reset)
        type(Particle), intent(in out) :: particleList(:)
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
        ! Re-use rho array since it doesn't get used after first Poisson
        type(Particle), intent(in out) :: particleList(:)
        type(Domain), intent(in) :: world
        integer(int32), intent(in) :: CurrentDiagStep
        character(*), intent(in) :: dirName
        real(real64) :: densities(NumberXNodes)
        integer(int32) :: i, iThread, j, leftThreadIndx, rightThreadIndx
        logical, intent(in) :: boolAverage
        character(len=5) :: char_i
        
        do i=1, numberChargedParticles
            densities = 0.0d0
            !$OMP parallel private(iThread, j, leftThreadIndx, rightThreadIndx)
            iThread = omp_get_thread_num() + 1
            leftThreadIndx = world%threadNodeIndx(1,iThread)
            rightThreadIndx = world%threadNodeIndx(2,iThread)
            particleList(i)%densities(leftThreadIndx:rightThreadIndx, 1) = &
                SUM(particleList(i)%densities(leftThreadIndx:rightThreadIndx, :), DIM=2) * particleList(i)%w_p
            !$OMP barrier
            if (world%gridSmoothBool) then
                do j = leftThreadIndx, rightThreadIndx
                    SELECT CASE(world%boundaryConditions(j+1) - world%boundaryConditions(j))
                    CASE(0)
                        densities(j) = 0.25d0 * (particleList(i)%densities(j-1, 1) + 2.0d0 * particleList(i)%densities(j, 1) + particleList(i)%densities(j+1, 1))
                    CASE(-1,-4)
                        !Dirichlet to left
                        densities(j) = 0.25d0 * (2.0d0 * particleList(i)%densities(1, 1) + particleList(i)%densities(2, 1))
                    CASE(1,4)
                        ! Dirichlet to right
                        densities(j) = 0.25d0 * (2.0d0 * particleList(i)%densities(NumberXNodes, 1) + particleList(i)%densities(NumberXNodes-1, 1))
                    CASE(-2)
                        !Neumann to left
                        densities(j) = 0.25d0 * (particleList(i)%densities(2, 1) + 3.0d0 * particleList(i)%densities(1, 1))
                    CASE(2)
                        !Neumann to right
                        densities(j) = 0.25d0 * (particleList(i)%densities(NumberXNodes-1, 1) + 3.0d0 * particleList(i)%densities(NumberXNodes, 1))
                    CASE(-3)
                        !periodic to left
                        densities(j) = 0.25d0 * (particleList(i)%densities(NumberXNodes, 1) + 2.0d0 * particleList(i)%densities(1, 1) + particleList(i)%densities(2, 1))
                    CASE(3)
                        !periodic to right
                        densities(j) = 0.25d0 * (particleList(i)%densities(1, 1) + 2.0d0 * particleList(i)%densities(NumberXNodes, 1) + particleList(i)%densities(NumberXNodes-1, 1))
                    END SELECT
                end do
            else
                densities(leftThreadIndx:rightThreadIndx) = particleList(i)%densities(leftThreadIndx:rightThreadIndx, 1)
            end if
            densities(leftThreadIndx:rightThreadIndx) = densities(leftThreadIndx:rightThreadIndx)/world%dx_dl(leftThreadIndx:rightThreadIndx)   
            !$OMP end parallel
            write(char_i, '(I3)'), CurrentDiagStep
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