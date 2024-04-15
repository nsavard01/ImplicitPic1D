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
            if (world%boundaryConditions(l_center) == 0 .and. world%boundaryConditions(l_right) == 0) then
                ! No Boundary either side
                part%densities(l_left, iThread) = part%densities(l_left, iThread) + 0.5d0 * (1.0d0 - d)**2
                part%densities(l_right, iThread) = part%densities(l_right, iThread) + 0.5d0 * d**2
            else if (world%boundaryConditions(l_right) == 1 .or. world%boundaryConditions(l_right) == 4) then
                ! Dirichlet to right
                part%densities(l_center, iThread) = part%densities(l_center, iThread) - 0.5d0 * d**2
                part%densities(l_left, iThread) = part%densities(l_left, iThread) + 0.5d0 * (1.0d0 - d)**2
            else if (world%boundaryConditions(l_center) == 1 .or. world%boundaryConditions(l_center) == 4) then
                !Dirichlet to left
                part%densities(l_center, iThread) = part%densities(l_center, iThread) - 0.5d0 * (1.0d0 - d)**2
                part%densities(l_right, iThread) = part%densities(l_right, iThread) + 0.5d0 * d**2
            else if (world%boundaryConditions(l_right) == 2) then
                !Neumann to right
                part%densities(l_center, iThread) = part%densities(l_center, iThread) + 0.5d0 * d**2
                part%densities(l_left, iThread) = part%densities(l_left, iThread) + 0.5d0 * (1.0d0 - d)**2
            else if (world%boundaryConditions(l_center) == 2) then
                !Neumann to left
                part%densities(l_center, iThread) = part%densities(l_center, iThread) + 0.5d0 * (1.0d0 - d)**2
                part%densities(l_right, iThread) = part%densities(l_right, iThread) + 0.5d0 * d**2
            end if
        end do
    end subroutine interpolateParticleToNodes

    subroutine loadParticleDensity(particleList, world, reset)
        type(Particle), intent(in out) :: particleList(:)
        type(Domain), intent(in) :: world
        logical, intent(in) :: reset
        integer(int32) :: i, iThread
        do i = 1, numberChargedParticles
            !$OMP parallel private(iThread)
            iThread = omp_get_thread_num() + 1 
            if (reset) then
                particleList(i)%densities(:,iThread) = 0.0d0
            end if
            call interpolateParticleToNodes(particleList(i), world, iThread)
            !$OMP end parallel
        end do
    end subroutine

    subroutine depositRho(rho, particleList, world) 
        real(real64), intent(in out) :: rho(NumberXNodes)
        type(Particle), intent(in out) :: particleList(:)
        type(Domain), intent(in) :: world
        integer(int32) :: i, iThread
        rho = 0.0d0
        do i = 1, numberChargedParticles
            !$OMP parallel private(iThread)
            iThread = omp_get_thread_num() + 1 
            particleList(i)%densities(:,iThread) = 0.0d0
            call interpolateParticleToNodes(particleList(i), world, iThread)
            !$OMP end parallel   
            rho = rho + SUM(particleList(i)%densities, DIM = 2) * particleList(i)%w_p * particleList(i)%q
        end do
    end subroutine depositRho

    subroutine WriteParticleDensity(particleList, world, CurrentDiagStep, boolAverage, dirName) 
        ! For diagnostics, deposit single particle density
        ! Re-use rho array since it doesn't get used after first Poisson
        type(Particle), intent(in) :: particleList(:)
        type(Domain), intent(in) :: world
        integer(int32), intent(in) :: CurrentDiagStep
        character(*), intent(in) :: dirName
        real(real64) :: densities(NumberXNodes)
        integer(int32) :: i
        logical, intent(in) :: boolAverage
        character(len=5) :: char_i
        
        do i=1, numberChargedParticles
            densities = (SUM(particleList(i)%densities, DIM=2)/world%dx_dl) * particleList(i)%w_p
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