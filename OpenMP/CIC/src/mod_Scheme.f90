module mod_Scheme
    use iso_fortran_env, only: int32, real64
    use constants
    use mod_BasicFunctions
    use mod_particle
    use mod_domain
    use omp_lib
    implicit none
    ! Scheme module for CIC
contains

    subroutine initializeScheme(schemeNum)
        integer(int32), intent(in out) :: schemeNum
        schemeNum = 1
        print *, "--Scheme---"
        print *, "CIC constant grid size between half nodes"
        print *, '----------'
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
            part%phaseSpace(1,i, iThread) = getLFromX(part%phaseSpace(1,i, iThread), world)
        end do
        !$OMP end parallel    
    end subroutine initialize_randUniform

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
        if (x >= world%grid(idxLower) + 0.5d0 * world%nodeVol(idxLower)) then
            l = idxHigher + (x - world%grid(idxHigher))/world%nodeVol(idxHigher)
        else
            l = idxLower + (x - world%grid(idxLower))/world%nodeVol(idxLower)
        end if   
    end function getLFromX



    subroutine depositRho(rho, particleList, world) 
        real(real64), intent(in out) :: rho(NumberXNodes)
        type(Particle), intent(in) :: particleList(:)
        type(Domain), intent(in) :: world
        integer(int32) :: i, j, l_center, l_left, l_right, iThread
        real(real64) :: d, rhoTemp(NumberXNodes, numThread)
        rho = 0.0d0
        do i=1, numberChargedParticles
            rhoTemp = 0.0d0
            !$OMP parallel private(iThread, j, l_center, l_left, l_right, d)
            iThread = omp_get_thread_num() + 1
            do j = 1, particleList(i)%N_p(iThread)
                l_center = NINT(particleList(i)%phaseSpace(1, j, iThread))
                d = particleList(i)%phaseSpace(1, j, iThread) - real(l_center)
                SELECT CASE (world%boundaryConditions(l_center))
                CASE(0)
                    ! Inside domain
                    l_left = l_center -1
                    l_right = l_center + 1
                    rhoTemp(l_center, iThread) = rhoTemp(l_center, iThread) + (0.75 - d**2)
                    rhoTemp(l_right, iThread) = rhoTemp(l_right, iThread) + 0.5d0 * (0.5d0 + d)**2
                    rhoTemp(l_left, iThread) = rhoTemp(l_left, iThread) + 0.5d0 * (0.5d0 - d)**2
                    ! SELECT CASE (world%boundaryConditions(l_right))
                    ! CASE(0,1)
                    !     rhoTemp(l_right, iThread) = rhoTemp(l_right, iThread) + 0.5d0 * (0.5d0 + d)**2   
                    ! CASE(2)
                    !     rhoTemp(l_right, iThread) = rhoTemp(l_right, iThread) + 0.5d0 * (0.5d0 + d)**2
                    ! CASE(3)
                    !     rhoTemp(l_right, iThread) = rhoTemp(l_right, iThread) + 0.5d0 * (0.5d0 + d)**2
                    !     rhoTemp(1, iThread) = rhoTemp(1, iThread) + 0.5d0 * (0.5d0 + d)**2
                    ! END SELECT
                    
                    ! SELECT CASE (world%boundaryConditions(l_left))
                    ! CASE(0,1)
                    !     rhoTemp(l_left, iThread) = rhoTemp(l_left, iThread) + 0.5d0 * (0.5d0 - d)**2
                    ! CASE(2)
                    !     rhoTemp(l_left, iThread) = rhoTemp(l_left, iThread) + 0.5d0 * (0.5d0 - d)**2
                    ! CASE(3)
                    !     rhoTemp(l_left, iThread) = rhoTemp(l_left, iThread) + 0.5d0 * (0.5d0 - d)**2
                    !     rhoTemp(NumberXNodes, iThread) = rhoTemp(NumberXNodes, iThread) + 0.5d0 * (0.5d0 - d)**2
                    ! END SELECT
                CASE(1)
                    !Dirichlet
                    rhoTemp(l_center, iThread) = rhoTemp(l_center, iThread) + (1.0d0-ABS(d))
                    rhoTemp(l_center + INT(SIGN(1.0, d)), iThread) = rhoTemp(l_center + INT(SIGN(1.0, d)), iThread) + ABS(d)
                CASE(2)
                    !Neumann symmetric
                    rhoTemp(l_center, iThread) = rhoTemp(l_center, iThread) + (0.75 - d**2)
                    rhoTemp(l_center + INT(SIGN(1.0, d)), iThread) = rhoTemp(l_center + INT(SIGN(1.0, d)), iThread) + (0.25d0 + d**2)
                CASE(3)
                    ! Periodic
                    rhoTemp(l_center, iThread) = rhoTemp(l_center, iThread) + (0.75 - d**2)
                    rhoTemp(ABS(l_center - NumberXNodes) + 1, iThread) = rhoTemp(ABS(l_center - NumberXNodes) + 1, iThread) + (0.75 - d**2)
                    ! towards domain
                    rhoTemp(l_center+INT(SIGN(1.0, d)), iThread) = rhoTemp(l_center+INT(SIGN(1.0, d)), iThread) + 0.5d0 * (0.5d0 + ABS(d))**2
                    ! across periodic boundary
                    rhoTemp(MODULO(l_center-2*INT(SIGN(1.0, d)),NumberXNodes), iThread) = rhoTemp(MODULO(l_center-2*INT(SIGN(1.0, d)),NumberXNodes), iThread) + 0.5d0 * (0.5d0 - ABS(d))**2
                END SELECT
            end do
            !$OMP end parallel
            rho = rho + SUM(rhoTemp, DIM = 2) * particleList(i)%w_p * particleList(i)%q
        end do
    end subroutine depositRho

    subroutine loadParticleDensity(densities, particleList, world)
        type(Particle), intent(in) :: particleList(numberChargedParticles)
        type(Domain), intent(in) :: world
        real(real64), intent(in out) :: densities(NumberXNodes,numberChargedParticles)
        integer(int32) :: i,j, l_left, l_right, l_center, iThread
        real(real64) :: d, tempDensity(NumberXNodes, numThread)
        do i=1, numberChargedParticles
            tempDensity = 0.0d0
            !$OMP parallel private(iThread, j, l_center, l_left, l_right, d)
            iThread = omp_get_thread_num() + 1
            do j = 1, particleList(i)%N_p(iThread)
                l_center = NINT(particleList(i)%phaseSpace(1, j, iThread))
                d = particleList(i)%phaseSpace(1, j, iThread) - real(l_center)
                SELECT CASE (world%boundaryConditions(l_center))
                CASE(0)
                    ! Inside domain
                    l_left = l_center -1
                    l_right = l_center + 1
                    tempDensity(l_center, iThread) = tempDensity(l_center, iThread) + (0.75 - d**2)
                    tempDensity(l_right, iThread) = tempDensity(l_right, iThread) + 0.5d0 * (0.5d0 + d)**2
                    tempDensity(l_left, iThread) = tempDensity(l_left, iThread) + 0.5d0 * (0.5d0 - d)**2
                    ! SELECT CASE (world%boundaryConditions(l_right))
                    ! CASE(0,1)
                    !     rhoTemp(l_right) = rhoTemp(l_right) + 0.5d0 * (0.5d0 + d)**2   
                    ! CASE(1)
                    !     rhoTemp(l_right) = rhoTemp(l_right) + 0.5d0 * (0.5d0 + d)**2
                    ! CASE(2)
                    !     rhoTemp(l_right) = rho(l_right) + particleList(i)%q * particleList(i)%w_p * (0.5d0 + d)**2
                    ! CASE(3)
                    !     rho(l_right) = rho(l_right) + particleList(i)%q * particleList(i)%w_p * 0.5d0 * (0.5d0 + d)**2
                    !     rho(1) = rho(1) + particleList(i)%q * particleList(i)%w_p * 0.5d0 * (0.5d0 + d)**2
                    ! END SELECT
                    ! if (l_left /= 1) then
                    !     rho(l_left) = rho(l_left) + particleList(i)%q * particleList(i)%w_p * 0.5d0 * (0.5d0 - d)**2
                    ! else
                    !     SELECT CASE (world%boundaryConditions(l_left))
                    !     CASE(1)
                    !         rho(l_left) = rho(l_left) + particleList(i)%q * particleList(i)%w_p * 0.5d0 * (0.5d0 - d)**2
                    !     CASE(2)
                    !         rho(l_left) = rho(l_left) + particleList(i)%q * particleList(i)%w_p * (0.5d0 - d)**2
                    !     CASE(3)
                    !         rho(l_left) = rho(l_left) + particleList(i)%q * particleList(i)%w_p * 0.5d0 * (0.5d0 - d)**2
                    !         rho(NumberXNodes) = rho(NumberXNodes) + particleList(i)%q * particleList(i)%w_p * 0.5d0 * (0.5d0 - d)**2
                    !     END SELECT
                    ! end if
                CASE(1)
                    !Dirichlet
                    tempDensity(l_center, iThread) = tempDensity(l_center, iThread) + (1.0d0-ABS(d))
                    tempDensity(l_center + INT(SIGN(1.0, d)), iThread) = tempDensity(l_center + INT(SIGN(1.0, d)), iThread) + ABS(d)
                CASE(2)
                    !Neumann symmetric
                    tempDensity(l_center, iThread) = tempDensity(l_center, iThread) + (0.75 - d**2)
                    tempDensity(l_center + INT(SIGN(1.0, d)), iThread) = tempDensity(l_center + INT(SIGN(1.0, d)), iThread) + (0.25d0 + d**2)
                CASE(3)
                    ! Periodic
                    tempDensity(l_center, iThread) = tempDensity(l_center, iThread) + (0.75 - d**2)
                    tempDensity(ABS(l_center - NumberXNodes) + 1, iThread) = tempDensity(ABS(l_center - NumberXNodes) + 1, iThread) + (0.75 - d**2)
                    ! towards domain
                    tempDensity(l_center+INT(SIGN(1.0, d)), iThread) = tempDensity(l_center+INT(SIGN(1.0, d)), iThread) + 0.5d0 * (0.5d0 + ABS(d))**2
                    ! across periodic boundary
                    tempDensity(MODULO(l_center-2*INT(SIGN(1.0, d)),NumberXNodes), iThread) = tempDensity(MODULO(l_center-2*INT(SIGN(1.0, d)),NumberXNodes), iThread) + 0.5d0 * (0.5d0 - ABS(d))**2
                END SELECT
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
        integer(int32) :: i, j
        logical, intent(in) :: boolAverage
        character(len=5) :: char_i
        
        do i=1, numberChargedParticles
            do j = 1, NumberXNodes
                if (world%boundaryConditions(j) == 0) then
                    densities(j, i) = densities(j,i)/world%nodeVol(j)
                else
                    densities(j, i) = 2.0d0 * densities(j,i)/world%nodeVol(j)
                end if
            end do
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