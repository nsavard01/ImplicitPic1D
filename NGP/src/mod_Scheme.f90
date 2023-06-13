module mod_Scheme
    ! module which has NGP specific procedures, typically having to do with the different domain grid structure
    use iso_fortran_env, only: int32, real64
    use constants
    use mod_BasicFunctions
    use mod_particle
    use mod_domain
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
        integer(int32), intent(in out) :: irand
        integer(int32) :: i, numInCell, idxLower, numPerCell(NumberXNodes-1)
        real(real64) :: sumDxDl, L_domain
        idxLower = 1
        sumDxDl = 0
        L_domain = world%grid(NumberXNodes) - world%grid(1)
        do i=1, NumberXNodes-1
            ! Use int to make sure always have a bit left over, otherwise will fill up before getting to end
            numInCell = INT(part%N_p * world%dx_dl(i)/L_domain)
            if (idxLower + numInCell > part % N_P + 1) then
                stop "You are putting too many particles for the uniform particle case"
            end if

            
            call getRandom(part%phaseSpace(1,idxLower:idxLower + numInCell-1), irand)
            part%phaseSpace(1, idxLower:idxLower + numInCell - 1) = part%phaseSpace(1, idxLower:idxLower + numInCell - 1) + i
            idxLower = idxLower + numInCell
            numPerCell(i) = numInCell
            
        end do
        if (idxLower < part%N_p + 1) then
            call getRandom(part%phaseSpace(1, idxLower:part%N_p), irand)
            part%phaseSpace(1, idxLower:part%N_p) = part%phaseSpace(1, idxLower:part%N_p) * (NumberXNodes - 1) + 1
        end if
        
        
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
        l = idxLower + (x - world%grid(idxLower))/world%dx_dl(idxLower)
        
    end function getLFromX

    subroutine addMaxwellianLostParticles(particleList, T_e, T_i, irand, delIdx, idxReFlux, reFluxMaxIdx, x_lower, world)
        ! Add power to all particles in domain
        type(Particle), intent(in out) :: particleList(2)
        type(Domain), intent(in) :: world
        integer(int32), intent(in out) :: irand
        real(real64), intent(in) :: T_e, T_i, x_lower
        integer(int32), intent(in) :: delIdx(2), idxReFlux(1000, 2), reFluxMaxIdx(2)
        integer(int32) :: i,j
        real(real64) :: x_random
        do i=1, delIdx(1)
            call getMaxwellianFluxSample(particleList(1)%phaseSpace(2:4, particleList(1)%N_p + i), particleList(1)%mass, T_e, irand)
            call getMaxwellianFluxSample(particleList(2)%phaseSpace(2:4, particleList(2)%N_p + i), particleList(2)%mass, T_i, irand)
            x_random = world%grid(NumberXNodes) - ran2(irand) * (world%grid(NumberXNodes) - x_lower)
            particleList(1)%phaseSpace(1, particleList(1)%N_p + i) = getLFromX(x_random, world)
            particleList(2)%phaseSpace(1, particleList(2)%N_p + i) = particleList(1)%phaseSpace(1, particleList(1)%N_p + i)
        end do
        particleList(2)%N_p = particleList(2)%N_p + delIdx(1)
        particleList(1)%N_p = particleList(1)%N_p + delIdx(1)
        do j = 1, 2
            do i = 1, reFluxMaxIdx(j)
                if (j == 1) then
                    call getMaxwellianFluxSample(particleList(j)%phaseSpace(2:4, idxReFlux(i, j)), particleList(j)%mass, T_e, irand)
                else
                    call getMaxwellianFluxSample(particleList(j)%phaseSpace(2:4, idxReFlux(i, j)), particleList(j)%mass, T_i, irand)
                end if
                particleList(j)%phaseSpace(2:4, idxReFlux(i, j)) = -ABS(particleList(j)%phaseSpace(2:4, idxReFlux(i, j)))
            end do
        end do

    end subroutine addMaxwellianLostParticles



    ! --------------------------- Diagnostics ------------------------------------

    subroutine depositRho(rho, particleList, world) 
        real(real64), intent(in out) :: rho(NumberXNodes)
        type(Particle), intent(in) :: particleList(:)
        type(Domain), intent(in) :: world
        integer(int32) :: i, j, l_left, l_right
        real(real64) :: d
        rho = 0.0d0
        do i=1, numberChargedParticles
            do j = 1, particleList(i)%N_p
                l_left = INT(particleList(i)%phaseSpace(1, j))
                d = MOD(particleList(i)%phaseSpace(1, j), 1.0d0)
                rho(l_left) = rho(l_left) + particleList(i)%q * particleList(i)%w_p * (1.0d0-d)
                rho(l_left + 1) = rho(l_left + 1) + particleList(i)%q * particleList(i)%w_p * d
            end do
        end do
        rho = rho / world%nodeVol
        if (world%boundaryConditions(1) == 3) then
            rho(1) = rho(1) + rho(NumberXNodes)
            rho(NumberXNodes) = rho(1)
        else if (world%boundaryConditions(1) == 2) then
            rho(1) = rho(1)*2.0d0
        end if

        if (world%boundaryConditions(NumberXNodes) == 2) rho(NumberXNodes) = rho(NumberXNodes)*2.0d0
    end subroutine depositRho

    subroutine loadParticleDensity(densities, particleList, world)
        type(Particle), intent(in) :: particleList(:)
        type(Domain), intent(in) :: world
        real(real64), intent(in out) :: densities(:,:)
        integer(int32) :: i,j, l_left
        real(real64) :: d
        do i=1, numberChargedParticles
            do j = 1, particleList(i)%N_p
                l_left = INT(particleList(i)%phaseSpace(1,j))
                d = MOD(particleList(i)%phaseSpace(1,j), 1.0d0)
                densities(l_left, i) = densities(l_left, i) + particleList(i)%w_p * (1.0d0-d)
                densities(l_left + 1, i) = densities(l_left + 1, i) + particleList(i)%w_p * d
            end do
        end do
        if (world%boundaryConditions(l_left) > 4) then
            print *, 'Just excuse to use world object since needed in CIC scheme'
        end if

    end subroutine loadParticleDensity

    subroutine WriteParticleDensity(densities, particleList, world, CurrentDiagStep, boolAverage) 
        ! For diagnostics, deposit single particle density
        ! Re-use rho array since it doesn't get used after first Poisson
        real(real64), intent(in out) :: densities(:,:)
        type(Particle), intent(in) :: particleList(:)
        type(Domain), intent(in) :: world
        integer(int32), intent(in) :: CurrentDiagStep
        integer(int32) :: i
        logical, intent(in) :: boolAverage
        character(len=5) :: char_i
        
        do i=1, numberChargedParticles
            densities(:,i) = densities(:,i)/world%nodeVol
            if (world%boundaryConditions(1) == 3) then
                densities(1,i) = densities(1,i) + densities(NumberXNodes,i)
                densities(NumberXNodes,i) = densities(1, i)
            else if (world%boundaryConditions(1) == 2) then
                densities(1,i) = densities(1,i)*2.0d0
            end if
    
            if (world%boundaryConditions(NumberXNodes) == 2) densities(NumberXNodes, i) = densities(NumberXNodes, i)*2.0d0
            write(char_i, '(I3)'), CurrentDiagStep
            if (boolAverage) then
                open(41,file='../Data/Density/density_'//particleList(i)%name//"_Average.dat", form='UNFORMATTED')
            else
                open(41,file='../Data/Density/density_'//particleList(i)%name//"_"//trim(adjustl(char_i))//".dat", form='UNFORMATTED')
            end if
            write(41) densities(:,i)
            close(41)
        end do
        
    end subroutine WriteParticleDensity

    !----------------------- Subroutines for specific examples-----------------------
    

end module mod_Scheme