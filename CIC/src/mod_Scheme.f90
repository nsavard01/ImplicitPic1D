module mod_Scheme
    use iso_fortran_env, only: int32, real64
    use constants
    use mod_BasicFunctions
    use mod_particle
    use mod_domain
    implicit none
    ! Scheme module for CIC
contains

    subroutine initializeScheme(schemeNum)
        integer(int32), intent(in out) :: schemeNum
        schemeNum = 1
    end subroutine initializeScheme


    subroutine initialize_randUniform(part, world, irand)
        ! place particles randomly in each dx_dl based on portion of volume it take up
        ! 
        type(Particle), intent(in out) :: part
        type(Domain), intent(in) :: world
        integer(int32), intent(in out) :: irand
        integer(int32) :: i, numInCell, idxLower, numPerCell(NumberXNodes)
        real(real64) :: L_domain
        idxLower = 1
        L_domain = world%grid(NumberXNodes) - world%grid(1)
        do i=1, NumberXNodes
            ! Use int to make sure always have a bit left over, otherwise will fill up before getting to end
            if (world%boundaryConditions(i) /= 0) then
                ! if near boundary, cell volume divided by two
                numInCell = INT(part%N_p * world%nodeVol(i)/L_domain/2.0d0)
                call getRandom(part%phaseSpace(1,idxLower:idxLower + numInCell-1), irand)
                numPerCell(i) = numInCell
                if (i == 1) then
                    part%phaseSpace(1, idxLower:idxLower + numInCell - 1) = part%phaseSpace(1, idxLower:idxLower + numInCell - 1) * 0.5d0 + 1.0d0
                else
                    part%phaseSpace(1, idxLower:idxLower + numInCell - 1) = NumberXNodes - part%phaseSpace(1, idxLower:idxLower + numInCell - 1) * 0.5d0
                end if
            else
                numInCell = INT(part%N_p * world%nodeVol(i)/L_domain)
                call getRandom(part%phaseSpace(1,idxLower:idxLower + numInCell-1), irand)
                part%phaseSpace(1, idxLower:idxLower + numInCell - 1) = part%phaseSpace(1, idxLower:idxLower + numInCell - 1) + i - 0.5d0
                numPerCell(i) = numInCell
            end if
            if (idxLower + numInCell> part % N_P + 1) then
                print *, "numInCell is:", numInCell
                print *, "on", i, "cell"
                stop "You are putting too many particles for the uniform particle case"
            end if  
            if (numInCell < 1) then
                print *, "You are placing 0 particles in a cell."
                print *, "Either the cell size is very small with respect to the domain, or too few particles."
            end if
            idxLower = idxLower + numInCell
            
        end do
        if (idxLower < part%N_p + 1) then
            call getRandom(part%phaseSpace(1, idxLower:part%N_p), irand)
            part%phaseSpace(1, idxLower:part%N_p) = part%phaseSpace(1, idxLower:part%N_p) * (NumberXNodes-1) + 1
        end if
        
    end subroutine initialize_randUniform



    subroutine depositRho(rho, particleList, world) 
        real(real64), intent(in out) :: rho(NumberXNodes)
        type(Particle), intent(in) :: particleList(:)
        type(Domain), intent(in) :: world
        integer(int32) :: i, j, l_center, l_left, l_right
        real(real64) :: d
        do i=1, numberChargedParticles
            do j = 1, particleList(i)%N_p
                l_center = NINT(particleList(i)%phaseSpace(1, j))
                d = particleList(i)%phaseSpace(1, j) - l_center
                l_left = l_center -1
                l_right = l_center + 1
                SELECT CASE (world%boundaryConditions(l_center))
                CASE(0)
                    ! Inside domain
                    rho(l_center) = rho(l_center) + particleList(i)%q * particleList(i)%w_p * (0.75 - d**2)
                    if (l_right /= NumberXNodes) then
                        rho(l_right) = rho(l_right) + particleList(i)%q * particleList(i)%w_p * 0.5d0 * (0.5d0 + d)**2
                    else
                        SELECT CASE (world%boundaryConditions(l_right))
                        CASE(1)
                            rho(l_right) = rho(l_right) + particleList(i)%q * particleList(i)%w_p * 0.5d0 * (0.5d0 + d)**2
                        CASE(2)
                            rho(l_right) = rho(l_right) + particleList(i)%q * particleList(i)%w_p * (0.5d0 + d)**2
                        CASE(3)
                            rho(l_right) = rho(l_right) + particleList(i)%q * particleList(i)%w_p * 0.5d0 * (0.5d0 + d)**2
                            rho(1) = rho(1) + particleList(i)%q * particleList(i)%w_p * 0.5d0 * (0.5d0 + d)**2
                        END SELECT
                    end if
                    if (l_left /= 1) then
                        rho(l_left) = rho(l_left) + particleList(i)%q * particleList(i)%w_p * 0.5d0 * (0.5d0 - d)**2
                    else
                        SELECT CASE (world%boundaryConditions(l_left))
                        CASE(1)
                            rho(l_left) = rho(l_left) + particleList(i)%q * particleList(i)%w_p * 0.5d0 * (0.5d0 - d)**2
                        CASE(2)
                            rho(l_left) = rho(l_left) + particleList(i)%q * particleList(i)%w_p * (0.5d0 - d)**2
                        CASE(3)
                            rho(l_left) = rho(l_left) + particleList(i)%q * particleList(i)%w_p * 0.5d0 * (0.5d0 - d)**2
                            rho(NumberXNodes) = rho(NumberXNodes) + particleList(i)%q * particleList(i)%w_p * 0.5d0 * (0.5d0 - d)**2
                        END SELECT
                    end if
                CASE(1)
                    !Dirichlet
                    rho(l_center) = rho(l_center) + particleList(i)%q * particleList(i)%w_p * (1.0d0-ABS(d))
                    rho(l_center + INT(SIGN(1.0, d))) = rho(l_center + INT(SIGN(1.0, d))) + particleList(i)%q * particleList(i)%w_p * ABS(d)
                CASE(2)
                    !Neumann symmetric
                    rho(l_center) = rho(l_center) + 2.0d0 * particleList(i)%q * particleList(i)%w_p * (0.75 - d**2)
                    rho(l_center + INT(SIGN(1.0, d))) = rho(l_center + INT(SIGN(1.0, d))) + particleList(i)%q * particleList(i)%w_p * (0.25d0 + d**2)
                CASE(3)
                    ! Periodic
                    rho(l_center) = rho(l_center) + particleList(i)%q * particleList(i)%w_p * (0.75 - d**2)
                    ! towards domain
                    rho(l_center+INT(SIGN(1.0, d))) = rho(l_center+INT(SIGN(1.0, d))) + particleList(i)%q * particleList(i)%w_p * 0.5d0 * (0.5d0 + ABS(d))**2
                    ! across periodic boundary
                    rho(MODULO(l_center-2*INT(SIGN(1.0, d)),NumberXNodes)) = rho(MODULO(l_center-2*INT(SIGN(1.0, d)),NumberXNodes)) + particleList(i)%q * particleList(i)%w_p * 0.5d0 * (0.5d0 - ABS(d))**2
                END SELECT
            end do
        end do
        rho = rho / world%nodeVol
        ! if (world%boundaryConditions(1) == 3) then
        !     rho(1) = rho(1) + rho(NumberXNodes)
        !     rho(NumberXNodes) = rho(1)
        ! end if
        ! if (world%boundaryConditions(1) == 2) then
        !     rho(1) = 2.0d0 * rho(1)
        ! else if (world%boundaryConditions(NumberXNodes) == 2) then
        !     rho(NumberXNodes) = 2.0d0 * rho(NumberXNodes)
        ! end if
    end subroutine depositRho

    subroutine loadParticleDensity(densities, particleList, world)
        type(Particle), intent(in) :: particleList(:)
        type(Domain), intent(in) :: world
        real(real64), intent(in out) :: densities(:,:)
        integer(int32) :: i,j, l_center, l_left, l_right
        real(real64) :: d
        do i=1, numberChargedParticles
            do j = 1, particleList(i)%N_p
                l_center = NINT(particleList(i)%phaseSpace(1, j))
                l_left = l_center -1 
                l_right = l_center + 1
                d = particleList(i)%phaseSpace(1, j) - l_center
                SELECT CASE (world%boundaryConditions(l_center))
                CASE(0)
                    ! Inside domain
                    densities(l_center, i) = densities(l_center, i) + particleList(i)%w_p * (0.75 - d**2)
                    if (l_right /= NumberXNodes) then
                        densities(l_right, i) = densities(l_right, i) + particleList(i)%w_p * 0.5d0 * (0.5d0 + d)**2
                    else
                        SELECT CASE (world%boundaryConditions(l_right))
                        CASE(1)
                            densities(l_right, i) = densities(l_right, i) + particleList(i)%w_p * 0.5d0 * (0.5d0 + d)**2
                        CASE(2)
                            densities(l_right, i) = densities(l_right, i) +particleList(i)%w_p * (0.5d0 + d)**2
                        CASE(3)
                            densities(l_right, i) = densities(l_right, i) + particleList(i)%w_p * 0.5d0 * (0.5d0 + d)**2
                            densities(1, i) = densities(1, i) +  particleList(i)%w_p * 0.5d0 * (0.5d0 + d)**2
                        END SELECT
                    end if
                    if (l_left /= 1) then
                        densities(l_left, i) = densities(l_left, i) +  particleList(i)%w_p * 0.5d0 * (0.5d0 - d)**2
                    else
                        SELECT CASE (world%boundaryConditions(l_left))
                        CASE(1)
                            densities(l_left, i) = densities(l_left, i) + particleList(i)%w_p * 0.5d0 * (0.5d0 - d)**2
                        CASE(2)
                            densities(l_left, i) = densities(l_left, i) +  particleList(i)%w_p * (0.5d0 - d)**2
                        CASE(3)
                            densities(l_left, i) = densities(l_left, i) + particleList(i)%w_p * 0.5d0 * (0.5d0 - d)**2
                            densities(NumberXNodes, i) = densities(NumberXNodes, i) + particleList(i)%w_p * 0.5d0 * (0.5d0 - d)**2
                        END SELECT
                    end if
                CASE(1)
                    !Dirichlet
                    densities(l_center, i) = densities(l_center, i) + particleList(i)%w_p * (1.0d0-ABS(d))
                    densities(l_center + INT(SIGN(1.0, d)), i) = densities(l_center + INT(SIGN(1.0, d)), i) + particleList(i)%w_p * ABS(d)
                CASE(2)
                    !Neumann symmetric
                    densities(l_center, i) = densities(l_center, i) + 2.0d0 * particleList(i)%w_p * (0.75 - d**2)
                    densities(l_center + INT(SIGN(1.0, d)), i) = densities(l_center + INT(SIGN(1.0, d)), i) + particleList(i)%w_p * (0.25d0 + d**2)
                CASE(3)
                    ! Periodic
                    densities(l_center, i) = densities(l_center, i) + particleList(i)%w_p * (0.75 - d**2)
                    ! towards domain
                    densities(l_center+INT(SIGN(1.0, d)), i) = densities(l_center+INT(SIGN(1.0, d)), i) + particleList(i)%w_p * 0.5d0 * (0.5d0 + ABS(d))**2
                    ! across periodic boundary
                    densities(MODULO(l_center-2*INT(SIGN(1.0, d)),NumberXNodes), i) = densities(MODULO(l_center-2*INT(SIGN(1.0, d)),NumberXNodes), i) + particleList(i)%w_p * 0.5d0 * (0.5d0 - ABS(d))**2
                END SELECT
            end do
        end do


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
            ! if (world%boundaryConditions(1) == 3) then
            !     densities(1,i) = densities(1,i) + densities(NumberXNodes,i)
            !     densities(NumberXNodes,i) = densities(1, i)
            ! else if (world%boundaryConditions(1) == 2) then
            !     densities(1,i) = densities(1,i)*2.0d0
            ! end if
    
            ! if (world%boundaryConditions(NumberXNodes) == 2) densities(NumberXNodes, i) = densities(NumberXNodes, i)*2.0d0
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

end module mod_Scheme