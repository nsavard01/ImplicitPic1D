module mod_Scheme
    use iso_fortran_env, only: int32, real64
    use constants
    use mod_BasicFunctions
    use mod_particle
    use mod_domain
    use mod_potentialSolver
    use mod_collisions
    implicit none
    ! Scheme module for CIC
contains

    subroutine solveSingleTimeStepDiagnostic(solver, particleList, world, del_t, maxIter, eps_r)
        ! Single time step solver with Divergence of ampere, followed by adding of power, followed by collisions
        type(Particle), intent(in out) :: particleList(:)
        type(potentialSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t, eps_r
        integer(int32), intent(in) :: maxIter
        integer(int32) :: j, k
        real(real64) :: KE_i, KE_f, PE_i, PE_f, rho_f(NumberXNodes)

        ! Get charge/energy conservation error
        solver%particleEnergyLoss = 0.0d0
        PE_i = solver%getTotalPE(world, .false.)
        KE_i = 0.0d0
        do j=1, numberChargedParticles
            KE_i = KE_i + particleList(j)%getTotalKE()
        end do
        call solver%depositRho(particleList, world) 
        call solver%solveDivAmpereAnderson(particleList, world, del_t, maxIter, eps_r)
        KE_f = solver%particleEnergyLoss
        do j=1, numberChargedParticles
            KE_f = KE_f + particleList(j)%getTotalKE()
        end do
        PE_f = solver%getTotalPE(world, .false.)
        solver%energyError = ABS((KE_i + PE_i - KE_f - PE_f)/(KE_i + PE_i))
        call depositRhoDiag(rho_f, particleList, world)
        solver%chargeError = 0.0d0
        j = 0
        if (world%boundaryConditions(1) == 3) then
            j = j + 1
            solver%chargeError = solver%chargeError + (1 + (solver%J(1) - solver%J(NumberXNodes-1)) *del_t/ world%dx_dl(1)/(rho_f(1) - solver%rho(1)))**2
        end if
        do k = 1, NumberXNodes -2
            if ((rho_f(k+1) - solver%rho(k+1)) /= 0) then
                solver%chargeError = solver%chargeError + (1 + (solver%J(k + 1) - solver%J(k)) *del_t/ world%dx_dl(k+1)/(rho_f(k+1) - solver%rho(k+1)))**2
                j = j + 1
            end if
        end do
        solver%chargeError = SQRT(solver%chargeError/j)

    end subroutine solveSingleTimeStepDiagnostic

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
                numInCell = INT(part%N_p * world%dx_dl(i)/L_domain/2.0d0)
                call getRandom(part%phaseSpace(1,idxLower:idxLower + numInCell-1), irand)
                numPerCell(i) = numInCell
                if (i == 1) then
                    part%phaseSpace(1, idxLower:idxLower + numInCell - 1) = part%phaseSpace(1, idxLower:idxLower + numInCell - 1) * 0.5d0 + 1.0d0
                else
                    part%phaseSpace(1, idxLower:idxLower + numInCell - 1) = NumberXNodes - part%phaseSpace(1, idxLower:idxLower + numInCell - 1) * 0.5d0
                end if
            else
                numInCell = INT(part%N_p * world%dx_dl(i)/L_domain)
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


    subroutine depositRhoDiag(rho, particleList, world) 
        real(real64), intent(in out) :: rho(:)
        type(Particle), intent(in) :: particleList(:)
        type(Domain), intent(in) :: world
        integer(int32) :: i, j, l_center
        real(real64) :: d
        rho = 0.0d0
        do i=1, numberChargedParticles
            do j = 1, particleList(i)%N_p
                l_center = NINT(particleList(i)%phaseSpace(1, j))
                d = particleList(i)%phaseSpace(1, j) - l_center
                if (world%boundaryConditions(l_center) == 0) then
                    ! Inside domain
                    rho(l_center) = rho(l_center) + particleList(i)%q * particleList(i)%w_p * (0.75 - d**2)
                    rho(l_center + 1) = rho(l_center + 1) + particleList(i)%q * particleList(i)%w_p * 0.5d0 * (0.5d0 + d)**2
                    rho(l_center - 1) = rho(l_center - 1) + particleList(i)%q * particleList(i)%w_p * 0.5d0 * (0.5d0 - d)**2
                else if (world%boundaryConditions(l_center) == 1) then
                    !Dirichlet
                    rho(l_center) = rho(l_center) + particleList(i)%q * particleList(i)%w_p * (1.0d0-ABS(d))
                    rho(l_center + INT(SIGN(1.0, d))) = rho(l_center + INT(SIGN(1.0, d))) + particleList(i)%q * particleList(i)%w_p * ABS(d)
                
                else if (world%boundaryConditions(l_center) == 3) then
                    ! Periodic
                    rho(l_center) = rho(l_center) + particleList(i)%q * particleList(i)%w_p * (0.75 - d**2)
                    ! towards domain
                    rho(l_center+INT(SIGN(1.0, d))) = rho(l_center+INT(SIGN(1.0, d))) + particleList(i)%q * particleList(i)%w_p * 0.5d0 * (0.5d0 + ABS(d))**2
                    ! across periodic boundary
                    rho(MODULO(l_center-2*INT(SIGN(1.0, d)),NumberXNodes)) = rho(MODULO(l_center-2*INT(SIGN(1.0, d)),NumberXNodes)) + particleList(i)%q * particleList(i)%w_p * 0.5d0 * (0.5d0 - ABS(d))**2
                end if
            end do
        end do
        if (world%boundaryConditions(1) == 3) then
            rho(1) = rho(1) + rho(NumberXNodes)
            rho(NumberXNodes) = rho(1)
        end if
        rho = rho / world%dx_dl
    end subroutine depositRhoDiag

    subroutine loadParticleDensity(densities, particleList, world)
        type(Particle), intent(in) :: particleList(:)
        type(Domain), intent(in) :: world
        real(real64), intent(in out) :: densities(:,:)
        integer(int32) :: i,j, l_center
        real(real64) :: d
        do i=1, numberChargedParticles
            do j = 1, particleList(i)%N_p
                l_center = NINT(particleList(i)%phaseSpace(1, j))
                d = particleList(i)%phaseSpace(1, j) - l_center
                if (world%boundaryConditions(l_center) == 0) then
                    ! Inside domain
                    densities(l_center, i) = densities(l_center, i) + particleList(i)%w_p * (0.75 - d**2)
                    densities(l_center + 1, i) = densities(l_center + 1, i) + particleList(i)%w_p * 0.5d0 * (0.5d0 + d)**2
                    densities(l_center - 1, i) = densities(l_center - 1, i) + particleList(i)%w_p * 0.5d0 * (0.5d0 - d)**2
                else if (world%boundaryConditions(l_center) == 1) then
                    !Dirichlet
                    densities(l_center, i) = densities(l_center, i) + particleList(i)%w_p * (1.0d0-ABS(d))
                    densities(l_center + INT(SIGN(1.0, d)), i) = densities(l_center + INT(SIGN(1.0, d)), i) + particleList(i)%w_p * ABS(d)
                
                else if (world%boundaryConditions(l_center) == 3) then
                    ! Periodic
                    densities(l_center, i) = densities(l_center, i) + particleList(i)%w_p * (0.75 - d**2)
                    ! towards domain
                    densities(l_center+INT(SIGN(1.0, d)), i) = densities(l_center+INT(SIGN(1.0, d)), i) + particleList(i)%w_p * 0.5d0 * (0.5d0 + ABS(d))**2
                    ! across periodic boundary
                    densities(MODULO(l_center-2*INT(SIGN(1.0, d)),NumberXNodes), i) = densities(MODULO(l_center-2*INT(SIGN(1.0, d)),NumberXNodes), i) + particleList(i)%w_p * 0.5d0 * (0.5d0 - ABS(d))**2
                end if
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
            if (world%boundaryConditions(1) == 3) then
                densities(1,i) = densities(1,i) + densities(NumberXNodes, i)
                densities(NumberXNodes, i) = densities(1, i)
            end if
            densities(:,i) = densities(:,i)/world%dx_dl
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