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

    function getChargeContinuityError(rho_i, rho_f, J_total, world, del_t) result(chargeError)
        real(real64), intent(in) :: del_t, rho_i(NumberXNodes), rho_f(NumberXNodes), J_total(NumberXNodes-1, numThread)
        type(Domain), intent(in) :: world
        integer(int32) :: i, k
        real(real64) :: chargeError, J(NumberXNodes-1), del_Rho
        J = SUM(J_total, DIM=2)
        chargeError = 0.0d0
        k = 0
        do i = 1, NumberXNodes
            del_Rho = rho_f(i) - rho_i(i)
            if (del_Rho /= 0) then
                SELECT CASE (world%boundaryConditions(i))
                CASE(0)
                    chargeError = chargeError + (1.0d0 + del_t * (J(i) - J(i-1))/del_Rho)**2
                    k = k + 1
                CASE(1)
                    continue
                CASE(2)
                    if (i == 1) then
                        chargeError = chargeError + (1.0d0 + del_t * J(1)/del_Rho)**2
                    else
                        chargeError = chargeError + (1.0d0 - del_t * J(NumberXNodes-1)/del_Rho)**2
                    end if
                    k = k + 1
                END SELECT
            end if
        end do
        chargeError = SQRT(chargeError/k)
    end function getChargeContinuityError


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
            densities = (SUM(particleList(i)%densities, DIM=2)/world%nodeVol) * particleList(i)%w_p
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