module mod_simulation
    ! module which actually uses others to run a simulation (single or multiple time steps)
    ! will contain diagnostics as well
    use iso_fortran_env, only: int32, real64
    use constants
    use mod_BasicFunctions
    use mod_particle
    use mod_domain
    use mod_potentialSolver
    use mod_collisions

    integer(int32) :: maxIter = 50, numDiagnosticSteps = 50, numTimeSteps = 500
    real(real64) :: eps_r = 1e-8, del_t, fractionFreq = 0.5d0

contains
    
    ! -------------------------- Simulation ------------------------------------------

    subroutine solveSingleTimeStep(solver, particleList, world, del_t, maxIter, eps_r, irand, boolDiagnostic)
        ! Single time step solver with Divergence of ampere, followed by adding of power, followed by collisions
        type(Particle), intent(in out) :: particleList(:)
        type(potSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t, eps_r
        logical, intent(in) :: boolDiagnostic
        integer(int32), intent(in) :: maxIter
        integer(int32), intent(in out) :: irand
        !real(real64) :: P_before(3), P_after(3), E_before, E_after

        call solver%solveDivAmperePicard(particleList, world, del_t, maxIter, eps_r, boolDiagnostic)

        call addUniformPowerMaxwellianNicolas(particleList(1), Power, 0.05d0, irand, del_t)
        !call addUniformPowerMaxwellian(particleList(1), Power, nu_h, irand, del_t)
        ! P_before = m_e * SUM(particleList(1)%phaseSpace(2:4,:), DIM = 2) + particleList(2)%mass *SUM(particleList(2)%phaseSpace(2:4,:), DIM = 2)
        ! E_before = m_e * SUM(particleList(1)%phaseSpace(2:4,:)**2) * 0.5d0 / e + particleList(2)%mass * SUM(particleList(2)%phaseSpace(2:4, :)**2) * 0.5d0/e
        call ionizationCollisionIsotropic(particleList(1), particleList(2), 1.0d20, 1.0d-20, del_t, 15.8d0, irand)
        ! P_after = m_e * SUM(particleList(1)%phaseSpace(2:4,:), DIM = 2) + particleList(2)%mass *SUM(particleList(2)%phaseSpace(2:4,:), DIM = 2)
        ! E_after = m_e * SUM(particleList(1)%phaseSpace(2:4,:)**2) * 0.5d0 / e + particleList(2)%mass * SUM(particleList(2)%phaseSpace(2:4, :)**2) * 0.5d0/e + inelasticEnergyLoss
        ! print *, "P_before is:", P_before
        ! print *, "P_after is:", P_after
        ! print *, "E_before is:", E_before
        ! print *, "E_after is:", E_after
        ! print *, "Inelastic Energy loss is:", inelasticEnergyLoss
        ! stop "Just to check collision"


    end subroutine solveSingleTimeStep


    subroutine solveSimulation(solver, particleList, world, del_t, maxIter, eps_r, irand, numTimeSteps)
        ! Perform certain amount of timesteps, with diagnostics taken at first and last time step
        type(Particle), intent(in out) :: particleList(:)
        type(potSolver), intent(in out) :: solver
        type(Domain), intent(in) :: world
        real(real64), intent(in) :: del_t, eps_r
        integer(int32), intent(in) :: maxIter, numTimeSteps
        integer(int32), intent(in out) :: irand
        integer(int32) :: numSkipSteps, i, j,CurrentDiagStep = 1
        real(real64) :: currentTime
        open(15,file='../Data/InitialConditions.dat',access='APPEND')
        write(15,'("Number Grid Nodes, Final Expected Time(s), Delta t(s)")')
        write(15,"((I3.3, 1x), 2(es16.8,1x))") n_x, numTimeSteps*del_t, del_t
        close(15)
        close(22)
        numSkipSteps = numTimeSteps/(numDiagnosticSteps)
        101 format(20(1x,es16.8))
        open(22,file='../Data/ScalarDiagnosticData.dat',access='APPEND')
        write(22,'("Time (s), Collision Loss (W/m^2), ParticleCurrentLoss (A/m^2), ParticlePowerLoss(W/m^2), chargeError (a.u), energyError(a.u)")')

        !Save initial particle/field data, along with domain
        call solver%writeParticleDensity(particleList, world, 0) 
        call solver%writePhi(0)
        call particleList(1)%writeLocalTemperature(0)
        call world%writeDomain()
        do j=1, size(particleList)
            call particleList(j)%writePhaseSpace(0)
        end do

        do i = 1, numTimeSteps-1
            if (MOD((i-1), numSkipSteps) /= 0) then
                call solveSingleTimeStep(solver, particleList, world, del_t, maxIter, eps_r, irand, .false.)
            else  
                print *, "Simulation is", real(i)/numTimeSteps * 100.0, "percent done"
                currentTime = (i) * del_t
                call solveSingleTimeStep(solver, particleList, world, del_t, maxIter, eps_r, irand, .true.)
                call solver%writeParticleDensity(particleList, world, CurrentDiagStep) 
                call solver%writePhi(CurrentDiagStep)
                call particleList(1)%writeLocalTemperature(CurrentDiagStep)
                do j=1, size(particleList)
                    call particleList(j)%writePhaseSpace(CurrentDiagStep)
                end do
                write(22,101) currentTime, inelasticEnergyLoss*e/del_t, solver%particleCurrentLoss, solver%particlePowerLoss, solver%chargeError, solver%energyError
                CurrentDiagStep = CurrentDiagStep + 1
            end if
        end do
        currentTime = (numTimeSteps) * del_t
        call solveSingleTimeStep(solver, particleList, world, del_t, maxIter, eps_r, irand, .true.)
        call solver%writeParticleDensity(particleList, world, CurrentDiagStep) 
        call solver%writePhi(CurrentDiagStep)
        call particleList(1)%writeLocalTemperature(CurrentDiagStep)
        do j=1, size(particleList)
            call particleList(j)%writePhaseSpace(CurrentDiagStep)
        end do
        write(22,101) currentTime, inelasticEnergyLoss*e/del_t, solver%particleCurrentLoss, solver%particlePowerLoss, solver%chargeError, solver%energyError
        close(22)






    end subroutine solveSimulation







end module mod_simulation