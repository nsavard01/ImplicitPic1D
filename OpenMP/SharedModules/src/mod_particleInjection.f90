module mod_particleInjection
    use iso_fortran_env, only: int32, real64
    use constants
    use mod_BasicFunctions
    use mod_domain
    use mod_particle
    use mod_Scheme
    use omp_lib
    implicit none
    real(real64), allocatable :: energyAddColl(:)

contains

    subroutine addMaxwellianLostParticles(particleList, T_e, T_i, irand, world)
        ! Add power to all particles in domain
        type(Particle), intent(in out) :: particleList(2)
        type(Domain), intent(in) :: world
        integer(int32), intent(in out) :: irand(numThread)
        real(real64), intent(in) :: T_e, T_i
        integer(int32) :: i, iThread
        real(real64) :: x_random
        !$OMP parallel private(iThread, i, x_random)
        iThread = omp_get_thread_num() + 1
        do i=1, particleList(1)%delIdx(iThread)
            call getMaxwellianSample(particleList(1)%phaseSpace(2:4, particleList(1)%N_p(iThread) + i, iThread), particleList(1)%mass, T_e, irand(iThread))
            call getMaxwellianFluxSample(particleList(2)%phaseSpace(2:4, particleList(2)%N_p(iThread) + i, iThread), particleList(2)%mass, T_i, irand(iThread))
            energyAddColl(iThread) = energyAddColl(iThread) + (SUM(particleList(1)%phaseSpace(2:4, particleList(1)%N_p(iThread) + i, iThread)**2) * particleList(1)%mass * particleList(1)%w_p &
            + SUM(particleList(2)%phaseSpace(2:4, particleList(2)%N_p(iThread) + i, iThread)**2) * particleList(2)%mass * particleList(2)%w_p) * 0.5d0
            x_random = world%grid(NumberXNodes) * ran2(irand(iThread))
            particleList(1)%phaseSpace(1, particleList(1)%N_p(iThread) + i, iThread) = getLFromX(x_random, world)
            particleList(2)%phaseSpace(1, particleList(2)%N_p(iThread) + i, iThread) = particleList(1)%phaseSpace(1, particleList(1)%N_p(iThread) + i, iThread)
        end do
        particleList(2)%N_p(iThread) = particleList(2)%N_p(iThread) + particleList(1)%delIdx(iThread)
        particleList(1)%N_p(iThread) = particleList(1)%N_p(iThread) + particleList(1)%delIdx(iThread)
        !$OMP end parallel
    end subroutine addMaxwellianLostParticles

    subroutine refluxParticles(particleList, T_e, T_i, irand, world)
        type(Particle), intent(in out) :: particleList(2)
        type(Domain), intent(in) :: world
        integer(int32), intent(in out) :: irand(numThread)
        real(real64), intent(in) :: T_e, T_i
        integer(int32) :: i,j, iThread
        do j = 1, 2
            !$OMP parallel private(iThread, i)
            iThread = omp_get_thread_num() + 1
            do i = 1, particleList(j)%refIdx(iThread)
                energyAddColl(iThread) = energyAddColl(iThread) - (SUM(particleList(j)%phaseSpace(2:4, particleList(j)%refRecordIdx(i, iThread), iThread)**2) * particleList(j)%mass * particleList(j)%w_p) * 0.5d0
                if (j == 1) then
                    call getMaxwellianFluxSample(particleList(j)%phaseSpace(2:4, particleList(j)%refRecordIdx(i, iThread), iThread), particleList(j)%mass, T_e, irand(iThread))
                else
                    call getMaxwellianFluxSample(particleList(j)%phaseSpace(2:4, particleList(j)%refRecordIdx(i, iThread), iThread), particleList(j)%mass, T_i, irand(iThread))
                end if
                if (world%boundaryConditions(1) == 2) then
                    particleList(j)%phaseSpace(2, particleList(j)%refRecordIdx(i, iThread), iThread) = ABS(particleList(j)%phaseSpace(2, particleList(j)%refRecordIdx(i, iThread), iThread))
                else
                    particleList(j)%phaseSpace(2, particleList(j)%refRecordIdx(i, iThread), iThread) = -ABS(particleList(j)%phaseSpace(2, particleList(j)%refRecordIdx(i, iThread), iThread))
                end if
                energyAddColl(iThread) = energyAddColl(iThread) + (SUM(particleList(j)%phaseSpace(2:4, particleList(j)%refRecordIdx(i, iThread), iThread)**2) * particleList(j)%mass * particleList(j)%w_p) * 0.5d0
            end do
            !$OMP end parallel
        end do
    end subroutine refluxParticles

end module mod_particleInjection