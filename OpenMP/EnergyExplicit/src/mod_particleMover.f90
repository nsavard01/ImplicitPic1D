module mod_particleMover
    use iso_fortran_env, only: int32, real64, output_unit
    use constants
    use mod_BasicFunctions
    use mod_particle
    use mod_domain
    use mod_potentialSolver
    use mod_particleInjection
    use omp_lib
    implicit none
    ! Procedures for moving particles and depositing J 
    integer(int32), parameter, private :: m_Anderson_Particle = 2
    integer(int32), parameter, private :: maxPartIter = 50

contains

   


    ! -------------------------------------------- Particle mover without boolean checks for depositing J ------------------------------------------------------------

    ! subroutine moveParticles(solver, particleList, world, del_t)
    ! !particle mover 1
    !     ! particle mover to avoid the boolean checks which mostly don't happen when depositing J
    !     class(potentialSolver), intent(in out) :: solver
    !     type(Domain), intent(in) :: world
    !     type(Particle), intent(in out) :: particleList(:)
    !     real(real64), intent(in) :: del_t
    !     !a and c correspond to quadratic equations | l_alongV is nearest integer boundary along velocity component, away is opposite
    !     real(real64) :: l_f, l_sub, v_sub, v_f, x_i, x_f
    !     integer(int32) :: j, i, l_cell, iThread, delIdx, l_boundary, numIter, refIdx, N_p, start_cell
    !     logical :: BoundaryBool

    !     !$OMP parallel private(iThread, i, j, l_f, l_sub, v_sub, v_f, l_cell, delIdx, l_boundary, numIter, &
    !             refIdx, BoundaryBool, N_p, start_cell)

    !     iThread = omp_get_thread_num() + 1 
    !     loopSpecies: do j = 1, numberChargedParticles
    !         delIdx = 0
    !         refIdx = 0
    !         particleList(j)%wallLoss(:, iThread) = 0
    !         particleList(j)%energyLoss(:, iThread) = 0.0d0
    !         N_p = particleList(j)%N_p(iThread)
    !         loopParticles: do i = 1, N_p


    !             ! Initial phase coords
    !             v_sub = particleList(j)%phaseSpace(2,i,iThread)
    !             l_sub = particleList(j)%phaseSpace(1,i,iThread)

    !             ! cell rounded down
    !             l_cell = INT(l_sub)
                    
    !             ! Final velocity coords
    !             v_f = v_sub + particleList(j)%q_over_m * solver%EField(l_cell) * del_t
    !             ! logical to phisical conversion
    !             x_i = world%grid(l_cell) + (l_sub - l_cell) * world%dx_dl(l_cell)
    !             x_f = x_i + del_t * v_f

    !             ! positive motion
    !             if(v_f >= 0) then
    !                 ! boundary check
    !                 if(x_f > world%grid(NumberXNodes)) then
    !                     ! periodic
    !                     if(world%boundaryConditions(NumberXNodes) == 3) then
    !                         x_f = MODULO(x_f, world%L_domain) + world%grid(1)
    !                         do start_cell = 2, NumberXNodes
    !                             if(world%grid(start_cell) > x_f) then
    !                                 l_f = (start_cell - 1) + (x_f - world%grid(start_cell - 1)) / world%dx_dl(start_cell-1)
    !                                 particleList(j)%phaseSpace(1, i-delIdx, iThread) = l_f
    !                                 particleList(j)%phaseSpace(2,i-delIdx, iThread) = v_f
    !                                 particleList(j)%phaseSpace(3:4,i-delIdx, iThread) = particleList(j)%phaseSpace(3:4,i,iThread)
    !                                 exit
    !                             endif
    !                         enddo
    !                     else
    !                         ! Dirichlet
    !                         if(world%boundaryConditions(NumberXNodes) == 1 .or. world%boundaryConditions(NumberXNodes) == 4) then
    !                             delIdx = delIdx + 1
    !                             particleList(j)%energyLoss(2, iThread) = particleList(j)%energyLoss(2, iThread) + v_f**2 + (SUM(particleList(j)%phaseSpace(3:4,i,iThread)**2)) !J/m^2 in 1D
    !                             particleList(j)%wallLoss(2, iThread) = particleList(j)%wallLoss(2, iThread) + 1 !C/m^2 in 1D
    !                             particleList(j)%momentumLoss(2, iThread) = particleList(j)%momentumLoss(2, iThread) + v_f
    !                         else
    !                         ! nuemann
    !                             x_f = -x_f + 2 * world%grid(NumberXNodes)
    !                             v_f = -v_f
    !                             ! dirichelt on other side
    !                             if(x_f <= world%grid(1)) then
    !                                 delIdx = delIdx + 1
    !                                 particleList(j)%energyLoss(1, iThread) = particleList(j)%energyLoss(1, iThread) + v_f**2 + (SUM(particleList(j)%phaseSpace(3:4,i,iThread)**2))
    !                                 particleList(j)%wallLoss(1, iThread) = particleList(j)%wallLoss(1, iThread) + 1 !C/m^2 in 1D
    !                                 particleList(j)%momentumLoss(1, iThread) = particleList(j)%momentumLoss(1, iThread) + v_f
    !                             else
    !                                 do start_cell = l_cell, 1, -1
    !                                     if(world%grid(start_cell) < x_f) then
    !                                         l_f = (start_cell) + (x_f - world%grid(start_cell)) / world%dx_dl(start_cell)
    !                                         particleList(j)%phaseSpace(1, i-delIdx, iThread) = l_f
    !                                         particleList(j)%phaseSpace(2,i-delIdx, iThread) = v_f
    !                                         particleList(j)%phaseSpace(3:4,i-delIdx, iThread) = particleList(j)%phaseSpace(3:4,i,iThread)
    !                                         exit
    !                                     end if
    !                                 end do
    !                             end if
    !                         end if
    !                     end if
    !                 else
    !                     do start_cell = l_cell + 1, NumberXNodes
    !                         if(world%grid(start_cell) > x_f) then
    !                             l_f = (start_cell - 1) + (x_f - world%grid(start_cell - 1)) / world%dx_dl(start_cell-1)
    !                             particleList(j)%phaseSpace(1, i-delIdx, iThread) = l_f
    !                             particleList(j)%phaseSpace(2,i-delIdx, iThread) = v_f
    !                             particleList(j)%phaseSpace(3:4,i-delIdx, iThread) = particleList(j)%phaseSpace(3:4,i,iThread)
    !                             exit
    !                         endif
    !                     enddo
    !                 end if
    !             else
    !                 if(x_f < world%grid(1) ) then
    !                 ! periodic
    !                     if(world%boundaryConditions(1) == 3) then
    !                         x_f = MODULO(x_f, world%L_domain) + world%grid(1)
    !                         do start_cell = NumberXHalfNodes, 1, -1
    !                             if(world%grid(start_cell) < x_f) then
    !                                 l_f = (start_cell) + (x_f - world%grid(start_cell)) / world%dx_dl(start_cell)
    !                                 particleList(j)%phaseSpace(1, i-delIdx, iThread) = l_f
    !                                 particleList(j)%phaseSpace(2,i-delIdx, iThread) = v_f
    !                                 particleList(j)%phaseSpace(3:4,i-delIdx, iThread) = particleList(j)%phaseSpace(3:4,i,iThread)
    !                                 exit
    !                             endif
    !                         enddo
    !                     else
    !                         ! Dirichlet
    !                         if(world%boundaryConditions(1) == 1 .or. world%boundaryConditions(1) == 4) then
    !                             delIdx = delIdx + 1
    !                             particleList(j)%energyLoss(1, iThread) = particleList(j)%energyLoss(1, iThread) + v_f**2 + (SUM(particleList(j)%phaseSpace(3:4,i,iThread)**2))
    !                             particleList(j)%wallLoss(1, iThread) = particleList(j)%wallLoss(1, iThread) + 1 !C/m^2 in 1D
    !                             particleList(j)%momentumLoss(1, iThread) = particleList(j)%momentumLoss(1, iThread) + v_f
    !                         else
    !                         ! nuemann
    !                             x_f = -x_f + 2 * world%grid(1)
    !                             v_f = -v_f
    !                             ! dirichelt on other side
    !                             if(x_f >= world%grid(NumberXNodes)) then
    !                                 delIdx = delIdx + 1
    !                                 particleList(j)%energyLoss(2, iThread) = particleList(j)%energyLoss(2, iThread) + v_f**2 + (SUM(particleList(j)%phaseSpace(3:4,i,iThread)**2)) !J/m^2 in 1D
    !                                 particleList(j)%wallLoss(2, iThread) = particleList(j)%wallLoss(2, iThread) + 1 !C/m^2 in 1D
    !                                 particleList(j)%momentumLoss(2, iThread) = particleList(j)%momentumLoss(2, iThread) + v_f
    !                             else
    !                                 do start_cell = l_cell + 1, NumberXHalfNodes
    !                                     if(world%grid(start_cell) > x_f) then
    !                                         l_f = (start_cell - 1) + (x_f - world%grid(start_cell - 1)) / world%dx_dl(start_cell-1)
    !                                         particleList(j)%phaseSpace(1, i-delIdx, iThread) = l_f
    !                                         particleList(j)%phaseSpace(2,i-delIdx, iThread) = v_f
    !                                         particleList(j)%phaseSpace(3:4,i-delIdx, iThread) = particleList(j)%phaseSpace(3:4,i,iThread)
    !                                         exit
    !                                     endif
    !                                 enddo
    !                             end if
    !                         end if
    !                     end if
    !                 else
    !                     do start_cell = l_cell, 1, -1
    !                         if(world%grid(start_cell) < x_f) then
    !                             l_f = (start_cell) + (x_f - world%grid(start_cell)) / world%dx_dl(start_cell)
    !                             particleList(j)%phaseSpace(1, i-delIdx, iThread) = l_f
    !                             particleList(j)%phaseSpace(2,i-delIdx, iThread) = v_f
    !                             particleList(j)%phaseSpace(3:4,i-delIdx, iThread) = particleList(j)%phaseSpace(3:4,i,iThread)
    !                             exit
    !                         endif
    !                     enddo
    !                 endif
    !             endif

    !         end do loopParticles
    !         particleList(j)%N_p(iThread) = N_p - delIdx
    !         particleList(j)%delIdx(iThread) = delIdx
    !         particleList(j)%refIdx(iThread) = refIdx
            
    !     end do loopSpecies
        
    !     !$OMP end parallel
    !     do j = 1, numberChargedParticles
    !         particleList(j)%numToCollide = particleList(j)%N_p
    !         particleList(j)%accumEnergyLoss = particleList(j)%accumEnergyLoss + SUM(particleList(j)%energyLoss, DIM = 2)
    !         particleList(j)%accumWallLoss = particleList(j)%accumWallLoss + SUM(particleList(j)%wallLoss, DIM=2)
    !     end do
    ! end subroutine moveParticles


    ! subroutine moveParticles(solver, particleList, world, del_t)
    ! ! particle mover 2
    !     ! particle mover to avoid the boolean checks which mostly don't happen when depositing J
    !     class(potentialSolver), intent(in out) :: solver
    !     type(Domain), intent(in) :: world
    !     type(Particle), intent(in out) :: particleList(:)
    !     real(real64), intent(in) :: del_t
    !     !a and c correspond to quadratic equations | l_alongV is nearest integer boundary along velocity component, away is opposite
    !     real(real64) :: l_f, l_sub, v_sub, v_f, x_i, x_f
    !     integer(int32) :: j, i, l_cell, iThread, delIdx, l_boundary, numIter, refIdx, N_p, start_cell, end_cell, loop_var, step
    !     logical :: BoundaryBool

    !     !$OMP parallel private(iThread, i, j, l_f, l_sub, v_sub, v_f, l_cell, delIdx, l_boundary, numIter, &
    !             refIdx, BoundaryBool, N_p, start_cell, end_cell, loop_var, step)
    !     BoundaryBool = .true.
    !     iThread = omp_get_thread_num() + 1 
    !     loopSpecies: do j = 1, numberChargedParticles
    !         delIdx = 0
    !         refIdx = 0
    !         particleList(j)%wallLoss(:, iThread) = 0
    !         particleList(j)%energyLoss(:, iThread) = 0.0d0
    !         N_p = particleList(j)%N_p(iThread)
    !         loopParticles: do i = 1, N_p


    !             ! Initial phase coords
    !             v_sub = particleList(j)%phaseSpace(2,i,iThread)
    !             l_sub = particleList(j)%phaseSpace(1,i,iThread)

    !             ! cell rounded down
    !             l_cell = INT(l_sub)
                    
    !             ! Final velocity coords
    !             v_f = v_sub + particleList(j)%q_over_m * solver%EField(l_cell) * del_t
    !             ! logical to phisical conversion
    !             x_i = world%grid(l_cell) + (l_sub - l_cell) * world%dx_dl(l_cell)
    !             x_f = x_i + del_t * v_f

    !             ! positive motion
    !             if(v_f >= 0) then
    !                 ! boundary check
    !                 if(x_f >= world%grid(NumberXNodes)) then
    !                     ! periodic
    !                     if(world%boundaryConditions(NumberXNodes) == 3) then
    !                         x_f = MODULO(x_f, world%L_domain) + world%grid(1)
    !                         start_cell = 2
    !                         end_cell = NumberXNodes
    !                         step = 1
    !                     else if (world%boundaryConditions(NumberXNodes) == 1 .or. world%boundaryConditions(NumberXNodes) == 4) then
    !                     ! Dirichlet
    !                         delIdx = delIdx + 1
    !                         particleList(j)%energyLoss(2, iThread) = particleList(j)%energyLoss(2, iThread) + v_f**2 + (SUM(particleList(j)%phaseSpace(3:4,i,iThread)**2)) !J/m^2 in 1D
    !                         particleList(j)%wallLoss(2, iThread) = particleList(j)%wallLoss(2, iThread) + 1 !C/m^2 in 1D
    !                         particleList(j)%momentumLoss(2, iThread) = particleList(j)%momentumLoss(2, iThread) + v_f
    !                         !BoundaryBool = .false.
    !                         cycle
    !                     else
    !                     ! nuemann
    !                         refIdx = refIdx +1
    !                         x_f = -x_f + 2 * world%grid(NumberXNodes)
    !                         v_f = -v_f
    !                         ! dirichelt on other side
    !                         if(x_f <= world%grid(1)) then
    !                             delIdx = delIdx + 1
    !                             particleList(j)%energyLoss(1, iThread) = particleList(j)%energyLoss(1, iThread) + v_f**2 + (SUM(particleList(j)%phaseSpace(3:4,i,iThread)**2))
    !                             particleList(j)%wallLoss(1, iThread) = particleList(j)%wallLoss(1, iThread) + 1 !C/m^2 in 1D
    !                             particleList(j)%momentumLoss(1, iThread) = particleList(j)%momentumLoss(1, iThread) + v_f
    !                             !BoundaryBool = .false.
    !                             cycle
    !                         else
    !                             start_cell = l_cell
    !                             end_cell = 1
    !                             step = -1
    !                         end if
    !                     end if
    !                 else
    !                     start_cell = l_cell + 1
    !                     end_cell = NumberXNodes
    !                     step = 1
    !                 end if
    !             else
    !                 if(x_f <= world%grid(1) ) then
    !                 ! periodic
    !                     if(world%boundaryConditions(1) == 3) then
    !                         x_f = MODULO(x_f, world%L_domain) + world%grid(1)
    !                         start_cell = NumberXHalfNodes
    !                         end_cell = 1
    !                         step = -1
    !                     else if (world%boundaryConditions(1) == 1 .or. world%boundaryConditions(1) == 4) then
    !                         ! Dirichlet
    !                         delIdx = delIdx + 1
    !                         particleList(j)%energyLoss(1, iThread) = particleList(j)%energyLoss(1, iThread) + v_f**2 + (SUM(particleList(j)%phaseSpace(3:4,i,iThread)**2))
    !                         particleList(j)%wallLoss(1, iThread) = particleList(j)%wallLoss(1, iThread) + 1 !C/m^2 in 1D
    !                         particleList(j)%momentumLoss(1, iThread) = particleList(j)%momentumLoss(1, iThread) + v_f
    !                         !BoundaryBool = .false.
    !                         cycle
    !                     else
    !                         ! nuemann
    !                         refIdx = refIdx +1
    !                         x_f = -x_f + 2 * world%grid(1)
    !                         v_f = -v_f
    !                         ! dirichelt on other side
    !                         if(x_f >= world%grid(NumberXNodes)) then
    !                             delIdx = delIdx + 1
    !                             particleList(j)%energyLoss(2, iThread) = particleList(j)%energyLoss(2, iThread) + v_f**2 + (SUM(particleList(j)%phaseSpace(3:4,i,iThread)**2)) !J/m^2 in 1D
    !                             particleList(j)%wallLoss(2, iThread) = particleList(j)%wallLoss(2, iThread) + 1 !C/m^2 in 1D
    !                             particleList(j)%momentumLoss(2, iThread) = particleList(j)%momentumLoss(2, iThread) + v_f
    !                             !BoundaryBool = .false.
    !                             cycle
    !                         else
    !                             start_cell = l_cell + 1
    !                             end_cell = NumberXNodes
    !                             step = 1
    !                         end if
    !                     end if
    !                 else
    !                     start_cell = l_cell
    !                     end_cell = 1
    !                     step = -1
    !                 endif
    !             endif

    !             ! do while(step * x_f > step * world%grid(start_cell))
    !             !     start_cell = start_cell + step
    !             ! end do
    !             ! end_cell = start_cell - step
    !             ! if (v_f >=0) then
    !             !     l_f = end_cell + ((x_f - world%grid(end_cell)) / world%dx_dl(end_cell))
    !             ! else
    !             !     l_f = end_cell - ((world%grid(end_cell) - x_f) / world%dx_dl(end_cell-1))
    !             ! end if

    !             if(v_f >= 0) then
    !                 l_f_loop1: do loop_var = start_cell, end_cell, step
    !                     if(world%grid(loop_var) > x_f) then
    !                         l_f = (loop_var - 1) + (x_f - world%grid(loop_var - 1)) / world%dx_dl(loop_var-1)
    !                         particleList(j)%phaseSpace(1, i-delIdx, iThread) = l_f
    !                         particleList(j)%phaseSpace(2,i-delIdx, iThread) = v_f
    !                         particleList(j)%phaseSpace(3:4,i-delIdx, iThread) = particleList(j)%phaseSpace(3:4,i,iThread)
    !                         exit
    !                     endif
    !                 enddo l_f_loop1
    !             else
    !                 l_f_loop2: do loop_var = start_cell, end_cell, step
    !                     if(world%grid(loop_var) < x_f) then
    !                         l_f = (loop_var) + (x_f - world%grid(loop_var)) / world%dx_dl(loop_var)
    !                         particleList(j)%phaseSpace(1, i-delIdx, iThread) = l_f
    !                         particleList(j)%phaseSpace(2,i-delIdx, iThread) = v_f
    !                         particleList(j)%phaseSpace(3:4,i-delIdx, iThread) = particleList(j)%phaseSpace(3:4,i,iThread)
    !                         exit
    !                     endif
    !                 enddo l_f_loop2
    !             end if

    !         end do loopParticles
    !         particleList(j)%N_p(iThread) = N_p - delIdx
    !         particleList(j)%delIdx(iThread) = delIdx
    !         particleList(j)%refIdx(iThread) = refIdx

            
    !     end do loopSpecies
        
    !     !$OMP end parallel
    !     do j = 1, numberChargedParticles
    !         particleList(j)%numToCollide = particleList(j)%N_p
    !         particleList(j)%accumEnergyLoss = particleList(j)%accumEnergyLoss + SUM(particleList(j)%energyLoss, DIM = 2)
    !         particleList(j)%accumWallLoss = particleList(j)%accumWallLoss + SUM(particleList(j)%wallLoss, DIM=2)
    !     end do
    ! end subroutine moveParticles

    ! subroutine moveParticles(solver, particleList, world, del_t)
    !     ! Particle mover version 3
    !     ! initialize classes and variables
    !     class(potentialSolver), intent(in out) :: solver
    !     type(Domain), intent(in) :: world
    !     type(Particle), intent(in out) :: particleList(:)
    !     real(real64), intent(in) :: del_t
    !     real(real64) :: l_f, l_sub, v_sub, v_f, distance, dx, pos
    !     integer(int32) :: j, i, l_cell, iThread, delIdx, l_boundary, numIter, refIdx, N_p, v_flag

    !     ! OMP parallel start
    !     !$OMP parallel private(iThread, i, j, l_f, l_sub, v_sub, v_f, l_cell, delIdx, l_boundary, numIter, pos, &
    !     !$OMP& refIdx, N_p, distance, v_flag, dx)
    !     ! get thread number for indexing
    !     iThread = omp_get_thread_num() + 1

    !     ! loop over particle types
    !     loopSpecies: do j = 1, numberChargedParticles

    !         ! count of deleted particles (dirichlet boundary conditions)
    !         delIdx = 0
    !         ! count of reflected particles (neuman boundary conditions)
    !         refIdx = 0

    !         ! reset loss counters for dirichlet walls
    !         particleList(j)%wallLoss(:, iThread) = 0
    !         particleList(j)%energyLoss(:, iThread) = 0
    !         particleList(j)%momentumLoss(:,iThread) = 0

    !         ! number of particles
    !         N_p = particleList(j)%N_p(iThread)

    !         ! loop over number of particles
    !         loopParticles: do i = 1, N_p

    !             ! Initial phase coords
    !             v_sub = particleList(j)%phaseSpace(2,i,iThread)
    !             l_sub = particleList(j)%phaseSpace(1,i,iThread)
   
    !             ! left side node
    !             l_cell = INT(l_sub)
                    
    !             ! Final velocity
    !             v_f = v_sub + particleList(j)%q_over_m * solver%EField(l_cell) * del_t

    !             ! distance to travel in real space
                
        
    !             ! set flags for left or right movment
    !             if(v_f >= 0) then
    !                 ! positive velocity flag
    !                 v_flag = 1
    !                 ! distance to right edge of cell
    !                 dx = world%dx_dl(l_cell) * (l_cell + 1 + - l_sub)

    !             else
    !                 ! negative velocity flag
    !                 v_flag = -1
    !                 ! distance to left egde of cell
    !                 dx = world%dx_dl(l_cell) * (l_sub - l_cell)
    !                 ! start at right side node
    !                 l_cell = l_cell + 1
                    
    !             end if
    !             distance = v_flag * del_t * v_f
    !             ! current possition
    !             pos = l_sub

    !             ! loop to find final position and check boundary conditions
    !             ! move the particle to the boundary of the next cell, if it wouldn't make it then end loop
    !             do while(distance > dx)
    !                 ! set cell to next cell by adding the flag (1 or -1)
    !                 l_cell = l_cell + v_flag

    !                 ! left boundary check
    !                 if(l_cell == 1 .and. v_flag == -1) then
    !                     ! periodic
    !                     if(world%boundaryConditions(1) == 3) then
    !                         ! move cell to other boundary minus one of left side
    !                         l_cell = NumberXHalfNodes
    !                     ! nueman
    !                     elseif(world%boundaryConditions(1) == 2) then
    !                         ! add to reflected particle count
    !                         refIdx = refIdx +1
    !                         ! flip velocity flag
    !                         v_flag = -v_flag
    !                     ! Dirichlet   
    !                     else
    !                         ! add to deleted paticle count and update loss
    !                         delIdx = delIdx + 1
    !                         particleList(j)%energyLoss(1, iThread) = particleList(j)%energyLoss(1, iThread) + v_f**2 + (SUM(particleList(j)%phaseSpace(3:4,i,iThread)**2))
    !                         particleList(j)%wallLoss(1, iThread) = particleList(j)%wallLoss(1, iThread) + 1 !C/m^2 in 1D
    !                         particleList(j)%momentumLoss(1, iThread) = particleList(j)%momentumLoss(1, iThread) + v_f
    !                         ! skip to next particle
    !                         cycle loopParticles
    !                     end if
    !                 endif
    !                 ! right boundary check
    !                 if(l_cell == NumberXNodes .and. v_flag == 1) then
    !                     ! periodic
    !                     if(world%boundaryConditions(NumberXNodes) == 3) then
    !                         ! move cell to other boundary
    !                         l_cell = 1
    !                     ! nueman
    !                     elseif(world%boundaryConditions(NumberXNodes) == 2) then
    !                         ! add to reflected particle count
    !                         refIdx = refIdx +1
    !                         ! flip velocity flag
    !                         v_flag = -v_flag
    !                     ! Dirichlet   
    !                     else
    !                         ! add to deleted paticle count and update loss
    !                         delIdx = delIdx + 1
    !                         particleList(j)%energyLoss(2, iThread) = particleList(j)%energyLoss(2, iThread) + v_f**2 + (SUM(particleList(j)%phaseSpace(3:4,i,iThread)**2)) !J/m^2 in 1D
    !                         particleList(j)%wallLoss(2, iThread) = particleList(j)%wallLoss(2, iThread) + 1 !C/m^2 in 1D
    !                         particleList(j)%momentumLoss(2, iThread) = particleList(j)%momentumLoss(2, iThread) + v_f
    !                         ! skip to next particle
    !                         cycle loopParticles
    !                     endif
    !                 endif

    !                 ! update distance to travel
    !                 distance = distance - dx

    !                 ! update dx to next cell width
    !                 dx = world%dx_dl(l_cell)

    !                 ! update curent position
    !                 pos = l_cell
    !             enddo

    !             ! move the remaining amount of distance depending on travel direction
    !             ! positive direction
    !             if(v_f < 0) then
    !                 l_cell = l_cell -1
    !             end if
                
    !             l_f = pos + v_flag * distance/world%dx_dl(l_cell)

    !             ! write new positions
    !             particleList(j)%phaseSpace(1, i-delIdx, iThread) = l_f
    !             particleList(j)%phaseSpace(2,i-delIdx, iThread) = v_f
    !             particleList(j)%phaseSpace(3:4,i-delIdx, iThread) = particleList(j)%phaseSpace(3:4,i,iThread)
                
    !         end do loopParticles
    !         ! update particle list length and reflect and deleteted particles
    !         particleList(j)%N_p(iThread) = N_p - delIdx
    !         particleList(j)%delIdx(iThread) = delIdx
    !         particleList(j)%refIdx(iThread) = refIdx

    !     end do loopSpecies
        
    !     !$OMP end parallel
        
    !     ! update loss values
    !     do j = 1, numberChargedParticles
    !         particleList(j)%numToCollide = particleList(j)%N_p
    !         particleList(j)%accumEnergyLoss = particleList(j)%accumEnergyLoss + SUM(particleList(j)%energyLoss, DIM = 2)
    !         particleList(j)%accumWallLoss = particleList(j)%accumWallLoss + SUM(particleList(j)%wallLoss, DIM=2)
    !     end do
    ! end subroutine moveParticles

subroutine moveParticles(solver, particleList, world, del_t)
    ! Particle mover version 4
    ! initialize classes and variables
    class(potentialSolver), intent(in out) :: solver
    type(Domain), intent(in) :: world
    type(Particle), intent(in out) :: particleList(:)
    real(real64), intent(in) :: del_t
    real(real64) :: l_f, l_sub, v_sub, v_f, x_i, x_f
    integer(int32) :: j, i, l_cell, iThread, delIdx, l_boundary, numIter, refIdx, N_p, v_flag
    logical :: refluxedBool, bool

    ! OMP parallel start
    !$OMP parallel private(iThread, i, j, l_f, l_sub, v_sub, v_f, l_cell, delIdx, l_boundary, numIter, &
    !$OMP& refIdx, N_p, v_flag, x_i, x_f, refluxedBool, bool)
    ! get thread number for indexing
    iThread = omp_get_thread_num() + 1

    ! loop over particle types
    loopSpecies: do j = 1, numberChargedParticles

        ! count of deleted particles (dirichlet boundary conditions)
        delIdx = 0
        ! count of reflected particles (neuman boundary conditions)
        refIdx = 0

        ! reset loss counters for dirichlet walls
        particleList(j)%wallLoss(:, iThread) = 0
        particleList(j)%energyLoss(:, iThread) = 0
        particleList(j)%momentumLoss(:,iThread) = 0

        ! number of particles
        N_p = particleList(j)%N_p(iThread)

        ! loop over number of particles
        loopParticles: do i = 1, N_p
            refluxedBool = .False.
            ! Initial phase coords
            v_sub = particleList(j)%phaseSpace(2,i,iThread)
            l_sub = particleList(j)%phaseSpace(1,i,iThread)

            ! left side node
            l_cell = INT(l_sub) 

            
            ! Final velocity
            v_f = v_sub + particleList(j)%q_over_m * solver%EField(l_cell) * del_t

            ! real postions
            x_i = world%grid(l_cell) + (l_sub - l_cell) * world%dx_dl(l_cell)
            x_f = x_i + del_t * v_f

            ! set flags for left or right movment
            if(v_f >= 0) then
                ! positive velocity flag
                v_flag = 1
                l_cell = l_cell + 1
            else
                ! negative velocity flag
                v_flag = -1
            end if
            
            bool = v_flag * x_f > v_flag * world%grid(l_cell)
            ! loop to find final position and check boundary conditions
            ! move the particle to the boundary of the next cell, if it wouldn't make it then end loop
            do while(bool)
                
                ! left boundary check
                if(l_cell == 1 .and. v_flag == -1) then
                    ! periodic
                    if(world%boundaryConditions(1) == 3) then
                        ! move cell to other boundary minus one of left side
                        l_cell = NumberXHalfNodes
                        x_f = x_f + world%L_domain
                    ! nueman
                    elseif(world%boundaryConditions(1) == 2) then
                        
                        ! flip velocity flag
                        v_flag = -v_flag
                        x_f = 2 * world%grid(1) - x_f

                        ! add to reflected particle count
                        refIdx = refIdx +1
                        particleList(j)%refRecordIdx(refIdx, iThread) = i - delIdx
                        refluxedBool = .true.
                    ! Dirichlet   
                    else
                        ! add to deleted paticle count and update loss
                        delIdx = delIdx + 1
                        particleList(j)%energyLoss(1, iThread) = particleList(j)%energyLoss(1, iThread) + v_f**2 + (SUM(particleList(j)%phaseSpace(3:4,i,iThread)**2))
                        particleList(j)%wallLoss(1, iThread) = particleList(j)%wallLoss(1, iThread) + 1 !C/m^2 in 1D
                        particleList(j)%momentumLoss(1, iThread) = particleList(j)%momentumLoss(1, iThread) + v_f
                        ! skip to next particle
                        cycle loopParticles
                    end if
                endif
                ! right boundary check
                if(l_cell == NumberXNodes .and. v_flag == 1) then
                    ! periodic
                    if(world%boundaryConditions(NumberXNodes) == 3) then
                        ! move cell to other boundary
                        l_cell = 1
                        x_f = x_f - world%L_domain
                    ! nueman
                    elseif(world%boundaryConditions(NumberXNodes) == 2) then
                        
                        ! flip velocity flag
                        v_flag = -v_flag
                        x_f = 2 * world%grid(NumberXNodes) - x_f

                        ! add to reflected particle count
                        refIdx = refIdx +1
                        particleList(j)%refRecordIdx(refIdx, iThread) = i - delIdx
                        refluxedBool = .true.
                    ! Dirichlet   
                    else
                        ! add to deleted paticle count and update loss
                        delIdx = delIdx + 1
                        particleList(j)%energyLoss(2, iThread) = particleList(j)%energyLoss(2, iThread) + v_f**2 + (SUM(particleList(j)%phaseSpace(3:4,i,iThread)**2)) !J/m^2 in 1D
                        particleList(j)%wallLoss(2, iThread) = particleList(j)%wallLoss(2, iThread) + 1 !C/m^2 in 1D
                        particleList(j)%momentumLoss(2, iThread) = particleList(j)%momentumLoss(2, iThread) + v_f
                        ! skip to next particle
                        cycle loopParticles
                    endif
                endif

                ! set cell to next cell by adding the flag (1 or -1)
                l_cell = l_cell + v_flag

                ! update distance to travel
                bool = v_flag * x_f > v_flag * world%grid(l_cell)
            enddo
            
            ! move the remaining amount of distance depending on travel direction
            ! positive direction
            if(v_f >= 0) then
                l_cell = l_cell -1
            end if
            l_f = l_cell + (x_f - world%grid(l_cell))/world%dx_dl(l_cell)

            ! write new positions
            particleList(j)%phaseSpace(1, i-delIdx, iThread) = l_f
            particleList(j)%phaseSpace(2,i-delIdx, iThread) = v_f
            particleList(j)%phaseSpace(3:4,i-delIdx, iThread) = particleList(j)%phaseSpace(3:4,i,iThread)
            
        end do loopParticles
        ! update particle list length and reflect and deleteted particles
        particleList(j)%N_p(iThread) = N_p - delIdx
        particleList(j)%delIdx(iThread) = delIdx
        particleList(j)%refIdx(iThread) = refIdx
    end do loopSpecies
    
    !$OMP end parallel
    
    ! update loss values
    do j = 1, numberChargedParticles
        particleList(j)%numToCollide = particleList(j)%N_p
        particleList(j)%accumEnergyLoss = particleList(j)%accumEnergyLoss + SUM(particleList(j)%energyLoss, DIM = 2)
        particleList(j)%accumWallLoss = particleList(j)%accumWallLoss + SUM(particleList(j)%wallLoss, DIM=2)
    end do
end subroutine moveParticles

subroutine moveParticlesEvenGrid(self, particleList, world, del_t)
    ! particle mover to avoid the boolean checks which mostly don't happen when depositing J
    class(potentialSolver), intent(in out) :: self
    type(Domain), intent(in) :: world
    type(Particle), intent(in out) :: particleList(:)
    real(real64), intent(in) :: del_t
    !a and c correspond to quadratic equations | l_alongV is nearest integer boundary along velocity component, away is opposite
    integer(int32) :: j, i, delIdx, refIdx, iThread, N_p
    real(real64) :: v_prime, q_over_m, partLoc, Del_x
    Del_x = world%dx_dl(1)
    !$OMP parallel private(iThread, i, j,delIdx, v_prime,partLoc, refIdx, N_p)
    iThread = omp_get_thread_num() + 1
    loopSpecies: do j = 1, numberChargedParticles
        delIdx = 0
        refIdx = 0
        particleList(j)%wallLoss(:, iThread) = 0
        particleList(j)%energyLoss(:, iThread) = 0.0d0
        N_p = particleList(j)%N_p(iThread)
        loopParticles: do i = 1, N_p
            ! First velocity change
            v_prime = particleList(j)%phaseSpace(2, i, iThread) + particleList(j)%q_over_m * self%getEField(particleList(j)%phaseSpace(1, i, iThread)) * del_t
            ! Get new position
            partLoc = particleList(j)%phaseSpace(1, i, iThread) + v_prime * del_t/Del_x
            if (partLoc <= 1) then
                SELECT CASE (world%boundaryConditions(1))
                CASE(1,4)
                    particleList(j)%energyLoss(1, iThread) = particleList(j)%energyLoss(1, iThread) + v_prime**2 + SUM(particleList(j)%phaseSpace(3:4, i, iThread)**2)
                    particleList(j)%wallLoss(1, iThread) = particleList(j)%wallLoss(1, iThread) + 1 !C/m^2 in 1D
                    particleList(j)%momentumLoss(1,iThread) = particleList(j)%momentumLoss(1,iThread) + v_prime
                    delIdx = delIdx + 1
                CASE(2)
                    particleList(j)%phaseSpace(1, i-delIdx, iThread) = 2.0d0 - partLoc
                    particleList(j)%phaseSpace(2, i-delIdx, iThread) = -v_prime
                    refIdx = refIdx + 1
                    particleList(j)%refRecordIdx(refIdx, iThread) = i - delIdx
                CASE(3)
                    particleList(j)%phaseSpace(1, i-delIdx, iThread) = MODULO(partLoc - 2.0d0, real(NumberXNodes, kind = real64)) + 1
                ! CASE default
                !     print *, 'no case, moveParticles'
                !     stop
                END SELECT
            else if ((partLoc >= NumberXNodes)) then
                SELECT CASE (world%boundaryConditions(NumberXNodes))
                CASE(1,4)
                    particleList(j)%energyLoss(2, iThread) = particleList(j)%energyLoss(2, iThread) + v_prime**2 + SUM(particleList(j)%phaseSpace(3:4, i, iThread)**2)
                    particleList(j)%wallLoss(2, iThread) = particleList(j)%wallLoss(2, iThread) + 1 !C/m^2 in 1D
                    particleList(j)%momentumLoss(2,iThread) = particleList(j)%momentumLoss(2,iThread) + v_prime
                    delIdx = delIdx + 1
                CASE(2)
                    particleList(j)%phaseSpace(1, i-delIdx, iThread) = 2.0d0 * NumberXNodes - partLoc
                    particleList(j)%phaseSpace(2, i-delIdx, iThread) = -v_prime
                    refIdx = refIdx + 1
                    particleList(j)%refRecordIdx(refIdx, iThread) = i - delIdx
                CASE(3)
                    particleList(j)%phaseSpace(1, i-delIdx, iThread) = MODULO(partLoc, real(NumberXNodes, kind = real64)) + 1
                ! CASE default
                !     print *, 'no case, moveParticles'
                !     stop
                END SELECT
            else
                particleList(j)%phaseSpace(1, i-delIdx, iThread) = partLoc
                particleList(j)%phaseSpace(2, i-delIdx, iThread) = v_prime
                particleList(j)%phaseSpace(3:4, i-delIdx, iThread) = particleList(j)%phaseSpace(3:4, i, iThread)
            end if
        end do loopParticles
        particleList(j)%N_p(iThread) = N_p - delIdx
        particleList(j)%delIdx(iThread) = delIdx
        particleList(j)%refIdx(iThread) = refIdx
    end do loopSpecies
    !$OMP end parallel
    do j = 1, numberChargedParticles
        particleList(j)%numToCollide = particleList(j)%N_p
        particleList(j)%accumEnergyLoss = particleList(j)%accumEnergyLoss + SUM(particleList(j)%energyLoss, DIM=2)
        particleList(j)%accumWallLoss = particleList(j)%accumWallLoss + SUM(particleList(j)%wallLoss, DIM=2)
    end do
end subroutine moveParticlesEvenGrid


end module mod_particleMover