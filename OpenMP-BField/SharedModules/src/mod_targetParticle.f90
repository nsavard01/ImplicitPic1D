module mod_targetParticle

    use iso_fortran_env, only: int32, real64, output_unit
    use constants
    use mod_BasicFunctions
    implicit none

    private
    public :: targetParticle, readNeutralParticleInputs
    integer(int32), public, protected :: numberNeutralParticles = 0

    ! target particle contains basic properties of basic neutral particle for target null collision
    type :: targetParticle
        character(:), allocatable :: name !name of the particle
        real(real64) :: mass, density, temperature, v_therm
    contains
        procedure, public, pass(self) :: generate3DMaxwellianVelocity
    end type targetParticle


    interface targetParticle
        module procedure :: targetParticle_constructor
    end interface targetParticle

contains

    type(targetParticle) function targetParticle_constructor(particleName, mass, density, temperature) result(self)
        ! Construct particle object, sizeIncrease is fraction larger stored array compared to initial amount of particles
        ! In future, use hash function for possible k = 1 .. Nx, m amount of boundaries, p = prime number  m < p < N_x. h(k) = (k%p)%m
        real(real64), intent(in) :: mass, density, temperature
        character(*), intent(in) :: particleName
        self % name = particleName
        self % mass = mass
        self % density = density ! #N/m^3
        self % temperature = temperature !temperature
        self % v_therm = SQRT(self%temperature*k_B/ self%mass)
    end function targetParticle_constructor


    ! ------------------------ Generating and initializing velocity with Maxwellian --------------------------


    function generate3DMaxwellianVelocity(self, irand) result(res)
        ! random velocity generator for the particle for temperature T (eV)
        ! Use box-muller method for random guassian variable, same as gwenael but doesn't have factor 2? Maybe factored into v_th
        class(targetParticle), intent(in) :: self
        integer(int32), intent(in out) :: irand
        real(real64) :: U1, U2, U3, U4, res(3)
        U1 = ran2(irand)
        U2 = ran2(irand)
        U3 = ran2(irand)
        U4 = ran2(irand)
        res(1) = self%v_therm * SQRT(-2 * LOG(U1)) * COS(2 * pi * U2)
        res(2) = self%v_therm * SQRT(-2 * LOG(U1)) * SIN(2 * pi * U2)
        res(3) = self%v_therm * SQRT(-2 * LOG(U3)) * SIN(2 * pi * U4)
    end function generate3DMaxwellianVelocity

    ! ------------------------ read in neutrals ------------------------------------------

    subroutine readNeutralParticleInputs(filename, targetParticleList)
        type(targetParticle), allocatable, intent(out) :: targetParticleList(:)
        character(len=*), intent(in) :: filename
        integer(int32) :: j, numNeutral = 0
        character(len=15) :: name
        character(len=8) :: neutralName
        real(real64) :: mass_neutral, Temp_neutral, density_neutral

        print *, "Reading particle inputs:"
        open(10,file='../InputData/'//filename, action = 'read')
        do j=1, 10000
            read(10,*) name


            if( name(1:4).eq.'NEUT' .or. name(1:4).eq.'Neut' .or. name(1:4).eq.'neut' ) then
                do while(name(1:4).ne.'----')
                    read(10,*) name
                end do
                read(10,'(A6)', ADVANCE = 'NO') name
                do while (name(1:4).ne.'----')
                    numNeutral = numNeutral + 1
                    read(10,*) mass_neutral, Temp_neutral, density_neutral
                    mass_neutral = mass_neutral * m_amu
                    neutralName = trim(name)
                    read(10,'(A6)', ADVANCE = 'NO') name
                end do
            endif       

            if (name(1:7) == 'ENDFILE') then
                close(10)
                exit
            end if

        end do


        numberNeutralParticles = numNeutral
        if (numberNeutralParticles > 0) then
            allocate(targetParticleList(numberNeutralParticles))
            do j = 1, numberNeutralParticles
                targetParticleList(j) = targetParticle(neutralName, mass_neutral, density_neutral, Temp_neutral)
                print *, 'Initializing target particle:', targetParticleList(j)%name
                print *, 'Particle mass is:', targetParticleList(j)%mass
                print *, 'Particle temperature(K) is:', targetParticleList(j)%temperature
                print *, 'Particle density is:', targetParticleList(j)%density
            end do
        end if

        
        print *, "---------------"
        print *, ""



    end subroutine readNeutralParticleInputs


end module mod_targetParticle