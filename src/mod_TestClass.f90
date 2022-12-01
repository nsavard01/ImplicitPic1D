module mod_TestClass

    use iso_fortran_env, only: int32, real64, output_unit
    use constants
    use mod_particle
    implicit none

    private
    public :: testClass

    ! Particle contains particle properties and stored values in phase space
    type :: testClass
        type(Particle), allocatable :: particleList(:) !name of the particle

    contains
        procedure, public, pass(self) :: testFunction
    end type testClass


    interface testClass
        module procedure :: testClass_constructor
    end interface testClass

contains

    type(testClass) function testClass_constructor(pList) result(self)
        ! Construct particle object, sizeIncrease is fraction larger stored array compared to initial amount of particles
        type(Particle), intent(in) :: pList(:)
        allocate(self%particleList(size(pList)))
        self%particleList = pList

    end function testClass_constructor

    subroutine testFunction(self, index, T_e, irand)
        class(testClass), intent(in out) :: self
        integer(int32), intent(in out) :: irand
        integer(int32), intent(in) :: index
        real(real64), intent(in) :: T_e
        integer(int32) :: i
        do i = 1, 1000
            call self % particleList(index) % generate3DMaxwellian(T_e, irand)
        end do


    end subroutine testFunction


end module mod_TestClass