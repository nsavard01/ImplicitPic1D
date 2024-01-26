module mod_NullCollision

    use iso_fortran_env, only: int32, real64, output_unit
    use constants
    use mod_BasicFunctions
    use mod_particle
    use mod_targetParticle
    implicit none

    private
    public :: nullCollision

    ! Particle contains particle properties and stored values in phase space df
    type :: nullCollision
        integer(int32) :: numberCollisions, numberReactants, lengthArrays
        real(real64), allocatable :: energyArray(:), sigmaVArray(:, :), energyThreshold(:)
        real(real64) :: sigmaVMax, reducedMass
        integer(int32), allocatable :: collisionType(:), reactantsIndx(:), numberProducts(:), productsIndx(:,:)

    end type nullCollision


    interface nullCollision 
        procedure :: nullCollision_constructor
    end interface nullCollision

contains

    type(nullCollision) function nullCollision_constructor(numberReactants, numberCollisions, lengthArrays, red_mass, energyArray, sigmaVArray, energyThreshold, collisionType, reactantsIndx, numberProducts, productsIndx) result(self)
        ! Construct particle object, sizeIncrease is fraction larger stored array compared to initial amount of particles
        ! In future, use hash function for possible k = 1 .. Nx, m amount of boundaries, p = prime number  m < p < N_x. h(k) = (k%p)%m
        integer(int32), intent(in) :: numberCollisions, numberReactants, lengthArrays
        real(real64), intent(in) :: energyArray(lengthArrays), sigmaVArray(lengthArrays, numberCollisions), energyThreshold(numberCollisions), red_mass
        integer(int32), intent(in) :: collisionType(numberCollisions), reactantsIndx(numberReactants), numberProducts(numberCollisions), productsIndx(3,numberCollisions)
        self%numberCollisions = numberCollisions
        self%numberReactants = numberReactants
        self%lengthArrays = lengthArrays
        self%reducedMass = red_mass
        allocate(self%energyArray(lengthArrays), self%sigmaVArray(lengthArrays, numberCollisions), self%energyThreshold(numberCollisions), self%collisionType(numberCollisions), &
            self%reactantsIndx(numberReactants), self%numberProducts(numberCollisions), self%productsIndx(3, numberCollisions))
        self%energyArray = energyArray
        self%sigmaVArray = sigmaVArray
        self%energyThreshold = energyThreshold
        self%collisionType = collisionType
        self%reactantsIndx = reactantsIndx
        self%numberProducts = numberProducts
        self%productsIndx = productsIndx
        self%sigmaVMax = MAXVAL(SUM(self%sigmaVArray, DIM=2))
    end function nullCollision_constructor


end module mod_NullCollision