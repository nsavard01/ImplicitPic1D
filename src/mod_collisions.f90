module mod_collisions
    use iso_fortran_env, only: int32, real64
    use constants
    use mod_BasicFunctions
    use mod_particle
    use mod_domain
    implicit none

    ! Type of collisions, will likely need arrays which you can loop through which has source, target particle, choose from particle list
    ! Will likely have construct method which takes in collision data and turns into array

end module mod_collisions