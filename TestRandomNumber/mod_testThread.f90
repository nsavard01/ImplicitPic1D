module mod_testThread
    use omp_lib
    implicit none

    type,public :: testThread
  
      !! main class for random number generator
  
        integer :: testThreadVar

    end type testThread

    interface testThread
        module procedure :: init
    end interface testThread

contains

    type(testThread) function init(num) result(self)
        !! Initializes `me%mt(nn)` with a seed
        integer, intent(in) :: num
        integer :: iThread
        !$OMP threadprivate(self%testThreadVar)
        !$OMP parallel private(iThread)
        iThread = omp_get_thread_num() + 1
        testThreadVar = iThread + num
        !$OMP end parallel
        
    
    end function init
    
end module mod_testThread