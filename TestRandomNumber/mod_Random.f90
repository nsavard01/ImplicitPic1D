!*****************************************************************************************
!>
!  64-bit version of the Mersenne Twister pseudorandom number generator.
!
!### History
!
!  * Contributors: Rémi Piatek, Takuji Nishimura, Makoto Matsumoto, Jacob Williams
!    See LICENSE file for details.
!
!### References
!
!  * T. Nishimura, "Tables of 64-bit Mersenne Twisters" ACM Transactions on Modeling and
!    Computer Simulation 10. (2000) 348--357.
!  * M. Matsumoto and T. Nishimura,
!    "Mersenne Twister: a 623-dimensionally equidistributed uniform pseudorandom number generator"
!    ACM Transactions on Modeling and Computer Simulation 8. (Jan. 1998) 3--30.
!  * Original source: http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/VERSIONS/FORTRAN/mt19937-64.f95

module mod_Random

  use, intrinsic :: iso_fortran_env

  implicit none

  private


  integer(int32),parameter :: IA = 16807, IM = 2147483647, IQ = 127773, IR = 2836
  real(real64), parameter :: am = 1.0d0/IM
  
  
  type,public :: randType

    !! main class for random number generator

    private

    integer(int32) :: ix, iy, irand

  contains

    private

    
    procedure, public :: initialize

    procedure, public :: getRand

  end type randType
  
contains

  subroutine initialize(self,seed)
    !! Initializes `me%mt(nn)` with a seed

    class(randType),intent(inout) :: self
    integer(int32), intent(in) :: seed
    self%ix = -1
    self%iy = -1
    self%irand = seed
    self%iy=ior(ieor(888889999,abs(self%irand)),1) 
    self%ix=ieor(777755555,abs(self%irand))
    self%irand=abs(self%irand)+1  

  end subroutine initialize
 
  function getRand(self) result(res)
    class(randType), intent(inout) :: self
    integer(int32) :: k
    real(real64) :: res
    self%ix=ieor(self%ix,ishft(self%ix,13))
    self%ix=ieor(self%ix,ishft(self%ix,-17))
    self%ix=ieor(self%ix,ishft(self%ix,5))
    k=self%iy/IQ 
    self%iy=IA*(self%iy-k*IQ)-IR*k
    if (self%iy < 0) self%iy=self%iy+IM
    res=am*ior(iand(IM,ieor(self%ix,self%iy)),1)
  end function getRand
end module mod_Random