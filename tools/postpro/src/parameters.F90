module parameters
!
use, intrinsic :: iso_fortran_env
implicit none

#ifdef SINGLE_PRECISION
integer, parameter :: rkind = REAL32
#else
integer, parameter :: rkind = REAL64
#endif
!
end module parameters
