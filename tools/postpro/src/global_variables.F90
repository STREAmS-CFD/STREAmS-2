module global_variables

use parameters
!
implicit none
integer, parameter         :: nv = 68
!
integer                    :: nx,ny
integer                    :: ng
integer                    :: flow_init
integer                    :: nstatloc
!
real(rkind)               :: Reynolds,Prandtl,rfac,Mach
real(rkind)               :: u0,rho0,p0,t0,Twall,mu0,k0,cp0,cv0,l0,c0,gam
real(rkind), dimension(:), allocatable  :: x, y, dxg, dyg
real(rkind), dimension(:,:,:), allocatable  :: wstat
integer, dimension(:),     allocatable      :: ixstat 
!
end module
