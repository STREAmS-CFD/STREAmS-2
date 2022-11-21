module derivatives
 use parameters
 contains

 subroutine ddyw_staggered(ny,u,y,ywall,dudyw)
! Conputes derivative using cetral differences and upwind closures up to 6th order
  implicit none
  integer, intent(in) :: ny
  real(rkind), dimension(ny), intent(in) :: u,y
  real(rkind), intent(in) :: ywall
  real(rkind), intent(out) :: dudyw
  real(rkind) :: d1,d2
!
  d1 =  y(1) - ywall 
  d2 =  y(2) - ywall 
  dudyw = (u(1)*d2**2 -u(2)*d1**2)/(d1*d2*(d2-d1))
!
 end subroutine ddyw_staggered
!
 subroutine ddy(ny,u,y,iorder,dudy)
  implicit none
  integer, intent(in) :: ny, iorder
  real(rkind) , dimension(ny), intent(in) :: u,y
  real(rkind) , dimension(ny), intent(out):: dudy
  real(rkind) , dimension(3,3):: cc
  real(rkind) , dimension(3,7,3):: cl,cr
  integer :: lmax,j,l,ll
  real(rkind) :: du,dy,du_l,dy_l,du_r,dy_r
! 
  lmax = iorder/2
  if (lmax>3) then
   write(*,*) 'Error: order of accuray for derivatives cannot be larger than 6'
   stop
  endif
!
! 
 cc(1,1) = 1._rkind/2._rkind
!
 cc(1,2)   =  2._rkind/3._rkind
 cc(2,2)   = -1._rkind/12._rkind
!
 cc(1,3)   =  3._rkind/4._rkind
 cc(2,3)   = -3._rkind/20._rkind
 cc(3,3)   =  1._rkind/60._rkind
!
! Left closure coefficients cl(j,l,lmax)
! lmax = 1, j=1
  cl(1,1,1) = -3._rkind/2._rkind
  cl(1,2,1) =  2._rkind
  cl(1,3,1) = -1._rkind/2._rkind
!
! lmax = 2, j=1
  cl(1,1,2) = -25._rkind/12._rkind
  cl(1,2,2) =  4._rkind
  cl(1,3,2) = -3._rkind
  cl(1,4,2) =  4._rkind/3._rkind 
  cl(1,5,2) = -1._rkind/4._rkind
!
! lmax = 2, j=2
  cl(2,1,2) = -1._rkind/4._rkind
  cl(2,2,2) = -5._rkind/6._rkind
  cl(2,3,2) =  3._rkind/2._rkind 
  cl(2,4,2) = -1._rkind/2._rkind
  cl(2,5,2) =  1._rkind/12._rkind
!  
! lmax = 3, j=1
  cl(1,1,3) = -49._rkind/20._rkind
  cl(1,2,3) =  6._rkind
  cl(1,3,3) = -15._rkind/2._rkind
  cl(1,4,3) =  20._rkind/3._rkind 
  cl(1,5,3) = -15._rkind/4._rkind
  cl(1,6,3) =  6._rkind/5._rkind 
  cl(1,7,3) = -1._rkind/6._rkind
! lmax = 3, j=2
  cl(2,1,3) = -1._rkind/6._rkind
  cl(2,2,3) = -77._rkind/60._rkind
  cl(2,3,3) =  5._rkind/2._rkind
  cl(2,4,3) = -5._rkind/3._rkind
  cl(2,5,3) =  5._rkind/6._rkind 
  cl(2,6,3) = -1._rkind/4._rkind
  cl(2,7,3) =  1._rkind/30._rkind 
! lmax = 3, j=3
  cl(3,1,3) =  1._rkind/30._rkind 
  cl(3,2,3) = -2._rkind/5._rkind
  cl(3,3,3) = -7._rkind/12._rkind
  cl(3,4,3) =  4._rkind/3._rkind
  cl(3,5,3) = -1._rkind/2._rkind
  cl(3,6,3) =  2._rkind/15._rkind
  cl(3,7,3) = -1._rkind/60._rkind
!
  do j=1+lmax,ny-lmax
   du = 0._rkind
   dy = 0._rkind
   do l=1,lmax
    du = du + cc(l,lmax)*(u(j+l)-u(j-l))
    dy = dy + cc(l,lmax)*(y(j+l)-y(j-l))
   enddo
   dudy(j) = du/dy
  enddo
!
  do j=1,lmax
   du_l = 0._rkind
   dy_l = 0._rkind
   du_r = 0._rkind
   dy_r = 0._rkind
   do l=-lmax,lmax
    ll = 1+l+lmax 
    du_l = du_l + cl(j,ll,lmax)*u(ll)
    dy_l = dy_l + cl(j,ll,lmax)*y(ll)
    du_r = du_r - cl(j,ll,lmax)*u(ny-(l+lmax))
    dy_r = dy_r - cl(j,ll,lmax)*y(ny-(l+lmax))
   enddo
   dudy(j)      = du_l/dy_l
   dudy(ny-j+1) = du_r/dy_r
  enddo

 end subroutine ddy

 subroutine interpolate(np,fjpl,xjpl,x0,f0,df0,d2f0)
!Interpolate at x0,f0 using points xjpl,fjpl
!In output returns f0 and its first and second derivatives
  use uty,              only : invmat
  implicit none
  integer, intent(in) :: np
  real(rkind), dimension(np), intent(in)    :: fjpl,xjpl
  real(rkind), intent(in)                   :: x0
  real(rkind), dimension(np)                :: sol
  real(rkind), dimension(np,np) :: amat
  real(rkind) :: f0, df0, d2f0, dx, den
  integer :: j,l
!
  if (np<3) then
   write(*,*) "Number of poits for interpolation must be larger or equal to 3"
   stop
  endif
!
  do j=1,np
   dx = xjpl(j) - x0 
   amat(j,1) = 1._rkind
   den = 1._rkind
   do l=2,np
    den = den*(l-1)
    amat(j,l) = dx**(l-1)/(den)
   enddo
  enddo
!
  call invmat(amat,np)
  sol = matmul(amat,fjpl)
  f0   = sol(1)
  df0  = sol(2)
  d2f0 = sol(3)
!
 end subroutine interpolate
end module
