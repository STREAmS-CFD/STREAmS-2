module postpro_cha
 use parameters
 use global_variables 
 use reader
 use derivatives
 use comp_transform

 contains 

 subroutine stats1d
 implicit none
 integer, parameter :: nasymm = 16 
 integer, parameter :: np     = 7  ! Number of points for interpolation
 integer, parameter :: naux   = 3  ! Number of auxiliary variables for compressibility transformations 
 real(rkind), dimension(ny)   :: u, dudy
 real(rkind), dimension(nv,ny) :: wstat1d
 real(rkind), dimension(nv,ny/2) :: w1dh
 real(rkind), dimension(np) :: fvec,yvec
 real(rkind), dimension(ny/2) :: yvd,yt,yv,ut,uvd,uv
 real(rkind), dimension(naux,ny/2) :: vaux
 integer, dimension(nasymm) :: asymm 
 real(rkind) :: yy,rho,dudyw,rhowall,drhodyw,d2rhodyw
 real(rkind) :: uwall,d2udyw
 real(rkind) :: qw,Twall,dTdyw,d2Tdyw,Ttau
 real(rkind) :: ufav(ny/2),vfav(ny/2),wfav(ny/2) 
 real(rkind) :: uu,muwall,nuwall,retau,deltav,utau,tauw
 real(rkind) :: pp,tt,tt2,pp2,rho2,rhou2,rhov2,rhow2,rhouv
 integer :: i,j,m
 character(20) :: tname


 call read_grid_cha

 wstat1d = 0._rkind
 do i=1,nx
  do j=1,ny
   do m=1,nv
    wstat1d(m,j) = wstat1d(m,j) + wstat(i,j,m)
   enddo
  enddo
 enddo
 wstat1d = wstat1d/nx
!
 asymm = (/3,14,19,31,36,38,40,42,44,50,52,54,55,60,64,65/) 
 call avg_walls(ny,nv,nasymm, asymm, wstat1d, w1dh)
! 
 do j=1,ny/2
  ufav(j) = w1dh(13,j)/w1dh(1,j)
 enddo
!
 call ddy(ny/2,ufav,y(1:ny/2),6,dudy)
!
!Interpolate rho at the wall
!
 fvec = w1dh(1,1:np)
 yvec   = y(1:np)
 call interpolate(np,fvec,yvec,-1._rkind,rhowall,drhodyw,d2rhodyw)
!
!Compute dudy at the wall
 fvec = ufav(1:np)
 call interpolate(np,fvec,yvec,-1._rkind,uwall,dudyw,d2udyw)
! 
!Compute dTdy at the wall
 fvec = w1dh(6,1:np) 
 call interpolate(np,fvec,yvec,-1._rkind,Twall,dTdyw,d2Tdyw)
!Uncoment for second order accurate approximation of dudyw
!call ddyw_staggered(ny/2,ufav,y(1:ny/2),-1._rkind,dudyw)
!
!
 muwall = mu0
 nuwall = muwall/rhowall
 tauw   = muwall*dudyw
 utau   = sqrt(tauw/rhowall)
 deltav = nuwall/utau
 retau  = l0/deltav
!
 Twall = t0
 qw    =-k0*dTdyw
 Ttau  = qw/(rhowall*cp0*utau)
!
! Compute compressibility transformations
!
 vaux(1,:) = w1dh(1,:)/rhowall 
 vaux(2,:) = w1dh(21,:)/nuwall 
 vaux(3,:) = w1dh(20,:)/muwall 
!
 tname = 'TrettelLarsson'
 call transform(ny/2,1._rkind+y(1:ny/2),ufav(1:ny/2),naux,vaux,yt,ut,tname)
!
 tname = 'vanDriest'
 call transform(ny/2,1._rkind+y(1:ny/2),ufav(1:ny/2),naux,vaux,yvd,uvd,tname)
!
 tname = 'Volpiani'
 call transform(ny/2,1._rkind+y(1:ny/2),ufav(1:ny/2),naux,vaux,yv,uv,tname)
!
 open(10,file='POSTPRO/channinfo.dat',form='formatted')
 write(10,*)'Retau,          rho_wall,         utau/u0,         Cf'
 write(10,*)retau,rhowall,utau/u0,2._rkind*tauw/(rho0*u0**2)
 close(10)
!
 open(10,file='POSTPRO/channstat.prof',form='formatted')
 do j=1,ny/2
  yy   = 1.+y(j)
  rho  = w1dh(1,j) 
  tt   = w1dh(6,j) 
  pp   = w1dh(5,j) 
  uu   = ufav(j)
  rhou2  = w1dh(16,j) - w1dh(13,j)**2/w1dh(1,j) 
  rhov2  = w1dh(17,j) - w1dh(14,j)**2/w1dh(1,j)
  rhow2  = w1dh(18,j) - w1dh(15,j)**2/w1dh(1,j) 
  rhouv  = w1dh(19,j) - w1dh(13,j)*w1dh(14,j)/w1dh(1,j)
  rho2   = w1dh(7,j)  - rho**2
  pp2    = w1dh(11,j) - pp**2 
  tt2    = w1dh(12,j) - tt**2 
  write(10,100)yy,          &         !1  y/h
               yy/deltav,   &         !2  y^+
               yt(j)/deltav,&         !3  y_T^+
               uu/u0,       &         !4  u/u_b
               uu/utau,     &         !5  u^+
               uvd(j)/utau, &         !6  u_VD^+
               ut(j)/utau,  &         !7  u_TL^+
               rho/rhowall, &         !8  rho/rhowall
               rhou2/tauw,  &         !9  rho/rhowall*u2^+   
               rhov2/tauw,  &         !10 rho/rhowall*v2^+  
               rhow2/tauw,  &         !11 rho/rhowall*w2^+ 
               rhouv/tauw,  &         !12 rho/rhowall*uv^+
               tt/Twall,    &         !13 Temp/Twall
               rho2/rhowall,&         !14 rho/rhowall
               tt2/Twall**2,&         !15 T2/Twall
               (tt-Twall)/Ttau,     & !16 T^+
               tt2/Ttau**2, &         !17 T2^+
               yv(j)/deltav,&         !18 y_V^+
               uv(j)/utau             !19 u_v^+
              

 enddo
 close(10)
 100  format(200ES20.10)
 end subroutine stats1d

 subroutine avg_walls(ny,nv,nasymm,asymm,wstat1d,w1dh)
!Average channel flow statistics on top and bottom wall
  implicit none
  integer, intent(in) :: ny,nv,nasymm
  real(rkind), dimension(nv,ny), intent(in) :: wstat1d
  real(rkind), dimension(nv,ny), intent(out):: w1dh
  integer, dimension(nasymm), intent(in) :: asymm
  integer :: m,j
!  
  do j=1,ny/2
   do m=1,nv
    if (any(m==asymm(:))) then
     w1dh(m,j) = 0.5_rkind*(wstat1d(m,j)-wstat1d(m,ny-j+1))
    else
     w1dh(m,j) = 0.5_rkind*(wstat1d(m,j)+wstat1d(m,ny-j+1))
    endif
   enddo
  enddo
!
 end subroutine avg_walls
end module
