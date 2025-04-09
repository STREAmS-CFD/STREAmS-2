module postpro_cha
 use parameters
 use global_variables 
 use reader
 use derivatives
 use comp_transform
!
 contains 
!
 subroutine stats1d
 implicit none
 integer, parameter :: nasymm = 15 
 integer, parameter :: naux   =  6 ! Number of auxiliary variables for compressibility transformations 
 real(rkind), dimension(ny)   :: u
 real(rkind), dimension(nv,ny) :: wstat1d
 real(rkind), dimension(nv,ny/2) :: w1dh
 real(rkind), dimension(ny/2) :: yvd,yt,yv,yh,ut,uvd,uv,ug,uh
 real(rkind), dimension(naux,ny/2) :: vaux
 integer, dimension(nasymm) :: asymm 
 real(rkind) :: yy,rho,dudyw,rhow
 real(rkind) :: uwall
 real(rkind) :: qw,dtdyw,Ttau,ystar,Mtau
 real(rkind), dimension(ny/2) :: ufav,vfav,tfav,ydist,ybig,duvddy,dutdy,duvdy,duhdy,dugdy
 real(rkind) :: uu,vv,rmuw,rnuw,retau,deltav,utau,tauw,ttw,dtw,duw,dyw
 real(rkind) :: pp,tt,tt2,pp2,rho2,u2p,v2p,w2p,uvp
 real(rkind) :: aad,bbd,ccd,aai,bbi,cci
 real(rkind) :: intj,intjm,integrand
 real(rkind) :: seq,seqp,setl,setlp
 real(rkind) :: rhormsp,trmsp,prmsp
 real(rkind) :: taup
 integer :: i,j,m
 character(20) :: tname
!
 call read_grid_cha
!
 do j=1,ny/2
  ydist(j) = 1._rkind+y(j)
 enddo
!
 wstat1d = 0._rkind
 do j=1,ny
  do i=1,nx
   do m=1,nv
    wstat1d(m,j) = wstat1d(m,j) + wstat(i,j,m)
   enddo
  enddo
 enddo
 wstat1d = wstat1d/nx
!
 asymm = (/3,14,19,31,36,38,40,42,44,50,52,54,55,60,64/) 
 call avg_walls(ny,nv,nasymm,asymm,wstat1d,w1dh)
! 
 do j=1,ny/2
  ufav(j) = w1dh(13,j)/w1dh(1,j)
  vfav(j) = w1dh(14,j)/w1dh(1,j)
  tfav(j) = w1dh(25,j)/w1dh(1,j)
 enddo
!
 aad =  9._rkind/8._rkind
 bbd = -1._rkind/8._rkind
 ccd =  0._rkind/8._rkind
 duw = 2._rkind*aad*ufav(1)+2._rkind/3._rkind*bbd*ufav(2)+2._rkind/5._rkind*ccd*ufav(3)
 dyw = 2._rkind*aad*ydist(1)+2._rkind/3._rkind*bbd*ydist(2)+2._rkind/5._rkind*ccd*ydist(3)
 dtw = 2._rkind*aad*tfav(1)+2._rkind/3._rkind*bbd*tfav(2)+2._rkind/5._rkind*ccd*tfav(3)
 dtw = dtw-2._rkind*Twall*(aad+bbd/3._rkind+ccd/5._rkind)
 dudyw  = duw/dyw
 dtdyw  = dtw/dyw
!Coefficients from Lele pag 39
 aai = 150._rkind/128._rkind
 bbi = -25._rkind/128._rkind
 cci =   3._rkind/128._rkind
 rhow  = aai*w1dh(1,1)+bbi*w1dh(1,2)+cci*w1dh(1,3)
!ttw   = aai*w1dh(6 ,1)+bbi*w1dh(6 ,2)+cci*w1dh(6 ,3)
!rmuw  = aai*w1dh(20,1)+bbi*w1dh(20,2)+cci*w1dh(20,3)
 rmuw  = mu0
 rnuw  = rmuw/rhow
 ttw   = Twall ! overwrite
!
 tauw   = rmuw*dudyw
 utau   = sqrt(tauw/rhow)
 deltav = rnuw/utau
 retau  = 1._rkind/deltav
!
 qw     = -k0*dTdyw
 Ttau   = qw/(rhow*cp0*utau)
 Mtau   = utau/sqrt(gam*ttw)
!
! Compute compressibility transformations
!
 vaux(1,:) = w1dh(1 ,:)/rhow
 vaux(2,:) = w1dh(21,:)/rnuw
 vaux(3,:) = vaux(1,:)*vaux(2,:)
 vaux(4,:) = ydist/deltav ! yplus
 vaux(5,:) = ydist*sqrt(tauw*w1dh(1,:))/w1dh(20,:) ! ystar
 vaux(6,:) = utau/sqrt(gam*ttw) ! Mach_tau constant
 do j=1,ny/2
  ybig(j) = ydist(j)*sqrt(vaux(1,j))/vaux(3,j) ! when normalized by deltav corresponds to ystar
 enddo
!
 tname = 'vanDriest'
 call transform(1,ny/2,ydist,ufav,naux,vaux,yvd,uvd,tname)
 call ddy(ny/2,uvd(:),yvd,2,duvddy(:))
!
 tname = 'TrettelLarsson'
 call transform(1,ny/2,ydist,ufav,naux,vaux,yt,ut,tname)
 call ddy(ny/2,ut(:),yt,2,dutdy(:))
!
 tname = 'Volpiani'
 call transform(1,ny/2,ydist,ufav,naux,vaux,yv,uv,tname)
 call ddy(ny/2,uv(:),yv,2,duvdy(:))
!
 tname = 'Hasan'
 call transform(1,ny/2,ydist,ufav,naux,vaux,yh,uh,tname)
 call ddy(ny/2,uh(:),yv,2,duhdy(:))
!
 tname = 'Griffin'
 seq   = ufav(1)/ybig(1)
 seqp  = seq/utau*deltav
 setl  = ufav(1)/y(1)
 setlp = setl/utau*deltav
 integrand  = seq/(1._rkind+seqp-setlp)
 ug(1) = integrand*ybig(1)
 do j=2,ny/2
  taup  = 1._rkind-ydist(j)
  seq   = 1._rkind/vaux(3,j)*(ufav(j)-ufav(j-1))/(ybig(j)-ybig(j-1))
  seqp  = seq/utau*deltav
  setl  = vaux(3,j)*(ufav(j)-ufav(j-1))/(ydist(j)-ydist(j-1))
  setlp = setl/utau*deltav
  integrand  = taup*seq/(taup+seqp-setlp)
  ug(j) = ug(j-1)+integrand*(ybig(j)-ybig(j-1))
 enddo
 call ddy(ny/2,ug(:),yh,2,dugdy(:)) ! yg = yh
!
 open(10,file='POSTPRO/channinfo.dat',form='formatted')
 if (theta_wall>=-1._rkind) then
  write(10,*) 'Turbulent channel flow with constant bulk temperature and theta_wall =', theta_wall
  write(10,*) 'Mach number based on bulk temperature:', Mach_input
 else
  write(10,*) 'Turbulent channel flow with free bulk temperature'
 endif
 write(10,*) 'Mach number based on wall temperature:', Mach
 write(10,*) 'Bulk Reynolds number based on wall viscosity:', Reynolds
 write(10,*) 'Friction Reynolds number:', retau
 write(10,*) 'Skin friction coefficient:', 2._rkind*tauw/(rho0*u0**2)
 write(10,*) 'Wall density / bulk density:', rhow/rho0
 write(10,*) 'Friction velocity/ bulk velocity:', utau/u0
 write(10,*) 'Friction temperature/ wall temperature:', ttau/Twall
 close(10)
!
 open(10,file='POSTPRO/channstat.prof',form='formatted')
 do j=1,ny/2
  yy    = ydist(j)
  ystar = vaux(5,j)
  rho   = w1dh(1,j) 
  uu    = ufav(j)
  vv    = vfav(j)
  tt    = tfav(j)
  pp    = w1dh(5,j) 
  duvddy(j) = yvd(j)*duvddy(j)/utau
  dutdy(j)  = yt(j)*dutdy(j)/utau
  duvdy(j)  = yv(j)*duvdy(j)/utau
  duhdy(j)  = yh(j)*duhdy(j)/utau
  dugdy(j)  = yh(j)*dugdy(j)/utau
  u2p = (w1dh(16,j)-w1dh(13,j)*w1dh(13,j)/w1dh(1,j))/tauw
  v2p = (w1dh(17,j)-w1dh(14,j)*w1dh(14,j)/w1dh(1,j))/tauw
  w2p = (w1dh(18,j)-w1dh(15,j)*w1dh(15,j)/w1dh(1,j))/tauw
  uvp = (w1dh(19,j)-w1dh(13,j)*w1dh(14,j)/w1dh(1,j))/tauw
  rho2  = w1dh(7,j)  - rho**2
  pp2   = w1dh(11,j) - pp**2 
  tt2   = w1dh(12,j) - tt**2 
  rhormsp = sqrt(abs(w1dh( 7,j)-w1dh(1,j)**2))/(rhow*gam*Mtau**2)
  trmsp   = sqrt(abs(w1dh(26,j)/w1dh(1,j)-tfav(j)**2))/(ttw*gam*Mtau**2)
  prmsp   = sqrt(abs(w1dh(11,j)-w1dh(5,j)**2))/tauw

  write(10,100) yy,              &         !1  y/h
                rho/rho0,        &         !2  rho/rho bulk
                uu/u0,           &         !3  uu/u bulk
                tt/t0,           &         !4  tt/Twall
                pp/rho0/u0**2,   &         !5  pp/rho bulk/u bulk**2
                w1dh(20,j)/rmuw, &         !6  mu/mu wall
                yy/deltav,       &         !7  y+
                ystar,           &         !8  y*
                yv(j)/deltav,    &         !9  y_V^+
                uu/utau,         &         !10 u^+
                uvd(j)/utau,     &         !11 u_VD^+
                ut(j)/utau,      &         !12 u_TL^+
                uv(j)/utau,      &         !13 u_v^+
                ug(j)/utau,      &         !14 u_g^+
                uh(j)/utau,      &         !15 u_h^+
                duvddy(j),       &         !16
                dutdy(j),        &         !17
                duvdy(j),        &         !18
                dugdy(j),        &         !19
                duhdy(j),        &         !20
                u2p,             &         !21
                v2p,             &         !22
                w2p,             &         !23
                uvp,             &         !24
                rhormsp,         &         !25
                trmsp,           &         !26
                prmsp                      !27
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
