module postpro_chacurv
 use parameters
 use global_variables 
 use reader
 use derivatives
 use comp_transform

 contains 

 subroutine stats1d_curv
 implicit none
 integer, parameter :: naux   = 6  ! Number of auxiliary variables for compressibility transformations 
 real(rkind), dimension(ny)   :: u, dudy
 real(rkind), dimension(:,:), allocatable :: wstat1d
 real(rkind), dimension(:,:), allocatable :: w1dh
 real(rkind), dimension(1-ng:ny/2+ng) :: yint,yext
 real(rkind), dimension(ny/2) :: yvd,yt,yv,ut,uvd,uv,yh,uh
 real(rkind), dimension(naux,ny/2) :: vaux
 real(rkind) :: yy,rho,dudyw,rhowall,drhodyw,d2rhodyw
 real(rkind) :: rhowall_i,rhowall_o,muwall_i,muwall_o,utau_i,utau_o
 real(rkind) :: uwall,d2udyw
 real(rkind) :: qw,dTdyw,d2Tdyw,Ttau
 real(rkind) :: ufav(ny/2),vfav(ny/2),wfav(ny/2),tfav(ny/2),uuu(ny)
 real(rkind), dimension(ny) :: ufav_tot,vfav_tot,wfav_tot,tfav_tot,uuu_tot,ulam
 real(rkind) :: uu,muwall,nuwall,retau,deltav,utau,tauw
 real(rkind) :: pp,tt,tt2,pp2,rho2,rhou2,rhov2,rhow2,rhouv,rhow,rmuw,ttw,rnuw
 real(rkind) :: rhormsp,prmsp,trmsp,dtw,duw,dyw,ystar,Mtau,vv
 real(rkind) :: aad,bbd,ccd,aai,bbi,cci
 real(rkind) :: ri2, ro2, rc2
 real(rkind) :: aa,alfa,bb,beta,rr,ubulk,umax1,umax2,umaxlam,umean,usum,uu2,vlam,vpoi,vv2,ww2
 integer :: i,j,m
 character(20) :: tname
 real(rkind) :: d1, d2, x10, xc1, b1, c1, u1, dudyw_streams1, dudyw_simple

 allocate(wstat1d(nv,ny))
 allocate(w1dh(nv,ny/2))

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
!-----------------------------------------------------------------------------------
! CONVEX
!-----------------------------------------------------------------------------------
 w1dh(:,:) = wstat1d(:,1:ny/2)
 yint(1-ng:ny/2+ng) = y(1-ng:ny/2+ng) - R_curv + 1._rkind
 print*,'yint: ',yint
 !print*,'yint fix: ',yint + R_curv/2.

 do j=1,ny/2
  ufav(j) = w1dh(13,j)/w1dh(1,j)
  vfav(j) = w1dh(14,j)/w1dh(1,j)
  wfav(j) = w1dh(15,j)/w1dh(1,j)
  tfav(j) = w1dh(25,j)/w1dh(1,j)
 enddo
!
 aad =  9._rkind/8._rkind
 bbd = -1._rkind/8._rkind
 ccd =  0._rkind/8._rkind
 duw = 2._rkind*aad*ufav(1)+2._rkind/3._rkind*bbd*ufav(2)+2._rkind/5._rkind*ccd*ufav(3)
 dyw = 2._rkind*aad*yint(1)+2._rkind/3._rkind*bbd*yint(2)+2._rkind/5._rkind*ccd*yint(3)
 dtw = 2._rkind*aad*tfav(1)+2._rkind/3._rkind*bbd*tfav(2)+2._rkind/5._rkind*ccd*tfav(3)
 print*,'dtw prima,Twall: ',dtw,Twall
 dtw = dtw-2._rkind*Twall*(aad+bbd/3._rkind+ccd/5._rkind)
 print*,'dtw, dyw: ',dtw,dyw
 dudyw  = duw/dyw
 dtdyw  = dtw/dyw
!Coefficients from Lele pag 39
 aai = 150._rkind/128._rkind
 bbi = -25._rkind/128._rkind
 cci =   3._rkind/128._rkind
 rhow  = aai*w1dh(1,1)+bbi*w1dh(1,2)+cci*w1dh(1,3)
!ttw   = aai*w1dh(6 ,1)+bbi*w1dh(6 ,2)+cci*w1dh(6 ,3)
!rmuw  = aai*w1dh(20,1)+bbi*w1dh(20,2)+cci*w1dh(20,3)
!
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
 rhowall_i = rhow
 muwall_i  = rmuw
 utau_i    = utau
!
 print*, 'Convex side'
 print*,'*************************************************************************'
 print*, 'rho_w =',rhow, '   ','mu_w      =',rmuw
 print*, 'dudyw =',dudyw,'   ','dTdyw     =',dtdyw
 print*, 'tau_w =',tauw, '   ','u_tau     =',utau
 print*, 'retau =',retau,'   ','u+ (ny/2) =',ufav(ny/2)/utau
 print*,'*************************************************************************'
!
! Compute compressibility transformations
!
 vaux(1,:) = w1dh(1 ,:)/rhow
 vaux(2,:) = w1dh(21,:)/rnuw
 vaux(3,:) = vaux(1,:)*vaux(2,:)
 vaux(4,:) = yint/deltav ! yplus
 vaux(5,:) = yint*sqrt(tauw*w1dh(1,:))/w1dh(20,:) ! ystar
 vaux(6,:) = utau/sqrt(gam*ttw) ! Mach_tau constant
!
 tname = 'vanDriest'
 call transform(1,ny/2,yint,ufav,naux,vaux,yvd,uvd,tname)
!
 tname = 'TrettelLarsson'
 call transform(1,ny/2,yint,ufav,naux,vaux,yt,ut,tname)
!
 tname = 'Volpiani'
 call transform(1,ny/2,yint,ufav,naux,vaux,yv,uv,tname)
!
 tname = 'Hasan'
 call transform(1,ny/2,yint,ufav,naux,vaux,yh,uh,tname)
!
 open(10,file='POSTPRO/channinfo_convex.dat',form='formatted')
 write(10,*)'Retau,          rho_wall,         utau/u0,         Cf'
 write(10,*)retau,           rhow,             utau/u0,  2._rkind*tauw/(rho0*u0**2)
 close(10)
!
 open(10,file='POSTPRO/channstat_convex.prof',form='formatted')
 do j=1,ny/2
  yy   = yint(j)+1.
  ystar = vaux(5,j)
  rho  = w1dh(1,j) 
  uu   = ufav(j)
  vv   = vfav(j)
  tt   = tfav(j)
  pp   = w1dh(5,j) 
  rhou2  = w1dh(16,j) - w1dh(13,j)**2/w1dh(1,j) 
  rhov2  = w1dh(17,j) - w1dh(14,j)**2/w1dh(1,j)
  rhow2  = w1dh(18,j) - w1dh(15,j)**2/w1dh(1,j) 
  rhouv  = w1dh(19,j) - w1dh(13,j)*w1dh(14,j)/w1dh(1,j)
  rho2   = w1dh(7,j)  - rho**2
  pp2    = w1dh(11,j) - pp**2 
  tt2    = w1dh(12,j) - tt**2 
  rhormsp = sqrt(abs(w1dh( 7,j)-w1dh(1,j)**2))/(rhow*gam*Mtau**2)
  trmsp   = sqrt(abs(w1dh(26,j)/w1dh(1,j)-tfav(j)**2))/(ttw*gam*Mtau**2)
  prmsp   = sqrt(abs(w1dh(11,j)-w1dh(5,j)**2))/tauw
  write(10,100)yy,              &         !1  y/h
               rho/rho0,        &         !2  rho/rho bulk
               uu/u0,           &         !3  u/u_b
               tt/t0,           &         !4  tt/Twall
               pp/rho0/u0**2,   &         !5  pp/rho bulk/u bulk**2
               w1dh(20,j)/rmuw, &         !6  mu/mu wall
               yy/deltav,       &         !7  y^+
               ystar,           &         !8  y*
               yv(j)/deltav,    &         !9  y_V^+
               uu/utau,         &         !10 u^+
               uvd(j)/utau,     &         !11 u_VD^+
               ut(j)/utau,      &         !12 u_TL^+
               uv(j)/utau,      &         !13 u_v^+
               uh(j)/utau,      &         !14 u_h^+
               rhou2/tauw,      &         !15 rho/rhowall*u2^+   
               rhov2/tauw,      &         !16 rho/rhowall*v2^+  
               rhow2/tauw,      &         !17 rho/rhowall*w2^+ 
               rhouv/tauw,      &         !18 rho/rhowall*uv^+
               rho2/rhowall,    &         !19 rho/rhowall
               rhormsp,         &         !20
               trmsp,           &         !21
               prmsp                      !22
 enddo
 close(10)

!-----------------------------------------------------------------------------------
! CONCAVE
!-----------------------------------------------------------------------------------
 w1dh(:,:) = wstat1d(:,ny:ny/2+1:-1)
 yext(1-ng:ny/2+ng) = -y(ny+ng:ny/2-ng:-1)+R_curv

 do j=1,ny/2
  ufav(j) = w1dh(13,j)/w1dh(1,j)
  vfav(j) = w1dh(14,j)/w1dh(1,j)
  wfav(j) = w1dh(15,j)/w1dh(1,j)
  tfav(j) = w1dh(25,j)/w1dh(1,j)
 enddo
!
 aad =  9._rkind/8._rkind
 bbd = -1._rkind/8._rkind
 ccd =  0._rkind/8._rkind
 duw = 2._rkind*aad*ufav(1)+2._rkind/3._rkind*bbd*ufav(2)+2._rkind/5._rkind*ccd*ufav(3)
 dyw = 2._rkind*aad*yint(1)+2._rkind/3._rkind*bbd*yint(2)+2._rkind/5._rkind*ccd*yint(3)
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
!
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
 rhowall_o = rhow
 muwall_o  = rmuw
 utau_o    = utau
!
 print*, 'Concave side'
 print*,'*************************************************************************'
 print*, 'rho_w =',rhow, '   ','mu_w      =',rmuw
 print*, 'dudyw =',dudyw,'   ','dTdyw     =',dtdyw
 print*, 'tau_w =',tauw, '   ','u_tau     =',utau
 print*, 'retau =',retau,'   ','u+ (ny/2) =',ufav(ny/2)/utau
 print*,'*************************************************************************'

 vaux(1,:) = w1dh(1 ,:)/rhow
 vaux(2,:) = w1dh(21,:)/rnuw
 vaux(3,:) = vaux(1,:)*vaux(2,:)
 vaux(4,:) = yext/deltav ! yplus
 vaux(5,:) = yext*sqrt(tauw*w1dh(1,:))/w1dh(20,:) ! ystar
 vaux(6,:) = utau/sqrt(gam*ttw) ! Mach_tau constant
!
 tname = 'vanDriest'
 call transform(1,ny/2,yext,ufav,naux,vaux,yvd,uvd,tname)
!
 tname = 'TrettelLarsson'
 call transform(1,ny/2,yext,ufav,naux,vaux,yt,ut,tname)
!
 tname = 'Volpiani'
 call transform(1,ny/2,yext,ufav,naux,vaux,yv,uv,tname)
!
 tname = 'Hasan'
 call transform(1,ny/2,yext,ufav,naux,vaux,yh,uh,tname)
!
 open(10,file='POSTPRO/channinfo_concave.dat',form='formatted')
 write(10,*)'Retau,          rho_wall,         utau/u0,         Cf'
 write(10,*)retau,           rhow,             utau/u0,  2._rkind*tauw/(rho0*u0**2)
 close(10)

 open(10,file='POSTPRO/channstat_concave.prof',form='formatted')
 do j=1,ny/2
  yy   = yext(j)+1.
  ystar = vaux(5,j)
  rho  = w1dh(1,j)
  uu   = ufav(j)
  vv   = vfav(j)
  tt   = tfav(j)
  pp   = w1dh(5,j)
  rhou2  = w1dh(16,j) - w1dh(13,j)**2/w1dh(1,j)
  rhov2  = w1dh(17,j) - w1dh(14,j)**2/w1dh(1,j)
  rhow2  = w1dh(18,j) - w1dh(15,j)**2/w1dh(1,j)
  rhouv  = w1dh(19,j) - w1dh(13,j)*w1dh(14,j)/w1dh(1,j)
  rho2   = w1dh(7,j)  - rho**2
  pp2    = w1dh(11,j) - pp**2
  tt2    = w1dh(12,j) - tt**2
  rhormsp = sqrt(abs(w1dh( 7,j)-w1dh(1,j)**2))/(rhow*gam*Mtau**2)
  trmsp   = sqrt(abs(w1dh(26,j)/w1dh(1,j)-tfav(j)**2))/(ttw*gam*Mtau**2)
  prmsp   = sqrt(abs(w1dh(11,j)-w1dh(5,j)**2))/tauw
  write(10,100)yy,              &         !1  y/h
               rho/rho0,        &         !2  rho/rho bulk
               uu/u0,           &         !3  u/u_b
               tt/t0,           &         !4  tt/Twall
               pp/rho0/u0**2,   &         !5  pp/rho bulk/u bulk**2
               w1dh(20,j)/rmuw, &         !6  mu/mu wall
               yy/deltav,       &         !7  y^+
               ystar,           &         !8  y*
               yv(j)/deltav,    &         !9  y_V^+
               uu/utau,         &         !10 u^+
               uvd(j)/utau,     &         !11 u_VD^+
               ut(j)/utau,      &         !12 u_TL^+
               uv(j)/utau,      &         !13 u_v^+
               uh(j)/utau,      &         !14 u_h^+
               rhou2/tauw,      &         !15 rho/rhowall*u2^+
               rhov2/tauw,      &         !16 rho/rhowall*v2^+
               rhow2/tauw,      &         !17 rho/rhowall*w2^+
               rhouv/tauw,      &         !18 rho/rhowall*uv^+
               rho2/rhowall,    &         !19 rho/rhowall
               rhormsp,         &         !20
               trmsp,           &         !21
               prmsp                      !22
 enddo

!-----------------------------------------------------------------------------------
! GLOBAL
!-----------------------------------------------------------------------------------
  ri2 = (R_curv-1._rkind)*(R_curv-1._rkind)
  ro2 = (R_curv+1._rkind)*(R_curv+1._rkind)
  rc2 = (R_curv*R_curv)*2._rkind

  do j=1,ny
   uuu_tot(j) = wstat1d(2,j)
   ufav_tot(j) = wstat1d(13,j)/wstat1d(1,j)
   vfav_tot(j) = wstat1d(14,j)/wstat1d(1,j)
   wfav_tot(j) = wstat1d(15,j)/wstat1d(1,j)
   tfav_tot(j) = wstat1d(25,j)/wstat1d(1,j)
  enddo

  rhowall = 0.5_rkind*(rhowall_i+rhowall_o)
  muwall  = 0.5_rkind*(muwall_i+muwall_o)
  utau    = sqrt((ri2*(utau_i*utau_i)+ro2*(utau_o*utau_o))/rc2)
  retau   = rhowall * utau / muwall
  print*, 'Global re_tau =', retau
  umax1  = maxval(ufav_tot)/u0
  umax2  = maxval(uuu)/u0
  aa = R_curv - 1._rkind
  bb = R_curv + 1._rkind
  alfa = (bb**2*log(bb)-aa**2*log(aa))/(bb**2-aa**2)
  beta = (aa**2*bb**2) * log(bb/aa)   /(bb**2-aa**2)
  umean = (0.25_rkind*(bb**2-aa**2)-(aa**2*bb**2)*(log(bb/aa))**2/(bb**2-aa**2))/(bb-aa)
  usum = 0._rkind
  do j=1,ny
    rr = yn(j) + R_curv       
    ulam(j) = alfa*rr-beta/rr-rr*log(rr)
    usum = usum + ufav_tot(j)
  enddo
  umaxlam  = maxval(ulam)/umean
  ubulk = usum/ny
  !print*, 'u_bulk', u0, ubulk
  !print*, 'u_max =', maxval(ufav_tot), 'u_center =', ufav_tot(ny/2)
!
  open(unit=10,file='POSTPRO/channstat_global.prof',form='formatted')
  do j=1,ny
   yy   = y(j) - R_curv
   vpoi = 1.5_rkind*(1._rkind-yy*yy)
   !vlam = ulam(j)/umaxlam 
   vlam = (alfa*y(j)-beta/y(j)-y(j)*log(y(j)))
   uu2  = wstat1d(16,j)/wstat1d(1,j) - ufav_tot(j)*ufav_tot(j)
   vv2  = wstat1d(17,j)/wstat1d(1,j) - vfav_tot(j)*vfav_tot(j)
   ww2  = wstat1d(18,j)/wstat1d(1,j) - wfav_tot(j)*wfav_tot(j)
   uv   = wstat1d(19,j)/wstat1d(1,j) - ufav_tot(j)*vfav_tot(j)
   !tauv = w_av(20,j)*w_av(21,j)
   !omz2 = w_av(25,j)-w_av(22,j)**2
   !omx2 = w_av(26,j)-w_av(23,j)**2
   !omy2 = w_av(27,j)-w_av(24,j)**2
!
   write(10,100) yy, &                !y/h
                 !vpoi,&
                 !vlam/umean,&
                 !ulam(j)/umean,&
                 ufav_tot(j)/u0,&         !u/u_b
                 uuu(j)/u0,&          !u/u_b
                 !ufav_tot(j)/ubulk,&     !u/u_b
                 uu2/utau**2,&        !u2^+
                 vv2/utau**2,&        !v2^+
                 ww2/utau**2,&        !w2^+
                 uv /utau**2!,&       !uv^+
                 !tauv/tauw,&         !tau_visc
                 !omz2,&
                 !omx2,&
                 !omy2
  enddo
  close(10)

 100  format(200ES20.10)
 end subroutine stats1d_curv

end module postpro_chacurv
