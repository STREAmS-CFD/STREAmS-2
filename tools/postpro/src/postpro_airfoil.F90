module postpro_airfoil
 use parameters
 use global_variables 
 use reader
 use derivatives
 use comp_transform

 contains 

 subroutine forces_airfoil
 implicit none

 integer     :: i
 real(rkind) :: cfs,costh,sinth,dudy,dudy1,dudyw,dudyw1,ut1,ut2,ut3,ut4
 real(rkind) :: al,ds,dN,dA,dpN,dpA,dtN,dtA,N,A,pN,pA,tN,tA
 real(rkind) :: lift,drag,pres,fric,pdyn,tauw,tauw1,tauwf,cf,cp,pw,pw1,pwf

 call read_grid_airfoil_ramp()

 open(unit=12,file='POSTPRO/avg_coeff.dat',form='formatted')
  N  = 0._rkind
  A  = 0._rkind
  pN = 0._rkind
  pA = 0._rkind
  tN = 0._rkind
  tA = 0._rkind
  do i=ite,itu
   ds    = .5_rkind*(csimod(i,1)+csimod(i+1,1))  ! csimod  id the module of the wall-tangent
   costh = .5_rkind*(dxdcsi(i,1)+dxdcsi(i+1,1))/ds ! dxdcsin is the x-component of the wall-ta
   sinth = .5_rkind*(dydcsi(i,1)+dydcsi(i+1,1))/ds ! dydcsin is the y-component of the wall-ta

   ut1    = wstat(i,1,2)*costh+wstat(i,1,3)*sinth
   ut2    = wstat(i,2,2)*costh+wstat(i,2,3)*sinth
   ut3    = wstat(i,3,2)*costh+wstat(i,3,3)*sinth
   ut4    = wstat(i,4,2)*costh+wstat(i,4,3)*sinth
   dudy   = -22._rkind*ut1+36._rkind*ut2-18._rkind*ut3+4._rkind*ut4
   dudyw  = dudy*meta(i,1)/12._rkind
   tauw   = wstat(i,1,20)*dudyw
   pw     = wstat(i,1,5)

   ut1    = wstat(i+1,1,2)*costh+wstat(i+1,1,3)*sinth
   ut2    = wstat(i+1,2,2)*costh+wstat(i+1,2,3)*sinth
   ut3    = wstat(i+1,3,2)*costh+wstat(i+1,3,3)*sinth
   ut4    = wstat(i+1,4,2)*costh+wstat(i+1,4,3)*sinth
   dudy1  = -22._rkind*ut1+36._rkind*ut2-18._rkind*ut3+4._rkind*ut4
   dudyw1 = dudy1*meta(i+1,1)/12._rkind
   tauw1  = wstat(i+1,1,20)*dudyw1
   pw1    = wstat(i+1,1,5)

   ! cf without sign
   cf = tauw/(0.5_rkind*u0*u0)
   ! cf with sign
   if (i< ile) cfs =-tauw/(0.5_rkind*u0*u0)
   if (i>=ile) cfs = tauw/(0.5_rkind*u0*u0)

   cp = (pw-p0)/(0.5_rkind*u0*u0)
   pwf   = .5_rkind*(pw+pw1)
   tauwf = .5_rkind*(tauw+tauw1)
   if (i<itu) then
    dN  = -costh*pwf+sinth*tauwf
    dA  = +sinth*pwf+costh*tauwf
    dpN = -costh*pwf
    dpA = +sinth*pwf
    dtN = sinth*tauwf
    dtA = costh*tauwf
    N  = N + dN*ds
    A  = A + dA*ds
    pN = pN + dpN*ds
    pA = pA + dpA*ds
    tN = tN + dtN*ds
    tA = tA + dtA*ds
   endif
   write(12,100) xg(i,1),cf,cfs,cp,pw,tauw
  enddo
 close(12)
!
 al   = aoa*pi/180._rkind
 pdyn = .5_rkind*u0*u0 
 lift = ( N*cos(al)- A*sin(al))/pdyn 
 drag = ( N*sin(al)+ A*cos(al))/pdyn
 pres = (pN*sin(al)+pA*cos(al))/pdyn
 fric = (tN*sin(al)+tA*cos(al))/pdyn
 print*, 'Lift     =',lift,'Drag     =',drag
 print*, 'Pressure =',pres,'Friction =',fric
 open(unit=30,file='POSTPRO/avg_forces.dat')
  write(30,100) lift,drag,pres,fric
 close(30)
!  
 100  format(200ES20.10)
 end subroutine forces_airfoil

 subroutine stats2d_airfoil
 implicit none
 integer, parameter :: naux = 6
 integer :: i,j,m,j99,ii,jvort
 real(rkind), dimension(3,3) :: sig
 real(rkind), dimension(npoints_bl) :: ufav_bl,vfav_bl,wfav_bl,tfav_bl,utfav_bl,vtfav_bl,rhout2_bl,rhovt2_bl,rhoutvt_bl
 real(rkind), dimension(nx) :: rmutil,d99_vec,deltav_vec,utau_vec,rhow_vec,ttw_vec,dudyw_vec,rmuw_vec,jvortical
 real(rkind), dimension(naux,npoints_bl) :: vaux
 real(rkind), dimension(npoints_bl) :: yt,yvd,yv,uv,ut,uvd,yh,uh,uinv,dutfav_bl_dy
 real(rkind) :: uwall,dudyw,d2udyw,rmuw,rnuw,tauw,rhow,ttw,deltav,pdyn,dthre
 real(rkind) :: cf,dely,utau,uu,unum,uden,delta99,retau,prms,rmu,chi3,fv1,epsi 
 real(rkind) :: dstarinc,thetainc,rhoe,pe,ue,tte,sqgmr,cv13,rnutil,rnut 
 real(rkind) :: rho,rhop,up,dy,dyh,dyw,dstar,theta,uup,gg,ry,py,udel,d_bl,uref,vref,al
 real(rkind) :: shapef,shapefinc,rmue,ff,aa,bb,alpha,beta,prod
 real(rkind) :: ucsi,ueta,vcsi,veta,ux,uy,uz,vx,vy,vz,wx,wy,wz,eps11,eps22,eps33,eps12
 real(rkind) :: fc,ftheta,cfinc,rethetainc,rethetawall,redelta99,redelta2,retheta
 real(rkind) :: y99,uvdp,u2p,v2p,w2p,uvp,rhofac,utp,uvpl,yvp,ytp,Mtau,mum,deltaj99,ch,ppm,ppw,qw
 real(rkind) :: yp,cpcoeff,pwall,yyy,Bq,ttm,trmsp,uhp,ystar,rlamw,trec,rhormsp,rhom
 real(rkind) :: dtdyw,qwall,vortx2,vorty2,vortz2,omag,omag_factor,eta
 character(20) :: tname
 character(6)  :: chstat
 logical filecheck
!
! Mean boundary layer properties
!
  udel = 0.99_rkind*u0 
  gg    = 2*gam/(gam-1)
  al    = aoa*pi/180._rkind
  uref  = u0*cos(al)
  vref  = u0*sin(al)
  sqgmr = u0/Reynolds
  cv13  = 7.1**3
!
  inquire(file='eddy.dat', exist=filecheck)
  if (filecheck) then
     open(11,file='eddy.dat',form='formatted')
     do i=1,nx
      read(11,*) rmutil(i)
     enddo
     close(11)

     open(13,file='POSTPRO/wake_rans.dat',form='formatted')
     do i=itu+1,nx,ix_out
      call extract_bl(i)
      do j=1,4
          ufav_bl(j)  = wstat_bl(j,13)/wstat_bl(j,1)
          vfav_bl(j)  = wstat_bl(j,14)/wstat_bl(j,1)
          utfav_bl(j) = (ufav_bl(j)*dxdcsi(i,1)+vfav_bl(j)*dydcsi(i,1))/csimod(i,1) ! csi tangent wall
          vtfav_bl(j) = (ufav_bl(j)*detadx(i,1)+vfav_bl(j)*detady(i,1))/meta(i,1) ! csi normal wall
      enddo
      dudyw  = -22._rkind*utfav_bl(1)+36._rkind*utfav_bl(2)-18._rkind*utfav_bl(3)+ 4._rkind*utfav_bl(4)
      dyw    = -22._rkind*0.000_rkind+36._rkind*delta_bl(1)-18._rkind*delta_bl(2)+ 4._rkind*delta_bl(3)
      dudyw  = dudyw / dyw

      rnutil = rmutil(i)/wstat_bl(1,1)
      rnuw   = wstat_bl(1,20)/wstat_bl(1,1)
      chi3   = (rnutil/rnuw)**3
      fv1    = chi3/(chi3+cv13)
      rnut   = rnutil*fv1
      epsi   = rnut*dudyw**2 ! dissipation ≈ production
      eta    =(rnuw**3/epsi)**0.25
      write(13,100) xg(i,1),eta,rnuw,epsi,rnut,dudyw
     enddo
     close(13)
  endif

  open(13,file='POSTPRO/wake.dat',form='formatted')
  do i=itu+1,nx,ix_out
   call extract_bl(i)
   do j=1,4
       ufav_bl(j)    = wstat_bl(j,13)/wstat_bl(j,1)
       vfav_bl(j)    = wstat_bl(j,14)/wstat_bl(j,1)
       wfav_bl(j)    = wstat_bl(j,15)/wstat_bl(j,1)
       tfav_bl(j)    = wstat_bl(j,25)/wstat_bl(j,1)
       utfav_bl(j)   = (ufav_bl(j)*dxdcsi(i,1)+vfav_bl(j)*dydcsi(i,1))/csimod(i,1) ! csi tangent wall
       vtfav_bl(j)   = (ufav_bl(j)*detadx(i,1)+vfav_bl(j)*detady(i,1))/meta(i,1)   ! csi normal
       rhout2_bl(j)  = (wstat_bl(j,16)*dxdcsi(i,1)**2+wstat_bl(j,17)*dydcsi(i,1)**2+ &
                       2._rkind*wstat_bl(j,19)*dxdcsi(i,1)*dydcsi(i,1))/csimod(i,1)**2
       rhovt2_bl(j)  = (wstat_bl(j,16)*detadx(i,1)**2+wstat_bl(j,17)*detady(i,1)**2+ &
                       2._rkind*wstat_bl(j,19)*detadx(i,1)*detady(i,1))/meta(i,1)**2
       rhoutvt_bl(j) = (wstat_bl(j,16)*dxdcsi(i,1)*detadx(i,1)+wstat_bl(j,17)*dydcsi(i,1)*detady(i,1)+ &
                       wstat_bl(j,19)*(dxdcsi(i,1)*detady(i,1)+dydcsi(i,1)*detadx(i,1)))/ &
                       csimod(i,1)/meta(i,1)
   enddo

   dudyw  = -22._rkind*utfav_bl(1)+36._rkind*utfav_bl(2)-18._rkind*utfav_bl(3)+ 4._rkind*utfav_bl(4)
   dyw    = -22._rkind*0.000_rkind+36._rkind*delta_bl(1)-18._rkind*delta_bl(2)+ 4._rkind*delta_bl(3)
   dudyw  = dudyw / dyw

   uvp  = rhoutvt_bl(1) -wstat_bl(1,1)*utfav_bl(1)*vtfav_bl(1)
   rnuw = wstat_bl(1,20)/wstat_bl(1,1)
   prod = abs(uvp*dudyw) 
   epsi = wstat_bl(1,57) + wstat_bl(1,58) + wstat_bl(1,59) + wstat_bl(1,60)

   ucsi = 0.5_rkind*(wstat(i+1,1,2)-wstat(i-1,1,2)) ! using wstat since j=1
   vcsi = 0.5_rkind*(wstat(i+1,1,3)-wstat(i-1,1,3)) ! using wstat since j=1
   ueta =-1.5_rkind*wstat_bl(1,2)+2._rkind*wstat_bl(2,2)-0.5_rkind*wstat_bl(3,2)
   veta =-1.5_rkind*wstat_bl(1,3)+2._rkind*wstat_bl(2,3)-0.5_rkind*wstat_bl(3,3)
   ux = ucsi*dcsidx(i,1) + ueta*detadx(i,1)
   vx = vcsi*dcsidx(i,1) + veta*detadx(i,1)
   uy = ucsi*dcsidy(i,1) + ueta*detady(i,1)
   vy = vcsi*dcsidy(i,1) + veta*detady(i,1)

   sig(1,1) = wstat_bl(1,43)
   sig(1,2) = wstat_bl(1,44)
   sig(1,3) = wstat_bl(1,45)
   sig(2,2) = wstat_bl(1,46)
   sig(2,3) = wstat_bl(1,47)
   sig(3,3) = wstat_bl(1,48)

   wx=0.; wy=0.; wz=0.; uz=0.; vz=0.;
   eps11 = sig(1,1)*ux+sig(1,2)*uy+sig(1,3)*uz
   eps22 = sig(2,1)*vx+sig(2,2)*vy+sig(2,3)*vz
   eps33 = sig(3,1)*wx+sig(3,2)*wy+sig(3,3)*wz 
   eps12 = sig(1,1)*vx+sig(1,2)*(ux+vy)+sig(2,2)*uy+sig(1,3)*vz+sig(2,3)*uz

   epsi = epsi - eps11 - eps22 - eps33 - eps12

   write(13,100) xg(i,1),yg(i,1),epsi,(rnuw**3/epsi)**0.25
  enddo
  close(13)

  open(10,file='POSTPRO/bl_pressure.dat',form='formatted')
  do i=ite+1,ile-1,ix_out
   print*,'Analyzing wall i-index: ',i
   ! get wstat_bl(j,iv)
   call extract_bl(i)

   do j=1,npoints_bl
      ufav_bl(j)  = wstat_bl(j,13)/wstat_bl(j,1)
      vfav_bl(j)  = wstat_bl(j,14)/wstat_bl(j,1)
      wfav_bl(j)  = wstat_bl(j,15)/wstat_bl(j,1)
      tfav_bl(j)  = wstat_bl(j,25)/wstat_bl(j,1)
      utfav_bl(j) = (ufav_bl(j)*dxdcsi(i,1)+vfav_bl(j)*dydcsi(i,1))/csimod(i,1) ! csi tangent wall
      vtfav_bl(j) = (ufav_bl(j)*detadx(i,1)+vfav_bl(j)*detady(i,1))/meta(i,1)   ! csi normal wall
   enddo

   dudyw  = -22._rkind*utfav_bl(1)+36._rkind*utfav_bl(2)-18._rkind*utfav_bl(3)+4._rkind*utfav_bl(4)
   dyw    = -22._rkind*0.000_rkind+36._rkind*delta_bl(1)-18._rkind*delta_bl(2)+4._rkind*delta_bl(3)
   dtdyw  = -22._rkind* tfav_bl(1)+36._rkind* tfav_bl(2)-18._rkind* tfav_bl(3)+4._rkind* tfav_bl(4)
   dudyw  = dudyw / dyw
   dtdyw  = dtdyw / dyw

   rhow   = wstat_bl(1,1)
   ttw    = tfav_bl(1)
   rmuw   = wstat_bl(1,20)
   rnuw   = rmuw/rhow

   ttw    = Twall
   rlamw  = rmuw*gam/(gam-1._rkind)/Prandtl
   tauw   = rmuw*dudyw
   qw     = -rlamw*dtdyw
   utau   = sqrt(abs(tauw)/rhow)
   deltav = rnuw/utau
   pdyn   = 0.5_rkind*rho0*u0**2
   cf     = tauw/pdyn
   trec   = t0 * (1.+0.5*(gam-1.)*rfac*Mach**2)
   Bq     = qw/((gam/(gam-1.))*rhow*utau*ttw)
   if (theta_wall==1._rkind) then
    ch = 0._rkind
   else
    ch = qw/u0/(gam/(gam-1.))/(ttw-trec)
   endif

   prms  = sqrt(abs(wstat_bl(1,11)-wstat_bl(1,5)**2))/pdyn
   ppw   = rhow*ttw 

   cpcoeff = (ppw-p0)/pdyn

   !-----------------------------------------------------------
   ! VORTICAL CRTIERIUM
   !-----------------------------------------------------------
   !TEORICO omag_factor = 0.05_rkind
   omag_factor = 0.25_rkind
   !omag_factor = 0.95_rkind
   jvort = 1
   do j=1,npoints_bl-1
    vortx2 = wstat_bl(j,22)
    vorty2 = wstat_bl(j,23)
    vortz2 = wstat_bl(j,24)
    omag = sqrt(vortx2+vorty2+vortz2)
    !print*,'omag: ',omag,vortx2,vorty2,vortz2,omag_factor*u0/l0
    if (omag < omag_factor*u0/l0) then
     jvort = j
     exit
    endif
   enddo

   udel = 0.99_rkind*abs(utfav_bl(jvort))

   j99 = 1
   do j=1,npoints_bl-1
    if (abs(utfav_bl(j))>udel) then
     j99 = j-1
     exit
    endif
   enddo

   print*,'jvort, j99: ',jvort, j99
   !-----------------------------------------------------------
   !-----------------------------------------------------------
   ! INVISCID CRTIERIUM
   !-----------------------------------------------------------
   ! Computation of delta_99 (see Griffin et al, PRF 2021)
   !uinv=0._rkind
   !gg = 2*gam/(gam-1)
   !do j=1,npoints_bl-1
   ! ry=wstat_bl(j,1)
   ! py=wstat_bl(j,5)
   ! uinv(j)=sqrt(gg*(p0/rho0-py/ry)+uref**2+vref**2-vtfav_bl(j)**2) ! inviscid velocity profile
   !enddo
   !j99   = 1
   !do j=1,npoints_bl-1
   ! uu = utfav_bl(j)
   ! udel = 0.99_rkind*uinv(j)
   ! if (uu>udel) then
   !  j99 = j-1
   !  exit
   ! endif
   !enddo
   !jvort = j99
   !print*,'j: ',jvort, j99
   !-----------------------------------------------------------
   jvortical(i) = jvort

   deltaj99 = delta_bl(j99-1)
   dely = delta_bl(j99) - delta_bl(j99-1)
   unum = udel-abs(utfav_bl(j99))
   uden = abs(utfav_bl(j99+1))-abs(utfav_bl(j99))
   delta99 = deltaj99 + dely*(unum/uden) ! b._rkindl._rkind thickness

   retau = delta99/deltav
   d99_vec(i)    = delta99
   deltav_vec(i) = deltav
   utau_vec(i)   = utau
   rhow_vec(i)   = rhow
   ttw_vec(i)    = ttw
   rmuw_vec(i)   = rmuw
!
!  Integral boundary layer thicknesses
!
   dstar     = 0._rkind
   theta     = 0._rkind
   dstarinc  = 0._rkind
   thetainc  = 0._rkind
   rhoe      = wstat_bl(j99,1)
   pe        = wstat_bl(j99,5)
   ue        = utfav_bl(j99)
   do j=1,j99
    rho  = wstat_bl(j,1)/rhoe
    rhop = wstat_bl(j+1,1)/rhoe
    uu   = utfav_bl(j)/ue
    uup  = utfav_bl(j+1)/ue
    if (j==1) then
     dy = delta_bl(1)
    else
     dy = delta_bl(j) - delta_bl(j-1)
    endif 
    dyh  = 0.5_rkind*dy
!   Trapezoidal rule
    dstar = dstar       + dyh*((1._rkind-rho*uu)+(1._rkind-rhop*uup))
    theta = theta       + dyh*((rho*uu*(1._rkind-uu))+(rhop*uup*(1._rkind-uup)))
    dstarinc = dstarinc + dyh*((1._rkind-uu)+(1._rkind-uup))
    thetainc = thetainc + dyh*((uu*(1._rkind-uu))+(uup*(1._rkind-uup)))
   enddo
   shapef    = dstar/theta ! Shape factor H
   shapefinc = dstarinc/thetainc ! Incompressible Shape factor H_i
!
!  VAN DRIEST II
!
   tte  = pe/rhoe
   rmue = wstat_bl(j99,20) 
   ff   = ttw/tte
   aa   = ((rfac*0.2_rkind*Mach**2)/ff)**0.5_rkind
   bb   = (1+rfac*0.2_rkind*Mach**2-ff)/ff
   alpha = 2*aa**2-bb
   alpha = alpha/(4*aa**2+bb**2)**0.5_rkind
   beta = bb/(4*aa**2+bb**2)**0.5_rkind
   fc   = rfac*0.2_rkind*Mach**2
   fc   = fc/(asin(alpha)+asin(beta))**2
   ftheta = rmue/rmuw
   cfinc = cf*fc

   redelta99  = ue*rhoe*delta99/rmue
   retheta    = ue*rhoe*theta/rmue
   rethetainc = ue*rhoe*thetainc/rmue
   redelta2   = ue*rhoe*theta/rmuw ! = retheta*ftheta

   !              1        2        3      4              5       6       7       8     9       10
   write(10,100) xg(i,1),yg(i,1),delta99,delta_bl(jvort),dstar,theta,dstarinc,thetainc,shapef,shapefinc, &
   !                   11   12 13    14   15    16    17  18  19     20      21   22
                      rhow,ttw,ppw,prms,deltav,utau,tauw,cf,cfinc,cpcoeff,ue/u0,ufav_bl(jvort)/u0, &
   !                     23        24      25      26  27 28
                      redelta99,retheta,redelta2,retau,Bq,ch

  enddo
  close(10)
!
  open(10,file='POSTPRO/bl_suction.dat',form='formatted')
  do i=ile+2,itu-1,ix_out
   print*,'Analyzing wall i-index: ',i
   ! get wstat_bl(j,iv)
   call extract_bl(i)

   do j=1,npoints_bl
       ufav_bl(j)  = wstat_bl(j,13)/wstat_bl(j,1)
       vfav_bl(j)  = wstat_bl(j,14)/wstat_bl(j,1)
       wfav_bl(j)  = wstat_bl(j,15)/wstat_bl(j,1)
       tfav_bl(j)  = wstat_bl(j,25)/wstat_bl(j,1)
       utfav_bl(j) = (ufav_bl(j)*dxdcsi(i,1)+vfav_bl(j)*dydcsi(i,1))/csimod(i,1) ! csi tangent wall
       vtfav_bl(j) = (ufav_bl(j)*detadx(i,1)+vfav_bl(j)*detady(i,1))/meta(i,1) ! csi normal wall
   enddo

   dudyw  = -22._rkind*utfav_bl(1)+36._rkind*utfav_bl(2)-18._rkind*utfav_bl(3)+4._rkind*utfav_bl(4)
   dyw    = -22._rkind*0.000_rkind+36._rkind*delta_bl(1)-18._rkind*delta_bl(2)+4._rkind*delta_bl(3)
   dtdyw  = -22._rkind* tfav_bl(1)+36._rkind* tfav_bl(2)-18._rkind* tfav_bl(3)+4._rkind* tfav_bl(4)
   dudyw  = dudyw / dyw
   dtdyw  = dtdyw / dyw

   rhow   = wstat_bl(1,1)
   ttw    = tfav_bl(1)
   rmuw   = wstat_bl(1,20)
   rnuw   = rmuw/rhow

   ttw    = Twall
   rlamw  = rmuw*gam/(gam-1._rkind)/Prandtl
   tauw   = rmuw*dudyw
   qw     = -rlamw*dtdyw
   utau   = sqrt(abs(tauw)/rhow)
   deltav = rnuw/utau
   pdyn   = 0.5_rkind*rho0*u0**2
   cf     = tauw/pdyn
   trec   = t0 * (1.+0.5*(gam-1.)*rfac*Mach**2)
   Bq     = qw/((gam/(gam-1.))*rhow*utau*ttw)
   if (theta_wall==1._rkind) then
    ch = 0._rkind
   else
    ch = qw/u0/(gam/(gam-1.))/(ttw-trec)
   endif

   prms = sqrt(abs(wstat_bl(1,11)-wstat_bl(1,5)**2))/pdyn
   ppw   = rhow*ttw 

   cpcoeff = (ppw-p0)/pdyn

   !-----------------------------------------------------------
   ! VORTICAL CRTIERIUM
   !-----------------------------------------------------------
   !TEORICO omag_factor = 0.05_rkind
   omag_factor = 0.25_rkind
   !omag_factor = 0.95_rkind
   jvort = 1
   do j=1,npoints_bl-1
    vortx2 = wstat_bl(j,22)
    vorty2 = wstat_bl(j,23)
    vortz2 = wstat_bl(j,24)
    omag = sqrt(vortx2+vorty2+vortz2)
    !print*,'omag: ',omag,vortx2,vorty2,vortz2,omag_factor*u0/l0
    if (omag < omag_factor*u0/l0) then
     jvort = j
     exit
    endif
   enddo

   udel = 0.99_rkind*utfav_bl(jvort)

   j99 = 1
   do j=1,npoints_bl-1
    if (utfav_bl(j)>udel) then
     j99 = j-1
     exit
    endif
   enddo

   print*,'u0, l0: ',u0, l0
   print*,'jvort, j99: ',jvort, j99
   !-----------------------------------------------------------
   !-----------------------------------------------------------
   ! INVISCID CRTIERIUM
   !-----------------------------------------------------------
   ! Computation of delta_99 (see Griffin et al, PRF 2021)
   !uinv=0._rkind
   !gg = 2*gam/(gam-1)
   !do j=1,npoints_bl-1
   ! ry=wstat_bl(j,1)
   ! py=wstat_bl(j,5)
   ! uinv(j)=sqrt(gg*(p0/rho0-py/ry)+uref**2+vref**2-vtfav_bl(j)**2) ! inviscid velocity profile
   !enddo
   !j99   = 1
   !do j=1,npoints_bl-1
   ! uu = utfav_bl(j)
   ! udel = 0.99_rkind*uinv(j)
   ! if (uu>udel) then
   !  j99 = j-1
   !  exit
   ! endif
   !enddo
   !jvort = j99
   !print*,'j: ',jvort, j99
   !-----------------------------------------------------------
   jvortical(i) = jvort

   deltaj99 = delta_bl(j99-1)
   dely = delta_bl(j99) - delta_bl(j99-1)
   unum = udel-utfav_bl(j99)
   uden = utfav_bl(j99+1)-utfav_bl(j99)
   if(abs(uden) > 1e-6) then
       delta99 = deltaj99 + dely*(unum/uden) ! b._rkindl._rkind thickness
   else
       delta99 = deltaj99
   endif
   retau = delta99/deltav
   d99_vec(i)    = delta99
   deltav_vec(i) = deltav
   utau_vec(i)   = utau
   rhow_vec(i)   = rhow
   ttw_vec(i)    = ttw
   rmuw_vec(i)   = rmuw
!
!  Integral boundary layer thicknesses
!
   dstar     = 0._rkind
   theta     = 0._rkind
   dstarinc  = 0._rkind
   thetainc  = 0._rkind
   rhoe      = wstat_bl(j99,1)
   pe        = wstat_bl(j99,5)
   ue        = utfav_bl(j99)
   do j=1,j99
    rho  = wstat_bl(j,1)/rhoe
    rhop = wstat_bl(j+1,1)/rhoe
    uu   = utfav_bl(j)/ue
    uup  = utfav_bl(j+1)/ue
    if (j==1) then
     dy = delta_bl(1)
    else
     dy = delta_bl(j) - delta_bl(j-1)
    endif
    dyh  = 0.5_rkind*dy
!   Trapezoidal rule
    dstar = dstar       + dyh*((1._rkind-rho*uu)+(1._rkind-rhop*uup))
    theta = theta       + dyh*((rho*uu*(1._rkind-uu))+(rhop*uup*(1._rkind-uup)))
    dstarinc = dstarinc + dyh*((1._rkind-uu)+(1._rkind-uup))
    thetainc = thetainc + dyh*((uu*(1._rkind-uu))+(uup*(1._rkind-uup)))
   enddo
   shapef    = dstar/theta ! Shape factor H
   shapefinc = dstarinc/thetainc ! Incompressible Shape factor H_i
!
!  VAN DRIEST II
!
   tte  = pe/rhoe
   rmue = wstat_bl(j99,20) 
   ff   = ttw/tte
   aa   = ((rfac*0.2_rkind*Mach**2)/ff)**0.5_rkind
   bb   = (1+rfac*0.2_rkind*Mach**2-ff)/ff
   alpha = 2*aa**2-bb
   alpha = alpha/(4*aa**2+bb**2)**0.5_rkind
   beta = bb/(4*aa**2+bb**2)**0.5_rkind
   fc   = rfac*0.2_rkind*Mach**2
   fc   = fc/(asin(alpha)+asin(beta))**2
   ftheta = rmue/rmuw
   cfinc = cf*fc

   redelta99  = ue*rhoe*delta99/rmue
   retheta    = ue*rhoe*theta/rmue
   rethetainc = ue*rhoe*thetainc/rmue
   redelta2   = ue*rhoe*theta/rmuw ! = retheta*ftheta

   !              1        2        3      4              5       6       7       8     9       10
   write(10,100) xg(i,1),yg(i,1),delta99,delta_bl(jvort),dstar,theta,dstarinc,thetainc,shapef,shapefinc, &
   !                   11   12 13    14   15    16    17  18  19     20      21   22
                      rhow,ttw,ppw,prms,deltav,utau,tauw,cf,cfinc,cpcoeff,ue/u0,ufav_bl(jvort)/u0, &
   !                     23        24      25      26  27 28
                      redelta99,retheta,redelta2,retau,Bq,ch, real(jvort), real(j99)

  enddo
  close(10)

  print*, 'Start/end indices of suction side =',ile,itu
!
  do ii = 1,nstatloc
   i = ixstat(ii)
   write(chstat,1006) i
   write(*,*)'analyzing profile,',i
   call extract_bl(i)

   do j=1,npoints_bl
       ufav_bl(j)    = wstat_bl(j,13)/wstat_bl(j,1)
       vfav_bl(j)    = wstat_bl(j,14)/wstat_bl(j,1)
       wfav_bl(j)    = wstat_bl(j,15)/wstat_bl(j,1)
       tfav_bl(j)    = wstat_bl(j,25)/wstat_bl(j,1)
       utfav_bl(j)   = (ufav_bl(j)*dxdcsi(i,1)+vfav_bl(j)*dydcsi(i,1))/csimod(i,1) ! csi tangent wall
       vtfav_bl(j)   = (ufav_bl(j)*detadx(i,1)+vfav_bl(j)*detady(i,1))/meta(i,1)   ! csi normal
       rhout2_bl(j)  = (wstat_bl(j,16)*dxdcsi(i,1)**2+wstat_bl(j,17)*dydcsi(i,1)**2+ &
                       2._rkind*wstat_bl(j,19)*dxdcsi(i,1)*dydcsi(i,1))/csimod(i,1)**2
       rhovt2_bl(j)  = (wstat_bl(j,16)*detadx(i,1)**2+wstat_bl(j,17)*detady(i,1)**2+ &
                       2._rkind*wstat_bl(j,19)*detadx(i,1)*detady(i,1))/meta(i,1)**2
       rhoutvt_bl(j) = (wstat_bl(j,16)*dxdcsi(i,1)*detadx(i,1)+wstat_bl(j,17)*dydcsi(i,1)*detady(i,1)+ &
                       wstat_bl(j,19)*(dxdcsi(i,1)*detady(i,1)+dydcsi(i,1)*detadx(i,1)))/ &
                       csimod(i,1)/meta(i,1)
   enddo
!
   jvort = jvortical(i)

   rhow    = rhow_vec(i)
   rmuw    = rmuw_vec(i)
   rnuw    = rmuw/rhow
   delta99 = d99_vec(i)
   deltav  = deltav_vec(i)
   utau    = utau_vec(i)
   tauw    = rhow*utau**2
   ttw     = ttw_vec(i)
!
! Compute compressibility transformations
!
   vaux(1,:) = wstat_bl(:,1)/rhow
   vaux(2,:) = wstat_bl(:,21)/rnuw
   vaux(3,:) = wstat_bl(:,20)/rmuw
   vaux(4,:) = delta_bl(1:npoints_bl)/deltav ! yplus
   vaux(5,:) = delta_bl(1:npoints_bl)*sqrt(tauw*wstat_bl(:,1))/wstat_bl(:,20) ! ystar
   vaux(6,:) = utau/sqrt(gam*ttw) ! Mach_tau constant
!
   tname = 'TrettelLarsson'
   call transform(0,npoints_bl,delta_bl,utfav_bl,naux,vaux,yt,ut,tname)
!
   tname = 'vanDriest'
   call transform(0,npoints_bl,delta_bl,utfav_bl,naux,vaux,yvd,uvd,tname)
!
   tname = 'Volpiani'
   call transform(0,npoints_bl,delta_bl,utfav_bl,naux,vaux,yv,uv,tname)

   tname = 'Hasan'
   call transform(0,npoints_bl,delta_bl,utfav_bl,naux,vaux,yh,uh,tname)
!
   open(unit=15,file='POSTPRO/stat_'//chstat//'.prof')
   do j=1,npoints_bl
    if(j == 1) then
      d_bl = 0._rkind
    else
      d_bl = delta_bl(j-1)
    endif
    y99    = d_bl/delta99
    yp     = d_bl/deltav
    ystar  = vaux(5,j)
    up     = utfav_bl(j)/utau
    uvdp   = uvd(j)/utau

    u2p  = abs(rhout2_bl(j)-wstat_bl(j,1)*utfav_bl(j)*utfav_bl(j))/tauw
    v2p  = abs(rhovt2_bl(j)-wstat_bl(j,1)*vtfav_bl(j)*vtfav_bl(j))/tauw
    w2p  = abs(wstat_bl(j,18)-wstat_bl(j,1)*wfav_bl(j)*wfav_bl(j))/tauw
    uvp  = abs(rhoutvt_bl(j)-wstat_bl(j,1)*utfav_bl(j)*vtfav_bl(j))/tauw

    rhofac = sqrt(wstat_bl(j,1)/rhow)
    ytp    = yt(j)/deltav
    yvp    = yv(j)/deltav
    utp    = ut(j)/utau
    uhp    = uh(j)/utau
    uvpl   = uv(j)/utau

    rhom    = wstat_bl(j,1)
    ttm     = tfav_bl(j)
    ppm     = ttm*rhom
    Mtau    = utau/sqrt(gam*ttw)
    rhormsp = (sqrt(wstat_bl(j, 7)-wstat_bl(j,1)**2))/(rhow*gam*Mtau**2)
    trmsp   = (sqrt(wstat_bl(j,26)/wstat_bl(j,1)-tfav_bl(j)**2))/(ttw*gam*Mtau**2)
    prms    = sqrt(abs(wstat_bl(j,11)-wstat_bl(j,5)**2))
    mum     = wstat_bl(j,20)

    !             1     2    3            4          5           6         7   8   9
    write(15,100) d_bl,rhom,ufav_bl(j),vfav_bl(j),utfav_bl(j),vtfav_bl(j),ttm,ppm,mum, &
    !             10  11 12    13 
                  y99,yp,ystar,yvp, &
    !             14  15   16  17  18
                  up,uvdp,utp,uvpl,uhp, &
    !              19  20  21  22
                  u2p,v2p,w2p,uvp, &
    !              23      24     25
                  rhormsp,trmsp,prms/p0
   enddo
   close(15)
  enddo

 100  format(200ES20.10)
 1006 format(I6.6)

 end subroutine stats2d_airfoil

end module postpro_airfoil
