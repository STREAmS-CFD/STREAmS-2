module postpro_bl
 use parameters
 use global_variables
 use reader
 use derivatives
 use comp_transform

 contains
!
 subroutine stats2d
 implicit none
 integer, parameter :: np = 7
 integer, parameter :: naux = 3
 integer, parameter :: nv_print = 34 !print less statistics than total (see compute_stats in singleideal.F90 for indices)
 integer :: i,j,m,j99,ii
 real(rkind), dimension(nx,ny) :: ufav, vfav,wfav
 real(rkind), dimension(np)    :: fvec, yvec
 real(rkind), dimension(ny)    :: uvd
 real(rkind), dimension(nx) :: d99_vec,deltav_vec,utau_vec
 real(rkind), dimension(naux,ny) :: vaux
 real(rkind), dimension(ny) :: yt,yvd,yv,uv,ut
 real(rkind) :: uwall, dudyw, d2udyw, rmuw, rnuw,tauw,rhow,ttw,deltav,pdyn
 real(rkind) :: cf,udel,dely,utau,uu,unum,uden,delta99,retau,prmsp 
 real(rkind) :: dstarinc,thetainc,rhoe,pe,ue,tte 
 real(rkind) :: rho,rhop,up,dy,dyh,dstar,theta,uup
 real(rkind) :: shapef,shapefinc,rmue,ff,aa,bb,alpha,beta 
 real(rkind) :: fc,ftheta,cfinc,rethetainc,rethetawall,reout 
 real(rkind) :: y99,uvdp,u2p,v2p,w2p,uvp,rhofac,utp,uvpl,yvp,ytp
 real(rkind) :: yp 
 character(20) :: tname
 character(6)  :: chstat
!
  reout = u0/mu0
!
  call read_grid_bl
!
! open(unit=10,file='wavplot.dat',form='formatted')
! write(10,*) 'zone i=',nx,', j=',ny
  do j=1,ny
   do i=1,nx
    ufav(i,j) = wstat(i,j,13)/wstat(i,j,1)
    vfav(i,j) = wstat(i,j,14)/wstat(i,j,1)
    wfav(i,j) = wstat(i,j,15)/wstat(i,j,1)
!   write(10,100) x(i),y(j),(wstat(i,j,m),m=1,nv_print)
   enddo
  enddo
! close(10)
!
! Mean boundary layer properties
!
  udel = 0.99_rkind*u0
!
  open(10,file='POSTPRO/cf.dat',form='formatted')
  do i=1,nx
   fvec  = ufav(i,1:np)
   yvec  = y(1:np)
   call interpolate(np,fvec,yvec,0._rkind,uwall,dudyw,d2udyw)
   rhow  = wstat(i,1,1)
   ttw   = wstat(i,1,6)
   rmuw  = wstat(i,1,20)
   tauw  = rmuw*dudyw
   utau  = sqrt(abs(tauw)/rhow)
   rnuw  = rmuw/rhow
   deltav= rnuw/utau
   pdyn  = 0.5_rkind*u0**2
   cf    = tauw/pdyn
   j99   = 1
   do j=1,ny-1
    uu = ufav(i,j)
    if (uu>udel) then
     j99 = j-1
     exit
    endif
   enddo
   dely = y(j99+1)-y(j99)
   unum = udel-ufav(i,j99)
   uden = ufav(i,j99+1)-ufav(i,j99)
   delta99 = y(j99)+dely*(unum/uden) ! b._rkindl._rkind thickness
   retau = delta99/deltav
   prmsp = sqrt(abs(wstat(i,1,11)-wstat(i,1,5)**2))/tauw
   d99_vec(i)    = delta99
   deltav_vec(i) = deltav
   utau_vec(i)   = utau
!
!  Integral boundary layer thicknesses
!
   dstar     = 0._rkind
   theta     = 0._rkind
   dstarinc  = 0._rkind
   thetainc  = 0._rkind
   rhoe      = wstat(i,j99,1)
   pe        = wstat(i,j99,5)
   ue        = ufav(i,j99)
   do j=1,j99
    rho  = wstat(i,j,1)/rhoe
    rhop = wstat(i,j+1,1)/rhoe
    uu   = ufav(i,j)/ue
    uup  = ufav(i,j+1)/ue
    dy   = y(j+1)-y(j)
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
   rmue = wstat(i,j99,20) 
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
   rethetainc = reout*thetainc*ftheta
   rethetawall = reout*theta*rmue/rmuw
!
   write(10,100) x(i),cf,retau,shapef,shapefinc,delta99,dstar,theta,utau/u0,rethetainc,cfinc,rethetawall,prmsp
  enddo
  close(10)
!
  do ii = 1,nstatloc
   i = ixstat(ii)
   write(chstat,1006) i
   write(*,*)'writing stat,',i,chstat
!
   rhow    = wstat(i,1,1)
   delta99 = d99_vec(i)
   deltav  = deltav_vec(i)
   utau    = utau_vec(i)
!
!
! Compute compressibility transformations
!
   vaux(1,:) = wstat(i,:,1)/rhow
   vaux(2,:) = wstat(i,:,21)/rnuw
   vaux(3,:) = wstat(i,:,20)/rmuw
!
   tname = 'TrettelLarsson'
   call transform(ny,y(1:ny),ufav(i,1:ny),naux,vaux,yt,ut,tname)
!
   tname = 'vanDriest'
   call transform(ny,y(1:ny),ufav(i,1:ny),naux,vaux,yvd,uvd,tname)
!
   tname = 'Volpiani'
   call transform(ny,y(1:ny),ufav(i,1:ny),naux,vaux,yv,uv,tname)
!
   open(unit=15,file='POSTPRO/stat_'//chstat//'.prof')
   do j=1,ny
    y99    = y(j)/delta99
    yp     = y(j)/deltav
    up     = ufav(i,j)/utau
    uvdp   = uvd(j)/utau
    u2p  = abs(wstat(i,j,16)-wstat(i,j,1)*ufav(i,j)*ufav(i,j))/tauw
    v2p  = abs(wstat(i,j,17)-wstat(i,j,1)*vfav(i,j)*vfav(i,j))/tauw
    w2p  = abs(wstat(i,j,18)-wstat(i,j,1)*wfav(i,j)*wfav(i,j))/tauw
    uvp    =        (wstat(i,j,19)/wstat(i,j,1)-ufav(i,j)*vfav(i,j)) /utau**2
    rhofac = sqrt(wstat(i,j,1)/rhow)
    prmsp  = sqrt(abs(wstat(i,j,11)-wstat(i,j,5)**2))/tauw
    ytp    = yt(j)/deltav
    yvp    = yt(j)/deltav
    utp    = ut(j)/utau
    uvpl   = uv(j)/utau
    write(15,100) y99,yp,up,uvdp,u2p,v2p,w2p,uvp,rhofac,prmsp,ytp,utp,yvp,uvpl
   enddo
   close(15)
!
  enddo
 100  format(200ES20.10)
 1006 format(I6.6)
 end subroutine stats2d
!
end module
