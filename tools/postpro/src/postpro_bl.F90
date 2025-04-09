module postpro_bl
 use parameters
 use global_variables
 use reader
 use derivatives
 use comp_transform
!
 contains
!
 subroutine stats2d
 implicit none
 integer, parameter :: naux = 6
 integer :: i,j,m,j99,ii,jint
 real(rkind), dimension(nx,ny) :: ufav,vfav,wfav,tfav,dufavdy
 real(rkind), dimension(nx) :: d99_vec,deltav_vec,utau_vec,rhow_vec,ttw_vec,dudyw_vec,rmuw_vec
 real(rkind), dimension(naux,ny) :: vaux
 real(rkind), dimension(ny) :: yt,yh,yvd,yv,uv,ut,uvd,uh,duvddy,dutdy
 real(rkind), dimension(ny) :: duvdy,duhdy,dugdy,ybig,ut_2,ug,dybigdy
 real(rkind) :: duw,dtw,dyw,dudyw,rmuw,rnuw,tauw,rhow,ttw,deltav,pdyn
 real(rkind) :: cf,u99,dely,utau,uu,unum,uden,delta99,retau,prmsp,rhormsp,trmsp
 real(rkind) :: dstarinc,thetainc,rhoe,pe,ue,tte,trec
 real(rkind) :: rho,rhop,up,dy,dyh,dstar,theta,uup
 real(rkind) :: shapef,shapefinc,rmue,ff,aa,bb,aai,bbi,cci,alpha,beta,aad,bbd,ccd
 real(rkind) :: fc,ftheta,cfinc,rethetainc,redelta2,retheta
 real(rkind) :: y99,uvdp,u2p,v2p,w2p,uvp,rhofac,utp,uvpl,yvp,ytp,uhp,rhom,ystar,ugp
 real(rkind) :: yp,ttm,ppw,Mtau,mum,ppm,prms
 real(rkind) :: qw,ch,dtdyw,rlamw,det,Bq,redelta99
 real(rkind) :: duthr,mach_edge,rmfac
 real(rkind) :: intj,intjm,integrand
 real(rkind) :: seq,seqp,setl,setlp
 character(20) :: tname
 character(6)  :: chstat
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
    tfav(i,j) = wstat(i,j,25)/wstat(i,j,1)
!   write(10,100) x(i),y(j),(wstat(i,j,m),m=1,nv)
   enddo
  enddo
! close(10)
  do i=1,nx
   call ddy(ny,ufav(i,:),y,2,dufavdy(i,:))
  enddo
!
! Mean boundary layer properties
!
  u99 = 0.99_rkind*u0
!
  open(10,file='POSTPRO/cf.dat',form='formatted')
  do i=1,nx
   if (ystag>0) then ! staggered grid
!   Coefficients from Lele pag 37 
!   aad = 450._rkind/384._rkind
!   bbd = -75._rkind/384._rkind
!   ccd =   9._rkind/384._rkind
    aad =  9._rkind/8._rkind
    bbd = -1._rkind/8._rkind
    ccd =  0._rkind/8._rkind
    duw = 2._rkind*aad*ufav(i,1)+2._rkind/3._rkind*bbd*ufav(i,2)+2._rkind/5._rkind*ccd*ufav(i,3)
    dyw = 2._rkind*aad*y(1)     +2._rkind/3._rkind*bbd*y(2)     +2._rkind/5._rkind*ccd*y(3)
    dtw = 2._rkind*aad*tfav(i,1)+2._rkind/3._rkind*bbd*tfav(i,2)+2._rkind/5._rkind*ccd*tfav(i,3)
    dtw = dtw-2._rkind*Twall*(aad+bbd/3._rkind+ccd/5._rkind)
    dudyw  = duw/dyw
    dtdyw  = dtw/dyw
!    Coefficients from Lele pag 39
     aai = 150._rkind/128._rkind
     bbi = -25._rkind/128._rkind
     cci =   3._rkind/128._rkind
     rhow  = aai*wstat(i,1,1 )+bbi*wstat(i,2,1 )+cci*wstat(i,3,1 )
     ttw   = aai*wstat(i,1,6 )+bbi*wstat(i,2,6 )+cci*wstat(i,3,6 )
     rmuw  = aai*wstat(i,1,20)+bbi*wstat(i,2,20)+cci*wstat(i,3,20)
     rnuw  = rmuw/rhow
   else
    duw = (-22.*ufav(i,1)+36.*ufav(i,2)-18.*ufav(i,3)+4.*ufav(i,4))!/12
    dtw = (-22.*tfav(i,1)+36.*tfav(i,2)-18.*tfav(i,3)+4.*tfav(i,4))!/12.
    dyw = (-22.*y(1)+36.*y(2)-18.*y(3)+ 4.*y(4))!/12.
    dudyw = duw/dyw
    dtdyw = dtw/dyw
    rhow  = wstat(i,1,1)
    ttw   = tfav(i,1)
    rmuw  = wstat(i,1,20)
    rnuw  = rmuw/rhow
   endif
   ttw    = Twall ! overwrite Twall with exact boundary condition
   rlamw  = rmuw*gam/(gam-1._rkind)/Prandtl
   tauw   = rmuw*dudyw
   qw     = -rlamw*dtdyw
   utau   = sqrt(abs(tauw)/rhow)
   deltav = rnuw/utau
   pdyn  = 0.5_rkind*rho0*u0**2
   cf    = tauw/pdyn
   trec  = t0 * (1.+0.5*(gam-1.)*rfac*Mach**2)
   Bq    = qw/((gam/(gam-1.))*rhow*utau*ttw)
   if (theta_wall==1._rkind) then
    ch = 0._rkind
   else
    ch = qw/u0/(gam/(gam-1.))/(ttw-trec)
   endif
   j99   = 1
   do j=1,ny-1
    uu = ufav(i,j)
    if (uu>u99) then
     j99 = j-1
     exit
    endif
   enddo
   dely = y(j99+1)-y(j99)
   unum = u99-ufav(i,j99)
   uden = ufav(i,j99+1)-ufav(i,j99)
   delta99 = y(j99)+dely*(unum/uden) ! bl thickness
   retau = delta99/deltav
   prmsp = sqrt(abs(wstat(i,1,11)-wstat(i,1,5)**2))/tauw
   ppw           = rhow*ttw 
   d99_vec(i)    = delta99
   deltav_vec(i) = deltav
   utau_vec(i)   = utau
   rhow_vec(i)   = rhow
   ttw_vec(i)    = ttw
   rmuw_vec(i)   = rmuw
   dudyw_vec(i)  = dudyw 
!
!  Integral boundary layer thicknesses
!
!  jint = ny
!  if (ystag==0) jint = ny-1
   duthr = 0.01_rkind*dufavdy(i,j99)
   jint  = 1
   do j=j99+1,ny-1
    uu = ufav(i,j)
    if (dufavdy(i,j)<duthr) then
     jint = j
     exit
    endif
   enddo
   dstar      = 0._rkind ! displacement thickness
   theta      = 0._rkind ! momentum thickness
   dstarinc   = 0._rkind ! inc displacement
   thetainc   = 0._rkind ! inc momentum
   rhoe       = wstat(i,jint,1)
   pe         = wstat(i,jint,5)
   ue         = ufav(i,jint)
   if (ystag>0) then
    do j=1,jint
     rho  = wstat(i,j,1)/rhoe
     uu   = ufav(i,j)/ue
     dy   = yn(j+1)-yn(j)
!    Trapezoidal rule
     dstar = dstar       + dy*(1._rkind-rho*uu)
     theta = theta       + dy*(rho*uu*(1._rkind-uu))
     dstarinc = dstarinc + dy*(1._rkind-uu)
     thetainc = thetainc + dy*(uu*(1._rkind-uu))
    enddo
   else
    do j=1,jint
     rho  = wstat(i,j, 1)/rhoe
     rhop = wstat(i,j+1,1)/rhoe
     uu   = ufav(i,j)/ue
     uup  = ufav(i,j+1)/ue
     dy   = y(j+1)-y(j)
     dyh  = 0.5_rkind*dy
!    Trapezoidal rule
     dstar = dstar       + dyh*((1._rkind-rho*uu)+(1._rkind-rhop*uup))
     theta = theta       + dyh*((rho*uu*(1._rkind-uu))+(rhop*uup*(1._rkind-uup)))
     dstarinc = dstarinc + dyh*((1._rkind-uu)+(1._rkind-uup))
     thetainc = thetainc + dyh*((uu*(1._rkind-uu))+(uup*(1._rkind-uup)))
    enddo
   endif
   shapef    = dstar/theta ! Shape factor H
   shapefinc = dstarinc/thetainc ! Incompressible Shape factor H_i
!
!  VAN DRIEST II pag 996 Hopkins and Inouye AIAA J 1971
!
   tte    = pe/rhoe
   rmue   = wstat(i,jint,20) 
   mach_edge = ue/sqrt(gam*tte)
   !print *, mach_edge, wstat(33,i,jint)
   rmfac  = rfac*0.2_rkind*mach_edge**2
   ff     = ttw/tte
   aa     = sqrt(rmfac/ff)
   bb     = (1._rkind+rmfac-ff)/ff
   alpha  = 2*aa**2-bb
   alpha  = alpha/sqrt(4*aa**2+bb**2)
   beta   = bb   /sqrt(4*aa**2+bb**2)
   fc     = rmfac/(asin(alpha)+asin(beta))**2
   ftheta = rmue/rmuw
   cfinc  = cf*fc
   redelta99  = ue*rhoe*delta99/rmue
   retheta    = ue*rhoe*theta/rmue
   rethetainc = ue*rhoe*thetainc/rmue
   redelta2   = ue*rhoe*theta/rmuw ! = retheta*ftheta
!
   !write(10,100) x(i),cf,retau,shapef,shapefinc,delta99,dstar,theta,utau/u0,rethetainc,&
   !              cfinc,redelta2,prmsp,ch,retheta
   !              1      2      3     4       5        6       7       8
   write(10,100) x(i),delta99,dstar,theta,dstarinc,thetainc,shapef,shapefinc, &
   !                    9   10  11   12    13    14   15  
                      rhow,ttw,ppw,prmsp,utau/u0,cf,cfinc, &
   !                      16      17       18      19  20 21
                      redelta99,retheta,redelta2,retau,Bq,ch
  enddo
  close(10)
!
  do ii = 1,nstatloc
   i = ixstat(ii)
   write(chstat,1006) i
   write(*,*)'writing stat,',i,chstat
!
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
   wstat(i,:,21) = wstat(i,:,20)/wstat(i,:,1)
!   
   vaux(1,:) = wstat(i,:,1 )/rhow
   vaux(2,:) = wstat(i,:,21)/rnuw
!  vaux(3,:) = wstat(i,:,20)/rmuw
   vaux(3,:) = vaux(1,:)*vaux(2,:)
   vaux(4,:) = y(1:ny)/deltav ! yplus
   vaux(5,:) = y(1:ny)*sqrt(tauw*wstat(i,:,1))/wstat(i,:,20) ! ystar
   vaux(6,:) = utau/sqrt(gam*ttw) ! Mach_tau constant
!   
   do j=1,ny
    ybig(j) = y(j)*sqrt(vaux(1,j))/vaux(3,j) ! when normalized by deltav corresponds to ystar
   enddo
!
!  TL (direct integration)
!  ut_2(1) = vaux(3,1)*(ufav(i,1))/y(1)*ybig(1)
!  do j=2,ny
!   ut_2(j) = ut_2(j-1)+vaux(3,j)*(ufav(i,j)-ufav(i,j-1))/(y(j)-y(j-1))*(ybig(j)-ybig(j-1))
!  enddo
!
   tname = 'vanDriest'
   call transform(ystag,ny,y(1:ny),ufav(i,1:ny),naux,vaux,yvd,uvd,tname)
   call ddy(ny,uvd(:),yvd,2,duvddy(:)) ! yvd = yp

   tname = 'TrettelLarsson'
   call transform(ystag,ny,y(1:ny),ufav(i,1:ny),naux,vaux,yt,ut,tname)
   call ddy(ny,ut(:),yt,2,dutdy(:))

   tname = 'Volpiani'
   call transform(ystag,ny,y(1:ny),ufav(i,1:ny),naux,vaux,yv,uv,tname)
   call ddy(ny,uv(:),yv,2,duvdy(:))

   tname = 'Hasan'
   call transform(ystag,ny,y(1:ny),ufav(i,1:ny),naux,vaux,yh,uh,tname)
   call ddy(ny,uh(:),yh,2,duhdy(:))
!
   tname = 'Griffin'
   seq   = ufav(i,1)/ybig(1)
   seqp  = seq/utau*deltav
   setl  = ufav(i,1)/y(1)
   setlp = setl/utau*deltav
   integrand  = seq/(1._rkind+seqp-setlp)
   ug(1) = integrand*ybig(1)
   do j=2,ny
    seq   = 1._rkind/vaux(3,j)*(ufav(i,j)-ufav(i,j-1))/(ybig(j)-ybig(j-1))
    seqp  = seq/utau*deltav
    setl  = vaux(3,j)*(ufav(i,j)-ufav(i,j-1))/(y(j)-y(j-1))
    setlp = setl/utau*deltav
    integrand  = seq/(1._rkind+seqp-setlp)
    ug(j) = ug(j-1)+integrand*(ybig(j)-ybig(j-1))
   enddo
   call ddy(ny,ug(:),yh,2,dugdy(:))
!
   open(unit=15,file='POSTPRO/stat_'//chstat//'.prof')
   do j=1,ny
    y99       = y(j)/delta99
    yp        = y(j)/deltav
    ystar     = vaux(5,j)
    up        = ufav(i,j)/utau
    uvdp      = uvd(j)/utau
    duvddy(j) = yvd(j)*duvddy(j)/utau
    dutdy(j)  = yt(j)*dutdy(j)/utau
    duvdy(j)  = yv(j)*duvdy(j)/utau
    duhdy(j)  = yh(j)*duhdy(j)/utau
    dugdy(j)  = yh(j)*dugdy(j)/utau
    u2p     = (wstat(i,j,16)-wstat(i,j,1)*ufav(i,j)*ufav(i,j))/tauw
    v2p     = (wstat(i,j,17)-wstat(i,j,1)*vfav(i,j)*vfav(i,j))/tauw
    w2p     = (wstat(i,j,18)-wstat(i,j,1)*wfav(i,j)*wfav(i,j))/tauw
    uvp     = (wstat(i,j,19)-wstat(i,j,1)*ufav(i,j)*vfav(i,j))/tauw
    rhofac  = sqrt(wstat(i,j,1)/rhow)
    ytp     = yt(j)/deltav
    yvp     = yv(j)/deltav
    utp     = ut(j)/utau
    uhp     = uh(j)/utau
    ugp     = ug(j)/utau
    uvpl    = uv(j)/utau
    rhom    = wstat(i,j,1)
    ttm     = tfav(i,j)
    ppm     = ttm*rhom
    Mtau    = utau/sqrt(gam*ttw)
    rhormsp = (sqrt(wstat(i,j, 7)-wstat(i,j,1)**2))/(rhow*gam*Mtau**2)
    trmsp   = (sqrt(wstat(i,j,26)/wstat(i,j,1)-tfav(i,j)**2))/(ttw*gam*Mtau**2)
    prmsp    = sqrt(abs(wstat(i,j,11)-wstat(i,j,5)**2))/tauw
    mum     = wstat(i,j,20)/wstat(i,1,20)
    ! note: ytp=ystar, yh=ystar, yvd=yp
    !              1     2    3         4        5   6   7
    write(15,100) y(j),rhom,ufav(i,j)/u0,vfav(i,j)/u0,ttm,ppm,mum, &
    !              8   9   10  11
                  y99,yp,ystar,yvp, &
    !             12  13   14  15   16  17
                  up,uvdp,utp,uvpl,ugp,uhp, &
    !             18          19       20       21      22
                  duvddy(j),dutdy(j),duvdy(j),dugdy(j),duhdy(j), &
    !             23  24  25  26
                  u2p,v2p,w2p,uvp, &
    !             27       28    29
                  rhormsp,trmsp,prmsp
   enddo
   close(15)
!
  enddo
 100  format(200ES20.10)
 1006 format(I6.6)
 end subroutine stats2d
!
end module
