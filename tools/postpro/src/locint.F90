      !{{{ int_quasibil
      ! mi metto in un riferimento che è quasi quello csi,eta, tramite
      ! l'angolo theta. è in questo piano che faccio l'interpolazione
      ! bilineare
      subroutine int_quasibil(v,p,xst,pint)
!
!     FS: secondo me la prossima è sbagliata perché pint non deve essere solo reale ma double precision
      !implicit double precision(a-h,o,q-z) ! p is real by default
      implicit double precision(a-h,o-z) ! p is real by default
!
      dimension v(4,2) ! vertici del quadrato
      dimension p(4)   ! valore in quei vertici della variabile da interpolare
      dimension xst(2) ! è dove voglio la funzione

      idebug = 0
!
!     assumo che i lati a i cost siano quelli più diritti, e trovo theta
!     come angolo medio tra theta1 e theta2 formati dai singoli lati
!     "diritti"
!
!     trovo sin e cos degli angoli dei lati diritti
      xx = v(3,1)-v(2,1)
      yy = v(3,2)-v(2,2)
      dl32 = sqrt(xx*xx+yy*yy)
      costh1 = yy/dl32
      sinth1 = -xx/dl32
      xx = v(4,1)-v(1,1)
      yy = v(4,2)-v(1,2)
      dl41 = sqrt(xx*xx+yy*yy)
      costh2 = yy/dl41
      sinth2 = -xx/dl41
!
!     stimo l'angolo medio, assumendo sia vicino agli altri due
      dth = 0.5*(sinth2*costh1 - sinth1*costh2)
      costh = costh1 - sinth1*dth
      sinth = sinth1 + costh1*dth
!
!     mi metto nel riferimento csi,eta (v1 è l'origine)
      csi1 = 0.
      eta1 = 0.
      xx = (v(2,1)-v(1,1))
      yy = (v(2,2)-v(1,2))
      csi2 =  xx*costh + yy*sinth
      eta2 = -xx*sinth + yy*costh
      xx = (v(3,1)-v(1,1))
      yy = (v(3,2)-v(1,2))
      csi3 =  xx*costh + yy*sinth
      eta3 = -xx*sinth + yy*costh
      xx = (v(4,1)-v(1,1))
      yy = (v(4,2)-v(1,2))
      csi4 =  xx*costh + yy*sinth
      eta4 = -xx*sinth + yy*costh
      xx = (xst(1)-v(1,1))
      yy = (xst(2)-v(1,2))
      csi5 =  xx*costh + yy*sinth
      eta5 = -xx*sinth + yy*costh
!
!     DEBUG
      if (idebug.eq.1) then
       open(14,file='debug_csieta.dat')
       write(14,*) csi1, eta1
       write(14,*) csi2, eta2
       write(14,*) csi3, eta3
       write(14,*) csi4, eta4
       write(14,*) csi1, eta1
       close(14)
      endif
!     stop
!
!     approssimazione bilineare: medio v1 e v2, poi v3 e v4, e usando le
!     loro eta medie (mediate allo stesso modo), medio quello che esce
      p1 = p(1)
      p2 = p(2)
      p3 = p(3)
      p4 = p(4)
!
!     medie in csi
      csi21 = csi2-csi1
      csi21i = 1./csi21
      csi51 = csi5-csi1
      csi52 = csi5-csi2
      p12 = csi51*csi21i*p2 - csi52*csi21i*p1
      eta12 = csi51*csi21i*eta2 - csi52*csi21i*eta1
      if (idebug.eq.1) then
       print*, csi21, csi21i, csi51,csi52,p12,p1,p2,eta12,'primo' ! debug
       print*, csi1, csi2, csi5, p1,p2, p12, 'prima media' ! debug
      endif
!
      csi34 = csi3-csi4
      csi34i = 1./csi34
      csi54 = csi5-csi4
      csi53 = csi5-csi3
      p34 = csi54*csi34i*p3 - csi53*csi34i*p4
      eta34 = csi54*csi34i*eta3 - csi53*csi34i*eta4
      if (idebug.eq.1) then
       print*, csi34, csi34i, csi54,csi53,p34,p3,p4,eta34,'secondo' ! debug
       print*, csi3, csi4, csi5, p3,p4, p34, 'seconda media' ! debug
      endif
!
!     finalmente, calcolo il valore interpolato
      eta_ba = eta34 - eta12
      eta_bai = 1./eta_ba
      eta_5a = eta5 - eta12
      eta_5b = eta5 - eta34
      pint = eta_5a*eta_bai*p34 - eta_5b*eta_bai*p12
      if (idebug.eq.1) then
       print*, eta_ba, eta_bai, eta_5a, eta_5b, pint, 'last' ! debug
       print*, eta12, eta34, eta5, p12, p34, pint, 'licenza' ! debug
      endif
!
      return
      end
      !}}}
!     !{{{ distmin
      ! imin t.c. minima distanza tra il punto xst e la curva xvec
      subroutine distmin(xst,xvec,istart,iend,ng,i1,i2,imin)
!
      implicit double precision(a-h,o-z)
!
      dimension xvec(2,istart-ng:iend+ng)
      dimension xst(2)
!
      dist_min = 1.d9 ! inizializzo con grande distanza
      do i=i1,i2
       xx = xst(1) - xvec(1,i)
       yy = xst(2) - xvec(2,i)
       dist = sqrt(xx*xx+yy*yy)
       dist_min = min(dist,dist_min)
       if (dist_min.eq.dist) imin = i
      enddo
!
      return
      end
      !}}}
!     !{{{ alpha
      ! assume angolo piccolo, nel senso che le frontiere le cammino
      ! cella per cella (ergo mi basta asin)
      ! ritorna anche info sui moduli
      subroutine alpha(xst,x1,x2,a,imf)
!
      implicit double precision(a-h,o-z)
!
      dimension xst(2), x1(2), x2(2), xx1(2), xx2(2)
!
      xx1(:) = x1(:) - xst(:)
      xx2(:) = x2(:) - xst(:)
!
      vec = xx1(1)*xx2(2) - xx1(2)*xx2(1) ! third component of the vector product
      sca = xx1(1)*xx2(1) + xx1(2)*xx2(2) ! scalar product
!
      xx1m = xx1(1)*xx1(1) + xx1(2)*xx1(2) ! amplitude squared
      xx2m = xx2(1)*xx2(1) + xx2(2)*xx2(2) ! amplitude squared
      rmod = sqrt(xx1m*xx2m)
!
      imf = 1 ! 1 means none of the points coincide with xst
      if (rmod.eq.0.) then
       imf = 0 ! one of the points coincide with xst
      else
       pi = 4.*atan(1.)
       c = sca/rmod ! cos
       if (abs(c).gt.1.) c = c/abs(c)
       s = vec/rmod ! sin
       if (abs(s).gt.1.) s = s/abs(s)
       al = asin(abs(s)) ! 0. < al < pi/2
       ! -pi < a < pi (always, since angle of a segment seen by a point)
       if ((s.ge.0.).and.(c.ge.0.)) then
        a = al
       elseif ((s.ge.0.).and.(c.lt.0.)) then
        a = pi-al
       elseif ((s.lt.0.).and.(c.lt.0.)) then
        a = -pi+al
       elseif ((s.lt.0.).and.(c.ge.0.)) then
        a = -al
       endif
!
!      print*, 'a', a*180./pi
!
      endif
!
      return
      end
      !}}}
      !{{{ robust_search
      subroutine robust_search(i1,i2,j1,j2, &
                                        xcg,istart,iend,ny, xst, &
                                        ii,jj,ierr)
!
      implicit double precision(a-h,o-z)
!
      dimension xcg(2,istart:iend,1:ny)
      dimension xst(2)
!
      pi = 4.*atan(1.)
!
!     REFINED SEARCH (slow, but works always)
!
!     checking all the cells to find the one with alpha=360
      ierr = 1
      do j=j1,j2-1
       do i=i1,i2-1
!       print*, i,j
!
        call alpha(xst,xcg(:,i+1,j  ),xcg(:,i+1,j+1),a1,imf1)
        call alpha(xst,xcg(:,i+1,j+1),xcg(:,i  ,j+1),a2,imf2)
        call alpha(xst,xcg(:,i  ,j+1),xcg(:,i  ,j  ),a3,imf3)
        call alpha(xst,xcg(:,i  ,j  ),xcg(:,i+1,j  ),a4,imf4)
        if (imf1*imf2*imf3*imf4.ne.0) then ! angoli sempre definiti
         atot = a1+a2+a3+a4
         if (abs(atot-2.*pi).lt.1.d-6) then
          ii = i
          jj = j
          ierr = 0
          exit
         endif
        else ! punto da trovare coincidente con un centro cella
         if ((xst(1).eq.xcg(1,i+1,j)).and.(xst(2).eq.xcg(2,i+1,j))) then
          ii = i+1
          jj = j
          ierr = 0
          exit
         elseif ((xst(1).eq.xcg(1,i+1,j+1)).and. &
           (xst(2).eq.xcg(2,i+1,j+1))) then
          ii = i+1
          jj = j+1
          ierr = 0
          exit
         elseif ((xst(1).eq.xcg(1,i,j+1)).and. &
           (xst(2).eq.xcg(2,i,j+1))) then
          ii = i
          jj = j+1
          ierr = 0
          exit
         elseif ((xst(1).eq.xcg(1,i,j)).and. &
           (xst(2).eq.xcg(2,i,j))) then
          ii = i
          jj = j
          ierr = 0
          exit
         endif
        endif
!
       enddo
      enddo
!
      return
      end
      !}}}
      !{{{ solid_search
      subroutine solid_search(i1,i2,j1,j2, &
                           xcg,istart,iend,ny,ng, xst, &
                           ii,jj,ierr,iprint) 
!
      implicit double precision(a-h,o-z)
!
      dimension xcg(2,istart-ng:iend+ng,1-ng:ny+ng)
      dimension xvec(2,istart-ng:iend+ng)
      dimension yvec(2,1-ng:ny+ng)
      dimension xst(2)
      character chit*3
!
      pi = 4.*atan(1.)
!
!     REFINED SEARCH (it should be a good compromise)
!
!     checking edge of half domain, to see where it is confined the
!     point, then proceeding by bisection
      ierr = 1
!
!     step 1: cerco la i t.c. la distanza del punto dai centricella sul
!     profilo sia minima; questa è solo un'inizializzazione
      xvec(:,:) = xcg(:,istart-ng:iend+ng,j1)
      call distmin(xst,xvec(:,:),istart,iend,ng,i1,i2,imin)
      im = imin
      if (iprint.eq.1) print*, "starting im=", im
!
!     step 2: è la parte iterativa:
!     cerco di capire da quale lato sto della linea a i costante
!     trovata. lo faccio guardando al putno a distanza minima.
!     sostituisco i1 o i2 e continuo a bisezionare: trovo la fila i di
!     celle che contengono il punto
!
      it=0
      if (iprint.eq.1) then
       open(66,file='debug_point.dat')
       write(66,*) xst
       close(66)
      endif
 100  continue
      it=it+1
      if (iprint.eq.1) then
       write(chit,'(I3.3)') it
       open(66,file='debug_'//chit//'.dat')
       do i=i1,i2
        write(66,*) xcg(:,i,j1)
       enddo
       do j=j1,j2
        write(66,*) xcg(:,i2,j)
       enddo
       do i=i2,i1,-1
        write(66,*) xcg(:,i,j2)
       enddo
       do j=j2,j1,-1
        write(66,*) xcg(:,i1,j)
       enddo
       close(66)
       print*, 'i1,i2,im, j1,j2', i1,i2,im, j1,j2
      endif
!
!     considero solo la metà destra. se non sta qua, allora sta a
!     sinistra
!
!     fila destra (i+1, j, j+1)
      alpha1 = 0.
      do j=j1,j2-1
       call alpha(xst,xcg(:,i2,j),xcg(:,i2,j+1),a1,imf1)
       if (imf1.ne.0) then ! angoli sempre definiti
        alpha1 = alpha1+a1
       else ! punto da trovare coincidente con un centro cella
        if ((xst(1).eq.xcg(1,i2,j)).and. &
      (xst(2).eq.xcg(2,i2,j))) then
         ii = i2
         jj = j
         ierr = 0
         goto 666
        elseif ((xst(1).eq.xcg(1,i2,j+1)).and. &
          (xst(2).eq.xcg(2,i2,j+1))) then
         ii = i2
         jj = j+1
         ierr = 0
         goto 666
        endif
       endif
      enddo
!
!     fila sopra (j+1, i+1, i)
      alpha2 = 0.
      do i=i2-1,im,-1
       call alpha(xst,xcg(:,i+1,j2),xcg(:,i,j2),a2,imf2)
       if (imf2.ne.0) then ! angoli sempre definiti
        alpha2 = alpha2+a2
       else ! punto da trovare coincidente con un centro cella
        if ((xst(1).eq.xcg(1,i,j2)).and. &
      (xst(2).eq.xcg(2,i,j2))) then
         ii = i
         jj = j2
         ierr = 0
         goto 666
        elseif ((xst(1).eq.xcg(1,i+1,j2)).and. &
          (xst(2).eq.xcg(2,i+1,j2))) then
         ii = i+1
         jj = j2
         ierr = 0
         goto 666
        endif
       endif
      enddo
!
!     fila sinistra (im, j+1, j)
      alpha3 = 0.
      do j=j2-1,j1,-1
       call alpha(xst,xcg(:,im,j+1),xcg(:,im,j),a3,imf3)
       if (imf3.ne.0) then ! angoli sempre definiti
        alpha3 = alpha3+a3
       else ! punto da trovare coincidente con un centro cella
        if ((xst(1).eq.xcg(1,im,j)).and. &
      (xst(2).eq.xcg(2,im,j))) then
         ii = im
         jj = j
         ierr = 0
         goto 666
        elseif ((xst(1).eq.xcg(1,im,j+1)).and. &
          (xst(2).eq.xcg(2,im,j+1))) then
         ii = im
         jj = j+1
         ierr = 0
         goto 666
        endif
       endif
      enddo
!
!     fila sotto (j1, i, i+1)
      alpha4 = 0.
      do i=im,i2-1
       call alpha(xst,xcg(:,i,j1),xcg(:,i+1,j1),a4,imf4)
       if (imf4.ne.0) then ! angoli sempre definiti
        alpha4 = alpha4+a4
       else ! punto da trovare coincidente con un centro cella
        if ((xst(1).eq.xcg(1,i,j1)).and. &
      (xst(2).eq.xcg(2,i,j1))) then
         ii = i
         jj = j1
         ierr = 0
         goto 666
        elseif ((xst(1).eq.xcg(1,i+1,j1)).and. &
          (xst(2).eq.xcg(2,i+1,j1))) then
         ii = i+1
         jj = j1
         ierr = 0
         goto 666
        endif
       endif
      enddo
!
!     se l'angolo è 360, allora im sta a sinistra, senno a destra
      atot = alpha1 + alpha2 + alpha3 + alpha4
      if (abs(atot-2.*pi).lt.1.d-6) then ! im sta a sinistra
       i1 = im
      else ! im sta a destra
       i2 = im
      endif
      im = (i1+i2)/2
      if (i2-i1.gt.1) goto 100
!
!     step 3: trovo la j, con la ricerca super fine, ma tra gli indici
!     nuovi
      call robust_search(i1,i2,j1,j2, &
                            xcg(:,istart:iend,1:ny), &
                            istart,iend,ny,xst,ii,jj,ierr)
!
      if ((ierr.eq.1).and.(iprint.eq.1)) then
       print*, '================================'
       print*, 'from locint:'
       print*, 'i1, i2', i1, i2
       print*, 'j1, j2', j1, j2
       print*, 'nx, ny', iend-istart+1, ny
       print*, 'ii, jj', ii, jj
       print*, 'xst', xst
       print*, 'xcg(:,ii,jj+1)  xcg(:,ii+1,jj+1)'
       print*, 'xcg(:,ii,jj  )  xcg(:,ii+1,jj  )'
       print*,  xcg(:,ii,jj+1), xcg(:,ii+1,jj+1)
       print*,  xcg(:,ii,jj  ), xcg(:,ii+1,jj  )
       print*, '================================'
       print*, 'warning, point inside the domain'
      endif
!
!     check
      call alpha(xst,xcg(:,ii+1,jj  ),xcg(:,ii+1,jj+1),a1,imf1)
      call alpha(xst,xcg(:,ii+1,jj+1),xcg(:,ii  ,jj+1),a2,imf2)
      call alpha(xst,xcg(:,ii  ,jj+1),xcg(:,ii  ,jj  ),a3,imf3)
      call alpha(xst,xcg(:,ii  ,jj  ),xcg(:,ii+1,jj  ),a4,imf4)
      if (imf1*imf2*imf3*imf4.ne.0) then ! angoli sempre definiti
       pi = 4.*atan(1.)
       atot = a1+a2+a3+a4
       if (abs(atot-2.*pi).lt.1.d-6) then
        ierr = 0
!       print*, 'a refined search is required'
!       print*, 2.*pi, atot
       endif
      endif
!
 666  continue
 1003 format(I3.3)
!
      return
      end
      !}}}
