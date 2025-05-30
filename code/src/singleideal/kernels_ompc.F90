module streams_kernels_ompc

  use streams_parameters, only : rkind, ikind, real64

  implicit none

contains

  subroutine zero_flux_subroutine(nx,ny,nz,nv,fl_ompc)
    integer :: nx, ny, nz, nv
    real(rkind), dimension(1:,1:,1:,1:), intent(inout) :: fl_ompc
    integer :: i,j,k,m,iercuda

    !$omp parallel do default(firstprivate) shared(fl_ompc)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          do m=1,nv
            fl_ompc(i,j,k,m) = 0._rkind
          enddo
        enddo
      enddo
    enddo


  endsubroutine zero_flux_subroutine



  subroutine init_flux_subroutine(nx,ny,nz,nv,fl_ompc,fln_ompc,rhodt)
    integer :: nx, ny, nz, nv
    real(rkind) :: rhodt
    real(rkind), dimension(1:,1:,1:,1:), intent(inout) :: fl_ompc, fln_ompc
    integer :: i,j,k,m,iercuda

    !$omp parallel do default(firstprivate) shared(fl_ompc, &
    !$omp& fln_ompc)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          do m=1,nv
            fln_ompc(i,j,k,m) = - rhodt * fl_ompc(i,j,k,m)
            fl_ompc(i,j,k,m) = 0._rkind
          enddo
        enddo
      enddo
    enddo


  endsubroutine init_flux_subroutine




  subroutine count_weno_subroutine(nv,nv_aux,nx,ny,nz,ng,eul_imin,eul_imax,eul_jmin,eul_jmax,&
  &eul_kmin,eul_kmax,weno_scheme,sensor_threshold,w_aux_ompc,ep_ord_change_ompc,count_weno_x,&
  &count_weno_y,count_weno_z)
    integer, intent(in) :: nv, nv_aux, nx, ny, nz, ng
    integer, intent(in) :: eul_imin, eul_imax, eul_jmin, eul_jmax, eul_kmin, eul_kmax
    integer, intent(in) :: weno_scheme
    real(rkind) :: sensor_threshold
    real(rkind), intent(out) :: count_weno_x, count_weno_y, count_weno_z
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv_aux) :: w_aux_ompc
    integer, dimension(0:nx,0:ny,0:nz,1:3) :: ep_ord_change_ompc
    integer :: i,j,k,ishk,ii,jj,kk,iercuda
    count_weno_x = 0._rkind
    !$omp parallel do default(firstprivate) shared(w_aux_ompc, &
    !$omp& ep_ord_change_ompc) &
    !$omp& reduction(+:count_weno_x)
    do k = 1,nz
      do j = 1,ny
        do i = eul_imin,eul_imax
          ishk = 0
          do ii=i-weno_scheme+1,i+weno_scheme
            if (w_aux_ompc(ii,j,k,8) > sensor_threshold) then
              ishk = 1
            endif
          enddo
          if (ishk > 0) then
            count_weno_x = count_weno_x + 1
          endif
        enddo
      enddo
    enddo
    count_weno_y = 0._rkind
    !$omp parallel do default(firstprivate) shared(w_aux_ompc, &
    !$omp& ep_ord_change_ompc) &
    !$omp& reduction(+:count_weno_y)
    do k = 1,nz
      do j = eul_jmin,eul_jmax
        do i = 1,nx
          ishk = 0
          do jj=j-weno_scheme+1,j+weno_scheme
            if (w_aux_ompc(i,jj,k,8) > sensor_threshold) then
              ishk = 1
            endif
          enddo
          if (ishk > 0) then
            count_weno_y = count_weno_y + 1
          endif
        enddo
      enddo
    enddo
    count_weno_z = 0._rkind
    !$omp parallel do default(firstprivate) shared(w_aux_ompc, &
    !$omp& ep_ord_change_ompc) &
    !$omp& reduction(+:count_weno_z)
    do k = eul_kmin,eul_kmax
      do j = 1,ny
        do i = 1,nx
          ishk = 0
          do kk=k-weno_scheme+1,k+weno_scheme
            if (w_aux_ompc(i,j,kk,8) > sensor_threshold) then
              ishk = 1
            endif
          enddo
          if (ishk > 0) then
            count_weno_z = count_weno_z + 1
          endif
        enddo
      enddo
    enddo


  endsubroutine count_weno_subroutine



  subroutine count_weno_c2_subroutine(nv,nv_aux,nx,ny,nz,ng,eul_imin,eul_imax,eul_jmin,eul_jmax,&
  &eul_kmin,eul_kmax,l_base,weno_scheme,sensor_threshold,w_aux_ompc,ep_ord_change_ompc,lmax_tag_ompc,&
  &wall_tag_ompc,count_weno_x,count_weno_y,count_weno_z)
    integer, intent(in) :: nv, nv_aux, nx, ny, nz, ng
    integer, intent(in) :: eul_imin, eul_imax, eul_jmin, eul_jmax, eul_kmin, eul_kmax
    integer, intent(in) :: l_base, weno_scheme
    real(rkind) :: sensor_threshold
    real(rkind), intent(out) :: count_weno_x, count_weno_y, count_weno_z
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv_aux) :: w_aux_ompc
    integer, dimension(0:nx,0:ny,0:nz,1:3) :: ep_ord_change_ompc
    integer, dimension(1-ng:nx+ng) :: lmax_tag_ompc
    integer, dimension(1-ng:nx+ng), intent(in) :: wall_tag_ompc
    integer :: i, j, k, ishk, weno_scheme_i, weno_scheme_j, weno_scheme_k, ii, jj, kk, iercuda
    count_weno_x = 0._rkind
    !$omp parallel do default(firstprivate) shared(w_aux_ompc, &
    !$omp& ep_ord_change_ompc, &
    !$omp& lmax_tag_ompc, &
    !$omp& wall_tag_ompc) &
    !$omp& reduction(+:count_weno_x)
    do k = 1,nz
      do j = 1,ny
        do i = eul_imin,eul_imax
          weno_scheme_i = weno_scheme+ep_ord_change_ompc(i,j,k,1)
          if(weno_scheme_i < 1) then
            weno_scheme_i = 1
          endif
          if(j == 1) then
            weno_scheme_i = lmax_tag_ompc(i)+ep_ord_change_ompc(i,j,k,1)
            if(weno_scheme_i < 1) then
              weno_scheme_i = 1
            endif
          endif
          ishk = 0
          do ii=i-weno_scheme_i+1,i+weno_scheme_i
            if (w_aux_ompc(ii,j,k,8) > sensor_threshold) then
              ishk = 1
            endif
          enddo
          if (ishk > 0) then
            count_weno_x = count_weno_x + 1
          endif
        enddo
      enddo
    enddo
    count_weno_y = 0._rkind
    !$omp parallel do default(firstprivate) shared(w_aux_ompc, &
    !$omp& ep_ord_change_ompc, &
    !$omp& lmax_tag_ompc, &
    !$omp& wall_tag_ompc) &
    !$omp& reduction(+:count_weno_y)
    do k = 1,nz
      do j = eul_jmin,eul_jmax
        do i = 1,nx
          weno_scheme_j = weno_scheme+ep_ord_change_ompc(i,j,k,2)
          if(weno_scheme_j < 1) then
            weno_scheme_j = 1
          endif
          if (j <= l_base .and. wall_tag_ompc(i) > 0) then
            weno_scheme_j = weno_scheme
          endif
          ishk = 0
          do jj=j-weno_scheme_j+1,j+weno_scheme_j
            if (w_aux_ompc(i,jj,k,8) > sensor_threshold) then
              ishk = 1
            endif
          enddo
          if (ishk > 0) then
            count_weno_y = count_weno_y + 1
          endif
        enddo
      enddo
    enddo
    count_weno_z = 0._rkind
    !$omp parallel do default(firstprivate) shared(w_aux_ompc, &
    !$omp& ep_ord_change_ompc, &
    !$omp& lmax_tag_ompc, &
    !$omp& wall_tag_ompc) &
    !$omp& reduction(+:count_weno_z)
    do k = eul_kmin,eul_kmax
      do j = 1,ny
        do i = 1,nx
          ishk = 0
          do kk=k-weno_scheme+1,k+weno_scheme
            if (w_aux_ompc(i,j,kk,8) > sensor_threshold) then
              ishk = 1
            endif
          enddo
          if (ishk > 0) then
            count_weno_z = count_weno_z + 1
          endif
        enddo
      enddo
    enddo


  endsubroutine count_weno_c2_subroutine



  subroutine euler_x_fluxes_hybrid_c2_kernel(nv,nv_aux,nx,ny,nz,ng,eul_imin,eul_imax,lmax_base,&
  &nkeep,rgas0,w_aux_ompc,coeff_deriv1_ompc,fhat_ompc,dcsidxc2_ompc,dcsidxnc2_ompc,dcsidyc2_ompc,&
  &dcsidync2_ompc,jac_ompc,mcsijac1_ompc,lmax_tag_ompc,force_zero_flux_min,force_zero_flux_max,&
  &weno_scheme,weno_version,sensor_threshold,weno_size,cp_coeff_ompc,indx_cp_l,indx_cp_r,&
  &ep_ord_change_ompc,calorically_perfect,tol_iter_nr,rho0,u0,t0)
    integer, value :: nv, nx, ny, nz, ng, nv_aux
    integer, value :: eul_imin, eul_imax, lmax_base, nkeep, indx_cp_l, indx_cp_r, calorically_perfect
    integer, value :: weno_scheme, weno_size, weno_version
    integer, dimension(0:nx,0:ny,0:nz,1:3) :: ep_ord_change_ompc
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv_aux) :: w_aux_ompc
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv) :: fhat_ompc
    real(rkind), dimension(indx_cp_l:indx_cp_r+1) :: cp_coeff_ompc
    real(rkind), dimension(4,4) :: coeff_deriv1_ompc
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng) :: dcsidxc2_ompc
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng) :: dcsidxnc2_ompc
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng) :: dcsidyc2_ompc
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng) :: dcsidync2_ompc
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng) :: jac_ompc
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng) :: mcsijac1_ompc
    integer, dimension(1-ng:nx+ng) :: lmax_tag_ompc
    integer, value :: force_zero_flux_min, force_zero_flux_max
    real(rkind), value :: sensor_threshold, rgas0, tol_iter_nr,rho0,u0,t0
    integer :: i, j, k, m, l
    real(rkind) :: fh1, fh2, fh3, fh4, fh5
    real(rkind) :: rhom, uui, uti, vvi, wwi, ppi, enti, rhoi, tti, dcsidxi, dcsidyi
    real(rkind) :: uuip, utip, vvip, wwip, ppip, entip, rhoip, ttip, dcsidxip, dcsidyip
    real(rkind) :: ft1, ft2, ft3, ft4, ft5, ft6, ft7
    real(rkind) :: uvs1, uvs2, uvs3, uvs4, uvs5, uvs6, uvs7, uv_part
    integer :: ii, lmax, wenorec_ord, lmaxi, weno_scheme_i, weno_size_i
    integer :: ishk
    real(rkind) :: b1, b2, b3, c, ci, h, uu, vv, ww
    real(rkind), dimension(5,5) :: el, er
    real(rkind), dimension(5) :: evmax, fi
    real(rkind), dimension(5,8) :: gp,gm
    integer :: ll, mm
    real(rkind) :: rho, pp, wc, gc, rhou, ut, rhov, rhow, rhoh
    real(rkind) :: dcsidxav, dcsidyav, gcsixav, gcsiyav, dcsimod, qq
    real(rkind) :: tt, gamloc
    real(rkind) :: uvs5_i,uvs5_k,uvs5_p,eei,eeip
    real(rkind) :: drho, dee, eem
    real(rkind) :: drhof, deef
    real(rkind) :: sumnumrho,sumnumee,sumdenrho,sumdenee
    integer :: n,n2
    real(rkind) :: t_sumdenrho, t_sumdenee, t2_sumdenrho, t2_sumdenee
    integer :: iii,jjj,kkk
    real(rkind), dimension(5,5) :: eler

    !$omp parallel do default(firstprivate) shared(w_aux_ompc, &
    !$omp& fhat_ompc, &
    !$omp& cp_coeff_ompc, &
    !$omp& coeff_deriv1_ompc, &
    !$omp& dcsidxc2_ompc, &
    !$omp& dcsidxnc2_ompc, &
    !$omp& dcsidyc2_ompc, &
    !$omp& dcsidync2_ompc, &
    !$omp& jac_ompc, &
    !$omp& mcsijac1_ompc, &
    !$omp& ep_ord_change_ompc, &
    !$omp& lmax_tag_ompc)
    do k = 1,nz
      do j = 1,ny
        do i = +eul_imin-2+1,eul_imax
          weno_scheme_i = max(weno_scheme+ep_ord_change_ompc(i,j,k,1),1)
          lmax = max(lmax_base+ep_ord_change_ompc(i,j,k,1),1)
          if(j == 1) then
            weno_scheme_i = max(lmax_tag_ompc(i)+ep_ord_change_ompc(i,j,k,1),1)
            lmax = weno_scheme_i
          endif
          weno_size = 2*weno_scheme_i

          ishk = 0
          do ii=i-weno_scheme_i+1,i+weno_scheme_i
            if (w_aux_ompc(ii,j,k,8) > sensor_threshold) ishk = 1
          enddo

          if (ishk == 0) then

            ft1 = 0._rkind
            ft2 = 0._rkind
            ft3 = 0._rkind
            ft4 = 0._rkind
            ft5 = 0._rkind
            ft6 = 0._rkind
            ft7 = 0._rkind
            do l=1,lmax
              uvs1 = 0._rkind
              uvs2 = 0._rkind
              uvs3 = 0._rkind
              uvs4 = 0._rkind
              uvs5 = 0._rkind
              uvs6 = 0._rkind
              uvs7 = 0._rkind
              do m=0,l-1

                rhoi = w_aux_ompc(i-m,j,k,1)
                uui = w_aux_ompc(i-m,j,k,2)
                vvi = w_aux_ompc(i-m,j,k,3)
                wwi = w_aux_ompc(i-m,j,k,4)
                enti = w_aux_ompc(i-m,j,k,5)
                tti = w_aux_ompc(i-m,j,k,6)
                ppi = tti*rhoi*rgas0
                uti = (uui*dcsidxc2_ompc(i-m,j) + vvi*dcsidyc2_ompc(i-m,j))/jac_ompc(i-m,j)
                dcsidxi = dcsidxc2_ompc(i-m,j)/jac_ompc(i-m,j)
                dcsidyi = dcsidyc2_ompc(i-m,j)/jac_ompc(i-m,j)

                rhoip = w_aux_ompc(i-m+l,j,k,1)
                uuip = w_aux_ompc(i-m+l,j,k,2)
                vvip = w_aux_ompc(i-m+l,j,k,3)
                wwip = w_aux_ompc(i-m+l,j,k,4)
                entip = w_aux_ompc(i-m+l,j,k,5)
                ttip = w_aux_ompc(i-m+l,j,k,6)
                ppip = ttip*rhoip*rgas0
                utip = (uuip*dcsidxc2_ompc(i-m+l,j) + vvip*dcsidyc2_ompc(i-m+l,j))/jac_ompc(i-m+l,j)
                dcsidxip = dcsidxc2_ompc(i-m+l,j)/jac_ompc(i-m+l,j)
                dcsidyip = dcsidyc2_ompc(i-m+l,j)/jac_ompc(i-m+l,j)

                rhom = rhoi+rhoip
                uv_part = (uti+utip) * rhom
                uvs1 = uvs1 + uv_part * (2._rkind)
                uvs2 = uvs2 + uv_part * (uui+uuip)
                uvs3 = uvs3 + uv_part * (vvi+vvip)
                uvs4 = uvs4 + uv_part * (wwi+wwip)
                uvs5 = uvs5 + uv_part * (enti+entip)
                uvs6 = uvs6 + (2._rkind)*(ppi+ppip)*(dcsidxi+dcsidxip)
                uvs7 = uvs7 + (2._rkind)*(ppi+ppip)*(dcsidyi+dcsidyip)
              enddo
              ft1 = ft1 + coeff_deriv1_ompc(l,lmax)*uvs1
              ft2 = ft2 + coeff_deriv1_ompc(l,lmax)*uvs2
              ft3 = ft3 + coeff_deriv1_ompc(l,lmax)*uvs3
              ft4 = ft4 + coeff_deriv1_ompc(l,lmax)*uvs4
              ft5 = ft5 + coeff_deriv1_ompc(l,lmax)*uvs5
              ft6 = ft6 + coeff_deriv1_ompc(l,lmax)*uvs6
              ft7 = ft7 + coeff_deriv1_ompc(l,lmax)*uvs7
            enddo

            fh1 = 0.25_rkind*ft1
            fh2 = 0.25_rkind*ft2
            fh3 = 0.25_rkind*ft3
            fh4 = 0.25_rkind*ft4
            fh5 = 0.25_rkind*ft5

            if ((i==0 .and. force_zero_flux_min == 1).or.(i==nx .and. force_zero_flux_max == 1)) then
              fh1 = 0._rkind
              fh2 = 0._rkind
              fh3 = 0._rkind
              fh4 = 0._rkind
              fh5 = 0._rkind
            endif
            fh2 = fh2 + 0.25_rkind*ft6
            fh3 = fh3 + 0.25_rkind*ft7

            fhat_ompc(i,j,k,1) = fh1
            fhat_ompc(i,j,k,2) = fh2
            fhat_ompc(i,j,k,3) = fh3
            fhat_ompc(i,j,k,4) = fh4
            fhat_ompc(i,j,k,5) = fh5
          else
            call compute_roe_average(nx, ny, nz, ng, i, i+1, j, j, k, k, w_aux_ompc, rgas0, b1, b2,&
            & b3, c, ci, h, uu, vv, ww, cp_coeff_ompc, indx_cp_l, indx_cp_r, calorically_perfect, tol_iter_nr,&
            &t0)
            dcsidxav = 0.5_rkind * (dcsidxnc2_ompc(i,j) + dcsidxnc2_ompc(i+1,j))
            dcsidyav = 0.5_rkind * (dcsidync2_ompc(i,j) + dcsidync2_ompc(i+1,j))
            dcsimod = 1._rkind/sqrt(dcsidxav**2+dcsidyav**2)
            dcsidxav = dcsidxav * dcsimod
            dcsidyav = dcsidyav * dcsimod
            ut = dcsidxav * uu + dcsidyav * vv

            call eigenvectors_x_c2(b1, b2, b3, uu, vv, ww, c, ci, h, el, er, ut, dcsidxav, dcsidyav)

            do m=1,5
              evmax(m) = -1._rkind
            enddo
            do l=1,weno_size
              ll = i + l - weno_scheme_i
              uu = w_aux_ompc(ll,j,k,2)
              vv = w_aux_ompc(ll,j,k,3)
              ut = dcsidxnc2_ompc(ll,j)*uu+dcsidync2_ompc(ll,j)*vv
              tt = w_aux_ompc(ll,j,k,6)
              gamloc = get_gamloc_dev(indx_cp_l,indx_cp_r,tt,t0,cp_coeff_ompc,calorically_perfect,rgas0)
              c = sqrt (gamloc*rgas0*tt)

              evmax(1) = max(abs(ut-c),evmax(1))
              evmax(2) = max(abs(ut ),evmax(2))
              evmax(3) = evmax(2)
              evmax(4) = evmax(2)
              evmax(5) = max(abs(ut+c),evmax(5))
            enddo
            do l=1,weno_size
              ll = i + l - weno_scheme_i

              rho = w_aux_ompc(ll,j,k,1)
              uu = w_aux_ompc(ll,j,k,2)
              vv = w_aux_ompc(ll,j,k,3)
              ut = dcsidxnc2_ompc(ll,j)*uu+dcsidync2_ompc(ll,j)*vv
              ww = w_aux_ompc(ll,j,k,4)
              h = w_aux_ompc(ll,j,k,5)
              pp = rho*w_aux_ompc(ll,j,k,6)*rgas0
              fi(1) = rho * ut
              fi(2) = (rho * uu * ut + pp * dcsidxnc2_ompc(ll,j))
              fi(3) = (rho * vv * ut + pp * dcsidync2_ompc(ll,j))
              fi(4) = rho * ww * ut
              fi(5) = rho * h * ut
              do m=1,5
                wc = 0._rkind
                gc = 0._rkind

                wc = wc + el(1,m) * rho * mcsijac1_ompc(ll,j)
                gc = gc + el(1,m) * fi(1) * mcsijac1_ompc(ll,j)
                wc = wc + el(2,m) * rho*uu * mcsijac1_ompc(ll,j)
                gc = gc + el(2,m) * fi(2) * mcsijac1_ompc(ll,j)
                wc = wc + el(3,m) * rho*vv * mcsijac1_ompc(ll,j)
                gc = gc + el(3,m) * fi(3) * mcsijac1_ompc(ll,j)
                wc = wc + el(4,m) * rho*ww * mcsijac1_ompc(ll,j)
                gc = gc + el(4,m) * fi(4) * mcsijac1_ompc(ll,j)
                wc = wc + el(5,m) * (rho*h-pp) * mcsijac1_ompc(ll,j)
                gc = gc + el(5,m) * fi(5) * mcsijac1_ompc(ll,j)

                c = 0.5_rkind * (gc + evmax(m) * wc)
                gp(m,l) = c
                gm(m,l) = gc - c
              enddo
            enddo
            call wenorec_1d(nv,gp,gm,fi,weno_scheme_i,weno_scheme_i,weno_version,rho0,u0)
            do m=1,5
              fhat_ompc(i,j,k,m) = 0._rkind
              do mm=1,5
                fhat_ompc(i,j,k,m) = fhat_ompc(i,j,k,m) + er(mm,m) * fi(mm)
              enddo
            enddo

          endif
        enddo
      enddo
    enddo


  endsubroutine euler_x_fluxes_hybrid_c2_kernel


  subroutine euler_x_fluxes_hybrid_kernel(nv,nv_aux,nx,ny,nz,ng,istart_face,iend_face,lmax_base,&
  &nkeep,rgas0,w_aux_ompc,coeff_deriv1_ompc,dcsidx_ompc,fhat_ompc,force_zero_flux_min,&
  &force_zero_flux_max,weno_scheme,weno_version,sensor_threshold,weno_size,cp_coeff_ompc,indx_cp_l,&
  &indx_cp_r,ep_ord_change_ompc,calorically_perfect,tol_iter_nr,rho0,u0,t0)
    integer, value :: nv, nx, ny, nz, ng, nv_aux
    integer, value :: istart_face, iend_face, lmax_base, nkeep, indx_cp_l, indx_cp_r, calorically_perfect
    integer, value :: weno_scheme, weno_size, weno_version
    integer, dimension(0:nx,0:ny,0:nz,1:3) :: ep_ord_change_ompc
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv_aux) :: w_aux_ompc
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv) :: fhat_ompc
    real(rkind), dimension(indx_cp_l:indx_cp_r+1) :: cp_coeff_ompc
    real(rkind), dimension(4,4) :: coeff_deriv1_ompc
    real(rkind), dimension(nx) :: dcsidx_ompc
    integer, value :: force_zero_flux_min, force_zero_flux_max
    real(rkind), value :: sensor_threshold, rgas0, tol_iter_nr,rho0,u0,t0
    integer :: i, j, k, m, l
    real(rkind) :: fh1, fh2, fh3, fh4, fh5
    real(rkind) :: rhom, uui, vvi, wwi, ppi, enti, rhoi, tti
    real(rkind) :: uuip, vvip, wwip, ppip, entip, rhoip, ttip
    real(rkind) :: ft1, ft2, ft3, ft4, ft5, ft6
    real(rkind) :: uvs1, uvs2, uvs3, uvs4, uvs6, uv_part
    real(rkind) :: uvs5
    integer :: ii, lmax, wenorec_ord
    integer :: ishk
    real(rkind) :: b1, b2, b3, c, ci, h, uu, vv, ww
    real(rkind), dimension(5,5) :: el, er
    real(rkind), dimension(5) :: evmax, fi
    real(rkind), dimension(5,8) :: gp,gm
    integer :: ll, mm
    real(rkind) :: rho, pp, wc, gc, rhou
    real(rkind) :: tt, gamloc
    real(rkind) :: uvs5_i,uvs5_k,uvs5_p,eei,eeip
    real(rkind) :: drho, dee, eem
    real(rkind) :: drhof, deef
    real(rkind) :: sumnumrho,sumnumee,sumdenrho,sumdenee
    integer :: n,n2
    real(rkind) :: t_sumdenrho, t_sumdenee, t2_sumdenrho, t2_sumdenee

    !$omp parallel do default(firstprivate) shared(w_aux_ompc, &
    !$omp& fhat_ompc, &
    !$omp& cp_coeff_ompc, &
    !$omp& coeff_deriv1_ompc, &
    !$omp& dcsidx_ompc, &
    !$omp& ep_ord_change_ompc)
    do k = 1,nz
      do j = 1,ny
        do i = +istart_face-1+1,iend_face
          ishk = 0
          do ii=i-weno_scheme+1,i+weno_scheme
            if (w_aux_ompc(ii,j,k,8) > sensor_threshold) ishk = 1
          enddo

          if (ishk == 0) then

            ft1 = 0._rkind
            ft2 = 0._rkind
            ft3 = 0._rkind
            ft4 = 0._rkind
            ft5 = 0._rkind
            ft6 = 0._rkind
            lmax = max(lmax_base+ep_ord_change_ompc(i,j,k,1),1)
            if (nkeep>=0) then
              do l=1,lmax
                uvs1 = 0._rkind
                uvs2 = 0._rkind
                uvs3 = 0._rkind
                uvs4 = 0._rkind
                uvs5_i = 0._rkind
                uvs5_k = 0._rkind
                uvs5_p = 0._rkind
                uvs6 = 0._rkind
                do m=0,l-1

                  rhoi = w_aux_ompc(i-m,j,k,1)
                  uui = w_aux_ompc(i-m,j,k,2)
                  vvi = w_aux_ompc(i-m,j,k,3)
                  wwi = w_aux_ompc(i-m,j,k,4)
                  enti = w_aux_ompc(i-m,j,k,5)
                  tti = w_aux_ompc(i-m,j,k,6)
                  ppi = tti*rhoi*rgas0
                  eei = enti-ppi/rhoi-0.5_rkind*(uui*uui+vvi*vvi+wwi*wwi)

                  rhoip = w_aux_ompc(i-m+l,j,k,1)
                  uuip = w_aux_ompc(i-m+l,j,k,2)
                  vvip = w_aux_ompc(i-m+l,j,k,3)
                  wwip = w_aux_ompc(i-m+l,j,k,4)
                  entip = w_aux_ompc(i-m+l,j,k,5)
                  ttip = w_aux_ompc(i-m+l,j,k,6)
                  ppip = ttip*rhoip*rgas0
                  eeip = entip-ppip/rhoip-0.5_rkind*(uuip*uuip+vvip*vvip+wwip*wwip)

                  rhom = rhoi+rhoip
                  eem = eei + eeip

                  if(nkeep == 0) then
                    drhof = 1._rkind
                    deef = 1._rkind
                  else
                    sumnumrho = 1._rkind
                    drho = 2._rkind*(rhoip-rhoi)/rhom
                    dee = 2._rkind*(eeip - eei)/eem
                    t_sumdenrho = (0.5_rkind*drho)*(0.5_rkind*drho)
                    t_sumdenee = (0.5_rkind*dee )*(0.5_rkind*dee )
                    t2_sumdenrho = t_sumdenrho
                    t2_sumdenee = t_sumdenee
                    sumdenrho = 1._rkind + t_sumdenrho / (3._rkind)
                    sumdenee = 1._rkind + t_sumdenee
                    sumnumee = 1._rkind + t_sumdenee / (3._rkind)
                    do n = 2, nkeep
                      n2 = 2*n
                      t_sumdenrho = t2_sumdenrho * t_sumdenrho
                      t_sumdenee = t2_sumdenee * t_sumdenee
                      sumdenrho = sumdenrho + t_sumdenrho / (1._rkind+n2)
                      sumdenee = sumdenee + t_sumdenee
                      sumnumee = sumnumee + t_sumdenee / (1._rkind+n2)
                    enddo
                    drhof = sumnumrho/sumdenrho
                    deef = sumnumee /sumdenee
                  endif

                  uv_part = (uui+uuip) * rhom * drhof
                  uvs1 = uvs1 + uv_part * (2._rkind)
                  uvs2 = uvs2 + uv_part * (uui+uuip)
                  uvs3 = uvs3 + uv_part * (vvi+vvip)
                  uvs4 = uvs4 + uv_part * (wwi+wwip)
                  uvs5_i = uvs5_i + uv_part * eem * deef
                  uvs5_k = uvs5_k + uv_part * (uui*uuip+vvi*vvip+wwi*wwip)
                  uvs5_p = uvs5_p + 4._rkind*(uui*ppip+uuip*ppi)
                  uvs6 = uvs6 + (2._rkind)*(ppi+ppip)
                enddo
                ft1 = ft1 + coeff_deriv1_ompc(l,lmax)*uvs1
                ft2 = ft2 + coeff_deriv1_ompc(l,lmax)*uvs2
                ft3 = ft3 + coeff_deriv1_ompc(l,lmax)*uvs3
                ft4 = ft4 + coeff_deriv1_ompc(l,lmax)*uvs4
                ft5 = ft5 + coeff_deriv1_ompc(l,lmax)*(uvs5_i+uvs5_k+uvs5_p)
                ft6 = ft6 + coeff_deriv1_ompc(l,lmax)*uvs6
              enddo
            else
              do l=1,lmax
                uvs1 = 0._rkind
                uvs2 = 0._rkind
                uvs3 = 0._rkind
                uvs4 = 0._rkind
                uvs5 = 0._rkind
                uvs6 = 0._rkind
                do m=0,l-1

                  rhoi = w_aux_ompc(i-m,j,k,1)
                  uui = w_aux_ompc(i-m,j,k,2)
                  vvi = w_aux_ompc(i-m,j,k,3)
                  wwi = w_aux_ompc(i-m,j,k,4)
                  enti = w_aux_ompc(i-m,j,k,5)
                  tti = w_aux_ompc(i-m,j,k,6)
                  ppi = tti*rhoi*rgas0

                  rhoip = w_aux_ompc(i-m+l,j,k,1)
                  uuip = w_aux_ompc(i-m+l,j,k,2)
                  vvip = w_aux_ompc(i-m+l,j,k,3)
                  wwip = w_aux_ompc(i-m+l,j,k,4)
                  entip = w_aux_ompc(i-m+l,j,k,5)
                  ttip = w_aux_ompc(i-m+l,j,k,6)
                  ppip = ttip*rhoip*rgas0

                  rhom = rhoi+rhoip
                  uv_part = (uui+uuip) * rhom
                  uvs1 = uvs1 + uv_part * (2._rkind)
                  uvs2 = uvs2 + uv_part * (uui+uuip)
                  uvs3 = uvs3 + uv_part * (vvi+vvip)
                  uvs4 = uvs4 + uv_part * (wwi+wwip)
                  uvs5 = uvs5 + uv_part * (enti+entip)
                  uvs6 = uvs6 + (2._rkind)*(ppi+ppip)
                enddo
                ft1 = ft1 + coeff_deriv1_ompc(l,lmax)*uvs1
                ft2 = ft2 + coeff_deriv1_ompc(l,lmax)*uvs2
                ft3 = ft3 + coeff_deriv1_ompc(l,lmax)*uvs3
                ft4 = ft4 + coeff_deriv1_ompc(l,lmax)*uvs4
                ft5 = ft5 + coeff_deriv1_ompc(l,lmax)*uvs5
                ft6 = ft6 + coeff_deriv1_ompc(l,lmax)*uvs6
              enddo
            endif

            fh1 = 0.25_rkind*ft1
            fh2 = 0.25_rkind*ft2
            fh3 = 0.25_rkind*ft3
            fh4 = 0.25_rkind*ft4
            fh5 = 0.25_rkind*ft5

            if ((i==0 .and. force_zero_flux_min == 1).or.(i==nx .and. force_zero_flux_max == 1)) then
              fh1 = 0._rkind
              fh2 = 0._rkind
              fh3 = 0._rkind
              fh4 = 0._rkind
              fh5 = 0._rkind
            endif
            fh2 = fh2 + 0.5_rkind*ft6

            fhat_ompc(i,j,k,1) = fh1
            fhat_ompc(i,j,k,2) = fh2
            fhat_ompc(i,j,k,3) = fh3
            fhat_ompc(i,j,k,4) = fh4
            fhat_ompc(i,j,k,5) = fh5
          else
            call compute_roe_average(nx, ny, nz, ng, i, i+1, j, j, k, k, w_aux_ompc, rgas0, b1, b2,&
            & b3, c, ci, h, uu, vv, ww, cp_coeff_ompc, indx_cp_l, indx_cp_r, calorically_perfect, tol_iter_nr,&
            &t0)

            call eigenvectors_x(b1, b2, b3, uu, vv, ww, c, ci, h, el, er)

            do m=1,5
              evmax(m) = -1._rkind
            enddo
            do l=1,weno_size
              ll = i + l - weno_scheme
              uu = w_aux_ompc(ll,j,k,2)
              tt = w_aux_ompc(ll,j,k,6)
              gamloc = get_gamloc_dev(indx_cp_l,indx_cp_r,tt,t0,cp_coeff_ompc,calorically_perfect,rgas0)
              c = sqrt (gamloc*rgas0*tt)
              evmax(1) = max(abs(uu-c),evmax(1))
              evmax(2) = max(abs(uu ),evmax(2))
              evmax(3) = evmax(2)
              evmax(4) = evmax(2)
              evmax(5) = max(abs(uu+c),evmax(5))
            enddo
            do l=1,weno_size
              ll = i + l - weno_scheme

              rho = w_aux_ompc(ll,j,k,1)
              uu = w_aux_ompc(ll,j,k,2)
              vv = w_aux_ompc(ll,j,k,3)
              ww = w_aux_ompc(ll,j,k,4)
              h = w_aux_ompc(ll,j,k,5)
              rhou = rho*uu
              pp = rho*w_aux_ompc(ll,j,k,6)*rgas0
              fi(1) = rhou
              fi(2) = uu * rhou + pp
              fi(3) = vv * rhou
              fi(4) = ww * rhou
              fi(5) = h * rhou
              do m=1,5
                wc = 0._rkind
                gc = 0._rkind

                wc = wc + el(1,m) * rho
                gc = gc + el(1,m) * fi(1)
                wc = wc + el(2,m) * rho*uu
                gc = gc + el(2,m) * fi(2)
                wc = wc + el(3,m) * rho*vv
                gc = gc + el(3,m) * fi(3)
                wc = wc + el(4,m) * rho*ww
                gc = gc + el(4,m) * fi(4)
                wc = wc + el(5,m) * (rho*h-pp)
                gc = gc + el(5,m) * fi(5)

                c = 0.5_rkind * (gc + evmax(m) * wc)
                gp(m,l) = c
                gm(m,l) = gc - c
              enddo
            enddo
            wenorec_ord = max(weno_scheme+ep_ord_change_ompc(i,j,k,1),1)
            call wenorec_1d(nv,gp,gm,fi,weno_scheme,wenorec_ord,weno_version,rho0,u0)
            do m=1,5
              fhat_ompc(i,j,k,m) = 0._rkind
              do mm=1,5
                fhat_ompc(i,j,k,m) = fhat_ompc(i,j,k,m) + er(mm,m) * fi(mm)
              enddo
            enddo

          endif
        enddo
      enddo
    enddo


  endsubroutine euler_x_fluxes_hybrid_kernel


  subroutine euler_x_fluxes_hybrid_rusanov_kernel(nv,nv_aux,nx,ny,nz,ng,istart_face,iend_face,&
  &lmax_base,nkeep,rgas0,w_aux_ompc,coeff_deriv1_ompc,dcsidx_ompc,fhat_ompc,force_zero_flux_min,&
  &force_zero_flux_max,weno_scheme,weno_version,sensor_threshold,weno_size,cp_coeff_ompc,indx_cp_l,&
  &indx_cp_r,ep_ord_change_ompc,calorically_perfect,tol_iter_nr,rho0,u0,t0)
    integer, value :: nv, nx, ny, nz, ng, nv_aux
    integer, value :: istart_face, iend_face, lmax_base, nkeep, indx_cp_l, indx_cp_r, calorically_perfect
    integer, value :: weno_scheme, weno_size, weno_version
    integer, dimension(0:nx,0:ny,0:nz,1:3) :: ep_ord_change_ompc
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv_aux) :: w_aux_ompc
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv) :: fhat_ompc
    real(rkind), dimension(indx_cp_l:indx_cp_r+1) :: cp_coeff_ompc
    real(rkind), dimension(4,4) :: coeff_deriv1_ompc
    real(rkind), dimension(nx) :: dcsidx_ompc
    integer, value :: force_zero_flux_min, force_zero_flux_max
    real(rkind), value :: sensor_threshold, rgas0, tol_iter_nr,rho0,u0,t0
    integer :: i, j, k, m, l
    real(rkind) :: fh1, fh2, fh3, fh4, fh5
    real(rkind) :: rhom, uui, vvi, wwi, ppi, enti, rhoi, tti
    real(rkind) :: uuip, vvip, wwip, ppip, entip, rhoip, ttip
    real(rkind) :: ft1, ft2, ft3, ft4, ft5, ft6
    real(rkind) :: uvs1, uvs2, uvs3, uvs4, uvs6, uv_part
    real(rkind) :: uvs5
    integer :: ii, lmax, wenorec_ord
    integer :: ishk
    real(rkind) :: b1, b2, b3, c, ci, h, uu, vv, ww
    real(rkind), dimension(5) :: fi
    real(rkind), dimension(5,8) :: gp,gm
    integer :: ll, mm
    real(rkind) :: evm,evmax,rhoevm
    real(rkind) :: rho, pp, wc, gc, rhou
    real(rkind) :: tt, gamloc
    real(rkind) :: uvs5_i,uvs5_k,uvs5_p,eei,eeip
    real(rkind) :: drho, dee, eem
    real(rkind) :: drhof, deef
    real(rkind) :: sumnumrho,sumnumee,sumdenrho,sumdenee
    integer :: n,n2
    real(rkind) :: t_sumdenrho, t_sumdenee, t2_sumdenrho, t2_sumdenee

    !$omp parallel do default(firstprivate) shared(w_aux_ompc, &
    !$omp& fhat_ompc, &
    !$omp& cp_coeff_ompc, &
    !$omp& coeff_deriv1_ompc, &
    !$omp& dcsidx_ompc, &
    !$omp& ep_ord_change_ompc)
    do k = 1,nz
      do j = 1,ny
        do i = +istart_face-1+1,iend_face
          ishk = 0
          do ii=i-weno_scheme+1,i+weno_scheme
            if (w_aux_ompc(ii,j,k,8) > sensor_threshold) ishk = 1
          enddo

          if (ishk == 0) then

            ft1 = 0._rkind
            ft2 = 0._rkind
            ft3 = 0._rkind
            ft4 = 0._rkind
            ft5 = 0._rkind
            ft6 = 0._rkind
            lmax = max(lmax_base+ep_ord_change_ompc(i,j,k,1),1)
            if (nkeep>=0) then
              do l=1,lmax
                uvs1 = 0._rkind
                uvs2 = 0._rkind
                uvs3 = 0._rkind
                uvs4 = 0._rkind
                uvs5_i = 0._rkind
                uvs5_k = 0._rkind
                uvs5_p = 0._rkind
                uvs6 = 0._rkind
                do m=0,l-1

                  rhoi = w_aux_ompc(i-m,j,k,1)
                  uui = w_aux_ompc(i-m,j,k,2)
                  vvi = w_aux_ompc(i-m,j,k,3)
                  wwi = w_aux_ompc(i-m,j,k,4)
                  enti = w_aux_ompc(i-m,j,k,5)
                  tti = w_aux_ompc(i-m,j,k,6)
                  ppi = tti*rhoi*rgas0
                  eei = enti-ppi/rhoi-0.5_rkind*(uui*uui+vvi*vvi+wwi*wwi)

                  rhoip = w_aux_ompc(i-m+l,j,k,1)
                  uuip = w_aux_ompc(i-m+l,j,k,2)
                  vvip = w_aux_ompc(i-m+l,j,k,3)
                  wwip = w_aux_ompc(i-m+l,j,k,4)
                  entip = w_aux_ompc(i-m+l,j,k,5)
                  ttip = w_aux_ompc(i-m+l,j,k,6)
                  ppip = ttip*rhoip*rgas0
                  eeip = entip-ppip/rhoip-0.5_rkind*(uuip*uuip+vvip*vvip+wwip*wwip)

                  rhom = rhoi+rhoip
                  eem = eei + eeip

                  if(nkeep == 0) then
                    drhof = 1._rkind
                    deef = 1._rkind
                  else
                    sumnumrho = 1._rkind
                    drho = 2._rkind*(rhoip-rhoi)/rhom
                    dee = 2._rkind*(eeip - eei)/eem
                    t_sumdenrho = (0.5_rkind*drho)*(0.5_rkind*drho)
                    t_sumdenee = (0.5_rkind*dee )*(0.5_rkind*dee )
                    t2_sumdenrho = t_sumdenrho
                    t2_sumdenee = t_sumdenee
                    sumdenrho = 1._rkind + t_sumdenrho / (3._rkind)
                    sumdenee = 1._rkind + t_sumdenee
                    sumnumee = 1._rkind + t_sumdenee / (3._rkind)
                    do n = 2, nkeep
                      n2 = 2*n
                      t_sumdenrho = t2_sumdenrho * t_sumdenrho
                      t_sumdenee = t2_sumdenee * t_sumdenee
                      sumdenrho = sumdenrho + t_sumdenrho / (1._rkind+n2)
                      sumdenee = sumdenee + t_sumdenee
                      sumnumee = sumnumee + t_sumdenee / (1._rkind+n2)
                    enddo
                    drhof = sumnumrho/sumdenrho
                    deef = sumnumee /sumdenee
                  endif

                  uv_part = (uui+uuip) * rhom * drhof
                  uvs1 = uvs1 + uv_part * (2._rkind)
                  uvs2 = uvs2 + uv_part * (uui+uuip)
                  uvs3 = uvs3 + uv_part * (vvi+vvip)
                  uvs4 = uvs4 + uv_part * (wwi+wwip)
                  uvs5_i = uvs5_i + uv_part * eem * deef
                  uvs5_k = uvs5_k + uv_part * (uui*uuip+vvi*vvip+wwi*wwip)
                  uvs5_p = uvs5_p + 4._rkind*(uui*ppip+uuip*ppi)
                  uvs6 = uvs6 + (2._rkind)*(ppi+ppip)
                enddo
                ft1 = ft1 + coeff_deriv1_ompc(l,lmax)*uvs1
                ft2 = ft2 + coeff_deriv1_ompc(l,lmax)*uvs2
                ft3 = ft3 + coeff_deriv1_ompc(l,lmax)*uvs3
                ft4 = ft4 + coeff_deriv1_ompc(l,lmax)*uvs4
                ft5 = ft5 + coeff_deriv1_ompc(l,lmax)*(uvs5_i+uvs5_k+uvs5_p)
                ft6 = ft6 + coeff_deriv1_ompc(l,lmax)*uvs6
              enddo
            else
              do l=1,lmax
                uvs1 = 0._rkind
                uvs2 = 0._rkind
                uvs3 = 0._rkind
                uvs4 = 0._rkind
                uvs5 = 0._rkind
                uvs6 = 0._rkind
                do m=0,l-1

                  rhoi = w_aux_ompc(i-m,j,k,1)
                  uui = w_aux_ompc(i-m,j,k,2)
                  vvi = w_aux_ompc(i-m,j,k,3)
                  wwi = w_aux_ompc(i-m,j,k,4)
                  enti = w_aux_ompc(i-m,j,k,5)
                  tti = w_aux_ompc(i-m,j,k,6)
                  ppi = tti*rhoi*rgas0

                  rhoip = w_aux_ompc(i-m+l,j,k,1)
                  uuip = w_aux_ompc(i-m+l,j,k,2)
                  vvip = w_aux_ompc(i-m+l,j,k,3)
                  wwip = w_aux_ompc(i-m+l,j,k,4)
                  entip = w_aux_ompc(i-m+l,j,k,5)
                  ttip = w_aux_ompc(i-m+l,j,k,6)
                  ppip = ttip*rhoip*rgas0

                  rhom = rhoi+rhoip
                  uv_part = (uui+uuip) * rhom
                  uvs1 = uvs1 + uv_part * (2._rkind)
                  uvs2 = uvs2 + uv_part * (uui+uuip)
                  uvs3 = uvs3 + uv_part * (vvi+vvip)
                  uvs4 = uvs4 + uv_part * (wwi+wwip)
                  uvs5 = uvs5 + uv_part * (enti+entip)
                  uvs6 = uvs6 + (2._rkind)*(ppi+ppip)
                enddo
                ft1 = ft1 + coeff_deriv1_ompc(l,lmax)*uvs1
                ft2 = ft2 + coeff_deriv1_ompc(l,lmax)*uvs2
                ft3 = ft3 + coeff_deriv1_ompc(l,lmax)*uvs3
                ft4 = ft4 + coeff_deriv1_ompc(l,lmax)*uvs4
                ft5 = ft5 + coeff_deriv1_ompc(l,lmax)*uvs5
                ft6 = ft6 + coeff_deriv1_ompc(l,lmax)*uvs6
              enddo
            endif

            fh1 = 0.25_rkind*ft1
            fh2 = 0.25_rkind*ft2
            fh3 = 0.25_rkind*ft3
            fh4 = 0.25_rkind*ft4
            fh5 = 0.25_rkind*ft5

            if ((i==0 .and. force_zero_flux_min == 1).or.(i==nx .and. force_zero_flux_max == 1)) then
              fh1 = 0._rkind
              fh2 = 0._rkind
              fh3 = 0._rkind
              fh4 = 0._rkind
              fh5 = 0._rkind
            endif
            fh2 = fh2 + 0.5_rkind*ft6

            fhat_ompc(i,j,k,1) = fh1
            fhat_ompc(i,j,k,2) = fh2
            fhat_ompc(i,j,k,3) = fh3
            fhat_ompc(i,j,k,4) = fh4
            fhat_ompc(i,j,k,5) = fh5
          else
            evmax = -1._rkind
            do l=1,weno_size
              ll = i + l - weno_scheme
              uu = w_aux_ompc(ll,j,k,2)
              tt = w_aux_ompc(ll,j,k,6)
              gamloc = get_gamloc_dev(indx_cp_l,indx_cp_r,tt,t0,cp_coeff_ompc,calorically_perfect,rgas0)
              c = sqrt (gamloc*rgas0*tt)
              evm = max(abs(uu-c),abs(uu+c))
              evmax = max(evm,evmax)
            enddo
            do l=1,weno_size
              ll = i + l - weno_scheme

              rho = w_aux_ompc(ll,j,k,1)
              uu = w_aux_ompc(ll,j,k,2)
              vv = w_aux_ompc(ll,j,k,3)
              ww = w_aux_ompc(ll,j,k,4)
              h = w_aux_ompc(ll,j,k,5)
              rhou = rho*uu
              pp = rho*w_aux_ompc(ll,j,k,6)*rgas0

              rhoevm = rho*evmax
              evm = rhou
              c = 0.5_rkind * (evm + rhoevm)
              gp(1,l) = c
              gm(1,l) = evm - c
              evm = uu * rhou + pp
              c = 0.5_rkind * (evm + rhoevm * uu)
              gp(2,l) = c
              gm(2,l) = evm - c
              evm = vv * rhou
              c = 0.5_rkind * (evm + rhoevm * vv)
              gp(3,l) = c
              gm(3,l) = evm - c
              evm = ww * rhou
              c = 0.5_rkind * (evm + rhoevm * ww)
              gp(4,l) = c
              gm(4,l) = evm - c
              evm = h * rhou
              c = 0.5_rkind * (evm + evmax * (rho*h-pp))
              gp(5,l) = c
              gm(5,l) = evm - c
            enddo
            wenorec_ord = max(weno_scheme+ep_ord_change_ompc(i,j,k,1),1)
            call wenorec_1d_rusanov(nv,gp,gm,fi,weno_scheme,wenorec_ord,weno_version,rho0,u0)
            do m=1,5
              fhat_ompc(i,j,k,m) = fi(m)
            enddo

          endif
        enddo
      enddo
    enddo


  endsubroutine euler_x_fluxes_hybrid_rusanov_kernel


  subroutine euler_x_update_subroutine(nx,ny,nz,ng,nv,eul_imin,eul_imax,fhat_ompc,fl_ompc,dcsidx_ompc,stream_id)
    integer :: nx,ny,nz,ng,nv,eul_imin,eul_imax
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv), intent(in) :: fhat_ompc
    real(rkind), dimension(1:nx,1:ny,1:nz,1:nv), intent(inout) :: fl_ompc
    real(rkind), dimension(1:nx), intent(in) :: dcsidx_ompc
    integer :: stream_id
    integer :: i,j,k,m,iv,iercuda

    !$omp parallel do default(firstprivate) shared(fhat_ompc, &
    !$omp& fl_ompc, &
    !$omp& dcsidx_ompc)
    do k = 1,nz
      do j = 1,ny
        do i = eul_imin,eul_imax
          do iv=1,nv
            fl_ompc(i,j,k,iv) = fl_ompc(i,j,k,iv) + (fhat_ompc(i,j,k,iv)-fhat_ompc(i-1,j,k,iv))*dcsidx_ompc(i)
          enddo
        enddo
      enddo
    enddo


  endsubroutine euler_x_update_subroutine



  subroutine euler_x_update_c2_subroutine(nx,ny,nz,ng,nv,eul_imin,eul_imax,fhat_ompc,fl_ompc,jac_ompc,stream_id)
    integer :: nx,ny,nz,ng,nv,eul_imin,eul_imax
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv), intent(in) :: fhat_ompc
    real(rkind), dimension(1:nx,1:ny,1:nz,1:nv), intent(inout) :: fl_ompc
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng), intent(in) :: jac_ompc
    integer :: stream_id
    integer :: i,j,k,m,iv,iercuda

    !$omp parallel do default(firstprivate) shared(fhat_ompc, &
    !$omp& fl_ompc, &
    !$omp& jac_ompc)
    do k = 1,nz
      do j = 1,ny
        do i = eul_imin,eul_imax
          do iv=1,nv
            fl_ompc(i,j,k,iv) = fl_ompc(i,j,k,iv) + (fhat_ompc(i,j,k,iv)-fhat_ompc(i-1,j,k,iv))*jac_ompc(i,j)
          enddo
        enddo
      enddo
    enddo


  endsubroutine euler_x_update_c2_subroutine




  subroutine euler_y_update_subroutine(nx,ny,nz,ng,nv,eul_jmin,eul_jmax,fhat_ompc,fl_ompc,detady_ompc,stream_id)
    integer :: nx,ny,nz,ng,nv,eul_jmin,eul_jmax
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv), intent(in) :: fhat_ompc
    real(rkind), dimension(1:nx,1:ny,1:nz,1:nv), intent(inout) :: fl_ompc
    real(rkind), dimension(1:ny), intent(in) :: detady_ompc
    integer :: stream_id
    integer :: i,j,k,m,iv,iercuda

    !$omp parallel do default(firstprivate) shared(fhat_ompc, &
    !$omp& fl_ompc, &
    !$omp& detady_ompc)
    do k = 1,nz
      do j = eul_jmin,eul_jmax
        do i = 1,nx
          do iv=1,nv
            fl_ompc(i,j,k,iv) = fl_ompc(i,j,k,iv) + (fhat_ompc(i,j,k,iv)-fhat_ompc(i,j-1,k,iv))*detady_ompc(j)
          enddo
        enddo
      enddo
    enddo


  endsubroutine euler_y_update_subroutine



  subroutine euler_y_update_c2_subroutine(nx,ny,nz,ng,nv,eul_jmin,eul_jmax,fhat_ompc,fl_ompc,jac_ompc,wall_tag_ompc,stream_id)
    integer :: nx,ny,nz,ng,nv,eul_jmin,eul_jmax
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv), intent(in) :: fhat_ompc
    real(rkind), dimension(1:nx,1:ny,1:nz,1:nv), intent(inout) :: fl_ompc
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng), intent(in) :: jac_ompc
    integer, dimension(1-ng:nx+ng), intent(in) :: wall_tag_ompc
    integer :: stream_id
    integer :: i,j,k,m,iv,iercuda

    !$omp parallel do default(firstprivate) shared(fhat_ompc, &
    !$omp& fl_ompc, &
    !$omp& jac_ompc, &
    !$omp& wall_tag_ompc)
    do k = 1,nz
      do j = 1,eul_jmax
        do i = 1,nx
          if(j > 1 .or. eul_jmin == 1 .or. wall_tag_ompc(i) > 0) then
            do iv=1,nv
              fl_ompc(i,j,k,iv) = fl_ompc(i,j,k,iv) + (fhat_ompc(i,j,k,iv)-fhat_ompc(i,j-1,k,iv))*jac_ompc(i,j)
            enddo
          endif
        enddo
      enddo
    enddo


  endsubroutine euler_y_update_c2_subroutine


  subroutine euler_z_update_subroutine(nx,ny,nz,ng,nv,eul_kmin,eul_kmax,fhat_ompc,fl_ompc,dzitdz_ompc,stream_id)
    integer :: nx,ny,nz,ng,nv,eul_kmin,eul_kmax
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv), intent(in) :: fhat_ompc
    real(rkind), dimension(1:nx,1:ny,1:nz,1:nv), intent(inout) :: fl_ompc
    real(rkind), dimension(1:nz), intent(in) :: dzitdz_ompc
    integer :: stream_id
    integer :: i,j,k,m,iv,iercuda

    !$omp parallel do default(firstprivate) shared(fhat_ompc, &
    !$omp& fl_ompc, &
    !$omp& dzitdz_ompc)
    do k = eul_kmin,eul_kmax
      do j = 1,ny
        do i = 1,nx
          do iv=1,nv
            fl_ompc(i,j,k,iv) = fl_ompc(i,j,k,iv) + (fhat_ompc(i,j,k,iv)-fhat_ompc(i,j,k-1,iv))*dzitdz_ompc(k)
          enddo
        enddo
      enddo
    enddo


  endsubroutine euler_z_update_subroutine




  subroutine euler_z_hybrid_kernel(nv,nv_aux,nx,ny,nz,ng,eul_kmin,eul_kmax,lmax_base,nkeep,rgas0,&
  &w_aux_ompc,fl_ompc,coeff_deriv1_ompc,dzitdz_ompc,fhat_ompc,force_zero_flux_min,force_zero_flux_max,&
  &weno_scheme,weno_version,sensor_threshold,weno_size,cp_coeff_ompc,indx_cp_l,indx_cp_r,&
  &ep_ord_change_ompc,calorically_perfect,tol_iter_nr,rho0,u0,t0)
    integer, value :: nv, nx, ny, nz, ng, nv_aux
    integer, value :: eul_kmin, eul_kmax, lmax_base, nkeep, indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv_aux) :: w_aux_ompc
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv) :: fhat_ompc
    integer, dimension(0:nx,0:ny,0:nz,1:3) :: ep_ord_change_ompc
    real(rkind), dimension(indx_cp_l:indx_cp_r+1) :: cp_coeff_ompc
    real(rkind), dimension(nx,ny,nz,nv) :: fl_ompc
    real(rkind), dimension(4,4) :: coeff_deriv1_ompc
    real(rkind), dimension(nz) :: dzitdz_ompc
    integer, value :: force_zero_flux_min, force_zero_flux_max
    real(rkind), value :: sensor_threshold, rgas0, tol_iter_nr, rho0, u0, t0
    integer, value :: weno_scheme, weno_size, weno_version
    integer :: i, j, k, m, l
    real(rkind) :: fh1, fh2, fh3, fh4, fh5
    real(rkind) :: rhom, uui, vvi, wwi, ppi, enti, rhoi, tti
    real(rkind) :: uuip, vvip, wwip, ppip, entip, rhoip, ttip
    real(rkind) :: ft1, ft2, ft3, ft4, ft5, ft6
    real(rkind) :: uvs1, uvs2, uvs3, uvs4, uvs6, uv_part
    real(rkind) :: uvs5
    integer :: kk, lmax, wenorec_ord
    integer :: ishk
    real(rkind) :: b1, b2, b3, c, ci, h, uu, vv, ww
    real(rkind), dimension(5,5) :: el, er
    real(rkind), dimension(5) :: evmax, fk
    real(rkind), dimension(5,8) :: gp, gm
    integer :: ll, mm
    real(rkind) :: rho, pp, wc, gc, rhow
    real(rkind) :: tt, gamloc
    real(rkind) :: uvs5_i,uvs5_k,uvs5_p,eei,eeip
    real(rkind) :: drho, dee, eem
    real(rkind) :: drhof, deef
    real(rkind) :: sumnumrho,sumnumee,sumdenrho,sumdenee
    integer :: n,n2
    real(rkind) :: t_sumdenrho, t_sumdenee, t2_sumdenrho, t2_sumdenee

    !$omp parallel do default(firstprivate) shared(w_aux_ompc, &
    !$omp& fhat_ompc, &
    !$omp& cp_coeff_ompc, &
    !$omp& fl_ompc, &
    !$omp& coeff_deriv1_ompc, &
    !$omp& dzitdz_ompc, &
    !$omp& ep_ord_change_ompc)
    do k = eul_kmin-1,eul_kmax
      do j = 1,ny
        do i = 1,nx
          ishk = 0
          do kk=k-weno_scheme+1,k+weno_scheme
            if (w_aux_ompc(i,j,kk,8) > sensor_threshold) ishk = 1
          enddo

          if (ishk == 0) then
            ft1 = 0._rkind
            ft2 = 0._rkind
            ft3 = 0._rkind
            ft4 = 0._rkind
            ft5 = 0._rkind
            ft6 = 0._rkind
            lmax = max(lmax_base+ep_ord_change_ompc(i,j,k,3),1)
            if (nkeep>=0) then
              do l=1,lmax
                uvs1 = 0._rkind
                uvs2 = 0._rkind
                uvs3 = 0._rkind
                uvs4 = 0._rkind
                uvs5_i = 0._rkind
                uvs5_k = 0._rkind
                uvs5_p = 0._rkind
                uvs6 = 0._rkind
                do m=0,l-1

                  rhoi = w_aux_ompc(i,j,k-m,1)
                  uui = w_aux_ompc(i,j,k-m,2)
                  vvi = w_aux_ompc(i,j,k-m,3)
                  wwi = w_aux_ompc(i,j,k-m,4)
                  enti = w_aux_ompc(i,j,k-m,5)
                  tti = w_aux_ompc(i,j,k-m,6)
                  ppi = tti*rhoi*rgas0
                  eei = enti-ppi/rhoi-0.5_rkind*(uui*uui+vvi*vvi+wwi*wwi)

                  rhoip = w_aux_ompc(i,j,k-m+l,1)
                  uuip = w_aux_ompc(i,j,k-m+l,2)
                  vvip = w_aux_ompc(i,j,k-m+l,3)
                  wwip = w_aux_ompc(i,j,k-m+l,4)
                  entip = w_aux_ompc(i,j,k-m+l,5)
                  ttip = w_aux_ompc(i,j,k-m+l,6)
                  ppip = ttip*rhoip*rgas0
                  eeip = entip-ppip/rhoip-0.5_rkind*(uuip*uuip+vvip*vvip+wwip*wwip)

                  rhom = rhoi + rhoip
                  eem = eei + eeip

                  if(nkeep == 0) then
                    drhof = 1._rkind
                    deef = 1._rkind
                  else
                    sumnumrho = 1._rkind
                    drho = 2._rkind*(rhoip-rhoi)/rhom
                    dee = 2._rkind*(eeip - eei)/eem
                    t_sumdenrho = (0.5_rkind*drho)*(0.5_rkind*drho)
                    t_sumdenee = (0.5_rkind*dee )*(0.5_rkind*dee )
                    t2_sumdenrho = t_sumdenrho
                    t2_sumdenee = t_sumdenee
                    sumdenrho = 1._rkind + t_sumdenrho / (3._rkind)
                    sumdenee = 1._rkind + t_sumdenee
                    sumnumee = 1._rkind + t_sumdenee / (3._rkind)
                    do n = 2, nkeep
                      n2 = 2*n
                      t_sumdenrho = t2_sumdenrho * t_sumdenrho
                      t_sumdenee = t2_sumdenee * t_sumdenee
                      sumdenrho = sumdenrho + t_sumdenrho / (1._rkind+n2)
                      sumdenee = sumdenee + t_sumdenee
                      sumnumee = sumnumee + t_sumdenee / (1._rkind+n2)
                    enddo
                    drhof = sumnumrho/sumdenrho
                    deef = sumnumee /sumdenee
                  endif

                  uv_part = (wwi+wwip) * rhom * drhof
                  uvs1 = uvs1 + uv_part * (2._rkind)
                  uvs2 = uvs2 + uv_part * (uui+uuip)
                  uvs3 = uvs3 + uv_part * (vvi+vvip)
                  uvs4 = uvs4 + uv_part * (wwi+wwip)
                  uvs5_i = uvs5_i + uv_part * eem * deef
                  uvs5_k = uvs5_k + uv_part * (uui*uuip+vvi*vvip+wwi*wwip)
                  uvs5_p = uvs5_p + 4._rkind*(wwi*ppip+wwip*ppi)
                  uvs6 = uvs6 + (2._rkind)*(ppi+ppip)
                enddo
                ft1 = ft1 + coeff_deriv1_ompc(l,lmax)*uvs1
                ft2 = ft2 + coeff_deriv1_ompc(l,lmax)*uvs2
                ft3 = ft3 + coeff_deriv1_ompc(l,lmax)*uvs3
                ft4 = ft4 + coeff_deriv1_ompc(l,lmax)*uvs4
                ft5 = ft5 + coeff_deriv1_ompc(l,lmax)*(uvs5_i+uvs5_k+uvs5_p)
                ft6 = ft6 + coeff_deriv1_ompc(l,lmax)*uvs6
              enddo
            else
              do l=1,lmax
                uvs1 = 0._rkind
                uvs2 = 0._rkind
                uvs3 = 0._rkind
                uvs4 = 0._rkind
                uvs5 = 0._rkind
                uvs6 = 0._rkind
                do m=0,l-1

                  rhoi = w_aux_ompc(i,j,k-m,1)
                  uui = w_aux_ompc(i,j,k-m,2)
                  vvi = w_aux_ompc(i,j,k-m,3)
                  wwi = w_aux_ompc(i,j,k-m,4)
                  enti = w_aux_ompc(i,j,k-m,5)
                  tti = w_aux_ompc(i,j,k-m,6)
                  ppi = tti*rhoi*rgas0

                  rhoip = w_aux_ompc(i,j,k-m+l,1)
                  uuip = w_aux_ompc(i,j,k-m+l,2)
                  vvip = w_aux_ompc(i,j,k-m+l,3)
                  wwip = w_aux_ompc(i,j,k-m+l,4)
                  entip = w_aux_ompc(i,j,k-m+l,5)
                  ttip = w_aux_ompc(i,j,k-m+l,6)
                  ppip = ttip*rhoip*rgas0

                  rhom = rhoi+rhoip
                  uv_part = (wwi+wwip) * rhom
                  uvs1 = uvs1 + uv_part * (2._rkind)
                  uvs2 = uvs2 + uv_part * (uui+uuip)
                  uvs3 = uvs3 + uv_part * (vvi+vvip)
                  uvs4 = uvs4 + uv_part * (wwi+wwip)
                  uvs5 = uvs5 + uv_part * (enti+entip)
                  uvs6 = uvs6 + (2._rkind)*(ppi+ppip)
                enddo
                ft1 = ft1 + coeff_deriv1_ompc(l,lmax)*uvs1
                ft2 = ft2 + coeff_deriv1_ompc(l,lmax)*uvs2
                ft3 = ft3 + coeff_deriv1_ompc(l,lmax)*uvs3
                ft4 = ft4 + coeff_deriv1_ompc(l,lmax)*uvs4
                ft5 = ft5 + coeff_deriv1_ompc(l,lmax)*uvs5
                ft6 = ft6 + coeff_deriv1_ompc(l,lmax)*uvs6
              enddo
            endif
            fh1 = 0.25_rkind*ft1
            fh2 = 0.25_rkind*ft2
            fh3 = 0.25_rkind*ft3
            fh4 = 0.25_rkind*ft4
            fh5 = 0.25_rkind*ft5
            if ((k==0 .and. force_zero_flux_min == 1).or.(k==nz .and. force_zero_flux_max == 1)) then
              fh1 = 0._rkind
              fh2 = 0._rkind
              fh3 = 0._rkind
              fh4 = 0._rkind
              fh5 = 0._rkind
            endif
            fh4 = fh4 + 0.5_rkind*ft6

            fhat_ompc(i,j,k,1) = fh1
            fhat_ompc(i,j,k,2) = fh2
            fhat_ompc(i,j,k,3) = fh3
            fhat_ompc(i,j,k,4) = fh4
            fhat_ompc(i,j,k,5) = fh5
          else
            call compute_roe_average(nx, ny, nz, ng, i, i, j, j, k, k+1, w_aux_ompc, rgas0, b1, b2,&
            & b3, c, ci, h, uu, vv, ww, cp_coeff_ompc, indx_cp_l, indx_cp_r, calorically_perfect, tol_iter_nr,&
            &t0)

            call eigenvectors_z(b1, b2, b3, uu, vv, ww, c, ci, h, el, er)

            do m=1,5
              evmax(m) = -1._rkind
            enddo
            do l=1,weno_size
              ll = k + l - weno_scheme
              ww = w_aux_ompc(i,j,ll,4)
              tt = w_aux_ompc(i,j,ll,6)
              gamloc = get_gamloc_dev(indx_cp_l,indx_cp_r,tt,t0,cp_coeff_ompc,calorically_perfect,rgas0)
              c = sqrt (gamloc*rgas0*tt)
              evmax(1) = max(abs(ww-c),evmax(1))
              evmax(2) = max(abs(ww ),evmax(2))
              evmax(3) = evmax(2)
              evmax(4) = evmax(2)
              evmax(5) = max(abs(ww+c),evmax(5))
            enddo
            do l=1,weno_size
              ll = k + l - weno_scheme

              rho = w_aux_ompc(i,j,ll,1)
              uu = w_aux_ompc(i,j,ll,2)
              vv = w_aux_ompc(i,j,ll,3)
              ww = w_aux_ompc(i,j,ll,4)
              h = w_aux_ompc(i,j,ll,5)
              rhow = rho*ww
              pp = rho*w_aux_ompc(i,j,ll,6)*rgas0
              fk(1) = rhow
              fk(2) = uu * rhow
              fk(3) = vv * rhow
              fk(4) = ww * rhow + pp
              fk(5) = h * rhow
              do m=1,5
                wc = 0._rkind
                gc = 0._rkind

                wc = wc + el(1,m) * rho
                gc = gc + el(1,m) * fk(1)
                wc = wc + el(2,m) * rho*uu
                gc = gc + el(2,m) * fk(2)
                wc = wc + el(3,m) * rho*vv
                gc = gc + el(3,m) * fk(3)
                wc = wc + el(4,m) * rho*ww
                gc = gc + el(4,m) * fk(4)
                wc = wc + el(5,m) * (rho*h-pp)
                gc = gc + el(5,m) * fk(5)

                c = 0.5_rkind * (gc + evmax(m) * wc)
                gp(m,l) = c
                gm(m,l) = gc - c
              enddo
            enddo
            wenorec_ord = max(weno_scheme+ep_ord_change_ompc(i,j,k,3),1)
            call wenorec_1d(nv,gp,gm,fk,weno_scheme,wenorec_ord,weno_version,rho0,u0)
            do m=1,5
              fhat_ompc(i,j,k,m) = 0._rkind
              do mm=1,5
                fhat_ompc(i,j,k,m) = fhat_ompc(i,j,k,m) + er(mm,m) * fk(mm)
              enddo
            enddo

          endif
        enddo
      enddo
    enddo


  endsubroutine euler_z_hybrid_kernel


  subroutine euler_z_hybrid_rusanov_kernel(nv,nv_aux,nx,ny,nz,ng,eul_kmin,eul_kmax,lmax_base,nkeep,&
  &rgas0,w_aux_ompc,fl_ompc,coeff_deriv1_ompc,dzitdz_ompc,fhat_ompc,force_zero_flux_min,&
  &force_zero_flux_max,weno_scheme,weno_version,sensor_threshold,weno_size,cp_coeff_ompc,indx_cp_l,&
  &indx_cp_r,ep_ord_change_ompc,calorically_perfect,tol_iter_nr,rho0,u0,t0)
    integer, value :: nv, nx, ny, nz, ng, nv_aux
    integer, value :: eul_kmin, eul_kmax, lmax_base, nkeep, indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv_aux) :: w_aux_ompc
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv) :: fhat_ompc
    integer, dimension(0:nx,0:ny,0:nz,1:3) :: ep_ord_change_ompc
    real(rkind), dimension(indx_cp_l:indx_cp_r+1) :: cp_coeff_ompc
    real(rkind), dimension(nx,ny,nz,nv) :: fl_ompc
    real(rkind), dimension(4,4) :: coeff_deriv1_ompc
    real(rkind), dimension(nz) :: dzitdz_ompc
    integer, value :: force_zero_flux_min, force_zero_flux_max
    real(rkind), value :: sensor_threshold, rgas0, tol_iter_nr, rho0, u0, t0
    integer, value :: weno_scheme, weno_size, weno_version
    integer :: i, j, k, m, l
    real(rkind) :: fh1, fh2, fh3, fh4, fh5
    real(rkind) :: rhom, uui, vvi, wwi, ppi, enti, rhoi, tti
    real(rkind) :: uuip, vvip, wwip, ppip, entip, rhoip, ttip
    real(rkind) :: ft1, ft2, ft3, ft4, ft5, ft6
    real(rkind) :: uvs1, uvs2, uvs3, uvs4, uvs6, uv_part
    real(rkind) :: uvs5
    integer :: kk, lmax, wenorec_ord
    integer :: ishk
    real(rkind) :: b1, b2, b3, c, ci, h, uu, vv, ww
    real(rkind), dimension(5) :: fk
    real(rkind), dimension(5,8) :: gp, gm
    integer :: ll, mm
    real(rkind) :: evm, evmax, rhoevm
    real(rkind) :: rho, pp, wc, gc, rhow
    real(rkind) :: tt, gamloc
    real(rkind) :: uvs5_i,uvs5_k,uvs5_p,eei,eeip
    real(rkind) :: drho, dee, eem
    real(rkind) :: drhof, deef
    real(rkind) :: sumnumrho,sumnumee,sumdenrho,sumdenee
    integer :: n,n2
    real(rkind) :: t_sumdenrho, t_sumdenee, t2_sumdenrho, t2_sumdenee

    !$omp parallel do default(firstprivate) shared(w_aux_ompc, &
    !$omp& fhat_ompc, &
    !$omp& cp_coeff_ompc, &
    !$omp& fl_ompc, &
    !$omp& coeff_deriv1_ompc, &
    !$omp& dzitdz_ompc, &
    !$omp& ep_ord_change_ompc)
    do k = eul_kmin-1,eul_kmax
      do j = 1,ny
        do i = 1,nx
          ishk = 0
          do kk=k-weno_scheme+1,k+weno_scheme
            if (w_aux_ompc(i,j,kk,8) > sensor_threshold) ishk = 1
          enddo

          if (ishk == 0) then
            ft1 = 0._rkind
            ft2 = 0._rkind
            ft3 = 0._rkind
            ft4 = 0._rkind
            ft5 = 0._rkind
            ft6 = 0._rkind
            lmax = max(lmax_base+ep_ord_change_ompc(i,j,k,3),1)
            if (nkeep>=0) then
              do l=1,lmax
                uvs1 = 0._rkind
                uvs2 = 0._rkind
                uvs3 = 0._rkind
                uvs4 = 0._rkind
                uvs5_i = 0._rkind
                uvs5_k = 0._rkind
                uvs5_p = 0._rkind
                uvs6 = 0._rkind
                do m=0,l-1

                  rhoi = w_aux_ompc(i,j,k-m,1)
                  uui = w_aux_ompc(i,j,k-m,2)
                  vvi = w_aux_ompc(i,j,k-m,3)
                  wwi = w_aux_ompc(i,j,k-m,4)
                  enti = w_aux_ompc(i,j,k-m,5)
                  tti = w_aux_ompc(i,j,k-m,6)
                  ppi = tti*rhoi*rgas0
                  eei = enti-ppi/rhoi-0.5_rkind*(uui*uui+vvi*vvi+wwi*wwi)

                  rhoip = w_aux_ompc(i,j,k-m+l,1)
                  uuip = w_aux_ompc(i,j,k-m+l,2)
                  vvip = w_aux_ompc(i,j,k-m+l,3)
                  wwip = w_aux_ompc(i,j,k-m+l,4)
                  entip = w_aux_ompc(i,j,k-m+l,5)
                  ttip = w_aux_ompc(i,j,k-m+l,6)
                  ppip = ttip*rhoip*rgas0
                  eeip = entip-ppip/rhoip-0.5_rkind*(uuip*uuip+vvip*vvip+wwip*wwip)

                  rhom = rhoi + rhoip
                  eem = eei + eeip

                  if(nkeep == 0) then
                    drhof = 1._rkind
                    deef = 1._rkind
                  else
                    sumnumrho = 1._rkind
                    drho = 2._rkind*(rhoip-rhoi)/rhom
                    dee = 2._rkind*(eeip - eei)/eem
                    t_sumdenrho = (0.5_rkind*drho)*(0.5_rkind*drho)
                    t_sumdenee = (0.5_rkind*dee )*(0.5_rkind*dee )
                    t2_sumdenrho = t_sumdenrho
                    t2_sumdenee = t_sumdenee
                    sumdenrho = 1._rkind + t_sumdenrho / (3._rkind)
                    sumdenee = 1._rkind + t_sumdenee
                    sumnumee = 1._rkind + t_sumdenee / (3._rkind)
                    do n = 2, nkeep
                      n2 = 2*n
                      t_sumdenrho = t2_sumdenrho * t_sumdenrho
                      t_sumdenee = t2_sumdenee * t_sumdenee
                      sumdenrho = sumdenrho + t_sumdenrho / (1._rkind+n2)
                      sumdenee = sumdenee + t_sumdenee
                      sumnumee = sumnumee + t_sumdenee / (1._rkind+n2)
                    enddo
                    drhof = sumnumrho/sumdenrho
                    deef = sumnumee /sumdenee
                  endif

                  uv_part = (wwi+wwip) * rhom * drhof
                  uvs1 = uvs1 + uv_part * (2._rkind)
                  uvs2 = uvs2 + uv_part * (uui+uuip)
                  uvs3 = uvs3 + uv_part * (vvi+vvip)
                  uvs4 = uvs4 + uv_part * (wwi+wwip)
                  uvs5_i = uvs5_i + uv_part * eem * deef
                  uvs5_k = uvs5_k + uv_part * (uui*uuip+vvi*vvip+wwi*wwip)
                  uvs5_p = uvs5_p + 4._rkind*(wwi*ppip+wwip*ppi)
                  uvs6 = uvs6 + (2._rkind)*(ppi+ppip)
                enddo
                ft1 = ft1 + coeff_deriv1_ompc(l,lmax)*uvs1
                ft2 = ft2 + coeff_deriv1_ompc(l,lmax)*uvs2
                ft3 = ft3 + coeff_deriv1_ompc(l,lmax)*uvs3
                ft4 = ft4 + coeff_deriv1_ompc(l,lmax)*uvs4
                ft5 = ft5 + coeff_deriv1_ompc(l,lmax)*(uvs5_i+uvs5_k+uvs5_p)
                ft6 = ft6 + coeff_deriv1_ompc(l,lmax)*uvs6
              enddo
            else
              do l=1,lmax
                uvs1 = 0._rkind
                uvs2 = 0._rkind
                uvs3 = 0._rkind
                uvs4 = 0._rkind
                uvs5 = 0._rkind
                uvs6 = 0._rkind
                do m=0,l-1

                  rhoi = w_aux_ompc(i,j,k-m,1)
                  uui = w_aux_ompc(i,j,k-m,2)
                  vvi = w_aux_ompc(i,j,k-m,3)
                  wwi = w_aux_ompc(i,j,k-m,4)
                  enti = w_aux_ompc(i,j,k-m,5)
                  tti = w_aux_ompc(i,j,k-m,6)
                  ppi = tti*rhoi*rgas0

                  rhoip = w_aux_ompc(i,j,k-m+l,1)
                  uuip = w_aux_ompc(i,j,k-m+l,2)
                  vvip = w_aux_ompc(i,j,k-m+l,3)
                  wwip = w_aux_ompc(i,j,k-m+l,4)
                  entip = w_aux_ompc(i,j,k-m+l,5)
                  ttip = w_aux_ompc(i,j,k-m+l,6)
                  ppip = ttip*rhoip*rgas0

                  rhom = rhoi+rhoip
                  uv_part = (wwi+wwip) * rhom
                  uvs1 = uvs1 + uv_part * (2._rkind)
                  uvs2 = uvs2 + uv_part * (uui+uuip)
                  uvs3 = uvs3 + uv_part * (vvi+vvip)
                  uvs4 = uvs4 + uv_part * (wwi+wwip)
                  uvs5 = uvs5 + uv_part * (enti+entip)
                  uvs6 = uvs6 + (2._rkind)*(ppi+ppip)
                enddo
                ft1 = ft1 + coeff_deriv1_ompc(l,lmax)*uvs1
                ft2 = ft2 + coeff_deriv1_ompc(l,lmax)*uvs2
                ft3 = ft3 + coeff_deriv1_ompc(l,lmax)*uvs3
                ft4 = ft4 + coeff_deriv1_ompc(l,lmax)*uvs4
                ft5 = ft5 + coeff_deriv1_ompc(l,lmax)*uvs5
                ft6 = ft6 + coeff_deriv1_ompc(l,lmax)*uvs6
              enddo
            endif
            fh1 = 0.25_rkind*ft1
            fh2 = 0.25_rkind*ft2
            fh3 = 0.25_rkind*ft3
            fh4 = 0.25_rkind*ft4
            fh5 = 0.25_rkind*ft5
            if ((k==0 .and. force_zero_flux_min == 1).or.(k==nz .and. force_zero_flux_max == 1)) then
              fh1 = 0._rkind
              fh2 = 0._rkind
              fh3 = 0._rkind
              fh4 = 0._rkind
              fh5 = 0._rkind
            endif
            fh4 = fh4 + 0.5_rkind*ft6

            fhat_ompc(i,j,k,1) = fh1
            fhat_ompc(i,j,k,2) = fh2
            fhat_ompc(i,j,k,3) = fh3
            fhat_ompc(i,j,k,4) = fh4
            fhat_ompc(i,j,k,5) = fh5
          else
            evmax = -1._rkind
            do l=1,weno_size
              ll = k + l - weno_scheme
              ww = w_aux_ompc(i,j,ll,4)
              tt = w_aux_ompc(i,j,ll,6)
              gamloc = get_gamloc_dev(indx_cp_l,indx_cp_r,tt,t0,cp_coeff_ompc,calorically_perfect,rgas0)
              c = sqrt (gamloc*rgas0*tt)
              evm = max(abs(ww-c),abs(ww+c))
              evmax = max(evm,evmax)
            enddo
            do l=1,weno_size
              ll = k + l - weno_scheme

              rho = w_aux_ompc(i,j,ll,1)
              uu = w_aux_ompc(i,j,ll,2)
              vv = w_aux_ompc(i,j,ll,3)
              ww = w_aux_ompc(i,j,ll,4)
              h = w_aux_ompc(i,j,ll,5)
              rhow = rho*ww
              pp = rho*w_aux_ompc(i,j,ll,6)*rgas0

              rhoevm = rho*evmax
              evm = rhow
              c = 0.5_rkind * (evm + rhoevm)
              gp(1,l) = c
              gm(1,l) = evm-c
              evm = uu * rhow
              c = 0.5_rkind * (evm + rhoevm * uu)
              gp(2,l) = c
              gm(2,l) = evm-c
              evm = vv * rhow
              c = 0.5_rkind * (evm + rhoevm * vv)
              gp(3,l) = c
              gm(3,l) = evm-c
              evm = ww * rhow + pp
              c = 0.5_rkind * (evm + rhoevm * ww)
              gp(4,l) = c
              gm(4,l) = evm-c
              evm = h * rhow
              c = 0.5_rkind * (evm + evmax * (rho*h-pp))
              gp(5,l) = c
              gm(5,l) = evm-c
            enddo
            wenorec_ord = max(weno_scheme+ep_ord_change_ompc(i,j,k,3),1)
            call wenorec_1d_rusanov(nv,gp,gm,fk,weno_scheme,wenorec_ord,weno_version,rho0,u0)
            do m=1,5
              fhat_ompc(i,j,k,m) = fk(m)
            enddo

          endif
        enddo
      enddo
    enddo


  endsubroutine euler_z_hybrid_rusanov_kernel


  subroutine euler_y_hybrid_c2_kernel(nv,nv_aux,nx,ny,nz,ng,eul_jmin,eul_jmax,lmax_base,nkeep,rgas0,&
  &w_aux_ompc,fl_ompc,coeff_deriv1_ompc,fhat_ompc,detadxc2_ompc,detadxnc2_ompc,detadyc2_ompc,&
  &detadync2_ompc,jac_ompc,metajac1_ompc,wall_tag_ompc,force_zero_flux_min,force_zero_flux_max,&
  &weno_scheme,weno_version,sensor_threshold,weno_size,cp_coeff_ompc,indx_cp_l,indx_cp_r,&
  &ep_ord_change_ompc,calorically_perfect,tol_iter_nr,rho0,u0,t0)
    integer, value :: nv, nx, ny, nz, ng, nv_aux
    integer, value :: eul_jmin, eul_jmax, lmax_base, nkeep, indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv_aux) :: w_aux_ompc
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv) :: fhat_ompc
    integer, dimension(1-ng:nx+ng) :: wall_tag_ompc
    integer, dimension(0:nx,0:ny,0:nz,1:3) :: ep_ord_change_ompc
    real(rkind), dimension(indx_cp_l:indx_cp_r+1) :: cp_coeff_ompc
    real(rkind), dimension(nx,ny,nz,nv) :: fl_ompc
    real(rkind), dimension(4,4) :: coeff_deriv1_ompc
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng) :: detadxc2_ompc
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng) :: detadxnc2_ompc
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng) :: detadyc2_ompc
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng) :: detadync2_ompc
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng) :: jac_ompc
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng) :: metajac1_ompc
    integer, value :: force_zero_flux_min, force_zero_flux_max
    real(rkind), value :: sensor_threshold, rgas0, tol_iter_nr,rho0,u0,t0
    integer, value :: weno_scheme, weno_size, weno_version
    integer :: i, j, k, m, l
    real(rkind) :: fh1, fh2, fh3, fh4, fh5
    real(rkind) :: rhom, uui, uti, vvi, wwi, ppi, enti, rhoi, tti, detadxi, detadyi
    real(rkind) :: uuip, utip, vvip, wwip, ppip, entip, rhoip, ttip, detadxip, detadyip
    real(rkind) :: ft1, ft2, ft3, ft4, ft5, ft6, ft7
    real(rkind) :: uvs1, uvs2, uvs3, uvs4, uvs5, uvs6, uvs7, uv_part
    integer :: jj
    integer :: ishk
    real(rkind) :: b1, b2, b3, c, ci, h, uu, vv, ww
    real(rkind), dimension(5,5) :: el, er
    real(rkind), dimension(5) :: evmax, fj
    real(rkind), dimension(5,8) :: gp, gm
    integer :: ll, mm, lmax, wenorec_ord
    real(rkind) :: rho, pp, wc, gc, rhou, rhov, rhow, rhoe, ut
    real(rkind) :: detadxav, detadyav, detamod
    real(rkind) :: tt, gamloc
    real(rkind) :: uvs5_i,uvs5_k,uvs5_p,eei,eeip
    real(rkind) :: drho, dee, eem
    real(rkind) :: drhof, deef
    real(rkind) :: sumnumrho,sumnumee,sumdenrho,sumdenee
    integer :: n,n2,weno_scheme_j
    real(rkind) :: t_sumdenrho, t_sumdenee, t2_sumdenrho, t2_sumdenee
    integer :: iii,jjj,kkk
    real(rkind), dimension(5,5) :: eler

    !$omp parallel do default(firstprivate) shared(w_aux_ompc, &
    !$omp& fhat_ompc, &
    !$omp& cp_coeff_ompc, &
    !$omp& fl_ompc, &
    !$omp& coeff_deriv1_ompc, &
    !$omp& detadxc2_ompc, &
    !$omp& detadxnc2_ompc, &
    !$omp& detadyc2_ompc, &
    !$omp& detadync2_ompc, &
    !$omp& jac_ompc, &
    !$omp& metajac1_ompc, &
    !$omp& wall_tag_ompc, &
    !$omp& ep_ord_change_ompc)
    do k = 1,nz
      do j = 0,eul_jmax
        do i = 1,nx
          lmax = lmax_base
          weno_scheme_j = weno_scheme
          if (j <= lmax) then
            if (wall_tag_ompc(i) < 1) then
              lmax = max(lmax_base+ep_ord_change_ompc(i,j,k,2),1)
              weno_scheme_j = max(weno_scheme+ep_ord_change_ompc(i,j,k,2),1)
            endif
          else
            lmax = max(lmax_base+ep_ord_change_ompc(i,j,k,2),1)
            weno_scheme_j = max(weno_scheme+ep_ord_change_ompc(i,j,k,2),1)
          endif
          weno_size = 2*weno_scheme_j

          ishk = 0
          do jj=j-weno_scheme_j+1,j+weno_scheme_j
            if (w_aux_ompc(i,jj,k,8) > sensor_threshold) ishk = 1
          enddo

          if (ishk == 0) then
            ft1 = 0._rkind
            ft2 = 0._rkind
            ft3 = 0._rkind
            ft4 = 0._rkind
            ft5 = 0._rkind
            ft6 = 0._rkind
            ft7 = 0._rkind
            do l=1,lmax
              uvs1 = 0._rkind
              uvs2 = 0._rkind
              uvs3 = 0._rkind
              uvs4 = 0._rkind
              uvs5 = 0._rkind
              uvs6 = 0._rkind
              uvs7 = 0._rkind
              do m=0,l-1

                rhoi = w_aux_ompc(i,j-m,k,1)
                uui = w_aux_ompc(i,j-m,k,2)
                vvi = w_aux_ompc(i,j-m,k,3)
                wwi = w_aux_ompc(i,j-m,k,4)
                enti = w_aux_ompc(i,j-m,k,5)
                tti = w_aux_ompc(i,j-m,k,6)
                ppi = tti*rhoi*rgas0
                uti = (uui*detadxc2_ompc(i,j-m)+vvi*detadyc2_ompc(i,j-m))/jac_ompc(i,j-m)
                detadxi = detadxc2_ompc(i,j-m)/jac_ompc(i,j-m)
                detadyi = detadyc2_ompc(i,j-m)/jac_ompc(i,j-m)

                rhoip = w_aux_ompc(i,j-m+l,k,1)
                uuip = w_aux_ompc(i,j-m+l,k,2)
                vvip = w_aux_ompc(i,j-m+l,k,3)
                wwip = w_aux_ompc(i,j-m+l,k,4)
                entip = w_aux_ompc(i,j-m+l,k,5)
                ttip = w_aux_ompc(i,j-m+l,k,6)
                ppip = ttip*rhoip*rgas0
                utip = (uuip*detadxc2_ompc(i,j-m+l)+vvip*detadyc2_ompc(i,j-m+l))/jac_ompc(i,j-m+l)
                detadxip = detadxc2_ompc(i,j-m+l)/jac_ompc(i,j-m+l)
                detadyip = detadyc2_ompc(i,j-m+l)/jac_ompc(i,j-m+l)

                rhom = rhoi + rhoip
                uv_part = (uti+utip) * rhom
                uvs1 = uvs1 + uv_part * (2._rkind)
                uvs2 = uvs2 + uv_part * (uui+uuip)
                uvs3 = uvs3 + uv_part * (vvi+vvip)
                uvs4 = uvs4 + uv_part * (wwi+wwip)
                uvs5 = uvs5 + uv_part * (enti+entip)
                uvs6 = uvs6 + (2._rkind)*(ppi+ppip)*(detadxi+detadxip)
                uvs7 = uvs7 + (2._rkind)*(ppi+ppip)*(detadyi+detadyip)
              enddo
              ft1 = ft1 + coeff_deriv1_ompc(l,lmax)*uvs1
              ft2 = ft2 + coeff_deriv1_ompc(l,lmax)*uvs2
              ft3 = ft3 + coeff_deriv1_ompc(l,lmax)*uvs3
              ft4 = ft4 + coeff_deriv1_ompc(l,lmax)*uvs4
              ft5 = ft5 + coeff_deriv1_ompc(l,lmax)*uvs5
              ft6 = ft6 + coeff_deriv1_ompc(l,lmax)*uvs6
              ft7 = ft7 + coeff_deriv1_ompc(l,lmax)*uvs7
            enddo

            fh1 = 0.25_rkind*ft1
            fh2 = 0.25_rkind*ft2
            fh3 = 0.25_rkind*ft3
            fh4 = 0.25_rkind*ft4
            fh5 = 0.25_rkind*ft5
            if ((j==0 .and. force_zero_flux_min == 1).or.(j==ny .and. force_zero_flux_max == 1)) then
              fh1 = 0._rkind
              fh2 = 0._rkind
              fh3 = 0._rkind
              fh4 = 0._rkind
              fh5 = 0._rkind
            endif

            fh2 = fh2 + 0.25_rkind*ft6
            fh3 = fh3 + 0.25_rkind*ft7

            fhat_ompc(i,j,k,1) = fh1
            fhat_ompc(i,j,k,2) = fh2
            fhat_ompc(i,j,k,3) = fh3
            fhat_ompc(i,j,k,4) = fh4
            fhat_ompc(i,j,k,5) = fh5
          else
            call compute_roe_average(nx, ny, nz, ng, i, i, j, j+1, k, k, w_aux_ompc, rgas0, b1, b2,&
            & b3, c, ci, h, uu, vv, ww, cp_coeff_ompc, indx_cp_l, indx_cp_r, calorically_perfect, tol_iter_nr,&
            &t0)
            detadxav = 0.5_rkind * (detadxnc2_ompc(i,j) + detadxnc2_ompc(i,j+1))
            detadyav = 0.5_rkind * (detadync2_ompc(i,j) + detadync2_ompc(i,j+1))
            detamod = 1._rkind / sqrt(detadxav**2+detadyav**2)
            detadxav = detadxav * detamod
            detadyav = detadyav * detamod
            ut = detadxav * uu + detadyav * vv

            call eigenvectors_y_c2(b1, b2, b3, uu, vv, ww, c, ci, h, el, er, ut, detadxav, detadyav)

            do m=1,5
              evmax(m) = -1._rkind
            enddo
            do l=1,weno_size
              ll = j + l - weno_scheme_j
              uu = w_aux_ompc(i,ll,k,2)
              vv = w_aux_ompc(i,ll,k,3)
              ut = detadxnc2_ompc(i,ll)*uu+detadync2_ompc(i,ll)*vv
              tt = w_aux_ompc(i,ll,k,6)
              gamloc = get_gamloc_dev(indx_cp_l,indx_cp_r,tt,t0,cp_coeff_ompc,calorically_perfect,rgas0)
              c = sqrt (gamloc*rgas0*tt)
              evmax(1) = max(abs(ut-c),evmax(1))
              evmax(2) = max(abs(ut ),evmax(2))
              evmax(3) = evmax(2)
              evmax(4) = evmax(2)
              evmax(5) = max(abs(ut+c),evmax(5))
            enddo
            do l=1,weno_size
              ll = j + l - weno_scheme_j

              rho = w_aux_ompc(i,ll,k,1)
              uu = w_aux_ompc(i,ll,k,2)
              vv = w_aux_ompc(i,ll,k,3)
              ut = detadxnc2_ompc(i,ll)*uu+detadync2_ompc(i,ll)*vv
              ww = w_aux_ompc(i,ll,k,4)
              h = w_aux_ompc(i,ll,k,5)
              rhov = rho*vv
              pp = rho*w_aux_ompc(i,ll,k,6)*rgas0
              fj(1) = rho * ut
              fj(2) = (rho * uu * ut + pp * detadxnc2_ompc(i,ll))
              fj(3) = (rho * vv * ut + pp * detadync2_ompc(i,ll))
              fj(4) = rho * ww * ut
              fj(5) = rho * h * ut
              do m=1,5
                wc = 0._rkind
                gc = 0._rkind

                wc = wc + el(1,m) * rho * metajac1_ompc(i,ll)
                gc = gc + el(1,m) * fj(1) * metajac1_ompc(i,ll)
                wc = wc + el(2,m) * rho*uu * metajac1_ompc(i,ll)
                gc = gc + el(2,m) * fj(2) * metajac1_ompc(i,ll)
                wc = wc + el(3,m) * rho*vv * metajac1_ompc(i,ll)
                gc = gc + el(3,m) * fj(3) * metajac1_ompc(i,ll)
                wc = wc + el(4,m) * rho*ww * metajac1_ompc(i,ll)
                gc = gc + el(4,m) * fj(4) * metajac1_ompc(i,ll)
                wc = wc + el(5,m) * (rho*h-pp) * metajac1_ompc(i,ll)
                gc = gc + el(5,m) * fj(5) * metajac1_ompc(i,ll)

                c = 0.5_rkind * (gc + evmax(m) * wc)
                gp(m,l) = c
                gm(m,l) = gc - c
              enddo
            enddo
            call wenorec_1d(nv,gp,gm,fj,weno_scheme_j,weno_scheme_j,weno_version,rho0,u0)
            do m=1,5
              fhat_ompc(i,j,k,m) = 0._rkind
              do mm=1,5
                fhat_ompc(i,j,k,m) = fhat_ompc(i,j,k,m) + er(mm,m) * fj(mm)
              enddo
            enddo

          endif
        enddo
      enddo
    enddo


  endsubroutine euler_y_hybrid_c2_kernel



  subroutine euler_y_hybrid_kernel(nv,nv_aux,nx,ny,nz,ng,eul_jmin,eul_jmax,lmax_base,nkeep,rgas0,&
  &w_aux_ompc,fl_ompc,coeff_deriv1_ompc,detady_ompc,fhat_ompc,force_zero_flux_min,force_zero_flux_max,&
  &weno_scheme,weno_version,sensor_threshold,weno_size,cp_coeff_ompc,indx_cp_l,indx_cp_r,&
  &ep_ord_change_ompc,calorically_perfect,tol_iter_nr,rho0,u0,t0)
    integer, value :: nv, nx, ny, nz, ng, nv_aux
    integer, value :: eul_jmin, eul_jmax, lmax_base, nkeep, indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv_aux) :: w_aux_ompc
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv) :: fhat_ompc
    integer, dimension(0:nx,0:ny,0:nz,1:3) :: ep_ord_change_ompc
    real(rkind), dimension(indx_cp_l:indx_cp_r+1) :: cp_coeff_ompc
    real(rkind), dimension(nx,ny,nz,nv) :: fl_ompc
    real(rkind), dimension(4,4) :: coeff_deriv1_ompc
    real(rkind), dimension(ny) :: detady_ompc
    integer, value :: force_zero_flux_min, force_zero_flux_max
    real(rkind), value :: sensor_threshold, rgas0, tol_iter_nr,rho0,u0,t0
    integer, value :: weno_scheme, weno_size, weno_version
    integer :: i, j, k, m, l
    real(rkind) :: fh1, fh2, fh3, fh4, fh5
    real(rkind) :: rhom, uui, vvi, wwi, ppi, enti, rhoi, tti
    real(rkind) :: uuip, vvip, wwip, ppip, entip, rhoip, ttip
    real(rkind) :: ft1, ft2, ft3, ft4, ft5, ft6
    real(rkind) :: uvs1, uvs2, uvs3, uvs4, uvs6, uv_part
    real(rkind) :: uvs5
    integer :: jj
    integer :: ishk
    real(rkind) :: b1, b2, b3, c, ci, h, uu, vv, ww
    real(rkind), dimension(5,5) :: el, er
    real(rkind), dimension(5) :: evmax, fj
    real(rkind), dimension(5,8) :: gp, gm
    integer :: ll, mm, lmax, wenorec_ord
    real(rkind) :: rho, pp, wc, gc, rhov
    real(rkind) :: tt, gamloc
    real(rkind) :: uvs5_i,uvs5_k,uvs5_p,eei,eeip
    real(rkind) :: drho, dee, eem
    real(rkind) :: drhof, deef
    real(rkind) :: sumnumrho,sumnumee,sumdenrho,sumdenee
    integer :: n,n2
    real(rkind) :: t_sumdenrho, t_sumdenee, t2_sumdenrho, t2_sumdenee

    !$omp parallel do default(firstprivate) shared(w_aux_ompc, &
    !$omp& fhat_ompc, &
    !$omp& cp_coeff_ompc, &
    !$omp& fl_ompc, &
    !$omp& coeff_deriv1_ompc, &
    !$omp& detady_ompc, &
    !$omp& ep_ord_change_ompc)
    do k = 1,nz
      do j = eul_jmin-1,eul_jmax
        do i = 1,nx
          ishk = 0
          do jj=j-weno_scheme+1,j+weno_scheme
            if (w_aux_ompc(i,jj,k,8) > sensor_threshold) ishk = 1
          enddo

          if (ishk == 0) then
            ft1 = 0._rkind
            ft2 = 0._rkind
            ft3 = 0._rkind
            ft4 = 0._rkind
            ft5 = 0._rkind
            ft6 = 0._rkind
            lmax = max(lmax_base+ep_ord_change_ompc(i,j,k,2),1)
            if (nkeep>=0) then
              do l=1,lmax
                uvs1 = 0._rkind
                uvs2 = 0._rkind
                uvs3 = 0._rkind
                uvs4 = 0._rkind
                uvs5_i = 0._rkind
                uvs5_k = 0._rkind
                uvs5_p = 0._rkind
                uvs6 = 0._rkind
                do m=0,l-1

                  rhoi = w_aux_ompc(i,j-m,k,1)
                  uui = w_aux_ompc(i,j-m,k,2)
                  vvi = w_aux_ompc(i,j-m,k,3)
                  wwi = w_aux_ompc(i,j-m,k,4)
                  enti = w_aux_ompc(i,j-m,k,5)
                  tti = w_aux_ompc(i,j-m,k,6)
                  ppi = tti*rhoi*rgas0
                  eei = enti-ppi/rhoi-0.5_rkind*(uui*uui+vvi*vvi+wwi*wwi)

                  rhoip = w_aux_ompc(i,j-m+l,k,1)
                  uuip = w_aux_ompc(i,j-m+l,k,2)
                  vvip = w_aux_ompc(i,j-m+l,k,3)
                  wwip = w_aux_ompc(i,j-m+l,k,4)
                  entip = w_aux_ompc(i,j-m+l,k,5)
                  ttip = w_aux_ompc(i,j-m+l,k,6)
                  ppip = ttip*rhoip*rgas0
                  eeip = entip-ppip/rhoip-0.5_rkind*(uuip*uuip+vvip*vvip+wwip*wwip)

                  rhom = rhoi + rhoip
                  eem = eei + eeip

                  if(nkeep == 0) then
                    drhof = 1._rkind
                    deef = 1._rkind
                  else
                    sumnumrho = 1._rkind
                    drho = 2._rkind*(rhoip-rhoi)/rhom
                    dee = 2._rkind*(eeip - eei)/eem
                    t_sumdenrho = (0.5_rkind*drho)*(0.5_rkind*drho)
                    t_sumdenee = (0.5_rkind*dee )*(0.5_rkind*dee )
                    t2_sumdenrho = t_sumdenrho
                    t2_sumdenee = t_sumdenee
                    sumdenrho = 1._rkind + t_sumdenrho / (3._rkind)
                    sumdenee = 1._rkind + t_sumdenee
                    sumnumee = 1._rkind + t_sumdenee / (3._rkind)
                    do n = 2, nkeep
                      n2 = 2*n
                      t_sumdenrho = t2_sumdenrho * t_sumdenrho
                      t_sumdenee = t2_sumdenee * t_sumdenee
                      sumdenrho = sumdenrho + t_sumdenrho / (1._rkind+n2)
                      sumdenee = sumdenee + t_sumdenee
                      sumnumee = sumnumee + t_sumdenee / (1._rkind+n2)
                    enddo
                    drhof = sumnumrho/sumdenrho
                    deef = sumnumee /sumdenee
                  endif

                  uv_part = (vvi+vvip) * rhom * drhof
                  uvs1 = uvs1 + uv_part * (2._rkind)
                  uvs2 = uvs2 + uv_part * (uui+uuip)
                  uvs3 = uvs3 + uv_part * (vvi+vvip)
                  uvs4 = uvs4 + uv_part * (wwi+wwip)
                  uvs5_i = uvs5_i + uv_part * eem * deef
                  uvs5_k = uvs5_k + uv_part * (uui*uuip+vvi*vvip+wwi*wwip)
                  uvs5_p = uvs5_p + 4._rkind*(vvi*ppip+vvip*ppi)
                  uvs6 = uvs6 + (2._rkind)*(ppi+ppip)
                enddo
                ft1 = ft1 + coeff_deriv1_ompc(l,lmax)*uvs1
                ft2 = ft2 + coeff_deriv1_ompc(l,lmax)*uvs2
                ft3 = ft3 + coeff_deriv1_ompc(l,lmax)*uvs3
                ft4 = ft4 + coeff_deriv1_ompc(l,lmax)*uvs4
                ft5 = ft5 + coeff_deriv1_ompc(l,lmax)*(uvs5_i+uvs5_k+uvs5_p)
                ft6 = ft6 + coeff_deriv1_ompc(l,lmax)*uvs6
              enddo
            else
              do l=1,lmax
                uvs1 = 0._rkind
                uvs2 = 0._rkind
                uvs3 = 0._rkind
                uvs4 = 0._rkind
                uvs5 = 0._rkind
                uvs6 = 0._rkind
                do m=0,l-1

                  rhoi = w_aux_ompc(i,j-m,k,1)
                  uui = w_aux_ompc(i,j-m,k,2)
                  vvi = w_aux_ompc(i,j-m,k,3)
                  wwi = w_aux_ompc(i,j-m,k,4)
                  enti = w_aux_ompc(i,j-m,k,5)
                  tti = w_aux_ompc(i,j-m,k,6)
                  ppi = tti*rhoi*rgas0

                  rhoip = w_aux_ompc(i,j-m+l,k,1)
                  uuip = w_aux_ompc(i,j-m+l,k,2)
                  vvip = w_aux_ompc(i,j-m+l,k,3)
                  wwip = w_aux_ompc(i,j-m+l,k,4)
                  entip = w_aux_ompc(i,j-m+l,k,5)
                  ttip = w_aux_ompc(i,j-m+l,k,6)
                  ppip = ttip*rhoip*rgas0

                  rhom = rhoi + rhoip
                  uv_part = (vvi+vvip) * rhom
                  uvs1 = uvs1 + uv_part * (2._rkind)
                  uvs2 = uvs2 + uv_part * (uui+uuip)
                  uvs3 = uvs3 + uv_part * (vvi+vvip)
                  uvs4 = uvs4 + uv_part * (wwi+wwip)
                  uvs5 = uvs5 + uv_part * (enti+entip)
                  uvs6 = uvs6 + (2._rkind)*(ppi+ppip)
                enddo
                ft1 = ft1 + coeff_deriv1_ompc(l,lmax)*uvs1
                ft2 = ft2 + coeff_deriv1_ompc(l,lmax)*uvs2
                ft3 = ft3 + coeff_deriv1_ompc(l,lmax)*uvs3
                ft4 = ft4 + coeff_deriv1_ompc(l,lmax)*uvs4
                ft5 = ft5 + coeff_deriv1_ompc(l,lmax)*uvs5
                ft6 = ft6 + coeff_deriv1_ompc(l,lmax)*uvs6
              enddo
            endif
            fh1 = 0.25_rkind*ft1
            fh2 = 0.25_rkind*ft2
            fh3 = 0.25_rkind*ft3
            fh4 = 0.25_rkind*ft4
            fh5 = 0.25_rkind*ft5
            if ((j==0 .and. force_zero_flux_min == 1).or.(j==ny .and. force_zero_flux_max == 1)) then
              fh1 = 0._rkind
              fh2 = 0._rkind
              fh3 = 0._rkind
              fh4 = 0._rkind
              fh5 = 0._rkind
            endif

            fh3 = fh3 + 0.5_rkind*ft6

            fhat_ompc(i,j,k,1) = fh1
            fhat_ompc(i,j,k,2) = fh2
            fhat_ompc(i,j,k,3) = fh3
            fhat_ompc(i,j,k,4) = fh4
            fhat_ompc(i,j,k,5) = fh5
          else
            call compute_roe_average(nx, ny, nz, ng, i, i, j, j+1, k, k, w_aux_ompc, rgas0, b1, b2,&
            & b3, c, ci, h, uu, vv, ww, cp_coeff_ompc, indx_cp_l, indx_cp_r, calorically_perfect, tol_iter_nr,&
            &t0)

            call eigenvectors_y(b1, b2, b3, uu, vv, ww, c, ci, h, el, er)

            do m=1,5
              evmax(m) = -1._rkind
            enddo
            do l=1,weno_size
              ll = j + l - weno_scheme
              vv = w_aux_ompc(i,ll,k,3)
              tt = w_aux_ompc(i,ll,k,6)
              gamloc = get_gamloc_dev(indx_cp_l,indx_cp_r,tt,t0,cp_coeff_ompc,calorically_perfect,rgas0)
              c = sqrt (gamloc*rgas0*tt)
              evmax(1) = max(abs(vv-c),evmax(1))
              evmax(2) = max(abs(vv ),evmax(2))
              evmax(3) = evmax(2)
              evmax(4) = evmax(2)
              evmax(5) = max(abs(vv+c),evmax(5))
            enddo
            do l=1,weno_size
              ll = j + l - weno_scheme

              rho = w_aux_ompc(i,ll,k,1)
              uu = w_aux_ompc(i,ll,k,2)
              vv = w_aux_ompc(i,ll,k,3)
              ww = w_aux_ompc(i,ll,k,4)
              h = w_aux_ompc(i,ll,k,5)
              rhov = rho*vv
              pp = rho*w_aux_ompc(i,ll,k,6)*rgas0
              fj(1) = rhov
              fj(2) = uu * rhov
              fj(3) = vv * rhov + pp
              fj(4) = ww * rhov
              fj(5) = h * rhov
              do m=1,5
                wc = 0._rkind
                gc = 0._rkind

                wc = wc + el(1,m) * rho
                gc = gc + el(1,m) * fj(1)
                wc = wc + el(2,m) * rho*uu
                gc = gc + el(2,m) * fj(2)
                wc = wc + el(3,m) * rho*vv
                gc = gc + el(3,m) * fj(3)
                wc = wc + el(4,m) * rho*ww
                gc = gc + el(4,m) * fj(4)
                wc = wc + el(5,m) * (rho*h-pp)
                gc = gc + el(5,m) * fj(5)

                c = 0.5_rkind * (gc + evmax(m) * wc)
                gp(m,l) = c
                gm(m,l) = gc - c
              enddo
            enddo
            wenorec_ord = max(weno_scheme+ep_ord_change_ompc(i,j,k,2),1)
            call wenorec_1d(nv,gp,gm,fj,weno_scheme,wenorec_ord,weno_version,rho0,u0)
            do m=1,5
              fhat_ompc(i,j,k,m) = 0._rkind
              do mm=1,5
                fhat_ompc(i,j,k,m) = fhat_ompc(i,j,k,m) + er(mm,m) * fj(mm)
              enddo
            enddo

          endif
        enddo
      enddo
    enddo


  endsubroutine euler_y_hybrid_kernel


  subroutine euler_y_hybrid_rusanov_kernel(nv,nv_aux,nx,ny,nz,ng,eul_jmin,eul_jmax,lmax_base,nkeep,&
  &rgas0,w_aux_ompc,fl_ompc,coeff_deriv1_ompc,detady_ompc,fhat_ompc,force_zero_flux_min,&
  &force_zero_flux_max,weno_scheme,weno_version,sensor_threshold,weno_size,cp_coeff_ompc,indx_cp_l,&
  &indx_cp_r,ep_ord_change_ompc,calorically_perfect,tol_iter_nr,rho0,u0,t0)
    integer, value :: nv, nx, ny, nz, ng, nv_aux
    integer, value :: eul_jmin, eul_jmax, lmax_base, nkeep, indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv_aux) :: w_aux_ompc
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv) :: fhat_ompc
    integer, dimension(0:nx,0:ny,0:nz,1:3) :: ep_ord_change_ompc
    real(rkind), dimension(indx_cp_l:indx_cp_r+1) :: cp_coeff_ompc
    real(rkind), dimension(nx,ny,nz,nv) :: fl_ompc
    real(rkind), dimension(4,4) :: coeff_deriv1_ompc
    real(rkind), dimension(ny) :: detady_ompc
    integer, value :: force_zero_flux_min, force_zero_flux_max
    real(rkind), value :: sensor_threshold, rgas0, tol_iter_nr,rho0,u0,t0
    integer, value :: weno_scheme, weno_size, weno_version
    integer :: i, j, k, m, l
    real(rkind) :: fh1, fh2, fh3, fh4, fh5
    real(rkind) :: rhom, uui, vvi, wwi, ppi, enti, rhoi, tti
    real(rkind) :: uuip, vvip, wwip, ppip, entip, rhoip, ttip
    real(rkind) :: ft1, ft2, ft3, ft4, ft5, ft6
    real(rkind) :: uvs1, uvs2, uvs3, uvs4, uvs6, uv_part
    real(rkind) :: uvs5
    integer :: jj
    integer :: ishk
    real(rkind) :: b1, b2, b3, c, ci, h, uu, vv, ww
    real(rkind), dimension(5) :: fj
    real(rkind), dimension(5,8) :: gp, gm
    integer :: ll, mm, lmax, wenorec_ord
    real(rkind) :: evm,evmax,rhoevm
    real(rkind) :: rho, pp, wc, gc, rhov
    real(rkind) :: tt, gamloc
    real(rkind) :: uvs5_i,uvs5_k,uvs5_p,eei,eeip
    real(rkind) :: drho, dee, eem
    real(rkind) :: drhof, deef
    real(rkind) :: sumnumrho,sumnumee,sumdenrho,sumdenee
    integer :: n,n2
    real(rkind) :: t_sumdenrho, t_sumdenee, t2_sumdenrho, t2_sumdenee

    !$omp parallel do default(firstprivate) shared(w_aux_ompc, &
    !$omp& fhat_ompc, &
    !$omp& cp_coeff_ompc, &
    !$omp& fl_ompc, &
    !$omp& coeff_deriv1_ompc, &
    !$omp& detady_ompc, &
    !$omp& ep_ord_change_ompc)
    do k = 1,nz
      do j = eul_jmin-1,eul_jmax
        do i = 1,nx
          ishk = 0
          do jj=j-weno_scheme+1,j+weno_scheme
            if (w_aux_ompc(i,jj,k,8) > sensor_threshold) ishk = 1
          enddo

          if (ishk == 0) then
            ft1 = 0._rkind
            ft2 = 0._rkind
            ft3 = 0._rkind
            ft4 = 0._rkind
            ft5 = 0._rkind
            ft6 = 0._rkind
            lmax = max(lmax_base+ep_ord_change_ompc(i,j,k,2),1)
            if (nkeep>=0) then
              do l=1,lmax
                uvs1 = 0._rkind
                uvs2 = 0._rkind
                uvs3 = 0._rkind
                uvs4 = 0._rkind
                uvs5_i = 0._rkind
                uvs5_k = 0._rkind
                uvs5_p = 0._rkind
                uvs6 = 0._rkind
                do m=0,l-1

                  rhoi = w_aux_ompc(i,j-m,k,1)
                  uui = w_aux_ompc(i,j-m,k,2)
                  vvi = w_aux_ompc(i,j-m,k,3)
                  wwi = w_aux_ompc(i,j-m,k,4)
                  enti = w_aux_ompc(i,j-m,k,5)
                  tti = w_aux_ompc(i,j-m,k,6)
                  ppi = tti*rhoi*rgas0
                  eei = enti-ppi/rhoi-0.5_rkind*(uui*uui+vvi*vvi+wwi*wwi)

                  rhoip = w_aux_ompc(i,j-m+l,k,1)
                  uuip = w_aux_ompc(i,j-m+l,k,2)
                  vvip = w_aux_ompc(i,j-m+l,k,3)
                  wwip = w_aux_ompc(i,j-m+l,k,4)
                  entip = w_aux_ompc(i,j-m+l,k,5)
                  ttip = w_aux_ompc(i,j-m+l,k,6)
                  ppip = ttip*rhoip*rgas0
                  eeip = entip-ppip/rhoip-0.5_rkind*(uuip*uuip+vvip*vvip+wwip*wwip)

                  rhom = rhoi + rhoip
                  eem = eei + eeip

                  if(nkeep == 0) then
                    drhof = 1._rkind
                    deef = 1._rkind
                  else
                    sumnumrho = 1._rkind
                    drho = 2._rkind*(rhoip-rhoi)/rhom
                    dee = 2._rkind*(eeip - eei)/eem
                    t_sumdenrho = (0.5_rkind*drho)*(0.5_rkind*drho)
                    t_sumdenee = (0.5_rkind*dee )*(0.5_rkind*dee )
                    t2_sumdenrho = t_sumdenrho
                    t2_sumdenee = t_sumdenee
                    sumdenrho = 1._rkind + t_sumdenrho / (3._rkind)
                    sumdenee = 1._rkind + t_sumdenee
                    sumnumee = 1._rkind + t_sumdenee / (3._rkind)
                    do n = 2, nkeep
                      n2 = 2*n
                      t_sumdenrho = t2_sumdenrho * t_sumdenrho
                      t_sumdenee = t2_sumdenee * t_sumdenee
                      sumdenrho = sumdenrho + t_sumdenrho / (1._rkind+n2)
                      sumdenee = sumdenee + t_sumdenee
                      sumnumee = sumnumee + t_sumdenee / (1._rkind+n2)
                    enddo
                    drhof = sumnumrho/sumdenrho
                    deef = sumnumee /sumdenee
                  endif

                  uv_part = (vvi+vvip) * rhom * drhof
                  uvs1 = uvs1 + uv_part * (2._rkind)
                  uvs2 = uvs2 + uv_part * (uui+uuip)
                  uvs3 = uvs3 + uv_part * (vvi+vvip)
                  uvs4 = uvs4 + uv_part * (wwi+wwip)
                  uvs5_i = uvs5_i + uv_part * eem * deef
                  uvs5_k = uvs5_k + uv_part * (uui*uuip+vvi*vvip+wwi*wwip)
                  uvs5_p = uvs5_p + 4._rkind*(vvi*ppip+vvip*ppi)
                  uvs6 = uvs6 + (2._rkind)*(ppi+ppip)
                enddo
                ft1 = ft1 + coeff_deriv1_ompc(l,lmax)*uvs1
                ft2 = ft2 + coeff_deriv1_ompc(l,lmax)*uvs2
                ft3 = ft3 + coeff_deriv1_ompc(l,lmax)*uvs3
                ft4 = ft4 + coeff_deriv1_ompc(l,lmax)*uvs4
                ft5 = ft5 + coeff_deriv1_ompc(l,lmax)*(uvs5_i+uvs5_k+uvs5_p)
                ft6 = ft6 + coeff_deriv1_ompc(l,lmax)*uvs6
              enddo
            else
              do l=1,lmax
                uvs1 = 0._rkind
                uvs2 = 0._rkind
                uvs3 = 0._rkind
                uvs4 = 0._rkind
                uvs5 = 0._rkind
                uvs6 = 0._rkind
                do m=0,l-1

                  rhoi = w_aux_ompc(i,j-m,k,1)
                  uui = w_aux_ompc(i,j-m,k,2)
                  vvi = w_aux_ompc(i,j-m,k,3)
                  wwi = w_aux_ompc(i,j-m,k,4)
                  enti = w_aux_ompc(i,j-m,k,5)
                  tti = w_aux_ompc(i,j-m,k,6)
                  ppi = tti*rhoi*rgas0

                  rhoip = w_aux_ompc(i,j-m+l,k,1)
                  uuip = w_aux_ompc(i,j-m+l,k,2)
                  vvip = w_aux_ompc(i,j-m+l,k,3)
                  wwip = w_aux_ompc(i,j-m+l,k,4)
                  entip = w_aux_ompc(i,j-m+l,k,5)
                  ttip = w_aux_ompc(i,j-m+l,k,6)
                  ppip = ttip*rhoip*rgas0

                  rhom = rhoi + rhoip
                  uv_part = (vvi+vvip) * rhom
                  uvs1 = uvs1 + uv_part * (2._rkind)
                  uvs2 = uvs2 + uv_part * (uui+uuip)
                  uvs3 = uvs3 + uv_part * (vvi+vvip)
                  uvs4 = uvs4 + uv_part * (wwi+wwip)
                  uvs5 = uvs5 + uv_part * (enti+entip)
                  uvs6 = uvs6 + (2._rkind)*(ppi+ppip)
                enddo
                ft1 = ft1 + coeff_deriv1_ompc(l,lmax)*uvs1
                ft2 = ft2 + coeff_deriv1_ompc(l,lmax)*uvs2
                ft3 = ft3 + coeff_deriv1_ompc(l,lmax)*uvs3
                ft4 = ft4 + coeff_deriv1_ompc(l,lmax)*uvs4
                ft5 = ft5 + coeff_deriv1_ompc(l,lmax)*uvs5
                ft6 = ft6 + coeff_deriv1_ompc(l,lmax)*uvs6
              enddo
            endif
            fh1 = 0.25_rkind*ft1
            fh2 = 0.25_rkind*ft2
            fh3 = 0.25_rkind*ft3
            fh4 = 0.25_rkind*ft4
            fh5 = 0.25_rkind*ft5
            if ((j==0 .and. force_zero_flux_min == 1).or.(j==ny .and. force_zero_flux_max == 1)) then
              fh1 = 0._rkind
              fh2 = 0._rkind
              fh3 = 0._rkind
              fh4 = 0._rkind
              fh5 = 0._rkind
            endif

            fh3 = fh3 + 0.5_rkind*ft6

            fhat_ompc(i,j,k,1) = fh1
            fhat_ompc(i,j,k,2) = fh2
            fhat_ompc(i,j,k,3) = fh3
            fhat_ompc(i,j,k,4) = fh4
            fhat_ompc(i,j,k,5) = fh5
          else
            evmax = -1._rkind
            do l=1,weno_size
              ll = j + l - weno_scheme
              vv = w_aux_ompc(i,ll,k,3)
              tt = w_aux_ompc(i,ll,k,6)
              gamloc = get_gamloc_dev(indx_cp_l,indx_cp_r,tt,t0,cp_coeff_ompc,calorically_perfect,rgas0)
              c = sqrt (gamloc*rgas0*tt)
              evm = max(abs(vv-c),abs(vv+c))
              evmax = max(evm,evmax)
            enddo
            do l=1,weno_size
              ll = j + l - weno_scheme

              rho = w_aux_ompc(i,ll,k,1)
              uu = w_aux_ompc(i,ll,k,2)
              vv = w_aux_ompc(i,ll,k,3)
              ww = w_aux_ompc(i,ll,k,4)
              h = w_aux_ompc(i,ll,k,5)
              rhov = rho*vv
              pp = rho*w_aux_ompc(i,ll,k,6)*rgas0

              rhoevm = rho*evmax
              evm = rhov
              c = 0.5_rkind * (evm + rhoevm)
              gp(1,l) = c
              gm(1,l) = evm-c
              evm = uu * rhov
              c = 0.5_rkind * (evm + rhoevm * uu)
              gp(2,l) = c
              gm(2,l) = evm-c
              evm = vv * rhov + pp
              c = 0.5_rkind * (evm + rhoevm * vv)
              gp(3,l) = c
              gm(3,l) = evm-c
              evm = ww * rhov
              c = 0.5_rkind * (evm + rhoevm * ww)
              gp(4,l) = c
              gm(4,l) = evm-c
              evm = h * rhov
              c = 0.5_rkind * (evm + evmax * (rho*h-pp))
              gp(5,l) = c
              gm(5,l) = evm-c
            enddo
            wenorec_ord = max(weno_scheme+ep_ord_change_ompc(i,j,k,2),1)
            call wenorec_1d_rusanov(nv,gp,gm,fj,weno_scheme,wenorec_ord,weno_version,rho0,u0)
            do m=1,5
              fhat_ompc(i,j,k,m) = fj(m)
            enddo

          endif
        enddo
      enddo
    enddo


  endsubroutine euler_y_hybrid_rusanov_kernel


  subroutine wenorec_1d_rusanov(nvar,vp,vm,vhat,iweno,wenorec_ord,weno_version,rho0,u0)
    integer :: nvar, iweno, wenorec_ord, weno_version
    real(rkind), dimension(5,8) :: vm,vp
    real(rkind), dimension(5) :: vhat
    real(rkind) :: rho0, u0
    real(rkind), dimension(-1:4) :: dwe
    real(rkind), dimension(-1:4) :: betap,betam
    real(rkind), dimension( 5) :: betascale
    real(rkind) :: vminus, vplus
    integer :: i,l,m
    real(rkind) :: c0,c1,c2,c3,c4,d0,d1,d2,d3,summ,sump
    real(rkind) :: tau5p,tau5m,eps40
    real(rkind) :: u0_2, rho0_2u0_2, rho0_2u0_4

    u0_2 = u0*u0
    rho0_2u0_2 = rho0*rho0*u0_2
    rho0_2u0_4 = rho0_2u0_2*u0_2
    betascale(1) = 1._rkind/rho0_2u0_2
    betascale(2) = 1._rkind/rho0_2u0_4
    betascale(3) = betascale(2)
    betascale(4) = betascale(2)
    betascale(5) = betascale(2)/u0_2
    if (wenorec_ord==1) then
      i = iweno
      do m=1,5
        vminus = vp(m,i)
        vplus = vm(m,i+1)
        vhat(m) = vminus+vplus
      enddo
    elseif (wenorec_ord==2) then
      i = iweno
      dwe(1) = 2._rkind/3._rkind
      dwe(0) = 1._rkind/3._rkind
      do m=1,5
        betap(0) = (vp(m,i )-vp(m,i-1))**2
        betap(1) = (vp(m,i+1)-vp(m,i ))**2
        betap(1) = betascale(m)*betap(1)
        betap(0) = betascale(m)*betap(0)

        betam(0) = (vm(m,i+2)-vm(m,i+1))**2
        betam(1) = (vm(m,i+1)-vm(m,i ))**2
        betam(1) = betascale(m)*betam(1)
        betam(0) = betascale(m)*betam(0)
        sump = 0._rkind
        summ = 0._rkind
        do l=0,1
          betap(l) = dwe(l)/(0.000001_rkind+betap(l))**2
          betam(l) = dwe(l)/(0.000001_rkind+betam(l))**2
          sump = sump + betap(l)
          summ = summ + betam(l)
        enddo
        do l=0,1
          betap(l) = betap(l)/sump
          betam(l) = betam(l)/summ
        enddo
        vminus = betap(0) *(-vp(m,i-1)+3*vp(m,i )) + betap(1) *( vp(m,i )+ vp(m,i+1))
        vplus = betam(0) *(-vm(m,i+2)+3*vm(m,i+1)) + betam(1) *( vm(m,i )+ vm(m,i+1))
        vhat(m) = 0.5_rkind*(vminus+vplus)
      enddo
    elseif (wenorec_ord==3) then
      i = iweno
      dwe( 0) = 1._rkind/10._rkind
      dwe( 1) = 6._rkind/10._rkind
      dwe( 2) = 3._rkind/10._rkind
      d0 = 13._rkind/12._rkind
      d1 = 1._rkind/4._rkind
      c0 = 1._rkind/3._rkind
      c1 = 5._rkind/6._rkind
      c2 =-1._rkind/6._rkind
      c3 =-7._rkind/6._rkind
      c4 =11._rkind/6._rkind
      if (weno_version==0) then
        do m=1,5
          betap(2) = d0*(vp(m,i)-2._rkind*vp(m,i+1)+vp(m,i+2))**2+d1*(3._rkind*vp(m,i)-4._rkind*vp(m,i+1)+vp(m,i+2))**2
          betap(1) = d0*(vp(m,i-1)-2._rkind*vp(m,i)+vp(m,i+1))**2+d1*( vp(m,i-1)-vp(m,i+1) )**2
          betap(0) = d0*(vp(m,i)-2._rkind*vp(m,i-1)+vp(m,i-2))**2+d1*(3._rkind*vp(m,i)-4._rkind*vp(m,i-1)+vp(m,i-2))**2
          betap(2) = betascale(m)*betap(2)
          betap(1) = betascale(m)*betap(1)
          betap(0) = betascale(m)*betap(0)
          betam(2) = d0*(vm(m,i+1)-2._rkind*vm(m,i)+vm(m,i-1))**2+d1*(3._rkind*vm(m,i+1)-4._rkind*vm(m,i)+vm(m,i-1))**2
          betam(1) = d0*(vm(m,i+2)-2._rkind*vm(m,i+1)+vm(m,i))**2+d1*( vm(m,i+2)-vm(m,i) )**2
          betam(0) = d0*(vm(m,i+1)-2._rkind*vm(m,i+2)+vm(m,i+3))**2+d1*(3._rkind*vm(m,i+1)-4._rkind*vm(m,i+2)+vm(m,i+3))**2
          betam(2) = betascale(m)*betam(2)
          betam(1) = betascale(m)*betam(1)
          betam(0) = betascale(m)*betam(0)
          sump = 0._rkind
          summ = 0._rkind
          do l=0,2
            betap(l) = dwe( l)/(0.000001_rkind+betap(l))**2
            betam(l) = dwe( l)/(0.000001_rkind+betam(l))**2
            sump = sump + betap(l)
            summ = summ + betam(l)
          enddo
          do l=0,2
            betap(l) = betap(l)/sump
            betam(l) = betam(l)/summ
          enddo
          vminus = betap(2)*(c0*vp(m,i )+c1*vp(m,i+1)+c2*vp(m,i+2)) + betap(1)*(c2*vp(m,i-1)+c1*&
          &vp(m,i )+c0*vp(m,i+1)) + betap(0)*(c0*vp(m,i-2)+c3*vp(m,i-1)+c4*vp(m,i ))
          vplus = betam(2)*(c0*vm(m,i+1)+c1*vm(m,i )+c2*vm(m,i-1)) + betam(1)*(c2*vm(m,i+2)+c1*vm(m,&
          &i+1)+c0*vm(m,i )) + betam(0)*(c0*vm(m,i+3)+c3*vm(m,i+2)+c4*vm(m,i+1))
          vhat(m) = vminus+vplus
        enddo
      else
        do m=1,5
          betap(2) = d0*(vp(m,i)-2._rkind*vp(m,i+1)+vp(m,i+2))**2+d1*(3._rkind*vp(m,i)-4._rkind*vp(m,i+1)+vp(m,i+2))**2
          betap(1) = d0*(vp(m,i-1)-2._rkind*vp(m,i)+vp(m,i+1))**2+d1*( vp(m,i-1)-vp(m,i+1) )**2
          betap(0) = d0*(vp(m,i)-2._rkind*vp(m,i-1)+vp(m,i-2))**2+d1*(3._rkind*vp(m,i)-4._rkind*vp(m,i-1)+vp(m,i-2))**2
          betap(2) = betascale(m)*betap(2)
          betap(1) = betascale(m)*betap(1)
          betap(0) = betascale(m)*betap(0)
          betam(2) = d0*(vm(m,i+1)-2._rkind*vm(m,i)+vm(m,i-1))**2+d1*(3._rkind*vm(m,i+1)-4._rkind*vm(m,i)+vm(m,i-1))**2
          betam(1) = d0*(vm(m,i+2)-2._rkind*vm(m,i+1)+vm(m,i))**2+d1*( vm(m,i+2)-vm(m,i) )**2
          betam(0) = d0*(vm(m,i+1)-2._rkind*vm(m,i+2)+vm(m,i+3))**2+d1*(3._rkind*vm(m,i+1)-4._rkind*vm(m,i+2)+vm(m,i+3))**2
          betam(2) = betascale(m)*betam(2)
          betam(1) = betascale(m)*betam(1)
          betam(0) = betascale(m)*betam(0)
          eps40 = 1.e-35
          tau5p = abs(betap(0)-betap(2))+eps40
          tau5m = abs(betam(0)-betam(2))+eps40
          do l=0,2
            betap(l) = (betap(l)+eps40)/(betap(l)+tau5p)
            betam(l) = (betam(l)+eps40)/(betam(l)+tau5m)
          enddo
          sump = 0._rkind
          summ = 0._rkind
          do l=0,2
            betap(l) = dwe(l)/betap(l)
            betam(l) = dwe(l)/betam(l)
            sump = sump + betap(l)
            summ = summ + betam(l)
          enddo
          do l=0,2
            betap(l) = betap(l)/sump
            betam(l) = betam(l)/summ
          enddo
          vminus = betap(2)*(c0*vp(m,i )+c1*vp(m,i+1)+c2*vp(m,i+2)) + betap(1)*(c2*vp(m,i-1)+c1*&
          &vp(m,i )+c0*vp(m,i+1)) + betap(0)*(c0*vp(m,i-2)+c3*vp(m,i-1)+c4*vp(m,i ))
          vplus = betam(2)*(c0*vm(m,i+1)+c1*vm(m,i )+c2*vm(m,i-1)) + betam(1)*(c2*vm(m,i+2)+c1*vm(m,&
          &i+1)+c0*vm(m,i )) + betam(0)*(c0*vm(m,i+3)+c3*vm(m,i+2)+c4*vm(m,i+1))
          vhat(m) = vminus+vplus
        enddo

      endif
    elseif (wenorec_ord==4) then
      i = iweno
      dwe( 0) = 1._rkind/35._rkind
      dwe( 1) = 12._rkind/35._rkind
      dwe( 2) = 18._rkind/35._rkind
      dwe( 3) = 4._rkind/35._rkind
      d1 = 1._rkind/36._rkind
      d2 = 13._rkind/12._rkind
      d3 = 781._rkind/720._rkind
      do m=1,5
        betap(3)= d1*(-11*vp(m, i)+18*vp(m,i+1)- 9*vp(m,i+2)+ 2*vp(m,i+3))**2+d2*( 2*vp(m, i)- 5*&
        &vp(m,i+1)+ 4*vp(m,i+2)- vp(m,i+3))**2+ d3*( -vp(m, i)+ 3*vp(m,i+1)- 3*vp(m,i+2)+ vp(m,i+3))**2
        betap(2)= d1*(- 2*vp(m,i-1)- 3*vp(m,i )+ 6*vp(m,i+1)- vp(m,i+2))**2+d2*( vp(m,i-1)- 2*vp(m,&
        &i )+ vp(m,i+1) )**2+d3*( -vp(m,i-1)+ 3*vp(m,i )- 3*vp(m,i+1)+ vp(m,i+2))**2
        betap(1)= d1*( vp(m,i-2)- 6*vp(m,i-1)+ 3*vp(m,i )+ 2*vp(m,i+1))**2+d2*( vp(m,i-1)- 2*vp(m,&
        &i )+ vp(m,i+1))**2+ d3*( -vp(m,i-2)+ 3*vp(m,i-1)- 3*vp(m,i )+ vp(m,i+1))**2
        betap(0)= d1*(- 2*vp(m,i-3)+ 9*vp(m,i-2)-18*vp(m,i-1)+11*vp(m,i ))**2+d2*(- vp(m,i-3)+ 4*&
        &vp(m,i-2)- 5*vp(m,i-1)+ 2*vp(m,i ))**2+d3*( -vp(m,i-3)+ 3*vp(m,i-2)- 3*vp(m,i-1)+ vp(m,i ))**2
        betap(3) = betascale(m)*betap(3)
        betap(2) = betascale(m)*betap(2)
        betap(1) = betascale(m)*betap(1)
        betap(0) = betascale(m)*betap(0)
        betam(3)= d1*(-11*vm(m,i+1)+18*vm(m,i )- 9*vm(m,i-1)+ 2*vm(m,i-2))**2+d2*( 2*vm(m,i+1)- 5*&
        &vm(m,i )+ 4*vm(m,i-1)- vm(m,i-2))**2+d3*( -vm(m,i+1)+ 3*vm(m,i )- 3*vm(m,i-1)+ vm(m,i-2))**2
        betam(2)= d1*(- 2*vm(m,i+2)- 3*vm(m,i+1)+ 6*vm(m,i )- vm(m,i-1))**2+d2*( vm(m,i+2)- 2*vm(m,&
        &i+1)+ vm(m,i ) )**2+d3*( -vm(m,i+2)+ 3*vm(m,i+1)- 3*vm(m,i )+ vm(m,i-1))**2
        betam(1)= d1*( vm(m,i+3)- 6*vm(m,i+2)+ 3*vm(m,i+1)+ 2*vm(m,i ))**2+d2*( vm(m,i+2)- 2*vm(m,i+&
        &1)+ vm(m,i ))**2+d3*( -vm(m,i+3)+ 3*vm(m,i+2)- 3*vm(m,i+1)+ vm(m,i ))**2
        betam(0)= d1*(- 2*vm(m,i+4)+ 9*vm(m,i+3)-18*vm(m,i+2)+11*vm(m,i+1))**2+d2*(- vm(m,i+4)+ 4*&
        &vm(m,i+3)- 5*vm(m,i+2)+ 2*vm(m,i+1))**2+d3*( -vm(m,i+4)+ 3*vm(m,i+3)- 3*vm(m,i+2)+ vm(m,i+1))**2
        betam(3) = betascale(m)*betam(3)
        betam(2) = betascale(m)*betam(2)
        betam(1) = betascale(m)*betam(1)
        betam(0) = betascale(m)*betam(0)
        sump = 0._rkind
        summ = 0._rkind
        do l=0,3
          betap(l) = dwe( l)/(0.000001_rkind+betap(l))**2
          betam(l) = dwe( l)/(0.000001_rkind+betam(l))**2
          sump = sump + betap(l)
          summ = summ + betam(l)
        enddo
        do l=0,3
          betap(l) = betap(l)/sump
          betam(l) = betam(l)/summ
        enddo
        vminus = betap(3)*( 6*vp(m,i )+26*vp(m,i+1)-10*vp(m,i+2)+ 2*vp(m,i+3))+betap(2)*(-2*vp(m,i-&
        &1)+14*vp(m,i )+14*vp(m,i+1)- 2*vp(m,i+2))+betap(1)*( 2*vp(m,i-2)-10*vp(m,i-1)+26*vp(m,i )+ 6*vp(m,i+&
        &1))+betap(0)*(-6*vp(m,i-3)+26*vp(m,i-2)-46*vp(m,i-1)+50*vp(m,i ))
        vplus = betam(3)*( 6*vm(m,i+1)+26*vm(m,i )-10*vm(m,i-1)+ 2*vm(m,i-2))+betam(2)*(-2*vm(m,i+&
        &2)+14*vm(m,i+1)+14*vm(m,i )- 2*vm(m,i-1))+betam(1)*( 2*vm(m,i+3)-10*vm(m,i+2)+26*vm(m,i+1)+ 6*vm(m,&
        &i ))+betam(0)*(-6*vm(m,i+4)+26*vm(m,i+3)-46*vm(m,i+2)+50*vm(m,i+1))
        vhat(m) = (vminus+vplus)/24._rkind
      enddo
    endif
  endsubroutine wenorec_1d_rusanov


  subroutine wenorec_1d(nvar,vp,vm,vhat,iweno,wenorec_ord,weno_version,rho0,u0)
    integer :: nvar, iweno, wenorec_ord, weno_version
    real(rkind), dimension(5,8) :: vm,vp
    real(rkind), dimension(5) :: vhat
    real(rkind) :: rho0, u0
    real(rkind), dimension(-1:4) :: dwe
    real(rkind), dimension(-1:4) :: betap,betam
    real(rkind), dimension( 5) :: betascale
    real(rkind) :: vminus, vplus
    integer :: i,l,m
    real(rkind) :: c0,c1,c2,c3,c4,d0,d1,d2,d3,summ,sump
    real(rkind) :: tau5p,tau5m,eps40
    real(rkind) :: u0_2, rho0_2u0_2, rho0_2u0_4

    u0_2 = u0*u0
    rho0_2u0_2 = rho0*rho0*u0_2
    rho0_2u0_4 = rho0_2u0_2*u0_2
    betascale(1) = 1._rkind/rho0_2u0_2
    betascale(2) = betascale(1)
    betascale(3) = 1._rkind/rho0_2u0_4
    betascale(4) = betascale(3)
    betascale(5) = betascale(1)
    if (wenorec_ord==1) then
      i = iweno
      do m=1,5
        vminus = vp(m,i)
        vplus = vm(m,i+1)
        vhat(m) = vminus+vplus
      enddo
    elseif (wenorec_ord==2) then
      i = iweno
      dwe(1) = 2._rkind/3._rkind
      dwe(0) = 1._rkind/3._rkind
      do m=1,5
        betap(0) = (vp(m,i )-vp(m,i-1))**2
        betap(1) = (vp(m,i+1)-vp(m,i ))**2
        betap(1) = betascale(m)*betap(1)
        betap(0) = betascale(m)*betap(0)

        betam(0) = (vm(m,i+2)-vm(m,i+1))**2
        betam(1) = (vm(m,i+1)-vm(m,i ))**2
        betam(1) = betascale(m)*betam(1)
        betam(0) = betascale(m)*betam(0)
        sump = 0._rkind
        summ = 0._rkind
        do l=0,1
          betap(l) = dwe(l)/(0.000001_rkind+betap(l))**2
          betam(l) = dwe(l)/(0.000001_rkind+betam(l))**2
          sump = sump + betap(l)
          summ = summ + betam(l)
        enddo
        do l=0,1
          betap(l) = betap(l)/sump
          betam(l) = betam(l)/summ
        enddo
        vminus = betap(0) *(-vp(m,i-1)+3*vp(m,i )) + betap(1) *( vp(m,i )+ vp(m,i+1))
        vplus = betam(0) *(-vm(m,i+2)+3*vm(m,i+1)) + betam(1) *( vm(m,i )+ vm(m,i+1))
        vhat(m) = 0.5_rkind*(vminus+vplus)
      enddo
    elseif (wenorec_ord==3) then
      i = iweno
      dwe( 0) = 1._rkind/10._rkind
      dwe( 1) = 6._rkind/10._rkind
      dwe( 2) = 3._rkind/10._rkind
      d0 = 13._rkind/12._rkind
      d1 = 1._rkind/4._rkind
      c0 = 1._rkind/3._rkind
      c1 = 5._rkind/6._rkind
      c2 =-1._rkind/6._rkind
      c3 =-7._rkind/6._rkind
      c4 =11._rkind/6._rkind
      if (weno_version==0) then
        do m=1,5
          betap(2) = d0*(vp(m,i)-2._rkind*vp(m,i+1)+vp(m,i+2))**2+d1*(3._rkind*vp(m,i)-4._rkind*vp(m,i+1)+vp(m,i+2))**2
          betap(1) = d0*(vp(m,i-1)-2._rkind*vp(m,i)+vp(m,i+1))**2+d1*( vp(m,i-1)-vp(m,i+1) )**2
          betap(0) = d0*(vp(m,i)-2._rkind*vp(m,i-1)+vp(m,i-2))**2+d1*(3._rkind*vp(m,i)-4._rkind*vp(m,i-1)+vp(m,i-2))**2
          betap(2) = betascale(m)*betap(2)
          betap(1) = betascale(m)*betap(1)
          betap(0) = betascale(m)*betap(0)
          betam(2) = d0*(vm(m,i+1)-2._rkind*vm(m,i)+vm(m,i-1))**2+d1*(3._rkind*vm(m,i+1)-4._rkind*vm(m,i)+vm(m,i-1))**2
          betam(1) = d0*(vm(m,i+2)-2._rkind*vm(m,i+1)+vm(m,i))**2+d1*( vm(m,i+2)-vm(m,i) )**2
          betam(0) = d0*(vm(m,i+1)-2._rkind*vm(m,i+2)+vm(m,i+3))**2+d1*(3._rkind*vm(m,i+1)-4._rkind*vm(m,i+2)+vm(m,i+3))**2
          betam(2) = betascale(m)*betam(2)
          betam(1) = betascale(m)*betam(1)
          betam(0) = betascale(m)*betam(0)
          sump = 0._rkind
          summ = 0._rkind
          do l=0,2
            betap(l) = dwe( l)/(0.000001_rkind+betap(l))**2
            betam(l) = dwe( l)/(0.000001_rkind+betam(l))**2
            sump = sump + betap(l)
            summ = summ + betam(l)
          enddo
          do l=0,2
            betap(l) = betap(l)/sump
            betam(l) = betam(l)/summ
          enddo
          vminus = betap(2)*(c0*vp(m,i )+c1*vp(m,i+1)+c2*vp(m,i+2)) + betap(1)*(c2*vp(m,i-1)+c1*&
          &vp(m,i )+c0*vp(m,i+1)) + betap(0)*(c0*vp(m,i-2)+c3*vp(m,i-1)+c4*vp(m,i ))
          vplus = betam(2)*(c0*vm(m,i+1)+c1*vm(m,i )+c2*vm(m,i-1)) + betam(1)*(c2*vm(m,i+2)+c1*vm(m,&
          &i+1)+c0*vm(m,i )) + betam(0)*(c0*vm(m,i+3)+c3*vm(m,i+2)+c4*vm(m,i+1))
          vhat(m) = vminus+vplus
        enddo
      else
        do m=1,5
          betap(2) = d0*(vp(m,i)-2._rkind*vp(m,i+1)+vp(m,i+2))**2+d1*(3._rkind*vp(m,i)-4._rkind*vp(m,i+1)+vp(m,i+2))**2
          betap(1) = d0*(vp(m,i-1)-2._rkind*vp(m,i)+vp(m,i+1))**2+d1*( vp(m,i-1)-vp(m,i+1) )**2
          betap(0) = d0*(vp(m,i)-2._rkind*vp(m,i-1)+vp(m,i-2))**2+d1*(3._rkind*vp(m,i)-4._rkind*vp(m,i-1)+vp(m,i-2))**2
          betap(2) = betascale(m)*betap(2)
          betap(1) = betascale(m)*betap(1)
          betap(0) = betascale(m)*betap(0)
          betam(2) = d0*(vm(m,i+1)-2._rkind*vm(m,i)+vm(m,i-1))**2+d1*(3._rkind*vm(m,i+1)-4._rkind*vm(m,i)+vm(m,i-1))**2
          betam(1) = d0*(vm(m,i+2)-2._rkind*vm(m,i+1)+vm(m,i))**2+d1*( vm(m,i+2)-vm(m,i) )**2
          betam(0) = d0*(vm(m,i+1)-2._rkind*vm(m,i+2)+vm(m,i+3))**2+d1*(3._rkind*vm(m,i+1)-4._rkind*vm(m,i+2)+vm(m,i+3))**2
          betam(2) = betascale(m)*betam(2)
          betam(1) = betascale(m)*betam(1)
          betam(0) = betascale(m)*betam(0)
          eps40 = 1.e-35
          tau5p = abs(betap(0)-betap(2))+eps40
          tau5m = abs(betam(0)-betam(2))+eps40
          do l=0,2
            betap(l) = (betap(l)+eps40)/(betap(l)+tau5p)
            betam(l) = (betam(l)+eps40)/(betam(l)+tau5m)
          enddo
          sump = 0._rkind
          summ = 0._rkind
          do l=0,2
            betap(l) = dwe(l)/betap(l)
            betam(l) = dwe(l)/betam(l)
            sump = sump + betap(l)
            summ = summ + betam(l)
          enddo
          do l=0,2
            betap(l) = betap(l)/sump
            betam(l) = betam(l)/summ
          enddo
          vminus = betap(2)*(c0*vp(m,i )+c1*vp(m,i+1)+c2*vp(m,i+2)) + betap(1)*(c2*vp(m,i-1)+c1*&
          &vp(m,i )+c0*vp(m,i+1)) + betap(0)*(c0*vp(m,i-2)+c3*vp(m,i-1)+c4*vp(m,i ))
          vplus = betam(2)*(c0*vm(m,i+1)+c1*vm(m,i )+c2*vm(m,i-1)) + betam(1)*(c2*vm(m,i+2)+c1*vm(m,&
          &i+1)+c0*vm(m,i )) + betam(0)*(c0*vm(m,i+3)+c3*vm(m,i+2)+c4*vm(m,i+1))
          vhat(m) = vminus+vplus
        enddo

      endif
    elseif (wenorec_ord==4) then
      i = iweno
      dwe( 0) = 1._rkind/35._rkind
      dwe( 1) = 12._rkind/35._rkind
      dwe( 2) = 18._rkind/35._rkind
      dwe( 3) = 4._rkind/35._rkind
      d1 = 1._rkind/36._rkind
      d2 = 13._rkind/12._rkind
      d3 = 781._rkind/720._rkind
      do m=1,5
        betap(3)= d1*(-11*vp(m, i)+18*vp(m,i+1)- 9*vp(m,i+2)+ 2*vp(m,i+3))**2+d2*( 2*vp(m, i)- 5*&
        &vp(m,i+1)+ 4*vp(m,i+2)- vp(m,i+3))**2+ d3*( -vp(m, i)+ 3*vp(m,i+1)- 3*vp(m,i+2)+ vp(m,i+3))**2
        betap(2)= d1*(- 2*vp(m,i-1)- 3*vp(m,i )+ 6*vp(m,i+1)- vp(m,i+2))**2+d2*( vp(m,i-1)- 2*vp(m,&
        &i )+ vp(m,i+1) )**2+d3*( -vp(m,i-1)+ 3*vp(m,i )- 3*vp(m,i+1)+ vp(m,i+2))**2
        betap(1)= d1*( vp(m,i-2)- 6*vp(m,i-1)+ 3*vp(m,i )+ 2*vp(m,i+1))**2+d2*( vp(m,i-1)- 2*vp(m,&
        &i )+ vp(m,i+1))**2+ d3*( -vp(m,i-2)+ 3*vp(m,i-1)- 3*vp(m,i )+ vp(m,i+1))**2
        betap(0)= d1*(- 2*vp(m,i-3)+ 9*vp(m,i-2)-18*vp(m,i-1)+11*vp(m,i ))**2+d2*(- vp(m,i-3)+ 4*&
        &vp(m,i-2)- 5*vp(m,i-1)+ 2*vp(m,i ))**2+d3*( -vp(m,i-3)+ 3*vp(m,i-2)- 3*vp(m,i-1)+ vp(m,i ))**2
        betap(3) = betascale(m)*betap(3)
        betap(2) = betascale(m)*betap(2)
        betap(1) = betascale(m)*betap(1)
        betap(0) = betascale(m)*betap(0)
        betam(3)= d1*(-11*vm(m,i+1)+18*vm(m,i )- 9*vm(m,i-1)+ 2*vm(m,i-2))**2+d2*( 2*vm(m,i+1)- 5*&
        &vm(m,i )+ 4*vm(m,i-1)- vm(m,i-2))**2+d3*( -vm(m,i+1)+ 3*vm(m,i )- 3*vm(m,i-1)+ vm(m,i-2))**2
        betam(2)= d1*(- 2*vm(m,i+2)- 3*vm(m,i+1)+ 6*vm(m,i )- vm(m,i-1))**2+d2*( vm(m,i+2)- 2*vm(m,&
        &i+1)+ vm(m,i ) )**2+d3*( -vm(m,i+2)+ 3*vm(m,i+1)- 3*vm(m,i )+ vm(m,i-1))**2
        betam(1)= d1*( vm(m,i+3)- 6*vm(m,i+2)+ 3*vm(m,i+1)+ 2*vm(m,i ))**2+d2*( vm(m,i+2)- 2*vm(m,i+&
        &1)+ vm(m,i ))**2+d3*( -vm(m,i+3)+ 3*vm(m,i+2)- 3*vm(m,i+1)+ vm(m,i ))**2
        betam(0)= d1*(- 2*vm(m,i+4)+ 9*vm(m,i+3)-18*vm(m,i+2)+11*vm(m,i+1))**2+d2*(- vm(m,i+4)+ 4*&
        &vm(m,i+3)- 5*vm(m,i+2)+ 2*vm(m,i+1))**2+d3*( -vm(m,i+4)+ 3*vm(m,i+3)- 3*vm(m,i+2)+ vm(m,i+1))**2
        betam(3) = betascale(m)*betam(3)
        betam(2) = betascale(m)*betam(2)
        betam(1) = betascale(m)*betam(1)
        betam(0) = betascale(m)*betam(0)
        sump = 0._rkind
        summ = 0._rkind
        do l=0,3
          betap(l) = dwe( l)/(0.000001_rkind+betap(l))**2
          betam(l) = dwe( l)/(0.000001_rkind+betam(l))**2
          sump = sump + betap(l)
          summ = summ + betam(l)
        enddo
        do l=0,3
          betap(l) = betap(l)/sump
          betam(l) = betam(l)/summ
        enddo
        vminus = betap(3)*( 6*vp(m,i )+26*vp(m,i+1)-10*vp(m,i+2)+ 2*vp(m,i+3))+betap(2)*(-2*vp(m,i-&
        &1)+14*vp(m,i )+14*vp(m,i+1)- 2*vp(m,i+2))+betap(1)*( 2*vp(m,i-2)-10*vp(m,i-1)+26*vp(m,i )+ 6*vp(m,i+&
        &1))+betap(0)*(-6*vp(m,i-3)+26*vp(m,i-2)-46*vp(m,i-1)+50*vp(m,i ))
        vplus = betam(3)*( 6*vm(m,i+1)+26*vm(m,i )-10*vm(m,i-1)+ 2*vm(m,i-2))+betam(2)*(-2*vm(m,i+&
        &2)+14*vm(m,i+1)+14*vm(m,i )- 2*vm(m,i-1))+betam(1)*( 2*vm(m,i+3)-10*vm(m,i+2)+26*vm(m,i+1)+ 6*vm(m,&
        &i ))+betam(0)*(-6*vm(m,i+4)+26*vm(m,i+3)-46*vm(m,i+2)+50*vm(m,i+1))
        vhat(m) = (vminus+vplus)/24._rkind
      enddo
    endif
  endsubroutine wenorec_1d



  subroutine eigenvectors_x(b1,b2,b3,uu,vv,ww,c,ci,h,el,er)
    real(rkind), intent(in) :: b1, b2, b3, uu, vv, ww, c, ci, h
    real(rkind), dimension(:,:), intent(out) :: el, er

    el(1,1) = 0.5_rkind * (b1 + uu * ci)
    el(2,1) = -0.5_rkind * (b2 * uu + ci)
    el(3,1) = -0.5_rkind * (b2 * vv )
    el(4,1) = -0.5_rkind * (b2 * ww )
    el(5,1) = 0.5_rkind * b2
    el(1,2) = 1._rkind - b1
    el(2,2) = b2*uu
    el(3,2) = b2*vv
    el(4,2) = b2*ww
    el(5,2) = -b2
    el(1,3) = -vv
    el(2,3) = 0._rkind
    el(3,3) = 1._rkind
    el(4,3) = 0._rkind
    el(5,3) = 0._rkind
    el(1,4) = -ww
    el(2,4) = 0._rkind
    el(3,4) = 0._rkind
    el(4,4) = 1._rkind
    el(5,4) = 0._rkind
    el(1,5) = 0.5_rkind * (b1 - uu * ci)
    el(2,5) = -0.5_rkind * (b2 * uu - ci)
    el(3,5) = -0.5_rkind * (b2 * vv )
    el(4,5) = -0.5_rkind * (b2 * ww )
    el(5,5) = 0.5_rkind * b2

    er(1,1) = 1._rkind
    er(2,1) = 1._rkind
    er(3,1) = 0._rkind
    er(4,1) = 0._rkind
    er(5,1) = 1._rkind
    er(1,2) = uu - c
    er(2,2) = uu
    er(3,2) = 0._rkind
    er(4,2) = 0._rkind
    er(5,2) = uu + c
    er(1,3) = vv
    er(2,3) = vv
    er(3,3) = 1._rkind
    er(4,3) = 0._rkind
    er(5,3) = vv
    er(1,4) = ww
    er(2,4) = ww
    er(3,4) = 0._rkind
    er(4,4) = 1._rkind
    er(5,4) = ww
    er(1,5) = h - uu * c
    er(2,5) = b3
    er(3,5) = vv
    er(4,5) = ww
    er(5,5) = h + uu * c
  endsubroutine eigenvectors_x


  subroutine eigenvectors_x_c2(b1,b2,b3,uu,vv,ww,c,ci,h,el,er,ut,dcsidxav,dcsidyav)
    real(rkind), intent(in) :: b1, b2, b3, uu, vv, ww, c, ci, h, ut, dcsidxav, dcsidyav
    real(rkind), dimension(:,:), intent(out) :: el, er

    el(1,1) = 0.5_rkind * (b1 + ut * ci)
    el(2,1) = -0.5_rkind * (b2 * uu + dcsidxav * ci)
    el(3,1) = -0.5_rkind * (b2 * vv + dcsidyav * ci)
    el(4,1) = -0.5_rkind * (b2 * ww )
    el(5,1) = 0.5_rkind * b2
    el(1,2) = 1._rkind - b1
    el(2,2) = b2*uu
    el(3,2) = b2*vv
    el(4,2) = b2*ww
    el(5,2) = -b2
    el(1,3) = dcsidyav * uu - dcsidxav * vv
    el(2,3) = - dcsidyav
    el(3,3) = dcsidxav
    el(4,3) = 0._rkind
    el(5,3) = 0._rkind
    el(1,4) = -ww
    el(2,4) = 0._rkind
    el(3,4) = 0._rkind
    el(4,4) = 1._rkind
    el(5,4) = 0._rkind
    el(1,5) = 0.5_rkind * (b1 - ut * ci)
    el(2,5) = -0.5_rkind * (b2 * uu - dcsidxav * ci)
    el(3,5) = -0.5_rkind * (b2 * vv - dcsidyav * ci)
    el(4,5) = -0.5_rkind * (b2 * ww )
    el(5,5) = 0.5_rkind * b2

    er(1,1) = 1._rkind
    er(2,1) = 1._rkind
    er(3,1) = 0._rkind
    er(4,1) = 0._rkind
    er(5,1) = 1._rkind
    er(1,2) = uu - dcsidxav * c
    er(2,2) = uu
    er(3,2) = - dcsidyav
    er(4,2) = 0._rkind
    er(5,2) = uu + dcsidxav * c
    er(1,3) = vv - dcsidyav * c
    er(2,3) = vv
    er(3,3) = dcsidxav
    er(4,3) = 0._rkind
    er(5,3) = vv + dcsidyav * c
    er(1,4) = ww
    er(2,4) = ww
    er(3,4) = 0._rkind
    er(4,4) = 1._rkind
    er(5,4) = ww
    er(1,5) = h - ut * c
    er(2,5) = b3
    er(3,5) = - dcsidyav * uu + dcsidxav * vv
    er(4,5) = ww
    er(5,5) = h + ut * c
  endsubroutine eigenvectors_x_c2


  subroutine eigenvectors_y(b1,b2,b3,uu,vv,ww,c,ci,h,el,er)
    real(rkind), intent(in) :: b1, b2, b3, uu, vv, ww, c, ci, h
    real(rkind), dimension(:,:), intent(out) :: el, er

    el(1,1) = 0.5_rkind * (b1 + vv * ci)
    el(2,1) = -0.5_rkind * (b2 * uu )
    el(3,1) = -0.5_rkind * (b2 * vv + ci)
    el(4,1) = -0.5_rkind * (b2 * ww )
    el(5,1) = 0.5_rkind * b2
    el(1,2) = 1._rkind - b1
    el(2,2) = b2 * uu
    el(3,2) = b2 * vv
    el(4,2) = b2 * ww
    el(5,2) = -b2
    el(1,3) = -uu
    el(2,3) = 1._rkind
    el(3,3) = 0._rkind
    el(4,3) = 0._rkind
    el(5,3) = 0._rkind
    el(1,4) = -ww
    el(2,4) = 0._rkind
    el(3,4) = 0._rkind
    el(4,4) = 1._rkind
    el(5,4) = 0._rkind
    el(1,5) = 0.5_rkind * (b1 - vv * ci)
    el(2,5) = -0.5_rkind * (b2 * uu )
    el(3,5) = -0.5_rkind * (b2 * vv - ci)
    el(4,5) = -0.5_rkind * (b2 * ww )
    el(5,5) = 0.5_rkind * b2

    er(1,1) = 1._rkind
    er(2,1) = 1._rkind
    er(3,1) = 0._rkind
    er(4,1) = 0._rkind
    er(5,1) = 1._rkind
    er(1,2) = uu
    er(2,2) = uu
    er(3,2) = 1._rkind
    er(4,2) = 0._rkind
    er(5,2) = uu
    er(1,3) = vv - c
    er(2,3) = vv
    er(3,3) = 0._rkind
    er(4,3) = 0._rkind
    er(5,3) = vv + c
    er(1,4) = ww
    er(2,4) = ww
    er(3,4) = 0._rkind
    er(4,4) = 1._rkind
    er(5,4) = ww
    er(1,5) = h - vv * c
    er(2,5) = b3
    er(3,5) = uu
    er(4,5) = ww
    er(5,5) = h + vv * c
  endsubroutine eigenvectors_y


  subroutine eigenvectors_y_c2(b1,b2,b3,uu,vv,ww,c,ci,h,el,er,ut,detadxav,detadyav)
    real(rkind), intent(in) :: b1, b2, b3, uu, vv, ww, c, ci, h, ut, detadxav, detadyav
    real(rkind), dimension(:,:), intent(out) :: el, er

    el(1,1) = 0.5_rkind * (b1 + ut * ci)
    el(2,1) = -0.5_rkind * (b2 * uu + detadxav * ci)
    el(3,1) = -0.5_rkind * (b2 * vv + detadyav * ci)
    el(4,1) = -0.5_rkind * (b2 * ww )
    el(5,1) = 0.5_rkind * b2
    el(1,2) = 1._rkind - b1
    el(2,2) = b2 * uu
    el(3,2) = b2 * vv
    el(4,2) = b2 * ww
    el(5,2) = -b2
    el(1,3) = -detadyav * uu + detadxav * vv
    el(2,3) = +detadyav
    el(3,3) = -detadxav
    el(4,3) = 0._rkind
    el(5,3) = 0._rkind
    el(1,4) = -ww
    el(2,4) = 0._rkind
    el(3,4) = 0._rkind
    el(4,4) = 1._rkind
    el(5,4) = 0._rkind
    el(1,5) = 0.5_rkind * (b1 - ut * ci)
    el(2,5) = -0.5_rkind * (b2 * uu - detadxav * ci)
    el(3,5) = -0.5_rkind * (b2 * vv - detadyav * ci)
    el(4,5) = -0.5_rkind * (b2 * ww )
    el(5,5) = 0.5_rkind * b2

    er(1,1) = 1._rkind
    er(2,1) = 1._rkind
    er(3,1) = 0._rkind
    er(4,1) = 0._rkind
    er(5,1) = 1._rkind
    er(1,2) = uu - detadxav * c
    er(2,2) = uu
    er(3,2) = detadyav
    er(4,2) = 0._rkind
    er(5,2) = uu + detadxav * c
    er(1,3) = vv - detadyav * c
    er(2,3) = vv
    er(3,3) = - detadxav
    er(4,3) = 0._rkind
    er(5,3) = vv + detadyav * c
    er(1,4) = ww
    er(2,4) = ww
    er(3,4) = 0._rkind
    er(4,4) = 1._rkind
    er(5,4) = ww
    er(1,5) = h - ut * c
    er(2,5) = b3
    er(3,5) = detadyav * uu - detadxav * vv
    er(4,5) = ww
    er(5,5) = h + ut * c
  endsubroutine eigenvectors_y_c2


  subroutine eigenvectors_z(b1,b2,b3,uu,vv,ww,c,ci,h,el,er)
    real(rkind), intent(in) :: b1, b2, b3, uu, vv, ww, c, ci, h
    real(rkind), dimension(:,:), intent(out) :: el, er

    el(1,1) = 0.5_rkind * (b1 + ww * ci)
    el(2,1) = -0.5_rkind * (b2 * uu )
    el(3,1) = -0.5_rkind * (b2 * vv )
    el(4,1) = -0.5_rkind * (b2 * ww + ci)
    el(5,1) = 0.5_rkind * b2
    el(1,2) = 1._rkind - b1
    el(2,2) = b2 * uu
    el(3,2) = b2 * vv
    el(4,2) = b2 * ww
    el(5,2) = -b2
    el(1,3) = -uu
    el(2,3) = 1._rkind
    el(3,3) = 0._rkind
    el(4,3) = 0._rkind
    el(5,3) = 0._rkind
    el(1,4) = -vv
    el(2,4) = 0._rkind
    el(3,4) = 1._rkind
    el(4,4) = 0._rkind
    el(5,4) = 0._rkind
    el(1,5) = 0.5_rkind * (b1 - ww * ci)
    el(2,5) = -0.5_rkind * (b2 * uu )
    el(3,5) = -0.5_rkind * (b2 * vv )
    el(4,5) = -0.5_rkind * (b2 * ww - ci)
    el(5,5) = 0.5_rkind * b2

    er(1,1) = 1._rkind
    er(2,1) = 1._rkind
    er(3,1) = 0._rkind
    er(4,1) = 0._rkind
    er(5,1) = 1._rkind
    er(1,2) = uu
    er(2,2) = uu
    er(3,2) = 1._rkind
    er(4,2) = 0._rkind
    er(5,2) = uu
    er(1,3) = vv
    er(2,3) = vv
    er(3,3) = 0._rkind
    er(4,3) = 1._rkind
    er(5,3) = vv
    er(1,4) = ww - c
    er(2,4) = ww
    er(3,4) = 0._rkind
    er(4,4) = 0._rkind
    er(5,4) = ww + c
    er(1,5) = h - ww * c
    er(2,5) = b3
    er(3,5) = uu
    er(4,5) = vv
    er(5,5) = h + ww * c
  endsubroutine eigenvectors_z


  subroutine compute_roe_average(nx,ny,nz,ng,i,ip,j,jp,k,kp,w_aux_ompc,rgas0,b1,b2,b3,c,ci,h,uu,vv,&
  &ww,cp_coeff_ompc,indx_cp_l,indx_cp_r,calorically_perfect,tol_iter_nr,t0)
    integer :: ng,i,ip,j,jp,k,kp,indx_cp_l,indx_cp_r,calorically_perfect
    integer :: nx,ny,nz
    real(rkind), intent(in) :: rgas0, tol_iter_nr,t0
    real(rkind), intent(out) :: b1, b2, b3, uu, vv, ww, ci, h, c
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in) :: w_aux_ompc
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in) :: cp_coeff_ompc
    real(rkind) :: up, vp, wp, qqp, hp, r, rp1, cc, qq, gam, gm1
    integer :: ll,iter,max_iter
    real(rkind) :: tt,gm1loc,hbar,ttp,told,num,den,gamloc,cploc,tpow,tpowp,p_rho,p_e,etot,rho

    max_iter = 50
    uu = w_aux_ompc(i,j,k,2)
    vv = w_aux_ompc(i,j,k,3)
    ww = w_aux_ompc(i,j,k,4)
    qq = 0.5_rkind * (uu*uu +vv*vv + ww*ww)
    h = w_aux_ompc(i,j,k,5)
    up = w_aux_ompc(ip,jp,kp,2)
    vp = w_aux_ompc(ip,jp,kp,3)
    wp = w_aux_ompc(ip,jp,kp,4)
    qqp = 0.5_rkind * (up*up +vp*vp +wp*wp)
    hp = w_aux_ompc(ip,jp,kp,5)
    r = w_aux_ompc(ip,jp,kp,1)/w_aux_ompc(i,j,k,1)
    r = sqrt(r)
    rho = r*w_aux_ompc(i,j,k,1)
    rp1 = 1._rkind/(r +1._rkind)
    uu = (r*up +uu)*rp1
    vv = (r*vp +vv)*rp1
    ww = (r*wp +ww)*rp1
    h = (r*hp +h)*rp1
    qq = 0.5_rkind * (uu*uu +vv*vv + ww*ww)
    hbar = h - qq - cp_coeff_ompc(indx_cp_r+1)*t0
    if (calorically_perfect==1) then
      tt = t0+hbar/cp_coeff_ompc(0)
      gamloc = cp_coeff_ompc(0)/(cp_coeff_ompc(0)-rgas0)
    else
      tt = w_aux_ompc(i ,j ,k ,6)
      ttp = w_aux_ompc(ip,jp,kp,6)
      tt = (r*ttp +tt)*rp1
      told = tt
      do iter=1,max_iter
        num = 0._rkind
        den = 0._rkind
        do ll=indx_cp_l,indx_cp_r
          if (ll==-1) then
            tpow = (told/t0)**ll
            den = den+cp_coeff_ompc(ll)*tpow
            num = num+cp_coeff_ompc(ll)*log(told/t0)
          else
            tpow = (told/t0)**ll
            tpowp = (told/t0)*tpow
            den = den+cp_coeff_ompc(ll)*tpow
            num = num+cp_coeff_ompc(ll)*(tpowp-1._rkind)/(ll+1._rkind)
          endif
        enddo
        num = num*t0
        tt = told+(hbar-num)/den
        if (abs(tt-told) < tol_iter_nr) exit
        told = tt
      enddo
      cploc = 0._rkind
      do ll=indx_cp_l,indx_cp_r
        cploc = cploc+cp_coeff_ompc(ll)*(tt/t0)**ll
      enddo
      gamloc = cploc/(cploc-rgas0)
    endif
    cc = gamloc*tt*rgas0
    gm1loc = gamloc-1._rkind
    c = sqrt(cc)
    ci = 1._rkind/c
    p_rho = tt*rgas0
    p_e = rho*gm1loc
    etot = h - tt*rgas0
    b3 = etot - rho * p_rho/p_e
    b2 = p_e/(rho*cc)
    b1 = p_rho/cc - b2*(etot - 2._rkind*qq)
  endsubroutine compute_roe_average


  subroutine force_rhs_2_subroutine(nx,ny,nz,ng,fln_ompc,w_aux_ompc,bulk5g,fluid_mask_ompc)
    integer :: nx, ny, nz, ng
    real(rkind), dimension(5), intent(in) :: bulk5g
    real(rkind), dimension(1:,1:,1:,1:), intent(inout) :: fln_ompc
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in) :: w_aux_ompc
    integer, dimension(1-ng:,1-ng:,1-ng:), intent(in) :: fluid_mask_ompc
    real(rkind) :: bulk_1, bulk_2, uu
    integer :: i,j,k,iercuda
    bulk_1 = bulk5g(1)
    bulk_2 = bulk5g(2)
    !$omp parallel do default(firstprivate) shared(fln_ompc, &
    !$omp& w_aux_ompc, &
    !$omp& fluid_mask_ompc)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          if (fluid_mask_ompc(i,j,k)==0) then
            uu = w_aux_ompc(i,j,k,2)
            fln_ompc(i,j,k,1) = fln_ompc(i,j,k,1) - bulk_1
            fln_ompc(i,j,k,2) = fln_ompc(i,j,k,2) - bulk_2
            fln_ompc(i,j,k,5) = fln_ompc(i,j,k,5) - uu*bulk_2
          endif
        enddo
      enddo
    enddo


  endsubroutine force_rhs_2_subroutine



  subroutine force_rhs_2_c2_subroutine(nx,ny,nz,ng,fln_ompc,w_aux_ompc,bulk5g,fluid_mask_ompc,&
  &r_curv,yn_ompc,dxdcsinc2_ompc,dydcsinc2_ompc)
    integer :: nx, ny, nz, ng
    real(rkind), dimension(5), intent(in) :: bulk5g
    real(rkind), intent(in) :: r_curv
    real(rkind), dimension(1:,1:,1:,1:), intent(inout) :: fln_ompc
    real(rkind), dimension(1:), intent(in) :: yn_ompc
    real(rkind), dimension(1-ng:,1-ng:), intent(in) :: dxdcsinc2_ompc, dydcsinc2_ompc
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in) :: w_aux_ompc
    integer, dimension(1-ng:,1-ng:,1-ng:), intent(in) :: fluid_mask_ompc
    real(rkind) :: bulk_1, bulk_5, uu, vv, r, dpdcsi, dpdx, dpdy
    integer :: i,j,k,iercuda
    bulk_1 = bulk5g(1)
    bulk_5 = bulk5g(5)
    !$omp parallel do default(firstprivate) shared(fln_ompc, &
    !$omp& yn_ompc, &
    !$omp& dxdcsinc2_ompc, &
    !$omp& dydcsinc2_ompc, &
    !$omp& w_aux_ompc, &
    !$omp& fluid_mask_ompc)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          if (fluid_mask_ompc(i,j,k)==0) then
            uu = w_aux_ompc(i,j,k,2)
            vv = w_aux_ompc(i,j,k,3)
            r = r_curv+0.5_rkind*(yn_ompc(j)+yn_ompc(j+1))
            dpdcsi = bulk_5/r
            dpdx = dpdcsi*dxdcsinc2_ompc(i,j)
            dpdy = dpdcsi*dydcsinc2_ompc(i,j)
            fln_ompc(i,j,k,1) = fln_ompc(i,j,k,1) - bulk_1
            fln_ompc(i,j,k,2) = fln_ompc(i,j,k,2) - dpdx
            fln_ompc(i,j,k,3) = fln_ompc(i,j,k,3) - dpdy
            fln_ompc(i,j,k,5) = fln_ompc(i,j,k,5) - uu*dpdx - vv*dpdy
          endif
        enddo
      enddo
    enddo


  endsubroutine force_rhs_2_c2_subroutine


  subroutine force_rhs_1_subroutine(nx,ny,nz,ng,yn_ompc,fln_ompc,w_ompc,w_aux_ompc,bulk5,fluid_mask_ompc)
    integer, intent(in) :: nx, ny, nz, ng
    real(rkind), dimension(1:), intent(in) :: yn_ompc
    real(rkind), dimension(1:,1:,1:,1:), intent(in) :: fln_ompc
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in) :: w_ompc
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in) :: w_aux_ompc
    integer, dimension(1-ng:,1-ng:,1-ng:), intent(in) :: fluid_mask_ompc
    real(rkind), dimension(5), intent(out) :: bulk5
    real(rkind) :: bulk_1, bulk_2, bulk_3, bulk_4, bulk_5
    real(rkind) :: dy
    integer :: i,j,k,iercuda
    bulk_1 = 0._rkind
    bulk_2 = 0._rkind
    bulk_3 = 0._rkind
    bulk_4 = 0._rkind
    bulk_5 = 0._rkind
    !$omp parallel do default(firstprivate) shared(yn_ompc, &
    !$omp& fln_ompc, &
    !$omp& w_ompc, &
    !$omp& w_aux_ompc, &
    !$omp& fluid_mask_ompc) &
    !$omp& reduction(+:bulk_1,bulk_2,bulk_3,bulk_4,bulk_5)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          if (fluid_mask_ompc(i,j,k)==0) then
            dy = yn_ompc(j+1)-yn_ompc(j)
            bulk_1 = bulk_1 + fln_ompc(i,j,k,1)*dy
            bulk_2 = bulk_2 + fln_ompc(i,j,k,2)*dy
            bulk_3 = bulk_3 + w_ompc (i,j,k,1)*dy
            bulk_4 = bulk_4 + w_ompc (i,j,k,2)*dy
            bulk_5 = bulk_5 + w_ompc (i,j,k,2)*dy*w_aux_ompc(i,j,k,6)
          endif
        enddo
      enddo
    enddo
    bulk5(1) = bulk_1
    bulk5(2) = bulk_2
    bulk5(3) = bulk_3
    bulk5(4) = bulk_4
    bulk5(5) = bulk_5

  endsubroutine force_rhs_1_subroutine



  subroutine force_rhs_1_c2_subroutine(nx,ny,nz,ng,yn_ompc,fln_ompc,w_ompc,w_aux_ompc,bulk5,&
  &fluid_mask_ompc,dcsidxnc2_ompc,dcsidync2_ompc,jac_ompc)
    integer, intent(in) :: nx, ny, nz, ng
    real(rkind), dimension(1:), intent(in) :: yn_ompc
    real(rkind), dimension(1:,1:,1:,1:), intent(in) :: fln_ompc
    real(rkind), dimension(1-ng:,1-ng:), intent(in) :: dcsidxnc2_ompc, dcsidync2_ompc, jac_ompc
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in) :: w_ompc
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in) :: w_aux_ompc
    integer, dimension(1-ng:,1-ng:,1-ng:), intent(in) :: fluid_mask_ompc
    real(rkind), dimension(5), intent(out) :: bulk5
    real(rkind) :: bulk_1, bulk_2, bulk_3, bulk_4, bulk_5
    real(rkind) :: dv, rhou, f_rhou
    integer :: i,j,k,iercuda
    bulk_1 = 0._rkind
    bulk_2 = 0._rkind
    bulk_3 = 0._rkind
    bulk_4 = 0._rkind
    bulk_5 = 0._rkind
    !$omp parallel do default(firstprivate) shared(yn_ompc, &
    !$omp& fln_ompc, &
    !$omp& dcsidxnc2_ompc, &
    !$omp& dcsidync2_ompc, &
    !$omp& jac_ompc, &
    !$omp& w_ompc, &
    !$omp& w_aux_ompc, &
    !$omp& fluid_mask_ompc) &
    !$omp& reduction(+:bulk_1,bulk_2,bulk_3,bulk_4,bulk_5)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          if (fluid_mask_ompc(i,j,k)==0) then
            dv = 1._rkind/jac_ompc(i,j)
            rhou = w_ompc(i,j,k,2)*dcsidxnc2_ompc(i,j)+ w_ompc(i,j,k,3)*dcsidync2_ompc(i,j)
            f_rhou = fln_ompc(i,j,k,2)*dcsidxnc2_ompc(i,j)+fln_ompc(i,j,k,3)*dcsidync2_ompc(i,j)
            bulk_1 = bulk_1 + fln_ompc(i,j,k,1)*dv
            bulk_2 = bulk_2 + w_ompc(i,j,k,1)*dv
            bulk_3 = bulk_3 + rhou*dv
            bulk_4 = bulk_4 + rhou*w_aux_ompc(i,j,k,6)*dv
            bulk_5 = bulk_5 + f_rhou*dv
          endif
        enddo
      enddo
    enddo
    bulk5(1) = bulk_1
    bulk5(2) = bulk_2
    bulk5(3) = bulk_3
    bulk5(4) = bulk_4
    bulk5(5) = bulk_5

  endsubroutine force_rhs_1_c2_subroutine



  subroutine force_var_1_subroutine(nx,ny,nz,ng,yn_ompc,fln_ompc,w_ompc,w_aux_ompc,bulkt,&
  &fluid_mask_ompc,cv_coeff_ompc,indx_cp_l,indx_cp_r,t0,calorically_perfect,tol_iter_nr)
    integer, intent(in) :: nx, ny, nz, ng, indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind), intent(in) :: t0, tol_iter_nr
    real(rkind), dimension(1:), intent(in) :: yn_ompc
    real(rkind), dimension(1:,1:,1:,1:), intent(in) :: fln_ompc
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in) :: w_ompc
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(inout) :: w_aux_ompc
    integer, dimension(1-ng:,1-ng:,1-ng:), intent(in) :: fluid_mask_ompc
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in) :: cv_coeff_ompc
    real(rkind), intent(out) :: bulkt
    real(rkind) :: bulk_5,rho,rhou,rhov,rhow,rhoe,ri,uu,vv,ww,qq,ee,tt
    real(rkind) :: dy
    integer :: i,j,k,iercuda
    bulk_5 = 0._rkind
    !$omp parallel do default(firstprivate) shared(yn_ompc, &
    !$omp& fln_ompc, &
    !$omp& w_ompc, &
    !$omp& w_aux_ompc, &
    !$omp& cv_coeff_ompc, &
    !$omp& fluid_mask_ompc) &
    !$omp& reduction(+:bulk_5)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          if (fluid_mask_ompc(i,j,k)==0) then
            dy = yn_ompc(j+1)-yn_ompc(j)
            rho = w_ompc(i,j,k,1)
            rhou = w_ompc(i,j,k,2)
            rhov = w_ompc(i,j,k,3)
            rhow = w_ompc(i,j,k,4)
            rhoe = w_ompc(i,j,k,5)
            ri = 1._rkind/rho
            uu = rhou*ri
            vv = rhov*ri
            ww = rhow*ri
            qq = 0.5_rkind*(uu*uu+vv*vv+ww*ww)
            ee = rhoe/rho-qq
            tt = get_temperature_from_e_dev(ee,w_aux_ompc(i,j,k,6),t0,cv_coeff_ompc,indx_cp_l,&
            &indx_cp_r, calorically_perfect,tol_iter_nr)
            w_aux_ompc(i,j,k,6) = tt
            bulk_5 = bulk_5 + rhou*tt*dy
          endif
        enddo
      enddo
    enddo
    bulkt = bulk_5

  endsubroutine force_var_1_subroutine


  subroutine force_var_2_subroutine(nx,ny,nz,ng,w_ompc,w_aux_ompc,tbdiff,fluid_mask_ompc,&
  &cv_coeff_ompc,indx_cp_l,indx_cp_r,t0,calorically_perfect)
    integer :: nx, ny, nz, ng, indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind), intent(in) :: tbdiff, t0
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(inout) :: w_ompc
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in) :: w_aux_ompc
    integer, dimension(1-ng:,1-ng:,1-ng:), intent(in) :: fluid_mask_ompc
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in) :: cv_coeff_ompc
    real(rkind) :: rho,rhou,rhov,rhow,tt,ttnew,ee
    integer :: i,j,k,iercuda

    !$omp parallel do default(firstprivate) shared(w_ompc, &
    !$omp& w_aux_ompc, &
    !$omp& cv_coeff_ompc, &
    !$omp& fluid_mask_ompc)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          if (fluid_mask_ompc(i,j,k)==0) then
            rho = w_ompc(i,j,k,1)
            rhou = w_ompc(i,j,k,2)
            rhov = w_ompc(i,j,k,3)
            rhow = w_ompc(i,j,k,4)
            tt = w_aux_ompc(i,j,k,6)
            ttnew = tt+tbdiff
            ee = get_e_from_temperature_dev(ttnew, t0, indx_cp_l, indx_cp_r, cv_coeff_ompc, calorically_perfect)
            w_ompc(i,j,k,5) = rho*ee + 0.5_rkind*(rhou**2+rhov**2+rhow**2)/rho
          endif
        enddo
      enddo
    enddo


  endsubroutine force_var_2_subroutine



  subroutine force_var_1_c2_subroutine(nx,ny,nz,ng,yn_ompc,fln_ompc,w_ompc,w_aux_ompc,bulkt,&
  &fluid_mask_ompc,cv_coeff_ompc,indx_cp_l,indx_cp_r,t0,calorically_perfect,tol_iter_nr,jac_ompc,&
  &dcsidxnc2_ompc,dcsidync2_ompc)
    integer, intent(in) :: nx, ny, nz, ng, indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind), intent(in) :: t0, tol_iter_nr
    real(rkind), dimension(1:), intent(in) :: yn_ompc
    real(rkind), dimension(1:,1:,1:,1:), intent(in) :: fln_ompc
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in) :: w_ompc
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(inout) :: w_aux_ompc
    integer, dimension(1-ng:,1-ng:,1-ng:), intent(in) :: fluid_mask_ompc
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in) :: cv_coeff_ompc
    real(rkind), dimension(1-ng:,1-ng:), intent(in) :: jac_ompc, dcsidxnc2_ompc, dcsidync2_ompc
    real(rkind), intent(out) :: bulkt
    real(rkind) :: bulk_5,rho,rhou,rhov,rhow,rhoe,ri,uu,vv,ww,qq,ee,tt
    real(rkind) :: dy, dv, rhout
    integer :: i,j,k,iercuda
    bulk_5 = 0._rkind
    !$omp parallel do default(firstprivate) shared(yn_ompc, &
    !$omp& fln_ompc, &
    !$omp& w_ompc, &
    !$omp& w_aux_ompc, &
    !$omp& cv_coeff_ompc, &
    !$omp& jac_ompc, &
    !$omp& dcsidxnc2_ompc, &
    !$omp& dcsidync2_ompc, &
    !$omp& fluid_mask_ompc) &
    !$omp& reduction(+:bulk_5)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          if (fluid_mask_ompc(i,j,k)==0) then
            dv = 1._rkind/jac_ompc(i,j)
            rho = w_ompc(i,j,k,1)
            rhou = w_ompc(i,j,k,2)
            rhov = w_ompc(i,j,k,3)
            rhow = w_ompc(i,j,k,4)
            rhoe = w_ompc(i,j,k,5)
            rhout = rhou*dcsidxnc2_ompc(i,j)+rhov*dcsidync2_ompc(i,j)
            ri = 1._rkind/rho
            uu = rhou*ri
            vv = rhov*ri
            ww = rhow*ri
            qq = 0.5_rkind*(uu*uu+vv*vv+ww*ww)
            ee = rhoe/rho-qq
            tt = get_temperature_from_e_dev(ee,w_aux_ompc(i,j,k,6),t0,cv_coeff_ompc,indx_cp_l,&
            &indx_cp_r, calorically_perfect,tol_iter_nr)
            w_aux_ompc(i,j,k,6) = tt
            bulk_5 = bulk_5 + rhou*tt*dv
          endif
        enddo
      enddo
    enddo
    bulkt = bulk_5

  endsubroutine force_var_1_c2_subroutine


  subroutine force_var_2_c2_subroutine(nx,ny,nz,ng,w_ompc,w_aux_ompc,tbdiff,fluid_mask_ompc,&
  &cv_coeff_ompc,indx_cp_l,indx_cp_r,t0,calorically_perfect,jac_ompc)
    integer :: nx, ny, nz, ng, indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind), intent(in) :: tbdiff, t0
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(inout) :: w_ompc
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in) :: w_aux_ompc
    integer, dimension(1-ng:,1-ng:,1-ng:), intent(in) :: fluid_mask_ompc
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in) :: cv_coeff_ompc
    real(rkind), dimension(1-ng:,1-ng:), intent(in) :: jac_ompc
    real(rkind) :: rho,rhou,rhov,rhow,tt,ttnew,ee,dv
    integer :: i,j,k,iercuda

    !$omp parallel do default(firstprivate) shared(w_ompc, &
    !$omp& w_aux_ompc, &
    !$omp& cv_coeff_ompc, &
    !$omp& jac_ompc, &
    !$omp& fluid_mask_ompc)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          if (fluid_mask_ompc(i,j,k)==0) then
            dv = 1._rkind/jac_ompc(i,j)
            rho = w_ompc(i,j,k,1)
            rhou = w_ompc(i,j,k,2)
            rhov = w_ompc(i,j,k,3)
            rhow = w_ompc(i,j,k,4)
            tt = w_aux_ompc(i,j,k,6)
            ttnew = tt+tbdiff
            ee = get_e_from_temperature_dev(ttnew, t0, indx_cp_l, indx_cp_r, cv_coeff_ompc, calorically_perfect)
            w_ompc(i,j,k,5) = rho*ee + 0.5_rkind*(rhou**2+rhov**2+rhow**2)/rho
          endif
        enddo
      enddo
    enddo


  endsubroutine force_var_2_c2_subroutine



  subroutine update_flux_subroutine(nx,ny,nz,nv,fl_ompc,fln_ompc,gamdt)
    integer :: nx, ny, nz, nv
    real(rkind) :: gamdt
    real(rkind), dimension(1:,1:,1:,1:), intent(in) :: fl_ompc
    real(rkind), dimension(1:,1:,1:,1:), intent(inout) :: fln_ompc
    integer :: i,j,k,m,iercuda

    !$omp parallel do default(firstprivate) shared(fl_ompc, &
    !$omp& fln_ompc)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          do m=1,nv
            fln_ompc(i,j,k,m) = fln_ompc(i,j,k,m)-gamdt*fl_ompc(i,j,k,m)
          enddo
        enddo
      enddo
    enddo


  endsubroutine update_flux_subroutine



  subroutine update_field_subroutine(nx,ny,nz,ng,nv,w_ompc,fln_ompc,fluid_mask_ompc)
    integer :: nx, ny, nz, nv, ng
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(inout) :: w_ompc
    real(rkind), dimension(1:,1:,1:,1:), intent(in) :: fln_ompc
    integer, dimension(1-ng:,1-ng:,1-ng:), intent(in) :: fluid_mask_ompc
    integer :: i,j,k,m,iercuda

    !$omp parallel do default(firstprivate) shared(w_ompc, &
    !$omp& fln_ompc, &
    !$omp& fluid_mask_ompc)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          if (fluid_mask_ompc(i,j,k)==0) then
            do m=1,nv
              w_ompc(i,j,k,m) = w_ompc(i,j,k,m)+fln_ompc(i,j,k,m)
            enddo
          endif
        enddo
      enddo
    enddo


  endsubroutine update_field_subroutine



  subroutine visflx_div_ord2_subroutine(nx,ny,nz,ng,w_aux_ompc,fl_ompc,x_ompc,y_ompc,z_ompc,stream_id)
    integer, intent(in) :: nx, ny, nz, ng
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(inout) :: w_aux_ompc
    real(rkind), dimension(1:,1:,1:, 1:), intent(inout) :: fl_ompc
    real(rkind), dimension(1-ng:), intent(in) :: x_ompc, y_ompc, z_ompc
    integer :: i,j,k,iercuda
    real(rkind) :: dxl,dyl,dzl
    real(rkind) :: uu,vv,ww,mu
    real(rkind) :: sigq,sigx,sigy,sigz
    integer, intent(in) :: stream_id

    !$omp parallel do default(firstprivate) shared(w_aux_ompc, &
    !$omp& fl_ompc, &
    !$omp& x_ompc, &
    !$omp& y_ompc, &
    !$omp& z_ompc)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          uu = w_aux_ompc(i,j,k,2)
          vv = w_aux_ompc(i,j,k,3)
          ww = w_aux_ompc(i,j,k,4)
          mu = w_aux_ompc(i,j,k,7)
          dxl = mu/(x_ompc(i+1)-x_ompc(i-1))
          dyl = mu/(y_ompc(j+1)-y_ompc(j-1))
          dzl = mu/(z_ompc(k+1)-z_ompc(k-1))
          sigx = dxl*(w_aux_ompc(i+1,j,k,10)-w_aux_ompc(i-1,j,k,10))
          sigy = dyl*(w_aux_ompc(i,j+1,k,10)-w_aux_ompc(i,j-1,k,10))
          sigz = dzl*(w_aux_ompc(i,j,k+1,10)-w_aux_ompc(i,j,k-1,10))
          sigq = sigx*uu+sigy*vv+sigz*ww
          fl_ompc(i,j,k,2) = fl_ompc(i,j,k,2) - sigx
          fl_ompc(i,j,k,3) = fl_ompc(i,j,k,3) - sigy
          fl_ompc(i,j,k,4) = fl_ompc(i,j,k,4) - sigz
          fl_ompc(i,j,k,5) = fl_ompc(i,j,k,5) - sigq
        enddo
      enddo
    enddo


  endsubroutine visflx_div_ord2_subroutine



  subroutine visflx_div_subroutine(nx,ny,nz,ng,visc_order,w_aux_ompc,fl_ompc,coeff_deriv1_ompc,&
  &dcsidx_ompc,detady_ompc,dzitdz_ompc,stream_id)
    integer, intent(in) :: nx, ny, nz, ng, visc_order
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(inout) :: w_aux_ompc
    real(rkind), dimension(1:,1:,1:, 1:), intent(inout) :: fl_ompc
    real(rkind), dimension(1:,1:), intent(in) :: coeff_deriv1_ompc
    real(rkind), dimension(1:), intent(in) :: dcsidx_ompc, detady_ompc, dzitdz_ompc
    integer :: i,j,k,l,iercuda
    real(rkind) :: ccl
    real(rkind) :: uu,vv,ww,tt,mu
    real(rkind) :: sigq,sigx,sigy,sigz
    real(rkind) :: divx3l, divy3l, divz3l
    integer :: lmax
    integer, intent(in) :: stream_id
    lmax = visc_order/2
    !$omp parallel do default(firstprivate) shared(w_aux_ompc, &
    !$omp& fl_ompc, &
    !$omp& coeff_deriv1_ompc, &
    !$omp& dcsidx_ompc, &
    !$omp& detady_ompc, &
    !$omp& dzitdz_ompc)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          uu = w_aux_ompc(i,j,k,2)
          vv = w_aux_ompc(i,j,k,3)
          ww = w_aux_ompc(i,j,k,4)
          tt = w_aux_ompc(i,j,k,6)
          mu = w_aux_ompc(i,j,k,7)
          divx3l = 0._rkind
          divy3l = 0._rkind
          divz3l = 0._rkind
          do l=1,lmax
            ccl = coeff_deriv1_ompc(l,lmax)
            divx3l = divx3l+ccl*(w_aux_ompc(i+l,j,k,10)-w_aux_ompc(i-l,j,k,10))
            divy3l = divy3l+ccl*(w_aux_ompc(i,j+l,k,10)-w_aux_ompc(i,j-l,k,10))
            divz3l = divz3l+ccl*(w_aux_ompc(i,j,k+l,10)-w_aux_ompc(i,j,k-l,10))
          enddo
          divx3l = divx3l*dcsidx_ompc(i)
          divy3l = divy3l*detady_ompc(j)
          divz3l = divz3l*dzitdz_ompc(k)
          sigx = mu*divx3l
          sigy = mu*divy3l
          sigz = mu*divz3l
          sigq = sigx*uu+sigy*vv+sigz*ww
          fl_ompc(i,j,k,2) = fl_ompc(i,j,k,2) - sigx
          fl_ompc(i,j,k,3) = fl_ompc(i,j,k,3) - sigy
          fl_ompc(i,j,k,4) = fl_ompc(i,j,k,4) - sigz
          fl_ompc(i,j,k,5) = fl_ompc(i,j,k,5) - sigq
        enddo
      enddo
    enddo


  endsubroutine visflx_div_subroutine



  subroutine visflx_div_c2_subroutine(nx,ny,nz,ng,visc_order,w_aux_ompc,fl_ompc,coeff_deriv1_ompc,&
  &dcsidxc2_ompc,dcsidyc2_ompc,detadxc2_ompc,detadyc2_ompc,dzitdz_ompc,vis_tag_ompc,wall_tag_ompc,&
  &stream_id)
    integer, intent(in) :: nx, ny, nz, ng, visc_order
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(inout) :: w_aux_ompc
    real(rkind), dimension(1:,1:,1:, 1:), intent(inout) :: fl_ompc
    real(rkind), dimension(1:,1:), intent(in) :: coeff_deriv1_ompc
    real(rkind), dimension(1-ng:,1-ng:), intent(in) :: dcsidxc2_ompc, detadyc2_ompc
    real(rkind), dimension(1-ng:,1-ng:), intent(in) :: dcsidyc2_ompc, detadxc2_ompc
    real(rkind), dimension(1:), intent(in) :: dzitdz_ompc
    integer, dimension(1-ng:), intent(in) :: vis_tag_ompc, wall_tag_ompc
    integer :: i,j,k,l,iercuda
    real(rkind) :: cli, clj, clk
    real(rkind) :: uu,vv,ww,tt,mu
    real(rkind) :: sigq,sigx,sigy,sigz
    real(rkind) :: divcsi, diveta, divzit
    real(rkind) :: divx3l, divy3l, divz3l
    integer :: lmax,lmaxi,lmaxj
    integer, intent(in) :: stream_id
    lmax = visc_order/2
    !$omp parallel do default(firstprivate) shared(w_aux_ompc, &
    !$omp& fl_ompc, &
    !$omp& coeff_deriv1_ompc, &
    !$omp& dcsidxc2_ompc, &
    !$omp& detadyc2_ompc, &
    !$omp& dcsidyc2_ompc, &
    !$omp& detadxc2_ompc, &
    !$omp& dzitdz_ompc, &
    !$omp& vis_tag_ompc, &
    !$omp& wall_tag_ompc)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          uu = w_aux_ompc(i,j,k,2)
          vv = w_aux_ompc(i,j,k,3)
          ww = w_aux_ompc(i,j,k,4)
          tt = w_aux_ompc(i,j,k,6)
          mu = w_aux_ompc(i,j,k,7)
          lmaxi = lmax
          lmaxj = lmax
          if (j==1) lmaxi = vis_tag_ompc(i)
          if (wall_tag_ompc(i) < 1) lmaxj = min(j,lmax)
          divcsi = 0._rkind
          diveta = 0._rkind
          divzit = 0._rkind
          do l=1,lmax
            cli = coeff_deriv1_ompc(l,lmaxi)
            clj = coeff_deriv1_ompc(l,lmaxj)
            clk = coeff_deriv1_ompc(l,lmax )
            divcsi = divcsi+cli*(w_aux_ompc(i+l,j,k,10)-w_aux_ompc(i-l,j,k,10))
            diveta = diveta+clj*(w_aux_ompc(i,j+l,k,10)-w_aux_ompc(i,j-l,k,10))
            divzit = divzit+clk*(w_aux_ompc(i,j,k+l,10)-w_aux_ompc(i,j,k-l,10))
          enddo
          divx3l = divcsi*dcsidxc2_ompc(i,j) + diveta*detadxc2_ompc(i,j)
          divy3l = divcsi*dcsidyc2_ompc(i,j) + diveta*detadyc2_ompc(i,j)
          divz3l = divzit*dzitdz_ompc(k)
          sigx = mu*divx3l
          sigy = mu*divy3l
          sigz = mu*divz3l
          sigq = sigx*uu+sigy*vv+sigz*ww
          fl_ompc(i,j,k,2) = fl_ompc(i,j,k,2) - sigx
          fl_ompc(i,j,k,3) = fl_ompc(i,j,k,3) - sigy
          fl_ompc(i,j,k,4) = fl_ompc(i,j,k,4) - sigz
          fl_ompc(i,j,k,5) = fl_ompc(i,j,k,5) - sigq
        enddo
      enddo
    enddo


  endsubroutine visflx_div_c2_subroutine



  subroutine visflx_subroutine(nx,ny,nz,ng,visc_order,prandtl,t0,indx_cp_l,indx_cp_r,cp_coeff_ompc,&
  &calorically_perfect,u0,l0,w_ompc,w_aux_ompc,fl_ompc,coeff_deriv1_ompc,coeff_deriv2_ompc,dcsidx_ompc,&
  &detady_ompc,dzitdz_ompc,dcsidxs_ompc,detadys_ompc,dzitdzs_ompc,dcsidx2_ompc,detady2_ompc,&
  &dzitdz2_ompc,wallprop_ompc)
    integer, intent(in) :: nx, ny, nz, ng, visc_order, calorically_perfect
    integer, intent(in) :: indx_cp_l, indx_cp_r
    real(rkind), intent(in) :: prandtl, u0, l0, t0
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(in) :: w_ompc
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(inout) :: w_aux_ompc
    real(rkind), dimension(1-ng:,1-ng:, 2:), intent(inout) :: wallprop_ompc
    real(rkind), dimension(1:,1:,1:, 1:), intent(inout) :: fl_ompc
    real(rkind), dimension(1:,1:), intent(in) :: coeff_deriv1_ompc
    real(rkind), dimension(0:,1:), intent(in) :: coeff_deriv2_ompc
    real(rkind), dimension(1:), intent(in) :: dcsidx_ompc, detady_ompc, dzitdz_ompc
    real(rkind), dimension(1:), intent(in) :: dcsidxs_ompc, detadys_ompc, dzitdzs_ompc
    real(rkind), dimension(1:), intent(in) :: dcsidx2_ompc, detady2_ompc, dzitdz2_ompc
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in) :: cp_coeff_ompc
    integer :: i,j,k,l,ll,iercuda
    real(rkind) :: ccl,clapl
    real(rkind) :: sig11,sig12,sig13
    real(rkind) :: sig22,sig23
    real(rkind) :: sig33
    real(rkind) :: uu,vv,ww,tt,mu
    real(rkind) :: ux,uy,uz
    real(rkind) :: vx,vy,vz
    real(rkind) :: wx,wy,wz
    real(rkind) :: tx,ty,tz
    real(rkind) :: mux, muy, muz
    real(rkind) :: ulap,ulapx,ulapy,ulapz
    real(rkind) :: vlap,vlapx,vlapy,vlapz
    real(rkind) :: wlap,wlapx,wlapy,wlapz
    real(rkind) :: tlap,tlapx,tlapy,tlapz
    real(rkind) :: sigq,sigx,sigy,sigz,sigqt,sigah
    real(rkind) :: div, div3l, omegax, omegay, omegaz, omod2
    real(rkind) :: cploc
    integer :: lmax
    lmax = visc_order/2
    !$omp parallel do default(firstprivate) shared(w_ompc, &
    !$omp& w_aux_ompc, &
    !$omp& wallprop_ompc, &
    !$omp& fl_ompc, &
    !$omp& coeff_deriv1_ompc, &
    !$omp& coeff_deriv2_ompc, &
    !$omp& dcsidx_ompc, &
    !$omp& detady_ompc, &
    !$omp& dzitdz_ompc, &
    !$omp& dcsidxs_ompc, &
    !$omp& detadys_ompc, &
    !$omp& dzitdzs_ompc, &
    !$omp& dcsidx2_ompc, &
    !$omp& detady2_ompc, &
    !$omp& dzitdz2_ompc, &
    !$omp& cp_coeff_ompc)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          uu = w_aux_ompc(i,j,k,2)
          vv = w_aux_ompc(i,j,k,3)
          ww = w_aux_ompc(i,j,k,4)
          tt = w_aux_ompc(i,j,k,6)
          mu = w_aux_ompc(i,j,k,7)
          if (calorically_perfect==1) then
            cploc = cp_coeff_ompc(0)
          else
            cploc = 0._rkind
            do ll=indx_cp_l,indx_cp_r
              cploc = cploc+cp_coeff_ompc(ll)*(tt/t0)**ll
            enddo
          endif
          ux = 0._rkind
          vx = 0._rkind
          wx = 0._rkind
          tx = 0._rkind
          mux = 0._rkind
          uy = 0._rkind
          vy = 0._rkind
          wy = 0._rkind
          ty = 0._rkind
          muy = 0._rkind
          uz = 0._rkind
          vz = 0._rkind
          wz = 0._rkind
          tz = 0._rkind
          muz = 0._rkind
          do l=1,lmax
            ccl = coeff_deriv1_ompc(l,lmax)
            ux = ux+ccl*(w_aux_ompc(i+l,j,k,2)-w_aux_ompc(i-l,j,k,2))
            vx = vx+ccl*(w_aux_ompc(i+l,j,k,3)-w_aux_ompc(i-l,j,k,3))
            wx = wx+ccl*(w_aux_ompc(i+l,j,k,4)-w_aux_ompc(i-l,j,k,4))
            tx = tx+ccl*(w_aux_ompc(i+l,j,k,6)-w_aux_ompc(i-l,j,k,6))
            mux = mux+ccl*(w_aux_ompc(i+l,j,k,7)-w_aux_ompc(i-l,j,k,7))
            uy = uy+ccl*(w_aux_ompc(i,j+l,k,2)-w_aux_ompc(i,j-l,k,2))
            vy = vy+ccl*(w_aux_ompc(i,j+l,k,3)-w_aux_ompc(i,j-l,k,3))
            wy = wy+ccl*(w_aux_ompc(i,j+l,k,4)-w_aux_ompc(i,j-l,k,4))
            ty = ty+ccl*(w_aux_ompc(i,j+l,k,6)-w_aux_ompc(i,j-l,k,6))
            muy = muy+ccl*(w_aux_ompc(i,j+l,k,7)-w_aux_ompc(i,j-l,k,7))
            uz = uz+ccl*(w_aux_ompc(i,j,k+l,2)-w_aux_ompc(i,j,k-l,2))
            vz = vz+ccl*(w_aux_ompc(i,j,k+l,3)-w_aux_ompc(i,j,k-l,3))
            wz = wz+ccl*(w_aux_ompc(i,j,k+l,4)-w_aux_ompc(i,j,k-l,4))
            tz = tz+ccl*(w_aux_ompc(i,j,k+l,6)-w_aux_ompc(i,j,k-l,6))
            muz = muz+ccl*(w_aux_ompc(i,j,k+l,7)-w_aux_ompc(i,j,k-l,7))
          enddo
          ux = ux*dcsidx_ompc(i)
          vx = vx*dcsidx_ompc(i)
          wx = wx*dcsidx_ompc(i)
          tx = tx*dcsidx_ompc(i)
          mux = mux*dcsidx_ompc(i)
          uy = uy*detady_ompc(j)
          vy = vy*detady_ompc(j)
          wy = wy*detady_ompc(j)
          ty = ty*detady_ompc(j)
          muy = muy*detady_ompc(j)
          uz = uz*dzitdz_ompc(k)
          vz = vz*dzitdz_ompc(k)
          wz = wz*dzitdz_ompc(k)
          tz = tz*dzitdz_ompc(k)
          muz = muz*dzitdz_ompc(k)
          if (j==1) then
            wallprop_ompc(i,k,2) = mu*uy
            wallprop_ompc(i,k,3) = mu*wy
            wallprop_ompc(i,k,4) = mu*ty*cploc/prandtl
          endif
          ulapx = coeff_deriv2_ompc(0,lmax)*uu
          ulapy = ulapx
          ulapz = ulapx
          vlapx = coeff_deriv2_ompc(0,lmax)*vv
          vlapy = vlapx
          vlapz = vlapx
          wlapx = coeff_deriv2_ompc(0,lmax)*ww
          wlapy = wlapx
          wlapz = wlapx
          tlapx = coeff_deriv2_ompc(0,lmax)*tt
          tlapy = tlapx
          tlapz = tlapx
          do l=1,lmax
            clapl = coeff_deriv2_ompc(l,lmax)
            ulapx = ulapx + clapl*(w_aux_ompc(i+l,j,k,2)+w_aux_ompc(i-l,j,k,2))
            ulapy = ulapy + clapl*(w_aux_ompc(i,j+l,k,2)+w_aux_ompc(i,j-l,k,2))
            ulapz = ulapz + clapl*(w_aux_ompc(i,j,k+l,2)+w_aux_ompc(i,j,k-l,2))
            vlapx = vlapx + clapl*(w_aux_ompc(i+l,j,k,3)+w_aux_ompc(i-l,j,k,3))
            vlapy = vlapy + clapl*(w_aux_ompc(i,j+l,k,3)+w_aux_ompc(i,j-l,k,3))
            vlapz = vlapz + clapl*(w_aux_ompc(i,j,k+l,3)+w_aux_ompc(i,j,k-l,3))
            wlapx = wlapx + clapl*(w_aux_ompc(i+l,j,k,4)+w_aux_ompc(i-l,j,k,4))
            wlapy = wlapy + clapl*(w_aux_ompc(i,j+l,k,4)+w_aux_ompc(i,j-l,k,4))
            wlapz = wlapz + clapl*(w_aux_ompc(i,j,k+l,4)+w_aux_ompc(i,j,k-l,4))
            tlapx = tlapx + clapl*(w_aux_ompc(i+l,j,k,6)+w_aux_ompc(i-l,j,k,6))
            tlapy = tlapy + clapl*(w_aux_ompc(i,j+l,k,6)+w_aux_ompc(i,j-l,k,6))
            tlapz = tlapz + clapl*(w_aux_ompc(i,j,k+l,6)+w_aux_ompc(i,j,k-l,6))
          enddo
          ulapx = ulapx*dcsidxs_ompc(i)+ux*dcsidx2_ompc(i)
          vlapx = vlapx*dcsidxs_ompc(i)+vx*dcsidx2_ompc(i)
          wlapx = wlapx*dcsidxs_ompc(i)+wx*dcsidx2_ompc(i)
          tlapx = tlapx*dcsidxs_ompc(i)+tx*dcsidx2_ompc(i)
          ulapy = ulapy*detadys_ompc(j)+uy*detady2_ompc(j)
          vlapy = vlapy*detadys_ompc(j)+vy*detady2_ompc(j)
          wlapy = wlapy*detadys_ompc(j)+wy*detady2_ompc(j)
          tlapy = tlapy*detadys_ompc(j)+ty*detady2_ompc(j)
          ulapz = ulapz*dzitdzs_ompc(k)+uz*dzitdz2_ompc(k)
          vlapz = vlapz*dzitdzs_ompc(k)+vz*dzitdz2_ompc(k)
          wlapz = wlapz*dzitdzs_ompc(k)+wz*dzitdz2_ompc(k)
          tlapz = tlapz*dzitdzs_ompc(k)+tz*dzitdz2_ompc(k)
          ulap = ulapx+ulapy+ulapz
          vlap = vlapx+vlapy+vlapz
          wlap = wlapx+wlapy+wlapz
          tlap = tlapx+tlapy+tlapz
          div = ux+vy+wz
          div3l = div/3._rkind
          w_aux_ompc(i,j,k,10) = div3l
          omegax = wy-vz
          omegay = uz-wx
          omegaz = vx-uy
          omod2 = omegax*omegax+omegay*omegay+omegaz*omegaz
          w_aux_ompc(i,j,k,9) = sqrt(omod2)
          w_aux_ompc(i,j,k,8) = (max(-div/sqrt(omod2+div**2+(u0/l0)**2),0._rkind))**2
          sig11 = 2._rkind*(ux-div3l)
          sig12 = uy+vx
          sig13 = uz+wx
          sig22 = 2._rkind*(vy-div3l)
          sig23 = vz+wy
          sig33 = 2._rkind*(wz-div3l)
          sigx = mux*sig11 + muy*sig12 + muz*sig13 + mu*ulap
          sigy = mux*sig12 + muy*sig22 + muz*sig23 + mu*vlap
          sigz = mux*sig13 + muy*sig23 + muz*sig33 + mu*wlap
          sigqt = (mux*tx+muy*ty+muz*tz+mu*tlap)*cploc/prandtl
          sigah = (sig11*ux+sig12*uy+sig13*uz+sig12*vx+sig22*vy+sig23*vz+sig13*wx+sig23*wy+sig33*wz)*mu
          sigq = sigx*uu+sigy*vv+sigz*ww+sigah+sigqt
          fl_ompc(i,j,k,2) = fl_ompc(i,j,k,2) - sigx
          fl_ompc(i,j,k,3) = fl_ompc(i,j,k,3) - sigy
          fl_ompc(i,j,k,4) = fl_ompc(i,j,k,4) - sigz
          fl_ompc(i,j,k,5) = fl_ompc(i,j,k,5) - sigq
        enddo
      enddo
    enddo


  endsubroutine visflx_subroutine



  subroutine visflx_c2_subroutine(nx,ny,nz,ng,visc_order,prandtl,t0,indx_cp_l,indx_cp_r,&
  &cp_coeff_ompc,calorically_perfect,u0,l0,w_ompc,w_aux_ompc,fl_ompc,coeff_deriv1_ompc,&
  &coeff_deriv2_ompc,dcsidxc2_ompc,detadyc2_ompc,detadxc2_ompc,dcsidyc2_ompc,dzitdz_ompc,dzitdzs_ompc,&
  &dzitdz2_ompc,g1_ompc,g2_ompc,g12_ompc,jac_ompc,wallprop_ompc,vis_tag_ompc,wall_tag_ompc,iblock,&
  &ite_rank_x,itu_rank_x,ite_l,itu_l,teshk)
    integer, intent(in) :: nx, ny, nz, ng, visc_order, calorically_perfect
    integer, intent(in) :: indx_cp_l, indx_cp_r
    real(rkind), intent(in) :: prandtl, u0, l0, t0
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(in) :: w_ompc
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(inout) :: w_aux_ompc
    real(rkind), dimension(1-ng:,1-ng:, 2:), intent(inout) :: wallprop_ompc
    real(rkind), dimension(1:,1:,1:, 1:), intent(inout) :: fl_ompc
    real(rkind), dimension(1:,1:), intent(in) :: coeff_deriv1_ompc
    real(rkind), dimension(0:,1:), intent(in) :: coeff_deriv2_ompc
    real(rkind), dimension(1:), intent(in) :: dzitdz_ompc, dzitdzs_ompc, dzitdz2_ompc
    real(rkind), dimension(1-ng:,1-ng:), intent(in) :: dcsidxc2_ompc, detadyc2_ompc
    real(rkind), dimension(1-ng:,1-ng:), intent(in) :: detadxc2_ompc, dcsidyc2_ompc
    real(rkind), dimension(1-ng:,1-ng:), intent(in) :: g1_ompc, g2_ompc, g12_ompc, jac_ompc
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in) :: cp_coeff_ompc
    integer, dimension(1-ng:), intent(in) :: vis_tag_ompc, wall_tag_ompc
    integer, intent(in) :: iblock, ite_rank_x, itu_rank_x, ite_l, itu_l
    real(rkind), intent(in) :: teshk
    integer :: i,j,k,l,ll,iercuda
    real(rkind) :: ccl,clapl,cli,clj,clk,clapi,clapj,clapk
    real(rkind) :: sig11,sig12,sig13
    real(rkind) :: sig22,sig23
    real(rkind) :: sig33
    real(rkind) :: uu,vv,ww,tt,mu
    real(rkind) :: ux,uy,uz
    real(rkind) :: vx,vy,vz
    real(rkind) :: wx,wy,wz
    real(rkind) :: tx,ty,tz
    real(rkind) :: mux, muy, muz
    real(rkind) :: ucsi,ueta,uzit,vcsi,veta,vzit,wcsi,weta,wzit
    real(rkind) :: tcsi,teta,tzit,mucsi,mueta,muzit
    real(rkind) :: ulapcsi,ulapeta,ulapzit,vlapcsi,vlapeta,vlapzit,wlapcsi,wlapeta,wlapzit
    real(rkind) :: ulapcsieta,vlapcsieta,wlapcsieta,tlapcsieta
    real(rkind) :: dg1,dg2,dg12csi,dg12eta
    real(rkind) :: tlapcsi,tlapeta,tlapzit
    real(rkind) :: ulap,ulapx,ulapy,ulapz
    real(rkind) :: vlap,vlapx,vlapy,vlapz
    real(rkind) :: wlap,wlapx,wlapy,wlapz
    real(rkind) :: tlap,tlapx,tlapy,tlapz
    real(rkind) :: sigq,sigx,sigy,sigz,sigqt,sigah
    real(rkind) :: div, div3l, omegax, omegay, omegaz, omod2
    real(rkind) :: cploc
    integer :: lmax,lmaxi,lmaxj
    lmax = visc_order/2
    !$omp parallel do default(firstprivate) shared(w_ompc, &
    !$omp& w_aux_ompc, &
    !$omp& wallprop_ompc, &
    !$omp& fl_ompc, &
    !$omp& coeff_deriv1_ompc, &
    !$omp& coeff_deriv2_ompc, &
    !$omp& dzitdz_ompc, &
    !$omp& dzitdzs_ompc, &
    !$omp& dzitdz2_ompc, &
    !$omp& dcsidxc2_ompc, &
    !$omp& detadyc2_ompc, &
    !$omp& detadxc2_ompc, &
    !$omp& dcsidyc2_ompc, &
    !$omp& g1_ompc, &
    !$omp& g2_ompc, &
    !$omp& g12_ompc, &
    !$omp& jac_ompc, &
    !$omp& cp_coeff_ompc, &
    !$omp& vis_tag_ompc, &
    !$omp& wall_tag_ompc)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          uu = w_aux_ompc(i,j,k,2)
          vv = w_aux_ompc(i,j,k,3)
          ww = w_aux_ompc(i,j,k,4)
          tt = w_aux_ompc(i,j,k,6)
          mu = w_aux_ompc(i,j,k,7)
          lmaxi = lmax
          lmaxj = lmax
          if (j == 1) lmaxi = vis_tag_ompc(i)
          if (wall_tag_ompc(i) < 1) lmaxj = min(j,lmax)
          if (calorically_perfect==1) then
            cploc = cp_coeff_ompc(0)
          else
            cploc = 0._rkind
            do ll=indx_cp_l,indx_cp_r
              cploc = cploc+cp_coeff_ompc(ll)*(tt/t0)**ll
            enddo
          endif
          ucsi = 0._rkind
          vcsi = 0._rkind
          wcsi = 0._rkind
          tcsi = 0._rkind
          mucsi = 0._rkind
          ueta = 0._rkind
          veta = 0._rkind
          weta = 0._rkind
          teta = 0._rkind
          mueta = 0._rkind
          uzit = 0._rkind
          vzit = 0._rkind
          wzit = 0._rkind
          tzit = 0._rkind
          muzit = 0._rkind
          dg1 = 0._rkind
          dg2 = 0._rkind
          dg12csi = 0._rkind
          dg12eta = 0._rkind
          do l=1,lmax
            cli = coeff_deriv1_ompc(l,lmaxi)
            clj = coeff_deriv1_ompc(l,lmaxj)
            clk = coeff_deriv1_ompc(l,lmax )
            ucsi = ucsi +cli*(w_aux_ompc(i+l,j,k,2)-w_aux_ompc(i-l,j,k,2))
            vcsi = vcsi +cli*(w_aux_ompc(i+l,j,k,3)-w_aux_ompc(i-l,j,k,3))
            wcsi = wcsi +cli*(w_aux_ompc(i+l,j,k,4)-w_aux_ompc(i-l,j,k,4))
            tcsi = tcsi +cli*(w_aux_ompc(i+l,j,k,6)-w_aux_ompc(i-l,j,k,6))
            mucsi = mucsi+cli*(w_aux_ompc(i+l,j,k,7)-w_aux_ompc(i-l,j,k,7))
            ueta = ueta +clj*(w_aux_ompc(i,j+l,k,2)-w_aux_ompc(i,j-l,k,2))
            veta = veta +clj*(w_aux_ompc(i,j+l,k,3)-w_aux_ompc(i,j-l,k,3))
            weta = weta +clj*(w_aux_ompc(i,j+l,k,4)-w_aux_ompc(i,j-l,k,4))
            teta = teta +clj*(w_aux_ompc(i,j+l,k,6)-w_aux_ompc(i,j-l,k,6))
            mueta = mueta+clj*(w_aux_ompc(i,j+l,k,7)-w_aux_ompc(i,j-l,k,7))
            uzit = uzit +clk*(w_aux_ompc(i,j,k+l,2)-w_aux_ompc(i,j,k-l,2))
            vzit = vzit +clk*(w_aux_ompc(i,j,k+l,3)-w_aux_ompc(i,j,k-l,3))
            wzit = wzit +clk*(w_aux_ompc(i,j,k+l,4)-w_aux_ompc(i,j,k-l,4))
            tzit = tzit +clk*(w_aux_ompc(i,j,k+l,6)-w_aux_ompc(i,j,k-l,6))
            muzit = muzit+clk*(w_aux_ompc(i,j,k+l,7)-w_aux_ompc(i,j,k-l,7))
            dg1 = dg1 +cli*(g1_ompc (i+l,j)-g1_ompc (i-l,j))
            dg12csi=dg12csi+cli*(g12_ompc(i+l,j)-g12_ompc(i-l,j))
            dg2 = dg2 +clj*(g2_ompc (i,j+l)-g2_ompc (i,j-l))
            dg12eta=dg12eta+clj*(g12_ompc(i,j+l)-g12_ompc(i,j-l))
          enddo
          ux = ucsi *dcsidxc2_ompc(i,j) + ueta *detadxc2_ompc(i,j)
          vx = vcsi *dcsidxc2_ompc(i,j) + veta *detadxc2_ompc(i,j)
          wx = wcsi *dcsidxc2_ompc(i,j) + weta *detadxc2_ompc(i,j)
          tx = tcsi *dcsidxc2_ompc(i,j) + teta *detadxc2_ompc(i,j)
          mux = mucsi*dcsidxc2_ompc(i,j) + mueta*detadxc2_ompc(i,j)
          uy = ucsi *dcsidyc2_ompc(i,j) + ueta *detadyc2_ompc(i,j)
          vy = vcsi *dcsidyc2_ompc(i,j) + veta *detadyc2_ompc(i,j)
          wy = wcsi *dcsidyc2_ompc(i,j) + weta *detadyc2_ompc(i,j)
          ty = tcsi *dcsidyc2_ompc(i,j) + teta *detadyc2_ompc(i,j)
          muy = mucsi*dcsidyc2_ompc(i,j) + mueta*detadyc2_ompc(i,j)
          uz = uzit *dzitdz_ompc(k)
          vz = vzit *dzitdz_ompc(k)
          wz = wzit *dzitdz_ompc(k)
          tz = tzit *dzitdz_ompc(k)
          muz = muzit*dzitdz_ompc(k)
          if (j==1) then
            wallprop_ompc(i,k,2) = mu*ueta
            wallprop_ompc(i,k,3) = mu*weta
            wallprop_ompc(i,k,4) = mu*teta*cploc/prandtl
          endif
          ulapcsi = coeff_deriv2_ompc(0,lmaxi)*uu
          ulapeta = coeff_deriv2_ompc(0,lmaxj)*uu
          ulapzit = coeff_deriv2_ompc(0,lmax)*uu
          vlapcsi = coeff_deriv2_ompc(0,lmaxi)*vv
          vlapeta = coeff_deriv2_ompc(0,lmaxj)*vv
          vlapzit = coeff_deriv2_ompc(0,lmax)*vv
          wlapcsi = coeff_deriv2_ompc(0,lmaxi)*ww
          wlapeta = coeff_deriv2_ompc(0,lmaxj)*ww
          wlapzit = coeff_deriv2_ompc(0,lmax)*ww
          tlapcsi = coeff_deriv2_ompc(0,lmaxi)*tt
          tlapeta = coeff_deriv2_ompc(0,lmaxj)*tt
          tlapzit = coeff_deriv2_ompc(0,lmax)*tt
          do l=1,lmax
            clapi = coeff_deriv2_ompc(l,lmaxi)
            clapj = coeff_deriv2_ompc(l,lmaxj)
            clapk = coeff_deriv2_ompc(l,lmax )
            ulapcsi = ulapcsi + clapi*(w_aux_ompc(i+l,j,k,2)+w_aux_ompc(i-l,j,k,2))
            ulapeta = ulapeta + clapj*(w_aux_ompc(i,j+l,k,2)+w_aux_ompc(i,j-l,k,2))
            ulapzit = ulapzit + clapk*(w_aux_ompc(i,j,k+l,2)+w_aux_ompc(i,j,k-l,2))
            vlapcsi = vlapcsi + clapi*(w_aux_ompc(i+l,j,k,3)+w_aux_ompc(i-l,j,k,3))
            vlapeta = vlapeta + clapj*(w_aux_ompc(i,j+l,k,3)+w_aux_ompc(i,j-l,k,3))
            vlapzit = vlapzit + clapk*(w_aux_ompc(i,j,k+l,3)+w_aux_ompc(i,j,k-l,3))
            wlapcsi = wlapcsi + clapi*(w_aux_ompc(i+l,j,k,4)+w_aux_ompc(i-l,j,k,4))
            wlapeta = wlapeta + clapj*(w_aux_ompc(i,j+l,k,4)+w_aux_ompc(i,j-l,k,4))
            wlapzit = wlapzit + clapk*(w_aux_ompc(i,j,k+l,4)+w_aux_ompc(i,j,k-l,4))
            tlapcsi = tlapcsi + clapi*(w_aux_ompc(i+l,j,k,6)+w_aux_ompc(i-l,j,k,6))
            tlapeta = tlapeta + clapj*(w_aux_ompc(i,j+l,k,6)+w_aux_ompc(i,j-l,k,6))
            tlapzit = tlapzit + clapk*(w_aux_ompc(i,j,k+l,6)+w_aux_ompc(i,j,k-l,6))
          enddo
          ulap = g1_ompc(i,j)*ulapcsi+ucsi*dg1+g2_ompc(i,j)*ulapeta+ueta*dg2
          vlap = g1_ompc(i,j)*vlapcsi+vcsi*dg1+g2_ompc(i,j)*vlapeta+veta*dg2
          wlap = g1_ompc(i,j)*wlapcsi+wcsi*dg1+g2_ompc(i,j)*wlapeta+weta*dg2
          tlap = g1_ompc(i,j)*tlapcsi+tcsi*dg1+g2_ompc(i,j)*tlapeta+teta*dg2
          ulap = ulap * jac_ompc(i,j)
          vlap = vlap * jac_ompc(i,j)
          wlap = wlap * jac_ompc(i,j)
          tlap = tlap * jac_ompc(i,j)
          ulapz = ulapzit*dzitdzs_ompc(k)+uz*dzitdz2_ompc(k)
          vlapz = vlapzit*dzitdzs_ompc(k)+vz*dzitdz2_ompc(k)
          wlapz = wlapzit*dzitdzs_ompc(k)+wz*dzitdz2_ompc(k)
          tlapz = tlapzit*dzitdzs_ompc(k)+tz*dzitdz2_ompc(k)
          ulap = ulap+ulapz
          vlap = vlap+vlapz
          wlap = wlap+wlapz
          tlap = tlap+tlapz
          div = ux+vy+wz
          div3l = div/3._rkind
          w_aux_ompc(i,j,k,10) = div3l
          omegax = wy-vz
          omegay = uz-wx
          omegaz = vx-uy
          omod2 = omegax*omegax+omegay*omegay+omegaz*omegaz
          w_aux_ompc(i,j,k,9) = sqrt(omod2)
          w_aux_ompc(i,j,k,8) = max(-div/sqrt(omod2+div**2+(u0/l0)**2),0._rkind)
          if( j == 1) then
            if(iblock == ite_rank_x .and. i == ite_l) w_aux_ompc(i,j,k,8) = teshk
            if(iblock == itu_rank_x .and. i == itu_l) w_aux_ompc(i,j,k,8) = teshk
          endif
          sig11 = 2._rkind*(ux-div3l)
          sig12 = uy+vx
          sig13 = uz+wx
          sig22 = 2._rkind*(vy-div3l)
          sig23 = vz+wy
          sig33 = 2._rkind*(wz-div3l)
          sigx = mux*sig11 + muy*sig12 + muz*sig13 + mu*ulap
          sigy = mux*sig12 + muy*sig22 + muz*sig23 + mu*vlap
          sigz = mux*sig13 + muy*sig23 + muz*sig33 + mu*wlap
          sigqt = (mux*tx+muy*ty+muz*tz+mu*tlap)*cploc/prandtl
          sigah = (sig11*ux+sig12*uy+sig13*uz+sig12*vx+sig22*vy+sig23*vz+sig13*wx+sig23*wy+sig33*wz)*mu
          sigq = sigx*uu+sigy*vv+sigz*ww+sigah+sigqt
          fl_ompc(i,j,k,2) = fl_ompc(i,j,k,2) - sigx
          fl_ompc(i,j,k,3) = fl_ompc(i,j,k,3) - sigy
          fl_ompc(i,j,k,4) = fl_ompc(i,j,k,4) - sigz
          fl_ompc(i,j,k,5) = fl_ompc(i,j,k,5) - sigq
        enddo
      enddo
    enddo


  endsubroutine visflx_c2_subroutine




  subroutine visflx_nosensor_subroutine(nx,ny,nz,ng,visc_order,prandtl,t0,indx_cp_l,indx_cp_r,&
  &cp_coeff_ompc,calorically_perfect,u0,l0,w_ompc,w_aux_ompc,fl_ompc,coeff_deriv1_ompc,&
  &coeff_deriv2_ompc,dcsidx_ompc,detady_ompc,dzitdz_ompc,dcsidxs_ompc,detadys_ompc,dzitdzs_ompc,&
  &dcsidx2_ompc,detady2_ompc,dzitdz2_ompc,wallprop_ompc)
    integer, intent(in) :: nx, ny, nz, ng, visc_order, calorically_perfect
    integer, intent(in) :: indx_cp_l, indx_cp_r
    real(rkind), intent(in) :: prandtl, u0, l0, t0
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(in) :: w_ompc
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(inout) :: w_aux_ompc
    real(rkind), dimension(1-ng:,1-ng:, 2:), intent(inout) :: wallprop_ompc
    real(rkind), dimension(1:,1:,1:, 1:), intent(inout) :: fl_ompc
    real(rkind), dimension(1:,1:), intent(in) :: coeff_deriv1_ompc
    real(rkind), dimension(0:,1:), intent(in) :: coeff_deriv2_ompc
    real(rkind), dimension(1:), intent(in) :: dcsidx_ompc, detady_ompc, dzitdz_ompc
    real(rkind), dimension(1:), intent(in) :: dcsidxs_ompc, detadys_ompc, dzitdzs_ompc
    real(rkind), dimension(1:), intent(in) :: dcsidx2_ompc, detady2_ompc, dzitdz2_ompc
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in) :: cp_coeff_ompc
    integer :: i,j,k,l,ll,iercuda
    real(rkind) :: ccl,clapl
    real(rkind) :: sig11,sig12,sig13
    real(rkind) :: sig22,sig23
    real(rkind) :: sig33
    real(rkind) :: uu,vv,ww,tt,mu
    real(rkind) :: ux,uy,uz
    real(rkind) :: vx,vy,vz
    real(rkind) :: wx,wy,wz
    real(rkind) :: tx,ty,tz
    real(rkind) :: mux, muy, muz
    real(rkind) :: ulap,ulapx,ulapy,ulapz
    real(rkind) :: vlap,vlapx,vlapy,vlapz
    real(rkind) :: wlap,wlapx,wlapy,wlapz
    real(rkind) :: tlap,tlapx,tlapy,tlapz
    real(rkind) :: sigq,sigx,sigy,sigz,sigqt,sigah
    real(rkind) :: div, div3l, omegax, omegay, omegaz, omod2
    real(rkind) :: cploc
    integer :: lmax
    lmax = visc_order/2
    !$omp parallel do default(firstprivate) shared(w_ompc, &
    !$omp& w_aux_ompc, &
    !$omp& wallprop_ompc, &
    !$omp& fl_ompc, &
    !$omp& coeff_deriv1_ompc, &
    !$omp& coeff_deriv2_ompc, &
    !$omp& dcsidx_ompc, &
    !$omp& detady_ompc, &
    !$omp& dzitdz_ompc, &
    !$omp& dcsidxs_ompc, &
    !$omp& detadys_ompc, &
    !$omp& dzitdzs_ompc, &
    !$omp& dcsidx2_ompc, &
    !$omp& detady2_ompc, &
    !$omp& dzitdz2_ompc, &
    !$omp& cp_coeff_ompc)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          uu = w_aux_ompc(i,j,k,2)
          vv = w_aux_ompc(i,j,k,3)
          ww = w_aux_ompc(i,j,k,4)
          tt = w_aux_ompc(i,j,k,6)
          mu = w_aux_ompc(i,j,k,7)
          ux = 0._rkind
          vx = 0._rkind
          wx = 0._rkind
          tx = 0._rkind
          mux = 0._rkind
          uy = 0._rkind
          vy = 0._rkind
          wy = 0._rkind
          ty = 0._rkind
          muy = 0._rkind
          uz = 0._rkind
          vz = 0._rkind
          wz = 0._rkind
          tz = 0._rkind
          muz = 0._rkind
          ulapx = coeff_deriv2_ompc(0,lmax)*uu
          ulapy = ulapx
          ulapz = ulapx
          vlapx = coeff_deriv2_ompc(0,lmax)*vv
          vlapy = vlapx
          vlapz = vlapx
          wlapx = coeff_deriv2_ompc(0,lmax)*ww
          wlapy = wlapx
          wlapz = wlapx
          tlapx = coeff_deriv2_ompc(0,lmax)*tt
          tlapy = tlapx
          tlapz = tlapx
          do l=1,lmax
            ccl = coeff_deriv1_ompc(l,lmax)
            clapl = coeff_deriv2_ompc(l,lmax)
            ux = ux+ccl*(w_aux_ompc(i+l,j,k,2)-w_aux_ompc(i-l,j,k,2))
            ulapx = ulapx + clapl*(w_aux_ompc(i+l,j,k,2)+w_aux_ompc(i-l,j,k,2))
            vx = vx+ccl*(w_aux_ompc(i+l,j,k,3)-w_aux_ompc(i-l,j,k,3))
            vlapx = vlapx + clapl*(w_aux_ompc(i+l,j,k,3)+w_aux_ompc(i-l,j,k,3))
            wx = wx+ccl*(w_aux_ompc(i+l,j,k,4)-w_aux_ompc(i-l,j,k,4))
            wlapx = wlapx + clapl*(w_aux_ompc(i+l,j,k,4)+w_aux_ompc(i-l,j,k,4))
            tx = tx+ccl*(w_aux_ompc(i+l,j,k,6)-w_aux_ompc(i-l,j,k,6))
            tlapx = tlapx + clapl*(w_aux_ompc(i+l,j,k,6)+w_aux_ompc(i-l,j,k,6))
            mux = mux+ccl*(w_aux_ompc(i+l,j,k,7)-w_aux_ompc(i-l,j,k,7))
          enddo
          do l=1,lmax
            ccl = coeff_deriv1_ompc(l,lmax)
            clapl = coeff_deriv2_ompc(l,lmax)
            uy = uy+ccl*(w_aux_ompc(i,j+l,k,2)-w_aux_ompc(i,j-l,k,2))
            ulapy = ulapy + clapl*(w_aux_ompc(i,j+l,k,2)+w_aux_ompc(i,j-l,k,2))
            vy = vy+ccl*(w_aux_ompc(i,j+l,k,3)-w_aux_ompc(i,j-l,k,3))
            vlapy = vlapy + clapl*(w_aux_ompc(i,j+l,k,3)+w_aux_ompc(i,j-l,k,3))
            wy = wy+ccl*(w_aux_ompc(i,j+l,k,4)-w_aux_ompc(i,j-l,k,4))
            wlapy = wlapy + clapl*(w_aux_ompc(i,j+l,k,4)+w_aux_ompc(i,j-l,k,4))
            ty = ty+ccl*(w_aux_ompc(i,j+l,k,6)-w_aux_ompc(i,j-l,k,6))
            tlapy = tlapy + clapl*(w_aux_ompc(i,j+l,k,6)+w_aux_ompc(i,j-l,k,6))
            muy = muy+ccl*(w_aux_ompc(i,j+l,k,7)-w_aux_ompc(i,j-l,k,7))
          enddo
          do l=1,lmax
            ccl = coeff_deriv1_ompc(l,lmax)
            clapl = coeff_deriv2_ompc(l,lmax)
            uz = uz+ccl*(w_aux_ompc(i,j,k+l,2)-w_aux_ompc(i,j,k-l,2))
            ulapz = ulapz + clapl*(w_aux_ompc(i,j,k+l,2)+w_aux_ompc(i,j,k-l,2))
            vz = vz+ccl*(w_aux_ompc(i,j,k+l,3)-w_aux_ompc(i,j,k-l,3))
            vlapz = vlapz + clapl*(w_aux_ompc(i,j,k+l,3)+w_aux_ompc(i,j,k-l,3))
            wz = wz+ccl*(w_aux_ompc(i,j,k+l,4)-w_aux_ompc(i,j,k-l,4))
            wlapz = wlapz + clapl*(w_aux_ompc(i,j,k+l,4)+w_aux_ompc(i,j,k-l,4))
            tz = tz+ccl*(w_aux_ompc(i,j,k+l,6)-w_aux_ompc(i,j,k-l,6))
            tlapz = tlapz + clapl*(w_aux_ompc(i,j,k+l,6)+w_aux_ompc(i,j,k-l,6))
            muz = muz+ccl*(w_aux_ompc(i,j,k+l,7)-w_aux_ompc(i,j,k-l,7))
          enddo
          ux = ux*dcsidx_ompc(i)
          vx = vx*dcsidx_ompc(i)
          wx = wx*dcsidx_ompc(i)
          tx = tx*dcsidx_ompc(i)
          mux = mux*dcsidx_ompc(i)
          uy = uy*detady_ompc(j)
          vy = vy*detady_ompc(j)
          wy = wy*detady_ompc(j)
          ty = ty*detady_ompc(j)
          muy = muy*detady_ompc(j)
          uz = uz*dzitdz_ompc(k)
          vz = vz*dzitdz_ompc(k)
          wz = wz*dzitdz_ompc(k)
          tz = tz*dzitdz_ompc(k)
          muz = muz*dzitdz_ompc(k)
          ulapx = ulapx*dcsidxs_ompc(i)+ux*dcsidx2_ompc(i)
          vlapx = vlapx*dcsidxs_ompc(i)+vx*dcsidx2_ompc(i)
          wlapx = wlapx*dcsidxs_ompc(i)+wx*dcsidx2_ompc(i)
          tlapx = tlapx*dcsidxs_ompc(i)+tx*dcsidx2_ompc(i)
          ulapy = ulapy*detadys_ompc(j)+uy*detady2_ompc(j)
          vlapy = vlapy*detadys_ompc(j)+vy*detady2_ompc(j)
          wlapy = wlapy*detadys_ompc(j)+wy*detady2_ompc(j)
          tlapy = tlapy*detadys_ompc(j)+ty*detady2_ompc(j)
          ulapz = ulapz*dzitdzs_ompc(k)+uz*dzitdz2_ompc(k)
          vlapz = vlapz*dzitdzs_ompc(k)+vz*dzitdz2_ompc(k)
          wlapz = wlapz*dzitdzs_ompc(k)+wz*dzitdz2_ompc(k)
          tlapz = tlapz*dzitdzs_ompc(k)+tz*dzitdz2_ompc(k)
          ulap = ulapx+ulapy+ulapz
          vlap = vlapx+vlapy+vlapz
          wlap = wlapx+wlapy+wlapz
          tlap = tlapx+tlapy+tlapz
          div = ux+vy+wz
          div3l = div/3._rkind
          w_aux_ompc(i,j,k,10) = div3l
          omegax = wy-vz
          omegay = uz-wx
          omegaz = vx-uy
          omod2 = omegax*omegax+omegay*omegay+omegaz*omegaz
          w_aux_ompc(i,j,k,9) = sqrt(omod2)
          if (calorically_perfect==1) then
            cploc = cp_coeff_ompc(0)
          else
            cploc = 0._rkind
            do ll=indx_cp_l,indx_cp_r
              cploc = cploc+cp_coeff_ompc(ll)*(tt/t0)**ll
            enddo
          endif
          sig11 = 2._rkind*(ux-div3l)
          sig12 = uy+vx
          sig13 = uz+wx
          sig22 = 2._rkind*(vy-div3l)
          sig23 = vz+wy
          sig33 = 2._rkind*(wz-div3l)
          sigx = mux*sig11 + muy*sig12 + muz*sig13 + mu*ulap
          sigy = mux*sig12 + muy*sig22 + muz*sig23 + mu*vlap
          sigz = mux*sig13 + muy*sig23 + muz*sig33 + mu*wlap
          sigqt = (mux*tx+muy*ty+muz*tz+mu*tlap)*cploc/prandtl
          sigah = (sig11*ux+sig12*uy+sig13*uz+sig12*vx+sig22*vy+sig23*vz+sig13*wx+sig23*wy+sig33*wz)*mu
          sigq = sigx*uu+sigy*vv+sigz*ww+sigah+sigqt
          fl_ompc(i,j,k,2) = fl_ompc(i,j,k,2) - sigx
          fl_ompc(i,j,k,3) = fl_ompc(i,j,k,3) - sigy
          fl_ompc(i,j,k,4) = fl_ompc(i,j,k,4) - sigz
          fl_ompc(i,j,k,5) = fl_ompc(i,j,k,5) - sigq
          if (j==1) then
            wallprop_ompc(i,k,2) = mu*uy
            wallprop_ompc(i,k,3) = mu*wy
            wallprop_ompc(i,k,4) = mu*ty*cploc/prandtl
          endif
        enddo
      enddo
    enddo


  endsubroutine visflx_nosensor_subroutine



  subroutine visflx_nosensor_c2_subroutine(nx,ny,nz,ng,visc_order,prandtl,t0,indx_cp_l,indx_cp_r,&
  &cp_coeff_ompc,calorically_perfect,iblock,ite_rank_x,itu_rank_x,ite_l,itu_l,u0,l0,w_ompc,w_aux_ompc,&
  &fl_ompc,coeff_deriv1_ompc,coeff_deriv2_ompc,dcsidxc2_ompc,detadyc2_ompc,detadxc2_ompc,dcsidyc2_ompc,&
  &dzitdz_ompc,dzitdzs_ompc,dzitdz2_ompc,g1_ompc,g2_ompc,g12_ompc,jac_ompc,wallprop_ompc,vis_tag_ompc,&
  &wall_tag_ompc,ortho)
    integer, intent(in) :: nx, ny, nz, ng, visc_order, calorically_perfect
    integer, intent(in) :: indx_cp_l, indx_cp_r, ortho
    integer, intent(in) :: iblock, ite_rank_x, itu_rank_x, ite_l, itu_l
    real(rkind), intent(in) :: prandtl, u0, l0, t0
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(in) :: w_ompc
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(inout) :: w_aux_ompc
    real(rkind), dimension(1-ng:,1-ng:, 2:), intent(inout) :: wallprop_ompc
    real(rkind), dimension(1:,1:,1:, 1:), intent(inout) :: fl_ompc
    real(rkind), dimension(1:,1:), intent(in) :: coeff_deriv1_ompc
    real(rkind), dimension(0:,1:), intent(in) :: coeff_deriv2_ompc
    real(rkind), dimension(1:), intent(in) :: dzitdz_ompc, dzitdzs_ompc, dzitdz2_ompc
    real(rkind), dimension(1-ng:,1-ng:), intent(in) :: dcsidxc2_ompc, detadyc2_ompc
    real(rkind), dimension(1-ng:,1-ng:), intent(in) :: detadxc2_ompc, dcsidyc2_ompc
    real(rkind), dimension(1-ng:,1-ng:), intent(in) :: g1_ompc, g2_ompc, g12_ompc, jac_ompc
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in) :: cp_coeff_ompc
    integer, dimension(1-ng:), intent(in) :: vis_tag_ompc, wall_tag_ompc
    integer :: i,j,k,l,ll,iercuda
    real(rkind) :: ccl,clapl,cli,clj,clk,clapi,clapj,clapk
    real(rkind) :: sig11,sig12,sig13
    real(rkind) :: sig22,sig23
    real(rkind) :: sig33
    real(rkind) :: uu,vv,ww,tt,mu
    real(rkind) :: ux,uy,uz
    real(rkind) :: vx,vy,vz
    real(rkind) :: wx,wy,wz
    real(rkind) :: tx,ty,tz
    real(rkind) :: mux, muy, muz
    real(rkind) :: ucsi,ueta,uzit,vcsi,veta,vzit,wcsi,weta,wzit
    real(rkind) :: tcsi,teta,tzit,mucsi,mueta,muzit
    real(rkind) :: ulapcsi,ulapeta,ulapzit,vlapcsi,vlapeta,vlapzit,wlapcsi,wlapeta,wlapzit
    real(rkind) :: ulapcsieta,vlapcsieta,wlapcsieta,tlapcsieta
    real(rkind) :: dg1,dg2,dg12csi,dg12eta
    real(rkind) :: tlapcsi,tlapeta,tlapzit
    real(rkind) :: ulap,ulapx,ulapy,ulapz
    real(rkind) :: vlap,vlapx,vlapy,vlapz
    real(rkind) :: wlap,wlapx,wlapy,wlapz
    real(rkind) :: tlap,tlapx,tlapy,tlapz
    real(rkind) :: sigq,sigx,sigy,sigz,sigqt,sigah
    real(rkind) :: div, div3l, omegax, omegay, omegaz, omod2
    real(rkind) :: cploc
    integer :: lmax,lmaxi,lmaxj
    lmax = visc_order/2
    !$omp parallel do default(firstprivate) shared(w_ompc, &
    !$omp& w_aux_ompc, &
    !$omp& wallprop_ompc, &
    !$omp& fl_ompc, &
    !$omp& coeff_deriv1_ompc, &
    !$omp& coeff_deriv2_ompc, &
    !$omp& dzitdz_ompc, &
    !$omp& dzitdzs_ompc, &
    !$omp& dzitdz2_ompc, &
    !$omp& dcsidxc2_ompc, &
    !$omp& detadyc2_ompc, &
    !$omp& detadxc2_ompc, &
    !$omp& dcsidyc2_ompc, &
    !$omp& g1_ompc, &
    !$omp& g2_ompc, &
    !$omp& g12_ompc, &
    !$omp& jac_ompc, &
    !$omp& cp_coeff_ompc, &
    !$omp& vis_tag_ompc, &
    !$omp& wall_tag_ompc)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          uu = w_aux_ompc(i,j,k,2)
          vv = w_aux_ompc(i,j,k,3)
          ww = w_aux_ompc(i,j,k,4)
          tt = w_aux_ompc(i,j,k,6)
          mu = w_aux_ompc(i,j,k,7)
          lmaxi = lmax
          lmaxj = lmax
          if (j == 1) lmaxi = vis_tag_ompc(i)
          if (wall_tag_ompc(i) < 1) lmaxj = min(j,lmax)
          if (calorically_perfect==1) then
            cploc = cp_coeff_ompc(0)
          else
            cploc = 0._rkind
            do ll=indx_cp_l,indx_cp_r
              cploc = cploc+cp_coeff_ompc(ll)*(tt/t0)**ll
            enddo
          endif
          ucsi = 0._rkind
          vcsi = 0._rkind
          wcsi = 0._rkind
          tcsi = 0._rkind
          mucsi = 0._rkind
          ueta = 0._rkind
          veta = 0._rkind
          weta = 0._rkind
          teta = 0._rkind
          mueta = 0._rkind
          uzit = 0._rkind
          vzit = 0._rkind
          wzit = 0._rkind
          tzit = 0._rkind
          muzit = 0._rkind
          dg1 = 0._rkind
          dg2 = 0._rkind
          dg12csi = 0._rkind
          dg12eta = 0._rkind
          do l=1,lmax
            cli = coeff_deriv1_ompc(l,lmaxi)
            clj = coeff_deriv1_ompc(l,lmaxj)
            clk = coeff_deriv1_ompc(l,lmax )
            ucsi = ucsi +cli*(w_aux_ompc(i+l,j,k,2)-w_aux_ompc(i-l,j,k,2))
            vcsi = vcsi +cli*(w_aux_ompc(i+l,j,k,3)-w_aux_ompc(i-l,j,k,3))
            wcsi = wcsi +cli*(w_aux_ompc(i+l,j,k,4)-w_aux_ompc(i-l,j,k,4))
            tcsi = tcsi +cli*(w_aux_ompc(i+l,j,k,6)-w_aux_ompc(i-l,j,k,6))
            mucsi = mucsi+cli*(w_aux_ompc(i+l,j,k,7)-w_aux_ompc(i-l,j,k,7))
            ueta = ueta +clj*(w_aux_ompc(i,j+l,k,2)-w_aux_ompc(i,j-l,k,2))
            veta = veta +clj*(w_aux_ompc(i,j+l,k,3)-w_aux_ompc(i,j-l,k,3))
            weta = weta +clj*(w_aux_ompc(i,j+l,k,4)-w_aux_ompc(i,j-l,k,4))
            teta = teta +clj*(w_aux_ompc(i,j+l,k,6)-w_aux_ompc(i,j-l,k,6))
            mueta = mueta+clj*(w_aux_ompc(i,j+l,k,7)-w_aux_ompc(i,j-l,k,7))
            uzit = uzit +clk*(w_aux_ompc(i,j,k+l,2)-w_aux_ompc(i,j,k-l,2))
            vzit = vzit +clk*(w_aux_ompc(i,j,k+l,3)-w_aux_ompc(i,j,k-l,3))
            wzit = wzit +clk*(w_aux_ompc(i,j,k+l,4)-w_aux_ompc(i,j,k-l,4))
            tzit = tzit +clk*(w_aux_ompc(i,j,k+l,6)-w_aux_ompc(i,j,k-l,6))
            muzit = muzit+clk*(w_aux_ompc(i,j,k+l,7)-w_aux_ompc(i,j,k-l,7))
            dg1 = dg1 +cli*(g1_ompc (i+l,j)-g1_ompc (i-l,j))
            dg12csi=dg12csi+cli*(g12_ompc(i+l,j)-g12_ompc(i-l,j))
            dg2 = dg2 +clj*(g2_ompc (i,j+l)-g2_ompc (i,j-l))
            dg12eta=dg12eta+clj*(g12_ompc(i,j+l)-g12_ompc(i,j-l))
          enddo
          ux = ucsi *dcsidxc2_ompc(i,j) + ueta *detadxc2_ompc(i,j)
          vx = vcsi *dcsidxc2_ompc(i,j) + veta *detadxc2_ompc(i,j)
          wx = wcsi *dcsidxc2_ompc(i,j) + weta *detadxc2_ompc(i,j)
          tx = tcsi *dcsidxc2_ompc(i,j) + teta *detadxc2_ompc(i,j)
          mux = mucsi*dcsidxc2_ompc(i,j) + mueta*detadxc2_ompc(i,j)
          uy = ucsi *dcsidyc2_ompc(i,j) + ueta *detadyc2_ompc(i,j)
          vy = vcsi *dcsidyc2_ompc(i,j) + veta *detadyc2_ompc(i,j)
          wy = wcsi *dcsidyc2_ompc(i,j) + weta *detadyc2_ompc(i,j)
          ty = tcsi *dcsidyc2_ompc(i,j) + teta *detadyc2_ompc(i,j)
          muy = mucsi*dcsidyc2_ompc(i,j) + mueta*detadyc2_ompc(i,j)
          uz = uzit *dzitdz_ompc(k)
          vz = vzit *dzitdz_ompc(k)
          wz = wzit *dzitdz_ompc(k)
          tz = tzit *dzitdz_ompc(k)
          muz = muzit*dzitdz_ompc(k)
          if (j==1) then
            wallprop_ompc(i,k,3) = mu*weta
            wallprop_ompc(i,k,4) = mu*teta*cploc/prandtl
          endif
          ulapcsi = coeff_deriv2_ompc(0,lmaxi)*uu
          ulapeta = coeff_deriv2_ompc(0,lmaxj)*uu
          ulapzit = coeff_deriv2_ompc(0,lmax)*uu
          vlapcsi = coeff_deriv2_ompc(0,lmaxi)*vv
          vlapeta = coeff_deriv2_ompc(0,lmaxj)*vv
          vlapzit = coeff_deriv2_ompc(0,lmax)*vv
          wlapcsi = coeff_deriv2_ompc(0,lmaxi)*ww
          wlapeta = coeff_deriv2_ompc(0,lmaxj)*ww
          wlapzit = coeff_deriv2_ompc(0,lmax)*ww
          tlapcsi = coeff_deriv2_ompc(0,lmaxi)*tt
          tlapeta = coeff_deriv2_ompc(0,lmaxj)*tt
          tlapzit = coeff_deriv2_ompc(0,lmax)*tt
          do l=1,lmax
            clapi = coeff_deriv2_ompc(l,lmaxi)
            clapj = coeff_deriv2_ompc(l,lmaxj)
            clapk = coeff_deriv2_ompc(l,lmax )
            ulapcsi = ulapcsi + clapi*(w_aux_ompc(i+l,j,k,2)+w_aux_ompc(i-l,j,k,2))
            ulapeta = ulapeta + clapj*(w_aux_ompc(i,j+l,k,2)+w_aux_ompc(i,j-l,k,2))
            ulapzit = ulapzit + clapk*(w_aux_ompc(i,j,k+l,2)+w_aux_ompc(i,j,k-l,2))
            vlapcsi = vlapcsi + clapi*(w_aux_ompc(i+l,j,k,3)+w_aux_ompc(i-l,j,k,3))
            vlapeta = vlapeta + clapj*(w_aux_ompc(i,j+l,k,3)+w_aux_ompc(i,j-l,k,3))
            vlapzit = vlapzit + clapk*(w_aux_ompc(i,j,k+l,3)+w_aux_ompc(i,j,k-l,3))
            wlapcsi = wlapcsi + clapi*(w_aux_ompc(i+l,j,k,4)+w_aux_ompc(i-l,j,k,4))
            wlapeta = wlapeta + clapj*(w_aux_ompc(i,j+l,k,4)+w_aux_ompc(i,j-l,k,4))
            wlapzit = wlapzit + clapk*(w_aux_ompc(i,j,k+l,4)+w_aux_ompc(i,j,k-l,4))
            tlapcsi = tlapcsi + clapi*(w_aux_ompc(i+l,j,k,6)+w_aux_ompc(i-l,j,k,6))
            tlapeta = tlapeta + clapj*(w_aux_ompc(i,j+l,k,6)+w_aux_ompc(i,j-l,k,6))
            tlapzit = tlapzit + clapk*(w_aux_ompc(i,j,k+l,6)+w_aux_ompc(i,j,k-l,6))
          enddo
          if(ortho == 1) then
            ulap = g1_ompc(i,j)*ulapcsi+ucsi*dg1+g2_ompc(i,j)*ulapeta+ueta*dg2
            vlap = g1_ompc(i,j)*vlapcsi+vcsi*dg1+g2_ompc(i,j)*vlapeta+veta*dg2
            wlap = g1_ompc(i,j)*wlapcsi+wcsi*dg1+g2_ompc(i,j)*wlapeta+weta*dg2
            tlap = g1_ompc(i,j)*tlapcsi+tcsi*dg1+g2_ompc(i,j)*tlapeta+teta*dg2
          else
            ulapcsieta = 0.25_rkind*(w_aux_ompc(i+1,j+1,k,2)-w_aux_ompc(i-1,j+1,k,2)-w_aux_ompc(i+1,&
            &j-1,k,2)+w_aux_ompc(i-1,j-1,k,2))
            vlapcsieta = 0.25_rkind*(w_aux_ompc(i+1,j+1,k,3)-w_aux_ompc(i-1,j+1,k,3)-w_aux_ompc(i+1,&
            &j-1,k,3)+w_aux_ompc(i-1,j-1,k,3))
            wlapcsieta = 0.25_rkind*(w_aux_ompc(i+1,j+1,k,4)-w_aux_ompc(i-1,j+1,k,4)-w_aux_ompc(i+1,&
            &j-1,k,4)+w_aux_ompc(i-1,j-1,k,4))
            tlapcsieta = 0.25_rkind*(w_aux_ompc(i+1,j+1,k,6)-w_aux_ompc(i-1,j+1,k,6)-w_aux_ompc(i+1,&
            &j-1,k,6)+w_aux_ompc(i-1,j-1,k,6))
            if (iblock==ite_rank_x.and.i==ite_l-1.and.j==1) then
              ulapcsieta = 0.5_rkind*(w_aux_ompc(i,j+1,k,2)-w_aux_ompc(i-1,j+1,k,2)-w_aux_ompc(i,j-1,k,2)+w_aux_ompc(i-1,j-1,k,2))
              vlapcsieta = 0.5_rkind*(w_aux_ompc(i,j+1,k,3)-w_aux_ompc(i-1,j+1,k,3)-w_aux_ompc(i,j-1,k,3)+w_aux_ompc(i-1,j-1,k,3))
              wlapcsieta = 0.5_rkind*(w_aux_ompc(i,j+1,k,4)-w_aux_ompc(i-1,j+1,k,4)-w_aux_ompc(i,j-1,k,4)+w_aux_ompc(i-1,j-1,k,4))
              tlapcsieta = 0.5_rkind*(w_aux_ompc(i,j+1,k,6)-w_aux_ompc(i-1,j+1,k,6)-w_aux_ompc(i,j-1,k,6)+w_aux_ompc(i-1,j-1,k,6))
            endif
            if (iblock==itu_rank_x.and.i==itu_l+1.and.j==1) then
              ulapcsieta = 0.5_rkind*(w_aux_ompc(i+1,j+1,k,2)-w_aux_ompc(i,j+1,k,2)-w_aux_ompc(i+1,j-1,k,2)+w_aux_ompc(i,j-1,k,2))
              vlapcsieta = 0.5_rkind*(w_aux_ompc(i+1,j+1,k,3)-w_aux_ompc(i,j+1,k,3)-w_aux_ompc(i+1,j-1,k,3)+w_aux_ompc(i,j-1,k,3))
              wlapcsieta = 0.5_rkind*(w_aux_ompc(i+1,j+1,k,4)-w_aux_ompc(i,j+1,k,4)-w_aux_ompc(i+1,j-1,k,4)+w_aux_ompc(i,j-1,k,4))
              tlapcsieta = 0.5_rkind*(w_aux_ompc(i+1,j+1,k,6)-w_aux_ompc(i,j+1,k,6)-w_aux_ompc(i+1,j-1,k,6)+w_aux_ompc(i,j-1,k,6))
            endif
            ulap = g1_ompc(i,j)*ulapcsi+ucsi*dg1+2._rkind*g12_ompc(i,j)*ulapcsieta+dg12csi*ueta+&
            &dg12eta*ucsi+g2_ompc(i,j)*ulapeta+ueta*dg2
            vlap = g1_ompc(i,j)*vlapcsi+vcsi*dg1+2._rkind*g12_ompc(i,j)*vlapcsieta+dg12csi*veta+&
            &dg12eta*vcsi+g2_ompc(i,j)*vlapeta+veta*dg2
            wlap = g1_ompc(i,j)*wlapcsi+wcsi*dg1+2._rkind*g12_ompc(i,j)*wlapcsieta+dg12csi*weta+&
            &dg12eta*wcsi+g2_ompc(i,j)*wlapeta+weta*dg2
            tlap = g1_ompc(i,j)*tlapcsi+tcsi*dg1+2._rkind*g12_ompc(i,j)*tlapcsieta+dg12csi*teta+&
            &dg12eta*tcsi+g2_ompc(i,j)*tlapeta+teta*dg2
          endif
          ulap = ulap * jac_ompc(i,j)
          vlap = vlap * jac_ompc(i,j)
          wlap = wlap * jac_ompc(i,j)
          tlap = tlap * jac_ompc(i,j)
          ulapz = ulapzit*dzitdzs_ompc(k)+uz*dzitdz2_ompc(k)
          vlapz = vlapzit*dzitdzs_ompc(k)+vz*dzitdz2_ompc(k)
          wlapz = wlapzit*dzitdzs_ompc(k)+wz*dzitdz2_ompc(k)
          tlapz = tlapzit*dzitdzs_ompc(k)+tz*dzitdz2_ompc(k)
          ulap = ulap+ulapz
          vlap = vlap+vlapz
          wlap = wlap+wlapz
          tlap = tlap+tlapz
          div = ux+vy+wz
          div3l = div/3._rkind
          w_aux_ompc(i,j,k,10) = div3l
          omegax = wy-vz
          omegay = uz-wx
          omegaz = vx-uy
          omod2 = omegax*omegax+omegay*omegay+omegaz*omegaz
          w_aux_ompc(i,j,k,9) = sqrt(omod2)
          sig11 = 2._rkind*(ux-div3l)
          sig12 = uy+vx
          sig13 = uz+wx
          sig22 = 2._rkind*(vy-div3l)
          sig23 = vz+wy
          sig33 = 2._rkind*(wz-div3l)
          sigx = mux*sig11 + muy*sig12 + muz*sig13 + mu*ulap
          sigy = mux*sig12 + muy*sig22 + muz*sig23 + mu*vlap
          sigz = mux*sig13 + muy*sig23 + muz*sig33 + mu*wlap
          sigqt = (mux*tx+muy*ty+muz*tz+mu*tlap)*cploc/prandtl
          sigah = (sig11*ux+sig12*uy+sig13*uz+sig12*vx+sig22*vy+sig23*vz+sig13*wx+sig23*wy+sig33*wz)*mu
          sigq = sigx*uu+sigy*vv+sigz*ww+sigah+sigqt
          fl_ompc(i,j,k,2) = fl_ompc(i,j,k,2) - sigx
          fl_ompc(i,j,k,3) = fl_ompc(i,j,k,3) - sigy
          fl_ompc(i,j,k,4) = fl_ompc(i,j,k,4) - sigz
          fl_ompc(i,j,k,5) = fl_ompc(i,j,k,5) - sigq
        enddo
      enddo
    enddo


  endsubroutine visflx_nosensor_c2_subroutine



  subroutine sponge_subroutine(nx,ny,nz,ng,nv,w_ompc,wfar_ompc,fln_ompc,f_sponge_ompc,j_sponge)
    integer, intent(in) :: nx,ny,nz,ng,nv,j_sponge
    real(rkind), dimension(1-ng:, 1-ng:, 1-ng:, 1:), intent(in) :: w_ompc
    real(rkind), dimension(1:,1:,1:,1:), intent(inout) :: fln_ompc
    real(rkind), dimension(nx,nv), intent(in) :: wfar_ompc
    real(rkind), dimension(ny), intent(in) :: f_sponge_ompc
    integer :: i,j,k,m,iercuda

    !$omp parallel do default(firstprivate) shared(w_ompc, &
    !$omp& fln_ompc, &
    !$omp& wfar_ompc, &
    !$omp& f_sponge_ompc)
    do k = 1,nz
      do j = j_sponge,ny
        do i = 1,nx
          do m=1,nv
            fln_ompc(i,j,k,m) = fln_ompc(i,j,k,m) - f_sponge_ompc(j)*(w_ompc(i,j,k,m) - wfar_ompc(i,m))
          enddo
        enddo
      enddo
    enddo


  endsubroutine sponge_subroutine



  subroutine limiter_subroutine(nx,ny,nz,ng,w_ompc,w_aux_ompc,iblock,kblock,indx_cp_l,indx_cp_r,&
  &cv_coeff_ompc,calorically_perfect,t0,tol_iter_nr,rho_lim,tem_lim,rho_lim_rescale,tem_lim_rescale)
    integer, intent(in) :: nx,ny,nz,ng,iblock,kblock,calorically_perfect,indx_cp_l,indx_cp_r
    real(rkind), intent(in) :: t0, tol_iter_nr
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(inout) :: w_ompc, w_aux_ompc
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in) :: cv_coeff_ompc
    real(rkind), intent(in) :: rho_lim, tem_lim, rho_lim_rescale, tem_lim_rescale
    real(rkind) :: rho, tem, eei, uu, vv, ww, rhoe, qq, ee
    integer :: i,j,k,iercuda
    real(rkind) :: n_limited_rho, n_limited_tem
    n_limited_rho = 0._rkind
    n_limited_tem = 0._rkind
    !$omp parallel do default(firstprivate) shared(w_ompc, &
    !$omp& w_aux_ompc, &
    !$omp& cv_coeff_ompc) &
    !$omp& reduction(+:n_limited_rho,n_limited_tem)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          rho = w_ompc(i,j,k,1)
          uu = w_ompc(i,j,k,2)/rho
          vv = w_ompc(i,j,k,3)/rho
          ww = w_ompc(i,j,k,4)/rho
          rhoe = w_ompc(i,j,k,5)
          qq = 0.5_rkind*(uu*uu+vv*vv+ww*ww)
          ee = rhoe/rho-qq
          tem = get_temperature_from_e_dev(ee,w_aux_ompc(i,j,k,6),t0,cv_coeff_ompc,indx_cp_l,&
          &indx_cp_r,calorically_perfect,tol_iter_nr)
          if((rho < rho_lim).or.(tem < tem_lim)) then
            if(rho < rho_lim) then
              n_limited_rho = n_limited_rho + 1._rkind
            endif
            if(tem < tem_lim) then
              n_limited_tem = n_limited_tem + 1._rkind
            endif
          endif
        enddo
      enddo
    enddo
    if(n_limited_rho > 0.) print*,'warning! n_limited_rho :',n_limited_rho
    if(n_limited_tem > 0.) print*,'warning! n_limited_tem :',n_limited_tem

    if(n_limited_rho > 0. .or. n_limited_tem > 0.) then
      !$omp parallel do default(firstprivate) shared(w_ompc, &
      !$omp& w_aux_ompc, &
      !$omp& cv_coeff_ompc)
      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            rho = w_ompc(i,j,k,1)
            uu = w_ompc(i,j,k,2)/rho
            vv = w_ompc(i,j,k,3)/rho
            ww = w_ompc(i,j,k,4)/rho
            rhoe = w_ompc(i,j,k,5)
            qq = 0.5_rkind*(uu*uu+vv*vv+ww*ww)
            ee = rhoe/rho-qq
            tem = get_temperature_from_e_dev(ee,w_aux_ompc(i,j,k,6),t0,cv_coeff_ompc,indx_cp_l,&
            &indx_cp_r,calorically_perfect,tol_iter_nr)
            if((rho < rho_lim).or.(tem < tem_lim)) then
              if(rho < rho_lim) then
                rho = rho_lim*rho_lim_rescale
                print*,'limiter rhofix: ',iblock,kblock,i,j,k
              endif
              if(tem < tem_lim) then
                tem = tem_lim*tem_lim_rescale
                print*,'limiter temfix: ',iblock,kblock,i,j,k
              endif
              w_ompc(i,j,k,1) = rho
              w_ompc(i,j,k,2) = rho*uu
              w_ompc(i,j,k,3) = rho*vv
              w_ompc(i,j,k,4) = rho*ww
              eei = get_e_from_temperature_dev(tem, t0, indx_cp_l, indx_cp_r, cv_coeff_ompc, calorically_perfect)
              w_ompc(i,j,k,5) = rho*eei + 0.5_rkind*(w_ompc(i,j,k,2)**2+w_ompc(i,j,k,3)**2+w_ompc(i,j,k,4)**2)/rho
            endif
          enddo
        enddo
      enddo
    endif

  endsubroutine limiter_subroutine



  subroutine filter_subroutine(nx,ny,nz,ng,w_ompc,w_aux_ompc,indx_cp_l,indx_cp_r,cv_coeff_ompc,&
  &calorically_perfect,t0,tol_iter_nr,coeff_filter_ompc,jfilter,wall_tag_ompc)
    integer, intent(in) :: nx,ny,nz,ng,calorically_perfect,indx_cp_l,indx_cp_r,jfilter
    integer, dimension(1-ng:nx+ng), intent(in) :: wall_tag_ompc
    real(rkind), intent(in) :: t0, tol_iter_nr
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(inout) :: w_ompc, w_aux_ompc
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in) :: cv_coeff_ompc
    real(rkind), dimension(0:4), intent(in) :: coeff_filter_ompc
    real(rkind) :: cfilt,rho,eei,uu,vv,ww,rhoe,qq,ee,rho_p,rho_m,rho_f,tem_p,tem_m,tem_f
    integer :: i,j,k,l,iercuda

    !$omp parallel do default(firstprivate) shared(w_ompc, &
    !$omp& w_aux_ompc, &
    !$omp& cv_coeff_ompc, &
    !$omp& coeff_filter_ompc, &
    !$omp& wall_tag_ompc)
    do k = 1,nz
      do j = 1+4,ny-4
        do i = 1,nx
          if (wall_tag_ompc(i) < 1 .and. j < jfilter) then
          else
            rho_f = 0._rkind
            tem_f = 0._rkind
            do l=0,4
              rho_p = w_ompc(i+l,j,k,1)
              uu = w_ompc(i+l,j,k,2)/rho_p
              vv = w_ompc(i+l,j,k,3)/rho_p
              ww = w_ompc(i+l,j,k,4)/rho_p
              rhoe = w_ompc(i+l,j,k,5)
              qq = 0.5_rkind*(uu*uu+vv*vv+ww*ww)
              ee = rhoe/rho_p-qq
              tem_p = get_temperature_from_e_dev(ee,w_aux_ompc(i+l,j,k,6),t0,cv_coeff_ompc,&
              &indx_cp_l,indx_cp_r,calorically_perfect,tol_iter_nr)
              rho_m = w_ompc(i-l,j,k,1)
              uu = w_ompc(i-l,j,k,2)/rho_m
              vv = w_ompc(i-l,j,k,3)/rho_m
              ww = w_ompc(i-l,j,k,4)/rho_m
              rhoe = w_ompc(i-l,j,k,5)
              qq = 0.5_rkind*(uu*uu+vv*vv+ww*ww)
              ee = rhoe/rho_m-qq
              tem_m = get_temperature_from_e_dev(ee,w_aux_ompc(i-l,j,k,6),t0,cv_coeff_ompc,&
              &indx_cp_l,indx_cp_r,calorically_perfect,tol_iter_nr)
              cfilt = coeff_filter_ompc(l)
              rho_f = rho_f + cfilt*(rho_p+rho_m)
              tem_f = tem_f + cfilt*(tem_p+tem_m)
            enddo
            uu = w_ompc(i,j,k,2)/w_ompc(i,j,k,1)
            vv = w_ompc(i,j,k,3)/w_ompc(i,j,k,1)
            ww = w_ompc(i,j,k,4)/w_ompc(i,j,k,1)
            w_ompc(i,j,k,1) = rho_f
            w_ompc(i,j,k,2) = rho_f*uu
            w_ompc(i,j,k,3) = rho_f*vv
            w_ompc(i,j,k,4) = rho_f*ww
            eei = get_e_from_temperature_dev(tem_f,t0,indx_cp_l,indx_cp_r,cv_coeff_ompc,calorically_perfect)
            w_ompc(i,j,k,5) = rho_f*eei + 0.5_rkind*(w_ompc(i,j,k,2)**2+w_ompc(i,j,k,3)**2+w_ompc(i,j,k,4)**2)/rho_f
            w_aux_ompc(i,j,k,6) = tem_f
          endif
        enddo
      enddo
    enddo

    !$omp parallel do default(firstprivate) shared(w_ompc, &
    !$omp& w_aux_ompc, &
    !$omp& cv_coeff_ompc, &
    !$omp& coeff_filter_ompc, &
    !$omp& wall_tag_ompc)
    do k = 1,nz
      do j = 1+4,ny-4
        do i = 1,nx
          if (wall_tag_ompc(i) < 1 .and. j < jfilter) then
          else
            rho_f = 0._rkind
            tem_f = 0._rkind
            do l=0,4
              rho_p = w_ompc(i,j+l,k,1)
              uu = w_ompc(i,j+l,k,2)/rho_p
              vv = w_ompc(i,j+l,k,3)/rho_p
              ww = w_ompc(i,j+l,k,4)/rho_p
              rhoe = w_ompc(i,j+l,k,5)
              qq = 0.5_rkind*(uu*uu+vv*vv+ww*ww)
              ee = rhoe/rho_p-qq
              tem_p = get_temperature_from_e_dev(ee,w_aux_ompc(i+l,j,k,6),t0,cv_coeff_ompc,&
              &indx_cp_l,indx_cp_r,calorically_perfect,tol_iter_nr)
              rho_m = w_ompc(i,j-l,k,1)
              uu = w_ompc(i,j-l,k,2)/rho_m
              vv = w_ompc(i,j-l,k,3)/rho_m
              ww = w_ompc(i,j-l,k,4)/rho_m
              rhoe = w_ompc(i,j-l,k,5)
              qq = 0.5_rkind*(uu*uu+vv*vv+ww*ww)
              ee = rhoe/rho_m-qq
              tem_m = get_temperature_from_e_dev(ee,w_aux_ompc(i-l,j,k,6),t0,cv_coeff_ompc,&
              &indx_cp_l,indx_cp_r,calorically_perfect,tol_iter_nr)
              cfilt = coeff_filter_ompc(l)
              rho_f = rho_f + cfilt*(rho_p+rho_m)
              tem_f = tem_f + cfilt*(tem_p+tem_m)
            enddo
            uu = w_ompc(i,j,k,2)/w_ompc(i,j,k,1)
            vv = w_ompc(i,j,k,3)/w_ompc(i,j,k,1)
            ww = w_ompc(i,j,k,4)/w_ompc(i,j,k,1)
            w_ompc(i,j,k,1) = rho_f
            w_ompc(i,j,k,2) = rho_f*uu
            w_ompc(i,j,k,3) = rho_f*vv
            w_ompc(i,j,k,4) = rho_f*ww
            eei = get_e_from_temperature_dev(tem_f,t0,indx_cp_l,indx_cp_r,cv_coeff_ompc,calorically_perfect)
            w_ompc(i,j,k,5) = rho_f*eei + 0.5_rkind*(w_ompc(i,j,k,2)**2+w_ompc(i,j,k,3)**2+w_ompc(i,j,k,4)**2)/rho_f
            w_aux_ompc(i,j,k,6) = tem_f
          endif
        enddo
      enddo
    enddo


  endsubroutine filter_subroutine



  subroutine visflx_reduced_ord2_subroutine(nx,ny,nz,ng,prandtl,t0,indx_cp_l,indx_cp_r,&
  &cp_coeff_ompc,calorically_perfect,u0,l0,w_ompc,w_aux_ompc,fl_ompc,x_ompc,y_ompc,z_ompc,&
  &wallprop_ompc,update_sensor)
    integer, intent(in) :: nx, ny, nz, ng, calorically_perfect
    integer, intent(in) :: indx_cp_l, indx_cp_r
    real(rkind), intent(in) :: prandtl, u0, l0, t0
    integer, intent(in) :: update_sensor
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(in) :: w_ompc
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(inout) :: w_aux_ompc
    real(rkind), dimension(1-ng:,1-ng:, 2:), intent(inout) :: wallprop_ompc
    real(rkind), dimension(1:,1:,1:, 1:), intent(inout) :: fl_ompc
    real(rkind), dimension(1-ng:), intent(in) :: x_ompc, y_ompc, z_ompc
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in) :: cp_coeff_ompc
    integer :: i,j,k,l,ll,iercuda
    real(rkind) :: dxl,dyl,dzl
    real(rkind) :: sig11,sig12,sig13
    real(rkind) :: sig21,sig22,sig23
    real(rkind) :: sig31,sig32,sig33
    real(rkind) :: uu,vv,ww,tt,mu
    real(rkind) :: ux,uy,uz
    real(rkind) :: vx,vy,vz
    real(rkind) :: wx,wy,wz
    real(rkind) :: tx,ty,tz
    real(rkind) :: mux, muy, muz
    real(rkind) :: sigq,sigx,sigy,sigz,sigqt,sigah
    real(rkind) :: div, div3l, omegax, omegay, omegaz, omod2
    real(rkind) :: cploc

    !$omp parallel do default(firstprivate) shared(w_ompc, &
    !$omp& w_aux_ompc, &
    !$omp& wallprop_ompc, &
    !$omp& fl_ompc, &
    !$omp& x_ompc, &
    !$omp& y_ompc, &
    !$omp& z_ompc, &
    !$omp& cp_coeff_ompc)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          uu = w_aux_ompc(i,j,k,2)
          vv = w_aux_ompc(i,j,k,3)
          ww = w_aux_ompc(i,j,k,4)
          tt = w_aux_ompc(i,j,k,6)
          mu = w_aux_ompc(i,j,k,7)
          if (calorically_perfect==1) then
            cploc = cp_coeff_ompc(0)
          else
            cploc = 0._rkind
            do ll=indx_cp_l,indx_cp_r
              cploc = cploc+cp_coeff_ompc(ll)*(tt/t0)**ll
            enddo
          endif
          dxl = 1._rkind/(x_ompc(i+1)-x_ompc(i-1))
          dyl = 1._rkind/(y_ompc(j+1)-y_ompc(j-1))
          dzl = 1._rkind/(z_ompc(k+1)-z_ompc(k-1))
          ux = dxl*(w_aux_ompc(i+1,j,k,2)-w_aux_ompc(i-1,j,k,2))
          vx = dxl*(w_aux_ompc(i+1,j,k,3)-w_aux_ompc(i-1,j,k,3))
          wx = dxl*(w_aux_ompc(i+1,j,k,4)-w_aux_ompc(i-1,j,k,4))
          mux = dxl*(w_aux_ompc(i+1,j,k,7)-w_aux_ompc(i-1,j,k,7))
          uy = dyl*(w_aux_ompc(i,j+1,k,2)-w_aux_ompc(i,j-1,k,2))
          vy = dyl*(w_aux_ompc(i,j+1,k,3)-w_aux_ompc(i,j-1,k,3))
          wy = dyl*(w_aux_ompc(i,j+1,k,4)-w_aux_ompc(i,j-1,k,4))
          ty = dyl*(w_aux_ompc(i,j+1,k,6)-w_aux_ompc(i,j-1,k,6))
          muy = dyl*(w_aux_ompc(i,j+1,k,7)-w_aux_ompc(i,j-1,k,7))
          uz = dzl*(w_aux_ompc(i,j,k+1,2)-w_aux_ompc(i,j,k-1,2))
          vz = dzl*(w_aux_ompc(i,j,k+1,3)-w_aux_ompc(i,j,k-1,3))
          wz = dzl*(w_aux_ompc(i,j,k+1,4)-w_aux_ompc(i,j,k-1,4))
          muz = dzl*(w_aux_ompc(i,j,k+1,7)-w_aux_ompc(i,j,k-1,7))
          if (j==1) then
            wallprop_ompc(i,k,2) = mu*uy
            wallprop_ompc(i,k,3) = mu*wy
            wallprop_ompc(i,k,4) = mu*ty*cploc/prandtl
          endif
          div = ux+vy+wz
          div3l = div/3._rkind
          w_aux_ompc(i,j,k,10) = div3l
          omegax = wy-vz
          omegay = uz-wx
          omegaz = vx-uy
          omod2 = omegax*omegax+omegay*omegay+omegaz*omegaz
          w_aux_ompc(i,j,k,9) = sqrt(omod2)
          if (update_sensor == 1) then
            w_aux_ompc(i,j,k,8) = (max(-div/sqrt(omod2+div**2+(u0/l0)**2),0._rkind))**2
          endif
          sig11 = ux-2._rkind*div3l
          sig12 = vx
          sig13 = wx
          sig21 = uy
          sig22 = vy-2._rkind*div3l
          sig23 = wy
          sig31 = uz
          sig32 = vz
          sig33 = wz-2._rkind*div3l
          sigx = mux*sig11 + muy*sig12 + muz*sig13
          sigy = mux*sig21 + muy*sig22 + muz*sig23
          sigz = mux*sig31 + muy*sig32 + muz*sig33
          sigah = (sig11*ux+sig12*uy+sig13*uz+sig21*vx+sig22*vy+sig23*vz+sig31*wx+sig32*wy+sig33*wz)*mu
          sigq = sigx*uu+sigy*vv+sigz*ww+sigah
          fl_ompc(i,j,k,2) = fl_ompc(i,j,k,2) - sigx
          fl_ompc(i,j,k,3) = fl_ompc(i,j,k,3) - sigy
          fl_ompc(i,j,k,4) = fl_ompc(i,j,k,4) - sigz
          fl_ompc(i,j,k,5) = fl_ompc(i,j,k,5) - sigq
        enddo
      enddo
    enddo


  endsubroutine visflx_reduced_ord2_subroutine




  subroutine sensor_subroutine(nx,ny,nz,ng,u0,l0,w_aux_ompc)
    integer, intent(in) :: nx, ny, nz, ng
    real(rkind), intent(in) :: u0, l0
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(inout) :: w_aux_ompc
    integer :: i,j,k,iercuda
    real(rkind) :: div, omod2

    !$omp parallel do default(firstprivate) shared(w_aux_ompc)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          omod2 = w_aux_ompc(i,j,k,9)**2
          div = 3._rkind*w_aux_ompc(i,j,k,10)
          w_aux_ompc(i,j,k,8) = (max(-div/sqrt(omod2+div**2+(u0/l0)**2),0._rkind))**2
        enddo
      enddo
    enddo


  endsubroutine sensor_subroutine



  subroutine sensor_c2_subroutine(nx,ny,nz,ng,u0,l0_ducros,w_aux_ompc,iblock,ite_rank_x,itu_rank_x,&
  &ite_l,itu_l,teshk,theta_ij_ompc,theta_threshold,jweno)
    integer, intent(in) :: nx, ny, nz, ng, iblock, ite_rank_x, itu_rank_x, ite_l, itu_l, jweno
    real(rkind), intent(in) :: u0, l0_ducros, teshk, theta_threshold
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(inout) :: w_aux_ompc
    real(rkind), dimension(1:,1:), intent(inout) :: theta_ij_ompc
    integer :: i,j,k,iercuda
    real(rkind) :: div, omod2

    !$omp parallel do default(firstprivate) shared(w_aux_ompc, &
    !$omp& theta_ij_ompc)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          omod2 = w_aux_ompc(i,j,k,9)**2
          div = 3._rkind*w_aux_ompc(i,j,k,10)
          w_aux_ompc(i,j,k,8) = max(-div/sqrt(omod2+div**2+(u0/l0_ducros)**2),0._rkind)
          if( j == 1) then
            if(iblock == ite_rank_x .and. i == ite_l) w_aux_ompc(i,j,k,8) = teshk
            if(iblock == itu_rank_x .and. i == itu_l) w_aux_ompc(i,j,k,8) = teshk
          endif
          if(theta_ij_ompc(i,j) > theta_threshold .or. j > jweno) then
            w_aux_ompc(i,j,k,8) = 1000._rkind
          endif
        enddo
      enddo
    enddo


  endsubroutine sensor_c2_subroutine



  subroutine visflx_x_subroutine(nx,ny,nz,nv,ng,prandtl,t0,indx_cp_l,indx_cp_r,cp_coeff_ompc,&
  &calorically_perfect,x_ompc,w_aux_ompc,fl_ompc,fhat_ompc)
    integer, intent(in) :: nx, ny, nz, nv, ng, calorically_perfect
    integer, intent(in) :: indx_cp_l, indx_cp_r
    real(rkind), intent(in) :: prandtl, t0
    real(rkind), dimension(1:nx,1:ny,1:nz,1:nv), intent(inout) :: fl_ompc
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in) :: w_aux_ompc
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv), intent(out) :: fhat_ompc
    real(rkind), dimension(1-ng:), intent(in) :: x_ompc
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in) :: cp_coeff_ompc
    integer :: i,j,k,iv,ll,iercuda
    real(rkind) :: uu,vv,ww,tt,mu,qq
    real(rkind) :: uup,vvp,wwp,ttp,mup,qqp
    real(rkind) :: sigq,sigx,sigy,sigz,sigq_tt,sigq_qq
    real(rkind) :: dxhl,fl2o,fl3o,fl4o,fl5o
    real(rkind) :: ttf,muf
    real(rkind) :: cploc

    !$omp parallel do default(firstprivate) shared(fl_ompc, &
    !$omp& w_aux_ompc, &
    !$omp& fhat_ompc, &
    !$omp& x_ompc, &
    !$omp& cp_coeff_ompc)
    do k = 1,nz
      do j = 1,ny
        do i = 0,nx
          uu = w_aux_ompc(i ,j,k,2)
          uup = w_aux_ompc(i+1,j,k,2)
          vv = w_aux_ompc(i ,j,k,3)
          vvp = w_aux_ompc(i+1,j,k,3)
          ww = w_aux_ompc(i ,j,k,4)
          wwp = w_aux_ompc(i+1,j,k,4)
          tt = w_aux_ompc(i ,j,k,6)
          ttp = w_aux_ompc(i+1,j,k,6)
          mu = w_aux_ompc(i ,j,k,7)
          mup = w_aux_ompc(i+1,j,k,7)
          qq = 0.5_rkind*(uu*uu+vv*vv+ww*ww)
          qqp = 0.5_rkind*(uup*uup+vvp*vvp+wwp*wwp)
          if (calorically_perfect==1) then
            cploc = cp_coeff_ompc(0)
          else
            ttf = 0.5_rkind*(tt+ttp)
            cploc = 0._rkind
            do ll=indx_cp_l,indx_cp_r
              cploc = cploc+cp_coeff_ompc(ll)*(ttf/t0)**ll
            enddo
          endif
          sigx = uup-uu
          sigy = vvp-vv
          sigz = wwp-ww
          sigq_tt = ttp-tt
          sigq_qq = qqp-qq
          muf = mu+mup
          muf = 0.5_rkind*muf/(x_ompc(i+1)-x_ompc(i))
          sigx = sigx*muf
          sigy = sigy*muf
          sigz = sigz*muf
          sigq = (sigq_tt*cploc/prandtl+sigq_qq)*muf
          fhat_ompc(i,j,k,2) = - sigx
          fhat_ompc(i,j,k,3) = - sigy
          fhat_ompc(i,j,k,4) = - sigz
          fhat_ompc(i,j,k,5) = - sigq
        enddo
      enddo
    enddo

    !$omp parallel do default(firstprivate) shared(fl_ompc, &
    !$omp& w_aux_ompc, &
    !$omp& fhat_ompc, &
    !$omp& x_ompc, &
    !$omp& cp_coeff_ompc)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          dxhl = 2._rkind/(x_ompc(i+1)-x_ompc(i-1))
          do iv=2,nv
            fl_ompc(i,j,k,iv) = fl_ompc(i,j,k,iv) + dxhl*(fhat_ompc(i,j,k,iv)-fhat_ompc(i-1,j,k,iv))
          enddo
        enddo
      enddo
    enddo


  endsubroutine visflx_x_subroutine



  subroutine visflx_y_subroutine(nx,ny,nz,nv,ng,prandtl,t0,indx_cp_l,indx_cp_r,cp_coeff_ompc,&
  &calorically_perfect,y_ompc,w_aux_ompc,fl_ompc)
    integer, intent(in) :: nx, ny, nz, nv, ng, calorically_perfect
    integer, intent(in) :: indx_cp_l, indx_cp_r
    real(rkind), intent(in) :: prandtl, t0
    real(rkind), dimension(1:nx,1:ny,1:nz,1:nv), intent(inout) :: fl_ompc
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in) :: w_aux_ompc
    real(rkind), dimension(1-ng:), intent(in) :: y_ompc
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in) :: cp_coeff_ompc
    integer :: i,j,k,iv,ll,iercuda
    real(rkind) :: uu,vv,ww,tt,mu,qq
    real(rkind) :: uup,vvp,wwp,ttp,mup,qqp
    real(rkind) :: sigq,sigx,sigy,sigz,sigq_tt,sigq_qq
    real(rkind) :: dyhl,fl2o,fl3o,fl4o,fl5o
    real(rkind) :: ttf,muf
    real(rkind) :: cploc

    !$omp parallel do default(firstprivate) shared(fl_ompc, &
    !$omp& w_aux_ompc, &
    !$omp& y_ompc, &
    !$omp& cp_coeff_ompc)
    do k = 1,nz
      do i = 1,nx
        do j=0,ny
          uu = w_aux_ompc(i,j ,k,2)
          uup = w_aux_ompc(i,j+1,k,2)
          vv = w_aux_ompc(i,j ,k,3)
          vvp = w_aux_ompc(i,j+1,k,3)
          ww = w_aux_ompc(i,j ,k,4)
          wwp = w_aux_ompc(i,j+1,k,4)
          tt = w_aux_ompc(i,j ,k,6)
          ttp = w_aux_ompc(i,j+1,k,6)
          mu = w_aux_ompc(i,j ,k,7)
          mup = w_aux_ompc(i,j+1,k,7)
          qq = 0.5_rkind*(uu*uu+vv*vv+ww*ww)
          qqp = 0.5_rkind*(uup*uup+vvp*vvp+wwp*wwp)
          if (calorically_perfect==1) then
            cploc = cp_coeff_ompc(0)
          else
            ttf = 0.5_rkind*(tt+ttp)
            cploc = 0._rkind
            do ll=indx_cp_l,indx_cp_r
              cploc = cploc+cp_coeff_ompc(ll)*(ttf/t0)**ll
            enddo
          endif
          sigx = uup-uu
          sigy = vvp-vv
          sigz = wwp-ww
          sigq_tt = ttp-tt
          sigq_qq = qqp-qq
          muf = mu+mup
          muf = 0.5_rkind*muf/(y_ompc(j+1)-y_ompc(j))
          sigx = sigx*muf
          sigy = sigy*muf
          sigz = sigz*muf
          sigq = (sigq_tt*cploc/prandtl+sigq_qq)*muf
          if (j>0) then
            fl_ompc(i,j,k,2) = fl_ompc(i,j,k,2) + fl2o-sigx*dyhl
            fl_ompc(i,j,k,3) = fl_ompc(i,j,k,3) + fl3o-sigy*dyhl
            fl_ompc(i,j,k,4) = fl_ompc(i,j,k,4) + fl4o-sigz*dyhl
            fl_ompc(i,j,k,5) = fl_ompc(i,j,k,5) + fl5o-sigq*dyhl
          endif
          if (j<ny) then
            dyhl = 2._rkind/(y_ompc(j+2)-y_ompc(j))
            fl2o = sigx*dyhl
            fl3o = sigy*dyhl
            fl4o = sigz*dyhl
            fl5o = sigq*dyhl
          endif
        enddo
      enddo
    enddo


  endsubroutine visflx_y_subroutine


  subroutine visflx_z_subroutine(nx,ny,nz,nv,ng,prandtl,t0,indx_cp_l,indx_cp_r,cp_coeff_ompc,&
  &calorically_perfect,z_ompc,w_aux_ompc,fl_ompc)
    integer, intent(in) :: nx, ny, nz, nv, ng, calorically_perfect
    integer, intent(in) :: indx_cp_l, indx_cp_r
    real(rkind), intent(in) :: prandtl, t0
    real(rkind), dimension(1:nx,1:ny,1:nz,1:nv), intent(inout) :: fl_ompc
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in) :: w_aux_ompc
    real(rkind), dimension(1-ng:), intent(in) :: z_ompc
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in) :: cp_coeff_ompc
    integer :: i,j,k,iv,ll,iercuda
    real(rkind) :: uu,vv,ww,tt,mu,qq
    real(rkind) :: uup,vvp,wwp,ttp,mup,qqp
    real(rkind) :: sigq,sigx,sigy,sigz,sigq_tt,sigq_qq
    real(rkind) :: dzhl,fl2o,fl3o,fl4o,fl5o
    real(rkind) :: ttf,muf
    real(rkind) :: cploc

    !$omp parallel do default(firstprivate) shared(fl_ompc, &
    !$omp& w_aux_ompc, &
    !$omp& z_ompc, &
    !$omp& cp_coeff_ompc)
    do j = 1,ny
      do i = 1,nx
        do k=0,nz
          uu = w_aux_ompc(i,j,k ,2)
          uup = w_aux_ompc(i,j,k+1,2)
          vv = w_aux_ompc(i,j,k ,3)
          vvp = w_aux_ompc(i,j,k+1,3)
          ww = w_aux_ompc(i,j,k ,4)
          wwp = w_aux_ompc(i,j,k+1,4)
          tt = w_aux_ompc(i,j,k ,6)
          ttp = w_aux_ompc(i,j,k+1,6)
          mu = w_aux_ompc(i,j,k ,7)
          mup = w_aux_ompc(i,j,k+1,7)
          qq = 0.5_rkind*(uu*uu+vv*vv+ww*ww)
          qqp = 0.5_rkind*(uup*uup+vvp*vvp+wwp*wwp)
          if (calorically_perfect==1) then
            cploc = cp_coeff_ompc(0)
          else
            ttf = 0.5_rkind*(tt+ttp)
            cploc = 0._rkind
            do ll=indx_cp_l,indx_cp_r
              cploc = cploc+cp_coeff_ompc(ll)*(ttf/t0)**ll
            enddo
          endif
          sigx = uup-uu
          sigy = vvp-vv
          sigz = wwp-ww
          sigq_tt = ttp-tt
          sigq_qq = qqp-qq
          muf = mu+mup
          muf = 0.5_rkind*muf/(z_ompc(k+1)-z_ompc(k))
          sigx = sigx*muf
          sigy = sigy*muf
          sigz = sigz*muf
          sigq = (sigq_tt*cploc/prandtl+sigq_qq)*muf
          if (k>0) then
            fl_ompc(i,j,k,2) = fl_ompc(i,j,k,2) + fl2o-sigx*dzhl
            fl_ompc(i,j,k,3) = fl_ompc(i,j,k,3) + fl3o-sigy*dzhl
            fl_ompc(i,j,k,4) = fl_ompc(i,j,k,4) + fl4o-sigz*dzhl
            fl_ompc(i,j,k,5) = fl_ompc(i,j,k,5) + fl5o-sigq*dzhl
          endif
          if (k<nz) then
            dzhl = 2._rkind/(z_ompc(k+2)-z_ompc(k))
            fl2o = sigx*dzhl
            fl3o = sigy*dzhl
            fl4o = sigz*dzhl
            fl5o = sigq*dzhl
          endif
        enddo
      enddo
    enddo


  endsubroutine visflx_z_subroutine



  subroutine recyc_exchange_subroutine_1(irecyc,w_ompc,wbuf1s_ompc,nx,ny,nz,ng,nv)
    integer, intent(in) :: irecyc, nx, ny, nz, ng, nv
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in) :: w_ompc
    real(rkind), dimension(:,:,:,:), intent(inout) :: wbuf1s_ompc
    integer :: i,j,k,m, iercuda

    !$omp parallel do default(firstprivate) shared(w_ompc, &
    !$omp& wbuf1s_ompc)
    do k = 1,nz
      do j = 1,ny
        do i = 1,ng
          do m=1,nv
            wbuf1s_ompc(i,j,k,m) = w_ompc(irecyc+1-i,j,k,m)
          enddo
        enddo
      enddo
    enddo


  endsubroutine recyc_exchange_subroutine_1


  subroutine recyc_exchange_subroutine_2(n1_start_recv,n1_start_send,n1_end_recv,wrecyc_ompc,wbuf1r_ompc,nx,ny,nz,ng,nv)
    integer, intent(in) :: n1_start_recv, n1_start_send, n1_end_recv, nx, ny, nz, ng, nv
    real(rkind), dimension(:,:,:,:), intent(in) :: wbuf1r_ompc
    real(rkind), dimension(:,:,:,:), intent(inout) :: wrecyc_ompc
    integer :: i,j,k,m, iercuda

    !$omp parallel do default(firstprivate) shared(wbuf1r_ompc, &
    !$omp& wrecyc_ompc)
    do k = n1_start_recv,n1_end_recv
      do j = 1,ny
        do i = 1,ng
          do m=1,nv
            wrecyc_ompc(i,j,k,m) = wbuf1r_ompc(i,j,k-n1_start_recv+n1_start_send,m)
          enddo
        enddo
      enddo
    enddo


  endsubroutine recyc_exchange_subroutine_2


  subroutine recyc_exchange_subroutine_3(n2_start_recv,n2_start_send,n2_end_recv,wrecyc_ompc,wbuf2r_ompc,nx,ny,nz,ng,nv)
    integer, intent(in) :: n2_start_recv, n2_start_send, n2_end_recv, nx, ny, nz, ng, nv
    real(rkind), dimension(:,:,:,:), intent(in) :: wbuf2r_ompc
    real(rkind), dimension(:,:,:,:), intent(inout) :: wrecyc_ompc
    integer :: i,j,k,m, iercuda

    !$omp parallel do default(firstprivate) shared(wbuf2r_ompc, &
    !$omp& wrecyc_ompc)
    do k = n2_start_recv,n2_end_recv
      do j = 1,ny
        do i = 1,ng
          do m=1,nv
            wrecyc_ompc(i,j,k,m) = wbuf2r_ompc(i,j,k-n2_start_recv+n2_start_send,m)
          enddo
        enddo
      enddo
    enddo


  endsubroutine recyc_exchange_subroutine_3



  subroutine bcextr_sub_subroutine(ilat,nx,ny,nz,ng,p0,rgas0,w_ompc,indx_cp_l,indx_cp_r,cv_coeff_ompc,t0,calorically_perfect)
    integer, intent(in) :: ilat, nx, ny, nz, ng, indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind), intent(in) :: p0, t0, rgas0
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in) :: cv_coeff_ompc
    real(rkind), dimension(1-ng:, 1-ng:, 1-ng:, 1:), intent(inout) :: w_ompc
    real(rkind) :: rho,rhou,rhov,rhow,tt,ee
    integer :: i,j,k,l,m,iercuda
    if (ilat==1) then
      !$omp parallel do default(firstprivate) shared(cv_coeff_ompc, &
      !$omp& w_ompc)
      do k = 1,nz
        do j = 1,ny
          do l = 1,ng
            rho = w_ompc(1,j,k,1)
            rhou = w_ompc(1,j,k,2)
            rhov = w_ompc(1,j,k,3)
            rhow = w_ompc(1,j,k,4)
            tt = p0/rho/rgas0
            ee = get_e_from_temperature_dev(tt,t0,indx_cp_l,indx_cp_r,cv_coeff_ompc,calorically_perfect)
            w_ompc(1-l,j,k,1) = rho
            w_ompc(1-l,j,k,2) = rhou
            w_ompc(1-l,j,k,3) = rhov
            w_ompc(1-l,j,k,4) = rhow
            w_ompc(1-l,j,k,5) = rho*ee + 0.5_rkind*(rhou**2+rhov**2+rhow**2)/rho
          enddo
        enddo
      enddo
    elseif (ilat==2) then
      !$omp parallel do default(firstprivate) shared(cv_coeff_ompc, &
      !$omp& w_ompc)
      do k = 1,nz
        do j = 1,ny
          do l = 1,ng
            rho = w_ompc(nx,j,k,1)
            rhou = w_ompc(nx,j,k,2)
            rhov = w_ompc(nx,j,k,3)
            rhow = w_ompc(nx,j,k,4)
            tt = p0/rho/rgas0
            ee = get_e_from_temperature_dev(tt,t0,indx_cp_l,indx_cp_r,cv_coeff_ompc,calorically_perfect)
            w_ompc(nx+l,j,k,1) = rho
            w_ompc(nx+l,j,k,2) = rhou
            w_ompc(nx+l,j,k,3) = rhov
            w_ompc(nx+l,j,k,4) = rhow
            w_ompc(nx+l,j,k,5) = rho*ee + 0.5_rkind*(rhou**2+rhov**2+rhow**2)/rho
          enddo
        enddo
      enddo
    elseif (ilat==3) then
    elseif (ilat==4) then
      !$omp parallel do default(firstprivate) shared(cv_coeff_ompc, &
      !$omp& w_ompc)
      do k = 1,nz
        do i = 1,nx
          do l = 1,ng
            rho = w_ompc(i,ny,k,1)
            rhou = w_ompc(i,ny,k,2)
            rhov = w_ompc(i,ny,k,3)
            rhow = w_ompc(i,ny,k,4)
            tt = p0/rho/rgas0
            ee = get_e_from_temperature_dev(tt, t0, indx_cp_l, indx_cp_r, cv_coeff_ompc, calorically_perfect)
            w_ompc(i,ny+l,k,1) = rho
            w_ompc(i,ny+l,k,2) = rhou
            w_ompc(i,ny+l,k,3) = rhov
            w_ompc(i,ny+l,k,4) = rhow
            w_ompc(i,ny+l,k,5) = rho*ee + 0.5_rkind*(rhou**2+rhov**2+rhow**2)/rho
          enddo
        enddo
      enddo
    elseif (ilat==5) then
      !$omp parallel do default(firstprivate) shared(cv_coeff_ompc, &
      !$omp& w_ompc)
      do l = 1,ng
        do j = 1,ny
          do i = 1,nx
            rho = w_ompc(i,j,1,1)
            rhou = w_ompc(i,j,1,2)
            rhov = w_ompc(i,j,1,3)
            rhow = w_ompc(i,j,1,4)
            tt = p0/rho/rgas0
            ee = get_e_from_temperature_dev(tt, t0, indx_cp_l, indx_cp_r, cv_coeff_ompc, calorically_perfect)
            w_ompc(i,j,1-l,1) = rho
            w_ompc(i,j,1-l,2) = rhou
            w_ompc(i,j,1-l,3) = rhov
            w_ompc(i,j,1-l,4) = rhow
            w_ompc(i,j,1-l,5) = rho*ee + 0.5_rkind*(rhou**2+rhov**2+rhow**2)/rho
          enddo
        enddo
      enddo
    elseif (ilat==6) then
      !$omp parallel do default(firstprivate) shared(cv_coeff_ompc, &
      !$omp& w_ompc)
      do l = 1,ng
        do j = 1,ny
          do i = 1,nx
            rho = w_ompc(i,j,nz,1)
            rhou = w_ompc(i,j,nz,2)
            rhov = w_ompc(i,j,nz,3)
            rhow = w_ompc(i,j,nz,4)
            tt = p0/rho/rgas0
            ee = get_e_from_temperature_dev(tt, t0, indx_cp_l, indx_cp_r, cv_coeff_ompc, calorically_perfect)
            w_ompc(i,j,nz+l,1) = rho
            w_ompc(i,j,nz+l,2) = rhou
            w_ompc(i,j,nz+l,3) = rhov
            w_ompc(i,j,nz+l,4) = rhow
            w_ompc(i,j,nz+l,5) = rho*ee + 0.5_rkind*(rhou**2+rhov**2+rhow**2)/rho
          enddo
        enddo
      enddo
    endif

  endsubroutine bcextr_sub_subroutine




  subroutine bc_nr_lat_x_kernel(start_or_end,nr_type,nx,ny,nz,ng,nv,w_aux_ompc,w_ompc,fl_ompc,&
  &dcsidx_ompc,indx_cp_l,indx_cp_r,cp_coeff_ompc,winf_ompc,calorically_perfect,rgas0,t0)
    integer, intent(in), value :: start_or_end, nr_type, nx, ny, nz, ng, nv, indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind), intent(in), value :: rgas0,t0
    real(rkind), dimension(3) :: c_one
    real(rkind), dimension(5) :: dw_dn, dwc_dn, ev, dw_dn_outer, dwc_dn_outer
    real(rkind), dimension(5,5) :: el, er
    real(rkind) :: w_target
    integer :: i, j, k, l, m, mm, sgn_dw
    real(rkind) :: df, uu, vv, ww, h, qq, cc, c, ci, b2, b1
    real(rkind), dimension(indx_cp_l:indx_cp_r+1) :: cp_coeff_ompc
    real(rkind), dimension(:,:,:,:) :: fl_ompc
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,:) :: w_aux_ompc, w_ompc
    real(rkind), dimension(:) :: dcsidx_ompc
    real(rkind), dimension(:) :: winf_ompc
    real(rkind) :: rho,tt,gamloc,p_rho,p_e,etot,b3

    !$omp parallel do default(firstprivate) shared(cp_coeff_ompc, &
    !$omp& fl_ompc, &
    !$omp& w_aux_ompc, &
    !$omp& w_ompc, &
    !$omp& dcsidx_ompc, &
    !$omp& winf_ompc)
    do k = 1,nz
      do j = 1,ny
        c_one = [-1.5_rkind, 2._rkind, -0.5_rkind]
        if(start_or_end == 1) then
          i = 1
          sgn_dw = 1
        elseif(start_or_end == 2) then
          i = nx
          sgn_dw = -1
        endif

        do m=1,5
          dw_dn(m) = 0._rkind
          do l=1,3
            dw_dn(m) = dw_dn(m) + sgn_dw * c_one(l)*w_ompc(i+sgn_dw*(l-1),j,k,m)
          enddo
          w_target = 0._rkind
          if (nr_type == 2) then
            w_target = w_ompc(i-sgn_dw,j,k,m)
          elseif(nr_type == 3) then
            w_target = winf_ompc(m)
          endif
          dw_dn_outer(m) = sgn_dw * (w_ompc(i,j,k,m)-w_target)
        enddo

        rho = w_aux_ompc(i,j,k,1)
        uu = w_aux_ompc(i,j,k,2)
        vv = w_aux_ompc(i,j,k,3)
        ww = w_aux_ompc(i,j,k,4)
        h = w_aux_ompc(i,j,k,5)
        tt = w_aux_ompc(i,j,k,6)
        qq = 0.5_rkind * (uu*uu +vv*vv + ww*ww)
        gamloc = get_gamloc_dev(indx_cp_l,indx_cp_r,tt,t0,cp_coeff_ompc,calorically_perfect,rgas0)
        cc = gamloc * tt * rgas0
        c = sqrt(cc)
        ci = 1._rkind/c
        p_rho = tt*rgas0
        p_e = rho*(gamloc-1._rkind)
        etot = h - tt*rgas0
        b3 = etot - rho * p_rho/p_e
        b2 = p_e/(rho*cc)
        b1 = p_rho/cc - b2*(etot - 2._rkind*qq)
        call eigenvectors_x(b1, b2, b3, uu, vv, ww, c, ci, h, el, er)

        do m=1,5
          dwc_dn(m) = 0._rkind
          dwc_dn_outer(m) = 0._rkind
          do mm=1,5
            dwc_dn(m) = dwc_dn(m) + el(mm,m) * dw_dn(mm)
            dwc_dn_outer(m) = dwc_dn_outer(m) + el(mm,m) * dw_dn_outer(mm)
          enddo
        enddo

        ev(1) = uu-c
        ev(2) = uu
        ev(3) = ev(2)
        ev(4) = ev(2)
        ev(5) = uu+c

        if (nr_type == 1) then
          do l=1,5
            ev(l) = sgn_dw*min(sgn_dw*ev(l) ,0._rkind)
          enddo
        endif

        if (nr_type == 2 .or. nr_type == 3) then
          do m=1,5
            if(sgn_dw*ev(m) > 0._rkind) then
              dwc_dn(m) = dwc_dn_outer(m)
            endif
          enddo
        endif

        do m=1,5
          dwc_dn(m) = ev(m) * dwc_dn(m)
        enddo

        if (nr_type == 6) then
          dwc_dn(2) = 0._rkind
          dwc_dn(3) = 0._rkind
          dwc_dn(4) = 0._rkind
          if(start_or_end == 1) then
            dwc_dn(5) = dwc_dn(1)
          elseif(start_or_end == 2) then
            dwc_dn(1) = dwc_dn(5)
          endif
        endif

        do m=1,5
          df = 0._rkind
          do mm=1,5
            df = df + er(mm,m) * dwc_dn(mm)
          enddo
          fl_ompc(i,j,k,m) = fl_ompc(i,j,k,m) + df * dcsidx_ompc(i)
        enddo
      enddo
    enddo


  endsubroutine bc_nr_lat_x_kernel


  subroutine bc_nr_lat_x_c2_kernel(start_or_end,nr_type,nx,ny,nz,ng,nv,w_aux_ompc,w_ompc,fl_ompc,&
  &dcsidx_ompc,indx_cp_l,indx_cp_r,cp_coeff_ompc,winf_ompc,jac_ompc,mcsi_ompc,dcsidxnc2_ompc,&
  &dcsidync2_ompc,dcsidxc2_ompc,dcsidyc2_ompc,calorically_perfect,rgas0,t0)
    integer, intent(in), value :: start_or_end, nr_type, nx, ny, nz, ng, nv, indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind), intent(in), value :: rgas0,t0
    real(rkind), dimension(3) :: c_one
    real(rkind), dimension(5) :: dw_dn, dwc_dn, ev, dw_dn_outer, dwc_dn_outer
    real(rkind), dimension(5,5) :: el, er
    real(rkind) :: w_target
    integer :: i, j, k, l, m, mm, sgn_dw
    real(rkind) :: df, uu, vv, ww, h, qq, cc, c, ci, b2, b1
    real(rkind), dimension(indx_cp_l:indx_cp_r+1) :: cp_coeff_ompc
    real(rkind), dimension(:,:,:,:) :: fl_ompc
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,:) :: w_aux_ompc, w_ompc
    real(rkind), dimension(:) :: dcsidx_ompc
    real(rkind), dimension(:) :: winf_ompc
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng), intent(in) :: jac_ompc, mcsi_ompc, dcsidxnc2_ompc, dcsidync2_ompc
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng), intent(in) :: dcsidxc2_ompc, dcsidyc2_ompc
    real(rkind) :: dcsixjdcsi, dcsiyjdcsi, pp, ut
    real(rkind) :: rho,tt,gamloc,p_rho,p_e,etot,b3

    !$omp parallel do default(firstprivate) shared(cp_coeff_ompc, &
    !$omp& fl_ompc, &
    !$omp& w_aux_ompc, &
    !$omp& w_ompc, &
    !$omp& dcsidx_ompc, &
    !$omp& winf_ompc, &
    !$omp& jac_ompc, &
    !$omp& mcsi_ompc, &
    !$omp& dcsidxnc2_ompc, &
    !$omp& dcsidync2_ompc, &
    !$omp& dcsidxc2_ompc, &
    !$omp& dcsidyc2_ompc)
    do k = 1,nz
      do j = 1,ny
        c_one = [-1.5_rkind, 2._rkind, -0.5_rkind]
        if(start_or_end == 1) then
          i = 1
          sgn_dw = 1
        elseif(start_or_end == 2) then
          i = nx
          sgn_dw = -1
        endif

        do m=1,nv
          dw_dn(m) = 0._rkind
          do l=1,3
            dw_dn(m) = dw_dn(m) + sgn_dw * c_one(l)*w_ompc(i+sgn_dw*(l-1),j,k,m)
          enddo
          w_target = 0._rkind
          if (nr_type == 2) then
            w_target = w_ompc(i-sgn_dw,j,k,m)
          elseif(nr_type == 3) then
            w_target = winf_ompc(m)
          endif
          dw_dn_outer(m) = sgn_dw * (w_ompc(i,j,k,m)-w_target)
        enddo

        rho = w_aux_ompc(i,j,k,1)
        uu = w_aux_ompc(i,j,k,2)
        vv = w_aux_ompc(i,j,k,3)
        ww = w_aux_ompc(i,j,k,4)
        h = w_aux_ompc(i,j,k,5)
        tt = w_aux_ompc(i,j,k,6)
        qq = 0.5_rkind * (uu*uu +vv*vv + ww*ww)
        gamloc = get_gamloc_dev(indx_cp_l,indx_cp_r,tt,t0,cp_coeff_ompc,calorically_perfect,rgas0)
        cc = gamloc * tt * rgas0
        c = sqrt(cc)
        ci = 1._rkind/c
        p_rho = tt*rgas0
        p_e = rho*(gamloc-1._rkind)
        etot = h - tt*rgas0
        b3 = etot - rho * p_rho/p_e
        b2 = p_e/(rho*cc)
        b1 = p_rho/cc - b2*(etot - 2._rkind*qq)
        ut = dcsidxnc2_ompc(i,j) * uu + dcsidync2_ompc(i,j) * vv
        call eigenvectors_x_c2(b1, b2, b3, uu, vv, ww, c, ci, h, el, er, ut, dcsidxnc2_ompc(i,j), dcsidync2_ompc(i,j))

        do m=1,5
          dwc_dn(m) = 0._rkind
          dwc_dn_outer(m) = 0._rkind
          do mm=1,5
            dwc_dn(m) = dwc_dn(m) + el(mm,m) * dw_dn(mm)
            dwc_dn_outer(m) = dwc_dn_outer(m) + el(mm,m) * dw_dn_outer(mm)
          enddo
        enddo

        ev(1) = ut-c
        ev(2) = ut
        ev(3) = ev(2)
        ev(4) = ev(2)
        ev(5) = ut+c

        if (nr_type == 1) then
          do l=1,5
            ev(l) = sgn_dw*min(sgn_dw*ev(l) ,0._rkind)
          enddo
        endif

        if (nr_type == 2 .or. nr_type == 3) then
          do m=1,5
            if(sgn_dw*ev(m) > 0._rkind) then
              dwc_dn(m) = dwc_dn_outer(m)
            endif
          enddo
        endif

        do m=1,5
          dwc_dn(m) = ev(m) * dwc_dn(m)
        enddo

        if (nr_type == 6) then
          dwc_dn(2) = 0._rkind
          dwc_dn(3) = 0._rkind
          dwc_dn(4) = 0._rkind
          if (start_or_end == 1) then
            dwc_dn(5) = dwc_dn(1)
          elseif(start_or_end == 2) then
            dwc_dn(1) = dwc_dn(5)
          endif
        endif

        do m=1,5
          df = 0._rkind
          do mm=1,5
            df = df + er(mm,m) * dwc_dn(mm)
          enddo
          fl_ompc(i,j,k,m) = fl_ompc(i,j,k,m) + df * mcsi_ompc(i,j)
        enddo


        dcsixjdcsi = sgn_dw * (c_one(1)*dcsidxc2_ompc(i,j)/jac_ompc(i,j)+c_one(2)*dcsidxc2_ompc(i+1*&
        &sgn_dw,j)/jac_ompc(i+1*sgn_dw,j)+ c_one(3)*dcsidxc2_ompc(i+2*sgn_dw,j)/jac_ompc(i+2*sgn_dw,j))
        dcsiyjdcsi = sgn_dw * (c_one(1)*dcsidyc2_ompc(i,j)/jac_ompc(i,j)+c_one(2)*dcsidyc2_ompc(i+1*&
        &sgn_dw,j)/jac_ompc(i+1*sgn_dw,j)+ c_one(3)*dcsidyc2_ompc(i+2*sgn_dw,j)/jac_ompc(i+2*sgn_dw,j))
        pp = p_rho * rho
        fl_ompc(i,j,k,1) = fl_ompc(i,j,k,1)+( w_ompc(i,j,k,1)*uu *dcsixjdcsi+ w_ompc(i,j,k,1)*vv *dcsiyjdcsi)*jac_ompc(i,j)
        fl_ompc(i,j,k,2) = fl_ompc(i,j,k,2)+((w_ompc(i,j,k,2)*uu+pp)*dcsixjdcsi+ w_ompc(i,j,k,2)*vv *dcsiyjdcsi)*jac_ompc(i,j)
        fl_ompc(i,j,k,3) = fl_ompc(i,j,k,3)+( w_ompc(i,j,k,3)*uu *dcsixjdcsi+(w_ompc(i,j,k,3)*vv+pp)*dcsiyjdcsi)*jac_ompc(i,j)
        fl_ompc(i,j,k,4) = fl_ompc(i,j,k,4)+( w_ompc(i,j,k,4)*uu *dcsixjdcsi+ w_ompc(i,j,k,4)*vv *dcsiyjdcsi)*jac_ompc(i,j)
        fl_ompc(i,j,k,5) = fl_ompc(i,j,k,5)+((w_ompc(i,j,k,5)+pp)*uu*dcsixjdcsi+(w_ompc(i,j,k,5)+pp)*vv*dcsiyjdcsi)*jac_ompc(i,j)
      enddo
    enddo


  endsubroutine bc_nr_lat_x_c2_kernel


  subroutine bc_nr_lat_y_kernel(start_or_end,nr_type,nx,ny,nz,ng,nv,w_aux_ompc,w_ompc,fl_ompc,&
  &detady_ompc,indx_cp_l,indx_cp_r,cp_coeff_ompc,winf_ompc,calorically_perfect,rgas0,t0)
    integer, intent(in), value :: start_or_end, nr_type, nx, ny, nz, ng, nv, indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind), intent(in), value :: rgas0,t0
    real(rkind), dimension(3) :: c_one
    real(rkind), dimension(5) :: dw_dn, dwc_dn, ev, dw_dn_outer, dwc_dn_outer
    real(rkind), dimension(5,5) :: el, er
    real(rkind) :: w_target
    integer :: i, j, k, l, m, mm, sgn_dw
    real(rkind) :: df, uu, vv, ww, h, qq, cc, c, ci, b2, b1
    real(rkind), dimension(indx_cp_l:indx_cp_r+1) :: cp_coeff_ompc
    real(rkind), dimension(:,:,:,:) :: fl_ompc
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: w_aux_ompc, w_ompc
    real(rkind), dimension(:) :: detady_ompc
    real(rkind), dimension(:) :: winf_ompc
    real(rkind) :: rho,tt,gamloc,p_rho,p_e,etot,b3

    !$omp parallel do default(firstprivate) shared(cp_coeff_ompc, &
    !$omp& fl_ompc, &
    !$omp& w_aux_ompc, &
    !$omp& w_ompc, &
    !$omp& detady_ompc, &
    !$omp& winf_ompc)
    do k = 1,nz
      do i = 1,nx
        c_one = [-1.5_rkind, 2._rkind, -0.5_rkind]
        if(start_or_end == 1) then
          j = 1
          sgn_dw = 1
        elseif(start_or_end == 2) then
          j = ny
          sgn_dw = -1
        endif

        do m=1,5
          dw_dn(m) = 0._rkind
          do l=1,3
            dw_dn(m) = dw_dn(m) + sgn_dw * c_one(l)*w_ompc(i,j+sgn_dw*(l-1),k,m)
          enddo
          w_target = 0._rkind
          if (nr_type == 2) then
            w_target = w_ompc(i,j-sgn_dw,k,m)
          elseif(nr_type == 3) then
            w_target = winf_ompc(m)
          endif
          dw_dn_outer(m) = sgn_dw * (w_ompc(i,j,k,m)-w_target)
        enddo

        rho = w_aux_ompc(i,j,k,1)
        uu = w_aux_ompc(i,j,k,2)
        vv = w_aux_ompc(i,j,k,3)
        ww = w_aux_ompc(i,j,k,4)
        h = w_aux_ompc(i,j,k,5)
        tt = w_aux_ompc(i,j,k,6)
        qq = 0.5_rkind * (uu*uu +vv*vv + ww*ww)
        gamloc = get_gamloc_dev(indx_cp_l,indx_cp_r,tt,t0,cp_coeff_ompc,calorically_perfect,rgas0)
        cc = gamloc * tt * rgas0
        c = sqrt(cc)
        ci = 1._rkind/c
        p_rho = tt*rgas0
        p_e = rho*(gamloc-1._rkind)
        etot = h - tt*rgas0
        b3 = etot - rho * p_rho/p_e
        b2 = p_e/(rho*cc)
        b1 = p_rho/cc - b2*(etot - 2._rkind*qq)
        call eigenvectors_y(b1, b2, b3, uu, vv, ww, c, ci, h, el, er)

        do m=1,5
          dwc_dn(m) = 0._rkind
          dwc_dn_outer(m) = 0._rkind
          do mm=1,5
            dwc_dn(m) = dwc_dn(m) + el(mm,m) * dw_dn(mm)
            dwc_dn_outer(m) = dwc_dn_outer(m) + el(mm,m) * dw_dn_outer(mm)
          enddo
        enddo

        ev(1) = vv-c
        ev(2) = vv
        ev(3) = ev(2)
        ev(4) = ev(2)
        ev(5) = vv+c

        if (nr_type == 1) then
          if(start_or_end == 1) then
            do l=1,5
              ev(l) = min(ev(l) ,0._rkind)
            enddo
          elseif(start_or_end == 2) then
            do l=1,5
              ev(l) = max(ev(l) ,0._rkind)
            enddo
          endif
        endif

        if (nr_type == 2 .or. nr_type == 3) then
          do m=1,5
            if(sgn_dw*ev(m) > 0._rkind) then
              dwc_dn(m) = dwc_dn_outer(m)
            endif
          enddo
        endif

        do m=1,5
          dwc_dn(m) = ev(m) * dwc_dn(m)
        enddo

        if (nr_type == 6) then
          dwc_dn(2) = 0._rkind
          dwc_dn(3) = 0._rkind
          dwc_dn(4) = 0._rkind
          if(start_or_end == 1) then
            dwc_dn(5) = dwc_dn(1)
          elseif(start_or_end == 2) then
            dwc_dn(1) = dwc_dn(5)
          endif
        endif

        do m=1,5
          df = 0._rkind
          do mm=1,5
            df = df + er(mm,m) * dwc_dn(mm)
          enddo
          fl_ompc(i,j,k,m) = fl_ompc(i,j,k,m) + df * detady_ompc(j)
        enddo
      enddo
    enddo


  endsubroutine bc_nr_lat_y_kernel


  subroutine bc_nr_lat_y_c2_kernel(start_or_end,nr_type,nx,ny,nz,ng,nv,w_aux_ompc,w_ompc,fl_ompc,&
  &detady_ompc,indx_cp_l,indx_cp_r,cp_coeff_ompc,wfar_ompc,jac_ompc,meta_ompc,detadxnc2_ompc,&
  &detadync2_ompc,detadxc2_ompc,detadyc2_ompc,wall_tag_ompc,calorically_perfect,rgas0,t0)
    integer, intent(in), value :: start_or_end, nr_type, nx, ny, nz, ng, nv, indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind), intent(in), value :: rgas0,t0
    real(rkind), dimension(3) :: c_one
    real(rkind), dimension(5) :: dw_dn, dwc_dn, ev, dw_dn_outer, dwc_dn_outer
    real(rkind), dimension(5,5) :: el, er
    real(rkind) :: w_target
    integer :: i, j, k, l, m, mm, sgn_dw
    real(rkind) :: df, uu, vv, ww, h, qq, cc, c, ci, b2, b1
    real(rkind), dimension(indx_cp_l:indx_cp_r+1) :: cp_coeff_ompc
    real(rkind), dimension(:,:,:,:) :: fl_ompc
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: w_aux_ompc, w_ompc
    real(rkind), dimension(:) :: detady_ompc
    real(rkind), dimension(nx,nv) :: wfar_ompc
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng), intent(in) :: jac_ompc, meta_ompc, detadxnc2_ompc, detadync2_ompc
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng), intent(in) :: detadxc2_ompc, detadyc2_ompc
    integer, dimension(1-ng:nx+ng), intent(in) :: wall_tag_ompc
    real(rkind) :: rho,tt,gamloc,p_rho,p_e,etot,b3
    real(rkind) :: detaxjdeta, detayjdeta, pp, ut

    !$omp parallel do default(firstprivate) shared(cp_coeff_ompc, &
    !$omp& fl_ompc, &
    !$omp& w_aux_ompc, &
    !$omp& w_ompc, &
    !$omp& detady_ompc, &
    !$omp& wfar_ompc, &
    !$omp& jac_ompc, &
    !$omp& meta_ompc, &
    !$omp& detadxnc2_ompc, &
    !$omp& detadync2_ompc, &
    !$omp& detadxc2_ompc, &
    !$omp& detadyc2_ompc, &
    !$omp& wall_tag_ompc)
    do k = 1,nz
      do i = 1,nx
        c_one = [-1.5_rkind, 2._rkind, -0.5_rkind]
        if(start_or_end == 1) then
          j = 1
          sgn_dw = 1
        elseif(start_or_end == 2) then
          j = ny
          sgn_dw = -1
        endif

        if(((start_or_end == 1) .and. (wall_tag_ompc(i) > 0)) ) then
          sgn_dw = 1
        else
          do m=1,nv
            dw_dn(m) = 0._rkind
            do l=1,3
              dw_dn(m) = dw_dn(m) + sgn_dw * c_one(l)*w_ompc(i,j+sgn_dw*(l-1),k,m)
            enddo
            w_target = 0._rkind
            if (nr_type == 2) then
              w_target = w_ompc(i,j-sgn_dw,k,m)
            elseif(nr_type == 3) then
              w_target = wfar_ompc(i,m)
            endif
            dw_dn_outer(m) = sgn_dw * (w_ompc(i,j,k,m)-w_target)
          enddo

          rho = w_aux_ompc(i,j,k,1)
          uu = w_aux_ompc(i,j,k,2)
          vv = w_aux_ompc(i,j,k,3)
          ww = w_aux_ompc(i,j,k,4)
          h = w_aux_ompc(i,j,k,5)
          tt = w_aux_ompc(i,j,k,6)
          qq = 0.5_rkind * (uu*uu +vv*vv + ww*ww)
          gamloc = get_gamloc_dev(indx_cp_l,indx_cp_r,tt,t0,cp_coeff_ompc,calorically_perfect,rgas0)
          cc = gamloc * tt * rgas0
          c = sqrt(cc)
          ci = 1._rkind/c
          p_rho = tt*rgas0
          p_e = rho*(gamloc-1._rkind)
          etot = h - tt*rgas0
          b3 = etot - rho * p_rho/p_e
          b2 = p_e/(rho*cc)
          b1 = p_rho/cc - b2*(etot - 2._rkind*qq)
          ut = detadxnc2_ompc(i,j) * uu + detadync2_ompc(i,j) * vv
          call eigenvectors_y_c2(b1, b2, b3, uu, vv, ww, c, ci, h, el, er, ut, detadxnc2_ompc(i,j), detadync2_ompc(i,j))

          do m=1,5
            dwc_dn(m) = 0._rkind
            dwc_dn_outer(m) = 0._rkind
            do mm=1,5
              dwc_dn(m) = dwc_dn(m) + el(mm,m) * dw_dn(mm)
              dwc_dn_outer(m) = dwc_dn_outer(m) + el(mm,m) * dw_dn_outer(mm)
            enddo
          enddo

          ev(1) = ut-c
          ev(2) = ut
          ev(3) = ev(2)
          ev(4) = ev(2)
          ev(5) = ut+c

          if(nr_type == 1) then
            if(start_or_end == 1) then
              do l=1,5
                ev(l) = min(ev(l) ,0._rkind)
              enddo
            elseif(start_or_end == 2) then
              do l=1,5
                ev(l) = max(ev(l) ,0._rkind)
              enddo
            endif
          endif

          if (nr_type == 2 .or. nr_type == 3) then
            do m=1,5
              if(sgn_dw*ev(m) > 0._rkind) then
                dwc_dn(m) = dwc_dn_outer(m)
              endif
            enddo
          endif

          do m=1,5
            dwc_dn(m) = ev(m) * dwc_dn(m)
          enddo

          if (nr_type == 6) then
            dwc_dn(2) = 0._rkind
            dwc_dn(3) = 0._rkind
            dwc_dn(4) = 0._rkind
            if(start_or_end == 1) then
              dwc_dn(5) = dwc_dn(1)
            elseif(start_or_end == 2) then
              dwc_dn(1) = dwc_dn(5)
            endif
          endif

          do m=1,5
            df = 0._rkind
            do mm=1,5
              df = df + er(mm,m) * dwc_dn(mm)
            enddo
            fl_ompc(i,j,k,m) = fl_ompc(i,j,k,m) + df * meta_ompc(i,j)
          enddo

          detaxjdeta=sgn_dw*(c_one(1)*detadxc2_ompc(i,j)/jac_ompc(i,j)+c_one(2)*detadxc2_ompc(i,j+1*&
          &sgn_dw)/jac_ompc(i,j+1*sgn_dw)+ c_one(3)*detadxc2_ompc(i,j+2*sgn_dw)/jac_ompc(i,j+2*sgn_dw))
          detayjdeta=sgn_dw*(c_one(1)*detadyc2_ompc(i,j)/jac_ompc(i,j)+c_one(2)*detadyc2_ompc(i,j+1*&
          &sgn_dw)/jac_ompc(i,j+1*sgn_dw)+ c_one(3)*detadyc2_ompc(i,j+2*sgn_dw)/jac_ompc(i,j+2*sgn_dw))

          pp = p_rho * rho
          fl_ompc(i,j,k,1) = fl_ompc(i,j,k,1) +( w_ompc(i,j,k,1)*uu *detaxjdeta + w_ompc(i,j,k,1)*vv *detayjdeta)*jac_ompc(i,j)
          fl_ompc(i,j,k,2) = fl_ompc(i,j,k,2) +((w_ompc(i,j,k,2)*uu+pp)*detaxjdeta + w_ompc(i,j,k,2)*vv *detayjdeta)*jac_ompc(i,j)
          fl_ompc(i,j,k,3) = fl_ompc(i,j,k,3) +( w_ompc(i,j,k,3)*uu *detaxjdeta + (w_ompc(i,j,k,3)*vv+pp)*detayjdeta)*jac_ompc(i,j)
          fl_ompc(i,j,k,4) = fl_ompc(i,j,k,4) +( w_ompc(i,j,k,4)*uu *detaxjdeta + w_ompc(i,j,k,4)*vv *detayjdeta)*jac_ompc(i,j)
          fl_ompc(i,j,k,5) = fl_ompc(i,j,k,5) +((w_ompc(i,j,k,5)+pp)*uu*detaxjdeta + (w_ompc(i,j,k,&
          &5)+pp)*vv*detayjdeta)*jac_ompc(i,j)
        endif
      enddo
    enddo


  endsubroutine bc_nr_lat_y_c2_kernel


  subroutine bc_nr_lat_z_kernel(start_or_end,nr_type,nx,ny,nz,ng,nv,w_aux_ompc,w_ompc,fl_ompc,&
  &dzitdz_ompc,indx_cp_l,indx_cp_r,cp_coeff_ompc,winf_ompc,calorically_perfect,rgas0,t0)
    integer, intent(in), value :: start_or_end, nr_type, nx, ny, nz, ng, nv, indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind), intent(in), value :: rgas0,t0
    real(rkind), dimension(3) :: c_one
    real(rkind), dimension(5) :: dw_dn, dwc_dn, ev, dw_dn_outer, dwc_dn_outer
    real(rkind), dimension(5,5) :: el, er
    real(rkind) :: w_target
    integer :: i, j, k, l, m, mm, sgn_dw
    real(rkind) :: df, uu, vv, ww, h, qq, cc, c, ci, b2, b1
    real(rkind), dimension(indx_cp_l:indx_cp_r+1) :: cp_coeff_ompc
    real(rkind), dimension(:,:,:,:) :: fl_ompc
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: w_aux_ompc, w_ompc
    real(rkind), dimension(:) :: dzitdz_ompc
    real(rkind), dimension(:) :: winf_ompc
    real(rkind) :: rho,tt,gamloc,p_rho,p_e,etot,b3

    !$omp parallel do default(firstprivate) shared(cp_coeff_ompc, &
    !$omp& fl_ompc, &
    !$omp& w_aux_ompc, &
    !$omp& w_ompc, &
    !$omp& dzitdz_ompc, &
    !$omp& winf_ompc)
    do j = 1,ny
      do i = 1,nx
        c_one = [-1.5_rkind, 2._rkind, -0.5_rkind]
        if(start_or_end == 1) then
          k = 1
          sgn_dw = 1
        elseif(start_or_end == 2) then
          k = nz
          sgn_dw = -1
        endif

        do m=1,5
          dw_dn(m) = 0._rkind
          do l=1,3
            dw_dn(m) = dw_dn(m) + sgn_dw * c_one(l)*w_ompc(i,j,k+sgn_dw*(l-1),m)
          enddo
          if (nr_type == 2) then
            w_target = w_ompc(i,j,k-sgn_dw,m)
          elseif(nr_type == 3) then
            w_target = winf_ompc(m)
          endif
          dw_dn_outer(m) = sgn_dw * (w_ompc(i,j,k,m)-w_target)
        enddo

        rho = w_aux_ompc(i,j,k,1)
        uu = w_aux_ompc(i,j,k,2)
        vv = w_aux_ompc(i,j,k,3)
        ww = w_aux_ompc(i,j,k,4)
        h = w_aux_ompc(i,j,k,5)
        tt = w_aux_ompc(i,j,k,6)
        qq = 0.5_rkind * (uu*uu +vv*vv + ww*ww)
        gamloc = get_gamloc_dev(indx_cp_l,indx_cp_r,tt,t0,cp_coeff_ompc,calorically_perfect,rgas0)
        cc = gamloc * tt * rgas0
        c = sqrt(cc)
        ci = 1._rkind/c
        p_rho = tt*rgas0
        p_e = rho*(gamloc-1._rkind)
        etot = h - tt*rgas0
        b3 = etot - rho * p_rho/p_e
        b2 = p_e/(rho*cc)
        b1 = p_rho/cc - b2*(etot - 2._rkind*qq)
        call eigenvectors_z(b1, b2, b3, uu, vv, ww, c, ci, h, el, er)

        do m=1,5
          dwc_dn(m) = 0._rkind
          dwc_dn_outer(m) = 0._rkind
          do mm=1,5
            dwc_dn(m) = dwc_dn(m) + el(mm,m) * dw_dn(mm)
            dwc_dn_outer(m) = dwc_dn_outer(m) + el(mm,m) * dw_dn_outer(mm)
          enddo
        enddo

        ev(1) = ww-c
        ev(2) = ww
        ev(3) = ev(2)
        ev(4) = ev(2)
        ev(5) = ww+c

        if (nr_type == 1) then
          if(start_or_end == 1) then
            do l=1,5
              ev(l) = min(ev(l) ,0._rkind)
            enddo
          elseif(start_or_end == 2) then
            do l=1,5
              ev(l) = max(ev(l) ,0._rkind)
            enddo
          endif
        endif

        if (nr_type == 2 .or. nr_type == 3) then
          do m=1,5
            if(sgn_dw*ev(m) > 0._rkind) then
              dwc_dn(m) = dwc_dn_outer(m)
            endif
          enddo
        endif

        do m=1,5
          dwc_dn(m) = ev(m) * dwc_dn(m)
        enddo

        if (nr_type == 6) then
          dwc_dn(2) = 0._rkind
          dwc_dn(3) = 0._rkind
          dwc_dn(4) = 0._rkind
          if(start_or_end == 1) then
            dwc_dn(5) = dwc_dn(1)
          elseif(start_or_end == 2) then
            dwc_dn(1) = dwc_dn(5)
          endif
        endif

        do m=1,5
          df = 0._rkind
          do mm=1,5
            df = df + er(mm,m) * dwc_dn(mm)
          enddo
          fl_ompc(i,j,k,m) = fl_ompc(i,j,k,m) + df * dzitdz_ompc(k)
        enddo
      enddo
    enddo


  endsubroutine bc_nr_lat_z_kernel


  subroutine bcrecyc_subroutine_1(nx,ny,nz,ng,nv,wrecycav_ompc,wrecyc_ompc)
    integer, intent(in) :: nx, ny, nz, ng, nv
    real(rkind), dimension(:,:,:), intent(inout) :: wrecycav_ompc
    real(rkind), dimension(:,:,:,:), intent(in) :: wrecyc_ompc
    integer :: i,j,k,m,iercuda

    !$omp parallel do default(firstprivate) shared(wrecycav_ompc, &
    !$omp& wrecyc_ompc)
    do j = 1,ny
      do i = 1,ng
        do m=1,nv
          wrecycav_ompc(i,j,m) = 0._rkind
          do k=1,nz
            wrecycav_ompc(i,j,m) = wrecycav_ompc(i,j,m)+wrecyc_ompc(i,j,k,m)
          enddo
        enddo
      enddo
    enddo


  endsubroutine bcrecyc_subroutine_1



  subroutine bcrecyc_subroutine_2(nx,ny,nz,nzmax,ng,wrecycav_ompc,wrecyc_ompc)
    integer, intent(in) :: nx, ny, nz, nzmax, ng
    real(rkind), dimension(:,:,:), intent(in) :: wrecycav_ompc
    real(rkind), dimension(:,:,:,:), intent(inout) :: wrecyc_ompc
    real(rkind) :: ufav, vfav, wfav, rhom
    integer :: i,j,k,iercuda

    !$omp parallel do default(firstprivate) shared(wrecycav_ompc, &
    !$omp& wrecyc_ompc)
    do j = 1,ny
      do i = 1,ng
        ufav = wrecycav_ompc(i,j,2)/wrecycav_ompc(i,j,1)
        vfav = wrecycav_ompc(i,j,3)/wrecycav_ompc(i,j,1)
        wfav = wrecycav_ompc(i,j,4)/wrecycav_ompc(i,j,1)
        rhom = wrecycav_ompc(i,j,1)/nzmax
        do k=1,nz
          wrecyc_ompc(i,j,k,2) = wrecyc_ompc(i,j,k,2)/wrecyc_ompc(i,j,k,1)-ufav
          wrecyc_ompc(i,j,k,3) = wrecyc_ompc(i,j,k,3)/wrecyc_ompc(i,j,k,1)-vfav
          wrecyc_ompc(i,j,k,4) = wrecyc_ompc(i,j,k,4)/wrecyc_ompc(i,j,k,1)-wfav
          wrecyc_ompc(i,j,k,1) = wrecyc_ompc(i,j,k,1)-rhom
        enddo
      enddo
    enddo


  endsubroutine bcrecyc_subroutine_2




  subroutine bcrecyc_subroutine_3(nx,ny,nz,ng,p0,u0,rgas0,w_ompc,wmean_ompc,wrecyc_ompc,&
  &weta_inflow_ompc,map_j_inn_ompc,map_j_out_ompc,map_j_out_blend_ompc,yplus_inflow_ompc,&
  &eta_inflow_ompc,yplus_recyc_ompc,eta_recyc_ompc,eta_recyc_blend_ompc,betarecyc,glund1,&
  &inflow_random_plane_ompc,indx_cp_l,indx_cp_r,cv_coeff_ompc,cp_coeff_ompc,t0,calorically_perfect,&
  &rand_type)
    integer, intent(in) :: nx, ny, nz, ng, indx_cp_l, indx_cp_r, calorically_perfect, rand_type
    real(rkind) :: p0, rgas0, betarecyc, glund1, t0, u0
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in) :: cv_coeff_ompc, cp_coeff_ompc
    real(rkind), dimension(1-ng:,:,:), intent(in) :: wmean_ompc
    real(rkind), dimension(:,:,:,:), intent(in) :: wrecyc_ompc
    real(rkind), dimension(:,:,:), intent(in) :: inflow_random_plane_ompc
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,:), intent(inout) :: w_ompc
    real(rkind), dimension(:), intent(in) :: weta_inflow_ompc
    real(rkind), dimension(1-ng:), intent(in) :: yplus_inflow_ompc, eta_inflow_ompc, yplus_recyc_ompc, eta_recyc_ompc
    real(rkind), dimension(1-ng:), intent(in) :: eta_recyc_blend_ompc
    integer, dimension(:), intent(in) :: map_j_inn_ompc, map_j_out_ompc, map_j_out_blend_ompc
    integer :: i,j,k,iercuda
    real(rkind) :: eta, weta, weta1, bdamp, disty_inn, disty_out, rhofluc, ufluc, vfluc, wfluc, rhof_inn, rhof_out
    real(rkind) :: uf_inn, uf_out, vf_inn, vf_out, wf_inn, wf_out, etamin
    real(rkind) :: rhomean, uumean , vvmean , wwmean , tmean , rho , uu , vv , ww , rhou , rhov , rhow, tt, ee
    real(rkind) :: u0_02, tfluc
    integer :: j_inn, j_out
    u0_02 = 0.02_rkind*u0
    if (rand_type==0) u0_02 = 0._rkind
    !$omp parallel do default(firstprivate) shared(cv_coeff_ompc, &
    !$omp& cp_coeff_ompc, &
    !$omp& wmean_ompc, &
    !$omp& wrecyc_ompc, &
    !$omp& inflow_random_plane_ompc, &
    !$omp& w_ompc, &
    !$omp& weta_inflow_ompc, &
    !$omp& yplus_inflow_ompc, &
    !$omp& eta_inflow_ompc, &
    !$omp& yplus_recyc_ompc, &
    !$omp& eta_recyc_ompc, &
    !$omp& eta_recyc_blend_ompc, &
    !$omp& map_j_inn_ompc, &
    !$omp& map_j_out_ompc, &
    !$omp& map_j_out_blend_ompc)
    do k = 1,nz
      do j = 1,ny
        eta = eta_inflow_ompc(j)
        etamin = min(eta,1._rkind)
        weta = weta_inflow_ompc(j)
        weta1 = 1._rkind-weta
        bdamp = 0.5_rkind*(1._rkind-tanh(4._rkind*(eta_inflow_ompc(j)-2._rkind)))
        j_inn = map_j_inn_ompc(j)
        j_out = map_j_out_ompc(j)
        disty_inn = (yplus_inflow_ompc(j)-yplus_recyc_ompc(j_inn))/(yplus_recyc_ompc(j_inn+1)-yplus_recyc_ompc(j_inn))
        disty_out = (eta_inflow_ompc(j)-eta_recyc_ompc(j_out))/(eta_recyc_ompc(j_out+1)-eta_recyc_ompc(j_out))
        do i=1,ng
          rhomean = wmean_ompc(1-i,j,1)
          uumean = wmean_ompc(1-i,j,2)/rhomean
          vvmean = wmean_ompc(1-i,j,3)/rhomean
          wwmean = wmean_ompc(1-i,j,4)/rhomean
          tmean = p0/rhomean/rgas0
          if (j==1.or.j_inn>=ny.or.j_out>=ny) then
            rhofluc = 0._rkind
            ufluc = 0._rkind
            vfluc = 0._rkind
            wfluc = 0._rkind
          else
            rhof_inn = wrecyc_ompc(i,j_inn,k,1)*(1._rkind-disty_inn)+wrecyc_ompc(i,j_inn+1,k,1)*disty_inn
            rhof_out = wrecyc_ompc(i,j_out,k,1)*(1._rkind-disty_out)+wrecyc_ompc(i,j_out+1,k,1)*disty_out
            uf_inn = wrecyc_ompc(i,j_inn,k,2)*(1._rkind-disty_inn)+wrecyc_ompc(i,j_inn+1,k,2)*disty_inn
            uf_out = wrecyc_ompc(i,j_out,k,2)*(1._rkind-disty_out)+wrecyc_ompc(i,j_out+1,k,2)*disty_out
            vf_inn = wrecyc_ompc(i,j_inn,k,3)*(1._rkind-disty_inn)+wrecyc_ompc(i,j_inn+1,k,3)*disty_inn
            vf_out = wrecyc_ompc(i,j_out,k,3)*(1._rkind-disty_out)+wrecyc_ompc(i,j_out+1,k,3)*disty_out
            wf_inn = wrecyc_ompc(i,j_inn,k,4)*(1._rkind-disty_inn)+wrecyc_ompc(i,j_inn+1,k,4)*disty_inn
            wf_out = wrecyc_ompc(i,j_out,k,4)*(1._rkind-disty_out)+wrecyc_ompc(i,j_out+1,k,4)*disty_out
            rhofluc = rhof_inn*weta1+rhof_out*weta
            ufluc = uf_inn*weta1+ uf_out*weta
            vfluc = vf_inn*weta1+ vf_out*weta
            wfluc = wf_inn*weta1+ wf_out*weta
            rhofluc = rhofluc*bdamp
            ufluc = ufluc *bdamp*betarecyc
            vfluc = vfluc *bdamp*betarecyc
            wfluc = wfluc *bdamp*betarecyc
            ufluc = ufluc+u0_02*(inflow_random_plane_ompc(j,k,1)-0.5_rkind)*etamin
            vfluc = vfluc+u0_02*(inflow_random_plane_ompc(j,k,2)-0.5_rkind)*etamin
            wfluc = wfluc+u0_02*(inflow_random_plane_ompc(j,k,3)-0.5_rkind)*etamin
          endif
          rhofluc = max(-0.75_rkind*rhomean,rhofluc)
          rhofluc = min( 4.0_rkind*rhomean,rhofluc)
          rho = rhomean + rhofluc
          tt = p0/rho/rgas0
          uu = uumean + ufluc
          vv = vvmean + vfluc
          ww = wwmean + wfluc
          rhou = rho*uu
          rhov = rho*vv
          rhow = rho*ww
          ee = get_e_from_temperature_dev(tt, t0, indx_cp_l, indx_cp_r, cv_coeff_ompc, calorically_perfect)
          w_ompc(1-i,j,k,1) = rho
          w_ompc(1-i,j,k,2) = rhou
          w_ompc(1-i,j,k,3) = rhov
          w_ompc(1-i,j,k,4) = rhow
          w_ompc(1-i,j,k,5) = rho*ee + 0.5_rkind*(rhou**2+rhov**2+rhow**2)/rho
        enddo
      enddo
    enddo


  endsubroutine bcrecyc_subroutine_3



  subroutine bclam_subroutine(ilat,nx,ny,nz,ng,nv,w_ompc,wmean_ompc,p0,rgas0,indx_cp_l,indx_cp_r,&
  &cv_coeff_ompc,t0,calorically_perfect)
    integer, intent(in) :: nx, ny, nz, ng, nv, indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind), intent(in) :: p0, rgas0, t0
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in) :: cv_coeff_ompc
    real(rkind), dimension(1-ng:, 1-ng:, 1-ng:, 1:), intent(inout) :: w_ompc
    real(rkind), dimension(1-ng:,:,:), intent(in) :: wmean_ompc
    integer :: ilat
    integer :: j,k,l,iercuda
    real(rkind) :: rho,rhou,rhov,rhow,tt,ee
    if (ilat==1) then
      !$omp parallel do default(firstprivate) shared(cv_coeff_ompc, &
      !$omp& w_ompc, &
      !$omp& wmean_ompc)
      do k = 1,nz
        do j = 1,ny
          do l = 1,ng
            rho = wmean_ompc(1-l,j,1)
            rhou = wmean_ompc(1-l,j,2)
            rhov = wmean_ompc(1-l,j,3)
            rhow = wmean_ompc(1-l,j,4)
            tt = p0/rho/rgas0
            w_ompc(1-l,j,k,1) = rho
            w_ompc(1-l,j,k,2) = rhou
            w_ompc(1-l,j,k,3) = rhov
            w_ompc(1-l,j,k,4) = rhow
            ee = get_e_from_temperature_dev(tt, t0, indx_cp_l, indx_cp_r, cv_coeff_ompc, calorically_perfect)
            w_ompc(1-l,j,k,5) = rho*ee + 0.5_rkind*(rhou**2+rhov**2+rhow**2)/rho
          enddo
        enddo
      enddo
    elseif (ilat==2) then
    elseif (ilat==3) then
    elseif (ilat==4) then
    elseif (ilat==5) then
    elseif (ilat==6) then
    endif

  endsubroutine bclam_subroutine



  subroutine bcfree_subroutine(ilat,nx,ny,nz,ng,nv,winf_ompc,w_ompc)
    integer, intent(in) :: ilat,nx,ny,nz,ng,nv
    real(rkind), dimension(1:), intent(in) :: winf_ompc
    real(rkind), dimension(1-ng:, 1-ng:, 1-ng:, 1:), intent(inout) :: w_ompc
    integer :: i,j,k,l,m,iercuda
    if (ilat==1) then
      !$omp parallel do default(firstprivate) shared(winf_ompc, &
      !$omp& w_ompc)
      do k = 1,nz
        do j = 1,ny
          do l = 1,ng
            do m=1,nv
              w_ompc(1-l,j,k,m) = winf_ompc(m)
            enddo
          enddo
        enddo
      enddo
    elseif (ilat==2) then
      !$omp parallel do default(firstprivate) shared(winf_ompc, &
      !$omp& w_ompc)
      do k = 1,nz
        do j = 1,ny
          do l=1,ng
            do m=1,nv
              w_ompc(nx+l,j,k,m) = winf_ompc(m)
            enddo
          enddo
        enddo
      enddo
    elseif (ilat==3) then
      !$omp parallel do default(firstprivate) shared(winf_ompc, &
      !$omp& w_ompc)
      do k = 1,nz
        do i = 1,nx
          do l=1,ng
            do m=1,nv
              w_ompc(i,1-l,k,m) = winf_ompc(m)
            enddo
          enddo
        enddo
      enddo
    elseif (ilat==4) then
      !$omp parallel do default(firstprivate) shared(winf_ompc, &
      !$omp& w_ompc)
      do k = 1,nz
        do i = 1,nx
          do l=1,ng
            do m=1,nv
              w_ompc(i,ny+l,k,m) = winf_ompc(m)
            enddo
          enddo
        enddo
      enddo
    elseif (ilat==5) then
      !$omp parallel do default(firstprivate) shared(winf_ompc, &
      !$omp& w_ompc)
      do j = 1,ny
        do i = 1,nx
          do l=1,ng
            do m=1,nv
              w_ompc(i,j,1-l,m) = winf_ompc(m)
            enddo
          enddo
        enddo
      enddo
    elseif (ilat==6) then
      !$omp parallel do default(firstprivate) shared(winf_ompc, &
      !$omp& w_ompc)
      do j = 1,ny
        do i = 1,nx
          do l=1,ng
            do m=1,nv
              w_ompc(i,j,nz+l,m) = winf_ompc(m)
            enddo
          enddo
        enddo
      enddo
    endif

  endsubroutine bcfree_subroutine



  subroutine bcfree_sub_subroutine(ilat,nx,ny,nz,ng,w_ompc,w_aux_ompc,aoa,t0,ptot0,ttot0,rgas0,&
  &indx_cp_l,indx_cp_r,cv_coeff_ompc,cp_coeff_ompc,calorically_perfect,tol_iter_nr)
    integer, intent(in) :: ilat,nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect
    real(rkind), intent(in) :: ptot0, ttot0, rgas0, aoa, t0, tol_iter_nr
    real(rkind), dimension(1-ng:, 1-ng:, 1-ng:, 1:), intent(inout) :: w_ompc
    real(rkind), dimension(1-ng:, 1-ng:, 1-ng:, 1:), intent(in) :: w_aux_ompc
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in) :: cv_coeff_ompc, cp_coeff_ompc
    real(rkind) :: rho,rhou,rhov,rhow,rhoe,ri,uu,vv,ww,qq,ee,tt,pp,del,gamloc,rmaf,rml,vel_mod
    real(rkind) :: rmfac,sinangle,cosangle
    integer :: i,j,k,l,m,iercuda
    cosangle = cos(aoa)
    sinangle = sin(aoa)

    if (ilat==1) then
      !$omp parallel do default(firstprivate) shared(w_ompc, &
      !$omp& w_aux_ompc, &
      !$omp& cv_coeff_ompc, &
      !$omp& cp_coeff_ompc)
      do k = 1,nz
        do j = 1,ny
          rho = w_ompc(1,j,k,1)
          rhou = w_ompc(1,j,k,2)
          rhov = w_ompc(1,j,k,3)
          rhow = w_ompc(1,j,k,4)
          rhoe = w_ompc(1,j,k,5)
          ri = 1._rkind/rho
          uu = rhou*ri
          vv = rhov*ri
          ww = rhow*ri
          qq = 0.5_rkind*(uu*uu+vv*vv+ww*ww)
          ee = rhoe/rho-qq
          tt = get_temperature_from_e_dev(ee,w_aux_ompc(1,j,k,6),t0,cv_coeff_ompc,indx_cp_l,&
          &indx_cp_r, calorically_perfect, tol_iter_nr)
          pp = rho*tt*rgas0
          gamloc = get_gamloc_dev(indx_cp_l,indx_cp_r,tt,t0,cp_coeff_ompc,calorically_perfect,rgas0)
          del = 0.5_rkind*(gamloc-1._rkind)
          rmfac = (pp/ptot0)**(-(gamloc-1._rkind)/gamloc)
          rml = sqrt((rmfac-1._rkind)/del)
          tt = ttot0/rmfac
          vel_mod = rml*sqrt(gamloc*rgas0*tt)
          uu = vel_mod*cosangle
          vv = vel_mod*sinangle
          ww = 0._rkind
          rho = pp/tt/rgas0
          rhou = rho*uu
          rhov = rho*vv
          rhow = rho*ww
          ee = get_e_from_temperature_dev(tt,t0,indx_cp_l,indx_cp_r,cv_coeff_ompc,calorically_perfect)
          do l=1,ng
            w_ompc(1-l,j,k,1) = rho
            w_ompc(1-l,j,k,2) = rhou
            w_ompc(1-l,j,k,3) = rhov
            w_ompc(1-l,j,k,4) = rhow
            w_ompc(1-l,j,k,5) = rho*ee+0.5_rkind*(rhou*rhou+rhov*rhov+rhow*rhow)/rho
          enddo
        enddo
      enddo
    elseif (ilat==4) then
      !$omp parallel do default(firstprivate) shared(w_ompc, &
      !$omp& w_aux_ompc, &
      !$omp& cv_coeff_ompc, &
      !$omp& cp_coeff_ompc)
      do k = 1,nz
        do i = 1,nx
          rho = w_ompc(i,ny,k,1)
          rhou = w_ompc(i,ny,k,2)
          rhov = w_ompc(i,ny,k,3)
          rhow = w_ompc(i,ny,k,4)
          rhoe = w_ompc(i,ny,k,5)
          ri = 1._rkind/rho
          uu = rhou*ri
          vv = rhov*ri
          ww = rhow*ri
          qq = 0.5_rkind*(uu*uu+vv*vv+ww*ww)
          ee = rhoe/rho-qq
          tt = get_temperature_from_e_dev(ee,w_aux_ompc(i,ny,k,6),t0,cv_coeff_ompc,indx_cp_l,&
          &indx_cp_r, calorically_perfect,tol_iter_nr)
          pp = rho*tt*rgas0
          gamloc = get_gamloc_dev(indx_cp_l,indx_cp_r,tt,t0,cp_coeff_ompc,calorically_perfect,rgas0)
          del = 0.5_rkind*(gamloc-1._rkind)
          rmfac = (pp/ptot0)**(-(gamloc-1._rkind)/gamloc)
          rml = sqrt((rmfac-1._rkind)/del)
          tt = ttot0/rmfac
          vel_mod = rml*sqrt(gamloc*rgas0*tt)
          uu = vel_mod*cosangle
          vv = vel_mod*sinangle
          ww = 0._rkind
          rho = pp/tt/rgas0
          rhou = rho*uu
          rhov = rho*vv
          rhow = rho*ww
          ee = get_e_from_temperature_dev(tt,t0,indx_cp_l,indx_cp_r,cv_coeff_ompc,calorically_perfect)
          do l=1,ng
            w_ompc(i,ny+l,k,1) = rho
            w_ompc(i,ny+l,k,2) = rhou
            w_ompc(i,ny+l,k,3) = rhov
            w_ompc(i,ny+l,k,4) = rhow
            w_ompc(i,ny+l,k,5) = rho*ee+0.5_rkind*(rhou*rhou+rhov*rhov+rhow*rhow)/rho
          enddo
        enddo
      enddo
    endif

  endsubroutine bcfree_sub_subroutine



  subroutine bcshock_subroutine(ilat,nx,ny,nz,ng,nv,w_ompc,winf_ompc,winf_past_shock_ompc,&
  &xshock_imp,shock_angle,x_ompc,y_ompc,tanhfacs)
    integer, intent(in) :: nx,ny,nz,ng,nv,ilat
    real(rkind), intent(in) :: xshock_imp, shock_angle, tanhfacs
    real(rkind), dimension(1:), intent(in) :: winf_ompc, winf_past_shock_ompc
    real(rkind), dimension(1-ng:, 1-ng:, 1-ng:, 1:), intent(inout) :: w_ompc
    real(rkind), dimension(1-ng:), intent(in) :: x_ompc, y_ompc
    integer :: i,k,l,m,iercuda
    real(rkind) :: xsh, xx, tanhlen, tanhf, dwinf
    tanhlen = 8._rkind*tanhfacs

    if (ilat==1) then
    elseif (ilat==2) then
    elseif (ilat==3) then
    elseif (ilat==4) then
      !$omp parallel do default(firstprivate) shared(winf_ompc, &
      !$omp& winf_past_shock_ompc, &
      !$omp& w_ompc, &
      !$omp& x_ompc, &
      !$omp& y_ompc)
      do k = 1,nz
        do i = 1,nx
          do l = 0,ng
            xsh = xshock_imp-y_ompc(ny+l)/tan(shock_angle)
            xx = x_ompc(i)-xsh
            if (abs(xx)>tanhlen) then
              do m=1,nv
                w_ompc(i,ny+l,k,m) = w_ompc(i,ny,k,m)
              enddo
            else
              tanhf = 0.5_rkind*(1._rkind+tanh(xx/tanhfacs))
              do m=1,nv
                dwinf = winf_past_shock_ompc(m)-winf_ompc(m)
                w_ompc(i,ny+l,k,m) = winf_ompc(m)+dwinf*tanhf
              enddo
            endif
          enddo
        enddo
      enddo
    elseif (ilat==5) then
    elseif (ilat==6) then
    endif

  endsubroutine bcshock_subroutine


  subroutine bcextr_var_subroutine(nx,ny,nz,ng,w_ompc)
    integer :: nx,ny,nz,ng
    integer :: ilat
    integer :: i,j,k,l,m,iercuda
    real(rkind), dimension(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng, 1:1), intent(inout) :: w_ompc

    !$omp parallel do default(firstprivate) shared(w_ompc)
    do k = 1,nz
      do j = 1,ny
        do l=1,ng
          w_ompc(1-l,j,k,1) = w_ompc(1,j,k,1)
        enddo
      enddo
    enddo

    !$omp parallel do default(firstprivate) shared(w_ompc)
    do k = 1,nz
      do j = 1,ny
        do l=1,ng
          w_ompc(nx+l,j,k,1) = w_ompc(nx,j,k,1)
        enddo
      enddo
    enddo

    !$omp parallel do default(firstprivate) shared(w_ompc)
    do k = 1,nz
      do i = 1,nx
        do l=1,ng
          w_ompc(i,1-l,k,1) = w_ompc(i,1,k,1)
        enddo
      enddo
    enddo

    !$omp parallel do default(firstprivate) shared(w_ompc)
    do k = 1,nz
      do i = 1,nx
        do l=1,ng
          w_ompc(i,ny+l,k,1) = w_ompc(i,ny,k,1)
        enddo
      enddo
    enddo

    !$omp parallel do default(firstprivate) shared(w_ompc)
    do j = 1,ny
      do i = 1,nx
        do l=1,ng
          w_ompc(i,j,1-l,1) = w_ompc(i,j,1,1)
        enddo
      enddo
    enddo

    !$omp parallel do default(firstprivate) shared(w_ompc)
    do j = 1,ny
      do i = 1,nx
        do l=1,ng
          w_ompc(i,j,nz+l,1) = w_ompc(i,j,nz,1)
        enddo
      enddo
    enddo


  endsubroutine bcextr_var_subroutine



  subroutine bcextr_airfoil_var_subroutine(nx,ny,nz,ng,ndim,wall_tag_ompc,w_ompc,ileftx,irightx,ileftz,irightz)
    integer :: nx,ny,nz,ng,ndim, ileftx, irightx, ileftz, irightz
    integer :: ilat
    integer :: i,j,k,l,m,iercuda
    real(rkind), dimension(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng, 1:1), intent(inout) :: w_ompc
    integer, dimension(1-ng:nx+ng), intent(in) :: wall_tag_ompc
    if (ileftx < 0) then
      !$omp parallel do default(firstprivate) shared(w_ompc, &
      !$omp& wall_tag_ompc)
      do k = 1,nz
        do j = 1,ny
          do l=1,ng
            w_ompc(1-l,j,k,1) = 2._rkind*w_ompc(2-l,j,k,1)-w_ompc(3-l,j,k,1)
          enddo
        enddo
      enddo
    endif
    if (irightx < 0) then
      !$omp parallel do default(firstprivate) shared(w_ompc, &
      !$omp& wall_tag_ompc)
      do k = 1,nz
        do j = 1,ny
          do l=1,ng
            w_ompc(nx+l,j,k,1) = 2._rkind*w_ompc(nx+l-1,j,k,1)-w_ompc(nx+l-2,j,k,1)
          enddo
        enddo
      enddo
    endif
    !$omp parallel do default(firstprivate) shared(w_ompc, &
    !$omp& wall_tag_ompc)
    do k = 1,nz
      do i = 1,nx
        do l=1,ng
          if(wall_tag_ompc(i) < 1) then
            w_ompc(i,1-l,k,1) = 2._rkind*w_ompc(i,2-l,k,1)-w_ompc(i,3-l,k,1)
          endif
        enddo
      enddo
    enddo

    !$omp parallel do default(firstprivate) shared(w_ompc, &
    !$omp& wall_tag_ompc)
    do k = 1,nz
      do i = 1,nx
        do l=1,ng
          w_ompc(i,ny+l,k,1) = 2._rkind*w_ompc(i,ny+l-1,k,1)-w_ompc(i,ny+l-2,k,1)
        enddo
      enddo
    enddo
    if(ndim == 3) then
      if (ileftz < 0) then
        !$omp parallel do default(firstprivate) shared(w_ompc, &
        !$omp& wall_tag_ompc)
        do j = 1,ny
          do i = 1,nx
            do l=1,ng
              w_ompc(i,j,1-l,1) = 2._rkind*w_ompc(i,j,2-l,1)-w_ompc(i,j,3-l,1)
            enddo
          enddo
        enddo
      endif
      if (irightz < 0) then
        !$omp parallel do default(firstprivate) shared(w_ompc, &
        !$omp& wall_tag_ompc)
        do j = 1,ny
          do i = 1,nx
            do l=1,ng
              w_ompc(i,j,nz+l,1) = 2._rkind*w_ompc(i,j,nz+l-1,1)-w_ompc(i,j,nz+l-2,1)
            enddo
          enddo
        enddo
      endif
    endif

  endsubroutine bcextr_airfoil_var_subroutine


  subroutine bcextr_subroutine(ilat,nx,ny,nz,ng,nv,w_ompc)
    integer :: nx,ny,nz,ng,nv
    integer :: ilat
    integer :: i,j,k,l,m,iercuda
    real(rkind), dimension(1-ng:, 1-ng:, 1-ng:, 1:), intent(inout) :: w_ompc
    if (ilat==1) then
      !$omp parallel do default(firstprivate) shared(w_ompc)
      do k = 1,nz
        do j = 1,ny
          do l=1,ng
            do m=1,nv
              w_ompc(1-l,j,k,m) = w_ompc(1,j,k,m)
            enddo
          enddo
        enddo
      enddo
    elseif (ilat==2) then
      !$omp parallel do default(firstprivate) shared(w_ompc)
      do k = 1,nz
        do j = 1,ny
          do l=1,ng
            do m=1,nv
              w_ompc(nx+l,j,k,m) = w_ompc(nx,j,k,m)
            enddo
          enddo
        enddo
      enddo
    elseif (ilat==3) then
      !$omp parallel do default(firstprivate) shared(w_ompc)
      do k = 1,nz
        do i = 1,nx
          do l=1,ng
            do m=1,nv
              w_ompc(i,1-l,k,m) = w_ompc(i,1,k,m)
            enddo
          enddo
        enddo
      enddo
    elseif (ilat==4) then
      !$omp parallel do default(firstprivate) shared(w_ompc)
      do k = 1,nz
        do i = 1,nx
          do l=1,ng
            do m=1,nv
              w_ompc(i,ny+l,k,m) = w_ompc(i,ny,k,m)
            enddo
          enddo
        enddo
      enddo
    elseif (ilat==5) then
      !$omp parallel do default(firstprivate) shared(w_ompc)
      do j = 1,ny
        do i = 1,nx
          do l=1,ng
            do m=1,nv
              w_ompc(i,j,1-l,m) = w_ompc(i,j,1,m)
            enddo
          enddo
        enddo
      enddo
    elseif (ilat==6) then
      !$omp parallel do default(firstprivate) shared(w_ompc)
      do j = 1,ny
        do i = 1,nx
          do l=1,ng
            do m=1,nv
              w_ompc(i,j,nz+l,m) = w_ompc(i,j,nz,m)
            enddo
          enddo
        enddo
      enddo
    endif

  endsubroutine bcextr_subroutine



  subroutine bcsym_subroutine(ilat,nx,ny,nz,ng,w_ompc)
    integer :: nx,ny,nz,ng
    integer :: ilat
    integer :: i,k,l, iercuda
    real(rkind), dimension(1-ng:, 1-ng:, 1-ng:, 1:), intent(inout) :: w_ompc
    if (ilat==1) then
    elseif (ilat==2) then
    elseif (ilat==3) then
      !$omp parallel do default(firstprivate) shared(w_ompc)
      do k = 1,nz
        do i = 1,nx
          do l=1,ng
            w_ompc(i,1-l,k,1) = w_ompc(i,1+l,k,1)
            w_ompc(i,1-l,k,2) = w_ompc(i,1+l,k,2)
            w_ompc(i,1-l,k,3) = w_ompc(i,1+l,k,3)
            w_ompc(i,1-l,k,4) = w_ompc(i,1+l,k,4)
            w_ompc(i,1-l,k,5) = w_ompc(i,1+l,k,5)
          enddo
        enddo
      enddo
    elseif (ilat==4) then
    elseif (ilat==5) then
    elseif (ilat==6) then
    endif

  endsubroutine bcsym_subroutine



  subroutine bcsym_c2_subroutine(ilat,nx,ny,nz,ng,w_ompc,dxdcsic2_ompc,dydcsic2_ompc)
    integer :: nx,ny,nz,ng
    integer :: ilat
    integer :: i,k,l, iercuda
    real(rkind), dimension(1-ng:, 1-ng:, 1-ng:, 1:), intent(inout) :: w_ompc
    real(rkind), dimension(1-ng:, 1-ng:), intent(in) :: dxdcsic2_ompc
    real(rkind), dimension(1-ng:, 1-ng:), intent(in) :: dydcsic2_ompc
    real(rkind) :: abcsym_i31,abcsym_i32
    real(rkind) :: abcsym_i41,abcsym_i42
    real(rkind) :: tauvers_i31,tauvers_i32
    real(rkind) :: tauvers_i41,tauvers_i42
    real(rkind) :: taumod
    if (ilat==1) then
    elseif (ilat==2) then
    elseif (ilat==3) then
      !$omp parallel do default(firstprivate) shared(w_ompc, &
      !$omp& dxdcsic2_ompc, &
      !$omp& dydcsic2_ompc)
      do k = 1,nz
        do i = 1,nx
          do l=1,ng
            taumod = sqrt(dxdcsic2_ompc(i,1)**2+dydcsic2_ompc(i,1)**2)
            tauvers_i31 = dxdcsic2_ompc(i,1)/taumod
            tauvers_i32 = dydcsic2_ompc(i,1)/taumod
            abcsym_i31 = tauvers_i31**2 - tauvers_i32**2
            abcsym_i32 = 2._rkind*tauvers_i31*tauvers_i32
            w_ompc(i,1-l,k,1) = w_ompc(i,1+l,k,1)
            w_ompc(i,1-l,k,4) = w_ompc(i,1+l,k,4)
            w_ompc(i,1-l,k,5) = w_ompc(i,1+l,k,5)
            w_ompc(i,1-l,k,2) = abcsym_i31*w_ompc(i,1+l,k,2) + abcsym_i32*w_ompc(i,1+l,k,3)
            w_ompc(i,1-l,k,3) = abcsym_i32*w_ompc(i,1+l,k,2) - abcsym_i31*w_ompc(i,1+l,k,3)
          enddo
        enddo
      enddo
    elseif (ilat==4) then
      !$omp parallel do default(firstprivate) shared(w_ompc, &
      !$omp& dxdcsic2_ompc, &
      !$omp& dydcsic2_ompc)
      do k = 1,nz
        do i = 1,nx
          do l=1,ng
            taumod = sqrt(dxdcsic2_ompc(i,ny)**2+dydcsic2_ompc(i,ny)**2)
            tauvers_i41 = -dxdcsic2_ompc(i,ny)/taumod
            tauvers_i42 = -dydcsic2_ompc(i,ny)/taumod
            abcsym_i41 = tauvers_i41**2 - tauvers_i42**2
            abcsym_i42 = 2._rkind*tauvers_i41*tauvers_i42
            w_ompc(i,ny+l,k,1) = w_ompc(i,ny-l,k,1)
            w_ompc(i,ny+l,k,4) = w_ompc(i,ny-l,k,4)
            w_ompc(i,ny+l,k,5) = w_ompc(i,ny-l,k,5)
            w_ompc(i,ny+l,k,2) = abcsym_i41*w_ompc(i,ny-l,k,2) + abcsym_i42*w_ompc(i,ny-l,k,3)
            w_ompc(i,ny+l,k,3) = abcsym_i42*w_ompc(i,ny-l,k,2) - abcsym_i41*w_ompc(i,ny-l,k,3)
          enddo
        enddo
      enddo
    elseif (ilat==5) then
    elseif (ilat==6) then
    endif

  endsubroutine bcsym_c2_subroutine



  subroutine bcwall_subroutine(ilat,nx,ny,nz,ng,twall,w_ompc,w_aux_ompc,indx_cp_l,indx_cp_r,&
  &cv_coeff_ompc,t0,rgas0,calorically_perfect,tol_iter_nr)
    integer :: nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect
    integer :: ilat
    real(rkind) :: twall, t0, rgas0, tol_iter_nr
    integer :: i,k,l, iercuda
    real(rkind) :: rho,uu,vv,ww,qq,pp,tt,rhoe,ee
    real(rkind), dimension(1-ng:, 1-ng:, 1-ng:, 1:), intent(inout) :: w_ompc
    real(rkind), dimension(1-ng:, 1-ng:, 1-ng:, 1:), intent(in) :: w_aux_ompc
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in) :: cv_coeff_ompc
    if (ilat==1) then
    elseif (ilat==2) then
    elseif (ilat==3) then
      !$omp parallel do default(firstprivate) shared(w_ompc, &
      !$omp& w_aux_ompc, &
      !$omp& cv_coeff_ompc)
      do k = 1,nz
        do i = 1,nx
          w_ompc(i,1,k,2) = 0._rkind
          w_ompc(i,1,k,3) = 0._rkind
          w_ompc(i,1,k,4) = 0._rkind
          ee = get_e_from_temperature_dev(twall, t0, indx_cp_l, indx_cp_r, cv_coeff_ompc, calorically_perfect)
          w_ompc(i,1,k,5) = w_ompc(i,1,k,1)*ee
          do l=1,ng
            rho = w_ompc(i,1+l,k,1)
            uu = w_ompc(i,1+l,k,2)/rho
            vv = w_ompc(i,1+l,k,3)/rho
            ww = w_ompc(i,1+l,k,4)/rho
            rhoe = w_ompc(i,1+l,k,5)
            qq = 0.5_rkind*(uu*uu+vv*vv+ww*ww)
            ee = rhoe/rho-qq
            tt = get_temperature_from_e_dev(ee,w_aux_ompc(i,1+l,k,6),t0,cv_coeff_ompc,indx_cp_l,&
            &indx_cp_r,calorically_perfect,tol_iter_nr)
            pp = rho*tt*rgas0
            tt = 2._rkind*twall-tt
            ee = get_e_from_temperature_dev(tt, t0, indx_cp_l, indx_cp_r, cv_coeff_ompc,calorically_perfect)
            rho = pp/tt/rgas0
            w_ompc(i,1-l,k,1) = rho
            w_ompc(i,1-l,k,2) = -rho*uu
            w_ompc(i,1-l,k,3) = -rho*vv
            w_ompc(i,1-l,k,4) = -rho*ww
            w_ompc(i,1-l,k,5) = rho*(ee+qq)
          enddo
        enddo
      enddo
    elseif (ilat==4) then
    elseif (ilat==5) then
    elseif (ilat==6) then
    endif

  endsubroutine bcwall_subroutine



  subroutine bcwall_airfoil_subroutine(ilat,nx,ny,nz,ng,twall,w_ompc,w_aux_ompc,wall_tag_ompc,&
  &indx_cp_l,indx_cp_r,cv_coeff_ompc,t0,rgas0,calorically_perfect,tol_iter_nr,xc2_ompc,yc2_ompc,a_tw,&
  &v_bs,thic,kx_tw,om_tw,xtw1,xtw2,time,dxdetanc2_ompc,dydetanc2_ompc,u0)
    integer :: nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect
    integer :: ilat
    integer :: i,k,l, iercuda
    integer, dimension(1-ng:nx+ng), intent(in) :: wall_tag_ompc
    real(rkind) :: u0, twall, t0, rgas0, tol_iter_nr
    real(rkind) :: rho,uu,vv,ww,qq,pp,tt,rhoe,ee
    real(rkind), dimension(1-ng:, 1-ng:, 1-ng:, 1:), intent(inout) :: w_ompc
    real(rkind), dimension(1-ng:, 1-ng:, 1-ng:, 1:), intent(in) :: w_aux_ompc
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in) :: cv_coeff_ompc
    real(rkind), intent(in), dimension(1-ng:nx+ng,1-ng:ny+ng) :: xc2_ompc, yc2_ompc
    real(rkind), intent(in), dimension(1-ng:, 1-ng:) :: dxdetanc2_ompc, dydetanc2_ompc
    real(rkind), intent(in) :: a_tw, v_bs, thic, kx_tw, om_tw, time, xtw1, xtw2
    real(rkind) :: dx1, dx2, ftanh, tww, ug, vg, wg
    if (ilat==1) then
    elseif (ilat==2) then
    elseif (ilat==3) then
      !$omp parallel do default(firstprivate) shared(w_ompc, &
      !$omp& w_aux_ompc, &
      !$omp& cv_coeff_ompc, &
      !$omp& xc2_ompc, &
      !$omp& yc2_ompc, &
      !$omp& dxdetanc2_ompc, &
      !$omp& dydetanc2_ompc, &
      !$omp& wall_tag_ompc)
      do k = 1,nz
        do i = 1,nx
          if(wall_tag_ompc(i) < 1) then
            dx1 = (xc2_ompc(i,1)-xtw1)
            dx2 = (xc2_ompc(i,1)-xtw2)
            if(a_tw > 0._rkind .and. wall_tag_ompc(i) == 0) then
              ftanh = 0.5_rkind*(tanh(dx1/thic)-tanh(dx2/thic))
              tww = ftanh*(a_tw*sin(kx_tw*xc2_ompc(i,1)-om_tw*time))
              w_ompc(i,1,k,2) = 0._rkind
              w_ompc(i,1,k,3) = 0._rkind
              w_ompc(i,1,k,4) = w_ompc(i,1,k,1)*tww
              ee = get_e_from_temperature_dev(twall, t0, indx_cp_l, indx_cp_r, cv_coeff_ompc, calorically_perfect)
              w_ompc(i,1,k,5) = w_ompc(i,1,k,1)*ee+w_ompc(i,1,k,1)*0.5_rkind*(tww*tww)
            elseif(abs(v_bs) > 0._rkind .and. wall_tag_ompc(i) == 0) then
              ftanh = 0.5_rkind*(tanh(dx1/thic)-tanh(dx2/thic))
              tww = ftanh*v_bs*u0
              w_ompc(i,1,k,2) = w_ompc(i,1,k,1)*tww*dxdetanc2_ompc(i,1)
              w_ompc(i,1,k,3) = w_ompc(i,1,k,1)*tww*dydetanc2_ompc(i,1)
              w_ompc(i,1,k,4) = 0._rkind
              ee = get_e_from_temperature_dev(twall, t0, indx_cp_l, indx_cp_r, cv_coeff_ompc, calorically_perfect)
              w_ompc(i,1,k,5) = w_ompc(i,1,k,1)*ee+w_ompc(i,1,k,1)*0.5_rkind*(tww*tww)
            else
              w_ompc(i,1,k,2) = 0._rkind
              w_ompc(i,1,k,3) = 0._rkind
              w_ompc(i,1,k,4) = 0._rkind
              ee = get_e_from_temperature_dev(twall, t0, indx_cp_l, indx_cp_r, cv_coeff_ompc, calorically_perfect)
              w_ompc(i,1,k,5) = w_ompc(i,1,k,1)*ee
            endif
            do l=1,ng
              rho = w_ompc(i,1+l,k,1)
              uu = w_ompc(i,1+l,k,2)/rho
              vv = w_ompc(i,1+l,k,3)/rho
              ww = w_ompc(i,1+l,k,4)/rho
              rhoe = w_ompc(i,1+l,k,5)
              qq = 0.5_rkind*(uu*uu+vv*vv+ww*ww)
              ee = rhoe/rho-qq
              tt = get_temperature_from_e_dev(ee,w_aux_ompc(i,1+l,k,6),t0,cv_coeff_ompc,indx_cp_l,&
              &indx_cp_r,calorically_perfect,tol_iter_nr)
              pp = rho*tt*rgas0
              tt = 2._rkind*twall-tt
              ee = get_e_from_temperature_dev(tt, t0, indx_cp_l, indx_cp_r, cv_coeff_ompc,calorically_perfect)
              rho = pp/tt/rgas0
              ug = 2._rkind*w_ompc(i,1,k,2)/w_ompc(i,1,k,1) - uu
              vg = 2._rkind*w_ompc(i,1,k,3)/w_ompc(i,1,k,1) - vv
              wg = 2._rkind*w_ompc(i,1,k,4)/w_ompc(i,1,k,1) - ww
              qq = 0.5_rkind*(ug*ug+vg*vg+wg*wg)
              w_ompc(i,1-l,k,1) = rho
              w_ompc(i,1-l,k,2) = rho*ug
              w_ompc(i,1-l,k,3) = rho*vg
              w_ompc(i,1-l,k,4) = rho*wg
              w_ompc(i,1-l,k,5) = rho*(ee+qq)
            enddo
          endif
        enddo
      enddo
    elseif (ilat==4) then
    elseif (ilat==5) then
    elseif (ilat==6) then
    endif

  endsubroutine bcwall_airfoil_subroutine



  subroutine bcwall_staggered_subroutine(ilat,nx,ny,nz,ng,twall,w_ompc,w_aux_ompc,indx_cp_l,&
  &indx_cp_r,cv_coeff_ompc,t0,rgas0,calorically_perfect,tol_iter_nr)
    integer :: nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect
    integer :: ilat
    real(rkind) :: twall, t0, rgas0, tol_iter_nr
    integer :: i,k,l, iercuda
    real(rkind) :: rho,uu,vv,ww,qq,pp,tt,rhoe,ee
    real(rkind), dimension(1-ng:, 1-ng:, 1-ng:, 1:), intent(inout) :: w_ompc
    real(rkind), dimension(1-ng:, 1-ng:, 1-ng:, 1:), intent(in) :: w_aux_ompc
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in) :: cv_coeff_ompc
    if (ilat==1) then
    elseif (ilat==2) then
    elseif (ilat==3) then
      !$omp parallel do default(firstprivate) shared(w_ompc, &
      !$omp& w_aux_ompc, &
      !$omp& cv_coeff_ompc)
      do k = 1,nz
        do l = 1,ng
          do i = 1,nx
            rho = w_ompc(i,l,k,1)
            uu = w_ompc(i,l,k,2)/rho
            vv = w_ompc(i,l,k,3)/rho
            ww = w_ompc(i,l,k,4)/rho
            rhoe = w_ompc(i,l,k,5)
            qq = 0.5_rkind*(uu*uu+vv*vv+ww*ww)
            ee = rhoe/rho-qq
            tt = get_temperature_from_e_dev(ee,w_aux_ompc(i,l,k,6),t0,cv_coeff_ompc,indx_cp_l,&
            &indx_cp_r,calorically_perfect,tol_iter_nr)
            pp = rho*tt*rgas0
            tt = 2._rkind*twall-tt
            ee = get_e_from_temperature_dev(tt, t0, indx_cp_l, indx_cp_r, cv_coeff_ompc,calorically_perfect)
            rho = pp/tt/rgas0
            w_ompc(i,1-l,k,1) = rho
            w_ompc(i,1-l,k,2) = -rho*uu
            w_ompc(i,1-l,k,3) = -rho*vv
            w_ompc(i,1-l,k,4) = -rho*ww
            w_ompc(i,1-l,k,5) = rho*(ee+qq)
          enddo
        enddo
      enddo
    elseif (ilat==4) then
      !$omp parallel do default(firstprivate) shared(w_ompc, &
      !$omp& w_aux_ompc, &
      !$omp& cv_coeff_ompc)
      do k = 1,nz
        do l = 1,ng
          do i = 1,nx
            rho = w_ompc(i,ny+1-l,k,1)
            uu = w_ompc(i,ny+1-l,k,2)/rho
            vv = w_ompc(i,ny+1-l,k,3)/rho
            ww = w_ompc(i,ny+1-l,k,4)/rho
            rhoe = w_ompc(i,ny+1-l,k,5)
            qq = 0.5_rkind*(uu*uu+vv*vv+ww*ww)
            ee = rhoe/rho-qq
            tt = get_temperature_from_e_dev(ee,w_aux_ompc(i,ny+1-l,k,6),t0,cv_coeff_ompc,indx_cp_l,&
            &indx_cp_r,calorically_perfect,tol_iter_nr)
            pp = rho*tt*rgas0
            tt = 2._rkind*twall-tt
            ee = get_e_from_temperature_dev(tt, t0, indx_cp_l, indx_cp_r, cv_coeff_ompc,calorically_perfect)
            rho = pp/tt/rgas0
            w_ompc(i,ny+l,k,1) = rho
            w_ompc(i,ny+l,k,2) = -rho*uu
            w_ompc(i,ny+l,k,3) = -rho*vv
            w_ompc(i,ny+l,k,4) = -rho*ww
            w_ompc(i,ny+l,k,5) = rho*(ee+qq)
          enddo
        enddo
      enddo
    elseif (ilat==5) then
    elseif (ilat==6) then
    endif

  endsubroutine bcwall_staggered_subroutine



  subroutine compute_residual_subroutine(nx,ny,nz,ng,nv,fln_ompc,dt,residual_rhou,fluid_mask_ompc)
    integer :: nx, ny, nz, ng, nv
    real(rkind), intent(out) :: residual_rhou
    real(rkind), intent(in) :: dt
    real(rkind), dimension(1:nx, 1:ny, 1:nz, nv), intent(in) :: fln_ompc
    integer, dimension(1-ng:,1-ng:,1-ng:), intent(in) :: fluid_mask_ompc
    integer :: i,j,k,iercuda
    residual_rhou = 0._rkind
    !$omp parallel do default(firstprivate) shared(fln_ompc, &
    !$omp& fluid_mask_ompc) &
    !$omp& reduction(+:residual_rhou)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          if (fluid_mask_ompc(i,j,k)==0) then
            residual_rhou = residual_rhou + (fln_ompc(i,j,k,2)/dt)**2
          endif
        enddo
      enddo
    enddo


  endsubroutine compute_residual_subroutine



  subroutine compute_airfoil_forces_runtime_subroutine(nx,ny,nz,ng,nv,p0,u0,rgas0,w_aux_ompc,&
  &meta_ompc,csimod_ompc,wall_tag_ompc,dxdcsic2_ompc,dydcsic2_ompc,n,a,pn,pa,tn,ta)
    integer :: nx, ny, nz, ng, nv
    real(rkind) :: p0, u0, rgas0
    real(rkind) :: n, a, pn, pa, tn, ta
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:) :: w_aux_ompc
    real(rkind), dimension(1-ng:,1-ng:), intent(in) :: meta_ompc, csimod_ompc
    integer, dimension(1-ng:nx+ng), intent(in) :: wall_tag_ompc
    real(rkind), dimension(1-ng:,1-ng:), intent(in) :: dxdcsic2_ompc, dydcsic2_ompc
    real(rkind) :: dudy, dudyw, tauw, pw, cf, cp, pwf, tauwf, ds, costh, sinth
    real(rkind) :: dn, da, dpn, dpa, dtn, dta, ut1, ut2, ut3, ut4
    integer :: i,j,k,iercuda
    n = 0._rkind
    a = 0._rkind
    pn = 0._rkind
    pa = 0._rkind
    tn = 0._rkind
    ta = 0._rkind
    !$omp parallel do default(firstprivate) shared(w_aux_ompc, &
    !$omp& meta_ompc, &
    !$omp& csimod_ompc, &
    !$omp& dxdcsic2_ompc, &
    !$omp& dydcsic2_ompc, &
    !$omp& wall_tag_ompc) &
    !$omp& reduction(+:n,a,pn,pa,tn,ta)
    do k = 1,nz
      do j = 1,1
        do i = 1,nx
          if(wall_tag_ompc(i) < 1) then
            ds = csimod_ompc(i,1)
            costh = dxdcsic2_ompc(i,1)/ds
            sinth = dydcsic2_ompc(i,1)/ds
            ut1 = w_aux_ompc(i,1,k,2)*costh+w_aux_ompc(i,1,k,3)*sinth
            ut2 = w_aux_ompc(i,2,k,2)*costh+w_aux_ompc(i,2,k,3)*sinth
            ut3 = w_aux_ompc(i,3,k,2)*costh+w_aux_ompc(i,3,k,3)*sinth
            ut4 = w_aux_ompc(i,4,k,2)*costh+w_aux_ompc(i,4,k,3)*sinth
            dudy = -22._rkind*ut1 + 36._rkind*ut2 - 18._rkind*ut3 + 4._rkind*ut4
            dudyw = dudy*meta_ompc(i,1)/12._rkind
            tauw = w_aux_ompc(i,1,k,7)*dudyw
            pw = w_aux_ompc(i,1,k,1)*w_aux_ompc(i,1,k,6)*rgas0
            cf = tauw/(0.5_rkind*u0*u0)
            cp = (pw-p0)/(0.5_rkind*u0*u0)
            pwf = pw
            tauwf = tauw
            dn = -costh*pwf+sinth*tauwf
            da = +sinth*pwf+costh*tauwf
            dpn = -costh*pwf
            dpa = +sinth*pwf
            dtn = sinth*tauwf
            dta = costh*tauwf
            if(wall_tag_ompc(i) < -1) ds = ds/2._rkind
            n = n + dn*ds
            a = a + da*ds
            pn = pn + dpn*ds
            pa = pa + dpa*ds
            tn = tn + dtn*ds
            ta = ta + dta*ds
          endif
        enddo
      enddo
    enddo


  endsubroutine compute_airfoil_forces_runtime_subroutine



  subroutine compute_rho_t_p_minmax_subroutine(nx,ny,nz,ng,rgas0,w_aux_ompc,rhomin,rhomax,tmin,tmax,pmin,pmax,fluid_mask_ompc)
    integer, intent(in) :: nx, ny, nz, ng
    real(rkind), intent(in) :: rgas0
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(in) :: w_aux_ompc
    integer, dimension(1-ng:,1-ng:,1-ng:), intent(in) :: fluid_mask_ompc
    real(rkind), intent(out) :: rhomin, rhomax, tmin, tmax, pmin, pmax
    integer :: i,j,k,iercuda
    real(rkind) :: rho,tt,pp
    rhomin = huge(1._rkind)
    rhomax = -100._rkind
    tmin = huge(1._rkind)
    tmax = -100._rkind
    pmin = huge(1._rkind)
    pmax = -100._rkind
    !$omp parallel do default(firstprivate) shared(w_aux_ompc, &
    !$omp& fluid_mask_ompc) &
    !$omp& reduction(min:rhomin,tmin,pmin)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          if (fluid_mask_ompc(i,j,k)==0) then
            rho = w_aux_ompc(i,j,k,1)
            tt = w_aux_ompc(i,j,k,6)
            pp = rho*tt*rgas0
            rhomin = min(rhomin,rho)
            tmin = min(tmin ,tt )
            pmin = min(pmin ,pp )
          endif
        enddo
      enddo
    enddo

    !$omp parallel do default(firstprivate) shared(w_aux_ompc, &
    !$omp& fluid_mask_ompc) &
    !$omp& reduction(max:rhomax,tmax,pmax)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          if (fluid_mask_ompc(i,j,k)==0) then
            rho = w_aux_ompc(i,j,k,1)
            tt = w_aux_ompc(i,j,k,6)
            pp = rho*tt*rgas0
            rhomax = max(rhomax,rho)
            tmax = max(tmax ,tt )
            pmax = max(pmax ,pp )
          endif
        enddo
      enddo
    enddo


  endsubroutine compute_rho_t_p_minmax_subroutine



  subroutine compute_dt_subroutine(nx,ny,nz,ng,rgas0,prandtl,dcsidx_ompc,detady_ompc,dzitdz_ompc,&
  &dcsidxs_ompc,detadys_ompc,dzitdzs_ompc,w_ompc,w_aux_ompc,dtxi_max,dtyi_max,dtzi_max,dtxv_max,&
  &dtyv_max,dtzv_max,dtxk_max,dtyk_max,dtzk_max,indx_cp_l,indx_cp_r,cp_coeff_ompc,fluid_mask_ompc,&
  &calorically_perfect,t0)
    integer :: nx, ny, nz, ng, indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind) :: rgas0, t0
    real(rkind) :: prandtl
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:) , intent(in) :: w_ompc
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(in) :: w_aux_ompc
    integer, dimension(1-ng:,1-ng:,1-ng:), intent(in) :: fluid_mask_ompc
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in) :: cp_coeff_ompc
    real(rkind), dimension(1:), intent(in) :: dcsidx_ompc, detady_ompc, dzitdz_ompc
    real(rkind), dimension(1:), intent(in) :: dcsidxs_ompc, detadys_ompc, dzitdzs_ompc
    real(rkind) :: dtxi, dtyi, dtzi, dtxv, dtyv, dtzv, dtxk, dtyk, dtzk
    integer :: i,j,k,ll,iercuda
    real(rkind) :: rho, ri, uu, vv, ww, tt, mu, nu, k_over_rhocp, c
    real(rkind) :: dtxi_max, dtyi_max, dtzi_max, dtxv_max, dtyv_max, dtzv_max, dtxk_max, dtyk_max, dtzk_max
    real(rkind) :: gamloc
    dtxi_max = 0._rkind
    dtyi_max = 0._rkind
    dtzi_max = 0._rkind
    dtxv_max = 0._rkind
    dtyv_max = 0._rkind
    dtzv_max = 0._rkind
    dtxk_max = 0._rkind
    dtyk_max = 0._rkind
    dtzk_max = 0._rkind
    !$omp parallel do default(firstprivate) shared(w_ompc, &
    !$omp& w_aux_ompc, &
    !$omp& cp_coeff_ompc, &
    !$omp& dcsidx_ompc, &
    !$omp& detady_ompc, &
    !$omp& dzitdz_ompc, &
    !$omp& dcsidxs_ompc, &
    !$omp& detadys_ompc, &
    !$omp& dzitdzs_ompc, &
    !$omp& fluid_mask_ompc) &
    !$omp& reduction(max:dtxi_max,dtyi_max,dtzi_max,dtxv_max,dtyv_max,dtzv_max,dtxk_max,dtyk_max,dtzk_max)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          if (fluid_mask_ompc(i,j,k)==0) then
            rho = w_ompc(i,j,k,1)
            ri = 1._rkind/rho
            uu = w_aux_ompc(i,j,k,2)
            vv = w_aux_ompc(i,j,k,3)
            ww = w_aux_ompc(i,j,k,4)
            tt = w_aux_ompc(i,j,k,6)
            mu = w_aux_ompc(i,j,k,7)
            nu = ri*mu
            k_over_rhocp = nu/prandtl
            gamloc = get_gamloc_dev(indx_cp_l,indx_cp_r,tt,t0,cp_coeff_ompc,calorically_perfect,rgas0)
            c = sqrt (gamloc*rgas0*tt)
            dtxi = (abs(uu)+c)*dcsidx_ompc(i)
            dtyi = (abs(vv)+c)*detady_ompc(j)
            dtzi = (abs(ww)+c)*dzitdz_ompc(k)
            dtxv = nu*dcsidxs_ompc(i)
            dtyv = nu*detadys_ompc(j)
            dtzv = nu*dzitdzs_ompc(k)
            dtxk = k_over_rhocp*dcsidxs_ompc(i)
            dtyk = k_over_rhocp*detadys_ompc(j)
            dtzk = k_over_rhocp*dzitdzs_ompc(k)
            dtxi_max = max(dtxi_max, dtxi)
            dtyi_max = max(dtyi_max, dtyi)
            dtzi_max = max(dtzi_max, dtzi)
            dtxv_max = max(dtxv_max, dtxv)
            dtyv_max = max(dtyv_max, dtyv)
            dtzv_max = max(dtzv_max, dtzv)
            dtxk_max = max(dtxk_max, dtxk)
            dtyk_max = max(dtyk_max, dtyk)
            dtzk_max = max(dtzk_max, dtzk)
          endif
        enddo
      enddo
    enddo


  endsubroutine compute_dt_subroutine



  subroutine compute_dt_c2_subroutine(nx,ny,nz,ng,rgas0,prandtl,dcsidxnc2_ompc,dcsidync2_ompc,&
  &detadxnc2_ompc,detadync2_ompc,mcsi_ompc,meta_ompc,dzitdz_ompc,dzitdzs_ompc,w_ompc,w_aux_ompc,&
  &dtxi_max,dtyi_max,dtzi_max,dtxv_max,dtyv_max,dtzv_max,dtxk_max,dtyk_max,dtzk_max,indx_cp_l,&
  &indx_cp_r,cp_coeff_ompc,fluid_mask_ompc,calorically_perfect,t0)
    integer :: nx, ny, nz, ng, indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind) :: rgas0, t0
    real(rkind) :: prandtl
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:) , intent(in) :: w_ompc
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(in) :: w_aux_ompc
    integer, dimension(1-ng:,1-ng:,1-ng:), intent(in) :: fluid_mask_ompc
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in) :: cp_coeff_ompc
    real(rkind), dimension(1:), intent(in) :: dzitdz_ompc, dzitdzs_ompc
    real(rkind), dimension(1-ng:,1-ng:), intent(in) :: dcsidxnc2_ompc, dcsidync2_ompc, detadxnc2_ompc, detadync2_ompc
    real(rkind), dimension(1-ng:,1-ng:), intent(in) :: mcsi_ompc, meta_ompc
    real(rkind) :: dtxi, dtyi, dtzi, dtxv, dtyv, dtzv, dtxk, dtyk, dtzk
    integer :: i,j,k,ll,iercuda
    real(rkind) :: rho, ri, uu, vv, ww, tt, mu, nu, k_over_rhocp, c
    real(rkind) :: dtxi_max, dtyi_max, dtzi_max, dtxv_max, dtyv_max, dtzv_max, dtxk_max, dtyk_max, dtzk_max
    real(rkind) :: gamloc, csis1, etas1, util, vtil
    dtxi_max = 0._rkind
    dtyi_max = 0._rkind
    dtzi_max = 0._rkind
    dtxv_max = 0._rkind
    dtyv_max = 0._rkind
    dtzv_max = 0._rkind
    dtxk_max = 0._rkind
    dtyk_max = 0._rkind
    dtzk_max = 0._rkind
    !$omp parallel do default(firstprivate) shared(w_ompc, &
    !$omp& w_aux_ompc, &
    !$omp& cp_coeff_ompc, &
    !$omp& dzitdz_ompc, &
    !$omp& dzitdzs_ompc, &
    !$omp& dcsidxnc2_ompc, &
    !$omp& dcsidync2_ompc, &
    !$omp& detadxnc2_ompc, &
    !$omp& detadync2_ompc, &
    !$omp& mcsi_ompc, &
    !$omp& meta_ompc, &
    !$omp& fluid_mask_ompc) &
    !$omp& reduction(max:dtxi_max,dtyi_max,dtzi_max,dtxv_max,dtyv_max,dtzv_max,dtxk_max,dtyk_max,dtzk_max)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          if (fluid_mask_ompc(i,j,k)==0) then
            rho = w_ompc(i,j,k,1)
            ri = 1._rkind/rho
            uu = w_aux_ompc(i,j,k,2)
            vv = w_aux_ompc(i,j,k,3)
            ww = w_aux_ompc(i,j,k,4)
            tt = w_aux_ompc(i,j,k,6)
            mu = w_aux_ompc(i,j,k,7)
            nu = ri*mu
            k_over_rhocp = nu/prandtl
            gamloc = get_gamloc_dev(indx_cp_l,indx_cp_r,tt,t0,cp_coeff_ompc,calorically_perfect,rgas0)
            k_over_rhocp = nu/prandtl*(gamloc)/(gamloc-1)
            c = sqrt (gamloc*rgas0*tt)
            csis1= mcsi_ompc(i,j)**2
            etas1= meta_ompc(i,j)**2
            util = uu*dcsidxnc2_ompc(i,j)+vv*dcsidync2_ompc(i,j)
            vtil = uu*detadxnc2_ompc(i,j)+vv*detadync2_ompc(i,j)
            dtxi = (abs(util)+c)*mcsi_ompc(i,j)
            dtyi = (abs(vtil)+c)*meta_ompc(i,j)
            dtzi = (abs(ww) +c)*dzitdz_ompc(k)
            dtxv = nu*csis1
            dtyv = nu*etas1
            dtzv = nu*dzitdzs_ompc(k)
            dtxk = k_over_rhocp*csis1
            dtyk = k_over_rhocp*etas1
            dtzk = k_over_rhocp*dzitdzs_ompc(k)
            dtxi_max = max(dtxi_max, dtxi)
            dtyi_max = max(dtyi_max, dtyi)
            dtzi_max = max(dtzi_max, dtzi)
            dtxv_max = max(dtxv_max, dtxv)
            dtyv_max = max(dtyv_max, dtyv)
            dtzv_max = max(dtzv_max, dtzv)
            dtxk_max = max(dtxk_max, dtxk)
            dtyk_max = max(dtyk_max, dtyk)
            dtzk_max = max(dtzk_max, dtzk)
          endif
        enddo
      enddo
    enddo


  endsubroutine compute_dt_c2_subroutine



  subroutine eval_aux_subroutine(nx,ny,nz,ng,istart,iend,jstart,jend,kstart,kend,w_ompc,w_aux_ompc,&
  &visc_model,mu0,t0,sutherland_s,t_ref_dim,powerlaw_vtexp,visc_power,visc_sutherland,visc_no,&
  &cv_coeff_ompc,indx_cp_l,indx_cp_r,rgas0,calorically_perfect,tol_iter_nr,stream_id)
    integer(ikind), intent(in) :: nx, ny, nz, ng, visc_model
    integer(ikind), intent(in) :: istart, iend, jstart, jend, kstart, kend
    integer(ikind), intent(in) :: visc_power, visc_sutherland, visc_no, indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind), intent(in) :: mu0, t0, sutherland_s, t_ref_dim, powerlaw_vtexp, rgas0, tol_iter_nr
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in) :: cv_coeff_ompc
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(in) :: w_ompc
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(inout) :: w_aux_ompc
    integer, intent(in) :: stream_id
    integer(ikind) :: i, j, k
    real(rkind) :: rho, rhou, rhov, rhow, rhoe, ri, uu, vv, ww, qq, pp, tt, mu, ee
    integer(ikind) :: iercuda

    !$omp parallel do default(firstprivate) shared(cv_coeff_ompc, &
    !$omp& w_ompc, &
    !$omp& w_aux_ompc)
    do k = kstart,kend
      do j = jstart,jend
        do i = istart,iend
          rho = w_ompc(i,j,k,1)
          rhou = w_ompc(i,j,k,2)
          rhov = w_ompc(i,j,k,3)
          rhow = w_ompc(i,j,k,4)
          rhoe = w_ompc(i,j,k,5)
          ri = 1._rkind/rho
          uu = rhou*ri
          vv = rhov*ri
          ww = rhow*ri
          qq = 0.5_rkind*(uu*uu+vv*vv+ww*ww)
          ee = rhoe/rho-qq
          tt = get_temperature_from_e_dev(ee, w_aux_ompc(i,j,k,6), t0, cv_coeff_ompc, indx_cp_l,&
          & indx_cp_r, calorically_perfect, tol_iter_nr)
          pp = rho*tt*rgas0
          w_aux_ompc(i,j,k,1) = rho
          w_aux_ompc(i,j,k,2) = uu
          w_aux_ompc(i,j,k,3) = vv
          w_aux_ompc(i,j,k,4) = ww
          w_aux_ompc(i,j,k,5) = (rhoe+pp)/rho
          w_aux_ompc(i,j,k,6) = tt
          if (visc_model == visc_power) then
            mu = mu0 * (tt/t0)**powerlaw_vtexp
          elseif (visc_model == visc_sutherland) then
            mu = mu0 * (tt/t0)**1.5_rkind * (1._rkind+sutherland_s/t_ref_dim)/(tt/t0 + sutherland_s/t_ref_dim)
          elseif (visc_model == visc_no) then
            mu = 0._rkind
          endif
          w_aux_ompc(i,j,k,7) = mu
        enddo
      enddo
    enddo


  endsubroutine eval_aux_subroutine



  subroutine tripping_pressure_subroutine(nx,ny,nz,ng,nv,i_rank_start,pi,itr1,itr2,x0tr,y0tr,x0ts,&
  &y0ts,lamx,lamy,lamz,lamz1,phiz,phiz1,asl,bt,xc2_ompc,yc2_ompc,z_ompc,w_ompc,fl_ompc,wall_tag_ompc,&
  &dxdetanc2_ompc,dydetanc2_ompc)
    real(rkind),intent(in) :: pi,x0tr,y0tr,x0ts,y0ts,lamx,lamy,lamz,lamz1,phiz,phiz1,asl,bt
    integer, intent(in) :: nx,ny,nz,nv,ng,itr1,itr2,i_rank_start
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in) :: w_ompc
    real(rkind), dimension(1:,1:,1:,1:), intent(inout) :: fl_ompc
    integer, dimension(1-ng:nx+ng), intent(in) :: wall_tag_ompc
    real(rkind), intent(in), dimension(1-ng:, 1-ng:) :: dxdetanc2_ompc, dydetanc2_ompc
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng), intent(in) :: xc2_ompc, yc2_ompc
    real(rkind), dimension(1-ng:), intent(in) :: z_ompc
    real(rkind) :: xx,yy,zz,hzi,hzi1,gzt,fz,fzx,fzy
    integer :: i,j,k,iercuda,ig

    !$omp parallel do default(firstprivate) shared(w_ompc, &
    !$omp& fl_ompc, &
    !$omp& dxdetanc2_ompc, &
    !$omp& dydetanc2_ompc, &
    !$omp& xc2_ompc, &
    !$omp& yc2_ompc, &
    !$omp& z_ompc, &
    !$omp& wall_tag_ompc)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          ig = i + i_rank_start
          if (ig > itr1 .and. ig < itr2) then
            xx = xc2_ompc(i,j)
            yy = yc2_ompc(i,j)
            zz = z_ompc(k)
            hzi = sin(2._rkind*pi*lamz *zz+phiz )
            hzi1 = sin(2._rkind*pi*lamz1*zz+phiz1)
            gzt = asl*((1-bt)*hzi+bt*hzi1)
            fz = gzt*exp(-((xx-x0tr)/lamx)**2-((yy-y0tr)/lamy)**2)
            fzx = 0._rkind
            fzy = fz
            fl_ompc(i,j,k,2) = fl_ompc(i,j,k,2) - fzx*w_ompc(i,j,k,1)
            fl_ompc(i,j,k,3) = fl_ompc(i,j,k,3) - fzy*w_ompc(i,j,k,1)
            fl_ompc(i,j,k,5) = fl_ompc(i,j,k,5) -(fzx*w_ompc(i,j,k,2)+fzy*w_ompc(i,j,k,3))
          endif
        enddo
      enddo
    enddo


  endsubroutine tripping_pressure_subroutine


  subroutine tripping_suction_subroutine(nx,ny,nz,ng,nv,i_rank_start,pi,its1,its2,x0tr,y0tr,x0ts,&
  &y0ts,lamx,lamy,lams,lams1,phis,phis1,asl,bt,xc2_ompc,yc2_ompc,z_ompc,w_ompc,fl_ompc,wall_tag_ompc,&
  &dxdetanc2_ompc,dydetanc2_ompc)
    real(rkind),intent(in) :: pi,x0tr,y0tr,x0ts,y0ts,lamx,lamy,lams,lams1,phis,phis1,asl,bt
    integer, intent(in) :: nx,ny,nz,nv,ng,its1,its2,i_rank_start
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in) :: w_ompc
    real(rkind), dimension(1:,1:,1:,1:), intent(inout) :: fl_ompc
    integer, dimension(1-ng:nx+ng), intent(in) :: wall_tag_ompc
    real(rkind), intent(in), dimension(1-ng:, 1-ng:) :: dxdetanc2_ompc, dydetanc2_ompc
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng), intent(in) :: xc2_ompc, yc2_ompc
    real(rkind), dimension(1-ng:), intent(in) :: z_ompc
    real(rkind) :: xx,yy,zz,hzi,hzi1,gzt,fz,fzx,fzy
    integer :: i,j,k,iercuda,ig

    !$omp parallel do default(firstprivate) shared(w_ompc, &
    !$omp& fl_ompc, &
    !$omp& dxdetanc2_ompc, &
    !$omp& dydetanc2_ompc, &
    !$omp& xc2_ompc, &
    !$omp& yc2_ompc, &
    !$omp& z_ompc, &
    !$omp& wall_tag_ompc)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          ig = i + i_rank_start
          if (ig > its1 .and. ig < its2) then
            xx = xc2_ompc(i,j)
            yy = yc2_ompc(i,j)
            zz = z_ompc(k)
            hzi = sin(2._rkind*pi*lams *zz+phis)
            hzi1 = sin(2._rkind*pi*lams1*zz+phis1)
            gzt = asl*((1-bt)*hzi+bt*hzi1)
            fz = gzt*exp(-((xx-x0ts)/lamx)**2-((yy-y0ts)/lamy)**2)
            fzx = 0._rkind
            fzy = fz
            fl_ompc(i,j,k,2) = fl_ompc(i,j,k,2) - fzx*w_ompc(i,j,k,1)
            fl_ompc(i,j,k,3) = fl_ompc(i,j,k,3) - fzy*w_ompc(i,j,k,1)
            fl_ompc(i,j,k,5) = fl_ompc(i,j,k,5) -(fzx*w_ompc(i,j,k,2)+fzy*w_ompc(i,j,k,3))
          endif
        enddo
      enddo
    enddo


  endsubroutine tripping_suction_subroutine



  function get_gamloc_dev(indx_cp_l,indx_cp_r,tt,t0,cp_coeff_ompc,calorically_perfect,rgas0)
    real(rkind) :: get_gamloc_dev
    integer, value :: indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind), value :: tt, rgas0, t0
    real(rkind), dimension(indx_cp_l:indx_cp_r+1) :: cp_coeff_ompc
    real(rkind) :: cploc, gamloc
    integer :: l

    if (calorically_perfect==1) then
      cploc = cp_coeff_ompc(0)
    else
      cploc = 0._rkind
      do l=indx_cp_l,indx_cp_r
        cploc = cploc+cp_coeff_ompc(l)*(tt/t0)**l
      enddo
    endif
    gamloc = cploc/(cploc-rgas0)
    get_gamloc_dev = gamloc
  endfunction get_gamloc_dev


  function get_temperature_from_e_dev(ee,t_start,t0,cv_coeff_ompc,indx_cp_l,indx_cp_r,calorically_perfect,tol_iter_nr)
    real(rkind) :: get_temperature_from_e_dev
    integer, value :: indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind), value :: ee, t_start, t0, tol_iter_nr
    real(rkind), dimension(indx_cp_l:indx_cp_r+1) :: cv_coeff_ompc
    real(rkind) :: tt, t_old, ebar, den, num, t_pow, t_powp
    integer :: l,iter,max_iter

    max_iter = 50

    ebar = ee - cv_coeff_ompc(indx_cp_r+1)*t0
    if (calorically_perfect==1) then
      tt = t0+ebar/cv_coeff_ompc(0)
    else
      t_old = t_start
      do iter=1,max_iter
        den = 0._rkind
        num = 0._rkind
        do l=indx_cp_l,indx_cp_r
          if (l==-1) then
            t_pow = (t_old/t0)**l
            den = den+cv_coeff_ompc(l)*t_pow
            num = num+cv_coeff_ompc(l)*log(t_old/t0)
          else
            t_pow = (t_old/t0)**l
            t_powp = (t_old/t0)*t_pow
            den = den+cv_coeff_ompc(l)*t_pow
            num = num+cv_coeff_ompc(l)*(t_powp-1._rkind)/(l+1._rkind)
          endif
        enddo
        num = num*t0
        tt = t_old+(ebar-num)/den
        if (abs(tt-t_old) < tol_iter_nr) exit
        t_old = tt
      enddo
    endif
    get_temperature_from_e_dev = tt
  endfunction get_temperature_from_e_dev


  function get_e_from_temperature_dev(tt,t0,indx_cp_l,indx_cp_r,cv_coeff_ompc,calorically_perfect)
    real(rkind) :: get_e_from_temperature_dev
    real(rkind),value :: tt, t0
    integer,value :: indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind), dimension(indx_cp_l:indx_cp_r+1) :: cv_coeff_ompc
    real(rkind) :: ee
    integer :: l

    ee = cv_coeff_ompc(indx_cp_r+1)
    if (calorically_perfect==1) then
      ee = ee+cv_coeff_ompc(0)*(tt/t0-1._rkind)
    else
      do l=indx_cp_l,indx_cp_r
        if (l==-1) then
          ee = ee+cv_coeff_ompc(l)*log(tt/t0)
        else
          ee = ee+cv_coeff_ompc(l)/(l+1._rkind)*((tt/t0)**(l+1)-1._rkind)
        endif
      enddo
    endif
    ee = ee*t0
    get_e_from_temperature_dev = ee
  endfunction get_e_from_temperature_dev

  function vel_law_of_the_wall(yp)
    real(rkind) :: vel_law_of_the_wall
    real(rkind), value :: yp
    real(rkind) :: vkc, vkci, b, blogk, up

    vkc = 0.41_rkind
    vkci = 1._rkind/vkc
    b = 5.25_rkind
    blogk = b-vkci*log(vkc)
    up = vkci*log(1._rkind+vkc*yp)
    up = up+blogk*(1._rkind-exp(-yp/11._rkind)-yp/11._rkind*exp(-yp/3._rkind))
    vel_law_of_the_wall = max(11._rkind,up)
  endfunction vel_law_of_the_wall


  function tem_law_of_the_wall(yp,pr)
    real(rkind) :: tem_law_of_the_wall
    real(rkind), value :: yp,pr
    real(rkind) :: vkt,vkti,big_gam,beta,onethird,tp

    vkt = 1._rkind/2.12_rkind
    vkti = 2.12_rkind
    big_gam = 0.01_rkind*(yp*pr)**4
    big_gam = big_gam/(1._rkind+5._rkind*pr**3*yp)
    onethird = 1._rkind/3._rkind
    beta = vkti*log(pr)+(3.85_rkind*pr**onethird-1.3_rkind)**2

    tp = pr*yp*exp(-big_gam)+exp(-1._rkind/big_gam)*(beta+vkti*log(1+yp))
    tem_law_of_the_wall = tp
  endfunction tem_law_of_the_wall





  subroutine insitu_swirling_kernel(nv,nx,ny,nz,visc_order,ng,npsi,mpsi,dcsidx_ompc,detady_ompc,&
  &dzitdz_ompc,w_aux_ompc,coeff_deriv1_ompc,psi_ompc,u0,x_ompc)
    integer, value :: nv, nx, ny, nz, visc_order, ng, npsi,mpsi
    real(rkind), value :: u0
    real(rkind), dimension(nx) :: dcsidx_ompc
    real(rkind), dimension(ny) :: detady_ompc
    real(rkind), dimension(nz) :: dzitdz_ompc
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: psi_ompc
    real(rkind), dimension(1-ng:, 1-ng:,1-ng:,1:) :: w_aux_ompc
    real(rkind), dimension(1:,1:), intent(in) :: coeff_deriv1_ompc
    real(rkind), dimension(1-ng:nx+ng) :: x_ompc
    real(rkind) :: uu,vv,ww
    real(rkind) :: ux,uy,uz
    real(rkind) :: vx,vy,vz
    real(rkind) :: wx,wy,wz
    real(rkind), dimension(3) :: eigr_a, eigi_a
    real(rkind), dimension(3,3) :: astar
    real(rkind) :: ccl, div3l
    real(rkind) :: epsi2, omx, omy, omz, omod2, div2, div
    integer :: i,j,k,l

    !$omp parallel do default(firstprivate) shared(dcsidx_ompc, &
    !$omp& detady_ompc, &
    !$omp& dzitdz_ompc, &
    !$omp& psi_ompc, &
    !$omp& w_aux_ompc, &
    !$omp& coeff_deriv1_ompc, &
    !$omp& x_ompc)
    do k = 1,nz
      do i = 1,nx
        do j=1,ny
          uu = w_aux_ompc(i,j,k,2)
          vv = w_aux_ompc(i,j,k,3)
          ww = w_aux_ompc(i,j,k,4)
          ux = 0._rkind
          vx = 0._rkind
          wx = 0._rkind
          uy = 0._rkind
          vy = 0._rkind
          wy = 0._rkind
          uz = 0._rkind
          vz = 0._rkind
          wz = 0._rkind
          do l=1,visc_order/2
            ccl = coeff_deriv1_ompc(l,visc_order/2)
            ux = ux+ccl*(w_aux_ompc(i+l,j,k,2)-w_aux_ompc(i-l,j,k,2))
            vx = vx+ccl*(w_aux_ompc(i+l,j,k,3)-w_aux_ompc(i-l,j,k,3))
            wx = wx+ccl*(w_aux_ompc(i+l,j,k,4)-w_aux_ompc(i-l,j,k,4))
            uy = uy+ccl*(w_aux_ompc(i,j+l,k,2)-w_aux_ompc(i,j-l,k,2))
            vy = vy+ccl*(w_aux_ompc(i,j+l,k,3)-w_aux_ompc(i,j-l,k,3))
            wy = wy+ccl*(w_aux_ompc(i,j+l,k,4)-w_aux_ompc(i,j-l,k,4))
            uz = uz+ccl*(w_aux_ompc(i,j,k+l,2)-w_aux_ompc(i,j,k-l,2))
            vz = vz+ccl*(w_aux_ompc(i,j,k+l,3)-w_aux_ompc(i,j,k-l,3))
            wz = wz+ccl*(w_aux_ompc(i,j,k+l,4)-w_aux_ompc(i,j,k-l,4))
          enddo
          ux = ux*dcsidx_ompc(i)
          vx = vx*dcsidx_ompc(i)
          wx = wx*dcsidx_ompc(i)
          uy = uy*detady_ompc(j)
          vy = vy*detady_ompc(j)
          wy = wy*detady_ompc(j)
          uz = uz*dzitdz_ompc(k)
          vz = vz*dzitdz_ompc(k)
          wz = wz*dzitdz_ompc(k)
          div = ux+vy+wz
          div3l = div/3._rkind
          epsi2 = u0**2
          omz = vx-uy
          omx = wy-vz
          omy = uz-wx
          omod2 = omx*omx+omy*omy+omz*omz
          div2 = div*div
          astar(1,1) = ux-div3l
          astar(1,2) = uy
          astar(1,3) = uz
          astar(2,1) = vx
          astar(2,2) = vy-div3l
          astar(2,3) = vz
          astar(3,1) = wx
          astar(3,2) = wy
          astar(3,3) = wz-div3l
          call eigs33(astar,eigr_a,eigi_a)
          psi_ompc(i,j,k,mpsi) = 2._rkind*max(0._rkind,eigi_a(2))
        enddo
      enddo
    enddo


  endsubroutine insitu_swirling_kernel


  subroutine insitu_swirling_c2_kernel(nv,nx,ny,nz,visc_order,ng,npsi,mpsi,dzitdz_ompc,w_aux_ompc,&
  &coeff_deriv1_ompc,dcsidxc2_ompc,detadyc2_ompc,detadxc2_ompc,dcsidyc2_ompc,psi_ompc,u0,x_ompc,&
  &vis_tag_ompc,wall_tag_ompc)
    integer, value :: nv, nx, ny, nz, visc_order, ng, npsi,mpsi
    real(rkind), value :: u0
    real(rkind), dimension(nz) :: dzitdz_ompc
    real(rkind), dimension(1-ng:,1-ng:), intent(in) :: dcsidxc2_ompc, detadyc2_ompc
    real(rkind), dimension(1-ng:,1-ng:), intent(in) :: detadxc2_ompc, dcsidyc2_ompc
    integer, dimension(1-ng:), intent(in) :: vis_tag_ompc, wall_tag_ompc
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: psi_ompc
    real(rkind), dimension(1-ng:, 1-ng:,1-ng:,1:) :: w_aux_ompc
    real(rkind), dimension(1:,1:), intent(in) :: coeff_deriv1_ompc
    real(rkind), dimension(1-ng:nx+ng) :: x_ompc
    real(rkind) :: uu,vv,ww
    real(rkind) :: ux,uy,uz
    real(rkind) :: vx,vy,vz
    real(rkind) :: wx,wy,wz
    real(rkind) :: ucsi,ueta,uzit,vcsi,veta,vzit,wcsi,weta,wzit
    real(rkind), dimension(3) :: eigr_a, eigi_a
    real(rkind), dimension(3,3) :: astar
    real(rkind) :: ccl, div3l, cli, clj, clk
    real(rkind) :: epsi2, omx, omy, omz, omod2, div2, div
    integer :: i,j,k,l,lmaxi,lmaxj

    !$omp parallel do default(firstprivate) shared(dzitdz_ompc, &
    !$omp& dcsidxc2_ompc, &
    !$omp& detadyc2_ompc, &
    !$omp& detadxc2_ompc, &
    !$omp& dcsidyc2_ompc, &
    !$omp& psi_ompc, &
    !$omp& w_aux_ompc, &
    !$omp& coeff_deriv1_ompc, &
    !$omp& x_ompc, &
    !$omp& vis_tag_ompc, &
    !$omp& wall_tag_ompc)
    do k = 1,nz
      do i = 1,nx
        do j=1,ny
          ucsi = 0._rkind
          vcsi = 0._rkind
          wcsi = 0._rkind
          ueta = 0._rkind
          veta = 0._rkind
          weta = 0._rkind
          uzit = 0._rkind
          vzit = 0._rkind
          wzit = 0._rkind

          lmaxi = visc_order/2
          lmaxj = visc_order/2
          if (j == 1) lmaxi = vis_tag_ompc(i)
          if (wall_tag_ompc(i) < 1) lmaxj = min(j,visc_order/2)

          do l=1,visc_order/2
            cli = coeff_deriv1_ompc(l,lmaxi)
            clj = coeff_deriv1_ompc(l,lmaxj)
            clk = coeff_deriv1_ompc(l,visc_order/2 )

            ucsi = ucsi +cli*(w_aux_ompc(i+l,j,k,2)-w_aux_ompc(i-l,j,k,2))
            vcsi = vcsi +cli*(w_aux_ompc(i+l,j,k,3)-w_aux_ompc(i-l,j,k,3))
            wcsi = wcsi +cli*(w_aux_ompc(i+l,j,k,4)-w_aux_ompc(i-l,j,k,4))

            ueta = ueta +clj*(w_aux_ompc(i,j+l,k,2)-w_aux_ompc(i,j-l,k,2))
            veta = veta +clj*(w_aux_ompc(i,j+l,k,3)-w_aux_ompc(i,j-l,k,3))
            weta = weta +clj*(w_aux_ompc(i,j+l,k,4)-w_aux_ompc(i,j-l,k,4))

            uzit = uzit +clk*(w_aux_ompc(i,j,k+l,2)-w_aux_ompc(i,j,k-l,2))
            vzit = vzit +clk*(w_aux_ompc(i,j,k+l,3)-w_aux_ompc(i,j,k-l,3))
            wzit = wzit +clk*(w_aux_ompc(i,j,k+l,4)-w_aux_ompc(i,j,k-l,4))
          enddo

          ux = ucsi *dcsidxc2_ompc(i,j) + ueta *detadxc2_ompc(i,j)
          vx = vcsi *dcsidxc2_ompc(i,j) + veta *detadxc2_ompc(i,j)
          wx = wcsi *dcsidxc2_ompc(i,j) + weta *detadxc2_ompc(i,j)

          uy = ucsi *dcsidyc2_ompc(i,j) + ueta *detadyc2_ompc(i,j)
          vy = vcsi *dcsidyc2_ompc(i,j) + veta *detadyc2_ompc(i,j)
          wy = wcsi *dcsidyc2_ompc(i,j) + weta *detadyc2_ompc(i,j)
          uz = uzit *dzitdz_ompc(k)
          vz = vzit *dzitdz_ompc(k)
          wz = wzit *dzitdz_ompc(k)
          div = ux+vy+wz
          div3l = div/3._rkind
          epsi2 = u0**2
          omz = vx-uy
          omx = wy-vz
          omy = uz-wx
          omod2 = omx*omx+omy*omy+omz*omz
          div2 = div*div
          astar(1,1) = ux-div3l
          astar(1,2) = uy
          astar(1,3) = uz
          astar(2,1) = vx
          astar(2,2) = vy-div3l
          astar(2,3) = vz
          astar(3,1) = wx
          astar(3,2) = wy
          astar(3,3) = wz-div3l
          call eigs33(astar,eigr_a,eigi_a)
          psi_ompc(i,j,k,mpsi) = 2._rkind*max(0._rkind,eigi_a(2))
        enddo
      enddo
    enddo


  endsubroutine insitu_swirling_c2_kernel


  subroutine insitu_schlieren_kernel(nv,nx,ny,nz,visc_order,ng,npsi,mpsi,dcsidx_ompc,detady_ompc,&
  &dzitdz_ompc,w_aux_ompc,coeff_deriv1_ompc,psi_ompc,u0,x_ompc)
    integer, value :: nv, nx, ny, nz, visc_order, ng, npsi,mpsi
    real(rkind), value :: u0
    real(rkind), dimension(nx) :: dcsidx_ompc
    real(rkind), dimension(ny) :: detady_ompc
    real(rkind), dimension(nz) :: dzitdz_ompc
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: psi_ompc
    real(rkind), dimension(1-ng:, 1-ng:,1-ng:,1:) :: w_aux_ompc
    real(rkind), dimension(1:,1:), intent(in) :: coeff_deriv1_ompc
    real(rkind), dimension(1-ng:nx+ng) :: x_ompc
    real(rkind) :: rhox, rhoy, rhoz
    real(rkind) :: ccl
    integer :: i,j,k,l
    !$omp parallel do default(firstprivate) shared(dcsidx_ompc, &
    !$omp& detady_ompc, &
    !$omp& dzitdz_ompc, &
    !$omp& psi_ompc, &
    !$omp& w_aux_ompc, &
    !$omp& coeff_deriv1_ompc, &
    !$omp& x_ompc)
    do k = 1,nz
      do i = 1,nx
        do j=1,ny
          rhox = 0._rkind
          rhoy = 0._rkind
          rhoz = 0._rkind
          do l=1,visc_order/2
            ccl = coeff_deriv1_ompc(l,visc_order/2)
            rhox = rhox+ccl*(w_aux_ompc(i+l,j,k,1)-w_aux_ompc(i-l,j,k,1))
            rhoy = rhoy+ccl*(w_aux_ompc(i,j+l,k,1)-w_aux_ompc(i,j-l,k,1))
            rhoz = rhoz+ccl*(w_aux_ompc(i,j,k+l,1)-w_aux_ompc(i,j,k-l,1))
          enddo
          rhox = rhox*dcsidx_ompc(i)
          rhoy = rhoy*detady_ompc(j)
          rhoz = rhoz*dzitdz_ompc(k)
          psi_ompc(i,j,k,mpsi) = exp(-sqrt((rhox)**2+(rhoy)**2+(rhoz)**2))
        enddo
      enddo
    enddo

  endsubroutine insitu_schlieren_kernel


  subroutine insitu_schlieren_c2_kernel(nv,nx,ny,nz,visc_order,ng,npsi,mpsi,dzitdz_ompc,w_aux_ompc,&
  &coeff_deriv1_ompc,dcsidxc2_ompc,detadyc2_ompc,detadxc2_ompc,dcsidyc2_ompc,psi_ompc,u0,x_ompc,&
  &vis_tag_ompc,wall_tag_ompc)
    integer, value :: nv, nx, ny, nz, visc_order, ng, npsi,mpsi
    real(rkind), value :: u0
    real(rkind), dimension(nz) :: dzitdz_ompc
    real(rkind), dimension(1-ng:,1-ng:), intent(in) :: dcsidxc2_ompc, detadyc2_ompc
    real(rkind), dimension(1-ng:,1-ng:), intent(in) :: detadxc2_ompc, dcsidyc2_ompc
    integer, dimension(1-ng:), intent(in) :: vis_tag_ompc, wall_tag_ompc
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: psi_ompc
    real(rkind), dimension(1-ng:, 1-ng:,1-ng:,1:) :: w_aux_ompc
    real(rkind), dimension(1:,1:), intent(in) :: coeff_deriv1_ompc
    real(rkind), dimension(1-ng:nx+ng) :: x_ompc
    real(rkind) :: rhocsi,rhoeta,rhozit,rhox,rhoy,rhoz
    real(rkind) :: cli,clj,clk
    integer :: i,j,k,l,lmaxi,lmaxj
    !$omp parallel do default(firstprivate) shared(dzitdz_ompc, &
    !$omp& dcsidxc2_ompc, &
    !$omp& detadyc2_ompc, &
    !$omp& detadxc2_ompc, &
    !$omp& dcsidyc2_ompc, &
    !$omp& psi_ompc, &
    !$omp& w_aux_ompc, &
    !$omp& coeff_deriv1_ompc, &
    !$omp& x_ompc, &
    !$omp& vis_tag_ompc, &
    !$omp& wall_tag_ompc)
    do k = 1,nz
      do i = 1,nx
        do j=1,ny
          rhocsi = 0._rkind
          rhoeta = 0._rkind
          rhozit = 0._rkind

          lmaxi = visc_order/2
          lmaxj = visc_order/2
          if (j == 1) lmaxi = vis_tag_ompc(i)
          if (wall_tag_ompc(i) < 1) lmaxj = min(j,visc_order/2)

          do l=1,visc_order/2
            cli = coeff_deriv1_ompc(l,lmaxi)
            clj = coeff_deriv1_ompc(l,lmaxj)
            clk = coeff_deriv1_ompc(l,visc_order/2 )

            rhocsi = rhocsi +cli*(w_aux_ompc(i+l,j,k,1)-w_aux_ompc(i-l,j,k,1))
            rhoeta = rhoeta +clj*(w_aux_ompc(i,j+l,k,1)-w_aux_ompc(i,j-l,k,1))
            rhozit = rhozit +clk*(w_aux_ompc(i,j,k+l,1)-w_aux_ompc(i,j,k-l,1))
          enddo

          rhox = rhocsi *dcsidxc2_ompc(i,j) + rhoeta *detadxc2_ompc(i,j)
          rhoy = rhocsi *dcsidyc2_ompc(i,j) + rhoeta *detadyc2_ompc(i,j)
          rhoz = rhozit *dzitdz_ompc(k)
          psi_ompc(i,j,k,mpsi) = exp(-sqrt((rhox)**2+(rhoy)**2+(rhoz)**2))
        enddo
      enddo
    enddo

  endsubroutine insitu_schlieren_c2_kernel


  subroutine eigs33(rmat,rex,rimx)
    integer, parameter :: doubtype = real64
    real(rkind), dimension(3) :: rex,rimx
    real(rkind), dimension(3) :: rey,rimy
    real(rkind), dimension(3) :: reu,rimu
    real(rkind), dimension(3) :: rev,rimv
    real(rkind), dimension(3,3) :: rmat,s,ww,temp
    real(doubtype) :: ddbl,pdbl,qdbl
    real(rkind) :: otrd, pi, a, temps, tempw, b, somma, c, p, q, sqdel, sqp, teta
    integer :: ii, jj, k

    otrd = 1./3._rkind
    pi = acos(-1._rkind)
    do ii = 1,3
      do jj = 1,3
        s(ii,jj) = 0.5*(rmat(ii,jj)+rmat(jj,ii))
        ww(ii,jj) = 0.5*(rmat(ii,jj)-rmat(jj,ii))
      enddo
    enddo
    a = -(s(1,1)+s(2,2)+s(3,3))
    temps = 0.
    tempw = 0.
    do ii = 1,3
      do jj = 1,3
        temps = temps + s(ii,jj)*s(ii,jj)
        tempw = tempw + ww(ii,jj)*ww(ii,jj)
      enddo
    enddo
    b = 0.5*(a**2-temps+tempw)
    do ii = 1,3
      do jj = 1,3
        temp(ii,jj)=0.
        do k = 1,3
          temp(ii,jj) = temp(ii,jj)+s(ii,k)*s(k,jj)+3.*ww(ii,k)*ww(k,jj)
        enddo
      enddo
    enddo
    somma=0.
    do ii=1,3
      do jj=1,3
        somma= somma+temp(ii,jj)*s(jj,ii)
      enddo
    enddo
    c = -1./3.*(a**3-3.*a*b+somma)
    pdbl = real(b-a**2/3., doubtype)
    qdbl = real(c-a*b/3.+2*a**3/27.,doubtype)
    ddbl = (qdbl**2)/4.D0 + (pdbl**3)/27.D0
    p = real(pdbl,rkind)
    q = real(qdbl,rkind)
    if(ddbl.gt.0.D0) then
      sqdel = real(sqrt(ddbl), rkind)
      reu(1) =-0.5*q+sqdel
      rev(1) =-0.5*q-sqdel
      reu(1) = sign(1._rkind,reu(1))*(abs(reu(1)))**otrd
      rev(1) = sign(1._rkind,rev(1))*(abs(rev(1)))**otrd
      reu(2) = -0.5*reu(1)
      rev(2) = -0.5*rev(1)
      reu(3) = reu(2)
      rev(3) = rev(2)
      rimu(1) = 0.
      rimv(1) = 0.
      rimu(2) = sqrt(3._rkind)/2.*reu(1)
      rimv(2) = sqrt(3._rkind)/2.*rev(1)
      rimu(3) = -rimu(2)
      rimv(3) = -rimv(2)
      rey(1) = reu(1)+rev(1)
      rimy(1) = rimu(1)+rimv(1)
      rey(2) = reu(2)+rev(3)
      rimy(2) = rimu(2)+rimv(3)
      rey(3) = reu(3)+rev(2)
      rimy(3) = rimu(3)+rimv(2)
    else
      if (q.eq.0.) then
        rey(1) = 0.
        rey(2) = sqrt(-p)
        rey(3) = -rey(2)
      else
        sqp = 2.*sqrt(-p/3.)
        sqdel = real(sqrt(-ddbl), rkind)
        if (q.lt.0.) then
          teta = atan(-2.*sqdel/q)
        else
          teta = pi+atan(-2.*sqdel/q)
        endif
        rey(1) = sqp*cos(teta/3.)
        rey(2) = sqp*cos((teta+2*pi)/3.)
        rey(3) = sqp*cos((teta+4*pi)/3.)
      endif
      rimy(1) = 0.
      rimy(2) = 0.
      rimy(3) = 0.
    endif
    rex(1) = rey(1)-(a/3.)
    rimx(1) = rimy(1)
    rex(2) = rey(2)-(a/3.)
    rimx(2) = rimy(2)
    rex(3) = rey(3)-(a/3.)
    rimx(3) = rimy(3)
    if (rimy(2).lt.0.) then
      rex(2) = rey(3)-(a/3.)
      rimx(2) = rimy(3)
      rex(3) = rey(2)-(a/3.)
      rimx(3) = rimy(2)
    endif
  endsubroutine eigs33


  subroutine insitu_div_subroutine(nx,ny,nz,ng,visc_order,npsi,mpsi,w_aux_ompc,coeff_deriv1_ompc,&
  &dcsidx_ompc,detady_ompc,dzitdz_ompc,psi_ompc)
    integer, intent(in) :: nx, ny, nz, ng, visc_order, npsi, mpsi
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(inout) :: w_aux_ompc
    real(rkind), dimension(1:,1:), intent(in) :: coeff_deriv1_ompc
    real(rkind), dimension(1:), intent(in) :: dcsidx_ompc, detady_ompc, dzitdz_ompc
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(inout) :: psi_ompc
    real(rkind) :: ccl
    real(rkind) :: uu,vv,ww
    real(rkind) :: ux,vy,wz
    real(rkind) :: div
    integer :: lmax
    integer :: i,j,k,l,iercuda
    lmax = visc_order / 2
    !$omp parallel do default(firstprivate) shared(w_aux_ompc, &
    !$omp& coeff_deriv1_ompc, &
    !$omp& dcsidx_ompc, &
    !$omp& detady_ompc, &
    !$omp& dzitdz_ompc, &
    !$omp& psi_ompc)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          uu = w_aux_ompc(i,j,k,2)
          vv = w_aux_ompc(i,j,k,3)
          ww = w_aux_ompc(i,j,k,4)
          ux = 0._rkind
          vy = 0._rkind
          wz = 0._rkind
          do l=1,lmax
            ccl = coeff_deriv1_ompc(l,lmax)
            ux = ux+ccl*(w_aux_ompc(i+l,j,k,2)-w_aux_ompc(i-l,j,k,2))
            vy = vy+ccl*(w_aux_ompc(i,j+l,k,3)-w_aux_ompc(i,j-l,k,3))
            wz = wz+ccl*(w_aux_ompc(i,j,k+l,4)-w_aux_ompc(i,j,k-l,4))
          enddo
          ux = ux*dcsidx_ompc(i)
          vy = vy*detady_ompc(j)
          wz = wz*dzitdz_ompc(k)
          div = ux+vy+wz
          psi_ompc(i,j,k,mpsi) = div
        enddo
      enddo
    enddo


  endsubroutine insitu_div_subroutine



  subroutine insitu_omega_subroutine(nx,ny,nz,ng,visc_order,npsi,mpsi,w_aux_ompc,coeff_deriv1_ompc,&
  &dcsidx_ompc,detady_ompc,dzitdz_ompc,psi_ompc)
    integer, intent(in) :: nx, ny, nz, ng, visc_order, npsi, mpsi
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(inout) :: w_aux_ompc
    real(rkind), dimension(1:,1:), intent(in) :: coeff_deriv1_ompc
    real(rkind), dimension(1:), intent(in) :: dcsidx_ompc, detady_ompc, dzitdz_ompc
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(inout) :: psi_ompc
    integer :: i,j,k,l,iercuda
    real(rkind) :: ccl
    real(rkind) :: uu,vv,ww
    real(rkind) :: uy,uz
    real(rkind) :: vx,vz
    real(rkind) :: wx,wy
    real(rkind) :: omegax, omegay, omegaz, omod2
    integer :: lmax
    lmax = visc_order/2
    !$omp parallel do default(firstprivate) shared(w_aux_ompc, &
    !$omp& coeff_deriv1_ompc, &
    !$omp& dcsidx_ompc, &
    !$omp& detady_ompc, &
    !$omp& dzitdz_ompc, &
    !$omp& psi_ompc)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          uu = w_aux_ompc(i,j,k,2)
          vv = w_aux_ompc(i,j,k,3)
          ww = w_aux_ompc(i,j,k,4)
          vx = 0._rkind
          wx = 0._rkind
          uy = 0._rkind
          wy = 0._rkind
          uz = 0._rkind
          vz = 0._rkind
          do l=1,lmax
            ccl = coeff_deriv1_ompc(l,lmax)
            vx = vx+ccl*(w_aux_ompc(i+l,j,k,3)-w_aux_ompc(i-l,j,k,3))
            wx = wx+ccl*(w_aux_ompc(i+l,j,k,4)-w_aux_ompc(i-l,j,k,4))
            uy = uy+ccl*(w_aux_ompc(i,j+l,k,2)-w_aux_ompc(i,j-l,k,2))
            wy = wy+ccl*(w_aux_ompc(i,j+l,k,4)-w_aux_ompc(i,j-l,k,4))
            uz = uz+ccl*(w_aux_ompc(i,j,k+l,2)-w_aux_ompc(i,j,k-l,2))
            vz = vz+ccl*(w_aux_ompc(i,j,k+l,3)-w_aux_ompc(i,j,k-l,3))
          enddo
          vx = vx*dcsidx_ompc(i)
          wx = wx*dcsidx_ompc(i)
          uy = uy*detady_ompc(j)
          wy = wy*detady_ompc(j)
          uz = uz*dzitdz_ompc(k)
          vz = vz*dzitdz_ompc(k)
          omegax = wy-vz
          omegay = uz-wx
          omegaz = vx-uy
          omod2 = omegax*omegax+omegay*omegay+omegaz*omegaz
          psi_ompc(i,j,k,mpsi) = sqrt(omod2)
        enddo
      enddo
    enddo


  endsubroutine insitu_omega_subroutine



  subroutine insitu_ducros_subroutine(nx,ny,nz,ng,visc_order,npsi,mpsi,u0,l0,w_aux_ompc,&
  &coeff_deriv1_ompc,dcsidx_ompc,detady_ompc,dzitdz_ompc,psi_ompc)
    integer, intent(in) :: nx, ny, nz, ng, visc_order, npsi, mpsi
    real(rkind), intent(in) :: u0, l0
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(inout) :: w_aux_ompc
    real(rkind), dimension(1:,1:), intent(in) :: coeff_deriv1_ompc
    real(rkind), dimension(1:), intent(in) :: dcsidx_ompc, detady_ompc, dzitdz_ompc
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(inout) :: psi_ompc
    integer :: i,j,k,l,iercuda
    real(rkind) :: ccl
    real(rkind) :: uu,vv,ww
    real(rkind) :: ux,uy,uz
    real(rkind) :: vx,vy,vz
    real(rkind) :: wx,wy,wz
    real(rkind) :: div, omegax, omegay, omegaz, omod2
    integer :: lmax
    lmax = visc_order/2
    !$omp parallel do default(firstprivate) shared(w_aux_ompc, &
    !$omp& coeff_deriv1_ompc, &
    !$omp& dcsidx_ompc, &
    !$omp& detady_ompc, &
    !$omp& dzitdz_ompc, &
    !$omp& psi_ompc)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          uu = w_aux_ompc(i,j,k,2)
          vv = w_aux_ompc(i,j,k,3)
          ww = w_aux_ompc(i,j,k,4)
          ux = 0._rkind
          vx = 0._rkind
          wx = 0._rkind
          uy = 0._rkind
          vy = 0._rkind
          wy = 0._rkind
          uz = 0._rkind
          vz = 0._rkind
          wz = 0._rkind
          do l=1,lmax
            ccl = coeff_deriv1_ompc(l,lmax)
            ux = ux+ccl*(w_aux_ompc(i+l,j,k,2)-w_aux_ompc(i-l,j,k,2))
            vx = vx+ccl*(w_aux_ompc(i+l,j,k,3)-w_aux_ompc(i-l,j,k,3))
            wx = wx+ccl*(w_aux_ompc(i+l,j,k,4)-w_aux_ompc(i-l,j,k,4))
            uy = uy+ccl*(w_aux_ompc(i,j+l,k,2)-w_aux_ompc(i,j-l,k,2))
            vy = vy+ccl*(w_aux_ompc(i,j+l,k,3)-w_aux_ompc(i,j-l,k,3))
            wy = wy+ccl*(w_aux_ompc(i,j+l,k,4)-w_aux_ompc(i,j-l,k,4))
            uz = uz+ccl*(w_aux_ompc(i,j,k+l,2)-w_aux_ompc(i,j,k-l,2))
            vz = vz+ccl*(w_aux_ompc(i,j,k+l,3)-w_aux_ompc(i,j,k-l,3))
            wz = wz+ccl*(w_aux_ompc(i,j,k+l,4)-w_aux_ompc(i,j,k-l,4))
          enddo
          ux = ux*dcsidx_ompc(i)
          vx = vx*dcsidx_ompc(i)
          wx = wx*dcsidx_ompc(i)
          uy = uy*detady_ompc(j)
          vy = vy*detady_ompc(j)
          wy = wy*detady_ompc(j)
          uz = uz*dzitdz_ompc(k)
          vz = vz*dzitdz_ompc(k)
          wz = wz*dzitdz_ompc(k)
          div = ux+vy+wz
          omegax = wy-vz
          omegay = uz-wx
          omegaz = vx-uy
          omod2 = omegax*omegax+omegay*omegay+omegaz*omegaz
          psi_ompc(i,j,k,mpsi) = max(-div/sqrt(omod2+div**2+(u0/l0)**2),0._rkind)
        enddo
      enddo
    enddo


  endsubroutine insitu_ducros_subroutine





  subroutine probe_interpolation_subroutine(num_probe,nx,ny,nz,ng,nv_aux,ijk_probe_ompc,&
  &w_aux_probe_ompc,w_aux_ompc,probe_coeff_ompc)
    integer, intent(in) :: num_probe,nx,ny,nz,ng,nv_aux
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv_aux), intent(in) :: w_aux_ompc
    real(rkind), dimension(6,num_probe), intent(inout) :: w_aux_probe_ompc
    real(rkind), dimension(2,2,2,num_probe), intent(in) :: probe_coeff_ompc
    integer, dimension(3,num_probe), intent(in) :: ijk_probe_ompc
    integer :: i,j,k,ii,jj,kk,l,iercuda
    real(rkind) :: w1,w2,w3,w4,w5,w6

    !$omp parallel do default(firstprivate) shared(w_aux_ompc, &
    !$omp& w_aux_probe_ompc, &
    !$omp& probe_coeff_ompc, &
    !$omp& ijk_probe_ompc)
    do l = 1,num_probe
      ii = ijk_probe_ompc(1,l)
      jj = ijk_probe_ompc(2,l)
      kk = ijk_probe_ompc(3,l)
      w1 = 0._rkind
      w2 = 0._rkind
      w3 = 0._rkind
      w4 = 0._rkind
      w5 = 0._rkind
      w6 = 0._rkind
      do k=1,2
        do j=1,2
          do i=1,2
            w1 = w1 + probe_coeff_ompc(i,j,k,l)*w_aux_ompc(i+ii-1,j+jj-1,k+kk-1,1)
            w2 = w2 + probe_coeff_ompc(i,j,k,l)*w_aux_ompc(i+ii-1,j+jj-1,k+kk-1,2)
            w3 = w3 + probe_coeff_ompc(i,j,k,l)*w_aux_ompc(i+ii-1,j+jj-1,k+kk-1,3)
            w4 = w4 + probe_coeff_ompc(i,j,k,l)*w_aux_ompc(i+ii-1,j+jj-1,k+kk-1,4)
            w5 = w5 + probe_coeff_ompc(i,j,k,l)*w_aux_ompc(i+ii-1,j+jj-1,k+kk-1,5)
            w6 = w6 + probe_coeff_ompc(i,j,k,l)*w_aux_ompc(i+ii-1,j+jj-1,k+kk-1,6)
          enddo
        enddo
      enddo
      w_aux_probe_ompc(1,l) = w1
      w_aux_probe_ompc(2,l) = w2
      w_aux_probe_ompc(3,l) = w3
      w_aux_probe_ompc(4,l) = w4
      w_aux_probe_ompc(5,l) = w5
      w_aux_probe_ompc(6,l) = w6
    enddo


  endsubroutine probe_interpolation_subroutine



  subroutine compute_tspec_subroutine(nx,ny,nz,ng,ndft,j_slice,i_win,it_win,w_aux_ompc,w_tspec_ompc)
    integer, intent(in) :: nx, ny, nz, ng, ndft, j_slice, i_win, it_win
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in) :: w_aux_ompc
    real(rkind), dimension(1:,1:,1:,1:,1:), intent(inout) :: w_tspec_ompc
    real(rkind) :: wei_real, wei_imag
    real(rkind) :: pp, win_scale
    integer :: i,j,k,iercuda, ll, l, nn, n
    real(rkind) :: pi
    pi = acos(-1._rkind)
    !$omp parallel do default(firstprivate) shared(w_aux_ompc, &
    !$omp& w_tspec_ompc)
    do k = 1,nz
      do i = 1,nx
        do l=1,ndft
          n=it_win
          ll = l-1
          nn = n-1
          wei_real = cos((-2.0_rkind*pi*ll*nn)/ndft)
          wei_imag = sin((-2.0_rkind*pi*ll*nn)/ndft)
          win_scale = 1._rkind
          pp = win_scale*w_aux_ompc(i,j_slice,k,1)*w_aux_ompc(i,j_slice,k,6)
          if(nn==0) then
            w_tspec_ompc(i,k,l,i_win,1) = pp*wei_real
            w_tspec_ompc(i,k,l,i_win,2) = pp*wei_imag
          else
            w_tspec_ompc(i,k,l,i_win,1) = w_tspec_ompc(i,k,l,i_win,1) + pp*wei_real
            w_tspec_ompc(i,k,l,i_win,2) = w_tspec_ompc(i,k,l,i_win,2) + pp*wei_imag
          endif
        enddo
      enddo
    enddo


  endsubroutine compute_tspec_subroutine



  subroutine compute_psd_tspec_subroutine(nx,ny,nz,ndft,i_win,dt_tspec,w_tspec_ompc,w_psd_tspec_ompc)
    integer, intent(in) :: nx, ny, nz, ndft, i_win
    real(rkind), intent(in) :: dt_tspec
    real(rkind), dimension(1:,1:,1:,1:,1:), intent(inout) :: w_tspec_ompc
    real(rkind), dimension(1:,1:), intent(inout) :: w_psd_tspec_ompc
    integer :: i,j,k,iercuda, ll, l, nn, n
    real(rkind) :: pi
    pi = acos(-1._rkind)
    !$omp parallel do default(firstprivate) shared(w_tspec_ompc, &
    !$omp& w_psd_tspec_ompc)
    do l = 1,ndft
      do i = 1,nx
        w_psd_tspec_ompc(i,l) = 0._rkind
        do k=1,nz
          w_psd_tspec_ompc(i,l) = w_psd_tspec_ompc(i,l)+(w_tspec_ompc(i,k,l,i_win,1)**2+w_tspec_ompc(i,k,l,i_win,2)**2)
        enddo
      enddo
    enddo


  endsubroutine compute_psd_tspec_subroutine



  subroutine compute_wallprop_c2_subroutine(nx,ny,nz,ng,w_aux_ompc,wallprop_ompc,dxdcsic2_ompc,&
  &dydcsic2_ompc,csimod_ompc,meta_ompc,wall_tag_ompc)
    integer, intent(in) :: nx,ny,nz,ng
    integer, dimension(1-ng:nx+ng), intent(in) :: wall_tag_ompc
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(in) :: w_aux_ompc
    real(rkind), dimension(1-ng:,1-ng:, 2:), intent(inout) :: wallprop_ompc
    real(rkind), dimension(1-ng:,1-ng:), intent(in) :: meta_ompc, csimod_ompc
    real(rkind), dimension(1-ng:,1-ng:), intent(in) :: dxdcsic2_ompc, dydcsic2_ompc
    integer :: i,j,k,iercuda
    real(rkind) :: ut1,ut2,ut3,ut4,ds,costh,sinth,dudy,dudyw

    !$omp parallel do default(firstprivate) shared(w_aux_ompc, &
    !$omp& wallprop_ompc, &
    !$omp& meta_ompc, &
    !$omp& csimod_ompc, &
    !$omp& dxdcsic2_ompc, &
    !$omp& dydcsic2_ompc, &
    !$omp& wall_tag_ompc)
    do k = 1,nz
      do i = 1,nx
        if(wall_tag_ompc(i) < 1) then
          ds = csimod_ompc(i,1)
          costh = dxdcsic2_ompc(i,1)/ds
          sinth = dydcsic2_ompc(i,1)/ds
          ut1 = w_aux_ompc(i,1,k,2)*costh+w_aux_ompc(i,1,k,3)*sinth
          ut2 = w_aux_ompc(i,2,k,2)*costh+w_aux_ompc(i,2,k,3)*sinth
          ut3 = w_aux_ompc(i,3,k,2)*costh+w_aux_ompc(i,3,k,3)*sinth
          ut4 = w_aux_ompc(i,4,k,2)*costh+w_aux_ompc(i,4,k,3)*sinth
          dudy = -22._rkind*ut1+36._rkind*ut2-18._rkind*ut3+4._rkind*ut4
          dudyw = dudy*meta_ompc(i,1)/12._rkind
          wallprop_ompc(i,k,2) = w_aux_ompc(i,1,k,7)*dudyw
        else
          wallprop_ompc(i,k,2) = w_aux_ompc(i,1,k,2)
        endif
      enddo
    enddo


  endsubroutine compute_wallprop_c2_subroutine



endmodule streams_kernels_ompc


