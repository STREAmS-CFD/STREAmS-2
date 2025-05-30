module streams_kernels_omp

  use streams_parameters, only : rkind, ikind, real64
  use utils_omp
  use omp_lib
  implicit none

contains

  subroutine zero_flux_kernel(nx,ny,nz,nv,fl_gpu)
    integer :: nx, ny, nz, nv
    real(rkind), dimension(1:,1:,1:,1:), intent(inout) :: fl_gpu
    integer :: i,j,k,m,iercuda

    !$omp target teams distribute parallel do collapse(3) has_device_addr(fl_gpu)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          do m=1,nv
            fl_gpu(i,j,k,m) = 0._rkind
          enddo
        enddo
      enddo
    enddo


  endsubroutine zero_flux_kernel



  subroutine init_flux_kernel(nx,ny,nz,nv,fl_gpu,fln_gpu,rhodt)
    integer :: nx, ny, nz, nv
    real(rkind) :: rhodt
    real(rkind), dimension(1:,1:,1:,1:), intent(inout) :: fl_gpu, fln_gpu
    integer :: i,j,k,m,iercuda

    !$omp target teams distribute parallel do collapse(3) has_device_addr(fl_gpu,fln_gpu)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          do m=1,nv
            fln_gpu(i,j,k,m) = - rhodt * fl_gpu(i,j,k,m)
            fl_gpu(i,j,k,m) = 0._rkind
          enddo
        enddo
      enddo
    enddo


  endsubroutine init_flux_kernel




  subroutine count_weno_kernel(nv,nv_aux,nx,ny,nz,ng,eul_imin,eul_imax,eul_jmin,eul_jmax,eul_kmin,&
  &eul_kmax,weno_scheme,sensor_threshold,w_aux_gpu,ep_ord_change_gpu,count_weno_x,count_weno_y,&
  &count_weno_z)
    integer, intent(in) :: nv, nv_aux, nx, ny, nz, ng
    integer, intent(in) :: eul_imin, eul_imax, eul_jmin, eul_jmax, eul_kmin, eul_kmax
    integer, intent(in) :: weno_scheme
    real(rkind) :: sensor_threshold
    real(rkind), intent(out) :: count_weno_x, count_weno_y, count_weno_z
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv_aux) :: w_aux_gpu
    integer, dimension(0:nx,0:ny,0:nz,1:3) :: ep_ord_change_gpu
    integer :: i,j,k,ishk,ii,jj,kk,iercuda
    count_weno_x = 0._rkind
    !$omp target teams distribute parallel do collapse(3) has_device_addr(w_aux_gpu,ep_ord_change_gpu) &
    !$omp& reduction(+:count_weno_x)
    do k = 1,nz
      do j = 1,ny
        do i = eul_imin,eul_imax
          ishk = 0
          do ii=i-weno_scheme+1,i+weno_scheme
            if (w_aux_gpu(ii,j,k,8) > sensor_threshold) then
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
    !$omp target teams distribute parallel do collapse(3) has_device_addr(w_aux_gpu,ep_ord_change_gpu) &
    !$omp& reduction(+:count_weno_y)
    do k = 1,nz
      do j = eul_jmin,eul_jmax
        do i = 1,nx
          ishk = 0
          do jj=j-weno_scheme+1,j+weno_scheme
            if (w_aux_gpu(i,jj,k,8) > sensor_threshold) then
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
    !$omp target teams distribute parallel do collapse(3) has_device_addr(w_aux_gpu,ep_ord_change_gpu) &
    !$omp& reduction(+:count_weno_z)
    do k = eul_kmin,eul_kmax
      do j = 1,ny
        do i = 1,nx
          ishk = 0
          do kk=k-weno_scheme+1,k+weno_scheme
            if (w_aux_gpu(i,j,kk,8) > sensor_threshold) then
              ishk = 1
            endif
          enddo
          if (ishk > 0) then
            count_weno_z = count_weno_z + 1
          endif
        enddo
      enddo
    enddo


  endsubroutine count_weno_kernel



  subroutine count_weno_c2_kernel(nv,nv_aux,nx,ny,nz,ng,eul_imin,eul_imax,eul_jmin,eul_jmax,&
  &eul_kmin,eul_kmax,l_base,weno_scheme,sensor_threshold,w_aux_gpu,ep_ord_change_gpu,lmax_tag_gpu,&
  &wall_tag_gpu,count_weno_x,count_weno_y,count_weno_z)
    integer, intent(in) :: nv, nv_aux, nx, ny, nz, ng
    integer, intent(in) :: eul_imin, eul_imax, eul_jmin, eul_jmax, eul_kmin, eul_kmax
    integer, intent(in) :: l_base, weno_scheme
    real(rkind) :: sensor_threshold
    real(rkind), intent(out) :: count_weno_x, count_weno_y, count_weno_z
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv_aux) :: w_aux_gpu
    integer, dimension(0:nx,0:ny,0:nz,1:3) :: ep_ord_change_gpu
    integer, dimension(1-ng:nx+ng) :: lmax_tag_gpu
    integer, dimension(1-ng:nx+ng), intent(in) :: wall_tag_gpu
    integer :: i, j, k, ishk, weno_scheme_i, weno_scheme_j, weno_scheme_k, ii, jj, kk, iercuda
    count_weno_x = 0._rkind
    !$omp target teams distribute parallel do collapse(3) has_device_addr(w_aux_gpu,ep_ord_change_gpu,lmax_tag_gpu,wall_tag_gpu) &
    !$omp& reduction(+:count_weno_x)
    do k = 1,nz
      do j = 1,ny
        do i = eul_imin,eul_imax
          weno_scheme_i = weno_scheme+ep_ord_change_gpu(i,j,k,1)
          if(weno_scheme_i < 1) then
            weno_scheme_i = 1
          endif
          if(j == 1) then
            weno_scheme_i = lmax_tag_gpu(i)+ep_ord_change_gpu(i,j,k,1)
            if(weno_scheme_i < 1) then
              weno_scheme_i = 1
            endif
          endif
          ishk = 0
          do ii=i-weno_scheme_i+1,i+weno_scheme_i
            if (w_aux_gpu(ii,j,k,8) > sensor_threshold) then
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
    !$omp target teams distribute parallel do collapse(3) has_device_addr(w_aux_gpu,ep_ord_change_gpu,lmax_tag_gpu,wall_tag_gpu) &
    !$omp& reduction(+:count_weno_y)
    do k = 1,nz
      do j = eul_jmin,eul_jmax
        do i = 1,nx
          weno_scheme_j = weno_scheme+ep_ord_change_gpu(i,j,k,2)
          if(weno_scheme_j < 1) then
            weno_scheme_j = 1
          endif
          if (j <= l_base .and. wall_tag_gpu(i) > 0) then
            weno_scheme_j = weno_scheme
          endif
          ishk = 0
          do jj=j-weno_scheme_j+1,j+weno_scheme_j
            if (w_aux_gpu(i,jj,k,8) > sensor_threshold) then
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
    !$omp target teams distribute parallel do collapse(3) has_device_addr(w_aux_gpu,ep_ord_change_gpu,lmax_tag_gpu,wall_tag_gpu) &
    !$omp& reduction(+:count_weno_z)
    do k = eul_kmin,eul_kmax
      do j = 1,ny
        do i = 1,nx
          ishk = 0
          do kk=k-weno_scheme+1,k+weno_scheme
            if (w_aux_gpu(i,j,kk,8) > sensor_threshold) then
              ishk = 1
            endif
          enddo
          if (ishk > 0) then
            count_weno_z = count_weno_z + 1
          endif
        enddo
      enddo
    enddo


  endsubroutine count_weno_c2_kernel



  subroutine euler_x_fluxes_hybrid_c2_kernel(nv,nv_aux,nx,ny,nz,ng,eul_imin,eul_imax,lmax_base,&
  &nkeep,rgas0,w_aux_gpu,coeff_deriv1_gpu,fhat_gpu,dcsidxc2_gpu,dcsidxnc2_gpu,dcsidyc2_gpu,&
  &dcsidync2_gpu,jac_gpu,mcsijac1_gpu,lmax_tag_gpu,force_zero_flux_min,force_zero_flux_max,weno_scheme,&
  &weno_version,sensor_threshold,weno_size,cp_coeff_gpu,indx_cp_l,indx_cp_r,ep_ord_change_gpu,&
  &calorically_perfect,tol_iter_nr,rho0,u0,t0)
    integer, value :: nv, nx, ny, nz, ng, nv_aux
    integer, value :: eul_imin, eul_imax, lmax_base, nkeep, indx_cp_l, indx_cp_r, calorically_perfect
    integer, value :: weno_scheme, weno_size, weno_version
    integer, dimension(0:nx,0:ny,0:nz,1:3) :: ep_ord_change_gpu
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv_aux) :: w_aux_gpu
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv) :: fhat_gpu
    real(rkind), dimension(indx_cp_l:indx_cp_r+1) :: cp_coeff_gpu
    real(rkind), dimension(4,4) :: coeff_deriv1_gpu
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng) :: dcsidxc2_gpu
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng) :: dcsidxnc2_gpu
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng) :: dcsidyc2_gpu
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng) :: dcsidync2_gpu
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng) :: jac_gpu
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng) :: mcsijac1_gpu
    integer, dimension(1-ng:nx+ng) :: lmax_tag_gpu
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
    !$omp target data map(alloc:el,er,evmax,fi,gp,gm,eler)
    !$omp target teams distribute parallel do collapse(3) has_device_addr(w_aux_gpu,fhat_gpu,cp_coeff_gpu,coeff_deriv1_gpu,dcsidxc2_gpu,dcsidxnc2_gpu,dcsidyc2_gpu,dcsidync2_gpu,jac_gpu,mcsijac1_gpu,ep_ord_change_gpu,lmax_tag_gpu) private(el,er,evmax,fi,gp,gm,eler)
    do k = 1,nz
      do j = 1,ny
        do i = +eul_imin-2+1,eul_imax
          weno_scheme_i = max(weno_scheme+ep_ord_change_gpu(i,j,k,1),1)
          lmax = max(lmax_base+ep_ord_change_gpu(i,j,k,1),1)
          if(j == 1) then
            weno_scheme_i = max(lmax_tag_gpu(i)+ep_ord_change_gpu(i,j,k,1),1)
            lmax = weno_scheme_i
          endif
          weno_size = 2*weno_scheme_i

          ishk = 0
          do ii=i-weno_scheme_i+1,i+weno_scheme_i
            if (w_aux_gpu(ii,j,k,8) > sensor_threshold) ishk = 1
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

                rhoi = w_aux_gpu(i-m,j,k,1)
                uui = w_aux_gpu(i-m,j,k,2)
                vvi = w_aux_gpu(i-m,j,k,3)
                wwi = w_aux_gpu(i-m,j,k,4)
                enti = w_aux_gpu(i-m,j,k,5)
                tti = w_aux_gpu(i-m,j,k,6)
                ppi = tti*rhoi*rgas0
                uti = (uui*dcsidxc2_gpu(i-m,j) + vvi*dcsidyc2_gpu(i-m,j))/jac_gpu(i-m,j)
                dcsidxi = dcsidxc2_gpu(i-m,j)/jac_gpu(i-m,j)
                dcsidyi = dcsidyc2_gpu(i-m,j)/jac_gpu(i-m,j)

                rhoip = w_aux_gpu(i-m+l,j,k,1)
                uuip = w_aux_gpu(i-m+l,j,k,2)
                vvip = w_aux_gpu(i-m+l,j,k,3)
                wwip = w_aux_gpu(i-m+l,j,k,4)
                entip = w_aux_gpu(i-m+l,j,k,5)
                ttip = w_aux_gpu(i-m+l,j,k,6)
                ppip = ttip*rhoip*rgas0
                utip = (uuip*dcsidxc2_gpu(i-m+l,j) + vvip*dcsidyc2_gpu(i-m+l,j))/jac_gpu(i-m+l,j)
                dcsidxip = dcsidxc2_gpu(i-m+l,j)/jac_gpu(i-m+l,j)
                dcsidyip = dcsidyc2_gpu(i-m+l,j)/jac_gpu(i-m+l,j)

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
              ft1 = ft1 + coeff_deriv1_gpu(l,lmax)*uvs1
              ft2 = ft2 + coeff_deriv1_gpu(l,lmax)*uvs2
              ft3 = ft3 + coeff_deriv1_gpu(l,lmax)*uvs3
              ft4 = ft4 + coeff_deriv1_gpu(l,lmax)*uvs4
              ft5 = ft5 + coeff_deriv1_gpu(l,lmax)*uvs5
              ft6 = ft6 + coeff_deriv1_gpu(l,lmax)*uvs6
              ft7 = ft7 + coeff_deriv1_gpu(l,lmax)*uvs7
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

            fhat_gpu(i,j,k,1) = fh1
            fhat_gpu(i,j,k,2) = fh2
            fhat_gpu(i,j,k,3) = fh3
            fhat_gpu(i,j,k,4) = fh4
            fhat_gpu(i,j,k,5) = fh5
          else
            call compute_roe_average(nx, ny, nz, ng, i, i+1, j, j, k, k, w_aux_gpu, rgas0, b1, b2,&
            & b3, c, ci, h, uu, vv, ww, cp_coeff_gpu, indx_cp_l, indx_cp_r, calorically_perfect, tol_iter_nr,t0)
            dcsidxav = 0.5_rkind * (dcsidxnc2_gpu(i,j) + dcsidxnc2_gpu(i+1,j))
            dcsidyav = 0.5_rkind * (dcsidync2_gpu(i,j) + dcsidync2_gpu(i+1,j))
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
              uu = w_aux_gpu(ll,j,k,2)
              vv = w_aux_gpu(ll,j,k,3)
              ut = dcsidxnc2_gpu(ll,j)*uu+dcsidync2_gpu(ll,j)*vv
              tt = w_aux_gpu(ll,j,k,6)
              gamloc = get_gamloc_dev(indx_cp_l,indx_cp_r,tt,t0,cp_coeff_gpu,calorically_perfect,rgas0)
              c = sqrt (gamloc*rgas0*tt)

              evmax(1) = max(abs(ut-c),evmax(1))
              evmax(2) = max(abs(ut ),evmax(2))
              evmax(3) = evmax(2)
              evmax(4) = evmax(2)
              evmax(5) = max(abs(ut+c),evmax(5))
            enddo
            do l=1,weno_size
              ll = i + l - weno_scheme_i

              rho = w_aux_gpu(ll,j,k,1)
              uu = w_aux_gpu(ll,j,k,2)
              vv = w_aux_gpu(ll,j,k,3)
              ut = dcsidxnc2_gpu(ll,j)*uu+dcsidync2_gpu(ll,j)*vv
              ww = w_aux_gpu(ll,j,k,4)
              h = w_aux_gpu(ll,j,k,5)
              pp = rho*w_aux_gpu(ll,j,k,6)*rgas0
              fi(1) = rho * ut
              fi(2) = (rho * uu * ut + pp * dcsidxnc2_gpu(ll,j))
              fi(3) = (rho * vv * ut + pp * dcsidync2_gpu(ll,j))
              fi(4) = rho * ww * ut
              fi(5) = rho * h * ut
              do m=1,5
                wc = 0._rkind
                gc = 0._rkind

                wc = wc + el(1,m) * rho * mcsijac1_gpu(ll,j)
                gc = gc + el(1,m) * fi(1) * mcsijac1_gpu(ll,j)
                wc = wc + el(2,m) * rho*uu * mcsijac1_gpu(ll,j)
                gc = gc + el(2,m) * fi(2) * mcsijac1_gpu(ll,j)
                wc = wc + el(3,m) * rho*vv * mcsijac1_gpu(ll,j)
                gc = gc + el(3,m) * fi(3) * mcsijac1_gpu(ll,j)
                wc = wc + el(4,m) * rho*ww * mcsijac1_gpu(ll,j)
                gc = gc + el(4,m) * fi(4) * mcsijac1_gpu(ll,j)
                wc = wc + el(5,m) * (rho*h-pp) * mcsijac1_gpu(ll,j)
                gc = gc + el(5,m) * fi(5) * mcsijac1_gpu(ll,j)

                c = 0.5_rkind * (gc + evmax(m) * wc)
                gp(m,l) = c
                gm(m,l) = gc - c
              enddo
            enddo
            call wenorec_1d(nv,gp,gm,fi,weno_scheme_i,weno_scheme_i,weno_version,rho0,u0)
            do m=1,5
              fhat_gpu(i,j,k,m) = 0._rkind
              do mm=1,5
                fhat_gpu(i,j,k,m) = fhat_gpu(i,j,k,m) + er(mm,m) * fi(mm)
              enddo
            enddo

          endif
        enddo
      enddo
    enddo

    !$omp end target data
  endsubroutine euler_x_fluxes_hybrid_c2_kernel


  subroutine euler_x_fluxes_hybrid_kernel(nv,nv_aux,nx,ny,nz,ng,istart_face,iend_face,lmax_base,&
  &nkeep,rgas0,w_aux_gpu,coeff_deriv1_gpu,dcsidx_gpu,fhat_gpu,force_zero_flux_min,force_zero_flux_max,&
  &weno_scheme,weno_version,sensor_threshold,weno_size,cp_coeff_gpu,indx_cp_l,indx_cp_r,&
  &ep_ord_change_gpu,calorically_perfect,tol_iter_nr,rho0,u0,t0)
    integer, value :: nv, nx, ny, nz, ng, nv_aux
    integer, value :: istart_face, iend_face, lmax_base, nkeep, indx_cp_l, indx_cp_r, calorically_perfect
    integer, value :: weno_scheme, weno_size, weno_version
    integer, dimension(0:nx,0:ny,0:nz,1:3) :: ep_ord_change_gpu
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv_aux) :: w_aux_gpu
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv) :: fhat_gpu
    real(rkind), dimension(indx_cp_l:indx_cp_r+1) :: cp_coeff_gpu
    real(rkind), dimension(4,4) :: coeff_deriv1_gpu
    real(rkind), dimension(nx) :: dcsidx_gpu
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
    !$omp target data map(alloc:el,er,evmax,fi,gp,gm)
    !$omp target teams distribute parallel do collapse(2) has_device_addr(w_aux_gpu,fhat_gpu,cp_coeff_gpu,coeff_deriv1_gpu,dcsidx_gpu,ep_ord_change_gpu) private(el,er,evmax,fi,gp,gm)
    do j = 1,ny
      do i = +istart_face-1+1,iend_face
        do k = 1,nz

          ishk = 0
          do ii=i-weno_scheme+1,i+weno_scheme
            if (w_aux_gpu(ii,j,k,8) > sensor_threshold) ishk = 1
          enddo

          if (ishk == 0) then

            ft1 = 0._rkind
            ft2 = 0._rkind
            ft3 = 0._rkind
            ft4 = 0._rkind
            ft5 = 0._rkind
            ft6 = 0._rkind
            lmax = max(lmax_base+ep_ord_change_gpu(i,j,k,1),1)
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

                  rhoi = w_aux_gpu(i-m,j,k,1)
                  uui = w_aux_gpu(i-m,j,k,2)
                  vvi = w_aux_gpu(i-m,j,k,3)
                  wwi = w_aux_gpu(i-m,j,k,4)
                  enti = w_aux_gpu(i-m,j,k,5)
                  tti = w_aux_gpu(i-m,j,k,6)
                  ppi = tti*rhoi*rgas0
                  eei = enti-ppi/rhoi-0.5_rkind*(uui*uui+vvi*vvi+wwi*wwi)

                  rhoip = w_aux_gpu(i-m+l,j,k,1)
                  uuip = w_aux_gpu(i-m+l,j,k,2)
                  vvip = w_aux_gpu(i-m+l,j,k,3)
                  wwip = w_aux_gpu(i-m+l,j,k,4)
                  entip = w_aux_gpu(i-m+l,j,k,5)
                  ttip = w_aux_gpu(i-m+l,j,k,6)
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
                ft1 = ft1 + coeff_deriv1_gpu(l,lmax)*uvs1
                ft2 = ft2 + coeff_deriv1_gpu(l,lmax)*uvs2
                ft3 = ft3 + coeff_deriv1_gpu(l,lmax)*uvs3
                ft4 = ft4 + coeff_deriv1_gpu(l,lmax)*uvs4
                ft5 = ft5 + coeff_deriv1_gpu(l,lmax)*(uvs5_i+uvs5_k+uvs5_p)
                ft6 = ft6 + coeff_deriv1_gpu(l,lmax)*uvs6
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

                  rhoi = w_aux_gpu(i-m,j,k,1)
                  uui = w_aux_gpu(i-m,j,k,2)
                  vvi = w_aux_gpu(i-m,j,k,3)
                  wwi = w_aux_gpu(i-m,j,k,4)
                  enti = w_aux_gpu(i-m,j,k,5)
                  tti = w_aux_gpu(i-m,j,k,6)
                  ppi = tti*rhoi*rgas0

                  rhoip = w_aux_gpu(i-m+l,j,k,1)
                  uuip = w_aux_gpu(i-m+l,j,k,2)
                  vvip = w_aux_gpu(i-m+l,j,k,3)
                  wwip = w_aux_gpu(i-m+l,j,k,4)
                  entip = w_aux_gpu(i-m+l,j,k,5)
                  ttip = w_aux_gpu(i-m+l,j,k,6)
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
                ft1 = ft1 + coeff_deriv1_gpu(l,lmax)*uvs1
                ft2 = ft2 + coeff_deriv1_gpu(l,lmax)*uvs2
                ft3 = ft3 + coeff_deriv1_gpu(l,lmax)*uvs3
                ft4 = ft4 + coeff_deriv1_gpu(l,lmax)*uvs4
                ft5 = ft5 + coeff_deriv1_gpu(l,lmax)*uvs5
                ft6 = ft6 + coeff_deriv1_gpu(l,lmax)*uvs6
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

            fhat_gpu(i,j,k,1) = fh1
            fhat_gpu(i,j,k,2) = fh2
            fhat_gpu(i,j,k,3) = fh3
            fhat_gpu(i,j,k,4) = fh4
            fhat_gpu(i,j,k,5) = fh5
          else
            call compute_roe_average(nx, ny, nz, ng, i, i+1, j, j, k, k, w_aux_gpu, rgas0, b1, b2,&
            & b3, c, ci, h, uu, vv, ww, cp_coeff_gpu, indx_cp_l, indx_cp_r, calorically_perfect, tol_iter_nr,t0)

            call eigenvectors_x(b1, b2, b3, uu, vv, ww, c, ci, h, el, er)

            do m=1,5
              evmax(m) = -1._rkind
            enddo
            do l=1,weno_size
              ll = i + l - weno_scheme
              uu = w_aux_gpu(ll,j,k,2)
              tt = w_aux_gpu(ll,j,k,6)
              gamloc = get_gamloc_dev(indx_cp_l,indx_cp_r,tt,t0,cp_coeff_gpu,calorically_perfect,rgas0)
              c = sqrt (gamloc*rgas0*tt)
              evmax(1) = max(abs(uu-c),evmax(1))
              evmax(2) = max(abs(uu ),evmax(2))
              evmax(3) = evmax(2)
              evmax(4) = evmax(2)
              evmax(5) = max(abs(uu+c),evmax(5))
            enddo
            do l=1,weno_size
              ll = i + l - weno_scheme

              rho = w_aux_gpu(ll,j,k,1)
              uu = w_aux_gpu(ll,j,k,2)
              vv = w_aux_gpu(ll,j,k,3)
              ww = w_aux_gpu(ll,j,k,4)
              h = w_aux_gpu(ll,j,k,5)
              rhou = rho*uu
              pp = rho*w_aux_gpu(ll,j,k,6)*rgas0
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
            wenorec_ord = max(weno_scheme+ep_ord_change_gpu(i,j,k,1),1)
            call wenorec_1d(nv,gp,gm,fi,weno_scheme,wenorec_ord,weno_version,rho0,u0)
            do m=1,5
              fhat_gpu(i,j,k,m) = 0._rkind
              do mm=1,5
                fhat_gpu(i,j,k,m) = fhat_gpu(i,j,k,m) + er(mm,m) * fi(mm)
              enddo
            enddo

          endif
        enddo
      enddo
    enddo

    !$omp end target data
  endsubroutine euler_x_fluxes_hybrid_kernel


  subroutine euler_x_fluxes_hybrid_rusanov_kernel(nv,nv_aux,nx,ny,nz,ng,istart_face,iend_face,&
  &lmax_base,nkeep,rgas0,w_aux_gpu,coeff_deriv1_gpu,dcsidx_gpu,fhat_gpu,force_zero_flux_min,&
  &force_zero_flux_max,weno_scheme,weno_version,sensor_threshold,weno_size,cp_coeff_gpu,indx_cp_l,&
  &indx_cp_r,ep_ord_change_gpu,calorically_perfect,tol_iter_nr,rho0,u0,t0)
    integer, value :: nv, nx, ny, nz, ng, nv_aux
    integer, value :: istart_face, iend_face, lmax_base, nkeep, indx_cp_l, indx_cp_r, calorically_perfect
    integer, value :: weno_scheme, weno_size, weno_version
    integer, dimension(0:nx,0:ny,0:nz,1:3) :: ep_ord_change_gpu
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv_aux) :: w_aux_gpu
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv) :: fhat_gpu
    real(rkind), dimension(indx_cp_l:indx_cp_r+1) :: cp_coeff_gpu
    real(rkind), dimension(4,4) :: coeff_deriv1_gpu
    real(rkind), dimension(nx) :: dcsidx_gpu
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
    !$omp target data map(alloc:fi,gp,gm)
    !$omp target teams distribute parallel do collapse(3) has_device_addr(w_aux_gpu,fhat_gpu,cp_coeff_gpu,coeff_deriv1_gpu,dcsidx_gpu,ep_ord_change_gpu) private(fi,gp,gm)
    do k = 1,nz
      do j = 1,ny
        do i = +istart_face-1+1,iend_face
          ishk = 0
          do ii=i-weno_scheme+1,i+weno_scheme
            if (w_aux_gpu(ii,j,k,8) > sensor_threshold) ishk = 1
          enddo

          if (ishk == 0) then

            ft1 = 0._rkind
            ft2 = 0._rkind
            ft3 = 0._rkind
            ft4 = 0._rkind
            ft5 = 0._rkind
            ft6 = 0._rkind
            lmax = max(lmax_base+ep_ord_change_gpu(i,j,k,1),1)
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

                  rhoi = w_aux_gpu(i-m,j,k,1)
                  uui = w_aux_gpu(i-m,j,k,2)
                  vvi = w_aux_gpu(i-m,j,k,3)
                  wwi = w_aux_gpu(i-m,j,k,4)
                  enti = w_aux_gpu(i-m,j,k,5)
                  tti = w_aux_gpu(i-m,j,k,6)
                  ppi = tti*rhoi*rgas0
                  eei = enti-ppi/rhoi-0.5_rkind*(uui*uui+vvi*vvi+wwi*wwi)

                  rhoip = w_aux_gpu(i-m+l,j,k,1)
                  uuip = w_aux_gpu(i-m+l,j,k,2)
                  vvip = w_aux_gpu(i-m+l,j,k,3)
                  wwip = w_aux_gpu(i-m+l,j,k,4)
                  entip = w_aux_gpu(i-m+l,j,k,5)
                  ttip = w_aux_gpu(i-m+l,j,k,6)
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
                ft1 = ft1 + coeff_deriv1_gpu(l,lmax)*uvs1
                ft2 = ft2 + coeff_deriv1_gpu(l,lmax)*uvs2
                ft3 = ft3 + coeff_deriv1_gpu(l,lmax)*uvs3
                ft4 = ft4 + coeff_deriv1_gpu(l,lmax)*uvs4
                ft5 = ft5 + coeff_deriv1_gpu(l,lmax)*(uvs5_i+uvs5_k+uvs5_p)
                ft6 = ft6 + coeff_deriv1_gpu(l,lmax)*uvs6
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

                  rhoi = w_aux_gpu(i-m,j,k,1)
                  uui = w_aux_gpu(i-m,j,k,2)
                  vvi = w_aux_gpu(i-m,j,k,3)
                  wwi = w_aux_gpu(i-m,j,k,4)
                  enti = w_aux_gpu(i-m,j,k,5)
                  tti = w_aux_gpu(i-m,j,k,6)
                  ppi = tti*rhoi*rgas0

                  rhoip = w_aux_gpu(i-m+l,j,k,1)
                  uuip = w_aux_gpu(i-m+l,j,k,2)
                  vvip = w_aux_gpu(i-m+l,j,k,3)
                  wwip = w_aux_gpu(i-m+l,j,k,4)
                  entip = w_aux_gpu(i-m+l,j,k,5)
                  ttip = w_aux_gpu(i-m+l,j,k,6)
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
                ft1 = ft1 + coeff_deriv1_gpu(l,lmax)*uvs1
                ft2 = ft2 + coeff_deriv1_gpu(l,lmax)*uvs2
                ft3 = ft3 + coeff_deriv1_gpu(l,lmax)*uvs3
                ft4 = ft4 + coeff_deriv1_gpu(l,lmax)*uvs4
                ft5 = ft5 + coeff_deriv1_gpu(l,lmax)*uvs5
                ft6 = ft6 + coeff_deriv1_gpu(l,lmax)*uvs6
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

            fhat_gpu(i,j,k,1) = fh1
            fhat_gpu(i,j,k,2) = fh2
            fhat_gpu(i,j,k,3) = fh3
            fhat_gpu(i,j,k,4) = fh4
            fhat_gpu(i,j,k,5) = fh5
          else
            evmax = -1._rkind
            do l=1,weno_size
              ll = i + l - weno_scheme
              uu = w_aux_gpu(ll,j,k,2)
              tt = w_aux_gpu(ll,j,k,6)
              gamloc = get_gamloc_dev(indx_cp_l,indx_cp_r,tt,t0,cp_coeff_gpu,calorically_perfect,rgas0)
              c = sqrt (gamloc*rgas0*tt)
              evm = max(abs(uu-c),abs(uu+c))
              evmax = max(evm,evmax)
            enddo
            do l=1,weno_size
              ll = i + l - weno_scheme

              rho = w_aux_gpu(ll,j,k,1)
              uu = w_aux_gpu(ll,j,k,2)
              vv = w_aux_gpu(ll,j,k,3)
              ww = w_aux_gpu(ll,j,k,4)
              h = w_aux_gpu(ll,j,k,5)
              rhou = rho*uu
              pp = rho*w_aux_gpu(ll,j,k,6)*rgas0

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
            wenorec_ord = max(weno_scheme+ep_ord_change_gpu(i,j,k,1),1)
            call wenorec_1d_rusanov(nv,gp,gm,fi,weno_scheme,wenorec_ord,weno_version,rho0,u0)
            do m=1,5
              fhat_gpu(i,j,k,m) = fi(m)
            enddo

          endif
        enddo
      enddo
    enddo

    !$omp end target data
  endsubroutine euler_x_fluxes_hybrid_rusanov_kernel


  subroutine euler_x_update_kernel(nx,ny,nz,ng,nv,eul_imin,eul_imax,fhat_gpu,fl_gpu,dcsidx_gpu,stream_id)
    integer :: nx,ny,nz,ng,nv,eul_imin,eul_imax
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv), intent(in) :: fhat_gpu
    real(rkind), dimension(1:nx,1:ny,1:nz,1:nv), intent(inout) :: fl_gpu
    real(rkind), dimension(1:nx), intent(in) :: dcsidx_gpu
    integer :: stream_id
    integer :: i,j,k,m,iv,iercuda

    !$omp target teams distribute parallel do collapse(3) has_device_addr(fhat_gpu,fl_gpu,dcsidx_gpu)
    do k = 1,nz
      do j = 1,ny
        do i = eul_imin,eul_imax
          do iv=1,nv
            fl_gpu(i,j,k,iv) = fl_gpu(i,j,k,iv) + (fhat_gpu(i,j,k,iv)-fhat_gpu(i-1,j,k,iv))*dcsidx_gpu(i)
          enddo
        enddo
      enddo
    enddo


  endsubroutine euler_x_update_kernel



  subroutine euler_x_update_c2_kernel(nx,ny,nz,ng,nv,eul_imin,eul_imax,fhat_gpu,fl_gpu,jac_gpu,stream_id)
    integer :: nx,ny,nz,ng,nv,eul_imin,eul_imax
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv), intent(in) :: fhat_gpu
    real(rkind), dimension(1:nx,1:ny,1:nz,1:nv), intent(inout) :: fl_gpu
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng), intent(in) :: jac_gpu
    integer :: stream_id
    integer :: i,j,k,m,iv,iercuda

    !$omp target teams distribute parallel do collapse(3) has_device_addr(fhat_gpu,fl_gpu,jac_gpu)
    do k = 1,nz
      do j = 1,ny
        do i = eul_imin,eul_imax
          do iv=1,nv
            fl_gpu(i,j,k,iv) = fl_gpu(i,j,k,iv) + (fhat_gpu(i,j,k,iv)-fhat_gpu(i-1,j,k,iv))*jac_gpu(i,j)
          enddo
        enddo
      enddo
    enddo


  endsubroutine euler_x_update_c2_kernel




  subroutine euler_y_update_kernel(nx,ny,nz,ng,nv,eul_jmin,eul_jmax,fhat_gpu,fl_gpu,detady_gpu,stream_id)
    integer :: nx,ny,nz,ng,nv,eul_jmin,eul_jmax
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv), intent(in) :: fhat_gpu
    real(rkind), dimension(1:nx,1:ny,1:nz,1:nv), intent(inout) :: fl_gpu
    real(rkind), dimension(1:ny), intent(in) :: detady_gpu
    integer :: stream_id
    integer :: i,j,k,m,iv,iercuda

    !$omp target teams distribute parallel do collapse(3) has_device_addr(fhat_gpu,fl_gpu,detady_gpu)
    do k = 1,nz
      do j = eul_jmin,eul_jmax
        do i = 1,nx
          do iv=1,nv
            fl_gpu(i,j,k,iv) = fl_gpu(i,j,k,iv) + (fhat_gpu(i,j,k,iv)-fhat_gpu(i,j-1,k,iv))*detady_gpu(j)
          enddo
        enddo
      enddo
    enddo


  endsubroutine euler_y_update_kernel



  subroutine euler_y_update_c2_kernel(nx,ny,nz,ng,nv,eul_jmin,eul_jmax,fhat_gpu,fl_gpu,jac_gpu,wall_tag_gpu,stream_id)
    integer :: nx,ny,nz,ng,nv,eul_jmin,eul_jmax
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv), intent(in) :: fhat_gpu
    real(rkind), dimension(1:nx,1:ny,1:nz,1:nv), intent(inout) :: fl_gpu
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng), intent(in) :: jac_gpu
    integer, dimension(1-ng:nx+ng), intent(in) :: wall_tag_gpu
    integer :: stream_id
    integer :: i,j,k,m,iv,iercuda

    !$omp target teams distribute parallel do collapse(3) has_device_addr(fhat_gpu,fl_gpu,jac_gpu,wall_tag_gpu)
    do k = 1,nz
      do j = 1,eul_jmax
        do i = 1,nx
          if(j > 1 .or. eul_jmin == 1 .or. wall_tag_gpu(i) > 0) then
            do iv=1,nv
              fl_gpu(i,j,k,iv) = fl_gpu(i,j,k,iv) + (fhat_gpu(i,j,k,iv)-fhat_gpu(i,j-1,k,iv))*jac_gpu(i,j)
            enddo
          endif
        enddo
      enddo
    enddo


  endsubroutine euler_y_update_c2_kernel


  subroutine euler_z_update_kernel(nx,ny,nz,ng,nv,eul_kmin,eul_kmax,fhat_gpu,fl_gpu,dzitdz_gpu,stream_id)
    integer :: nx,ny,nz,ng,nv,eul_kmin,eul_kmax
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv), intent(in) :: fhat_gpu
    real(rkind), dimension(1:nx,1:ny,1:nz,1:nv), intent(inout) :: fl_gpu
    real(rkind), dimension(1:nz), intent(in) :: dzitdz_gpu
    integer :: stream_id
    integer :: i,j,k,m,iv,iercuda

    !$omp target teams distribute parallel do collapse(3) has_device_addr(fhat_gpu,fl_gpu,dzitdz_gpu)
    do k = eul_kmin,eul_kmax
      do j = 1,ny
        do i = 1,nx
          do iv=1,nv
            fl_gpu(i,j,k,iv) = fl_gpu(i,j,k,iv) + (fhat_gpu(i,j,k,iv)-fhat_gpu(i,j,k-1,iv))*dzitdz_gpu(k)
          enddo
        enddo
      enddo
    enddo


  endsubroutine euler_z_update_kernel




  subroutine euler_z_hybrid_kernel(nv,nv_aux,nx,ny,nz,ng,eul_kmin,eul_kmax,lmax_base,nkeep,rgas0,&
  &w_aux_gpu,fl_gpu,coeff_deriv1_gpu,dzitdz_gpu,fhat_gpu,force_zero_flux_min,force_zero_flux_max,&
  &weno_scheme,weno_version,sensor_threshold,weno_size,cp_coeff_gpu,indx_cp_l,indx_cp_r,&
  &ep_ord_change_gpu,calorically_perfect,tol_iter_nr,rho0,u0,t0)
    integer, value :: nv, nx, ny, nz, ng, nv_aux
    integer, value :: eul_kmin, eul_kmax, lmax_base, nkeep, indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv_aux) :: w_aux_gpu
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv) :: fhat_gpu
    integer, dimension(0:nx,0:ny,0:nz,1:3) :: ep_ord_change_gpu
    real(rkind), dimension(indx_cp_l:indx_cp_r+1) :: cp_coeff_gpu
    real(rkind), dimension(nx,ny,nz,nv) :: fl_gpu
    real(rkind), dimension(4,4) :: coeff_deriv1_gpu
    real(rkind), dimension(nz) :: dzitdz_gpu
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
    !$omp target data map(alloc:el,er,evmax,fk,gp,gm)
    !$omp target teams distribute parallel do collapse(3) has_device_addr(w_aux_gpu,fhat_gpu,cp_coeff_gpu,fl_gpu,coeff_deriv1_gpu,dzitdz_gpu,ep_ord_change_gpu) private(el,er,evmax,fk,gp,gm)
    do k = eul_kmin-1,eul_kmax
      do j = 1,ny
        do i = 1,nx
          ishk = 0
          do kk=k-weno_scheme+1,k+weno_scheme
            if (w_aux_gpu(i,j,kk,8) > sensor_threshold) ishk = 1
          enddo

          if (ishk == 0) then
            ft1 = 0._rkind
            ft2 = 0._rkind
            ft3 = 0._rkind
            ft4 = 0._rkind
            ft5 = 0._rkind
            ft6 = 0._rkind
            lmax = max(lmax_base+ep_ord_change_gpu(i,j,k,3),1)
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

                  rhoi = w_aux_gpu(i,j,k-m,1)
                  uui = w_aux_gpu(i,j,k-m,2)
                  vvi = w_aux_gpu(i,j,k-m,3)
                  wwi = w_aux_gpu(i,j,k-m,4)
                  enti = w_aux_gpu(i,j,k-m,5)
                  tti = w_aux_gpu(i,j,k-m,6)
                  ppi = tti*rhoi*rgas0
                  eei = enti-ppi/rhoi-0.5_rkind*(uui*uui+vvi*vvi+wwi*wwi)

                  rhoip = w_aux_gpu(i,j,k-m+l,1)
                  uuip = w_aux_gpu(i,j,k-m+l,2)
                  vvip = w_aux_gpu(i,j,k-m+l,3)
                  wwip = w_aux_gpu(i,j,k-m+l,4)
                  entip = w_aux_gpu(i,j,k-m+l,5)
                  ttip = w_aux_gpu(i,j,k-m+l,6)
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
                ft1 = ft1 + coeff_deriv1_gpu(l,lmax)*uvs1
                ft2 = ft2 + coeff_deriv1_gpu(l,lmax)*uvs2
                ft3 = ft3 + coeff_deriv1_gpu(l,lmax)*uvs3
                ft4 = ft4 + coeff_deriv1_gpu(l,lmax)*uvs4
                ft5 = ft5 + coeff_deriv1_gpu(l,lmax)*(uvs5_i+uvs5_k+uvs5_p)
                ft6 = ft6 + coeff_deriv1_gpu(l,lmax)*uvs6
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

                  rhoi = w_aux_gpu(i,j,k-m,1)
                  uui = w_aux_gpu(i,j,k-m,2)
                  vvi = w_aux_gpu(i,j,k-m,3)
                  wwi = w_aux_gpu(i,j,k-m,4)
                  enti = w_aux_gpu(i,j,k-m,5)
                  tti = w_aux_gpu(i,j,k-m,6)
                  ppi = tti*rhoi*rgas0

                  rhoip = w_aux_gpu(i,j,k-m+l,1)
                  uuip = w_aux_gpu(i,j,k-m+l,2)
                  vvip = w_aux_gpu(i,j,k-m+l,3)
                  wwip = w_aux_gpu(i,j,k-m+l,4)
                  entip = w_aux_gpu(i,j,k-m+l,5)
                  ttip = w_aux_gpu(i,j,k-m+l,6)
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
                ft1 = ft1 + coeff_deriv1_gpu(l,lmax)*uvs1
                ft2 = ft2 + coeff_deriv1_gpu(l,lmax)*uvs2
                ft3 = ft3 + coeff_deriv1_gpu(l,lmax)*uvs3
                ft4 = ft4 + coeff_deriv1_gpu(l,lmax)*uvs4
                ft5 = ft5 + coeff_deriv1_gpu(l,lmax)*uvs5
                ft6 = ft6 + coeff_deriv1_gpu(l,lmax)*uvs6
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

            fhat_gpu(i,j,k,1) = fh1
            fhat_gpu(i,j,k,2) = fh2
            fhat_gpu(i,j,k,3) = fh3
            fhat_gpu(i,j,k,4) = fh4
            fhat_gpu(i,j,k,5) = fh5
          else
            call compute_roe_average(nx, ny, nz, ng, i, i, j, j, k, k+1, w_aux_gpu, rgas0, b1, b2,&
            & b3, c, ci, h, uu, vv, ww, cp_coeff_gpu, indx_cp_l, indx_cp_r, calorically_perfect, tol_iter_nr,t0)

            call eigenvectors_z(b1, b2, b3, uu, vv, ww, c, ci, h, el, er)

            do m=1,5
              evmax(m) = -1._rkind
            enddo
            do l=1,weno_size
              ll = k + l - weno_scheme
              ww = w_aux_gpu(i,j,ll,4)
              tt = w_aux_gpu(i,j,ll,6)
              gamloc = get_gamloc_dev(indx_cp_l,indx_cp_r,tt,t0,cp_coeff_gpu,calorically_perfect,rgas0)
              c = sqrt (gamloc*rgas0*tt)
              evmax(1) = max(abs(ww-c),evmax(1))
              evmax(2) = max(abs(ww ),evmax(2))
              evmax(3) = evmax(2)
              evmax(4) = evmax(2)
              evmax(5) = max(abs(ww+c),evmax(5))
            enddo
            do l=1,weno_size
              ll = k + l - weno_scheme

              rho = w_aux_gpu(i,j,ll,1)
              uu = w_aux_gpu(i,j,ll,2)
              vv = w_aux_gpu(i,j,ll,3)
              ww = w_aux_gpu(i,j,ll,4)
              h = w_aux_gpu(i,j,ll,5)
              rhow = rho*ww
              pp = rho*w_aux_gpu(i,j,ll,6)*rgas0
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
            wenorec_ord = max(weno_scheme+ep_ord_change_gpu(i,j,k,3),1)
            call wenorec_1d(nv,gp,gm,fk,weno_scheme,wenorec_ord,weno_version,rho0,u0)
            do m=1,5
              fhat_gpu(i,j,k,m) = 0._rkind
              do mm=1,5
                fhat_gpu(i,j,k,m) = fhat_gpu(i,j,k,m) + er(mm,m) * fk(mm)
              enddo
            enddo

          endif
        enddo
      enddo
    enddo

    !$omp end target data
  endsubroutine euler_z_hybrid_kernel


  subroutine euler_z_hybrid_rusanov_kernel(nv,nv_aux,nx,ny,nz,ng,eul_kmin,eul_kmax,lmax_base,nkeep,&
  &rgas0,w_aux_gpu,fl_gpu,coeff_deriv1_gpu,dzitdz_gpu,fhat_gpu,force_zero_flux_min,force_zero_flux_max,&
  &weno_scheme,weno_version,sensor_threshold,weno_size,cp_coeff_gpu,indx_cp_l,indx_cp_r,&
  &ep_ord_change_gpu,calorically_perfect,tol_iter_nr,rho0,u0,t0)
    integer, value :: nv, nx, ny, nz, ng, nv_aux
    integer, value :: eul_kmin, eul_kmax, lmax_base, nkeep, indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv_aux) :: w_aux_gpu
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv) :: fhat_gpu
    integer, dimension(0:nx,0:ny,0:nz,1:3) :: ep_ord_change_gpu
    real(rkind), dimension(indx_cp_l:indx_cp_r+1) :: cp_coeff_gpu
    real(rkind), dimension(nx,ny,nz,nv) :: fl_gpu
    real(rkind), dimension(4,4) :: coeff_deriv1_gpu
    real(rkind), dimension(nz) :: dzitdz_gpu
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
    !$omp target data map(alloc:fk,gp,gm)
    !$omp target teams distribute parallel do collapse(3) has_device_addr(w_aux_gpu,fhat_gpu,cp_coeff_gpu,fl_gpu,coeff_deriv1_gpu,dzitdz_gpu,ep_ord_change_gpu) private(fk,gp,gm)
    do k = eul_kmin-1,eul_kmax
      do j = 1,ny
        do i = 1,nx
          ishk = 0
          do kk=k-weno_scheme+1,k+weno_scheme
            if (w_aux_gpu(i,j,kk,8) > sensor_threshold) ishk = 1
          enddo

          if (ishk == 0) then
            ft1 = 0._rkind
            ft2 = 0._rkind
            ft3 = 0._rkind
            ft4 = 0._rkind
            ft5 = 0._rkind
            ft6 = 0._rkind
            lmax = max(lmax_base+ep_ord_change_gpu(i,j,k,3),1)
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

                  rhoi = w_aux_gpu(i,j,k-m,1)
                  uui = w_aux_gpu(i,j,k-m,2)
                  vvi = w_aux_gpu(i,j,k-m,3)
                  wwi = w_aux_gpu(i,j,k-m,4)
                  enti = w_aux_gpu(i,j,k-m,5)
                  tti = w_aux_gpu(i,j,k-m,6)
                  ppi = tti*rhoi*rgas0
                  eei = enti-ppi/rhoi-0.5_rkind*(uui*uui+vvi*vvi+wwi*wwi)

                  rhoip = w_aux_gpu(i,j,k-m+l,1)
                  uuip = w_aux_gpu(i,j,k-m+l,2)
                  vvip = w_aux_gpu(i,j,k-m+l,3)
                  wwip = w_aux_gpu(i,j,k-m+l,4)
                  entip = w_aux_gpu(i,j,k-m+l,5)
                  ttip = w_aux_gpu(i,j,k-m+l,6)
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
                ft1 = ft1 + coeff_deriv1_gpu(l,lmax)*uvs1
                ft2 = ft2 + coeff_deriv1_gpu(l,lmax)*uvs2
                ft3 = ft3 + coeff_deriv1_gpu(l,lmax)*uvs3
                ft4 = ft4 + coeff_deriv1_gpu(l,lmax)*uvs4
                ft5 = ft5 + coeff_deriv1_gpu(l,lmax)*(uvs5_i+uvs5_k+uvs5_p)
                ft6 = ft6 + coeff_deriv1_gpu(l,lmax)*uvs6
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

                  rhoi = w_aux_gpu(i,j,k-m,1)
                  uui = w_aux_gpu(i,j,k-m,2)
                  vvi = w_aux_gpu(i,j,k-m,3)
                  wwi = w_aux_gpu(i,j,k-m,4)
                  enti = w_aux_gpu(i,j,k-m,5)
                  tti = w_aux_gpu(i,j,k-m,6)
                  ppi = tti*rhoi*rgas0

                  rhoip = w_aux_gpu(i,j,k-m+l,1)
                  uuip = w_aux_gpu(i,j,k-m+l,2)
                  vvip = w_aux_gpu(i,j,k-m+l,3)
                  wwip = w_aux_gpu(i,j,k-m+l,4)
                  entip = w_aux_gpu(i,j,k-m+l,5)
                  ttip = w_aux_gpu(i,j,k-m+l,6)
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
                ft1 = ft1 + coeff_deriv1_gpu(l,lmax)*uvs1
                ft2 = ft2 + coeff_deriv1_gpu(l,lmax)*uvs2
                ft3 = ft3 + coeff_deriv1_gpu(l,lmax)*uvs3
                ft4 = ft4 + coeff_deriv1_gpu(l,lmax)*uvs4
                ft5 = ft5 + coeff_deriv1_gpu(l,lmax)*uvs5
                ft6 = ft6 + coeff_deriv1_gpu(l,lmax)*uvs6
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

            fhat_gpu(i,j,k,1) = fh1
            fhat_gpu(i,j,k,2) = fh2
            fhat_gpu(i,j,k,3) = fh3
            fhat_gpu(i,j,k,4) = fh4
            fhat_gpu(i,j,k,5) = fh5
          else
            evmax = -1._rkind
            do l=1,weno_size
              ll = k + l - weno_scheme
              ww = w_aux_gpu(i,j,ll,4)
              tt = w_aux_gpu(i,j,ll,6)
              gamloc = get_gamloc_dev(indx_cp_l,indx_cp_r,tt,t0,cp_coeff_gpu,calorically_perfect,rgas0)
              c = sqrt (gamloc*rgas0*tt)
              evm = max(abs(ww-c),abs(ww+c))
              evmax = max(evm,evmax)
            enddo
            do l=1,weno_size
              ll = k + l - weno_scheme

              rho = w_aux_gpu(i,j,ll,1)
              uu = w_aux_gpu(i,j,ll,2)
              vv = w_aux_gpu(i,j,ll,3)
              ww = w_aux_gpu(i,j,ll,4)
              h = w_aux_gpu(i,j,ll,5)
              rhow = rho*ww
              pp = rho*w_aux_gpu(i,j,ll,6)*rgas0

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
            wenorec_ord = max(weno_scheme+ep_ord_change_gpu(i,j,k,3),1)
            call wenorec_1d_rusanov(nv,gp,gm,fk,weno_scheme,wenorec_ord,weno_version,rho0,u0)
            do m=1,5
              fhat_gpu(i,j,k,m) = fk(m)
            enddo

          endif
        enddo
      enddo
    enddo

    !$omp end target data
  endsubroutine euler_z_hybrid_rusanov_kernel


  subroutine euler_y_hybrid_c2_kernel(nv,nv_aux,nx,ny,nz,ng,eul_jmin,eul_jmax,lmax_base,nkeep,rgas0,&
  &w_aux_gpu,fl_gpu,coeff_deriv1_gpu,fhat_gpu,detadxc2_gpu,detadxnc2_gpu,detadyc2_gpu,detadync2_gpu,&
  &jac_gpu,metajac1_gpu,wall_tag_gpu,force_zero_flux_min,force_zero_flux_max,weno_scheme,weno_version,&
  &sensor_threshold,weno_size,cp_coeff_gpu,indx_cp_l,indx_cp_r,ep_ord_change_gpu,calorically_perfect,&
  &tol_iter_nr,rho0,u0,t0)
    integer, value :: nv, nx, ny, nz, ng, nv_aux
    integer, value :: eul_jmin, eul_jmax, lmax_base, nkeep, indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv_aux) :: w_aux_gpu
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv) :: fhat_gpu
    integer, dimension(1-ng:nx+ng) :: wall_tag_gpu
    integer, dimension(0:nx,0:ny,0:nz,1:3) :: ep_ord_change_gpu
    real(rkind), dimension(indx_cp_l:indx_cp_r+1) :: cp_coeff_gpu
    real(rkind), dimension(nx,ny,nz,nv) :: fl_gpu
    real(rkind), dimension(4,4) :: coeff_deriv1_gpu
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng) :: detadxc2_gpu
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng) :: detadxnc2_gpu
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng) :: detadyc2_gpu
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng) :: detadync2_gpu
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng) :: jac_gpu
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng) :: metajac1_gpu
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
    !$omp target data map(alloc:el,er,evmax,fj,gp,gm,eler)
    !$omp target teams distribute parallel do collapse(3) has_device_addr(w_aux_gpu,fhat_gpu,cp_coeff_gpu,fl_gpu,coeff_deriv1_gpu,detadxc2_gpu,detadxnc2_gpu,detadyc2_gpu,detadync2_gpu,jac_gpu,metajac1_gpu,wall_tag_gpu,ep_ord_change_gpu) private(el,er,evmax,fj,gp,gm,eler)
    do k = 1,nz
      do j = 0,eul_jmax
        do i = 1,nx
          lmax = lmax_base
          weno_scheme_j = weno_scheme
          if (j <= lmax) then
            if (wall_tag_gpu(i) < 1) then
              lmax = max(lmax_base+ep_ord_change_gpu(i,j,k,2),1)
              weno_scheme_j = max(weno_scheme+ep_ord_change_gpu(i,j,k,2),1)
            endif
          else
            lmax = max(lmax_base+ep_ord_change_gpu(i,j,k,2),1)
            weno_scheme_j = max(weno_scheme+ep_ord_change_gpu(i,j,k,2),1)
          endif
          weno_size = 2*weno_scheme_j

          ishk = 0
          do jj=j-weno_scheme_j+1,j+weno_scheme_j
            if (w_aux_gpu(i,jj,k,8) > sensor_threshold) ishk = 1
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

                rhoi = w_aux_gpu(i,j-m,k,1)
                uui = w_aux_gpu(i,j-m,k,2)
                vvi = w_aux_gpu(i,j-m,k,3)
                wwi = w_aux_gpu(i,j-m,k,4)
                enti = w_aux_gpu(i,j-m,k,5)
                tti = w_aux_gpu(i,j-m,k,6)
                ppi = tti*rhoi*rgas0
                uti = (uui*detadxc2_gpu(i,j-m)+vvi*detadyc2_gpu(i,j-m))/jac_gpu(i,j-m)
                detadxi = detadxc2_gpu(i,j-m)/jac_gpu(i,j-m)
                detadyi = detadyc2_gpu(i,j-m)/jac_gpu(i,j-m)

                rhoip = w_aux_gpu(i,j-m+l,k,1)
                uuip = w_aux_gpu(i,j-m+l,k,2)
                vvip = w_aux_gpu(i,j-m+l,k,3)
                wwip = w_aux_gpu(i,j-m+l,k,4)
                entip = w_aux_gpu(i,j-m+l,k,5)
                ttip = w_aux_gpu(i,j-m+l,k,6)
                ppip = ttip*rhoip*rgas0
                utip = (uuip*detadxc2_gpu(i,j-m+l)+vvip*detadyc2_gpu(i,j-m+l))/jac_gpu(i,j-m+l)
                detadxip = detadxc2_gpu(i,j-m+l)/jac_gpu(i,j-m+l)
                detadyip = detadyc2_gpu(i,j-m+l)/jac_gpu(i,j-m+l)

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
              ft1 = ft1 + coeff_deriv1_gpu(l,lmax)*uvs1
              ft2 = ft2 + coeff_deriv1_gpu(l,lmax)*uvs2
              ft3 = ft3 + coeff_deriv1_gpu(l,lmax)*uvs3
              ft4 = ft4 + coeff_deriv1_gpu(l,lmax)*uvs4
              ft5 = ft5 + coeff_deriv1_gpu(l,lmax)*uvs5
              ft6 = ft6 + coeff_deriv1_gpu(l,lmax)*uvs6
              ft7 = ft7 + coeff_deriv1_gpu(l,lmax)*uvs7
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

            fhat_gpu(i,j,k,1) = fh1
            fhat_gpu(i,j,k,2) = fh2
            fhat_gpu(i,j,k,3) = fh3
            fhat_gpu(i,j,k,4) = fh4
            fhat_gpu(i,j,k,5) = fh5
          else
            call compute_roe_average(nx, ny, nz, ng, i, i, j, j+1, k, k, w_aux_gpu, rgas0, b1, b2,&
            & b3, c, ci, h, uu, vv, ww, cp_coeff_gpu, indx_cp_l, indx_cp_r, calorically_perfect, tol_iter_nr,t0)
            detadxav = 0.5_rkind * (detadxnc2_gpu(i,j) + detadxnc2_gpu(i,j+1))
            detadyav = 0.5_rkind * (detadync2_gpu(i,j) + detadync2_gpu(i,j+1))
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
              uu = w_aux_gpu(i,ll,k,2)
              vv = w_aux_gpu(i,ll,k,3)
              ut = detadxnc2_gpu(i,ll)*uu+detadync2_gpu(i,ll)*vv
              tt = w_aux_gpu(i,ll,k,6)
              gamloc = get_gamloc_dev(indx_cp_l,indx_cp_r,tt,t0,cp_coeff_gpu,calorically_perfect,rgas0)
              c = sqrt (gamloc*rgas0*tt)
              evmax(1) = max(abs(ut-c),evmax(1))
              evmax(2) = max(abs(ut ),evmax(2))
              evmax(3) = evmax(2)
              evmax(4) = evmax(2)
              evmax(5) = max(abs(ut+c),evmax(5))
            enddo
            do l=1,weno_size
              ll = j + l - weno_scheme_j

              rho = w_aux_gpu(i,ll,k,1)
              uu = w_aux_gpu(i,ll,k,2)
              vv = w_aux_gpu(i,ll,k,3)
              ut = detadxnc2_gpu(i,ll)*uu+detadync2_gpu(i,ll)*vv
              ww = w_aux_gpu(i,ll,k,4)
              h = w_aux_gpu(i,ll,k,5)
              rhov = rho*vv
              pp = rho*w_aux_gpu(i,ll,k,6)*rgas0
              fj(1) = rho * ut
              fj(2) = (rho * uu * ut + pp * detadxnc2_gpu(i,ll))
              fj(3) = (rho * vv * ut + pp * detadync2_gpu(i,ll))
              fj(4) = rho * ww * ut
              fj(5) = rho * h * ut
              do m=1,5
                wc = 0._rkind
                gc = 0._rkind

                wc = wc + el(1,m) * rho * metajac1_gpu(i,ll)
                gc = gc + el(1,m) * fj(1) * metajac1_gpu(i,ll)
                wc = wc + el(2,m) * rho*uu * metajac1_gpu(i,ll)
                gc = gc + el(2,m) * fj(2) * metajac1_gpu(i,ll)
                wc = wc + el(3,m) * rho*vv * metajac1_gpu(i,ll)
                gc = gc + el(3,m) * fj(3) * metajac1_gpu(i,ll)
                wc = wc + el(4,m) * rho*ww * metajac1_gpu(i,ll)
                gc = gc + el(4,m) * fj(4) * metajac1_gpu(i,ll)
                wc = wc + el(5,m) * (rho*h-pp) * metajac1_gpu(i,ll)
                gc = gc + el(5,m) * fj(5) * metajac1_gpu(i,ll)

                c = 0.5_rkind * (gc + evmax(m) * wc)
                gp(m,l) = c
                gm(m,l) = gc - c
              enddo
            enddo
            call wenorec_1d(nv,gp,gm,fj,weno_scheme_j,weno_scheme_j,weno_version,rho0,u0)
            do m=1,5
              fhat_gpu(i,j,k,m) = 0._rkind
              do mm=1,5
                fhat_gpu(i,j,k,m) = fhat_gpu(i,j,k,m) + er(mm,m) * fj(mm)
              enddo
            enddo

          endif
        enddo
      enddo
    enddo

    !$omp end target data
  endsubroutine euler_y_hybrid_c2_kernel



  subroutine euler_y_hybrid_kernel(nv,nv_aux,nx,ny,nz,ng,eul_jmin,eul_jmax,lmax_base,nkeep,rgas0,&
  &w_aux_gpu,fl_gpu,coeff_deriv1_gpu,detady_gpu,fhat_gpu,force_zero_flux_min,force_zero_flux_max,&
  &weno_scheme,weno_version,sensor_threshold,weno_size,cp_coeff_gpu,indx_cp_l,indx_cp_r,&
  &ep_ord_change_gpu,calorically_perfect,tol_iter_nr,rho0,u0,t0)
    integer, value :: nv, nx, ny, nz, ng, nv_aux
    integer, value :: eul_jmin, eul_jmax, lmax_base, nkeep, indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv_aux) :: w_aux_gpu
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv) :: fhat_gpu
    integer, dimension(0:nx,0:ny,0:nz,1:3) :: ep_ord_change_gpu
    real(rkind), dimension(indx_cp_l:indx_cp_r+1) :: cp_coeff_gpu
    real(rkind), dimension(nx,ny,nz,nv) :: fl_gpu
    real(rkind), dimension(4,4) :: coeff_deriv1_gpu
    real(rkind), dimension(ny) :: detady_gpu
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
    !$omp target data map(alloc:el,er,evmax,fj,gp,gm)
    !$omp target teams distribute parallel do collapse(3) has_device_addr(w_aux_gpu,fhat_gpu,cp_coeff_gpu,fl_gpu,coeff_deriv1_gpu,detady_gpu,ep_ord_change_gpu) private(el,er,evmax,fj,gp,gm)
    do k = 1,nz
      do j = eul_jmin-1,eul_jmax
        do i = 1,nx
          ishk = 0
          do jj=j-weno_scheme+1,j+weno_scheme
            if (w_aux_gpu(i,jj,k,8) > sensor_threshold) ishk = 1
          enddo

          if (ishk == 0) then
            ft1 = 0._rkind
            ft2 = 0._rkind
            ft3 = 0._rkind
            ft4 = 0._rkind
            ft5 = 0._rkind
            ft6 = 0._rkind
            lmax = max(lmax_base+ep_ord_change_gpu(i,j,k,2),1)
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

                  rhoi = w_aux_gpu(i,j-m,k,1)
                  uui = w_aux_gpu(i,j-m,k,2)
                  vvi = w_aux_gpu(i,j-m,k,3)
                  wwi = w_aux_gpu(i,j-m,k,4)
                  enti = w_aux_gpu(i,j-m,k,5)
                  tti = w_aux_gpu(i,j-m,k,6)
                  ppi = tti*rhoi*rgas0
                  eei = enti-ppi/rhoi-0.5_rkind*(uui*uui+vvi*vvi+wwi*wwi)

                  rhoip = w_aux_gpu(i,j-m+l,k,1)
                  uuip = w_aux_gpu(i,j-m+l,k,2)
                  vvip = w_aux_gpu(i,j-m+l,k,3)
                  wwip = w_aux_gpu(i,j-m+l,k,4)
                  entip = w_aux_gpu(i,j-m+l,k,5)
                  ttip = w_aux_gpu(i,j-m+l,k,6)
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
                ft1 = ft1 + coeff_deriv1_gpu(l,lmax)*uvs1
                ft2 = ft2 + coeff_deriv1_gpu(l,lmax)*uvs2
                ft3 = ft3 + coeff_deriv1_gpu(l,lmax)*uvs3
                ft4 = ft4 + coeff_deriv1_gpu(l,lmax)*uvs4
                ft5 = ft5 + coeff_deriv1_gpu(l,lmax)*(uvs5_i+uvs5_k+uvs5_p)
                ft6 = ft6 + coeff_deriv1_gpu(l,lmax)*uvs6
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

                  rhoi = w_aux_gpu(i,j-m,k,1)
                  uui = w_aux_gpu(i,j-m,k,2)
                  vvi = w_aux_gpu(i,j-m,k,3)
                  wwi = w_aux_gpu(i,j-m,k,4)
                  enti = w_aux_gpu(i,j-m,k,5)
                  tti = w_aux_gpu(i,j-m,k,6)
                  ppi = tti*rhoi*rgas0

                  rhoip = w_aux_gpu(i,j-m+l,k,1)
                  uuip = w_aux_gpu(i,j-m+l,k,2)
                  vvip = w_aux_gpu(i,j-m+l,k,3)
                  wwip = w_aux_gpu(i,j-m+l,k,4)
                  entip = w_aux_gpu(i,j-m+l,k,5)
                  ttip = w_aux_gpu(i,j-m+l,k,6)
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
                ft1 = ft1 + coeff_deriv1_gpu(l,lmax)*uvs1
                ft2 = ft2 + coeff_deriv1_gpu(l,lmax)*uvs2
                ft3 = ft3 + coeff_deriv1_gpu(l,lmax)*uvs3
                ft4 = ft4 + coeff_deriv1_gpu(l,lmax)*uvs4
                ft5 = ft5 + coeff_deriv1_gpu(l,lmax)*uvs5
                ft6 = ft6 + coeff_deriv1_gpu(l,lmax)*uvs6
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

            fhat_gpu(i,j,k,1) = fh1
            fhat_gpu(i,j,k,2) = fh2
            fhat_gpu(i,j,k,3) = fh3
            fhat_gpu(i,j,k,4) = fh4
            fhat_gpu(i,j,k,5) = fh5
          else
            call compute_roe_average(nx, ny, nz, ng, i, i, j, j+1, k, k, w_aux_gpu, rgas0, b1, b2,&
            & b3, c, ci, h, uu, vv, ww, cp_coeff_gpu, indx_cp_l, indx_cp_r, calorically_perfect, tol_iter_nr,t0)

            call eigenvectors_y(b1, b2, b3, uu, vv, ww, c, ci, h, el, er)

            do m=1,5
              evmax(m) = -1._rkind
            enddo
            do l=1,weno_size
              ll = j + l - weno_scheme
              vv = w_aux_gpu(i,ll,k,3)
              tt = w_aux_gpu(i,ll,k,6)
              gamloc = get_gamloc_dev(indx_cp_l,indx_cp_r,tt,t0,cp_coeff_gpu,calorically_perfect,rgas0)
              c = sqrt (gamloc*rgas0*tt)
              evmax(1) = max(abs(vv-c),evmax(1))
              evmax(2) = max(abs(vv ),evmax(2))
              evmax(3) = evmax(2)
              evmax(4) = evmax(2)
              evmax(5) = max(abs(vv+c),evmax(5))
            enddo
            do l=1,weno_size
              ll = j + l - weno_scheme

              rho = w_aux_gpu(i,ll,k,1)
              uu = w_aux_gpu(i,ll,k,2)
              vv = w_aux_gpu(i,ll,k,3)
              ww = w_aux_gpu(i,ll,k,4)
              h = w_aux_gpu(i,ll,k,5)
              rhov = rho*vv
              pp = rho*w_aux_gpu(i,ll,k,6)*rgas0
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
            wenorec_ord = max(weno_scheme+ep_ord_change_gpu(i,j,k,2),1)
            call wenorec_1d(nv,gp,gm,fj,weno_scheme,wenorec_ord,weno_version,rho0,u0)
            do m=1,5
              fhat_gpu(i,j,k,m) = 0._rkind
              do mm=1,5
                fhat_gpu(i,j,k,m) = fhat_gpu(i,j,k,m) + er(mm,m) * fj(mm)
              enddo
            enddo

          endif
        enddo
      enddo
    enddo

    !$omp end target data
  endsubroutine euler_y_hybrid_kernel


  subroutine euler_y_hybrid_rusanov_kernel(nv,nv_aux,nx,ny,nz,ng,eul_jmin,eul_jmax,lmax_base,nkeep,&
  &rgas0,w_aux_gpu,fl_gpu,coeff_deriv1_gpu,detady_gpu,fhat_gpu,force_zero_flux_min,force_zero_flux_max,&
  &weno_scheme,weno_version,sensor_threshold,weno_size,cp_coeff_gpu,indx_cp_l,indx_cp_r,&
  &ep_ord_change_gpu,calorically_perfect,tol_iter_nr,rho0,u0,t0)
    integer, value :: nv, nx, ny, nz, ng, nv_aux
    integer, value :: eul_jmin, eul_jmax, lmax_base, nkeep, indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv_aux) :: w_aux_gpu
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv) :: fhat_gpu
    integer, dimension(0:nx,0:ny,0:nz,1:3) :: ep_ord_change_gpu
    real(rkind), dimension(indx_cp_l:indx_cp_r+1) :: cp_coeff_gpu
    real(rkind), dimension(nx,ny,nz,nv) :: fl_gpu
    real(rkind), dimension(4,4) :: coeff_deriv1_gpu
    real(rkind), dimension(ny) :: detady_gpu
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
    !$omp target data map(alloc:fj,gp,gm)
    !$omp target teams distribute parallel do collapse(3) has_device_addr(w_aux_gpu,fhat_gpu,cp_coeff_gpu,fl_gpu,coeff_deriv1_gpu,detady_gpu,ep_ord_change_gpu) private(fj,gp,gm)
    do k = 1,nz
      do j = eul_jmin-1,eul_jmax
        do i = 1,nx
          ishk = 0
          do jj=j-weno_scheme+1,j+weno_scheme
            if (w_aux_gpu(i,jj,k,8) > sensor_threshold) ishk = 1
          enddo

          if (ishk == 0) then
            ft1 = 0._rkind
            ft2 = 0._rkind
            ft3 = 0._rkind
            ft4 = 0._rkind
            ft5 = 0._rkind
            ft6 = 0._rkind
            lmax = max(lmax_base+ep_ord_change_gpu(i,j,k,2),1)
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

                  rhoi = w_aux_gpu(i,j-m,k,1)
                  uui = w_aux_gpu(i,j-m,k,2)
                  vvi = w_aux_gpu(i,j-m,k,3)
                  wwi = w_aux_gpu(i,j-m,k,4)
                  enti = w_aux_gpu(i,j-m,k,5)
                  tti = w_aux_gpu(i,j-m,k,6)
                  ppi = tti*rhoi*rgas0
                  eei = enti-ppi/rhoi-0.5_rkind*(uui*uui+vvi*vvi+wwi*wwi)

                  rhoip = w_aux_gpu(i,j-m+l,k,1)
                  uuip = w_aux_gpu(i,j-m+l,k,2)
                  vvip = w_aux_gpu(i,j-m+l,k,3)
                  wwip = w_aux_gpu(i,j-m+l,k,4)
                  entip = w_aux_gpu(i,j-m+l,k,5)
                  ttip = w_aux_gpu(i,j-m+l,k,6)
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
                ft1 = ft1 + coeff_deriv1_gpu(l,lmax)*uvs1
                ft2 = ft2 + coeff_deriv1_gpu(l,lmax)*uvs2
                ft3 = ft3 + coeff_deriv1_gpu(l,lmax)*uvs3
                ft4 = ft4 + coeff_deriv1_gpu(l,lmax)*uvs4
                ft5 = ft5 + coeff_deriv1_gpu(l,lmax)*(uvs5_i+uvs5_k+uvs5_p)
                ft6 = ft6 + coeff_deriv1_gpu(l,lmax)*uvs6
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

                  rhoi = w_aux_gpu(i,j-m,k,1)
                  uui = w_aux_gpu(i,j-m,k,2)
                  vvi = w_aux_gpu(i,j-m,k,3)
                  wwi = w_aux_gpu(i,j-m,k,4)
                  enti = w_aux_gpu(i,j-m,k,5)
                  tti = w_aux_gpu(i,j-m,k,6)
                  ppi = tti*rhoi*rgas0

                  rhoip = w_aux_gpu(i,j-m+l,k,1)
                  uuip = w_aux_gpu(i,j-m+l,k,2)
                  vvip = w_aux_gpu(i,j-m+l,k,3)
                  wwip = w_aux_gpu(i,j-m+l,k,4)
                  entip = w_aux_gpu(i,j-m+l,k,5)
                  ttip = w_aux_gpu(i,j-m+l,k,6)
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
                ft1 = ft1 + coeff_deriv1_gpu(l,lmax)*uvs1
                ft2 = ft2 + coeff_deriv1_gpu(l,lmax)*uvs2
                ft3 = ft3 + coeff_deriv1_gpu(l,lmax)*uvs3
                ft4 = ft4 + coeff_deriv1_gpu(l,lmax)*uvs4
                ft5 = ft5 + coeff_deriv1_gpu(l,lmax)*uvs5
                ft6 = ft6 + coeff_deriv1_gpu(l,lmax)*uvs6
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

            fhat_gpu(i,j,k,1) = fh1
            fhat_gpu(i,j,k,2) = fh2
            fhat_gpu(i,j,k,3) = fh3
            fhat_gpu(i,j,k,4) = fh4
            fhat_gpu(i,j,k,5) = fh5
          else
            evmax = -1._rkind
            do l=1,weno_size
              ll = j + l - weno_scheme
              vv = w_aux_gpu(i,ll,k,3)
              tt = w_aux_gpu(i,ll,k,6)
              gamloc = get_gamloc_dev(indx_cp_l,indx_cp_r,tt,t0,cp_coeff_gpu,calorically_perfect,rgas0)
              c = sqrt (gamloc*rgas0*tt)
              evm = max(abs(vv-c),abs(vv+c))
              evmax = max(evm,evmax)
            enddo
            do l=1,weno_size
              ll = j + l - weno_scheme

              rho = w_aux_gpu(i,ll,k,1)
              uu = w_aux_gpu(i,ll,k,2)
              vv = w_aux_gpu(i,ll,k,3)
              ww = w_aux_gpu(i,ll,k,4)
              h = w_aux_gpu(i,ll,k,5)
              rhov = rho*vv
              pp = rho*w_aux_gpu(i,ll,k,6)*rgas0

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
            wenorec_ord = max(weno_scheme+ep_ord_change_gpu(i,j,k,2),1)
            call wenorec_1d_rusanov(nv,gp,gm,fj,weno_scheme,wenorec_ord,weno_version,rho0,u0)
            do m=1,5
              fhat_gpu(i,j,k,m) = fj(m)
            enddo

          endif
        enddo
      enddo
    enddo

    !$omp end target data
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
    !$omp declare target
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
    !$omp declare target
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
    !$omp declare target
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
    !$omp declare target
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
    !$omp declare target
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
    !$omp declare target
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
    !$omp declare target
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


  subroutine compute_roe_average(nx,ny,nz,ng,i,ip,j,jp,k,kp,w_aux_gpu,rgas0,b1,b2,b3,c,ci,h,uu,vv,&
  &ww,cp_coeff_gpu,indx_cp_l,indx_cp_r,calorically_perfect,tol_iter_nr,t0)
    integer :: ng,i,ip,j,jp,k,kp,indx_cp_l,indx_cp_r,calorically_perfect
    integer :: nx,ny,nz
    real(rkind), intent(in) :: rgas0, tol_iter_nr,t0
    real(rkind), intent(out) :: b1, b2, b3, uu, vv, ww, ci, h, c
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in) :: w_aux_gpu
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in) :: cp_coeff_gpu
    real(rkind) :: up, vp, wp, qqp, hp, r, rp1, cc, qq, gam, gm1
    integer :: ll,iter,max_iter
    real(rkind) :: tt,gm1loc,hbar,ttp,told,num,den,gamloc,cploc,tpow,tpowp,p_rho,p_e,etot,rho
    !$omp declare target
    max_iter = 50
    uu = w_aux_gpu(i,j,k,2)
    vv = w_aux_gpu(i,j,k,3)
    ww = w_aux_gpu(i,j,k,4)
    qq = 0.5_rkind * (uu*uu +vv*vv + ww*ww)
    h = w_aux_gpu(i,j,k,5)
    up = w_aux_gpu(ip,jp,kp,2)
    vp = w_aux_gpu(ip,jp,kp,3)
    wp = w_aux_gpu(ip,jp,kp,4)
    qqp = 0.5_rkind * (up*up +vp*vp +wp*wp)
    hp = w_aux_gpu(ip,jp,kp,5)
    r = w_aux_gpu(ip,jp,kp,1)/w_aux_gpu(i,j,k,1)
    r = sqrt(r)
    rho = r*w_aux_gpu(i,j,k,1)
    rp1 = 1._rkind/(r +1._rkind)
    uu = (r*up +uu)*rp1
    vv = (r*vp +vv)*rp1
    ww = (r*wp +ww)*rp1
    h = (r*hp +h)*rp1
    qq = 0.5_rkind * (uu*uu +vv*vv + ww*ww)
    hbar = h - qq - cp_coeff_gpu(indx_cp_r+1)*t0
    if (calorically_perfect==1) then
      tt = t0+hbar/cp_coeff_gpu(0)
      gamloc = cp_coeff_gpu(0)/(cp_coeff_gpu(0)-rgas0)
    else
      tt = w_aux_gpu(i ,j ,k ,6)
      ttp = w_aux_gpu(ip,jp,kp,6)
      tt = (r*ttp +tt)*rp1
      told = tt
      do iter=1,max_iter
        num = 0._rkind
        den = 0._rkind
        do ll=indx_cp_l,indx_cp_r
          if (ll==-1) then
            tpow = (told/t0)**ll
            den = den+cp_coeff_gpu(ll)*tpow
            num = num+cp_coeff_gpu(ll)*log(told/t0)
          else
            tpow = (told/t0)**ll
            tpowp = (told/t0)*tpow
            den = den+cp_coeff_gpu(ll)*tpow
            num = num+cp_coeff_gpu(ll)*(tpowp-1._rkind)/(ll+1._rkind)
          endif
        enddo
        num = num*t0
        tt = told+(hbar-num)/den
        if (abs(tt-told) < tol_iter_nr) exit
        told = tt
      enddo
      cploc = 0._rkind
      do ll=indx_cp_l,indx_cp_r
        cploc = cploc+cp_coeff_gpu(ll)*(tt/t0)**ll
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


  subroutine force_rhs_2_kernel(nx,ny,nz,ng,fln_gpu,w_aux_gpu,bulk5g,fluid_mask_gpu)
    integer :: nx, ny, nz, ng
    real(rkind), dimension(5), intent(in) :: bulk5g
    real(rkind), dimension(1:,1:,1:,1:), intent(inout) :: fln_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in) :: w_aux_gpu
    integer, dimension(1-ng:,1-ng:,1-ng:), intent(in) :: fluid_mask_gpu
    real(rkind) :: bulk_1, bulk_2, uu
    integer :: i,j,k,iercuda
    bulk_1 = bulk5g(1)
    bulk_2 = bulk5g(2)
    !$omp target teams distribute parallel do collapse(3) has_device_addr(fln_gpu,w_aux_gpu,fluid_mask_gpu)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          if (fluid_mask_gpu(i,j,k)==0) then
            uu = w_aux_gpu(i,j,k,2)
            fln_gpu(i,j,k,1) = fln_gpu(i,j,k,1) - bulk_1
            fln_gpu(i,j,k,2) = fln_gpu(i,j,k,2) - bulk_2
            fln_gpu(i,j,k,5) = fln_gpu(i,j,k,5) - uu*bulk_2
          endif
        enddo
      enddo
    enddo


  endsubroutine force_rhs_2_kernel



  subroutine force_rhs_2_c2_kernel(nx,ny,nz,ng,fln_gpu,w_aux_gpu,bulk5g,fluid_mask_gpu,r_curv,yn_gpu,dxdcsinc2_gpu,dydcsinc2_gpu)
    integer :: nx, ny, nz, ng
    real(rkind), dimension(5), intent(in) :: bulk5g
    real(rkind), intent(in) :: r_curv
    real(rkind), dimension(1:,1:,1:,1:), intent(inout) :: fln_gpu
    real(rkind), dimension(1:), intent(in) :: yn_gpu
    real(rkind), dimension(1-ng:,1-ng:), intent(in) :: dxdcsinc2_gpu, dydcsinc2_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in) :: w_aux_gpu
    integer, dimension(1-ng:,1-ng:,1-ng:), intent(in) :: fluid_mask_gpu
    real(rkind) :: bulk_1, bulk_5, uu, vv, r, dpdcsi, dpdx, dpdy
    integer :: i,j,k,iercuda
    bulk_1 = bulk5g(1)
    bulk_5 = bulk5g(5)
    !$omp target teams distribute parallel do collapse(3) has_device_addr(fln_gpu,yn_gpu,dxdcsinc2_gpu,dydcsinc2_gpu,w_aux_gpu,fluid_mask_gpu)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          if (fluid_mask_gpu(i,j,k)==0) then
            uu = w_aux_gpu(i,j,k,2)
            vv = w_aux_gpu(i,j,k,3)
            r = r_curv+0.5_rkind*(yn_gpu(j)+yn_gpu(j+1))
            dpdcsi = bulk_5/r
            dpdx = dpdcsi*dxdcsinc2_gpu(i,j)
            dpdy = dpdcsi*dydcsinc2_gpu(i,j)
            fln_gpu(i,j,k,1) = fln_gpu(i,j,k,1) - bulk_1
            fln_gpu(i,j,k,2) = fln_gpu(i,j,k,2) - dpdx
            fln_gpu(i,j,k,3) = fln_gpu(i,j,k,3) - dpdy
            fln_gpu(i,j,k,5) = fln_gpu(i,j,k,5) - uu*dpdx - vv*dpdy
          endif
        enddo
      enddo
    enddo


  endsubroutine force_rhs_2_c2_kernel


  subroutine force_rhs_1_kernel(nx,ny,nz,ng,yn_gpu,fln_gpu,w_gpu,w_aux_gpu,bulk5,fluid_mask_gpu)
    integer, intent(in) :: nx, ny, nz, ng
    real(rkind), dimension(1:), intent(in) :: yn_gpu
    real(rkind), dimension(1:,1:,1:,1:), intent(in) :: fln_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in) :: w_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in) :: w_aux_gpu
    integer, dimension(1-ng:,1-ng:,1-ng:), intent(in) :: fluid_mask_gpu
    real(rkind), dimension(5), intent(out) :: bulk5
    real(rkind) :: bulk_1, bulk_2, bulk_3, bulk_4, bulk_5
    real(rkind) :: dy
    integer :: i,j,k,iercuda
    bulk_1 = 0._rkind
    bulk_2 = 0._rkind
    bulk_3 = 0._rkind
    bulk_4 = 0._rkind
    bulk_5 = 0._rkind
    !$omp target teams distribute parallel do collapse(3) has_device_addr(yn_gpu,fln_gpu,w_gpu,w_aux_gpu,fluid_mask_gpu) &
    !$omp& reduction(+:bulk_1,bulk_2,bulk_3,bulk_4,bulk_5)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          if (fluid_mask_gpu(i,j,k)==0) then
            dy = yn_gpu(j+1)-yn_gpu(j)
            bulk_1 = bulk_1 + fln_gpu(i,j,k,1)*dy
            bulk_2 = bulk_2 + fln_gpu(i,j,k,2)*dy
            bulk_3 = bulk_3 + w_gpu (i,j,k,1)*dy
            bulk_4 = bulk_4 + w_gpu (i,j,k,2)*dy
            bulk_5 = bulk_5 + w_gpu (i,j,k,2)*dy*w_aux_gpu(i,j,k,6)
          endif
        enddo
      enddo
    enddo
    bulk5(1) = bulk_1
    bulk5(2) = bulk_2
    bulk5(3) = bulk_3
    bulk5(4) = bulk_4
    bulk5(5) = bulk_5

  endsubroutine force_rhs_1_kernel



  subroutine force_rhs_1_c2_kernel(nx,ny,nz,ng,yn_gpu,fln_gpu,w_gpu,w_aux_gpu,bulk5,fluid_mask_gpu,&
  &dcsidxnc2_gpu,dcsidync2_gpu,jac_gpu)
    integer, intent(in) :: nx, ny, nz, ng
    real(rkind), dimension(1:), intent(in) :: yn_gpu
    real(rkind), dimension(1:,1:,1:,1:), intent(in) :: fln_gpu
    real(rkind), dimension(1-ng:,1-ng:), intent(in) :: dcsidxnc2_gpu, dcsidync2_gpu, jac_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in) :: w_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in) :: w_aux_gpu
    integer, dimension(1-ng:,1-ng:,1-ng:), intent(in) :: fluid_mask_gpu
    real(rkind), dimension(5), intent(out) :: bulk5
    real(rkind) :: bulk_1, bulk_2, bulk_3, bulk_4, bulk_5
    real(rkind) :: dv, rhou, f_rhou
    integer :: i,j,k,iercuda
    bulk_1 = 0._rkind
    bulk_2 = 0._rkind
    bulk_3 = 0._rkind
    bulk_4 = 0._rkind
    bulk_5 = 0._rkind
    !$omp target teams distribute parallel do collapse(3) has_device_addr(yn_gpu,fln_gpu,dcsidxnc2_gpu,dcsidync2_gpu,jac_gpu,w_gpu,w_aux_gpu,fluid_mask_gpu) &
    !$omp& reduction(+:bulk_1,bulk_2,bulk_3,bulk_4,bulk_5)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          if (fluid_mask_gpu(i,j,k)==0) then
            dv = 1._rkind/jac_gpu(i,j)
            rhou = w_gpu(i,j,k,2)*dcsidxnc2_gpu(i,j)+ w_gpu(i,j,k,3)*dcsidync2_gpu(i,j)
            f_rhou = fln_gpu(i,j,k,2)*dcsidxnc2_gpu(i,j)+fln_gpu(i,j,k,3)*dcsidync2_gpu(i,j)
            bulk_1 = bulk_1 + fln_gpu(i,j,k,1)*dv
            bulk_2 = bulk_2 + w_gpu(i,j,k,1)*dv
            bulk_3 = bulk_3 + rhou*dv
            bulk_4 = bulk_4 + rhou*w_aux_gpu(i,j,k,6)*dv
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

  endsubroutine force_rhs_1_c2_kernel



  subroutine force_var_1_kernel(nx,ny,nz,ng,yn_gpu,fln_gpu,w_gpu,w_aux_gpu,bulkt,fluid_mask_gpu,&
  &cv_coeff_gpu,indx_cp_l,indx_cp_r,t0,calorically_perfect,tol_iter_nr)
    integer, intent(in) :: nx, ny, nz, ng, indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind), intent(in) :: t0, tol_iter_nr
    real(rkind), dimension(1:), intent(in) :: yn_gpu
    real(rkind), dimension(1:,1:,1:,1:), intent(in) :: fln_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in) :: w_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(inout) :: w_aux_gpu
    integer, dimension(1-ng:,1-ng:,1-ng:), intent(in) :: fluid_mask_gpu
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in) :: cv_coeff_gpu
    real(rkind), intent(out) :: bulkt
    real(rkind) :: bulk_5,rho,rhou,rhov,rhow,rhoe,ri,uu,vv,ww,qq,ee,tt
    real(rkind) :: dy
    integer :: i,j,k,iercuda
    bulk_5 = 0._rkind
    !$omp target teams distribute parallel do collapse(2) has_device_addr(yn_gpu,fln_gpu,w_gpu,w_aux_gpu,cv_coeff_gpu,fluid_mask_gpu) &
    !$omp& reduction(+:bulk_5)
    do j = 1,ny
      do i = 1,nx
        do k = 1,nz
          if (fluid_mask_gpu(i,j,k)==0) then
            dy = yn_gpu(j+1)-yn_gpu(j)
            rho = w_gpu(i,j,k,1)
            rhou = w_gpu(i,j,k,2)
            rhov = w_gpu(i,j,k,3)
            rhow = w_gpu(i,j,k,4)
            rhoe = w_gpu(i,j,k,5)
            ri = 1._rkind/rho
            uu = rhou*ri
            vv = rhov*ri
            ww = rhow*ri
            qq = 0.5_rkind*(uu*uu+vv*vv+ww*ww)
            ee = rhoe/rho-qq
            tt = get_temperature_from_e_dev(ee,w_aux_gpu(i,j,k,6),t0,cv_coeff_gpu,indx_cp_l,&
            &indx_cp_r, calorically_perfect,tol_iter_nr)
            w_aux_gpu(i,j,k,6) = tt
            bulk_5 = bulk_5 + rhou*tt*dy
          endif
        enddo
      enddo
    enddo
    bulkt = bulk_5

  endsubroutine force_var_1_kernel


  subroutine force_var_2_kernel(nx,ny,nz,ng,w_gpu,w_aux_gpu,tbdiff,fluid_mask_gpu,cv_coeff_gpu,&
  &indx_cp_l,indx_cp_r,t0,calorically_perfect)
    integer :: nx, ny, nz, ng, indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind), intent(in) :: tbdiff, t0
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(inout) :: w_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in) :: w_aux_gpu
    integer, dimension(1-ng:,1-ng:,1-ng:), intent(in) :: fluid_mask_gpu
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in) :: cv_coeff_gpu
    real(rkind) :: rho,rhou,rhov,rhow,tt,ttnew,ee
    integer :: i,j,k,iercuda

    !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu,w_aux_gpu,cv_coeff_gpu,fluid_mask_gpu)
    do j = 1,ny
      do k = 1,nz
        do i = 1,nx
          if (fluid_mask_gpu(i,j,k)==0) then
            rho = w_gpu(i,j,k,1)
            rhou = w_gpu(i,j,k,2)
            rhov = w_gpu(i,j,k,3)
            rhow = w_gpu(i,j,k,4)
            tt = w_aux_gpu(i,j,k,6)
            ttnew = tt+tbdiff
            ee = get_e_from_temperature_dev(ttnew, t0, indx_cp_l, indx_cp_r, cv_coeff_gpu, calorically_perfect)
            w_gpu(i,j,k,5) = rho*ee + 0.5_rkind*(rhou**2+rhov**2+rhow**2)/rho
          endif
        enddo
      enddo
    enddo


  endsubroutine force_var_2_kernel



  subroutine force_var_1_c2_kernel(nx,ny,nz,ng,yn_gpu,fln_gpu,w_gpu,w_aux_gpu,bulkt,fluid_mask_gpu,&
  &cv_coeff_gpu,indx_cp_l,indx_cp_r,t0,calorically_perfect,tol_iter_nr,jac_gpu,dcsidxnc2_gpu,&
  &dcsidync2_gpu)
    integer, intent(in) :: nx, ny, nz, ng, indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind), intent(in) :: t0, tol_iter_nr
    real(rkind), dimension(1:), intent(in) :: yn_gpu
    real(rkind), dimension(1:,1:,1:,1:), intent(in) :: fln_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in) :: w_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(inout) :: w_aux_gpu
    integer, dimension(1-ng:,1-ng:,1-ng:), intent(in) :: fluid_mask_gpu
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in) :: cv_coeff_gpu
    real(rkind), dimension(1-ng:,1-ng:), intent(in) :: jac_gpu, dcsidxnc2_gpu, dcsidync2_gpu
    real(rkind), intent(out) :: bulkt
    real(rkind) :: bulk_5,rho,rhou,rhov,rhow,rhoe,ri,uu,vv,ww,qq,ee,tt
    real(rkind) :: dy, dv, rhout
    integer :: i,j,k,iercuda
    bulk_5 = 0._rkind
    !$omp target teams distribute parallel do collapse(2) has_device_addr(yn_gpu,fln_gpu,w_gpu,w_aux_gpu,cv_coeff_gpu,jac_gpu,dcsidxnc2_gpu,dcsidync2_gpu,fluid_mask_gpu) &
    !$omp& reduction(+:bulk_5)
    do j = 1,ny
      do i = 1,nx
        do k = 1,nz
          if (fluid_mask_gpu(i,j,k)==0) then
            dv = 1._rkind/jac_gpu(i,j)
            rho = w_gpu(i,j,k,1)
            rhou = w_gpu(i,j,k,2)
            rhov = w_gpu(i,j,k,3)
            rhow = w_gpu(i,j,k,4)
            rhoe = w_gpu(i,j,k,5)
            rhout = rhou*dcsidxnc2_gpu(i,j)+rhov*dcsidync2_gpu(i,j)
            ri = 1._rkind/rho
            uu = rhou*ri
            vv = rhov*ri
            ww = rhow*ri
            qq = 0.5_rkind*(uu*uu+vv*vv+ww*ww)
            ee = rhoe/rho-qq
            tt = get_temperature_from_e_dev(ee,w_aux_gpu(i,j,k,6),t0,cv_coeff_gpu,indx_cp_l,&
            &indx_cp_r, calorically_perfect,tol_iter_nr)
            w_aux_gpu(i,j,k,6) = tt
            bulk_5 = bulk_5 + rhou*tt*dv
          endif
        enddo
      enddo
    enddo
    bulkt = bulk_5

  endsubroutine force_var_1_c2_kernel


  subroutine force_var_2_c2_kernel(nx,ny,nz,ng,w_gpu,w_aux_gpu,tbdiff,fluid_mask_gpu,cv_coeff_gpu,&
  &indx_cp_l,indx_cp_r,t0,calorically_perfect,jac_gpu)
    integer :: nx, ny, nz, ng, indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind), intent(in) :: tbdiff, t0
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(inout) :: w_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in) :: w_aux_gpu
    integer, dimension(1-ng:,1-ng:,1-ng:), intent(in) :: fluid_mask_gpu
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in) :: cv_coeff_gpu
    real(rkind), dimension(1-ng:,1-ng:), intent(in) :: jac_gpu
    real(rkind) :: rho,rhou,rhov,rhow,tt,ttnew,ee,dv
    integer :: i,j,k,iercuda

    !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu,w_aux_gpu,cv_coeff_gpu,jac_gpu,fluid_mask_gpu)
    do j = 1,ny
      do k = 1,nz
        do i = 1,nx
          if (fluid_mask_gpu(i,j,k)==0) then
            dv = 1._rkind/jac_gpu(i,j)
            rho = w_gpu(i,j,k,1)
            rhou = w_gpu(i,j,k,2)
            rhov = w_gpu(i,j,k,3)
            rhow = w_gpu(i,j,k,4)
            tt = w_aux_gpu(i,j,k,6)
            ttnew = tt+tbdiff
            ee = get_e_from_temperature_dev(ttnew, t0, indx_cp_l, indx_cp_r, cv_coeff_gpu, calorically_perfect)
            w_gpu(i,j,k,5) = rho*ee + 0.5_rkind*(rhou**2+rhov**2+rhow**2)/rho
          endif
        enddo
      enddo
    enddo


  endsubroutine force_var_2_c2_kernel



  subroutine update_flux_kernel(nx,ny,nz,nv,fl_gpu,fln_gpu,gamdt)
    integer :: nx, ny, nz, nv
    real(rkind) :: gamdt
    real(rkind), dimension(1:,1:,1:,1:), intent(in) :: fl_gpu
    real(rkind), dimension(1:,1:,1:,1:), intent(inout) :: fln_gpu
    integer :: i,j,k,m,iercuda

    !$omp target teams distribute parallel do collapse(3) has_device_addr(fl_gpu,fln_gpu)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          do m=1,nv
            fln_gpu(i,j,k,m) = fln_gpu(i,j,k,m)-gamdt*fl_gpu(i,j,k,m)
          enddo
        enddo
      enddo
    enddo


  endsubroutine update_flux_kernel



  subroutine update_field_kernel(nx,ny,nz,ng,nv,w_gpu,fln_gpu,fluid_mask_gpu)
    integer :: nx, ny, nz, nv, ng
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(inout) :: w_gpu
    real(rkind), dimension(1:,1:,1:,1:), intent(in) :: fln_gpu
    integer, dimension(1-ng:,1-ng:,1-ng:), intent(in) :: fluid_mask_gpu
    integer :: i,j,k,m,iercuda

    !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu,fln_gpu,fluid_mask_gpu)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          if (fluid_mask_gpu(i,j,k)==0) then
            do m=1,nv
              w_gpu(i,j,k,m) = w_gpu(i,j,k,m)+fln_gpu(i,j,k,m)
            enddo
          endif
        enddo
      enddo
    enddo


  endsubroutine update_field_kernel



  subroutine visflx_div_ord2_kernel(nx,ny,nz,ng,w_aux_gpu,fl_gpu,x_gpu,y_gpu,z_gpu,stream_id)
    integer, intent(in) :: nx, ny, nz, ng
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(inout) :: w_aux_gpu
    real(rkind), dimension(1:,1:,1:, 1:), intent(inout) :: fl_gpu
    real(rkind), dimension(1-ng:), intent(in) :: x_gpu, y_gpu, z_gpu
    integer :: i,j,k,iercuda
    real(rkind) :: dxl,dyl,dzl
    real(rkind) :: uu,vv,ww,mu
    real(rkind) :: sigq,sigx,sigy,sigz
    integer, intent(in) :: stream_id

    !$omp target teams distribute parallel do collapse(3) has_device_addr(w_aux_gpu,fl_gpu,x_gpu,y_gpu,z_gpu)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          uu = w_aux_gpu(i,j,k,2)
          vv = w_aux_gpu(i,j,k,3)
          ww = w_aux_gpu(i,j,k,4)
          mu = w_aux_gpu(i,j,k,7)
          dxl = mu/(x_gpu(i+1)-x_gpu(i-1))
          dyl = mu/(y_gpu(j+1)-y_gpu(j-1))
          dzl = mu/(z_gpu(k+1)-z_gpu(k-1))
          sigx = dxl*(w_aux_gpu(i+1,j,k,10)-w_aux_gpu(i-1,j,k,10))
          sigy = dyl*(w_aux_gpu(i,j+1,k,10)-w_aux_gpu(i,j-1,k,10))
          sigz = dzl*(w_aux_gpu(i,j,k+1,10)-w_aux_gpu(i,j,k-1,10))
          sigq = sigx*uu+sigy*vv+sigz*ww
          fl_gpu(i,j,k,2) = fl_gpu(i,j,k,2) - sigx
          fl_gpu(i,j,k,3) = fl_gpu(i,j,k,3) - sigy
          fl_gpu(i,j,k,4) = fl_gpu(i,j,k,4) - sigz
          fl_gpu(i,j,k,5) = fl_gpu(i,j,k,5) - sigq
        enddo
      enddo
    enddo


  endsubroutine visflx_div_ord2_kernel



  subroutine visflx_div_kernel(nx,ny,nz,ng,visc_order,w_aux_gpu,fl_gpu,coeff_deriv1_gpu,dcsidx_gpu,detady_gpu,dzitdz_gpu,stream_id)
    integer, intent(in) :: nx, ny, nz, ng, visc_order
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(inout) :: w_aux_gpu
    real(rkind), dimension(1:,1:,1:, 1:), intent(inout) :: fl_gpu
    real(rkind), dimension(1:,1:), intent(in) :: coeff_deriv1_gpu
    real(rkind), dimension(1:), intent(in) :: dcsidx_gpu, detady_gpu, dzitdz_gpu
    integer :: i,j,k,l,iercuda
    real(rkind) :: ccl
    real(rkind) :: uu,vv,ww,tt,mu
    real(rkind) :: sigq,sigx,sigy,sigz
    real(rkind) :: divx3l, divy3l, divz3l
    integer :: lmax
    integer, intent(in) :: stream_id
    lmax = visc_order/2
    !$omp target teams distribute parallel do collapse(2) has_device_addr(w_aux_gpu,fl_gpu,coeff_deriv1_gpu,dcsidx_gpu,detady_gpu,dzitdz_gpu)
    do k = 1,nz
      do i = 1,nx
        do j = 1,ny
          uu = w_aux_gpu(i,j,k,2)
          vv = w_aux_gpu(i,j,k,3)
          ww = w_aux_gpu(i,j,k,4)
          tt = w_aux_gpu(i,j,k,6)
          mu = w_aux_gpu(i,j,k,7)
          divx3l = 0._rkind
          divy3l = 0._rkind
          divz3l = 0._rkind
          do l=1,lmax
            ccl = coeff_deriv1_gpu(l,lmax)
            divx3l = divx3l+ccl*(w_aux_gpu(i+l,j,k,10)-w_aux_gpu(i-l,j,k,10))
            divy3l = divy3l+ccl*(w_aux_gpu(i,j+l,k,10)-w_aux_gpu(i,j-l,k,10))
            divz3l = divz3l+ccl*(w_aux_gpu(i,j,k+l,10)-w_aux_gpu(i,j,k-l,10))
          enddo
          divx3l = divx3l*dcsidx_gpu(i)
          divy3l = divy3l*detady_gpu(j)
          divz3l = divz3l*dzitdz_gpu(k)
          sigx = mu*divx3l
          sigy = mu*divy3l
          sigz = mu*divz3l
          sigq = sigx*uu+sigy*vv+sigz*ww
          fl_gpu(i,j,k,2) = fl_gpu(i,j,k,2) - sigx
          fl_gpu(i,j,k,3) = fl_gpu(i,j,k,3) - sigy
          fl_gpu(i,j,k,4) = fl_gpu(i,j,k,4) - sigz
          fl_gpu(i,j,k,5) = fl_gpu(i,j,k,5) - sigq
        enddo
      enddo
    enddo


  endsubroutine visflx_div_kernel



  subroutine visflx_div_c2_kernel(nx,ny,nz,ng,visc_order,w_aux_gpu,fl_gpu,coeff_deriv1_gpu,&
  &dcsidxc2_gpu,dcsidyc2_gpu,detadxc2_gpu,detadyc2_gpu,dzitdz_gpu,vis_tag_gpu,wall_tag_gpu,stream_id)
    integer, intent(in) :: nx, ny, nz, ng, visc_order
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(inout) :: w_aux_gpu
    real(rkind), dimension(1:,1:,1:, 1:), intent(inout) :: fl_gpu
    real(rkind), dimension(1:,1:), intent(in) :: coeff_deriv1_gpu
    real(rkind), dimension(1-ng:,1-ng:), intent(in) :: dcsidxc2_gpu, detadyc2_gpu
    real(rkind), dimension(1-ng:,1-ng:), intent(in) :: dcsidyc2_gpu, detadxc2_gpu
    real(rkind), dimension(1:), intent(in) :: dzitdz_gpu
    integer, dimension(1-ng:), intent(in) :: vis_tag_gpu, wall_tag_gpu
    integer :: i,j,k,l,iercuda
    real(rkind) :: cli, clj, clk
    real(rkind) :: uu,vv,ww,tt,mu
    real(rkind) :: sigq,sigx,sigy,sigz
    real(rkind) :: divcsi, diveta, divzit
    real(rkind) :: divx3l, divy3l, divz3l
    integer :: lmax,lmaxi,lmaxj
    integer, intent(in) :: stream_id
    lmax = visc_order/2
    !$omp target teams distribute parallel do collapse(3) has_device_addr(w_aux_gpu,fl_gpu,coeff_deriv1_gpu,dcsidxc2_gpu,detadyc2_gpu,dcsidyc2_gpu,detadxc2_gpu,dzitdz_gpu,vis_tag_gpu,wall_tag_gpu)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          uu = w_aux_gpu(i,j,k,2)
          vv = w_aux_gpu(i,j,k,3)
          ww = w_aux_gpu(i,j,k,4)
          tt = w_aux_gpu(i,j,k,6)
          mu = w_aux_gpu(i,j,k,7)
          lmaxi = lmax
          lmaxj = lmax
          if (j==1) lmaxi = vis_tag_gpu(i)
          if (wall_tag_gpu(i) < 1) lmaxj = min(j,lmax)
          divcsi = 0._rkind
          diveta = 0._rkind
          divzit = 0._rkind
          do l=1,lmax
            cli = coeff_deriv1_gpu(l,lmaxi)
            clj = coeff_deriv1_gpu(l,lmaxj)
            clk = coeff_deriv1_gpu(l,lmax )
            divcsi = divcsi+cli*(w_aux_gpu(i+l,j,k,10)-w_aux_gpu(i-l,j,k,10))
            diveta = diveta+clj*(w_aux_gpu(i,j+l,k,10)-w_aux_gpu(i,j-l,k,10))
            divzit = divzit+clk*(w_aux_gpu(i,j,k+l,10)-w_aux_gpu(i,j,k-l,10))
          enddo
          divx3l = divcsi*dcsidxc2_gpu(i,j) + diveta*detadxc2_gpu(i,j)
          divy3l = divcsi*dcsidyc2_gpu(i,j) + diveta*detadyc2_gpu(i,j)
          divz3l = divzit*dzitdz_gpu(k)
          sigx = mu*divx3l
          sigy = mu*divy3l
          sigz = mu*divz3l
          sigq = sigx*uu+sigy*vv+sigz*ww
          fl_gpu(i,j,k,2) = fl_gpu(i,j,k,2) - sigx
          fl_gpu(i,j,k,3) = fl_gpu(i,j,k,3) - sigy
          fl_gpu(i,j,k,4) = fl_gpu(i,j,k,4) - sigz
          fl_gpu(i,j,k,5) = fl_gpu(i,j,k,5) - sigq
        enddo
      enddo
    enddo


  endsubroutine visflx_div_c2_kernel



  subroutine visflx_kernel(nx,ny,nz,ng,visc_order,prandtl,t0,indx_cp_l,indx_cp_r,cp_coeff_gpu,&
  &calorically_perfect,u0,l0,w_gpu,w_aux_gpu,fl_gpu,coeff_deriv1_gpu,coeff_deriv2_gpu,dcsidx_gpu,&
  &detady_gpu,dzitdz_gpu,dcsidxs_gpu,detadys_gpu,dzitdzs_gpu,dcsidx2_gpu,detady2_gpu,dzitdz2_gpu,&
  &wallprop_gpu)
    integer, intent(in) :: nx, ny, nz, ng, visc_order, calorically_perfect
    integer, intent(in) :: indx_cp_l, indx_cp_r
    real(rkind), intent(in) :: prandtl, u0, l0, t0
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(in) :: w_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(inout) :: w_aux_gpu
    real(rkind), dimension(1-ng:,1-ng:, 2:), intent(inout) :: wallprop_gpu
    real(rkind), dimension(1:,1:,1:, 1:), intent(inout) :: fl_gpu
    real(rkind), dimension(1:,1:), intent(in) :: coeff_deriv1_gpu
    real(rkind), dimension(0:,1:), intent(in) :: coeff_deriv2_gpu
    real(rkind), dimension(1:), intent(in) :: dcsidx_gpu, detady_gpu, dzitdz_gpu
    real(rkind), dimension(1:), intent(in) :: dcsidxs_gpu, detadys_gpu, dzitdzs_gpu
    real(rkind), dimension(1:), intent(in) :: dcsidx2_gpu, detady2_gpu, dzitdz2_gpu
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in) :: cp_coeff_gpu
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
    !$omp target teams distribute parallel do collapse(2) has_device_addr(w_gpu,w_aux_gpu,wallprop_gpu,fl_gpu,coeff_deriv1_gpu,coeff_deriv2_gpu,dcsidx_gpu,detady_gpu,dzitdz_gpu,dcsidxs_gpu,detadys_gpu,dzitdzs_gpu,dcsidx2_gpu,detady2_gpu,dzitdz2_gpu,cp_coeff_gpu)
    do j = 1,ny
      do i = 1,nx
        do k = 1,nz
          uu = w_aux_gpu(i,j,k,2)
          vv = w_aux_gpu(i,j,k,3)
          ww = w_aux_gpu(i,j,k,4)
          tt = w_aux_gpu(i,j,k,6)
          mu = w_aux_gpu(i,j,k,7)
          if (calorically_perfect==1) then
            cploc = cp_coeff_gpu(0)
          else
            cploc = 0._rkind
            do ll=indx_cp_l,indx_cp_r
              cploc = cploc+cp_coeff_gpu(ll)*(tt/t0)**ll
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
            ccl = coeff_deriv1_gpu(l,lmax)
            ux = ux+ccl*(w_aux_gpu(i+l,j,k,2)-w_aux_gpu(i-l,j,k,2))
            vx = vx+ccl*(w_aux_gpu(i+l,j,k,3)-w_aux_gpu(i-l,j,k,3))
            wx = wx+ccl*(w_aux_gpu(i+l,j,k,4)-w_aux_gpu(i-l,j,k,4))
            tx = tx+ccl*(w_aux_gpu(i+l,j,k,6)-w_aux_gpu(i-l,j,k,6))
            mux = mux+ccl*(w_aux_gpu(i+l,j,k,7)-w_aux_gpu(i-l,j,k,7))
            uy = uy+ccl*(w_aux_gpu(i,j+l,k,2)-w_aux_gpu(i,j-l,k,2))
            vy = vy+ccl*(w_aux_gpu(i,j+l,k,3)-w_aux_gpu(i,j-l,k,3))
            wy = wy+ccl*(w_aux_gpu(i,j+l,k,4)-w_aux_gpu(i,j-l,k,4))
            ty = ty+ccl*(w_aux_gpu(i,j+l,k,6)-w_aux_gpu(i,j-l,k,6))
            muy = muy+ccl*(w_aux_gpu(i,j+l,k,7)-w_aux_gpu(i,j-l,k,7))
            uz = uz+ccl*(w_aux_gpu(i,j,k+l,2)-w_aux_gpu(i,j,k-l,2))
            vz = vz+ccl*(w_aux_gpu(i,j,k+l,3)-w_aux_gpu(i,j,k-l,3))
            wz = wz+ccl*(w_aux_gpu(i,j,k+l,4)-w_aux_gpu(i,j,k-l,4))
            tz = tz+ccl*(w_aux_gpu(i,j,k+l,6)-w_aux_gpu(i,j,k-l,6))
            muz = muz+ccl*(w_aux_gpu(i,j,k+l,7)-w_aux_gpu(i,j,k-l,7))
          enddo
          ux = ux*dcsidx_gpu(i)
          vx = vx*dcsidx_gpu(i)
          wx = wx*dcsidx_gpu(i)
          tx = tx*dcsidx_gpu(i)
          mux = mux*dcsidx_gpu(i)
          uy = uy*detady_gpu(j)
          vy = vy*detady_gpu(j)
          wy = wy*detady_gpu(j)
          ty = ty*detady_gpu(j)
          muy = muy*detady_gpu(j)
          uz = uz*dzitdz_gpu(k)
          vz = vz*dzitdz_gpu(k)
          wz = wz*dzitdz_gpu(k)
          tz = tz*dzitdz_gpu(k)
          muz = muz*dzitdz_gpu(k)
          if (j==1) then
            wallprop_gpu(i,k,2) = mu*uy
            wallprop_gpu(i,k,3) = mu*wy
            wallprop_gpu(i,k,4) = mu*ty*cploc/prandtl
          endif
          ulapx = coeff_deriv2_gpu(0,lmax)*uu
          ulapy = ulapx
          ulapz = ulapx
          vlapx = coeff_deriv2_gpu(0,lmax)*vv
          vlapy = vlapx
          vlapz = vlapx
          wlapx = coeff_deriv2_gpu(0,lmax)*ww
          wlapy = wlapx
          wlapz = wlapx
          tlapx = coeff_deriv2_gpu(0,lmax)*tt
          tlapy = tlapx
          tlapz = tlapx
          do l=1,lmax
            clapl = coeff_deriv2_gpu(l,lmax)
            ulapx = ulapx + clapl*(w_aux_gpu(i+l,j,k,2)+w_aux_gpu(i-l,j,k,2))
            ulapy = ulapy + clapl*(w_aux_gpu(i,j+l,k,2)+w_aux_gpu(i,j-l,k,2))
            ulapz = ulapz + clapl*(w_aux_gpu(i,j,k+l,2)+w_aux_gpu(i,j,k-l,2))
            vlapx = vlapx + clapl*(w_aux_gpu(i+l,j,k,3)+w_aux_gpu(i-l,j,k,3))
            vlapy = vlapy + clapl*(w_aux_gpu(i,j+l,k,3)+w_aux_gpu(i,j-l,k,3))
            vlapz = vlapz + clapl*(w_aux_gpu(i,j,k+l,3)+w_aux_gpu(i,j,k-l,3))
            wlapx = wlapx + clapl*(w_aux_gpu(i+l,j,k,4)+w_aux_gpu(i-l,j,k,4))
            wlapy = wlapy + clapl*(w_aux_gpu(i,j+l,k,4)+w_aux_gpu(i,j-l,k,4))
            wlapz = wlapz + clapl*(w_aux_gpu(i,j,k+l,4)+w_aux_gpu(i,j,k-l,4))
            tlapx = tlapx + clapl*(w_aux_gpu(i+l,j,k,6)+w_aux_gpu(i-l,j,k,6))
            tlapy = tlapy + clapl*(w_aux_gpu(i,j+l,k,6)+w_aux_gpu(i,j-l,k,6))
            tlapz = tlapz + clapl*(w_aux_gpu(i,j,k+l,6)+w_aux_gpu(i,j,k-l,6))
          enddo
          ulapx = ulapx*dcsidxs_gpu(i)+ux*dcsidx2_gpu(i)
          vlapx = vlapx*dcsidxs_gpu(i)+vx*dcsidx2_gpu(i)
          wlapx = wlapx*dcsidxs_gpu(i)+wx*dcsidx2_gpu(i)
          tlapx = tlapx*dcsidxs_gpu(i)+tx*dcsidx2_gpu(i)
          ulapy = ulapy*detadys_gpu(j)+uy*detady2_gpu(j)
          vlapy = vlapy*detadys_gpu(j)+vy*detady2_gpu(j)
          wlapy = wlapy*detadys_gpu(j)+wy*detady2_gpu(j)
          tlapy = tlapy*detadys_gpu(j)+ty*detady2_gpu(j)
          ulapz = ulapz*dzitdzs_gpu(k)+uz*dzitdz2_gpu(k)
          vlapz = vlapz*dzitdzs_gpu(k)+vz*dzitdz2_gpu(k)
          wlapz = wlapz*dzitdzs_gpu(k)+wz*dzitdz2_gpu(k)
          tlapz = tlapz*dzitdzs_gpu(k)+tz*dzitdz2_gpu(k)
          ulap = ulapx+ulapy+ulapz
          vlap = vlapx+vlapy+vlapz
          wlap = wlapx+wlapy+wlapz
          tlap = tlapx+tlapy+tlapz
          div = ux+vy+wz
          div3l = div/3._rkind
          w_aux_gpu(i,j,k,10) = div3l
          omegax = wy-vz
          omegay = uz-wx
          omegaz = vx-uy
          omod2 = omegax*omegax+omegay*omegay+omegaz*omegaz
          w_aux_gpu(i,j,k,9) = sqrt(omod2)
          w_aux_gpu(i,j,k,8) = (max(-div/sqrt(omod2+div**2+(u0/l0)**2),0._rkind))**2
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
          fl_gpu(i,j,k,2) = fl_gpu(i,j,k,2) - sigx
          fl_gpu(i,j,k,3) = fl_gpu(i,j,k,3) - sigy
          fl_gpu(i,j,k,4) = fl_gpu(i,j,k,4) - sigz
          fl_gpu(i,j,k,5) = fl_gpu(i,j,k,5) - sigq
        enddo
      enddo
    enddo


  endsubroutine visflx_kernel



  subroutine visflx_c2_kernel(nx,ny,nz,ng,visc_order,prandtl,t0,indx_cp_l,indx_cp_r,cp_coeff_gpu,&
  &calorically_perfect,u0,l0,w_gpu,w_aux_gpu,fl_gpu,coeff_deriv1_gpu,coeff_deriv2_gpu,dcsidxc2_gpu,&
  &detadyc2_gpu,detadxc2_gpu,dcsidyc2_gpu,dzitdz_gpu,dzitdzs_gpu,dzitdz2_gpu,g1_gpu,g2_gpu,g12_gpu,&
  &jac_gpu,wallprop_gpu,vis_tag_gpu,wall_tag_gpu,iblock,ite_rank_x,itu_rank_x,ite_l,itu_l,teshk)
    integer, intent(in) :: nx, ny, nz, ng, visc_order, calorically_perfect
    integer, intent(in) :: indx_cp_l, indx_cp_r
    real(rkind), intent(in) :: prandtl, u0, l0, t0
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(in) :: w_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(inout) :: w_aux_gpu
    real(rkind), dimension(1-ng:,1-ng:, 2:), intent(inout) :: wallprop_gpu
    real(rkind), dimension(1:,1:,1:, 1:), intent(inout) :: fl_gpu
    real(rkind), dimension(1:,1:), intent(in) :: coeff_deriv1_gpu
    real(rkind), dimension(0:,1:), intent(in) :: coeff_deriv2_gpu
    real(rkind), dimension(1:), intent(in) :: dzitdz_gpu, dzitdzs_gpu, dzitdz2_gpu
    real(rkind), dimension(1-ng:,1-ng:), intent(in) :: dcsidxc2_gpu, detadyc2_gpu
    real(rkind), dimension(1-ng:,1-ng:), intent(in) :: detadxc2_gpu, dcsidyc2_gpu
    real(rkind), dimension(1-ng:,1-ng:), intent(in) :: g1_gpu, g2_gpu, g12_gpu, jac_gpu
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in) :: cp_coeff_gpu
    integer, dimension(1-ng:), intent(in) :: vis_tag_gpu, wall_tag_gpu
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
    !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu,w_aux_gpu,wallprop_gpu,fl_gpu,coeff_deriv1_gpu,coeff_deriv2_gpu,dzitdz_gpu,dzitdzs_gpu,dzitdz2_gpu,dcsidxc2_gpu,detadyc2_gpu,detadxc2_gpu,dcsidyc2_gpu,g1_gpu,g2_gpu,g12_gpu,jac_gpu,cp_coeff_gpu,vis_tag_gpu,wall_tag_gpu)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          uu = w_aux_gpu(i,j,k,2)
          vv = w_aux_gpu(i,j,k,3)
          ww = w_aux_gpu(i,j,k,4)
          tt = w_aux_gpu(i,j,k,6)
          mu = w_aux_gpu(i,j,k,7)
          lmaxi = lmax
          lmaxj = lmax
          if (j == 1) lmaxi = vis_tag_gpu(i)
          if (wall_tag_gpu(i) < 1) lmaxj = min(j,lmax)
          if (calorically_perfect==1) then
            cploc = cp_coeff_gpu(0)
          else
            cploc = 0._rkind
            do ll=indx_cp_l,indx_cp_r
              cploc = cploc+cp_coeff_gpu(ll)*(tt/t0)**ll
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
            cli = coeff_deriv1_gpu(l,lmaxi)
            clj = coeff_deriv1_gpu(l,lmaxj)
            clk = coeff_deriv1_gpu(l,lmax )
            ucsi = ucsi +cli*(w_aux_gpu(i+l,j,k,2)-w_aux_gpu(i-l,j,k,2))
            vcsi = vcsi +cli*(w_aux_gpu(i+l,j,k,3)-w_aux_gpu(i-l,j,k,3))
            wcsi = wcsi +cli*(w_aux_gpu(i+l,j,k,4)-w_aux_gpu(i-l,j,k,4))
            tcsi = tcsi +cli*(w_aux_gpu(i+l,j,k,6)-w_aux_gpu(i-l,j,k,6))
            mucsi = mucsi+cli*(w_aux_gpu(i+l,j,k,7)-w_aux_gpu(i-l,j,k,7))
            ueta = ueta +clj*(w_aux_gpu(i,j+l,k,2)-w_aux_gpu(i,j-l,k,2))
            veta = veta +clj*(w_aux_gpu(i,j+l,k,3)-w_aux_gpu(i,j-l,k,3))
            weta = weta +clj*(w_aux_gpu(i,j+l,k,4)-w_aux_gpu(i,j-l,k,4))
            teta = teta +clj*(w_aux_gpu(i,j+l,k,6)-w_aux_gpu(i,j-l,k,6))
            mueta = mueta+clj*(w_aux_gpu(i,j+l,k,7)-w_aux_gpu(i,j-l,k,7))
            uzit = uzit +clk*(w_aux_gpu(i,j,k+l,2)-w_aux_gpu(i,j,k-l,2))
            vzit = vzit +clk*(w_aux_gpu(i,j,k+l,3)-w_aux_gpu(i,j,k-l,3))
            wzit = wzit +clk*(w_aux_gpu(i,j,k+l,4)-w_aux_gpu(i,j,k-l,4))
            tzit = tzit +clk*(w_aux_gpu(i,j,k+l,6)-w_aux_gpu(i,j,k-l,6))
            muzit = muzit+clk*(w_aux_gpu(i,j,k+l,7)-w_aux_gpu(i,j,k-l,7))
            dg1 = dg1 +cli*(g1_gpu (i+l,j)-g1_gpu (i-l,j))
            dg12csi=dg12csi+cli*(g12_gpu(i+l,j)-g12_gpu(i-l,j))
            dg2 = dg2 +clj*(g2_gpu (i,j+l)-g2_gpu (i,j-l))
            dg12eta=dg12eta+clj*(g12_gpu(i,j+l)-g12_gpu(i,j-l))
          enddo
          ux = ucsi *dcsidxc2_gpu(i,j) + ueta *detadxc2_gpu(i,j)
          vx = vcsi *dcsidxc2_gpu(i,j) + veta *detadxc2_gpu(i,j)
          wx = wcsi *dcsidxc2_gpu(i,j) + weta *detadxc2_gpu(i,j)
          tx = tcsi *dcsidxc2_gpu(i,j) + teta *detadxc2_gpu(i,j)
          mux = mucsi*dcsidxc2_gpu(i,j) + mueta*detadxc2_gpu(i,j)
          uy = ucsi *dcsidyc2_gpu(i,j) + ueta *detadyc2_gpu(i,j)
          vy = vcsi *dcsidyc2_gpu(i,j) + veta *detadyc2_gpu(i,j)
          wy = wcsi *dcsidyc2_gpu(i,j) + weta *detadyc2_gpu(i,j)
          ty = tcsi *dcsidyc2_gpu(i,j) + teta *detadyc2_gpu(i,j)
          muy = mucsi*dcsidyc2_gpu(i,j) + mueta*detadyc2_gpu(i,j)
          uz = uzit *dzitdz_gpu(k)
          vz = vzit *dzitdz_gpu(k)
          wz = wzit *dzitdz_gpu(k)
          tz = tzit *dzitdz_gpu(k)
          muz = muzit*dzitdz_gpu(k)
          if (j==1) then
            wallprop_gpu(i,k,2) = mu*ueta
            wallprop_gpu(i,k,3) = mu*weta
            wallprop_gpu(i,k,4) = mu*teta*cploc/prandtl
          endif
          ulapcsi = coeff_deriv2_gpu(0,lmaxi)*uu
          ulapeta = coeff_deriv2_gpu(0,lmaxj)*uu
          ulapzit = coeff_deriv2_gpu(0,lmax)*uu
          vlapcsi = coeff_deriv2_gpu(0,lmaxi)*vv
          vlapeta = coeff_deriv2_gpu(0,lmaxj)*vv
          vlapzit = coeff_deriv2_gpu(0,lmax)*vv
          wlapcsi = coeff_deriv2_gpu(0,lmaxi)*ww
          wlapeta = coeff_deriv2_gpu(0,lmaxj)*ww
          wlapzit = coeff_deriv2_gpu(0,lmax)*ww
          tlapcsi = coeff_deriv2_gpu(0,lmaxi)*tt
          tlapeta = coeff_deriv2_gpu(0,lmaxj)*tt
          tlapzit = coeff_deriv2_gpu(0,lmax)*tt
          do l=1,lmax
            clapi = coeff_deriv2_gpu(l,lmaxi)
            clapj = coeff_deriv2_gpu(l,lmaxj)
            clapk = coeff_deriv2_gpu(l,lmax )
            ulapcsi = ulapcsi + clapi*(w_aux_gpu(i+l,j,k,2)+w_aux_gpu(i-l,j,k,2))
            ulapeta = ulapeta + clapj*(w_aux_gpu(i,j+l,k,2)+w_aux_gpu(i,j-l,k,2))
            ulapzit = ulapzit + clapk*(w_aux_gpu(i,j,k+l,2)+w_aux_gpu(i,j,k-l,2))
            vlapcsi = vlapcsi + clapi*(w_aux_gpu(i+l,j,k,3)+w_aux_gpu(i-l,j,k,3))
            vlapeta = vlapeta + clapj*(w_aux_gpu(i,j+l,k,3)+w_aux_gpu(i,j-l,k,3))
            vlapzit = vlapzit + clapk*(w_aux_gpu(i,j,k+l,3)+w_aux_gpu(i,j,k-l,3))
            wlapcsi = wlapcsi + clapi*(w_aux_gpu(i+l,j,k,4)+w_aux_gpu(i-l,j,k,4))
            wlapeta = wlapeta + clapj*(w_aux_gpu(i,j+l,k,4)+w_aux_gpu(i,j-l,k,4))
            wlapzit = wlapzit + clapk*(w_aux_gpu(i,j,k+l,4)+w_aux_gpu(i,j,k-l,4))
            tlapcsi = tlapcsi + clapi*(w_aux_gpu(i+l,j,k,6)+w_aux_gpu(i-l,j,k,6))
            tlapeta = tlapeta + clapj*(w_aux_gpu(i,j+l,k,6)+w_aux_gpu(i,j-l,k,6))
            tlapzit = tlapzit + clapk*(w_aux_gpu(i,j,k+l,6)+w_aux_gpu(i,j,k-l,6))
          enddo
          ulap = g1_gpu(i,j)*ulapcsi+ucsi*dg1+g2_gpu(i,j)*ulapeta+ueta*dg2
          vlap = g1_gpu(i,j)*vlapcsi+vcsi*dg1+g2_gpu(i,j)*vlapeta+veta*dg2
          wlap = g1_gpu(i,j)*wlapcsi+wcsi*dg1+g2_gpu(i,j)*wlapeta+weta*dg2
          tlap = g1_gpu(i,j)*tlapcsi+tcsi*dg1+g2_gpu(i,j)*tlapeta+teta*dg2
          ulap = ulap * jac_gpu(i,j)
          vlap = vlap * jac_gpu(i,j)
          wlap = wlap * jac_gpu(i,j)
          tlap = tlap * jac_gpu(i,j)
          ulapz = ulapzit*dzitdzs_gpu(k)+uz*dzitdz2_gpu(k)
          vlapz = vlapzit*dzitdzs_gpu(k)+vz*dzitdz2_gpu(k)
          wlapz = wlapzit*dzitdzs_gpu(k)+wz*dzitdz2_gpu(k)
          tlapz = tlapzit*dzitdzs_gpu(k)+tz*dzitdz2_gpu(k)
          ulap = ulap+ulapz
          vlap = vlap+vlapz
          wlap = wlap+wlapz
          tlap = tlap+tlapz
          div = ux+vy+wz
          div3l = div/3._rkind
          w_aux_gpu(i,j,k,10) = div3l
          omegax = wy-vz
          omegay = uz-wx
          omegaz = vx-uy
          omod2 = omegax*omegax+omegay*omegay+omegaz*omegaz
          w_aux_gpu(i,j,k,9) = sqrt(omod2)
          w_aux_gpu(i,j,k,8) = max(-div/sqrt(omod2+div**2+(u0/l0)**2),0._rkind)
          if( j == 1) then
            if(iblock == ite_rank_x .and. i == ite_l) w_aux_gpu(i,j,k,8) = teshk
            if(iblock == itu_rank_x .and. i == itu_l) w_aux_gpu(i,j,k,8) = teshk
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
          fl_gpu(i,j,k,2) = fl_gpu(i,j,k,2) - sigx
          fl_gpu(i,j,k,3) = fl_gpu(i,j,k,3) - sigy
          fl_gpu(i,j,k,4) = fl_gpu(i,j,k,4) - sigz
          fl_gpu(i,j,k,5) = fl_gpu(i,j,k,5) - sigq
        enddo
      enddo
    enddo


  endsubroutine visflx_c2_kernel




  subroutine visflx_nosensor_kernel(nx,ny,nz,ng,visc_order,prandtl,t0,indx_cp_l,indx_cp_r,&
  &cp_coeff_gpu,calorically_perfect,u0,l0,w_gpu,w_aux_gpu,fl_gpu,coeff_deriv1_gpu,coeff_deriv2_gpu,&
  &dcsidx_gpu,detady_gpu,dzitdz_gpu,dcsidxs_gpu,detadys_gpu,dzitdzs_gpu,dcsidx2_gpu,detady2_gpu,&
  &dzitdz2_gpu,wallprop_gpu)
    integer, intent(in) :: nx, ny, nz, ng, visc_order, calorically_perfect
    integer, intent(in) :: indx_cp_l, indx_cp_r
    real(rkind), intent(in) :: prandtl, u0, l0, t0
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(in) :: w_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(inout) :: w_aux_gpu
    real(rkind), dimension(1-ng:,1-ng:, 2:), intent(inout) :: wallprop_gpu
    real(rkind), dimension(1:,1:,1:, 1:), intent(inout) :: fl_gpu
    real(rkind), dimension(1:,1:), intent(in) :: coeff_deriv1_gpu
    real(rkind), dimension(0:,1:), intent(in) :: coeff_deriv2_gpu
    real(rkind), dimension(1:), intent(in) :: dcsidx_gpu, detady_gpu, dzitdz_gpu
    real(rkind), dimension(1:), intent(in) :: dcsidxs_gpu, detadys_gpu, dzitdzs_gpu
    real(rkind), dimension(1:), intent(in) :: dcsidx2_gpu, detady2_gpu, dzitdz2_gpu
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in) :: cp_coeff_gpu
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
    !$omp target teams distribute parallel do collapse(2) has_device_addr(w_gpu,w_aux_gpu,wallprop_gpu,fl_gpu,coeff_deriv1_gpu,coeff_deriv2_gpu,dcsidx_gpu,detady_gpu,dzitdz_gpu,dcsidxs_gpu,detadys_gpu,dzitdzs_gpu,dcsidx2_gpu,detady2_gpu,dzitdz2_gpu,cp_coeff_gpu)
    do j = 1,ny
      do i = 1,nx
        do k = 1,nz
          uu = w_aux_gpu(i,j,k,2)
          vv = w_aux_gpu(i,j,k,3)
          ww = w_aux_gpu(i,j,k,4)
          tt = w_aux_gpu(i,j,k,6)
          mu = w_aux_gpu(i,j,k,7)
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
          ulapx = coeff_deriv2_gpu(0,lmax)*uu
          ulapy = ulapx
          ulapz = ulapx
          vlapx = coeff_deriv2_gpu(0,lmax)*vv
          vlapy = vlapx
          vlapz = vlapx
          wlapx = coeff_deriv2_gpu(0,lmax)*ww
          wlapy = wlapx
          wlapz = wlapx
          tlapx = coeff_deriv2_gpu(0,lmax)*tt
          tlapy = tlapx
          tlapz = tlapx
          do l=1,lmax
            ccl = coeff_deriv1_gpu(l,lmax)
            clapl = coeff_deriv2_gpu(l,lmax)
            ux = ux+ccl*(w_aux_gpu(i+l,j,k,2)-w_aux_gpu(i-l,j,k,2))
            ulapx = ulapx + clapl*(w_aux_gpu(i+l,j,k,2)+w_aux_gpu(i-l,j,k,2))
            vx = vx+ccl*(w_aux_gpu(i+l,j,k,3)-w_aux_gpu(i-l,j,k,3))
            vlapx = vlapx + clapl*(w_aux_gpu(i+l,j,k,3)+w_aux_gpu(i-l,j,k,3))
            wx = wx+ccl*(w_aux_gpu(i+l,j,k,4)-w_aux_gpu(i-l,j,k,4))
            wlapx = wlapx + clapl*(w_aux_gpu(i+l,j,k,4)+w_aux_gpu(i-l,j,k,4))
            tx = tx+ccl*(w_aux_gpu(i+l,j,k,6)-w_aux_gpu(i-l,j,k,6))
            tlapx = tlapx + clapl*(w_aux_gpu(i+l,j,k,6)+w_aux_gpu(i-l,j,k,6))
            mux = mux+ccl*(w_aux_gpu(i+l,j,k,7)-w_aux_gpu(i-l,j,k,7))
          enddo
          do l=1,lmax
            ccl = coeff_deriv1_gpu(l,lmax)
            clapl = coeff_deriv2_gpu(l,lmax)
            uy = uy+ccl*(w_aux_gpu(i,j+l,k,2)-w_aux_gpu(i,j-l,k,2))
            ulapy = ulapy + clapl*(w_aux_gpu(i,j+l,k,2)+w_aux_gpu(i,j-l,k,2))
            vy = vy+ccl*(w_aux_gpu(i,j+l,k,3)-w_aux_gpu(i,j-l,k,3))
            vlapy = vlapy + clapl*(w_aux_gpu(i,j+l,k,3)+w_aux_gpu(i,j-l,k,3))
            wy = wy+ccl*(w_aux_gpu(i,j+l,k,4)-w_aux_gpu(i,j-l,k,4))
            wlapy = wlapy + clapl*(w_aux_gpu(i,j+l,k,4)+w_aux_gpu(i,j-l,k,4))
            ty = ty+ccl*(w_aux_gpu(i,j+l,k,6)-w_aux_gpu(i,j-l,k,6))
            tlapy = tlapy + clapl*(w_aux_gpu(i,j+l,k,6)+w_aux_gpu(i,j-l,k,6))
            muy = muy+ccl*(w_aux_gpu(i,j+l,k,7)-w_aux_gpu(i,j-l,k,7))
          enddo
          do l=1,lmax
            ccl = coeff_deriv1_gpu(l,lmax)
            clapl = coeff_deriv2_gpu(l,lmax)
            uz = uz+ccl*(w_aux_gpu(i,j,k+l,2)-w_aux_gpu(i,j,k-l,2))
            ulapz = ulapz + clapl*(w_aux_gpu(i,j,k+l,2)+w_aux_gpu(i,j,k-l,2))
            vz = vz+ccl*(w_aux_gpu(i,j,k+l,3)-w_aux_gpu(i,j,k-l,3))
            vlapz = vlapz + clapl*(w_aux_gpu(i,j,k+l,3)+w_aux_gpu(i,j,k-l,3))
            wz = wz+ccl*(w_aux_gpu(i,j,k+l,4)-w_aux_gpu(i,j,k-l,4))
            wlapz = wlapz + clapl*(w_aux_gpu(i,j,k+l,4)+w_aux_gpu(i,j,k-l,4))
            tz = tz+ccl*(w_aux_gpu(i,j,k+l,6)-w_aux_gpu(i,j,k-l,6))
            tlapz = tlapz + clapl*(w_aux_gpu(i,j,k+l,6)+w_aux_gpu(i,j,k-l,6))
            muz = muz+ccl*(w_aux_gpu(i,j,k+l,7)-w_aux_gpu(i,j,k-l,7))
          enddo
          ux = ux*dcsidx_gpu(i)
          vx = vx*dcsidx_gpu(i)
          wx = wx*dcsidx_gpu(i)
          tx = tx*dcsidx_gpu(i)
          mux = mux*dcsidx_gpu(i)
          uy = uy*detady_gpu(j)
          vy = vy*detady_gpu(j)
          wy = wy*detady_gpu(j)
          ty = ty*detady_gpu(j)
          muy = muy*detady_gpu(j)
          uz = uz*dzitdz_gpu(k)
          vz = vz*dzitdz_gpu(k)
          wz = wz*dzitdz_gpu(k)
          tz = tz*dzitdz_gpu(k)
          muz = muz*dzitdz_gpu(k)
          ulapx = ulapx*dcsidxs_gpu(i)+ux*dcsidx2_gpu(i)
          vlapx = vlapx*dcsidxs_gpu(i)+vx*dcsidx2_gpu(i)
          wlapx = wlapx*dcsidxs_gpu(i)+wx*dcsidx2_gpu(i)
          tlapx = tlapx*dcsidxs_gpu(i)+tx*dcsidx2_gpu(i)
          ulapy = ulapy*detadys_gpu(j)+uy*detady2_gpu(j)
          vlapy = vlapy*detadys_gpu(j)+vy*detady2_gpu(j)
          wlapy = wlapy*detadys_gpu(j)+wy*detady2_gpu(j)
          tlapy = tlapy*detadys_gpu(j)+ty*detady2_gpu(j)
          ulapz = ulapz*dzitdzs_gpu(k)+uz*dzitdz2_gpu(k)
          vlapz = vlapz*dzitdzs_gpu(k)+vz*dzitdz2_gpu(k)
          wlapz = wlapz*dzitdzs_gpu(k)+wz*dzitdz2_gpu(k)
          tlapz = tlapz*dzitdzs_gpu(k)+tz*dzitdz2_gpu(k)
          ulap = ulapx+ulapy+ulapz
          vlap = vlapx+vlapy+vlapz
          wlap = wlapx+wlapy+wlapz
          tlap = tlapx+tlapy+tlapz
          div = ux+vy+wz
          div3l = div/3._rkind
          w_aux_gpu(i,j,k,10) = div3l
          omegax = wy-vz
          omegay = uz-wx
          omegaz = vx-uy
          omod2 = omegax*omegax+omegay*omegay+omegaz*omegaz
          w_aux_gpu(i,j,k,9) = sqrt(omod2)
          if (calorically_perfect==1) then
            cploc = cp_coeff_gpu(0)
          else
            cploc = 0._rkind
            do ll=indx_cp_l,indx_cp_r
              cploc = cploc+cp_coeff_gpu(ll)*(tt/t0)**ll
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
          fl_gpu(i,j,k,2) = fl_gpu(i,j,k,2) - sigx
          fl_gpu(i,j,k,3) = fl_gpu(i,j,k,3) - sigy
          fl_gpu(i,j,k,4) = fl_gpu(i,j,k,4) - sigz
          fl_gpu(i,j,k,5) = fl_gpu(i,j,k,5) - sigq
          if (j==1) then
            wallprop_gpu(i,k,2) = mu*uy
            wallprop_gpu(i,k,3) = mu*wy
            wallprop_gpu(i,k,4) = mu*ty*cploc/prandtl
          endif
        enddo
      enddo
    enddo


  endsubroutine visflx_nosensor_kernel



  subroutine visflx_nosensor_c2_kernel(nx,ny,nz,ng,visc_order,prandtl,t0,indx_cp_l,indx_cp_r,&
  &cp_coeff_gpu,calorically_perfect,iblock,ite_rank_x,itu_rank_x,ite_l,itu_l,u0,l0,w_gpu,w_aux_gpu,&
  &fl_gpu,coeff_deriv1_gpu,coeff_deriv2_gpu,dcsidxc2_gpu,detadyc2_gpu,detadxc2_gpu,dcsidyc2_gpu,&
  &dzitdz_gpu,dzitdzs_gpu,dzitdz2_gpu,g1_gpu,g2_gpu,g12_gpu,jac_gpu,wallprop_gpu,vis_tag_gpu,&
  &wall_tag_gpu,ortho)
    integer, intent(in) :: nx, ny, nz, ng, visc_order, calorically_perfect
    integer, intent(in) :: indx_cp_l, indx_cp_r, ortho
    integer, intent(in) :: iblock, ite_rank_x, itu_rank_x, ite_l, itu_l
    real(rkind), intent(in) :: prandtl, u0, l0, t0
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(in) :: w_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(inout) :: w_aux_gpu
    real(rkind), dimension(1-ng:,1-ng:, 2:), intent(inout) :: wallprop_gpu
    real(rkind), dimension(1:,1:,1:, 1:), intent(inout) :: fl_gpu
    real(rkind), dimension(1:,1:), intent(in) :: coeff_deriv1_gpu
    real(rkind), dimension(0:,1:), intent(in) :: coeff_deriv2_gpu
    real(rkind), dimension(1:), intent(in) :: dzitdz_gpu, dzitdzs_gpu, dzitdz2_gpu
    real(rkind), dimension(1-ng:,1-ng:), intent(in) :: dcsidxc2_gpu, detadyc2_gpu
    real(rkind), dimension(1-ng:,1-ng:), intent(in) :: detadxc2_gpu, dcsidyc2_gpu
    real(rkind), dimension(1-ng:,1-ng:), intent(in) :: g1_gpu, g2_gpu, g12_gpu, jac_gpu
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in) :: cp_coeff_gpu
    integer, dimension(1-ng:), intent(in) :: vis_tag_gpu, wall_tag_gpu
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
    !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu,w_aux_gpu,wallprop_gpu,fl_gpu,coeff_deriv1_gpu,coeff_deriv2_gpu,dzitdz_gpu,dzitdzs_gpu,dzitdz2_gpu,dcsidxc2_gpu,detadyc2_gpu,detadxc2_gpu,dcsidyc2_gpu,g1_gpu,g2_gpu,g12_gpu,jac_gpu,cp_coeff_gpu,vis_tag_gpu,wall_tag_gpu)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          uu = w_aux_gpu(i,j,k,2)
          vv = w_aux_gpu(i,j,k,3)
          ww = w_aux_gpu(i,j,k,4)
          tt = w_aux_gpu(i,j,k,6)
          mu = w_aux_gpu(i,j,k,7)
          lmaxi = lmax
          lmaxj = lmax
          if (j == 1) lmaxi = vis_tag_gpu(i)
          if (wall_tag_gpu(i) < 1) lmaxj = min(j,lmax)
          if (calorically_perfect==1) then
            cploc = cp_coeff_gpu(0)
          else
            cploc = 0._rkind
            do ll=indx_cp_l,indx_cp_r
              cploc = cploc+cp_coeff_gpu(ll)*(tt/t0)**ll
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
            cli = coeff_deriv1_gpu(l,lmaxi)
            clj = coeff_deriv1_gpu(l,lmaxj)
            clk = coeff_deriv1_gpu(l,lmax )
            ucsi = ucsi +cli*(w_aux_gpu(i+l,j,k,2)-w_aux_gpu(i-l,j,k,2))
            vcsi = vcsi +cli*(w_aux_gpu(i+l,j,k,3)-w_aux_gpu(i-l,j,k,3))
            wcsi = wcsi +cli*(w_aux_gpu(i+l,j,k,4)-w_aux_gpu(i-l,j,k,4))
            tcsi = tcsi +cli*(w_aux_gpu(i+l,j,k,6)-w_aux_gpu(i-l,j,k,6))
            mucsi = mucsi+cli*(w_aux_gpu(i+l,j,k,7)-w_aux_gpu(i-l,j,k,7))
            ueta = ueta +clj*(w_aux_gpu(i,j+l,k,2)-w_aux_gpu(i,j-l,k,2))
            veta = veta +clj*(w_aux_gpu(i,j+l,k,3)-w_aux_gpu(i,j-l,k,3))
            weta = weta +clj*(w_aux_gpu(i,j+l,k,4)-w_aux_gpu(i,j-l,k,4))
            teta = teta +clj*(w_aux_gpu(i,j+l,k,6)-w_aux_gpu(i,j-l,k,6))
            mueta = mueta+clj*(w_aux_gpu(i,j+l,k,7)-w_aux_gpu(i,j-l,k,7))
            uzit = uzit +clk*(w_aux_gpu(i,j,k+l,2)-w_aux_gpu(i,j,k-l,2))
            vzit = vzit +clk*(w_aux_gpu(i,j,k+l,3)-w_aux_gpu(i,j,k-l,3))
            wzit = wzit +clk*(w_aux_gpu(i,j,k+l,4)-w_aux_gpu(i,j,k-l,4))
            tzit = tzit +clk*(w_aux_gpu(i,j,k+l,6)-w_aux_gpu(i,j,k-l,6))
            muzit = muzit+clk*(w_aux_gpu(i,j,k+l,7)-w_aux_gpu(i,j,k-l,7))
            dg1 = dg1 +cli*(g1_gpu (i+l,j)-g1_gpu (i-l,j))
            dg12csi=dg12csi+cli*(g12_gpu(i+l,j)-g12_gpu(i-l,j))
            dg2 = dg2 +clj*(g2_gpu (i,j+l)-g2_gpu (i,j-l))
            dg12eta=dg12eta+clj*(g12_gpu(i,j+l)-g12_gpu(i,j-l))
          enddo
          ux = ucsi *dcsidxc2_gpu(i,j) + ueta *detadxc2_gpu(i,j)
          vx = vcsi *dcsidxc2_gpu(i,j) + veta *detadxc2_gpu(i,j)
          wx = wcsi *dcsidxc2_gpu(i,j) + weta *detadxc2_gpu(i,j)
          tx = tcsi *dcsidxc2_gpu(i,j) + teta *detadxc2_gpu(i,j)
          mux = mucsi*dcsidxc2_gpu(i,j) + mueta*detadxc2_gpu(i,j)
          uy = ucsi *dcsidyc2_gpu(i,j) + ueta *detadyc2_gpu(i,j)
          vy = vcsi *dcsidyc2_gpu(i,j) + veta *detadyc2_gpu(i,j)
          wy = wcsi *dcsidyc2_gpu(i,j) + weta *detadyc2_gpu(i,j)
          ty = tcsi *dcsidyc2_gpu(i,j) + teta *detadyc2_gpu(i,j)
          muy = mucsi*dcsidyc2_gpu(i,j) + mueta*detadyc2_gpu(i,j)
          uz = uzit *dzitdz_gpu(k)
          vz = vzit *dzitdz_gpu(k)
          wz = wzit *dzitdz_gpu(k)
          tz = tzit *dzitdz_gpu(k)
          muz = muzit*dzitdz_gpu(k)
          if (j==1) then
            wallprop_gpu(i,k,3) = mu*weta
            wallprop_gpu(i,k,4) = mu*teta*cploc/prandtl
          endif
          ulapcsi = coeff_deriv2_gpu(0,lmaxi)*uu
          ulapeta = coeff_deriv2_gpu(0,lmaxj)*uu
          ulapzit = coeff_deriv2_gpu(0,lmax)*uu
          vlapcsi = coeff_deriv2_gpu(0,lmaxi)*vv
          vlapeta = coeff_deriv2_gpu(0,lmaxj)*vv
          vlapzit = coeff_deriv2_gpu(0,lmax)*vv
          wlapcsi = coeff_deriv2_gpu(0,lmaxi)*ww
          wlapeta = coeff_deriv2_gpu(0,lmaxj)*ww
          wlapzit = coeff_deriv2_gpu(0,lmax)*ww
          tlapcsi = coeff_deriv2_gpu(0,lmaxi)*tt
          tlapeta = coeff_deriv2_gpu(0,lmaxj)*tt
          tlapzit = coeff_deriv2_gpu(0,lmax)*tt
          do l=1,lmax
            clapi = coeff_deriv2_gpu(l,lmaxi)
            clapj = coeff_deriv2_gpu(l,lmaxj)
            clapk = coeff_deriv2_gpu(l,lmax )
            ulapcsi = ulapcsi + clapi*(w_aux_gpu(i+l,j,k,2)+w_aux_gpu(i-l,j,k,2))
            ulapeta = ulapeta + clapj*(w_aux_gpu(i,j+l,k,2)+w_aux_gpu(i,j-l,k,2))
            ulapzit = ulapzit + clapk*(w_aux_gpu(i,j,k+l,2)+w_aux_gpu(i,j,k-l,2))
            vlapcsi = vlapcsi + clapi*(w_aux_gpu(i+l,j,k,3)+w_aux_gpu(i-l,j,k,3))
            vlapeta = vlapeta + clapj*(w_aux_gpu(i,j+l,k,3)+w_aux_gpu(i,j-l,k,3))
            vlapzit = vlapzit + clapk*(w_aux_gpu(i,j,k+l,3)+w_aux_gpu(i,j,k-l,3))
            wlapcsi = wlapcsi + clapi*(w_aux_gpu(i+l,j,k,4)+w_aux_gpu(i-l,j,k,4))
            wlapeta = wlapeta + clapj*(w_aux_gpu(i,j+l,k,4)+w_aux_gpu(i,j-l,k,4))
            wlapzit = wlapzit + clapk*(w_aux_gpu(i,j,k+l,4)+w_aux_gpu(i,j,k-l,4))
            tlapcsi = tlapcsi + clapi*(w_aux_gpu(i+l,j,k,6)+w_aux_gpu(i-l,j,k,6))
            tlapeta = tlapeta + clapj*(w_aux_gpu(i,j+l,k,6)+w_aux_gpu(i,j-l,k,6))
            tlapzit = tlapzit + clapk*(w_aux_gpu(i,j,k+l,6)+w_aux_gpu(i,j,k-l,6))
          enddo
          if(ortho == 1) then
            ulap = g1_gpu(i,j)*ulapcsi+ucsi*dg1+g2_gpu(i,j)*ulapeta+ueta*dg2
            vlap = g1_gpu(i,j)*vlapcsi+vcsi*dg1+g2_gpu(i,j)*vlapeta+veta*dg2
            wlap = g1_gpu(i,j)*wlapcsi+wcsi*dg1+g2_gpu(i,j)*wlapeta+weta*dg2
            tlap = g1_gpu(i,j)*tlapcsi+tcsi*dg1+g2_gpu(i,j)*tlapeta+teta*dg2
          else
            ulapcsieta = 0.25_rkind*(w_aux_gpu(i+1,j+1,k,2)-w_aux_gpu(i-1,j+1,k,2)-w_aux_gpu(i+1,j-1,k,2)+w_aux_gpu(i-1,j-1,k,2))
            vlapcsieta = 0.25_rkind*(w_aux_gpu(i+1,j+1,k,3)-w_aux_gpu(i-1,j+1,k,3)-w_aux_gpu(i+1,j-1,k,3)+w_aux_gpu(i-1,j-1,k,3))
            wlapcsieta = 0.25_rkind*(w_aux_gpu(i+1,j+1,k,4)-w_aux_gpu(i-1,j+1,k,4)-w_aux_gpu(i+1,j-1,k,4)+w_aux_gpu(i-1,j-1,k,4))
            tlapcsieta = 0.25_rkind*(w_aux_gpu(i+1,j+1,k,6)-w_aux_gpu(i-1,j+1,k,6)-w_aux_gpu(i+1,j-1,k,6)+w_aux_gpu(i-1,j-1,k,6))
            if (iblock==ite_rank_x.and.i==ite_l-1.and.j==1) then
              ulapcsieta = 0.5_rkind*(w_aux_gpu(i,j+1,k,2)-w_aux_gpu(i-1,j+1,k,2)-w_aux_gpu(i,j-1,k,2)+w_aux_gpu(i-1,j-1,k,2))
              vlapcsieta = 0.5_rkind*(w_aux_gpu(i,j+1,k,3)-w_aux_gpu(i-1,j+1,k,3)-w_aux_gpu(i,j-1,k,3)+w_aux_gpu(i-1,j-1,k,3))
              wlapcsieta = 0.5_rkind*(w_aux_gpu(i,j+1,k,4)-w_aux_gpu(i-1,j+1,k,4)-w_aux_gpu(i,j-1,k,4)+w_aux_gpu(i-1,j-1,k,4))
              tlapcsieta = 0.5_rkind*(w_aux_gpu(i,j+1,k,6)-w_aux_gpu(i-1,j+1,k,6)-w_aux_gpu(i,j-1,k,6)+w_aux_gpu(i-1,j-1,k,6))
            endif
            if (iblock==itu_rank_x.and.i==itu_l+1.and.j==1) then
              ulapcsieta = 0.5_rkind*(w_aux_gpu(i+1,j+1,k,2)-w_aux_gpu(i,j+1,k,2)-w_aux_gpu(i+1,j-1,k,2)+w_aux_gpu(i,j-1,k,2))
              vlapcsieta = 0.5_rkind*(w_aux_gpu(i+1,j+1,k,3)-w_aux_gpu(i,j+1,k,3)-w_aux_gpu(i+1,j-1,k,3)+w_aux_gpu(i,j-1,k,3))
              wlapcsieta = 0.5_rkind*(w_aux_gpu(i+1,j+1,k,4)-w_aux_gpu(i,j+1,k,4)-w_aux_gpu(i+1,j-1,k,4)+w_aux_gpu(i,j-1,k,4))
              tlapcsieta = 0.5_rkind*(w_aux_gpu(i+1,j+1,k,6)-w_aux_gpu(i,j+1,k,6)-w_aux_gpu(i+1,j-1,k,6)+w_aux_gpu(i,j-1,k,6))
            endif
            ulap = g1_gpu(i,j)*ulapcsi+ucsi*dg1+2._rkind*g12_gpu(i,j)*ulapcsieta+dg12csi*ueta+&
            &dg12eta*ucsi+g2_gpu(i,j)*ulapeta+ueta*dg2
            vlap = g1_gpu(i,j)*vlapcsi+vcsi*dg1+2._rkind*g12_gpu(i,j)*vlapcsieta+dg12csi*veta+&
            &dg12eta*vcsi+g2_gpu(i,j)*vlapeta+veta*dg2
            wlap = g1_gpu(i,j)*wlapcsi+wcsi*dg1+2._rkind*g12_gpu(i,j)*wlapcsieta+dg12csi*weta+&
            &dg12eta*wcsi+g2_gpu(i,j)*wlapeta+weta*dg2
            tlap = g1_gpu(i,j)*tlapcsi+tcsi*dg1+2._rkind*g12_gpu(i,j)*tlapcsieta+dg12csi*teta+&
            &dg12eta*tcsi+g2_gpu(i,j)*tlapeta+teta*dg2
          endif
          ulap = ulap * jac_gpu(i,j)
          vlap = vlap * jac_gpu(i,j)
          wlap = wlap * jac_gpu(i,j)
          tlap = tlap * jac_gpu(i,j)
          ulapz = ulapzit*dzitdzs_gpu(k)+uz*dzitdz2_gpu(k)
          vlapz = vlapzit*dzitdzs_gpu(k)+vz*dzitdz2_gpu(k)
          wlapz = wlapzit*dzitdzs_gpu(k)+wz*dzitdz2_gpu(k)
          tlapz = tlapzit*dzitdzs_gpu(k)+tz*dzitdz2_gpu(k)
          ulap = ulap+ulapz
          vlap = vlap+vlapz
          wlap = wlap+wlapz
          tlap = tlap+tlapz
          div = ux+vy+wz
          div3l = div/3._rkind
          w_aux_gpu(i,j,k,10) = div3l
          omegax = wy-vz
          omegay = uz-wx
          omegaz = vx-uy
          omod2 = omegax*omegax+omegay*omegay+omegaz*omegaz
          w_aux_gpu(i,j,k,9) = sqrt(omod2)
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
          fl_gpu(i,j,k,2) = fl_gpu(i,j,k,2) - sigx
          fl_gpu(i,j,k,3) = fl_gpu(i,j,k,3) - sigy
          fl_gpu(i,j,k,4) = fl_gpu(i,j,k,4) - sigz
          fl_gpu(i,j,k,5) = fl_gpu(i,j,k,5) - sigq
        enddo
      enddo
    enddo


  endsubroutine visflx_nosensor_c2_kernel



  subroutine sponge_kernel(nx,ny,nz,ng,nv,w_gpu,wfar_gpu,fln_gpu,f_sponge_gpu,j_sponge)
    integer, intent(in) :: nx,ny,nz,ng,nv,j_sponge
    real(rkind), dimension(1-ng:, 1-ng:, 1-ng:, 1:), intent(in) :: w_gpu
    real(rkind), dimension(1:,1:,1:,1:), intent(inout) :: fln_gpu
    real(rkind), dimension(nx,nv), intent(in) :: wfar_gpu
    real(rkind), dimension(ny), intent(in) :: f_sponge_gpu
    integer :: i,j,k,m,iercuda

    !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu,fln_gpu,wfar_gpu,f_sponge_gpu)
    do k = 1,nz
      do j = j_sponge,ny
        do i = 1,nx
          do m=1,nv
            fln_gpu(i,j,k,m) = fln_gpu(i,j,k,m) - f_sponge_gpu(j)*(w_gpu(i,j,k,m) - wfar_gpu(i,m))
          enddo
        enddo
      enddo
    enddo


  endsubroutine sponge_kernel



  subroutine limiter_kernel(nx,ny,nz,ng,w_gpu,w_aux_gpu,iblock,kblock,indx_cp_l,indx_cp_r,&
  &cv_coeff_gpu,calorically_perfect,t0,tol_iter_nr,rho_lim,tem_lim,rho_lim_rescale,tem_lim_rescale)
    integer, intent(in) :: nx,ny,nz,ng,iblock,kblock,calorically_perfect,indx_cp_l,indx_cp_r
    real(rkind), intent(in) :: t0, tol_iter_nr
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(inout) :: w_gpu, w_aux_gpu
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in) :: cv_coeff_gpu
    real(rkind), intent(in) :: rho_lim, tem_lim, rho_lim_rescale, tem_lim_rescale
    real(rkind) :: rho, tem, eei, uu, vv, ww, rhoe, qq, ee
    integer :: i,j,k,iercuda
    real(rkind) :: n_limited_rho, n_limited_tem
    n_limited_rho = 0._rkind
    n_limited_tem = 0._rkind
    !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu,w_aux_gpu,cv_coeff_gpu) &
    !$omp& reduction(+:n_limited_rho,n_limited_tem)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          rho = w_gpu(i,j,k,1)
          uu = w_gpu(i,j,k,2)/rho
          vv = w_gpu(i,j,k,3)/rho
          ww = w_gpu(i,j,k,4)/rho
          rhoe = w_gpu(i,j,k,5)
          qq = 0.5_rkind*(uu*uu+vv*vv+ww*ww)
          ee = rhoe/rho-qq
          tem = get_temperature_from_e_dev(ee,w_aux_gpu(i,j,k,6),t0,cv_coeff_gpu,indx_cp_l,&
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
      !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu,w_aux_gpu,cv_coeff_gpu)
      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            rho = w_gpu(i,j,k,1)
            uu = w_gpu(i,j,k,2)/rho
            vv = w_gpu(i,j,k,3)/rho
            ww = w_gpu(i,j,k,4)/rho
            rhoe = w_gpu(i,j,k,5)
            qq = 0.5_rkind*(uu*uu+vv*vv+ww*ww)
            ee = rhoe/rho-qq
            tem = get_temperature_from_e_dev(ee,w_aux_gpu(i,j,k,6),t0,cv_coeff_gpu,indx_cp_l,&
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
              w_gpu(i,j,k,1) = rho
              w_gpu(i,j,k,2) = rho*uu
              w_gpu(i,j,k,3) = rho*vv
              w_gpu(i,j,k,4) = rho*ww
              eei = get_e_from_temperature_dev(tem, t0, indx_cp_l, indx_cp_r, cv_coeff_gpu, calorically_perfect)
              w_gpu(i,j,k,5) = rho*eei + 0.5_rkind*(w_gpu(i,j,k,2)**2+w_gpu(i,j,k,3)**2+w_gpu(i,j,k,4)**2)/rho
            endif
          enddo
        enddo
      enddo
    endif

  endsubroutine limiter_kernel



  subroutine filter_kernel(nx,ny,nz,ng,w_gpu,w_aux_gpu,indx_cp_l,indx_cp_r,cv_coeff_gpu,&
  &calorically_perfect,t0,tol_iter_nr,coeff_filter_gpu,jfilter,wall_tag_gpu)
    integer, intent(in) :: nx,ny,nz,ng,calorically_perfect,indx_cp_l,indx_cp_r,jfilter
    integer, dimension(1-ng:nx+ng), intent(in) :: wall_tag_gpu
    real(rkind), intent(in) :: t0, tol_iter_nr
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(inout) :: w_gpu, w_aux_gpu
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in) :: cv_coeff_gpu
    real(rkind), dimension(0:4), intent(in) :: coeff_filter_gpu
    real(rkind) :: cfilt,rho,eei,uu,vv,ww,rhoe,qq,ee,rho_p,rho_m,rho_f,tem_p,tem_m,tem_f
    integer :: i,j,k,l,iercuda

    !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu,w_aux_gpu,cv_coeff_gpu,coeff_filter_gpu,wall_tag_gpu)
    do k = 1,nz
      do j = 1+4,ny-4
        do i = 1,nx
          if (wall_tag_gpu(i) < 1 .and. j < jfilter) then
          else
            rho_f = 0._rkind
            tem_f = 0._rkind
            do l=0,4
              rho_p = w_gpu(i+l,j,k,1)
              uu = w_gpu(i+l,j,k,2)/rho_p
              vv = w_gpu(i+l,j,k,3)/rho_p
              ww = w_gpu(i+l,j,k,4)/rho_p
              rhoe = w_gpu(i+l,j,k,5)
              qq = 0.5_rkind*(uu*uu+vv*vv+ww*ww)
              ee = rhoe/rho_p-qq
              tem_p = get_temperature_from_e_dev(ee,w_aux_gpu(i+l,j,k,6),t0,cv_coeff_gpu,indx_cp_l,&
              &indx_cp_r,calorically_perfect,tol_iter_nr)
              rho_m = w_gpu(i-l,j,k,1)
              uu = w_gpu(i-l,j,k,2)/rho_m
              vv = w_gpu(i-l,j,k,3)/rho_m
              ww = w_gpu(i-l,j,k,4)/rho_m
              rhoe = w_gpu(i-l,j,k,5)
              qq = 0.5_rkind*(uu*uu+vv*vv+ww*ww)
              ee = rhoe/rho_m-qq
              tem_m = get_temperature_from_e_dev(ee,w_aux_gpu(i-l,j,k,6),t0,cv_coeff_gpu,indx_cp_l,&
              &indx_cp_r,calorically_perfect,tol_iter_nr)
              cfilt = coeff_filter_gpu(l)
              rho_f = rho_f + cfilt*(rho_p+rho_m)
              tem_f = tem_f + cfilt*(tem_p+tem_m)
            enddo
            uu = w_gpu(i,j,k,2)/w_gpu(i,j,k,1)
            vv = w_gpu(i,j,k,3)/w_gpu(i,j,k,1)
            ww = w_gpu(i,j,k,4)/w_gpu(i,j,k,1)
            w_gpu(i,j,k,1) = rho_f
            w_gpu(i,j,k,2) = rho_f*uu
            w_gpu(i,j,k,3) = rho_f*vv
            w_gpu(i,j,k,4) = rho_f*ww
            eei = get_e_from_temperature_dev(tem_f,t0,indx_cp_l,indx_cp_r,cv_coeff_gpu,calorically_perfect)
            w_gpu(i,j,k,5) = rho_f*eei + 0.5_rkind*(w_gpu(i,j,k,2)**2+w_gpu(i,j,k,3)**2+w_gpu(i,j,k,4)**2)/rho_f
            w_aux_gpu(i,j,k,6) = tem_f
          endif
        enddo
      enddo
    enddo

    !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu,w_aux_gpu,cv_coeff_gpu,coeff_filter_gpu,wall_tag_gpu)
    do k = 1,nz
      do j = 1+4,ny-4
        do i = 1,nx
          if (wall_tag_gpu(i) < 1 .and. j < jfilter) then
          else
            rho_f = 0._rkind
            tem_f = 0._rkind
            do l=0,4
              rho_p = w_gpu(i,j+l,k,1)
              uu = w_gpu(i,j+l,k,2)/rho_p
              vv = w_gpu(i,j+l,k,3)/rho_p
              ww = w_gpu(i,j+l,k,4)/rho_p
              rhoe = w_gpu(i,j+l,k,5)
              qq = 0.5_rkind*(uu*uu+vv*vv+ww*ww)
              ee = rhoe/rho_p-qq
              tem_p = get_temperature_from_e_dev(ee,w_aux_gpu(i+l,j,k,6),t0,cv_coeff_gpu,indx_cp_l,&
              &indx_cp_r,calorically_perfect,tol_iter_nr)
              rho_m = w_gpu(i,j-l,k,1)
              uu = w_gpu(i,j-l,k,2)/rho_m
              vv = w_gpu(i,j-l,k,3)/rho_m
              ww = w_gpu(i,j-l,k,4)/rho_m
              rhoe = w_gpu(i,j-l,k,5)
              qq = 0.5_rkind*(uu*uu+vv*vv+ww*ww)
              ee = rhoe/rho_m-qq
              tem_m = get_temperature_from_e_dev(ee,w_aux_gpu(i-l,j,k,6),t0,cv_coeff_gpu,indx_cp_l,&
              &indx_cp_r,calorically_perfect,tol_iter_nr)
              cfilt = coeff_filter_gpu(l)
              rho_f = rho_f + cfilt*(rho_p+rho_m)
              tem_f = tem_f + cfilt*(tem_p+tem_m)
            enddo
            uu = w_gpu(i,j,k,2)/w_gpu(i,j,k,1)
            vv = w_gpu(i,j,k,3)/w_gpu(i,j,k,1)
            ww = w_gpu(i,j,k,4)/w_gpu(i,j,k,1)
            w_gpu(i,j,k,1) = rho_f
            w_gpu(i,j,k,2) = rho_f*uu
            w_gpu(i,j,k,3) = rho_f*vv
            w_gpu(i,j,k,4) = rho_f*ww
            eei = get_e_from_temperature_dev(tem_f,t0,indx_cp_l,indx_cp_r,cv_coeff_gpu,calorically_perfect)
            w_gpu(i,j,k,5) = rho_f*eei + 0.5_rkind*(w_gpu(i,j,k,2)**2+w_gpu(i,j,k,3)**2+w_gpu(i,j,k,4)**2)/rho_f
            w_aux_gpu(i,j,k,6) = tem_f
          endif
        enddo
      enddo
    enddo


  endsubroutine filter_kernel



  subroutine visflx_reduced_ord2_kernel(nx,ny,nz,ng,prandtl,t0,indx_cp_l,indx_cp_r,cp_coeff_gpu,&
  &calorically_perfect,u0,l0,w_gpu,w_aux_gpu,fl_gpu,x_gpu,y_gpu,z_gpu,wallprop_gpu,update_sensor)
    integer, intent(in) :: nx, ny, nz, ng, calorically_perfect
    integer, intent(in) :: indx_cp_l, indx_cp_r
    real(rkind), intent(in) :: prandtl, u0, l0, t0
    integer, intent(in) :: update_sensor
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(in) :: w_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(inout) :: w_aux_gpu
    real(rkind), dimension(1-ng:,1-ng:, 2:), intent(inout) :: wallprop_gpu
    real(rkind), dimension(1:,1:,1:, 1:), intent(inout) :: fl_gpu
    real(rkind), dimension(1-ng:), intent(in) :: x_gpu, y_gpu, z_gpu
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in) :: cp_coeff_gpu
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

    !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu,w_aux_gpu,wallprop_gpu,fl_gpu,x_gpu,y_gpu,z_gpu,cp_coeff_gpu)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          uu = w_aux_gpu(i,j,k,2)
          vv = w_aux_gpu(i,j,k,3)
          ww = w_aux_gpu(i,j,k,4)
          tt = w_aux_gpu(i,j,k,6)
          mu = w_aux_gpu(i,j,k,7)
          if (calorically_perfect==1) then
            cploc = cp_coeff_gpu(0)
          else
            cploc = 0._rkind
            do ll=indx_cp_l,indx_cp_r
              cploc = cploc+cp_coeff_gpu(ll)*(tt/t0)**ll
            enddo
          endif
          dxl = 1._rkind/(x_gpu(i+1)-x_gpu(i-1))
          dyl = 1._rkind/(y_gpu(j+1)-y_gpu(j-1))
          dzl = 1._rkind/(z_gpu(k+1)-z_gpu(k-1))
          ux = dxl*(w_aux_gpu(i+1,j,k,2)-w_aux_gpu(i-1,j,k,2))
          vx = dxl*(w_aux_gpu(i+1,j,k,3)-w_aux_gpu(i-1,j,k,3))
          wx = dxl*(w_aux_gpu(i+1,j,k,4)-w_aux_gpu(i-1,j,k,4))
          mux = dxl*(w_aux_gpu(i+1,j,k,7)-w_aux_gpu(i-1,j,k,7))
          uy = dyl*(w_aux_gpu(i,j+1,k,2)-w_aux_gpu(i,j-1,k,2))
          vy = dyl*(w_aux_gpu(i,j+1,k,3)-w_aux_gpu(i,j-1,k,3))
          wy = dyl*(w_aux_gpu(i,j+1,k,4)-w_aux_gpu(i,j-1,k,4))
          ty = dyl*(w_aux_gpu(i,j+1,k,6)-w_aux_gpu(i,j-1,k,6))
          muy = dyl*(w_aux_gpu(i,j+1,k,7)-w_aux_gpu(i,j-1,k,7))
          uz = dzl*(w_aux_gpu(i,j,k+1,2)-w_aux_gpu(i,j,k-1,2))
          vz = dzl*(w_aux_gpu(i,j,k+1,3)-w_aux_gpu(i,j,k-1,3))
          wz = dzl*(w_aux_gpu(i,j,k+1,4)-w_aux_gpu(i,j,k-1,4))
          muz = dzl*(w_aux_gpu(i,j,k+1,7)-w_aux_gpu(i,j,k-1,7))
          if (j==1) then
            wallprop_gpu(i,k,2) = mu*uy
            wallprop_gpu(i,k,3) = mu*wy
            wallprop_gpu(i,k,4) = mu*ty*cploc/prandtl
          endif
          div = ux+vy+wz
          div3l = div/3._rkind
          w_aux_gpu(i,j,k,10) = div3l
          omegax = wy-vz
          omegay = uz-wx
          omegaz = vx-uy
          omod2 = omegax*omegax+omegay*omegay+omegaz*omegaz
          w_aux_gpu(i,j,k,9) = sqrt(omod2)
          if (update_sensor == 1) then
            w_aux_gpu(i,j,k,8) = (max(-div/sqrt(omod2+div**2+(u0/l0)**2),0._rkind))**2
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
          fl_gpu(i,j,k,2) = fl_gpu(i,j,k,2) - sigx
          fl_gpu(i,j,k,3) = fl_gpu(i,j,k,3) - sigy
          fl_gpu(i,j,k,4) = fl_gpu(i,j,k,4) - sigz
          fl_gpu(i,j,k,5) = fl_gpu(i,j,k,5) - sigq
        enddo
      enddo
    enddo


  endsubroutine visflx_reduced_ord2_kernel




  subroutine sensor_kernel(nx,ny,nz,ng,u0,l0,w_aux_gpu)
    integer, intent(in) :: nx, ny, nz, ng
    real(rkind), intent(in) :: u0, l0
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(inout) :: w_aux_gpu
    integer :: i,j,k,iercuda
    real(rkind) :: div, omod2

    !$omp target teams distribute parallel do collapse(3) has_device_addr(w_aux_gpu)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          omod2 = w_aux_gpu(i,j,k,9)**2
          div = 3._rkind*w_aux_gpu(i,j,k,10)
          w_aux_gpu(i,j,k,8) = (max(-div/sqrt(omod2+div**2+(u0/l0)**2),0._rkind))**2
        enddo
      enddo
    enddo


  endsubroutine sensor_kernel



  subroutine sensor_c2_kernel(nx,ny,nz,ng,u0,l0_ducros,w_aux_gpu,iblock,ite_rank_x,itu_rank_x,ite_l,&
  &itu_l,teshk,theta_ij_gpu,theta_threshold,jweno)
    integer, intent(in) :: nx, ny, nz, ng, iblock, ite_rank_x, itu_rank_x, ite_l, itu_l, jweno
    real(rkind), intent(in) :: u0, l0_ducros, teshk, theta_threshold
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(inout) :: w_aux_gpu
    real(rkind), dimension(1:,1:), intent(inout) :: theta_ij_gpu
    integer :: i,j,k,iercuda
    real(rkind) :: div, omod2

    !$omp target teams distribute parallel do collapse(3) has_device_addr(w_aux_gpu,theta_ij_gpu)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          omod2 = w_aux_gpu(i,j,k,9)**2
          div = 3._rkind*w_aux_gpu(i,j,k,10)
          w_aux_gpu(i,j,k,8) = max(-div/sqrt(omod2+div**2+(u0/l0_ducros)**2),0._rkind)
          if( j == 1) then
            if(iblock == ite_rank_x .and. i == ite_l) w_aux_gpu(i,j,k,8) = teshk
            if(iblock == itu_rank_x .and. i == itu_l) w_aux_gpu(i,j,k,8) = teshk
          endif
          if(theta_ij_gpu(i,j) > theta_threshold .or. j > jweno) then
            w_aux_gpu(i,j,k,8) = 1000._rkind
          endif
        enddo
      enddo
    enddo


  endsubroutine sensor_c2_kernel



  subroutine visflx_x_kernel(nx,ny,nz,nv,ng,prandtl,t0,indx_cp_l,indx_cp_r,cp_coeff_gpu,&
  &calorically_perfect,x_gpu,w_aux_gpu,fl_gpu,fhat_gpu)
    integer, intent(in) :: nx, ny, nz, nv, ng, calorically_perfect
    integer, intent(in) :: indx_cp_l, indx_cp_r
    real(rkind), intent(in) :: prandtl, t0
    real(rkind), dimension(1:nx,1:ny,1:nz,1:nv), intent(inout) :: fl_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in) :: w_aux_gpu
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv), intent(out) :: fhat_gpu
    real(rkind), dimension(1-ng:), intent(in) :: x_gpu
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in) :: cp_coeff_gpu
    integer :: i,j,k,iv,ll,iercuda
    real(rkind) :: uu,vv,ww,tt,mu,qq
    real(rkind) :: uup,vvp,wwp,ttp,mup,qqp
    real(rkind) :: sigq,sigx,sigy,sigz,sigq_tt,sigq_qq
    real(rkind) :: dxhl,fl2o,fl3o,fl4o,fl5o
    real(rkind) :: ttf,muf
    real(rkind) :: cploc

    !$omp target teams distribute parallel do collapse(3) has_device_addr(fl_gpu,w_aux_gpu,fhat_gpu,x_gpu,cp_coeff_gpu)
    do k = 1,nz
      do j = 1,ny
        do i = 0,nx
          uu = w_aux_gpu(i ,j,k,2)
          uup = w_aux_gpu(i+1,j,k,2)
          vv = w_aux_gpu(i ,j,k,3)
          vvp = w_aux_gpu(i+1,j,k,3)
          ww = w_aux_gpu(i ,j,k,4)
          wwp = w_aux_gpu(i+1,j,k,4)
          tt = w_aux_gpu(i ,j,k,6)
          ttp = w_aux_gpu(i+1,j,k,6)
          mu = w_aux_gpu(i ,j,k,7)
          mup = w_aux_gpu(i+1,j,k,7)
          qq = 0.5_rkind*(uu*uu+vv*vv+ww*ww)
          qqp = 0.5_rkind*(uup*uup+vvp*vvp+wwp*wwp)
          if (calorically_perfect==1) then
            cploc = cp_coeff_gpu(0)
          else
            ttf = 0.5_rkind*(tt+ttp)
            cploc = 0._rkind
            do ll=indx_cp_l,indx_cp_r
              cploc = cploc+cp_coeff_gpu(ll)*(ttf/t0)**ll
            enddo
          endif
          sigx = uup-uu
          sigy = vvp-vv
          sigz = wwp-ww
          sigq_tt = ttp-tt
          sigq_qq = qqp-qq
          muf = mu+mup
          muf = 0.5_rkind*muf/(x_gpu(i+1)-x_gpu(i))
          sigx = sigx*muf
          sigy = sigy*muf
          sigz = sigz*muf
          sigq = (sigq_tt*cploc/prandtl+sigq_qq)*muf
          fhat_gpu(i,j,k,2) = - sigx
          fhat_gpu(i,j,k,3) = - sigy
          fhat_gpu(i,j,k,4) = - sigz
          fhat_gpu(i,j,k,5) = - sigq
        enddo
      enddo
    enddo

    !$omp target teams distribute parallel do collapse(3) has_device_addr(fl_gpu,w_aux_gpu,fhat_gpu,x_gpu,cp_coeff_gpu)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          dxhl = 2._rkind/(x_gpu(i+1)-x_gpu(i-1))
          do iv=2,nv
            fl_gpu(i,j,k,iv) = fl_gpu(i,j,k,iv) + dxhl*(fhat_gpu(i,j,k,iv)-fhat_gpu(i-1,j,k,iv))
          enddo
        enddo
      enddo
    enddo


  endsubroutine visflx_x_kernel



  subroutine visflx_y_kernel(nx,ny,nz,nv,ng,prandtl,t0,indx_cp_l,indx_cp_r,cp_coeff_gpu,calorically_perfect,y_gpu,w_aux_gpu,fl_gpu)
    integer, intent(in) :: nx, ny, nz, nv, ng, calorically_perfect
    integer, intent(in) :: indx_cp_l, indx_cp_r
    real(rkind), intent(in) :: prandtl, t0
    real(rkind), dimension(1:nx,1:ny,1:nz,1:nv), intent(inout) :: fl_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in) :: w_aux_gpu
    real(rkind), dimension(1-ng:), intent(in) :: y_gpu
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in) :: cp_coeff_gpu
    integer :: i,j,k,iv,ll,iercuda
    real(rkind) :: uu,vv,ww,tt,mu,qq
    real(rkind) :: uup,vvp,wwp,ttp,mup,qqp
    real(rkind) :: sigq,sigx,sigy,sigz,sigq_tt,sigq_qq
    real(rkind) :: dyhl,fl2o,fl3o,fl4o,fl5o
    real(rkind) :: ttf,muf
    real(rkind) :: cploc

    !$omp target teams distribute parallel do collapse(2) has_device_addr(fl_gpu,w_aux_gpu,y_gpu,cp_coeff_gpu)
    do k = 1,nz
      do i = 1,nx
        do j=0,ny
          uu = w_aux_gpu(i,j ,k,2)
          uup = w_aux_gpu(i,j+1,k,2)
          vv = w_aux_gpu(i,j ,k,3)
          vvp = w_aux_gpu(i,j+1,k,3)
          ww = w_aux_gpu(i,j ,k,4)
          wwp = w_aux_gpu(i,j+1,k,4)
          tt = w_aux_gpu(i,j ,k,6)
          ttp = w_aux_gpu(i,j+1,k,6)
          mu = w_aux_gpu(i,j ,k,7)
          mup = w_aux_gpu(i,j+1,k,7)
          qq = 0.5_rkind*(uu*uu+vv*vv+ww*ww)
          qqp = 0.5_rkind*(uup*uup+vvp*vvp+wwp*wwp)
          if (calorically_perfect==1) then
            cploc = cp_coeff_gpu(0)
          else
            ttf = 0.5_rkind*(tt+ttp)
            cploc = 0._rkind
            do ll=indx_cp_l,indx_cp_r
              cploc = cploc+cp_coeff_gpu(ll)*(ttf/t0)**ll
            enddo
          endif
          sigx = uup-uu
          sigy = vvp-vv
          sigz = wwp-ww
          sigq_tt = ttp-tt
          sigq_qq = qqp-qq
          muf = mu+mup
          muf = 0.5_rkind*muf/(y_gpu(j+1)-y_gpu(j))
          sigx = sigx*muf
          sigy = sigy*muf
          sigz = sigz*muf
          sigq = (sigq_tt*cploc/prandtl+sigq_qq)*muf
          if (j>0) then
            fl_gpu(i,j,k,2) = fl_gpu(i,j,k,2) + fl2o-sigx*dyhl
            fl_gpu(i,j,k,3) = fl_gpu(i,j,k,3) + fl3o-sigy*dyhl
            fl_gpu(i,j,k,4) = fl_gpu(i,j,k,4) + fl4o-sigz*dyhl
            fl_gpu(i,j,k,5) = fl_gpu(i,j,k,5) + fl5o-sigq*dyhl
          endif
          if (j<ny) then
            dyhl = 2._rkind/(y_gpu(j+2)-y_gpu(j))
            fl2o = sigx*dyhl
            fl3o = sigy*dyhl
            fl4o = sigz*dyhl
            fl5o = sigq*dyhl
          endif
        enddo
      enddo
    enddo


  endsubroutine visflx_y_kernel


  subroutine visflx_z_kernel(nx,ny,nz,nv,ng,prandtl,t0,indx_cp_l,indx_cp_r,cp_coeff_gpu,calorically_perfect,z_gpu,w_aux_gpu,fl_gpu)
    integer, intent(in) :: nx, ny, nz, nv, ng, calorically_perfect
    integer, intent(in) :: indx_cp_l, indx_cp_r
    real(rkind), intent(in) :: prandtl, t0
    real(rkind), dimension(1:nx,1:ny,1:nz,1:nv), intent(inout) :: fl_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in) :: w_aux_gpu
    real(rkind), dimension(1-ng:), intent(in) :: z_gpu
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in) :: cp_coeff_gpu
    integer :: i,j,k,iv,ll,iercuda
    real(rkind) :: uu,vv,ww,tt,mu,qq
    real(rkind) :: uup,vvp,wwp,ttp,mup,qqp
    real(rkind) :: sigq,sigx,sigy,sigz,sigq_tt,sigq_qq
    real(rkind) :: dzhl,fl2o,fl3o,fl4o,fl5o
    real(rkind) :: ttf,muf
    real(rkind) :: cploc

    !$omp target teams distribute parallel do collapse(2) has_device_addr(fl_gpu,w_aux_gpu,z_gpu,cp_coeff_gpu)
    do j = 1,ny
      do i = 1,nx
        do k=0,nz
          uu = w_aux_gpu(i,j,k ,2)
          uup = w_aux_gpu(i,j,k+1,2)
          vv = w_aux_gpu(i,j,k ,3)
          vvp = w_aux_gpu(i,j,k+1,3)
          ww = w_aux_gpu(i,j,k ,4)
          wwp = w_aux_gpu(i,j,k+1,4)
          tt = w_aux_gpu(i,j,k ,6)
          ttp = w_aux_gpu(i,j,k+1,6)
          mu = w_aux_gpu(i,j,k ,7)
          mup = w_aux_gpu(i,j,k+1,7)
          qq = 0.5_rkind*(uu*uu+vv*vv+ww*ww)
          qqp = 0.5_rkind*(uup*uup+vvp*vvp+wwp*wwp)
          if (calorically_perfect==1) then
            cploc = cp_coeff_gpu(0)
          else
            ttf = 0.5_rkind*(tt+ttp)
            cploc = 0._rkind
            do ll=indx_cp_l,indx_cp_r
              cploc = cploc+cp_coeff_gpu(ll)*(ttf/t0)**ll
            enddo
          endif
          sigx = uup-uu
          sigy = vvp-vv
          sigz = wwp-ww
          sigq_tt = ttp-tt
          sigq_qq = qqp-qq
          muf = mu+mup
          muf = 0.5_rkind*muf/(z_gpu(k+1)-z_gpu(k))
          sigx = sigx*muf
          sigy = sigy*muf
          sigz = sigz*muf
          sigq = (sigq_tt*cploc/prandtl+sigq_qq)*muf
          if (k>0) then
            fl_gpu(i,j,k,2) = fl_gpu(i,j,k,2) + fl2o-sigx*dzhl
            fl_gpu(i,j,k,3) = fl_gpu(i,j,k,3) + fl3o-sigy*dzhl
            fl_gpu(i,j,k,4) = fl_gpu(i,j,k,4) + fl4o-sigz*dzhl
            fl_gpu(i,j,k,5) = fl_gpu(i,j,k,5) + fl5o-sigq*dzhl
          endif
          if (k<nz) then
            dzhl = 2._rkind/(z_gpu(k+2)-z_gpu(k))
            fl2o = sigx*dzhl
            fl3o = sigy*dzhl
            fl4o = sigz*dzhl
            fl5o = sigq*dzhl
          endif
        enddo
      enddo
    enddo


  endsubroutine visflx_z_kernel



  subroutine recyc_exchange_kernel_1(irecyc,w_gpu,wbuf1s_gpu,nx,ny,nz,ng,nv)
    integer, intent(in) :: irecyc, nx, ny, nz, ng, nv
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in) :: w_gpu
    real(rkind), dimension(:,:,:,:), intent(inout) :: wbuf1s_gpu
    integer :: i,j,k,m, iercuda

    !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu,wbuf1s_gpu)
    do k = 1,nz
      do j = 1,ny
        do i = 1,ng
          do m=1,nv
            wbuf1s_gpu(i,j,k,m) = w_gpu(irecyc+1-i,j,k,m)
          enddo
        enddo
      enddo
    enddo


  endsubroutine recyc_exchange_kernel_1


  subroutine recyc_exchange_kernel_2(n1_start_recv,n1_start_send,n1_end_recv,wrecyc_gpu,wbuf1r_gpu,nx,ny,nz,ng,nv)
    integer, intent(in) :: n1_start_recv, n1_start_send, n1_end_recv, nx, ny, nz, ng, nv
    real(rkind), dimension(:,:,:,:), intent(in) :: wbuf1r_gpu
    real(rkind), dimension(:,:,:,:), intent(inout) :: wrecyc_gpu
    integer :: i,j,k,m, iercuda

    !$omp target teams distribute parallel do collapse(3) has_device_addr(wbuf1r_gpu,wrecyc_gpu)
    do k = n1_start_recv,n1_end_recv
      do j = 1,ny
        do i = 1,ng
          do m=1,nv
            wrecyc_gpu(i,j,k,m) = wbuf1r_gpu(i,j,k-n1_start_recv+n1_start_send,m)
          enddo
        enddo
      enddo
    enddo


  endsubroutine recyc_exchange_kernel_2


  subroutine recyc_exchange_kernel_3(n2_start_recv,n2_start_send,n2_end_recv,wrecyc_gpu,wbuf2r_gpu,nx,ny,nz,ng,nv)
    integer, intent(in) :: n2_start_recv, n2_start_send, n2_end_recv, nx, ny, nz, ng, nv
    real(rkind), dimension(:,:,:,:), intent(in) :: wbuf2r_gpu
    real(rkind), dimension(:,:,:,:), intent(inout) :: wrecyc_gpu
    integer :: i,j,k,m, iercuda

    !$omp target teams distribute parallel do collapse(3) has_device_addr(wbuf2r_gpu,wrecyc_gpu)
    do k = n2_start_recv,n2_end_recv
      do j = 1,ny
        do i = 1,ng
          do m=1,nv
            wrecyc_gpu(i,j,k,m) = wbuf2r_gpu(i,j,k-n2_start_recv+n2_start_send,m)
          enddo
        enddo
      enddo
    enddo


  endsubroutine recyc_exchange_kernel_3



  subroutine bcextr_sub_kernel(ilat,nx,ny,nz,ng,p0,rgas0,w_gpu,indx_cp_l,indx_cp_r,cv_coeff_gpu,t0,calorically_perfect)
    integer, intent(in) :: ilat, nx, ny, nz, ng, indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind), intent(in) :: p0, t0, rgas0
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in) :: cv_coeff_gpu
    real(rkind), dimension(1-ng:, 1-ng:, 1-ng:, 1:), intent(inout) :: w_gpu
    real(rkind) :: rho,rhou,rhov,rhow,tt,ee
    integer :: i,j,k,l,m,iercuda
    if (ilat==1) then
      !$omp target teams distribute parallel do collapse(3) has_device_addr(cv_coeff_gpu,w_gpu)
      do k = 1,nz
        do j = 1,ny
          do l = 1,ng
            rho = w_gpu(1,j,k,1)
            rhou = w_gpu(1,j,k,2)
            rhov = w_gpu(1,j,k,3)
            rhow = w_gpu(1,j,k,4)
            tt = p0/rho/rgas0
            ee = get_e_from_temperature_dev(tt,t0,indx_cp_l,indx_cp_r,cv_coeff_gpu,calorically_perfect)
            w_gpu(1-l,j,k,1) = rho
            w_gpu(1-l,j,k,2) = rhou
            w_gpu(1-l,j,k,3) = rhov
            w_gpu(1-l,j,k,4) = rhow
            w_gpu(1-l,j,k,5) = rho*ee + 0.5_rkind*(rhou**2+rhov**2+rhow**2)/rho
          enddo
        enddo
      enddo
    elseif (ilat==2) then
      !$omp target teams distribute parallel do collapse(3) has_device_addr(cv_coeff_gpu,w_gpu)
      do k = 1,nz
        do j = 1,ny
          do l = 1,ng
            rho = w_gpu(nx,j,k,1)
            rhou = w_gpu(nx,j,k,2)
            rhov = w_gpu(nx,j,k,3)
            rhow = w_gpu(nx,j,k,4)
            tt = p0/rho/rgas0
            ee = get_e_from_temperature_dev(tt,t0,indx_cp_l,indx_cp_r,cv_coeff_gpu,calorically_perfect)
            w_gpu(nx+l,j,k,1) = rho
            w_gpu(nx+l,j,k,2) = rhou
            w_gpu(nx+l,j,k,3) = rhov
            w_gpu(nx+l,j,k,4) = rhow
            w_gpu(nx+l,j,k,5) = rho*ee + 0.5_rkind*(rhou**2+rhov**2+rhow**2)/rho
          enddo
        enddo
      enddo
    elseif (ilat==3) then
    elseif (ilat==4) then
      !$omp target teams distribute parallel do collapse(3) has_device_addr(cv_coeff_gpu,w_gpu)
      do k = 1,nz
        do i = 1,nx
          do l = 1,ng
            rho = w_gpu(i,ny,k,1)
            rhou = w_gpu(i,ny,k,2)
            rhov = w_gpu(i,ny,k,3)
            rhow = w_gpu(i,ny,k,4)
            tt = p0/rho/rgas0
            ee = get_e_from_temperature_dev(tt, t0, indx_cp_l, indx_cp_r, cv_coeff_gpu, calorically_perfect)
            w_gpu(i,ny+l,k,1) = rho
            w_gpu(i,ny+l,k,2) = rhou
            w_gpu(i,ny+l,k,3) = rhov
            w_gpu(i,ny+l,k,4) = rhow
            w_gpu(i,ny+l,k,5) = rho*ee + 0.5_rkind*(rhou**2+rhov**2+rhow**2)/rho
          enddo
        enddo
      enddo
    elseif (ilat==5) then
      !$omp target teams distribute parallel do collapse(3) has_device_addr(cv_coeff_gpu,w_gpu)
      do l = 1,ng
        do j = 1,ny
          do i = 1,nx
            rho = w_gpu(i,j,1,1)
            rhou = w_gpu(i,j,1,2)
            rhov = w_gpu(i,j,1,3)
            rhow = w_gpu(i,j,1,4)
            tt = p0/rho/rgas0
            ee = get_e_from_temperature_dev(tt, t0, indx_cp_l, indx_cp_r, cv_coeff_gpu, calorically_perfect)
            w_gpu(i,j,1-l,1) = rho
            w_gpu(i,j,1-l,2) = rhou
            w_gpu(i,j,1-l,3) = rhov
            w_gpu(i,j,1-l,4) = rhow
            w_gpu(i,j,1-l,5) = rho*ee + 0.5_rkind*(rhou**2+rhov**2+rhow**2)/rho
          enddo
        enddo
      enddo
    elseif (ilat==6) then
      !$omp target teams distribute parallel do collapse(3) has_device_addr(cv_coeff_gpu,w_gpu)
      do l = 1,ng
        do j = 1,ny
          do i = 1,nx
            rho = w_gpu(i,j,nz,1)
            rhou = w_gpu(i,j,nz,2)
            rhov = w_gpu(i,j,nz,3)
            rhow = w_gpu(i,j,nz,4)
            tt = p0/rho/rgas0
            ee = get_e_from_temperature_dev(tt, t0, indx_cp_l, indx_cp_r, cv_coeff_gpu, calorically_perfect)
            w_gpu(i,j,nz+l,1) = rho
            w_gpu(i,j,nz+l,2) = rhou
            w_gpu(i,j,nz+l,3) = rhov
            w_gpu(i,j,nz+l,4) = rhow
            w_gpu(i,j,nz+l,5) = rho*ee + 0.5_rkind*(rhou**2+rhov**2+rhow**2)/rho
          enddo
        enddo
      enddo
    endif

  endsubroutine bcextr_sub_kernel




  subroutine bc_nr_lat_x_kernel(start_or_end,nr_type,nx,ny,nz,ng,nv,w_aux_gpu,w_gpu,fl_gpu,&
  &dcsidx_gpu,indx_cp_l,indx_cp_r,cp_coeff_gpu,winf_gpu,calorically_perfect,rgas0,t0)
    integer, intent(in), value :: start_or_end, nr_type, nx, ny, nz, ng, nv, indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind), intent(in), value :: rgas0,t0
    real(rkind), dimension(3) :: c_one
    real(rkind), dimension(5) :: dw_dn, dwc_dn, ev, dw_dn_outer, dwc_dn_outer
    real(rkind), dimension(5,5) :: el, er
    real(rkind) :: w_target
    integer :: i, j, k, l, m, mm, sgn_dw
    real(rkind) :: df, uu, vv, ww, h, qq, cc, c, ci, b2, b1
    real(rkind), dimension(indx_cp_l:indx_cp_r+1) :: cp_coeff_gpu
    real(rkind), dimension(:,:,:,:) :: fl_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,:) :: w_aux_gpu, w_gpu
    real(rkind), dimension(:) :: dcsidx_gpu
    real(rkind), dimension(:) :: winf_gpu
    real(rkind) :: rho,tt,gamloc,p_rho,p_e,etot,b3
    !$omp target data map(alloc:c_one,dw_dn,dwc_dn,ev,dw_dn_outer,dwc_dn_outer,el,er)
    !$omp target teams distribute parallel do collapse(2) has_device_addr(cp_coeff_gpu,fl_gpu,w_aux_gpu,w_gpu,dcsidx_gpu,winf_gpu) private(c_one,dw_dn,dwc_dn,ev,dw_dn_outer,dwc_dn_outer,el,er)
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
            dw_dn(m) = dw_dn(m) + sgn_dw * c_one(l)*w_gpu(i+sgn_dw*(l-1),j,k,m)
          enddo
          w_target = 0._rkind
          if (nr_type == 2) then
            w_target = w_gpu(i-sgn_dw,j,k,m)
          elseif(nr_type == 3) then
            w_target = winf_gpu(m)
          endif
          dw_dn_outer(m) = sgn_dw * (w_gpu(i,j,k,m)-w_target)
        enddo

        rho = w_aux_gpu(i,j,k,1)
        uu = w_aux_gpu(i,j,k,2)
        vv = w_aux_gpu(i,j,k,3)
        ww = w_aux_gpu(i,j,k,4)
        h = w_aux_gpu(i,j,k,5)
        tt = w_aux_gpu(i,j,k,6)
        qq = 0.5_rkind * (uu*uu +vv*vv + ww*ww)
        gamloc = get_gamloc_dev(indx_cp_l,indx_cp_r,tt,t0,cp_coeff_gpu,calorically_perfect,rgas0)
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
          fl_gpu(i,j,k,m) = fl_gpu(i,j,k,m) + df * dcsidx_gpu(i)
        enddo
      enddo
    enddo

    !$omp end target data
  endsubroutine bc_nr_lat_x_kernel


  subroutine bc_nr_lat_x_c2_kernel(start_or_end,nr_type,nx,ny,nz,ng,nv,w_aux_gpu,w_gpu,fl_gpu,&
  &dcsidx_gpu,indx_cp_l,indx_cp_r,cp_coeff_gpu,winf_gpu,jac_gpu,mcsi_gpu,dcsidxnc2_gpu,dcsidync2_gpu,&
  &dcsidxc2_gpu,dcsidyc2_gpu,calorically_perfect,rgas0,t0)
    integer, intent(in), value :: start_or_end, nr_type, nx, ny, nz, ng, nv, indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind), intent(in), value :: rgas0,t0
    real(rkind), dimension(3) :: c_one
    real(rkind), dimension(5) :: dw_dn, dwc_dn, ev, dw_dn_outer, dwc_dn_outer
    real(rkind), dimension(5,5) :: el, er
    real(rkind) :: w_target
    integer :: i, j, k, l, m, mm, sgn_dw
    real(rkind) :: df, uu, vv, ww, h, qq, cc, c, ci, b2, b1
    real(rkind), dimension(indx_cp_l:indx_cp_r+1) :: cp_coeff_gpu
    real(rkind), dimension(:,:,:,:) :: fl_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,:) :: w_aux_gpu, w_gpu
    real(rkind), dimension(:) :: dcsidx_gpu
    real(rkind), dimension(:) :: winf_gpu
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng), intent(in) :: jac_gpu, mcsi_gpu, dcsidxnc2_gpu, dcsidync2_gpu
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng), intent(in) :: dcsidxc2_gpu, dcsidyc2_gpu
    real(rkind) :: dcsixjdcsi, dcsiyjdcsi, pp, ut
    real(rkind) :: rho,tt,gamloc,p_rho,p_e,etot,b3
    !$omp target data map(alloc:c_one,dw_dn,dwc_dn,ev,dw_dn_outer,dwc_dn_outer,el,er)
    !$omp target teams distribute parallel do collapse(2) has_device_addr(cp_coeff_gpu,fl_gpu,w_aux_gpu,w_gpu,dcsidx_gpu,winf_gpu,jac_gpu,mcsi_gpu,dcsidxnc2_gpu,dcsidync2_gpu,dcsidxc2_gpu,dcsidyc2_gpu) private(c_one,dw_dn,dwc_dn,ev,dw_dn_outer,dwc_dn_outer,el,er)
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
            dw_dn(m) = dw_dn(m) + sgn_dw * c_one(l)*w_gpu(i+sgn_dw*(l-1),j,k,m)
          enddo
          w_target = 0._rkind
          if (nr_type == 2) then
            w_target = w_gpu(i-sgn_dw,j,k,m)
          elseif(nr_type == 3) then
            w_target = winf_gpu(m)
          endif
          dw_dn_outer(m) = sgn_dw * (w_gpu(i,j,k,m)-w_target)
        enddo

        rho = w_aux_gpu(i,j,k,1)
        uu = w_aux_gpu(i,j,k,2)
        vv = w_aux_gpu(i,j,k,3)
        ww = w_aux_gpu(i,j,k,4)
        h = w_aux_gpu(i,j,k,5)
        tt = w_aux_gpu(i,j,k,6)
        qq = 0.5_rkind * (uu*uu +vv*vv + ww*ww)
        gamloc = get_gamloc_dev(indx_cp_l,indx_cp_r,tt,t0,cp_coeff_gpu,calorically_perfect,rgas0)
        cc = gamloc * tt * rgas0
        c = sqrt(cc)
        ci = 1._rkind/c
        p_rho = tt*rgas0
        p_e = rho*(gamloc-1._rkind)
        etot = h - tt*rgas0
        b3 = etot - rho * p_rho/p_e
        b2 = p_e/(rho*cc)
        b1 = p_rho/cc - b2*(etot - 2._rkind*qq)
        ut = dcsidxnc2_gpu(i,j) * uu + dcsidync2_gpu(i,j) * vv
        call eigenvectors_x_c2(b1, b2, b3, uu, vv, ww, c, ci, h, el, er, ut, dcsidxnc2_gpu(i,j), dcsidync2_gpu(i,j))

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
          fl_gpu(i,j,k,m) = fl_gpu(i,j,k,m) + df * mcsi_gpu(i,j)
        enddo


        dcsixjdcsi = sgn_dw * (c_one(1)*dcsidxc2_gpu(i,j)/jac_gpu(i,j)+c_one(2)*dcsidxc2_gpu(i+1*&
        &sgn_dw,j)/jac_gpu(i+1*sgn_dw,j)+ c_one(3)*dcsidxc2_gpu(i+2*sgn_dw,j)/jac_gpu(i+2*sgn_dw,j))
        dcsiyjdcsi = sgn_dw * (c_one(1)*dcsidyc2_gpu(i,j)/jac_gpu(i,j)+c_one(2)*dcsidyc2_gpu(i+1*&
        &sgn_dw,j)/jac_gpu(i+1*sgn_dw,j)+ c_one(3)*dcsidyc2_gpu(i+2*sgn_dw,j)/jac_gpu(i+2*sgn_dw,j))
        pp = p_rho * rho
        fl_gpu(i,j,k,1) = fl_gpu(i,j,k,1)+( w_gpu(i,j,k,1)*uu *dcsixjdcsi+ w_gpu(i,j,k,1)*vv *dcsiyjdcsi)*jac_gpu(i,j)
        fl_gpu(i,j,k,2) = fl_gpu(i,j,k,2)+((w_gpu(i,j,k,2)*uu+pp)*dcsixjdcsi+ w_gpu(i,j,k,2)*vv *dcsiyjdcsi)*jac_gpu(i,j)
        fl_gpu(i,j,k,3) = fl_gpu(i,j,k,3)+( w_gpu(i,j,k,3)*uu *dcsixjdcsi+(w_gpu(i,j,k,3)*vv+pp)*dcsiyjdcsi)*jac_gpu(i,j)
        fl_gpu(i,j,k,4) = fl_gpu(i,j,k,4)+( w_gpu(i,j,k,4)*uu *dcsixjdcsi+ w_gpu(i,j,k,4)*vv *dcsiyjdcsi)*jac_gpu(i,j)
        fl_gpu(i,j,k,5) = fl_gpu(i,j,k,5)+((w_gpu(i,j,k,5)+pp)*uu*dcsixjdcsi+(w_gpu(i,j,k,5)+pp)*vv*dcsiyjdcsi)*jac_gpu(i,j)
      enddo
    enddo

    !$omp end target data
  endsubroutine bc_nr_lat_x_c2_kernel


  subroutine bc_nr_lat_y_kernel(start_or_end,nr_type,nx,ny,nz,ng,nv,w_aux_gpu,w_gpu,fl_gpu,&
  &detady_gpu,indx_cp_l,indx_cp_r,cp_coeff_gpu,winf_gpu,calorically_perfect,rgas0,t0)
    integer, intent(in), value :: start_or_end, nr_type, nx, ny, nz, ng, nv, indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind), intent(in), value :: rgas0,t0
    real(rkind), dimension(3) :: c_one
    real(rkind), dimension(5) :: dw_dn, dwc_dn, ev, dw_dn_outer, dwc_dn_outer
    real(rkind), dimension(5,5) :: el, er
    real(rkind) :: w_target
    integer :: i, j, k, l, m, mm, sgn_dw
    real(rkind) :: df, uu, vv, ww, h, qq, cc, c, ci, b2, b1
    real(rkind), dimension(indx_cp_l:indx_cp_r+1) :: cp_coeff_gpu
    real(rkind), dimension(:,:,:,:) :: fl_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: w_aux_gpu, w_gpu
    real(rkind), dimension(:) :: detady_gpu
    real(rkind), dimension(:) :: winf_gpu
    real(rkind) :: rho,tt,gamloc,p_rho,p_e,etot,b3
    !$omp target data map(alloc:c_one,dw_dn,dwc_dn,ev,dw_dn_outer,dwc_dn_outer,el,er)
    !$omp target teams distribute parallel do collapse(2) has_device_addr(cp_coeff_gpu,fl_gpu,w_aux_gpu,w_gpu,detady_gpu,winf_gpu) private(c_one,dw_dn,dwc_dn,ev,dw_dn_outer,dwc_dn_outer,el,er)
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
            dw_dn(m) = dw_dn(m) + sgn_dw * c_one(l)*w_gpu(i,j+sgn_dw*(l-1),k,m)
          enddo
          w_target = 0._rkind
          if (nr_type == 2) then
            w_target = w_gpu(i,j-sgn_dw,k,m)
          elseif(nr_type == 3) then
            w_target = winf_gpu(m)
          endif
          dw_dn_outer(m) = sgn_dw * (w_gpu(i,j,k,m)-w_target)
        enddo

        rho = w_aux_gpu(i,j,k,1)
        uu = w_aux_gpu(i,j,k,2)
        vv = w_aux_gpu(i,j,k,3)
        ww = w_aux_gpu(i,j,k,4)
        h = w_aux_gpu(i,j,k,5)
        tt = w_aux_gpu(i,j,k,6)
        qq = 0.5_rkind * (uu*uu +vv*vv + ww*ww)
        gamloc = get_gamloc_dev(indx_cp_l,indx_cp_r,tt,t0,cp_coeff_gpu,calorically_perfect,rgas0)
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
          fl_gpu(i,j,k,m) = fl_gpu(i,j,k,m) + df * detady_gpu(j)
        enddo
      enddo
    enddo

    !$omp end target data
  endsubroutine bc_nr_lat_y_kernel


  subroutine bc_nr_lat_y_c2_kernel(start_or_end,nr_type,nx,ny,nz,ng,nv,w_aux_gpu,w_gpu,fl_gpu,&
  &detady_gpu,indx_cp_l,indx_cp_r,cp_coeff_gpu,wfar_gpu,jac_gpu,meta_gpu,detadxnc2_gpu,detadync2_gpu,&
  &detadxc2_gpu,detadyc2_gpu,wall_tag_gpu,calorically_perfect,rgas0,t0)
    integer, intent(in), value :: start_or_end, nr_type, nx, ny, nz, ng, nv, indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind), intent(in), value :: rgas0,t0
    real(rkind), dimension(3) :: c_one
    real(rkind), dimension(5) :: dw_dn, dwc_dn, ev, dw_dn_outer, dwc_dn_outer
    real(rkind), dimension(5,5) :: el, er
    real(rkind) :: w_target
    integer :: i, j, k, l, m, mm, sgn_dw
    real(rkind) :: df, uu, vv, ww, h, qq, cc, c, ci, b2, b1
    real(rkind), dimension(indx_cp_l:indx_cp_r+1) :: cp_coeff_gpu
    real(rkind), dimension(:,:,:,:) :: fl_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: w_aux_gpu, w_gpu
    real(rkind), dimension(:) :: detady_gpu
    real(rkind), dimension(nx,nv) :: wfar_gpu
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng), intent(in) :: jac_gpu, meta_gpu, detadxnc2_gpu, detadync2_gpu
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng), intent(in) :: detadxc2_gpu, detadyc2_gpu
    integer, dimension(1-ng:nx+ng), intent(in) :: wall_tag_gpu
    real(rkind) :: rho,tt,gamloc,p_rho,p_e,etot,b3
    real(rkind) :: detaxjdeta, detayjdeta, pp, ut
    !$omp target data map(alloc:c_one,dw_dn,dwc_dn,ev,dw_dn_outer,dwc_dn_outer,el,er)
    !$omp target teams distribute parallel do collapse(2) has_device_addr(cp_coeff_gpu,fl_gpu,w_aux_gpu,w_gpu,detady_gpu,wfar_gpu,jac_gpu,meta_gpu,detadxnc2_gpu,detadync2_gpu,detadxc2_gpu,detadyc2_gpu,wall_tag_gpu) private(c_one,dw_dn,dwc_dn,ev,dw_dn_outer,dwc_dn_outer,el,er)
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

        if(((start_or_end == 1) .and. (wall_tag_gpu(i) > 0)) ) then
          sgn_dw = 1
        else
          do m=1,nv
            dw_dn(m) = 0._rkind
            do l=1,3
              dw_dn(m) = dw_dn(m) + sgn_dw * c_one(l)*w_gpu(i,j+sgn_dw*(l-1),k,m)
            enddo
            w_target = 0._rkind
            if (nr_type == 2) then
              w_target = w_gpu(i,j-sgn_dw,k,m)
            elseif(nr_type == 3) then
              w_target = wfar_gpu(i,m)
            endif
            dw_dn_outer(m) = sgn_dw * (w_gpu(i,j,k,m)-w_target)
          enddo

          rho = w_aux_gpu(i,j,k,1)
          uu = w_aux_gpu(i,j,k,2)
          vv = w_aux_gpu(i,j,k,3)
          ww = w_aux_gpu(i,j,k,4)
          h = w_aux_gpu(i,j,k,5)
          tt = w_aux_gpu(i,j,k,6)
          qq = 0.5_rkind * (uu*uu +vv*vv + ww*ww)
          gamloc = get_gamloc_dev(indx_cp_l,indx_cp_r,tt,t0,cp_coeff_gpu,calorically_perfect,rgas0)
          cc = gamloc * tt * rgas0
          c = sqrt(cc)
          ci = 1._rkind/c
          p_rho = tt*rgas0
          p_e = rho*(gamloc-1._rkind)
          etot = h - tt*rgas0
          b3 = etot - rho * p_rho/p_e
          b2 = p_e/(rho*cc)
          b1 = p_rho/cc - b2*(etot - 2._rkind*qq)
          ut = detadxnc2_gpu(i,j) * uu + detadync2_gpu(i,j) * vv
          call eigenvectors_y_c2(b1, b2, b3, uu, vv, ww, c, ci, h, el, er, ut, detadxnc2_gpu(i,j), detadync2_gpu(i,j))

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
            fl_gpu(i,j,k,m) = fl_gpu(i,j,k,m) + df * meta_gpu(i,j)
          enddo

          detaxjdeta=sgn_dw*(c_one(1)*detadxc2_gpu(i,j)/jac_gpu(i,j)+c_one(2)*detadxc2_gpu(i,j+1*&
          &sgn_dw)/jac_gpu(i,j+1*sgn_dw)+ c_one(3)*detadxc2_gpu(i,j+2*sgn_dw)/jac_gpu(i,j+2*sgn_dw))
          detayjdeta=sgn_dw*(c_one(1)*detadyc2_gpu(i,j)/jac_gpu(i,j)+c_one(2)*detadyc2_gpu(i,j+1*&
          &sgn_dw)/jac_gpu(i,j+1*sgn_dw)+ c_one(3)*detadyc2_gpu(i,j+2*sgn_dw)/jac_gpu(i,j+2*sgn_dw))

          pp = p_rho * rho
          fl_gpu(i,j,k,1) = fl_gpu(i,j,k,1) +( w_gpu(i,j,k,1)*uu *detaxjdeta + w_gpu(i,j,k,1)*vv *detayjdeta)*jac_gpu(i,j)
          fl_gpu(i,j,k,2) = fl_gpu(i,j,k,2) +((w_gpu(i,j,k,2)*uu+pp)*detaxjdeta + w_gpu(i,j,k,2)*vv *detayjdeta)*jac_gpu(i,j)
          fl_gpu(i,j,k,3) = fl_gpu(i,j,k,3) +( w_gpu(i,j,k,3)*uu *detaxjdeta + (w_gpu(i,j,k,3)*vv+pp)*detayjdeta)*jac_gpu(i,j)
          fl_gpu(i,j,k,4) = fl_gpu(i,j,k,4) +( w_gpu(i,j,k,4)*uu *detaxjdeta + w_gpu(i,j,k,4)*vv *detayjdeta)*jac_gpu(i,j)
          fl_gpu(i,j,k,5) = fl_gpu(i,j,k,5) +((w_gpu(i,j,k,5)+pp)*uu*detaxjdeta + (w_gpu(i,j,k,5)+pp)*vv*detayjdeta)*jac_gpu(i,j)
        endif
      enddo
    enddo

    !$omp end target data
  endsubroutine bc_nr_lat_y_c2_kernel


  subroutine bc_nr_lat_z_kernel(start_or_end,nr_type,nx,ny,nz,ng,nv,w_aux_gpu,w_gpu,fl_gpu,&
  &dzitdz_gpu,indx_cp_l,indx_cp_r,cp_coeff_gpu,winf_gpu,calorically_perfect,rgas0,t0)
    integer, intent(in), value :: start_or_end, nr_type, nx, ny, nz, ng, nv, indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind), intent(in), value :: rgas0,t0
    real(rkind), dimension(3) :: c_one
    real(rkind), dimension(5) :: dw_dn, dwc_dn, ev, dw_dn_outer, dwc_dn_outer
    real(rkind), dimension(5,5) :: el, er
    real(rkind) :: w_target
    integer :: i, j, k, l, m, mm, sgn_dw
    real(rkind) :: df, uu, vv, ww, h, qq, cc, c, ci, b2, b1
    real(rkind), dimension(indx_cp_l:indx_cp_r+1) :: cp_coeff_gpu
    real(rkind), dimension(:,:,:,:) :: fl_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: w_aux_gpu, w_gpu
    real(rkind), dimension(:) :: dzitdz_gpu
    real(rkind), dimension(:) :: winf_gpu
    real(rkind) :: rho,tt,gamloc,p_rho,p_e,etot,b3
    !$omp target data map(alloc:c_one,dw_dn,dwc_dn,ev,dw_dn_outer,dwc_dn_outer,el,er)
    !$omp target teams distribute parallel do collapse(2) has_device_addr(cp_coeff_gpu,fl_gpu,w_aux_gpu,w_gpu,dzitdz_gpu,winf_gpu) private(c_one,dw_dn,dwc_dn,ev,dw_dn_outer,dwc_dn_outer,el,er)
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
            dw_dn(m) = dw_dn(m) + sgn_dw * c_one(l)*w_gpu(i,j,k+sgn_dw*(l-1),m)
          enddo
          if (nr_type == 2) then
            w_target = w_gpu(i,j,k-sgn_dw,m)
          elseif(nr_type == 3) then
            w_target = winf_gpu(m)
          endif
          dw_dn_outer(m) = sgn_dw * (w_gpu(i,j,k,m)-w_target)
        enddo

        rho = w_aux_gpu(i,j,k,1)
        uu = w_aux_gpu(i,j,k,2)
        vv = w_aux_gpu(i,j,k,3)
        ww = w_aux_gpu(i,j,k,4)
        h = w_aux_gpu(i,j,k,5)
        tt = w_aux_gpu(i,j,k,6)
        qq = 0.5_rkind * (uu*uu +vv*vv + ww*ww)
        gamloc = get_gamloc_dev(indx_cp_l,indx_cp_r,tt,t0,cp_coeff_gpu,calorically_perfect,rgas0)
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
          fl_gpu(i,j,k,m) = fl_gpu(i,j,k,m) + df * dzitdz_gpu(k)
        enddo
      enddo
    enddo

    !$omp end target data
  endsubroutine bc_nr_lat_z_kernel


  subroutine bcrecyc_kernel_1(nx,ny,nz,ng,nv,wrecycav_gpu,wrecyc_gpu)
    integer, intent(in) :: nx, ny, nz, ng, nv
    real(rkind), dimension(:,:,:), intent(inout) :: wrecycav_gpu
    real(rkind), dimension(:,:,:,:), intent(in) :: wrecyc_gpu
    integer :: i,j,k,m,iercuda

    !$omp target teams distribute parallel do collapse(2) has_device_addr(wrecycav_gpu,wrecyc_gpu)
    do j = 1,ny
      do i = 1,ng
        do m=1,nv
          wrecycav_gpu(i,j,m) = 0._rkind
          do k=1,nz
            wrecycav_gpu(i,j,m) = wrecycav_gpu(i,j,m)+wrecyc_gpu(i,j,k,m)
          enddo
        enddo
      enddo
    enddo


  endsubroutine bcrecyc_kernel_1



  subroutine bcrecyc_kernel_2(nx,ny,nz,nzmax,ng,wrecycav_gpu,wrecyc_gpu)
    integer, intent(in) :: nx, ny, nz, nzmax, ng
    real(rkind), dimension(:,:,:), intent(in) :: wrecycav_gpu
    real(rkind), dimension(:,:,:,:), intent(inout) :: wrecyc_gpu
    real(rkind) :: ufav, vfav, wfav, rhom
    integer :: i,j,k,iercuda

    !$omp target teams distribute parallel do collapse(2) has_device_addr(wrecycav_gpu,wrecyc_gpu)
    do j = 1,ny
      do i = 1,ng
        ufav = wrecycav_gpu(i,j,2)/wrecycav_gpu(i,j,1)
        vfav = wrecycav_gpu(i,j,3)/wrecycav_gpu(i,j,1)
        wfav = wrecycav_gpu(i,j,4)/wrecycav_gpu(i,j,1)
        rhom = wrecycav_gpu(i,j,1)/nzmax
        do k=1,nz
          wrecyc_gpu(i,j,k,2) = wrecyc_gpu(i,j,k,2)/wrecyc_gpu(i,j,k,1)-ufav
          wrecyc_gpu(i,j,k,3) = wrecyc_gpu(i,j,k,3)/wrecyc_gpu(i,j,k,1)-vfav
          wrecyc_gpu(i,j,k,4) = wrecyc_gpu(i,j,k,4)/wrecyc_gpu(i,j,k,1)-wfav
          wrecyc_gpu(i,j,k,1) = wrecyc_gpu(i,j,k,1)-rhom
        enddo
      enddo
    enddo


  endsubroutine bcrecyc_kernel_2




  subroutine bcrecyc_kernel_3(nx,ny,nz,ng,p0,u0,rgas0,w_gpu,wmean_gpu,wrecyc_gpu,weta_inflow_gpu,&
  &map_j_inn_gpu,map_j_out_gpu,map_j_out_blend_gpu,yplus_inflow_gpu,eta_inflow_gpu,yplus_recyc_gpu,&
  &eta_recyc_gpu,eta_recyc_blend_gpu,betarecyc,glund1,inflow_random_plane_gpu,indx_cp_l,indx_cp_r,&
  &cv_coeff_gpu,cp_coeff_gpu,t0,calorically_perfect,rand_type)
    integer, intent(in) :: nx, ny, nz, ng, indx_cp_l, indx_cp_r, calorically_perfect, rand_type
    real(rkind) :: p0, rgas0, betarecyc, glund1, t0, u0
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in) :: cv_coeff_gpu, cp_coeff_gpu
    real(rkind), dimension(1-ng:,:,:), intent(in) :: wmean_gpu
    real(rkind), dimension(:,:,:,:), intent(in) :: wrecyc_gpu
    real(rkind), dimension(:,:,:), intent(in) :: inflow_random_plane_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,:), intent(inout) :: w_gpu
    real(rkind), dimension(:), intent(in) :: weta_inflow_gpu
    real(rkind), dimension(1-ng:), intent(in) :: yplus_inflow_gpu, eta_inflow_gpu, yplus_recyc_gpu, eta_recyc_gpu
    real(rkind), dimension(1-ng:), intent(in) :: eta_recyc_blend_gpu
    integer, dimension(:), intent(in) :: map_j_inn_gpu, map_j_out_gpu, map_j_out_blend_gpu
    integer :: i,j,k,iercuda
    real(rkind) :: eta, weta, weta1, bdamp, disty_inn, disty_out, rhofluc, ufluc, vfluc, wfluc, rhof_inn, rhof_out
    real(rkind) :: uf_inn, uf_out, vf_inn, vf_out, wf_inn, wf_out, etamin
    real(rkind) :: rhomean, uumean , vvmean , wwmean , tmean , rho , uu , vv , ww , rhou , rhov , rhow, tt, ee
    real(rkind) :: u0_02, tfluc
    integer :: j_inn, j_out
    u0_02 = 0.02_rkind*u0
    if (rand_type==0) u0_02 = 0._rkind
    !$omp target teams distribute parallel do collapse(2) has_device_addr(cv_coeff_gpu,cp_coeff_gpu,wmean_gpu,wrecyc_gpu,inflow_random_plane_gpu,w_gpu,weta_inflow_gpu,yplus_inflow_gpu,eta_inflow_gpu,yplus_recyc_gpu,eta_recyc_gpu,eta_recyc_blend_gpu,map_j_inn_gpu,map_j_out_gpu,map_j_out_blend_gpu)
    do k = 1,nz
      do j = 1,ny
        eta = eta_inflow_gpu(j)
        etamin = min(eta,1._rkind)
        weta = weta_inflow_gpu(j)
        weta1 = 1._rkind-weta
        bdamp = 0.5_rkind*(1._rkind-tanh(4._rkind*(eta_inflow_gpu(j)-2._rkind)))
        j_inn = map_j_inn_gpu(j)
        j_out = map_j_out_gpu(j)
        disty_inn = (yplus_inflow_gpu(j)-yplus_recyc_gpu(j_inn))/(yplus_recyc_gpu(j_inn+1)-yplus_recyc_gpu(j_inn))
        disty_out = (eta_inflow_gpu(j)-eta_recyc_gpu(j_out))/(eta_recyc_gpu(j_out+1)-eta_recyc_gpu(j_out))
        do i=1,ng
          rhomean = wmean_gpu(1-i,j,1)
          uumean = wmean_gpu(1-i,j,2)/rhomean
          vvmean = wmean_gpu(1-i,j,3)/rhomean
          wwmean = wmean_gpu(1-i,j,4)/rhomean
          tmean = p0/rhomean/rgas0
          if (j==1.or.j_inn>=ny.or.j_out>=ny) then
            rhofluc = 0._rkind
            ufluc = 0._rkind
            vfluc = 0._rkind
            wfluc = 0._rkind
          else
            rhof_inn = wrecyc_gpu(i,j_inn,k,1)*(1._rkind-disty_inn)+wrecyc_gpu(i,j_inn+1,k,1)*disty_inn
            rhof_out = wrecyc_gpu(i,j_out,k,1)*(1._rkind-disty_out)+wrecyc_gpu(i,j_out+1,k,1)*disty_out
            uf_inn = wrecyc_gpu(i,j_inn,k,2)*(1._rkind-disty_inn)+wrecyc_gpu(i,j_inn+1,k,2)*disty_inn
            uf_out = wrecyc_gpu(i,j_out,k,2)*(1._rkind-disty_out)+wrecyc_gpu(i,j_out+1,k,2)*disty_out
            vf_inn = wrecyc_gpu(i,j_inn,k,3)*(1._rkind-disty_inn)+wrecyc_gpu(i,j_inn+1,k,3)*disty_inn
            vf_out = wrecyc_gpu(i,j_out,k,3)*(1._rkind-disty_out)+wrecyc_gpu(i,j_out+1,k,3)*disty_out
            wf_inn = wrecyc_gpu(i,j_inn,k,4)*(1._rkind-disty_inn)+wrecyc_gpu(i,j_inn+1,k,4)*disty_inn
            wf_out = wrecyc_gpu(i,j_out,k,4)*(1._rkind-disty_out)+wrecyc_gpu(i,j_out+1,k,4)*disty_out
            rhofluc = rhof_inn*weta1+rhof_out*weta
            ufluc = uf_inn*weta1+ uf_out*weta
            vfluc = vf_inn*weta1+ vf_out*weta
            wfluc = wf_inn*weta1+ wf_out*weta
            rhofluc = rhofluc*bdamp
            ufluc = ufluc *bdamp*betarecyc
            vfluc = vfluc *bdamp*betarecyc
            wfluc = wfluc *bdamp*betarecyc
            ufluc = ufluc+u0_02*(inflow_random_plane_gpu(j,k,1)-0.5_rkind)*etamin
            vfluc = vfluc+u0_02*(inflow_random_plane_gpu(j,k,2)-0.5_rkind)*etamin
            wfluc = wfluc+u0_02*(inflow_random_plane_gpu(j,k,3)-0.5_rkind)*etamin
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
          ee = get_e_from_temperature_dev(tt, t0, indx_cp_l, indx_cp_r, cv_coeff_gpu, calorically_perfect)
          w_gpu(1-i,j,k,1) = rho
          w_gpu(1-i,j,k,2) = rhou
          w_gpu(1-i,j,k,3) = rhov
          w_gpu(1-i,j,k,4) = rhow
          w_gpu(1-i,j,k,5) = rho*ee + 0.5_rkind*(rhou**2+rhov**2+rhow**2)/rho
        enddo
      enddo
    enddo


  endsubroutine bcrecyc_kernel_3



  subroutine bclam_kernel(ilat,nx,ny,nz,ng,nv,w_gpu,wmean_gpu,p0,rgas0,indx_cp_l,indx_cp_r,cv_coeff_gpu,t0,calorically_perfect)
    integer, intent(in) :: nx, ny, nz, ng, nv, indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind), intent(in) :: p0, rgas0, t0
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in) :: cv_coeff_gpu
    real(rkind), dimension(1-ng:, 1-ng:, 1-ng:, 1:), intent(inout) :: w_gpu
    real(rkind), dimension(1-ng:,:,:), intent(in) :: wmean_gpu
    integer :: ilat
    integer :: j,k,l,iercuda
    real(rkind) :: rho,rhou,rhov,rhow,tt,ee
    if (ilat==1) then
      !$omp target teams distribute parallel do collapse(3) has_device_addr(cv_coeff_gpu,w_gpu,wmean_gpu)
      do k = 1,nz
        do j = 1,ny
          do l = 1,ng
            rho = wmean_gpu(1-l,j,1)
            rhou = wmean_gpu(1-l,j,2)
            rhov = wmean_gpu(1-l,j,3)
            rhow = wmean_gpu(1-l,j,4)
            tt = p0/rho/rgas0
            w_gpu(1-l,j,k,1) = rho
            w_gpu(1-l,j,k,2) = rhou
            w_gpu(1-l,j,k,3) = rhov
            w_gpu(1-l,j,k,4) = rhow
            ee = get_e_from_temperature_dev(tt, t0, indx_cp_l, indx_cp_r, cv_coeff_gpu, calorically_perfect)
            w_gpu(1-l,j,k,5) = rho*ee + 0.5_rkind*(rhou**2+rhov**2+rhow**2)/rho
          enddo
        enddo
      enddo
    elseif (ilat==2) then
    elseif (ilat==3) then
    elseif (ilat==4) then
    elseif (ilat==5) then
    elseif (ilat==6) then
    endif

  endsubroutine bclam_kernel



  subroutine bcfree_kernel(ilat,nx,ny,nz,ng,nv,winf_gpu,w_gpu)
    integer, intent(in) :: ilat,nx,ny,nz,ng,nv
    real(rkind), dimension(1:), intent(in) :: winf_gpu
    real(rkind), dimension(1-ng:, 1-ng:, 1-ng:, 1:), intent(inout) :: w_gpu
    integer :: i,j,k,l,m,iercuda
    if (ilat==1) then
      !$omp target teams distribute parallel do collapse(3) has_device_addr(winf_gpu,w_gpu)
      do k = 1,nz
        do j = 1,ny
          do l = 1,ng
            do m=1,nv
              w_gpu(1-l,j,k,m) = winf_gpu(m)
            enddo
          enddo
        enddo
      enddo
    elseif (ilat==2) then
      !$omp target teams distribute parallel do collapse(2) has_device_addr(winf_gpu,w_gpu)
      do k = 1,nz
        do j = 1,ny
          do l=1,ng
            do m=1,nv
              w_gpu(nx+l,j,k,m) = winf_gpu(m)
            enddo
          enddo
        enddo
      enddo
    elseif (ilat==3) then
      !$omp target teams distribute parallel do collapse(2) has_device_addr(winf_gpu,w_gpu)
      do k = 1,nz
        do i = 1,nx
          do l=1,ng
            do m=1,nv
              w_gpu(i,1-l,k,m) = winf_gpu(m)
            enddo
          enddo
        enddo
      enddo
    elseif (ilat==4) then
      !$omp target teams distribute parallel do collapse(2) has_device_addr(winf_gpu,w_gpu)
      do k = 1,nz
        do i = 1,nx
          do l=1,ng
            do m=1,nv
              w_gpu(i,ny+l,k,m) = winf_gpu(m)
            enddo
          enddo
        enddo
      enddo
    elseif (ilat==5) then
      !$omp target teams distribute parallel do collapse(2) has_device_addr(winf_gpu,w_gpu)
      do j = 1,ny
        do i = 1,nx
          do l=1,ng
            do m=1,nv
              w_gpu(i,j,1-l,m) = winf_gpu(m)
            enddo
          enddo
        enddo
      enddo
    elseif (ilat==6) then
      !$omp target teams distribute parallel do collapse(2) has_device_addr(winf_gpu,w_gpu)
      do j = 1,ny
        do i = 1,nx
          do l=1,ng
            do m=1,nv
              w_gpu(i,j,nz+l,m) = winf_gpu(m)
            enddo
          enddo
        enddo
      enddo
    endif

  endsubroutine bcfree_kernel



  subroutine bcfree_sub_kernel(ilat,nx,ny,nz,ng,w_gpu,w_aux_gpu,aoa,t0,ptot0,ttot0,rgas0,indx_cp_l,&
  &indx_cp_r,cv_coeff_gpu,cp_coeff_gpu,calorically_perfect,tol_iter_nr)
    integer, intent(in) :: ilat,nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect
    real(rkind), intent(in) :: ptot0, ttot0, rgas0, aoa, t0, tol_iter_nr
    real(rkind), dimension(1-ng:, 1-ng:, 1-ng:, 1:), intent(inout) :: w_gpu
    real(rkind), dimension(1-ng:, 1-ng:, 1-ng:, 1:), intent(in) :: w_aux_gpu
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in) :: cv_coeff_gpu, cp_coeff_gpu
    real(rkind) :: rho,rhou,rhov,rhow,rhoe,ri,uu,vv,ww,qq,ee,tt,pp,del,gamloc,rmaf,rml,vel_mod
    real(rkind) :: rmfac,sinangle,cosangle
    integer :: i,j,k,l,m,iercuda
    cosangle = cos(aoa)
    sinangle = sin(aoa)

    if (ilat==1) then
      !$omp target teams distribute parallel do collapse(2) has_device_addr(w_gpu,w_aux_gpu,cv_coeff_gpu,cp_coeff_gpu)
      do k = 1,nz
        do j = 1,ny
          rho = w_gpu(1,j,k,1)
          rhou = w_gpu(1,j,k,2)
          rhov = w_gpu(1,j,k,3)
          rhow = w_gpu(1,j,k,4)
          rhoe = w_gpu(1,j,k,5)
          ri = 1._rkind/rho
          uu = rhou*ri
          vv = rhov*ri
          ww = rhow*ri
          qq = 0.5_rkind*(uu*uu+vv*vv+ww*ww)
          ee = rhoe/rho-qq
          tt = get_temperature_from_e_dev(ee,w_aux_gpu(1,j,k,6),t0,cv_coeff_gpu,indx_cp_l,indx_cp_r,&
          & calorically_perfect, tol_iter_nr)
          pp = rho*tt*rgas0
          gamloc = get_gamloc_dev(indx_cp_l,indx_cp_r,tt,t0,cp_coeff_gpu,calorically_perfect,rgas0)
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
          ee = get_e_from_temperature_dev(tt,t0,indx_cp_l,indx_cp_r,cv_coeff_gpu,calorically_perfect)
          do l=1,ng
            w_gpu(1-l,j,k,1) = rho
            w_gpu(1-l,j,k,2) = rhou
            w_gpu(1-l,j,k,3) = rhov
            w_gpu(1-l,j,k,4) = rhow
            w_gpu(1-l,j,k,5) = rho*ee+0.5_rkind*(rhou*rhou+rhov*rhov+rhow*rhow)/rho
          enddo
        enddo
      enddo
    elseif (ilat==4) then
      !$omp target teams distribute parallel do collapse(2) has_device_addr(w_gpu,w_aux_gpu,cv_coeff_gpu,cp_coeff_gpu)
      do k = 1,nz
        do i = 1,nx
          rho = w_gpu(i,ny,k,1)
          rhou = w_gpu(i,ny,k,2)
          rhov = w_gpu(i,ny,k,3)
          rhow = w_gpu(i,ny,k,4)
          rhoe = w_gpu(i,ny,k,5)
          ri = 1._rkind/rho
          uu = rhou*ri
          vv = rhov*ri
          ww = rhow*ri
          qq = 0.5_rkind*(uu*uu+vv*vv+ww*ww)
          ee = rhoe/rho-qq
          tt = get_temperature_from_e_dev(ee,w_aux_gpu(i,ny,k,6),t0,cv_coeff_gpu,indx_cp_l,&
          &indx_cp_r, calorically_perfect,tol_iter_nr)
          pp = rho*tt*rgas0
          gamloc = get_gamloc_dev(indx_cp_l,indx_cp_r,tt,t0,cp_coeff_gpu,calorically_perfect,rgas0)
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
          ee = get_e_from_temperature_dev(tt,t0,indx_cp_l,indx_cp_r,cv_coeff_gpu,calorically_perfect)
          do l=1,ng
            w_gpu(i,ny+l,k,1) = rho
            w_gpu(i,ny+l,k,2) = rhou
            w_gpu(i,ny+l,k,3) = rhov
            w_gpu(i,ny+l,k,4) = rhow
            w_gpu(i,ny+l,k,5) = rho*ee+0.5_rkind*(rhou*rhou+rhov*rhov+rhow*rhow)/rho
          enddo
        enddo
      enddo
    endif

  endsubroutine bcfree_sub_kernel



  subroutine bcshock_kernel(ilat,nx,ny,nz,ng,nv,w_gpu,winf_gpu,winf_past_shock_gpu,xshock_imp,shock_angle,x_gpu,y_gpu,tanhfacs)
    integer, intent(in) :: nx,ny,nz,ng,nv,ilat
    real(rkind), intent(in) :: xshock_imp, shock_angle, tanhfacs
    real(rkind), dimension(1:), intent(in) :: winf_gpu, winf_past_shock_gpu
    real(rkind), dimension(1-ng:, 1-ng:, 1-ng:, 1:), intent(inout) :: w_gpu
    real(rkind), dimension(1-ng:), intent(in) :: x_gpu, y_gpu
    integer :: i,k,l,m,iercuda
    real(rkind) :: xsh, xx, tanhlen, tanhf, dwinf
    tanhlen = 8._rkind*tanhfacs

    if (ilat==1) then
    elseif (ilat==2) then
    elseif (ilat==3) then
    elseif (ilat==4) then
      !$omp target teams distribute parallel do collapse(3) has_device_addr(winf_gpu,winf_past_shock_gpu,w_gpu,x_gpu,y_gpu)
      do k = 1,nz
        do i = 1,nx
          do l = 0,ng
            xsh = xshock_imp-y_gpu(ny+l)/tan(shock_angle)
            xx = x_gpu(i)-xsh
            if (abs(xx)>tanhlen) then
              do m=1,nv
                w_gpu(i,ny+l,k,m) = w_gpu(i,ny,k,m)
              enddo
            else
              tanhf = 0.5_rkind*(1._rkind+tanh(xx/tanhfacs))
              do m=1,nv
                dwinf = winf_past_shock_gpu(m)-winf_gpu(m)
                w_gpu(i,ny+l,k,m) = winf_gpu(m)+dwinf*tanhf
              enddo
            endif
          enddo
        enddo
      enddo
    elseif (ilat==5) then
    elseif (ilat==6) then
    endif

  endsubroutine bcshock_kernel


  subroutine bcextr_var_kernel(nx,ny,nz,ng,w_gpu)
    integer :: nx,ny,nz,ng
    integer :: ilat
    integer :: i,j,k,l,m,iercuda
    real(rkind), dimension(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng, 1:1), intent(inout) :: w_gpu

    !$omp target teams distribute parallel do collapse(2) has_device_addr(w_gpu)
    do k = 1,nz
      do j = 1,ny
        do l=1,ng
          w_gpu(1-l,j,k,1) = w_gpu(1,j,k,1)
        enddo
      enddo
    enddo

    !$omp target teams distribute parallel do collapse(2) has_device_addr(w_gpu)
    do k = 1,nz
      do j = 1,ny
        do l=1,ng
          w_gpu(nx+l,j,k,1) = w_gpu(nx,j,k,1)
        enddo
      enddo
    enddo

    !$omp target teams distribute parallel do collapse(2) has_device_addr(w_gpu)
    do k = 1,nz
      do i = 1,nx
        do l=1,ng
          w_gpu(i,1-l,k,1) = w_gpu(i,1,k,1)
        enddo
      enddo
    enddo

    !$omp target teams distribute parallel do collapse(2) has_device_addr(w_gpu)
    do k = 1,nz
      do i = 1,nx
        do l=1,ng
          w_gpu(i,ny+l,k,1) = w_gpu(i,ny,k,1)
        enddo
      enddo
    enddo

    !$omp target teams distribute parallel do collapse(2) has_device_addr(w_gpu)
    do j = 1,ny
      do i = 1,nx
        do l=1,ng
          w_gpu(i,j,1-l,1) = w_gpu(i,j,1,1)
        enddo
      enddo
    enddo

    !$omp target teams distribute parallel do collapse(2) has_device_addr(w_gpu)
    do j = 1,ny
      do i = 1,nx
        do l=1,ng
          w_gpu(i,j,nz+l,1) = w_gpu(i,j,nz,1)
        enddo
      enddo
    enddo


  endsubroutine bcextr_var_kernel



  subroutine bcextr_airfoil_var_kernel(nx,ny,nz,ng,ndim,wall_tag_gpu,w_gpu,ileftx,irightx,ileftz,irightz)
    integer :: nx,ny,nz,ng,ndim, ileftx, irightx, ileftz, irightz
    integer :: ilat
    integer :: i,j,k,l,m,iercuda
    real(rkind), dimension(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng, 1:1), intent(inout) :: w_gpu
    integer, dimension(1-ng:nx+ng), intent(in) :: wall_tag_gpu
    if (ileftx < 0) then
      !$omp target teams distribute parallel do collapse(2) has_device_addr(w_gpu,wall_tag_gpu)
      do k = 1,nz
        do j = 1,ny
          do l=1,ng
            w_gpu(1-l,j,k,1) = 2._rkind*w_gpu(2-l,j,k,1)-w_gpu(3-l,j,k,1)
          enddo
        enddo
      enddo
    endif
    if (irightx < 0) then
      !$omp target teams distribute parallel do collapse(2) has_device_addr(w_gpu,wall_tag_gpu)
      do k = 1,nz
        do j = 1,ny
          do l=1,ng
            w_gpu(nx+l,j,k,1) = 2._rkind*w_gpu(nx+l-1,j,k,1)-w_gpu(nx+l-2,j,k,1)
          enddo
        enddo
      enddo
    endif
    !$omp target teams distribute parallel do collapse(2) has_device_addr(w_gpu,wall_tag_gpu)
    do k = 1,nz
      do i = 1,nx
        do l=1,ng
          if(wall_tag_gpu(i) < 1) then
            w_gpu(i,1-l,k,1) = 2._rkind*w_gpu(i,2-l,k,1)-w_gpu(i,3-l,k,1)
          endif
        enddo
      enddo
    enddo

    !$omp target teams distribute parallel do collapse(2) has_device_addr(w_gpu,wall_tag_gpu)
    do k = 1,nz
      do i = 1,nx
        do l=1,ng
          w_gpu(i,ny+l,k,1) = 2._rkind*w_gpu(i,ny+l-1,k,1)-w_gpu(i,ny+l-2,k,1)
        enddo
      enddo
    enddo
    if(ndim == 3) then
      if (ileftz < 0) then
        !$omp target teams distribute parallel do collapse(2) has_device_addr(w_gpu,wall_tag_gpu)
        do j = 1,ny
          do i = 1,nx
            do l=1,ng
              w_gpu(i,j,1-l,1) = 2._rkind*w_gpu(i,j,2-l,1)-w_gpu(i,j,3-l,1)
            enddo
          enddo
        enddo
      endif
      if (irightz < 0) then
        !$omp target teams distribute parallel do collapse(2) has_device_addr(w_gpu,wall_tag_gpu)
        do j = 1,ny
          do i = 1,nx
            do l=1,ng
              w_gpu(i,j,nz+l,1) = 2._rkind*w_gpu(i,j,nz+l-1,1)-w_gpu(i,j,nz+l-2,1)
            enddo
          enddo
        enddo
      endif
    endif

  endsubroutine bcextr_airfoil_var_kernel


  subroutine bcextr_kernel(ilat,nx,ny,nz,ng,nv,w_gpu)
    integer :: nx,ny,nz,ng,nv
    integer :: ilat
    integer :: i,j,k,l,m,iercuda
    real(rkind), dimension(1-ng:, 1-ng:, 1-ng:, 1:), intent(inout) :: w_gpu
    if (ilat==1) then
      !$omp target teams distribute parallel do collapse(2) has_device_addr(w_gpu)
      do k = 1,nz
        do j = 1,ny
          do l=1,ng
            do m=1,nv
              w_gpu(1-l,j,k,m) = w_gpu(1,j,k,m)
            enddo
          enddo
        enddo
      enddo
    elseif (ilat==2) then
      !$omp target teams distribute parallel do collapse(2) has_device_addr(w_gpu)
      do k = 1,nz
        do j = 1,ny
          do l=1,ng
            do m=1,nv
              w_gpu(nx+l,j,k,m) = w_gpu(nx,j,k,m)
            enddo
          enddo
        enddo
      enddo
    elseif (ilat==3) then
      !$omp target teams distribute parallel do collapse(2) has_device_addr(w_gpu)
      do k = 1,nz
        do i = 1,nx
          do l=1,ng
            do m=1,nv
              w_gpu(i,1-l,k,m) = w_gpu(i,1,k,m)
            enddo
          enddo
        enddo
      enddo
    elseif (ilat==4) then
      !$omp target teams distribute parallel do collapse(2) has_device_addr(w_gpu)
      do k = 1,nz
        do i = 1,nx
          do l=1,ng
            do m=1,nv
              w_gpu(i,ny+l,k,m) = w_gpu(i,ny,k,m)
            enddo
          enddo
        enddo
      enddo
    elseif (ilat==5) then
      !$omp target teams distribute parallel do collapse(2) has_device_addr(w_gpu)
      do j = 1,ny
        do i = 1,nx
          do l=1,ng
            do m=1,nv
              w_gpu(i,j,1-l,m) = w_gpu(i,j,1,m)
            enddo
          enddo
        enddo
      enddo
    elseif (ilat==6) then
      !$omp target teams distribute parallel do collapse(2) has_device_addr(w_gpu)
      do j = 1,ny
        do i = 1,nx
          do l=1,ng
            do m=1,nv
              w_gpu(i,j,nz+l,m) = w_gpu(i,j,nz,m)
            enddo
          enddo
        enddo
      enddo
    endif

  endsubroutine bcextr_kernel



  subroutine bcsym_kernel(ilat,nx,ny,nz,ng,w_gpu)
    integer :: nx,ny,nz,ng
    integer :: ilat
    integer :: i,k,l, iercuda
    real(rkind), dimension(1-ng:, 1-ng:, 1-ng:, 1:), intent(inout) :: w_gpu
    if (ilat==1) then
    elseif (ilat==2) then
    elseif (ilat==3) then
      !$omp target teams distribute parallel do collapse(2) has_device_addr(w_gpu)
      do k = 1,nz
        do i = 1,nx
          do l=1,ng
            w_gpu(i,1-l,k,1) = w_gpu(i,1+l,k,1)
            w_gpu(i,1-l,k,2) = w_gpu(i,1+l,k,2)
            w_gpu(i,1-l,k,3) = w_gpu(i,1+l,k,3)
            w_gpu(i,1-l,k,4) = w_gpu(i,1+l,k,4)
            w_gpu(i,1-l,k,5) = w_gpu(i,1+l,k,5)
          enddo
        enddo
      enddo
    elseif (ilat==4) then
    elseif (ilat==5) then
    elseif (ilat==6) then
    endif

  endsubroutine bcsym_kernel



  subroutine bcsym_c2_kernel(ilat,nx,ny,nz,ng,w_gpu,dxdcsic2_gpu,dydcsic2_gpu)
    integer :: nx,ny,nz,ng
    integer :: ilat
    integer :: i,k,l, iercuda
    real(rkind), dimension(1-ng:, 1-ng:, 1-ng:, 1:), intent(inout) :: w_gpu
    real(rkind), dimension(1-ng:, 1-ng:), intent(in) :: dxdcsic2_gpu
    real(rkind), dimension(1-ng:, 1-ng:), intent(in) :: dydcsic2_gpu
    real(rkind) :: abcsym_i31,abcsym_i32
    real(rkind) :: abcsym_i41,abcsym_i42
    real(rkind) :: tauvers_i31,tauvers_i32
    real(rkind) :: tauvers_i41,tauvers_i42
    real(rkind) :: taumod
    if (ilat==1) then
    elseif (ilat==2) then
    elseif (ilat==3) then
      !$omp target teams distribute parallel do collapse(2) has_device_addr(w_gpu,dxdcsic2_gpu,dydcsic2_gpu)
      do k = 1,nz
        do i = 1,nx
          do l=1,ng
            taumod = sqrt(dxdcsic2_gpu(i,1)**2+dydcsic2_gpu(i,1)**2)
            tauvers_i31 = dxdcsic2_gpu(i,1)/taumod
            tauvers_i32 = dydcsic2_gpu(i,1)/taumod
            abcsym_i31 = tauvers_i31**2 - tauvers_i32**2
            abcsym_i32 = 2._rkind*tauvers_i31*tauvers_i32
            w_gpu(i,1-l,k,1) = w_gpu(i,1+l,k,1)
            w_gpu(i,1-l,k,4) = w_gpu(i,1+l,k,4)
            w_gpu(i,1-l,k,5) = w_gpu(i,1+l,k,5)
            w_gpu(i,1-l,k,2) = abcsym_i31*w_gpu(i,1+l,k,2) + abcsym_i32*w_gpu(i,1+l,k,3)
            w_gpu(i,1-l,k,3) = abcsym_i32*w_gpu(i,1+l,k,2) - abcsym_i31*w_gpu(i,1+l,k,3)
          enddo
        enddo
      enddo
    elseif (ilat==4) then
      !$omp target teams distribute parallel do collapse(2) has_device_addr(w_gpu,dxdcsic2_gpu,dydcsic2_gpu)
      do k = 1,nz
        do i = 1,nx
          do l=1,ng
            taumod = sqrt(dxdcsic2_gpu(i,ny)**2+dydcsic2_gpu(i,ny)**2)
            tauvers_i41 = -dxdcsic2_gpu(i,ny)/taumod
            tauvers_i42 = -dydcsic2_gpu(i,ny)/taumod
            abcsym_i41 = tauvers_i41**2 - tauvers_i42**2
            abcsym_i42 = 2._rkind*tauvers_i41*tauvers_i42
            w_gpu(i,ny+l,k,1) = w_gpu(i,ny-l,k,1)
            w_gpu(i,ny+l,k,4) = w_gpu(i,ny-l,k,4)
            w_gpu(i,ny+l,k,5) = w_gpu(i,ny-l,k,5)
            w_gpu(i,ny+l,k,2) = abcsym_i41*w_gpu(i,ny-l,k,2) + abcsym_i42*w_gpu(i,ny-l,k,3)
            w_gpu(i,ny+l,k,3) = abcsym_i42*w_gpu(i,ny-l,k,2) - abcsym_i41*w_gpu(i,ny-l,k,3)
          enddo
        enddo
      enddo
    elseif (ilat==5) then
    elseif (ilat==6) then
    endif

  endsubroutine bcsym_c2_kernel



  subroutine bcwall_kernel(ilat,nx,ny,nz,ng,twall,w_gpu,w_aux_gpu,indx_cp_l,indx_cp_r,cv_coeff_gpu,&
  &t0,rgas0,calorically_perfect,tol_iter_nr)
    integer :: nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect
    integer :: ilat
    real(rkind) :: twall, t0, rgas0, tol_iter_nr
    integer :: i,k,l, iercuda
    real(rkind) :: rho,uu,vv,ww,qq,pp,tt,rhoe,ee
    real(rkind), dimension(1-ng:, 1-ng:, 1-ng:, 1:), intent(inout) :: w_gpu
    real(rkind), dimension(1-ng:, 1-ng:, 1-ng:, 1:), intent(in) :: w_aux_gpu
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in) :: cv_coeff_gpu
    if (ilat==1) then
    elseif (ilat==2) then
    elseif (ilat==3) then
      !$omp target teams distribute parallel do collapse(2) has_device_addr(w_gpu,w_aux_gpu,cv_coeff_gpu)
      do k = 1,nz
        do i = 1,nx
          w_gpu(i,1,k,2) = 0._rkind
          w_gpu(i,1,k,3) = 0._rkind
          w_gpu(i,1,k,4) = 0._rkind
          ee = get_e_from_temperature_dev(twall, t0, indx_cp_l, indx_cp_r, cv_coeff_gpu, calorically_perfect)
          w_gpu(i,1,k,5) = w_gpu(i,1,k,1)*ee
          do l=1,ng
            rho = w_gpu(i,1+l,k,1)
            uu = w_gpu(i,1+l,k,2)/rho
            vv = w_gpu(i,1+l,k,3)/rho
            ww = w_gpu(i,1+l,k,4)/rho
            rhoe = w_gpu(i,1+l,k,5)
            qq = 0.5_rkind*(uu*uu+vv*vv+ww*ww)
            ee = rhoe/rho-qq
            tt = get_temperature_from_e_dev(ee,w_aux_gpu(i,1+l,k,6),t0,cv_coeff_gpu,indx_cp_l,&
            &indx_cp_r,calorically_perfect,tol_iter_nr)
            pp = rho*tt*rgas0
            tt = 2._rkind*twall-tt
            ee = get_e_from_temperature_dev(tt, t0, indx_cp_l, indx_cp_r, cv_coeff_gpu,calorically_perfect)
            rho = pp/tt/rgas0
            w_gpu(i,1-l,k,1) = rho
            w_gpu(i,1-l,k,2) = -rho*uu
            w_gpu(i,1-l,k,3) = -rho*vv
            w_gpu(i,1-l,k,4) = -rho*ww
            w_gpu(i,1-l,k,5) = rho*(ee+qq)
          enddo
        enddo
      enddo
    elseif (ilat==4) then
    elseif (ilat==5) then
    elseif (ilat==6) then
    endif

  endsubroutine bcwall_kernel



  subroutine bcwall_airfoil_kernel(ilat,nx,ny,nz,ng,twall,w_gpu,w_aux_gpu,wall_tag_gpu,indx_cp_l,&
  &indx_cp_r,cv_coeff_gpu,t0,rgas0,calorically_perfect,tol_iter_nr,xc2_gpu,yc2_gpu,a_tw,v_bs,thic,&
  &kx_tw,om_tw,xtw1,xtw2,time,dxdetanc2_gpu,dydetanc2_gpu,u0)
    integer :: nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect
    integer :: ilat
    integer :: i,k,l, iercuda
    integer, dimension(1-ng:nx+ng), intent(in) :: wall_tag_gpu
    real(rkind) :: u0, twall, t0, rgas0, tol_iter_nr
    real(rkind) :: rho,uu,vv,ww,qq,pp,tt,rhoe,ee
    real(rkind), dimension(1-ng:, 1-ng:, 1-ng:, 1:), intent(inout) :: w_gpu
    real(rkind), dimension(1-ng:, 1-ng:, 1-ng:, 1:), intent(in) :: w_aux_gpu
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in) :: cv_coeff_gpu
    real(rkind), intent(in), dimension(1-ng:nx+ng,1-ng:ny+ng) :: xc2_gpu, yc2_gpu
    real(rkind), intent(in), dimension(1-ng:, 1-ng:) :: dxdetanc2_gpu, dydetanc2_gpu
    real(rkind), intent(in) :: a_tw, v_bs, thic, kx_tw, om_tw, time, xtw1, xtw2
    real(rkind) :: dx1, dx2, ftanh, tww, ug, vg, wg
    if (ilat==1) then
    elseif (ilat==2) then
    elseif (ilat==3) then
      !$omp target teams distribute parallel do collapse(2) has_device_addr(w_gpu,w_aux_gpu,cv_coeff_gpu,xc2_gpu,yc2_gpu,dxdetanc2_gpu,dydetanc2_gpu,wall_tag_gpu)
      do k = 1,nz
        do i = 1,nx
          if(wall_tag_gpu(i) < 1) then
            dx1 = (xc2_gpu(i,1)-xtw1)
            dx2 = (xc2_gpu(i,1)-xtw2)
            if(a_tw > 0._rkind .and. wall_tag_gpu(i) == 0) then
              ftanh = 0.5_rkind*(tanh(dx1/thic)-tanh(dx2/thic))
              tww = ftanh*(a_tw*sin(kx_tw*xc2_gpu(i,1)-om_tw*time))
              w_gpu(i,1,k,2) = 0._rkind
              w_gpu(i,1,k,3) = 0._rkind
              w_gpu(i,1,k,4) = w_gpu(i,1,k,1)*tww
              ee = get_e_from_temperature_dev(twall, t0, indx_cp_l, indx_cp_r, cv_coeff_gpu, calorically_perfect)
              w_gpu(i,1,k,5) = w_gpu(i,1,k,1)*ee+w_gpu(i,1,k,1)*0.5_rkind*(tww*tww)
            elseif(abs(v_bs) > 0._rkind .and. wall_tag_gpu(i) == 0) then
              ftanh = 0.5_rkind*(tanh(dx1/thic)-tanh(dx2/thic))
              tww = ftanh*v_bs*u0
              w_gpu(i,1,k,2) = w_gpu(i,1,k,1)*tww*dxdetanc2_gpu(i,1)
              w_gpu(i,1,k,3) = w_gpu(i,1,k,1)*tww*dydetanc2_gpu(i,1)
              w_gpu(i,1,k,4) = 0._rkind
              ee = get_e_from_temperature_dev(twall, t0, indx_cp_l, indx_cp_r, cv_coeff_gpu, calorically_perfect)
              w_gpu(i,1,k,5) = w_gpu(i,1,k,1)*ee+w_gpu(i,1,k,1)*0.5_rkind*(tww*tww)
            else
              w_gpu(i,1,k,2) = 0._rkind
              w_gpu(i,1,k,3) = 0._rkind
              w_gpu(i,1,k,4) = 0._rkind
              ee = get_e_from_temperature_dev(twall, t0, indx_cp_l, indx_cp_r, cv_coeff_gpu, calorically_perfect)
              w_gpu(i,1,k,5) = w_gpu(i,1,k,1)*ee
            endif
            do l=1,ng
              rho = w_gpu(i,1+l,k,1)
              uu = w_gpu(i,1+l,k,2)/rho
              vv = w_gpu(i,1+l,k,3)/rho
              ww = w_gpu(i,1+l,k,4)/rho
              rhoe = w_gpu(i,1+l,k,5)
              qq = 0.5_rkind*(uu*uu+vv*vv+ww*ww)
              ee = rhoe/rho-qq
              tt = get_temperature_from_e_dev(ee,w_aux_gpu(i,1+l,k,6),t0,cv_coeff_gpu,indx_cp_l,&
              &indx_cp_r,calorically_perfect,tol_iter_nr)
              pp = rho*tt*rgas0
              tt = 2._rkind*twall-tt
              ee = get_e_from_temperature_dev(tt, t0, indx_cp_l, indx_cp_r, cv_coeff_gpu,calorically_perfect)
              rho = pp/tt/rgas0
              ug = 2._rkind*w_gpu(i,1,k,2)/w_gpu(i,1,k,1) - uu
              vg = 2._rkind*w_gpu(i,1,k,3)/w_gpu(i,1,k,1) - vv
              wg = 2._rkind*w_gpu(i,1,k,4)/w_gpu(i,1,k,1) - ww
              qq = 0.5_rkind*(ug*ug+vg*vg+wg*wg)
              w_gpu(i,1-l,k,1) = rho
              w_gpu(i,1-l,k,2) = rho*ug
              w_gpu(i,1-l,k,3) = rho*vg
              w_gpu(i,1-l,k,4) = rho*wg
              w_gpu(i,1-l,k,5) = rho*(ee+qq)
            enddo
          endif
        enddo
      enddo
    elseif (ilat==4) then
    elseif (ilat==5) then
    elseif (ilat==6) then
    endif

  endsubroutine bcwall_airfoil_kernel



  subroutine bcwall_staggered_kernel(ilat,nx,ny,nz,ng,twall,w_gpu,w_aux_gpu,indx_cp_l,indx_cp_r,&
  &cv_coeff_gpu,t0,rgas0,calorically_perfect,tol_iter_nr)
    integer :: nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect
    integer :: ilat
    real(rkind) :: twall, t0, rgas0, tol_iter_nr
    integer :: i,k,l, iercuda
    real(rkind) :: rho,uu,vv,ww,qq,pp,tt,rhoe,ee
    real(rkind), dimension(1-ng:, 1-ng:, 1-ng:, 1:), intent(inout) :: w_gpu
    real(rkind), dimension(1-ng:, 1-ng:, 1-ng:, 1:), intent(in) :: w_aux_gpu
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in) :: cv_coeff_gpu
    if (ilat==1) then
    elseif (ilat==2) then
    elseif (ilat==3) then
      !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu,w_aux_gpu,cv_coeff_gpu)
      do k = 1,nz
        do l = 1,ng
          do i = 1,nx
            rho = w_gpu(i,l,k,1)
            uu = w_gpu(i,l,k,2)/rho
            vv = w_gpu(i,l,k,3)/rho
            ww = w_gpu(i,l,k,4)/rho
            rhoe = w_gpu(i,l,k,5)
            qq = 0.5_rkind*(uu*uu+vv*vv+ww*ww)
            ee = rhoe/rho-qq
            tt = get_temperature_from_e_dev(ee,w_aux_gpu(i,l,k,6),t0,cv_coeff_gpu,indx_cp_l,&
            &indx_cp_r,calorically_perfect,tol_iter_nr)
            pp = rho*tt*rgas0
            tt = 2._rkind*twall-tt
            ee = get_e_from_temperature_dev(tt, t0, indx_cp_l, indx_cp_r, cv_coeff_gpu,calorically_perfect)
            rho = pp/tt/rgas0
            w_gpu(i,1-l,k,1) = rho
            w_gpu(i,1-l,k,2) = -rho*uu
            w_gpu(i,1-l,k,3) = -rho*vv
            w_gpu(i,1-l,k,4) = -rho*ww
            w_gpu(i,1-l,k,5) = rho*(ee+qq)
          enddo
        enddo
      enddo
    elseif (ilat==4) then
      !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu,w_aux_gpu,cv_coeff_gpu)
      do k = 1,nz
        do l = 1,ng
          do i = 1,nx
            rho = w_gpu(i,ny+1-l,k,1)
            uu = w_gpu(i,ny+1-l,k,2)/rho
            vv = w_gpu(i,ny+1-l,k,3)/rho
            ww = w_gpu(i,ny+1-l,k,4)/rho
            rhoe = w_gpu(i,ny+1-l,k,5)
            qq = 0.5_rkind*(uu*uu+vv*vv+ww*ww)
            ee = rhoe/rho-qq
            tt = get_temperature_from_e_dev(ee,w_aux_gpu(i,ny+1-l,k,6),t0,cv_coeff_gpu,indx_cp_l,&
            &indx_cp_r,calorically_perfect,tol_iter_nr)
            pp = rho*tt*rgas0
            tt = 2._rkind*twall-tt
            ee = get_e_from_temperature_dev(tt, t0, indx_cp_l, indx_cp_r, cv_coeff_gpu,calorically_perfect)
            rho = pp/tt/rgas0
            w_gpu(i,ny+l,k,1) = rho
            w_gpu(i,ny+l,k,2) = -rho*uu
            w_gpu(i,ny+l,k,3) = -rho*vv
            w_gpu(i,ny+l,k,4) = -rho*ww
            w_gpu(i,ny+l,k,5) = rho*(ee+qq)
          enddo
        enddo
      enddo
    elseif (ilat==5) then
    elseif (ilat==6) then
    endif

  endsubroutine bcwall_staggered_kernel



  subroutine compute_residual_kernel(nx,ny,nz,ng,nv,fln_gpu,dt,residual_rhou,fluid_mask_gpu)
    integer :: nx, ny, nz, ng, nv
    real(rkind), intent(out) :: residual_rhou
    real(rkind), intent(in) :: dt
    real(rkind), dimension(1:nx, 1:ny, 1:nz, nv), intent(in) :: fln_gpu
    integer, dimension(1-ng:,1-ng:,1-ng:), intent(in) :: fluid_mask_gpu
    integer :: i,j,k,iercuda
    residual_rhou = 0._rkind
    !$omp target teams distribute parallel do collapse(3) has_device_addr(fln_gpu,fluid_mask_gpu) &
    !$omp& reduction(+:residual_rhou)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          if (fluid_mask_gpu(i,j,k)==0) then
            residual_rhou = residual_rhou + (fln_gpu(i,j,k,2)/dt)**2
          endif
        enddo
      enddo
    enddo


  endsubroutine compute_residual_kernel



  subroutine compute_airfoil_forces_runtime_kernel(nx,ny,nz,ng,nv,p0,u0,rgas0,w_aux_gpu,meta_gpu,&
  &csimod_gpu,wall_tag_gpu,dxdcsic2_gpu,dydcsic2_gpu,n,a,pn,pa,tn,ta)
    integer :: nx, ny, nz, ng, nv
    real(rkind) :: p0, u0, rgas0
    real(rkind) :: n, a, pn, pa, tn, ta
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:) :: w_aux_gpu
    real(rkind), dimension(1-ng:,1-ng:), intent(in) :: meta_gpu, csimod_gpu
    integer, dimension(1-ng:nx+ng), intent(in) :: wall_tag_gpu
    real(rkind), dimension(1-ng:,1-ng:), intent(in) :: dxdcsic2_gpu, dydcsic2_gpu
    real(rkind) :: dudy, dudyw, tauw, pw, cf, cp, pwf, tauwf, ds, costh, sinth
    real(rkind) :: dn, da, dpn, dpa, dtn, dta, ut1, ut2, ut3, ut4
    integer :: i,j,k,iercuda
    n = 0._rkind
    a = 0._rkind
    pn = 0._rkind
    pa = 0._rkind
    tn = 0._rkind
    ta = 0._rkind
    !$omp target teams distribute parallel do collapse(3) has_device_addr(w_aux_gpu,meta_gpu,csimod_gpu,dxdcsic2_gpu,dydcsic2_gpu,wall_tag_gpu) &
    !$omp& reduction(+:n,a,pn,pa,tn,ta)
    do k = 1,nz
      do j = 1,1
        do i = 1,nx
          if(wall_tag_gpu(i) < 1) then
            ds = csimod_gpu(i,1)
            costh = dxdcsic2_gpu(i,1)/ds
            sinth = dydcsic2_gpu(i,1)/ds
            ut1 = w_aux_gpu(i,1,k,2)*costh+w_aux_gpu(i,1,k,3)*sinth
            ut2 = w_aux_gpu(i,2,k,2)*costh+w_aux_gpu(i,2,k,3)*sinth
            ut3 = w_aux_gpu(i,3,k,2)*costh+w_aux_gpu(i,3,k,3)*sinth
            ut4 = w_aux_gpu(i,4,k,2)*costh+w_aux_gpu(i,4,k,3)*sinth
            dudy = -22._rkind*ut1 + 36._rkind*ut2 - 18._rkind*ut3 + 4._rkind*ut4
            dudyw = dudy*meta_gpu(i,1)/12._rkind
            tauw = w_aux_gpu(i,1,k,7)*dudyw
            pw = w_aux_gpu(i,1,k,1)*w_aux_gpu(i,1,k,6)*rgas0
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
            if(wall_tag_gpu(i) < -1) ds = ds/2._rkind
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


  endsubroutine compute_airfoil_forces_runtime_kernel



  subroutine compute_rho_t_p_minmax_kernel(nx,ny,nz,ng,rgas0,w_aux_gpu,rhomin,rhomax,tmin,tmax,pmin,pmax,fluid_mask_gpu)
    integer, intent(in) :: nx, ny, nz, ng
    real(rkind), intent(in) :: rgas0
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(in) :: w_aux_gpu
    integer, dimension(1-ng:,1-ng:,1-ng:), intent(in) :: fluid_mask_gpu
    real(rkind), intent(out) :: rhomin, rhomax, tmin, tmax, pmin, pmax
    integer :: i,j,k,iercuda
    real(rkind) :: rho,tt,pp
    rhomin = huge(1._rkind)
    rhomax = -100._rkind
    tmin = huge(1._rkind)
    tmax = -100._rkind
    pmin = huge(1._rkind)
    pmax = -100._rkind
    !$omp target teams distribute parallel do collapse(3) has_device_addr(w_aux_gpu,fluid_mask_gpu) &
    !$omp& reduction(min:rhomin,tmin,pmin)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          if (fluid_mask_gpu(i,j,k)==0) then
            rho = w_aux_gpu(i,j,k,1)
            tt = w_aux_gpu(i,j,k,6)
            pp = rho*tt*rgas0
            rhomin = min(rhomin,rho)
            tmin = min(tmin ,tt )
            pmin = min(pmin ,pp )
          endif
        enddo
      enddo
    enddo

    !$omp target teams distribute parallel do collapse(3) has_device_addr(w_aux_gpu,fluid_mask_gpu) &
    !$omp& reduction(max:rhomax,tmax,pmax)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          if (fluid_mask_gpu(i,j,k)==0) then
            rho = w_aux_gpu(i,j,k,1)
            tt = w_aux_gpu(i,j,k,6)
            pp = rho*tt*rgas0
            rhomax = max(rhomax,rho)
            tmax = max(tmax ,tt )
            pmax = max(pmax ,pp )
          endif
        enddo
      enddo
    enddo


  endsubroutine compute_rho_t_p_minmax_kernel



  subroutine compute_dt_kernel(nx,ny,nz,ng,rgas0,prandtl,dcsidx_gpu,detady_gpu,dzitdz_gpu,&
  &dcsidxs_gpu,detadys_gpu,dzitdzs_gpu,w_gpu,w_aux_gpu,dtxi_max,dtyi_max,dtzi_max,dtxv_max,dtyv_max,&
  &dtzv_max,dtxk_max,dtyk_max,dtzk_max,indx_cp_l,indx_cp_r,cp_coeff_gpu,fluid_mask_gpu,&
  &calorically_perfect,t0)
    integer :: nx, ny, nz, ng, indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind) :: rgas0, t0
    real(rkind) :: prandtl
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:) , intent(in) :: w_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(in) :: w_aux_gpu
    integer, dimension(1-ng:,1-ng:,1-ng:), intent(in) :: fluid_mask_gpu
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in) :: cp_coeff_gpu
    real(rkind), dimension(1:), intent(in) :: dcsidx_gpu, detady_gpu, dzitdz_gpu
    real(rkind), dimension(1:), intent(in) :: dcsidxs_gpu, detadys_gpu, dzitdzs_gpu
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
    !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu,w_aux_gpu,cp_coeff_gpu,dcsidx_gpu,detady_gpu,dzitdz_gpu,dcsidxs_gpu,detadys_gpu,dzitdzs_gpu,fluid_mask_gpu) &
    !$omp& reduction(max:dtxi_max,dtyi_max,dtzi_max,dtxv_max,dtyv_max,dtzv_max,dtxk_max,dtyk_max,dtzk_max)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          if (fluid_mask_gpu(i,j,k)==0) then
            rho = w_gpu(i,j,k,1)
            ri = 1._rkind/rho
            uu = w_aux_gpu(i,j,k,2)
            vv = w_aux_gpu(i,j,k,3)
            ww = w_aux_gpu(i,j,k,4)
            tt = w_aux_gpu(i,j,k,6)
            mu = w_aux_gpu(i,j,k,7)
            nu = ri*mu
            k_over_rhocp = nu/prandtl
            gamloc = get_gamloc_dev(indx_cp_l,indx_cp_r,tt,t0,cp_coeff_gpu,calorically_perfect,rgas0)
            c = sqrt (gamloc*rgas0*tt)
            dtxi = (abs(uu)+c)*dcsidx_gpu(i)
            dtyi = (abs(vv)+c)*detady_gpu(j)
            dtzi = (abs(ww)+c)*dzitdz_gpu(k)
            dtxv = nu*dcsidxs_gpu(i)
            dtyv = nu*detadys_gpu(j)
            dtzv = nu*dzitdzs_gpu(k)
            dtxk = k_over_rhocp*dcsidxs_gpu(i)
            dtyk = k_over_rhocp*detadys_gpu(j)
            dtzk = k_over_rhocp*dzitdzs_gpu(k)
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


  endsubroutine compute_dt_kernel



  subroutine compute_dt_c2_kernel(nx,ny,nz,ng,rgas0,prandtl,dcsidxnc2_gpu,dcsidync2_gpu,&
  &detadxnc2_gpu,detadync2_gpu,mcsi_gpu,meta_gpu,dzitdz_gpu,dzitdzs_gpu,w_gpu,w_aux_gpu,dtxi_max,&
  &dtyi_max,dtzi_max,dtxv_max,dtyv_max,dtzv_max,dtxk_max,dtyk_max,dtzk_max,indx_cp_l,indx_cp_r,&
  &cp_coeff_gpu,fluid_mask_gpu,calorically_perfect,t0)
    integer :: nx, ny, nz, ng, indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind) :: rgas0, t0
    real(rkind) :: prandtl
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:) , intent(in) :: w_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(in) :: w_aux_gpu
    integer, dimension(1-ng:,1-ng:,1-ng:), intent(in) :: fluid_mask_gpu
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in) :: cp_coeff_gpu
    real(rkind), dimension(1:), intent(in) :: dzitdz_gpu, dzitdzs_gpu
    real(rkind), dimension(1-ng:,1-ng:), intent(in) :: dcsidxnc2_gpu, dcsidync2_gpu, detadxnc2_gpu, detadync2_gpu
    real(rkind), dimension(1-ng:,1-ng:), intent(in) :: mcsi_gpu, meta_gpu
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
    !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu,w_aux_gpu,cp_coeff_gpu,dzitdz_gpu,dzitdzs_gpu,dcsidxnc2_gpu,dcsidync2_gpu,detadxnc2_gpu,detadync2_gpu,mcsi_gpu,meta_gpu,fluid_mask_gpu) &
    !$omp& reduction(max:dtxi_max,dtyi_max,dtzi_max,dtxv_max,dtyv_max,dtzv_max,dtxk_max,dtyk_max,dtzk_max)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          if (fluid_mask_gpu(i,j,k)==0) then
            rho = w_gpu(i,j,k,1)
            ri = 1._rkind/rho
            uu = w_aux_gpu(i,j,k,2)
            vv = w_aux_gpu(i,j,k,3)
            ww = w_aux_gpu(i,j,k,4)
            tt = w_aux_gpu(i,j,k,6)
            mu = w_aux_gpu(i,j,k,7)
            nu = ri*mu
            k_over_rhocp = nu/prandtl
            gamloc = get_gamloc_dev(indx_cp_l,indx_cp_r,tt,t0,cp_coeff_gpu,calorically_perfect,rgas0)
            k_over_rhocp = nu/prandtl*(gamloc)/(gamloc-1)
            c = sqrt (gamloc*rgas0*tt)
            csis1= mcsi_gpu(i,j)**2
            etas1= meta_gpu(i,j)**2
            util = uu*dcsidxnc2_gpu(i,j)+vv*dcsidync2_gpu(i,j)
            vtil = uu*detadxnc2_gpu(i,j)+vv*detadync2_gpu(i,j)
            dtxi = (abs(util)+c)*mcsi_gpu(i,j)
            dtyi = (abs(vtil)+c)*meta_gpu(i,j)
            dtzi = (abs(ww) +c)*dzitdz_gpu(k)
            dtxv = nu*csis1
            dtyv = nu*etas1
            dtzv = nu*dzitdzs_gpu(k)
            dtxk = k_over_rhocp*csis1
            dtyk = k_over_rhocp*etas1
            dtzk = k_over_rhocp*dzitdzs_gpu(k)
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


  endsubroutine compute_dt_c2_kernel



  subroutine eval_aux_kernel(nx,ny,nz,ng,istart,iend,jstart,jend,kstart,kend,w_gpu,w_aux_gpu,&
  &visc_model,mu0,t0,sutherland_s,t_ref_dim,powerlaw_vtexp,visc_power,visc_sutherland,visc_no,&
  &cv_coeff_gpu,indx_cp_l,indx_cp_r,rgas0,calorically_perfect,tol_iter_nr,stream_id)
    integer(ikind), intent(in) :: nx, ny, nz, ng, visc_model
    integer(ikind), intent(in) :: istart, iend, jstart, jend, kstart, kend
    integer(ikind), intent(in) :: visc_power, visc_sutherland, visc_no, indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind), intent(in) :: mu0, t0, sutherland_s, t_ref_dim, powerlaw_vtexp, rgas0, tol_iter_nr
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in) :: cv_coeff_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(in) :: w_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(inout) :: w_aux_gpu
    integer, intent(in) :: stream_id
    integer(ikind) :: i, j, k
    real(rkind) :: rho, rhou, rhov, rhow, rhoe, ri, uu, vv, ww, qq, pp, tt, mu, ee
    integer(ikind) :: iercuda

    !$omp target teams distribute parallel do collapse(3) has_device_addr(cv_coeff_gpu,w_gpu,w_aux_gpu)
    do j = jstart,jend
      do k = kstart,kend
        do i = istart,iend
          rho = w_gpu(i,j,k,1)
          rhou = w_gpu(i,j,k,2)
          rhov = w_gpu(i,j,k,3)
          rhow = w_gpu(i,j,k,4)
          rhoe = w_gpu(i,j,k,5)
          ri = 1._rkind/rho
          uu = rhou*ri
          vv = rhov*ri
          ww = rhow*ri
          qq = 0.5_rkind*(uu*uu+vv*vv+ww*ww)
          ee = rhoe/rho-qq
          tt = get_temperature_from_e_dev(ee, w_aux_gpu(i,j,k,6), t0, cv_coeff_gpu, indx_cp_l,&
          & indx_cp_r, calorically_perfect, tol_iter_nr)
          pp = rho*tt*rgas0
          w_aux_gpu(i,j,k,1) = rho
          w_aux_gpu(i,j,k,2) = uu
          w_aux_gpu(i,j,k,3) = vv
          w_aux_gpu(i,j,k,4) = ww
          w_aux_gpu(i,j,k,5) = (rhoe+pp)/rho
          w_aux_gpu(i,j,k,6) = tt
          if (visc_model == visc_power) then
            mu = mu0 * (tt/t0)**powerlaw_vtexp
          elseif (visc_model == visc_sutherland) then
            mu = mu0 * (tt/t0)**1.5_rkind * (1._rkind+sutherland_s/t_ref_dim)/(tt/t0 + sutherland_s/t_ref_dim)
          elseif (visc_model == visc_no) then
            mu = 0._rkind
          endif
          w_aux_gpu(i,j,k,7) = mu
        enddo
      enddo
    enddo


  endsubroutine eval_aux_kernel



  subroutine tripping_pressure_kernel(nx,ny,nz,ng,nv,i_rank_start,pi,itr1,itr2,x0tr,y0tr,x0ts,y0ts,&
  &lamx,lamy,lamz,lamz1,phiz,phiz1,asl,bt,xc2_gpu,yc2_gpu,z_gpu,w_gpu,fl_gpu,wall_tag_gpu,&
  &dxdetanc2_gpu,dydetanc2_gpu)
    real(rkind),intent(in) :: pi,x0tr,y0tr,x0ts,y0ts,lamx,lamy,lamz,lamz1,phiz,phiz1,asl,bt
    integer, intent(in) :: nx,ny,nz,nv,ng,itr1,itr2,i_rank_start
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in) :: w_gpu
    real(rkind), dimension(1:,1:,1:,1:), intent(inout) :: fl_gpu
    integer, dimension(1-ng:nx+ng), intent(in) :: wall_tag_gpu
    real(rkind), intent(in), dimension(1-ng:, 1-ng:) :: dxdetanc2_gpu, dydetanc2_gpu
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng), intent(in) :: xc2_gpu, yc2_gpu
    real(rkind), dimension(1-ng:), intent(in) :: z_gpu
    real(rkind) :: xx,yy,zz,hzi,hzi1,gzt,fz,fzx,fzy
    integer :: i,j,k,iercuda,ig

    !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu,fl_gpu,dxdetanc2_gpu,dydetanc2_gpu,xc2_gpu,yc2_gpu,z_gpu,wall_tag_gpu)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          ig = i + i_rank_start
          if (ig > itr1 .and. ig < itr2) then
            xx = xc2_gpu(i,j)
            yy = yc2_gpu(i,j)
            zz = z_gpu(k)
            hzi = sin(2._rkind*pi*lamz *zz+phiz )
            hzi1 = sin(2._rkind*pi*lamz1*zz+phiz1)
            gzt = asl*((1-bt)*hzi+bt*hzi1)
            fz = gzt*exp(-((xx-x0tr)/lamx)**2-((yy-y0tr)/lamy)**2)
            fzx = 0._rkind
            fzy = fz
            fl_gpu(i,j,k,2) = fl_gpu(i,j,k,2) - fzx*w_gpu(i,j,k,1)
            fl_gpu(i,j,k,3) = fl_gpu(i,j,k,3) - fzy*w_gpu(i,j,k,1)
            fl_gpu(i,j,k,5) = fl_gpu(i,j,k,5) -(fzx*w_gpu(i,j,k,2)+fzy*w_gpu(i,j,k,3))
          endif
        enddo
      enddo
    enddo


  endsubroutine tripping_pressure_kernel


  subroutine tripping_suction_kernel(nx,ny,nz,ng,nv,i_rank_start,pi,its1,its2,x0tr,y0tr,x0ts,y0ts,&
  &lamx,lamy,lams,lams1,phis,phis1,asl,bt,xc2_gpu,yc2_gpu,z_gpu,w_gpu,fl_gpu,wall_tag_gpu,&
  &dxdetanc2_gpu,dydetanc2_gpu)
    real(rkind),intent(in) :: pi,x0tr,y0tr,x0ts,y0ts,lamx,lamy,lams,lams1,phis,phis1,asl,bt
    integer, intent(in) :: nx,ny,nz,nv,ng,its1,its2,i_rank_start
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in) :: w_gpu
    real(rkind), dimension(1:,1:,1:,1:), intent(inout) :: fl_gpu
    integer, dimension(1-ng:nx+ng), intent(in) :: wall_tag_gpu
    real(rkind), intent(in), dimension(1-ng:, 1-ng:) :: dxdetanc2_gpu, dydetanc2_gpu
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng), intent(in) :: xc2_gpu, yc2_gpu
    real(rkind), dimension(1-ng:), intent(in) :: z_gpu
    real(rkind) :: xx,yy,zz,hzi,hzi1,gzt,fz,fzx,fzy
    integer :: i,j,k,iercuda,ig

    !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu,fl_gpu,dxdetanc2_gpu,dydetanc2_gpu,xc2_gpu,yc2_gpu,z_gpu,wall_tag_gpu)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          ig = i + i_rank_start
          if (ig > its1 .and. ig < its2) then
            xx = xc2_gpu(i,j)
            yy = yc2_gpu(i,j)
            zz = z_gpu(k)
            hzi = sin(2._rkind*pi*lams *zz+phis)
            hzi1 = sin(2._rkind*pi*lams1*zz+phis1)
            gzt = asl*((1-bt)*hzi+bt*hzi1)
            fz = gzt*exp(-((xx-x0ts)/lamx)**2-((yy-y0ts)/lamy)**2)
            fzx = 0._rkind
            fzy = fz
            fl_gpu(i,j,k,2) = fl_gpu(i,j,k,2) - fzx*w_gpu(i,j,k,1)
            fl_gpu(i,j,k,3) = fl_gpu(i,j,k,3) - fzy*w_gpu(i,j,k,1)
            fl_gpu(i,j,k,5) = fl_gpu(i,j,k,5) -(fzx*w_gpu(i,j,k,2)+fzy*w_gpu(i,j,k,3))
          endif
        enddo
      enddo
    enddo


  endsubroutine tripping_suction_kernel



  function get_gamloc_dev(indx_cp_l,indx_cp_r,tt,t0,cp_coeff_gpu,calorically_perfect,rgas0)
    real(rkind) :: get_gamloc_dev
    integer, value :: indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind), value :: tt, rgas0, t0
    real(rkind), dimension(indx_cp_l:indx_cp_r+1) :: cp_coeff_gpu
    real(rkind) :: cploc, gamloc
    integer :: l
    !$omp declare target
    if (calorically_perfect==1) then
      cploc = cp_coeff_gpu(0)
    else
      cploc = 0._rkind
      do l=indx_cp_l,indx_cp_r
        cploc = cploc+cp_coeff_gpu(l)*(tt/t0)**l
      enddo
    endif
    gamloc = cploc/(cploc-rgas0)
    get_gamloc_dev = gamloc
  endfunction get_gamloc_dev


  function get_temperature_from_e_dev(ee,t_start,t0,cv_coeff_gpu,indx_cp_l,indx_cp_r,calorically_perfect,tol_iter_nr)
    real(rkind) :: get_temperature_from_e_dev
    integer, value :: indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind), value :: ee, t_start, t0, tol_iter_nr
    real(rkind), dimension(indx_cp_l:indx_cp_r+1) :: cv_coeff_gpu
    real(rkind) :: tt, t_old, ebar, den, num, t_pow, t_powp
    integer :: l,iter,max_iter
    !$omp declare target
    max_iter = 50

    ebar = ee - cv_coeff_gpu(indx_cp_r+1)*t0
    if (calorically_perfect==1) then
      tt = t0+ebar/cv_coeff_gpu(0)
    else
      t_old = t_start
      do iter=1,max_iter
        den = 0._rkind
        num = 0._rkind
        do l=indx_cp_l,indx_cp_r
          if (l==-1) then
            t_pow = (t_old/t0)**l
            den = den+cv_coeff_gpu(l)*t_pow
            num = num+cv_coeff_gpu(l)*log(t_old/t0)
          else
            t_pow = (t_old/t0)**l
            t_powp = (t_old/t0)*t_pow
            den = den+cv_coeff_gpu(l)*t_pow
            num = num+cv_coeff_gpu(l)*(t_powp-1._rkind)/(l+1._rkind)
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


  function get_e_from_temperature_dev(tt,t0,indx_cp_l,indx_cp_r,cv_coeff_gpu,calorically_perfect)
    real(rkind) :: get_e_from_temperature_dev
    real(rkind),value :: tt, t0
    integer,value :: indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind), dimension(indx_cp_l:indx_cp_r+1) :: cv_coeff_gpu
    real(rkind) :: ee
    integer :: l
    !$omp declare target
    ee = cv_coeff_gpu(indx_cp_r+1)
    if (calorically_perfect==1) then
      ee = ee+cv_coeff_gpu(0)*(tt/t0-1._rkind)
    else
      do l=indx_cp_l,indx_cp_r
        if (l==-1) then
          ee = ee+cv_coeff_gpu(l)*log(tt/t0)
        else
          ee = ee+cv_coeff_gpu(l)/(l+1._rkind)*((tt/t0)**(l+1)-1._rkind)
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
    !$omp declare target
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
    !$omp declare target
    vkt = 1._rkind/2.12_rkind
    vkti = 2.12_rkind
    big_gam = 0.01_rkind*(yp*pr)**4
    big_gam = big_gam/(1._rkind+5._rkind*pr**3*yp)
    onethird = 1._rkind/3._rkind
    beta = vkti*log(pr)+(3.85_rkind*pr**onethird-1.3_rkind)**2

    tp = pr*yp*exp(-big_gam)+exp(-1._rkind/big_gam)*(beta+vkti*log(1+yp))
    tem_law_of_the_wall = tp
  endfunction tem_law_of_the_wall





  subroutine insitu_swirling_kernel(nv,nx,ny,nz,visc_order,ng,npsi,mpsi,dcsidx_gpu,detady_gpu,&
  &dzitdz_gpu,w_aux_gpu,coeff_deriv1_gpu,psi_gpu,u0,x_gpu)
    integer, value :: nv, nx, ny, nz, visc_order, ng, npsi,mpsi
    real(rkind), value :: u0
    real(rkind), dimension(nx) :: dcsidx_gpu
    real(rkind), dimension(ny) :: detady_gpu
    real(rkind), dimension(nz) :: dzitdz_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: psi_gpu
    real(rkind), dimension(1-ng:, 1-ng:,1-ng:,1:) :: w_aux_gpu
    real(rkind), dimension(1:,1:), intent(in) :: coeff_deriv1_gpu
    real(rkind), dimension(1-ng:nx+ng) :: x_gpu
    real(rkind) :: uu,vv,ww
    real(rkind) :: ux,uy,uz
    real(rkind) :: vx,vy,vz
    real(rkind) :: wx,wy,wz
    real(rkind), dimension(3) :: eigr_a, eigi_a
    real(rkind), dimension(3,3) :: astar
    real(rkind) :: ccl, div3l
    real(rkind) :: epsi2, omx, omy, omz, omod2, div2, div
    integer :: i,j,k,l
    !$omp target data map(alloc:eigr_a,eigi_a,astar)
    !$omp target teams distribute parallel do collapse(2) has_device_addr(dcsidx_gpu,detady_gpu,dzitdz_gpu,psi_gpu,w_aux_gpu,coeff_deriv1_gpu,x_gpu) private(eigr_a,eigi_a,astar)
    do k = 1,nz
      do i = 1,nx
        do j=1,ny
          uu = w_aux_gpu(i,j,k,2)
          vv = w_aux_gpu(i,j,k,3)
          ww = w_aux_gpu(i,j,k,4)
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
            ccl = coeff_deriv1_gpu(l,visc_order/2)
            ux = ux+ccl*(w_aux_gpu(i+l,j,k,2)-w_aux_gpu(i-l,j,k,2))
            vx = vx+ccl*(w_aux_gpu(i+l,j,k,3)-w_aux_gpu(i-l,j,k,3))
            wx = wx+ccl*(w_aux_gpu(i+l,j,k,4)-w_aux_gpu(i-l,j,k,4))
            uy = uy+ccl*(w_aux_gpu(i,j+l,k,2)-w_aux_gpu(i,j-l,k,2))
            vy = vy+ccl*(w_aux_gpu(i,j+l,k,3)-w_aux_gpu(i,j-l,k,3))
            wy = wy+ccl*(w_aux_gpu(i,j+l,k,4)-w_aux_gpu(i,j-l,k,4))
            uz = uz+ccl*(w_aux_gpu(i,j,k+l,2)-w_aux_gpu(i,j,k-l,2))
            vz = vz+ccl*(w_aux_gpu(i,j,k+l,3)-w_aux_gpu(i,j,k-l,3))
            wz = wz+ccl*(w_aux_gpu(i,j,k+l,4)-w_aux_gpu(i,j,k-l,4))
          enddo
          ux = ux*dcsidx_gpu(i)
          vx = vx*dcsidx_gpu(i)
          wx = wx*dcsidx_gpu(i)
          uy = uy*detady_gpu(j)
          vy = vy*detady_gpu(j)
          wy = wy*detady_gpu(j)
          uz = uz*dzitdz_gpu(k)
          vz = vz*dzitdz_gpu(k)
          wz = wz*dzitdz_gpu(k)
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
          psi_gpu(i,j,k,mpsi) = 2._rkind*max(0._rkind,eigi_a(2))
        enddo
      enddo
    enddo

    !$omp end target data
  endsubroutine insitu_swirling_kernel


  subroutine insitu_swirling_c2_kernel(nv,nx,ny,nz,visc_order,ng,npsi,mpsi,dzitdz_gpu,w_aux_gpu,&
  &coeff_deriv1_gpu,dcsidxc2_gpu,detadyc2_gpu,detadxc2_gpu,dcsidyc2_gpu,psi_gpu,u0,x_gpu,vis_tag_gpu,&
  &wall_tag_gpu)
    integer, value :: nv, nx, ny, nz, visc_order, ng, npsi,mpsi
    real(rkind), value :: u0
    real(rkind), dimension(nz) :: dzitdz_gpu
    real(rkind), dimension(1-ng:,1-ng:), intent(in) :: dcsidxc2_gpu, detadyc2_gpu
    real(rkind), dimension(1-ng:,1-ng:), intent(in) :: detadxc2_gpu, dcsidyc2_gpu
    integer, dimension(1-ng:), intent(in) :: vis_tag_gpu, wall_tag_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: psi_gpu
    real(rkind), dimension(1-ng:, 1-ng:,1-ng:,1:) :: w_aux_gpu
    real(rkind), dimension(1:,1:), intent(in) :: coeff_deriv1_gpu
    real(rkind), dimension(1-ng:nx+ng) :: x_gpu
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
    !$omp target data map(alloc:eigr_a,eigi_a,astar)
    !$omp target teams distribute parallel do collapse(2) has_device_addr(dzitdz_gpu,dcsidxc2_gpu,detadyc2_gpu,detadxc2_gpu,dcsidyc2_gpu,psi_gpu,w_aux_gpu,coeff_deriv1_gpu,x_gpu,vis_tag_gpu,wall_tag_gpu) private(eigr_a,eigi_a,astar)
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
          if (j == 1) lmaxi = vis_tag_gpu(i)
          if (wall_tag_gpu(i) < 1) lmaxj = min(j,visc_order/2)

          do l=1,visc_order/2
            cli = coeff_deriv1_gpu(l,lmaxi)
            clj = coeff_deriv1_gpu(l,lmaxj)
            clk = coeff_deriv1_gpu(l,visc_order/2 )

            ucsi = ucsi +cli*(w_aux_gpu(i+l,j,k,2)-w_aux_gpu(i-l,j,k,2))
            vcsi = vcsi +cli*(w_aux_gpu(i+l,j,k,3)-w_aux_gpu(i-l,j,k,3))
            wcsi = wcsi +cli*(w_aux_gpu(i+l,j,k,4)-w_aux_gpu(i-l,j,k,4))

            ueta = ueta +clj*(w_aux_gpu(i,j+l,k,2)-w_aux_gpu(i,j-l,k,2))
            veta = veta +clj*(w_aux_gpu(i,j+l,k,3)-w_aux_gpu(i,j-l,k,3))
            weta = weta +clj*(w_aux_gpu(i,j+l,k,4)-w_aux_gpu(i,j-l,k,4))

            uzit = uzit +clk*(w_aux_gpu(i,j,k+l,2)-w_aux_gpu(i,j,k-l,2))
            vzit = vzit +clk*(w_aux_gpu(i,j,k+l,3)-w_aux_gpu(i,j,k-l,3))
            wzit = wzit +clk*(w_aux_gpu(i,j,k+l,4)-w_aux_gpu(i,j,k-l,4))
          enddo

          ux = ucsi *dcsidxc2_gpu(i,j) + ueta *detadxc2_gpu(i,j)
          vx = vcsi *dcsidxc2_gpu(i,j) + veta *detadxc2_gpu(i,j)
          wx = wcsi *dcsidxc2_gpu(i,j) + weta *detadxc2_gpu(i,j)

          uy = ucsi *dcsidyc2_gpu(i,j) + ueta *detadyc2_gpu(i,j)
          vy = vcsi *dcsidyc2_gpu(i,j) + veta *detadyc2_gpu(i,j)
          wy = wcsi *dcsidyc2_gpu(i,j) + weta *detadyc2_gpu(i,j)
          uz = uzit *dzitdz_gpu(k)
          vz = vzit *dzitdz_gpu(k)
          wz = wzit *dzitdz_gpu(k)
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
          psi_gpu(i,j,k,mpsi) = 2._rkind*max(0._rkind,eigi_a(2))
        enddo
      enddo
    enddo

    !$omp end target data
  endsubroutine insitu_swirling_c2_kernel


  subroutine insitu_schlieren_kernel(nv,nx,ny,nz,visc_order,ng,npsi,mpsi,dcsidx_gpu,detady_gpu,&
  &dzitdz_gpu,w_aux_gpu,coeff_deriv1_gpu,psi_gpu,u0,x_gpu)
    integer, value :: nv, nx, ny, nz, visc_order, ng, npsi,mpsi
    real(rkind), value :: u0
    real(rkind), dimension(nx) :: dcsidx_gpu
    real(rkind), dimension(ny) :: detady_gpu
    real(rkind), dimension(nz) :: dzitdz_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: psi_gpu
    real(rkind), dimension(1-ng:, 1-ng:,1-ng:,1:) :: w_aux_gpu
    real(rkind), dimension(1:,1:), intent(in) :: coeff_deriv1_gpu
    real(rkind), dimension(1-ng:nx+ng) :: x_gpu
    real(rkind) :: rhox, rhoy, rhoz
    real(rkind) :: ccl
    integer :: i,j,k,l
    !$omp target teams distribute parallel do collapse(2) has_device_addr(dcsidx_gpu,detady_gpu,dzitdz_gpu,psi_gpu,w_aux_gpu,coeff_deriv1_gpu,x_gpu)
    do k = 1,nz
      do i = 1,nx
        do j=1,ny
          rhox = 0._rkind
          rhoy = 0._rkind
          rhoz = 0._rkind
          do l=1,visc_order/2
            ccl = coeff_deriv1_gpu(l,visc_order/2)
            rhox = rhox+ccl*(w_aux_gpu(i+l,j,k,1)-w_aux_gpu(i-l,j,k,1))
            rhoy = rhoy+ccl*(w_aux_gpu(i,j+l,k,1)-w_aux_gpu(i,j-l,k,1))
            rhoz = rhoz+ccl*(w_aux_gpu(i,j,k+l,1)-w_aux_gpu(i,j,k-l,1))
          enddo
          rhox = rhox*dcsidx_gpu(i)
          rhoy = rhoy*detady_gpu(j)
          rhoz = rhoz*dzitdz_gpu(k)
          psi_gpu(i,j,k,mpsi) = exp(-sqrt((rhox)**2+(rhoy)**2+(rhoz)**2))
        enddo
      enddo
    enddo

  endsubroutine insitu_schlieren_kernel


  subroutine insitu_schlieren_c2_kernel(nv,nx,ny,nz,visc_order,ng,npsi,mpsi,dzitdz_gpu,w_aux_gpu,&
  &coeff_deriv1_gpu,dcsidxc2_gpu,detadyc2_gpu,detadxc2_gpu,dcsidyc2_gpu,psi_gpu,u0,x_gpu,vis_tag_gpu,&
  &wall_tag_gpu)
    integer, value :: nv, nx, ny, nz, visc_order, ng, npsi,mpsi
    real(rkind), value :: u0
    real(rkind), dimension(nz) :: dzitdz_gpu
    real(rkind), dimension(1-ng:,1-ng:), intent(in) :: dcsidxc2_gpu, detadyc2_gpu
    real(rkind), dimension(1-ng:,1-ng:), intent(in) :: detadxc2_gpu, dcsidyc2_gpu
    integer, dimension(1-ng:), intent(in) :: vis_tag_gpu, wall_tag_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: psi_gpu
    real(rkind), dimension(1-ng:, 1-ng:,1-ng:,1:) :: w_aux_gpu
    real(rkind), dimension(1:,1:), intent(in) :: coeff_deriv1_gpu
    real(rkind), dimension(1-ng:nx+ng) :: x_gpu
    real(rkind) :: rhocsi,rhoeta,rhozit,rhox,rhoy,rhoz
    real(rkind) :: cli,clj,clk
    integer :: i,j,k,l,lmaxi,lmaxj
    !$omp target teams distribute parallel do collapse(2) has_device_addr(dzitdz_gpu,dcsidxc2_gpu,detadyc2_gpu,detadxc2_gpu,dcsidyc2_gpu,psi_gpu,w_aux_gpu,coeff_deriv1_gpu,x_gpu,vis_tag_gpu,wall_tag_gpu)
    do k = 1,nz
      do i = 1,nx
        do j=1,ny
          rhocsi = 0._rkind
          rhoeta = 0._rkind
          rhozit = 0._rkind

          lmaxi = visc_order/2
          lmaxj = visc_order/2
          if (j == 1) lmaxi = vis_tag_gpu(i)
          if (wall_tag_gpu(i) < 1) lmaxj = min(j,visc_order/2)

          do l=1,visc_order/2
            cli = coeff_deriv1_gpu(l,lmaxi)
            clj = coeff_deriv1_gpu(l,lmaxj)
            clk = coeff_deriv1_gpu(l,visc_order/2 )

            rhocsi = rhocsi +cli*(w_aux_gpu(i+l,j,k,1)-w_aux_gpu(i-l,j,k,1))
            rhoeta = rhoeta +clj*(w_aux_gpu(i,j+l,k,1)-w_aux_gpu(i,j-l,k,1))
            rhozit = rhozit +clk*(w_aux_gpu(i,j,k+l,1)-w_aux_gpu(i,j,k-l,1))
          enddo

          rhox = rhocsi *dcsidxc2_gpu(i,j) + rhoeta *detadxc2_gpu(i,j)
          rhoy = rhocsi *dcsidyc2_gpu(i,j) + rhoeta *detadyc2_gpu(i,j)
          rhoz = rhozit *dzitdz_gpu(k)
          psi_gpu(i,j,k,mpsi) = exp(-sqrt((rhox)**2+(rhoy)**2+(rhoz)**2))
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
    !$omp declare target
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


  subroutine insitu_div_kernel(nx,ny,nz,ng,visc_order,npsi,mpsi,w_aux_gpu,coeff_deriv1_gpu,&
  &dcsidx_gpu,detady_gpu,dzitdz_gpu,psi_gpu)
    integer, intent(in) :: nx, ny, nz, ng, visc_order, npsi, mpsi
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(inout) :: w_aux_gpu
    real(rkind), dimension(1:,1:), intent(in) :: coeff_deriv1_gpu
    real(rkind), dimension(1:), intent(in) :: dcsidx_gpu, detady_gpu, dzitdz_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(inout) :: psi_gpu
    real(rkind) :: ccl
    real(rkind) :: uu,vv,ww
    real(rkind) :: ux,vy,wz
    real(rkind) :: div
    integer :: lmax
    integer :: i,j,k,l,iercuda
    lmax = visc_order / 2
    !$omp target teams distribute parallel do collapse(3) has_device_addr(w_aux_gpu,coeff_deriv1_gpu,dcsidx_gpu,detady_gpu,dzitdz_gpu,psi_gpu)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          uu = w_aux_gpu(i,j,k,2)
          vv = w_aux_gpu(i,j,k,3)
          ww = w_aux_gpu(i,j,k,4)
          ux = 0._rkind
          vy = 0._rkind
          wz = 0._rkind
          do l=1,lmax
            ccl = coeff_deriv1_gpu(l,lmax)
            ux = ux+ccl*(w_aux_gpu(i+l,j,k,2)-w_aux_gpu(i-l,j,k,2))
            vy = vy+ccl*(w_aux_gpu(i,j+l,k,3)-w_aux_gpu(i,j-l,k,3))
            wz = wz+ccl*(w_aux_gpu(i,j,k+l,4)-w_aux_gpu(i,j,k-l,4))
          enddo
          ux = ux*dcsidx_gpu(i)
          vy = vy*detady_gpu(j)
          wz = wz*dzitdz_gpu(k)
          div = ux+vy+wz
          psi_gpu(i,j,k,mpsi) = div
        enddo
      enddo
    enddo


  endsubroutine insitu_div_kernel



  subroutine insitu_omega_kernel(nx,ny,nz,ng,visc_order,npsi,mpsi,w_aux_gpu,coeff_deriv1_gpu,&
  &dcsidx_gpu,detady_gpu,dzitdz_gpu,psi_gpu)
    integer, intent(in) :: nx, ny, nz, ng, visc_order, npsi, mpsi
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(inout) :: w_aux_gpu
    real(rkind), dimension(1:,1:), intent(in) :: coeff_deriv1_gpu
    real(rkind), dimension(1:), intent(in) :: dcsidx_gpu, detady_gpu, dzitdz_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(inout) :: psi_gpu
    integer :: i,j,k,l,iercuda
    real(rkind) :: ccl
    real(rkind) :: uu,vv,ww
    real(rkind) :: uy,uz
    real(rkind) :: vx,vz
    real(rkind) :: wx,wy
    real(rkind) :: omegax, omegay, omegaz, omod2
    integer :: lmax
    lmax = visc_order/2
    !$omp target teams distribute parallel do collapse(3) has_device_addr(w_aux_gpu,coeff_deriv1_gpu,dcsidx_gpu,detady_gpu,dzitdz_gpu,psi_gpu)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          uu = w_aux_gpu(i,j,k,2)
          vv = w_aux_gpu(i,j,k,3)
          ww = w_aux_gpu(i,j,k,4)
          vx = 0._rkind
          wx = 0._rkind
          uy = 0._rkind
          wy = 0._rkind
          uz = 0._rkind
          vz = 0._rkind
          do l=1,lmax
            ccl = coeff_deriv1_gpu(l,lmax)
            vx = vx+ccl*(w_aux_gpu(i+l,j,k,3)-w_aux_gpu(i-l,j,k,3))
            wx = wx+ccl*(w_aux_gpu(i+l,j,k,4)-w_aux_gpu(i-l,j,k,4))
            uy = uy+ccl*(w_aux_gpu(i,j+l,k,2)-w_aux_gpu(i,j-l,k,2))
            wy = wy+ccl*(w_aux_gpu(i,j+l,k,4)-w_aux_gpu(i,j-l,k,4))
            uz = uz+ccl*(w_aux_gpu(i,j,k+l,2)-w_aux_gpu(i,j,k-l,2))
            vz = vz+ccl*(w_aux_gpu(i,j,k+l,3)-w_aux_gpu(i,j,k-l,3))
          enddo
          vx = vx*dcsidx_gpu(i)
          wx = wx*dcsidx_gpu(i)
          uy = uy*detady_gpu(j)
          wy = wy*detady_gpu(j)
          uz = uz*dzitdz_gpu(k)
          vz = vz*dzitdz_gpu(k)
          omegax = wy-vz
          omegay = uz-wx
          omegaz = vx-uy
          omod2 = omegax*omegax+omegay*omegay+omegaz*omegaz
          psi_gpu(i,j,k,mpsi) = sqrt(omod2)
        enddo
      enddo
    enddo


  endsubroutine insitu_omega_kernel



  subroutine insitu_ducros_kernel(nx,ny,nz,ng,visc_order,npsi,mpsi,u0,l0,w_aux_gpu,coeff_deriv1_gpu,&
  &dcsidx_gpu,detady_gpu,dzitdz_gpu,psi_gpu)
    integer, intent(in) :: nx, ny, nz, ng, visc_order, npsi, mpsi
    real(rkind), intent(in) :: u0, l0
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(inout) :: w_aux_gpu
    real(rkind), dimension(1:,1:), intent(in) :: coeff_deriv1_gpu
    real(rkind), dimension(1:), intent(in) :: dcsidx_gpu, detady_gpu, dzitdz_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(inout) :: psi_gpu
    integer :: i,j,k,l,iercuda
    real(rkind) :: ccl
    real(rkind) :: uu,vv,ww
    real(rkind) :: ux,uy,uz
    real(rkind) :: vx,vy,vz
    real(rkind) :: wx,wy,wz
    real(rkind) :: div, omegax, omegay, omegaz, omod2
    integer :: lmax
    lmax = visc_order/2
    !$omp target teams distribute parallel do collapse(3) has_device_addr(w_aux_gpu,coeff_deriv1_gpu,dcsidx_gpu,detady_gpu,dzitdz_gpu,psi_gpu)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          uu = w_aux_gpu(i,j,k,2)
          vv = w_aux_gpu(i,j,k,3)
          ww = w_aux_gpu(i,j,k,4)
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
            ccl = coeff_deriv1_gpu(l,lmax)
            ux = ux+ccl*(w_aux_gpu(i+l,j,k,2)-w_aux_gpu(i-l,j,k,2))
            vx = vx+ccl*(w_aux_gpu(i+l,j,k,3)-w_aux_gpu(i-l,j,k,3))
            wx = wx+ccl*(w_aux_gpu(i+l,j,k,4)-w_aux_gpu(i-l,j,k,4))
            uy = uy+ccl*(w_aux_gpu(i,j+l,k,2)-w_aux_gpu(i,j-l,k,2))
            vy = vy+ccl*(w_aux_gpu(i,j+l,k,3)-w_aux_gpu(i,j-l,k,3))
            wy = wy+ccl*(w_aux_gpu(i,j+l,k,4)-w_aux_gpu(i,j-l,k,4))
            uz = uz+ccl*(w_aux_gpu(i,j,k+l,2)-w_aux_gpu(i,j,k-l,2))
            vz = vz+ccl*(w_aux_gpu(i,j,k+l,3)-w_aux_gpu(i,j,k-l,3))
            wz = wz+ccl*(w_aux_gpu(i,j,k+l,4)-w_aux_gpu(i,j,k-l,4))
          enddo
          ux = ux*dcsidx_gpu(i)
          vx = vx*dcsidx_gpu(i)
          wx = wx*dcsidx_gpu(i)
          uy = uy*detady_gpu(j)
          vy = vy*detady_gpu(j)
          wy = wy*detady_gpu(j)
          uz = uz*dzitdz_gpu(k)
          vz = vz*dzitdz_gpu(k)
          wz = wz*dzitdz_gpu(k)
          div = ux+vy+wz
          omegax = wy-vz
          omegay = uz-wx
          omegaz = vx-uy
          omod2 = omegax*omegax+omegay*omegay+omegaz*omegaz
          psi_gpu(i,j,k,mpsi) = max(-div/sqrt(omod2+div**2+(u0/l0)**2),0._rkind)
        enddo
      enddo
    enddo


  endsubroutine insitu_ducros_kernel





  subroutine probe_interpolation_kernel(num_probe,nx,ny,nz,ng,nv_aux,ijk_probe_gpu,w_aux_probe_gpu,w_aux_gpu,probe_coeff_gpu)
    integer, intent(in) :: num_probe,nx,ny,nz,ng,nv_aux
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv_aux), intent(in) :: w_aux_gpu
    real(rkind), dimension(6,num_probe), intent(inout) :: w_aux_probe_gpu
    real(rkind), dimension(2,2,2,num_probe), intent(in) :: probe_coeff_gpu
    integer, dimension(3,num_probe), intent(in) :: ijk_probe_gpu
    integer :: i,j,k,ii,jj,kk,l,iercuda
    real(rkind) :: w1,w2,w3,w4,w5,w6

    !$omp target teams distribute parallel do collapse(1) has_device_addr(w_aux_gpu,w_aux_probe_gpu,probe_coeff_gpu,ijk_probe_gpu)
    do l = 1,num_probe
      ii = ijk_probe_gpu(1,l)
      jj = ijk_probe_gpu(2,l)
      kk = ijk_probe_gpu(3,l)
      w1 = 0._rkind
      w2 = 0._rkind
      w3 = 0._rkind
      w4 = 0._rkind
      w5 = 0._rkind
      w6 = 0._rkind
      do k=1,2
        do j=1,2
          do i=1,2
            w1 = w1 + probe_coeff_gpu(i,j,k,l)*w_aux_gpu(i+ii-1,j+jj-1,k+kk-1,1)
            w2 = w2 + probe_coeff_gpu(i,j,k,l)*w_aux_gpu(i+ii-1,j+jj-1,k+kk-1,2)
            w3 = w3 + probe_coeff_gpu(i,j,k,l)*w_aux_gpu(i+ii-1,j+jj-1,k+kk-1,3)
            w4 = w4 + probe_coeff_gpu(i,j,k,l)*w_aux_gpu(i+ii-1,j+jj-1,k+kk-1,4)
            w5 = w5 + probe_coeff_gpu(i,j,k,l)*w_aux_gpu(i+ii-1,j+jj-1,k+kk-1,5)
            w6 = w6 + probe_coeff_gpu(i,j,k,l)*w_aux_gpu(i+ii-1,j+jj-1,k+kk-1,6)
          enddo
        enddo
      enddo
      w_aux_probe_gpu(1,l) = w1
      w_aux_probe_gpu(2,l) = w2
      w_aux_probe_gpu(3,l) = w3
      w_aux_probe_gpu(4,l) = w4
      w_aux_probe_gpu(5,l) = w5
      w_aux_probe_gpu(6,l) = w6
    enddo


  endsubroutine probe_interpolation_kernel



  subroutine compute_tspec_kernel(nx,ny,nz,ng,ndft,j_slice,i_win,it_win,w_aux_gpu,w_tspec_gpu)
    integer, intent(in) :: nx, ny, nz, ng, ndft, j_slice, i_win, it_win
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in) :: w_aux_gpu
    real(rkind), dimension(1:,1:,1:,1:,1:), intent(inout) :: w_tspec_gpu
    real(rkind) :: wei_real, wei_imag
    real(rkind) :: pp, win_scale
    integer :: i,j,k,iercuda, ll, l, nn, n
    real(rkind) :: pi
    pi = acos(-1._rkind)
    !$omp target teams distribute parallel do collapse(2) has_device_addr(w_aux_gpu,w_tspec_gpu)
    do k = 1,nz
      do i = 1,nx
        do l=1,ndft
          n=it_win
          ll = l-1
          nn = n-1
          wei_real = cos((-2.0_rkind*pi*ll*nn)/ndft)
          wei_imag = sin((-2.0_rkind*pi*ll*nn)/ndft)
          win_scale = 1._rkind
          pp = win_scale*w_aux_gpu(i,j_slice,k,1)*w_aux_gpu(i,j_slice,k,6)
          if(nn==0) then
            w_tspec_gpu(i,k,l,i_win,1) = pp*wei_real
            w_tspec_gpu(i,k,l,i_win,2) = pp*wei_imag
          else
            w_tspec_gpu(i,k,l,i_win,1) = w_tspec_gpu(i,k,l,i_win,1) + pp*wei_real
            w_tspec_gpu(i,k,l,i_win,2) = w_tspec_gpu(i,k,l,i_win,2) + pp*wei_imag
          endif
        enddo
      enddo
    enddo


  endsubroutine compute_tspec_kernel



  subroutine compute_psd_tspec_kernel(nx,ny,nz,ndft,i_win,dt_tspec,w_tspec_gpu,w_psd_tspec_gpu)
    integer, intent(in) :: nx, ny, nz, ndft, i_win
    real(rkind), intent(in) :: dt_tspec
    real(rkind), dimension(1:,1:,1:,1:,1:), intent(inout) :: w_tspec_gpu
    real(rkind), dimension(1:,1:), intent(inout) :: w_psd_tspec_gpu
    integer :: i,j,k,iercuda, ll, l, nn, n
    real(rkind) :: pi
    pi = acos(-1._rkind)
    !$omp target teams distribute parallel do collapse(2) has_device_addr(w_tspec_gpu,w_psd_tspec_gpu)
    do l = 1,ndft
      do i = 1,nx
        w_psd_tspec_gpu(i,l) = 0._rkind
        do k=1,nz
          w_psd_tspec_gpu(i,l) = w_psd_tspec_gpu(i,l)+(w_tspec_gpu(i,k,l,i_win,1)**2+w_tspec_gpu(i,k,l,i_win,2)**2)
        enddo
      enddo
    enddo


  endsubroutine compute_psd_tspec_kernel



  subroutine compute_wallprop_c2_kernel(nx,ny,nz,ng,w_aux_gpu,wallprop_gpu,dxdcsic2_gpu,&
  &dydcsic2_gpu,csimod_gpu,meta_gpu,wall_tag_gpu)
    integer, intent(in) :: nx,ny,nz,ng
    integer, dimension(1-ng:nx+ng), intent(in) :: wall_tag_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(in) :: w_aux_gpu
    real(rkind), dimension(1-ng:,1-ng:, 2:), intent(inout) :: wallprop_gpu
    real(rkind), dimension(1-ng:,1-ng:), intent(in) :: meta_gpu, csimod_gpu
    real(rkind), dimension(1-ng:,1-ng:), intent(in) :: dxdcsic2_gpu, dydcsic2_gpu
    integer :: i,j,k,iercuda
    real(rkind) :: ut1,ut2,ut3,ut4,ds,costh,sinth,dudy,dudyw

    !$omp target teams distribute parallel do collapse(2) has_device_addr(w_aux_gpu,wallprop_gpu,meta_gpu,csimod_gpu,dxdcsic2_gpu,dydcsic2_gpu,wall_tag_gpu)
    do k = 1,nz
      do i = 1,nx
        if(wall_tag_gpu(i) < 1) then
          ds = csimod_gpu(i,1)
          costh = dxdcsic2_gpu(i,1)/ds
          sinth = dydcsic2_gpu(i,1)/ds
          ut1 = w_aux_gpu(i,1,k,2)*costh+w_aux_gpu(i,1,k,3)*sinth
          ut2 = w_aux_gpu(i,2,k,2)*costh+w_aux_gpu(i,2,k,3)*sinth
          ut3 = w_aux_gpu(i,3,k,2)*costh+w_aux_gpu(i,3,k,3)*sinth
          ut4 = w_aux_gpu(i,4,k,2)*costh+w_aux_gpu(i,4,k,3)*sinth
          dudy = -22._rkind*ut1+36._rkind*ut2-18._rkind*ut3+4._rkind*ut4
          dudyw = dudy*meta_gpu(i,1)/12._rkind
          wallprop_gpu(i,k,2) = w_aux_gpu(i,1,k,7)*dudyw
        else
          wallprop_gpu(i,k,2) = w_aux_gpu(i,1,k,2)
        endif
      enddo
    enddo


  endsubroutine compute_wallprop_c2_kernel



endmodule streams_kernels_omp


