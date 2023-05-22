module streams_kernels_gpu
!
  use streams_parameters, only : rkind, ikind, REAL64
  use CUDAFOR
  implicit none
!
contains
!
  subroutine zero_flux_cuf(nx, ny, nz, nv, fl_gpu)
    integer :: nx, ny, nz, nv
    real(rkind), dimension(1:,1:,1:,1:), intent(inout), device :: fl_gpu
    integer :: i,j,k,m,iercuda
    !$cuf kernel do(3) <<<*,*>>>
    do k=1,nz
      do j=1,ny
        do i=1,nx
          do m=1,nv
            fl_gpu(i,j,k,m)  = 0._rkind
          enddo
        enddo
      enddo
    enddo
    !@cuf iercuda=cudaDeviceSynchronize()
  endsubroutine zero_flux_cuf
!
  subroutine init_flux_cuf(nx, ny, nz, nv, fl_gpu, fln_gpu, rhodt)
    integer :: nx, ny, nz, nv
    real(rkind) :: rhodt
    real(rkind), dimension(1:,1:,1:,1:), intent(inout), device :: fl_gpu, fln_gpu
    integer :: i,j,k,m,iercuda
    !$cuf kernel do(3) <<<*,*>>>
    do k=1,nz
      do j=1,ny
        do i=1,nx
          do m=1,nv
            fln_gpu(i,j,k,m) = - rhodt * fl_gpu(i,j,k,m)
            fl_gpu(i,j,k,m)  = 0._rkind
          enddo
        enddo
      enddo
    enddo
    !@cuf iercuda=cudaDeviceSynchronize()
  endsubroutine init_flux_cuf
!
  subroutine euler_x_transp_cuf(nx, ny, nz, ng, istart, iend, nv_aux, w_aux_gpu, w_aux_trans_gpu, stream_id)
    integer :: nx,ny,nz,ng,nv_aux
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv_aux), intent(in), device :: w_aux_gpu
    real(rkind), dimension(1-ng:ny+ng,1-ng:nx+ng,1-ng:nz+ng,1:8), intent(inout), device :: w_aux_trans_gpu
    integer, intent(in) :: istart, iend
    integer(kind=cuda_stream_kind), intent(in) :: stream_id
    integer :: i,j,k,iv,iercuda
    !$cuf kernel do(3) <<<*,*,stream=stream_id>>>
    do k=1,nz
      do i=istart,iend
        do j=1,ny
          do iv=1,8
            w_aux_trans_gpu(j,i,k,iv) = w_aux_gpu(i,j,k,iv)
          enddo
        enddo
      enddo
    enddo
  endsubroutine euler_x_transp_cuf
!
  attributes(global) launch_bounds(256) subroutine euler_x_fluxes_hybrid_kernel(nv, nv_aux, nx, ny, nz, ng, &
    eul_imin, eul_imax, lmax_base, nkeep, rgas0, coeff_deriv1_gpu, dcsidx_gpu, w_aux_trans_gpu, fhat_trans_gpu, &
    force_zero_flux_min, force_zero_flux_max, &
    weno_scheme, weno_version, sensor_threshold, weno_size, gplus_x_gpu, gminus_x_gpu, cp_coeff_gpu, &
    indx_cp_l, indx_cp_r, ep_ord_change_x_gpu, calorically_perfect, tol_iter_nr,rho0,u0,t0)
!
    implicit none
!   Passed arguments
    integer, value :: nv, nx, ny, nz, ng, nv_aux
    integer, value :: eul_imin, eul_imax, lmax_base, nkeep, indx_cp_l, indx_cp_r, calorically_perfect
    integer, value :: weno_scheme, weno_size, weno_version
    integer, dimension(0:ny,0:nx,0:nz) :: ep_ord_change_x_gpu
    real(rkind), dimension(1-ng:ny+ng,1-ng:nx+ng,1-ng:nz+ng,1:8) :: w_aux_trans_gpu
    real(rkind), dimension(1-ng:ny+ng,1-ng:nx+ng,1-ng:nz+ng,1:nv) :: fhat_trans_gpu
    real(rkind), dimension(ny,nz,nv,2*weno_scheme) :: gplus_x_gpu
    real(rkind), dimension(ny,nz,nv,2*weno_scheme) :: gminus_x_gpu
    real(rkind), dimension(indx_cp_l:indx_cp_r+1) :: cp_coeff_gpu
    real(rkind), dimension(4,4) :: coeff_deriv1_gpu
    real(rkind), dimension(nx) :: dcsidx_gpu
    integer, value :: force_zero_flux_min, force_zero_flux_max
    real(rkind), value :: sensor_threshold, rgas0, tol_iter_nr,rho0,u0,t0
!   Local variables
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
    integer :: ll, mm
    real(rkind) :: rho, pp, wc, gc, rhou
    real(rkind) :: tt, gamloc
    real(rkind) :: uvs5_i,uvs5_k,uvs5_p,eei,eeip
    real(rkind) :: drho, dee, eem
    real(rkind) :: drhof, deef
    real(rkind) :: sumnumrho,sumnumee,sumdenrho,sumdenee
    integer :: n,n2
!
    j = blockDim%x * (blockIdx%x - 1) + threadIdx%x
    k = blockDim%y * (blockIdx%y - 1) + threadIdx%y
    if (j > ny .or. k > nz) return
!
    do i=eul_imin-1,eul_imax
!
      ishk = 0
      do ii=i-weno_scheme+1,i+weno_scheme
        if (w_aux_trans_gpu(j,ii,k,8) > sensor_threshold) ishk = 1
      enddo
!
      if (ishk == 0) then
!
        ft1  = 0._rkind
        ft2  = 0._rkind
        ft3  = 0._rkind
        ft4  = 0._rkind
        ft5  = 0._rkind
        ft6  = 0._rkind
        lmax = max(lmax_base+ep_ord_change_x_gpu(j,i,k),1)
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
!
              rhoi  = w_aux_trans_gpu(j,i-m,k,1)
              uui   = w_aux_trans_gpu(j,i-m,k,2)
              vvi   = w_aux_trans_gpu(j,i-m,k,3)
              wwi   = w_aux_trans_gpu(j,i-m,k,4)
              enti  = w_aux_trans_gpu(j,i-m,k,5)
              tti   = w_aux_trans_gpu(j,i-m,k,6)
              ppi   = tti*rhoi*rgas0
              eei   = enti-ppi/rhoi-0.5_rkind*(uui*uui+vvi*vvi+wwi*wwi)
!
              rhoip = w_aux_trans_gpu(j,i-m+l,k,1)
              uuip  = w_aux_trans_gpu(j,i-m+l,k,2)
              vvip  = w_aux_trans_gpu(j,i-m+l,k,3)
              wwip  = w_aux_trans_gpu(j,i-m+l,k,4)
              entip = w_aux_trans_gpu(j,i-m+l,k,5)
              ttip  = w_aux_trans_gpu(j,i-m+l,k,6)
              ppip  = ttip*rhoip*rgas0
              eeip  = entip-ppip/rhoip-0.5_rkind*(uuip*uuip+vvip*vvip+wwip*wwip)
!
              rhom  = rhoi+rhoip
              eem   = eei + eeip
!
              drho   = 2._rkind*(rhoip-rhoi)/rhom
              dee    = 2._rkind*(eeip - eei)/eem
              sumnumrho = 1._rkind
              sumdenrho = 1._rkind
              sumnumee  = 1._rkind
              sumdenee  = 1._rkind
              do n = 1, nkeep
                n2 = 2*n
                sumdenrho = sumdenrho + (0.5_rkind*drho)**n2 / (1._rkind+n2)
                sumdenee  = sumdenee  + (0.5_rkind*dee )**n2
                sumnumee  = sumnumee  + (0.5_rkind*dee )**n2 / (1._rkind+n2)
              enddo
              drhof = sumnumrho/sumdenrho
              deef  = sumnumee /sumdenee
!
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
            ft1  = ft1  + coeff_deriv1_gpu(l,lmax)*uvs1
            ft2  = ft2  + coeff_deriv1_gpu(l,lmax)*uvs2
            ft3  = ft3  + coeff_deriv1_gpu(l,lmax)*uvs3
            ft4  = ft4  + coeff_deriv1_gpu(l,lmax)*uvs4
            ft5  = ft5  + coeff_deriv1_gpu(l,lmax)*(uvs5_i+uvs5_k+uvs5_p)
            ft6  = ft6  + coeff_deriv1_gpu(l,lmax)*uvs6
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
!
              rhoi  = w_aux_trans_gpu(j,i-m,k,1)
              uui   = w_aux_trans_gpu(j,i-m,k,2)
              vvi   = w_aux_trans_gpu(j,i-m,k,3)
              wwi   = w_aux_trans_gpu(j,i-m,k,4)
              enti  = w_aux_trans_gpu(j,i-m,k,5)
              tti   = w_aux_trans_gpu(j,i-m,k,6)
              ppi   = tti*rhoi*rgas0
!
              rhoip = w_aux_trans_gpu(j,i-m+l,k,1)
              uuip  = w_aux_trans_gpu(j,i-m+l,k,2)
              vvip  = w_aux_trans_gpu(j,i-m+l,k,3)
              wwip  = w_aux_trans_gpu(j,i-m+l,k,4)
              entip = w_aux_trans_gpu(j,i-m+l,k,5)
              ttip  = w_aux_trans_gpu(j,i-m+l,k,6)
              ppip  = ttip*rhoip*rgas0
!
              rhom  = rhoi+rhoip
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
!
        fh1 = 0.25_rkind*ft1
        fh2 = 0.25_rkind*ft2
        fh3 = 0.25_rkind*ft3
        fh4 = 0.25_rkind*ft4
        fh5 = 0.25_rkind*ft5
!
        if ((i==0 .and. force_zero_flux_min == 1).or.(i==nx .and. force_zero_flux_max == 1)) then
          fh1 = 0._rkind
          fh2 = 0._rkind
          fh3 = 0._rkind
          fh4 = 0._rkind
          fh5 = 0._rkind
        endif
        fh2 = fh2 + 0.5_rkind*ft6
!
        fhat_trans_gpu(j,i,k,1) = fh1
        fhat_trans_gpu(j,i,k,2) = fh2
        fhat_trans_gpu(j,i,k,3) = fh3
        fhat_trans_gpu(j,i,k,4) = fh4
        fhat_trans_gpu(j,i,k,5) = fh5
      else
        call compute_roe_average(nx, ny, nz, ng, j, j, i, i+1, k, k, w_aux_trans_gpu, rgas0, &
          b1, b2, b3, c, ci, h, uu, vv, ww, cp_coeff_gpu, indx_cp_l, indx_cp_r, &
          calorically_perfect, tol_iter_nr,t0)
!
        call eigenvectors_x(b1, b2, b3, uu, vv, ww, c, ci, h, el, er)
!
        do m=1,5 ! loop on characteristic fields
          evmax(m) = -1._rkind
        enddo
        do l=1,weno_size ! LLF
          ll   = i + l - weno_scheme
          uu   = w_aux_trans_gpu(j,ll,k,2)
          tt   = w_aux_trans_gpu(j,ll,k,6)
          gamloc = get_gamloc_dev(indx_cp_l,indx_cp_r,tt,t0,cp_coeff_gpu,calorically_perfect,rgas0)
          c    = sqrt (gamloc*rgas0*tt)
          evmax(1) = max(abs(uu-c),evmax(1))
          evmax(2) = max(abs(uu  ),evmax(2))
          evmax(3) = max(abs(uu+c),evmax(3))
          evmax(4) = evmax(2)
          evmax(5) = evmax(2)
        enddo
        do l=1,weno_size ! loop over the stencil centered at face i
          ll = i + l - weno_scheme
!
          rho    = w_aux_trans_gpu(j,ll,k,1)
          uu     = w_aux_trans_gpu(j,ll,k,2)
          vv     = w_aux_trans_gpu(j,ll,k,3)
          ww     = w_aux_trans_gpu(j,ll,k,4)
          h      = w_aux_trans_gpu(j,ll,k,5)
          rhou   = rho*uu
          pp     = rho*w_aux_trans_gpu(j,ll,k,6)*rgas0
          fi(1)  =      rhou
          fi(2)  = uu * rhou + pp
          fi(3)  = vv * rhou
          fi(4)  = ww * rhou
          fi(5)  = h  * rhou
          do m=1,5
            wc = 0._rkind
            gc = 0._rkind
!
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
!
            c = 0.5_rkind * (gc + evmax(m) * wc)
            gplus_x_gpu (j,k,m,l) = c
            gminus_x_gpu(j,k,m,l) = gc - c
          enddo
        enddo
!       
!       Reconstruction of the '+' and '-' fluxes
!       
        wenorec_ord = max(weno_scheme+ep_ord_change_x_gpu(j,i,k),1)
        call wenorec(j,k,nv,ny,nz,gplus_x_gpu,gminus_x_gpu,fi,weno_scheme,wenorec_ord,weno_version,rho0,u0)
!       
!       !Return to conservative fluxes
        do m=1,5
          fhat_trans_gpu(j,i,k,m) = 0._rkind
          do mm=1,5
            fhat_trans_gpu(j,i,k,m) = fhat_trans_gpu(j,i,k,m) + er(mm,m) * fi(mm)
          enddo
        enddo
!
      endif
    enddo
!
  endsubroutine euler_x_fluxes_hybrid_kernel
!
  attributes(global) launch_bounds(256) subroutine euler_x_fluxes_hybrid_rusanov_kernel(nv, nv_aux, nx, ny, nz, ng, &
    eul_imin, eul_imax, lmax_base, nkeep, rgas0, coeff_deriv1_gpu, dcsidx_gpu, w_aux_trans_gpu, fhat_trans_gpu, &
    force_zero_flux_min, force_zero_flux_max, &
    weno_scheme, weno_version, sensor_threshold, weno_size, gplus_x_gpu, gminus_x_gpu, cp_coeff_gpu, &
    indx_cp_l, indx_cp_r, ep_ord_change_x_gpu, calorically_perfect, tol_iter_nr,rho0,u0,t0)
!
    implicit none
!   Passed arguments
    integer, value :: nv, nx, ny, nz, ng, nv_aux
    integer, value :: eul_imin, eul_imax, lmax_base, nkeep, indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind), dimension(1-ng:ny+ng,1-ng:nx+ng,1-ng:nz+ng,1:8) :: w_aux_trans_gpu
    real(rkind), dimension(1-ng:ny+ng,1-ng:nx+ng,1-ng:nz+ng,1:nv) :: fhat_trans_gpu
    integer, dimension(0:ny,0:nx,0:nz) :: ep_ord_change_x_gpu
    real(rkind), dimension(indx_cp_l:indx_cp_r+1) :: cp_coeff_gpu
    real(rkind), dimension(4,4) :: coeff_deriv1_gpu
    real(rkind), dimension(nx) :: dcsidx_gpu
    integer, value :: force_zero_flux_min, force_zero_flux_max
    real(rkind), value :: sensor_threshold, rgas0, tol_iter_nr,rho0,u0,t0
    integer, value :: weno_scheme, weno_size, weno_version
!   Local variables
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
    real(rkind), dimension(ny,nz,nv,2*weno_scheme) :: gplus_x_gpu
    real(rkind), dimension(ny,nz,nv,2*weno_scheme) :: gminus_x_gpu
    integer :: ll, mm
    real(rkind) :: evm,evmax,rhoevm
    real(rkind) :: rho, pp, rhou
    real(rkind) :: tt, gamloc
    real(rkind) :: uvs5_i,uvs5_k,uvs5_p,eei,eeip
    real(rkind) :: drho, dee, eem
    real(rkind) :: drhof, deef
    real(rkind) :: sumnumrho,sumnumee,sumdenrho,sumdenee
    integer :: n,n2
!
    j = blockDim%x * (blockIdx%x - 1) + threadIdx%x
    k = blockDim%y * (blockIdx%y - 1) + threadIdx%y
    if (j > ny .or. k > nz) return
!
    do i=eul_imin-1,eul_imax
!
      ishk = 0
      do ii=i-weno_scheme+1,i+weno_scheme
        if (w_aux_trans_gpu(j,ii,k,8) > sensor_threshold) ishk = 1
      enddo
!
      if (ishk == 0) then
!
        ft1  = 0._rkind
        ft2  = 0._rkind
        ft3  = 0._rkind
        ft4  = 0._rkind
        ft5  = 0._rkind
        ft6  = 0._rkind
        lmax = max(lmax_base+ep_ord_change_x_gpu(j,i,k),1)
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
!
              rhoi  = w_aux_trans_gpu(j,i-m,k,1)
              uui   = w_aux_trans_gpu(j,i-m,k,2)
              vvi   = w_aux_trans_gpu(j,i-m,k,3)
              wwi   = w_aux_trans_gpu(j,i-m,k,4)
              enti  = w_aux_trans_gpu(j,i-m,k,5)
              tti   = w_aux_trans_gpu(j,i-m,k,6)
              ppi   = tti*rhoi*rgas0
              eei   = enti-ppi/rhoi-0.5_rkind*(uui*uui+vvi*vvi+wwi*wwi)
!
              rhoip = w_aux_trans_gpu(j,i-m+l,k,1)
              uuip  = w_aux_trans_gpu(j,i-m+l,k,2)
              vvip  = w_aux_trans_gpu(j,i-m+l,k,3)
              wwip  = w_aux_trans_gpu(j,i-m+l,k,4)
              entip = w_aux_trans_gpu(j,i-m+l,k,5)
              ttip  = w_aux_trans_gpu(j,i-m+l,k,6)
              ppip  = ttip*rhoip*rgas0
              eeip  = entip-ppip/rhoip-0.5_rkind*(uuip*uuip+vvip*vvip+wwip*wwip)
!
              rhom  = rhoi+rhoip
              eem   = eei + eeip
!
              drho   = 2._rkind*(rhoip-rhoi)/rhom
              dee    = 2._rkind*(eeip - eei)/eem
              sumnumrho = 1._rkind
              sumdenrho = 1._rkind
              sumnumee  = 1._rkind
              sumdenee  = 1._rkind
              do n = 1, nkeep
                n2 = 2*n
                sumdenrho = sumdenrho + (0.5_rkind*drho)**n2 / (1._rkind+n2)
                sumdenee  = sumdenee  + (0.5_rkind*dee )**n2
                sumnumee  = sumnumee  + (0.5_rkind*dee )**n2 / (1._rkind+n2)
              enddo
              drhof = sumnumrho/sumdenrho
              deef  = sumnumee /sumdenee
!
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
            ft1  = ft1  + coeff_deriv1_gpu(l,lmax)*uvs1
            ft2  = ft2  + coeff_deriv1_gpu(l,lmax)*uvs2
            ft3  = ft3  + coeff_deriv1_gpu(l,lmax)*uvs3
            ft4  = ft4  + coeff_deriv1_gpu(l,lmax)*uvs4
            ft5  = ft5  + coeff_deriv1_gpu(l,lmax)*(uvs5_i+uvs5_k+uvs5_p)
            ft6  = ft6  + coeff_deriv1_gpu(l,lmax)*uvs6
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
!
              rhoi  = w_aux_trans_gpu(j,i-m,k,1)
              uui   = w_aux_trans_gpu(j,i-m,k,2)
              vvi   = w_aux_trans_gpu(j,i-m,k,3)
              wwi   = w_aux_trans_gpu(j,i-m,k,4)
              enti  = w_aux_trans_gpu(j,i-m,k,5)
              tti   = w_aux_trans_gpu(j,i-m,k,6)
              ppi   = tti*rhoi*rgas0
!
              rhoip = w_aux_trans_gpu(j,i-m+l,k,1)
              uuip  = w_aux_trans_gpu(j,i-m+l,k,2)
              vvip  = w_aux_trans_gpu(j,i-m+l,k,3)
              wwip  = w_aux_trans_gpu(j,i-m+l,k,4)
              entip = w_aux_trans_gpu(j,i-m+l,k,5)
              ttip  = w_aux_trans_gpu(j,i-m+l,k,6)
              ppip  = ttip*rhoip*rgas0
!
              rhom  = rhoi+rhoip
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
!
        fh1 = 0.25_rkind*ft1
        fh2 = 0.25_rkind*ft2
        fh3 = 0.25_rkind*ft3
        fh4 = 0.25_rkind*ft4
        fh5 = 0.25_rkind*ft5
!
        if ((i==0 .and. force_zero_flux_min == 1).or.(i==nx .and. force_zero_flux_max == 1)) then
          fh1 = 0._rkind
          fh2 = 0._rkind
          fh3 = 0._rkind
          fh4 = 0._rkind
          fh5 = 0._rkind
        endif
        fh2 = fh2 + 0.5_rkind*ft6
!
        fhat_trans_gpu(j,i,k,1) = fh1
        fhat_trans_gpu(j,i,k,2) = fh2
        fhat_trans_gpu(j,i,k,3) = fh3
        fhat_trans_gpu(j,i,k,4) = fh4
        fhat_trans_gpu(j,i,k,5) = fh5
      else
!
        evmax = -1._rkind
        do l=1,weno_size ! LLF
          ll   = i + l - weno_scheme
          uu   = w_aux_trans_gpu(j,ll,k,2)
          tt   = w_aux_trans_gpu(j,ll,k,6)
          gamloc = get_gamloc_dev(indx_cp_l,indx_cp_r,tt,t0,cp_coeff_gpu,calorically_perfect,rgas0)
          c      = sqrt (gamloc*rgas0*tt)
          evm    = max(abs(uu-c),abs(uu+c))
          evmax  = max(evm,evmax)
        enddo
        do l=1,weno_size ! loop over the stencil centered at face i
          ll = i + l - weno_scheme
          rho    = w_aux_trans_gpu(j,ll,k,1)
          uu     = w_aux_trans_gpu(j,ll,k,2)
          vv     = w_aux_trans_gpu(j,ll,k,3)
          ww     = w_aux_trans_gpu(j,ll,k,4)
          h      = w_aux_trans_gpu(j,ll,k,5)
          rhou   = rho*uu
          pp     = rho*w_aux_trans_gpu(j,ll,k,6)*rgas0
          rhoevm = rho*evmax
!
          evm    = rhou
          c = 0.5_rkind * (evm + rhoevm)
          gplus_x_gpu (j,k,1,l) = c
          gminus_x_gpu(j,k,1,l) = evm - c
          evm    = uu * rhou + pp
          c = 0.5_rkind * (evm + rhoevm * uu)
          gplus_x_gpu (j,k,2,l) = c
          gminus_x_gpu(j,k,2,l) = evm - c
          evm    = vv * rhou
          c = 0.5_rkind * (evm + rhoevm * vv)
          gplus_x_gpu (j,k,3,l) = c
          gminus_x_gpu(j,k,3,l) = evm - c
          evm    = ww * rhou
          c = 0.5_rkind * (evm + rhoevm * ww)
          gplus_x_gpu (j,k,4,l) = c
          gminus_x_gpu(j,k,4,l) = evm - c
          evm    = h  * rhou
          c = 0.5_rkind * (evm + evmax * (rho*h-pp))
          gplus_x_gpu (j,k,5,l) = c
          gminus_x_gpu(j,k,5,l) = evm - c
        enddo
!       
!       Reconstruction of the '+' and '-' fluxes
!       
        wenorec_ord = max(weno_scheme+ep_ord_change_x_gpu(j,i,k),1)
        call wenorec_rusanov(j,k,nv,ny,nz,gplus_x_gpu,gminus_x_gpu,fhat_trans_gpu,&
          weno_scheme,wenorec_ord,weno_version,ng,ny,nx,nz,j,i,k)
!       
      endif
    enddo
!
  endsubroutine euler_x_fluxes_hybrid_rusanov_kernel
!
  subroutine euler_x_update_cuf(nx, ny, nz, ng, nv, eul_imin, eul_imax, &
    fhat_trans_gpu,fl_trans_gpu,fl_gpu,dcsidx_gpu,stream_id)
    integer :: nx,ny,nz,ng,nv,eul_imin,eul_imax
    real(rkind), dimension(1-ng:ny+ng,1-ng:nx+ng,1-ng:nz+ng,1:nv), intent(in), device :: fhat_trans_gpu
    real(rkind), dimension(1:ny,1:nx,1:nz,1:nv), intent(inout), device :: fl_trans_gpu
    real(rkind), dimension(1:nx,1:ny,1:nz,1:nv), intent(inout), device :: fl_gpu
    real(rkind), dimension(1:nx), intent(in), device :: dcsidx_gpu
    integer(kind=cuda_stream_kind) :: stream_id
    integer :: i,j,k,m,iv,iercuda
!
    !$cuf kernel do(3) <<<*,*,stream=stream_id>>>
    do k=1,nz
      do i=eul_imin,eul_imax ! loop on the inner nodes
        do j=1,ny
          do m=1,nv
            fl_trans_gpu(j,i,k,m) = (fhat_trans_gpu(j,i,k,m)-fhat_trans_gpu(j,i-1,k,m))*dcsidx_gpu(i)
          enddo
        enddo
      enddo
    enddo
!
!   Note: transposed is done now alongside flux update. Hope performance are good!
    !$cuf kernel do(3) <<<*,*,stream=stream_id>>>
    do k=1,nz
      do j=1,ny
        do i=eul_imin,eul_imax
          do iv=1,nv
            fl_gpu(i,j,k,iv) = fl_gpu(i,j,k,iv) + fl_trans_gpu(j,i,k,iv)
          enddo
        enddo
      enddo
    enddo
  endsubroutine euler_x_update_cuf
!
  attributes(global) launch_bounds(256) subroutine euler_z_hybrid_kernel(nv, nv_aux, nx, ny, nz, ng, &
    eul_kmin, eul_kmax, lmax_base, nkeep, rgas0, w_aux_gpu, fl_gpu, coeff_deriv1_gpu, dzitdz_gpu, fhat_gpu, &
    force_zero_flux_min, force_zero_flux_max, &
    weno_scheme, weno_version, sensor_threshold, weno_size, w_gpu, gplus_z_gpu, gminus_z_gpu, cp_coeff_gpu, indx_cp_l, indx_cp&
    &_r, &
    ep_ord_change_gpu, calorically_perfect, tol_iter_nr,rho0,u0,t0)
    implicit none
!   Passed arguments
    integer, value :: nv, nx, ny, nz, ng, nv_aux
    integer, value :: eul_kmin, eul_kmax, lmax_base, nkeep, indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv_aux) :: w_aux_gpu
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv) :: w_gpu
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv) :: fhat_gpu
    integer, dimension(0:nx,0:ny,0:nz,2:3) :: ep_ord_change_gpu
    real(rkind), dimension(indx_cp_l:indx_cp_r+1) :: cp_coeff_gpu
    real(rkind), dimension(nx,ny,nz,nv) :: fl_gpu
    real(rkind), dimension(4,4) :: coeff_deriv1_gpu
    real(rkind), dimension(nz) :: dzitdz_gpu
    integer, value :: force_zero_flux_min, force_zero_flux_max
    real(rkind), value :: sensor_threshold, rgas0, tol_iter_nr, rho0, u0, t0
    integer, value :: weno_scheme, weno_size, weno_version
!   Local variables
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
    real(rkind), dimension(nx,ny,nv,2*weno_scheme) :: gplus_z_gpu
    real(rkind), dimension(nx,ny,nv,2*weno_scheme) :: gminus_z_gpu
    integer :: ll, mm
    real(rkind) :: rho, pp, wc, gc, rhow
    real(rkind) :: tt, gamloc
    real(rkind) :: uvs5_i,uvs5_k,uvs5_p,eei,eeip
    real(rkind) :: drho, dee, eem
    real(rkind) :: drhof, deef
    real(rkind) :: sumnumrho,sumnumee,sumdenrho,sumdenee
    integer :: n,n2
!
    i = blockDim%x * (blockIdx%x - 1) + threadIdx%x
    j = blockDim%y * (blockIdx%y - 1) + threadIdx%y
    if(i > nx .or. j > ny) return
!
    do k=eul_kmin-1,eul_kmax
      ishk = 0
      do kk=k-weno_scheme+1,k+weno_scheme
        if (w_aux_gpu(i,j,kk,8) > sensor_threshold) ishk = 1
      enddo
!
      if (ishk == 0) then
        ft1  = 0._rkind
        ft2  = 0._rkind
        ft3  = 0._rkind
        ft4  = 0._rkind
        ft5  = 0._rkind
        ft6  = 0._rkind
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
!
              rhoi  = w_aux_gpu(i,j,k-m,1)
              uui   = w_aux_gpu(i,j,k-m,2)
              vvi   = w_aux_gpu(i,j,k-m,3)
              wwi   = w_aux_gpu(i,j,k-m,4)
              enti  = w_aux_gpu(i,j,k-m,5)
              tti   = w_aux_gpu(i,j,k-m,6)
              ppi   = tti*rhoi*rgas0
              eei   = enti-ppi/rhoi-0.5_rkind*(uui*uui+vvi*vvi+wwi*wwi)
!
              rhoip = w_aux_gpu(i,j,k-m+l,1)
              uuip  = w_aux_gpu(i,j,k-m+l,2)
              vvip  = w_aux_gpu(i,j,k-m+l,3)
              wwip  = w_aux_gpu(i,j,k-m+l,4)
              entip = w_aux_gpu(i,j,k-m+l,5)
              ttip  = w_aux_gpu(i,j,k-m+l,6)
              ppip  = ttip*rhoip*rgas0
              eeip  = entip-ppip/rhoip-0.5_rkind*(uuip*uuip+vvip*vvip+wwip*wwip)
!
              rhom  = rhoi + rhoip
              eem   = eei  + eeip
!
              drho   = 2._rkind*(rhoip-rhoi)/rhom
              dee    = 2._rkind*(eeip - eei)/eem
              sumnumrho = 1._rkind
              sumdenrho = 1._rkind
              sumnumee  = 1._rkind
              sumdenee  = 1._rkind
              do n = 1, nkeep
                n2 = 2*n
                sumdenrho = sumdenrho + (0.5_rkind*drho)**n2 / (1._rkind+n2)
                sumdenee  = sumdenee  + (0.5_rkind*dee )**n2
                sumnumee  = sumnumee  + (0.5_rkind*dee )**n2 / (1._rkind+n2)
              enddo
              drhof = sumnumrho/sumdenrho
              deef  = sumnumee /sumdenee
!
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
            ft1  = ft1  + coeff_deriv1_gpu(l,lmax)*uvs1
            ft2  = ft2  + coeff_deriv1_gpu(l,lmax)*uvs2
            ft3  = ft3  + coeff_deriv1_gpu(l,lmax)*uvs3
            ft4  = ft4  + coeff_deriv1_gpu(l,lmax)*uvs4
            ft5  = ft5  + coeff_deriv1_gpu(l,lmax)*(uvs5_i+uvs5_k+uvs5_p)
            ft6  = ft6  + coeff_deriv1_gpu(l,lmax)*uvs6
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
!
              rhoi  = w_aux_gpu(i,j,k-m,1)
              uui   = w_aux_gpu(i,j,k-m,2)
              vvi   = w_aux_gpu(i,j,k-m,3)
              wwi   = w_aux_gpu(i,j,k-m,4)
              enti  = w_aux_gpu(i,j,k-m,5)
              tti   = w_aux_gpu(i,j,k-m,6)
              ppi   = tti*rhoi*rgas0
!
              rhoip = w_aux_gpu(i,j,k-m+l,1)
              uuip  = w_aux_gpu(i,j,k-m+l,2)
              vvip  = w_aux_gpu(i,j,k-m+l,3)
              wwip  = w_aux_gpu(i,j,k-m+l,4)
              entip = w_aux_gpu(i,j,k-m+l,5)
              ttip  = w_aux_gpu(i,j,k-m+l,6)
              ppip  = ttip*rhoip*rgas0
!
              rhom  = rhoi+rhoip
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
!
        fhat_gpu(i,j,k,1) = fh1
        fhat_gpu(i,j,k,2) = fh2
        fhat_gpu(i,j,k,3) = fh3
        fhat_gpu(i,j,k,4) = fh4
        fhat_gpu(i,j,k,5) = fh5
      else
        call compute_roe_average(nx, ny, nz, ng, i, i, j, j, k, k+1, w_aux_gpu, rgas0, &
          b1, b2, b3, c, ci, h, uu, vv, ww, cp_coeff_gpu, indx_cp_l, indx_cp_r, &
          calorically_perfect, tol_iter_nr,t0)
!
        call eigenvectors_z(b1, b2, b3, uu, vv, ww, c, ci, h, el, er)
!
        do m=1,5 ! loop on characteristic fields
          evmax(m) = -1._rkind
        enddo
        do l=1,weno_size ! LLF
          ll = k + l - weno_scheme
          ww   = w_aux_gpu(i,j,ll,4)
          tt   = w_aux_gpu(i,j,ll,6)
          gamloc = get_gamloc_dev(indx_cp_l,indx_cp_r,tt,t0,cp_coeff_gpu,calorically_perfect,rgas0)
          c    = sqrt (gamloc*rgas0*tt)
          evmax(1) = max(abs(ww-c),evmax(1))
          evmax(2) = max(abs(ww  ),evmax(2))
          evmax(3) = max(abs(ww+c),evmax(3))
          evmax(4) = evmax(2)
          evmax(5) = evmax(2)
        enddo
        do l=1,weno_size ! loop over the stencil centered at face i
          ll = k + l - weno_scheme
!
          rho    = w_aux_gpu(i,j,ll,1)
          uu     = w_aux_gpu(i,j,ll,2)
          vv     = w_aux_gpu(i,j,ll,3)
          ww     = w_aux_gpu(i,j,ll,4)
          h      = w_aux_gpu(i,j,ll,5)
          rhow   = rho*ww
          pp     = rho*w_aux_gpu(i,j,ll,6)*rgas0
          fk(1)  =      rhow
          fk(2)  = uu * rhow
          fk(3)  = vv * rhow
          fk(4)  = ww * rhow + pp
          fk(5)  = h  * rhow
          do m=1,5
            wc = 0._rkind
            gc = 0._rkind
!
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
!           do mm=1,5
!             wc = wc + el(mm,m) * w_gpu(i,j,ll,mm)
!             gc = gc + el(mm,m) * fk(mm)
!           enddo
!
            c = 0.5_rkind * (gc + evmax(m) * wc)
            gplus_z_gpu (i,j,m,l) = c
            gminus_z_gpu(i,j,m,l) = gc - c
          enddo
        enddo
!       
!       Reconstruction of the '+' and '-' fluxes
!       
        wenorec_ord = max(weno_scheme+ep_ord_change_gpu(i,j,k,3),1)
        call wenorec(i,j,nv,nx,ny,gplus_z_gpu,gminus_z_gpu,fk,weno_scheme,wenorec_ord,weno_version,rho0,u0)
!       
!       !Return to conservative fluxes
        do m=1,5
          fhat_gpu(i,j,k,m) = 0._rkind
          do mm=1,5
            fhat_gpu(i,j,k,m) = fhat_gpu(i,j,k,m) + er(mm,m) * fk(mm)
          enddo
        enddo
!
      endif
    enddo
!
!   Update net flux
    do k=eul_kmin,eul_kmax ! loop on the inner nodes
      do m=1,5
        fl_gpu(i,j,k,m) = fl_gpu(i,j,k,m) + (fhat_gpu(i,j,k,m)-fhat_gpu(i,j,k-1,m))*dzitdz_gpu(k)
      enddo
    enddo
!
  endsubroutine euler_z_hybrid_kernel
!
  attributes(global) launch_bounds(256) subroutine euler_z_hybrid_rusanov_kernel(nv, nv_aux, nx, ny, nz, ng, &
    eul_kmin, eul_kmax, lmax_base, nkeep, rgas0, w_aux_gpu, fl_gpu, coeff_deriv1_gpu, dzitdz_gpu, fhat_gpu, &
    force_zero_flux_min, force_zero_flux_max, &
    weno_scheme, weno_version, sensor_threshold, weno_size, w_gpu, gplus_z_gpu, gminus_z_gpu, cp_coeff_gpu, indx_cp_l, indx_cp&
    &_r, &
    ep_ord_change_gpu, calorically_perfect, tol_iter_nr,rho0,u0,t0)
    implicit none
!   Passed arguments
    integer, value :: nv, nx, ny, nz, ng, nv_aux
    integer, value :: eul_kmin, eul_kmax, lmax_base, nkeep, indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv_aux) :: w_aux_gpu
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv) :: w_gpu
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv) :: fhat_gpu
    integer, dimension(0:nx,0:ny,0:nz,2:3) :: ep_ord_change_gpu
    real(rkind), dimension(indx_cp_l:indx_cp_r+1) :: cp_coeff_gpu
    real(rkind), dimension(nx,ny,nz,nv) :: fl_gpu
    real(rkind), dimension(4,4) :: coeff_deriv1_gpu
    real(rkind), dimension(nz) :: dzitdz_gpu
    integer, value :: force_zero_flux_min, force_zero_flux_max
    real(rkind), value :: sensor_threshold, rgas0, tol_iter_nr, rho0, u0, t0
    integer, value :: weno_scheme, weno_size, weno_version
!   Local variables
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
    real(rkind), dimension(nx,ny,nv,2*weno_scheme) :: gplus_z_gpu
    real(rkind), dimension(nx,ny,nv,2*weno_scheme) :: gminus_z_gpu
    integer :: ll, mm
    real(rkind) :: evm, evmax, rhoevm
    real(rkind) :: rho, pp, rhow
    real(rkind) :: tt, gamloc
    real(rkind) :: uvs5_i,uvs5_k,uvs5_p,eei,eeip
    real(rkind) :: drho, dee, eem
    real(rkind) :: drhof, deef
    real(rkind) :: sumnumrho,sumnumee,sumdenrho,sumdenee
    integer :: n,n2
!
    i = blockDim%x * (blockIdx%x - 1) + threadIdx%x
    j = blockDim%y * (blockIdx%y - 1) + threadIdx%y
    if(i > nx .or. j > ny) return
!
    do k=eul_kmin-1,eul_kmax
      ishk = 0
      do kk=k-weno_scheme+1,k+weno_scheme
        if (w_aux_gpu(i,j,kk,8) > sensor_threshold) ishk = 1
      enddo
!
      if (ishk == 0) then
        ft1  = 0._rkind
        ft2  = 0._rkind
        ft3  = 0._rkind
        ft4  = 0._rkind
        ft5  = 0._rkind
        ft6  = 0._rkind
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
!
              rhoi  = w_aux_gpu(i,j,k-m,1)
              uui   = w_aux_gpu(i,j,k-m,2)
              vvi   = w_aux_gpu(i,j,k-m,3)
              wwi   = w_aux_gpu(i,j,k-m,4)
              enti  = w_aux_gpu(i,j,k-m,5)
              tti   = w_aux_gpu(i,j,k-m,6)
              ppi   = tti*rhoi*rgas0
              eei   = enti-ppi/rhoi-0.5_rkind*(uui*uui+vvi*vvi+wwi*wwi)
!
              rhoip = w_aux_gpu(i,j,k-m+l,1)
              uuip  = w_aux_gpu(i,j,k-m+l,2)
              vvip  = w_aux_gpu(i,j,k-m+l,3)
              wwip  = w_aux_gpu(i,j,k-m+l,4)
              entip = w_aux_gpu(i,j,k-m+l,5)
              ttip  = w_aux_gpu(i,j,k-m+l,6)
              ppip  = ttip*rhoip*rgas0
              eeip  = entip-ppip/rhoip-0.5_rkind*(uuip*uuip+vvip*vvip+wwip*wwip)
!
              rhom  = rhoi + rhoip
              eem   = eei  + eeip
!
              drho   = 2._rkind*(rhoip-rhoi)/rhom
              dee    = 2._rkind*(eeip - eei)/eem
              sumnumrho = 1._rkind
              sumdenrho = 1._rkind
              sumnumee  = 1._rkind
              sumdenee  = 1._rkind
              do n = 1, nkeep
                n2 = 2*n
                sumdenrho = sumdenrho + (0.5_rkind*drho)**n2 / (1._rkind+n2)
                sumdenee  = sumdenee  + (0.5_rkind*dee )**n2
                sumnumee  = sumnumee  + (0.5_rkind*dee )**n2 / (1._rkind+n2)
              enddo
              drhof = sumnumrho/sumdenrho
              deef  = sumnumee /sumdenee
!
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
            ft1  = ft1  + coeff_deriv1_gpu(l,lmax)*uvs1
            ft2  = ft2  + coeff_deriv1_gpu(l,lmax)*uvs2
            ft3  = ft3  + coeff_deriv1_gpu(l,lmax)*uvs3
            ft4  = ft4  + coeff_deriv1_gpu(l,lmax)*uvs4
            ft5  = ft5  + coeff_deriv1_gpu(l,lmax)*(uvs5_i+uvs5_k+uvs5_p)
            ft6  = ft6  + coeff_deriv1_gpu(l,lmax)*uvs6
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
!
              rhoi  = w_aux_gpu(i,j,k-m,1)
              uui   = w_aux_gpu(i,j,k-m,2)
              vvi   = w_aux_gpu(i,j,k-m,3)
              wwi   = w_aux_gpu(i,j,k-m,4)
              enti  = w_aux_gpu(i,j,k-m,5)
              tti   = w_aux_gpu(i,j,k-m,6)
              ppi   = tti*rhoi*rgas0
!
              rhoip = w_aux_gpu(i,j,k-m+l,1)
              uuip  = w_aux_gpu(i,j,k-m+l,2)
              vvip  = w_aux_gpu(i,j,k-m+l,3)
              wwip  = w_aux_gpu(i,j,k-m+l,4)
              entip = w_aux_gpu(i,j,k-m+l,5)
              ttip  = w_aux_gpu(i,j,k-m+l,6)
              ppip  = ttip*rhoip*rgas0
!
              rhom  = rhoi+rhoip
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
!
        fhat_gpu(i,j,k,1) = fh1
        fhat_gpu(i,j,k,2) = fh2
        fhat_gpu(i,j,k,3) = fh3
        fhat_gpu(i,j,k,4) = fh4
        fhat_gpu(i,j,k,5) = fh5
      else
        evmax = -1._rkind
        do l=1,weno_size ! LLF
          ll = k + l - weno_scheme
          ww   = w_aux_gpu(i,j,ll,4)
          tt   = w_aux_gpu(i,j,ll,6)
          gamloc = get_gamloc_dev(indx_cp_l,indx_cp_r,tt,t0,cp_coeff_gpu,calorically_perfect,rgas0)
          c     = sqrt (gamloc*rgas0*tt)
          evm   = max(abs(ww-c),abs(ww+c))
          evmax = max(evm,evmax)
        enddo
        do l=1,weno_size ! loop over the stencil centered at face i
          ll = k + l - weno_scheme
!
          rho    = w_aux_gpu(i,j,ll,1)
          uu     = w_aux_gpu(i,j,ll,2)
          vv     = w_aux_gpu(i,j,ll,3)
          ww     = w_aux_gpu(i,j,ll,4)
          h      = w_aux_gpu(i,j,ll,5)
          rhow   = rho*ww
          pp     = rho*w_aux_gpu(i,j,ll,6)*rgas0
          rhoevm = rho*evmax
!
          evm = rhow
          c = 0.5_rkind * (evm + rhoevm)
          gplus_z_gpu (i,j,1,l) = c
          gminus_z_gpu(i,j,1,l) = evm-c
          evm  = uu * rhow
          c = 0.5_rkind * (evm + rhoevm * uu)
          gplus_z_gpu (i,j,2,l) = c
          gminus_z_gpu(i,j,2,l) = evm-c
          evm = vv * rhow
          c = 0.5_rkind * (evm + rhoevm * vv)
          gplus_z_gpu (i,j,3,l) = c
          gminus_z_gpu(i,j,3,l) = evm-c
          evm = ww * rhow + pp
          c = 0.5_rkind * (evm + rhoevm * ww)
          gplus_z_gpu (i,j,4,l) = c
          gminus_z_gpu(i,j,4,l) = evm-c
          evm = h  * rhow
          c = 0.5_rkind * (evm + evmax * (rho*h-pp))
          gplus_z_gpu (i,j,5,l) = c
          gminus_z_gpu(i,j,5,l) = evm-c
        enddo
!       
!       Reconstruction of the '+' and '-' fluxes
!       
        wenorec_ord = max(weno_scheme+ep_ord_change_gpu(i,j,k,3),1)
        call wenorec_rusanov(i,j,nv,nx,ny,gplus_z_gpu,gminus_z_gpu,fhat_gpu,&
          weno_scheme,wenorec_ord,weno_version,ng,nx,ny,nz,i,j,k)
!
      endif
    enddo
!
!   Update net flux
    do k=eul_kmin,eul_kmax ! loop on the inner nodes
      do m=1,5
        fl_gpu(i,j,k,m) = fl_gpu(i,j,k,m) + (fhat_gpu(i,j,k,m)-fhat_gpu(i,j,k-1,m))*dzitdz_gpu(k)
      enddo
    enddo
!
  endsubroutine euler_z_hybrid_rusanov_kernel
!
  attributes(global) launch_bounds(256) subroutine euler_y_hybrid_kernel(nv, nv_aux, nx, ny, nz, ng, &
    eul_jmin, eul_jmax, lmax_base, nkeep, rgas0, w_aux_gpu, fl_gpu, coeff_deriv1_gpu, detady_gpu, fhat_gpu, &
    force_zero_flux_min, force_zero_flux_max, &
    weno_scheme, weno_version, sensor_threshold, weno_size, w_gpu, gplus_y_gpu, gminus_y_gpu, cp_coeff_gpu, indx_cp_l, indx_cp&
    &_r, &
    ep_ord_change_gpu, calorically_perfect, tol_iter_nr,rho0,u0,t0)
    implicit none
!   Passed arguments
    integer, value :: nv, nx, ny, nz, ng, nv_aux
    integer, value :: eul_jmin, eul_jmax, lmax_base, nkeep, indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv_aux) :: w_aux_gpu
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv) :: w_gpu
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv) :: fhat_gpu
    integer, dimension(0:nx,0:ny,0:nz,2:3) :: ep_ord_change_gpu
    real(rkind), dimension(indx_cp_l:indx_cp_r+1) :: cp_coeff_gpu
    real(rkind), dimension(nx,ny,nz,nv) :: fl_gpu
    real(rkind), dimension(4,4) :: coeff_deriv1_gpu
    real(rkind), dimension(ny) :: detady_gpu
    integer, value :: force_zero_flux_min, force_zero_flux_max
    real(rkind), value :: sensor_threshold, rgas0, tol_iter_nr,rho0,u0,t0
    integer, value :: weno_scheme, weno_size, weno_version
!   Local variables
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
    real(rkind), dimension(nx,nz,nv,2*weno_scheme) :: gplus_y_gpu
    real(rkind), dimension(nx,nz,nv,2*weno_scheme) :: gminus_y_gpu
    integer :: ll, mm, lmax, wenorec_ord
    real(rkind) :: rho, pp, wc, gc, rhov
    real(rkind) :: tt, gamloc
    real(rkind) :: uvs5_i,uvs5_k,uvs5_p,eei,eeip
    real(rkind) :: drho, dee, eem
    real(rkind) :: drhof, deef
    real(rkind) :: sumnumrho,sumnumee,sumdenrho,sumdenee
    integer :: n,n2
!
    i = blockDim%x * (blockIdx%x - 1) + threadIdx%x
    k = blockDim%y * (blockIdx%y - 1) + threadIdx%y
    if (i > nx .or. k > nz) return
!
    do j=eul_jmin-1,eul_jmax
!
      ishk = 0
      do jj=j-weno_scheme+1,j+weno_scheme
        if (w_aux_gpu(i,jj,k,8) > sensor_threshold) ishk = 1
      enddo
!
      if (ishk == 0) then
        ft1  = 0._rkind
        ft2  = 0._rkind
        ft3  = 0._rkind
        ft4  = 0._rkind
        ft5  = 0._rkind
        ft6  = 0._rkind
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
!
              rhoi  = w_aux_gpu(i,j-m,k,1)
              uui   = w_aux_gpu(i,j-m,k,2)
              vvi   = w_aux_gpu(i,j-m,k,3)
              wwi   = w_aux_gpu(i,j-m,k,4)
              enti  = w_aux_gpu(i,j-m,k,5)
              tti   = w_aux_gpu(i,j-m,k,6)
              ppi   = tti*rhoi*rgas0
              eei   = enti-ppi/rhoi-0.5_rkind*(uui*uui+vvi*vvi+wwi*wwi)
!
              rhoip = w_aux_gpu(i,j-m+l,k,1)
              uuip  = w_aux_gpu(i,j-m+l,k,2)
              vvip  = w_aux_gpu(i,j-m+l,k,3)
              wwip  = w_aux_gpu(i,j-m+l,k,4)
              entip = w_aux_gpu(i,j-m+l,k,5)
              ttip  = w_aux_gpu(i,j-m+l,k,6)
              ppip  = ttip*rhoip*rgas0
              eeip  = entip-ppip/rhoip-0.5_rkind*(uuip*uuip+vvip*vvip+wwip*wwip)
!
              rhom  = rhoi + rhoip
              eem   = eei  + eeip
!
              drho   = 2._rkind*(rhoip-rhoi)/rhom
              dee    = 2._rkind*(eeip - eei)/eem
              sumnumrho = 1._rkind
              sumdenrho = 1._rkind
              sumnumee  = 1._rkind
              sumdenee  = 1._rkind
              do n = 1, nkeep
                n2 = 2*n
                sumdenrho = sumdenrho + (0.5_rkind*drho)**n2 / (1._rkind+n2)
                sumdenee  = sumdenee  + (0.5_rkind*dee )**n2
                sumnumee  = sumnumee  + (0.5_rkind*dee )**n2 / (1._rkind+n2)
              enddo
              drhof = sumnumrho/sumdenrho
              deef  = sumnumee /sumdenee
!
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
            ft1  = ft1  + coeff_deriv1_gpu(l,lmax)*uvs1
            ft2  = ft2  + coeff_deriv1_gpu(l,lmax)*uvs2
            ft3  = ft3  + coeff_deriv1_gpu(l,lmax)*uvs3
            ft4  = ft4  + coeff_deriv1_gpu(l,lmax)*uvs4
            ft5  = ft5  + coeff_deriv1_gpu(l,lmax)*(uvs5_i+uvs5_k+uvs5_p)
            ft6  = ft6  + coeff_deriv1_gpu(l,lmax)*uvs6
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
!
              rhoi  = w_aux_gpu(i,j-m,k,1)
              uui   = w_aux_gpu(i,j-m,k,2)
              vvi   = w_aux_gpu(i,j-m,k,3)
              wwi   = w_aux_gpu(i,j-m,k,4)
              enti  = w_aux_gpu(i,j-m,k,5)
              tti   = w_aux_gpu(i,j-m,k,6)
              ppi   = tti*rhoi*rgas0
!
              rhoip = w_aux_gpu(i,j-m+l,k,1)
              uuip  = w_aux_gpu(i,j-m+l,k,2)
              vvip  = w_aux_gpu(i,j-m+l,k,3)
              wwip  = w_aux_gpu(i,j-m+l,k,4)
              entip = w_aux_gpu(i,j-m+l,k,5)
              ttip  = w_aux_gpu(i,j-m+l,k,6)
              ppip  = ttip*rhoip*rgas0
!
              rhom  = rhoi + rhoip
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
!
        fh3 = fh3 + 0.5_rkind*ft6
!
        fhat_gpu(i,j,k,1) = fh1
        fhat_gpu(i,j,k,2) = fh2
        fhat_gpu(i,j,k,3) = fh3
        fhat_gpu(i,j,k,4) = fh4
        fhat_gpu(i,j,k,5) = fh5
      else
        call compute_roe_average(nx, ny, nz, ng, i, i, j, j+1, k, k, w_aux_gpu, rgas0, &
          b1, b2, b3, c, ci, h, uu, vv, ww, cp_coeff_gpu, indx_cp_l, indx_cp_r, &
          calorically_perfect, tol_iter_nr,t0)
!
        call eigenvectors_y(b1, b2, b3, uu, vv, ww, c, ci, h, el, er)
!
        do m=1,5 ! loop on characteristic fields
          evmax(m) = -1._rkind
        enddo
        do l=1,weno_size ! LLF
          ll = j + l - weno_scheme
          vv   = w_aux_gpu(i,ll,k,3)
          tt   = w_aux_gpu(i,ll,k,6)
          gamloc = get_gamloc_dev(indx_cp_l,indx_cp_r,tt,t0,cp_coeff_gpu,calorically_perfect,rgas0)
          c    = sqrt (gamloc*rgas0*tt)
          evmax(1) = max(abs(vv-c),evmax(1))
          evmax(2) = max(abs(vv  ),evmax(2))
          evmax(3) = max(abs(vv+c),evmax(3))
          evmax(4) = evmax(2)
          evmax(5) = evmax(2)
        enddo
        do l=1,weno_size ! loop over the stencil centered at face i
          ll = j + l - weno_scheme
!
          rho    = w_aux_gpu(i,ll,k,1)
          uu     = w_aux_gpu(i,ll,k,2)
          vv     = w_aux_gpu(i,ll,k,3)
          ww     = w_aux_gpu(i,ll,k,4)
          h      = w_aux_gpu(i,ll,k,5)
          rhov   = rho*vv
          pp     = rho*w_aux_gpu(i,ll,k,6)*rgas0
          fj(1)  =      rhov
          fj(2)  = uu * rhov
          fj(3)  = vv * rhov + pp
          fj(4)  = ww * rhov
          fj(5)  = h  * rhov
          do m=1,5
            wc = 0._rkind
            gc = 0._rkind
!
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
!
            c = 0.5_rkind * (gc + evmax(m) * wc)
            gplus_y_gpu (i,k,m,l) = c
            gminus_y_gpu(i,k,m,l) = gc - c
          enddo
        enddo
!       
!       Reconstruction of the '+' and '-' fluxes
!       
        wenorec_ord = max(weno_scheme+ep_ord_change_gpu(i,j,k,2),1)
        call wenorec(i,k,nv,nx,nz,gplus_y_gpu,gminus_y_gpu,fj,weno_scheme,wenorec_ord,weno_version,rho0,u0)
!       
!       !Return to conservative fluxes
        do m=1,5
          fhat_gpu(i,j,k,m) = 0._rkind
          do mm=1,5
            fhat_gpu(i,j,k,m) = fhat_gpu(i,j,k,m) + er(mm,m) * fj(mm)
          enddo
        enddo
!
      endif
    enddo
!
!   Update net flux
    do j=eul_jmin,eul_jmax ! loop on the inner nodes
      do m=1,5
        fl_gpu(i,j,k,m) = fl_gpu(i,j,k,m) + (fhat_gpu(i,j,k,m)-fhat_gpu(i,j-1,k,m))*detady_gpu(j)
      enddo
    enddo
!
  endsubroutine euler_y_hybrid_kernel
!
  attributes(global) launch_bounds(256) subroutine euler_y_hybrid_rusanov_kernel(nv, nv_aux, nx, ny, nz, ng, &
    eul_jmin, eul_jmax, lmax_base, nkeep, rgas0, w_aux_gpu, fl_gpu, coeff_deriv1_gpu, detady_gpu, fhat_gpu, &
    force_zero_flux_min, force_zero_flux_max, &
    weno_scheme, weno_version, sensor_threshold, weno_size, w_gpu, gplus_y_gpu, gminus_y_gpu, cp_coeff_gpu, indx_cp_l, indx_cp&
    &_r, &
    ep_ord_change_gpu, calorically_perfect, tol_iter_nr,rho0,u0,t0)
    implicit none
!   Passed arguments
    integer, value :: nv, nx, ny, nz, ng, nv_aux
    integer, value :: eul_jmin, eul_jmax, lmax_base, nkeep, indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv_aux) :: w_aux_gpu
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv) :: w_gpu
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv) :: fhat_gpu
    integer, dimension(0:nx,0:ny,0:nz,2:3) :: ep_ord_change_gpu
    real(rkind), dimension(indx_cp_l:indx_cp_r+1) :: cp_coeff_gpu
    real(rkind), dimension(nx,ny,nz,nv) :: fl_gpu
    real(rkind), dimension(4,4) :: coeff_deriv1_gpu
    real(rkind), dimension(ny) :: detady_gpu
    integer, value :: force_zero_flux_min, force_zero_flux_max
    real(rkind), value :: sensor_threshold, rgas0, tol_iter_nr,rho0,u0,t0
    integer, value :: weno_scheme, weno_size, weno_version
!   Local variables
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
    real(rkind), dimension(nx,nz,nv,2*weno_scheme) :: gplus_y_gpu
    real(rkind), dimension(nx,nz,nv,2*weno_scheme) :: gminus_y_gpu
    integer :: ll, mm, lmax, wenorec_ord
    real(rkind) :: evm, evmax, rhoevm
    real(rkind) :: rho, pp, rhov
    real(rkind) :: tt, gamloc
    real(rkind) :: uvs5_i,uvs5_k,uvs5_p,eei,eeip
    real(rkind) :: drho, dee, eem
    real(rkind) :: drhof, deef
    real(rkind) :: sumnumrho,sumnumee,sumdenrho,sumdenee
    integer :: n,n2
!
    i = blockDim%x * (blockIdx%x - 1) + threadIdx%x
    k = blockDim%y * (blockIdx%y - 1) + threadIdx%y
    if (i > nx .or. k > nz) return
!
    do j=eul_jmin-1,eul_jmax
!
      ishk = 0
      do jj=j-weno_scheme+1,j+weno_scheme
        if (w_aux_gpu(i,jj,k,8) > sensor_threshold) ishk = 1
      enddo
!
      if (ishk == 0) then
        ft1  = 0._rkind
        ft2  = 0._rkind
        ft3  = 0._rkind
        ft4  = 0._rkind
        ft5  = 0._rkind
        ft6  = 0._rkind
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
!
              rhoi  = w_aux_gpu(i,j-m,k,1)
              uui   = w_aux_gpu(i,j-m,k,2)
              vvi   = w_aux_gpu(i,j-m,k,3)
              wwi   = w_aux_gpu(i,j-m,k,4)
              enti  = w_aux_gpu(i,j-m,k,5)
              tti   = w_aux_gpu(i,j-m,k,6)
              ppi   = tti*rhoi*rgas0
              eei   = enti-ppi/rhoi-0.5_rkind*(uui*uui+vvi*vvi+wwi*wwi)
!
              rhoip = w_aux_gpu(i,j-m+l,k,1)
              uuip  = w_aux_gpu(i,j-m+l,k,2)
              vvip  = w_aux_gpu(i,j-m+l,k,3)
              wwip  = w_aux_gpu(i,j-m+l,k,4)
              entip = w_aux_gpu(i,j-m+l,k,5)
              ttip  = w_aux_gpu(i,j-m+l,k,6)
              ppip  = ttip*rhoip*rgas0
              eeip  = entip-ppip/rhoip-0.5_rkind*(uuip*uuip+vvip*vvip+wwip*wwip)
!
              rhom  = rhoi + rhoip
              eem   = eei  + eeip
!
              drho   = 2._rkind*(rhoip-rhoi)/rhom
              dee    = 2._rkind*(eeip - eei)/eem
              sumnumrho = 1._rkind
              sumdenrho = 1._rkind
              sumnumee  = 1._rkind
              sumdenee  = 1._rkind
              do n = 1, nkeep
                n2 = 2*n
                sumdenrho = sumdenrho + (0.5_rkind*drho)**n2 / (1._rkind+n2)
                sumdenee  = sumdenee  + (0.5_rkind*dee )**n2
                sumnumee  = sumnumee  + (0.5_rkind*dee )**n2 / (1._rkind+n2)
              enddo
              drhof = sumnumrho/sumdenrho
              deef  = sumnumee /sumdenee
!
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
            ft1  = ft1  + coeff_deriv1_gpu(l,lmax)*uvs1
            ft2  = ft2  + coeff_deriv1_gpu(l,lmax)*uvs2
            ft3  = ft3  + coeff_deriv1_gpu(l,lmax)*uvs3
            ft4  = ft4  + coeff_deriv1_gpu(l,lmax)*uvs4
            ft5  = ft5  + coeff_deriv1_gpu(l,lmax)*(uvs5_i+uvs5_k+uvs5_p)
            ft6  = ft6  + coeff_deriv1_gpu(l,lmax)*uvs6
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
!
              rhoi  = w_aux_gpu(i,j-m,k,1)
              uui   = w_aux_gpu(i,j-m,k,2)
              vvi   = w_aux_gpu(i,j-m,k,3)
              wwi   = w_aux_gpu(i,j-m,k,4)
              enti  = w_aux_gpu(i,j-m,k,5)
              tti   = w_aux_gpu(i,j-m,k,6)
              ppi   = tti*rhoi*rgas0
!
              rhoip = w_aux_gpu(i,j-m+l,k,1)
              uuip  = w_aux_gpu(i,j-m+l,k,2)
              vvip  = w_aux_gpu(i,j-m+l,k,3)
              wwip  = w_aux_gpu(i,j-m+l,k,4)
              entip = w_aux_gpu(i,j-m+l,k,5)
              ttip  = w_aux_gpu(i,j-m+l,k,6)
              ppip  = ttip*rhoip*rgas0
!
              rhom  = rhoi + rhoip
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
!
        fh3 = fh3 + 0.5_rkind*ft6
!
        fhat_gpu(i,j,k,1) = fh1
        fhat_gpu(i,j,k,2) = fh2
        fhat_gpu(i,j,k,3) = fh3
        fhat_gpu(i,j,k,4) = fh4
        fhat_gpu(i,j,k,5) = fh5
      else
        evmax = -1._rkind
        do l=1,weno_size ! LLF
          ll = j + l - weno_scheme
          vv   = w_aux_gpu(i,ll,k,3)
          tt   = w_aux_gpu(i,ll,k,6)
          gamloc = get_gamloc_dev(indx_cp_l,indx_cp_r,tt,t0,cp_coeff_gpu,calorically_perfect,rgas0)
          c     = sqrt (gamloc*rgas0*tt)
          evm   = max(abs(vv-c),abs(vv+c))
          evmax = max(evm,evmax)
        enddo
        do l=1,weno_size ! loop over the stencil centered at face i
          ll = j + l - weno_scheme
!
          rho    = w_aux_gpu(i,ll,k,1)
          uu     = w_aux_gpu(i,ll,k,2)
          vv     = w_aux_gpu(i,ll,k,3)
          ww     = w_aux_gpu(i,ll,k,4)
          h      = w_aux_gpu(i,ll,k,5)
          rhov   = rho*vv
          pp     = rho*w_aux_gpu(i,ll,k,6)*rgas0
!
          rhoevm = rho*evmax
          evm = rhov
          c = 0.5_rkind * (evm + rhoevm)
          gplus_y_gpu (i,k,1,l) = c
          gminus_y_gpu(i,k,1,l) = evm-c
          evm = uu * rhov
          c = 0.5_rkind * (evm + rhoevm * uu)
          gplus_y_gpu (i,k,2,l) = c
          gminus_y_gpu(i,k,2,l) = evm-c
          evm = vv * rhov + pp
          c = 0.5_rkind * (evm + rhoevm * vv)
          gplus_y_gpu (i,k,3,l) = c
          gminus_y_gpu(i,k,3,l) = evm-c
          evm = ww * rhov
          c = 0.5_rkind * (evm + rhoevm * ww)
          gplus_y_gpu (i,k,4,l) = c
          gminus_y_gpu(i,k,4,l) = evm-c
          evm =  h  * rhov
          c = 0.5_rkind * (evm + evmax * (rho*h-pp))
          gplus_y_gpu (i,k,5,l) = c
          gminus_y_gpu(i,k,5,l) = evm-c
        enddo
!       
!       Reconstruction of the '+' and '-' fluxes
!       
        wenorec_ord = max(weno_scheme+ep_ord_change_gpu(i,j,k,2),1)
        call wenorec_rusanov(i,k,nv,nx,nz,gplus_y_gpu,gminus_y_gpu,fhat_gpu,&
          weno_scheme,wenorec_ord,weno_version,ng,nx,ny,nz,i,j,k)
!       
      endif
    enddo
!
!   Update net flux
    do j=eul_jmin,eul_jmax ! loop on the inner nodes
      do m=1,5
        fl_gpu(i,j,k,m) = fl_gpu(i,j,k,m) + (fhat_gpu(i,j,k,m)-fhat_gpu(i,j-1,k,m))*detady_gpu(j)
      enddo
    enddo
!
  endsubroutine euler_y_hybrid_rusanov_kernel
!
  attributes(device) subroutine wenorec_rusanov(ii,jj,nvar,n1,n2,vp,vm,vhat,&
    iweno,wenorec_ord,weno_version,ng,mx,my,mz,iii,jjj,kkk)
!
!   Passed arguments
    integer :: nvar, iweno, wenorec_ord, weno_version, ii, jj, n1, n2
    integer :: ng,mx,my,mz,iii,jjj,kkk
    real(rkind), dimension(n1,n2,nvar,2*iweno) :: vm,vp
    real(rkind), dimension(1-ng:mx+ng,1-ng:my+ng,1-ng:mz+ng,nvar) :: vhat
!
!   Local variables
    real(rkind), dimension(-1:4) :: dwe           ! linear weights
    real(rkind), dimension(-1:4) :: betap,betam   ! beta_l
    real(rkind) :: vminus, vplus
!   
    integer :: i,l,m
    real(rkind) :: c0,c1,c2,c3,c4,d0,d1,d2,d3,summ,sump
    real(rkind) :: tau5p,tau5m,eps40
!   
    if (wenorec_ord==1) then ! Godunov
!     
      i = iweno ! index of intermediate node to perform reconstruction
!     
      do m=1,nvar
        vminus  = vp(ii,jj,m,i)
        vplus   = vm(ii,jj,m,i+1)
        vhat(iii,jjj,kkk,m) = vminus+vplus
      enddo
!     
    elseif (wenorec_ord==2) then ! WENO-3
!     
      i = iweno ! index of intermediate node to perform reconstruction
!     
      dwe(1)   = 2._rkind/3._rkind
      dwe(0)   = 1._rkind/3._rkind
!     
      do m=1,nvar
!       
        betap(0)  = (vp(ii,jj,m,i  )-vp(ii,jj,m,i-1))**2
        betap(1)  = (vp(ii,jj,m,i+1)-vp(ii,jj,m,i  ))**2
!
        betam(0)  = (vm(ii,jj,m,i+2)-vm(ii,jj,m,i+1))**2
        betam(1)  = (vm(ii,jj,m,i+1)-vm(ii,jj,m,i  ))**2
!       
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
!       
        vminus = betap(0) *(-vp(ii,jj,m,i-1)+3*vp(ii,jj,m,i  )) + betap(1) *( vp(ii,jj,m,i  )+ vp(ii,jj,m,i+1))
        vplus  = betam(0) *(-vm(ii,jj,m,i+2)+3*vm(ii,jj,m,i+1)) + betam(1) *( vm(ii,jj,m,i  )+ vm(ii,jj,m,i+1))
        vhat(iii,jjj,kkk,m) = 0.5_rkind*(vminus+vplus)
!       
      enddo ! end of m-loop
!     
    elseif (wenorec_ord==3) then ! WENO-5
!     
      i = iweno ! index of intermediate node to perform reconstruction
!     
      dwe( 0) = 1._rkind/10._rkind
      dwe( 1) = 6._rkind/10._rkind
      dwe( 2) = 3._rkind/10._rkind
!     
!     JS
      d0 = 13._rkind/12._rkind
      d1 = 1._rkind/4._rkind
!     Weights for polynomial reconstructions
      c0 = 1._rkind/3._rkind
      c1 = 5._rkind/6._rkind
      c2 =-1._rkind/6._rkind
      c3 =-7._rkind/6._rkind
      c4 =11._rkind/6._rkind
!     
      if (weno_version==0) then ! Standard JS WENO 5
!       
        do m=1,nvar
!         
          betap(2) = d0*(vp(ii,jj,m,i)-2._rkind*vp(ii,jj,m,i+1)+vp(ii,jj,m,i+2))**2+&
          d1*(3._rkind*vp(ii,jj,m,i)-4._rkind*vp(ii,jj,m,i+1)+vp(ii,jj,m,i+2))**2
          betap(1) = d0*(vp(ii,jj,m,i-1)-2._rkind*vp(ii,jj,m,i)+vp(ii,jj,m,i+1))**2+&
          d1*(     vp(ii,jj,m,i-1)-vp(ii,jj,m,i+1) )**2
          betap(0) = d0*(vp(ii,jj,m,i)-2._rkind*vp(ii,jj,m,i-1)+vp(ii,jj,m,i-2))**2+&
          d1*(3._rkind*vp(ii,jj,m,i)-4._rkind*vp(ii,jj,m,i-1)+vp(ii,jj,m,i-2))**2
!         
          betam(2) = d0*(vm(ii,jj,m,i+1)-2._rkind*vm(ii,jj,m,i)+vm(ii,jj,m,i-1))**2+&
          d1*(3._rkind*vm(ii,jj,m,i+1)-4._rkind*vm(ii,jj,m,i)+vm(ii,jj,m,i-1))**2
          betam(1) = d0*(vm(ii,jj,m,i+2)-2._rkind*vm(ii,jj,m,i+1)+vm(ii,jj,m,i))**2+&
          d1*(     vm(ii,jj,m,i+2)-vm(ii,jj,m,i) )**2
          betam(0) = d0*(vm(ii,jj,m,i+1)-2._rkind*vm(ii,jj,m,i+2)+vm(ii,jj,m,i+3))**2+&
          d1*(3._rkind*vm(ii,jj,m,i+1)-4._rkind*vm(ii,jj,m,i+2)+vm(ii,jj,m,i+3))**2
!         
          sump = 0._rkind
          summ = 0._rkind
          do l=0,2
            betap(l) = dwe(  l)/(0.000001_rkind+betap(l))**2
            betam(l) = dwe(  l)/(0.000001_rkind+betam(l))**2
            sump = sump + betap(l)
            summ = summ + betam(l)
          enddo
          do l=0,2
            betap(l) = betap(l)/sump
            betam(l) = betam(l)/summ
          enddo
!         
          vminus = betap(2)*(c0*vp(ii,jj,m,i  )+c1*vp(ii,jj,m,i+1)+c2*vp(ii,jj,m,i+2)) + &
          betap(1)*(c2*vp(ii,jj,m,i-1)+c1*vp(ii,jj,m,i  )+c0*vp(ii,jj,m,i+1)) + &
          betap(0)*(c0*vp(ii,jj,m,i-2)+c3*vp(ii,jj,m,i-1)+c4*vp(ii,jj,m,i  ))
          vplus  = betam(2)*(c0*vm(ii,jj,m,i+1)+c1*vm(ii,jj,m,i  )+c2*vm(ii,jj,m,i-1)) + &
          betam(1)*(c2*vm(ii,jj,m,i+2)+c1*vm(ii,jj,m,i+1)+c0*vm(ii,jj,m,i  )) + &
          betam(0)*(c0*vm(ii,jj,m,i+3)+c3*vm(ii,jj,m,i+2)+c4*vm(ii,jj,m,i+1))
!         
          vhat(iii,jjj,kkk,m) = vminus+vplus
!         
        enddo ! end of m-loop
!       
        else ! WENO 5Z
!       
        do m=1,nvar
!         
          betap(2) = d0*(vp(ii,jj,m,i)-2._rkind*vp(ii,jj,m,i+1)+vp(ii,jj,m,i+2))**2+&
          d1*(3._rkind*vp(ii,jj,m,i)-4._rkind*vp(ii,jj,m,i+1)+vp(ii,jj,m,i+2))**2
          betap(1) = d0*(vp(ii,jj,m,i-1)-2._rkind*vp(ii,jj,m,i)+vp(ii,jj,m,i+1))**2+&
          d1*(     vp(ii,jj,m,i-1)-vp(ii,jj,m,i+1) )**2
          betap(0) = d0*(vp(ii,jj,m,i)-2._rkind*vp(ii,jj,m,i-1)+vp(ii,jj,m,i-2))**2+&
          d1*(3._rkind*vp(ii,jj,m,i)-4._rkind*vp(ii,jj,m,i-1)+vp(ii,jj,m,i-2))**2
!         
          betam(2) = d0*(vm(ii,jj,m,i+1)-2._rkind*vm(ii,jj,m,i)+vm(ii,jj,m,i-1))**2+&
          d1*(3._rkind*vm(ii,jj,m,i+1)-4._rkind*vm(ii,jj,m,i)+vm(ii,jj,m,i-1))**2
          betam(1) = d0*(vm(ii,jj,m,i+2)-2._rkind*vm(ii,jj,m,i+1)+vm(ii,jj,m,i))**2+&
          d1*(     vm(ii,jj,m,i+2)-vm(ii,jj,m,i) )**2
          betam(0) = d0*(vm(ii,jj,m,i+1)-2._rkind*vm(ii,jj,m,i+2)+vm(ii,jj,m,i+3))**2+&
          d1*(3._rkind*vm(ii,jj,m,i+1)-4._rkind*vm(ii,jj,m,i+2)+vm(ii,jj,m,i+3))**2
!         
          eps40 = 1.D-40
          tau5p = abs(betap(0)-betap(2))+eps40
          tau5m = abs(betam(0)-betam(2))+eps40
!         
          do l=0,2
            betap(l) = (betap(l)+eps40)/(betap(l)+tau5p)
            betam(l) = (betam(l)+eps40)/(betam(l)+tau5m)
          enddo
!         
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
!         
          vminus = betap(2)*(c0*vp(ii,jj,m,i  )+c1*vp(ii,jj,m,i+1)+c2*vp(ii,jj,m,i+2)) + &
          betap(1)*(c2*vp(ii,jj,m,i-1)+c1*vp(ii,jj,m,i  )+c0*vp(ii,jj,m,i+1)) + &
          betap(0)*(c0*vp(ii,jj,m,i-2)+c3*vp(ii,jj,m,i-1)+c4*vp(ii,jj,m,i  ))
          vplus  = betam(2)*(c0*vm(ii,jj,m,i+1)+c1*vm(ii,jj,m,i  )+c2*vm(ii,jj,m,i-1)) + &
          betam(1)*(c2*vm(ii,jj,m,i+2)+c1*vm(ii,jj,m,i+1)+c0*vm(ii,jj,m,i  )) + &
          betam(0)*(c0*vm(ii,jj,m,i+3)+c3*vm(ii,jj,m,i+2)+c4*vm(ii,jj,m,i+1))
!         
          vhat(iii,jjj,kkk,m) = vminus+vplus
!         
        enddo ! end of m-loop
!
      endif
!     
    elseif (wenorec_ord==4) then ! WENO-7
!     
      i = iweno ! index of intermediate node to perform reconstruction
!     
      dwe( 0) = 1._rkind/35._rkind
      dwe( 1) = 12._rkind/35._rkind
      dwe( 2) = 18._rkind/35._rkind
      dwe( 3) = 4._rkind/35._rkind
!     
!     JS weights
      d1 = 1._rkind/36._rkind
      d2 = 13._rkind/12._rkind
      d3 = 781._rkind/720._rkind
!     
      do m=1,nvar
!       
        betap(3)= d1*(-11*vp(ii,jj,m,  i)+18*vp(ii,jj,m,i+1)- 9*vp(ii,jj,m,i+2)+ 2*vp(ii,jj,m,i+3))**2+&
        d2*(  2*vp(ii,jj,m,  i)- 5*vp(ii,jj,m,i+1)+ 4*vp(ii,jj,m,i+2)-   vp(ii,jj,m,i+3))**2+ &
        d3*(   -vp(ii,jj,m,  i)+ 3*vp(ii,jj,m,i+1)- 3*vp(ii,jj,m,i+2)+   vp(ii,jj,m,i+3))**2
        betap(2)= d1*(- 2*vp(ii,jj,m,i-1)- 3*vp(ii,jj,m,i  )+ 6*vp(ii,jj,m,i+1)-   vp(ii,jj,m,i+2))**2+&
        d2*(    vp(ii,jj,m,i-1)- 2*vp(ii,jj,m,i  )+   vp(ii,jj,m,i+1)             )**2+&
        d3*(   -vp(ii,jj,m,i-1)+ 3*vp(ii,jj,m,i  )- 3*vp(ii,jj,m,i+1)+   vp(ii,jj,m,i+2))**2
        betap(1)= d1*(    vp(ii,jj,m,i-2)- 6*vp(ii,jj,m,i-1)+ 3*vp(ii,jj,m,i  )+ 2*vp(ii,jj,m,i+1))**2+&
        d2*( vp(ii,jj,m,i-1)- 2*vp(ii,jj,m,i  )+   vp(ii,jj,m,i+1))**2+ &
        d3*(   -vp(ii,jj,m,i-2)+ 3*vp(ii,jj,m,i-1)- 3*vp(ii,jj,m,i  )+   vp(ii,jj,m,i+1))**2
        betap(0)= d1*(- 2*vp(ii,jj,m,i-3)+ 9*vp(ii,jj,m,i-2)-18*vp(ii,jj,m,i-1)+11*vp(ii,jj,m,i  ))**2+&
        d2*(-   vp(ii,jj,m,i-3)+ 4*vp(ii,jj,m,i-2)- 5*vp(ii,jj,m,i-1)+ 2*vp(ii,jj,m,i  ))**2+&
        d3*(   -vp(ii,jj,m,i-3)+ 3*vp(ii,jj,m,i-2)- 3*vp(ii,jj,m,i-1)+   vp(ii,jj,m,i  ))**2
!       
        betam(3)= d1*(-11*vm(ii,jj,m,i+1)+18*vm(ii,jj,m,i  )- 9*vm(ii,jj,m,i-1)+ 2*vm(ii,jj,m,i-2))**2+&
        d2*(  2*vm(ii,jj,m,i+1)- 5*vm(ii,jj,m,i  )+ 4*vm(ii,jj,m,i-1)-   vm(ii,jj,m,i-2))**2+&
        d3*(   -vm(ii,jj,m,i+1)+ 3*vm(ii,jj,m,i  )- 3*vm(ii,jj,m,i-1)+   vm(ii,jj,m,i-2))**2
        betam(2)= d1*(- 2*vm(ii,jj,m,i+2)- 3*vm(ii,jj,m,i+1)+ 6*vm(ii,jj,m,i  )-   vm(ii,jj,m,i-1))**2+&
        d2*(    vm(ii,jj,m,i+2)- 2*vm(ii,jj,m,i+1)+   vm(ii,jj,m,i  )             )**2+&
        d3*(   -vm(ii,jj,m,i+2)+ 3*vm(ii,jj,m,i+1)- 3*vm(ii,jj,m,i  )+   vm(ii,jj,m,i-1))**2
        betam(1)= d1*(    vm(ii,jj,m,i+3)- 6*vm(ii,jj,m,i+2)+ 3*vm(ii,jj,m,i+1)+ 2*vm(ii,jj,m,i  ))**2+&
        d2*(                 vm(ii,jj,m,i+2)- 2*vm(ii,jj,m,i+1)+   vm(ii,jj,m,i  ))**2+&
        d3*(   -vm(ii,jj,m,i+3)+ 3*vm(ii,jj,m,i+2)- 3*vm(ii,jj,m,i+1)+   vm(ii,jj,m,i  ))**2
        betam(0)= d1*(- 2*vm(ii,jj,m,i+4)+ 9*vm(ii,jj,m,i+3)-18*vm(ii,jj,m,i+2)+11*vm(ii,jj,m,i+1))**2+&
        d2*(-   vm(ii,jj,m,i+4)+ 4*vm(ii,jj,m,i+3)- 5*vm(ii,jj,m,i+2)+ 2*vm(ii,jj,m,i+1))**2+&
        d3*(   -vm(ii,jj,m,i+4)+ 3*vm(ii,jj,m,i+3)- 3*vm(ii,jj,m,i+2)+   vm(ii,jj,m,i+1))**2
!       
        sump = 0._rkind
        summ = 0._rkind
        do l=0,3
          betap(l) = dwe(  l)/(0.000001_rkind+betap(l))**2
          betam(l) = dwe(  l)/(0.000001_rkind+betam(l))**2
          sump = sump + betap(l)
          summ = summ + betam(l)
        enddo
        do l=0,3
          betap(l) = betap(l)/sump
          betam(l) = betam(l)/summ
        enddo
!       
        vminus = betap(3)*( 6*vp(ii,jj,m,i  )+26*vp(ii,jj,m,i+1)-10*vp(ii,jj,m,i+2)+ 2*vp(ii,jj,m,i+3))+&
        betap(2)*(-2*vp(ii,jj,m,i-1)+14*vp(ii,jj,m,i  )+14*vp(ii,jj,m,i+1)- 2*vp(ii,jj,m,i+2))+&
        betap(1)*( 2*vp(ii,jj,m,i-2)-10*vp(ii,jj,m,i-1)+26*vp(ii,jj,m,i  )+ 6*vp(ii,jj,m,i+1))+&
        betap(0)*(-6*vp(ii,jj,m,i-3)+26*vp(ii,jj,m,i-2)-46*vp(ii,jj,m,i-1)+50*vp(ii,jj,m,i  ))
        vplus  =  betam(3)*( 6*vm(ii,jj,m,i+1)+26*vm(ii,jj,m,i  )-10*vm(ii,jj,m,i-1)+ 2*vm(ii,jj,m,i-2))+&
        betam(2)*(-2*vm(ii,jj,m,i+2)+14*vm(ii,jj,m,i+1)+14*vm(ii,jj,m,i  )- 2*vm(ii,jj,m,i-1))+&
        betam(1)*( 2*vm(ii,jj,m,i+3)-10*vm(ii,jj,m,i+2)+26*vm(ii,jj,m,i+1)+ 6*vm(ii,jj,m,i  ))+&
        betam(0)*(-6*vm(ii,jj,m,i+4)+26*vm(ii,jj,m,i+3)-46*vm(ii,jj,m,i+2)+50*vm(ii,jj,m,i+1))
!       
        vhat(iii,jjj,kkk,m) = (vminus+vplus)/24._rkind
!       
      enddo ! end of m-loop
!     
    else
      write(*,*) 'Error! WENO scheme not implemented'
      stop
    endif
!
  endsubroutine wenorec_rusanov
!
  attributes(device) subroutine wenorec(ii,jj,nvar,n1,n2,vp,vm,vhat,iweno,wenorec_ord,weno_version,rho0,u0)
!
!   Passed arguments
    integer :: nvar, iweno, wenorec_ord, weno_version, ii, jj, n1, n2
    real(rkind), dimension(n1,n2,nvar,2*iweno) :: vm,vp
    real(rkind), dimension(nvar) :: vhat
    real(rkind) :: rho0, u0
!
!   Local variables
    real(rkind), dimension(-1:4) :: dwe           ! linear weights
    real(rkind), dimension(-1:4) :: betap,betam   ! beta_l
    real(rkind), dimension(   5) :: betascale
    real(rkind) :: vminus, vplus
!   
    integer :: i,l,m
    real(rkind) :: c0,c1,c2,c3,c4,d0,d1,d2,d3,summ,sump
    real(rkind) :: tau5p,tau5m,eps40
    real(rkind) :: u0_2, rho0_2u0_2, rho0_2u0_4
!   
    u0_2       = u0*u0
    rho0_2u0_2 = rho0*rho0*u0_2
    rho0_2u0_4 = rho0_2u0_2*u0_2
    betascale(1) = 1._rkind/rho0_2u0_2
    betascale(2) = betascale(1)
    betascale(3) = betascale(1)
    betascale(4) = 1._rkind/rho0_2u0_4
    betascale(5) = betascale(4)
!   
    if (wenorec_ord==1) then ! Godunov
!     
      i = iweno ! index of intermediate node to perform reconstruction
!     
      do m=1,nvar
        vminus  = vp(ii,jj,m,i)
        vplus   = vm(ii,jj,m,i+1)
        vhat(m) = vminus+vplus
      enddo
!     
    elseif (wenorec_ord==2) then ! WENO-3
!     
      i = iweno ! index of intermediate node to perform reconstruction
!     
      dwe(1)   = 2._rkind/3._rkind
      dwe(0)   = 1._rkind/3._rkind
!     
      do m=1,nvar
!       
        betap(0)  = (vp(ii,jj,m,i  )-vp(ii,jj,m,i-1))**2
        betap(1)  = (vp(ii,jj,m,i+1)-vp(ii,jj,m,i  ))**2
        betap(1) = betascale(m)*betap(1)
        betap(0) = betascale(m)*betap(0)
!
        betam(0)  = (vm(ii,jj,m,i+2)-vm(ii,jj,m,i+1))**2
        betam(1)  = (vm(ii,jj,m,i+1)-vm(ii,jj,m,i  ))**2
        betam(1) = betascale(m)*betam(1)
        betam(0) = betascale(m)*betam(0)
!       
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
!       
        vminus = betap(0) *(-vp(ii,jj,m,i-1)+3*vp(ii,jj,m,i  )) + betap(1) *( vp(ii,jj,m,i  )+ vp(ii,jj,m,i+1))
        vplus  = betam(0) *(-vm(ii,jj,m,i+2)+3*vm(ii,jj,m,i+1)) + betam(1) *( vm(ii,jj,m,i  )+ vm(ii,jj,m,i+1))
        vhat(m) = 0.5_rkind*(vminus+vplus)
!       
      enddo ! end of m-loop
!     
    elseif (wenorec_ord==3) then ! WENO-5
!     
      i = iweno ! index of intermediate node to perform reconstruction
!     
      dwe( 0) = 1._rkind/10._rkind
      dwe( 1) = 6._rkind/10._rkind
      dwe( 2) = 3._rkind/10._rkind
!     
!     JS
      d0 = 13._rkind/12._rkind
      d1 = 1._rkind/4._rkind
!     Weights for polynomial reconstructions
      c0 = 1._rkind/3._rkind
      c1 = 5._rkind/6._rkind
      c2 =-1._rkind/6._rkind
      c3 =-7._rkind/6._rkind
      c4 =11._rkind/6._rkind
!     
      if (weno_version==0) then ! Standard JS WENO 5
!       
        do m=1,nvar
!         
          betap(2) = d0*(vp(ii,jj,m,i)-2._rkind*vp(ii,jj,m,i+1)+vp(ii,jj,m,i+2))**2+&
          d1*(3._rkind*vp(ii,jj,m,i)-4._rkind*vp(ii,jj,m,i+1)+vp(ii,jj,m,i+2))**2
          betap(1) = d0*(vp(ii,jj,m,i-1)-2._rkind*vp(ii,jj,m,i)+vp(ii,jj,m,i+1))**2+&
          d1*(     vp(ii,jj,m,i-1)-vp(ii,jj,m,i+1) )**2
          betap(0) = d0*(vp(ii,jj,m,i)-2._rkind*vp(ii,jj,m,i-1)+vp(ii,jj,m,i-2))**2+&
          d1*(3._rkind*vp(ii,jj,m,i)-4._rkind*vp(ii,jj,m,i-1)+vp(ii,jj,m,i-2))**2
          betap(2) = betascale(m)*betap(2)
          betap(1) = betascale(m)*betap(1)
          betap(0) = betascale(m)*betap(0)
!         
          betam(2) = d0*(vm(ii,jj,m,i+1)-2._rkind*vm(ii,jj,m,i)+vm(ii,jj,m,i-1))**2+&
          d1*(3._rkind*vm(ii,jj,m,i+1)-4._rkind*vm(ii,jj,m,i)+vm(ii,jj,m,i-1))**2
          betam(1) = d0*(vm(ii,jj,m,i+2)-2._rkind*vm(ii,jj,m,i+1)+vm(ii,jj,m,i))**2+&
          d1*(     vm(ii,jj,m,i+2)-vm(ii,jj,m,i) )**2
          betam(0) = d0*(vm(ii,jj,m,i+1)-2._rkind*vm(ii,jj,m,i+2)+vm(ii,jj,m,i+3))**2+&
          d1*(3._rkind*vm(ii,jj,m,i+1)-4._rkind*vm(ii,jj,m,i+2)+vm(ii,jj,m,i+3))**2
          betam(2) = betascale(m)*betam(2)
          betam(1) = betascale(m)*betam(1)
          betam(0) = betascale(m)*betam(0)
!         
          sump = 0._rkind
          summ = 0._rkind
          do l=0,2
            betap(l) = dwe(  l)/(0.000001_rkind+betap(l))**2
            betam(l) = dwe(  l)/(0.000001_rkind+betam(l))**2
            sump = sump + betap(l)
            summ = summ + betam(l)
          enddo
          do l=0,2
            betap(l) = betap(l)/sump
            betam(l) = betam(l)/summ
          enddo
!         
          vminus = betap(2)*(c0*vp(ii,jj,m,i  )+c1*vp(ii,jj,m,i+1)+c2*vp(ii,jj,m,i+2)) + &
          betap(1)*(c2*vp(ii,jj,m,i-1)+c1*vp(ii,jj,m,i  )+c0*vp(ii,jj,m,i+1)) + &
          betap(0)*(c0*vp(ii,jj,m,i-2)+c3*vp(ii,jj,m,i-1)+c4*vp(ii,jj,m,i  ))
          vplus  = betam(2)*(c0*vm(ii,jj,m,i+1)+c1*vm(ii,jj,m,i  )+c2*vm(ii,jj,m,i-1)) + &
          betam(1)*(c2*vm(ii,jj,m,i+2)+c1*vm(ii,jj,m,i+1)+c0*vm(ii,jj,m,i  )) + &
          betam(0)*(c0*vm(ii,jj,m,i+3)+c3*vm(ii,jj,m,i+2)+c4*vm(ii,jj,m,i+1))
!         
          vhat(m) = vminus+vplus
!         
        enddo ! end of m-loop
!       
        else ! WENO 5Z
!       
        do m=1,nvar
!         
          betap(2) = d0*(vp(ii,jj,m,i)-2._rkind*vp(ii,jj,m,i+1)+vp(ii,jj,m,i+2))**2+&
          d1*(3._rkind*vp(ii,jj,m,i)-4._rkind*vp(ii,jj,m,i+1)+vp(ii,jj,m,i+2))**2
          betap(1) = d0*(vp(ii,jj,m,i-1)-2._rkind*vp(ii,jj,m,i)+vp(ii,jj,m,i+1))**2+&
          d1*(     vp(ii,jj,m,i-1)-vp(ii,jj,m,i+1) )**2
          betap(0) = d0*(vp(ii,jj,m,i)-2._rkind*vp(ii,jj,m,i-1)+vp(ii,jj,m,i-2))**2+&
          d1*(3._rkind*vp(ii,jj,m,i)-4._rkind*vp(ii,jj,m,i-1)+vp(ii,jj,m,i-2))**2
          betap(2) = betascale(m)*betap(2)
          betap(1) = betascale(m)*betap(1)
          betap(0) = betascale(m)*betap(0)
!         
          betam(2) = d0*(vm(ii,jj,m,i+1)-2._rkind*vm(ii,jj,m,i)+vm(ii,jj,m,i-1))**2+&
          d1*(3._rkind*vm(ii,jj,m,i+1)-4._rkind*vm(ii,jj,m,i)+vm(ii,jj,m,i-1))**2
          betam(1) = d0*(vm(ii,jj,m,i+2)-2._rkind*vm(ii,jj,m,i+1)+vm(ii,jj,m,i))**2+&
          d1*(     vm(ii,jj,m,i+2)-vm(ii,jj,m,i) )**2
          betam(0) = d0*(vm(ii,jj,m,i+1)-2._rkind*vm(ii,jj,m,i+2)+vm(ii,jj,m,i+3))**2+&
          d1*(3._rkind*vm(ii,jj,m,i+1)-4._rkind*vm(ii,jj,m,i+2)+vm(ii,jj,m,i+3))**2
          betam(2) = betascale(m)*betam(2)
          betam(1) = betascale(m)*betam(1)
          betam(0) = betascale(m)*betam(0)
!         
          eps40 = 1.D-40
          tau5p = abs(betap(0)-betap(2))+eps40
          tau5m = abs(betam(0)-betam(2))+eps40
!         
          do l=0,2
            betap(l) = (betap(l)+eps40)/(betap(l)+tau5p)
            betam(l) = (betam(l)+eps40)/(betam(l)+tau5m)
          enddo
!         
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
!         
          vminus = betap(2)*(c0*vp(ii,jj,m,i  )+c1*vp(ii,jj,m,i+1)+c2*vp(ii,jj,m,i+2)) + &
          betap(1)*(c2*vp(ii,jj,m,i-1)+c1*vp(ii,jj,m,i  )+c0*vp(ii,jj,m,i+1)) + &
          betap(0)*(c0*vp(ii,jj,m,i-2)+c3*vp(ii,jj,m,i-1)+c4*vp(ii,jj,m,i  ))
          vplus  = betam(2)*(c0*vm(ii,jj,m,i+1)+c1*vm(ii,jj,m,i  )+c2*vm(ii,jj,m,i-1)) + &
          betam(1)*(c2*vm(ii,jj,m,i+2)+c1*vm(ii,jj,m,i+1)+c0*vm(ii,jj,m,i  )) + &
          betam(0)*(c0*vm(ii,jj,m,i+3)+c3*vm(ii,jj,m,i+2)+c4*vm(ii,jj,m,i+1))
!         
          vhat(m) = vminus+vplus
!         
        enddo ! end of m-loop
!
      endif
!     
    elseif (wenorec_ord==4) then ! WENO-7
!     
      i = iweno ! index of intermediate node to perform reconstruction
!     
      dwe( 0) = 1._rkind/35._rkind
      dwe( 1) = 12._rkind/35._rkind
      dwe( 2) = 18._rkind/35._rkind
      dwe( 3) = 4._rkind/35._rkind
!     
!     JS weights
      d1 = 1._rkind/36._rkind
      d2 = 13._rkind/12._rkind
      d3 = 781._rkind/720._rkind
!     
      do m=1,nvar
!       
        betap(3)= d1*(-11*vp(ii,jj,m,  i)+18*vp(ii,jj,m,i+1)- 9*vp(ii,jj,m,i+2)+ 2*vp(ii,jj,m,i+3))**2+&
        d2*(  2*vp(ii,jj,m,  i)- 5*vp(ii,jj,m,i+1)+ 4*vp(ii,jj,m,i+2)-   vp(ii,jj,m,i+3))**2+ &
        d3*(   -vp(ii,jj,m,  i)+ 3*vp(ii,jj,m,i+1)- 3*vp(ii,jj,m,i+2)+   vp(ii,jj,m,i+3))**2
        betap(2)= d1*(- 2*vp(ii,jj,m,i-1)- 3*vp(ii,jj,m,i  )+ 6*vp(ii,jj,m,i+1)-   vp(ii,jj,m,i+2))**2+&
        d2*(    vp(ii,jj,m,i-1)- 2*vp(ii,jj,m,i  )+   vp(ii,jj,m,i+1)             )**2+&
        d3*(   -vp(ii,jj,m,i-1)+ 3*vp(ii,jj,m,i  )- 3*vp(ii,jj,m,i+1)+   vp(ii,jj,m,i+2))**2
        betap(1)= d1*(    vp(ii,jj,m,i-2)- 6*vp(ii,jj,m,i-1)+ 3*vp(ii,jj,m,i  )+ 2*vp(ii,jj,m,i+1))**2+&
        d2*( vp(ii,jj,m,i-1)- 2*vp(ii,jj,m,i  )+   vp(ii,jj,m,i+1))**2+ &
        d3*(   -vp(ii,jj,m,i-2)+ 3*vp(ii,jj,m,i-1)- 3*vp(ii,jj,m,i  )+   vp(ii,jj,m,i+1))**2
        betap(0)= d1*(- 2*vp(ii,jj,m,i-3)+ 9*vp(ii,jj,m,i-2)-18*vp(ii,jj,m,i-1)+11*vp(ii,jj,m,i  ))**2+&
        d2*(-   vp(ii,jj,m,i-3)+ 4*vp(ii,jj,m,i-2)- 5*vp(ii,jj,m,i-1)+ 2*vp(ii,jj,m,i  ))**2+&
        d3*(   -vp(ii,jj,m,i-3)+ 3*vp(ii,jj,m,i-2)- 3*vp(ii,jj,m,i-1)+   vp(ii,jj,m,i  ))**2
        betap(3) = betascale(m)*betap(3)
        betap(2) = betascale(m)*betap(2)
        betap(1) = betascale(m)*betap(1)
        betap(0) = betascale(m)*betap(0)
!       
        betam(3)= d1*(-11*vm(ii,jj,m,i+1)+18*vm(ii,jj,m,i  )- 9*vm(ii,jj,m,i-1)+ 2*vm(ii,jj,m,i-2))**2+&
        d2*(  2*vm(ii,jj,m,i+1)- 5*vm(ii,jj,m,i  )+ 4*vm(ii,jj,m,i-1)-   vm(ii,jj,m,i-2))**2+&
        d3*(   -vm(ii,jj,m,i+1)+ 3*vm(ii,jj,m,i  )- 3*vm(ii,jj,m,i-1)+   vm(ii,jj,m,i-2))**2
        betam(2)= d1*(- 2*vm(ii,jj,m,i+2)- 3*vm(ii,jj,m,i+1)+ 6*vm(ii,jj,m,i  )-   vm(ii,jj,m,i-1))**2+&
        d2*(    vm(ii,jj,m,i+2)- 2*vm(ii,jj,m,i+1)+   vm(ii,jj,m,i  )             )**2+&
        d3*(   -vm(ii,jj,m,i+2)+ 3*vm(ii,jj,m,i+1)- 3*vm(ii,jj,m,i  )+   vm(ii,jj,m,i-1))**2
        betam(1)= d1*(    vm(ii,jj,m,i+3)- 6*vm(ii,jj,m,i+2)+ 3*vm(ii,jj,m,i+1)+ 2*vm(ii,jj,m,i  ))**2+&
        d2*(                 vm(ii,jj,m,i+2)- 2*vm(ii,jj,m,i+1)+   vm(ii,jj,m,i  ))**2+&
        d3*(   -vm(ii,jj,m,i+3)+ 3*vm(ii,jj,m,i+2)- 3*vm(ii,jj,m,i+1)+   vm(ii,jj,m,i  ))**2
        betam(0)= d1*(- 2*vm(ii,jj,m,i+4)+ 9*vm(ii,jj,m,i+3)-18*vm(ii,jj,m,i+2)+11*vm(ii,jj,m,i+1))**2+&
        d2*(-   vm(ii,jj,m,i+4)+ 4*vm(ii,jj,m,i+3)- 5*vm(ii,jj,m,i+2)+ 2*vm(ii,jj,m,i+1))**2+&
        d3*(   -vm(ii,jj,m,i+4)+ 3*vm(ii,jj,m,i+3)- 3*vm(ii,jj,m,i+2)+   vm(ii,jj,m,i+1))**2
        betam(3) = betascale(m)*betam(3)
        betam(2) = betascale(m)*betam(2)
        betam(1) = betascale(m)*betam(1)
        betam(0) = betascale(m)*betam(0)
!       
        sump = 0._rkind
        summ = 0._rkind
        do l=0,3
          betap(l) = dwe(  l)/(0.000001_rkind+betap(l))**2
          betam(l) = dwe(  l)/(0.000001_rkind+betam(l))**2
          sump = sump + betap(l)
          summ = summ + betam(l)
        enddo
        do l=0,3
          betap(l) = betap(l)/sump
          betam(l) = betam(l)/summ
        enddo
!       
        vminus = betap(3)*( 6*vp(ii,jj,m,i  )+26*vp(ii,jj,m,i+1)-10*vp(ii,jj,m,i+2)+ 2*vp(ii,jj,m,i+3))+&
        betap(2)*(-2*vp(ii,jj,m,i-1)+14*vp(ii,jj,m,i  )+14*vp(ii,jj,m,i+1)- 2*vp(ii,jj,m,i+2))+&
        betap(1)*( 2*vp(ii,jj,m,i-2)-10*vp(ii,jj,m,i-1)+26*vp(ii,jj,m,i  )+ 6*vp(ii,jj,m,i+1))+&
        betap(0)*(-6*vp(ii,jj,m,i-3)+26*vp(ii,jj,m,i-2)-46*vp(ii,jj,m,i-1)+50*vp(ii,jj,m,i  ))
        vplus  =  betam(3)*( 6*vm(ii,jj,m,i+1)+26*vm(ii,jj,m,i  )-10*vm(ii,jj,m,i-1)+ 2*vm(ii,jj,m,i-2))+&
        betam(2)*(-2*vm(ii,jj,m,i+2)+14*vm(ii,jj,m,i+1)+14*vm(ii,jj,m,i  )- 2*vm(ii,jj,m,i-1))+&
        betam(1)*( 2*vm(ii,jj,m,i+3)-10*vm(ii,jj,m,i+2)+26*vm(ii,jj,m,i+1)+ 6*vm(ii,jj,m,i  ))+&
        betam(0)*(-6*vm(ii,jj,m,i+4)+26*vm(ii,jj,m,i+3)-46*vm(ii,jj,m,i+2)+50*vm(ii,jj,m,i+1))
!       
        vhat(m) = (vminus+vplus)/24._rkind
!       
      enddo ! end of m-loop
!     
    else
      write(*,*) 'Error! WENO scheme not implemented'
      stop
    endif
!
  endsubroutine wenorec
! 
  attributes(device) subroutine eigenvectors_x(b1, b2, b3, uu, vv, ww, c, ci, h, el, er)
    real(rkind), intent(in) :: b1, b2, b3, uu, vv, ww, c, ci, h
    real(rkind), dimension(:,:), intent(out) ::  el, er
!
!   left eigenvectors matrix (at Roe state)
!   matrix L-1 AIAA2001-2609 - Rohde after transposition
    el(1,1)   =   0.5_rkind * (b1     + uu * ci)
    el(2,1)   =  -0.5_rkind * (b2 * uu +     ci)
    el(3,1)   =  -0.5_rkind * (b2 * vv         )
    el(4,1)   =  -0.5_rkind * (b2 * ww         )
    el(5,1)   =   0.5_rkind * b2
    el(1,2)   =   1._rkind - b1
    el(2,2)   =   b2*uu
    el(3,2)   =   b2*vv
    el(4,2)   =   b2*ww
    el(5,2)   =  -b2
    el(1,3)   =   0.5_rkind * (b1     - uu * ci)
    el(2,3)   =  -0.5_rkind * (b2 * uu -     ci)
    el(3,3)   =  -0.5_rkind * (b2 * vv         )
    el(4,3)   =  -0.5_rkind * (b2 * ww         )
    el(5,3)   =   0.5_rkind * b2
    el(1,4)   =   -vv ! vv
    el(2,4)   =   0._rkind
    el(3,4)   =   1._rkind ! -1._rkind
    el(4,4)   =   0._rkind
    el(5,4)   =   0._rkind
    el(1,5)   =  -ww
    el(2,5)   =   0._rkind
    el(3,5)   =   0._rkind
    el(4,5)   =   1._rkind
    el(5,5)   =   0._rkind
!
!   right eigenvectors matrix (at Roe state)
!   matrix R-1 AIAA2001-2609 - Rohde after transposition
    er(1,1)   =  1._rkind
    er(2,1)   =  1._rkind
    er(3,1)   =  1._rkind
    er(4,1)   =  0._rkind
    er(5,1)   =  0._rkind
    er(1,2)   =  uu -  c
    er(2,2)   =  uu
    er(3,2)   =  uu +  c
    er(4,2)   =  0._rkind
    er(5,2)   =  0._rkind
    er(1,3)   =  vv
    er(2,3)   =  vv
    er(3,3)   =  vv
    er(4,3)   =  1._rkind ! -1._rkind
    er(5,3)   =   0._rkind
    er(1,4)   =  ww
    er(2,4)   =  ww
    er(3,4)   =  ww
    er(4,4)   =  0._rkind
    er(5,4)   =  1._rkind
    er(1,5)   =  h  - uu * c
    er(2,5)   =  b3 ! etot - rho * p_rho/p_e
    er(3,5)   =  h  + uu * c
    er(4,5)   =  vv ! -vv
    er(5,5)   =  ww
!
  endsubroutine eigenvectors_x
!
  attributes(device) subroutine eigenvectors_y(b1, b2, b3, uu, vv, ww, c, ci, h, el, er)
    real(rkind), intent(in) :: b1, b2, b3, uu, vv, ww, c, ci, h
    real(rkind), dimension(:,:), intent(out) ::  el, er
!
!   left eigenvectors matrix (at Roe state)
!   matrix L-2 AIAA2001-2609 - Rohde after transposition
    el(1,1) =   0.5_rkind * (b1     + vv * ci)
    el(2,1) =  -0.5_rkind * (b2 * uu         )
    el(3,1) =  -0.5_rkind * (b2 * vv +     ci)
    el(4,1) =  -0.5_rkind * (b2 * ww         )
    el(5,1) =   0.5_rkind * b2
    el(1,2) =   1._rkind - b1
    el(2,2) =   b2 * uu
    el(3,2) =   b2 * vv
    el(4,2) =   b2 * ww
    el(5,2) =  -b2
    el(1,3) =   0.5_rkind * (b1     - vv * ci)
    el(2,3) =  -0.5_rkind * (b2 * uu         )
    el(3,3) =  -0.5_rkind * (b2 * vv -     ci)
    el(4,3) =  -0.5_rkind * (b2 * ww         )
    el(5,3) =   0.5_rkind * b2
    el(1,4) =  -ww ! -uu
    el(2,4) =   0._rkind !  1._rkind
    el(3,4) =   0._rkind
    el(4,4) =   1._rkind !  0._rkind
    el(5,4) =   0._rkind
    el(1,5) =   uu ! ww
    el(2,5) =  -1._rkind !  0._rkind
    el(3,5) =   0._rkind
    el(4,5) =   0._rkind ! -1._rkind
    el(5,5) =   0._rkind
!
!   right eigenvectors matrix (at Roe state)
!   matrix R-2 AIAA2001-2609 - Rohde after transposition
    er(1,1) =  1._rkind
    er(2,1) =  1._rkind
    er(3,1) =  1._rkind
    er(4,1) =  0._rkind
    er(5,1) =  0._rkind
    er(1,2) =  uu
    er(2,2) =  uu
    er(3,2) =  uu
    er(4,2) =  0._rkind !  1._rkind
    er(5,2) = -1._rkind !  0._rkind
    er(1,3) =  vv - c
    er(2,3) =  vv
    er(3,3) =  vv + c
    er(4,3) =  0._rkind
    er(5,3) =  0._rkind
    er(1,4) =  ww
    er(2,4) =  ww
    er(3,4) =  ww
    er(4,4) =  1._rkind  !  0._rkind
    er(5,4) =  0._rkind ! -1._rkind
    er(1,5) =  h  - vv * c
    er(2,5) =  b3 ! etot - rho * p_rho/p_e !qq
    er(3,5) =  h  + vv * c
    er(4,5) =  ww     ! uu
    er(5,5) =  -uu    ! -ww
  endsubroutine eigenvectors_y
!
  attributes(device) subroutine eigenvectors_z(b1, b2, b3, uu, vv, ww, c, ci, h, el, er)
    real(rkind), intent(in) :: b1, b2, b3, uu, vv, ww, c, ci, h
    real(rkind), dimension(:,:), intent(out) ::  el, er
!
!   left eigenvectors matrix (at Roe state)
!   matrix L-3 AIAA2001-2609 - Rohde after transposition
    el(1,1)   =   0.5_rkind * (b1     + ww * ci)
    el(2,1)   =  -0.5_rkind * (b2 * uu         )
    el(3,1)   =  -0.5_rkind * (b2 * vv         )
    el(4,1)   =  -0.5_rkind * (b2 * ww +     ci)
    el(5,1)   =   0.5_rkind * b2
    el(1,2)   =   1._rkind - b1
    el(2,2)   =   b2 * uu
    el(3,2)   =   b2 * vv
    el(4,2)   =   b2 * ww
    el(5,2)   =  -b2
    el(1,3)   =   0.5_rkind * (b1     - ww * ci)
    el(2,3)   =  -0.5_rkind * (b2 * uu         )
    el(3,3)   =  -0.5_rkind * (b2 * vv         )
    el(4,3)   =  -0.5_rkind * (b2 * ww -     ci)
    el(5,3)   =   0.5_rkind * b2
    el(1,4)   =  -uu ! uu
    el(2,4)   =   1._rkind ! -1._rkind
    el(3,4)   =   0._rkind
    el(4,4)   =   0._rkind
    el(5,4)   =   0._rkind
    el(1,5)   =  -vv
    el(2,5)   =   0._rkind
    el(3,5)   =   1._rkind
    el(4,5)   =   0._rkind
    el(5,5)   =   0._rkind
!
!   right eigenvectors matrix (at Roe state)
!   matrix R-3 AIAA2001-2609 - Rohde after transposition
    er(1,1)   =  1._rkind
    er(2,1)   =  1._rkind
    er(3,1)   =  1._rkind
    er(4,1)   =  0._rkind
    er(5,1)   =  0._rkind
    er(1,2)   =  uu
    er(2,2)   =  uu
    er(3,2)   =  uu
    er(4,2)   =  1._rkind ! -1._rkind
    er(5,2)   =  0._rkind
    er(1,3)   =  vv
    er(2,3)   =  vv
    er(3,3)   =  vv
    er(4,3)   =  0._rkind
    er(5,3)   =  1._rkind
    er(1,4)   =  ww - c
    er(2,4)   =  ww
    er(3,4)   =  ww + c
    er(4,4)   =  0._rkind
    er(5,4)   =  0._rkind
    er(1,5)   =  h  - ww * c
    er(2,5)   =  b3 ! etot - rho * p_rho/p_e
    er(3,5)   =  h  + ww * c
    er(4,5)   =  uu ! -uu
    er(5,5)   =  vv
!
  endsubroutine eigenvectors_z
!
  attributes(device) subroutine compute_roe_average(nx, ny, nz, ng, i, ip, j, jp, k, kp, w_aux_gpu, rgas0, &
    b1, b2, b3, c, ci, h, uu, vv, ww, cp_coeff_gpu, indx_cp_l, indx_cp_r,&
    calorically_perfect,tol_iter_nr,t0)
    integer :: ng,i,ip,j,jp,k,kp,indx_cp_l,indx_cp_r,calorically_perfect
    integer :: nx,ny,nz
    real(rkind), intent(in) :: rgas0, tol_iter_nr,t0
    real(rkind), intent(out) :: b1, b2, b3, uu, vv, ww, ci, h, c
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in) :: w_aux_gpu
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in) :: cp_coeff_gpu
    real(rkind) :: up, vp, wp, qqp, hp, r, rp1, cc, qq, gam, gm1
    integer :: ll,iter,max_iter
    real(rkind) :: tt,gm1loc,hbar,ttp,told,num,den,gamloc,cploc,tpow,tpowp,p_rho,p_e,etot,rho,cp0,cv0
!   
    max_iter = 50
!   
    cp0 = cp_coeff_gpu(0)
    cv0 = cp_coeff_gpu(0)-rgas0
    gam = cp0/cv0
    gm1 = gam-1._rkind
!   Compute Roe average
!   Left state (node i)
!   ri        =  1._rkind/w_aux_gpu(im,jm,km,1)
    uu        =  w_aux_gpu(i,j,k,2)
    vv        =  w_aux_gpu(i,j,k,3)
    ww        =  w_aux_gpu(i,j,k,4)
    qq        =  0.5_rkind * (uu*uu  +vv*vv + ww*ww)
!   pp        =  gm1 * (w_gpu(i,j,k,5) - w_gpu(i,j,k,1) * qq)
!   h         =  (w_gpu(i,j,k,5)  + pp) * ri
    h         =  w_aux_gpu(i,j,k,5) ! gamma*w_gpu(im,jm,km,5)*ri-gm1*qq
!   Right state (node i+1)
!   rip       =  1._rkind/w_aux_gpu(ip,jp,kp,1)
    up        =  w_aux_gpu(ip,jp,kp,2)
    vp        =  w_aux_gpu(ip,jp,kp,3)
    wp        =  w_aux_gpu(ip,jp,kp,4)
    qqp       =  0.5_rkind * (up*up  +vp*vp +wp*wp)
!   ppp       =  gm1 * (w_gpu(i,jp,k,5) - w_gpu(i,jp,k,1) * qqp)
!   hp        =  (w_gpu(i,jp,k,5)  + ppp) * rip
    hp        =  w_aux_gpu(ip,jp,kp,5) ! gamma*w_gpu(ip,jp,kp,5)*rip-gm1*qqp
!   average state
    r         =  w_aux_gpu(ip,jp,kp,1)/w_aux_gpu(i,j,k,1)
    r         =  sqrt(r)
    rho       =  r*w_aux_gpu(i,j,k,1)
    rp1       =  1._rkind/(r  +1._rkind)
    uu        =  (r*up  +uu)*rp1
    vv        =  (r*vp  +vv)*rp1
    ww        =  (r*wp  +ww)*rp1
    h         =  (r*hp  +h)*rp1
    qq        =  0.5_rkind * (uu*uu  +vv*vv + ww*ww)
!   
    if (calorically_perfect==1) then
!     cc       = gm1 * (h-qq+cv0*t0) ! cv0*t0 needed because e = cv*(tt-t0)
      cc       = gm1 * (h-qq       )
      tt       = cc/gam/rgas0
      gm1loc   = gm1
    else
      hbar      = h - qq - cp_coeff_gpu(indx_cp_r+1)*t0
      tt        = w_aux_gpu(i ,j ,k ,6)
      ttp       = w_aux_gpu(ip,jp,kp,6)
      tt        = (r*ttp +tt)*rp1
      told      = tt ! First attempt computed as Roe average of temperature
      do iter=1,max_iter
        num = 0._rkind
        den = 0._rkind
        do ll=indx_cp_l,indx_cp_r
          if (ll==-1) then
            tpow  = (told/t0)**ll
            den = den+cp_coeff_gpu(ll)*tpow
            num = num+cp_coeff_gpu(ll)*log(told/t0)
          else
            tpow  = (told/t0)**ll
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
      gm1loc = gamloc-1._rkind
      cc     = gamloc*tt*rgas0
    endif
!   
    c         =  sqrt(cc)
    ci        =  1._rkind/c
!   
    p_rho     = tt*rgas0
    p_e       = rho*gm1loc   ! p_e   = rho * R_gas / Cv
    etot      = h - tt*rgas0
!   
    b3        = etot - rho * p_rho/p_e ! qq in the case of calorically perfect
    b2        = p_e/(rho*cc)
    b1        = p_rho/cc - b2*(etot - 2._rkind*qq)
!   
!   b2        = gm1/cc  ! 1/(cp*T)
!   b1        = b2 * qq
!   b3        = qq
!   b2        = gm1/cc
!   b1        = b2 * b3
!
  endsubroutine compute_roe_average
!
  subroutine force_rhs_2_cuf(nx, ny, nz, ng, fln_gpu, w_aux_gpu, bulk5g, fluid_mask_gpu)
    integer :: nx, ny, nz, ng
    real(rkind), dimension(5), intent(in) :: bulk5g
    real(rkind), dimension(1:,1:,1:,1:), intent(inout), device :: fln_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in), device :: w_aux_gpu
    integer, dimension(1-ng:,1-ng:,1-ng:), intent(in), device :: fluid_mask_gpu
    real(rkind) :: bulk_1, bulk_2, uu
    integer :: i,j,k,iercuda
    bulk_1 = bulk5g(1)
    bulk_2 = bulk5g(2)
    !$cuf kernel do(3) <<<*,*>>>
    do k=1,nz
      do j=1,ny
        do i=1,nx
          if (fluid_mask_gpu(i,j,k)==0) then
            uu = w_aux_gpu(i,j,k,2)
            fln_gpu(i,j,k,1) = fln_gpu(i,j,k,1) -    bulk_1
            fln_gpu(i,j,k,2) = fln_gpu(i,j,k,2) -    bulk_2
            fln_gpu(i,j,k,5) = fln_gpu(i,j,k,5) - uu*bulk_2
          endif
        enddo
      enddo
    enddo
    !@cuf iercuda=cudaDeviceSynchronize()
  endsubroutine force_rhs_2_cuf
!
  subroutine force_rhs_1_cuf(nx, ny, nz, ng, yn_gpu, fln_gpu, w_gpu, w_aux_gpu, bulk5, fluid_mask_gpu)
    integer, intent(in) :: nx, ny, nz, ng
    real(rkind), dimension(1:), intent(in), device :: yn_gpu
    real(rkind), dimension(1:,1:,1:,1:), intent(in), device :: fln_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in), device :: w_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in), device :: w_aux_gpu
    integer, dimension(1-ng:,1-ng:,1-ng:), intent(in), device :: fluid_mask_gpu
    real(rkind), dimension(5), intent(out) :: bulk5
    real(rkind) :: bulk_1, bulk_2, bulk_3, bulk_4, bulk_5
    real(rkind) :: dy
    integer :: i,j,k,iercuda
    bulk_1 = 0._rkind
    bulk_2 = 0._rkind
    bulk_3 = 0._rkind
    bulk_4 = 0._rkind
    bulk_5 = 0._rkind
!
    !$cuf kernel do(3) <<<*,*>>> reduce(+:bulk_1,bulk_2,bulk_3,bulk_4,bulk_5)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          if (fluid_mask_gpu(i,j,k)==0) then
            dy = yn_gpu(j+1)-yn_gpu(j)
            bulk_1 = bulk_1 + fln_gpu(i,j,k,1)*dy
            bulk_2 = bulk_2 + fln_gpu(i,j,k,2)*dy
            bulk_3 = bulk_3 + w_gpu  (i,j,k,1)*dy
            bulk_4 = bulk_4 + w_gpu  (i,j,k,2)*dy
            bulk_5 = bulk_5 + w_gpu  (i,j,k,2)*dy*w_aux_gpu(i,j,k,6)
          endif
        enddo
      enddo
    enddo
    !@cuf iercuda=cudaDeviceSynchronize()
!
    bulk5(1) = bulk_1
    bulk5(2) = bulk_2
    bulk5(3) = bulk_3
    bulk5(4) = bulk_4
    bulk5(5) = bulk_5
  endsubroutine force_rhs_1_cuf
!
  subroutine force_var_1_cuf(nx, ny, nz, ng, yn_gpu, fln_gpu, w_gpu, w_aux_gpu, bulkt, fluid_mask_gpu, cv_coeff_gpu, &
    indx_cp_l, indx_cp_r, t0, calorically_perfect, tol_iter_nr)
    integer, intent(in) :: nx, ny, nz, ng, indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind), intent(in) :: t0, tol_iter_nr
    real(rkind), dimension(1:), intent(in), device :: yn_gpu
    real(rkind), dimension(1:,1:,1:,1:), intent(in), device :: fln_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in), device :: w_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(inout), device :: w_aux_gpu
    integer, dimension(1-ng:,1-ng:,1-ng:), intent(in), device :: fluid_mask_gpu
    real(rkind),    dimension(indx_cp_l:indx_cp_r+1), intent(in),  device :: cv_coeff_gpu
    real(rkind), intent(out) :: bulkt
    real(rkind) :: bulk_5,rho,rhou,rhov,rhow,rhoe,ri,uu,vv,ww,qq,ee,tt
    real(rkind) :: dy
    integer :: i,j,k,iercuda
!   
    bulk_5 = 0._rkind
!   
    !$cuf kernel do(3) <<<*,*>>> reduce(+:bulk_5)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          if (fluid_mask_gpu(i,j,k)==0) then
            dy = yn_gpu(j+1)-yn_gpu(j)
            rho  = w_gpu(i,j,k,1)
            rhou = w_gpu(i,j,k,2)
            rhov = w_gpu(i,j,k,3)
            rhow = w_gpu(i,j,k,4)
            rhoe = w_gpu(i,j,k,5)
            ri   = 1._rkind/rho
            uu   = rhou*ri
            vv   = rhov*ri
            ww   = rhow*ri
            qq   = 0.5_rkind*(uu*uu+vv*vv+ww*ww)
            ee = rhoe/rho-qq
            tt = get_temperature_from_e_dev(ee,w_aux_gpu(i,j,k,6),t0,cv_coeff_gpu,indx_cp_l,indx_cp_r, &
            calorically_perfect,tol_iter_nr)
            w_aux_gpu(i,j,k,6) = tt
            bulk_5 = bulk_5 + rhou*tt*dy
          endif
        enddo
      enddo
    enddo
    !@cuf iercuda=cudaDeviceSynchronize()
!   
    bulkt = bulk_5
!   
  endsubroutine force_var_1_cuf
! 
  subroutine force_var_2_cuf(nx, ny, nz, ng, w_gpu, w_aux_gpu, tbdiff, fluid_mask_gpu, cv_coeff_gpu, indx_cp_l, indx_cp_r, &
    t0, calorically_perfect)
    integer :: nx, ny, nz, ng, indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind), intent(in) :: tbdiff, t0
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(inout), device :: w_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in), device :: w_aux_gpu
    integer, dimension(1-ng:,1-ng:,1-ng:), intent(in), device :: fluid_mask_gpu
    real(rkind),    dimension(indx_cp_l:indx_cp_r+1), intent(in),  device :: cv_coeff_gpu
    real(rkind) :: rho,rhou,rhov,rhow,tt,ttnew,ee
    integer :: i,j,k,iercuda
!   
    !$cuf kernel do(3) <<<*,*>>>
    do k=1,nz
      do j=1,ny
        do i=1,nx
          if (fluid_mask_gpu(i,j,k)==0) then
            rho  = w_gpu(i,j,k,1)
            rhou = w_gpu(i,j,k,2)
            rhov = w_gpu(i,j,k,3)
            rhow = w_gpu(i,j,k,4)
            tt   = w_aux_gpu(i,j,k,6)
            ttnew = tt+tbdiff
            ee = get_e_from_temperature_dev(ttnew, t0, indx_cp_l, indx_cp_r, cv_coeff_gpu, calorically_perfect)
            w_gpu(i,j,k,5) = rho*ee + 0.5_rkind*(rhou**2+rhov**2+rhow**2)/rho
          endif
        enddo
      enddo
    enddo
    !@cuf iercuda=cudaDeviceSynchronize()
!   
  endsubroutine force_var_2_cuf
!
  subroutine update_flux_cuf(nx, ny, nz, nv, fl_gpu, fln_gpu, gamdt)
    integer :: nx, ny, nz, nv
    real(rkind) :: gamdt
    real(rkind), dimension(1:,1:,1:,1:), intent(in), device :: fl_gpu
    real(rkind), dimension(1:,1:,1:,1:), intent(inout), device :: fln_gpu
    integer :: i,j,k,m,iercuda
    !$cuf kernel do(3) <<<*,*>>>
    do k=1,nz
      do j=1,ny
        do i=1,nx
          do m=1,nv
            fln_gpu(i,j,k,m) = fln_gpu(i,j,k,m)-gamdt*fl_gpu(i,j,k,m)
          enddo
        enddo
      enddo
    enddo
    !@cuf iercuda=cudaDeviceSynchronize()
  endsubroutine update_flux_cuf
!
  subroutine update_field_cuf(nx, ny, nz, ng, nv, w_gpu, fln_gpu, fluid_mask_gpu)
    integer :: nx, ny, nz, nv, ng
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(inout), device :: w_gpu
    real(rkind), dimension(1:,1:,1:,1:), intent(in), device :: fln_gpu
    integer, dimension(1-ng:,1-ng:,1-ng:), intent(in), device :: fluid_mask_gpu
    integer :: i,j,k,m,iercuda
    !$cuf kernel do(3) <<<*,*>>>
    do k=1,nz
      do j=1,ny
        do i=1,nx
          if (fluid_mask_gpu(i,j,k)==0) then
            do m=1,nv
              w_gpu(i,j,k,m) = w_gpu(i,j,k,m)+fln_gpu(i,j,k,m)
            enddo
          endif
        enddo
      enddo
    enddo
    !@cuf iercuda=cudaDeviceSynchronize()
  endsubroutine update_field_cuf
!
  subroutine visflx_div_ord2_cuf(nx, ny, nz, ng, &
    w_aux_gpu, fl_gpu, &
    x_gpu, y_gpu, z_gpu, stream_id)
!
    integer, intent(in) :: nx, ny, nz, ng
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(inout), device :: w_aux_gpu
    real(rkind), dimension(1:,1:,1:, 1:), intent(inout), device :: fl_gpu
    real(rkind), dimension(1-ng:), intent(in), device :: x_gpu, y_gpu, z_gpu
    integer     :: i,j,k,iercuda
    real(rkind) :: dxl,dyl,dzl
    real(rkind) :: uu,vv,ww,mu
    real(rkind) :: sigq,sigx,sigy,sigz
    integer(kind=cuda_stream_kind), intent(in) :: stream_id
!
    !$cuf kernel do(3) <<<*,*,stream=stream_id>>>
    do k=1,nz
      do j=1,ny
        do i=1,nx
!
          uu = w_aux_gpu(i,j,k,2)
          vv = w_aux_gpu(i,j,k,3)
          ww = w_aux_gpu(i,j,k,4)
          mu = w_aux_gpu(i,j,k,7)
!
          dxl  = mu/(x_gpu(i+1)-x_gpu(i-1))
          dyl  = mu/(y_gpu(j+1)-y_gpu(j-1))
          dzl  = mu/(z_gpu(k+1)-z_gpu(k-1))
          sigx = dxl*(w_aux_gpu(i+1,j,k,10)-w_aux_gpu(i-1,j,k,10))
          sigy = dyl*(w_aux_gpu(i,j+1,k,10)-w_aux_gpu(i,j-1,k,10))
          sigz = dzl*(w_aux_gpu(i,j,k+1,10)-w_aux_gpu(i,j,k-1,10))
          sigq = sigx*uu+sigy*vv+sigz*ww
!
          fl_gpu(i,j,k,2) = fl_gpu(i,j,k,2) - sigx
          fl_gpu(i,j,k,3) = fl_gpu(i,j,k,3) - sigy
          fl_gpu(i,j,k,4) = fl_gpu(i,j,k,4) - sigz
          fl_gpu(i,j,k,5) = fl_gpu(i,j,k,5) - sigq
!
        enddo
      enddo
    enddo
!
  endsubroutine visflx_div_ord2_cuf
!
  subroutine visflx_div_cuf(nx, ny, nz, ng, visc_order, &
    w_aux_gpu, fl_gpu, coeff_deriv1_gpu, &
    dcsidx_gpu, detady_gpu, dzitdz_gpu, stream_id)
!
    integer, intent(in) :: nx, ny, nz, ng, visc_order
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(inout), device :: w_aux_gpu
    real(rkind), dimension(1:,1:,1:, 1:), intent(inout), device :: fl_gpu
    real(rkind), dimension(1:,1:), intent(in), device :: coeff_deriv1_gpu
    real(rkind), dimension(1:), intent(in), device :: dcsidx_gpu, detady_gpu, dzitdz_gpu
    integer     :: i,j,k,l,iercuda
    real(rkind) :: ccl
    real(rkind) :: uu,vv,ww,tt,mu
    real(rkind) :: sigq,sigx,sigy,sigz
    real(rkind) :: divx3l, divy3l, divz3l
    integer     :: lmax
    integer(kind=cuda_stream_kind), intent(in) :: stream_id
!
    lmax = visc_order/2
!
    !$cuf kernel do(3) <<<*,*,stream=stream_id>>>
    do k=1,nz
      do j=1,ny
        do i=1,nx
!
          uu = w_aux_gpu(i,j,k,2)
          vv = w_aux_gpu(i,j,k,3)
          ww = w_aux_gpu(i,j,k,4)
          tt = w_aux_gpu(i,j,k,6)
          mu = w_aux_gpu(i,j,k,7)
!
          divx3l = 0._rkind
          divy3l = 0._rkind
          divz3l = 0._rkind
!
          do l=1,lmax
            ccl = coeff_deriv1_gpu(l,lmax)
            divx3l = divx3l+ccl*(w_aux_gpu(i+l,j,k,10)-w_aux_gpu(i-l,j,k,10))
            divy3l = divy3l+ccl*(w_aux_gpu(i,j+l,k,10)-w_aux_gpu(i,j-l,k,10))
            divz3l = divz3l+ccl*(w_aux_gpu(i,j,k+l,10)-w_aux_gpu(i,j,k-l,10))
          enddo
!
          divx3l = divx3l*dcsidx_gpu(i)
          divy3l = divy3l*detady_gpu(j)
          divz3l = divz3l*dzitdz_gpu(k)
!
          sigx  = mu*divx3l
          sigy  = mu*divy3l
          sigz  = mu*divz3l
          sigq  = sigx*uu+sigy*vv+sigz*ww
!
          fl_gpu(i,j,k,2) = fl_gpu(i,j,k,2) - sigx
          fl_gpu(i,j,k,3) = fl_gpu(i,j,k,3) - sigy
          fl_gpu(i,j,k,4) = fl_gpu(i,j,k,4) - sigz
          fl_gpu(i,j,k,5) = fl_gpu(i,j,k,5) - sigq
!
        enddo
      enddo
    enddo
!
  endsubroutine visflx_div_cuf
!
  subroutine visflx_cuf(nx, ny, nz, ng, visc_order, &
    Prandtl, t0, indx_cp_l, indx_cp_r, cp_coeff_gpu, calorically_perfect, &
    u0, l0, w_gpu, w_aux_gpu, fl_gpu, &
    coeff_deriv1_gpu, coeff_deriv2_gpu, &
    dcsidx_gpu, detady_gpu, dzitdz_gpu,  &
    dcsidxs_gpu, detadys_gpu, dzitdzs_gpu, &
    dcsidx2_gpu, detady2_gpu, dzitdz2_gpu, wallprop_gpu)
!
    integer, intent(in) :: nx, ny, nz, ng, visc_order, calorically_perfect
    integer, intent(in) :: indx_cp_l, indx_cp_r
    real(rkind), intent(in) :: Prandtl, u0, l0, t0
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(in), device :: w_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(inout), device :: w_aux_gpu
    real(rkind), dimension(1-ng:,1-ng:, 2:), intent(inout), device :: wallprop_gpu
    real(rkind), dimension(1:,1:,1:, 1:), intent(inout), device :: fl_gpu
    real(rkind), dimension(1:,1:), intent(in), device :: coeff_deriv1_gpu
    real(rkind), dimension(0:,1:), intent(in), device  :: coeff_deriv2_gpu
    real(rkind), dimension(1:), intent(in), device :: dcsidx_gpu, detady_gpu, dzitdz_gpu
    real(rkind), dimension(1:), intent(in), device :: dcsidxs_gpu, detadys_gpu, dzitdzs_gpu
    real(rkind), dimension(1:), intent(in), device :: dcsidx2_gpu, detady2_gpu, dzitdz2_gpu
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in), device :: cp_coeff_gpu
    integer     :: i,j,k,l,ll,iercuda
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
    integer     :: lmax
!
    lmax = visc_order/2
!
    !$cuf kernel do(2) <<<*,*>>>
    do j=1,ny
      do i=1,nx
        do k=1,nz
!
          uu = w_aux_gpu(i,j,k,2)
          vv = w_aux_gpu(i,j,k,3)
          ww = w_aux_gpu(i,j,k,4)
          tt = w_aux_gpu(i,j,k,6)
          mu = w_aux_gpu(i,j,k,7)
!
          if (calorically_perfect==1) then
            cploc = cp_coeff_gpu(0)
          else
            cploc = 0._rkind
            do ll=indx_cp_l,indx_cp_r
              cploc = cploc+cp_coeff_gpu(ll)*(tt/t0)**ll
            enddo
          endif
!
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
!
          do l=1,lmax
            ccl = coeff_deriv1_gpu(l,lmax)
            ux = ux+ccl*(w_aux_gpu(i+l,j,k,2)-w_aux_gpu(i-l,j,k,2))
            vx = vx+ccl*(w_aux_gpu(i+l,j,k,3)-w_aux_gpu(i-l,j,k,3))
            wx = wx+ccl*(w_aux_gpu(i+l,j,k,4)-w_aux_gpu(i-l,j,k,4))
            tx = tx+ccl*(w_aux_gpu(i+l,j,k,6)-w_aux_gpu(i-l,j,k,6))
            mux = mux+ccl*(w_aux_gpu(i+l,j,k,7)-w_aux_gpu(i-l,j,k,7))
!
            uy = uy+ccl*(w_aux_gpu(i,j+l,k,2)-w_aux_gpu(i,j-l,k,2))
            vy = vy+ccl*(w_aux_gpu(i,j+l,k,3)-w_aux_gpu(i,j-l,k,3))
            wy = wy+ccl*(w_aux_gpu(i,j+l,k,4)-w_aux_gpu(i,j-l,k,4))
            ty = ty+ccl*(w_aux_gpu(i,j+l,k,6)-w_aux_gpu(i,j-l,k,6))
            muy = muy+ccl*(w_aux_gpu(i,j+l,k,7)-w_aux_gpu(i,j-l,k,7))
!
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
!         
          if (j==1) then
            wallprop_gpu(i,k,2) = mu*uy
            wallprop_gpu(i,k,3) = mu*wy
            wallprop_gpu(i,k,4) = mu*ty*cploc/Prandtl
          endif
!         
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
!
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
!
          ulap  = ulapx+ulapy+ulapz
          vlap  = vlapx+vlapy+vlapz
          wlap  = wlapx+wlapy+wlapz
          tlap  = tlapx+tlapy+tlapz
!
          div     = ux+vy+wz
          div3l   = div/3._rkind
          w_aux_gpu(i,j,k,10) = div3l
          omegax = wy-vz
          omegay = uz-wx
          omegaz = vx-uy
          omod2 = omegax*omegax+omegay*omegay+omegaz*omegaz
          w_aux_gpu(i,j,k,9) = sqrt(omod2)
!
!         Original Ducros shock sensor
!         w_aux_gpu(i,j,k,8) = div2/(omod2+div2+1.D-12)
!         
!         Modified Ducros shock sensor
!         w_aux_gpu(i,j,k,8) = div2/(omod2+div2+(u0/l0)**2)
!         
!         Modified Ducros shock sensor (remove expansions)
          w_aux_gpu(i,j,k,8) = (max(-div/sqrt(omod2+div**2+(u0/l0)**2),0._rkind))**2
!         
!         ducfilter = div2/(div2+0.1_rkind*sqrt(uu*uu+vv*vv+ww*ww)*max(dcsidx_gpu(i),detady_gpu(j),dzitdz_gpu(k))+1.D-12)
!         w_aux_gpu(i,j,k,8) = (max(-div/sqrt(omod2+div**2+1.D-12),0._rkind))**2*ducfilter
!
          sig11 = 2._rkind*(ux-div3l)
          sig12 = uy+vx
          sig13 = uz+wx
          sig22 = 2._rkind*(vy-div3l)
          sig23 = vz+wy
          sig33 = 2._rkind*(wz-div3l)
          sigx  = mux*sig11 + muy*sig12 + muz*sig13 + mu*ulap
          sigy  = mux*sig12 + muy*sig22 + muz*sig23 + mu*vlap
          sigz  = mux*sig13 + muy*sig23 + muz*sig33 + mu*wlap
          sigqt = (mux*tx+muy*ty+muz*tz+mu*tlap)*cploc/Prandtl ! Conduction
!
          sigah = (sig11*ux+sig12*uy+sig13*uz+sig12*vx+sig22*vy+sig23*vz+sig13*wx+sig23*wy+sig33*wz)*mu ! Aerodynamic heating
!
          sigq  = sigx*uu+sigy*vv+sigz*ww+sigah+sigqt
!
          fl_gpu(i,j,k,2) = fl_gpu(i,j,k,2) - sigx
          fl_gpu(i,j,k,3) = fl_gpu(i,j,k,3) - sigy
          fl_gpu(i,j,k,4) = fl_gpu(i,j,k,4) - sigz
          fl_gpu(i,j,k,5) = fl_gpu(i,j,k,5) - sigq
!
        enddo
      enddo
    enddo
    !@cuf iercuda=cudaDeviceSynchronize()
  endsubroutine visflx_cuf
!
  subroutine visflx_nosensor_cuf(nx, ny, nz, ng, visc_order, &
    Prandtl, t0, indx_cp_l, indx_cp_r, cp_coeff_gpu, calorically_perfect, &
    u0, l0, w_gpu, w_aux_gpu, fl_gpu, &
    coeff_deriv1_gpu, coeff_deriv2_gpu, &
    dcsidx_gpu, detady_gpu, dzitdz_gpu,  &
    dcsidxs_gpu, detadys_gpu, dzitdzs_gpu, &
    dcsidx2_gpu, detady2_gpu, dzitdz2_gpu, wallprop_gpu)
!
    integer, intent(in) :: nx, ny, nz, ng, visc_order, calorically_perfect
    integer, intent(in) :: indx_cp_l, indx_cp_r
    real(rkind), intent(in) :: Prandtl, u0, l0, t0
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(in), device :: w_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(inout), device :: w_aux_gpu
    real(rkind), dimension(1-ng:,1-ng:, 2:), intent(inout), device :: wallprop_gpu
    real(rkind), dimension(1:,1:,1:, 1:), intent(inout), device :: fl_gpu
    real(rkind), dimension(1:,1:), intent(in), device :: coeff_deriv1_gpu
    real(rkind), dimension(0:,1:), intent(in), device  :: coeff_deriv2_gpu
    real(rkind), dimension(1:), intent(in), device :: dcsidx_gpu, detady_gpu, dzitdz_gpu
    real(rkind), dimension(1:), intent(in), device :: dcsidxs_gpu, detadys_gpu, dzitdzs_gpu
    real(rkind), dimension(1:), intent(in), device :: dcsidx2_gpu, detady2_gpu, dzitdz2_gpu
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in), device :: cp_coeff_gpu
    integer     :: i,j,k,l,ll,iercuda
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
    integer     :: lmax
!
    lmax = visc_order/2
!
    !$cuf kernel do(2) <<<*,*>>>
    do j=1,ny
      do i=1,nx
        do k=1,nz
!
          uu = w_aux_gpu(i,j,k,2)
          vv = w_aux_gpu(i,j,k,3)
          ww = w_aux_gpu(i,j,k,4)
          tt = w_aux_gpu(i,j,k,6)
          mu = w_aux_gpu(i,j,k,7)
!
          if (calorically_perfect==1) then
            cploc = cp_coeff_gpu(0)
          else
            cploc = 0._rkind
            do ll=indx_cp_l,indx_cp_r
              cploc = cploc+cp_coeff_gpu(ll)*(tt/t0)**ll
            enddo
          endif
!
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
!
          do l=1,lmax
            ccl = coeff_deriv1_gpu(l,lmax)
            ux = ux+ccl*(w_aux_gpu(i+l,j,k,2)-w_aux_gpu(i-l,j,k,2))
            vx = vx+ccl*(w_aux_gpu(i+l,j,k,3)-w_aux_gpu(i-l,j,k,3))
            wx = wx+ccl*(w_aux_gpu(i+l,j,k,4)-w_aux_gpu(i-l,j,k,4))
            tx = tx+ccl*(w_aux_gpu(i+l,j,k,6)-w_aux_gpu(i-l,j,k,6))
            mux = mux+ccl*(w_aux_gpu(i+l,j,k,7)-w_aux_gpu(i-l,j,k,7))
!
            uy = uy+ccl*(w_aux_gpu(i,j+l,k,2)-w_aux_gpu(i,j-l,k,2))
            vy = vy+ccl*(w_aux_gpu(i,j+l,k,3)-w_aux_gpu(i,j-l,k,3))
            wy = wy+ccl*(w_aux_gpu(i,j+l,k,4)-w_aux_gpu(i,j-l,k,4))
            ty = ty+ccl*(w_aux_gpu(i,j+l,k,6)-w_aux_gpu(i,j-l,k,6))
            muy = muy+ccl*(w_aux_gpu(i,j+l,k,7)-w_aux_gpu(i,j-l,k,7))
!
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
!         
          if (j==1) then
            wallprop_gpu(i,k,2) = mu*uy
            wallprop_gpu(i,k,3) = mu*wy
            wallprop_gpu(i,k,4) = mu*ty*cploc/Prandtl
          endif
!         
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
!
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
!
          ulap  = ulapx+ulapy+ulapz
          vlap  = vlapx+vlapy+vlapz
          wlap  = wlapx+wlapy+wlapz
          tlap  = tlapx+tlapy+tlapz
!
          div     = ux+vy+wz
          div3l   = div/3._rkind
          w_aux_gpu(i,j,k,10) = div3l
          omegax = wy-vz
          omegay = uz-wx
          omegaz = vx-uy
          omod2 = omegax*omegax+omegay*omegay+omegaz*omegaz
          w_aux_gpu(i,j,k,9) = sqrt(omod2)
!
          sig11 = 2._rkind*(ux-div3l)
          sig12 = uy+vx
          sig13 = uz+wx
          sig22 = 2._rkind*(vy-div3l)
          sig23 = vz+wy
          sig33 = 2._rkind*(wz-div3l)
          sigx  = mux*sig11 + muy*sig12 + muz*sig13 + mu*ulap
          sigy  = mux*sig12 + muy*sig22 + muz*sig23 + mu*vlap
          sigz  = mux*sig13 + muy*sig23 + muz*sig33 + mu*wlap
          sigqt = (mux*tx+muy*ty+muz*tz+mu*tlap)*cploc/Prandtl ! Conduction
!
          sigah = (sig11*ux+sig12*uy+sig13*uz+sig12*vx+sig22*vy+sig23*vz+sig13*wx+sig23*wy+sig33*wz)*mu ! Aerodynamic heating
!
          sigq  = sigx*uu+sigy*vv+sigz*ww+sigah+sigqt
!
          fl_gpu(i,j,k,2) = fl_gpu(i,j,k,2) - sigx
          fl_gpu(i,j,k,3) = fl_gpu(i,j,k,3) - sigy
          fl_gpu(i,j,k,4) = fl_gpu(i,j,k,4) - sigz
          fl_gpu(i,j,k,5) = fl_gpu(i,j,k,5) - sigq
!
        enddo
      enddo
    enddo
    !@cuf iercuda=cudaDeviceSynchronize()
  endsubroutine visflx_nosensor_cuf
!
  subroutine visflx_reduced_ord2_cuf(nx, ny, nz, ng, &
    Prandtl, t0, indx_cp_l, indx_cp_r, cp_coeff_gpu, calorically_perfect, &
    u0, l0, w_gpu, w_aux_gpu, fl_gpu, &
    x_gpu, y_gpu, z_gpu, wallprop_gpu, update_sensor)
!
    integer, intent(in) :: nx, ny, nz, ng, calorically_perfect
    integer, intent(in) :: indx_cp_l, indx_cp_r
    real(rkind), intent(in) :: Prandtl, u0, l0, t0
    integer, intent(in) :: update_sensor
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(in), device :: w_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(inout), device :: w_aux_gpu
    real(rkind), dimension(1-ng:,1-ng:, 2:), intent(inout), device :: wallprop_gpu
    real(rkind), dimension(1:,1:,1:, 1:), intent(inout), device :: fl_gpu
    real(rkind), dimension(1-ng:), intent(in), device :: x_gpu, y_gpu, z_gpu
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in), device :: cp_coeff_gpu
    integer     :: i,j,k,l,ll,iercuda
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
!
    !$cuf kernel do(3) <<<*,*>>>
    do k=1,nz
      do j=1,ny
        do i=1,nx
!
          uu = w_aux_gpu(i,j,k,2)
          vv = w_aux_gpu(i,j,k,3)
          ww = w_aux_gpu(i,j,k,4)
          tt = w_aux_gpu(i,j,k,6)
          mu = w_aux_gpu(i,j,k,7)
!
          if (calorically_perfect==1) then
            cploc = cp_coeff_gpu(0)
          else
            cploc = 0._rkind
            do ll=indx_cp_l,indx_cp_r
              cploc = cploc+cp_coeff_gpu(ll)*(tt/t0)**ll
            enddo
          endif
!
          dxl = 1._rkind/(x_gpu(i+1)-x_gpu(i-1))
          dyl = 1._rkind/(y_gpu(j+1)-y_gpu(j-1))
          dzl = 1._rkind/(z_gpu(k+1)-z_gpu(k-1))
!
          ux  = dxl*(w_aux_gpu(i+1,j,k,2)-w_aux_gpu(i-1,j,k,2))
          vx  = dxl*(w_aux_gpu(i+1,j,k,3)-w_aux_gpu(i-1,j,k,3))
          wx  = dxl*(w_aux_gpu(i+1,j,k,4)-w_aux_gpu(i-1,j,k,4))
          mux = dxl*(w_aux_gpu(i+1,j,k,7)-w_aux_gpu(i-1,j,k,7))
!
          uy  = dyl*(w_aux_gpu(i,j+1,k,2)-w_aux_gpu(i,j-1,k,2))
          vy  = dyl*(w_aux_gpu(i,j+1,k,3)-w_aux_gpu(i,j-1,k,3))
          wy  = dyl*(w_aux_gpu(i,j+1,k,4)-w_aux_gpu(i,j-1,k,4))
          ty  = dyl*(w_aux_gpu(i,j+1,k,6)-w_aux_gpu(i,j-1,k,6))
          muy = dyl*(w_aux_gpu(i,j+1,k,7)-w_aux_gpu(i,j-1,k,7))
!
          uz  = dzl*(w_aux_gpu(i,j,k+1,2)-w_aux_gpu(i,j,k-1,2))
          vz  = dzl*(w_aux_gpu(i,j,k+1,3)-w_aux_gpu(i,j,k-1,3))
          wz  = dzl*(w_aux_gpu(i,j,k+1,4)-w_aux_gpu(i,j,k-1,4))
          muz = dzl*(w_aux_gpu(i,j,k+1,7)-w_aux_gpu(i,j,k-1,7))
!         
          if (j==1) then
            wallprop_gpu(i,k,2) = mu*uy
            wallprop_gpu(i,k,3) = mu*wy
            wallprop_gpu(i,k,4) = mu*ty*cploc/Prandtl
          endif
!         
          div     = ux+vy+wz
          div3l   = div/3._rkind
          w_aux_gpu(i,j,k,10) = div3l
          omegax = wy-vz
          omegay = uz-wx
          omegaz = vx-uy
          omod2 = omegax*omegax+omegay*omegay+omegaz*omegaz
          w_aux_gpu(i,j,k,9) = sqrt(omod2)
!
          if (update_sensor == 1) then
!           Original Ducros shock sensor
!           w_aux_gpu(i,j,k,8) = div2/(omod2+div2+1.D-12)
!           
!           Modified Ducros shock sensor
!           w_aux_gpu(i,j,k,8) = div2/(omod2+div2+(u0/l0)**2)
!           
!           Modified Ducros shock sensor (remove expansions)
            w_aux_gpu(i,j,k,8) = (max(-div/sqrt(omod2+div**2+(u0/l0)**2),0._rkind))**2
!           
!           ducfilter = div2/(div2+0.1_rkind*sqrt(uu*uu+vv*vv+ww*ww)*max(dcsidx_gpu(i),detady_gpu(j),dzitdz_gpu(k))+1.D-12)
!           w_aux_gpu(i,j,k,8) = (max(-div/sqrt(omod2+div**2+1.D-12),0._rkind))**2*ducfilter
          endif
!
          sig11 = ux-2._rkind*div3l
          sig12 = vx
          sig13 = wx
          sig21 = uy
          sig22 = vy-2._rkind*div3l
          sig23 = wy
          sig31 = uz
          sig32 = vz
          sig33 = wz-2._rkind*div3l
          sigx  = mux*sig11 + muy*sig12 + muz*sig13
          sigy  = mux*sig12 + muy*sig22 + muz*sig23
          sigz  = mux*sig13 + muy*sig23 + muz*sig33
!         
          sigah = (sig11*ux+sig12*uy+sig13*uz+sig12*vx+sig22*vy+sig23*vz+sig13*wx+sig23*wy+sig33*wz)*mu ! Aerodynamic heating
!         
          sigq  = sigx*uu+sigy*vv+sigz*ww+sigah
!         
          fl_gpu(i,j,k,2) = fl_gpu(i,j,k,2) - sigx
          fl_gpu(i,j,k,3) = fl_gpu(i,j,k,3) - sigy
          fl_gpu(i,j,k,4) = fl_gpu(i,j,k,4) - sigz
          fl_gpu(i,j,k,5) = fl_gpu(i,j,k,5) - sigq
!         
        enddo
      enddo
    enddo
    !@cuf iercuda=cudaDeviceSynchronize()
  endsubroutine visflx_reduced_ord2_cuf
!
  subroutine visflx_reduced_cuf(nx, ny, nz, ng, visc_order, &
    Prandtl, t0, indx_cp_l, indx_cp_r, cp_coeff_gpu, calorically_perfect, &
    u0, l0, w_gpu, w_aux_gpu, fl_gpu, &
    coeff_deriv1_gpu, coeff_deriv2_gpu, &
    dcsidx_gpu, detady_gpu, dzitdz_gpu,  &
    dcsidxs_gpu, detadys_gpu, dzitdzs_gpu, &
    dcsidx2_gpu, detady2_gpu, dzitdz2_gpu, wallprop_gpu, update_sensor)
!
    integer, intent(in) :: nx, ny, nz, ng, visc_order, calorically_perfect
    integer, intent(in) :: indx_cp_l, indx_cp_r
    real(rkind), intent(in) :: Prandtl, u0, l0, t0
    integer, intent(in) :: update_sensor
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(in), device :: w_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(inout), device :: w_aux_gpu
    real(rkind), dimension(1-ng:,1-ng:, 2:), intent(inout), device :: wallprop_gpu
    real(rkind), dimension(1:,1:,1:, 1:), intent(inout), device :: fl_gpu
    real(rkind), dimension(1:,1:), intent(in), device :: coeff_deriv1_gpu
    real(rkind), dimension(0:,1:), intent(in), device  :: coeff_deriv2_gpu
    real(rkind), dimension(1:), intent(in), device :: dcsidx_gpu, detady_gpu, dzitdz_gpu
    real(rkind), dimension(1:), intent(in), device :: dcsidxs_gpu, detadys_gpu, dzitdzs_gpu
    real(rkind), dimension(1:), intent(in), device :: dcsidx2_gpu, detady2_gpu, dzitdz2_gpu
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in), device :: cp_coeff_gpu
    integer     :: i,j,k,l,ll,iercuda
    real(rkind) :: ccl,clapl
    real(rkind) :: sig11,sig12,sig13
    real(rkind) :: sig21,sig22,sig23
    real(rkind) :: sig31,sig32,sig33
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
    integer     :: lmax
!
    lmax = visc_order/2
!
    !$cuf kernel do(3) <<<*,*>>>
    do k=1,nz
      do j=1,ny
        do i=1,nx
!
          uu = w_aux_gpu(i,j,k,2)
          vv = w_aux_gpu(i,j,k,3)
          ww = w_aux_gpu(i,j,k,4)
          tt = w_aux_gpu(i,j,k,6)
          mu = w_aux_gpu(i,j,k,7)
!
          if (calorically_perfect==1) then
            cploc = cp_coeff_gpu(0)
          else
            cploc = 0._rkind
            do ll=indx_cp_l,indx_cp_r
              cploc = cploc+cp_coeff_gpu(ll)*(tt/t0)**ll
            enddo
          endif
!
          ux = 0._rkind
          vx = 0._rkind
          wx = 0._rkind
!         tx = 0._rkind
          mux = 0._rkind
          uy = 0._rkind
          vy = 0._rkind
          wy = 0._rkind
          ty = 0._rkind
          muy = 0._rkind
          uz = 0._rkind
          vz = 0._rkind
          wz = 0._rkind
!         tz = 0._rkind
          muz = 0._rkind
!
          do l=1,lmax
            ccl = coeff_deriv1_gpu(l,lmax)
            ux = ux+ccl*(w_aux_gpu(i+l,j,k,2)-w_aux_gpu(i-l,j,k,2))
            vx = vx+ccl*(w_aux_gpu(i+l,j,k,3)-w_aux_gpu(i-l,j,k,3))
            wx = wx+ccl*(w_aux_gpu(i+l,j,k,4)-w_aux_gpu(i-l,j,k,4))
!           tx = tx+ccl*(w_aux_gpu(i+l,j,k,6)-w_aux_gpu(i-l,j,k,6))
            mux = mux+ccl*(w_aux_gpu(i+l,j,k,7)-w_aux_gpu(i-l,j,k,7))
!
            uy = uy+ccl*(w_aux_gpu(i,j+l,k,2)-w_aux_gpu(i,j-l,k,2))
            vy = vy+ccl*(w_aux_gpu(i,j+l,k,3)-w_aux_gpu(i,j-l,k,3))
            wy = wy+ccl*(w_aux_gpu(i,j+l,k,4)-w_aux_gpu(i,j-l,k,4))
            ty = ty+ccl*(w_aux_gpu(i,j+l,k,6)-w_aux_gpu(i,j-l,k,6))
            muy = muy+ccl*(w_aux_gpu(i,j+l,k,7)-w_aux_gpu(i,j-l,k,7))
!
            uz = uz+ccl*(w_aux_gpu(i,j,k+l,2)-w_aux_gpu(i,j,k-l,2))
            vz = vz+ccl*(w_aux_gpu(i,j,k+l,3)-w_aux_gpu(i,j,k-l,3))
            wz = wz+ccl*(w_aux_gpu(i,j,k+l,4)-w_aux_gpu(i,j,k-l,4))
!           tz = tz+ccl*(w_aux_gpu(i,j,k+l,6)-w_aux_gpu(i,j,k-l,6))
            muz = muz+ccl*(w_aux_gpu(i,j,k+l,7)-w_aux_gpu(i,j,k-l,7))
          enddo
          ux = ux*dcsidx_gpu(i)
          vx = vx*dcsidx_gpu(i)
          wx = wx*dcsidx_gpu(i)
!         tx = tx*dcsidx_gpu(i)
          mux = mux*dcsidx_gpu(i)
          uy = uy*detady_gpu(j)
          vy = vy*detady_gpu(j)
          wy = wy*detady_gpu(j)
          ty = ty*detady_gpu(j)
          muy = muy*detady_gpu(j)
          uz = uz*dzitdz_gpu(k)
          vz = vz*dzitdz_gpu(k)
          wz = wz*dzitdz_gpu(k)
!         tz = tz*dzitdz_gpu(k)
          muz = muz*dzitdz_gpu(k)
!         
          if (j==1) then
            wallprop_gpu(i,k,2) = mu*uy
            wallprop_gpu(i,k,3) = mu*wy
            wallprop_gpu(i,k,4) = mu*ty*cploc/Prandtl
          endif
!         
          div     = ux+vy+wz
          div3l   = div/3._rkind
          w_aux_gpu(i,j,k,10) = div3l
          omegax = wy-vz
          omegay = uz-wx
          omegaz = vx-uy
          omod2 = omegax*omegax+omegay*omegay+omegaz*omegaz
          w_aux_gpu(i,j,k,9) = sqrt(omod2)
!
          if (update_sensor == 1) then
!           Original Ducros shock sensor
!           w_aux_gpu(i,j,k,8) = div2/(omod2+div2+1.D-12)
!           
!           Modified Ducros shock sensor
!           w_aux_gpu(i,j,k,8) = div2/(omod2+div2+(u0/l0)**2)
!           
!           Modified Ducros shock sensor (remove expansions)
            w_aux_gpu(i,j,k,8) = (max(-div/sqrt(omod2+div**2+(u0/l0)**2),0._rkind))**2
!           
!           ducfilter = div2/(div2+0.1_rkind*sqrt(uu*uu+vv*vv+ww*ww)*max(dcsidx_gpu(i),detady_gpu(j),dzitdz_gpu(k))+1.D-12)
!           w_aux_gpu(i,j,k,8) = (max(-div/sqrt(omod2+div**2+1.D-12),0._rkind))**2*ducfilter
          endif
!
          sig11 = ux-2._rkind*div3l
          sig12 = vx
          sig13 = wx
          sig21 = uy
          sig22 = vy-2._rkind*div3l
          sig23 = wy
          sig31 = uz
          sig32 = vz
          sig33 = wz-2._rkind*div3l
          sigx  = mux*sig11 + muy*sig12 + muz*sig13
          sigy  = mux*sig12 + muy*sig22 + muz*sig23
          sigz  = mux*sig13 + muy*sig23 + muz*sig33
!         
          sigah = (sig11*ux+sig12*uy+sig13*uz+sig12*vx+sig22*vy+sig23*vz+sig13*wx+sig23*wy+sig33*wz)*mu ! Aerodynamic heating
!         
          sigq  = sigx*uu+sigy*vv+sigz*ww+sigah
!         
          fl_gpu(i,j,k,2) = fl_gpu(i,j,k,2) - sigx
          fl_gpu(i,j,k,3) = fl_gpu(i,j,k,3) - sigy
          fl_gpu(i,j,k,4) = fl_gpu(i,j,k,4) - sigz
          fl_gpu(i,j,k,5) = fl_gpu(i,j,k,5) - sigq
!         
        enddo
      enddo
    enddo
    !@cuf iercuda=cudaDeviceSynchronize()
  endsubroutine visflx_reduced_cuf
!
  subroutine sensor_cuf(nx, ny, nz, ng, u0, l0, w_aux_gpu)
    integer, intent(in) :: nx, ny, nz, ng
    real(rkind), intent(in) :: u0, l0
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(inout), device :: w_aux_gpu
    integer     :: i,j,k,iercuda
    real(rkind) :: div, omod2
!
    !$cuf kernel do(3) <<<*,*>>>
    do k=1,nz
      do j=1,ny
        do i=1,nx
!
          omod2 = w_aux_gpu(i,j,k,9)**2
          div   = 3._rkind*w_aux_gpu(i,j,k,10)
!
!         Original Ducros shock sensor
!         w_aux_gpu(i,j,k,8) = div2/(omod2+div2+1.D-12)
!         
!         Modified Ducros shock sensor
!         w_aux_gpu(i,j,k,8) = div2/(omod2+div2+(u0/l0)**2)
!         
!         Modified Ducros shock sensor (remove expansions)
          w_aux_gpu(i,j,k,8) = (max(-div/sqrt(omod2+div**2+(u0/l0)**2),0._rkind))**2
!         
!         ducfilter = div2/(div2+0.1_rkind*sqrt(uu*uu+vv*vv+ww*ww)*max(dcsidx_gpu(i),detady_gpu(j),dzitdz_gpu(k))+1.D-12)
!         w_aux_gpu(i,j,k,8) = (max(-div/sqrt(omod2+div**2+1.D-12),0._rkind))**2*ducfilter
!
        enddo
      enddo
    enddo
    !@cuf iercuda=cudaDeviceSynchronize()
  endsubroutine sensor_cuf
!
  subroutine visflx_x_cuf(nx, ny, nz, nv, ng, &
    Prandtl, t0, indx_cp_l, indx_cp_r, cp_coeff_gpu, calorically_perfect, &
    x_gpu, w_aux_trans_gpu, fl_trans_gpu, fl_gpu)
!
    integer, intent(in) :: nx, ny, nz, nv, ng, calorically_perfect
    integer, intent(in) :: indx_cp_l, indx_cp_r
    real(rkind), intent(in) :: Prandtl, t0
    real(rkind), dimension(1:ny,1:nx,1:nz,1:nv), intent(inout), device :: fl_trans_gpu
    real(rkind), dimension(1:nx,1:ny,1:nz,1:nv), intent(inout), device :: fl_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in), device :: w_aux_trans_gpu
    real(rkind), dimension(1-ng:), intent(in), device :: x_gpu
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in), device :: cp_coeff_gpu
    integer     :: i,j,k,iv,ll,iercuda
    real(rkind) :: uu,vv,ww,tt,mu,qq
    real(rkind) :: uup,vvp,wwp,ttp,mup,qqp
    real(rkind) :: sigq,sigx,sigy,sigz,sigq_tt,sigq_qq
    real(rkind) :: dxhl,fl2o,fl3o,fl4o,fl5o
    real(rkind) :: ttf,muf
    real(rkind) :: cploc
!   
    !$cuf kernel do(2) <<<*,*>>>
    do k=1,nz
      do j=1,ny
        do i=0,nx
!         
          uu  = w_aux_trans_gpu(j,i  ,k,2)
          uup = w_aux_trans_gpu(j,i+1,k,2)
          vv  = w_aux_trans_gpu(j,i  ,k,3)
          vvp = w_aux_trans_gpu(j,i+1,k,3)
          ww  = w_aux_trans_gpu(j,i  ,k,4)
          wwp = w_aux_trans_gpu(j,i+1,k,4)
          tt  = w_aux_trans_gpu(j,i  ,k,6)
          ttp = w_aux_trans_gpu(j,i+1,k,6)
          mu  = w_aux_trans_gpu(j,i  ,k,7)
          mup = w_aux_trans_gpu(j,i+1,k,7)
          qq  = 0.5_rkind*(uu*uu+vv*vv+ww*ww)
          qqp = 0.5_rkind*(uup*uup+vvp*vvp+wwp*wwp)
!         
          if (calorically_perfect==1) then
            cploc = cp_coeff_gpu(0)
          else
            ttf = 0.5_rkind*(tt+ttp)
            cploc = 0._rkind
            do ll=indx_cp_l,indx_cp_r
              cploc = cploc+cp_coeff_gpu(ll)*(ttf/t0)**ll
            enddo
          endif
!         
          sigx    = uup-uu
          sigy    = vvp-vv
          sigz    = wwp-ww
          sigq_tt = ttp-tt
          sigq_qq = qqp-qq
          muf     = mu+mup
          muf     = 0.5_rkind*muf/(x_gpu(i+1)-x_gpu(i))
!         
          sigx = sigx*muf
          sigy = sigy*muf
          sigz = sigz*muf
          sigq = (sigq_tt*cploc/Prandtl+sigq_qq)*muf
!         
          if (i>0) then
            fl_trans_gpu(j,i,k,2) = fl2o-sigx*dxhl
            fl_trans_gpu(j,i,k,3) = fl3o-sigy*dxhl
            fl_trans_gpu(j,i,k,4) = fl4o-sigz*dxhl
            fl_trans_gpu(j,i,k,5) = fl5o-sigq*dxhl
          endif
          if (i<nx) then
            dxhl = 2._rkind/(x_gpu(i+2)-x_gpu(i))
            fl2o = sigx*dxhl
            fl3o = sigy*dxhl
            fl4o = sigz*dxhl
            fl5o = sigq*dxhl
          endif
!         
        enddo
      enddo
    enddo
    !@cuf iercuda=cudaDeviceSynchronize()
!   
    !$cuf kernel do(3) <<<*,*>>>
    do k=1,nz
      do j=1,ny
        do i=1,nx
          do iv=2,nv
            fl_gpu(i,j,k,iv) = fl_gpu(i,j,k,iv) + fl_trans_gpu(j,i,k,iv)
          enddo
        enddo
      enddo
    enddo
    !@cuf iercuda=cudaDeviceSynchronize()
!
  endsubroutine visflx_x_cuf
!
  subroutine visflx_y_cuf(nx, ny, nz, nv, ng, &
    Prandtl, t0, indx_cp_l, indx_cp_r, cp_coeff_gpu, calorically_perfect, &
    y_gpu, w_aux_gpu, fl_gpu)
!
    integer, intent(in) :: nx, ny, nz, nv, ng, calorically_perfect
    integer, intent(in) :: indx_cp_l, indx_cp_r
    real(rkind), intent(in) :: Prandtl, t0
    real(rkind), dimension(1:nx,1:ny,1:nz,1:nv), intent(inout), device :: fl_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in), device :: w_aux_gpu
    real(rkind), dimension(1-ng:), intent(in), device :: y_gpu
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in), device :: cp_coeff_gpu
    integer     :: i,j,k,iv,ll,iercuda
    real(rkind) :: uu,vv,ww,tt,mu,qq
    real(rkind) :: uup,vvp,wwp,ttp,mup,qqp
    real(rkind) :: sigq,sigx,sigy,sigz,sigq_tt,sigq_qq
    real(rkind) :: dyhl,fl2o,fl3o,fl4o,fl5o
    real(rkind) :: ttf,muf
    real(rkind) :: cploc
!   
    !$cuf kernel do(2) <<<*,*>>>
    do k=1,nz
      do i=1,nx
        do j=0,ny
!         
          uu  = w_aux_gpu(i,j  ,k,2)
          uup = w_aux_gpu(i,j+1,k,2)
          vv  = w_aux_gpu(i,j  ,k,3)
          vvp = w_aux_gpu(i,j+1,k,3)
          ww  = w_aux_gpu(i,j  ,k,4)
          wwp = w_aux_gpu(i,j+1,k,4)
          tt  = w_aux_gpu(i,j  ,k,6)
          ttp = w_aux_gpu(i,j+1,k,6)
          mu  = w_aux_gpu(i,j  ,k,7)
          mup = w_aux_gpu(i,j+1,k,7)
          qq  = 0.5_rkind*(uu*uu+vv*vv+ww*ww)
          qqp = 0.5_rkind*(uup*uup+vvp*vvp+wwp*wwp)
!         
          if (calorically_perfect==1) then
            cploc = cp_coeff_gpu(0)
          else
            ttf = 0.5_rkind*(tt+ttp)
            cploc = 0._rkind
            do ll=indx_cp_l,indx_cp_r
              cploc = cploc+cp_coeff_gpu(ll)*(ttf/t0)**ll
            enddo
          endif
!         
          sigx    = uup-uu
          sigy    = vvp-vv
          sigz    = wwp-ww
          sigq_tt = ttp-tt
          sigq_qq = qqp-qq
          muf     = mu+mup
          muf     = 0.5_rkind*muf/(y_gpu(j+1)-y_gpu(j))
!         
          sigx = sigx*muf
          sigy = sigy*muf
          sigz = sigz*muf
          sigq = (sigq_tt*cploc/Prandtl+sigq_qq)*muf
!         
          if (j>0) then
            fl_gpu(i,j,k,2) = fl_gpu(i,j,k,2) + fl2o-sigx*dyhl
            fl_gpu(i,j,k,3) = fl_gpu(i,j,k,3) + fl3o-sigy*dyhl
            fl_gpu(i,j,k,4) = fl_gpu(i,j,k,4) + fl4o-sigz*dyhl
            fl_gpu(i,j,k,5) = fl_gpu(i,j,k,5) + fl5o-sigq*dyhl
          endif
          if (j<ny) then
            dyhl  = 2._rkind/(y_gpu(j+2)-y_gpu(j))
            fl2o = sigx*dyhl
            fl3o = sigy*dyhl
            fl4o = sigz*dyhl
            fl5o = sigq*dyhl
          endif
!         
        enddo
      enddo
    enddo
    !@cuf iercuda=cudaDeviceSynchronize()
!   
  endsubroutine visflx_y_cuf
! 
  subroutine visflx_z_cuf(nx, ny, nz, nv, ng, &
    Prandtl, t0, indx_cp_l, indx_cp_r, cp_coeff_gpu, calorically_perfect, &
    z_gpu, w_aux_gpu, fl_gpu)
!
    integer, intent(in) :: nx, ny, nz, nv, ng, calorically_perfect
    integer, intent(in) :: indx_cp_l, indx_cp_r
    real(rkind), intent(in) :: Prandtl, t0
    real(rkind), dimension(1:nx,1:ny,1:nz,1:nv), intent(inout), device :: fl_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in), device :: w_aux_gpu
    real(rkind), dimension(1-ng:), intent(in), device :: z_gpu
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in), device :: cp_coeff_gpu
    integer     :: i,j,k,iv,ll,iercuda
    real(rkind) :: uu,vv,ww,tt,mu,qq
    real(rkind) :: uup,vvp,wwp,ttp,mup,qqp
    real(rkind) :: sigq,sigx,sigy,sigz,sigq_tt,sigq_qq
    real(rkind) :: dzhl,fl2o,fl3o,fl4o,fl5o
    real(rkind) :: ttf,muf
    real(rkind) :: cploc
!   
    !$cuf kernel do(2) <<<*,*>>>
    do j=1,ny
      do i=1,nx
        do k=0,nz
!         
          uu  = w_aux_gpu(i,j,k  ,2)
          uup = w_aux_gpu(i,j,k+1,2)
          vv  = w_aux_gpu(i,j,k  ,3)
          vvp = w_aux_gpu(i,j,k+1,3)
          ww  = w_aux_gpu(i,j,k  ,4)
          wwp = w_aux_gpu(i,j,k+1,4)
          tt  = w_aux_gpu(i,j,k  ,6)
          ttp = w_aux_gpu(i,j,k+1,6)
          mu  = w_aux_gpu(i,j,k  ,7)
          mup = w_aux_gpu(i,j,k+1,7)
          qq  = 0.5_rkind*(uu*uu+vv*vv+ww*ww)
          qqp = 0.5_rkind*(uup*uup+vvp*vvp+wwp*wwp)
!         
          if (calorically_perfect==1) then
            cploc = cp_coeff_gpu(0)
          else
            ttf = 0.5_rkind*(tt+ttp)
            cploc = 0._rkind
            do ll=indx_cp_l,indx_cp_r
              cploc = cploc+cp_coeff_gpu(ll)*(ttf/t0)**ll
            enddo
          endif
!         
          sigx    = uup-uu
          sigy    = vvp-vv
          sigz    = wwp-ww
          sigq_tt = ttp-tt
          sigq_qq = qqp-qq
          muf     = mu+mup
          muf     = 0.5_rkind*muf/(z_gpu(k+1)-z_gpu(k))
!         
          sigx = sigx*muf
          sigy = sigy*muf
          sigz = sigz*muf
          sigq = (sigq_tt*cploc/Prandtl+sigq_qq)*muf
!         
          if (k>0) then
            fl_gpu(i,j,k,2) = fl_gpu(i,j,k,2) + fl2o-sigx*dzhl
            fl_gpu(i,j,k,3) = fl_gpu(i,j,k,3) + fl3o-sigy*dzhl
            fl_gpu(i,j,k,4) = fl_gpu(i,j,k,4) + fl4o-sigz*dzhl
            fl_gpu(i,j,k,5) = fl_gpu(i,j,k,5) + fl5o-sigq*dzhl
          endif
          if (k<nz) then
            dzhl  = 2._rkind/(z_gpu(k+2)-z_gpu(k))
            fl2o = sigx*dzhl
            fl3o = sigy*dzhl
            fl4o = sigz*dzhl
            fl5o = sigq*dzhl
          endif
!         
        enddo
      enddo
    enddo
    !@cuf iercuda=cudaDeviceSynchronize()
!   
  endsubroutine visflx_z_cuf
!
  subroutine recyc_exchange_cuf_1(irecyc, w_gpu, wbuf1s_gpu, nx, ny, nz, ng, nv)
    integer, intent(in) :: irecyc, nx, ny, nz, ng, nv
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in), device :: w_gpu
    real(rkind), dimension(:,:,:,:), intent(inout), device :: wbuf1s_gpu
    integer :: i,j,k,m, iercuda
    !$cuf kernel do(2) <<<*,*>>>
    do k=1,nz
      do j=1,ny
        do i=1,ng
          do m=1,nv
            wbuf1s_gpu(i,j,k,m) = w_gpu(irecyc+1-i,j,k,m)
          enddo
        enddo
      enddo
    enddo
    !@cuf iercuda=cudaDeviceSynchronize()
  endsubroutine recyc_exchange_cuf_1
  subroutine recyc_exchange_cuf_2(n1_start_recv, n1_start_send, n1_end_recv, wrecyc_gpu, wbuf1r_gpu, nx, ny, nz, ng, nv)
    integer, intent(in) :: n1_start_recv, n1_start_send, n1_end_recv, nx, ny, nz, ng, nv
    real(rkind), dimension(:,:,:,:), intent(in), device :: wbuf1r_gpu
    real(rkind), dimension(:,:,:,:), intent(inout), device :: wrecyc_gpu
    integer :: i,j,k,m, iercuda
    !$cuf kernel do(2) <<<*,*>>>
    do k=n1_start_recv,n1_end_recv
      do j=1,ny
        do i=1,ng
          do m=1,nv
            wrecyc_gpu(i,j,k,m) = wbuf1r_gpu(i,j,k-n1_start_recv+n1_start_send,m)
          enddo
        enddo
      enddo
    enddo
    !@cuf iercuda=cudaDeviceSynchronize()
  endsubroutine recyc_exchange_cuf_2
  subroutine recyc_exchange_cuf_3(n2_start_recv, n2_start_send, n2_end_recv, wrecyc_gpu, wbuf2r_gpu, nx, ny, nz, ng, nv)
    integer, intent(in) :: n2_start_recv, n2_start_send, n2_end_recv, nx, ny, nz, ng, nv
    real(rkind), dimension(:,:,:,:), intent(in), device :: wbuf2r_gpu
    real(rkind), dimension(:,:,:,:), intent(inout), device :: wrecyc_gpu
    integer :: i,j,k,m, iercuda
    !$cuf kernel do(2) <<<*,*>>>
    do k=n2_start_recv,n2_end_recv
      do j=1,ny
        do i=1,ng
          do m=1,nv
            wrecyc_gpu(i,j,k,m) = wbuf2r_gpu(i,j,k-n2_start_recv+n2_start_send,m)
          enddo
        enddo
      enddo
    enddo
    !@cuf iercuda=cudaDeviceSynchronize()
  endsubroutine recyc_exchange_cuf_3
!
  subroutine bcextr_sub_cuf(ilat, nx, ny, nz, ng, p0, rgas0, w_gpu, indx_cp_l, indx_cp_r, cv_coeff_gpu, t0, calorically_perfect)
!
    integer, intent(in) :: ilat, nx, ny, nz, ng, indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind), intent(in) :: p0, t0, rgas0
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in), device :: cv_coeff_gpu
    real(rkind), dimension(1-ng:, 1-ng:, 1-ng:, 1:), intent(inout), device :: w_gpu
!
    real(rkind) :: rho,rhou,rhov,rhow,tt,ee
    integer :: i,j,k,l,m,iercuda
!
    if (ilat==1) then     ! left side
    elseif (ilat==2) then ! right side
      !$cuf kernel do(2) <<<*,*>>>
      do k=1,nz
        do j=1,ny
          do l=1,ng
            rho  = w_gpu(nx,j,k,1)
            rhou = w_gpu(nx,j,k,2)
            rhov = w_gpu(nx,j,k,3)
            rhow = w_gpu(nx,j,k,4)
            tt   = p0/rho/rgas0
            ee   = get_e_from_temperature_dev(tt, t0, indx_cp_l, indx_cp_r, cv_coeff_gpu, calorically_perfect)
            w_gpu(nx+l,j,k,1) = rho
            w_gpu(nx+l,j,k,2) = rhou
            w_gpu(nx+l,j,k,3) = rhov
            w_gpu(nx+l,j,k,4) = rhow
            w_gpu(nx+l,j,k,5) = rho*ee + 0.5_rkind*(rhou**2+rhov**2+rhow**2)/rho
          enddo
        enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()
    elseif (ilat==3) then ! lower side
    elseif (ilat==4) then  ! upper side
      !$cuf kernel do(2) <<<*,*>>>
      do k=1,nz
        do i=1,nx
          do l=1,ng
            rho  = w_gpu(i,ny,k,1)
            rhou = w_gpu(i,ny,k,2)
            rhov = w_gpu(i,ny,k,3)
            rhow = w_gpu(i,ny,k,4)
            tt   = p0/rho/rgas0
            ee   = get_e_from_temperature_dev(tt, t0, indx_cp_l, indx_cp_r, cv_coeff_gpu, calorically_perfect)
            w_gpu(i,ny+l,k,1) = rho
            w_gpu(i,ny+l,k,2) = rhou
            w_gpu(i,ny+l,k,3) = rhov
            w_gpu(i,ny+l,k,4) = rhow
            w_gpu(i,ny+l,k,5) = rho*ee + 0.5_rkind*(rhou**2+rhov**2+rhow**2)/rho
          enddo
        enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()
!
    elseif (ilat==5) then  ! back side
    elseif (ilat==6) then  ! fore side
    endif
  endsubroutine bcextr_sub_cuf
!
!
  attributes(global) subroutine bc_nr_lat_x_kernel(start_or_end, nr_type, &
    nx, ny, nz, ng, nv, w_aux_gpu, w_gpu, fl_gpu, &
    dcsidx_gpu, indx_cp_l, indx_cp_r, cp_coeff_gpu, winf_gpu, calorically_perfect, rgas0,t0)
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
!
    j = blockDim%x * (blockIdx%x - 1) + threadIdx%x
    k = blockDim%y * (blockIdx%y - 1) + threadIdx%y
    if (j > ny .or. k > nz) return
!
!   Setup min or max boundary
    c_one = [-1.5_rkind, 2._rkind, -0.5_rkind]
    if(start_or_end == 1) then
      i      = 1
      sgn_dw = 1
    elseif(start_or_end == 2) then
      i      = nx
      sgn_dw = -1
    endif
!
!   Compute d(U_cons)/dx: inner, and outer for relaxation
    do m=1,nv
      dw_dn(m) = 0._rkind
      do l=1,3
        dw_dn(m) = dw_dn(m) + sgn_dw * c_one(l)*w_gpu(i+sgn_dw*(l-1),j,k,m)
      enddo
!     Relax to w_gpu which if imposed by recycing works. Problems
!     could arise at the exit, but we do not relax there.
!     Another possible choice would be relaxing to w_inf.
!     w_target       = winf_gpu(m)
      w_target       = w_gpu(i-sgn_dw,j,k,m)
      dw_dn_outer(m) = sgn_dw * (w_gpu(i,j,k,m)-w_target)
    enddo
!
!   Compute eigenvectors
    rho       = w_aux_gpu(i,j,k,1)
    uu        = w_aux_gpu(i,j,k,2)
    vv        = w_aux_gpu(i,j,k,3)
    ww        = w_aux_gpu(i,j,k,4)
    h         = w_aux_gpu(i,j,k,5)
    tt        = w_aux_gpu(i,j,k,6)
    qq        = 0.5_rkind * (uu*uu  +vv*vv + ww*ww)
!   
    gamloc    = get_gamloc_dev(indx_cp_l,indx_cp_r,tt,t0,cp_coeff_gpu,calorically_perfect,rgas0)
    cc        = gamloc * tt * rgas0
    c         =  sqrt(cc)
    ci        =  1._rkind/c
!   
    p_rho     = tt*rgas0
    p_e       = rho*(gamloc-1._rkind)   ! p_e   = rho * R_gas / Cv
    etot      = h - tt*rgas0
!   
    b3        = etot - rho * p_rho/p_e ! qq in the case of calorically perfect
    b2        = p_e/(rho*cc)
    b1        = p_rho/cc - b2*(etot - 2._rkind*qq)
!   
!   b2        =  gm1/cc  ! 1/(cp*T)
!   b1        =  b2 * qq
!   
    call eigenvectors_x(b1, b2, b3, uu, vv, ww, c, ci, h, el, er)
!
!   Pre-multiply to L to get derivative of characteristic variables
    do m=1,5
      dwc_dn(m)       = 0._rkind
      dwc_dn_outer(m) = 0._rkind
      do mm=1,5
        dwc_dn(m)       = dwc_dn(m)       + el(mm,m) * dw_dn(mm)
        dwc_dn_outer(m) = dwc_dn_outer(m) + el(mm,m) * dw_dn_outer(mm)
      enddo
    enddo
!
!   Compute eigenvalues
    ev(1) = uu-c
    ev(2) = uu
    ev(3) = uu+c
    ev(4) = ev(2)
    ev(5) = ev(2)
!
!   If nr_type=1, kill acousting ingoing waves
    if(nr_type == 1) then
      do l=1,5
        ev(l) = sgn_dw*min(sgn_dw*ev(l) ,0._rkind)
      enddo
    endif
!
!   If nr_type=2, exiting waves are kept, incoming waves are assigned from outer derivatives
    if(nr_type == 2) then
      do m=1,5
        if(sgn_dw*ev(m) > 0._rkind) then
          dwc_dn(m) = dwc_dn_outer(m)
        endif
      enddo
    endif
!
!   Compute wave amplitude vector (lambda*L*dw)
    do m=1,5
      dwc_dn(m) = ev(m) * dwc_dn(m)
    enddo
!
!   If nr_type=6, enforce wave reflection to simulate purely reflective wall (LODI relations)
    if(nr_type == 6) then
      dwc_dn(2) = 0._rkind
      if(start_or_end == 1) then
        dwc_dn(3) = dwc_dn(1) ! exiting wave value is used to impose entering wave
      elseif(start_or_end == 2) then
        dwc_dn(1) = dwc_dn(3) ! exiting wave value is used to impose entering wave
      endif
      dwc_dn(4) = 0._rkind
      dwc_dn(5) = 0._rkind
    endif
!
!   Pre-multiply to R to return to conservative variables and assign result to fl
    do m=1,5
      df = 0._rkind
      do mm=1,5
        df = df + er(mm,m) * dwc_dn(mm)
      enddo
      fl_gpu(i,j,k,m) = fl_gpu(i,j,k,m) + df * dcsidx_gpu(i)
    enddo
  endsubroutine bc_nr_lat_x_kernel
!
  attributes(global) subroutine bc_nr_lat_y_kernel(start_or_end, nr_type, &
    nx, ny, nz, ng, nv, w_aux_gpu, w_gpu, fl_gpu, detady_gpu, &
    indx_cp_l, indx_cp_r, cp_coeff_gpu, calorically_perfect,rgas0,t0)
    integer, intent(in), value :: start_or_end, nr_type, nx, ny, nz, ng, nv, indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind), intent(in), value :: rgas0,t0
    real(rkind), dimension(3) :: c_one
    real(rkind), dimension(5) :: dw_dn, dwc_dn, ev
    real(rkind), dimension(5,5) :: el, er
    integer :: i, j, k, l, m, mm, sgn_dw
    real(rkind) :: df, uu, vv, ww, h, qq, cc, c, ci, b2, b1
    real(rkind), dimension(indx_cp_l:indx_cp_r+1) :: cp_coeff_gpu
    real(rkind), dimension(:,:,:,:) :: fl_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: w_aux_gpu, w_gpu
    real(rkind), dimension(:) :: detady_gpu
    real(rkind) :: rho,tt,gamloc,p_rho,p_e,etot,b3
!
    i = blockDim%x * (blockIdx%x - 1) + threadIdx%x
    k = blockDim%y * (blockIdx%y - 1) + threadIdx%y
    if (i > nx .or. k > nz) return
!
!   Setup min or max boundary
    c_one = [-1.5_rkind, 2._rkind, -0.5_rkind]
    if(start_or_end == 1) then
      j      = 1
      sgn_dw = 1
    elseif(start_or_end == 2) then
      j      = ny
      sgn_dw = -1
    endif
!
!   Compute d(U_cons)/dx
    do m=1,nv
      dw_dn(m) = 0._rkind
      do l=1,3
        dw_dn(m) = dw_dn(m) + sgn_dw * c_one(l)*w_gpu(i,j+sgn_dw*(l-1),k,m)
      enddo
    enddo
!
!   Compute eigenvectors
    rho       =  w_aux_gpu(i,j,k,1)
    uu        =  w_aux_gpu(i,j,k,2)
    vv        =  w_aux_gpu(i,j,k,3)
    ww        =  w_aux_gpu(i,j,k,4)
    h         =  w_aux_gpu(i,j,k,5)
    tt        =  w_aux_gpu(i,j,k,6)
    qq        =  0.5_rkind * (uu*uu  +vv*vv + ww*ww)
!   
    gamloc    = get_gamloc_dev(indx_cp_l,indx_cp_r,tt,t0,cp_coeff_gpu,calorically_perfect,rgas0)
    cc        = gamloc * tt * rgas0
    c         =  sqrt(cc)
    ci        =  1._rkind/c
!   
    p_rho     = tt*rgas0
    p_e       = rho*(gamloc-1._rkind)   ! p_e   = rho * R_gas / Cv
    etot      = h - tt*rgas0
!   
    b3        = etot - rho * p_rho/p_e ! qq in the case of calorically perfect
    b2        = p_e/(rho*cc)
    b1        = p_rho/cc - b2*(etot - 2._rkind*qq)
!   
!   b2        =  gm1/cc  ! 1/(cp*T)
!   b1        =  b2 * qq
!   
    call eigenvectors_y(b1, b2, b3, uu, vv, ww, c, ci, h, el, er)
!
!   Pre-multiply to L to get derivative of characteristic variables
    do m=1,5
      dwc_dn(m) = 0._rkind
      do mm=1,5
        dwc_dn(m) = dwc_dn(m) + el(mm,m) * dw_dn(mm)
      enddo
    enddo
!
!   Compute eigenvalues
    ev(1) = vv-c
    ev(2) = vv
    ev(3) = vv+c
    ev(4) = ev(2)
    ev(5) = ev(2)
!
!   If nr_type=1, kill acousting ingoing waves
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
!
!   Compute wave amplitude vector (lambda*L*dw)
    do m=1,5
      dwc_dn(m) = ev(m) * dwc_dn(m)
    enddo
!
!   If nr_type=6, enforce wave reflection to simulate purely reflective wall (LODI relations)
    if(nr_type == 6) then
      dwc_dn(2) = 0._rkind
      if(start_or_end == 1) then
        dwc_dn(3) = dwc_dn(1) ! exiting wave value is used to impose entering wave
      elseif(start_or_end == 2) then
        dwc_dn(1) = dwc_dn(3) ! exiting wave value is used to impose entering wave
      endif
      dwc_dn(4) = 0._rkind
      dwc_dn(5) = 0._rkind
    endif
!
!   Pre-multiply to R to return to conservative variables and assign result to fl
    do m=1,5
      df = 0._rkind
      do mm=1,5
        df = df + er(mm,m) * dwc_dn(mm)
      enddo
      fl_gpu(i,j,k,m) = fl_gpu(i,j,k,m) + df * detady_gpu(j)
    enddo
  endsubroutine bc_nr_lat_y_kernel
!
  attributes(global) subroutine bc_nr_lat_z_kernel(start_or_end, nr_type, &
    nx, ny, nz, ng, nv, w_aux_gpu, w_gpu, fl_gpu, dzitdz_gpu, &
    indx_cp_l, indx_cp_r, cp_coeff_gpu, calorically_perfect,rgas0,t0)
    integer, intent(in), value :: start_or_end, nr_type, nx, ny, nz, ng, nv, indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind), intent(in), value :: rgas0,t0
    real(rkind), dimension(3) :: c_one
    real(rkind), dimension(5) :: dw_dn, dwc_dn, ev
    real(rkind), dimension(5,5) :: el, er
    integer :: i, j, k, l, m, mm, sgn_dw
    real(rkind) :: df, uu, vv, ww, h, qq, cc, c, ci, b2, b1
    real(rkind), dimension(indx_cp_l:indx_cp_r+1) :: cp_coeff_gpu
    real(rkind), dimension(:,:,:,:) :: fl_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: w_aux_gpu, w_gpu
    real(rkind), dimension(:) :: dzitdz_gpu
    real(rkind) :: rho,tt,gamloc,p_rho,p_e,etot,b3
!
    i = blockDim%x * (blockIdx%x - 1) + threadIdx%x
    j = blockDim%y * (blockIdx%y - 1) + threadIdx%y
    if (i > nx .or. j > ny) return
!
!   Setup min or max boundary
    c_one = [-1.5_rkind, 2._rkind, -0.5_rkind]
    if(start_or_end == 1) then
      k      = 1
      sgn_dw = 1
    elseif(start_or_end == 2) then
      k      = nz
      sgn_dw = -1
    endif
!
!   Compute d(U_cons)/dx
    do m=1,nv
      dw_dn(m) = 0._rkind
      do l=1,3
        dw_dn(m) = dw_dn(m) + sgn_dw * c_one(l)*w_gpu(i,j,k+sgn_dw*(l-1),m)
      enddo
    enddo
!
!   Compute eigenvectors
    rho       =  w_aux_gpu(i,j,k,1)
    uu        =  w_aux_gpu(i,j,k,2)
    vv        =  w_aux_gpu(i,j,k,3)
    ww        =  w_aux_gpu(i,j,k,4)
    h         =  w_aux_gpu(i,j,k,5)
    tt        =  w_aux_gpu(i,j,k,6)
    qq        =  0.5_rkind * (uu*uu  +vv*vv + ww*ww)
!   
    gamloc    = get_gamloc_dev(indx_cp_l,indx_cp_r,tt,t0,cp_coeff_gpu,calorically_perfect,rgas0)
    cc        = gamloc * tt * rgas0
    c         =  sqrt(cc)
    ci        =  1._rkind/c
!   
    p_rho     = tt*rgas0
    p_e       = rho*(gamloc-1._rkind)   ! p_e   = rho * R_gas / Cv
    etot      = h - tt*rgas0
!   
    b3        = etot - rho * p_rho/p_e ! qq in the case of calorically perfect
    b2        = p_e/(rho*cc)
    b1        = p_rho/cc - b2*(etot - 2._rkind*qq)
!   
!   b2        =  gm1/cc  ! 1/(cp*T)
!   b1        =  b2 * qq
!   
    call eigenvectors_z(b1, b2, b3, uu, vv, ww, c, ci, h, el, er)
!
!   Pre-multiply to L to get derivative of characteristic variables
    do m=1,5
      dwc_dn(m) = 0._rkind
      do mm=1,5
        dwc_dn(m) = dwc_dn(m) + el(mm,m) * dw_dn(mm)
      enddo
    enddo
!
!   Compute eigenvalues
    ev(1) = ww-c
    ev(2) = ww
    ev(3) = ww+c
    ev(4) = ev(2)
    ev(5) = ev(2)
!
!   If nr_type=1, kill acousting ingoing waves
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
!
!   Compute wave amplitude vector (lambda*L*dw)
    do m=1,5
      dwc_dn(m) = ev(m) * dwc_dn(m)
    enddo
!
!   If nr_type=6, enforce wave reflection to simulate purely reflective wall (LODI relations)
    if(nr_type == 6) then
      dwc_dn(2) = 0._rkind
      if(start_or_end == 1) then
        dwc_dn(3) = dwc_dn(1) ! exiting wave value is used to impose entering wave
      elseif(start_or_end == 2) then
        dwc_dn(1) = dwc_dn(3) ! exiting wave value is used to impose entering wave
      endif
      dwc_dn(4) = 0._rkind
      dwc_dn(5) = 0._rkind
    endif
!
!   Pre-multiply to R to return to conservative variables and assign result to fl
    do m=1,5
      df = 0._rkind
      do mm=1,5
        df = df + er(mm,m) * dwc_dn(mm)
      enddo
      fl_gpu(i,j,k,m) = fl_gpu(i,j,k,m) + df * dzitdz_gpu(k)
    enddo
  endsubroutine bc_nr_lat_z_kernel
!
  subroutine  bcrecyc_cuf_1(nx, ny, nz, ng, nv, wrecycav_gpu, wrecyc_gpu)
    integer, intent(in) :: nx, ny, nz, ng, nv
    real(rkind), dimension(:,:,:), intent(inout), device :: wrecycav_gpu
    real(rkind), dimension(:,:,:,:), intent(in), device :: wrecyc_gpu
    integer :: i,j,k,m,iercuda
    !$cuf kernel do(2) <<<*,*>>>
    do j=1,ny
      do i=1,ng
        do m=1,nv
          wrecycav_gpu(i,j,m) = 0._rkind
          do k=1,nz
            wrecycav_gpu(i,j,m) = wrecycav_gpu(i,j,m)+wrecyc_gpu(i,j,k,m)
          enddo
        enddo
      enddo
    enddo
    !@cuf iercuda=cudaDeviceSynchronize()
  endsubroutine  bcrecyc_cuf_1
!
  subroutine  bcrecyc_cuf_2(nx, ny, nz, nzmax, ng, wrecycav_gpu, wrecyc_gpu)
    integer, intent(in) :: nx, ny, nz, nzmax, ng
    real(rkind), dimension(:,:,:), intent(in), device :: wrecycav_gpu
    real(rkind), dimension(:,:,:,:), intent(inout), device :: wrecyc_gpu
    real(rkind) :: ufav, vfav, wfav, rhom
    integer :: i,j,k,iercuda
!
    !$cuf kernel do(2) <<<*,*>>>
    do j=1,ny
      do i=1,ng
        ufav = wrecycav_gpu(i,j,2)/wrecycav_gpu(i,j,1)
        vfav = wrecycav_gpu(i,j,3)/wrecycav_gpu(i,j,1)
        wfav = wrecycav_gpu(i,j,4)/wrecycav_gpu(i,j,1)
        rhom = wrecycav_gpu(i,j,1)/nzmax
        do k=1,nz
          wrecyc_gpu(i,j,k,2) = wrecyc_gpu(i,j,k,2)/wrecyc_gpu(i,j,k,1)-ufav ! Velocity fluctuations
          wrecyc_gpu(i,j,k,3) = wrecyc_gpu(i,j,k,3)/wrecyc_gpu(i,j,k,1)-vfav
          wrecyc_gpu(i,j,k,4) = wrecyc_gpu(i,j,k,4)/wrecyc_gpu(i,j,k,1)-wfav
          wrecyc_gpu(i,j,k,1) = wrecyc_gpu(i,j,k,1)-rhom                     ! Density fluctuations
        enddo
      enddo
    enddo
    !@cuf iercuda=cudaDeviceSynchronize()
!
  endsubroutine  bcrecyc_cuf_2
!
  subroutine bcrecyc_cuf_3(nx, ny, nz, ng, p0, u0, rgas0, w_gpu, wmean_gpu, wrecyc_gpu, &
    weta_inflow_gpu, map_j_inn_gpu, map_j_out_gpu, &
    yplus_inflow_gpu, eta_inflow_gpu, yplus_recyc_gpu, eta_recyc_gpu, betarecyc, inflow_random_plane_gpu, &
    indx_cp_l, indx_cp_r, cv_coeff_gpu, t0, calorically_perfect)
    integer, intent(in) :: nx, ny, nz, ng, indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind) :: p0, rgas0, betarecyc, t0, u0
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in), device :: cv_coeff_gpu
    real(rkind), dimension(1-ng:,:,:), intent(in), device :: wmean_gpu
    real(rkind), dimension(:,:,:,:), intent(in), device :: wrecyc_gpu
    real(rkind), dimension(:,:,:), intent(in), device :: inflow_random_plane_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,:), intent(inout), device :: w_gpu
    real(rkind), dimension(:), intent(in), device :: weta_inflow_gpu
    real(rkind), dimension(1-ng:), intent(in), device :: yplus_inflow_gpu, eta_inflow_gpu, yplus_recyc_gpu, eta_recyc_gpu
    integer, dimension(:), intent(in), device :: map_j_inn_gpu, map_j_out_gpu
    integer :: i,j,k,iercuda
    real(rkind) :: eta, weta, weta1, bdamp, disty_inn, disty_out, rhofluc, ufluc, vfluc, wfluc, rhof_inn, rhof_out
    real(rkind) :: uf_inn, uf_out, vf_inn, vf_out, wf_inn, wf_out
    real(rkind) :: rhomean, uumean , vvmean , wwmean , tmean  , rho , uu , vv , ww , rhou , rhov , rhow, tt, ee
    integer :: j_inn, j_out
!
    !$cuf kernel do(2) <<<*,*>>>
    do k=1,nz
      do j=1,ny
        eta   = eta_inflow_gpu(j)
        weta  = weta_inflow_gpu(j)
        weta1 = 1._rkind-weta
        bdamp = 0.5_rkind*(1._rkind-tanh(4._rkind*(eta_inflow_gpu(j)-2._rkind)))
        j_inn = map_j_inn_gpu(j)
        j_out = map_j_out_gpu(j)
        disty_inn = (yplus_inflow_gpu(j)-yplus_recyc_gpu(j_inn))/(yplus_recyc_gpu(j_inn+1)-yplus_recyc_gpu(j_inn))
        disty_out = (eta_inflow_gpu(j)-eta_recyc_gpu(j_out))/(eta_recyc_gpu(j_out+1)-eta_recyc_gpu(j_out))
!
        do i=1,ng
!
          if (j==1.or.j_inn>=ny.or.j_out>=ny) then
            rhofluc = 0._rkind
            ufluc   = 0._rkind
            vfluc   = 0._rkind
            wfluc   = 0._rkind
          else
            rhof_inn = wrecyc_gpu(i,j_inn,k,1)*(1._rkind-disty_inn)+wrecyc_gpu(i,j_inn+1,k,1)*disty_inn
            rhof_out = wrecyc_gpu(i,j_out,k,1)*(1._rkind-disty_out)+wrecyc_gpu(i,j_out+1,k,1)*disty_out
            uf_inn   = wrecyc_gpu(i,j_inn,k,2)*(1._rkind-disty_inn)+wrecyc_gpu(i,j_inn+1,k,2)*disty_inn
            uf_out   = wrecyc_gpu(i,j_out,k,2)*(1._rkind-disty_out)+wrecyc_gpu(i,j_out+1,k,2)*disty_out
            vf_inn   = wrecyc_gpu(i,j_inn,k,3)*(1._rkind-disty_inn)+wrecyc_gpu(i,j_inn+1,k,3)*disty_inn
            vf_out   = wrecyc_gpu(i,j_out,k,3)*(1._rkind-disty_out)+wrecyc_gpu(i,j_out+1,k,3)*disty_out
            wf_inn   = wrecyc_gpu(i,j_inn,k,4)*(1._rkind-disty_inn)+wrecyc_gpu(i,j_inn+1,k,4)*disty_inn
            wf_out   = wrecyc_gpu(i,j_out,k,4)*(1._rkind-disty_out)+wrecyc_gpu(i,j_out+1,k,4)*disty_out
!           
            rhofluc = rhof_inn*weta1+rhof_out*weta
            ufluc   =   uf_inn*weta1+  uf_out*weta
            vfluc   =   vf_inn*weta1+  vf_out*weta
            wfluc   =   wf_inn*weta1+  wf_out*weta
            rhofluc = rhofluc*bdamp
            ufluc   = ufluc  *bdamp*betarecyc
            vfluc   = vfluc  *bdamp*betarecyc
            wfluc   = wfluc  *bdamp*betarecyc
            ufluc   = ufluc+0.05_rkind*u0*(inflow_random_plane_gpu(j,k,1)-0.5_rkind)*eta
            vfluc   = vfluc+0.05_rkind*u0*(inflow_random_plane_gpu(j,k,2)-0.5_rkind)*eta
            wfluc   = wfluc+0.05_rkind*u0*(inflow_random_plane_gpu(j,k,3)-0.5_rkind)*eta
          endif
!
          rhomean = wmean_gpu(1-i,j,1)
          uumean  = wmean_gpu(1-i,j,2)/rhomean
          vvmean  = wmean_gpu(1-i,j,3)/rhomean
          wwmean  = wmean_gpu(1-i,j,4)/rhomean
          tmean   = p0/rhomean/rgas0
          rho     = rhomean + rhofluc
          uu      = uumean  + ufluc
          vv      = vvmean  + vfluc
          ww      = wwmean  + wfluc
          rhou    = rho*uu
          rhov    = rho*vv
          rhow    = rho*ww
!
          w_gpu(1-i,j,k,1) = rho
          w_gpu(1-i,j,k,2) = rhou
          w_gpu(1-i,j,k,3) = rhov
          w_gpu(1-i,j,k,4) = rhow
          tt               = p0/rho/rgas0
          ee = get_e_from_temperature_dev(tt, t0, indx_cp_l, indx_cp_r, cv_coeff_gpu, calorically_perfect)
          w_gpu(1-i,j,k,5) = rho*ee + 0.5_rkind*(rhou**2+rhov**2+rhow**2)/rho
        enddo
      enddo
    enddo
    !@cuf iercuda=cudaDeviceSynchronize()
  endsubroutine bcrecyc_cuf_3
!
  subroutine bclam_cuf(ilat, nx, ny, nz, ng, nv, w_gpu, wmean_gpu, p0, rgas0, &
    indx_cp_l, indx_cp_r, cv_coeff_gpu, t0, calorically_perfect)
!
    integer, intent(in) :: nx, ny, nz, ng, nv, indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind), intent(in) :: p0, rgas0, t0
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in), device :: cv_coeff_gpu
    real(rkind), dimension(1-ng:, 1-ng:, 1-ng:, 1:), intent(inout), device :: w_gpu
    real(rkind), dimension(1-ng:,:,:), intent(in), device :: wmean_gpu
!
    integer :: ilat
    integer :: j,k,l,iercuda
    real(rkind) :: rho,rhou,rhov,rhow,tt,ee
!
    if (ilat==1) then     ! left side
      !$cuf kernel do(2) <<<*,*>>>
      do k=1,nz
        do j=1,ny
          do l=1,ng
            rho   = wmean_gpu(1-l,j,1)
            rhou  = wmean_gpu(1-l,j,2)
            rhov  = wmean_gpu(1-l,j,3)
            rhow  = wmean_gpu(1-l,j,4)
            tt    = p0/rho/rgas0
            w_gpu(1-l,j,k,1) = rho
            w_gpu(1-l,j,k,2) = rhou
            w_gpu(1-l,j,k,3) = rhov
            w_gpu(1-l,j,k,4) = rhow
            ee = get_e_from_temperature_dev(tt, t0, indx_cp_l, indx_cp_r, cv_coeff_gpu, calorically_perfect)
            w_gpu(1-l,j,k,5) = rho*ee + 0.5_rkind*(rhou**2+rhov**2+rhow**2)/rho
          enddo
        enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()
    elseif (ilat==2) then ! right side
    elseif (ilat==3) then ! lower side
    elseif (ilat==4) then  ! upper side
    elseif (ilat==5) then  ! back side
    elseif (ilat==6) then  ! fore side
    endif
  endsubroutine bclam_cuf
!
  subroutine bcfree_cuf(ilat, nx, ny, nz, ng, nv, winf_gpu, w_gpu)
    integer, intent(in) :: ilat,nx,ny,nz,ng,nv
    real(rkind), dimension(1:), intent(in), device :: winf_gpu
    real(rkind), dimension(1-ng:, 1-ng:, 1-ng:, 1:), intent(inout), device :: w_gpu
    integer :: j,k,l,m,iercuda
!
    if (ilat==1) then     ! left side
      !$cuf kernel do(2) <<<*,*>>>
      do k=1,nz
        do j=1,ny
          do l=1,ng
            do m=1,nv
              w_gpu(1-l,j,k,m) = winf_gpu(m)
            enddo
          enddo
        enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()
    elseif (ilat==2) then ! right side
    elseif (ilat==3) then ! lower side
    elseif (ilat==4) then  ! upper side
    elseif (ilat==5) then  ! back side
    elseif (ilat==6) then  ! fore side
    endif
  endsubroutine bcfree_cuf
!
  subroutine bcshock_cuf(ilat, nx, ny, nz, ng, nv, w_gpu, winf_gpu, winf_past_shock_gpu, xshock_imp, shock_angle, &
    x_gpu, y_gpu, tanhfacs)
    integer, intent(in) :: nx,ny,nz,ng,nv,ilat
    real(rkind), intent(in) :: xshock_imp, shock_angle, tanhfacs
    real(rkind), dimension(1:), intent(in), device :: winf_gpu, winf_past_shock_gpu
    real(rkind), dimension(1-ng:, 1-ng:, 1-ng:, 1:), intent(inout), device :: w_gpu
    real(rkind), dimension(1-ng:), intent(in), device :: x_gpu, y_gpu
    integer :: i,k,l,m,iercuda
    real(rkind) :: xsh, xx, tanhlen, tanhf, dwinf
!
    tanhlen = 8._rkind*tanhfacs
!
    if (ilat==1) then     ! left side
    elseif (ilat==2) then ! right side
    elseif (ilat==3) then ! lower side
    elseif (ilat==4) then  ! upper side
      !$cuf kernel do(2) <<<*,*>>>
      do k=1,nz
        do i=1,nx
          do l=0,ng
            xsh = xshock_imp-y_gpu(ny+l)/tan(shock_angle)
            xx  = x_gpu(i)-xsh
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
      !@cuf iercuda=cudaDeviceSynchronize()
    elseif (ilat==5) then  ! back side
    elseif (ilat==6) then  ! fore side
    endif
!   
  endsubroutine bcshock_cuf
! 
  subroutine bcextr_var_cuf(nx, ny, nz, ng, w_var_gpu)
    integer :: nx,ny,nz,ng
    integer :: ilat
    integer :: i,j,k,l,m,iercuda
    real(rkind), dimension(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng, 1:1), intent(inout), device :: w_var_gpu
!
    !$cuf kernel do(2) <<<*,*>>>
    do k=1,nz
      do j=1,ny
        do l=1,ng
          w_var_gpu(1-l,j,k,1) = w_var_gpu(1,j,k,1)
        enddo
      enddo
    enddo
    !@cuf iercuda=cudaDeviceSynchronize()
    !$cuf kernel do(2) <<<*,*>>>
    do k=1,nz
      do j=1,ny
        do l=1,ng
          w_var_gpu(nx+l,j,k,1) = w_var_gpu(nx,j,k,1)
        enddo
      enddo
    enddo
    !@cuf iercuda=cudaDeviceSynchronize()
    !$cuf kernel do(2) <<<*,*>>>
    do k=1,nz
      do i=1,nx
        do l=1,ng
          w_var_gpu(i,1-l,k,1) = w_var_gpu(i,1,k,1)
        enddo
      enddo
    enddo
    !@cuf iercuda=cudaDeviceSynchronize()
    !$cuf kernel do(2) <<<*,*>>>
    do k=1,nz
      do i=1,nx
        do l=1,ng
          w_var_gpu(i,ny+l,k,1) = w_var_gpu(i,ny,k,1)
        enddo
      enddo
    enddo
    !@cuf iercuda=cudaDeviceSynchronize()
    !$cuf kernel do(2) <<<*,*>>>
    do j=1,ny
      do i=1,nx
        do l=1,ng
          w_var_gpu(i,j,1-l,1) = w_var_gpu(i,j,1,1)
        enddo
      enddo
    enddo
    !@cuf iercuda=cudaDeviceSynchronize()
    !$cuf kernel do(2) <<<*,*>>>
    do j=1,ny
      do i=1,nx
        do l=1,ng
          w_var_gpu(i,j,nz+l,1) = w_var_gpu(i,j,nz,1)
        enddo
      enddo
    enddo
    !@cuf iercuda=cudaDeviceSynchronize()
  endsubroutine bcextr_var_cuf
!
  subroutine bcextr_cuf(ilat, nx, ny, nz, ng, nv, w_gpu)
    integer :: nx,ny,nz,ng,nv
    integer :: ilat
    integer :: i,j,k,l,m,iercuda
    real(rkind), dimension(1-ng:, 1-ng:, 1-ng:, 1:), intent(inout), device :: w_gpu
!
    if (ilat==1) then     ! left side
    elseif (ilat==2) then ! right side
      !$cuf kernel do(2) <<<*,*>>>
      do k=1,nz
        do j=1,ny
          do l=1,ng
            do m=1,nv
              w_gpu(nx+l,j,k,m) = w_gpu(nx,j,k,m)
            enddo
          enddo
        enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()
    elseif (ilat==3) then ! lower side
      !$cuf kernel do(2) <<<*,*>>>
      do k=1,nz
        do i=1,nx
          do l=1,ng
            do m=1,nv
              w_gpu(i,1-l,k,m) = w_gpu(i,1,k,m)
            enddo
          enddo
        enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()
    elseif (ilat==4) then  ! upper side
      !$cuf kernel do(2) <<<*,*>>>
      do k=1,nz
        do i=1,nx
          do l=1,ng
            do m=1,nv
              w_gpu(i,ny+l,k,m) = w_gpu(i,ny,k,m)
            enddo
          enddo
        enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()
    elseif (ilat==5) then  ! back side
      !$cuf kernel do(2) <<<*,*>>>
      do j=1,ny
        do i=1,nx
          do l=1,ng
            do m=1,nv
              w_gpu(i,j,1-l,m) = w_gpu(i,j,1,m)
            enddo
          enddo
        enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()
    elseif (ilat==6) then  ! fore side
      !$cuf kernel do(2) <<<*,*>>>
      do j=1,ny
        do i=1,nx
          do l=1,ng
            do m=1,nv
              w_gpu(i,j,nz+l,m) = w_gpu(i,j,nz,m)
            enddo
          enddo
        enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()
    endif
  endsubroutine bcextr_cuf
!
  subroutine bcsym_cuf(ilat, nx, ny, nz, ng, twall, w_gpu, w_aux_gpu, indx_cp_l, indx_cp_r, cv_coeff_gpu, calorically_perfect)
    integer :: nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect
    integer :: ilat
    real(rkind) :: twall
    integer :: i,k,l, iercuda
    real(rkind), dimension(1-ng:, 1-ng:, 1-ng:, 1:), intent(inout), device :: w_gpu
    real(rkind), dimension(1-ng:, 1-ng:, 1-ng:, 1:), intent(in), device :: w_aux_gpu
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in), device :: cv_coeff_gpu
!
    if (ilat==1) then     ! left side
    elseif (ilat==2) then ! right side
    elseif (ilat==3) then ! lower side
      !$cuf kernel do(2) <<<*,*>>>
      do k=1,nz
        do i=1,nx
          do l=1,ng
            w_gpu(i,1-l,k,1) = w_gpu(i,1+l,k,1)
            w_gpu(i,1-l,k,2) = w_gpu(i,1+l,k,2)
            w_gpu(i,1-l,k,3) = w_gpu(i,1+l,k,3)
            w_gpu(i,1-l,k,4) = w_gpu(i,1+l,k,4)
            w_gpu(i,1-l,k,5) = w_gpu(i,1+l,k,5)
          enddo
        enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()
    elseif (ilat==4) then  ! upper side
    elseif (ilat==5) then  ! back side
    elseif (ilat==6) then  ! fore side
    endif
  endsubroutine bcsym_cuf
!
  subroutine bcwall_cuf(ilat, nx, ny, nz, ng, twall, w_gpu, w_aux_gpu, indx_cp_l, indx_cp_r, cv_coeff_gpu, t0, rgas0, &
    calorically_perfect, tol_iter_nr)
    integer :: nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect
    integer :: ilat
    real(rkind) :: twall, t0, rgas0, tol_iter_nr
    integer :: i,k,l, iercuda
    real(rkind) :: rho,uu,vv,ww,qq,pp,tt,rhoe,ee
    real(rkind), dimension(1-ng:, 1-ng:, 1-ng:, 1:), intent(inout), device :: w_gpu
    real(rkind), dimension(1-ng:, 1-ng:, 1-ng:, 1:), intent(in), device :: w_aux_gpu
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in), device :: cv_coeff_gpu
!
    if (ilat==1) then     ! left side
    elseif (ilat==2) then ! right side
    elseif (ilat==3) then ! lower side
      !$cuf kernel do(2) <<<*,*>>>
      do k=1,nz
        do i=1,nx
          w_gpu(i,1,k,2) = 0._rkind
          w_gpu(i,1,k,3) = 0._rkind
          w_gpu(i,1,k,4) = 0._rkind
!         w_gpu(i,1,k,5) = w_gpu(i,1,k,1)*gm*twall
          ee = get_e_from_temperature_dev(twall, t0, indx_cp_l, indx_cp_r, cv_coeff_gpu, calorically_perfect)
          w_gpu(i,1,k,5) = w_gpu(i,1,k,1)*ee
          do l=1,ng
            rho  = w_gpu(i,1+l,k,1)
            uu   = w_gpu(i,1+l,k,2)/rho
            vv   = w_gpu(i,1+l,k,3)/rho
            ww   = w_gpu(i,1+l,k,4)/rho
            rhoe = w_gpu(i,1+l,k,5)
            qq   = 0.5_rkind*(uu*uu+vv*vv+ww*ww)
            ee   = rhoe/rho-qq
            tt   = get_temperature_from_e_dev(ee,w_aux_gpu(i,1+l,k,6),t0,cv_coeff_gpu,indx_cp_l,indx_cp_r,calorically_perfect,&
            tol_iter_nr)
            pp   = rho*tt*rgas0
            tt   = 2._rkind*twall-tt ! bc
            ee   = get_e_from_temperature_dev(tt, t0, indx_cp_l, indx_cp_r, cv_coeff_gpu,calorically_perfect)
            rho  = pp/tt/rgas0
            w_gpu(i,1-l,k,1) =  rho
            w_gpu(i,1-l,k,2) = -rho*uu
            w_gpu(i,1-l,k,3) = -rho*vv
            w_gpu(i,1-l,k,4) = -rho*ww
            w_gpu(i,1-l,k,5) =  rho*(ee+qq)
          enddo
        enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()
    elseif (ilat==4) then  ! upper side
    elseif (ilat==5) then  ! back side
    elseif (ilat==6) then  ! fore side
    endif
  endsubroutine bcwall_cuf
!
  subroutine bcwall_staggered_cuf(ilat, nx, ny, nz, ng, twall, w_gpu, w_aux_gpu, indx_cp_l, indx_cp_r, cv_coeff_gpu, t0, rgas0, &
    calorically_perfect, tol_iter_nr)
    integer :: nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect
    integer :: ilat
    real(rkind) :: twall, t0, rgas0, tol_iter_nr
    integer :: i,k,l, iercuda
    real(rkind) :: rho,uu,vv,ww,qq,pp,tt,rhoe,ee
    real(rkind), dimension(1-ng:, 1-ng:, 1-ng:, 1:), intent(inout), device :: w_gpu
    real(rkind), dimension(1-ng:, 1-ng:, 1-ng:, 1:), intent(in), device :: w_aux_gpu
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in), device :: cv_coeff_gpu
!
    if (ilat==1) then     ! left side
    elseif (ilat==2) then ! right side
    elseif (ilat==3) then ! lower side
      !$cuf kernel do(2) <<<*,*>>>
      do k=1,nz
        do i=1,nx
          do l=1,ng
            rho  = w_gpu(i,l,k,1)
            uu   = w_gpu(i,l,k,2)/rho
            vv   = w_gpu(i,l,k,3)/rho
            ww   = w_gpu(i,l,k,4)/rho
            rhoe = w_gpu(i,l,k,5)
            qq   = 0.5_rkind*(uu*uu+vv*vv+ww*ww)
            ee   = rhoe/rho-qq
            tt   = get_temperature_from_e_dev(ee,w_aux_gpu(i,l,k,6),t0,cv_coeff_gpu,indx_cp_l,indx_cp_r,calorically_perfect,&
            tol_iter_nr)
            pp   = rho*tt*rgas0
            tt   = 2._rkind*twall-tt ! bc
            ee   = get_e_from_temperature_dev(tt, t0, indx_cp_l, indx_cp_r, cv_coeff_gpu,calorically_perfect)
            rho  = pp/tt/rgas0
            w_gpu(i,1-l,k,1) =  rho
            w_gpu(i,1-l,k,2) = -rho*uu
            w_gpu(i,1-l,k,3) = -rho*vv
            w_gpu(i,1-l,k,4) = -rho*ww
            w_gpu(i,1-l,k,5) =  rho*(ee+qq)
          enddo
        enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()
    elseif (ilat==4) then  ! upper side
      !$cuf kernel do(2) <<<*,*>>>
      do k=1,nz
        do i=1,nx
          do l=1,ng
            rho  = w_gpu(i,ny+1-l,k,1)
            uu   = w_gpu(i,ny+1-l,k,2)/rho
            vv   = w_gpu(i,ny+1-l,k,3)/rho
            ww   = w_gpu(i,ny+1-l,k,4)/rho
            rhoe = w_gpu(i,ny+1-l,k,5)
            qq   = 0.5_rkind*(uu*uu+vv*vv+ww*ww)
            ee   = rhoe/rho-qq
            tt   = get_temperature_from_e_dev(ee,w_aux_gpu(i,ny+1-l,k,6),t0,cv_coeff_gpu,indx_cp_l,indx_cp_r,calorically_perfect,&
            tol_iter_nr)
            pp   = rho*tt*rgas0
            tt   = 2._rkind*twall-tt ! bc
            ee   = get_e_from_temperature_dev(tt, t0, indx_cp_l, indx_cp_r, cv_coeff_gpu,calorically_perfect)
            rho  = pp/tt/rgas0
            w_gpu(i,ny+l,k,1) =  rho
            w_gpu(i,ny+l,k,2) = -rho*uu
            w_gpu(i,ny+l,k,3) = -rho*vv
            w_gpu(i,ny+l,k,4) = -rho*ww
            w_gpu(i,ny+l,k,5) =  rho*(ee+qq)
          enddo
        enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()
    elseif (ilat==5) then  ! back side
    elseif (ilat==6) then  ! fore side
    endif
  endsubroutine bcwall_staggered_cuf
!
  subroutine compute_residual_cuf(nx, ny, nz, ng, nv, fln_gpu, dt, residual_rhou, fluid_mask_gpu)
    integer :: nx, ny, nz, ng, nv
    real(rkind), intent(out) :: residual_rhou
    real(rkind), intent(in) :: dt
    real(rkind), dimension(1:nx, 1:ny, 1:nz, nv), intent(in), device :: fln_gpu
    integer, dimension(1-ng:,1-ng:,1-ng:), intent(in), device :: fluid_mask_gpu
    integer :: i,j,k,iercuda
!   Note: should be modified to include metrics
!
    residual_rhou = 0._rkind
    !$cuf kernel do(2) <<<*,*>>> reduce(+:residual_rhou)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          if (fluid_mask_gpu(i,j,k)==0) then
            residual_rhou = residual_rhou + (fln_gpu(i,j,k,2)/dt)**2
          endif
        enddo
      enddo
    enddo
    !@cuf iercuda=cudaDeviceSynchronize()
  endsubroutine compute_residual_cuf
!
  subroutine compute_dt_cuf(nx, ny, nz, ng, rgas0, Prandtl, &
    dcsidx_gpu, detady_gpu, dzitdz_gpu, dcsidxs_gpu, detadys_gpu, dzitdzs_gpu, w_gpu, w_aux_gpu, &
    dtxi_max, dtyi_max, dtzi_max, dtxv_max, dtyv_max, dtzv_max, dtxk_max, dtyk_max, dtzk_max,    &
    indx_cp_l, indx_cp_r, cp_coeff_gpu,fluid_mask_gpu,calorically_perfect,t0)
    integer :: nx, ny, nz, ng, indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind) :: rgas0, t0
    real(rkind) :: Prandtl
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:) , intent(in), device :: w_gpu
    real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(in), device :: w_aux_gpu
    integer, dimension(1-ng:,1-ng:,1-ng:), intent(in), device :: fluid_mask_gpu
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in), device :: cp_coeff_gpu
    real(rkind), dimension(1:), intent(in), device :: dcsidx_gpu, detady_gpu, dzitdz_gpu
    real(rkind), dimension(1:), intent(in), device :: dcsidxs_gpu, detadys_gpu, dzitdzs_gpu
    real(rkind) :: dtxi, dtyi, dtzi, dtxv, dtyv, dtzv, dtxk, dtyk, dtzk
    integer     :: i,j,k,ll,iercuda
    real(rkind) :: rho, ri, uu, vv, ww, tt, mu, nu, k_over_rhocp, c
    real(rkind) :: dtxi_max, dtyi_max, dtzi_max, dtxv_max, dtyv_max, dtzv_max, dtxk_max, dtyk_max, dtzk_max
    real(rkind) :: gamloc
!
    dtxi_max = 0._rkind
    dtyi_max = 0._rkind
    dtzi_max = 0._rkind
    dtxv_max = 0._rkind
    dtyv_max = 0._rkind
    dtzv_max = 0._rkind
    dtxk_max = 0._rkind
    dtyk_max = 0._rkind
    dtzk_max = 0._rkind
!
    !$cuf kernel do(2) <<<*,*>>> reduce(max:dtxi_max,dtyi_max,dtzi_max,dtxv_max,dtyv_max,dtzv_max,dtxk_max,dtyk_max,dtzk_max)
    do k=1,nz
      do j=1,ny
        do i=1,nx
!         
          if (fluid_mask_gpu(i,j,k)==0) then
!           
            rho  = w_gpu(i,j,k,1)
            ri   = 1._rkind/rho
            uu   = w_aux_gpu(i,j,k,2)
            vv   = w_aux_gpu(i,j,k,3)
            ww   = w_aux_gpu(i,j,k,4)
            tt   = w_aux_gpu(i,j,k,6)
            mu   = w_aux_gpu(i,j,k,7)
!           
            nu  = ri*mu
            k_over_rhocp = nu/Prandtl
            gamloc = get_gamloc_dev(indx_cp_l,indx_cp_r,tt,t0,cp_coeff_gpu,calorically_perfect,rgas0)
            c      = sqrt (gamloc*rgas0*tt)
            dtxi   = (abs(uu)+c)*dcsidx_gpu(i)
            dtyi   = (abs(vv)+c)*detady_gpu(j)
            dtzi   = (abs(ww)+c)*dzitdz_gpu(k)
            dtxv  = nu*dcsidxs_gpu(i)
            dtyv  = nu*detadys_gpu(j)
            dtzv  = nu*dzitdzs_gpu(k)
            dtxk  = k_over_rhocp*dcsidxs_gpu(i)
            dtyk  = k_over_rhocp*detadys_gpu(j)
            dtzk  = k_over_rhocp*dzitdzs_gpu(k)
!
            dtxi_max = max(dtxi_max, dtxi)
            dtyi_max = max(dtyi_max, dtyi)
            dtzi_max = max(dtzi_max, dtzi)
            dtxv_max = max(dtxv_max, dtxv)
            dtyv_max = max(dtyv_max, dtyv)
            dtzv_max = max(dtzv_max, dtzv)
            dtxk_max = max(dtxk_max, dtxk)
            dtyk_max = max(dtyk_max, dtyk)
            dtzk_max = max(dtzk_max, dtzk)
!           
          endif
!         
        enddo
      enddo
    enddo
    !@cuf iercuda=cudaDeviceSynchronize()
  endsubroutine compute_dt_cuf
!
  subroutine eval_aux_cuf(nx, ny, nz, ng, istart, iend, jstart, jend, kstart, kend, w_gpu, w_aux_gpu, &
    visc_model, mu0, t0, sutherland_S, T_ref_dim, &
    powerlaw_vtexp, VISC_POWER, VISC_SUTHERLAND, VISC_NO, &
    cv_coeff_gpu, indx_cp_l, indx_cp_r, rgas0, calorically_perfect, tol_iter_nr, stream_id)
    integer(ikind), intent(in) :: nx, ny, nz, ng, visc_model
    integer(ikind), intent(in) :: istart, iend, jstart, jend, kstart, kend
    integer(ikind), intent(in) :: VISC_POWER, VISC_SUTHERLAND, VISC_NO, indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind),    intent(in) :: mu0, t0, sutherland_S, T_ref_dim, powerlaw_vtexp, rgas0, tol_iter_nr
    real(rkind),    dimension(indx_cp_l:indx_cp_r+1), intent(in),  device   :: cv_coeff_gpu
    real(rkind),    dimension(1-ng:,1-ng:,1-ng:, 1:), intent(in),  device   :: w_gpu
    real(rkind),    dimension(1-ng:,1-ng:,1-ng:, 1:), intent(inout), device :: w_aux_gpu
    integer(kind=cuda_stream_kind), intent(in) :: stream_id
    integer(ikind)                        :: i, j, k
    real(rkind)                           :: rho, rhou, rhov, rhow, rhoe, ri, uu, vv, ww, qq, pp, tt, mu, ee
    integer(ikind)                        :: iercuda
!
!   real(rkind) :: T_start, T_old, ee, ebar, den, num, T_pow, T_powp, tt2
!   integer :: l
!
    !$cuf kernel do(3) <<<*,*,stream=stream_id>>>
    do k=kstart,kend
      do j=jstart,jend
        do i=istart,iend
!         w_aux(:) : rho, u, v, w, h, T, viscosity, div, |omega|, ducros
          rho  = w_gpu(i,j,k,1)
          rhou = w_gpu(i,j,k,2)
          rhov = w_gpu(i,j,k,3)
          rhow = w_gpu(i,j,k,4)
          rhoe = w_gpu(i,j,k,5)
          ri   = 1._rkind/rho
          uu   = rhou*ri
          vv   = rhov*ri
          ww   = rhow*ri
          qq   = 0.5_rkind*(uu*uu+vv*vv+ww*ww)
          ee = rhoe/rho-qq
          tt = get_temperature_from_e_dev(ee, w_aux_gpu(i,j,k,6), t0, cv_coeff_gpu, indx_cp_l, indx_cp_r, &
          calorically_perfect, tol_iter_nr)
          pp = rho*tt*rgas0
!
          w_aux_gpu(i,j,k,1) = rho
          w_aux_gpu(i,j,k,2) = uu
          w_aux_gpu(i,j,k,3) = vv
          w_aux_gpu(i,j,k,4) = ww
          w_aux_gpu(i,j,k,5) = (rhoe+pp)/rho
          w_aux_gpu(i,j,k,6) = tt
!
          if (visc_model == VISC_POWER) then
            mu = mu0 * (tt/t0)**powerlaw_vtexp
          elseif (visc_model == VISC_SUTHERLAND) then
            mu = mu0 * (tt/t0)**1.5_rkind * &
            (1._rkind+sutherland_S/T_ref_dim)/(tt/t0 + sutherland_S/T_ref_dim)
          elseif (visc_model == VISC_NO) then
            mu = 0._rkind
          endif
          w_aux_gpu(i,j,k,7) = mu
!         STREAMS v1.0 : mu0 = sqgmr  ;  ggmopr = cp/Pr ; k = sqgmr * ggmopr
        enddo
      enddo
    enddo
  endsubroutine eval_aux_cuf
!
  attributes(device) function get_gamloc_dev(indx_cp_l,indx_cp_r,tt,t0,cp_coeff_gpu,calorically_perfect,rgas0)
    real(rkind) :: get_gamloc_dev
    integer, value :: indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind), value :: tt, rgas0, t0
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), device :: cp_coeff_gpu
    real(rkind) :: cploc, gamloc
    integer :: l
!
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
!
  endfunction get_gamloc_dev
!
  attributes(device) function get_temperature_from_e_dev(ee, T_start, t0, cv_coeff_gpu, indx_cp_l, indx_cp_r, &
    calorically_perfect,tol_iter_nr)
    real(rkind) :: get_temperature_from_e_dev
    integer, value :: indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind), value :: ee, T_start, t0, tol_iter_nr
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), device :: cv_coeff_gpu
    real(rkind) :: tt, T_old, ebar, den, num, T_pow, T_powp
    integer :: l,iter,max_iter
!
    max_iter = 50
!
    if (calorically_perfect==1) then
!     tt  = t0+ee/cv_coeff_gpu(0)
      tt  =    ee/cv_coeff_gpu(0)
    else
      T_old = T_start
      ebar  = ee - cv_coeff_gpu(indx_cp_r+1)*t0
      do iter=1,max_iter
        den = 0._rkind
        num = 0._rkind
        do l=indx_cp_l,indx_cp_r
          if (l==-1) then
            T_pow  = (T_old/t0)**l
            den    = den+cv_coeff_gpu(l)*T_pow
            num    = num+cv_coeff_gpu(l)*log(T_old/t0)
          else
            T_pow  = (T_old/t0)**l
            T_powp = (T_old/t0)*T_pow
            den    = den+cv_coeff_gpu(l)*T_pow
            num    = num+cv_coeff_gpu(l)*(T_powp-1._rkind)/(l+1._rkind)
          endif
        enddo
        num = num*t0
        tt = T_old+(ebar-num)/den
        if (abs(tt-T_old) < tol_iter_nr) exit
        T_old = tt
      enddo
    endif
    get_temperature_from_e_dev = tt
  endfunction get_temperature_from_e_dev
!
  attributes(device) function get_e_from_temperature_dev(tt, t0, indx_cp_l, indx_cp_r, cv_coeff_gpu, calorically_perfect)
    real(rkind) :: get_e_from_temperature_dev
    real(rkind),value :: tt, t0
    integer,value     :: indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind), dimension(indx_cp_l:indx_cp_r+1) :: cv_coeff_gpu
    real(rkind) :: ee
    integer :: l
!
    if (calorically_perfect==1) then
!     ee = cv_coeff_gpu(0)*(tt-t0)
      ee = cv_coeff_gpu(0)* tt
    else
      ee = cv_coeff_gpu(indx_cp_r+1)
      do l=indx_cp_l,indx_cp_r
        if (l==-1) then
          ee = ee+cv_coeff_gpu(l)*log(tt/t0)
        else
          ee = ee+cv_coeff_gpu(l)/(l+1._rkind)*((tt/t0)**(l+1)-1._rkind)
        endif
      enddo
      ee = ee*t0
    endif
    get_e_from_temperature_dev = ee
  endfunction get_e_from_temperature_dev
!
! 
!
!
!
! eigs33
!
!
!
!
! NOMANAGEDsubroutine copy_to_psi_pv_managed_cuf(nxsl_ins,nxel_ins,nysl_ins,nyel_ins,nzsl_ins,nzel_ins, &
! NOMANAGED         ng,nx,ny,nz,npsi, npsi_pv, nv_aux, n_aux_list, &
! NOMANAGED         psi_gpu,psi_pv_managed,w_aux_gpu,aux_list_gpu)
! NOMANAGED integer :: i,j,k,l,ll, iercuda
! NOMANAGED real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv_aux), device :: w_aux_gpu
!
! NOMANAGED !$cuf kernel do(3) <<<*,*>>>
! NOMANAGED do k=nzsl_ins,nzel_ins
! NOMANAGED  do j=nysl_ins,nyel_ins
! NOMANAGED   do i=nxsl_ins,nxel_ins
! NOMANAGED    do l=1,n_aux_list
! NOMANAGED     ll = aux_list_gpu(l)
! NOMANAGED     psi_pv_managed(i,j,k,l)  = w_aux_gpu(i,j,k,ll)
! NOMANAGED    enddo
! NOMANAGED    do l=1,npsi
! NOMANAGED     psi_pv_managed(i,j,k,n_aux_list+l) = psi_gpu(i,j,k,l)
! NOMANAGED    enddo
! NOMANAGED   enddo
! NOMANAGED  enddo
! NOMANAGED enddo
! NOMANAGED !@cuf iercuda=cudaDeviceSynchronize()
! NOMANAGEDend subroutine copy_to_psi_pv_managed_cuf
!
  subroutine probe_interpolation_cuf(num_probe,nx,ny,nz,ng,nv_aux,ijk_probe_gpu,w_aux_probe_gpu,w_aux_gpu,probe_coeff_gpu)
    integer, intent(in) :: num_probe,nx,ny,nz,ng,nv_aux
    real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv_aux), intent(in), device :: w_aux_gpu
    real(rkind), dimension(6,num_probe), intent(inout), device :: w_aux_probe_gpu
    real(rkind), dimension(2,2,2,num_probe), intent(in), device :: probe_coeff_gpu
    integer, dimension(3,num_probe), intent(in), device :: ijk_probe_gpu
    integer :: i,j,k,ii,jj,kk,l,iercuda
    real(rkind) :: w1,w2,w3,w4,w5,w6
!
    !$cuf kernel do(1) <<<*,*>>>
    do l=1,num_probe
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
    !@cuf iercuda=cudaDeviceSynchronize()
!   
  endsubroutine probe_interpolation_cuf
!
endmodule streams_kernels_gpu
!
