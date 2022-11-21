module streams_kernels_cpu

    use streams_parameters, only : rkind, ikind, REAL64
    implicit none

contains

    subroutine init_flux_cpu(nx, ny, nz, nv, fl_cpu, fln_cpu, rhodt) 
        integer :: nx, ny, nz, nv
        real(rkind) :: rhodt
        real(rkind), dimension(1:,1:,1:,1:), intent(inout) :: fl_cpu, fln_cpu
        integer :: i,j,k,m
         do k=1,nz
          do j=1,ny
           do i=1,nx
            do m=1,nv
             fln_cpu(i,j,k,m) = - rhodt * fl_cpu(i,j,k,m)
             fl_cpu(i,j,k,m)  = 0._rkind
            enddo
           enddo
          enddo
         enddo
    endsubroutine init_flux_cpu

    subroutine euler_x_transpose_cpu(nx, ny, nz, ng, nv_aux, w_aux_cpu, w_aux_trans_cpu)
        integer :: nx,ny,nz,ng,nv_aux
        real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv_aux), intent(in) :: w_aux_cpu
        real(rkind), dimension(1-ng:ny+ng,1-ng:nx+ng,1-ng:nz+ng,1:8), intent(inout) :: w_aux_trans_cpu
        integer :: i,j,k,iv
        do k=1,nz
         do i=1-ng,nx+ng
          do j=1,ny
           do iv=1,8
            w_aux_trans_cpu(j,i,k,iv) = w_aux_cpu(i,j,k,iv)
           enddo
          enddo
         enddo
        enddo
    endsubroutine euler_x_transpose_cpu

     subroutine euler_x_fluxes_hybrid_kernel(nv, nv_aux, nx, ny, nz, ng, &
        eul_imin, eul_imax, lmax_base, nkeep, cp0, cv0, coeff_deriv1_cpu, dcsidx_cpu, w_aux_trans_cpu, fhat_trans_cpu, &
        force_zero_flux_min, force_zero_flux_max, &
        weno_scheme, sensor_threshold, weno_size, gplus_x_cpu, gminus_x_cpu, cp_coeff_cpu, &
        indx_cp_l, indx_cp_r, order_modify_x_cpu, calorically_perfect, tol_iter_nr)

        implicit none
        ! Passed arguments
        integer :: nv, nx, ny, nz, ng, nv_aux
        integer :: eul_imin, eul_imax, lmax_base, nkeep, indx_cp_l, indx_cp_r, calorically_perfect
        real(rkind), dimension(1-ng:ny+ng,1-ng:nx+ng,1-ng:nz+ng,1:8) :: w_aux_trans_cpu
        real(rkind), dimension(1-ng:ny+ng,1-ng:nx+ng,1-ng:nz+ng,1:nv) :: fhat_trans_cpu
        integer, dimension(0:ny,0:nx,0:nz) :: order_modify_x_cpu
        real(rkind), dimension(indx_cp_l:indx_cp_r) :: cp_coeff_cpu
        real(rkind), dimension(4,4) :: coeff_deriv1_cpu
        real(rkind), dimension(nx) :: dcsidx_cpu
        integer :: force_zero_flux_min, force_zero_flux_max
        real(rkind) :: sensor_threshold, cp0, cv0, tol_iter_nr
        integer :: weno_scheme, weno_size
        ! Local variables
        integer :: i, j, k, m, l
        real(rkind) :: fh1, fh2, fh3, fh4, fh5
        real(rkind) :: rhom, uui, vvi, wwi, ppi, enti, rhoi, tti
        real(rkind) :: uuip, vvip, wwip, ppip, entip, rhoip, ttip
        real(rkind) :: ft1, ft2, ft3, ft4, ft5, ft6
        real(rkind) :: uvs1, uvs2, uvs3, uvs4, uvs6, uv_part
!       real(rkind) :: uvs5
        integer :: ii, lmax, wenorec_ord
        integer :: ishk
        real(rkind) :: b1, b2, b3, c, ci, h, uu, vv, ww
        real(rkind), dimension(5,5) :: el, er
        real(rkind), dimension(5) :: ev, evmax, fi, ghat, gl, gr
        real(rkind), dimension(nv,2*weno_scheme,ny,nz) :: gplus_x_cpu
        real(rkind), dimension(nv,2*weno_scheme,ny,nz) :: gminus_x_cpu
        integer :: ll, mm
        real(rkind) :: rho, pp, wc, gc, rhou
        real(rkind) :: tt, gamloc
        real(rkind) :: uvs5_i,uvs5_k,uvs5_p,eei,eeip
        real(rkind) :: drho, dee, eem
        real(rkind) :: drhof, deef
        real(rkind) :: sumnumrho,sumnumee,sumdenrho,sumdenee
        integer :: n,n2



        do j=1,ny
do k=1,nz


        do i=eul_imin-1,eul_imax

            ishk = 0
            do ii=i-weno_scheme+1,i+weno_scheme
                if (w_aux_trans_cpu(j,ii,k,8) > sensor_threshold) ishk = 1
            enddo

            if (ishk == 0) then

                ft1 = 0._rkind
                ft2 = 0._rkind
                ft3 = 0._rkind
                ft4 = 0._rkind
                ft5 = 0._rkind
                ft6 = 0._rkind
                lmax = lmax_base+order_modify_x_cpu(j,i,k)
                do l=1,lmax
                    uvs1 = 0._rkind
                    uvs2 = 0._rkind
                    uvs3 = 0._rkind
                    uvs4 = 0._rkind
!                   uvs5 = 0._rkind
                    uvs5_i = 0._rkind
                    uvs5_k = 0._rkind
                    uvs5_p = 0._rkind
                    uvs6 = 0._rkind
                    do m=0,l-1

                        rhoi  = w_aux_trans_cpu(j,i-m,k,1)
                        uui   = w_aux_trans_cpu(j,i-m,k,2)
                        vvi   = w_aux_trans_cpu(j,i-m,k,3)
                        wwi   = w_aux_trans_cpu(j,i-m,k,4)
                        enti  = w_aux_trans_cpu(j,i-m,k,5)
                        tti   = w_aux_trans_cpu(j,i-m,k,6)
                        ppi   = tti*rhoi
                        eei   = enti-tti-0.5_rkind*(uui*uui+vvi*vvi+wwi*wwi)

                        rhoip = w_aux_trans_cpu(j,i-m+l,k,1)
                        uuip  = w_aux_trans_cpu(j,i-m+l,k,2)
                        vvip  = w_aux_trans_cpu(j,i-m+l,k,3)
                        wwip  = w_aux_trans_cpu(j,i-m+l,k,4)
                        entip = w_aux_trans_cpu(j,i-m+l,k,5)
                        ttip  = w_aux_trans_cpu(j,i-m+l,k,6)
                        ppip  = ttip*rhoip
                        eeip  = entip-ttip-0.5_rkind*(uuip*uuip+vvip*vvip+wwip*wwip)

                        rhom  = rhoi+rhoip
                        eem   = eei + eeip

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

                        uv_part = (uui+uuip) * rhom * drhof
                        uvs1 = uvs1 + uv_part * (2._rkind)
                        uvs2 = uvs2 + uv_part * (uui+uuip)
                        uvs3 = uvs3 + uv_part * (vvi+vvip)
                        uvs4 = uvs4 + uv_part * (wwi+wwip)
!                       uvs5 = uvs5 + uv_part * (enti+entip)
                        uvs5_i = uvs5_i + uv_part * eem * deef
                        uvs5_k = uvs5_k + uv_part * (uui*uuip+vvi*vvip+wwi*wwip)
                        uvs5_p = uvs5_p + 4._rkind*(uui*ppip+uuip*ppi)
                        uvs6 = uvs6 + (2._rkind)*(ppi+ppip)
                    enddo
                    ft1 = ft1 + coeff_deriv1_cpu(l,lmax)*uvs1
                    ft2 = ft2 + coeff_deriv1_cpu(l,lmax)*uvs2
                    ft3 = ft3 + coeff_deriv1_cpu(l,lmax)*uvs3
                    ft4 = ft4 + coeff_deriv1_cpu(l,lmax)*uvs4
                    !ft5 = ft5 + coeff_deriv1_cpu(l,lmax)*uvs5
                    ft5 = ft5 + coeff_deriv1_cpu(l,lmax)*(uvs5_i+uvs5_k+uvs5_p)
                    ft6 = ft6 + coeff_deriv1_cpu(l,lmax)*uvs6
                enddo
                fh1 = 0.25_rkind*ft1
                fh2 = 0.25_rkind*ft2 
                fh3 = 0.25_rkind*ft3
                fh4 = 0.25_rkind*ft4
                fh5 = 0.25_rkind*ft5
                if((i==0 .and. force_zero_flux_min == 1).or.(i==nx .and. force_zero_flux_max == 1)) then
                   fh1 = 0._rkind
                   fh2 = 0._rkind
                   fh3 = 0._rkind
                   fh4 = 0._rkind
                   fh5 = 0._rkind
                endif
                fh2 = fh2 + 0.5_rkind*ft6

                fhat_trans_cpu(j,i,k,1) = fh1
                fhat_trans_cpu(j,i,k,2) = fh2
                fhat_trans_cpu(j,i,k,3) = fh3
                fhat_trans_cpu(j,i,k,4) = fh4
                fhat_trans_cpu(j,i,k,5) = fh5
            else
                call compute_roe_average(nx, ny, nz, ng, j, j, i, i+1, k, k, w_aux_trans_cpu, cp0, cv0, &
                                         b1, b2, b3, c, ci, h, uu, vv, ww, cp_coeff_cpu, indx_cp_l, indx_cp_r, &
                                         calorically_perfect, tol_iter_nr)

                call eigenvectors_x(b1, b2, b3, uu, vv, ww, c, ci, h, el, er)

                do m=1,5 ! loop on characteristic fields
                    evmax(m) = -1._rkind
                enddo
                do l=1,weno_size ! LLF
                    ll   = i + l - weno_scheme
                    uu   = w_aux_trans_cpu(j,ll,k,2)
                    tt   = w_aux_trans_cpu(j,ll,k,6)
                    gamloc = get_gamloc_dev(indx_cp_l,indx_cp_r,tt,cp_coeff_cpu,calorically_perfect)
                    c    = sqrt (gamloc*tt)
                    ev(1) = abs(uu-c)
                    ev(2) = abs(uu)
                    ev(3) = abs(uu+c)
                    ev(4) = ev(2)
                    ev(5) = ev(2)
                    do m=1,5
                        evmax(m) = max(ev(m),evmax(m))
                    enddo
                enddo
                do l=1,weno_size ! loop over the stencil centered at face i
                    ll = i + l - weno_scheme

                    rho    = w_aux_trans_cpu(j,ll,k,1)
                    uu     = w_aux_trans_cpu(j,ll,k,2)
                    vv     = w_aux_trans_cpu(j,ll,k,3)
                    ww     = w_aux_trans_cpu(j,ll,k,4)
                    h      = w_aux_trans_cpu(j,ll,k,5) 
                    rhou   = rho*uu
                    pp     = rho*w_aux_trans_cpu(j,ll,k,6)
                    fi(1)  =      rhou
                    fi(2)  = uu * rhou + pp
                    fi(3)  = vv * rhou
                    fi(4)  = ww * rhou
                    fi(5)  = h  * rhou
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

                        gplus_x_cpu (m,l,j,k) = 0.5_rkind * (gc + evmax(m) * wc)
                        gminus_x_cpu(m,l,j,k) = gc - gplus_x_cpu(m,l,j,k)
                    enddo
                enddo
!
!               Reconstruction of the '+' and '-' fluxes
!
                wenorec_ord = weno_scheme+order_modify_x_cpu(j,i,k)
                call wenorec(nv,gplus_x_cpu(:,:,j,k),gminus_x_cpu(:,:,j,k),gl,gr,weno_scheme,wenorec_ord)
!
                do m=1,5
                    ghat(m) = gl(m) + gr(m) ! char. flux
                enddo

!               !Return to conservative fluxes
                do m=1,5
                    fhat_trans_cpu(j,i,k,m) = 0._rkind
                    do mm=1,5
                       fhat_trans_cpu(j,i,k,m) = fhat_trans_cpu(j,i,k,m) + er(mm,m) * ghat(mm)
                    enddo
                enddo

            endif
        enddo

    enddo
enddo
endsubroutine euler_x_fluxes_hybrid_kernel


     subroutine euler_x_fluxes_central_kernel(nv, nv_aux, nx, ny, nz, ng, &
                eul_imin, eul_imax, lmax, w_aux_trans_cpu, coeff_deriv1_cpu, dcsidx_cpu, fhat_trans_cpu, &
                force_zero_flux_min, force_zero_flux_max)
        implicit none
        integer :: nv, nv_aux, nx, ny, nz, ng
        integer :: eul_imin, eul_imax, lmax
        real(rkind), dimension(4,4) :: coeff_deriv1_cpu
        real(rkind), dimension(1-ng:ny+ng,1-ng:nx+ng,1-ng:nz+ng,1:8) :: w_aux_trans_cpu
        real(rkind), dimension(1-ng:ny+ng,1-ng:nx+ng,1-ng:nz+ng,1:nv) :: fhat_trans_cpu
        real(rkind), dimension(nx) :: dcsidx_cpu
        integer :: force_zero_flux_min, force_zero_flux_max

        ! Local variables
        integer :: i, j, k, m, l
        real(rkind) :: rhom, uui, vvi, wwi, ppi, enti, rhoi
        real(rkind) :: uuip, vvip, wwip, ppip, entip, rhoip
        real(rkind) :: ft1, ft2, ft3, ft4, ft5, ft6
        real(rkind) :: fh1, fh2, fh3, fh4, fh5
        real(rkind) :: uvs1, uvs2, uvs3, uvs4, uvs5, uvs6, uv_part
 


        do j=1,ny
do k=1,nz


        do i=eul_imin-1,eul_imax

            ft1 = 0._rkind
            ft2 = 0._rkind
            ft3 = 0._rkind
            ft4 = 0._rkind
            ft5 = 0._rkind
            ft6 = 0._rkind
            do l=1,lmax
                uvs1 = 0._rkind
                uvs2 = 0._rkind
                uvs3 = 0._rkind
                uvs4 = 0._rkind
                uvs5 = 0._rkind
                uvs6 = 0._rkind
                do m=0,l-1

                    rhoi  = w_aux_trans_cpu(j,i-m,k,1)
                    uui   = w_aux_trans_cpu(j,i-m,k,2)
                    vvi   = w_aux_trans_cpu(j,i-m,k,3)
                    wwi   = w_aux_trans_cpu(j,i-m,k,4)
                    enti  = w_aux_trans_cpu(j,i-m,k,5)
                    ppi   = w_aux_trans_cpu(j,i-m,k,6)*rhoi

                    rhoip = w_aux_trans_cpu(j,i-m+l,k,1)
                    uuip  = w_aux_trans_cpu(j,i-m+l,k,2)
                    vvip  = w_aux_trans_cpu(j,i-m+l,k,3)
                    wwip  = w_aux_trans_cpu(j,i-m+l,k,4)
                    entip = w_aux_trans_cpu(j,i-m+l,k,5)
                    ppip  = w_aux_trans_cpu(j,i-m+l,k,6)*rhoip

                    rhom  = rhoi+rhoip

                    uv_part = (uui+uuip) * rhom
                    uvs1 = uvs1 + uv_part * (2._rkind)
                    uvs2 = uvs2 + uv_part * (uui+uuip)
                    uvs3 = uvs3 + uv_part * (vvi+vvip)
                    uvs4 = uvs4 + uv_part * (wwi+wwip)
                    uvs5 = uvs5 + uv_part * (enti+entip)
                    uvs6 = uvs6 + (2._rkind)*(ppi+ppip)
                enddo
                ft1 = ft1 + coeff_deriv1_cpu(l,lmax)*uvs1
                ft2 = ft2 + coeff_deriv1_cpu(l,lmax)*uvs2
                ft3 = ft3 + coeff_deriv1_cpu(l,lmax)*uvs3
                ft4 = ft4 + coeff_deriv1_cpu(l,lmax)*uvs4
                ft5 = ft5 + coeff_deriv1_cpu(l,lmax)*uvs5
                ft6 = ft6 + coeff_deriv1_cpu(l,lmax)*uvs6
            enddo
            fh1 = 0.25_rkind*ft1
            fh2 = 0.25_rkind*ft2 
            fh3 = 0.25_rkind*ft3
            fh4 = 0.25_rkind*ft4
            fh5 = 0.25_rkind*ft5
            if((i==0 .and. force_zero_flux_min == 1).or.(i==nx .and. force_zero_flux_max == 1)) then
               fh1 = 0._rkind
               fh2 = 0._rkind
               fh3 = 0._rkind
               fh4 = 0._rkind
               fh5 = 0._rkind
            endif
            fh2 = fh2 + 0.5_rkind*ft6

            fhat_trans_cpu(j,i,k,1) = fh1
            fhat_trans_cpu(j,i,k,2) = fh2
            fhat_trans_cpu(j,i,k,3) = fh3
            fhat_trans_cpu(j,i,k,4) = fh4
            fhat_trans_cpu(j,i,k,5) = fh5
        enddo

    enddo
enddo
endsubroutine euler_x_fluxes_central_kernel


    subroutine euler_x_update_cpu(nx, ny, nz, ng, nv, eul_imin, eul_imax, fhat_trans_cpu, fl_trans_cpu, fl_cpu, dcsidx_cpu)
        integer :: nx,ny,nz,ng,nv,eul_imin,eul_imax
        real(rkind), dimension(1-ng:ny+ng,1-ng:nx+ng,1-ng:nz+ng,1:nv), intent(in) :: fhat_trans_cpu
        real(rkind), dimension(1:ny,1:nx,1:nz,1:nv), intent(inout) :: fl_trans_cpu
        real(rkind), dimension(1:nx,1:ny,1:nz,1:nv), intent(inout) :: fl_cpu
        real(rkind), dimension(1:nx), intent(in) :: dcsidx_cpu
        integer :: i,j,k,m,iv

        do k=1,nz
         do i=eul_imin,eul_imax ! loop on the inner nodes
          do j=1,ny
           do m=1,nv
            fl_trans_cpu(j,i,k,m) = (fhat_trans_cpu(j,i,k,m)-fhat_trans_cpu(j,i-1,k,m))*dcsidx_cpu(i)
           enddo
          enddo
         enddo
        enddo

        ! Note: transposed is done now alongside flux update. Hope performance are good!
        do k=1,nz
         do j=1,ny
          do i=eul_imin,eul_imax
           do iv=1,nv
            fl_cpu(i,j,k,iv) = fl_cpu(i,j,k,iv) + fl_trans_cpu(j,i,k,iv)
           enddo
          enddo
         enddo
        enddo
    endsubroutine euler_x_update_cpu

     subroutine euler_z_hybrid_kernel(nv, nv_aux, nx, ny, nz, ng, &
        eul_kmin, eul_kmax, lmax_base, nkeep, cp0, cv0, w_aux_cpu, fl_cpu, coeff_deriv1_cpu, dzitdz_cpu, fhat_cpu, &
        force_zero_flux_min, force_zero_flux_max, &
        weno_scheme, sensor_threshold, weno_size, w_cpu, gplus_z_cpu, gminus_z_cpu, cp_coeff_cpu, indx_cp_l, indx_cp_r, &
        order_modify_cpu, calorically_perfect, tol_iter_nr)
        implicit none
        ! Passed arguments
        integer :: nv, nx, ny, nz, ng, nv_aux
        integer :: eul_kmin, eul_kmax, lmax_base, nkeep, indx_cp_l, indx_cp_r, calorically_perfect
        real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv_aux) :: w_aux_cpu
        real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv) :: w_cpu
        real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv) :: fhat_cpu
        integer, dimension(0:nx,0:ny,0:nz,2:3) :: order_modify_cpu
        real(rkind), dimension(indx_cp_l:indx_cp_r) :: cp_coeff_cpu
        real(rkind), dimension(nx,ny,nz,nv) :: fl_cpu
        real(rkind), dimension(4,4) :: coeff_deriv1_cpu
        real(rkind), dimension(nz) :: dzitdz_cpu
        integer :: force_zero_flux_min, force_zero_flux_max
        real(rkind) :: sensor_threshold, cp0, cv0, tol_iter_nr
        integer :: weno_scheme, weno_size
        ! Local variables
        integer :: i, j, k, m, l
        real(rkind) :: fh1, fh2, fh3, fh4, fh5
        real(rkind) :: rhom, uui, vvi, wwi, ppi, enti, rhoi, tti
        real(rkind) :: uuip, vvip, wwip, ppip, entip, rhoip, ttip
        real(rkind) :: ft1, ft2, ft3, ft4, ft5, ft6
        real(rkind) :: uvs1, uvs2, uvs3, uvs4, uvs6, uv_part
!       real(rkind) :: uvs5
        integer :: kk, lmax, wenorec_ord
        integer :: ishk
        real(rkind) :: b1, b2, b3, c, ci, h, uu, vv, ww
        real(rkind), dimension(5,5) :: el, er
        real(rkind), dimension(5) :: ev, evmax, fk, ghat, gl, gr
        real(rkind), dimension(nv,2*weno_scheme,nx,ny) :: gplus_z_cpu
        real(rkind), dimension(nv,2*weno_scheme,nx,ny) :: gminus_z_cpu
        integer :: ll, mm
        real(rkind) :: rho, pp, wc, gc, rhow
        real(rkind) :: tt, gamloc
        real(rkind) :: uvs5_i,uvs5_k,uvs5_p,eei,eeip
        real(rkind) :: drho, dee, eem
        real(rkind) :: drhof, deef
        real(rkind) :: sumnumrho,sumnumee,sumdenrho,sumdenee
        integer :: n,n2
 


        do i=1,nx
do j=1,ny


        do k=eul_kmin-1,eul_kmax
            ishk = 0
            do kk=k-weno_scheme+1,k+weno_scheme
                if (w_aux_cpu(i,j,kk,8) > sensor_threshold) ishk = 1
            enddo

            if (ishk == 0) then
                ft1 = 0._rkind
                ft2 = 0._rkind
                ft3 = 0._rkind
                ft4 = 0._rkind
                ft5 = 0._rkind
                ft6 = 0._rkind
                lmax = lmax_base+order_modify_cpu(i,j,k,3)
                do l=1,lmax
                    uvs1 = 0._rkind
                    uvs2 = 0._rkind
                    uvs3 = 0._rkind
                    uvs4 = 0._rkind
!                   uvs5 = 0._rkind
                    uvs5_i = 0._rkind
                    uvs5_k = 0._rkind
                    uvs5_p = 0._rkind
                    uvs6 = 0._rkind
                    do m=0,l-1

                        rhoi  = w_aux_cpu(i,j,k-m,1)
                        uui   = w_aux_cpu(i,j,k-m,2)
                        vvi   = w_aux_cpu(i,j,k-m,3)
                        wwi   = w_aux_cpu(i,j,k-m,4)
                        enti  = w_aux_cpu(i,j,k-m,5)
                        tti   = w_aux_cpu(i,j,k-m,6)
                        ppi   = tti*rhoi
                        eei   = enti-tti-0.5_rkind*(uui*uui+vvi*vvi+wwi*wwi)

                        rhoip = w_aux_cpu(i,j,k-m+l,1)
                        uuip  = w_aux_cpu(i,j,k-m+l,2)
                        vvip  = w_aux_cpu(i,j,k-m+l,3)
                        wwip  = w_aux_cpu(i,j,k-m+l,4)
                        entip = w_aux_cpu(i,j,k-m+l,5)
                        ttip  = w_aux_cpu(i,j,k-m+l,6)
                        ppip  = ttip*rhoip
                        eeip  = entip-ttip-0.5_rkind*(uuip*uuip+vvip*vvip+wwip*wwip)

                        rhom  = rhoi + rhoip
                        eem   = eei  + eeip

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

                        uv_part = (wwi+wwip) * rhom * drhof
                        uvs1 = uvs1 + uv_part * (2._rkind)
                        uvs2 = uvs2 + uv_part * (uui+uuip)
                        uvs3 = uvs3 + uv_part * (vvi+vvip)
                        uvs4 = uvs4 + uv_part * (wwi+wwip)
!                       uvs5 = uvs5 + uv_part * (enti+entip)
                        uvs5_i = uvs5_i + uv_part * eem * deef
                        uvs5_k = uvs5_k + uv_part * (uui*uuip+vvi*vvip+wwi*wwip)
                        uvs5_p = uvs5_p + 4._rkind*(wwi*ppip+wwip*ppi)
                        uvs6 = uvs6 + (2._rkind)*(ppi+ppip)
                    enddo
                    ft1 = ft1 + coeff_deriv1_cpu(l,lmax)*uvs1
                    ft2 = ft2 + coeff_deriv1_cpu(l,lmax)*uvs2
                    ft3 = ft3 + coeff_deriv1_cpu(l,lmax)*uvs3
                    ft4 = ft4 + coeff_deriv1_cpu(l,lmax)*uvs4
!                   ft5 = ft5 + coeff_deriv1_cpu(l,lmax)*uvs5
                    ft5 = ft5 + coeff_deriv1_cpu(l,lmax)*(uvs5_i+uvs5_k+uvs5_p)
                    ft6 = ft6 + coeff_deriv1_cpu(l,lmax)*uvs6
                enddo
                fh1 = 0.25_rkind*ft1
                fh2 = 0.25_rkind*ft2
                fh3 = 0.25_rkind*ft3
                fh4 = 0.25_rkind*ft4
                fh5 = 0.25_rkind*ft5
                if((k==0 .and. force_zero_flux_min == 1).or.(k==nz .and. force_zero_flux_max == 1)) then
                   fh1 = 0._rkind
                   fh2 = 0._rkind
                   fh3 = 0._rkind
                   fh4 = 0._rkind
                   fh5 = 0._rkind
                endif
                fh4 = fh4 + 0.5_rkind*ft6

                fhat_cpu(i,j,k,1) = fh1
                fhat_cpu(i,j,k,2) = fh2
                fhat_cpu(i,j,k,3) = fh3
                fhat_cpu(i,j,k,4) = fh4
                fhat_cpu(i,j,k,5) = fh5
            else
                call compute_roe_average(nx, ny, nz, ng, i, i, j, j, k, k+1, w_aux_cpu, cp0, cv0, &
                                         b1, b2, b3, c, ci, h, uu, vv, ww, cp_coeff_cpu, indx_cp_l, indx_cp_r, &
                                         calorically_perfect, tol_iter_nr)

                call eigenvectors_z(b1, b2, b3, uu, vv, ww, c, ci, h, el, er)

                do m=1,5 ! loop on characteristic fields
                    evmax(m) = -1._rkind
                enddo
                do l=1,weno_size ! LLF
                    ll = k + l - weno_scheme
                    ww   = w_aux_cpu(i,j,ll,4)
                    tt   = w_aux_cpu(i,j,ll,6)
                    gamloc = get_gamloc_dev(indx_cp_l,indx_cp_r,tt,cp_coeff_cpu,calorically_perfect)
                    c    = sqrt (gamloc*tt)
                    ev(1) = abs(ww-c)
                    ev(2) = abs(ww)
                    ev(3) = abs(ww+c)
                    ev(4) = ev(2)
                    ev(5) = ev(2)
                    do m=1,5
                        evmax(m) = max(ev(m),evmax(m))
                    enddo
                enddo
                do l=1,weno_size ! loop over the stencil centered at face i
                    ll = k + l - weno_scheme

                    rho    = w_aux_cpu(i,j,ll,1)
                    uu     = w_aux_cpu(i,j,ll,2)
                    vv     = w_aux_cpu(i,j,ll,3)
                    ww     = w_aux_cpu(i,j,ll,4)
                    h      = w_aux_cpu(i,j,ll,5) 
                    rhow   = rho*ww
                    pp     = rho*w_aux_cpu(i,j,ll,6)
                    fk(1)  =      rhow
                    fk(2)  = uu * rhow
                    fk(3)  = vv * rhow 
                    fk(4)  = ww * rhow + pp
                    fk(5)  = h  * rhow
                    do m=1,5
                        wc = 0._rkind
                        gc = 0._rkind

                        do mm=1,5
                            wc = wc + el(mm,m) * w_cpu(i,j,ll,mm)
                            gc = gc + el(mm,m) * fk(mm)
                        enddo
                        gplus_z_cpu (m,l,i,j) = 0.5_rkind * (gc + evmax(m) * wc)
                        gminus_z_cpu(m,l,i,j) = gc - gplus_z_cpu(m,l,i,j)
                    enddo
                enddo
!
!               Reconstruction of the '+' and '-' fluxes
!
                wenorec_ord = weno_scheme+order_modify_cpu(i,j,k,3)
                call wenorec(nv,gplus_z_cpu(:,:,i,j),gminus_z_cpu(:,:,i,j),gl,gr,weno_scheme,wenorec_ord)
!
                do m=1,5
                    ghat(m) = gl(m) + gr(m) ! char. flux
                enddo

!               !Return to conservative fluxes
                do m=1,5
                    fhat_cpu(i,j,k,m) = 0._rkind
                    do mm=1,5
                       fhat_cpu(i,j,k,m) = fhat_cpu(i,j,k,m) + er(mm,m) * ghat(mm)
                    enddo
                enddo

            endif
        enddo

!       Update net flux 
        do k=eul_kmin,eul_kmax ! loop on the inner nodes
            do m=1,5
                fl_cpu(i,j,k,m) = fl_cpu(i,j,k,m) + (fhat_cpu(i,j,k,m)-fhat_cpu(i,j,k-1,m))*dzitdz_cpu(k)
            enddo
        enddo

    enddo
enddo
endsubroutine euler_z_hybrid_kernel


     subroutine euler_z_central_kernel(nv, nv_aux, nx, ny, nz, ng, &
        eul_kmin, eul_kmax, lmax, w_aux_cpu, fl_cpu, coeff_deriv1_cpu, dzitdz_cpu, fhat_cpu, &
        force_zero_flux_min, force_zero_flux_max)
        implicit none
        ! Passed arguments
        integer :: nv, nx, ny, nz, ng, nv_aux
        integer :: eul_kmin, eul_kmax, lmax
        real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv_aux) :: w_aux_cpu
        real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv) :: fhat_cpu
        real(rkind), dimension(nx,ny,nz,nv) :: fl_cpu
        real(rkind), dimension(4,4) :: coeff_deriv1_cpu
        real(rkind), dimension(nz) :: dzitdz_cpu
        integer :: force_zero_flux_min, force_zero_flux_max
        ! Local variables
        integer :: i, j, k, m, l
        real(rkind) :: fh1, fh2, fh3, fh4, fh5
        real(rkind) :: rhom, uui, vvi, wwi, ppi, enti, rhoi
        real(rkind) :: uuip, vvip, wwip, ppip, entip, rhoip
        real(rkind) :: ft1, ft2, ft3, ft4, ft5, ft6
        real(rkind) :: uvs1, uvs2, uvs3, uvs4, uvs5, uvs6, uv_part
 


        do i=1,nx
do j=1,ny


        do k=eul_kmin-1,eul_kmax
            ft1 = 0._rkind
            ft2 = 0._rkind
            ft3 = 0._rkind
            ft4 = 0._rkind
            ft5 = 0._rkind
            ft6 = 0._rkind
            do l=1,lmax
                uvs1 = 0._rkind
                uvs2 = 0._rkind
                uvs3 = 0._rkind
                uvs4 = 0._rkind
                uvs5 = 0._rkind
                uvs6 = 0._rkind
                do m=0,l-1

                    rhoi  = w_aux_cpu(i,j,k-m,1)
                    uui   = w_aux_cpu(i,j,k-m,2)
                    vvi   = w_aux_cpu(i,j,k-m,3)
                    wwi   = w_aux_cpu(i,j,k-m,4)
                    enti  = w_aux_cpu(i,j,k-m,5)
                    ppi   = w_aux_cpu(i,j,k-m,6)*rhoi

                    rhoip = w_aux_cpu(i,j,k-m+l,1)
                    uuip  = w_aux_cpu(i,j,k-m+l,2)
                    vvip  = w_aux_cpu(i,j,k-m+l,3)
                    wwip  = w_aux_cpu(i,j,k-m+l,4)
                    entip = w_aux_cpu(i,j,k-m+l,5)
                    ppip  = w_aux_cpu(i,j,k-m+l,6)*rhoip

                    rhom  = rhoi+rhoip

                    uv_part = (wwi+wwip) * rhom
                    uvs1 = uvs1 + uv_part * (2._rkind)
                    uvs2 = uvs2 + uv_part * (uui+uuip)
                    uvs3 = uvs3 + uv_part * (vvi+vvip)
                    uvs4 = uvs4 + uv_part * (wwi+wwip)
                    uvs5 = uvs5 + uv_part * (enti+entip)
                    uvs6 = uvs6 + (2._rkind)*(ppi+ppip)
                enddo
                ft1 = ft1 + coeff_deriv1_cpu(l,lmax)*uvs1
                ft2 = ft2 + coeff_deriv1_cpu(l,lmax)*uvs2
                ft3 = ft3 + coeff_deriv1_cpu(l,lmax)*uvs3
                ft4 = ft4 + coeff_deriv1_cpu(l,lmax)*uvs4
                ft5 = ft5 + coeff_deriv1_cpu(l,lmax)*uvs5
                ft6 = ft6 + coeff_deriv1_cpu(l,lmax)*uvs6
            enddo
            fh1 = 0.25_rkind*ft1
            fh2 = 0.25_rkind*ft2
            fh3 = 0.25_rkind*ft3
            fh4 = 0.25_rkind*ft4
            fh5 = 0.25_rkind*ft5
            if((k==0 .and. force_zero_flux_min == 1).or.(k==nz .and. force_zero_flux_max == 1)) then
               fh1 = 0._rkind
               fh2 = 0._rkind
               fh3 = 0._rkind
               fh4 = 0._rkind
               fh5 = 0._rkind
            endif
            fh4 = fh4 + 0.5_rkind*ft6

            fhat_cpu(i,j,k,1) = fh1
            fhat_cpu(i,j,k,2) = fh2
            fhat_cpu(i,j,k,3) = fh3
            fhat_cpu(i,j,k,4) = fh4
            fhat_cpu(i,j,k,5) = fh5
        enddo

!       Update net flux 
        do k=eul_kmin,eul_kmax ! loop on the inner nodes
            do m=1,5
                fl_cpu(i,j,k,m) = fl_cpu(i,j,k,m) + (fhat_cpu(i,j,k,m)-fhat_cpu(i,j,k-1,m))*dzitdz_cpu(k)
            enddo
        enddo

    enddo
enddo
endsubroutine euler_z_central_kernel


     subroutine euler_y_central_kernel(nv, nv_aux, nx, ny, nz, ng, &
        eul_jmin, eul_jmax, lmax, w_aux_cpu, fl_cpu, coeff_deriv1_cpu, detady_cpu, fhat_cpu, &
        force_zero_flux_min, force_zero_flux_max)
        implicit none
        ! Passed arguments
        integer :: nv, nx, ny, nz, ng, nv_aux
        integer :: eul_jmin, eul_jmax, lmax
        real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv_aux) :: w_aux_cpu
        real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv) :: fhat_cpu
        real(rkind), dimension(nx,ny,nz,nv) :: fl_cpu
        real(rkind), dimension(4,4) :: coeff_deriv1_cpu
        real(rkind), dimension(ny) :: detady_cpu
        integer :: force_zero_flux_min, force_zero_flux_max
        ! Local variables
        integer :: i, j, k, m, l
        real(rkind) :: fh1, fh2, fh3, fh4, fh5
        real(rkind) :: rhom, uui, vvi, wwi, ppi, enti, rhoi
        real(rkind) :: uuip, vvip, wwip, ppip, entip, rhoip
        real(rkind) :: ft1, ft2, ft3, ft4, ft5, ft6
        real(rkind) :: uvs1, uvs2, uvs3, uvs4, uvs5, uvs6, uv_part
 


        do i=1,nx
do k=1,nz


        do j=eul_jmin-1,eul_jmax
            ft1 = 0._rkind
            ft2 = 0._rkind
            ft3 = 0._rkind
            ft4 = 0._rkind
            ft5 = 0._rkind
            ft6 = 0._rkind

            !lmax = scheme_cpu(i,j,k)

            do l=1,lmax
                uvs1 = 0._rkind
                uvs2 = 0._rkind
                uvs3 = 0._rkind
                uvs4 = 0._rkind
                uvs5 = 0._rkind
                uvs6 = 0._rkind
                do m=0,l-1

                    rhoi  = w_aux_cpu(i,j-m,k,1)
                    uui   = w_aux_cpu(i,j-m,k,2)
                    vvi   = w_aux_cpu(i,j-m,k,3)
                    wwi   = w_aux_cpu(i,j-m,k,4)
                    enti  = w_aux_cpu(i,j-m,k,5)
                    ppi   = w_aux_cpu(i,j-m,k,6)*rhoi

                    rhoip = w_aux_cpu(i,j-m+l,k,1)
                    uuip  = w_aux_cpu(i,j-m+l,k,2)
                    vvip  = w_aux_cpu(i,j-m+l,k,3)
                    wwip  = w_aux_cpu(i,j-m+l,k,4)
                    entip = w_aux_cpu(i,j-m+l,k,5)
                    ppip  = w_aux_cpu(i,j-m+l,k,6)*rhoip

                    rhom  = rhoi+rhoip

                    uv_part = (vvi+vvip) * rhom
                    uvs1 = uvs1 + uv_part * (2._rkind)
                    uvs2 = uvs2 + uv_part * (uui+uuip)
                    uvs3 = uvs3 + uv_part * (vvi+vvip)
                    uvs4 = uvs4 + uv_part * (wwi+wwip)
                    uvs5 = uvs5 + uv_part * (enti+entip)
                    uvs6 = uvs6 + (2._rkind)*(ppi+ppip)
                enddo
                ft1 = ft1 + coeff_deriv1_cpu(l,lmax)*uvs1
                ft2 = ft2 + coeff_deriv1_cpu(l,lmax)*uvs2
                ft3 = ft3 + coeff_deriv1_cpu(l,lmax)*uvs3
                ft4 = ft4 + coeff_deriv1_cpu(l,lmax)*uvs4
                ft5 = ft5 + coeff_deriv1_cpu(l,lmax)*uvs5
                ft6 = ft6 + coeff_deriv1_cpu(l,lmax)*uvs6
            enddo
            fh1 = 0.25_rkind*ft1
            fh2 = 0.25_rkind*ft2
            fh3 = 0.25_rkind*ft3
            fh4 = 0.25_rkind*ft4
            fh5 = 0.25_rkind*ft5
            if((j==0 .and. force_zero_flux_min == 1).or.(j==ny .and. force_zero_flux_max == 1)) then
               fh1 = 0._rkind
               fh2 = 0._rkind
               fh3 = 0._rkind
               fh4 = 0._rkind
               fh5 = 0._rkind
            endif

            fh3 = fh3 + 0.5_rkind*ft6

            fhat_cpu(i,j,k,1) = fh1
            fhat_cpu(i,j,k,2) = fh2
            fhat_cpu(i,j,k,3) = fh3
            fhat_cpu(i,j,k,4) = fh4
            fhat_cpu(i,j,k,5) = fh5
        enddo

!       Update net flux 
        do j=eul_jmin,eul_jmax ! loop on the inner nodes
            do m=1,5
                fl_cpu(i,j,k,m) = fl_cpu(i,j,k,m) + (fhat_cpu(i,j,k,m)-fhat_cpu(i,j-1,k,m))*detady_cpu(j)
            enddo
        enddo

    enddo
enddo
endsubroutine euler_y_central_kernel


     subroutine euler_y_hybrid_kernel(nv, nv_aux, nx, ny, nz, ng, &
        eul_jmin, eul_jmax, lmax_base, nkeep, cp0, cv0, w_aux_cpu, fl_cpu, coeff_deriv1_cpu, detady_cpu, fhat_cpu, &
        force_zero_flux_min, force_zero_flux_max, &
        weno_scheme, sensor_threshold, weno_size, w_cpu, gplus_y_cpu, gminus_y_cpu, cp_coeff_cpu, indx_cp_l, indx_cp_r, &
        order_modify_cpu, calorically_perfect, tol_iter_nr)
        implicit none
        ! Passed arguments
        integer :: nv, nx, ny, nz, ng, nv_aux
        integer :: eul_jmin, eul_jmax, lmax_base, nkeep, indx_cp_l, indx_cp_r, calorically_perfect
        real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv_aux) :: w_aux_cpu
        real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv) :: w_cpu
        real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv) :: fhat_cpu
        integer, dimension(0:nx,0:ny,0:nz,2:3) :: order_modify_cpu
        real(rkind), dimension(indx_cp_l:indx_cp_r) :: cp_coeff_cpu
        real(rkind), dimension(nx,ny,nz,nv) :: fl_cpu
        real(rkind), dimension(4,4) :: coeff_deriv1_cpu
        real(rkind), dimension(ny) :: detady_cpu
        integer :: force_zero_flux_min, force_zero_flux_max
        real(rkind) :: sensor_threshold, cp0, cv0, tol_iter_nr
        integer :: weno_scheme, weno_size
        ! Local variables
        integer :: i, j, k, m, l
        real(rkind) :: fh1, fh2, fh3, fh4, fh5
        real(rkind) :: rhom, uui, vvi, wwi, ppi, enti, rhoi, tti
        real(rkind) :: uuip, vvip, wwip, ppip, entip, rhoip, ttip
        real(rkind) :: ft1, ft2, ft3, ft4, ft5, ft6
        real(rkind) :: uvs1, uvs2, uvs3, uvs4, uvs6, uv_part
!       real(rkind) :: uvs5
        integer :: jj
        integer :: ishk
        real(rkind) :: b1, b2, b3, c, ci, h, uu, vv, ww
        real(rkind), dimension(5,5) :: el, er
        real(rkind), dimension(5) :: ev, evmax, fj, ghat, gl, gr
        real(rkind), dimension(nv,2*weno_scheme,nx,nz) :: gplus_y_cpu
        real(rkind), dimension(nv,2*weno_scheme,nx,nz) :: gminus_y_cpu
        integer :: ll, mm, lmax, wenorec_ord
        real(rkind) :: rho, pp, wc, gc, rhov
        real(rkind) :: tt, gamloc
        real(rkind) :: uvs5_i,uvs5_k,uvs5_p,eei,eeip
        real(rkind) :: drho, dee, eem
        real(rkind) :: drhof, deef
        real(rkind) :: sumnumrho,sumnumee,sumdenrho,sumdenee
        integer :: n,n2
 


        do i=1,nx
do k=1,nz


        do j=eul_jmin-1,eul_jmax

            ishk = 0
            do jj=j-weno_scheme+1,j+weno_scheme
                if (w_aux_cpu(i,jj,k,8) > sensor_threshold) ishk = 1
            enddo

            if (ishk == 0) then
                ft1 = 0._rkind
                ft2 = 0._rkind
                ft3 = 0._rkind
                ft4 = 0._rkind
                ft5 = 0._rkind
                ft6 = 0._rkind

                lmax = lmax_base+order_modify_cpu(i,j,k,2)

                do l=1,lmax
                    uvs1 = 0._rkind
                    uvs2 = 0._rkind
                    uvs3 = 0._rkind
                    uvs4 = 0._rkind
!                   uvs5 = 0._rkind
                    uvs5_i = 0._rkind
                    uvs5_k = 0._rkind
                    uvs5_p = 0._rkind
                    uvs6 = 0._rkind
                    do m=0,l-1

                        rhoi  = w_aux_cpu(i,j-m,k,1)
                        uui   = w_aux_cpu(i,j-m,k,2)
                        vvi   = w_aux_cpu(i,j-m,k,3)
                        wwi   = w_aux_cpu(i,j-m,k,4)
                        enti  = w_aux_cpu(i,j-m,k,5)
                        tti   = w_aux_cpu(i,j-m,k,6)
                        ppi   = tti*rhoi
                        eei   = enti-tti-0.5_rkind*(uui*uui+vvi*vvi+wwi*wwi)

                        rhoip = w_aux_cpu(i,j-m+l,k,1)
                        uuip  = w_aux_cpu(i,j-m+l,k,2)
                        vvip  = w_aux_cpu(i,j-m+l,k,3)
                        wwip  = w_aux_cpu(i,j-m+l,k,4)
                        entip = w_aux_cpu(i,j-m+l,k,5)
                        ttip  = w_aux_cpu(i,j-m+l,k,6)
                        ppip  = ttip*rhoip
                        eeip  = entip-ttip-0.5_rkind*(uuip*uuip+vvip*vvip+wwip*wwip)

                        rhom  = rhoi + rhoip
                        eem   = eei  + eeip

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

                        uv_part = (vvi+vvip) * rhom * drhof
                        uvs1 = uvs1 + uv_part * (2._rkind)
                        uvs2 = uvs2 + uv_part * (uui+uuip)
                        uvs3 = uvs3 + uv_part * (vvi+vvip)
                        uvs4 = uvs4 + uv_part * (wwi+wwip)
!                       uvs5 = uvs5 + uv_part * (enti+entip)
                        uvs5_i = uvs5_i + uv_part * eem * deef
                        uvs5_k = uvs5_k + uv_part * (uui*uuip+vvi*vvip+wwi*wwip)
                        uvs5_p = uvs5_p + 4._rkind*(vvi*ppip+vvip*ppi)
                        uvs6 = uvs6 + (2._rkind)*(ppi+ppip)
                    enddo
                    ft1 = ft1 + coeff_deriv1_cpu(l,lmax)*uvs1
                    ft2 = ft2 + coeff_deriv1_cpu(l,lmax)*uvs2
                    ft3 = ft3 + coeff_deriv1_cpu(l,lmax)*uvs3
                    ft4 = ft4 + coeff_deriv1_cpu(l,lmax)*uvs4
!                   ft5 = ft5 + coeff_deriv1_cpu(l,lmax)*uvs5
                    ft5 = ft5 + coeff_deriv1_cpu(l,lmax)*(uvs5_i+uvs5_k+uvs5_p)
                    ft6 = ft6 + coeff_deriv1_cpu(l,lmax)*uvs6
                enddo
                fh1 = 0.25_rkind*ft1
                fh2 = 0.25_rkind*ft2
                fh3 = 0.25_rkind*ft3
                fh4 = 0.25_rkind*ft4
                fh5 = 0.25_rkind*ft5
                if((j==0 .and. force_zero_flux_min == 1).or.(j==ny .and. force_zero_flux_max == 1)) then
                   fh1 = 0._rkind
                   fh2 = 0._rkind
                   fh3 = 0._rkind
                   fh4 = 0._rkind
                   fh5 = 0._rkind
                endif

                fh3 = fh3 + 0.5_rkind*ft6

                fhat_cpu(i,j,k,1) = fh1
                fhat_cpu(i,j,k,2) = fh2
                fhat_cpu(i,j,k,3) = fh3
                fhat_cpu(i,j,k,4) = fh4
                fhat_cpu(i,j,k,5) = fh5
            else
                call compute_roe_average(nx, ny, nz, ng, i, i, j, j+1, k, k, w_aux_cpu, cp0, cv0, &
                                         b1, b2, b3, c, ci, h, uu, vv, ww, cp_coeff_cpu, indx_cp_l, indx_cp_r, &
                                         calorically_perfect, tol_iter_nr)

                call eigenvectors_y(b1, b2, b3, uu, vv, ww, c, ci, h, el, er)

                do m=1,5 ! loop on characteristic fields
                    evmax(m) = -1._rkind
                enddo
                do l=1,weno_size ! LLF
                    ll = j + l - weno_scheme
                    vv   = w_aux_cpu(i,ll,k,3)
                    tt   = w_aux_cpu(i,ll,k,6)
                    gamloc = get_gamloc_dev(indx_cp_l,indx_cp_r,tt,cp_coeff_cpu,calorically_perfect)
                    c    = sqrt (gamloc*tt)
                    ev(1) = abs(vv-c)
                    ev(2) = abs(vv)
                    ev(3) = abs(vv+c)
                    ev(4) = ev(2)
                    ev(5) = ev(2)
                    do m=1,5
                        evmax(m) = max(ev(m),evmax(m))
                    enddo
                enddo
                do l=1,weno_size ! loop over the stencil centered at face i
                    ll = j + l - weno_scheme

                    rho    = w_aux_cpu(i,ll,k,1)
                    uu     = w_aux_cpu(i,ll,k,2)
                    vv     = w_aux_cpu(i,ll,k,3)
                    ww     = w_aux_cpu(i,ll,k,4)
                    h      = w_aux_cpu(i,ll,k,5) 
                    rhov   = rho*vv
                    pp     = rho*w_aux_cpu(i,ll,k,6)
                    fj(1)  =      rhov
                    fj(2)  = uu * rhov
                    fj(3)  = vv * rhov + pp
                    fj(4)  = ww * rhov
                    fj(5)  = h  * rhov
                    do m=1,5
                        wc = 0._rkind
                        gc = 0._rkind

                        do mm=1,5
                            wc = wc + el(mm,m) * w_cpu(i,ll,k,mm)
                            gc = gc + el(mm,m) * fj(mm)
                        enddo
                        gplus_y_cpu (m,l,i,k) = 0.5_rkind * (gc + evmax(m) * wc)
                        gminus_y_cpu(m,l,i,k) = gc - gplus_y_cpu(m,l,i,k)
                    enddo
                enddo
!
!               Reconstruction of the '+' and '-' fluxes
!
                wenorec_ord = weno_scheme+order_modify_cpu(i,j,k,2)
                call wenorec(nv,gplus_y_cpu(:,:,i,k),gminus_y_cpu(:,:,i,k),gl,gr,weno_scheme,wenorec_ord)
!
                do m=1,5
                    ghat(m) = gl(m) + gr(m) ! char. flux
                enddo

!               !Return to conservative fluxes
                do m=1,5
                    fhat_cpu(i,j,k,m) = 0._rkind
                    do mm=1,5
                       fhat_cpu(i,j,k,m) = fhat_cpu(i,j,k,m) + er(mm,m) * ghat(mm)
                    enddo
                enddo

            endif
        enddo

!       Update net flux 
        do j=eul_jmin,eul_jmax ! loop on the inner nodes
            do m=1,5
                fl_cpu(i,j,k,m) = fl_cpu(i,j,k,m) + (fhat_cpu(i,j,k,m)-fhat_cpu(i,j-1,k,m))*detady_cpu(j)
            enddo
        enddo

    enddo
enddo
endsubroutine euler_y_hybrid_kernel


     subroutine wenorec(nvar,vp,vm,vminus,vplus,iweno,wenorec_ord)
    
    !    Passed arguments
         integer :: nvar, iweno, wenorec_ord
         real(rkind), dimension(nvar,2*iweno) :: vm,vp
         real(rkind), dimension(nvar) :: vminus,vplus
    
    !    Local variables
         real(rkind), dimension(-1:4) :: dwe           ! linear weights
         real(rkind), dimension(-1:4) :: alfp,alfm     ! alpha_l
!        real(rkind), dimension(-1:4) :: alfp_map,alfm_map ! alpha_l
         real(rkind), dimension(-1:4) :: betap,betam   ! beta_l
         real(rkind), dimension(-1:4) :: omp,omm       ! WENO weights
    !    
         integer :: i,l,m
         real(rkind) :: c0,c1,c2,c3,c4,d0,d1,d2,d3,summ,sump
    !    
         if (wenorec_ord==1) then ! Godunov
    !    
             i = iweno ! index of intermediate node to perform reconstruction
    !    
             vminus(1:nvar) = vp(1:nvar,i)
             vplus (1:nvar) = vm(1:nvar,i+1)
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
                 betap(0)  = (vp(m,i  )-vp(m,i-1))**2
                 betap(1)  = (vp(m,i+1)-vp(m,i  ))**2
                 betam(0)  = (vm(m,i+2)-vm(m,i+1))**2
                 betam(1)  = (vm(m,i+1)-vm(m,i  ))**2
    !    
                 sump = 0._rkind
                 summ = 0._rkind
                 do l=0,1
                     alfp(l) = dwe(l)/(0.000001_rkind+betap(l))**2
                     alfm(l) = dwe(l)/(0.000001_rkind+betam(l))**2
                     sump = sump + alfp(l)
                     summ = summ + alfm(l)
                 enddo
                 do l=0,1
                     omp(l) = alfp(l)/sump
                     omm(l) = alfm(l)/summ
                 enddo
    !    
                 vminus(m) = omp(0) *(-vp(m,i-1)+3*vp(m,i  )) + omp(1) *( vp(m,i  )+ vp(m,i+1))
                 vplus(m)  = omm(0) *(-vm(m,i+2)+3*vm(m,i+1)) + omm(1) *( vm(m,i  )+ vm(m,i+1))
    !    
             enddo ! end of m-loop
    !    
             do m=1,nvar
                 vminus(m) = 0.5_rkind*vminus(m)
                 vplus(m)  = 0.5_rkind*vplus(m)
             enddo
    !    
         elseif (wenorec_ord==3) then ! WENO-5
    !    
          i = iweno ! index of intermediate node to perform reconstruction
    !    
          dwe( 0) = 1._rkind/10._rkind
          dwe( 1) = 6._rkind/10._rkind
          dwe( 2) = 3._rkind/10._rkind
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
          do m=1,nvar
    !    
           betap(2) = d0*(     vp(m,i)-2._rkind*vp(m,i+1)+vp(m,i+2))**2+d1*(3._rkind*vp(m,i)-4._rkind*vp(m,i+1)+vp(m,i+2))**2
           betap(1) = d0*(     vp(m,i-1)-2._rkind*vp(m,i)+vp(m,i+1))**2+d1*(     vp(m,i-1)-vp(m,i+1) )**2
           betap(0) = d0*(     vp(m,i)-2._rkind*vp(m,i-1)+vp(m,i-2))**2+d1*(3._rkind*vp(m,i)-4._rkind*vp(m,i-1)+vp(m,i-2))**2
    !    
           betam(2) = d0*(     vm(m,i+1)-2._rkind*vm(m,i)+vm(m,i-1))**2+d1*(3._rkind*vm(m,i+1)-4._rkind*vm(m,i)+vm(m,i-1))**2
           betam(1) = d0*(     vm(m,i+2)-2._rkind*vm(m,i+1)+vm(m,i))**2+d1*(     vm(m,i+2)-vm(m,i) )**2
           betam(0) = d0*(     vm(m,i+1)-2._rkind*vm(m,i+2)+vm(m,i+3))**2+d1*(3._rkind*vm(m,i+1)-4._rkind*vm(m,i+2)+vm(m,i+3))**2
    !    
           sump = 0._rkind
           summ = 0._rkind
           do l=0,2
            alfp(l) = dwe(  l)/(0.000001_rkind+betap(l))**2
            alfm(l) = dwe(  l)/(0.000001_rkind+betam(l))**2
            sump = sump + alfp(l)
            summ = summ + alfm(l)
           enddo
           do l=0,2
            omp(l) = alfp(l)/sump
            omm(l) = alfm(l)/summ
           enddo
    !
           vminus(m)   = omp(2)*(c0*vp(m,i  )+c1*vp(m,i+1)+c2*vp(m,i+2)) + &
             & omp(1)*(c2*vp(m,i-1)+c1*vp(m,i  )+c0*vp(m,i+1)) + omp(0)*(c0*vp(m,i-2)+c3*vp(m,i-1)+c4*vp(m,i  ))
           vplus(m)   = omm(2)*(c0*vm(m,i+1)+c1*vm(m,i  )+c2*vm(m,i-1)) +  &
             & omm(1)*(c2*vm(m,i+2)+c1*vm(m,i+1)+c0*vm(m,i  )) + omm(0)*(c0*vm(m,i+3)+c3*vm(m,i+2)+c4*vm(m,i+1))
    !    
          enddo ! end of m-loop 
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
           betap(3)= d1*(-11*vp(m,  i)+18*vp(m,i+1)- 9*vp(m,i+2)+ 2*vp(m,i+3))**2+&
           &  d2*(  2*vp(m,  i)- 5*vp(m,i+1)+ 4*vp(m,i+2)-   vp(m,i+3))**2+ &
           & d3*(   -vp(m,  i)+ 3*vp(m,i+1)- 3*vp(m,i+2)+   vp(m,i+3))**2
           betap(2)= d1*(- 2*vp(m,i-1)- 3*vp(m,i  )+ 6*vp(m,i+1)-   vp(m,i+2))**2+&
           &  d2*(    vp(m,i-1)- 2*vp(m,i  )+   vp(m,i+1)             )**2+&
           &  d3*(   -vp(m,i-1)+ 3*vp(m,i  )- 3*vp(m,i+1)+   vp(m,i+2))**2
           betap(1)= d1*(    vp(m,i-2)- 6*vp(m,i-1)+ 3*vp(m,i  )+ 2*vp(m,i+1))**2+&
           &  d2*( vp(m,i-1)- 2*vp(m,i  )+   vp(m,i+1))**2+ &
           &  d3*(   -vp(m,i-2)+ 3*vp(m,i-1)- 3*vp(m,i  )+   vp(m,i+1))**2
           betap(0)= d1*(- 2*vp(m,i-3)+ 9*vp(m,i-2)-18*vp(m,i-1)+11*vp(m,i  ))**2+&
           &  d2*(-   vp(m,i-3)+ 4*vp(m,i-2)- 5*vp(m,i-1)+ 2*vp(m,i  ))**2+&
           &  d3*(   -vp(m,i-3)+ 3*vp(m,i-2)- 3*vp(m,i-1)+   vp(m,i  ))**2
    !    
           betam(3)= d1*(-11*vm(m,i+1)+18*vm(m,i  )- 9*vm(m,i-1)+ 2*vm(m,i-2))**2+&
           &  d2*(  2*vm(m,i+1)- 5*vm(m,i  )+ 4*vm(m,i-1)-   vm(m,i-2))**2+&
           &  d3*(   -vm(m,i+1)+ 3*vm(m,i  )- 3*vm(m,i-1)+   vm(m,i-2))**2
           betam(2)= d1*(- 2*vm(m,i+2)- 3*vm(m,i+1)+ 6*vm(m,i  )-   vm(m,i-1))**2+&
           &  d2*(    vm(m,i+2)- 2*vm(m,i+1)+   vm(m,i  )             )**2+&
           &  d3*(   -vm(m,i+2)+ 3*vm(m,i+1)- 3*vm(m,i  )+   vm(m,i-1))**2
           betam(1)= d1*(    vm(m,i+3)- 6*vm(m,i+2)+ 3*vm(m,i+1)+ 2*vm(m,i  ))**2+&
           &  d2*(                 vm(m,i+2)- 2*vm(m,i+1)+   vm(m,i  ))**2+&
           &  d3*(   -vm(m,i+3)+ 3*vm(m,i+2)- 3*vm(m,i+1)+   vm(m,i  ))**2
           betam(0)= d1*(- 2*vm(m,i+4)+ 9*vm(m,i+3)-18*vm(m,i+2)+11*vm(m,i+1))**2+&
           &  d2*(-   vm(m,i+4)+ 4*vm(m,i+3)- 5*vm(m,i+2)+ 2*vm(m,i+1))**2+&
           &  d3*(   -vm(m,i+4)+ 3*vm(m,i+3)- 3*vm(m,i+2)+   vm(m,i+1))**2 
    !    
           sump = 0._rkind
           summ = 0._rkind
           do l=0,3
            alfp(l) = dwe(  l)/(0.000001_rkind+betap(l))**2
            alfm(l) = dwe(  l)/(0.000001_rkind+betam(l))**2
            sump = sump + alfp(l)
            summ = summ + alfm(l)
           enddo
           do l=0,3
            omp(l) = alfp(l)/sump
            omm(l) = alfm(l)/summ
           enddo
    !    
           vminus(m)   = omp(3)*( 6*vp(m,i  )+26*vp(m,i+1)-10*vp(m,i+2)+ 2*vp(m,i+3))+&
            omp(2)*(-2*vp(m,i-1)+14*vp(m,i  )+14*vp(m,i+1)- 2*vp(m,i+2))+&
            omp(1)*( 2*vp(m,i-2)-10*vp(m,i-1)+26*vp(m,i  )+ 6*vp(m,i+1))+&
            omp(0)*(-6*vp(m,i-3)+26*vp(m,i-2)-46*vp(m,i-1)+50*vp(m,i  ))
           vplus(m)   =  omm(3)*( 6*vm(m,i+1)+26*vm(m,i  )-10*vm(m,i-1)+ 2*vm(m,i-2))+&
            omm(2)*(-2*vm(m,i+2)+14*vm(m,i+1)+14*vm(m,i  )- 2*vm(m,i-1))+&
            omm(1)*( 2*vm(m,i+3)-10*vm(m,i+2)+26*vm(m,i+1)+ 6*vm(m,i  ))+&
            omm(0)*(-6*vm(m,i+4)+26*vm(m,i+3)-46*vm(m,i+2)+50*vm(m,i+1))
    !    
          enddo ! end of m-loop 
    !    
          vminus = vminus/24._rkind
          vplus  = vplus /24._rkind
    !    
         else
          write(*,*) 'Error! WENO scheme not implemented'
          stop
         endif
    
    endsubroutine wenorec
!
     subroutine eigenvectors_x(b1, b2, b3, uu, vv, ww, c, ci, h, el, er)
        real(rkind), intent(in) :: b1, b2, b3, uu, vv, ww, c, ci, h
        real(rkind), dimension(:,:), intent(out) ::  el, er

!       left eigenvectors matrix (at Roe state)
!       matrix L-1 AIAA2001-2609 - Rohde after transposition
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

!       right eigenvectors matrix (at Roe state)
!       matrix R-1 AIAA2001-2609 - Rohde after transposition
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

    endsubroutine eigenvectors_x

     subroutine eigenvectors_y(b1, b2, b3, uu, vv, ww, c, ci, h, el, er)
        real(rkind), intent(in) :: b1, b2, b3, uu, vv, ww, c, ci, h
        real(rkind), dimension(:,:), intent(out) ::  el, er

!       left eigenvectors matrix (at Roe state)
!       matrix L-2 AIAA2001-2609 - Rohde after transposition
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

!       right eigenvectors matrix (at Roe state)
!       matrix R-2 AIAA2001-2609 - Rohde after transposition
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

     subroutine eigenvectors_z(b1, b2, b3, uu, vv, ww, c, ci, h, el, er)
        real(rkind), intent(in) :: b1, b2, b3, uu, vv, ww, c, ci, h
        real(rkind), dimension(:,:), intent(out) ::  el, er

!       left eigenvectors matrix (at Roe state)
!       matrix L-3 AIAA2001-2609 - Rohde after transposition
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

!       right eigenvectors matrix (at Roe state)
!       matrix R-3 AIAA2001-2609 - Rohde after transposition
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


    endsubroutine eigenvectors_z

     subroutine compute_roe_average(nx, ny, nz, ng, i, ip, j, jp, k, kp, w_aux_cpu, cp0, cv0, &
                                                      b1, b2, b3, c, ci, h, uu, vv, ww, cp_coeff_cpu, indx_cp_l, indx_cp_r,&
                                                      calorically_perfect,tol_iter_nr)
        integer :: ng,i,ip,j,jp,k,kp,indx_cp_l,indx_cp_r,calorically_perfect
        integer :: nx,ny,nz
        real(rkind), intent(in) :: cp0, cv0, tol_iter_nr
        real(rkind), intent(out) :: b1, b2, b3, uu, vv, ww, ci, h, c
        real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in) :: w_aux_cpu 
        real(rkind), dimension(indx_cp_l:), intent(in) :: cp_coeff_cpu
        real(rkind) :: up, vp, wp, qqp, hp, r, rp1, cc, qq, gam, gm1
        integer :: ll,iter,max_iter
        real(rkind) :: tt,gm1loc,hbar,ttp,told,num,den,gamloc,cploc,tpow,tpowp,p_rho,p_e,etot,rho
!
        max_iter = 50
!
        gam = cp0/cv0
        gm1 = gam-1._rkind
!       Compute Roe average
!       Left state (node i)
        !ri        =  1._rkind/w_aux_cpu(im,jm,km,1)
        uu        =  w_aux_cpu(i,j,k,2)
        vv        =  w_aux_cpu(i,j,k,3)
        ww        =  w_aux_cpu(i,j,k,4)
        qq        =  0.5_rkind * (uu*uu  +vv*vv + ww*ww)
!       pp        =  gm1 * (w_cpu(i,j,k,5) - w_cpu(i,j,k,1) * qq)
!       h         =  (w_cpu(i,j,k,5)  + pp) * ri
        h         =  w_aux_cpu(i,j,k,5) ! gamma*w_cpu(im,jm,km,5)*ri-gm1*qq
!       Right state (node i+1)
        !rip       =  1._rkind/w_aux_cpu(ip,jp,kp,1)
        up        =  w_aux_cpu(ip,jp,kp,2)
        vp        =  w_aux_cpu(ip,jp,kp,3)
        wp        =  w_aux_cpu(ip,jp,kp,4)
        qqp       =  0.5_rkind * (up*up  +vp*vp +wp*wp)
!       ppp       =  gm1 * (w_cpu(i,jp,k,5) - w_cpu(i,jp,k,1) * qqp)
!       hp        =  (w_cpu(i,jp,k,5)  + ppp) * rip
        hp        =  w_aux_cpu(ip,jp,kp,5) ! gamma*w_cpu(ip,jp,kp,5)*rip-gm1*qqp
!       average state
        r         =  w_aux_cpu(ip,jp,kp,1)/w_aux_cpu(i,j,k,1)
        r         =  sqrt(r)
        rho       =  r*w_aux_cpu(i,j,k,1)
        rp1       =  1._rkind/(r  +1._rkind)
        uu        =  (r*up  +uu)*rp1
        vv        =  (r*vp  +vv)*rp1
        ww        =  (r*wp  +ww)*rp1
        h         =  (r*hp  +h)*rp1
        qq        =  0.5_rkind * (uu*uu  +vv*vv + ww*ww)
!
        if (calorically_perfect==1) then
         cc       = gm1 * (h - qq) 
         tt       = cc/gam
         gm1loc   = gm1
        else
         hbar      =  (h - qq) - cp0
         tt        =  w_aux_cpu(i ,j ,k ,6)
         ttp       =  w_aux_cpu(ip,jp,kp,6)
         tt        =  (r*ttp +tt)*rp1
         told      = tt ! First attempt computed as Roe average of temperature
         do iter=1,max_iter
          num = 0._rkind
          den = 0._rkind
          do ll=indx_cp_l,indx_cp_r
           if (ll==-1) then
            tpow  = told**ll
            den = den+cp_coeff_cpu(ll)*tpow
            num = num+cp_coeff_cpu(ll)*log(told)
           else
            tpow  = told**ll
            tpowp = told*tpow
            den = den+cp_coeff_cpu(ll)*tpow
            num = num+cp_coeff_cpu(ll)*(tpowp-1._rkind)/(ll+1._rkind)
           endif
          enddo
          tt = told+(hbar-num)/den
          if (abs(tt-told)<tol_iter_nr) exit
           told = tt
         enddo
         cploc = 0._rkind
         do ll=indx_cp_l,indx_cp_r
          cploc = cploc+cp_coeff_cpu(ll)*tt**ll
         enddo
         gamloc = cploc/(cploc-1._rkind)
         gm1loc = gamloc-1._rkind
         cc     =  gamloc*tt
        endif
!
        c         =  sqrt(cc) 
        ci        =  1._rkind/c 
!
        p_rho     = tt
        p_e       = rho*gm1loc   ! p_e   = rho * R_gas / Cv
        etot      = h - tt
!
        b3        = etot - rho * p_rho/p_e ! qq in the case of calorically perfect
        b2        = p_e/(rho*cc)
        b1        = p_rho/cc - b2*(etot - 2._rkind*qq)
!
!       b2        = gm1/cc  ! 1/(cp*T)
!       b1        = b2 * qq
!       b3        = qq 
!       b2        = gm1/cc
!       b1        = b2 * b3

    endsubroutine compute_roe_average

    subroutine force_rhs_2_cpu(nx, ny, nz, ng, fln_cpu, w_aux_cpu, bulk5g, fluid_mask_cpu)
        integer :: nx, ny, nz, ng
        real(rkind), dimension(5), intent(in) :: bulk5g
        real(rkind), dimension(1:,1:,1:,1:), intent(inout) :: fln_cpu
        real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in) :: w_aux_cpu
        integer, dimension(1-ng:,1-ng:,1-ng:), intent(in) :: fluid_mask_cpu
        real(rkind) :: bulk_1, bulk_2, uu
        integer :: i,j,k
        bulk_1 = bulk5g(1)
        bulk_2 = bulk5g(2)
        do k=1,nz
         do j=1,ny
          do i=1,nx
           if (fluid_mask_cpu(i,j,k)==0) then
            uu = w_aux_cpu(i,j,k,2)
            fln_cpu(i,j,k,1) = fln_cpu(i,j,k,1) -    bulk_1
            fln_cpu(i,j,k,2) = fln_cpu(i,j,k,2) -    bulk_2
            fln_cpu(i,j,k,5) = fln_cpu(i,j,k,5) - uu*bulk_2
           endif
          enddo
         enddo
        enddo
    endsubroutine force_rhs_2_cpu

    subroutine force_rhs_1_cpu(nx, ny, nz, ng, yn_cpu, fln_cpu, w_cpu, w_aux_cpu, bulk5, fluid_mask_cpu)
        integer, intent(in) :: nx, ny, nz, ng
        real(rkind), dimension(1:), intent(in) :: yn_cpu
        real(rkind), dimension(1:,1:,1:,1:), intent(in) :: fln_cpu
        real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in) :: w_cpu
        real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in) :: w_aux_cpu
        integer, dimension(1-ng:,1-ng:,1-ng:), intent(in) :: fluid_mask_cpu
        real(rkind), dimension(5), intent(out) :: bulk5
        real(rkind) :: bulk_1, bulk_2, bulk_3, bulk_4, bulk_5
        real(rkind) :: dy
        integer :: i,j,k
        bulk_1 = 0._rkind
        bulk_2 = 0._rkind
        bulk_3 = 0._rkind
        bulk_4 = 0._rkind
        bulk_5 = 0._rkind
        
        do k=1,nz
         do j=1,ny
          do i=1,nx
           if (fluid_mask_cpu(i,j,k)==0) then
            dy = yn_cpu(j+1)-yn_cpu(j)
            bulk_1 = bulk_1 + fln_cpu(i,j,k,1)*dy
            bulk_2 = bulk_2 + fln_cpu(i,j,k,2)*dy
            bulk_3 = bulk_3 + w_cpu  (i,j,k,1)*dy
            bulk_4 = bulk_4 + w_cpu  (i,j,k,2)*dy
            bulk_5 = bulk_5 + w_cpu  (i,j,k,2)*dy*w_aux_cpu(i,j,k,6)
           endif
          enddo
         enddo
        enddo
       
        bulk5(1) = bulk_1
        bulk5(2) = bulk_2
        bulk5(3) = bulk_3
        bulk5(4) = bulk_4
        bulk5(5) = bulk_5
    endsubroutine force_rhs_1_cpu

    subroutine force_var_1_cpu(nx, ny, nz, ng, yn_cpu, fln_cpu, w_cpu, w_aux_cpu, bulkt, fluid_mask_cpu, cv_coeff_cpu, &
                               indx_cp_l, indx_cp_r, cv0, calorically_perfect, tol_iter_nr)
        integer, intent(in) :: nx, ny, nz, ng, indx_cp_l, indx_cp_r, calorically_perfect
        real(rkind), intent(in) :: cv0, tol_iter_nr
        real(rkind), dimension(1:), intent(in) :: yn_cpu
        real(rkind), dimension(1:,1:,1:,1:), intent(in) :: fln_cpu
        real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in) :: w_cpu
        real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(inout) :: w_aux_cpu
        integer, dimension(1-ng:,1-ng:,1-ng:), intent(in) :: fluid_mask_cpu
        real(rkind),    dimension(indx_cp_l:), intent(in) :: cv_coeff_cpu
        real(rkind), intent(out) :: bulkt
        real(rkind) :: bulk_5,rho,rhou,rhov,rhow,rhoe,ri,uu,vv,ww,qq,ee,tt
        real(rkind) :: dy
        integer :: i,j,k
!
        bulk_5 = 0._rkind
!        
        do k=1,nz
         do j=1,ny
          do i=1,nx
           if (fluid_mask_cpu(i,j,k)==0) then
            dy = yn_cpu(j+1)-yn_cpu(j)
            rho  = w_cpu(i,j,k,1)
            rhou = w_cpu(i,j,k,2)
            rhov = w_cpu(i,j,k,3)
            rhow = w_cpu(i,j,k,4)
            rhoe = w_cpu(i,j,k,5)
            ri   = 1._rkind/rho
            uu   = rhou*ri
            vv   = rhov*ri
            ww   = rhow*ri
            qq   = 0.5_rkind*(uu*uu+vv*vv+ww*ww)
            ee = rhoe/rho-qq
            tt = get_temperature_from_e_dev(ee,w_aux_cpu(i,j,k,6),cv0,cv_coeff_cpu,indx_cp_l,indx_cp_r, &
                                            calorically_perfect,tol_iter_nr)
            w_aux_cpu(i,j,k,6) = tt
            bulk_5 = bulk_5 + rhou*tt*dy
           endif
          enddo
         enddo
        enddo
!       
        bulkt = bulk_5
!
    endsubroutine force_var_1_cpu
!
    subroutine force_var_2_cpu(nx, ny, nz, ng, w_cpu, w_aux_cpu, tbdiff, fluid_mask_cpu, cv_coeff_cpu, indx_cp_l, indx_cp_r, &
                               cv0, calorically_perfect)
        integer :: nx, ny, nz, ng, indx_cp_l, indx_cp_r, calorically_perfect
        real(rkind), intent(in) :: tbdiff, cv0
        real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(inout) :: w_cpu
        real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in) :: w_aux_cpu
        integer, dimension(1-ng:,1-ng:,1-ng:), intent(in) :: fluid_mask_cpu
        real(rkind),    dimension(indx_cp_l:), intent(in) :: cv_coeff_cpu
        real(rkind) :: rho,rhou,rhov,rhow,tt,ttnew,ee
        integer :: i,j,k
!
        do k=1,nz
         do j=1,ny
          do i=1,nx
           if (fluid_mask_cpu(i,j,k)==0) then
            rho  = w_cpu(i,j,k,1)
            rhou = w_cpu(i,j,k,2)
            rhov = w_cpu(i,j,k,3)
            rhow = w_cpu(i,j,k,4)
            tt   = w_aux_cpu(i,j,k,6)
            ttnew = tt+tbdiff
            ee = get_e_from_temperature_dev(ttnew, cv0, indx_cp_l, indx_cp_r, cv_coeff_cpu, calorically_perfect)
            w_cpu(i,j,k,5) = rho*ee + 0.5_rkind*(rhou**2+rhov**2+rhow**2)/rho
           endif
          enddo
         enddo
        enddo
!
    endsubroutine force_var_2_cpu

    subroutine update_flux_cpu(nx, ny, nz, nv, fl_cpu, fln_cpu, gamdt) 
        integer :: nx, ny, nz, nv
        real(rkind) :: gamdt
        real(rkind), dimension(1:,1:,1:,1:), intent(in) :: fl_cpu
        real(rkind), dimension(1:,1:,1:,1:), intent(inout) :: fln_cpu
        integer :: i,j,k,m
         do k=1,nz
          do j=1,ny
           do i=1,nx
            do m=1,nv
             fln_cpu(i,j,k,m) = fln_cpu(i,j,k,m)-gamdt*fl_cpu(i,j,k,m)
            enddo
           enddo
          enddo
         enddo
    endsubroutine update_flux_cpu

    subroutine update_field_cpu(nx, ny, nz, ng, nv, w_cpu, fln_cpu, fluid_mask_cpu)
        integer :: nx, ny, nz, nv, ng
        real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(inout) :: w_cpu
        real(rkind), dimension(1:,1:,1:,1:), intent(in) :: fln_cpu
        integer, dimension(1-ng:,1-ng:,1-ng:), intent(in) :: fluid_mask_cpu
        integer :: i,j,k,m
         do k=1,nz
          do j=1,ny
           do i=1,nx
            if (fluid_mask_cpu(i,j,k)==0) then
             do m=1,nv
              w_cpu(i,j,k,m) = w_cpu(i,j,k,m)+fln_cpu(i,j,k,m)
             enddo
            endif
           enddo
          enddo
         enddo
    endsubroutine update_field_cpu

    subroutine visflx_div_cpu(nx, ny, nz, ng, visc_order, &
            w_aux_cpu, fl_cpu, coeff_deriv1_cpu, &
            dcsidx_cpu, detady_cpu, dzitdz_cpu)

        integer, intent(in) :: nx, ny, nz, ng, visc_order
        real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(inout) :: w_aux_cpu  
        real(rkind), dimension(1:,1:,1:, 1:), intent(inout) :: fl_cpu  
        real(rkind), dimension(1:,1:), intent(in) :: coeff_deriv1_cpu
        real(rkind), dimension(1:), intent(in) :: dcsidx_cpu, detady_cpu, dzitdz_cpu 
        integer     :: i,j,k,l
        real(rkind) :: ccl
        real(rkind) :: uu,vv,ww,tt,mu
        real(rkind) :: sigq,sigx,sigy,sigz
        real(rkind) :: divx3l, divy3l, divz3l
        integer     :: lmax

        lmax = visc_order/2
     
        do k=1,nz
         do j=1,ny
          do i=1,nx
      
           uu = w_aux_cpu(i,j,k,2)
           vv = w_aux_cpu(i,j,k,3)
           ww = w_aux_cpu(i,j,k,4)
           tt = w_aux_cpu(i,j,k,6)
           mu = w_aux_cpu(i,j,k,7)

           divx3l = 0._rkind
           divy3l = 0._rkind
           divz3l = 0._rkind

           do l=1,lmax
            ccl = coeff_deriv1_cpu(l,lmax)
            divx3l = divx3l+ccl*(w_aux_cpu(i+l,j,k,10)-w_aux_cpu(i-l,j,k,10))
            divy3l = divy3l+ccl*(w_aux_cpu(i,j+l,k,10)-w_aux_cpu(i,j-l,k,10))
            divz3l = divz3l+ccl*(w_aux_cpu(i,j,k+l,10)-w_aux_cpu(i,j,k-l,10))
           enddo

           divx3l = divx3l*dcsidx_cpu(i)
           divy3l = divy3l*detady_cpu(j)
           divz3l = divz3l*dzitdz_cpu(k)
       
           sigx  = mu*divx3l
           sigy  = mu*divy3l
           sigz  = mu*divz3l
           sigq  = sigx*uu+sigy*vv+sigz*ww
      
           fl_cpu(i,j,k,2) = fl_cpu(i,j,k,2) - sigx
           fl_cpu(i,j,k,3) = fl_cpu(i,j,k,3) - sigy
           fl_cpu(i,j,k,4) = fl_cpu(i,j,k,4) - sigz
           fl_cpu(i,j,k,5) = fl_cpu(i,j,k,5) - sigq 
     
          enddo
         enddo
        enddo

    endsubroutine visflx_div_cpu

    subroutine visflx_cpu(nx, ny, nz, ng, visc_order, &
            Prandtl, cp0, indx_cp_l, indx_cp_r, cp_coeff_cpu, calorically_perfect, &
            u0, l0, w_cpu, w_aux_cpu, fl_cpu, &
            coeff_deriv1_cpu, coeff_deriv2_cpu, &
            dcsidx_cpu, detady_cpu, dzitdz_cpu,  &
            dcsidxs_cpu, detadys_cpu, dzitdzs_cpu, & 
            dcsidx2_cpu, detady2_cpu, dzitdz2_cpu, wallprop_cpu)

        integer, intent(in) :: nx, ny, nz, ng, visc_order, calorically_perfect
        integer, intent(in) :: indx_cp_l, indx_cp_r
        real(rkind), intent(in) :: Prandtl, u0, l0, cp0
        real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(in) :: w_cpu 
        real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(inout) :: w_aux_cpu 
        real(rkind), dimension(1-ng:,1-ng:, 2:), intent(inout) :: wallprop_cpu 
        real(rkind), dimension(1:,1:,1:, 1:), intent(inout) :: fl_cpu
        real(rkind), dimension(1:,1:), intent(in) :: coeff_deriv1_cpu
        real(rkind), dimension(0:,1:), intent(in)  :: coeff_deriv2_cpu
        real(rkind), dimension(1:), intent(in) :: dcsidx_cpu, detady_cpu, dzitdz_cpu 
        real(rkind), dimension(1:), intent(in) :: dcsidxs_cpu, detadys_cpu, dzitdzs_cpu
        real(rkind), dimension(1:), intent(in) :: dcsidx2_cpu, detady2_cpu, dzitdz2_cpu
        real(rkind), dimension(indx_cp_l:), intent(in) :: cp_coeff_cpu
        integer     :: i,j,k,l,ll
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

        lmax = visc_order/2
     
        do j=1,ny
         do i=1,nx
          do k=1,nz
      
           uu = w_aux_cpu(i,j,k,2)
           vv = w_aux_cpu(i,j,k,3)
           ww = w_aux_cpu(i,j,k,4)
           tt = w_aux_cpu(i,j,k,6)
           mu = w_aux_cpu(i,j,k,7)

           if (calorically_perfect==1) then
            cploc = cp0
           else
            cploc = 0._rkind
            do ll=indx_cp_l,indx_cp_r
             cploc = cploc+cp_coeff_cpu(ll)*tt**ll
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
            ccl = coeff_deriv1_cpu(l,lmax)
            ux = ux+ccl*(w_aux_cpu(i+l,j,k,2)-w_aux_cpu(i-l,j,k,2))
            vx = vx+ccl*(w_aux_cpu(i+l,j,k,3)-w_aux_cpu(i-l,j,k,3))
            wx = wx+ccl*(w_aux_cpu(i+l,j,k,4)-w_aux_cpu(i-l,j,k,4))
            tx = tx+ccl*(w_aux_cpu(i+l,j,k,6)-w_aux_cpu(i-l,j,k,6))
            mux = mux+ccl*(w_aux_cpu(i+l,j,k,7)-w_aux_cpu(i-l,j,k,7))
      
            uy = uy+ccl*(w_aux_cpu(i,j+l,k,2)-w_aux_cpu(i,j-l,k,2))
            vy = vy+ccl*(w_aux_cpu(i,j+l,k,3)-w_aux_cpu(i,j-l,k,3))
            wy = wy+ccl*(w_aux_cpu(i,j+l,k,4)-w_aux_cpu(i,j-l,k,4))
            ty = ty+ccl*(w_aux_cpu(i,j+l,k,6)-w_aux_cpu(i,j-l,k,6))
            muy = muy+ccl*(w_aux_cpu(i,j+l,k,7)-w_aux_cpu(i,j-l,k,7))
     
            uz = uz+ccl*(w_aux_cpu(i,j,k+l,2)-w_aux_cpu(i,j,k-l,2))
            vz = vz+ccl*(w_aux_cpu(i,j,k+l,3)-w_aux_cpu(i,j,k-l,3))
            wz = wz+ccl*(w_aux_cpu(i,j,k+l,4)-w_aux_cpu(i,j,k-l,4))
            tz = tz+ccl*(w_aux_cpu(i,j,k+l,6)-w_aux_cpu(i,j,k-l,6))
            muz = muz+ccl*(w_aux_cpu(i,j,k+l,7)-w_aux_cpu(i,j,k-l,7))
           enddo
           ux = ux*dcsidx_cpu(i)
           vx = vx*dcsidx_cpu(i)
           wx = wx*dcsidx_cpu(i)
           tx = tx*dcsidx_cpu(i)
           mux = mux*dcsidx_cpu(i)
           uy = uy*detady_cpu(j)
           vy = vy*detady_cpu(j)
           wy = wy*detady_cpu(j)
           ty = ty*detady_cpu(j)
           muy = muy*detady_cpu(j)
           uz = uz*dzitdz_cpu(k)
           vz = vz*dzitdz_cpu(k)
           wz = wz*dzitdz_cpu(k)
           tz = tz*dzitdz_cpu(k)
           muz = muz*dzitdz_cpu(k)
!
           if (j==1) then
            wallprop_cpu(i,k,2) = mu*uy
            wallprop_cpu(i,k,3) = mu*wy
            wallprop_cpu(i,k,4) = mu*ty*cploc/Prandtl
           endif
!      
           ulapx = coeff_deriv2_cpu(0,lmax)*uu
           ulapy = ulapx
           ulapz = ulapx
           vlapx = coeff_deriv2_cpu(0,lmax)*vv
           vlapy = vlapx
           vlapz = vlapx
           wlapx = coeff_deriv2_cpu(0,lmax)*ww
           wlapy = wlapx
           wlapz = wlapx
           tlapx = coeff_deriv2_cpu(0,lmax)*tt
           tlapy = tlapx
           tlapz = tlapx
           do l=1,lmax
            clapl = coeff_deriv2_cpu(l,lmax)
            ulapx = ulapx + clapl*(w_aux_cpu(i+l,j,k,2)+w_aux_cpu(i-l,j,k,2))
            ulapy = ulapy + clapl*(w_aux_cpu(i,j+l,k,2)+w_aux_cpu(i,j-l,k,2))
            ulapz = ulapz + clapl*(w_aux_cpu(i,j,k+l,2)+w_aux_cpu(i,j,k-l,2))
            vlapx = vlapx + clapl*(w_aux_cpu(i+l,j,k,3)+w_aux_cpu(i-l,j,k,3))
            vlapy = vlapy + clapl*(w_aux_cpu(i,j+l,k,3)+w_aux_cpu(i,j-l,k,3))
            vlapz = vlapz + clapl*(w_aux_cpu(i,j,k+l,3)+w_aux_cpu(i,j,k-l,3))
            wlapx = wlapx + clapl*(w_aux_cpu(i+l,j,k,4)+w_aux_cpu(i-l,j,k,4))
            wlapy = wlapy + clapl*(w_aux_cpu(i,j+l,k,4)+w_aux_cpu(i,j-l,k,4))
            wlapz = wlapz + clapl*(w_aux_cpu(i,j,k+l,4)+w_aux_cpu(i,j,k-l,4))
            tlapx = tlapx + clapl*(w_aux_cpu(i+l,j,k,6)+w_aux_cpu(i-l,j,k,6))
            tlapy = tlapy + clapl*(w_aux_cpu(i,j+l,k,6)+w_aux_cpu(i,j-l,k,6))
            tlapz = tlapz + clapl*(w_aux_cpu(i,j,k+l,6)+w_aux_cpu(i,j,k-l,6))
           enddo
      
           ulapx = ulapx*dcsidxs_cpu(i)+ux*dcsidx2_cpu(i)
           vlapx = vlapx*dcsidxs_cpu(i)+vx*dcsidx2_cpu(i)
           wlapx = wlapx*dcsidxs_cpu(i)+wx*dcsidx2_cpu(i)
           tlapx = tlapx*dcsidxs_cpu(i)+tx*dcsidx2_cpu(i)
           ulapy = ulapy*detadys_cpu(j)+uy*detady2_cpu(j)
           vlapy = vlapy*detadys_cpu(j)+vy*detady2_cpu(j)
           wlapy = wlapy*detadys_cpu(j)+wy*detady2_cpu(j)
           tlapy = tlapy*detadys_cpu(j)+ty*detady2_cpu(j)
           ulapz = ulapz*dzitdzs_cpu(k)+uz*dzitdz2_cpu(k)
           vlapz = vlapz*dzitdzs_cpu(k)+vz*dzitdz2_cpu(k)
           wlapz = wlapz*dzitdzs_cpu(k)+wz*dzitdz2_cpu(k)
           tlapz = tlapz*dzitdzs_cpu(k)+tz*dzitdz2_cpu(k)
      
           ulap  = ulapx+ulapy+ulapz
           vlap  = vlapx+vlapy+vlapz
           wlap  = wlapx+wlapy+wlapz
           tlap  = tlapx+tlapy+tlapz
     
           div     = ux+vy+wz
           div3l   = div/3._rkind
           w_aux_cpu(i,j,k,10) = div3l
           omegax = wy-vz
           omegay = uz-wx
           omegaz = vx-uy
           omod2 = omegax*omegax+omegay*omegay+omegaz*omegaz
           w_aux_cpu(i,j,k,9) = sqrt(omod2)

           ! Original Ducros shock sensor
           ! w_aux_cpu(i,j,k,8) = div2/(omod2+div2+1.D-12)
!
           ! Modified Ducros shock sensor
           ! w_aux_cpu(i,j,k,8) = div2/(omod2+div2+(u0/l0)**2)
!
           ! Modified Ducros shock sensor (remove expansions)
            w_aux_cpu(i,j,k,8) = (max(-div/sqrt(omod2+div**2+(u0/l0)**2),0._rkind))**2
!
           ! ducfilter = div2/(div2+0.1_rkind*sqrt(uu*uu+vv*vv+ww*ww)*max(dcsidx_cpu(i),detady_cpu(j),dzitdz_cpu(k))+1.D-12)
           ! w_aux_cpu(i,j,k,8) = (max(-div/sqrt(omod2+div**2+1.D-12),0._rkind))**2*ducfilter

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
           
           sigah = (sig11*ux+sig12*uy+sig13*uz+sig12*vx+sig22*vy+sig23*vz+sig13*wx+sig23*wy+sig33*wz)*mu ! Aerodynamic heating
       
           sigq  = sigx*uu+sigy*vv+sigz*ww+sigah+sigqt
      
           fl_cpu(i,j,k,2) = fl_cpu(i,j,k,2) - sigx
           fl_cpu(i,j,k,3) = fl_cpu(i,j,k,3) - sigy
           fl_cpu(i,j,k,4) = fl_cpu(i,j,k,4) - sigz
           fl_cpu(i,j,k,5) = fl_cpu(i,j,k,5) - sigq 
     
          enddo
         enddo
        enddo
    endsubroutine visflx_cpu

    subroutine recyc_exchange_cpu_1(irecyc, w_cpu, wbuf1s_cpu, nx, ny, nz, ng, nv)
      integer, intent(in) :: irecyc, nx, ny, nz, ng, nv
      real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in) :: w_cpu
      real(rkind), dimension(:,:,:,:), intent(inout) :: wbuf1s_cpu
      integer :: i,j,k,m
      do k=1,nz
       do j=1,ny
        do i=1,ng
         do m=1,nv
          wbuf1s_cpu(i,j,k,m) = w_cpu(irecyc+1-i,j,k,m)
         enddo
        enddo
       enddo
      enddo
    endsubroutine recyc_exchange_cpu_1
    subroutine recyc_exchange_cpu_2(n1_start_recv, n1_start_send, n1_end_recv, wrecyc_cpu, wbuf1r_cpu, nx, ny, nz, ng, nv)
        integer, intent(in) :: n1_start_recv, n1_start_send, n1_end_recv, nx, ny, nz, ng, nv
        real(rkind), dimension(:,:,:,:), intent(in) :: wbuf1r_cpu
        real(rkind), dimension(:,:,:,:), intent(inout) :: wrecyc_cpu
        integer :: i,j,k,m
         do k=n1_start_recv,n1_end_recv
          do j=1,ny
           do i=1,ng
            do m=1,nv
             wrecyc_cpu(i,j,k,m) = wbuf1r_cpu(i,j,k-n1_start_recv+n1_start_send,m)
            enddo
           enddo
          enddo
         enddo
    endsubroutine recyc_exchange_cpu_2
    subroutine recyc_exchange_cpu_3(n2_start_recv, n2_start_send, n2_end_recv, wrecyc_cpu, wbuf2r_cpu, nx, ny, nz, ng, nv)
        integer, intent(in) :: n2_start_recv, n2_start_send, n2_end_recv, nx, ny, nz, ng, nv
        real(rkind), dimension(:,:,:,:), intent(in) :: wbuf2r_cpu
        real(rkind), dimension(:,:,:,:), intent(inout) :: wrecyc_cpu
        integer :: i,j,k,m
         do k=n2_start_recv,n2_end_recv
          do j=1,ny
           do i=1,ng
            do m=1,nv
             wrecyc_cpu(i,j,k,m) = wbuf2r_cpu(i,j,k-n2_start_recv+n2_start_send,m)
            enddo
           enddo
          enddo
         enddo
    endsubroutine recyc_exchange_cpu_3

     subroutine bc_nr_lat_x_kernel(start_or_end, nr_type, &
                                  nx, ny, nz, ng, nv, w_aux_cpu, w_cpu, fl_cpu, &
                                  dcsidx_cpu, indx_cp_l, indx_cp_r, cp_coeff_cpu, winf_cpu, calorically_perfect)
        integer, intent(in) :: start_or_end, nr_type, nx, ny, nz, ng, nv, indx_cp_l, indx_cp_r, calorically_perfect
        real(rkind), dimension(3) :: c_one
        real(rkind), dimension(5) :: dw_dn, dwc_dn, ev, dw_dn_outer, dwc_dn_outer
        real(rkind), dimension(5,5) :: el, er
        real(rkind) :: w_target
        integer :: i, j, k, l, m, mm, sgn_dw
        real(rkind) :: df, uu, vv, ww, h, qq, cc, c, ci, b2, b1
        real(rkind), dimension(indx_cp_l:) :: cp_coeff_cpu
        real(rkind), dimension(:,:,:,:) :: fl_cpu
        real(rkind), dimension(1-ng:,1-ng:,1-ng:,:) :: w_aux_cpu, w_cpu
        real(rkind), dimension(:) :: dcsidx_cpu
        real(rkind), dimension(:) :: winf_cpu
        real(rkind) :: rho,tt,gamloc,p_rho,p_e,etot,b3



        do j=1,ny
do k=1,nz


        ! Setup min or max boundary
        c_one = [-1.5_rkind, 2._rkind, -0.5_rkind]
        if(start_or_end == 1) then
            i      = 1
            sgn_dw = 1
        elseif(start_or_end == 2) then
            i      = nx
            sgn_dw = -1
        endif

        ! Compute d(U_cons)/dx: inner, and outer for relaxation
        do m=1,nv
            dw_dn(m) = 0._rkind
            do l=1,3
                dw_dn(m) = dw_dn(m) + sgn_dw * c_one(l)*w_cpu(i+sgn_dw*(l-1),j,k,m)
            enddo
            ! Relax to w_cpu which if imposed by recycing works. Problems
            ! could arise at the exit, but we do not relax there.
            ! Another possible choice would be relaxing to w_inf.
            !w_target       = winf_cpu(m)
            w_target       = w_cpu(i-sgn_dw,j,k,m)
            dw_dn_outer(m) = sgn_dw * (w_cpu(i,j,k,m)-w_target)
        enddo

        ! Compute eigenvectors
        rho       = w_aux_cpu(i,j,k,1)
        uu        = w_aux_cpu(i,j,k,2)
        vv        = w_aux_cpu(i,j,k,3)
        ww        = w_aux_cpu(i,j,k,4)
        h         = w_aux_cpu(i,j,k,5)
        tt        = w_aux_cpu(i,j,k,6)
        qq        = 0.5_rkind * (uu*uu  +vv*vv + ww*ww)
!
        gamloc    = get_gamloc_dev(indx_cp_l,indx_cp_r,tt,cp_coeff_cpu,calorically_perfect)
        cc        = gamloc * tt
        c         =  sqrt(cc) 
        ci        =  1._rkind/c 
!
        p_rho     = tt
        p_e       = rho*(gamloc-1._rkind)   ! p_e   = rho * R_gas / Cv
        etot      = h - tt
!
        b3        = etot - rho * p_rho/p_e ! qq in the case of calorically perfect
        b2        = p_e/(rho*cc)
        b1        = p_rho/cc - b2*(etot - 2._rkind*qq)
!
!       b2        =  gm1/cc  ! 1/(cp*T)
!       b1        =  b2 * qq
!
        call eigenvectors_x(b1, b2, b3, uu, vv, ww, c, ci, h, el, er)

        ! Pre-multiply to L to get derivative of characteristic variables
        do m=1,5
            dwc_dn(m)       = 0._rkind
            dwc_dn_outer(m) = 0._rkind
            do mm=1,5
                dwc_dn(m)       = dwc_dn(m)       + el(mm,m) * dw_dn(mm)
                dwc_dn_outer(m) = dwc_dn_outer(m) + el(mm,m) * dw_dn_outer(mm)
            enddo
        enddo

        ! Compute eigenvalues
        ev(1) = uu-c 
        ev(2) = uu   
        ev(3) = uu+c 
        ev(4) = ev(2)
        ev(5) = ev(2)

        ! If nr_type=1, kill acousting ingoing waves
        if(nr_type == 1) then
            do l=1,5
                ev(l) = sgn_dw*min(sgn_dw*ev(l) ,0._rkind)
            enddo
        endif

        ! If nr_type=2, exiting waves are kept, incoming waves are assigned from outer derivatives
        if(nr_type == 2) then
            do m=1,5
                if(sgn_dw*ev(m) > 0._rkind) then
                    dwc_dn(m) = dwc_dn_outer(m)
                endif
            enddo
        endif

        ! Compute wave amplitude vector (lambda*L*dw)
        do m=1,5
            dwc_dn(m) = ev(m) * dwc_dn(m)
        enddo

        ! If nr_type=6, enforce wave reflection to simulate purely reflective wall (LODI relations)
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

        ! Pre-multiply to R to return to conservative variables and assign result to fl
        do m=1,5
            df = 0._rkind
            do mm=1,5
                df = df + er(mm,m) * dwc_dn(mm)
            enddo
            fl_cpu(i,j,k,m) = fl_cpu(i,j,k,m) + df * dcsidx_cpu(i)
        enddo
    enddo
enddo
endsubroutine bc_nr_lat_x_kernel


     subroutine bc_nr_lat_y_kernel(start_or_end, nr_type, &
                                  nx, ny, nz, ng, nv, w_aux_cpu, w_cpu, fl_cpu, detady_cpu, &
                                  indx_cp_l, indx_cp_r, cp_coeff_cpu, calorically_perfect)
        integer, intent(in) :: start_or_end, nr_type, nx, ny, nz, ng, nv, indx_cp_l, indx_cp_r, calorically_perfect
        real(rkind), dimension(3) :: c_one
        real(rkind), dimension(5) :: dw_dn, dwc_dn, ev
        real(rkind), dimension(5,5) :: el, er
        integer :: i, j, k, l, m, mm, sgn_dw
        real(rkind) :: df, uu, vv, ww, h, qq, cc, c, ci, b2, b1
        real(rkind), dimension(indx_cp_l:) :: cp_coeff_cpu
        real(rkind), dimension(:,:,:,:) :: fl_cpu
        real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: w_aux_cpu, w_cpu
        real(rkind), dimension(:) :: detady_cpu
        real(rkind) :: rho,tt,gamloc,p_rho,p_e,etot,b3



        do i=1,nx
do k=1,nz


        ! Setup min or max boundary
        c_one = [-1.5_rkind, 2._rkind, -0.5_rkind]
        if(start_or_end == 1) then
            j      = 1
            sgn_dw = 1
        elseif(start_or_end == 2) then
            j      = ny
            sgn_dw = -1
        endif

        ! Compute d(U_cons)/dx
        do m=1,nv
            dw_dn(m) = 0._rkind
            do l=1,3
                dw_dn(m) = dw_dn(m) + sgn_dw * c_one(l)*w_cpu(i,j+sgn_dw*(l-1),k,m)
            enddo
        enddo

        ! Compute eigenvectors
        rho       =  w_aux_cpu(i,j,k,1)
        uu        =  w_aux_cpu(i,j,k,2)
        vv        =  w_aux_cpu(i,j,k,3)
        ww        =  w_aux_cpu(i,j,k,4)
        h         =  w_aux_cpu(i,j,k,5)
        tt        =  w_aux_cpu(i,j,k,6)
        qq        =  0.5_rkind * (uu*uu  +vv*vv + ww*ww)
!
        gamloc    = get_gamloc_dev(indx_cp_l,indx_cp_r,tt,cp_coeff_cpu,calorically_perfect)
        cc        = gamloc * tt
        c         =  sqrt(cc) 
        ci        =  1._rkind/c 
!
        p_rho     = tt
        p_e       = rho*(gamloc-1._rkind)   ! p_e   = rho * R_gas / Cv
        etot      = h - tt
!
        b3        = etot - rho * p_rho/p_e ! qq in the case of calorically perfect
        b2        = p_e/(rho*cc)
        b1        = p_rho/cc - b2*(etot - 2._rkind*qq)
!
!       b2        =  gm1/cc  ! 1/(cp*T)
!       b1        =  b2 * qq
!
        call eigenvectors_y(b1, b2, b3, uu, vv, ww, c, ci, h, el, er)

        ! Pre-multiply to L to get derivative of characteristic variables
        do m=1,5
            dwc_dn(m) = 0._rkind
            do mm=1,5
                dwc_dn(m) = dwc_dn(m) + el(mm,m) * dw_dn(mm)
            enddo
        enddo

        ! Compute eigenvalues
        ev(1) = vv-c 
        ev(2) = vv   
        ev(3) = vv+c 
        ev(4) = ev(2)
        ev(5) = ev(2)

        ! If nr_type=1, kill acousting ingoing waves
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

        ! Compute wave amplitude vector (lambda*L*dw)
        do m=1,5
            dwc_dn(m) = ev(m) * dwc_dn(m)
        enddo

        ! If nr_type=6, enforce wave reflection to simulate purely reflective wall (LODI relations)
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

        ! Pre-multiply to R to return to conservative variables and assign result to fl
        do m=1,5
            df = 0._rkind
            do mm=1,5
                df = df + er(mm,m) * dwc_dn(mm)
            enddo
            fl_cpu(i,j,k,m) = fl_cpu(i,j,k,m) + df * detady_cpu(j)
        enddo
    enddo
enddo
endsubroutine bc_nr_lat_y_kernel


     subroutine bc_nr_lat_z_kernel(start_or_end, nr_type, &
                                  nx, ny, nz, ng, nv, w_aux_cpu, w_cpu, fl_cpu, dzitdz_cpu, &
                                  indx_cp_l, indx_cp_r, cp_coeff_cpu, calorically_perfect)
        integer, intent(in) :: start_or_end, nr_type, nx, ny, nz, ng, nv, indx_cp_l, indx_cp_r, calorically_perfect
        real(rkind), dimension(3) :: c_one
        real(rkind), dimension(5) :: dw_dn, dwc_dn, ev
        real(rkind), dimension(5,5) :: el, er
        integer :: i, j, k, l, m, mm, sgn_dw
        real(rkind) :: df, uu, vv, ww, h, qq, cc, c, ci, b2, b1
        real(rkind), dimension(indx_cp_l:) :: cp_coeff_cpu
        real(rkind), dimension(:,:,:,:) :: fl_cpu
        real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: w_aux_cpu, w_cpu
        real(rkind), dimension(:) :: dzitdz_cpu
        real(rkind) :: rho,tt,gamloc,p_rho,p_e,etot,b3



        do i=1,nx
do j=1,ny


        ! Setup min or max boundary
        c_one = [-1.5_rkind, 2._rkind, -0.5_rkind]
        if(start_or_end == 1) then
            k      = 1
            sgn_dw = 1
        elseif(start_or_end == 2) then
            k      = nz
            sgn_dw = -1
        endif

        ! Compute d(U_cons)/dx
        do m=1,nv
            dw_dn(m) = 0._rkind
            do l=1,3
                dw_dn(m) = dw_dn(m) + sgn_dw * c_one(l)*w_cpu(i,j,k+sgn_dw*(l-1),m)
            enddo
        enddo

        ! Compute eigenvectors
        rho       =  w_aux_cpu(i,j,k,1)
        uu        =  w_aux_cpu(i,j,k,2)
        vv        =  w_aux_cpu(i,j,k,3)
        ww        =  w_aux_cpu(i,j,k,4)
        h         =  w_aux_cpu(i,j,k,5)
        tt        =  w_aux_cpu(i,j,k,6)
        qq        =  0.5_rkind * (uu*uu  +vv*vv + ww*ww)
!
        gamloc    = get_gamloc_dev(indx_cp_l,indx_cp_r,tt,cp_coeff_cpu,calorically_perfect)
        cc        = gamloc * tt
        c         =  sqrt(cc) 
        ci        =  1._rkind/c 
!
        p_rho     = tt
        p_e       = rho*(gamloc-1._rkind)   ! p_e   = rho * R_gas / Cv
        etot      = h - tt
!
        b3        = etot - rho * p_rho/p_e ! qq in the case of calorically perfect
        b2        = p_e/(rho*cc)
        b1        = p_rho/cc - b2*(etot - 2._rkind*qq)
!
!       b2        =  gm1/cc  ! 1/(cp*T)
!       b1        =  b2 * qq
!
        call eigenvectors_z(b1, b2, b3, uu, vv, ww, c, ci, h, el, er)

        ! Pre-multiply to L to get derivative of characteristic variables
        do m=1,5
            dwc_dn(m) = 0._rkind
            do mm=1,5
                dwc_dn(m) = dwc_dn(m) + el(mm,m) * dw_dn(mm)
            enddo
        enddo

        ! Compute eigenvalues
        ev(1) = ww-c 
        ev(2) = ww   
        ev(3) = ww+c 
        ev(4) = ev(2)
        ev(5) = ev(2)

        ! If nr_type=1, kill acousting ingoing waves
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

        ! Compute wave amplitude vector (lambda*L*dw)
        do m=1,5
            dwc_dn(m) = ev(m) * dwc_dn(m)
        enddo

        ! If nr_type=6, enforce wave reflection to simulate purely reflective wall (LODI relations)
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

        ! Pre-multiply to R to return to conservative variables and assign result to fl
        do m=1,5
            df = 0._rkind
            do mm=1,5
                df = df + er(mm,m) * dwc_dn(mm)
            enddo
            fl_cpu(i,j,k,m) = fl_cpu(i,j,k,m) + df * dzitdz_cpu(k)
        enddo
    enddo
enddo
endsubroutine bc_nr_lat_z_kernel


    subroutine  bcrecyc_cpu_1(nx, ny, nz, ng, nv, wrecycav_cpu, wrecyc_cpu)
        integer, intent(in) :: nx, ny, nz, ng, nv
        real(rkind), dimension(:,:,:), intent(inout) :: wrecycav_cpu
        real(rkind), dimension(:,:,:,:), intent(in) :: wrecyc_cpu
        integer :: i,j,k,m
        do j=1,ny
         do i=1,ng
          do m=1,nv
           wrecycav_cpu(i,j,m) = 0._rkind
           do k=1,nz
            wrecycav_cpu(i,j,m) = wrecycav_cpu(i,j,m)+wrecyc_cpu(i,j,k,m)
           enddo
          enddo
         enddo
        enddo
    endsubroutine  bcrecyc_cpu_1

 subroutine  bcrecyc_cpu_2(nx, ny, nz, nzmax, ng, wrecycav_cpu, wrecyc_cpu)
     integer, intent(in) :: nx, ny, nz, nzmax, ng
     real(rkind), dimension(:,:,:), intent(in) :: wrecycav_cpu
     real(rkind), dimension(:,:,:,:), intent(inout) :: wrecyc_cpu
     real(rkind) :: ufav, vfav, wfav, rhom
     integer :: i,j,k

     do j=1,ny
      do i=1,ng
       ufav = wrecycav_cpu(i,j,2)/wrecycav_cpu(i,j,1)
       vfav = wrecycav_cpu(i,j,3)/wrecycav_cpu(i,j,1)
       wfav = wrecycav_cpu(i,j,4)/wrecycav_cpu(i,j,1)
       rhom = wrecycav_cpu(i,j,1)/nzmax
       do k=1,nz
        wrecyc_cpu(i,j,k,2) = wrecyc_cpu(i,j,k,2)/wrecyc_cpu(i,j,k,1)-ufav ! Velocity fluctuations
        wrecyc_cpu(i,j,k,3) = wrecyc_cpu(i,j,k,3)/wrecyc_cpu(i,j,k,1)-vfav
        wrecyc_cpu(i,j,k,4) = wrecyc_cpu(i,j,k,4)/wrecyc_cpu(i,j,k,1)-wfav
        wrecyc_cpu(i,j,k,1) = wrecyc_cpu(i,j,k,1)-rhom                     ! Density fluctuations
       enddo
      enddo
     enddo

 endsubroutine  bcrecyc_cpu_2

 subroutine bcrecyc_cpu_3(nx, ny, nz, ng, p0, w_cpu, wmean_cpu, wrecyc_cpu, &
     weta_inflow_cpu, map_j_inn_cpu, map_j_out_cpu, &
     yplus_inflow_cpu, eta_inflow_cpu, yplus_recyc_cpu, eta_recyc_cpu, betarecyc, &
     indx_cp_l, indx_cp_r, cv_coeff_cpu, cv0, calorically_perfect)
     integer, intent(in) :: nx, ny, nz, ng, indx_cp_l, indx_cp_r, calorically_perfect
     real(rkind) :: p0, betarecyc, cv0
     real(rkind), dimension(indx_cp_l:), intent(in) :: cv_coeff_cpu
     real(rkind), dimension(1-ng:,:,:), intent(in) :: wmean_cpu
     real(rkind), dimension(:,:,:,:), intent(in) :: wrecyc_cpu
     real(rkind), dimension(1-ng:,1-ng:,1-ng:,:), intent(inout) :: w_cpu
     real(rkind), dimension(:), intent(in) :: weta_inflow_cpu 
     real(rkind), dimension(1-ng:), intent(in) :: yplus_inflow_cpu, eta_inflow_cpu, yplus_recyc_cpu, eta_recyc_cpu
     integer, dimension(:), intent(in) :: map_j_inn_cpu, map_j_out_cpu
     integer :: i,j,k
     real(rkind) :: weta, weta1, bdamp, disty_inn, disty_out, rhofluc, ufluc, vfluc, wfluc, rhof_inn, rhof_out
     real(rkind) :: uf_inn, uf_out, vf_inn, vf_out, wf_inn, wf_out
     real(rkind) :: rhomean, uumean , vvmean , wwmean , tmean  , rho , uu , vv , ww , rhou , rhov , rhow, tt, ee
     integer :: j_inn, j_out

      do k=1,nz
       do j=1,ny
        weta  = weta_inflow_cpu(j)
        weta1 = 1._rkind-weta
        bdamp = 0.5_rkind*(1._rkind-tanh(4._rkind*(eta_inflow_cpu(j)-2._rkind))) 
        j_inn = map_j_inn_cpu(j) 
        j_out = map_j_out_cpu(j) 
        disty_inn = (yplus_inflow_cpu(j)-yplus_recyc_cpu(j_inn))/(yplus_recyc_cpu(j_inn+1)-yplus_recyc_cpu(j_inn))    
        disty_out = (eta_inflow_cpu(j)-eta_recyc_cpu(j_out))/(eta_recyc_cpu(j_out+1)-eta_recyc_cpu(j_out))    
     
        do i=1,ng
     
         if (j==1.or.j_inn>=ny.or.j_out>=ny) then
          rhofluc = 0._rkind
          ufluc   = 0._rkind
          vfluc   = 0._rkind
          wfluc   = 0._rkind
         else
          rhof_inn = wrecyc_cpu(i,j_inn,k,1)*(1._rkind-disty_inn)+wrecyc_cpu(i,j_inn+1,k,1)*disty_inn
          rhof_out = wrecyc_cpu(i,j_out,k,1)*(1._rkind-disty_out)+wrecyc_cpu(i,j_out+1,k,1)*disty_out
          uf_inn   = wrecyc_cpu(i,j_inn,k,2)*(1._rkind-disty_inn)+wrecyc_cpu(i,j_inn+1,k,2)*disty_inn
          uf_out   = wrecyc_cpu(i,j_out,k,2)*(1._rkind-disty_out)+wrecyc_cpu(i,j_out+1,k,2)*disty_out
          vf_inn   = wrecyc_cpu(i,j_inn,k,3)*(1._rkind-disty_inn)+wrecyc_cpu(i,j_inn+1,k,3)*disty_inn
          vf_out   = wrecyc_cpu(i,j_out,k,3)*(1._rkind-disty_out)+wrecyc_cpu(i,j_out+1,k,3)*disty_out
          wf_inn   = wrecyc_cpu(i,j_inn,k,4)*(1._rkind-disty_inn)+wrecyc_cpu(i,j_inn+1,k,4)*disty_inn
          wf_out   = wrecyc_cpu(i,j_out,k,4)*(1._rkind-disty_out)+wrecyc_cpu(i,j_out+1,k,4)*disty_out
     !
          rhofluc = rhof_inn*weta1+rhof_out*weta
          ufluc   =   uf_inn*weta1+  uf_out*weta
          vfluc   =   vf_inn*weta1+  vf_out*weta
          wfluc   =   wf_inn*weta1+  wf_out*weta
          rhofluc = rhofluc*bdamp
          ufluc   = ufluc  *bdamp*betarecyc 
          vfluc   = vfluc  *bdamp*betarecyc 
          wfluc   = wfluc  *bdamp*betarecyc 
         endif
     
         rhomean = wmean_cpu(1-i,j,1)
         uumean  = wmean_cpu(1-i,j,2)/rhomean
         vvmean  = wmean_cpu(1-i,j,3)/rhomean
         wwmean  = wmean_cpu(1-i,j,4)/rhomean
         tmean   = p0/rhomean
         rho     = rhomean + rhofluc
         uu      = uumean  + ufluc
         vv      = vvmean  + vfluc
         ww      = wwmean  + wfluc
         rhou    = rho*uu
         rhov    = rho*vv
         rhow    = rho*ww
         
         w_cpu(1-i,j,k,1) = rho
         w_cpu(1-i,j,k,2) = rhou 
         w_cpu(1-i,j,k,3) = rhov 
         w_cpu(1-i,j,k,4) = rhow 
         tt               = p0/rho
         ee = get_e_from_temperature_dev(tt, cv0, indx_cp_l, indx_cp_r, cv_coeff_cpu, calorically_perfect)
         w_cpu(1-i,j,k,5) = rho*ee + 0.5_rkind*(rhou**2+rhov**2+rhow**2)/rho
         !w_cpu(1-i,j,k,5) = p0*gm + 0.5_rkind*(rhou**2+rhov**2+rhow**2)/rho
        enddo
       enddo
      enddo
    endsubroutine bcrecyc_cpu_3

    subroutine bcfree_cpu(ilat, nx, ny, nz, ng, nv, rho0, u0, e0, w_cpu)
       integer :: nx,ny,nz,ng,nv
       integer :: ilat
       real(rkind) :: rho0,u0,e0
       integer :: j,k,l
       real(rkind), dimension(1-ng:, 1-ng:, 1-ng:, 1:), intent(inout) :: w_cpu

       if (ilat==1) then     ! left side
        do k=1,nz
         do j=1,ny
          do l=1,ng
           w_cpu(1-l,j,k,1) = rho0
           w_cpu(1-l,j,k,2) = rho0*u0
           w_cpu(1-l,j,k,3) = 0._rkind
           w_cpu(1-l,j,k,4) = 0._rkind
           w_cpu(1-l,j,k,5) = rho0*e0+0.5_rkind*rho0*u0**2
          enddo
         enddo
        enddo
       elseif (ilat==2) then ! right side
       elseif (ilat==3) then ! lower side
       elseif (ilat==4) then  ! upper side
       elseif (ilat==5) then  ! back side
       elseif (ilat==6) then  ! fore side
       endif
    endsubroutine bcfree_cpu

    subroutine bcshock_cpu(ilat, nx, ny, nz, ng, nv, w_cpu, winf_cpu, winf_past_shock_cpu, xshock_imp, shock_angle, &
                           x_cpu, y_cpu, tanhfacs)
       integer, intent(in) :: nx,ny,nz,ng,nv,ilat
       real(rkind), intent(in) :: xshock_imp, shock_angle, tanhfacs
       real(rkind), dimension(1:), intent(in) :: winf_cpu, winf_past_shock_cpu
       real(rkind), dimension(1-ng:, 1-ng:, 1-ng:, 1:), intent(inout) :: w_cpu
       real(rkind), dimension(1-ng:), intent(in) :: x_cpu, y_cpu
       integer :: i,k,l,m
       real(rkind) :: xsh, xx, tanhlen, tanhf, dwinf

       tanhlen = 8._rkind*tanhfacs

       if (ilat==1) then     ! left side
       elseif (ilat==2) then ! right side
       elseif (ilat==3) then ! lower side
       elseif (ilat==4) then  ! upper side
        do k=1,nz
         do i=1,nx
          do l=0,ng
           xsh = xshock_imp-y_cpu(ny+l)/tan(shock_angle)
           xx  = x_cpu(i)-xsh
           if (abs(xx)>tanhlen) then
            do m=1,nv
             w_cpu(i,ny+l,k,m) = w_cpu(i,ny,k,m)
            enddo
           else
            tanhf = 0.5_rkind*(1._rkind+tanh(xx/tanhfacs))
            do m=1,nv
             dwinf = winf_past_shock_cpu(m)-winf_cpu(m)
             w_cpu(i,ny+l,k,m) = winf_cpu(m)+dwinf*tanhf
            enddo
           endif
          enddo
         enddo
        enddo
       elseif (ilat==5) then  ! back side
       elseif (ilat==6) then  ! fore side
       endif
!
    endsubroutine bcshock_cpu
!
    subroutine bcextr_var_cpu(nx, ny, nz, ng, w_var_cpu)
       integer :: nx,ny,nz,ng
       integer :: ilat
       integer :: i,j,k,l,m
       real(rkind), dimension(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng, 1:1), intent(inout) :: w_var_cpu

        do k=1,nz
         do j=1,ny
          do l=1,ng
           w_var_cpu(1-l,j,k,1) = w_var_cpu(1,j,k,1)
          enddo
         enddo
        enddo
        do k=1,nz
         do j=1,ny
          do l=1,ng
           w_var_cpu(nx+l,j,k,1) = w_var_cpu(nx,j,k,1)
          enddo
         enddo
        enddo
        do k=1,nz
         do i=1,nx
          do l=1,ng
           w_var_cpu(i,1-l,k,1) = w_var_cpu(i,1,k,1)
          enddo
         enddo
        enddo
        do k=1,nz
         do i=1,nx
          do l=1,ng
           w_var_cpu(i,ny+l,k,1) = w_var_cpu(i,ny,k,1)
          enddo
         enddo
        enddo
        do j=1,ny
         do i=1,nx
          do l=1,ng
           w_var_cpu(i,j,1-l,1) = w_var_cpu(i,j,1,1)
          enddo
         enddo
        enddo
        do j=1,ny
         do i=1,nx
          do l=1,ng
           w_var_cpu(i,j,nz+l,1) = w_var_cpu(i,j,nz,1)
          enddo
         enddo
        enddo
    endsubroutine bcextr_var_cpu

    subroutine bcextr_cpu(ilat, nx, ny, nz, ng, nv, w_cpu)
       integer :: nx,ny,nz,ng,nv
       integer :: ilat
       integer :: i,j,k,l,m
       real(rkind), dimension(1-ng:, 1-ng:, 1-ng:, 1:), intent(inout) :: w_cpu

       if (ilat==1) then     ! left side
       elseif (ilat==2) then ! right side
        do k=1,nz
         do j=1,ny
          do l=1,ng
           do m=1,nv
            w_cpu(nx+l,j,k,m) = w_cpu(nx,j,k,m)
           enddo
          enddo
         enddo
        enddo
       elseif (ilat==3) then ! lower side
        do k=1,nz
         do i=1,nx
          do l=1,ng
           do m=1,nv
            w_cpu(i,1-l,k,m) = w_cpu(i,1,k,m)
           enddo
          enddo
         enddo
        enddo
       elseif (ilat==4) then  ! upper side
        do k=1,nz
         do i=1,nx
          do l=1,ng
           do m=1,nv
            w_cpu(i,ny+l,k,m) = w_cpu(i,ny,k,m)
           enddo
          enddo
         enddo
        enddo
       elseif (ilat==5) then  ! back side
        do j=1,ny
         do i=1,nx
          do l=1,ng
           do m=1,nv
            w_cpu(i,j,1-l,m) = w_cpu(i,j,1,m)
           enddo
          enddo
         enddo
        enddo
       elseif (ilat==6) then  ! fore side
        do j=1,ny
         do i=1,nx
          do l=1,ng
           do m=1,nv
            w_cpu(i,j,nz+l,m) = w_cpu(i,j,nz,m)
           enddo
          enddo
         enddo
        enddo
       endif
    endsubroutine bcextr_cpu

    subroutine bcsym_cpu(ilat, nx, ny, nz, ng, twall, w_cpu, w_aux_cpu, indx_cp_l, indx_cp_r, cv_coeff_cpu, calorically_perfect)
       integer :: nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect
       integer :: ilat
       real(rkind) :: twall
       integer :: i,k,l
       real(rkind), dimension(1-ng:, 1-ng:, 1-ng:, 1:), intent(inout) :: w_cpu
       real(rkind), dimension(1-ng:, 1-ng:, 1-ng:, 1:), intent(in) :: w_aux_cpu
       real(rkind), dimension(indx_cp_l:), intent(in) :: cv_coeff_cpu

       if (ilat==1) then     ! left side
       elseif (ilat==2) then ! right side
       elseif (ilat==3) then ! lower side
        do k=1,nz
         do i=1,nx
          do l=1,ng
           w_cpu(i,1-l,k,1) = w_cpu(i,1+l,k,1)
           w_cpu(i,1-l,k,2) = w_cpu(i,1+l,k,2)
           w_cpu(i,1-l,k,3) = w_cpu(i,1+l,k,3)
           w_cpu(i,1-l,k,4) = w_cpu(i,1+l,k,4)
           w_cpu(i,1-l,k,5) = w_cpu(i,1+l,k,5)
          enddo
         enddo
        enddo
       elseif (ilat==4) then  ! upper side
       elseif (ilat==5) then  ! back side
       elseif (ilat==6) then  ! fore side
       endif
    endsubroutine bcsym_cpu

    subroutine bcwall_cpu(ilat, nx, ny, nz, ng, twall, w_cpu, w_aux_cpu, indx_cp_l, indx_cp_r, cv_coeff_cpu, cv0, &
                          calorically_perfect, tol_iter_nr)
       integer :: nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect
       integer :: ilat
       real(rkind) :: twall, cv0, tol_iter_nr
       integer :: i,k,l
       real(rkind) :: rho,uu,vv,ww,qq,pp,tt,rhoe,ee
       real(rkind), dimension(1-ng:, 1-ng:, 1-ng:, 1:), intent(inout) :: w_cpu
       real(rkind), dimension(1-ng:, 1-ng:, 1-ng:, 1:), intent(in) :: w_aux_cpu
       real(rkind), dimension(indx_cp_l:), intent(in) :: cv_coeff_cpu

       if (ilat==1) then     ! left side
       elseif (ilat==2) then ! right side
       elseif (ilat==3) then ! lower side
        do k=1,nz
         do i=1,nx
          w_cpu(i,1,k,2) = 0._rkind
          w_cpu(i,1,k,3) = 0._rkind
          w_cpu(i,1,k,4) = 0._rkind
          !w_cpu(i,1,k,5) = w_cpu(i,1,k,1)*gm*twall
          ee = get_e_from_temperature_dev(twall, cv0, indx_cp_l, indx_cp_r, cv_coeff_cpu, calorically_perfect)
          w_cpu(i,1,k,5) = w_cpu(i,1,k,1)*ee
          do l=1,ng
           rho  = w_cpu(i,1+l,k,1)
           uu   = w_cpu(i,1+l,k,2)/rho
           vv   = w_cpu(i,1+l,k,3)/rho
           ww   = w_cpu(i,1+l,k,4)/rho
           rhoe = w_cpu(i,1+l,k,5)
           qq   = 0.5_rkind*(uu*uu+vv*vv+ww*ww)
           ee   = rhoe/rho-qq
           tt   = get_temperature_from_e_dev(ee,w_aux_cpu(i,1+l,k,6),cv0,cv_coeff_cpu,indx_cp_l,indx_cp_r,calorically_perfect,&
                                             tol_iter_nr)
           pp   = rho*tt
           tt   = 2._rkind*twall-tt ! bc
           ee   = get_e_from_temperature_dev(tt, cv0, indx_cp_l, indx_cp_r, cv_coeff_cpu,calorically_perfect)
           rho  = pp/tt
           w_cpu(i,1-l,k,1) =  rho
           w_cpu(i,1-l,k,2) = -rho*uu
           w_cpu(i,1-l,k,3) = -rho*vv
           w_cpu(i,1-l,k,4) = -rho*ww
           w_cpu(i,1-l,k,5) =  rho*(ee+qq)
          enddo
         enddo
        enddo
       elseif (ilat==4) then  ! upper side
       elseif (ilat==5) then  ! back side
       elseif (ilat==6) then  ! fore side
       endif
    endsubroutine bcwall_cpu

    subroutine bcwall_staggered_cpu(ilat, nx, ny, nz, ng, twall, w_cpu, w_aux_cpu, indx_cp_l, indx_cp_r, cv_coeff_cpu, cv0, &
                                    calorically_perfect, tol_iter_nr)
       integer :: nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect
       integer :: ilat
       real(rkind) :: twall, cv0, tol_iter_nr
       integer :: i,k,l
       real(rkind) :: rho,uu,vv,ww,qq,pp,tt,rhoe,ee
       real(rkind), dimension(1-ng:, 1-ng:, 1-ng:, 1:), intent(inout) :: w_cpu
       real(rkind), dimension(1-ng:, 1-ng:, 1-ng:, 1:), intent(in) :: w_aux_cpu
       real(rkind), dimension(indx_cp_l:), intent(in) :: cv_coeff_cpu

       if (ilat==1) then     ! left side
       elseif (ilat==2) then ! right side
       elseif (ilat==3) then ! lower side
        do k=1,nz
         do i=1,nx
          do l=1,ng
           rho  = w_cpu(i,l,k,1)
           uu   = w_cpu(i,l,k,2)/rho
           vv   = w_cpu(i,l,k,3)/rho
           ww   = w_cpu(i,l,k,4)/rho
           rhoe = w_cpu(i,l,k,5)
           qq   = 0.5_rkind*(uu*uu+vv*vv+ww*ww)
           ee   = rhoe/rho-qq
           tt   = get_temperature_from_e_dev(ee,w_aux_cpu(i,l,k,6),cv0,cv_coeff_cpu,indx_cp_l,indx_cp_r,calorically_perfect,&
                                             tol_iter_nr)
           pp   = rho*tt
           tt   = 2._rkind*twall-tt ! bc
           ee   = get_e_from_temperature_dev(tt, cv0, indx_cp_l, indx_cp_r, cv_coeff_cpu,calorically_perfect)
           rho  = pp/tt
           w_cpu(i,1-l,k,1) =  rho
           w_cpu(i,1-l,k,2) = -rho*uu
           w_cpu(i,1-l,k,3) = -rho*vv
           w_cpu(i,1-l,k,4) = -rho*ww
           w_cpu(i,1-l,k,5) =  rho*(ee+qq)
          enddo
         enddo
        enddo
       elseif (ilat==4) then  ! upper side
        do k=1,nz
         do i=1,nx
          do l=1,ng
           rho  = w_cpu(i,ny+1-l,k,1)
           uu   = w_cpu(i,ny+1-l,k,2)/rho
           vv   = w_cpu(i,ny+1-l,k,3)/rho
           ww   = w_cpu(i,ny+1-l,k,4)/rho
           rhoe = w_cpu(i,ny+1-l,k,5)
           qq   = 0.5_rkind*(uu*uu+vv*vv+ww*ww)
           ee   = rhoe/rho-qq
           tt   = get_temperature_from_e_dev(ee,w_aux_cpu(i,ny+1-l,k,6),cv0,cv_coeff_cpu,indx_cp_l,indx_cp_r,calorically_perfect,&
                                             tol_iter_nr)
           pp   = rho*tt
           tt   = 2._rkind*twall-tt ! bc
           ee   = get_e_from_temperature_dev(tt, cv0, indx_cp_l, indx_cp_r, cv_coeff_cpu,calorically_perfect)
           rho  = pp/tt
           w_cpu(i,ny+l,k,1) =  rho
           w_cpu(i,ny+l,k,2) = -rho*uu
           w_cpu(i,ny+l,k,3) = -rho*vv
           w_cpu(i,ny+l,k,4) = -rho*ww
           w_cpu(i,ny+l,k,5) =  rho*(ee+qq)
          enddo
         enddo
        enddo
       elseif (ilat==5) then  ! back side
       elseif (ilat==6) then  ! fore side
       endif
    endsubroutine bcwall_staggered_cpu

    subroutine compute_residual_cpu(nx, ny, nz, ng, nv, fln_cpu, dt, residual_rhou, fluid_mask_cpu)
        integer :: nx, ny, nz, ng, nv
        real(rkind), intent(out) :: residual_rhou
        real(rkind), intent(in) :: dt
        real(rkind), dimension(1:nx, 1:ny, 1:nz, nv), intent(in) :: fln_cpu
        integer, dimension(1-ng:,1-ng:,1-ng:), intent(in) :: fluid_mask_cpu
        integer :: i,j,k
        ! Note: should be modified to include metrics

        residual_rhou = 0._rkind
        do k=1,nz
         do j=1,ny
          do i=1,nx
           if (fluid_mask_cpu(i,j,k)==0) then
            residual_rhou = residual_rhou + (fln_cpu(i,j,k,2)/dt)**2
           endif
          enddo
         enddo
       enddo
    endsubroutine compute_residual_cpu

    subroutine compute_dt_cpu(nx, ny, nz, ng, mu0, visc_model, &
                              VISC_POWER, VISC_SUTHERLAND, powerlaw_vtexp, sutherland_S, T_ref_dim, &
                              t0, cp0, Prandtl, &
                              dcsidx_cpu, detady_cpu, dzitdz_cpu, dcsidxs_cpu, detadys_cpu, dzitdzs_cpu, w_cpu, w_aux_cpu, &
                              dtxi_max, dtyi_max, dtzi_max, dtxv_max, dtyv_max, dtzv_max, dtxk_max, dtyk_max, dtzk_max,    &
                              indx_cp_l, indx_cp_r, cp_coeff_cpu,fluid_mask_cpu,calorically_perfect)
     integer :: nx, ny, nz, ng, indx_cp_l, indx_cp_r, calorically_perfect
     real(rkind) :: mu0, cp0
     integer  :: visc_model
     real(rkind) :: powerlaw_vtexp, sutherland_S, T_ref_dim
     integer :: VISC_POWER, VISC_SUTHERLAND
     real(rkind) :: t0, Prandtl
     real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:) , intent(in) :: w_cpu 
     real(rkind), dimension(1-ng:,1-ng:,1-ng:, 1:), intent(in) :: w_aux_cpu  
     integer, dimension(1-ng:,1-ng:,1-ng:), intent(in) :: fluid_mask_cpu
     real(rkind), dimension(indx_cp_l:), intent(in) :: cp_coeff_cpu
     real(rkind), dimension(1:), intent(in) :: dcsidx_cpu, detady_cpu, dzitdz_cpu 
     real(rkind), dimension(1:), intent(in) :: dcsidxs_cpu, detadys_cpu, dzitdzs_cpu
     real(rkind) :: dtxi, dtyi, dtzi, dtxv, dtyv, dtzv, dtxk, dtyk, dtzk
     integer     :: i,j,k,ll
     real(rkind) :: rho, ri, uu, vv, ww, tt, mu, nu, k_over_rho, c
     real(rkind) :: dtxi_max, dtyi_max, dtzi_max, dtxv_max, dtyv_max, dtzv_max, dtxk_max, dtyk_max, dtzk_max
     real(rkind) :: gamloc, cploc

     dtxi_max = 0._rkind
     dtyi_max = 0._rkind
     dtzi_max = 0._rkind
     dtxv_max = 0._rkind
     dtyv_max = 0._rkind
     dtzv_max = 0._rkind
     dtxk_max = 0._rkind
     dtyk_max = 0._rkind
     dtzk_max = 0._rkind

     do k=1,nz
      do j=1,ny
       do i=1,nx
!
        if (fluid_mask_cpu(i,j,k)==0) then
!
!        rho  = w_cpu(i,j,k,1)
!        ri   = 1._rkind/rho
!        uu   = w_cpu(i,j,k,2)*ri
!        vv   = w_cpu(i,j,k,3)*ri
!        ww   = w_cpu(i,j,k,4)*ri
!        rhoe = w_cpu(i,j,k,5)
!        qq   = 0.5_rkind*(uu*uu+vv*vv+ww*ww)
!!       pp   = gm1*(rhoe-rho*qq)
!!       tt   = pp*ri
!        tt   = w_aux_cpu(i,j,k,6)
!        pp   = rho*tt
!        if(visc_model == VISC_POWER) then
!            mu = mu0 * (tt / t0)**powerlaw_vtexp
!        elseif(visc_model == VISC_SUTHERLAND) then
!            mu = mu0 * (tt / t0)**1.5_rkind * &
!                (1._rkind+sutherland_S/T_ref_dim)/(tt/t0 +
!                sutherland_S/T_ref_dim)
!        endif
!
         rho  = w_cpu(i,j,k,1)
         ri   = 1._rkind/rho
         uu   = w_aux_cpu(i,j,k,2)
         vv   = w_aux_cpu(i,j,k,3)
         ww   = w_aux_cpu(i,j,k,4)
         tt   = w_aux_cpu(i,j,k,6)
         mu   = w_aux_cpu(i,j,k,7)
!
         if (calorically_perfect==1) then
          cploc = cp0
         else
          cploc = 0._rkind
          do ll=indx_cp_l,indx_cp_r
           cploc = cploc+cp_coeff_cpu(ll)*tt**ll
          enddo
         endif
!
         nu  = ri*mu
         k_over_rho = nu*cploc/Prandtl
         gamloc = get_gamloc_dev(indx_cp_l,indx_cp_r,tt,cp_coeff_cpu,calorically_perfect)
         c      = sqrt (gamloc*tt)
         dtxi = (abs(uu)+c)*dcsidx_cpu(i) 
         dtyi = (abs(vv)+c)*detady_cpu(j) 
         dtzi = (abs(ww)+c)*dzitdz_cpu(k) 
         dtxv  = nu*dcsidxs_cpu(i)
         dtyv  = nu*detadys_cpu(j)
         dtzv  = nu*dzitdzs_cpu(k)
         dtxk  = k_over_rho*dcsidxs_cpu(i)
         dtyk  = k_over_rho*detadys_cpu(j)
         dtzk  = k_over_rho*dzitdzs_cpu(k)

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
    endsubroutine compute_dt_cpu

    subroutine compute_aux_cpu(nx, ny, nz, ng, w_cpu, w_aux_cpu, &
            visc_model, mu0, t0, sutherland_S, T_ref_dim, &
            powerlaw_vtexp, VISC_POWER, VISC_SUTHERLAND, VISC_NO, &
            cv_coeff_cpu, indx_cp_l, indx_cp_r, cv0, calorically_perfect, tol_iter_nr)
        integer(ikind), intent(in) :: nx, ny, nz, ng, visc_model
        integer(ikind), intent(in) :: VISC_POWER, VISC_SUTHERLAND, VISC_NO, indx_cp_l, indx_cp_r, calorically_perfect
        real(rkind),    intent(in) :: mu0, t0, sutherland_S, T_ref_dim, powerlaw_vtexp, cv0, tol_iter_nr
        real(rkind),    dimension(indx_cp_l:), intent(in)   :: cv_coeff_cpu
        real(rkind),    dimension(1-ng:,1-ng:,1-ng:, 1:), intent(in)   :: w_cpu  
        real(rkind),    dimension(1-ng:,1-ng:,1-ng:, 1:), intent(inout) :: w_aux_cpu
        integer(ikind)                        :: i, j, k
        real(rkind)                           :: rho, rhou, rhov, rhow, rhoe, ri, uu, vv, ww, qq, pp, tt, mu, ee
        !real(rkind) :: T_start, T_old, ee, ebar, den, num, T_pow, T_powp, tt2
        !integer :: l

        do k=1-ng, nz+ng
         do j=1-ng, ny+ng
          do i=1-ng, nx+ng
              ! w_aux(:) : rho, u, v, w, h, T, viscosity, div, |omega|, ducros
              rho  = w_cpu(i,j,k,1)
              rhou = w_cpu(i,j,k,2)
              rhov = w_cpu(i,j,k,3)
              rhow = w_cpu(i,j,k,4)
              rhoe = w_cpu(i,j,k,5)
              ri   = 1._rkind/rho
              uu   = rhou*ri
              vv   = rhov*ri
              ww   = rhow*ri
              qq   = 0.5_rkind*(uu*uu+vv*vv+ww*ww)
              ee = rhoe/rho-qq
              tt = get_temperature_from_e_dev(ee, w_aux_cpu(i,j,k,6), cv0, cv_coeff_cpu, indx_cp_l, indx_cp_r, &
                                              calorically_perfect, tol_iter_nr)
              pp = rho*tt

              w_aux_cpu(i,j,k,1) = rho
              w_aux_cpu(i,j,k,2) = uu
              w_aux_cpu(i,j,k,3) = vv
              w_aux_cpu(i,j,k,4) = ww
              w_aux_cpu(i,j,k,5) = (rhoe+pp)/rho
              w_aux_cpu(i,j,k,6) = tt

              if (visc_model == VISC_POWER) then
                  mu = mu0 * (tt / t0)**powerlaw_vtexp
              elseif (visc_model == VISC_SUTHERLAND) then
                  mu = mu0 * (tt / t0)**1.5_rkind * &
                      (1._rkind+sutherland_S/T_ref_dim)/(tt/t0 + sutherland_S/T_ref_dim)
              elseif (visc_model == VISC_NO) then
                  mu = 0._rkind
              endif
              w_aux_cpu(i,j,k,7) = mu
              ! STREAMS v1.0 : mu0 = sqgmr  ;  ggmopr = cp/Pr ; k = sqgmr * ggmopr
          enddo
         enddo
        enddo
    endsubroutine compute_aux_cpu

     function get_gamloc_dev(indx_cp_l,indx_cp_r,tt,cp_coeff_cpu,calorically_perfect)
    real(rkind) :: get_gamloc_dev
    integer :: indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind) :: tt
    real(rkind), dimension(indx_cp_l:indx_cp_r) :: cp_coeff_cpu
    real(rkind) :: cploc, gamloc
    integer :: l
    
    if (calorically_perfect==1) then
     cploc = cp_coeff_cpu(0)
    else
     cploc = 0._rkind
     do l=indx_cp_l,indx_cp_r
      cploc = cploc+cp_coeff_cpu(l)*tt**l
     enddo
    endif
    gamloc = cploc/(cploc-1._rkind)
    get_gamloc_dev = gamloc

    endfunction get_gamloc_dev

     function get_temperature_from_e_dev(ee, T_start, cv0, cv_coeff_cpu, indx_cp_l, indx_cp_r, &
                                                           calorically_perfect,tol_iter_nr)
    real(rkind) :: get_temperature_from_e_dev
    integer :: indx_cp_l, indx_cp_r, calorically_perfect
    real(rkind) :: ee, T_start, cv0, tol_iter_nr
    real(rkind), dimension(indx_cp_l:indx_cp_r) :: cv_coeff_cpu
    real(rkind) :: tt, T_old, ebar, den, num, T_pow, T_powp
    integer :: l,iter,max_iter

    max_iter = 50

    if (calorically_perfect==1) then
        tt    = ee/cv0
    else
        T_old = T_start 
        ebar  = ee - cv0
        do iter=1,max_iter
            den = 0._rkind
            num = 0._rkind
            do l=indx_cp_l,indx_cp_r
             if (l==-1) then
              T_pow  = T_old**l
              den    = den+cv_coeff_cpu(l)*T_pow
              num    = num+cv_coeff_cpu(l)*log(T_old)
             else
              T_pow  = T_old**l
              T_powp = T_old*T_pow
              den    = den+cv_coeff_cpu(l)*T_pow
              num    = num+cv_coeff_cpu(l)*(T_powp-1._rkind)/(l+1._rkind)
             endif
            enddo
            tt = T_old+(ebar-num)/den
            if (abs(tt-T_old) < tol_iter_nr) exit
            T_old = tt
        enddo
    endif
    get_temperature_from_e_dev = tt
    endfunction get_temperature_from_e_dev

     function get_e_from_temperature_dev(tt, cv0, indx_cp_l, indx_cp_r, cv_coeff_cpu, calorically_perfect)
        real(rkind) :: get_e_from_temperature_dev
        real(rkind) :: tt, cv0
        integer     :: indx_cp_l, indx_cp_r, calorically_perfect
        real(rkind), dimension(indx_cp_l:indx_cp_r) :: cv_coeff_cpu
        real(rkind) :: ee
        integer :: l

        if (calorically_perfect==1) then
            ee = tt * cv0
        else
            ee = cv0
            do l=indx_cp_l,indx_cp_r
             if (l==-1) then
              ee = ee+cv_coeff_cpu(l)*log(tt)
             else
              ee = ee+cv_coeff_cpu(l)/(l+1._rkind)*(tt**(l+1)-1._rkind)
             endif
            enddo
        endif
        get_e_from_temperature_dev = ee
    endfunction get_e_from_temperature_dev

!


    
    ! eigs33




    !NOMANAGEDsubroutine copy_to_psi_pv_managed_cpu(nxsl_ins,nxel_ins,nysl_ins,nyel_ins,nzsl_ins,nzel_ins, &
    !NOMANAGED         ng,nx,ny,nz,npsi, npsi_pv, nv_aux, n_aux_list, &
    !NOMANAGED         psi_cpu,psi_pv_managed,w_aux_cpu,aux_list_cpu)
    !NOMANAGED integer :: i,j,k,l,ll
    !NOMANAGED real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv_aux) :: w_aux_cpu

    !NOMANAGED do k=nzsl_ins,nzel_ins
    !NOMANAGED  do j=nysl_ins,nyel_ins
    !NOMANAGED   do i=nxsl_ins,nxel_ins
    !NOMANAGED    do l=1,n_aux_list
    !NOMANAGED     ll = aux_list_cpu(l)
    !NOMANAGED     psi_pv_managed(i,j,k,l)  = w_aux_cpu(i,j,k,ll)
    !NOMANAGED    enddo
    !NOMANAGED    do l=1,npsi
    !NOMANAGED     psi_pv_managed(i,j,k,n_aux_list+l) = psi_cpu(i,j,k,l)
    !NOMANAGED    enddo
    !NOMANAGED   enddo
    !NOMANAGED  enddo
    !NOMANAGED enddo
    !NOMANAGEDend subroutine copy_to_psi_pv_managed_cpu

    subroutine probe_interpolation_cpu(num_probe,nx,ny,nz,ng,nv_aux,ijk_probe_cpu,w_aux_probe_cpu,w_aux_cpu,probe_coeff_cpu)
        integer, intent(in) :: num_probe,nx,ny,nz,ng,nv_aux
        real(rkind), dimension(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:nv_aux), intent(in) :: w_aux_cpu
        real(rkind), dimension(6,num_probe), intent(inout) :: w_aux_probe_cpu
        real(rkind), dimension(2,2,2,num_probe), intent(in) :: probe_coeff_cpu
        integer, dimension(3,num_probe), intent(in) :: ijk_probe_cpu
        integer :: i,j,k,ii,jj,kk,l
        real(rkind) :: w1,w2,w3,w4,w5,w6

        do l=1,num_probe
         ii = ijk_probe_cpu(1,l)
         jj = ijk_probe_cpu(2,l)
         kk = ijk_probe_cpu(3,l)
         w1 = 0._rkind
         w2 = 0._rkind
         w3 = 0._rkind
         w4 = 0._rkind
         w5 = 0._rkind
         w6 = 0._rkind
         do k=1,2
          do j=1,2
           do i=1,2
            w1 = w1 + probe_coeff_cpu(i,j,k,l)*w_aux_cpu(i+ii-1,j+jj-1,k+kk-1,1)
            w2 = w2 + probe_coeff_cpu(i,j,k,l)*w_aux_cpu(i+ii-1,j+jj-1,k+kk-1,2)
            w3 = w3 + probe_coeff_cpu(i,j,k,l)*w_aux_cpu(i+ii-1,j+jj-1,k+kk-1,3)
            w4 = w4 + probe_coeff_cpu(i,j,k,l)*w_aux_cpu(i+ii-1,j+jj-1,k+kk-1,4)
            w5 = w5 + probe_coeff_cpu(i,j,k,l)*w_aux_cpu(i+ii-1,j+jj-1,k+kk-1,5)
            w6 = w6 + probe_coeff_cpu(i,j,k,l)*w_aux_cpu(i+ii-1,j+jj-1,k+kk-1,6)
           enddo
          enddo
         enddo
         w_aux_probe_cpu(1,l) = w1
         w_aux_probe_cpu(2,l) = w2
         w_aux_probe_cpu(3,l) = w3
         w_aux_probe_cpu(4,l) = w4
         w_aux_probe_cpu(5,l) = w5
         w_aux_probe_cpu(6,l) = w6
        enddo
!
    endsubroutine probe_interpolation_cpu

endmodule streams_kernels_cpu


