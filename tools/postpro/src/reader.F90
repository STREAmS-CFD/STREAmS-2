module reader
    use cfgio_mod, only: cfg_t, parse_cfg
    use parameters
    use global_variables

contains

    subroutine read_input
        implicit none
        !
        character(20) :: ch
        character(100)            :: filename
        type(cfg_t)               :: cfg
        type(cfg_t)               :: flow_params_cfg
        type(cfg_t)               :: post_cfg

        filename = 'singleideal.ini'
        cfg = parse_cfg(filename)
        call cfg%get("grid","nxmax",nx)
        call cfg%get("grid","nymax",ny)
        call cfg%get("grid","ng",ng)
        ystag = 0
        if (cfg%has_key("grid","ystag")) then
         call cfg%get("grid","ystag",ystag)
        endif
        call cfg%get("flow","Reynolds",Reynolds_friction)
        call cfg%get("flow","Mach",Mach_input)
        call cfg%get("flow","flow_init",flow_init)
        call cfg%get("flow","theta_wall",theta_wall)
        call cfg%get("fluid","Prandtl",Prandtl)
        call cfg%get("grid","grid_dim",grid_dim)
        call cfg%get("output","io_type_w",io_type_w)
        call cfg%get("mpi","x_split",mpi_split_x)
        if (cfg%has_key("curvi","R_curv")) call cfg%get("curvi","R_curv",R_curv)
        aoa = 0._rkind
        if (cfg%has_key("flow","aoa")) call cfg%get("flow","aoa",aoa)
        nxb = nx/mpi_split_x
        !
        ! rfac  = Prandtl**(1._rkind/3._rkind)
        !
        filename = 'flow_params.dat'
        flow_params_cfg = parse_cfg(filename)
        call flow_params_cfg%get("flow_params","u0",u0)
        call flow_params_cfg%get("flow_params","p0",p0)
        call flow_params_cfg%get("flow_params","rho0",rho0)
        call flow_params_cfg%get("flow_params","t0",t0)
        call flow_params_cfg%get("flow_params","c0",c0)
        call flow_params_cfg%get("flow_params","l0",l0)
        call flow_params_cfg%get("flow_params","cp0",cp0)
        call flow_params_cfg%get("flow_params","cv0",cv0)
        call flow_params_cfg%get("flow_params","mu0",mu0)
        call flow_params_cfg%get("flow_params","k0",k0)
        call flow_params_cfg%get("flow_params","gam",gam)
        call flow_params_cfg%get("flow_params","T_wall",Twall)
        call flow_params_cfg%get("flow_params","rfac",rfac)
        call flow_params_cfg%get("flow_params","Mach",Mach)
        call flow_params_cfg%get("flow_params","Reynolds",Reynolds) ! Bulk for channel, re delta for BL

        write(*,*) 'Mesh size, nx= ',nx,'ny= ',ny
        write(*,*) 'Friction Reynolds number= ', Reynolds_friction

        filename = 'postpro.ini'
        post_cfg = parse_cfg(filename)
        call post_cfg%get("postpro","recompute_avg",recompute_avg)
        call post_cfg%get("postpro","it_start",it_start)
        call post_cfg%get("postpro","it_end",it_end)
        call post_cfg%get("postpro","it_out",it_out)
        call post_cfg%get("postpro","save_plot3d",save_plot3d)
        call post_cfg%get("postpro","plot3d_vars",plot3d_vars)
        call post_cfg%get("postpro","ixstat",ixstat)
        nstatloc = size(ixstat)
        call post_cfg%get("postpro","stat_0_1",stat_0_1)
        call post_cfg%get("postpro","npoints_bl",npoints_bl)
        allocate(delta_bl(npoints_bl-1))
        call post_cfg%get("postpro","ix_ramp_skip",ix_ramp_skip)
        call post_cfg%get("postpro","ix_out",ix_out)

    end subroutine read_input

    subroutine read_grid_cha
        implicit none
        integer :: j

        allocate(y(1-ng:ny+ng))
        allocate(yn(1:ny+1))

        open(10,file='y.dat',form='formatted')
        do j=1-ng,ny+ng
        read(10,*) y(j)
        enddo 
        close(10)

        open(10,file='yn.dat',form='formatted')
        do j=1,ny+1
        read(10,*) yn(j)
        enddo 
        close(10)
    end subroutine read_grid_cha

    subroutine read_grid_airfoil_ramp()
        implicit none
        integer :: i,j,i1,i2,i3,i4
        real(rkind) :: detj, rmod, deltay2, deltax2
        logical filecheck

        if(flow_init == 5) then
            open(18,file='wall_indices.dat')
            read(18,*) ite,ile,itu
            close(18)
        endif
        ! 
        allocate(xg(nx,ny))
        allocate(yg(nx,ny))
        allocate(dcsidxn(nx,ny))
        allocate(dcsidyn(nx,ny))
        allocate(detadxn(nx,ny))
        allocate(detadyn(nx,ny))
        allocate(dcsidx (nx,ny))
        allocate(dcsidy (nx,ny))
        allocate(detadx (nx,ny))
        allocate(detady (nx,ny))
        allocate(dxdcsin(nx,ny))
        allocate(dydcsin(nx,ny))
        allocate(dxdetan(nx,ny))
        allocate(dydetan(nx,ny))
        allocate(dxdcsi (nx,ny))
        allocate(dydcsi (nx,ny))
        allocate(dxdeta (nx,ny))
        allocate(dydeta (nx,ny))
        allocate(mcsi   (nx,ny))
        allocate(meta   (nx,ny))
        allocate(csimod (nx,ny))
        allocate(etamod (nx,ny))
        allocate(jac    (nx,ny))

        inquire(file='metrics.q', exist=filecheck)
        if (filecheck) then
         print*,'Binary metrics file found'
         open(10,file='grid2d.xyz', form="unformatted", access="stream")
          read(10) i1,i2,i3
          read(10) xg,yg
         close(10)
         open(10,file='metrics.q', form="unformatted", access="stream")
          read(10) i1,i2,i3,i4
          read(10) dxdcsi,dydcsi,dxdeta,dydeta
         close(10)
        else
         print*,'Binary metrics file not found. Try to load ascii version'
         open(10,file='metrics.dat')
         do j=1,ny
         do i=1,nx
         read(10,100) xg(i,j),yg(i,j),dxdcsi(i,j),dydcsi(i,j),dxdeta(i,j),dydeta(i,j) 
         enddo
         enddo
         close(10)
        endif

        100  format(200ES20.10)
        ! Derived metrics quantities: normalized metrics, jacobian, metric tensor components, and others.
        ! See Hung 2002
        do j=1,ny
        do i=1,nx

        !          Inverse metrics 
        detj        = 1._rkind/(dxdcsi(i,j)*dydeta(i,j)-dxdeta(i,j)*dydcsi(i,j))
        dcsidx(i,j) =  dydeta(i,j)*detj
        dcsidy(i,j) = -dxdeta(i,j)*detj
        detadx(i,j) = -dydcsi(i,j)*detj
        detady(i,j) =  dxdcsi(i,j)*detj

        !          Determinant of the jacobian matrix of the coordinate transformation
        jac(i,j)    = dcsidx(i,j)*detady(i,j)-dcsidy(i,j)*detadx(i,j)
        !        
        csimod(i,j) = sqrt(dxdcsi(i,j)**2+dydcsi(i,j)**2) ! |g_1|, module of the csi-tangent vector 
        mcsi(i,j)   = sqrt(dcsidx(i,j)**2+dcsidy(i,j)**2) ! |g^1|, module of the eta-normal vector  
        etamod(i,j) = sqrt(dxdeta(i,j)**2+dydeta(i,j)**2) ! |g_2|, module of the eta-tangent vector  
        meta(i,j)   = sqrt(detadx(i,j)**2+detady(i,j)**2) ! |g^2|, module of the csi-normal vector  

        dxdcsin(i,j) = dxdcsi(i,j)*mcsi(i,j) ! x-component of the csi-tangent unit vector
        dydcsin(i,j) = dydcsi(i,j)*mcsi(i,j) ! y-component of the csi-tangent unit vector
        dxdetan(i,j) = dxdeta(i,j)*meta(i,j) ! x-component of the eta-tangent unit vector
        dydetan(i,j) = dydeta(i,j)*meta(i,j) ! y-component of the eta-tangent unit vector

        dcsidxn(i,j) = dcsidx(i,j)*csimod(i,j) ! x-component of the eta-normal unit vector 
        dcsidyn(i,j) = dcsidy(i,j)*csimod(i,j) ! y-component of the eta-normal unit vector
        detadxn(i,j) = detadx(i,j)*etamod(i,j) ! x-component of the csi-normal unit vector 
        detadyn(i,j) = detady(i,j)*etamod(i,j) ! y-component of the csi-normal unit vector 

        enddo
        enddo

        !delta_bl(1:npoints_bl-1) = abs(yg(1,2:npoints_bl)-yg(1,1))
        delta_bl = 0._rkind
        do j=1,npoints_bl-1
         deltax2     = (xg(1,j+1)-xg(1,j))**2 
         deltay2     = (yg(1,j+1)-yg(1,j))**2 
         if (j==1) then
          delta_bl(j) = sqrt(deltax2+deltay2)
         else
          delta_bl(j) = delta_bl(j-1) + sqrt(deltax2+deltay2)
         endif
        enddo

        allocate(points_bl(npoints_bl,2))
        allocate(wstat_bl(npoints_bl,nv))
        allocate(xnf(2,nx,ny))
        xnf(1,:,:) = xg
        xnf(2,:,:) = yg

    end subroutine read_grid_airfoil_ramp

    subroutine read_grid_bl
        implicit none
        integer :: i,j

        allocate(x(1-ng:nx+ng+1))
        allocate(y(1-ng:ny+ng))
        allocate(yn(1:ny+1))

        open(10,file='x.dat',form='formatted')
        do i=1-ng,nx+ng+1
        read(10,*) x(i)
        enddo 
        close(10)

        open(10,file='y.dat',form='formatted')
        do j=1-ng,ny+ng
        read(10,*) y(j)
        enddo 
        close(10)
        if (ystag>0) then
         open(10,file='yn.dat',form='formatted')
         do j=1,ny+1
          read(10,*) yn(j)
         enddo 
         close(10)
        endif

    end subroutine read_grid_bl

    subroutine read_stat
        implicit none
        logical :: file_exist
        integer :: l, istart, iend, it, i, j, m, ib, ii
        character(4) :: chstore, chx, chy

        if (flow_init==0) nv = 70
        if (flow_init==1 .or. flow_init==2) nv = 70
        if (flow_init==4 .or. flow_init == 5 .or. &
            (flow_init == 1 .and. grid_dim==2)) then
            nv = 70
        endif

        allocate(wstat(nx,ny,nv))

        if(recompute_avg == 1) then
            print*,'Reading AVGZ/stat_z_YYY and averaging'
            allocate(wstatz(nx,ny,nv))
            wstat = 0._rkind
            do it=it_start,it_end,it_out
                print*, it
                write(chstore,1004) it
                open(11,file='AVGZ/stat_z_'//chstore//'.bin',form='unformatted',access="stream")
                read(11) wstatz
                close(11)
                wstat = wstat + wstatz
            enddo
            deallocate(wstatz)
            1004 format(I4.4)

            wstat = wstat/((real(it_end,rkind)-real(it_start,rkind))/it_out+1._rkind)

            open(11,file='z_t_averages.bin',form='unformatted')
            write(11) wstat
            close(11)

        elseif(recompute_avg == 2) then
            print*,'Time averages already performed, reading .bin'
            open(11,file='z_t_averages.bin',form='unformatted')
            read(11) wstat
            close(11)

        else
            if(io_type_w == 2) then
                print*,'Reading stat.bin'
                open(10,file='stat.bin',form='unformatted',access='stream')
                read(10) wstat
                close(10)
            elseif(io_type_w == 1) then
                print*,'Reading and merging stat?_YYY.bin'
                allocate(wstatb(nv,nxb,ny))
                ii = 0
                do ib=0,mpi_split_x-1
                    write(chx,1004) ib
                    write(chy,1004) 0
                    print *, 'stat?_'//chx//'_'//chy//'.bin'
                    if(stat_0_1 == 0) then
                        open (11,file='stat0_'//chx//'_'//chy//'.bin',form='unformatted')
                    else
                        open (11,file='stat1_'//chx//'_'//chy//'.bin',form='unformatted')
                    endif
                        read(11) wstatb(1:nv,1:nxb,1:ny)
                    close(11)
                    do j=1,ny
                    do i=1,nxb
                    do m=1,nv
                        wstat(ii+i,j,m) = wstatb(m,i,j)
                    enddo
                    enddo
                    enddo
                    ii = ii+nxb
                enddo
            endif
        endif
        !
    end subroutine read_stat

    subroutine save_plot3d_fields()
        integer :: n_vars, i, j, ii, i_var
        real(rkind) :: num, den
        n_vars = size(plot3d_vars)
        open(unit=123, file="POSTPRO/stat_fields.q", form="unformatted", access="stream")
        write(123) nx,ny,1,n_vars
        allocate(wstattemp(nx,ny))
        allocate(vtemp(nx,ny))
        do ii=1,n_vars
            i_var = plot3d_vars(ii)
            wstattemp(:,:) = wstat(:,:,i_var)
            !REMOVE if(grid_dim == 1) then
            !REMOVE     wstattemp(:,:) = wstat(i_var,:,:)
            !REMOVE elseif(grid_dim == 2) then
            !REMOVE     wstattemp(:,:) = wstat(:,:,i_var)
            !REMOVE     if(i_var == 2) then
            !REMOVE         do j=1,ny
            !REMOVE         do i=1,nx
            !REMOVE             num = wstat(i,j,3)*(dxdetan(i,j)/dcsidxn(i,j)-1._rkind)
            !REMOVE             den = dcsidyn(i,j)*dxdetan(i,j)/dcsidxn(i,j)-dydetan(i,j)
            !REMOVE             vtemp(i,j) = num/den
            !REMOVE             wstattemp(i,j) = (wstat(i,j,3)-vtemp(i,j)*dydetan(i,j))/dxdetan(i,j)
            !REMOVE         enddo
            !REMOVE         enddo
            !REMOVE     endif
            !REMOVE     if(i_var == 3) then
            !REMOVE         do j=1,ny
            !REMOVE         do i=1,nx
            !REMOVE             num = wstat(i,j,3)*(dxdetan(i,j)/dcsidxn(i,j)-1._rkind)
            !REMOVE             den = dcsidyn(i,j)*dxdetan(i,j)/dcsidxn(i,j)-dydetan(i,j)
            !REMOVE             vtemp(i,j) = num/den
            !REMOVE             wstattemp(i,j) = vtemp(i,j)
            !REMOVE         enddo
            !REMOVE         enddo
            !REMOVE     endif
            !REMOVE endif
            write(123) wstattemp
        enddo
        close(123)
    endsubroutine save_plot3d_fields

    subroutine extract_bl(i)
        real(rkind) :: vx,vy,normx,normy,nmod
        real(rkind) :: tol=1e-9
        integer :: i
        !delta_bl(:) = 0.00001_rkind
        ! find normal (counterclockwise so it is external)
        vx = xg(min(i+1,nx),1)-xg(max(i-1,1),1)
        vy = yg(min(i+1,nx),1)-yg(max(i-1,1),1)
        normx = - vy
        normy = vx 
        nmod = sqrt(normx**2+normy**2)
        normx = normx/nmod
        normy = normy/nmod
        print*,'point = ',xg(i,1),yg(i,1), vx, vy
        print*,'manual normal = ',normx,normy,normx**2+normy**2
        !normx = dxdetan(i,j)/(dxdetan(i,j)**2+dydetan(i,j)**2)**0.5
        !normy = dydetan(i,j)/(dxdetan(i,j)**2+dydetan(i,j)**2)**0.5
        !print*,'metrics normal = ',normx,normy,normx**2+normy**2

        ! sets points along normal
        points_bl(1,:) = [xg(i,1),yg(i,1)]
        do i_bl=2,npoints_bl
            points_bl(i_bl,:) = [xg(i,1),yg(i,1)] + delta_bl(i_bl-1)*[normx,normy]
            !print*,'PBL: ',i_bl, points_bl(i_bl,:)
        enddo
        !print*,'points_bl(:,1): ',points_bl(:,1)
        !print*,'points_bl(:,2): ',points_bl(:,2)

        ! interpolate quantities along bl line
        wstat_bl(1,:) = wstat(i,1,:)
        !print*,'wstat wall: ',wstat(i,1,1:5)
        ii = i
        jj = 1
        do i_bl=2,npoints_bl
            x_bl   = points_bl(i_bl,1)
            y_bl   = points_bl(i_bl,2)
            xst    = [x_bl,y_bl]
            istart = max(ii-10,1)
            iend   = min(ii+10,nx)
            i1     = istart
            i2     = iend
            j1     = max(jj-10,1)
            j2     = min(jj+10,ny)
            call solid_search(i1,i2,j1,j2,xnf(:,istart:iend,1:ny), &
                                 istart,iend,ny,0,xst,ii,jj,ierr,iprint)
            if (ierr.ne.0) then
                print*, 'Fast approach failed for i_bl,xst()', i_bl,xst
                istart = 1
                iend   = nx
                i1     = istart
                i2     = iend
                j1     = 1
                j2     = ny 
                call solid_search(i1,i2,j1,j2,xnf(:,1:nx,1:ny), &
                                  istart,iend,ny,0,xst,ii,jj,ierr,iprint)
            else
                !print*,'tank xst,ii,jj:',xst,ii,jj
            endif
            i_loc = ii
            j_loc = jj
            v_loc(1,:) = [xg(i_loc  ,j_loc  ), yg(i_loc  ,j_loc  )]
            v_loc(2,:) = [xg(i_loc+1,j_loc  ), yg(i_loc+1,j_loc  )]
            v_loc(3,:) = [xg(i_loc+1,j_loc+1), yg(i_loc+1,j_loc+1)]
            v_loc(4,:) = [xg(i_loc  ,j_loc+1), yg(i_loc  ,j_loc+1)]
            do iv=1,nv
                p_loc(1) = wstat(i_loc  ,j_loc  ,iv)
                p_loc(2) = wstat(i_loc+1,j_loc  ,iv)
                p_loc(3) = wstat(i_loc+1,j_loc+1,iv)
                p_loc(4) = wstat(i_loc  ,j_loc+1,iv)
                call int_quasibil(v_loc,p_loc,xst,p_int)
                wstat_bl(i_bl,iv) = p_int
            enddo
        enddo
   endsubroutine extract_bl

end module reader
