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

   filename = 'singleideal.ini'
   cfg = parse_cfg(filename)
   call cfg%get("grid","nxmax",nx)
   call cfg%get("grid","nymax",ny)
   call cfg%get("grid","ng",ng)
   call cfg%get("flow","Reynolds",Reynolds)
   call cfg%get("flow","Mach",Mach)
   call cfg%get("flow","flow_init",flow_init)
   call cfg%get("fluid","Prandtl",Prandtl)
!
   rfac  = Prandtl**(1._rkind/3._rkind)
!
   filename = 'flow_params.dat'
   flow_params_cfg = parse_cfg(filename)
   call flow_params_cfg%get("flow_params","u0",u0)
   call flow_params_cfg%get("flow_params","rho0",rho0)
   call flow_params_cfg%get("flow_params","p0",p0)
   call flow_params_cfg%get("flow_params","t0",t0)
   call flow_params_cfg%get("flow_params","c0",c0)
   call flow_params_cfg%get("flow_params","l0",l0)
   call flow_params_cfg%get("flow_params","cp0",cp0)
   call flow_params_cfg%get("flow_params","cv0",cv0)
   call flow_params_cfg%get("flow_params","mu0",mu0)
   call flow_params_cfg%get("flow_params","k0",k0)
   call flow_params_cfg%get("flow_params","gam",gam)
   call flow_params_cfg%get("flow_params","T_wall",Twall)

   !open(10,file='flow_params.dat',form='formatted')
   !read(10,*)ch,ch,u0
   !read(10,*)ch,ch,rho0
   !read(10,*)ch,ch,p0
   !read(10,*)ch,ch,t0
   !read(10,*)ch,ch,c0
   !read(10,*)ch,ch,l0
   !read(10,*)ch,ch,cp0
   !read(10,*)ch,ch,cv0
   !read(10,*)ch,ch,mu0
   !read(10,*)ch,ch,k0
   !read(10,*)ch,ch,gam
   !read(10,*)ch,ch,Twall
   !close(10)
 
   write(*,*)'Mesh size, nx= ',nx,'ny= ',ny
   write(*,*)'Reynolds = ', Reynolds

  end subroutine read_input

  subroutine read_grid_cha
  implicit none
  integer :: j

  allocate(y(1-ng:ny+ng))

  open(10,file='y.dat',form='formatted')
  do j=1-ng,ny+ng
   read(10,*)y(j)
  enddo 

  end subroutine read_grid_cha

  subroutine read_grid_bl
  implicit none
  integer :: i,j

  allocate(x(1-ng:nx+ng+1))
  allocate(y(1-ng:ny+ng))

  open(10,file='x.dat',form='formatted')
  do i=1-ng,nx+ng+1
   read(10,*)x(i)
  enddo 
!
  open(10,file='y.dat',form='formatted')
  do j=1-ng,ny+ng
   read(10,*)y(j)
  enddo 

  end subroutine read_grid_bl


  subroutine read_stat
   implicit none
   logical :: file_exist
   integer :: l

   allocate(wstat(nx,ny,nv))

   open(10,file='stat.bin',form='unformatted',access='stream')
   read(10)wstat
   close(10)
!
   if (flow_init>=1) then

    inquire(file = 'bl_profiles.dat', exist= file_exist)

    if (file_exist) then

     open(10,file='bl_profiles.dat',form='formatted')
     read(10,*)nstatloc
     allocate(ixstat(nstatloc))
     read(10,*)(ixstat(l),l=1,nstatloc)
     close(10)

    else

     write(*,*)'Error: file bl_profiles.dat does not exist'
     stop

    endif
   endif
  end subroutine read_stat

 end module reader
