!< STREAmS, field class definition.
module streams_field_object

  use streams_grid_object, only : grid_object, GRID_AIRFOIL, GRID_CHA
  use streams_parameters
  use MPI
  use iso_c_binding
  use, intrinsic :: iso_fortran_env, only : error_unit
  use cpumem_mod

  implicit none
  private
  public :: field_object

  type :: field_object
    !< Field class definition.
    type(grid_object), pointer :: grid=>null()
    integer(ikind) :: nv=1_ikind
    integer(ikind) :: nx, ny, nz
    integer :: save_metrics_2d
    real(rkind), allocatable, dimension(:) :: x, y, z, yn
    real(rkind), allocatable, dimension(:) :: dcsidx,dcsidx2,dcsidxs
    real(rkind), allocatable, dimension(:) :: detady,detady2,detadys
    real(rkind), allocatable, dimension(:) :: dzitdz,dzitdz2,dzitdzs
    !
    !Coordinates and metric (curvilinear) related quantities
    real(rkind), allocatable, dimension(:,:) :: xc2, yc2
    real(rkind), allocatable, dimension(:,:) :: dcsidxc2,dcsidyc2
    real(rkind), allocatable, dimension(:,:) :: detadxc2,detadyc2
    !
    real(rkind), allocatable, dimension(:,:) :: dcsidxnc2,dcsidync2
    real(rkind), allocatable, dimension(:,:) :: detadxnc2,detadync2
    !
    real(rkind), allocatable, dimension(:,:) :: dxdcsic2,dydcsic2
    real(rkind), allocatable, dimension(:,:) :: dxdetac2,dydetac2

    real(rkind), allocatable, dimension(:,:) :: dxdcsinc2,dydcsinc2
    real(rkind), allocatable, dimension(:,:) :: dxdetanc2,dydetanc2
    !
    real(rkind), allocatable, dimension(:,:) :: jac,mcsijac1,metajac1
    real(rkind), allocatable, dimension(:,:) :: mcsi,meta,csimod,etamod
    real(rkind), allocatable, dimension(:,:) :: g1,g2,g12
    real(rkind) :: volume
    !
    real(rkind), allocatable, dimension(:,:) :: theta_ij
    !
    integer :: ite_l, itu_l, ite_rank_x, itu_rank_x
    integer :: icl_l, icu_l, icl_rank_x, icu_rank_x
    integer :: ile_l, ile_rank_x
    integer, dimension(:), allocatable :: wall_tag
    !Field variable
    real(rkind), allocatable, dimension(:,:,:,:) :: w

    real(rkind), dimension(:,:,:,:), allocatable :: grid3d

    !MPI data
    integer(ikind) :: mpi_err
    integer(ikind) :: myrank
    integer(ikind) :: nprocs
    logical :: masterproc
    integer, dimension(3) :: nblocks, ncoords
    integer :: mp_cart,mp_cartx,mp_carty,mp_cartz
    integer :: nproc,nrank_x, nrank_y, nrank_z, nrank_wake
    integer :: ileftx,irightx,ilefty,irighty,ileftz,irightz
    integer :: ileftbottom, irightbottom, ilefttop, irighttop
    !ranks for edges
    integer :: lxly, lxry, rxly, rxry
    integer :: lxlz, lxrz, rxlz, rxrz
    integer :: lylz, lyrz, rylz, ryrz
    !ranks for corners
    integer :: lxlylz, lxlyrz, lxrylz, lxryrz
    integer :: rxlylz, rxlyrz, rxrylz, rxryrz


  contains
    !public methods
    procedure, pass(self) :: initialize
    procedure, pass(self) :: alloc
    procedure, pass(self) :: cartesian_mpi
    procedure, pass(self) :: correct_bc
    procedure, pass(self) :: define_corners
    procedure, pass(self) :: define_edges_corners
    procedure, pass(self) :: compute_local_grid_metrics
    procedure, pass(self) :: set_local_grid
    procedure, pass(self) :: process_grid_par_c2
    procedure, pass(self) :: process_grid_c2
    procedure, pass(self) :: compute_metrics_c2
    procedure, pass(self) :: read_grid_c2
    procedure, pass(self) :: init_grid3d
    procedure, pass(self) :: metric_i
    procedure, pass(self) :: metric_j
    procedure, pass(self) :: save_metrics_p3d
    procedure, pass(self) :: bcswap_met
    procedure, pass(self) :: bcswap_xy
    procedure, pass(self) :: bcwake_met
    procedure, pass(self) :: bcwake_xy
    procedure, pass(self) :: bcte_met
    procedure, pass(self) :: read_field
    procedure, pass(self) :: read_field_serial
    procedure, pass(self) :: write_field
    procedure, pass(self) :: write_field_serial
    procedure, pass(self) :: write_plot3d
    procedure, pass(self) :: read_plot3d
    procedure, pass(self) :: write_vtk
    procedure, pass(self) :: write_vtk_c1
    procedure, pass(self) :: write_vtk_c2
    procedure, pass(self) :: check_cpu_mem
    procedure, pass(self) :: gen_grid_cha_curv
  endtype field_object

contains

  subroutine check_cpu_mem(self, description)
    class(field_object), intent(inout) :: self !< The field backend.
    character(*) :: description
    integer :: ierr
    character(128) :: proc_name
    integer :: resultlen

    integer(c_long) :: mem_total, mem_av
    real(rkind) :: memtotal, memav

    !call mpi_barrier(mpi_comm_world, ierr)
    call mpi_get_processor_name(proc_name, resultlen, ierr)

    call getmemory(mem_total,mem_av)
    memtotal = real(mem_total,rkind)/(1024_rkind**3)
    memav = real(mem_av,rkind)/(1024_rkind**3)
    !call mpi_barrier(mpi_comm_world, ierr)
    write(error_unit, "(A,2x,A,2x,A,2x,I0,2x,I0,2x,I0)") 'CPU rank,mems: ', description,&
    &proc_name(1:resultlen),self%myrank, mem_av, mem_total
  endsubroutine check_cpu_mem

  subroutine initialize(self, grid, nv, mpi_splits, save_metrics_2d)
    class(field_object), intent(inout) :: self !< The field.
    type(grid_object), intent(in), target :: grid !< Grid data.
    integer(ikind), intent(in), optional :: nv !< Number of field variables.
    integer, dimension(3) :: mpi_splits
    integer :: i, j, save_metrics_2d

    self%grid => grid
    self%nv = nv
    self%save_metrics_2d = save_metrics_2d

    call get_mpi_basic_info(self%nprocs, self%myrank, self%masterproc, self%mpi_err)

    call self%cartesian_mpi(mpi_splits)

    call self%alloc()

    if(self%grid%grid_dim == 2) then
      if (self%grid%grid2d_par == 0) then
        call self%process_grid_c2()
        else ! build grid2d if not already built (i.e., grid2d_par>0), set boundaries and set x/y-wall and wall_tag
        if (self%grid%grid_type == GRID_CHA) then
          call self%gen_grid_cha_curv(self%grid%dyp_target, self%grid%retaucha, self%grid%angle, self%grid%R_curv)
        else
          call self%read_grid_c2() ! read xc2, yc2, grid%xwall, grid%ywall, grid%xg(:), grid%yg(:)
          call self%process_grid_par_c2()
        endif
        call self%grid%compute_metrics()
      endif
    endif

    !Decompose 1d-grid and compute metrics (when grid_dim=2 some 1d grids are used: often z, x/y for init or other tasks)
    call self%set_local_grid()

    !do i=1-self%grid%ng,self%grid%nxmax+self%grid%ng+1
      !write(170+self%myrank,*) self%grid%xg(i)
    !enddo
    !do j=1-self%grid%ng,self%grid%nymax+self%grid%ng
      !write(180+self%myrank,*) self%grid%yg(j)
    !enddo
    !call mpi_barrier(mpi_comm_world, i)

  endsubroutine initialize

  subroutine init_grid3d(self)
    class(field_object), intent(inout) :: self !< The field.
    integer :: i,j,k
    associate(grid3d => self%grid3d)
      if (self%grid%grid_dim == 1) then
        do k=1,self%nz
          do j=1,self%ny
            do i=1,self%nx
              grid3d(1,i,j,k) = self%x(i)
              grid3d(2,i,j,k) = self%y(j)
              grid3d(3,i,j,k) = self%z(k)
            enddo
          enddo
        enddo
      elseif (self%grid%grid_dim == 2) then
        do k=1,self%nz
          do j=1,self%ny
            do i=1,self%nx
              grid3d(1,i,j,k) = self%xc2(i,j)
              grid3d(2,i,j,k) = self%yc2(i,j)
              grid3d(3,i,j,k) = self%z(k)
            enddo
          enddo
        enddo
        !print *, 'minmax grid1', MINVAL(grid3d(1,:,:,:)), MAXVAL(grid3d(1,:,:,:))
        !print *, 'minmax grid2', MINVAL(grid3d(2,:,:,:)), MAXVAL(grid3d(2,:,:,:))
        !print *, 'minmax grid3', MINVAL(grid3d(3,:,:,:)), MAXVAL(grid3d(3,:,:,:))
      endif
    endassociate
  endsubroutine init_grid3d

  subroutine compute_local_grid_metrics(self, coeff_deriv1, ep_ord_change, lmax_tag)
    class(field_object), intent(inout) :: self !< The field.
    real(rkind), dimension(1:4,4) :: coeff_deriv1
    integer, dimension(0:,0:,0:,1:) :: ep_ord_change
    integer, dimension(1-self%grid%ng:) :: lmax_tag

    if(self%grid%grid_dim == 2) then
      call self%compute_metrics_c2(coeff_deriv1, ep_ord_change, lmax_tag, save_metrics=self%save_metrics_2d)
    endif

    call self%init_grid3d()

  endsubroutine compute_local_grid_metrics

  subroutine cartesian_mpi(self, mpi_splits)
    class(field_object), intent(inout) :: self !< The field.
    integer, dimension(3) :: mpi_splits
    integer, parameter :: ndims=3
    logical :: remain_dims(ndims)
    logical :: reord

    associate(nblocks => self%nblocks, ncoords => self%ncoords, nprocs => self%nprocs,&
    & mpi_err => self%mpi_err, pbc => self%grid%is_xyz_periodic, masterproc => self%masterproc,&
    & mp_cart=>self%mp_cart,mp_cartx=>self%mp_cartx,mp_carty=>self%mp_carty,mp_cartz=> self%mp_cartz,&
    & nrank=>self%myrank,nproc=>self%nprocs,nrank_x=>self%nrank_x, nrank_y=>self%nrank_y,&
    & nrank_z=>self%nrank_z, ileftx=>self%ileftx,irightx=>self%irightx,ilefty=>self%ilefty,&
    &irighty=>self%irighty,ileftz=>self%ileftz,irightz=>self%irightz, nrank_wake => self%nrank_wake)

      nblocks(:) = mpi_splits(:)

      self%nx = self%grid%nxmax/nblocks(1)
      self%ny = self%grid%nymax/nblocks(2)
      self%nz = self%grid%nzmax/nblocks(3)

      !Create 3D topology
      reord = .false.
      call mpi_cart_create(mpi_comm_world,ndims,nblocks,pbc,reord,mp_cart,mpi_err)
      call mpi_cart_coords(mp_cart,nrank,ndims,ncoords,mpi_err)

      !Create 1D communicators
      remain_dims(1) = .true.
      remain_dims(2) = .false.
      remain_dims(3) = .false.
      call mpi_cart_sub(mp_cart,remain_dims,mp_cartx,mpi_err)
      call mpi_comm_rank(mp_cartx,nrank_x,mpi_err)
      call mpi_cart_shift(mp_cartx,0,1,ileftx,irightx,mpi_err)
      remain_dims(2) = .true.
      remain_dims(1) = .false.
      remain_dims(3) = .false.
      call mpi_cart_sub(mp_cart,remain_dims,mp_carty,mpi_err)
      call mpi_comm_rank(mp_carty,nrank_y,mpi_err)
      call mpi_cart_shift(mp_carty,0,1,ilefty,irighty,mpi_err)
      remain_dims(3) = .true.
      remain_dims(1) = .false.
      remain_dims(2) = .false.
      call mpi_cart_sub(mp_cart,remain_dims,mp_cartz,mpi_err)
      call mpi_comm_rank(mp_cartz,nrank_z,mpi_err)
      call mpi_cart_shift(mp_cartz,0,1,ileftz,irightz,mpi_err)

      !Find neighbor process in the airfoil C-grid logic
      call mpi_cart_shift_general(mp_cart, [nblocks(1)-2*ncoords(1)-1,0,0], nrank, nrank_wake)
    endassociate

    call self%define_corners()

    call self%define_edges_corners()

  endsubroutine cartesian_mpi

  subroutine define_corners(self)
    class(field_object), intent(inout) :: self
    !
    integer, dimension(3) :: ncoords_nei
    !
    associate(ileftbottom => self%ileftbottom, irightbottom => self%irightbottom,&
    & ilefttop => self%ilefttop, irighttop => self%irighttop, ncoords => self%ncoords,&
    & nblocks => self%nblocks, iermpi => self%mpi_err, mp_cart => self%mp_cart, pbc => self%grid%is_xyz_p&
    &eriodic)

      ileftbottom = mpi_proc_null
      irightbottom = mpi_proc_null
      ilefttop = mpi_proc_null
      irighttop = mpi_proc_null
      !
      !ileftbottom
      ncoords_nei(1) = ncoords(1)-1
      if (pbc(1)) then
        if (ncoords_nei(1)==-1) ncoords_nei(1) = nblocks(1)-1
        if (ncoords_nei(1)==nblocks(1)) ncoords_nei(1) = 0
      endif
      ncoords_nei(2) = ncoords(2)
      ncoords_nei(3) = ncoords(3)-1
      if (pbc(3)) then
        if (ncoords_nei(3)==-1) ncoords_nei(3) = nblocks(3)-1
        if (ncoords_nei(3)==nblocks(3)) ncoords_nei(3) = 0
      endif
      if (ncoords_nei(1)>=0.and.ncoords_nei(1)<nblocks(1)) then
        if (ncoords_nei(3)>=0.and.ncoords_nei(3)<nblocks(3)) then
          call mpi_cart_rank(mp_cart,ncoords_nei,ileftbottom,iermpi)
        endif
      endif
      !ilefttop
      ncoords_nei(1) = ncoords(1)-1
      if (pbc(1)) then
        if (ncoords_nei(1)==-1) ncoords_nei(1) = nblocks(1)-1
        if (ncoords_nei(1)==nblocks(1)) ncoords_nei(1) = 0
      endif
      ncoords_nei(2) = ncoords(2)
      ncoords_nei(3) = ncoords(3)+1
      if (pbc(3)) then
        if (ncoords_nei(3)==-1) ncoords_nei(3) = nblocks(3)-1
        if (ncoords_nei(3)==nblocks(3)) ncoords_nei(3) = 0
      endif
      if (ncoords_nei(1)>=0.and.ncoords_nei(1)<nblocks(1)) then
        if (ncoords_nei(3)>=0.and.ncoords_nei(3)<nblocks(3)) then
          call mpi_cart_rank(mp_cart,ncoords_nei,ilefttop,iermpi)
        endif
      endif
      !
      !irightbottom
      ncoords_nei(1) = ncoords(1)+1
      if (pbc(1)) then
        if (ncoords_nei(1)==-1) ncoords_nei(1) = nblocks(1)-1
        if (ncoords_nei(1)==nblocks(1)) ncoords_nei(1) = 0
      endif
      ncoords_nei(2) = ncoords(2)
      ncoords_nei(3) = ncoords(3)-1
      if (pbc(3)) then
        if (ncoords_nei(3)==-1) ncoords_nei(3) = nblocks(3)-1
        if (ncoords_nei(3)==nblocks(3)) ncoords_nei(3) = 0
      endif
      if (ncoords_nei(1)>=0.and.ncoords_nei(1)<nblocks(1)) then
        if (ncoords_nei(3)>=0.and.ncoords_nei(3)<nblocks(3)) then
          call mpi_cart_rank(mp_cart,ncoords_nei,irightbottom,iermpi)
        endif
      endif
      !irighttop
      ncoords_nei(1) = ncoords(1)+1
      if (pbc(1)) then
        if (ncoords_nei(1)==-1) ncoords_nei(1) = nblocks(1)-1
        if (ncoords_nei(1)==nblocks(1)) ncoords_nei(1) = 0
      endif
      ncoords_nei(2) = ncoords(2)
      ncoords_nei(3) = ncoords(3)+1
      if (pbc(3)) then
        if (ncoords_nei(3)==-1) ncoords_nei(3) = nblocks(3)-1
        if (ncoords_nei(3)==nblocks(3)) ncoords_nei(3) = 0
      endif
      if (ncoords_nei(1)>=0.and.ncoords_nei(1)<nblocks(1)) then
        if (ncoords_nei(3)>=0.and.ncoords_nei(3)<nblocks(3)) then
          call mpi_cart_rank(mp_cart,ncoords_nei,irighttop,iermpi)
        endif
      endif
      !
    endassociate
  end subroutine define_corners

  subroutine define_edges_corners(self)
    class(field_object), intent(inout) :: self
    !
    integer, dimension(3) :: ncoords_nei
    !
    associate(myrank => self%myrank, mp_cart => self%mp_cart, lxly => self%lxly, lxry => self%lxry,&
    & rxly => self%rxly, rxry => self%rxry, lxlz => self%lxlz, lxrz => self%lxrz, rxlz => self%rxlz,&
    & rxrz => self%rxrz, lylz => self%lylz, lyrz => self%lyrz, rylz => self%rylz, ryrz => self%ryrz,&
    & lxlylz => self%lxlylz, lxlyrz => self%lxlyrz, lxrylz => self%lxrylz, lxryrz => self%lxryrz,&
    & rxlylz => self%rxlylz, rxlyrz => self%rxlyrz, rxrylz => self%rxrylz, rxryrz => self%rxryrz)

      call mpi_cart_shift_general(mp_cart, [-1,-1, 0], myrank, lxly)
      call mpi_cart_shift_general(mp_cart, [-1, 1, 0], myrank, lxry)
      call mpi_cart_shift_general(mp_cart, [ 1,-1, 0], myrank, rxly)
      call mpi_cart_shift_general(mp_cart, [ 1, 1, 0], myrank, rxry)

      call mpi_cart_shift_general(mp_cart, [-1, 0, -1], myrank, lxlz)
      call mpi_cart_shift_general(mp_cart, [-1, 0, 1], myrank, lxrz)
      call mpi_cart_shift_general(mp_cart, [ 1, 0, -1], myrank, rxlz)
      call mpi_cart_shift_general(mp_cart, [ 1, 0, 1], myrank, rxrz)

      call mpi_cart_shift_general(mp_cart, [0 ,-1, -1], myrank, lylz)
      call mpi_cart_shift_general(mp_cart, [0 ,-1, 1], myrank, lyrz)
      call mpi_cart_shift_general(mp_cart, [0 , 1, -1], myrank, rylz)
      call mpi_cart_shift_general(mp_cart, [0 , 1, 1], myrank, ryrz)

      call mpi_cart_shift_general(mp_cart, [-1 ,-1, -1], myrank, lxlylz)
      call mpi_cart_shift_general(mp_cart, [-1 ,-1, 1], myrank, lxlyrz)
      call mpi_cart_shift_general(mp_cart, [-1 , 1, -1], myrank, lxrylz)
      call mpi_cart_shift_general(mp_cart, [-1 , 1, 1], myrank, lxryrz)
      call mpi_cart_shift_general(mp_cart, [ 1 ,-1, -1], myrank, rxlylz)
      call mpi_cart_shift_general(mp_cart, [ 1 ,-1, 1], myrank, rxlyrz)
      call mpi_cart_shift_general(mp_cart, [ 1 , 1, -1], myrank, rxrylz)
      call mpi_cart_shift_general(mp_cart, [ 1 , 1, 1], myrank, rxryrz)

    endassociate
  end subroutine define_edges_corners

  subroutine alloc(self)
    class(field_object), intent(inout) :: self !< The field.
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%grid%ng, nv => self%nv)
      allocate(self%w(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,nv))
      allocate(self%x(1-ng:nx+ng))
      allocate(self%y(1-ng:ny+ng))
      allocate(self%yn(1:ny+1))
      allocate(self%z(1-ng:nz+ng))
      allocate(self%dcsidx(nx), self%dcsidx2(nx), self%dcsidxs(nx))
      allocate(self%detady(ny), self%detady2(ny), self%detadys(ny))
      allocate(self%dzitdz(nz), self%dzitdz2(nz), self%dzitdzs(nz))

      !Coordinates and metric (curvilinear) related quantities
      if (self%grid%grid_dim == 2) then
        allocate(self%xc2(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%yc2(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%dcsidxnc2(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%dcsidync2(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%detadxnc2(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%detadync2(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%dcsidxc2(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%dcsidyc2(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%detadxc2(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%detadyc2(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%dxdcsic2(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%dydcsic2(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%dxdetac2(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%dydetac2(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%dxdcsinc2(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%dydcsinc2(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%dxdetanc2(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%dydetanc2(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%mcsi(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%meta(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%mcsijac1(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%metajac1(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%jac(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%g1(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%g2(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%g12(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%csimod(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%etamod(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%wall_tag(1-ng:nx+ng)) ! used only for flow_init=5
        self%wall_tag(:) = 0
        allocate(self%theta_ij(1:nx,1:ny))
      endif

      allocate(self%grid3d(3,nx,ny,nz))

    endassociate
  endsubroutine alloc

  subroutine set_local_grid(self)
    class(field_object), intent(inout) :: self !< The field.
    integer :: i, j, k
    integer :: ii, jj, kk

    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%grid%ng, ncoords => self%ncoor&
    &ds, x => self%x, xg => self%grid%xg, y => self%y, yg => self%grid%yg, z => self%z,&
    & zg => self%grid%zg, dcsidx => self%dcsidx, dcsidxs => self%dcsidxs, dcsidx2 => self%dcsidx2,&
    & detady => self%detady, detadys => self%detadys, detady2 => self%detady2, dzitdz => self%dzitdz,&
    & dzitdzs => self%dzitdzs, dzitdz2 => self%dzitdz2, dxg => self%grid%dxg, dyg => self%grid%dyg,&
    & dzg => self%grid%dzg, d2xg => self%grid%d2xg, d2yg => self%grid%d2yg, d2zg => self%grid%d2zg,&
    & ndim => self%grid%ndim)
      !local coordinates (nodes)
      ii = nx*ncoords(1)
      do i=1-ng,nx+ng
        x(i) = xg(ii+i)
      enddo
      jj = ny*ncoords(2)
      do j=1-ng,ny+ng
        y(j) = yg(jj+j)
      enddo
      if(ndim == 3) then
        kk = nz*ncoords(3)
        do k=1-ng,nz+ng
          z(k) = zg(kk+k)
        enddo
      else
        z(:) = 0._rkind
      endif
      !
      jj = ny*ncoords(2)
      do j=1,ny+1
        self%yn(j) = self%grid%yn(jj+j)
      enddo
      !
      ii = nx*ncoords(1)
      do i=1,nx
        dcsidx (i) = 1._rkind/(dxg(ii+i))
        dcsidxs(i) = dcsidx(i)*dcsidx(i)
        dcsidx2(i) = -d2xg(ii+i)*dcsidxs(i)
      enddo
      jj = ny*ncoords(2)
      do j=1,ny
        detady (j) = 1._rkind/dyg(jj+j)
        detadys(j) = detady(j)*detady(j)
        detady2(j) = -d2yg(jj+j)*detadys(j)
      enddo
      if(ndim == 3) then
        kk = nz*ncoords(3)
        do k=1,nz
          dzitdz (k) = 1._rkind/dzg(kk+k)
          dzitdzs(k) = dzitdz(k)*dzitdz(k)
          dzitdz2(k) = -d2zg(kk+k)*dzitdzs(k)
        enddo
      else
        dzitdz (:) = 0._rkind
        dzitdzs(:) = 0._rkind
        dzitdz2(:) = 0._rkind
      endif
    endassociate
  endsubroutine set_local_grid
  !
  subroutine process_grid_c2(self)
    class(field_object), intent(inout) :: self
    integer :: i, j, k, ii, jj, kk

    associate(nx => self%nx, ny => self%ny, ng => self%grid%ng, nxmax => self%grid%nxmax,&
    & xc2 => self%xc2, yc2 => self%yc2, xc2g => self%grid%xc2g, yc2g => self%grid%yc2g,&
    & wall_tagg => self%grid%wall_tagg, ite => self%grid%ite, itu => self%grid%itu, wall_tag => self%wall&
    &_tag, xwall => self%grid%xwall, ywall => self%grid%ywall, ile => self%grid%ile, icl => self%grid%icl&
    &, icu => self%grid%icu, ncoords => self%ncoords, ite_rank_x => self%ite_rank_x, itu_rank_x => self%i&
    &tu_rank_x, icl_rank_x => self%icl_rank_x, icu_rank_x => self%icu_rank_x, ile_l => self%ile_l,&
    & ile_rank_x => self%ile_rank_x, ite_l => self%ite_l, itu_l => self%itu_l, nrank => self%myrank,&
    & icl_l => self%icl_l, icu_l => self%icu_l, masterproc => self%masterproc)

      !Decompose global grid
      ii = nx*ncoords(1)
      jj = ny*ncoords(2)
      do j=1-ng,ny+ng
        do i=1-ng,nx+ng
          xc2(i,j) = xc2g(ii+i,jj+j)
          yc2(i,j) = yc2g(ii+i,jj+j)
        enddo
      enddo

      xwall(1-ng:nxmax+ng) = xc2g(1-ng:nxmax+ng, 1)
      ywall(1-ng:nxmax+ng) = yc2g(1-ng:nxmax+ng, 1)

      !Global to local for airfoil C-grid
      if (self%grid%grid_type == GRID_AIRFOIL) then
        !From global to local vectors
        ii = nx*ncoords(1)
        do i=1-ng,nx+ng
          wall_tag(i) = wall_tagg(ii+i)
        enddo

        !Find the local position of TE (should be put in field probably)
        ite_rank_x = (ite-1)/nx ! x-rank of the proc containing TE (pressure side)
        itu_rank_x = (itu-1)/nx ! x-rank of the proc containing TE (suction side)
        ile_rank_x = (ile-1)/nx ! x-rank of the proc containing LE
        icl_rank_x = (icl-1)/nx ! x-rank of the proc containing c/2 (pressure side)
        icu_rank_x = (icu-1)/nx ! x-rank of the proc containing c/2 (suction side)
        ite_l = ite-nx*ite_rank_x ! local position of the TE (pressure side)
        itu_l = itu-nx*itu_rank_x ! local position of the TE (suction side)
        ile_l = ile-nx*ile_rank_x ! local position of the LE
        icl_l = icl-nx*icl_rank_x ! local position of the c/2 (pressure side)
        icu_l = icu-nx*icu_rank_x ! local position of the c/2 (suction side)

        if (masterproc) print*, '======================================'
        if (masterproc) print*, 'X-ranks containing TE ',ite_rank_x,ile_rank_x,itu_rank_x
        if (masterproc) print*, 'Local nodes of the TE ',ite_l,ile_l,itu_l
      else
        ite_rank_x = MPI_PROC_NULL
        itu_rank_x = MPI_PROC_NULL
        ite_l = -100
        itu_l = -100
      endif
    endassociate

  endsubroutine process_grid_c2

  subroutine read_grid_c2(self)
    !read xc2, yc2, xwall, ywall, grid%xg(wall), grid%yg(inflow)
    class(field_object), intent(inout) :: self !< The field.
    integer :: mpi_io_file, m
    integer :: filetype
    integer,dimension(3) :: sizes ! Dimensions of the total grid
    integer,dimension(3) :: subsizes ! Dimensions of grid local to a procs
    integer,dimension(3) :: starts ! Starting coordinates
    integer :: size_real, size_integer
    integer :: ntot, iermpi
    integer (kind=mpi_offset_kind) :: offset
    integer, dimension(mpi_status_size) :: istatus
    integer :: nxmax_read,nymax_read,nzmax_one
    integer :: imin, imax, i, j
    real(rkind) :: rtemp
    character(len=5) :: i_str, k_str

    associate(mp_cartx => self%mp_cartx, nblocks => self%nblocks, nx => self%nx, ny => self%ny,&
    & ng => self%grid%ng, xc2 => self%xc2, yc2 => self%yc2, xwall => self%grid%xwall, ywall => self%grid%&
    &ywall, ncoords => self%ncoords, xg => self%grid%xg, yg => self%grid%yg, nxmax => self%grid%nxmax,&
    & nymax => self%grid%nymax)

      if(self%grid%grid2d_par == 2) then
        size_integer = storage_size(1_ikind)/8
        offset = 3*size_integer
        size_real = storage_size(1._rkind)/8
        call mpi_barrier(mp_cartx,iermpi)
        sizes(1) = nblocks(1)*nx
        sizes(2) = nblocks(2)*ny
        sizes(3) = 1
        subsizes(1) = nx
        subsizes(2) = ny
        subsizes(3) = 1
        starts(1) = 0 + self%ncoords(1)*subsizes(1)
        starts(2) = 0 + self%ncoords(2)*subsizes(2)
        starts(3) = 0
        ntot = nx*ny
        !
        call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,mpi_prec,filetype,iermpi)
        call MPI_TYPE_COMMIT(filetype,iermpi)
        call MPI_FILE_OPEN(mp_cartx,'grid2d.xyz',MPI_MODE_RDWR,MPI_INFO_NULL,mpi_io_file,iermpi)

        call MPI_FILE_SET_VIEW(mpi_io_file,offset,mpi_prec,filetype,"native",MPI_INFO_NULL,iermpi)
        call MPI_FILE_READ_ALL(mpi_io_file,self%xc2(1:nx, 1:ny),ntot,mpi_prec,istatus,iermpi)
        do m=1,nblocks(1)*nblocks(2)
          offset = offset+size_real*ntot
        enddo

        call MPI_FILE_SET_VIEW(mpi_io_file,offset,mpi_prec,filetype,"native",MPI_INFO_NULL,iermpi)
        call MPI_FILE_READ_ALL(mpi_io_file,self%yc2(1:nx, 1:ny),ntot,mpi_prec,istatus,iermpi)
        !do m=1,nblocks(1)*nblocks(2)
          !offset = offset+size_real*ntot
        !enddo

        call MPI_FILE_CLOSE(mpi_io_file,iermpi)
        call MPI_TYPE_FREE(filetype,iermpi)

      elseif(self%grid%grid2d_par == 1) then

        open(10,file='grid2d.xyz',form='unformatted',access='stream')
        read(10) nxmax_read,nymax_read,nzmax_one
        if(nxmax_read /= nxmax) then
          print*,'Error! Grid size along x does not match: ',nxmax_read,nxmax
          call mpi_abort(mpi_comm_world,88,self%mpi_err)
        endif
        if(nymax_read /= nymax) then
          print*,'Error! Grid size along y does not match: ',nymax_read,nymax
          call mpi_abort(mpi_comm_world,88,self%mpi_err)
        endif
        if(nzmax_one /= 1) then
          print*,'Error! Grid along z from file must be one: ',nzmax_one
          call mpi_abort(mpi_comm_world,88,self%mpi_err)
        endif
        imin = nx*ncoords(1)+1
        imax = imin+nx-1
        do j=1,nymax
          do i=1,nxmax
            read(10) rtemp
            if(i>=imin .and. i<=imax) then
              xc2(i-nx*ncoords(1),j) = rtemp
            endif
          enddo
        enddo
        do j=1,nymax
          do i=1,nxmax
            read(10) rtemp
            if(i>=imin .and. i<=imax) then
              yc2(i-nx*ncoords(1),j) = rtemp
            endif
          enddo
        enddo
        close(10)

      elseif(self%grid%grid2d_par == 3) then
        write(i_str(1:5),'(I5.5)') self%ncoords(1)
        write(k_str(1:5),'(I5.5)') self%ncoords(3)
        print*,'Reading grid: ',i_str," x ",k_str
        open(unit=10, file="grid_"//i_str//"/grid_"//i_str//"_"//k_str//".bin", form="unformatted", action="read")
        read(10) xc2(1:nx,1:ny)
        read(10) yc2(1:nx,1:ny)
        close(10)
        print*,'Completed reading grid: ',i_str," x ",k_str
      endif

    endassociate
  endsubroutine read_grid_c2

  subroutine process_grid_par_c2(self)
    class(field_object), intent(inout) :: self !< The field.
    integer :: i, j, k, ii, jj, kk
    integer :: im,jm,l,m,n,iend,jend,imin,imax,jmin,jmax,lmax
    integer :: ind,ind1,ind2
    real(rkind) :: theta,theta_max,theta_maxl
    real(rkind), dimension(self%nx) :: xcsil,ycsil,xetal,yetal
    real(rkind), dimension(self%grid%nxmax) :: xcsig,ycsig,xetag,yetag
    real(rkind) :: detj,rmod,vol
    real(rkind) :: dxdcsic2_sca, dydcsic2_sca, dxdetac2_sca, dydetac2_sca
    real(rkind) :: phi1, phi2, theta_ij_max, tol, errmin, error
    real(rkind), dimension(self%nx) :: etaxl, etaml, metal
    real(rkind), dimension(self%grid%nxmax) :: etaxg, etamg, metag

    tol = 1E-10

    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%grid%ng, nxmax => self%grid%nx&
    &max, volume => self%volume, c => self%grid%metrics_cd1, cc => self%grid%metrics_cd2,&
    & xg => self%grid%xg, yg => self%grid%yg, xc2 => self%xc2, yc2 => self%yc2, wall_tagg => self%grid%wa&
    &ll_tagg, ite => self%grid%ite, itu => self%grid%itu, wall_tag => self%wall_tag, nymax => self%grid%n&
    &ymax, ile => self%grid%ile, icl => self%grid%icl, icu => self%grid%icu, xwall => self%grid%xwall,&
    & ywall => self%grid%ywall, ileftx => self%ileftx, irightx=>self%irightx, ncoords => self%ncoords,&
    & ite_rank_x => self%ite_rank_x, itu_rank_x => self%itu_rank_x, icl_rank_x => self%icl_rank_x,&
    & icu_rank_x => self%icu_rank_x, ile_l => self%ile_l, ile_rank_x => self%ile_rank_x,&
    & ite_l => self%ite_l, itu_l => self%itu_l, nrank => self%myrank, icl_l => self%icl_l,&
    & icu_l => self%icu_l, mp_cart => self%mp_cart, R_curv => self%grid%R_curv, nblocks => self%nblocks,&
    & angle => self%grid%angle, mp_cartx => self%mp_cartx, iermpi => self%mpi_err, masterproc => self%mas&
    &terproc)

      !xwall, ywall, yg (without ghosts)
      call mpi_allgather(xc2(1:nx,1),nx,mpi_prec,xwall(1:nxmax),nx,mpi_prec,mp_cartx,iermpi)
      call mpi_allgather(yc2(1:nx,1),nx,mpi_prec,ywall(1:nxmax),nx,mpi_prec,mp_cartx,iermpi)
      if(self%masterproc) yg(1:nymax) = yc2(1,1:nymax)
      call MPI_Bcast(yg(1:nymax), nymax, mpi_prec, 0, mpi_comm_world, iermpi)

      !Swap x-y along x-direction
      call self%bcswap_xy()

      !Set ghost nodes in i direction (extrapolation)
      if (.not. self%grid%is_xyz_periodic(1)) then
        if(ncoords(1) == 0) then
          do j=1,nymax
            do i=1,ng
              xc2(1-i,j) = 2._rkind*xc2(2-i,j)-xc2(3-i,j)
              yc2(1-i,j) = 2._rkind*yc2(2-i,j)-yc2(3-i,j)
              !mirror xc2(1-i,j) = 2._rkind*xc2(1,j)-xc2(1+i,j)
              !mirror yc2(1-i,j) = 2._rkind*yc2(1,j)-yc2(1+i,j)
            enddo
          enddo
        endif
        if(ncoords(1) == nblocks(1)-1) then
          do j=1,nymax
            do i=1,ng
              xc2(nx+i,j) = 2._rkind*xc2(nx+i-1,j)-xc2(nx+i-2,j)
              yc2(nx+i,j) = 2._rkind*yc2(nx+i-1,j)-yc2(nx+i-2,j)
            enddo
          enddo
        endif
      endif

      !xwall, ywall, xg: extend i-ghosts
      !TODO: extrapolation should be avoided for x-periodic cases (OGRID)
      do i=1,ng
        xwall(1-i) = 2._rkind*xwall(2-i)-xwall(3-i)
        ywall(1-i) = 2._rkind*ywall(2-i)-ywall(3-i)
        xwall(nxmax+i) = 2._rkind*xwall(nxmax+i-1)-xwall(nxmax+i-2)
        ywall(nxmax+i) = 2._rkind*ywall(nxmax+i-1)-ywall(nxmax+i-2)
      enddo
      xg(1-ng:nxmax+ng) = xwall(1-ng:nxmax+ng)
      !xg has one point more (nxmax+ng+1) to init bl
      xg(nxmax+ng+1) = 2._rkind*xg(nxmax+ng) - xg(nxmax+ng-1)

      !Set ghost nodes in j direction at outer boundary
      do j=1,ng
        do i=1-ng,nx+ng
          xc2(i,nymax+j) = 2._rkind*xc2(i,nymax+j-1)-xc2(i,nymax+j-2)
          yc2(i,nymax+j) = 2._rkind*yc2(i,nymax+j-1)-yc2(i,nymax+j-2)
        enddo
      enddo

      !Setting ghost nodes in j direction (mirroring to be corrected for wake)
      do j=1,ng
        do i=1-ng,nx+ng
          !extrapolation xc2g(i,1-j) = 2._rkind*xc2g(i,2-j)-xc2g(i,3-j)
          !extrapolation yc2g(i,1-j) = 2._rkind*yc2g(i,2-j)-yc2g(i,3-j)
          xc2(i,1-j) = 2._rkind*xc2(i,1)-xc2(i,1+j)
          yc2(i,1-j) = 2._rkind*yc2(i,1)-yc2(i,1+j)
        enddo
      enddo
      !yg with ghosts (assuming no wake, may be corrected in the next part)
      do j=1,ng
        yg(1-j) = 2._rkind*yg(1)-yg(1+j)
        yg(nymax+j) = 2._rkind*yg(nymax+j-1)-yg(nymax+j-2)
      enddo

      !Init default, valid for bl (ramp)
      ite_rank_x = MPI_PROC_NULL
      itu_rank_x = MPI_PROC_NULL
      ite_l = -100
      itu_l = -100

      !do i=1-ng,nxmax+ng+1
        !write(555,*) xg(i)
      !enddo
      !do i=1-ng,nymax+ng
        !write(556,*) yg(i)
      !enddo

      if_airfoil: if (self%grid%grid_type == GRID_AIRFOIL) then

        !Find the location of the trailing edge
        i=1
        do while (i<=nxmax)
          if (abs(ywall(i)-ywall(nxmax-i+1))<tol) then
            i=i+1
          else
            exit
          endif
        enddo
        ite = i-1
        itu = nxmax-ite+1

        !Find the location of the leading edge
        ile = minloc(xwall(1:nxmax),1)

        if (masterproc) then
          open(18,file='wall_indices.dat')
          write(18,*) ite,ile,itu
          close(18)
        endif

        !wall_tag is used to distinguish airfoil wall from wake
        wall_tagg = 1
        do i=ite,ile
          wall_tagg(i) = -1 ! pressure side
        enddo
        do i=ile+1,itu
          wall_tagg(i) = 0 ! suction side
        enddo
        wall_tagg(ite) = -2
        wall_tagg(itu) = -3
        !
        if (masterproc) then
          open(18,file='tagged_wall.dat')
          do i=1-ng,nxmax+ng
            write(18,*) i,xwall(i),ywall(i),wall_tagg(i)
          enddo
          close(18)
        endif
        !
        if (masterproc) then
          print*, '======================================'
          print*, ' LE, TE, Tripping, Actuation '
          print*, '======================================'
          print*, 'Global leading edge node ', ile
          print*, 'Global trailing edge nodes ', ite,itu
        endif

        !Exchange j-min ghost nodes along wake assuming all is wake
        call self%bcwake_xy()

        !Correct j-min ghosts of yg according to wake condition
        if(self%masterproc) yg(1-ng:0) = yc2(1,1-ng:0)
        call MPI_Bcast(yg(1-ng:0), ng, mpi_prec, 0, mpi_comm_world, iermpi)

        !Setting ghost nodes in j direction (airfoil-->extrapolation)
        do j=1,ng
          do i=1-ng,nx+ng
            ii = i + nx*ncoords(1)
            if(ii > ite .and. ii < itu) then
              !extrapolation xc2g(i,1-j) = 2._rkind*xc2g(i,2-j)-xc2g(i,3-j)
              !extrapolation yc2g(i,1-j) = 2._rkind*yc2g(i,2-j)-yc2g(i,3-j)
              xc2(i,1-j) = 2._rkind*xc2(i,1)-xc2(i,1+j)
              yc2(i,1-j) = 2._rkind*yc2(i,1)-yc2(i,1+j)
            endif
          enddo
        enddo

        !open(11,file='grid_streams.xyz',form='unformatted')
        !write(11) nx+2*ng,ny+2*ng
        !write(11) (((xc2(i,j), i=1-ng,nx+ng), j=1-ng,ny+ng), l=1,1), ! (((yc2(i,j), i=1-ng,nx+ng), j=1-ng,ny+ng), l=1,1)
        !close(11)

        errmin = 1000._rkind
        do i=ite,ile ! pressure side of the airfoil
          error = abs(xwall(i)-.5_rkind)
          if (error < errmin) then
            errmin = error
            icl = i
          endif
        enddo
        errmin = 1000._rkind
        do i=ile,itu ! suction side of the airfoil
          error = abs(xwall(i)-.5_rkind)
          if (error < errmin) then
            errmin = error
            icu = i
          endif
        enddo
        !print*,'icl, icu V2 :',icl, icu

        !From global to local vectors
        ii = nx*ncoords(1)
        do i=1-ng,nx+ng
          wall_tag(i) = wall_tagg(ii+i)
        enddo

        !Find the local position of TE (should be put in field probably)
        ite_rank_x = (ite-1)/nx ! x-rank of the proc containing TE (pressure side)
        itu_rank_x = (itu-1)/nx ! x-rank of the proc containing TE (suction side)
        ile_rank_x = (ile-1)/nx ! x-rank of the proc containing LE
        icl_rank_x = (icl-1)/nx ! x-rank of the proc containing c/2 (pressure side)
        icu_rank_x = (icu-1)/nx ! x-rank of the proc containing c/2 (suction side)
        ite_l = ite-nx*ite_rank_x ! local position of the TE (pressure side)
        itu_l = itu-nx*itu_rank_x ! local position of the TE (suction side)
        ile_l = ile-nx*ile_rank_x ! local position of the LE
        icl_l = icl-nx*icl_rank_x ! local position of the c/2 (pressure side)
        icu_l = icu-nx*icu_rank_x ! local position of the c/2 (suction side)

        if (masterproc) print*, '======================================'
        if (masterproc) print*, 'X-ranks containing TE ',ite_rank_x,ile_rank_x,itu_rank_x
        if (masterproc) print*, 'Local nodes of the TE ',ite_l,ile_l,itu_l

        !Global to local for airfoil C-grid
      endif if_airfoil

    endassociate

  endsubroutine process_grid_par_c2

  subroutine compute_metrics_c2(self, coeff_deriv1, ep_ord_change, lmax_tag, save_metrics)
    class(field_object), intent(inout) :: self !< The field.
    real(rkind), dimension(1:4,4), intent(in) :: coeff_deriv1
    integer, dimension(0:,0:,0:,1:), intent(in) :: ep_ord_change
    integer, dimension(1-self%grid%ng:), intent(in) :: lmax_tag
    integer, intent(in) :: save_metrics
    integer :: i, j, k, ii, jj, kk
    integer :: im,jm,l,m,n,iend,jend,imin,imax,jmin,jmax,lmax
    integer :: ind,ind1,ind2
    real(rkind) :: theta,theta_max,theta_maxl
    real(rkind), dimension(self%nx) :: xcsil,ycsil,xetal,yetal
    real(rkind), dimension(self%grid%nxmax) :: xcsig,ycsig,xetag,yetag
    real(rkind) :: detj,rmod,vol
    real(rkind) :: dxdcsic2_sca, dydcsic2_sca, dxdetac2_sca, dydetac2_sca
    real(rkind) :: phi1, phi2, theta_ij_max, tol, errmin, error
    real(rkind), dimension(self%nx) :: etaxl, etaml, metal
    real(rkind), dimension(self%grid%nxmax) :: etaxg, etamg, metag

    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%grid%ng, nxmax => self%grid%nx&
    &max, volume => self%volume, c => self%grid%metrics_cd1, cc => self%grid%metrics_cd2,&
    & xc2 => self%xc2, yc2 => self%yc2, xc2g => self%grid%xc2g, yc2g => self%grid%yc2g,&
    & wall_tagg => self%grid%wall_tagg, ite => self%grid%ite, itu => self%grid%itu, wall_tag => self%wall&
    &_tag, nymax => self%grid%nymax, ile => self%grid%ile, icl => self%grid%icl, icu => self%grid%icu,&
    & xwall => self%grid%xwall, ywall => self%grid%ywall, dcsidxnc2 => self%dcsidxnc2,&
    & dcsidync2 => self%dcsidync2, detadxnc2 => self%detadxnc2, detadync2 => self%detadync2,&
    & dcsidxc2 => self%dcsidxc2, dcsidyc2 => self%dcsidyc2, detadxc2 => self%detadxc2,&
    & detadyc2 => self%detadyc2, dxdcsic2 => self%dxdcsic2, dydcsic2 => self%dydcsic2,&
    & dxdetac2 => self%dxdetac2, dydetac2 => self%dydetac2, dxdcsinc2 => self%dxdcsinc2,&
    & dydcsinc2 => self%dydcsinc2, dxdetanc2 => self%dxdetanc2, dydetanc2 => self%dydetanc2,&
    & mcsi => self%mcsi, meta => self%meta, mcsijac1 => self%mcsijac1, metajac1 => self%metajac1,&
    & jac => self%jac, g1 => self%g1, g2 => self%g2, g12 => self%g12, csimod => self%csimod,&
    & etamod => self%etamod, ileftx => self%ileftx, irightx=>self%irightx, ncoords => self%ncoords,&
    & ite_rank_x => self%ite_rank_x, itu_rank_x => self%itu_rank_x, icl_rank_x => self%icl_rank_x,&
    & icu_rank_x => self%icu_rank_x, ile_l => self%ile_l, ile_rank_x => self%ile_rank_x,&
    & ite_l => self%ite_l, itu_l => self%itu_l, nrank => self%myrank, icl_l => self%icl_l,&
    & icu_l => self%icu_l, mp_cart => self%mp_cart, R_curv => self%grid%R_curv, nblocks => self%nblocks,&
    & angle => self%grid%angle, theta_ij => self%theta_ij, mp_cartx => self%mp_cartx, iermpi => self%mpi_&
    &err, masterproc => self%masterproc)

      !Initialize metrics
      dxdcsic2(:,:) = 400._rkind
      dxdetac2(:,:) = 200._rkind
      dydcsic2(:,:) = 200._rkind
      dydetac2(:,:) = 400._rkind

      !Metrics (computed with the same scheme of convective terms, to guarantee freestream preservation)
      call self%metric_i(coeff_deriv1, ep_ord_change, lmax_tag)
      call self%metric_j(coeff_deriv1, ep_ord_change)

      !Ghost (BCs and inter-comm) along csi
      call self%bcswap_met() ! swap metrics between blocks (only along i)

      !Fill x-ghost metrics:
      !if periodic(x) (channel) mirror metrics, i.e., keep exchanged value but change sign to mixed metrics
      if(self%grid%is_xyz_periodic(1)) then
        if (ncoords(1).eq.0) then
          do j=1,ny
            do l=1,ng
              dxdcsic2_sca = dxdcsic2(1-l,j)
              dydcsic2_sca = dydcsic2(1-l,j)
              dxdetac2_sca = dxdetac2(1-l,j)
              dydetac2_sca = dydetac2(1-l,j)
              dxdcsic2(1-l,j) = cos(angle) * dxdcsic2_sca - sin(angle) * dydcsic2_sca
              dydcsic2(1-l,j) = sin(angle) * dxdcsic2_sca + cos(angle) * dydcsic2_sca
              dxdetac2(1-l,j) = cos(angle) * dxdetac2_sca - sin(angle) * dydetac2_sca
              dydetac2(1-l,j) = sin(angle) * dxdetac2_sca + cos(angle) * dydetac2_sca
            enddo
          enddo
        endif
        iend = nblocks(1)-1
        if (ncoords(1).eq.iend) then
          do j=1,ny
            do l=1,ng
              dxdcsic2_sca = dxdcsic2(nx+l,j)
              dydcsic2_sca = dydcsic2(nx+l,j)
              dxdetac2_sca = dxdetac2(nx+l,j)
              dydetac2_sca = dydetac2(nx+l,j)
              dxdcsic2(nx+l,j) = cos(angle) * dxdcsic2_sca + sin(angle) * dydcsic2_sca
              dydcsic2(nx+l,j) = - sin(angle) * dxdcsic2_sca + cos(angle) * dydcsic2_sca
              dxdetac2(nx+l,j) = cos(angle) * dxdetac2_sca + sin(angle) * dydetac2_sca
              dydetac2(nx+l,j) = - sin(angle) * dxdetac2_sca + cos(angle) * dydetac2_sca
            enddo
          enddo
        endif
        else ! if not periodic(x), extrapolate metrics along i (only initial & final blocks)
        if (ncoords(1).eq.0) then
          do j=1,ny
            do l=1,ng
              dxdcsic2(1-l,j) = 2._rkind*dxdcsic2(2-l,j)-dxdcsic2(3-l,j)
              dxdetac2(1-l,j) = 2._rkind*dxdetac2(2-l,j)-dxdetac2(3-l,j)
              dydcsic2(1-l,j) = 2._rkind*dydcsic2(2-l,j)-dydcsic2(3-l,j)
              dydetac2(1-l,j) = 2._rkind*dydetac2(2-l,j)-dydetac2(3-l,j)
            enddo
          enddo
        endif
        iend = nblocks(1)-1
        if (ncoords(1).eq.iend) then
          do j=1,ny
            do l=1,ng
              dxdcsic2(nx+l,j) = 2._rkind*dxdcsic2(nx+l-1,j)-dxdcsic2(nx+l-2,j)
              dxdetac2(nx+l,j) = 2._rkind*dxdetac2(nx+l-1,j)-dxdetac2(nx+l-2,j)
              dydcsic2(nx+l,j) = 2._rkind*dydcsic2(nx+l-1,j)-dydcsic2(nx+l-2,j)
              dydetac2(nx+l,j) = 2._rkind*dydetac2(nx+l-1,j)-dydetac2(nx+l-2,j)
            enddo
          enddo
        endif
      endif
      !print*,dxdcsic2(nx+1:nx+ng,3)
      !print*,dxdetac2(nx+1:nx+ng,3)
      !print*,dydcsic2(nx+1:nx+ng,3)
      !print*,dydetac2(nx+1:nx+ng,3)
      !print*,dxdcsic2(1:ng,3)
      !print*,dxdetac2(1:ng,3)
      !print*,dydcsic2(1:ng,3)
      !print*,dydetac2(1:ng,3)
      !STOP

      !Ghost along eta (only BCs since no block-y partitioning)
      !Extrapolate metrics along j
      do i=1-ng,nx+ng
        do l=1,ng
          if (wall_tag(i)<1) then ! y<0 extrapolation only on the airfoil wall
            !extrapolation dxdcsic2(i,1-l) = 2._rkind*dxdcsic2(i,2-l)-dxdcsic2(i,3-l)
            !extrapolation dxdetac2(i,1-l) = 2._rkind*dxdetac2(i,2-l)-dxdetac2(i,3-l)
            !extrapolation dydcsic2(i,1-l) = 2._rkind*dydcsic2(i,2-l)-dydcsic2(i,3-l)
            !extrapolation dydetac2(i,1-l) = 2._rkind*dydetac2(i,2-l)-dydetac2(i,3-l)

            dxdcsic2(i,1-l) = 2._rkind*dxdcsic2(i,1)-dxdcsic2(i,1+l)
            dxdetac2(i,1-l) = dxdetac2(i,1+l)
            dydcsic2(i,1-l) = 2._rkind*dydcsic2(i,1)-dydcsic2(i,1+l)
            dydetac2(i,1-l) = dydetac2(i,1+l)
          endif ! y>ly extrapolation always
          dxdcsic2(i,ny+l) = 2._rkind*dxdcsic2(i,ny+l-1)-dxdcsic2(i,ny+l-2)
          dxdetac2(i,ny+l) = 2._rkind*dxdetac2(i,ny+l-1)-dxdetac2(i,ny+l-2)
          dydcsic2(i,ny+l) = 2._rkind*dydcsic2(i,ny+l-1)-dydcsic2(i,ny+l-2)
          dydetac2(i,ny+l) = 2._rkind*dydetac2(i,ny+l-1)-dydetac2(i,ny+l-2)
        enddo
      enddo
      if (self%grid%grid_type == GRID_AIRFOIL) then
        call self%bcwake_met() ! adjust metrics in y-ghost nodes of the wake
        call self%bcte_met() ! metrics on TE must coincide
      endif

      !Derived metrics quantities: normalized metrics, jacobian, metric tensor components, and others.
      !See Hung 2002
      do j=1-ng,ny+ng
        do i=1-ng,nx+ng

          !Inverse metrics
          detj = 1._rkind/(dxdcsic2(i,j)*dydetac2(i,j)-dxdetac2(i,j)*dydcsic2(i,j))
          dcsidxc2(i,j) = dydetac2(i,j)*detj
          dcsidyc2(i,j) = -dxdetac2(i,j)*detj
          detadxc2(i,j) = -dydcsic2(i,j)*detj
          detadyc2(i,j) = dxdcsic2(i,j)*detj

          !Determinant of the jacobian matrix of the coordinate transformation
          jac(i,j) = dcsidxc2(i,j)*detadyc2(i,j)-dcsidyc2(i,j)*detadxc2(i,j)
          if (jac(i,j).le.0._rkind) then
            print *,'Negative volume',jac(i,j),'at node',i,j,'from rank', nrank
            stop
          endif
          !
          csimod(i,j) = sqrt(dxdcsic2(i,j)**2+dydcsic2(i,j)**2) ! |g_1|, module of the csi-tangent vector
          mcsi(i,j) = sqrt(dcsidxc2(i,j)**2+dcsidyc2(i,j)**2) ! |g^1|, module of the eta-normal vector
          etamod(i,j) = sqrt(dxdetac2(i,j)**2+dydetac2(i,j)**2) ! |g_2|, module of the eta-tangent vector
          meta(i,j) = sqrt(detadxc2(i,j)**2+detadyc2(i,j)**2) ! |g^2|, module of the csi-normal vector

          g1(i,j) = mcsi(i,j)**2/jac(i,j)
          g2(i,j) = meta(i,j)**2/jac(i,j)
          g12(i,j) = (dcsidxc2(i,j)*detadxc2(i,j)+dcsidyc2(i,j)*detadyc2(i,j))/jac(i,j)

          !----------------------------------------------------------------------------------------
          !normalization only for orthogonal grids
          !----------------------------------------------------------------------------------------
          !dxdcsinc2(i,j) = dxdcsic2(i,j)/csimod(i,j) ! x-component of the csi-tangent unit vector
          !dydcsinc2(i,j) = dydcsic2(i,j)/csimod(i,j) ! y-component of the csi-tangent unit vector
          !dxdetanc2(i,j) = dxdetac2(i,j)/etamod(i,j) ! x-component of the eta-tangent unit vector
          !dydetanc2(i,j) = dydetac2(i,j)/etamod(i,j) ! y-component of the eta-tangent unit vector
          !dcsidxnc2(i,j) = dcsidxc2(i,j)/mcsi(i,j) ! x-component of the eta-normal unit vector
          !dcsidync2(i,j) = dcsidyc2(i,j)/mcsi(i,j) ! y-component of the eta-normal unit vector
          !detadxnc2(i,j) = detadxc2(i,j)/meta(i,j) ! x-component of the csi-normal unit vector
          !detadync2(i,j) = detadyc2(i,j)/meta(i,j) ! y-component of the csi-normal unit vector
          !mcsijac1(i,j) = mcsi(i,j)/jac(i,j)
          !metajac1(i,j) = meta(i,j)/jac(i,j)
          !----------------------------------------------------------------------------------------
          !normalization also for non orthogonal grids
          !----------------------------------------------------------------------------------------
          dxdcsinc2(i,j) = dxdcsic2(i,j)*mcsi(i,j) ! x-component of the csi-tangent unit vector
          dydcsinc2(i,j) = dydcsic2(i,j)*mcsi(i,j) ! y-component of the csi-tangent unit vector
          dxdetanc2(i,j) = dxdetac2(i,j)*meta(i,j) ! x-component of the eta-tangent unit vector
          dydetanc2(i,j) = dydetac2(i,j)*meta(i,j) ! y-component of the eta-tangent unit vector
          dcsidxnc2(i,j) = dcsidxc2(i,j)*csimod(i,j) ! x-component of the eta-normal unit vector
          dcsidync2(i,j) = dcsidyc2(i,j)*csimod(i,j) ! y-component of the eta-normal unit vector
          detadxnc2(i,j) = detadxc2(i,j)*etamod(i,j) ! x-component of the csi-normal unit vector
          detadync2(i,j) = detadyc2(i,j)*etamod(i,j) ! y-component of the csi-normal unit vector
          mcsijac1(i,j) = 1._rkind/(jac(i,j)*csimod(i,j))
          metajac1(i,j) = 1._rkind/(jac(i,j)*etamod(i,j))
          !----------------------------------------------------------------------------------------

        enddo
      enddo

      !Save deviation of ortogonality of cells to discriminate bad cells later
      theta_ij_max = 0.
      do j=1,ny
        do i=1,nx
          theta_ij(i,j) = abs(acos(dxdetac2(i,j)/etamod(i,j))-acos(detadxc2(i,j)/meta(i,j))) * 180.d0/pi
          theta_ij_max = max(theta_ij_max,theta_ij(i,j))
        enddo
      enddo
      call mpi_allreduce(MPI_IN_PLACE,theta_ij_max,1,mpi_prec,mpi_max,mpi_comm_world,iermpi)
      if (masterproc) print*, 'Max deviation from orthogonality (first eval) [deg] =', theta_ij_max
      !
      !Write a global vector of metrics and check grid conformity
      if (save_metrics > 0) then

        call self%save_metrics_p3d()

      endif

      !Total volume, especially useful to compute bulk quantities for channel flow
      vol = 0._rkind
      do k=1,nz
        do j=1,ny
          do i=1,nx
            vol = vol + 1/jac(i,j)
          enddo
        enddo
      enddo
      call mpi_allreduce(vol,volume,1,mpi_prec,mpi_sum,mp_cart,iermpi)
      if (masterproc) print*,'Volume =', volume

      if (self%grid%grid_type == GRID_AIRFOIL.and.ncoords(3)==0) then
        if (masterproc) print*, '======================================'
        if (masterproc) print*, ' Wall spacings '
        if (masterproc) print*, '======================================'
      endif
      call mpi_barrier(mp_cart,iermpi)
      if (self%grid%grid_type == GRID_AIRFOIL.and.ncoords(3)==0) then
        if (ncoords(1)==ite_rank_x) print*, 'Delta_y at TE, pressure =', etamod(ite_l,1)
        if (ncoords(1)==ite_rank_x) print*, 'Delta_x at TE, pressure =', csimod(ite_l,1)
        if (ncoords(1)==itu_rank_x) print*, 'Delta_y at TE, suction =', etamod(itu_l,1)
        if (ncoords(1)==itu_rank_x) print*, 'Delta_x at TE, suction =', csimod(itu_l,1)
        if (ncoords(1)==icl_rank_x) print*, 'Delta_y at c/2, pressure =', etamod(icl_l,1)
        if (ncoords(1)==icl_rank_x) print*, 'Delta_x at c/2, pressure =', csimod(icl_l,1)
        if (ncoords(1)==icu_rank_x) print*, 'Delta_y at c/2, suction =', etamod(icu_l,1)
        if (ncoords(1)==icu_rank_x) print*, 'Delta_x at c/2, suction =', csimod(icu_l,1)
        if (ncoords(1)==ile_rank_x) print*, 'Delta_y at LE =', etamod(ile_l,1)
        if (ncoords(1)==ile_rank_x) print*, 'Delta_x at LE =', csimod(ile_l,1)
      endif
      call mpi_barrier(mp_cart,iermpi)

    endassociate
    !
    100 format(200ES20.10)
    200 format(1I10,1I10,40ES20.10)
    300 format(1I10,40ES20.10)

  endsubroutine compute_metrics_c2

  subroutine gen_grid_cha_curv(self, dyp_target, retaucha, angle, R_curv)
    class(field_object), intent(inout) :: self !< The grid.
    real(rkind), intent(in) :: dyp_target, retaucha, angle, R_curv
    real(rkind), dimension(1-self%grid%ng:self%grid%nxmax+self%grid%ng) :: theta
    real(rkind), dimension(1-self%grid%ng:self%grid%nymax+self%grid%ng) :: ys
    real(rkind) :: dt, theta1, deltat, dcsi, csi, b, deta, bold, eta
    integer :: i, j, ii
    !
    associate(nxmax => self%grid%nxmax,nymax=>self%grid%nymax, ng => self%grid%ng, xc2 => self%xc2,&
    & yc2 => self%yc2, yn => self%grid%yn, domain_size => self%grid%domain_size, nx => self%nx,&
    & ny=>self%ny, ncoords => self%ncoords, xg => self%grid%xg, yg => self%grid%yg, xwall => self%grid%xw&
    &all, ywall => self%grid%ywall, mp_cartx => self%mp_cartx, iermpi => self%mpi_err)

      !Arc of circumference in csi (stream-wise)
      dt = angle/nxmax
      theta1 = 0.5_rkind*pi+0.5_rkind*angle
      deltat = angle + dt*2*ng
      dcsi = 1._rkind/(nxmax+2*ng)
      do i=1-ng,nx+ng
        ii = i+ncoords(1)*nx
        csi = real(ii-1,rkind)*dcsi
        theta(i) = theta1 - csi * deltat
      enddo
      !
      !Tanh mapping function in eta (radial direction)
      b = 1._rkind
      deta = 1._rkind/nymax ! not nymax-1, because yn goes up to nymax+1
      do
        bold = b
        b = atanh(tanh(0.5_rkind*b)*(dyp_target/retaucha-1._rkind))/(deta-0.5_rkind)
        if (abs(b-bold) < tol_iter) exit
      enddo
      if (self%masterproc) write(*,*) 'Stretching parameter b =', b
      do j=1,nymax+1 ! so that ys(1) and ys(nymax) are both inside the domain
        eta = (real(j,rkind)-1._rkind)*deta
        yn(j) = tanh(b*(eta-0.5_rkind))/tanh(b*0.5_rkind) ! (from -1 to +1)
      enddo
      do j=1,nymax
        ys(j) = R_curv + 0.5_rkind*(yn(j)+yn(j+1)) ! staggered radius (from rly-1 to rly+1)
      enddo
      do j=1,ng ! symmetry conditions to fill ghost nodes
        ys(1-j) = 2._rkind*ys(1) - ys(1+j)
        ys(nymax+j) = 2._rkind*ys(nymax) - ys(nymax-j)
      enddo
      !
      !Composition of the 2D grid (including ghosts)
      do j=1-ng,nymax+ng
        do i=1-ng,nx+ng
          xc2(i,j) = ys(j)*cos(theta(i))
          yc2(i,j) = ys(j)*sin(theta(i))
        enddo
      enddo

      if (self%masterproc) then
        write(*,*) 'Stretching parameter b =', b
        write(*,*) 'Delta x+ (centerline ) =', dt*R_curv*retaucha
        write(*,*) 'L_x (centerline) =', angle*R_curv
        write(*,*) 'theta_i =', theta(1)*180/pi, 'theta_f =', theta(nx)*180/pi
        open(18,file='y.dat')
        do j=1-ng,nymax+ng
          write(18,*) ys(j) !cannot stay here, it would be out of bounds, yn(j)
        enddo
        close(18)
        open(18,file='yn.dat') ! consistent with grid_dim=1 writings
        do j=1,nymax+1
          write(18,*) yn(j)
        enddo
        close(18)
      endif

      do i=1-ng,nxmax+ng
        csi = real(i-1,rkind)*dcsi
        xwall(i) = ys(1)*cos(theta1 - csi * deltat)
        ywall(i) = ys(1)*sin(theta1 - csi * deltat)
      enddo
      xg(1-ng:nxmax+ng) = xwall(1-ng:nxmax+ng) ! the last point of xg is not assigned, not used for channel

      !yg(:) = ys(:)
      do j=1-ng,nymax+ng
        yg(j) = ys(j)*sin(theta1)
      enddo

      !do j=1-ng,ny+ng
        !do i=1-ng,nx+ng
          !write(21,*) xc2(i,j), yc2(i,j)
        !enddo
      !enddo

    endassociate
  endsubroutine gen_grid_cha_curv

  subroutine save_metrics_p3d(self)
    class(field_object), intent(inout) :: self !< The field.
    integer :: mpi_io_file, m
    integer :: filetype
    integer,dimension(3) :: sizes ! Dimensions of the total grid
    integer,dimension(3) :: subsizes ! Dimensions of grid local to a procs
    integer,dimension(3) :: starts ! Starting coordinates
    integer :: size_real, size_integer
    integer :: ntot, iermpi
    integer (kind=mpi_offset_kind) :: offset
    integer, dimension(mpi_status_size) :: istatus

    associate(mp_cartx => self%mp_cartx, nblocks => self%nblocks, nx => self%nx, ny => self%ny)

      if(self%ncoords(3) == 0) then

        size_integer = storage_size(1_ikind)/8
        size_real = storage_size(1._rkind)/8
        if (self%masterproc) then
          open(unit=123, file='metrics.q', form="unformatted", access="stream")
          write(123) self%grid%nxmax,self%grid%nymax,1,7
          flush(123)
          close(123)
        endif
        call mpi_barrier(mp_cartx,iermpi)
        sizes(1) = nblocks(1)*nx
        sizes(2) = nblocks(2)*ny
        sizes(3) = 1
        subsizes(1) = nx
        subsizes(2) = ny
        subsizes(3) = 1
        starts(1) = 0 + self%ncoords(1)*subsizes(1)
        starts(2) = 0 + self%ncoords(2)*subsizes(2)
        starts(3) = 0
        ntot = nx*ny
        !
        call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,mpi_prec,filetype,iermpi)
        call MPI_TYPE_COMMIT(filetype,iermpi)
        call MPI_FILE_OPEN(mp_cartx,'metrics.q',MPI_MODE_RDWR,MPI_INFO_NULL,mpi_io_file,iermpi)
        offset = 4*size_integer

        call MPI_FILE_SET_VIEW(mpi_io_file,offset,mpi_prec,filetype,"native",MPI_INFO_NULL,iermpi)
        call MPI_FILE_WRITE_ALL(mpi_io_file,self%dxdcsic2(1:nx, 1:ny),ntot,mpi_prec,istatus,iermpi)
        do m=1,nblocks(1)*nblocks(2)
          offset = offset+size_real*ntot
        enddo

        call MPI_FILE_SET_VIEW(mpi_io_file,offset,mpi_prec,filetype,"native",MPI_INFO_NULL,iermpi)
        call MPI_FILE_WRITE_ALL(mpi_io_file,self%dydcsic2(1:nx, 1:ny),ntot,mpi_prec,istatus,iermpi)
        do m=1,nblocks(1)*nblocks(2)
          offset = offset+size_real*ntot
        enddo

        call MPI_FILE_SET_VIEW(mpi_io_file,offset,mpi_prec,filetype,"native",MPI_INFO_NULL,iermpi)
        call MPI_FILE_WRITE_ALL(mpi_io_file,self%dxdetac2(1:nx, 1:ny),ntot,mpi_prec,istatus,iermpi)
        do m=1,nblocks(1)*nblocks(2)
          offset = offset+size_real*ntot
        enddo

        call MPI_FILE_SET_VIEW(mpi_io_file,offset,mpi_prec,filetype,"native",MPI_INFO_NULL,iermpi)
        call MPI_FILE_WRITE_ALL(mpi_io_file,self%dydetac2(1:nx, 1:ny),ntot,mpi_prec,istatus,iermpi)
        do m=1,nblocks(1)*nblocks(2)
          offset = offset+size_real*ntot
        enddo

        call MPI_FILE_SET_VIEW(mpi_io_file,offset,mpi_prec,filetype,"native",MPI_INFO_NULL,iermpi)
        call MPI_FILE_WRITE_ALL(mpi_io_file,self%etamod(1:nx, 1:ny),ntot,mpi_prec,istatus,iermpi)
        do m=1,nblocks(1)*nblocks(2)
          offset = offset+size_real*ntot
        enddo

        call MPI_FILE_SET_VIEW(mpi_io_file,offset,mpi_prec,filetype,"native",MPI_INFO_NULL,iermpi)
        call MPI_FILE_WRITE_ALL(mpi_io_file,self%meta(1:nx, 1:ny),ntot,mpi_prec,istatus,iermpi)
        do m=1,nblocks(1)*nblocks(2)
          offset = offset+size_real*ntot
        enddo

        call MPI_FILE_SET_VIEW(mpi_io_file,offset,mpi_prec,filetype,"native",MPI_INFO_NULL,iermpi)
        call MPI_FILE_WRITE_ALL(mpi_io_file,self%theta_ij,ntot,mpi_prec,istatus,iermpi)

        call MPI_FILE_CLOSE(mpi_io_file,iermpi)
        call MPI_TYPE_FREE(filetype,iermpi)

      endif
    endassociate
  endsubroutine save_metrics_p3d

  subroutine bcswap_met(self)
    !=====================================
    !Virtual interface BC for metric terms
    !=====================================
    class(field_object), intent(inout) :: self !< The field.

    integer :: i,j,ind,ind1,ind2,iermpi
    integer, dimension(mpi_status_size) :: istatus
    !
    real(rkind), dimension(self%ny*self%grid%ng) :: wbuf1s_dxdcsi
    real(rkind), dimension(self%ny*self%grid%ng) :: wbuf2s_dxdcsi
    real(rkind), dimension(self%ny*self%grid%ng) :: wbuf1r_dxdcsi
    real(rkind), dimension(self%ny*self%grid%ng) :: wbuf2r_dxdcsi
    !
    real(rkind), dimension(self%ny*self%grid%ng) :: wbuf1s_dydcsi
    real(rkind), dimension(self%ny*self%grid%ng) :: wbuf2s_dydcsi
    real(rkind), dimension(self%ny*self%grid%ng) :: wbuf1r_dydcsi
    real(rkind), dimension(self%ny*self%grid%ng) :: wbuf2r_dydcsi
    !
    real(rkind), dimension(self%ny*self%grid%ng) :: wbuf1s_dxdeta
    real(rkind), dimension(self%ny*self%grid%ng) :: wbuf2s_dxdeta
    real(rkind), dimension(self%ny*self%grid%ng) :: wbuf1r_dxdeta
    real(rkind), dimension(self%ny*self%grid%ng) :: wbuf2r_dxdeta
    !
    real(rkind), dimension(self%ny*self%grid%ng) :: wbuf1s_dydeta
    real(rkind), dimension(self%ny*self%grid%ng) :: wbuf2s_dydeta
    real(rkind), dimension(self%ny*self%grid%ng) :: wbuf1r_dydeta
    real(rkind), dimension(self%ny*self%grid%ng) :: wbuf2r_dydeta

    associate(nx => self%nx, ny => self%ny, ng => self%grid%ng, ileftx => self%ileftx,&
    & irightx => self%irightx, nrank_x => self%nrank_x, mp_cartx => self%mp_cartx, dxdcsic2 => self%dxdcs&
    &ic2, dydcsic2 => self%dydcsic2, dxdetac2 => self%dxdetac2, dydetac2 => self%dydetac2)
      !
      ind1 = 0
      do j=1,ny
        do i=1,ng
          ind1 = ind1 + 1
          wbuf1s_dxdcsi(ind1) = dxdcsic2(i,j)
          wbuf1s_dydcsi(ind1) = dydcsic2(i,j)
          wbuf1s_dxdeta(ind1) = dxdetac2(i,j)
          wbuf1s_dydeta(ind1) = dydetac2(i,j)
        enddo
      enddo
      ind2 = 0
      do j=1,ny
        do i=nx-ng+1,nx
          ind2 = ind2 + 1
          wbuf2s_dxdcsi(ind2) = dxdcsic2(i,j)
          wbuf2s_dydcsi(ind2) = dydcsic2(i,j)
          wbuf2s_dxdeta(ind2) = dxdetac2(i,j)
          wbuf2s_dydeta(ind2) = dydetac2(i,j)
        enddo
      enddo
      !
      call mpi_sendrecv(wbuf1s_dxdcsi,ind1,mpi_prec,ileftx ,1,wbuf2r_dxdcsi,ind2,mpi_prec,irightx,1,mp_cartx,istatus,iermpi)
      call mpi_sendrecv(wbuf1s_dydcsi,ind1,mpi_prec,ileftx ,1,wbuf2r_dydcsi,ind2,mpi_prec,irightx,1,mp_cartx,istatus,iermpi)
      call mpi_sendrecv(wbuf1s_dxdeta,ind1,mpi_prec,ileftx ,1,wbuf2r_dxdeta,ind2,mpi_prec,irightx,1,mp_cartx,istatus,iermpi)
      call mpi_sendrecv(wbuf1s_dydeta,ind1,mpi_prec,ileftx ,1,wbuf2r_dydeta,ind2,mpi_prec,irightx,1,mp_cartx,istatus,iermpi)
      !
      call mpi_sendrecv(wbuf2s_dxdcsi,ind2,mpi_prec,irightx,2,wbuf1r_dxdcsi,ind1,mpi_prec,ileftx ,2,mp_cartx,istatus,iermpi)
      call mpi_sendrecv(wbuf2s_dydcsi,ind2,mpi_prec,irightx,2,wbuf1r_dydcsi,ind1,mpi_prec,ileftx ,2,mp_cartx,istatus,iermpi)
      call mpi_sendrecv(wbuf2s_dxdeta,ind2,mpi_prec,irightx,2,wbuf1r_dxdeta,ind1,mpi_prec,ileftx ,2,mp_cartx,istatus,iermpi)
      call mpi_sendrecv(wbuf2s_dydeta,ind2,mpi_prec,irightx,2,wbuf1r_dydeta,ind1,mpi_prec,ileftx ,2,mp_cartx,istatus,iermpi)
      !
      ind = 0
      do j=1,ny
        do i=1,ng
          ind = ind + 1
          dxdcsic2(i-ng,j) = wbuf1r_dxdcsi(ind)
          dydcsic2(i-ng,j) = wbuf1r_dydcsi(ind)
          dxdetac2(i-ng,j) = wbuf1r_dxdeta(ind)
          dydetac2(i-ng,j) = wbuf1r_dydeta(ind)
        enddo
      enddo
      ind = 0
      do j=1,ny
        do i=1,ng
          ind = ind + 1
          dxdcsic2(nx+i,j) = wbuf2r_dxdcsi(ind)
          dydcsic2(nx+i,j) = wbuf2r_dydcsi(ind)
          dxdetac2(nx+i,j) = wbuf2r_dxdeta(ind)
          dydetac2(nx+i,j) = wbuf2r_dydeta(ind)
        enddo
      enddo

    endassociate
    !
  end subroutine bcswap_met

  subroutine bcswap_xy(self)
    !=====================================
    !Virtual interface BC for metric terms
    !=====================================
    class(field_object), intent(inout) :: self !< The field.

    integer :: i,j,ind,ind1,ind2,iermpi
    integer, dimension(mpi_status_size) :: istatus
    !
    real(rkind), dimension(self%ny*self%grid%ng) :: wbuf1s_x
    real(rkind), dimension(self%ny*self%grid%ng) :: wbuf2s_x
    real(rkind), dimension(self%ny*self%grid%ng) :: wbuf1r_x
    real(rkind), dimension(self%ny*self%grid%ng) :: wbuf2r_x
    !
    real(rkind), dimension(self%ny*self%grid%ng) :: wbuf1s_y
    real(rkind), dimension(self%ny*self%grid%ng) :: wbuf2s_y
    real(rkind), dimension(self%ny*self%grid%ng) :: wbuf1r_y
    real(rkind), dimension(self%ny*self%grid%ng) :: wbuf2r_y
    !
    associate(nx => self%nx, ny => self%ny, ng => self%grid%ng, ileftx => self%ileftx,&
    & irightx => self%irightx, nrank_x => self%nrank_x, mp_cartx => self%mp_cartx, xc2 => self%xc2,&
    & yc2 => self%yc2)
      !
      ind1 = 0
      do j=1,ny
        do i=1,ng
          ind1 = ind1 + 1
          wbuf1s_x(ind1) = xc2(i,j)
          wbuf1s_y(ind1) = yc2(i,j)
        enddo
      enddo
      ind2 = 0
      do j=1,ny
        do i=nx-ng+1,nx
          ind2 = ind2 + 1
          wbuf2s_x(ind2) = xc2(i,j)
          wbuf2s_y(ind2) = yc2(i,j)
        enddo
      enddo
      !
      call mpi_sendrecv(wbuf1s_x,ind1,mpi_prec,ileftx ,1,wbuf2r_x,ind2,mpi_prec,irightx,1,mp_cartx,istatus,iermpi)
      call mpi_sendrecv(wbuf1s_y,ind1,mpi_prec,ileftx ,1,wbuf2r_y,ind2,mpi_prec,irightx,1,mp_cartx,istatus,iermpi)
      !
      call mpi_sendrecv(wbuf2s_x,ind2,mpi_prec,irightx,2,wbuf1r_x,ind1,mpi_prec,ileftx ,2,mp_cartx,istatus,iermpi)
      call mpi_sendrecv(wbuf2s_y,ind2,mpi_prec,irightx,2,wbuf1r_y,ind1,mpi_prec,ileftx ,2,mp_cartx,istatus,iermpi)
      !
      ind = 0
      do j=1,ny
        do i=1,ng
          ind = ind + 1
          xc2(i-ng,j) = wbuf1r_x(ind)
          yc2(i-ng,j) = wbuf1r_y(ind)
        enddo
      enddo
      ind = 0
      do j=1,ny
        do i=1,ng
          ind = ind + 1
          xc2(nx+i,j) = wbuf2r_x(ind)
          yc2(nx+i,j) = wbuf2r_y(ind)
        enddo
      enddo

    endassociate
    !
  end subroutine bcswap_xy

  subroutine bcwake_met(self)
    !===================================================
    !Virtual interface BC for metrics in the wake region
    !===================================================
    !
    class(field_object), intent(inout) :: self !< The field.

    integer :: i,j,l, ind, iermpi
    integer, dimension(mpi_status_size) :: istatus
    !
    real(rkind), dimension(self%nx,self%grid%ng) :: wbufs_dxdcsi
    real(rkind), dimension(self%nx,self%grid%ng) :: wbufs_dydcsi
    real(rkind), dimension(self%nx,self%grid%ng) :: wbufs_dxdeta
    real(rkind), dimension(self%nx,self%grid%ng) :: wbufs_dydeta
    real(rkind), dimension(self%nx,self%grid%ng) :: wbufr_dxdcsi
    real(rkind), dimension(self%nx,self%grid%ng) :: wbufr_dydcsi
    real(rkind), dimension(self%nx,self%grid%ng) :: wbufr_dxdeta
    real(rkind), dimension(self%nx,self%grid%ng) :: wbufr_dydeta
    !
    associate(nx => self%nx, ny => self%ny, ng => self%grid%ng, wall_tag => self%wall_tag,&
    & ileftx => self%ileftx, irightx => self%irightx, nrank_x => self%nrank_x, nrank => self%myrank,&
    & nrank_wake => self%nrank_wake, mp_cart => self%mp_cart, dxdcsic2 => self%dxdcsic2,&
    & dydcsic2 => self%dydcsic2, dxdetac2 => self%dxdetac2, dydetac2 => self%dydetac2)

      ind = nx*ng
      do j=1,ng
        do i=1,nx
          wbufs_dxdcsi(i,j) = -dxdcsic2(nx-i+1,1+j)
          wbufs_dydcsi(i,j) = -dydcsic2(nx-i+1,1+j)
          wbufs_dxdeta(i,j) = -dxdetac2(nx-i+1,1+j)
          wbufs_dydeta(i,j) = -dydetac2(nx-i+1,1+j)
        enddo
      enddo
      !
      if (nrank_wake==nrank) then
        do j=1,ng
          do i=1,nx
            wbufr_dxdcsi(i,j) = wbufs_dxdcsi(i,j)
            wbufr_dydcsi(i,j) = wbufs_dydcsi(i,j)
            wbufr_dxdeta(i,j) = wbufs_dxdeta(i,j)
            wbufr_dydeta(i,j) = wbufs_dydeta(i,j)
          enddo
        enddo
      else
        call mpi_sendrecv(wbufs_dxdcsi,ind,mpi_prec,nrank_wake,1,wbufr_dxdcsi,ind,mpi_prec,nrank_wake,1,mp_cart,istatus,iermpi)
        call mpi_sendrecv(wbufs_dydcsi,ind,mpi_prec,nrank_wake,2,wbufr_dydcsi,ind,mpi_prec,nrank_wake,2,mp_cart,istatus,iermpi)
        call mpi_sendrecv(wbufs_dxdeta,ind,mpi_prec,nrank_wake,3,wbufr_dxdeta,ind,mpi_prec,nrank_wake,3,mp_cart,istatus,iermpi)
        call mpi_sendrecv(wbufs_dydeta,ind,mpi_prec,nrank_wake,4,wbufr_dydeta,ind,mpi_prec,nrank_wake,4,mp_cart,istatus,iermpi)
      endif
      !
      do j=1,ng
        do i=1,nx
          if (wall_tag(i)>0) then
            dxdcsic2(i,1-j) = wbufr_dxdcsi(i,j)
            dydcsic2(i,1-j) = wbufr_dydcsi(i,j)
            dxdetac2(i,1-j) = wbufr_dxdeta(i,j)
            dydetac2(i,1-j) = wbufr_dydeta(i,j)
          endif
        enddo
      enddo
      !
    endassociate
    !
  endsubroutine bcwake_met

  subroutine bcwake_xy(self)
    !
    class(field_object), intent(inout) :: self !< The field.

    integer :: i,j,l, ind, iermpi
    integer, dimension(mpi_status_size) :: istatus
    !
    real(rkind), dimension(1-self%grid%ng:self%nx+self%grid%ng,self%grid%ng) :: wbufs_x
    real(rkind), dimension(1-self%grid%ng:self%nx+self%grid%ng,self%grid%ng) :: wbufs_y
    real(rkind), dimension(1-self%grid%ng:self%nx+self%grid%ng,self%grid%ng) :: wbufr_x
    real(rkind), dimension(1-self%grid%ng:self%nx+self%grid%ng,self%grid%ng) :: wbufr_y
    !
    associate(nx => self%nx, ny => self%ny, ng => self%grid%ng, wall_tag => self%wall_tag,&
    & ileftx => self%ileftx, irightx => self%irightx, nrank_x => self%nrank_x, nrank => self%myrank,&
    & nrank_wake => self%nrank_wake, mp_cart => self%mp_cart, xc2 => self%xc2, yc2 => self%yc2,&
    &dxdcsic2 => self%dxdcsic2, dydcsic2 => self%dydcsic2, dxdetac2 => self%dxdetac2, dydetac2 => self%dy&
    &detac2)

      ind = (nx+2*ng)*ng
      do j=1,ng
        do i=1-ng,nx+ng
          wbufs_x(i,j) = xc2(nx-i+1,1+j)
          wbufs_y(i,j) = yc2(nx-i+1,1+j)
        enddo
      enddo
      !
      if (nrank_wake==nrank) then
        do j=1,ng
          do i=1-ng,nx+ng
            wbufr_x(i,j) = wbufs_x(i,j)
            wbufr_y(i,j) = wbufs_y(i,j)
          enddo
        enddo
      else
        call mpi_sendrecv(wbufs_x,ind,mpi_prec,nrank_wake,1,wbufr_x,ind,mpi_prec,nrank_wake,1,mp_cart,istatus,iermpi)
        call mpi_sendrecv(wbufs_y,ind,mpi_prec,nrank_wake,2,wbufr_y,ind,mpi_prec,nrank_wake,2,mp_cart,istatus,iermpi)
      endif
      !
      do j=1,ng
        do i=1-ng,nx+ng
          xc2(i,1-j) = wbufr_x(i,j)
          yc2(i,1-j) = wbufr_y(i,j)
        enddo
      enddo
      !
    endassociate
    !
  end subroutine bcwake_xy
  !
  subroutine bcte_met(self)
    class(field_object), intent(inout) :: self !< The field.
    !==================================
    !Enforcing TE treatment for metrics
    !==================================
    !
    integer :: i,j,l
    integer :: ind, iermpi
    integer, dimension(mpi_status_size) :: istatus
    !
    real(rkind), dimension(self%nx) :: tes_dxdcsi
    real(rkind), dimension(self%nx) :: tes_dydcsi
    real(rkind), dimension(self%nx) :: tes_dxdeta
    real(rkind), dimension(self%nx) :: tes_dydeta
    real(rkind), dimension(self%nx) :: ter_dxdcsi
    real(rkind), dimension(self%nx) :: ter_dydcsi
    real(rkind), dimension(self%nx) :: ter_dxdeta
    real(rkind), dimension(self%nx) :: ter_dydeta
    !
    associate(nx => self%nx, ny => self%ny, ng => self%grid%ng, ncoords => self%ncoords,&
    & mp_cart => self%mp_cart, ileftx => self%ileftx, irightx => self%irightx, nrank_x => self%nrank_x,&
    & itu_l => self%itu_l, itu_rank_x => self%itu_rank_x, nrank => self%myrank, nrank_wake => self%nrank_&
    &wake, dxdcsic2 => self%dxdcsic2, dydcsic2 => self%dydcsic2, dxdetac2 => self%dxdetac2,&
    & dydetac2 => self%dydetac2)

      ind = nx
      do i=1,nx
        tes_dxdcsi(i) =-dxdcsic2(nx-i+1,1)
        tes_dydcsi(i) =-dydcsic2(nx-i+1,1)
        tes_dxdeta(i) =-dxdetac2(nx-i+1,1)
        tes_dydeta(i) =-dydetac2(nx-i+1,1)
      enddo
      !
      if (nrank_wake==nrank) then
        do i=1,nx
          ter_dxdcsi(i) = tes_dxdcsi(i)
          ter_dydcsi(i) = tes_dydcsi(i)
          ter_dxdeta(i) = tes_dxdeta(i)
          ter_dydeta(i) = tes_dydeta(i)
        enddo
      else
        call mpi_sendrecv(tes_dxdcsi,ind,mpi_prec,nrank_wake,1,ter_dxdcsi,ind,mpi_prec,nrank_wake,1,mp_cart,istatus,iermpi)
        call mpi_sendrecv(tes_dydcsi,ind,mpi_prec,nrank_wake,2,ter_dydcsi,ind,mpi_prec,nrank_wake,2,mp_cart,istatus,iermpi)
        call mpi_sendrecv(tes_dxdeta,ind,mpi_prec,nrank_wake,3,ter_dxdeta,ind,mpi_prec,nrank_wake,3,mp_cart,istatus,iermpi)
        call mpi_sendrecv(tes_dydeta,ind,mpi_prec,nrank_wake,4,ter_dydeta,ind,mpi_prec,nrank_wake,4,mp_cart,istatus,iermpi)
      endif
      !
      if (ncoords(1)==itu_rank_x) then
        dxdcsic2(itu_l,1) = ter_dxdcsi(itu_l)
        dydcsic2(itu_l,1) = ter_dydcsi(itu_l)
        dxdetac2(itu_l,1) = ter_dxdeta(itu_l)
        dydetac2(itu_l,1) = ter_dydeta(itu_l)
      endif

    endassociate
    !
  end subroutine bcte_met

  subroutine metric_i(self, coeff_deriv1, ep_ord_change, lmax_tag)
    !=======================================
    !Evaluate metric terms along i direction
    !=======================================
    !
    class(field_object), intent(inout) :: self !< The field.
    real(rkind), dimension(1:4,4), intent(in) :: coeff_deriv1
    integer, dimension(0:,0:,0:,1:) , intent(in) :: ep_ord_change
    integer, dimension(1-self%grid%ng:), intent(in) :: lmax_tag
    integer :: i,j,l,m,n,im,istart,iend,imin,imax,lmax,lmaxi
    real(rkind), dimension(2,1-self%grid%ng:self%nx+self%grid%ng) :: u,v,r
    real(rkind), dimension(2,1-self%grid%ng:self%nx+self%grid%ng) :: fhat_met
    real(rkind), dimension(2,self%grid%ng,1-self%grid%ng:self%nx+self%grid%ng) :: uv
    real(rkind), dimension(2) :: fh,ft,uvs

    associate(nx => self%nx, ny => self%ny, metrics_order => self%grid%metrics_order,&
    & xc2 => self%xc2, yc2 => self%yc2, c => self%grid%metrics_cd1, dxdcsic2 => self%dxdcsic2,&
    & dydcsic2 => self%dydcsic2)
      !
      istart = 1
      iend = nx
      lmax = metrics_order/2
      imin = istart - lmax
      imax = iend + lmax
      !
      do j=1,ny ! loop in the j-direction
        !
        !Sweep in the i-direction
        do i=imin,imax
          u(1,i) = xc2(i,j)
          u(2,i) = yc2(i,j)
          v(1,i) = 1._rkind
          v(2,i) = 1._rkind
          r(1,i) = 1._rkind
          r(2,i) = 1._rkind
        enddo
        !
        !Compute u*v averages
        do i=imin,iend
          do l=1,lmax
            do m=1,2
              uv(m,l,i) = (u(m,i)+u(m,i+l))*(v(m,i)+v(m,i+l))*(r(m,i)+r(m,i+l))
            enddo
          enddo
        enddo
        !
        !Evaluation of numerical fluxes
        do i=istart-1,iend ! i is the index of intermediate nodes
          ft = 0._rkind
          !
          lmaxi = max(lmax + ep_ord_change(i,j,1,1), 1 ) ! consider k=1 for order reduction
          if (j==1) lmaxi = lmax_tag(i) ! overwrite lmax for j=1
          !
          do l=1,lmaxi
            uvs = 0._rkind
            do n=0,l-1
              uvs = uvs + uv(:,l,i-n)
            enddo
            ft = ft + coeff_deriv1(l,lmaxi)*uvs
          enddo
          ft = 0.25_rkind*ft
          do m=1,2
            fh(m) = ft(m)
          enddo
          fhat_met(:,i) = fh(:)
        enddo ! end of loop on the cell faces
        !
        !Evaluation of f(u)_csi
        do i=istart,iend ! loop on the inner nodes
          im = i - 1
          dxdcsic2(i,j) = fhat_met(1,i)-fhat_met(1,im)
          dydcsic2(i,j) = fhat_met(2,i)-fhat_met(2,im)
        enddo
      enddo ! end of j-loop
      !
    endassociate
  end subroutine metric_i

  subroutine metric_j(self, coeff_deriv1, ep_ord_change)
    !=======================================
    !Evaluate metric terms along j direction
    !=======================================
    class(field_object), intent(inout) :: self !< The field.
    real(rkind), dimension(1:4,4), intent(in) :: coeff_deriv1
    integer, dimension(0:,0:,0:,1:), intent(in) :: ep_ord_change
    !
    integer :: i,j,jm,l,m,n,jstart,jend,jmin,jmax,lmax,lmaxj
    real(rkind), dimension(2,1-self%grid%ng:self%ny+self%grid%ng) :: u,v,r
    real(rkind), dimension(2,1-self%grid%ng:self%ny+self%grid%ng) :: fhat_met
    real(rkind), dimension(2,self%grid%ng,1-self%grid%ng:self%ny+self%grid%ng) :: uv
    real(rkind), dimension(2) :: fh,ft,uvs
    !
    associate(nx => self%nx, ny => self%ny, metrics_order => self%grid%metrics_order,&
    & xc2 => self%xc2, yc2 => self%yc2, c => self%grid%metrics_cd1, dxdcsic2 => self%dxdcsic2,&
    & dydcsic2 => self%dydcsic2, dxdetac2 => self%dxdetac2, dydetac2 => self%dydetac2,&
    & wall_tag => self%wall_tag)

      jstart = 1
      jend = ny
      lmax = metrics_order/2
      jmin = jstart - lmax
      jmax = jend + lmax
      !
      do i=1,nx ! loop in the i-direction
        !
        !Sweep in the i-direction
        do j=jmin,jmax
          u(1,j) = xc2(i,j)
          u(2,j) = yc2(i,j)
          v(1,j) = 1._rkind
          v(2,j) = 1._rkind
          r(1,j) = 1._rkind
          r(2,j) = 1._rkind
        enddo
        !
        !Three-points averages
        do j=jmin,jend
          do l=1,lmax
            do m=1,2
              uv(m,l,j) = (u(m,j)+u(m,j+l))*(v(m,j)+v(m,j+l))*(r(m,j)+r(m,j+l))
            enddo
          enddo
        enddo
        !
        !Evaluation of numerical fluxes
        do j=jstart-1,jend ! i is the index of intermediate nodes
          ft = 0._rkind
          !
          lmaxj = max(lmax + ep_ord_change(i,j,1,2), 1 ) ! consider k=1 for order reduction
          if (wall_tag(i)<1.and.j<=lmax) lmaxj = max(min(j,lmax),1) ! only for wall (do not reduce wake)
          !
          do l=1,lmaxj
            uvs = 0._rkind
            do n=0,l-1
              uvs = uvs + uv(:,l,j-n)
            enddo
            ft = ft + coeff_deriv1(l,lmaxj)*uvs
          enddo
          ft = 0.25_rkind*ft
          do m=1,2
            fh(m) = ft(m)
          enddo
          fhat_met(:,j) = fh(:)
        enddo ! end of loop on the cell faces
        !
        !Evaluation of g(u)_eta
        do j=jstart,jend ! loop on the inner nodes
          jm = j - 1
          dxdetac2(i,j) = fhat_met(1,j)-fhat_met(1,jm)
          dydetac2(i,j) = fhat_met(2,j)-fhat_met(2,jm)
        enddo
      enddo ! end of i-loop

    endassociate

  end subroutine metric_j
  !

  subroutine read_field_serial(self)
    class(field_object), intent(inout) :: self !< The field.
    character(4) :: chx,chy,chz

    if (self%masterproc) write(*,*) 'Reading rst0_XXX_XXX_XXX.bin'
    1004 format(I4.4)
    write(chx,1004) self%ncoords(1)
    write(chy,1004) self%ncoords(2)
    write(chz,1004) self%ncoords(3)

    open (11,file='rst0_'//chx//'_'//chy//'_'//chz//'.bin',form='unformatted')
    read(11) self%w(1:self%nx,1:self%ny,1:self%nz,1:self%nv)
    close(11)
  endsubroutine read_field_serial

  subroutine read_field(self)
    class(field_object), intent(inout) :: self !< The field.
    !Writing rst.bin
    integer :: mpi_io_file
    integer :: filetype
    integer, dimension(3) :: sizes ! Dimensions of the total grid
    integer, dimension(3) :: subsizes ! Dimensions of grid local to a procs
    integer, dimension(3) :: starts ! Starting coordinates
    integer :: size_real
    integer :: ntot
    integer (kind=mpi_offset_kind) :: offset
    integer :: l,m
    integer, dimension(mpi_status_size) :: istatus
    character(len=256) :: oldname, newname
    !
    if (self%masterproc) write(*,*) 'Reading rst.bin'
    !
    sizes(1) = self%nblocks(1)*self%nx
    sizes(2) = self%nblocks(2)*self%ny
    sizes(3) = self%nblocks(3)*self%nz
    subsizes(1) = self%nx
    subsizes(2) = self%ny
    subsizes(3) = self%nz
    starts(1) = 0 + self%ncoords(1)*subsizes(1)
    starts(2) = 0 + self%ncoords(2)*subsizes(2)
    starts(3) = 0 + self%ncoords(3)*subsizes(3)
    ntot = self%nx*self%ny*self%nz
    !
    call mpi_type_create_subarray(3,sizes,subsizes,starts,mpi_order_fortran,mpi_prec,filetype,self%mpi_err)
    call mpi_type_commit(filetype,self%mpi_err)
    call mpi_file_open(self%mp_cart,'rst.bin',mpi_mode_rdonly,mpi_info_null,mpi_io_file,self%mpi_err)
    offset = 0
    do l=1,self%nv
      call mpi_file_set_view(mpi_io_file,offset,mpi_prec,filetype,"native",mpi_info_null,self%mpi_err)
      call mpi_file_read_all(mpi_io_file,self%w(1:self%nx,1:self%ny,1:self%nz,l),ntot,mpi_prec,istatus,self%mpi_err)
      call mpi_type_size(mpi_prec,size_real,self%mpi_err)
      do m=1,self%nblocks(1)*self%nblocks(2)*self%nblocks(3)
        offset = offset+size_real*ntot
      enddo
    enddo
    !
    call mpi_file_close(mpi_io_file,self%mpi_err)
    call mpi_type_free(filetype,self%mpi_err)
    !
    call mpi_barrier(mpi_comm_world, self%mpi_err)
    if (self%masterproc) then
      oldname = c_char_"rst.bin"//c_null_char
      newname = c_char_"rst.bak"//c_null_char
      self%mpi_err = rename_wrapper(oldname, newname)
      if (self%mpi_err /= 0) write(error_unit,*) "Warning! Cannot rename file rst.bin to rst.bak"
    endif
    call mpi_barrier(mpi_comm_world, self%mpi_err)
  endsubroutine read_field

  subroutine write_field_serial(self)
    class(field_object), intent(inout) :: self !< The field.
    character(4) :: chx,chy,chz

    if (self%masterproc) write(*,*) 'Writing rst1_XXX_XXX_XXX.bin'
    1004 format(I4.4)
    write(chx,1004) self%ncoords(1)
    write(chy,1004) self%ncoords(2)
    write(chz,1004) self%ncoords(3)

    open (11,file='rst1_'//chx//'_'//chy//'_'//chz//'.bin',form='unformatted')
    write(11) self%w(1:self%nx,1:self%ny,1:self%nz,1:self%nv)
    close(11)

  endsubroutine write_field_serial

  subroutine write_field(self)
    class(field_object), intent(inout) :: self !< The field.
    !Writing rst.bin and finaltime.dat
    integer :: mpi_io_file
    integer :: filetype
    integer, dimension(3) :: sizes ! Dimensions of the total grid
    integer, dimension(3) :: subsizes ! Dimensions of grid local to a procs
    integer, dimension(3) :: starts ! Starting coordinates
    integer :: size_real
    integer :: ntot
    integer (kind=mpi_offset_kind) :: offset
    integer :: l,m
    integer, dimension(mpi_status_size) :: istatus
    !
    if (self%masterproc) write(*,*) 'Writing rst.bin'
    !
    sizes(1) = self%nblocks(1)*self%nx
    sizes(2) = self%nblocks(2)*self%ny
    sizes(3) = self%nblocks(3)*self%nz
    subsizes(1) = self%nx
    subsizes(2) = self%ny
    subsizes(3) = self%nz
    starts(1) = 0 + self%ncoords(1)*subsizes(1)
    starts(2) = 0 + self%ncoords(2)*subsizes(2)
    starts(3) = 0 + self%ncoords(3)*subsizes(3)
    ntot = self%nx*self%ny*self%nz
    !
    call mpi_type_create_subarray(3,sizes,subsizes,starts,mpi_order_fortran,mpi_prec,filetype,self%mpi_err)
    call mpi_type_commit(filetype,self%mpi_err)
    call mpi_file_open(self%mp_cart,'rst.bin',mpi_mode_create+mpi_mode_wronly,mpi_info_null,mpi_io_file,self%mpi_err)
    offset = 0
    do l=1,self%nv
      call mpi_file_set_view(mpi_io_file,offset,mpi_prec,filetype,"native",mpi_info_null,self%mpi_err)
      call mpi_file_write_all(mpi_io_file,self%w(1:self%nx,1:self%ny,1:self%nz,l),ntot,mpi_prec,istatus,self%mpi_err)
      call mpi_type_size(mpi_prec,size_real,self%mpi_err)
      do m=1,self%nblocks(1)*self%nblocks(2)*self%nblocks(3)
        offset = offset+size_real*ntot
      enddo
    enddo
    !
    call mpi_file_close(mpi_io_file,self%mpi_err)
    call mpi_type_free(filetype,self%mpi_err)
    !
  endsubroutine write_field

  subroutine write_plot3d(self, mach, reynolds, time, istore, w_io, plot3dgrid, plot3dfield, file_prefix, w_aux_io)
    class(field_object), intent(inout) :: self !< The field.
    !Writing field (MPI I/O)
    real(rkind), optional :: mach, reynolds, time
    real(rkind) :: mach_, reynolds_, time_
    integer :: istore
    character(*), optional :: file_prefix
    character(32) :: file_prefix_
    logical, optional :: plot3dgrid
    logical, optional :: plot3dfield
    logical :: plot3dgrid_, plot3dfield_
    real(rkind), dimension(:,:,:,:), optional :: w_io
    real(rkind), dimension(:,:,:,:), allocatable :: w_io_
    real(rkind), dimension(:,:,:,:) :: w_aux_io
    integer :: nv_io
    integer :: i,j,k,l,m
    integer :: mpi_io_file
    integer :: filetype
    integer,dimension(3) :: sizes ! Dimensions of the total grid
    integer,dimension(3) :: subsizes ! Dimensions of grid local to a procs
    integer,dimension(3) :: starts ! Starting coordinates
    integer :: size_real,size_integer
    integer :: ntot
    integer (kind=mpi_offset_kind) :: offset
    character(4) :: chstore
    integer, dimension(mpi_status_size) :: istatus
    integer :: istat

    write(chstore(1:4), '(I4.4)') istore

    mach_ = 1._rkind ; if(present(mach)) mach_ = mach
    reynolds_ = 1._rkind ; if(present(reynolds)) reynolds_ = reynolds
    time_ = 1._rkind ; if(present(time)) time_ = time
    plot3dfield_ = .true. ; if(present(plot3dfield)) plot3dfield_ = plot3dfield
    plot3dgrid_ = .false. ; if(present(plot3dgrid)) plot3dgrid_ = plot3dgrid

    if(plot3dgrid_) then
      istat = create_folder("FIELDS")

      if (self%masterproc) write(*,*) 'Storing grid in plot3D format'
      !
      if (self%masterproc) then
        open(unit=123, file='FIELDS/plot3dgrid.xyz', access="stream", form="unformatted")
        write(123) self%grid%nxmax,self%grid%nymax,self%grid%nzmax
        flush(123)
        close(123)
      endif
      call mpi_barrier(mpi_comm_world,self%mpi_err)
      !
      sizes(1) = self%nblocks(1)*self%nx
      sizes(2) = self%nblocks(2)*self%ny
      sizes(3) = self%nblocks(3)*self%nz
      subsizes(1) = self%nx
      subsizes(2) = self%ny
      subsizes(3) = self%nz
      starts(1) = 0 + self%ncoords(1)*subsizes(1)
      starts(2) = 0 + self%ncoords(2)*subsizes(2)
      starts(3) = 0 + self%ncoords(3)*subsizes(3)
      ntot = self%nx*self%ny*self%nz
      !
      call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,mpi_prec,filetype,self%mpi_err)
      call MPI_TYPE_COMMIT(filetype,self%mpi_err)
      call MPI_FILE_OPEN(self%mp_cart,'FIELDS/plot3dgrid.xyz',MPI_MODE_RDWR,MPI_INFO_NULL,mpi_io_file,self%mpi_err)
      size_integer = storage_size(1_ikind)/8
      offset = 3*size_integer
      do l=1,3
        call MPI_FILE_SET_VIEW(mpi_io_file,offset,mpi_prec,filetype,"native",MPI_INFO_NULL,self%mpi_err)
        call MPI_FILE_WRITE_ALL(mpi_io_file,self%grid3d(l,:,:,:),ntot,mpi_prec,istatus,self%mpi_err)
        call MPI_TYPE_SIZE(mpi_prec,size_real,self%mpi_err)
        do m=1,self%nblocks(1)*self%nblocks(2)*self%nblocks(3)
          offset = offset+size_real*ntot
        enddo
      enddo
      !
      call MPI_FILE_CLOSE(mpi_io_file,self%mpi_err)
      call MPI_TYPE_FREE(filetype,self%mpi_err)

    endif

    if (plot3dfield_) then
      if (present(w_io)) then
        nv_io = 1
        w_io_ = w_io
        file_prefix_ = file_prefix
      else
        nv_io = self%nv
        !w_io_ = self%w(1:self%nx,1:self%ny,1:self%nz,:)
        allocate(w_io_(1:self%nx,1:self%ny,1:self%nz,1:self%nv))
        w_io_(:,:,:,1) = w_aux_io(:,:,:,1)
        w_io_(:,:,:,2) = w_aux_io(:,:,:,2)
        w_io_(:,:,:,3) = w_aux_io(:,:,:,3)
        w_io_(:,:,:,4) = w_aux_io(:,:,:,4)
        w_io_(:,:,:,5) = w_aux_io(:,:,:,6)
        file_prefix_ = "flow_"
      endif

      size_real = storage_size(1.0_rkind)/8
      size_integer = storage_size(1_ikind)/8
      !
      if (self%masterproc) then
        open(unit=123, file="FIELDS/"//trim(file_prefix_)//chstore//'.q', access="stream", form="unformatted")
        if(present(w_io)) then
          write(123) self%grid%nxmax, self%grid%nymax, self%grid%nzmax, nv_io
        else
          write(123) self%grid%nxmax, self%grid%nymax, self%grid%nzmax
          write(123) mach_, 0._rkind, reynolds_, time_
        endif
        flush(123)
        close(123)
      endif
      call mpi_barrier(self%mp_cart,self%mpi_err)
      !
      sizes(1) = self%nblocks(1)*self%nx
      sizes(2) = self%nblocks(2)*self%ny
      sizes(3) = self%nblocks(3)*self%nz
      subsizes(1) = self%nx
      subsizes(2) = self%ny
      subsizes(3) = self%nz
      starts(1) = 0 + self%ncoords(1)*subsizes(1)
      starts(2) = 0 + self%ncoords(2)*subsizes(2)
      starts(3) = 0 + self%ncoords(3)*subsizes(3)
      ntot = self%nx*self%ny*self%nz
      !
      call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,mpi_prec,filetype,self%mpi_err)
      call MPI_TYPE_COMMIT(filetype,self%mpi_err)
      call MPI_FILE_OPEN(self%mp_cart,"FIELDS/"//trim(file_prefix_)//chstore//'.q',MPI_MODE_RDWR,&
      &MPI_INFO_NULL,mpi_io_file,self%mpi_err)
      if(present(w_io)) then
        offset = 4*size_integer
      else
        offset = 3*size_integer+4*size_real
      endif
      do l=1,nv_io
        call MPI_FILE_SET_VIEW(mpi_io_file,offset,mpi_prec,filetype,"native",MPI_INFO_NULL,self%mpi_err)
        call MPI_FILE_WRITE_ALL(mpi_io_file,w_io_(:,:,:,l),ntot,mpi_prec,istatus,self%mpi_err)
        do m=1,self%nblocks(1)*self%nblocks(2)*self%nblocks(3)
          offset = offset+size_real*ntot
        enddo
      enddo
      !
      call MPI_FILE_CLOSE(mpi_io_file,self%mpi_err)
      call MPI_TYPE_FREE(filetype,self%mpi_err)
    endif
  end subroutine write_plot3d

  subroutine read_plot3d(self, istore)
    class(field_object), intent(inout) :: self !< The field.
    !Writing field (MPI I/O)
    integer :: istore
    character(32) :: file_prefix
    integer :: i,j,k,l,m
    integer :: mpi_io_file
    integer :: filetype
    integer,dimension(3) :: sizes ! Dimensions of the total grid
    integer,dimension(3) :: subsizes ! Dimensions of grid local to a procs
    integer,dimension(3) :: starts ! Starting coordinates
    integer :: size_real,size_integer
    integer :: ntot
    integer (kind=mpi_offset_kind) :: offset
    character(4) :: chstore
    integer, dimension(mpi_status_size) :: istatus
    integer :: istat

    write(chstore(1:4), '(I4.4)') istore

    file_prefix = "flow_"
    if (self%masterproc) write(*,*) trim(file_prefix)//chstore//'.q'

    size_real = storage_size(1.0_rkind)/8
    size_integer = storage_size(1_ikind)/8
    !
    sizes(1) = self%nblocks(1)*self%nx
    sizes(2) = self%nblocks(2)*self%ny
    sizes(3) = self%nblocks(3)*self%nz
    subsizes(1) = self%nx
    subsizes(2) = self%ny
    subsizes(3) = self%nz
    starts(1) = 0 + self%ncoords(1)*subsizes(1)
    starts(2) = 0 + self%ncoords(2)*subsizes(2)
    starts(3) = 0 + self%ncoords(3)*subsizes(3)
    ntot = self%nx*self%ny*self%nz
    !
    call mpi_type_create_subarray(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,mpi_prec,filetype,self%mpi_err)
    call mpi_type_commit(filetype,self%mpi_err)
    call mpi_file_open(self%mp_cart,"FIELDS/"//trim(file_prefix)//chstore//'.q',mpi_mode_rdonly,&
    &mpi_info_null,mpi_io_file,self%mpi_err)
    offset = 3*size_integer+4*size_real
    do l=1,self%nv
      call mpi_file_set_view(mpi_io_file,offset,mpi_prec,filetype,"native",mpi_info_null,self%mpi_err)
      call mpi_file_read_all(mpi_io_file,self%w(1:self%nx,1:self%ny,1:self%nz,l),ntot,mpi_prec,istatus,self%mpi_err)
      do m=1,self%nblocks(1)*self%nblocks(2)*self%nblocks(3)
        offset = offset+size_real*ntot
      enddo
    enddo
    !
    call mpi_file_close(mpi_io_file,self%mpi_err)
    call mpi_type_free(filetype,self%mpi_err)
    !
  end subroutine read_plot3d

  subroutine correct_bc(self, ibc)
    class(field_object), intent(inout) :: self !< The field.
    integer, dimension(6), intent(inout) :: ibc
    integer, dimension(6) :: itag
    integer :: iend, jend, kend
    itag = ibc
    ibc = 0
    iend = self%nblocks(1)-1
    jend = self%nblocks(2)-1
    kend = self%nblocks(3)-1

    if (self%ncoords(1)== 0) ibc(1) = itag(1)
    if (self%ncoords(1)==iend) ibc(2) = itag(2)
    if (self%ncoords(2)== 0) ibc(3) = itag(3)
    if (self%ncoords(2)==jend) ibc(4) = itag(4)
    if (self%ncoords(3)== 0) ibc(5) = itag(5)
    if (self%ncoords(3)==kend) ibc(6) = itag(6)
  endsubroutine correct_bc

  subroutine write_vtk(self, time, istore, w_aux_io)
    class(field_object), intent(inout) :: self !< The field.
    real(rkind), optional :: time
    integer :: istore
    real(rkind), dimension(:,:,:,:) :: w_aux_io
    integer :: istat
    istat = create_folder("FIELDS")
    if(self%grid%grid_dim == 1) then
      call write_vtk_c1(self, time, istore, w_aux_io)
    elseif(self%grid%grid_dim == 2) then
      call write_vtk_c2(self, time, istore, w_aux_io)
    endif
  endsubroutine write_vtk

  subroutine write_vtk_c1(self, time, istore, w_aux_io)
    class(field_object), intent(inout) :: self !< The field.
    !Writing field (MPI I/O)
    real(rkind), optional :: time
    real(rkind) :: time_
    integer :: istore
    character(32) :: file_prefix_
    real(rkind), dimension(:,:,:,:), allocatable :: w_io_
    real(rkind), dimension(:,:,:,:) :: w_aux_io
    integer :: nv_io
    integer :: l
    integer :: mpi_io_file
    integer :: filetype
    integer,dimension(3) :: sizes ! Dimensions of the total grid
    integer,dimension(3) :: subsizes ! Dimensions of grid local to a procs
    integer,dimension(3) :: starts ! Starting coordinates
    integer :: size_real
    integer :: ntot
    integer (kind=mpi_offset_kind) :: offset
    integer (kind=mpi_offset_kind) :: offset_x,offset_y,offset_z,delta_offset_w
    character(4) :: chstore
    integer, dimension(mpi_status_size) :: istatus
    integer, parameter :: int64_kind = selected_int_kind(2*range(1))
    integer(int64_kind) :: gridsize_64
    !character(len=16) :: int2str
    !character(len=32) :: int2str_o
    character(len=65536) :: xml_part
    character(len=7) :: vtk_float
    character(len=12), dimension(5) :: names


    write(chstore(1:4), '(I4.4)') istore

    time_ = 1._rkind ; if(present(time)) time_ = time
    if (self%masterproc) print *, 'Storing VTK sol', istore,'at time', time

    call MPI_TYPE_SIZE(mpi_prec,size_real,self%mpi_err)
    if(size_real == 4) then
      vtk_float = "Float32"
    elseif(size_real == 8) then
      vtk_float = "Float64"
    else
      if(self%masterproc) write(*,*) "Error on VTK write! size_real must be either 4 or 8"
      call MPI_ABORT(MPI_COMM_WORLD,self%mpi_err,self%mpi_err)
    endif
    gridsize_64 = int(size_real,int64_kind)*int(self%grid%nxmax,int64_kind)*int(self%grid%nymax,&
    &int64_kind)*int(self%grid%nzmax,int64_kind)

    if(storage_size(gridsize_64) /= 64) then
      if(self%masterproc) write(*,*) "Error on VTK write! Size of int64_kind integers is not 8 bytes!"
      call MPI_ABORT(MPI_COMM_WORLD, self%mpi_err, self%mpi_err)
    endif

    nv_io = self%nv
    file_prefix_ = "flow_"
    !w_io_ = self%w(1:self%nx,1:self%ny,1:self%nz,:)
    allocate(w_io_(1:self%nx,1:self%ny,1:self%nz,1:self%nv))
    w_io_(:,:,:,1) = w_aux_io(:,:,:,1)
    w_io_(:,:,:,2) = w_aux_io(:,:,:,2)
    w_io_(:,:,:,3) = w_aux_io(:,:,:,3)
    w_io_(:,:,:,4) = w_aux_io(:,:,:,4)
    w_io_(:,:,:,5) = w_aux_io(:,:,:,6)

    names= [character(len=12) :: "density", "velocity_x", "velocity_y", "velocity_z", "temperature"]


    if (self%masterproc) then
      offset_x = 0
      offset_y = size_real*self%grid%nxmax + storage_size(gridsize_64)/8
      offset_z = offset_y + size_real*self%grid%nymax + storage_size(gridsize_64)/8
      delta_offset_w = gridsize_64 + storage_size(gridsize_64)/8 ! the second part is because of the header of bytes before data

      open(unit=123, file=trim(file_prefix_)//chstore//'.vtr', access="stream", form="unformatted", status="replace")
      xml_part = ' <?xml version="1.0"?> <VTKFile type="RectilinearGrid" version="1.0" byte_order="L&
      &ittleEndian" header_type="UInt64"> <RectilinearGrid WholeExtent="+1 +' //int2str(self%grid%nxmax)//&
      &' +1 +'//int2str(self%grid%nymax)//' +1 +'//int2str(self%grid%nzmax)//'"> <Piece Extent="+1 +'//&
      &int2str(self%grid%nxmax)// ' +1 +'//int2str(self%grid%nymax)//' +1 +'//int2str(self%grid%nzmax)//&
      &'"> <Coordinates> <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="X" format="appended&
      &" offset="'// int2str_o(offset_x)//'"/> <DataArray type="'//vtk_float//'" NumberOfComponents="1" Nam&
      &e="Y" format="appended" offset="'// int2str_o(offset_y)//'"/> <DataArray type="'//vtk_float//&
      &'" NumberOfComponents="1" Name="Z" format="appended" offset="'// int2str_o(offset_z)//'"/> </&
      &Coordinates> <PointData> '
      offset = offset_z + size_real*self%grid%nzmax + storage_size(gridsize_64)/8
      do l=1,nv_io
        xml_part = trim(adjustl(xml_part)) // ' <DataArray type="'//vtk_float//'" NumberOfComponents&
        &="1" Name="'//trim(names(l))//'" format="appended" offset="'//int2str_o(offset)//'"/>'
        offset = offset + delta_offset_w
      enddo
      xml_part = trim(adjustl(xml_part)) // ' </PointData> </Piece> </RectilinearGrid> <AppendedData encoding="raw"> '
      write(123) trim(adjustl(xml_part))

      write(123) "_"
      !write(123) storage_size(gridsize_64)/8*int(nxmax,int64_kind) , xg(1:nxmax)
      !write(123) storage_size(gridsize_64)/8*int(nymax,int64_kind) , yg(1:nymax)
      !write(123) storage_size(gridsize_64)/8*int(nzmax,int64_kind) , zg(1:nzmax)
      write(123) size_real*int(self%grid%nxmax,int64_kind) , self%grid%xg(1:self%grid%nxmax)
      write(123) size_real*int(self%grid%nymax,int64_kind) , self%grid%yg(1:self%grid%nymax)
      write(123) size_real*int(self%grid%nzmax,int64_kind) , self%grid%zg(1:self%grid%nzmax)
      flush(123)
      close(123)
    endif
    !
    sizes(1) = self%nblocks(1)*self%nx
    sizes(2) = self%nblocks(2)*self%ny
    sizes(3) = self%nblocks(3)*self%nz
    subsizes(1) = self%nx
    subsizes(2) = self%ny
    subsizes(3) = self%nz
    starts(1) = 0 + self%ncoords(1)*subsizes(1)
    starts(2) = 0 + self%ncoords(2)*subsizes(2)
    starts(3) = 0 + self%ncoords(3)*subsizes(3)
    ntot = self%nx*self%ny*self%nz
    !
    call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,mpi_prec,filetype,self%mpi_err)
    call MPI_TYPE_COMMIT(filetype,self%mpi_err)
    !
    do l=1,nv_io
      call MPI_BARRIER(self%mp_cart,self%mpi_err)
      if (self%masterproc) then
        open(unit=123, file=trim(file_prefix_)//chstore//'.vtr', access="stream", form="unformatted", position="append")
        write(123) gridsize_64
        flush(123)
        close(123)
      endif
      call MPI_BARRIER(self%mp_cart, self%mpi_err)
      call MPI_FILE_OPEN(self%mp_cart,trim(file_prefix_)//chstore//'.vtr',MPI_MODE_RDWR,MPI_INFO_NULL,mpi_io_file,self%mpi_err)
      call MPI_FILE_GET_SIZE(mpi_io_file, offset, self%mpi_err)
      call MPI_BARRIER(self%mp_cart, self%mpi_err)
      call MPI_FILE_SET_VIEW(mpi_io_file,offset,mpi_prec,filetype,"native",MPI_INFO_NULL,self%mpi_err)
      call MPI_FILE_WRITE_ALL(mpi_io_file,w_io_(1:self%nx,1:self%ny,1:self%nz,l),ntot,mpi_prec,istatus,self%mpi_err)
      call MPI_FILE_CLOSE(mpi_io_file,self%mpi_err)
    enddo

    call MPI_TYPE_FREE(filetype,self%mpi_err)
    if (self%masterproc) then
      open(unit=123, file=trim(file_prefix_)//chstore//'.vtr', access="stream", position="append", form="unformatted")
      write(123) ' </AppendedData> </VTKFile>'
      close(123)
    endif
  end subroutine write_vtk_c1

  subroutine write_vtk_c2(self, time, istore, w_aux_io)
    class(field_object), intent(inout) :: self !< The field.
    !Writing field (MPI I/O)
    real(rkind), optional :: time
    real(rkind) :: time_
    integer :: istore
    character(32) :: file_prefix_
    real(rkind), dimension(:,:,:,:), allocatable :: w_io_
    real(rkind), dimension(:,:,:,:) :: w_aux_io
    integer :: nv_io
    integer :: l
    integer :: mpi_io_file
    integer :: filetype
    integer,dimension(3) :: sizes ! Dimensions of the total grid
    integer,dimension(3) :: subsizes ! Dimensions of grid local to a procs
    integer,dimension(3) :: starts ! Starting coordinates
    integer :: size_real
    integer :: ntot
    integer (kind=mpi_offset_kind) :: offset
    integer (kind=mpi_offset_kind) :: offset_x,offset_y,offset_z,delta_offset_w
    character(4) :: chstore
    integer, dimension(mpi_status_size) :: istatus
    integer, parameter :: int64_kind = selected_int_kind(2*range(1))
    integer(int64_kind) :: gridsize_64
    !character(len=16) :: int2str
    !character(len=32) :: int2str_o
    character(len=65536) :: xml_part
    character(len=7) :: vtk_float
    character(len=12), dimension(5) :: names
    integer :: mpi_triplet

    write(chstore(1:4), '(I4.4)') istore

    time_ = 1._rkind ; if(present(time)) time_ = time
    if (self%masterproc) print *, 'Storing VTK sol', istore,'at time', time

    call MPI_TYPE_SIZE(mpi_prec,size_real,self%mpi_err)
    if(size_real == 4) then
      vtk_float = "Float32"
    elseif(size_real == 8) then
      vtk_float = "Float64"
    else
      if(self%masterproc) write(*,*) "Error on VTK write! size_real must be either 4 or 8"
      call MPI_ABORT(MPI_COMM_WORLD,self%mpi_err,self%mpi_err)
    endif
    gridsize_64 = int(size_real,int64_kind)*int(self%grid%nxmax,int64_kind)*int(self%grid%nymax,&
    &int64_kind)*int(self%grid%nzmax,int64_kind)

    if(storage_size(gridsize_64) /= 64) then
      if(self%masterproc) write(*,*) "Error on VTK write! Size of int64_kind integers is not 8 bytes!"
      call MPI_ABORT(MPI_COMM_WORLD, self%mpi_err, self%mpi_err)
    endif

    nv_io = self%nv
    file_prefix_ = "flow_"
    !w_io_ = self%w(1:self%nx,1:self%ny,1:self%nz,:)
    allocate(w_io_(1:self%nx,1:self%ny,1:self%nz,1:self%nv))
    w_io_(:,:,:,1) = w_aux_io(:,:,:,1)
    w_io_(:,:,:,2) = w_aux_io(:,:,:,2)
    w_io_(:,:,:,3) = w_aux_io(:,:,:,3)
    w_io_(:,:,:,4) = w_aux_io(:,:,:,4)
    w_io_(:,:,:,5) = w_aux_io(:,:,:,6)

    names= [character(len=12) :: "density", "velocity_x", "velocity_y", "velocity_z", "temperature"]

    if (self%masterproc) then
      offset_x = 0
      offset_y = gridsize_64 + storage_size(gridsize_64)/8
      offset_z = offset_y + gridsize_64 !+ storage_size(gridsize_64)/8
      delta_offset_w = gridsize_64 + storage_size(gridsize_64)/8 ! the second part is because of the header of bytes before data

      open(unit=123, file="FIELDS/"//trim(file_prefix_)//chstore//'.vts', access="stream", form="unformatted", status="replace")
      xml_part = ' <?xml version="1.0"?> <VTKFile type="StructuredGrid" version="1.0" byte_order="Li&
      &ttleEndian" header_type="UInt64"> <StructuredGrid WholeExtent="+1 +' //int2str(self%grid%nxmax)//' +&
      &1 +'//int2str(self%grid%nymax)//' +1 +'//int2str(self%grid%nzmax)//'"> <Piece Extent="+1 +'//&
      &int2str(self%grid%nxmax)// ' +1 +'//int2str(self%grid%nymax)//' +1 +'//int2str(self%grid%nzmax)//&
      &'"> <Points> <DataArray type="'//vtk_float//'" NumberOfComponents="3" Name="XYZ" format="appended" o&
      &ffset="'// int2str_o(offset_x)//'"/> </Points> <PointData> '
      offset = offset_z + gridsize_64 !+ storage_size(gridsize_64)/8
      do l=1,nv_io
        xml_part = trim(adjustl(xml_part)) // ' <DataArray type="'//vtk_float//'" NumberOfComponents&
        &="1" Name="'//trim(names(l))//'" format="appended" offset="'//int2str_o(offset)//'"/>'
        offset = offset + delta_offset_w
      enddo
      xml_part = trim(adjustl(xml_part)) // ' </PointData> </Piece> </StructuredGrid> <AppendedData encoding="raw"> '
      write(123) trim(adjustl(xml_part))

      write(123) "_"
      flush(123)
      close(123)
    endif
    !
    sizes(1) = self%nblocks(1)*self%nx
    sizes(2) = self%nblocks(2)*self%ny
    sizes(3) = self%nblocks(3)*self%nz
    subsizes(1) = self%nx
    subsizes(2) = self%ny
    subsizes(3) = self%nz
    starts(1) = 0 + self%ncoords(1)*subsizes(1)
    starts(2) = 0 + self%ncoords(2)*subsizes(2)
    starts(3) = 0 + self%ncoords(3)*subsizes(3)
    ntot = self%nx*self%ny*self%nz
    !
    call MPI_TYPE_CONTIGUOUS(3, mpi_prec, mpi_triplet, self%mpi_err)
    call MPI_TYPE_COMMIT(mpi_triplet,self%mpi_err)
    call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,mpi_triplet,filetype,self%mpi_err)
    call MPI_TYPE_COMMIT(filetype,self%mpi_err)

    call MPI_BARRIER(self%mp_cart,self%mpi_err)
    if (self%masterproc) then
      open(unit=123, file="FIELDS/"//trim(file_prefix_)//chstore//'.vts', access="stream", form="unformatted", position="append")
      write(123) gridsize_64*3
      flush(123)
      close(123)
    endif
    call MPI_BARRIER(self%mp_cart, self%mpi_err)
    call MPI_FILE_OPEN(self%mp_cart,"FIELDS/"//trim(file_prefix_)//chstore//'.vts',MPI_MODE_RDWR,&
    &MPI_INFO_NULL,mpi_io_file,self%mpi_err)
    call MPI_FILE_GET_SIZE(mpi_io_file, offset, self%mpi_err)
    call MPI_BARRIER(self%mp_cart, self%mpi_err)
    call MPI_FILE_SET_VIEW(mpi_io_file,offset,mpi_prec,filetype,"native",MPI_INFO_NULL,self%mpi_err)
    call MPI_FILE_WRITE_ALL(mpi_io_file,self%grid3d,ntot,mpi_triplet,istatus,self%mpi_err)
    call MPI_FILE_CLOSE(mpi_io_file,self%mpi_err)
    call MPI_TYPE_FREE(filetype,self%mpi_err)

    call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,mpi_prec,filetype,self%mpi_err)
    call MPI_TYPE_COMMIT(filetype,self%mpi_err)
    !
    call MPI_BARRIER(self%mp_cart,self%mpi_err)
    do l=1,nv_io
      call MPI_BARRIER(self%mp_cart,self%mpi_err)
      if (self%masterproc) then
        open(unit=123, file="FIELDS/"//trim(file_prefix_)//chstore//'.vts', access="stream", form="unformatted", position="append")
        write(123) gridsize_64
        flush(123)
        close(123)
      endif
      call MPI_BARRIER(self%mp_cart, self%mpi_err)
      call MPI_FILE_OPEN(self%mp_cart,"FIELDS/"//trim(file_prefix_)//chstore//'.vts',MPI_MODE_RDWR,&
      &MPI_INFO_NULL,mpi_io_file,self%mpi_err)
      call MPI_FILE_GET_SIZE(mpi_io_file, offset, self%mpi_err)
      call MPI_BARRIER(self%mp_cart, self%mpi_err)
      call MPI_FILE_SET_VIEW(mpi_io_file,offset,mpi_prec,filetype,"native",MPI_INFO_NULL,self%mpi_err)
      call MPI_FILE_WRITE_ALL(mpi_io_file,w_io_(1:self%nx,1:self%ny,1:self%nz,l),ntot,mpi_prec,istatus,self%mpi_err)
      call MPI_FILE_CLOSE(mpi_io_file,self%mpi_err)
    enddo

    call MPI_TYPE_FREE(filetype,self%mpi_err)
    if (self%masterproc) then
      open(unit=123, file="FIELDS/"//trim(file_prefix_)//chstore//'.vts', access="stream", position="append", form="unformatted")
      write(123) ' </AppendedData> </VTKFile>'
      close(123)
    endif
  end subroutine write_vtk_c2

endmodule streams_field_object

