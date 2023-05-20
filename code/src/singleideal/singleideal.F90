module streams_equation_singleideal_object
! < STREAmS, Euler equations backend-independent variables and routines.
!
  use streams_field_object
  use streams_grid_object
  use streams_parameters
  use cfgio_mod, only: cfg_t, parse_cfg
  use, intrinsic :: iso_fortran_env, only : error_unit
  use crandom_f_mod
  use MPI
  use ISO_C_BINDING
!
!
!
  implicit none
  private
  public :: equation_singleideal_object
  public :: VISC_POWER, VISC_SUTHERLAND, VISC_NO
  public :: RK_WRAY, RK_JAMESON, RK_SHU
  public :: slicexy_aux, slicexz_aux, sliceyz_aux
  public :: ibm_MAX_PARBC
!
! public constants
  integer(ikind), parameter :: VISC_NO            = 0_ikind
  integer(ikind), parameter :: VISC_POWER         = 1_ikind
  integer(ikind), parameter :: VISC_SUTHERLAND    = 2_ikind
!
  integer(ikind), parameter :: RK_WRAY            = 1_ikind
  integer(ikind), parameter :: RK_JAMESON         = 2_ikind
  integer(ikind), parameter :: RK_SHU             = 3_ikind
!
  integer(ikind), parameter :: ibm_MAX_PARBC      = 3_ikind
!
  real(rkind), dimension(:,:,:,:), allocatable :: slicexy_aux, slicexz_aux, sliceyz_aux
!
  type :: equation_singleideal_object
    type(grid_object)         :: grid
    type(field_object)        :: field
    type(cfg_t)               :: cfg
    type(cfg_t)               :: flow_params_cfg
    type(cfg_t)               :: field_info_cfg
!
    integer(ikind)            :: nx, ny, nz        ! replica from field
    integer(ikind)            :: ng                ! replica from grid
    integer(ikind)            :: mpi_err
    integer(ikind)            :: myrank
    integer(ikind)            :: nprocs
    logical                   :: masterproc
!
    integer(ikind)            :: nv
    integer(ikind)            :: nv_aux
    integer(ikind)            :: nv_stat
    integer(ikind)            :: nv_stat_3d
    integer(ikind)            :: enable_stat_3d
!
    real(rkind)               :: dt
    real(rkind)               :: cfl
    real(rkind)               :: time0, time
    integer                   :: num_iter
    integer                   :: icyc0, icyc
!
    integer(ikind)            :: visc_order, conservative_viscous
    integer(ikind)            :: ep_order, nkeep
    integer(ikind)            :: weno_scheme, weno_version, flux_splitting
    real(rkind)               :: sensor_threshold
    real(rkind)               :: xshock_imp, shock_angle, tanhfacs
    integer                   :: rand_type
    integer                   :: eul_imin, eul_imax, eul_jmin, eul_jmax, eul_kmin, eul_kmax
    integer, dimension(6)     :: force_zero_flux
    integer(ikind), dimension(:), allocatable :: supported_orders
!
    real(rkind), dimension(1:4,4) :: coeff_deriv1
    real(rkind), dimension(0:4,4) :: coeff_deriv2
!
    integer(ikind)            :: rk_type
    integer(ikind)            :: nrk
    real(rkind), allocatable, dimension(:) :: rhork,gamrk,alprk
    real(rkind), allocatable, dimension(:) :: ark,brk,crk
!
    integer                   :: iter_dt_recompute, print_control
!
    real(rkind)               :: dtsave, dtsave_restart, dtstat, dtslice
    real(rkind)               :: time_from_last_rst, time_from_last_write, time_from_last_stat, time_from_last_slice
    integer                   :: restart_type
    integer                   :: io_type
    integer                   :: istore
    integer                   :: itav
    integer                   :: enable_plot3d, enable_vtk
    real(rkind)               :: residual_rhou
    real(rkind)               :: entropy
!
    integer                   :: flow_init
    real(rkind)               :: Mach, Reynolds, Prandtl, theta_wall
    real(rkind)               :: Reynolds_friction
    real(rkind)               :: rgas0
    real(rkind)               :: rho0
    real(rkind)               :: t0
    real(rkind)               :: l0
    real(rkind)               :: u0, p0, mu0
    real(rkind)               :: gam, gm, gm1, rfac
    real(rkind)               :: T_wall, T_bulk_target, T_recovery
    real(rkind)               :: powerlaw_vtexp, T_ref_dim, sutherland_S
    integer                   :: visc_model
    real(rkind)               :: dpdx,rhobulk,ubulk,tbulk,volchan   ! for channel flow
    logical                   :: channel_case, bl_case, bl_laminar
    real(rkind)               :: x_recyc
    logical                   :: recyc
    integer                   :: i_recyc, ib_recyc
    real(rkind), dimension(:), allocatable :: winf, winf_past_shock
!
    real(rkind), dimension(:,:,:,:), allocatable :: w_aux
    real(rkind), dimension(:,:,:), allocatable   :: wmean
    integer, dimension(:,:,:), allocatable       :: fluid_mask
    integer, dimension(:,:,:,:), allocatable     :: ep_ord_change
    integer, dimension(:,:,:), allocatable       :: ep_ord_change_x
    integer :: correct_bound_ord
!
    real(rkind), dimension(:,:,:), allocatable   :: w_stat
    real(rkind), dimension(:,:,:,:), allocatable :: w_stat_3d
!
    integer, dimension(6)     :: bctags, bctags_nr
!
    real(rkind), dimension(:), allocatable :: deltavec, deltavvec, cfvec
    real(rkind) :: betarecyc
    real(rkind), dimension(:), allocatable :: yplus_inflow
    real(rkind), dimension(:), allocatable :: eta_inflow
    real(rkind), dimension(:), allocatable :: yplus_recyc
    real(rkind), dimension(:), allocatable :: eta_recyc
    integer, dimension(:), allocatable :: map_j_inn, map_j_out
    real(rkind), dimension(:), allocatable :: weta_inflow
!
    integer :: calorically_perfect, indx_cp_l, indx_cp_r
    real(rkind), allocatable, dimension(:) :: cp_coeff, cv_coeff
!   
    integer(ikind), dimension(:), allocatable :: igslice, jgslice, kgslice
    integer(ikind), dimension(:), allocatable :: islice, jslice, kslice
    integer(ikind) :: num_probe, num_probe_tot
    real(rkind), dimension(:,:), allocatable :: probe_coord, w_aux_probe
    real(rkind), dimension(:,:,:,:), allocatable :: probe_coeff
    integer, dimension(:,:), allocatable :: ijk_probe
!
    integer :: debug_memory = 0
    integer :: mode_async = 0
!   
!
!
!   
!
    logical :: time_is_freezed=.false.
!
!   
!
  contains
!   public methods
    procedure, pass(self) :: initialize
    procedure, pass(self) :: read_input
    procedure, pass(self) :: runge_kutta_initialize
    procedure, pass(self) :: initial_conditions
    procedure, pass(self) :: read_field_info
    procedure, pass(self) :: write_field_info
    procedure, pass(self) :: set_oblique_shock
    procedure, pass(self) :: set_fluid_prop
    procedure, pass(self) :: set_flow_params
    procedure, pass(self) :: set_bl_prop
    procedure, pass(self) :: set_chan_prop
    procedure, pass(self) :: init_channel
    procedure, pass(self) :: init_wind_tunnel
    procedure, pass(self) :: init_bl
    procedure, pass(self) :: init_bl_lam
    procedure, pass(self) :: init_tgv
    procedure, pass(self) :: alloc
    procedure, pass(self) :: bc_preproc
    procedure, pass(self) :: compute_stats
    procedure, pass(self) :: write_stats
    procedure, pass(self) :: write_stats_serial
    procedure, pass(self) :: read_stats
    procedure, pass(self) :: read_stats_serial
    procedure, pass(self) :: compute_stats_3d
    procedure, pass(self) :: write_stats_3d
    procedure, pass(self) :: read_stats_3d
    procedure, pass(self) :: write_stats_3d_serial
    procedure, pass(self) :: read_stats_3d_serial
    procedure, pass(self) :: recyc_prepare
    procedure, pass(self) :: add_synthetic_perturbations
    procedure, pass(self) :: sliceprobe_prepare
    procedure, pass(self) :: correct_bc_order
!
!
!   
!
  endtype equation_singleideal_object
!
contains
!
  subroutine initialize(self, filename)
!   < Initialize the equation.
    class(equation_singleideal_object), intent(inout) :: self              !< The equation.
    character(*)                    , intent(in)      :: filename          !< Input file name.
    logical, dimension(3) :: periodic
    integer, dimension(3) :: mpi_splits
    integer :: mpi_split_x ,mpi_split_y ,mpi_split_z
    integer :: nxmax, nymax, nzmax, ng
    integer :: grid_type, metrics_order
    real(rkind) :: domain_size_x,  domain_size_y,  domain_size_z, rtemp
    real(rkind), dimension(:), allocatable :: grid_vars
    logical :: rebuild_ghost, ystaggering
    integer :: bctag_temp, order
    real(rkind) :: d1_temp(1:4), d2_temp(0:4)
    integer :: i
!
!   Get MPI basic quantities
    call get_mpi_basic_info(self%nprocs, self%myrank, self%masterproc, self%mpi_err)
!
!   Read input
    call self%read_input(filename)
!
!   Enable memory debugger
    if (self%cfg%has_key("output","debug_memory")) call self%cfg%get("output","debug_memory",self%debug_memory)
!
!   Enable async pattern
    if (self%cfg%has_key("numerics","mode_async")) call self%cfg%get("numerics","mode_async",self%mode_async)
!
!   Random numbers' configuration
    call self%cfg%get("numerics","rand_type",self%rand_type)
    if(self%rand_type < 0) then
      call init_crandom_f(self%myrank,reproducible=.false.)
      if (self%masterproc) write(*,*) 'Random numbers NOT reproducible'
    else
      call init_crandom_f(self%myrank,reproducible=.true.)
      if (self%masterproc) write(*,*) 'Random numbers reproducible'
    endif
!
!   Normalization
    self%l0    = 1._rkind
    self%t0    = 1._rkind
    self%rho0  = 1._rkind
    self%rgas0 = 1._rkind
    if (self%cfg%has_key("ref_prop","l0")) call self%cfg%get("ref_prop","l0",self%l0)
    if (self%cfg%has_key("ref_prop","t0")) call self%cfg%get("ref_prop","t0",self%t0)
    if (self%cfg%has_key("ref_prop","rho0")) call self%cfg%get("ref_prop","rho0",self%rho0)
    if (self%cfg%has_key("ref_prop","rgas0")) call self%cfg%get("ref_prop","rgas0",self%rgas0)
!
!   Flow init
    call self%cfg%get("flow","flow_init",self%flow_init)
    if (self%flow_init<-1 .or.self%flow_init>2) call fail_input_any("flow_init not implemented")
    self%channel_case = .false.
    if (self%flow_init == 0) self%channel_case = .true.
    self%bl_case = .false.
    if (self%flow_init == 1) then
      self%bl_case    = .true.
      self%bl_laminar = .false. ! Default is turbulent BL
    endif
!
    ystaggering = .false.
    if (self%channel_case) ystaggering = .true.
!
!   Boundary conditions
    call self%cfg%get("bc","xmin",bctag_temp) ; self%bctags(1) = bctag_temp
    if (self%bctags(1)==9) self%bl_laminar = .true.
    call self%cfg%get("bc","xmax",bctag_temp) ; self%bctags(2) = bctag_temp
    call self%cfg%get("bc","ymin",bctag_temp) ; self%bctags(3) = bctag_temp
    call self%cfg%get("bc","ymax",bctag_temp) ; self%bctags(4) = bctag_temp
    call self%cfg%get("bc","zmin",bctag_temp) ; self%bctags(5) = bctag_temp
    call self%cfg%get("bc","zmax",bctag_temp) ; self%bctags(6) = bctag_temp
    periodic(1:3) = .false.
    if(any(self%bctags(1:2) == 0)) then
      self%bctags(1:2) = 0
      periodic(1) = .true.
    endif
    if(any(self%bctags(3:4) == 0)) then
      self%bctags(3:4) = 0
      periodic(2) = .true.
    endif
    if(any(self%bctags(5:6) == 0)) then
      self%bctags(5:6) = 0
      periodic(3) = .true.
    endif
    call self%cfg%get("bc","xmin_nr",bctag_temp) ; self%bctags_nr(1) = bctag_temp
    call self%cfg%get("bc","xmax_nr",bctag_temp) ; self%bctags_nr(2) = bctag_temp
    call self%cfg%get("bc","ymin_nr",bctag_temp) ; self%bctags_nr(3) = bctag_temp
    call self%cfg%get("bc","ymax_nr",bctag_temp) ; self%bctags_nr(4) = bctag_temp
    call self%cfg%get("bc","zmin_nr",bctag_temp) ; self%bctags_nr(5) = bctag_temp
    call self%cfg%get("bc","zmax_nr",bctag_temp) ; self%bctags_nr(6) = bctag_temp
!   Restart
    call self%cfg%get("controls","restart_type",self%restart_type)
    call self%cfg%get("grid","grid_type",grid_type)
    rebuild_ghost = .false.
    if (grid_type == GRID_FROMFILE) rebuild_ghost = .true.
    if (self%restart_type > 0) grid_type = GRID_FROMFILE
!
!   Grid
    call self%cfg%get("grid","nxmax",nxmax)
    call self%cfg%get("grid","nymax",nymax)
    call self%cfg%get("grid","nzmax",nzmax)
    call self%cfg%get("grid","ng",ng)
    call self%cfg%get("grid","domain_size_x",domain_size_x)
    call self%cfg%get("grid","domain_size_y",domain_size_y)
    call self%cfg%get("grid","domain_size_z",domain_size_z)
    call self%cfg%get("grid","metrics_order",metrics_order)
!
    select case (grid_type)
    case(GRID_FROMFILE)
    case(GRID_UNIFORM)
    case(GRID_CHA)
      allocate(grid_vars(4))
      call self%cfg%get("grid","jbgrid",rtemp)     ; grid_vars(1) = rtemp
      call self%cfg%get("grid","dyptarget",rtemp)  ; grid_vars(2) = rtemp
      call self%cfg%get("grid","ysmoosteps",rtemp) ; grid_vars(3) = rtemp
      call self%cfg%get("flow","Reynolds",rtemp)   ; grid_vars(4) = rtemp
    case(GRID_BL)
      allocate(grid_vars(6))
      call self%cfg%get("grid","jbgrid",rtemp)     ; grid_vars(1) = rtemp
      call self%cfg%get("grid","dyptarget",rtemp)  ; grid_vars(2) = rtemp
      call self%cfg%get("grid","nywr",rtemp)       ; grid_vars(3) = rtemp
      call self%cfg%get("grid","lywr",rtemp)       ; grid_vars(4) = rtemp
      call self%cfg%get("flow","Reynolds",rtemp)   ; grid_vars(5) = rtemp
      call self%cfg%get("grid","ysmoosteps",rtemp) ; grid_vars(6) = rtemp
    end select
!
    call self%grid%initialize(periodic, nxmax, nymax, nzmax, ng, grid_type, &
      domain_size_x, domain_size_y, domain_size_z, self%l0, &
      grid_vars, metrics_order, rebuild_ghost, ystaggering)
!
!   Number of variables and auxiliary variables
    self%nv         = 5
    self%nv_aux     = 10
    self%nv_stat    = 70
    self%nv_stat_3d = 14
!
!   Set fluid properties
    call self%set_fluid_prop()
!   Set flow parameters
    call self%set_flow_params()
!
!   Numerics
    call self%cfg%get("numerics","ep_order",self%ep_order)
    if (self%cfg%has_key("numerics","nkeep")) then
      call self%cfg%get("numerics","nkeep",self%nkeep)
    else
      self%nkeep = 0
    endif
    if (self%calorically_perfect/=1) then
      if (self%nkeep>0) then
        if (self%masterproc) write(*,*) "nkeep forced equal to zero for non calorically perfect gas"
        self%nkeep = 0
      endif
    endif
    call self%cfg%get("numerics","weno_scheme",self%weno_scheme)
    call self%cfg%get("numerics","weno_version",self%weno_version)
    if (self%cfg%has_key("numerics","flux_splitting")) then
      call self%cfg%get("numerics","flux_splitting",self%flux_splitting)
    else
      self%flux_splitting = 0
    endif
    call self%cfg%get("numerics","visc_order",self%visc_order)
    call self%cfg%get("numerics","conservative_viscous",self%conservative_viscous)
    allocate(self%supported_orders(1:4))
    self%supported_orders = [2,4,6,8]
    do i=1,size(self%supported_orders)
      order = self%supported_orders(i)
      call self%grid%get_deriv_coeffs(d1_temp, d2_temp, order, order)
      self%coeff_deriv1(:,i) = d1_temp
      self%coeff_deriv2(:,i) = d2_temp
    enddo
    call self%cfg%get("numerics","sensor_threshold",self%sensor_threshold)
!
!   Stats 3D
    call self%cfg%get("output","enable_plot3d",self%enable_plot3d)
    if (self%cfg%has_key("output","enable_stat_3d")) then
      call self%cfg%get("output","enable_stat_3d",self%enable_stat_3d)
    else
      self%enable_stat_3d = 0
    endif
!
!   MPI split
    call self%cfg%get("mpi","x_split",mpi_split_x) ; mpi_splits(1) = mpi_split_x
    call self%cfg%get("mpi","y_split",mpi_split_y) ; mpi_splits(2) = mpi_split_y
    call self%cfg%get("mpi","z_split",mpi_split_z) ; mpi_splits(3) = mpi_split_z
!
!   Field initialize
    call self%field%initialize(self%grid, self%nv, mpi_splits)
!
!   Allocate variables
    call self%alloc()
!
!   Fix boundary conditions considering MPI neighbours
!   print*,'before rank= ',self%myrank, ' - bc: ',self%bctags
    call self%field%correct_bc(self%bctags)
    call self%field%correct_bc(self%bctags_nr)
!
!   Pre-process boundary conditions
    call self%bc_preproc()
!   print*,'rank= ',self%myrank, ' - bc: ',self%bctags,' - bc_nr: ',self%bctags_nr
!
!   Initialize temporal integrator
    call self%cfg%get("numerics","rk_type",self%rk_type)
    call self%runge_kutta_initialize()
    call self%cfg%get("controls","cfl",self%cfl)
    call self%cfg%get("controls","num_iter",self%num_iter)
!
!
    call self%initial_conditions()
    select case(self%restart_type)
    case(0)
      self%time0    = 0._rkind
      self%icyc0    = 0
      self%itav     = 0
      self%istore   = 1
      self%time_from_last_rst    = 0._rkind
      self%time_from_last_write  = 0._rkind
      self%time_from_last_stat   = 0._rkind
      self%time_from_last_slice  = 0._rkind
      self%w_stat = 0._rkind
      if (self%enable_stat_3d>0)  self%w_stat_3d = 0._rkind
    case(1)
      if (self%io_type==1) call self%field%read_field_serial()
      if (self%io_type==2) call self%field%read_field()
      call self%read_field_info()
      self%itav   = 0
      self%w_stat = 0._rkind
      if (self%enable_stat_3d>0)  self%w_stat_3d = 0._rkind
    case(2)
      if (self%io_type==1) then
        call self%field%read_field_serial()
        call self%read_stats_serial()
        if (self%enable_stat_3d>0) call self%read_stats_3d_serial()
      endif
      if (self%io_type==2) then
        call self%field%read_field()
        call self%read_stats()
        if (self%enable_stat_3d>0) call self%read_stats_3d()
      endif
      call self%read_field_info()
    endselect
!
    self%recyc = .false.
!   only ncoords(1) procs have bctags 10, otherwise MPI tag replaced it
    if (self%field%ncoords(1) == 0) then
      if(self%bctags(1) == 10) then
        self%recyc = .true.
      endif
    endif
    call mpi_bcast(self%recyc,1,mpi_logical,0,self%field%mp_cartx,self%mpi_err)
    if (self%recyc) call self%recyc_prepare()
!
    self%time     = self%time0
    self%icyc     = self%icyc0
!
!   Save field and restart intervals
    call self%cfg%get("output","dtsave",self%dtsave)
    call self%cfg%get("output","dtsave_restart",self%dtsave_restart)
    call self%cfg%get("output","dtstat",self%dtstat)
    call self%cfg%get("output","enable_plot3d",self%enable_plot3d)
    call self%cfg%get("output","enable_vtk",self%enable_vtk)
    call self%cfg%get("output","dtslice",self%dtslice)
    call self%cfg%get("output","igslice",self%igslice)
    call self%cfg%get("output","jgslice",self%jgslice)
    call self%cfg%get("output","kgslice",self%kgslice)
    self%io_type = 2 ! (0 => no IO, 1 => serial, 2 => parallel)
    if (self%cfg%has_key("output","io_type")) then
      call self%cfg%get("output","io_type",self%io_type)
    endif
    call self%sliceprobe_prepare()
!
    call self%cfg%get("output","print_control",self%print_control)
    call self%cfg%get("controls","iter_dt_recompute",self%iter_dt_recompute)
!
!   Correct order of accuracy at the borders
    self%correct_bound_ord = 0
    if (self%cfg%has_key("bc","correct_bound_ord")) then
      call self%cfg%get("bc","correct_bound_ord",self%correct_bound_ord)
    endif
    if (self%correct_bound_ord>0) call self%correct_bc_order()
!
!
!
!   oldelse
!
!   oldendif
!
!
!
!   put as type default self%i_freeze = 0
!
  endsubroutine initialize
!
  subroutine correct_bc_order(self)
    class(equation_singleideal_object), intent(inout) :: self
!   
    integer :: i,j,k
    integer :: stencil_size
!   
    associate(nx => self%field%nx, ny => self%field%ny,nz => self%field%nz, &
      ng => self%grid%ng, ep_order => self%ep_order, weno_scheme => self%weno_scheme, &
      ncoords => self%field%ncoords, bctags => self%bctags, nblocks => self%field%nblocks)
!
      stencil_size = max(ep_order/2, weno_scheme)
!     
!     IMIN
      if ((ncoords(1)==0).and.(any(bctags(1:2)/=0))) then
        self%ep_ord_change_x(:,0,:)  = -stencil_size+1
        do i=1,stencil_size-1
          self%ep_ord_change_x(:,i,:) = -stencil_size+i
        enddo
      endif
!     
!     IMAX
      if ((ncoords(1)==(nblocks(1)-1)).and.(any(bctags(1:2)/=0))) then
        self%ep_ord_change_x(:,nx,:) = -stencil_size+1
        do i=1,stencil_size-1
          self%ep_ord_change_x(:,nx-i,:) = -stencil_size+i
        enddo
      endif
!     
!     JMIN
      if ((ncoords(2)==0).and.(any(bctags(3:4)/=0))) then
        self%ep_ord_change(:,0,:,2)  = -stencil_size+1
        do j=1,stencil_size-1
          self%ep_ord_change(:,j,:,2) = -stencil_size+j
        enddo
      endif
!     
!     JMAX
      if ((ncoords(2)==(nblocks(2)-1)).and.(any(bctags(3:4)/=0))) then
        self%ep_ord_change(:,ny,:,2) = -stencil_size+1
        do j=1,stencil_size-1
          self%ep_ord_change(:,ny-j,:,2) = -stencil_size+j
        enddo
      endif
!     
!     KMIN
      if ((ncoords(3)==0).and.(any(bctags(5:6)/=0))) then
        self%ep_ord_change(:,:,0,3)  = -stencil_size+1
        do k=1,stencil_size-1
          self%ep_ord_change(:,:,k,3) = -stencil_size+k
        enddo
      endif
!     
!     KMAX
      if ((ncoords(3)==(nblocks(3)-1)).and.(any(bctags(5:6)/=0))) then
        self%ep_ord_change(:,:,nz,3) = -stencil_size+1
        do k=1,stencil_size-1
          self%ep_ord_change(:,:,nz-k,3) = -stencil_size+k
        enddo
      endif
!
    endassociate
  end subroutine correct_bc_order
!
!
!
!
!
!
!
  subroutine sliceprobe_prepare(self)
    class(equation_singleideal_object), intent(inout) :: self
!
    integer :: i, j, k, l, ip, jp, kp, n
    integer :: ii, jj, kk, ll
    integer :: icord, jcord, kcord
    integer :: inum, jnum, knum
    real(rkind) :: xp, yp, zp
    real(rkind) :: x0, y0, z0, dxloc, dyloc, dzloc
    real(rkind) :: xyz1, xyz2, xyz3
    logical :: probe_exists
    logical :: in_i, in_j, in_k
    logical, dimension(:), allocatable :: in_ijk
    real(rkind), dimension(8,8) :: amat3d
    real(rkind), dimension(1,8) :: xtrasp3d,alftrasp3d
!
    associate(ng => self%grid%ng, nx => self%field%nx, ny => self%field%ny, nz => self%field%nz, &
      x => self%field%x, y => self%field%y, z => self%field%z)
!
      inum = 0
      do i = 1,size(self%igslice)
        if (self%igslice(i)>0) then
          icord = (self%igslice(i)-1)/nx
          if (self%field%ncoords(1)==icord) inum = inum + 1
        endif
      enddo
      jnum = 0
      do j = 1,size(self%jgslice)
        if (self%jgslice(j)>0) then
          jcord = (self%jgslice(j)-1)/ny
          if (self%field%ncoords(2)==jcord) jnum = jnum + 1
        endif
      enddo
      knum = 0
      do k = 1,size(self%kgslice)
        if (self%kgslice(k)>0) then
          kcord = (self%kgslice(k)-1)/nz
          if (self%field%ncoords(3)==kcord) knum = knum + 1
        endif
      enddo
!     
      if (inum>0) then
        allocate(self%islice(inum))
        allocate(sliceyz_aux(inum,1-ng:ny+ng,1-ng:nz+ng,6))
      endif
      if (jnum>0) then
        allocate(self%jslice(jnum))
        allocate(slicexz_aux(1-ng:nx+ng,jnum,1-ng:nz+ng,6))
      endif
      if (knum>0) then
        allocate(self%kslice(knum))
        allocate(slicexy_aux(1-ng:nx+ng,1-ng:ny+ng,knum,6))
      endif
!     
      inum = 0
      do i = 1,size(self%igslice)
        if (self%igslice(i)>0) then
          icord = (self%igslice(i)-1)/nx
          if (self%field%ncoords(1)==icord) then
            inum = inum+1
            self%islice(inum) = self%igslice(i)-self%field%ncoords(1)*nx
          endif
        endif
      enddo
      jnum = 0
      do j = 1,size(self%jgslice)
        if (self%jgslice(j)>0) then
          jcord = (self%jgslice(j)-1)/ny
          if (self%field%ncoords(2)==jcord) then
            jnum = jnum + 1
            self%jslice(jnum) = self%jgslice(j)-self%field%ncoords(2)*ny
          endif
        endif
      enddo
      knum = 0
      do k = 1,size(self%kgslice)
        if (self%kgslice(k)>0) then
          kcord = (self%kgslice(k)-1)/nz
          if (self%field%ncoords(3)==kcord) then
            knum = knum + 1
            self%kslice(knum) = self%kgslice(k)-self%field%ncoords(3)*nz
          endif
        endif
      enddo
!
      inquire(file='probe_list.dat',exist=probe_exists)
      self%num_probe = 0
      if (probe_exists) then
        open(18,file='probe_list.dat')
        read(18,*) self%num_probe_tot
        allocate(in_ijk(self%num_probe_tot))
        in_ijk = .false.
        do l=1,self%num_probe_tot
          read(18,*) xp, yp, zp
          in_i   = .false.
          in_j   = .false.
          in_k   = .false.
          if (xp>=self%field%x(1).and.xp<self%field%x(nx+1)) in_i = .true.
          if (yp>=self%field%y(1).and.yp<self%field%y(ny+1)) in_j = .true.
          if (zp>=self%field%z(1).and.zp<self%field%z(nz+1)) in_k = .true.
          if (in_i.and.in_j.and.in_k) then
            in_ijk(l) = .true.
            self%num_probe = self%num_probe+1
          endif
        enddo
        close(18)
!
        if (self%num_probe>0) then
          allocate(self%probe_coord(3,self%num_probe))
          allocate(self%w_aux_probe(6,self%num_probe))
          allocate(self%ijk_probe(3,self%num_probe))
          allocate(self%probe_coeff(2,2,2,self%num_probe))
        endif
!
        open(18,file='probe_list.dat')
        read(18,*) self%num_probe_tot
        n = 0
        do l=1,self%num_probe_tot
          read(18,*) xp, yp, zp
          if (in_ijk(l)) then
            n = n+1
            self%probe_coord(1,n) = xp
            self%probe_coord(2,n) = yp
            self%probe_coord(3,n) = zp
            call locateval(x(1:nx),nx,xp,ip)
            call locateval(y(1:ny),ny,yp,jp)
            call locateval(z(1:nz),nz,zp,kp)
            self%ijk_probe(1,n) = ip
            self%ijk_probe(2,n) = jp
            self%ijk_probe(3,n) = kp
!           
!           Find coefficients for trilinear interpolation
!           
            x0 = x(ip)
            y0 = y(jp)
            z0 = z(kp)
            dxloc = x(ip+1)-x(ip)
            dyloc = y(jp+1)-y(jp)
            dzloc = z(kp+1)-z(kp)
            xp = (xp-x0)/dxloc
            yp = (yp-y0)/dyloc
            zp = (zp-z0)/dzloc
!           
            xtrasp3d(1,1) = xp*yp*zp
            xtrasp3d(1,2) = xp*yp
            xtrasp3d(1,3) = xp*zp
            xtrasp3d(1,4) = yp*zp
            xtrasp3d(1,5) = xp
            xtrasp3d(1,6) = yp
            xtrasp3d(1,7) = zp
            xtrasp3d(1,8) = 1._rkind
!           
!           Dirichlet
            ll = 0
            do kk=0,1
              do jj=0,1
                do ii=0,1
                  ll = ll+1
                  xyz1 = x(ip+ii)
                  xyz2 = y(jp+jj)
                  xyz3 = z(kp+kk)
                  xyz1 = (xyz1-x0)/dxloc
                  xyz2 = (xyz2-y0)/dyloc
                  xyz3 = (xyz3-z0)/dzloc
                  amat3d(ll,:) = [xyz1*xyz2*xyz3,xyz1*xyz2,xyz1*xyz3,xyz2*xyz3,xyz1,xyz2,xyz3,1._rkind]
                enddo
              enddo
            enddo
!           
            call invmat(amat3d,8)
            alftrasp3d = matmul(xtrasp3d,amat3d)
            self%probe_coeff(1,1,1,n) = alftrasp3d(1,1)
            self%probe_coeff(2,1,1,n) = alftrasp3d(1,2)
            self%probe_coeff(1,2,1,n) = alftrasp3d(1,3)
            self%probe_coeff(2,2,1,n) = alftrasp3d(1,4)
            self%probe_coeff(1,1,2,n) = alftrasp3d(1,5)
            self%probe_coeff(2,1,2,n) = alftrasp3d(1,6)
            self%probe_coeff(1,2,2,n) = alftrasp3d(1,7)
            self%probe_coeff(2,2,2,n) = alftrasp3d(1,8)
          endif
!
        enddo
        close(18)
      endif
!
!
    endassociate
!
  endsubroutine sliceprobe_prepare
!
  subroutine recyc_prepare(self)
    class(equation_singleideal_object), intent(inout) :: self
    integer :: ig_recyc, j
!
    call self%cfg%get("bc","x_recyc",self%x_recyc)
    self%x_recyc = self%x_recyc*self%l0
!
    associate(ng => self%grid%ng, ny => self%field%ny)
      allocate(self%yplus_inflow(1-ng:ny+ng))
      allocate(self%eta_inflow(1-ng:ny+ng))
      allocate(self%yplus_recyc(1-ng:ny+ng))
      allocate(self%eta_recyc(1-ng:ny+ng))
      allocate(self%map_j_inn(1:ny))
      allocate(self%map_j_out(1:ny))
      allocate(self%weta_inflow(1:ny))
    endassociate
!
    associate(xg => self%grid%xg, nxmax => self%grid%nxmax, nx => self%field%nx, ny => self%field%ny, &
      ng => self%grid%ng, xrecyc => self%x_recyc, i_recyc => self%i_recyc, ib_recyc => self%ib_recyc, &
      deltavec => self%deltavec, deltavvec => self%deltavvec, betarecyc => self%betarecyc, &
      y => self%field%y, yplus_inflow => self%yplus_inflow, eta_inflow => self%eta_inflow, &
      yplus_recyc  => self%yplus_recyc,  eta_recyc  => self%eta_recyc, l0 => self%l0, &
      map_j_inn => self%map_j_inn, map_j_out => self%map_j_out, weta_inflow => self%weta_inflow &
      )
      call locateval(xg(1:nxmax),nxmax,xrecyc,ig_recyc) ! xrecyc is between xg(ii) and xg(ii+1), ii is between 0 and nxmax
      ib_recyc = (ig_recyc-1)/nx
      i_recyc  = ig_recyc-nx*ib_recyc
      if (i_recyc<ng) then
        i_recyc  = ng
        ig_recyc = i_recyc+nx*ib_recyc
      endif
      if (self%masterproc) write(*,*) 'Recycling station exactly at = ', xg(ig_recyc)
!
      betarecyc = deltavvec(ig_recyc)/deltavvec(1)
      if (self%masterproc) write(*,*) 'Recycling beta factor = ', betarecyc
      if (self%masterproc) write(*,*) 'Pirozzoli-Ceci beta factor = ', (deltavec(ig_recyc)/deltavec(1))**0.13_rkind
!     check which Reynolds is the next in the formula
!     if (self%masterproc) write(*,*) 'xrecyc: ',xrecyc, self%Reynolds
      if (self%masterproc) then
        write(*,*) 'Urbin-Knight beta factor = ', &
          ((1._rkind+(xrecyc/l0)*0.27_rkind**1.2_rkind/self%Reynolds**0.2_rkind)**(5._rkind/6._rkind))**0.1_rkind
      endif
!
      do j=1-ng,ny+ng
        yplus_inflow(j) = y(j)/deltavvec(1)
        eta_inflow(j)   = y(j)/deltavec (1)
        yplus_recyc(j)  = y(j)/deltavvec(ig_recyc)
        eta_recyc(j)    = y(j)/deltavec (ig_recyc)
      enddo
!     
      do j=1,ny
        call locateval(yplus_recyc(1:ny),ny,yplus_inflow(j),map_j_inn(j))
        call locateval(eta_recyc(1:ny),ny,eta_inflow(j),map_j_out(j))
        weta_inflow(j) = 0.5_rkind*(1._rkind+tanh((4._rkind*(eta_inflow(j)-0.2_rkind))&
        /(0.2_rkind+eta_inflow(j)*0.6_rkind))/tanh(4._rkind))
      enddo
    endassociate
  endsubroutine recyc_prepare
!
  subroutine alloc(self)
    class(equation_singleideal_object), intent(inout) :: self
    associate(nx => self%field%nx, ny => self%field%ny, nz => self%field%nz,  &
      ng => self%grid%ng, nv => self%nv, nv_aux => self%nv_aux, nv_stat => self%nv_stat, &
      nv_stat_3d => self%nv_stat_3d, enable_stat_3d => self%enable_stat_3d)
      allocate(self%w_aux(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng, nv_aux))
      allocate(self%w_stat(nv_stat, 1:nx, 1:ny))
      allocate(self%fluid_mask(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng))
      allocate(self%ep_ord_change(0:nx, 0:ny, 0:nz, 2:3))
      allocate(self%ep_ord_change_x(0:ny, 0:nx, 0:nz))
      if (enable_stat_3d>0) allocate(self%w_stat_3d(nv_stat_3d, 1:nx, 1:ny, 1:nz))
      self%fluid_mask     = 0
      self%ep_ord_change   = 0
      self%ep_ord_change_x = 0
    endassociate
  endsubroutine alloc
!
  subroutine read_stats(self)
    class(equation_singleideal_object), intent(inout) :: self
    integer, dimension(3) :: sizes    ! Dimensions of the total grid
    integer, dimension(3) :: subsizes ! Dimensions of grid local to a procs
    integer, dimension(3) :: starts   ! Starting coordinates
    integer :: ntotxy
    integer, dimension(mpi_status_size) :: istatus
    integer :: mpi_io_file
    integer :: filetype, l, m
    integer :: size_real
    integer (kind=mpi_offset_kind) :: offset
    character(len=256) :: oldname, newname
    associate(nx => self%field%nx, ny => self%field%ny, nv_stat => self%nv_stat, &
      w_stat => self%w_stat, mp_cartx => self%field%mp_cartx, mp_cartz => self%field%mp_cartz, &
      ncoords => self%field%ncoords, nblocks => self%field%nblocks, iermpi => self%mpi_err)
      if (ncoords(3)==0) then
        sizes(1) = nblocks(1)*nx
        sizes(2) = nblocks(2)*ny
        sizes(3) = 1
        subsizes(1) = nx
        subsizes(2) = ny
        subsizes(3) = 1
        starts(1) = 0 + ncoords(1)*subsizes(1)
        starts(2) = 0 + ncoords(2)*subsizes(2)
        starts(3) = 0
        ntotxy = nx*ny
!
        call mpi_type_create_subarray(3,sizes,subsizes,starts,mpi_order_fortran,mpi_prec,filetype,iermpi)
        call mpi_type_commit(filetype,iermpi)
        call mpi_file_open(mp_cartx,'stat.bin',mpi_mode_rdonly,mpi_info_null,mpi_io_file,iermpi)
        offset = 0
        do l=1,nv_stat
          call mpi_file_set_view(mpi_io_file,offset,mpi_prec,filetype,"native",mpi_info_null,iermpi)
          call mpi_file_read_all(mpi_io_file,w_stat(l,1:nx,1:ny),ntotxy,mpi_prec,istatus,iermpi)
          call mpi_type_size(mpi_prec,size_real,iermpi)
          do m=1,nblocks(1)*nblocks(2)
            offset = offset+size_real*ntotxy
          enddo
        enddo
!
        call mpi_file_close(mpi_io_file,iermpi)
        call mpi_type_free(filetype,iermpi)
      endif
!     
      call mpi_bcast(w_stat,nv_stat*nx*ny,mpi_prec,0,mp_cartz,iermpi)
!     
      call mpi_barrier(mpi_comm_world, iermpi)
      if (self%masterproc) then
        oldname = c_char_"stat.bin"//c_null_char
        newname = c_char_"stat.bak"//c_null_char
        iermpi = rename_wrapper(oldname, newname)
        if (iermpi /= 0) write(error_unit,*) "Warning! Cannot rename file stat.bin to stat.bak"
      endif
      call mpi_barrier(mpi_comm_world, iermpi)
    endassociate
!
  endsubroutine read_stats
!
  subroutine read_stats_serial(self)
    class(equation_singleideal_object), intent(inout) :: self
    character(4) :: chx,chy
    associate(nx => self%field%nx, ny => self%field%ny, nz => self%field%nz, nv_stat => self%nv_stat, &
      w_stat => self%w_stat, ncoords => self%field%ncoords, mp_cartz => self%field%mp_cartz, &
      iermpi => self%mpi_err)
!
      if (ncoords(3)==0) then
!
        if (self%masterproc) write(*,*) 'Reading stat0_XXX_XXX.bin'
        1004 format(I4.4)
        write(chx,1004) ncoords(1)
        write(chy,1004) ncoords(2)
!
        open (11,file='stat0_'//chx//'_'//chy//'.bin',form='unformatted')
        read(11) w_stat(1:nv_stat,1:nx,1:ny)
        close(11)
!
      endif
!     
      call mpi_bcast(w_stat,nv_stat*nx*ny,mpi_prec,0,mp_cartz,iermpi)
!     
    endassociate
!
  endsubroutine read_stats_serial
!
  subroutine read_stats_3d(self)
    class(equation_singleideal_object), intent(inout) :: self
    integer, dimension(3) :: sizes    ! Dimensions of the total grid
    integer, dimension(3) :: subsizes ! Dimensions of grid local to a procs
    integer, dimension(3) :: starts   ! Starting coordinates
    integer :: ntot3d
    integer, dimension(mpi_status_size) :: istatus
    integer :: mpi_io_file
    integer :: filetype, l, m
    integer :: size_real
    integer (kind=mpi_offset_kind) :: offset
    character(len=256) :: oldname, newname
!
    associate(nx => self%field%nx, ny => self%field%ny, nz => self%field%nz, nv_stat_3d => self%nv_stat_3d, &
      w_stat_3d => self%w_stat_3d, mp_cart => self%field%mp_cart, &
      ncoords => self%field%ncoords, nblocks => self%field%nblocks, iermpi => self%mpi_err)
!
      sizes(1) = nblocks(1)*nx
      sizes(2) = nblocks(2)*ny
      sizes(3) = nblocks(3)*nz
      subsizes(1) = nx
      subsizes(2) = ny
      subsizes(3) = nz
      starts(1) = 0 + ncoords(1)*subsizes(1)
      starts(2) = 0 + ncoords(2)*subsizes(2)
      starts(3) = 0 + ncoords(3)*subsizes(3)
      ntot3d = nx*ny*nz
!
      call mpi_type_create_subarray(3,sizes,subsizes,starts,mpi_order_fortran,mpi_prec,filetype,iermpi)
      call mpi_type_commit(filetype,iermpi)
      call mpi_file_open(mp_cart,'stat3d.bin',mpi_mode_rdonly,mpi_info_null,mpi_io_file,iermpi)
      offset = 0
      do l=1,nv_stat_3d
        call mpi_file_set_view(mpi_io_file,offset,mpi_prec,filetype,"native",mpi_info_null,iermpi)
        call mpi_file_read_all(mpi_io_file,w_stat_3d(l,1:nx,1:ny,1:nz),ntot3d,mpi_prec,istatus,iermpi)
        call mpi_type_size(mpi_prec,size_real,iermpi)
        do m=1,nblocks(1)*nblocks(2)*nblocks(3)
          offset = offset+size_real*ntot3d
        enddo
      enddo
!
      call mpi_file_close(mpi_io_file,iermpi)
      call mpi_type_free(filetype,iermpi)
!     
      if (self%masterproc) then
        oldname = c_char_"stat3d.bin"//c_null_char
        newname = c_char_"stat3d.bak"//c_null_char
        iermpi = rename_wrapper(oldname, newname)
        if (iermpi /= 0) write(error_unit,*) "Warning! Cannot rename file stat3d.bin to stat3d.bak"
      endif
      call mpi_barrier(mpi_comm_world, iermpi)
!     
    endassociate
!
  endsubroutine read_stats_3d
!
  subroutine write_stats(self)
    class(equation_singleideal_object), intent(inout) :: self
    integer, dimension(3) :: sizes    ! Dimensions of the total grid
    integer, dimension(3) :: subsizes ! Dimensions of grid local to a procs
    integer, dimension(3) :: starts   ! Starting coordinates
    integer :: ntotxy
    integer, dimension(mpi_status_size) :: istatus
    integer :: mpi_io_file
    integer :: filetype, l, m
    integer :: size_real
    integer (kind=mpi_offset_kind) :: offset
    associate(nx => self%field%nx, ny => self%field%ny, nv_stat => self%nv_stat, &
      w_stat => self%w_stat, mp_cartx => self%field%mp_cartx, &
      ncoords => self%field%ncoords, nblocks => self%field%nblocks, iermpi => self%mpi_err)
      if (ncoords(3)==0) then
        sizes(1) = nblocks(1)*nx
        sizes(2) = nblocks(2)*ny
        sizes(3) = 1
        subsizes(1) = nx
        subsizes(2) = ny
        subsizes(3) = 1
        starts(1) = 0 + ncoords(1)*subsizes(1)
        starts(2) = 0 + ncoords(2)*subsizes(2)
        starts(3) = 0
        ntotxy = nx*ny
!
        call mpi_type_create_subarray(3,sizes,subsizes,starts,mpi_order_fortran,mpi_prec,filetype,iermpi)
        call mpi_type_commit(filetype,iermpi)
        call mpi_file_open(mp_cartx,'stat.bin',mpi_mode_create+mpi_mode_wronly,mpi_info_null,mpi_io_file,iermpi)
        offset = 0
        do l=1,nv_stat
          call mpi_file_set_view(mpi_io_file,offset,mpi_prec,filetype,"native",mpi_info_null,iermpi)
          call mpi_file_write_all(mpi_io_file,w_stat(l,1:nx,1:ny),ntotxy,mpi_prec,istatus,iermpi)
          call mpi_type_size(mpi_prec,size_real,iermpi)
          do m=1,nblocks(1)*nblocks(2)
            offset = offset+size_real*ntotxy
          enddo
        enddo
!
        call mpi_file_close(mpi_io_file,iermpi)
        call mpi_type_free(filetype,iermpi)
      endif
    endassociate
  endsubroutine write_stats
!
  subroutine write_stats_serial(self)
    class(equation_singleideal_object), intent(inout) :: self
    character(4) :: chx,chy
    associate(nx => self%field%nx, ny => self%field%ny, nz => self%field%nz, nv_stat => self%nv_stat, &
      w_stat => self%w_stat, ncoords => self%field%ncoords, mp_cartz => self%field%mp_cartz)
!
      if (ncoords(3)==0) then
!
        if (self%masterproc) write(*,*) 'Writing stat1_XXX_XXX.bin'
        1004 format(I4.4)
        write(chx,1004) ncoords(1)
        write(chy,1004) ncoords(2)
!
        open (11,file='stat1_'//chx//'_'//chy//'.bin',form='unformatted')
        write(11) w_stat(1:nv_stat,1:nx,1:ny)
        close(11)
!
      endif
!     
    endassociate
!
  endsubroutine write_stats_serial
!
  subroutine write_stats_3d(self)
    class(equation_singleideal_object), intent(inout) :: self
    integer, dimension(3) :: sizes    ! Dimensions of the total grid
    integer, dimension(3) :: subsizes ! Dimensions of grid local to a procs
    integer, dimension(3) :: starts   ! Starting coordinates
    integer :: ntot3d
    integer, dimension(mpi_status_size) :: istatus
    integer :: mpi_io_file
    integer :: filetype, l, m
    integer :: size_real
    integer (kind=mpi_offset_kind) :: offset
    associate(nx => self%field%nx, ny => self%field%ny, nz => self%field%nz, nv_stat_3d => self%nv_stat_3d, &
      w_stat_3d => self%w_stat_3d, mp_cart => self%field%mp_cart, &
      ncoords => self%field%ncoords, nblocks => self%field%nblocks, iermpi => self%mpi_err)
!
      sizes(1) = nblocks(1)*nx
      sizes(2) = nblocks(2)*ny
      sizes(3) = nblocks(3)*nz
      subsizes(1) = nx
      subsizes(2) = ny
      subsizes(3) = nz
      starts(1) = 0 + ncoords(1)*subsizes(1)
      starts(2) = 0 + ncoords(2)*subsizes(2)
      starts(3) = 0 + ncoords(3)*subsizes(3)
      ntot3d = nx*ny*nz
!
      call mpi_type_create_subarray(3,sizes,subsizes,starts,mpi_order_fortran,mpi_prec,filetype,iermpi)
      call mpi_type_commit(filetype,iermpi)
      call mpi_file_open(mp_cart,'stat3d.bin',mpi_mode_create+mpi_mode_wronly,mpi_info_null,mpi_io_file,iermpi)
      offset = 0
      do l=1,nv_stat_3d
        call mpi_file_set_view(mpi_io_file,offset,mpi_prec,filetype,"native",mpi_info_null,iermpi)
        call mpi_file_write_all(mpi_io_file,w_stat_3d(l,1:nx,1:ny,1:nz),ntot3d,mpi_prec,istatus,iermpi)
        call mpi_type_size(mpi_prec,size_real,iermpi)
        do m=1,nblocks(1)*nblocks(2)*nblocks(3)
          offset = offset+size_real*ntot3d
        enddo
      enddo
!
      call mpi_file_close(mpi_io_file,iermpi)
      call mpi_type_free(filetype,iermpi)
!
    endassociate
  endsubroutine write_stats_3d
! 
  subroutine write_stats_3d_serial(self)
    class(equation_singleideal_object), intent(inout) :: self
    character(4) :: chx,chy,chz
    associate(nx => self%field%nx, ny => self%field%ny, nz => self%field%nz, nv_stat_3d => self%nv_stat_3d, &
      w_stat_3d => self%w_stat_3d, ncoords => self%field%ncoords)
!
      if (self%masterproc) write(*,*) 'Writing stat3d1_XXX_XXX_XXX.bin'
      1004 format(I4.4)
      write(chx,1004) ncoords(1)
      write(chy,1004) ncoords(2)
      write(chz,1004) ncoords(3)
!
      open (11,file='stat3d1_'//chx//'_'//chy//'_'//chz//'.bin',form='unformatted')
      write(11) w_stat_3d(1:nv_stat_3d,1:nx,1:ny,1:nz)
      close(11)
!
    endassociate
  endsubroutine write_stats_3d_serial
! 
  subroutine read_stats_3d_serial(self)
    class(equation_singleideal_object), intent(inout) :: self
    character(4) :: chx,chy,chz
    associate(nx => self%field%nx, ny => self%field%ny, nz => self%field%nz, nv_stat_3d => self%nv_stat_3d, &
      w_stat_3d => self%w_stat_3d, ncoords => self%field%ncoords)
!
      if (self%masterproc) write(*,*) 'Reading stat3d0_XXX_XXX_XXX.bin'
      1004 format(I4.4)
      write(chx,1004) ncoords(1)
      write(chy,1004) ncoords(2)
      write(chz,1004) ncoords(3)
!
      open (11,file='stat3d0_'//chx//'_'//chy//'_'//chz//'.bin',form='unformatted')
      read(11) w_stat_3d(1:nv_stat_3d,1:nx,1:ny,1:nz)
      close(11)
!
    endassociate
  endsubroutine read_stats_3d_serial
! 
  subroutine compute_stats(self)
    class(equation_singleideal_object), intent(inout) :: self
    real(rkind), dimension(self%nv_stat, self%field%nx, self%field%ny) :: w_stat_z
    real(rkind) :: rho, rhou, rhov, rhow, rhoe, ri, uu, vv, ww, qq, pp, tt, rho2, uu2, vv2, ww2, pp2, tt2, mu, nu, uv
    real(rkind) :: ux,uy,uz,vx,vy,vz,wx,wy,wz,omx,omy,omz,divl,div3l,machlocal,gamloc,c,t_tot,ccl,cploc
    real(rkind), dimension(3,3) :: sig
    integer :: i,j,l,k,npt,lmax
    associate(nx => self%field%nx, ny => self%field%ny, nz => self%field%nz, nv_stat => self%nv_stat, &
      ng => self%grid%ng, nv_aux => self%nv_aux, w_stat => self%w_stat, nzmax => self%grid%nzmax, &
      w => self%field%w, visc_model => self%visc_model, &
      powerlaw_vtexp => self%powerlaw_vtexp, mu0 => self%mu0, &
      T_ref_dim => self%T_ref_dim, itav => self%itav, mp_cartz => self%field%mp_cartz, &
      sutherland_S => self%sutherland_S, &
      mpi_err => self%mpi_err, cv_coeff => self%cv_coeff, w_aux => self%w_aux, &
      dcsidx => self%field%dcsidx, detady => self%field%detady, dzitdz => self%field%dzitdz, &
      indx_cp_l => self%indx_cp_l, indx_cp_r => self%indx_cp_r, cp_coeff => self%cp_coeff, &
      calorically_perfect => self%calorically_perfect, rgas0 => self%rgas0, t0 => self%t0)
!
      lmax = self%visc_order/2
!
      w_stat_z = 0._rkind
      do k=1,nz
        do j=1,ny
          do i=1,nx
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
              ccl = self%coeff_deriv1(l,lmax)
              ux = ux+ccl*(w(2,i+l,j,k)/w(1,i+l,j,k)-w(2,i-l,j,k)/w(1,i-l,j,k))
              vx = vx+ccl*(w(3,i+l,j,k)/w(1,i+l,j,k)-w(3,i-l,j,k)/w(1,i-l,j,k))
              wx = wx+ccl*(w(4,i+l,j,k)/w(1,i+l,j,k)-w(4,i-l,j,k)/w(1,i-l,j,k))
              uy = uy+ccl*(w(2,i,j+l,k)/w(1,i,j+l,k)-w(2,i,j-l,k)/w(1,i,j-l,k))
              vy = vy+ccl*(w(3,i,j+l,k)/w(1,i,j+l,k)-w(3,i,j-l,k)/w(1,i,j-l,k))
              wy = wy+ccl*(w(4,i,j+l,k)/w(1,i,j+l,k)-w(4,i,j-l,k)/w(1,i,j-l,k))
              uz = uz+ccl*(w(2,i,j,k+l)/w(1,i,j,k+l)-w(2,i,j,k-l)/w(1,i,j,k-l))
              vz = vz+ccl*(w(3,i,j,k+l)/w(1,i,j,k+l)-w(3,i,j,k-l)/w(1,i,j,k-l))
              wz = wz+ccl*(w(4,i,j,k+l)/w(1,i,j,k+l)-w(4,i,j,k-l)/w(1,i,j,k-l))
            enddo
            ux = ux*dcsidx(i)
            vx = vx*dcsidx(i)
            wx = wx*dcsidx(i)
            uy = uy*detady(j)
            vy = vy*detady(j)
            wy = wy*detady(j)
            uz = uz*dzitdz(k)
            vz = vz*dzitdz(k)
            wz = wz*dzitdz(k)
            omx = wy-vz
            omy = uz-wx
            omz = vx-uy
            rho  = w(1,i,j,k)
            rhou = w(2,i,j,k)
            rhov = w(3,i,j,k)
            rhow = w(4,i,j,k)
            rhoe = w(5,i,j,k)
            ri   = 1._rkind/rho
            uu   = rhou*ri
            vv   = rhov*ri
            ww   = rhow*ri
            qq   = 0.5_rkind*(uu*uu+vv*vv+ww*ww)
!
!           pp   = gm1*(rhoe-rho*qq)
!           tt   = pp*ri
!
            tt = w_aux(i,j,k,6)
            pp = rho*rgas0*tt
            mu = w_aux(i,j,k,7)
            nu = mu/rho
!           
            gamloc    = get_gamloc(tt,t0,indx_cp_l,indx_cp_r,cp_coeff,calorically_perfect,rgas0)
            cploc     = gamloc/(gamloc-1._rkind)*rgas0
            c         = sqrt(gamloc*rgas0*tt)
            machlocal = sqrt(2._rkind*qq)/c
            t_tot     = tt+qq*(gamloc-1._rkind)/gamloc
!           
            rho2 = rho*rho
            uu2  = uu*uu
            vv2  = vv*vv
            ww2  = ww*ww
            uv   = uu*vv
            pp2  = pp*pp
            tt2  = tt*tt
!           
            divl  = (ux+vy+wz)
            div3l = divl/3._rkind
            sig(1,1) = 2._rkind*(ux-div3l)
            sig(1,2) = uy+vx
            sig(1,3) = uz+wx
            sig(2,1) = sig(1,2)
            sig(2,2) = 2._rkind*(vy-div3l)
            sig(2,3) = vz+wy
            sig(3,1) = sig(1,3)
            sig(3,2) = sig(2,3)
            sig(3,3) = 2._rkind*(wz-div3l)
            sig      = sig*mu
!           
            w_stat_z(1,i,j)  = w_stat_z(1,i,j)+rho
            w_stat_z(2,i,j)  = w_stat_z(2,i,j)+uu
            w_stat_z(3,i,j)  = w_stat_z(3,i,j)+vv
            w_stat_z(4,i,j)  = w_stat_z(4,i,j)+ww
            w_stat_z(5,i,j)  = w_stat_z(5,i,j)+pp
            w_stat_z(6,i,j)  = w_stat_z(6,i,j)+tt
            w_stat_z(7,i,j)  = w_stat_z(7,i,j)+rho2
            w_stat_z(8,i,j)  = w_stat_z(8,i,j)+uu2
            w_stat_z(9,i,j)  = w_stat_z(9,i,j)+vv2
            w_stat_z(10,i,j) = w_stat_z(10,i,j)+ww2
            w_stat_z(11,i,j) = w_stat_z(11,i,j)+pp2
            w_stat_z(12,i,j) = w_stat_z(12,i,j)+tt2
            w_stat_z(13,i,j) = w_stat_z(13,i,j)+rhou
            w_stat_z(14,i,j) = w_stat_z(14,i,j)+rhov
            w_stat_z(15,i,j) = w_stat_z(15,i,j)+rhow
            w_stat_z(16,i,j) = w_stat_z(16,i,j)+rhou*uu
            w_stat_z(17,i,j) = w_stat_z(17,i,j)+rhov*vv
            w_stat_z(18,i,j) = w_stat_z(18,i,j)+rhow*ww
            w_stat_z(19,i,j) = w_stat_z(19,i,j)+rhou*vv
            w_stat_z(20,i,j) = w_stat_z(20,i,j)+mu
            w_stat_z(21,i,j) = w_stat_z(21,i,j)+nu
!           
!           ! Vorticity fluctuations terms
            w_stat_z(22,i,j) = w_stat_z(22,i,j)+omx**2
            w_stat_z(23,i,j) = w_stat_z(23,i,j)+omy**2
            w_stat_z(24,i,j) = w_stat_z(24,i,j)+omz**2
!           ! Temperature fluctuations terms
            w_stat_z(25,i,j) = w_stat_z(25,i,j)+rho*tt
            w_stat_z(26,i,j) = w_stat_z(26,i,j)+rho*tt2
            w_stat_z(27,i,j) = w_stat_z(27,i,j)+t_tot
            w_stat_z(28,i,j) = w_stat_z(28,i,j)+rho*t_tot
            w_stat_z(29,i,j) = w_stat_z(29,i,j)+t_tot**2
            w_stat_z(30,i,j) = w_stat_z(30,i,j)+rhou*tt
            w_stat_z(31,i,j) = w_stat_z(31,i,j)+rhov*tt
            w_stat_z(32,i,j) = w_stat_z(32,i,j)+rhow*tt
!           ! Fluctuating Mach number
            w_stat_z(33,i,j) = w_stat_z(33,i,j)+machlocal
            w_stat_z(34,i,j) = w_stat_z(34,i,j)+machlocal**2
!           
!           ! Turbulent transport
!           w_stat_z(35,i,j) = w_stat_z(35,i,j)+rhou*(uu2+vv2+ww2)
!           w_stat_z(36,i,j) = w_stat_z(36,i,j)+rhov*(uu2+vv2+ww2)
!           w_stat_z(37,i,j) = w_stat_z(37,i,j)+pp*uu
!           w_stat_z(38,i,j) = w_stat_z(38,i,j)+pp*vv
!           ! Viscous diffusion
!           w_stat_z(39,i,j) = w_stat_z(39,i,j)+sig(1,1)*uu+sig(2,1)*vv+sig(3,1)*ww
!           w_stat_z(40,i,j) = w_stat_z(40,i,j)+sig(1,2)*uu+sig(2,2)*vv+sig(3,2)*ww
!           w_stat_z(41,i,j) = w_stat_z(41,i,j)+sig(1,1)
!           w_stat_z(42,i,j) = w_stat_z(42,i,j)+sig(1,2)
!           w_stat_z(43,i,j) = w_stat_z(43,i,j)+sig(1,3)
!           w_stat_z(44,i,j) = w_stat_z(44,i,j)+sig(2,2)
!           w_stat_z(45,i,j) = w_stat_z(45,i,j)+sig(2,3)
!           w_stat_z(46,i,j) = w_stat_z(46,i,j)+sig(3,3)
!           ! Dissipation
!           w_stat_z(47,i,j) = w_stat_z(47,i,j)+sig(1,1)*ux+sig(1,2)*uy+sig(1,3)*uz  &
!           +sig(2,1)*vx+sig(2,2)*vy+sig(2,3)*vz  &
!           +sig(3,1)*wx+sig(3,2)*wy+sig(3,3)*wz
!           !! K (Other terms due to compressibility)
!           w_stat_z(48,i,j) = w_stat_z(48,i,j)+pp*(ux+vy+wz)
!           
!           Turbulent transport
            w_stat_z(35,i,j) = w_stat_z(35,i,j)+rhou*uu2
            w_stat_z(36,i,j) = w_stat_z(36,i,j)+rhov*uu2
            w_stat_z(37,i,j) = w_stat_z(37,i,j)+rhou*vv2
            w_stat_z(38,i,j) = w_stat_z(38,i,j)+rhov*vv2
            w_stat_z(39,i,j) = w_stat_z(39,i,j)+rhou*ww2
            w_stat_z(40,i,j) = w_stat_z(40,i,j)+rhov*ww2
!           Pressure transport
            w_stat_z(41,i,j) = w_stat_z(41,i,j)+pp*uu
            w_stat_z(42,i,j) = w_stat_z(42,i,j)+pp*vv
!           Useful for compressibility terms
            w_stat_z(43,i,j) = w_stat_z(43,i,j)+sig(1,1)
            w_stat_z(44,i,j) = w_stat_z(44,i,j)+sig(1,2)
            w_stat_z(45,i,j) = w_stat_z(45,i,j)+sig(1,3)
            w_stat_z(46,i,j) = w_stat_z(46,i,j)+sig(2,2)
            w_stat_z(47,i,j) = w_stat_z(47,i,j)+sig(2,3)
            w_stat_z(48,i,j) = w_stat_z(48,i,j)+sig(3,3)
!           ! Viscous transport (diffusion)
!           11
            w_stat_z(49,i,j) = w_stat_z(49,i,j)+sig(1,1)*uu
            w_stat_z(50,i,j) = w_stat_z(50,i,j)+sig(1,2)*uu
!           22
            w_stat_z(51,i,j) = w_stat_z(51,i,j)+sig(2,1)*vv
            w_stat_z(52,i,j) = w_stat_z(52,i,j)+sig(2,2)*vv
!           33
            w_stat_z(53,i,j) = w_stat_z(53,i,j)+sig(3,1)*ww
            w_stat_z(54,i,j) = w_stat_z(54,i,j)+sig(3,2)*ww
!           12
            w_stat_z(55,i,j) = w_stat_z(55,i,j)+sig(1,1)*vv+sig(2,1)*uu
            w_stat_z(56,i,j) = w_stat_z(56,i,j)+sig(1,2)*vv+sig(2,2)*uu
!           Dissipation
!           ! 11
            w_stat_z(57,i,j) = w_stat_z(57,i,j)+sig(1,1)*ux+sig(1,2)*uy+sig(1,3)*uz
!           ! 22
            w_stat_z(58,i,j) = w_stat_z(58,i,j)+sig(2,1)*vx+sig(2,2)*vy+sig(2,3)*vz
!           ! 33
            w_stat_z(59,i,j) = w_stat_z(59,i,j)+sig(3,1)*wx+sig(3,2)*wy+sig(3,3)*wz
!           ! 12
            w_stat_z(60,i,j) = w_stat_z(60,i,j)+sig(1,1)*vx+sig(1,2)*(ux+vy)+sig(2,2)*uy+sig(1,3)*vz+sig(2,3)*uz
!           Pressure-strain redistribution
            w_stat_z(61,i,j) = w_stat_z(61,i,j)+pp*ux
            w_stat_z(62,i,j) = w_stat_z(62,i,j)+pp*vy
            w_stat_z(63,i,j) = w_stat_z(63,i,j)+pp*wz
            w_stat_z(64,i,j) = w_stat_z(64,i,j)+pp*(uy+vx)
!           Compressible dissipation
            w_stat_z(65,i,j)  = w_stat_z(65,i,j)+divl*divl
!           
            w_stat_z(66,i,j) = w_stat_z(66,i,j)+rho*tt2*tt
            w_stat_z(67,i,j) = w_stat_z(67,i,j)+rho*tt2*tt2
            w_stat_z(68,i,j) = w_stat_z(68,i,j)+rhou*uu2*uu
            w_stat_z(69,i,j) = w_stat_z(69,i,j)+cploc
            w_stat_z(70,i,j) = w_stat_z(70,i,j)+gamloc
!           
          enddo
        enddo
      enddo
!
      npt      = nv_stat*nx*ny
      call mpi_allreduce(MPI_IN_PLACE,w_stat_z,npt,mpi_prec,mpi_sum,mp_cartz,mpi_err)
      w_stat_z = w_stat_z/nzmax
!
      w_stat = w_stat*itav + w_stat_z
      w_stat = w_stat/(itav+1)
    endassociate
  endsubroutine compute_stats
!
  subroutine compute_stats_3d(self)
    class(equation_singleideal_object), intent(inout) :: self
    integer :: i,j,k
    associate(nx => self%field%nx, ny => self%field%ny, nz => self%field%nz, itav => self%itav, &
      w_stat_3d => self%w_stat_3d, w => self%field%w, w_aux => self%w_aux)
!
      do k=1,nz
        do j=1,ny
          do i=1,nx
            w_stat_3d(1,i,j,k)  = w_stat_3d(1,i,j,k) *itav + w(1,i,j,k)
            w_stat_3d(2,i,j,k)  = w_stat_3d(2,i,j,k) *itav + w(2,i,j,k)
            w_stat_3d(3,i,j,k)  = w_stat_3d(3,i,j,k) *itav + w(3,i,j,k)
            w_stat_3d(4,i,j,k)  = w_stat_3d(4,i,j,k) *itav + w(4,i,j,k)
            w_stat_3d(5,i,j,k)  = w_stat_3d(5,i,j,k) *itav + w(1,i,j,k)*w_aux(i,j,k,6)
            w_stat_3d(6,i,j,k)  = w_stat_3d(6,i,j,k) *itav + w(1,i,j,k)**2
            w_stat_3d(7,i,j,k)  = w_stat_3d(7,i,j,k) *itav + w(2,i,j,k)**2/w(1,i,j,k)
            w_stat_3d(8,i,j,k)  = w_stat_3d(8,i,j,k) *itav + w(3,i,j,k)**2/w(1,i,j,k)
            w_stat_3d(9,i,j,k)  = w_stat_3d(9,i,j,k) *itav + w(4,i,j,k)**2/w(1,i,j,k)
            w_stat_3d(10,i,j,k) = w_stat_3d(10,i,j,k)*itav + w(2,i,j,k)*w(3,i,j,k)/w(1,i,j,k)
            w_stat_3d(11,i,j,k) = w_stat_3d(11,i,j,k)*itav + w(2,i,j,k)*w(4,i,j,k)/w(1,i,j,k)
            w_stat_3d(12,i,j,k) = w_stat_3d(12,i,j,k)*itav + w(3,i,j,k)*w(4,i,j,k)/w(1,i,j,k)
            w_stat_3d(13,i,j,k) = w_stat_3d(13,i,j,k)*itav + w(1,i,j,k)*w_aux(i,j,k,6)**2
            w_stat_3d(14,i,j,k) = w_stat_3d(14,i,j,k)*itav + (w(1,i,j,k)*w_aux(i,j,k,6))**2
          enddo
        enddo
      enddo
      w_stat_3d = w_stat_3d/(itav+1)
!
    end associate
  end subroutine compute_stats_3d
!
  subroutine read_field_info(self)
    class(equation_singleideal_object), intent(inout) :: self
    character(len=256) :: oldname, newname
!
    self%field_info_cfg = parse_cfg("field_info.dat")
    call self%field_info_cfg%get("field_info","icyc0",                self%icyc0)
    call self%field_info_cfg%get("field_info","time0",                self%time0)
    call self%field_info_cfg%get("field_info","itav",                 self%itav)
    call self%field_info_cfg%get("field_info","time_from_last_rst",   self%time_from_last_rst)
    call self%field_info_cfg%get("field_info","time_from_last_write", self%time_from_last_write)
    call self%field_info_cfg%get("field_info","time_from_last_stat",  self%time_from_last_stat)
    call self%field_info_cfg%get("field_info","time_from_last_slice", self%time_from_last_slice)
    call self%field_info_cfg%get("field_info","istore",               self%istore)
!
!   open(unit=15, file="finaltime.dat")
!   read(15,*) self%icyc0
!   read(15,*) self%time0
!   read(15,*) self%itav
!   read(15,*) self%time_from_last_rst
!   read(15,*) self%time_from_last_write
!   read(15,*) self%time_from_last_stat
!   read(15,*) self%time_from_last_slice
!   read(15,*) self%istore
!   close(15)
!
    call mpi_barrier(mpi_comm_world, self%mpi_err)
    if (self%masterproc) then
      oldname = c_char_"field_info.dat"//c_null_char
      newname = c_char_"field_info.bak"//c_null_char
      self%mpi_err = rename_wrapper(oldname, newname)
      if (self%mpi_err /= 0) write(error_unit,*) "Warning! Cannot rename file field_info.dat to field_info.bak"
    endif
    call mpi_barrier(mpi_comm_world, self%mpi_err)
  endsubroutine read_field_info
!
  subroutine write_field_info(self)
    class(equation_singleideal_object), intent(inout) :: self
    if(self%masterproc) then
      call self%field_info_cfg%set("field_info","icyc0",                self%icyc)
      call self%field_info_cfg%set("field_info","time0",                self%time)
      call self%field_info_cfg%set("field_info","itav",                 self%itav)
      call self%field_info_cfg%set("field_info","time_from_last_rst",   self%time_from_last_rst)
      call self%field_info_cfg%set("field_info","time_from_last_write", self%time_from_last_write)
      call self%field_info_cfg%set("field_info","time_from_last_stat",  self%time_from_last_stat)
      call self%field_info_cfg%set("field_info","time_from_last_slice", self%time_from_last_slice)
      call self%field_info_cfg%set("field_info","istore",               self%istore)
      call self%field_info_cfg%write("field_info.dat")
!     open(unit=15, file="finaltime.dat")
!     write(15,*) self%icyc
!     write(15,*) self%time
!     write(15,*) self%itav
!     write(15,*) self%time_from_last_rst
!     write(15,*) self%time_from_last_write
!     write(15,*) self%time_from_last_stat
!     write(15,*) self%time_from_last_slice
!     write(15,*) self%istore
!     close(15)
    endif
  endsubroutine write_field_info
!
  function get_e_from_temperature(tt, t0, indx_cp_l, indx_cp_r, cv_coeff, calorically_perfect)
    real(rkind), intent(in) :: tt, t0
    integer, intent(in) :: indx_cp_l, indx_cp_r,calorically_perfect
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in) :: cv_coeff
    real(rkind) :: get_e_from_temperature
    real(rkind) :: ee
    integer :: l
!
    if (calorically_perfect==1) then
!     ee = cv_coeff(0)*(tt-t0)
      ee = cv_coeff(0)* tt
    else
      ee = cv_coeff(indx_cp_r+1)
      do l=indx_cp_l, indx_cp_r
        if (l==-1) then
          ee = ee+cv_coeff(l)*log(tt/t0)
        else
          ee = ee+cv_coeff(l)/(l+1._rkind)*((tt/t0)**(l+1)-1._rkind)
        endif
      enddo
      ee = ee*t0
    endif
    get_e_from_temperature = ee
  endfunction get_e_from_temperature
!
  subroutine set_fluid_prop(self)
    class(equation_singleideal_object), intent(inout) :: self
!   
    integer :: l, num_coeff_cp
    real(rkind) :: cp_tref, cv_tref, hstar_over_rgas_tref_dim, bcoeff
    real(rkind), allocatable, dimension(:) :: cp_temp
!   
    associate(Prandtl => self%Prandtl, gam => self%gam, visc_model => self%visc_model, &
      sutherland_S => self%sutherland_S, calorically_perfect => self%calorically_perfect, &
      powerlaw_vtexp => self%powerlaw_vtexp, indx_cp_l => self%indx_cp_l, &
      indx_cp_r => self%indx_cp_r, T_ref_dim => self%T_ref_dim, &
      gm => self%gm, gm1 => self%gm1, rfac => self%rfac, rgas0 => self%rgas0, &
      cfg => self%cfg,t0 => self%t0)
!
      call cfg%get("flow","T_ref",T_ref_dim)
      call cfg%get("fluid","Prandtl",Prandtl)
      call cfg%get("fluid","gam",gam)
      call cfg%get("fluid","visc_model",visc_model)
      call cfg%get("fluid","s_suth",sutherland_S)
      call cfg%get("fluid","vt_exp",powerlaw_vtexp)
      call cfg%get("fluid","calorically_perfect",calorically_perfect)
      if (calorically_perfect==0) then
        call cfg%get("fluid","indx_cp_l",indx_cp_l)
        call cfg%get("fluid","indx_cp_r",indx_cp_r)
      elseif (calorically_perfect==1) then
        indx_cp_l = 0
        indx_cp_r = 0
      elseif (calorically_perfect/=1) then
        call fail_input_any("calorically_perfect must be 0 or 1")
      endif
!     
      rfac = Prandtl**(1._rkind/3._rkind)
!     
      allocate(self%cp_coeff(indx_cp_l:indx_cp_r+1))
      allocate(self%cv_coeff(indx_cp_l:indx_cp_r+1))
!     
      self%cp_coeff = 0._rkind
      self%cv_coeff = 0._rkind
      bcoeff = 1._rkind
      if (calorically_perfect==0) then
        call cfg%get("fluid","cp_coeff",cp_temp)
        num_coeff_cp = size(cp_temp)
        if (num_coeff_cp==(indx_cp_r-indx_cp_l+1)) then
          self%cp_coeff(indx_cp_l:indx_cp_r)  = cp_temp(1:num_coeff_cp) ! Assume constant zero
        elseif (num_coeff_cp==(indx_cp_r-indx_cp_l+2)) then
          self%cp_coeff(indx_cp_l:indx_cp_r+1)  = cp_temp(1:num_coeff_cp) ! Read also NASA b1 coefficient
          hstar_over_rgas_tref_dim = self%cp_coeff(indx_cp_r+1)
          do l = indx_cp_l, indx_cp_r
            if (l==-1) then
              hstar_over_rgas_tref_dim = hstar_over_rgas_tref_dim+self%cp_coeff(l)*log(T_ref_dim)
            else
              hstar_over_rgas_tref_dim = hstar_over_rgas_tref_dim+self%cp_coeff(l)/(l+1)*T_ref_dim**(l+1)
            endif
          enddo
          bcoeff = hstar_over_rgas_tref_dim/T_ref_dim
          if (self%masterproc) write(*,*) 'h(T_ref_dim) (kJ/mol):', hstar_over_rgas_tref_dim*8.314510_rkind/1000._rkind
          if (self%masterproc) write(*,*) 'bcoeff:', bcoeff
        else
          call fail_input_any("Error! Check number of cp coefficients")
        endif
        cp_tref = 0._rkind
        do l=indx_cp_l,indx_cp_r
          cp_tref = cp_tref + self%cp_coeff(l)*T_ref_dim**l
        enddo
        cp_tref = cp_tref*rgas0
        cv_tref = cp_tref-rgas0
        gam = cp_tref/cv_tref
        do l=indx_cp_l,indx_cp_r
          self%cp_coeff(l) = self%cp_coeff(l)*T_ref_dim**l
          self%cv_coeff(l) = self%cp_coeff(l)
        enddo
        self%cv_coeff(0)  = self%cp_coeff(0) - 1._rkind
      endif
!     
      gm1 = gam-1._rkind
      gm  = 1._rkind/gm1
      if (self%masterproc) write(*,*) 'Gamma: ', gam
!     
      if (calorically_perfect==1) then
        self%cv_coeff(0) = gm
        self%cp_coeff(0) = gm*gam
      endif
      self%cp_coeff(indx_cp_r+1) =  bcoeff
      self%cv_coeff(indx_cp_r+1) = (bcoeff-1._rkind)
!
      self%cv_coeff = self%cv_coeff*rgas0
      self%cp_coeff = self%cp_coeff*rgas0
!     
    endassociate
!
  endsubroutine set_fluid_prop
!
  subroutine set_flow_params(self)
    class(equation_singleideal_object), intent(inout) :: self
!   
    real(rkind) :: gm1h
    real(rkind) :: c0, e0, etot0, k0
!   
    allocate(self%winf(self%nv))
    allocate(self%winf_past_shock(self%nv))
!   
    associate(Mach => self%Mach, Reynolds => self%Reynolds, theta_wall => self%theta_wall, &
      t0 => self%t0, rho0 => self%rho0, l0 => self%l0, p0 => self%p0, &
      u0 => self%u0, mu0 => self%mu0, &
      powerlaw_vtexp => self%powerlaw_vtexp, T_ref_dim => self%T_ref_dim, &
      rgas0 => self%rgas0, rfac => self%rfac, gm1 => self%gm1, gam => self%gam, &
      T_recovery => self%T_recovery, T_wall => self%T_wall, &
      winf => self%winf, flow_params_cfg => self%flow_params_cfg, &
      indx_cp_l => self%indx_cp_l, indx_cp_r => self%indx_cp_r, Prandtl => self%Prandtl, &
      cv_coeff => self%cv_coeff, calorically_perfect => self%calorically_perfect, &
      cfg => self%cfg, winf_past_shock => self%winf_past_shock, cp_coeff => self%cp_coeff)
!     
      call cfg%get("flow","Reynolds",Reynolds)
      call cfg%get("flow","Mach",Mach)
      call cfg%get("flow","theta_wall",theta_wall)
      if (self%bctags(4)==7) call self%set_oblique_shock()
!
      gm1h       = 0.5_rkind*gm1
      T_recovery = t0 * (1._rkind+gm1h*rfac*Mach**2)
      T_wall     = theta_wall * (T_recovery-t0) + t0
!
      select case(self%flow_init)
      case(0)
        call self%set_chan_prop()
      case(1)
        call self%set_bl_prop()
      end select
!
      p0    = rho0*rgas0*t0
!     
      c0    = sqrt(gam*rgas0*t0)
      u0    = Mach*c0
      e0    = get_e_from_temperature(t0, t0, indx_cp_l, indx_cp_r, cv_coeff, calorically_perfect)
      etot0 = e0+0.5_rkind*u0**2
!
      winf(1) = rho0
      winf(2) = rho0*u0
      winf(3) = 0._rkind
      winf(4) = 0._rkind
      winf(5) = rho0*etot0
!     
      if (abs(u0)<tol_iter) then
        mu0  = c0*rho0*l0/Reynolds ! Reynolds number based on speed of sound
        u0   = c0
        if (self%masterproc) write(*,*) 'Reynolds number based on speed of sound'
      else
        mu0 = u0*rho0*l0/Reynolds
      endif
      k0 = mu0*cp_coeff(0)/Prandtl
!     
      if(self%masterproc) then
        call flow_params_cfg%set("flow_params","u0",             u0)
        call flow_params_cfg%set("flow_params","p0",             p0)
        call flow_params_cfg%set("flow_params","rho0",           rho0)
        call flow_params_cfg%set("flow_params","t0",             t0)
        call flow_params_cfg%set("flow_params","c0",             c0)
        call flow_params_cfg%set("flow_params","l0",             l0)
        call flow_params_cfg%set("flow_params","cp0",            cp_coeff(0))
        call flow_params_cfg%set("flow_params","cv0",            cv_coeff(0))
        call flow_params_cfg%set("flow_params","mu0",            mu0)
        call flow_params_cfg%set("flow_params","k0",             k0)
        call flow_params_cfg%set("flow_params","gam",            gam)
        call flow_params_cfg%set("flow_params","T_wall",         T_wall)
        call flow_params_cfg%set("flow_params","theta_wall",     theta_wall)
        call flow_params_cfg%set("flow_params","rfac",           rfac)
        call flow_params_cfg%set("flow_params","powerlaw_vtexp", powerlaw_vtexp)
        call flow_params_cfg%set("flow_params","T_ref_dim",      T_ref_dim)
        call flow_params_cfg%set("flow_params","Reynolds",       Reynolds)
        call flow_params_cfg%set("flow_params","Mach",           Mach)
        call flow_params_cfg%set("flow_params","Prandtl",        Prandtl)
!
        call flow_params_cfg%write("flow_params.dat")
      endif
!
    endassociate
!
  endsubroutine set_flow_params
!
  subroutine set_bl_prop(self)
    class(equation_singleideal_object), intent(inout) :: self
!   
    real(rkind) :: Re_out, s2tinf, Trat
    real(rkind), allocatable, dimension(:) :: uvec, rhovec, tvec, viscvec, yl0
    real(rkind) :: cf, thrat
    integer :: l, i, m
    logical :: visc_exp
!
    associate(nymax => self%grid%nymax, Reynolds_friction => self%Reynolds_friction, &
      Reynolds => self%Reynolds, l0 => self%l0, &
      Mach => self%Mach, powerlaw_vtexp => self%powerlaw_vtexp, &
      gam => self%gam, rfac => self%rfac, yg => self%grid%yg)
!     
      allocate(uvec(nymax),tvec(nymax),rhovec(nymax),viscvec(nymax),yl0(nymax))
!     
      visc_exp = .false.
      if (self%visc_model == VISC_POWER) visc_exp = .true.
!     
      if (self%bl_case) then
        s2tinf = self%sutherland_S/self%T_ref_dim
        if (.not.self%bl_laminar) then
          if (self%masterproc) write(*,*) 'Input friction Reynolds number: ', Reynolds
          Reynolds_friction  = Reynolds
          Trat = self%T_wall/self%T_recovery
          yl0  = yg(1:nymax)/l0
          call meanvelocity_bl(Reynolds_friction,Mach,Trat,s2tinf,powerlaw_vtexp,visc_exp,gam,rfac,nymax,yl0(1:nymax),uvec,&
            rhovec,tvec,viscvec,Re_out,cf,thrat)
          Reynolds = Re_out
        endif
        if (self%masterproc) write(*,*) 'Reynolds based on free-stream properties: ', Reynolds
      endif
!     
    endassociate
!   
  endsubroutine set_bl_prop
!
  subroutine set_chan_prop(self)
    class(equation_singleideal_object), intent(inout) :: self
!   
    real(rkind) :: s2tinf, Re_out
    real(rkind) :: trec_over_tb, tw_over_tb, tb_over_tw, tw_over_tr
    real(rkind) :: thtmp, rmb, rmtw
    real(rkind) :: gamloc,gamlocold,gamlocgm1
    logical :: visc_exp
    integer :: l, max_iter
!   
    associate(T_wall => self%T_wall, t0 => self%t0, theta_wall => self%theta_wall, &
      T_bulk_target => self%T_bulk_target, &
      indx_cp_l => self%indx_cp_l, indx_cp_r => self%indx_cp_r, &
      cp_coeff => self%cp_coeff, calorically_perfect => self%calorically_perfect, &
      Mach => self%Mach, gam => self%gam, rfac => self%rfac, &
      Reynolds => self%Reynolds, Reynolds_friction => self%Reynolds_friction, &
      powerlaw_vtexp => self%powerlaw_vtexp, gm1 => self%gm1, rgas0 => self%rgas0, &
      yg => self%grid%yg, nymax => self%grid%nymax, yn => self%grid%yn)
!     
      self%T_wall = self%t0 ! Reference temperature in the channel is twall
!     
      if (theta_wall>=-1._rkind) then ! T bulk fixed
        thtmp     = theta_wall
        rmb       = Mach ! Bulk Mach number based on T_bulk (input)
!       
        gamloc    = gam
        max_iter = 50
!       Iterative loop to find gamma at bulk temperature
        do l=1,max_iter
          gamlocold = gamloc
          gamlocgm1 = gamloc-1._rkind
          trec_over_tb  = 1._rkind+gamlocgm1/2._rkind*rfac*rmb**2
          tw_over_tb    = 1._rkind+thtmp*(trec_over_tb-1._rkind)
          tb_over_tw    = 1._rkind/tw_over_tb
          tw_over_tr    = 1._rkind/(trec_over_tb*tb_over_tw)
          T_bulk_target = tb_over_tw*T_wall
          gamloc = get_gamloc(T_bulk_target,t0,indx_cp_l,indx_cp_r,cp_coeff,calorically_perfect,rgas0)
          if (abs(gamloc-gamlocold)<1.D-9) exit
        enddo
        if (self%masterproc) write(*,*) 'Target bulk temperature: ', T_bulk_target
        if (self%masterproc) write(*,*) 'Gamma at bulk temperature: ', gamloc
        rmtw = rmb*sqrt(tb_over_tw) ! Bulk Mach number based on T_wall
        Mach = rmtw
      else
        thtmp = -1._rkind ! For estimation of Re_out
        rmtw  = Mach      ! Bulk Mach number based on T_wall
        rmb   = sqrt(rmtw**2/(1._rkind-gm1/2._rkind*rfac*thtmp*rmtw**2))
        trec_over_tb = 1._rkind+gm1/2._rkind*rfac*rmb**2
        tw_over_tb   = 1._rkind+thtmp*(trec_over_tb-1._rkind)
        tb_over_tw   = 1._rkind/tw_over_tb
        tw_over_tr   = 1._rkind/(trec_over_tb*tb_over_tw)
      endif
!     
      if (self%masterproc) write(*,*) 'Mach number (based on T bulk): ', rmb
      if (self%masterproc) write(*,*) 'Mach number (based on T wall): ', rmtw
!     
      visc_exp = .false.
      if (self%visc_model == VISC_POWER) visc_exp = .true.
!     
      s2tinf = self%sutherland_S/self%T_ref_dim
      Reynolds_friction = Reynolds
      call get_reynolds_cha(Reynolds_friction,rmb,tw_over_tr,s2tinf,powerlaw_vtexp,visc_exp,gam,rfac, &
        nymax/2,yn(1:nymax/2+1),Re_out)
!     
      Reynolds = Re_out
      if (self%masterproc) write(*,*) 'Re bulk (viscosity evaluated at Twall) = ', Reynolds
!     
    endassociate
!
  endsubroutine set_chan_prop
! 
  subroutine set_oblique_shock(self)
    class(equation_singleideal_object), intent(inout) :: self
!   
    integer :: i, m
    real(rkind) :: sins,coss,Mach_normal,Mach_normal2,vel_ratio,rho_ratio,p_ratio,deflection_angle_tan,deflection_angle
    real(rkind) :: Mach_past_shock, Mach_normal_past_shock, t_past_shock, c_past_shock, u_past_shock, v_past_shock, t_ratio
    real(rkind) :: rho_past_shock, e_past_shock, etot_past_shock, xshock_top
    real(rkind) :: xx, tanhf
    integer :: ig_shock
!
    associate(Mach => self%Mach, u0 => self%u0, rho0 => self%rho0, &
      gam => self%gam, gm => self%gm, gm1 => self%gm1, &
      nymax => self%grid%nymax, yg => self%grid%yg, &
      winf => self%winf, indx_cp_l => self%indx_cp_l, indx_cp_r => self%indx_cp_r, &
      xshock_imp => self%xshock_imp, shock_angle => self%shock_angle, winf_past_shock => self%winf_past_shock, &
      nxmax => self%grid%nxmax, tanhfacs => self%tanhfacs, cfg => self%cfg, rgas0 => self%rgas0, &
      calorically_perfect => self%calorically_perfect, t0 => self%t0, l0 => self%l0)
!     
!     Oblique shock deflection is assumed to be negative
!     
      call cfg%get("flow","xshock_imp",xshock_imp)
      xshock_imp = xshock_imp*l0
      call cfg%get("flow","shock_angle",shock_angle)
!     
      shock_angle = shock_angle*pi/180._rkind
      xshock_top  = xshock_imp-yg(nymax)/tan(shock_angle)
      call locateval(self%grid%xg(1:nxmax),nxmax,xshock_top,ig_shock)
      tanhfacs = self%grid%dxg(ig_shock)
      if (self%masterproc) write(*,*) 'xshock_top, tanhfacs:', xshock_top
!     
      sins = sin(shock_angle)
      coss = cos(shock_angle)
      Mach_normal  = Mach*sins
      Mach_normal2 = Mach_normal**2
      Mach_normal_past_shock = (1._rkind+0.5_rkind*gm1*Mach_normal2)/(gam*Mach_normal2-0.5_rkind*gm1)
      Mach_normal_past_shock = sqrt(Mach_normal_past_shock)
      vel_ratio = (2._rkind+gm1*Mach_normal2)/((gam+1._rkind)*Mach_normal2)
      rho_ratio = 1._rkind/vel_ratio
      p_ratio   = 1._rkind+2._rkind*gam/(gam+1._rkind)*(Mach_normal2-1._rkind)
      t_ratio   = p_ratio/rho_ratio
      deflection_angle_tan = 2._rkind*coss/sins*(Mach**2*sins**2-1._rkind)&
      /(2._rkind+Mach**2*(gam+cos(2._rkind*shock_angle)))
      deflection_angle = atan(deflection_angle_tan)
      if (self%masterproc) write(*,*) 'Deflection angle:', deflection_angle*180._rkind/pi
      Mach_past_shock = Mach_normal_past_shock/sin(shock_angle-deflection_angle)
      t_past_shock = t0*t_ratio
      c_past_shock = sqrt(gam*rgas0*t_past_shock)
      u_past_shock = Mach_past_shock*c_past_shock*cos(deflection_angle)
      v_past_shock = Mach_past_shock*c_past_shock*sin(deflection_angle)
      if (self%masterproc) write(*,*) 'Streamwise velocity past shock:', u_past_shock
      if (self%masterproc) write(*,*) 'Vertical velocity past shock:', v_past_shock
      rho_past_shock = rho0*rho_ratio
      e_past_shock    = get_e_from_temperature(t_past_shock, t0, indx_cp_l, indx_cp_r, self%cv_coeff, calorically_perfect)
      etot_past_shock = e_past_shock+0.5_rkind*(u_past_shock**2+v_past_shock**2)
      winf_past_shock(1) =  rho_past_shock
      winf_past_shock(2) =  rho_past_shock*u_past_shock
      winf_past_shock(3) = -rho_past_shock*v_past_shock
      winf_past_shock(4) = 0._rkind
      winf_past_shock(5) =  rho_past_shock*etot_past_shock
      if (self%masterproc) then
        open(18,file='shock_profile.dat')
        do i=1,nxmax
          xx = self%grid%xg(i)-xshock_top
          tanhf = 0.5_rkind*(1._rkind+tanh(xx/tanhfacs))
          write(18,100) self%grid%xg(i), ((winf(m)+tanhf*(winf_past_shock(m)-winf(m))),m=1,5)
          100  format(20ES20.10)
        enddo
        close(18)
      endif
!     
    endassociate
!
  endsubroutine set_oblique_shock
!
  subroutine initial_conditions(self)
    class(equation_singleideal_object), intent(inout) :: self
!   
    integer :: l
!   
    do l=1,self%nv
      self%field%w(l,:,:,:) = self%winf(l)
    enddo
!   ! Temperature initialized so that compute_aux_gpu is able to have the iteration first value (for thermally perfect)
    self%w_aux(:,:,:,6)   = self%t0
!   
    select case (self%flow_init)
    case(-1) ! wind tunnel
      call self%init_wind_tunnel()
    case(0) ! CHANNEL
      call self%init_channel()
    case(1) ! BL, SBLI
      if (self%bl_laminar) then
        call self%init_bl_lam()
      else
        call self%init_bl()
      endif
    case(2) ! TGV
      call self%init_tgv()
    endselect
  endsubroutine initial_conditions
! 
  subroutine get_reynolds_cha(retau,rm,trat,s2tinf,vtexp,visc_exp,gam,rfac,ny,yn,rebulk)
!   
    logical, intent(in) :: visc_exp
    integer, intent(in) :: ny
    real(rkind), intent(in)  :: retau,rm,trat,s2tinf,gam,rfac,vtexp
    real(rkind), intent(out) :: rebulk
    real(rkind), dimension(ny+1), intent(in)  :: yn
!   
    integer :: j,m,itransf,jp99,jj99,jj,l
    real(rkind) :: gm1h, tr, tw, alf, retau_inc, retau_target, retau_correction, retau_musker
    real(rkind) :: retau_inc_target, retau_inc_correction, up99, d99plus_inc
    real(rkind) :: fuu, du, rhow_over_rho, rhow_over_rhojm, rmu_over_rmuw, rmu_over_rmuwjm
    real(rkind) :: ucij, ucijm, res, dy, ycij, ycijm, u99, d99, d99plus, yy, ui
    real(rkind) :: ue, uuu, uuum, thi
    real(rkind) :: rhosum,rhousum,ubulk,ysum
    real(rkind) :: vkc,bcost,c1,b_rich,eta_rich
!   
    real(rkind), dimension(ny+1) :: yplus       ! Incompressible (transformed) y
    real(rkind), dimension(ny+1) :: uplus       ! Incompressible (transformed )velocity profile
    real(rkind), dimension(ny+1) :: ycplus      ! Compressible y
    real(rkind), dimension(ny+1) :: ucplus      ! Compressible velocity profile
    real(rkind), dimension(ny+1) :: yplusold
    real(rkind), dimension(ny+1) :: ucplusold
    real(rkind), dimension(ny+1) :: density     ! Density profile
    real(rkind), dimension(ny+1) :: viscosity   ! Viscosity profile
    real(rkind), dimension(ny+1) :: temperature ! Temperature profile
    real(rkind), dimension(ny+1) :: uu,uci,uc
!   
    itransf = 0 ! Volpiani transformation
    gm1h    = 0.5_rkind*(gam-1._rkind)
    tr   = 1._rkind+gm1h*rfac*rm**2
    tw   = trat*tr
    alf = 0.8259_rkind
    vkc   = 0.39_rkind
    bcost = 5.5_rkind
!   Parameters for Reichardt's inner layer and Finley's wall law
    c1       = -1._rkind/vkc*log(vkc)+bcost
    b_rich   = 0.33_rkind
    eta_rich = 11._rkind
!   
    ycplus = 0._rkind
    do j=2,ny+1
      ycplus(j) = (yn(j)-yn(1))/(yn(ny+1)-yn(1))*retau ! y^+ (compressible)
!     ycplus(j) = (1._rkind+yn(j))*retau ! y^+ (compressible)
    enddo
    yplus    = ycplus  ! Assume as initial guess incompressible (transformed) yplus equal to ycplus
    yplusold = yplus
    do ! outer loop to find yplus
!     
      uplus = 0._rkind
      do j=2,ny+1
        uplus(j) = velpiros_channel(yplus(j),retau)
!       uplus(j) = 5.5_rkind+log(yplus(j))/0.4_rkind
      enddo
!     
      ucplus = uplus ! Assume as initial guess compressible velocity profile equal to incompressible
      do
!       
        ucplusold = ucplus
        do j=1,ny+1
          uu(j)          = ucplus(j)/ucplus(ny)
!         fuu            = alf*uu(j)+(1._rkind-alf)*uu(j)**2 ! Duan and Martin (see Zhang JFM 2014)
          fuu            = uu(j)                             ! Walz
          temperature(j) = tw+(tr-tw)*fuu+(1._rkind-tr)*uu(j)**2
          density(j)     = 1._rkind/temperature(j)
          if (visc_exp) then
            viscosity(j) = temperature(j)**vtexp ! Power-law
          else
            viscosity(j)   = temperature(j)**(1.5_rkind)*(1._rkind+s2tinf)/(temperature(j)+s2tinf)
          endif
        enddo
!       
!       Inverse Van Driest or Volpiani transformation for u
!       
        do j=2,ny+1
          du              = uplus(j)-uplus(j-1)
          rhow_over_rho   = density(1)/density(j)
          rhow_over_rhojm = density(1)/density(j-1)
          rmu_over_rmuw   = viscosity(j)/viscosity(1)
          rmu_over_rmuwjm = viscosity(j-1)/viscosity(1)
          if (itransf==1) then
!           VD
            ucij  = sqrt(rhow_over_rho)
            ucijm = sqrt(rhow_over_rhojm)
          else
!           Volpiani
            ucij  = sqrt(rhow_over_rho)*sqrt(rmu_over_rmuw)
            ucijm = sqrt(rhow_over_rhojm)*sqrt(rmu_over_rmuwjm)
          endif
!         
          ucplus(j) = ucplus(j-1)+0.5_rkind*(ucij+ucijm)*du
!         
        enddo
!       
        res = 0._rkind
        do j=1,ny+1
          res = res+abs(ucplus(j)-ucplusold(j))
        enddo
        res = res/(ny+1)
        if (res<tol_iter) exit
!       
      enddo ! End of iterative loop to find ucplus
!     
!     At this stage we have ucplus and density, temperature, viscosity distributions
!     
!     Compute incompressible yplus with transformation
!     
      yplus = 0._rkind
      do j=2,ny+1
        dy = ycplus(j)-ycplus(j-1)
!       
        rhow_over_rho   = density(1)/density(j)
        rhow_over_rhojm = density(1)/density(j-1)
        rmu_over_rmuw   = viscosity(j)/viscosity(1)
        rmu_over_rmuwjm = viscosity(j-1)/viscosity(1)
!       
!       ! VD
        if (itransf==1) then
          ycij  = 1._rkind
          ycijm = 1._rkind
        else
!         ! Volpiani
          ycij  = sqrt(rhow_over_rho)*sqrt(rmu_over_rmuw**3)
          ycij  = 1._rkind/ycij
          ycijm = sqrt(rhow_over_rhojm)*sqrt(rmu_over_rmuwjm**3)
          ycijm  = 1._rkind/ycijm
        endif
!       
        yplus(j) = yplus(j-1)+0.5_rkind*(ycij+ycijm)*dy
      enddo
!     
      res  = 0._rkind
      do j=1,ny+1
        res  = abs(yplus(j)-yplusold(j))
      enddo
      res  = res/(ny+1)
      if (res<tol_iter) exit
      yplusold = yplus
    enddo
!   
    rhosum   = 0._rkind
    rhousum  = 0._rkind
    ysum     = 0._rkind
    do j=1,ny
      dy      = yn(j+1)-yn(j)
      ysum    = ysum+dy
      rhosum  = rhosum+0.5*(density(j)+density(j+1))*dy
      rhousum = rhousum+0.5*(density(j)*ucplus(j)+density(j+1)*ucplus(j+1))*dy
    enddo
    rhosum  = rhosum/ysum
    rhousum = rhousum/ysum
    ubulk   = rhousum/rhosum
    rebulk  = retau*ubulk*rhosum/density(1)
!   
    return
  end subroutine get_reynolds_cha
! 
  subroutine init_bl_lam(self)
    class(equation_singleideal_object), intent(inout) :: self
!
    integer :: i,j,k,n,nst,ii
    integer :: ne
    real(rkind) :: etad, deta, Tinf_dim, Twall_dim, etaedge, xbl
    real(rkind) :: xx, yy, etast, wl, wr, ust, vst, tst, rst
    real(rkind) :: rho,uu,vv,ww,rhouu,rhovv,rhoww,ee,tt
    real(rkind), allocatable, dimension(:) :: ubl,tbl,vbl,eta
!
    allocate(self%wmean(1-self%grid%ng:self%field%nx+self%grid%ng+1, 1:self%field%ny, 4))
!
    associate(nx => self%field%nx, ny => self%field%ny, nz => self%field%nz, ng => self%grid%ng, &
      x => self%field%x, y => self%field%y, z => self%field%z, Reynolds => self%Reynolds, &
      xg => self%grid%xg, nxmax => self%grid%nxmax, &
      rho0 => self%rho0, u0 => self%u0, p0 => self%p0, gm => self%gm, &
      wmean => self%wmean, w => self%field%w, Mach => self%Mach, gam => self%gam, &
      rfac => self%rfac, Prandtl => self%Prandtl, w_aux => self%w_aux,  &
      deltavec => self%deltavec ,deltavvec => self%deltavvec, cfvec => self%cfvec, &
      t0 => self%t0, T_ref_dim => self%T_ref_dim, T_wall => self%T_wall, l0 => self%l0, &
      indx_cp_l => self%indx_cp_l, indx_cp_r => self%indx_cp_r, cv_coeff => self%cv_coeff, &
      calorically_perfect => self%calorically_perfect, rgas0 => self%rgas0)
!     
      etad = 20._rkind
      deta = 0.02_rkind
      ne = nint(etad/deta)
!     
      allocate(ubl(0:ne))
      allocate(vbl(0:ne))
      allocate(tbl(0:ne))
      allocate(eta(0:ne))
!     
      Tinf_dim  = t0*T_ref_dim/t0
      Twall_dim = T_wall*T_ref_dim/t0
!     
      call compressible_blasius(ne,etad,deta,gam,Mach,Prandtl,Tinf_dim,Twall_dim,eta,ubl,vbl,tbl,etaedge)
!     
      xbl = Reynolds/(etaedge**2)*l0
      if (self%masterproc) write(*,*) 'XBL', xbl
!     
      wmean = 0._rkind
      do i=1-ng,nx+ng+1
        ii = self%field%ncoords(1)*nx+i
        do j=1,ny
          xx = xg(ii)+xbl
          yy = y(j)
          etast = (yy/l0)*etaedge*sqrt(xbl/xx)
          nst = 1
          do n=1,ne-1
            if (eta(n+1)<etast) nst = n
          enddo
          wl = etast - eta(nst)
          wr = eta(nst+1) - etast
          ust = (wr*ubl(nst)+wl*ubl(nst+1))/deta * u0
          vst = (wr*vbl(nst)+wl*vbl(nst+1))/deta/sqrt(Reynolds*xx/l0)*u0
          tst = (wr*tbl(nst)+wl*tbl(nst+1))/deta/Tinf_dim*t0
          rst = p0/rgas0/tst
          wmean(i,j,1) = rst
          wmean(i,j,2) = rst*ust
          wmean(i,j,3) = rst*vst
        enddo
      enddo
!     
      do k=1,nz
        do j=1,ny
          do i=1,nx
            rho = wmean(i,j,1)
            uu  = wmean(i,j,2)/rho
            vv  = wmean(i,j,3)/rho
            ww  = wmean(i,j,4)/rho
            rhouu = rho*uu
            rhovv = rho*vv
            rhoww = rho*ww
            w(1,i,j,k) = rho
            w(2,i,j,k) = rhouu
            w(3,i,j,k) = rhovv
            w(4,i,j,k) = rhoww
            tt         = p0/rgas0/rho
            w_aux(i,j,k,6) = tt
            ee = get_e_from_temperature(tt, t0, indx_cp_l, indx_cp_r, cv_coeff, calorically_perfect)
            w(5,i,j,k) = rho*ee + 0.5_rkind*(rhouu**2+rhovv**2+rhoww**2)/rho
          enddo
        enddo
      enddo
!
    endassociate
  endsubroutine init_bl_lam
! 
  subroutine init_bl(self)
    class(equation_singleideal_object), intent(inout) :: self
!
    real(rkind), dimension(1-self%grid%ng:self%grid%nxmax+self%grid%ng+1) :: thvec
    real(rkind), dimension(self%field%ny) :: uvec, rhovec, tvec, viscvec, yl0
    real(rkind), dimension(3) :: rr
    real(rkind) :: Trat,thrat,cf
    real(rkind) :: retau,retauold,delta,deltav,th
    real(rkind) :: rhoi,ui,ti,yl
    real(rkind) :: vi,vi_j,vi_jm
    real(rkind) :: rho,uu,vv,ww,rhouu,rhovv,rhoww,ee,tt
    real(rkind) :: u0_02
    real(rkind) :: s2tinf, vtexp, redelta
    integer :: i,j,k,ii,jj,jjj,m
    logical :: visc_exp
    logical :: file_exists
!
!   ghost only on x dir
    allocate(self%wmean(1-self%grid%ng:self%field%nx+self%grid%ng+1, 1:self%field%ny, 4))
    allocate(self%deltavec(1-self%grid%ng:self%grid%nxmax+self%grid%ng+1))
    allocate(self%deltavvec(1-self%grid%ng:self%grid%nxmax+self%grid%ng+1))
    allocate(self%cfvec(1-self%grid%ng:self%grid%nxmax+self%grid%ng+1))
!
    associate(nx => self%field%nx, ny => self%field%ny, nz => self%field%nz, ng => self%grid%ng, &
      y => self%field%y, z => self%field%z, Reynolds_friction => self%Reynolds_friction, &
      xg => self%grid%xg, nxmax => self%grid%nxmax, &
      rho0 => self%rho0, u0 => self%u0, p0 => self%p0, gm => self%gm, l0 => self%l0, &
      wmean => self%wmean, w => self%field%w, Mach => self%Mach, gam => self%gam, &
      rfac => self%rfac, Prandtl => self%Prandtl, w_aux => self%w_aux,  &
      deltavec => self%deltavec ,deltavvec => self%deltavvec, cfvec => self%cfvec, &
      indx_cp_l => self%indx_cp_l, indx_cp_r => self%indx_cp_r, cv_coeff => self%cv_coeff, &
      calorically_perfect => self%calorically_perfect, rgas0 => self%rgas0, t0 => self%t0)
!
      Trat = self%T_wall/self%T_recovery
!
      deltavec(1)  = l0
      deltavvec(1) = l0/Reynolds_friction
      yl0          = y(1:ny)/l0
!
      s2tinf = self%sutherland_S/self%T_ref_dim
      vtexp  = self%powerlaw_vtexp
      visc_exp = .false.
      if (self%visc_model == VISC_POWER) visc_exp = .true.
!     
      inquire(file='blvec.bin',exist=file_exists)
      if (file_exists) then
        open(183,file='blvec.bin',form='unformatted')
        read(183) cfvec,thvec,deltavec,deltavvec
        close(183)
      else
        call meanvelocity_bl(Reynolds_friction,Mach,Trat,s2tinf,vtexp,visc_exp,gam,rfac,ny,yl0(1:ny),uvec,&
          rhovec,tvec,viscvec,redelta,cf,thrat)
        cfvec(1)     = cf
        thvec(1)     = thrat*deltavec(1)
        retauold     = Reynolds_friction
!
        do i=1,ng
          retau = retauold
          do
            call meanvelocity_bl(retau,Mach,Trat,s2tinf,vtexp,visc_exp,gam,rfac,ny,yl0(1:ny),uvec,&
              rhovec,tvec,viscvec,redelta,cf,thrat)
            th = thvec(2-i)-0.25_rkind*abs((xg(1-i)-xg(2-i)))*(cf+cfvec(2-i)) ! dthdx evaluated with second order accuracy
            delta  = th/thrat
            deltav = deltavvec(2-i)*sqrt(cfvec(2-i)/cf)
            retau  = delta/deltav
            if (abs(retau-retauold) < 0.01_rkind) exit
            retauold = retau
          enddo
          thvec    (1-i) = th
          cfvec    (1-i) = cf
          deltavec (1-i) = delta
          deltavvec(1-i) = deltav
        enddo
!
        retauold = Reynolds_friction
!
        if (self%masterproc) open(182,file='cfstart.dat')
        do i=2,nxmax+ng+1
          retau = retauold
          do
            call meanvelocity_bl(retau,Mach,Trat,s2tinf,vtexp,visc_exp,gam,rfac,ny,yl0(1:ny),uvec,&
              rhovec,tvec,viscvec,redelta,cf,thrat)
            th = thvec(i-1)+0.25_rkind*(xg(i)-xg(i-1))*(cf+cfvec(i-1))
            delta  = th/thrat
            deltav = deltavvec(i-1)*sqrt(cfvec(i-1)/cf)
            retau  = delta/deltav
            if (abs(retau-retauold)<0.01_rkind) exit
            retauold = retau
          enddo
          thvec    (i) = th
          cfvec    (i) = cf
          deltavec (i) = delta
          deltavvec(i) = deltav
          if (self%masterproc) write(182,100) xg(i),delta,deltav,cf,th
          100  format(20ES20.10)
        enddo
        if (self%masterproc) close(182)
        if (self%masterproc) then
          open(183,file='blvec.bin',form='unformatted')
          write(183) cfvec,thvec,deltavec,deltavvec
          close(183)
        endif
      endif
!     
!     Compute locally wmean from 1-ng to nx+ng+1
!     
      wmean = 0._rkind
      do i=1-ng,nx+ng+1
        ii = self%field%ncoords(1)*nx+i
        delta  = deltavec(ii)
        deltav = deltavvec(ii)
        retau  = delta/deltav
        call meanvelocity_bl(retau,Mach,Trat,s2tinf,vtexp,visc_exp,gam,rfac,ny,yl0(1:ny),uvec,&
          rhovec,tvec,viscvec,redelta,cf,thrat)
        do j=1,ny
          yl = y(j)/delta
          call locateval(yl0(1:ny),ny,yl,jj) ! yl is between yvec(jj) and yvec(jj+1)
          m = 4
          jjj = min(max(jj-(m-1)/2,1),ny+1-m)
          call pol_int(yl0(jjj),rhovec(jjj),m,yl,rhoi)
          call pol_int(yl0(jjj),uvec(jjj),m,yl,ui  )
          call pol_int(yl0(jjj),tvec(jjj),m,yl,ti  )
          wmean(i,j,1) = rhoi*self%rho0
          wmean(i,j,2) = rhoi*ui*self%u0*self%rho0
        enddo
      enddo
!     
      do i=1-ng,nx+ng
        ii = self%field%ncoords(1)*nx+i
        do j=2,ny
          vi_j  = -(wmean(i+1,j,2)-wmean(i,j,2))/(xg(ii+1)-xg(ii))
          vi_jm = -(wmean(i+1,j-1,2)-wmean(i,j-1,2))/(xg(ii+1)-xg(ii))
          vi    = 0.5_rkind*(vi_j+vi_jm)
          wmean(i,j,3) = wmean(i,j-1,3)+vi*(y(j)-y(j-1))
        enddo
      enddo
!
      u0_02 = 0._rkind !0.02_rkind*u0
      do k=1,nz
        do j=1,ny
          do i=1,nx
            rho = wmean(i,j,1)
            call get_crandom_f(rr)
            rr = rr-0.5_rkind
            uu  = wmean(i,j,2)/rho+u0_02*rr(1)
            vv  = wmean(i,j,3)/rho+u0_02*rr(2)
            ww  = wmean(i,j,4)/rho+u0_02*rr(3)
!
            rhouu = rho*uu
            rhovv = rho*vv
            rhoww = rho*ww
            w(1,i,j,k) = rho
            w(2,i,j,k) = rhouu
            w(3,i,j,k) = rhovv
            w(4,i,j,k) = rhoww
            tt         = p0/rgas0/rho
            w_aux(i,j,k,6) = tt
            ee = get_e_from_temperature(tt, t0, indx_cp_l, indx_cp_r, cv_coeff, calorically_perfect)
            w(5,i,j,k) = rho*ee + 0.5_rkind*(rhouu**2+rhovv**2+rhoww**2)/rho
          enddo
        enddo
      enddo
!
      call self%add_synthetic_perturbations()
!
    endassociate
  endsubroutine init_bl
!
  subroutine add_synthetic_perturbations(self)
    class(equation_singleideal_object), intent(inout) :: self
    real(rkind), dimension(:,:), allocatable :: synth_params
    real(rkind) :: rho_wall, tau_wall, u_tau, lz_plus, rr, rhofac, up, vp, wp
    real(rkind) :: arg_sin, arg_cos, ys , ufy , vfy , dup_dx , dvfy_dy, dvp_dy
    real(rkind) :: rho,ee,tt
    real(rkind) :: delta
    real(rkind), dimension(3) :: rr3
    integer :: i,j,k,l,n_streaks,ii
!
    allocate(synth_params(5,7))
!
    synth_params(:,1) = [  12._rkind,   0.3_rkind,  0.45_rkind,   0.6_rkind,   0.5_rkind]
    synth_params(:,2) = [  1.2_rkind,   0.3_rkind,   0.2_rkind,  0.08_rkind,  0.04_rkind]
    synth_params(:,3) = [-0.25_rkind, -0.06_rkind, -0.05_rkind, -0.04_rkind, -0.03_rkind]
    synth_params(:,4) = [ 0.12_rkind,   1.2_rkind,   0.6_rkind,   0.4_rkind,   0.2_rkind]
    synth_params(:,5) = [ 10.0_rkind,   0.9_rkind,   0.9_rkind,   0.9_rkind,   0.9_rkind]
    synth_params(:,6) = [120.0_rkind, 0.333_rkind,  0.25_rkind,   0.2_rkind, 0.166_rkind]
    synth_params(:,7) = [  0.0_rkind,    1._rkind,    1._rkind,    1._rkind,    1._rkind]
!
    associate(nx => self%field%nx, ny => self%field%ny, nz => self%field%nz, ng => self%grid%ng, &
      deltavec => self%deltavec, deltavvec => self%deltavvec, wmean => self%wmean, cfvec => self%cfvec, &
      rho0 => self%rho0, u0 => self%u0, lz => self%grid%domain_size(3), w => self%field%w, l0 => self%l0, &
      x => self%field%x, y => self%field%y, z => self%field%z, p0 => self%p0, gm => self%gm, &
      w_aux => self%w_aux, indx_cp_l => self%indx_cp_l, indx_cp_r => self%indx_cp_r, cv_coeff => self%cv_coeff, &
      calorically_perfect => self%calorically_perfect, rgas0 => self%rgas0, t0 => self%t0)
!
      rho_wall = wmean(1,1,1) ! inflow, wall, density
      tau_wall = cfvec(1)*rho0*u0**2*0.5_rkind
      u_tau    = sqrt(tau_wall/rho_wall)
!
      synth_params(1,1) = synth_params(1,1) * deltavvec(1)
      synth_params(1,4) = synth_params(1,4) * u_tau / deltavvec(1)
      synth_params(1,5) = synth_params(1,5) * u_tau
!
      lz_plus = lz/deltavvec(1)
      n_streaks = nint(lz_plus / synth_params(1,6))
      synth_params(1,6) = lz / n_streaks
!
      do l=2,5
        synth_params(l,1) = synth_params(l,1) * deltavec(1)
        synth_params(l,4) = synth_params(l,4) * u0 / deltavec(1)
        synth_params(l,5) = synth_params(l,5) * u0
        synth_params(l,6) = synth_params(l,6) * lz
        call get_crandom_f(rr)
!       rr = 0._rkind
        synth_params(l,7) = synth_params(l,7) * 2._rkind*pi *rr
      enddo
!     random must be synced across processes
      call mpi_bcast(synth_params(2:5,7),4,mpi_prec,0,self%field%mp_cart,self%mpi_err)
!
      do k=1,nz
        do j=1,ny
          do i=1,nx
!           
            rhofac  = wmean(i,1,1)/wmean(i,j,1)
            rhofac  = sqrt(rhofac)*u0
            up      = 0._rkind
            vp      = 0._rkind
            wp      = 0._rkind
            do l=1,5
              arg_sin = synth_params(l,4)*x(i)/synth_params(l,5)
              arg_cos = 2._rkind*pi*z(k)/synth_params(l,6) + synth_params(l,7)
              ys      = y(j)/synth_params(l,1)
              ufy     = ys    * exp(-ys)
              vfy     = ys**2 * exp(-(ys**2))
              up      = up+synth_params(l,2) * ufy * sin(arg_sin) * cos(arg_cos)
              vp      = vp+synth_params(l,3) * vfy * sin(arg_sin) * cos(arg_cos)
!
              dup_dx  = synth_params(l,2) * ufy * synth_params(l,4)/synth_params(l,5) * cos(arg_sin)
              dvfy_dy = 2._rkind*ys/synth_params(l,1)*exp(-(ys**2))*(1._rkind-ys**2)
              dvp_dy  = synth_params(l,3) * dvfy_dy * sin(arg_sin)
              wp      = wp-(dup_dx+dvp_dy)*sin(arg_cos)*synth_params(l,6)/(2._rkind*pi)
            enddo
            up = up * rhofac
            vp = vp * rhofac
            wp = wp * rhofac
!           
            ii     = self%field%ncoords(1)*nx+i
            delta  = deltavec(ii)
            if (y(j)<delta) then
              call get_crandom_f(rr3)
              rr3 = rr3-0.5_rkind
              up = up+0.03_rkind*u0*rr3(1)*(y(j)/l0)
              vp = vp+0.03_rkind*u0*rr3(2)*(y(j)/l0)
              wp = wp+0.03_rkind*u0*rr3(3)*(y(j)/l0)
            endif
!           
            rho        = w(1,i,j,k)
            w(2,i,j,k) = w(2,i,j,k) + up*w(1,i,j,k)
            w(3,i,j,k) = w(3,i,j,k) + vp*w(1,i,j,k)
            w(4,i,j,k) = w(4,i,j,k) + wp*w(1,i,j,k)
            tt         = p0/rgas0/rho
            w_aux(i,j,k,6) = tt
            ee = get_e_from_temperature(tt, t0, indx_cp_l, indx_cp_r, cv_coeff, calorically_perfect)
            w(5,i,j,k) = rho*ee + 0.5_rkind*(w(2,i,j,k)**2 + w(3,i,j,k)**2 + w(4,i,j,k)**2)/w(1,i,j,k)
!           
          enddo
        enddo
      enddo
!
    endassociate
  endsubroutine add_synthetic_perturbations
! 
  function velmusker(yplus,retau)
    real(rkind), intent(in) :: yplus,retau
    real(rkind) :: yp,eta,pi_wake
    real(rkind) :: velmusker
    yp        = yplus
    eta       = yp/retau
    eta       = min(1._rkind,eta)
    yp        = eta*retau
    pi_wake   = 0.434_rkind
    velmusker = 5.424_rkind*atan((2*yp-8.15_rkind)/16.7_rkind)+ &
    log10((yp+10.6_rkind)**9.6_rkind/(yp**2-8.15_rkind*yp+86)**2)-3.51132976630723_rkind &
    +2.44_rkind*(pi_wake*(6*eta**2-4*eta**3)+(eta**2*(1-eta)))
  end function velmusker
! 
  function velpiros_channel(yplus,retau)
    real(rkind), intent(in) :: yplus,retau
    real(rkind) :: velpiros_channel
    real(rkind) :: yp,eta,vkc,vkci,c,c1,pi_wake,etastar,b,b_rich,eta_rich,uinn,uout
!   
    yp  = yplus
    eta = yp/retau
    eta = min(1._rkind,eta)
    yp  = eta*retau
!   
!   c       = 4.17_rkind
!   vkc     = 0.383_rkind
    c       = 5.3_rkind
    vkc     = 0.39_rkind
    vkci    = 1._rkind/vkc
    c1      = -1._rkind/vkc*log(vkc)+c
    pi_wake = 0.00857_rkind
    etastar = 0.275_rkind
    b       = 0.0447_rkind
    b_rich   = 0.33_rkind
    eta_rich = 11._rkind
!   
!   uinn    = c+log(yp)*vkci
    uinn    = vkci*log(1._rkind+vkc*yp)+c1*(1._rkind-exp(-yp/eta_rich)-yp/eta_rich*exp(-b_rich*yp))
    uout    = c+vkci*log(retau)+vkci*log(etastar)
    uout    = uout+(1-etastar)/(2*vkc*etastar)*(1-(1-eta)**2/(1-etastar)**2)
    if (eta.le.etastar) then
      velpiros_channel = uinn
    else
      velpiros_channel = uout
    endif
!   
    return
  end function velpiros_channel
! 
  subroutine compressible_blasius(n,etad,deta,gam,Mach,Prandtl,Te,Twall,eta,u,v,tbl,delta1)
    integer, intent(in) :: n
    real(rkind), intent(in) :: etad, deta, Mach, Prandtl, Te, Twall, gam
    real(rkind), intent(out) :: delta1
    real(rkind), dimension(0:n), intent(inout) :: u,v,tbl,eta
    real(rkind), dimension(0:n) :: f,t,g,h,a1,a2,a3,a4,a5,a6,a7,a8,s,r,visc,tcr
!   
    integer :: iflag,nm,i,j,kk
    real(rkind) :: Rg,s_suth,eps,z1,z2,z3,gm,adiab,T0e,Ae,Ue,Cp,He,H0e,visce,s0,r0,tt,vis
    real(rkind) :: u0,f0,t0,g0,h0
    real(rkind) :: u1,f1,t1,g1,h1
    real(rkind) :: u2,f2,t2,g2,h2
    real(rkind) :: u3,f3,t3,g3,h3
    real(rkind) :: a10,a20,a30,a40,a50,a60,a70,a80
    real(rkind) :: a11,a21,a31,a41,a51,a61,a71,a81
    real(rkind) :: a12,a22,a32,a42,a52,a62,a72,a82
    real(rkind) :: a13,a23,a33,a43,a53,a63,a73,a83
    real(rkind) :: eq1,eq2,b1,b2,b3,b4,det,c1,c2,c3,c4,da,db,dd,rho,delta2
!   
    Rg = 287.15_rkind
    s_suth = 110.4_rkind
    eps    = 0.000001_rkind
    z1     = 0.334_rkind
    z2     = 0.82_rkind
    z3     = 0.22_rkind
    iflag  =  1
!   
    nm = n-1
!   
!   *************************************  fluid conditions
!   
    gm    = gam / ( gam -1._rkind)
    adiab = (1._rkind+ (gam - 1._rkind) / 2._rkind *Mach*Mach)
    T0e   = Te * adiab
    Ae    = sqrt (gam * Rg * Te)
    Ue    = Mach * ae
    Cp    = gm * Rg
    He    = Cp * Te
    H0e   = Cp * T0e
!   
    if (Te.lt.110.4_rkind) then
      visce = .693873D-6*te
    else
      visce = 1.458D-5 * te**1.5_rkind/( te + s_suth )
    end if
!   
!   ************************************* initial conditions
!   
    j = 0
    eta(0)  = 0._rkind
    u(0)    = 0._rkind
    h(0)    = 0._rkind
    a1(0)   = 0._rkind
    a2(0)   = 0._rkind
    a3(0)   = 0._rkind
    a5(0)   = 1._rkind
    a6(0)   = 0._rkind
    a7(0)   = 0._rkind
    s0      = z1
    if(iflag.eq.0) then
      g(0)  = 0._rkind
      a4(0) = 1._rkind
      a8(0) = 0._rkind
      r0    = z2
    else
      t(0)  = (Twall - Te )/(T0e-te)
      a4(0) = 0._rkind
      a8(0) = 1._rkind
      r0    = z3
    end if
!   
    do
!
      f(0) = s0
      if ( iflag.eq.0)  t(0) = r0
      if ( iflag.eq.1)  g(0) = r0
!     
      Tbl(0)  = Te + (T0e - Te) * t(0)
      tt      = tbl(0)
      if (tt<110.4_rkind) then
        visc(0) = .693873D-6*tt
      else
        visc(0) = 1.458D-5 * tt**1.5_rkind  / ( tt + s_suth )
      end if
      vis     = visc(0) / visce
!     
!     ******************************** Runge-Kutta integration
!     
      do i = 0,nm
!       
        u0 =   f(i) / vis * deta
        f0 = - h(i) * f(i) / vis* deta
        t0 =   g(i) * Prandtl / vis * deta
        g0 = -(h(i) * g(i) * Prandtl + 2._rkind * f(i) **2._rkind) / vis * deta
        h0 =  .5_rkind*u(i) / (1 + (T0e/Te - 1) * t(i)) * deta
        u1 =   (f(i) + .5_rkind*f0) / vis * deta
        f1 = - (h(i) + .5_rkind*h0) * (f(i) + .5_rkind*f0) /vis * deta
        t1 =   (g(i) + .5_rkind*g0) * Prandtl / vis * deta
        g1 = -((h(i) + .5_rkind*h0) * (g(i) + .5_rkind*g0) * Prandtl + 2._rkind * (f(i) + .5_rkind*f0)**2._rkind) /vis*deta
        h1 =  .5_rkind*(u(i)+.5_rkind*u0) / (1 + (T0e/Te - 1) * (t(i)+.5_rkind*t0)) *deta
        u2 =   (f(i) + .5_rkind*f1) / vis * deta
        f2 = - (h(i) + .5_rkind*h1) * (f(i) + .5_rkind*f1) /vis * deta
        t2 =   (g(i) + .5_rkind*g1) * Prandtl / vis * deta
        g2 = -((h(i) + .5_rkind*h1) * (g(i) + .5_rkind*g1) * Prandtl + 2._rkind * (f(i) + .5_rkind*f1)**2._rkind) /vis*deta
        h2 =  .5_rkind*(u(i)+.5_rkind*u1) / (1 + (T0e/Te - 1) * (t(i)+.5*t1)) *deta
        u3 =   (f(i) + f2) / vis * deta
        f3 = - (h(i) + h2) * (f(i) + f2) /vis * deta
        t3 =   (g(i) + g2) * Prandtl / vis * deta
        g3 = -((h(i) + h2) * (g(i) + g2) * Prandtl + 2._rkind * (f(i) + f2)**2._rkind) /vis*deta
        h3 = .5_rkind*(u(i)+u2) / (1 + (T0e/Te - 1) *(t(i)+t2))*deta
!       
        a10 =   a5(i)/vis *deta
        a20 =   a6(i) / vis *deta
        a30 =   a7(i) * Prandtl / vis *deta
        a40 =   a8(i) *Prandtl / vis *deta
        a50 = - a5(i) * h(i) / vis *deta
        a60 = - a6(i) * h(i) / vis *deta
        a70 = -(4*f(i)*a5(i)+ Prandtl * h(i) *a7(i))/vis *deta
        a80 = -(4*f(i)*a6(i)+ Prandtl * h(i)*a8(i))/vis *deta
!       
        a11 =   (a5(i) + .5_rkind*a50) / vis *deta
        a21 =   (a6(i) + .5_rkind*a60) / vis *deta
        a31 =   (a7(i) + .5_rkind*a70) * Prandtl / vis *deta
        a41 =   (a8(i) + .5_rkind*a80) * Prandtl / vis *deta
        a51 = - (a5(i) + .5_rkind*a50) * (h(i) + .5_rkind*h0) / vis *deta
        a61 = - (a6(i) + .5_rkind*a60) * (h(i) + .5_rkind*h0) / vis *deta
        a71 = -(4* (f(i) + .5_rkind*f0) * (a5(i) + .5_rkind*a50)+Prandtl * (h(i) + .5_rkind*h0) *(a7(i) + .5_rkind*a70))/vis *deta
        a81 = -(4* (f(i) + .5_rkind*f0) * (a6(i) + .5_rkind*a60)+Prandtl * (h(i) + .5_rkind*h0) *(a8(i) + .5_rkind*a80))/vis *deta
!       
        a12 =   (a5(i) + .5_rkind*a51) / vis *deta
        a22 =   (a6(i) + .5_rkind*a61) / vis *deta
        a32 =   (a7(i) + .5_rkind*a71) * Prandtl / vis *deta
        a42 =   (a8(i) + .5_rkind*a81) * Prandtl / vis *deta
        a52 = - (a5(i) + .5_rkind*a51) * (h(i) + .5_rkind*h1) / vis *deta
        a62 = - (a6(i) + .5_rkind*a61) * (h(i) + .5_rkind*h1) / vis *deta
        a72 = -(4* (f(i) + .5_rkind*f1) * (a5(i) + .5_rkind*a51)+Prandtl * (h(i) + .5_rkind*h1) *(a7(i) + .5_rkind*a71))/vis *deta
        a82 = -(4* (f(i) + .5_rkind*f1) * (a6(i) + .5_rkind*a61)+Prandtl * (h(i) + .5_rkind*h1) *(a8(i) + .5_rkind*a81))/vis *deta
!       
        a13 =   (a5(i) + a52) / vis *deta
        a23 =   (a6(i) + a62) / vis *deta
        a33 =   (a7(i) + a72) * Prandtl / vis *deta
        a43 =   (a8(i) + a82) * Prandtl / vis *deta
        a53 = - (a5(i) + a52) * (h(i) + .5_rkind*h2) / vis *deta
        a63 = - (a6(i) + a62) * (h(i) + .5_rkind*h2) / vis *deta
        a73 = -(4* (f(i) + f2) * (a5(i) + a52) + Prandtl * (h(i) + h2) *(a7(i) + a72))/vis *deta
        a83 = -(4* (f(i) + f2) * (a6(i) + a62) + Prandtl * (h(i) + h2) *(a8(i) + a82))/vis *deta
!       
        f(i+1) = f(i) + (f0 + 2._rkind*f1 + 2*f2 + f3) / 6._rkind
        u(i+1) = u(i) + (u0 + 2._rkind*u1 + 2*u2 + u3) / 6._rkind
        t(i+1) = t(i) + (t0 + 2._rkind*t1 + 2*t2 + t3) / 6._rkind
        g(i+1) = g(i) + (g0 + 2._rkind*g1 + 2*g2 + g3) / 6._rkind
        h(i+1) = h(i) + (h0 + 2._rkind*h1 + 2*h2 + h3) / 6._rkind
!       
        a1(i+1) = a1(i) + (a10 + 2._rkind*a11 + 2._rkind*a12 + a13) / 6._rkind
        a2(i+1) = a2(i) + (a20 + 2._rkind*a21 + 2._rkind*a22 + a23) / 6._rkind
        a3(i+1) = a3(i) + (a30 + 2._rkind*a31 + 2._rkind*a32 + a33) / 6._rkind
        a4(i+1) = a4(i) + (a40 + 2._rkind*a41 + 2._rkind*a42 + a43) / 6._rkind
        a5(i+1) = a5(i) + (a50 + 2._rkind*a51 + 2._rkind*a52 + a53) / 6._rkind
        a6(i+1) = a6(i) + (a60 + 2._rkind*a61 + 2._rkind*a62 + a63) / 6._rkind
        a7(i+1) = a7(i) + (a70 + 2._rkind*a71 + 2._rkind*a72 + a73) / 6._rkind
        a8(i+1) = a8(i) + (a80 + 2._rkind*a81 + 2._rkind*a82 + a83) / 6._rkind
        eta(i+1)  = eta(i) + deta
!       
!       **************************************   new value of visc
!       
        Tbl(i+1)  = Te + (T0e - Te) * t(i+1)
        tt      = Tbl(i+1)
        if (tt<110.4_rkind) then
          visc(i+1) = .693873D-6*tt
        else
          visc(i+1) = 1.458D-5 * tt**1.5_rkind  / ( tt + s_suth )
        end if
        vis     = visc(i+1) / visce
!       
      end do
!     
!     ******************************************* shooting method
!     ******************************************* with Newton Raphson
!     
      eq1 = 1._rkind - u(n)
      eq2 = - t(n)
      b1  =   a1(n)
      b2  =   a1(n)
      b3  =   a3(n)
      b4  =   a4(n)
      det =   b1 * b4 - b2 * b3
      c1  =   b4 / det
      c2  = - b2 / det
      c3  = - b3 / det
      c4  =   b1 / det
      da  =   c1 * eq1 + c2 * eq2
      db  =   c3 * eq1 + c4 * eq2
      s0  =   s0 + da
      r0  =   r0 + db
      j = j + 1
      if (abs(u(n)-1._rkind)<eps.and.abs(t(n))<eps) exit
    enddo
!   
!   ********************************* b.l. thickness
!   
    dd = .99_rkind
    do i = 0, nm
      if (u(i).ge.dd) exit
    end do
    delta1 = eta(i-1)
    kk = i - 1
!   
!   ******************************* displacement thickness
!   
    delta2 = delta1 -2*h(kk)
!   
!   ********************************* Crocco's temperature profile
!   
    do i = 0,n
      tcr(i) = tbl(0)/Te  + (1._rkind - tbl(0)/Te) * u(i) + (gam - 1._rkind) / 2._rkind *Mach**2. * u(i) * (1._rkind - u(i))
    end do
!   
!   ********************************************* printing
!   
    do i=0,n
      g1  = .5_rkind*u(i) / tbl(i) * Te
      rho = Te/tbl(i)
      v(i) =  eta(i)*g1 - h(i)
    end do
!   
    return
  end subroutine compressible_blasius
! 
  subroutine meanvelocity_bl(retau,rm,trat,s2tinf,vtexp,visc_exp,gam,rfac,ny,y,u,rho,t,visc,redelta,cf,th)
!   
    logical, intent(in) :: visc_exp
    integer, intent(in) :: ny
    real(rkind), intent(in)  :: retau,rm,trat,s2tinf,gam,rfac,vtexp
    real(rkind), intent(out) :: redelta,cf,th
    real(rkind), dimension(ny), intent(in)  :: y
    real(rkind), dimension(ny), intent(out) :: u,rho,t,visc
!   
    integer :: j,m,itransf,jp99,jj99,jj
    real(rkind) :: gm1h, tr, tw, alf, retau_inc, retau_target, retau_correction, retau_musker
    real(rkind) :: retau_inc_target, retau_inc_correction, up99, d99plus_inc
    real(rkind) :: fuu, du, rhow_over_rho, rhow_over_rhojm, rmu_over_rmuw, rmu_over_rmuwjm
    real(rkind) :: ucij, ucijm, res, dy, ycij, ycijm, u99, d99, d99plus, yy, ui
    real(rkind) :: ue, uuu, uuum, thi
    real(rkind) :: tol_iter_loc
!   
    real(rkind), dimension(ny) :: yplus       ! Incompressible (transformed) y
    real(rkind), dimension(ny) :: uplus       ! Incompressible (transformed )velocity profile
    real(rkind), dimension(ny) :: ycplus      ! Compressible y
    real(rkind), dimension(ny) :: ucplus      ! Compressible velocity profile
    real(rkind), dimension(ny) :: yplusold
    real(rkind), dimension(ny) :: ucplusold
    real(rkind), dimension(ny) :: density     ! Density profile
    real(rkind), dimension(ny) :: viscosity   ! Viscosity profile
    real(rkind), dimension(ny) :: temperature ! Temperature profile
    real(rkind), dimension(ny) :: uu,uci,uc
!   
    itransf = 0 ! Volpiani transformation
    gm1h    = 0.5_rkind*(gam-1._rkind)
    tr      = 1._rkind+gm1h*rfac*rm**2
    tw      = trat*tr
    alf     = 0.8259_rkind
    tol_iter_loc = 0.000001_rkind
!   
    m = 4
    retau_inc        = retau
    retau_target     = retau
    retau_correction = 1._rkind
    do ! Outer loop to find retau_inc
      retau_inc            = retau_inc*retau_correction
!     
      retau_musker         = retau_inc
      retau_inc_target     = retau_inc
      retau_inc_correction = 1._rkind
      do ! Local loop to find retau_musker
        retau_musker = retau_musker*retau_inc_correction
        uplus = 0._rkind
        do j=2,ny
          yplus(j) = y(j)*retau_musker
          uplus(j) = velmusker(yplus(j),retau_musker) ! Incompressible velocity profile (Musker, AIAA J 1979)
        enddo
        up99 = 0.99_rkind*uplus(ny)
        call locateval(uplus,ny,up99,jp99)
        jp99 = min(max(jp99-(m-1)/2,1),ny+1-m)
        call pol_int(uplus(jp99),yplus(jp99),m,up99,d99plus_inc)
        retau_inc_correction = retau_inc_target/d99plus_inc
        if (abs(retau_inc_correction-1._rkind)<tol_iter_loc) exit
      enddo
!     
      ycplus   = y*retau_musker
      yplus    = ycplus  ! Assume as initial guess incompressible (transformed) yplus equal to ycplus
      yplusold = yplus
!     
      do ! outer loop to find yplus
!       
        uplus = 0._rkind
        do j=2,ny
          uplus(j) = velmusker(yplus(j),retau_musker) ! Incompressible velocity profile (Musker, AIAA J 1979)
        enddo
!       
        ucplus = uplus ! Assume as initial guess compressible velocity profile equal to incompressible
!       
        do
!         
          ucplusold = ucplus
          do j=1,ny
            uu(j)          = ucplus(j)/ucplus(ny)
!           fuu            = alf*uu(j)+(1._rkind-alf)*uu(j)**2 ! Duan and Martin (see Zhang JFM 2014)
            fuu            = uu(j)                             ! Walz
            temperature(j) = tw+(tr-tw)*fuu+(1._rkind-tr)*uu(j)**2
            density(j)     = 1._rkind/temperature(j)
            if (visc_exp) then
              viscosity(j) = temperature(j)**vtexp ! Power-law
            else
              viscosity(j)   = temperature(j)**(1.5_rkind)*(1._rkind+s2tinf)/(temperature(j)+s2tinf)
            endif
          enddo
!         
!         Inverse Van Driest or Volpiani transformation for u
!         
          do j=2,ny
            du              = uplus(j)-uplus(j-1)
            rhow_over_rho   = density(1)/density(j)
            rhow_over_rhojm = density(1)/density(j-1)
            rmu_over_rmuw   = viscosity(j)/viscosity(1)
            rmu_over_rmuwjm = viscosity(j-1)/viscosity(1)
            if (itransf==1) then
!             VD
              ucij  = sqrt(rhow_over_rho)
              ucijm = sqrt(rhow_over_rhojm)
            else
!             Volpiani
              ucij  = sqrt(rhow_over_rho)*sqrt(rmu_over_rmuw)
              ucijm = sqrt(rhow_over_rhojm)*sqrt(rmu_over_rmuwjm)
            endif
!           
            ucplus(j) = ucplus(j-1)+0.5_rkind*(ucij+ucijm)*du
!           
          enddo
!         
          res = 0._rkind
          do j=1,ny
            res = res+abs(ucplus(j)-ucplusold(j))
          enddo
          res = res/ny
          if (res<tol_iter_loc) exit
!         
        enddo ! End of iterative loop to find ucplus
!       
!       At this stage we have ucplus and density, temperature, viscosity distributions
!       
!       Compute incompressible yplus with transformation
!       
        yplus = 0._rkind
        do j=2,ny
          dy = ycplus(j)-ycplus(j-1)
!         
          rhow_over_rho   = density(1)/density(j)
          rhow_over_rhojm = density(1)/density(j-1)
          rmu_over_rmuw   = viscosity(j)/viscosity(1)
          rmu_over_rmuwjm = viscosity(j-1)/viscosity(1)
!         
!         ! VD
          if (itransf==1) then
            ycij  = 1._rkind
            ycijm = 1._rkind
          else
!           ! Volpiani
            ycij  = sqrt(rhow_over_rho)*sqrt(rmu_over_rmuw**3)
            ycij  = 1._rkind/ycij
            ycijm = sqrt(rhow_over_rhojm)*sqrt(rmu_over_rmuwjm**3)
            ycijm  = 1._rkind/ycijm
          endif
!         
          yplus(j) = yplus(j-1)+0.5_rkind*(ycij+ycijm)*dy
        enddo
!       
        res  = 0._rkind
        do j=1,ny
          res  = abs(yplus(j)-yplusold(j))
        enddo
        res  = res/ny
        if (res<tol_iter_loc) exit
        yplusold = yplus
      enddo
!     
      do j=1,ny
        uc(j) = ucplus(j)/ucplus(ny)
      enddo
      u99 = 0.99_rkind
      call locateval(uc,ny,u99,jj99)
      jj99 = min(max(jj99-(m-1)/2,1),ny+1-m)
      call pol_int(uc(jj99),y(jj99),m,u99,d99)
      d99plus = d99*retau_musker
      retau_correction = retau_target/d99plus
      if (abs(retau_correction-1._rkind) < tol_iter_loc) exit
!     
    enddo ! Outer loop
!   
    do j=1,ny
      uc(j) = ucplus(j)/ucplus(ny)
    enddo
    u99 = 0.99_rkind
    call locateval(uc,ny,u99,jj99)
    jj99 = min(max(jj99-(m-1)/2,1),ny+1-m)
    call pol_int(uc(jj99),y(jj99),m,u99,d99)
!   
    do j=1,ny
      yy = y(j)*d99
      call locateval(y,ny,yy,jj)
      jj = min(max(jj-(m-1)/2,1),ny+1-m)
      call pol_int(y(jj),uc(jj),m,yy,ui)
      uci(j) = ui
    enddo
    uc = uci
    do j=1,ny
      uu(j) = uc(j)/uc(ny)
!     fuu   = alf*uu(j)+(1._rkind-alf)*uu(j)**2 ! Duan and Martin (see Zhang JFM 2014)
      fuu   = uu(j)                             ! Walz
      temperature(j) = tw+(tr-tw)*fuu+(1.-tr)*uu(j)**2
      density(j)     = 1._rkind/temperature(j)
      if (visc_exp) then
        viscosity(j) = temperature(j)**vtexp ! Power-law
      else
        viscosity(j) = temperature(j)**(1.5_rkind)*(1._rkind+s2tinf)/(temperature(j)+s2tinf)
      endif
    enddo
!   
    retau_musker = retau_musker*d99
    redelta      = retau_musker*ucplus(ny)/density(1)*viscosity(1)
    cf           = 2._rkind*density(1)/ucplus(ny)**2
!   dstar        = 0._rkind
    th           = 0._rkind
    ue           = uc(ny)
!   
    do j=2,ny
      uuu     = uc(j)
      uuum    = uc(j-1)
      dy     = y(j)-y(j-1)
!     dstari = 0.5*((1.-density(j)*uuu)+(1.-density(j-1)*uuum))
!     dstar  = dstar+dstari*dy
      thi    = 0.5_rkind*(density(j)*uuu*(1.-uuu)+density(j-1)*uuum*(1.-uuum))
      th     = th+thi*dy
    enddo
!   
    do j=1,ny
      u(j)    = uc(j)
      rho(j)  = density(j)
      t(j)    = temperature(j)
      visc(j) = viscosity(j)
    enddo
!   
    return
  end subroutine meanvelocity_bl
! 
  subroutine init_wind_tunnel(self)
    class(equation_singleideal_object), intent(inout) :: self
!
    integer :: i,j,k,l
!   only ghost in x to be similar to bl needs
    allocate(self%wmean(1-self%grid%ng:self%field%nx+self%grid%ng+1, 1:self%field%ny, 4))
!
    associate(nx => self%field%nx, ny => self%field%ny, nz => self%field%nz, ng => self%grid%ng, nv => self%nv, &
      wmean => self%wmean, w => self%field%w, winf => self%winf)
!
      wmean = 0._rkind
      do j=1,ny
        do i=1-ng,nx+ng+1
          wmean(i,j,1) = winf(1)
          wmean(i,j,2) = winf(2)
          wmean(i,j,3) = winf(3)
          wmean(i,j,4) = winf(4)
        enddo
      enddo
!
      do k=1,nz
        do j=1,ny
          do i=1,nx
            do l=1,nv
              w(l,i,j,k) = winf(l)
            enddo
          enddo
        enddo
      enddo
!
    endassociate
  endsubroutine init_wind_tunnel
!
  subroutine init_tgv(self)
    class(equation_singleideal_object), intent(inout) :: self
!
    integer :: i,j,k
    real(rkind) :: rho, uu, vv, ww
    real(rkind) :: rhouu, rhovv, rhoww
    real(rkind) :: tt, ee, pp
!   only ghost in x to be similar to bl needs
    allocate(self%wmean(1-self%grid%ng:self%field%nx+self%grid%ng+1, 1:self%field%ny, 4))
!
    associate(nx => self%field%nx, ny => self%field%ny, nz => self%field%nz, ng => self%grid%ng, nv => self%nv, &
      wmean => self%wmean, w => self%field%w, winf => self%winf, w_aux => self%w_aux, &
      rho0 => self%rho0, u0 => self%u0, p0 => self%p0, gam => self%gam, y => self%field%y, z => self%field%z, &
      x => self%field%x, indx_cp_l => self%indx_cp_l, indx_cp_r => self%indx_cp_r, cv_coeff => self%cv_coeff, &
      calorically_perfect => self%calorically_perfect, rgas0 => self%rgas0, t0 => self%t0)
!
      wmean = 0._rkind
      do j=1,ny
        do i=1-ng,nx+ng+1
          wmean(i,j,1) = winf(1)
          wmean(i,j,2) = winf(2)
          wmean(i,j,3) = winf(3)
          wmean(i,j,4) = winf(4)
        enddo
      enddo
!
      do k=1,nz
        do j=1,ny
          do i=1,nx
            rho   = rho0
            uu    =  u0*sin(x(i))*cos(y(j))*cos(z(k))
            vv    = -u0*cos(x(i))*sin(y(j))*cos(z(k))
            ww    =  0._rkind
            pp    = 1._rkind+u0**2/16._rkind*(cos(2._rkind*x(i))+cos(2._rkind*y(j)))*(2._rkind+cos(2.*z(k)))
            rhouu = rho*uu
            rhovv = rho*vv
            rhoww = rho*ww
            w(1,i,j,k) = rho
            w(2,i,j,k) = rhouu
            w(3,i,j,k) = rhovv
            w(4,i,j,k) = rhoww
            tt         = pp/rgas0/rho
            w_aux(i,j,k,6) = tt
            ee = get_e_from_temperature(tt, t0, indx_cp_l, indx_cp_r, cv_coeff, calorically_perfect)
            w(5,i,j,k) = rho*ee + 0.5_rkind*(rhouu**2+rhovv**2+rhoww**2)/rho
          enddo
        enddo
      enddo
!
    endassociate
  endsubroutine init_tgv
!
  subroutine init_channel(self)
    class(equation_singleideal_object), intent(inout) :: self
    integer :: i, j, k
    real(rkind), dimension(3) :: rr
    real(rkind) :: u0_02, u0_05, rho, uu, vv, ww
    real(rkind) :: ufluc, vfluc, wfluc
    real(rkind) :: rhouu, rhovv, rhoww
    real(rkind) :: tt, ee
    integer :: nroll
!
!   only ghost in x to be similar to bl needs
    allocate(self%wmean(1-self%grid%ng:self%field%nx+self%grid%ng+1, 1:self%field%ny, 4))
!
    associate(nx => self%field%nx, ny => self%field%ny, nz => self%field%nz, ng => self%grid%ng, &
      y => self%field%y, z => self%field%z, &
      rho0 => self%rho0, u0 => self%u0, p0 => self%p0, gm => self%gm, t0 => self%t0, &
      wmean => self%wmean, w => self%field%w, w_aux => self%w_aux, l0 => self%l0, &
      indx_cp_l => self%indx_cp_l, indx_cp_r => self%indx_cp_r, cv_coeff => self%cv_coeff, &
      volchan => self%volchan, nxmax => self%grid%nxmax, nzmax => self%grid%nzmax, yn => self%grid%yn, &
      lz => self%grid%domain_size(3), calorically_perfect => self%calorically_perfect, rgas0 => self%rgas0)
!
      volchan = yn(ny+1)-yn(1)
      volchan = volchan*nxmax
      volchan = volchan*nzmax
!
      wmean = 0._rkind
      do j=1,ny
        do i=1-ng,nx+ng+1
          wmean(i,j,1) = rho0
          wmean(i,j,2) = 1.5_rkind*rho0*u0*(1._rkind-(y(j)/l0)**2)
          wmean(i,j,3) = 0._rkind
          wmean(i,j,4) = 0._rkind
        enddo
      enddo
!
      u0_02 = 0.02_rkind*u0
      u0_05 = 0.05_rkind*u0
      nroll = int(lz/yn(ny+1))
!
      do k=1,nz
        do j=1,ny
          do i=1,nx
!
            rho = wmean(i,j,1)
            uu  = wmean(i,j,2)/rho
            vv  = wmean(i,j,3)/rho
            ww  = wmean(i,j,4)/rho
!
            call get_crandom_f(rr)
            rr = rr-0.5_rkind
!
            ufluc   = u0_02*rr(1)
            vfluc   = u0_02*rr(2)
            wfluc   = u0_02*rr(3)
            vfluc   = vfluc+u0_05*sin(0.5_rkind*pi*y(j)/l0)*cos(2*pi*z(k)/lz*nroll)
            wfluc   = wfluc+u0_05*sin(0.5_rkind*pi*y(j)/l0)*sin(2*pi*z(k)/lz*nroll)
            uu      = uu+ufluc
            vv      = vv+vfluc
            ww      = ww+wfluc
            rhouu   = rho*uu
            rhovv   = rho*vv
            rhoww   = rho*ww
            w(1,i,j,k) = rho
            w(2,i,j,k) = rhouu
            w(3,i,j,k) = rhovv
            w(4,i,j,k) = rhoww
            tt        = p0/rgas0/rho
            w_aux(i,j,k,6) = tt
            ee = get_e_from_temperature(tt, t0, indx_cp_l, indx_cp_r, cv_coeff, calorically_perfect)
            w(5,i,j,k) = rho*ee + 0.5_rkind*(rhouu**2+rhovv**2+rhoww**2)/rho
          enddo
        enddo
      enddo
!
    endassociate
  endsubroutine init_channel
!
  subroutine bc_preproc(self)
    class(equation_singleideal_object), intent(inout) :: self
    integer :: ilat, offset
!
    self%force_zero_flux = [0,0,0,0,0,0]
    self%eul_imin = 1
    self%eul_imax = self%field%nx
    self%eul_jmin = 1
    self%eul_jmax = self%field%ny
    self%eul_kmin = 1
    self%eul_kmax = self%field%nz
    do ilat=1,6 ! loop on all sides of the boundary (3D -> 6)
      select case(self%bctags(ilat))
      case(0)
      case(1)
      case(2)
      case(4)
      case(5)
      case(6)
        if (self%channel_case) self%force_zero_flux(ilat) = 1
      case(7)
      case(8)
      case(9)
      endselect
      if(self%bctags_nr(ilat) /= 0) then
        offset = 1
      else
        offset = 0
      endif
      select case(ilat)
      case(1)
        self%eul_imin = self%eul_imin + offset
      case(2)
        self%eul_imax = self%eul_imax - offset
      case(3)
        self%eul_jmin = self%eul_jmin + offset
      case(4)
        self%eul_jmax = self%eul_jmax - offset
      case(5)
        self%eul_kmin = self%eul_kmin + offset
      case(6)
        self%eul_kmax = self%eul_kmax - offset
      endselect
    enddo
  endsubroutine bc_preproc
!
  function get_gamloc(tt,t0,indx_cp_l,indx_cp_r,cp_coeff,calorically_perfect,rgas0)
    real(rkind) :: get_gamloc
    integer :: indx_cp_l,indx_cp_r,calorically_perfect
    real(rkind) :: tt, rgas0, t0
    real(rkind), dimension(indx_cp_l:indx_cp_r+1) :: cp_coeff
    real(rkind) :: cploc, gamloc
    integer :: l
!
    if (calorically_perfect==1) then
      cploc = cp_coeff(0)
    else
      cploc = 0._rkind
      do l=indx_cp_l,indx_cp_r
        cploc = cploc+cp_coeff(l)*(tt/t0)**l
      enddo
    endif
    gamloc = cploc/(cploc-rgas0)
    get_gamloc = gamloc
!
  endfunction get_gamloc
!
  subroutine runge_kutta_initialize(self)
!   < Initialize Runge-Kutta data.
    class(equation_singleideal_object), intent(inout) :: self !< The equation.
!
    select case(self%rk_type)
    case(RK_WRAY)
      self%nrk = 3
      allocate(self%rhork(self%nrk),self%gamrk(self%nrk),self%alprk(self%nrk))
      self%rhork(:) = [0._rkind, -17._rkind/60._rkind , -5._rkind /12._rkind]
      self%gamrk(:) = [8._rkind  /15._rkind, 5._rkind  /12._rkind, 3._rkind  /4._rkind]
      self%alprk(:) = self%rhork(:) + self%gamrk(:)
    case(RK_JAMESON)
      self%nrk = 4
      allocate(self%alprk(self%nrk))
      self%alprk(:) = [1._rkind/4._rkind, 1._rkind/3._rkind, 1._rkind/2._rkind, 1._rkind]
    case(RK_SHU)
      self%nrk = 3
      allocate(self%ark(self%nrk),self%brk(self%nrk),self%crk(self%nrk))
      self%ark(:) = [1._rkind, 0.75_rkind, 1._rkind /3._rkind]
      self%brk(:) = [0._rkind, 0.25_rkind, 2._rkind /3._rkind]
      self%crk(:) = [1._rkind, 0.25_rkind, 2._rkind /3._rkind]
    endselect
!
  endsubroutine runge_kutta_initialize
!
  subroutine read_input(self, filename)
!   < Initialize the equation.
    class(equation_singleideal_object), intent(inout) :: self              !< The equation.
    character(*)                    , intent(in)      :: filename          !< Input file name.
!
!   Test and possibly use .ini file as input format
!   https://github.com/pkgpl/cfgio
    self%cfg=parse_cfg(filename)
!
  endsubroutine read_input
!
! 
! 
! 
!
! 
!
!
! 
!
! 
! 
! 
! 
endmodule streams_equation_singleideal_object
!
!
