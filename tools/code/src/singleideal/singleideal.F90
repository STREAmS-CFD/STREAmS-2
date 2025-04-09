module streams_equation_singleideal_object

  use streams_field_object
  use streams_grid_object
  use streams_parameters
  use cfgio_mod, only: cfg_t, parse_cfg
  use, intrinsic :: iso_fortran_env, only : error_unit
  use crandom_f_mod
  use MPI
  use ISO_C_BINDING


  use catalyst_api
  use catalyst_conduit

  implicit none
  private
  public :: equation_singleideal_object
  public :: VISC_POWER, VISC_SUTHERLAND, VISC_NO
  public :: RK_WRAY, RK_JAMESON, RK_SHU
  public :: slicexy_aux, slicexz_aux, sliceyz_aux
  public :: IBM_MAX_PARBC

  integer(ikind), parameter :: VISC_NO = 0_ikind
  integer(ikind), parameter :: VISC_POWER = 1_ikind
  integer(ikind), parameter :: VISC_SUTHERLAND = 2_ikind

  integer(ikind), parameter :: RK_WRAY = 1_ikind
  integer(ikind), parameter :: RK_JAMESON = 2_ikind
  integer(ikind), parameter :: RK_SHU = 3_ikind

  integer(ikind), parameter :: IBM_MAX_PARBC = 3_ikind

  real(rkind), dimension(:,:,:,:), allocatable :: slicexy_aux, slicexz_aux, sliceyz_aux

  type :: equation_singleideal_object
    type(grid_object) :: grid
    type(field_object) :: field
    type(cfg_t) :: cfg
    type(cfg_t) :: flow_params_cfg
    type(cfg_t) :: field_info_cfg

    integer :: gpu_bind
    integer(ikind) :: nx, ny, nz
    integer(ikind) :: ng
    integer(ikind) :: mpi_err
    integer(ikind) :: myrank
    integer(ikind) :: nprocs
    logical :: masterproc

    integer(ikind) :: nv
    integer(ikind) :: nv_aux
    integer(ikind) :: nv_stat
    integer(ikind) :: nv_stat_3d
    integer(ikind) :: enable_stat_time
    integer(ikind) :: enable_stat_3d

    real(rkind) :: dt
    real(rkind) :: cfl
    real(rkind) :: time0, time
    integer :: num_iter
    integer :: icyc0, icyc

    integer(ikind) :: visc_order, conservative_viscous
    integer(ikind) :: ep_order, nkeep
    integer(ikind) :: weno_scheme, weno_version, flux_splitting
    real(rkind) :: sensor_threshold
    real(rkind) :: xshock_imp, shock_angle, tanhfacs
    integer :: rand_type
    integer :: eul_imin, eul_imax, eul_jmin, eul_jmax, eul_kmin, eul_kmax
    integer, dimension(6) :: force_zero_flux
    integer(ikind), dimension(:), allocatable :: supported_orders

    real(rkind), dimension(1:4,4) :: coeff_deriv1
    real(rkind), dimension(0:4,4) :: coeff_deriv2
    real(rkind), dimension(0:4) :: coeff_filter

    integer(ikind) :: rk_type
    integer(ikind) :: nrk
    real(rkind), allocatable, dimension(:) :: rhork,gamrk,alprk
    real(rkind), allocatable, dimension(:) :: ark,brk,crk

    integer :: iter_dt_recompute, print_control

    real(rkind) :: dtsave, dtsave_restart, dtstat, dtslice, dtprobe, dtslice_p3d
    real(rkind) :: time_from_last_rst, time_from_last_write, time_from_last_stat
    real(rkind) :: time_from_last_slice, time_from_last_probe, time_from_last_slice_p3d
    real(rkind) :: time_from_last_tspec
    integer :: restart_type
    integer :: io_type_r
    integer :: io_type_w
    integer :: restart_from_flow
    integer :: istore
    integer :: itav, itslice_p3d
    integer :: enable_plot3d, enable_vtk
    integer :: enable_slice_plot3d, slice_plot3d_j1, slice_plot3d_j2
    real(rkind) :: residual_rhou,vmax

    integer :: flow_init
    real(rkind) :: Mach, Reynolds, Prandtl, theta_wall
    real(rkind) :: aoa
    real(rkind) :: Reynolds_friction
    real(rkind) :: rgas0
    real(rkind) :: rho0, u0, t0, s0, p0, mu0, ptot0, ttot0
    real(rkind) :: l0, l0_ducros
    real(rkind) :: gam, gm, gm1, rfac
    integer :: rfac_type
    real(rkind) :: T_wall, T_bulk_target, T_recovery
    real(rkind) :: powerlaw_vtexp, T_ref_dim, sutherland_S
    real(rkind) :: rhomin, rhomax, tmin, tmax, pmin, pmax
    integer :: visc_model
    real(rkind) :: dpdx,rhobulk,ubulk,tbulk,volchan,dpth
    logical :: channel_case, bl_case, bl_laminar
    real(rkind) :: x_recyc
    logical :: recyc
    integer :: i_recyc, ib_recyc
    real(rkind), dimension(:), allocatable :: winf, winf_past_shock
    real(rkind), dimension(:,:), allocatable :: wfar, wfar_glob

    real(rkind), dimension(:,:,:,:), allocatable :: w_aux
    real(rkind), dimension(:,:,:), allocatable :: wmean
    integer, dimension(:,:,:), allocatable :: fluid_mask
    integer, dimension(:,:,:,:), allocatable :: ep_ord_change
    integer :: correct_bound_ord

    real(rkind), dimension(:,:,:), allocatable :: w_stat, w_stat_z
    real(rkind), dimension(:,:,:,:), allocatable :: w_stat_3d
    logical :: save_stat_z

    integer, dimension(6) :: bctags, bctags_nr

    real(rkind), dimension(:), allocatable :: deltavec, deltavvec, cfvec
    real(rkind) :: betarecyc, glund1
    real(rkind), dimension(:), allocatable :: yplus_inflow
    real(rkind), dimension(:), allocatable :: eta_inflow
    real(rkind), dimension(:), allocatable :: yplus_recyc
    real(rkind), dimension(:), allocatable :: eta_recyc,eta_recyc_blend
    integer, dimension(:), allocatable :: map_j_inn,map_j_out,map_j_out_blend
    integer :: jbl_inflow
    real(rkind), dimension(:), allocatable :: weta_inflow
    real(rkind), dimension(:,:,:), allocatable :: inflow_random_plane

    integer :: calorically_perfect, indx_cp_l, indx_cp_r
    real(rkind), allocatable, dimension(:) :: cp_coeff, cv_coeff
    integer(ikind), dimension(:), allocatable :: igslice, jgslice, kgslice
    integer(ikind), dimension(:), allocatable :: islice, jslice, kslice
    integer(ikind), dimension(:), allocatable :: list_aux_slice
    integer(ikind) :: num_probe, num_probe_tot, num_aux_slice
    real(rkind), dimension(:,:), allocatable :: probe_coord, w_aux_probe
    real(rkind), dimension(:,:,:,:), allocatable :: probe_coeff
    integer, dimension(:,:), allocatable :: ijk_probe
    integer, dimension(:), allocatable :: moving_probe, id_probe

    integer(ikind), dimension(:), allocatable :: igslice_p3d, jgslice_p3d, kgslice_p3d
    integer(ikind), dimension(:), allocatable :: islice_p3d, jslice_p3d, kslice_p3d
    real(rkind), allocatable, dimension(:,:,:) :: wallprop

    integer :: enable_tspec
    real(rkind) :: dt_tspec, wss_tspec
    integer :: ndft_tspec, jg_tspec, igstart_tspec, igend_tspec
    integer :: it_win1_tspec, it_win2_tspec, it_nwin_tspec
    real(rkind), dimension(:,:,:,:,:), allocatable :: w_tspec
    real(rkind), dimension(:,:), allocatable :: w_psd_tspec
    real(rkind), dimension(:,:), allocatable :: w_psd_avg_tspec

    integer :: enable_sponge, j_sponge
    real(rkind), dimension(:), allocatable :: f_sponge
    real(rkind) :: k_sponge, theta_threshold

    integer :: enable_forces_runtime, count_weno_control
    real(rkind) :: count_weno_x, count_weno_y, count_weno_z
    real(rkind) :: xtr1, xtr2, x0tr, y0tr, x0ts, y0ts, xtrip, xtw1, xtw2, a_tr, a_tw, v_bs, kx_tw, om_tw, thic, teshk
    integer :: itr, its, is, is_old, itr1, itr2, its1, its2, jtr, jts
    real(rkind) :: lamz , lamz1 , lams , lams1 , phiz , phiz1 , phis , phis1, del0_tr
    real(rkind) :: lamz_old , lamz1_old , lams_old , lams1_old , phiz_old , phiz1_old , phis_old , phis1_old

    integer, dimension(:), allocatable :: lmax_tag
    integer, dimension(:), allocatable :: vis_tag
    real(rkind) :: lift_af, drag_af, pres_af, fric_af
    real(rkind) :: rho_lim, tem_lim, rho_lim_rescale, tem_lim_rescale
    integer :: enable_limiter, ifilter, jfilter, jweno
    integer :: ortho

    integer :: debug_memory = 0
    integer :: mode_async = 0

    logical :: time_is_freezed=.false.
    !insitu_var_start
    integer :: enable_insitu = 0
    integer :: i_freeze = 0
    integer :: i_insitu
    real(rkind) :: dt_insitu
    real(rkind) :: time_from_last_insitu
    real(rkind) :: time_insitu
    character(len=128) :: vtkpipeline
    integer, allocatable, dimension(:) :: aux_list
    integer, allocatable, dimension(:) :: add_list
    integer, allocatable, dimension(:) :: freeze_intervals
    character(len=64), allocatable, dimension(:) :: aux_list_name
    character(len=64), allocatable, dimension(:) :: add_list_name
    logical :: fcoproc
    logical :: enable_freeze_intervals
    real(rkind) :: perc_ny_cut
    integer :: ny_cut
    integer :: ngm, ngp
    integer :: nxsl_ins, nxel_ins, nysl_ins, nyel_ins, nzsl_ins, nzel_ins
    integer :: nxs_ins, nxe_ins, nys_ins, nye_ins, nzs_ins, nze_ins
    integer :: nxstartg, nxendg, nystartg, nyendg, nzstartg, nzendg
    integer :: npsi, npsi_pv, n_aux_list, n_add_list
    real(rkind), allocatable, dimension(:,:,:,:) :: psi, psi_pv
    real(rkind), allocatable, dimension(:) :: xyzc
    integer :: insitu_flag
    character(64) :: insitu_platform
    real(rkind), allocatable, dimension(:) :: points_x, points_y, points_z
    integer :: n_points_x, n_points_y, n_points_z
    integer(int64) :: n_points
    !insitu_var_end

  contains
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
    procedure, pass(self) :: init_channel_c2
    procedure, pass(self) :: init_airfoil
    procedure, pass(self) :: init_wind_tunnel
    procedure, pass(self) :: init_bl_lam
    procedure, pass(self) :: init_bl_old
    procedure, pass(self) :: init_bl
    procedure, pass(self) :: init_cv
    procedure, pass(self) :: init_sponge
    procedure, pass(self) :: alloc
    procedure, pass(self) :: bc_preproc
    procedure, pass(self) :: compute_stats
    procedure, pass(self) :: compute_stats_c2
    procedure, pass(self) :: compute_airfoil_forces
    procedure, pass(self) :: write_stats
    procedure, pass(self) :: write_stats_z
    procedure, pass(self) :: write_stats_serial
    procedure, pass(self) :: write_stats_z_serial
    procedure, pass(self) :: read_stats
    procedure, pass(self) :: read_stats_serial
    procedure, pass(self) :: compute_stats_3d
    procedure, pass(self) :: write_stats_3d
    procedure, pass(self) :: read_stats_3d
    procedure, pass(self) :: write_stats_3d_serial
    procedure, pass(self) :: read_stats_3d_serial
    procedure, pass(self) :: recyc_prepare
    procedure, pass(self) :: add_synthetic_perturbations
    procedure, pass(self) :: slice_prepare
    procedure, pass(self) :: probe_prepare
    procedure, pass(self) :: probe_compute_coeff
    procedure, pass(self) :: correct_bc_order
    procedure, pass(self) :: correct_airfoil_order
    procedure, pass(self) :: write_slice_p3d
    procedure, pass(self) :: fromaux_to_w
    procedure, pass(self) :: write_tspec
    procedure, pass(self) :: write_psd_tspec
    procedure, pass(self) :: write_psd_avg_tspec
    procedure, pass(self) :: read_tspec
    procedure, pass(self) :: time_is_freezed_fun

    !insitu_proc_start
    procedure, pass(self) :: insitu_initialize
    procedure, pass(self) :: insitu_finalize
    procedure, pass(self) :: insitu_allocate
    procedure, pass(self) :: insitu_define_limits
    procedure, pass(self) :: insitu_do_catalyst_initialization
    procedure, pass(self) :: insitu_do_catalyst_finalization
    !insitu_proc_end
  endtype equation_singleideal_object

contains

  subroutine initialize(self, filename)
    class(equation_singleideal_object), intent(inout) :: self
    character(*) , intent(in) :: filename
    logical, dimension(3) :: periodic
    integer, dimension(3) :: mpi_splits
    integer :: mpi_split_x ,mpi_split_y ,mpi_split_z
    integer :: nxmax, nymax, nzmax, ng, l
    integer :: grid_type, metrics_order, ystag, grid_dim
    real(rkind) :: domain_size_x, domain_size_y, domain_size_z, rtemp, wind
    real(rkind), dimension(:), allocatable :: grid_vars
    logical :: rebuild_ghost, ystaggering
    integer :: bctag_temp, order, itemp
    real(rkind) :: d1_temp(1:4), d2_temp(0:4)
    integer :: i, grid2d_par, save_metrics_2d, istat


    call get_mpi_basic_info(self%nprocs, self%myrank, self%masterproc, self%mpi_err)

    call self%read_input(filename)

    if (self%cfg%has_key("output","debug_memory")) call self%cfg%get("output","debug_memory",self%debug_memory)
    if(.not. any(self%debug_memory == [0,1])) call fail_input_any("debug_memory must be 0 or 1")

    if (self%cfg%has_key("numerics","mode_async")) call self%cfg%get("numerics","mode_async",self%mode_async)
    if(.not. any(self%mode_async == [0,1])) call fail_input_any("mode_async must be 0 or 1")

    self%gpu_bind = 0
    if (self%cfg%has_key("grid","gpu_bind")) call self%cfg%get("grid","gpu_bind",self%gpu_bind)
    if(.not. any(self%gpu_bind == [0,1])) call fail_input_any("gpu_bind must be 0 or 1")

    call self%cfg%get("numerics","rand_type",self%rand_type)
    if(self%rand_type == 0) then
      call init_crandom_f(0,reproducible=.true.)
      if (self%masterproc) write(*,*) 'Random numbers disabled'
    elseif(self%rand_type < 0) then
      call init_crandom_f(self%myrank+1,reproducible=.false.)
      if (self%masterproc) write(*,*) 'Random numbers NOT reproducible'
    else
      call init_crandom_f(self%myrank+1,reproducible=.true.)
      if (self%masterproc) write(*,*) 'Random numbers reproducible'
    endif

    self%l0 = 1._rkind
    self%t0 = 1._rkind
    self%rho0 = 1._rkind
    self%rgas0 = 1._rkind
    if (self%cfg%has_key("ref_prop","l0")) call self%cfg%get("ref_prop","l0",self%l0)
    if (self%cfg%has_key("ref_prop","t0")) call self%cfg%get("ref_prop","t0",self%t0)
    if (self%cfg%has_key("ref_prop","rho0")) call self%cfg%get("ref_prop","rho0",self%rho0)
    if (self%cfg%has_key("ref_prop","rgas0")) call self%cfg%get("ref_prop","rgas0",self%rgas0)

    grid_dim = 1
    if (self%cfg%has_key("grid","grid_dim")) call self%cfg%get("grid","grid_dim",grid_dim)
    if(.not. any(grid_dim == [1,2])) call fail_input_any("grid_dim must be 1 or 2")
    if (self%masterproc) print*,'Grid Cartesian (1) / Curvilinear (2) : ',grid_dim

    grid2d_par = 2
    if (grid_dim == 2) then
      if (self%cfg%has_key("grid","grid2d_par")) then
        call self%cfg%get("grid","grid2d_par",grid2d_par)
      endif
    endif
    if(.not. any(grid2d_par == [0,1,2,3])) call fail_input_any("grid2d_par must be 0, 1, 2, or 3")

    save_metrics_2d = 1
    if (grid_dim == 2) then
      if (self%cfg%has_key("grid","save_metrics_2d")) then
        call self%cfg%get("grid","save_metrics_2d",save_metrics_2d)
      endif
    endif
    if(.not. any(save_metrics_2d == [0,1])) call fail_input_any("save_metrics_2d must be 0 or 1")

    call self%cfg%get("flow","flow_init",self%flow_init)
    if (self%flow_init<-1 .or.self%flow_init>6) call fail_input_any("flow_init not implemented")
    self%channel_case = .false.
    if (self%flow_init == 0 .or. self%flow_init == 4) self%channel_case = .true.
    self%bl_case = .false.
    if (self%flow_init == 1 .or. self%flow_init == 3) then
      self%bl_case = .true.
      self%bl_laminar = .false.
    endif

    self%aoa = 0._rkind
    if (self%cfg%has_key("flow","aoa")) then
      call self%cfg%get("flow","aoa",self%aoa)
    endif
    self%l0_ducros = self%l0
    if (self%flow_init == 5) then
      if (self%cfg%has_key("curvi","l0_ducros")) then
        call self%cfg%get("curvi","l0_ducros",self%l0_ducros)
      endif
    endif

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
    call self%cfg%get("controls","restart_type",self%restart_type)
    if(.not. any(self%restart_type == [0,1,2])) call fail_input_any("restart_type must be 0, 1, or 2")
    call self%cfg%get("grid","grid_type",grid_type)
    if(.not. any(grid_type == [1,2,3,4,5,6])) call fail_input_any("grid_type must be 1, 2, 3, 4, 5 or 6")
    rebuild_ghost = .false.
    if (grid_type == GRID_FROMFILE) rebuild_ghost = .true.
    if (self%restart_type > 0 .and. grid_dim == 1) grid_type = GRID_FROMFILE

    call self%cfg%get("grid","nxmax",nxmax)
    call self%cfg%get("grid","nymax",nymax)
    call self%cfg%get("grid","nzmax",nzmax)
    call self%cfg%get("grid","ng",ng)
    call self%cfg%get("grid","domain_size_x",domain_size_x)
    call self%cfg%get("grid","domain_size_y",domain_size_y)
    call self%cfg%get("grid","domain_size_z",domain_size_z)
    call self%cfg%get("grid","metrics_order",metrics_order)
    if(.not. any(metrics_order == [2,4,6,8])) call fail_input_any("metrics_order must be 2, 4, 6 or 8")
    ystag = 0
    if (self%cfg%has_key("grid","ystag")) then
      call self%cfg%get("grid","ystag",ystag)
    endif
    if(.not. any(ystag == [0,1])) call fail_input_any("ystag must be 0 or 1")
    ystaggering = .false.
    if (ystag>0) ystaggering = .true.
    if (self%channel_case) ystaggering = .true.

    call self%cfg%get("numerics","ep_order",self%ep_order)
    if (self%ep_order /= metrics_order) call fail_input_any("ep_order must be equal to metrics_order")

    select case (grid_type)
    case (GRID_FROMFILE)
    case (GRID_UNIFORM)
    case (GRID_CHA)
      if (grid_dim == 1) then
        allocate(grid_vars(4))
        call self%cfg%get("grid","jbgrid",rtemp) ; grid_vars(1) = rtemp
        call self%cfg%get("grid","dyptarget",rtemp) ; grid_vars(2) = rtemp
        call self%cfg%get("grid","ysmoosteps",rtemp) ; grid_vars(3) = rtemp
        call self%cfg%get("flow","Reynolds",rtemp) ; grid_vars(4) = rtemp
      endif
      if (grid_dim == 2) then
        allocate(grid_vars(4))
        call self%cfg%get("grid","dyptarget",rtemp) ; grid_vars(1) = rtemp
        call self%cfg%get("curvi","Reynolds_friction",rtemp) ; grid_vars(2) = rtemp
        call self%cfg%get("curvi","angle",rtemp) ; grid_vars(3) = rtemp
        call self%cfg%get("curvi","R_curv",rtemp) ; grid_vars(4) = rtemp
      endif
    case (GRID_BL,GRID_SBLI)
      if(grid_dim == 1) then
        allocate(grid_vars(6))
        call self%cfg%get("grid","jbgrid",rtemp) ; grid_vars(1) = rtemp
        call self%cfg%get("grid","dyptarget",rtemp) ; grid_vars(2) = rtemp
        call self%cfg%get("grid","nywr",rtemp) ; grid_vars(3) = rtemp
        call self%cfg%get("grid","lywr",rtemp) ; grid_vars(4) = rtemp
        call self%cfg%get("flow","Reynolds",rtemp) ; grid_vars(5) = rtemp
        call self%cfg%get("grid","ysmoosteps",rtemp) ; grid_vars(6) = rtemp
      endif
      if(grid_dim == 2) then
      endif
    case(GRID_AIRFOIL)
    end select

    call self%grid%initialize(periodic, nxmax, nymax, nzmax, ng, grid_type, domain_size_x,&
    & domain_size_y, domain_size_z, self%l0, grid_vars, metrics_order, rebuild_ghost, ystaggering,&
    & grid_dim, grid2d_par)

    self%nv = 5
    self%nv_aux = 10
    if(self%grid%grid_dim == 1) then
      self%nv_stat = 70
    elseif(self%grid%grid_dim == 2) then
      self%nv_stat = 70
    endif
    self%nv_stat_3d = 14

    self%save_stat_z = .false.
    if (self%cfg%has_key("output","save_stat_z")) then
      call self%cfg%get("output","save_stat_z",itemp)
      if (itemp > 0) self%save_stat_z = .true.
    endif
    istat = create_folder("AVGZ")

    call self%set_fluid_prop()
    call self%set_flow_params()

    self%nkeep = 0
    if (self%cfg%has_key("numerics","nkeep")) then
      call self%cfg%get("numerics","nkeep",self%nkeep)
      if(self%grid%grid_dim == 2 .and. self%nkeep /= -1) call fail_input_any("nkeep >= 0 not support&
      &ed for curvilinear grids. Only -1 supported")
    endif
    call self%cfg%get("numerics","weno_scheme",self%weno_scheme)
    if(.not. any(self%weno_scheme == [1,2,3,4])) call fail_input_any("weno_scheme must be 1, 2, 3 or 4")
    call self%cfg%get("numerics","weno_version",self%weno_version)
    if(.not. any(self%weno_version == [0,1])) call fail_input_any("weno_version must be 0, or 1")
    if(self%weno_scheme /= 3 .and. self%weno_version == 1) call fail_input_any("weno_version 1 supported only for weno_scheme = 3")
    self%flux_splitting = 0
    if (self%cfg%has_key("numerics","flux_splitting")) then
      call self%cfg%get("numerics","flux_splitting",self%flux_splitting)
    endif
    if(.not. any(self%flux_splitting == [0,1])) call fail_input_any("flux_splitting must be 0, or 1")
    call self%cfg%get("numerics","visc_order",self%visc_order)
    if (self%visc_order /= self%ep_order) call fail_input_any("viscous order must be equal to ep_order")
    self%count_weno_control = 0
    if (self%cfg%has_key("numerics","count_weno_control")) then
      call self%cfg%get("numerics","count_weno_control",self%count_weno_control)
    endif

    call self%cfg%get("numerics","conservative_viscous",self%conservative_viscous)
    if(.not. any(self%conservative_viscous == [0,1])) call fail_input_any("conservative_viscous must be 0, or 1")
    allocate(self%supported_orders(1:4))
    self%supported_orders = [2,4,6,8]
    do i=1,size(self%supported_orders)
      order = self%supported_orders(i)
      call self%grid%get_deriv_coeffs(d1_temp, d2_temp, order, order)
      self%coeff_deriv1(:,i) = d1_temp
      self%coeff_deriv2(:,i) = d2_temp
    enddo
    call self%cfg%get("numerics","sensor_threshold",self%sensor_threshold)
    if(self%masterproc) then
      if(self%sensor_threshold < 0._rkind) print*,'Running full WENO'
      if(self%sensor_threshold > 1._rkind) print*,'Running full central'
    endif
    if (grid_dim==2) then
      self%ortho = 1
      if (self%cfg%has_key("numerics","ortho")) then
        call self%cfg%get("numerics","ortho",self%ortho)
      endif
      if (self%masterproc) write(*,*) 'Ortogonality of grid: ',self%ortho
      if(.not. any(self%ortho == [0,1])) call fail_input_any("ortho of grid must be 0, or 1")
    endif

    call self%cfg%get("output","enable_plot3d",self%enable_plot3d)

    if (self%cfg%has_key("output","enable_stat_time")) then
      call self%cfg%get("output","enable_stat_time",self%enable_stat_time)
    else
      self%enable_stat_time = 1
    endif

    if (self%cfg%has_key("output","enable_stat_3d")) then
      call self%cfg%get("output","enable_stat_3d",self%enable_stat_3d)
    else
      self%enable_stat_3d = 0
    endif

    self%enable_tspec = 0
    if (self%cfg%has_key("output","enable_tspec")) then
      call self%cfg%get("output","enable_tspec",self%enable_tspec)
    endif
    if (self%enable_tspec>0) then
      call self%cfg%get("output","dt_tspec",self%dt_tspec)
      call self%cfg%get("output","ndft_tspec",self%ndft_tspec)
      call self%cfg%get("output","jg_tspec",self%jg_tspec)
      if (self%masterproc) write(*,*) 'Time spectra enabled dt_tspec =',self%dt_tspec, ' - ndft_tspec=',self%ndft_tspec
      self%wss_tspec = 0._rkind
      do l=1,self%ndft_tspec
        wind = 1._rkind
        self%wss_tspec = self%wss_tspec+wind*wind
      enddo
    endif

    call self%cfg%get("mpi","x_split",mpi_split_x) ; mpi_splits(1) = mpi_split_x
    call self%cfg%get("mpi","y_split",mpi_split_y) ; mpi_splits(2) = mpi_split_y
    call self%cfg%get("mpi","z_split",mpi_split_z) ; mpi_splits(3) = mpi_split_z
    if(mod(nxmax, mpi_split_x) /= 0) call fail_input_any("nxmax not muplitple of mpi_split_x")
    if(mpi_split_y /= 1) call fail_input_any("mpi_split_y must be 1")
    if(mod(nzmax, mpi_split_z) /= 0) call fail_input_any("nzmax not muplitple of mpi_split_z")

    call self%field%initialize(self%grid, self%nv, mpi_splits, save_metrics_2d)

    call self%alloc()

    call self%field%correct_bc(self%bctags)
    call self%field%correct_bc(self%bctags_nr)

    call self%bc_preproc()

    self%correct_bound_ord = 0
    if (self%cfg%has_key("bc","correct_bound_ord")) then
      call self%cfg%get("bc","correct_bound_ord",self%correct_bound_ord)
    endif
    if(.not. any(self%correct_bound_ord == [0,1])) call fail_input_any("correct_bound_ord must be 0, or 1")
    if (self%correct_bound_ord>0) call self%correct_bc_order()

    if (grid_dim==2) then

      if (self%flow_init == 5) then
        call self%correct_airfoil_order()
      else
        self%lmax_tag(:) = max(self%ep_order/2, self%weno_scheme)
        self%vis_tag(:) = self%visc_order/2
      endif


      self%jweno = self%grid%nymax
      if (self%cfg%has_key("curvi","jweno")) then
        call self%cfg%get("curvi","jweno",self%jweno)
        if (self%jweno==0) self%jweno = self%grid%nymax
        if (self%jweno < self%grid%nymax) then
          if (self%masterproc) write(*,*) 'Weno always active starting from j = ', self%jweno
        endif
      endif

      self%jfilter = self%grid%nymax
      self%ifilter = 0
      if (self%cfg%has_key("curvi","ifilter")) then
        call self%cfg%get("curvi","ifilter",self%ifilter)
        if (self%ifilter > 0) then
          if (self%cfg%has_key("curvi","jfilter")) then
            call self%cfg%get("curvi","jfilter",self%jfilter)
            if (self%jfilter==0) self%jfilter = self%grid%nymax
          endif
          if (self%grid%ng < 4) then
            if (self%masterproc) print*, 'Error: filter supported with al least 4 ghost nodes'
            call mpi_abort(mpi_comm_world,99,self%mpi_err)
          endif
          if (self%masterproc) write(*,*) 'Filter enabled starting from j = ', self%jfilter
          self%coeff_filter(0) = 3565._rkind/20736._rkind
          self%coeff_filter(1) = 3091._rkind/12960._rkind
          self%coeff_filter(2) = 1997._rkind/25920._rkind
          self%coeff_filter(3) = 149._rkind /12960._rkind
          self%coeff_filter(4) = 107._rkind/103680._rkind
        endif
      endif

      self%enable_limiter = 0
      if (self%cfg%has_key("curvi","enable_limiter")) then
        call self%cfg%get("curvi","enable_limiter",self%enable_limiter)
        if (self%enable_limiter > 0) then
          if (self%masterproc) write(*,*) 'Limiter enabled'
        endif
      endif
      self%rho_lim = 0.0_rkind
      if (self%cfg%has_key("curvi","rho_lim")) then
        call self%cfg%get("curvi","rho_lim",self%rho_lim)
      endif
      self%tem_lim = 0.0_rkind
      if (self%cfg%has_key("curvi","tem_lim")) then
        call self%cfg%get("curvi","tem_lim",self%tem_lim)
      endif
      self%rho_lim_rescale = 1.0_rkind
      if (self%cfg%has_key("curvi","rho_lim_rescale")) then
        call self%cfg%get("curvi","rho_lim_rescale",self%rho_lim_rescale)
      endif
      self%tem_lim_rescale = 1.0_rkind
      if (self%cfg%has_key("curvi","tem_lim_rescale")) then
        call self%cfg%get("curvi","tem_lim_rescale",self%tem_lim_rescale)
      endif
      if (self%enable_limiter > 0) then
        if (self%masterproc) then
          write(*,*) 'Limited minimum rho :',self%rho_lim
          write(*,*) 'Limited minimum tem :',self%tem_lim
          write(*,*) 'Limited minimum rho rescale :',self%rho_lim_rescale
          write(*,*) 'Limited minimum tem rescale :',self%tem_lim_rescale
        endif
      endif

      self%enable_forces_runtime = 0
      if (self%cfg%has_key("curvi","enable_forces_runtime")) then
        call self%cfg%get("curvi","enable_forces_runtime",self%enable_forces_runtime)
      endif

      self%theta_threshold = 100000._rkind
      if (self%cfg%has_key("curvi","theta_threshold")) then
        call self%cfg%get("curvi","theta_threshold",self%theta_threshold)
        if (self%masterproc) write(*,*) 'theta_threshold (bad cell theta level) :',self%theta_threshold
      endif
    endif

    call self%field%compute_local_grid_metrics(self%coeff_deriv1, self%ep_ord_change, self%lmax_tag)

    call self%cfg%get("numerics","rk_type",self%rk_type)
    call self%runge_kutta_initialize()
    call self%cfg%get("controls","cfl",self%cfl)
    call self%cfg%get("controls","num_iter",self%num_iter)

    self%io_type_r = 2
    if (self%cfg%has_key("output","io_type_r")) then
      call self%cfg%get("output","io_type_r",self%io_type_r)
    endif
    self%io_type_w = 2
    if (self%cfg%has_key("output","io_type_w")) then
      call self%cfg%get("output","io_type_w",self%io_type_w)
    endif
    if(.not. any(self%io_type_r == [0,1,2])) call fail_input_any("io_type_r must be 0, 1 or 2")
    if(.not. any(self%io_type_w == [0,1,2])) call fail_input_any("io_type_w must be 0, 1 or 2")
    self%restart_from_flow = -1
    if (self%cfg%has_key("controls","restart_from_flow")) then
      call self%cfg%get("controls","restart_from_flow",self%restart_from_flow)
    endif

    if (self%cfg%has_key("insitu","enable_insitu")) then
      call self%cfg%get("insitu","enable_insitu",self%enable_insitu)
      self%insitu_platform = "catalyst-v2"
      if (self%cfg%has_key("insitu","insitu_platform")) then
        call self%cfg%get("insitu","insitu_platform",self%insitu_platform)
      endif
      if (self%enable_insitu>0) then
        if (self%masterproc) write(*,*) 'Insitu-platform :',self%insitu_platform
      endif
    endif
    call self%initial_conditions()
    select case(self%restart_type)
    case(0)
      self%time0 = 0._rkind
      self%icyc0 = 0
      self%itav = 0
      self%itslice_p3d = 0
      self%istore = 1
      self%time_from_last_rst = 0._rkind
      self%time_from_last_write = 0._rkind
      self%time_from_last_stat = 0._rkind
      self%time_from_last_slice = 0._rkind
      self%time_from_last_probe = 0._rkind
      self%time_from_last_slice_p3d = 0._rkind
      if(self%enable_insitu > 0) then
        self%i_insitu = 1
        self%time_from_last_insitu = 0._rkind
        self%time_insitu = 0._rkind
      endif
      self%w_stat = 0._rkind
      if (self%enable_stat_3d>0) self%w_stat_3d = 0._rkind
      if (self%enable_tspec>0) then
        self%it_win1_tspec = 1
        self%it_win2_tspec = 1-self%ndft_tspec/2
        self%it_nwin_tspec = 0
        self%w_psd_avg_tspec = 0._rkind
        self%w_tspec = 0._rkind
        self%time_from_last_tspec = 0._rkind
      endif
    case(1)
      if (self%restart_from_flow>=0) then
        call self%field%read_plot3d(istore=self%restart_from_flow)
        call self%fromaux_to_w()
      else
        if (self%io_type_r==1) call self%field%read_field_serial()
        if (self%io_type_r==2) call self%field%read_field()
      endif
      call self%read_field_info()
      self%itav = 0
      self%itslice_p3d = 0
      self%w_stat = 0._rkind
      if (self%enable_stat_3d>0) self%w_stat_3d = 0._rkind
      if (self%enable_tspec>0) then
        self%it_win1_tspec = 1
        self%it_win2_tspec = 1-self%ndft_tspec/2
        self%it_nwin_tspec = 0
        self%w_psd_avg_tspec = 0._rkind
        self%w_tspec = 0._rkind
        self%time_from_last_tspec = 0._rkind
      endif
    case(2)
      if (self%io_type_r==1) then
        if (self%restart_from_flow>=0) then
          call self%field%read_plot3d(istore=self%restart_from_flow)
          call self%fromaux_to_w()
        else
          call self%field%read_field_serial()
        endif
        if (self%enable_stat_time>0) call self%read_stats_serial()
        if (self%enable_stat_3d>0) call self%read_stats_3d_serial()
      endif
      if (self%io_type_r==2) then
        if (self%restart_from_flow>=0) then
          call self%field%read_plot3d(istore=self%restart_from_flow)
          call self%fromaux_to_w()
        else
          call self%field%read_field()
        endif
        if (self%enable_stat_time>0) call self%read_stats()
        if (self%enable_stat_3d>0) call self%read_stats_3d()
      endif

      call self%read_field_info()

      if (self%enable_tspec == 1) then
        self%it_win1_tspec = 1
        self%it_win2_tspec = 1-self%ndft_tspec/2
        self%it_nwin_tspec = 0
        self%w_psd_avg_tspec = 0._rkind
        self%w_tspec = 0._rkind
        self%time_from_last_tspec = 0._rkind
      elseif (self%enable_tspec > 1) then
        call self%read_tspec()
      endif
    endselect

    if (self%cfg%has_key("curvi","enable_sponge")) then
      call self%cfg%get("curvi","enable_sponge",itemp) ; self%enable_sponge = itemp
    else
      self%enable_sponge = 0
    endif
    if (self%enable_sponge > 0) call self%init_sponge

    self%recyc = .false.
    if (self%field%ncoords(1) == 0) then
      if(self%bctags(1) == 10) then
        self%recyc = .true.
      endif
    endif
    call mpi_bcast(self%recyc,1,mpi_logical,0,self%field%mp_cartx,self%mpi_err)
    if (self%recyc) call self%recyc_prepare()

    self%time = self%time0
    self%icyc = self%icyc0

    call self%cfg%get("output","dtsave",self%dtsave)
    call self%cfg%get("output","dtsave_restart",self%dtsave_restart)
    call self%cfg%get("output","dtstat",self%dtstat)
    call self%cfg%get("output","enable_plot3d",self%enable_plot3d)
    call self%cfg%get("output","enable_vtk",self%enable_vtk)
    call self%cfg%get("output","dtslice",self%dtslice)
    if (self%cfg%has_key("output","list_aux_slice")) then
      call self%cfg%get("output","list_aux_slice",self%list_aux_slice)
      self%num_aux_slice = size(self%list_aux_slice)
    else
      self%list_aux_slice = [1, 2, 3, 4, 5, 6]
      self%num_aux_slice = 6
    endif
    if (self%cfg%has_key("output","dtprobe")) then
      call self%cfg%get("output","dtprobe",self%dtprobe)
    else
      self%dtprobe = self%dtslice
    endif
    if (self%cfg%has_key("output","dtslice_p3d")) then
      call self%cfg%get("output","dtslice_p3d",self%dtslice_p3d)
    else
      self%dtslice_p3d = self%dtslice
    endif
    call self%cfg%get("output","igslice",self%igslice)
    call self%cfg%get("output","jgslice",self%jgslice)
    call self%cfg%get("output","kgslice",self%kgslice)
    call self%slice_prepare(self%igslice, self%jgslice, self%kgslice, self%islice, self%jslice, self%kslice, alloc_aux=.true.)
    call self%probe_prepare()
    call self%probe_compute_coeff()

    self%enable_slice_plot3d = 0
    if (self%cfg%has_key("output","enable_slice_plot3d")) then
      call self%cfg%get("output","enable_slice_plot3d",self%enable_slice_plot3d)
    endif
    if(self%enable_slice_plot3d > 0) then
      call self%cfg%get("output","igslice_p3d",self%igslice_p3d)
      call self%cfg%get("output","jgslice_p3d",self%jgslice_p3d)
      call self%cfg%get("output","kgslice_p3d",self%kgslice_p3d)
      call self%slice_prepare(self%igslice_p3d, self%jgslice_p3d, self%kgslice_p3d, self%islice_p3d,&
      & self%jslice_p3d, self%kslice_p3d, alloc_aux = .false.)
    endif

    call self%cfg%get("output","print_control",self%print_control)
    call self%cfg%get("controls","iter_dt_recompute",self%iter_dt_recompute)


    if (self%enable_insitu>0) call self%insitu_initialize()

  endsubroutine initialize

  subroutine init_sponge(self)
    class(equation_singleideal_object), intent(inout) :: self

    integer :: i, ii, j, l, itemp, ile_rank
    real(rkind) :: rad, radmax, radmin
    logical :: wfar_check

    allocate(self%f_sponge(self%field%ny))

    associate(nx => self%field%nx, ny => self%field%ny, ncoords => self%field%ncoords,&
    & ierr => self%mpi_err, nxmax => self%grid%nxmax, xc2 => self%field%xc2, yc2 => self%field%yc2,&
    & f_sponge => self%f_sponge, j_sponge => self%j_sponge, wfar => self%wfar, ile_l => self%field%ile_l,&
    & ile_rank_x => self%field%ile_rank_x)

      if (self%cfg%has_key("curvi","j_sponge")) then
        call self%cfg%get("curvi","j_sponge",j_sponge)
      else
        j_sponge = nint(0.9375_rkind*real(ny))
      endif
      if (self%cfg%has_key("curvi","k_sponge")) then
        call self%cfg%get("curvi","k_sponge",self%k_sponge)
      else
        self%k_sponge = .1_rkind
      endif

      f_sponge = .0_rkind
      if (ncoords(1) == ile_rank_x .and. ncoords(3)==0) then
        radmax = ((xc2(ile_l,ny )-xc2(ile_l,1))**2 + (yc2(ile_l,ny )-yc2(ile_l,1))**2)**0.5
        radmin = ((xc2(ile_l,j_sponge)-xc2(ile_l,1))**2 + (yc2(ile_l,j_sponge)-yc2(ile_l,1))**2)**0.5
        rad = radmin
        write(*,*) 'C-mesh radius =', radmax
        write(*,*) 'Absorbing layer is active on ymax starting from radius =', radmin
        do j=j_sponge,ny
          rad = rad + self%field%etamod(ile_l,j)
          f_sponge(j) = self%k_sponge*((rad-radmin)/(radmax-radmin))**3
          if (f_sponge(j) > 1._rkind) f_sponge(j) = 1._rkind
        enddo
        open(unit=11,file='sponge_strength.dat',form='formatted')
        do j=1,ny
          write(11,*) j, abs(xc2(ile_l,j)), f_sponge(j)
        enddo
        close(11)
      endif
      call mpi_cart_rank(self%field%mp_cart, [ile_rank_x, 0, 0], ile_rank, ierr)
      call mpi_bcast(f_sponge, ny, mpi_prec, ile_rank, self%field%mp_cart, ierr)

      if (self%enable_sponge > 1) then
        allocate(self%wfar_glob(self%grid%nxmax,self%nv))
        inquire(file='wfar.bin', exist=wfar_check)
        if (wfar_check) then
          if (self%masterproc) write(*,*) 'Reading far-field solution from file'
          open(11,file='wfar.bin',form='unformatted')
          read(11) self%wfar_glob(1:nxmax,1:self%nv)
          close(11)
          do l=1,self%nv
            ii = nx*ncoords(1)
            do i=1,nx
              wfar(i,l) = self%wfar_glob(ii+i,l)
            enddo
          enddo
        else
          if (self%masterproc) write(*,*) 'Storing current far-field solution to file'
          do l=1,self%nv
            do i=1,nx
              wfar(i,l) = self%field%w(i,ny,1,l)
            enddo
            call mpi_gather(wfar(:,l),nx,mpi_prec,self%wfar_glob(:,l),nx,mpi_prec,0,self%field%mp_cartx,self%mpi_err)
          enddo
          if (self%masterproc) then
            open(11,file='wfar.bin',form='unformatted')
            write(11) self%wfar_glob(1:nxmax,1:self%nv)
            close(11)
          endif
        endif
      endif

    endassociate
  endsubroutine init_sponge

  subroutine correct_bc_order(self)
    class(equation_singleideal_object), intent(inout) :: self
    integer :: i,j,k
    integer :: stencil_size
    associate(nx => self%field%nx, ny => self%field%ny,nz => self%field%nz, ng => self%grid%ng,&
    & ep_order => self%ep_order, weno_scheme => self%weno_scheme, ncoords => self%field%ncoords,&
    & bctags => self%bctags, nblocks => self%field%nblocks)

      stencil_size = max(ep_order/2, weno_scheme)
      if ((ncoords(1)==0).and.(any(bctags(1:2)/=0))) then
        self%ep_ord_change(0,:,:,1) = -stencil_size+1
        do i=1,stencil_size-1
          self%ep_ord_change(i,:,:,1) = -stencil_size+i
        enddo
      endif
      if ((ncoords(1)==(nblocks(1)-1)).and.(any(bctags(1:2)/=0))) then
        self%ep_ord_change(nx,:,:,1) = -stencil_size+1
        do i=1,stencil_size-1
          self%ep_ord_change(nx-i,:,:,1) = -stencil_size+i
        enddo
      endif
      if ((ncoords(2)==0).and.(any(bctags(3:4)/=0))) then
        self%ep_ord_change(:,0,:,2) = -stencil_size+1
        do j=1,stencil_size-1
          self%ep_ord_change(:,j,:,2) = -stencil_size+j
        enddo
      endif
      if ((ncoords(2)==(nblocks(2)-1)).and.(any(bctags(3:4)/=0))) then
        self%ep_ord_change(:,ny,:,2) = -stencil_size+1
        do j=1,stencil_size-1
          self%ep_ord_change(:,ny-j,:,2) = -stencil_size+j
        enddo
      endif
      if ((ncoords(3)==0).and.(any(bctags(5:6)/=0))) then
        self%ep_ord_change(:,:,0,3) = -stencil_size+1
        do k=1,stencil_size-1
          self%ep_ord_change(:,:,k,3) = -stencil_size+k
        enddo
      endif
      if ((ncoords(3)==(nblocks(3)-1)).and.(any(bctags(5:6)/=0))) then
        self%ep_ord_change(:,:,nz,3) = -stencil_size+1
        do k=1,stencil_size-1
          self%ep_ord_change(:,:,nz-k,3) = -stencil_size+k
        enddo
      endif

    endassociate
  end subroutine correct_bc_order

  subroutine insitu_do_catalyst_initialization(self)
    class(equation_singleideal_object), intent(inout) :: self
    type(C_PTR) :: params, about, exec, final, results
    integer, parameter :: f64 = selected_real_kind(8)
    integer(kind(catalyst_status)) :: code
    integer :: exit_code
    character(128) :: PARAVIEW_IMPL_DIR=""

    params = catalyst_conduit_node_create()
    call catalyst_conduit_node_set_path_char8_str(params, "catalyst/scripts/script0", trim(adjustl(self%vtkpipeline)))
    call catalyst_conduit_node_set_path_char8_str(params, "catalyst_load/implementation", "paraview")
    call catalyst_conduit_node_set_path_char8_str(params, "catalyst_load/search_paths/paraview", PARAVIEW_IMPL_DIR)
    code = c_catalyst_initialize(params)
    if (code /= catalyst_status_ok) then
      write (error_unit, *) "failed to initialize: ", code
      exit_code = 1
    end if
    call catalyst_conduit_node_destroy(params)
  endsubroutine insitu_do_catalyst_initialization

  subroutine correct_airfoil_order(self)
    class(equation_singleideal_object), intent(inout) :: self
    integer :: i,j,k,l, ii, lmax, mm
    integer, dimension(1-self%grid%ng:self%grid%nxmax+self%grid%ng) :: lmax_tagg,vis_tagg
    associate(masterproc => self%masterproc, nx => self%field%nx, nxmax => self%grid%nxmax,&
    & ng => self%grid%ng, ep_order => self%ep_order, weno_scheme => self%weno_scheme, ncoords => self%fie&
    &ld%ncoords, visc_order => self%visc_order, xwall => self%grid%xwall, ywall => self%grid%ywall,&
    & wall_tagg => self%grid%wall_tagg, ite => self%grid%ite, itu => self%grid%itu, lmax_tag => self%lmax&
    &_tag, vis_tag => self%vis_tag)

      lmax = max(ep_order/2, weno_scheme)
      mm = visc_order/2

      lmax_tagg = lmax
      do l=1,lmax
        lmax_tagg(ite-l) = l
        lmax_tagg(itu+l-1) = l
      enddo

      vis_tagg = mm
      do l=1,mm
        vis_tagg(ite-l) = l
        vis_tagg(itu+l) = l
      enddo

      if (masterproc) then
        open(18,file='tagged_wall_orders.dat')
        do i=1-ng,nxmax+ng
          write(18,'(I0,2X,ES24.14,2X,ES24.14,2X,3(I0,2X))') i,xwall(i),ywall(i),wall_tagg(i),lmax_tagg(i),vis_tagg(i)
        enddo
        close(18)
      endif

      ii = nx*ncoords(1)
      do i=1-ng,nx+ng
        lmax_tag(i) = lmax_tagg(ii+i)
        vis_tag(i) = vis_tagg(ii+i)
      enddo

    endassociate
  end subroutine correct_airfoil_order

  subroutine insitu_initialize(self)
    use tcp
    class(equation_singleideal_object), intent(inout) :: self
    integer :: counter, i, j, k
    call self%cfg%get("insitu","vtkpipeline",self%vtkpipeline)

    if(self%insitu_platform == "catalyst-v1") then
      call insitu_start(self%fcoproc,trim(adjustl(self%vtkpipeline)),self%masterproc)
    elseif(self%insitu_platform == "catalyst-v2") then
      call self%insitu_do_catalyst_initialization()
    else
      write(*,*) 'insitu_platform not available. Aborting!'
      call MPI_ABORT(mpi_comm_world,-40,self%mpi_err)
    endif

    call self%cfg%get("insitu","dt_insitu",self%dt_insitu)
    call self%cfg%get("insitu","perc_ny_cut",self%perc_ny_cut)
    call self%cfg%get("insitu","freeze_intervals",self%freeze_intervals)
    if (allocated(self%freeze_intervals)) then
      if (size(self%freeze_intervals) > 0) then
        self%enable_freeze_intervals = .true.
      endif
    endif
    call self%insitu_define_limits()
    call self%insitu_allocate()
    associate( x => self%field%x, y => self%field%y, z => self%field%z, xc2 => self%field%xc2,&
    & yc2 => self%field%yc2, nrank => self%myrank, nproc => self%nprocs, nxsl_ins => self%nxsl_ins,&
    & nxel_ins => self%nxel_ins, nysl_ins => self%nysl_ins, nyel_ins => self%nyel_ins,&
    & nzsl_ins => self%nzsl_ins, nzel_ins => self%nzel_ins, nxs_ins => self%nxs_ins, nxe_ins => self%nxe_&
    &ins, nys_ins => self%nys_ins, nye_ins => self%nye_ins, nzs_ins => self%nzs_ins, nze_ins => self%nze_&
    &ins, nxstartg => self%nxstartg, nxendg => self%nxendg, nystartg => self%nystartg,&
    & nyendg => self%nyendg, nzstartg => self%nzstartg, nzendg => self%nzendg, xyzc => self%xyzc,&
    & points_x => self%points_x, points_y => self%points_y, points_z => self%points_z,&
    & n_points_x => self%n_points_x, n_points_y => self%n_points_y, n_points_z => self%n_points_z,&
    & n_points => self%n_points)
      if(self%insitu_platform == "catalyst-v1") then
        call createcpstructureddata(nxel_ins-nxsl_ins+1,nyel_ins-nysl_ins+1,nzel_ins-nzsl_ins+1,&
        &nxs_ins,nys_ins,nzs_ins,nxe_ins,nye_ins,nze_ins,nxendg-nxstartg+1,nyendg-nystartg+1,nzendg-nzstartg+&
        &1,x(nxsl_ins:nxel_ins),y(nysl_ins:nyel_ins),z(nzsl_ins:nzel_ins) ,nrank,nproc ,"3d_struct"//&
        &c_null_char,xyzc)
      elseif(self%insitu_platform == "catalyst-v2") then
        counter = 1
        n_points_x = nxel_ins-nxsl_ins+1
        n_points_y = nyel_ins-nysl_ins+1
        n_points_z = nzel_ins-nzsl_ins+1
        n_points = n_points_x * n_points_y * n_points_z
        if(self%grid%grid_dim == 1) then
          do k=nzsl_ins,nzel_ins
            do j=nysl_ins,nyel_ins
              do i=nxsl_ins,nxel_ins
                points_x(counter) = x(i)
                points_y(counter) = y(j)
                points_z(counter) = z(k)
                counter = counter + 1
              enddo
            enddo
          enddo
        elseif(self%grid%grid_dim == 2) then
          do k=nzsl_ins,nzel_ins
            do j=nysl_ins,nyel_ins
              do i=nxsl_ins,nxel_ins
                points_x(counter) = xc2(i,j)
                points_y(counter) = yc2(i,j)
                points_z(counter) = z(k)
                counter = counter + 1
              enddo
            enddo
          enddo
        endif
      endif
    endassociate
  endsubroutine insitu_initialize

  subroutine insitu_define_limits(self)
    class(equation_singleideal_object), intent(inout) :: self
    integer :: ngg, j
    integer :: i, ii
    real(rkind) :: ycut
    call self%cfg%get("insitu","aux_list",self%aux_list)
    self%n_aux_list = size(self%aux_list)
    allocate(self%aux_list_name(self%n_aux_list))
    call self%cfg%get("insitu","add_list",self%add_list)
    self%n_add_list = size(self%add_list)
    allocate(self%add_list_name(self%n_add_list))
    associate( aux_list => self%aux_list, aux_list_name => self%aux_list_name, mpi_err => self%mpi_e&
    &rr, n_aux_list => self%n_aux_list, add_list => self%add_list, add_list_name => self%add_list_name,&
    & n_add_list => self%n_add_list, ngm => self%ngm, ngp => self%ngp, yg => self%grid%yg,&
    & nxmax => self%grid%nxmax, nymax => self%grid%nymax, nzmax => self%grid%nzmax, perc_ny_cut => self%p&
    &erc_ny_cut, ny_cut => self%ny_cut, nxsl_ins => self%nxsl_ins, nxel_ins => self%nxel_ins,&
    & nysl_ins => self%nysl_ins, nyel_ins => self%nyel_ins, nzsl_ins => self%nzsl_ins,&
    & nzel_ins => self%nzel_ins, nxs_ins => self%nxs_ins, nxe_ins => self%nxe_ins, nys_ins => self%nys_in&
    &s, nye_ins => self%nye_ins, nzs_ins => self%nzs_ins, nze_ins => self%nze_ins, nxstartg => self%nxsta&
    &rtg, nxendg => self%nxendg, nystartg => self%nystartg, nyendg => self%nyendg, nzstartg => self%nzsta&
    &rtg, nzendg => self%nzendg, npsi => self%npsi, npsi_pv => self%npsi_pv, nx => self%field%nx,&
    & ny => self%field%ny, nz => self%field%nz, ncoords => self%field%ncoords )
      do i=1,n_aux_list
        ii = aux_list(i)
        if (ii == 1) aux_list_name(i) = "density"
        if (ii == 2) aux_list_name(i) = "u"
        if (ii == 3) aux_list_name(i) = "v"
        if (ii == 4) aux_list_name(i) = "w"
        if (ii == 5) aux_list_name(i) = "h"
        if (ii == 6) aux_list_name(i) = "T"
        if (ii == 7) aux_list_name(i) = "viscosity"
        if (ii > 7) call MPI_ABORT(mpi_comm_world,-3,mpi_err)
      enddo
      do i=1,n_add_list
        ii = add_list(i)
        if (ii == 1) add_list_name(i) = "div"
        if (ii == 2) add_list_name(i) = "abs_omega"
        if (ii == 3) add_list_name(i) = "ducros"
        if (ii == 4) add_list_name(i) = "swirling_strength"
        if (ii == 5) add_list_name(i) = "schlieren"
        if (ii > 5) call MPI_ABORT(mpi_comm_world,-3,mpi_err)
      enddo
      ngg = 1
      if (ngg == 0) then
        ngm = 0
        ngp = 1
      else
        ngm = ngg
        ngp = ngg
      endif
      ny_cut=nymax
      if(perc_ny_cut > 0._rkind) then
        ycut = yg(nymax)*(1.-perc_ny_cut/100.)
        do j=1,nymax
          if (yg(j)>ycut) then
            ny_cut = j
            exit
          endif
        enddo
      endif
      nxsl_ins = 1-ngm
      nxel_ins = nx+ngp
      nysl_ins = 1
      nyel_ins = ny_cut
      if(self%grid%ndim == 2) then
        nzsl_ins = 1
        nzel_ins = nz
      else
        nzsl_ins = 1-ngm
        nzel_ins = nz+ngp
      endif
      nxs_ins = nxsl_ins + ncoords(1)*nx
      nxe_ins = nxel_ins + ncoords(1)*nx
      nys_ins = nysl_ins + ncoords(2)*ny
      nye_ins = nyel_ins + ncoords(2)*ny
      nzs_ins = nzsl_ins + ncoords(3)*nz
      nze_ins = nzel_ins + ncoords(3)*nz
      nxstartg = 1-ngm
      nxendg = nxmax+ngp
      nystartg = 1
      nyendg = ny_cut
      nzstartg = 1-ngm
      nzendg = nzmax+ngp
      npsi = size(self%add_list)
      npsi_pv = npsi+size(self%aux_list)
      if (self%masterproc) write(*,*) "Insitu variable: ", aux_list_name, add_list_name
    endassociate
  endsubroutine insitu_define_limits

  subroutine insitu_allocate(self)
    class(equation_singleideal_object), intent(inout) :: self
    integer :: n_points_insitu
    associate( nxsl_ins => self%nxsl_ins, nxel_ins => self%nxel_ins, nysl_ins => self%nysl_ins,&
    & nyel_ins => self%nyel_ins, nzsl_ins => self%nzsl_ins, nzel_ins => self%nzel_ins, npsi => self%npsi,&
    & npsi_pv => self%npsi_pv, nx => self%field%nx, ny => self%field%ny, nz => self%field%nz,&
    & ng => self%grid%ng, insitu_platform => self%insitu_platform )
      if (npsi_pv > 0) allocate(self%psi_pv(nxsl_ins:nxel_ins,nysl_ins:nyel_ins,nzsl_ins:nzel_ins,npsi_pv))
      if (npsi > 0) allocate(self%psi(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,npsi))
      if(insitu_platform == "catalyst-v1") then
        allocate(self%xyzc(3*(nx+2)*(ny)*(nz+2)))
      elseif(insitu_platform == "catalyst-v2") then
        n_points_insitu = (nxel_ins-nxsl_ins+1)*(nyel_ins-nysl_ins+1)*(nzel_ins-nzsl_ins+1)
        allocate(self%points_x(n_points_insitu))
        allocate(self%points_y(n_points_insitu))
        allocate(self%points_z(n_points_insitu))
      endif
    endassociate
  endsubroutine insitu_allocate

  subroutine insitu_finalize(self)
    class(equation_singleideal_object), intent(inout) :: self
    if(self%insitu_platform == "catalyst-v1") then
      call insitu_end(self%fcoproc)
    elseif(self%insitu_platform == "catalyst-v2") then
      call self%insitu_do_catalyst_finalization()
    endif
  endsubroutine insitu_finalize

  subroutine insitu_do_catalyst_finalization(self)
    class(equation_singleideal_object), intent(inout) :: self
    type(C_PTR) :: params
    integer(kind(catalyst_status)) :: code
    integer :: exit_code
    print*,'start do_catalyst_finalization'
    params = catalyst_conduit_node_create()
    code = c_catalyst_finalize(params)
    if (code /= catalyst_status_ok) then
      write (error_unit, *) "failed to call finalize:", code
      exit_code = 1
    end if
    call catalyst_conduit_node_destroy(params)
    print*,'end do_catalyst_finalization'
  endsubroutine insitu_do_catalyst_finalization

  function time_is_freezed_fun(self)
    class(equation_singleideal_object), intent(inout) :: self
    integer :: i
    logical :: time_is_freezed_fun
    time_is_freezed_fun = .false.
    if (self%enable_insitu > 0) then
      do i=1,size(self%freeze_intervals),2
        if (self%icyc == self%freeze_intervals(i)) then
          self%i_freeze = self%i_freeze + 1
          time_is_freezed_fun = .true.
          if (self%i_freeze == self%freeze_intervals(i+1) + 1) then
            self%i_freeze = 0
            time_is_freezed_fun = .false.
            return
          endif
          return
        endif
      enddo
    endif
  endfunction time_is_freezed_fun

  subroutine slice_prepare(self, igslice, jgslice, kgslice, islice, jslice, kslice, alloc_aux)
    class(equation_singleideal_object), intent(inout) :: self

    integer, dimension(:), intent(in) :: igslice, jgslice, kgslice
    integer, allocatable, dimension(:), intent(out):: islice, jslice, kslice
    logical, intent(in) :: alloc_aux
    integer :: i, j, k, l, ip, jp, kp, n
    integer :: ii, jj, kk, ll
    integer :: icord, jcord, kcord
    integer :: inum, jnum, knum

    associate(ng => self%grid%ng, nx => self%field%nx, ny => self%field%ny, nz => self%field%nz,&
    & x => self%field%x, y => self%field%y, z => self%field%z, num_aux_slice => self%num_aux_slice)

      inum = 0
      do i = 1,size(igslice)
        if (igslice(i)>0) then
          icord = (igslice(i)-1)/nx
          if (self%field%ncoords(1)==icord) inum = inum + 1
        endif
      enddo
      jnum = 0
      do j = 1,size(jgslice)
        if (jgslice(j)>0) then
          jcord = (jgslice(j)-1)/ny
          if (self%field%ncoords(2)==jcord) jnum = jnum + 1
        endif
      enddo
      knum = 0
      do k = 1,size(kgslice)
        if (kgslice(k)>0) then
          kcord = (kgslice(k)-1)/nz
          if (self%field%ncoords(3)==kcord) knum = knum + 1
        endif
      enddo
      if (inum>0) then
        allocate(islice(inum))
        if (alloc_aux) allocate(sliceyz_aux(inum,1-ng:ny+ng,1-ng:nz+ng,num_aux_slice))
      endif
      if (jnum>0) then
        allocate(jslice(jnum))
        if (alloc_aux) allocate(slicexz_aux(1-ng:nx+ng,jnum,1-ng:nz+ng,num_aux_slice))
      endif
      if (knum>0) then
        allocate(kslice(knum))
        if (alloc_aux) allocate(slicexy_aux(1-ng:nx+ng,1-ng:ny+ng,knum,num_aux_slice))
      endif
      inum = 0
      do i = 1,size(igslice)
        if (igslice(i)>0) then
          icord = (igslice(i)-1)/nx
          if (self%field%ncoords(1)==icord) then
            inum = inum+1
            islice(inum) = igslice(i)-self%field%ncoords(1)*nx
          endif
        endif
      enddo
      jnum = 0
      do j = 1,size(jgslice)
        if (jgslice(j)>0) then
          jcord = (jgslice(j)-1)/ny
          if (self%field%ncoords(2)==jcord) then
            jnum = jnum + 1
            jslice(jnum) = jgslice(j)-self%field%ncoords(2)*ny
          endif
        endif
      enddo
      knum = 0
      do k = 1,size(kgslice)
        if (kgslice(k)>0) then
          kcord = (kgslice(k)-1)/nz
          if (self%field%ncoords(3)==kcord) then
            knum = knum + 1
            kslice(knum) = kgslice(k)-self%field%ncoords(3)*nz
          endif
        endif
      enddo

    endassociate

  endsubroutine slice_prepare

  subroutine probe_prepare(self)
    class(equation_singleideal_object), intent(inout) :: self

    integer :: i, j, k, l, ip, jp, kp, n
    integer :: ii, jj, kk, ll
    integer :: inum, jnum, knum
    real(rkind) :: xp, yp, zp
    logical :: probe_exists
    logical :: in_i, in_j, in_k
    logical, dimension(:), allocatable :: in_ijk
    integer :: mp
    integer :: idp

    associate(ng => self%grid%ng, nx => self%field%nx, ny => self%field%ny, nz => self%field%nz,&
    & x => self%field%x, y => self%field%y, z => self%field%z)

      inquire(file='probe_list.dat',exist=probe_exists)
      self%num_probe = 0
      if (probe_exists) then
        open(18,file='probe_list.dat')
        read(18,*) self%num_probe_tot
        allocate(in_ijk(self%num_probe_tot))
        in_ijk = .false.
        do l=1,self%num_probe_tot
          read(18,*) idp, xp, yp, zp
          in_i = .false.
          in_j = .false.
          in_k = .false.
          if (xp>=self%field%x(1).and.xp<self%field%x(nx+1)) in_i = .true.
          if (yp>=self%field%y(1).and.yp<self%field%y(ny+1)) in_j = .true.
          if (zp>=self%field%z(1).and.zp<self%field%z(nz+1)) in_k = .true.
          if (in_i.and.in_j.and.in_k) then
            in_ijk(l) = .true.
            self%num_probe = self%num_probe+1
          endif
        enddo
        close(18)

        if (self%num_probe>0) then
          allocate(self%moving_probe(self%num_probe))
          allocate(self%id_probe(self%num_probe))
          allocate(self%probe_coord(3,self%num_probe))
          allocate(self%w_aux_probe(6,self%num_probe))
          allocate(self%ijk_probe(3,self%num_probe))
          allocate(self%probe_coeff(2,2,2,self%num_probe))
        endif

        open(18,file='probe_list.dat')
        read(18,*) self%num_probe_tot
        n = 0
        do l=1,self%num_probe_tot
          read(18,*) idp, xp, yp, zp, mp
          if (in_ijk(l)) then
            n = n+1
            self%id_probe(n) = idp
            self%probe_coord(1,n) = xp
            self%probe_coord(2,n) = yp
            self%probe_coord(3,n) = zp
            self%moving_probe(n) = mp
          endif
        enddo
        close(18)
      endif
    endassociate

  endsubroutine probe_prepare

  subroutine probe_compute_coeff(self)
    class(equation_singleideal_object), intent(inout) :: self

    integer :: i, j, k, l, ip, jp, kp, n
    integer :: ii, jj, kk, ll
    integer :: inum, jnum, knum
    real(rkind) :: xp, yp, zp
    real(rkind) :: x0, y0, z0, dxloc, dyloc, dzloc
    real(rkind) :: xyz1, xyz2, xyz3
    logical :: probe_exists
    real(rkind), dimension(8,8) :: amat3d
    real(rkind), dimension(1,8) :: xtrasp3d,alftrasp3d
    integer :: moving_probe

    associate(ng => self%grid%ng, nx => self%field%nx, ny => self%field%ny, nz => self%field%nz,&
    & x => self%field%x, y => self%field%y, z => self%field%z)

      do n=1,self%num_probe
        xp = self%probe_coord(1,n)
        yp = self%probe_coord(2,n)
        zp = self%probe_coord(3,n)
        call locateval(x(1:nx),nx,xp,ip)
        call locateval(y(1:ny),ny,yp,jp)
        call locateval(z(1:nz),nz,zp,kp)
        self%ijk_probe(1,n) = ip
        self%ijk_probe(2,n) = jp
        self%ijk_probe(3,n) = kp
        x0 = x(ip)
        y0 = y(jp)
        z0 = z(kp)
        dxloc = x(ip+1)-x(ip)
        dyloc = y(jp+1)-y(jp)
        dzloc = z(kp+1)-z(kp)
        xp = (xp-x0)/dxloc
        yp = (yp-y0)/dyloc
        zp = (zp-z0)/dzloc
        xtrasp3d(1,1) = xp*yp*zp
        xtrasp3d(1,2) = xp*yp
        xtrasp3d(1,3) = xp*zp
        xtrasp3d(1,4) = yp*zp
        xtrasp3d(1,5) = xp
        xtrasp3d(1,6) = yp
        xtrasp3d(1,7) = zp
        xtrasp3d(1,8) = 1._rkind
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
      enddo

    endassociate

  endsubroutine probe_compute_coeff

  subroutine recyc_prepare(self)
    class(equation_singleideal_object), intent(inout) :: self
    integer :: ig_recyc, j, k
    real(rkind) :: bexprecyc_base, my_eta0, my_deltablend
    real(rkind) :: alfa_blend, beta_blend, c_blend

    call self%cfg%get("bc","x_recyc",self%x_recyc)
    self%x_recyc = self%x_recyc*self%l0

    associate(ng => self%grid%ng, ny => self%field%ny, nz => self%field%nz)
      allocate(self%yplus_inflow(1-ng:ny+ng))
      allocate(self%eta_inflow(1-ng:ny+ng))
      allocate(self%yplus_recyc(1-ng:ny+ng))
      allocate(self%eta_recyc(1-ng:ny+ng))
      allocate(self%eta_recyc_blend(1-ng:ny+ng))
      allocate(self%map_j_inn(1:ny))
      allocate(self%map_j_out(1:ny))
      allocate(self%map_j_out_blend(1:ny))
      allocate(self%weta_inflow(1:ny))
      allocate(self%inflow_random_plane(1:ny,1:nz,3))
    endassociate

    associate(xg => self%grid%xg, nxmax => self%grid%nxmax, nx => self%field%nx, ny => self%field%ny&
    &, ng => self%grid%ng, xrecyc => self%x_recyc, i_recyc => self%i_recyc, ib_recyc => self%ib_recyc,&
    & deltavec => self%deltavec, deltavvec => self%deltavvec, betarecyc => self%betarecyc,&
    & y => self%field%y, yplus_inflow => self%yplus_inflow, eta_inflow => self%eta_inflow,&
    & yplus_recyc => self%yplus_recyc, eta_recyc => self%eta_recyc, l0 => self%l0, eta_recyc_blend => sel&
    &f%eta_recyc_blend, map_j_out_blend => self%map_j_out_blend, map_j_inn => self%map_j_inn,&
    & map_j_out => self%map_j_out, weta_inflow => self%weta_inflow, inflow_random_plane => self%inflow_ra&
    &ndom_plane, glund1 => self%glund1)

      inflow_random_plane = 0.5_rkind

      bexprecyc_base = 0.13_rkind
      my_eta0 = 0.08_rkind
      my_deltablend = 1.10_rkind

      call locateval(xg(1:nxmax),nxmax,xrecyc,ig_recyc)
      ib_recyc = (ig_recyc-1)/nx
      i_recyc = ig_recyc-nx*ib_recyc
      if (i_recyc<ng) then
        i_recyc = ng
        ig_recyc = i_recyc+nx*ib_recyc
      endif
      if (self%masterproc) write(*,*) 'Recycling station exactly at = ', xg(ig_recyc)

      betarecyc = deltavvec(ig_recyc)/deltavvec(1)
      glund1 = (deltavec(ig_recyc)/deltavec(1))**bexprecyc_base
      if (self%masterproc) write(*,*) 'Recycling beta factor = ', betarecyc
      if (self%masterproc) write(*,*) 'Pirozzoli-Ceci beta factor = ', glund1
      if (self%masterproc) then
        write(*,*) 'Urbin-Knight beta factor = ', ((1._rkind+(xrecyc/l0)*0.27_rkind**1.2_rkind/&
        &self%Reynolds**0.2_rkind)**(5._rkind/6._rkind))**0.1_rkind
      endif

      do j=1-ng,ny+ng
        yplus_inflow(j) = y(j)/deltavvec(1)
        eta_inflow(j) = y(j)/deltavec (1)
        yplus_recyc(j) = y(j)/deltavvec(ig_recyc)
        eta_recyc(j) = y(j)/deltavec (ig_recyc)
        eta_recyc_blend(j) = y(j)*((deltavec(1)/deltavec(ig_recyc))**bexprecyc_base+0.5_rkind*&
        &(1._rkind+tanh(log(abs(eta_recyc(j))/my_eta0)/my_deltablend))*((deltavec(1)/deltavec(ig_recyc))-&
        &(deltavec(1)/deltavec(ig_recyc))**bexprecyc_base))

      enddo
      alfa_blend = 4._rkind
      beta_blend = .2_rkind
      c_blend = 1._rkind-2._rkind*beta_blend
      do j=1,ny
        call locateval(yplus_recyc(1:ny),ny,yplus_inflow(j),map_j_inn(j))
        call locateval(eta_recyc(1:ny),ny,eta_inflow(j),map_j_out(j))
        call locateval(eta_recyc_blend(1:ny),ny,eta_inflow(j),map_j_out_blend(j))
        weta_inflow(j) = 0.5_rkind*(1._rkind+tanh((alfa_blend*(eta_inflow(j)-beta_blend))/&
        &(beta_blend+eta_inflow(j)*c_blend))/tanh(alfa_blend))
      enddo
    endassociate
  endsubroutine recyc_prepare

  subroutine alloc(self)
    class(equation_singleideal_object), intent(inout) :: self
    associate(nx => self%field%nx, ny => self%field%ny, nz => self%field%nz, ng => self%grid%ng,&
    & nv => self%nv, nv_aux => self%nv_aux, nv_stat => self%nv_stat, nv_stat_3d => self%nv_stat_3d,&
    & enable_stat_3d => self%enable_stat_3d, ndft => self%ndft_tspec)
      allocate(self%w_aux(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng, nv_aux))
      self%w_aux = 0
      allocate(self%w_stat(nv_stat, 1:nx, 1:ny))
      allocate(self%w_stat_z(nv_stat, 1:nx, 1:ny))
      allocate(self%fluid_mask(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng))
      allocate(self%ep_ord_change(0:nx, 0:ny, 0:nz, 1:3))
      allocate(self%wallprop(1-ng:nx+ng, 1-ng:nz+ng, 2:4))
      if (enable_stat_3d>0) allocate(self%w_stat_3d(nv_stat_3d, 1:nx, 1:ny, 1:nz))
      self%fluid_mask = 0
      self%ep_ord_change = 0

      allocate(self%lmax_tag(1-ng:nx+ng))
      allocate(self%vis_tag(1-ng:nx+ng))

      if (self%enable_tspec > 0) then
        allocate(self%w_tspec(nx, nz, ndft, 2, 2))
        allocate(self%w_psd_tspec(nx, ndft))
        allocate(self%w_psd_avg_tspec(nx, ndft))
      endif
    endassociate
  endsubroutine alloc

  subroutine fromaux_to_w(self)
    class(equation_singleideal_object), intent(inout) :: self
    integer :: i,j,k
    real(rkind) :: rho,uu,vv,ww,qq,tt,ee
    associate(nx => self%field%nx, ny => self%field%ny, nz => self%field%nz, t0 => self%t0,&
    & indx_cp_l => self%indx_cp_l, indx_cp_r => self%indx_cp_r, cv_coeff => self%cv_coeff,&
    & calorically_perfect => self%calorically_perfect, w => self%field%w)
      do k=1,nz
        do j=1,ny
          do i=1,nx
            rho = w(i,j,k,1)
            uu = w(i,j,k,2)
            vv = w(i,j,k,3)
            ww = w(i,j,k,4)
            tt = w(i,j,k,5)
            qq = 0.5_rkind*(uu*uu+vv*vv+ww*ww)
            ee = get_e_from_temperature(tt,t0,indx_cp_l,indx_cp_r,cv_coeff,calorically_perfect)
            w(i,j,k,1) = rho
            w(i,j,k,2) = rho*uu
            w(i,j,k,3) = rho*vv
            w(i,j,k,4) = rho*ww
            w(i,j,k,5) = rho*(ee+qq)
          enddo
        enddo
      enddo
    endassociate
  endsubroutine fromaux_to_w

  subroutine read_stats(self)
    class(equation_singleideal_object), intent(inout) :: self
    integer, dimension(3) :: sizes
    integer, dimension(3) :: subsizes
    integer, dimension(3) :: starts
    integer :: ntotxy
    integer, dimension(mpi_status_size) :: istatus
    integer :: mpi_io_file
    integer :: filetype, l, m
    integer :: size_real
    integer (kind=mpi_offset_kind) :: offset
    character(len=256) :: oldname, newname
    associate(nx => self%field%nx, ny => self%field%ny, nv_stat => self%nv_stat, w_stat => self%w_st&
    &at, mp_cartx => self%field%mp_cartx, mp_cartz => self%field%mp_cartz, ncoords => self%field%ncoords,&
    & nblocks => self%field%nblocks, iermpi => self%mpi_err)
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

        call mpi_file_close(mpi_io_file,iermpi)
        call mpi_type_free(filetype,iermpi)
      endif
      call mpi_bcast(w_stat,nv_stat*nx*ny,mpi_prec,0,mp_cartz,iermpi)
      call mpi_barrier(mpi_comm_world, iermpi)
      if (self%masterproc) then
        oldname = c_char_"stat.bin"//c_null_char
        newname = c_char_"stat.bak"//c_null_char
        iermpi = rename_wrapper(oldname, newname)
        if (iermpi /= 0) write(error_unit,*) "Warning! Cannot rename file stat.bin to stat.bak"
      endif
      call mpi_barrier(mpi_comm_world, iermpi)
    endassociate

  endsubroutine read_stats

  subroutine read_stats_serial(self)
    class(equation_singleideal_object), intent(inout) :: self
    character(4) :: chx,chy
    associate(nx => self%field%nx, ny => self%field%ny, nz => self%field%nz, nv_stat => self%nv_stat&
    &, w_stat => self%w_stat, ncoords => self%field%ncoords, mp_cartz => self%field%mp_cartz,&
    & iermpi => self%mpi_err)

      if (ncoords(3)==0) then

        if (self%masterproc) write(*,*) 'Reading stat0_XXX_XXX.bin'
        1004 format(I4.4)
        write(chx,1004) ncoords(1)
        write(chy,1004) ncoords(2)

        open (11,file='stat0_'//chx//'_'//chy//'.bin',form='unformatted')
        read(11) w_stat(1:nv_stat,1:nx,1:ny)
        close(11)

      endif
      call mpi_bcast(w_stat,nv_stat*nx*ny,mpi_prec,0,mp_cartz,iermpi)
    endassociate

  endsubroutine read_stats_serial

  subroutine read_stats_3d(self)
    class(equation_singleideal_object), intent(inout) :: self
    integer, dimension(3) :: sizes
    integer, dimension(3) :: subsizes
    integer, dimension(3) :: starts
    integer :: ntot3d
    integer, dimension(mpi_status_size) :: istatus
    integer :: mpi_io_file
    integer :: filetype, l, m
    integer :: size_real
    integer (kind=mpi_offset_kind) :: offset
    character(len=256) :: oldname, newname

    associate(nx => self%field%nx, ny => self%field%ny, nz => self%field%nz, nv_stat_3d => self%nv_s&
    &tat_3d, w_stat_3d => self%w_stat_3d, mp_cart => self%field%mp_cart, ncoords => self%field%ncoords,&
    & nblocks => self%field%nblocks, iermpi => self%mpi_err)

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

      call mpi_file_close(mpi_io_file,iermpi)
      call mpi_type_free(filetype,iermpi)
      if (self%masterproc) then
        oldname = c_char_"stat3d.bin"//c_null_char
        newname = c_char_"stat3d.bak"//c_null_char
        iermpi = rename_wrapper(oldname, newname)
        if (iermpi /= 0) write(error_unit,*) "Warning! Cannot rename file stat3d.bin to stat3d.bak"
      endif
      call mpi_barrier(mpi_comm_world, iermpi)
    endassociate

  endsubroutine read_stats_3d

  subroutine write_stats(self)
    class(equation_singleideal_object), intent(inout) :: self
    integer, dimension(3) :: sizes
    integer, dimension(3) :: subsizes
    integer, dimension(3) :: starts
    integer :: ntotxy
    integer, dimension(mpi_status_size) :: istatus
    integer :: mpi_io_file
    integer :: filetype, l, m
    integer :: size_real
    integer (kind=mpi_offset_kind) :: offset
    associate(nx => self%field%nx, ny => self%field%ny, nv_stat => self%nv_stat, w_stat => self%w_st&
    &at, mp_cartx => self%field%mp_cartx, ncoords => self%field%ncoords, nblocks => self%field%nblocks,&
    & iermpi => self%mpi_err)
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

        call mpi_file_close(mpi_io_file,iermpi)
        call mpi_type_free(filetype,iermpi)
      endif
    endassociate
  endsubroutine write_stats

  subroutine write_stats_z(self)
    class(equation_singleideal_object), intent(inout) :: self
    integer, dimension(3) :: sizes
    integer, dimension(3) :: subsizes
    integer, dimension(3) :: starts
    integer :: ntotxy
    integer, dimension(mpi_status_size) :: istatus
    integer :: mpi_io_file
    integer :: filetype, l, m
    integer :: size_real, istat
    integer (kind=mpi_offset_kind) :: offset
    character(4) :: chstore

    associate(nx => self%field%nx, ny => self%field%ny, nv_stat => self%nv_stat, w_stat => self%w_st&
    &at, mp_cartx => self%field%mp_cartx, w_stat_z => self%w_stat_z, ncoords => self%field%ncoords,&
    & nblocks => self%field%nblocks, iermpi => self%mpi_err)

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

        write(chstore(1:4), '(I4.4)') self%itav
        if (self%masterproc) write(*,*) 'saving stat z: ',self%itav, chstore
        call mpi_type_create_subarray(3,sizes,subsizes,starts,mpi_order_fortran,mpi_prec,filetype,iermpi)
        call mpi_type_commit(filetype,iermpi)
        call mpi_file_open(mp_cartx,'AVGZ/stat_z_'//chstore//'.bin',mpi_mode_create+mpi_mode_wronly,&
        &mpi_info_null, mpi_io_file,iermpi)
        offset = 0
        do l=1,nv_stat
          call mpi_file_set_view(mpi_io_file,offset,mpi_prec,filetype,"native",mpi_info_null,iermpi)
          call mpi_file_write_all(mpi_io_file,w_stat_z(l,1:nx,1:ny),ntotxy,mpi_prec,istatus,iermpi)
          call mpi_type_size(mpi_prec,size_real,iermpi)
          do m=1,nblocks(1)*nblocks(2)
            offset = offset+size_real*ntotxy
          enddo
        enddo
        call mpi_file_close(mpi_io_file,iermpi)
        call mpi_type_free(filetype,iermpi)
      endif
    endassociate
  endsubroutine write_stats_z

  subroutine write_stats_serial(self)
    class(equation_singleideal_object), intent(inout) :: self
    character(4) :: chx,chy
    associate(nx => self%field%nx, ny => self%field%ny, nz => self%field%nz, nv_stat => self%nv_stat&
    &, w_stat => self%w_stat, ncoords => self%field%ncoords, mp_cartz => self%field%mp_cartz)

      if (ncoords(3)==0) then

        if (self%masterproc) write(*,*) 'Writing stat1_XXX_XXX.bin'
        1004 format(I4.4)
        write(chx,1004) ncoords(1)
        write(chy,1004) ncoords(2)

        open (11,file='stat1_'//chx//'_'//chy//'.bin',form='unformatted')
        write(11) w_stat(1:nv_stat,1:nx,1:ny)
        close(11)

      endif
    endassociate

  endsubroutine write_stats_serial

  subroutine write_stats_z_serial(self)
    class(equation_singleideal_object), intent(inout) :: self
    character(4) :: chx,chy,chstore
    integer :: istat

    associate(nx => self%field%nx, ny => self%field%ny, nz => self%field%nz, nv_stat => self%nv_stat&
    &, w_stat => self%w_stat, w_stat_z => self%w_stat_z, ncoords => self%field%ncoords,&
    & mp_cartz => self%field%mp_cartz)

      if (ncoords(3)==0) then

        if (self%masterproc) write(*,*) 'Writing stat_z_YYYY_XXX_XXX.bin'
        1004 format(I4.4)
        write(chx,1004) ncoords(1)
        write(chy,1004) ncoords(2)

        if(self%save_stat_z) then
          write(chstore(1:4), '(I4.4)') self%itav
          open (11,file='AVGZ/stat_z_'//chstore//'_'//chx//'_'//chy//'.bin',form='unformatted')
          write(11) w_stat_z(1:nv_stat,1:nx,1:ny)
          close(11)
        endif
      endif
    endassociate

  endsubroutine write_stats_z_serial

  subroutine write_stats_3d(self)
    class(equation_singleideal_object), intent(inout) :: self
    integer, dimension(3) :: sizes
    integer, dimension(3) :: subsizes
    integer, dimension(3) :: starts
    integer :: ntot3d
    integer, dimension(mpi_status_size) :: istatus
    integer :: mpi_io_file
    integer :: filetype, l, m
    integer :: size_real
    integer (kind=mpi_offset_kind) :: offset
    associate(nx => self%field%nx, ny => self%field%ny, nz => self%field%nz, nv_stat_3d => self%nv_s&
    &tat_3d, w_stat_3d => self%w_stat_3d, mp_cart => self%field%mp_cart, ncoords => self%field%ncoords,&
    & nblocks => self%field%nblocks, iermpi => self%mpi_err)

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

      call mpi_file_close(mpi_io_file,iermpi)
      call mpi_type_free(filetype,iermpi)

    endassociate
  endsubroutine write_stats_3d
  subroutine write_stats_3d_serial(self)
    class(equation_singleideal_object), intent(inout) :: self
    character(4) :: chx,chy,chz
    associate(nx => self%field%nx, ny => self%field%ny, nz => self%field%nz, nv_stat_3d => self%nv_s&
    &tat_3d, w_stat_3d => self%w_stat_3d, ncoords => self%field%ncoords)

      if (self%masterproc) write(*,*) 'Writing stat3d1_XXX_XXX_XXX.bin'
      1004 format(I4.4)
      write(chx,1004) ncoords(1)
      write(chy,1004) ncoords(2)
      write(chz,1004) ncoords(3)

      open (11,file='stat3d1_'//chx//'_'//chy//'_'//chz//'.bin',form='unformatted')
      write(11) w_stat_3d(1:nv_stat_3d,1:nx,1:ny,1:nz)
      close(11)

    endassociate
  endsubroutine write_stats_3d_serial
  subroutine read_stats_3d_serial(self)
    class(equation_singleideal_object), intent(inout) :: self
    character(4) :: chx,chy,chz
    associate(nx => self%field%nx, ny => self%field%ny, nz => self%field%nz, nv_stat_3d => self%nv_s&
    &tat_3d, w_stat_3d => self%w_stat_3d, ncoords => self%field%ncoords)

      if (self%masterproc) write(*,*) 'Reading stat3d0_XXX_XXX_XXX.bin'
      1004 format(I4.4)
      write(chx,1004) ncoords(1)
      write(chy,1004) ncoords(2)
      write(chz,1004) ncoords(3)

      open (11,file='stat3d0_'//chx//'_'//chy//'_'//chz//'.bin',form='unformatted')
      read(11) w_stat_3d(1:nv_stat_3d,1:nx,1:ny,1:nz)
      close(11)

    endassociate
  endsubroutine read_stats_3d_serial
  subroutine compute_stats(self)
    class(equation_singleideal_object), intent(inout) :: self
    real(rkind) :: rho, rhou, rhov, rhow, rhoe, ri, uu, vv, ww, qq, pp, tt, rho2, uu2, vv2, ww2, pp2, tt2, mu, nu, uv
    real(rkind) :: ux,uy,uz,vx,vy,vz,wx,wy,wz,omx,omy,omz,divl,div3l,machlocal,gamloc,c,t_tot,ccl,cploc
    real(rkind), dimension(3,3) :: sig
    integer :: i,j,l,k,npt,lmax
    associate(nx => self%field%nx, ny => self%field%ny, nz => self%field%nz, nv_stat => self%nv_stat&
    &, ng => self%grid%ng, nv_aux => self%nv_aux, w_stat => self%w_stat, nzmax => self%grid%nzmax,&
    & w => self%field%w, visc_model => self%visc_model, w_stat_z => self%w_stat_z, powerlaw_vtexp => self&
    &%powerlaw_vtexp, mu0 => self%mu0, T_ref_dim => self%T_ref_dim, itav => self%itav,&
    & mp_cartz => self%field%mp_cartz, sutherland_S => self%sutherland_S, mpi_err => self%mpi_err,&
    & cv_coeff => self%cv_coeff, w_aux => self%w_aux, dcsidx => self%field%dcsidx, detady => self%field%d&
    &etady, dzitdz => self%field%dzitdz, indx_cp_l => self%indx_cp_l, indx_cp_r => self%indx_cp_r,&
    & cp_coeff => self%cp_coeff, calorically_perfect => self%calorically_perfect, rgas0 => self%rgas0,&
    & t0 => self%t0)

      lmax = self%visc_order/2

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
              ux = ux+ccl*(w(i+l,j,k,2)/w(i+l,j,k,1)-w(i-l,j,k,2)/w(i-l,j,k,1))
              vx = vx+ccl*(w(i+l,j,k,3)/w(i+l,j,k,1)-w(i-l,j,k,3)/w(i-l,j,k,1))
              wx = wx+ccl*(w(i+l,j,k,4)/w(i+l,j,k,1)-w(i-l,j,k,4)/w(i-l,j,k,1))
              uy = uy+ccl*(w(i,j+l,k,2)/w(i,j+l,k,1)-w(i,j-l,k,2)/w(i,j-l,k,1))
              vy = vy+ccl*(w(i,j+l,k,3)/w(i,j+l,k,1)-w(i,j-l,k,3)/w(i,j-l,k,1))
              wy = wy+ccl*(w(i,j+l,k,4)/w(i,j+l,k,1)-w(i,j-l,k,4)/w(i,j-l,k,1))
              uz = uz+ccl*(w(i,j,k+l,2)/w(i,j,k+l,1)-w(i,j,k-l,2)/w(i,j,k-l,1))
              vz = vz+ccl*(w(i,j,k+l,3)/w(i,j,k+l,1)-w(i,j,k-l,3)/w(i,j,k-l,1))
              wz = wz+ccl*(w(i,j,k+l,4)/w(i,j,k+l,1)-w(i,j,k-l,4)/w(i,j,k-l,1))
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
            rho = w(i,j,k,1)
            rhou = w(i,j,k,2)
            rhov = w(i,j,k,3)
            rhow = w(i,j,k,4)
            rhoe = w(i,j,k,5)
            ri = 1._rkind/rho
            uu = rhou*ri
            vv = rhov*ri
            ww = rhow*ri
            qq = 0.5_rkind*(uu*uu+vv*vv+ww*ww)


            tt = w_aux(i,j,k,6)
            pp = rho*rgas0*tt
            mu = w_aux(i,j,k,7)
            nu = mu/rho
            gamloc = get_gamloc(tt,t0,indx_cp_l,indx_cp_r,cp_coeff,calorically_perfect,rgas0)
            cploc = gamloc/(gamloc-1._rkind)*rgas0
            c = sqrt(gamloc*rgas0*tt)
            machlocal = sqrt(2._rkind*qq)/c
            t_tot = tt+qq*(gamloc-1._rkind)/gamloc
            rho2 = rho*rho
            uu2 = uu*uu
            vv2 = vv*vv
            ww2 = ww*ww
            uv = uu*vv
            pp2 = pp*pp
            tt2 = tt*tt
            divl = (ux+vy+wz)
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
            sig = sig*mu
            w_stat_z(1,i,j) = w_stat_z(1,i,j)+rho
            w_stat_z(2,i,j) = w_stat_z(2,i,j)+uu
            w_stat_z(3,i,j) = w_stat_z(3,i,j)+vv
            w_stat_z(4,i,j) = w_stat_z(4,i,j)+ww
            w_stat_z(5,i,j) = w_stat_z(5,i,j)+pp
            w_stat_z(6,i,j) = w_stat_z(6,i,j)+tt
            w_stat_z(7,i,j) = w_stat_z(7,i,j)+rho2
            w_stat_z(8,i,j) = w_stat_z(8,i,j)+uu2
            w_stat_z(9,i,j) = w_stat_z(9,i,j)+vv2
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
            w_stat_z(22,i,j) = w_stat_z(22,i,j)+omx**2
            w_stat_z(23,i,j) = w_stat_z(23,i,j)+omy**2
            w_stat_z(24,i,j) = w_stat_z(24,i,j)+omz**2
            w_stat_z(25,i,j) = w_stat_z(25,i,j)+rho*tt
            w_stat_z(26,i,j) = w_stat_z(26,i,j)+rho*tt2
            w_stat_z(27,i,j) = w_stat_z(27,i,j)+t_tot
            w_stat_z(28,i,j) = w_stat_z(28,i,j)+rho*t_tot
            w_stat_z(29,i,j) = w_stat_z(29,i,j)+t_tot**2
            w_stat_z(30,i,j) = w_stat_z(30,i,j)+rhou*tt
            w_stat_z(31,i,j) = w_stat_z(31,i,j)+rhov*tt
            w_stat_z(32,i,j) = w_stat_z(32,i,j)+rhow*tt
            w_stat_z(33,i,j) = w_stat_z(33,i,j)+machlocal
            w_stat_z(34,i,j) = w_stat_z(34,i,j)+machlocal**2
            w_stat_z(35,i,j) = w_stat_z(35,i,j)+rhou*uu2
            w_stat_z(36,i,j) = w_stat_z(36,i,j)+rhov*uu2
            w_stat_z(37,i,j) = w_stat_z(37,i,j)+rhou*vv2
            w_stat_z(38,i,j) = w_stat_z(38,i,j)+rhov*vv2
            w_stat_z(39,i,j) = w_stat_z(39,i,j)+rhou*ww2
            w_stat_z(40,i,j) = w_stat_z(40,i,j)+rhov*ww2
            w_stat_z(41,i,j) = w_stat_z(41,i,j)+pp*uu
            w_stat_z(42,i,j) = w_stat_z(42,i,j)+pp*vv
            w_stat_z(43,i,j) = w_stat_z(43,i,j)+sig(1,1)
            w_stat_z(44,i,j) = w_stat_z(44,i,j)+sig(1,2)
            w_stat_z(45,i,j) = w_stat_z(45,i,j)+sig(1,3)
            w_stat_z(46,i,j) = w_stat_z(46,i,j)+sig(2,2)
            w_stat_z(47,i,j) = w_stat_z(47,i,j)+sig(2,3)
            w_stat_z(48,i,j) = w_stat_z(48,i,j)+sig(3,3)
            w_stat_z(49,i,j) = w_stat_z(49,i,j)+sig(1,1)*uu
            w_stat_z(50,i,j) = w_stat_z(50,i,j)+sig(1,2)*uu
            w_stat_z(51,i,j) = w_stat_z(51,i,j)+sig(2,1)*vv
            w_stat_z(52,i,j) = w_stat_z(52,i,j)+sig(2,2)*vv
            w_stat_z(53,i,j) = w_stat_z(53,i,j)+sig(3,1)*ww
            w_stat_z(54,i,j) = w_stat_z(54,i,j)+sig(3,2)*ww
            w_stat_z(55,i,j) = w_stat_z(55,i,j)+sig(1,1)*vv+sig(2,1)*uu
            w_stat_z(56,i,j) = w_stat_z(56,i,j)+sig(1,2)*vv+sig(2,2)*uu
            w_stat_z(57,i,j) = w_stat_z(57,i,j)+sig(1,1)*ux+sig(1,2)*uy+sig(1,3)*uz
            w_stat_z(58,i,j) = w_stat_z(58,i,j)+sig(2,1)*vx+sig(2,2)*vy+sig(2,3)*vz
            w_stat_z(59,i,j) = w_stat_z(59,i,j)+sig(3,1)*wx+sig(3,2)*wy+sig(3,3)*wz
            w_stat_z(60,i,j) = w_stat_z(60,i,j)+sig(1,1)*vx+sig(1,2)*(ux+vy)+sig(2,2)*uy+sig(1,3)*vz+sig(2,3)*uz
            w_stat_z(61,i,j) = w_stat_z(61,i,j)+pp*ux
            w_stat_z(62,i,j) = w_stat_z(62,i,j)+pp*vy
            w_stat_z(63,i,j) = w_stat_z(63,i,j)+pp*wz
            w_stat_z(64,i,j) = w_stat_z(64,i,j)+pp*(uy+vx)
            w_stat_z(65,i,j) = w_stat_z(65,i,j)+divl*divl
            w_stat_z(66,i,j) = w_stat_z(66,i,j)+rho*tt2*tt
            w_stat_z(67,i,j) = w_stat_z(67,i,j)+rho*tt2*tt2
            w_stat_z(68,i,j) = w_stat_z(68,i,j)+rhou*uu2*uu
            w_stat_z(69,i,j) = w_stat_z(69,i,j)+cploc
            w_stat_z(70,i,j) = w_stat_z(70,i,j)+gamloc
          enddo
        enddo
      enddo

      npt = nv_stat*nx*ny
      call mpi_allreduce(MPI_IN_PLACE,w_stat_z,npt,mpi_prec,mpi_sum,mp_cartz,mpi_err)
      w_stat_z = w_stat_z/nzmax

      w_stat = w_stat*itav + w_stat_z
      w_stat = w_stat/(itav+1)
    endassociate
  endsubroutine compute_stats

  subroutine compute_airfoil_forces(self)
    class(equation_singleideal_object), intent(inout) :: self
    real(rkind) :: dudy, dudyw, tauw, pw, dudy1, dudyw1, tauw1, pw1, cf, cp, pwf, tauwf
    real(rkind) :: ds, costh, sinth, dN, dA, dpN, dpA, dtN, dtA, ut1, ut2, ut3, ut4
    character(4) :: chstore
    real(rkind) :: N, A, pN, pA, tN, tA
    real(rkind) :: Nglo, Aglo, pNglo, pAglo, tNglo, tAglo
    real(rkind), dimension(self%field%nx) :: xc2_vec,cf_vec,cp_vec,pw_vec,tauw_vec
    real(rkind), dimension(self%grid%nxmax) :: xc2_vecg,cf_vecg,cp_vecg,pw_vecg,tauw_vecg
    integer :: i, i_airfoil_start, i_airfoil_end
    real(rkind) :: al , pdyn , lift , drag , pres , fric

    if (self%field%ncoords(3)==0) then
      associate(nx => self%field%nx, ny => self%field%ny, nz => self%field%nz, nv_stat => self%nv_st&
      &at, ng => self%grid%ng, w_stat_z => self%w_stat_z, nxmax => self%grid%nxmax, ite_rank_x => self%fiel&
      &d%ite_rank_x, itu_rank_x => self%field%itu_rank_x, nrank => self%myrank, ite => self%grid%ite,&
      & itu => self%grid%itu, ite_l => self%field%ite_l, itu_l => self%field%itu_l, xc2 => self%field%xc2,&
      & meta => self%field%meta, csimod => self%field%csimod, dxdcsic2 => self%field%dxdcsic2,&
      & dydcsic2 => self%field%dydcsic2, iermpi => self%mpi_err, mp_cartx => self%field%mp_cartx,&
      & time => self%time, aoa => self%aoa, itav => self%itav, masterproc => self%masterproc,&
      & u0 => self%u0, p0 => self%p0, wall_tag => self%field%wall_tag)

        N = 0._rkind
        A = 0._rkind
        pN = 0._rkind
        pA = 0._rkind
        tN = 0._rkind
        tA = 0._rkind

        xc2_vec(:) = 0._rkind
        cf_vec(:) = 0._rkind
        cp_vec(:) = 0._rkind
        pw_vec(:) = 0._rkind
        tauw_vec(:) = 0._rkind

        do i=1,nx
          if(wall_tag(i) < 1) then
            ds = csimod(i,1)
            costh = dxdcsic2(i,1)/ds
            sinth = dydcsic2(i,1)/ds

            ut1 = (w_stat_z(2,i,1)*costh+w_stat_z(3,i,1)*sinth)/w_stat_z(1,i,1)
            ut2 = (w_stat_z(2,i,2)*costh+w_stat_z(3,i,2)*sinth)/w_stat_z(1,i,2)
            ut3 = (w_stat_z(2,i,3)*costh+w_stat_z(3,i,3)*sinth)/w_stat_z(1,i,3)
            ut4 = (w_stat_z(2,i,4)*costh+w_stat_z(3,i,4)*sinth)/w_stat_z(1,i,4)
            dudy = -22._rkind*ut1+36._rkind*ut2-18._rkind*ut3+4._rkind*ut4
            dudyw = dudy*meta(i,1)/12._rkind
            tauw = w_stat_z(20,i,1)*dudyw
            pw = w_stat_z(5,i,1)

            cf = tauw/(0.5_rkind*u0*u0)
            cp = (pw-p0)/(0.5_rkind*u0*u0)
            pwf = pw
            tauwf = tauw

            dN = -costh*pwf+sinth*tauwf
            dA = +sinth*pwf+costh*tauwf
            dpN = -costh*pwf
            dpA = +sinth*pwf
            dtN = sinth*tauwf
            dtA = costh*tauwf

            if(wall_tag(i) < -1) ds = ds/2._rkind

            N = N + dN*ds
            A = A + dA*ds
            pN = pN + dpN*ds
            pA = pA + dpA*ds
            tN = tN + dtN*ds
            tA = tA + dtA*ds

            xc2_vec(i) = xc2(i,1)
            cf_vec(i) = cf
            cp_vec(i) = cp
            pw_vec(i) = pw
            tauw_vec(i) = tauw
          endif
        enddo
        call mpi_gather(xc2_vec,nx,mpi_prec,xc2_vecg,nx,mpi_prec,0,mp_cartx,iermpi)
        call mpi_gather(cf_vec,nx,mpi_prec,cf_vecg,nx,mpi_prec,0,mp_cartx,iermpi)
        call mpi_gather(cp_vec,nx,mpi_prec,cp_vecg,nx,mpi_prec,0,mp_cartx,iermpi)
        call mpi_gather(pw_vec,nx,mpi_prec,pw_vecg,nx,mpi_prec,0,mp_cartx,iermpi)
        call mpi_gather(tauw_vec,nx,mpi_prec,tauw_vecg,nx,mpi_prec,0,mp_cartx,iermpi)
        if(masterproc) then
          write(chstore(1:4), '(I4.4)') itav
          open(unit=12,file='AVGZ/airfoil_coeffs_'//chstore//'.dat',form='formatted')
          do i=ite,itu
            write(12,100) xc2_vecg(i),cf_vecg(i),cp_vecg(i),pw_vecg(i),tauw_vecg(i)
          enddo
          close(12)
        endif
        Nglo = 0.
        Aglo = 0.
        pNglo = 0.
        pAglo = 0.
        tNglo = 0.
        tAglo = 0.

        call MPI_REDUCE(N, Nglo, 1, mpi_prec, MPI_SUM, 0, mp_cartx, iermpi)
        call MPI_REDUCE(A, Aglo, 1, mpi_prec, MPI_SUM, 0, mp_cartx, iermpi)
        call MPI_REDUCE(pN, pNglo, 1, mpi_prec, MPI_SUM, 0, mp_cartx, iermpi)
        call MPI_REDUCE(pA, pAglo, 1, mpi_prec, MPI_SUM, 0, mp_cartx, iermpi)
        call MPI_REDUCE(tN, tNglo, 1, mpi_prec, MPI_SUM, 0, mp_cartx, iermpi)
        call MPI_REDUCE(tA, tAglo, 1, mpi_prec, MPI_SUM, 0, mp_cartx, iermpi)

        if(masterproc) then
          al = aoa*pi/180._rkind
          pdyn = 0.5_rkind*u0*u0
          lift = ( Nglo*cos(al)- Aglo*sin(al))/pdyn
          drag = ( Nglo*sin(al)+ Aglo*cos(al))/pdyn
          pres = (pNglo*sin(al)+pAglo*cos(al))/pdyn
          fric = (tNglo*sin(al)+tAglo*cos(al))/pdyn
          write(*,*) 'Lift =',lift,' Drag =',drag
          write(*,*) 'Pressure =',pres,' Friction =',fric
          open(unit=30,file='airfoil_forces.dat',position='append')
          write(30,100) self%icyc,time,lift,drag,pres,fric,Nglo,Aglo,pNglo,pAglo,tNglo,tAglo
          close(30)
        endif
      endassociate
    endif
    100 format(I0,2X,200ES20.10)
  endsubroutine compute_airfoil_forces

  subroutine compute_stats_c2(self)
    class(equation_singleideal_object), intent(inout) :: self
    real(rkind) :: rho, rhou, rhov, rhow, rhoe, ri, uu, vv, ww, qq, pp, tt, rho2, uu2, vv2, ww2, pp2, tt2, mu, nu, uv
    real(rkind) :: ux,uy,uz,vx,vy,vz,wx,wy,wz,omx,omy,omz,divl,div3l,machlocal,gamloc,c,t_tot,ccl,cploc
    real(rkind), dimension(3,3) :: sig
    real(rkind) :: ucsi, vcsi, wcsi, ueta, veta, weta, uzit, vzit, wzit, ut, vt, rhout, rhovt,cli,clj,clk
    integer :: i,j,l,k,npt,lmax,lmaxi,lmaxj

    associate(nx => self%field%nx, ny => self%field%ny, nz => self%field%nz, nv_stat => self%nv_stat&
    &, ng => self%grid%ng, nv_aux => self%nv_aux, w_stat => self%w_stat, nzmax => self%grid%nzmax,&
    & w => self%field%w, visc_model => self%visc_model, w_stat_z => self%w_stat_z, powerlaw_vtexp => self&
    &%powerlaw_vtexp, mu0 => self%mu0, T_ref_dim => self%T_ref_dim, itav => self%itav,&
    & mp_cartz => self%field%mp_cartz, sutherland_S => self%sutherland_S, mpi_err => self%mpi_err,&
    & cv_coeff => self%cv_coeff, w_aux => self%w_aux, dcsidx => self%field%dcsidx, detady => self%field%d&
    &etady, dzitdz => self%field%dzitdz, indx_cp_l => self%indx_cp_l, indx_cp_r => self%indx_cp_r,&
    & cp_coeff => self%cp_coeff, calorically_perfect => self%calorically_perfect, rgas0 => self%rgas0,&
    & t0 => self%t0, wall_tag => self%field%wall_tag, vis_tag => self%vis_tag, dcsidxc2 => self%field%dcs&
    &idxc2, detadxc2 => self%field%detadxc2, dcsidyc2 => self%field%dcsidyc2, detadyc2 => self%field%deta&
    &dyc2, dxdcsinc2 => self%field%dxdcsinc2, dydcsinc2 => self%field%dydcsinc2, dxdetanc2 => self%field%&
    &dxdetanc2, dydetanc2 => self%field%dydetanc2, dcsidxnc2 => self%field%dcsidxnc2, dcsidync2 => self%f&
    &ield%dcsidync2)

      lmax = self%visc_order/2
      w_stat_z = 0._rkind

      do k=1,nz
        do j=1,ny
          do i=1,nx
            lmaxi = lmax
            lmaxj = lmax
            if (j==1) lmaxi = vis_tag(i)
            if (wall_tag(i)<1) lmaxj = min(j,lmax)

            ucsi = 0._rkind
            vcsi = 0._rkind
            wcsi = 0._rkind
            ueta = 0._rkind
            veta = 0._rkind
            weta = 0._rkind
            uzit = 0._rkind
            vzit = 0._rkind
            wzit = 0._rkind
            do l=1,lmax
              cli = self%coeff_deriv1(l,lmaxi)
              clj = self%coeff_deriv1(l,lmaxj)
              clk = self%coeff_deriv1(l,lmax )
              ucsi = ucsi+cli*(w(i+l,j,k,2)/w(i+l,j,k,1)-w(i-l,j,k,2)/w(i-l,j,k,1))
              vcsi = vcsi+cli*(w(i+l,j,k,3)/w(i+l,j,k,1)-w(i-l,j,k,3)/w(i-l,j,k,1))
              wcsi = wcsi+cli*(w(i+l,j,k,4)/w(i+l,j,k,1)-w(i-l,j,k,4)/w(i-l,j,k,1))
              ueta = ueta+clj*(w(i,j+l,k,2)/w(i,j+l,k,1)-w(i,j-l,k,2)/w(i,j-l,k,1))
              veta = veta+clj*(w(i,j+l,k,3)/w(i,j+l,k,1)-w(i,j-l,k,3)/w(i,j-l,k,1))
              weta = weta+clj*(w(i,j+l,k,4)/w(i,j+l,k,1)-w(i,j-l,k,4)/w(i,j-l,k,1))
              uzit = uzit+clk*(w(i,j,k+l,2)/w(i,j,k+l,1)-w(i,j,k-l,2)/w(i,j,k-l,1))
              vzit = vzit+clk*(w(i,j,k+l,3)/w(i,j,k+l,1)-w(i,j,k-l,3)/w(i,j,k-l,1))
              wzit = wzit+clk*(w(i,j,k+l,4)/w(i,j,k+l,1)-w(i,j,k-l,4)/w(i,j,k-l,1))
            enddo
            ux = ucsi*dcsidxc2(i,j) + ueta*detadxc2(i,j)
            vx = vcsi*dcsidxc2(i,j) + veta*detadxc2(i,j)
            wx = wcsi*dcsidxc2(i,j) + weta*detadxc2(i,j)
            uy = ucsi*dcsidyc2(i,j) + ueta*detadyc2(i,j)
            vy = vcsi*dcsidyc2(i,j) + veta*detadyc2(i,j)
            wy = wcsi*dcsidyc2(i,j) + weta*detadyc2(i,j)
            uz = dzitdz(k)*uzit
            vz = dzitdz(k)*vzit
            wz = dzitdz(k)*wzit

            omx = wy-vz
            omy = uz-wx
            omz = vx-uy

            rho = w(i,j,k,1)


            if(self%flow_init == 4) then
              rhout = w(i,j,k,2)*dxdcsinc2(i,j)+w(i,j,k,3)*dydcsinc2(i,j)
              rhovt = w(i,j,k,2)*dxdetanc2(i,j)+w(i,j,k,3)*dydetanc2(i,j)
              rhou = rhout
              rhov = rhovt
            else
              rhou = w(i,j,k,2)
              rhov = w(i,j,k,3)
            endif



            rhow = w(i,j,k,4)
            rhoe = w(i,j,k,5)

            ri = 1._rkind/rho
            uu = rhou*ri
            vv = rhov*ri
            ww = rhow*ri
            qq = 0.5_rkind*(uu*uu+vv*vv+ww*ww)


            tt = w_aux(i,j,k,6)
            pp = rho*rgas0*tt
            mu = w_aux(i,j,k,7)
            nu = mu/rho

            gamloc = get_gamloc(tt,t0,indx_cp_l,indx_cp_r,cp_coeff,calorically_perfect,rgas0)
            cploc = gamloc/(gamloc-1._rkind)*rgas0
            c = sqrt(gamloc*rgas0*tt)
            machlocal = sqrt(2._rkind*qq)/c
            t_tot = tt+qq*(gamloc-1._rkind)/gamloc
            rho2 = rho*rho
            uu2 = uu*uu
            vv2 = vv*vv
            ww2 = ww*ww
            uv = uu*vv
            pp2 = pp*pp
            tt2 = tt*tt
            divl = (ux+vy+wz)
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
            sig = sig*mu
            w_stat_z(1,i,j) = w_stat_z(1,i,j)+rho
            w_stat_z(2,i,j) = w_stat_z(2,i,j)+uu
            w_stat_z(3,i,j) = w_stat_z(3,i,j)+vv
            w_stat_z(4,i,j) = w_stat_z(4,i,j)+ww
            w_stat_z(5,i,j) = w_stat_z(5,i,j)+pp
            w_stat_z(6,i,j) = w_stat_z(6,i,j)+tt
            w_stat_z(7,i,j) = w_stat_z(7,i,j)+rho2
            w_stat_z(8,i,j) = w_stat_z(8,i,j)+uu2
            w_stat_z(9,i,j) = w_stat_z(9,i,j)+vv2
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
            w_stat_z(22,i,j) = w_stat_z(22,i,j)+omx**2
            w_stat_z(23,i,j) = w_stat_z(23,i,j)+omy**2
            w_stat_z(24,i,j) = w_stat_z(24,i,j)+omz**2
            w_stat_z(25,i,j) = w_stat_z(25,i,j)+rho*tt
            w_stat_z(26,i,j) = w_stat_z(26,i,j)+rho*tt2
            w_stat_z(27,i,j) = w_stat_z(27,i,j)+t_tot
            w_stat_z(28,i,j) = w_stat_z(28,i,j)+rho*t_tot
            w_stat_z(29,i,j) = w_stat_z(29,i,j)+t_tot**2
            w_stat_z(30,i,j) = w_stat_z(30,i,j)+rhou*tt
            w_stat_z(31,i,j) = w_stat_z(31,i,j)+rhov*tt
            w_stat_z(32,i,j) = w_stat_z(32,i,j)+rhow*tt
            w_stat_z(33,i,j) = w_stat_z(33,i,j)+machlocal
            w_stat_z(34,i,j) = w_stat_z(34,i,j)+machlocal**2
            w_stat_z(35,i,j) = w_stat_z(35,i,j)+rhou*uu2
            w_stat_z(36,i,j) = w_stat_z(36,i,j)+rhov*uu2
            w_stat_z(37,i,j) = w_stat_z(37,i,j)+rhou*vv2
            w_stat_z(38,i,j) = w_stat_z(38,i,j)+rhov*vv2
            w_stat_z(39,i,j) = w_stat_z(39,i,j)+rhou*ww2
            w_stat_z(40,i,j) = w_stat_z(40,i,j)+rhov*ww2
            w_stat_z(41,i,j) = w_stat_z(41,i,j)+pp*uu
            w_stat_z(42,i,j) = w_stat_z(42,i,j)+pp*vv
            w_stat_z(43,i,j) = w_stat_z(43,i,j)+sig(1,1)
            w_stat_z(44,i,j) = w_stat_z(44,i,j)+sig(1,2)
            w_stat_z(45,i,j) = w_stat_z(45,i,j)+sig(1,3)
            w_stat_z(46,i,j) = w_stat_z(46,i,j)+sig(2,2)
            w_stat_z(47,i,j) = w_stat_z(47,i,j)+sig(2,3)
            w_stat_z(48,i,j) = w_stat_z(48,i,j)+sig(3,3)
            w_stat_z(49,i,j) = w_stat_z(49,i,j)+sig(1,1)*uu
            w_stat_z(50,i,j) = w_stat_z(50,i,j)+sig(1,2)*uu
            w_stat_z(51,i,j) = w_stat_z(51,i,j)+sig(2,1)*vv
            w_stat_z(52,i,j) = w_stat_z(52,i,j)+sig(2,2)*vv
            w_stat_z(53,i,j) = w_stat_z(53,i,j)+sig(3,1)*ww
            w_stat_z(54,i,j) = w_stat_z(54,i,j)+sig(3,2)*ww
            w_stat_z(55,i,j) = w_stat_z(55,i,j)+sig(1,1)*vv+sig(2,1)*uu
            w_stat_z(56,i,j) = w_stat_z(56,i,j)+sig(1,2)*vv+sig(2,2)*uu
            w_stat_z(57,i,j) = w_stat_z(57,i,j)+sig(1,1)*ux+sig(1,2)*uy+sig(1,3)*uz
            w_stat_z(58,i,j) = w_stat_z(58,i,j)+sig(2,1)*vx+sig(2,2)*vy+sig(2,3)*vz
            w_stat_z(59,i,j) = w_stat_z(59,i,j)+sig(3,1)*wx+sig(3,2)*wy+sig(3,3)*wz
            w_stat_z(60,i,j) = w_stat_z(60,i,j)+sig(1,1)*vx+sig(1,2)*(ux+vy)+sig(2,2)*uy+sig(1,3)*vz+sig(2,3)*uz
            w_stat_z(61,i,j) = w_stat_z(61,i,j)+pp*ux
            w_stat_z(62,i,j) = w_stat_z(62,i,j)+pp*vy
            w_stat_z(63,i,j) = w_stat_z(63,i,j)+pp*wz
            w_stat_z(64,i,j) = w_stat_z(64,i,j)+pp*(uy+vx)
            w_stat_z(65,i,j) = w_stat_z(65,i,j)+divl*divl
            w_stat_z(66,i,j) = w_stat_z(66,i,j)+rho*tt2*tt
            w_stat_z(67,i,j) = w_stat_z(67,i,j)+rho*tt2*tt2
            w_stat_z(68,i,j) = w_stat_z(68,i,j)+rhou*uu2*uu
            w_stat_z(69,i,j) = w_stat_z(69,i,j)+cploc
            w_stat_z(70,i,j) = w_stat_z(70,i,j)+gamloc
          enddo
        enddo
      enddo

      npt = nv_stat*nx*ny
      call mpi_allreduce(MPI_IN_PLACE,w_stat_z,npt,mpi_prec,mpi_sum,mp_cartz,mpi_err)
      w_stat_z = w_stat_z/nzmax

      w_stat = w_stat*itav + w_stat_z
      w_stat = w_stat/(itav+1)
    endassociate
  endsubroutine compute_stats_c2

  subroutine compute_stats_3d(self)
    class(equation_singleideal_object), intent(inout) :: self
    integer :: i,j,k
    real(rkind) :: pp
    associate(nx => self%field%nx, ny => self%field%ny, nz => self%field%nz, itav => self%itav,&
    & w_stat_3d => self%w_stat_3d, w => self%field%w, w_aux => self%w_aux)

      do k=1,nz
        do j=1,ny
          do i=1,nx
            w_stat_3d(1,i,j,k) = w_stat_3d(1,i,j,k) *itav + w(i,j,k,1)
            w_stat_3d(2,i,j,k) = w_stat_3d(2,i,j,k) *itav + w(i,j,k,2)
            w_stat_3d(3,i,j,k) = w_stat_3d(3,i,j,k) *itav + w(i,j,k,3)
            w_stat_3d(4,i,j,k) = w_stat_3d(4,i,j,k) *itav + w(i,j,k,4)
            w_stat_3d(5,i,j,k) = w_stat_3d(5,i,j,k) *itav + w(i,j,k,1)*w_aux(i,j,k,6)
            w_stat_3d(6,i,j,k) = w_stat_3d(6,i,j,k) *itav + w(i,j,k,1)**2
            w_stat_3d(7,i,j,k) = w_stat_3d(7,i,j,k) *itav + w(i,j,k,2)**2/w(i,j,k,1)
            w_stat_3d(8,i,j,k) = w_stat_3d(8,i,j,k) *itav + w(i,j,k,3)**2/w(i,j,k,1)
            w_stat_3d(9,i,j,k) = w_stat_3d(9,i,j,k) *itav + w(i,j,k,4)**2/w(i,j,k,1)
            w_stat_3d(10,i,j,k) = w_stat_3d(10,i,j,k)*itav + w(i,j,k,2)*w(i,j,k,3)/w(i,j,k,1)
            w_stat_3d(11,i,j,k) = w_stat_3d(11,i,j,k)*itav + w(i,j,k,2)*w(i,j,k,4)/w(i,j,k,1)
            w_stat_3d(12,i,j,k) = w_stat_3d(12,i,j,k)*itav + w(i,j,k,3)*w(i,j,k,4)/w(i,j,k,1)
            w_stat_3d(13,i,j,k) = w_stat_3d(13,i,j,k)*itav + w(i,j,k,1)*w_aux(i,j,k,6)**2
            w_stat_3d(14,i,j,k) = w_stat_3d(14,i,j,k)*itav + (w(i,j,k,1)*w_aux(i,j,k,6))**2
          enddo
        enddo
      enddo
      w_stat_3d = w_stat_3d/(itav+1)

    end associate
  end subroutine compute_stats_3d


  subroutine read_field_info(self)
    class(equation_singleideal_object), intent(inout) :: self
    character(len=256) :: oldname, newname

    if (self%io_type_r==1) self%field_info_cfg = parse_cfg("field_info0.dat")
    if (self%io_type_r==2) self%field_info_cfg = parse_cfg("field_info0.dat")
    call self%field_info_cfg%get("field_info","icyc0", self%icyc0)
    call self%field_info_cfg%get("field_info","time0", self%time0)
    call self%field_info_cfg%get("field_info","itav", self%itav)
    call self%field_info_cfg%get("field_info","itslice_p3d", self%itslice_p3d)
    call self%field_info_cfg%get("field_info","time_from_last_rst", self%time_from_last_rst)
    call self%field_info_cfg%get("field_info","time_from_last_write", self%time_from_last_write)
    call self%field_info_cfg%get("field_info","time_from_last_stat", self%time_from_last_stat)
    call self%field_info_cfg%get("field_info","time_from_last_slice", self%time_from_last_slice)
    call self%field_info_cfg%get("field_info","time_from_last_probe", self%time_from_last_probe)
    if (self%field_info_cfg%has_key("field_info","time_from_last_slice_p3d")) then
      call self%field_info_cfg%get("field_info","time_from_last_slice_p3d", self%time_from_last_slice_p3d)
    else
      self%time_from_last_slice_p3d = self%time_from_last_slice
    endif
    call self%field_info_cfg%get("field_info","istore", self%istore)
    if (self%enable_insitu > 0) then
      if (self%field_info_cfg%has_key("field_info","time_from_last_insitu")) then
        call self%field_info_cfg%get("field_info","time_from_last_insitu",self%time_from_last_insitu)
      else
        self%time_from_last_insitu = 0._rkind
      endif
      if (self%field_info_cfg%has_key("field_info","i_insitu")) then
        call self%field_info_cfg%get("field_info","i_insitu",self%i_insitu)
      else
        self%i_insitu = 1
      endif
      if (self%field_info_cfg%has_key("field_info","time_insitu")) then
        call self%field_info_cfg%get("field_info","time_insitu",self%time_insitu)
      else
        self%time_insitu = 0._rkind
      endif
    endif
    if (self%flow_init == 5) then
      call self%field_info_cfg%get("field_info","is",self%is_old)
      call self%field_info_cfg%get("field_info","lamz",self%lamz_old)
      call self%field_info_cfg%get("field_info","lamz1",self%lamz1_old)
      call self%field_info_cfg%get("field_info","lams",self%lams_old)
      call self%field_info_cfg%get("field_info","lams1",self%lams1_old)
      call self%field_info_cfg%get("field_info","phiz",self%phiz_old)
      call self%field_info_cfg%get("field_info","phiz1",self%phiz1_old)
      call self%field_info_cfg%get("field_info","phis",self%phis_old)
      call self%field_info_cfg%get("field_info","phis1",self%phis1_old)
    endif
    if (self%enable_tspec > 0) then
      if (self%field_info_cfg%has_key("field_info","it_win1_tspec")) then
        call self%field_info_cfg%get("field_info","it_win1_tspec",self%it_win1_tspec)
        call self%field_info_cfg%get("field_info","it_win2_tspec",self%it_win2_tspec)
        call self%field_info_cfg%get("field_info","it_nwin_tspec",self%it_nwin_tspec)
        call self%field_info_cfg%get("field_info","time_from_last_tspec", self%time_from_last_tspec)
      else
        self%it_win1_tspec = 1
        self%it_win2_tspec = 1-self%ndft_tspec/2
        self%it_nwin_tspec = 0
        self%time_from_last_tspec = 0._rkind
      endif
    endif
    call mpi_barrier(mpi_comm_world, self%mpi_err)
    if (self%io_type_r==2) then
      if (self%masterproc) then
        oldname = c_char_"field_info0.dat"//c_null_char
        newname = c_char_"field_info0.bak"//c_null_char
        self%mpi_err = rename_wrapper(oldname, newname)
        if (self%mpi_err /= 0) write(error_unit,*) "Warning! Cannot rename file field_info.dat to field_info.bak"
      endif
      call mpi_barrier(mpi_comm_world, self%mpi_err)
    endif
  endsubroutine read_field_info

  subroutine write_field_info(self)
    class(equation_singleideal_object), intent(inout) :: self
    if(self%masterproc) then
      call self%field_info_cfg%set("field_info","icyc0", self%icyc)
      call self%field_info_cfg%set("field_info","time0", self%time)
      call self%field_info_cfg%set("field_info","itav", self%itav)
      call self%field_info_cfg%set("field_info","itslice_p3d", self%itslice_p3d)
      call self%field_info_cfg%set("field_info","time_from_last_rst", self%time_from_last_rst)
      call self%field_info_cfg%set("field_info","time_from_last_write", self%time_from_last_write)
      call self%field_info_cfg%set("field_info","time_from_last_stat", self%time_from_last_stat)
      call self%field_info_cfg%set("field_info","time_from_last_slice", self%time_from_last_slice)
      call self%field_info_cfg%set("field_info","time_from_last_probe", self%time_from_last_probe)
      call self%field_info_cfg%set("field_info","time_from_last_slice_p3d", self%time_from_last_slice_p3d)
      call self%field_info_cfg%set("field_info","istore", self%istore)
      if (self%enable_insitu > 0) then
        call self%field_info_cfg%set("field_info","time_from_last_insitu",self%time_from_last_insitu)
        call self%field_info_cfg%set("field_info","i_insitu",self%i_insitu)
        call self%field_info_cfg%set("field_info","time_insitu",self%time_insitu)
      endif
      if (self%flow_init==5) then
        call self%field_info_cfg%set("field_info","is",self%is)
        call self%field_info_cfg%set("field_info","lamz",self%lamz)
        call self%field_info_cfg%set("field_info","lamz1",self%lamz1)
        call self%field_info_cfg%set("field_info","lams",self%lams)
        call self%field_info_cfg%set("field_info","lams1",self%lams1)
        call self%field_info_cfg%set("field_info","phiz",self%phiz)
        call self%field_info_cfg%set("field_info","phiz1",self%phiz1)
        call self%field_info_cfg%set("field_info","phis",self%phis)
        call self%field_info_cfg%set("field_info","phis1",self%phis1)
      endif
      if (self%enable_tspec > 0) then
        call self%field_info_cfg%set("field_info","it_win1_tspec",self%it_win1_tspec)
        call self%field_info_cfg%set("field_info","it_win2_tspec",self%it_win2_tspec)
        call self%field_info_cfg%set("field_info","it_nwin_tspec",self%it_nwin_tspec)
        call self%field_info_cfg%set("field_info","time_from_last_tspec", self%time_from_last_tspec)
      endif
      if (self%io_type_w==1) call self%field_info_cfg%write("field_info1.dat")
      if (self%io_type_w==2) call self%field_info_cfg%write("field_info0.dat")
    endif
  endsubroutine write_field_info

  function get_e_from_temperature(tt, t0, indx_cp_l, indx_cp_r, cv_coeff, calorically_perfect)
    real(rkind), intent(in) :: tt, t0
    integer, intent(in) :: indx_cp_l, indx_cp_r,calorically_perfect
    real(rkind), dimension(indx_cp_l:indx_cp_r+1), intent(in) :: cv_coeff
    real(rkind) :: get_e_from_temperature
    real(rkind) :: ee
    integer :: l

    ee = cv_coeff(indx_cp_r+1)
    if (calorically_perfect==1) then
      ee = ee+cv_coeff(0)*(tt/t0-1._rkind)
    else
      do l=indx_cp_l, indx_cp_r
        if (l==-1) then
          ee = ee+cv_coeff(l)*log(tt/t0)
        else
          ee = ee+cv_coeff(l)/(l+1._rkind)*((tt/t0)**(l+1)-1._rkind)
        endif
      enddo
    endif
    ee = ee*t0
    get_e_from_temperature = ee
  endfunction get_e_from_temperature

  subroutine set_fluid_prop(self)
    class(equation_singleideal_object), intent(inout) :: self
    integer :: l, num_coeff_cp
    real(rkind) :: cp_tref, cv_tref, hstar_over_rgas_tref_dim, bcoeff
    real(rkind), allocatable, dimension(:) :: cp_temp
    associate(Prandtl => self%Prandtl, gam => self%gam, visc_model => self%visc_model,&
    & sutherland_S => self%sutherland_S, calorically_perfect => self%calorically_perfect,&
    & powerlaw_vtexp => self%powerlaw_vtexp, indx_cp_l => self%indx_cp_l, indx_cp_r => self%indx_cp_r,&
    & T_ref_dim => self%T_ref_dim, gm => self%gm, gm1 => self%gm1, rfac => self%rfac, rgas0 => self%rgas0&
    &, cfg => self%cfg,t0 => self%t0, rfac_type => self%rfac_type)

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
      if (self%cfg%has_key("fluid","rfac_type")) then
        call self%cfg%get("fluid","rfac_type",rfac_type)
      else
        rfac_type = 0
      endif

      if(rfac_type == 0) then
        rfac = Prandtl**(1._rkind/3._rkind)
      else
        rfac = Prandtl**(1._rkind/2._rkind)
      endif
      allocate(self%cp_coeff(indx_cp_l:indx_cp_r+1))
      allocate(self%cv_coeff(indx_cp_l:indx_cp_r+1))
      self%cp_coeff = 0._rkind
      self%cv_coeff = 0._rkind
      if (calorically_perfect==0) then
        call cfg%get("fluid","cp_coeff",cp_temp)
        num_coeff_cp = size(cp_temp)
        if (num_coeff_cp==(indx_cp_r-indx_cp_l+2)) then
          self%cp_coeff(indx_cp_l:indx_cp_r+1) = cp_temp(1:num_coeff_cp)
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
        self%cv_coeff(0) = self%cp_coeff(0) - 1._rkind
      endif
      gm1 = gam-1._rkind
      gm = 1._rkind/gm1
      if (self%masterproc) write(*,*) 'Gamma: ', gam
      if (calorically_perfect==1) then
        self%cv_coeff(0) = gm
        self%cp_coeff(0) = gm*gam
        bcoeff = self%cv_coeff(0)+1._rkind
      endif
      self%cp_coeff(indx_cp_r+1) = bcoeff
      self%cv_coeff(indx_cp_r+1) = (bcoeff-1._rkind)

      self%cv_coeff = self%cv_coeff*rgas0
      self%cp_coeff = self%cp_coeff*rgas0
    endassociate

  endsubroutine set_fluid_prop

  subroutine set_flow_params(self)
    class(equation_singleideal_object), intent(inout) :: self
    real(rkind) :: gm1h
    real(rkind) :: c0, e0, etot0, k0, rtemp
    allocate(self%winf(self%nv))
    allocate(self%winf_past_shock(self%nv))
    associate(Mach => self%Mach, Reynolds => self%Reynolds, theta_wall => self%theta_wall,&
    & t0 => self%t0, rho0 => self%rho0, l0 => self%l0, p0 => self%p0, u0 => self%u0, mu0 => self%mu0,&
    & s0 => self%s0, ptot0 => self%ptot0, ttot0 => self%ttot0, powerlaw_vtexp => self%powerlaw_vtexp,&
    & T_ref_dim => self%T_ref_dim, rgas0 => self%rgas0, rfac => self%rfac, gm1 => self%gm1,&
    & gam => self%gam, T_recovery => self%T_recovery, T_wall => self%T_wall, winf => self%winf,&
    & flow_params_cfg => self%flow_params_cfg, indx_cp_l => self%indx_cp_l, indx_cp_r => self%indx_cp_r,&
    & Prandtl => self%Prandtl, cv_coeff => self%cv_coeff, calorically_perfect => self%calorically_perfect&
    &, aoa => self%aoa, cfg => self%cfg, winf_past_shock => self%winf_past_shock, cp_coeff => self%cp_coe&
    &ff)
      call cfg%get("flow","Reynolds",Reynolds)
      call cfg%get("flow","Mach",Mach)
      call cfg%get("flow","theta_wall",theta_wall)
      if (self%bctags(4)==7) call self%set_oblique_shock()

      gm1h = 0.5_rkind*gm1
      T_recovery = t0 * (1._rkind+gm1h*rfac*Mach**2)
      T_wall = theta_wall * (T_recovery-t0) + t0

      select case(self%flow_init)
      case(0,4)
        call self%set_chan_prop()
      case(1)
        call self%set_bl_prop()
      end select

      p0 = rho0*rgas0*t0
      c0 = sqrt(gam*rgas0*t0)
      u0 = Mach*c0
      e0 = get_e_from_temperature(t0, t0, indx_cp_l, indx_cp_r, cv_coeff, calorically_perfect)
      etot0 = e0+0.5_rkind*u0**2
      s0 = log(p0/rho0**gam)
      ttot0 = t0 * (1._rkind+gm1h*Mach**2)
      ptot0 = p0 * (1._rkind+gm1h*Mach**2)**(gam/gm1)

      winf(1) = rho0
      winf(2) = rho0*u0*cos(aoa*pi/180._rkind)
      winf(3) = rho0*u0*sin(aoa*pi/180._rkind)
      winf(4) = 0._rkind
      winf(5) = rho0*etot0
      if(self%masterproc) write(*,*) 'Angle of attack:', aoa
      if (abs(u0)<tol_iter) then
        mu0 = c0*rho0*l0/Reynolds
        u0 = c0
        if (self%masterproc) write(*,*) 'Reynolds number based on speed of sound'
      else
        mu0 = u0*rho0*l0/Reynolds
      endif
      k0 = mu0*cp_coeff(0)/Prandtl
      if(self%masterproc) then
        call flow_params_cfg%set("flow_params","u0", u0)
        call flow_params_cfg%set("flow_params","p0", p0)
        call flow_params_cfg%set("flow_params","rho0", rho0)
        call flow_params_cfg%set("flow_params","t0", t0)
        call flow_params_cfg%set("flow_params","c0", c0)
        call flow_params_cfg%set("flow_params","l0", l0)
        call flow_params_cfg%set("flow_params","cp0", cp_coeff(0))
        call flow_params_cfg%set("flow_params","cv0", cv_coeff(0))
        call flow_params_cfg%set("flow_params","mu0", mu0)
        call flow_params_cfg%set("flow_params","k0", k0)
        call flow_params_cfg%set("flow_params","gam", gam)
        call flow_params_cfg%set("flow_params","T_wall", T_wall)
        call flow_params_cfg%set("flow_params","theta_wall", theta_wall)
        call flow_params_cfg%set("flow_params","rfac", rfac)
        call flow_params_cfg%set("flow_params","powerlaw_vtexp", powerlaw_vtexp)
        call flow_params_cfg%set("flow_params","T_ref_dim", T_ref_dim)
        call flow_params_cfg%set("flow_params","Reynolds", Reynolds)
        call flow_params_cfg%set("flow_params","Mach", Mach)
        call flow_params_cfg%set("flow_params","Prandtl", Prandtl)

        call flow_params_cfg%write("flow_params.dat")
      endif

    endassociate

  endsubroutine set_flow_params

  subroutine set_bl_prop(self)
    class(equation_singleideal_object), intent(inout) :: self
    real(rkind) :: Re_out, s2tinf, Trat
    real(rkind), allocatable, dimension(:) :: uvec, rhovec, tvec, viscvec, yl0, yvec
    real(rkind) :: cf,thrat,retheta,redelta,th,ch,mtau,deltav,spr
    integer :: l,i,m,imode,icompute
    logical :: visc_exp

    associate(nymax => self%grid%nymax, Reynolds_friction => self%Reynolds_friction,&
    & Reynolds => self%Reynolds, l0 => self%l0, Mach => self%Mach, powerlaw_vtexp => self%powerlaw_vtexp,&
    & gam => self%gam, rfac => self%rfac, yg => self%grid%yg, theta_wall => self%theta_wall,&
    & Prandtl => self%Prandtl)
      allocate(uvec(nymax),tvec(nymax),rhovec(nymax),viscvec(nymax),yl0(nymax),yvec(nymax))
      visc_exp = .false.
      if (self%visc_model == VISC_POWER) visc_exp = .true.
      if (self%bl_case) then
        s2tinf = self%sutherland_S/self%T_ref_dim
        if (.not.self%bl_laminar) then
          if (self%masterproc) write(*,*) 'Input friction Reynolds number: ', Reynolds
          Reynolds_friction = Reynolds
          Trat = self%T_wall/self%T_recovery
          yl0 = yg(1:nymax)/l0
          imode = 0
          icompute = 0
          yvec = yg(1:nymax)
          spr = 0.8_rkind
          call hasan_meanprofile(nymax,l0,Mach,theta_wall,Reynolds_friction,retheta,redelta,Prandtl,&
          &spr,rfac,gam,powerlaw_vtexp,visc_exp,s2tinf,yvec,uvec,tvec,rhovec,viscvec,th,cf,ch,mtau,deltav,&
          &imode,icompute)
          Reynolds = redelta
        endif
        if (self%masterproc) write(*,*) 'Reynolds based on free-stream properties: ', Reynolds
      endif
    endassociate
  endsubroutine set_bl_prop

  subroutine set_chan_prop(self)
    class(equation_singleideal_object), intent(inout) :: self
    real(rkind) :: s2tinf, Re_out
    real(rkind) :: trec_over_tb, tw_over_tb, tb_over_tw, tw_over_tr
    real(rkind) :: thtmp, rmb, rmtw
    real(rkind) :: gamloc,gamlocold,gamlocgm1
    logical :: visc_exp
    integer :: l, max_iter
    associate(T_wall => self%T_wall, t0 => self%t0, theta_wall => self%theta_wall,&
    & T_bulk_target => self%T_bulk_target, indx_cp_l => self%indx_cp_l, indx_cp_r => self%indx_cp_r,&
    & cp_coeff => self%cp_coeff, calorically_perfect => self%calorically_perfect, Mach => self%Mach,&
    & gam => self%gam, rfac => self%rfac, Reynolds => self%Reynolds, Reynolds_friction => self%Reynolds_f&
    &riction, powerlaw_vtexp => self%powerlaw_vtexp, gm1 => self%gm1, rgas0 => self%rgas0,&
    & yg => self%grid%yg, nymax => self%grid%nymax, yn => self%grid%yn)
      self%T_wall = self%t0
      if (theta_wall>=-1._rkind) then
        thtmp = theta_wall
        rmb = Mach
        gamloc = gam
        max_iter = 50
        do l=1,max_iter
          gamlocold = gamloc
          gamlocgm1 = gamloc-1._rkind
          trec_over_tb = 1._rkind+gamlocgm1/2._rkind*rfac*rmb**2
          tw_over_tb = 1._rkind+thtmp*(trec_over_tb-1._rkind)
          tb_over_tw = 1._rkind/tw_over_tb
          tw_over_tr = 1._rkind/(trec_over_tb*tb_over_tw)
          T_bulk_target = tb_over_tw*T_wall
          gamloc = get_gamloc(T_bulk_target,t0,indx_cp_l,indx_cp_r,cp_coeff,calorically_perfect,rgas0)
          if (abs(gamloc-gamlocold)<1.D-9) exit
        enddo
        if (self%masterproc) write(*,*) 'Target bulk temperature: ', T_bulk_target
        if (self%masterproc) write(*,*) 'Gamma at bulk temperature: ', gamloc
        rmtw = rmb*sqrt(tb_over_tw)
        Mach = rmtw
      else
        thtmp = -1._rkind
        rmtw = Mach
        rmb = sqrt(rmtw**2/(1._rkind-gm1/2._rkind*rfac*thtmp*rmtw**2))
        trec_over_tb = 1._rkind+gm1/2._rkind*rfac*rmb**2
        tw_over_tb = 1._rkind+thtmp*(trec_over_tb-1._rkind)
        tb_over_tw = 1._rkind/tw_over_tb
        tw_over_tr = 1._rkind/(trec_over_tb*tb_over_tw)
      endif
      if (self%masterproc) write(*,*) 'Mach number (based on T bulk): ', rmb
      if (self%masterproc) write(*,*) 'Mach number (based on T wall): ', rmtw
      visc_exp = .false.
      if (self%visc_model == VISC_POWER) visc_exp = .true.
      s2tinf = self%sutherland_S/self%T_ref_dim
      if(self%grid%grid_dim == 1) then
        Reynolds_friction = Reynolds
        call get_reynolds_cha(Reynolds_friction,rmb,tw_over_tr,s2tinf,powerlaw_vtexp,visc_exp,gam,&
        &rfac, nymax/2,yn(1:nymax/2+1),Re_out)
        Reynolds = Re_out
      endif
      if (self%masterproc) write(*,*) 'Re bulk (viscosity evaluated at Twall) = ', Reynolds
    endassociate

  endsubroutine set_chan_prop
  subroutine set_oblique_shock(self)
    class(equation_singleideal_object), intent(inout) :: self
    integer :: i, m
    real(rkind) :: sins,coss,Mach_normal,Mach_normal2,vel_ratio,rho_ratio,p_ratio,deflection_angle_tan,deflection_angle
    real(rkind) :: Mach_past_shock, Mach_normal_past_shock, t_past_shock, c_past_shock, u_past_shock, v_past_shock, t_ratio
    real(rkind) :: rho_past_shock, e_past_shock, etot_past_shock, xshock_top
    real(rkind) :: xx, tanhf
    integer :: ig_shock

    associate(Mach => self%Mach, u0 => self%u0, rho0 => self%rho0, gam => self%gam, gm => self%gm,&
    & gm1 => self%gm1, nymax => self%grid%nymax, yg => self%grid%yg, winf => self%winf,&
    & indx_cp_l => self%indx_cp_l, indx_cp_r => self%indx_cp_r, xshock_imp => self%xshock_imp,&
    & shock_angle => self%shock_angle, winf_past_shock => self%winf_past_shock, nxmax => self%grid%nxmax,&
    & tanhfacs => self%tanhfacs, cfg => self%cfg, rgas0 => self%rgas0, calorically_perfect => self%calori&
    &cally_perfect, t0 => self%t0, l0 => self%l0)
      call cfg%get("flow","xshock_imp",xshock_imp)
      xshock_imp = xshock_imp*l0
      call cfg%get("flow","shock_angle",shock_angle)
      shock_angle = shock_angle*pi/180._rkind
      xshock_top = xshock_imp-yg(nymax)/tan(shock_angle)
      call locateval(self%grid%xg(1:nxmax),nxmax,xshock_top,ig_shock)
      tanhfacs = self%grid%dxg(ig_shock)
      if (self%masterproc) write(*,*) 'xshock_top, tanhfacs:', xshock_top
      sins = sin(shock_angle)
      coss = cos(shock_angle)
      Mach_normal = Mach*sins
      Mach_normal2 = Mach_normal**2
      Mach_normal_past_shock = (1._rkind+0.5_rkind*gm1*Mach_normal2)/(gam*Mach_normal2-0.5_rkind*gm1)
      Mach_normal_past_shock = sqrt(Mach_normal_past_shock)
      vel_ratio = (2._rkind+gm1*Mach_normal2)/((gam+1._rkind)*Mach_normal2)
      rho_ratio = 1._rkind/vel_ratio
      p_ratio = 1._rkind+2._rkind*gam/(gam+1._rkind)*(Mach_normal2-1._rkind)
      t_ratio = p_ratio/rho_ratio
      deflection_angle_tan = 2._rkind*coss/sins*(Mach**2*sins**2-1._rkind)/(2._rkind+Mach**2*(gam+cos(2._rkind*shock_angle)))
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
      e_past_shock = get_e_from_temperature(t_past_shock, t0, indx_cp_l, indx_cp_r, self%cv_coeff, calorically_perfect)
      etot_past_shock = e_past_shock+0.5_rkind*(u_past_shock**2+v_past_shock**2)
      winf_past_shock(1) = rho_past_shock
      winf_past_shock(2) = rho_past_shock*u_past_shock
      winf_past_shock(3) = -rho_past_shock*v_past_shock
      winf_past_shock(4) = 0._rkind
      winf_past_shock(5) = rho_past_shock*etot_past_shock
      if (self%masterproc) then
        open(18,file='shock_profile.dat')
        do i=1,nxmax
          xx = self%grid%xg(i)-xshock_top
          tanhf = 0.5_rkind*(1._rkind+tanh(xx/tanhfacs))
          write(18,100) self%grid%xg(i), ((winf(m)+tanhf*(winf_past_shock(m)-winf(m))),m=1,5)
          100 format(20ES20.10)
        enddo
        close(18)
      endif
    endassociate

  endsubroutine set_oblique_shock

  subroutine initial_conditions(self)
    class(equation_singleideal_object), intent(inout) :: self
    integer :: l
    if(self%grid%grid_dim == 2) then
      allocate(self%wfar(self%field%nx,self%nv))
      do l=1,self%nv
        self%wfar(:,l) = self%winf(l)
      enddo
    endif

    do l=1,self%nv
      self%field%w(:,:,:,l) = self%winf(l)
    enddo

    self%w_aux(:,:,:,6) = self%t0
    select case (self%flow_init)
    case(-1)
      call self%init_wind_tunnel()
    case(0)
      call self%init_channel()
    case(1)
      if (self%bl_laminar) then
        call self%init_bl_lam()
      else
        call self%init_bl()
      endif
    case(2)
      call self%init_cv()
    case(4)
      call self%init_channel_c2()
    case(5)
      call self%init_airfoil()
    endselect
  endsubroutine initial_conditions
  subroutine get_reynolds_cha(retau,rm,trat,s2tinf,vtexp,visc_exp,gam,rfac,ny,yn,rebulk)
    logical, intent(in) :: visc_exp
    integer, intent(in) :: ny
    real(rkind), intent(in) :: retau,rm,trat,s2tinf,gam,rfac,vtexp
    real(rkind), intent(out) :: rebulk
    real(rkind), dimension(ny+1), intent(in) :: yn
    integer :: j,m,itransf,jp99,jj99,jj,l
    real(rkind) :: gm1h, tr, tw, alf, retau_inc, retau_target, retau_correction, retau_musker
    real(rkind) :: retau_inc_target, retau_inc_correction, up99, d99plus_inc
    real(rkind) :: fuu, du, rhow_over_rho, rhow_over_rhojm, rmu_over_rmuw, rmu_over_rmuwjm
    real(rkind) :: ucij, ucijm, res, dy, ycij, ycijm, u99, d99, d99plus, yy, ui
    real(rkind) :: ue, uuu, uuum, thi
    real(rkind) :: rhosum,rhousum,ubulk,ysum
    real(rkind) :: vkc,bcost,c1,b_rich,eta_rich
    real(rkind), dimension(ny+1) :: yplus
    real(rkind), dimension(ny+1) :: uplus
    real(rkind), dimension(ny+1) :: ycplus
    real(rkind), dimension(ny+1) :: ucplus
    real(rkind), dimension(ny+1) :: yplusold
    real(rkind), dimension(ny+1) :: ucplusold
    real(rkind), dimension(ny+1) :: density
    real(rkind), dimension(ny+1) :: viscosity
    real(rkind), dimension(ny+1) :: temperature
    real(rkind), dimension(ny+1) :: uu,uci,uc
    itransf = 0
    gm1h = 0.5_rkind*(gam-1._rkind)
    tr = 1._rkind+gm1h*rfac*rm**2
    tw = trat*tr
    alf = 0.8259_rkind
    vkc = 0.39_rkind
    bcost = 5.5_rkind
    c1 = -1._rkind/vkc*log(vkc)+bcost
    b_rich = 0.33_rkind
    eta_rich = 11._rkind
    ycplus = 0._rkind
    do j=2,ny+1
      ycplus(j) = (yn(j)-yn(1))/(yn(ny+1)-yn(1))*retau
    enddo
    yplus = ycplus
    yplusold = yplus
    do
      uplus = 0._rkind
      do j=2,ny+1
        uplus(j) = velpiros_channel(yplus(j),retau)
      enddo
      ucplus = uplus
      do
        ucplusold = ucplus
        do j=1,ny+1
          uu(j) = ucplus(j)/ucplus(ny)
          fuu = uu(j)
          temperature(j) = tw+(tr-tw)*fuu+(1._rkind-tr)*uu(j)**2
          density(j) = 1._rkind/temperature(j)
          if (visc_exp) then
            viscosity(j) = temperature(j)**vtexp
          else
            viscosity(j) = temperature(j)**(1.5_rkind)*(1._rkind+s2tinf)/(temperature(j)+s2tinf)
          endif
        enddo
        do j=2,ny+1
          du = uplus(j)-uplus(j-1)
          rhow_over_rho = density(1)/density(j)
          rhow_over_rhojm = density(1)/density(j-1)
          rmu_over_rmuw = viscosity(j)/viscosity(1)
          rmu_over_rmuwjm = viscosity(j-1)/viscosity(1)
          if (itransf==1) then
            ucij = sqrt(rhow_over_rho)
            ucijm = sqrt(rhow_over_rhojm)
          else
            ucij = sqrt(rhow_over_rho)*sqrt(rmu_over_rmuw)
            ucijm = sqrt(rhow_over_rhojm)*sqrt(rmu_over_rmuwjm)
          endif
          ucplus(j) = ucplus(j-1)+0.5_rkind*(ucij+ucijm)*du
        enddo
        res = 0._rkind
        do j=1,ny+1
          res = res+abs(ucplus(j)-ucplusold(j))
        enddo
        res = res/(ny+1)
        if (res<tol_iter) exit
      enddo
      yplus = 0._rkind
      do j=2,ny+1
        dy = ycplus(j)-ycplus(j-1)
        rhow_over_rho = density(1)/density(j)
        rhow_over_rhojm = density(1)/density(j-1)
        rmu_over_rmuw = viscosity(j)/viscosity(1)
        rmu_over_rmuwjm = viscosity(j-1)/viscosity(1)
        if (itransf==1) then
          ycij = 1._rkind
          ycijm = 1._rkind
        else
          ycij = sqrt(rhow_over_rho)*sqrt(rmu_over_rmuw**3)
          ycij = 1._rkind/ycij
          ycijm = sqrt(rhow_over_rhojm)*sqrt(rmu_over_rmuwjm**3)
          ycijm = 1._rkind/ycijm
        endif
        yplus(j) = yplus(j-1)+0.5_rkind*(ycij+ycijm)*dy
      enddo
      res = 0._rkind
      do j=1,ny+1
        res = abs(yplus(j)-yplusold(j))
      enddo
      res = res/(ny+1)
      if (res<tol_iter) exit
      yplusold = yplus
    enddo
    rhosum = 0._rkind
    rhousum = 0._rkind
    ysum = 0._rkind
    do j=1,ny
      dy = yn(j+1)-yn(j)
      ysum = ysum+dy
      rhosum = rhosum+0.5*(density(j)+density(j+1))*dy
      rhousum = rhousum+0.5*(density(j)*ucplus(j)+density(j+1)*ucplus(j+1))*dy
    enddo
    rhosum = rhosum/ysum
    rhousum = rhousum/ysum
    ubulk = rhousum/rhosum
    rebulk = retau*ubulk*rhosum/density(1)
    return
  end subroutine get_reynolds_cha
  subroutine init_bl_lam(self)
    class(equation_singleideal_object), intent(inout) :: self

    integer :: i,j,k,n,nst,ii
    integer :: ne
    real(rkind) :: etad, deta, Tinf_dim, Twall_dim, etaedge, xbl
    real(rkind) :: xx, yy, etast, wl, wr, ust, vst, tst, rst
    real(rkind) :: rho,uu,vv,ww,rhouu,rhovv,rhoww,ee,tt
    real(rkind), allocatable, dimension(:) :: ubl,tbl,vbl,eta

    allocate(self%wmean(1-self%grid%ng:self%field%nx+self%grid%ng+1, 1:self%field%ny, 4))

    associate(nx => self%field%nx, ny => self%field%ny, nz => self%field%nz, ng => self%grid%ng,&
    & x => self%field%x, y => self%field%y, z => self%field%z, Reynolds => self%Reynolds,&
    & xg => self%grid%xg, nxmax => self%grid%nxmax, rho0 => self%rho0, u0 => self%u0, p0 => self%p0,&
    & gm => self%gm, wmean => self%wmean, w => self%field%w, Mach => self%Mach, gam => self%gam,&
    & rfac => self%rfac, Prandtl => self%Prandtl, w_aux => self%w_aux, deltavec => self%deltavec ,&
    &deltavvec => self%deltavvec, cfvec => self%cfvec, t0 => self%t0, T_ref_dim => self%T_ref_dim,&
    & T_wall => self%T_wall, l0 => self%l0, indx_cp_l => self%indx_cp_l, indx_cp_r => self%indx_cp_r,&
    & cv_coeff => self%cv_coeff, calorically_perfect => self%calorically_perfect, rgas0 => self%rgas0)
      etad = 20._rkind
      deta = 0.02_rkind
      ne = nint(etad/deta)
      allocate(ubl(0:ne))
      allocate(vbl(0:ne))
      allocate(tbl(0:ne))
      allocate(eta(0:ne))
      Tinf_dim = t0*T_ref_dim/t0
      Twall_dim = T_wall*T_ref_dim/t0
      call compressible_blasius(ne,etad,deta,gam,Mach,Prandtl,Tinf_dim,Twall_dim,eta,ubl,vbl,tbl,etaedge)
      xbl = Reynolds/(etaedge**2)*l0
      if (self%masterproc) write(*,*) 'XBL', xbl
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
      do k=1,nz
        do j=1,ny
          do i=1,nx
            rho = wmean(i,j,1)
            uu = wmean(i,j,2)/rho
            vv = wmean(i,j,3)/rho
            ww = wmean(i,j,4)/rho
            rhouu = rho*uu
            rhovv = rho*vv
            rhoww = rho*ww
            w(i,j,k,1) = rho
            w(i,j,k,2) = rhouu
            w(i,j,k,3) = rhovv
            w(i,j,k,4) = rhoww
            tt = p0/rgas0/rho
            w_aux(i,j,k,6) = tt
            ee = get_e_from_temperature(tt, t0, indx_cp_l, indx_cp_r, cv_coeff, calorically_perfect)
            w(i,j,k,5) = rho*ee + 0.5_rkind*(rhouu**2+rhovv**2+rhoww**2)/rho
          enddo
        enddo
      enddo

    endassociate
  endsubroutine init_bl_lam
  subroutine init_bl(self)
    class(equation_singleideal_object), intent(inout) :: self

    real(rkind), dimension(1-self%grid%ng:self%grid%nxmax+self%grid%ng+1) :: thvec, retauvec
    real(rkind), dimension(self%field%ny) :: yvec, uvec, rhovec, tvec, viscvec
    real(rkind), dimension(3) :: rr
    real(rkind) :: spr
    real(rkind) :: cf,ch,mtau,deltav
    real(rkind) :: retheta,retau,redelta,retauold,delta,th,retheta_inflow,deltaold
    real(rkind) :: vi,vi_j,vi_jm
    real(rkind) :: rho,uu,vv,ww,rhouu,rhovv,rhoww,ee,tt
    real(rkind) :: u0_02
    real(rkind) :: s2tinf,vtexp
    integer :: i,j,k,ii,imode,icompute
    logical :: visc_exp
    logical :: file_exists

    allocate(self%wmean(1-self%grid%ng:self%field%nx+self%grid%ng+1, 1:self%field%ny, 4))
    allocate(self%deltavec(1-self%grid%ng:self%grid%nxmax+self%grid%ng+1))
    allocate(self%deltavvec(1-self%grid%ng:self%grid%nxmax+self%grid%ng+1))
    allocate(self%cfvec(1-self%grid%ng:self%grid%nxmax+self%grid%ng+1))

    associate(nx => self%field%nx, ny => self%field%ny, nz => self%field%nz, ng => self%grid%ng,&
    & y => self%field%y, z => self%field%z, Reynolds_friction => self%Reynolds_friction,&
    & xg => self%grid%xg, nxmax => self%grid%nxmax, rho0 => self%rho0, u0 => self%u0, p0 => self%p0,&
    & gm => self%gm, l0 => self%l0, wmean => self%wmean, w => self%field%w, Mach => self%Mach,&
    & gam => self%gam, rfac => self%rfac, Prandtl => self%Prandtl, w_aux => self%w_aux,&
    & deltavec => self%deltavec ,deltavvec => self%deltavvec, cfvec => self%cfvec, indx_cp_l => self%indx&
    &_cp_l, indx_cp_r => self%indx_cp_r, cv_coeff => self%cv_coeff, calorically_perfect => self%calorical&
    &ly_perfect, rgas0 => self%rgas0, t0 => self%t0, jbl_inflow => self%jbl_inflow, theta_wall => self%th&
    &eta_wall)

      call locateval(y(1:ny),ny,l0,jbl_inflow)

      spr = 0.8_rkind
      yvec = y(1:ny)

      s2tinf = self%sutherland_S/self%T_ref_dim
      vtexp = self%powerlaw_vtexp
      visc_exp = .false.
      if (self%visc_model == VISC_POWER) visc_exp = .true.
      inquire(file='blvec.bin',exist=file_exists)
      if (file_exists) then
        open(183,file='blvec.bin',form='unformatted')
        read(183) cfvec,thvec,deltavec,deltavvec
        close(183)
      else
        imode = 0
        icompute = 0
        call hasan_meanprofile(ny,l0,Mach,theta_wall,Reynolds_friction,retheta,redelta,Prandtl,spr,&
        &rfac,gam,vtexp,visc_exp,s2tinf,yvec,uvec,tvec,rhovec,viscvec,th,cf,ch,mtau,deltav,imode,icompute)
        retheta_inflow = retheta
        deltavec(1) = l0
        deltavvec(1) = deltav
        cfvec(1) = cf
        thvec(1) = th
        retauvec(1) = Reynolds_friction

        retheta = retheta_inflow
        do i=1,ng
          thvec(1-i) = thvec(2-i)-0.5_rkind*abs((xg(1-i)-xg(2-i)))*cfvec(2-i)
          retheta = retheta/thvec(2-i)*thvec(1-i)
          delta = deltavec(2-i)/thvec(2-i)*thvec(1-i)
          do
            deltaold = delta
            imode = 1
            icompute = 0
            call hasan_meanprofile(ny,delta,Mach,theta_wall,retau,retheta,redelta,Prandtl,spr,rfac,&
            &gam,vtexp,visc_exp,s2tinf,yvec,uvec,tvec,rhovec,viscvec,th,cf,ch,mtau,deltav,imode,icompute)
            delta = deltavec(2-i)*retau/retauvec(2-i)*sqrt(cfvec(2-i)/cf)
            if (abs(delta-deltaold)<0.000000001_rkind) exit
          enddo
          deltavec (1-i) = delta
          cfvec (1-i) = cf
          deltavvec(1-i) = deltav
          retauvec (1-i) = retau
        enddo

        if (self%masterproc) open(182,file='cfstart.dat')
        if (self%masterproc) write(182,100) xg(1),deltavec(1),deltavvec(1),cfvec(1),thvec(1),retau,retheta,redelta
        retheta = retheta_inflow
        do i=2,nxmax+ng+1
          thvec(i) = thvec(i-1)+0.5_rkind*abs((xg(i)-xg(i-1)))*cfvec(i-1)
          retheta = retheta/thvec(i-1)*thvec(i)
          delta = deltavec(i-1)/thvec(i-1)*thvec(i)
          do
            deltaold = delta
            imode = 1
            icompute = 0
            call hasan_meanprofile(ny,delta,Mach,theta_wall,retau,retheta,redelta,Prandtl,spr,rfac,&
            &gam,vtexp,visc_exp,s2tinf,yvec,uvec,tvec,rhovec,viscvec,th,cf,ch,mtau,deltav,imode,icompute)
            delta = deltavec(i-1)*retau/retauvec(i-1)*sqrt(cfvec(i-1)/cf)
            if (abs(delta-deltaold)<0.000000001_rkind) exit
          enddo
          deltavec(i) = delta
          cfvec (i) = cf
          deltavvec(i) = deltav
          retauvec (i) = retau
          if (self%masterproc) write(182,100) xg(i),delta,deltav,cf,th,retau,retheta,redelta
          100 format(20ES20.10)
        enddo

        if (self%masterproc) close(182)
        if (self%masterproc) then
          open(183,file='blvec.bin',form='unformatted')
          write(183) cfvec,thvec,deltavec,deltavvec
          close(183)
        endif
      endif
      wmean = 0._rkind
      do i=1-ng,nx+ng+1
        ii = self%field%ncoords(1)*nx+i
        delta = deltavec(ii)
        deltav = deltavvec(ii)
        retau = delta/deltav
        imode = 0
        icompute = 1
        call hasan_meanprofile(ny,delta,Mach,theta_wall,retau,retheta,redelta,Prandtl,spr,rfac,gam,&
        &vtexp,visc_exp,s2tinf,yvec,uvec,tvec,rhovec,viscvec,th,cf,ch,mtau,deltav,imode,icompute)
        do j=1,ny
          wmean(i,j,1) = rhovec(j)*self%rho0
          wmean(i,j,2) = rhovec(j)*uvec(j)*self%u0*self%rho0
        enddo
      enddo
      do i=1-ng,nx+ng
        ii = self%field%ncoords(1)*nx+i
        do j=2,ny
          vi_j = -(wmean(i+1,j,2)-wmean(i,j,2))/(xg(ii+1)-xg(ii))
          vi_jm = -(wmean(i+1,j-1,2)-wmean(i,j-1,2))/(xg(ii+1)-xg(ii))
          vi = 0.5_rkind*(vi_j+vi_jm)
          wmean(i,j,3) = wmean(i,j-1,3)+vi*(y(j)-y(j-1))
        enddo
      enddo

      u0_02 = 0.02_rkind*u0
      if (self%rand_type==0) u0_02 = 0._rkind
      do k=1,nz
        do j=1,ny
          do i=1,nx
            rho = wmean(i,j,1)
            call get_crandom_f(rr)
            rr = rr-0.5_rkind
            uu = wmean(i,j,2)/rho+u0_02*rr(1)
            vv = wmean(i,j,3)/rho+u0_02*rr(2)
            ww = wmean(i,j,4)/rho+u0_02*rr(3)

            rhouu = rho*uu
            rhovv = rho*vv
            rhoww = rho*ww
            w(i,j,k,1) = rho
            w(i,j,k,2) = rhouu
            w(i,j,k,3) = rhovv
            w(i,j,k,4) = rhoww
            tt = p0/rgas0/rho
            w_aux(i,j,k,6) = tt
            ee = get_e_from_temperature(tt, t0, indx_cp_l, indx_cp_r, cv_coeff, calorically_perfect)
            w(i,j,k,5) = rho*ee + 0.5_rkind*(rhouu**2+rhovv**2+rhoww**2)/rho
          enddo
        enddo
      enddo
      if (self%masterproc) then
        open(183,file='blinflow.dat',form='formatted')
        do j=1,ny
          write(183,*) y(j),wmean(1,j,1),wmean(1,j,2)/wmean(1,j,1)
        enddo
        close(183)
      endif

      call self%add_synthetic_perturbations()

    endassociate
  endsubroutine init_bl
  subroutine init_bl_old(self)
    class(equation_singleideal_object), intent(inout) :: self

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

    allocate(self%wmean(1-self%grid%ng:self%field%nx+self%grid%ng+1, 1:self%field%ny, 4))
    allocate(self%deltavec(1-self%grid%ng:self%grid%nxmax+self%grid%ng+1))
    allocate(self%deltavvec(1-self%grid%ng:self%grid%nxmax+self%grid%ng+1))
    allocate(self%cfvec(1-self%grid%ng:self%grid%nxmax+self%grid%ng+1))

    associate(nx => self%field%nx, ny => self%field%ny, nz => self%field%nz, ng => self%grid%ng,&
    & y => self%field%y, z => self%field%z, Reynolds_friction => self%Reynolds_friction,&
    & xg => self%grid%xg, nxmax => self%grid%nxmax, rho0 => self%rho0, u0 => self%u0, p0 => self%p0,&
    & gm => self%gm, l0 => self%l0, wmean => self%wmean, w => self%field%w, Mach => self%Mach,&
    & gam => self%gam, rfac => self%rfac, Prandtl => self%Prandtl, w_aux => self%w_aux,&
    & deltavec => self%deltavec ,deltavvec => self%deltavvec, cfvec => self%cfvec, indx_cp_l => self%indx&
    &_cp_l, indx_cp_r => self%indx_cp_r, cv_coeff => self%cv_coeff, calorically_perfect => self%calorical&
    &ly_perfect, rgas0 => self%rgas0, t0 => self%t0, jbl_inflow => self%jbl_inflow)

      call locateval(y(1:ny),ny,l0,jbl_inflow)

      Trat = self%T_wall/self%T_recovery

      deltavec(1) = l0
      deltavvec(1) = l0/Reynolds_friction
      yl0 = y(1:ny)/l0

      s2tinf = self%sutherland_S/self%T_ref_dim
      vtexp = self%powerlaw_vtexp
      visc_exp = .false.
      if (self%visc_model == VISC_POWER) visc_exp = .true.
      inquire(file='blvec.bin',exist=file_exists)
      if (file_exists) then
        open(183,file='blvec.bin',form='unformatted')
        read(183) cfvec,thvec,deltavec,deltavvec
        close(183)
      else
        call meanvelocity_bl(Reynolds_friction,Mach,Trat,s2tinf,vtexp,visc_exp,gam,rfac,ny,&
        &yl0(1:ny),uvec,rhovec,tvec,viscvec,redelta,cf,thrat)
        cfvec(1) = cf
        thvec(1) = thrat*deltavec(1)
        retauold = Reynolds_friction

        do i=1,ng
          retau = retauold
          do
            call meanvelocity_bl(retau,Mach,Trat,s2tinf,vtexp,visc_exp,gam,rfac,ny,yl0(1:ny),uvec,&
            &rhovec,tvec,viscvec,redelta,cf,thrat)
            th = thvec(2-i)-0.25_rkind*abs((xg(1-i)-xg(2-i)))*(cf+cfvec(2-i))
            delta = th/thrat
            deltav = deltavvec(2-i)*sqrt(cfvec(2-i)/cf)
            retau = delta/deltav
            if (abs(retau-retauold) < 0.01_rkind) exit
            retauold = retau
          enddo
          thvec (1-i) = th
          cfvec (1-i) = cf
          deltavec (1-i) = delta
          deltavvec(1-i) = deltav
        enddo

        retauold = Reynolds_friction

        if (self%masterproc) open(182,file='cfstart.dat')
        do i=2,nxmax+ng+1
          retau = retauold
          do
            call meanvelocity_bl(retau,Mach,Trat,s2tinf,vtexp,visc_exp,gam,rfac,ny,yl0(1:ny),uvec,&
            &rhovec,tvec,viscvec,redelta,cf,thrat)
            th = thvec(i-1)+0.25_rkind*(xg(i)-xg(i-1))*(cf+cfvec(i-1))
            delta = th/thrat
            deltav = deltavvec(i-1)*sqrt(cfvec(i-1)/cf)
            retau = delta/deltav
            if (abs(retau-retauold)<0.01_rkind) exit
            retauold = retau
          enddo
          thvec (i) = th
          cfvec (i) = cf
          deltavec (i) = delta
          deltavvec(i) = deltav
          if (self%masterproc) write(182,100) xg(i),delta,deltav,cf,th
          100 format(20ES20.10)
        enddo
        if (self%masterproc) close(182)
        if (self%masterproc) then
          open(183,file='blvec.bin',form='unformatted')
          write(183) cfvec,thvec,deltavec,deltavvec
          close(183)
        endif
      endif
      wmean = 0._rkind
      do i=1-ng,nx+ng+1
        ii = self%field%ncoords(1)*nx+i
        delta = deltavec(ii)
        deltav = deltavvec(ii)
        retau = delta/deltav
        call meanvelocity_bl(retau,Mach,Trat,s2tinf,vtexp,visc_exp,gam,rfac,ny,yl0(1:ny),uvec,rhovec,tvec,viscvec,redelta,cf,thrat)
        do j=1,ny
          yl = y(j)/delta
          call locateval(yl0(1:ny),ny,yl,jj)
          m = 4
          jjj = min(max(jj-(m-1)/2,1),ny+1-m)
          call pol_int(yl0(jjj),rhovec(jjj),m,yl,rhoi)
          call pol_int(yl0(jjj),uvec(jjj),m,yl,ui )
          call pol_int(yl0(jjj),tvec(jjj),m,yl,ti )
          wmean(i,j,1) = rhoi*self%rho0
          wmean(i,j,2) = rhoi*ui*self%u0*self%rho0
        enddo
      enddo
      do i=1-ng,nx+ng
        ii = self%field%ncoords(1)*nx+i
        do j=2,ny
          vi_j = -(wmean(i+1,j,2)-wmean(i,j,2))/(xg(ii+1)-xg(ii))
          vi_jm = -(wmean(i+1,j-1,2)-wmean(i,j-1,2))/(xg(ii+1)-xg(ii))
          vi = 0.5_rkind*(vi_j+vi_jm)
          wmean(i,j,3) = wmean(i,j-1,3)+vi*(y(j)-y(j-1))
        enddo
      enddo

      u0_02 = 0.02_rkind*u0
      if (self%rand_type==0) u0_02 = 0._rkind
      do k=1,nz
        do j=1,ny
          do i=1,nx
            rho = wmean(i,j,1)
            call get_crandom_f(rr)
            rr = rr-0.5_rkind
            uu = wmean(i,j,2)/rho+u0_02*rr(1)
            vv = wmean(i,j,3)/rho+u0_02*rr(2)
            ww = wmean(i,j,4)/rho+u0_02*rr(3)

            rhouu = rho*uu
            rhovv = rho*vv
            rhoww = rho*ww
            w(i,j,k,1) = rho
            w(i,j,k,2) = rhouu
            w(i,j,k,3) = rhovv
            w(i,j,k,4) = rhoww
            tt = p0/rgas0/rho
            w_aux(i,j,k,6) = tt
            ee = get_e_from_temperature(tt, t0, indx_cp_l, indx_cp_r, cv_coeff, calorically_perfect)
            w(i,j,k,5) = rho*ee + 0.5_rkind*(rhouu**2+rhovv**2+rhoww**2)/rho
          enddo
        enddo
      enddo

      call self%add_synthetic_perturbations()

    endassociate
  endsubroutine init_bl_old

  subroutine add_synthetic_perturbations(self)
    class(equation_singleideal_object), intent(inout) :: self
    real(rkind), dimension(:,:), allocatable :: synth_params
    real(rkind) :: rho_wall, tau_wall, u_tau, lz_plus, rr, rhofac, up, vp, wp
    real(rkind) :: arg_sin, arg_cos, ys , ufy , vfy , dup_dx , dvfy_dy, dvp_dy
    real(rkind) :: rho,ee,tt
    real(rkind) :: delta, u0_03
    real(rkind), dimension(3) :: rr3
    integer :: i,j,k,l,n_streaks,ii

    allocate(synth_params(5,7))

    synth_params(:,1) = [ 12._rkind, 0.3_rkind, 0.45_rkind, 0.6_rkind, 0.5_rkind]
    synth_params(:,2) = [ 1.2_rkind, 0.3_rkind, 0.2_rkind, 0.08_rkind, 0.04_rkind]
    synth_params(:,3) = [-0.25_rkind, -0.06_rkind, -0.05_rkind, -0.04_rkind, -0.03_rkind]
    synth_params(:,4) = [ 0.12_rkind, 1.2_rkind, 0.6_rkind, 0.4_rkind, 0.2_rkind]
    synth_params(:,5) = [ 10.0_rkind, 0.9_rkind, 0.9_rkind, 0.9_rkind, 0.9_rkind]
    synth_params(:,6) = [120.0_rkind, 0.333_rkind, 0.25_rkind, 0.2_rkind, 0.166_rkind]
    synth_params(:,7) = [ 0.0_rkind, 1._rkind, 1._rkind, 1._rkind, 1._rkind]

    associate(nx => self%field%nx, ny => self%field%ny, nz => self%field%nz, ng => self%grid%ng,&
    & deltavec => self%deltavec, deltavvec => self%deltavvec, wmean => self%wmean, cfvec => self%cfvec,&
    & rho0 => self%rho0, u0 => self%u0, lz => self%grid%domain_size(3), w => self%field%w, l0 => self%l0,&
    & x => self%field%x, y => self%field%y, z => self%field%z, p0 => self%p0, gm => self%gm,&
    & w_aux => self%w_aux, indx_cp_l => self%indx_cp_l, indx_cp_r => self%indx_cp_r, cv_coeff => self%cv_&
    &coeff, calorically_perfect => self%calorically_perfect, rgas0 => self%rgas0, t0 => self%t0)

      rho_wall = wmean(1,1,1)
      tau_wall = cfvec(1)*rho0*u0**2*0.5_rkind
      u_tau = sqrt(tau_wall/rho_wall)

      synth_params(1,1) = synth_params(1,1) * deltavvec(1)
      synth_params(1,4) = synth_params(1,4) * u_tau / deltavvec(1)
      synth_params(1,5) = synth_params(1,5) * u_tau

      lz_plus = lz/deltavvec(1)
      n_streaks = nint(lz_plus / synth_params(1,6))
      synth_params(1,6) = lz / n_streaks

      do l=2,5
        synth_params(l,1) = synth_params(l,1) * deltavec(1)
        synth_params(l,4) = synth_params(l,4) * u0 / deltavec(1)
        synth_params(l,5) = synth_params(l,5) * u0
        synth_params(l,6) = synth_params(l,6) * lz
        call get_crandom_f(rr)
        synth_params(l,7) = synth_params(l,7) * 2._rkind*pi *rr
      enddo
      call mpi_bcast(synth_params(2:5,7),4,mpi_prec,0,self%field%mp_cart,self%mpi_err)

      u0_03 = 0.03_rkind*u0
      if (self%rand_type==0) u0_03 = 0._rkind
      do k=1,nz
        do j=1,ny
          do i=1,nx
            rhofac = wmean(i,1,1)/wmean(i,j,1)
            rhofac = sqrt(rhofac)*u0
            up = 0._rkind
            vp = 0._rkind
            wp = 0._rkind
            do l=2,5
              arg_sin = synth_params(l,4)*x(i)/synth_params(l,5)
              arg_cos = 2._rkind*pi*z(k)/synth_params(l,6) + synth_params(l,7)
              ys = y(j)/synth_params(l,1)
              ufy = ys * exp(-ys)
              vfy = ys**2 * exp(-(ys**2))
              up = up+synth_params(l,2) * ufy * sin(arg_sin) * cos(arg_cos)
              vp = vp+synth_params(l,3) * vfy * sin(arg_sin) * cos(arg_cos)

              dup_dx = synth_params(l,2) * ufy * synth_params(l,4)/synth_params(l,5) * cos(arg_sin)
              dvfy_dy = 2._rkind*ys/synth_params(l,1)*exp(-(ys**2))*(1._rkind-ys**2)
              dvp_dy = synth_params(l,3) * dvfy_dy * sin(arg_sin)
              wp = wp-(dup_dx+dvp_dy)*sin(arg_cos)*synth_params(l,6)/(2._rkind*pi)
            enddo
            up = up * rhofac
            vp = vp * rhofac
            wp = wp * rhofac
            ii = self%field%ncoords(1)*nx+i
            delta = deltavec(ii)
            if (y(j)<delta) then
              call get_crandom_f(rr3)
              rr3 = rr3-0.5_rkind
              up = up+u0_03*rr3(1)*(y(j)/l0)
              vp = vp+u0_03*rr3(2)*(y(j)/l0)
              wp = wp+u0_03*rr3(3)*(y(j)/l0)
            endif
            rho = w(i,j,k,1)
            w(i,j,k,2) = w(i,j,k,2) + up*w(i,j,k,1)
            w(i,j,k,3) = w(i,j,k,3) + vp*w(i,j,k,1)
            w(i,j,k,4) = w(i,j,k,4) + wp*w(i,j,k,1)
            tt = p0/rgas0/rho
            w_aux(i,j,k,6) = tt
            ee = get_e_from_temperature(tt, t0, indx_cp_l, indx_cp_r, cv_coeff, calorically_perfect)
            w(i,j,k,5) = rho*ee + 0.5_rkind*(w(i,j,k,2)**2 + w(i,j,k,3)**2 + w(i,j,k,4)**2)/w(i,j,k,1)
          enddo
        enddo
      enddo

    endassociate
  endsubroutine add_synthetic_perturbations
  function velmusker(yplus,retau)
    real(rkind), intent(in) :: yplus,retau
    real(rkind) :: yp,eta,pi_wake
    real(rkind) :: velmusker
    yp = yplus
    eta = yp/retau
    eta = min(1._rkind,eta)
    yp = eta*retau
    pi_wake = 0.434_rkind
    velmusker = 5.424_rkind*atan((2*yp-8.15_rkind)/16.7_rkind)+ log10((yp+10.6_rkind)**9.6_rkind/&
    &(yp**2-8.15_rkind*yp+86)**2)-3.51132976630723_rkind +2.44_rkind*(pi_wake*(6*eta**2-4*eta**3)+(eta**&
    &2*(1-eta)))
  end function velmusker
  function velpiros_channel(yplus,retau)
    real(rkind), intent(in) :: yplus,retau
    real(rkind) :: velpiros_channel
    real(rkind) :: yp,eta,vkc,vkci,c,c1,pi_wake,etastar,b,b_rich,eta_rich,uinn,uout
    yp = yplus
    eta = yp/retau
    eta = min(1._rkind,eta)
    yp = eta*retau
    c = 5.3_rkind
    vkc = 0.39_rkind
    vkci = 1._rkind/vkc
    c1 = -1._rkind/vkc*log(vkc)+c
    pi_wake = 0.00857_rkind
    etastar = 0.275_rkind
    b = 0.0447_rkind
    b_rich = 0.33_rkind
    eta_rich = 11._rkind
    uinn = vkci*log(1._rkind+vkc*yp)+c1*(1._rkind-exp(-yp/eta_rich)-yp/eta_rich*exp(-b_rich*yp))
    uout = c+vkci*log(retau)+vkci*log(etastar)
    uout = uout+(1-etastar)/(2*vkc*etastar)*(1-(1-eta)**2/(1-etastar)**2)
    if (eta.le.etastar) then
      velpiros_channel = uinn
    else
      velpiros_channel = uout
    endif
    return
  end function velpiros_channel
  subroutine compressible_blasius(n,etad,deta,gam,Mach,Prandtl,Te,Twall,eta,u,v,tbl,delta1)
    integer, intent(in) :: n
    real(rkind), intent(in) :: etad, deta, Mach, Prandtl, Te, Twall, gam
    real(rkind), intent(out) :: delta1
    real(rkind), dimension(0:n), intent(inout) :: u,v,tbl,eta
    real(rkind), dimension(0:n) :: f,t,g,h,a1,a2,a3,a4,a5,a6,a7,a8,s,r,visc,tcr
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
    Rg = 287.15_rkind
    s_suth = 110.4_rkind
    eps = 0.000001_rkind
    z1 = 0.334_rkind
    z2 = 0.82_rkind
    z3 = 0.22_rkind
    iflag = 1
    nm = n-1
    gm = gam / ( gam -1._rkind)
    adiab = (1._rkind+ (gam - 1._rkind) / 2._rkind *Mach*Mach)
    T0e = Te * adiab
    Ae = sqrt (gam * Rg * Te)
    Ue = Mach * ae
    Cp = gm * Rg
    He = Cp * Te
    H0e = Cp * T0e
    if (Te.lt.110.4_rkind) then
      visce = .693873D-6*te
    else
      visce = 1.458D-5 * te**1.5_rkind/( te + s_suth )
    end if
    j = 0
    eta(0) = 0._rkind
    u(0) = 0._rkind
    h(0) = 0._rkind
    a1(0) = 0._rkind
    a2(0) = 0._rkind
    a3(0) = 0._rkind
    a5(0) = 1._rkind
    a6(0) = 0._rkind
    a7(0) = 0._rkind
    s0 = z1
    if(iflag.eq.0) then
      g(0) = 0._rkind
      a4(0) = 1._rkind
      a8(0) = 0._rkind
      r0 = z2
    else
      t(0) = (Twall - Te )/(T0e-te)
      a4(0) = 0._rkind
      a8(0) = 1._rkind
      r0 = z3
    end if
    do

      f(0) = s0
      if ( iflag.eq.0) t(0) = r0
      if ( iflag.eq.1) g(0) = r0
      Tbl(0) = Te + (T0e - Te) * t(0)
      tt = tbl(0)
      if (tt<110.4_rkind) then
        visc(0) = .693873D-6*tt
      else
        visc(0) = 1.458D-5 * tt**1.5_rkind / ( tt + s_suth )
      end if
      vis = visc(0) / visce
      do i = 0,nm
        u0 = f(i) / vis * deta
        f0 = - h(i) * f(i) / vis* deta
        t0 = g(i) * Prandtl / vis * deta
        g0 = -(h(i) * g(i) * Prandtl + 2._rkind * f(i) **2._rkind) / vis * deta
        h0 = .5_rkind*u(i) / (1 + (T0e/Te - 1) * t(i)) * deta
        u1 = (f(i) + .5_rkind*f0) / vis * deta
        f1 = - (h(i) + .5_rkind*h0) * (f(i) + .5_rkind*f0) /vis * deta
        t1 = (g(i) + .5_rkind*g0) * Prandtl / vis * deta
        g1 = -((h(i) + .5_rkind*h0) * (g(i) + .5_rkind*g0) * Prandtl + 2._rkind * (f(i) + .5_rkind*f0)**2._rkind) /vis*deta
        h1 = .5_rkind*(u(i)+.5_rkind*u0) / (1 + (T0e/Te - 1) * (t(i)+.5_rkind*t0)) *deta
        u2 = (f(i) + .5_rkind*f1) / vis * deta
        f2 = - (h(i) + .5_rkind*h1) * (f(i) + .5_rkind*f1) /vis * deta
        t2 = (g(i) + .5_rkind*g1) * Prandtl / vis * deta
        g2 = -((h(i) + .5_rkind*h1) * (g(i) + .5_rkind*g1) * Prandtl + 2._rkind * (f(i) + .5_rkind*f1)**2._rkind) /vis*deta
        h2 = .5_rkind*(u(i)+.5_rkind*u1) / (1 + (T0e/Te - 1) * (t(i)+.5*t1)) *deta
        u3 = (f(i) + f2) / vis * deta
        f3 = - (h(i) + h2) * (f(i) + f2) /vis * deta
        t3 = (g(i) + g2) * Prandtl / vis * deta
        g3 = -((h(i) + h2) * (g(i) + g2) * Prandtl + 2._rkind * (f(i) + f2)**2._rkind) /vis*deta
        h3 = .5_rkind*(u(i)+u2) / (1 + (T0e/Te - 1) *(t(i)+t2))*deta
        a10 = a5(i)/vis *deta
        a20 = a6(i) / vis *deta
        a30 = a7(i) * Prandtl / vis *deta
        a40 = a8(i) *Prandtl / vis *deta
        a50 = - a5(i) * h(i) / vis *deta
        a60 = - a6(i) * h(i) / vis *deta
        a70 = -(4*f(i)*a5(i)+ Prandtl * h(i) *a7(i))/vis *deta
        a80 = -(4*f(i)*a6(i)+ Prandtl * h(i)*a8(i))/vis *deta
        a11 = (a5(i) + .5_rkind*a50) / vis *deta
        a21 = (a6(i) + .5_rkind*a60) / vis *deta
        a31 = (a7(i) + .5_rkind*a70) * Prandtl / vis *deta
        a41 = (a8(i) + .5_rkind*a80) * Prandtl / vis *deta
        a51 = - (a5(i) + .5_rkind*a50) * (h(i) + .5_rkind*h0) / vis *deta
        a61 = - (a6(i) + .5_rkind*a60) * (h(i) + .5_rkind*h0) / vis *deta
        a71 = -(4* (f(i) + .5_rkind*f0) * (a5(i) + .5_rkind*a50)+Prandtl * (h(i) + .5_rkind*h0) *(a7(i) + .5_rkind*a70))/vis *deta
        a81 = -(4* (f(i) + .5_rkind*f0) * (a6(i) + .5_rkind*a60)+Prandtl * (h(i) + .5_rkind*h0) *(a8(i) + .5_rkind*a80))/vis *deta
        a12 = (a5(i) + .5_rkind*a51) / vis *deta
        a22 = (a6(i) + .5_rkind*a61) / vis *deta
        a32 = (a7(i) + .5_rkind*a71) * Prandtl / vis *deta
        a42 = (a8(i) + .5_rkind*a81) * Prandtl / vis *deta
        a52 = - (a5(i) + .5_rkind*a51) * (h(i) + .5_rkind*h1) / vis *deta
        a62 = - (a6(i) + .5_rkind*a61) * (h(i) + .5_rkind*h1) / vis *deta
        a72 = -(4* (f(i) + .5_rkind*f1) * (a5(i) + .5_rkind*a51)+Prandtl * (h(i) + .5_rkind*h1) *(a7(i) + .5_rkind*a71))/vis *deta
        a82 = -(4* (f(i) + .5_rkind*f1) * (a6(i) + .5_rkind*a61)+Prandtl * (h(i) + .5_rkind*h1) *(a8(i) + .5_rkind*a81))/vis *deta
        a13 = (a5(i) + a52) / vis *deta
        a23 = (a6(i) + a62) / vis *deta
        a33 = (a7(i) + a72) * Prandtl / vis *deta
        a43 = (a8(i) + a82) * Prandtl / vis *deta
        a53 = - (a5(i) + a52) * (h(i) + .5_rkind*h2) / vis *deta
        a63 = - (a6(i) + a62) * (h(i) + .5_rkind*h2) / vis *deta
        a73 = -(4* (f(i) + f2) * (a5(i) + a52) + Prandtl * (h(i) + h2) *(a7(i) + a72))/vis *deta
        a83 = -(4* (f(i) + f2) * (a6(i) + a62) + Prandtl * (h(i) + h2) *(a8(i) + a82))/vis *deta
        f(i+1) = f(i) + (f0 + 2._rkind*f1 + 2*f2 + f3) / 6._rkind
        u(i+1) = u(i) + (u0 + 2._rkind*u1 + 2*u2 + u3) / 6._rkind
        t(i+1) = t(i) + (t0 + 2._rkind*t1 + 2*t2 + t3) / 6._rkind
        g(i+1) = g(i) + (g0 + 2._rkind*g1 + 2*g2 + g3) / 6._rkind
        h(i+1) = h(i) + (h0 + 2._rkind*h1 + 2*h2 + h3) / 6._rkind
        a1(i+1) = a1(i) + (a10 + 2._rkind*a11 + 2._rkind*a12 + a13) / 6._rkind
        a2(i+1) = a2(i) + (a20 + 2._rkind*a21 + 2._rkind*a22 + a23) / 6._rkind
        a3(i+1) = a3(i) + (a30 + 2._rkind*a31 + 2._rkind*a32 + a33) / 6._rkind
        a4(i+1) = a4(i) + (a40 + 2._rkind*a41 + 2._rkind*a42 + a43) / 6._rkind
        a5(i+1) = a5(i) + (a50 + 2._rkind*a51 + 2._rkind*a52 + a53) / 6._rkind
        a6(i+1) = a6(i) + (a60 + 2._rkind*a61 + 2._rkind*a62 + a63) / 6._rkind
        a7(i+1) = a7(i) + (a70 + 2._rkind*a71 + 2._rkind*a72 + a73) / 6._rkind
        a8(i+1) = a8(i) + (a80 + 2._rkind*a81 + 2._rkind*a82 + a83) / 6._rkind
        eta(i+1) = eta(i) + deta
        Tbl(i+1) = Te + (T0e - Te) * t(i+1)
        tt = Tbl(i+1)
        if (tt<110.4_rkind) then
          visc(i+1) = .693873D-6*tt
        else
          visc(i+1) = 1.458D-5 * tt**1.5_rkind / ( tt + s_suth )
        end if
        vis = visc(i+1) / visce
      end do
      eq1 = 1._rkind - u(n)
      eq2 = - t(n)
      b1 = a1(n)
      b2 = a1(n)
      b3 = a3(n)
      b4 = a4(n)
      det = b1 * b4 - b2 * b3
      c1 = b4 / det
      c2 = - b2 / det
      c3 = - b3 / det
      c4 = b1 / det
      da = c1 * eq1 + c2 * eq2
      db = c3 * eq1 + c4 * eq2
      s0 = s0 + da
      r0 = r0 + db
      j = j + 1
      if (abs(u(n)-1._rkind)<eps.and.abs(t(n))<eps) exit
    enddo
    dd = .99_rkind
    do i = 0, nm
      if (u(i).ge.dd) exit
    end do
    delta1 = eta(i-1)
    kk = i - 1
    delta2 = delta1 -2*h(kk)
    do i = 0,n
      tcr(i) = tbl(0)/Te + (1._rkind - tbl(0)/Te) * u(i) + (gam - 1._rkind) / 2._rkind *Mach**2. * u(i) * (1._rkind - u(i))
    end do
    do i=0,n
      g1 = .5_rkind*u(i) / tbl(i) * Te
      rho = Te/tbl(i)
      v(i) = eta(i)*g1 - h(i)
    end do
    return
  end subroutine compressible_blasius
  subroutine meanvelocity_bl(retau,rm,trat,s2tinf,vtexp,visc_exp,gam,rfac,ny,y,u,rho,t,visc,redelta,cf,th)
    logical, intent(in) :: visc_exp
    integer, intent(in) :: ny
    real(rkind), intent(in) :: retau,rm,trat,s2tinf,gam,rfac,vtexp
    real(rkind), intent(out) :: redelta,cf,th
    real(rkind), dimension(ny), intent(in) :: y
    real(rkind), dimension(ny), intent(out) :: u,rho,t,visc
    integer :: j,m,itransf,jp99,jj99,jj
    real(rkind) :: gm1h, tr, tw, alf, retau_inc, retau_target, retau_correction, retau_musker
    real(rkind) :: retau_inc_target, retau_inc_correction, up99, d99plus_inc
    real(rkind) :: fuu, du, rhow_over_rho, rhow_over_rhojm, rmu_over_rmuw, rmu_over_rmuwjm
    real(rkind) :: ucij, ucijm, res, dy, ycij, ycijm, u99, d99, d99plus, yy, ui
    real(rkind) :: ue, uuu, uuum, thi
    real(rkind) :: tol_iter_loc
    real(rkind), dimension(ny) :: yplus
    real(rkind), dimension(ny) :: uplus
    real(rkind), dimension(ny) :: ycplus
    real(rkind), dimension(ny) :: ucplus
    real(rkind), dimension(ny) :: yplusold
    real(rkind), dimension(ny) :: ucplusold
    real(rkind), dimension(ny) :: density
    real(rkind), dimension(ny) :: viscosity
    real(rkind), dimension(ny) :: temperature
    real(rkind), dimension(ny) :: uu,uci,uc
    itransf = 1
    gm1h = 0.5_rkind*(gam-1._rkind)
    tr = 1._rkind+gm1h*rfac*rm**2
    tw = trat*tr
    alf = 0.8259_rkind
    tol_iter_loc = 0.000001_rkind
    m = 4
    retau_inc = retau
    retau_target = retau
    retau_correction = 1._rkind
    do
      retau_inc = retau_inc*retau_correction
      retau_musker = retau_inc
      retau_inc_target = retau_inc
      retau_inc_correction = 1._rkind
      do
        retau_musker = retau_musker*retau_inc_correction
        uplus = 0._rkind
        do j=2,ny
          yplus(j) = y(j)*retau_musker
          uplus(j) = velmusker(yplus(j),retau_musker)
        enddo
        up99 = 0.99_rkind*uplus(ny)
        call locateval(uplus,ny,up99,jp99)
        jp99 = min(max(jp99-(m-1)/2,1),ny+1-m)
        call pol_int(uplus(jp99),yplus(jp99),m,up99,d99plus_inc)
        retau_inc_correction = retau_inc_target/d99plus_inc
        if (abs(retau_inc_correction-1._rkind)<tol_iter_loc) exit
      enddo
      ycplus = y*retau_musker
      yplus = ycplus
      yplusold = yplus
      do
        uplus = 0._rkind
        do j=2,ny
          uplus(j) = velmusker(yplus(j),retau_musker)
        enddo
        ucplus = uplus
        do
          ucplusold = ucplus
          do j=1,ny
            uu(j) = ucplus(j)/ucplus(ny)
            fuu = uu(j)
            temperature(j) = tw+(tr-tw)*fuu+(1._rkind-tr)*uu(j)**2
            density(j) = 1._rkind/temperature(j)
            if (visc_exp) then
              viscosity(j) = temperature(j)**vtexp
            else
              viscosity(j) = temperature(j)**(1.5_rkind)*(1._rkind+s2tinf)/(temperature(j)+s2tinf)
            endif
          enddo
          do j=2,ny
            du = uplus(j)-uplus(j-1)
            rhow_over_rho = density(1)/density(j)
            rhow_over_rhojm = density(1)/density(j-1)
            rmu_over_rmuw = viscosity(j)/viscosity(1)
            rmu_over_rmuwjm = viscosity(j-1)/viscosity(1)
            if (itransf==1) then
              ucij = sqrt(rhow_over_rho)
              ucijm = sqrt(rhow_over_rhojm)
            else
              ucij = sqrt(rhow_over_rho)*sqrt(rmu_over_rmuw)
              ucijm = sqrt(rhow_over_rhojm)*sqrt(rmu_over_rmuwjm)
            endif
            ucplus(j) = ucplus(j-1)+0.5_rkind*(ucij+ucijm)*du
          enddo
          res = 0._rkind
          do j=1,ny
            res = res+abs(ucplus(j)-ucplusold(j))
          enddo
          res = res/ny
          if (res<tol_iter_loc) exit
        enddo
        yplus = 0._rkind
        do j=2,ny
          dy = ycplus(j)-ycplus(j-1)
          rhow_over_rho = density(1)/density(j)
          rhow_over_rhojm = density(1)/density(j-1)
          rmu_over_rmuw = viscosity(j)/viscosity(1)
          rmu_over_rmuwjm = viscosity(j-1)/viscosity(1)
          if (itransf==1) then
            ycij = 1._rkind
            ycijm = 1._rkind
          else
            ycij = sqrt(rhow_over_rho)*sqrt(rmu_over_rmuw**3)
            ycij = 1._rkind/ycij
            ycijm = sqrt(rhow_over_rhojm)*sqrt(rmu_over_rmuwjm**3)
            ycijm = 1._rkind/ycijm
          endif
          yplus(j) = yplus(j-1)+0.5_rkind*(ycij+ycijm)*dy
        enddo
        res = 0._rkind
        do j=1,ny
          res = abs(yplus(j)-yplusold(j))
        enddo
        res = res/ny
        if (res<tol_iter_loc) exit
        yplusold = yplus
      enddo
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
    enddo
    do j=1,ny
      uc(j) = ucplus(j)/ucplus(ny)
    enddo
    u99 = 0.99_rkind
    call locateval(uc,ny,u99,jj99)
    jj99 = min(max(jj99-(m-1)/2,1),ny+1-m)
    call pol_int(uc(jj99),y(jj99),m,u99,d99)
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
      fuu = uu(j)
      temperature(j) = tw+(tr-tw)*fuu+(1.-tr)*uu(j)**2
      density(j) = 1._rkind/temperature(j)
      if (visc_exp) then
        viscosity(j) = temperature(j)**vtexp
      else
        viscosity(j) = temperature(j)**(1.5_rkind)*(1._rkind+s2tinf)/(temperature(j)+s2tinf)
      endif
    enddo
    retau_musker = retau_musker*d99
    redelta = retau_musker*ucplus(ny)/density(1)*viscosity(1)
    cf = 2._rkind*density(1)/ucplus(ny)**2
    th = 0._rkind
    ue = uc(ny)
    do j=2,ny
      uuu = uc(j)
      uuum = uc(j-1)
      dy = y(j)-y(j-1)
      thi = 0.5_rkind*(density(j)*uuu*(1.-uuu)+density(j-1)*uuum*(1.-uuum))
      th = th+thi*dy
    enddo
    do j=1,ny
      u(j) = uc(j)
      rho(j) = density(j)
      t(j) = temperature(j)
      visc(j) = viscosity(j)
    enddo
    return
  end subroutine meanvelocity_bl
  subroutine hasan_meanprofile(n,d99,mach,theta,retau,retheta,redelta,pr,spr,rfac,gam,vtexp,&
  &visc_exp,s2tinf,ye,ue,te,rhoe,mue,th,cf,ch,mtau,deltav,imode,icompute)
    implicit none

    integer, intent(in) :: n,imode,icompute
    real(rkind), intent(in) :: d99,mach,theta,pr,spr,rfac,gam,s2tinf,vtexp
    real(rkind), intent(inout) :: retau,retheta
    real(rkind), intent(inout) :: cf,ch,th,mtau,deltav,redelta
    real(rkind), dimension(n+1), intent(in) :: ye
    real(rkind), dimension(n+1), intent(out) :: ue,te,rhoe,mue
    logical, intent(in) :: visc_exp

    integer, parameter :: ny = 10000
    integer, parameter :: nyext = 2*ny
    integer :: j,niter,nitermax,l
    real(rkind) :: pi
    real(rkind) :: gm1h
    real(rkind) :: rethetaold,retauold
    real(rkind) :: tr,tw,picole,apl,vkc,z1
    real(rkind) :: absdiff,tol
    real(rkind) :: uinf,uinf_star,uinf_plus,ufac
    real(rkind) :: dy,dudy
    real(rkind) :: thint,thint_j,thint_jm
    real(rkind) :: uint,uint_j,uint_jm,uint_inn_j,uint_inn_jm,uint_out_j,uint_out_jm
    real(rkind), dimension(ny+1) :: y,u,t,rho,mu
    real(rkind), dimension(ny+1) :: ypl,upl
    real(rkind), dimension(ny+1) :: yst,d,mut,wake
    real(rkind), dimension(nyext+1) :: yext,uext,text,rhoext,muext
    integer :: jj,jjj,m
    real(rkind) :: yy
    pi = 4._rkind*atan(1._rkind)
    gm1h = 0.5_rkind*(gam-1._rkind)
    tr = 1._rkind+gm1h*rfac*mach**2
    tw = 1._rkind+theta*(tr-1._rkind)
    apl = 17._rkind
    vkc = 0.41_rkind
    nitermax = 10000
    tol = 0.00001_rkind
    dy = d99/ny
    do j=1,ny+1
      y(j) = (j-1)*dy
    enddo
    do j=1,nyext+1
      yext(j) = (j-1)*dy
    enddo
    ufac = 0.99_rkind
    do j=1,ny+1
      u(j) = y(j)/d99
      u(j) = min(1._rkind,u(j))
    enddo
    uinf = u(ny+1)/ufac
    u = u/uinf
    niter = 0
    if (imode==0) then
      retheta = 1000._rkind
    else
      retau = 200._rkind
    endif
    mtau = 0.1_rkind
    absdiff = huge(1._rkind)
    do while (absdiff > tol)
      niter = niter+1
      if (niter>nitermax) exit
      if (imode==0) then
        rethetaold = retheta
      else
        retauold = retau
      endif
      do j=1,ny+1
        t(j) = tw+(tr-tw)*(spr*u(j)+(1._rkind-spr)*u(j)**2)+(1._rkind-tr)*u(j)**2
        rho(j) = 1._rkind/t(j)
        if (visc_exp) then
          mu(j) = t(j)**vtexp
        else
          mu(j) = t(j)**1.5_rkind*(1._rkind+s2tinf)/(t(j)+s2tinf)
        endif
        ypl(j) = (y(j)/d99)*retau
        yst(j) = ypl(j)*sqrt(rho(j)/rho(1))*mu(1)/mu(j)
        d(j) = (1._rkind-exp(-yst(j)/(apl+19.3_rkind*mtau)))**2
        mut(j) = mu(j)/mu(1)*vkc*yst(j)*d(j)
        z1 = retheta/425._rkind-1._rkind
        picole = 0.69_rkind*(1._rkind-exp(-0.243_rkind*sqrt(z1)-0.15_rkind*z1))
        wake(j) = picole/vkc*pi*sin(pi*(y(j)/d99))
      enddo
      upl = 0._rkind
      do j=2,ny+1
        uint_inn_j = 1._rkind/(mu(j)/mu(1)+mut(j))
        uint_inn_jm = 1._rkind/(mu(j-1)/mu(1)+mut(j-1))
        uint_out_j = sqrt(rho(1)/rho(j))/retau*wake(j)
        uint_out_jm = sqrt(rho(1)/rho(j-1))/retau*wake(j-1)
        uint_j = uint_inn_j +uint_out_j
        uint_jm = uint_inn_jm+uint_out_jm
        uint = 0.5_rkind*(uint_j+uint_jm)
        upl(j) = upl(j-1)+uint*(ypl(j)-ypl(j-1))
      enddo
      uinf_plus = upl(ny+1)/ufac
      u = upl/uinf_plus
      th = 0._rkind
      do j=2,ny+1
        thint_j = rho(j)*u(j)*(1._rkind-u(j))
        thint_jm = rho(j-1)*u(j-1)*(1._rkind-u(j-1))
        thint = 0.5_rkind*(thint_j+thint_jm)
        th = th+thint*(y(j)-y(j-1))
      enddo
      cf = 2._rkind/uinf_plus**2*rho(1)
      ch = 0.5_rkind*cf*spr/pr
      mtau = mach*sqrt(0.5_rkind*cf)
      if (imode==0) then
        retheta = retau*th/d99*mu(1)/rho(1)*uinf_plus
        absdiff = abs(retheta-rethetaold)
      else
        retau = retheta*d99/th/mu(1)*rho(1)/uinf_plus
        absdiff = abs(retau-retauold)
      endif
      redelta = retheta*d99/th
    enddo
    deltav = d99/retau
    if (niter>nitermax) write(*,*) 'Maximum number of iterations in Hasan!',absdiff
    do j=1,ny+1
      uext(j) = u(j)
      text(j) = t(j)
    enddo
    dudy = (uext(ny+1)-uext(ny))/(yext(ny+1)-yext(ny))
    do j=ny+2,nyext+1
      uext(j) = uext(ny+1)+dudy*(yext(j)-yext(ny+1))
      uext(j) = min(uext(j),1._rkind)
      text(j) = tw+(tr-tw)*(spr*uext(j)+(1._rkind-spr)*uext(j)**2)+(1._rkind-tr)*uext(j)**2
    enddo
    do j=1,nyext+1
      rhoext(j) = 1._rkind/text(j)
      if (visc_exp) then
        muext(j) = text(j)**vtexp
      else
        muext(j) = text(j)**1.5_rkind*(1._rkind+s2tinf)/(text(j)+s2tinf)
      endif
    enddo
    ue = 1._rkind
    rhoe = 1._rkind
    te = 1._rkind
    mue = 1._rkind
    if (icompute/=0) then
      ue(1) = uext(1)
      rhoe(1) = rhoext(1)
      te(1) = text(1)
      mue(1) = muext(1)
      do j=1,n
        yy = ye(j)
        if (yy<yext(nyext+1)) then
          call locateval(yext,nyext+1,yy,jj)
          m = 2
          jjj = min(max(jj-(m-1)/2,1),nyext+1+1-m)
          call pol_int(yext(jjj),uext(jjj),m,yy,ue(j))
          call pol_int(yext(jjj),rhoext(jjj),m,yy,rhoe(j))
          call pol_int(yext(jjj),text(jjj),m,yy,te(j))
          call pol_int(yext(jjj),muext(jjj),m,yy,mue(j))
        endif
      enddo
    endif
    return
  endsubroutine hasan_meanprofile
  subroutine init_wind_tunnel(self)
    class(equation_singleideal_object), intent(inout) :: self

    integer :: i,j,k,l
    allocate(self%wmean(1-self%grid%ng:self%field%nx+self%grid%ng+1, 1:self%field%ny, 4))

    associate(nx => self%field%nx, ny => self%field%ny, nz => self%field%nz, ng => self%grid%ng,&
    & nv => self%nv, wmean => self%wmean, w => self%field%w, winf => self%winf)

      wmean = 0._rkind
      do j=1,ny
        do i=1-ng,nx+ng+1
          wmean(i,j,1) = winf(1)
          wmean(i,j,2) = winf(2)
          wmean(i,j,3) = winf(3)
          wmean(i,j,4) = winf(4)
        enddo
      enddo

      do k=1,nz
        do j=1,ny
          do i=1,nx
            do l=1,nv
              w(i,j,k,l) = winf(l)
            enddo
          enddo
        enddo
      enddo

    endassociate
  endsubroutine init_wind_tunnel

  subroutine init_cv(self)
    class(equation_singleideal_object), intent(inout) :: self

    integer :: i,j,k
    real(rkind) :: rmv,xx,yy,rr2,gm1h,gm1i,exparg,uvtx,vvtx
    real(rkind) :: rho,pp,tt,uu,vv,ww,rhouu,rhovv,rhoww,ee
    allocate(self%wmean(1-self%grid%ng:self%field%nx+self%grid%ng+1, 1:self%field%ny, 4))

    associate(nx => self%field%nx, ny => self%field%ny, nz => self%field%nz, ng => self%grid%ng,&
    & nv => self%nv, wmean => self%wmean, w => self%field%w, winf => self%winf, w_aux => self%w_aux,&
    & Mach => self%Mach, rho0 => self%rho0, u0 => self%u0, gam => self%gam, y => self%field%y,&
    & z => self%field%z, x => self%field%x, indx_cp_l => self%indx_cp_l, indx_cp_r => self%indx_cp_r,&
    & cv_coeff => self%cv_coeff, calorically_perfect => self%calorically_perfect, rgas0 => self%rgas0,&
    & t0 => self%t0, l0 => self%l0, p0 => self%p0)

      wmean = 0._rkind
      do j=1,ny
        do i=1-ng,nx+ng+1
          wmean(i,j,1) = winf(1)
          wmean(i,j,2) = winf(2)
          wmean(i,j,3) = winf(3)
          wmean(i,j,4) = winf(4)
        enddo
      enddo

      rmv = 0.1_rkind
      gm1h = 0.5_rkind*(gam-1._rkind)
      gm1i = 1._rkind/(gam-1._rkind)
      do k=1,nz
        do j=1,ny
          do i=1,nx
            xx = x(i)/l0-20._rkind
            yy = y(j)/l0-10._rkind
            rr2 = xx*xx+yy*yy
            exparg = 0.5_rkind*(1._rkind-rr2)
            rho = rho0*(1._rkind-gm1h*rmv**2*exp(1._rkind-rr2))**gm1i
            pp = p0*(1._rkind-gm1h*rmv**2*exp(1._rkind-rr2))**gm1i
            uvtx = -rmv/Mach*yy*exp(exparg)
            vvtx = rmv/Mach*xx*exp(exparg)
            uu = 1._rkind+uvtx
            vv = vvtx
            uu = uu*u0
            vv = vv*u0
            ww = 0._rkind
            rhouu = rho*uu
            rhovv = rho*vv
            rhoww = rho*ww
            w(i,j,k,1) = rho
            w(i,j,k,2) = rhouu
            w(i,j,k,3) = rhovv
            w(i,j,k,4) = rhoww
            tt = pp/rgas0/rho
            w_aux(i,j,k,6) = tt
            ee = get_e_from_temperature(tt, t0, indx_cp_l, indx_cp_r, cv_coeff, calorically_perfect)
            w(i,j,k,5) = rho*ee + 0.5_rkind*(rhouu**2+rhovv**2+rhoww**2)/rho
          enddo
        enddo
      enddo

    endassociate
  endsubroutine init_cv

  subroutine init_channel(self)
    class(equation_singleideal_object), intent(inout) :: self
    integer :: i, j, k
    real(rkind), dimension(3) :: rr
    real(rkind) :: u0_02, u0_05, rho, uu, vv, ww
    real(rkind) :: ufluc, vfluc, wfluc
    real(rkind) :: rhouu, rhovv, rhoww
    real(rkind) :: tt, ee
    integer :: nroll

    allocate(self%wmean(1-self%grid%ng:self%field%nx+self%grid%ng+1, 1:self%field%ny, 4))

    associate(nx => self%field%nx, ny => self%field%ny, nz => self%field%nz, ng => self%grid%ng,&
    & y => self%field%y, z => self%field%z, rho0 => self%rho0, u0 => self%u0, p0 => self%p0,&
    & gm => self%gm, t0 => self%t0, wmean => self%wmean, w => self%field%w, w_aux => self%w_aux,&
    & l0 => self%l0, indx_cp_l => self%indx_cp_l, indx_cp_r => self%indx_cp_r, cv_coeff => self%cv_coeff,&
    & volchan => self%volchan, nxmax => self%grid%nxmax, nzmax => self%grid%nzmax, yn => self%grid%yn,&
    & lz => self%grid%domain_size(3), calorically_perfect => self%calorically_perfect,&
    & rgas0 => self%rgas0)

      volchan = yn(ny+1)-yn(1)
      volchan = volchan*nxmax
      volchan = volchan*nzmax

      wmean = 0._rkind
      do j=1,ny
        do i=1-ng,nx+ng+1
          wmean(i,j,1) = rho0
          wmean(i,j,2) = 1.5_rkind*rho0*u0*(1._rkind-(y(j)/l0)**2)
          wmean(i,j,3) = 0._rkind
          wmean(i,j,4) = 0._rkind
        enddo
      enddo

      u0_02 = 0.02_rkind*u0
      if (self%rand_type==0) u0_02 = 0._rkind
      u0_05 = 0.05_rkind*u0
      nroll = int(lz/yn(ny+1))

      do k=1,nz
        do j=1,ny
          do i=1,nx

            rho = wmean(i,j,1)
            uu = wmean(i,j,2)/rho
            vv = wmean(i,j,3)/rho
            ww = wmean(i,j,4)/rho

            call get_crandom_f(rr)
            rr = rr-0.5_rkind

            ufluc = u0_02*rr(1)
            vfluc = u0_02*rr(2)
            wfluc = u0_02*rr(3)
            vfluc = vfluc+u0_05*sin(0.5_rkind*pi*y(j)/l0)*cos(2*pi*z(k)/lz*nroll)
            wfluc = wfluc+u0_05*sin(0.5_rkind*pi*y(j)/l0)*sin(2*pi*z(k)/lz*nroll)
            uu = uu+ufluc
            vv = vv+vfluc
            ww = ww+wfluc
            rhouu = rho*uu
            rhovv = rho*vv
            rhoww = rho*ww
            w(i,j,k,1) = rho
            w(i,j,k,2) = rhouu
            w(i,j,k,3) = rhovv
            w(i,j,k,4) = rhoww
            tt = p0/rgas0/rho
            w_aux(i,j,k,6) = tt
            ee = get_e_from_temperature(tt, t0, indx_cp_l, indx_cp_r, cv_coeff, calorically_perfect)
            w(i,j,k,5) = rho*ee + 0.5_rkind*(rhouu**2+rhovv**2+rhoww**2)/rho
          enddo
        enddo
      enddo

    endassociate
  endsubroutine init_channel

  subroutine init_airfoil(self)
    class(equation_singleideal_object), intent(inout) :: self
    real(rkind) :: errmin, error, xtrip, rtemp
    real(rkind) :: ds, ds_tot
    integer :: i, j, rank_tr, rank_ts, itr_l, its_l, rank_tr_x, rank_ts_x, itemp, ii, l
    character(1) :: var

    call self%cfg%get("curvi","xtr1",rtemp) ; self%xtr1 = rtemp
    call self%cfg%get("curvi","xtr2",rtemp) ; self%xtr2 = rtemp
    call self%cfg%get("curvi","xtw1",rtemp) ; self%xtw1 = rtemp
    call self%cfg%get("curvi","xtw2",rtemp) ; self%xtw2 = rtemp
    call self%cfg%get("curvi","a_tr",rtemp) ; self%a_tr = rtemp
    call self%cfg%get("curvi","a_tw",rtemp) ; self%a_tw = rtemp
    if (self%cfg%has_key("curvi","del0_tr")) then
      call self%cfg%get("curvi","del0_tr",rtemp) ; self%del0_tr = rtemp
    else
      self%del0_tr = 0.001_rkind
    endif
    if (self%cfg%has_key("curvi","v_bs")) then
      call self%cfg%get("curvi","v_bs",rtemp) ; self%v_bs = rtemp
    else
      self%v_bs = 0._rkind
    endif
    call self%cfg%get("curvi","kx_tw",rtemp) ; self%kx_tw = rtemp
    call self%cfg%get("curvi","om_tw",rtemp) ; self%om_tw = rtemp
    call self%cfg%get("curvi","thic",rtemp) ; self%thic = rtemp
    call self%cfg%get("curvi","teshk",rtemp) ; self%teshk = rtemp

    associate(ite => self%grid%ite, ile => self%grid%ile, itu => self%grid%itu, itr => self%itr,&
    & its => self%its, jtr => self%jtr, jts => self%jts, nx => self%field%nx, ny => self%field%ny,&
    & ncoords => self%field%ncoords, ierr => self%mpi_err, itr1 => self%itr1, itr2 => self%itr2,&
    & its1 => self%its1, its2 => self%its2, nxmax => self%grid%nxmax, xwall => self%grid%xwall,&
    & ywall => self%grid%ywall, x0tr => self%x0tr, y0tr => self%y0tr, x0ts => self%x0ts,&
    & y0ts => self%y0ts, xc2 => self%field%xc2, yc2 => self%field%yc2, mp_cart => self%field%mp_cart,&
    & masterproc => self%masterproc, xtr1 => self%xtr1, xtr2 => self%xtr2, xtrip => self%xtrip,&
    & xtw1 => self%xtw1, xtw2 => self%xtw2, a_tr => self%a_tr, a_tw => self%a_tw, v_bs => self%v_bs,&
    &lamz_old => self%lamz_old , lamz1_old => self%lamz1_old , lams_old => self%lams_old ,&
    & lams1_old => self%lams1_old , phiz_old => self%phiz_old , phiz1_old => self%phiz1_old ,&
    & phis_old => self%phis_old , phis1_old => self%phis1_old, is_old => self%is_old)

      xtrip = xtr1+0.5_rkind*abs(xtr2-xtr1)

      errmin = 100.
      do i=ite,ile
        error = abs(xwall(i)-xtrip)
        if (error < errmin) then
          errmin = error
          itr = i
        endif
      enddo
      errmin = 100.
      do i=ile,itu
        error = abs(xwall(i)-xtrip)
        if (error < errmin) then
          errmin = error
          its = i
        endif
      enddo

      if (self%cfg%has_key("curvi","itr")) then
        call self%cfg%get("curvi","itr",itr)
        if(self%masterproc) then
          print*,'Tripping index pressure found: ',itr
        endif
      endif
      if (self%cfg%has_key("curvi","its")) then
        call self%cfg%get("curvi","its",its)
        if(self%masterproc) then
          print*,'Tripping index suction found: ',its
        endif
      endif
      jtr = 1
      if (self%cfg%has_key("curvi","jtr")) then
        call self%cfg%get("curvi","jtr",jtr)
        if(self%masterproc) then
          print*,'Tripping j-index pressure found: ',jtr
        endif
      endif
      jts = 1
      if (self%cfg%has_key("curvi","jts")) then
        call self%cfg%get("curvi","jts",jts)
        if(self%masterproc) then
          print*,'Tripping j-index suction found: ',jts
        endif
      endif

      rank_tr_x = int((itr-1)/nx)
      rank_ts_x = int((its-1)/nx)
      if(ncoords(1) == rank_tr_x) then
        itr_l = itr - rank_tr_x*nx
        x0tr = xc2(itr_l,jtr)
        y0tr = yc2(itr_l,jtr)
      endif
      if(ncoords(1) == rank_ts_x) then
        its_l = its - rank_ts_x*nx
        x0ts = xc2(its_l,jts)
        y0ts = yc2(its_l,jts)
        print*,'local x0ts, y0ts: ',x0ts, y0ts
      endif
      call MPI_Cart_rank(self%field%mp_cart, [rank_tr_x, 0, 0], rank_tr, ierr)
      call MPI_Cart_rank(self%field%mp_cart, [rank_ts_x, 0, 0], rank_ts, ierr)

      call MPI_Bcast(x0tr, 1, mpi_prec, rank_tr, self%field%mp_cart, ierr)
      call MPI_Bcast(y0tr, 1, mpi_prec, rank_tr, self%field%mp_cart, ierr)
      call MPI_Bcast(x0ts, 1, mpi_prec, rank_ts, self%field%mp_cart, ierr)
      call MPI_Bcast(y0ts, 1, mpi_prec, rank_ts, self%field%mp_cart, ierr)
      if(self%masterproc) print*,'global: x0ts, y0ts: ',x0ts, y0ts

      ds_tot = 0.
      do i=itr-1,1,-1
        ds = ((xwall(i)-xwall(i+1))**2+(ywall(i)-ywall(i+1))**2)**0.5
        ds_tot = ds_tot + ds
        if(ds_tot > 0.5_rkind*abs(xtr2-xtr1)) then
          itr1 = i
          exit
        endif
      enddo
      ds_tot = 0.
      do i=itr+1,nxmax
        ds = ((xwall(i)-xwall(i-1))**2+(ywall(i)-ywall(i-1))**2)**0.5
        ds_tot = ds_tot + ds
        if(ds_tot > 0.5_rkind*abs(xtr2-xtr1)) then
          itr2 = i
          exit
        endif
      enddo
      ds_tot = 0.
      do i=its-1,1,-1
        ds = ((xwall(i)-xwall(i+1))**2+(ywall(i)-ywall(i+1))**2)**0.5
        ds_tot = ds_tot + ds
        if(ds_tot > 0.5_rkind*abs(xtr2-xtr1)) then
          its1 = i
          exit
        endif
      enddo
      ds_tot = 0.
      do i=its+1,nxmax
        ds = ((xwall(i)-xwall(i-1))**2+(ywall(i)-ywall(i-1))**2)**0.5
        ds_tot = ds_tot + ds
        if(ds_tot > 0.5_rkind*abs(xtr2-xtr1)) then
          its2 = i
          exit
        endif
      enddo

      if (masterproc) then
        if (a_tr > 0.) then
          print*, 'Tripping center on pressure side', x0tr,y0tr
          print*, 'Tripping center on suction side', x0ts,y0ts
          print*, 'Global tripping center nodes pressure ', itr,jtr
          print*, 'Global tripping center nodes suction ', its,jts
          print*, 'Global start-end pressure tripping ', itr1,itr2
          print*, 'Global start-end suction tripping ', its1,its2
        else
          print*, 'No tripping'
        endif
        if (a_tw > 0. .and. abs(v_bs) > 0.) then
          print*, 'Error: traveling waves and blowing/suction cannot be both active'
          call mpi_abort(mpi_comm_world,99,self%mpi_err)
        endif
        if (a_tw > 0.) then
          print*, 'Traveling waves active with amplitude =', a_tw
        elseif (v_bs > 0.) then
          print*, 'Blowing active with velocity =', v_bs
        elseif (v_bs < 0.) then
          print*, 'Suction active with velocity =', v_bs
        endif
        if (a_tw > 0. .or. abs(v_bs) > 0.) then
          print*, 'Start - end of actuated region ', xtw1,xtw2
        else
          print*, 'No actuation'
        endif
      endif

      lamz_old = 0._rkind
      lamz1_old = 0._rkind
      lams_old = 0._rkind
      lams1_old = 0._rkind
      phiz_old = 0._rkind
      phiz1_old = 0._rkind
      phis_old = 0._rkind
      phis1_old = 0._rkind
      is_old = 0
    endassociate

    if(self%restart_type == 0) then
      open(unit=30,file='airfoil_forces.dat')
      write(30,*) '# time,lift,drag,pres,fric'
      close(30)
    endif

  endsubroutine init_airfoil

  subroutine init_channel_c2(self)
    class(equation_singleideal_object), intent(inout) :: self
    integer :: i, j, k
    real(rkind), dimension(3) :: rr
    real(rkind) :: u0_02, u0_05, rho, uu, vv, ww
    real(rkind) :: ufluc, vfluc, wfluc
    real(rkind) :: rhouu, rhovv, rhoww
    real(rkind) :: tt, ee, eta
    real(rkind) :: aa, bb, alfa, beta, ulam, ulam_mean, y_temp, ulam_mean_integ, area
    integer :: nroll

    allocate(self%wmean(1-self%grid%ng:self%field%nx+self%grid%ng+1, 1:self%field%ny, 4))

    associate(nx => self%field%nx, ny => self%field%ny, nz => self%field%nz, ng => self%grid%ng,&
    & y => self%field%y, z => self%field%z, rho0 => self%rho0, u0 => self%u0, p0 => self%p0,&
    & gm => self%gm, t0 => self%t0, wmean => self%wmean, w => self%field%w, w_aux => self%w_aux,&
    & l0 => self%l0, indx_cp_l => self%indx_cp_l, indx_cp_r => self%indx_cp_r, cv_coeff => self%cv_coeff,&
    & yn => self%grid%yn, dxdcsinc2 => self%field%dxdcsinc2, dydcsinc2 => self%field%dydcsinc2,&
    & lz => self%grid%domain_size(3), rgas0 => self%rgas0, calorically_perfect => self%calorically_perfec&
    &t, R_curv => self%grid%R_curv)

      aa = R_curv - 1._rkind
      bb = R_curv + 1._rkind
      alfa = (bb**2*log(bb)-aa**2*log(aa))/(bb**2-aa**2)
      beta = (aa**2*bb**2) * log(bb/aa) /(bb**2-aa**2)
      ulam_mean = (0.25_rkind*(bb**2-aa**2)-(aa**2*bb**2)*(log(bb/aa))**2/(bb**2-aa**2))/(bb-aa)

      wmean = 0._rkind
      ulam_mean_integ = 0._rkind
      area = 0._rkind
      do i=1,nx
        do j=1,ny
          wmean(i,j,1) = rho0



          y_temp = R_curv + 0.5_rkind*(yn(j)+yn(j+1))
          ulam = (alfa*y_temp-beta/y_temp-y_temp*log(y_temp))

          ulam_mean_integ = ulam_mean_integ + ulam/self%field%jac(i,j)
          area = area + 1._rkind/self%field%jac(i,j)
        enddo
      enddo
      call mpi_allreduce(MPI_IN_PLACE,ulam_mean_integ,1,mpi_prec,mpi_sum,self%field%mp_cartx,self%mpi_err)
      call mpi_allreduce(MPI_IN_PLACE,area,1,mpi_prec,mpi_sum,self%field%mp_cartx,self%mpi_err)
      ulam_mean_integ = ulam_mean_integ/area

      do i=1-ng,nx+ng
        do j=1,ny
          y_temp = R_curv + 0.5_rkind*(yn(j)+yn(j+1))
          ulam = (alfa*y_temp-beta/y_temp-y_temp*log(y_temp))
          ulam = ulam/ulam_mean_integ
          wmean(i,j,2) = rho0*u0*ulam*dxdcsinc2(i,j)
          wmean(i,j,3) = rho0*u0*ulam*dydcsinc2(i,j)
          wmean(i,j,4) = 0._rkind
        enddo
      enddo

      u0_05 = 0.5_rkind*u0

      nroll = int(lz/yn(ny+1))
      do k=1,nz
        do j=1,ny
          do i=1,nx
            rho = wmean(i,j,1)
            uu = wmean(i,j,2)/rho
            vv = wmean(i,j,3)/rho
            ww = wmean(i,j,4)/rho
            uu = uu + u0_05*sin(0.5_rkind*pi*y(j)/l0)*sin(2*pi*z(k)/lz*nroll)
            vv = vv + u0_05*sin(0.5_rkind*pi*y(j)/l0)*cos(2*pi*z(k)/lz*nroll)
            ww = ww + u0_05*sin(0.5_rkind*pi*y(j)/l0)*sin(2*pi*z(k)/lz*nroll)
            rhouu = rho*uu
            rhovv = rho*vv
            rhoww = rho*ww
            w(i,j,k,1) = rho
            w(i,j,k,2) = rhouu
            w(i,j,k,3) = rhovv
            w(i,j,k,4) = rhoww
            tt = p0/rgas0/rho
            w_aux(i,j,k,6) = tt
            ee = get_e_from_temperature(tt, t0, indx_cp_l, indx_cp_r, cv_coeff, calorically_perfect)
            w(i,j,k,5) = rho*ee + 0.5_rkind*(rhouu**2+rhovv**2+rhoww**2)/rho
          enddo
        enddo
      enddo

    endassociate
  endsubroutine init_channel_c2

  subroutine bc_preproc(self)
    class(equation_singleideal_object), intent(inout) :: self
    integer :: ilat, offset

    self%force_zero_flux = [0,0,0,0,0,0]
    self%eul_imin = 1
    self%eul_imax = self%field%nx
    self%eul_jmin = 1
    self%eul_jmax = self%field%ny
    self%eul_kmin = 1
    self%eul_kmax = self%field%nz
    do ilat=1,6
      select case(self%bctags(ilat))
      case(0)
      case(1)
      case(2)
      case(4)
      case(5)
      case(6)
        if (self%grid%is_y_staggered) then
          self%force_zero_flux(ilat) = 1
          self%bctags_nr(ilat) = 0
        endif
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

  function get_gamloc(tt,t0,indx_cp_l,indx_cp_r,cp_coeff,calorically_perfect,rgas0)
    real(rkind) :: get_gamloc
    integer :: indx_cp_l,indx_cp_r,calorically_perfect
    real(rkind) :: tt, rgas0, t0
    real(rkind), dimension(indx_cp_l:indx_cp_r+1) :: cp_coeff
    real(rkind) :: cploc, gamloc
    integer :: l

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

  endfunction get_gamloc

  subroutine runge_kutta_initialize(self)
    class(equation_singleideal_object), intent(inout) :: self

    select case(self%rk_type)
    case(RK_WRAY)
      self%nrk = 3
      allocate(self%rhork(self%nrk),self%gamrk(self%nrk),self%alprk(self%nrk))
      self%rhork(:) = [0._rkind, -17._rkind/60._rkind , -5._rkind /12._rkind]
      self%gamrk(:) = [8._rkind /15._rkind, 5._rkind /12._rkind, 3._rkind /4._rkind]
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

  endsubroutine runge_kutta_initialize
  subroutine read_input(self, filename)
    class(equation_singleideal_object), intent(inout) :: self
    character(*) , intent(in) :: filename

    self%cfg=parse_cfg(filename)

  endsubroutine read_input
  subroutine write_slice_p3d(self, plot3dgrid, plot3dfield)
    class(equation_singleideal_object), intent(inout) :: self
    integer :: i,j,k,l,m,ll
    integer :: mpi_io_file
    integer :: filetype
    integer,dimension(3) :: sizes
    integer,dimension(3) :: subsizes
    integer,dimension(3) :: starts
    integer :: size_real, size_integer
    integer :: ntot
    integer (kind=mpi_offset_kind) :: offset
    integer, dimension(mpi_status_size) :: istatus
    integer :: ii, jj, kk, i_slice, j_slice, k_slice, isize, jsize, ksize, istat
    character(len=14) :: i_folder, j_folder, k_folder
    logical, optional :: plot3dgrid, plot3dfield
    logical :: plot3dgrid_, plot3dfield_
    character(6) :: nastore

    plot3dfield_ = .true. ; if(present(plot3dfield)) plot3dfield_ = plot3dfield
    plot3dgrid_ = .false. ; if(present(plot3dgrid)) plot3dgrid_ = plot3dgrid

    size_real = storage_size(1._rkind)/8
    size_integer = storage_size(1)/8


    associate(masterproc => self%masterproc, nx => self%field%nx, ny => self%field%ny,&
    & nz => self%field%nz, ndim => self%grid%ndim, nxmax => self%grid%nxmax, nymax => self%grid%nymax,&
    & nzmax => self%grid%nzmax, mp_cart => self%field%mp_cart, mp_cartx => self%field%mp_cartx,&
    & mp_cartz => self%field%mp_cartz, nblocks => self%field%nblocks, ncoords => self%field%ncoords,&
    & iermpi => self%mpi_err, time => self%time, itslice_p3d => self%itslice_p3d, flow_init => self%flow_&
    &init, num_aux_slice => self%num_aux_slice, x => self%field%xc2, y => self%field%yc2,&
    & z => self%field%z, grid => self%field%grid3d, w => self%field%w, w_aux => self%w_aux)

      p3d_grid: if(plot3dgrid_) then

        if (allocated(self%jslice_p3d)) then
          jsize = size(self%jslice_p3d)
          do jj=1,jsize
            j_slice = self%jslice_p3d(jj)
            j_folder = "SLICEXZ_??????"
            write(j_folder(9:14),'(I6.6)') j_slice
            if (masterproc) then
              istat = create_folder(j_folder)
              open(unit=123, file=trim(j_folder)//"/grid2d.xyz", form="unformatted", access="stream")
              write(123) nxmax,1,nzmax
              flush(123)
              close(123)
            endif
            call mpi_barrier(mpi_comm_world,iermpi)
            sizes(1) = nblocks(1)*nx
            sizes(2) = 1
            sizes(3) = nblocks(3)*nz
            subsizes(1) = nx
            subsizes(2) = 1
            subsizes(3) = nz
            starts(1) = 0 + ncoords(1)*subsizes(1)
            starts(2) = 0
            starts(3) = 0 + ncoords(3)*subsizes(3)
            ntot = nx*nz
            call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,mpi_prec,filetype,iermpi)
            call MPI_TYPE_COMMIT(filetype,iermpi)
            call MPI_FILE_OPEN(mp_cart,trim(j_folder)//"/grid2d.xyz",MPI_MODE_RDWR,MPI_INFO_NULL,mpi_io_file,iermpi)
            offset = 3*size_integer
            do l=1,3
              call MPI_FILE_SET_VIEW(mpi_io_file,offset,mpi_prec,filetype,"native",MPI_INFO_NULL,iermpi)
              call MPI_FILE_WRITE_ALL(mpi_io_file,grid(l,1:nx,j_slice:j_slice,1:nz),ntot,mpi_prec,istatus,iermpi)
              call MPI_TYPE_SIZE(mpi_prec,size_real,iermpi)
              do m=1,nblocks(1)*nblocks(3)
                offset = offset+size_real*ntot
              enddo
            enddo
            call MPI_FILE_CLOSE(mpi_io_file,iermpi)
            call MPI_TYPE_FREE(filetype,iermpi)
          enddo
        endif

        if (allocated(self%kslice_p3d)) then
          ksize = size(self%kslice_p3d)
          do kk=1,ksize
            k_slice = self%kslice_p3d(kk)
            k_folder = "SLICEXY_??????"
            write(k_folder(9:14),'(I6.6)') k_slice + ncoords(3)*nz
            if (ncoords(1)==0) then
              istat = create_folder(k_folder)
              open(unit=123, file=trim(k_folder)//"/grid2d.xyz", form="unformatted", access="stream")
              write(123) nxmax,nymax,1
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
            starts(1) = 0 + ncoords(1)*subsizes(1)
            starts(2) = 0 + ncoords(2)*subsizes(2)
            starts(3) = 0
            ntot = nx*ny
            call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,mpi_prec,filetype,iermpi)
            call MPI_TYPE_COMMIT(filetype,iermpi)
            call MPI_FILE_OPEN(mp_cartx,trim(k_folder)//'/grid2d.xyz',MPI_MODE_RDWR,MPI_INFO_NULL,mpi_io_file,iermpi)
            offset = 3*size_integer
            do l=1,3
              call MPI_FILE_SET_VIEW(mpi_io_file,offset,mpi_prec,filetype,"native",MPI_INFO_NULL,iermpi)
              call MPI_FILE_WRITE_ALL(mpi_io_file,grid(l,1:nx,1:ny,k_slice:k_slice),ntot,mpi_prec,istatus,iermpi)
              call MPI_TYPE_SIZE(mpi_prec,size_real,iermpi)
              do m=1,nblocks(1)*nblocks(2)
                offset = offset+size_real*ntot
              enddo
            enddo
            call MPI_FILE_CLOSE(mpi_io_file,iermpi)
            call MPI_TYPE_FREE(filetype,iermpi)
          enddo
        endif

        if (allocated(self%islice_p3d)) then
          isize = size(self%islice_p3d)
          do ii=1,isize
            i_slice = self%islice_p3d(ii)
            i_folder = "SLICEYZ_??????"
            write(i_folder(9:14),'(I6.6)') i_slice + ncoords(1)*nx
            if (ncoords(3)==0) then
              istat = create_folder(i_folder)
              open(unit=123, file=trim(i_folder)//"/grid2d.xyz", form="unformatted", access="stream")
              write(123) 1,nymax,nzmax
              flush(123)
              close(123)
            endif
            call mpi_barrier(mp_cartz,iermpi)
            sizes(1) = 1
            sizes(2) = nblocks(2)*ny
            sizes(3) = nblocks(3)*nz
            subsizes(1) = 1
            subsizes(2) = ny
            subsizes(3) = nz
            starts(1) = 0
            starts(2) = 0 + ncoords(2)*subsizes(2)
            starts(3) = 0 + ncoords(3)*subsizes(3)
            ntot = ny*nz
            call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,mpi_prec,filetype,iermpi)
            call MPI_TYPE_COMMIT(filetype,iermpi)
            call MPI_FILE_OPEN(mp_cartz,trim(i_folder)//'/grid2d.xyz',MPI_MODE_RDWR,MPI_INFO_NULL,mpi_io_file,iermpi)
            offset = 3*size_integer
            do l=1,3
              call MPI_FILE_SET_VIEW(mpi_io_file,offset,mpi_prec,filetype,"native",MPI_INFO_NULL,iermpi)
              call MPI_FILE_WRITE_ALL(mpi_io_file,grid(l,i_slice:i_slice,1:ny,1:nz),ntot,mpi_prec,istatus,iermpi)
              call MPI_TYPE_SIZE(mpi_prec,size_real,iermpi)
              do m=1,nblocks(2)*nblocks(3)
                offset = offset+size_real*ntot
              enddo
            enddo
            call MPI_FILE_CLOSE(mpi_io_file,iermpi)
            call MPI_TYPE_FREE(filetype,iermpi)
          enddo
        endif

      endif p3d_grid

      p3d_field: if (plot3dfield_) then

        write(nastore(1:6), '(I6.6)') itslice_p3d
        if (masterproc) write(*,*) 'Storing slices at time', time

        if (allocated(self%jslice_p3d)) then
          jsize = size(self%jslice_p3d)
          do jj=1,jsize
            j_slice = self%jslice_p3d(jj)
            j_folder = "SLICEXZ_??????"
            write(j_folder(9:14),'(I6.6)') j_slice
            if (masterproc) then
              open(unit=123, file=trim(j_folder)//'/slicexz_'//nastore//'.q', form="unformatted", access="stream")
              write(123) nxmax,1,nzmax,num_aux_slice
              flush(123)
              close(123)
            endif
            call mpi_barrier(mpi_comm_world,iermpi)
            sizes(1) = nblocks(1)*nx
            sizes(2) = 1
            sizes(3) = nblocks(3)*nz
            subsizes(1) = nx
            subsizes(2) = 1
            subsizes(3) = nz
            starts(1) = 0 + ncoords(1)*subsizes(1)
            starts(2) = 0
            starts(3) = 0 + ncoords(3)*subsizes(3)
            ntot = nx*nz
            call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,mpi_prec,filetype,iermpi)
            call MPI_TYPE_COMMIT(filetype,iermpi)
            call MPI_FILE_OPEN(mp_cart,trim(j_folder)//'/slicexz_'//nastore//'.q',MPI_MODE_RDWR,MPI_INFO_NULL,mpi_io_file,iermpi)
            offset = 4*size_integer
            do l=1,num_aux_slice
              ll = self%list_aux_slice(l)
              call MPI_FILE_SET_VIEW(mpi_io_file,offset,mpi_prec,filetype,"native",MPI_INFO_NULL,iermpi)
              if (.not.(self%grid%is_y_staggered) .and. j_slice == 1 .and. ll >=2 .and. ll <=4) then
                call MPI_FILE_WRITE_ALL(mpi_io_file,self%wallprop(1:nx,1:nz,ll),ntot,mpi_prec,istatus,iermpi)
              else
                call MPI_FILE_WRITE_ALL(mpi_io_file,w_aux(1:nx,j_slice:j_slice,1:nz,ll),ntot,mpi_prec,istatus,iermpi)
              endif
              do m=1,nblocks(1)*nblocks(3)
                offset = offset+size_real*ntot
              enddo
            enddo
            call MPI_FILE_CLOSE(mpi_io_file,iermpi)
            call MPI_TYPE_FREE(filetype,iermpi)
          enddo
        endif
        if (allocated(self%kslice_p3d)) then
          ksize = size(self%kslice_p3d)
          do kk=1,ksize
            k_slice = self%kslice_p3d(kk)
            k_folder = "SLICEXY_??????"
            write(k_folder(9:14),'(I6.6)') k_slice + ncoords(3)*nz
            if (ncoords(1)==0) then
              istat = create_folder(k_folder)
              open(unit=123, file=trim(k_folder)//'/slicexy_'//nastore//'.q', form="unformatted", access="stream")
              write(123) nxmax,nymax,1,num_aux_slice
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
            starts(1) = 0 + ncoords(1)*subsizes(1)
            starts(2) = 0 + ncoords(2)*subsizes(2)
            starts(3) = 0
            ntot = nx*ny
            call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,mpi_prec,filetype,iermpi)
            call MPI_TYPE_COMMIT(filetype,iermpi)
            call MPI_FILE_OPEN(mp_cartx,trim(k_folder)//'/slicexy_'//nastore//'.q',MPI_MODE_RDWR,MPI_INFO_NULL, mpi_io_file,iermpi)
            offset = 4*size_integer
            do l=1,num_aux_slice
              ll = self%list_aux_slice(l)
              call MPI_FILE_SET_VIEW(mpi_io_file,offset,mpi_prec,filetype,"native",MPI_INFO_NULL,iermpi)
              call MPI_FILE_WRITE_ALL(mpi_io_file,w_aux(1:nx,1:ny,k_slice:k_slice,ll),ntot,mpi_prec,istatus,iermpi)
              do m=1,nblocks(1)*nblocks(2)
                offset = offset+size_real*ntot
              enddo
            enddo
            call MPI_FILE_CLOSE(mpi_io_file,iermpi)
            call MPI_TYPE_FREE(filetype,iermpi)
          enddo
        endif

        if (allocated(self%islice_p3d)) then
          isize = size(self%islice_p3d)
          do ii=1,isize
            i_slice = self%islice_p3d(ii)
            i_folder = "SLICEYZ_??????"
            write(i_folder(9:14),'(I6.6)') i_slice + ncoords(1)*nx
            if (ncoords(3)==0) then
              istat = create_folder(i_folder)
              open(unit=123, file=trim(i_folder)//'/sliceyz_'//nastore//'.q', form="unformatted", access="stream")
              write(123) 1,nymax,nzmax,num_aux_slice
              flush(123)
              close(123)
            endif
            call mpi_barrier(mp_cartz,iermpi)
            sizes(1) = 1
            sizes(2) = nblocks(2)*ny
            sizes(3) = nblocks(3)*nz
            subsizes(1) = 1
            subsizes(2) = ny
            subsizes(3) = nz
            starts(1) = 0
            starts(2) = 0 + ncoords(2)*subsizes(2)
            starts(3) = 0 + ncoords(3)*subsizes(3)
            ntot = ny*nz

            call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,mpi_prec,filetype,iermpi)
            call MPI_TYPE_COMMIT(filetype,iermpi)
            call MPI_FILE_OPEN(mp_cartz,trim(i_folder)//'/sliceyz_'//nastore//'.q',MPI_MODE_RDWR,MPI_INFO_NULL, mpi_io_file,iermpi)
            offset = 4*size_integer
            do l=1,num_aux_slice
              ll = self%list_aux_slice(l)
              call MPI_FILE_SET_VIEW(mpi_io_file,offset,mpi_prec,filetype,"native",MPI_INFO_NULL,iermpi)
              call MPI_FILE_WRITE_ALL(mpi_io_file,w_aux(i_slice:i_slice,1:ny,1:nz,ll),ntot,mpi_prec,istatus,iermpi)
              call MPI_TYPE_SIZE(mpi_prec,size_real,iermpi)
              do m=1,nblocks(2)*nblocks(3)
                offset = offset+size_real*ntot
              enddo
            enddo
            call MPI_FILE_CLOSE(mpi_io_file,iermpi)
            call MPI_TYPE_FREE(filetype,iermpi)
          enddo
        endif

      endif p3d_field

    endassociate
  endsubroutine write_slice_p3d
  subroutine write_tspec(self)
    class(equation_singleideal_object), intent(inout) :: self
    character(4) :: chx,chz
    associate(nx => self%field%nx, ny => self%field%ny, nz => self%field%nz, nv_stat => self%nv_stat&
    &, w_tspec => self%w_tspec, w_psd_avg_tspec => self%w_psd_avg_tspec, ncoords => self%field%ncoords)

      if (self%masterproc) write(*,*) 'Writing rst_tspec'
      write(chx,'(I4.4)') ncoords(1)
      write(chz,'(I4.4)') ncoords(3)

      open (11,file='rst1_tspec_'//chx//'_'//chz//'.bin',form='unformatted')
      write(11) w_tspec
      write(11) w_psd_avg_tspec
      close(11)

    endassociate
  endsubroutine write_tspec
  subroutine read_tspec(self)
    class(equation_singleideal_object), intent(inout) :: self
    character(4) :: chx,chz
    associate(nx => self%field%nx, ny => self%field%ny, nz => self%field%nz, nv_stat => self%nv_stat&
    &, w_tspec => self%w_tspec, w_psd_avg_tspec => self%w_psd_avg_tspec, ncoords => self%field%ncoords)

      if (self%masterproc) write(*,*) 'Reading rst_tspec'
      write(chx,'(I4.4)') ncoords(1)
      write(chz,'(I4.4)') ncoords(3)

      open (11,file='rst0_tspec_'//chx//'_'//chz//'.bin',form='unformatted')
      read(11) w_tspec
      read(11) w_psd_avg_tspec
      close(11)

    endassociate
  endsubroutine read_tspec
  subroutine write_psd_avg_tspec(self)
    class(equation_singleideal_object), intent(inout) :: self
    character(4) :: chx,chz
    integer :: i, idft, nfreq, l
    real(rkind) :: freq, df, prms_tspec_cha
    real(rkind), dimension(self%field%nx) :: prms_tspec
    real(rkind), dimension(self%ndft_tspec) :: w_psd_avg_tspec_cha
    associate(nx => self%field%nx, ny => self%field%ny, nz => self%field%nz, nv_stat => self%nv_stat&
    &, w_tspec => self%w_tspec, w_psd_avg_tspec => self%w_psd_avg_tspec, nxmax => self%grid%nxmax,&
    & ncoords => self%field%ncoords, ndft => self%ndft_tspec, dt_tspec => self%dt_tspec,&
    & mp_cartx => self%field%mp_cartx, mpi_err => self%mpi_err)

      if (self%masterproc) write(*,*) 'Writing psd_tspec'
      write(chx,'(I4.4)') ncoords(1)

      if(ncoords(3) == 0) then
        nfreq = ndft/2
        df = 1._rkind/(ndft*dt_tspec)

        if (self%channel_case) then
          prms_tspec_cha = 0._rkind
          do i=1,nx
            prms_tspec_cha = w_psd_avg_tspec(i,nfreq+1)*df
            do l=1,nfreq-1
              prms_tspec_cha = prms_tspec_cha+2._rkind*w_psd_avg_tspec(i,l+1)*df
            enddo
          enddo
          call mpi_allreduce(MPI_IN_PLACE,prms_tspec_cha,1,mpi_prec,mpi_sum,mp_cartx,mpi_err)
          prms_tspec_cha = prms_tspec_cha/nxmax

          if(self%masterproc) then
            open(10,file='prms_tspec_online_avg.dat')
            write(10,*) sqrt(prms_tspec_cha)
            close(10)
          endif

          w_psd_avg_tspec_cha = 0._rkind
          do i=1,nx
            w_psd_avg_tspec_cha(:) = w_psd_avg_tspec_cha(:)+w_psd_avg_tspec(i,:)
          enddo
          call mpi_allreduce(MPI_IN_PLACE,w_psd_avg_tspec_cha,ndft,mpi_prec,mpi_sum,mp_cartx,mpi_err)
          w_psd_avg_tspec_cha = w_psd_avg_tspec_cha/nxmax

          if(self%masterproc) then
            open (11,file='psd_tspec_online_avg.dat',form='formatted')
            do idft=1+1,nfreq+1
              freq = 0.5_rkind*(idft-1)/dt_tspec/nfreq
              write(11,'(1(I0,2X),3(E15.8,2X))') idft,freq,w_psd_avg_tspec_cha(idft),w_psd_avg_tspec_cha(idft)/prms_tspec_cha
            enddo
            close(11)
          endif

        else
          open(10,file='prms_tspec_online_avg_'//chx//'.dat')
          prms_tspec = 0._rkind
          do i=1,nx
            prms_tspec(i) = w_psd_avg_tspec(i,nfreq+1)*df
            do l=1,nfreq-1
              prms_tspec(i) = prms_tspec(i)+2._rkind*w_psd_avg_tspec(i,l+1)*df
            enddo
            write(10,*) i, sqrt(prms_tspec(i))
          enddo
          close(10)

          open (11,file='psd_tspec_online_avg_'//chx//'.dat',form='formatted')
          do i=1,nx
            do idft=1+1,ndft/2+1
              freq = 0.5_rkind*(idft-1)/dt_tspec/nfreq
              write(11,'(2(I0,2X),3(E15.8,2X))') i,idft,freq,w_psd_avg_tspec(i,idft),w_psd_avg_tspec(i,idft)/prms_tspec(i)
            enddo
          enddo
          close(11)

        endif
      endif

    endassociate
  endsubroutine write_psd_avg_tspec
  subroutine write_psd_tspec(self)
    class(equation_singleideal_object), intent(inout) :: self
    character(4) :: chx,chz
    character(128) :: filename
    integer :: i, idft, nfreq, l
    real(rkind) :: freq, df, prms_tspec_cha
    real(rkind), dimension(self%field%nx) :: prms_tspec
    real(rkind), dimension(self%ndft_tspec) :: w_psd_tspec_cha
    associate(nx => self%field%nx, ny => self%field%ny, nz => self%field%nz, nv_stat => self%nv_stat&
    &, w_tspec => self%w_tspec, w_psd_tspec => self%w_psd_tspec, nxmax => self%grid%nxmax,&
    & it_nwin_tspec => self%it_nwin_tspec, ncoords => self%field%ncoords, ndft => self%ndft_tspec,&
    & dt_tspec => self%dt_tspec, mp_cartx => self%field%mp_cartx, mpi_err => self%mpi_err)

      if (self%masterproc) write(*,*) 'Writing psd_tspec'
      write(chx,'(I4.4)') ncoords(1)

      if(ncoords(3) == 0) then
        nfreq = ndft/2
        df = 1._rkind/(ndft*dt_tspec)

        if (self%channel_case) then
          prms_tspec_cha = 0._rkind
          do i=1,nx
            prms_tspec_cha = w_psd_tspec(i,nfreq+1)*df
            do l=1,nfreq-1
              prms_tspec_cha = prms_tspec_cha+2._rkind*w_psd_tspec(i,l+1)*df
            enddo
          enddo
          call mpi_allreduce(MPI_IN_PLACE,prms_tspec_cha,1,mpi_prec,mpi_sum,mp_cartx,mpi_err)
          prms_tspec_cha = prms_tspec_cha/nxmax

          if(self%masterproc) then
            filename = 'prms_tspec_online_YYYY.dat'
            write(filename(19:22),'(I4.4)') it_nwin_tspec
            open(10,file=filename)
            write(10,*) sqrt(prms_tspec_cha)
            close(10)
          endif

          w_psd_tspec_cha = 0._rkind
          do i=1,nx
            w_psd_tspec_cha(:) = w_psd_tspec_cha(:)+w_psd_tspec(i,:)
          enddo
          call mpi_allreduce(MPI_IN_PLACE,w_psd_tspec_cha,ndft,mpi_prec,mpi_sum,mp_cartx,mpi_err)
          w_psd_tspec_cha = w_psd_tspec_cha/nxmax

          if(self%masterproc) then
            filename = 'psd_tspec_online_YYYY.dat'
            write(filename(18:21),'(I4.4)') it_nwin_tspec
            open (11,file=filename,form='formatted')
            do idft=1+1,nfreq+1
              freq = 0.5_rkind*(idft-1)/dt_tspec/nfreq
              write(11,'(1(I0,2X),3(E15.8,2X))') idft,freq,w_psd_tspec_cha(idft),w_psd_tspec_cha(idft)/prms_tspec_cha
            enddo
            close(11)
          endif

        else
          filename = 'prms_tspec_online_YYYY_'//chx//'.dat'
          write(filename(19:22),'(I4.4)') it_nwin_tspec
          open(10,file=filename)
          prms_tspec = 0._rkind
          do i=1,nx
            prms_tspec(i) = w_psd_tspec(i,nfreq+1)*df
            do l=1,nfreq-1
              prms_tspec(i) = prms_tspec(i)+2._rkind*w_psd_tspec(i,l+1)*df
            enddo
            write(10,*) i, sqrt(prms_tspec(i))
          enddo
          close(10)

          filename = 'psd_tspec_online_YYYY_'//chx//'.dat'
          write(filename(18:21),'(I4.4)') it_nwin_tspec
          open (11,file=filename,form='formatted')
          do i=1,nx
            do idft=1+1,ndft/2+1
              freq = 0.5_rkind*(idft-1)/dt_tspec/nfreq
              write(11,'(2(I0,2X),3(E15.8,2X))') i,idft,freq,w_psd_tspec(i,idft),w_psd_tspec(i,idft)/prms_tspec(i)
            enddo
          enddo
          close(11)

        endif
      endif

    endassociate
  endsubroutine write_psd_tspec









endmodule streams_equation_singleideal_object

