module streams_equation_singleideal_gpu_object

  use streams_base_gpu_object
  use streams_field_object
  use streams_grid_object
  use streams_kernels_gpu
  use streams_parameters
  use crandom_f_mod
  use streams_equation_singleideal_object
  use mpi
  use cudafor
  use iso_c_binding
  use, intrinsic :: iso_fortran_env
  use tcp

  use catalyst_api
  use catalyst_conduit
  use catalyst2_conduit_wrappers

  implicit none
  private
  public :: equation_singleideal_gpu_object

  integer(kind=cuda_stream_kind) :: stream1

  integer, parameter :: eulercentral_threads_x=128,eulercentral_threads_y=3
  integer, parameter :: eulerweno_threads_x=128 ,eulerweno_threads_y=2

  type :: equation_singleideal_gpu_object
    type(base_gpu_object) :: base_gpu
    type(equation_singleideal_object) :: equation_base
    type(field_object), pointer :: field=>null()
    type(grid_object), pointer :: grid=>null()
    integer(ikind) :: ng
    integer(ikind) :: nx
    integer(ikind) :: ny
    integer(ikind) :: nz
    integer(ikind) :: nv
    integer(ikind) :: nv_aux
    integer(ikind) :: nprocs
    integer(ikind) :: myrank
    integer(ikind) :: error
    real(rkind) :: time0
    real(rkind), allocatable, dimension(:,:), device :: coeff_deriv1_gpu
    real(rkind), allocatable, dimension(:,:), device :: coeff_deriv2_gpu
    real(rkind), allocatable, dimension(:) , device :: coeff_filter_gpu

    integer :: ierr
    integer :: icyc0, num_iter
    integer :: visc_model
    real(rkind) :: mu0, sutherland_s , t_ref_dim, powerlaw_vtexp
    logical :: masterproc
    integer :: mpi_err

    real(rkind), allocatable, dimension(:,:,:,:), device :: w_aux_gpu
    integer, allocatable, dimension(:,:,:), device :: fluid_mask_gpu
    integer, allocatable, dimension(:,:,:,:), device :: ep_ord_change_gpu

    real(rkind), allocatable, dimension(:), device :: winf_gpu, winf_past_shock_gpu
    real(rkind), allocatable, dimension(:,:,:,:) :: w_aux_debug
    real(rkind), allocatable, dimension(:,:,:,:) :: fl_debug
    real(rkind), allocatable, dimension(:,:,:,:), device :: fhat_gpu
    real(rkind), allocatable, dimension(:,:,:,:), device :: w_aux_trans_gpu, fl_trans_gpu, fhat_trans_gpu
    real(rkind), allocatable, dimension(:,:,:,:), device :: fl_gpu, fln_gpu
    real(rkind), allocatable, dimension(:,:), device :: dcoe_gpu

    real(rkind), allocatable, dimension(:,:,:,:), device :: wrecyc_gpu
    real(rkind), allocatable, dimension(:,:,:), device :: wrecycav_gpu
    real(rkind), dimension(:,:,:), allocatable, device :: wmean_gpu

    real(rkind), dimension(:,:,:), allocatable, device :: inflow_random_plane_gpu
    real(rkind), dimension(:), allocatable, device :: weta_inflow_gpu
    real(rkind), dimension(:), allocatable, device :: yplus_inflow_gpu, eta_inflow_gpu
    real(rkind), dimension(:), allocatable, device :: yplus_recyc_gpu, eta_recyc_gpu, eta_recyc_blend_gpu
    integer, dimension(:), allocatable, device :: map_j_inn_gpu, map_j_out_gpu, map_j_out_blend_gpu

    real(rkind), dimension(:), allocatable, device :: cv_coeff_gpu, cp_coeff_gpu

    real(rkind), dimension(:,:), allocatable, device :: w_aux_probe_gpu
    real(rkind), dimension(:,:,:,:), allocatable, device :: probe_coeff_gpu
    integer, dimension(:,:), allocatable, device :: ijk_probe_gpu
    real(rkind), allocatable, dimension(:,:,:), device :: wallprop_gpu

    integer, dimension(:), allocatable, device :: lmax_tag_gpu
    integer, dimension(:), allocatable, device :: vis_tag_gpu

    real(rkind), dimension(:,:,:,:,:), allocatable, device :: w_tspec_gpu
    real(rkind), dimension(:,:), allocatable, device :: w_psd_tspec_gpu
    real(rkind), dimension(:,:), allocatable, device :: wfar_gpu
    real(rkind), dimension(:), allocatable, device :: f_sponge_gpu


    !insitu_var_start
    real(rkind), dimension(:,:,:,:), allocatable, device :: psi_gpu
    integer, dimension(:), allocatable, device :: aux_list_gpu, add_list_gpu
    !insitu_var_end

  contains
    procedure, pass(self) :: compute_dt
    procedure, pass(self) :: initialize
    procedure, pass(self) :: compute_residual
    procedure, pass(self) :: compute_airfoil_forces_runtime
    procedure, pass(self) :: print_progress
    procedure, pass(self) :: run
    procedure, pass(self) :: rk_sync
    procedure, pass(self) :: rk_async
    procedure, pass(self) :: point_to_field
    procedure, pass(self) :: point_to_grid
    procedure, pass(self) :: alloc
    procedure, pass(self) :: update_ghost
    procedure, pass(self) :: force_rhs
    procedure, pass(self) :: force_var
    procedure, pass(self) :: euler_x
    procedure, pass(self) :: euler_y
    procedure, pass(self) :: euler_z
    procedure, pass(self) :: visflx
    procedure, pass(self) :: compute_aux
    procedure, pass(self) :: recyc_exchange
    procedure, pass(self) :: bcrecyc
    procedure, pass(self) :: bc_nr
    procedure, pass(self) :: manage_output
    procedure, pass(self) :: debug_airfoil
    procedure, pass(self) :: limiter
    procedure, pass(self) :: filter
    procedure, pass(self) :: tripping
    procedure, pass(self) :: bcextr_var
    procedure, pass(self) :: analyze_runtime
    procedure, pass(self) :: count_weno
    procedure, pass(self) :: sponge
    !insitu_proc_start
    procedure, pass(self) :: insitu_alloc_gpu
    procedure, pass(self) :: insitu_coprocess
    procedure, pass(self) :: insitu_compute_psi
    procedure, pass(self) :: insitu_do_catalyst_execute
    !insitu_proc_end

  endtype equation_singleideal_gpu_object

contains


  subroutine point_to_field(self, field)
    class(equation_singleideal_gpu_object), intent(inout) :: self
    class(field_object), target :: field
    self%field => field
  endsubroutine point_to_field

  subroutine point_to_grid(self, grid)
    class(equation_singleideal_gpu_object), intent(inout) :: self
    class(grid_object), target :: grid
    self%grid => grid
  endsubroutine point_to_grid




  subroutine insitu_alloc_gpu(self)
    class(equation_singleideal_gpu_object), intent(inout) :: self
    associate( nxsl_ins => self%equation_base%nxsl_ins, nxel_ins => self%equation_base%nxel_ins,&
    & nysl_ins => self%equation_base%nysl_ins, nyel_ins => self%equation_base%nyel_ins,&
    & nzsl_ins => self%equation_base%nzsl_ins, nzel_ins => self%equation_base%nzel_ins,&
    & npsi => self%equation_base%npsi, npsi_pv => self%equation_base%npsi_pv, nx => self%field%nx,&
    & ny => self%field%ny, nz => self%field%nz, ng => self%grid%ng, n_aux_list => self%equation_base%n_au&
    &x_list, n_add_list => self%equation_base%n_add_list )
      if (npsi > 0) allocate(self%psi_gpu(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,npsi))
      if (n_aux_list > 0) then
        allocate(self%aux_list_gpu(1:n_aux_list))
        self%aux_list_gpu = self%equation_base%aux_list
      endif
      if (n_add_list > 0) then
        allocate(self%add_list_gpu(1:n_add_list))
        self%add_list_gpu = self%equation_base%add_list
      endif
    endassociate
  endsubroutine insitu_alloc_gpu

  subroutine visflx(self, mode)
    class(equation_singleideal_gpu_object), intent(inout) :: self
    integer, intent(in) :: mode
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv,&
    & nv_aux => self%nv_aux, dt => self%equation_base%dt, visc_order => self%equation_base%visc_order,&
    & prandtl => self%equation_base%prandtl, visc_model => self%visc_model, mu0 => self%mu0,&
    & u0 => self%equation_base%u0, l0 => self%equation_base%l0, l0_ducros => self%equation_base%l0_ducros&
    &, t_ref_dim => self%t_ref_dim, sutherland_s => self%sutherland_s, powerlaw_vtexp => self%powerlaw_vt&
    &exp, coeff_deriv1_gpu => self%coeff_deriv1_gpu, coeff_deriv2_gpu => self%coeff_deriv2_gpu,&
    & fhat_trans_gpu => self%fhat_trans_gpu, fl_trans_gpu => self%fl_trans_gpu, fl_gpu => self%fl_gpu,&
    & w_aux_gpu => self%w_aux_gpu, w_aux_trans_gpu => self%w_aux_trans_gpu, dcsidx_gpu => self%base_gpu%d&
    &csidx_gpu, detady_gpu => self%base_gpu%detady_gpu, dzitdz_gpu => self%base_gpu%dzitdz_gpu,&
    & dcsidxs_gpu => self%base_gpu%dcsidxs_gpu, detadys_gpu => self%base_gpu%detadys_gpu,&
    & dzitdzs_gpu => self%base_gpu%dzitdzs_gpu, dcsidx2_gpu => self%base_gpu%dcsidx2_gpu,&
    & detady2_gpu => self%base_gpu%detady2_gpu, dzitdz2_gpu => self%base_gpu%dzitdz2_gpu,&
    & x_gpu => self%base_gpu%x_gpu, y_gpu => self%base_gpu%y_gpu, z_gpu => self%base_gpu%z_gpu,&
    & eul_imin => self%equation_base%eul_imin, eul_imax => self%equation_base%eul_imax,&
    & eul_jmin => self%equation_base%eul_jmin, eul_jmax => self%equation_base%eul_jmax,&
    & eul_kmin => self%equation_base%eul_kmin, eul_kmax => self%equation_base%eul_kmax,&
    & cv_coeff_gpu => self%cv_coeff_gpu, cp_coeff_gpu => self%cp_coeff_gpu, indx_cp_l => self%equation_ba&
    &se%indx_cp_l, indx_cp_r => self%equation_base%indx_cp_r, calorically_perfect => self%equation_base%c&
    &alorically_perfect, t0 => self%equation_base%t0, channel_case => self%equation_base%channel_case,&
    & grid_dim => self%grid%grid_dim, dcsidxc2_gpu => self%base_gpu%dcsidxc2_gpu, dcsidyc2_gpu => self%ba&
    &se_gpu%dcsidyc2_gpu, detadxc2_gpu => self%base_gpu%detadxc2_gpu, detadyc2_gpu => self%base_gpu%detad&
    &yc2_gpu, g1_gpu => self%base_gpu%g1_gpu, g2_gpu => self%base_gpu%g2_gpu, g12_gpu => self%base_gpu%g1&
    &2_gpu, jac_gpu => self%base_gpu%jac_gpu, vis_tag_gpu => self%vis_tag_gpu, wall_tag_gpu => self%base_&
    &gpu%wall_tag_gpu, iblock => self%field%ncoords(1), ite_rank_x => self%field%ite_rank_x,&
    & itu_rank_x => self%field%itu_rank_x, ite_l => self%field%ite_l, itu_l => self%field%itu_l,&
    & teshk => self%equation_base%teshk, ortho => self%equation_base%ortho, jweno => self%equation_base%j&
    &weno, theta_ij_gpu => self%base_gpu%theta_ij_gpu, theta_threshold => self%equation_base%theta_thresh&
    &old)

      if (mode == 0) then
        if (grid_dim == 1) then
          call visflx_cuf(nx, ny, nz, ng, visc_order, prandtl, t0, indx_cp_l, indx_cp_r,&
          & cp_coeff_gpu, calorically_perfect, u0, l0, self%base_gpu%w_gpu, self%w_aux_gpu, self%fl_gpu,&
          & coeff_deriv1_gpu, coeff_deriv2_gpu, dcsidx_gpu, detady_gpu, dzitdz_gpu, dcsidxs_gpu, detadys_gpu,&
          & dzitdzs_gpu, dcsidx2_gpu, detady2_gpu, dzitdz2_gpu, self%wallprop_gpu)
        elseif (grid_dim == 2) then
          call visflx_c2_cuf(nx, ny, nz, ng, visc_order, prandtl, t0, indx_cp_l, indx_cp_r,&
          & cp_coeff_gpu, calorically_perfect, u0, l0, self%base_gpu%w_gpu, self%w_aux_gpu, self%fl_gpu,&
          & coeff_deriv1_gpu, coeff_deriv2_gpu, dcsidxc2_gpu, detadyc2_gpu, detadxc2_gpu, dcsidyc2_gpu,&
          & dzitdz_gpu, dzitdzs_gpu, dzitdz2_gpu, g1_gpu, g2_gpu, g12_gpu, jac_gpu, self%wallprop_gpu,&
          & vis_tag_gpu, wall_tag_gpu, iblock, ite_rank_x, itu_rank_x, ite_l, itu_l, teshk )
        endif
      elseif (mode == 1) then
        if (grid_dim == 1) then
          call visflx_div_cuf(nx, ny, nz, ng, visc_order, self%w_aux_gpu, self%fl_gpu,&
          & coeff_deriv1_gpu, dcsidx_gpu, detady_gpu, dzitdz_gpu, stream1)
        elseif (grid_dim == 2) then
          call visflx_div_c2_cuf(nx, ny, nz, ng, visc_order, self%w_aux_gpu, self%fl_gpu,&
          & coeff_deriv1_gpu, dcsidxc2_gpu, dcsidyc2_gpu, detadxc2_gpu, detadyc2_gpu, dzitdz_gpu, vis_tag_gpu,&
          & wall_tag_gpu, stream1)
        endif
      elseif (mode == 2) then
        call visflx_reduced_ord2_cuf(nx, ny, nz, ng, prandtl, t0, indx_cp_l, indx_cp_r,&
        & cp_coeff_gpu, calorically_perfect, u0, l0, self%base_gpu%w_gpu, self%w_aux_gpu, self%fl_gpu, x_gpu,&
        & y_gpu, z_gpu, self%wallprop_gpu, 1)
      elseif (mode == 3) then
        if (grid_dim == 1) then
          call sensor_cuf(nx, ny, nz, ng, u0, l0, self%w_aux_gpu)
        elseif (grid_dim == 2) then
          call sensor_c2_cuf(nx, ny, nz, ng, u0, l0_ducros, self%w_aux_gpu, iblock, ite_rank_x,&
          & itu_rank_x, ite_l, itu_l, teshk, theta_ij_gpu, theta_threshold, jweno)
        endif
      elseif (mode == 4) then
        if (grid_dim == 1) then
          call visflx_nosensor_cuf(nx, ny, nz, ng, visc_order, prandtl, t0, indx_cp_l, indx_cp_r,&
          & cp_coeff_gpu, calorically_perfect, u0, l0, self%base_gpu%w_gpu, self%w_aux_gpu, self%fl_gpu,&
          & coeff_deriv1_gpu, coeff_deriv2_gpu, dcsidx_gpu, detady_gpu, dzitdz_gpu, dcsidxs_gpu, detadys_gpu,&
          & dzitdzs_gpu, dcsidx2_gpu, detady2_gpu, dzitdz2_gpu, self%wallprop_gpu)
        elseif (grid_dim == 2) then
          call visflx_nosensor_c2_cuf(nx, ny, nz, ng, visc_order, prandtl, t0, indx_cp_l, indx_cp_r,&
          & cp_coeff_gpu, calorically_perfect, iblock, ite_rank_x, itu_rank_x, ite_l, itu_l, u0, l0,&
          & self%base_gpu%w_gpu, self%w_aux_gpu, self%fl_gpu, coeff_deriv1_gpu, coeff_deriv2_gpu, dcsidxc2_gpu,&
          & detadyc2_gpu, detadxc2_gpu, dcsidyc2_gpu, dzitdz_gpu, dzitdzs_gpu, dzitdz2_gpu, g1_gpu, g2_gpu,&
          & g12_gpu, jac_gpu, self%wallprop_gpu, vis_tag_gpu, wall_tag_gpu, ortho)
        endif
      elseif (mode == 5) then
        call visflx_x_cuf(nx, ny, nz, nv, ng, prandtl, t0, indx_cp_l, indx_cp_r, cp_coeff_gpu,&
        & calorically_perfect, self%base_gpu%x_gpu, w_aux_gpu, self%fl_gpu, self%fhat_gpu)
        call visflx_y_cuf(nx, ny, nz, nv, ng, prandtl, t0, indx_cp_l, indx_cp_r, cp_coeff_gpu,&
        & calorically_perfect, self%base_gpu%y_gpu, self%w_aux_gpu, self%fl_gpu)
        call visflx_z_cuf(nx, ny, nz, nv, ng, prandtl, t0, indx_cp_l, indx_cp_r, cp_coeff_gpu,&
        & calorically_perfect, self%base_gpu%z_gpu, self%w_aux_gpu, self%fl_gpu)
      elseif (mode == 6) then
        call visflx_reduced_ord2_cuf(nx, ny, nz, ng, prandtl, t0, indx_cp_l, indx_cp_r,&
        & cp_coeff_gpu, calorically_perfect, u0, l0, self%base_gpu%w_gpu, self%w_aux_gpu, self%fl_gpu, x_gpu,&
        & y_gpu, z_gpu, self%wallprop_gpu, 0)
      elseif (mode == 7) then
        call visflx_div_ord2_cuf(nx, ny, nz, ng, self%w_aux_gpu, self%fl_gpu, x_gpu, y_gpu, z_gpu, stream1)
      endif
    endassociate
  endsubroutine visflx

  subroutine compute_aux(self, central, ghost)
    class(equation_singleideal_gpu_object), intent(inout) :: self
    integer, intent(in), optional :: central, ghost
    integer :: central_, ghost_
    central_ = 1 ; if(present(central)) central_ = central
    ghost_ = 1 ; if(present(ghost)) ghost_ = ghost
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, visc_model => self%visc_mo&
    &del, mu0 => self%mu0, t0 => self%equation_base%t0, t_ref_dim => self%t_ref_dim, sutherland_s => self&
    &%sutherland_s, powerlaw_vtexp => self%powerlaw_vtexp, cv_coeff_gpu => self%cv_coeff_gpu,&
    & cp_coeff_gpu => self%cp_coeff_gpu, indx_cp_l => self%equation_base%indx_cp_l, indx_cp_r => self%equ&
    &ation_base%indx_cp_r, calorically_perfect => self%equation_base%calorically_perfect,&
    & rgas0 => self%equation_base%rgas0)
      if(central_ == 1 .and. ghost_ == 1) then
        call eval_aux_cuf(nx, ny, nz, ng, 1-ng, nx+ng, 1-ng, ny+ng, 1-ng, nz+ng, self%base_gpu%w_gpu&
        &, self%w_aux_gpu, visc_model, mu0, t0, sutherland_s, t_ref_dim, powerlaw_vtexp, visc_power,&
        & visc_sutherland, visc_no, cv_coeff_gpu, indx_cp_l, indx_cp_r, rgas0, calorically_perfect,&
        & tol_iter_nr,stream1)
      elseif(central_ == 1 .and. ghost_ == 0) then
        call eval_aux_cuf(nx, ny, nz, ng, 1, nx, 1, ny, 1, nz, self%base_gpu%w_gpu, self%w_aux_gpu,&
        & visc_model, mu0, t0, sutherland_s, t_ref_dim, powerlaw_vtexp, visc_power, visc_sutherland, visc_no,&
        & cv_coeff_gpu, indx_cp_l, indx_cp_r, rgas0, calorically_perfect, tol_iter_nr,stream1)
      elseif(central_ == 0 .and. ghost_ == 1) then
        call eval_aux_cuf(nx, ny, nz, ng, 1-ng, 0, 1-ng, ny+ng, 1-ng, nz+ng, self%base_gpu%w_gpu,&
        & self%w_aux_gpu, visc_model, mu0, t0, sutherland_s, t_ref_dim, powerlaw_vtexp, visc_power,&
        & visc_sutherland, visc_no, cv_coeff_gpu, indx_cp_l, indx_cp_r, rgas0, calorically_perfect,&
        & tol_iter_nr,stream1)
        call eval_aux_cuf(nx, ny, nz, ng, nx+1, nx+ng, 1-ng, ny+ng, 1-ng, nz+ng, self%base_gpu%w_gpu&
        &, self%w_aux_gpu, visc_model, mu0, t0, sutherland_s, t_ref_dim, powerlaw_vtexp, visc_power,&
        & visc_sutherland, visc_no, cv_coeff_gpu, indx_cp_l, indx_cp_r, rgas0, calorically_perfect,&
        & tol_iter_nr,stream1)
        call eval_aux_cuf(nx, ny, nz, ng, 1-ng, nx+ng, 1-ng, 0, 1-ng, nz+ng, self%base_gpu%w_gpu,&
        & self%w_aux_gpu, visc_model, mu0, t0, sutherland_s, t_ref_dim, powerlaw_vtexp, visc_power,&
        & visc_sutherland, visc_no, cv_coeff_gpu, indx_cp_l, indx_cp_r, rgas0, calorically_perfect,&
        & tol_iter_nr,stream1)
        call eval_aux_cuf(nx, ny, nz, ng, 1-ng, nx+ng, ny+1, ny+ng, 1-ng, nz+ng, self%base_gpu%w_gpu&
        &, self%w_aux_gpu, visc_model, mu0, t0, sutherland_s, t_ref_dim, powerlaw_vtexp, visc_power,&
        & visc_sutherland, visc_no, cv_coeff_gpu, indx_cp_l, indx_cp_r, rgas0, calorically_perfect,&
        & tol_iter_nr,stream1)
        call eval_aux_cuf(nx, ny, nz, ng, 1-ng, nx+ng, 1-ng, ny+ng, 1-ng, 0, self%base_gpu%w_gpu,&
        & self%w_aux_gpu, visc_model, mu0, t0, sutherland_s, t_ref_dim, powerlaw_vtexp, visc_power,&
        & visc_sutherland, visc_no, cv_coeff_gpu, indx_cp_l, indx_cp_r, rgas0, calorically_perfect,&
        & tol_iter_nr,stream1)
        call eval_aux_cuf(nx, ny, nz, ng, 1-ng, nx+ng, 1-ng, ny+ng, nz+1, nz+ng, self%base_gpu%w_gpu&
        &, self%w_aux_gpu, visc_model, mu0, t0, sutherland_s, t_ref_dim, powerlaw_vtexp, visc_power,&
        & visc_sutherland, visc_no, cv_coeff_gpu, indx_cp_l, indx_cp_r, rgas0, calorically_perfect,&
        & tol_iter_nr,stream1)
      endif
    endassociate
  endsubroutine compute_aux

  subroutine rk_sync(self)
    class(equation_singleideal_gpu_object), intent(inout) :: self
    integer :: istep, lmax, iercuda
    real(rkind) :: rhodt, gamdt, alpdt
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv,&
    & dt => self%equation_base%dt, ep_order => self%equation_base%ep_order, weno_scheme => self%equation_&
    &base%weno_scheme, flow_init => self%equation_base%flow_init, a_tr => self%equation_base%a_tr,&
    & ileftx => self%field%ileftx, irightx => self%field%irightx, ileftz => self%field%ileftz,&
    & irightz => self%field%irightz, conservative_viscous => self%equation_base%conservative_viscous,&
    & eul_imin => self%equation_base%eul_imin, eul_imax => self%equation_base%eul_imax,&
    & eul_jmin => self%equation_base%eul_jmin, eul_jmax => self%equation_base%eul_jmax,&
    & eul_kmin => self%equation_base%eul_kmin, eul_kmax => self%equation_base%eul_kmax,&
    & channel_case => self%equation_base%channel_case, nv_aux => self%nv_aux, enable_limiter => self%equa&
    &tion_base%enable_limiter, enable_sponge => self%equation_base%enable_sponge, ndim => self%grid%ndim,&
    & grid_dim => self%grid%grid_dim)
      if (channel_case) then
        self%equation_base%dpdx = 0._rkind
        self%equation_base%dpth = 0._rkind
      endif

      do istep=1,self%equation_base%nrk
        rhodt = self%equation_base%rhork(istep)*dt
        gamdt = self%equation_base%gamrk(istep)*dt
        alpdt = self%equation_base%alprk(istep)*dt

        call init_flux_cuf(nx, ny, nz, nv, self%fl_gpu, self%fln_gpu, rhodt)
        call self%base_gpu%bcswap()
        call self%compute_aux()
        !@cuf iercuda=cudadevicesynchronize()
        call self%euler_x(eul_imin, eul_imax)
        !@cuf iercuda=cudadevicesynchronize()
        if (conservative_viscous== 1) then
          call self%visflx(mode=6)
          call self%visflx(mode=5)
        else
          call self%visflx(mode=4)
        endif
        if(flow_init == 5) then
          call self%base_gpu%bcswap_wake_var(self%w_aux_gpu(:,:,:,10:10))
          call self%base_gpu%bcswap_var(self%w_aux_gpu(:,:,:,10:10))
          call self%bcextr_var(self%w_aux_gpu(:,:,:,10:10), mode=2)
          call self%base_gpu%bcte_var(self%w_aux_gpu(:,:,:,10:10))
        else
          call self%bcextr_var(self%w_aux_gpu(:,:,:,10:10), mode=1)
          call self%base_gpu%bcswap_var(self%w_aux_gpu(:,:,:,10:10))
        endif
        call self%euler_y(eul_jmin,eul_jmax)
        !@cuf iercuda=cudadevicesynchronize()
        if (ndim == 3) call self%euler_z(eul_kmin,eul_kmax)
        if (istep == 3) then
          call self%visflx(mode=3)
          call self%base_gpu%bcswap_var(self%w_aux_gpu(:,:,:,8:8))
          if(flow_init == 5) then
            call self%base_gpu%bcswap_wake_var(self%w_aux_gpu(:,:,:,8:8))
          endif
        endif
        if (conservative_viscous==1) then
          call self%visflx(mode=7)
        else
          call self%visflx(mode=1)
        endif
        !@cuf iercuda=cudadevicesynchronize()
        call self%bc_nr()

        if (flow_init == 5 .and. a_tr>0._rkind) call self%tripping()

        call update_flux_cuf(nx, ny, nz, nv, self%fl_gpu, self%fln_gpu, gamdt)
        if (channel_case) call self%force_rhs()
        if (grid_dim == 2 .and. enable_sponge > 0) call self%sponge()
        call update_field_cuf(nx, ny, nz, ng, nv, self%base_gpu%w_gpu, self%fln_gpu, self%fluid_mask_gpu)
        if (grid_dim == 2 .and. enable_limiter > 0) call self%limiter()
        if (channel_case) call self%force_var()
        call self%update_ghost(do_swap=0)

      enddo

    endassociate
  endsubroutine rk_sync

  subroutine bcextr_var(self, w_swap, mode)
    class(equation_singleideal_gpu_object), intent(inout) :: self
    integer :: mode
    real(rkind), dimension(1-self%ng:self%nx+self%ng, 1-self%ng:self%ny+self%ng, 1-self%ng:self%nz+self%ng,1:1), device :: w_swap
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, ndim => self%grid%ndim,&
    & ileftx => self%field%ileftx, irightx => self%field%irightx, ileftz => self%field%ileftz,&
    & irightz => self%field%irightz)
      if(mode == 1) then
        call bcextr_var_cuf(nx, ny, nz, ng, w_swap)
      elseif(mode == 2) then
        call bcextr_airfoil_var_cuf(nx, ny, nz, ng, ndim, self%base_gpu%wall_tag_gpu, w_swap, ileftx, irightx, ileftz, irightz)
      endif
    endassociate
  endsubroutine bcextr_var

  subroutine rk_async(self)
    class(equation_singleideal_gpu_object), intent(inout) :: self
    integer :: istep, lmax, iercuda
    real(rkind) :: rhodt, gamdt, alpdt
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv,&
    & dt => self%equation_base%dt, ep_order => self%equation_base%ep_order, weno_scheme => self%equation_&
    &base%weno_scheme, flow_init => self%equation_base%flow_init, a_tr => self%equation_base%a_tr,&
    & ileftx => self%field%ileftx, irightx => self%field%irightx, ileftz => self%field%ileftz,&
    & irightz => self%field%irightz, conservative_viscous => self%equation_base%conservative_viscous,&
    & eul_imin => self%equation_base%eul_imin, eul_imax => self%equation_base%eul_imax,&
    & eul_jmin => self%equation_base%eul_jmin, eul_jmax => self%equation_base%eul_jmax,&
    & eul_kmin => self%equation_base%eul_kmin, eul_kmax => self%equation_base%eul_kmax,&
    & channel_case => self%equation_base%channel_case, nv_aux => self%nv_aux, enable_limiter => self%equa&
    &tion_base%enable_limiter, enable_sponge => self%equation_base%enable_sponge, ndim => self%grid%ndim,&
    & grid_dim => self%grid%grid_dim)
      if (channel_case) then
        self%equation_base%dpdx = 0._rkind
        self%equation_base%dpth = 0._rkind
      endif
      lmax = max(ep_order/2, weno_scheme)

      do istep=1,self%equation_base%nrk
        rhodt = self%equation_base%rhork(istep)*dt
        gamdt = self%equation_base%gamrk(istep)*dt
        alpdt = self%equation_base%alprk(istep)*dt

        call init_flux_cuf(nx, ny, nz, nv, self%fl_gpu, self%fln_gpu, rhodt)
        call self%compute_aux(central=1, ghost=0)
        !@cuf iercuda=cudadevicesynchronize()
        call self%base_gpu%bcswap(steps=[.true.,.false.,.false.])
        call self%euler_x(lmax+1,nx-lmax,lmax,nx-lmax,do_update=.false.)
        call self%base_gpu%bcswap(steps=[.false.,.true.,.true.])
        call self%compute_aux(central=0, ghost=1)
        call self%euler_x(eul_imin,lmax,eul_imin-1,lmax-1,do_update=.false.)
        call self%euler_x(nx-lmax+1,eul_imax,nx-lmax+1,eul_imax,do_update=.true.)
        !@cuf iercuda=cudadevicesynchronize()
        if (conservative_viscous== 1) then
          call self%visflx(mode=6)
          call self%visflx(mode=5)
        else
          call self%visflx(mode=4)
        endif
        if(flow_init == 5) then
          call self%base_gpu%bcswap_wake_var(self%w_aux_gpu(:,:,:,10:10))
          call self%base_gpu%bcswap_var(self%w_aux_gpu(:,:,:,10:10), steps=[.true.,.false.,.false.])
          call self%euler_y(eul_jmin,eul_jmax)
          call self%base_gpu%bcswap_var(self%w_aux_gpu(:,:,:,10:10), steps=[.false.,.true.,.true.])
          call self%bcextr_var(self%w_aux_gpu(:,:,:,10:10), mode=2)
          call self%base_gpu%bcte_var(self%w_aux_gpu(:,:,:,10:10))
        else
          call self%bcextr_var(self%w_aux_gpu(:,:,:,10:10), mode=1)
          call self%base_gpu%bcswap_var(self%w_aux_gpu(:,:,:,10:10), steps=[.true.,.false.,.false.])
          call self%euler_y(eul_jmin,eul_jmax)
          call self%base_gpu%bcswap_var(self%w_aux_gpu(:,:,:,10:10), steps=[.false.,.true.,.true.])
        endif
        if (ndim == 3) call self%euler_z(eul_kmin,eul_kmax)
        if (istep == 3) then
          call self%visflx(mode=3)
          if(flow_init == 5) then
            call self%base_gpu%bcswap_var(self%w_aux_gpu(:,:,:,8:8))
            call self%base_gpu%bcswap_wake_var(self%w_aux_gpu(:,:,:,8:8), steps=[.true.,.false.,.false.])
            if (conservative_viscous==1) then
              call self%visflx(mode=7)
            else
              call self%visflx(mode=1)
            endif
            call self%base_gpu%bcswap_wake_var(self%w_aux_gpu(:,:,:,8:8), steps=[.false.,.true.,.true.])
          else
            call self%base_gpu%bcswap_var(self%w_aux_gpu(:,:,:,8:8), steps=[.true.,.false.,.false.])
            if (conservative_viscous==1) then
              call self%visflx(mode=7)
            else
              call self%visflx(mode=1)
            endif
            call self%base_gpu%bcswap_var(self%w_aux_gpu(:,:,:,8:8), steps=[.false.,.true.,.true.])
          endif
        else
          if (conservative_viscous==1) then
            call self%visflx(mode=7)
          else
            call self%visflx(mode=1)
          endif
          !@cuf iercuda=cudadevicesynchronize()
        endif
        call self%bc_nr()

        if (flow_init == 5 .and. a_tr>0._rkind) call self%tripping()

        call update_flux_cuf(nx, ny, nz, nv, self%fl_gpu, self%fln_gpu, gamdt)
        if (channel_case) call self%force_rhs()
        if (grid_dim == 2 .and. enable_sponge > 0) call self%sponge()
        call update_field_cuf(nx, ny, nz, ng, nv, self%base_gpu%w_gpu, self%fln_gpu, self%fluid_mask_gpu)
        if (grid_dim == 2 .and. enable_limiter > 0) call self%limiter()
        if (channel_case) call self%force_var()
        call self%update_ghost(do_swap=0)

      enddo

    endassociate
  endsubroutine rk_async

  subroutine sponge(self)
    class(equation_singleideal_gpu_object), intent(inout) :: self
    associate(nx=>self%nx,ny=>self%ny,nz=>self%nz,ng=>self%ng,nv=>self%nv,w_gpu=>self%base_gpu%w_gpu&
    &, wfar_gpu=>self%wfar_gpu,fln_gpu=>self%fln_gpu,f_sponge_gpu=>self%f_sponge_gpu, j_sponge=>self%equa&
    &tion_base%j_sponge)

      call sponge_cuf(nx,ny,nz,ng,nv,w_gpu,wfar_gpu,fln_gpu,f_sponge_gpu,j_sponge)
    endassociate
  endsubroutine sponge

  subroutine limiter(self)
    class(equation_singleideal_gpu_object), intent(inout) :: self
    integer :: iblock,kblock
    associate(nx=>self%nx,ny=>self%ny,nz=>self%nz,ng=>self%ng,w_gpu=>self%base_gpu%w_gpu,&
    &w_aux_gpu=>self%w_aux_gpu, ncoords=>self%field%ncoords, indx_cp_l=>self%equation_base%indx_cp_l,&
    & indx_cp_r=>self%equation_base%indx_cp_r, cv_coeff_gpu => self%cv_coeff_gpu, calorically_perfect => &
    &self%equation_base%calorically_perfect, t0 => self%equation_base%t0, flow_init => self%equation_base&
    &%flow_init, rho_lim => self%equation_base%rho_lim, tem_lim => self%equation_base%tem_lim,&
    & rho_lim_rescale => self%equation_base%rho_lim_rescale, tem_lim_rescale => self%equation_base%tem_li&
    &m_rescale)
      iblock = ncoords(1)
      kblock = ncoords(3)

      call limiter_cuf(nx,ny,nz,ng,w_gpu,w_aux_gpu,iblock,kblock,indx_cp_l,indx_cp_r,cv_coeff_gpu,&
      &calorically_perfect,t0, tol_iter_nr, rho_lim, tem_lim, rho_lim_rescale, tem_lim_rescale)
    endassociate
  endsubroutine limiter

  subroutine filter(self)
    class(equation_singleideal_gpu_object), intent(inout) :: self

    associate(nx=>self%nx,ny=>self%ny,nz=>self%nz,ng=>self%ng,w_gpu=>self%base_gpu%w_gpu,&
    &w_aux_gpu=>self%w_aux_gpu,indx_cp_l=>self%equation_base%indx_cp_l,indx_cp_r=>self%equation_base%indx&
    &_cp_r,cv_coeff_gpu => self%cv_coeff_gpu,calorically_perfect => self%equation_base%calorically_perfec&
    &t,t0 => self%equation_base%t0, coeff_filter_gpu => self%coeff_filter_gpu,jfilter => self%equation_ba&
    &se%jfilter, wall_tag_gpu => self%base_gpu%wall_tag_gpu)


      call filter_cuf(nx,ny,nz,ng,w_gpu,w_aux_gpu,indx_cp_l,indx_cp_r,cv_coeff_gpu,&
      &calorically_perfect,t0,tol_iter_nr,coeff_filter_gpu,jfilter,wall_tag_gpu)
    endassociate
  endsubroutine filter



  subroutine tripping(self)
    class(equation_singleideal_gpu_object), intent(inout) :: self
    real(rkind) :: asl, del0,zita,lamx,lamy,tau,pp,bt,rand
    integer :: i_rank_start, ierr
    real(rkind), dimension(2) :: rsend
    associate(nx=>self%nx,ny=>self%ny,nz=>self%nz,ng=>self%ng,nv=>self%nv,w_gpu=>self%base_gpu%w_gpu&
    &,fl_gpu=>self%fl_gpu, xc2_gpu => self%base_gpu%xc2_gpu, yc2_gpu => self%base_gpu%yc2_gpu,&
    & z_gpu => self%base_gpu%z_gpu, ncoords=>self%field%ncoords, rlz => self%grid%domain_size(3),&
    & wall_tag_gpu => self%base_gpu%wall_tag_gpu, dxdetanc2_gpu => self%base_gpu%dxdetanc2_gpu,&
    & dydetanc2_gpu => self%base_gpu%dydetanc2_gpu, a_tr => self%equation_base%a_tr, u0 => self%equation_&
    &base%u0, time => self%equation_base%time, del0_tr => self%equation_base%del0_tr, masterproc => self%&
    &masterproc, itr1 => self%equation_base%itr1, itr2 => self%equation_base%itr2, its1 => self%equation_&
    &base%its1, its2 => self%equation_base%its2, x0tr => self%equation_base%x0tr, y0tr => self%equation_b&
    &ase%y0tr, x0ts => self%equation_base%x0ts, y0ts => self%equation_base%y0ts, lamz => self%equation_ba&
    &se%lamz , lamz1 => self%equation_base%lamz1 , lams => self%equation_base%lams , lams1 => self%equati&
    &on_base%lams1 , phiz => self%equation_base%phiz , phiz1 => self%equation_base%phiz1 ,&
    & phis => self%equation_base%phis , phis1 => self%equation_base%phis1, lamz_old => self%equation_base&
    &%lamz_old , lamz1_old => self%equation_base%lamz1_old , lams_old => self%equation_base%lams_old ,&
    & lams1_old => self%equation_base%lams1_old , phiz_old => self%equation_base%phiz_old ,&
    & phiz1_old => self%equation_base%phiz1_old , phis_old => self%equation_base%phis_old ,&
    & phis1_old => self%equation_base%phis1_old, is => self%equation_base%is, is_old => self%equation_bas&
    &e%is_old)
      asl = a_tr*u0*u0
      del0 = del0_tr
      zita = 1.7_rkind*del0
      lamx = 4._rkind*del0
      lamy = del0
      tau = lamx/u0
      is = int(time/tau)
      pp = time/tau-is
      bt = 3._rkind*pp*pp-2._rkind*pp*pp*pp
      if (is==is_old) then
        lamz = lamz_old
        lamz1 = lamz1_old
        phiz = phiz_old
        phiz1 = phiz1_old
        lams = lams_old
        lams1 = lams1_old
        phis = phis_old
        phis1 = phis1_old
      else
        lamz = lamz1_old
        phiz = phiz1_old
        if(masterproc) then
          call get_crandom_f(rand)
          lamz1 = 1._rkind/rlz*int(rand*rlz/zita)
          call get_crandom_f(rand)
          phiz1 = 2._rkind*pi*rand
          rsend(1) = lamz1
          rsend(2) = phiz1
        endif
        call mpi_bcast(rsend,2,mpi_prec,0,mpi_comm_world,ierr)
        if(.not.masterproc) then
          lamz1 = rsend(1)
          phiz1 = rsend(2)
        endif
        lamz_old = lamz
        lamz1_old = lamz1
        phiz_old = phiz
        phiz1_old = phiz1
        lams = lams1_old
        phis = phis1_old
        if(masterproc) then
          call get_crandom_f(rand)
          lams1 = 1._rkind/rlz*int(rand*rlz/zita)
          call get_crandom_f(rand)
          phis1 = 2._rkind*pi*rand
          rsend(1) = lams1
          rsend(2) = phis1
        endif
        call mpi_bcast(rsend,2,mpi_prec,0,mpi_comm_world,ierr)
        if(.not.masterproc) then
          lams1 = rsend(1)
          phis1 = rsend(2)
        endif
        lams_old = lams
        lams1_old = lams1
        phis_old = phis
        phis1_old = phis1
        is_old = is
      endif

      i_rank_start = nx*ncoords(1)
      call tripping_pressure_cuf(nx,ny,nz,ng,nv,i_rank_start,pi,itr1,itr2,x0tr,y0tr,x0ts,y0ts, lamx,&
      &lamy,lamz,lamz1,phiz,phiz1,asl,bt,xc2_gpu,yc2_gpu,z_gpu,w_gpu,fl_gpu,wall_tag_gpu, dxdetanc2_gpu,&
      & dydetanc2_gpu)

      call tripping_suction_cuf(nx,ny,nz,ng,nv,i_rank_start,pi,its1,its2,x0tr,y0tr,x0ts,y0ts,lamx,&
      &lamy,lams,lams1,phis,phis1,asl,bt,xc2_gpu,yc2_gpu,z_gpu,w_gpu,fl_gpu,wall_tag_gpu, dxdetanc2_gpu,&
      & dydetanc2_gpu)

    endassociate
  endsubroutine tripping

  subroutine debug_airfoil(self)
    class(equation_singleideal_gpu_object), intent(inout) :: self
    integer :: j,m,l,mm,ll,iv,ii
    real(rkind), parameter :: tol_check=1e-14
  endsubroutine debug_airfoil

  subroutine euler_x(self, istart, iend, istart_face, iend_face, do_update)
    class(equation_singleideal_gpu_object), intent(inout) :: self
    integer, intent(in) :: istart, iend
    integer :: lmax, weno_size, iercuda, ierror
    logical, optional :: do_update
    integer, optional :: istart_face, iend_face
    integer :: istart_face_, iend_face_
    logical :: do_update_
    type(dim3) :: grid, tblock
    integer :: force_zero_flux_min,force_zero_flux_max

    do_update_ = .true. ; if(present(do_update)) do_update_ = do_update
    istart_face_ = istart-1 ; if(present(istart_face)) istart_face_ = istart_face
    iend_face_ = iend ; if(present(iend_face)) iend_face_ = iend_face

    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv,&
    & nv_aux => self%nv_aux, ep_order => self%equation_base%ep_order, force_zero_flux => self%equation_ba&
    &se%force_zero_flux, coeff_deriv1_gpu => self%coeff_deriv1_gpu, dcsidx_gpu => self%base_gpu%dcsidx_gp&
    &u, fhat_gpu => self%fhat_gpu, w_gpu => self%base_gpu%w_gpu, w_aux_gpu => self%w_aux_gpu,&
    & fl_gpu => self%fl_gpu, sensor_threshold => self%equation_base%sensor_threshold, weno_scheme => self&
    &%equation_base%weno_scheme, weno_version => self%equation_base%weno_version, cp_coeff_gpu => self%cp&
    &_coeff_gpu, indx_cp_l => self%equation_base%indx_cp_l, indx_cp_r => self%equation_base%indx_cp_r,&
    & calorically_perfect => self%equation_base%calorically_perfect, ep_ord_change_gpu => self%ep_ord_cha&
    &nge_gpu, nkeep => self%equation_base%nkeep, flux_splitting => self%equation_base%flux_splitting,&
    & dcsidxc2_gpu => self%base_gpu%dcsidxc2_gpu, dcsidxnc2_gpu => self%base_gpu%dcsidxnc2_gpu,&
    & dcsidyc2_gpu => self%base_gpu%dcsidyc2_gpu, dcsidync2_gpu => self%base_gpu%dcsidync2_gpu,&
    & jac_gpu => self%base_gpu%jac_gpu, mcsijac1_gpu => self%base_gpu%mcsijac1_gpu, lmax_tag_gpu => self%&
    &lmax_tag_gpu, eul_imin => self%equation_base%eul_imin, eul_imax => self%equation_base%eul_imax,&
    & rgas0 => self%equation_base%rgas0)

      if(iend - istart >= 0) then
        weno_size = 2*weno_scheme
        lmax = ep_order/2
        force_zero_flux_min = force_zero_flux(1)
        force_zero_flux_max = force_zero_flux(2)
        tblock = dim3(eulerweno_threads_x,eulerweno_threads_y,1)

        grid = dim3(ceiling(real(iend_face_-istart_face_+1)/tblock%x),ceiling(real(ny)/tblock%y),1)

        if (self%grid%grid_dim == 1) then
          if (flux_splitting==1) then

            call euler_x_fluxes_hybrid_rusanov_kernel<<<grid, tblock, 0, stream1>>>(nv, nv_aux, nx,&
            & ny, nz, ng, istart_face_, iend_face_, lmax, nkeep, rgas0, w_aux_gpu, coeff_deriv1_gpu, dcsidx_gpu,&
            & fhat_gpu, force_zero_flux_min, force_zero_flux_max, weno_scheme, weno_version, sensor_threshold,&
            & weno_size, cp_coeff_gpu, indx_cp_l, indx_cp_r, ep_ord_change_gpu, calorically_perfect, tol_iter_nr,&
            &self%equation_base%rho0,self%equation_base%u0, self%equation_base%t0)
          else

            call euler_x_fluxes_hybrid_kernel<<<grid, tblock, 0, stream1>>>(nv, nv_aux, nx, ny, nz,&
            & ng, istart_face_, iend_face_, lmax, nkeep, rgas0, w_aux_gpu, coeff_deriv1_gpu, dcsidx_gpu,&
            & fhat_gpu, force_zero_flux_min, force_zero_flux_max, weno_scheme, weno_version, sensor_threshold,&
            & weno_size, cp_coeff_gpu, indx_cp_l, indx_cp_r, ep_ord_change_gpu, calorically_perfect, tol_iter_nr,&
            &self%equation_base%rho0,self%equation_base%u0, self%equation_base%t0)
          endif

          if(do_update_) then
            call euler_x_update_cuf(nx, ny, nz, ng, nv, eul_imin, eul_imax, fhat_gpu, fl_gpu, dcsidx_gpu, stream1)
          endif

        elseif (self%grid%grid_dim == 2) then
          call euler_x_fluxes_hybrid_c2_kernel<<<grid, tblock, 0, stream1>>>(nv, nv_aux, nx, ny, nz,&
          & ng, istart, iend, lmax, nkeep, rgas0, w_aux_gpu, coeff_deriv1_gpu, fhat_gpu, dcsidxc2_gpu,&
          & dcsidxnc2_gpu, dcsidyc2_gpu, dcsidync2_gpu, jac_gpu, mcsijac1_gpu, lmax_tag_gpu,&
          & force_zero_flux_min, force_zero_flux_max, weno_scheme, weno_version, sensor_threshold, weno_size,&
          & cp_coeff_gpu, indx_cp_l, indx_cp_r, ep_ord_change_gpu, calorically_perfect, tol_iter_nr,&
          &self%equation_base%rho0, self%equation_base%u0, self%equation_base%t0)
          call euler_x_update_c2_cuf(nx, ny, nz, ng, nv, istart, iend, fhat_gpu, fl_gpu, jac_gpu, stream1)
        endif

      endif

    endassociate
  endsubroutine euler_x

  subroutine euler_y(self, eul_jmin, eul_jmax)
    class(equation_singleideal_gpu_object), intent(inout) :: self
    integer, intent(in) :: eul_jmin, eul_jmax
    integer :: lmax, weno_size
    type(dim3) :: grid, tblock
    integer :: force_zero_flux_min,force_zero_flux_max

    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv,&
    & nv_aux => self%nv_aux, ep_order => self%equation_base%ep_order, force_zero_flux => self%equation_ba&
    &se%force_zero_flux, coeff_deriv1_gpu => self%coeff_deriv1_gpu, detady_gpu => self%base_gpu%detady_gp&
    &u, fhat_gpu => self%fhat_gpu, w_aux_gpu => self%w_aux_gpu, fl_gpu => self%fl_gpu,&
    & sensor_threshold => self%equation_base%sensor_threshold, weno_scheme => self%equation_base%weno_sch&
    &eme, weno_version => self%equation_base%weno_version, indx_cp_l => self%equation_base%indx_cp_l,&
    & indx_cp_r => self%equation_base%indx_cp_r, calorically_perfect => self%equation_base%calorically_pe&
    &rfect, cp_coeff_gpu => self%cp_coeff_gpu, ep_ord_change_gpu => self%ep_ord_change_gpu,&
    & nkeep => self%equation_base%nkeep, wall_tag_gpu => self%base_gpu%wall_tag_gpu, flux_splitting => se&
    &lf%equation_base%flux_splitting, detadxc2_gpu => self%base_gpu%detadxc2_gpu, detadxnc2_gpu => self%b&
    &ase_gpu%detadxnc2_gpu, detadyc2_gpu => self%base_gpu%detadyc2_gpu, detadync2_gpu => self%base_gpu%de&
    &tadync2_gpu, jac_gpu => self%base_gpu%jac_gpu, metajac1_gpu => self%base_gpu%metajac1_gpu,&
    & rgas0 => self%equation_base%rgas0)
      weno_size = 2*weno_scheme
      lmax = ep_order/2
      force_zero_flux_min = force_zero_flux(3)
      force_zero_flux_max = force_zero_flux(4)

      tblock = dim3(eulerweno_threads_x,eulerweno_threads_y,1)
      grid = dim3(ceiling(real(nx)/tblock%x),ceiling(real(nz)/tblock%y),1)
      if (self%grid%grid_dim == 1) then
        if (flux_splitting==1) then
          call euler_y_hybrid_rusanov_kernel<<<grid, tblock, 0, stream1>>>(nv, nv_aux, nx, ny, nz,&
          & ng, eul_jmin, eul_jmax, lmax, nkeep, rgas0, w_aux_gpu, fl_gpu, coeff_deriv1_gpu, detady_gpu,&
          & fhat_gpu, force_zero_flux_min, force_zero_flux_max, weno_scheme, weno_version, sensor_threshold,&
          & weno_size, cp_coeff_gpu, indx_cp_l, indx_cp_r, ep_ord_change_gpu, calorically_perfect, tol_iter_nr,&
          &self%equation_base%rho0,self%equation_base%u0, self%equation_base%t0)
        else
          call euler_y_hybrid_kernel<<<grid, tblock, 0, stream1>>>(nv, nv_aux, nx, ny, nz, ng,&
          & eul_jmin, eul_jmax, lmax, nkeep, rgas0, w_aux_gpu, fl_gpu, coeff_deriv1_gpu, detady_gpu, fhat_gpu,&
          & force_zero_flux_min, force_zero_flux_max, weno_scheme, weno_version, sensor_threshold, weno_size,&
          & cp_coeff_gpu, indx_cp_l, indx_cp_r, ep_ord_change_gpu, calorically_perfect, tol_iter_nr,&
          &self%equation_base%rho0,self%equation_base%u0, self%equation_base%t0)
        endif
        call euler_y_update_cuf(nx, ny, nz, ng, nv, eul_jmin, eul_jmax, fhat_gpu, fl_gpu, detady_gpu, stream1)
      elseif (self%grid%grid_dim == 2) then
        call euler_y_hybrid_c2_kernel<<<grid, tblock, 0, stream1>>>(nv, nv_aux, nx, ny, nz, ng,&
        & eul_jmin, eul_jmax, lmax, nkeep, rgas0, w_aux_gpu, fl_gpu, coeff_deriv1_gpu, fhat_gpu,&
        & detadxc2_gpu, detadxnc2_gpu, detadyc2_gpu, detadync2_gpu, jac_gpu, metajac1_gpu, wall_tag_gpu,&
        & force_zero_flux_min, force_zero_flux_max, weno_scheme, weno_version, sensor_threshold, weno_size,&
        & cp_coeff_gpu, indx_cp_l, indx_cp_r, ep_ord_change_gpu, calorically_perfect, tol_iter_nr,&
        &self%equation_base%rho0,self%equation_base%u0, self%equation_base%t0)
        call euler_y_update_c2_cuf(nx, ny, nz, ng, nv, eul_jmin, eul_jmax, fhat_gpu, fl_gpu, jac_gpu, wall_tag_gpu, stream1)
      endif
    endassociate
  endsubroutine euler_y

  subroutine euler_z(self, eul_kmin, eul_kmax)
    class(equation_singleideal_gpu_object), intent(inout) :: self
    integer, intent(in) :: eul_kmin, eul_kmax
    integer :: lmax, weno_size
    type(dim3) :: grid, tblock
    integer :: force_zero_flux_min, force_zero_flux_max

    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv,&
    & nv_aux => self%nv_aux, ep_order => self%equation_base%ep_order, force_zero_flux => self%equation_ba&
    &se%force_zero_flux, coeff_deriv1_gpu => self%coeff_deriv1_gpu, dzitdz_gpu => self%base_gpu%dzitdz_gp&
    &u, fhat_gpu => self%fhat_gpu, w_aux_gpu => self%w_aux_gpu, fl_gpu => self%fl_gpu,&
    & sensor_threshold => self%equation_base%sensor_threshold, weno_scheme => self%equation_base%weno_sch&
    &eme, weno_version => self%equation_base%weno_version, cp_coeff_gpu => self%cp_coeff_gpu,&
    & indx_cp_l => self%equation_base%indx_cp_l, indx_cp_r => self%equation_base%indx_cp_r,&
    & calorically_perfect => self%equation_base%calorically_perfect, ep_ord_change_gpu => self%ep_ord_cha&
    &nge_gpu, nkeep => self%equation_base%nkeep, flux_splitting => self%equation_base%flux_splitting,&
    & rgas0 => self%equation_base%rgas0, inflow_random_plane => self%equation_base%inflow_random_plane,&
    & inflow_random_plane_gpu => self%inflow_random_plane_gpu)
      weno_size = 2*weno_scheme
      lmax = ep_order/2
      force_zero_flux_min = force_zero_flux(5)
      force_zero_flux_max = force_zero_flux(6)
      tblock = dim3(eulerweno_threads_x,eulerweno_threads_y,1)
      grid = dim3(ceiling(real(nx)/tblock%x),ceiling(real(ny)/tblock%y),1)
      if (flux_splitting==1) then
        call euler_z_hybrid_rusanov_kernel<<<grid, tblock, 0, 0>>>(nv, nv_aux, nx, ny, nz, ng,&
        & eul_kmin, eul_kmax, lmax, nkeep, rgas0, w_aux_gpu, fl_gpu, coeff_deriv1_gpu, dzitdz_gpu, fhat_gpu,&
        & force_zero_flux_min, force_zero_flux_max, weno_scheme, weno_version, sensor_threshold, weno_size,&
        & cp_coeff_gpu, indx_cp_l, indx_cp_r, ep_ord_change_gpu, calorically_perfect, tol_iter_nr,&
        &self%equation_base%rho0,self%equation_base%u0, self%equation_base%t0)
      else
        call euler_z_hybrid_kernel<<<grid, tblock, 0, 0>>>(nv, nv_aux, nx, ny, nz, ng, eul_kmin,&
        & eul_kmax, lmax, nkeep, rgas0, w_aux_gpu, fl_gpu, coeff_deriv1_gpu, dzitdz_gpu, fhat_gpu,&
        & force_zero_flux_min, force_zero_flux_max, weno_scheme, weno_version, sensor_threshold, weno_size,&
        & cp_coeff_gpu, indx_cp_l, indx_cp_r, ep_ord_change_gpu, calorically_perfect, tol_iter_nr,&
        &self%equation_base%rho0,self%equation_base%u0, self%equation_base%t0)
      endif
      call euler_z_update_cuf(nx, ny, nz, ng, nv, eul_kmin, eul_kmax, fhat_gpu, fl_gpu, dzitdz_gpu, 0_cuda_stream_kind)
      if (self%equation_base%recyc) call get_crandom_f(inflow_random_plane(2:self%equation_base%jbl_inflow,1:nz,1:3))
    endassociate
  endsubroutine euler_z

  subroutine force_rhs(self)
    class(equation_singleideal_gpu_object), intent(inout) :: self
    real(rkind), dimension(5) :: bulk5, bulk5g

    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nxmax => self%equation_bas&
    &e%grid%nxmax, nzmax => self%equation_base%grid%nzmax, yn => self%equation_base%field%yn,&
    & yn_gpu => self%base_gpu%yn_gpu, w_gpu => self%base_gpu%w_gpu, w_aux_gpu => self%w_aux_gpu,&
    & fln_gpu => self%fln_gpu, volchan => self%equation_base%volchan, fluid_mask_gpu => self%fluid_mask_g&
    &pu, r_curv => self%grid%r_curv, volume => self%field%volume, dxdcsinc2_gpu => self%base_gpu%dxdcsinc&
    &2_gpu, dydcsinc2_gpu => self%base_gpu%dydcsinc2_gpu, dcsidxnc2_gpu => self%base_gpu%dcsidxnc2_gpu,&
    & dcsidync2_gpu => self%base_gpu%dcsidync2_gpu, jac_gpu => self%base_gpu%jac_gpu)

      if(self%grid%grid_dim == 1) then
        call force_rhs_1_cuf(nx, ny, nz, ng, yn_gpu, fln_gpu, w_gpu, w_aux_gpu, bulk5, fluid_mask_gpu)
      elseif(self%grid%grid_dim == 2) then
        call force_rhs_1_c2_cuf(nx, ny, nz, ng, yn_gpu, fln_gpu, w_gpu, w_aux_gpu, bulk5,&
        & fluid_mask_gpu, dcsidxnc2_gpu, dcsidync2_gpu, jac_gpu)
      endif


      call mpi_allreduce(bulk5,bulk5g,5,mpi_prec,mpi_sum,self%equation_base%field%mp_cart,self%mpi_err)

      if(self%grid%grid_dim == 1) then
        bulk5g = bulk5g/volchan
        self%equation_base%dpdx = self%equation_base%dpdx+bulk5g(2)
        self%equation_base%rhobulk = bulk5g(3)
        self%equation_base%ubulk = bulk5g(4)/self%equation_base%rhobulk
        self%equation_base%tbulk = bulk5g(5)/self%equation_base%rhobulk/self%equation_base%ubulk
      elseif(self%grid%grid_dim == 2) then
        bulk5g = bulk5g/volume
        bulk5g(5) = bulk5g(5)*r_curv
        self%equation_base%rhobulk = bulk5g(2)
        self%equation_base%ubulk = bulk5g(3)/self%equation_base%rhobulk
        self%equation_base%tbulk = bulk5g(4)/self%equation_base%rhobulk/self%equation_base%ubulk
        self%equation_base%dpth = bulk5g(5)
      endif

      if(self%grid%grid_dim == 1) then
        call force_rhs_2_cuf(nx, ny, nz, ng, fln_gpu, w_aux_gpu, bulk5g, fluid_mask_gpu)
      elseif(self%grid%grid_dim == 2) then
        call force_rhs_2_c2_cuf(nx, ny, nz, ng, fln_gpu, w_aux_gpu, bulk5g, fluid_mask_gpu, r_curv,&
        & yn_gpu, dxdcsinc2_gpu, dydcsinc2_gpu)
      endif
    endassociate
  endsubroutine force_rhs
  subroutine force_var(self)
    class(equation_singleideal_gpu_object), intent(inout) :: self
    real(rkind) :: bulkt, bulktg, tbtarget, tbdiff, trelax

    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nxmax => self%equation_bas&
    &e%grid%nxmax, nzmax => self%equation_base%grid%nzmax, yn => self%equation_base%field%yn,&
    & yn_gpu => self%base_gpu%yn_gpu, w_gpu => self%base_gpu%w_gpu, w_aux_gpu => self%w_aux_gpu,&
    & fln_gpu => self%fln_gpu, volchan => self%equation_base%volchan, fluid_mask_gpu => self%fluid_mask_g&
    &pu, cv_coeff_gpu => self%cv_coeff_gpu, volume => self%field%volume, t0 => self%equation_base%t0,&
    & time => self%equation_base%time, calorically_perfect => self%equation_base%calorically_perfect,&
    & jac_gpu => self%base_gpu%jac_gpu, dcsidxnc2_gpu => self%base_gpu%dcsidxnc2_gpu, dcsidync2_gpu => se&
    &lf%base_gpu%dcsidync2_gpu, indx_cp_l => self%equation_base%indx_cp_l, indx_cp_r => self%equation_bas&
    &e%indx_cp_r)

      if (self%equation_base%theta_wall>=-1._rkind) then
        if(self%grid%grid_dim == 1) then
          call force_var_1_cuf(nx, ny, nz, ng, yn_gpu, fln_gpu, w_gpu, w_aux_gpu, bulkt,&
          & fluid_mask_gpu, cv_coeff_gpu, indx_cp_l, indx_cp_r, t0, calorically_perfect, tol_iter_nr)
          call mpi_allreduce(bulkt,bulktg,1,mpi_prec,mpi_sum,self%equation_base%field%mp_cart,self%mpi_err)

          bulktg = bulktg/volchan
          bulktg = bulktg/self%equation_base%rhobulk/self%equation_base%ubulk
          tbtarget = self%equation_base%t_bulk_target
          tbdiff = tbtarget-bulktg
          call force_var_2_cuf(nx, ny, nz, ng, w_gpu, w_aux_gpu, tbdiff, fluid_mask_gpu,&
          & cv_coeff_gpu, indx_cp_l, indx_cp_r, t0, calorically_perfect)
        elseif(self%grid%grid_dim == 2) then
          call force_var_1_c2_cuf(nx, ny, nz, ng, yn_gpu, fln_gpu, w_gpu, w_aux_gpu, bulkt,&
          & fluid_mask_gpu, cv_coeff_gpu, indx_cp_l, indx_cp_r, t0, calorically_perfect, tol_iter_nr, jac_gpu,&
          & dcsidxnc2_gpu, dcsidync2_gpu)
          call mpi_allreduce(bulkt,bulktg,1,mpi_prec,mpi_sum,self%equation_base%field%mp_cart,self%mpi_err)

          trelax = min(1._rkind, time/10._rkind)
          bulktg = bulktg/volume
          bulktg = bulktg/self%equation_base%rhobulk/self%equation_base%ubulk
          tbtarget = trelax*self%equation_base%t_bulk_target+(1._rkind-trelax)*bulktg
          tbdiff = tbtarget-bulktg

          call force_var_2_c2_cuf(nx, ny, nz, ng, w_gpu, w_aux_gpu, tbdiff, fluid_mask_gpu,&
          & cv_coeff_gpu, indx_cp_l, indx_cp_r, t0, calorically_perfect, jac_gpu)
        endif
      endif
    endassociate
  endsubroutine force_var

  subroutine recyc_exchange(self)
    class(equation_singleideal_gpu_object), intent(inout) :: self

    integer :: indx
    integer, dimension(4) :: req
    integer :: kcoordsendto1,kcoordsendto2
    integer :: kcoordrecvfrom1,kcoordrecvfrom2
    integer :: sendto1,sendto2
    integer :: recvfrom1,recvfrom2
    integer :: kshiftglob
    integer :: n1_start_send,n1_end_send,n2_start_send,n2_end_send
    integer :: n1_start_recv,n1_end_recv,n2_start_recv,n2_end_recv
    integer :: iercuda

    associate(nzmax => self%equation_base%grid%nzmax, ng => self%equation_base%grid%ng,&
    & nx => self%equation_base%field%nx, ny => self%equation_base%field%ny, nz => self%equation_base%fiel&
    &d%nz, ncoords => self%equation_base%field%ncoords, mp_cart => self%equation_base%field%mp_cart,&
    & nblocks => self%equation_base%field%nblocks, iermpi => self%mpi_err, w_gpu => self%base_gpu%w_gpu,&
    & nv => self%equation_base%nv, wbuf1s_gpu => self%base_gpu%wbuf1s_gpu, wbuf2r_gpu => self%base_gpu%wb&
    &uf2r_gpu, wbuf1r_gpu => self%base_gpu%wbuf1r_gpu, wrecyc_gpu => self%wrecyc_gpu, ibrecyc => self%equ&
    &ation_base%ib_recyc, irecyc => self%equation_base%i_recyc)

      kshiftglob = nzmax/2
      n1_start_send = 1
      n1_end_send = nz-mod(kshiftglob,nz)
      n2_start_send = n1_end_send+1
      n2_end_send = nz
      n1_start_recv = 1+mod(kshiftglob,nz)
      n1_end_recv = nz
      n2_start_recv = 1
      n2_end_recv = mod(kshiftglob,nz)

      req = mpi_request_null

      if (ncoords(1)==ibrecyc) then
        kcoordsendto1 = ncoords(3)+kshiftglob/nz
        kcoordsendto2 = kcoordsendto1+1
        kcoordsendto1 = mod(kcoordsendto1,nblocks(3))
        kcoordsendto2 = mod(kcoordsendto2,nblocks(3))
        call mpi_cart_rank(mp_cart,[0,0,kcoordsendto1],sendto1,iermpi)
        call mpi_cart_rank(mp_cart,[0,0,kcoordsendto2],sendto2,iermpi)
        call recyc_exchange_cuf_1(irecyc, w_gpu, wbuf1s_gpu, nx, ny, nz, ng, nv)
        indx = nv*ng*ny*nz
        call mpi_isend(wbuf1s_gpu,indx,mpi_prec,sendto1,2000,mp_cart,req(1),iermpi)
        call mpi_isend(wbuf1s_gpu,indx,mpi_prec,sendto2,3000,mp_cart,req(2),iermpi)
      endif
      if (ncoords(1)==0) then
        kcoordrecvfrom1 = ncoords(3)-kshiftglob/nz+nblocks(3)
        kcoordrecvfrom2 = kcoordrecvfrom1-1
        kcoordrecvfrom1 = mod(kcoordrecvfrom1,nblocks(3))
        kcoordrecvfrom2 = mod(kcoordrecvfrom2,nblocks(3))
        call mpi_cart_rank(mp_cart,[ibrecyc,0,kcoordrecvfrom1],recvfrom1,iermpi)
        call mpi_cart_rank(mp_cart,[ibrecyc,0,kcoordrecvfrom2],recvfrom2,iermpi)
        indx = nv*ng*ny*nz
        call mpi_irecv(wbuf1r_gpu,indx,mpi_prec,recvfrom1,2000,mp_cart,req(3),iermpi)
        call mpi_irecv(wbuf2r_gpu,indx,mpi_prec,recvfrom2,3000,mp_cart,req(4),iermpi)
      endif
      call mpi_waitall(4,req,mpi_statuses_ignore,iermpi)
      if (ncoords(1)==0) then
        call recyc_exchange_cuf_2(n1_start_recv, n1_start_send, n1_end_recv, wrecyc_gpu, wbuf1r_gpu, nx, ny, nz, ng, nv)
        call recyc_exchange_cuf_3(n2_start_recv, n2_start_send, n2_end_recv, wrecyc_gpu, wbuf2r_gpu, nx, ny, nz, ng, nv)
      endif
    endassociate

  end subroutine recyc_exchange

  subroutine bc_nr(self)
    class(equation_singleideal_gpu_object), intent(inout) :: self
    integer :: ilat, dir, start_or_end
    type(dim3) :: grid, tblock

    associate(bctags_nr => self%equation_base%bctags_nr, nx => self%nx, ny => self%ny,&
    & nz => self%nz, nv => self%nv, ng => self%ng, w_aux_gpu => self%w_aux_gpu, w_gpu => self%base_gpu%w_&
    &gpu, fl_gpu => self%fl_gpu, dcsidx_gpu => self%base_gpu%dcsidx_gpu, detady_gpu => self%base_gpu%deta&
    &dy_gpu, jac_gpu => self%base_gpu%jac_gpu , dcsidxnc2_gpu => self%base_gpu%dcsidxnc2_gpu ,&
    & dcsidync2_gpu => self%base_gpu%dcsidync2_gpu , dcsidxc2_gpu => self%base_gpu%dcsidxc2_gpu ,&
    & dcsidyc2_gpu => self%base_gpu%dcsidyc2_gpu , detadxnc2_gpu => self%base_gpu%detadxnc2_gpu ,&
    & detadync2_gpu => self%base_gpu%detadync2_gpu , detadxc2_gpu => self%base_gpu%detadxc2_gpu ,&
    & detadyc2_gpu => self%base_gpu%detadyc2_gpu , wall_tag_gpu => self%base_gpu%wall_tag_gpu ,&
    & dzitdz_gpu => self%base_gpu%dzitdz_gpu, wfar_gpu => self%wfar_gpu, mcsi_gpu => self%base_gpu%mcsi_g&
    &pu, meta_gpu => self%base_gpu%meta_gpu, indx_cp_l => self%equation_base%indx_cp_l,&
    & indx_cp_r => self%equation_base%indx_cp_r, calorically_perfect => self%equation_base%calorically_pe&
    &rfect, rgas0 => self%equation_base%rgas0, t0 => self%equation_base%t0, cp_coeff_gpu => self%cp_coeff&
    &_gpu, winf_gpu => self%winf_gpu)

      do ilat=1,2*self%grid%ndim
        if(bctags_nr(ilat) > 0) then
          dir = (ilat-1)/2 +1
          start_or_end = mod(ilat-1,2)+1
          if(dir == 1) then
            tblock = dim3(eulerweno_threads_x,eulerweno_threads_y,1)
            grid = dim3(ceiling(real(ny)/tblock%x),ceiling(real(nz)/tblock%y),1)
            if(self%grid%grid_dim == 1) then
              call bc_nr_lat_x_kernel<<<grid, tblock, 0, 0>>>(start_or_end, bctags_nr(ilat), nx, ny,&
              & nz, ng, nv, w_aux_gpu, w_gpu, fl_gpu, dcsidx_gpu, indx_cp_l, indx_cp_r, cp_coeff_gpu, winf_gpu,&
              & calorically_perfect,rgas0,t0)
            elseif(self%grid%grid_dim == 2) then
              call bc_nr_lat_x_c2_kernel<<<grid, tblock, 0, 0>>>(start_or_end, bctags_nr(ilat), nx,&
              & ny, nz, ng, nv, w_aux_gpu, w_gpu, fl_gpu, dcsidx_gpu, indx_cp_l, indx_cp_r, cp_coeff_gpu, winf_gpu,&
              & jac_gpu, mcsi_gpu, dcsidxnc2_gpu, dcsidync2_gpu, dcsidxc2_gpu, dcsidyc2_gpu, calorically_perfect,&
              &rgas0,t0)
            endif
          endif
          if(dir == 2) then
            tblock = dim3(eulerweno_threads_x,eulerweno_threads_y,1)
            grid = dim3(ceiling(real(nx)/tblock%x),ceiling(real(nz)/tblock%y),1)
            if(self%grid%grid_dim == 1) then
              call bc_nr_lat_y_kernel<<<grid, tblock, 0, 0>>>(start_or_end, bctags_nr(ilat), nx, ny,&
              & nz, ng, nv, w_aux_gpu, w_gpu, fl_gpu, detady_gpu, indx_cp_l, indx_cp_r, cp_coeff_gpu, winf_gpu,&
              & calorically_perfect,rgas0,t0)
            elseif(self%grid%grid_dim == 2) then
              call bc_nr_lat_y_c2_kernel<<<grid, tblock, 0, 0>>>(start_or_end, bctags_nr(ilat), nx,&
              & ny, nz, ng, nv, w_aux_gpu, w_gpu, fl_gpu, detady_gpu, indx_cp_l, indx_cp_r, cp_coeff_gpu, wfar_gpu,&
              & jac_gpu, meta_gpu, detadxnc2_gpu, detadync2_gpu, detadxc2_gpu, detadyc2_gpu, wall_tag_gpu,&
              & calorically_perfect,rgas0,t0)
            endif
          endif
          if(dir == 3) then
            tblock = dim3(eulerweno_threads_x,eulerweno_threads_y,1)
            grid = dim3(ceiling(real(nx)/tblock%x),ceiling(real(ny)/tblock%y),1)
            call bc_nr_lat_z_kernel<<<grid, tblock, 0, 0>>>(start_or_end, bctags_nr(ilat), nx, ny,&
            & nz, ng, nv, w_aux_gpu, w_gpu, fl_gpu, dzitdz_gpu, indx_cp_l, indx_cp_r, cp_coeff_gpu, winf_gpu,&
            & calorically_perfect,rgas0,t0)
          endif
        endif
      enddo
    endassociate

  endsubroutine bc_nr

  subroutine update_ghost(self, do_swap)
    class(equation_singleideal_gpu_object), intent(inout) :: self
    integer, intent(in), optional :: do_swap
    integer :: do_swap_
    integer :: ilat

    do_swap_ = 1 ; if (present(do_swap)) do_swap_ = do_swap

    if (self%equation_base%recyc) call self%recyc_exchange()

    do ilat=1,2*self%grid%ndim
      select case(self%equation_base%bctags(ilat))
      case(0)
      case(1)
        call bcfree_cuf(ilat, self%nx, self%ny, self%nz, self%ng, self%nv, self%winf_gpu, self%base_gpu%w_gpu)
      case(3)
        call bcfree_sub_cuf(ilat, self%nx, self%ny, self%nz, self%ng, self%base_gpu%w_gpu,&
        & self%w_aux_gpu, self%equation_base%aoa, self%equation_base%t0, self%equation_base%ptot0,&
        & self%equation_base%ttot0, self%equation_base%rgas0, self%equation_base%indx_cp_l,&
        & self%equation_base%indx_cp_r, self%cv_coeff_gpu, self%cp_coeff_gpu, self%equation_base%calorically_&
        &perfect, tol_iter_nr)
      case(2)
        call bcextr_cuf(ilat, self%nx, self%ny, self%nz, self%ng, self%nv, self%base_gpu%w_gpu)
      case(4)
        call bcextr_sub_cuf(ilat, self%nx, self%ny, self%nz, self%ng, self%equation_base%p0,&
        & self%equation_base%rgas0, self%base_gpu%w_gpu, self%equation_base%indx_cp_l, self%equation_base%ind&
        &x_cp_r, self%cv_coeff_gpu, self%equation_base%t0, self%equation_base%calorically_perfect)
      case(5)
        if(self%grid%grid_dim == 1) then
          call bcsym_cuf(ilat, self%nx, self%ny, self%nz, self%ng, self%base_gpu%w_gpu)
        elseif(self%grid%grid_dim == 2) then
          call bcsym_c2_cuf(ilat, self%nx, self%ny, self%nz, self%ng, self%base_gpu%w_gpu,&
          & self%base_gpu%dxdcsic2_gpu, self%base_gpu%dydcsic2_gpu)
        endif
      case(6)
        if (self%grid%is_y_staggered) then
          call bcwall_staggered_cuf(ilat, self%nx, self%ny, self%nz, self%ng, self%equation_base%t_w&
          &all, self%base_gpu%w_gpu, self%w_aux_gpu, self%equation_base%indx_cp_l, self%equation_base%indx_cp_r&
          &, self%cv_coeff_gpu, self%equation_base%t0, self%equation_base%rgas0, self%equation_base%calorically&
          &_perfect,tol_iter_nr)
        else
          call bcwall_cuf(ilat, self%nx, self%ny, self%nz, self%ng, self%equation_base%t_wall,&
          & self%base_gpu%w_gpu, self%w_aux_gpu, self%equation_base%indx_cp_l, self%equation_base%indx_cp_r,&
          & self%cv_coeff_gpu, self%equation_base%t0, self%equation_base%rgas0, self%equation_base%calorically_&
          &perfect,tol_iter_nr)
        endif
      case(7)
        call bcshock_cuf(ilat, self%nx, self%ny, self%nz, self%ng, self%nv, self%base_gpu%w_gpu,&
        & self%winf_gpu, self%winf_past_shock_gpu, self%equation_base%xshock_imp, self%equation_base%shock_an&
        &gle, self%base_gpu%x_gpu, self%base_gpu%y_gpu, self%equation_base%tanhfacs)
      case(8)
      case(9)
        call bclam_cuf(ilat, self%nx, self%ny, self%nz, self%ng, self%nv, self%base_gpu%w_gpu,&
        & self%wmean_gpu, self%equation_base%p0, self%equation_base%rgas0, self%equation_base%indx_cp_l,&
        & self%equation_base%indx_cp_r, self%cv_coeff_gpu, self%equation_base%t0, self%equation_base%calorica&
        &lly_perfect)
      case(10)
        call self%bcrecyc(ilat)
      case(11)
        call self%base_gpu%bcswap_wake()
        call bcwall_airfoil_cuf(ilat, self%nx, self%ny, self%nz, self%ng, self%equation_base%t_wall,&
        & self%base_gpu%w_gpu, self%w_aux_gpu, self%base_gpu%wall_tag_gpu, self%equation_base%indx_cp_l,&
        & self%equation_base%indx_cp_r, self%cv_coeff_gpu, self%equation_base%t0, self%equation_base%rgas0,&
        & self%equation_base%calorically_perfect,tol_iter_nr, self%base_gpu%xc2_gpu, self%base_gpu%yc2_gpu,&
        & self%equation_base%a_tw, self%equation_base%v_bs, self%equation_base%thic, self%equation_base%kx_tw&
        &, self%equation_base%om_tw, self%equation_base%xtw1, self%equation_base%xtw2, self%equation_base%tim&
        &e, self%base_gpu%dxdetanc2_gpu, self%base_gpu%dydetanc2_gpu, self%equation_base%u0)
        call self%base_gpu%bcte()
      endselect
    enddo

    if (do_swap_ == 1) call self%base_gpu%bcswap()

  endsubroutine update_ghost

  subroutine bcrecyc(self, ilat)
    class(equation_singleideal_gpu_object), intent(inout) :: self
    integer, intent(in) :: ilat
    integer :: ntot
    if (ilat == 1) then
      associate(nx => self%nx, ny => self%ny, nz => self%nz, nzmax => self%grid%nzmax,&
      & ng => self%ng, nv => self%nv, wrecycav_gpu => self%wrecycav_gpu, wrecyc_gpu => self%wrecyc_gpu,&
      & mp_cartz => self%field%mp_cartz, iermpi => self%mpi_err, p0 => self%equation_base%p0,&
      & w_gpu => self%base_gpu%w_gpu, wmean_gpu => self%wmean_gpu, weta_inflow_gpu => self%weta_inflow_gpu,&
      & map_j_inn_gpu => self%map_j_inn_gpu, map_j_out_gpu => self%map_j_out_gpu, map_j_out_blend_gpu => se&
      &lf%map_j_out_blend_gpu, eta_recyc_blend_gpu => self%eta_recyc_blend_gpu, yplus_inflow_gpu => self%yp&
      &lus_inflow_gpu, eta_inflow_gpu => self%eta_inflow_gpu, yplus_recyc_gpu => self%yplus_recyc_gpu,&
      & eta_recyc_gpu => self%eta_recyc_gpu, betarecyc => self%equation_base%betarecyc, i_recyc => self%equ&
      &ation_base%i_recyc, glund1 => self%equation_base%glund1, indx_cp_l => self%equation_base%indx_cp_l,&
      & indx_cp_r => self%equation_base%indx_cp_r, rgas0 => self%equation_base%rgas0, calorically_perfect =&
      &> self%equation_base%calorically_perfect, rand_type => self%equation_base%rand_type,&
      & t0 => self%equation_base%t0, u0 => self%equation_base%u0, l0 => self%equation_base%l0,&
      & cv_coeff_gpu => self%cv_coeff_gpu, cp_coeff_gpu => self%cp_coeff_gpu, inflow_random_plane => self%e&
      &quation_base%inflow_random_plane, inflow_random_plane_gpu => self%inflow_random_plane_gpu)
        call bcrecyc_cuf_1(nx, ny, nz, ng, nv, wrecycav_gpu, wrecyc_gpu)

        inflow_random_plane_gpu = inflow_random_plane

        ntot = ng*ny*nv
        call mpi_allreduce(mpi_in_place,wrecycav_gpu,ntot,mpi_prec,mpi_sum,mp_cartz,iermpi)

        call bcrecyc_cuf_2(nx, ny, nz, nzmax, ng, wrecycav_gpu, wrecyc_gpu)

        call bcrecyc_cuf_3(nx, ny, nz, ng, p0, u0, rgas0, w_gpu, wmean_gpu, wrecyc_gpu,&
        & weta_inflow_gpu, map_j_inn_gpu, map_j_out_gpu, map_j_out_blend_gpu, yplus_inflow_gpu,&
        & eta_inflow_gpu, yplus_recyc_gpu, eta_recyc_gpu, eta_recyc_blend_gpu, betarecyc, glund1,&
        & inflow_random_plane_gpu, indx_cp_l, indx_cp_r, cv_coeff_gpu, cp_coeff_gpu, t0, calorically_perfect,&
        &rand_type)
      endassociate
    endif

  end subroutine bcrecyc

  subroutine initialize(self, filename)
    class(equation_singleideal_gpu_object), intent(inout) :: self
    character(*) , intent(in) :: filename

    call get_mpi_basic_info(self%nprocs, self%myrank, self%masterproc, self%mpi_err)

    call self%equation_base%initialize(filename)
    self%nx = self%equation_base%field%nx
    self%ny = self%equation_base%field%ny
    self%nz = self%equation_base%field%nz
    self%ng = self%equation_base%grid%ng
    self%nv = self%equation_base%nv
    self%nv_aux = self%equation_base%nv_aux

    self%visc_model = self%equation_base%visc_model
    self%mu0 = self%equation_base%mu0
    self%t_ref_dim = self%equation_base%t_ref_dim
    self%sutherland_s = self%equation_base%sutherland_s
    self%powerlaw_vtexp = self%equation_base%powerlaw_vtexp

    call self%point_to_field(self%equation_base%field)
    call self%point_to_grid(self%equation_base%grid)
    self%num_iter = self%equation_base%num_iter
    self%time0 = self%equation_base%time0
    self%icyc0 = self%equation_base%icyc0

    call self%base_gpu%initialize(self%equation_base%field, self%equation_base%gpu_bind)

    if(self%equation_base%debug_memory>0) then
      call self%field%check_cpu_mem(description="--initialize-start--")
      call self%base_gpu%check_gpu_mem(description="--initialize-start--")
    endif

    call self%base_gpu%copy_cpu_gpu()

    if(self%equation_base%debug_memory>0) then
      call self%field%check_cpu_mem(description="--initialize-first-gpu-usage--")
      call self%base_gpu%check_gpu_mem(description="--initialize-first-gpu-usage--")
    endif

    call self%alloc()

    self%w_aux_gpu = self%equation_base%w_aux
    self%fluid_mask_gpu = self%equation_base%fluid_mask
    self%ep_ord_change_gpu = self%equation_base%ep_ord_change
    self%winf_gpu = self%equation_base%winf
    self%winf_past_shock_gpu = self%equation_base%winf_past_shock
    if (self%grid%grid_dim == 2) self%wfar_gpu = self%equation_base%wfar

    self%lmax_tag_gpu = self%equation_base%lmax_tag
    self%vis_tag_gpu = self%equation_base%vis_tag

    self%cp_coeff_gpu = self%equation_base%cp_coeff
    self%cv_coeff_gpu = self%equation_base%cv_coeff

    if (self%equation_base%num_probe>0) then
      self%probe_coeff_gpu = self%equation_base%probe_coeff
      self%ijk_probe_gpu = self%equation_base%ijk_probe
    endif

    if(self%equation_base%flow_init /= 5) then
      self%wmean_gpu = self%equation_base%wmean
    endif

    if(self%equation_base%enable_tspec > 0) then
      self%w_tspec_gpu = self%equation_base%w_tspec
    endif

    allocate(self%coeff_deriv1_gpu(1:4,4))
    allocate(self%coeff_deriv2_gpu(0:4,4))
    self%coeff_deriv1_gpu = self%equation_base%coeff_deriv1
    self%coeff_deriv2_gpu = self%equation_base%coeff_deriv2

    if (self%equation_base%ifilter > 0) then
      self%coeff_filter_gpu = self%equation_base%coeff_filter
    endif

    if (self%equation_base%enable_sponge > 0) then
      self%f_sponge_gpu = self%equation_base%f_sponge
    endif

    if (self%equation_base%recyc) then
      self%yplus_inflow_gpu = self%equation_base%yplus_inflow
      self%eta_inflow_gpu = self%equation_base%eta_inflow
      self%yplus_recyc_gpu = self%equation_base%yplus_recyc
      self%eta_recyc_gpu = self%equation_base%eta_recyc
      self%eta_recyc_blend_gpu = self%equation_base%eta_recyc_blend
      self%map_j_inn_gpu = self%equation_base%map_j_inn
      self%map_j_out_gpu = self%equation_base%map_j_out
      self%map_j_out_blend_gpu = self%equation_base%map_j_out_blend
      self%weta_inflow_gpu = self%equation_base%weta_inflow
    endif


    if (self%equation_base%enable_insitu > 0) then
      call self%insitu_alloc_gpu()
    endif
    if (self%equation_base%debug_memory>0) then
      call self%field%check_cpu_mem(description="initialize-completed")
      call self%base_gpu%check_gpu_mem(description="initialize-completed")
    endif


    self%mpi_err = cudastreamcreate(stream1)
  endsubroutine initialize

  subroutine alloc(self)
    class(equation_singleideal_gpu_object), intent(inout) :: self
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv,&
    & nv_aux => self%nv_aux, weno_scheme => self%equation_base%weno_scheme,calorically_perfect => self%eq&
    &uation_base%calorically_perfect, indx_cp_l => self%equation_base%indx_cp_l ,indx_cp_r => self%equati&
    &on_base%indx_cp_r, ndft => self%equation_base%ndft_tspec, grid_dim => self%grid%grid_dim)

      allocate(self%w_aux_gpu(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng, nv_aux))
      allocate(self%fluid_mask_gpu(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng))
      allocate(self%ep_ord_change_gpu(0:nx, 0:ny, 0:nz, 1:3))

      self%w_aux_gpu = self%equation_base%w_aux
      allocate(self%winf_gpu(nv))
      allocate(self%winf_past_shock_gpu(nv))
      if (grid_dim == 2) allocate (self%wfar_gpu(nx,nv))


      allocate(self%fl_gpu(1:nx, 1:ny, 1:nz, nv))
      allocate(self%fln_gpu(1:nx, 1:ny, 1:nz, nv))
      allocate(self%wallprop_gpu(1-ng:nx+ng, 1-ng:nz+ng, 2:4))
      allocate(self%fhat_gpu(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng, nv))

      allocate(self%wrecyc_gpu(ng,ny,nz,nv))
      allocate(self%wrecycav_gpu(ng,ny,nv))
      allocate(self%wmean_gpu(1-self%grid%ng:self%field%nx+self%grid%ng+1, 1:self%field%ny,4))

      allocate(self%yplus_inflow_gpu(1-ng:ny+ng))
      allocate(self%eta_inflow_gpu(1-ng:ny+ng))
      allocate(self%yplus_recyc_gpu(1-ng:ny+ng))
      allocate(self%eta_recyc_gpu(1-ng:ny+ng))
      allocate(self%eta_recyc_blend_gpu(1-ng:ny+ng))
      allocate(self%map_j_inn_gpu(1:ny))
      allocate(self%map_j_out_gpu(1:ny))
      allocate(self%map_j_out_blend_gpu(1:ny))
      allocate(self%weta_inflow_gpu(1:ny))
      allocate(self%inflow_random_plane_gpu(1:ny,1:nz,3))

      allocate(self%cv_coeff_gpu(indx_cp_l:indx_cp_r+1))
      allocate(self%cp_coeff_gpu(indx_cp_l:indx_cp_r+1))

      if (self%equation_base%num_probe>0) then
        allocate(self%w_aux_probe_gpu(6,self%equation_base%num_probe))
        allocate(self%ijk_probe_gpu(3,self%equation_base%num_probe))
        allocate(self%probe_coeff_gpu(2,2,2,self%equation_base%num_probe))
      endif

      allocate(self%lmax_tag_gpu(1-ng:nx+ng))
      allocate(self%vis_tag_gpu(1-ng:nx+ng))

      if (self%equation_base%enable_tspec > 0) then
        allocate(self%w_psd_tspec_gpu(nx, ndft))
        allocate(self%w_tspec_gpu(nx, nz, ndft, 2, 2))
      endif

      if (self%equation_base%ifilter > 0) allocate(self%coeff_filter_gpu(0:4))
      if (self%equation_base%enable_sponge > 0) allocate(self%f_sponge_gpu(ny))

    endassociate
  endsubroutine alloc

  subroutine compute_residual(self)
    class(equation_singleideal_gpu_object), intent(inout) :: self
    real(rkind) :: restemp
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv,&
    & dt => self%equation_base%dt, fln_gpu => self%fln_gpu, residual_rhou => self%equation_base%residual_&
    &rhou, fluid_mask_gpu => self%fluid_mask_gpu, vmax => self%equation_base%vmax, w_aux_gpu => self%w_au&
    &x_gpu, rgas0 => self%equation_base%rgas0, rhomin => self%equation_base%rhomin, rhomax => self%equati&
    &on_base%rhomax, pmin => self%equation_base%pmin, pmax => self%equation_base%pmax,&
    & tmin => self%equation_base%tmin, tmax => self%equation_base%tmax)

      call compute_residual_cuf(nx, ny, nz, ng, nv, fln_gpu, dt, residual_rhou, fluid_mask_gpu)
      call mpi_allreduce(mpi_in_place,residual_rhou,1,mpi_prec,mpi_sum,self%equation_base%field%mp_cart,self%mpi_err)
      restemp = real(nx,rkind)*real(ny,rkind)*real(nz,rkind)*real(self%nprocs,rkind)
      residual_rhou = sqrt(residual_rhou/restemp)

      call compute_rho_t_p_minmax_cuf(nx, ny, nz, ng, rgas0, w_aux_gpu,rhomin,rhomax, tmin, tmax, pmin, pmax, fluid_mask_gpu)
      call mpi_allreduce(mpi_in_place,rhomin,1,mpi_prec,mpi_min,self%equation_base%field%mp_cart,self%mpi_err)
      call mpi_allreduce(mpi_in_place,rhomax,1,mpi_prec,mpi_max,self%equation_base%field%mp_cart,self%mpi_err)
      call mpi_allreduce(mpi_in_place,tmin,1,mpi_prec,mpi_min,self%equation_base%field%mp_cart,self%mpi_err)
      call mpi_allreduce(mpi_in_place,tmax,1,mpi_prec,mpi_max,self%equation_base%field%mp_cart,self%mpi_err)
      call mpi_allreduce(mpi_in_place,pmin,1,mpi_prec,mpi_min,self%equation_base%field%mp_cart,self%mpi_err)
      call mpi_allreduce(mpi_in_place,pmax,1,mpi_prec,mpi_max,self%equation_base%field%mp_cart,self%mpi_err)

    endassociate
  endsubroutine compute_residual

  subroutine analyze_runtime(self)
    class(equation_singleideal_gpu_object), intent(inout) :: self
    associate(flow_init => self%equation_base%flow_init, enable_forces_runtime => self%equation_base&
    &%enable_forces_runtime , grid_dim => self%grid%grid_dim , count_weno_control => self%equation_base%c&
    &ount_weno_control , icyc => self%equation_base%icyc , icyc0 => self%equation_base%icyc0,&
    & ifilter => self%equation_base%ifilter)

      if(grid_dim == 2 .and. enable_forces_runtime > 0) then
        call self%compute_airfoil_forces_runtime()
      endif
      if (count_weno_control >0) then
        if (mod(icyc-icyc0,count_weno_control) == 0) then
          call self%count_weno()
        endif
      endif
      if (ifilter > 0 .and. grid_dim == 2) then
        if (mod(icyc-icyc0,ifilter) == 0) then
          call self%filter()
        endif
      endif
    endassociate
  endsubroutine analyze_runtime

  subroutine count_weno(self)
    class(equation_singleideal_gpu_object), intent(inout) :: self
    real(rkind) :: ntot, count_weno_x_frac, count_weno_y_frac, count_weno_z_frac
    integer :: lmax
    associate(nv_aux => self%equation_base%nv_aux, nx => self%field%nx, ny => self%field%ny,&
    & nz => self%field%nz, ng => self%grid%ng, nv => self%field%nv, nprocs => self%nprocs,&
    & eul_imin => self%equation_base%eul_imin, eul_imax => self%equation_base%eul_imax,&
    & eul_jmin => self%equation_base%eul_jmin, eul_jmax => self%equation_base%eul_jmax,&
    & eul_kmin => self%equation_base%eul_kmin, eul_kmax => self%equation_base%eul_kmax,&
    & ep_order => self%equation_base%ep_order, grid_dim => self%grid%grid_dim , weno_scheme => self%equat&
    &ion_base%weno_scheme, sensor_threshold => self%equation_base%sensor_threshold, w_aux_gpu => self%w_a&
    &ux_gpu, ep_ord_change_gpu => self%ep_ord_change_gpu, lmax_tag_gpu => self%lmax_tag_gpu,&
    & wall_tag_gpu => self%base_gpu%wall_tag_gpu, count_weno_x => self%equation_base%count_weno_x,&
    & count_weno_y => self%equation_base%count_weno_y, count_weno_z => self%equation_base%count_weno_z)

      if (grid_dim==1) then
        call count_weno_cuf(nv, nv_aux, nx, ny, nz, ng, eul_imin, eul_imax, eul_jmin, eul_jmax,&
        & eul_kmin, eul_kmax, weno_scheme, sensor_threshold, w_aux_gpu, ep_ord_change_gpu, count_weno_x,&
        & count_weno_y, count_weno_z)
      else
        lmax = ep_order/2
        call count_weno_c2_cuf(nv, nv_aux, nx, ny, nz, ng, eul_imin, eul_imax, eul_jmin, eul_jmax,&
        & eul_kmin, eul_kmax, lmax, weno_scheme, sensor_threshold, w_aux_gpu, ep_ord_change_gpu ,&
        & lmax_tag_gpu, wall_tag_gpu, count_weno_x, count_weno_y, count_weno_z)
      endif

      call mpi_allreduce(mpi_in_place, count_weno_x, 1, mpi_prec, mpi_sum, mpi_comm_world, self%mpi_err)
      call mpi_allreduce(mpi_in_place, count_weno_y, 1, mpi_prec, mpi_sum, mpi_comm_world, self%mpi_err)
      call mpi_allreduce(mpi_in_place, count_weno_z, 1, mpi_prec, mpi_sum, mpi_comm_world, self%mpi_err)
      ntot = real(eul_imax-eul_imin+1,rkind)*real(ny,rkind)*real(nz,rkind)*real(nprocs,rkind)
      count_weno_x_frac = count_weno_x/ntot
      ntot = real(nx,rkind)*real(eul_jmax-eul_jmin+1,rkind)*real(nz,rkind)*real(nprocs,rkind)
      count_weno_y_frac = count_weno_y/ntot
      ntot = real(nx,rkind)*real(ny,rkind)*real(eul_kmax-eul_kmin+1,rkind)*real(nprocs,rkind)
      count_weno_z_frac = count_weno_z/ntot
      open(unit=88,file='count_weno_xyz.dat',position='append')
      write(88,'(3e15.8)') count_weno_x_frac, count_weno_y_frac, count_weno_z_frac
      close(88)

    endassociate
  endsubroutine count_weno

  subroutine compute_airfoil_forces_runtime(self)
    class(equation_singleideal_gpu_object), intent(inout) :: self
    real(rkind) :: n, a, pn, pa, tn, ta, al, pdyn

    associate(nx => self%field%nx, ny => self%field%ny, nz => self%field%nz, ng => self%grid%ng,&
    & nv => self%field%nv, nzmax => self%grid%nzmax, dzitdz => self%field%dzitdz, icyc => self%equation_b&
    &ase%icyc, time => self%equation_base%time, aoa => self%equation_base%aoa, p0 => self%equation_base%p&
    &0, u0 => self%equation_base%u0, rgas0 => self%equation_base%rgas0, w_aux_gpu => self%w_aux_gpu,&
    & wall_tag_gpu => self%base_gpu%wall_tag_gpu, meta_gpu => self%base_gpu%meta_gpu, csimod_gpu => self%&
    &base_gpu%csimod_gpu, dxdcsic2_gpu => self%base_gpu%dxdcsic2_gpu, dydcsic2_gpu => self%base_gpu%dydcs&
    &ic2_gpu, lift_af => self%equation_base%lift_af, drag_af => self%equation_base%drag_af,&
    & pres_af => self%equation_base%pres_af, fric_af => self%equation_base%fric_af)

      call compute_airfoil_forces_runtime_cuf(nx, ny, nz, ng, nv, p0, u0, rgas0, w_aux_gpu,&
      & meta_gpu, csimod_gpu, wall_tag_gpu, dxdcsic2_gpu, dydcsic2_gpu, n, a, pn, pa, tn, ta)
      call mpi_allreduce(mpi_in_place, n, 1, mpi_prec, mpi_sum, self%equation_base%field%mp_cart, self%mpi_err)
      call mpi_allreduce(mpi_in_place, a, 1, mpi_prec, mpi_sum, self%equation_base%field%mp_cart, self%mpi_err)
      call mpi_allreduce(mpi_in_place, pn, 1, mpi_prec, mpi_sum, self%equation_base%field%mp_cart, self%mpi_err)
      call mpi_allreduce(mpi_in_place, pa, 1, mpi_prec, mpi_sum, self%equation_base%field%mp_cart, self%mpi_err)
      call mpi_allreduce(mpi_in_place, tn, 1, mpi_prec, mpi_sum, self%equation_base%field%mp_cart, self%mpi_err)
      call mpi_allreduce(mpi_in_place, ta, 1, mpi_prec, mpi_sum, self%equation_base%field%mp_cart, self%mpi_err)
      n = n/nzmax
      a = a/nzmax
      pn = pn/nzmax
      pa = pa/nzmax
      tn = tn/nzmax
      ta = ta/nzmax

      if(self%masterproc) then
        al = aoa
        pdyn = 0.5_rkind*u0*u0
        lift_af = ( n*cos(al)- a*sin(al))/pdyn
        drag_af = ( n*sin(al)+ a*cos(al))/pdyn
        pres_af = (pn*sin(al)+pa*cos(al))/pdyn
        fric_af = (tn*sin(al)+ta*cos(al))/pdyn
        open(unit=30,file='airfoil_forces_runtime.dat',position='append')
        write(30,100) icyc,time,lift_af,drag_af,pres_af,fric_af,n,a,pn,pa,tn,ta
        close(30)
      endif

    endassociate

    100 format(i0,2x,20es20.10)
  endsubroutine compute_airfoil_forces_runtime

  subroutine compute_dt(self)
    class(equation_singleideal_gpu_object), intent(inout) :: self
    real(rkind) :: dt_min, dtinv_max
    real(rkind) :: dtxi_max, dtyi_max, dtzi_max, dtxv_max, dtyv_max, dtzv_max, dtxk_max, dtyk_max, dtzk_max
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, dt => self%equation_base%d&
    &t, cfl => self%equation_base%cfl, prandtl => self%equation_base%prandtl, visc_model => self%visc_mod&
    &el, mu0 => self%mu0, t_ref_dim => self%t_ref_dim, powerlaw_vtexp => self%powerlaw_vtexp,&
    & sutherland_s => self%sutherland_s, w_gpu => self%base_gpu%w_gpu, w_aux_gpu => self%w_aux_gpu,&
    & dcsidx_gpu => self%base_gpu%dcsidx_gpu, detady_gpu => self%base_gpu%detady_gpu, dzitdz_gpu => self%&
    &base_gpu%dzitdz_gpu, dcsidxs_gpu => self%base_gpu%dcsidxs_gpu, detadys_gpu => self%base_gpu%detadys_&
    &gpu, dzitdzs_gpu => self%base_gpu%dzitdzs_gpu, dcsidxnc2_gpu => self%base_gpu%dcsidxnc2_gpu,&
    & dcsidync2_gpu => self%base_gpu%dcsidync2_gpu, detadxnc2_gpu => self%base_gpu%detadxnc2_gpu,&
    & detadync2_gpu => self%base_gpu%detadync2_gpu, mcsi_gpu => self%base_gpu%mcsi_gpu,&
    & meta_gpu => self%base_gpu%meta_gpu, calorically_perfect => self%equation_base%calorically_perfect,&
    & indx_cp_l => self%equation_base%indx_cp_l, indx_cp_r => self%equation_base%indx_cp_r,&
    & cp_coeff_gpu => self%cp_coeff_gpu, fluid_mask_gpu => self%fluid_mask_gpu, rgas0 => self%equation_ba&
    &se%rgas0, t0 => self%equation_base%t0, ndim => self%grid%ndim)
      if (cfl < 0) then
        dt = -cfl
      else
        if(self%grid%grid_dim == 1) then
          call compute_dt_cuf(nx, ny, nz, ng, rgas0, prandtl, dcsidx_gpu, detady_gpu, dzitdz_gpu,&
          & dcsidxs_gpu, detadys_gpu, dzitdzs_gpu, w_gpu, w_aux_gpu, dtxi_max, dtyi_max, dtzi_max, dtxv_max,&
          & dtyv_max, dtzv_max, dtxk_max, dtyk_max, dtzk_max, indx_cp_l, indx_cp_r, cp_coeff_gpu,&
          &fluid_mask_gpu,calorically_perfect,t0)
        elseif(self%grid%grid_dim == 2) then
          call compute_dt_c2_cuf(nx, ny, nz, ng, rgas0, prandtl, dcsidxnc2_gpu, dcsidync2_gpu,&
          & detadxnc2_gpu, detadync2_gpu, mcsi_gpu, meta_gpu, dzitdz_gpu, dzitdzs_gpu, w_gpu, w_aux_gpu,&
          & dtxi_max, dtyi_max, dtzi_max, dtxv_max, dtyv_max, dtzv_max, dtxk_max, dtyk_max, dtzk_max,&
          & indx_cp_l, indx_cp_r, cp_coeff_gpu,fluid_mask_gpu,calorically_perfect,t0)
        endif
        if(ndim == 2) then
          dtinv_max = maxval([dtxi_max, dtyi_max, dtzi_max, dtxv_max, dtyv_max, dtzv_max, dtxk_max, dtyk_max, dtzk_max])
        elseif(ndim == 3) then
          dtinv_max = maxval([dtxi_max, dtyi_max, dtzi_max, dtxv_max, dtyv_max, dtzv_max, dtxk_max, dtyk_max, dtzk_max])
        endif
        call mpi_allreduce(mpi_in_place,dtinv_max,1,mpi_prec,mpi_max,self%equation_base%field%mp_cart,self%mpi_err)
        dt_min = 1._rkind/dtinv_max
        dt = self%equation_base%cfl*dt_min
      endif
    endassociate
  endsubroutine compute_dt

  subroutine run(self, filename)

    class(equation_singleideal_gpu_object), intent(inout) :: self
    character(*) , intent(in) :: filename
    real(rkind) :: timing(1:2)
    real(rkind) :: timing_step(1:2)
    integer :: icyc_loop, iercuda, istat

    call self%initialize(filename=filename)

    associate(icyc0 => self%equation_base%icyc0, icyc => self%equation_base%icyc,&
    & time => self%equation_base%time, iter_dt_recompute => self%equation_base%iter_dt_recompute,&
    & residual_rhou => self%equation_base%residual_rhou, dpdx => self%equation_base%dpdx,&
    & rhobulk => self%equation_base%rhobulk, ubulk => self%equation_base%ubulk, tbulk => self%equation_ba&
    &se%tbulk, nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, visc_model => self%visc_model,&
    & mu0 => self%mu0, t_ref_dim => self%t_ref_dim, sutherland_s => self%sutherland_s,&
    & powerlaw_vtexp => self%powerlaw_vtexp, cv_coeff_gpu => self%cv_coeff_gpu, calorically_perfect => se&
    &lf%equation_base%calorically_perfect, indx_cp_l => self%equation_base%indx_cp_l, indx_cp_r => self%e&
    &quation_base%indx_cp_r, mode_async => self%equation_base%mode_async, flow_init => self%equation_base&
    &%flow_init, enable_forces_runtime => self%equation_base%enable_forces_runtime, count_weno_control =>&
    & self%equation_base%count_weno_control, time_from_last_rst => self%equation_base%time_from_last_rst,&
    & time_from_last_write => self%equation_base%time_from_last_write, time_from_last_stat => self%equati&
    &on_base%time_from_last_stat, time_from_last_slice => self%equation_base%time_from_last_slice,&
    & time_from_last_probe => self%equation_base%time_from_last_probe, time_from_last_slice_p3d => self%e&
    &quation_base%time_from_last_slice_p3d, time_from_last_insitu => self%equation_base%time_from_last_in&
    &situ, time_from_last_tspec => self%equation_base%time_from_last_tspec, time_is_freezed => self%equat&
    &ion_base%time_is_freezed)

      call self%update_ghost()
      call self%compute_aux()
      !@cuf iercuda=cudadevicesynchronize()

      if (mode_async >= 0) then
        if (self%equation_base%conservative_viscous==1) then
          call self%visflx(mode=2)
        else
          call self%visflx(mode=0)
        endif
        call self%base_gpu%bcswap_var(self%w_aux_gpu(:,:,:,8:8 ))
      endif
      call zero_flux_cuf(self%nx,self%ny,self%nz,self%nv,self%fl_gpu)
      if (self%equation_base%enable_plot3d>0) then
        call self%field%write_plot3d(mach=self%equation_base%mach, reynolds=self%equation_base%reyno&
        &lds, time=0._rkind, istore=0, plot3dgrid=.true., plot3dfield=.false., w_aux_io=self%equation_base%w_&
        &aux(1:self%nx,1:self%ny,1:self%nz,1:6))
      endif
      if(self%equation_base%enable_slice_plot3d > 0) then
        call self%equation_base%write_slice_p3d(plot3dgrid=.true.,plot3dfield=.false.)
      endif
      if (self%equation_base%restart_type==0) then
        self%equation_base%w_aux = self%w_aux_gpu
        if (self%equation_base%enable_plot3d>0) then
          call self%field%write_plot3d(mach=self%equation_base%mach, reynolds=self%equation_base%rey&
          &nolds, time=0._rkind, istore=0, plot3dgrid=.false., plot3dfield=.true., w_aux_io=self%equation_base%&
          &w_aux(1:self%nx,1:self%ny,1:self%nz,1:6))
        endif
        if (self%equation_base%enable_vtk>0) then
          call self%field%write_vtk(time=0._rkind, istore=0, w_aux_io=self%equation_base%w_aux(1:self%nx,1:self%ny,1:self%nz,1:6))
        endif
        if(self%equation_base%enable_slice_plot3d > 0) then
          call self%equation_base%write_slice_p3d(plot3dgrid=.false.,plot3dfield=.true.)
        endif
      endif
      call self%compute_dt()
      if (self%masterproc) write(*,*) 'dt =', self%equation_base%dt
      if (self%masterproc) write(*,*) 'dt*u0/l0 =', self%equation_base%dt*self%equation_base%u0/self%equation_base%l0
      if(self%equation_base%rand_type == 0) then
        call init_crandom_f(0,reproducible=.true.)
        if (self%masterproc) write(*,*) 'random numbers disabled'
      elseif(self%equation_base%rand_type < 0) then
        call init_crandom_f(self%myrank+1,reproducible=.false.)
        if (self%masterproc) write(*,*) 'random numbers not reproducible'
      else
        call init_crandom_f(self%myrank+1,reproducible=.true.)
        if (self%masterproc) write(*,*) 'random numbers reproducible'
      endif

      call mpi_barrier(mpi_comm_world, self%error) ; timing(1) = mpi_wtime()

      icyc_loop = icyc
      integration: do

        icyc_loop = icyc_loop + 1
        time_is_freezed = self%equation_base%time_is_freezed_fun()

        if ( time_is_freezed ) then
          self%equation_base%dt = 0._rkind
        else
          icyc = icyc + 1
          if(self%equation_base%cfl > 0 .and. mod(icyc-icyc0, iter_dt_recompute)==0) then
            call self%compute_aux(central=1, ghost=0)
            !@cuf iercuda=cudadevicesynchronize()
            call self%compute_dt()
          endif

          if(mode_async == 0) call self%rk_sync()
          if(mode_async == 1) call self%rk_async()
          if (mod(icyc-icyc0, self%equation_base%print_control)==0) then
            call self%compute_residual()
            if (ieee_is_nan(self%equation_base%residual_rhou)) then
              if (self%masterproc) write(*,*) 'boom!!!'
              call self%base_gpu%copy_gpu_cpu()
              if (self%equation_base%enable_plot3d>0) then
                call self%field%write_plot3d(mach=self%equation_base%mach, reynolds=self%equation_ba&
                &se%reynolds, time=0._rkind, istore=0, plot3dgrid=.true., plot3dfield=.true., w_aux_io=self%equation_&
                &base%w_aux(1:self%nx,1:self%ny,1:self%nz,1:6))
              endif
              if (self%equation_base%enable_vtk>0) then
                call self%field%write_vtk(time=0._rkind, istore=0, w_aux_io=self%equation_base%w_aux&
                &(1:self%nx,1:self%ny,1:self%nz,1:6))
              endif
              call mpi_barrier(mpi_comm_world,self%mpi_err)
              call mpi_abort(mpi_comm_world,99,self%mpi_err)
            endif
          endif
          call self%analyze_runtime()
        endif
        call self%manage_output()

        self%equation_base%time = self%equation_base%time + self%equation_base%dt

        if(self%masterproc.and.mod(icyc-icyc0, self%equation_base%print_control)==0) then
          call self%print_progress()
        endif
        if ((self%equation_base%icyc-self%equation_base%icyc0) >= self%num_iter) exit integration
        call mpi_barrier(mpi_comm_world, self%error) ; timing_step(2) = mpi_wtime()
      enddo integration

      if (allocated(self%equation_base%islice)) then
        wait(133)
        close(133)
      endif
      if (allocated(self%equation_base%jslice)) then
        wait(134)
        close(134)
      endif
      if (allocated(self%equation_base%kslice)) then
        wait(135)
        close(135)
      endif
      if (self%equation_base%num_probe>0) close(136)

      call mpi_barrier(mpi_comm_world, self%error) ; timing(2) = mpi_wtime()
      if(self%num_iter > 0) then
        if (self%masterproc) then
          write(*,'(a, f18.10)') 'averaged timing: ', (timing(2) - timing(1))/(self%equation_base%icyc-self%equation_base%icyc0)
        endif
      endif

      call self%base_gpu%copy_gpu_cpu()
      if (self%equation_base%io_type_w==1) then
        call self%field%write_field_serial()
        if (self%equation_base%enable_stat_time>0) call self%equation_base%write_stats_serial()
        if (self%equation_base%enable_stat_3d>0) call self%equation_base%write_stats_3d_serial()
        if (self%equation_base%enable_tspec>0) then
          self%equation_base%w_tspec = self%w_tspec_gpu
          call self%equation_base%write_tspec()
          call self%equation_base%write_psd_avg_tspec()
        endif
      endif
      if (self%equation_base%io_type_w==2) then
        call self%field%write_field()
        if (self%equation_base%enable_stat_time>0) call self%equation_base%write_stats()
        if (self%equation_base%enable_stat_3d>0) call self%equation_base%write_stats_3d()
        if (self%equation_base%enable_tspec>0) then
          self%equation_base%w_tspec = self%w_tspec_gpu
          call self%equation_base%write_tspec()
          call self%equation_base%write_psd_avg_tspec()
        endif
      endif

      call self%equation_base%write_field_info()

    endassociate
  endsubroutine run

  subroutine print_progress(self)
    class(equation_singleideal_gpu_object), intent(inout) :: self
    character(6) :: pos_io
    associate(icyc => self%equation_base%icyc, time => self%equation_base%time, dt => self%equation_&
    &base%dt, residual_rhou => self%equation_base%residual_rhou, dpdx => self%equation_base%dpdx,&
    & dpth => self%equation_base%dpth, vmax => self%equation_base%vmax, rhobulk => self%equation_base%rho&
    &bulk, ubulk => self%equation_base%ubulk, tbulk => self%equation_base%tbulk, rhomin => self%equation_&
    &base%rhomin, rhomax => self%equation_base%rhomax, pmin => self%equation_base%pmin,&
    & pmax => self%equation_base%pmax, tmin => self%equation_base%tmin, tmax => self%equation_base%tmax)
      residual_rhou = residual_rhou/(self%equation_base%u0**2)/self%equation_base%rho0*self%equation_base%l0
      dpdx = dpdx /(self%equation_base%u0**2)/self%equation_base%rho0*self%equation_base%l0
      pos_io = 'append'
      if (self%equation_base%icyc==1) pos_io = 'rewind'
      if (self%masterproc) then
        open(unit=15,file='progress.out',position=pos_io)
        if (self%equation_base%channel_case) then
          if (self%grid%grid_dim == 1) then
            write(* ,'(1i10,20es20.10)') icyc, dt, time, residual_rhou, -dpdx/dt, rhobulk, ubulk,&
            & tbulk, rhomin, rhomax, tmin, tmax, pmin, pmax
            write(15,'(1i10,20es20.10)') icyc, dt, time, residual_rhou, -dpdx/dt, rhobulk, ubulk,&
            & tbulk, rhomin, rhomax, tmin, tmax, pmin, pmax
          elseif (self%grid%grid_dim == 2) then
            write(* ,'(1i10,20es20.10)') icyc, dt, time, residual_rhou, -dpth/dt, rhobulk, ubulk,&
            & tbulk, rhomin, rhomax, tmin, tmax, pmin, pmax
            write(15,'(1i10,20es20.10)') icyc, dt, time, residual_rhou, -dpth/dt, rhobulk, ubulk,&
            & tbulk, rhomin, rhomax, tmin, tmax, pmin, pmax
          endif
        else if(self%equation_base%flow_init == 5) then
          write(* ,'(1i10,20es20.10)') icyc, dt, time, residual_rhou, rhomin,rhomax,tmin,tmax,pmin,pmax
          write(15,'(1i10,20es20.10)') icyc, dt, time, residual_rhou, rhomin,rhomax,tmin,tmax,pmin,pmax
        else
          write(* ,'(1i10,20es20.10)') icyc, dt, time, residual_rhou, rhomin, rhomax, tmin, tmax, pmin, pmax
          write(15,'(1i10,20es20.10)') icyc, dt, time, residual_rhou, rhomin, rhomax, tmin, tmax, pmin, pmax
        endif
        close(15)
      endif
    endassociate
  endsubroutine print_progress




  subroutine insitu_coprocess(self)
    class(equation_singleideal_gpu_object), intent(inout) :: self
    integer :: l,i,j,k,ll
    call self%update_ghost()
    call self%base_gpu%bcswap_edges_corners()
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, visc_model => self%visc_mo&
    &del, mu0 => self%mu0, t_ref_dim => self%t_ref_dim, sutherland_s => self%sutherland_s,&
    & powerlaw_vtexp => self%powerlaw_vtexp, w_aux_gpu => self%w_aux_gpu, cv_coeff_gpu => self%cv_coeff_g&
    &pu, indx_cp_l => self%equation_base%indx_cp_l, indx_cp_r => self%equation_base%indx_cp_r,&
    & calorically_perfect => self%equation_base%calorically_perfect)
      call self%compute_aux()
    endassociate
    if (self%equation_base%npsi > 0) then
      call self%insitu_compute_psi()
      do l=1,self%equation_base%npsi
        call self%base_gpu%bcswap_var(self%psi_gpu(:,:,:,l:l ))
      enddo
      do l=1,self%equation_base%npsi
        call self%base_gpu%bcswap_edges_corners_var(self%psi_gpu(:,:,:,l:l ))
      enddo
    endif
    associate( nxsl_ins => self%equation_base%nxsl_ins, nxel_ins => self%equation_base%nxel_ins,&
    & nysl_ins => self%equation_base%nysl_ins, nyel_ins => self%equation_base%nyel_ins,&
    & nzsl_ins => self%equation_base%nzsl_ins, nzel_ins => self%equation_base%nzel_ins,&
    & npsi => self%equation_base%npsi, npsi_pv => self%equation_base%npsi_pv, nv_aux => self%equation_bas&
    &e%nv_aux, n_aux_list => self%equation_base%n_aux_list, psi_gpu => self%psi_gpu, n_add_list => self%e&
    &quation_base%n_add_list, w_aux_gpu => self%w_aux_gpu, aux_list_gpu => self%aux_list_gpu,&
    & add_list_gpu => self%add_list_gpu, nx => self%field%nx, ny => self%field%ny, nz => self%field%nz,&
    & icyc => self%equation_base%icyc, time => self%equation_base%time, i_insitu => self%equation_base%i_&
    &insitu, flag => self%equation_base%insitu_flag, nrank => self%myrank, ng => self%grid%ng,&
    & aux_list_name => self%equation_base%aux_list_name, add_list_name => self%equation_base%add_list_nam&
    &e, time_insitu => self%equation_base%time_insitu)

      self%equation_base%w_aux = w_aux_gpu
      self%equation_base%psi = psi_gpu

      do k=nzsl_ins,nzel_ins
        do j=nysl_ins,nyel_ins
          do i=nxsl_ins,nxel_ins
            do l=1,n_aux_list
              ll = self%equation_base%aux_list(l)
              self%equation_base%psi_pv(i,j,k,l) = self%equation_base%w_aux(i,j,k,ll)
            enddo
            do l=1,npsi
              self%equation_base%psi_pv(i,j,k,n_aux_list+l) = self%equation_base%psi(i,j,k,l)
            enddo
          enddo
        enddo
      enddo

      if(self%equation_base%insitu_platform == "catalyst-v1") then
        call requestdatadescription(i_insitu,time_insitu,flag)
        if (flag.ne.0) then
          call needtocreategrid(flag)
          if (flag.ne.0) then
            do l=1,n_aux_list
              call addfieldtostructured("3d_struct"//c_null_char,self%equation_base%psi_pv(:,:,:,l),&
              & trim(adjustl(aux_list_name(l)))//c_null_char,nrank)
            enddo
            do l=1,n_add_list
              call addfieldtostructured("3d_struct"//c_null_char,self%equation_base%psi_pv(:,:,:,&
              &n_aux_list+l), trim(adjustl(add_list_name(l)))//c_null_char,nrank)
            enddo
            call coprocess()
          end if
        end if
      elseif(self%equation_base%insitu_platform == "catalyst-v2") then
        call self%insitu_do_catalyst_execute()
      endif
    endassociate

  end subroutine insitu_coprocess

  subroutine insitu_do_catalyst_execute(self)
    class(equation_singleideal_gpu_object), intent(inout) :: self

    integer :: cycle
    real(real64) :: time
    type(c_ptr) :: catalyst_exec_params, mesh, info
    type(c_ptr) :: xt, yt, zt
    type(c_ptr), dimension(:), allocatable :: vx
    integer(kind(catalyst_status)) :: code
    integer :: exit_code
    integer :: l

    catalyst_exec_params = catalyst_conduit_node_create()
    call catalyst_conduit_node_set_path_int64(catalyst_exec_params, "catalyst/state/timestep",&
    & int(self%equation_base%i_insitu,int64))
    call catalyst_conduit_node_set_path_float_32_64(catalyst_exec_params, "catalyst/state/time", self%equation_base%time_insitu)


    call catalyst_conduit_node_set_path_char8_str(catalyst_exec_params, "catalyst/channels/grid/type", "mesh")

    mesh = catalyst_conduit_node_create()

    associate(points_x => self%equation_base%points_x, points_y => self%equation_base%points_y,&
    & points_z => self%equation_base%points_z, n_points_x => self%equation_base%n_points_x,&
    & n_points_y => self%equation_base%n_points_y, n_points_z => self%equation_base%n_points_z,&
    & n_points => self%equation_base%n_points, n_aux_list => self%equation_base%n_aux_list,&
    & n_add_list => self%equation_base%n_add_list, aux_list_name => self%equation_base%aux_list_name,&
    & add_list_name => self%equation_base%add_list_name)

      allocate(vx(n_aux_list+n_add_list))

      call catalyst_conduit_node_set_path_char8_str(mesh,"coordsets/coords/type","explicit")
      xt = catalyst_conduit_node_create()
      yt = catalyst_conduit_node_create()
      if(n_points_z > 1) zt = catalyst_conduit_node_create()
      call catalyst_conduit_node_set_external_float_32_64_ptr(xt, points_x, n_points)
      call catalyst_conduit_node_set_external_float_32_64_ptr(yt, points_y, n_points)
      if(n_points_z > 1) call catalyst_conduit_node_set_external_float_32_64_ptr(zt, points_z, n_points)
      call catalyst_conduit_node_set_path_external_node(mesh, "coordsets/coords/values/x", xt)
      call catalyst_conduit_node_set_path_external_node(mesh, "coordsets/coords/values/y", yt)
      if(n_points_z > 1) call catalyst_conduit_node_set_path_external_node(mesh, "coordsets/coords/values/z", zt)
      call catalyst_conduit_node_set_path_char8_str(mesh, "topologies/mesh/type", "structured")
      call catalyst_conduit_node_set_path_char8_str(mesh, "topologies/mesh/coordset", "coords")
      call catalyst_conduit_node_set_path_int32(mesh, "topologies/mesh/elements/dims/i", n_points_x-1)
      call catalyst_conduit_node_set_path_int32(mesh, "topologies/mesh/elements/dims/j", n_points_y-1)
      if(n_points_z > 1) call catalyst_conduit_node_set_path_int32(mesh, "topologies/mesh/elements/dims/k", n_points_z-1)

      do l=1,n_aux_list
        vx(l) = catalyst_conduit_node_create()
        call catalyst_conduit_node_set_external_float_32_64_ptr(vx(l), self%equation_base%psi_pv(:,:,:,l), n_points)
        call catalyst_conduit_node_set_path_char8_str(mesh, "fields/"//trim(adjustl(aux_list_name(l)))//"/association", "vertex")
        call catalyst_conduit_node_set_path_char8_str(mesh, "fields/"//trim(adjustl(aux_list_name(l)))//"/topology", "mesh")
        call catalyst_conduit_node_set_path_char8_str(mesh, "fields/"//trim(adjustl(aux_list_name(l)&
        &))//"/volume_dependent", "false")
        call catalyst_conduit_node_set_path_external_node(mesh, "fields/"//trim(adjustl(aux_list_name(l)))//"/values", vx(l))
      enddo
      do l=1,n_add_list
        vx(n_aux_list+l) = catalyst_conduit_node_create()
        call catalyst_conduit_node_set_external_float_32_64_ptr(vx(n_aux_list+l),&
        & self%equation_base%psi_pv(:,:,:,n_aux_list+l), n_points)
        call catalyst_conduit_node_set_path_char8_str(mesh, "fields/"//trim(adjustl(add_list_name(l)))//"/association", "vertex")
        call catalyst_conduit_node_set_path_char8_str(mesh, "fields/"//trim(adjustl(add_list_name(l)))//"/topology", "mesh")
        call catalyst_conduit_node_set_path_char8_str(mesh, "fields/"//trim(adjustl(add_list_name(l)&
        &))//"/volume_dependent", "false")
        call catalyst_conduit_node_set_path_external_node(mesh, "fields/"//trim(adjustl(add_list_nam&
        &e(l)))//"/values", vx(n_aux_list+l))
      enddo


      call catalyst_conduit_node_set_path_external_node(catalyst_exec_params, "catalyst/channels/grid/data", mesh)

#if 1
      call catalyst_conduit_node_print(catalyst_exec_params)

#endif

      code = c_catalyst_execute(catalyst_exec_params)
      if (code /= catalyst_status_ok) then
        write (error_unit, *) "failed to call `execute`:", code
        exit_code = 1
      end if
      call catalyst_conduit_node_destroy(catalyst_exec_params)
      call catalyst_conduit_node_destroy(mesh)
      call catalyst_conduit_node_destroy(xt)
      call catalyst_conduit_node_destroy(yt)
      if(n_points_z > 1) call catalyst_conduit_node_destroy(zt)
      do l=1,n_aux_list+n_add_list
        call catalyst_conduit_node_destroy(vx(l))
      enddo

      deallocate(vx)

    endassociate

  endsubroutine insitu_do_catalyst_execute

  subroutine insitu_compute_psi(self)
    class(equation_singleideal_gpu_object), intent(inout) :: self
    integer :: l
    type(dim3) :: grid, tblock
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv,&
    & nv_aux => self%nv_aux, dt => self%equation_base%dt, visc_order => self%equation_base%visc_order,&
    & prandtl => self%equation_base%prandtl, visc_model => self%visc_model, mu0 => self%mu0,&
    & u0 => self%equation_base%u0, l0 => self%equation_base%l0, t_ref_dim => self%t_ref_dim,&
    & sutherland_s => self%sutherland_s, powerlaw_vtexp => self%powerlaw_vtexp, coeff_deriv1_gpu => self%&
    &coeff_deriv1_gpu, coeff_deriv2_gpu => self%coeff_deriv2_gpu, fhat_trans_gpu => self%fhat_trans_gpu,&
    & fl_trans_gpu => self%fl_trans_gpu, fl_gpu => self%fl_gpu, w_aux_gpu => self%w_aux_gpu,&
    & w_aux_trans_gpu => self%w_aux_trans_gpu, dcsidx_gpu => self%base_gpu%dcsidx_gpu,&
    & detady_gpu => self%base_gpu%detady_gpu, dzitdz_gpu => self%base_gpu%dzitdz_gpu, dcsidxs_gpu => self&
    &%base_gpu%dcsidxs_gpu, detadys_gpu => self%base_gpu%detadys_gpu, dzitdzs_gpu => self%base_gpu%dzitdz&
    &s_gpu, dcsidx2_gpu => self%base_gpu%dcsidx2_gpu, detady2_gpu => self%base_gpu%detady2_gpu,&
    & dzitdz2_gpu => self%base_gpu%dzitdz2_gpu, dcsidxc2_gpu => self%base_gpu%dcsidxc2_gpu,&
    & detadyc2_gpu => self%base_gpu%detadyc2_gpu, detadxc2_gpu => self%base_gpu%detadxc2_gpu,&
    & dcsidyc2_gpu => self%base_gpu%dcsidyc2_gpu, vis_tag_gpu => self%vis_tag_gpu, wall_tag_gpu => self%b&
    &ase_gpu%wall_tag_gpu, eul_imin => self%equation_base%eul_imin, eul_imax => self%equation_base%eul_im&
    &ax, eul_jmin => self%equation_base%eul_jmin, eul_jmax => self%equation_base%eul_jmax,&
    & eul_kmin => self%equation_base%eul_kmin, eul_kmax => self%equation_base%eul_kmax,&
    & cv_coeff_gpu => self%cv_coeff_gpu, indx_cp_l => self%equation_base%indx_cp_l, indx_cp_r => self%equ&
    &ation_base%indx_cp_r, calorically_perfect => self%equation_base%calorically_perfect,&
    & npsi => self%equation_base%npsi, psi_gpu => self%psi_gpu, add_list => self%equation_base%add_list,&
    & x_gpu => self%base_gpu%x_gpu )
      do l=1,npsi
        if (add_list(l) == 1) then
          call insitu_div_cuf(nx, ny, nz, ng, visc_order, npsi,l, w_aux_gpu, coeff_deriv1_gpu,&
          & dcsidx_gpu, detady_gpu, dzitdz_gpu, psi_gpu )
        endif
        if (add_list(l) == 2) then
          call insitu_omega_cuf(nx, ny, nz, ng, visc_order, npsi,l, w_aux_gpu, coeff_deriv1_gpu,&
          & dcsidx_gpu, detady_gpu, dzitdz_gpu, psi_gpu )
        endif
        if (add_list(l) == 3) then
          call insitu_ducros_cuf(nx, ny, nz, ng, visc_order, npsi,l, u0,l0,w_aux_gpu,&
          & coeff_deriv1_gpu, dcsidx_gpu, detady_gpu, dzitdz_gpu, psi_gpu )
        endif
        if (add_list(l) == 4) then
          tblock = dim3(eulerweno_threads_x,eulerweno_threads_y,1)
          grid = dim3(ceiling(real(nx)/tblock%x),ceiling(real(nz)/tblock%y),1)
          if(self%grid%grid_dim == 1) then
            call insitu_swirling_kernel<<<grid, tblock>>>(nv, nx, ny, nz, visc_order, ng, npsi,l,&
            &dcsidx_gpu, detady_gpu, dzitdz_gpu, w_aux_gpu, coeff_deriv1_gpu, psi_gpu, u0, x_gpu )
          elseif(self%grid%grid_dim == 2) then
            call insitu_swirling_c2_kernel<<<grid, tblock>>>(nv, nx, ny, nz, visc_order, ng, npsi,l,&
            &dzitdz_gpu, w_aux_gpu, coeff_deriv1_gpu, dcsidxc2_gpu, detadyc2_gpu, detadxc2_gpu, dcsidyc2_gpu,&
            & psi_gpu, u0, x_gpu, vis_tag_gpu, wall_tag_gpu)
          endif
        endif
        if (add_list(l) == 5) then
          tblock = dim3(eulerweno_threads_x,eulerweno_threads_y,1)
          grid = dim3(ceiling(real(nx)/tblock%x),ceiling(real(nz)/tblock%y),1)
          if(self%grid%grid_dim == 1) then
            call insitu_schlieren_kernel<<<grid, tblock>>>(nv, nx, ny, nz, visc_order, ng, npsi,l,&
            &dcsidx_gpu, detady_gpu, dzitdz_gpu, w_aux_gpu, coeff_deriv1_gpu, psi_gpu, u0, x_gpu )
          elseif(self%grid%grid_dim == 2) then
            call insitu_schlieren_c2_kernel<<<grid, tblock>>>(nv, nx, ny, nz, visc_order, ng, npsi,&
            &l,dzitdz_gpu, w_aux_gpu, coeff_deriv1_gpu, dcsidxc2_gpu, detadyc2_gpu, detadxc2_gpu, dcsidyc2_gpu,&
            & psi_gpu, u0, x_gpu, vis_tag_gpu, wall_tag_gpu)
          endif
        endif
      enddo
    endassociate

  end subroutine insitu_compute_psi

  subroutine manage_output(self)
    class(equation_singleideal_gpu_object), intent(inout) :: self
    integer :: i,j,k,ii,jj,kk,l,n,m,mm
    integer :: isize, jsize, ksize
    character(3) :: chx, chy, chz
    logical :: sliceyz_exist, slicexz_exist, slicexy_exist, probe_exist
    integer :: iercuda, istat
    associate(time_from_last_rst => self%equation_base%time_from_last_rst, time_from_last_write => s&
    &elf%equation_base%time_from_last_write, time_from_last_stat => self%equation_base%time_from_last_sta&
    &t, time_from_last_slice => self%equation_base%time_from_last_slice, time_from_last_probe => self%equ&
    &ation_base%time_from_last_probe, time_from_last_slice_p3d => self%equation_base%time_from_last_slice&
    &_p3d, time_from_last_insitu => self%equation_base%time_from_last_insitu, time_from_last_tspec => sel&
    &f%equation_base%time_from_last_tspec, icyc0 => self%equation_base%icyc0, icyc => self%equation_base%&
    &icyc, time => self%equation_base%time, w_aux_gpu => self%w_aux_gpu, grid_dim => self%grid%grid_dim,&
    & ijk_probe_gpu => self%ijk_probe_gpu, w_aux_probe_gpu => self%w_aux_probe_gpu, probe_coeff_gpu => se&
    &lf%probe_coeff_gpu, w_aux_probe => self%equation_base%w_aux_probe, wallprop => self%equation_base%wa&
    &llprop, wallprop_gpu => self%wallprop_gpu, nx => self%nx, ny => self%ny, nz => self%nz,&
    & ng => self%ng, nv_aux => self%nv_aux, nzmax => self%grid%nzmax, mp_cartz => self%field%mp_cartz,&
    & mpi_err => self%mpi_err, time_insitu => self%equation_base%time_insitu, ndft => self%equation_base%&
    &ndft_tspec, jg_tspec => self%equation_base%jg_tspec, dt_tspec => self%equation_base%dt_tspec,&
    & it_win1_tspec => self%equation_base%it_win1_tspec, it_win2_tspec => self%equation_base%it_win2_tspe&
    &c, it_nwin_tspec => self%equation_base%it_nwin_tspec, w_tspec_gpu => self%w_tspec_gpu,&
    & w_psd_tspec_gpu => self%w_psd_tspec_gpu, w_psd_tspec => self%equation_base%w_psd_tspec,&
    & w_psd_avg_tspec => self%equation_base%w_psd_avg_tspec, wss => self%equation_base%wss_tspec,&
    & wall_tag_gpu => self%base_gpu%wall_tag_gpu, meta_gpu => self%base_gpu%meta_gpu, csimod_gpu => self%&
    &base_gpu%csimod_gpu, dxdcsic2_gpu => self%base_gpu%dxdcsic2_gpu, dydcsic2_gpu => self%base_gpu%dydcs&
    &ic2_gpu )

      time_from_last_write = time_from_last_write + self%equation_base%dt
      if (time_from_last_write >= self%equation_base%dtsave) then
        if (self%masterproc) write(*,*) 'time_from_last_write=',time_from_last_write
        if (self%masterproc) write(*,*) 'istore =', self%equation_base%istore
        call self%base_gpu%copy_gpu_cpu()
        self%equation_base%w_aux = self%w_aux_gpu
        if (self%equation_base%enable_plot3d>0) then
          call self%field%write_plot3d(mach=self%equation_base%mach, reynolds=self%equation_base%rey&
          &nolds, time=time, istore=self%equation_base%istore, plot3dgrid=.false., plot3dfield=.true.,&
          & w_aux_io=self%equation_base%w_aux(1:self%nx,1:self%ny,1:self%nz,1:6))
        endif
        if (self%equation_base%enable_vtk>0) then
          call self%field%write_vtk(time=time, istore=self%equation_base%istore, w_aux_io=self%equat&
          &ion_base%w_aux(1:self%nx,1:self%ny,1:self%nz,1:6))
        endif
        time_from_last_write = 0._rkind
        self%equation_base%istore = self%equation_base%istore + 1
      endif

      time_from_last_stat = time_from_last_stat + self%equation_base%dt
      if (time_from_last_stat >= self%equation_base%dtstat) then
        if (self%masterproc) write(*,*) 'time_from_last_stat=',time_from_last_stat
        if (self%masterproc) write(*,*) 'itav =',self%equation_base%itav
        call self%base_gpu%copy_gpu_cpu()


        self%equation_base%w_aux(:,:,:,6) = self%w_aux_gpu(:,:,:,6)
        self%equation_base%w_aux(:,:,:,7) = self%w_aux_gpu(:,:,:,7)

        if (self%grid%grid_dim == 1) then
          call self%equation_base%compute_stats()
        elseif(self%grid%grid_dim == 2) then
          call self%equation_base%compute_stats_c2()
          if(self%equation_base%flow_init == 5) then
            call self%equation_base%compute_airfoil_forces()
          endif
        endif
        if (self%equation_base%save_stat_z) call self%equation_base%write_stats_z()
        if (self%equation_base%enable_stat_3d>0) call self%equation_base%compute_stats_3d()
        time_from_last_stat = 0._rkind
        self%equation_base%itav = self%equation_base%itav + 1
      endif

      if (self%equation_base%enable_insitu > 0) then
        if (self%equation_base%time_is_freezed) then
          time_from_last_insitu = time_from_last_insitu + self%equation_base%dt_insitu
          time_insitu = time_insitu+self%equation_base%dt_insitu
        else
          time_from_last_insitu = time_from_last_insitu + self%equation_base%dt
          time_insitu = time_insitu+self%equation_base%dt
        endif
        if (time_from_last_insitu >= self%equation_base%dt_insitu) then
          if (self%masterproc) write(*,*) 'time_from_last_insitu=',time_from_last_insitu
          if (self%masterproc) write(*,*) 'i_insitu =', self%equation_base%i_insitu
          call self%insitu_coprocess()
          time_from_last_insitu = time_from_last_insitu - self%equation_base%dt_insitu
          self%equation_base%i_insitu = self%equation_base%i_insitu + 1
        endif
      endif

      time_from_last_slice = time_from_last_slice + self%equation_base%dt
      time_from_last_probe = time_from_last_probe + self%equation_base%dt
      write(chx,'(i3.3)') self%field%ncoords(1)
      write(chy,'(i3.3)') self%field%ncoords(2)
      write(chz,'(i3.3)') self%field%ncoords(3)
      if (icyc-icyc0 == 1) then
        if (allocated(self%equation_base%islice)) then
          inquire(file='sliceyz_'//chx//'_'//chy//'_'//chz//'.bin', exist=sliceyz_exist)
          if(.not.(sliceyz_exist)) then
            open(133,file='sliceyz_'//chx//'_'//chy//'_'//chz//'.bin',form='unformatted', asynchronous="yes")
          else
            open(133,file='sliceyz_'//chx//'_'//chy//'_'//chz//'.bin',form='unformatted', position="append",asynchronous="yes")
          endif
        endif
        if (allocated(self%equation_base%jslice)) then
          inquire(file='slicexz_'//chx//'_'//chy//'_'//chz//'.bin', exist=slicexz_exist)
          if (.not.(slicexz_exist)) then
            open(134,file='slicexz_'//chx//'_'//chy//'_'//chz//'.bin',form='unformatted', asynchronous="yes")
          else
            open(134,file='slicexz_'//chx//'_'//chy//'_'//chz//'.bin',form='unformatted', position="append",asynchronous="yes")
          endif
        endif
        if (allocated(self%equation_base%kslice)) then
          inquire(file='slicexy_'//chx//'_'//chy//'_'//chz//'.bin', exist=slicexy_exist)
          if (.not.(slicexy_exist)) then
            open(135,file='slicexy_'//chx//'_'//chy//'_'//chz//'.bin',form='unformatted', asynchronous="yes")
          else
            open(135,file='slicexy_'//chx//'_'//chy//'_'//chz//'.bin',form='unformatted', position="append",asynchronous="yes")
          endif
        endif
        if (self%equation_base%num_probe>0) then
          inquire(file='probe_'//chx//'_'//chy//'_'//chz//'.bin', exist=probe_exist)
          if (.not.(probe_exist)) then
            open(136,file='probe_'//chx//'_'//chy//'_'//chz//'.bin',form='unformatted')
            write(136) self%equation_base%num_probe
            do i=1,self%equation_base%num_probe
              write(136) self%equation_base%id_probe(i),self%equation_base%probe_coord(1:3,i)
            enddo
          else
            open(136,file='probe_'//chx//'_'//chy//'_'//chz//'.bin',form='unformatted', position="append")
          endif
        endif
      endif
      if (time_from_last_slice >= self%equation_base%dtslice) then
        if (self%masterproc) write(*,*) 'time_from_last_slice=',time_from_last_slice
        if (allocated(self%equation_base%islice)) then
          isize = size(self%equation_base%islice)
          do i=1,size(self%equation_base%islice)
            ii = self%equation_base%islice(i)
            do m=1,self%equation_base%num_aux_slice
              mm = self%equation_base%list_aux_slice(m)
              sliceyz_aux(i,:,:,m) = self%w_aux_gpu(ii,:,:,mm)
            enddo
          enddo
          write(133,asynchronous="no") icyc,time
          write(133,asynchronous="no") self%grid%nxmax,self%grid%nymax,self%grid%nzmax
          write(133,asynchronous="no") isize,(self%equation_base%islice(i),i=1,isize)
          write(133,asynchronous="no") isize, self%field%ny, self%field%nz, self%equation_base%num_aux_slice
          write(133,asynchronous="yes") sliceyz_aux !(1:isize,1:self%field%ny,1:self%field%nz,1:6)
        endif
        if (allocated(self%equation_base%jslice)) then
          jsize = size(self%equation_base%jslice)
          do j=1,size(self%equation_base%jslice)
            jj = self%equation_base%jslice(j)
            do m=1,self%equation_base%num_aux_slice
              mm = self%equation_base%list_aux_slice(m)
              if(jj == 1 .and. mm >= 2 .and. mm <= 4) then
                if (grid_dim == 2 .and. mm == 2) then
                  call compute_wallprop_c2_cuf(nx,ny,nz,ng,w_aux_gpu,wallprop_gpu,dxdcsic2_gpu,&
                  &dydcsic2_gpu,csimod_gpu,meta_gpu,wall_tag_gpu)
                endif
                slicexz_aux(:,j,:,m) = wallprop_gpu(:,:,mm)
              else
                slicexz_aux(:,j,:,m) = self%w_aux_gpu(:,jj,:,mm)
              endif
            enddo
          enddo
          write(134,asynchronous="no") icyc,time
          write(134,asynchronous="no") self%grid%nxmax,self%grid%nymax,self%grid%nzmax
          write(134,asynchronous="no") jsize,(self%equation_base%jslice(j),j=1,jsize)
          write(134,asynchronous="no") self%field%nx, jsize, self%field%nz, self%equation_base%num_aux_slice
          write(134,asynchronous="yes") slicexz_aux !(1:self%field%nx,1:jsize,1:self%field%nz,1:6)
        endif
        if (allocated(self%equation_base%kslice)) then
          ksize = size(self%equation_base%kslice)
          do k=1,size(self%equation_base%kslice)
            kk = self%equation_base%kslice(k)
            do m=1,self%equation_base%num_aux_slice
              mm = self%equation_base%list_aux_slice(m)
              slicexy_aux(:,:,k,m) = self%w_aux_gpu(:,:,kk,mm)
            enddo
          enddo
          write(135,asynchronous="no") icyc,time
          write(135,asynchronous="no") self%grid%nxmax,self%grid%nymax,self%grid%nzmax
          write(135,asynchronous="no") ksize,(self%equation_base%kslice(k),k=1,ksize)
          write(135,asynchronous="no") self%field%nx, self%field%ny, ksize, self%equation_base%num_aux_slice
          write(135,asynchronous="yes") slicexy_aux !(1:self%field%nx,1:self%field%ny,1:ksize,1:6)
        endif
        time_from_last_slice = 0._rkind
      endif

      if (time_from_last_probe >= self%equation_base%dtprobe) then
        if (self%masterproc) write(*,*) 'time_from_last_probe=',time_from_last_probe
        if (self%equation_base%num_probe>0) then
          call probe_interpolation_cuf(self%equation_base%num_probe,nx,ny,nz,ng,nv_aux,&
          &ijk_probe_gpu,w_aux_probe_gpu,w_aux_gpu,probe_coeff_gpu)
          w_aux_probe = w_aux_probe_gpu
          write(136) time, ((w_aux_probe(l,n),l=1,6),n=1,self%equation_base%num_probe)
        endif
        time_from_last_probe = 0._rkind
      endif

      time_from_last_slice_p3d = time_from_last_slice_p3d + self%equation_base%dt
      if (time_from_last_slice_p3d >= self%equation_base%dtslice_p3d) then
        self%equation_base%itslice_p3d = self%equation_base%itslice_p3d + 1
        if(self%equation_base%enable_slice_plot3d > 0) then
          if (self%masterproc) write(*,*) 'time_from_last_slice_p3d=',time_from_last_slice_p3d

          if (allocated(self%equation_base%islice_p3d)) then
            do i=1,size(self%equation_base%islice_p3d)
              ii = self%equation_base%islice_p3d(i)
              do m=1,self%equation_base%num_aux_slice
                mm = self%equation_base%list_aux_slice(m)
                self%equation_base%w_aux(ii,:,:,mm) = self%w_aux_gpu(ii,:,:,mm)
              enddo
            enddo
          endif
          if (allocated(self%equation_base%jslice_p3d)) then
            do j=1,size(self%equation_base%jslice_p3d)
              jj = self%equation_base%jslice_p3d(j)
              do m=1,self%equation_base%num_aux_slice
                mm = self%equation_base%list_aux_slice(m)
                self%equation_base%w_aux(:,jj,:,mm) = self%w_aux_gpu(:,jj,:,mm)
              enddo
              if(jj == 1) then
                if (grid_dim == 2) then
                  call compute_wallprop_c2_cuf(nx,ny,nz,ng,w_aux_gpu,wallprop_gpu,dxdcsic2_gpu,&
                  &dydcsic2_gpu,csimod_gpu,meta_gpu,wall_tag_gpu)
                endif
                wallprop(:,:,:) = wallprop_gpu(:,:,:)
              endif
            enddo
          endif
          if (allocated(self%equation_base%kslice_p3d)) then
            do k=1,size(self%equation_base%kslice_p3d)
              kk = self%equation_base%kslice_p3d(k)
              do m=1,self%equation_base%num_aux_slice
                mm = self%equation_base%list_aux_slice(m)
                self%equation_base%w_aux(:,:,kk,mm) = self%w_aux_gpu(:,:,kk,mm)
              enddo
            enddo
          endif
          call self%equation_base%write_slice_p3d(plot3dgrid=.false.,plot3dfield=.true.)
        endif
        time_from_last_slice_p3d = time_from_last_slice_p3d - self%equation_base%dtslice_p3d
      endif

      time_from_last_tspec = time_from_last_tspec + self%equation_base%dt
      if (time_from_last_tspec >= self%equation_base%dt_tspec) then
        if(self%equation_base%enable_tspec > 0) then
          if (self%masterproc) write(*,*) 'time_from_last_tspec=',time_from_last_tspec
          call compute_tspec_cuf(nx, ny, nz, ng, ndft, jg_tspec, 1, it_win1_tspec, w_aux_gpu, w_tspec_gpu)
          if(self%equation_base%it_win2_tspec > 0) then
            call compute_tspec_cuf(nx, ny, nz, ng, ndft, jg_tspec, 2, it_win2_tspec, w_aux_gpu, w_tspec_gpu)
          endif
          if(it_win1_tspec == ndft .or. it_win2_tspec == ndft) then
            if(it_win1_tspec == ndft) then
              call compute_psd_tspec_cuf(nx, ny, nz, ndft, 1, dt_tspec, w_tspec_gpu, w_psd_tspec_gpu)
            endif
            if(it_win2_tspec == ndft) then
              call compute_psd_tspec_cuf(nx, ny, nz, ndft, 2, dt_tspec, w_tspec_gpu, w_psd_tspec_gpu)
            endif

            w_psd_tspec = w_psd_tspec_gpu
            call mpi_allreduce(mpi_in_place,w_psd_tspec,nx*ndft,mpi_prec,mpi_sum,mp_cartz,mpi_err)
            w_psd_tspec = dt_tspec*w_psd_tspec/nzmax/wss

            call self%equation_base%write_psd_tspec()

            w_psd_avg_tspec = w_psd_avg_tspec*it_nwin_tspec + w_psd_tspec
            w_psd_avg_tspec = w_psd_avg_tspec/(it_nwin_tspec+1)

            if(it_win1_tspec == ndft) it_win1_tspec = 0
            if(it_win2_tspec == ndft) it_win2_tspec = 0
            it_nwin_tspec = it_nwin_tspec + 1
          endif
          it_win1_tspec = it_win1_tspec + 1
          it_win2_tspec = it_win2_tspec + 1
        endif
        time_from_last_tspec = time_from_last_tspec - self%equation_base%dt_tspec
      endif

      time_from_last_rst = time_from_last_rst + self%equation_base%dt
      if (time_from_last_rst >= self%equation_base%dtsave_restart) then
        if (self%masterproc) write(*,*) 'time_from_last_rst=',time_from_last_rst
        call self%base_gpu%copy_gpu_cpu()

        if (self%equation_base%io_type_w==1) then
          call self%field%write_field_serial()
          if (self%equation_base%enable_stat_time>0) call self%equation_base%write_stats_serial()
          if (self%equation_base%enable_stat_3d>0) call self%equation_base%write_stats_3d_serial()
          if (self%equation_base%enable_tspec>0) then
            self%equation_base%w_tspec = self%w_tspec_gpu
            call self%equation_base%write_tspec()
            call self%equation_base%write_psd_avg_tspec()
          endif
        endif
        if (self%equation_base%io_type_w==2) then
          call self%field%write_field()
          if (self%equation_base%enable_stat_time>0) call self%equation_base%write_stats()
          if (self%equation_base%enable_stat_3d>0) call self%equation_base%write_stats_3d()
          if (self%equation_base%enable_tspec>0) then
            self%equation_base%w_tspec = self%w_tspec_gpu
            call self%equation_base%write_tspec()
            call self%equation_base%write_psd_avg_tspec()
          endif
        endif
        call self%equation_base%write_field_info()
        time_from_last_rst = 0._rkind
      endif

    endassociate
  end subroutine manage_output

endmodule streams_equation_singleideal_gpu_object

