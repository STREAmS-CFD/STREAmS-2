module streams_equation_singleideal_cpu_object

  use streams_base_cpu_object
  use streams_field_object
  use streams_grid_object
  use streams_kernels_cpu
  use streams_parameters
  use crandom_f_mod
  use streams_equation_singleideal_object
  use mpi

  use iso_c_binding
  use, intrinsic :: iso_fortran_env


  implicit none
  private
  public :: equation_singleideal_cpu_object

  integer :: stream1

  integer, parameter :: eulercentral_threads_x=128,eulercentral_threads_y=3
  integer, parameter :: eulerweno_threads_x=128 ,eulerweno_threads_y=2

  type :: equation_singleideal_cpu_object
    type(base_cpu_object) :: base_cpu
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
    real(rkind), allocatable, dimension(:,:) :: coeff_deriv1_cpu
    real(rkind), allocatable, dimension(:,:) :: coeff_deriv2_cpu

    integer :: ierr
    integer :: icyc0, num_iter
    integer :: visc_model
    real(rkind) :: mu0, sutherland_s , t_ref_dim, powerlaw_vtexp
    logical :: masterproc
    integer :: mpi_err
    real(rkind), allocatable, dimension(:) :: winf_cpu, winf_past_shock_cpu
    real(rkind), allocatable, dimension(:,:,:,:) :: w_aux_cpu
    real(rkind), allocatable, dimension(:,:,:,:) :: fhat_cpu
    real(rkind), allocatable, dimension(:,:,:,:) :: w_aux_trans_cpu, fl_trans_cpu, fhat_trans_cpu
    real(rkind), allocatable, dimension(:,:,:,:) :: fl_cpu, fln_cpu
    real(rkind), allocatable, dimension(:,:) :: dcoe_cpu
    real(rkind), allocatable, dimension(:,:,:,:) :: w_var, w_var_t
    integer, allocatable, dimension(:,:,:) :: fluid_mask_cpu
    integer, allocatable, dimension(:,:,:) :: fluid_mask_ini_cpu
    integer, allocatable, dimension(:,:,:,:) :: ep_ord_change_cpu

    real(rkind), allocatable, dimension(:,:,:,:) :: wrecyc_cpu
    real(rkind), allocatable, dimension(:,:,:) :: wrecycav_cpu
    real(rkind), dimension(:,:,:), allocatable :: wmean_cpu

    real(rkind), dimension(:,:,:), allocatable :: inflow_random_plane_cpu
    real(rkind), dimension(:), allocatable :: weta_inflow_cpu
    real(rkind), dimension(:), allocatable :: yplus_inflow_cpu, eta_inflow_cpu
    real(rkind), dimension(:), allocatable :: yplus_recyc_cpu, eta_recyc_cpu, eta_recyc_blend_cpu
    integer, dimension(:), allocatable :: map_j_inn_cpu, map_j_out_cpu, map_j_out_blend_cpu

    real(rkind), dimension(:), allocatable :: cv_coeff_cpu, cp_coeff_cpu

    real(rkind), dimension(:,:), allocatable :: w_aux_probe_cpu
    real(rkind), dimension(:,:,:,:), allocatable :: probe_coeff_cpu
    integer, dimension(:,:), allocatable :: ijk_probe_cpu
    real(rkind), allocatable, dimension(:,:,:) :: wallprop_cpu



  contains
    procedure, pass(self) :: compute_dt
    procedure, pass(self) :: initialize
    procedure, pass(self) :: compute_residual
    procedure, pass(self) :: print_progress
    procedure, pass(self) :: run
    procedure, pass(self) :: rk_sync_old
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

  endtype equation_singleideal_cpu_object

contains


  subroutine point_to_field(self, field)
    class(equation_singleideal_cpu_object), intent(inout) :: self
    class(field_object), target :: field
    self%field => field
  endsubroutine point_to_field

  subroutine point_to_grid(self, grid)
    class(equation_singleideal_cpu_object), intent(inout) :: self
    class(grid_object), target :: grid
    self%grid => grid
  endsubroutine point_to_grid





  subroutine rk_sync_old(self)
    class(equation_singleideal_cpu_object), intent(inout) :: self
    integer :: istep, iercuda
    real(rkind) :: rhodt, gamdt, alpdt
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv,&
    & nv_aux => self%nv_aux, dt => self%equation_base%dt, eul_imin => self%equation_base%eul_imin,&
    & eul_imax => self%equation_base%eul_imax, eul_jmin => self%equation_base%eul_jmin,&
    & eul_jmax => self%equation_base%eul_jmax, eul_kmin => self%equation_base%eul_kmin,&
    & eul_kmax => self%equation_base%eul_kmax, channel_case => self%equation_base%channel_case)

      if (channel_case) self%equation_base%dpdx = 0._rkind
      do istep=1,self%equation_base%nrk
        rhodt = self%equation_base%rhork(istep)*dt
        gamdt = self%equation_base%gamrk(istep)*dt
        alpdt = self%equation_base%alprk(istep)*dt

        call init_flux_subroutine(nx, ny, nz, nv, self%fl_cpu, self%fln_cpu, rhodt)
        if (self%equation_base%conservative_viscous==1) then
          call self%visflx(mode=2)
        else
          call self%visflx(mode=0)
        endif
        call bcextr_var_subroutine(nx, ny, nz, ng, self%w_aux_cpu(:,:,:,10:10))
        call self%base_cpu%bcswap_var(self%w_aux_cpu(:,:,:,10:10))
        call self%base_cpu%bcswap_var(self%w_aux_cpu(:,:,:,8:8 ))
        if (self%equation_base%conservative_viscous==1) then
          call self%visflx(mode=7)
        else
          call self%visflx(mode=1)
        endif

        call self%euler_x(eul_imin, eul_imax)

        call self%euler_y(eul_jmin,eul_jmax)

        call self%euler_z(eul_kmin,eul_kmax)
        if (self%equation_base%conservative_viscous==1) then
          call self%visflx(mode=5)
        endif
        call self%bc_nr()
        call update_flux_subroutine(nx, ny, nz, nv, self%fl_cpu, self%fln_cpu, gamdt)
        if (channel_case) call self%force_rhs()
        call update_field_subroutine(nx, ny, nz, ng, nv, self%base_cpu%w_cpu, self%fln_cpu, self%fluid_mask_cpu)
        if (channel_case) call self%force_var()
        call self%update_ghost()
        call self%compute_aux()

      enddo
    endassociate
  endsubroutine rk_sync_old

  subroutine visflx(self, mode)
    class(equation_singleideal_cpu_object), intent(inout) :: self
    integer, intent(in) :: mode
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv,&
    & nv_aux => self%nv_aux, dt => self%equation_base%dt, visc_order => self%equation_base%visc_order,&
    & prandtl => self%equation_base%prandtl, visc_model => self%visc_model, mu0 => self%mu0,&
    & u0 => self%equation_base%u0, l0 => self%equation_base%l0, t_ref_dim => self%t_ref_dim,&
    & sutherland_s => self%sutherland_s, powerlaw_vtexp => self%powerlaw_vtexp, coeff_deriv1_cpu => self%&
    &coeff_deriv1_cpu, coeff_deriv2_cpu => self%coeff_deriv2_cpu, fhat_trans_cpu => self%fhat_trans_cpu,&
    & fl_trans_cpu => self%fl_trans_cpu, fl_cpu => self%fl_cpu, w_aux_cpu => self%w_aux_cpu,&
    & w_aux_trans_cpu => self%w_aux_trans_cpu, dcsidx_cpu => self%base_cpu%dcsidx_cpu,&
    & detady_cpu => self%base_cpu%detady_cpu, dzitdz_cpu => self%base_cpu%dzitdz_cpu, dcsidxs_cpu => self&
    &%base_cpu%dcsidxs_cpu, detadys_cpu => self%base_cpu%detadys_cpu, dzitdzs_cpu => self%base_cpu%dzitdz&
    &s_cpu, dcsidx2_cpu => self%base_cpu%dcsidx2_cpu, detady2_cpu => self%base_cpu%detady2_cpu,&
    & dzitdz2_cpu => self%base_cpu%dzitdz2_cpu, x_cpu => self%base_cpu%x_cpu, y_cpu => self%base_cpu%y_cp&
    &u, z_cpu => self%base_cpu%z_cpu, eul_imin => self%equation_base%eul_imin, eul_imax => self%equation_&
    &base%eul_imax, eul_jmin => self%equation_base%eul_jmin, eul_jmax => self%equation_base%eul_jmax,&
    & eul_kmin => self%equation_base%eul_kmin, eul_kmax => self%equation_base%eul_kmax,&
    & cv_coeff_cpu => self%cv_coeff_cpu, cp_coeff_cpu => self%cp_coeff_cpu, indx_cp_l => self%equation_ba&
    &se%indx_cp_l, indx_cp_r => self%equation_base%indx_cp_r, calorically_perfect => self%equation_base%c&
    &alorically_perfect, t0 => self%equation_base%t0, channel_case => self%equation_base%channel_case)

      if (mode == 0) then
        call visflx_subroutine(nx, ny, nz, ng, visc_order, prandtl, t0, indx_cp_l, indx_cp_r,&
        & cp_coeff_cpu, calorically_perfect, u0, l0, self%base_cpu%w_cpu, self%w_aux_cpu, self%fl_cpu,&
        & coeff_deriv1_cpu, coeff_deriv2_cpu, dcsidx_cpu, detady_cpu, dzitdz_cpu, dcsidxs_cpu, detadys_cpu,&
        & dzitdzs_cpu, dcsidx2_cpu, detady2_cpu, dzitdz2_cpu, self%wallprop_cpu)
      elseif (mode == 1) then
        call visflx_div_subroutine(nx, ny, nz, ng, visc_order, self%w_aux_cpu, self%fl_cpu,&
        & coeff_deriv1_cpu, dcsidx_cpu, detady_cpu, dzitdz_cpu, stream1)
      elseif (mode == 2) then
        call visflx_reduced_ord2_subroutine(nx, ny, nz, ng, prandtl, t0, indx_cp_l, indx_cp_r,&
        & cp_coeff_cpu, calorically_perfect, u0, l0, self%base_cpu%w_cpu, self%w_aux_cpu, self%fl_cpu, x_cpu,&
        & y_cpu, z_cpu, self%wallprop_cpu, 1)
      elseif (mode == 3) then
        call sensor_subroutine(nx, ny, nz, ng, u0, l0, self%w_aux_cpu)
      elseif (mode == 4) then
        call visflx_nosensor_subroutine(nx, ny, nz, ng, visc_order, prandtl, t0, indx_cp_l,&
        & indx_cp_r, cp_coeff_cpu, calorically_perfect, u0, l0, self%base_cpu%w_cpu, self%w_aux_cpu,&
        & self%fl_cpu, coeff_deriv1_cpu, coeff_deriv2_cpu, dcsidx_cpu, detady_cpu, dzitdz_cpu, dcsidxs_cpu,&
        & detadys_cpu, dzitdzs_cpu, dcsidx2_cpu, detady2_cpu, dzitdz2_cpu, self%wallprop_cpu)
      elseif (mode == 5) then
        call visflx_x_subroutine(nx, ny, nz, nv, ng, prandtl, t0, indx_cp_l, indx_cp_r,&
        & cp_coeff_cpu, calorically_perfect, self%base_cpu%x_cpu, w_aux_cpu, self%fl_cpu, self%fhat_cpu)
        call visflx_y_subroutine(nx, ny, nz, nv, ng, prandtl, t0, indx_cp_l, indx_cp_r,&
        & cp_coeff_cpu, calorically_perfect, self%base_cpu%y_cpu, self%w_aux_cpu, self%fl_cpu)
        call visflx_z_subroutine(nx, ny, nz, nv, ng, prandtl, t0, indx_cp_l, indx_cp_r,&
        & cp_coeff_cpu, calorically_perfect, self%base_cpu%z_cpu, self%w_aux_cpu, self%fl_cpu)
      elseif (mode == 6) then
        call visflx_reduced_ord2_subroutine(nx, ny, nz, ng, prandtl, t0, indx_cp_l, indx_cp_r,&
        & cp_coeff_cpu, calorically_perfect, u0, l0, self%base_cpu%w_cpu, self%w_aux_cpu, self%fl_cpu, x_cpu,&
        & y_cpu, z_cpu, self%wallprop_cpu, 0)
      elseif (mode == 7) then
        call visflx_div_ord2_subroutine(nx, ny, nz, ng, self%w_aux_cpu, self%fl_cpu, x_cpu, y_cpu, z_cpu, stream1)
      endif
    endassociate
  endsubroutine visflx

  subroutine compute_aux(self, central, ghost)
    class(equation_singleideal_cpu_object), intent(inout) :: self
    integer, intent(in), optional :: central, ghost
    integer :: central_, ghost_
    central_ = 1 ; if(present(central)) central_ = central
    ghost_ = 1 ; if(present(ghost)) ghost_ = ghost
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, visc_model => self%visc_mo&
    &del, mu0 => self%mu0, t0 => self%equation_base%t0, t_ref_dim => self%t_ref_dim, sutherland_s => self&
    &%sutherland_s, powerlaw_vtexp => self%powerlaw_vtexp, cv_coeff_cpu => self%cv_coeff_cpu,&
    & cp_coeff_cpu => self%cp_coeff_cpu, indx_cp_l => self%equation_base%indx_cp_l, indx_cp_r => self%equ&
    &ation_base%indx_cp_r, calorically_perfect => self%equation_base%calorically_perfect,&
    & rgas0 => self%equation_base%rgas0)
      if(central_ == 1 .and. ghost_ == 1) then
        call eval_aux_subroutine(nx, ny, nz, ng, 1-ng, nx+ng, 1-ng, ny+ng, 1-ng, nz+ng,&
        & self%base_cpu%w_cpu, self%w_aux_cpu, visc_model, mu0, t0, sutherland_s, t_ref_dim, powerlaw_vtexp,&
        & visc_power, visc_sutherland, visc_no, cv_coeff_cpu, indx_cp_l, indx_cp_r, rgas0,&
        & calorically_perfect, tol_iter_nr,stream1)
      elseif(central_ == 1 .and. ghost_ == 0) then
        call eval_aux_subroutine(nx, ny, nz, ng, 1, nx, 1, ny, 1, nz, self%base_cpu%w_cpu,&
        & self%w_aux_cpu, visc_model, mu0, t0, sutherland_s, t_ref_dim, powerlaw_vtexp, visc_power,&
        & visc_sutherland, visc_no, cv_coeff_cpu, indx_cp_l, indx_cp_r, rgas0, calorically_perfect,&
        & tol_iter_nr,stream1)
      elseif(central_ == 0 .and. ghost_ == 1) then
        call eval_aux_subroutine(nx, ny, nz, ng, 1-ng, 0, 1-ng, ny+ng, 1-ng, nz+ng,&
        & self%base_cpu%w_cpu, self%w_aux_cpu, visc_model, mu0, t0, sutherland_s, t_ref_dim, powerlaw_vtexp,&
        & visc_power, visc_sutherland, visc_no, cv_coeff_cpu, indx_cp_l, indx_cp_r, rgas0,&
        & calorically_perfect, tol_iter_nr,stream1)
        call eval_aux_subroutine(nx, ny, nz, ng, nx+1, nx+ng, 1-ng, ny+ng, 1-ng, nz+ng,&
        & self%base_cpu%w_cpu, self%w_aux_cpu, visc_model, mu0, t0, sutherland_s, t_ref_dim, powerlaw_vtexp,&
        & visc_power, visc_sutherland, visc_no, cv_coeff_cpu, indx_cp_l, indx_cp_r, rgas0,&
        & calorically_perfect, tol_iter_nr,stream1)
        call eval_aux_subroutine(nx, ny, nz, ng, 1-ng, nx+ng, 1-ng, 0, 1-ng, nz+ng,&
        & self%base_cpu%w_cpu, self%w_aux_cpu, visc_model, mu0, t0, sutherland_s, t_ref_dim, powerlaw_vtexp,&
        & visc_power, visc_sutherland, visc_no, cv_coeff_cpu, indx_cp_l, indx_cp_r, rgas0,&
        & calorically_perfect, tol_iter_nr,stream1)
        call eval_aux_subroutine(nx, ny, nz, ng, 1-ng, nx+ng, ny+1, ny+ng, 1-ng, nz+ng,&
        & self%base_cpu%w_cpu, self%w_aux_cpu, visc_model, mu0, t0, sutherland_s, t_ref_dim, powerlaw_vtexp,&
        & visc_power, visc_sutherland, visc_no, cv_coeff_cpu, indx_cp_l, indx_cp_r, rgas0,&
        & calorically_perfect, tol_iter_nr,stream1)
        call eval_aux_subroutine(nx, ny, nz, ng, 1-ng, nx+ng, 1-ng, ny+ng, 1-ng, 0,&
        & self%base_cpu%w_cpu, self%w_aux_cpu, visc_model, mu0, t0, sutherland_s, t_ref_dim, powerlaw_vtexp,&
        & visc_power, visc_sutherland, visc_no, cv_coeff_cpu, indx_cp_l, indx_cp_r, rgas0,&
        & calorically_perfect, tol_iter_nr,stream1)
        call eval_aux_subroutine(nx, ny, nz, ng, 1-ng, nx+ng, 1-ng, ny+ng, nz+1, nz+ng,&
        & self%base_cpu%w_cpu, self%w_aux_cpu, visc_model, mu0, t0, sutherland_s, t_ref_dim, powerlaw_vtexp,&
        & visc_power, visc_sutherland, visc_no, cv_coeff_cpu, indx_cp_l, indx_cp_r, rgas0,&
        & calorically_perfect, tol_iter_nr,stream1)
      endif
    endassociate
  endsubroutine compute_aux

  subroutine rk_sync(self)
    class(equation_singleideal_cpu_object), intent(inout) :: self
    integer :: istep, lmax, iercuda
    real(rkind) :: rhodt, gamdt, alpdt
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv,&
    & dt => self%equation_base%dt, ep_order => self%equation_base%ep_order, weno_scheme => self%equation_&
    &base%weno_scheme, conservative_viscous => self%equation_base%conservative_viscous,&
    & eul_imin => self%equation_base%eul_imin, eul_imax => self%equation_base%eul_imax,&
    & eul_jmin => self%equation_base%eul_jmin, eul_jmax => self%equation_base%eul_jmax,&
    & eul_kmin => self%equation_base%eul_kmin, eul_kmax => self%equation_base%eul_kmax,&
    & channel_case => self%equation_base%channel_case, nv_aux => self%nv_aux)


      if (channel_case) self%equation_base%dpdx = 0._rkind

      do istep=1,self%equation_base%nrk
        rhodt = self%equation_base%rhork(istep)*dt
        gamdt = self%equation_base%gamrk(istep)*dt
        alpdt = self%equation_base%alprk(istep)*dt

        call init_flux_subroutine(nx, ny, nz, nv, self%fl_cpu, self%fln_cpu, rhodt)
        call self%base_cpu%bcswap()
        call self%compute_aux()

        call self%euler_x(eul_imin, eul_imax)

        if (conservative_viscous== 1) then
          call self%visflx(mode=6)
          call self%visflx(mode=5)
        else
          call self%visflx(mode=4)
        endif
        call bcextr_var_subroutine(nx, ny, nz, ng, self%w_aux_cpu(:,:,:,10:10))
        call self%base_cpu%bcswap_var(self%w_aux_cpu(:,:,:,10:10))
        call self%euler_y(eul_jmin,eul_jmax)

        call self%euler_z(eul_kmin,eul_kmax)
        if (istep == 3) then
          call self%visflx(mode=3)
          call self%base_cpu%bcswap_var(self%w_aux_cpu(:,:,:,8:8))
        endif
        if (conservative_viscous==1) then
          call self%visflx(mode=7)
        else
          call self%visflx(mode=1)
        endif

        call self%bc_nr()
        call update_flux_subroutine(nx, ny, nz, nv, self%fl_cpu, self%fln_cpu, gamdt)
        if (channel_case) call self%force_rhs()
        call update_field_subroutine(nx, ny, nz, ng, nv, self%base_cpu%w_cpu, self%fln_cpu, self%fluid_mask_cpu)
        if (channel_case) call self%force_var()
        call self%update_ghost(do_swap=0)
      enddo

    endassociate
  endsubroutine rk_sync

  subroutine rk_async(self)
    class(equation_singleideal_cpu_object), intent(inout) :: self
    integer :: istep, lmax, iercuda
    real(rkind) :: rhodt, gamdt, alpdt
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv,&
    & dt => self%equation_base%dt, ep_order => self%equation_base%ep_order, weno_scheme => self%equation_&
    &base%weno_scheme, conservative_viscous => self%equation_base%conservative_viscous,&
    & eul_imin => self%equation_base%eul_imin, eul_imax => self%equation_base%eul_imax,&
    & eul_jmin => self%equation_base%eul_jmin, eul_jmax => self%equation_base%eul_jmax,&
    & eul_kmin => self%equation_base%eul_kmin, eul_kmax => self%equation_base%eul_kmax,&
    & channel_case => self%equation_base%channel_case, nv_aux => self%nv_aux)



      if (channel_case) self%equation_base%dpdx = 0._rkind
      lmax = max(ep_order/2, weno_scheme)

      do istep=1,self%equation_base%nrk
        rhodt = self%equation_base%rhork(istep)*dt
        gamdt = self%equation_base%gamrk(istep)*dt
        alpdt = self%equation_base%alprk(istep)*dt

        call init_flux_subroutine(nx, ny, nz, nv, self%fl_cpu, self%fln_cpu, rhodt)


        call self%compute_aux(central=1, ghost=0)

        call self%base_cpu%bcswap(steps=[.true.,.false.,.false.])
        call self%euler_x(lmax+1,nx-lmax,lmax,nx-lmax,do_update=.false.)
        call self%base_cpu%bcswap(steps=[.false.,.true.,.true.])
        call self%compute_aux(central=0, ghost=1)
        call self%euler_x(eul_imin,lmax,eul_imin-1,lmax-1,do_update=.false.)
        call self%euler_x(nx-lmax+1,eul_imax,nx-lmax+1,eul_imax,do_update=.true.)
        if (conservative_viscous == 1) then
          call self%visflx(mode=6)
          call self%visflx(mode=5)
        else
          call self%visflx(mode=4)
        endif
        call bcextr_var_subroutine(nx, ny, nz, ng, self%w_aux_cpu(:,:,:,10:10))
        call self%base_cpu%bcswap_var(self%w_aux_cpu(:,:,:,10:10), steps=[.true.,.false.,.false.])
        call self%euler_y(eul_jmin,eul_jmax)
        call self%base_cpu%bcswap_var(self%w_aux_cpu(:,:,:,10:10), steps=[.false.,.true.,.true.])
        call self%euler_z(eul_kmin,eul_kmax)
        if (istep == 3) then
          call self%visflx(mode=3)
          call self%base_cpu%bcswap_var(self%w_aux_cpu(:,:,:,8:8), steps=[.true.,.false.,.false.])
          if (conservative_viscous==1) then
            call self%visflx(mode=7)
          else
            call self%visflx(mode=1)
          endif
          call self%base_cpu%bcswap_var(self%w_aux_cpu(:,:,:,8:8), steps=[.false.,.true.,.true.])
        else
          if (conservative_viscous==1) then
            call self%visflx(mode=7)
          else
            call self%visflx(mode=1)
          endif

        endif
        call self%bc_nr()
        call update_flux_subroutine(nx, ny, nz, nv, self%fl_cpu, self%fln_cpu, gamdt)
        if (channel_case) call self%force_rhs()
        call update_field_subroutine(nx, ny, nz, ng, nv, self%base_cpu%w_cpu, self%fln_cpu, self%fluid_mask_cpu)
        if (channel_case) call self%force_var()
        call self%update_ghost(do_swap=0)
      enddo
    endassociate
  endsubroutine rk_async

  subroutine euler_x(self, istart, iend, istart_face, iend_face, do_update)
    class(equation_singleideal_cpu_object), intent(inout) :: self
    integer, intent(in) :: istart, iend
    integer :: lmax, weno_size, iercuda, ierror
    logical, optional :: do_update
    integer, optional :: istart_face, iend_face
    integer :: istart_face_, iend_face_
    logical :: do_update_

    integer :: force_zero_flux_min,force_zero_flux_max

    do_update_ = .true. ; if(present(do_update)) do_update_ = do_update
    istart_face_ = istart-1 ; if(present(istart_face)) istart_face_ = istart_face
    iend_face_ = iend ; if(present(iend_face)) iend_face_ = iend_face

    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv,&
    & nv_aux => self%nv_aux, ep_order => self%equation_base%ep_order, force_zero_flux => self%equation_ba&
    &se%force_zero_flux, coeff_deriv1_cpu => self%coeff_deriv1_cpu, dcsidx_cpu => self%base_cpu%dcsidx_cp&
    &u, fhat_cpu => self%fhat_cpu, w_cpu => self%base_cpu%w_cpu, w_aux_cpu => self%w_aux_cpu,&
    & fl_cpu => self%fl_cpu, sensor_threshold => self%equation_base%sensor_threshold, weno_scheme => self&
    &%equation_base%weno_scheme, weno_version => self%equation_base%weno_version, cp_coeff_cpu => self%cp&
    &_coeff_cpu, indx_cp_l => self%equation_base%indx_cp_l, indx_cp_r => self%equation_base%indx_cp_r,&
    & calorically_perfect => self%equation_base%calorically_perfect, ep_ord_change_cpu => self%ep_ord_cha&
    &nge_cpu, nkeep => self%equation_base%nkeep, flux_splitting => self%equation_base%flux_splitting,&
    & eul_imin => self%equation_base%eul_imin, eul_imax => self%equation_base%eul_imax,&
    & rgas0 => self%equation_base%rgas0)

      if(iend - istart >= 0) then
        weno_size = 2*weno_scheme
        lmax = ep_order/2
        force_zero_flux_min = force_zero_flux(1)
        force_zero_flux_max = force_zero_flux(2)




        if (flux_splitting==1) then

          call euler_x_fluxes_hybrid_rusanov_kernel(nv, nv_aux, nx, ny, nz, ng, istart_face_,&
          & iend_face_, lmax, nkeep, rgas0, w_aux_cpu, coeff_deriv1_cpu, dcsidx_cpu, fhat_cpu,&
          & force_zero_flux_min, force_zero_flux_max, weno_scheme, weno_version, sensor_threshold, weno_size,&
          & cp_coeff_cpu, indx_cp_l, indx_cp_r, ep_ord_change_cpu, calorically_perfect, tol_iter_nr,&
          &self%equation_base%rho0,self%equation_base%u0, self%equation_base%t0)
        else

          call euler_x_fluxes_hybrid_kernel(nv, nv_aux, nx, ny, nz, ng, istart_face_, iend_face_,&
          & lmax, nkeep, rgas0, w_aux_cpu, coeff_deriv1_cpu, dcsidx_cpu, fhat_cpu, force_zero_flux_min,&
          & force_zero_flux_max, weno_scheme, weno_version, sensor_threshold, weno_size, cp_coeff_cpu,&
          & indx_cp_l, indx_cp_r, ep_ord_change_cpu, calorically_perfect, tol_iter_nr,self%equation_base%rho0,&
          &self%equation_base%u0, self%equation_base%t0)
        endif
      endif

      if(do_update_) then
        call euler_x_update_subroutine(nx, ny, nz, ng, nv, eul_imin, eul_imax, fhat_cpu, fl_cpu, dcsidx_cpu, stream1)
      endif

    endassociate
  endsubroutine euler_x

  subroutine euler_y(self, eul_jmin, eul_jmax)
    class(equation_singleideal_cpu_object), intent(inout) :: self
    integer, intent(in) :: eul_jmin, eul_jmax
    integer :: lmax, weno_size

    integer :: force_zero_flux_min,force_zero_flux_max

    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv,&
    & nv_aux => self%nv_aux, ep_order => self%equation_base%ep_order, force_zero_flux => self%equation_ba&
    &se%force_zero_flux, coeff_deriv1_cpu => self%coeff_deriv1_cpu, detady_cpu => self%base_cpu%detady_cp&
    &u, fhat_cpu => self%fhat_cpu, w_aux_cpu => self%w_aux_cpu, fl_cpu => self%fl_cpu,&
    & sensor_threshold => self%equation_base%sensor_threshold, weno_scheme => self%equation_base%weno_sch&
    &eme, weno_version => self%equation_base%weno_version, indx_cp_l => self%equation_base%indx_cp_l,&
    & indx_cp_r => self%equation_base%indx_cp_r, calorically_perfect => self%equation_base%calorically_pe&
    &rfect, cp_coeff_cpu => self%cp_coeff_cpu, ep_ord_change_cpu => self%ep_ord_change_cpu,&
    & nkeep => self%equation_base%nkeep, flux_splitting => self%equation_base%flux_splitting,&
    & rgas0 => self%equation_base%rgas0)
      weno_size = 2*weno_scheme
      lmax = ep_order/2
      force_zero_flux_min = force_zero_flux(3)
      force_zero_flux_max = force_zero_flux(4)



      if (flux_splitting==1) then
        call euler_y_hybrid_rusanov_kernel(nv, nv_aux, nx, ny, nz, ng, eul_jmin, eul_jmax, lmax,&
        & nkeep, rgas0, w_aux_cpu, fl_cpu, coeff_deriv1_cpu, detady_cpu, fhat_cpu, force_zero_flux_min,&
        & force_zero_flux_max, weno_scheme, weno_version, sensor_threshold, weno_size, cp_coeff_cpu,&
        & indx_cp_l, indx_cp_r, ep_ord_change_cpu, calorically_perfect, tol_iter_nr,self%equation_base%rho0,&
        &self%equation_base%u0, self%equation_base%t0)
      else
        call euler_y_hybrid_kernel(nv, nv_aux, nx, ny, nz, ng, eul_jmin, eul_jmax, lmax, nkeep,&
        & rgas0, w_aux_cpu, fl_cpu, coeff_deriv1_cpu, detady_cpu, fhat_cpu, force_zero_flux_min,&
        & force_zero_flux_max, weno_scheme, weno_version, sensor_threshold, weno_size, cp_coeff_cpu,&
        & indx_cp_l, indx_cp_r, ep_ord_change_cpu, calorically_perfect, tol_iter_nr,self%equation_base%rho0,&
        &self%equation_base%u0, self%equation_base%t0)
      endif
      call euler_y_update_subroutine(nx, ny, nz, ng, nv, eul_jmin, eul_jmax, fhat_cpu, fl_cpu, detady_cpu, stream1)
    endassociate
  endsubroutine euler_y

  subroutine euler_z(self, eul_kmin, eul_kmax)
    class(equation_singleideal_cpu_object), intent(inout) :: self
    integer, intent(in) :: eul_kmin, eul_kmax
    integer :: lmax, weno_size

    integer :: force_zero_flux_min, force_zero_flux_max

    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv,&
    & nv_aux => self%nv_aux, ep_order => self%equation_base%ep_order, force_zero_flux => self%equation_ba&
    &se%force_zero_flux, coeff_deriv1_cpu => self%coeff_deriv1_cpu, dzitdz_cpu => self%base_cpu%dzitdz_cp&
    &u, fhat_cpu => self%fhat_cpu, w_aux_cpu => self%w_aux_cpu, fl_cpu => self%fl_cpu,&
    & sensor_threshold => self%equation_base%sensor_threshold, weno_scheme => self%equation_base%weno_sch&
    &eme, weno_version => self%equation_base%weno_version, cp_coeff_cpu => self%cp_coeff_cpu,&
    & indx_cp_l => self%equation_base%indx_cp_l, indx_cp_r => self%equation_base%indx_cp_r,&
    & calorically_perfect => self%equation_base%calorically_perfect, ep_ord_change_cpu => self%ep_ord_cha&
    &nge_cpu, nkeep => self%equation_base%nkeep, flux_splitting => self%equation_base%flux_splitting,&
    & rgas0 => self%equation_base%rgas0, inflow_random_plane => self%equation_base%inflow_random_plane,&
    & inflow_random_plane_cpu => self%inflow_random_plane_cpu)
      weno_size = 2*weno_scheme
      lmax = ep_order/2
      force_zero_flux_min = force_zero_flux(5)
      force_zero_flux_max = force_zero_flux(6)


      if (flux_splitting==1) then
        call euler_z_hybrid_rusanov_kernel(nv, nv_aux, nx, ny, nz, ng, eul_kmin, eul_kmax, lmax,&
        & nkeep, rgas0, w_aux_cpu, fl_cpu, coeff_deriv1_cpu, dzitdz_cpu, fhat_cpu, force_zero_flux_min,&
        & force_zero_flux_max, weno_scheme, weno_version, sensor_threshold, weno_size, cp_coeff_cpu,&
        & indx_cp_l, indx_cp_r, ep_ord_change_cpu, calorically_perfect, tol_iter_nr,self%equation_base%rho0,&
        &self%equation_base%u0, self%equation_base%t0)
      else
        call euler_z_hybrid_kernel(nv, nv_aux, nx, ny, nz, ng, eul_kmin, eul_kmax, lmax, nkeep,&
        & rgas0, w_aux_cpu, fl_cpu, coeff_deriv1_cpu, dzitdz_cpu, fhat_cpu, force_zero_flux_min,&
        & force_zero_flux_max, weno_scheme, weno_version, sensor_threshold, weno_size, cp_coeff_cpu,&
        & indx_cp_l, indx_cp_r, ep_ord_change_cpu, calorically_perfect, tol_iter_nr,self%equation_base%rho0,&
        &self%equation_base%u0, self%equation_base%t0)
      endif
      call euler_z_update_subroutine(nx, ny, nz, ng, nv, eul_kmin, eul_kmax, fhat_cpu, fl_cpu, dzitdz_cpu, 0)
      if (self%equation_base%recyc) call get_crandom_f(inflow_random_plane(2:self%equation_base%jbl_inflow,1:nz,1:3))
    endassociate
  endsubroutine euler_z

  subroutine force_rhs(self)
    class(equation_singleideal_cpu_object), intent(inout) :: self
    real(rkind), dimension(5) :: bulk5, bulk5g

    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nxmax => self%equation_bas&
    &e%grid%nxmax, nzmax => self%equation_base%grid%nzmax, yn => self%equation_base%field%yn,&
    & yn_cpu => self%base_cpu%yn_cpu, w_cpu => self%base_cpu%w_cpu, w_aux_cpu => self%w_aux_cpu,&
    & fln_cpu => self%fln_cpu, volchan => self%equation_base%volchan, fluid_mask_cpu => self%fluid_mask_c&
    &pu)

      call force_rhs_1_subroutine(nx, ny, nz, ng, yn_cpu, fln_cpu, w_cpu, w_aux_cpu, bulk5, fluid_mask_cpu)

      call mpi_allreduce(bulk5,bulk5g,5,mpi_prec,mpi_sum,self%equation_base%field%mp_cart,self%mpi_err)

      bulk5g = bulk5g/volchan
      self%equation_base%dpdx = self%equation_base%dpdx+bulk5g(2)
      self%equation_base%rhobulk = bulk5g(3)
      self%equation_base%ubulk = bulk5g(4)/self%equation_base%rhobulk
      self%equation_base%tbulk = bulk5g(5)/self%equation_base%rhobulk/self%equation_base%ubulk

      call force_rhs_2_subroutine(nx, ny, nz, ng, fln_cpu, w_aux_cpu, bulk5g, fluid_mask_cpu)
    endassociate
  endsubroutine force_rhs
  subroutine force_var(self)
    class(equation_singleideal_cpu_object), intent(inout) :: self
    real(rkind) :: bulkt, bulktg, tbtarget, tbdiff

    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nxmax => self%equation_bas&
    &e%grid%nxmax, nzmax => self%equation_base%grid%nzmax, yn => self%equation_base%field%yn,&
    & yn_cpu => self%base_cpu%yn_cpu, w_cpu => self%base_cpu%w_cpu, w_aux_cpu => self%w_aux_cpu,&
    & fln_cpu => self%fln_cpu, volchan => self%equation_base%volchan, fluid_mask_cpu => self%fluid_mask_c&
    &pu, cv_coeff_cpu => self%cv_coeff_cpu, t0 => self%equation_base%t0, calorically_perfect => self%equa&
    &tion_base%calorically_perfect, indx_cp_l => self%equation_base%indx_cp_l, indx_cp_r => self%equation&
    &_base%indx_cp_r)

      if (self%equation_base%theta_wall>=-1._rkind) then
        call force_var_1_subroutine(nx, ny, nz, ng, yn_cpu, fln_cpu, w_cpu, w_aux_cpu, bulkt,&
        & fluid_mask_cpu, cv_coeff_cpu, indx_cp_l, indx_cp_r, t0, calorically_perfect, tol_iter_nr)
        call mpi_allreduce(bulkt,bulktg,1,mpi_prec,mpi_sum,self%equation_base%field%mp_cart,self%mpi_err)

        bulktg = bulktg/volchan
        bulktg = bulktg/self%equation_base%rhobulk/self%equation_base%ubulk
        tbtarget = self%equation_base%t_bulk_target
        tbdiff = tbtarget-bulktg
        call force_var_2_subroutine(nx, ny, nz, ng, w_cpu, w_aux_cpu, tbdiff, fluid_mask_cpu,&
        & cv_coeff_cpu, indx_cp_l, indx_cp_r, t0, calorically_perfect)
      endif
    endassociate
  endsubroutine force_var

  subroutine recyc_exchange(self)
    class(equation_singleideal_cpu_object), intent(inout) :: self

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
    & nblocks => self%equation_base%field%nblocks, iermpi => self%mpi_err, w_cpu => self%base_cpu%w_cpu,&
    & nv => self%equation_base%nv, wbuf1s_cpu => self%base_cpu%wbuf1s_cpu, wbuf2r_cpu => self%base_cpu%wb&
    &uf2r_cpu, wbuf1r_cpu => self%base_cpu%wbuf1r_cpu, wrecyc_cpu => self%wrecyc_cpu, ibrecyc => self%equ&
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
        call recyc_exchange_subroutine_1(irecyc, w_cpu, wbuf1s_cpu, nx, ny, nz, ng, nv)
        indx = nv*ng*ny*nz
        call mpi_isend(wbuf1s_cpu,indx,mpi_prec,sendto1,2000,mp_cart,req(1),iermpi)
        call mpi_isend(wbuf1s_cpu,indx,mpi_prec,sendto2,3000,mp_cart,req(2),iermpi)
      endif
      if (ncoords(1)==0) then
        kcoordrecvfrom1 = ncoords(3)-kshiftglob/nz+nblocks(3)
        kcoordrecvfrom2 = kcoordrecvfrom1-1
        kcoordrecvfrom1 = mod(kcoordrecvfrom1,nblocks(3))
        kcoordrecvfrom2 = mod(kcoordrecvfrom2,nblocks(3))
        call mpi_cart_rank(mp_cart,[ibrecyc,0,kcoordrecvfrom1],recvfrom1,iermpi)
        call mpi_cart_rank(mp_cart,[ibrecyc,0,kcoordrecvfrom2],recvfrom2,iermpi)
        indx = nv*ng*ny*nz
        call mpi_irecv(wbuf1r_cpu,indx,mpi_prec,recvfrom1,2000,mp_cart,req(3),iermpi)
        call mpi_irecv(wbuf2r_cpu,indx,mpi_prec,recvfrom2,3000,mp_cart,req(4),iermpi)
      endif
      call mpi_waitall(4,req,mpi_statuses_ignore,iermpi)
      if (ncoords(1)==0) then
        call recyc_exchange_subroutine_2(n1_start_recv, n1_start_send, n1_end_recv, wrecyc_cpu, wbuf1r_cpu, nx, ny, nz, ng, nv)
        call recyc_exchange_subroutine_3(n2_start_recv, n2_start_send, n2_end_recv, wrecyc_cpu, wbuf2r_cpu, nx, ny, nz, ng, nv)
      endif
    endassociate

  end subroutine recyc_exchange

  subroutine bc_nr(self)
    class(equation_singleideal_cpu_object), intent(inout) :: self
    integer :: ilat, dir, start_or_end


    associate(bctags_nr => self%equation_base%bctags_nr, nx => self%nx, ny => self%ny,&
    & nz => self%nz, nv => self%nv, ng => self%ng, w_aux_cpu => self%w_aux_cpu, w_cpu => self%base_cpu%w_&
    &cpu, fl_cpu => self%fl_cpu, dcsidx_cpu => self%base_cpu%dcsidx_cpu, detady_cpu => self%base_cpu%deta&
    &dy_cpu, dzitdz_cpu => self%base_cpu%dzitdz_cpu, indx_cp_l => self%equation_base%indx_cp_l,&
    & indx_cp_r => self%equation_base%indx_cp_r, calorically_perfect => self%equation_base%calorically_pe&
    &rfect, rgas0 => self%equation_base%rgas0, t0 => self%equation_base%t0, cp_coeff_cpu => self%cp_coeff&
    &_cpu, winf_cpu => self%winf_cpu)

      do ilat=1,6
        if(bctags_nr(ilat) > 0) then
          dir = (ilat-1)/2 +1
          start_or_end = mod(ilat-1,2)+1
          if(dir == 1) then


            call bc_nr_lat_x_kernel(start_or_end, bctags_nr(ilat), nx, ny, nz, ng, nv, w_aux_cpu,&
            & w_cpu, fl_cpu, dcsidx_cpu, indx_cp_l, indx_cp_r, cp_coeff_cpu, winf_cpu, calorically_perfect,rgas0,&
            &t0)
          endif
          if(dir == 2) then


            call bc_nr_lat_y_kernel(start_or_end, bctags_nr(ilat), nx, ny, nz, ng, nv, w_aux_cpu,&
            & w_cpu, fl_cpu, detady_cpu, indx_cp_l, indx_cp_r, cp_coeff_cpu, winf_cpu, calorically_perfect,rgas0,&
            &t0)
          endif
          if(dir == 3) then


            call bc_nr_lat_z_kernel(start_or_end, bctags_nr(ilat), nx, ny, nz, ng, nv, w_aux_cpu,&
            & w_cpu, fl_cpu, dzitdz_cpu, indx_cp_l, indx_cp_r, cp_coeff_cpu, winf_cpu, calorically_perfect,rgas0,&
            &t0)
          endif
        endif
      enddo
    endassociate

  endsubroutine bc_nr

  subroutine update_ghost(self, do_swap)
    class(equation_singleideal_cpu_object), intent(inout) :: self
    integer, intent(in), optional :: do_swap
    integer :: do_swap_
    integer :: ilat

    do_swap_ = 1 ; if (present(do_swap)) do_swap_ = do_swap

    if (self%equation_base%recyc) call self%recyc_exchange()

    do ilat=1,6
      select case(self%equation_base%bctags(ilat))
      case(0)
      case(1)
        call bcfree_subroutine(ilat, self%nx, self%ny, self%nz, self%ng, self%nv, self%winf_cpu, self%base_cpu%w_cpu)
      case(2)
        call bcextr_subroutine(ilat, self%nx, self%ny, self%nz, self%ng, self%nv, self%base_cpu%w_cpu)
      case(4)
        call bcextr_sub_subroutine(ilat, self%nx, self%ny, self%nz, self%ng, self%equation_base%p0,&
        & self%equation_base%rgas0, self%base_cpu%w_cpu, self%equation_base%indx_cp_l, self%equation_base%ind&
        &x_cp_r, self%cv_coeff_cpu, self%equation_base%t0, self%equation_base%calorically_perfect)
      case(5)
        call bcsym_subroutine(ilat, self%nx, self%ny, self%nz, self%ng, self%equation_base%t_wall,&
        & self%base_cpu%w_cpu, self%w_aux_cpu, self%equation_base%indx_cp_l, self%equation_base%indx_cp_r,&
        & self%cv_coeff_cpu, self%equation_base%calorically_perfect)
      case(6)
        if (self%equation_base%channel_case) then
          call bcwall_staggered_subroutine(ilat, self%nx, self%ny, self%nz, self%ng,&
          & self%equation_base%t_wall, self%base_cpu%w_cpu, self%w_aux_cpu, self%equation_base%indx_cp_l,&
          & self%equation_base%indx_cp_r, self%cv_coeff_cpu, self%equation_base%t0, self%equation_base%rgas0,&
          & self%equation_base%calorically_perfect,tol_iter_nr)
        else
          call bcwall_subroutine(ilat, self%nx, self%ny, self%nz, self%ng, self%equation_base%t_wall&
          &, self%base_cpu%w_cpu, self%w_aux_cpu, self%equation_base%indx_cp_l, self%equation_base%indx_cp_r,&
          & self%cv_coeff_cpu, self%equation_base%t0, self%equation_base%rgas0, self%equation_base%calorically_&
          &perfect,tol_iter_nr)
        endif
      case(7)
        call bcshock_subroutine(ilat, self%nx, self%ny, self%nz, self%ng, self%nv,&
        & self%base_cpu%w_cpu, self%winf_cpu, self%winf_past_shock_cpu, self%equation_base%xshock_imp,&
        & self%equation_base%shock_angle, self%base_cpu%x_cpu, self%base_cpu%y_cpu, self%equation_base%tanhfa&
        &cs)
      case(8)
      case(9)
        call bclam_subroutine(ilat, self%nx, self%ny, self%nz, self%ng, self%nv, self%base_cpu%w_cpu&
        &, self%wmean_cpu, self%equation_base%p0, self%equation_base%rgas0, self%equation_base%indx_cp_l,&
        & self%equation_base%indx_cp_r, self%cv_coeff_cpu, self%equation_base%t0, self%equation_base%calorica&
        &lly_perfect)
      case(10)
        call self%bcrecyc(ilat)
      endselect
    enddo

    if (do_swap_ == 1) call self%base_cpu%bcswap()

  endsubroutine update_ghost

  subroutine bcrecyc(self, ilat)
    class(equation_singleideal_cpu_object), intent(inout) :: self
    integer, intent(in) :: ilat
    integer :: ntot
    if (ilat == 1) then
      associate(nx => self%nx, ny => self%ny, nz => self%nz, nzmax => self%grid%nzmax,&
      & ng => self%ng, nv => self%nv, wrecycav_cpu => self%wrecycav_cpu, wrecyc_cpu => self%wrecyc_cpu,&
      & mp_cartz => self%field%mp_cartz, iermpi => self%mpi_err, p0 => self%equation_base%p0,&
      & w_cpu => self%base_cpu%w_cpu, wmean_cpu => self%wmean_cpu, weta_inflow_cpu => self%weta_inflow_cpu,&
      & map_j_inn_cpu => self%map_j_inn_cpu, map_j_out_cpu => self%map_j_out_cpu, map_j_out_blend_cpu => se&
      &lf%map_j_out_blend_cpu, eta_recyc_blend_cpu => self%eta_recyc_blend_cpu, yplus_inflow_cpu => self%yp&
      &lus_inflow_cpu, eta_inflow_cpu => self%eta_inflow_cpu, yplus_recyc_cpu => self%yplus_recyc_cpu,&
      & eta_recyc_cpu => self%eta_recyc_cpu, betarecyc => self%equation_base%betarecyc, i_recyc => self%equ&
      &ation_base%i_recyc, glund1 => self%equation_base%glund1, indx_cp_l => self%equation_base%indx_cp_l,&
      & indx_cp_r => self%equation_base%indx_cp_r, rgas0 => self%equation_base%rgas0, calorically_perfect =&
      &> self%equation_base%calorically_perfect, rand_type => self%equation_base%rand_type,&
      & t0 => self%equation_base%t0, u0 => self%equation_base%u0, l0 => self%equation_base%l0,&
      & cv_coeff_cpu => self%cv_coeff_cpu, inflow_random_plane => self%equation_base%inflow_random_plane,&
      & inflow_random_plane_cpu => self%inflow_random_plane_cpu)
        call bcrecyc_subroutine_1(nx, ny, nz, ng, nv, wrecycav_cpu, wrecyc_cpu)

        inflow_random_plane_cpu = inflow_random_plane

        ntot = ng*ny*nv
        call mpi_allreduce(mpi_in_place,wrecycav_cpu,ntot,mpi_prec,mpi_sum,mp_cartz,iermpi)

        call bcrecyc_subroutine_2(nx, ny, nz, nzmax, ng, wrecycav_cpu, wrecyc_cpu)

        call bcrecyc_subroutine_3(nx, ny, nz, ng, p0, u0, rgas0, w_cpu, wmean_cpu, wrecyc_cpu,&
        & weta_inflow_cpu, map_j_inn_cpu, map_j_out_cpu, map_j_out_blend_cpu, yplus_inflow_cpu,&
        & eta_inflow_cpu, yplus_recyc_cpu, eta_recyc_cpu, eta_recyc_blend_cpu, betarecyc, glund1,&
        & inflow_random_plane_cpu, indx_cp_l, indx_cp_r, cv_coeff_cpu, t0, calorically_perfect,rand_type)
      endassociate
    endif

  end subroutine bcrecyc

  subroutine initialize(self, filename)
    class(equation_singleideal_cpu_object), intent(inout) :: self
    character(*) , intent(in) :: filename

    call get_mpi_basic_info(self%nprocs, self%myrank, self%masterproc, self%mpi_err)

    call self%equation_base%initialize(filename)
    self%nx = self%equation_base%field%nx
    self%ny = self%equation_base%field%ny
    self%nz = self%equation_base%field%nz
    self%ng = self%equation_base%grid%ng
    self%nv = self%equation_base%nv

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

    call self%base_cpu%initialize(self%equation_base%field)

    if(self%equation_base%debug_memory>0) then
      call self%field%check_cpu_mem(description="--initialize-start--")
    endif

    call self%base_cpu%copy_from_field()

    if(self%equation_base%debug_memory>0) then
      call self%field%check_cpu_mem(description="--initialize-first-gpu-usage--")
    endif

    self%nv_aux = 10

    call self%alloc()

    self%w_aux_cpu = self%equation_base%w_aux
    self%fluid_mask_cpu = self%equation_base%fluid_mask
    self%fluid_mask_ini_cpu = self%fluid_mask_cpu
    self%ep_ord_change_cpu = self%equation_base%ep_ord_change
    self%winf_cpu = self%equation_base%winf
    self%winf_past_shock_cpu = self%equation_base%winf_past_shock

    self%cp_coeff_cpu = self%equation_base%cp_coeff
    self%cv_coeff_cpu = self%equation_base%cv_coeff

    if (self%equation_base%num_probe>0) then
      self%probe_coeff_cpu = self%equation_base%probe_coeff
      self%ijk_probe_cpu = self%equation_base%ijk_probe
    endif

    self%wmean_cpu = self%equation_base%wmean

    allocate(self%coeff_deriv1_cpu(1:4,4))
    allocate(self%coeff_deriv2_cpu(0:4,4))
    self%coeff_deriv1_cpu = self%equation_base%coeff_deriv1
    self%coeff_deriv2_cpu = self%equation_base%coeff_deriv2

    if (self%equation_base%recyc) then
      self%yplus_inflow_cpu = self%equation_base%yplus_inflow
      self%eta_inflow_cpu = self%equation_base%eta_inflow
      self%yplus_recyc_cpu = self%equation_base%yplus_recyc
      self%eta_recyc_cpu = self%equation_base%eta_recyc
      self%eta_recyc_blend_cpu = self%equation_base%eta_recyc_blend
      self%map_j_inn_cpu = self%equation_base%map_j_inn
      self%map_j_out_cpu = self%equation_base%map_j_out
      self%map_j_out_blend_cpu = self%equation_base%map_j_out_blend
      self%weta_inflow_cpu = self%equation_base%weta_inflow
    endif


    if (self%equation_base%debug_memory>0) then
      call self%field%check_cpu_mem(description="initialize-completed")
    endif



  endsubroutine initialize

  subroutine alloc(self)
    class(equation_singleideal_cpu_object), intent(inout) :: self
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv,&
    & nv_aux => self%nv_aux, weno_scheme => self%equation_base%weno_scheme,calorically_perfect => self%eq&
    &uation_base%calorically_perfect, indx_cp_l => self%equation_base%indx_cp_l ,indx_cp_r => self%equati&
    &on_base%indx_cp_r)

      allocate(self%winf_cpu(nv))
      allocate(self%winf_past_shock_cpu(nv))
      allocate(self%w_aux_cpu(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng, nv_aux))
      allocate(self%fl_cpu(1:nx, 1:ny, 1:nz, nv))
      allocate(self%fln_cpu(1:nx, 1:ny, 1:nz, nv))
      allocate(self%wallprop_cpu(1-ng:nx+ng, 1-ng:nz+ng, 2:4))
      allocate(self%w_var(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng, 1))
      allocate(self%w_var_t(1, 1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng))
      allocate(self%fhat_cpu(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng, nv))
      allocate(self%fluid_mask_cpu(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng))
      allocate(self%fluid_mask_ini_cpu(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng))
      allocate(self%ep_ord_change_cpu(0:nx, 0:ny, 0:nz, 1:3))

      allocate(self%wrecyc_cpu(ng,ny,nz,nv))
      allocate(self%wrecycav_cpu(ng,ny,nv))
      allocate(self%wmean_cpu(1-self%grid%ng:self%field%nx+self%grid%ng+1, 1:self%field%ny,4))

      allocate(self%yplus_inflow_cpu(1-ng:ny+ng))
      allocate(self%eta_inflow_cpu(1-ng:ny+ng))
      allocate(self%yplus_recyc_cpu(1-ng:ny+ng))
      allocate(self%eta_recyc_cpu(1-ng:ny+ng))
      allocate(self%eta_recyc_blend_cpu(1-ng:ny+ng))
      allocate(self%map_j_inn_cpu(1:ny))
      allocate(self%map_j_out_cpu(1:ny))
      allocate(self%map_j_out_blend_cpu(1:ny))
      allocate(self%weta_inflow_cpu(1:ny))
      allocate(self%inflow_random_plane_cpu(1:ny,1:nz,3))

      allocate(self%cv_coeff_cpu(indx_cp_l:indx_cp_r+1))
      allocate(self%cp_coeff_cpu(indx_cp_l:indx_cp_r+1))

      if (self%equation_base%num_probe>0) then
        allocate(self%w_aux_probe_cpu(6,self%equation_base%num_probe))
        allocate(self%ijk_probe_cpu(3,self%equation_base%num_probe))
        allocate(self%probe_coeff_cpu(2,2,2,self%equation_base%num_probe))
      endif

    endassociate
  endsubroutine alloc

  subroutine compute_residual(self)
    class(equation_singleideal_cpu_object), intent(inout) :: self
    real(rkind) :: restemp
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv,&
    & dt => self%equation_base%dt, fln_cpu => self%fln_cpu, residual_rhou => self%equation_base%residual_&
    &rhou, fluid_mask_cpu => self%fluid_mask_cpu, vmax => self%equation_base%vmax, w_aux_cpu => self%w_au&
    &x_cpu, rgas0 => self%equation_base%rgas0, rhomin => self%equation_base%rhomin, rhomax => self%equati&
    &on_base%rhomax, pmin => self%equation_base%pmin, pmax => self%equation_base%pmax,&
    & tmin => self%equation_base%tmin, tmax => self%equation_base%tmax)

      call compute_residual_subroutine(nx, ny, nz, ng, nv, fln_cpu, dt, residual_rhou, fluid_mask_cpu)
      call mpi_allreduce(mpi_in_place,residual_rhou,1,mpi_prec,mpi_sum,self%equation_base%field%mp_cart,self%mpi_err)
      restemp = real(nx,rkind)*real(ny,rkind)*real(nz,rkind)*real(self%nprocs,rkind)
      residual_rhou = sqrt(residual_rhou/restemp)

      call compute_rho_t_p_minmax_subroutine(nx, ny, nz, ng, rgas0, w_aux_cpu,rhomin,rhomax, tmin,&
      & tmax, pmin, pmax, fluid_mask_cpu)
      call mpi_allreduce(mpi_in_place,rhomin,1,mpi_prec,mpi_min,self%equation_base%field%mp_cart,self%mpi_err)
      call mpi_allreduce(mpi_in_place,rhomax,1,mpi_prec,mpi_max,self%equation_base%field%mp_cart,self%mpi_err)
      call mpi_allreduce(mpi_in_place,tmin,1,mpi_prec,mpi_min,self%equation_base%field%mp_cart,self%mpi_err)
      call mpi_allreduce(mpi_in_place,tmax,1,mpi_prec,mpi_max,self%equation_base%field%mp_cart,self%mpi_err)
      call mpi_allreduce(mpi_in_place,pmin,1,mpi_prec,mpi_min,self%equation_base%field%mp_cart,self%mpi_err)
      call mpi_allreduce(mpi_in_place,pmax,1,mpi_prec,mpi_max,self%equation_base%field%mp_cart,self%mpi_err)

    endassociate
  endsubroutine compute_residual

  subroutine compute_dt(self)
    class(equation_singleideal_cpu_object), intent(inout) :: self
    real(rkind) :: dt_min, dtinv_max
    real(rkind) :: dtxi_max, dtyi_max, dtzi_max, dtxv_max, dtyv_max, dtzv_max, dtxk_max, dtyk_max, dtzk_max
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, dt => self%equation_base%d&
    &t, cfl => self%equation_base%cfl, prandtl => self%equation_base%prandtl, visc_model => self%visc_mod&
    &el, mu0 => self%mu0, t_ref_dim => self%t_ref_dim, powerlaw_vtexp => self%powerlaw_vtexp,&
    & sutherland_s => self%sutherland_s, w_cpu => self%base_cpu%w_cpu, w_aux_cpu => self%w_aux_cpu,&
    & dcsidx_cpu => self%base_cpu%dcsidx_cpu, detady_cpu => self%base_cpu%detady_cpu, dzitdz_cpu => self%&
    &base_cpu%dzitdz_cpu, dcsidxs_cpu => self%base_cpu%dcsidxs_cpu, detadys_cpu => self%base_cpu%detadys_&
    &cpu, dzitdzs_cpu => self%base_cpu%dzitdzs_cpu, calorically_perfect => self%equation_base%calorically&
    &_perfect, indx_cp_l => self%equation_base%indx_cp_l, indx_cp_r => self%equation_base%indx_cp_r,&
    & cp_coeff_cpu => self%cp_coeff_cpu, fluid_mask_cpu => self%fluid_mask_cpu, rgas0 => self%equation_ba&
    &se%rgas0, t0 => self%equation_base%t0)
      if (cfl < 0) then
        dt = -cfl
      else
        call compute_dt_subroutine(nx, ny, nz, ng, rgas0, prandtl, dcsidx_cpu, detady_cpu,&
        & dzitdz_cpu, dcsidxs_cpu, detadys_cpu, dzitdzs_cpu, w_cpu, w_aux_cpu, dtxi_max, dtyi_max, dtzi_max,&
        & dtxv_max, dtyv_max, dtzv_max, dtxk_max, dtyk_max, dtzk_max, indx_cp_l, indx_cp_r, cp_coeff_cpu,&
        &fluid_mask_cpu,calorically_perfect,t0)
        dtinv_max = maxval([dtxi_max, dtyi_max, dtzi_max, dtxv_max, dtyv_max, dtzv_max, dtxk_max, dtyk_max, dtzk_max])
        call mpi_allreduce(mpi_in_place,dtinv_max,1,mpi_prec,mpi_max,self%equation_base%field%mp_cart,self%mpi_err)
        dt_min = 1._rkind/dtinv_max
        dt = self%equation_base%cfl*dt_min
      endif
    endassociate
  endsubroutine compute_dt

  subroutine run(self, filename)

    class(equation_singleideal_cpu_object), intent(inout) :: self
    character(*) , intent(in) :: filename
    real(rkind) :: timing(1:2)
    real(rkind) :: timing_step(1:2)
    integer :: icyc_loop, iercuda

    call self%initialize(filename=filename)

    associate(icyc0 => self%equation_base%icyc0, icyc => self%equation_base%icyc,&
    & time => self%equation_base%time, iter_dt_recompute => self%equation_base%iter_dt_recompute,&
    & residual_rhou => self%equation_base%residual_rhou, dpdx => self%equation_base%dpdx,&
    & rhobulk => self%equation_base%rhobulk, ubulk => self%equation_base%ubulk, tbulk => self%equation_ba&
    &se%tbulk, nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, visc_model => self%visc_model,&
    & mu0 => self%mu0, t_ref_dim => self%t_ref_dim, sutherland_s => self%sutherland_s,&
    & powerlaw_vtexp => self%powerlaw_vtexp, cv_coeff_cpu => self%cv_coeff_cpu, calorically_perfect => se&
    &lf%equation_base%calorically_perfect, indx_cp_l => self%equation_base%indx_cp_l, indx_cp_r => self%e&
    &quation_base%indx_cp_r, mode_async => self%equation_base%mode_async, time_from_last_rst => self%equa&
    &tion_base%time_from_last_rst, time_from_last_write => self%equation_base%time_from_last_write,&
    & time_from_last_stat => self%equation_base%time_from_last_stat, time_from_last_slice => self%equatio&
    &n_base%time_from_last_slice, time_from_last_probe => self%equation_base%time_from_last_probe,&
    & time_is_freezed => self%equation_base%time_is_freezed)

      call self%update_ghost()
      call self%compute_aux()


      if (mode_async >= 0) then
        if (self%equation_base%conservative_viscous==1) then
          call self%visflx(mode=2)
        else
          call self%visflx(mode=0)
        endif
        call self%base_cpu%bcswap_var(self%w_aux_cpu(:,:,:,8:8 ))
      endif
      call zero_flux_subroutine(self%nx,self%ny,self%nz,self%nv,self%fl_cpu)

      if (self%equation_base%restart_type==0) then
        self%equation_base%w_aux = self%w_aux_cpu
        if (self%equation_base%enable_plot3d>0) then
          call self%field%write_plot3d(mach=self%equation_base%mach, reynolds=self%equation_base%rey&
          &nolds, time=0._rkind, istore=0, plot3dgrid=.true., plot3dfield=.true., w_aux_io=self%equation_base%w&
          &_aux(1:self%nx,1:self%ny,1:self%nz,1:6))
        endif
        if (self%equation_base%enable_vtk>0) then
          call self%field%write_vtk(time=0._rkind, istore=0, w_aux_io=self%equation_base%w_aux(1:self%nx,1:self%ny,1:self%nz,1:6))
        endif
      endif
      call self%compute_dt()
      if (self%masterproc) write(*,*) 'dt =', self%equation_base%dt
      if (self%masterproc) write(*,*) 'dt*u0/l0 =', self%equation_base%dt*self%equation_base%u0/self%equation_base%l0

      call mpi_barrier(mpi_comm_world, self%error) ; timing(1) = mpi_wtime()

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
      icyc_loop = icyc
      integration: do

        icyc_loop = icyc_loop + 1
        time_is_freezed = self%equation_base%time_is_freezed_fun()

        if ( time_is_freezed ) then
          self%equation_base%dt = 0._rkind
        else
          icyc = icyc + 1

          if(mod(icyc-icyc0, iter_dt_recompute)==0) then
            call self%compute_aux(central=1, ghost=0)

            call self%compute_dt()
          endif

          if(mode_async == -1) call self%rk_sync_old()
          if(mode_async == 0) call self%rk_sync()
          if(mode_async == 1) call self%rk_async()
          if (mod(icyc-icyc0, self%equation_base%print_control)==0) call self%compute_residual()
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

      call self%base_cpu%copy_to_field()
      if (self%equation_base%io_type_w==1) then
        call self%field%write_field_serial()
        call self%equation_base%write_stats_serial()
        if (self%equation_base%enable_stat_3d>0) call self%equation_base%write_stats_3d_serial()
      endif
      if (self%equation_base%io_type_w==2) then
        call self%field%write_field()
        call self%equation_base%write_stats()
        if (self%equation_base%enable_stat_3d>0) call self%equation_base%write_stats_3d()
      endif
      call self%equation_base%write_field_info()

    endassociate
  endsubroutine run

  subroutine print_progress(self)
    class(equation_singleideal_cpu_object), intent(inout) :: self
    character(6) :: pos_io
    associate(icyc => self%equation_base%icyc, time => self%equation_base%time, dt => self%equation_&
    &base%dt, residual_rhou => self%equation_base%residual_rhou, dpdx => self%equation_base%dpdx,&
    & vmax => self%equation_base%vmax, rhobulk => self%equation_base%rhobulk, ubulk => self%equation_base&
    &%ubulk, tbulk => self%equation_base%tbulk, rhomin => self%equation_base%rhomin, rhomax => self%equat&
    &ion_base%rhomax, pmin => self%equation_base%pmin, pmax => self%equation_base%pmax,&
    & tmin => self%equation_base%tmin, tmax => self%equation_base%tmax)
      residual_rhou = residual_rhou/(self%equation_base%u0**2)/self%equation_base%rho0*self%equation_base%l0
      dpdx = dpdx /(self%equation_base%u0**2)/self%equation_base%rho0*self%equation_base%l0
      pos_io = 'append'
      if (self%equation_base%icyc==1) pos_io = 'rewind'
      if (self%masterproc) then
        open(unit=15,file='progress.out',position=pos_io)
        if (self%equation_base%channel_case) then
          write(* ,'(1i10,20es20.10)') icyc, dt, time, residual_rhou, -dpdx/dt, rhobulk, ubulk,&
          & tbulk, rhomin, rhomax, tmin, tmax, pmin, pmax
          write(15,'(1i10,20es20.10)') icyc, dt, time, residual_rhou, -dpdx/dt, rhobulk, ubulk,&
          & tbulk, rhomin, rhomax, tmin, tmax, pmin, pmax
        else
          write(* ,'(1i10,20es20.10)') icyc, dt, time, residual_rhou, rhomin, rhomax, tmin, tmax, pmin, pmax
          write(15,'(1i10,20es20.10)') icyc, dt, time, residual_rhou, rhomin, rhomax, tmin, tmax, pmin, pmax
        endif
        close(15)
      endif
    endassociate
  endsubroutine print_progress







  subroutine manage_output(self)
    class(equation_singleideal_cpu_object), intent(inout) :: self
    integer :: i,j,k,ii,jj,kk,l,n
    integer :: isize, jsize, ksize
    character(3) :: chx, chy, chz
    logical :: sliceyz_exist, slicexz_exist, slicexy_exist, probe_exist
    associate(time_from_last_rst => self%equation_base%time_from_last_rst, time_from_last_write => s&
    &elf%equation_base%time_from_last_write, time_from_last_stat => self%equation_base%time_from_last_sta&
    &t, time_from_last_slice => self%equation_base%time_from_last_slice, time_from_last_probe => self%equ&
    &ation_base%time_from_last_probe, icyc0 => self%equation_base%icyc0, icyc => self%equation_base%icyc,&
    & time => self%equation_base%time, w_aux_cpu => self%w_aux_cpu, ijk_probe_cpu => self%ijk_probe_cpu,&
    & w_aux_probe_cpu => self%w_aux_probe_cpu, probe_coeff_cpu => self%probe_coeff_cpu,&
    & w_aux_probe => self%equation_base%w_aux_probe, nx => self%nx, ny => self%ny, nz => self%nz,&
    & ng => self%ng, nv_aux => self%nv_aux)


      time_from_last_write = time_from_last_write + self%equation_base%dt
      if (time_from_last_write >= self%equation_base%dtsave) then
        if (self%masterproc) write(*,*) 'time_from_last_write=',time_from_last_write
        if (self%masterproc) write(*,*) 'istore =', self%equation_base%istore
        call self%base_cpu%copy_to_field()
        self%equation_base%w_aux = self%w_aux_cpu
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
        call self%base_cpu%copy_to_field()
        self%equation_base%w_aux(:,:,:,6) = self%w_aux_cpu(:,:,:,6)
        self%equation_base%w_aux(:,:,:,7) = self%w_aux_cpu(:,:,:,7)
        call self%equation_base%compute_stats()
        if (self%equation_base%enable_stat_3d>0) call self%equation_base%compute_stats_3d()
        time_from_last_stat = 0._rkind
        self%equation_base%itav = self%equation_base%itav + 1
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
            sliceyz_aux(i,:,:,1:6) = self%w_aux_cpu(ii,:,:,1:6)
          enddo
          write(133,asynchronous="no") icyc,time
          write(133,asynchronous="no") self%grid%nxmax,self%grid%nymax,self%grid%nzmax
          write(133,asynchronous="no") isize,(self%equation_base%islice(i),i=1,isize)
          write(133,asynchronous="no") isize, self%field%ny, self%field%nz, 6
          write(133,asynchronous="yes") sliceyz_aux !(1:isize,1:self%field%ny,1:self%field%nz,1:6)
        endif
        if (allocated(self%equation_base%jslice)) then
          jsize = size(self%equation_base%jslice)
          do j=1,size(self%equation_base%jslice)
            jj = self%equation_base%jslice(j)
            slicexz_aux(:,j,:,1:6) = self%w_aux_cpu(:,jj,:,1:6)
          enddo
          if (self%equation_base%jslice(1)==1) slicexz_aux(:,1,:,2:4) = self%wallprop_cpu(:,:,2:4)
          write(134,asynchronous="no") icyc,time
          write(134,asynchronous="no") self%grid%nxmax,self%grid%nymax,self%grid%nzmax
          write(134,asynchronous="no") jsize,(self%equation_base%jslice(j),j=1,jsize)
          write(134,asynchronous="no") self%field%nx, jsize, self%field%nz, 6
          write(134,asynchronous="yes") slicexz_aux !(1:self%field%nx,1:jsize,1:self%field%nz,1:6)
        endif
        if (allocated(self%equation_base%kslice)) then
          ksize = size(self%equation_base%kslice)
          do k=1,size(self%equation_base%kslice)
            kk = self%equation_base%kslice(k)
            slicexy_aux(:,:,k,1:6) = self%w_aux_cpu(:,:,kk,1:6)
          enddo
          write(135,asynchronous="no") icyc,time
          write(135,asynchronous="no") self%grid%nxmax,self%grid%nymax,self%grid%nzmax
          write(135,asynchronous="no") ksize,(self%equation_base%kslice(k),k=1,ksize)
          write(135,asynchronous="no") self%field%nx, self%field%ny, ksize, 6
          write(135,asynchronous="yes") slicexy_aux !(1:self%field%nx,1:self%field%ny,1:ksize,1:6)
        endif
        time_from_last_slice = 0._rkind
      endif

      if (time_from_last_probe >= self%equation_base%dtprobe) then
        if (self%masterproc) write(*,*) 'time_from_last_probe=',time_from_last_probe
        if (self%equation_base%num_probe>0) then
          call probe_interpolation_subroutine(self%equation_base%num_probe,nx,ny,nz,ng,nv_aux,&
          &ijk_probe_cpu,w_aux_probe_cpu,w_aux_cpu,probe_coeff_cpu)
          w_aux_probe = w_aux_probe_cpu
          write(136) time, ((w_aux_probe(l,n),l=1,6),n=1,self%equation_base%num_probe)
        endif
        time_from_last_probe = 0._rkind
      endif

      time_from_last_rst = time_from_last_rst + self%equation_base%dt
      if (time_from_last_rst >= self%equation_base%dtsave_restart) then
        if (self%masterproc) write(*,*) 'time_from_last_rst=',time_from_last_rst
        call self%base_cpu%copy_to_field()

        if (self%equation_base%io_type_w==1) then
          call self%field%write_field_serial()
          call self%equation_base%write_stats_serial()
          if (self%equation_base%enable_stat_3d>0) call self%equation_base%write_stats_3d_serial()
        endif
        if (self%equation_base%io_type_w==2) then
          call self%field%write_field()
          call self%equation_base%write_stats()
          if (self%equation_base%enable_stat_3d>0) call self%equation_base%write_stats_3d()
        endif
        call self%equation_base%write_field_info()
        time_from_last_rst = 0._rkind
      endif

    endassociate
  end subroutine manage_output

endmodule streams_equation_singleideal_cpu_object


