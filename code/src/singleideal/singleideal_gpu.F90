module streams_equation_singleideal_gpu_object
! < STREAmS, Navier-Stokes equations, GPU backend.
!
  use streams_base_gpu_object
  use streams_field_object
  use streams_grid_object
  use streams_kernels_gpu
  use streams_parameters
  use crandom_f_mod
  use streams_equation_singleideal_object
  use MPI
  use CUDAFOR
  use ISO_C_BINDING
  use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
!
  implicit none
  private
  public :: equation_singleideal_gpu_object
!
  integer(kind=cuda_stream_kind) :: stream1
!
  integer, parameter :: EULERCENTRAL_THREADS_X=128,EULERCENTRAL_THREADS_Y=3
  integer, parameter :: EULERWENO_THREADS_X=128   ,EULERWENO_THREADS_Y=2
!
  type :: equation_singleideal_gpu_object
!   <
!   < w(1): rho
!   < w(2): rho * u
!   < w(3): rho * v
!   < w(4): rho * w
!   < w(5): rho * E
!   <```
!   < w_aux(1) : rho
!   < w_aux(2) : u
!   < w_aux(3) : v
!   < w_aux(4) : w
!   < w_aux(5) : h
!   < w_aux(6) : T
!   < w_aux(7) : viscosity
!   < w_aux(8) : ducros
!   < w_aux(9) : |omega|
!   < w_aux(10): div
!   <```
    type(base_gpu_object)       :: base_gpu               !< The base GPU handler.
    type(equation_singleideal_object)    :: equation_base  !< The equation base.
    type(field_object), pointer :: field=>null()          !< The field.
    type(grid_object),  pointer :: grid=>null()           !< The grid.
    integer(ikind)              :: ng                     !< Number of ghost cells.
    integer(ikind)              :: nx                     !< Number of cells in i direction.
    integer(ikind)              :: ny                     !< Number of cells in j direction.
    integer(ikind)              :: nz                     !< Number of cells in k direction.
    integer(ikind)              :: nv                     !< Number of variables.
    integer(ikind)              :: nv_aux                 !< Number of auxiliary variables.
    integer(ikind)              :: nprocs                 !< Number of auxiliary variables.
    integer(ikind)              :: myrank                 !< Number of auxiliary variables.
    integer(ikind)              :: error                  !< Number of auxiliary variables.
    real(rkind)                 :: time0
    real(rkind), allocatable, dimension(:,:), device :: coeff_deriv1_gpu
    real(rkind), allocatable, dimension(:,:), device :: coeff_deriv2_gpu
!
    integer :: icyc0, num_iter
    integer :: visc_model
    real(rkind) :: mu0, sutherland_S , T_ref_dim, powerlaw_vtexp
    logical :: masterproc
    integer :: mpi_err
    real(rkind), allocatable, dimension(:), device :: winf_gpu, winf_past_shock_gpu
    real(rkind), allocatable, dimension(:,:,:,:), device :: w_aux_gpu
    real(rkind), allocatable, dimension(:,:,:,:), device :: fhat_gpu
    real(rkind), allocatable, dimension(:,:,:,:), device :: w_aux_trans_gpu, fl_trans_gpu, fhat_trans_gpu
    real(rkind), allocatable, dimension(:,:,:,:), device :: fl_gpu, fln_gpu
    real(rkind), allocatable, dimension(:,:), device :: dcoe_gpu
    real(rkind), allocatable, dimension(:,:,:,:) :: w_var, w_var_t
    real(rkind), allocatable, dimension(:,:,:,:), device :: gplus_x_gpu, gminus_x_gpu
    real(rkind), allocatable, dimension(:,:,:,:), device :: gplus_y_gpu, gminus_y_gpu
    real(rkind), allocatable, dimension(:,:,:,:), device :: gplus_z_gpu, gminus_z_gpu
    integer, allocatable, dimension(:,:,:), device :: fluid_mask_gpu
    integer, allocatable, dimension(:,:,:,:), device :: ep_ord_change_gpu
    integer, allocatable, dimension(:,:,:), device :: ep_ord_change_x_gpu
!
    real(rkind), allocatable, dimension(:,:,:,:), device :: wrecyc_gpu
    real(rkind), allocatable, dimension(:,:,:), device :: wrecycav_gpu
    real(rkind), dimension(:,:,:), allocatable, device   :: wmean_gpu
!
    real(rkind), dimension(:,:,:), allocatable, device :: inflow_random_plane_gpu
    real(rkind), dimension(:), allocatable, device :: weta_inflow_gpu
    real(rkind), dimension(:), allocatable, device :: yplus_inflow_gpu, eta_inflow_gpu
    real(rkind), dimension(:), allocatable, device :: yplus_recyc_gpu, eta_recyc_gpu
    integer, dimension(:), allocatable, device :: map_j_inn_gpu, map_j_out_gpu
!
    real(rkind), dimension(:), allocatable, device :: cv_coeff_gpu, cp_coeff_gpu
!
    real(rkind), dimension(:,:), allocatable, device :: w_aux_probe_gpu
    real(rkind), dimension(:,:,:,:), allocatable, device :: probe_coeff_gpu
    integer, dimension(:,:), allocatable, device :: ijk_probe_gpu
    real(rkind), allocatable, dimension(:,:,:), device :: wallprop_gpu
!
!
!
!
!
!
!   
!
  contains
!   public methods
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
!
!
!
  endtype equation_singleideal_gpu_object
!
contains
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Utilities
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
  subroutine point_to_field(self, field)
    class(equation_singleideal_gpu_object), intent(inout) :: self              !< The equation.
    class(field_object), target :: field              !< The equation.
    self%field => field
  endsubroutine point_to_field
!
  subroutine point_to_grid(self, grid)
    class(equation_singleideal_gpu_object), intent(inout) :: self              !< The equation.
    class(grid_object), target :: grid              !< The equation.
    self%grid => grid
  endsubroutine point_to_grid
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
!
!
  subroutine rk_sync_old(self)
    class(equation_singleideal_gpu_object), intent(inout) :: self              !< The equation.
    integer :: istep, iercuda
    real(rkind) :: rhodt, gamdt, alpdt
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, &
      nv => self%nv,  nv_aux => self%nv_aux, &
      dt => self%equation_base%dt, &
      eul_imin => self%equation_base%eul_imin, eul_imax => self%equation_base%eul_imax, &
      eul_jmin => self%equation_base%eul_jmin, eul_jmax => self%equation_base%eul_jmax, &
      eul_kmin => self%equation_base%eul_kmin, eul_kmax => self%equation_base%eul_kmax, &
      channel_case => self%equation_base%channel_case)
!     
!     
      if (channel_case) self%equation_base%dpdx = 0._rkind
!     
      do istep=1,self%equation_base%nrk
        rhodt = self%equation_base%rhork(istep)*dt
        gamdt = self%equation_base%gamrk(istep)*dt
        alpdt = self%equation_base%alprk(istep)*dt
!
        call init_flux_cuf(nx, ny, nz, nv, self%fl_gpu, self%fln_gpu, rhodt)
        if (self%equation_base%conservative_viscous==1) then
          call self%visflx(mode=2)
        else
          call self%visflx(mode=0)
        endif
        call bcextr_var_cuf(nx, ny, nz, ng, self%w_aux_gpu(:,:,:,10:10))
        call self%base_gpu%bcswap_var(self%w_aux_gpu(:,:,:,10:10)) ! div/3
        call self%base_gpu%bcswap_var(self%w_aux_gpu(:,:,:,8:8 )) ! ducros
        if (self%equation_base%conservative_viscous==1) then
          call self%visflx(mode=7)
        else
          call self%visflx(mode=1)
        endif
        !@cuf iercuda=cudaDeviceSynchronize()
        call self%euler_x(eul_imin, eul_imax, 1-ng, nx+ng)
        !@cuf iercuda=cudaDeviceSynchronize()
        call self%euler_y(eul_jmin,eul_jmax)
        !@cuf iercuda=cudaDeviceSynchronize()
        call self%euler_z(eul_kmin,eul_kmax)
        if (self%equation_base%conservative_viscous==1) then
          call self%visflx(mode=5)
        endif
        call self%bc_nr()
        call update_flux_cuf(nx, ny, nz, nv, self%fl_gpu, self%fln_gpu, gamdt)
        if (channel_case) call self%force_rhs()
        call update_field_cuf(nx, ny, nz, ng, nv, self%base_gpu%w_gpu, self%fln_gpu, self%fluid_mask_gpu)
        if (channel_case) call self%force_var()
        call self%update_ghost()
        call self%compute_aux()
        !@cuf iercuda=cudaDeviceSynchronize()
      enddo
    endassociate
  endsubroutine rk_sync_old
!
  subroutine visflx(self, mode)
    class(equation_singleideal_gpu_object), intent(inout) :: self              !< The equation.
    integer, intent(in) :: mode
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, &
      nv => self%nv,  nv_aux => self%nv_aux, &
      dt => self%equation_base%dt, &
      visc_order => self%equation_base%visc_order, Prandtl => self%equation_base%Prandtl, &
      visc_model => self%visc_model, mu0 => self%mu0, &
      u0 => self%equation_base%u0, l0 => self%equation_base%l0, &
      T_ref_dim => self%T_ref_dim, &
      sutherland_S => self%sutherland_S, &
      powerlaw_vtexp => self%powerlaw_vtexp, &
      coeff_deriv1_gpu => self%coeff_deriv1_gpu, coeff_deriv2_gpu => self%coeff_deriv2_gpu, &
      fhat_trans_gpu => self%fhat_trans_gpu, fl_trans_gpu => self%fl_trans_gpu, fl_gpu => self%fl_gpu, &
      w_aux_gpu  => self%w_aux_gpu, w_aux_trans_gpu => self%w_aux_trans_gpu, &
      dcsidx_gpu => self%base_gpu%dcsidx_gpu,   &
      detady_gpu => self%base_gpu%detady_gpu,   &
      dzitdz_gpu => self%base_gpu%dzitdz_gpu,   &
      dcsidxs_gpu => self%base_gpu%dcsidxs_gpu, &
      detadys_gpu => self%base_gpu%detadys_gpu, &
      dzitdzs_gpu => self%base_gpu%dzitdzs_gpu, &
      dcsidx2_gpu => self%base_gpu%dcsidx2_gpu, &
      detady2_gpu => self%base_gpu%detady2_gpu, &
      dzitdz2_gpu => self%base_gpu%dzitdz2_gpu, &
      x_gpu => self%base_gpu%x_gpu, &
      y_gpu => self%base_gpu%y_gpu, &
      z_gpu => self%base_gpu%z_gpu, &
      eul_imin => self%equation_base%eul_imin, eul_imax => self%equation_base%eul_imax, &
      eul_jmin => self%equation_base%eul_jmin, eul_jmax => self%equation_base%eul_jmax, &
      eul_kmin => self%equation_base%eul_kmin, eul_kmax => self%equation_base%eul_kmax, &
      cv_coeff_gpu => self%cv_coeff_gpu, &
      cp_coeff_gpu => self%cp_coeff_gpu, &
      indx_cp_l => self%equation_base%indx_cp_l, &
      indx_cp_r => self%equation_base%indx_cp_r, &
      calorically_perfect => self%equation_base%calorically_perfect, &
      t0 => self%equation_base%t0, &
      channel_case => self%equation_base%channel_case)
      if (mode == 0) then ! laplacian
        call visflx_cuf(nx, ny, nz, ng, visc_order, &
          Prandtl, t0, indx_cp_l, indx_cp_r, cp_coeff_gpu, calorically_perfect, &
          u0, l0, self%base_gpu%w_gpu, self%w_aux_gpu, self%fl_gpu, &
          coeff_deriv1_gpu, coeff_deriv2_gpu, &
          dcsidx_gpu, detady_gpu, dzitdz_gpu,  &
          dcsidxs_gpu, detadys_gpu, dzitdzs_gpu,  &
          dcsidx2_gpu, detady2_gpu, dzitdz2_gpu, self%wallprop_gpu)
      elseif (mode == 1) then ! div
        call visflx_div_cuf(nx, ny, nz, ng, visc_order, &
          self%w_aux_gpu, self%fl_gpu, coeff_deriv1_gpu, &
          dcsidx_gpu, detady_gpu, dzitdz_gpu, stream1)
      elseif (mode == 2) then ! reduced
!       call visflx_reduced_cuf(nx, ny, nz, ng, visc_order, &
!         Prandtl, t0, indx_cp_l, indx_cp_r, cp_coeff_gpu, calorically_perfect, &
!         u0, l0, self%base_gpu%w_gpu, self%w_aux_gpu, self%fl_gpu, &
!         coeff_deriv1_gpu, coeff_deriv2_gpu, &
!         dcsidx_gpu, detady_gpu, dzitdz_gpu,  &
!         dcsidxs_gpu, detadys_gpu, dzitdzs_gpu,  &
!         dcsidx2_gpu, detady2_gpu, dzitdz2_gpu, self%wallprop_gpu, 1)
        call visflx_reduced_ord2_cuf(nx, ny, nz, ng, &
          Prandtl, t0, indx_cp_l, indx_cp_r, cp_coeff_gpu, calorically_perfect, &
          u0, l0, self%base_gpu%w_gpu, self%w_aux_gpu, self%fl_gpu, &
          x_gpu, y_gpu, z_gpu, self%wallprop_gpu, 1)
      elseif (mode == 3) then ! only sensor
        call sensor_cuf(nx, ny, nz, ng, u0, l0, self%w_aux_gpu)
      elseif (mode == 4) then ! laplacian
        call visflx_nosensor_cuf(nx, ny, nz, ng, visc_order, &
          Prandtl, t0, indx_cp_l, indx_cp_r, cp_coeff_gpu, calorically_perfect, &
          u0, l0, self%base_gpu%w_gpu, self%w_aux_gpu, self%fl_gpu, &
          coeff_deriv1_gpu, coeff_deriv2_gpu, &
          dcsidx_gpu, detady_gpu, dzitdz_gpu,  &
          dcsidxs_gpu, detadys_gpu, dzitdzs_gpu,  &
          dcsidx2_gpu, detady2_gpu, dzitdz2_gpu, self%wallprop_gpu)
      elseif (mode == 5) then ! staggered
        call visflx_x_cuf(nx, ny, nz, nv, ng, &
          Prandtl, t0, indx_cp_l, indx_cp_r, cp_coeff_gpu, calorically_perfect, &
          self%base_gpu%x_gpu, w_aux_trans_gpu, fl_trans_gpu, self%fl_gpu)
        call visflx_y_cuf(nx, ny, nz, nv, ng, &
          Prandtl, t0, indx_cp_l, indx_cp_r, cp_coeff_gpu, calorically_perfect, &
          self%base_gpu%y_gpu, self%w_aux_gpu, self%fl_gpu)
        call visflx_z_cuf(nx, ny, nz, nv, ng, &
          Prandtl, t0, indx_cp_l, indx_cp_r, cp_coeff_gpu, calorically_perfect, &
          self%base_gpu%z_gpu, self%w_aux_gpu, self%fl_gpu)
      elseif (mode == 6) then ! reduce_nosensor
!       call visflx_reduced_cuf(nx, ny, nz, ng, visc_order, &
!         Prandtl, t0, indx_cp_l, indx_cp_r, cp_coeff_gpu, calorically_perfect, &
!         u0, l0, self%base_gpu%w_gpu, self%w_aux_gpu, self%fl_gpu, &
!         coeff_deriv1_gpu, coeff_deriv2_gpu, &
!         dcsidx_gpu, detady_gpu, dzitdz_gpu,  &
!         dcsidxs_gpu, detadys_gpu, dzitdzs_gpu,  &
!         dcsidx2_gpu, detady2_gpu, dzitdz2_gpu, self%wallprop_gpu, 0)
        call visflx_reduced_ord2_cuf(nx, ny, nz, ng, &
          Prandtl, t0, indx_cp_l, indx_cp_r, cp_coeff_gpu, calorically_perfect, &
          u0, l0, self%base_gpu%w_gpu, self%w_aux_gpu, self%fl_gpu, &
          x_gpu, y_gpu, z_gpu, self%wallprop_gpu, 0)
      elseif (mode == 7) then ! div_ord2
        call visflx_div_ord2_cuf(nx, ny, nz, ng, &
          self%w_aux_gpu, self%fl_gpu, &
          x_gpu, y_gpu, z_gpu, stream1)
      endif
    endassociate
  endsubroutine visflx
!
  subroutine compute_aux(self, central, ghost)
    class(equation_singleideal_gpu_object), intent(inout) :: self              !< The equation.
    integer, intent(in), optional :: central, ghost
    integer :: central_, ghost_
    central_ = 1 ; if(present(central)) central_ = central
    ghost_ = 1   ; if(present(ghost))   ghost_ = ghost
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, &
      visc_model => self%visc_model, mu0 => self%mu0, &
      t0 => self%equation_base%t0, T_ref_dim => self%T_ref_dim, &
      sutherland_S => self%sutherland_S, &
      powerlaw_vtexp => self%powerlaw_vtexp, &
      cv_coeff_gpu => self%cv_coeff_gpu, &
      cp_coeff_gpu => self%cp_coeff_gpu, &
      indx_cp_l => self%equation_base%indx_cp_l, &
      indx_cp_r => self%equation_base%indx_cp_r, &
      calorically_perfect => self%equation_base%calorically_perfect, &
      rgas0 => self%equation_base%rgas0)
      if(central_ == 1 .and. ghost_ == 1) then
        call eval_aux_cuf(nx, ny, nz, ng, 1-ng, nx+ng, 1-ng, ny+ng, 1-ng, nz+ng, self%base_gpu%w_gpu, self%w_aux_gpu, &
          visc_model, mu0, t0, sutherland_S, T_ref_dim, &
          powerlaw_vtexp, VISC_POWER, VISC_SUTHERLAND, VISC_NO, cv_coeff_gpu,  &
          indx_cp_l, indx_cp_r, rgas0, calorically_perfect, tol_iter_nr,stream1)
      elseif(central_ == 1 .and. ghost_ == 0) then
        call eval_aux_cuf(nx, ny, nz, ng, 1, nx, 1, ny, 1, nz, self%base_gpu%w_gpu, self%w_aux_gpu, &
          visc_model, mu0, t0, sutherland_S, T_ref_dim, &
          powerlaw_vtexp, VISC_POWER, VISC_SUTHERLAND, VISC_NO, cv_coeff_gpu,  &
          indx_cp_l, indx_cp_r, rgas0, calorically_perfect, tol_iter_nr,stream1)
      elseif(central_ == 0 .and. ghost_ == 1) then
        call eval_aux_cuf(nx, ny, nz, ng, 1-ng, 0, 1-ng, ny+ng, 1-ng, nz+ng, self%base_gpu%w_gpu, self%w_aux_gpu, &
          visc_model, mu0, t0, sutherland_S, T_ref_dim, &
          powerlaw_vtexp, VISC_POWER, VISC_SUTHERLAND, VISC_NO, cv_coeff_gpu,  &
          indx_cp_l, indx_cp_r, rgas0, calorically_perfect, tol_iter_nr,stream1)
        call eval_aux_cuf(nx, ny, nz, ng, nx+1, nx+ng, 1-ng, ny+ng, 1-ng, nz+ng, self%base_gpu%w_gpu, self%w_aux_gpu, &
          visc_model, mu0, t0, sutherland_S, T_ref_dim, &
          powerlaw_vtexp, VISC_POWER, VISC_SUTHERLAND, VISC_NO, cv_coeff_gpu,  &
          indx_cp_l, indx_cp_r, rgas0, calorically_perfect, tol_iter_nr,stream1)
        call eval_aux_cuf(nx, ny, nz, ng, 1-ng, nx+ng, 1-ng, 0, 1-ng, nz+ng, self%base_gpu%w_gpu, self%w_aux_gpu, &
          visc_model, mu0, t0, sutherland_S, T_ref_dim, &
          powerlaw_vtexp, VISC_POWER, VISC_SUTHERLAND, VISC_NO, cv_coeff_gpu,  &
          indx_cp_l, indx_cp_r, rgas0, calorically_perfect, tol_iter_nr,stream1)
        call eval_aux_cuf(nx, ny, nz, ng, 1-ng, nx+ng, ny+1, ny+ng, 1-ng, nz+ng, self%base_gpu%w_gpu, self%w_aux_gpu, &
          visc_model, mu0, t0, sutherland_S, T_ref_dim, &
          powerlaw_vtexp, VISC_POWER, VISC_SUTHERLAND, VISC_NO, cv_coeff_gpu,  &
          indx_cp_l, indx_cp_r, rgas0, calorically_perfect, tol_iter_nr,stream1)
        call eval_aux_cuf(nx, ny, nz, ng, 1-ng, nx+ng, 1-ng, ny+ng, 1-ng, 0, self%base_gpu%w_gpu, self%w_aux_gpu, &
          visc_model, mu0, t0, sutherland_S, T_ref_dim, &
          powerlaw_vtexp, VISC_POWER, VISC_SUTHERLAND, VISC_NO, cv_coeff_gpu,  &
          indx_cp_l, indx_cp_r, rgas0, calorically_perfect, tol_iter_nr,stream1)
        call eval_aux_cuf(nx, ny, nz, ng, 1-ng, nx+ng, 1-ng, ny+ng, nz+1, nz+ng, self%base_gpu%w_gpu, self%w_aux_gpu, &
          visc_model, mu0, t0, sutherland_S, T_ref_dim, &
          powerlaw_vtexp, VISC_POWER, VISC_SUTHERLAND, VISC_NO, cv_coeff_gpu,  &
          indx_cp_l, indx_cp_r, rgas0, calorically_perfect, tol_iter_nr,stream1)
      endif
    endassociate
  endsubroutine compute_aux
!
  subroutine rk_sync(self)
    class(equation_singleideal_gpu_object), intent(inout) :: self              !< The equation.
    integer :: istep, lmax, iercuda
    real(rkind) :: rhodt, gamdt, alpdt
!   
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv, &
      dt => self%equation_base%dt, ep_order => self%equation_base%ep_order, &
      weno_scheme => self%equation_base%weno_scheme, &
      conservative_viscous => self%equation_base%conservative_viscous, &
      eul_imin => self%equation_base%eul_imin, eul_imax => self%equation_base%eul_imax, &
      eul_jmin => self%equation_base%eul_jmin, eul_jmax => self%equation_base%eul_jmax, &
      eul_kmin => self%equation_base%eul_kmin, eul_kmax => self%equation_base%eul_kmax, &
      channel_case => self%equation_base%channel_case, &
      nv_aux => self%nv_aux)
!
!     
!     
      if (channel_case) self%equation_base%dpdx = 0._rkind
!
      do istep=1,self%equation_base%nrk
        rhodt = self%equation_base%rhork(istep)*dt
        gamdt = self%equation_base%gamrk(istep)*dt
        alpdt = self%equation_base%alprk(istep)*dt
!
        call init_flux_cuf(nx, ny, nz, nv, self%fl_gpu, self%fln_gpu, rhodt)
        call self%base_gpu%bcswap()
        call self%compute_aux()
        !@cuf iercuda=cudaDeviceSynchronize()
        call self%euler_x(eul_imin, eul_imax, 1-ng, nx+ng)
        !@cuf iercuda=cudaDeviceSynchronize()
        if (conservative_viscous== 1) then
          call self%visflx(mode=6) ! 0=lapl, 1=div, 2=reduced, 3=sensor, 4=lapl_nosensor, 5=stag, 6=reduce_nosens
          call self%visflx(mode=5) ! 0=lapl, 1=div, 2=reduced, 3=sensor, 4=lapl_nosensor, 5=stag, 6=reduce_nosens
        else
          call self%visflx(mode=4) ! 0=lapl, 1=div, 2=reduced, 3=sensor, 4=lapl_nosensor, 5=stag, 6=reduce_nosens
        endif
        call bcextr_var_cuf(nx, ny, nz, ng, self%w_aux_gpu(:,:,:,10:10))
        call self%base_gpu%bcswap_var(self%w_aux_gpu(:,:,:,10:10)) ! div/3
        call self%euler_y(eul_jmin,eul_jmax)
        !@cuf iercuda=cudaDeviceSynchronize()
        call self%euler_z(eul_kmin,eul_kmax)
        if(istep == 3) then
          call self%visflx(mode=3) ! 0=lapl, 1=div, 2=reduced, 3=sensor, 4=lapl_nosensor, 5=stag, 6=reduce_nosens
          call self%base_gpu%bcswap_var(self%w_aux_gpu(:,:,:,8:8)) ! ducros
        endif
        if (conservative_viscous==1) then
          call self%visflx(mode=7)  ! 0=lapl, 1=div, 2=reduced, 3=sensor, 4=lapl_nosensor, 5=stag, 6=reduce_nosens
        else
          call self%visflx(mode=1)
        endif
        !@cuf iercuda=cudaDeviceSynchronize()
        call self%bc_nr()
        call update_flux_cuf(nx, ny, nz, nv, self%fl_gpu, self%fln_gpu, gamdt)
        if (channel_case) call self%force_rhs()
        call update_field_cuf(nx, ny, nz, ng, nv, self%base_gpu%w_gpu, self%fln_gpu, self%fluid_mask_gpu)
        if (channel_case) call self%force_var()
        call self%update_ghost(do_swap=0)
      enddo
!
    endassociate
  endsubroutine rk_sync
!
  subroutine rk_async(self)
    class(equation_singleideal_gpu_object), intent(inout) :: self              !< The equation.
    integer :: istep, lmax, iercuda
    real(rkind) :: rhodt, gamdt, alpdt
!   
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv, &
      dt => self%equation_base%dt, ep_order => self%equation_base%ep_order, &
      weno_scheme => self%equation_base%weno_scheme, &
      conservative_viscous => self%equation_base%conservative_viscous, &
      eul_imin => self%equation_base%eul_imin, eul_imax => self%equation_base%eul_imax, &
      eul_jmin => self%equation_base%eul_jmin, eul_jmax => self%equation_base%eul_jmax, &
      eul_kmin => self%equation_base%eul_kmin, eul_kmax => self%equation_base%eul_kmax, &
      channel_case => self%equation_base%channel_case, &
      nv_aux => self%nv_aux)
!
!
      if (channel_case) self%equation_base%dpdx = 0._rkind
      lmax = max(ep_order/2, weno_scheme) ! max stencil width
!
      do istep=1,self%equation_base%nrk
        rhodt = self%equation_base%rhork(istep)*dt
        gamdt = self%equation_base%gamrk(istep)*dt
        alpdt = self%equation_base%alprk(istep)*dt
!
        call init_flux_cuf(nx, ny, nz, nv, self%fl_gpu, self%fln_gpu, rhodt)
!
!
        call self%compute_aux(central=1, ghost=0)
        !@cuf iercuda=cudaDeviceSynchronize()
        call self%base_gpu%bcswap(steps=[.true.,.false.,.false.])
        call self%euler_x(lmax+1,nx-lmax,1,nx)
        call self%base_gpu%bcswap(steps=[.false.,.true.,.true.])
        call self%compute_aux(central=0, ghost=1)
        if (lmax-eul_imin >= 0)      call self%euler_x(eul_imin,lmax,1-ng,0)
        if (eul_imax-nx+lmax-1 >= 0) call self%euler_x(nx-lmax+1,eul_imax,nx+1,nx+ng)
        if (conservative_viscous == 1) then
          call self%visflx(mode=6) ! 0=lapl, 1=div, 2=reduced, 3=sensor, 4=lapl_nosensor, 5=stag, 6=reduce_nosens
          call self%visflx(mode=5) ! 0=lapl, 1=div, 2=reduced, 3=sensor, 4=lapl_nosensor, 5=stag, 6=reduce_nosens
        else
          call self%visflx(mode=4) ! 0=lapl, 1=div, 2=reduced, 3=sensor, 4=lapl_nosensor, 5=stag, 6=reduce_nosens
        endif
        call bcextr_var_cuf(nx, ny, nz, ng, self%w_aux_gpu(:,:,:,10:10))
        call self%base_gpu%bcswap_var(self%w_aux_gpu(:,:,:,10:10), steps=[.true.,.false.,.false.]) ! div/3
        call self%euler_y(eul_jmin,eul_jmax)
        call self%base_gpu%bcswap_var(self%w_aux_gpu(:,:,:,10:10), steps=[.false.,.true.,.true.]) ! div/3
        call self%euler_z(eul_kmin,eul_kmax)
        if(istep == 3) then
          call self%visflx(mode=3) ! 0=lapl, 1=div, 2=reduced, 3=sensor, 4=lapl_nosensor, 5=stag
          call self%base_gpu%bcswap_var(self%w_aux_gpu(:,:,:,8:8), steps=[.true.,.false.,.false.]) ! ducros
          if (conservative_viscous==1) then
            call self%visflx(mode=7)  ! 0=lapl, 1=div, 2=reduced, 3=sensor, 4=lapl_nosensor, 5=stag, 6=reduce_nosens
          else
            call self%visflx(mode=1)
          endif
          call self%base_gpu%bcswap_var(self%w_aux_gpu(:,:,:,8:8), steps=[.false.,.true.,.true.]) ! ducros
        else
          if (conservative_viscous==1) then
            call self%visflx(mode=7)  ! 0=lapl, 1=div, 2=reduced, 3=sensor, 4=lapl_nosensor, 5=stag, 6=reduce_nosens
          else
            call self%visflx(mode=1)
          endif
          !@cuf iercuda=cudaDeviceSynchronize()
        endif
        call self%bc_nr()
        call update_flux_cuf(nx, ny, nz, nv, self%fl_gpu, self%fln_gpu, gamdt)
        if (channel_case) call self%force_rhs()
        call update_field_cuf(nx, ny, nz, ng, nv, self%base_gpu%w_gpu, self%fln_gpu, self%fluid_mask_gpu)
        if (channel_case) call self%force_var()
        call self%update_ghost(do_swap=0)
      enddo
    endassociate
  endsubroutine rk_async
!
  subroutine euler_x(self, istart, iend, istart_trans, iend_trans)
    class(equation_singleideal_gpu_object), intent(inout) :: self
    integer, intent(in) :: istart, iend, istart_trans, iend_trans
    integer :: lmax, weno_size
    type(dim3) :: grid, tBlock
    integer :: force_zero_flux_min,force_zero_flux_max
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, &
      nv => self%nv,  nv_aux => self%nv_aux, &
      ep_order => self%equation_base%ep_order, &
      force_zero_flux => self%equation_base%force_zero_flux, &
      coeff_deriv1_gpu => self%coeff_deriv1_gpu,  dcsidx_gpu => self%base_gpu%dcsidx_gpu, &
      fhat_gpu => self%fhat_gpu, w_gpu => self%base_gpu%w_gpu, w_aux_gpu => self%w_aux_gpu, &
      fl_gpu => self%fl_gpu, sensor_threshold => self%equation_base%sensor_threshold, &
      fl_trans_gpu => self%fl_trans_gpu, &
      weno_scheme => self%equation_base%weno_scheme, &
      weno_version => self%equation_base%weno_version, &
      gplus_x_gpu => self%gplus_x_gpu, gminus_x_gpu => self%gminus_x_gpu, &
      w_aux_trans_gpu => self%w_aux_trans_gpu, fhat_trans_gpu => self%fhat_trans_gpu, &
      cp_coeff_gpu => self%cp_coeff_gpu, &
      indx_cp_l => self%equation_base%indx_cp_l, indx_cp_r => self%equation_base%indx_cp_r, &
      calorically_perfect => self%equation_base%calorically_perfect, &
      ep_ord_change_x_gpu => self%ep_ord_change_x_gpu, nkeep => self%equation_base%nkeep, &
      flux_splitting => self%equation_base%flux_splitting, &
      rgas0 => self%equation_base%rgas0)
!
      call euler_x_transp_cuf(nx, ny, nz, ng, istart_trans, iend_trans, nv_aux, w_aux_gpu, w_aux_trans_gpu, stream1)
!
      weno_size  = 2*weno_scheme
      lmax = ep_order/2 ! max stencil width
      force_zero_flux_min = force_zero_flux(1)
      force_zero_flux_max = force_zero_flux(2)
      tBlock = dim3(EULERWENO_THREADS_X,EULERWENO_THREADS_Y,1)
      grid = dim3(ceiling(real(ny)/tBlock%x),ceiling(real(nz)/tBlock%y),1)
!
      if (flux_splitting==1) then
        call euler_x_fluxes_hybrid_rusanov_kernel<<<grid, tBlock, 0, stream1>>>(nv, nv_aux, nx, ny, nz, ng, &
          istart, iend, lmax, nkeep, rgas0, coeff_deriv1_gpu, dcsidx_gpu, w_aux_trans_gpu, fhat_trans_gpu, &
          force_zero_flux_min, force_zero_flux_max, weno_scheme, weno_version, &
          sensor_threshold, weno_size, gplus_x_gpu, gminus_x_gpu, cp_coeff_gpu, indx_cp_l, indx_cp_r, &
          ep_ord_change_x_gpu, calorically_perfect, tol_iter_nr,self%equation_base%rho0,self%equation_base%u0, &
          self%equation_base%t0)
      else
        call euler_x_fluxes_hybrid_kernel<<<grid, tBlock, 0, stream1>>>(nv, nv_aux, nx, ny, nz, ng, &
          istart, iend, lmax, nkeep, rgas0, coeff_deriv1_gpu, dcsidx_gpu, w_aux_trans_gpu, fhat_trans_gpu, &
          force_zero_flux_min, force_zero_flux_max, weno_scheme, weno_version, &
          sensor_threshold, weno_size, gplus_x_gpu, gminus_x_gpu, cp_coeff_gpu, indx_cp_l, indx_cp_r, &
          ep_ord_change_x_gpu, calorically_perfect, tol_iter_nr,self%equation_base%rho0,self%equation_base%u0, &
          self%equation_base%t0)
      endif
!
      call euler_x_update_cuf(nx, ny, nz, ng, nv, istart, iend, fhat_trans_gpu, fl_trans_gpu, fl_gpu, dcsidx_gpu, stream1)
    endassociate
  endsubroutine euler_x
!
  subroutine euler_y(self, eul_jmin, eul_jmax)
    class(equation_singleideal_gpu_object), intent(inout) :: self
    integer, intent(in) :: eul_jmin, eul_jmax
    integer :: lmax, weno_size
    type(dim3) :: grid, tBlock
    integer :: force_zero_flux_min,force_zero_flux_max
!
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, &
      nv => self%nv,  nv_aux => self%nv_aux, &
      ep_order => self%equation_base%ep_order, &
      force_zero_flux => self%equation_base%force_zero_flux, &
      coeff_deriv1_gpu => self%coeff_deriv1_gpu, detady_gpu => self%base_gpu%detady_gpu, &
      fhat_gpu => self%fhat_gpu, w_gpu => self%base_gpu%w_gpu, w_aux_gpu => self%w_aux_gpu, &
      fl_gpu => self%fl_gpu, sensor_threshold => self%equation_base%sensor_threshold, &
      weno_scheme => self%equation_base%weno_scheme, &
      weno_version => self%equation_base%weno_version, &
      gplus_y_gpu => self%gplus_y_gpu, gminus_y_gpu => self%gminus_y_gpu, &
      indx_cp_l => self%equation_base%indx_cp_l, indx_cp_r => self%equation_base%indx_cp_r, &
      calorically_perfect => self%equation_base%calorically_perfect, &
      cp_coeff_gpu => self%cp_coeff_gpu, &
      ep_ord_change_gpu => self%ep_ord_change_gpu, nkeep => self%equation_base%nkeep, &
      flux_splitting => self%equation_base%flux_splitting, &
      rgas0 => self%equation_base%rgas0)
      weno_size = 2*weno_scheme
      lmax  = ep_order/2 ! max stencil width
      force_zero_flux_min = force_zero_flux(3)
      force_zero_flux_max = force_zero_flux(4)
!
      tBlock = dim3(EULERWENO_THREADS_X,EULERWENO_THREADS_Y,1)
      grid = dim3(ceiling(real(nx)/tBlock%x),ceiling(real(nz)/tBlock%y),1)
      if (flux_splitting==1) then
        call euler_y_hybrid_rusanov_kernel<<<grid, tBlock, 0, stream1>>>(nv, nv_aux, nx, ny, nz, ng, &
          eul_jmin, eul_jmax, lmax, nkeep, rgas0, w_aux_gpu, fl_gpu, coeff_deriv1_gpu, detady_gpu, fhat_gpu, &
          force_zero_flux_min, force_zero_flux_max, weno_scheme, weno_version, &
          sensor_threshold, weno_size, w_gpu, gplus_y_gpu, gminus_y_gpu, cp_coeff_gpu, indx_cp_l, indx_cp_r, &
          ep_ord_change_gpu, calorically_perfect, tol_iter_nr,self%equation_base%rho0,self%equation_base%u0, &
          self%equation_base%t0)
      else
        call euler_y_hybrid_kernel<<<grid, tBlock, 0, stream1>>>(nv, nv_aux, nx, ny, nz, ng, &
          eul_jmin, eul_jmax, lmax, nkeep, rgas0, w_aux_gpu, fl_gpu, coeff_deriv1_gpu, detady_gpu, fhat_gpu, &
          force_zero_flux_min, force_zero_flux_max, weno_scheme, weno_version, &
          sensor_threshold, weno_size, w_gpu, gplus_y_gpu, gminus_y_gpu, cp_coeff_gpu, indx_cp_l, indx_cp_r, &
          ep_ord_change_gpu, calorically_perfect, tol_iter_nr,self%equation_base%rho0,self%equation_base%u0, &
          self%equation_base%t0)
      endif
    endassociate
  endsubroutine euler_y
!
  subroutine euler_z(self, eul_kmin, eul_kmax)
    class(equation_singleideal_gpu_object), intent(inout) :: self
    integer, intent(in) :: eul_kmin, eul_kmax
    integer :: lmax, weno_size
    type(dim3) :: grid, tBlock
    integer :: force_zero_flux_min, force_zero_flux_max
!
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, &
      nv => self%nv,  nv_aux => self%nv_aux, &
      ep_order => self%equation_base%ep_order, &
      force_zero_flux => self%equation_base%force_zero_flux, &
      coeff_deriv1_gpu => self%coeff_deriv1_gpu, dzitdz_gpu => self%base_gpu%dzitdz_gpu, &
      fhat_gpu => self%fhat_gpu, w_gpu => self%base_gpu%w_gpu, w_aux_gpu => self%w_aux_gpu, &
      fl_gpu => self%fl_gpu, sensor_threshold => self%equation_base%sensor_threshold, &
      weno_scheme => self%equation_base%weno_scheme, &
      weno_version => self%equation_base%weno_version, &
      gplus_z_gpu => self%gplus_z_gpu, gminus_z_gpu => self%gminus_z_gpu, &
      cp_coeff_gpu => self%cp_coeff_gpu, &
      indx_cp_l => self%equation_base%indx_cp_l, indx_cp_r => self%equation_base%indx_cp_r, &
      calorically_perfect => self%equation_base%calorically_perfect, &
      ep_ord_change_gpu => self%ep_ord_change_gpu, nkeep => self%equation_base%nkeep, &
      flux_splitting => self%equation_base%flux_splitting, &
      rgas0 => self%equation_base%rgas0, &
      inflow_random_plane => self%equation_base%inflow_random_plane, &
      inflow_random_plane_gpu => self%inflow_random_plane_gpu)
      weno_size = 2*weno_scheme
      lmax = ep_order/2 ! max stencil width
      force_zero_flux_min = force_zero_flux(5)
      force_zero_flux_max = force_zero_flux(6)
      tBlock = dim3(EULERWENO_THREADS_X,EULERWENO_THREADS_Y,1)
      grid = dim3(ceiling(real(nx)/tBlock%x),ceiling(real(ny)/tBlock%y),1)
      if (flux_splitting==1) then
        call euler_z_hybrid_rusanov_kernel<<<grid, tBlock, 0, 0>>>(nv, nv_aux, nx, ny, nz, ng, &
          eul_kmin, eul_kmax, lmax, nkeep, rgas0, w_aux_gpu, fl_gpu, coeff_deriv1_gpu, dzitdz_gpu, fhat_gpu, &
          force_zero_flux_min, force_zero_flux_max, weno_scheme, weno_version, &
          sensor_threshold, weno_size, w_gpu, gplus_z_gpu, gminus_z_gpu, cp_coeff_gpu, indx_cp_l, indx_cp_r, &
          ep_ord_change_gpu, calorically_perfect, tol_iter_nr,self%equation_base%rho0,self%equation_base%u0, &
          self%equation_base%t0)
      else
        call euler_z_hybrid_kernel<<<grid, tBlock, 0, 0>>>(nv, nv_aux, nx, ny, nz, ng, &
          eul_kmin, eul_kmax, lmax, nkeep, rgas0, w_aux_gpu, fl_gpu, coeff_deriv1_gpu, dzitdz_gpu, fhat_gpu, &
          force_zero_flux_min, force_zero_flux_max, weno_scheme, weno_version, &
          sensor_threshold, weno_size, w_gpu, gplus_z_gpu, gminus_z_gpu, cp_coeff_gpu, indx_cp_l, indx_cp_r, &
          ep_ord_change_gpu, calorically_perfect, tol_iter_nr,self%equation_base%rho0,self%equation_base%u0, &
          self%equation_base%t0)
      endif
      if (self%equation_base%recyc) call get_crandom_f(inflow_random_plane(2:self%equation_base%jbl_inflow,1:nz,1:3))
!
    endassociate
  endsubroutine euler_z
!
  subroutine force_rhs(self)
    class(equation_singleideal_gpu_object), intent(inout) :: self
    real(rkind), dimension(5) :: bulk5, bulk5g
!
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, &
      nxmax => self%equation_base%grid%nxmax, nzmax => self%equation_base%grid%nzmax, &
      yn => self%equation_base%field%yn, yn_gpu => self%base_gpu%yn_gpu, &
      w_gpu => self%base_gpu%w_gpu, w_aux_gpu => self%w_aux_gpu, fln_gpu => self%fln_gpu, &
      volchan => self%equation_base%volchan, fluid_mask_gpu => self%fluid_mask_gpu)
!
      call force_rhs_1_cuf(nx, ny, nz, ng, yn_gpu, fln_gpu, w_gpu, w_aux_gpu, bulk5, fluid_mask_gpu)
!
      call mpi_allreduce(bulk5,bulk5g,5,mpi_prec,mpi_sum,self%equation_base%field%mp_cart,self%mpi_err)
!
      bulk5g = bulk5g/volchan !/nxmax/nzmax/(yn(ny+1)-yn(1))
      self%equation_base%dpdx    = self%equation_base%dpdx+bulk5g(2)
      self%equation_base%rhobulk = bulk5g(3)
      self%equation_base%ubulk   = bulk5g(4)/self%equation_base%rhobulk
      self%equation_base%tbulk   = bulk5g(5)/self%equation_base%rhobulk/self%equation_base%ubulk
!
!     Add forcing terms in momentum and energy equation
      call force_rhs_2_cuf(nx, ny, nz, ng, fln_gpu, w_aux_gpu, bulk5g, fluid_mask_gpu)
!     
    endassociate
  endsubroutine force_rhs
! 
  subroutine force_var(self)
    class(equation_singleideal_gpu_object), intent(inout) :: self
    real(rkind) :: bulkt, bulktg, tbtarget, tbdiff
!
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, &
      nxmax => self%equation_base%grid%nxmax, nzmax => self%equation_base%grid%nzmax, &
      yn => self%equation_base%field%yn, yn_gpu => self%base_gpu%yn_gpu, &
      w_gpu => self%base_gpu%w_gpu, w_aux_gpu => self%w_aux_gpu, fln_gpu => self%fln_gpu, &
      volchan => self%equation_base%volchan, fluid_mask_gpu => self%fluid_mask_gpu, &
      cv_coeff_gpu => self%cv_coeff_gpu, &
      t0 => self%equation_base%t0, &
      calorically_perfect => self%equation_base%calorically_perfect, &
      indx_cp_l => self%equation_base%indx_cp_l, indx_cp_r => self%equation_base%indx_cp_r)
!
      if (self%equation_base%theta_wall>=-1._rkind) then
!       
        call force_var_1_cuf(nx, ny, nz, ng, yn_gpu, fln_gpu, w_gpu, w_aux_gpu, bulkt, fluid_mask_gpu, cv_coeff_gpu, &
          indx_cp_l, indx_cp_r, t0, calorically_perfect, tol_iter_nr)
!       
        call mpi_allreduce(bulkt,bulktg,1,mpi_prec,mpi_sum,self%equation_base%field%mp_cart,self%mpi_err)
!
        bulktg = bulktg/volchan !/nxmax/nzmax/(yn(ny+1)-yn(1))
        bulktg = bulktg/self%equation_base%rhobulk/self%equation_base%ubulk
        tbtarget = self%equation_base%T_bulk_target
        tbdiff   = tbtarget-bulktg
!       
        call force_var_2_cuf(nx, ny, nz, ng, w_gpu, w_aux_gpu, tbdiff, fluid_mask_gpu, cv_coeff_gpu, indx_cp_l, indx_cp_r, &
          t0, calorically_perfect)
!       
      endif
!     
    endassociate
  endsubroutine force_var
!
  subroutine recyc_exchange(self)
    class(equation_singleideal_gpu_object), intent(inout) :: self
!
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
!
    associate(nzmax => self%equation_base%grid%nzmax, ng => self%equation_base%grid%ng, &
      nx => self%equation_base%field%nx, ny => self%equation_base%field%ny, nz => self%equation_base%field%nz, &
      ncoords => self%equation_base%field%ncoords, mp_cart => self%equation_base%field%mp_cart, &
      nblocks => self%equation_base%field%nblocks, &
      iermpi => self%mpi_err, w_gpu => self%base_gpu%w_gpu, nv => self%equation_base%nv, &
      wbuf1s_gpu => self%base_gpu%wbuf1s_gpu, wbuf2r_gpu => self%base_gpu%wbuf2r_gpu, &
      wbuf1r_gpu => self%base_gpu%wbuf1r_gpu, wrecyc_gpu => self%wrecyc_gpu, &
      ibrecyc => self%equation_base%ib_recyc, irecyc => self%equation_base%i_recyc)
!
      kshiftglob    = nzmax/2 ! global shift in the spanwise direction (between 0 and nzmax-1)
      n1_start_send = 1
      n1_end_send   = nz-mod(kshiftglob,nz)
      n2_start_send = n1_end_send+1
      n2_end_send   = nz
      n1_start_recv = 1+mod(kshiftglob,nz)
      n1_end_recv   = nz
      n2_start_recv = 1
      n2_end_recv   = mod(kshiftglob,nz)
!
      req = mpi_request_null
!
      if (ncoords(1)==ibrecyc) then ! Send data
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
!       ! call mpi_ssend(wbuf1s_gpu,indx,mpi_prec,0,2000,mp_cartx,iermpi)
      endif
      if (ncoords(1)==0) then ! Receive data
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
!
  end subroutine recyc_exchange
!
  subroutine bc_nr(self)
    class(equation_singleideal_gpu_object), intent(inout) :: self
    integer :: ilat, dir, start_or_end
    type(dim3) :: grid, tBlock
!
    associate(bctags_nr => self%equation_base%bctags_nr, nx => self%nx, ny => self%ny, nz => self%nz, nv => self%nv, &
      ng => self%ng, w_aux_gpu => self%w_aux_gpu, w_gpu => self%base_gpu%w_gpu, &
      fl_gpu => self%fl_gpu, dcsidx_gpu => self%base_gpu%dcsidx_gpu,  detady_gpu => self%base_gpu%detady_gpu, &
      dzitdz_gpu => self%base_gpu%dzitdz_gpu, &
      indx_cp_l => self%equation_base%indx_cp_l, indx_cp_r => self%equation_base%indx_cp_r, &
      calorically_perfect => self%equation_base%calorically_perfect, rgas0 => self%equation_base%rgas0, &
      t0 => self%equation_base%t0, cp_coeff_gpu => self%cp_coeff_gpu, winf_gpu => self%winf_gpu)
!
      do ilat=1,6! loop on all sides of the boundary (3D -> 6)
        if(bctags_nr(ilat) > 0) then
          dir          = (ilat-1)/2   +1
          start_or_end = mod(ilat-1,2)+1
!         1 - NR
!         2 - relax
!         6 - reflective wall
          if(dir == 1) then
            tBlock = dim3(EULERWENO_THREADS_X,EULERWENO_THREADS_Y,1)
            grid = dim3(ceiling(real(ny)/tBlock%x),ceiling(real(nz)/tBlock%y),1)
            call bc_nr_lat_x_kernel<<<grid, tBlock, 0, 0>>>(start_or_end, bctags_nr(ilat), &
              nx, ny, nz, ng, nv, w_aux_gpu, w_gpu, fl_gpu, dcsidx_gpu, &
              indx_cp_l, indx_cp_r, cp_coeff_gpu, winf_gpu, calorically_perfect,rgas0,t0)
          endif
          if(dir == 2) then
            tBlock = dim3(EULERWENO_THREADS_X,EULERWENO_THREADS_Y,1)
            grid = dim3(ceiling(real(nx)/tBlock%x),ceiling(real(nz)/tBlock%y),1)
            call bc_nr_lat_y_kernel<<<grid, tBlock, 0, 0>>>(start_or_end, bctags_nr(ilat), &
              nx, ny, nz, ng, nv, w_aux_gpu, w_gpu, fl_gpu, detady_gpu, &
              indx_cp_l, indx_cp_r, cp_coeff_gpu, calorically_perfect,rgas0,t0)
          endif
          if(dir == 3) then
            tBlock = dim3(EULERWENO_THREADS_X,EULERWENO_THREADS_Y,1)
            grid = dim3(ceiling(real(nx)/tBlock%x),ceiling(real(ny)/tBlock%y),1)
            call bc_nr_lat_z_kernel<<<grid, tBlock, 0, 0>>>(start_or_end, bctags_nr(ilat), &
              nx, ny, nz, ng, nv, w_aux_gpu, w_gpu, fl_gpu, dzitdz_gpu, &
              indx_cp_l, indx_cp_r, cp_coeff_gpu, calorically_perfect,rgas0,t0)
          endif
        endif
      enddo
    endassociate
!
  endsubroutine bc_nr
!
  subroutine update_ghost(self, do_swap)
    class(equation_singleideal_gpu_object), intent(inout) :: self
    integer, intent(in), optional :: do_swap
    integer :: do_swap_
    integer :: ilat
!
    do_swap_ = 1 ; if (present(do_swap)) do_swap_ = do_swap
!
    if (self%equation_base%recyc) call self%recyc_exchange()
!
    do ilat=1,6! loop on all sides of the boundary (3D -> 6)
      select case(self%equation_base%bctags(ilat))
      case(0)
      case(1)
        call bcfree_cuf(ilat, self%nx, self%ny, self%nz, self%ng, self%nv, self%winf_gpu, self%base_gpu%w_gpu)
      case(2)
        call bcextr_cuf(ilat, self%nx, self%ny, self%nz, self%ng, self%nv, self%base_gpu%w_gpu)
      case(4)
        call bcextr_sub_cuf(ilat, self%nx, self%ny, self%nz, self%ng, &
          self%equation_base%p0, self%equation_base%rgas0, self%base_gpu%w_gpu, self%equation_base%indx_cp_l, &
          self%equation_base%indx_cp_r, self%cv_coeff_gpu, self%equation_base%t0, &
          self%equation_base%calorically_perfect)
      case(5)
        call bcsym_cuf(ilat, self%nx, self%ny, self%nz, self%ng, &
          self%equation_base%T_wall, self%base_gpu%w_gpu, self%w_aux_gpu, &
          self%equation_base%indx_cp_l, self%equation_base%indx_cp_r, self%cv_coeff_gpu, &
          self%equation_base%calorically_perfect)
      case(6)
        if (self%equation_base%channel_case) then
          call bcwall_staggered_cuf(ilat, self%nx, self%ny, self%nz, self%ng, &
            self%equation_base%T_wall, self%base_gpu%w_gpu, self%w_aux_gpu, &
            self%equation_base%indx_cp_l, self%equation_base%indx_cp_r, self%cv_coeff_gpu, self%equation_base%t0, &
            self%equation_base%rgas0, self%equation_base%calorically_perfect,tol_iter_nr)
        else
          call bcwall_cuf(ilat, self%nx, self%ny, self%nz, self%ng, &
            self%equation_base%T_wall, self%base_gpu%w_gpu, self%w_aux_gpu, &
            self%equation_base%indx_cp_l, self%equation_base%indx_cp_r, self%cv_coeff_gpu, self%equation_base%t0, &
            self%equation_base%rgas0, self%equation_base%calorically_perfect,tol_iter_nr)
        endif
      case(7)
        call bcshock_cuf(ilat, self%nx, self%ny, self%nz, self%ng, self%nv, self%base_gpu%w_gpu, &
          self%winf_gpu, self%winf_past_shock_gpu, self%equation_base%xshock_imp, &
          self%equation_base%shock_angle, self%base_gpu%x_gpu, self%base_gpu%y_gpu, &
          self%equation_base%tanhfacs)
      case(8)
!
      case(9)
        call bclam_cuf(ilat, self%nx, self%ny, self%nz, self%ng, self%nv, self%base_gpu%w_gpu, self%wmean_gpu, &
          self%equation_base%p0, self%equation_base%rgas0, self%equation_base%indx_cp_l, &
          self%equation_base%indx_cp_r, self%cv_coeff_gpu, self%equation_base%t0, &
          self%equation_base%calorically_perfect)
      case(10)
        call self%bcrecyc(ilat)
      endselect
    enddo
!
    if (do_swap_ == 1) call self%base_gpu%bcswap()
!
  endsubroutine update_ghost
!
  subroutine bcrecyc(self, ilat)
!   
    class(equation_singleideal_gpu_object), intent(inout) :: self              !< The equation.
!   Apply recycling-rescaling boundary condition
!   
    integer, intent(in) :: ilat
    integer :: ntot
!   
    if (ilat == 1) then
      associate(nx => self%nx, ny => self%ny, nz => self%nz, nzmax => self%grid%nzmax, ng => self%ng, nv => self%nv, &
        wrecycav_gpu => self%wrecycav_gpu, wrecyc_gpu => self%wrecyc_gpu, &
        mp_cartz => self%field%mp_cartz, iermpi => self%mpi_err, &
        p0 => self%equation_base%p0, w_gpu => self%base_gpu%w_gpu, &
        wmean_gpu => self%wmean_gpu, weta_inflow_gpu  => self%weta_inflow_gpu, &
        map_j_inn_gpu => self%map_j_inn_gpu, map_j_out_gpu => self%map_j_out_gpu, &
        yplus_inflow_gpu => self%yplus_inflow_gpu, eta_inflow_gpu => self%eta_inflow_gpu, &
        yplus_recyc_gpu => self%yplus_recyc_gpu, eta_recyc_gpu => self%eta_recyc_gpu, &
        betarecyc => self%equation_base%betarecyc, i_recyc => self%equation_base%i_recyc, &
        indx_cp_l => self%equation_base%indx_cp_l, &
        indx_cp_r => self%equation_base%indx_cp_r, &
        rgas0 => self%equation_base%rgas0, &
        calorically_perfect => self%equation_base%calorically_perfect, &
        t0 => self%equation_base%t0, &
        u0 => self%equation_base%u0, &
        l0 => self%equation_base%l0, &
        cv_coeff_gpu => self%cv_coeff_gpu, &
        inflow_random_plane => self%equation_base%inflow_random_plane, &
        inflow_random_plane_gpu => self%inflow_random_plane_gpu)
!       Compute spanwise averages at the recycling station
        call bcrecyc_cuf_1(nx, ny, nz, ng, nv, wrecycav_gpu, wrecyc_gpu)
!
!       call get_crandom_f(inflow_random_plane(2:ny,1:nz,1:3)) ! this is done in euler z to overlap
        inflow_random_plane_gpu = inflow_random_plane
!
        ntot = ng*ny*nv
        call mpi_allreduce(MPI_IN_PLACE,wrecycav_gpu,ntot,mpi_prec,mpi_sum,mp_cartz,iermpi)
!
!       Remove average
        call bcrecyc_cuf_2(nx, ny, nz, nzmax, ng, wrecycav_gpu, wrecyc_gpu)
!
!       Apply bc recycling
        call bcrecyc_cuf_3(nx, ny, nz, ng, p0, u0, rgas0, w_gpu, wmean_gpu, wrecyc_gpu, &
          weta_inflow_gpu, map_j_inn_gpu, map_j_out_gpu, &
          yplus_inflow_gpu, eta_inflow_gpu, yplus_recyc_gpu, eta_recyc_gpu, betarecyc, inflow_random_plane_gpu, &
          indx_cp_l, indx_cp_r, cv_coeff_gpu, t0, calorically_perfect)
      endassociate
    endif
!
  end subroutine bcrecyc
!
  subroutine initialize(self, filename)
!   < Initialize the equation.
    class(equation_singleideal_gpu_object), intent(inout) :: self              !< The equation.
    character(*) , intent(in) :: filename !< Input file name.
!
    call get_mpi_basic_info(self%nprocs, self%myrank, self%masterproc, self%mpi_err)
!
    call self%equation_base%initialize(filename)
    self%nx = self%equation_base%field%nx
    self%ny = self%equation_base%field%ny
    self%nz = self%equation_base%field%nz
    self%ng = self%equation_base%grid%ng
    self%nv = self%equation_base%nv
!
    self%visc_model       = self%equation_base%visc_model
    self%mu0              = self%equation_base%mu0
    self%T_ref_dim        = self%equation_base%T_ref_dim
    self%sutherland_S     = self%equation_base%sutherland_S
    self%powerlaw_vtexp   = self%equation_base%powerlaw_vtexp
!
    call self%point_to_field(self%equation_base%field)
    call self%point_to_grid(self%equation_base%grid)
    self%num_iter = self%equation_base%num_iter
    self%time0    = self%equation_base%time0
    self%icyc0    = self%equation_base%icyc0
!
    call self%base_gpu%initialize(self%equation_base%field)
!
    if(self%equation_base%debug_memory>0) then
      call self%field%check_cpu_mem(description="--Initialize-start--")
      call self%base_gpu%check_gpu_mem(description="--Initialize-start--")
    endif
!
    call self%base_gpu%copy_cpu_gpu()
!
    if(self%equation_base%debug_memory>0) then
      call self%field%check_cpu_mem(description="--Initialize-first-GPU-usage--")
      call self%base_gpu%check_gpu_mem(description="--Initialize-first-GPU-usage--")
    endif
!
!   self%nv = self%equation_base%field%nv
    self%nv_aux        =  10
!
    call self%alloc()
!
    self%w_aux_gpu = self%equation_base%w_aux
    self%fluid_mask_gpu = self%equation_base%fluid_mask
    self%ep_ord_change_gpu = self%equation_base%ep_ord_change
    self%ep_ord_change_x_gpu = self%equation_base%ep_ord_change_x
    self%winf_gpu       = self%equation_base%winf
    self%winf_past_shock_gpu = self%equation_base%winf_past_shock
!
    self%cp_coeff_gpu = self%equation_base%cp_coeff
    self%cv_coeff_gpu = self%equation_base%cv_coeff
!
    if (self%equation_base%num_probe>0) then
      self%probe_coeff_gpu = self%equation_base%probe_coeff
      self%ijk_probe_gpu   = self%equation_base%ijk_probe
    endif
!
    self%wmean_gpu = self%equation_base%wmean
!
    allocate(self%coeff_deriv1_gpu(1:4,4))
    allocate(self%coeff_deriv2_gpu(0:4,4))
    self%coeff_deriv1_gpu = self%equation_base%coeff_deriv1
    self%coeff_deriv2_gpu = self%equation_base%coeff_deriv2
!
    if (self%equation_base%recyc) then
      self%yplus_inflow_gpu = self%equation_base%yplus_inflow
      self%eta_inflow_gpu   = self%equation_base%eta_inflow
      self%yplus_recyc_gpu  = self%equation_base%yplus_recyc
      self%eta_recyc_gpu    = self%equation_base%eta_recyc
      self%map_j_inn_gpu    = self%equation_base%map_j_inn
      self%map_j_out_gpu    = self%equation_base%map_j_out
      self%weta_inflow_gpu  = self%equation_base%weta_inflow
    endif
!
!
!   
    if(self%equation_base%debug_memory>0) then
      call self%field%check_cpu_mem(description="Initialize-completed")
      call self%base_gpu%check_gpu_mem(description="Initialize-completed")
    endif
!   
!   Allocate field_gpu variables
!   call self%base_gpu%alloc(field=self%field, nv_aux=self%nv_aux)
!
!   ! Use base_gpu as pointee
!   self%field         => self%base_gpu%field
!   self%grid          => self%base_gpu%field%grid
!
    self%mpi_err = cudaStreamCreate(stream1)
  endsubroutine initialize
!
  subroutine alloc(self)
    class(equation_singleideal_gpu_object), intent(inout) :: self
    associate(nx => self%nx, ny => self%ny, nz => self%nz,  &
      ng => self%ng, nv => self%nv, nv_aux => self%nv_aux, &
      weno_scheme => self%equation_base%weno_scheme,&
      calorically_perfect => self%equation_base%calorically_perfect, &
      indx_cp_l => self%equation_base%indx_cp_l ,&
      indx_cp_r => self%equation_base%indx_cp_r)
!
      allocate(self%winf_gpu(nv))
      allocate(self%winf_past_shock_gpu(nv))
      allocate(self%w_aux_gpu(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng, nv_aux))
      allocate(self%fl_gpu(1:nx, 1:ny, 1:nz, nv))
      allocate(self%fln_gpu(1:nx, 1:ny, 1:nz, nv))
      allocate(self%wallprop_gpu(1-ng:nx+ng, 1-ng:nz+ng, 2:4))
      allocate(self%w_var(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng, 1))
      allocate(self%w_var_t(1, 1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng))
      allocate(self%fhat_gpu(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng, nv))
      allocate(self%w_aux_trans_gpu(1-ng:ny+ng, 1-ng:nx+ng, 1-ng:nz+ng, 8))
      allocate(self%fhat_trans_gpu(1-ng:ny+ng, 1-ng:nx+ng, 1-ng:nz+ng, nv))
      allocate(self%fl_trans_gpu(1:ny, 1:nx, 1:nz, nv))
      allocate(self%gplus_x_gpu (ny,nz,nv,2*weno_scheme))
      allocate(self%gminus_x_gpu(ny,nz,nv,2*weno_scheme))
      allocate(self%gplus_y_gpu (nx,nz,nv,2*weno_scheme))
      allocate(self%gminus_y_gpu(nx,nz,nv,2*weno_scheme))
      allocate(self%gplus_z_gpu (nx,ny,nv,2*weno_scheme))
      allocate(self%gminus_z_gpu(nx,ny,nv,2*weno_scheme))
      allocate(self%fluid_mask_gpu(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng))
      allocate(self%ep_ord_change_gpu(0:nx, 0:ny, 0:nz, 2:3))
      allocate(self%ep_ord_change_x_gpu(0:ny, 0:nx, 0:nz))
!
      allocate(self%wrecyc_gpu(ng,ny,nz,nv))
      allocate(self%wrecycav_gpu(ng,ny,nv))
      allocate(self%wmean_gpu(1-self%grid%ng:self%field%nx+self%grid%ng+1, 1:self%field%ny,4))
!
      allocate(self%yplus_inflow_gpu(1-ng:ny+ng))
      allocate(self%eta_inflow_gpu(1-ng:ny+ng))
      allocate(self%yplus_recyc_gpu(1-ng:ny+ng))
      allocate(self%eta_recyc_gpu(1-ng:ny+ng))
      allocate(self%map_j_inn_gpu(1:ny))
      allocate(self%map_j_out_gpu(1:ny))
      allocate(self%weta_inflow_gpu(1:ny))
      allocate(self%inflow_random_plane_gpu(1:ny,1:nz,3))
!
      allocate(self%cv_coeff_gpu(indx_cp_l:indx_cp_r+1))
      allocate(self%cp_coeff_gpu(indx_cp_l:indx_cp_r+1))
!
      if (self%equation_base%num_probe>0) then
        allocate(self%w_aux_probe_gpu(6,self%equation_base%num_probe))
        allocate(self%ijk_probe_gpu(3,self%equation_base%num_probe))
        allocate(self%probe_coeff_gpu(2,2,2,self%equation_base%num_probe))
      endif
!
    endassociate
  endsubroutine alloc
!
  subroutine compute_residual(self)
    class(equation_singleideal_gpu_object), intent(inout) :: self              !< The equation.
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv, &
      dt => self%equation_base%dt, fln_gpu => self%fln_gpu, residual_rhou => self%equation_base%residual_rhou, &
      fluid_mask_gpu => self%fluid_mask_gpu)
      call compute_residual_cuf(nx, ny, nz, ng, nv, fln_gpu, dt, residual_rhou, fluid_mask_gpu)
      call mpi_allreduce(MPI_IN_PLACE,residual_rhou,1,mpi_prec,mpi_sum,self%equation_base%field%mp_cart,self%mpi_err)
      residual_rhou = residual_rhou / (real(nx,rkind)*real(ny,rkind)*real(nz,rkind)*real(self%nprocs,rkind))
      residual_rhou = sqrt(residual_rhou)
    endassociate
  endsubroutine compute_residual
!
  subroutine compute_dt(self)
!   < Initialize the equation.
    class(equation_singleideal_gpu_object), intent(inout) :: self              !< The equation.
    real(rkind) :: dt_min, dtinv_max
    real(rkind) :: dtxi_max, dtyi_max, dtzi_max, dtxv_max, dtyv_max, dtzv_max, dtxk_max, dtyk_max, dtzk_max
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, &
      dt => self%equation_base%dt, CFL => self%equation_base%CFL, &
      Prandtl => self%equation_base%Prandtl, &
      visc_model => self%visc_model, mu0 => self%mu0, &
      T_ref_dim => self%T_ref_dim, powerlaw_vtexp => self%powerlaw_vtexp, &
      sutherland_S => self%sutherland_S, &
      w_gpu => self%base_gpu%w_gpu, &
      w_aux_gpu  => self%w_aux_gpu, &
      dcsidx_gpu => self%base_gpu%dcsidx_gpu,   &
      detady_gpu => self%base_gpu%detady_gpu,   &
      dzitdz_gpu => self%base_gpu%dzitdz_gpu,   &
      dcsidxs_gpu  => self%base_gpu%dcsidxs_gpu, &
      detadys_gpu  => self%base_gpu%detadys_gpu, &
      dzitdzs_gpu  => self%base_gpu%dzitdzs_gpu, &
      calorically_perfect => self%equation_base%calorically_perfect, &
      indx_cp_l    => self%equation_base%indx_cp_l, &
      indx_cp_r    => self%equation_base%indx_cp_r, &
      cp_coeff_gpu => self%cp_coeff_gpu, &
      fluid_mask_gpu => self%fluid_mask_gpu, &
      rgas0 => self%equation_base%rgas0, &
      t0 => self%equation_base%t0)
      if (CFL < 0) then
        dt = -CFL
      else
        call compute_dt_cuf(nx, ny, nz, ng, rgas0, Prandtl, &
          dcsidx_gpu, detady_gpu, dzitdz_gpu, dcsidxs_gpu, detadys_gpu, dzitdzs_gpu, w_gpu, w_aux_gpu, &
          dtxi_max, dtyi_max, dtzi_max, dtxv_max, dtyv_max, dtzv_max, dtxk_max, dtyk_max, dtzk_max,    &
          indx_cp_l, indx_cp_r, cp_coeff_gpu,fluid_mask_gpu,calorically_perfect,t0)
!       open(unit=116, file="dt_values.dat", position="append")
!       write(116,'(100(f16.8,2x))') dtxi_max, dtyi_max, dtzi_max, dtxv_max, dtyv_max, dtzv_max, dtxk_max, dtyk_max, dtzk_max
!       close(116)
        dtinv_max = maxval([dtxi_max, dtyi_max, dtzi_max, dtxv_max, dtyv_max, dtzv_max, dtxk_max, dtyk_max, dtzk_max])
        call mpi_allreduce(MPI_IN_PLACE,dtinv_max,1,mpi_prec,mpi_max,self%equation_base%field%mp_cart,self%mpi_err)
        dt_min = 1._rkind/dtinv_max
        dt = self%equation_base%CFL*dt_min
      endif
    endassociate
  endsubroutine compute_dt
!
  subroutine run(self, filename)
!
    class(equation_singleideal_gpu_object), intent(inout) :: self              !< The equation.
    character(*)                    , intent(in)          :: filename          !< Input file name.
    real(rkind)                                           :: timing(1:2)       !< Tic toc timing.
    real(rkind)                                           :: timing_step(1:2)  !< Tic toc timing.
    integer :: icyc_loop, iercuda
!
    call self%initialize(filename=filename)
!
    associate(icyc0 => self%equation_base%icyc0, icyc => self%equation_base%icyc, &
      time => self%equation_base%time, iter_dt_recompute => self%equation_base%iter_dt_recompute, &
      residual_rhou => self%equation_base%residual_rhou, dpdx => self%equation_base%dpdx, &
      rhobulk => self%equation_base%rhobulk, ubulk => self%equation_base%ubulk, &
      tbulk => self%equation_base%tbulk, nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, &
      visc_model => self%visc_model, mu0 => self%mu0, &
      T_ref_dim => self%T_ref_dim, &
      sutherland_S => self%sutherland_S, &
      powerlaw_vtexp => self%powerlaw_vtexp, &
      cv_coeff_gpu => self%cv_coeff_gpu, &
      calorically_perfect => self%equation_base%calorically_perfect, &
      indx_cp_l    => self%equation_base%indx_cp_l, &
      indx_cp_r    => self%equation_base%indx_cp_r, &
      mode_async => self%equation_base%mode_async, &
      time_from_last_rst    => self%equation_base%time_from_last_rst,    &
      time_from_last_write  => self%equation_base%time_from_last_write,  &
      time_from_last_stat   => self%equation_base%time_from_last_stat,   &
      time_from_last_slice  => self%equation_base%time_from_last_slice,  &
      time_is_freezed => self%equation_base%time_is_freezed)
!
      call self%update_ghost()
      call self%compute_aux()
      !@cuf iercuda=cudaDeviceSynchronize()
!
      if (mode_async >= 0) then
        if (self%equation_base%conservative_viscous==1) then
          call self%visflx(mode=2)
        else
          call self%visflx(mode=0)
        endif
        call self%base_gpu%bcswap_var(self%w_aux_gpu(:,:,:,8:8 )) ! ducros
      endif
      call zero_flux_cuf(self%nx,self%ny,self%nz,self%nv,self%fl_gpu)
!
      if (self%equation_base%restart_type==0) then
        self%equation_base%w_aux = self%w_aux_gpu
        if (self%equation_base%enable_plot3d>0) then
          call self%field%write_plot3d(mach=self%equation_base%Mach, reynolds=self%equation_base%Reynolds, &
            time=0._rkind, istore=0, plot3dgrid=.true., plot3dfield=.true., &
            w_aux_io=self%equation_base%w_aux(1:self%nx,1:self%ny,1:self%nz,1:6))
        endif
        if (self%equation_base%enable_vtk>0) then
          call self%field%write_vtk(time=0._rkind, istore=0, w_aux_io=self%equation_base%w_aux(1:self%nx,1:self%ny,1:self%nz,1:6))
        endif
      endif
      call self%compute_dt()
      if (self%masterproc) write(*,*) 'dt =', self%equation_base%dt
      if (self%masterproc) write(*,*) 'dt*u0/l0 =', self%equation_base%dt*self%equation_base%u0/self%equation_base%l0
!
      call MPI_BARRIER(MPI_COMM_WORLD, self%error) ; timing(1) = MPI_Wtime()
!     
      icyc_loop = icyc
      integration: do
!
        icyc_loop = icyc_loop + 1
!
        if ( time_is_freezed ) then
          self%equation_base%dt = 0._rkind
        else
!         call MPI_BARRIER(MPI_COMM_WORLD, self%error) ; timing_step(1) = MPI_Wtime()
          icyc = icyc + 1
!
          if(mod(icyc-icyc0, iter_dt_recompute)==0) then
            call self%compute_aux(central=1, ghost=0)
            !@cuf iercuda=cudaDeviceSynchronize()
            call self%compute_dt()
          endif
!
!         select case(self%equation_base%rk_type)
!         case(RK_WRAY,RK_JAMESON)
            if(mode_async == -1) call self%rk_sync_old()
            if(mode_async ==  0) call self%rk_sync()
            if(mode_async ==  1) call self%rk_async()
!         case(RK_SHU)
!           call self%rk()
!         end select
          if (mod(icyc-icyc0, self%equation_base%print_control)==0) call self%compute_residual()
          if (ieee_is_nan(self%equation_base%residual_rhou)) then
            if (self%masterproc) write(*,*) 'BOOM!!!'
            call mpi_barrier(mpi_comm_world,self%mpi_err)
            call mpi_abort(mpi_comm_world,99,self%mpi_err)
          endif
        endif
!       
        call self%manage_output()
!
        self%equation_base%time = self%equation_base%time + self%equation_base%dt
!
        if(self%masterproc.and.mod(icyc-icyc0, self%equation_base%print_control)==0) then
          call self%print_progress()
        endif
        if ((self%equation_base%icyc-self%equation_base%icyc0) >= self%num_iter) exit integration
        call MPI_BARRIER(MPI_COMM_WORLD, self%error) ; timing_step(2) = MPI_Wtime()
!       print '(A, F18.10)', 'step timing: ', timing_step(2) - timing_step(1)
      enddo integration
!
      if (allocated(self%equation_base%islice)) close(133)
      if (allocated(self%equation_base%jslice)) close(134)
      if (allocated(self%equation_base%kslice)) close(135)
      if (self%equation_base%num_probe>0) close(136)
!
      call MPI_BARRIER(MPI_COMM_WORLD, self%error) ; timing(2) = MPI_Wtime()
      if(self%num_iter > 0) then
        if (self%masterproc) then
          write(*,'(A, F18.10)') 'averaged timing: ', &
            (timing(2) - timing(1))/(self%equation_base%icyc-self%equation_base%icyc0)
        endif
      endif
!
      call self%base_gpu%copy_gpu_cpu()
      if (self%equation_base%io_type==1) then
        call self%field%write_field_serial()
        call self%equation_base%write_stats_serial()
        if (self%equation_base%enable_stat_3d>0) call self%equation_base%write_stats_3d_serial()
      endif
      if (self%equation_base%io_type==2) then
        call self%field%write_field()
        call self%equation_base%write_stats()
        if (self%equation_base%enable_stat_3d>0) call self%equation_base%write_stats_3d()
      endif
      call self%equation_base%write_field_info()
!
    endassociate
  endsubroutine run
!
  subroutine print_progress(self)
    class(equation_singleideal_gpu_object), intent(inout) :: self
    character(6) :: pos_io
    associate(icyc => self%equation_base%icyc, time => self%equation_base%time, dt => self%equation_base%dt, &
      residual_rhou => self%equation_base%residual_rhou, dpdx => self%equation_base%dpdx, &
      rhobulk => self%equation_base%rhobulk, ubulk => self%equation_base%ubulk, &
      tbulk => self%equation_base%tbulk)
!     
      residual_rhou = residual_rhou/(self%equation_base%u0**2)/self%equation_base%rho0*self%equation_base%l0
      dpdx          = dpdx         /(self%equation_base%u0**2)/self%equation_base%rho0*self%equation_base%l0
      pos_io = 'append'
      if (self%equation_base%icyc==1) pos_io = 'rewind'
      if (self%masterproc) then
        open(unit=15,file='progress.out',position=pos_io)
        if (self%equation_base%channel_case) then
          write(* ,'(1I10,20ES20.10)') icyc, dt, time, residual_rhou, -dpdx/dt, rhobulk, ubulk, tbulk
          write(15,'(1I10,20ES20.10)') icyc, dt, time, residual_rhou, -dpdx/dt, rhobulk, ubulk, tbulk
        else
          write(* ,'(1I10,20ES20.10)') icyc, dt, time, residual_rhou
          write(15,'(1I10,20ES20.10)') icyc, dt, time, residual_rhou
        endif
        close(15)
      endif
    endassociate
  endsubroutine print_progress
!
!
!
!
!
  subroutine manage_output(self)
    class(equation_singleideal_gpu_object), intent(inout) :: self
!   
    integer :: i,j,k,ii,jj,kk,l,n
    integer :: isize, jsize, ksize
    character(3) :: chx, chy, chz
    logical :: sliceyz_exist, slicexz_exist, slicexy_exist, probe_exist
!   
    associate(time_from_last_rst    => self%equation_base%time_from_last_rst,    &
      time_from_last_write  => self%equation_base%time_from_last_write,  &
      time_from_last_stat   => self%equation_base%time_from_last_stat,   &
      time_from_last_slice  => self%equation_base%time_from_last_slice,  &
      icyc0 => self%equation_base%icyc0, icyc => self%equation_base%icyc, &
      time => self%equation_base%time, &
      w_aux_gpu => self%w_aux_gpu, &
      ijk_probe_gpu => self%ijk_probe_gpu, &
      w_aux_probe_gpu => self%w_aux_probe_gpu, &
      probe_coeff_gpu => self%probe_coeff_gpu, &
      nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv_aux => self%nv_aux, &
      w_aux_probe => self%equation_base%w_aux_probe)
!
!     Save flow samples
      time_from_last_write = time_from_last_write + self%equation_base%dt
      if (time_from_last_write >= self%equation_base%dtsave) then
        if (self%masterproc) write(*,*) 'time_from_last_write=',time_from_last_write
        if (self%masterproc) write(*,*) 'istore =', self%equation_base%istore
        call self%base_gpu%copy_gpu_cpu()
        self%equation_base%w_aux = self%w_aux_gpu
        if (self%equation_base%enable_plot3d>0) then
          call self%field%write_plot3d(mach=self%equation_base%Mach, reynolds=self%equation_base%Reynolds, &
            time=time, istore=self%equation_base%istore, plot3dgrid=.false., plot3dfield=.true., &
            w_aux_io=self%equation_base%w_aux(1:self%nx,1:self%ny,1:self%nz,1:6))
        endif
        if (self%equation_base%enable_vtk>0) then
          call self%field%write_vtk(time=time, istore=self%equation_base%istore, &
            w_aux_io=self%equation_base%w_aux(1:self%nx,1:self%ny,1:self%nz,1:6))
        endif
        time_from_last_write = 0._rkind !time_from_last_write - self%equation_base%dtsave
        self%equation_base%istore = self%equation_base%istore + 1
      endif
!
!     Compute stats
      time_from_last_stat = time_from_last_stat + self%equation_base%dt
      if (time_from_last_stat >= self%equation_base%dtstat) then
        if (self%masterproc) write(*,*) 'time_from_last_stat=',time_from_last_stat
        if (self%masterproc) write(*,*) 'itav =',self%equation_base%itav
        call self%base_gpu%copy_gpu_cpu()
        self%equation_base%w_aux(:,:,:,6) = self%w_aux_gpu(:,:,:,6)
        self%equation_base%w_aux(:,:,:,7) = self%w_aux_gpu(:,:,:,7)
        call self%equation_base%compute_stats()
        if (self%equation_base%enable_stat_3d>0) call self%equation_base%compute_stats_3d()
        time_from_last_stat = 0._rkind ! time_from_last_stat - self%equation_base%dtstat
        self%equation_base%itav = self%equation_base%itav + 1
      endif
!
!
!
!     Write slice
      time_from_last_slice = time_from_last_slice + self%equation_base%dt
      write(chx,'(I3.3)') self%field%ncoords(1)
      write(chy,'(I3.3)') self%field%ncoords(2)
      write(chz,'(I3.3)') self%field%ncoords(3)
      if (icyc-icyc0 == 1) then
        if (allocated(self%equation_base%islice)) then
          inquire(file='sliceyz_'//chx//'_'//chy//'_'//chz//'.bin', exist=sliceyz_exist)
          if(.not.(sliceyz_exist)) then
            open(133,file='sliceyz_'//chx//'_'//chy//'_'//chz//'.bin',form='unformatted',                   asynchronous="yes")
          else
            open(133,file='sliceyz_'//chx//'_'//chy//'_'//chz//'.bin',form='unformatted', position="append",asynchronous="yes")
          endif
        endif
        if (allocated(self%equation_base%jslice)) then
          inquire(file='slicexz_'//chx//'_'//chy//'_'//chz//'.bin', exist=slicexz_exist)
          if(.not.(slicexz_exist)) then
            open(134,file='slicexz_'//chx//'_'//chy//'_'//chz//'.bin',form='unformatted',                   asynchronous="yes")
          else
            open(134,file='slicexz_'//chx//'_'//chy//'_'//chz//'.bin',form='unformatted', position="append",asynchronous="yes")
          endif
        endif
        if (allocated(self%equation_base%kslice)) then
          inquire(file='slicexy_'//chx//'_'//chy//'_'//chz//'.bin', exist=slicexy_exist)
          if(.not.(slicexy_exist)) then
            open(135,file='slicexy_'//chx//'_'//chy//'_'//chz//'.bin',form='unformatted',                   asynchronous="yes")
          else
            open(135,file='slicexy_'//chx//'_'//chy//'_'//chz//'.bin',form='unformatted', position="append",asynchronous="yes")
          endif
        endif
        if (self%equation_base%num_probe>0) then
          inquire(file='probe_'//chx//'_'//chy//'_'//chz//'.bin', exist=probe_exist)
          if(.not.(probe_exist)) then
            open(136,file='probe_'//chx//'_'//chy//'_'//chz//'.bin',form='unformatted')
            write(136) self%equation_base%num_probe
          else
            open(136,file='probe_'//chx//'_'//chy//'_'//chz//'.bin',form='unformatted', position="append")
          endif
        endif
      endif
!     
      if (time_from_last_slice >= self%equation_base%dtslice) then
        if (self%masterproc) write(*,*) 'time_from_last_slice=',time_from_last_slice
        if (allocated(self%equation_base%islice)) then
          isize = size(self%equation_base%islice)
          do i=1,size(self%equation_base%islice)
            ii = self%equation_base%islice(i)
            sliceyz_aux(i,:,:,1:6) = self%w_aux_gpu(ii,:,:,1:6)
          enddo
!         wait(133)
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
            slicexz_aux(:,j,:,1:6) = self%w_aux_gpu(:,jj,:,1:6)
          enddo
          if (self%equation_base%jslice(1)==1) slicexz_aux(:,1,:,2:4) = self%wallprop_gpu(:,:,2:4)
!         wait(134)
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
            slicexy_aux(:,:,k,1:6) = self%w_aux_gpu(:,:,kk,1:6)
          enddo
!         wait(135)
          write(135,asynchronous="no") icyc,time
          write(135,asynchronous="no") self%grid%nxmax,self%grid%nymax,self%grid%nzmax
          write(135,asynchronous="no") ksize,(self%equation_base%kslice(k),k=1,ksize)
          write(135,asynchronous="no") self%field%nx, self%field%ny, ksize, 6
          write(135,asynchronous="yes") slicexy_aux !(1:self%field%nx,1:self%field%ny,1:ksize,1:6)
        endif
!       
        if (self%equation_base%num_probe>0) then
          call probe_interpolation_cuf(self%equation_base%num_probe,nx,ny,nz,ng,nv_aux,ijk_probe_gpu,&
            w_aux_probe_gpu,w_aux_gpu,probe_coeff_gpu)
!         
          w_aux_probe = w_aux_probe_gpu
!         write(136,100) time, ((w_aux_probe(l,n),l=1,6),n=1,self%equation_base%num_probe)
          write(136) time, ((w_aux_probe(l,n),l=1,6),n=1,self%equation_base%num_probe)
        endif
!       
        time_from_last_slice = 0._rkind ! time_from_last_slice - self%equation_base%dtslice
      endif
!
!     Save restart & stats
      time_from_last_rst = time_from_last_rst + self%equation_base%dt
      if (time_from_last_rst >= self%equation_base%dtsave_restart) then
        if (self%masterproc) write(*,*) 'time_from_last_rst=',time_from_last_rst
        call self%base_gpu%copy_gpu_cpu()
!
        if (self%equation_base%io_type==1) then
          call self%field%write_field_serial()
          call self%equation_base%write_stats_serial()
          if (self%equation_base%enable_stat_3d>0) call self%equation_base%write_stats_3d_serial()
        endif
        if (self%equation_base%io_type==2) then
          call self%field%write_field()
          call self%equation_base%write_stats()
          if (self%equation_base%enable_stat_3d>0) call self%equation_base%write_stats_3d()
        endif
        call self%equation_base%write_field_info()
        time_from_last_rst = 0._rkind ! time_from_last_rst - self%equation_base%dtsave_restart
      endif
!
    endassociate
!   
  end subroutine manage_output
!
endmodule streams_equation_singleideal_gpu_object
!
