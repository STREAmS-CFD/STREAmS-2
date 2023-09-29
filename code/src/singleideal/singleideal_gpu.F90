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
  use, intrinsic :: iso_fortran_env
  use tcp
!
! INSITU CATALYST2
  use catalyst_api
  use catalyst_conduit
!
  implicit none
  private
  public :: equation_singleideal_gpu_object
!
  integer(kind=cuda_stream_kind) :: stream1
!
  integer, parameter :: EULERCENTRAL_THREADS_X=128,EULERCENTRAL_THREADS_Y=3
  integer, parameter :: EULERWENO_THREADS_X=128 ,EULERWENO_THREADS_Y=2
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
    type(base_gpu_object) :: base_gpu !< The base GPU handler.
    type(equation_singleideal_object) :: equation_base !< The equation base.
    type(field_object), pointer :: field=>null() !< The field.
    type(grid_object), pointer :: grid=>null() !< The grid.
    integer(ikind) :: ng !< Number of ghost cells.
    integer(ikind) :: nx !< Number of cells in i direction.
    integer(ikind) :: ny !< Number of cells in j direction.
    integer(ikind) :: nz !< Number of cells in k direction.
    integer(ikind) :: nv !< Number of variables.
    integer(ikind) :: nv_aux !< Number of auxiliary variables.
    integer(ikind) :: nprocs !< Number of auxiliary variables.
    integer(ikind) :: myrank !< Number of auxiliary variables.
    integer(ikind) :: error !< Number of auxiliary variables.
    real(rkind) :: time0
    real(rkind), allocatable, dimension(:,:), device :: coeff_deriv1_gpu
    real(rkind), allocatable, dimension(:,:), device :: coeff_deriv2_gpu
!
    integer :: ierr
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
!   real(rkind), allocatable, dimension(:,:,:,:), device :: gplus_x_gpu, gminus_x_gpu
!   real(rkind), allocatable, dimension(:,:,:,:), device :: gplus_y_gpu, gminus_y_gpu
!   real(rkind), allocatable, dimension(:,:,:,:), device :: gplus_z_gpu, gminus_z_gpu
    integer, allocatable, dimension(:,:,:), device :: fluid_mask_gpu
    integer, allocatable, dimension(:,:,:,:), device :: ep_ord_change_gpu
!
    real(rkind), allocatable, dimension(:,:,:,:), device :: wrecyc_gpu
    real(rkind), allocatable, dimension(:,:,:), device :: wrecycav_gpu
    real(rkind), dimension(:,:,:), allocatable, device :: wmean_gpu
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
!   ibm_var_start
    integer, dimension(:,:,:), allocatable, device :: ibm_sbody_gpu
    integer, dimension(:,:,:), allocatable, device :: ibm_is_interface_node_gpu
    integer, dimension(:,:,:), allocatable, device :: ibm_inside_moving_gpu
    integer, dimension(:,: ), allocatable, device :: ibm_ijk_interface_gpu
    integer, dimension(:,: ), allocatable, device :: ibm_ijk_refl_gpu
    integer, dimension(:,: ), allocatable, device :: ibm_ijk_wall_gpu
    real(rkind), dimension(:,:), allocatable, device :: ibm_nxyz_interface_gpu
    integer, dimension(:,:), allocatable, device :: ibm_bc_gpu
    real(rkind), dimension(:,:), allocatable, device :: ibm_dist_gpu
    real(rkind), dimension(:,:,:,:), allocatable, device :: ibm_coeff_d_gpu
    real(rkind), dimension(:,:,:,:), allocatable, device :: ibm_coeff_n_gpu
    integer, dimension(: ), allocatable, device :: ibm_refl_type_gpu
    real(rkind), dimension(:,: ), allocatable, device :: ibm_w_refl_gpu
    real(rkind), dimension(:,: ), allocatable, device :: ibm_parbc_gpu
    real(rkind), dimension(: ), allocatable, device :: ibm_vega_y_gpu, ibm_vega_r_gpu
!   ibm_var_end
!
!   insitu_var_start
!   real(rkind), dimension(:,:,:,:), allocatable, managed :: psi_pv_managed
    real(rkind), dimension(:,:,:,:), allocatable, device :: psi_gpu
    integer, dimension(:), allocatable, device :: aux_list_gpu, add_list_gpu
!   insitu_var_end
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
!   ibm
    procedure, pass(self) :: ibm_apply
    procedure, pass(self) :: ibm_inside
    procedure, pass(self) :: ibm_compute_force
    procedure, pass(self) :: ibm_alloc_gpu
!   insitu
    procedure, pass(self) :: insitu_alloc_gpu
    procedure, pass(self) :: insitu_coprocess
    procedure, pass(self) :: insitu_compute_psi
    procedure, pass(self) :: insitu_do_catalyst_execute
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
    class(equation_singleideal_gpu_object), intent(inout) :: self !< The equation.
    class(field_object), target :: field !< The equation.
    self%field => field
  endsubroutine point_to_field
!
  subroutine point_to_grid(self, grid)
    class(equation_singleideal_gpu_object), intent(inout) :: self !< The equation.
    class(grid_object), target :: grid !< The equation.
    self%grid => grid
  endsubroutine point_to_grid
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
  subroutine ibm_alloc_gpu(self)
    class(equation_singleideal_gpu_object), intent(inout) :: self
    associate(ibm_num_interface => self%equation_base%ibm_num_interface, ibm_num_bc => self%equation&
    &_base%ibm_num_bc, nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv)
      allocate(self%ibm_sbody_gpu(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
      allocate(self%ibm_is_interface_node_gpu(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
      if (self%equation_base%ibm_vega_moving>0) then
        allocate(self%ibm_inside_moving_gpu(nx,ny,nz))
!       notneeded self%ibm_inside_moving_gpu = 0
        allocate(self%ibm_vega_y_gpu(1:self%equation_base%ibm_vega_ny+1))
        self%ibm_vega_y_gpu = self%equation_base%ibm_vega_y
        allocate(self%ibm_vega_r_gpu(1:self%equation_base%ibm_vega_ny+1))
        self%ibm_vega_r_gpu = self%equation_base%ibm_vega_r
      endif
      allocate(self%ibm_parbc_gpu(ibm_num_bc,ibm_MAX_PARBC))
      if (ibm_num_interface>0) then
        allocate(self%ibm_ijk_interface_gpu (3,ibm_num_interface)) ! Local values of i,j,k for the interface node
        allocate(self%ibm_ijk_refl_gpu (3,ibm_num_interface)) ! Reflected node bwtween i,i+1 and j,j+1 and k,k+1
        allocate(self%ibm_ijk_wall_gpu (3,ibm_num_interface)) ! Wall node between i,i+1 and j,j+1 and k,k+1
        allocate(self%ibm_nxyz_interface_gpu(3,ibm_num_interface)) ! Wall-normal components
        allocate(self%ibm_bc_gpu (2,ibm_num_interface)) ! Bc tag for interface nodes
!       Distance between interface node and wall point (1) and reflected point and wall point (2)
        allocate(self%ibm_dist_gpu (2,ibm_num_interface))
        allocate(self%ibm_coeff_d_gpu (2,2,2,ibm_num_interface)) ! Coefficients for trilin interpolation (Dirichlet)
        allocate(self%ibm_coeff_n_gpu (2,2,2,ibm_num_interface)) ! Coefficients for trilin interpolation (Neumann)
        allocate(self%ibm_refl_type_gpu(ibm_num_interface))
        allocate(self%ibm_w_refl_gpu(ibm_num_interface,nv))
      endif
      self%ibm_sbody_gpu = self%equation_base%ibm_sbody
      self%ibm_is_interface_node_gpu = self%equation_base%ibm_is_interface_node
      self%ibm_parbc_gpu = self%equation_base%ibm_parbc
      if (ibm_num_interface>0) then
        self%ibm_ijk_interface_gpu = self%equation_base%ibm_ijk_interface
        self%ibm_ijk_refl_gpu = self%equation_base%ibm_ijk_refl
        self%ibm_ijk_wall_gpu = self%equation_base%ibm_ijk_wall
        self%ibm_nxyz_interface_gpu = self%equation_base%ibm_nxyz_interface
        self%ibm_bc_gpu = self%equation_base%ibm_bc
        self%ibm_dist_gpu = self%equation_base%ibm_dist
        self%ibm_coeff_d_gpu = self%equation_base%ibm_coeff_d
        self%ibm_coeff_n_gpu = self%equation_base%ibm_coeff_n
        self%ibm_refl_type_gpu = self%equation_base%ibm_refl_type
      endif
    endassociate
  endsubroutine ibm_alloc_gpu
!
  subroutine insitu_alloc_gpu(self)
    class(equation_singleideal_gpu_object), intent(inout) :: self
    associate( nxsl_ins => self%equation_base%nxsl_ins, nxel_ins => self%equation_base%nxel_ins, nys&
    &l_ins => self%equation_base%nysl_ins, nyel_ins => self%equation_base%nyel_ins, nzsl_ins => self%equa&
    &tion_base%nzsl_ins, nzel_ins => self%equation_base%nzel_ins, npsi => self%equation_base%npsi, npsi_p&
    &v => self%equation_base%npsi_pv, nx => self%field%nx, ny => self%field%ny, nz => self%field%nz, ng =&
    &> self%grid%ng, n_aux_list => self%equation_base%n_aux_list, n_add_list => self%equation_base%n_add_&
    &list )
!     NOMANAGED if (npsi_pv > 0) allocate(self%psi_pv_managed(nxsl_ins:nxel_ins,nysl_ins:nyel_ins,nzsl_ins:nzel_ins,npsi_pv))
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
!
  subroutine rk_sync_old(self)
    class(equation_singleideal_gpu_object), intent(inout) :: self !< The equation.
    integer :: istep, iercuda
    real(rkind) :: rhodt, gamdt, alpdt
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv, nv_aux => s&
    &elf%nv_aux, dt => self%equation_base%dt, eul_imin => self%equation_base%eul_imin, eul_imax => self%e&
    &quation_base%eul_imax, eul_jmin => self%equation_base%eul_jmin, eul_jmax => self%equation_base%eul_j&
    &max, eul_kmin => self%equation_base%eul_kmin, eul_kmax => self%equation_base%eul_kmax, channel_case &
    &=> self%equation_base%channel_case, enable_ibm => self%equation_base%enable_ibm)
!     
      if (enable_ibm>0) then
        self%equation_base%ibm_force_x = 0._rkind
        self%equation_base%ibm_force_y = 0._rkind
        self%equation_base%ibm_force_z = 0._rkind
      endif
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
        call self%euler_x(eul_imin, eul_imax)
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
        if (enable_ibm>0) then
          call self%base_gpu%bcswap_corner()
          call self%compute_aux()
          !@cuf iercuda=cudaDeviceSynchronize()
          if (istep==1) then
            call self%equation_base%ibm_bc_prepare()
            if (self%equation_base%ibm_trajectory_points>0) self%ibm_vega_y_gpu = self%equation_base%ibm_vega_y
            if (self%equation_base%ibm_vega_moving>0) call self%ibm_inside()
          endif
          call self%ibm_apply()
          call self%ibm_compute_force(istep)
          call self%update_ghost() ! needed after application of ibm
        endif
        call self%compute_aux()
        !@cuf iercuda=cudaDeviceSynchronize()
      enddo
    endassociate
  endsubroutine rk_sync_old
!
  subroutine visflx(self, mode)
    class(equation_singleideal_gpu_object), intent(inout) :: self !< The equation.
    integer, intent(in) :: mode
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv, nv_aux => s&
    &elf%nv_aux, dt => self%equation_base%dt, visc_order => self%equation_base%visc_order, Prandtl => sel&
    &f%equation_base%Prandtl, visc_model => self%visc_model, mu0 => self%mu0, u0 => self%equation_base%u0&
    &, l0 => self%equation_base%l0, T_ref_dim => self%T_ref_dim, sutherland_S => self%sutherland_S, power&
    &law_vtexp => self%powerlaw_vtexp, coeff_deriv1_gpu => self%coeff_deriv1_gpu, coeff_deriv2_gpu => sel&
    &f%coeff_deriv2_gpu, fhat_trans_gpu => self%fhat_trans_gpu, fl_trans_gpu => self%fl_trans_gpu, fl_gpu&
    & => self%fl_gpu, w_aux_gpu => self%w_aux_gpu, w_aux_trans_gpu => self%w_aux_trans_gpu, dcsidx_gpu =>&
    & self%base_gpu%dcsidx_gpu, detady_gpu => self%base_gpu%detady_gpu, dzitdz_gpu => self%base_gpu%dzitd&
    &z_gpu, dcsidxs_gpu => self%base_gpu%dcsidxs_gpu, detadys_gpu => self%base_gpu%detadys_gpu, dzitdzs_g&
    &pu => self%base_gpu%dzitdzs_gpu, dcsidx2_gpu => self%base_gpu%dcsidx2_gpu, detady2_gpu => self%base_&
    &gpu%detady2_gpu, dzitdz2_gpu => self%base_gpu%dzitdz2_gpu, x_gpu => self%base_gpu%x_gpu, y_gpu => se&
    &lf%base_gpu%y_gpu, z_gpu => self%base_gpu%z_gpu, eul_imin => self%equation_base%eul_imin, eul_imax =&
    &> self%equation_base%eul_imax, eul_jmin => self%equation_base%eul_jmin, eul_jmax => self%equation_ba&
    &se%eul_jmax, eul_kmin => self%equation_base%eul_kmin, eul_kmax => self%equation_base%eul_kmax, cv_co&
    &eff_gpu => self%cv_coeff_gpu, cp_coeff_gpu => self%cp_coeff_gpu, indx_cp_l => self%equation_base%ind&
    &x_cp_l, indx_cp_r => self%equation_base%indx_cp_r, calorically_perfect => self%equation_base%caloric&
    &ally_perfect, t0 => self%equation_base%t0, channel_case => self%equation_base%channel_case, enable_i&
    &bm => self%equation_base%enable_ibm)
      if (mode == 0) then ! laplacian
        call visflx_cuf(nx, ny, nz, ng, visc_order, Prandtl, t0, indx_cp_l, indx_cp_r, cp_coeff_gpu,&
        & calorically_perfect, u0, l0, self%base_gpu%w_gpu, self%w_aux_gpu, self%fl_gpu, coeff_deriv1_gpu, co&
        &eff_deriv2_gpu, dcsidx_gpu, detady_gpu, dzitdz_gpu, dcsidxs_gpu, detadys_gpu, dzitdzs_gpu, dcsidx2_g&
        &pu, detady2_gpu, dzitdz2_gpu, self%wallprop_gpu)
      elseif (mode == 1) then ! div
        call visflx_div_cuf(nx, ny, nz, ng, visc_order, self%w_aux_gpu, self%fl_gpu, coeff_deriv1_gp&
        &u, dcsidx_gpu, detady_gpu, dzitdz_gpu, stream1)
      elseif (mode == 2) then ! reduced
!       call visflx_reduced_cuf(nx, ny, nz, ng, visc_order, ! Prandtl, t0, indx_cp_l, indx_cp_r, cp_
!coeff_gpu, calorically_perfect, ! u0, l0, self%base_gpu%w_gpu, self%w_aux_gpu, self%fl_gpu, ! coeff_
!deriv1_gpu, coeff_deriv2_gpu, ! dcsidx_gpu, detady_gpu, dzitdz_gpu, ! dcsidxs_gpu, detadys_gpu, dzit
!dzs_gpu, ! dcsidx2_gpu, detady2_gpu, dzitdz2_gpu, self%wallprop_gpu, 1)
        call visflx_reduced_ord2_cuf(nx, ny, nz, ng, Prandtl, t0, indx_cp_l, indx_cp_r, cp_coeff_gpu&
        &, calorically_perfect, u0, l0, self%base_gpu%w_gpu, self%w_aux_gpu, self%fl_gpu, x_gpu, y_gpu, z_gpu&
        &, self%wallprop_gpu, 1)
      elseif (mode == 3) then ! only sensor
        call sensor_cuf(nx, ny, nz, ng, u0, l0, self%w_aux_gpu)
      elseif (mode == 4) then ! laplacian
        call visflx_nosensor_cuf(nx, ny, nz, ng, visc_order, Prandtl, t0, indx_cp_l, indx_cp_r, cp_c&
        &oeff_gpu, calorically_perfect, u0, l0, self%base_gpu%w_gpu, self%w_aux_gpu, self%fl_gpu, coeff_deriv&
        &1_gpu, coeff_deriv2_gpu, dcsidx_gpu, detady_gpu, dzitdz_gpu, dcsidxs_gpu, detadys_gpu, dzitdzs_gpu, &
        &dcsidx2_gpu, detady2_gpu, dzitdz2_gpu, self%wallprop_gpu)
      elseif (mode == 5) then ! staggered
        call visflx_x_cuf(nx, ny, nz, nv, ng, Prandtl, t0, indx_cp_l, indx_cp_r, cp_coeff_gpu, calor&
        &ically_perfect, self%base_gpu%x_gpu, w_aux_gpu, self%fl_gpu, self%fhat_gpu)
        call visflx_y_cuf(nx, ny, nz, nv, ng, Prandtl, t0, indx_cp_l, indx_cp_r, cp_coeff_gpu, calor&
        &ically_perfect, self%base_gpu%y_gpu, self%w_aux_gpu, self%fl_gpu)
        call visflx_z_cuf(nx, ny, nz, nv, ng, Prandtl, t0, indx_cp_l, indx_cp_r, cp_coeff_gpu, calor&
        &ically_perfect, self%base_gpu%z_gpu, self%w_aux_gpu, self%fl_gpu)
      elseif (mode == 6) then ! reduce_nosensor
!       call visflx_reduced_cuf(nx, ny, nz, ng, visc_order, ! Prandtl, t0, indx_cp_l, indx_cp_r, cp_
!coeff_gpu, calorically_perfect, ! u0, l0, self%base_gpu%w_gpu, self%w_aux_gpu, self%fl_gpu, ! coeff_
!deriv1_gpu, coeff_deriv2_gpu, ! dcsidx_gpu, detady_gpu, dzitdz_gpu, ! dcsidxs_gpu, detadys_gpu, dzit
!dzs_gpu, ! dcsidx2_gpu, detady2_gpu, dzitdz2_gpu, self%wallprop_gpu, 0)
        call visflx_reduced_ord2_cuf(nx, ny, nz, ng, Prandtl, t0, indx_cp_l, indx_cp_r, cp_coeff_gpu&
        &, calorically_perfect, u0, l0, self%base_gpu%w_gpu, self%w_aux_gpu, self%fl_gpu, x_gpu, y_gpu, z_gpu&
        &, self%wallprop_gpu, 0)
      elseif (mode == 7) then ! div_ord2
        call visflx_div_ord2_cuf(nx, ny, nz, ng, self%w_aux_gpu, self%fl_gpu, x_gpu, y_gpu, z_gpu, stream1)
      endif
    endassociate
  endsubroutine visflx
!
  subroutine compute_aux(self, central, ghost)
    class(equation_singleideal_gpu_object), intent(inout) :: self !< The equation.
    integer, intent(in), optional :: central, ghost
    integer :: central_, ghost_
    central_ = 1 ; if(present(central)) central_ = central
    ghost_ = 1 ; if(present(ghost)) ghost_ = ghost
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, visc_model => self%visc_mo&
    &del, mu0 => self%mu0, t0 => self%equation_base%t0, T_ref_dim => self%T_ref_dim, sutherland_S => self&
    &%sutherland_S, powerlaw_vtexp => self%powerlaw_vtexp, cv_coeff_gpu => self%cv_coeff_gpu, cp_coeff_gp&
    &u => self%cp_coeff_gpu, indx_cp_l => self%equation_base%indx_cp_l, indx_cp_r => self%equation_base%i&
    &ndx_cp_r, calorically_perfect => self%equation_base%calorically_perfect, rgas0 => self%equation_base&
    &%rgas0)
      if(central_ == 1 .and. ghost_ == 1) then
        call eval_aux_cuf(nx, ny, nz, ng, 1-ng, nx+ng, 1-ng, ny+ng, 1-ng, nz+ng, self%base_gpu%w_gpu&
        &, self%w_aux_gpu, visc_model, mu0, t0, sutherland_S, T_ref_dim, powerlaw_vtexp, VISC_POWER, VISC_SUT&
        &HERLAND, VISC_NO, cv_coeff_gpu, indx_cp_l, indx_cp_r, rgas0, calorically_perfect, tol_iter_nr,stream&
        &1)
      elseif(central_ == 1 .and. ghost_ == 0) then
        call eval_aux_cuf(nx, ny, nz, ng, 1, nx, 1, ny, 1, nz, self%base_gpu%w_gpu, self%w_aux_gpu, &
        &visc_model, mu0, t0, sutherland_S, T_ref_dim, powerlaw_vtexp, VISC_POWER, VISC_SUTHERLAND, VISC_NO, &
        &cv_coeff_gpu, indx_cp_l, indx_cp_r, rgas0, calorically_perfect, tol_iter_nr,stream1)
      elseif(central_ == 0 .and. ghost_ == 1) then
        call eval_aux_cuf(nx, ny, nz, ng, 1-ng, 0, 1-ng, ny+ng, 1-ng, nz+ng, self%base_gpu%w_gpu, se&
        &lf%w_aux_gpu, visc_model, mu0, t0, sutherland_S, T_ref_dim, powerlaw_vtexp, VISC_POWER, VISC_SUTHERL&
        &AND, VISC_NO, cv_coeff_gpu, indx_cp_l, indx_cp_r, rgas0, calorically_perfect, tol_iter_nr,stream1)
        call eval_aux_cuf(nx, ny, nz, ng, nx+1, nx+ng, 1-ng, ny+ng, 1-ng, nz+ng, self%base_gpu%w_gpu&
        &, self%w_aux_gpu, visc_model, mu0, t0, sutherland_S, T_ref_dim, powerlaw_vtexp, VISC_POWER, VISC_SUT&
        &HERLAND, VISC_NO, cv_coeff_gpu, indx_cp_l, indx_cp_r, rgas0, calorically_perfect, tol_iter_nr,stream&
        &1)
        call eval_aux_cuf(nx, ny, nz, ng, 1-ng, nx+ng, 1-ng, 0, 1-ng, nz+ng, self%base_gpu%w_gpu, se&
        &lf%w_aux_gpu, visc_model, mu0, t0, sutherland_S, T_ref_dim, powerlaw_vtexp, VISC_POWER, VISC_SUTHERL&
        &AND, VISC_NO, cv_coeff_gpu, indx_cp_l, indx_cp_r, rgas0, calorically_perfect, tol_iter_nr,stream1)
        call eval_aux_cuf(nx, ny, nz, ng, 1-ng, nx+ng, ny+1, ny+ng, 1-ng, nz+ng, self%base_gpu%w_gpu&
        &, self%w_aux_gpu, visc_model, mu0, t0, sutherland_S, T_ref_dim, powerlaw_vtexp, VISC_POWER, VISC_SUT&
        &HERLAND, VISC_NO, cv_coeff_gpu, indx_cp_l, indx_cp_r, rgas0, calorically_perfect, tol_iter_nr,stream&
        &1)
        call eval_aux_cuf(nx, ny, nz, ng, 1-ng, nx+ng, 1-ng, ny+ng, 1-ng, 0, self%base_gpu%w_gpu, se&
        &lf%w_aux_gpu, visc_model, mu0, t0, sutherland_S, T_ref_dim, powerlaw_vtexp, VISC_POWER, VISC_SUTHERL&
        &AND, VISC_NO, cv_coeff_gpu, indx_cp_l, indx_cp_r, rgas0, calorically_perfect, tol_iter_nr,stream1)
        call eval_aux_cuf(nx, ny, nz, ng, 1-ng, nx+ng, 1-ng, ny+ng, nz+1, nz+ng, self%base_gpu%w_gpu&
        &, self%w_aux_gpu, visc_model, mu0, t0, sutherland_S, T_ref_dim, powerlaw_vtexp, VISC_POWER, VISC_SUT&
        &HERLAND, VISC_NO, cv_coeff_gpu, indx_cp_l, indx_cp_r, rgas0, calorically_perfect, tol_iter_nr,stream&
        &1)
      endif
    endassociate
  endsubroutine compute_aux
!
  subroutine rk_sync(self)
    class(equation_singleideal_gpu_object), intent(inout) :: self !< The equation.
    integer :: istep, lmax, iercuda
    real(rkind) :: rhodt, gamdt, alpdt
!   
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv, dt => self%&
    &equation_base%dt, ep_order => self%equation_base%ep_order, weno_scheme => self%equation_base%weno_sc&
    &heme, conservative_viscous => self%equation_base%conservative_viscous, eul_imin => self%equation_bas&
    &e%eul_imin, eul_imax => self%equation_base%eul_imax, eul_jmin => self%equation_base%eul_jmin, eul_jm&
    &ax => self%equation_base%eul_jmax, eul_kmin => self%equation_base%eul_kmin, eul_kmax => self%equatio&
    &n_base%eul_kmax, channel_case => self%equation_base%channel_case, nv_aux => self%nv_aux, enable_ibm &
    &=> self%equation_base%enable_ibm)
!
!     
      if (enable_ibm>0) then
        self%equation_base%ibm_force_x = 0._rkind
        self%equation_base%ibm_force_y = 0._rkind
        self%equation_base%ibm_force_z = 0._rkind
      endif
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
        if (enable_ibm>0) then
          call self%base_gpu%bcswap_corner()
          call self%compute_aux()
          !@cuf iercuda=cudaDeviceSynchronize()
          if (istep==1) then
            call self%equation_base%ibm_bc_prepare()
            if (self%equation_base%ibm_vega_moving>0) then
              if (self%equation_base%ibm_trajectory_points>0) self%ibm_vega_y_gpu = self%equation_base%ibm_vega_y
              call self%ibm_inside()
            endif
          endif
          call self%ibm_apply()
          call self%update_ghost() ! needed after application of ibm
        endif
        call self%compute_aux()
        !@cuf iercuda=cudaDeviceSynchronize()
        call self%euler_x(eul_imin, eul_imax)
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
        if (istep == 3) then
          call self%visflx(mode=3) ! 0=lapl, 1=div, 2=reduced, 3=sensor, 4=lapl_nosensor, 5=stag, 6=reduce_nosens
          call self%base_gpu%bcswap_var(self%w_aux_gpu(:,:,:,8:8)) ! ducros
        endif
        if (conservative_viscous==1) then
          call self%visflx(mode=7) ! 0=lapl, 1=div, 2=reduced, 3=sensor, 4=lapl_nosensor, 5=stag, 6=reduce_nosens
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
        if (enable_ibm>0) then
          call self%ibm_compute_force(istep)
        endif
      enddo
!
    endassociate
  endsubroutine rk_sync
!
  subroutine rk_async(self)
    class(equation_singleideal_gpu_object), intent(inout) :: self !< The equation.
    integer :: istep, lmax, iercuda
    real(rkind) :: rhodt, gamdt, alpdt
!   
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv, dt => self%&
    &equation_base%dt, ep_order => self%equation_base%ep_order, weno_scheme => self%equation_base%weno_sc&
    &heme, conservative_viscous => self%equation_base%conservative_viscous, eul_imin => self%equation_bas&
    &e%eul_imin, eul_imax => self%equation_base%eul_imax, eul_jmin => self%equation_base%eul_jmin, eul_jm&
    &ax => self%equation_base%eul_jmax, eul_kmin => self%equation_base%eul_kmin, eul_kmax => self%equatio&
    &n_base%eul_kmax, channel_case => self%equation_base%channel_case, nv_aux => self%nv_aux, enable_ibm &
    &=> self%equation_base%enable_ibm)
!
      if (enable_ibm>0) then
        self%equation_base%ibm_force_x = 0._rkind
        self%equation_base%ibm_force_y = 0._rkind
        self%equation_base%ibm_force_z = 0._rkind
      endif
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
        if (enable_ibm>0) then
          call self%base_gpu%bcswap(steps=[.true.,.false.,.false.])
          call self%base_gpu%bcswap_corner(steps=[.true.,.false.,.false.])
          call self%compute_aux(central=1,ghost=0)
          call self%base_gpu%bcswap(steps=[.false.,.true.,.false.])
          call self%base_gpu%bcswap_corner(steps=[.false.,.true.,.false.])
          call self%base_gpu%bcswap(steps=[.false.,.false.,.true.])
          call self%base_gpu%bcswap_corner(steps=[.false.,.false.,.true.])
          call self%compute_aux(central=0,ghost=1)
          if (istep==1) then
            call self%equation_base%ibm_bc_prepare()
            if (self%equation_base%ibm_vega_moving>0) then
              if (self%equation_base%ibm_trajectory_points>0) self%ibm_vega_y_gpu = self%equation_base%ibm_vega_y
              call self%ibm_inside()
            endif
          endif
          call self%ibm_apply()
          call self%update_ghost(do_swap=0) ! needed after application of ibm
        endif
!
        call self%compute_aux(central=1, ghost=0)
        !@cuf iercuda=cudaDeviceSynchronize()
        call self%base_gpu%bcswap(steps=[.true.,.false.,.false.])
        call self%euler_x(lmax+1,nx-lmax)
        call self%base_gpu%bcswap(steps=[.false.,.true.,.true.])
        call self%compute_aux(central=0, ghost=1)
        if (lmax-eul_imin >= 0) call self%euler_x(eul_imin,lmax)
        if (eul_imax-nx+lmax-1 >= 0) call self%euler_x(nx-lmax+1,eul_imax)
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
        if (istep == 3) then
          call self%visflx(mode=3) ! 0=lapl, 1=div, 2=reduced, 3=sensor, 4=lapl_nosensor, 5=stag
          call self%base_gpu%bcswap_var(self%w_aux_gpu(:,:,:,8:8), steps=[.true.,.false.,.false.]) ! ducros
          if (conservative_viscous==1) then
            call self%visflx(mode=7) ! 0=lapl, 1=div, 2=reduced, 3=sensor, 4=lapl_nosensor, 5=stag, 6=reduce_nosens
          else
            call self%visflx(mode=1)
          endif
          call self%base_gpu%bcswap_var(self%w_aux_gpu(:,:,:,8:8), steps=[.false.,.true.,.true.]) ! ducros
        else
          if (conservative_viscous==1) then
            call self%visflx(mode=7) ! 0=lapl, 1=div, 2=reduced, 3=sensor, 4=lapl_nosensor, 5=stag, 6=reduce_nosens
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
        if (enable_ibm>0) then
          call self%ibm_compute_force(istep)
        endif
      enddo
    endassociate
  endsubroutine rk_async
!
  subroutine euler_x(self, istart, iend)
    class(equation_singleideal_gpu_object), intent(inout) :: self
    integer, intent(in) :: istart, iend
    integer :: lmax, weno_size, iercuda, ierror
    type(dim3) :: grid, tBlock
    integer :: force_zero_flux_min,force_zero_flux_max
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv, nv_aux => s&
    &elf%nv_aux, ep_order => self%equation_base%ep_order, force_zero_flux => self%equation_base%force_zer&
    &o_flux, coeff_deriv1_gpu => self%coeff_deriv1_gpu, dcsidx_gpu => self%base_gpu%dcsidx_gpu, fhat_gpu &
    &=> self%fhat_gpu, w_gpu => self%base_gpu%w_gpu, w_aux_gpu => self%w_aux_gpu, fl_gpu => self%fl_gpu, &
    &sensor_threshold => self%equation_base%sensor_threshold, weno_scheme => self%equation_base%weno_sche&
    &me, weno_version => self%equation_base%weno_version, cp_coeff_gpu => self%cp_coeff_gpu, indx_cp_l =>&
    & self%equation_base%indx_cp_l, indx_cp_r => self%equation_base%indx_cp_r, calorically_perfect => sel&
    &f%equation_base%calorically_perfect, ep_ord_change_gpu => self%ep_ord_change_gpu, nkeep => self%equa&
    &tion_base%nkeep, flux_splitting => self%equation_base%flux_splitting, rgas0 => self%equation_base%rg&
    &as0)
!
      weno_size = 2*weno_scheme
      lmax = ep_order/2 ! max stencil width
      force_zero_flux_min = force_zero_flux(1)
      force_zero_flux_max = force_zero_flux(2)
      tBlock = dim3(EULERWENO_THREADS_X,EULERWENO_THREADS_Y,1)
!
      grid = dim3(ceiling(real(iend-istart+2)/tBlock%x),ceiling(real(ny)/tBlock%y),1)
!
      if (flux_splitting==1) then
!
        call euler_x_fluxes_hybrid_rusanov_kernel<<<grid, tBlock, 0, stream1>>>(nv, nv_aux, nx, ny, &
        &nz, ng, istart, iend, lmax, nkeep, rgas0, w_aux_gpu, coeff_deriv1_gpu, dcsidx_gpu, fhat_gpu, force_z&
        &ero_flux_min, force_zero_flux_max, weno_scheme, weno_version, sensor_threshold, weno_size, cp_coeff_&
        &gpu, indx_cp_l, indx_cp_r, ep_ord_change_gpu, calorically_perfect, tol_iter_nr,self%equation_base%rh&
        &o0,self%equation_base%u0, self%equation_base%t0)
      else
!
        call euler_x_fluxes_hybrid_kernel<<<grid, tBlock, 0, stream1>>>(nv, nv_aux, nx, ny, nz, ng, &
        &istart, iend, lmax, nkeep, rgas0, w_aux_gpu, coeff_deriv1_gpu, dcsidx_gpu, fhat_gpu, force_zero_flux&
        &_min, force_zero_flux_max, weno_scheme, weno_version, sensor_threshold, weno_size, cp_coeff_gpu, ind&
        &x_cp_l, indx_cp_r, ep_ord_change_gpu, calorically_perfect, tol_iter_nr,self%equation_base%rho0,self%&
        &equation_base%u0, self%equation_base%t0)
      endif
!
      call euler_x_update_cuf(nx, ny, nz, ng, nv, istart, iend, fhat_gpu, fl_gpu, dcsidx_gpu, stream1)
!
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
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv, nv_aux => s&
    &elf%nv_aux, ep_order => self%equation_base%ep_order, force_zero_flux => self%equation_base%force_zer&
    &o_flux, coeff_deriv1_gpu => self%coeff_deriv1_gpu, detady_gpu => self%base_gpu%detady_gpu, fhat_gpu &
    &=> self%fhat_gpu, w_aux_gpu => self%w_aux_gpu, fl_gpu => self%fl_gpu, sensor_threshold => self%equat&
    &ion_base%sensor_threshold, weno_scheme => self%equation_base%weno_scheme, weno_version => self%equat&
    &ion_base%weno_version, indx_cp_l => self%equation_base%indx_cp_l, indx_cp_r => self%equation_base%in&
    &dx_cp_r, calorically_perfect => self%equation_base%calorically_perfect, cp_coeff_gpu => self%cp_coef&
    &f_gpu, ep_ord_change_gpu => self%ep_ord_change_gpu, nkeep => self%equation_base%nkeep, flux_splittin&
    &g => self%equation_base%flux_splitting, rgas0 => self%equation_base%rgas0)
      weno_size = 2*weno_scheme
      lmax = ep_order/2 ! max stencil width
      force_zero_flux_min = force_zero_flux(3)
      force_zero_flux_max = force_zero_flux(4)
!
      tBlock = dim3(EULERWENO_THREADS_X,EULERWENO_THREADS_Y,1)
      grid = dim3(ceiling(real(nx)/tBlock%x),ceiling(real(nz)/tBlock%y),1)
      if (flux_splitting==1) then
        call euler_y_hybrid_rusanov_kernel<<<grid, tBlock, 0, stream1>>>(nv, nv_aux, nx, ny, nz, ng,&
        & eul_jmin, eul_jmax, lmax, nkeep, rgas0, w_aux_gpu, fl_gpu, coeff_deriv1_gpu, detady_gpu, fhat_gpu, &
        &force_zero_flux_min, force_zero_flux_max, weno_scheme, weno_version, sensor_threshold, weno_size, cp&
        &_coeff_gpu, indx_cp_l, indx_cp_r, ep_ord_change_gpu, calorically_perfect, tol_iter_nr,self%equation_&
        &base%rho0,self%equation_base%u0, self%equation_base%t0)
      else
        call euler_y_hybrid_kernel<<<grid, tBlock, 0, stream1>>>(nv, nv_aux, nx, ny, nz, ng, eul_jmi&
        &n, eul_jmax, lmax, nkeep, rgas0, w_aux_gpu, fl_gpu, coeff_deriv1_gpu, detady_gpu, fhat_gpu, force_ze&
        &ro_flux_min, force_zero_flux_max, weno_scheme, weno_version, sensor_threshold, weno_size, cp_coeff_g&
        &pu, indx_cp_l, indx_cp_r, ep_ord_change_gpu, calorically_perfect, tol_iter_nr,self%equation_base%rho&
        &0,self%equation_base%u0, self%equation_base%t0)
      endif
      call euler_y_update_cuf(nx, ny, nz, ng, nv, eul_jmin, eul_jmax, fhat_gpu, fl_gpu, detady_gpu, stream1)
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
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv, nv_aux => s&
    &elf%nv_aux, ep_order => self%equation_base%ep_order, force_zero_flux => self%equation_base%force_zer&
    &o_flux, coeff_deriv1_gpu => self%coeff_deriv1_gpu, dzitdz_gpu => self%base_gpu%dzitdz_gpu, fhat_gpu &
    &=> self%fhat_gpu, w_aux_gpu => self%w_aux_gpu, fl_gpu => self%fl_gpu, sensor_threshold => self%equat&
    &ion_base%sensor_threshold, weno_scheme => self%equation_base%weno_scheme, weno_version => self%equat&
    &ion_base%weno_version, cp_coeff_gpu => self%cp_coeff_gpu, indx_cp_l => self%equation_base%indx_cp_l,&
    & indx_cp_r => self%equation_base%indx_cp_r, calorically_perfect => self%equation_base%calorically_pe&
    &rfect, ep_ord_change_gpu => self%ep_ord_change_gpu, nkeep => self%equation_base%nkeep, flux_splittin&
    &g => self%equation_base%flux_splitting, rgas0 => self%equation_base%rgas0, inflow_random_plane => se&
    &lf%equation_base%inflow_random_plane, inflow_random_plane_gpu => self%inflow_random_plane_gpu)
      weno_size = 2*weno_scheme
      lmax = ep_order/2 ! max stencil width
      force_zero_flux_min = force_zero_flux(5)
      force_zero_flux_max = force_zero_flux(6)
      tBlock = dim3(EULERWENO_THREADS_X,EULERWENO_THREADS_Y,1)
      grid = dim3(ceiling(real(nx)/tBlock%x),ceiling(real(ny)/tBlock%y),1)
      if (flux_splitting==1) then
        call euler_z_hybrid_rusanov_kernel<<<grid, tBlock, 0, 0>>>(nv, nv_aux, nx, ny, nz, ng, eul_k&
        &min, eul_kmax, lmax, nkeep, rgas0, w_aux_gpu, fl_gpu, coeff_deriv1_gpu, dzitdz_gpu, fhat_gpu, force_&
        &zero_flux_min, force_zero_flux_max, weno_scheme, weno_version, sensor_threshold, weno_size, cp_coeff&
        &_gpu, indx_cp_l, indx_cp_r, ep_ord_change_gpu, calorically_perfect, tol_iter_nr,self%equation_base%r&
        &ho0,self%equation_base%u0, self%equation_base%t0)
      else
        call euler_z_hybrid_kernel<<<grid, tBlock, 0, 0>>>(nv, nv_aux, nx, ny, nz, ng, eul_kmin, eul&
        &_kmax, lmax, nkeep, rgas0, w_aux_gpu, fl_gpu, coeff_deriv1_gpu, dzitdz_gpu, fhat_gpu, force_zero_flu&
        &x_min, force_zero_flux_max, weno_scheme, weno_version, sensor_threshold, weno_size, cp_coeff_gpu, in&
        &dx_cp_l, indx_cp_r, ep_ord_change_gpu, calorically_perfect, tol_iter_nr,self%equation_base%rho0,self&
        &%equation_base%u0, self%equation_base%t0)
      endif
      call euler_z_update_cuf(nx, ny, nz, ng, nv, eul_kmin, eul_kmax, fhat_gpu, fl_gpu, dzitdz_gpu, 0_cuda_stream_kind)
      if (self%equation_base%recyc) call get_crandom_f(inflow_random_plane(2:self%equation_base%jbl_inflow,1:nz,1:3))
    endassociate
  endsubroutine euler_z
!
  subroutine force_rhs(self)
    class(equation_singleideal_gpu_object), intent(inout) :: self
    real(rkind), dimension(5) :: bulk5, bulk5g
!
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nxmax => self%equation_bas&
    &e%grid%nxmax, nzmax => self%equation_base%grid%nzmax, yn => self%equation_base%field%yn, yn_gpu => s&
    &elf%base_gpu%yn_gpu, w_gpu => self%base_gpu%w_gpu, w_aux_gpu => self%w_aux_gpu, fln_gpu => self%fln_&
    &gpu, volchan => self%equation_base%volchan, fluid_mask_gpu => self%fluid_mask_gpu)
!
      call force_rhs_1_cuf(nx, ny, nz, ng, yn_gpu, fln_gpu, w_gpu, w_aux_gpu, bulk5, fluid_mask_gpu)
!
      call mpi_allreduce(bulk5,bulk5g,5,mpi_prec,mpi_sum,self%equation_base%field%mp_cart,self%mpi_err)
!
      bulk5g = bulk5g/volchan !/nxmax/nzmax/(yn(ny+1)-yn(1))
      self%equation_base%dpdx = self%equation_base%dpdx+bulk5g(2)
      self%equation_base%rhobulk = bulk5g(3)
      self%equation_base%ubulk = bulk5g(4)/self%equation_base%rhobulk
      self%equation_base%tbulk = bulk5g(5)/self%equation_base%rhobulk/self%equation_base%ubulk
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
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nxmax => self%equation_bas&
    &e%grid%nxmax, nzmax => self%equation_base%grid%nzmax, yn => self%equation_base%field%yn, yn_gpu => s&
    &elf%base_gpu%yn_gpu, w_gpu => self%base_gpu%w_gpu, w_aux_gpu => self%w_aux_gpu, fln_gpu => self%fln_&
    &gpu, volchan => self%equation_base%volchan, fluid_mask_gpu => self%fluid_mask_gpu, cv_coeff_gpu => s&
    &elf%cv_coeff_gpu, t0 => self%equation_base%t0, calorically_perfect => self%equation_base%calorically&
    &_perfect, indx_cp_l => self%equation_base%indx_cp_l, indx_cp_r => self%equation_base%indx_cp_r)
!
      if (self%equation_base%theta_wall>=-1._rkind) then
!       
        call force_var_1_cuf(nx, ny, nz, ng, yn_gpu, fln_gpu, w_gpu, w_aux_gpu, bulkt, fluid_mask_gp&
        &u, cv_coeff_gpu, indx_cp_l, indx_cp_r, t0, calorically_perfect, tol_iter_nr)
!       
        call mpi_allreduce(bulkt,bulktg,1,mpi_prec,mpi_sum,self%equation_base%field%mp_cart,self%mpi_err)
!
        bulktg = bulktg/volchan !/nxmax/nzmax/(yn(ny+1)-yn(1))
        bulktg = bulktg/self%equation_base%rhobulk/self%equation_base%ubulk
        tbtarget = self%equation_base%T_bulk_target
        tbdiff = tbtarget-bulktg
!       
        call force_var_2_cuf(nx, ny, nz, ng, w_gpu, w_aux_gpu, tbdiff, fluid_mask_gpu, cv_coeff_gpu,&
        & indx_cp_l, indx_cp_r, t0, calorically_perfect)
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
    associate(nzmax => self%equation_base%grid%nzmax, ng => self%equation_base%grid%ng, nx => self%e&
    &quation_base%field%nx, ny => self%equation_base%field%ny, nz => self%equation_base%field%nz, ncoords&
    & => self%equation_base%field%ncoords, mp_cart => self%equation_base%field%mp_cart, nblocks => self%e&
    &quation_base%field%nblocks, iermpi => self%mpi_err, w_gpu => self%base_gpu%w_gpu, nv => self%equatio&
    &n_base%nv, wbuf1s_gpu => self%base_gpu%wbuf1s_gpu, wbuf2r_gpu => self%base_gpu%wbuf2r_gpu, wbuf1r_gp&
    &u => self%base_gpu%wbuf1r_gpu, wrecyc_gpu => self%wrecyc_gpu, ibrecyc => self%equation_base%ib_recyc&
    &, irecyc => self%equation_base%i_recyc)
!
      kshiftglob = nzmax/2 ! global shift in the spanwise direction (between 0 and nzmax-1)
      n1_start_send = 1
      n1_end_send = nz-mod(kshiftglob,nz)
      n2_start_send = n1_end_send+1
      n2_end_send = nz
      n1_start_recv = 1+mod(kshiftglob,nz)
      n1_end_recv = nz
      n2_start_recv = 1
      n2_end_recv = mod(kshiftglob,nz)
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
    associate(bctags_nr => self%equation_base%bctags_nr, nx => self%nx, ny => self%ny, nz => self%nz&
    &, nv => self%nv, ng => self%ng, w_aux_gpu => self%w_aux_gpu, w_gpu => self%base_gpu%w_gpu, fl_gpu =>&
    & self%fl_gpu, dcsidx_gpu => self%base_gpu%dcsidx_gpu, detady_gpu => self%base_gpu%detady_gpu, dzitdz&
    &_gpu => self%base_gpu%dzitdz_gpu, indx_cp_l => self%equation_base%indx_cp_l, indx_cp_r => self%equat&
    &ion_base%indx_cp_r, calorically_perfect => self%equation_base%calorically_perfect, rgas0 => self%equ&
    &ation_base%rgas0, t0 => self%equation_base%t0, cp_coeff_gpu => self%cp_coeff_gpu, winf_gpu => self%w&
    &inf_gpu)
!
      do ilat=1,6! loop on all sides of the boundary (3D -> 6)
        if(bctags_nr(ilat) > 0) then
          dir = (ilat-1)/2 +1
          start_or_end = mod(ilat-1,2)+1
!         1 - NR
!         2 - relax
!         6 - reflective wall
          if(dir == 1) then
            tBlock = dim3(EULERWENO_THREADS_X,EULERWENO_THREADS_Y,1)
            grid = dim3(ceiling(real(ny)/tBlock%x),ceiling(real(nz)/tBlock%y),1)
            call bc_nr_lat_x_kernel<<<grid, tBlock, 0, 0>>>(start_or_end, bctags_nr(ilat), nx, ny, n&
            &z, ng, nv, w_aux_gpu, w_gpu, fl_gpu, dcsidx_gpu, indx_cp_l, indx_cp_r, cp_coeff_gpu, winf_gpu, calor&
            &ically_perfect,rgas0,t0)
          endif
          if(dir == 2) then
            tBlock = dim3(EULERWENO_THREADS_X,EULERWENO_THREADS_Y,1)
            grid = dim3(ceiling(real(nx)/tBlock%x),ceiling(real(nz)/tBlock%y),1)
            call bc_nr_lat_y_kernel<<<grid, tBlock, 0, 0>>>(start_or_end, bctags_nr(ilat), nx, ny, n&
            &z, ng, nv, w_aux_gpu, w_gpu, fl_gpu, detady_gpu, indx_cp_l, indx_cp_r, cp_coeff_gpu, winf_gpu, calor&
            &ically_perfect,rgas0,t0)
          endif
          if(dir == 3) then
            tBlock = dim3(EULERWENO_THREADS_X,EULERWENO_THREADS_Y,1)
            grid = dim3(ceiling(real(nx)/tBlock%x),ceiling(real(ny)/tBlock%y),1)
            call bc_nr_lat_z_kernel<<<grid, tBlock, 0, 0>>>(start_or_end, bctags_nr(ilat), nx, ny, n&
            &z, ng, nv, w_aux_gpu, w_gpu, fl_gpu, dzitdz_gpu, indx_cp_l, indx_cp_r, cp_coeff_gpu, winf_gpu, calor&
            &ically_perfect,rgas0,t0)
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
        call bcextr_sub_cuf(ilat, self%nx, self%ny, self%nz, self%ng, self%equation_base%p0, self%eq&
        &uation_base%rgas0, self%base_gpu%w_gpu, self%equation_base%indx_cp_l, self%equation_base%indx_cp_r, &
        &self%cv_coeff_gpu, self%equation_base%t0, self%equation_base%calorically_perfect)
      case(5)
        call bcsym_cuf(ilat, self%nx, self%ny, self%nz, self%ng, self%equation_base%T_wall, self%bas&
        &e_gpu%w_gpu, self%w_aux_gpu, self%equation_base%indx_cp_l, self%equation_base%indx_cp_r, self%cv_coe&
        &ff_gpu, self%equation_base%calorically_perfect)
      case(6)
        if (self%equation_base%channel_case) then
          call bcwall_staggered_cuf(ilat, self%nx, self%ny, self%nz, self%ng, self%equation_base%T_w&
          &all, self%base_gpu%w_gpu, self%w_aux_gpu, self%equation_base%indx_cp_l, self%equation_base%indx_cp_r&
          &, self%cv_coeff_gpu, self%equation_base%t0, self%equation_base%rgas0, self%equation_base%calorically&
          &_perfect,tol_iter_nr)
        else
          call bcwall_cuf(ilat, self%nx, self%ny, self%nz, self%ng, self%equation_base%T_wall, self%&
          &base_gpu%w_gpu, self%w_aux_gpu, self%equation_base%indx_cp_l, self%equation_base%indx_cp_r, self%cv_&
          &coeff_gpu, self%equation_base%t0, self%equation_base%rgas0, self%equation_base%calorically_perfect,t&
          &ol_iter_nr)
        endif
      case(7)
        call bcshock_cuf(ilat, self%nx, self%ny, self%nz, self%ng, self%nv, self%base_gpu%w_gpu, sel&
        &f%winf_gpu, self%winf_past_shock_gpu, self%equation_base%xshock_imp, self%equation_base%shock_angle,&
        & self%base_gpu%x_gpu, self%base_gpu%y_gpu, self%equation_base%tanhfacs)
      case(8)
!       adiabatic wall (not implemented here, but implemented for IBM
      case(9)
        call bclam_cuf(ilat, self%nx, self%ny, self%nz, self%ng, self%nv, self%base_gpu%w_gpu, self%&
        &wmean_gpu, self%equation_base%p0, self%equation_base%rgas0, self%equation_base%indx_cp_l, self%equat&
        &ion_base%indx_cp_r, self%cv_coeff_gpu, self%equation_base%t0, self%equation_base%calorically_perfect&
        &)
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
    class(equation_singleideal_gpu_object), intent(inout) :: self !< The equation.
!   Apply recycling-rescaling boundary condition
!   
    integer, intent(in) :: ilat
    integer :: ntot
!   
    if (ilat == 1) then
      associate(nx => self%nx, ny => self%ny, nz => self%nz, nzmax => self%grid%nzmax, ng => self%ng&
      &, nv => self%nv, wrecycav_gpu => self%wrecycav_gpu, wrecyc_gpu => self%wrecyc_gpu, mp_cartz => self%&
      &field%mp_cartz, iermpi => self%mpi_err, p0 => self%equation_base%p0, w_gpu => self%base_gpu%w_gpu, w&
      &mean_gpu => self%wmean_gpu, weta_inflow_gpu => self%weta_inflow_gpu, map_j_inn_gpu => self%map_j_inn&
      &_gpu, map_j_out_gpu => self%map_j_out_gpu, yplus_inflow_gpu => self%yplus_inflow_gpu, eta_inflow_gpu&
      & => self%eta_inflow_gpu, yplus_recyc_gpu => self%yplus_recyc_gpu, eta_recyc_gpu => self%eta_recyc_gp&
      &u, betarecyc => self%equation_base%betarecyc, i_recyc => self%equation_base%i_recyc, indx_cp_l => se&
      &lf%equation_base%indx_cp_l, indx_cp_r => self%equation_base%indx_cp_r, rgas0 => self%equation_base%r&
      &gas0, calorically_perfect => self%equation_base%calorically_perfect, rand_type => self%equation_base&
      &%rand_type, t0 => self%equation_base%t0, u0 => self%equation_base%u0, l0 => self%equation_base%l0, c&
      &v_coeff_gpu => self%cv_coeff_gpu, inflow_random_plane => self%equation_base%inflow_random_plane, inf&
      &low_random_plane_gpu => self%inflow_random_plane_gpu)
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
        call bcrecyc_cuf_3(nx, ny, nz, ng, p0, u0, rgas0, w_gpu, wmean_gpu, wrecyc_gpu, weta_inflow_&
        &gpu, map_j_inn_gpu, map_j_out_gpu, yplus_inflow_gpu, eta_inflow_gpu, yplus_recyc_gpu, eta_recyc_gpu,&
        & betarecyc, inflow_random_plane_gpu, indx_cp_l, indx_cp_r, cv_coeff_gpu, t0, calorically_perfect,ran&
        &d_type)
      endassociate
    endif
!
  end subroutine bcrecyc
!
  subroutine initialize(self, filename)
!   < Initialize the equation.
    class(equation_singleideal_gpu_object), intent(inout) :: self !< The equation.
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
    self%visc_model = self%equation_base%visc_model
    self%mu0 = self%equation_base%mu0
    self%T_ref_dim = self%equation_base%T_ref_dim
    self%sutherland_S = self%equation_base%sutherland_S
    self%powerlaw_vtexp = self%equation_base%powerlaw_vtexp
!
    call self%point_to_field(self%equation_base%field)
    call self%point_to_grid(self%equation_base%grid)
    self%num_iter = self%equation_base%num_iter
    self%time0 = self%equation_base%time0
    self%icyc0 = self%equation_base%icyc0
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
    self%nv_aux = 10
!
    call self%alloc()
!
    self%w_aux_gpu = self%equation_base%w_aux
    self%fluid_mask_gpu = self%equation_base%fluid_mask
    self%ep_ord_change_gpu = self%equation_base%ep_ord_change
    self%winf_gpu = self%equation_base%winf
    self%winf_past_shock_gpu = self%equation_base%winf_past_shock
!
    self%cp_coeff_gpu = self%equation_base%cp_coeff
    self%cv_coeff_gpu = self%equation_base%cv_coeff
!
    if (self%equation_base%num_probe>0) then
      self%probe_coeff_gpu = self%equation_base%probe_coeff
      self%ijk_probe_gpu = self%equation_base%ijk_probe
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
      self%eta_inflow_gpu = self%equation_base%eta_inflow
      self%yplus_recyc_gpu = self%equation_base%yplus_recyc
      self%eta_recyc_gpu = self%equation_base%eta_recyc
      self%map_j_inn_gpu = self%equation_base%map_j_inn
      self%map_j_out_gpu = self%equation_base%map_j_out
      self%weta_inflow_gpu = self%equation_base%weta_inflow
    endif
!
    if (self%equation_base%enable_ibm > 0) call self%ibm_alloc_gpu()
!
    if (self%equation_base%enable_insitu > 0) then
      call self%insitu_alloc_gpu()
    endif
!   
    if (self%equation_base%debug_memory>0) then
      call self%field%check_cpu_mem(description="Initialize-completed")
      call self%base_gpu%check_gpu_mem(description="Initialize-completed")
    endif
!   
!   Allocate field_gpu variables
!   call self%base_gpu%alloc(field=self%field, nv_aux=self%nv_aux)
!
!   ! Use base_gpu as pointee
!   self%field => self%base_gpu%field
!   self%grid => self%base_gpu%field%grid
!
    self%mpi_err = cudaStreamCreate(stream1)
  endsubroutine initialize
!
  subroutine alloc(self)
    class(equation_singleideal_gpu_object), intent(inout) :: self
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv, nv_aux => s&
    &elf%nv_aux, weno_scheme => self%equation_base%weno_scheme,calorically_perfect => self%equation_base%&
    &calorically_perfect, indx_cp_l => self%equation_base%indx_cp_l ,indx_cp_r => self%equation_base%indx&
    &_cp_r)
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
!     allocate(self%w_aux_trans_gpu(1-ng:ny+ng, 1-ng:nx+ng, 1-ng:nz+ng, 8))
!     allocate(self%fhat_trans_gpu(1-ng:ny+ng, 1-ng:nx+ng, 1-ng:nz+ng, nv))
!     allocate(self%fl_trans_gpu(1:ny, 1:nx, 1:nz, nv))
!     allocate(self%gplus_x_gpu (0:nx,ny,nv,2*weno_scheme))
!     allocate(self%gminus_x_gpu(0:nx,ny,nv,2*weno_scheme))
!     allocate(self%gplus_y_gpu (nx,nz,nv,2*weno_scheme))
!     allocate(self%gminus_y_gpu(nx,nz,nv,2*weno_scheme))
!     allocate(self%gplus_z_gpu (nx,ny,nv,2*weno_scheme))
!     allocate(self%gminus_z_gpu(nx,ny,nv,2*weno_scheme))
      allocate(self%fluid_mask_gpu(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng))
      allocate(self%ep_ord_change_gpu(0:nx, 0:ny, 0:nz, 1:3))
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
    class(equation_singleideal_gpu_object), intent(inout) :: self !< The equation.
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv, dt => self%&
    &equation_base%dt, fln_gpu => self%fln_gpu, residual_rhou => self%equation_base%residual_rhou, fluid_&
    &mask_gpu => self%fluid_mask_gpu, vmax => self%equation_base%vmax, w_gpu => self%base_gpu%w_gpu)
      call compute_residual_cuf(nx, ny, nz, ng, nv, fln_gpu, dt, residual_rhou, fluid_mask_gpu)
      call mpi_allreduce(MPI_IN_PLACE,residual_rhou,1,mpi_prec,mpi_sum,self%equation_base%field%mp_cart,self%mpi_err)
      residual_rhou = residual_rhou / (real(nx,rkind)*real(ny,rkind)*real(nz,rkind)*real(self%nprocs,rkind))
      residual_rhou = sqrt(residual_rhou)
!     call compute_vmax_cuf(nx, ny, nz, ng, nv, w_gpu, vmax, fluid_mask_gpu)
    endassociate
  endsubroutine compute_residual
!
  subroutine compute_dt(self)
!   < Initialize the equation.
    class(equation_singleideal_gpu_object), intent(inout) :: self !< The equation.
    real(rkind) :: dt_min, dtinv_max
    real(rkind) :: dtxi_max, dtyi_max, dtzi_max, dtxv_max, dtyv_max, dtzv_max, dtxk_max, dtyk_max, dtzk_max
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, dt => self%equation_base%d&
    &t, CFL => self%equation_base%CFL, Prandtl => self%equation_base%Prandtl, visc_model => self%visc_mod&
    &el, mu0 => self%mu0, T_ref_dim => self%T_ref_dim, powerlaw_vtexp => self%powerlaw_vtexp, sutherland_&
    &S => self%sutherland_S, w_gpu => self%base_gpu%w_gpu, w_aux_gpu => self%w_aux_gpu, dcsidx_gpu => sel&
    &f%base_gpu%dcsidx_gpu, detady_gpu => self%base_gpu%detady_gpu, dzitdz_gpu => self%base_gpu%dzitdz_gp&
    &u, dcsidxs_gpu => self%base_gpu%dcsidxs_gpu, detadys_gpu => self%base_gpu%detadys_gpu, dzitdzs_gpu =&
    &> self%base_gpu%dzitdzs_gpu, calorically_perfect => self%equation_base%calorically_perfect, indx_cp_&
    &l => self%equation_base%indx_cp_l, indx_cp_r => self%equation_base%indx_cp_r, cp_coeff_gpu => self%c&
    &p_coeff_gpu, fluid_mask_gpu => self%fluid_mask_gpu, rgas0 => self%equation_base%rgas0, t0 => self%eq&
    &uation_base%t0)
      if (CFL < 0) then
        dt = -CFL
      else
        call compute_dt_cuf(nx, ny, nz, ng, rgas0, Prandtl, dcsidx_gpu, detady_gpu, dzitdz_gpu, dcsi&
        &dxs_gpu, detadys_gpu, dzitdzs_gpu, w_gpu, w_aux_gpu, dtxi_max, dtyi_max, dtzi_max, dtxv_max, dtyv_ma&
        &x, dtzv_max, dtxk_max, dtyk_max, dtzk_max, indx_cp_l, indx_cp_r, cp_coeff_gpu,fluid_mask_gpu,caloric&
        &ally_perfect,t0)
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
    class(equation_singleideal_gpu_object), intent(inout) :: self !< The equation.
    character(*) , intent(in) :: filename !< Input file name.
    real(rkind) :: timing(1:2) !< Tic toc timing.
    real(rkind) :: timing_step(1:2) !< Tic toc timing.
    integer :: icyc_loop, iercuda
!
    call self%initialize(filename=filename)
!
    associate(icyc0 => self%equation_base%icyc0, icyc => self%equation_base%icyc, time => self%equat&
    &ion_base%time, iter_dt_recompute => self%equation_base%iter_dt_recompute, residual_rhou => self%equa&
    &tion_base%residual_rhou, dpdx => self%equation_base%dpdx, rhobulk => self%equation_base%rhobulk, ubu&
    &lk => self%equation_base%ubulk, tbulk => self%equation_base%tbulk, nx => self%nx, ny => self%ny, nz &
    &=> self%nz, ng => self%ng, visc_model => self%visc_model, mu0 => self%mu0, T_ref_dim => self%T_ref_d&
    &im, sutherland_S => self%sutherland_S, powerlaw_vtexp => self%powerlaw_vtexp, cv_coeff_gpu => self%c&
    &v_coeff_gpu, calorically_perfect => self%equation_base%calorically_perfect, indx_cp_l => self%equati&
    &on_base%indx_cp_l, indx_cp_r => self%equation_base%indx_cp_r, enable_ibm => self%equation_base%enabl&
    &e_ibm, mode_async => self%equation_base%mode_async, time_from_last_rst => self%equation_base%time_fr&
    &om_last_rst, time_from_last_write => self%equation_base%time_from_last_write, time_from_last_stat =>&
    & self%equation_base%time_from_last_stat, time_from_last_slice => self%equation_base%time_from_last_s&
    &lice, time_from_last_insitu => self%equation_base%time_from_last_insitu, time_is_freezed => self%equ&
    &ation_base%time_is_freezed)
!
      call self%update_ghost()
      if (enable_ibm>0) then
        call self%base_gpu%bcswap_corner()
        call self%compute_aux()
        !@cuf iercuda=cudaDeviceSynchronize()
        call self%equation_base%ibm_bc_prepare()
        if (self%equation_base%ibm_vega_moving>0) then
          if (self%equation_base%ibm_trajectory_points>0) self%ibm_vega_y_gpu = self%equation_base%ibm_vega_y
          call self%ibm_inside()
        endif
        call self%ibm_apply()
        call self%update_ghost() ! needed after application of ibm
      endif
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
          call self%field%write_plot3d(mach=self%equation_base%Mach, reynolds=self%equation_base%Rey&
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
!
      call MPI_BARRIER(MPI_COMM_WORLD, self%error) ; timing(1) = MPI_Wtime()
!
!     Reinit random to avoid missing reproducibility with inflow random numbers
      if(self%equation_base%rand_type == 0) then
        call init_crandom_f(0,reproducible=.true.)
        if (self%masterproc) write(*,*) 'Random numbers disabled'
      elseif(self%equation_base%rand_type < 0) then
        call init_crandom_f(self%myrank,reproducible=.false.)
        if (self%masterproc) write(*,*) 'Random numbers NOT reproducible'
      else
        call init_crandom_f(self%myrank,reproducible=.true.)
        if (self%masterproc) write(*,*) 'Random numbers reproducible'
      endif
!     
      icyc_loop = icyc
      integration: do
!
        icyc_loop = icyc_loop + 1
        time_is_freezed = self%equation_base%time_is_freezed_fun()
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
            if(mode_async == 0) call self%rk_sync()
            if(mode_async == 1) call self%rk_async()
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
          write(*,'(A, F18.10)') 'averaged timing: ', (timing(2) - timing(1))/(self%equation_base%icyc-self%equation_base%icyc0)
        endif
      endif
!
      call self%base_gpu%copy_gpu_cpu()
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
!
    endassociate
  endsubroutine run
!
  subroutine print_progress(self)
    class(equation_singleideal_gpu_object), intent(inout) :: self
    character(6) :: pos_io
    associate(icyc => self%equation_base%icyc, time => self%equation_base%time, dt => self%equation_&
    &base%dt, residual_rhou => self%equation_base%residual_rhou, dpdx => self%equation_base%dpdx, vmax =>&
    & self%equation_base%vmax, rhobulk => self%equation_base%rhobulk, ubulk => self%equation_base%ubulk, &
    &tbulk => self%equation_base%tbulk, ibm_force_x => self%equation_base%ibm_force_x, ibm_force_y => sel&
    &f%equation_base%ibm_force_y, ibm_force_z => self%equation_base%ibm_force_z)
!     
      residual_rhou = residual_rhou/(self%equation_base%u0**2)/self%equation_base%rho0*self%equation_base%l0
      dpdx = dpdx /(self%equation_base%u0**2)/self%equation_base%rho0*self%equation_base%l0
      pos_io = 'append'
      if (self%equation_base%icyc==1) pos_io = 'rewind'
      if (self%masterproc) then
        open(unit=15,file='progress.out',position=pos_io)
        if (self%equation_base%channel_case) then
          write(* ,'(1I10,20ES20.10)') icyc, dt, time, residual_rhou, -dpdx/dt, rhobulk, ubulk, tbulk
          write(15,'(1I10,20ES20.10)') icyc, dt, time, residual_rhou, -dpdx/dt, rhobulk, ubulk, tbulk
        else
          if (self%equation_base%enable_ibm>0) then
            write(* ,'(1I10,20ES20.10)') icyc, dt, time, residual_rhou, ibm_force_x, ibm_force_y, ibm_force_z
            write(15,'(1I10,20ES20.10)') icyc, dt, time, residual_rhou, ibm_force_x, ibm_force_y, ibm_force_z
          else
!           write(* ,'(1I10,20ES20.10)') icyc, dt, time, residual_rhou,vmax
            write(* ,'(1I10,20ES20.10)') icyc, dt, time, residual_rhou
            write(15,'(1I10,20ES20.10)') icyc, dt, time, residual_rhou
          endif
        endif
        close(15)
      endif
    endassociate
  endsubroutine print_progress
!
  subroutine ibm_inside(self)
    class(equation_singleideal_gpu_object), intent(inout) :: self
!   
    associate(ng => self%ng, nx => self%nx, ny => self%ny, nz => self%nz, x_gpu => self%base_gpu%x_g&
    &pu, y_gpu => self%base_gpu%y_gpu, z_gpu => self%base_gpu%z_gpu, ibm_inside_moving_gpu => self%ibm_in&
    &side_moving_gpu, ibm_vega_y_gpu => self%ibm_vega_y_gpu, ibm_vega_r_gpu => self%ibm_vega_r_gpu, ibm_v&
    &ega_ny => self%equation_base%ibm_vega_ny, ibm_vega_dy => self%equation_base%ibm_vega_dy, ep_ord_chan&
    &ge_gpu => self%ep_ord_change_gpu, ibm_order_reduce => self%equation_base%ibm_order_reduce)
!
      call ibm_inside_moving_cuf(nx,ny,nz,ng,ibm_inside_moving_gpu, x_gpu,y_gpu,z_gpu,ibm_vega_ny, i&
      &bm_vega_dy,ibm_vega_y_gpu,ibm_vega_r_gpu, ep_ord_change_gpu,ibm_order_reduce)
!
    endassociate
!   
  end subroutine ibm_inside
!
  subroutine ibm_apply(self)
    class(equation_singleideal_gpu_object), intent(inout) :: self
!   
    associate(ibm_num_interface => self%equation_base%ibm_num_interface, ibm_ijk_refl_gpu => self%ib&
    &m_ijk_refl_gpu, ibm_refl_type_gpu => self%ibm_refl_type_gpu, ibm_coeff_d_gpu => self%ibm_coeff_d_gpu&
    &, ibm_coeff_n_gpu => self%ibm_coeff_n_gpu, ibm_is_interface_node_gpu => self%ibm_is_interface_node_g&
    &pu, ibm_bc_gpu => self%ibm_bc_gpu, w_gpu => self%base_gpu%w_gpu, w_aux_gpu => self%w_aux_gpu, cv_coe&
    &ff_gpu => self%cv_coeff_gpu, indx_cp_l => self%equation_base%indx_cp_l, ibm_w_refl_gpu => self%ibm_w&
    &_refl_gpu, ng => self%ng, nx => self%nx, ny => self%ny, nz => self%nz, indx_cp_r => self%equation_ba&
    &se%indx_cp_r, calorically_perfect => self%equation_base%calorically_perfect, t0 => self%equation_bas&
    &e%t0, ibm_ijk_interface_gpu => self%ibm_ijk_interface_gpu, ibm_nxyz_interface_gpu => self%ibm_nxyz_i&
    &nterface_gpu, ibm_aero_rad => self%equation_base%ibm_aero_rad, ibm_aero_pp => self%equation_base%ibm&
    &_aero_pp, ibm_aero_tt => self%equation_base%ibm_aero_tt, ibm_aero_modvel => self%equation_base%ibm_a&
    &ero_modvel, x_gpu => self%base_gpu%x_gpu, y_gpu => self%base_gpu%y_gpu, z_gpu => self%base_gpu%z_gpu&
    &, ibm_twall => self%equation_base%T_wall, ibm_bc_relax_factor => self%equation_base%ibm_bc_relax_fac&
    &tor, ibm_parbc_gpu => self%ibm_parbc_gpu, ibm_inside_moving_gpu => self%ibm_inside_moving_gpu, rgas0&
    & => self%equation_base%rgas0, fluid_mask_gpu => self%fluid_mask_gpu, ibm_num_bc => self%equation_bas&
    &e%ibm_num_bc, ibm_vega_vel => self%equation_base%ibm_vega_vel)
!
      if (self%equation_base%ibm_vega_moving>0) then
        call ibm_forcing1_cuf(nx,ny,nz,ng,indx_cp_l,indx_cp_r,w_gpu,ibm_inside_moving_gpu, cv_coeff_&
        &gpu,ibm_aero_rad,ibm_aero_pp,ibm_aero_tt, ibm_vega_vel, ibm_aero_modvel,x_gpu,y_gpu,z_gpu,ibm_bc_rel&
        &ax_factor,t0,calorically_perfect,rgas0)
      endif
      call ibm_interpolation_cuf(ibm_num_interface,nx,ny,nz,ng,indx_cp_l,indx_cp_r,ibm_ijk_refl_gpu,&
      &ibm_refl_type_gpu,w_gpu,w_aux_gpu,ibm_is_interface_node_gpu, ibm_coeff_d_gpu,ibm_coeff_n_gpu,ibm_bc_&
      &gpu,cv_coeff_gpu,ibm_w_refl_gpu,ibm_twall,ibm_nxyz_interface_gpu,calorically_perfect,rgas0)
      call ibm_forcing_cuf(ibm_num_interface,nx,ny,nz,ng,indx_cp_l,indx_cp_r,ibm_ijk_interface_gpu,w&
      &_gpu,w_aux_gpu,ibm_bc_gpu, cv_coeff_gpu,ibm_w_refl_gpu,ibm_nxyz_interface_gpu,ibm_aero_rad,ibm_aero_&
      &pp,ibm_aero_tt,ibm_aero_modvel, x_gpu,y_gpu,z_gpu,ibm_twall,ibm_bc_relax_factor,t0,calorically_perfe&
      &ct,rgas0,ibm_parbc_gpu, ibm_num_bc)
    endassociate
!   
  end subroutine ibm_apply
!
  subroutine ibm_compute_force(self,istep)
    class(equation_singleideal_gpu_object), intent(inout) :: self
!   
    integer, intent(in) :: istep
    real(rkind) :: ibm_force_x_s, ibm_force_y_s, ibm_force_z_s
!   
    associate(ibm_num_interface => self%equation_base%ibm_num_interface, ibm_bc_gpu => self%ibm_bc_g&
    &pu, w_gpu => self%base_gpu%w_gpu, w_aux_gpu => self%w_aux_gpu, ng => self%ng, ibm_ijk_interface_gpu &
    &=> self%ibm_ijk_interface_gpu, ibm_force_x => self%equation_base%ibm_force_x, ibm_force_y => self%eq&
    &uation_base%ibm_force_y, ibm_force_z => self%equation_base%ibm_force_z, dcsidx_gpu => self%base_gpu%&
    &dcsidx_gpu, detady_gpu => self%base_gpu%detady_gpu, dzitdz_gpu => self%base_gpu%dzitdz_gpu, fluid_ma&
    &sk_gpu => self%fluid_mask_gpu, fln_gpu => self%fln_gpu, nx => self%nx, ny => self%ny, nz => self%nz,&
    & iermpi => self%mpi_err, dt => self%equation_base%dt)
!
      call ibm_compute_force_cuf(ibm_num_interface,nx,ny,nz,ng,ibm_ijk_interface_gpu,fln_gpu,w_gpu,w&
      &_aux_gpu,ibm_bc_gpu,dcsidx_gpu,detady_gpu,dzitdz_gpu,ibm_force_x_s,ibm_force_y_s,ibm_force_z_s,fluid&
      &_mask_gpu)
      ibm_force_x = ibm_force_x + ibm_force_x_s
      ibm_force_y = ibm_force_y + ibm_force_y_s
      ibm_force_z = ibm_force_z + ibm_force_z_s
!     
      if (istep==self%equation_base%nrk) then
        call mpi_allreduce(mpi_in_place,ibm_force_x,1,mpi_prec,mpi_sum,mpi_comm_world,iermpi)
        call mpi_allreduce(mpi_in_place,ibm_force_y,1,mpi_prec,mpi_sum,mpi_comm_world,iermpi)
        call mpi_allreduce(mpi_in_place,ibm_force_z,1,mpi_prec,mpi_sum,mpi_comm_world,iermpi)
        ibm_force_x = -ibm_force_x/dt
        ibm_force_y = -ibm_force_y/dt
        ibm_force_z = -ibm_force_z/dt
      endif
!     
    endassociate
!   
  end subroutine ibm_compute_force
!
  subroutine insitu_coprocess(self)
    class(equation_singleideal_gpu_object), intent(inout) :: self
    integer :: l,i,j,k,ll
    call self%update_ghost()
!   call self%base_gpu%bcswap()
    call self%base_gpu%bcswap_corner()
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, visc_model => self%visc_mo&
    &del, mu0 => self%mu0, T_ref_dim => self%T_ref_dim, sutherland_S => self%sutherland_S, powerlaw_vtexp&
    & => self%powerlaw_vtexp, w_aux_gpu => self%w_aux_gpu, cv_coeff_gpu => self%cv_coeff_gpu, indx_cp_l =&
    &> self%equation_base%indx_cp_l, indx_cp_r => self%equation_base%indx_cp_r, calorically_perfect => se&
    &lf%equation_base%calorically_perfect)
      call self%compute_aux()
    endassociate
!   filling psi (optionally we can fill ghosts and corners only of psi_pv_managed)
    if (self%equation_base%npsi > 0) then
      call self%insitu_compute_psi()
      do l=1,self%equation_base%npsi
        call self%base_gpu%bcswap_var(self%psi_gpu(:,:,:,l:l ))
      enddo
      do l=1,self%equation_base%npsi
        call self%base_gpu%bcswap_corner_var(self%psi_gpu(:,:,:,l:l ))
      enddo
    endif
    associate( nxsl_ins => self%equation_base%nxsl_ins, nxel_ins => self%equation_base%nxel_ins, nys&
    &l_ins => self%equation_base%nysl_ins, nyel_ins => self%equation_base%nyel_ins, nzsl_ins => self%equa&
    &tion_base%nzsl_ins, nzel_ins => self%equation_base%nzel_ins, npsi => self%equation_base%npsi, npsi_p&
    &v => self%equation_base%npsi_pv, nv_aux => self%equation_base%nv_aux, n_aux_list => self%equation_ba&
    &se%n_aux_list, psi_gpu => self%psi_gpu, n_add_list => self%equation_base%n_add_list, w_aux_gpu => se&
    &lf%w_aux_gpu, aux_list_gpu => self%aux_list_gpu, add_list_gpu => self%add_list_gpu, nx => self%field&
    &%nx, ny => self%field%ny, nz => self%field%nz, icyc => self%equation_base%icyc, time => self%equatio&
    &n_base%time, i_insitu => self%equation_base%i_insitu, flag => self%equation_base%insitu_flag, nrank &
    &=> self%myrank, ng => self%grid%ng, aux_list_name => self%equation_base%aux_list_name, add_list_name&
    & => self%equation_base%add_list_name, time_insitu => self%equation_base%time_insitu)
!
      self%equation_base%w_aux = w_aux_gpu
      self%equation_base%psi = psi_gpu
!
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
!
!     NOMANAGEDcall copy_to_psi_pv_managed_cuf(nxsl_ins,nxel_ins,nysl_ins,nyel_ins,nzsl_ins,nzel_ins
!, !NOMANAGED ng,nx,ny,nz,npsi, npsi_pv, nv_aux, n_aux_list, !NOMANAGED psi_gpu,self%psi_pv_managed,w
!_aux_gpu,aux_list_gpu)
!     call requestdatadescription(icyc,time,flag)
      if(self%equation_base%insitu_platform == "catalyst-v1") then
        call requestdatadescription(i_insitu,time_insitu,flag)
        if (flag.ne.0) then
          call needtocreategrid(flag)
          if (flag.ne.0) then
            do l=1,n_aux_list
              call addfieldtostructured("3d_struct"//c_null_char,self%equation_base%psi_pv(:,:,:,l),&
              & trim(adjustl(aux_list_name(l)))//c_null_char,nrank)
!             NOMANAGEDcall addfieldtostructured("3d_struct"//c_null_char,self%psi_pv_managed(:,:,:,
!l), !NOMANAGED trim(adjustl(aux_list_name(l)))//c_null_char,nrank)
            enddo
            do l=1,n_add_list
              call addfieldtostructured("3d_struct"//c_null_char,self%equation_base%psi_pv(:,:,:,n_a&
              &ux_list+l), trim(adjustl(add_list_name(l)))//c_null_char,nrank)
!             NOMANAGEDcall addfieldtostructured("3d_struct"//c_null_char,self%psi_pv_managed(:,:,:,
!n_aux_list+l), !NOMANAGED trim(adjustl(add_list_name(l)))//c_null_char,nrank)
            enddo
            call coprocess()
          end if
        end if
      elseif(self%equation_base%insitu_platform == "catalyst-v2") then
        call self%insitu_do_catalyst_execute()
      endif
    endassociate
!
  end subroutine insitu_coprocess
!
  subroutine insitu_do_catalyst_execute(self)
    class(equation_singleideal_gpu_object), intent(inout) :: self
!
    integer :: cycle
    real(real64) :: time
    type(C_PTR) :: catalyst_exec_params, mesh, info
    type(C_PTR) :: xt, yt, zt
    type(C_PTR), dimension(:), allocatable :: vx
    integer(kind(catalyst_status)) :: code
    integer :: exit_code
    integer :: l
!
    catalyst_exec_params = catalyst_conduit_node_create()
    call catalyst_conduit_node_set_path_int64(catalyst_exec_params, "catalyst/state/timestep", int(s&
    &elf%equation_base%i_insitu,int64))
!   one can also use "catalyst/cycle" for the same purpose.
!   conduit_node_set_path_int64(catalyst_exec_params, "catalyst/state/cycle", cycle)
    call catalyst_conduit_node_set_path_float64(catalyst_exec_params, "catalyst/state/time", self%equation_base%time_insitu)
!
!   the data must be provided on a named channel. the name is determined by the
!   simulation. for this one, we're calling it "grid".
!
!   declare the type of the channel; we're using Conduit Mesh Blueprint
!   to describe the mesh and fields.
    call catalyst_conduit_node_set_path_char8_str(catalyst_exec_params, "catalyst/channels/grid/type", "mesh")
!
!   now, create the mesh.
    mesh = catalyst_conduit_node_create()
!
    associate(points_x => self%equation_base%points_x, points_y => self%equation_base%points_y, poin&
    &ts_z => self%equation_base%points_z, n_points_x => self%equation_base%n_points_x, n_points_y => self&
    &%equation_base%n_points_y, n_points_z => self%equation_base%n_points_z, n_points => self%equation_ba&
    &se%n_points, n_aux_list => self%equation_base%n_aux_list, n_add_list => self%equation_base%n_add_lis&
    &t, aux_list_name => self%equation_base%aux_list_name, add_list_name => self%equation_base%add_list_n&
    &ame)
!
      allocate(vx(n_aux_list+n_add_list))
!
!     *************** STUCTURED START
!     add coordsets
      call catalyst_conduit_node_set_path_char8_str(mesh,"coordsets/coords/type","explicit")
      xt = catalyst_conduit_node_create()
      yt = catalyst_conduit_node_create()
      zt = catalyst_conduit_node_create()
      call catalyst_conduit_node_set_external_float64_ptr(xt, points_x, n_points)
      call catalyst_conduit_node_set_external_float64_ptr(yt, points_y, n_points)
      call catalyst_conduit_node_set_external_float64_ptr(zt, points_z, n_points)
      call catalyst_conduit_node_set_path_external_node(mesh, "coordsets/coords/values/x", xt)
      call catalyst_conduit_node_set_path_external_node(mesh, "coordsets/coords/values/y", yt)
      call catalyst_conduit_node_set_path_external_node(mesh, "coordsets/coords/values/z", zt)
!     call c_catalyst_conduit_node_set_path_external_float64_ptr(mesh, "coordsets/coords/values/x"//
!C_NULL_CHAR, ! points_x, n_points)
!     call c_catalyst_conduit_node_set_path_external_float64_ptr(mesh, "coordsets/coords/values/y"//
!C_NULL_CHAR, ! points_y, n_points)
!     call c_catalyst_conduit_node_set_path_external_float64_ptr(mesh, "coordsets/coords/values/z"//
!C_NULL_CHAR, ! points_z, n_points)
      call catalyst_conduit_node_set_path_char8_str(mesh, "topologies/mesh/type", "structured")
      call catalyst_conduit_node_set_path_char8_str(mesh, "topologies/mesh/coordset", "coords")
      call catalyst_conduit_node_set_path_int32(mesh, "topologies/mesh/elements/dims/i", n_points_x-1)
      call catalyst_conduit_node_set_path_int32(mesh, "topologies/mesh/elements/dims/j", n_points_y-1)
      call catalyst_conduit_node_set_path_int32(mesh, "topologies/mesh/elements/dims/k", n_points_z-1)
!
      do l=1,n_aux_list
        vx(l) = catalyst_conduit_node_create()
        call catalyst_conduit_node_set_external_float64_ptr(vx(l), self%equation_base%psi_pv(:,:,:,l), n_points)
        call catalyst_conduit_node_set_path_char8_str(mesh, "fields/"//trim(adjustl(aux_list_name(l)))//"/association", "vertex")
        call catalyst_conduit_node_set_path_char8_str(mesh, "fields/"//trim(adjustl(aux_list_name(l)))//"/topology", "mesh")
        call catalyst_conduit_node_set_path_char8_str(mesh, "fields/"//trim(adjustl(aux_list_name(l)&
        &))//"/volume_dependent", "false")
        call catalyst_conduit_node_set_path_external_node(mesh, "fields/"//trim(adjustl(aux_list_name(l)))//"/values", vx(l))
      enddo
      do l=1,n_add_list
        vx(n_aux_list+l) = catalyst_conduit_node_create()
        call catalyst_conduit_node_set_external_float64_ptr(vx(n_aux_list+l), self%equation_base%psi&
        &_pv(:,:,:,n_aux_list+l), n_points)
        call catalyst_conduit_node_set_path_char8_str(mesh, "fields/"//trim(adjustl(add_list_name(l)))//"/association", "vertex")
        call catalyst_conduit_node_set_path_char8_str(mesh, "fields/"//trim(adjustl(add_list_name(l)))//"/topology", "mesh")
        call catalyst_conduit_node_set_path_char8_str(mesh, "fields/"//trim(adjustl(add_list_name(l)&
        &))//"/volume_dependent", "false")
        call catalyst_conduit_node_set_path_external_node(mesh, "fields/"//trim(adjustl(add_list_nam&
        &e(l)))//"/values", vx(n_aux_list+l))
      enddo
!
!     *************** STUCTURED END
!
      call catalyst_conduit_node_set_path_external_node(catalyst_exec_params, "catalyst/channels/grid/data", mesh)
!
#if 1
!     print for debugging purposes, if needed
      call catalyst_conduit_node_print(catalyst_exec_params)
!
!     print information with details about memory allocation
!     info = catalyst_conduit_node_create()
!     call catalyst_conduit_node_info(catalyst_exec_params, info)
!     call catalyst_conduit_node_print(info)
!     call catalyst_conduit_node_destroy(info)
#endif
!
      code = c_catalyst_execute(catalyst_exec_params)
      if (code /= catalyst_status_ok) then
        write (error_unit, *) "failed to call `execute`:", code
        exit_code = 1
      end if
      call catalyst_conduit_node_destroy(catalyst_exec_params)
      call catalyst_conduit_node_destroy(mesh)
      call catalyst_conduit_node_destroy(xt)
      call catalyst_conduit_node_destroy(yt)
      call catalyst_conduit_node_destroy(zt)
      do l=1,n_aux_list+n_add_list
        call catalyst_conduit_node_destroy(vx(l))
      enddo
!
      deallocate(vx)
!
    endassociate
!
  endsubroutine insitu_do_catalyst_execute
!
  subroutine insitu_compute_psi(self)
    class(equation_singleideal_gpu_object), intent(inout) :: self
    integer :: l
    type(dim3) :: grid, tBlock
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv, nv_aux => s&
    &elf%nv_aux, dt => self%equation_base%dt, visc_order => self%equation_base%visc_order, Prandtl => sel&
    &f%equation_base%Prandtl, visc_model => self%visc_model, mu0 => self%mu0, u0 => self%equation_base%u0&
    &, l0 => self%equation_base%l0, T_ref_dim => self%T_ref_dim, sutherland_S => self%sutherland_S, power&
    &law_vtexp => self%powerlaw_vtexp, coeff_deriv1_gpu => self%coeff_deriv1_gpu, coeff_deriv2_gpu => sel&
    &f%coeff_deriv2_gpu, fhat_trans_gpu => self%fhat_trans_gpu, fl_trans_gpu => self%fl_trans_gpu, fl_gpu&
    & => self%fl_gpu, w_aux_gpu => self%w_aux_gpu, w_aux_trans_gpu => self%w_aux_trans_gpu, dcsidx_gpu =>&
    & self%base_gpu%dcsidx_gpu, detady_gpu => self%base_gpu%detady_gpu, dzitdz_gpu => self%base_gpu%dzitd&
    &z_gpu, dcsidxs_gpu => self%base_gpu%dcsidxs_gpu, detadys_gpu => self%base_gpu%detadys_gpu, dzitdzs_g&
    &pu => self%base_gpu%dzitdzs_gpu, dcsidx2_gpu => self%base_gpu%dcsidx2_gpu, detady2_gpu => self%base_&
    &gpu%detady2_gpu, dzitdz2_gpu => self%base_gpu%dzitdz2_gpu, eul_imin => self%equation_base%eul_imin, &
    &eul_imax => self%equation_base%eul_imax, eul_jmin => self%equation_base%eul_jmin, eul_jmax => self%e&
    &quation_base%eul_jmax, eul_kmin => self%equation_base%eul_kmin, eul_kmax => self%equation_base%eul_k&
    &max, cv_coeff_gpu => self%cv_coeff_gpu, indx_cp_l => self%equation_base%indx_cp_l, indx_cp_r => self&
    &%equation_base%indx_cp_r, calorically_perfect => self%equation_base%calorically_perfect, ! ibm => se&
    &lf%equation_base%ibm, npsi => self%equation_base%npsi, psi_gpu => self%psi_gpu, add_list => self%equ&
    &ation_base%add_list, x_gpu => self%base_gpu%x_gpu )
      do l=1,npsi
        if (add_list(l) == 1) then
          call insitu_div_cuf(nx, ny, nz, ng, visc_order, npsi,l, w_aux_gpu, coeff_deriv1_gpu, dcsid&
          &x_gpu, detady_gpu, dzitdz_gpu, psi_gpu )
        endif
        if (add_list(l) == 2) then
          call insitu_omega_cuf(nx, ny, nz, ng, visc_order, npsi,l, w_aux_gpu, coeff_deriv1_gpu, dcs&
          &idx_gpu, detady_gpu, dzitdz_gpu, psi_gpu )
        endif
        if (add_list(l) == 3) then
          call insitu_ducros_cuf(nx, ny, nz, ng, visc_order, npsi,l, u0,l0,w_aux_gpu, coeff_deriv1_g&
          &pu, dcsidx_gpu, detady_gpu, dzitdz_gpu, psi_gpu )
        endif
        if (add_list(l) == 4) then
          tBlock = dim3(EULERWENO_THREADS_X,EULERWENO_THREADS_Y,1)
          grid = dim3(ceiling(real(nx)/tBlock%x),ceiling(real(nz)/tBlock%y),1)
          call insitu_swirling_kernel<<<grid, tBlock>>>(nv, nx, ny, nz, visc_order, ng, npsi,l,dcsid&
          &x_gpu, detady_gpu, dzitdz_gpu, w_aux_gpu, coeff_deriv1_gpu, psi_gpu, u0, x_gpu )
        endif
      enddo
    endassociate
!
  end subroutine insitu_compute_psi
!
  subroutine manage_output(self)
    class(equation_singleideal_gpu_object), intent(inout) :: self
!   
    integer :: i,j,k,ii,jj,kk,l,n
    integer :: isize, jsize, ksize
    character(3) :: chx, chy, chz
    logical :: sliceyz_exist, slicexz_exist, slicexy_exist, probe_exist
!   
    associate(time_from_last_rst => self%equation_base%time_from_last_rst, time_from_last_write => s&
    &elf%equation_base%time_from_last_write, time_from_last_stat => self%equation_base%time_from_last_sta&
    &t, time_from_last_slice => self%equation_base%time_from_last_slice, time_from_last_insitu => self%eq&
    &uation_base%time_from_last_insitu, icyc0 => self%equation_base%icyc0, icyc => self%equation_base%icy&
    &c, time => self%equation_base%time, w_aux_gpu => self%w_aux_gpu, ijk_probe_gpu => self%ijk_probe_gpu&
    &, w_aux_probe_gpu => self%w_aux_probe_gpu, probe_coeff_gpu => self%probe_coeff_gpu, nx => self%nx, n&
    &y => self%ny, nz => self%nz, ng => self%ng, nv_aux => self%nv_aux, w_aux_probe => self%equation_base&
    &%w_aux_probe, time_insitu => self%equation_base%time_insitu)
!
!     Save flow samples
      time_from_last_write = time_from_last_write + self%equation_base%dt
      if (time_from_last_write >= self%equation_base%dtsave) then
        if (self%masterproc) write(*,*) 'time_from_last_write=',time_from_last_write
        if (self%masterproc) write(*,*) 'istore =', self%equation_base%istore
        call self%base_gpu%copy_gpu_cpu()
        self%equation_base%w_aux = self%w_aux_gpu
        if (self%equation_base%enable_plot3d>0) then
          call self%field%write_plot3d(mach=self%equation_base%Mach, reynolds=self%equation_base%Rey&
          &nolds, time=time, istore=self%equation_base%istore, plot3dgrid=.false., plot3dfield=.true., w_aux_io&
          &=self%equation_base%w_aux(1:self%nx,1:self%ny,1:self%nz,1:6))
        endif
        if (self%equation_base%enable_vtk>0) then
          call self%field%write_vtk(time=time, istore=self%equation_base%istore, w_aux_io=self%equat&
          &ion_base%w_aux(1:self%nx,1:self%ny,1:self%nz,1:6))
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
!     Save insitu
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
!         time_from_last_insitu = 0._rkind !time_from_last_insitu - self%equation_base%dt_insitu
          time_from_last_insitu = time_from_last_insitu - self%equation_base%dt_insitu
          self%equation_base%i_insitu = self%equation_base%i_insitu + 1
        endif
      endif
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
            open(133,file='sliceyz_'//chx//'_'//chy//'_'//chz//'.bin',form='unformatted', asynchronous="yes")
          else
            open(133,file='sliceyz_'//chx//'_'//chy//'_'//chz//'.bin',form='unformatted', position="append",asynchronous="yes")
          endif
        endif
        if (allocated(self%equation_base%jslice)) then
          inquire(file='slicexz_'//chx//'_'//chy//'_'//chz//'.bin', exist=slicexz_exist)
          if(.not.(slicexz_exist)) then
            open(134,file='slicexz_'//chx//'_'//chy//'_'//chz//'.bin',form='unformatted', asynchronous="yes")
          else
            open(134,file='slicexz_'//chx//'_'//chy//'_'//chz//'.bin',form='unformatted', position="append",asynchronous="yes")
          endif
        endif
        if (allocated(self%equation_base%kslice)) then
          inquire(file='slicexy_'//chx//'_'//chy//'_'//chz//'.bin', exist=slicexy_exist)
          if(.not.(slicexy_exist)) then
            open(135,file='slicexy_'//chx//'_'//chy//'_'//chz//'.bin',form='unformatted', asynchronous="yes")
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
          call probe_interpolation_cuf(self%equation_base%num_probe,nx,ny,nz,ng,nv_aux,ijk_probe_gpu&
          &,w_aux_probe_gpu,w_aux_gpu,probe_coeff_gpu)
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
        time_from_last_rst = 0._rkind ! time_from_last_rst - self%equation_base%dtsave_restart
      endif
!
    endassociate
!   
  end subroutine manage_output
!
endmodule streams_equation_singleideal_gpu_object

