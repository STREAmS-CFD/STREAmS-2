module streams_base_omp_object

  use streams_field_object, only : field_object
  use streams_parameters
  use mpi
  use utils_omp
  use omp_lib

  implicit none
  private
  public :: base_omp_object

  type :: base_omp_object
    type(field_object), pointer :: field=>null()
    integer :: nx, ny, nz, ng, nv
    integer(ikind) :: myrank=0_ikind
    integer(ikind) :: nprocs=1_ikind
    logical :: masterproc
    integer(ikind) :: mpi_err=0_ikind
    integer(ikind) :: mydev=0_ikind
    integer(ikind) :: myhost=0_ikind
    integer(ikind) :: ierr
    integer(ikind) :: local_comm=0_ikind
    integer(ikind) :: gpu_bind=0_ikind
    real(rkind), dimension(:,:,:,:), pointer :: w_gpu

    real(rkind), dimension(:,:,:,:), pointer :: wbuf1s_gpu, wbuf2s_gpu, wbuf3s_gpu, wbuf4s_gpu, wbuf5s_gpu, wbuf6s_gpu
    real(rkind), dimension(:,:,:,:), pointer :: wbuf1r_gpu, wbuf2r_gpu, wbuf3r_gpu, wbuf4r_gpu, wbuf5r_gpu, wbuf6r_gpu
    real(rkind), dimension(:,:,:,:), pointer :: wbuf1xybcs_gpu, wbuf2xybcs_gpu, wbuf1xybcr_gpu, wbuf2xybcr_gpu
    real(rkind), dimension(:,:,:,:), pointer :: wbuf1s_c_gpu, wbuf2s_c_gpu, wbuf3s_c_gpu, wbuf4s_c_gpu
    real(rkind), dimension(:,:,:,:), pointer :: wbuf1r_c_gpu, wbuf2r_c_gpu, wbuf3r_c_gpu, wbuf4r_c_gpu
    real(rkind), dimension(:), pointer :: dcsidx_gpu, dcsidx2_gpu, dcsidxs_gpu
    real(rkind), dimension(:), pointer :: detady_gpu, detady2_gpu, detadys_gpu
    real(rkind), dimension(:), pointer :: dzitdz_gpu, dzitdz2_gpu, dzitdzs_gpu
    real(rkind), dimension(:), pointer :: x_gpu, y_gpu, z_gpu, yn_gpu
    real(rkind), dimension(:,:), pointer :: xc2_gpu,yc2_gpu
    real(rkind), dimension(:,:), pointer :: dcsidxc2_gpu,dcsidyc2_gpu
    real(rkind), dimension(:,:), pointer :: detadxc2_gpu,detadyc2_gpu
    real(rkind), dimension(:,:), pointer :: dcsidxnc2_gpu,dcsidync2_gpu
    real(rkind), dimension(:,:), pointer :: detadxnc2_gpu,detadync2_gpu
    real(rkind), dimension(:,:), pointer :: dxdcsic2_gpu,dydcsic2_gpu
    real(rkind), dimension(:,:), pointer :: dxdetac2_gpu,dydetac2_gpu
    real(rkind), dimension(:,:), pointer :: dxdcsinc2_gpu,dydcsinc2_gpu
    real(rkind), dimension(:,:), pointer :: dxdetanc2_gpu,dydetanc2_gpu
    real(rkind), dimension(:,:), pointer :: jac_gpu,mcsijac1_gpu,metajac1_gpu
    real(rkind), dimension(:,:), pointer :: mcsi_gpu,meta_gpu,csimod_gpu,etamod_gpu
    real(rkind), dimension(:,:), pointer :: g1_gpu,g2_gpu,g12_gpu
    real(rkind), dimension(:,:), pointer :: wbuftus_gpu, wbuftes_gpu , wbuftur_gpu, wbufter_gpu
    real(rkind), dimension(:,:), pointer :: theta_ij_gpu

    integer, dimension(:), pointer :: wall_tag_gpu

    real(rkind),dimension(:,:,:,:),pointer :: wbuf_lxly_s_gpu , wbuf_lxry_s_gpu , wbuf_rxly_s_gpu , wbuf_rxry_s_gpu
    real(rkind),dimension(:,:,:,:),pointer :: wbuf_lxlz_s_gpu , wbuf_lxrz_s_gpu , wbuf_rxlz_s_gpu , wbuf_rxrz_s_gpu
    real(rkind),dimension(:,:,:,:),pointer :: wbuf_lylz_s_gpu , wbuf_lyrz_s_gpu , wbuf_rylz_s_gpu , wbuf_ryrz_s_gpu
    real(rkind),dimension(:,:,:,:),pointer :: wbuf_lxlylz_s_gpu,wbuf_lxlyrz_s_gpu,wbuf_lxrylz_s_gpu,wbuf_lxryrz_s_gpu
    real(rkind),dimension(:,:,:,:),pointer :: wbuf_rxlylz_s_gpu,wbuf_rxlyrz_s_gpu,wbuf_rxrylz_s_gpu,wbuf_rxryrz_s_gpu
    real(rkind),dimension(:,:,:,:),pointer :: wbuf_lxly_r_gpu , wbuf_lxry_r_gpu , wbuf_rxly_r_gpu , wbuf_rxry_r_gpu
    real(rkind),dimension(:,:,:,:),pointer :: wbuf_lxlz_r_gpu , wbuf_lxrz_r_gpu , wbuf_rxlz_r_gpu , wbuf_rxrz_r_gpu
    real(rkind),dimension(:,:,:,:),pointer :: wbuf_lylz_r_gpu , wbuf_lyrz_r_gpu , wbuf_rylz_r_gpu , wbuf_ryrz_r_gpu
    real(rkind),dimension(:,:,:,:),pointer :: wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu
    real(rkind),dimension(:,:,:,:),pointer :: wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu

  contains
    procedure, pass(self) :: alloc
    procedure, pass(self) :: copy_cpu_gpu
    procedure, pass(self) :: copy_gpu_cpu
    procedure, pass(self) :: initialize
    procedure, pass(self) :: bcswap
    procedure, pass(self) :: bcswap_var
    procedure, pass(self) :: bcswap_corner
    procedure, pass(self) :: bcswap_corner_var
    procedure, pass(self) :: bcswap_wake
    procedure, pass(self) :: bcswap_wake_var
    procedure, pass(self) :: bcte
    procedure, pass(self) :: bcte_var
    procedure, pass(self) :: check_gpu_mem
    procedure, pass(self) :: bcswap_edges_corners
    procedure, pass(self) :: bcswap_edges_corners_var
  endtype base_omp_object

contains
  subroutine alloc(self)
    class(base_omp_object), intent(inout) :: self

    associate(nx => self%field%nx, ny => self%field%ny, nz => self%field%nz, ng => self%field%grid%ng, nv => self%field%nv)

      call omp_target_alloc_f(fptr_dev=self%w_gpu,ubounds=[nx+ng,ny+ng,nz+ng,nv],lbounds=[1-ng,1-ng,&
      &1-ng,1],omp_dev=self%mydev,ierr=self%ierr)

      call omp_target_alloc_f(fptr_dev=self%wbuf1s_gpu,ubounds=[ng,ny,nz,nv],lbounds=[1,1,1,1],omp_dev=self%mydev,ierr=self%ierr)
      call omp_target_alloc_f(fptr_dev=self%wbuf2s_gpu,ubounds=[ng,ny,nz,nv],lbounds=[1,1,1,1],omp_dev=self%mydev,ierr=self%ierr)
      call omp_target_alloc_f(fptr_dev=self%wbuf3s_gpu,ubounds=[nx,ng,nz,nv],lbounds=[1,1,1,1],omp_dev=self%mydev,ierr=self%ierr)
      call omp_target_alloc_f(fptr_dev=self%wbuf4s_gpu,ubounds=[nx,ng,nz,nv],lbounds=[1,1,1,1],omp_dev=self%mydev,ierr=self%ierr)
      call omp_target_alloc_f(fptr_dev=self%wbuf5s_gpu,ubounds=[nx,ny,ng,nv],lbounds=[1,1,1,1],omp_dev=self%mydev,ierr=self%ierr)
      call omp_target_alloc_f(fptr_dev=self%wbuf6s_gpu,ubounds=[nx,ny,ng,nv],lbounds=[1,1,1,1],omp_dev=self%mydev,ierr=self%ierr)
      call omp_target_alloc_f(fptr_dev=self%wbuf1r_gpu,ubounds=[ng,ny,nz,nv],lbounds=[1,1,1,1],omp_dev=self%mydev,ierr=self%ierr)
      call omp_target_alloc_f(fptr_dev=self%wbuf2r_gpu,ubounds=[ng,ny,nz,nv],lbounds=[1,1,1,1],omp_dev=self%mydev,ierr=self%ierr)
      call omp_target_alloc_f(fptr_dev=self%wbuf3r_gpu,ubounds=[nx,ng,nz,nv],lbounds=[1,1,1,1],omp_dev=self%mydev,ierr=self%ierr)
      call omp_target_alloc_f(fptr_dev=self%wbuf4r_gpu,ubounds=[nx,ng,nz,nv],lbounds=[1,1,1,1],omp_dev=self%mydev,ierr=self%ierr)
      call omp_target_alloc_f(fptr_dev=self%wbuf5r_gpu,ubounds=[nx,ny,ng,nv],lbounds=[1,1,1,1],omp_dev=self%mydev,ierr=self%ierr)
      call omp_target_alloc_f(fptr_dev=self%wbuf6r_gpu,ubounds=[nx,ny,ng,nv],lbounds=[1,1,1,1],omp_dev=self%mydev,ierr=self%ierr)

      call omp_target_alloc_f(fptr_dev=self%wbuf1xybcs_gpu,ubounds=[ng,ny+ng,nz,nv],lbounds=[1,1-ng,&
      &1,1],omp_dev=self%mydev,ierr=self%ierr)
      call omp_target_alloc_f(fptr_dev=self%wbuf2xybcs_gpu,ubounds=[ng,ny+ng,nz,nv],lbounds=[1,1-ng,&
      &1,1],omp_dev=self%mydev,ierr=self%ierr)
      call omp_target_alloc_f(fptr_dev=self%wbuf1xybcr_gpu,ubounds=[ng,ny+ng,nz,nv],lbounds=[1,1-ng,&
      &1,1],omp_dev=self%mydev,ierr=self%ierr)
      call omp_target_alloc_f(fptr_dev=self%wbuf2xybcr_gpu,ubounds=[ng,ny+ng,nz,nv],lbounds=[1,1-ng,&
      &1,1],omp_dev=self%mydev,ierr=self%ierr)
      call omp_target_alloc_f(fptr_dev=self%x_gpu,ubounds=[nx+ng],lbounds=[1-ng],omp_dev=self%mydev,ierr=self%ierr)
      call omp_target_alloc_f(fptr_dev=self%y_gpu,ubounds=[ny+ng],lbounds=[1-ng],omp_dev=self%mydev,ierr=self%ierr)
      call omp_target_alloc_f(fptr_dev=self%z_gpu,ubounds=[nz+ng],lbounds=[1-ng],omp_dev=self%mydev,ierr=self%ierr)

      call omp_target_alloc_f(fptr_dev=self%yn_gpu,ubounds=[ny+1],lbounds=[1],omp_dev=self%mydev,ierr=self%ierr)

      call omp_target_alloc_f(fptr_dev=self%dcsidx_gpu,ubounds=[nx],lbounds=[1],omp_dev=self%mydev,ierr=self%ierr)
      call omp_target_alloc_f(fptr_dev=self%dcsidx2_gpu,ubounds=[nx],lbounds=[1],omp_dev=self%mydev,ierr=self%ierr)
      call omp_target_alloc_f(fptr_dev=self%dcsidxs_gpu,ubounds=[nx],lbounds=[1],omp_dev=self%mydev,ierr=self%ierr)

      call omp_target_alloc_f(fptr_dev=self%detady_gpu,ubounds=[ny],lbounds=[1],omp_dev=self%mydev,ierr=self%ierr)
      call omp_target_alloc_f(fptr_dev=self%detady2_gpu,ubounds=[ny],lbounds=[1],omp_dev=self%mydev,ierr=self%ierr)
      call omp_target_alloc_f(fptr_dev=self%detadys_gpu,ubounds=[ny],lbounds=[1],omp_dev=self%mydev,ierr=self%ierr)

      call omp_target_alloc_f(fptr_dev=self%dzitdz_gpu,ubounds=[nz],lbounds=[1],omp_dev=self%mydev,ierr=self%ierr)
      call omp_target_alloc_f(fptr_dev=self%dzitdz2_gpu,ubounds=[nz],lbounds=[1],omp_dev=self%mydev,ierr=self%ierr)
      call omp_target_alloc_f(fptr_dev=self%dzitdzs_gpu,ubounds=[nz],lbounds=[1],omp_dev=self%mydev,ierr=self%ierr)


      call omp_target_alloc_f(fptr_dev=self%wbuf1s_c_gpu,ubounds=[nv,ng,ny,ng],lbounds=[1,1,1,1],omp_dev=self%mydev,ierr=self%ierr)
      call omp_target_alloc_f(fptr_dev=self%wbuf2s_c_gpu,ubounds=[nv,ng,ny,ng],lbounds=[1,1,1,1],omp_dev=self%mydev,ierr=self%ierr)
      call omp_target_alloc_f(fptr_dev=self%wbuf3s_c_gpu,ubounds=[nv,ng,ny,ng],lbounds=[1,1,1,1],omp_dev=self%mydev,ierr=self%ierr)
      call omp_target_alloc_f(fptr_dev=self%wbuf4s_c_gpu,ubounds=[nv,ng,ny,ng],lbounds=[1,1,1,1],omp_dev=self%mydev,ierr=self%ierr)
      call omp_target_alloc_f(fptr_dev=self%wbuf1r_c_gpu,ubounds=[nv,ng,ny,ng],lbounds=[1,1,1,1],omp_dev=self%mydev,ierr=self%ierr)
      call omp_target_alloc_f(fptr_dev=self%wbuf2r_c_gpu,ubounds=[nv,ng,ny,ng],lbounds=[1,1,1,1],omp_dev=self%mydev,ierr=self%ierr)
      call omp_target_alloc_f(fptr_dev=self%wbuf3r_c_gpu,ubounds=[nv,ng,ny,ng],lbounds=[1,1,1,1],omp_dev=self%mydev,ierr=self%ierr)
      call omp_target_alloc_f(fptr_dev=self%wbuf4r_c_gpu,ubounds=[nv,ng,ny,ng],lbounds=[1,1,1,1],omp_dev=self%mydev,ierr=self%ierr)

      call omp_target_alloc_f(fptr_dev=self%wbuf_lxly_s_gpu,ubounds=[nv,ng,ng,nz],lbounds=[1,1,1,1],&
      &omp_dev=self%mydev,ierr=self%ierr)
      call omp_target_alloc_f(fptr_dev=self%wbuf_lxry_s_gpu,ubounds=[nv,ng,ng,nz],lbounds=[1,1,1,1],&
      &omp_dev=self%mydev,ierr=self%ierr)

      call omp_target_alloc_f(fptr_dev=self%wbuf_rxly_s_gpu,ubounds=[nv,ng,ng,nz],lbounds=[1,1,1,1],&
      &omp_dev=self%mydev,ierr=self%ierr)
      call omp_target_alloc_f(fptr_dev=self%wbuf_rxry_s_gpu,ubounds=[nv,ng,ng,nz],lbounds=[1,1,1,1],&
      &omp_dev=self%mydev,ierr=self%ierr)

      call omp_target_alloc_f(fptr_dev=self%wbuf_lxlz_s_gpu,ubounds=[nv,ng,ny,ng],lbounds=[1,1,1,1],&
      &omp_dev=self%mydev,ierr=self%ierr)
      call omp_target_alloc_f(fptr_dev=self%wbuf_lxrz_s_gpu,ubounds=[nv,ng,ny,ng],lbounds=[1,1,1,1],&
      &omp_dev=self%mydev,ierr=self%ierr)

      call omp_target_alloc_f(fptr_dev=self%wbuf_rxlz_s_gpu,ubounds=[nv,ng,ny,ng],lbounds=[1,1,1,1],&
      &omp_dev=self%mydev,ierr=self%ierr)
      call omp_target_alloc_f(fptr_dev=self%wbuf_rxrz_s_gpu,ubounds=[nv,ng,ny,ng],lbounds=[1,1,1,1],&
      &omp_dev=self%mydev,ierr=self%ierr)

      call omp_target_alloc_f(fptr_dev=self%wbuf_lylz_s_gpu,ubounds=[nv,nx,ng,ng],lbounds=[1,1,1,1],&
      &omp_dev=self%mydev,ierr=self%ierr)
      call omp_target_alloc_f(fptr_dev=self%wbuf_lyrz_s_gpu,ubounds=[nv,nx,ng,ng],lbounds=[1,1,1,1],&
      &omp_dev=self%mydev,ierr=self%ierr)

      call omp_target_alloc_f(fptr_dev=self%wbuf_rylz_s_gpu,ubounds=[nv,nx,ng,ng],lbounds=[1,1,1,1],&
      &omp_dev=self%mydev,ierr=self%ierr)
      call omp_target_alloc_f(fptr_dev=self%wbuf_ryrz_s_gpu,ubounds=[nv,nx,ng,ng],lbounds=[1,1,1,1],&
      &omp_dev=self%mydev,ierr=self%ierr)

      call omp_target_alloc_f(fptr_dev=self%wbuf_lxlylz_s_gpu,ubounds=[nv,ng,ng,ng],lbounds=[1,1,1,&
      &1],omp_dev=self%mydev,ierr=self%ierr)
      call omp_target_alloc_f(fptr_dev=self%wbuf_lxlyrz_s_gpu,ubounds=[nv,ng,ng,ng],lbounds=[1,1,1,&
      &1],omp_dev=self%mydev,ierr=self%ierr)

      call omp_target_alloc_f(fptr_dev=self%wbuf_lxrylz_s_gpu,ubounds=[nv,ng,ng,ng],lbounds=[1,1,1,&
      &1],omp_dev=self%mydev,ierr=self%ierr)
      call omp_target_alloc_f(fptr_dev=self%wbuf_lxryrz_s_gpu,ubounds=[nv,ng,ng,ng],lbounds=[1,1,1,&
      &1],omp_dev=self%mydev,ierr=self%ierr)

      call omp_target_alloc_f(fptr_dev=self%wbuf_rxlylz_s_gpu,ubounds=[nv,ng,ng,ng],lbounds=[1,1,1,&
      &1],omp_dev=self%mydev,ierr=self%ierr)
      call omp_target_alloc_f(fptr_dev=self%wbuf_rxlyrz_s_gpu,ubounds=[nv,ng,ng,ng],lbounds=[1,1,1,&
      &1],omp_dev=self%mydev,ierr=self%ierr)

      call omp_target_alloc_f(fptr_dev=self%wbuf_rxrylz_s_gpu,ubounds=[nv,ng,ng,ng],lbounds=[1,1,1,&
      &1],omp_dev=self%mydev,ierr=self%ierr)
      call omp_target_alloc_f(fptr_dev=self%wbuf_rxryrz_s_gpu,ubounds=[nv,ng,ng,ng],lbounds=[1,1,1,&
      &1],omp_dev=self%mydev,ierr=self%ierr)


      call omp_target_alloc_f(fptr_dev=self%wbuf_lxly_r_gpu,ubounds=[nv,ng,ng,nz],lbounds=[1,1,1,1],&
      &omp_dev=self%mydev,ierr=self%ierr)
      call omp_target_alloc_f(fptr_dev=self%wbuf_lxry_r_gpu,ubounds=[nv,ng,ng,nz],lbounds=[1,1,1,1],&
      &omp_dev=self%mydev,ierr=self%ierr)

      call omp_target_alloc_f(fptr_dev=self%wbuf_rxly_r_gpu,ubounds=[nv,ng,ng,nz],lbounds=[1,1,1,1],&
      &omp_dev=self%mydev,ierr=self%ierr)
      call omp_target_alloc_f(fptr_dev=self%wbuf_rxry_r_gpu,ubounds=[nv,ng,ng,nz],lbounds=[1,1,1,1],&
      &omp_dev=self%mydev,ierr=self%ierr)

      call omp_target_alloc_f(fptr_dev=self%wbuf_lxlz_r_gpu,ubounds=[nv,ng,ny,ng],lbounds=[1,1,1,1],&
      &omp_dev=self%mydev,ierr=self%ierr)
      call omp_target_alloc_f(fptr_dev=self%wbuf_lxrz_r_gpu,ubounds=[nv,ng,ny,ng],lbounds=[1,1,1,1],&
      &omp_dev=self%mydev,ierr=self%ierr)

      call omp_target_alloc_f(fptr_dev=self%wbuf_rxlz_r_gpu,ubounds=[nv,ng,ny,ng],lbounds=[1,1,1,1],&
      &omp_dev=self%mydev,ierr=self%ierr)
      call omp_target_alloc_f(fptr_dev=self%wbuf_rxrz_r_gpu,ubounds=[nv,ng,ny,ng],lbounds=[1,1,1,1],&
      &omp_dev=self%mydev,ierr=self%ierr)

      call omp_target_alloc_f(fptr_dev=self%wbuf_lylz_r_gpu,ubounds=[nv,nx,ng,ng],lbounds=[1,1,1,1],&
      &omp_dev=self%mydev,ierr=self%ierr)
      call omp_target_alloc_f(fptr_dev=self%wbuf_lyrz_r_gpu,ubounds=[nv,nx,ng,ng],lbounds=[1,1,1,1],&
      &omp_dev=self%mydev,ierr=self%ierr)

      call omp_target_alloc_f(fptr_dev=self%wbuf_rylz_r_gpu,ubounds=[nv,nx,ng,ng],lbounds=[1,1,1,1],&
      &omp_dev=self%mydev,ierr=self%ierr)
      call omp_target_alloc_f(fptr_dev=self%wbuf_ryrz_r_gpu,ubounds=[nv,nx,ng,ng],lbounds=[1,1,1,1],&
      &omp_dev=self%mydev,ierr=self%ierr)

      call omp_target_alloc_f(fptr_dev=self%wbuf_lxlylz_r_gpu,ubounds=[nv,ng,ng,ng],lbounds=[1,1,1,&
      &1],omp_dev=self%mydev,ierr=self%ierr)
      call omp_target_alloc_f(fptr_dev=self%wbuf_lxlyrz_r_gpu,ubounds=[nv,ng,ng,ng],lbounds=[1,1,1,&
      &1],omp_dev=self%mydev,ierr=self%ierr)

      call omp_target_alloc_f(fptr_dev=self%wbuf_lxrylz_r_gpu,ubounds=[nv,ng,ng,ng],lbounds=[1,1,1,&
      &1],omp_dev=self%mydev,ierr=self%ierr)
      call omp_target_alloc_f(fptr_dev=self%wbuf_lxryrz_r_gpu,ubounds=[nv,ng,ng,ng],lbounds=[1,1,1,&
      &1],omp_dev=self%mydev,ierr=self%ierr)

      call omp_target_alloc_f(fptr_dev=self%wbuf_rxlylz_r_gpu,ubounds=[nv,ng,ng,ng],lbounds=[1,1,1,&
      &1],omp_dev=self%mydev,ierr=self%ierr)
      call omp_target_alloc_f(fptr_dev=self%wbuf_rxlyrz_r_gpu,ubounds=[nv,ng,ng,ng],lbounds=[1,1,1,&
      &1],omp_dev=self%mydev,ierr=self%ierr)

      call omp_target_alloc_f(fptr_dev=self%wbuf_rxrylz_r_gpu,ubounds=[nv,ng,ng,ng],lbounds=[1,1,1,&
      &1],omp_dev=self%mydev,ierr=self%ierr)
      call omp_target_alloc_f(fptr_dev=self%wbuf_rxryrz_r_gpu,ubounds=[nv,ng,ng,ng],lbounds=[1,1,1,&
      &1],omp_dev=self%mydev,ierr=self%ierr)


      if (self%field%grid%grid_dim == 2) then
        call omp_target_alloc_f(fptr_dev=self%xc2_gpu,ubounds=[nx+ng,ny+ng],lbounds=[1-ng,1-ng],omp_dev=self%mydev,ierr=self%ierr)
        call omp_target_alloc_f(fptr_dev=self%yc2_gpu,ubounds=[nx+ng,ny+ng],lbounds=[1-ng,1-ng],omp_dev=self%mydev,ierr=self%ierr)
        call omp_target_alloc_f(fptr_dev=self%dcsidxc2_gpu,ubounds=[nx+ng,ny+ng],lbounds=[1-ng,1-&
        &ng],omp_dev=self%mydev,ierr=self%ierr)
        call omp_target_alloc_f(fptr_dev=self%dcsidyc2_gpu,ubounds=[nx+ng,ny+ng],lbounds=[1-ng,1-&
        &ng],omp_dev=self%mydev,ierr=self%ierr)
        call omp_target_alloc_f(fptr_dev=self%detadxc2_gpu,ubounds=[nx+ng,ny+ng],lbounds=[1-ng,1-&
        &ng],omp_dev=self%mydev,ierr=self%ierr)
        call omp_target_alloc_f(fptr_dev=self%detadyc2_gpu,ubounds=[nx+ng,ny+ng],lbounds=[1-ng,1-&
        &ng],omp_dev=self%mydev,ierr=self%ierr)
        call omp_target_alloc_f(fptr_dev=self%dcsidxnc2_gpu,ubounds=[nx+ng,ny+ng],lbounds=[1-ng,1-&
        &ng],omp_dev=self%mydev,ierr=self%ierr)
        call omp_target_alloc_f(fptr_dev=self%dcsidync2_gpu,ubounds=[nx+ng,ny+ng],lbounds=[1-ng,1-&
        &ng],omp_dev=self%mydev,ierr=self%ierr)
        call omp_target_alloc_f(fptr_dev=self%detadxnc2_gpu,ubounds=[nx+ng,ny+ng],lbounds=[1-ng,1-&
        &ng],omp_dev=self%mydev,ierr=self%ierr)
        call omp_target_alloc_f(fptr_dev=self%detadync2_gpu,ubounds=[nx+ng,ny+ng],lbounds=[1-ng,1-&
        &ng],omp_dev=self%mydev,ierr=self%ierr)
        call omp_target_alloc_f(fptr_dev=self%dxdcsic2_gpu,ubounds=[nx+ng,ny+ng],lbounds=[1-ng,1-&
        &ng],omp_dev=self%mydev,ierr=self%ierr)
        call omp_target_alloc_f(fptr_dev=self%dydcsic2_gpu,ubounds=[nx+ng,ny+ng],lbounds=[1-ng,1-&
        &ng],omp_dev=self%mydev,ierr=self%ierr)
        call omp_target_alloc_f(fptr_dev=self%dxdetac2_gpu,ubounds=[nx+ng,ny+ng],lbounds=[1-ng,1-&
        &ng],omp_dev=self%mydev,ierr=self%ierr)
        call omp_target_alloc_f(fptr_dev=self%dydetac2_gpu,ubounds=[nx+ng,ny+ng],lbounds=[1-ng,1-&
        &ng],omp_dev=self%mydev,ierr=self%ierr)
        call omp_target_alloc_f(fptr_dev=self%dxdcsinc2_gpu,ubounds=[nx+ng,ny+ng],lbounds=[1-ng,1-&
        &ng],omp_dev=self%mydev,ierr=self%ierr)
        call omp_target_alloc_f(fptr_dev=self%dydcsinc2_gpu,ubounds=[nx+ng,ny+ng],lbounds=[1-ng,1-&
        &ng],omp_dev=self%mydev,ierr=self%ierr)
        call omp_target_alloc_f(fptr_dev=self%dxdetanc2_gpu,ubounds=[nx+ng,ny+ng],lbounds=[1-ng,1-&
        &ng],omp_dev=self%mydev,ierr=self%ierr)
        call omp_target_alloc_f(fptr_dev=self%dydetanc2_gpu,ubounds=[nx+ng,ny+ng],lbounds=[1-ng,1-&
        &ng],omp_dev=self%mydev,ierr=self%ierr)
        call omp_target_alloc_f(fptr_dev=self%mcsi_gpu,ubounds=[nx+ng,ny+ng],lbounds=[1-ng,1-ng],omp_dev=self%mydev,ierr=self%ierr)
        call omp_target_alloc_f(fptr_dev=self%meta_gpu,ubounds=[nx+ng,ny+ng],lbounds=[1-ng,1-ng],omp_dev=self%mydev,ierr=self%ierr)
        call omp_target_alloc_f(fptr_dev=self%mcsijac1_gpu,ubounds=[nx+ng,ny+ng],lbounds=[1-ng,1-&
        &ng],omp_dev=self%mydev,ierr=self%ierr)
        call omp_target_alloc_f(fptr_dev=self%metajac1_gpu,ubounds=[nx+ng,ny+ng],lbounds=[1-ng,1-&
        &ng],omp_dev=self%mydev,ierr=self%ierr)
        call omp_target_alloc_f(fptr_dev=self%jac_gpu,ubounds=[nx+ng,ny+ng],lbounds=[1-ng,1-ng],omp_dev=self%mydev,ierr=self%ierr)
        call omp_target_alloc_f(fptr_dev=self%g1_gpu,ubounds=[nx+ng,ny+ng],lbounds=[1-ng,1-ng],omp_dev=self%mydev,ierr=self%ierr)
        call omp_target_alloc_f(fptr_dev=self%g2_gpu,ubounds=[nx+ng,ny+ng],lbounds=[1-ng,1-ng],omp_dev=self%mydev,ierr=self%ierr)
        call omp_target_alloc_f(fptr_dev=self%g12_gpu,ubounds=[nx+ng,ny+ng],lbounds=[1-ng,1-ng],omp_dev=self%mydev,ierr=self%ierr)
        call omp_target_alloc_f(fptr_dev=self%csimod_gpu,ubounds=[nx+ng,ny+ng],lbounds=[1-ng,1-ng],&
        &omp_dev=self%mydev,ierr=self%ierr)
        call omp_target_alloc_f(fptr_dev=self%etamod_gpu,ubounds=[nx+ng,ny+ng],lbounds=[1-ng,1-ng],&
        &omp_dev=self%mydev,ierr=self%ierr)
        call omp_target_alloc_f(fptr_dev=self%theta_ij_gpu,ubounds=[nx,ny],lbounds=[1,1],omp_dev=self%mydev,ierr=self%ierr)

        call omp_target_alloc_f(fptr_dev=self%wall_tag_gpu,ubounds=[nx+ng],lbounds=[1-ng],omp_dev=self%mydev,ierr=self%ierr)
        call omp_target_alloc_f(fptr_dev=self%wbuftus_gpu,ubounds=[nz,2],lbounds=[1,1],omp_dev=self%mydev,ierr=self%ierr)
        call omp_target_alloc_f(fptr_dev=self%wbuftes_gpu,ubounds=[nz,2],lbounds=[1,1],omp_dev=self%mydev,ierr=self%ierr)
        call omp_target_alloc_f(fptr_dev=self%wbuftur_gpu,ubounds=[nz,2],lbounds=[1,1],omp_dev=self%mydev,ierr=self%ierr)
        call omp_target_alloc_f(fptr_dev=self%wbufter_gpu,ubounds=[nz,2],lbounds=[1,1],omp_dev=self%mydev,ierr=self%ierr)
      endif

    endassociate

    self%ierr = omp_target_memcpy_f(self%x_gpu , self%field%x,0,0,self%mydev,self%myhost)
    self%ierr = omp_target_memcpy_f(self%y_gpu , self%field%y,0,0,self%mydev,self%myhost)
    self%ierr = omp_target_memcpy_f(self%z_gpu , self%field%z,0,0,self%mydev,self%myhost)
    self%ierr = omp_target_memcpy_f(self%yn_gpu , self%field%yn,0,0,self%mydev,self%myhost)

    self%ierr = omp_target_memcpy_f(self%dcsidx_gpu , self%field%dcsidx,0,0,self%mydev,self%myhost)
    self%ierr = omp_target_memcpy_f(self%dcsidxs_gpu , self%field%dcsidxs,0,0,self%mydev,self%myhost)
    self%ierr = omp_target_memcpy_f(self%dcsidx2_gpu , self%field%dcsidx2,0,0,self%mydev,self%myhost)
    self%ierr = omp_target_memcpy_f(self%detady_gpu , self%field%detady,0,0,self%mydev,self%myhost)
    self%ierr = omp_target_memcpy_f(self%detadys_gpu , self%field%detadys,0,0,self%mydev,self%myhost)
    self%ierr = omp_target_memcpy_f(self%detady2_gpu , self%field%detady2,0,0,self%mydev,self%myhost)
    self%ierr = omp_target_memcpy_f(self%dzitdz_gpu , self%field%dzitdz,0,0,self%mydev,self%myhost)
    self%ierr = omp_target_memcpy_f(self%dzitdzs_gpu , self%field%dzitdzs,0,0,self%mydev,self%myhost)
    self%ierr = omp_target_memcpy_f(self%dzitdz2_gpu , self%field%dzitdz2,0,0,self%mydev,self%myhost)

    if (self%field%grid%grid_dim == 2) then
      self%ierr = omp_target_memcpy_f(self%xc2_gpu , self%field%xc2,0,0,self%mydev,self%myhost)
      self%ierr = omp_target_memcpy_f(self%yc2_gpu , self%field%yc2,0,0,self%mydev,self%myhost)
      self%ierr = omp_target_memcpy_f(self%dcsidxc2_gpu , self%field%dcsidxc2,0,0,self%mydev,self%myhost)
      self%ierr = omp_target_memcpy_f(self%dcsidyc2_gpu , self%field%dcsidyc2,0,0,self%mydev,self%myhost)
      self%ierr = omp_target_memcpy_f(self%detadxc2_gpu , self%field%detadxc2,0,0,self%mydev,self%myhost)
      self%ierr = omp_target_memcpy_f(self%detadyc2_gpu , self%field%detadyc2,0,0,self%mydev,self%myhost)
      self%ierr = omp_target_memcpy_f(self%dcsidxnc2_gpu , self%field%dcsidxnc2,0,0,self%mydev,self%myhost)
      self%ierr = omp_target_memcpy_f(self%dcsidync2_gpu , self%field%dcsidync2,0,0,self%mydev,self%myhost)
      self%ierr = omp_target_memcpy_f(self%detadxnc2_gpu , self%field%detadxnc2,0,0,self%mydev,self%myhost)
      self%ierr = omp_target_memcpy_f(self%detadync2_gpu , self%field%detadync2,0,0,self%mydev,self%myhost)
      self%ierr = omp_target_memcpy_f(self%dxdcsic2_gpu , self%field%dxdcsic2,0,0,self%mydev,self%myhost)
      self%ierr = omp_target_memcpy_f(self%dydcsic2_gpu , self%field%dydcsic2,0,0,self%mydev,self%myhost)
      self%ierr = omp_target_memcpy_f(self%dxdetac2_gpu , self%field%dxdetac2,0,0,self%mydev,self%myhost)
      self%ierr = omp_target_memcpy_f(self%dydetac2_gpu , self%field%dydetac2,0,0,self%mydev,self%myhost)
      self%ierr = omp_target_memcpy_f(self%dxdcsinc2_gpu , self%field%dxdcsinc2,0,0,self%mydev,self%myhost)
      self%ierr = omp_target_memcpy_f(self%dydcsinc2_gpu , self%field%dydcsinc2,0,0,self%mydev,self%myhost)
      self%ierr = omp_target_memcpy_f(self%dxdetanc2_gpu , self%field%dxdetanc2,0,0,self%mydev,self%myhost)
      self%ierr = omp_target_memcpy_f(self%dydetanc2_gpu , self%field%dydetanc2,0,0,self%mydev,self%myhost)
      self%ierr = omp_target_memcpy_f(self%mcsi_gpu , self%field%mcsi,0,0,self%mydev,self%myhost)
      self%ierr = omp_target_memcpy_f(self%meta_gpu , self%field%meta,0,0,self%mydev,self%myhost)
      self%ierr = omp_target_memcpy_f(self%mcsijac1_gpu , self%field%mcsijac1,0,0,self%mydev,self%myhost)
      self%ierr = omp_target_memcpy_f(self%metajac1_gpu , self%field%metajac1,0,0,self%mydev,self%myhost)
      self%ierr = omp_target_memcpy_f(self%jac_gpu , self%field%jac,0,0,self%mydev,self%myhost)
      self%ierr = omp_target_memcpy_f(self%g1_gpu , self%field%g1,0,0,self%mydev,self%myhost)
      self%ierr = omp_target_memcpy_f(self%g2_gpu , self%field%g2,0,0,self%mydev,self%myhost)
      self%ierr = omp_target_memcpy_f(self%g12_gpu , self%field%g12,0,0,self%mydev,self%myhost)
      self%ierr = omp_target_memcpy_f(self%csimod_gpu , self%field%csimod,0,0,self%mydev,self%myhost)
      self%ierr = omp_target_memcpy_f(self%etamod_gpu , self%field%etamod,0,0,self%mydev,self%myhost)
      self%ierr = omp_target_memcpy_f(self%wall_tag_gpu , self%field%wall_tag,0,0,self%mydev,self%myhost)
      self%ierr = omp_target_memcpy_f(self%theta_ij_gpu , self%field%theta_ij,0,0,self%mydev,self%myhost)
    endif

  endsubroutine alloc

  subroutine bcte_step_1_kernel(nx,ny,nz,nv,ng,nrank_x,ite_rank_x,itu_rank_x,ite_l,itu_l,w_gpu,wbuftus_gpu,wbuftes_gpu)
    integer, intent(in) :: nx, ny, nz, nv, ng
    integer, intent(in) :: nrank_x, ite_rank_x, itu_rank_x, ite_l, itu_l
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in) :: w_gpu
    real(rkind), dimension(nz, 2), intent(out) :: wbuftus_gpu, wbuftes_gpu
    integer :: i,j,k,m,iercuda
    if(nrank_x == ite_rank_x) then
      !$omp target teams distribute parallel do collapse(1) has_device_addr(w_gpu,wbuftus_gpu,wbuftes_gpu)
      do k = 1,nz
        wbuftes_gpu(k,1) = w_gpu(ite_l,1,k,1)
        if(nv > 1) wbuftes_gpu(k,2) = w_gpu(ite_l,1,k,5)
      enddo
    endif
    if(nrank_x == itu_rank_x) then
      !$omp target teams distribute parallel do collapse(1) has_device_addr(w_gpu,wbuftus_gpu,wbuftes_gpu)
      do k = 1,nz
        wbuftus_gpu(k,1) = w_gpu(itu_l,1,k,1)
        if(nv > 1) wbuftus_gpu(k,2) = w_gpu(itu_l,1,k,5)
      enddo
    endif

  endsubroutine bcte_step_1_kernel



  subroutine bcte(self, steps)
    class(base_omp_object), intent(inout) :: self
    logical, dimension(3), optional :: steps
    logical, dimension(3) :: steps_
    integer, dimension(4) :: requests

    steps_ = .true. ; if(present(steps)) steps_ = steps
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv,&
    & w_gpu => self%w_gpu, wbuftes_gpu => self%wbuftes_gpu, wbuftus_gpu => self%wbuftus_gpu,&
    & wbufter_gpu => self%wbufter_gpu, wbuftur_gpu => self%wbuftur_gpu, ite_rank_x => self%field%ite_rank&
    &_x, itu_rank_x => self%field%itu_rank_x, ite_l => self%field%ite_l, itu_l => self%field%itu_l,&
    & nrank_x => self%field%nrank_x, mp_cartx => self%field%mp_cartx, iermpi => self%mpi_err)

      if(steps_(1)) then
        call bcte_step_1_kernel(nx, ny, nz, nv, ng, nrank_x, ite_rank_x, itu_rank_x, ite_l, itu_l, w_gpu, wbuftus_gpu, wbuftes_gpu)
      endif

      if(steps_(2)) then
        if (nrank_x == ite_rank_x .and. nrank_x == itu_rank_x) then
          self%ierr = omp_target_memcpy_f(wbufter_gpu , wbuftes_gpu,0,0,self%mydev,self%mydev)
          self%ierr = omp_target_memcpy_f(wbuftur_gpu , wbuftus_gpu,0,0,self%mydev,self%mydev)
        else
          requests(:) = mpi_request_null
          if (nrank_x == ite_rank_x) then
            call mpi_isend(wbuftes_gpu,2*nz,mpi_prec,itu_rank_x,1,mp_cartx,requests(1),iermpi)
            call mpi_irecv(wbuftur_gpu,2*nz,mpi_prec,itu_rank_x,1,mp_cartx,requests(2),iermpi)
          elseif (nrank_x == itu_rank_x) then
            call mpi_isend(wbuftus_gpu,2*nz,mpi_prec,ite_rank_x,1,mp_cartx,requests(3),iermpi)
            call mpi_irecv(wbufter_gpu,2*nz,mpi_prec,ite_rank_x,1,mp_cartx,requests(4),iermpi)
          endif
          call mpi_waitall(4,requests,mpi_statuses_ignore,iermpi)
        endif
      endif

      if(steps_(3)) then
        call bcte_step_3_kernel(nx, ny, nz, nv, ng, nrank_x, ite_rank_x, itu_rank_x, ite_l, itu_l, w_gpu, wbuftur_gpu, wbufter_gpu)
      endif
    endassociate
  endsubroutine bcte

  subroutine bcte_var(self, w_swap, steps)
    class(base_omp_object), intent(inout) :: self
    real(rkind), dimension(1-self%ng:self%nx+self%ng, 1-self%ng:self%ny+self%ng, 1-self%ng:self%nz+self%ng,1:1) :: w_swap
    logical, dimension(3), optional :: steps
    logical, dimension(3) :: steps_
    integer, dimension(4) :: requests

    steps_ = .true. ; if(present(steps)) steps_ = steps
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv,&
    & w_gpu => self%w_gpu, wbuftes_gpu => self%wbuftes_gpu, wbuftus_gpu => self%wbuftus_gpu,&
    & wbufter_gpu => self%wbufter_gpu, wbuftur_gpu => self%wbuftur_gpu, ite_rank_x => self%field%ite_rank&
    &_x, itu_rank_x => self%field%itu_rank_x, ite_l => self%field%ite_l, itu_l => self%field%itu_l,&
    & nrank_x => self%field%nrank_x, mp_cartx => self%field%mp_cartx, iermpi => self%mpi_err)

      if(steps_(1)) then
        call bcte_step_1_kernel(nx, ny, nz, 1, ng, nrank_x, ite_rank_x, itu_rank_x, ite_l, itu_l, w_swap, wbuftus_gpu, wbuftes_gpu)
      endif

      if(steps_(2)) then
        if (nrank_x == ite_rank_x .and. nrank_x == itu_rank_x) then
          self%ierr = omp_target_memcpy_f(wbufter_gpu , wbuftes_gpu,0,0,self%mydev,self%mydev)
          self%ierr = omp_target_memcpy_f(wbuftur_gpu , wbuftus_gpu,0,0,self%mydev,self%mydev)
        else
          requests(:) = mpi_request_null
          if (nrank_x == ite_rank_x) then
            call mpi_isend(wbuftes_gpu,1*nz,mpi_prec,itu_rank_x,3,mp_cartx,requests(1),iermpi)
            call mpi_irecv(wbuftur_gpu,1*nz,mpi_prec,itu_rank_x,3,mp_cartx,requests(2),iermpi)
          elseif (nrank_x == itu_rank_x) then
            call mpi_isend(wbuftus_gpu,1*nz,mpi_prec,ite_rank_x,3,mp_cartx,requests(3),iermpi)
            call mpi_irecv(wbufter_gpu,1*nz,mpi_prec,ite_rank_x,3,mp_cartx,requests(4),iermpi)
          endif
          call mpi_waitall(4,requests,mpi_statuses_ignore,iermpi)
        endif
      endif

      if(steps_(3)) then
        call bcte_step_3_kernel(nx, ny, nz, 1, ng, nrank_x, ite_rank_x, itu_rank_x, ite_l, itu_l, w_swap, wbuftur_gpu, wbufter_gpu)
      endif
    endassociate
  endsubroutine bcte_var

  subroutine bcte_step_3_kernel(nx,ny,nz,nv,ng,nrank_x,ite_rank_x,itu_rank_x,ite_l,itu_l,w_gpu,wbuftur_gpu,wbufter_gpu)
    integer, intent(in) :: nx, ny, nz, nv, ng
    integer, intent(in) :: nrank_x, ite_rank_x, itu_rank_x, ite_l, itu_l
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(inout) :: w_gpu
    real(rkind), dimension(nz, 2), intent(in) :: wbuftur_gpu, wbufter_gpu
    integer :: k, iercuda
    if(nrank_x == ite_rank_x) then
      !$omp target teams distribute parallel do collapse(1) has_device_addr(w_gpu,wbuftur_gpu,wbufter_gpu)
      do k = 1,nz
        w_gpu(ite_l,1,k,1) = 0.5_rkind*(w_gpu(ite_l,1,k,1)+wbuftur_gpu(k,1))
        if(nv > 1) w_gpu(ite_l,1,k,5) = 0.5_rkind*(w_gpu(ite_l,1,k,5)+wbuftur_gpu(k,2))
      enddo
    endif

    if (nrank_x==itu_rank_x) then
      !$omp target teams distribute parallel do collapse(1) has_device_addr(w_gpu,wbuftur_gpu,wbufter_gpu)
      do k = 1,nz
        w_gpu(itu_l,1,k,1) = 0.5_rkind*(w_gpu(itu_l,1,k,1)+wbufter_gpu(k,1))
        if(nv > 1) w_gpu(itu_l,1,k,5) = 0.5_rkind*(w_gpu(itu_l,1,k,5)+wbufter_gpu(k,2))
      enddo
    endif

  endsubroutine bcte_step_3_kernel



  subroutine bcswap_wake_step_1_kernel(nx,ny,nz,nv,ng,w_gpu,wbuf4s_gpu)
    integer :: nx, ny, nz, ng, nv
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: w_gpu
    real(rkind), dimension(1:,1:,1:,1:) :: wbuf4s_gpu
    integer :: i,j,k,m,iercuda

    !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu,wbuf4s_gpu)
    do k = 1,nz
      do j = 1,ng
        do i = 1,nx
          do m=1,nv
            wbuf4s_gpu(i,j,k,m) = w_gpu(nx-i+1,1+j,k,m)
          enddo
        enddo
      enddo
    enddo


  endsubroutine bcswap_wake_step_1_kernel



  subroutine bcswap_wake_step_3_kernel(nx,ny,nz,nv,ng,w_gpu,wbuf3r_gpu,wall_tag_gpu)
    integer :: nx, ny, nz, ng, nv
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: w_gpu
    real(rkind), dimension(1:,1:,1:,1:) :: wbuf3r_gpu
    integer, dimension(1-ng:nx+ng) :: wall_tag_gpu
    integer :: i,j,k,m, iercuda

    !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu,wbuf3r_gpu,wall_tag_gpu)
    do k = 1,nz
      do j = 1,ng
        do i = 1,nx
          do m=1,nv
            if (wall_tag_gpu(i) > 0) then
              w_gpu(i,1-j,k,m) = wbuf3r_gpu(i,j,k,m)
            endif
          enddo
        enddo
      enddo
    enddo


  endsubroutine bcswap_wake_step_3_kernel



  subroutine bcswap_step_1_kernel(nx,ny,nz,nv,ng,ndim,w_gpu,wbuf1s_gpu,wbuf2s_gpu,wbuf3s_gpu,wbuf4s_gpu,wbuf5s_gpu,wbuf6s_gpu)
    integer :: nx, ny, nz, ng, nv, ndim
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: w_gpu
    real(rkind), dimension(1:,1:,1:,1:) :: wbuf1s_gpu, wbuf2s_gpu, wbuf3s_gpu, wbuf4s_gpu, wbuf5s_gpu, wbuf6s_gpu
    integer :: i,j,k,m,iercuda

    !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu,wbuf1s_gpu,wbuf2s_gpu,wbuf3s_gpu,wbuf4s_gpu,wbuf5s_gpu,wbuf6s_gpu)
    do k = 1,nz
      do j = 1,ny
        do i = 1,ng
          do m=1,nv
            wbuf1s_gpu(i,j,k,m) = w_gpu(i,j,k,m)
            wbuf2s_gpu(i,j,k,m) = w_gpu(nx-ng+i,j,k,m)
          enddo
        enddo
      enddo
    enddo

    !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu,wbuf1s_gpu,wbuf2s_gpu,wbuf3s_gpu,wbuf4s_gpu,wbuf5s_gpu,wbuf6s_gpu)
    do k = 1,nz
      do j = 1,ng
        do i = 1,nx
          do m=1,nv
            wbuf3s_gpu(i,j,k,m) = w_gpu(i,j,k,m)
            wbuf4s_gpu(i,j,k,m) = w_gpu(i,ny-ng+j,k,m)
          enddo
        enddo
      enddo
    enddo
    if (ndim==3) then
      !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu,wbuf1s_gpu,wbuf2s_gpu,wbuf3s_gpu,wbuf4s_gpu,wbuf5s_gpu,wbuf6s_gpu)
      do k = 1,ng
        do j = 1,ny
          do i = 1,nx
            do m=1,nv
              wbuf5s_gpu(i,j,k,m) = w_gpu(i,j,k,m)
              wbuf6s_gpu(i,j,k,m) = w_gpu(i,j,nz-ng+k,m)
            enddo
          enddo
        enddo
      enddo
    endif

  endsubroutine bcswap_step_1_kernel



  subroutine bcswap_c2_step_1_kernel(nx,ny,nz,nv,ng,ndim,w_gpu,wbuf1s_gpu,wbuf2s_gpu,wbuf3s_gpu,&
  &wbuf4s_gpu,wbuf5s_gpu,wbuf6s_gpu,dcsidxnc2_gpu,dcsidync2_gpu,detadxnc2_gpu,detadync2_gpu,&
  &is_periodic)
    integer :: nx, ny, nz, ng, nv, ndim
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: w_gpu
    real(rkind), dimension(1:,1:,1:,1:) :: wbuf1s_gpu, wbuf2s_gpu, wbuf3s_gpu, wbuf4s_gpu, wbuf5s_gpu, wbuf6s_gpu
    real(rkind), dimension(1-ng:,1-ng:) :: dcsidxnc2_gpu, dcsidync2_gpu, detadxnc2_gpu, detadync2_gpu
    integer, dimension(3) :: is_periodic
    integer :: i,j,k,m,iercuda
    if (is_periodic(1) == 1) then
      !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu,wbuf1s_gpu,wbuf2s_gpu,wbuf3s_gpu,wbuf4s_gpu,wbuf5s_gpu,wbuf6s_gpu,dcsidxnc2_gpu,dcsidync2_gpu,detadxnc2_gpu,detadync2_gpu)
      do k = 1,nz
        do j = 1,ny
          do i = 1,ng
            wbuf1s_gpu(i,j,k,1) = w_gpu(i,j,k,1)
            wbuf1s_gpu(i,j,k,2) = w_gpu(i,j,k,2)*dcsidxnc2_gpu(i,j)+w_gpu(i,j,k,3)*dcsidync2_gpu(i,j)
            wbuf1s_gpu(i,j,k,3) = w_gpu(i,j,k,2)*detadxnc2_gpu(i,j)+w_gpu(i,j,k,3)*detadync2_gpu(i,j)
            wbuf1s_gpu(i,j,k,4) = w_gpu(i,j,k,4)
            wbuf1s_gpu(i,j,k,5) = w_gpu(i,j,k,5)
            wbuf2s_gpu(i,j,k,1) = w_gpu(nx-ng+i,j,k,1)
            wbuf2s_gpu(i,j,k,2) = w_gpu(nx-ng+i,j,k,2)*dcsidxnc2_gpu(nx-ng+i,j)+w_gpu(nx-ng+i,j,k,3)*dcsidync2_gpu(nx-ng+i,j)
            wbuf2s_gpu(i,j,k,3) = w_gpu(nx-ng+i,j,k,2)*detadxnc2_gpu(nx-ng+i,j)+w_gpu(nx-ng+i,j,k,3)*detadync2_gpu(nx-ng+i,j)
            wbuf2s_gpu(i,j,k,4) = w_gpu(nx-ng+i,j,k,4)
            wbuf2s_gpu(i,j,k,5) = w_gpu(nx-ng+i,j,k,5)
          enddo
        enddo
      enddo
    else
      !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu,wbuf1s_gpu,wbuf2s_gpu,wbuf3s_gpu,wbuf4s_gpu,wbuf5s_gpu,wbuf6s_gpu,dcsidxnc2_gpu,dcsidync2_gpu,detadxnc2_gpu,detadync2_gpu)
      do k = 1,nz
        do j = 1,ny
          do i = 1,ng
            do m=1,nv
              wbuf1s_gpu(i,j,k,m) = w_gpu(i,j,k,m)
              wbuf2s_gpu(i,j,k,m) = w_gpu(nx-ng+i,j,k,m)
            enddo
          enddo
        enddo
      enddo
    endif
    !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu,wbuf1s_gpu,wbuf2s_gpu,wbuf3s_gpu,wbuf4s_gpu,wbuf5s_gpu,wbuf6s_gpu,dcsidxnc2_gpu,dcsidync2_gpu,detadxnc2_gpu,detadync2_gpu)
    do k = 1,nz
      do j = 1,ng
        do i = 1,nx
          do m=1,nv
            wbuf3s_gpu(i,j,k,m) = w_gpu(i,j,k,m)
            wbuf4s_gpu(i,j,k,m) = w_gpu(i,ny-ng+j,k,m)
          enddo
        enddo
      enddo
    enddo
    if (ndim==3) then
      !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu,wbuf1s_gpu,wbuf2s_gpu,wbuf3s_gpu,wbuf4s_gpu,wbuf5s_gpu,wbuf6s_gpu,dcsidxnc2_gpu,dcsidync2_gpu,detadxnc2_gpu,detadync2_gpu)
      do k = 1,ng
        do j = 1,ny
          do i = 1,nx
            do m=1,nv
              wbuf5s_gpu(i,j,k,m) = w_gpu(i,j,k,m)
              wbuf6s_gpu(i,j,k,m) = w_gpu(i,j,nz-ng+k,m)
            enddo
          enddo
        enddo
      enddo
    endif

  endsubroutine bcswap_c2_step_1_kernel



  subroutine bcswap_c2xybc_step_1_kernel(nx,ny,nz,nv,ng,ndim,w_gpu,wbuf1xybcs_gpu,wbuf2xybcs_gpu,&
  &wbuf3s_gpu,wbuf4s_gpu,wbuf5s_gpu,wbuf6s_gpu,dcsidxnc2_gpu,dcsidync2_gpu,detadxnc2_gpu,detadync2_gpu,&
  &is_periodic)
    integer :: nx, ny, nz, ng, nv, ndim
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: w_gpu
    real(rkind), dimension(1:,1-ng:,1:,1:) :: wbuf1xybcs_gpu, wbuf2xybcs_gpu
    real(rkind), dimension(1:,1:,1:,1:) :: wbuf3s_gpu, wbuf4s_gpu, wbuf5s_gpu, wbuf6s_gpu
    real(rkind), dimension(1-ng:,1-ng:) :: dcsidxnc2_gpu, dcsidync2_gpu, detadxnc2_gpu, detadync2_gpu
    integer, dimension(3) :: is_periodic
    integer :: i,j,k,m,iercuda
    if (is_periodic(1) == 1) then
      !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu,wbuf1xybcs_gpu,wbuf2xybcs_gpu,wbuf3s_gpu,wbuf4s_gpu,wbuf5s_gpu,wbuf6s_gpu,dcsidxnc2_gpu,dcsidync2_gpu,detadxnc2_gpu,detadync2_gpu)
      do k = 1,nz
        do j = 1-ng,ny+ng
          do i = 1,ng
            wbuf1xybcs_gpu(i,j,k,1) = w_gpu(i,j,k,1)
            wbuf1xybcs_gpu(i,j,k,2) = w_gpu(i,j,k,2)*dcsidxnc2_gpu(i,j)+w_gpu(i,j,k,3)*dcsidync2_gpu(i,j)
            wbuf1xybcs_gpu(i,j,k,3) = w_gpu(i,j,k,2)*detadxnc2_gpu(i,j)+w_gpu(i,j,k,3)*detadync2_gpu(i,j)
            wbuf1xybcs_gpu(i,j,k,4) = w_gpu(i,j,k,4)
            wbuf1xybcs_gpu(i,j,k,5) = w_gpu(i,j,k,5)
            wbuf2xybcs_gpu(i,j,k,1) = w_gpu(nx-ng+i,j,k,1)
            wbuf2xybcs_gpu(i,j,k,2) = w_gpu(nx-ng+i,j,k,2)*dcsidxnc2_gpu(nx-ng+i,j)+w_gpu(nx-ng+i,j,k,3)*dcsidync2_gpu(nx-ng+i,j)
            wbuf2xybcs_gpu(i,j,k,3) = w_gpu(nx-ng+i,j,k,2)*detadxnc2_gpu(nx-ng+i,j)+w_gpu(nx-ng+i,j,k,3)*detadync2_gpu(nx-ng+i,j)
            wbuf2xybcs_gpu(i,j,k,4) = w_gpu(nx-ng+i,j,k,4)
            wbuf2xybcs_gpu(i,j,k,5) = w_gpu(nx-ng+i,j,k,5)
          enddo
        enddo
      enddo
    else
      !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu,wbuf1xybcs_gpu,wbuf2xybcs_gpu,wbuf3s_gpu,wbuf4s_gpu,wbuf5s_gpu,wbuf6s_gpu,dcsidxnc2_gpu,dcsidync2_gpu,detadxnc2_gpu,detadync2_gpu)
      do k = 1,nz
        do j = 1-ng,ny+ng
          do i = 1,ng
            do m=1,nv
              wbuf1xybcs_gpu(i,j,k,m) = w_gpu(i,j,k,m)
              wbuf2xybcs_gpu(i,j,k,m) = w_gpu(nx-ng+i,j,k,m)
            enddo
          enddo
        enddo
      enddo
    endif
    !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu,wbuf1xybcs_gpu,wbuf2xybcs_gpu,wbuf3s_gpu,wbuf4s_gpu,wbuf5s_gpu,wbuf6s_gpu,dcsidxnc2_gpu,dcsidync2_gpu,detadxnc2_gpu,detadync2_gpu)
    do k = 1,nz
      do j = 1,ng
        do i = 1,nx
          do m=1,nv
            wbuf3s_gpu(i,j,k,m) = w_gpu(i,j,k,m)
            wbuf4s_gpu(i,j,k,m) = w_gpu(i,ny-ng+j,k,m)
          enddo
        enddo
      enddo
    enddo
    if (ndim==3) then
      !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu,wbuf1xybcs_gpu,wbuf2xybcs_gpu,wbuf3s_gpu,wbuf4s_gpu,wbuf5s_gpu,wbuf6s_gpu,dcsidxnc2_gpu,dcsidync2_gpu,detadxnc2_gpu,detadync2_gpu)
      do k = 1,ng
        do j = 1,ny
          do i = 1,nx
            do m=1,nv
              wbuf5s_gpu(i,j,k,m) = w_gpu(i,j,k,m)
              wbuf6s_gpu(i,j,k,m) = w_gpu(i,j,nz-ng+k,m)
            enddo
          enddo
        enddo
      enddo
    endif

  endsubroutine bcswap_c2xybc_step_1_kernel



  subroutine bcswap_corner_step_1_kernel(nx,ny,nz,nv,ng,w_gpu,wbuf1s_c_gpu,wbuf2s_c_gpu,wbuf3s_c_gpu,wbuf4s_c_gpu)
    integer :: nx, ny, nz, ng, nv
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: w_gpu
    real(rkind), dimension(1:,1:,1:,1:) :: wbuf1s_c_gpu, wbuf2s_c_gpu, wbuf3s_c_gpu, wbuf4s_c_gpu
    integer :: i,j,k,m,iercuda

    !$omp target teams distribute parallel do collapse(2) has_device_addr(w_gpu,wbuf1s_c_gpu,wbuf2s_c_gpu,wbuf3s_c_gpu,wbuf4s_c_gpu)
    do j = 1,ny
      do m = 1,nv
        do k=1,ng
          do i=1,ng
            wbuf1s_c_gpu(m,i,j,k) = w_gpu(i,j,k,m)
            wbuf2s_c_gpu(m,i,j,k) = w_gpu(nx-ng+i,j,nz-ng+k,m)
            wbuf3s_c_gpu(m,i,j,k) = w_gpu(i,j,nz-ng+k,m)
            wbuf4s_c_gpu(m,i,j,k) = w_gpu(nx-ng+i,j,k,m)
          enddo
        enddo
      enddo
    enddo


  endsubroutine bcswap_corner_step_1_kernel



  subroutine bcswap_corner_step_1b_kernel(nx,ny,nz,nv,ng,w_gpu,wbuf_lxly_s_gpu,wbuf_lxry_s_gpu,&
  &wbuf_rxly_s_gpu,wbuf_rxry_s_gpu,wbuf_lylz_s_gpu,wbuf_lyrz_s_gpu,wbuf_rylz_s_gpu,wbuf_ryrz_s_gpu,&
  &wbuf_lxlylz_s_gpu,wbuf_lxlyrz_s_gpu,wbuf_lxrylz_s_gpu,wbuf_lxryrz_s_gpu,wbuf_rxlylz_s_gpu,&
  &wbuf_rxlyrz_s_gpu,wbuf_rxrylz_s_gpu,wbuf_rxryrz_s_gpu)
    integer :: nx, ny, nz, ng, nv
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: w_gpu
    real(rkind), dimension(1:,1:,1:,1:) :: wbuf_lxly_s_gpu, wbuf_lxry_s_gpu, wbuf_rxly_s_gpu, wbuf_rxry_s_gpu
    real(rkind), dimension(1:,1:,1:,1:) :: wbuf_lylz_s_gpu, wbuf_lyrz_s_gpu, wbuf_rylz_s_gpu, wbuf_ryrz_s_gpu
    real(rkind), dimension(1:,1:,1:,1:) :: wbuf_lxlylz_s_gpu, wbuf_lxlyrz_s_gpu, wbuf_lxrylz_s_gpu, wbuf_lxryrz_s_gpu
    real(rkind), dimension(1:,1:,1:,1:) :: wbuf_rxlylz_s_gpu, wbuf_rxlyrz_s_gpu, wbuf_rxrylz_s_gpu, wbuf_rxryrz_s_gpu
    integer :: i,j,k,m,iercuda

    !$omp target teams distribute parallel do collapse(2) has_device_addr(w_gpu,wbuf_lxly_s_gpu,wbuf_lxry_s_gpu,wbuf_rxly_s_gpu,wbuf_rxry_s_gpu,wbuf_lylz_s_gpu,wbuf_lyrz_s_gpu,wbuf_rylz_s_gpu,wbuf_ryrz_s_gpu,wbuf_lxlylz_s_gpu,wbuf_lxlyrz_s_gpu,wbuf_lxrylz_s_gpu,wbuf_lxryrz_s_gpu,wbuf_rxlylz_s_gpu,wbuf_rxlyrz_s_gpu,wbuf_rxrylz_s_gpu,wbuf_rxryrz_s_gpu)
    do k = 1,nz
      do m = 1,nv
        do j=1,ng
          do i=1,ng
            wbuf_lxly_s_gpu(m,i,j,k) = w_gpu(i,j,k,m)
            wbuf_lxry_s_gpu(m,i,j,k) = w_gpu(i,ny-ng+j,k,m)
            wbuf_rxly_s_gpu(m,i,j,k) = w_gpu(nx-ng+i,j,k,m)
            wbuf_rxry_s_gpu(m,i,j,k) = w_gpu(nx-ng+i,ny-ng+j,k,m)
          enddo
        enddo
      enddo
    enddo

    !$omp target teams distribute parallel do collapse(2) has_device_addr(w_gpu,wbuf_lxly_s_gpu,wbuf_lxry_s_gpu,wbuf_rxly_s_gpu,wbuf_rxry_s_gpu,wbuf_lylz_s_gpu,wbuf_lyrz_s_gpu,wbuf_rylz_s_gpu,wbuf_ryrz_s_gpu,wbuf_lxlylz_s_gpu,wbuf_lxlyrz_s_gpu,wbuf_lxrylz_s_gpu,wbuf_lxryrz_s_gpu,wbuf_rxlylz_s_gpu,wbuf_rxlyrz_s_gpu,wbuf_rxrylz_s_gpu,wbuf_rxryrz_s_gpu)
    do i = 1,nx
      do m = 1,nv
        do k=1,ng
          do j=1,ng
            wbuf_lylz_s_gpu(m,i,j,k) = w_gpu(i,j,k,m)
            wbuf_lyrz_s_gpu(m,i,j,k) = w_gpu(i,j,nz-ng+k,m)
            wbuf_rylz_s_gpu(m,i,j,k) = w_gpu(i,ny-ng+j,k,m)
            wbuf_ryrz_s_gpu(m,i,j,k) = w_gpu(i,ny-ng+j,nz-ng+k,m)
          enddo
        enddo
      enddo
    enddo

    !$omp target teams distribute parallel do collapse(2) has_device_addr(w_gpu,wbuf_lxly_s_gpu,wbuf_lxry_s_gpu,wbuf_rxly_s_gpu,wbuf_rxry_s_gpu,wbuf_lylz_s_gpu,wbuf_lyrz_s_gpu,wbuf_rylz_s_gpu,wbuf_ryrz_s_gpu,wbuf_lxlylz_s_gpu,wbuf_lxlyrz_s_gpu,wbuf_lxrylz_s_gpu,wbuf_lxryrz_s_gpu,wbuf_rxlylz_s_gpu,wbuf_rxlyrz_s_gpu,wbuf_rxrylz_s_gpu,wbuf_rxryrz_s_gpu)
    do k = 1,ng
      do m = 1,nv
        do j=1,ng
          do i=1,ng
            wbuf_lxlylz_s_gpu(m,i,j,k) = w_gpu(i,j,k,m)
            wbuf_lxlyrz_s_gpu(m,i,j,k) = w_gpu(i,j,nz-ng+k,m)
            wbuf_lxrylz_s_gpu(m,i,j,k) = w_gpu(i,ny-ng+j,k,m)
            wbuf_lxryrz_s_gpu(m,i,j,k) = w_gpu(i,ny-ng+j,nz-ng+k,m)
            wbuf_rxlylz_s_gpu(m,i,j,k) = w_gpu(nx-ng+i,j,k,m)
            wbuf_rxlyrz_s_gpu(m,i,j,k) = w_gpu(nx-ng+i,j,nz-ng+k,m)
            wbuf_rxrylz_s_gpu(m,i,j,k) = w_gpu(nx-ng+i,ny-ng+j,k,m)
            wbuf_rxryrz_s_gpu(m,i,j,k) = w_gpu(nx-ng+i,ny-ng+j,nz-ng+k,m)
          enddo
        enddo
      enddo
    enddo


  endsubroutine bcswap_corner_step_1b_kernel



  subroutine bcswap_step_3_kernel(nx,ny,nz,nv,ng,ndim,w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,&
  &wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,ileftx,ilefty,ileftz,irightx,irighty,irightz)
    integer :: nx, ny, nz, ng, nv, ndim
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: w_gpu
    real(rkind), dimension(1:,1:,1:,1:) :: wbuf1r_gpu, wbuf2r_gpu, wbuf3r_gpu, wbuf4r_gpu, wbuf5r_gpu, wbuf6r_gpu
    integer :: i,j,k,m, iercuda
    integer :: ileftx, ilefty, ileftz, irightx, irighty, irightz
    if (ileftx/=mpi_proc_null) then
      !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu)
      do k = 1,nz
        do j = 1,ny
          do i = 1,ng
            do m=1,nv
              w_gpu(i-ng,j,k,m) = wbuf1r_gpu(i,j,k,m)
            enddo
          enddo
        enddo
      enddo
    endif
    if (irightx/=mpi_proc_null) then
      !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu)
      do k = 1,nz
        do j = 1,ny
          do i = 1,ng
            do m=1,nv
              w_gpu(nx+i,j,k,m) = wbuf2r_gpu(i,j,k,m)
            enddo
          enddo
        enddo
      enddo
    endif
    if (ilefty/=mpi_proc_null) then
      !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu)
      do k = 1,nz
        do j = 1,ng
          do i = 1,nx
            do m=1,nv
              w_gpu(i,j-ng,k,m) = wbuf3r_gpu(i,j,k,m)
            enddo
          enddo
        enddo
      enddo
    endif
    if (irighty/=mpi_proc_null) then
      !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu)
      do k = 1,nz
        do j = 1,ng
          do i = 1,nx
            do m=1,nv
              w_gpu(i,ny+j,k,m) = wbuf4r_gpu(i,j,k,m)
            enddo
          enddo
        enddo
      enddo
    endif
    if(ndim == 3) then
      if (ileftz/=mpi_proc_null) then
        !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu)
        do k = 1,ng
          do j = 1,ny
            do i = 1,nx
              do m=1,nv
                w_gpu(i,j,k-ng,m) = wbuf5r_gpu(i,j,k,m)
              enddo
            enddo
          enddo
        enddo
      endif
      if (irightz/=mpi_proc_null) then
        !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu)
        do k = 1,ng
          do j = 1,ny
            do i = 1,nx
              do m=1,nv
                w_gpu(i,j,nz+k,m) = wbuf6r_gpu(i,j,k,m)
              enddo
            enddo
          enddo
        enddo
      endif
    endif

  endsubroutine bcswap_step_3_kernel



  subroutine bcswap_c2_step_3_kernel(nx,ny,nz,nv,ng,ndim,w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,&
  &wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,ileftx,ilefty,ileftz,irightx,irighty,irightz,dxdcsinc2_gpu,&
  &dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu,is_periodic)
    integer :: nx, ny, nz, ng, nv, ndim
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: w_gpu
    real(rkind), dimension(1:,1:,1:,1:) :: wbuf1r_gpu, wbuf2r_gpu, wbuf3r_gpu, wbuf4r_gpu, wbuf5r_gpu, wbuf6r_gpu
    real(rkind), dimension(1-ng:,1-ng:) :: dxdcsinc2_gpu, dydcsinc2_gpu, dxdetanc2_gpu, dydetanc2_gpu
    integer :: i,j,k,m, iercuda
    integer :: ileftx, ilefty, ileftz, irightx, irighty, irightz
    integer, dimension(3) :: is_periodic
    if (is_periodic(1) == 1) then
      if (ileftx/=mpi_proc_null) then
        !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu)
        do k = 1,nz
          do j = 1,ny
            do i = 1,ng
              w_gpu(i-ng,j,k,1) = wbuf1r_gpu(i,j,k,1)
              w_gpu(i-ng,j,k,2) = wbuf1r_gpu(i,j,k,2)*dxdcsinc2_gpu(i-ng,j)+wbuf1r_gpu(i,j,k,3)*dxdetanc2_gpu(i-ng,j)
              w_gpu(i-ng,j,k,3) = wbuf1r_gpu(i,j,k,2)*dydcsinc2_gpu(i-ng,j)+wbuf1r_gpu(i,j,k,3)*dydetanc2_gpu(i-ng,j)
              w_gpu(i-ng,j,k,4) = wbuf1r_gpu(i,j,k,4)
              w_gpu(i-ng,j,k,5) = wbuf1r_gpu(i,j,k,5)
            enddo
          enddo
        enddo
      endif
      if (irightx/=mpi_proc_null) then
        !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu)
        do k = 1,nz
          do j = 1,ny
            do i = 1,ng
              w_gpu(nx+i,j,k,1) = wbuf2r_gpu(i,j,k,1)
              w_gpu(nx+i,j,k,2) = wbuf2r_gpu(i,j,k,2)*dxdcsinc2_gpu(nx+i,j)+wbuf2r_gpu(i,j,k,3)*dxdetanc2_gpu(nx+i,j)
              w_gpu(nx+i,j,k,3) = wbuf2r_gpu(i,j,k,2)*dydcsinc2_gpu(nx+i,j)+wbuf2r_gpu(i,j,k,3)*dydetanc2_gpu(nx+i,j)
              w_gpu(nx+i,j,k,4) = wbuf2r_gpu(i,j,k,4)
              w_gpu(nx+i,j,k,5) = wbuf2r_gpu(i,j,k,5)
            enddo
          enddo
        enddo
      endif
    else
      if (ileftx/=mpi_proc_null) then
        !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu)
        do k = 1,nz
          do j = 1,ny
            do i = 1,ng
              do m=1,nv
                w_gpu(i-ng,j,k,m) = wbuf1r_gpu(i,j,k,m)
              enddo
            enddo
          enddo
        enddo
      endif
      if (irightx/=mpi_proc_null) then
        !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu)
        do k = 1,nz
          do j = 1,ny
            do i = 1,ng
              do m=1,nv
                w_gpu(nx+i,j,k,m) = wbuf2r_gpu(i,j,k,m)
              enddo
            enddo
          enddo
        enddo
      endif
    endif
    if (ilefty/=mpi_proc_null) then
      !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu)
      do k = 1,nz
        do j = 1,ng
          do i = 1,nx
            do m=1,nv
              w_gpu(i,j-ng,k,m) = wbuf3r_gpu(i,j,k,m)
            enddo
          enddo
        enddo
      enddo
    endif
    if (irighty/=mpi_proc_null) then
      !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu)
      do k = 1,nz
        do j = 1,ng
          do i = 1,nx
            do m=1,nv
              w_gpu(i,ny+j,k,m) = wbuf4r_gpu(i,j,k,m)
            enddo
          enddo
        enddo
      enddo
    endif
    if(ndim == 3) then
      if (ileftz/=mpi_proc_null) then
        !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu)
        do k = 1,ng
          do j = 1,ny
            do i = 1,nx
              do m=1,nv
                w_gpu(i,j,k-ng,m) = wbuf5r_gpu(i,j,k,m)
              enddo
            enddo
          enddo
        enddo
      endif
      if (irightz/=mpi_proc_null) then
        !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu)
        do k = 1,ng
          do j = 1,ny
            do i = 1,nx
              do m=1,nv
                w_gpu(i,j,nz+k,m) = wbuf6r_gpu(i,j,k,m)
              enddo
            enddo
          enddo
        enddo
      endif
    endif

  endsubroutine bcswap_c2_step_3_kernel



  subroutine bcswap_c2xybc_step_3_kernel(nx,ny,nz,nv,ng,ndim,w_gpu,wbuf1xybcr_gpu,wbuf2xybcr_gpu,&
  &wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,ileftx,ilefty,ileftz,irightx,irighty,irightz,&
  &dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu,is_periodic)
    integer :: nx, ny, nz, ng, nv, ndim
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: w_gpu
    real(rkind), dimension(1:,1-ng:,1:,1:) :: wbuf1xybcr_gpu, wbuf2xybcr_gpu
    real(rkind), dimension(1:,1:,1:,1:) :: wbuf3r_gpu, wbuf4r_gpu, wbuf5r_gpu, wbuf6r_gpu
    real(rkind), dimension(1-ng:,1-ng:) :: dxdcsinc2_gpu, dydcsinc2_gpu, dxdetanc2_gpu, dydetanc2_gpu
    integer :: i,j,k,m, iercuda
    integer :: ileftx, ilefty, ileftz, irightx, irighty, irightz
    integer, dimension(3) :: is_periodic
    if (is_periodic(1) == 1) then
      if (ileftx/=mpi_proc_null) then
        !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu,wbuf1xybcr_gpu,wbuf2xybcr_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu)
        do k = 1,nz
          do j = 1-ng,ny+ng
            do i = 1,ng
              w_gpu(i-ng,j,k,1) = wbuf1xybcr_gpu(i,j,k,1)
              w_gpu(i-ng,j,k,2) = wbuf1xybcr_gpu(i,j,k,2)*dxdcsinc2_gpu(i-ng,j)+wbuf1xybcr_gpu(i,j,k,3)*dxdetanc2_gpu(i-ng,j)
              w_gpu(i-ng,j,k,3) = wbuf1xybcr_gpu(i,j,k,2)*dydcsinc2_gpu(i-ng,j)+wbuf1xybcr_gpu(i,j,k,3)*dydetanc2_gpu(i-ng,j)
              w_gpu(i-ng,j,k,4) = wbuf1xybcr_gpu(i,j,k,4)
              w_gpu(i-ng,j,k,5) = wbuf1xybcr_gpu(i,j,k,5)
            enddo
          enddo
        enddo
      endif
      if (irightx/=mpi_proc_null) then
        !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu,wbuf1xybcr_gpu,wbuf2xybcr_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu)
        do k = 1,nz
          do j = 1-ng,ny+ng
            do i = 1,ng
              w_gpu(nx+i,j,k,1) = wbuf2xybcr_gpu(i,j,k,1)
              w_gpu(nx+i,j,k,2) = wbuf2xybcr_gpu(i,j,k,2)*dxdcsinc2_gpu(nx+i,j)+wbuf2xybcr_gpu(i,j,k,3)*dxdetanc2_gpu(nx+i,j)
              w_gpu(nx+i,j,k,3) = wbuf2xybcr_gpu(i,j,k,2)*dydcsinc2_gpu(nx+i,j)+wbuf2xybcr_gpu(i,j,k,3)*dydetanc2_gpu(nx+i,j)
              w_gpu(nx+i,j,k,4) = wbuf2xybcr_gpu(i,j,k,4)
              w_gpu(nx+i,j,k,5) = wbuf2xybcr_gpu(i,j,k,5)
            enddo
          enddo
        enddo
      endif
    else
      if (ileftx/=mpi_proc_null) then
        !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu,wbuf1xybcr_gpu,wbuf2xybcr_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu)
        do k = 1,nz
          do j = 1-ng,ny+ng
            do i = 1,ng
              do m=1,nv
                w_gpu(i-ng,j,k,m) = wbuf1xybcr_gpu(i,j,k,m)
              enddo
            enddo
          enddo
        enddo
      endif
      if (irightx/=mpi_proc_null) then
        !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu,wbuf1xybcr_gpu,wbuf2xybcr_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu)
        do k = 1,nz
          do j = 1-ng,ny+ng
            do i = 1,ng
              do m=1,nv
                w_gpu(nx+i,j,k,m) = wbuf2xybcr_gpu(i,j,k,m)
              enddo
            enddo
          enddo
        enddo
      endif
    endif
    if (ilefty/=mpi_proc_null) then
      !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu,wbuf1xybcr_gpu,wbuf2xybcr_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu)
      do k = 1,nz
        do j = 1,ng
          do i = 1,nx
            do m=1,nv
              w_gpu(i,j-ng,k,m) = wbuf3r_gpu(i,j,k,m)
            enddo
          enddo
        enddo
      enddo
    endif
    if (irighty/=mpi_proc_null) then
      !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu,wbuf1xybcr_gpu,wbuf2xybcr_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu)
      do k = 1,nz
        do j = 1,ng
          do i = 1,nx
            do m=1,nv
              w_gpu(i,ny+j,k,m) = wbuf4r_gpu(i,j,k,m)
            enddo
          enddo
        enddo
      enddo
    endif
    if(ndim == 3) then
      if (ileftz/=mpi_proc_null) then
        !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu,wbuf1xybcr_gpu,wbuf2xybcr_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu)
        do k = 1,ng
          do j = 1,ny
            do i = 1,nx
              do m=1,nv
                w_gpu(i,j,k-ng,m) = wbuf5r_gpu(i,j,k,m)
              enddo
            enddo
          enddo
        enddo
      endif
      if (irightz/=mpi_proc_null) then
        !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu,wbuf1xybcr_gpu,wbuf2xybcr_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu)
        do k = 1,ng
          do j = 1,ny
            do i = 1,nx
              do m=1,nv
                w_gpu(i,j,nz+k,m) = wbuf6r_gpu(i,j,k,m)
              enddo
            enddo
          enddo
        enddo
      endif
    endif

  endsubroutine bcswap_c2xybc_step_3_kernel



  subroutine bcswap_corner_step_3_kernel(nx,ny,nz,nv,ng,w_gpu,wbuf1r_c_gpu,wbuf2r_c_gpu,&
  &wbuf3r_c_gpu,wbuf4r_c_gpu,ileftbottom,ilefttop,irightbottom,irighttop)
    integer :: nx, ny, nz, ng, nv
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: w_gpu
    real(rkind), dimension(1:,1:,1:,1:) :: wbuf1r_c_gpu, wbuf2r_c_gpu, wbuf3r_c_gpu, wbuf4r_c_gpu
    integer :: i,j,k,m, iercuda
    integer :: ileftbottom, ilefttop, irightbottom, irighttop
    if (ileftbottom/=mpi_proc_null) then
      !$omp target teams distribute parallel do collapse(2) has_device_addr(w_gpu,wbuf1r_c_gpu,wbuf2r_c_gpu,wbuf3r_c_gpu,wbuf4r_c_gpu)
      do j = 1,ny
        do m = 1,nv
          do k=1,ng
            do i=1,ng
              w_gpu(i-ng,j,k-ng,m) = wbuf1r_c_gpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif
    if (irighttop/=mpi_proc_null) then
      !$omp target teams distribute parallel do collapse(2) has_device_addr(w_gpu,wbuf1r_c_gpu,wbuf2r_c_gpu,wbuf3r_c_gpu,wbuf4r_c_gpu)
      do j = 1,ny
        do m = 1,nv
          do k=1,ng
            do i=1,ng
              w_gpu(i+nx,j,k+nz,m) = wbuf2r_c_gpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif
    if (ilefttop/=mpi_proc_null) then
      !$omp target teams distribute parallel do collapse(2) has_device_addr(w_gpu,wbuf1r_c_gpu,wbuf2r_c_gpu,wbuf3r_c_gpu,wbuf4r_c_gpu)
      do j = 1,ny
        do m = 1,nv
          do k=1,ng
            do i=1,ng
              w_gpu(i-ng,j,k+nz,m) = wbuf3r_c_gpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif
    if (irightbottom/=mpi_proc_null) then
      !$omp target teams distribute parallel do collapse(2) has_device_addr(w_gpu,wbuf1r_c_gpu,wbuf2r_c_gpu,wbuf3r_c_gpu,wbuf4r_c_gpu)
      do j = 1,ny
        do m = 1,nv
          do k=1,ng
            do i=1,ng
              w_gpu(i+nx,j,k-ng,m) = wbuf4r_c_gpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif

  endsubroutine bcswap_corner_step_3_kernel



  subroutine bcswap_corner_step_3b_kernel(nx,ny,nz,nv,ng,w_gpu,wbuf_lxly_r_gpu,wbuf_lxry_r_gpu,&
  &wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,wbuf_ryrz_r_gpu,&
  &wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,wbuf_rxlylz_r_gpu,&
  &wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu,lxly,lxry,rxly,rxry,lylz,lyrz,rylz,ryrz,&
  &lxlylz,lxlyrz,lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz)
    integer :: nx, ny, nz, ng, nv
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: w_gpu
    real(rkind), dimension(1:,1:,1:,1:) :: wbuf_lxly_r_gpu, wbuf_lxry_r_gpu, wbuf_rxly_r_gpu, wbuf_rxry_r_gpu
    real(rkind), dimension(1:,1:,1:,1:) :: wbuf_lylz_r_gpu, wbuf_lyrz_r_gpu, wbuf_rylz_r_gpu, wbuf_ryrz_r_gpu
    real(rkind), dimension(1:,1:,1:,1:) :: wbuf_lxlylz_r_gpu, wbuf_lxlyrz_r_gpu, wbuf_lxrylz_r_gpu, wbuf_lxryrz_r_gpu
    real(rkind), dimension(1:,1:,1:,1:) :: wbuf_rxlylz_r_gpu, wbuf_rxlyrz_r_gpu, wbuf_rxrylz_r_gpu, wbuf_rxryrz_r_gpu
    integer :: i,j,k,m, iercuda
    integer :: lxly, lxry, rxly, rxry
    integer :: lylz, lyrz, rylz, ryrz
    integer :: lxlylz, lxlyrz, lxrylz, lxryrz
    integer :: rxlylz, rxlyrz, rxrylz, rxryrz
    if (lxly/=mpi_proc_null) then
      !$omp target teams distribute parallel do collapse(2) has_device_addr(w_gpu,wbuf_lxly_r_gpu,wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu)
      do k = 1,nz
        do m = 1,nv
          do i=1,ng
            do j=1,ng
              w_gpu(i-ng,j-ng,k,m) = wbuf_lxly_r_gpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif
    if (rxry/=mpi_proc_null) then
      !$omp target teams distribute parallel do collapse(2) has_device_addr(w_gpu,wbuf_lxly_r_gpu,wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu)
      do k = 1,nz
        do m = 1,nv
          do i=1,ng
            do j=1,ng
              w_gpu(i+nx,j+ny,k,m) = wbuf_rxry_r_gpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif
    if (lxry/=mpi_proc_null) then
      !$omp target teams distribute parallel do collapse(2) has_device_addr(w_gpu,wbuf_lxly_r_gpu,wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu)
      do k = 1,nz
        do m = 1,nv
          do i=1,ng
            do j=1,ng
              w_gpu(i-ng,j+ny,k,m) = wbuf_lxry_r_gpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif
    if (rxly/=mpi_proc_null) then
      !$omp target teams distribute parallel do collapse(2) has_device_addr(w_gpu,wbuf_lxly_r_gpu,wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu)
      do k = 1,nz
        do m = 1,nv
          do i=1,ng
            do j=1,ng
              w_gpu(i+nx,j-ng,k,m) = wbuf_rxly_r_gpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif

    if (lylz/=mpi_proc_null) then
      !$omp target teams distribute parallel do collapse(2) has_device_addr(w_gpu,wbuf_lxly_r_gpu,wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu)
      do i = 1,nx
        do m = 1,nv
          do k=1,ng
            do j=1,ng
              w_gpu(i,j-ng,k-ng,m) = wbuf_lylz_r_gpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif
    if (ryrz/=mpi_proc_null) then
      !$omp target teams distribute parallel do collapse(2) has_device_addr(w_gpu,wbuf_lxly_r_gpu,wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu)
      do i = 1,nx
        do m = 1,nv
          do k=1,ng
            do j=1,ng
              w_gpu(i,j+ny,k+nz,m) = wbuf_ryrz_r_gpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif
    if (lyrz/=mpi_proc_null) then
      !$omp target teams distribute parallel do collapse(2) has_device_addr(w_gpu,wbuf_lxly_r_gpu,wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu)
      do i = 1,nx
        do m = 1,nv
          do k=1,ng
            do j=1,ng
              w_gpu(i,j-ng,k+nz,m) = wbuf_lyrz_r_gpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif
    if (rylz/=mpi_proc_null) then
      !$omp target teams distribute parallel do collapse(2) has_device_addr(w_gpu,wbuf_lxly_r_gpu,wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu)
      do i = 1,nx
        do m = 1,nv
          do k=1,ng
            do j=1,ng
              w_gpu(i,j+ny,k-ng,m) = wbuf_rylz_r_gpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif

    if (lxlylz/=mpi_proc_null) then
      !$omp target teams distribute parallel do collapse(2) has_device_addr(w_gpu,wbuf_lxly_r_gpu,wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu)
      do k = 1,ng
        do m = 1,nv
          do i=1,ng
            do j=1,ng
              w_gpu(i-ng,j-ng,k-ng,m) = wbuf_lxlylz_r_gpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif
    if (lxlyrz/=mpi_proc_null) then
      !$omp target teams distribute parallel do collapse(2) has_device_addr(w_gpu,wbuf_lxly_r_gpu,wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu)
      do k = 1,ng
        do m = 1,nv
          do i=1,ng
            do j=1,ng
              w_gpu(i-ng,j-ng,k+nz,m) = wbuf_lxlyrz_r_gpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif
    if (lxrylz/=mpi_proc_null) then
      !$omp target teams distribute parallel do collapse(2) has_device_addr(w_gpu,wbuf_lxly_r_gpu,wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu)
      do k = 1,ng
        do m = 1,nv
          do i=1,ng
            do j=1,ng
              w_gpu(i-ng,j+ny,k-ng,m) = wbuf_lxrylz_r_gpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif
    if (lxryrz/=mpi_proc_null) then
      !$omp target teams distribute parallel do collapse(2) has_device_addr(w_gpu,wbuf_lxly_r_gpu,wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu)
      do k = 1,ng
        do m = 1,nv
          do i=1,ng
            do j=1,ng
              w_gpu(i-ng,j+ny,k+nz,m) = wbuf_lxryrz_r_gpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif
    if (rxlylz/=mpi_proc_null) then
      !$omp target teams distribute parallel do collapse(2) has_device_addr(w_gpu,wbuf_lxly_r_gpu,wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu)
      do k = 1,ng
        do m = 1,nv
          do i=1,ng
            do j=1,ng
              w_gpu(i+nx,j-ng,k-ng,m) = wbuf_rxlylz_r_gpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif
    if (rxlyrz/=mpi_proc_null) then
      !$omp target teams distribute parallel do collapse(2) has_device_addr(w_gpu,wbuf_lxly_r_gpu,wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu)
      do k = 1,ng
        do m = 1,nv
          do i=1,ng
            do j=1,ng
              w_gpu(i+nx,j-ng,k+nz,m) = wbuf_rxlyrz_r_gpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif
    if (rxrylz/=mpi_proc_null) then
      !$omp target teams distribute parallel do collapse(2) has_device_addr(w_gpu,wbuf_lxly_r_gpu,wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu)
      do k = 1,ng
        do m = 1,nv
          do i=1,ng
            do j=1,ng
              w_gpu(i+nx,j+ny,k-ng,m) = wbuf_rxrylz_r_gpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif
    if (rxryrz/=mpi_proc_null) then
      !$omp target teams distribute parallel do collapse(2) has_device_addr(w_gpu,wbuf_lxly_r_gpu,wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu)
      do k = 1,ng
        do m = 1,nv
          do i=1,ng
            do j=1,ng
              w_gpu(i+nx,j+ny,k+nz,m) = wbuf_rxryrz_r_gpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif

  endsubroutine bcswap_corner_step_3b_kernel



  subroutine bcswap_corner(self, steps)
    class(base_omp_object), intent(inout) :: self
    logical, dimension(3), optional :: steps
    logical, dimension(3) :: steps_
    integer :: iercuda, indc
    integer, dimension(mpi_status_size) :: istatus

    steps_ = .true. ; if(present(steps)) steps_ = steps

    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv,&
    & w_gpu => self%w_gpu, wbuf1s_c_gpu => self%wbuf1s_c_gpu, wbuf2s_c_gpu => self%wbuf2s_c_gpu,&
    & wbuf3s_c_gpu => self%wbuf3s_c_gpu, wbuf4s_c_gpu => self%wbuf4s_c_gpu, wbuf1r_c_gpu => self%wbuf1r_c&
    &_gpu, wbuf2r_c_gpu => self%wbuf2r_c_gpu, wbuf3r_c_gpu => self%wbuf3r_c_gpu, wbuf4r_c_gpu => self%wbu&
    &f4r_c_gpu, ileftbottom => self%field%ileftbottom, irightbottom => self%field%irightbottom,&
    & ilefttop => self%field%ilefttop, irighttop => self%field%irighttop, mp_cart => self%field%mp_cart,&
    & iermpi => self%mpi_err, nrank => self%field%myrank)

      if(steps_(1)) then
        call bcswap_corner_step_1_kernel(nx, ny, nz, nv, ng, w_gpu, wbuf1s_c_gpu, wbuf2s_c_gpu, wbuf3s_c_gpu, wbuf4s_c_gpu)
      endif
      if(steps_(2)) then
        indc = nv*ny*ng*ng
        if (ileftbottom == nrank) then
          self%ierr = omp_target_memcpy_f(wbuf2r_c_gpu , wbuf1s_c_gpu,0,0,self%mydev,self%mydev)
          self%ierr = omp_target_memcpy_f(wbuf1r_c_gpu , wbuf2s_c_gpu,0,0,self%mydev,self%mydev)
        else
          call mpi_sendrecv(wbuf1s_c_gpu,indc,mpi_prec,ileftbottom ,1, wbuf2r_c_gpu,indc,mpi_prec,&
          &irighttop ,1,mp_cart,istatus,iermpi)
          call mpi_sendrecv(wbuf2s_c_gpu,indc,mpi_prec,irighttop ,2, wbuf1r_c_gpu,indc,mpi_prec,&
          &ileftbottom ,2,mp_cart,istatus,iermpi)
        endif
        if (ilefttop == nrank) then
          self%ierr = omp_target_memcpy_f(wbuf4r_c_gpu , wbuf3s_c_gpu,0,0,self%mydev,self%mydev)
          self%ierr = omp_target_memcpy_f(wbuf3r_c_gpu , wbuf4s_c_gpu,0,0,self%mydev,self%mydev)
        else
          call mpi_sendrecv(wbuf3s_c_gpu,indc,mpi_prec,ilefttop ,3, wbuf4r_c_gpu,indc,mpi_prec,&
          &irightbottom,3,mp_cart,istatus,iermpi)
          call mpi_sendrecv(wbuf4s_c_gpu,indc,mpi_prec,irightbottom,4, wbuf3r_c_gpu,indc,mpi_prec,&
          &ilefttop ,4,mp_cart,istatus,iermpi)
        endif
      endif
      if(steps_(3)) then
        call bcswap_corner_step_3_kernel(nx, ny, nz, nv, ng, w_gpu, wbuf1r_c_gpu, wbuf2r_c_gpu,&
        & wbuf3r_c_gpu, wbuf4r_c_gpu, ileftbottom, ilefttop, irightbottom, irighttop)
      endif
    endassociate
  endsubroutine bcswap_corner

  subroutine bcswap_corner_var(self, w_swap, steps)
    class(base_omp_object), intent(inout) :: self
    real(rkind), dimension(1-self%ng:self%nx+self%ng, 1-self%ng:self%ny+self%ng, 1-self%ng:self%nz+self%ng,1:1) :: w_swap
    logical, dimension(3), optional :: steps
    logical, dimension(3) :: steps_
    integer :: iercuda, indc
    integer, dimension(mpi_status_size) :: istatus

    steps_ = .true. ; if(present(steps)) steps_ = steps

    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv,&
    & wbuf1s_c_gpu => self%wbuf1s_c_gpu, wbuf2s_c_gpu => self%wbuf2s_c_gpu, wbuf3s_c_gpu => self%wbuf3s_c&
    &_gpu, wbuf4s_c_gpu => self%wbuf4s_c_gpu, wbuf1r_c_gpu => self%wbuf1r_c_gpu, wbuf2r_c_gpu => self%wbu&
    &f2r_c_gpu, wbuf3r_c_gpu => self%wbuf3r_c_gpu, wbuf4r_c_gpu => self%wbuf4r_c_gpu, ileftbottom => self&
    &%field%ileftbottom, irightbottom => self%field%irightbottom, ilefttop => self%field%ilefttop,&
    & irighttop => self%field%irighttop, mp_cart => self%field%mp_cart, iermpi => self%mpi_err,&
    & nrank => self%field%myrank)

      if(steps_(1)) then
        call bcswap_corner_step_1_kernel(nx, ny, nz, 1, ng, w_swap, wbuf1s_c_gpu, wbuf2s_c_gpu, wbuf3s_c_gpu, wbuf4s_c_gpu)
      endif
      if(steps_(2)) then
        indc = ny*ng*ng
        if (ileftbottom == nrank) then
          self%ierr = omp_target_memcpy_f(wbuf2r_c_gpu , wbuf1s_c_gpu,0,0,self%mydev,self%mydev)
          self%ierr = omp_target_memcpy_f(wbuf1r_c_gpu , wbuf2s_c_gpu,0,0,self%mydev,self%mydev)
        else
          call mpi_sendrecv(wbuf1s_c_gpu,indc,mpi_prec,ileftbottom ,1, wbuf2r_c_gpu,indc,mpi_prec,&
          &irighttop ,1,mp_cart,istatus,iermpi)
          call mpi_sendrecv(wbuf2s_c_gpu,indc,mpi_prec,irighttop ,2, wbuf1r_c_gpu,indc,mpi_prec,&
          &ileftbottom ,2,mp_cart,istatus,iermpi)
        endif
        if (ilefttop == nrank) then
          self%ierr = omp_target_memcpy_f(wbuf4r_c_gpu , wbuf3s_c_gpu,0,0,self%mydev,self%mydev)
          self%ierr = omp_target_memcpy_f(wbuf3r_c_gpu , wbuf4s_c_gpu,0,0,self%mydev,self%mydev)
        else
          call mpi_sendrecv(wbuf3s_c_gpu,indc,mpi_prec,ilefttop ,3, wbuf4r_c_gpu,indc,mpi_prec,&
          &irightbottom,3,mp_cart,istatus,iermpi)
          call mpi_sendrecv(wbuf4s_c_gpu,indc,mpi_prec,irightbottom,4, wbuf3r_c_gpu,indc,mpi_prec,&
          &ilefttop ,4,mp_cart,istatus,iermpi)
        endif
      endif
      if(steps_(3)) then
        call bcswap_corner_step_3_kernel(nx, ny, nz, 1, ng, w_swap, wbuf1r_c_gpu, wbuf2r_c_gpu,&
        & wbuf3r_c_gpu, wbuf4r_c_gpu, ileftbottom, ilefttop, irightbottom, irighttop)
      endif
    endassociate
  endsubroutine bcswap_corner_var

  subroutine bcswap(self, steps, swap_xy_corner_bc)
    class(base_omp_object), intent(inout) :: self
    logical, dimension(3), optional :: steps
    integer, optional :: swap_xy_corner_bc
    integer :: swap_xy_corner_bc_
    logical, dimension(3) :: steps_
    integer :: iercuda, indx, indy, indz
    integer, dimension(mpi_status_size) :: istatus
    integer, dimension(3) :: is_periodic

    steps_ = .true. ; if(present(steps)) steps_ = steps
    swap_xy_corner_bc_ = 1 ; if(present(swap_xy_corner_bc)) swap_xy_corner_bc_ = swap_xy_corner_bc

    is_periodic(:) = 0
    if(self%field%grid%is_xyz_periodic(1)) is_periodic(1) = 1
    if(self%field%grid%is_xyz_periodic(2)) is_periodic(2) = 1
    if(self%field%grid%is_xyz_periodic(3)) is_periodic(3) = 1

    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv,&
    & w_gpu => self%w_gpu, wbuf1s_gpu => self%wbuf1s_gpu, wbuf2s_gpu => self%wbuf2s_gpu,&
    & wbuf1xybcs_gpu => self%wbuf1xybcs_gpu, wbuf2xybcs_gpu => self%wbuf2xybcs_gpu, wbuf3s_gpu => self%wb&
    &uf3s_gpu, wbuf4s_gpu => self%wbuf4s_gpu, wbuf5s_gpu => self%wbuf5s_gpu, wbuf6s_gpu => self%wbuf6s_gp&
    &u, wbuf1r_gpu => self%wbuf1r_gpu, wbuf2r_gpu => self%wbuf2r_gpu, wbuf1xybcr_gpu => self%wbuf1xybcr_g&
    &pu, wbuf2xybcr_gpu => self%wbuf2xybcr_gpu, wbuf3r_gpu => self%wbuf3r_gpu, wbuf4r_gpu => self%wbuf4r_&
    &gpu, wbuf5r_gpu => self%wbuf5r_gpu, wbuf6r_gpu => self%wbuf6r_gpu, ileftx => self%field%ileftx,&
    &irightx => self%field%irightx, nrank_x => self%field%nrank_x, ilefty => self%field%ilefty,&
    &irighty => self%field%irighty, nrank_y => self%field%nrank_y, ileftz => self%field%ileftz,&
    &irightz => self%field%irightz, nrank_z => self%field%nrank_z, dcsidxnc2_gpu => self%dcsidxnc2_gpu,&
    & dcsidync2_gpu => self%dcsidync2_gpu, detadxnc2_gpu => self%detadxnc2_gpu, detadync2_gpu => self%det&
    &adync2_gpu, dxdcsinc2_gpu => self%dxdcsinc2_gpu, dydcsinc2_gpu => self%dydcsinc2_gpu,&
    & dxdetanc2_gpu => self%dxdetanc2_gpu, dydetanc2_gpu => self%dydetanc2_gpu, mp_cartx => self%field%mp&
    &_cartx, mp_carty => self%field%mp_carty, mp_cartz => self%field%mp_cartz, ncoords => self%field%ncoo&
    &rds, nblocks => self%field%nblocks, ndim => self%field%grid%ndim)

      if(steps_(1)) then
        if(self%field%grid%grid_dim == 1) then
          call bcswap_step_1_kernel(nx, ny, nz, nv, ng, ndim, w_gpu, wbuf1s_gpu, wbuf2s_gpu,&
          & wbuf3s_gpu, wbuf4s_gpu, wbuf5s_gpu, wbuf6s_gpu)
        elseif(self%field%grid%grid_dim == 2) then
          if(swap_xy_corner_bc_ == 0) then
            call bcswap_c2_step_1_kernel(nx, ny, nz, nv, ng, ndim, w_gpu, wbuf1s_gpu, wbuf2s_gpu,&
            & wbuf3s_gpu, wbuf4s_gpu, wbuf5s_gpu, wbuf6s_gpu, dcsidxnc2_gpu, dcsidync2_gpu, detadxnc2_gpu,&
            & detadync2_gpu, is_periodic)
          else
            call bcswap_c2xybc_step_1_kernel(nx, ny, nz, nv, ng, ndim, w_gpu, wbuf1xybcs_gpu,&
            & wbuf2xybcs_gpu, wbuf3s_gpu, wbuf4s_gpu, wbuf5s_gpu, wbuf6s_gpu, dcsidxnc2_gpu, dcsidync2_gpu,&
            & detadxnc2_gpu, detadync2_gpu, is_periodic)
          endif
        endif
      endif
      if(steps_(2)) then
        indy = nv*nx*ng*nz
        indz = nv*nx*ny*ng
        if(self%field%grid%grid_dim == 1 .or. (self%field%grid%grid_dim == 2 .and. swap_xy_corner_bc_ == 0)) then
          indx = nv*ng*ny*nz
          if(ileftx == nrank_x) then
            self%ierr = omp_target_memcpy_f(wbuf2r_gpu , wbuf1s_gpu,0,0,self%mydev,self%mydev)
            self%ierr = omp_target_memcpy_f(wbuf1r_gpu , wbuf2s_gpu,0,0,self%mydev,self%mydev)
          else
            call mpi_sendrecv(wbuf1s_gpu,indx,mpi_prec,ileftx ,1,wbuf2r_gpu,indx,mpi_prec,irightx,1,mp_cartx,istatus,self%mpi_err)
            call mpi_sendrecv(wbuf2s_gpu,indx,mpi_prec,irightx,2,wbuf1r_gpu,indx,mpi_prec,ileftx ,2,mp_cartx,istatus,self%mpi_err)
          endif
        else
          indx = nv*ng*(2*ng+ny)*nz
          if(ileftx == nrank_x) then
            self%ierr = omp_target_memcpy_f(wbuf2xybcr_gpu , wbuf1xybcs_gpu,0,0,self%mydev,self%mydev)
            self%ierr = omp_target_memcpy_f(wbuf1xybcr_gpu , wbuf2xybcs_gpu,0,0,self%mydev,self%mydev)
          else
            call mpi_sendrecv(wbuf1xybcs_gpu,indx,mpi_prec,ileftx ,1,wbuf2xybcr_gpu,indx,mpi_prec,&
            &irightx,1,mp_cartx,istatus,self%mpi_err)
            call mpi_sendrecv(wbuf2xybcs_gpu,indx,mpi_prec,irightx,2,wbuf1xybcr_gpu,indx,mpi_prec,&
            &ileftx ,2,mp_cartx,istatus,self%mpi_err)
          endif
        endif
        if(ilefty == nrank_y) then
          self%ierr = omp_target_memcpy_f(wbuf4r_gpu , wbuf3s_gpu,0,0,self%mydev,self%mydev)
          self%ierr = omp_target_memcpy_f(wbuf3r_gpu , wbuf4s_gpu,0,0,self%mydev,self%mydev)
        else
          call mpi_sendrecv(wbuf3s_gpu,indy,mpi_prec,ilefty ,3,wbuf4r_gpu,indy,mpi_prec,irighty,3,mp_carty,istatus,self%mpi_err)
          call mpi_sendrecv(wbuf4s_gpu,indy,mpi_prec,irighty,4,wbuf3r_gpu,indy,mpi_prec,ilefty ,4,mp_carty,istatus,self%mpi_err)
        endif
        if(ndim == 3) then
          if(ileftz == nrank_z) then
            self%ierr = omp_target_memcpy_f(wbuf6r_gpu , wbuf5s_gpu,0,0,self%mydev,self%mydev)
            self%ierr = omp_target_memcpy_f(wbuf5r_gpu , wbuf6s_gpu,0,0,self%mydev,self%mydev)
          else
            call mpi_sendrecv(wbuf5s_gpu,indz,mpi_prec,ileftz ,5,wbuf6r_gpu,indz,mpi_prec,irightz,5,mp_cartz,istatus,self%mpi_err)
            call mpi_sendrecv(wbuf6s_gpu,indz,mpi_prec,irightz,6,wbuf5r_gpu,indz,mpi_prec,ileftz ,6,mp_cartz,istatus,self%mpi_err)
          endif
        endif
      endif
      if(steps_(3)) then
        if(self%field%grid%grid_dim == 1) then
          call bcswap_step_3_kernel(nx, ny, nz, nv, ng, ndim, w_gpu, wbuf1r_gpu, wbuf2r_gpu,&
          & wbuf3r_gpu, wbuf4r_gpu, wbuf5r_gpu, wbuf6r_gpu, ileftx, ilefty, ileftz, irightx, irighty, irightz)
        elseif(self%field%grid%grid_dim == 2) then
          if(swap_xy_corner_bc_ == 0) then
            call bcswap_c2_step_3_kernel(nx, ny, nz, nv, ng, ndim, w_gpu, wbuf1r_gpu, wbuf2r_gpu,&
            & wbuf3r_gpu, wbuf4r_gpu, wbuf5r_gpu, wbuf6r_gpu, ileftx, ilefty, ileftz, irightx, irighty, irightz,&
            & dxdcsinc2_gpu, dydcsinc2_gpu, dxdetanc2_gpu, dydetanc2_gpu, is_periodic)
          else
            call bcswap_c2xybc_step_3_kernel(nx, ny, nz, nv, ng, ndim, w_gpu, wbuf1xybcr_gpu,&
            & wbuf2xybcr_gpu, wbuf3r_gpu, wbuf4r_gpu, wbuf5r_gpu, wbuf6r_gpu, ileftx, ilefty, ileftz, irightx,&
            & irighty, irightz, dxdcsinc2_gpu, dydcsinc2_gpu, dxdetanc2_gpu, dydetanc2_gpu, is_periodic)
            if(is_periodic(1) == 0) call extr_corner_ymin_kernel(ncoords, nblocks, nx, ny, nz, ng, nv, w_gpu)
          endif
        endif
      endif
    endassociate
  endsubroutine bcswap

  subroutine extr_corner_ymin_kernel(ncoords,nblocks,nx,ny,nz,ng,nv,w_gpu)
    integer, intent(in) :: nx, ny, nz, ng, nv
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(inout) :: w_gpu
    integer, dimension(3) :: ncoords, nblocks
    integer :: i,j,k,m,iercuda
    if(ncoords(1) == 0) then
      !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu)
      do k = 1,nz
        do j = 1,ng
          do i = 1,ng
            do m=1,nv
              w_gpu(1-i,1-j,k,m) = w_gpu(1-i,1,k,m)
            enddo
          enddo
        enddo
      enddo
    endif
    if(ncoords(1) == nblocks(1)-1) then
      !$omp target teams distribute parallel do collapse(3) has_device_addr(w_gpu)
      do k = 1,nz
        do j = 1,ng
          do i = 1,ng
            do m=1,nv
              w_gpu(nx+i,1-j,k,m) = w_gpu(nx+i,1,k,m)
            enddo
          enddo
        enddo
      enddo
    endif

  endsubroutine extr_corner_ymin_kernel



  subroutine bcswap_wake(self, steps)
    class(base_omp_object), intent(inout) :: self
    logical, dimension(3), optional :: steps
    logical, dimension(3) :: steps_
    integer :: iercuda, indx, indy, indz
    integer, dimension(mpi_status_size) :: istatus

    steps_ = .true. ; if(present(steps)) steps_ = steps

    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv,&
    & w_gpu => self%w_gpu, wbuf4s_gpu => self%wbuf4s_gpu, wbuf3r_gpu => self%wbuf3r_gpu,&
    & wall_tag_gpu => self%wall_tag_gpu, mp_cart => self%field%mp_cart, nrank_wake => self%field%nrank_wa&
    &ke, iermpi => self%mpi_err, nrank => self%myrank)

      if(steps_(1)) then
        call bcswap_wake_step_1_kernel(nx, ny, nz, nv, ng, w_gpu, wbuf4s_gpu)
      endif
      if(steps_(2)) then
        indy = nv*nx*ng*nz
        if (nrank_wake==nrank) then
          self%ierr = omp_target_memcpy_f(wbuf3r_gpu , wbuf4s_gpu,0,0,self%mydev,self%mydev)
        else
          call mpi_sendrecv(wbuf4s_gpu,indy,mpi_prec,nrank_wake,4,wbuf3r_gpu,indy,mpi_prec,nrank_wake,4,mp_cart,istatus,iermpi)
        endif
      endif
      if(steps_(3)) then
        call bcswap_wake_step_3_kernel(nx, ny, nz, nv, ng, w_gpu, wbuf3r_gpu, wall_tag_gpu)
      endif
    endassociate
  endsubroutine bcswap_wake

  subroutine bcswap_wake_var(self, w_swap, steps)
    class(base_omp_object), intent(inout) :: self
    real(rkind), dimension(1-self%ng:self%nx+self%ng, 1-self%ng:self%ny+self%ng, 1-self%ng:self%nz+self%ng,1:1) :: w_swap
    logical, dimension(3), optional :: steps
    logical, dimension(3) :: steps_
    integer :: iercuda, indx, indy, indz
    integer, dimension(mpi_status_size) :: istatus

    steps_ = .true. ; if(present(steps)) steps_ = steps

    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv,&
    & w_gpu => self%w_gpu, wbuf4s_gpu => self%wbuf4s_gpu, wbuf3r_gpu => self%wbuf3r_gpu,&
    & wall_tag_gpu => self%wall_tag_gpu, mp_cart => self%field%mp_cart, nrank_wake => self%field%nrank_wa&
    &ke, iermpi => self%mpi_err, nrank => self%myrank)

      if(steps_(1)) then
        call bcswap_wake_step_1_kernel(nx, ny, nz, 1, ng, w_swap, wbuf4s_gpu)
      endif
      if(steps_(2)) then
        indy = 1*nx*ng*nz
        if (nrank_wake==nrank) then
          self%ierr = omp_target_memcpy_f(wbuf3r_gpu , wbuf4s_gpu,0,0,self%mydev,self%mydev)
        else
          call mpi_sendrecv(wbuf4s_gpu,indy,mpi_prec,nrank_wake,4,wbuf3r_gpu,indy,mpi_prec,nrank_wake,4,mp_cart,istatus,iermpi)
        endif
      endif
      if(steps_(3)) then
        call bcswap_wake_step_3_kernel(nx, ny, nz, 1, ng, w_swap, wbuf3r_gpu, wall_tag_gpu)
      endif
    endassociate
  endsubroutine bcswap_wake_var

  subroutine bcswap_edges_corners(self, steps)
    class(base_omp_object), intent(inout) :: self
    logical, dimension(3), optional :: steps
    logical, dimension(3) :: steps_
    integer :: i,j,k,m, iercuda, indc
    integer, dimension(mpi_status_size) :: istatus

    steps_ = .true. ; if(present(steps)) steps_ = steps

    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv,&
    & w_gpu => self%w_gpu, wbuf1s_c_gpu => self%wbuf1s_c_gpu, wbuf2s_c_gpu => self%wbuf2s_c_gpu,&
    & wbuf3s_c_gpu => self%wbuf3s_c_gpu, wbuf4s_c_gpu => self%wbuf4s_c_gpu, wbuf1r_c_gpu => self%wbuf1r_c&
    &_gpu, wbuf2r_c_gpu => self%wbuf2r_c_gpu, wbuf3r_c_gpu => self%wbuf3r_c_gpu, wbuf4r_c_gpu => self%wbu&
    &f4r_c_gpu, ileftbottom => self%field%ileftbottom, irightbottom => self%field%irightbottom,&
    & ilefttop => self%field%ilefttop, irighttop => self%field%irighttop, pbc => self%field%grid%is_xyz_p&
    &eriodic, mp_cart => self%field%mp_cart, iermpi => self%mpi_err, nrank => self%field%myrank,&
    & lxly => self%field%lxly, lxry => self%field%lxry, rxly => self%field%rxly, rxry => self%field%rxry,&
    &lxlz => self%field%lxlz, lxrz => self%field%lxrz, rxlz => self%field%rxlz, rxrz => self%field%rxrz,&
    &lylz => self%field%lylz, lyrz => self%field%lyrz, rylz => self%field%rylz, ryrz => self%field%ryrz,&
    &lxlylz => self%field%lxlylz, lxlyrz => self%field%lxlyrz, lxrylz => self%field%lxrylz,&
    & rxlylz => self%field%rxlylz,lxryrz => self%field%lxryrz, rxlyrz => self%field%rxlyrz,&
    & rxrylz => self%field%rxrylz, rxryrz => self%field%rxryrz,wbuf_lxly_s_gpu => self%wbuf_lxly_s_gpu ,&
    & wbuf_lxry_s_gpu => self%wbuf_lxry_s_gpu, wbuf_rxly_s_gpu => self%wbuf_rxly_s_gpu ,&
    & wbuf_rxry_s_gpu => self%wbuf_rxry_s_gpu, wbuf_lxlz_s_gpu => self%wbuf_lxlz_s_gpu ,&
    & wbuf_lxrz_s_gpu => self%wbuf_lxrz_s_gpu, wbuf_rxlz_s_gpu => self%wbuf_rxlz_s_gpu ,&
    & wbuf_rxrz_s_gpu => self%wbuf_rxrz_s_gpu, wbuf_lylz_s_gpu => self%wbuf_lylz_s_gpu ,&
    & wbuf_lyrz_s_gpu => self%wbuf_lyrz_s_gpu, wbuf_rylz_s_gpu => self%wbuf_rylz_s_gpu ,&
    & wbuf_ryrz_s_gpu => self%wbuf_ryrz_s_gpu, wbuf_lxlylz_s_gpu => self%wbuf_lxlylz_s_gpu ,&
    & wbuf_lxlyrz_s_gpu => self%wbuf_lxlyrz_s_gpu, wbuf_lxrylz_s_gpu => self%wbuf_lxrylz_s_gpu ,&
    & wbuf_lxryrz_s_gpu => self%wbuf_lxryrz_s_gpu, wbuf_rxlylz_s_gpu => self%wbuf_rxlylz_s_gpu ,&
    & wbuf_rxlyrz_s_gpu => self%wbuf_rxlyrz_s_gpu, wbuf_rxrylz_s_gpu => self%wbuf_rxrylz_s_gpu ,&
    & wbuf_rxryrz_s_gpu => self%wbuf_rxryrz_s_gpu, wbuf_lxly_r_gpu => self%wbuf_lxly_r_gpu ,&
    & wbuf_lxry_r_gpu => self%wbuf_lxry_r_gpu, wbuf_rxly_r_gpu => self%wbuf_rxly_r_gpu ,&
    & wbuf_rxry_r_gpu => self%wbuf_rxry_r_gpu, wbuf_lxlz_r_gpu => self%wbuf_lxlz_r_gpu ,&
    & wbuf_lxrz_r_gpu => self%wbuf_lxrz_r_gpu, wbuf_rxlz_r_gpu => self%wbuf_rxlz_r_gpu ,&
    & wbuf_rxrz_r_gpu => self%wbuf_rxrz_r_gpu, wbuf_lylz_r_gpu => self%wbuf_lylz_r_gpu ,&
    & wbuf_lyrz_r_gpu => self%wbuf_lyrz_r_gpu, wbuf_rylz_r_gpu => self%wbuf_rylz_r_gpu ,&
    & wbuf_ryrz_r_gpu => self%wbuf_ryrz_r_gpu, wbuf_lxlylz_r_gpu => self%wbuf_lxlylz_r_gpu ,&
    & wbuf_lxlyrz_r_gpu => self%wbuf_lxlyrz_r_gpu, wbuf_lxrylz_r_gpu => self%wbuf_lxrylz_r_gpu ,&
    & wbuf_lxryrz_r_gpu => self%wbuf_lxryrz_r_gpu, wbuf_rxlylz_r_gpu => self%wbuf_rxlylz_r_gpu ,&
    & wbuf_rxlyrz_r_gpu => self%wbuf_rxlyrz_r_gpu, wbuf_rxrylz_r_gpu => self%wbuf_rxrylz_r_gpu ,&
    & wbuf_rxryrz_r_gpu => self%wbuf_rxryrz_r_gpu)

      if(steps_(1)) then
        call bcswap_corner_step_1_kernel(nx, ny, nz, nv, ng, w_gpu, wbuf1s_c_gpu, wbuf2s_c_gpu, wbuf3s_c_gpu, wbuf4s_c_gpu)
        if(pbc(2)) then
          call bcswap_corner_step_1b_kernel(nx, ny, nz, nv, ng, w_gpu, wbuf_lxly_s_gpu,&
          & wbuf_lxry_s_gpu, wbuf_rxly_s_gpu, wbuf_rxry_s_gpu, wbuf_lylz_s_gpu, wbuf_lyrz_s_gpu,&
          & wbuf_rylz_s_gpu, wbuf_ryrz_s_gpu, wbuf_lxlylz_s_gpu, wbuf_lxlyrz_s_gpu, wbuf_lxrylz_s_gpu,&
          & wbuf_lxryrz_s_gpu, wbuf_rxlylz_s_gpu, wbuf_rxlyrz_s_gpu, wbuf_rxrylz_s_gpu, wbuf_rxryrz_s_gpu)
        endif
      endif
      if(steps_(2)) then
        indc = nv*ny*ng*ng
        if (ileftbottom == nrank) then
          self%ierr = omp_target_memcpy_f(wbuf2r_c_gpu , wbuf1s_c_gpu,0,0,self%mydev,self%mydev)
          self%ierr = omp_target_memcpy_f(wbuf1r_c_gpu , wbuf2s_c_gpu,0,0,self%mydev,self%mydev)
        else
          call mpi_sendrecv(wbuf1s_c_gpu,indc,mpi_prec,ileftbottom ,1, wbuf2r_c_gpu,indc,mpi_prec,&
          &irighttop ,1,mp_cart,istatus,iermpi)
          call mpi_sendrecv(wbuf2s_c_gpu,indc,mpi_prec,irighttop ,2, wbuf1r_c_gpu,indc,mpi_prec,&
          &ileftbottom ,2,mp_cart,istatus,iermpi)
        endif
        if (ilefttop == nrank) then
          self%ierr = omp_target_memcpy_f(wbuf4r_c_gpu , wbuf3s_c_gpu,0,0,self%mydev,self%mydev)
          self%ierr = omp_target_memcpy_f(wbuf3r_c_gpu , wbuf4s_c_gpu,0,0,self%mydev,self%mydev)
        else
          call mpi_sendrecv(wbuf3s_c_gpu,indc,mpi_prec,ilefttop ,3, wbuf4r_c_gpu,indc,mpi_prec,&
          &irightbottom,3,mp_cart,istatus,iermpi)
          call mpi_sendrecv(wbuf4s_c_gpu,indc,mpi_prec,irightbottom,4, wbuf3r_c_gpu,indc,mpi_prec,&
          &ilefttop ,4,mp_cart,istatus,iermpi)
        endif
        if(pbc(2)) then
          indc = nv*ng*ng*nz
          if (lxly == nrank) then
            self%ierr = omp_target_memcpy_f(wbuf_rxry_r_gpu , wbuf_lxly_s_gpu,0,0,self%mydev,self%mydev)
            self%ierr = omp_target_memcpy_f(wbuf_lxly_r_gpu , wbuf_rxry_s_gpu,0,0,self%mydev,self%mydev)
          else
            call mpi_sendrecv(wbuf_lxly_s_gpu,indc,mpi_prec,lxly ,5, wbuf_rxry_r_gpu,indc,mpi_prec,rxry ,5,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rxry_s_gpu,indc,mpi_prec,rxry ,6, wbuf_lxly_r_gpu,indc,mpi_prec,lxly ,6,mp_cart,istatus,iermpi)
          endif
          if (lxry == nrank) then
            self%ierr = omp_target_memcpy_f(wbuf_rxly_r_gpu , wbuf_lxry_s_gpu,0,0,self%mydev,self%mydev)
            self%ierr = omp_target_memcpy_f(wbuf_lxry_r_gpu , wbuf_rxly_s_gpu,0,0,self%mydev,self%mydev)
          else
            call mpi_sendrecv(wbuf_lxry_s_gpu,indc,mpi_prec,lxry ,7, wbuf_rxly_r_gpu,indc,mpi_prec,rxly ,7,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rxly_s_gpu,indc,mpi_prec,rxly ,8, wbuf_lxry_r_gpu,indc,mpi_prec,lxry ,8,mp_cart,istatus,iermpi)
          endif
          indc = nv*nx*ng*ng
          if (lylz == nrank) then
            self%ierr = omp_target_memcpy_f(wbuf_ryrz_r_gpu , wbuf_lylz_s_gpu,0,0,self%mydev,self%mydev)
            self%ierr = omp_target_memcpy_f(wbuf_lylz_r_gpu , wbuf_ryrz_s_gpu,0,0,self%mydev,self%mydev)
          else
            call mpi_sendrecv(wbuf_lylz_s_gpu,indc,mpi_prec,lylz ,9, wbuf_ryrz_r_gpu,indc,mpi_prec,ryrz ,9,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_ryrz_s_gpu,indc,mpi_prec,ryrz ,10, wbuf_lylz_r_gpu,indc,mpi_prec,&
            &lylz ,10,mp_cart,istatus,iermpi)
          endif
          if (lyrz == nrank) then
            self%ierr = omp_target_memcpy_f(wbuf_rylz_r_gpu , wbuf_lyrz_s_gpu,0,0,self%mydev,self%mydev)
            self%ierr = omp_target_memcpy_f(wbuf_lyrz_r_gpu , wbuf_rylz_s_gpu,0,0,self%mydev,self%mydev)
          else
            call mpi_sendrecv(wbuf_lyrz_s_gpu,indc,mpi_prec,lyrz ,11, wbuf_rylz_r_gpu,indc,mpi_prec,&
            &rylz ,11,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rylz_s_gpu,indc,mpi_prec,rylz ,12, wbuf_lyrz_r_gpu,indc,mpi_prec,&
            &lyrz ,12,mp_cart,istatus,iermpi)
          endif
          indc = nv*ng*ng*ng
          if (lxlylz == nrank) then
            self%ierr = omp_target_memcpy_f(wbuf_rxryrz_r_gpu , wbuf_lxlylz_s_gpu,0,0,self%mydev,self%mydev)
            self%ierr = omp_target_memcpy_f(wbuf_rxrylz_r_gpu , wbuf_lxlyrz_s_gpu,0,0,self%mydev,self%mydev)
            self%ierr = omp_target_memcpy_f(wbuf_rxlyrz_r_gpu , wbuf_lxrylz_s_gpu,0,0,self%mydev,self%mydev)
            self%ierr = omp_target_memcpy_f(wbuf_rxlylz_r_gpu , wbuf_lxryrz_s_gpu,0,0,self%mydev,self%mydev)
            self%ierr = omp_target_memcpy_f(wbuf_lxryrz_r_gpu , wbuf_rxlylz_s_gpu,0,0,self%mydev,self%mydev)
            self%ierr = omp_target_memcpy_f(wbuf_lxrylz_r_gpu , wbuf_rxlyrz_s_gpu,0,0,self%mydev,self%mydev)
            self%ierr = omp_target_memcpy_f(wbuf_lxlyrz_r_gpu , wbuf_rxrylz_s_gpu,0,0,self%mydev,self%mydev)
            self%ierr = omp_target_memcpy_f(wbuf_lxlylz_r_gpu , wbuf_rxryrz_s_gpu,0,0,self%mydev,self%mydev)
          else
            call mpi_sendrecv(wbuf_lxlylz_s_gpu,indc,mpi_prec,lxlylz ,13, wbuf_rxryrz_r_gpu,indc,&
            &mpi_prec,rxryrz ,13,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rxryrz_s_gpu,indc,mpi_prec,rxryrz ,14, wbuf_lxlylz_r_gpu,indc,&
            &mpi_prec,lxlylz ,14,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_lxlyrz_s_gpu,indc,mpi_prec,lxlyrz ,15, wbuf_rxrylz_r_gpu,indc,&
            &mpi_prec,rxrylz ,15,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rxrylz_s_gpu,indc,mpi_prec,rxrylz ,16, wbuf_lxlyrz_r_gpu,indc,&
            &mpi_prec,lxlyrz ,16,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rxlylz_s_gpu,indc,mpi_prec,rxlylz ,17, wbuf_lxryrz_r_gpu,indc,&
            &mpi_prec,lxryrz ,17,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_lxryrz_s_gpu,indc,mpi_prec,lxryrz ,18, wbuf_rxlylz_r_gpu,indc,&
            &mpi_prec,rxlylz ,18,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rxlyrz_s_gpu,indc,mpi_prec,rxlyrz ,19, wbuf_lxrylz_r_gpu,indc,&
            &mpi_prec,lxrylz ,19,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_lxrylz_s_gpu,indc,mpi_prec,lxrylz ,20, wbuf_rxlyrz_r_gpu,indc,&
            &mpi_prec,rxlyrz ,20,mp_cart,istatus,iermpi)
          endif
        endif
      endif
      if(steps_(3)) then
        call bcswap_corner_step_3_kernel(nx, ny, nz, nv, ng, w_gpu, wbuf1r_c_gpu, wbuf2r_c_gpu,&
        & wbuf3r_c_gpu, wbuf4r_c_gpu, ileftbottom, ilefttop, irightbottom, irighttop)
        if(pbc(2)) then
          call bcswap_corner_step_3b_kernel(nx, ny, nz, nv, ng, w_gpu, wbuf_lxly_r_gpu,&
          & wbuf_lxry_r_gpu, wbuf_rxly_r_gpu, wbuf_rxry_r_gpu, wbuf_lylz_r_gpu, wbuf_lyrz_r_gpu,&
          & wbuf_rylz_r_gpu, wbuf_ryrz_r_gpu, wbuf_lxlylz_r_gpu, wbuf_lxlyrz_r_gpu, wbuf_lxrylz_r_gpu,&
          & wbuf_lxryrz_r_gpu, wbuf_rxlylz_r_gpu, wbuf_rxlyrz_r_gpu, wbuf_rxrylz_r_gpu, wbuf_rxryrz_r_gpu,&
          & lxly, lxry, rxly, rxry, lylz, lyrz, rylz, ryrz, lxlylz, lxlyrz, lxrylz, lxryrz, rxlylz, rxlyrz,&
          & rxrylz, rxryrz)
        endif
      endif
    endassociate
  endsubroutine bcswap_edges_corners

  subroutine bcswap_edges_corners_var(self, w_swap, steps)
    class(base_omp_object), intent(inout) :: self
    real(rkind), dimension(1-self%ng:self%nx+self%ng, 1-self%ng:self%ny+self%ng, 1-self%ng:self%nz+self%ng,1:1) :: w_swap
    logical, dimension(3), optional :: steps
    logical, dimension(3) :: steps_
    integer :: i,j,k,m, iercuda, indc
    integer, dimension(mpi_status_size) :: istatus

    steps_ = .true. ; if(present(steps)) steps_ = steps

    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv,&
    & wbuf1s_c_gpu => self%wbuf1s_c_gpu, wbuf2s_c_gpu => self%wbuf2s_c_gpu, wbuf3s_c_gpu => self%wbuf3s_c&
    &_gpu, wbuf4s_c_gpu => self%wbuf4s_c_gpu, wbuf1r_c_gpu => self%wbuf1r_c_gpu, wbuf2r_c_gpu => self%wbu&
    &f2r_c_gpu, wbuf3r_c_gpu => self%wbuf3r_c_gpu, wbuf4r_c_gpu => self%wbuf4r_c_gpu, ileftbottom => self&
    &%field%ileftbottom, irightbottom => self%field%irightbottom, ilefttop => self%field%ilefttop,&
    & irighttop => self%field%irighttop, pbc => self%field%grid%is_xyz_periodic, mp_cart => self%field%mp&
    &_cart, iermpi => self%mpi_err, nrank => self%field%myrank, lxly => self%field%lxly,&
    & lxry => self%field%lxry, rxly => self%field%rxly, rxry => self%field%rxry,lxlz => self%field%lxlz,&
    & lxrz => self%field%lxrz, rxlz => self%field%rxlz, rxrz => self%field%rxrz,lylz => self%field%lylz,&
    & lyrz => self%field%lyrz, rylz => self%field%rylz, ryrz => self%field%ryrz,lxlylz => self%field%lxly&
    &lz, lxlyrz => self%field%lxlyrz, lxrylz => self%field%lxrylz, rxlylz => self%field%rxlylz,&
    &lxryrz => self%field%lxryrz, rxlyrz => self%field%rxlyrz, rxrylz => self%field%rxrylz,&
    & rxryrz => self%field%rxryrz,wbuf_lxly_s_gpu => self%wbuf_lxly_s_gpu , wbuf_lxry_s_gpu => self%wbuf_&
    &lxry_s_gpu, wbuf_rxly_s_gpu => self%wbuf_rxly_s_gpu , wbuf_rxry_s_gpu => self%wbuf_rxry_s_gpu,&
    & wbuf_lxlz_s_gpu => self%wbuf_lxlz_s_gpu , wbuf_lxrz_s_gpu => self%wbuf_lxrz_s_gpu,&
    & wbuf_rxlz_s_gpu => self%wbuf_rxlz_s_gpu , wbuf_rxrz_s_gpu => self%wbuf_rxrz_s_gpu,&
    & wbuf_lylz_s_gpu => self%wbuf_lylz_s_gpu , wbuf_lyrz_s_gpu => self%wbuf_lyrz_s_gpu,&
    & wbuf_rylz_s_gpu => self%wbuf_rylz_s_gpu , wbuf_ryrz_s_gpu => self%wbuf_ryrz_s_gpu,&
    & wbuf_lxlylz_s_gpu => self%wbuf_lxlylz_s_gpu , wbuf_lxlyrz_s_gpu => self%wbuf_lxlyrz_s_gpu,&
    & wbuf_lxrylz_s_gpu => self%wbuf_lxrylz_s_gpu , wbuf_lxryrz_s_gpu => self%wbuf_lxryrz_s_gpu,&
    & wbuf_rxlylz_s_gpu => self%wbuf_rxlylz_s_gpu , wbuf_rxlyrz_s_gpu => self%wbuf_rxlyrz_s_gpu,&
    & wbuf_rxrylz_s_gpu => self%wbuf_rxrylz_s_gpu , wbuf_rxryrz_s_gpu => self%wbuf_rxryrz_s_gpu,&
    & wbuf_lxly_r_gpu => self%wbuf_lxly_r_gpu , wbuf_lxry_r_gpu => self%wbuf_lxry_r_gpu,&
    & wbuf_rxly_r_gpu => self%wbuf_rxly_r_gpu , wbuf_rxry_r_gpu => self%wbuf_rxry_r_gpu,&
    & wbuf_lxlz_r_gpu => self%wbuf_lxlz_r_gpu , wbuf_lxrz_r_gpu => self%wbuf_lxrz_r_gpu,&
    & wbuf_rxlz_r_gpu => self%wbuf_rxlz_r_gpu , wbuf_rxrz_r_gpu => self%wbuf_rxrz_r_gpu,&
    & wbuf_lylz_r_gpu => self%wbuf_lylz_r_gpu , wbuf_lyrz_r_gpu => self%wbuf_lyrz_r_gpu,&
    & wbuf_rylz_r_gpu => self%wbuf_rylz_r_gpu , wbuf_ryrz_r_gpu => self%wbuf_ryrz_r_gpu,&
    & wbuf_lxlylz_r_gpu => self%wbuf_lxlylz_r_gpu , wbuf_lxlyrz_r_gpu => self%wbuf_lxlyrz_r_gpu,&
    & wbuf_lxrylz_r_gpu => self%wbuf_lxrylz_r_gpu , wbuf_lxryrz_r_gpu => self%wbuf_lxryrz_r_gpu,&
    & wbuf_rxlylz_r_gpu => self%wbuf_rxlylz_r_gpu , wbuf_rxlyrz_r_gpu => self%wbuf_rxlyrz_r_gpu,&
    & wbuf_rxrylz_r_gpu => self%wbuf_rxrylz_r_gpu , wbuf_rxryrz_r_gpu => self%wbuf_rxryrz_r_gpu)

      if(steps_(1)) then
        call bcswap_corner_step_1_kernel(nx, ny, nz, 1, ng, w_swap, wbuf1s_c_gpu, wbuf2s_c_gpu, wbuf3s_c_gpu, wbuf4s_c_gpu)
        if(pbc(2)) then
          call bcswap_corner_step_1b_kernel(nx, ny, nz, 1, ng, w_swap, wbuf_lxly_s_gpu,&
          & wbuf_lxry_s_gpu, wbuf_rxly_s_gpu, wbuf_rxry_s_gpu, wbuf_lylz_s_gpu, wbuf_lyrz_s_gpu,&
          & wbuf_rylz_s_gpu, wbuf_ryrz_s_gpu, wbuf_lxlylz_s_gpu, wbuf_lxlyrz_s_gpu, wbuf_lxrylz_s_gpu,&
          & wbuf_lxryrz_s_gpu, wbuf_rxlylz_s_gpu, wbuf_rxlyrz_s_gpu, wbuf_rxrylz_s_gpu, wbuf_rxryrz_s_gpu)
        endif
      endif
      if(steps_(2)) then
        indc = 1*ny*ng*ng
        if (ileftbottom == nrank) then
          self%ierr = omp_target_memcpy_f(wbuf2r_c_gpu , wbuf1s_c_gpu,0,0,self%mydev,self%mydev)
          self%ierr = omp_target_memcpy_f(wbuf1r_c_gpu , wbuf2s_c_gpu,0,0,self%mydev,self%mydev)
        else
          call mpi_sendrecv(wbuf1s_c_gpu,indc,mpi_prec,ileftbottom ,1, wbuf2r_c_gpu,indc,mpi_prec,&
          &irighttop ,1,mp_cart,istatus,iermpi)
          call mpi_sendrecv(wbuf2s_c_gpu,indc,mpi_prec,irighttop ,2, wbuf1r_c_gpu,indc,mpi_prec,&
          &ileftbottom ,2,mp_cart,istatus,iermpi)
        endif
        if (ilefttop == nrank) then
          self%ierr = omp_target_memcpy_f(wbuf4r_c_gpu , wbuf3s_c_gpu,0,0,self%mydev,self%mydev)
          self%ierr = omp_target_memcpy_f(wbuf3r_c_gpu , wbuf4s_c_gpu,0,0,self%mydev,self%mydev)
        else
          call mpi_sendrecv(wbuf3s_c_gpu,indc,mpi_prec,ilefttop ,3, wbuf4r_c_gpu,indc,mpi_prec,&
          &irightbottom,3,mp_cart,istatus,iermpi)
          call mpi_sendrecv(wbuf4s_c_gpu,indc,mpi_prec,irightbottom,4, wbuf3r_c_gpu,indc,mpi_prec,&
          &ilefttop ,4,mp_cart,istatus,iermpi)
        endif
        if(pbc(2)) then
          indc = 1*ng*ng*nz
          if (lxly == nrank) then
            self%ierr = omp_target_memcpy_f(wbuf_rxry_r_gpu , wbuf_lxly_s_gpu,0,0,self%mydev,self%mydev)
            self%ierr = omp_target_memcpy_f(wbuf_lxly_r_gpu , wbuf_rxry_s_gpu,0,0,self%mydev,self%mydev)
          else
            call mpi_sendrecv(wbuf_lxly_s_gpu,indc,mpi_prec,lxly ,5, wbuf_rxry_r_gpu,indc,mpi_prec,rxry ,5,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rxry_s_gpu,indc,mpi_prec,rxry ,6, wbuf_lxly_r_gpu,indc,mpi_prec,lxly ,6,mp_cart,istatus,iermpi)
          endif
          if (lxry == nrank) then
            self%ierr = omp_target_memcpy_f(wbuf_rxly_r_gpu , wbuf_lxry_s_gpu,0,0,self%mydev,self%mydev)
            self%ierr = omp_target_memcpy_f(wbuf_lxry_r_gpu , wbuf_rxly_s_gpu,0,0,self%mydev,self%mydev)
          else
            call mpi_sendrecv(wbuf_lxry_s_gpu,indc,mpi_prec,lxry ,7, wbuf_rxly_r_gpu,indc,mpi_prec,rxly ,7,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rxly_s_gpu,indc,mpi_prec,rxly ,8, wbuf_lxry_r_gpu,indc,mpi_prec,lxry ,8,mp_cart,istatus,iermpi)
          endif
          indc = 1*nx*ng*ng
          if (lylz == nrank) then
            self%ierr = omp_target_memcpy_f(wbuf_ryrz_r_gpu , wbuf_lylz_s_gpu,0,0,self%mydev,self%mydev)
            self%ierr = omp_target_memcpy_f(wbuf_lylz_r_gpu , wbuf_ryrz_s_gpu,0,0,self%mydev,self%mydev)
          else
            call mpi_sendrecv(wbuf_lylz_s_gpu,indc,mpi_prec,lylz ,9, wbuf_ryrz_r_gpu,indc,mpi_prec,ryrz ,9,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_ryrz_s_gpu,indc,mpi_prec,ryrz ,10, wbuf_lylz_r_gpu,indc,mpi_prec,&
            &lylz ,10,mp_cart,istatus,iermpi)
          endif
          if (lyrz == nrank) then
            self%ierr = omp_target_memcpy_f(wbuf_rylz_r_gpu , wbuf_lyrz_s_gpu,0,0,self%mydev,self%mydev)
            self%ierr = omp_target_memcpy_f(wbuf_lyrz_r_gpu , wbuf_rylz_s_gpu,0,0,self%mydev,self%mydev)
          else
            call mpi_sendrecv(wbuf_lyrz_s_gpu,indc,mpi_prec,lyrz ,11, wbuf_rylz_r_gpu,indc,mpi_prec,&
            &rylz ,11,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rylz_s_gpu,indc,mpi_prec,rylz ,12, wbuf_lyrz_r_gpu,indc,mpi_prec,&
            &lyrz ,12,mp_cart,istatus,iermpi)
          endif
          indc = 1*ng*ng*ng
          if (lxlylz == nrank) then
            self%ierr = omp_target_memcpy_f(wbuf_rxryrz_r_gpu , wbuf_lxlylz_s_gpu,0,0,self%mydev,self%mydev)
            self%ierr = omp_target_memcpy_f(wbuf_rxrylz_r_gpu , wbuf_lxlyrz_s_gpu,0,0,self%mydev,self%mydev)
            self%ierr = omp_target_memcpy_f(wbuf_rxlyrz_r_gpu , wbuf_lxrylz_s_gpu,0,0,self%mydev,self%mydev)
            self%ierr = omp_target_memcpy_f(wbuf_rxlylz_r_gpu , wbuf_lxryrz_s_gpu,0,0,self%mydev,self%mydev)
            self%ierr = omp_target_memcpy_f(wbuf_lxryrz_r_gpu , wbuf_rxlylz_s_gpu,0,0,self%mydev,self%mydev)
            self%ierr = omp_target_memcpy_f(wbuf_lxrylz_r_gpu , wbuf_rxlyrz_s_gpu,0,0,self%mydev,self%mydev)
            self%ierr = omp_target_memcpy_f(wbuf_lxlyrz_r_gpu , wbuf_rxrylz_s_gpu,0,0,self%mydev,self%mydev)
            self%ierr = omp_target_memcpy_f(wbuf_lxlylz_r_gpu , wbuf_rxryrz_s_gpu,0,0,self%mydev,self%mydev)
          else
            call mpi_sendrecv(wbuf_lxlylz_s_gpu,indc,mpi_prec,lxlylz ,13, wbuf_rxryrz_r_gpu,indc,&
            &mpi_prec,rxryrz ,13,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rxryrz_s_gpu,indc,mpi_prec,rxryrz ,14, wbuf_lxlylz_r_gpu,indc,&
            &mpi_prec,lxlylz ,14,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_lxlyrz_s_gpu,indc,mpi_prec,lxlyrz ,15, wbuf_rxrylz_r_gpu,indc,&
            &mpi_prec,rxrylz ,15,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rxrylz_s_gpu,indc,mpi_prec,rxrylz ,16, wbuf_lxlyrz_r_gpu,indc,&
            &mpi_prec,lxlyrz ,16,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rxlylz_s_gpu,indc,mpi_prec,rxlylz ,17, wbuf_lxryrz_r_gpu,indc,&
            &mpi_prec,lxryrz ,17,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_lxryrz_s_gpu,indc,mpi_prec,lxryrz ,18, wbuf_rxlylz_r_gpu,indc,&
            &mpi_prec,rxlylz ,18,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rxlyrz_s_gpu,indc,mpi_prec,rxlyrz ,19, wbuf_lxrylz_r_gpu,indc,&
            &mpi_prec,lxrylz ,19,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_lxrylz_s_gpu,indc,mpi_prec,lxrylz ,20, wbuf_rxlyrz_r_gpu,indc,&
            &mpi_prec,rxlyrz ,20,mp_cart,istatus,iermpi)
          endif
        endif
      endif
      if(steps_(3)) then
        call bcswap_corner_step_3_kernel(nx, ny, nz, 1, ng, w_swap, wbuf1r_c_gpu, wbuf2r_c_gpu,&
        & wbuf3r_c_gpu, wbuf4r_c_gpu, ileftbottom, ilefttop, irightbottom, irighttop)
        if(pbc(2)) then
          call bcswap_corner_step_3b_kernel(nx, ny, nz, 1, ng, w_swap, wbuf_lxly_r_gpu,&
          & wbuf_lxry_r_gpu, wbuf_rxly_r_gpu, wbuf_rxry_r_gpu, wbuf_lylz_r_gpu, wbuf_lyrz_r_gpu,&
          & wbuf_rylz_r_gpu, wbuf_ryrz_r_gpu, wbuf_lxlylz_r_gpu, wbuf_lxlyrz_r_gpu, wbuf_lxrylz_r_gpu,&
          & wbuf_lxryrz_r_gpu, wbuf_rxlylz_r_gpu, wbuf_rxlyrz_r_gpu, wbuf_rxrylz_r_gpu, wbuf_rxryrz_r_gpu,&
          & lxly, lxry, rxly, rxry, lylz, lyrz, rylz, ryrz, lxlylz, lxlyrz, lxrylz, lxryrz, rxlylz, rxlyrz,&
          & rxrylz, rxryrz)
        endif
      endif
    endassociate
  endsubroutine bcswap_edges_corners_var


  subroutine bcswap_var(self, w_swap, steps)
    class(base_omp_object), intent(inout) :: self
    real(rkind), dimension(1-self%ng:self%nx+self%ng, 1-self%ng:self%ny+self%ng, 1-self%ng:self%nz+self%ng,1:1) :: w_swap
    logical, dimension(3), optional :: steps
    logical, dimension(3) :: steps_
    integer :: iercuda, indx, indy, indz
    integer, dimension(mpi_status_size) :: istatus
    integer, dimension(3) :: is_periodic

    steps_ = .true. ; if(present(steps)) steps_ = steps

    is_periodic(:) = 0
    if(self%field%grid%is_xyz_periodic(1)) is_periodic(1) = 1
    if(self%field%grid%is_xyz_periodic(2)) is_periodic(2) = 1
    if(self%field%grid%is_xyz_periodic(3)) is_periodic(3) = 1

    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, wbuf1s_gpu => self%wbuf1s_&
    &gpu, wbuf2s_gpu => self%wbuf2s_gpu, wbuf3s_gpu => self%wbuf3s_gpu, wbuf4s_gpu => self%wbuf4s_gpu,&
    & wbuf5s_gpu => self%wbuf5s_gpu, wbuf6s_gpu => self%wbuf6s_gpu, wbuf1r_gpu => self%wbuf1r_gpu,&
    & wbuf2r_gpu => self%wbuf2r_gpu, wbuf3r_gpu => self%wbuf3r_gpu, wbuf4r_gpu => self%wbuf4r_gpu,&
    & wbuf5r_gpu => self%wbuf5r_gpu, wbuf6r_gpu => self%wbuf6r_gpu, ileftx => self%field%ileftx,&
    &irightx => self%field%irightx, nrank_x => self%field%nrank_x, ilefty => self%field%ilefty,&
    &irighty => self%field%irighty, nrank_y => self%field%nrank_y, ileftz => self%field%ileftz,&
    &irightz => self%field%irightz, nrank_z => self%field%nrank_z, dcsidxnc2_gpu => self%dcsidxnc2_gpu,&
    & dcsidync2_gpu => self%dcsidync2_gpu, detadxnc2_gpu => self%detadxnc2_gpu, detadync2_gpu => self%det&
    &adync2_gpu, dxdcsinc2_gpu => self%dxdcsinc2_gpu, dydcsinc2_gpu => self%dydcsinc2_gpu,&
    & dxdetanc2_gpu => self%dxdetanc2_gpu, dydetanc2_gpu => self%dydetanc2_gpu, mp_cartx => self%field%mp&
    &_cartx, mp_carty => self%field%mp_carty, mp_cartz => self%field%mp_cartz, ndim => self%field%grid%nd&
    &im)

      if(steps_(1)) then
        if(self%field%grid%grid_dim == 1) then
          call bcswap_step_1_kernel(nx, ny, nz, 1, ng, ndim, w_swap, wbuf1s_gpu, wbuf2s_gpu,&
          & wbuf3s_gpu, wbuf4s_gpu, wbuf5s_gpu, wbuf6s_gpu)
        elseif(self%field%grid%grid_dim == 2) then
          call bcswap_step_1_kernel(nx, ny, nz, 1, ng, ndim, w_swap, wbuf1s_gpu, wbuf2s_gpu,&
          & wbuf3s_gpu, wbuf4s_gpu, wbuf5s_gpu, wbuf6s_gpu)
        endif
      endif
      if(steps_(2)) then
        indx = ng*ny*nz
        indy = nx*ng*nz
        indz = nx*ny*ng
        if(ileftx == nrank_x) then
          self%ierr = omp_target_memcpy_f(wbuf2r_gpu , wbuf1s_gpu,0,0,self%mydev,self%mydev)
          self%ierr = omp_target_memcpy_f(wbuf1r_gpu , wbuf2s_gpu,0,0,self%mydev,self%mydev)
        else
          call mpi_sendrecv(wbuf1s_gpu,indx,mpi_prec,ileftx ,1,wbuf2r_gpu,indx,mpi_prec,irightx,1,mp_cartx,istatus,self%mpi_err)
          call mpi_sendrecv(wbuf2s_gpu,indx,mpi_prec,irightx,2,wbuf1r_gpu,indx,mpi_prec,ileftx ,2,mp_cartx,istatus,self%mpi_err)
        endif
        if(ilefty == nrank_y) then
          self%ierr = omp_target_memcpy_f(wbuf4r_gpu , wbuf3s_gpu,0,0,self%mydev,self%mydev)
          self%ierr = omp_target_memcpy_f(wbuf3r_gpu , wbuf4s_gpu,0,0,self%mydev,self%mydev)
        else
          call mpi_sendrecv(wbuf3s_gpu,indy,mpi_prec,ilefty ,3,wbuf4r_gpu,indy,mpi_prec,irighty,3,mp_carty,istatus,self%mpi_err)
          call mpi_sendrecv(wbuf4s_gpu,indy,mpi_prec,irighty,4,wbuf3r_gpu,indy,mpi_prec,ilefty ,4,mp_carty,istatus,self%mpi_err)
        endif
        if(ndim == 3) then
          if(ileftz == nrank_z) then
            self%ierr = omp_target_memcpy_f(wbuf6r_gpu , wbuf5s_gpu,0,0,self%mydev,self%mydev)
            self%ierr = omp_target_memcpy_f(wbuf5r_gpu , wbuf6s_gpu,0,0,self%mydev,self%mydev)
          else
            call mpi_sendrecv(wbuf5s_gpu,indz,mpi_prec,ileftz ,5,wbuf6r_gpu,indz,mpi_prec,irightz,5,mp_cartz,istatus,self%mpi_err)
            call mpi_sendrecv(wbuf6s_gpu,indz,mpi_prec,irightz,6,wbuf5r_gpu,indz,mpi_prec,ileftz ,6,mp_cartz,istatus,self%mpi_err)
          endif
        endif
      endif
      if(steps_(3)) then
        if(self%field%grid%grid_dim == 1) then
          call bcswap_step_3_kernel(nx, ny, nz, 1, ng, ndim, w_swap, wbuf1r_gpu, wbuf2r_gpu,&
          & wbuf3r_gpu, wbuf4r_gpu, wbuf5r_gpu, wbuf6r_gpu, ileftx, ilefty, ileftz, irightx, irighty, irightz)
        elseif(self%field%grid%grid_dim == 2) then
          call bcswap_step_3_kernel(nx, ny, nz, 1, ng, ndim, w_swap, wbuf1r_gpu, wbuf2r_gpu,&
          & wbuf3r_gpu, wbuf4r_gpu, wbuf5r_gpu, wbuf6r_gpu, ileftx, ilefty, ileftz, irightx, irighty, irightz)
        endif
      endif
    endassociate
  endsubroutine bcswap_var

  subroutine copy_cpu_gpu(self)
    class(base_omp_object), intent(inout) :: self
    integer :: i,j,k,iv

    self%ierr = omp_target_memcpy_f(self%w_gpu , self%field%w,0,0,self%mydev,self%myhost)

  endsubroutine copy_cpu_gpu


  subroutine copy_gpu_cpu(self)
    class(base_omp_object), intent(inout) :: self
    integer :: i,j,k,iv

    self%ierr = omp_target_memcpy_f(self%field%w , self%w_gpu,0,0,self%myhost,self%mydev)

  endsubroutine copy_gpu_cpu

  subroutine initialize(self, field, gpu_bind)
    !< Initialize base backend.
    class(base_omp_object), intent(inout) :: self !< The base backend.
    class(field_object), target :: field
    integer, intent(in) :: gpu_bind

    self%field => field
    self%nx = self%field%nx
    self%ny = self%field%ny
    self%nz = self%field%nz
    self%ng = self%field%grid%ng
    self%nv = self%field%nv

    call get_mpi_basic_info(self%nprocs, self%myrank, self%masterproc, self%mpi_err)

    call self%field%check_cpu_mem(description="--Base-initialization--")

    call MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, self%local_comm, self%mpi_err)
    if(gpu_bind == 0) then
      call MPI_COMM_RANK(self%local_comm, self%mydev, self%mpi_err)
    else
      self%mydev = 0
    endif
    call omp_set_default_device(self%mydev)
    self%myhost = omp_get_initial_device()

    call self%alloc()

  endsubroutine initialize

  subroutine check_gpu_mem(self, description)
    class(base_omp_object), intent(inout) :: self !< The base backend.
    character(*) :: description
    integer :: ierr
    !integer(cuda_count_kind) :: mem_free, mem_total
    !character(128) :: proc_name
    !integer :: resultlen
    !!call mpi_barrier(mpi_comm_world, ierr)
    !call mpi_get_processor_name(proc_name, resultlen, ierr)
    !ierr = cudaMemGetInfo(mem_free, mem_total)
    !!call mpi_barrier(mpi_comm_world, ierr)
    !write(error_unit, "(A,2x,A,2x,A,2x,I0,2x,I0,2x,I0)") 'GPU rank,mems: ', ! description,
!    proc_name(1:resultlen),self%myrank, mem_free, mem_total
  endsubroutine check_gpu_mem

endmodule streams_base_omp_object


