module streams_base_ompc_object

  use streams_field_object, only : field_object
  use streams_parameters
  use mpi


  implicit none
  private
  public :: base_ompc_object

  type :: base_ompc_object
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
    real(rkind), pointer, contiguous, dimension(:,:,:,:) :: w_ompc

    real(rkind), dimension(:,:,:,:), allocatable :: wbuf1s_ompc, wbuf2s_ompc, wbuf3s_ompc, wbuf4s_ompc, wbuf5s_ompc, wbuf6s_ompc
    real(rkind), dimension(:,:,:,:), allocatable :: wbuf1r_ompc, wbuf2r_ompc, wbuf3r_ompc, wbuf4r_ompc, wbuf5r_ompc, wbuf6r_ompc
    real(rkind), dimension(:,:,:,:), allocatable :: wbuf1xybcs_ompc, wbuf2xybcs_ompc, wbuf1xybcr_ompc, wbuf2xybcr_ompc
    real(rkind), dimension(:,:,:,:), allocatable :: wbuf1s_c_ompc, wbuf2s_c_ompc, wbuf3s_c_ompc, wbuf4s_c_ompc
    real(rkind), dimension(:,:,:,:), allocatable :: wbuf1r_c_ompc, wbuf2r_c_ompc, wbuf3r_c_ompc, wbuf4r_c_ompc
    real(rkind), allocatable, dimension(:) :: dcsidx_ompc, dcsidx2_ompc, dcsidxs_ompc
    real(rkind), allocatable, dimension(:) :: detady_ompc, detady2_ompc, detadys_ompc
    real(rkind), allocatable, dimension(:) :: dzitdz_ompc, dzitdz2_ompc, dzitdzs_ompc
    real(rkind), allocatable, dimension(:) :: x_ompc, y_ompc, z_ompc, yn_ompc
    real(rkind), allocatable, dimension(:,:) :: xc2_ompc,yc2_ompc
    real(rkind), allocatable, dimension(:,:) :: dcsidxc2_ompc,dcsidyc2_ompc
    real(rkind), allocatable, dimension(:,:) :: detadxc2_ompc,detadyc2_ompc
    real(rkind), allocatable, dimension(:,:) :: dcsidxnc2_ompc,dcsidync2_ompc
    real(rkind), allocatable, dimension(:,:) :: detadxnc2_ompc,detadync2_ompc
    real(rkind), allocatable, dimension(:,:) :: dxdcsic2_ompc,dydcsic2_ompc
    real(rkind), allocatable, dimension(:,:) :: dxdetac2_ompc,dydetac2_ompc
    real(rkind), allocatable, dimension(:,:) :: dxdcsinc2_ompc,dydcsinc2_ompc
    real(rkind), allocatable, dimension(:,:) :: dxdetanc2_ompc,dydetanc2_ompc
    real(rkind), allocatable, dimension(:,:) :: jac_ompc,mcsijac1_ompc,metajac1_ompc
    real(rkind), allocatable, dimension(:,:) :: mcsi_ompc,meta_ompc,csimod_ompc,etamod_ompc
    real(rkind), allocatable, dimension(:,:) :: g1_ompc,g2_ompc,g12_ompc
    real(rkind), allocatable, dimension(:,:) :: wbuftus_ompc, wbuftes_ompc , wbuftur_ompc, wbufter_ompc
    real(rkind), allocatable, dimension(:,:) :: theta_ij_ompc

    integer, allocatable, dimension(:) :: wall_tag_ompc

    real(rkind), allocatable,dimension(:,:,:,:) :: wbuf_lxly_s_ompc , wbuf_lxry_s_ompc , wbuf_rxly_s_ompc , wbuf_rxry_s_ompc
    real(rkind), allocatable,dimension(:,:,:,:) :: wbuf_lxlz_s_ompc , wbuf_lxrz_s_ompc , wbuf_rxlz_s_ompc , wbuf_rxrz_s_ompc
    real(rkind), allocatable,dimension(:,:,:,:) :: wbuf_lylz_s_ompc , wbuf_lyrz_s_ompc , wbuf_rylz_s_ompc , wbuf_ryrz_s_ompc
    real(rkind), allocatable,dimension(:,:,:,:) :: wbuf_lxlylz_s_ompc,wbuf_lxlyrz_s_ompc,wbuf_lxrylz_s_ompc,wbuf_lxryrz_s_ompc
    real(rkind), allocatable,dimension(:,:,:,:) :: wbuf_rxlylz_s_ompc,wbuf_rxlyrz_s_ompc,wbuf_rxrylz_s_ompc,wbuf_rxryrz_s_ompc
    real(rkind), allocatable,dimension(:,:,:,:) :: wbuf_lxly_r_ompc , wbuf_lxry_r_ompc , wbuf_rxly_r_ompc , wbuf_rxry_r_ompc
    real(rkind), allocatable,dimension(:,:,:,:) :: wbuf_lxlz_r_ompc , wbuf_lxrz_r_ompc , wbuf_rxlz_r_ompc , wbuf_rxrz_r_ompc
    real(rkind), allocatable,dimension(:,:,:,:) :: wbuf_lylz_r_ompc , wbuf_lyrz_r_ompc , wbuf_rylz_r_ompc , wbuf_ryrz_r_ompc
    real(rkind), allocatable,dimension(:,:,:,:) :: wbuf_lxlylz_r_ompc,wbuf_lxlyrz_r_ompc,wbuf_lxrylz_r_ompc,wbuf_lxryrz_r_ompc
    real(rkind), allocatable,dimension(:,:,:,:) :: wbuf_rxlylz_r_ompc,wbuf_rxlyrz_r_ompc,wbuf_rxrylz_r_ompc,wbuf_rxryrz_r_ompc

  contains
    procedure, pass(self) :: alloc
    procedure, pass(self) :: copy_from_field
    procedure, pass(self) :: copy_to_field
    procedure, pass(self) :: initialize
    procedure, pass(self) :: bcswap
    procedure, pass(self) :: bcswap_var
    procedure, pass(self) :: bcswap_corner
    procedure, pass(self) :: bcswap_corner_var
    procedure, pass(self) :: bcswap_wake
    procedure, pass(self) :: bcswap_wake_var
    procedure, pass(self) :: bcte
    procedure, pass(self) :: bcte_var
    procedure, pass(self) :: bcswap_edges_corners
    procedure, pass(self) :: bcswap_edges_corners_var
  endtype base_ompc_object

contains
  subroutine alloc(self)
    class(base_ompc_object), intent(inout) :: self

    associate(nx => self%field%nx, ny => self%field%ny, nz => self%field%nz, ng => self%field%grid%ng, nv => self%field%nv)

      self%w_ompc => self%field%w 

      allocate(self%wbuf1s_ompc(ng,ny,nz,nv))
      allocate(self%wbuf2s_ompc(ng,ny,nz,nv))
      allocate(self%wbuf3s_ompc(nx,ng,nz,nv))
      allocate(self%wbuf4s_ompc(nx,ng,nz,nv))
      allocate(self%wbuf5s_ompc(nx,ny,ng,nv))
      allocate(self%wbuf6s_ompc(nx,ny,ng,nv))
      allocate(self%wbuf1r_ompc(ng,ny,nz,nv))
      allocate(self%wbuf2r_ompc(ng,ny,nz,nv))
      allocate(self%wbuf3r_ompc(nx,ng,nz,nv))
      allocate(self%wbuf4r_ompc(nx,ng,nz,nv))
      allocate(self%wbuf5r_ompc(nx,ny,ng,nv))
      allocate(self%wbuf6r_ompc(nx,ny,ng,nv))

      allocate(self%wbuf1xybcs_ompc(ng,1-ng:ny+ng,nz,nv))
      allocate(self%wbuf2xybcs_ompc(ng,1-ng:ny+ng,nz,nv))
      allocate(self%wbuf1xybcr_ompc(ng,1-ng:ny+ng,nz,nv))
      allocate(self%wbuf2xybcr_ompc(ng,1-ng:ny+ng,nz,nv))
      allocate(self%x_ompc(1-ng:nx+ng), self%y_ompc(1-ng:ny+ng), self%z_ompc(1-ng:nz+ng))
      allocate(self%yn_ompc(1:ny+1))

      allocate(self%dcsidx_ompc(nx), self%dcsidx2_ompc(nx), self%dcsidxs_ompc(nx))
      allocate(self%detady_ompc(ny), self%detady2_ompc(ny), self%detadys_ompc(ny))
      allocate(self%dzitdz_ompc(nz), self%dzitdz2_ompc(nz), self%dzitdzs_ompc(nz))

      allocate(self%wbuf1s_c_ompc(nv,ng,ny,ng))
      allocate(self%wbuf2s_c_ompc(nv,ng,ny,ng))
      allocate(self%wbuf3s_c_ompc(nv,ng,ny,ng))
      allocate(self%wbuf4s_c_ompc(nv,ng,ny,ng))
      allocate(self%wbuf1r_c_ompc(nv,ng,ny,ng))
      allocate(self%wbuf2r_c_ompc(nv,ng,ny,ng))
      allocate(self%wbuf3r_c_ompc(nv,ng,ny,ng))
      allocate(self%wbuf4r_c_ompc(nv,ng,ny,ng))

      allocate(self%wbuf_lxly_s_ompc(nv,ng,ng,nz), self%wbuf_lxry_s_ompc(nv,ng,ng,nz))
      allocate(self%wbuf_rxly_s_ompc(nv,ng,ng,nz), self%wbuf_rxry_s_ompc(nv,ng,ng,nz))
      allocate(self%wbuf_lxlz_s_ompc(nv,ng,ny,ng), self%wbuf_lxrz_s_ompc(nv,ng,ny,ng))
      allocate(self%wbuf_rxlz_s_ompc(nv,ng,ny,ng), self%wbuf_rxrz_s_ompc(nv,ng,ny,ng))
      allocate(self%wbuf_lylz_s_ompc(nv,nx,ng,ng), self%wbuf_lyrz_s_ompc(nv,nx,ng,ng))
      allocate(self%wbuf_rylz_s_ompc(nv,nx,ng,ng), self%wbuf_ryrz_s_ompc(nv,nx,ng,ng))
      allocate(self%wbuf_lxlylz_s_ompc(nv,ng,ng,ng),self%wbuf_lxlyrz_s_ompc(nv,ng,ng,ng))
      allocate(self%wbuf_lxrylz_s_ompc(nv,ng,ng,ng),self%wbuf_lxryrz_s_ompc(nv,ng,ng,ng))
      allocate(self%wbuf_rxlylz_s_ompc(nv,ng,ng,ng),self%wbuf_rxlyrz_s_ompc(nv,ng,ng,ng))
      allocate(self%wbuf_rxrylz_s_ompc(nv,ng,ng,ng),self%wbuf_rxryrz_s_ompc(nv,ng,ng,ng))

      allocate(self%wbuf_lxly_r_ompc(nv,ng,ng,nz), self%wbuf_lxry_r_ompc(nv,ng,ng,nz))
      allocate(self%wbuf_rxly_r_ompc(nv,ng,ng,nz), self%wbuf_rxry_r_ompc(nv,ng,ng,nz))
      allocate(self%wbuf_lxlz_r_ompc(nv,ng,ny,ng), self%wbuf_lxrz_r_ompc(nv,ng,ny,ng))
      allocate(self%wbuf_rxlz_r_ompc(nv,ng,ny,ng), self%wbuf_rxrz_r_ompc(nv,ng,ny,ng))
      allocate(self%wbuf_lylz_r_ompc(nv,nx,ng,ng), self%wbuf_lyrz_r_ompc(nv,nx,ng,ng))
      allocate(self%wbuf_rylz_r_ompc(nv,nx,ng,ng), self%wbuf_ryrz_r_ompc(nv,nx,ng,ng))
      allocate(self%wbuf_lxlylz_r_ompc(nv,ng,ng,ng),self%wbuf_lxlyrz_r_ompc(nv,ng,ng,ng))
      allocate(self%wbuf_lxrylz_r_ompc(nv,ng,ng,ng),self%wbuf_lxryrz_r_ompc(nv,ng,ng,ng))
      allocate(self%wbuf_rxlylz_r_ompc(nv,ng,ng,ng),self%wbuf_rxlyrz_r_ompc(nv,ng,ng,ng))
      allocate(self%wbuf_rxrylz_r_ompc(nv,ng,ng,ng),self%wbuf_rxryrz_r_ompc(nv,ng,ng,ng))

      if (self%field%grid%grid_dim == 2) then
        allocate(self%xc2_ompc(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%yc2_ompc(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%dcsidxc2_ompc(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%dcsidyc2_ompc(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%detadxc2_ompc(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%detadyc2_ompc(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%dcsidxnc2_ompc(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%dcsidync2_ompc(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%detadxnc2_ompc(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%detadync2_ompc(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%dxdcsic2_ompc(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%dydcsic2_ompc(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%dxdetac2_ompc(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%dydetac2_ompc(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%dxdcsinc2_ompc(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%dydcsinc2_ompc(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%dxdetanc2_ompc(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%dydetanc2_ompc(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%mcsi_ompc(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%meta_ompc(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%mcsijac1_ompc(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%metajac1_ompc(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%jac_ompc(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%g1_ompc(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%g2_ompc(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%g12_ompc(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%csimod_ompc(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%etamod_ompc(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%theta_ij_ompc(1:nx,1:ny))

        allocate(self%wall_tag_ompc(1-ng:nx+ng))
        allocate(self%wbuftus_ompc(nz,2))
        allocate(self%wbuftes_ompc(nz,2))
        allocate(self%wbuftur_ompc(nz,2))
        allocate(self%wbufter_ompc(nz,2))
      endif

    endassociate

    self%x_ompc = self%field%x
    self%y_ompc = self%field%y
    self%z_ompc = self%field%z
    self%yn_ompc = self%field%yn

    self%dcsidx_ompc = self%field%dcsidx
    self%dcsidxs_ompc = self%field%dcsidxs
    self%dcsidx2_ompc = self%field%dcsidx2
    self%detady_ompc = self%field%detady
    self%detadys_ompc = self%field%detadys
    self%detady2_ompc = self%field%detady2
    self%dzitdz_ompc = self%field%dzitdz
    self%dzitdzs_ompc = self%field%dzitdzs
    self%dzitdz2_ompc = self%field%dzitdz2

    if (self%field%grid%grid_dim == 2) then
      self%xc2_ompc = self%field%xc2
      self%yc2_ompc = self%field%yc2
      self%dcsidxc2_ompc = self%field%dcsidxc2
      self%dcsidyc2_ompc = self%field%dcsidyc2
      self%detadxc2_ompc = self%field%detadxc2
      self%detadyc2_ompc = self%field%detadyc2
      self%dcsidxnc2_ompc = self%field%dcsidxnc2
      self%dcsidync2_ompc = self%field%dcsidync2
      self%detadxnc2_ompc = self%field%detadxnc2
      self%detadync2_ompc = self%field%detadync2
      self%dxdcsic2_ompc = self%field%dxdcsic2
      self%dydcsic2_ompc = self%field%dydcsic2
      self%dxdetac2_ompc = self%field%dxdetac2
      self%dydetac2_ompc = self%field%dydetac2
      self%dxdcsinc2_ompc = self%field%dxdcsinc2
      self%dydcsinc2_ompc = self%field%dydcsinc2
      self%dxdetanc2_ompc = self%field%dxdetanc2
      self%dydetanc2_ompc = self%field%dydetanc2
      self%mcsi_ompc = self%field%mcsi
      self%meta_ompc = self%field%meta
      self%mcsijac1_ompc = self%field%mcsijac1
      self%metajac1_ompc = self%field%metajac1
      self%jac_ompc = self%field%jac
      self%g1_ompc = self%field%g1
      self%g2_ompc = self%field%g2
      self%g12_ompc = self%field%g12
      self%csimod_ompc = self%field%csimod
      self%etamod_ompc = self%field%etamod
      self%wall_tag_ompc = self%field%wall_tag
      self%theta_ij_ompc = self%field%theta_ij
    endif

  endsubroutine alloc

  subroutine bcte_step_1_subroutine(nx,ny,nz,nv,ng,nrank_x,ite_rank_x,itu_rank_x,ite_l,itu_l,w_ompc,wbuftus_ompc,wbuftes_ompc)
    integer, intent(in) :: nx, ny, nz, nv, ng
    integer, intent(in) :: nrank_x, ite_rank_x, itu_rank_x, ite_l, itu_l
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in) :: w_ompc
    real(rkind), dimension(nz, 2), intent(out) :: wbuftus_ompc, wbuftes_ompc
    integer :: i,j,k,m,iercuda
    if(nrank_x == ite_rank_x) then
      !$omp parallel do default(firstprivate) shared(w_ompc, &
      !$omp& wbuftus_ompc, &
      !$omp& wbuftes_ompc)
      do k = 1,nz
        wbuftes_ompc(k,1) = w_ompc(ite_l,1,k,1)
        if(nv > 1) wbuftes_ompc(k,2) = w_ompc(ite_l,1,k,5)
      enddo
    endif
    if(nrank_x == itu_rank_x) then
      !$omp parallel do default(firstprivate) shared(w_ompc, &
      !$omp& wbuftus_ompc, &
      !$omp& wbuftes_ompc)
      do k = 1,nz
        wbuftus_ompc(k,1) = w_ompc(itu_l,1,k,1)
        if(nv > 1) wbuftus_ompc(k,2) = w_ompc(itu_l,1,k,5)
      enddo
    endif

  endsubroutine bcte_step_1_subroutine



  subroutine bcte(self, steps)
    class(base_ompc_object), intent(inout) :: self
    logical, dimension(3), optional :: steps
    logical, dimension(3) :: steps_
    integer, dimension(4) :: requests

    steps_ = .true. ; if(present(steps)) steps_ = steps
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv,&
    & w_ompc => self%w_ompc, wbuftes_ompc => self%wbuftes_ompc, wbuftus_ompc => self%wbuftus_ompc,&
    & wbufter_ompc => self%wbufter_ompc, wbuftur_ompc => self%wbuftur_ompc, ite_rank_x => self%field%ite_&
    &rank_x, itu_rank_x => self%field%itu_rank_x, ite_l => self%field%ite_l, itu_l => self%field%itu_l,&
    & nrank_x => self%field%nrank_x, mp_cartx => self%field%mp_cartx, iermpi => self%mpi_err)

      if(steps_(1)) then
        call bcte_step_1_subroutine(nx, ny, nz, nv, ng, nrank_x, ite_rank_x, itu_rank_x, ite_l,&
        & itu_l, w_ompc, wbuftus_ompc, wbuftes_ompc)
      endif

      if(steps_(2)) then
        if (nrank_x == ite_rank_x .and. nrank_x == itu_rank_x) then
          wbufter_ompc = wbuftes_ompc
          wbuftur_ompc = wbuftus_ompc
        else
          requests(:) = mpi_request_null
          if (nrank_x == ite_rank_x) then
            call mpi_isend(wbuftes_ompc,2*nz,mpi_prec,itu_rank_x,1,mp_cartx,requests(1),iermpi)
            call mpi_irecv(wbuftur_ompc,2*nz,mpi_prec,itu_rank_x,1,mp_cartx,requests(2),iermpi)
          elseif (nrank_x == itu_rank_x) then
            call mpi_isend(wbuftus_ompc,2*nz,mpi_prec,ite_rank_x,1,mp_cartx,requests(3),iermpi)
            call mpi_irecv(wbufter_ompc,2*nz,mpi_prec,ite_rank_x,1,mp_cartx,requests(4),iermpi)
          endif
          call mpi_waitall(4,requests,mpi_statuses_ignore,iermpi)
        endif
      endif

      if(steps_(3)) then
        call bcte_step_3_subroutine(nx, ny, nz, nv, ng, nrank_x, ite_rank_x, itu_rank_x, ite_l,&
        & itu_l, w_ompc, wbuftur_ompc, wbufter_ompc)
      endif
    endassociate
  endsubroutine bcte

  subroutine bcte_var(self, w_swap, steps)
    class(base_ompc_object), intent(inout) :: self
    real(rkind), dimension(1-self%ng:self%nx+self%ng, 1-self%ng:self%ny+self%ng, 1-self%ng:self%nz+self%ng,1:1) :: w_swap
    logical, dimension(3), optional :: steps
    logical, dimension(3) :: steps_
    integer, dimension(4) :: requests

    steps_ = .true. ; if(present(steps)) steps_ = steps
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv,&
    & w_ompc => self%w_ompc, wbuftes_ompc => self%wbuftes_ompc, wbuftus_ompc => self%wbuftus_ompc,&
    & wbufter_ompc => self%wbufter_ompc, wbuftur_ompc => self%wbuftur_ompc, ite_rank_x => self%field%ite_&
    &rank_x, itu_rank_x => self%field%itu_rank_x, ite_l => self%field%ite_l, itu_l => self%field%itu_l,&
    & nrank_x => self%field%nrank_x, mp_cartx => self%field%mp_cartx, iermpi => self%mpi_err)

      if(steps_(1)) then
        call bcte_step_1_subroutine(nx, ny, nz, 1, ng, nrank_x, ite_rank_x, itu_rank_x, ite_l,&
        & itu_l, w_swap, wbuftus_ompc, wbuftes_ompc)
      endif

      if(steps_(2)) then
        if (nrank_x == ite_rank_x .and. nrank_x == itu_rank_x) then
          wbufter_ompc = wbuftes_ompc
          wbuftur_ompc = wbuftus_ompc
        else
          requests(:) = mpi_request_null
          if (nrank_x == ite_rank_x) then
            call mpi_isend(wbuftes_ompc,1*nz,mpi_prec,itu_rank_x,3,mp_cartx,requests(1),iermpi)
            call mpi_irecv(wbuftur_ompc,1*nz,mpi_prec,itu_rank_x,3,mp_cartx,requests(2),iermpi)
          elseif (nrank_x == itu_rank_x) then
            call mpi_isend(wbuftus_ompc,1*nz,mpi_prec,ite_rank_x,3,mp_cartx,requests(3),iermpi)
            call mpi_irecv(wbufter_ompc,1*nz,mpi_prec,ite_rank_x,3,mp_cartx,requests(4),iermpi)
          endif
          call mpi_waitall(4,requests,mpi_statuses_ignore,iermpi)
        endif
      endif

      if(steps_(3)) then
        call bcte_step_3_subroutine(nx, ny, nz, 1, ng, nrank_x, ite_rank_x, itu_rank_x, ite_l,&
        & itu_l, w_swap, wbuftur_ompc, wbufter_ompc)
      endif
    endassociate
  endsubroutine bcte_var

  subroutine bcte_step_3_subroutine(nx,ny,nz,nv,ng,nrank_x,ite_rank_x,itu_rank_x,ite_l,itu_l,w_ompc,wbuftur_ompc,wbufter_ompc)
    integer, intent(in) :: nx, ny, nz, nv, ng
    integer, intent(in) :: nrank_x, ite_rank_x, itu_rank_x, ite_l, itu_l
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(inout) :: w_ompc
    real(rkind), dimension(nz, 2), intent(in) :: wbuftur_ompc, wbufter_ompc
    integer :: k, iercuda
    if(nrank_x == ite_rank_x) then
      !$omp parallel do default(firstprivate) shared(w_ompc, &
      !$omp& wbuftur_ompc, &
      !$omp& wbufter_ompc)
      do k = 1,nz
        w_ompc(ite_l,1,k,1) = 0.5_rkind*(w_ompc(ite_l,1,k,1)+wbuftur_ompc(k,1))
        if(nv > 1) w_ompc(ite_l,1,k,5) = 0.5_rkind*(w_ompc(ite_l,1,k,5)+wbuftur_ompc(k,2))
      enddo
    endif

    if (nrank_x==itu_rank_x) then
      !$omp parallel do default(firstprivate) shared(w_ompc, &
      !$omp& wbuftur_ompc, &
      !$omp& wbufter_ompc)
      do k = 1,nz
        w_ompc(itu_l,1,k,1) = 0.5_rkind*(w_ompc(itu_l,1,k,1)+wbufter_ompc(k,1))
        if(nv > 1) w_ompc(itu_l,1,k,5) = 0.5_rkind*(w_ompc(itu_l,1,k,5)+wbufter_ompc(k,2))
      enddo
    endif

  endsubroutine bcte_step_3_subroutine



  subroutine bcswap_wake_step_1_subroutine(nx,ny,nz,nv,ng,w_ompc,wbuf4s_ompc)
    integer :: nx, ny, nz, ng, nv
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: w_ompc
    real(rkind), dimension(1:,1:,1:,1:) :: wbuf4s_ompc
    integer :: i,j,k,m,iercuda

    !$omp parallel do default(firstprivate) shared(w_ompc, &
    !$omp& wbuf4s_ompc)
    do k = 1,nz
      do j = 1,ng
        do i = 1,nx
          do m=1,nv
            wbuf4s_ompc(i,j,k,m) = w_ompc(nx-i+1,1+j,k,m)
          enddo
        enddo
      enddo
    enddo


  endsubroutine bcswap_wake_step_1_subroutine



  subroutine bcswap_wake_step_3_subroutine(nx,ny,nz,nv,ng,w_ompc,wbuf3r_ompc,wall_tag_ompc)
    integer :: nx, ny, nz, ng, nv
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: w_ompc
    real(rkind), dimension(1:,1:,1:,1:) :: wbuf3r_ompc
    integer, dimension(1-ng:nx+ng) :: wall_tag_ompc
    integer :: i,j,k,m, iercuda

    !$omp parallel do default(firstprivate) shared(w_ompc, &
    !$omp& wbuf3r_ompc, &
    !$omp& wall_tag_ompc)
    do k = 1,nz
      do j = 1,ng
        do i = 1,nx
          do m=1,nv
            if (wall_tag_ompc(i) > 0) then
              w_ompc(i,1-j,k,m) = wbuf3r_ompc(i,j,k,m)
            endif
          enddo
        enddo
      enddo
    enddo


  endsubroutine bcswap_wake_step_3_subroutine



  subroutine bcswap_step_1_subroutine(nx,ny,nz,nv,ng,ndim,w_ompc,wbuf1s_ompc,wbuf2s_ompc,&
  &wbuf3s_ompc,wbuf4s_ompc,wbuf5s_ompc,wbuf6s_ompc)
    integer :: nx, ny, nz, ng, nv, ndim
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: w_ompc
    real(rkind), dimension(1:,1:,1:,1:) :: wbuf1s_ompc, wbuf2s_ompc, wbuf3s_ompc, wbuf4s_ompc, wbuf5s_ompc, wbuf6s_ompc
    integer :: i,j,k,m,iercuda

    !$omp parallel do default(firstprivate) shared(w_ompc, &
    !$omp& wbuf1s_ompc, &
    !$omp& wbuf2s_ompc, &
    !$omp& wbuf3s_ompc, &
    !$omp& wbuf4s_ompc, &
    !$omp& wbuf5s_ompc, &
    !$omp& wbuf6s_ompc)
    do k = 1,nz
      do j = 1,ny
        do i = 1,ng
          do m=1,nv
            wbuf1s_ompc(i,j,k,m) = w_ompc(i,j,k,m)
            wbuf2s_ompc(i,j,k,m) = w_ompc(nx-ng+i,j,k,m)
          enddo
        enddo
      enddo
    enddo

    !$omp parallel do default(firstprivate) shared(w_ompc, &
    !$omp& wbuf1s_ompc, &
    !$omp& wbuf2s_ompc, &
    !$omp& wbuf3s_ompc, &
    !$omp& wbuf4s_ompc, &
    !$omp& wbuf5s_ompc, &
    !$omp& wbuf6s_ompc)
    do k = 1,nz
      do j = 1,ng
        do i = 1,nx
          do m=1,nv
            wbuf3s_ompc(i,j,k,m) = w_ompc(i,j,k,m)
            wbuf4s_ompc(i,j,k,m) = w_ompc(i,ny-ng+j,k,m)
          enddo
        enddo
      enddo
    enddo
    if (ndim==3) then
      !$omp parallel do default(firstprivate) shared(w_ompc, &
      !$omp& wbuf1s_ompc, &
      !$omp& wbuf2s_ompc, &
      !$omp& wbuf3s_ompc, &
      !$omp& wbuf4s_ompc, &
      !$omp& wbuf5s_ompc, &
      !$omp& wbuf6s_ompc)
      do k = 1,ng
        do j = 1,ny
          do i = 1,nx
            do m=1,nv
              wbuf5s_ompc(i,j,k,m) = w_ompc(i,j,k,m)
              wbuf6s_ompc(i,j,k,m) = w_ompc(i,j,nz-ng+k,m)
            enddo
          enddo
        enddo
      enddo
    endif

  endsubroutine bcswap_step_1_subroutine



  subroutine bcswap_c2_step_1_subroutine(nx,ny,nz,nv,ng,ndim,w_ompc,wbuf1s_ompc,wbuf2s_ompc,&
  &wbuf3s_ompc,wbuf4s_ompc,wbuf5s_ompc,wbuf6s_ompc,dcsidxnc2_ompc,dcsidync2_ompc,detadxnc2_ompc,&
  &detadync2_ompc,is_periodic)
    integer :: nx, ny, nz, ng, nv, ndim
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: w_ompc
    real(rkind), dimension(1:,1:,1:,1:) :: wbuf1s_ompc, wbuf2s_ompc, wbuf3s_ompc, wbuf4s_ompc, wbuf5s_ompc, wbuf6s_ompc
    real(rkind), dimension(1-ng:,1-ng:) :: dcsidxnc2_ompc, dcsidync2_ompc, detadxnc2_ompc, detadync2_ompc
    integer, dimension(3) :: is_periodic
    integer :: i,j,k,m,iercuda
    if (is_periodic(1) == 1) then
      !$omp parallel do default(firstprivate) shared(w_ompc, &
      !$omp& wbuf1s_ompc, &
      !$omp& wbuf2s_ompc, &
      !$omp& wbuf3s_ompc, &
      !$omp& wbuf4s_ompc, &
      !$omp& wbuf5s_ompc, &
      !$omp& wbuf6s_ompc, &
      !$omp& dcsidxnc2_ompc, &
      !$omp& dcsidync2_ompc, &
      !$omp& detadxnc2_ompc, &
      !$omp& detadync2_ompc)
      do k = 1,nz
        do j = 1,ny
          do i = 1,ng
            wbuf1s_ompc(i,j,k,1) = w_ompc(i,j,k,1)
            wbuf1s_ompc(i,j,k,2) = w_ompc(i,j,k,2)*dcsidxnc2_ompc(i,j)+w_ompc(i,j,k,3)*dcsidync2_ompc(i,j)
            wbuf1s_ompc(i,j,k,3) = w_ompc(i,j,k,2)*detadxnc2_ompc(i,j)+w_ompc(i,j,k,3)*detadync2_ompc(i,j)
            wbuf1s_ompc(i,j,k,4) = w_ompc(i,j,k,4)
            wbuf1s_ompc(i,j,k,5) = w_ompc(i,j,k,5)
            wbuf2s_ompc(i,j,k,1) = w_ompc(nx-ng+i,j,k,1)
            wbuf2s_ompc(i,j,k,2) = w_ompc(nx-ng+i,j,k,2)*dcsidxnc2_ompc(nx-ng+i,j)+w_ompc(nx-ng+i,j,k,3)*dcsidync2_ompc(nx-ng+i,j)
            wbuf2s_ompc(i,j,k,3) = w_ompc(nx-ng+i,j,k,2)*detadxnc2_ompc(nx-ng+i,j)+w_ompc(nx-ng+i,j,k,3)*detadync2_ompc(nx-ng+i,j)
            wbuf2s_ompc(i,j,k,4) = w_ompc(nx-ng+i,j,k,4)
            wbuf2s_ompc(i,j,k,5) = w_ompc(nx-ng+i,j,k,5)
          enddo
        enddo
      enddo
    else
      !$omp parallel do default(firstprivate) shared(w_ompc, &
      !$omp& wbuf1s_ompc, &
      !$omp& wbuf2s_ompc, &
      !$omp& wbuf3s_ompc, &
      !$omp& wbuf4s_ompc, &
      !$omp& wbuf5s_ompc, &
      !$omp& wbuf6s_ompc, &
      !$omp& dcsidxnc2_ompc, &
      !$omp& dcsidync2_ompc, &
      !$omp& detadxnc2_ompc, &
      !$omp& detadync2_ompc)
      do k = 1,nz
        do j = 1,ny
          do i = 1,ng
            do m=1,nv
              wbuf1s_ompc(i,j,k,m) = w_ompc(i,j,k,m)
              wbuf2s_ompc(i,j,k,m) = w_ompc(nx-ng+i,j,k,m)
            enddo
          enddo
        enddo
      enddo
    endif
    !$omp parallel do default(firstprivate) shared(w_ompc, &
    !$omp& wbuf1s_ompc, &
    !$omp& wbuf2s_ompc, &
    !$omp& wbuf3s_ompc, &
    !$omp& wbuf4s_ompc, &
    !$omp& wbuf5s_ompc, &
    !$omp& wbuf6s_ompc, &
    !$omp& dcsidxnc2_ompc, &
    !$omp& dcsidync2_ompc, &
    !$omp& detadxnc2_ompc, &
    !$omp& detadync2_ompc)
    do k = 1,nz
      do j = 1,ng
        do i = 1,nx
          do m=1,nv
            wbuf3s_ompc(i,j,k,m) = w_ompc(i,j,k,m)
            wbuf4s_ompc(i,j,k,m) = w_ompc(i,ny-ng+j,k,m)
          enddo
        enddo
      enddo
    enddo
    if (ndim==3) then
      !$omp parallel do default(firstprivate) shared(w_ompc, &
      !$omp& wbuf1s_ompc, &
      !$omp& wbuf2s_ompc, &
      !$omp& wbuf3s_ompc, &
      !$omp& wbuf4s_ompc, &
      !$omp& wbuf5s_ompc, &
      !$omp& wbuf6s_ompc, &
      !$omp& dcsidxnc2_ompc, &
      !$omp& dcsidync2_ompc, &
      !$omp& detadxnc2_ompc, &
      !$omp& detadync2_ompc)
      do k = 1,ng
        do j = 1,ny
          do i = 1,nx
            do m=1,nv
              wbuf5s_ompc(i,j,k,m) = w_ompc(i,j,k,m)
              wbuf6s_ompc(i,j,k,m) = w_ompc(i,j,nz-ng+k,m)
            enddo
          enddo
        enddo
      enddo
    endif

  endsubroutine bcswap_c2_step_1_subroutine



  subroutine bcswap_c2xybc_step_1_subroutine(nx,ny,nz,nv,ng,ndim,w_ompc,wbuf1xybcs_ompc,&
  &wbuf2xybcs_ompc,wbuf3s_ompc,wbuf4s_ompc,wbuf5s_ompc,wbuf6s_ompc,dcsidxnc2_ompc,dcsidync2_ompc,&
  &detadxnc2_ompc,detadync2_ompc,is_periodic)
    integer :: nx, ny, nz, ng, nv, ndim
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: w_ompc
    real(rkind), dimension(1:,1-ng:,1:,1:) :: wbuf1xybcs_ompc, wbuf2xybcs_ompc
    real(rkind), dimension(1:,1:,1:,1:) :: wbuf3s_ompc, wbuf4s_ompc, wbuf5s_ompc, wbuf6s_ompc
    real(rkind), dimension(1-ng:,1-ng:) :: dcsidxnc2_ompc, dcsidync2_ompc, detadxnc2_ompc, detadync2_ompc
    integer, dimension(3) :: is_periodic
    integer :: i,j,k,m,iercuda
    if (is_periodic(1) == 1) then
      !$omp parallel do default(firstprivate) shared(w_ompc, &
      !$omp& wbuf1xybcs_ompc, &
      !$omp& wbuf2xybcs_ompc, &
      !$omp& wbuf3s_ompc, &
      !$omp& wbuf4s_ompc, &
      !$omp& wbuf5s_ompc, &
      !$omp& wbuf6s_ompc, &
      !$omp& dcsidxnc2_ompc, &
      !$omp& dcsidync2_ompc, &
      !$omp& detadxnc2_ompc, &
      !$omp& detadync2_ompc)
      do k = 1,nz
        do j = 1-ng,ny+ng
          do i = 1,ng
            wbuf1xybcs_ompc(i,j,k,1) = w_ompc(i,j,k,1)
            wbuf1xybcs_ompc(i,j,k,2) = w_ompc(i,j,k,2)*dcsidxnc2_ompc(i,j)+w_ompc(i,j,k,3)*dcsidync2_ompc(i,j)
            wbuf1xybcs_ompc(i,j,k,3) = w_ompc(i,j,k,2)*detadxnc2_ompc(i,j)+w_ompc(i,j,k,3)*detadync2_ompc(i,j)
            wbuf1xybcs_ompc(i,j,k,4) = w_ompc(i,j,k,4)
            wbuf1xybcs_ompc(i,j,k,5) = w_ompc(i,j,k,5)
            wbuf2xybcs_ompc(i,j,k,1) = w_ompc(nx-ng+i,j,k,1)
            wbuf2xybcs_ompc(i,j,k,2) = w_ompc(nx-ng+i,j,k,2)*dcsidxnc2_ompc(nx-ng+i,j)+w_ompc(nx-ng+&
            &i,j,k,3)*dcsidync2_ompc(nx-ng+i,j)
            wbuf2xybcs_ompc(i,j,k,3) = w_ompc(nx-ng+i,j,k,2)*detadxnc2_ompc(nx-ng+i,j)+w_ompc(nx-ng+&
            &i,j,k,3)*detadync2_ompc(nx-ng+i,j)
            wbuf2xybcs_ompc(i,j,k,4) = w_ompc(nx-ng+i,j,k,4)
            wbuf2xybcs_ompc(i,j,k,5) = w_ompc(nx-ng+i,j,k,5)
          enddo
        enddo
      enddo
    else
      !$omp parallel do default(firstprivate) shared(w_ompc, &
      !$omp& wbuf1xybcs_ompc, &
      !$omp& wbuf2xybcs_ompc, &
      !$omp& wbuf3s_ompc, &
      !$omp& wbuf4s_ompc, &
      !$omp& wbuf5s_ompc, &
      !$omp& wbuf6s_ompc, &
      !$omp& dcsidxnc2_ompc, &
      !$omp& dcsidync2_ompc, &
      !$omp& detadxnc2_ompc, &
      !$omp& detadync2_ompc)
      do k = 1,nz
        do j = 1-ng,ny+ng
          do i = 1,ng
            do m=1,nv
              wbuf1xybcs_ompc(i,j,k,m) = w_ompc(i,j,k,m)
              wbuf2xybcs_ompc(i,j,k,m) = w_ompc(nx-ng+i,j,k,m)
            enddo
          enddo
        enddo
      enddo
    endif
    !$omp parallel do default(firstprivate) shared(w_ompc, &
    !$omp& wbuf1xybcs_ompc, &
    !$omp& wbuf2xybcs_ompc, &
    !$omp& wbuf3s_ompc, &
    !$omp& wbuf4s_ompc, &
    !$omp& wbuf5s_ompc, &
    !$omp& wbuf6s_ompc, &
    !$omp& dcsidxnc2_ompc, &
    !$omp& dcsidync2_ompc, &
    !$omp& detadxnc2_ompc, &
    !$omp& detadync2_ompc)
    do k = 1,nz
      do j = 1,ng
        do i = 1,nx
          do m=1,nv
            wbuf3s_ompc(i,j,k,m) = w_ompc(i,j,k,m)
            wbuf4s_ompc(i,j,k,m) = w_ompc(i,ny-ng+j,k,m)
          enddo
        enddo
      enddo
    enddo
    if (ndim==3) then
      !$omp parallel do default(firstprivate) shared(w_ompc, &
      !$omp& wbuf1xybcs_ompc, &
      !$omp& wbuf2xybcs_ompc, &
      !$omp& wbuf3s_ompc, &
      !$omp& wbuf4s_ompc, &
      !$omp& wbuf5s_ompc, &
      !$omp& wbuf6s_ompc, &
      !$omp& dcsidxnc2_ompc, &
      !$omp& dcsidync2_ompc, &
      !$omp& detadxnc2_ompc, &
      !$omp& detadync2_ompc)
      do k = 1,ng
        do j = 1,ny
          do i = 1,nx
            do m=1,nv
              wbuf5s_ompc(i,j,k,m) = w_ompc(i,j,k,m)
              wbuf6s_ompc(i,j,k,m) = w_ompc(i,j,nz-ng+k,m)
            enddo
          enddo
        enddo
      enddo
    endif

  endsubroutine bcswap_c2xybc_step_1_subroutine



  subroutine bcswap_corner_step_1_subroutine(nx,ny,nz,nv,ng,w_ompc,wbuf1s_c_ompc,wbuf2s_c_ompc,wbuf3s_c_ompc,wbuf4s_c_ompc)
    integer :: nx, ny, nz, ng, nv
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: w_ompc
    real(rkind), dimension(1:,1:,1:,1:) :: wbuf1s_c_ompc, wbuf2s_c_ompc, wbuf3s_c_ompc, wbuf4s_c_ompc
    integer :: i,j,k,m,iercuda

    !$omp parallel do default(firstprivate) shared(w_ompc, &
    !$omp& wbuf1s_c_ompc, &
    !$omp& wbuf2s_c_ompc, &
    !$omp& wbuf3s_c_ompc, &
    !$omp& wbuf4s_c_ompc)
    do j = 1,ny
      do m = 1,nv
        do k=1,ng
          do i=1,ng
            wbuf1s_c_ompc(m,i,j,k) = w_ompc(i,j,k,m)
            wbuf2s_c_ompc(m,i,j,k) = w_ompc(nx-ng+i,j,nz-ng+k,m)
            wbuf3s_c_ompc(m,i,j,k) = w_ompc(i,j,nz-ng+k,m)
            wbuf4s_c_ompc(m,i,j,k) = w_ompc(nx-ng+i,j,k,m)
          enddo
        enddo
      enddo
    enddo


  endsubroutine bcswap_corner_step_1_subroutine



  subroutine bcswap_corner_step_1b_subroutine(nx,ny,nz,nv,ng,w_ompc,wbuf_lxly_s_ompc,&
  &wbuf_lxry_s_ompc,wbuf_rxly_s_ompc,wbuf_rxry_s_ompc,wbuf_lylz_s_ompc,wbuf_lyrz_s_ompc,&
  &wbuf_rylz_s_ompc,wbuf_ryrz_s_ompc,wbuf_lxlylz_s_ompc,wbuf_lxlyrz_s_ompc,wbuf_lxrylz_s_ompc,&
  &wbuf_lxryrz_s_ompc,wbuf_rxlylz_s_ompc,wbuf_rxlyrz_s_ompc,wbuf_rxrylz_s_ompc,wbuf_rxryrz_s_ompc)
    integer :: nx, ny, nz, ng, nv
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: w_ompc
    real(rkind), dimension(1:,1:,1:,1:) :: wbuf_lxly_s_ompc, wbuf_lxry_s_ompc, wbuf_rxly_s_ompc, wbuf_rxry_s_ompc
    real(rkind), dimension(1:,1:,1:,1:) :: wbuf_lylz_s_ompc, wbuf_lyrz_s_ompc, wbuf_rylz_s_ompc, wbuf_ryrz_s_ompc
    real(rkind), dimension(1:,1:,1:,1:) :: wbuf_lxlylz_s_ompc, wbuf_lxlyrz_s_ompc, wbuf_lxrylz_s_ompc, wbuf_lxryrz_s_ompc
    real(rkind), dimension(1:,1:,1:,1:) :: wbuf_rxlylz_s_ompc, wbuf_rxlyrz_s_ompc, wbuf_rxrylz_s_ompc, wbuf_rxryrz_s_ompc
    integer :: i,j,k,m,iercuda

    !$omp parallel do default(firstprivate) shared(w_ompc, &
    !$omp& wbuf_lxly_s_ompc, &
    !$omp& wbuf_lxry_s_ompc, &
    !$omp& wbuf_rxly_s_ompc, &
    !$omp& wbuf_rxry_s_ompc, &
    !$omp& wbuf_lylz_s_ompc, &
    !$omp& wbuf_lyrz_s_ompc, &
    !$omp& wbuf_rylz_s_ompc, &
    !$omp& wbuf_ryrz_s_ompc, &
    !$omp& wbuf_lxlylz_s_ompc, &
    !$omp& wbuf_lxlyrz_s_ompc, &
    !$omp& wbuf_lxrylz_s_ompc, &
    !$omp& wbuf_lxryrz_s_ompc, &
    !$omp& wbuf_rxlylz_s_ompc, &
    !$omp& wbuf_rxlyrz_s_ompc, &
    !$omp& wbuf_rxrylz_s_ompc, &
    !$omp& wbuf_rxryrz_s_ompc)
    do k = 1,nz
      do m = 1,nv
        do j=1,ng
          do i=1,ng
            wbuf_lxly_s_ompc(m,i,j,k) = w_ompc(i,j,k,m)
            wbuf_lxry_s_ompc(m,i,j,k) = w_ompc(i,ny-ng+j,k,m)
            wbuf_rxly_s_ompc(m,i,j,k) = w_ompc(nx-ng+i,j,k,m)
            wbuf_rxry_s_ompc(m,i,j,k) = w_ompc(nx-ng+i,ny-ng+j,k,m)
          enddo
        enddo
      enddo
    enddo

    !$omp parallel do default(firstprivate) shared(w_ompc, &
    !$omp& wbuf_lxly_s_ompc, &
    !$omp& wbuf_lxry_s_ompc, &
    !$omp& wbuf_rxly_s_ompc, &
    !$omp& wbuf_rxry_s_ompc, &
    !$omp& wbuf_lylz_s_ompc, &
    !$omp& wbuf_lyrz_s_ompc, &
    !$omp& wbuf_rylz_s_ompc, &
    !$omp& wbuf_ryrz_s_ompc, &
    !$omp& wbuf_lxlylz_s_ompc, &
    !$omp& wbuf_lxlyrz_s_ompc, &
    !$omp& wbuf_lxrylz_s_ompc, &
    !$omp& wbuf_lxryrz_s_ompc, &
    !$omp& wbuf_rxlylz_s_ompc, &
    !$omp& wbuf_rxlyrz_s_ompc, &
    !$omp& wbuf_rxrylz_s_ompc, &
    !$omp& wbuf_rxryrz_s_ompc)
    do i = 1,nx
      do m = 1,nv
        do k=1,ng
          do j=1,ng
            wbuf_lylz_s_ompc(m,i,j,k) = w_ompc(i,j,k,m)
            wbuf_lyrz_s_ompc(m,i,j,k) = w_ompc(i,j,nz-ng+k,m)
            wbuf_rylz_s_ompc(m,i,j,k) = w_ompc(i,ny-ng+j,k,m)
            wbuf_ryrz_s_ompc(m,i,j,k) = w_ompc(i,ny-ng+j,nz-ng+k,m)
          enddo
        enddo
      enddo
    enddo

    !$omp parallel do default(firstprivate) shared(w_ompc, &
    !$omp& wbuf_lxly_s_ompc, &
    !$omp& wbuf_lxry_s_ompc, &
    !$omp& wbuf_rxly_s_ompc, &
    !$omp& wbuf_rxry_s_ompc, &
    !$omp& wbuf_lylz_s_ompc, &
    !$omp& wbuf_lyrz_s_ompc, &
    !$omp& wbuf_rylz_s_ompc, &
    !$omp& wbuf_ryrz_s_ompc, &
    !$omp& wbuf_lxlylz_s_ompc, &
    !$omp& wbuf_lxlyrz_s_ompc, &
    !$omp& wbuf_lxrylz_s_ompc, &
    !$omp& wbuf_lxryrz_s_ompc, &
    !$omp& wbuf_rxlylz_s_ompc, &
    !$omp& wbuf_rxlyrz_s_ompc, &
    !$omp& wbuf_rxrylz_s_ompc, &
    !$omp& wbuf_rxryrz_s_ompc)
    do k = 1,ng
      do m = 1,nv
        do j=1,ng
          do i=1,ng
            wbuf_lxlylz_s_ompc(m,i,j,k) = w_ompc(i,j,k,m)
            wbuf_lxlyrz_s_ompc(m,i,j,k) = w_ompc(i,j,nz-ng+k,m)
            wbuf_lxrylz_s_ompc(m,i,j,k) = w_ompc(i,ny-ng+j,k,m)
            wbuf_lxryrz_s_ompc(m,i,j,k) = w_ompc(i,ny-ng+j,nz-ng+k,m)
            wbuf_rxlylz_s_ompc(m,i,j,k) = w_ompc(nx-ng+i,j,k,m)
            wbuf_rxlyrz_s_ompc(m,i,j,k) = w_ompc(nx-ng+i,j,nz-ng+k,m)
            wbuf_rxrylz_s_ompc(m,i,j,k) = w_ompc(nx-ng+i,ny-ng+j,k,m)
            wbuf_rxryrz_s_ompc(m,i,j,k) = w_ompc(nx-ng+i,ny-ng+j,nz-ng+k,m)
          enddo
        enddo
      enddo
    enddo


  endsubroutine bcswap_corner_step_1b_subroutine



  subroutine bcswap_step_3_subroutine(nx,ny,nz,nv,ng,ndim,w_ompc,wbuf1r_ompc,wbuf2r_ompc,&
  &wbuf3r_ompc,wbuf4r_ompc,wbuf5r_ompc,wbuf6r_ompc,ileftx,ilefty,ileftz,irightx,irighty,irightz)
    integer :: nx, ny, nz, ng, nv, ndim
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: w_ompc
    real(rkind), dimension(1:,1:,1:,1:) :: wbuf1r_ompc, wbuf2r_ompc, wbuf3r_ompc, wbuf4r_ompc, wbuf5r_ompc, wbuf6r_ompc
    integer :: i,j,k,m, iercuda
    integer :: ileftx, ilefty, ileftz, irightx, irighty, irightz
    if (ileftx/=mpi_proc_null) then
      !$omp parallel do default(firstprivate) shared(w_ompc, &
      !$omp& wbuf1r_ompc, &
      !$omp& wbuf2r_ompc, &
      !$omp& wbuf3r_ompc, &
      !$omp& wbuf4r_ompc, &
      !$omp& wbuf5r_ompc, &
      !$omp& wbuf6r_ompc)
      do k = 1,nz
        do j = 1,ny
          do i = 1,ng
            do m=1,nv
              w_ompc(i-ng,j,k,m) = wbuf1r_ompc(i,j,k,m)
            enddo
          enddo
        enddo
      enddo
    endif
    if (irightx/=mpi_proc_null) then
      !$omp parallel do default(firstprivate) shared(w_ompc, &
      !$omp& wbuf1r_ompc, &
      !$omp& wbuf2r_ompc, &
      !$omp& wbuf3r_ompc, &
      !$omp& wbuf4r_ompc, &
      !$omp& wbuf5r_ompc, &
      !$omp& wbuf6r_ompc)
      do k = 1,nz
        do j = 1,ny
          do i = 1,ng
            do m=1,nv
              w_ompc(nx+i,j,k,m) = wbuf2r_ompc(i,j,k,m)
            enddo
          enddo
        enddo
      enddo
    endif
    if (ilefty/=mpi_proc_null) then
      !$omp parallel do default(firstprivate) shared(w_ompc, &
      !$omp& wbuf1r_ompc, &
      !$omp& wbuf2r_ompc, &
      !$omp& wbuf3r_ompc, &
      !$omp& wbuf4r_ompc, &
      !$omp& wbuf5r_ompc, &
      !$omp& wbuf6r_ompc)
      do k = 1,nz
        do j = 1,ng
          do i = 1,nx
            do m=1,nv
              w_ompc(i,j-ng,k,m) = wbuf3r_ompc(i,j,k,m)
            enddo
          enddo
        enddo
      enddo
    endif
    if (irighty/=mpi_proc_null) then
      !$omp parallel do default(firstprivate) shared(w_ompc, &
      !$omp& wbuf1r_ompc, &
      !$omp& wbuf2r_ompc, &
      !$omp& wbuf3r_ompc, &
      !$omp& wbuf4r_ompc, &
      !$omp& wbuf5r_ompc, &
      !$omp& wbuf6r_ompc)
      do k = 1,nz
        do j = 1,ng
          do i = 1,nx
            do m=1,nv
              w_ompc(i,ny+j,k,m) = wbuf4r_ompc(i,j,k,m)
            enddo
          enddo
        enddo
      enddo
    endif
    if(ndim == 3) then
      if (ileftz/=mpi_proc_null) then
        !$omp parallel do default(firstprivate) shared(w_ompc, &
        !$omp& wbuf1r_ompc, &
        !$omp& wbuf2r_ompc, &
        !$omp& wbuf3r_ompc, &
        !$omp& wbuf4r_ompc, &
        !$omp& wbuf5r_ompc, &
        !$omp& wbuf6r_ompc)
        do k = 1,ng
          do j = 1,ny
            do i = 1,nx
              do m=1,nv
                w_ompc(i,j,k-ng,m) = wbuf5r_ompc(i,j,k,m)
              enddo
            enddo
          enddo
        enddo
      endif
      if (irightz/=mpi_proc_null) then
        !$omp parallel do default(firstprivate) shared(w_ompc, &
        !$omp& wbuf1r_ompc, &
        !$omp& wbuf2r_ompc, &
        !$omp& wbuf3r_ompc, &
        !$omp& wbuf4r_ompc, &
        !$omp& wbuf5r_ompc, &
        !$omp& wbuf6r_ompc)
        do k = 1,ng
          do j = 1,ny
            do i = 1,nx
              do m=1,nv
                w_ompc(i,j,nz+k,m) = wbuf6r_ompc(i,j,k,m)
              enddo
            enddo
          enddo
        enddo
      endif
    endif

  endsubroutine bcswap_step_3_subroutine



  subroutine bcswap_c2_step_3_subroutine(nx,ny,nz,nv,ng,ndim,w_ompc,wbuf1r_ompc,wbuf2r_ompc,&
  &wbuf3r_ompc,wbuf4r_ompc,wbuf5r_ompc,wbuf6r_ompc,ileftx,ilefty,ileftz,irightx,irighty,irightz,&
  &dxdcsinc2_ompc,dydcsinc2_ompc,dxdetanc2_ompc,dydetanc2_ompc,is_periodic)
    integer :: nx, ny, nz, ng, nv, ndim
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: w_ompc
    real(rkind), dimension(1:,1:,1:,1:) :: wbuf1r_ompc, wbuf2r_ompc, wbuf3r_ompc, wbuf4r_ompc, wbuf5r_ompc, wbuf6r_ompc
    real(rkind), dimension(1-ng:,1-ng:) :: dxdcsinc2_ompc, dydcsinc2_ompc, dxdetanc2_ompc, dydetanc2_ompc
    integer :: i,j,k,m, iercuda
    integer :: ileftx, ilefty, ileftz, irightx, irighty, irightz
    integer, dimension(3) :: is_periodic
    if (is_periodic(1) == 1) then
      if (ileftx/=mpi_proc_null) then
        !$omp parallel do default(firstprivate) shared(w_ompc, &
        !$omp& wbuf1r_ompc, &
        !$omp& wbuf2r_ompc, &
        !$omp& wbuf3r_ompc, &
        !$omp& wbuf4r_ompc, &
        !$omp& wbuf5r_ompc, &
        !$omp& wbuf6r_ompc, &
        !$omp& dxdcsinc2_ompc, &
        !$omp& dydcsinc2_ompc, &
        !$omp& dxdetanc2_ompc, &
        !$omp& dydetanc2_ompc)
        do k = 1,nz
          do j = 1,ny
            do i = 1,ng
              w_ompc(i-ng,j,k,1) = wbuf1r_ompc(i,j,k,1)
              w_ompc(i-ng,j,k,2) = wbuf1r_ompc(i,j,k,2)*dxdcsinc2_ompc(i-ng,j)+wbuf1r_ompc(i,j,k,3)*dxdetanc2_ompc(i-ng,j)
              w_ompc(i-ng,j,k,3) = wbuf1r_ompc(i,j,k,2)*dydcsinc2_ompc(i-ng,j)+wbuf1r_ompc(i,j,k,3)*dydetanc2_ompc(i-ng,j)
              w_ompc(i-ng,j,k,4) = wbuf1r_ompc(i,j,k,4)
              w_ompc(i-ng,j,k,5) = wbuf1r_ompc(i,j,k,5)
            enddo
          enddo
        enddo
      endif
      if (irightx/=mpi_proc_null) then
        !$omp parallel do default(firstprivate) shared(w_ompc, &
        !$omp& wbuf1r_ompc, &
        !$omp& wbuf2r_ompc, &
        !$omp& wbuf3r_ompc, &
        !$omp& wbuf4r_ompc, &
        !$omp& wbuf5r_ompc, &
        !$omp& wbuf6r_ompc, &
        !$omp& dxdcsinc2_ompc, &
        !$omp& dydcsinc2_ompc, &
        !$omp& dxdetanc2_ompc, &
        !$omp& dydetanc2_ompc)
        do k = 1,nz
          do j = 1,ny
            do i = 1,ng
              w_ompc(nx+i,j,k,1) = wbuf2r_ompc(i,j,k,1)
              w_ompc(nx+i,j,k,2) = wbuf2r_ompc(i,j,k,2)*dxdcsinc2_ompc(nx+i,j)+wbuf2r_ompc(i,j,k,3)*dxdetanc2_ompc(nx+i,j)
              w_ompc(nx+i,j,k,3) = wbuf2r_ompc(i,j,k,2)*dydcsinc2_ompc(nx+i,j)+wbuf2r_ompc(i,j,k,3)*dydetanc2_ompc(nx+i,j)
              w_ompc(nx+i,j,k,4) = wbuf2r_ompc(i,j,k,4)
              w_ompc(nx+i,j,k,5) = wbuf2r_ompc(i,j,k,5)
            enddo
          enddo
        enddo
      endif
    else
      if (ileftx/=mpi_proc_null) then
        !$omp parallel do default(firstprivate) shared(w_ompc, &
        !$omp& wbuf1r_ompc, &
        !$omp& wbuf2r_ompc, &
        !$omp& wbuf3r_ompc, &
        !$omp& wbuf4r_ompc, &
        !$omp& wbuf5r_ompc, &
        !$omp& wbuf6r_ompc, &
        !$omp& dxdcsinc2_ompc, &
        !$omp& dydcsinc2_ompc, &
        !$omp& dxdetanc2_ompc, &
        !$omp& dydetanc2_ompc)
        do k = 1,nz
          do j = 1,ny
            do i = 1,ng
              do m=1,nv
                w_ompc(i-ng,j,k,m) = wbuf1r_ompc(i,j,k,m)
              enddo
            enddo
          enddo
        enddo
      endif
      if (irightx/=mpi_proc_null) then
        !$omp parallel do default(firstprivate) shared(w_ompc, &
        !$omp& wbuf1r_ompc, &
        !$omp& wbuf2r_ompc, &
        !$omp& wbuf3r_ompc, &
        !$omp& wbuf4r_ompc, &
        !$omp& wbuf5r_ompc, &
        !$omp& wbuf6r_ompc, &
        !$omp& dxdcsinc2_ompc, &
        !$omp& dydcsinc2_ompc, &
        !$omp& dxdetanc2_ompc, &
        !$omp& dydetanc2_ompc)
        do k = 1,nz
          do j = 1,ny
            do i = 1,ng
              do m=1,nv
                w_ompc(nx+i,j,k,m) = wbuf2r_ompc(i,j,k,m)
              enddo
            enddo
          enddo
        enddo
      endif
    endif
    if (ilefty/=mpi_proc_null) then
      !$omp parallel do default(firstprivate) shared(w_ompc, &
      !$omp& wbuf1r_ompc, &
      !$omp& wbuf2r_ompc, &
      !$omp& wbuf3r_ompc, &
      !$omp& wbuf4r_ompc, &
      !$omp& wbuf5r_ompc, &
      !$omp& wbuf6r_ompc, &
      !$omp& dxdcsinc2_ompc, &
      !$omp& dydcsinc2_ompc, &
      !$omp& dxdetanc2_ompc, &
      !$omp& dydetanc2_ompc)
      do k = 1,nz
        do j = 1,ng
          do i = 1,nx
            do m=1,nv
              w_ompc(i,j-ng,k,m) = wbuf3r_ompc(i,j,k,m)
            enddo
          enddo
        enddo
      enddo
    endif
    if (irighty/=mpi_proc_null) then
      !$omp parallel do default(firstprivate) shared(w_ompc, &
      !$omp& wbuf1r_ompc, &
      !$omp& wbuf2r_ompc, &
      !$omp& wbuf3r_ompc, &
      !$omp& wbuf4r_ompc, &
      !$omp& wbuf5r_ompc, &
      !$omp& wbuf6r_ompc, &
      !$omp& dxdcsinc2_ompc, &
      !$omp& dydcsinc2_ompc, &
      !$omp& dxdetanc2_ompc, &
      !$omp& dydetanc2_ompc)
      do k = 1,nz
        do j = 1,ng
          do i = 1,nx
            do m=1,nv
              w_ompc(i,ny+j,k,m) = wbuf4r_ompc(i,j,k,m)
            enddo
          enddo
        enddo
      enddo
    endif
    if(ndim == 3) then
      if (ileftz/=mpi_proc_null) then
        !$omp parallel do default(firstprivate) shared(w_ompc, &
        !$omp& wbuf1r_ompc, &
        !$omp& wbuf2r_ompc, &
        !$omp& wbuf3r_ompc, &
        !$omp& wbuf4r_ompc, &
        !$omp& wbuf5r_ompc, &
        !$omp& wbuf6r_ompc, &
        !$omp& dxdcsinc2_ompc, &
        !$omp& dydcsinc2_ompc, &
        !$omp& dxdetanc2_ompc, &
        !$omp& dydetanc2_ompc)
        do k = 1,ng
          do j = 1,ny
            do i = 1,nx
              do m=1,nv
                w_ompc(i,j,k-ng,m) = wbuf5r_ompc(i,j,k,m)
              enddo
            enddo
          enddo
        enddo
      endif
      if (irightz/=mpi_proc_null) then
        !$omp parallel do default(firstprivate) shared(w_ompc, &
        !$omp& wbuf1r_ompc, &
        !$omp& wbuf2r_ompc, &
        !$omp& wbuf3r_ompc, &
        !$omp& wbuf4r_ompc, &
        !$omp& wbuf5r_ompc, &
        !$omp& wbuf6r_ompc, &
        !$omp& dxdcsinc2_ompc, &
        !$omp& dydcsinc2_ompc, &
        !$omp& dxdetanc2_ompc, &
        !$omp& dydetanc2_ompc)
        do k = 1,ng
          do j = 1,ny
            do i = 1,nx
              do m=1,nv
                w_ompc(i,j,nz+k,m) = wbuf6r_ompc(i,j,k,m)
              enddo
            enddo
          enddo
        enddo
      endif
    endif

  endsubroutine bcswap_c2_step_3_subroutine



  subroutine bcswap_c2xybc_step_3_subroutine(nx,ny,nz,nv,ng,ndim,w_ompc,wbuf1xybcr_ompc,&
  &wbuf2xybcr_ompc,wbuf3r_ompc,wbuf4r_ompc,wbuf5r_ompc,wbuf6r_ompc,ileftx,ilefty,ileftz,irightx,&
  &irighty,irightz,dxdcsinc2_ompc,dydcsinc2_ompc,dxdetanc2_ompc,dydetanc2_ompc,is_periodic)
    integer :: nx, ny, nz, ng, nv, ndim
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: w_ompc
    real(rkind), dimension(1:,1-ng:,1:,1:) :: wbuf1xybcr_ompc, wbuf2xybcr_ompc
    real(rkind), dimension(1:,1:,1:,1:) :: wbuf3r_ompc, wbuf4r_ompc, wbuf5r_ompc, wbuf6r_ompc
    real(rkind), dimension(1-ng:,1-ng:) :: dxdcsinc2_ompc, dydcsinc2_ompc, dxdetanc2_ompc, dydetanc2_ompc
    integer :: i,j,k,m, iercuda
    integer :: ileftx, ilefty, ileftz, irightx, irighty, irightz
    integer, dimension(3) :: is_periodic
    if (is_periodic(1) == 1) then
      if (ileftx/=mpi_proc_null) then
        !$omp parallel do default(firstprivate) shared(w_ompc, &
        !$omp& wbuf1xybcr_ompc, &
        !$omp& wbuf2xybcr_ompc, &
        !$omp& wbuf3r_ompc, &
        !$omp& wbuf4r_ompc, &
        !$omp& wbuf5r_ompc, &
        !$omp& wbuf6r_ompc, &
        !$omp& dxdcsinc2_ompc, &
        !$omp& dydcsinc2_ompc, &
        !$omp& dxdetanc2_ompc, &
        !$omp& dydetanc2_ompc)
        do k = 1,nz
          do j = 1-ng,ny+ng
            do i = 1,ng
              w_ompc(i-ng,j,k,1) = wbuf1xybcr_ompc(i,j,k,1)
              w_ompc(i-ng,j,k,2) = wbuf1xybcr_ompc(i,j,k,2)*dxdcsinc2_ompc(i-ng,j)+wbuf1xybcr_ompc(i,j,k,3)*dxdetanc2_ompc(i-ng,j)
              w_ompc(i-ng,j,k,3) = wbuf1xybcr_ompc(i,j,k,2)*dydcsinc2_ompc(i-ng,j)+wbuf1xybcr_ompc(i,j,k,3)*dydetanc2_ompc(i-ng,j)
              w_ompc(i-ng,j,k,4) = wbuf1xybcr_ompc(i,j,k,4)
              w_ompc(i-ng,j,k,5) = wbuf1xybcr_ompc(i,j,k,5)
            enddo
          enddo
        enddo
      endif
      if (irightx/=mpi_proc_null) then
        !$omp parallel do default(firstprivate) shared(w_ompc, &
        !$omp& wbuf1xybcr_ompc, &
        !$omp& wbuf2xybcr_ompc, &
        !$omp& wbuf3r_ompc, &
        !$omp& wbuf4r_ompc, &
        !$omp& wbuf5r_ompc, &
        !$omp& wbuf6r_ompc, &
        !$omp& dxdcsinc2_ompc, &
        !$omp& dydcsinc2_ompc, &
        !$omp& dxdetanc2_ompc, &
        !$omp& dydetanc2_ompc)
        do k = 1,nz
          do j = 1-ng,ny+ng
            do i = 1,ng
              w_ompc(nx+i,j,k,1) = wbuf2xybcr_ompc(i,j,k,1)
              w_ompc(nx+i,j,k,2) = wbuf2xybcr_ompc(i,j,k,2)*dxdcsinc2_ompc(nx+i,j)+wbuf2xybcr_ompc(i,j,k,3)*dxdetanc2_ompc(nx+i,j)
              w_ompc(nx+i,j,k,3) = wbuf2xybcr_ompc(i,j,k,2)*dydcsinc2_ompc(nx+i,j)+wbuf2xybcr_ompc(i,j,k,3)*dydetanc2_ompc(nx+i,j)
              w_ompc(nx+i,j,k,4) = wbuf2xybcr_ompc(i,j,k,4)
              w_ompc(nx+i,j,k,5) = wbuf2xybcr_ompc(i,j,k,5)
            enddo
          enddo
        enddo
      endif
    else
      if (ileftx/=mpi_proc_null) then
        !$omp parallel do default(firstprivate) shared(w_ompc, &
        !$omp& wbuf1xybcr_ompc, &
        !$omp& wbuf2xybcr_ompc, &
        !$omp& wbuf3r_ompc, &
        !$omp& wbuf4r_ompc, &
        !$omp& wbuf5r_ompc, &
        !$omp& wbuf6r_ompc, &
        !$omp& dxdcsinc2_ompc, &
        !$omp& dydcsinc2_ompc, &
        !$omp& dxdetanc2_ompc, &
        !$omp& dydetanc2_ompc)
        do k = 1,nz
          do j = 1-ng,ny+ng
            do i = 1,ng
              do m=1,nv
                w_ompc(i-ng,j,k,m) = wbuf1xybcr_ompc(i,j,k,m)
              enddo
            enddo
          enddo
        enddo
      endif
      if (irightx/=mpi_proc_null) then
        !$omp parallel do default(firstprivate) shared(w_ompc, &
        !$omp& wbuf1xybcr_ompc, &
        !$omp& wbuf2xybcr_ompc, &
        !$omp& wbuf3r_ompc, &
        !$omp& wbuf4r_ompc, &
        !$omp& wbuf5r_ompc, &
        !$omp& wbuf6r_ompc, &
        !$omp& dxdcsinc2_ompc, &
        !$omp& dydcsinc2_ompc, &
        !$omp& dxdetanc2_ompc, &
        !$omp& dydetanc2_ompc)
        do k = 1,nz
          do j = 1-ng,ny+ng
            do i = 1,ng
              do m=1,nv
                w_ompc(nx+i,j,k,m) = wbuf2xybcr_ompc(i,j,k,m)
              enddo
            enddo
          enddo
        enddo
      endif
    endif
    if (ilefty/=mpi_proc_null) then
      !$omp parallel do default(firstprivate) shared(w_ompc, &
      !$omp& wbuf1xybcr_ompc, &
      !$omp& wbuf2xybcr_ompc, &
      !$omp& wbuf3r_ompc, &
      !$omp& wbuf4r_ompc, &
      !$omp& wbuf5r_ompc, &
      !$omp& wbuf6r_ompc, &
      !$omp& dxdcsinc2_ompc, &
      !$omp& dydcsinc2_ompc, &
      !$omp& dxdetanc2_ompc, &
      !$omp& dydetanc2_ompc)
      do k = 1,nz
        do j = 1,ng
          do i = 1,nx
            do m=1,nv
              w_ompc(i,j-ng,k,m) = wbuf3r_ompc(i,j,k,m)
            enddo
          enddo
        enddo
      enddo
    endif
    if (irighty/=mpi_proc_null) then
      !$omp parallel do default(firstprivate) shared(w_ompc, &
      !$omp& wbuf1xybcr_ompc, &
      !$omp& wbuf2xybcr_ompc, &
      !$omp& wbuf3r_ompc, &
      !$omp& wbuf4r_ompc, &
      !$omp& wbuf5r_ompc, &
      !$omp& wbuf6r_ompc, &
      !$omp& dxdcsinc2_ompc, &
      !$omp& dydcsinc2_ompc, &
      !$omp& dxdetanc2_ompc, &
      !$omp& dydetanc2_ompc)
      do k = 1,nz
        do j = 1,ng
          do i = 1,nx
            do m=1,nv
              w_ompc(i,ny+j,k,m) = wbuf4r_ompc(i,j,k,m)
            enddo
          enddo
        enddo
      enddo
    endif
    if(ndim == 3) then
      if (ileftz/=mpi_proc_null) then
        !$omp parallel do default(firstprivate) shared(w_ompc, &
        !$omp& wbuf1xybcr_ompc, &
        !$omp& wbuf2xybcr_ompc, &
        !$omp& wbuf3r_ompc, &
        !$omp& wbuf4r_ompc, &
        !$omp& wbuf5r_ompc, &
        !$omp& wbuf6r_ompc, &
        !$omp& dxdcsinc2_ompc, &
        !$omp& dydcsinc2_ompc, &
        !$omp& dxdetanc2_ompc, &
        !$omp& dydetanc2_ompc)
        do k = 1,ng
          do j = 1,ny
            do i = 1,nx
              do m=1,nv
                w_ompc(i,j,k-ng,m) = wbuf5r_ompc(i,j,k,m)
              enddo
            enddo
          enddo
        enddo
      endif
      if (irightz/=mpi_proc_null) then
        !$omp parallel do default(firstprivate) shared(w_ompc, &
        !$omp& wbuf1xybcr_ompc, &
        !$omp& wbuf2xybcr_ompc, &
        !$omp& wbuf3r_ompc, &
        !$omp& wbuf4r_ompc, &
        !$omp& wbuf5r_ompc, &
        !$omp& wbuf6r_ompc, &
        !$omp& dxdcsinc2_ompc, &
        !$omp& dydcsinc2_ompc, &
        !$omp& dxdetanc2_ompc, &
        !$omp& dydetanc2_ompc)
        do k = 1,ng
          do j = 1,ny
            do i = 1,nx
              do m=1,nv
                w_ompc(i,j,nz+k,m) = wbuf6r_ompc(i,j,k,m)
              enddo
            enddo
          enddo
        enddo
      endif
    endif

  endsubroutine bcswap_c2xybc_step_3_subroutine



  subroutine bcswap_corner_step_3_subroutine(nx,ny,nz,nv,ng,w_ompc,wbuf1r_c_ompc,wbuf2r_c_ompc,&
  &wbuf3r_c_ompc,wbuf4r_c_ompc,ileftbottom,ilefttop,irightbottom,irighttop)
    integer :: nx, ny, nz, ng, nv
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: w_ompc
    real(rkind), dimension(1:,1:,1:,1:) :: wbuf1r_c_ompc, wbuf2r_c_ompc, wbuf3r_c_ompc, wbuf4r_c_ompc
    integer :: i,j,k,m, iercuda
    integer :: ileftbottom, ilefttop, irightbottom, irighttop
    if (ileftbottom/=mpi_proc_null) then
      !$omp parallel do default(firstprivate) shared(w_ompc, &
      !$omp& wbuf1r_c_ompc, &
      !$omp& wbuf2r_c_ompc, &
      !$omp& wbuf3r_c_ompc, &
      !$omp& wbuf4r_c_ompc)
      do j = 1,ny
        do m = 1,nv
          do k=1,ng
            do i=1,ng
              w_ompc(i-ng,j,k-ng,m) = wbuf1r_c_ompc(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif
    if (irighttop/=mpi_proc_null) then
      !$omp parallel do default(firstprivate) shared(w_ompc, &
      !$omp& wbuf1r_c_ompc, &
      !$omp& wbuf2r_c_ompc, &
      !$omp& wbuf3r_c_ompc, &
      !$omp& wbuf4r_c_ompc)
      do j = 1,ny
        do m = 1,nv
          do k=1,ng
            do i=1,ng
              w_ompc(i+nx,j,k+nz,m) = wbuf2r_c_ompc(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif
    if (ilefttop/=mpi_proc_null) then
      !$omp parallel do default(firstprivate) shared(w_ompc, &
      !$omp& wbuf1r_c_ompc, &
      !$omp& wbuf2r_c_ompc, &
      !$omp& wbuf3r_c_ompc, &
      !$omp& wbuf4r_c_ompc)
      do j = 1,ny
        do m = 1,nv
          do k=1,ng
            do i=1,ng
              w_ompc(i-ng,j,k+nz,m) = wbuf3r_c_ompc(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif
    if (irightbottom/=mpi_proc_null) then
      !$omp parallel do default(firstprivate) shared(w_ompc, &
      !$omp& wbuf1r_c_ompc, &
      !$omp& wbuf2r_c_ompc, &
      !$omp& wbuf3r_c_ompc, &
      !$omp& wbuf4r_c_ompc)
      do j = 1,ny
        do m = 1,nv
          do k=1,ng
            do i=1,ng
              w_ompc(i+nx,j,k-ng,m) = wbuf4r_c_ompc(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif

  endsubroutine bcswap_corner_step_3_subroutine



  subroutine bcswap_corner_step_3b_subroutine(nx,ny,nz,nv,ng,w_ompc,wbuf_lxly_r_ompc,&
  &wbuf_lxry_r_ompc,wbuf_rxly_r_ompc,wbuf_rxry_r_ompc,wbuf_lylz_r_ompc,wbuf_lyrz_r_ompc,&
  &wbuf_rylz_r_ompc,wbuf_ryrz_r_ompc,wbuf_lxlylz_r_ompc,wbuf_lxlyrz_r_ompc,wbuf_lxrylz_r_ompc,&
  &wbuf_lxryrz_r_ompc,wbuf_rxlylz_r_ompc,wbuf_rxlyrz_r_ompc,wbuf_rxrylz_r_ompc,wbuf_rxryrz_r_ompc,lxly,&
  &lxry,rxly,rxry,lylz,lyrz,rylz,ryrz,lxlylz,lxlyrz,lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz)
    integer :: nx, ny, nz, ng, nv
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: w_ompc
    real(rkind), dimension(1:,1:,1:,1:) :: wbuf_lxly_r_ompc, wbuf_lxry_r_ompc, wbuf_rxly_r_ompc, wbuf_rxry_r_ompc
    real(rkind), dimension(1:,1:,1:,1:) :: wbuf_lylz_r_ompc, wbuf_lyrz_r_ompc, wbuf_rylz_r_ompc, wbuf_ryrz_r_ompc
    real(rkind), dimension(1:,1:,1:,1:) :: wbuf_lxlylz_r_ompc, wbuf_lxlyrz_r_ompc, wbuf_lxrylz_r_ompc, wbuf_lxryrz_r_ompc
    real(rkind), dimension(1:,1:,1:,1:) :: wbuf_rxlylz_r_ompc, wbuf_rxlyrz_r_ompc, wbuf_rxrylz_r_ompc, wbuf_rxryrz_r_ompc
    integer :: i,j,k,m, iercuda
    integer :: lxly, lxry, rxly, rxry
    integer :: lylz, lyrz, rylz, ryrz
    integer :: lxlylz, lxlyrz, lxrylz, lxryrz
    integer :: rxlylz, rxlyrz, rxrylz, rxryrz
    if (lxly/=mpi_proc_null) then
      !$omp parallel do default(firstprivate) shared(w_ompc, &
      !$omp& wbuf_lxly_r_ompc, &
      !$omp& wbuf_lxry_r_ompc, &
      !$omp& wbuf_rxly_r_ompc, &
      !$omp& wbuf_rxry_r_ompc, &
      !$omp& wbuf_lylz_r_ompc, &
      !$omp& wbuf_lyrz_r_ompc, &
      !$omp& wbuf_rylz_r_ompc, &
      !$omp& wbuf_ryrz_r_ompc, &
      !$omp& wbuf_lxlylz_r_ompc, &
      !$omp& wbuf_lxlyrz_r_ompc, &
      !$omp& wbuf_lxrylz_r_ompc, &
      !$omp& wbuf_lxryrz_r_ompc, &
      !$omp& wbuf_rxlylz_r_ompc, &
      !$omp& wbuf_rxlyrz_r_ompc, &
      !$omp& wbuf_rxrylz_r_ompc, &
      !$omp& wbuf_rxryrz_r_ompc)
      do k = 1,nz
        do m = 1,nv
          do i=1,ng
            do j=1,ng
              w_ompc(i-ng,j-ng,k,m) = wbuf_lxly_r_ompc(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif
    if (rxry/=mpi_proc_null) then
      !$omp parallel do default(firstprivate) shared(w_ompc, &
      !$omp& wbuf_lxly_r_ompc, &
      !$omp& wbuf_lxry_r_ompc, &
      !$omp& wbuf_rxly_r_ompc, &
      !$omp& wbuf_rxry_r_ompc, &
      !$omp& wbuf_lylz_r_ompc, &
      !$omp& wbuf_lyrz_r_ompc, &
      !$omp& wbuf_rylz_r_ompc, &
      !$omp& wbuf_ryrz_r_ompc, &
      !$omp& wbuf_lxlylz_r_ompc, &
      !$omp& wbuf_lxlyrz_r_ompc, &
      !$omp& wbuf_lxrylz_r_ompc, &
      !$omp& wbuf_lxryrz_r_ompc, &
      !$omp& wbuf_rxlylz_r_ompc, &
      !$omp& wbuf_rxlyrz_r_ompc, &
      !$omp& wbuf_rxrylz_r_ompc, &
      !$omp& wbuf_rxryrz_r_ompc)
      do k = 1,nz
        do m = 1,nv
          do i=1,ng
            do j=1,ng
              w_ompc(i+nx,j+ny,k,m) = wbuf_rxry_r_ompc(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif
    if (lxry/=mpi_proc_null) then
      !$omp parallel do default(firstprivate) shared(w_ompc, &
      !$omp& wbuf_lxly_r_ompc, &
      !$omp& wbuf_lxry_r_ompc, &
      !$omp& wbuf_rxly_r_ompc, &
      !$omp& wbuf_rxry_r_ompc, &
      !$omp& wbuf_lylz_r_ompc, &
      !$omp& wbuf_lyrz_r_ompc, &
      !$omp& wbuf_rylz_r_ompc, &
      !$omp& wbuf_ryrz_r_ompc, &
      !$omp& wbuf_lxlylz_r_ompc, &
      !$omp& wbuf_lxlyrz_r_ompc, &
      !$omp& wbuf_lxrylz_r_ompc, &
      !$omp& wbuf_lxryrz_r_ompc, &
      !$omp& wbuf_rxlylz_r_ompc, &
      !$omp& wbuf_rxlyrz_r_ompc, &
      !$omp& wbuf_rxrylz_r_ompc, &
      !$omp& wbuf_rxryrz_r_ompc)
      do k = 1,nz
        do m = 1,nv
          do i=1,ng
            do j=1,ng
              w_ompc(i-ng,j+ny,k,m) = wbuf_lxry_r_ompc(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif
    if (rxly/=mpi_proc_null) then
      !$omp parallel do default(firstprivate) shared(w_ompc, &
      !$omp& wbuf_lxly_r_ompc, &
      !$omp& wbuf_lxry_r_ompc, &
      !$omp& wbuf_rxly_r_ompc, &
      !$omp& wbuf_rxry_r_ompc, &
      !$omp& wbuf_lylz_r_ompc, &
      !$omp& wbuf_lyrz_r_ompc, &
      !$omp& wbuf_rylz_r_ompc, &
      !$omp& wbuf_ryrz_r_ompc, &
      !$omp& wbuf_lxlylz_r_ompc, &
      !$omp& wbuf_lxlyrz_r_ompc, &
      !$omp& wbuf_lxrylz_r_ompc, &
      !$omp& wbuf_lxryrz_r_ompc, &
      !$omp& wbuf_rxlylz_r_ompc, &
      !$omp& wbuf_rxlyrz_r_ompc, &
      !$omp& wbuf_rxrylz_r_ompc, &
      !$omp& wbuf_rxryrz_r_ompc)
      do k = 1,nz
        do m = 1,nv
          do i=1,ng
            do j=1,ng
              w_ompc(i+nx,j-ng,k,m) = wbuf_rxly_r_ompc(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif

    if (lylz/=mpi_proc_null) then
      !$omp parallel do default(firstprivate) shared(w_ompc, &
      !$omp& wbuf_lxly_r_ompc, &
      !$omp& wbuf_lxry_r_ompc, &
      !$omp& wbuf_rxly_r_ompc, &
      !$omp& wbuf_rxry_r_ompc, &
      !$omp& wbuf_lylz_r_ompc, &
      !$omp& wbuf_lyrz_r_ompc, &
      !$omp& wbuf_rylz_r_ompc, &
      !$omp& wbuf_ryrz_r_ompc, &
      !$omp& wbuf_lxlylz_r_ompc, &
      !$omp& wbuf_lxlyrz_r_ompc, &
      !$omp& wbuf_lxrylz_r_ompc, &
      !$omp& wbuf_lxryrz_r_ompc, &
      !$omp& wbuf_rxlylz_r_ompc, &
      !$omp& wbuf_rxlyrz_r_ompc, &
      !$omp& wbuf_rxrylz_r_ompc, &
      !$omp& wbuf_rxryrz_r_ompc)
      do i = 1,nx
        do m = 1,nv
          do k=1,ng
            do j=1,ng
              w_ompc(i,j-ng,k-ng,m) = wbuf_lylz_r_ompc(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif
    if (ryrz/=mpi_proc_null) then
      !$omp parallel do default(firstprivate) shared(w_ompc, &
      !$omp& wbuf_lxly_r_ompc, &
      !$omp& wbuf_lxry_r_ompc, &
      !$omp& wbuf_rxly_r_ompc, &
      !$omp& wbuf_rxry_r_ompc, &
      !$omp& wbuf_lylz_r_ompc, &
      !$omp& wbuf_lyrz_r_ompc, &
      !$omp& wbuf_rylz_r_ompc, &
      !$omp& wbuf_ryrz_r_ompc, &
      !$omp& wbuf_lxlylz_r_ompc, &
      !$omp& wbuf_lxlyrz_r_ompc, &
      !$omp& wbuf_lxrylz_r_ompc, &
      !$omp& wbuf_lxryrz_r_ompc, &
      !$omp& wbuf_rxlylz_r_ompc, &
      !$omp& wbuf_rxlyrz_r_ompc, &
      !$omp& wbuf_rxrylz_r_ompc, &
      !$omp& wbuf_rxryrz_r_ompc)
      do i = 1,nx
        do m = 1,nv
          do k=1,ng
            do j=1,ng
              w_ompc(i,j+ny,k+nz,m) = wbuf_ryrz_r_ompc(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif
    if (lyrz/=mpi_proc_null) then
      !$omp parallel do default(firstprivate) shared(w_ompc, &
      !$omp& wbuf_lxly_r_ompc, &
      !$omp& wbuf_lxry_r_ompc, &
      !$omp& wbuf_rxly_r_ompc, &
      !$omp& wbuf_rxry_r_ompc, &
      !$omp& wbuf_lylz_r_ompc, &
      !$omp& wbuf_lyrz_r_ompc, &
      !$omp& wbuf_rylz_r_ompc, &
      !$omp& wbuf_ryrz_r_ompc, &
      !$omp& wbuf_lxlylz_r_ompc, &
      !$omp& wbuf_lxlyrz_r_ompc, &
      !$omp& wbuf_lxrylz_r_ompc, &
      !$omp& wbuf_lxryrz_r_ompc, &
      !$omp& wbuf_rxlylz_r_ompc, &
      !$omp& wbuf_rxlyrz_r_ompc, &
      !$omp& wbuf_rxrylz_r_ompc, &
      !$omp& wbuf_rxryrz_r_ompc)
      do i = 1,nx
        do m = 1,nv
          do k=1,ng
            do j=1,ng
              w_ompc(i,j-ng,k+nz,m) = wbuf_lyrz_r_ompc(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif
    if (rylz/=mpi_proc_null) then
      !$omp parallel do default(firstprivate) shared(w_ompc, &
      !$omp& wbuf_lxly_r_ompc, &
      !$omp& wbuf_lxry_r_ompc, &
      !$omp& wbuf_rxly_r_ompc, &
      !$omp& wbuf_rxry_r_ompc, &
      !$omp& wbuf_lylz_r_ompc, &
      !$omp& wbuf_lyrz_r_ompc, &
      !$omp& wbuf_rylz_r_ompc, &
      !$omp& wbuf_ryrz_r_ompc, &
      !$omp& wbuf_lxlylz_r_ompc, &
      !$omp& wbuf_lxlyrz_r_ompc, &
      !$omp& wbuf_lxrylz_r_ompc, &
      !$omp& wbuf_lxryrz_r_ompc, &
      !$omp& wbuf_rxlylz_r_ompc, &
      !$omp& wbuf_rxlyrz_r_ompc, &
      !$omp& wbuf_rxrylz_r_ompc, &
      !$omp& wbuf_rxryrz_r_ompc)
      do i = 1,nx
        do m = 1,nv
          do k=1,ng
            do j=1,ng
              w_ompc(i,j+ny,k-ng,m) = wbuf_rylz_r_ompc(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif

    if (lxlylz/=mpi_proc_null) then
      !$omp parallel do default(firstprivate) shared(w_ompc, &
      !$omp& wbuf_lxly_r_ompc, &
      !$omp& wbuf_lxry_r_ompc, &
      !$omp& wbuf_rxly_r_ompc, &
      !$omp& wbuf_rxry_r_ompc, &
      !$omp& wbuf_lylz_r_ompc, &
      !$omp& wbuf_lyrz_r_ompc, &
      !$omp& wbuf_rylz_r_ompc, &
      !$omp& wbuf_ryrz_r_ompc, &
      !$omp& wbuf_lxlylz_r_ompc, &
      !$omp& wbuf_lxlyrz_r_ompc, &
      !$omp& wbuf_lxrylz_r_ompc, &
      !$omp& wbuf_lxryrz_r_ompc, &
      !$omp& wbuf_rxlylz_r_ompc, &
      !$omp& wbuf_rxlyrz_r_ompc, &
      !$omp& wbuf_rxrylz_r_ompc, &
      !$omp& wbuf_rxryrz_r_ompc)
      do k = 1,ng
        do m = 1,nv
          do i=1,ng
            do j=1,ng
              w_ompc(i-ng,j-ng,k-ng,m) = wbuf_lxlylz_r_ompc(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif
    if (lxlyrz/=mpi_proc_null) then
      !$omp parallel do default(firstprivate) shared(w_ompc, &
      !$omp& wbuf_lxly_r_ompc, &
      !$omp& wbuf_lxry_r_ompc, &
      !$omp& wbuf_rxly_r_ompc, &
      !$omp& wbuf_rxry_r_ompc, &
      !$omp& wbuf_lylz_r_ompc, &
      !$omp& wbuf_lyrz_r_ompc, &
      !$omp& wbuf_rylz_r_ompc, &
      !$omp& wbuf_ryrz_r_ompc, &
      !$omp& wbuf_lxlylz_r_ompc, &
      !$omp& wbuf_lxlyrz_r_ompc, &
      !$omp& wbuf_lxrylz_r_ompc, &
      !$omp& wbuf_lxryrz_r_ompc, &
      !$omp& wbuf_rxlylz_r_ompc, &
      !$omp& wbuf_rxlyrz_r_ompc, &
      !$omp& wbuf_rxrylz_r_ompc, &
      !$omp& wbuf_rxryrz_r_ompc)
      do k = 1,ng
        do m = 1,nv
          do i=1,ng
            do j=1,ng
              w_ompc(i-ng,j-ng,k+nz,m) = wbuf_lxlyrz_r_ompc(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif
    if (lxrylz/=mpi_proc_null) then
      !$omp parallel do default(firstprivate) shared(w_ompc, &
      !$omp& wbuf_lxly_r_ompc, &
      !$omp& wbuf_lxry_r_ompc, &
      !$omp& wbuf_rxly_r_ompc, &
      !$omp& wbuf_rxry_r_ompc, &
      !$omp& wbuf_lylz_r_ompc, &
      !$omp& wbuf_lyrz_r_ompc, &
      !$omp& wbuf_rylz_r_ompc, &
      !$omp& wbuf_ryrz_r_ompc, &
      !$omp& wbuf_lxlylz_r_ompc, &
      !$omp& wbuf_lxlyrz_r_ompc, &
      !$omp& wbuf_lxrylz_r_ompc, &
      !$omp& wbuf_lxryrz_r_ompc, &
      !$omp& wbuf_rxlylz_r_ompc, &
      !$omp& wbuf_rxlyrz_r_ompc, &
      !$omp& wbuf_rxrylz_r_ompc, &
      !$omp& wbuf_rxryrz_r_ompc)
      do k = 1,ng
        do m = 1,nv
          do i=1,ng
            do j=1,ng
              w_ompc(i-ng,j+ny,k-ng,m) = wbuf_lxrylz_r_ompc(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif
    if (lxryrz/=mpi_proc_null) then
      !$omp parallel do default(firstprivate) shared(w_ompc, &
      !$omp& wbuf_lxly_r_ompc, &
      !$omp& wbuf_lxry_r_ompc, &
      !$omp& wbuf_rxly_r_ompc, &
      !$omp& wbuf_rxry_r_ompc, &
      !$omp& wbuf_lylz_r_ompc, &
      !$omp& wbuf_lyrz_r_ompc, &
      !$omp& wbuf_rylz_r_ompc, &
      !$omp& wbuf_ryrz_r_ompc, &
      !$omp& wbuf_lxlylz_r_ompc, &
      !$omp& wbuf_lxlyrz_r_ompc, &
      !$omp& wbuf_lxrylz_r_ompc, &
      !$omp& wbuf_lxryrz_r_ompc, &
      !$omp& wbuf_rxlylz_r_ompc, &
      !$omp& wbuf_rxlyrz_r_ompc, &
      !$omp& wbuf_rxrylz_r_ompc, &
      !$omp& wbuf_rxryrz_r_ompc)
      do k = 1,ng
        do m = 1,nv
          do i=1,ng
            do j=1,ng
              w_ompc(i-ng,j+ny,k+nz,m) = wbuf_lxryrz_r_ompc(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif
    if (rxlylz/=mpi_proc_null) then
      !$omp parallel do default(firstprivate) shared(w_ompc, &
      !$omp& wbuf_lxly_r_ompc, &
      !$omp& wbuf_lxry_r_ompc, &
      !$omp& wbuf_rxly_r_ompc, &
      !$omp& wbuf_rxry_r_ompc, &
      !$omp& wbuf_lylz_r_ompc, &
      !$omp& wbuf_lyrz_r_ompc, &
      !$omp& wbuf_rylz_r_ompc, &
      !$omp& wbuf_ryrz_r_ompc, &
      !$omp& wbuf_lxlylz_r_ompc, &
      !$omp& wbuf_lxlyrz_r_ompc, &
      !$omp& wbuf_lxrylz_r_ompc, &
      !$omp& wbuf_lxryrz_r_ompc, &
      !$omp& wbuf_rxlylz_r_ompc, &
      !$omp& wbuf_rxlyrz_r_ompc, &
      !$omp& wbuf_rxrylz_r_ompc, &
      !$omp& wbuf_rxryrz_r_ompc)
      do k = 1,ng
        do m = 1,nv
          do i=1,ng
            do j=1,ng
              w_ompc(i+nx,j-ng,k-ng,m) = wbuf_rxlylz_r_ompc(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif
    if (rxlyrz/=mpi_proc_null) then
      !$omp parallel do default(firstprivate) shared(w_ompc, &
      !$omp& wbuf_lxly_r_ompc, &
      !$omp& wbuf_lxry_r_ompc, &
      !$omp& wbuf_rxly_r_ompc, &
      !$omp& wbuf_rxry_r_ompc, &
      !$omp& wbuf_lylz_r_ompc, &
      !$omp& wbuf_lyrz_r_ompc, &
      !$omp& wbuf_rylz_r_ompc, &
      !$omp& wbuf_ryrz_r_ompc, &
      !$omp& wbuf_lxlylz_r_ompc, &
      !$omp& wbuf_lxlyrz_r_ompc, &
      !$omp& wbuf_lxrylz_r_ompc, &
      !$omp& wbuf_lxryrz_r_ompc, &
      !$omp& wbuf_rxlylz_r_ompc, &
      !$omp& wbuf_rxlyrz_r_ompc, &
      !$omp& wbuf_rxrylz_r_ompc, &
      !$omp& wbuf_rxryrz_r_ompc)
      do k = 1,ng
        do m = 1,nv
          do i=1,ng
            do j=1,ng
              w_ompc(i+nx,j-ng,k+nz,m) = wbuf_rxlyrz_r_ompc(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif
    if (rxrylz/=mpi_proc_null) then
      !$omp parallel do default(firstprivate) shared(w_ompc, &
      !$omp& wbuf_lxly_r_ompc, &
      !$omp& wbuf_lxry_r_ompc, &
      !$omp& wbuf_rxly_r_ompc, &
      !$omp& wbuf_rxry_r_ompc, &
      !$omp& wbuf_lylz_r_ompc, &
      !$omp& wbuf_lyrz_r_ompc, &
      !$omp& wbuf_rylz_r_ompc, &
      !$omp& wbuf_ryrz_r_ompc, &
      !$omp& wbuf_lxlylz_r_ompc, &
      !$omp& wbuf_lxlyrz_r_ompc, &
      !$omp& wbuf_lxrylz_r_ompc, &
      !$omp& wbuf_lxryrz_r_ompc, &
      !$omp& wbuf_rxlylz_r_ompc, &
      !$omp& wbuf_rxlyrz_r_ompc, &
      !$omp& wbuf_rxrylz_r_ompc, &
      !$omp& wbuf_rxryrz_r_ompc)
      do k = 1,ng
        do m = 1,nv
          do i=1,ng
            do j=1,ng
              w_ompc(i+nx,j+ny,k-ng,m) = wbuf_rxrylz_r_ompc(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif
    if (rxryrz/=mpi_proc_null) then
      !$omp parallel do default(firstprivate) shared(w_ompc, &
      !$omp& wbuf_lxly_r_ompc, &
      !$omp& wbuf_lxry_r_ompc, &
      !$omp& wbuf_rxly_r_ompc, &
      !$omp& wbuf_rxry_r_ompc, &
      !$omp& wbuf_lylz_r_ompc, &
      !$omp& wbuf_lyrz_r_ompc, &
      !$omp& wbuf_rylz_r_ompc, &
      !$omp& wbuf_ryrz_r_ompc, &
      !$omp& wbuf_lxlylz_r_ompc, &
      !$omp& wbuf_lxlyrz_r_ompc, &
      !$omp& wbuf_lxrylz_r_ompc, &
      !$omp& wbuf_lxryrz_r_ompc, &
      !$omp& wbuf_rxlylz_r_ompc, &
      !$omp& wbuf_rxlyrz_r_ompc, &
      !$omp& wbuf_rxrylz_r_ompc, &
      !$omp& wbuf_rxryrz_r_ompc)
      do k = 1,ng
        do m = 1,nv
          do i=1,ng
            do j=1,ng
              w_ompc(i+nx,j+ny,k+nz,m) = wbuf_rxryrz_r_ompc(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif

  endsubroutine bcswap_corner_step_3b_subroutine



  subroutine bcswap_corner(self, steps)
    class(base_ompc_object), intent(inout) :: self
    logical, dimension(3), optional :: steps
    logical, dimension(3) :: steps_
    integer :: iercuda, indc
    integer, dimension(mpi_status_size) :: istatus

    steps_ = .true. ; if(present(steps)) steps_ = steps

    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv,&
    & w_ompc => self%w_ompc, wbuf1s_c_ompc => self%wbuf1s_c_ompc, wbuf2s_c_ompc => self%wbuf2s_c_ompc,&
    & wbuf3s_c_ompc => self%wbuf3s_c_ompc, wbuf4s_c_ompc => self%wbuf4s_c_ompc, wbuf1r_c_ompc => self%wbu&
    &f1r_c_ompc, wbuf2r_c_ompc => self%wbuf2r_c_ompc, wbuf3r_c_ompc => self%wbuf3r_c_ompc,&
    & wbuf4r_c_ompc => self%wbuf4r_c_ompc, ileftbottom => self%field%ileftbottom, irightbottom => self%fi&
    &eld%irightbottom, ilefttop => self%field%ilefttop, irighttop => self%field%irighttop,&
    & mp_cart => self%field%mp_cart, iermpi => self%mpi_err, nrank => self%field%myrank)

      if(steps_(1)) then
        call bcswap_corner_step_1_subroutine(nx, ny, nz, nv, ng, w_ompc, wbuf1s_c_ompc,&
        & wbuf2s_c_ompc, wbuf3s_c_ompc, wbuf4s_c_ompc)
      endif
      if(steps_(2)) then
        indc = nv*ny*ng*ng
        if (ileftbottom == nrank) then
          wbuf2r_c_ompc = wbuf1s_c_ompc
          wbuf1r_c_ompc = wbuf2s_c_ompc
        else
          call mpi_sendrecv(wbuf1s_c_ompc,indc,mpi_prec,ileftbottom ,1, wbuf2r_c_ompc,indc,mpi_prec,&
          &irighttop ,1,mp_cart,istatus,iermpi)
          call mpi_sendrecv(wbuf2s_c_ompc,indc,mpi_prec,irighttop ,2, wbuf1r_c_ompc,indc,mpi_prec,&
          &ileftbottom ,2,mp_cart,istatus,iermpi)
        endif
        if (ilefttop == nrank) then
          wbuf4r_c_ompc = wbuf3s_c_ompc
          wbuf3r_c_ompc = wbuf4s_c_ompc
        else
          call mpi_sendrecv(wbuf3s_c_ompc,indc,mpi_prec,ilefttop ,3, wbuf4r_c_ompc,indc,mpi_prec,&
          &irightbottom,3,mp_cart,istatus,iermpi)
          call mpi_sendrecv(wbuf4s_c_ompc,indc,mpi_prec,irightbottom,4, wbuf3r_c_ompc,indc,mpi_prec,&
          &ilefttop ,4,mp_cart,istatus,iermpi)
        endif
      endif
      if(steps_(3)) then
        call bcswap_corner_step_3_subroutine(nx, ny, nz, nv, ng, w_ompc, wbuf1r_c_ompc,&
        & wbuf2r_c_ompc, wbuf3r_c_ompc, wbuf4r_c_ompc, ileftbottom, ilefttop, irightbottom, irighttop)
      endif
    endassociate
  endsubroutine bcswap_corner

  subroutine bcswap_corner_var(self, w_swap, steps)
    class(base_ompc_object), intent(inout) :: self
    real(rkind), dimension(1-self%ng:self%nx+self%ng, 1-self%ng:self%ny+self%ng, 1-self%ng:self%nz+self%ng,1:1) :: w_swap
    logical, dimension(3), optional :: steps
    logical, dimension(3) :: steps_
    integer :: iercuda, indc
    integer, dimension(mpi_status_size) :: istatus

    steps_ = .true. ; if(present(steps)) steps_ = steps

    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv,&
    & wbuf1s_c_ompc => self%wbuf1s_c_ompc, wbuf2s_c_ompc => self%wbuf2s_c_ompc, wbuf3s_c_ompc => self%wbu&
    &f3s_c_ompc, wbuf4s_c_ompc => self%wbuf4s_c_ompc, wbuf1r_c_ompc => self%wbuf1r_c_ompc,&
    & wbuf2r_c_ompc => self%wbuf2r_c_ompc, wbuf3r_c_ompc => self%wbuf3r_c_ompc, wbuf4r_c_ompc => self%wbu&
    &f4r_c_ompc, ileftbottom => self%field%ileftbottom, irightbottom => self%field%irightbottom,&
    & ilefttop => self%field%ilefttop, irighttop => self%field%irighttop, mp_cart => self%field%mp_cart,&
    & iermpi => self%mpi_err, nrank => self%field%myrank)

      if(steps_(1)) then
        call bcswap_corner_step_1_subroutine(nx, ny, nz, 1, ng, w_swap, wbuf1s_c_ompc, wbuf2s_c_ompc, wbuf3s_c_ompc, wbuf4s_c_ompc)
      endif
      if(steps_(2)) then
        indc = ny*ng*ng
        if (ileftbottom == nrank) then
          wbuf2r_c_ompc = wbuf1s_c_ompc
          wbuf1r_c_ompc = wbuf2s_c_ompc
        else
          call mpi_sendrecv(wbuf1s_c_ompc,indc,mpi_prec,ileftbottom ,1, wbuf2r_c_ompc,indc,mpi_prec,&
          &irighttop ,1,mp_cart,istatus,iermpi)
          call mpi_sendrecv(wbuf2s_c_ompc,indc,mpi_prec,irighttop ,2, wbuf1r_c_ompc,indc,mpi_prec,&
          &ileftbottom ,2,mp_cart,istatus,iermpi)
        endif
        if (ilefttop == nrank) then
          wbuf4r_c_ompc = wbuf3s_c_ompc
          wbuf3r_c_ompc = wbuf4s_c_ompc
        else
          call mpi_sendrecv(wbuf3s_c_ompc,indc,mpi_prec,ilefttop ,3, wbuf4r_c_ompc,indc,mpi_prec,&
          &irightbottom,3,mp_cart,istatus,iermpi)
          call mpi_sendrecv(wbuf4s_c_ompc,indc,mpi_prec,irightbottom,4, wbuf3r_c_ompc,indc,mpi_prec,&
          &ilefttop ,4,mp_cart,istatus,iermpi)
        endif
      endif
      if(steps_(3)) then
        call bcswap_corner_step_3_subroutine(nx, ny, nz, 1, ng, w_swap, wbuf1r_c_ompc,&
        & wbuf2r_c_ompc, wbuf3r_c_ompc, wbuf4r_c_ompc, ileftbottom, ilefttop, irightbottom, irighttop)
      endif
    endassociate
  endsubroutine bcswap_corner_var

  subroutine bcswap(self, steps, swap_xy_corner_bc)
    class(base_ompc_object), intent(inout) :: self
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
    & w_ompc => self%w_ompc, wbuf1s_ompc => self%wbuf1s_ompc, wbuf2s_ompc => self%wbuf2s_ompc,&
    & wbuf1xybcs_ompc => self%wbuf1xybcs_ompc, wbuf2xybcs_ompc => self%wbuf2xybcs_ompc,&
    & wbuf3s_ompc => self%wbuf3s_ompc, wbuf4s_ompc => self%wbuf4s_ompc, wbuf5s_ompc => self%wbuf5s_ompc,&
    & wbuf6s_ompc => self%wbuf6s_ompc, wbuf1r_ompc => self%wbuf1r_ompc, wbuf2r_ompc => self%wbuf2r_ompc,&
    & wbuf1xybcr_ompc => self%wbuf1xybcr_ompc, wbuf2xybcr_ompc => self%wbuf2xybcr_ompc,&
    & wbuf3r_ompc => self%wbuf3r_ompc, wbuf4r_ompc => self%wbuf4r_ompc, wbuf5r_ompc => self%wbuf5r_ompc,&
    & wbuf6r_ompc => self%wbuf6r_ompc, ileftx => self%field%ileftx,irightx => self%field%irightx,&
    & nrank_x => self%field%nrank_x, ilefty => self%field%ilefty,irighty => self%field%irighty,&
    & nrank_y => self%field%nrank_y, ileftz => self%field%ileftz,irightz => self%field%irightz,&
    & nrank_z => self%field%nrank_z, dcsidxnc2_ompc => self%dcsidxnc2_ompc, dcsidync2_ompc => self%dcsidy&
    &nc2_ompc, detadxnc2_ompc => self%detadxnc2_ompc, detadync2_ompc => self%detadync2_ompc,&
    & dxdcsinc2_ompc => self%dxdcsinc2_ompc, dydcsinc2_ompc => self%dydcsinc2_ompc, dxdetanc2_ompc => sel&
    &f%dxdetanc2_ompc, dydetanc2_ompc => self%dydetanc2_ompc, mp_cartx => self%field%mp_cartx,&
    & mp_carty => self%field%mp_carty, mp_cartz => self%field%mp_cartz, ncoords => self%field%ncoords,&
    & nblocks => self%field%nblocks, ndim => self%field%grid%ndim)

      if(steps_(1)) then
        if(self%field%grid%grid_dim == 1) then
          call bcswap_step_1_subroutine(nx, ny, nz, nv, ng, ndim, w_ompc, wbuf1s_ompc, wbuf2s_ompc,&
          & wbuf3s_ompc, wbuf4s_ompc, wbuf5s_ompc, wbuf6s_ompc)
        elseif(self%field%grid%grid_dim == 2) then
          if(swap_xy_corner_bc_ == 0) then
            call bcswap_c2_step_1_subroutine(nx, ny, nz, nv, ng, ndim, w_ompc, wbuf1s_ompc,&
            & wbuf2s_ompc, wbuf3s_ompc, wbuf4s_ompc, wbuf5s_ompc, wbuf6s_ompc, dcsidxnc2_ompc, dcsidync2_ompc,&
            & detadxnc2_ompc, detadync2_ompc, is_periodic)
          else
            call bcswap_c2xybc_step_1_subroutine(nx, ny, nz, nv, ng, ndim, w_ompc, wbuf1xybcs_ompc,&
            & wbuf2xybcs_ompc, wbuf3s_ompc, wbuf4s_ompc, wbuf5s_ompc, wbuf6s_ompc, dcsidxnc2_ompc,&
            & dcsidync2_ompc, detadxnc2_ompc, detadync2_ompc, is_periodic)
          endif
        endif
      endif
      if(steps_(2)) then
        indy = nv*nx*ng*nz
        indz = nv*nx*ny*ng
        if(self%field%grid%grid_dim == 1 .or. (self%field%grid%grid_dim == 2 .and. swap_xy_corner_bc_ == 0)) then
          indx = nv*ng*ny*nz
          if(ileftx == nrank_x) then
            wbuf2r_ompc = wbuf1s_ompc
            wbuf1r_ompc = wbuf2s_ompc
          else
            call mpi_sendrecv(wbuf1s_ompc,indx,mpi_prec,ileftx ,1,wbuf2r_ompc,indx,mpi_prec,irightx,&
            &1,mp_cartx,istatus,self%mpi_err)
            call mpi_sendrecv(wbuf2s_ompc,indx,mpi_prec,irightx,2,wbuf1r_ompc,indx,mpi_prec,ileftx ,&
            &2,mp_cartx,istatus,self%mpi_err)
          endif
        else
          indx = nv*ng*(2*ng+ny)*nz
          if(ileftx == nrank_x) then
            wbuf2xybcr_ompc = wbuf1xybcs_ompc
            wbuf1xybcr_ompc = wbuf2xybcs_ompc
          else
            call mpi_sendrecv(wbuf1xybcs_ompc,indx,mpi_prec,ileftx ,1,wbuf2xybcr_ompc,indx,mpi_prec,&
            &irightx,1,mp_cartx,istatus,self%mpi_err)
            call mpi_sendrecv(wbuf2xybcs_ompc,indx,mpi_prec,irightx,2,wbuf1xybcr_ompc,indx,mpi_prec,&
            &ileftx ,2,mp_cartx,istatus,self%mpi_err)
          endif
        endif
        if(ilefty == nrank_y) then
          wbuf4r_ompc = wbuf3s_ompc
          wbuf3r_ompc = wbuf4s_ompc
        else
          call mpi_sendrecv(wbuf3s_ompc,indy,mpi_prec,ilefty ,3,wbuf4r_ompc,indy,mpi_prec,irighty,3,mp_carty,istatus,self%mpi_err)
          call mpi_sendrecv(wbuf4s_ompc,indy,mpi_prec,irighty,4,wbuf3r_ompc,indy,mpi_prec,ilefty ,4,mp_carty,istatus,self%mpi_err)
        endif
        if(ndim == 3) then
          if(ileftz == nrank_z) then
            wbuf6r_ompc = wbuf5s_ompc
            wbuf5r_ompc = wbuf6s_ompc
          else
            call mpi_sendrecv(wbuf5s_ompc,indz,mpi_prec,ileftz ,5,wbuf6r_ompc,indz,mpi_prec,irightz,&
            &5,mp_cartz,istatus,self%mpi_err)
            call mpi_sendrecv(wbuf6s_ompc,indz,mpi_prec,irightz,6,wbuf5r_ompc,indz,mpi_prec,ileftz ,&
            &6,mp_cartz,istatus,self%mpi_err)
          endif
        endif
      endif
      if(steps_(3)) then
        if(self%field%grid%grid_dim == 1) then
          call bcswap_step_3_subroutine(nx, ny, nz, nv, ng, ndim, w_ompc, wbuf1r_ompc, wbuf2r_ompc,&
          & wbuf3r_ompc, wbuf4r_ompc, wbuf5r_ompc, wbuf6r_ompc, ileftx, ilefty, ileftz, irightx, irighty,&
          & irightz)
        elseif(self%field%grid%grid_dim == 2) then
          if(swap_xy_corner_bc_ == 0) then
            call bcswap_c2_step_3_subroutine(nx, ny, nz, nv, ng, ndim, w_ompc, wbuf1r_ompc,&
            & wbuf2r_ompc, wbuf3r_ompc, wbuf4r_ompc, wbuf5r_ompc, wbuf6r_ompc, ileftx, ilefty, ileftz, irightx,&
            & irighty, irightz, dxdcsinc2_ompc, dydcsinc2_ompc, dxdetanc2_ompc, dydetanc2_ompc, is_periodic)
          else
            call bcswap_c2xybc_step_3_subroutine(nx, ny, nz, nv, ng, ndim, w_ompc, wbuf1xybcr_ompc,&
            & wbuf2xybcr_ompc, wbuf3r_ompc, wbuf4r_ompc, wbuf5r_ompc, wbuf6r_ompc, ileftx, ilefty, ileftz,&
            & irightx, irighty, irightz, dxdcsinc2_ompc, dydcsinc2_ompc, dxdetanc2_ompc, dydetanc2_ompc,&
            & is_periodic)
            if(is_periodic(1) == 0) call extr_corner_ymin_subroutine(ncoords, nblocks, nx, ny, nz, ng, nv, w_ompc)
          endif
        endif
      endif
    endassociate
  endsubroutine bcswap

  subroutine extr_corner_ymin_subroutine(ncoords,nblocks,nx,ny,nz,ng,nv,w_ompc)
    integer, intent(in) :: nx, ny, nz, ng, nv
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(inout) :: w_ompc
    integer, dimension(3) :: ncoords, nblocks
    integer :: i,j,k,m,iercuda
    if(ncoords(1) == 0) then
      !$omp parallel do default(firstprivate) shared(w_ompc)
      do k = 1,nz
        do j = 1,ng
          do i = 1,ng
            do m=1,nv
              w_ompc(1-i,1-j,k,m) = w_ompc(1-i,1,k,m)
            enddo
          enddo
        enddo
      enddo
    endif
    if(ncoords(1) == nblocks(1)-1) then
      !$omp parallel do default(firstprivate) shared(w_ompc)
      do k = 1,nz
        do j = 1,ng
          do i = 1,ng
            do m=1,nv
              w_ompc(nx+i,1-j,k,m) = w_ompc(nx+i,1,k,m)
            enddo
          enddo
        enddo
      enddo
    endif

  endsubroutine extr_corner_ymin_subroutine



  subroutine bcswap_wake(self, steps)
    class(base_ompc_object), intent(inout) :: self
    logical, dimension(3), optional :: steps
    logical, dimension(3) :: steps_
    integer :: iercuda, indx, indy, indz
    integer, dimension(mpi_status_size) :: istatus

    steps_ = .true. ; if(present(steps)) steps_ = steps

    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv,&
    & w_ompc => self%w_ompc, wbuf4s_ompc => self%wbuf4s_ompc, wbuf3r_ompc => self%wbuf3r_ompc,&
    & wall_tag_ompc => self%wall_tag_ompc, mp_cart => self%field%mp_cart, nrank_wake => self%field%nrank_&
    &wake, iermpi => self%mpi_err, nrank => self%myrank)

      if(steps_(1)) then
        call bcswap_wake_step_1_subroutine(nx, ny, nz, nv, ng, w_ompc, wbuf4s_ompc)
      endif
      if(steps_(2)) then
        indy = nv*nx*ng*nz
        if (nrank_wake==nrank) then
          wbuf3r_ompc = wbuf4s_ompc
        else
          call mpi_sendrecv(wbuf4s_ompc,indy,mpi_prec,nrank_wake,4,wbuf3r_ompc,indy,mpi_prec,nrank_wake,4,mp_cart,istatus,iermpi)
        endif
      endif
      if(steps_(3)) then
        call bcswap_wake_step_3_subroutine(nx, ny, nz, nv, ng, w_ompc, wbuf3r_ompc, wall_tag_ompc)
      endif
    endassociate
  endsubroutine bcswap_wake

  subroutine bcswap_wake_var(self, w_swap, steps)
    class(base_ompc_object), intent(inout) :: self
    real(rkind), dimension(1-self%ng:self%nx+self%ng, 1-self%ng:self%ny+self%ng, 1-self%ng:self%nz+self%ng,1:1) :: w_swap
    logical, dimension(3), optional :: steps
    logical, dimension(3) :: steps_
    integer :: iercuda, indx, indy, indz
    integer, dimension(mpi_status_size) :: istatus

    steps_ = .true. ; if(present(steps)) steps_ = steps

    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv,&
    & w_ompc => self%w_ompc, wbuf4s_ompc => self%wbuf4s_ompc, wbuf3r_ompc => self%wbuf3r_ompc,&
    & wall_tag_ompc => self%wall_tag_ompc, mp_cart => self%field%mp_cart, nrank_wake => self%field%nrank_&
    &wake, iermpi => self%mpi_err, nrank => self%myrank)

      if(steps_(1)) then
        call bcswap_wake_step_1_subroutine(nx, ny, nz, 1, ng, w_swap, wbuf4s_ompc)
      endif
      if(steps_(2)) then
        indy = 1*nx*ng*nz
        if (nrank_wake==nrank) then
          wbuf3r_ompc = wbuf4s_ompc
        else
          call mpi_sendrecv(wbuf4s_ompc,indy,mpi_prec,nrank_wake,4,wbuf3r_ompc,indy,mpi_prec,nrank_wake,4,mp_cart,istatus,iermpi)
        endif
      endif
      if(steps_(3)) then
        call bcswap_wake_step_3_subroutine(nx, ny, nz, 1, ng, w_swap, wbuf3r_ompc, wall_tag_ompc)
      endif
    endassociate
  endsubroutine bcswap_wake_var

  subroutine bcswap_edges_corners(self, steps)
    class(base_ompc_object), intent(inout) :: self
    logical, dimension(3), optional :: steps
    logical, dimension(3) :: steps_
    integer :: i,j,k,m, iercuda, indc
    integer, dimension(mpi_status_size) :: istatus

    steps_ = .true. ; if(present(steps)) steps_ = steps

    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv,&
    & w_ompc => self%w_ompc, wbuf1s_c_ompc => self%wbuf1s_c_ompc, wbuf2s_c_ompc => self%wbuf2s_c_ompc,&
    & wbuf3s_c_ompc => self%wbuf3s_c_ompc, wbuf4s_c_ompc => self%wbuf4s_c_ompc, wbuf1r_c_ompc => self%wbu&
    &f1r_c_ompc, wbuf2r_c_ompc => self%wbuf2r_c_ompc, wbuf3r_c_ompc => self%wbuf3r_c_ompc,&
    & wbuf4r_c_ompc => self%wbuf4r_c_ompc, ileftbottom => self%field%ileftbottom, irightbottom => self%fi&
    &eld%irightbottom, ilefttop => self%field%ilefttop, irighttop => self%field%irighttop,&
    & pbc => self%field%grid%is_xyz_periodic, mp_cart => self%field%mp_cart, iermpi => self%mpi_err,&
    & nrank => self%field%myrank, lxly => self%field%lxly, lxry => self%field%lxry, rxly => self%field%rx&
    &ly, rxry => self%field%rxry,lxlz => self%field%lxlz, lxrz => self%field%lxrz, rxlz => self%field%rxl&
    &z, rxrz => self%field%rxrz,lylz => self%field%lylz, lyrz => self%field%lyrz, rylz => self%field%rylz&
    &, ryrz => self%field%ryrz,lxlylz => self%field%lxlylz, lxlyrz => self%field%lxlyrz,&
    & lxrylz => self%field%lxrylz, rxlylz => self%field%rxlylz,lxryrz => self%field%lxryrz,&
    & rxlyrz => self%field%rxlyrz, rxrylz => self%field%rxrylz, rxryrz => self%field%rxryrz,&
    &wbuf_lxly_s_ompc => self%wbuf_lxly_s_ompc , wbuf_lxry_s_ompc => self%wbuf_lxry_s_ompc,&
    & wbuf_rxly_s_ompc => self%wbuf_rxly_s_ompc , wbuf_rxry_s_ompc => self%wbuf_rxry_s_ompc,&
    & wbuf_lxlz_s_ompc => self%wbuf_lxlz_s_ompc , wbuf_lxrz_s_ompc => self%wbuf_lxrz_s_ompc,&
    & wbuf_rxlz_s_ompc => self%wbuf_rxlz_s_ompc , wbuf_rxrz_s_ompc => self%wbuf_rxrz_s_ompc,&
    & wbuf_lylz_s_ompc => self%wbuf_lylz_s_ompc , wbuf_lyrz_s_ompc => self%wbuf_lyrz_s_ompc,&
    & wbuf_rylz_s_ompc => self%wbuf_rylz_s_ompc , wbuf_ryrz_s_ompc => self%wbuf_ryrz_s_ompc,&
    & wbuf_lxlylz_s_ompc => self%wbuf_lxlylz_s_ompc , wbuf_lxlyrz_s_ompc => self%wbuf_lxlyrz_s_ompc,&
    & wbuf_lxrylz_s_ompc => self%wbuf_lxrylz_s_ompc , wbuf_lxryrz_s_ompc => self%wbuf_lxryrz_s_ompc,&
    & wbuf_rxlylz_s_ompc => self%wbuf_rxlylz_s_ompc , wbuf_rxlyrz_s_ompc => self%wbuf_rxlyrz_s_ompc,&
    & wbuf_rxrylz_s_ompc => self%wbuf_rxrylz_s_ompc , wbuf_rxryrz_s_ompc => self%wbuf_rxryrz_s_ompc,&
    & wbuf_lxly_r_ompc => self%wbuf_lxly_r_ompc , wbuf_lxry_r_ompc => self%wbuf_lxry_r_ompc,&
    & wbuf_rxly_r_ompc => self%wbuf_rxly_r_ompc , wbuf_rxry_r_ompc => self%wbuf_rxry_r_ompc,&
    & wbuf_lxlz_r_ompc => self%wbuf_lxlz_r_ompc , wbuf_lxrz_r_ompc => self%wbuf_lxrz_r_ompc,&
    & wbuf_rxlz_r_ompc => self%wbuf_rxlz_r_ompc , wbuf_rxrz_r_ompc => self%wbuf_rxrz_r_ompc,&
    & wbuf_lylz_r_ompc => self%wbuf_lylz_r_ompc , wbuf_lyrz_r_ompc => self%wbuf_lyrz_r_ompc,&
    & wbuf_rylz_r_ompc => self%wbuf_rylz_r_ompc , wbuf_ryrz_r_ompc => self%wbuf_ryrz_r_ompc,&
    & wbuf_lxlylz_r_ompc => self%wbuf_lxlylz_r_ompc , wbuf_lxlyrz_r_ompc => self%wbuf_lxlyrz_r_ompc,&
    & wbuf_lxrylz_r_ompc => self%wbuf_lxrylz_r_ompc , wbuf_lxryrz_r_ompc => self%wbuf_lxryrz_r_ompc,&
    & wbuf_rxlylz_r_ompc => self%wbuf_rxlylz_r_ompc , wbuf_rxlyrz_r_ompc => self%wbuf_rxlyrz_r_ompc,&
    & wbuf_rxrylz_r_ompc => self%wbuf_rxrylz_r_ompc , wbuf_rxryrz_r_ompc => self%wbuf_rxryrz_r_ompc)

      if(steps_(1)) then
        call bcswap_corner_step_1_subroutine(nx, ny, nz, nv, ng, w_ompc, wbuf1s_c_ompc,&
        & wbuf2s_c_ompc, wbuf3s_c_ompc, wbuf4s_c_ompc)
        if(pbc(2)) then
          call bcswap_corner_step_1b_subroutine(nx, ny, nz, nv, ng, w_ompc, wbuf_lxly_s_ompc,&
          & wbuf_lxry_s_ompc, wbuf_rxly_s_ompc, wbuf_rxry_s_ompc, wbuf_lylz_s_ompc, wbuf_lyrz_s_ompc,&
          & wbuf_rylz_s_ompc, wbuf_ryrz_s_ompc, wbuf_lxlylz_s_ompc, wbuf_lxlyrz_s_ompc, wbuf_lxrylz_s_ompc,&
          & wbuf_lxryrz_s_ompc, wbuf_rxlylz_s_ompc, wbuf_rxlyrz_s_ompc, wbuf_rxrylz_s_ompc, wbuf_rxryrz_s_ompc)
        endif
      endif
      if(steps_(2)) then
        indc = nv*ny*ng*ng
        if (ileftbottom == nrank) then
          wbuf2r_c_ompc = wbuf1s_c_ompc
          wbuf1r_c_ompc = wbuf2s_c_ompc
        else
          call mpi_sendrecv(wbuf1s_c_ompc,indc,mpi_prec,ileftbottom ,1, wbuf2r_c_ompc,indc,mpi_prec,&
          &irighttop ,1,mp_cart,istatus,iermpi)
          call mpi_sendrecv(wbuf2s_c_ompc,indc,mpi_prec,irighttop ,2, wbuf1r_c_ompc,indc,mpi_prec,&
          &ileftbottom ,2,mp_cart,istatus,iermpi)
        endif
        if (ilefttop == nrank) then
          wbuf4r_c_ompc = wbuf3s_c_ompc
          wbuf3r_c_ompc = wbuf4s_c_ompc
        else
          call mpi_sendrecv(wbuf3s_c_ompc,indc,mpi_prec,ilefttop ,3, wbuf4r_c_ompc,indc,mpi_prec,&
          &irightbottom,3,mp_cart,istatus,iermpi)
          call mpi_sendrecv(wbuf4s_c_ompc,indc,mpi_prec,irightbottom,4, wbuf3r_c_ompc,indc,mpi_prec,&
          &ilefttop ,4,mp_cart,istatus,iermpi)
        endif
        if(pbc(2)) then
          indc = nv*ng*ng*nz
          if (lxly == nrank) then
            wbuf_rxry_r_ompc = wbuf_lxly_s_ompc
            wbuf_lxly_r_ompc = wbuf_rxry_s_ompc
          else
            call mpi_sendrecv(wbuf_lxly_s_ompc,indc,mpi_prec,lxly ,5, wbuf_rxry_r_ompc,indc,&
            &mpi_prec,rxry ,5,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rxry_s_ompc,indc,mpi_prec,rxry ,6, wbuf_lxly_r_ompc,indc,&
            &mpi_prec,lxly ,6,mp_cart,istatus,iermpi)
          endif
          if (lxry == nrank) then
            wbuf_rxly_r_ompc = wbuf_lxry_s_ompc
            wbuf_lxry_r_ompc = wbuf_rxly_s_ompc
          else
            call mpi_sendrecv(wbuf_lxry_s_ompc,indc,mpi_prec,lxry ,7, wbuf_rxly_r_ompc,indc,&
            &mpi_prec,rxly ,7,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rxly_s_ompc,indc,mpi_prec,rxly ,8, wbuf_lxry_r_ompc,indc,&
            &mpi_prec,lxry ,8,mp_cart,istatus,iermpi)
          endif
          indc = nv*nx*ng*ng
          if (lylz == nrank) then
            wbuf_ryrz_r_ompc = wbuf_lylz_s_ompc
            wbuf_lylz_r_ompc = wbuf_ryrz_s_ompc
          else
            call mpi_sendrecv(wbuf_lylz_s_ompc,indc,mpi_prec,lylz ,9, wbuf_ryrz_r_ompc,indc,&
            &mpi_prec,ryrz ,9,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_ryrz_s_ompc,indc,mpi_prec,ryrz ,10, wbuf_lylz_r_ompc,indc,&
            &mpi_prec,lylz ,10,mp_cart,istatus,iermpi)
          endif
          if (lyrz == nrank) then
            wbuf_rylz_r_ompc = wbuf_lyrz_s_ompc
            wbuf_lyrz_r_ompc = wbuf_rylz_s_ompc
          else
            call mpi_sendrecv(wbuf_lyrz_s_ompc,indc,mpi_prec,lyrz ,11, wbuf_rylz_r_ompc,indc,&
            &mpi_prec,rylz ,11,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rylz_s_ompc,indc,mpi_prec,rylz ,12, wbuf_lyrz_r_ompc,indc,&
            &mpi_prec,lyrz ,12,mp_cart,istatus,iermpi)
          endif
          indc = nv*ng*ng*ng
          if (lxlylz == nrank) then
            wbuf_rxryrz_r_ompc = wbuf_lxlylz_s_ompc
            wbuf_rxrylz_r_ompc = wbuf_lxlyrz_s_ompc
            wbuf_rxlyrz_r_ompc = wbuf_lxrylz_s_ompc
            wbuf_rxlylz_r_ompc = wbuf_lxryrz_s_ompc
            wbuf_lxryrz_r_ompc = wbuf_rxlylz_s_ompc
            wbuf_lxrylz_r_ompc = wbuf_rxlyrz_s_ompc
            wbuf_lxlyrz_r_ompc = wbuf_rxrylz_s_ompc
            wbuf_lxlylz_r_ompc = wbuf_rxryrz_s_ompc
          else
            call mpi_sendrecv(wbuf_lxlylz_s_ompc,indc,mpi_prec,lxlylz ,13, wbuf_rxryrz_r_ompc,indc,&
            &mpi_prec,rxryrz ,13,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rxryrz_s_ompc,indc,mpi_prec,rxryrz ,14, wbuf_lxlylz_r_ompc,indc,&
            &mpi_prec,lxlylz ,14,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_lxlyrz_s_ompc,indc,mpi_prec,lxlyrz ,15, wbuf_rxrylz_r_ompc,indc,&
            &mpi_prec,rxrylz ,15,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rxrylz_s_ompc,indc,mpi_prec,rxrylz ,16, wbuf_lxlyrz_r_ompc,indc,&
            &mpi_prec,lxlyrz ,16,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rxlylz_s_ompc,indc,mpi_prec,rxlylz ,17, wbuf_lxryrz_r_ompc,indc,&
            &mpi_prec,lxryrz ,17,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_lxryrz_s_ompc,indc,mpi_prec,lxryrz ,18, wbuf_rxlylz_r_ompc,indc,&
            &mpi_prec,rxlylz ,18,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rxlyrz_s_ompc,indc,mpi_prec,rxlyrz ,19, wbuf_lxrylz_r_ompc,indc,&
            &mpi_prec,lxrylz ,19,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_lxrylz_s_ompc,indc,mpi_prec,lxrylz ,20, wbuf_rxlyrz_r_ompc,indc,&
            &mpi_prec,rxlyrz ,20,mp_cart,istatus,iermpi)
          endif
        endif
      endif
      if(steps_(3)) then
        call bcswap_corner_step_3_subroutine(nx, ny, nz, nv, ng, w_ompc, wbuf1r_c_ompc,&
        & wbuf2r_c_ompc, wbuf3r_c_ompc, wbuf4r_c_ompc, ileftbottom, ilefttop, irightbottom, irighttop)
        if(pbc(2)) then
          call bcswap_corner_step_3b_subroutine(nx, ny, nz, nv, ng, w_ompc, wbuf_lxly_r_ompc,&
          & wbuf_lxry_r_ompc, wbuf_rxly_r_ompc, wbuf_rxry_r_ompc, wbuf_lylz_r_ompc, wbuf_lyrz_r_ompc,&
          & wbuf_rylz_r_ompc, wbuf_ryrz_r_ompc, wbuf_lxlylz_r_ompc, wbuf_lxlyrz_r_ompc, wbuf_lxrylz_r_ompc,&
          & wbuf_lxryrz_r_ompc, wbuf_rxlylz_r_ompc, wbuf_rxlyrz_r_ompc, wbuf_rxrylz_r_ompc, wbuf_rxryrz_r_ompc,&
          & lxly, lxry, rxly, rxry, lylz, lyrz, rylz, ryrz, lxlylz, lxlyrz, lxrylz, lxryrz, rxlylz, rxlyrz,&
          & rxrylz, rxryrz)
        endif
      endif
    endassociate
  endsubroutine bcswap_edges_corners

  subroutine bcswap_edges_corners_var(self, w_swap, steps)
    class(base_ompc_object), intent(inout) :: self
    real(rkind), dimension(1-self%ng:self%nx+self%ng, 1-self%ng:self%ny+self%ng, 1-self%ng:self%nz+self%ng,1:1) :: w_swap
    logical, dimension(3), optional :: steps
    logical, dimension(3) :: steps_
    integer :: i,j,k,m, iercuda, indc
    integer, dimension(mpi_status_size) :: istatus

    steps_ = .true. ; if(present(steps)) steps_ = steps

    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv,&
    & wbuf1s_c_ompc => self%wbuf1s_c_ompc, wbuf2s_c_ompc => self%wbuf2s_c_ompc, wbuf3s_c_ompc => self%wbu&
    &f3s_c_ompc, wbuf4s_c_ompc => self%wbuf4s_c_ompc, wbuf1r_c_ompc => self%wbuf1r_c_ompc,&
    & wbuf2r_c_ompc => self%wbuf2r_c_ompc, wbuf3r_c_ompc => self%wbuf3r_c_ompc, wbuf4r_c_ompc => self%wbu&
    &f4r_c_ompc, ileftbottom => self%field%ileftbottom, irightbottom => self%field%irightbottom,&
    & ilefttop => self%field%ilefttop, irighttop => self%field%irighttop, pbc => self%field%grid%is_xyz_p&
    &eriodic, mp_cart => self%field%mp_cart, iermpi => self%mpi_err, nrank => self%field%myrank,&
    & lxly => self%field%lxly, lxry => self%field%lxry, rxly => self%field%rxly, rxry => self%field%rxry,&
    &lxlz => self%field%lxlz, lxrz => self%field%lxrz, rxlz => self%field%rxlz, rxrz => self%field%rxrz,&
    &lylz => self%field%lylz, lyrz => self%field%lyrz, rylz => self%field%rylz, ryrz => self%field%ryrz,&
    &lxlylz => self%field%lxlylz, lxlyrz => self%field%lxlyrz, lxrylz => self%field%lxrylz,&
    & rxlylz => self%field%rxlylz,lxryrz => self%field%lxryrz, rxlyrz => self%field%rxlyrz,&
    & rxrylz => self%field%rxrylz, rxryrz => self%field%rxryrz,wbuf_lxly_s_ompc => self%wbuf_lxly_s_ompc &
    &, wbuf_lxry_s_ompc => self%wbuf_lxry_s_ompc, wbuf_rxly_s_ompc => self%wbuf_rxly_s_ompc ,&
    & wbuf_rxry_s_ompc => self%wbuf_rxry_s_ompc, wbuf_lxlz_s_ompc => self%wbuf_lxlz_s_ompc ,&
    & wbuf_lxrz_s_ompc => self%wbuf_lxrz_s_ompc, wbuf_rxlz_s_ompc => self%wbuf_rxlz_s_ompc ,&
    & wbuf_rxrz_s_ompc => self%wbuf_rxrz_s_ompc, wbuf_lylz_s_ompc => self%wbuf_lylz_s_ompc ,&
    & wbuf_lyrz_s_ompc => self%wbuf_lyrz_s_ompc, wbuf_rylz_s_ompc => self%wbuf_rylz_s_ompc ,&
    & wbuf_ryrz_s_ompc => self%wbuf_ryrz_s_ompc, wbuf_lxlylz_s_ompc => self%wbuf_lxlylz_s_ompc ,&
    & wbuf_lxlyrz_s_ompc => self%wbuf_lxlyrz_s_ompc, wbuf_lxrylz_s_ompc => self%wbuf_lxrylz_s_ompc ,&
    & wbuf_lxryrz_s_ompc => self%wbuf_lxryrz_s_ompc, wbuf_rxlylz_s_ompc => self%wbuf_rxlylz_s_ompc ,&
    & wbuf_rxlyrz_s_ompc => self%wbuf_rxlyrz_s_ompc, wbuf_rxrylz_s_ompc => self%wbuf_rxrylz_s_ompc ,&
    & wbuf_rxryrz_s_ompc => self%wbuf_rxryrz_s_ompc, wbuf_lxly_r_ompc => self%wbuf_lxly_r_ompc ,&
    & wbuf_lxry_r_ompc => self%wbuf_lxry_r_ompc, wbuf_rxly_r_ompc => self%wbuf_rxly_r_ompc ,&
    & wbuf_rxry_r_ompc => self%wbuf_rxry_r_ompc, wbuf_lxlz_r_ompc => self%wbuf_lxlz_r_ompc ,&
    & wbuf_lxrz_r_ompc => self%wbuf_lxrz_r_ompc, wbuf_rxlz_r_ompc => self%wbuf_rxlz_r_ompc ,&
    & wbuf_rxrz_r_ompc => self%wbuf_rxrz_r_ompc, wbuf_lylz_r_ompc => self%wbuf_lylz_r_ompc ,&
    & wbuf_lyrz_r_ompc => self%wbuf_lyrz_r_ompc, wbuf_rylz_r_ompc => self%wbuf_rylz_r_ompc ,&
    & wbuf_ryrz_r_ompc => self%wbuf_ryrz_r_ompc, wbuf_lxlylz_r_ompc => self%wbuf_lxlylz_r_ompc ,&
    & wbuf_lxlyrz_r_ompc => self%wbuf_lxlyrz_r_ompc, wbuf_lxrylz_r_ompc => self%wbuf_lxrylz_r_ompc ,&
    & wbuf_lxryrz_r_ompc => self%wbuf_lxryrz_r_ompc, wbuf_rxlylz_r_ompc => self%wbuf_rxlylz_r_ompc ,&
    & wbuf_rxlyrz_r_ompc => self%wbuf_rxlyrz_r_ompc, wbuf_rxrylz_r_ompc => self%wbuf_rxrylz_r_ompc ,&
    & wbuf_rxryrz_r_ompc => self%wbuf_rxryrz_r_ompc)

      if(steps_(1)) then
        call bcswap_corner_step_1_subroutine(nx, ny, nz, 1, ng, w_swap, wbuf1s_c_ompc, wbuf2s_c_ompc, wbuf3s_c_ompc, wbuf4s_c_ompc)
        if(pbc(2)) then
          call bcswap_corner_step_1b_subroutine(nx, ny, nz, 1, ng, w_swap, wbuf_lxly_s_ompc,&
          & wbuf_lxry_s_ompc, wbuf_rxly_s_ompc, wbuf_rxry_s_ompc, wbuf_lylz_s_ompc, wbuf_lyrz_s_ompc,&
          & wbuf_rylz_s_ompc, wbuf_ryrz_s_ompc, wbuf_lxlylz_s_ompc, wbuf_lxlyrz_s_ompc, wbuf_lxrylz_s_ompc,&
          & wbuf_lxryrz_s_ompc, wbuf_rxlylz_s_ompc, wbuf_rxlyrz_s_ompc, wbuf_rxrylz_s_ompc, wbuf_rxryrz_s_ompc)
        endif
      endif
      if(steps_(2)) then
        indc = 1*ny*ng*ng
        if (ileftbottom == nrank) then
          wbuf2r_c_ompc = wbuf1s_c_ompc
          wbuf1r_c_ompc = wbuf2s_c_ompc
        else
          call mpi_sendrecv(wbuf1s_c_ompc,indc,mpi_prec,ileftbottom ,1, wbuf2r_c_ompc,indc,mpi_prec,&
          &irighttop ,1,mp_cart,istatus,iermpi)
          call mpi_sendrecv(wbuf2s_c_ompc,indc,mpi_prec,irighttop ,2, wbuf1r_c_ompc,indc,mpi_prec,&
          &ileftbottom ,2,mp_cart,istatus,iermpi)
        endif
        if (ilefttop == nrank) then
          wbuf4r_c_ompc = wbuf3s_c_ompc
          wbuf3r_c_ompc = wbuf4s_c_ompc
        else
          call mpi_sendrecv(wbuf3s_c_ompc,indc,mpi_prec,ilefttop ,3, wbuf4r_c_ompc,indc,mpi_prec,&
          &irightbottom,3,mp_cart,istatus,iermpi)
          call mpi_sendrecv(wbuf4s_c_ompc,indc,mpi_prec,irightbottom,4, wbuf3r_c_ompc,indc,mpi_prec,&
          &ilefttop ,4,mp_cart,istatus,iermpi)
        endif
        if(pbc(2)) then
          indc = 1*ng*ng*nz
          if (lxly == nrank) then
            wbuf_rxry_r_ompc = wbuf_lxly_s_ompc
            wbuf_lxly_r_ompc = wbuf_rxry_s_ompc
          else
            call mpi_sendrecv(wbuf_lxly_s_ompc,indc,mpi_prec,lxly ,5, wbuf_rxry_r_ompc,indc,&
            &mpi_prec,rxry ,5,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rxry_s_ompc,indc,mpi_prec,rxry ,6, wbuf_lxly_r_ompc,indc,&
            &mpi_prec,lxly ,6,mp_cart,istatus,iermpi)
          endif
          if (lxry == nrank) then
            wbuf_rxly_r_ompc = wbuf_lxry_s_ompc
            wbuf_lxry_r_ompc = wbuf_rxly_s_ompc
          else
            call mpi_sendrecv(wbuf_lxry_s_ompc,indc,mpi_prec,lxry ,7, wbuf_rxly_r_ompc,indc,&
            &mpi_prec,rxly ,7,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rxly_s_ompc,indc,mpi_prec,rxly ,8, wbuf_lxry_r_ompc,indc,&
            &mpi_prec,lxry ,8,mp_cart,istatus,iermpi)
          endif
          indc = 1*nx*ng*ng
          if (lylz == nrank) then
            wbuf_ryrz_r_ompc = wbuf_lylz_s_ompc
            wbuf_lylz_r_ompc = wbuf_ryrz_s_ompc
          else
            call mpi_sendrecv(wbuf_lylz_s_ompc,indc,mpi_prec,lylz ,9, wbuf_ryrz_r_ompc,indc,&
            &mpi_prec,ryrz ,9,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_ryrz_s_ompc,indc,mpi_prec,ryrz ,10, wbuf_lylz_r_ompc,indc,&
            &mpi_prec,lylz ,10,mp_cart,istatus,iermpi)
          endif
          if (lyrz == nrank) then
            wbuf_rylz_r_ompc = wbuf_lyrz_s_ompc
            wbuf_lyrz_r_ompc = wbuf_rylz_s_ompc
          else
            call mpi_sendrecv(wbuf_lyrz_s_ompc,indc,mpi_prec,lyrz ,11, wbuf_rylz_r_ompc,indc,&
            &mpi_prec,rylz ,11,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rylz_s_ompc,indc,mpi_prec,rylz ,12, wbuf_lyrz_r_ompc,indc,&
            &mpi_prec,lyrz ,12,mp_cart,istatus,iermpi)
          endif
          indc = 1*ng*ng*ng
          if (lxlylz == nrank) then
            wbuf_rxryrz_r_ompc = wbuf_lxlylz_s_ompc
            wbuf_rxrylz_r_ompc = wbuf_lxlyrz_s_ompc
            wbuf_rxlyrz_r_ompc = wbuf_lxrylz_s_ompc
            wbuf_rxlylz_r_ompc = wbuf_lxryrz_s_ompc
            wbuf_lxryrz_r_ompc = wbuf_rxlylz_s_ompc
            wbuf_lxrylz_r_ompc = wbuf_rxlyrz_s_ompc
            wbuf_lxlyrz_r_ompc = wbuf_rxrylz_s_ompc
            wbuf_lxlylz_r_ompc = wbuf_rxryrz_s_ompc
          else
            call mpi_sendrecv(wbuf_lxlylz_s_ompc,indc,mpi_prec,lxlylz ,13, wbuf_rxryrz_r_ompc,indc,&
            &mpi_prec,rxryrz ,13,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rxryrz_s_ompc,indc,mpi_prec,rxryrz ,14, wbuf_lxlylz_r_ompc,indc,&
            &mpi_prec,lxlylz ,14,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_lxlyrz_s_ompc,indc,mpi_prec,lxlyrz ,15, wbuf_rxrylz_r_ompc,indc,&
            &mpi_prec,rxrylz ,15,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rxrylz_s_ompc,indc,mpi_prec,rxrylz ,16, wbuf_lxlyrz_r_ompc,indc,&
            &mpi_prec,lxlyrz ,16,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rxlylz_s_ompc,indc,mpi_prec,rxlylz ,17, wbuf_lxryrz_r_ompc,indc,&
            &mpi_prec,lxryrz ,17,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_lxryrz_s_ompc,indc,mpi_prec,lxryrz ,18, wbuf_rxlylz_r_ompc,indc,&
            &mpi_prec,rxlylz ,18,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rxlyrz_s_ompc,indc,mpi_prec,rxlyrz ,19, wbuf_lxrylz_r_ompc,indc,&
            &mpi_prec,lxrylz ,19,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_lxrylz_s_ompc,indc,mpi_prec,lxrylz ,20, wbuf_rxlyrz_r_ompc,indc,&
            &mpi_prec,rxlyrz ,20,mp_cart,istatus,iermpi)
          endif
        endif
      endif
      if(steps_(3)) then
        call bcswap_corner_step_3_subroutine(nx, ny, nz, 1, ng, w_swap, wbuf1r_c_ompc,&
        & wbuf2r_c_ompc, wbuf3r_c_ompc, wbuf4r_c_ompc, ileftbottom, ilefttop, irightbottom, irighttop)
        if(pbc(2)) then
          call bcswap_corner_step_3b_subroutine(nx, ny, nz, 1, ng, w_swap, wbuf_lxly_r_ompc,&
          & wbuf_lxry_r_ompc, wbuf_rxly_r_ompc, wbuf_rxry_r_ompc, wbuf_lylz_r_ompc, wbuf_lyrz_r_ompc,&
          & wbuf_rylz_r_ompc, wbuf_ryrz_r_ompc, wbuf_lxlylz_r_ompc, wbuf_lxlyrz_r_ompc, wbuf_lxrylz_r_ompc,&
          & wbuf_lxryrz_r_ompc, wbuf_rxlylz_r_ompc, wbuf_rxlyrz_r_ompc, wbuf_rxrylz_r_ompc, wbuf_rxryrz_r_ompc,&
          & lxly, lxry, rxly, rxry, lylz, lyrz, rylz, ryrz, lxlylz, lxlyrz, lxrylz, lxryrz, rxlylz, rxlyrz,&
          & rxrylz, rxryrz)
        endif
      endif
    endassociate
  endsubroutine bcswap_edges_corners_var


  subroutine bcswap_var(self, w_swap, steps)
    class(base_ompc_object), intent(inout) :: self
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

    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, wbuf1s_ompc => self%wbuf1s&
    &_ompc, wbuf2s_ompc => self%wbuf2s_ompc, wbuf3s_ompc => self%wbuf3s_ompc, wbuf4s_ompc => self%wbuf4s_&
    &ompc, wbuf5s_ompc => self%wbuf5s_ompc, wbuf6s_ompc => self%wbuf6s_ompc, wbuf1r_ompc => self%wbuf1r_o&
    &mpc, wbuf2r_ompc => self%wbuf2r_ompc, wbuf3r_ompc => self%wbuf3r_ompc, wbuf4r_ompc => self%wbuf4r_om&
    &pc, wbuf5r_ompc => self%wbuf5r_ompc, wbuf6r_ompc => self%wbuf6r_ompc, ileftx => self%field%ileftx,&
    &irightx => self%field%irightx, nrank_x => self%field%nrank_x, ilefty => self%field%ilefty,&
    &irighty => self%field%irighty, nrank_y => self%field%nrank_y, ileftz => self%field%ileftz,&
    &irightz => self%field%irightz, nrank_z => self%field%nrank_z, dcsidxnc2_ompc => self%dcsidxnc2_ompc,&
    & dcsidync2_ompc => self%dcsidync2_ompc, detadxnc2_ompc => self%detadxnc2_ompc, detadync2_ompc => sel&
    &f%detadync2_ompc, dxdcsinc2_ompc => self%dxdcsinc2_ompc, dydcsinc2_ompc => self%dydcsinc2_ompc,&
    & dxdetanc2_ompc => self%dxdetanc2_ompc, dydetanc2_ompc => self%dydetanc2_ompc, mp_cartx => self%fiel&
    &d%mp_cartx, mp_carty => self%field%mp_carty, mp_cartz => self%field%mp_cartz, ndim => self%field%gri&
    &d%ndim)

      if(steps_(1)) then
        if(self%field%grid%grid_dim == 1) then
          call bcswap_step_1_subroutine(nx, ny, nz, 1, ng, ndim, w_swap, wbuf1s_ompc, wbuf2s_ompc,&
          & wbuf3s_ompc, wbuf4s_ompc, wbuf5s_ompc, wbuf6s_ompc)
        elseif(self%field%grid%grid_dim == 2) then
          call bcswap_step_1_subroutine(nx, ny, nz, 1, ng, ndim, w_swap, wbuf1s_ompc, wbuf2s_ompc,&
          & wbuf3s_ompc, wbuf4s_ompc, wbuf5s_ompc, wbuf6s_ompc)
        endif
      endif
      if(steps_(2)) then
        indx = ng*ny*nz
        indy = nx*ng*nz
        indz = nx*ny*ng
        if(ileftx == nrank_x) then
          wbuf2r_ompc = wbuf1s_ompc
          wbuf1r_ompc = wbuf2s_ompc
        else
          call mpi_sendrecv(wbuf1s_ompc,indx,mpi_prec,ileftx ,1,wbuf2r_ompc,indx,mpi_prec,irightx,1,mp_cartx,istatus,self%mpi_err)
          call mpi_sendrecv(wbuf2s_ompc,indx,mpi_prec,irightx,2,wbuf1r_ompc,indx,mpi_prec,ileftx ,2,mp_cartx,istatus,self%mpi_err)
        endif
        if(ilefty == nrank_y) then
          wbuf4r_ompc = wbuf3s_ompc
          wbuf3r_ompc = wbuf4s_ompc
        else
          call mpi_sendrecv(wbuf3s_ompc,indy,mpi_prec,ilefty ,3,wbuf4r_ompc,indy,mpi_prec,irighty,3,mp_carty,istatus,self%mpi_err)
          call mpi_sendrecv(wbuf4s_ompc,indy,mpi_prec,irighty,4,wbuf3r_ompc,indy,mpi_prec,ilefty ,4,mp_carty,istatus,self%mpi_err)
        endif
        if(ndim == 3) then
          if(ileftz == nrank_z) then
            wbuf6r_ompc = wbuf5s_ompc
            wbuf5r_ompc = wbuf6s_ompc
          else
            call mpi_sendrecv(wbuf5s_ompc,indz,mpi_prec,ileftz ,5,wbuf6r_ompc,indz,mpi_prec,irightz,&
            &5,mp_cartz,istatus,self%mpi_err)
            call mpi_sendrecv(wbuf6s_ompc,indz,mpi_prec,irightz,6,wbuf5r_ompc,indz,mpi_prec,ileftz ,&
            &6,mp_cartz,istatus,self%mpi_err)
          endif
        endif
      endif
      if(steps_(3)) then
        if(self%field%grid%grid_dim == 1) then
          call bcswap_step_3_subroutine(nx, ny, nz, 1, ng, ndim, w_swap, wbuf1r_ompc, wbuf2r_ompc,&
          & wbuf3r_ompc, wbuf4r_ompc, wbuf5r_ompc, wbuf6r_ompc, ileftx, ilefty, ileftz, irightx, irighty,&
          & irightz)
        elseif(self%field%grid%grid_dim == 2) then
          call bcswap_step_3_subroutine(nx, ny, nz, 1, ng, ndim, w_swap, wbuf1r_ompc, wbuf2r_ompc,&
          & wbuf3r_ompc, wbuf4r_ompc, wbuf5r_ompc, wbuf6r_ompc, ileftx, ilefty, ileftz, irightx, irighty,&
          & irightz)
        endif
      endif
    endassociate
  endsubroutine bcswap_var

  subroutine copy_from_field(self)
    class(base_ompc_object), intent(inout) :: self
    integer :: i,j,k,iv

    self%w_ompc = self%field%w

  endsubroutine copy_from_field


  subroutine copy_to_field(self)
    class(base_ompc_object), intent(inout) :: self
    integer :: i,j,k,iv

    self%field%w = self%w_ompc

  endsubroutine copy_to_field

  subroutine initialize(self, field, gpu_bind)
    !< Initialize base backend.
    class(base_ompc_object), intent(inout) :: self !< The base backend.
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


    call self%alloc()

  endsubroutine initialize


endmodule streams_base_ompc_object


