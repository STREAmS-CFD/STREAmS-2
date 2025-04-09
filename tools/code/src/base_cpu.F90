module streams_base_cpu_object

  use streams_field_object, only : field_object
  use streams_parameters
  use mpi


  implicit none
  private
  public :: base_cpu_object

  type :: base_cpu_object
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
    real(rkind), pointer, contiguous, dimension(:,:,:,:) :: w_cpu

    real(rkind), dimension(:,:,:,:), allocatable :: wbuf1s_cpu, wbuf2s_cpu, wbuf3s_cpu, wbuf4s_cpu, wbuf5s_cpu, wbuf6s_cpu
    real(rkind), dimension(:,:,:,:), allocatable :: wbuf1r_cpu, wbuf2r_cpu, wbuf3r_cpu, wbuf4r_cpu, wbuf5r_cpu, wbuf6r_cpu
    real(rkind), dimension(:,:,:,:), allocatable :: wbuf1xybcs_cpu, wbuf2xybcs_cpu, wbuf1xybcr_cpu, wbuf2xybcr_cpu
    real(rkind), dimension(:,:,:,:), allocatable :: wbuf1s_c_cpu, wbuf2s_c_cpu, wbuf3s_c_cpu, wbuf4s_c_cpu
    real(rkind), dimension(:,:,:,:), allocatable :: wbuf1r_c_cpu, wbuf2r_c_cpu, wbuf3r_c_cpu, wbuf4r_c_cpu
    real(rkind), allocatable, dimension(:) :: dcsidx_cpu, dcsidx2_cpu, dcsidxs_cpu
    real(rkind), allocatable, dimension(:) :: detady_cpu, detady2_cpu, detadys_cpu
    real(rkind), allocatable, dimension(:) :: dzitdz_cpu, dzitdz2_cpu, dzitdzs_cpu
    real(rkind), allocatable, dimension(:) :: x_cpu, y_cpu, z_cpu, yn_cpu
    real(rkind), allocatable, dimension(:,:) :: xc2_cpu,yc2_cpu
    real(rkind), allocatable, dimension(:,:) :: dcsidxc2_cpu,dcsidyc2_cpu
    real(rkind), allocatable, dimension(:,:) :: detadxc2_cpu,detadyc2_cpu
    real(rkind), allocatable, dimension(:,:) :: dcsidxnc2_cpu,dcsidync2_cpu
    real(rkind), allocatable, dimension(:,:) :: detadxnc2_cpu,detadync2_cpu
    real(rkind), allocatable, dimension(:,:) :: dxdcsic2_cpu,dydcsic2_cpu
    real(rkind), allocatable, dimension(:,:) :: dxdetac2_cpu,dydetac2_cpu
    real(rkind), allocatable, dimension(:,:) :: dxdcsinc2_cpu,dydcsinc2_cpu
    real(rkind), allocatable, dimension(:,:) :: dxdetanc2_cpu,dydetanc2_cpu
    real(rkind), allocatable, dimension(:,:) :: jac_cpu,mcsijac1_cpu,metajac1_cpu
    real(rkind), allocatable, dimension(:,:) :: mcsi_cpu,meta_cpu,csimod_cpu,etamod_cpu
    real(rkind), allocatable, dimension(:,:) :: g1_cpu,g2_cpu,g12_cpu
    real(rkind), allocatable, dimension(:,:) :: wbuftus_cpu, wbuftes_cpu , wbuftur_cpu, wbufter_cpu
    real(rkind), allocatable, dimension(:,:) :: theta_ij_cpu

    integer, allocatable, dimension(:) :: wall_tag_cpu

    real(rkind), allocatable,dimension(:,:,:,:) :: wbuf_lxly_s_cpu , wbuf_lxry_s_cpu , wbuf_rxly_s_cpu , wbuf_rxry_s_cpu
    real(rkind), allocatable,dimension(:,:,:,:) :: wbuf_lxlz_s_cpu , wbuf_lxrz_s_cpu , wbuf_rxlz_s_cpu , wbuf_rxrz_s_cpu
    real(rkind), allocatable,dimension(:,:,:,:) :: wbuf_lylz_s_cpu , wbuf_lyrz_s_cpu , wbuf_rylz_s_cpu , wbuf_ryrz_s_cpu
    real(rkind), allocatable,dimension(:,:,:,:) :: wbuf_lxlylz_s_cpu,wbuf_lxlyrz_s_cpu,wbuf_lxrylz_s_cpu,wbuf_lxryrz_s_cpu
    real(rkind), allocatable,dimension(:,:,:,:) :: wbuf_rxlylz_s_cpu,wbuf_rxlyrz_s_cpu,wbuf_rxrylz_s_cpu,wbuf_rxryrz_s_cpu
    real(rkind), allocatable,dimension(:,:,:,:) :: wbuf_lxly_r_cpu , wbuf_lxry_r_cpu , wbuf_rxly_r_cpu , wbuf_rxry_r_cpu
    real(rkind), allocatable,dimension(:,:,:,:) :: wbuf_lxlz_r_cpu , wbuf_lxrz_r_cpu , wbuf_rxlz_r_cpu , wbuf_rxrz_r_cpu
    real(rkind), allocatable,dimension(:,:,:,:) :: wbuf_lylz_r_cpu , wbuf_lyrz_r_cpu , wbuf_rylz_r_cpu , wbuf_ryrz_r_cpu
    real(rkind), allocatable,dimension(:,:,:,:) :: wbuf_lxlylz_r_cpu,wbuf_lxlyrz_r_cpu,wbuf_lxrylz_r_cpu,wbuf_lxryrz_r_cpu
    real(rkind), allocatable,dimension(:,:,:,:) :: wbuf_rxlylz_r_cpu,wbuf_rxlyrz_r_cpu,wbuf_rxrylz_r_cpu,wbuf_rxryrz_r_cpu

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
  endtype base_cpu_object

contains
  subroutine alloc(self)
    class(base_cpu_object), intent(inout) :: self

    associate(nx => self%field%nx, ny => self%field%ny, nz => self%field%nz, ng => self%field%grid%ng, nv => self%field%nv)

      self%w_cpu => self%field%w 

      allocate(self%wbuf1s_cpu(ng,ny,nz,nv))
      allocate(self%wbuf2s_cpu(ng,ny,nz,nv))
      allocate(self%wbuf3s_cpu(nx,ng,nz,nv))
      allocate(self%wbuf4s_cpu(nx,ng,nz,nv))
      allocate(self%wbuf5s_cpu(nx,ny,ng,nv))
      allocate(self%wbuf6s_cpu(nx,ny,ng,nv))
      allocate(self%wbuf1r_cpu(ng,ny,nz,nv))
      allocate(self%wbuf2r_cpu(ng,ny,nz,nv))
      allocate(self%wbuf3r_cpu(nx,ng,nz,nv))
      allocate(self%wbuf4r_cpu(nx,ng,nz,nv))
      allocate(self%wbuf5r_cpu(nx,ny,ng,nv))
      allocate(self%wbuf6r_cpu(nx,ny,ng,nv))

      allocate(self%wbuf1xybcs_cpu(ng,1-ng:ny+ng,nz,nv))
      allocate(self%wbuf2xybcs_cpu(ng,1-ng:ny+ng,nz,nv))
      allocate(self%wbuf1xybcr_cpu(ng,1-ng:ny+ng,nz,nv))
      allocate(self%wbuf2xybcr_cpu(ng,1-ng:ny+ng,nz,nv))
      allocate(self%x_cpu(1-ng:nx+ng), self%y_cpu(1-ng:ny+ng), self%z_cpu(1-ng:nz+ng))
      allocate(self%yn_cpu(1:ny+1))

      allocate(self%dcsidx_cpu(nx), self%dcsidx2_cpu(nx), self%dcsidxs_cpu(nx))
      allocate(self%detady_cpu(ny), self%detady2_cpu(ny), self%detadys_cpu(ny))
      allocate(self%dzitdz_cpu(nz), self%dzitdz2_cpu(nz), self%dzitdzs_cpu(nz))

      allocate(self%wbuf1s_c_cpu(nv,ng,ny,ng))
      allocate(self%wbuf2s_c_cpu(nv,ng,ny,ng))
      allocate(self%wbuf3s_c_cpu(nv,ng,ny,ng))
      allocate(self%wbuf4s_c_cpu(nv,ng,ny,ng))
      allocate(self%wbuf1r_c_cpu(nv,ng,ny,ng))
      allocate(self%wbuf2r_c_cpu(nv,ng,ny,ng))
      allocate(self%wbuf3r_c_cpu(nv,ng,ny,ng))
      allocate(self%wbuf4r_c_cpu(nv,ng,ny,ng))

      allocate(self%wbuf_lxly_s_cpu(nv,ng,ng,nz), self%wbuf_lxry_s_cpu(nv,ng,ng,nz))
      allocate(self%wbuf_rxly_s_cpu(nv,ng,ng,nz), self%wbuf_rxry_s_cpu(nv,ng,ng,nz))
      allocate(self%wbuf_lxlz_s_cpu(nv,ng,ny,ng), self%wbuf_lxrz_s_cpu(nv,ng,ny,ng))
      allocate(self%wbuf_rxlz_s_cpu(nv,ng,ny,ng), self%wbuf_rxrz_s_cpu(nv,ng,ny,ng))
      allocate(self%wbuf_lylz_s_cpu(nv,nx,ng,ng), self%wbuf_lyrz_s_cpu(nv,nx,ng,ng))
      allocate(self%wbuf_rylz_s_cpu(nv,nx,ng,ng), self%wbuf_ryrz_s_cpu(nv,nx,ng,ng))
      allocate(self%wbuf_lxlylz_s_cpu(nv,ng,ng,ng),self%wbuf_lxlyrz_s_cpu(nv,ng,ng,ng))
      allocate(self%wbuf_lxrylz_s_cpu(nv,ng,ng,ng),self%wbuf_lxryrz_s_cpu(nv,ng,ng,ng))
      allocate(self%wbuf_rxlylz_s_cpu(nv,ng,ng,ng),self%wbuf_rxlyrz_s_cpu(nv,ng,ng,ng))
      allocate(self%wbuf_rxrylz_s_cpu(nv,ng,ng,ng),self%wbuf_rxryrz_s_cpu(nv,ng,ng,ng))

      allocate(self%wbuf_lxly_r_cpu(nv,ng,ng,nz), self%wbuf_lxry_r_cpu(nv,ng,ng,nz))
      allocate(self%wbuf_rxly_r_cpu(nv,ng,ng,nz), self%wbuf_rxry_r_cpu(nv,ng,ng,nz))
      allocate(self%wbuf_lxlz_r_cpu(nv,ng,ny,ng), self%wbuf_lxrz_r_cpu(nv,ng,ny,ng))
      allocate(self%wbuf_rxlz_r_cpu(nv,ng,ny,ng), self%wbuf_rxrz_r_cpu(nv,ng,ny,ng))
      allocate(self%wbuf_lylz_r_cpu(nv,nx,ng,ng), self%wbuf_lyrz_r_cpu(nv,nx,ng,ng))
      allocate(self%wbuf_rylz_r_cpu(nv,nx,ng,ng), self%wbuf_ryrz_r_cpu(nv,nx,ng,ng))
      allocate(self%wbuf_lxlylz_r_cpu(nv,ng,ng,ng),self%wbuf_lxlyrz_r_cpu(nv,ng,ng,ng))
      allocate(self%wbuf_lxrylz_r_cpu(nv,ng,ng,ng),self%wbuf_lxryrz_r_cpu(nv,ng,ng,ng))
      allocate(self%wbuf_rxlylz_r_cpu(nv,ng,ng,ng),self%wbuf_rxlyrz_r_cpu(nv,ng,ng,ng))
      allocate(self%wbuf_rxrylz_r_cpu(nv,ng,ng,ng),self%wbuf_rxryrz_r_cpu(nv,ng,ng,ng))

      if (self%field%grid%grid_dim == 2) then
        allocate(self%xc2_cpu(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%yc2_cpu(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%dcsidxc2_cpu(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%dcsidyc2_cpu(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%detadxc2_cpu(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%detadyc2_cpu(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%dcsidxnc2_cpu(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%dcsidync2_cpu(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%detadxnc2_cpu(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%detadync2_cpu(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%dxdcsic2_cpu(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%dydcsic2_cpu(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%dxdetac2_cpu(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%dydetac2_cpu(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%dxdcsinc2_cpu(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%dydcsinc2_cpu(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%dxdetanc2_cpu(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%dydetanc2_cpu(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%mcsi_cpu(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%meta_cpu(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%mcsijac1_cpu(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%metajac1_cpu(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%jac_cpu(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%g1_cpu(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%g2_cpu(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%g12_cpu(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%csimod_cpu(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%etamod_cpu(1-ng:nx+ng,1-ng:ny+ng))
        allocate(self%theta_ij_cpu(1:nx,1:ny))

        allocate(self%wall_tag_cpu(1-ng:nx+ng))
        allocate(self%wbuftus_cpu(nz,2))
        allocate(self%wbuftes_cpu(nz,2))
        allocate(self%wbuftur_cpu(nz,2))
        allocate(self%wbufter_cpu(nz,2))
      endif

    endassociate

    self%x_cpu = self%field%x
    self%y_cpu = self%field%y
    self%z_cpu = self%field%z
    self%yn_cpu = self%field%yn

    self%dcsidx_cpu = self%field%dcsidx
    self%dcsidxs_cpu = self%field%dcsidxs
    self%dcsidx2_cpu = self%field%dcsidx2
    self%detady_cpu = self%field%detady
    self%detadys_cpu = self%field%detadys
    self%detady2_cpu = self%field%detady2
    self%dzitdz_cpu = self%field%dzitdz
    self%dzitdzs_cpu = self%field%dzitdzs
    self%dzitdz2_cpu = self%field%dzitdz2

    if (self%field%grid%grid_dim == 2) then
      self%xc2_cpu = self%field%xc2
      self%yc2_cpu = self%field%yc2
      self%dcsidxc2_cpu = self%field%dcsidxc2
      self%dcsidyc2_cpu = self%field%dcsidyc2
      self%detadxc2_cpu = self%field%detadxc2
      self%detadyc2_cpu = self%field%detadyc2
      self%dcsidxnc2_cpu = self%field%dcsidxnc2
      self%dcsidync2_cpu = self%field%dcsidync2
      self%detadxnc2_cpu = self%field%detadxnc2
      self%detadync2_cpu = self%field%detadync2
      self%dxdcsic2_cpu = self%field%dxdcsic2
      self%dydcsic2_cpu = self%field%dydcsic2
      self%dxdetac2_cpu = self%field%dxdetac2
      self%dydetac2_cpu = self%field%dydetac2
      self%dxdcsinc2_cpu = self%field%dxdcsinc2
      self%dydcsinc2_cpu = self%field%dydcsinc2
      self%dxdetanc2_cpu = self%field%dxdetanc2
      self%dydetanc2_cpu = self%field%dydetanc2
      self%mcsi_cpu = self%field%mcsi
      self%meta_cpu = self%field%meta
      self%mcsijac1_cpu = self%field%mcsijac1
      self%metajac1_cpu = self%field%metajac1
      self%jac_cpu = self%field%jac
      self%g1_cpu = self%field%g1
      self%g2_cpu = self%field%g2
      self%g12_cpu = self%field%g12
      self%csimod_cpu = self%field%csimod
      self%etamod_cpu = self%field%etamod
      self%wall_tag_cpu = self%field%wall_tag
      self%theta_ij_cpu = self%field%theta_ij
    endif

  endsubroutine alloc

  subroutine bcte_step_1_subroutine(nx,ny,nz,nv,ng,nrank_x,ite_rank_x,itu_rank_x,ite_l,itu_l,w_cpu,wbuftus_cpu,wbuftes_cpu)
    integer, intent(in) :: nx, ny, nz, nv, ng
    integer, intent(in) :: nrank_x, ite_rank_x, itu_rank_x, ite_l, itu_l
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(in) :: w_cpu
    real(rkind), dimension(nz, 2), intent(out) :: wbuftus_cpu, wbuftes_cpu
    integer :: i,j,k,m,iercuda
    if(nrank_x == ite_rank_x) then
      do k = 1,nz
        wbuftes_cpu(k,1) = w_cpu(ite_l,1,k,1)
        if(nv > 1) wbuftes_cpu(k,2) = w_cpu(ite_l,1,k,5)
      enddo
    endif
    if(nrank_x == itu_rank_x) then
      do k = 1,nz
        wbuftus_cpu(k,1) = w_cpu(itu_l,1,k,1)
        if(nv > 1) wbuftus_cpu(k,2) = w_cpu(itu_l,1,k,5)
      enddo
    endif

  endsubroutine bcte_step_1_subroutine



  subroutine bcte(self, steps)
    class(base_cpu_object), intent(inout) :: self
    logical, dimension(3), optional :: steps
    logical, dimension(3) :: steps_
    integer, dimension(4) :: requests

    steps_ = .true. ; if(present(steps)) steps_ = steps
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv,&
    & w_cpu => self%w_cpu, wbuftes_cpu => self%wbuftes_cpu, wbuftus_cpu => self%wbuftus_cpu,&
    & wbufter_cpu => self%wbufter_cpu, wbuftur_cpu => self%wbuftur_cpu, ite_rank_x => self%field%ite_rank&
    &_x, itu_rank_x => self%field%itu_rank_x, ite_l => self%field%ite_l, itu_l => self%field%itu_l,&
    & nrank_x => self%field%nrank_x, mp_cartx => self%field%mp_cartx, iermpi => self%mpi_err)

      if(steps_(1)) then
        call bcte_step_1_subroutine(nx, ny, nz, nv, ng, nrank_x, ite_rank_x, itu_rank_x, ite_l,&
        & itu_l, w_cpu, wbuftus_cpu, wbuftes_cpu)
      endif

      if(steps_(2)) then
        if (nrank_x == ite_rank_x .and. nrank_x == itu_rank_x) then
          wbufter_cpu = wbuftes_cpu
          wbuftur_cpu = wbuftus_cpu
        else
          requests(:) = mpi_request_null
          if (nrank_x == ite_rank_x) then
            call mpi_isend(wbuftes_cpu,2*nz,mpi_prec,itu_rank_x,1,mp_cartx,requests(1),iermpi)
            call mpi_irecv(wbuftur_cpu,2*nz,mpi_prec,itu_rank_x,1,mp_cartx,requests(2),iermpi)
          elseif (nrank_x == itu_rank_x) then
            call mpi_isend(wbuftus_cpu,2*nz,mpi_prec,ite_rank_x,1,mp_cartx,requests(3),iermpi)
            call mpi_irecv(wbufter_cpu,2*nz,mpi_prec,ite_rank_x,1,mp_cartx,requests(4),iermpi)
          endif
          call mpi_waitall(4,requests,mpi_statuses_ignore,iermpi)
        endif
      endif

      if(steps_(3)) then
        call bcte_step_3_subroutine(nx, ny, nz, nv, ng, nrank_x, ite_rank_x, itu_rank_x, ite_l,&
        & itu_l, w_cpu, wbuftur_cpu, wbufter_cpu)
      endif
    endassociate
  endsubroutine bcte

  subroutine bcte_var(self, w_swap, steps)
    class(base_cpu_object), intent(inout) :: self
    real(rkind), dimension(1-self%ng:self%nx+self%ng, 1-self%ng:self%ny+self%ng, 1-self%ng:self%nz+self%ng,1:1) :: w_swap
    logical, dimension(3), optional :: steps
    logical, dimension(3) :: steps_
    integer, dimension(4) :: requests

    steps_ = .true. ; if(present(steps)) steps_ = steps
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv,&
    & w_cpu => self%w_cpu, wbuftes_cpu => self%wbuftes_cpu, wbuftus_cpu => self%wbuftus_cpu,&
    & wbufter_cpu => self%wbufter_cpu, wbuftur_cpu => self%wbuftur_cpu, ite_rank_x => self%field%ite_rank&
    &_x, itu_rank_x => self%field%itu_rank_x, ite_l => self%field%ite_l, itu_l => self%field%itu_l,&
    & nrank_x => self%field%nrank_x, mp_cartx => self%field%mp_cartx, iermpi => self%mpi_err)

      if(steps_(1)) then
        call bcte_step_1_subroutine(nx, ny, nz, 1, ng, nrank_x, ite_rank_x, itu_rank_x, ite_l,&
        & itu_l, w_swap, wbuftus_cpu, wbuftes_cpu)
      endif

      if(steps_(2)) then
        if (nrank_x == ite_rank_x .and. nrank_x == itu_rank_x) then
          wbufter_cpu = wbuftes_cpu
          wbuftur_cpu = wbuftus_cpu
        else
          requests(:) = mpi_request_null
          if (nrank_x == ite_rank_x) then
            call mpi_isend(wbuftes_cpu,1*nz,mpi_prec,itu_rank_x,3,mp_cartx,requests(1),iermpi)
            call mpi_irecv(wbuftur_cpu,1*nz,mpi_prec,itu_rank_x,3,mp_cartx,requests(2),iermpi)
          elseif (nrank_x == itu_rank_x) then
            call mpi_isend(wbuftus_cpu,1*nz,mpi_prec,ite_rank_x,3,mp_cartx,requests(3),iermpi)
            call mpi_irecv(wbufter_cpu,1*nz,mpi_prec,ite_rank_x,3,mp_cartx,requests(4),iermpi)
          endif
          call mpi_waitall(4,requests,mpi_statuses_ignore,iermpi)
        endif
      endif

      if(steps_(3)) then
        call bcte_step_3_subroutine(nx, ny, nz, 1, ng, nrank_x, ite_rank_x, itu_rank_x, ite_l,&
        & itu_l, w_swap, wbuftur_cpu, wbufter_cpu)
      endif
    endassociate
  endsubroutine bcte_var

  subroutine bcte_step_3_subroutine(nx,ny,nz,nv,ng,nrank_x,ite_rank_x,itu_rank_x,ite_l,itu_l,w_cpu,wbuftur_cpu,wbufter_cpu)
    integer, intent(in) :: nx, ny, nz, nv, ng
    integer, intent(in) :: nrank_x, ite_rank_x, itu_rank_x, ite_l, itu_l
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(inout) :: w_cpu
    real(rkind), dimension(nz, 2), intent(in) :: wbuftur_cpu, wbufter_cpu
    integer :: k, iercuda
    if(nrank_x == ite_rank_x) then
      do k = 1,nz
        w_cpu(ite_l,1,k,1) = 0.5_rkind*(w_cpu(ite_l,1,k,1)+wbuftur_cpu(k,1))
        if(nv > 1) w_cpu(ite_l,1,k,5) = 0.5_rkind*(w_cpu(ite_l,1,k,5)+wbuftur_cpu(k,2))
      enddo
    endif

    if (nrank_x==itu_rank_x) then
      do k = 1,nz
        w_cpu(itu_l,1,k,1) = 0.5_rkind*(w_cpu(itu_l,1,k,1)+wbufter_cpu(k,1))
        if(nv > 1) w_cpu(itu_l,1,k,5) = 0.5_rkind*(w_cpu(itu_l,1,k,5)+wbufter_cpu(k,2))
      enddo
    endif

  endsubroutine bcte_step_3_subroutine



  subroutine bcswap_wake_step_1_subroutine(nx,ny,nz,nv,ng,w_cpu,wbuf4s_cpu)
    integer :: nx, ny, nz, ng, nv
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: w_cpu
    real(rkind), dimension(1:,1:,1:,1:) :: wbuf4s_cpu
    integer :: i,j,k,m,iercuda

    do k = 1,nz
      do j = 1,ng
        do i = 1,nx
          do m=1,nv
            wbuf4s_cpu(i,j,k,m) = w_cpu(nx-i+1,1+j,k,m)
          enddo
        enddo
      enddo
    enddo


  endsubroutine bcswap_wake_step_1_subroutine



  subroutine bcswap_wake_step_3_subroutine(nx,ny,nz,nv,ng,w_cpu,wbuf3r_cpu,wall_tag_cpu)
    integer :: nx, ny, nz, ng, nv
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: w_cpu
    real(rkind), dimension(1:,1:,1:,1:) :: wbuf3r_cpu
    integer, dimension(1-ng:nx+ng) :: wall_tag_cpu
    integer :: i,j,k,m, iercuda

    do k = 1,nz
      do j = 1,ng
        do i = 1,nx
          do m=1,nv
            if (wall_tag_cpu(i) > 0) then
              w_cpu(i,1-j,k,m) = wbuf3r_cpu(i,j,k,m)
            endif
          enddo
        enddo
      enddo
    enddo


  endsubroutine bcswap_wake_step_3_subroutine



  subroutine bcswap_step_1_subroutine(nx,ny,nz,nv,ng,ndim,w_cpu,wbuf1s_cpu,wbuf2s_cpu,wbuf3s_cpu,wbuf4s_cpu,wbuf5s_cpu,wbuf6s_cpu)
    integer :: nx, ny, nz, ng, nv, ndim
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: w_cpu
    real(rkind), dimension(1:,1:,1:,1:) :: wbuf1s_cpu, wbuf2s_cpu, wbuf3s_cpu, wbuf4s_cpu, wbuf5s_cpu, wbuf6s_cpu
    integer :: i,j,k,m,iercuda

    do k = 1,nz
      do j = 1,ny
        do i = 1,ng
          do m=1,nv
            wbuf1s_cpu(i,j,k,m) = w_cpu(i,j,k,m)
            wbuf2s_cpu(i,j,k,m) = w_cpu(nx-ng+i,j,k,m)
          enddo
        enddo
      enddo
    enddo

    do k = 1,nz
      do j = 1,ng
        do i = 1,nx
          do m=1,nv
            wbuf3s_cpu(i,j,k,m) = w_cpu(i,j,k,m)
            wbuf4s_cpu(i,j,k,m) = w_cpu(i,ny-ng+j,k,m)
          enddo
        enddo
      enddo
    enddo
    if (ndim==3) then
      do k = 1,ng
        do j = 1,ny
          do i = 1,nx
            do m=1,nv
              wbuf5s_cpu(i,j,k,m) = w_cpu(i,j,k,m)
              wbuf6s_cpu(i,j,k,m) = w_cpu(i,j,nz-ng+k,m)
            enddo
          enddo
        enddo
      enddo
    endif

  endsubroutine bcswap_step_1_subroutine



  subroutine bcswap_c2_step_1_subroutine(nx,ny,nz,nv,ng,ndim,w_cpu,wbuf1s_cpu,wbuf2s_cpu,wbuf3s_cpu,&
  &wbuf4s_cpu,wbuf5s_cpu,wbuf6s_cpu,dcsidxnc2_cpu,dcsidync2_cpu,detadxnc2_cpu,detadync2_cpu,&
  &is_periodic)
    integer :: nx, ny, nz, ng, nv, ndim
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: w_cpu
    real(rkind), dimension(1:,1:,1:,1:) :: wbuf1s_cpu, wbuf2s_cpu, wbuf3s_cpu, wbuf4s_cpu, wbuf5s_cpu, wbuf6s_cpu
    real(rkind), dimension(1-ng:,1-ng:) :: dcsidxnc2_cpu, dcsidync2_cpu, detadxnc2_cpu, detadync2_cpu
    integer, dimension(3) :: is_periodic
    integer :: i,j,k,m,iercuda
    if (is_periodic(1) == 1) then
      do k = 1,nz
        do j = 1,ny
          do i = 1,ng
            wbuf1s_cpu(i,j,k,1) = w_cpu(i,j,k,1)
            wbuf1s_cpu(i,j,k,2) = w_cpu(i,j,k,2)*dcsidxnc2_cpu(i,j)+w_cpu(i,j,k,3)*dcsidync2_cpu(i,j)
            wbuf1s_cpu(i,j,k,3) = w_cpu(i,j,k,2)*detadxnc2_cpu(i,j)+w_cpu(i,j,k,3)*detadync2_cpu(i,j)
            wbuf1s_cpu(i,j,k,4) = w_cpu(i,j,k,4)
            wbuf1s_cpu(i,j,k,5) = w_cpu(i,j,k,5)
            wbuf2s_cpu(i,j,k,1) = w_cpu(nx-ng+i,j,k,1)
            wbuf2s_cpu(i,j,k,2) = w_cpu(nx-ng+i,j,k,2)*dcsidxnc2_cpu(nx-ng+i,j)+w_cpu(nx-ng+i,j,k,3)*dcsidync2_cpu(nx-ng+i,j)
            wbuf2s_cpu(i,j,k,3) = w_cpu(nx-ng+i,j,k,2)*detadxnc2_cpu(nx-ng+i,j)+w_cpu(nx-ng+i,j,k,3)*detadync2_cpu(nx-ng+i,j)
            wbuf2s_cpu(i,j,k,4) = w_cpu(nx-ng+i,j,k,4)
            wbuf2s_cpu(i,j,k,5) = w_cpu(nx-ng+i,j,k,5)
          enddo
        enddo
      enddo
    else
      do k = 1,nz
        do j = 1,ny
          do i = 1,ng
            do m=1,nv
              wbuf1s_cpu(i,j,k,m) = w_cpu(i,j,k,m)
              wbuf2s_cpu(i,j,k,m) = w_cpu(nx-ng+i,j,k,m)
            enddo
          enddo
        enddo
      enddo
    endif
    do k = 1,nz
      do j = 1,ng
        do i = 1,nx
          do m=1,nv
            wbuf3s_cpu(i,j,k,m) = w_cpu(i,j,k,m)
            wbuf4s_cpu(i,j,k,m) = w_cpu(i,ny-ng+j,k,m)
          enddo
        enddo
      enddo
    enddo
    if (ndim==3) then
      do k = 1,ng
        do j = 1,ny
          do i = 1,nx
            do m=1,nv
              wbuf5s_cpu(i,j,k,m) = w_cpu(i,j,k,m)
              wbuf6s_cpu(i,j,k,m) = w_cpu(i,j,nz-ng+k,m)
            enddo
          enddo
        enddo
      enddo
    endif

  endsubroutine bcswap_c2_step_1_subroutine



  subroutine bcswap_c2xybc_step_1_subroutine(nx,ny,nz,nv,ng,ndim,w_cpu,wbuf1xybcs_cpu,&
  &wbuf2xybcs_cpu,wbuf3s_cpu,wbuf4s_cpu,wbuf5s_cpu,wbuf6s_cpu,dcsidxnc2_cpu,dcsidync2_cpu,&
  &detadxnc2_cpu,detadync2_cpu,is_periodic)
    integer :: nx, ny, nz, ng, nv, ndim
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: w_cpu
    real(rkind), dimension(1:,1-ng:,1:,1:) :: wbuf1xybcs_cpu, wbuf2xybcs_cpu
    real(rkind), dimension(1:,1:,1:,1:) :: wbuf3s_cpu, wbuf4s_cpu, wbuf5s_cpu, wbuf6s_cpu
    real(rkind), dimension(1-ng:,1-ng:) :: dcsidxnc2_cpu, dcsidync2_cpu, detadxnc2_cpu, detadync2_cpu
    integer, dimension(3) :: is_periodic
    integer :: i,j,k,m,iercuda
    if (is_periodic(1) == 1) then
      do k = 1,nz
        do j = 1-ng,ny+ng
          do i = 1,ng
            wbuf1xybcs_cpu(i,j,k,1) = w_cpu(i,j,k,1)
            wbuf1xybcs_cpu(i,j,k,2) = w_cpu(i,j,k,2)*dcsidxnc2_cpu(i,j)+w_cpu(i,j,k,3)*dcsidync2_cpu(i,j)
            wbuf1xybcs_cpu(i,j,k,3) = w_cpu(i,j,k,2)*detadxnc2_cpu(i,j)+w_cpu(i,j,k,3)*detadync2_cpu(i,j)
            wbuf1xybcs_cpu(i,j,k,4) = w_cpu(i,j,k,4)
            wbuf1xybcs_cpu(i,j,k,5) = w_cpu(i,j,k,5)
            wbuf2xybcs_cpu(i,j,k,1) = w_cpu(nx-ng+i,j,k,1)
            wbuf2xybcs_cpu(i,j,k,2) = w_cpu(nx-ng+i,j,k,2)*dcsidxnc2_cpu(nx-ng+i,j)+w_cpu(nx-ng+i,j,k,3)*dcsidync2_cpu(nx-ng+i,j)
            wbuf2xybcs_cpu(i,j,k,3) = w_cpu(nx-ng+i,j,k,2)*detadxnc2_cpu(nx-ng+i,j)+w_cpu(nx-ng+i,j,k,3)*detadync2_cpu(nx-ng+i,j)
            wbuf2xybcs_cpu(i,j,k,4) = w_cpu(nx-ng+i,j,k,4)
            wbuf2xybcs_cpu(i,j,k,5) = w_cpu(nx-ng+i,j,k,5)
          enddo
        enddo
      enddo
    else
      do k = 1,nz
        do j = 1-ng,ny+ng
          do i = 1,ng
            do m=1,nv
              wbuf1xybcs_cpu(i,j,k,m) = w_cpu(i,j,k,m)
              wbuf2xybcs_cpu(i,j,k,m) = w_cpu(nx-ng+i,j,k,m)
            enddo
          enddo
        enddo
      enddo
    endif
    do k = 1,nz
      do j = 1,ng
        do i = 1,nx
          do m=1,nv
            wbuf3s_cpu(i,j,k,m) = w_cpu(i,j,k,m)
            wbuf4s_cpu(i,j,k,m) = w_cpu(i,ny-ng+j,k,m)
          enddo
        enddo
      enddo
    enddo
    if (ndim==3) then
      do k = 1,ng
        do j = 1,ny
          do i = 1,nx
            do m=1,nv
              wbuf5s_cpu(i,j,k,m) = w_cpu(i,j,k,m)
              wbuf6s_cpu(i,j,k,m) = w_cpu(i,j,nz-ng+k,m)
            enddo
          enddo
        enddo
      enddo
    endif

  endsubroutine bcswap_c2xybc_step_1_subroutine



  subroutine bcswap_corner_step_1_subroutine(nx,ny,nz,nv,ng,w_cpu,wbuf1s_c_cpu,wbuf2s_c_cpu,wbuf3s_c_cpu,wbuf4s_c_cpu)
    integer :: nx, ny, nz, ng, nv
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: w_cpu
    real(rkind), dimension(1:,1:,1:,1:) :: wbuf1s_c_cpu, wbuf2s_c_cpu, wbuf3s_c_cpu, wbuf4s_c_cpu
    integer :: i,j,k,m,iercuda

    do j = 1,ny
      do m = 1,nv
        do k=1,ng
          do i=1,ng
            wbuf1s_c_cpu(m,i,j,k) = w_cpu(i,j,k,m)
            wbuf2s_c_cpu(m,i,j,k) = w_cpu(nx-ng+i,j,nz-ng+k,m)
            wbuf3s_c_cpu(m,i,j,k) = w_cpu(i,j,nz-ng+k,m)
            wbuf4s_c_cpu(m,i,j,k) = w_cpu(nx-ng+i,j,k,m)
          enddo
        enddo
      enddo
    enddo


  endsubroutine bcswap_corner_step_1_subroutine



  subroutine bcswap_corner_step_1b_subroutine(nx,ny,nz,nv,ng,w_cpu,wbuf_lxly_s_cpu,wbuf_lxry_s_cpu,&
  &wbuf_rxly_s_cpu,wbuf_rxry_s_cpu,wbuf_lylz_s_cpu,wbuf_lyrz_s_cpu,wbuf_rylz_s_cpu,wbuf_ryrz_s_cpu,&
  &wbuf_lxlylz_s_cpu,wbuf_lxlyrz_s_cpu,wbuf_lxrylz_s_cpu,wbuf_lxryrz_s_cpu,wbuf_rxlylz_s_cpu,&
  &wbuf_rxlyrz_s_cpu,wbuf_rxrylz_s_cpu,wbuf_rxryrz_s_cpu)
    integer :: nx, ny, nz, ng, nv
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: w_cpu
    real(rkind), dimension(1:,1:,1:,1:) :: wbuf_lxly_s_cpu, wbuf_lxry_s_cpu, wbuf_rxly_s_cpu, wbuf_rxry_s_cpu
    real(rkind), dimension(1:,1:,1:,1:) :: wbuf_lylz_s_cpu, wbuf_lyrz_s_cpu, wbuf_rylz_s_cpu, wbuf_ryrz_s_cpu
    real(rkind), dimension(1:,1:,1:,1:) :: wbuf_lxlylz_s_cpu, wbuf_lxlyrz_s_cpu, wbuf_lxrylz_s_cpu, wbuf_lxryrz_s_cpu
    real(rkind), dimension(1:,1:,1:,1:) :: wbuf_rxlylz_s_cpu, wbuf_rxlyrz_s_cpu, wbuf_rxrylz_s_cpu, wbuf_rxryrz_s_cpu
    integer :: i,j,k,m,iercuda

    do k = 1,nz
      do m = 1,nv
        do j=1,ng
          do i=1,ng
            wbuf_lxly_s_cpu(m,i,j,k) = w_cpu(i,j,k,m)
            wbuf_lxry_s_cpu(m,i,j,k) = w_cpu(i,ny-ng+j,k,m)
            wbuf_rxly_s_cpu(m,i,j,k) = w_cpu(nx-ng+i,j,k,m)
            wbuf_rxry_s_cpu(m,i,j,k) = w_cpu(nx-ng+i,ny-ng+j,k,m)
          enddo
        enddo
      enddo
    enddo

    do i = 1,nx
      do m = 1,nv
        do k=1,ng
          do j=1,ng
            wbuf_lylz_s_cpu(m,i,j,k) = w_cpu(i,j,k,m)
            wbuf_lyrz_s_cpu(m,i,j,k) = w_cpu(i,j,nz-ng+k,m)
            wbuf_rylz_s_cpu(m,i,j,k) = w_cpu(i,ny-ng+j,k,m)
            wbuf_ryrz_s_cpu(m,i,j,k) = w_cpu(i,ny-ng+j,nz-ng+k,m)
          enddo
        enddo
      enddo
    enddo

    do k = 1,ng
      do m = 1,nv
        do j=1,ng
          do i=1,ng
            wbuf_lxlylz_s_cpu(m,i,j,k) = w_cpu(i,j,k,m)
            wbuf_lxlyrz_s_cpu(m,i,j,k) = w_cpu(i,j,nz-ng+k,m)
            wbuf_lxrylz_s_cpu(m,i,j,k) = w_cpu(i,ny-ng+j,k,m)
            wbuf_lxryrz_s_cpu(m,i,j,k) = w_cpu(i,ny-ng+j,nz-ng+k,m)
            wbuf_rxlylz_s_cpu(m,i,j,k) = w_cpu(nx-ng+i,j,k,m)
            wbuf_rxlyrz_s_cpu(m,i,j,k) = w_cpu(nx-ng+i,j,nz-ng+k,m)
            wbuf_rxrylz_s_cpu(m,i,j,k) = w_cpu(nx-ng+i,ny-ng+j,k,m)
            wbuf_rxryrz_s_cpu(m,i,j,k) = w_cpu(nx-ng+i,ny-ng+j,nz-ng+k,m)
          enddo
        enddo
      enddo
    enddo


  endsubroutine bcswap_corner_step_1b_subroutine



  subroutine bcswap_step_3_subroutine(nx,ny,nz,nv,ng,ndim,w_cpu,wbuf1r_cpu,wbuf2r_cpu,wbuf3r_cpu,&
  &wbuf4r_cpu,wbuf5r_cpu,wbuf6r_cpu,ileftx,ilefty,ileftz,irightx,irighty,irightz)
    integer :: nx, ny, nz, ng, nv, ndim
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: w_cpu
    real(rkind), dimension(1:,1:,1:,1:) :: wbuf1r_cpu, wbuf2r_cpu, wbuf3r_cpu, wbuf4r_cpu, wbuf5r_cpu, wbuf6r_cpu
    integer :: i,j,k,m, iercuda
    integer :: ileftx, ilefty, ileftz, irightx, irighty, irightz
    if (ileftx/=mpi_proc_null) then
      do k = 1,nz
        do j = 1,ny
          do i = 1,ng
            do m=1,nv
              w_cpu(i-ng,j,k,m) = wbuf1r_cpu(i,j,k,m)
            enddo
          enddo
        enddo
      enddo
    endif
    if (irightx/=mpi_proc_null) then
      do k = 1,nz
        do j = 1,ny
          do i = 1,ng
            do m=1,nv
              w_cpu(nx+i,j,k,m) = wbuf2r_cpu(i,j,k,m)
            enddo
          enddo
        enddo
      enddo
    endif
    if (ilefty/=mpi_proc_null) then
      do k = 1,nz
        do j = 1,ng
          do i = 1,nx
            do m=1,nv
              w_cpu(i,j-ng,k,m) = wbuf3r_cpu(i,j,k,m)
            enddo
          enddo
        enddo
      enddo
    endif
    if (irighty/=mpi_proc_null) then
      do k = 1,nz
        do j = 1,ng
          do i = 1,nx
            do m=1,nv
              w_cpu(i,ny+j,k,m) = wbuf4r_cpu(i,j,k,m)
            enddo
          enddo
        enddo
      enddo
    endif
    if(ndim == 3) then
      if (ileftz/=mpi_proc_null) then
        do k = 1,ng
          do j = 1,ny
            do i = 1,nx
              do m=1,nv
                w_cpu(i,j,k-ng,m) = wbuf5r_cpu(i,j,k,m)
              enddo
            enddo
          enddo
        enddo
      endif
      if (irightz/=mpi_proc_null) then
        do k = 1,ng
          do j = 1,ny
            do i = 1,nx
              do m=1,nv
                w_cpu(i,j,nz+k,m) = wbuf6r_cpu(i,j,k,m)
              enddo
            enddo
          enddo
        enddo
      endif
    endif

  endsubroutine bcswap_step_3_subroutine



  subroutine bcswap_c2_step_3_subroutine(nx,ny,nz,nv,ng,ndim,w_cpu,wbuf1r_cpu,wbuf2r_cpu,wbuf3r_cpu,&
  &wbuf4r_cpu,wbuf5r_cpu,wbuf6r_cpu,ileftx,ilefty,ileftz,irightx,irighty,irightz,dxdcsinc2_cpu,&
  &dydcsinc2_cpu,dxdetanc2_cpu,dydetanc2_cpu,is_periodic)
    integer :: nx, ny, nz, ng, nv, ndim
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: w_cpu
    real(rkind), dimension(1:,1:,1:,1:) :: wbuf1r_cpu, wbuf2r_cpu, wbuf3r_cpu, wbuf4r_cpu, wbuf5r_cpu, wbuf6r_cpu
    real(rkind), dimension(1-ng:,1-ng:) :: dxdcsinc2_cpu, dydcsinc2_cpu, dxdetanc2_cpu, dydetanc2_cpu
    integer :: i,j,k,m, iercuda
    integer :: ileftx, ilefty, ileftz, irightx, irighty, irightz
    integer, dimension(3) :: is_periodic
    if (is_periodic(1) == 1) then
      if (ileftx/=mpi_proc_null) then
        do k = 1,nz
          do j = 1,ny
            do i = 1,ng
              w_cpu(i-ng,j,k,1) = wbuf1r_cpu(i,j,k,1)
              w_cpu(i-ng,j,k,2) = wbuf1r_cpu(i,j,k,2)*dxdcsinc2_cpu(i-ng,j)+wbuf1r_cpu(i,j,k,3)*dxdetanc2_cpu(i-ng,j)
              w_cpu(i-ng,j,k,3) = wbuf1r_cpu(i,j,k,2)*dydcsinc2_cpu(i-ng,j)+wbuf1r_cpu(i,j,k,3)*dydetanc2_cpu(i-ng,j)
              w_cpu(i-ng,j,k,4) = wbuf1r_cpu(i,j,k,4)
              w_cpu(i-ng,j,k,5) = wbuf1r_cpu(i,j,k,5)
            enddo
          enddo
        enddo
      endif
      if (irightx/=mpi_proc_null) then
        do k = 1,nz
          do j = 1,ny
            do i = 1,ng
              w_cpu(nx+i,j,k,1) = wbuf2r_cpu(i,j,k,1)
              w_cpu(nx+i,j,k,2) = wbuf2r_cpu(i,j,k,2)*dxdcsinc2_cpu(nx+i,j)+wbuf2r_cpu(i,j,k,3)*dxdetanc2_cpu(nx+i,j)
              w_cpu(nx+i,j,k,3) = wbuf2r_cpu(i,j,k,2)*dydcsinc2_cpu(nx+i,j)+wbuf2r_cpu(i,j,k,3)*dydetanc2_cpu(nx+i,j)
              w_cpu(nx+i,j,k,4) = wbuf2r_cpu(i,j,k,4)
              w_cpu(nx+i,j,k,5) = wbuf2r_cpu(i,j,k,5)
            enddo
          enddo
        enddo
      endif
    else
      if (ileftx/=mpi_proc_null) then
        do k = 1,nz
          do j = 1,ny
            do i = 1,ng
              do m=1,nv
                w_cpu(i-ng,j,k,m) = wbuf1r_cpu(i,j,k,m)
              enddo
            enddo
          enddo
        enddo
      endif
      if (irightx/=mpi_proc_null) then
        do k = 1,nz
          do j = 1,ny
            do i = 1,ng
              do m=1,nv
                w_cpu(nx+i,j,k,m) = wbuf2r_cpu(i,j,k,m)
              enddo
            enddo
          enddo
        enddo
      endif
    endif
    if (ilefty/=mpi_proc_null) then
      do k = 1,nz
        do j = 1,ng
          do i = 1,nx
            do m=1,nv
              w_cpu(i,j-ng,k,m) = wbuf3r_cpu(i,j,k,m)
            enddo
          enddo
        enddo
      enddo
    endif
    if (irighty/=mpi_proc_null) then
      do k = 1,nz
        do j = 1,ng
          do i = 1,nx
            do m=1,nv
              w_cpu(i,ny+j,k,m) = wbuf4r_cpu(i,j,k,m)
            enddo
          enddo
        enddo
      enddo
    endif
    if(ndim == 3) then
      if (ileftz/=mpi_proc_null) then
        do k = 1,ng
          do j = 1,ny
            do i = 1,nx
              do m=1,nv
                w_cpu(i,j,k-ng,m) = wbuf5r_cpu(i,j,k,m)
              enddo
            enddo
          enddo
        enddo
      endif
      if (irightz/=mpi_proc_null) then
        do k = 1,ng
          do j = 1,ny
            do i = 1,nx
              do m=1,nv
                w_cpu(i,j,nz+k,m) = wbuf6r_cpu(i,j,k,m)
              enddo
            enddo
          enddo
        enddo
      endif
    endif

  endsubroutine bcswap_c2_step_3_subroutine



  subroutine bcswap_c2xybc_step_3_subroutine(nx,ny,nz,nv,ng,ndim,w_cpu,wbuf1xybcr_cpu,&
  &wbuf2xybcr_cpu,wbuf3r_cpu,wbuf4r_cpu,wbuf5r_cpu,wbuf6r_cpu,ileftx,ilefty,ileftz,irightx,irighty,&
  &irightz,dxdcsinc2_cpu,dydcsinc2_cpu,dxdetanc2_cpu,dydetanc2_cpu,is_periodic)
    integer :: nx, ny, nz, ng, nv, ndim
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: w_cpu
    real(rkind), dimension(1:,1-ng:,1:,1:) :: wbuf1xybcr_cpu, wbuf2xybcr_cpu
    real(rkind), dimension(1:,1:,1:,1:) :: wbuf3r_cpu, wbuf4r_cpu, wbuf5r_cpu, wbuf6r_cpu
    real(rkind), dimension(1-ng:,1-ng:) :: dxdcsinc2_cpu, dydcsinc2_cpu, dxdetanc2_cpu, dydetanc2_cpu
    integer :: i,j,k,m, iercuda
    integer :: ileftx, ilefty, ileftz, irightx, irighty, irightz
    integer, dimension(3) :: is_periodic
    if (is_periodic(1) == 1) then
      if (ileftx/=mpi_proc_null) then
        do k = 1,nz
          do j = 1-ng,ny+ng
            do i = 1,ng
              w_cpu(i-ng,j,k,1) = wbuf1xybcr_cpu(i,j,k,1)
              w_cpu(i-ng,j,k,2) = wbuf1xybcr_cpu(i,j,k,2)*dxdcsinc2_cpu(i-ng,j)+wbuf1xybcr_cpu(i,j,k,3)*dxdetanc2_cpu(i-ng,j)
              w_cpu(i-ng,j,k,3) = wbuf1xybcr_cpu(i,j,k,2)*dydcsinc2_cpu(i-ng,j)+wbuf1xybcr_cpu(i,j,k,3)*dydetanc2_cpu(i-ng,j)
              w_cpu(i-ng,j,k,4) = wbuf1xybcr_cpu(i,j,k,4)
              w_cpu(i-ng,j,k,5) = wbuf1xybcr_cpu(i,j,k,5)
            enddo
          enddo
        enddo
      endif
      if (irightx/=mpi_proc_null) then
        do k = 1,nz
          do j = 1-ng,ny+ng
            do i = 1,ng
              w_cpu(nx+i,j,k,1) = wbuf2xybcr_cpu(i,j,k,1)
              w_cpu(nx+i,j,k,2) = wbuf2xybcr_cpu(i,j,k,2)*dxdcsinc2_cpu(nx+i,j)+wbuf2xybcr_cpu(i,j,k,3)*dxdetanc2_cpu(nx+i,j)
              w_cpu(nx+i,j,k,3) = wbuf2xybcr_cpu(i,j,k,2)*dydcsinc2_cpu(nx+i,j)+wbuf2xybcr_cpu(i,j,k,3)*dydetanc2_cpu(nx+i,j)
              w_cpu(nx+i,j,k,4) = wbuf2xybcr_cpu(i,j,k,4)
              w_cpu(nx+i,j,k,5) = wbuf2xybcr_cpu(i,j,k,5)
            enddo
          enddo
        enddo
      endif
    else
      if (ileftx/=mpi_proc_null) then
        do k = 1,nz
          do j = 1-ng,ny+ng
            do i = 1,ng
              do m=1,nv
                w_cpu(i-ng,j,k,m) = wbuf1xybcr_cpu(i,j,k,m)
              enddo
            enddo
          enddo
        enddo
      endif
      if (irightx/=mpi_proc_null) then
        do k = 1,nz
          do j = 1-ng,ny+ng
            do i = 1,ng
              do m=1,nv
                w_cpu(nx+i,j,k,m) = wbuf2xybcr_cpu(i,j,k,m)
              enddo
            enddo
          enddo
        enddo
      endif
    endif
    if (ilefty/=mpi_proc_null) then
      do k = 1,nz
        do j = 1,ng
          do i = 1,nx
            do m=1,nv
              w_cpu(i,j-ng,k,m) = wbuf3r_cpu(i,j,k,m)
            enddo
          enddo
        enddo
      enddo
    endif
    if (irighty/=mpi_proc_null) then
      do k = 1,nz
        do j = 1,ng
          do i = 1,nx
            do m=1,nv
              w_cpu(i,ny+j,k,m) = wbuf4r_cpu(i,j,k,m)
            enddo
          enddo
        enddo
      enddo
    endif
    if(ndim == 3) then
      if (ileftz/=mpi_proc_null) then
        do k = 1,ng
          do j = 1,ny
            do i = 1,nx
              do m=1,nv
                w_cpu(i,j,k-ng,m) = wbuf5r_cpu(i,j,k,m)
              enddo
            enddo
          enddo
        enddo
      endif
      if (irightz/=mpi_proc_null) then
        do k = 1,ng
          do j = 1,ny
            do i = 1,nx
              do m=1,nv
                w_cpu(i,j,nz+k,m) = wbuf6r_cpu(i,j,k,m)
              enddo
            enddo
          enddo
        enddo
      endif
    endif

  endsubroutine bcswap_c2xybc_step_3_subroutine



  subroutine bcswap_corner_step_3_subroutine(nx,ny,nz,nv,ng,w_cpu,wbuf1r_c_cpu,wbuf2r_c_cpu,&
  &wbuf3r_c_cpu,wbuf4r_c_cpu,ileftbottom,ilefttop,irightbottom,irighttop)
    integer :: nx, ny, nz, ng, nv
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: w_cpu
    real(rkind), dimension(1:,1:,1:,1:) :: wbuf1r_c_cpu, wbuf2r_c_cpu, wbuf3r_c_cpu, wbuf4r_c_cpu
    integer :: i,j,k,m, iercuda
    integer :: ileftbottom, ilefttop, irightbottom, irighttop
    if (ileftbottom/=mpi_proc_null) then
      do j = 1,ny
        do m = 1,nv
          do k=1,ng
            do i=1,ng
              w_cpu(i-ng,j,k-ng,m) = wbuf1r_c_cpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif
    if (irighttop/=mpi_proc_null) then
      do j = 1,ny
        do m = 1,nv
          do k=1,ng
            do i=1,ng
              w_cpu(i+nx,j,k+nz,m) = wbuf2r_c_cpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif
    if (ilefttop/=mpi_proc_null) then
      do j = 1,ny
        do m = 1,nv
          do k=1,ng
            do i=1,ng
              w_cpu(i-ng,j,k+nz,m) = wbuf3r_c_cpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif
    if (irightbottom/=mpi_proc_null) then
      do j = 1,ny
        do m = 1,nv
          do k=1,ng
            do i=1,ng
              w_cpu(i+nx,j,k-ng,m) = wbuf4r_c_cpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif

  endsubroutine bcswap_corner_step_3_subroutine



  subroutine bcswap_corner_step_3b_subroutine(nx,ny,nz,nv,ng,w_cpu,wbuf_lxly_r_cpu,wbuf_lxry_r_cpu,&
  &wbuf_rxly_r_cpu,wbuf_rxry_r_cpu,wbuf_lylz_r_cpu,wbuf_lyrz_r_cpu,wbuf_rylz_r_cpu,wbuf_ryrz_r_cpu,&
  &wbuf_lxlylz_r_cpu,wbuf_lxlyrz_r_cpu,wbuf_lxrylz_r_cpu,wbuf_lxryrz_r_cpu,wbuf_rxlylz_r_cpu,&
  &wbuf_rxlyrz_r_cpu,wbuf_rxrylz_r_cpu,wbuf_rxryrz_r_cpu,lxly,lxry,rxly,rxry,lylz,lyrz,rylz,ryrz,&
  &lxlylz,lxlyrz,lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz)
    integer :: nx, ny, nz, ng, nv
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: w_cpu
    real(rkind), dimension(1:,1:,1:,1:) :: wbuf_lxly_r_cpu, wbuf_lxry_r_cpu, wbuf_rxly_r_cpu, wbuf_rxry_r_cpu
    real(rkind), dimension(1:,1:,1:,1:) :: wbuf_lylz_r_cpu, wbuf_lyrz_r_cpu, wbuf_rylz_r_cpu, wbuf_ryrz_r_cpu
    real(rkind), dimension(1:,1:,1:,1:) :: wbuf_lxlylz_r_cpu, wbuf_lxlyrz_r_cpu, wbuf_lxrylz_r_cpu, wbuf_lxryrz_r_cpu
    real(rkind), dimension(1:,1:,1:,1:) :: wbuf_rxlylz_r_cpu, wbuf_rxlyrz_r_cpu, wbuf_rxrylz_r_cpu, wbuf_rxryrz_r_cpu
    integer :: i,j,k,m, iercuda
    integer :: lxly, lxry, rxly, rxry
    integer :: lylz, lyrz, rylz, ryrz
    integer :: lxlylz, lxlyrz, lxrylz, lxryrz
    integer :: rxlylz, rxlyrz, rxrylz, rxryrz
    if (lxly/=mpi_proc_null) then
      do k = 1,nz
        do m = 1,nv
          do i=1,ng
            do j=1,ng
              w_cpu(i-ng,j-ng,k,m) = wbuf_lxly_r_cpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif
    if (rxry/=mpi_proc_null) then
      do k = 1,nz
        do m = 1,nv
          do i=1,ng
            do j=1,ng
              w_cpu(i+nx,j+ny,k,m) = wbuf_rxry_r_cpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif
    if (lxry/=mpi_proc_null) then
      do k = 1,nz
        do m = 1,nv
          do i=1,ng
            do j=1,ng
              w_cpu(i-ng,j+ny,k,m) = wbuf_lxry_r_cpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif
    if (rxly/=mpi_proc_null) then
      do k = 1,nz
        do m = 1,nv
          do i=1,ng
            do j=1,ng
              w_cpu(i+nx,j-ng,k,m) = wbuf_rxly_r_cpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif

    if (lylz/=mpi_proc_null) then
      do i = 1,nx
        do m = 1,nv
          do k=1,ng
            do j=1,ng
              w_cpu(i,j-ng,k-ng,m) = wbuf_lylz_r_cpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif
    if (ryrz/=mpi_proc_null) then
      do i = 1,nx
        do m = 1,nv
          do k=1,ng
            do j=1,ng
              w_cpu(i,j+ny,k+nz,m) = wbuf_ryrz_r_cpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif
    if (lyrz/=mpi_proc_null) then
      do i = 1,nx
        do m = 1,nv
          do k=1,ng
            do j=1,ng
              w_cpu(i,j-ng,k+nz,m) = wbuf_lyrz_r_cpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif
    if (rylz/=mpi_proc_null) then
      do i = 1,nx
        do m = 1,nv
          do k=1,ng
            do j=1,ng
              w_cpu(i,j+ny,k-ng,m) = wbuf_rylz_r_cpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif

    if (lxlylz/=mpi_proc_null) then
      do k = 1,ng
        do m = 1,nv
          do i=1,ng
            do j=1,ng
              w_cpu(i-ng,j-ng,k-ng,m) = wbuf_lxlylz_r_cpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif
    if (lxlyrz/=mpi_proc_null) then
      do k = 1,ng
        do m = 1,nv
          do i=1,ng
            do j=1,ng
              w_cpu(i-ng,j-ng,k+nz,m) = wbuf_lxlyrz_r_cpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif
    if (lxrylz/=mpi_proc_null) then
      do k = 1,ng
        do m = 1,nv
          do i=1,ng
            do j=1,ng
              w_cpu(i-ng,j+ny,k-ng,m) = wbuf_lxrylz_r_cpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif
    if (lxryrz/=mpi_proc_null) then
      do k = 1,ng
        do m = 1,nv
          do i=1,ng
            do j=1,ng
              w_cpu(i-ng,j+ny,k+nz,m) = wbuf_lxryrz_r_cpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif
    if (rxlylz/=mpi_proc_null) then
      do k = 1,ng
        do m = 1,nv
          do i=1,ng
            do j=1,ng
              w_cpu(i+nx,j-ng,k-ng,m) = wbuf_rxlylz_r_cpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif
    if (rxlyrz/=mpi_proc_null) then
      do k = 1,ng
        do m = 1,nv
          do i=1,ng
            do j=1,ng
              w_cpu(i+nx,j-ng,k+nz,m) = wbuf_rxlyrz_r_cpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif
    if (rxrylz/=mpi_proc_null) then
      do k = 1,ng
        do m = 1,nv
          do i=1,ng
            do j=1,ng
              w_cpu(i+nx,j+ny,k-ng,m) = wbuf_rxrylz_r_cpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif
    if (rxryrz/=mpi_proc_null) then
      do k = 1,ng
        do m = 1,nv
          do i=1,ng
            do j=1,ng
              w_cpu(i+nx,j+ny,k+nz,m) = wbuf_rxryrz_r_cpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif

  endsubroutine bcswap_corner_step_3b_subroutine



  subroutine bcswap_corner(self, steps)
    class(base_cpu_object), intent(inout) :: self
    logical, dimension(3), optional :: steps
    logical, dimension(3) :: steps_
    integer :: iercuda, indc
    integer, dimension(mpi_status_size) :: istatus

    steps_ = .true. ; if(present(steps)) steps_ = steps

    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv,&
    & w_cpu => self%w_cpu, wbuf1s_c_cpu => self%wbuf1s_c_cpu, wbuf2s_c_cpu => self%wbuf2s_c_cpu,&
    & wbuf3s_c_cpu => self%wbuf3s_c_cpu, wbuf4s_c_cpu => self%wbuf4s_c_cpu, wbuf1r_c_cpu => self%wbuf1r_c&
    &_cpu, wbuf2r_c_cpu => self%wbuf2r_c_cpu, wbuf3r_c_cpu => self%wbuf3r_c_cpu, wbuf4r_c_cpu => self%wbu&
    &f4r_c_cpu, ileftbottom => self%field%ileftbottom, irightbottom => self%field%irightbottom,&
    & ilefttop => self%field%ilefttop, irighttop => self%field%irighttop, mp_cart => self%field%mp_cart,&
    & iermpi => self%mpi_err, nrank => self%field%myrank)

      if(steps_(1)) then
        call bcswap_corner_step_1_subroutine(nx, ny, nz, nv, ng, w_cpu, wbuf1s_c_cpu, wbuf2s_c_cpu, wbuf3s_c_cpu, wbuf4s_c_cpu)
      endif
      if(steps_(2)) then
        indc = nv*ny*ng*ng
        if (ileftbottom == nrank) then
          wbuf2r_c_cpu = wbuf1s_c_cpu
          wbuf1r_c_cpu = wbuf2s_c_cpu
        else
          call mpi_sendrecv(wbuf1s_c_cpu,indc,mpi_prec,ileftbottom ,1, wbuf2r_c_cpu,indc,mpi_prec,&
          &irighttop ,1,mp_cart,istatus,iermpi)
          call mpi_sendrecv(wbuf2s_c_cpu,indc,mpi_prec,irighttop ,2, wbuf1r_c_cpu,indc,mpi_prec,&
          &ileftbottom ,2,mp_cart,istatus,iermpi)
        endif
        if (ilefttop == nrank) then
          wbuf4r_c_cpu = wbuf3s_c_cpu
          wbuf3r_c_cpu = wbuf4s_c_cpu
        else
          call mpi_sendrecv(wbuf3s_c_cpu,indc,mpi_prec,ilefttop ,3, wbuf4r_c_cpu,indc,mpi_prec,&
          &irightbottom,3,mp_cart,istatus,iermpi)
          call mpi_sendrecv(wbuf4s_c_cpu,indc,mpi_prec,irightbottom,4, wbuf3r_c_cpu,indc,mpi_prec,&
          &ilefttop ,4,mp_cart,istatus,iermpi)
        endif
      endif
      if(steps_(3)) then
        call bcswap_corner_step_3_subroutine(nx, ny, nz, nv, ng, w_cpu, wbuf1r_c_cpu, wbuf2r_c_cpu,&
        & wbuf3r_c_cpu, wbuf4r_c_cpu, ileftbottom, ilefttop, irightbottom, irighttop)
      endif
    endassociate
  endsubroutine bcswap_corner

  subroutine bcswap_corner_var(self, w_swap, steps)
    class(base_cpu_object), intent(inout) :: self
    real(rkind), dimension(1-self%ng:self%nx+self%ng, 1-self%ng:self%ny+self%ng, 1-self%ng:self%nz+self%ng,1:1) :: w_swap
    logical, dimension(3), optional :: steps
    logical, dimension(3) :: steps_
    integer :: iercuda, indc
    integer, dimension(mpi_status_size) :: istatus

    steps_ = .true. ; if(present(steps)) steps_ = steps

    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv,&
    & wbuf1s_c_cpu => self%wbuf1s_c_cpu, wbuf2s_c_cpu => self%wbuf2s_c_cpu, wbuf3s_c_cpu => self%wbuf3s_c&
    &_cpu, wbuf4s_c_cpu => self%wbuf4s_c_cpu, wbuf1r_c_cpu => self%wbuf1r_c_cpu, wbuf2r_c_cpu => self%wbu&
    &f2r_c_cpu, wbuf3r_c_cpu => self%wbuf3r_c_cpu, wbuf4r_c_cpu => self%wbuf4r_c_cpu, ileftbottom => self&
    &%field%ileftbottom, irightbottom => self%field%irightbottom, ilefttop => self%field%ilefttop,&
    & irighttop => self%field%irighttop, mp_cart => self%field%mp_cart, iermpi => self%mpi_err,&
    & nrank => self%field%myrank)

      if(steps_(1)) then
        call bcswap_corner_step_1_subroutine(nx, ny, nz, 1, ng, w_swap, wbuf1s_c_cpu, wbuf2s_c_cpu, wbuf3s_c_cpu, wbuf4s_c_cpu)
      endif
      if(steps_(2)) then
        indc = ny*ng*ng
        if (ileftbottom == nrank) then
          wbuf2r_c_cpu = wbuf1s_c_cpu
          wbuf1r_c_cpu = wbuf2s_c_cpu
        else
          call mpi_sendrecv(wbuf1s_c_cpu,indc,mpi_prec,ileftbottom ,1, wbuf2r_c_cpu,indc,mpi_prec,&
          &irighttop ,1,mp_cart,istatus,iermpi)
          call mpi_sendrecv(wbuf2s_c_cpu,indc,mpi_prec,irighttop ,2, wbuf1r_c_cpu,indc,mpi_prec,&
          &ileftbottom ,2,mp_cart,istatus,iermpi)
        endif
        if (ilefttop == nrank) then
          wbuf4r_c_cpu = wbuf3s_c_cpu
          wbuf3r_c_cpu = wbuf4s_c_cpu
        else
          call mpi_sendrecv(wbuf3s_c_cpu,indc,mpi_prec,ilefttop ,3, wbuf4r_c_cpu,indc,mpi_prec,&
          &irightbottom,3,mp_cart,istatus,iermpi)
          call mpi_sendrecv(wbuf4s_c_cpu,indc,mpi_prec,irightbottom,4, wbuf3r_c_cpu,indc,mpi_prec,&
          &ilefttop ,4,mp_cart,istatus,iermpi)
        endif
      endif
      if(steps_(3)) then
        call bcswap_corner_step_3_subroutine(nx, ny, nz, 1, ng, w_swap, wbuf1r_c_cpu, wbuf2r_c_cpu,&
        & wbuf3r_c_cpu, wbuf4r_c_cpu, ileftbottom, ilefttop, irightbottom, irighttop)
      endif
    endassociate
  endsubroutine bcswap_corner_var

  subroutine bcswap(self, steps, swap_xy_corner_bc)
    class(base_cpu_object), intent(inout) :: self
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
    & w_cpu => self%w_cpu, wbuf1s_cpu => self%wbuf1s_cpu, wbuf2s_cpu => self%wbuf2s_cpu,&
    & wbuf1xybcs_cpu => self%wbuf1xybcs_cpu, wbuf2xybcs_cpu => self%wbuf2xybcs_cpu, wbuf3s_cpu => self%wb&
    &uf3s_cpu, wbuf4s_cpu => self%wbuf4s_cpu, wbuf5s_cpu => self%wbuf5s_cpu, wbuf6s_cpu => self%wbuf6s_cp&
    &u, wbuf1r_cpu => self%wbuf1r_cpu, wbuf2r_cpu => self%wbuf2r_cpu, wbuf1xybcr_cpu => self%wbuf1xybcr_c&
    &pu, wbuf2xybcr_cpu => self%wbuf2xybcr_cpu, wbuf3r_cpu => self%wbuf3r_cpu, wbuf4r_cpu => self%wbuf4r_&
    &cpu, wbuf5r_cpu => self%wbuf5r_cpu, wbuf6r_cpu => self%wbuf6r_cpu, ileftx => self%field%ileftx,&
    &irightx => self%field%irightx, nrank_x => self%field%nrank_x, ilefty => self%field%ilefty,&
    &irighty => self%field%irighty, nrank_y => self%field%nrank_y, ileftz => self%field%ileftz,&
    &irightz => self%field%irightz, nrank_z => self%field%nrank_z, dcsidxnc2_cpu => self%dcsidxnc2_cpu,&
    & dcsidync2_cpu => self%dcsidync2_cpu, detadxnc2_cpu => self%detadxnc2_cpu, detadync2_cpu => self%det&
    &adync2_cpu, dxdcsinc2_cpu => self%dxdcsinc2_cpu, dydcsinc2_cpu => self%dydcsinc2_cpu,&
    & dxdetanc2_cpu => self%dxdetanc2_cpu, dydetanc2_cpu => self%dydetanc2_cpu, mp_cartx => self%field%mp&
    &_cartx, mp_carty => self%field%mp_carty, mp_cartz => self%field%mp_cartz, ncoords => self%field%ncoo&
    &rds, nblocks => self%field%nblocks, ndim => self%field%grid%ndim)

      if(steps_(1)) then
        if(self%field%grid%grid_dim == 1) then
          call bcswap_step_1_subroutine(nx, ny, nz, nv, ng, ndim, w_cpu, wbuf1s_cpu, wbuf2s_cpu,&
          & wbuf3s_cpu, wbuf4s_cpu, wbuf5s_cpu, wbuf6s_cpu)
        elseif(self%field%grid%grid_dim == 2) then
          if(swap_xy_corner_bc_ == 0) then
            call bcswap_c2_step_1_subroutine(nx, ny, nz, nv, ng, ndim, w_cpu, wbuf1s_cpu,&
            & wbuf2s_cpu, wbuf3s_cpu, wbuf4s_cpu, wbuf5s_cpu, wbuf6s_cpu, dcsidxnc2_cpu, dcsidync2_cpu,&
            & detadxnc2_cpu, detadync2_cpu, is_periodic)
          else
            call bcswap_c2xybc_step_1_subroutine(nx, ny, nz, nv, ng, ndim, w_cpu, wbuf1xybcs_cpu,&
            & wbuf2xybcs_cpu, wbuf3s_cpu, wbuf4s_cpu, wbuf5s_cpu, wbuf6s_cpu, dcsidxnc2_cpu, dcsidync2_cpu,&
            & detadxnc2_cpu, detadync2_cpu, is_periodic)
          endif
        endif
      endif
      if(steps_(2)) then
        indy = nv*nx*ng*nz
        indz = nv*nx*ny*ng
        if(self%field%grid%grid_dim == 1 .or. (self%field%grid%grid_dim == 2 .and. swap_xy_corner_bc_ == 0)) then
          indx = nv*ng*ny*nz
          if(ileftx == nrank_x) then
            wbuf2r_cpu = wbuf1s_cpu
            wbuf1r_cpu = wbuf2s_cpu
          else
            call mpi_sendrecv(wbuf1s_cpu,indx,mpi_prec,ileftx ,1,wbuf2r_cpu,indx,mpi_prec,irightx,1,mp_cartx,istatus,self%mpi_err)
            call mpi_sendrecv(wbuf2s_cpu,indx,mpi_prec,irightx,2,wbuf1r_cpu,indx,mpi_prec,ileftx ,2,mp_cartx,istatus,self%mpi_err)
          endif
        else
          indx = nv*ng*(2*ng+ny)*nz
          if(ileftx == nrank_x) then
            wbuf2xybcr_cpu = wbuf1xybcs_cpu
            wbuf1xybcr_cpu = wbuf2xybcs_cpu
          else
            call mpi_sendrecv(wbuf1xybcs_cpu,indx,mpi_prec,ileftx ,1,wbuf2xybcr_cpu,indx,mpi_prec,&
            &irightx,1,mp_cartx,istatus,self%mpi_err)
            call mpi_sendrecv(wbuf2xybcs_cpu,indx,mpi_prec,irightx,2,wbuf1xybcr_cpu,indx,mpi_prec,&
            &ileftx ,2,mp_cartx,istatus,self%mpi_err)
          endif
        endif
        if(ilefty == nrank_y) then
          wbuf4r_cpu = wbuf3s_cpu
          wbuf3r_cpu = wbuf4s_cpu
        else
          call mpi_sendrecv(wbuf3s_cpu,indy,mpi_prec,ilefty ,3,wbuf4r_cpu,indy,mpi_prec,irighty,3,mp_carty,istatus,self%mpi_err)
          call mpi_sendrecv(wbuf4s_cpu,indy,mpi_prec,irighty,4,wbuf3r_cpu,indy,mpi_prec,ilefty ,4,mp_carty,istatus,self%mpi_err)
        endif
        if(ndim == 3) then
          if(ileftz == nrank_z) then
            wbuf6r_cpu = wbuf5s_cpu
            wbuf5r_cpu = wbuf6s_cpu
          else
            call mpi_sendrecv(wbuf5s_cpu,indz,mpi_prec,ileftz ,5,wbuf6r_cpu,indz,mpi_prec,irightz,5,mp_cartz,istatus,self%mpi_err)
            call mpi_sendrecv(wbuf6s_cpu,indz,mpi_prec,irightz,6,wbuf5r_cpu,indz,mpi_prec,ileftz ,6,mp_cartz,istatus,self%mpi_err)
          endif
        endif
      endif
      if(steps_(3)) then
        if(self%field%grid%grid_dim == 1) then
          call bcswap_step_3_subroutine(nx, ny, nz, nv, ng, ndim, w_cpu, wbuf1r_cpu, wbuf2r_cpu,&
          & wbuf3r_cpu, wbuf4r_cpu, wbuf5r_cpu, wbuf6r_cpu, ileftx, ilefty, ileftz, irightx, irighty, irightz)
        elseif(self%field%grid%grid_dim == 2) then
          if(swap_xy_corner_bc_ == 0) then
            call bcswap_c2_step_3_subroutine(nx, ny, nz, nv, ng, ndim, w_cpu, wbuf1r_cpu,&
            & wbuf2r_cpu, wbuf3r_cpu, wbuf4r_cpu, wbuf5r_cpu, wbuf6r_cpu, ileftx, ilefty, ileftz, irightx,&
            & irighty, irightz, dxdcsinc2_cpu, dydcsinc2_cpu, dxdetanc2_cpu, dydetanc2_cpu, is_periodic)
          else
            call bcswap_c2xybc_step_3_subroutine(nx, ny, nz, nv, ng, ndim, w_cpu, wbuf1xybcr_cpu,&
            & wbuf2xybcr_cpu, wbuf3r_cpu, wbuf4r_cpu, wbuf5r_cpu, wbuf6r_cpu, ileftx, ilefty, ileftz, irightx,&
            & irighty, irightz, dxdcsinc2_cpu, dydcsinc2_cpu, dxdetanc2_cpu, dydetanc2_cpu, is_periodic)
            if(is_periodic(1) == 0) call extr_corner_ymin_subroutine(ncoords, nblocks, nx, ny, nz, ng, nv, w_cpu)
          endif
        endif
      endif
    endassociate
  endsubroutine bcswap

  subroutine extr_corner_ymin_subroutine(ncoords,nblocks,nx,ny,nz,ng,nv,w_cpu)
    integer, intent(in) :: nx, ny, nz, ng, nv
    real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:), intent(inout) :: w_cpu
    integer, dimension(3) :: ncoords, nblocks
    integer :: i,j,k,m,iercuda
    if(ncoords(1) == 0) then
      do k = 1,nz
        do j = 1,ng
          do i = 1,ng
            do m=1,nv
              w_cpu(1-i,1-j,k,m) = w_cpu(1-i,1,k,m)
            enddo
          enddo
        enddo
      enddo
    endif
    if(ncoords(1) == nblocks(1)-1) then
      do k = 1,nz
        do j = 1,ng
          do i = 1,ng
            do m=1,nv
              w_cpu(nx+i,1-j,k,m) = w_cpu(nx+i,1,k,m)
            enddo
          enddo
        enddo
      enddo
    endif

  endsubroutine extr_corner_ymin_subroutine



  subroutine bcswap_wake(self, steps)
    class(base_cpu_object), intent(inout) :: self
    logical, dimension(3), optional :: steps
    logical, dimension(3) :: steps_
    integer :: iercuda, indx, indy, indz
    integer, dimension(mpi_status_size) :: istatus

    steps_ = .true. ; if(present(steps)) steps_ = steps

    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv,&
    & w_cpu => self%w_cpu, wbuf4s_cpu => self%wbuf4s_cpu, wbuf3r_cpu => self%wbuf3r_cpu,&
    & wall_tag_cpu => self%wall_tag_cpu, mp_cart => self%field%mp_cart, nrank_wake => self%field%nrank_wa&
    &ke, iermpi => self%mpi_err, nrank => self%myrank)

      if(steps_(1)) then
        call bcswap_wake_step_1_subroutine(nx, ny, nz, nv, ng, w_cpu, wbuf4s_cpu)
      endif
      if(steps_(2)) then
        indy = nv*nx*ng*nz
        if (nrank_wake==nrank) then
          wbuf3r_cpu = wbuf4s_cpu
        else
          call mpi_sendrecv(wbuf4s_cpu,indy,mpi_prec,nrank_wake,4,wbuf3r_cpu,indy,mpi_prec,nrank_wake,4,mp_cart,istatus,iermpi)
        endif
      endif
      if(steps_(3)) then
        call bcswap_wake_step_3_subroutine(nx, ny, nz, nv, ng, w_cpu, wbuf3r_cpu, wall_tag_cpu)
      endif
    endassociate
  endsubroutine bcswap_wake

  subroutine bcswap_wake_var(self, w_swap, steps)
    class(base_cpu_object), intent(inout) :: self
    real(rkind), dimension(1-self%ng:self%nx+self%ng, 1-self%ng:self%ny+self%ng, 1-self%ng:self%nz+self%ng,1:1) :: w_swap
    logical, dimension(3), optional :: steps
    logical, dimension(3) :: steps_
    integer :: iercuda, indx, indy, indz
    integer, dimension(mpi_status_size) :: istatus

    steps_ = .true. ; if(present(steps)) steps_ = steps

    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv,&
    & w_cpu => self%w_cpu, wbuf4s_cpu => self%wbuf4s_cpu, wbuf3r_cpu => self%wbuf3r_cpu,&
    & wall_tag_cpu => self%wall_tag_cpu, mp_cart => self%field%mp_cart, nrank_wake => self%field%nrank_wa&
    &ke, iermpi => self%mpi_err, nrank => self%myrank)

      if(steps_(1)) then
        call bcswap_wake_step_1_subroutine(nx, ny, nz, 1, ng, w_swap, wbuf4s_cpu)
      endif
      if(steps_(2)) then
        indy = 1*nx*ng*nz
        if (nrank_wake==nrank) then
          wbuf3r_cpu = wbuf4s_cpu
        else
          call mpi_sendrecv(wbuf4s_cpu,indy,mpi_prec,nrank_wake,4,wbuf3r_cpu,indy,mpi_prec,nrank_wake,4,mp_cart,istatus,iermpi)
        endif
      endif
      if(steps_(3)) then
        call bcswap_wake_step_3_subroutine(nx, ny, nz, 1, ng, w_swap, wbuf3r_cpu, wall_tag_cpu)
      endif
    endassociate
  endsubroutine bcswap_wake_var

  subroutine bcswap_edges_corners(self, steps)
    class(base_cpu_object), intent(inout) :: self
    logical, dimension(3), optional :: steps
    logical, dimension(3) :: steps_
    integer :: i,j,k,m, iercuda, indc
    integer, dimension(mpi_status_size) :: istatus

    steps_ = .true. ; if(present(steps)) steps_ = steps

    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv,&
    & w_cpu => self%w_cpu, wbuf1s_c_cpu => self%wbuf1s_c_cpu, wbuf2s_c_cpu => self%wbuf2s_c_cpu,&
    & wbuf3s_c_cpu => self%wbuf3s_c_cpu, wbuf4s_c_cpu => self%wbuf4s_c_cpu, wbuf1r_c_cpu => self%wbuf1r_c&
    &_cpu, wbuf2r_c_cpu => self%wbuf2r_c_cpu, wbuf3r_c_cpu => self%wbuf3r_c_cpu, wbuf4r_c_cpu => self%wbu&
    &f4r_c_cpu, ileftbottom => self%field%ileftbottom, irightbottom => self%field%irightbottom,&
    & ilefttop => self%field%ilefttop, irighttop => self%field%irighttop, pbc => self%field%grid%is_xyz_p&
    &eriodic, mp_cart => self%field%mp_cart, iermpi => self%mpi_err, nrank => self%field%myrank,&
    & lxly => self%field%lxly, lxry => self%field%lxry, rxly => self%field%rxly, rxry => self%field%rxry,&
    &lxlz => self%field%lxlz, lxrz => self%field%lxrz, rxlz => self%field%rxlz, rxrz => self%field%rxrz,&
    &lylz => self%field%lylz, lyrz => self%field%lyrz, rylz => self%field%rylz, ryrz => self%field%ryrz,&
    &lxlylz => self%field%lxlylz, lxlyrz => self%field%lxlyrz, lxrylz => self%field%lxrylz,&
    & rxlylz => self%field%rxlylz,lxryrz => self%field%lxryrz, rxlyrz => self%field%rxlyrz,&
    & rxrylz => self%field%rxrylz, rxryrz => self%field%rxryrz,wbuf_lxly_s_cpu => self%wbuf_lxly_s_cpu ,&
    & wbuf_lxry_s_cpu => self%wbuf_lxry_s_cpu, wbuf_rxly_s_cpu => self%wbuf_rxly_s_cpu ,&
    & wbuf_rxry_s_cpu => self%wbuf_rxry_s_cpu, wbuf_lxlz_s_cpu => self%wbuf_lxlz_s_cpu ,&
    & wbuf_lxrz_s_cpu => self%wbuf_lxrz_s_cpu, wbuf_rxlz_s_cpu => self%wbuf_rxlz_s_cpu ,&
    & wbuf_rxrz_s_cpu => self%wbuf_rxrz_s_cpu, wbuf_lylz_s_cpu => self%wbuf_lylz_s_cpu ,&
    & wbuf_lyrz_s_cpu => self%wbuf_lyrz_s_cpu, wbuf_rylz_s_cpu => self%wbuf_rylz_s_cpu ,&
    & wbuf_ryrz_s_cpu => self%wbuf_ryrz_s_cpu, wbuf_lxlylz_s_cpu => self%wbuf_lxlylz_s_cpu ,&
    & wbuf_lxlyrz_s_cpu => self%wbuf_lxlyrz_s_cpu, wbuf_lxrylz_s_cpu => self%wbuf_lxrylz_s_cpu ,&
    & wbuf_lxryrz_s_cpu => self%wbuf_lxryrz_s_cpu, wbuf_rxlylz_s_cpu => self%wbuf_rxlylz_s_cpu ,&
    & wbuf_rxlyrz_s_cpu => self%wbuf_rxlyrz_s_cpu, wbuf_rxrylz_s_cpu => self%wbuf_rxrylz_s_cpu ,&
    & wbuf_rxryrz_s_cpu => self%wbuf_rxryrz_s_cpu, wbuf_lxly_r_cpu => self%wbuf_lxly_r_cpu ,&
    & wbuf_lxry_r_cpu => self%wbuf_lxry_r_cpu, wbuf_rxly_r_cpu => self%wbuf_rxly_r_cpu ,&
    & wbuf_rxry_r_cpu => self%wbuf_rxry_r_cpu, wbuf_lxlz_r_cpu => self%wbuf_lxlz_r_cpu ,&
    & wbuf_lxrz_r_cpu => self%wbuf_lxrz_r_cpu, wbuf_rxlz_r_cpu => self%wbuf_rxlz_r_cpu ,&
    & wbuf_rxrz_r_cpu => self%wbuf_rxrz_r_cpu, wbuf_lylz_r_cpu => self%wbuf_lylz_r_cpu ,&
    & wbuf_lyrz_r_cpu => self%wbuf_lyrz_r_cpu, wbuf_rylz_r_cpu => self%wbuf_rylz_r_cpu ,&
    & wbuf_ryrz_r_cpu => self%wbuf_ryrz_r_cpu, wbuf_lxlylz_r_cpu => self%wbuf_lxlylz_r_cpu ,&
    & wbuf_lxlyrz_r_cpu => self%wbuf_lxlyrz_r_cpu, wbuf_lxrylz_r_cpu => self%wbuf_lxrylz_r_cpu ,&
    & wbuf_lxryrz_r_cpu => self%wbuf_lxryrz_r_cpu, wbuf_rxlylz_r_cpu => self%wbuf_rxlylz_r_cpu ,&
    & wbuf_rxlyrz_r_cpu => self%wbuf_rxlyrz_r_cpu, wbuf_rxrylz_r_cpu => self%wbuf_rxrylz_r_cpu ,&
    & wbuf_rxryrz_r_cpu => self%wbuf_rxryrz_r_cpu)

      if(steps_(1)) then
        call bcswap_corner_step_1_subroutine(nx, ny, nz, nv, ng, w_cpu, wbuf1s_c_cpu, wbuf2s_c_cpu, wbuf3s_c_cpu, wbuf4s_c_cpu)
        if(pbc(2)) then
          call bcswap_corner_step_1b_subroutine(nx, ny, nz, nv, ng, w_cpu, wbuf_lxly_s_cpu,&
          & wbuf_lxry_s_cpu, wbuf_rxly_s_cpu, wbuf_rxry_s_cpu, wbuf_lylz_s_cpu, wbuf_lyrz_s_cpu,&
          & wbuf_rylz_s_cpu, wbuf_ryrz_s_cpu, wbuf_lxlylz_s_cpu, wbuf_lxlyrz_s_cpu, wbuf_lxrylz_s_cpu,&
          & wbuf_lxryrz_s_cpu, wbuf_rxlylz_s_cpu, wbuf_rxlyrz_s_cpu, wbuf_rxrylz_s_cpu, wbuf_rxryrz_s_cpu)
        endif
      endif
      if(steps_(2)) then
        indc = nv*ny*ng*ng
        if (ileftbottom == nrank) then
          wbuf2r_c_cpu = wbuf1s_c_cpu
          wbuf1r_c_cpu = wbuf2s_c_cpu
        else
          call mpi_sendrecv(wbuf1s_c_cpu,indc,mpi_prec,ileftbottom ,1, wbuf2r_c_cpu,indc,mpi_prec,&
          &irighttop ,1,mp_cart,istatus,iermpi)
          call mpi_sendrecv(wbuf2s_c_cpu,indc,mpi_prec,irighttop ,2, wbuf1r_c_cpu,indc,mpi_prec,&
          &ileftbottom ,2,mp_cart,istatus,iermpi)
        endif
        if (ilefttop == nrank) then
          wbuf4r_c_cpu = wbuf3s_c_cpu
          wbuf3r_c_cpu = wbuf4s_c_cpu
        else
          call mpi_sendrecv(wbuf3s_c_cpu,indc,mpi_prec,ilefttop ,3, wbuf4r_c_cpu,indc,mpi_prec,&
          &irightbottom,3,mp_cart,istatus,iermpi)
          call mpi_sendrecv(wbuf4s_c_cpu,indc,mpi_prec,irightbottom,4, wbuf3r_c_cpu,indc,mpi_prec,&
          &ilefttop ,4,mp_cart,istatus,iermpi)
        endif
        if(pbc(2)) then
          indc = nv*ng*ng*nz
          if (lxly == nrank) then
            wbuf_rxry_r_cpu = wbuf_lxly_s_cpu
            wbuf_lxly_r_cpu = wbuf_rxry_s_cpu
          else
            call mpi_sendrecv(wbuf_lxly_s_cpu,indc,mpi_prec,lxly ,5, wbuf_rxry_r_cpu,indc,mpi_prec,rxry ,5,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rxry_s_cpu,indc,mpi_prec,rxry ,6, wbuf_lxly_r_cpu,indc,mpi_prec,lxly ,6,mp_cart,istatus,iermpi)
          endif
          if (lxry == nrank) then
            wbuf_rxly_r_cpu = wbuf_lxry_s_cpu
            wbuf_lxry_r_cpu = wbuf_rxly_s_cpu
          else
            call mpi_sendrecv(wbuf_lxry_s_cpu,indc,mpi_prec,lxry ,7, wbuf_rxly_r_cpu,indc,mpi_prec,rxly ,7,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rxly_s_cpu,indc,mpi_prec,rxly ,8, wbuf_lxry_r_cpu,indc,mpi_prec,lxry ,8,mp_cart,istatus,iermpi)
          endif
          indc = nv*nx*ng*ng
          if (lylz == nrank) then
            wbuf_ryrz_r_cpu = wbuf_lylz_s_cpu
            wbuf_lylz_r_cpu = wbuf_ryrz_s_cpu
          else
            call mpi_sendrecv(wbuf_lylz_s_cpu,indc,mpi_prec,lylz ,9, wbuf_ryrz_r_cpu,indc,mpi_prec,ryrz ,9,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_ryrz_s_cpu,indc,mpi_prec,ryrz ,10, wbuf_lylz_r_cpu,indc,mpi_prec,&
            &lylz ,10,mp_cart,istatus,iermpi)
          endif
          if (lyrz == nrank) then
            wbuf_rylz_r_cpu = wbuf_lyrz_s_cpu
            wbuf_lyrz_r_cpu = wbuf_rylz_s_cpu
          else
            call mpi_sendrecv(wbuf_lyrz_s_cpu,indc,mpi_prec,lyrz ,11, wbuf_rylz_r_cpu,indc,mpi_prec,&
            &rylz ,11,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rylz_s_cpu,indc,mpi_prec,rylz ,12, wbuf_lyrz_r_cpu,indc,mpi_prec,&
            &lyrz ,12,mp_cart,istatus,iermpi)
          endif
          indc = nv*ng*ng*ng
          if (lxlylz == nrank) then
            wbuf_rxryrz_r_cpu = wbuf_lxlylz_s_cpu
            wbuf_rxrylz_r_cpu = wbuf_lxlyrz_s_cpu
            wbuf_rxlyrz_r_cpu = wbuf_lxrylz_s_cpu
            wbuf_rxlylz_r_cpu = wbuf_lxryrz_s_cpu
            wbuf_lxryrz_r_cpu = wbuf_rxlylz_s_cpu
            wbuf_lxrylz_r_cpu = wbuf_rxlyrz_s_cpu
            wbuf_lxlyrz_r_cpu = wbuf_rxrylz_s_cpu
            wbuf_lxlylz_r_cpu = wbuf_rxryrz_s_cpu
          else
            call mpi_sendrecv(wbuf_lxlylz_s_cpu,indc,mpi_prec,lxlylz ,13, wbuf_rxryrz_r_cpu,indc,&
            &mpi_prec,rxryrz ,13,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rxryrz_s_cpu,indc,mpi_prec,rxryrz ,14, wbuf_lxlylz_r_cpu,indc,&
            &mpi_prec,lxlylz ,14,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_lxlyrz_s_cpu,indc,mpi_prec,lxlyrz ,15, wbuf_rxrylz_r_cpu,indc,&
            &mpi_prec,rxrylz ,15,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rxrylz_s_cpu,indc,mpi_prec,rxrylz ,16, wbuf_lxlyrz_r_cpu,indc,&
            &mpi_prec,lxlyrz ,16,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rxlylz_s_cpu,indc,mpi_prec,rxlylz ,17, wbuf_lxryrz_r_cpu,indc,&
            &mpi_prec,lxryrz ,17,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_lxryrz_s_cpu,indc,mpi_prec,lxryrz ,18, wbuf_rxlylz_r_cpu,indc,&
            &mpi_prec,rxlylz ,18,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rxlyrz_s_cpu,indc,mpi_prec,rxlyrz ,19, wbuf_lxrylz_r_cpu,indc,&
            &mpi_prec,lxrylz ,19,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_lxrylz_s_cpu,indc,mpi_prec,lxrylz ,20, wbuf_rxlyrz_r_cpu,indc,&
            &mpi_prec,rxlyrz ,20,mp_cart,istatus,iermpi)
          endif
        endif
      endif
      if(steps_(3)) then
        call bcswap_corner_step_3_subroutine(nx, ny, nz, nv, ng, w_cpu, wbuf1r_c_cpu, wbuf2r_c_cpu,&
        & wbuf3r_c_cpu, wbuf4r_c_cpu, ileftbottom, ilefttop, irightbottom, irighttop)
        if(pbc(2)) then
          call bcswap_corner_step_3b_subroutine(nx, ny, nz, nv, ng, w_cpu, wbuf_lxly_r_cpu,&
          & wbuf_lxry_r_cpu, wbuf_rxly_r_cpu, wbuf_rxry_r_cpu, wbuf_lylz_r_cpu, wbuf_lyrz_r_cpu,&
          & wbuf_rylz_r_cpu, wbuf_ryrz_r_cpu, wbuf_lxlylz_r_cpu, wbuf_lxlyrz_r_cpu, wbuf_lxrylz_r_cpu,&
          & wbuf_lxryrz_r_cpu, wbuf_rxlylz_r_cpu, wbuf_rxlyrz_r_cpu, wbuf_rxrylz_r_cpu, wbuf_rxryrz_r_cpu,&
          & lxly, lxry, rxly, rxry, lylz, lyrz, rylz, ryrz, lxlylz, lxlyrz, lxrylz, lxryrz, rxlylz, rxlyrz,&
          & rxrylz, rxryrz)
        endif
      endif
    endassociate
  endsubroutine bcswap_edges_corners

  subroutine bcswap_edges_corners_var(self, w_swap, steps)
    class(base_cpu_object), intent(inout) :: self
    real(rkind), dimension(1-self%ng:self%nx+self%ng, 1-self%ng:self%ny+self%ng, 1-self%ng:self%nz+self%ng,1:1) :: w_swap
    logical, dimension(3), optional :: steps
    logical, dimension(3) :: steps_
    integer :: i,j,k,m, iercuda, indc
    integer, dimension(mpi_status_size) :: istatus

    steps_ = .true. ; if(present(steps)) steps_ = steps

    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv,&
    & wbuf1s_c_cpu => self%wbuf1s_c_cpu, wbuf2s_c_cpu => self%wbuf2s_c_cpu, wbuf3s_c_cpu => self%wbuf3s_c&
    &_cpu, wbuf4s_c_cpu => self%wbuf4s_c_cpu, wbuf1r_c_cpu => self%wbuf1r_c_cpu, wbuf2r_c_cpu => self%wbu&
    &f2r_c_cpu, wbuf3r_c_cpu => self%wbuf3r_c_cpu, wbuf4r_c_cpu => self%wbuf4r_c_cpu, ileftbottom => self&
    &%field%ileftbottom, irightbottom => self%field%irightbottom, ilefttop => self%field%ilefttop,&
    & irighttop => self%field%irighttop, pbc => self%field%grid%is_xyz_periodic, mp_cart => self%field%mp&
    &_cart, iermpi => self%mpi_err, nrank => self%field%myrank, lxly => self%field%lxly,&
    & lxry => self%field%lxry, rxly => self%field%rxly, rxry => self%field%rxry,lxlz => self%field%lxlz,&
    & lxrz => self%field%lxrz, rxlz => self%field%rxlz, rxrz => self%field%rxrz,lylz => self%field%lylz,&
    & lyrz => self%field%lyrz, rylz => self%field%rylz, ryrz => self%field%ryrz,lxlylz => self%field%lxly&
    &lz, lxlyrz => self%field%lxlyrz, lxrylz => self%field%lxrylz, rxlylz => self%field%rxlylz,&
    &lxryrz => self%field%lxryrz, rxlyrz => self%field%rxlyrz, rxrylz => self%field%rxrylz,&
    & rxryrz => self%field%rxryrz,wbuf_lxly_s_cpu => self%wbuf_lxly_s_cpu , wbuf_lxry_s_cpu => self%wbuf_&
    &lxry_s_cpu, wbuf_rxly_s_cpu => self%wbuf_rxly_s_cpu , wbuf_rxry_s_cpu => self%wbuf_rxry_s_cpu,&
    & wbuf_lxlz_s_cpu => self%wbuf_lxlz_s_cpu , wbuf_lxrz_s_cpu => self%wbuf_lxrz_s_cpu,&
    & wbuf_rxlz_s_cpu => self%wbuf_rxlz_s_cpu , wbuf_rxrz_s_cpu => self%wbuf_rxrz_s_cpu,&
    & wbuf_lylz_s_cpu => self%wbuf_lylz_s_cpu , wbuf_lyrz_s_cpu => self%wbuf_lyrz_s_cpu,&
    & wbuf_rylz_s_cpu => self%wbuf_rylz_s_cpu , wbuf_ryrz_s_cpu => self%wbuf_ryrz_s_cpu,&
    & wbuf_lxlylz_s_cpu => self%wbuf_lxlylz_s_cpu , wbuf_lxlyrz_s_cpu => self%wbuf_lxlyrz_s_cpu,&
    & wbuf_lxrylz_s_cpu => self%wbuf_lxrylz_s_cpu , wbuf_lxryrz_s_cpu => self%wbuf_lxryrz_s_cpu,&
    & wbuf_rxlylz_s_cpu => self%wbuf_rxlylz_s_cpu , wbuf_rxlyrz_s_cpu => self%wbuf_rxlyrz_s_cpu,&
    & wbuf_rxrylz_s_cpu => self%wbuf_rxrylz_s_cpu , wbuf_rxryrz_s_cpu => self%wbuf_rxryrz_s_cpu,&
    & wbuf_lxly_r_cpu => self%wbuf_lxly_r_cpu , wbuf_lxry_r_cpu => self%wbuf_lxry_r_cpu,&
    & wbuf_rxly_r_cpu => self%wbuf_rxly_r_cpu , wbuf_rxry_r_cpu => self%wbuf_rxry_r_cpu,&
    & wbuf_lxlz_r_cpu => self%wbuf_lxlz_r_cpu , wbuf_lxrz_r_cpu => self%wbuf_lxrz_r_cpu,&
    & wbuf_rxlz_r_cpu => self%wbuf_rxlz_r_cpu , wbuf_rxrz_r_cpu => self%wbuf_rxrz_r_cpu,&
    & wbuf_lylz_r_cpu => self%wbuf_lylz_r_cpu , wbuf_lyrz_r_cpu => self%wbuf_lyrz_r_cpu,&
    & wbuf_rylz_r_cpu => self%wbuf_rylz_r_cpu , wbuf_ryrz_r_cpu => self%wbuf_ryrz_r_cpu,&
    & wbuf_lxlylz_r_cpu => self%wbuf_lxlylz_r_cpu , wbuf_lxlyrz_r_cpu => self%wbuf_lxlyrz_r_cpu,&
    & wbuf_lxrylz_r_cpu => self%wbuf_lxrylz_r_cpu , wbuf_lxryrz_r_cpu => self%wbuf_lxryrz_r_cpu,&
    & wbuf_rxlylz_r_cpu => self%wbuf_rxlylz_r_cpu , wbuf_rxlyrz_r_cpu => self%wbuf_rxlyrz_r_cpu,&
    & wbuf_rxrylz_r_cpu => self%wbuf_rxrylz_r_cpu , wbuf_rxryrz_r_cpu => self%wbuf_rxryrz_r_cpu)

      if(steps_(1)) then
        call bcswap_corner_step_1_subroutine(nx, ny, nz, 1, ng, w_swap, wbuf1s_c_cpu, wbuf2s_c_cpu, wbuf3s_c_cpu, wbuf4s_c_cpu)
        if(pbc(2)) then
          call bcswap_corner_step_1b_subroutine(nx, ny, nz, 1, ng, w_swap, wbuf_lxly_s_cpu,&
          & wbuf_lxry_s_cpu, wbuf_rxly_s_cpu, wbuf_rxry_s_cpu, wbuf_lylz_s_cpu, wbuf_lyrz_s_cpu,&
          & wbuf_rylz_s_cpu, wbuf_ryrz_s_cpu, wbuf_lxlylz_s_cpu, wbuf_lxlyrz_s_cpu, wbuf_lxrylz_s_cpu,&
          & wbuf_lxryrz_s_cpu, wbuf_rxlylz_s_cpu, wbuf_rxlyrz_s_cpu, wbuf_rxrylz_s_cpu, wbuf_rxryrz_s_cpu)
        endif
      endif
      if(steps_(2)) then
        indc = 1*ny*ng*ng
        if (ileftbottom == nrank) then
          wbuf2r_c_cpu = wbuf1s_c_cpu
          wbuf1r_c_cpu = wbuf2s_c_cpu
        else
          call mpi_sendrecv(wbuf1s_c_cpu,indc,mpi_prec,ileftbottom ,1, wbuf2r_c_cpu,indc,mpi_prec,&
          &irighttop ,1,mp_cart,istatus,iermpi)
          call mpi_sendrecv(wbuf2s_c_cpu,indc,mpi_prec,irighttop ,2, wbuf1r_c_cpu,indc,mpi_prec,&
          &ileftbottom ,2,mp_cart,istatus,iermpi)
        endif
        if (ilefttop == nrank) then
          wbuf4r_c_cpu = wbuf3s_c_cpu
          wbuf3r_c_cpu = wbuf4s_c_cpu
        else
          call mpi_sendrecv(wbuf3s_c_cpu,indc,mpi_prec,ilefttop ,3, wbuf4r_c_cpu,indc,mpi_prec,&
          &irightbottom,3,mp_cart,istatus,iermpi)
          call mpi_sendrecv(wbuf4s_c_cpu,indc,mpi_prec,irightbottom,4, wbuf3r_c_cpu,indc,mpi_prec,&
          &ilefttop ,4,mp_cart,istatus,iermpi)
        endif
        if(pbc(2)) then
          indc = 1*ng*ng*nz
          if (lxly == nrank) then
            wbuf_rxry_r_cpu = wbuf_lxly_s_cpu
            wbuf_lxly_r_cpu = wbuf_rxry_s_cpu
          else
            call mpi_sendrecv(wbuf_lxly_s_cpu,indc,mpi_prec,lxly ,5, wbuf_rxry_r_cpu,indc,mpi_prec,rxry ,5,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rxry_s_cpu,indc,mpi_prec,rxry ,6, wbuf_lxly_r_cpu,indc,mpi_prec,lxly ,6,mp_cart,istatus,iermpi)
          endif
          if (lxry == nrank) then
            wbuf_rxly_r_cpu = wbuf_lxry_s_cpu
            wbuf_lxry_r_cpu = wbuf_rxly_s_cpu
          else
            call mpi_sendrecv(wbuf_lxry_s_cpu,indc,mpi_prec,lxry ,7, wbuf_rxly_r_cpu,indc,mpi_prec,rxly ,7,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rxly_s_cpu,indc,mpi_prec,rxly ,8, wbuf_lxry_r_cpu,indc,mpi_prec,lxry ,8,mp_cart,istatus,iermpi)
          endif
          indc = 1*nx*ng*ng
          if (lylz == nrank) then
            wbuf_ryrz_r_cpu = wbuf_lylz_s_cpu
            wbuf_lylz_r_cpu = wbuf_ryrz_s_cpu
          else
            call mpi_sendrecv(wbuf_lylz_s_cpu,indc,mpi_prec,lylz ,9, wbuf_ryrz_r_cpu,indc,mpi_prec,ryrz ,9,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_ryrz_s_cpu,indc,mpi_prec,ryrz ,10, wbuf_lylz_r_cpu,indc,mpi_prec,&
            &lylz ,10,mp_cart,istatus,iermpi)
          endif
          if (lyrz == nrank) then
            wbuf_rylz_r_cpu = wbuf_lyrz_s_cpu
            wbuf_lyrz_r_cpu = wbuf_rylz_s_cpu
          else
            call mpi_sendrecv(wbuf_lyrz_s_cpu,indc,mpi_prec,lyrz ,11, wbuf_rylz_r_cpu,indc,mpi_prec,&
            &rylz ,11,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rylz_s_cpu,indc,mpi_prec,rylz ,12, wbuf_lyrz_r_cpu,indc,mpi_prec,&
            &lyrz ,12,mp_cart,istatus,iermpi)
          endif
          indc = 1*ng*ng*ng
          if (lxlylz == nrank) then
            wbuf_rxryrz_r_cpu = wbuf_lxlylz_s_cpu
            wbuf_rxrylz_r_cpu = wbuf_lxlyrz_s_cpu
            wbuf_rxlyrz_r_cpu = wbuf_lxrylz_s_cpu
            wbuf_rxlylz_r_cpu = wbuf_lxryrz_s_cpu
            wbuf_lxryrz_r_cpu = wbuf_rxlylz_s_cpu
            wbuf_lxrylz_r_cpu = wbuf_rxlyrz_s_cpu
            wbuf_lxlyrz_r_cpu = wbuf_rxrylz_s_cpu
            wbuf_lxlylz_r_cpu = wbuf_rxryrz_s_cpu
          else
            call mpi_sendrecv(wbuf_lxlylz_s_cpu,indc,mpi_prec,lxlylz ,13, wbuf_rxryrz_r_cpu,indc,&
            &mpi_prec,rxryrz ,13,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rxryrz_s_cpu,indc,mpi_prec,rxryrz ,14, wbuf_lxlylz_r_cpu,indc,&
            &mpi_prec,lxlylz ,14,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_lxlyrz_s_cpu,indc,mpi_prec,lxlyrz ,15, wbuf_rxrylz_r_cpu,indc,&
            &mpi_prec,rxrylz ,15,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rxrylz_s_cpu,indc,mpi_prec,rxrylz ,16, wbuf_lxlyrz_r_cpu,indc,&
            &mpi_prec,lxlyrz ,16,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rxlylz_s_cpu,indc,mpi_prec,rxlylz ,17, wbuf_lxryrz_r_cpu,indc,&
            &mpi_prec,lxryrz ,17,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_lxryrz_s_cpu,indc,mpi_prec,lxryrz ,18, wbuf_rxlylz_r_cpu,indc,&
            &mpi_prec,rxlylz ,18,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rxlyrz_s_cpu,indc,mpi_prec,rxlyrz ,19, wbuf_lxrylz_r_cpu,indc,&
            &mpi_prec,lxrylz ,19,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_lxrylz_s_cpu,indc,mpi_prec,lxrylz ,20, wbuf_rxlyrz_r_cpu,indc,&
            &mpi_prec,rxlyrz ,20,mp_cart,istatus,iermpi)
          endif
        endif
      endif
      if(steps_(3)) then
        call bcswap_corner_step_3_subroutine(nx, ny, nz, 1, ng, w_swap, wbuf1r_c_cpu, wbuf2r_c_cpu,&
        & wbuf3r_c_cpu, wbuf4r_c_cpu, ileftbottom, ilefttop, irightbottom, irighttop)
        if(pbc(2)) then
          call bcswap_corner_step_3b_subroutine(nx, ny, nz, 1, ng, w_swap, wbuf_lxly_r_cpu,&
          & wbuf_lxry_r_cpu, wbuf_rxly_r_cpu, wbuf_rxry_r_cpu, wbuf_lylz_r_cpu, wbuf_lyrz_r_cpu,&
          & wbuf_rylz_r_cpu, wbuf_ryrz_r_cpu, wbuf_lxlylz_r_cpu, wbuf_lxlyrz_r_cpu, wbuf_lxrylz_r_cpu,&
          & wbuf_lxryrz_r_cpu, wbuf_rxlylz_r_cpu, wbuf_rxlyrz_r_cpu, wbuf_rxrylz_r_cpu, wbuf_rxryrz_r_cpu,&
          & lxly, lxry, rxly, rxry, lylz, lyrz, rylz, ryrz, lxlylz, lxlyrz, lxrylz, lxryrz, rxlylz, rxlyrz,&
          & rxrylz, rxryrz)
        endif
      endif
    endassociate
  endsubroutine bcswap_edges_corners_var


  subroutine bcswap_var(self, w_swap, steps)
    class(base_cpu_object), intent(inout) :: self
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

    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, wbuf1s_cpu => self%wbuf1s_&
    &cpu, wbuf2s_cpu => self%wbuf2s_cpu, wbuf3s_cpu => self%wbuf3s_cpu, wbuf4s_cpu => self%wbuf4s_cpu,&
    & wbuf5s_cpu => self%wbuf5s_cpu, wbuf6s_cpu => self%wbuf6s_cpu, wbuf1r_cpu => self%wbuf1r_cpu,&
    & wbuf2r_cpu => self%wbuf2r_cpu, wbuf3r_cpu => self%wbuf3r_cpu, wbuf4r_cpu => self%wbuf4r_cpu,&
    & wbuf5r_cpu => self%wbuf5r_cpu, wbuf6r_cpu => self%wbuf6r_cpu, ileftx => self%field%ileftx,&
    &irightx => self%field%irightx, nrank_x => self%field%nrank_x, ilefty => self%field%ilefty,&
    &irighty => self%field%irighty, nrank_y => self%field%nrank_y, ileftz => self%field%ileftz,&
    &irightz => self%field%irightz, nrank_z => self%field%nrank_z, dcsidxnc2_cpu => self%dcsidxnc2_cpu,&
    & dcsidync2_cpu => self%dcsidync2_cpu, detadxnc2_cpu => self%detadxnc2_cpu, detadync2_cpu => self%det&
    &adync2_cpu, dxdcsinc2_cpu => self%dxdcsinc2_cpu, dydcsinc2_cpu => self%dydcsinc2_cpu,&
    & dxdetanc2_cpu => self%dxdetanc2_cpu, dydetanc2_cpu => self%dydetanc2_cpu, mp_cartx => self%field%mp&
    &_cartx, mp_carty => self%field%mp_carty, mp_cartz => self%field%mp_cartz, ndim => self%field%grid%nd&
    &im)

      if(steps_(1)) then
        if(self%field%grid%grid_dim == 1) then
          call bcswap_step_1_subroutine(nx, ny, nz, 1, ng, ndim, w_swap, wbuf1s_cpu, wbuf2s_cpu,&
          & wbuf3s_cpu, wbuf4s_cpu, wbuf5s_cpu, wbuf6s_cpu)
        elseif(self%field%grid%grid_dim == 2) then
          call bcswap_step_1_subroutine(nx, ny, nz, 1, ng, ndim, w_swap, wbuf1s_cpu, wbuf2s_cpu,&
          & wbuf3s_cpu, wbuf4s_cpu, wbuf5s_cpu, wbuf6s_cpu)
        endif
      endif
      if(steps_(2)) then
        indx = ng*ny*nz
        indy = nx*ng*nz
        indz = nx*ny*ng
        if(ileftx == nrank_x) then
          wbuf2r_cpu = wbuf1s_cpu
          wbuf1r_cpu = wbuf2s_cpu
        else
          call mpi_sendrecv(wbuf1s_cpu,indx,mpi_prec,ileftx ,1,wbuf2r_cpu,indx,mpi_prec,irightx,1,mp_cartx,istatus,self%mpi_err)
          call mpi_sendrecv(wbuf2s_cpu,indx,mpi_prec,irightx,2,wbuf1r_cpu,indx,mpi_prec,ileftx ,2,mp_cartx,istatus,self%mpi_err)
        endif
        if(ilefty == nrank_y) then
          wbuf4r_cpu = wbuf3s_cpu
          wbuf3r_cpu = wbuf4s_cpu
        else
          call mpi_sendrecv(wbuf3s_cpu,indy,mpi_prec,ilefty ,3,wbuf4r_cpu,indy,mpi_prec,irighty,3,mp_carty,istatus,self%mpi_err)
          call mpi_sendrecv(wbuf4s_cpu,indy,mpi_prec,irighty,4,wbuf3r_cpu,indy,mpi_prec,ilefty ,4,mp_carty,istatus,self%mpi_err)
        endif
        if(ndim == 3) then
          if(ileftz == nrank_z) then
            wbuf6r_cpu = wbuf5s_cpu
            wbuf5r_cpu = wbuf6s_cpu
          else
            call mpi_sendrecv(wbuf5s_cpu,indz,mpi_prec,ileftz ,5,wbuf6r_cpu,indz,mpi_prec,irightz,5,mp_cartz,istatus,self%mpi_err)
            call mpi_sendrecv(wbuf6s_cpu,indz,mpi_prec,irightz,6,wbuf5r_cpu,indz,mpi_prec,ileftz ,6,mp_cartz,istatus,self%mpi_err)
          endif
        endif
      endif
      if(steps_(3)) then
        if(self%field%grid%grid_dim == 1) then
          call bcswap_step_3_subroutine(nx, ny, nz, 1, ng, ndim, w_swap, wbuf1r_cpu, wbuf2r_cpu,&
          & wbuf3r_cpu, wbuf4r_cpu, wbuf5r_cpu, wbuf6r_cpu, ileftx, ilefty, ileftz, irightx, irighty, irightz)
        elseif(self%field%grid%grid_dim == 2) then
          call bcswap_step_3_subroutine(nx, ny, nz, 1, ng, ndim, w_swap, wbuf1r_cpu, wbuf2r_cpu,&
          & wbuf3r_cpu, wbuf4r_cpu, wbuf5r_cpu, wbuf6r_cpu, ileftx, ilefty, ileftz, irightx, irighty, irightz)
        endif
      endif
    endassociate
  endsubroutine bcswap_var

  subroutine copy_from_field(self)
    class(base_cpu_object), intent(inout) :: self
    integer :: i,j,k,iv

    self%w_cpu = self%field%w

  endsubroutine copy_from_field


  subroutine copy_to_field(self)
    class(base_cpu_object), intent(inout) :: self
    integer :: i,j,k,iv

    self%field%w = self%w_cpu

  endsubroutine copy_to_field

  subroutine initialize(self, field, gpu_bind)
    !< Initialize base backend.
    class(base_cpu_object), intent(inout) :: self !< The base backend.
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


endmodule streams_base_cpu_object


