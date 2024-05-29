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
    real(rkind), allocatable, dimension(:,:,:,:) :: w_cpu
    real(rkind), allocatable, dimension(:,:,:,:) :: w_t

    real(rkind), dimension(:,:,:,:), allocatable :: wbuf1s_cpu, wbuf2s_cpu, wbuf3s_cpu, wbuf4s_cpu, wbuf5s_cpu, wbuf6s_cpu
    real(rkind), dimension(:,:,:,:), allocatable :: wbuf1r_cpu, wbuf2r_cpu, wbuf3r_cpu, wbuf4r_cpu, wbuf5r_cpu, wbuf6r_cpu
    real(rkind), dimension(:,:,:,:), allocatable :: wbuf1s_c_cpu, wbuf2s_c_cpu, wbuf3s_c_cpu, wbuf4s_c_cpu
    real(rkind), dimension(:,:,:,:), allocatable :: wbuf1r_c_cpu, wbuf2r_c_cpu, wbuf3r_c_cpu, wbuf4r_c_cpu
    real(rkind), allocatable, dimension(:) :: dcsidx_cpu, dcsidx2_cpu, dcsidxs_cpu
    real(rkind), allocatable, dimension(:) :: detady_cpu, detady2_cpu, detadys_cpu
    real(rkind), allocatable, dimension(:) :: dzitdz_cpu, dzitdz2_cpu, dzitdzs_cpu
    real(rkind), allocatable, dimension(:) :: x_cpu, y_cpu, z_cpu, yn_cpu

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
    procedure, pass(self) :: bcswap_edges_corners
    procedure, pass(self) :: bcswap_edges_corners_var
  endtype base_cpu_object


contains
  subroutine alloc(self)
    class(base_cpu_object), intent(inout) :: self

    associate(nx => self%field%nx, ny => self%field%ny, nz => self%field%nz, ng => self%field%grid%ng, nv => self%field%nv)

      allocate(self%w_cpu(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng, nv))
      allocate(self%w_t(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng, nv))

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

  endsubroutine alloc

  subroutine bcswap_step_1_subroutine(nx,ny,nz,nv,ng,w_cpu,wbuf1s_cpu,wbuf2s_cpu,wbuf3s_cpu,wbuf4s_cpu,wbuf5s_cpu,wbuf6s_cpu)
    integer :: nx, ny, nz, ng, nv
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


  endsubroutine bcswap_step_1_subroutine


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


  subroutine bcswap_step_3_subroutine(nx,ny,nz,nv,ng,w_cpu,wbuf1r_cpu,wbuf2r_cpu,wbuf3r_cpu,&
  &wbuf4r_cpu,wbuf5r_cpu,wbuf6r_cpu,ileftx,ilefty,ileftz,irightx,irighty,irightz)
    integer :: nx, ny, nz, ng, nv
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

  endsubroutine bcswap_step_3_subroutine


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

  subroutine bcswap(self, steps)
    class(base_cpu_object), intent(inout) :: self
    logical, dimension(3), optional :: steps
    logical, dimension(3) :: steps_
    integer :: iercuda, indx, indy, indz
    integer, dimension(mpi_status_size) :: istatus

    steps_ = .true. ; if(present(steps)) steps_ = steps

    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv,&
    & w_cpu => self%w_cpu, wbuf1s_cpu => self%wbuf1s_cpu, wbuf2s_cpu => self%wbuf2s_cpu,&
    & wbuf3s_cpu => self%wbuf3s_cpu, wbuf4s_cpu => self%wbuf4s_cpu, wbuf5s_cpu => self%wbuf5s_cpu,&
    & wbuf6s_cpu => self%wbuf6s_cpu, wbuf1r_cpu => self%wbuf1r_cpu, wbuf2r_cpu => self%wbuf2r_cpu,&
    & wbuf3r_cpu => self%wbuf3r_cpu, wbuf4r_cpu => self%wbuf4r_cpu, wbuf5r_cpu => self%wbuf5r_cpu,&
    & wbuf6r_cpu => self%wbuf6r_cpu, ileftx => self%field%ileftx,irightx => self%field%irightx,&
    & nrank_x => self%field%nrank_x, ilefty => self%field%ilefty,irighty => self%field%irighty,&
    & nrank_y => self%field%nrank_y, ileftz => self%field%ileftz,irightz => self%field%irightz,&
    & nrank_z => self%field%nrank_z, mp_cartx => self%field%mp_cartx, mp_carty => self%field%mp_carty,&
    & mp_cartz => self%field%mp_cartz)

      if(steps_(1)) then
        call bcswap_step_1_subroutine(nx, ny, nz, nv, ng, w_cpu, wbuf1s_cpu, wbuf2s_cpu, wbuf3s_cpu,&
        & wbuf4s_cpu, wbuf5s_cpu, wbuf6s_cpu)
      endif
      if(steps_(2)) then
        indx = nv*ng*ny*nz
        indy = nv*nx*ng*nz
        indz = nv*nx*ny*ng
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
        if(ileftz == nrank_z) then
          wbuf6r_cpu = wbuf5s_cpu
          wbuf5r_cpu = wbuf6s_cpu
        else
          call mpi_sendrecv(wbuf5s_cpu,indz,mpi_prec,ileftz ,5,wbuf6r_cpu,indz,mpi_prec,irightz,5,mp_cartz,istatus,self%mpi_err)
          call mpi_sendrecv(wbuf6s_cpu,indz,mpi_prec,irightz,6,wbuf5r_cpu,indz,mpi_prec,ileftz ,6,mp_cartz,istatus,self%mpi_err)
        endif
      endif
      if(steps_(3)) then
        call bcswap_step_3_subroutine(nx, ny, nz, nv, ng, w_cpu, wbuf1r_cpu, wbuf2r_cpu, wbuf3r_cpu,&
        & wbuf4r_cpu, wbuf5r_cpu, wbuf6r_cpu, ileftx, ilefty, ileftz, irightx, irighty, irightz)
      endif
    endassociate
  endsubroutine bcswap

  subroutine bcswap_var(self, w_swap, steps)
    class(base_cpu_object), intent(inout) :: self
    real(rkind), dimension(1-self%ng:self%nx+self%ng, 1-self%ng:self%ny+self%ng, 1-self%ng:self%nz+self%ng,1:1) :: w_swap
    logical, dimension(3), optional :: steps
    logical, dimension(3) :: steps_
    integer :: iercuda, indx, indy, indz
    integer, dimension(mpi_status_size) :: istatus

    steps_ = .true. ; if(present(steps)) steps_ = steps

    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, wbuf1s_cpu => self%wbuf1s_&
    &cpu, wbuf2s_cpu => self%wbuf2s_cpu, wbuf3s_cpu => self%wbuf3s_cpu, wbuf4s_cpu => self%wbuf4s_cpu,&
    & wbuf5s_cpu => self%wbuf5s_cpu, wbuf6s_cpu => self%wbuf6s_cpu, wbuf1r_cpu => self%wbuf1r_cpu,&
    & wbuf2r_cpu => self%wbuf2r_cpu, wbuf3r_cpu => self%wbuf3r_cpu, wbuf4r_cpu => self%wbuf4r_cpu,&
    & wbuf5r_cpu => self%wbuf5r_cpu, wbuf6r_cpu => self%wbuf6r_cpu, ileftx => self%field%ileftx,&
    &irightx => self%field%irightx, nrank_x => self%field%nrank_x, ilefty => self%field%ilefty,&
    &irighty => self%field%irighty, nrank_y => self%field%nrank_y, ileftz => self%field%ileftz,&
    &irightz => self%field%irightz, nrank_z => self%field%nrank_z, mp_cartx => self%field%mp_cartx,&
    & mp_carty => self%field%mp_carty, mp_cartz => self%field%mp_cartz)

      if(steps_(1)) then
        call bcswap_step_1_subroutine(nx, ny, nz, 1, ng, w_swap, wbuf1s_cpu, wbuf2s_cpu, wbuf3s_cpu,&
        & wbuf4s_cpu, wbuf5s_cpu, wbuf6s_cpu)
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
        if(ileftz == nrank_z) then
          wbuf6r_cpu = wbuf5s_cpu
          wbuf5r_cpu = wbuf6s_cpu
        else
          call mpi_sendrecv(wbuf5s_cpu,indz,mpi_prec,ileftz ,5,wbuf6r_cpu,indz,mpi_prec,irightz,5,mp_cartz,istatus,self%mpi_err)
          call mpi_sendrecv(wbuf6s_cpu,indz,mpi_prec,irightz,6,wbuf5r_cpu,indz,mpi_prec,ileftz ,6,mp_cartz,istatus,self%mpi_err)
        endif
      endif
      if(steps_(3)) then
        call bcswap_step_3_subroutine(nx, ny, nz, 1, ng, w_swap, wbuf1r_cpu, wbuf2r_cpu, wbuf3r_cpu,&
        & wbuf4r_cpu, wbuf5r_cpu, wbuf6r_cpu, ileftx, ilefty, ileftz, irightx, irighty, irightz)
      endif
    endassociate
  endsubroutine bcswap_var

  subroutine copy_from_field(self)
    class(base_cpu_object), intent(inout) :: self
    integer :: i,j,k,iv

    associate(nx => self%field%nx, ny => self%field%ny, nz => self%field%nz, ng => self%field%grid%ng, nv => self%field%nv)
      do k=1-ng,nz+ng
        do j=1-ng,ny+ng
          do i=1-ng,nx+ng
            do iv=1,nv
              self%w_t(i,j,k,iv) = self%field%w(iv,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endassociate

    self%w_cpu = self%w_t

  endsubroutine copy_from_field


  subroutine copy_to_field(self)
    class(base_cpu_object), intent(inout) :: self
    integer :: i,j,k,iv

    self%w_t = self%w_cpu

    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv)
      do k=1-ng,nz+ng
        do j=1-ng,ny+ng
          do i=1-ng,nx+ng
            do iv=1,nv
              self%field%w(iv,i,j,k) = self%w_t(i,j,k,iv)
            enddo
          enddo
        enddo
      enddo
    endassociate

  endsubroutine copy_to_field

  subroutine initialize(self, field)
    !< Initialize base backend.
    class(base_cpu_object), intent(inout) :: self !< The base backend.
    class(field_object), target :: field

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


