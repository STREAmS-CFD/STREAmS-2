module streams_base_gpu_object

  use streams_field_object, only : field_object
  use streams_parameters
  use mpi
  use cudafor

  implicit none
  private
  public :: base_gpu_object

  type :: base_gpu_object
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
    real(rkind), allocatable, dimension(:,:,:,:), device :: w_gpu
    real(rkind), allocatable, dimension(:,:,:,:) :: w_t

    real(rkind), dimension(:,:,:,:), device, allocatable :: wbuf1s_gpu, wbuf2s_gpu, wbuf3s_gpu, wbuf4s_gpu, wbuf5s_gpu, wbuf6s_gpu
    real(rkind), dimension(:,:,:,:), device, allocatable :: wbuf1r_gpu, wbuf2r_gpu, wbuf3r_gpu, wbuf4r_gpu, wbuf5r_gpu, wbuf6r_gpu
    real(rkind), dimension(:,:,:,:), allocatable, device :: wbuf1s_c_gpu, wbuf2s_c_gpu, wbuf3s_c_gpu, wbuf4s_c_gpu
    real(rkind), dimension(:,:,:,:), allocatable, device :: wbuf1r_c_gpu, wbuf2r_c_gpu, wbuf3r_c_gpu, wbuf4r_c_gpu
    real(rkind), allocatable, dimension(:), device :: dcsidx_gpu, dcsidx2_gpu, dcsidxs_gpu
    real(rkind), allocatable, dimension(:), device :: detady_gpu, detady2_gpu, detadys_gpu
    real(rkind), allocatable, dimension(:), device :: dzitdz_gpu, dzitdz2_gpu, dzitdzs_gpu
    real(rkind), allocatable, dimension(:), device :: x_gpu, y_gpu, z_gpu, yn_gpu

    real(rkind),allocatable,dimension(:,:,:,:),device :: wbuf_lxly_s_gpu , wbuf_lxry_s_gpu , wbuf_rxly_s_gpu , wbuf_rxry_s_gpu
    real(rkind),allocatable,dimension(:,:,:,:),device :: wbuf_lxlz_s_gpu , wbuf_lxrz_s_gpu , wbuf_rxlz_s_gpu , wbuf_rxrz_s_gpu
    real(rkind),allocatable,dimension(:,:,:,:),device :: wbuf_lylz_s_gpu , wbuf_lyrz_s_gpu , wbuf_rylz_s_gpu , wbuf_ryrz_s_gpu
    real(rkind),allocatable,dimension(:,:,:,:),device :: wbuf_lxlylz_s_gpu,wbuf_lxlyrz_s_gpu,wbuf_lxrylz_s_gpu,wbuf_lxryrz_s_gpu
    real(rkind),allocatable,dimension(:,:,:,:),device :: wbuf_rxlylz_s_gpu,wbuf_rxlyrz_s_gpu,wbuf_rxrylz_s_gpu,wbuf_rxryrz_s_gpu
    real(rkind),allocatable,dimension(:,:,:,:),device :: wbuf_lxly_r_gpu , wbuf_lxry_r_gpu , wbuf_rxly_r_gpu , wbuf_rxry_r_gpu
    real(rkind),allocatable,dimension(:,:,:,:),device :: wbuf_lxlz_r_gpu , wbuf_lxrz_r_gpu , wbuf_rxlz_r_gpu , wbuf_rxrz_r_gpu
    real(rkind),allocatable,dimension(:,:,:,:),device :: wbuf_lylz_r_gpu , wbuf_lyrz_r_gpu , wbuf_rylz_r_gpu , wbuf_ryrz_r_gpu
    real(rkind),allocatable,dimension(:,:,:,:),device :: wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu
    real(rkind),allocatable,dimension(:,:,:,:),device :: wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu

  contains
    procedure, pass(self) :: alloc
    procedure, pass(self) :: copy_cpu_gpu
    procedure, pass(self) :: copy_gpu_cpu
    procedure, pass(self) :: initialize
    procedure, pass(self) :: bcswap
    procedure, pass(self) :: bcswap_var
    procedure, pass(self) :: bcswap_corner
    procedure, pass(self) :: bcswap_corner_var
    procedure, pass(self) :: check_gpu_mem
    procedure, pass(self) :: bcswap_edges_corners
    procedure, pass(self) :: bcswap_edges_corners_var
  endtype base_gpu_object


contains
  subroutine alloc(self)
    class(base_gpu_object), intent(inout) :: self

    associate(nx => self%field%nx, ny => self%field%ny, nz => self%field%nz, ng => self%field%grid%ng, nv => self%field%nv)

      allocate(self%w_gpu(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng, nv))
      allocate(self%w_t(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng, nv))

      allocate(self%wbuf1s_gpu(ng,ny,nz,nv))
      allocate(self%wbuf2s_gpu(ng,ny,nz,nv))
      allocate(self%wbuf3s_gpu(nx,ng,nz,nv))
      allocate(self%wbuf4s_gpu(nx,ng,nz,nv))
      allocate(self%wbuf5s_gpu(nx,ny,ng,nv))
      allocate(self%wbuf6s_gpu(nx,ny,ng,nv))
      allocate(self%wbuf1r_gpu(ng,ny,nz,nv))
      allocate(self%wbuf2r_gpu(ng,ny,nz,nv))
      allocate(self%wbuf3r_gpu(nx,ng,nz,nv))
      allocate(self%wbuf4r_gpu(nx,ng,nz,nv))
      allocate(self%wbuf5r_gpu(nx,ny,ng,nv))
      allocate(self%wbuf6r_gpu(nx,ny,ng,nv))

      allocate(self%x_gpu(1-ng:nx+ng), self%y_gpu(1-ng:ny+ng), self%z_gpu(1-ng:nz+ng))
      allocate(self%yn_gpu(1:ny+1))

      allocate(self%dcsidx_gpu(nx), self%dcsidx2_gpu(nx), self%dcsidxs_gpu(nx))
      allocate(self%detady_gpu(ny), self%detady2_gpu(ny), self%detadys_gpu(ny))
      allocate(self%dzitdz_gpu(nz), self%dzitdz2_gpu(nz), self%dzitdzs_gpu(nz))

      allocate(self%wbuf1s_c_gpu(nv,ng,ny,ng))
      allocate(self%wbuf2s_c_gpu(nv,ng,ny,ng))
      allocate(self%wbuf3s_c_gpu(nv,ng,ny,ng))
      allocate(self%wbuf4s_c_gpu(nv,ng,ny,ng))
      allocate(self%wbuf1r_c_gpu(nv,ng,ny,ng))
      allocate(self%wbuf2r_c_gpu(nv,ng,ny,ng))
      allocate(self%wbuf3r_c_gpu(nv,ng,ny,ng))
      allocate(self%wbuf4r_c_gpu(nv,ng,ny,ng))

      allocate(self%wbuf_lxly_s_gpu(nv,ng,ng,nz), self%wbuf_lxry_s_gpu(nv,ng,ng,nz))
      allocate(self%wbuf_rxly_s_gpu(nv,ng,ng,nz), self%wbuf_rxry_s_gpu(nv,ng,ng,nz))
      allocate(self%wbuf_lxlz_s_gpu(nv,ng,ny,ng), self%wbuf_lxrz_s_gpu(nv,ng,ny,ng))
      allocate(self%wbuf_rxlz_s_gpu(nv,ng,ny,ng), self%wbuf_rxrz_s_gpu(nv,ng,ny,ng))
      allocate(self%wbuf_lylz_s_gpu(nv,nx,ng,ng), self%wbuf_lyrz_s_gpu(nv,nx,ng,ng))
      allocate(self%wbuf_rylz_s_gpu(nv,nx,ng,ng), self%wbuf_ryrz_s_gpu(nv,nx,ng,ng))
      allocate(self%wbuf_lxlylz_s_gpu(nv,ng,ng,ng),self%wbuf_lxlyrz_s_gpu(nv,ng,ng,ng))
      allocate(self%wbuf_lxrylz_s_gpu(nv,ng,ng,ng),self%wbuf_lxryrz_s_gpu(nv,ng,ng,ng))
      allocate(self%wbuf_rxlylz_s_gpu(nv,ng,ng,ng),self%wbuf_rxlyrz_s_gpu(nv,ng,ng,ng))
      allocate(self%wbuf_rxrylz_s_gpu(nv,ng,ng,ng),self%wbuf_rxryrz_s_gpu(nv,ng,ng,ng))

      allocate(self%wbuf_lxly_r_gpu(nv,ng,ng,nz), self%wbuf_lxry_r_gpu(nv,ng,ng,nz))
      allocate(self%wbuf_rxly_r_gpu(nv,ng,ng,nz), self%wbuf_rxry_r_gpu(nv,ng,ng,nz))
      allocate(self%wbuf_lxlz_r_gpu(nv,ng,ny,ng), self%wbuf_lxrz_r_gpu(nv,ng,ny,ng))
      allocate(self%wbuf_rxlz_r_gpu(nv,ng,ny,ng), self%wbuf_rxrz_r_gpu(nv,ng,ny,ng))
      allocate(self%wbuf_lylz_r_gpu(nv,nx,ng,ng), self%wbuf_lyrz_r_gpu(nv,nx,ng,ng))
      allocate(self%wbuf_rylz_r_gpu(nv,nx,ng,ng), self%wbuf_ryrz_r_gpu(nv,nx,ng,ng))
      allocate(self%wbuf_lxlylz_r_gpu(nv,ng,ng,ng),self%wbuf_lxlyrz_r_gpu(nv,ng,ng,ng))
      allocate(self%wbuf_lxrylz_r_gpu(nv,ng,ng,ng),self%wbuf_lxryrz_r_gpu(nv,ng,ng,ng))
      allocate(self%wbuf_rxlylz_r_gpu(nv,ng,ng,ng),self%wbuf_rxlyrz_r_gpu(nv,ng,ng,ng))
      allocate(self%wbuf_rxrylz_r_gpu(nv,ng,ng,ng),self%wbuf_rxryrz_r_gpu(nv,ng,ng,ng))

    endassociate

    self%x_gpu = self%field%x
    self%y_gpu = self%field%y
    self%z_gpu = self%field%z
    self%yn_gpu = self%field%yn

    self%dcsidx_gpu = self%field%dcsidx
    self%dcsidxs_gpu = self%field%dcsidxs
    self%dcsidx2_gpu = self%field%dcsidx2
    self%detady_gpu = self%field%detady
    self%detadys_gpu = self%field%detadys
    self%detady2_gpu = self%field%detady2
    self%dzitdz_gpu = self%field%dzitdz
    self%dzitdzs_gpu = self%field%dzitdzs
    self%dzitdz2_gpu = self%field%dzitdz2

  endsubroutine alloc

  subroutine bcswap_step_1_cuf(nx, ny, nz, nv, ng, w_gpu, wbuf1s_gpu, wbuf2s_gpu, wbuf3s_gpu, wbuf4s_gpu, wbuf5s_gpu, wbuf6s_gpu)
    integer :: nx, ny, nz, ng, nv
    real(rkind), device, dimension(1-ng:,1-ng:,1-ng:,1:) :: w_gpu
    real(rkind), device, dimension(1:,1:,1:,1:) :: wbuf1s_gpu, wbuf2s_gpu, wbuf3s_gpu, wbuf4s_gpu, wbuf5s_gpu, wbuf6s_gpu
    integer :: i,j,k,m,iercuda

    !$cuf kernel do(3) <<<*,*>>>
    do k=1,nz
      do j=1,ny
        do i=1,ng
          do m=1,nv
            wbuf1s_gpu(i,j,k,m) = w_gpu(i,j,k,m)
            wbuf2s_gpu(i,j,k,m) = w_gpu(nx-ng+i,j,k,m)
          enddo
        enddo
      enddo
    enddo
    !@cuf iercuda=cudadevicesynchronize()


    !$cuf kernel do(3) <<<*,*>>>
    do k=1,nz
      do j=1,ng
        do i=1,nx
          do m=1,nv
            wbuf3s_gpu(i,j,k,m) = w_gpu(i,j,k,m)
            wbuf4s_gpu(i,j,k,m) = w_gpu(i,ny-ng+j,k,m)
          enddo
        enddo
      enddo
    enddo
    !@cuf iercuda=cudadevicesynchronize()

    !$cuf kernel do(3) <<<*,*>>>
    do k=1,ng
      do j=1,ny
        do i=1,nx
          do m=1,nv
            wbuf5s_gpu(i,j,k,m) = w_gpu(i,j,k,m)
            wbuf6s_gpu(i,j,k,m) = w_gpu(i,j,nz-ng+k,m)
          enddo
        enddo
      enddo
    enddo
    !@cuf iercuda=cudadevicesynchronize()

  endsubroutine bcswap_step_1_cuf

  subroutine bcswap_corner_step_1_cuf(nx, ny, nz, nv, ng, w_gpu, wbuf1s_c_gpu, wbuf2s_c_gpu, wbuf3s_c_gpu, wbuf4s_c_gpu)
    integer :: nx, ny, nz, ng, nv
    real(rkind), device, dimension(1-ng:,1-ng:,1-ng:,1:) :: w_gpu
    real(rkind), device, dimension(1:,1:,1:,1:) :: wbuf1s_c_gpu, wbuf2s_c_gpu, wbuf3s_c_gpu, wbuf4s_c_gpu
    integer :: i,j,k,m,iercuda

    !$cuf kernel do(2) <<<*,*>>>
    do j=1,ny
      do m=1,nv
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
    !@cuf iercuda=cudadevicesynchronize()

  endsubroutine bcswap_corner_step_1_cuf

  subroutine bcswap_corner_step_1b_cuf(nx, ny, nz, nv, ng, w_gpu, wbuf_lxly_s_gpu, wbuf_lxry_s_gpu,&
  & wbuf_rxly_s_gpu, wbuf_rxry_s_gpu, wbuf_lylz_s_gpu, wbuf_lyrz_s_gpu, wbuf_rylz_s_gpu,&
  & wbuf_ryrz_s_gpu, wbuf_lxlylz_s_gpu, wbuf_lxlyrz_s_gpu, wbuf_lxrylz_s_gpu, wbuf_lxryrz_s_gpu,&
  & wbuf_rxlylz_s_gpu, wbuf_rxlyrz_s_gpu, wbuf_rxrylz_s_gpu, wbuf_rxryrz_s_gpu)
    integer :: nx, ny, nz, ng, nv
    real(rkind), device, dimension(1-ng:,1-ng:,1-ng:,1:) :: w_gpu
    real(rkind), device, dimension(1:,1:,1:,1:) :: wbuf_lxly_s_gpu, wbuf_lxry_s_gpu, wbuf_rxly_s_gpu, wbuf_rxry_s_gpu
    real(rkind), device, dimension(1:,1:,1:,1:) :: wbuf_lylz_s_gpu, wbuf_lyrz_s_gpu, wbuf_rylz_s_gpu, wbuf_ryrz_s_gpu
    real(rkind), device, dimension(1:,1:,1:,1:) :: wbuf_lxlylz_s_gpu, wbuf_lxlyrz_s_gpu, wbuf_lxrylz_s_gpu, wbuf_lxryrz_s_gpu
    real(rkind), device, dimension(1:,1:,1:,1:) :: wbuf_rxlylz_s_gpu, wbuf_rxlyrz_s_gpu, wbuf_rxrylz_s_gpu, wbuf_rxryrz_s_gpu
    integer :: i,j,k,m,iercuda

    !$cuf kernel do(2) <<<*,*>>>
    do k=1,nz
      do m=1,nv
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
    !@cuf iercuda=cudadevicesynchronize()

    !$cuf kernel do(2) <<<*,*>>>
    do i=1,nx
      do m=1,nv
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
    !@cuf iercuda=cudadevicesynchronize()

    !$cuf kernel do(2) <<<*,*>>>
    do k=1,ng
      do m=1,nv
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
    !@cuf iercuda=cudadevicesynchronize()

  endsubroutine bcswap_corner_step_1b_cuf

  subroutine bcswap_step_3_cuf(nx, ny, nz, nv, ng, w_gpu, wbuf1r_gpu, wbuf2r_gpu, wbuf3r_gpu,&
  & wbuf4r_gpu, wbuf5r_gpu, wbuf6r_gpu, ileftx, ilefty, ileftz, irightx, irighty, irightz)
    integer :: nx, ny, nz, ng, nv
    real(rkind), device, dimension(1-ng:,1-ng:,1-ng:,1:) :: w_gpu
    real(rkind), device, dimension(1:,1:,1:,1:) :: wbuf1r_gpu, wbuf2r_gpu, wbuf3r_gpu, wbuf4r_gpu, wbuf5r_gpu, wbuf6r_gpu
    integer :: i,j,k,m, iercuda
    integer :: ileftx, ilefty, ileftz, irightx, irighty, irightz
    if (ileftx/=mpi_proc_null) then
      !$cuf kernel do(3) <<<*,*>>>
      do k=1,nz
        do j=1,ny
          do i=1,ng
            do m=1,nv
              w_gpu(i-ng,j,k,m) = wbuf1r_gpu(i,j,k,m)
            enddo
          enddo
        enddo
      enddo
      !@cuf iercuda=cudadevicesynchronize()
    endif
    if (irightx/=mpi_proc_null) then
      !$cuf kernel do(3) <<<*,*>>>
      do k=1,nz
        do j=1,ny
          do i=1,ng
            do m=1,nv
              w_gpu(nx+i,j,k,m) = wbuf2r_gpu(i,j,k,m)
            enddo
          enddo
        enddo
      enddo
      !@cuf iercuda=cudadevicesynchronize()
    endif
    if (ilefty/=mpi_proc_null) then
      !$cuf kernel do(3) <<<*,*>>>
      do k=1,nz
        do j=1,ng
          do i=1,nx
            do m=1,nv
              w_gpu(i,j-ng,k,m) = wbuf3r_gpu(i,j,k,m)
            enddo
          enddo
        enddo
      enddo
      !@cuf iercuda=cudadevicesynchronize()
    endif
    if (irighty/=mpi_proc_null) then
      !$cuf kernel do(3) <<<*,*>>>
      do k=1,nz
        do j=1,ng
          do i=1,nx
            do m=1,nv
              w_gpu(i,ny+j,k,m) = wbuf4r_gpu(i,j,k,m)
            enddo
          enddo
        enddo
      enddo
      !@cuf iercuda=cudadevicesynchronize()
    endif
    if (ileftz/=mpi_proc_null) then
      !$cuf kernel do(3) <<<*,*>>>
      do k=1,ng
        do j=1,ny
          do i=1,nx
            do m=1,nv
              w_gpu(i,j,k-ng,m) = wbuf5r_gpu(i,j,k,m)
            enddo
          enddo
        enddo
      enddo
      !@cuf iercuda=cudadevicesynchronize()
    endif
    if (irightz/=mpi_proc_null) then
      !$cuf kernel do(3) <<<*,*>>>
      do k=1,ng
        do j=1,ny
          do i=1,nx
            do m=1,nv
              w_gpu(i,j,nz+k,m) = wbuf6r_gpu(i,j,k,m)
            enddo
          enddo
        enddo
      enddo
      !@cuf iercuda=cudadevicesynchronize()
    endif
  endsubroutine bcswap_step_3_cuf

  subroutine bcswap_corner_step_3_cuf(nx, ny, nz, nv, ng, w_gpu, wbuf1r_c_gpu, wbuf2r_c_gpu,&
  & wbuf3r_c_gpu, wbuf4r_c_gpu, ileftbottom, ilefttop, irightbottom, irighttop)
    integer :: nx, ny, nz, ng, nv
    real(rkind), device, dimension(1-ng:,1-ng:,1-ng:,1:) :: w_gpu
    real(rkind), device, dimension(1:,1:,1:,1:) :: wbuf1r_c_gpu, wbuf2r_c_gpu, wbuf3r_c_gpu, wbuf4r_c_gpu
    integer :: i,j,k,m, iercuda
    integer :: ileftbottom, ilefttop, irightbottom, irighttop

    if (ileftbottom/=mpi_proc_null) then
      !$cuf kernel do(2) <<<*,*>>>
      do j=1,ny
        do m=1,nv
          do k=1,ng
            do i=1,ng
              w_gpu(i-ng,j,k-ng,m) = wbuf1r_c_gpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
      !@cuf iercuda=cudadevicesynchronize()
    endif
    if (irighttop/=mpi_proc_null) then
      !$cuf kernel do(2) <<<*,*>>>
      do j=1,ny
        do m=1,nv
          do k=1,ng
            do i=1,ng
              w_gpu(i+nx,j,k+nz,m) = wbuf2r_c_gpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
      !@cuf iercuda=cudadevicesynchronize()
    endif
    if (ilefttop/=mpi_proc_null) then
      !$cuf kernel do(2) <<<*,*>>>
      do j=1,ny
        do m=1,nv
          do k=1,ng
            do i=1,ng
              w_gpu(i-ng,j,k+nz,m) = wbuf3r_c_gpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
      !@cuf iercuda=cudadevicesynchronize()
    endif
    if (irightbottom/=mpi_proc_null) then
      !$cuf kernel do(2) <<<*,*>>>
      do j=1,ny
        do m=1,nv
          do k=1,ng
            do i=1,ng
              w_gpu(i+nx,j,k-ng,m) = wbuf4r_c_gpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
      !@cuf iercuda=cudadevicesynchronize()
    endif
  endsubroutine bcswap_corner_step_3_cuf

  subroutine bcswap_corner_step_3b_cuf(nx, ny, nz, nv, ng, w_gpu, wbuf_lxly_r_gpu, wbuf_lxry_r_gpu,&
  & wbuf_rxly_r_gpu, wbuf_rxry_r_gpu, wbuf_lylz_r_gpu, wbuf_lyrz_r_gpu, wbuf_rylz_r_gpu,&
  & wbuf_ryrz_r_gpu, wbuf_lxlylz_r_gpu, wbuf_lxlyrz_r_gpu, wbuf_lxrylz_r_gpu, wbuf_lxryrz_r_gpu,&
  & wbuf_rxlylz_r_gpu, wbuf_rxlyrz_r_gpu, wbuf_rxrylz_r_gpu, wbuf_rxryrz_r_gpu, lxly, lxry, rxly, rxry,&
  & lylz, lyrz, rylz, ryrz, lxlylz, lxlyrz, lxrylz, lxryrz, rxlylz, rxlyrz, rxrylz, rxryrz)
    integer :: nx, ny, nz, ng, nv
    real(rkind), device, dimension(1-ng:,1-ng:,1-ng:,1:) :: w_gpu
    real(rkind), device, dimension(1:,1:,1:,1:) :: wbuf_lxly_r_gpu, wbuf_lxry_r_gpu, wbuf_rxly_r_gpu, wbuf_rxry_r_gpu
    real(rkind), device, dimension(1:,1:,1:,1:) :: wbuf_lylz_r_gpu, wbuf_lyrz_r_gpu, wbuf_rylz_r_gpu, wbuf_ryrz_r_gpu
    real(rkind), device, dimension(1:,1:,1:,1:) :: wbuf_lxlylz_r_gpu, wbuf_lxlyrz_r_gpu, wbuf_lxrylz_r_gpu, wbuf_lxryrz_r_gpu
    real(rkind), device, dimension(1:,1:,1:,1:) :: wbuf_rxlylz_r_gpu, wbuf_rxlyrz_r_gpu, wbuf_rxrylz_r_gpu, wbuf_rxryrz_r_gpu
    integer :: i,j,k,m, iercuda
    integer :: lxly, lxry, rxly, rxry
    integer :: lylz, lyrz, rylz, ryrz
    integer :: lxlylz, lxlyrz, lxrylz, lxryrz
    integer :: rxlylz, rxlyrz, rxrylz, rxryrz

    if (lxly/=mpi_proc_null) then
      !$cuf kernel do(2) <<<*,*>>>
      do k=1,nz
        do m=1,nv
          do i=1,ng
            do j=1,ng
              w_gpu(i-ng,j-ng,k,m) = wbuf_lxly_r_gpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
      !@cuf iercuda=cudadevicesynchronize()
    endif
    if (rxry/=mpi_proc_null) then
      !$cuf kernel do(2) <<<*,*>>>
      do k=1,nz
        do m=1,nv
          do i=1,ng
            do j=1,ng
              w_gpu(i+nx,j+ny,k,m) = wbuf_rxry_r_gpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
      !@cuf iercuda=cudadevicesynchronize()
    endif
    if (lxry/=mpi_proc_null) then
      !$cuf kernel do(2) <<<*,*>>>
      do k=1,nz
        do m=1,nv
          do i=1,ng
            do j=1,ng
              w_gpu(i-ng,j+ny,k,m) = wbuf_lxry_r_gpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
      !@cuf iercuda=cudadevicesynchronize()
    endif
    if (rxly/=mpi_proc_null) then
      !$cuf kernel do(2) <<<*,*>>>
      do k=1,nz
        do m=1,nv
          do i=1,ng
            do j=1,ng
              w_gpu(i+nx,j-ng,k,m) = wbuf_rxly_r_gpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
      !@cuf iercuda=cudadevicesynchronize()
    endif

    if (lylz/=mpi_proc_null) then
      !$cuf kernel do(2) <<<*,*>>>
      do i=1,nx
        do m=1,nv
          do k=1,ng
            do j=1,ng
              w_gpu(i,j-ng,k-ng,m) = wbuf_lylz_r_gpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
      !@cuf iercuda=cudadevicesynchronize()
    endif
    if (ryrz/=mpi_proc_null) then
      !$cuf kernel do(2) <<<*,*>>>
      do i=1,nx
        do m=1,nv
          do k=1,ng
            do j=1,ng
              w_gpu(i,j+ny,k+nz,m) = wbuf_ryrz_r_gpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
      !@cuf iercuda=cudadevicesynchronize()
    endif
    if (lyrz/=mpi_proc_null) then
      !$cuf kernel do(2) <<<*,*>>>
      do i=1,nx
        do m=1,nv
          do k=1,ng
            do j=1,ng
              w_gpu(i,j-ng,k+nz,m) = wbuf_lyrz_r_gpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
      !@cuf iercuda=cudadevicesynchronize()
    endif
    if (rylz/=mpi_proc_null) then
      !$cuf kernel do(2) <<<*,*>>>
      do i=1,nx
        do m=1,nv
          do k=1,ng
            do j=1,ng
              w_gpu(i,j+ny,k-ng,m) = wbuf_rylz_r_gpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
      !@cuf iercuda=cudadevicesynchronize()
    endif

    if (lxlylz/=mpi_proc_null) then
      !$cuf kernel do(2) <<<*,*>>>
      do k=1,ng
        do m=1,nv
          do i=1,ng
            do j=1,ng
              w_gpu(i-ng,j-ng,k-ng,m) = wbuf_lxlylz_r_gpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
      !@cuf iercuda=cudadevicesynchronize()
    endif
    if (lxlyrz/=mpi_proc_null) then
      !$cuf kernel do(2) <<<*,*>>>
      do k=1,ng
        do m=1,nv
          do i=1,ng
            do j=1,ng
              w_gpu(i-ng,j-ng,k+nz,m) = wbuf_lxlyrz_r_gpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
      !@cuf iercuda=cudadevicesynchronize()
    endif
    if (lxrylz/=mpi_proc_null) then
      !$cuf kernel do(2) <<<*,*>>>
      do k=1,ng
        do m=1,nv
          do i=1,ng
            do j=1,ng
              w_gpu(i-ng,j+ny,k-ng,m) = wbuf_lxrylz_r_gpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
      !@cuf iercuda=cudadevicesynchronize()
    endif
    if (lxryrz/=mpi_proc_null) then
      !$cuf kernel do(2) <<<*,*>>>
      do k=1,ng
        do m=1,nv
          do i=1,ng
            do j=1,ng
              w_gpu(i-ng,j+ny,k+nz,m) = wbuf_lxryrz_r_gpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
      !@cuf iercuda=cudadevicesynchronize()
    endif
    if (rxlylz/=mpi_proc_null) then
      !$cuf kernel do(2) <<<*,*>>>
      do k=1,ng
        do m=1,nv
          do i=1,ng
            do j=1,ng
              w_gpu(i+nx,j-ng,k-ng,m) = wbuf_rxlylz_r_gpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
      !@cuf iercuda=cudadevicesynchronize()
    endif
    if (rxlyrz/=mpi_proc_null) then
      !$cuf kernel do(2) <<<*,*>>>
      do k=1,ng
        do m=1,nv
          do i=1,ng
            do j=1,ng
              w_gpu(i+nx,j-ng,k+nz,m) = wbuf_rxlyrz_r_gpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
      !@cuf iercuda=cudadevicesynchronize()
    endif
    if (rxrylz/=mpi_proc_null) then
      !$cuf kernel do(2) <<<*,*>>>
      do k=1,ng
        do m=1,nv
          do i=1,ng
            do j=1,ng
              w_gpu(i+nx,j+ny,k-ng,m) = wbuf_rxrylz_r_gpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
      !@cuf iercuda=cudadevicesynchronize()
    endif
    if (rxryrz/=mpi_proc_null) then
      !$cuf kernel do(2) <<<*,*>>>
      do k=1,ng
        do m=1,nv
          do i=1,ng
            do j=1,ng
              w_gpu(i+nx,j+ny,k+nz,m) = wbuf_rxryrz_r_gpu(m,i,j,k)
            enddo
          enddo
        enddo
      enddo
      !@cuf iercuda=cudadevicesynchronize()
    endif
  endsubroutine bcswap_corner_step_3b_cuf

  subroutine bcswap_corner(self, steps)
    class(base_gpu_object), intent(inout) :: self
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
        call bcswap_corner_step_1_cuf(nx, ny, nz, nv, ng, w_gpu, wbuf1s_c_gpu, wbuf2s_c_gpu, wbuf3s_c_gpu, wbuf4s_c_gpu)
      endif
      if(steps_(2)) then
        indc = nv*ny*ng*ng
        if (ileftbottom == nrank) then
          wbuf2r_c_gpu = wbuf1s_c_gpu
          wbuf1r_c_gpu = wbuf2s_c_gpu
        else
          call mpi_sendrecv(wbuf1s_c_gpu,indc,mpi_prec,ileftbottom ,1, wbuf2r_c_gpu,indc,mpi_prec,&
          &irighttop ,1,mp_cart,istatus,iermpi)
          call mpi_sendrecv(wbuf2s_c_gpu,indc,mpi_prec,irighttop ,2, wbuf1r_c_gpu,indc,mpi_prec,&
          &ileftbottom ,2,mp_cart,istatus,iermpi)
        endif
        if (ilefttop == nrank) then
          wbuf4r_c_gpu = wbuf3s_c_gpu
          wbuf3r_c_gpu = wbuf4s_c_gpu
        else
          call mpi_sendrecv(wbuf3s_c_gpu,indc,mpi_prec,ilefttop ,3, wbuf4r_c_gpu,indc,mpi_prec,&
          &irightbottom,3,mp_cart,istatus,iermpi)
          call mpi_sendrecv(wbuf4s_c_gpu,indc,mpi_prec,irightbottom,4, wbuf3r_c_gpu,indc,mpi_prec,&
          &ilefttop ,4,mp_cart,istatus,iermpi)
        endif
      endif
      if(steps_(3)) then
        call bcswap_corner_step_3_cuf(nx, ny, nz, nv, ng, w_gpu, wbuf1r_c_gpu, wbuf2r_c_gpu,&
        & wbuf3r_c_gpu, wbuf4r_c_gpu, ileftbottom, ilefttop, irightbottom, irighttop)
      endif
    endassociate
  endsubroutine bcswap_corner

  subroutine bcswap_corner_var(self, w_swap, steps)
    class(base_gpu_object), intent(inout) :: self
    real(rkind), dimension(1-self%ng:self%nx+self%ng, 1-self%ng:self%ny+self%ng, 1-self%ng:self%nz+self%ng,1:1), device :: w_swap
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
        call bcswap_corner_step_1_cuf(nx, ny, nz, 1, ng, w_swap, wbuf1s_c_gpu, wbuf2s_c_gpu, wbuf3s_c_gpu, wbuf4s_c_gpu)
      endif
      if(steps_(2)) then
        indc = ny*ng*ng
        if (ileftbottom == nrank) then
          wbuf2r_c_gpu = wbuf1s_c_gpu
          wbuf1r_c_gpu = wbuf2s_c_gpu
        else
          call mpi_sendrecv(wbuf1s_c_gpu,indc,mpi_prec,ileftbottom ,1, wbuf2r_c_gpu,indc,mpi_prec,&
          &irighttop ,1,mp_cart,istatus,iermpi)
          call mpi_sendrecv(wbuf2s_c_gpu,indc,mpi_prec,irighttop ,2, wbuf1r_c_gpu,indc,mpi_prec,&
          &ileftbottom ,2,mp_cart,istatus,iermpi)
        endif
        if (ilefttop == nrank) then
          wbuf4r_c_gpu = wbuf3s_c_gpu
          wbuf3r_c_gpu = wbuf4s_c_gpu
        else
          call mpi_sendrecv(wbuf3s_c_gpu,indc,mpi_prec,ilefttop ,3, wbuf4r_c_gpu,indc,mpi_prec,&
          &irightbottom,3,mp_cart,istatus,iermpi)
          call mpi_sendrecv(wbuf4s_c_gpu,indc,mpi_prec,irightbottom,4, wbuf3r_c_gpu,indc,mpi_prec,&
          &ilefttop ,4,mp_cart,istatus,iermpi)
        endif
      endif
      if(steps_(3)) then
        call bcswap_corner_step_3_cuf(nx, ny, nz, 1, ng, w_swap, wbuf1r_c_gpu, wbuf2r_c_gpu,&
        & wbuf3r_c_gpu, wbuf4r_c_gpu, ileftbottom, ilefttop, irightbottom, irighttop)
      endif
    endassociate
  endsubroutine bcswap_corner_var

  subroutine bcswap_edges_corners(self, steps)
    class(base_gpu_object), intent(inout) :: self
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
        call bcswap_corner_step_1_cuf(nx, ny, nz, nv, ng, w_gpu, wbuf1s_c_gpu, wbuf2s_c_gpu, wbuf3s_c_gpu, wbuf4s_c_gpu)
        if(pbc(2)) then
          call bcswap_corner_step_1b_cuf(nx, ny, nz, nv, ng, w_gpu, wbuf_lxly_s_gpu,&
          & wbuf_lxry_s_gpu, wbuf_rxly_s_gpu, wbuf_rxry_s_gpu, wbuf_lylz_s_gpu, wbuf_lyrz_s_gpu,&
          & wbuf_rylz_s_gpu, wbuf_ryrz_s_gpu, wbuf_lxlylz_s_gpu, wbuf_lxlyrz_s_gpu, wbuf_lxrylz_s_gpu,&
          & wbuf_lxryrz_s_gpu, wbuf_rxlylz_s_gpu, wbuf_rxlyrz_s_gpu, wbuf_rxrylz_s_gpu, wbuf_rxryrz_s_gpu)
        endif
      endif
      if(steps_(2)) then
        indc = nv*ny*ng*ng
        if (ileftbottom == nrank) then
          wbuf2r_c_gpu = wbuf1s_c_gpu
          wbuf1r_c_gpu = wbuf2s_c_gpu
        else
          call mpi_sendrecv(wbuf1s_c_gpu,indc,mpi_prec,ileftbottom ,1, wbuf2r_c_gpu,indc,mpi_prec,&
          &irighttop ,1,mp_cart,istatus,iermpi)
          call mpi_sendrecv(wbuf2s_c_gpu,indc,mpi_prec,irighttop ,2, wbuf1r_c_gpu,indc,mpi_prec,&
          &ileftbottom ,2,mp_cart,istatus,iermpi)
        endif
        if (ilefttop == nrank) then
          wbuf4r_c_gpu = wbuf3s_c_gpu
          wbuf3r_c_gpu = wbuf4s_c_gpu
        else
          call mpi_sendrecv(wbuf3s_c_gpu,indc,mpi_prec,ilefttop ,3, wbuf4r_c_gpu,indc,mpi_prec,&
          &irightbottom,3,mp_cart,istatus,iermpi)
          call mpi_sendrecv(wbuf4s_c_gpu,indc,mpi_prec,irightbottom,4, wbuf3r_c_gpu,indc,mpi_prec,&
          &ilefttop ,4,mp_cart,istatus,iermpi)
        endif
        if(pbc(2)) then
          indc = nv*ng*ng*nz
          if (lxly == nrank) then
            wbuf_rxry_r_gpu = wbuf_lxly_s_gpu
            wbuf_lxly_r_gpu = wbuf_rxry_s_gpu
          else
            call mpi_sendrecv(wbuf_lxly_s_gpu,indc,mpi_prec,lxly ,5, wbuf_rxry_r_gpu,indc,mpi_prec,rxry ,5,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rxry_s_gpu,indc,mpi_prec,rxry ,6, wbuf_lxly_r_gpu,indc,mpi_prec,lxly ,6,mp_cart,istatus,iermpi)
          endif
          if (lxry == nrank) then
            wbuf_rxly_r_gpu = wbuf_lxry_s_gpu
            wbuf_lxry_r_gpu = wbuf_rxly_s_gpu
          else
            call mpi_sendrecv(wbuf_lxry_s_gpu,indc,mpi_prec,lxry ,7, wbuf_rxly_r_gpu,indc,mpi_prec,rxly ,7,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rxly_s_gpu,indc,mpi_prec,rxly ,8, wbuf_lxry_r_gpu,indc,mpi_prec,lxry ,8,mp_cart,istatus,iermpi)
          endif
          indc = nv*nx*ng*ng
          if (lylz == nrank) then
            wbuf_ryrz_r_gpu = wbuf_lylz_s_gpu
            wbuf_lylz_r_gpu = wbuf_ryrz_s_gpu
          else
            call mpi_sendrecv(wbuf_lylz_s_gpu,indc,mpi_prec,lylz ,9, wbuf_ryrz_r_gpu,indc,mpi_prec,ryrz ,9,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_ryrz_s_gpu,indc,mpi_prec,ryrz ,10, wbuf_lylz_r_gpu,indc,mpi_prec,&
            &lylz ,10,mp_cart,istatus,iermpi)
          endif
          if (lyrz == nrank) then
            wbuf_rylz_r_gpu = wbuf_lyrz_s_gpu
            wbuf_lyrz_r_gpu = wbuf_rylz_s_gpu
          else
            call mpi_sendrecv(wbuf_lyrz_s_gpu,indc,mpi_prec,lyrz ,11, wbuf_rylz_r_gpu,indc,mpi_prec,&
            &rylz ,11,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rylz_s_gpu,indc,mpi_prec,rylz ,12, wbuf_lyrz_r_gpu,indc,mpi_prec,&
            &lyrz ,12,mp_cart,istatus,iermpi)
          endif
          indc = nv*ng*ng*ng
          if (lxlylz == nrank) then
            wbuf_rxryrz_r_gpu = wbuf_lxlylz_s_gpu
            wbuf_rxrylz_r_gpu = wbuf_lxlyrz_s_gpu
            wbuf_rxlyrz_r_gpu = wbuf_lxrylz_s_gpu
            wbuf_rxlylz_r_gpu = wbuf_lxryrz_s_gpu
            wbuf_lxryrz_r_gpu = wbuf_rxlylz_s_gpu
            wbuf_lxrylz_r_gpu = wbuf_rxlyrz_s_gpu
            wbuf_lxlyrz_r_gpu = wbuf_rxrylz_s_gpu
            wbuf_lxlylz_r_gpu = wbuf_rxryrz_s_gpu
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
        call bcswap_corner_step_3_cuf(nx, ny, nz, nv, ng, w_gpu, wbuf1r_c_gpu, wbuf2r_c_gpu,&
        & wbuf3r_c_gpu, wbuf4r_c_gpu, ileftbottom, ilefttop, irightbottom, irighttop)
        if(pbc(2)) then
          call bcswap_corner_step_3b_cuf(nx, ny, nz, nv, ng, w_gpu, wbuf_lxly_r_gpu,&
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
    class(base_gpu_object), intent(inout) :: self
    real(rkind), dimension(1-self%ng:self%nx+self%ng, 1-self%ng:self%ny+self%ng, 1-self%ng:self%nz+self%ng,1:1), device :: w_swap
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
        call bcswap_corner_step_1_cuf(nx, ny, nz, 1, ng, w_swap, wbuf1s_c_gpu, wbuf2s_c_gpu, wbuf3s_c_gpu, wbuf4s_c_gpu)
        if(pbc(2)) then
          call bcswap_corner_step_1b_cuf(nx, ny, nz, 1, ng, w_swap, wbuf_lxly_s_gpu,&
          & wbuf_lxry_s_gpu, wbuf_rxly_s_gpu, wbuf_rxry_s_gpu, wbuf_lylz_s_gpu, wbuf_lyrz_s_gpu,&
          & wbuf_rylz_s_gpu, wbuf_ryrz_s_gpu, wbuf_lxlylz_s_gpu, wbuf_lxlyrz_s_gpu, wbuf_lxrylz_s_gpu,&
          & wbuf_lxryrz_s_gpu, wbuf_rxlylz_s_gpu, wbuf_rxlyrz_s_gpu, wbuf_rxrylz_s_gpu, wbuf_rxryrz_s_gpu)
        endif
      endif
      if(steps_(2)) then
        indc = 1*ny*ng*ng
        if (ileftbottom == nrank) then
          wbuf2r_c_gpu = wbuf1s_c_gpu
          wbuf1r_c_gpu = wbuf2s_c_gpu
        else
          call mpi_sendrecv(wbuf1s_c_gpu,indc,mpi_prec,ileftbottom ,1, wbuf2r_c_gpu,indc,mpi_prec,&
          &irighttop ,1,mp_cart,istatus,iermpi)
          call mpi_sendrecv(wbuf2s_c_gpu,indc,mpi_prec,irighttop ,2, wbuf1r_c_gpu,indc,mpi_prec,&
          &ileftbottom ,2,mp_cart,istatus,iermpi)
        endif
        if (ilefttop == nrank) then
          wbuf4r_c_gpu = wbuf3s_c_gpu
          wbuf3r_c_gpu = wbuf4s_c_gpu
        else
          call mpi_sendrecv(wbuf3s_c_gpu,indc,mpi_prec,ilefttop ,3, wbuf4r_c_gpu,indc,mpi_prec,&
          &irightbottom,3,mp_cart,istatus,iermpi)
          call mpi_sendrecv(wbuf4s_c_gpu,indc,mpi_prec,irightbottom,4, wbuf3r_c_gpu,indc,mpi_prec,&
          &ilefttop ,4,mp_cart,istatus,iermpi)
        endif
        if(pbc(2)) then
          indc = 1*ng*ng*nz
          if (lxly == nrank) then
            wbuf_rxry_r_gpu = wbuf_lxly_s_gpu
            wbuf_lxly_r_gpu = wbuf_rxry_s_gpu
          else
            call mpi_sendrecv(wbuf_lxly_s_gpu,indc,mpi_prec,lxly ,5, wbuf_rxry_r_gpu,indc,mpi_prec,rxry ,5,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rxry_s_gpu,indc,mpi_prec,rxry ,6, wbuf_lxly_r_gpu,indc,mpi_prec,lxly ,6,mp_cart,istatus,iermpi)
          endif
          if (lxry == nrank) then
            wbuf_rxly_r_gpu = wbuf_lxry_s_gpu
            wbuf_lxry_r_gpu = wbuf_rxly_s_gpu
          else
            call mpi_sendrecv(wbuf_lxry_s_gpu,indc,mpi_prec,lxry ,7, wbuf_rxly_r_gpu,indc,mpi_prec,rxly ,7,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rxly_s_gpu,indc,mpi_prec,rxly ,8, wbuf_lxry_r_gpu,indc,mpi_prec,lxry ,8,mp_cart,istatus,iermpi)
          endif
          indc = 1*nx*ng*ng
          if (lylz == nrank) then
            wbuf_ryrz_r_gpu = wbuf_lylz_s_gpu
            wbuf_lylz_r_gpu = wbuf_ryrz_s_gpu
          else
            call mpi_sendrecv(wbuf_lylz_s_gpu,indc,mpi_prec,lylz ,9, wbuf_ryrz_r_gpu,indc,mpi_prec,ryrz ,9,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_ryrz_s_gpu,indc,mpi_prec,ryrz ,10, wbuf_lylz_r_gpu,indc,mpi_prec,&
            &lylz ,10,mp_cart,istatus,iermpi)
          endif
          if (lyrz == nrank) then
            wbuf_rylz_r_gpu = wbuf_lyrz_s_gpu
            wbuf_lyrz_r_gpu = wbuf_rylz_s_gpu
          else
            call mpi_sendrecv(wbuf_lyrz_s_gpu,indc,mpi_prec,lyrz ,11, wbuf_rylz_r_gpu,indc,mpi_prec,&
            &rylz ,11,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rylz_s_gpu,indc,mpi_prec,rylz ,12, wbuf_lyrz_r_gpu,indc,mpi_prec,&
            &lyrz ,12,mp_cart,istatus,iermpi)
          endif
          indc = 1*ng*ng*ng
          if (lxlylz == nrank) then
            wbuf_rxryrz_r_gpu = wbuf_lxlylz_s_gpu
            wbuf_rxrylz_r_gpu = wbuf_lxlyrz_s_gpu
            wbuf_rxlyrz_r_gpu = wbuf_lxrylz_s_gpu
            wbuf_rxlylz_r_gpu = wbuf_lxryrz_s_gpu
            wbuf_lxryrz_r_gpu = wbuf_rxlylz_s_gpu
            wbuf_lxrylz_r_gpu = wbuf_rxlyrz_s_gpu
            wbuf_lxlyrz_r_gpu = wbuf_rxrylz_s_gpu
            wbuf_lxlylz_r_gpu = wbuf_rxryrz_s_gpu
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
        call bcswap_corner_step_3_cuf(nx, ny, nz, 1, ng, w_swap, wbuf1r_c_gpu, wbuf2r_c_gpu,&
        & wbuf3r_c_gpu, wbuf4r_c_gpu, ileftbottom, ilefttop, irightbottom, irighttop)
        if(pbc(2)) then
          call bcswap_corner_step_3b_cuf(nx, ny, nz, 1, ng, w_swap, wbuf_lxly_r_gpu,&
          & wbuf_lxry_r_gpu, wbuf_rxly_r_gpu, wbuf_rxry_r_gpu, wbuf_lylz_r_gpu, wbuf_lyrz_r_gpu,&
          & wbuf_rylz_r_gpu, wbuf_ryrz_r_gpu, wbuf_lxlylz_r_gpu, wbuf_lxlyrz_r_gpu, wbuf_lxrylz_r_gpu,&
          & wbuf_lxryrz_r_gpu, wbuf_rxlylz_r_gpu, wbuf_rxlyrz_r_gpu, wbuf_rxrylz_r_gpu, wbuf_rxryrz_r_gpu,&
          & lxly, lxry, rxly, rxry, lylz, lyrz, rylz, ryrz, lxlylz, lxlyrz, lxrylz, lxryrz, rxlylz, rxlyrz,&
          & rxrylz, rxryrz)
        endif
      endif
    endassociate
  endsubroutine bcswap_edges_corners_var

  subroutine bcswap(self, steps)
    class(base_gpu_object), intent(inout) :: self
    logical, dimension(3), optional :: steps
    logical, dimension(3) :: steps_
    integer :: iercuda, indx, indy, indz
    integer, dimension(mpi_status_size) :: istatus

    steps_ = .true. ; if(present(steps)) steps_ = steps

    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv,&
    & w_gpu => self%w_gpu, wbuf1s_gpu => self%wbuf1s_gpu, wbuf2s_gpu => self%wbuf2s_gpu,&
    & wbuf3s_gpu => self%wbuf3s_gpu, wbuf4s_gpu => self%wbuf4s_gpu, wbuf5s_gpu => self%wbuf5s_gpu,&
    & wbuf6s_gpu => self%wbuf6s_gpu, wbuf1r_gpu => self%wbuf1r_gpu, wbuf2r_gpu => self%wbuf2r_gpu,&
    & wbuf3r_gpu => self%wbuf3r_gpu, wbuf4r_gpu => self%wbuf4r_gpu, wbuf5r_gpu => self%wbuf5r_gpu,&
    & wbuf6r_gpu => self%wbuf6r_gpu, ileftx => self%field%ileftx,irightx => self%field%irightx,&
    & nrank_x => self%field%nrank_x, ilefty => self%field%ilefty,irighty => self%field%irighty,&
    & nrank_y => self%field%nrank_y, ileftz => self%field%ileftz,irightz => self%field%irightz,&
    & nrank_z => self%field%nrank_z, mp_cartx => self%field%mp_cartx, mp_carty => self%field%mp_carty,&
    & mp_cartz => self%field%mp_cartz)

      if(steps_(1)) then
        call bcswap_step_1_cuf(nx, ny, nz, nv, ng, w_gpu, wbuf1s_gpu, wbuf2s_gpu, wbuf3s_gpu, wbuf4s_gpu, wbuf5s_gpu, wbuf6s_gpu)
      endif
      if(steps_(2)) then
        indx = nv*ng*ny*nz
        indy = nv*nx*ng*nz
        indz = nv*nx*ny*ng
        if(ileftx == nrank_x) then
          wbuf2r_gpu = wbuf1s_gpu
          wbuf1r_gpu = wbuf2s_gpu
        else
          call mpi_sendrecv(wbuf1s_gpu,indx,mpi_prec,ileftx ,1,wbuf2r_gpu,indx,mpi_prec,irightx,1,mp_cartx,istatus,self%mpi_err)
          call mpi_sendrecv(wbuf2s_gpu,indx,mpi_prec,irightx,2,wbuf1r_gpu,indx,mpi_prec,ileftx ,2,mp_cartx,istatus,self%mpi_err)
        endif
        if(ilefty == nrank_y) then
          wbuf4r_gpu = wbuf3s_gpu
          wbuf3r_gpu = wbuf4s_gpu
        else
          call mpi_sendrecv(wbuf3s_gpu,indy,mpi_prec,ilefty ,3,wbuf4r_gpu,indy,mpi_prec,irighty,3,mp_carty,istatus,self%mpi_err)
          call mpi_sendrecv(wbuf4s_gpu,indy,mpi_prec,irighty,4,wbuf3r_gpu,indy,mpi_prec,ilefty ,4,mp_carty,istatus,self%mpi_err)
        endif
        if(ileftz == nrank_z) then
          wbuf6r_gpu = wbuf5s_gpu
          wbuf5r_gpu = wbuf6s_gpu
        else
          call mpi_sendrecv(wbuf5s_gpu,indz,mpi_prec,ileftz ,5,wbuf6r_gpu,indz,mpi_prec,irightz,5,mp_cartz,istatus,self%mpi_err)
          call mpi_sendrecv(wbuf6s_gpu,indz,mpi_prec,irightz,6,wbuf5r_gpu,indz,mpi_prec,ileftz ,6,mp_cartz,istatus,self%mpi_err)
        endif
      endif
      if(steps_(3)) then
        call bcswap_step_3_cuf(nx, ny, nz, nv, ng, w_gpu, wbuf1r_gpu, wbuf2r_gpu, wbuf3r_gpu,&
        & wbuf4r_gpu, wbuf5r_gpu, wbuf6r_gpu, ileftx, ilefty, ileftz, irightx, irighty, irightz)
      endif
    endassociate
  endsubroutine bcswap

  subroutine bcswap_var(self, w_swap, steps)
    class(base_gpu_object), intent(inout) :: self
    real(rkind), dimension(1-self%ng:self%nx+self%ng, 1-self%ng:self%ny+self%ng, 1-self%ng:self%nz+self%ng,1:1), device :: w_swap
    logical, dimension(3), optional :: steps
    logical, dimension(3) :: steps_
    integer :: iercuda, indx, indy, indz
    integer, dimension(mpi_status_size) :: istatus

    steps_ = .true. ; if(present(steps)) steps_ = steps

    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, wbuf1s_gpu => self%wbuf1s_&
    &gpu, wbuf2s_gpu => self%wbuf2s_gpu, wbuf3s_gpu => self%wbuf3s_gpu, wbuf4s_gpu => self%wbuf4s_gpu,&
    & wbuf5s_gpu => self%wbuf5s_gpu, wbuf6s_gpu => self%wbuf6s_gpu, wbuf1r_gpu => self%wbuf1r_gpu,&
    & wbuf2r_gpu => self%wbuf2r_gpu, wbuf3r_gpu => self%wbuf3r_gpu, wbuf4r_gpu => self%wbuf4r_gpu,&
    & wbuf5r_gpu => self%wbuf5r_gpu, wbuf6r_gpu => self%wbuf6r_gpu, ileftx => self%field%ileftx,&
    &irightx => self%field%irightx, nrank_x => self%field%nrank_x, ilefty => self%field%ilefty,&
    &irighty => self%field%irighty, nrank_y => self%field%nrank_y, ileftz => self%field%ileftz,&
    &irightz => self%field%irightz, nrank_z => self%field%nrank_z, mp_cartx => self%field%mp_cartx,&
    & mp_carty => self%field%mp_carty, mp_cartz => self%field%mp_cartz)

      if(steps_(1)) then
        call bcswap_step_1_cuf(nx, ny, nz, 1, ng, w_swap, wbuf1s_gpu, wbuf2s_gpu, wbuf3s_gpu, wbuf4s_gpu, wbuf5s_gpu, wbuf6s_gpu)
      endif
      if(steps_(2)) then
        indx = ng*ny*nz
        indy = nx*ng*nz
        indz = nx*ny*ng
        if(ileftx == nrank_x) then
          wbuf2r_gpu = wbuf1s_gpu
          wbuf1r_gpu = wbuf2s_gpu
        else
          call mpi_sendrecv(wbuf1s_gpu,indx,mpi_prec,ileftx ,1,wbuf2r_gpu,indx,mpi_prec,irightx,1,mp_cartx,istatus,self%mpi_err)
          call mpi_sendrecv(wbuf2s_gpu,indx,mpi_prec,irightx,2,wbuf1r_gpu,indx,mpi_prec,ileftx ,2,mp_cartx,istatus,self%mpi_err)
        endif
        if(ilefty == nrank_y) then
          wbuf4r_gpu = wbuf3s_gpu
          wbuf3r_gpu = wbuf4s_gpu
        else
          call mpi_sendrecv(wbuf3s_gpu,indy,mpi_prec,ilefty ,3,wbuf4r_gpu,indy,mpi_prec,irighty,3,mp_carty,istatus,self%mpi_err)
          call mpi_sendrecv(wbuf4s_gpu,indy,mpi_prec,irighty,4,wbuf3r_gpu,indy,mpi_prec,ilefty ,4,mp_carty,istatus,self%mpi_err)
        endif
        if(ileftz == nrank_z) then
          wbuf6r_gpu = wbuf5s_gpu
          wbuf5r_gpu = wbuf6s_gpu
        else
          call mpi_sendrecv(wbuf5s_gpu,indz,mpi_prec,ileftz ,5,wbuf6r_gpu,indz,mpi_prec,irightz,5,mp_cartz,istatus,self%mpi_err)
          call mpi_sendrecv(wbuf6s_gpu,indz,mpi_prec,irightz,6,wbuf5r_gpu,indz,mpi_prec,ileftz ,6,mp_cartz,istatus,self%mpi_err)
        endif
      endif
      if(steps_(3)) then
        call bcswap_step_3_cuf(nx, ny, nz, 1, ng, w_swap, wbuf1r_gpu, wbuf2r_gpu, wbuf3r_gpu,&
        & wbuf4r_gpu, wbuf5r_gpu, wbuf6r_gpu, ileftx, ilefty, ileftz, irightx, irighty, irightz)
      endif
    endassociate
  endsubroutine bcswap_var

  subroutine copy_cpu_gpu(self)
    class(base_gpu_object), intent(inout) :: self
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

    self%w_gpu = self%w_t

  endsubroutine copy_cpu_gpu


  subroutine copy_gpu_cpu(self)
    class(base_gpu_object), intent(inout) :: self
    integer :: i,j,k,iv

    self%w_t = self%w_gpu

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

  endsubroutine copy_gpu_cpu

  subroutine initialize(self, field)
    class(base_gpu_object), intent(inout) :: self
    class(field_object), target :: field
    type(cudadeviceprop) :: device_properties
    real(rkind) :: device_mem_avail

    self%field => field
    self%nx = self%field%nx
    self%ny = self%field%ny
    self%nz = self%field%nz
    self%ng = self%field%grid%ng
    self%nv = self%field%nv

    call get_mpi_basic_info(self%nprocs, self%myrank, self%masterproc, self%mpi_err)

    call self%field%check_cpu_mem(description="--base-initialization--")

    call mpi_comm_split_type(mpi_comm_world, mpi_comm_type_shared, 0, mpi_info_null, self%local_comm, self%mpi_err)
    call mpi_comm_rank(self%local_comm, self%mydev, self%mpi_err)
    self%mpi_err = cudasetdevice(self%mydev)
    self%mpi_err = cudagetdeviceproperties(device_properties, self%mydev)

    device_mem_avail = real(device_properties%totalglobalmem, rkind)/(1024_rkind**3)

    call self%alloc()

  endsubroutine initialize

  subroutine check_gpu_mem(self, description)
    class(base_gpu_object), intent(inout) :: self
    character(*) :: description
    integer :: ierr
    integer(cuda_count_kind) :: mem_free, mem_total
    character(128) :: proc_name
    integer :: resultlen
    call mpi_get_processor_name(proc_name, resultlen, ierr)
    ierr = cudamemgetinfo(mem_free, mem_total)
    write(error_unit, "(a,2x,a,2x,a,2x,i0,2x,i0,2x,i0)") 'gpu rank,mems: ', description,&
    &proc_name(1:resultlen),self%myrank, mem_free, mem_total
  endsubroutine check_gpu_mem

endmodule streams_base_gpu_object

