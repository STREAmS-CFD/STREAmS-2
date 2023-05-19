module streams_base_gpu_object
! < STREAmS, base GPU class definition
!
  use streams_field_object, only : field_object
  use streams_parameters
  use MPI
  use CUDAFOR
!
  implicit none
  private
  public :: base_gpu_object
!
  type :: base_gpu_object
    type(field_object), pointer :: field=>null()
!   Replica of field and grid sizes
    integer :: nx, ny, nz, ng, nv
!   MPI data
    integer(ikind)              :: myrank=0_ikind       !< MPI rank process.
    integer(ikind)              :: nprocs=1_ikind       !< Number of MPI processes.
    logical                     :: masterproc
    integer(ikind)              :: mpi_err=0_ikind      !< Error traping flag.
    integer(ikind)              :: mydev=0_ikind        !< My GPU rank.
    integer(ikind)              :: local_comm=0_ikind   !< Local communicator.
!   GPU data
    real(rkind), allocatable, dimension(:,:,:,:), device :: w_gpu
    real(rkind), allocatable, dimension(:,:,:,:)         :: w_t
!
    real(rkind), dimension(:,:,:,:), device, allocatable :: wbuf1s_gpu, wbuf2s_gpu, wbuf3s_gpu, &
    wbuf4s_gpu, wbuf5s_gpu, wbuf6s_gpu
    real(rkind), dimension(:,:,:,:), device, allocatable :: wbuf1r_gpu, wbuf2r_gpu, wbuf3r_gpu, &
    wbuf4r_gpu, wbuf5r_gpu, wbuf6r_gpu
    real(rkind), dimension(:,:,:,:), allocatable, device :: wbuf1s_c_gpu, wbuf2s_c_gpu, wbuf3s_c_gpu, wbuf4s_c_gpu
    real(rkind), dimension(:,:,:,:), allocatable, device :: wbuf1r_c_gpu, wbuf2r_c_gpu, wbuf3r_c_gpu, wbuf4r_c_gpu
    real(rkind), allocatable, dimension(:), device :: dcsidx_gpu, dcsidx2_gpu, dcsidxs_gpu
    real(rkind), allocatable, dimension(:), device :: detady_gpu, detady2_gpu, detadys_gpu
    real(rkind), allocatable, dimension(:), device :: dzitdz_gpu, dzitdz2_gpu, dzitdzs_gpu
    real(rkind), allocatable, dimension(:), device :: x_gpu, y_gpu, z_gpu, yn_gpu
  contains
!   public methods
    procedure, pass(self) :: alloc
    procedure, pass(self) :: copy_cpu_gpu
    procedure, pass(self) :: copy_gpu_cpu
    procedure, pass(self) :: initialize
    procedure, pass(self) :: bcswap
    procedure, pass(self) :: bcswap_var
    procedure, pass(self) :: bcswap_corner
    procedure, pass(self) :: bcswap_corner_var
    procedure, pass(self) :: check_gpu_mem
  endtype base_gpu_object
!
! interface assign_allocatable_gpu
!   !< Safe assign allocatable arrays (GPU), generic interface.
!   module procedure assign_allocatable_INT64_1D_gpu !< Safe assign allocatable arrays, I8P 1D type.
!   module procedure assign_allocatable_INT64_2D_gpu !< Safe assign allocatable arrays, I8P 2D type.
!   module procedure assign_allocatable_INT32_1D_gpu !< Safe assign allocatable arrays, ikind 1D type.
! endinterface assign_allocatable_gpu
!
contains
! public methods
  subroutine alloc(self)
    class(base_gpu_object), intent(inout)        :: self     !< The base backend.
!
    associate(nx => self%field%nx, ny => self%field%ny, nz => self%field%nz, &
      ng => self%field%grid%ng, nv => self%field%nv)
!
      allocate(self%w_gpu(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng, nv))
      allocate(self%w_t(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng, nv))
!
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
!
      allocate(self%x_gpu(1-ng:nx+ng), self%y_gpu(1-ng:ny+ng), self%z_gpu(1-ng:nz+ng))
      allocate(self%yn_gpu(1:ny+1))
!
      allocate(self%dcsidx_gpu(nx), self%dcsidx2_gpu(nx), self%dcsidxs_gpu(nx))
      allocate(self%detady_gpu(ny), self%detady2_gpu(ny), self%detadys_gpu(ny))
      allocate(self%dzitdz_gpu(nz), self%dzitdz2_gpu(nz), self%dzitdzs_gpu(nz))
!
      allocate(self%wbuf1s_c_gpu(nv,ng,ny,ng))
      allocate(self%wbuf2s_c_gpu(nv,ng,ny,ng))
      allocate(self%wbuf3s_c_gpu(nv,ng,ny,ng))
      allocate(self%wbuf4s_c_gpu(nv,ng,ny,ng))
      allocate(self%wbuf1r_c_gpu(nv,ng,ny,ng))
      allocate(self%wbuf2r_c_gpu(nv,ng,ny,ng))
      allocate(self%wbuf3r_c_gpu(nv,ng,ny,ng))
      allocate(self%wbuf4r_c_gpu(nv,ng,ny,ng))
!
    endassociate
!
    self%x_gpu  = self%field%x
    self%y_gpu  = self%field%y
    self%z_gpu  = self%field%z
    self%yn_gpu = self%field%yn
!
    self%dcsidx_gpu  = self%field%dcsidx
    self%dcsidxs_gpu = self%field%dcsidxs
    self%dcsidx2_gpu = self%field%dcsidx2
    self%detady_gpu  = self%field%detady
    self%detadys_gpu = self%field%detadys
    self%detady2_gpu = self%field%detady2
    self%dzitdz_gpu  = self%field%dzitdz
    self%dzitdzs_gpu = self%field%dzitdzs
    self%dzitdz2_gpu = self%field%dzitdz2
!
  endsubroutine alloc
!
  subroutine bcswap_step_1_cuf(nx, ny, nz, nv, ng, w_gpu, &
    wbuf1s_gpu, wbuf2s_gpu, wbuf3s_gpu, wbuf4s_gpu, wbuf5s_gpu, wbuf6s_gpu)
    integer :: nx, ny, nz, ng, nv
    real(rkind), device, dimension(1-ng:,1-ng:,1-ng:,1:) :: w_gpu
    real(rkind), device, dimension(1:,1:,1:,1:) :: wbuf1s_gpu, wbuf2s_gpu, wbuf3s_gpu, wbuf4s_gpu, wbuf5s_gpu, wbuf6s_gpu
    integer :: i,j,k,m,iercuda
!
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
    !@cuf iercuda=cudaDeviceSynchronize()
!
!   CUDA check (debugging purposes)
!   iercuda = cudaGetLastError()
!   if (iercuda /= cudaSuccess) then
!     print*,"CUDA ERROR! ",cudaGetErrorString(iercuda)
!   endif
!
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
    !@cuf iercuda=cudaDeviceSynchronize()
!
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
    !@cuf iercuda=cudaDeviceSynchronize()
!
  endsubroutine bcswap_step_1_cuf
!
  subroutine bcswap_corner_step_1_cuf(nx, ny, nz, nv, ng, w_gpu, &
    wbuf1s_c_gpu, wbuf2s_c_gpu, wbuf3s_c_gpu, wbuf4s_c_gpu)
    integer :: nx, ny, nz, ng, nv
    real(rkind), device, dimension(1-ng:,1-ng:,1-ng:,1:) :: w_gpu
    real(rkind), device, dimension(1:,1:,1:,1:) :: wbuf1s_c_gpu, wbuf2s_c_gpu, wbuf3s_c_gpu, wbuf4s_c_gpu
    integer :: i,j,k,m,iercuda
!
    !$cuf kernel do(2) <<<*,*>>>
    do j=1,ny
      do m=1,nv
        do k=1,ng
          do i=1,ng
            wbuf1s_c_gpu(m,i,j,k) = w_gpu(i,j,k,m)              ! ileftbottom
            wbuf2s_c_gpu(m,i,j,k) = w_gpu(nx-ng+i,j,nz-ng+k,m)  ! irighttop
            wbuf3s_c_gpu(m,i,j,k) = w_gpu(i,j,nz-ng+k,m)        ! ilefttop
            wbuf4s_c_gpu(m,i,j,k) = w_gpu(nx-ng+i,j,k,m)        ! irightbottom
          enddo
        enddo
      enddo
    enddo
    !@cuf iercuda=cudaDeviceSynchronize()
!
  endsubroutine bcswap_corner_step_1_cuf
!
  subroutine bcswap_step_3_cuf(nx, ny, nz, nv, ng, w_gpu, &
    wbuf1r_gpu, wbuf2r_gpu, wbuf3r_gpu, wbuf4r_gpu, wbuf5r_gpu, wbuf6r_gpu, &
    ileftx, ilefty, ileftz, irightx, irighty, irightz)
    integer :: nx, ny, nz, ng, nv
    real(rkind), device, dimension(1-ng:,1-ng:,1-ng:,1:) :: w_gpu
    real(rkind), device, dimension(1:,1:,1:,1:) :: wbuf1r_gpu, wbuf2r_gpu, wbuf3r_gpu, wbuf4r_gpu, wbuf5r_gpu, wbuf6r_gpu
    integer :: i,j,k,m, iercuda
    integer ::  ileftx, ilefty, ileftz, irightx, irighty, irightz
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
      !@cuf iercuda=cudaDeviceSynchronize()
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
      !@cuf iercuda=cudaDeviceSynchronize()
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
      !@cuf iercuda=cudaDeviceSynchronize()
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
      !@cuf iercuda=cudaDeviceSynchronize()
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
      !@cuf iercuda=cudaDeviceSynchronize()
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
      !@cuf iercuda=cudaDeviceSynchronize()
    endif
  endsubroutine bcswap_step_3_cuf
!
  subroutine bcswap_corner_step_3_cuf(nx, ny, nz, nv, ng, w_gpu, &
    wbuf1r_c_gpu, wbuf2r_c_gpu, wbuf3r_c_gpu, wbuf4r_c_gpu, &
    ileftbottom, ilefttop, irightbottom, irighttop)
    integer :: nx, ny, nz, ng, nv
    real(rkind), device, dimension(1-ng:,1-ng:,1-ng:,1:) :: w_gpu
    real(rkind), device, dimension(1:,1:,1:,1:) :: wbuf1r_c_gpu, wbuf2r_c_gpu, wbuf3r_c_gpu, wbuf4r_c_gpu
    integer :: i,j,k,m, iercuda
    integer :: ileftbottom, ilefttop, irightbottom, irighttop
!
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
      !@cuf iercuda=cudaDeviceSynchronize()
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
      !@cuf iercuda=cudaDeviceSynchronize()
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
      !@cuf iercuda=cudaDeviceSynchronize()
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
      !@cuf iercuda=cudaDeviceSynchronize()
    endif
  endsubroutine bcswap_corner_step_3_cuf
!
  subroutine bcswap_corner(self, steps)
    class(base_gpu_object), intent(inout) :: self !< The base backend.
    logical, dimension(3), optional :: steps
    logical, dimension(3) :: steps_
    integer :: iercuda, indc
    integer, dimension(mpi_status_size) :: istatus
!
    steps_ = .true. ; if(present(steps)) steps_ = steps
!
    associate(nx => self%nx, ny => self%ny, nz => self%nz, &
      ng => self%ng, nv => self%nv, &
      w_gpu => self%w_gpu, &
      wbuf1s_c_gpu => self%wbuf1s_c_gpu, wbuf2s_c_gpu => self%wbuf2s_c_gpu,            &
      wbuf3s_c_gpu => self%wbuf3s_c_gpu, wbuf4s_c_gpu => self%wbuf4s_c_gpu,            &
      wbuf1r_c_gpu => self%wbuf1r_c_gpu, wbuf2r_c_gpu => self%wbuf2r_c_gpu,            &
      wbuf3r_c_gpu => self%wbuf3r_c_gpu, wbuf4r_c_gpu => self%wbuf4r_c_gpu,            &
      ileftbottom  => self%field%ileftbottom, irightbottom => self%field%irightbottom, &
      ilefttop     => self%field%ilefttop,    irighttop => self%field%irighttop,       &
      mp_cart      => self%field%mp_cart, iermpi => self%mpi_err, nrank => self%field%myrank)
!
      if(steps_(1)) then
        call bcswap_corner_step_1_cuf(nx, ny, nz, nv, ng, w_gpu, &
          wbuf1s_c_gpu, wbuf2s_c_gpu, wbuf3s_c_gpu, wbuf4s_c_gpu)
      endif
      if(steps_(2)) then
        indc = nv*ny*ng*ng
        if (ileftbottom == nrank) then
          wbuf2r_c_gpu = wbuf1s_c_gpu
          wbuf1r_c_gpu = wbuf2s_c_gpu
        else
          call mpi_sendrecv(wbuf1s_c_gpu,indc,mpi_prec,ileftbottom ,1, &
            wbuf2r_c_gpu,indc,mpi_prec,irighttop   ,1,mp_cart,istatus,iermpi)
          call mpi_sendrecv(wbuf2s_c_gpu,indc,mpi_prec,irighttop   ,2, &
            wbuf1r_c_gpu,indc,mpi_prec,ileftbottom ,2,mp_cart,istatus,iermpi)
        endif
        if (ilefttop == nrank) then
          wbuf4r_c_gpu = wbuf3s_c_gpu
          wbuf3r_c_gpu = wbuf4s_c_gpu
        else
          call mpi_sendrecv(wbuf3s_c_gpu,indc,mpi_prec,ilefttop    ,3, &
            wbuf4r_c_gpu,indc,mpi_prec,irightbottom,3,mp_cart,istatus,iermpi)
          call mpi_sendrecv(wbuf4s_c_gpu,indc,mpi_prec,irightbottom,4, &
            wbuf3r_c_gpu,indc,mpi_prec,ilefttop    ,4,mp_cart,istatus,iermpi)
        endif
      endif
      if(steps_(3)) then
        call bcswap_corner_step_3_cuf(nx, ny, nz, nv, ng, w_gpu, &
          wbuf1r_c_gpu, wbuf2r_c_gpu, wbuf3r_c_gpu, wbuf4r_c_gpu, &
          ileftbottom, ilefttop, irightbottom, irighttop)
      endif
    endassociate
  endsubroutine bcswap_corner
!
  subroutine bcswap_corner_var(self, w_swap, steps)
    class(base_gpu_object), intent(inout) :: self !< The base backend.
    real(rkind), dimension(1-self%ng:self%nx+self%ng, &
    1-self%ng:self%ny+self%ng, &
    1-self%ng:self%nz+self%ng,1:1), device :: w_swap
    logical, dimension(3), optional :: steps
    logical, dimension(3) :: steps_
    integer :: iercuda, indc
    integer, dimension(mpi_status_size) :: istatus
!
    steps_ = .true. ; if(present(steps)) steps_ = steps
!
    associate(nx => self%nx, ny => self%ny, nz => self%nz, &
      ng => self%ng, nv => self%nv, &
      wbuf1s_c_gpu => self%wbuf1s_c_gpu, wbuf2s_c_gpu => self%wbuf2s_c_gpu,            &
      wbuf3s_c_gpu => self%wbuf3s_c_gpu, wbuf4s_c_gpu => self%wbuf4s_c_gpu,            &
      wbuf1r_c_gpu => self%wbuf1r_c_gpu, wbuf2r_c_gpu => self%wbuf2r_c_gpu,            &
      wbuf3r_c_gpu => self%wbuf3r_c_gpu, wbuf4r_c_gpu => self%wbuf4r_c_gpu,            &
      ileftbottom  => self%field%ileftbottom, irightbottom => self%field%irightbottom, &
      ilefttop     => self%field%ilefttop,    irighttop => self%field%irighttop,       &
      mp_cart      => self%field%mp_cart, iermpi => self%mpi_err, nrank => self%field%myrank)
!
      if(steps_(1)) then
        call bcswap_corner_step_1_cuf(nx, ny, nz, 1, ng, w_swap, &
          wbuf1s_c_gpu, wbuf2s_c_gpu, wbuf3s_c_gpu, wbuf4s_c_gpu)
      endif
      if(steps_(2)) then
        indc = ny*ng*ng
        if (ileftbottom == nrank) then
          wbuf2r_c_gpu = wbuf1s_c_gpu
          wbuf1r_c_gpu = wbuf2s_c_gpu
        else
          call mpi_sendrecv(wbuf1s_c_gpu,indc,mpi_prec,ileftbottom ,1, &
            wbuf2r_c_gpu,indc,mpi_prec,irighttop   ,1,mp_cart,istatus,iermpi)
          call mpi_sendrecv(wbuf2s_c_gpu,indc,mpi_prec,irighttop   ,2, &
            wbuf1r_c_gpu,indc,mpi_prec,ileftbottom ,2,mp_cart,istatus,iermpi)
        endif
        if (ilefttop == nrank) then
          wbuf4r_c_gpu = wbuf3s_c_gpu
          wbuf3r_c_gpu = wbuf4s_c_gpu
        else
          call mpi_sendrecv(wbuf3s_c_gpu,indc,mpi_prec,ilefttop    ,3, &
            wbuf4r_c_gpu,indc,mpi_prec,irightbottom,3,mp_cart,istatus,iermpi)
          call mpi_sendrecv(wbuf4s_c_gpu,indc,mpi_prec,irightbottom,4, &
            wbuf3r_c_gpu,indc,mpi_prec,ilefttop    ,4,mp_cart,istatus,iermpi)
        endif
      endif
      if(steps_(3)) then
        call bcswap_corner_step_3_cuf(nx, ny, nz, 1, ng, w_swap, &
          wbuf1r_c_gpu, wbuf2r_c_gpu, wbuf3r_c_gpu, wbuf4r_c_gpu, &
          ileftbottom, ilefttop, irightbottom, irighttop)
      endif
    endassociate
  endsubroutine bcswap_corner_var
!
  subroutine bcswap(self, steps)
    class(base_gpu_object), intent(inout) :: self !< The base backend.
    logical, dimension(3), optional :: steps
    logical, dimension(3) :: steps_
    integer :: iercuda, indx, indy, indz
    integer, dimension(mpi_status_size) :: istatus
!
    steps_ = .true. ; if(present(steps)) steps_ = steps
!
    associate(nx => self%nx, ny => self%ny, nz => self%nz, &
      ng => self%ng, nv => self%nv, &
      w_gpu => self%w_gpu, &
      wbuf1s_gpu => self%wbuf1s_gpu, wbuf2s_gpu => self%wbuf2s_gpu, &
      wbuf3s_gpu => self%wbuf3s_gpu, wbuf4s_gpu => self%wbuf4s_gpu, &
      wbuf5s_gpu => self%wbuf5s_gpu, wbuf6s_gpu => self%wbuf6s_gpu, &
      wbuf1r_gpu => self%wbuf1r_gpu, wbuf2r_gpu => self%wbuf2r_gpu, &
      wbuf3r_gpu => self%wbuf3r_gpu, wbuf4r_gpu => self%wbuf4r_gpu, &
      wbuf5r_gpu => self%wbuf5r_gpu, wbuf6r_gpu => self%wbuf6r_gpu, &
      ileftx => self%field%ileftx,irightx => self%field%irightx, nrank_x => self%field%nrank_x, &
      ilefty => self%field%ilefty,irighty => self%field%irighty, nrank_y => self%field%nrank_y, &
      ileftz => self%field%ileftz,irightz => self%field%irightz, nrank_z => self%field%nrank_z, &
      mp_cartx => self%field%mp_cartx, mp_carty => self%field%mp_carty, mp_cartz => self%field%mp_cartz)
!
      if(steps_(1)) then
        call bcswap_step_1_cuf(nx, ny, nz, nv, ng, w_gpu, &
          wbuf1s_gpu, wbuf2s_gpu, wbuf3s_gpu, wbuf4s_gpu, wbuf5s_gpu, wbuf6s_gpu)
      endif
      if(steps_(2)) then
        indx = nv*ng*ny*nz
        indy = nv*nx*ng*nz
        indz = nv*nx*ny*ng
!       http://developer.download.nvidia.com/compute/cuda/2_3/toolkit/...
!       docs/online/group__CUDART__MEMORY_ge4366f68c6fa8c85141448f187d2aa13.html
!       IMPORTANT NOTE: Copies with kind == cudaMemcpyDeviceToDevice are asynchronous
!       with respect to the host, but never overlap with kernel execution
!       For this reason here simple copies are used and not cudaMemcpyAsync D2D
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
        call bcswap_step_3_cuf(nx, ny, nz, nv, ng, w_gpu, &
          wbuf1r_gpu, wbuf2r_gpu, wbuf3r_gpu, wbuf4r_gpu, wbuf5r_gpu, wbuf6r_gpu, &
          ileftx, ilefty, ileftz, irightx, irighty, irightz)
      endif
    endassociate
  endsubroutine bcswap
!
  subroutine bcswap_var(self, w_swap, steps)
    class(base_gpu_object), intent(inout) :: self !< The base backend.
    real(rkind), dimension(1-self%ng:self%nx+self%ng, &
    1-self%ng:self%ny+self%ng, &
    1-self%ng:self%nz+self%ng,1:1), device :: w_swap
    logical, dimension(3), optional :: steps
    logical, dimension(3) :: steps_
    integer :: iercuda, indx, indy, indz
    integer, dimension(mpi_status_size) :: istatus
!
    steps_ = .true. ; if(present(steps)) steps_ = steps
!
    associate(nx => self%nx, ny => self%ny, nz => self%nz, &
      ng => self%ng, &
      wbuf1s_gpu => self%wbuf1s_gpu, wbuf2s_gpu => self%wbuf2s_gpu, &
      wbuf3s_gpu => self%wbuf3s_gpu, wbuf4s_gpu => self%wbuf4s_gpu, &
      wbuf5s_gpu => self%wbuf5s_gpu, wbuf6s_gpu => self%wbuf6s_gpu, &
      wbuf1r_gpu => self%wbuf1r_gpu, wbuf2r_gpu => self%wbuf2r_gpu, &
      wbuf3r_gpu => self%wbuf3r_gpu, wbuf4r_gpu => self%wbuf4r_gpu, &
      wbuf5r_gpu => self%wbuf5r_gpu, wbuf6r_gpu => self%wbuf6r_gpu, &
      ileftx => self%field%ileftx,irightx => self%field%irightx, nrank_x => self%field%nrank_x, &
      ilefty => self%field%ilefty,irighty => self%field%irighty, nrank_y => self%field%nrank_y, &
      ileftz => self%field%ileftz,irightz => self%field%irightz, nrank_z => self%field%nrank_z, &
      mp_cartx => self%field%mp_cartx, mp_carty => self%field%mp_carty, mp_cartz => self%field%mp_cartz)
!
      if(steps_(1)) then
        call bcswap_step_1_cuf(nx, ny, nz, 1, ng, w_swap, &
          wbuf1s_gpu, wbuf2s_gpu, wbuf3s_gpu, wbuf4s_gpu, wbuf5s_gpu, wbuf6s_gpu)
      endif
      if(steps_(2)) then
        indx = ng*ny*nz
        indy = nx*ng*nz
        indz = nx*ny*ng
!       http://developer.download.nvidia.com/compute/cuda/2_3/toolkit...
!       /docs/online/group__CUDART__MEMORY_ge4366f68c6fa8c85141448f187d2aa13.html
!       IMPORTANT NOTE: Copies with kind == cudaMemcpyDeviceToDevice are asynchronous
!       with respect to the host, but never overlap with kernel execution
!       For this reason here simple copies are used and not cudaMemcpyAsync D2D
        if(ileftx == nrank_x) then
          wbuf2r_gpu = wbuf1s_gpu ! inefficient, only 1 variable here to be copied
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
        call bcswap_step_3_cuf(nx, ny, nz, 1, ng, w_swap, &
          wbuf1r_gpu, wbuf2r_gpu, wbuf3r_gpu, wbuf4r_gpu, wbuf5r_gpu, wbuf6r_gpu, &
          ileftx, ilefty, ileftz, irightx, irighty, irightz)
      endif
    endassociate
  endsubroutine bcswap_var
!
  subroutine copy_cpu_gpu(self)
!   < Copy data from CPU to GPU.
    class(base_gpu_object), intent(inout) :: self !< The base backend.
    integer :: i,j,k,iv
!
    associate(nx => self%field%nx, ny => self%field%ny, nz => self%field%nz, &
      ng => self%field%grid%ng, nv => self%field%nv)
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
!
!   copy cpu to gpu
    self%w_gpu = self%w_t
!
  endsubroutine copy_cpu_gpu
!
! subroutine copy_gpu_cpu_var(self, w_var_gpu, w_var)
!   class(base_gpu_object), intent(inout) :: self
!   real(rkind), dimension(1-self%ng:, 1-self%ng, 1-self%ng, :), device :: w_io_gpu
!   ! copy gpu to cpu
!   self%w_var_t = self%w_var_gpu
!   associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng)
!     do k=1-ng,nz+ng
!       do j=1-ng,ny+ng
!         do i=1-ng,nx+ng
!           w_var(i,j,k) = self%w_var_t(i,j,k,iv)
!         enddo
!       enddo
!     enddo
!   endassociate
! endsubroutine copy_gpu_cpu_var
!
  subroutine copy_gpu_cpu(self)
!   < Copy data from CPU to GPU.
    class(base_gpu_object), intent(inout) :: self !< The base backend.
    integer :: i,j,k,iv
!
!   copy gpu to cpu
    self%w_t = self%w_gpu
!
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
!
  endsubroutine copy_gpu_cpu
!
  subroutine initialize(self, field)
!   < Initialize base backend.
    class(base_gpu_object), intent(inout)        :: self               !< The base backend.
    class(field_object), target :: field
    type(cudadeviceprop)                         :: device_properties  !< Device properties.
    real(rkind)                     :: device_mem_avail   !< Device memory available (Gb).
!
    self%field => field
    self%nx    = self%field%nx
    self%ny    = self%field%ny
    self%nz    = self%field%nz
    self%ng    = self%field%grid%ng
    self%nv    = self%field%nv
!
    call get_mpi_basic_info(self%nprocs, self%myrank, self%masterproc, self%mpi_err)
!
    call self%field%check_cpu_mem(description="--Base-initialization--")
!
    call MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, self%local_comm, self%mpi_err)
    call MPI_COMM_RANK(self%local_comm, self%mydev, self%mpi_err)
    self%mpi_err = CudaSetDevice(self%mydev)
    self%mpi_err = cudaGetDeviceProperties(device_properties, self%mydev)
!
    device_mem_avail = real(device_properties%totalGlobalMem, rkind)/(1024_rkind**3)
!   write(*,'(A,F5.2,A)') ' available device memory ', device_mem_avail, ' Gb'
!
    call self%alloc()
!
  endsubroutine initialize
!
  subroutine check_gpu_mem(self, description)
    class(base_gpu_object), intent(inout) :: self !< The base backend.
    character(*) :: description
    integer                       :: ierr
    integer(cuda_count_kind)      :: mem_free, mem_total
    character(128) :: proc_name
    integer :: resultlen
!   call mpi_barrier(mpi_comm_world, ierr)
    call mpi_get_processor_name(proc_name, resultlen, ierr)
    ierr = cudaMemGetInfo(mem_free, mem_total)
!   call mpi_barrier(mpi_comm_world, ierr)
    write(error_unit, "(A,2x,A,2x,A,2x,I0,2x,I0,2x,I0)") 'GPU rank,mems: ', &
      description,proc_name(1:resultlen),self%myrank, mem_free, mem_total
  endsubroutine check_gpu_mem
!
endmodule streams_base_gpu_object
