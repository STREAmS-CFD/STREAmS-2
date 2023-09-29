module streams_base_gpu_object
!
  use streams_field_object, only : field_object
  use streams_parameters
  use mpi
  use cudafor
!
  implicit none
  private
  public :: base_gpu_object
!
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
!
    real(rkind), dimension(:,:,:,:), device, allocatable :: wbuf1s_gpu, wbuf2s_gpu, wbuf3s_gpu, wbuf4s_gpu, wbuf5s_gpu, wbuf6s_gpu
    real(rkind), dimension(:,:,:,:), device, allocatable :: wbuf1r_gpu, wbuf2r_gpu, wbuf3r_gpu, wbuf4r_gpu, wbuf5r_gpu, wbuf6r_gpu
    real(rkind), dimension(:,:,:,:), allocatable, device :: wbuf1s_c_gpu, wbuf2s_c_gpu, wbuf3s_c_gpu, wbuf4s_c_gpu
    real(rkind), dimension(:,:,:,:), allocatable, device :: wbuf1r_c_gpu, wbuf2r_c_gpu, wbuf3r_c_gpu, wbuf4r_c_gpu
    real(rkind), allocatable, dimension(:), device :: dcsidx_gpu, dcsidx2_gpu, dcsidxs_gpu
    real(rkind), allocatable, dimension(:), device :: detady_gpu, detady2_gpu, detadys_gpu
    real(rkind), allocatable, dimension(:), device :: dzitdz_gpu, dzitdz2_gpu, dzitdzs_gpu
    real(rkind), allocatable, dimension(:), device :: x_gpu, y_gpu, z_gpu, yn_gpu
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
  endtype base_gpu_object
!
!
contains
  subroutine alloc(self)
    class(base_gpu_object), intent(inout) :: self
!
    associate(nx => self%field%nx, ny => self%field%ny, nz => self%field%nz, ng => self%field%grid%ng, nv => self%field%nv)
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
    self%x_gpu = self%field%x
    self%y_gpu = self%field%y
    self%z_gpu = self%field%z
    self%yn_gpu = self%field%yn
!
    self%dcsidx_gpu = self%field%dcsidx
    self%dcsidxs_gpu = self%field%dcsidxs
    self%dcsidx2_gpu = self%field%dcsidx2
    self%detady_gpu = self%field%detady
    self%detadys_gpu = self%field%detadys
    self%detady2_gpu = self%field%detady2
    self%dzitdz_gpu = self%field%dzitdz
    self%dzitdzs_gpu = self%field%dzitdzs
    self%dzitdz2_gpu = self%field%dzitdz2
!
  endsubroutine alloc
!
  subroutine bcswap_step_1_cuf(nx, ny, nz, nv, ng, w_gpu, wbuf1s_gpu, wbuf2s_gpu, wbuf3s_gpu, wbuf4s_gpu, wbuf5s_gpu, wbuf6s_gpu)
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
    !@cuf iercuda=cudadevicesynchronize()
!
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
    !@cuf iercuda=cudadevicesynchronize()
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
    !@cuf iercuda=cudadevicesynchronize()
!
  endsubroutine bcswap_step_1_cuf
!
  subroutine bcswap_corner_step_1_cuf(nx, ny, nz, nv, ng, w_gpu, wbuf1s_c_gpu, wbuf2s_c_gpu, wbuf3s_c_gpu, wbuf4s_c_gpu)
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
            wbuf1s_c_gpu(m,i,j,k) = w_gpu(i,j,k,m)
            wbuf2s_c_gpu(m,i,j,k) = w_gpu(nx-ng+i,j,nz-ng+k,m)
            wbuf3s_c_gpu(m,i,j,k) = w_gpu(i,j,nz-ng+k,m)
            wbuf4s_c_gpu(m,i,j,k) = w_gpu(nx-ng+i,j,k,m)
          enddo
        enddo
      enddo
    enddo
    !@cuf iercuda=cudadevicesynchronize()
!
  endsubroutine bcswap_corner_step_1_cuf
!
  subroutine bcswap_step_3_cuf(nx, ny, nz, nv, ng, w_gpu, wbuf1r_gpu, wbuf2r_gpu, wbuf3r_gpu, wbuf4r&
  &_gpu, wbuf5r_gpu, wbuf6r_gpu, ileftx, ilefty, ileftz, irightx, irighty, irightz)
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
!
  subroutine bcswap_corner_step_3_cuf(nx, ny, nz, nv, ng, w_gpu, wbuf1r_c_gpu, wbuf2r_c_gpu, wbuf3r_&
  &c_gpu, wbuf4r_c_gpu, ileftbottom, ilefttop, irightbottom, irighttop)
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
!
  subroutine bcswap_corner(self, steps)
    class(base_gpu_object), intent(inout) :: self
    logical, dimension(3), optional :: steps
    logical, dimension(3) :: steps_
    integer :: iercuda, indc
    integer, dimension(mpi_status_size) :: istatus
!
    steps_ = .true. ; if(present(steps)) steps_ = steps
!
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv, w_gpu => se&
    &lf%w_gpu, wbuf1s_c_gpu => self%wbuf1s_c_gpu, wbuf2s_c_gpu => self%wbuf2s_c_gpu, wbuf3s_c_gpu => self&
    &%wbuf3s_c_gpu, wbuf4s_c_gpu => self%wbuf4s_c_gpu, wbuf1r_c_gpu => self%wbuf1r_c_gpu, wbuf2r_c_gpu =>&
    & self%wbuf2r_c_gpu, wbuf3r_c_gpu => self%wbuf3r_c_gpu, wbuf4r_c_gpu => self%wbuf4r_c_gpu, ileftbotto&
    &m => self%field%ileftbottom, irightbottom => self%field%irightbottom, ilefttop => self%field%ileftto&
    &p, irighttop => self%field%irighttop, mp_cart => self%field%mp_cart, iermpi => self%mpi_err, nrank =&
    &> self%field%myrank)
!
      if(steps_(1)) then
        call bcswap_corner_step_1_cuf(nx, ny, nz, nv, ng, w_gpu, wbuf1s_c_gpu, wbuf2s_c_gpu, wbuf3s_c_gpu, wbuf4s_c_gpu)
      endif
      if(steps_(2)) then
        indc = nv*ny*ng*ng
        if (ileftbottom == nrank) then
          wbuf2r_c_gpu = wbuf1s_c_gpu
          wbuf1r_c_gpu = wbuf2s_c_gpu
        else
          call mpi_sendrecv(wbuf1s_c_gpu,indc,mpi_prec,ileftbottom ,1, wbuf2r_c_gpu,indc,mpi_prec,ir&
          &ighttop ,1,mp_cart,istatus,iermpi)
          call mpi_sendrecv(wbuf2s_c_gpu,indc,mpi_prec,irighttop ,2, wbuf1r_c_gpu,indc,mpi_prec,ilef&
          &tbottom ,2,mp_cart,istatus,iermpi)
        endif
        if (ilefttop == nrank) then
          wbuf4r_c_gpu = wbuf3s_c_gpu
          wbuf3r_c_gpu = wbuf4s_c_gpu
        else
          call mpi_sendrecv(wbuf3s_c_gpu,indc,mpi_prec,ilefttop ,3, wbuf4r_c_gpu,indc,mpi_prec,irigh&
          &tbottom,3,mp_cart,istatus,iermpi)
          call mpi_sendrecv(wbuf4s_c_gpu,indc,mpi_prec,irightbottom,4, wbuf3r_c_gpu,indc,mpi_prec,il&
          &efttop ,4,mp_cart,istatus,iermpi)
        endif
      endif
      if(steps_(3)) then
        call bcswap_corner_step_3_cuf(nx, ny, nz, nv, ng, w_gpu, wbuf1r_c_gpu, wbuf2r_c_gpu, wbuf3r_&
        &c_gpu, wbuf4r_c_gpu, ileftbottom, ilefttop, irightbottom, irighttop)
      endif
    endassociate
  endsubroutine bcswap_corner
!
  subroutine bcswap_corner_var(self, w_swap, steps)
    class(base_gpu_object), intent(inout) :: self
    real(rkind), dimension(1-self%ng:self%nx+self%ng, 1-self%ng:self%ny+self%ng, 1-self%ng:self%nz+self%ng,1:1), device :: w_swap
    logical, dimension(3), optional :: steps
    logical, dimension(3) :: steps_
    integer :: iercuda, indc
    integer, dimension(mpi_status_size) :: istatus
!
    steps_ = .true. ; if(present(steps)) steps_ = steps
!
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv, wbuf1s_c_gp&
    &u => self%wbuf1s_c_gpu, wbuf2s_c_gpu => self%wbuf2s_c_gpu, wbuf3s_c_gpu => self%wbuf3s_c_gpu, wbuf4s&
    &_c_gpu => self%wbuf4s_c_gpu, wbuf1r_c_gpu => self%wbuf1r_c_gpu, wbuf2r_c_gpu => self%wbuf2r_c_gpu, w&
    &buf3r_c_gpu => self%wbuf3r_c_gpu, wbuf4r_c_gpu => self%wbuf4r_c_gpu, ileftbottom => self%field%ileft&
    &bottom, irightbottom => self%field%irightbottom, ilefttop => self%field%ilefttop, irighttop => self%&
    &field%irighttop, mp_cart => self%field%mp_cart, iermpi => self%mpi_err, nrank => self%field%myrank)
!
      if(steps_(1)) then
        call bcswap_corner_step_1_cuf(nx, ny, nz, 1, ng, w_swap, wbuf1s_c_gpu, wbuf2s_c_gpu, wbuf3s_c_gpu, wbuf4s_c_gpu)
      endif
      if(steps_(2)) then
        indc = ny*ng*ng
        if (ileftbottom == nrank) then
          wbuf2r_c_gpu = wbuf1s_c_gpu
          wbuf1r_c_gpu = wbuf2s_c_gpu
        else
          call mpi_sendrecv(wbuf1s_c_gpu,indc,mpi_prec,ileftbottom ,1, wbuf2r_c_gpu,indc,mpi_prec,ir&
          &ighttop ,1,mp_cart,istatus,iermpi)
          call mpi_sendrecv(wbuf2s_c_gpu,indc,mpi_prec,irighttop ,2, wbuf1r_c_gpu,indc,mpi_prec,ilef&
          &tbottom ,2,mp_cart,istatus,iermpi)
        endif
        if (ilefttop == nrank) then
          wbuf4r_c_gpu = wbuf3s_c_gpu
          wbuf3r_c_gpu = wbuf4s_c_gpu
        else
          call mpi_sendrecv(wbuf3s_c_gpu,indc,mpi_prec,ilefttop ,3, wbuf4r_c_gpu,indc,mpi_prec,irigh&
          &tbottom,3,mp_cart,istatus,iermpi)
          call mpi_sendrecv(wbuf4s_c_gpu,indc,mpi_prec,irightbottom,4, wbuf3r_c_gpu,indc,mpi_prec,il&
          &efttop ,4,mp_cart,istatus,iermpi)
        endif
      endif
      if(steps_(3)) then
        call bcswap_corner_step_3_cuf(nx, ny, nz, 1, ng, w_swap, wbuf1r_c_gpu, wbuf2r_c_gpu, wbuf3r_&
        &c_gpu, wbuf4r_c_gpu, ileftbottom, ilefttop, irightbottom, irighttop)
      endif
    endassociate
  endsubroutine bcswap_corner_var
!
  subroutine bcswap(self, steps)
    class(base_gpu_object), intent(inout) :: self
    logical, dimension(3), optional :: steps
    logical, dimension(3) :: steps_
    integer :: iercuda, indx, indy, indz
    integer, dimension(mpi_status_size) :: istatus
!
    steps_ = .true. ; if(present(steps)) steps_ = steps
!
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv, w_gpu => se&
    &lf%w_gpu, wbuf1s_gpu => self%wbuf1s_gpu, wbuf2s_gpu => self%wbuf2s_gpu, wbuf3s_gpu => self%wbuf3s_gp&
    &u, wbuf4s_gpu => self%wbuf4s_gpu, wbuf5s_gpu => self%wbuf5s_gpu, wbuf6s_gpu => self%wbuf6s_gpu, wbuf&
    &1r_gpu => self%wbuf1r_gpu, wbuf2r_gpu => self%wbuf2r_gpu, wbuf3r_gpu => self%wbuf3r_gpu, wbuf4r_gpu &
    &=> self%wbuf4r_gpu, wbuf5r_gpu => self%wbuf5r_gpu, wbuf6r_gpu => self%wbuf6r_gpu, ileftx => self%fie&
    &ld%ileftx,irightx => self%field%irightx, nrank_x => self%field%nrank_x, ilefty => self%field%ilefty,&
    &irighty => self%field%irighty, nrank_y => self%field%nrank_y, ileftz => self%field%ileftz,irightz =>&
    & self%field%irightz, nrank_z => self%field%nrank_z, mp_cartx => self%field%mp_cartx, mp_carty => sel&
    &f%field%mp_carty, mp_cartz => self%field%mp_cartz)
!
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
        call bcswap_step_3_cuf(nx, ny, nz, nv, ng, w_gpu, wbuf1r_gpu, wbuf2r_gpu, wbuf3r_gpu, wbuf4r&
        &_gpu, wbuf5r_gpu, wbuf6r_gpu, ileftx, ilefty, ileftz, irightx, irighty, irightz)
      endif
    endassociate
  endsubroutine bcswap
!
  subroutine bcswap_var(self, w_swap, steps)
    class(base_gpu_object), intent(inout) :: self
    real(rkind), dimension(1-self%ng:self%nx+self%ng, 1-self%ng:self%ny+self%ng, 1-self%ng:self%nz+self%ng,1:1), device :: w_swap
    logical, dimension(3), optional :: steps
    logical, dimension(3) :: steps_
    integer :: iercuda, indx, indy, indz
    integer, dimension(mpi_status_size) :: istatus
!
    steps_ = .true. ; if(present(steps)) steps_ = steps
!
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, wbuf1s_gpu => self%wbuf1s_&
    &gpu, wbuf2s_gpu => self%wbuf2s_gpu, wbuf3s_gpu => self%wbuf3s_gpu, wbuf4s_gpu => self%wbuf4s_gpu, wb&
    &uf5s_gpu => self%wbuf5s_gpu, wbuf6s_gpu => self%wbuf6s_gpu, wbuf1r_gpu => self%wbuf1r_gpu, wbuf2r_gp&
    &u => self%wbuf2r_gpu, wbuf3r_gpu => self%wbuf3r_gpu, wbuf4r_gpu => self%wbuf4r_gpu, wbuf5r_gpu => se&
    &lf%wbuf5r_gpu, wbuf6r_gpu => self%wbuf6r_gpu, ileftx => self%field%ileftx,irightx => self%field%irig&
    &htx, nrank_x => self%field%nrank_x, ilefty => self%field%ilefty,irighty => self%field%irighty, nrank&
    &_y => self%field%nrank_y, ileftz => self%field%ileftz,irightz => self%field%irightz, nrank_z => self&
    &%field%nrank_z, mp_cartx => self%field%mp_cartx, mp_carty => self%field%mp_carty, mp_cartz => self%f&
    &ield%mp_cartz)
!
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
        call bcswap_step_3_cuf(nx, ny, nz, 1, ng, w_swap, wbuf1r_gpu, wbuf2r_gpu, wbuf3r_gpu, wbuf4r&
        &_gpu, wbuf5r_gpu, wbuf6r_gpu, ileftx, ilefty, ileftz, irightx, irighty, irightz)
      endif
    endassociate
  endsubroutine bcswap_var
!
  subroutine copy_cpu_gpu(self)
    class(base_gpu_object), intent(inout) :: self
    integer :: i,j,k,iv
!
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
!
    self%w_gpu = self%w_t
!
  endsubroutine copy_cpu_gpu
!
!
  subroutine copy_gpu_cpu(self)
    class(base_gpu_object), intent(inout) :: self
    integer :: i,j,k,iv
!
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
    class(base_gpu_object), intent(inout) :: self
    class(field_object), target :: field
    type(cudadeviceprop) :: device_properties
    real(rkind) :: device_mem_avail
!
    self%field => field
    self%nx = self%field%nx
    self%ny = self%field%ny
    self%nz = self%field%nz
    self%ng = self%field%grid%ng
    self%nv = self%field%nv
!
    call get_mpi_basic_info(self%nprocs, self%myrank, self%masterproc, self%mpi_err)
!
    call self%field%check_cpu_mem(description="--base-initialization--")
!
    call mpi_comm_split_type(mpi_comm_world, mpi_comm_type_shared, 0, mpi_info_null, self%local_comm, self%mpi_err)
    call mpi_comm_rank(self%local_comm, self%mydev, self%mpi_err)
    self%mpi_err = cudasetdevice(self%mydev)
    self%mpi_err = cudagetdeviceproperties(device_properties, self%mydev)
!
    device_mem_avail = real(device_properties%totalglobalmem, rkind)/(1024_rkind**3)
!
    call self%alloc()
!
  endsubroutine initialize
!
  subroutine check_gpu_mem(self, description)
    class(base_gpu_object), intent(inout) :: self
    character(*) :: description
    integer :: ierr
    integer(cuda_count_kind) :: mem_free, mem_total
    character(128) :: proc_name
    integer :: resultlen
    call mpi_get_processor_name(proc_name, resultlen, ierr)
    ierr = cudamemgetinfo(mem_free, mem_total)
    write(error_unit, "(a,2x,a,2x,a,2x,i0,2x,i0,2x,i0)") 'gpu rank,mems: ', description,proc_name(1:&
    &resultlen),self%myrank, mem_free, mem_total
  endsubroutine check_gpu_mem
!
endmodule streams_base_gpu_object

