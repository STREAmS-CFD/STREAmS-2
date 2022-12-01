module streams_base_cpu_object
    !< STREAmS, base CPU class definition

    use streams_field_object, only : field_object
    use streams_parameters
    use MPI

    implicit none
    private
    public :: base_cpu_object

    type :: base_cpu_object
        type(field_object), pointer :: field=>null() 
        ! Replica of field and grid sizes
        integer :: nx, ny, nz, ng, nv
        ! MPI data
        integer(ikind)              :: myrank=0_ikind       !< MPI rank process.
        integer(ikind)              :: nprocs=1_ikind       !< Number of MPI processes.
        logical                     :: masterproc  
        integer(ikind)              :: mpi_err=0_ikind      !< Error traping flag.
        integer(ikind)              :: local_comm=0_ikind   !< Local communicator.
        ! CPU data
        real(rkind), allocatable, dimension(:,:,:,:) :: w_cpu
        real(rkind), allocatable, dimension(:,:,:,:)         :: w_t

        real(rkind), dimension(:,:,:,:), allocatable :: wbuf1s_cpu, wbuf2s_cpu, wbuf3s_cpu, &
                                                                wbuf4s_cpu, wbuf5s_cpu, wbuf6s_cpu
        real(rkind), dimension(:,:,:,:), allocatable :: wbuf1r_cpu, wbuf2r_cpu, wbuf3r_cpu, &
                                                                wbuf4r_cpu, wbuf5r_cpu, wbuf6r_cpu
        real(rkind), dimension(:,:,:,:), allocatable :: wbuf1s_c_cpu, wbuf2s_c_cpu, wbuf3s_c_cpu, wbuf4s_c_cpu
        real(rkind), dimension(:,:,:,:), allocatable :: wbuf1r_c_cpu, wbuf2r_c_cpu, wbuf3r_c_cpu, wbuf4r_c_cpu
        real(rkind), allocatable, dimension(:) :: dcsidx_cpu, dcsidx2_cpu, dcsidxs_cpu
        real(rkind), allocatable, dimension(:) :: detady_cpu, detady2_cpu, detadys_cpu
        real(rkind), allocatable, dimension(:) :: dzitdz_cpu, dzitdz2_cpu, dzitdzs_cpu
        real(rkind), allocatable, dimension(:) :: x_cpu, y_cpu, z_cpu, yn_cpu
    contains
        ! public methods
        procedure, pass(self) :: alloc        
        procedure, pass(self) :: copy_from_field 
        procedure, pass(self) :: copy_to_field 
        procedure, pass(self) :: initialize   
        procedure, pass(self) :: bcswap       
        procedure, pass(self) :: bcswap_var       
        procedure, pass(self) :: bcswap_corner       
        procedure, pass(self) :: bcswap_corner_var       
            endtype base_cpu_object

    !interface assign_allocatable_cpu
    !    !< Safe assign allocatable arrays (CPU), generic interface.
    !    module procedure assign_allocatable_INT64_1D_cpu !< Safe assign allocatable arrays, I8P 1D type.
    !    module procedure assign_allocatable_INT64_2D_cpu !< Safe assign allocatable arrays, I8P 2D type.
    !    module procedure assign_allocatable_INT32_1D_cpu !< Safe assign allocatable arrays, ikind 1D type.
    !endinterface assign_allocatable_cpu

contains
    ! public methods
    subroutine alloc(self)
    class(base_cpu_object), intent(inout)        :: self     !< The base backend.

        associate(nx => self%field%nx, ny => self%field%ny, nz => self%field%nz, &
                  ng => self%field%grid%ng, nv => self%field%nv)

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

        endassociate

        self%x_cpu  = self%field%x
        self%y_cpu  = self%field%y
        self%z_cpu  = self%field%z
        self%yn_cpu = self%field%yn

        self%dcsidx_cpu  = self%field%dcsidx
        self%dcsidxs_cpu = self%field%dcsidxs
        self%dcsidx2_cpu = self%field%dcsidx2
        self%detady_cpu  = self%field%detady
        self%detadys_cpu = self%field%detadys
        self%detady2_cpu = self%field%detady2
        self%dzitdz_cpu  = self%field%dzitdz
        self%dzitdzs_cpu = self%field%dzitdzs
        self%dzitdz2_cpu = self%field%dzitdz2

    endsubroutine alloc

    subroutine bcswap_step_1_cpu(nx, ny, nz, nv, ng, w_cpu, &
                    wbuf1s_cpu, wbuf2s_cpu, wbuf3s_cpu, wbuf4s_cpu, wbuf5s_cpu, wbuf6s_cpu)
      integer :: nx, ny, nz, ng, nv
      real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: w_cpu
      real(rkind), dimension(1:,1:,1:,1:) :: wbuf1s_cpu, wbuf2s_cpu, wbuf3s_cpu, wbuf4s_cpu, wbuf5s_cpu, wbuf6s_cpu
      integer :: i,j,k,m

      do k=1,nz
       do j=1,ny
        do i=1,ng
         do m=1,nv
          wbuf1s_cpu(i,j,k,m) = w_cpu(i,j,k,m)
          wbuf2s_cpu(i,j,k,m) = w_cpu(nx-ng+i,j,k,m)
         enddo
        enddo
       enddo
      enddo

      ! CUDA check (debugging purposes)
      !iercuda = cudaGetLastError()
      !if (iercuda /= cudaSuccess) then
      ! print*,"CUDA ERROR! ",cudaGetErrorString(iercuda)
      !endif

      do k=1,nz
       do j=1,ng
        do i=1,nx
         do m=1,nv
          wbuf3s_cpu(i,j,k,m) = w_cpu(i,j,k,m)
          wbuf4s_cpu(i,j,k,m) = w_cpu(i,ny-ng+j,k,m)
         enddo
        enddo
       enddo
      enddo
  
     do k=1,ng
      do j=1,ny
       do i=1,nx
        do m=1,nv
         wbuf5s_cpu(i,j,k,m) = w_cpu(i,j,k,m)
         wbuf6s_cpu(i,j,k,m) = w_cpu(i,j,nz-ng+k,m)
        enddo
       enddo
      enddo
     enddo
      
    endsubroutine bcswap_step_1_cpu

    subroutine bcswap_corner_step_1_cpu(nx, ny, nz, nv, ng, w_cpu, &
                    wbuf1s_c_cpu, wbuf2s_c_cpu, wbuf3s_c_cpu, wbuf4s_c_cpu)
      integer :: nx, ny, nz, ng, nv
      real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: w_cpu
      real(rkind), dimension(1:,1:,1:,1:) :: wbuf1s_c_cpu, wbuf2s_c_cpu, wbuf3s_c_cpu, wbuf4s_c_cpu
      integer :: i,j,k,m

       do j=1,ny
        do m=1,nv
         do k=1,ng
          do i=1,ng
           wbuf1s_c_cpu(m,i,j,k) = w_cpu(i,j,k,m)              ! ileftbottom
           wbuf2s_c_cpu(m,i,j,k) = w_cpu(nx-ng+i,j,nz-ng+k,m)  ! irighttop
           wbuf3s_c_cpu(m,i,j,k) = w_cpu(i,j,nz-ng+k,m)        ! ilefttop
           wbuf4s_c_cpu(m,i,j,k) = w_cpu(nx-ng+i,j,k,m)        ! irightbottom
          enddo
         enddo
        enddo
       enddo
      
    endsubroutine bcswap_corner_step_1_cpu

    subroutine bcswap_step_3_cpu(nx, ny, nz, nv, ng, w_cpu, &
                    wbuf1r_cpu, wbuf2r_cpu, wbuf3r_cpu, wbuf4r_cpu, wbuf5r_cpu, wbuf6r_cpu, &
                    ileftx, ilefty, ileftz, irightx, irighty, irightz)
      integer :: nx, ny, nz, ng, nv
      real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: w_cpu
      real(rkind), dimension(1:,1:,1:,1:) :: wbuf1r_cpu, wbuf2r_cpu, wbuf3r_cpu, wbuf4r_cpu, wbuf5r_cpu, wbuf6r_cpu
      integer :: i,j,k,m
      integer ::  ileftx, ilefty, ileftz, irightx, irighty, irightz
       if (ileftx/=mpi_proc_null) then
       do k=1,nz
        do j=1,ny
         do i=1,ng
          do m=1,nv
           w_cpu(i-ng,j,k,m) = wbuf1r_cpu(i,j,k,m)
          enddo
         enddo
        enddo
       enddo
      endif
      if (irightx/=mpi_proc_null) then
       do k=1,nz
        do j=1,ny
         do i=1,ng
          do m=1,nv
           w_cpu(nx+i,j,k,m) = wbuf2r_cpu(i,j,k,m)
          enddo
         enddo
        enddo
       enddo
      endif
      if (ilefty/=mpi_proc_null) then
       do k=1,nz
        do j=1,ng
         do i=1,nx
          do m=1,nv
           w_cpu(i,j-ng,k,m) = wbuf3r_cpu(i,j,k,m)
          enddo
         enddo
        enddo
       enddo
      endif
      if (irighty/=mpi_proc_null) then
       do k=1,nz
        do j=1,ng
         do i=1,nx
          do m=1,nv
           w_cpu(i,ny+j,k,m) = wbuf4r_cpu(i,j,k,m)
          enddo
         enddo
        enddo
       enddo
      endif
      if (ileftz/=mpi_proc_null) then
       do k=1,ng
        do j=1,ny
         do i=1,nx
          do m=1,nv
           w_cpu(i,j,k-ng,m) = wbuf5r_cpu(i,j,k,m)
          enddo
         enddo
        enddo
       enddo
      endif
      if (irightz/=mpi_proc_null) then
       do k=1,ng
        do j=1,ny
         do i=1,nx
          do m=1,nv
           w_cpu(i,j,nz+k,m) = wbuf6r_cpu(i,j,k,m)
          enddo
         enddo
        enddo
       enddo
      endif
    endsubroutine bcswap_step_3_cpu

    subroutine bcswap_corner_step_3_cpu(nx, ny, nz, nv, ng, w_cpu, &
                    wbuf1r_c_cpu, wbuf2r_c_cpu, wbuf3r_c_cpu, wbuf4r_c_cpu, &
                    ileftbottom, ilefttop, irightbottom, irighttop)
      integer :: nx, ny, nz, ng, nv
      real(rkind), dimension(1-ng:,1-ng:,1-ng:,1:) :: w_cpu
      real(rkind), dimension(1:,1:,1:,1:) :: wbuf1r_c_cpu, wbuf2r_c_cpu, wbuf3r_c_cpu, wbuf4r_c_cpu
      integer :: i,j,k,m
      integer :: ileftbottom, ilefttop, irightbottom, irighttop

      if (ileftbottom/=mpi_proc_null) then
       do j=1,ny
        do m=1,nv
         do k=1,ng
          do i=1,ng
           w_cpu(i-ng,j,k-ng,m) = wbuf1r_c_cpu(m,i,j,k)
          enddo
         enddo
        enddo
       enddo
      endif
      if (irighttop/=mpi_proc_null) then
       do j=1,ny
        do m=1,nv
         do k=1,ng
          do i=1,ng
           w_cpu(i+nx,j,k+nz,m) = wbuf2r_c_cpu(m,i,j,k)
          enddo
         enddo
        enddo
       enddo
      endif
      if (ilefttop/=mpi_proc_null) then
       do j=1,ny
        do m=1,nv
         do k=1,ng
          do i=1,ng
           w_cpu(i-ng,j,k+nz,m) = wbuf3r_c_cpu(m,i,j,k)
          enddo
         enddo
        enddo
       enddo
      endif
      if (irightbottom/=mpi_proc_null) then
       do j=1,ny
        do m=1,nv
         do k=1,ng
          do i=1,ng
           w_cpu(i+nx,j,k-ng,m) = wbuf4r_c_cpu(m,i,j,k)
          enddo
         enddo
        enddo
       enddo
      endif
    endsubroutine bcswap_corner_step_3_cpu

    subroutine bcswap_corner(self, steps)
     class(base_cpu_object), intent(inout) :: self !< The base backend.
     logical, dimension(3), optional :: steps
     logical, dimension(3) :: steps_
     integer :: iercuda, indc
     integer, dimension(mpi_status_size) :: istatus

     steps_ = .true. ; if(present(steps)) steps_ = steps

     associate(nx => self%nx, ny => self%ny, nz => self%nz, &
               ng => self%ng, nv => self%nv, &
               w_cpu => self%w_cpu, &
               wbuf1s_c_cpu => self%wbuf1s_c_cpu, wbuf2s_c_cpu => self%wbuf2s_c_cpu,            &
               wbuf3s_c_cpu => self%wbuf3s_c_cpu, wbuf4s_c_cpu => self%wbuf4s_c_cpu,            &
               wbuf1r_c_cpu => self%wbuf1r_c_cpu, wbuf2r_c_cpu => self%wbuf2r_c_cpu,            &
               wbuf3r_c_cpu => self%wbuf3r_c_cpu, wbuf4r_c_cpu => self%wbuf4r_c_cpu,            &
               ileftbottom  => self%field%ileftbottom, irightbottom => self%field%irightbottom, &
               ilefttop     => self%field%ilefttop,    irighttop => self%field%irighttop,       &
               mp_cart      => self%field%mp_cart, iermpi => self%mpi_err, nrank => self%field%myrank)

     if(steps_(1)) then
      call bcswap_corner_step_1_cpu(nx, ny, nz, nv, ng, w_cpu, &
             wbuf1s_c_cpu, wbuf2s_c_cpu, wbuf3s_c_cpu, wbuf4s_c_cpu)
     endif
     if(steps_(2)) then
      indc = nv*ny*ng*ng
      if (ileftbottom == nrank) then
          wbuf2r_c_cpu = wbuf1s_c_cpu
          wbuf1r_c_cpu = wbuf2s_c_cpu
      else
          call mpi_sendrecv(wbuf1s_c_cpu,indc,mpi_prec,ileftbottom ,1, &
                            wbuf2r_c_cpu,indc,mpi_prec,irighttop   ,1,mp_cart,istatus,iermpi)
          call mpi_sendrecv(wbuf2s_c_cpu,indc,mpi_prec,irighttop   ,2, &
                            wbuf1r_c_cpu,indc,mpi_prec,ileftbottom ,2,mp_cart,istatus,iermpi)
      endif
      if (ilefttop == nrank) then
          wbuf4r_c_cpu = wbuf3s_c_cpu
          wbuf3r_c_cpu = wbuf4s_c_cpu
      else
          call mpi_sendrecv(wbuf3s_c_cpu,indc,mpi_prec,ilefttop    ,3, &
                            wbuf4r_c_cpu,indc,mpi_prec,irightbottom,3,mp_cart,istatus,iermpi)
          call mpi_sendrecv(wbuf4s_c_cpu,indc,mpi_prec,irightbottom,4, &
                            wbuf3r_c_cpu,indc,mpi_prec,ilefttop    ,4,mp_cart,istatus,iermpi)
      endif
     endif
     if(steps_(3)) then
      call bcswap_corner_step_3_cpu(nx, ny, nz, nv, ng, w_cpu, &
          wbuf1r_c_cpu, wbuf2r_c_cpu, wbuf3r_c_cpu, wbuf4r_c_cpu, &
          ileftbottom, ilefttop, irightbottom, irighttop)
     endif
     endassociate
    endsubroutine bcswap_corner

    subroutine bcswap_corner_var(self, w_swap, steps)
     class(base_cpu_object), intent(inout) :: self !< The base backend.
     real(rkind), dimension(1-self%ng:self%nx+self%ng, &
                            1-self%ng:self%ny+self%ng, &
                            1-self%ng:self%nz+self%ng,1:1) :: w_swap
     logical, dimension(3), optional :: steps
     logical, dimension(3) :: steps_
     integer :: iercuda, indc
     integer, dimension(mpi_status_size) :: istatus

     steps_ = .true. ; if(present(steps)) steps_ = steps

     associate(nx => self%nx, ny => self%ny, nz => self%nz, &
               ng => self%ng, nv => self%nv, &
               wbuf1s_c_cpu => self%wbuf1s_c_cpu, wbuf2s_c_cpu => self%wbuf2s_c_cpu,            &
               wbuf3s_c_cpu => self%wbuf3s_c_cpu, wbuf4s_c_cpu => self%wbuf4s_c_cpu,            &
               wbuf1r_c_cpu => self%wbuf1r_c_cpu, wbuf2r_c_cpu => self%wbuf2r_c_cpu,            &
               wbuf3r_c_cpu => self%wbuf3r_c_cpu, wbuf4r_c_cpu => self%wbuf4r_c_cpu,            &
               ileftbottom  => self%field%ileftbottom, irightbottom => self%field%irightbottom, &
               ilefttop     => self%field%ilefttop,    irighttop => self%field%irighttop,       &
               mp_cart      => self%field%mp_cart, iermpi => self%mpi_err, nrank => self%field%myrank)

     if(steps_(1)) then
      call bcswap_corner_step_1_cpu(nx, ny, nz, 1, ng, w_swap, &
             wbuf1s_c_cpu, wbuf2s_c_cpu, wbuf3s_c_cpu, wbuf4s_c_cpu)
     endif
     if(steps_(2)) then
      indc = ny*ng*ng
      if (ileftbottom == nrank) then
          wbuf2r_c_cpu = wbuf1s_c_cpu
          wbuf1r_c_cpu = wbuf2s_c_cpu
      else
          call mpi_sendrecv(wbuf1s_c_cpu,indc,mpi_prec,ileftbottom ,1, &
                            wbuf2r_c_cpu,indc,mpi_prec,irighttop   ,1,mp_cart,istatus,iermpi)
          call mpi_sendrecv(wbuf2s_c_cpu,indc,mpi_prec,irighttop   ,2, &
                            wbuf1r_c_cpu,indc,mpi_prec,ileftbottom ,2,mp_cart,istatus,iermpi)
      endif
      if (ilefttop == nrank) then
          wbuf4r_c_cpu = wbuf3s_c_cpu
          wbuf3r_c_cpu = wbuf4s_c_cpu
      else
          call mpi_sendrecv(wbuf3s_c_cpu,indc,mpi_prec,ilefttop    ,3, &
                            wbuf4r_c_cpu,indc,mpi_prec,irightbottom,3,mp_cart,istatus,iermpi)
          call mpi_sendrecv(wbuf4s_c_cpu,indc,mpi_prec,irightbottom,4, &
                            wbuf3r_c_cpu,indc,mpi_prec,ilefttop    ,4,mp_cart,istatus,iermpi)
      endif
     endif
     if(steps_(3)) then
      call bcswap_corner_step_3_cpu(nx, ny, nz, 1, ng, w_swap, &
          wbuf1r_c_cpu, wbuf2r_c_cpu, wbuf3r_c_cpu, wbuf4r_c_cpu, &
          ileftbottom, ilefttop, irightbottom, irighttop)
     endif
     endassociate
    endsubroutine bcswap_corner_var

    subroutine bcswap(self, steps)
     class(base_cpu_object), intent(inout) :: self !< The base backend.
     logical, dimension(3), optional :: steps
     logical, dimension(3) :: steps_
     integer :: iercuda, indx, indy, indz
     integer, dimension(mpi_status_size) :: istatus

     steps_ = .true. ; if(present(steps)) steps_ = steps

     associate(nx => self%nx, ny => self%ny, nz => self%nz, &
               ng => self%ng, nv => self%nv, &
               w_cpu => self%w_cpu, &
               wbuf1s_cpu => self%wbuf1s_cpu, wbuf2s_cpu => self%wbuf2s_cpu, &
               wbuf3s_cpu => self%wbuf3s_cpu, wbuf4s_cpu => self%wbuf4s_cpu, &
               wbuf5s_cpu => self%wbuf5s_cpu, wbuf6s_cpu => self%wbuf6s_cpu, &
               wbuf1r_cpu => self%wbuf1r_cpu, wbuf2r_cpu => self%wbuf2r_cpu, &
               wbuf3r_cpu => self%wbuf3r_cpu, wbuf4r_cpu => self%wbuf4r_cpu, &
               wbuf5r_cpu => self%wbuf5r_cpu, wbuf6r_cpu => self%wbuf6r_cpu, &
               ileftx => self%field%ileftx,irightx => self%field%irightx, nrank_x => self%field%nrank_x, &
               ilefty => self%field%ilefty,irighty => self%field%irighty, nrank_y => self%field%nrank_y, &
               ileftz => self%field%ileftz,irightz => self%field%irightz, nrank_z => self%field%nrank_z, &
               mp_cartx => self%field%mp_cartx, mp_carty => self%field%mp_carty, mp_cartz => self%field%mp_cartz)

     if(steps_(1)) then
      call bcswap_step_1_cpu(nx, ny, nz, nv, ng, w_cpu, &
             wbuf1s_cpu, wbuf2s_cpu, wbuf3s_cpu, wbuf4s_cpu, wbuf5s_cpu, wbuf6s_cpu)
     endif
     if(steps_(2)) then
      indx = nv*ng*ny*nz
      indy = nv*nx*ng*nz
      indz = nv*nx*ny*ng
      ! http://developer.download.nvidia.com/compute/cuda/2_3/toolkit/...
      !   docs/online/group__CUDART__MEMORY_ge4366f68c6fa8c85141448f187d2aa13.html
      ! IMPORTANT NOTE: Copies with kind == cudaMemcpyDeviceToDevice are asynchronous
      ! with respect to the host, but never overlap with kernel execution
      ! For this reason here simple copies are used and not cudaMemcpyAsync D2D
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
      call bcswap_step_3_cpu(nx, ny, nz, nv, ng, w_cpu, &
          wbuf1r_cpu, wbuf2r_cpu, wbuf3r_cpu, wbuf4r_cpu, wbuf5r_cpu, wbuf6r_cpu, &
          ileftx, ilefty, ileftz, irightx, irighty, irightz)
     endif
     endassociate
    endsubroutine bcswap

    subroutine bcswap_var(self, w_swap, steps)
     class(base_cpu_object), intent(inout) :: self !< The base backend.
     real(rkind), dimension(1-self%ng:self%nx+self%ng, &
                            1-self%ng:self%ny+self%ng, &
                            1-self%ng:self%nz+self%ng,1:1) :: w_swap
     logical, dimension(3), optional :: steps
     logical, dimension(3) :: steps_
     integer :: iercuda, indx, indy, indz
     integer, dimension(mpi_status_size) :: istatus

     steps_ = .true. ; if(present(steps)) steps_ = steps

     associate(nx => self%nx, ny => self%ny, nz => self%nz, &
               ng => self%ng, &
               wbuf1s_cpu => self%wbuf1s_cpu, wbuf2s_cpu => self%wbuf2s_cpu, &
               wbuf3s_cpu => self%wbuf3s_cpu, wbuf4s_cpu => self%wbuf4s_cpu, &
               wbuf5s_cpu => self%wbuf5s_cpu, wbuf6s_cpu => self%wbuf6s_cpu, &
               wbuf1r_cpu => self%wbuf1r_cpu, wbuf2r_cpu => self%wbuf2r_cpu, &
               wbuf3r_cpu => self%wbuf3r_cpu, wbuf4r_cpu => self%wbuf4r_cpu, &
               wbuf5r_cpu => self%wbuf5r_cpu, wbuf6r_cpu => self%wbuf6r_cpu, &
               ileftx => self%field%ileftx,irightx => self%field%irightx, nrank_x => self%field%nrank_x, &
               ilefty => self%field%ilefty,irighty => self%field%irighty, nrank_y => self%field%nrank_y, &
               ileftz => self%field%ileftz,irightz => self%field%irightz, nrank_z => self%field%nrank_z, &
               mp_cartx => self%field%mp_cartx, mp_carty => self%field%mp_carty, mp_cartz => self%field%mp_cartz)

     if(steps_(1)) then
      call bcswap_step_1_cpu(nx, ny, nz, 1, ng, w_swap, &
             wbuf1s_cpu, wbuf2s_cpu, wbuf3s_cpu, wbuf4s_cpu, wbuf5s_cpu, wbuf6s_cpu)
     endif
     if(steps_(2)) then
      indx = ng*ny*nz
      indy = nx*ng*nz
      indz = nx*ny*ng
      ! http://developer.download.nvidia.com/compute/cuda/2_3/toolkit...
      !   /docs/online/group__CUDART__MEMORY_ge4366f68c6fa8c85141448f187d2aa13.html
      ! IMPORTANT NOTE: Copies with kind == cudaMemcpyDeviceToDevice are asynchronous
      ! with respect to the host, but never overlap with kernel execution
      ! For this reason here simple copies are used and not cudaMemcpyAsync D2D
      if(ileftx == nrank_x) then
          wbuf2r_cpu = wbuf1s_cpu ! inefficient, only 1 variable here to be copied
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
      call bcswap_step_3_cpu(nx, ny, nz, 1, ng, w_swap, &
          wbuf1r_cpu, wbuf2r_cpu, wbuf3r_cpu, wbuf4r_cpu, wbuf5r_cpu, wbuf6r_cpu, &
          ileftx, ilefty, ileftz, irightx, irighty, irightz)
     endif
     endassociate
    endsubroutine bcswap_var

    subroutine copy_from_field(self)
     !< Copy data from CPU to field.

     class(base_cpu_object), intent(inout) :: self !< The base backend.
     integer :: i,j,k,iv

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

     ! copy cpu to field

     self%w_cpu = self%w_t

    endsubroutine copy_from_field

    !subroutine copy_to_field_var(self, w_var_cpu, w_var)
    ! class(base_cpu_object), intent(inout) :: self
    ! real(rkind), dimension(1-self%ng:, 1-self%ng, 1-self%ng, :) :: w_io_cpu
    ! ! copy field to cpu

    ! self%w_var_t = self%w_var_cpu 
    ! associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng)
    !  do k=1-ng,nz+ng
    !   do j=1-ng,ny+ng
    !    do i=1-ng,nx+ng
    !     w_var(i,j,k) = self%w_var_t(i,j,k,iv)
    !    enddo
    !   enddo
    !  enddo
    ! endassociate
    !endsubroutine copy_to_field_var

    subroutine copy_to_field(self)
     !< Copy data from CPU to field.

     class(base_cpu_object), intent(inout) :: self !< The base backend.
     integer :: i,j,k,iv

     ! copy field to cpu

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
        class(base_cpu_object), intent(inout)        :: self               !< The base backend.
        class(field_object), target :: field

        self%field => field
        self%nx    = self%field%nx
        self%ny    = self%field%ny
        self%nz    = self%field%nz
        self%ng    = self%field%grid%ng
        self%nv    = self%field%nv

        call get_mpi_basic_info(self%nprocs, self%myrank, self%masterproc, self%mpi_err)

        call self%field%check_cpu_mem(description="--Base-initialization--")


        call self%alloc()

    endsubroutine initialize


endmodule streams_base_cpu_object

