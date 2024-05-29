!< STREAmS, field class definition.
module streams_field_object

  use streams_grid_object, only : grid_object
  use streams_parameters
  use MPI
  use iso_c_binding
  use, intrinsic :: iso_fortran_env, only : error_unit
  use cpumem_mod

  implicit none
  private
  public :: field_object

  type :: field_object
    !< Field class definition.
    type(grid_object), pointer :: grid=>null()
    integer(ikind) :: nv=1_ikind
    integer(ikind) :: nx, ny, nz
    real(rkind), allocatable, dimension(:) :: x, y, z, yn
    real(rkind), allocatable, dimension(:) :: dcsidx,dcsidx2,dcsidxs
    real(rkind), allocatable, dimension(:) :: detady,detady2,detadys
    real(rkind), allocatable, dimension(:) :: dzitdz,dzitdz2,dzitdzs
    !Field variable
    real(rkind), allocatable, dimension(:,:,:,:) :: w

    !MPI data
    integer(ikind) :: mpi_err
    integer(ikind) :: myrank
    integer(ikind) :: nprocs
    logical :: masterproc
    integer, dimension(3) :: nblocks, ncoords
    integer :: mp_cart,mp_cartx,mp_carty,mp_cartz
    integer :: nproc,nrank_x, nrank_y, nrank_z
    integer :: ileftx,irightx,ilefty,irighty,ileftz,irightz
    integer :: ileftbottom, irightbottom, ilefttop, irighttop
    !ranks for edges
    integer :: lxly, lxry, rxly, rxry
    integer :: lxlz, lxrz, rxlz, rxrz
    integer :: lylz, lyrz, rylz, ryrz
    !ranks for corners
    integer :: lxlylz, lxlyrz, lxrylz, lxryrz
    integer :: rxlylz, rxlyrz, rxrylz, rxryrz


  contains
    !public methods
    procedure, pass(self) :: initialize
    procedure, pass(self) :: alloc
    procedure, pass(self) :: cartesian_mpi
    procedure, pass(self) :: correct_bc
    procedure, pass(self) :: define_corners
    procedure, pass(self) :: define_edges_corners
    procedure, pass(self) :: set_local_grid
    procedure, pass(self) :: read_field
    procedure, pass(self) :: read_field_serial
    procedure, pass(self) :: write_field
    procedure, pass(self) :: write_field_serial
    procedure, pass(self) :: write_plot3d
    procedure, pass(self) :: write_vtk
    procedure, pass(self) :: check_cpu_mem
  endtype field_object

contains

  subroutine check_cpu_mem(self, description)
    class(field_object), intent(inout) :: self !< The field backend.
    character(*) :: description
    integer :: ierr
    character(128) :: proc_name
    integer :: resultlen

    integer(c_long) :: mem_total, mem_av
    real(rkind) :: memtotal, memav

    !call mpi_barrier(mpi_comm_world, ierr)
    call mpi_get_processor_name(proc_name, resultlen, ierr)

    call getmemory(mem_total,mem_av)
    memtotal = real(mem_total,rkind)/(1024_rkind**3)
    memav = real(mem_av,rkind)/(1024_rkind**3)
    !call mpi_barrier(mpi_comm_world, ierr)
    write(error_unit, "(A,2x,A,2x,A,2x,I0,2x,I0,2x,I0)") 'CPU rank,mems: ', description,&
    &proc_name(1:resultlen),self%myrank, mem_av, mem_total
  endsubroutine check_cpu_mem

  subroutine initialize(self, grid, nv, mpi_splits)
    class(field_object), intent(inout) :: self !< The field.
    type(grid_object), intent(in), target :: grid !< Grid data.
    integer(ikind), intent(in), optional :: nv !< Number of field variables.
    integer, dimension(3) :: mpi_splits

    self%grid => grid
    self%nv = nv

    call get_mpi_basic_info(self%nprocs, self%myrank, self%masterproc, self%mpi_err)

    call self%cartesian_mpi(mpi_splits)

    call self%alloc()

    call self%set_local_grid()

  endsubroutine initialize

  subroutine cartesian_mpi(self, mpi_splits)
    class(field_object), intent(inout) :: self !< The field.
    integer, dimension(3) :: mpi_splits
    integer, parameter :: ndims=3
    logical :: remain_dims(ndims)
    logical :: reord

    associate(nblocks => self%nblocks, ncoords => self%ncoords, nprocs => self%nprocs,&
    & mpi_err => self%mpi_err, pbc => self%grid%is_xyz_periodic, masterproc => self%masterproc,&
    & mp_cart=>self%mp_cart,mp_cartx=>self%mp_cartx,mp_carty=>self%mp_carty,mp_cartz=> self%mp_cartz,&
    & nrank=>self%myrank,nproc=>self%nprocs,nrank_x=>self%nrank_x, nrank_y=>self%nrank_y,&
    & nrank_z=>self%nrank_z, ileftx=>self%ileftx,irightx=>self%irightx,ilefty=>self%ilefty,&
    &irighty=>self%irighty,ileftz=>self%ileftz,irightz=>self%irightz)

      nblocks(:) = mpi_splits(:)

      self%nx = self%grid%nxmax/nblocks(1)
      self%ny = self%grid%nymax/nblocks(2)
      self%nz = self%grid%nzmax/nblocks(3)

      !Create 3D topology
      reord = .false.
      call mpi_cart_create(mpi_comm_world,ndims,nblocks,pbc,reord,mp_cart,mpi_err)
      call mpi_cart_coords(mp_cart,nrank,ndims,ncoords,mpi_err)

      !Create 1D communicators
      remain_dims(1) = .true.
      remain_dims(2) = .false.
      remain_dims(3) = .false.
      call mpi_cart_sub(mp_cart,remain_dims,mp_cartx,mpi_err)
      call mpi_comm_rank(mp_cartx,nrank_x,mpi_err)
      call mpi_cart_shift(mp_cartx,0,1,ileftx,irightx,mpi_err)
      remain_dims(2) = .true.
      remain_dims(1) = .false.
      remain_dims(3) = .false.
      call mpi_cart_sub(mp_cart,remain_dims,mp_carty,mpi_err)
      call mpi_comm_rank(mp_carty,nrank_y,mpi_err)
      call mpi_cart_shift(mp_carty,0,1,ilefty,irighty,mpi_err)
      remain_dims(3) = .true.
      remain_dims(1) = .false.
      remain_dims(2) = .false.
      call mpi_cart_sub(mp_cart,remain_dims,mp_cartz,mpi_err)
      call mpi_comm_rank(mp_cartz,nrank_z,mpi_err)
      call mpi_cart_shift(mp_cartz,0,1,ileftz,irightz,mpi_err)

    endassociate

    call self%define_corners()

    call self%define_edges_corners()

  endsubroutine cartesian_mpi

  subroutine define_corners(self)
    class(field_object), intent(inout) :: self
    !
    integer, dimension(3) :: ncoords_nei
    !
    associate(ileftbottom => self%ileftbottom, irightbottom => self%irightbottom,&
    & ilefttop => self%ilefttop, irighttop => self%irighttop, ncoords => self%ncoords,&
    & nblocks => self%nblocks, iermpi => self%mpi_err, mp_cart => self%mp_cart, pbc => self%grid%is_xyz_p&
    &eriodic)

      ileftbottom = mpi_proc_null
      irightbottom = mpi_proc_null
      ilefttop = mpi_proc_null
      irighttop = mpi_proc_null
      !
      !ileftbottom
      ncoords_nei(1) = ncoords(1)-1
      if (pbc(1)) then
        if (ncoords_nei(1)==-1) ncoords_nei(1) = nblocks(1)-1
        if (ncoords_nei(1)==nblocks(1)) ncoords_nei(1) = 0
      endif
      ncoords_nei(2) = ncoords(2)
      ncoords_nei(3) = ncoords(3)-1
      if (pbc(3)) then
        if (ncoords_nei(3)==-1) ncoords_nei(3) = nblocks(3)-1
        if (ncoords_nei(3)==nblocks(3)) ncoords_nei(3) = 0
      endif
      if (ncoords_nei(1)>=0.and.ncoords_nei(1)<nblocks(1)) then
        if (ncoords_nei(3)>=0.and.ncoords_nei(3)<nblocks(3)) then
          call mpi_cart_rank(mp_cart,ncoords_nei,ileftbottom,iermpi)
        endif
      endif
      !ilefttop
      ncoords_nei(1) = ncoords(1)-1
      if (pbc(1)) then
        if (ncoords_nei(1)==-1) ncoords_nei(1) = nblocks(1)-1
        if (ncoords_nei(1)==nblocks(1)) ncoords_nei(1) = 0
      endif
      ncoords_nei(2) = ncoords(2)
      ncoords_nei(3) = ncoords(3)+1
      if (pbc(3)) then
        if (ncoords_nei(3)==-1) ncoords_nei(3) = nblocks(3)-1
        if (ncoords_nei(3)==nblocks(3)) ncoords_nei(3) = 0
      endif
      if (ncoords_nei(1)>=0.and.ncoords_nei(1)<nblocks(1)) then
        if (ncoords_nei(3)>=0.and.ncoords_nei(3)<nblocks(3)) then
          call mpi_cart_rank(mp_cart,ncoords_nei,ilefttop,iermpi)
        endif
      endif
      !
      !irightbottom
      ncoords_nei(1) = ncoords(1)+1
      if (pbc(1)) then
        if (ncoords_nei(1)==-1) ncoords_nei(1) = nblocks(1)-1
        if (ncoords_nei(1)==nblocks(1)) ncoords_nei(1) = 0
      endif
      ncoords_nei(2) = ncoords(2)
      ncoords_nei(3) = ncoords(3)-1
      if (pbc(3)) then
        if (ncoords_nei(3)==-1) ncoords_nei(3) = nblocks(3)-1
        if (ncoords_nei(3)==nblocks(3)) ncoords_nei(3) = 0
      endif
      if (ncoords_nei(1)>=0.and.ncoords_nei(1)<nblocks(1)) then
        if (ncoords_nei(3)>=0.and.ncoords_nei(3)<nblocks(3)) then
          call mpi_cart_rank(mp_cart,ncoords_nei,irightbottom,iermpi)
        endif
      endif
      !irighttop
      ncoords_nei(1) = ncoords(1)+1
      if (pbc(1)) then
        if (ncoords_nei(1)==-1) ncoords_nei(1) = nblocks(1)-1
        if (ncoords_nei(1)==nblocks(1)) ncoords_nei(1) = 0
      endif
      ncoords_nei(2) = ncoords(2)
      ncoords_nei(3) = ncoords(3)+1
      if (pbc(3)) then
        if (ncoords_nei(3)==-1) ncoords_nei(3) = nblocks(3)-1
        if (ncoords_nei(3)==nblocks(3)) ncoords_nei(3) = 0
      endif
      if (ncoords_nei(1)>=0.and.ncoords_nei(1)<nblocks(1)) then
        if (ncoords_nei(3)>=0.and.ncoords_nei(3)<nblocks(3)) then
          call mpi_cart_rank(mp_cart,ncoords_nei,irighttop,iermpi)
        endif
      endif
      !
    endassociate
  end subroutine define_corners

  subroutine define_edges_corners(self)
    class(field_object), intent(inout) :: self
    !
    integer, dimension(3) :: ncoords_nei
    !
    associate(myrank => self%myrank, mp_cart => self%mp_cart, lxly => self%lxly, lxry => self%lxry,&
    & rxly => self%rxly, rxry => self%rxry, lxlz => self%lxlz, lxrz => self%lxrz, rxlz => self%rxlz,&
    & rxrz => self%rxrz, lylz => self%lylz, lyrz => self%lyrz, rylz => self%rylz, ryrz => self%ryrz,&
    & lxlylz => self%lxlylz, lxlyrz => self%lxlyrz, lxrylz => self%lxrylz, lxryrz => self%lxryrz,&
    & rxlylz => self%rxlylz, rxlyrz => self%rxlyrz, rxrylz => self%rxrylz, rxryrz => self%rxryrz)

      call mpi_cart_shift_general(mp_cart, [-1,-1, 0], myrank, lxly)
      call mpi_cart_shift_general(mp_cart, [-1, 1, 0], myrank, lxry)
      call mpi_cart_shift_general(mp_cart, [ 1,-1, 0], myrank, rxly)
      call mpi_cart_shift_general(mp_cart, [ 1, 1, 0], myrank, rxry)

      call mpi_cart_shift_general(mp_cart, [-1, 0, -1], myrank, lxlz)
      call mpi_cart_shift_general(mp_cart, [-1, 0, 1], myrank, lxrz)
      call mpi_cart_shift_general(mp_cart, [ 1, 0, -1], myrank, rxlz)
      call mpi_cart_shift_general(mp_cart, [ 1, 0, 1], myrank, rxrz)

      call mpi_cart_shift_general(mp_cart, [0 ,-1, -1], myrank, lylz)
      call mpi_cart_shift_general(mp_cart, [0 ,-1, 1], myrank, lyrz)
      call mpi_cart_shift_general(mp_cart, [0 , 1, -1], myrank, rylz)
      call mpi_cart_shift_general(mp_cart, [0 , 1, 1], myrank, ryrz)

      call mpi_cart_shift_general(mp_cart, [-1 ,-1, -1], myrank, lxlylz)
      call mpi_cart_shift_general(mp_cart, [-1 ,-1, 1], myrank, lxlyrz)
      call mpi_cart_shift_general(mp_cart, [-1 , 1, -1], myrank, lxrylz)
      call mpi_cart_shift_general(mp_cart, [-1 , 1, 1], myrank, lxryrz)
      call mpi_cart_shift_general(mp_cart, [ 1 ,-1, -1], myrank, rxlylz)
      call mpi_cart_shift_general(mp_cart, [ 1 ,-1, 1], myrank, rxlyrz)
      call mpi_cart_shift_general(mp_cart, [ 1 , 1, -1], myrank, rxrylz)
      call mpi_cart_shift_general(mp_cart, [ 1 , 1, 1], myrank, rxryrz)

    endassociate
  end subroutine define_edges_corners

  subroutine alloc(self)
    class(field_object), intent(inout) :: self !< The field.
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%grid%ng, nv => self%nv)
      allocate(self%w(nv,1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
      allocate(self%x(1-ng:nx+ng))
      allocate(self%y(1-ng:ny+ng))
      allocate(self%yn(1:ny+1))
      allocate(self%z(1-ng:nz+ng))
      allocate(self%dcsidx(nx), self%dcsidx2(nx), self%dcsidxs(nx))
      allocate(self%detady(ny), self%detady2(ny), self%detadys(ny))
      allocate(self%dzitdz(nz), self%dzitdz2(nz), self%dzitdzs(nz))
    endassociate
  endsubroutine alloc

  subroutine set_local_grid(self)
    class(field_object), intent(inout) :: self !< The field.
    integer :: i, j, k
    integer :: ii, jj, kk

    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%grid%ng, ncoords => self%ncoor&
    &ds, x => self%x, xg => self%grid%xg, y => self%y, yg => self%grid%yg, z => self%z,&
    & zg => self%grid%zg, dcsidx => self%dcsidx, dcsidxs => self%dcsidxs, dcsidx2 => self%dcsidx2,&
    & detady => self%detady, detadys => self%detadys, detady2 => self%detady2, dzitdz => self%dzitdz,&
    & dzitdzs => self%dzitdzs, dzitdz2 => self%dzitdz2, dxg => self%grid%dxg, dyg => self%grid%dyg,&
    & dzg => self%grid%dzg, d2xg => self%grid%d2xg, d2yg => self%grid%d2yg, d2zg => self%grid%d2zg)
      !local coordinates (nodes)
      ii = nx*ncoords(1)
      do i=1-ng,nx+ng
        x(i) = xg(ii+i)
      enddo
      jj = ny*ncoords(2)
      do j=1-ng,ny+ng
        y(j) = yg(jj+j)
      enddo
      kk = nz*ncoords(3)
      do k=1-ng,nz+ng
        z(k) = zg(kk+k)
      enddo

      self%yn = self%grid%yn

      ii = nx*ncoords(1)
      do i=1,nx
        dcsidx (i) = 1._rkind/(dxg(ii+i))
        dcsidxs(i) = dcsidx(i)*dcsidx(i)
        dcsidx2(i) = -d2xg(ii+i)*dcsidxs(i)
      enddo
      jj = ny*ncoords(2)
      do j=1,ny
        detady (j) = 1._rkind/dyg(jj+j)
        detadys(j) = detady(j)*detady(j)
        detady2(j) = -d2yg(jj+j)*detadys(j)
      enddo
      kk = nz*ncoords(3)
      do k=1,nz
        dzitdz (k) = 1._rkind/dzg(kk+k)
        dzitdzs(k) = dzitdz(k)*dzitdz(k)
        dzitdz2(k) = -d2zg(kk+k)*dzitdzs(k)
      enddo
    endassociate
  endsubroutine set_local_grid

  subroutine read_field_serial(self)
    class(field_object), intent(inout) :: self !< The field.
    character(4) :: chx,chy,chz

    if (self%masterproc) write(*,*) 'Reading rst0_XXX_XXX_XXX.bin'
    1004 format(I4.4)
    write(chx,1004) self%ncoords(1)
    write(chy,1004) self%ncoords(2)
    write(chz,1004) self%ncoords(3)

    open (11,file='rst0_'//chx//'_'//chy//'_'//chz//'.bin',form='unformatted')
    read(11) self%w(1:self%nv,1:self%nx,1:self%ny,1:self%nz)
    close(11)
  endsubroutine read_field_serial

  subroutine read_field(self)
    class(field_object), intent(inout) :: self !< The field.
    !Writing rst.bin
    integer :: mpi_io_file
    integer :: filetype
    integer, dimension(3) :: sizes ! Dimensions of the total grid
    integer, dimension(3) :: subsizes ! Dimensions of grid local to a procs
    integer, dimension(3) :: starts ! Starting coordinates
    integer :: size_real
    integer :: ntot
    integer (kind=mpi_offset_kind) :: offset
    integer :: l,m
    integer, dimension(mpi_status_size) :: istatus
    character(len=256) :: oldname, newname
    !
    if (self%masterproc) write(*,*) 'Reading rst.bin'
    !
    sizes(1) = self%nblocks(1)*self%nx
    sizes(2) = self%nblocks(2)*self%ny
    sizes(3) = self%nblocks(3)*self%nz
    subsizes(1) = self%nx
    subsizes(2) = self%ny
    subsizes(3) = self%nz
    starts(1) = 0 + self%ncoords(1)*subsizes(1)
    starts(2) = 0 + self%ncoords(2)*subsizes(2)
    starts(3) = 0 + self%ncoords(3)*subsizes(3)
    ntot = self%nx*self%ny*self%nz
    !
    call mpi_type_create_subarray(3,sizes,subsizes,starts,mpi_order_fortran,mpi_prec,filetype,self%mpi_err)
    call mpi_type_commit(filetype,self%mpi_err)
    call mpi_file_open(self%mp_cart,'rst.bin',mpi_mode_rdonly,mpi_info_null,mpi_io_file,self%mpi_err)
    offset = 0
    do l=1,self%nv
      call mpi_file_set_view(mpi_io_file,offset,mpi_prec,filetype,"native",mpi_info_null,self%mpi_err)
      call mpi_file_read_all(mpi_io_file,self%w(l,1:self%nx,1:self%ny,1:self%nz),ntot,mpi_prec,istatus,self%mpi_err)
      call mpi_type_size(mpi_prec,size_real,self%mpi_err)
      do m=1,self%nblocks(1)*self%nblocks(2)*self%nblocks(3)
        offset = offset+size_real*ntot
      enddo
    enddo
    !
    call mpi_file_close(mpi_io_file,self%mpi_err)
    call mpi_type_free(filetype,self%mpi_err)
    !
    call mpi_barrier(mpi_comm_world, self%mpi_err)
    if (self%masterproc) then
      oldname = c_char_"rst.bin"//c_null_char
      newname = c_char_"rst.bak"//c_null_char
      self%mpi_err = rename_wrapper(oldname, newname)
      if (self%mpi_err /= 0) write(error_unit,*) "Warning! Cannot rename file rst.bin to rst.bak"
    endif
    call mpi_barrier(mpi_comm_world, self%mpi_err)
  endsubroutine read_field

  subroutine write_field_serial(self)
    class(field_object), intent(inout) :: self !< The field.
    character(4) :: chx,chy,chz

    if (self%masterproc) write(*,*) 'Writing rst1_XXX_XXX_XXX.bin'
    1004 format(I4.4)
    write(chx,1004) self%ncoords(1)
    write(chy,1004) self%ncoords(2)
    write(chz,1004) self%ncoords(3)

    open (11,file='rst1_'//chx//'_'//chy//'_'//chz//'.bin',form='unformatted')
    write(11) self%w(1:self%nv,1:self%nx,1:self%ny,1:self%nz)
    close(11)

  endsubroutine write_field_serial

  subroutine write_field(self)
    class(field_object), intent(inout) :: self !< The field.
    !Writing rst.bin and finaltime.dat
    integer :: mpi_io_file
    integer :: filetype
    integer, dimension(3) :: sizes ! Dimensions of the total grid
    integer, dimension(3) :: subsizes ! Dimensions of grid local to a procs
    integer, dimension(3) :: starts ! Starting coordinates
    integer :: size_real
    integer :: ntot
    integer (kind=mpi_offset_kind) :: offset
    integer :: l,m
    integer, dimension(mpi_status_size) :: istatus
    !
    if (self%masterproc) write(*,*) 'Writing rst.bin'
    !
    sizes(1) = self%nblocks(1)*self%nx
    sizes(2) = self%nblocks(2)*self%ny
    sizes(3) = self%nblocks(3)*self%nz
    subsizes(1) = self%nx
    subsizes(2) = self%ny
    subsizes(3) = self%nz
    starts(1) = 0 + self%ncoords(1)*subsizes(1)
    starts(2) = 0 + self%ncoords(2)*subsizes(2)
    starts(3) = 0 + self%ncoords(3)*subsizes(3)
    ntot = self%nx*self%ny*self%nz
    !
    call mpi_type_create_subarray(3,sizes,subsizes,starts,mpi_order_fortran,mpi_prec,filetype,self%mpi_err)
    call mpi_type_commit(filetype,self%mpi_err)
    call mpi_file_open(self%mp_cart,'rst.bin',mpi_mode_create+mpi_mode_wronly,mpi_info_null,mpi_io_file,self%mpi_err)
    offset = 0
    do l=1,self%nv
      call mpi_file_set_view(mpi_io_file,offset,mpi_prec,filetype,"native",mpi_info_null,self%mpi_err)
      call mpi_file_write_all(mpi_io_file,self%w(l,1:self%nx,1:self%ny,1:self%nz),ntot,mpi_prec,istatus,self%mpi_err)
      call mpi_type_size(mpi_prec,size_real,self%mpi_err)
      do m=1,self%nblocks(1)*self%nblocks(2)*self%nblocks(3)
        offset = offset+size_real*ntot
      enddo
    enddo
    !
    call mpi_file_close(mpi_io_file,self%mpi_err)
    call mpi_type_free(filetype,self%mpi_err)
    !
  endsubroutine write_field

  subroutine write_plot3d(self, mach, reynolds, time, istore, w_io, plot3dgrid, plot3dfield, file_prefix, w_aux_io)
    class(field_object), intent(inout) :: self !< The field.
    !Writing field (MPI I/O)
    real(rkind), optional :: mach, reynolds, time
    real(rkind) :: mach_, reynolds_, time_
    integer :: istore
    character(*), optional :: file_prefix
    character(32) :: file_prefix_
    logical, optional :: plot3dgrid
    logical, optional :: plot3dfield
    logical :: plot3dgrid_, plot3dfield_
    real(rkind), dimension(:,:,:,:), optional :: w_io
    real(rkind), dimension(:,:,:,:), allocatable :: w_io_
    real(rkind), dimension(:,:,:,:) :: w_aux_io
    integer :: nv_io
    integer :: i,j,k,l,m
    integer :: mpi_io_file
    integer :: filetype
    integer,dimension(3) :: sizes ! Dimensions of the total grid
    integer,dimension(3) :: subsizes ! Dimensions of grid local to a procs
    integer,dimension(3) :: starts ! Starting coordinates
    integer :: size_real,size_integer
    integer :: ntot
    integer (kind=mpi_offset_kind) :: offset
    character(4) :: chstore
    real(rkind), dimension(:,:,:,:), allocatable :: grid3d
    integer, dimension(mpi_status_size) :: istatus

    write(chstore(1:4), '(I4.4)') istore

    mach_ = 1._rkind ; if(present(mach)) mach_ = mach
    reynolds_ = 1._rkind ; if(present(reynolds)) reynolds_ = reynolds
    time_ = 1._rkind ; if(present(time)) time_ = time
    plot3dfield_ = .true. ; if(present(plot3dfield)) plot3dfield_ = plot3dfield
    plot3dgrid_ = .false. ; if(present(plot3dgrid)) plot3dgrid_ = plot3dgrid

    if(plot3dgrid_) then
      allocate(grid3d(3,self%nx,self%ny,self%nz))
      if (self%masterproc) write(*,*) 'Storing grid in plot3D format'
      !
      if (self%masterproc) then
        open(unit=123, file='plot3dgrid.xyz', access="stream", form="unformatted")
        write(123) self%grid%nxmax,self%grid%nymax,self%grid%nzmax
        flush(123)
        close(123)
      endif
      call mpi_barrier(mpi_comm_world,self%mpi_err)
      !
      do k=1,self%nz
        do j=1,self%ny
          do i=1,self%nx
            grid3d(1,i,j,k) = self%x(i)
            grid3d(2,i,j,k) = self%y(j)
            grid3d(3,i,j,k) = self%z(k)
          enddo
        enddo
      enddo
      !
      sizes(1) = self%nblocks(1)*self%nx
      sizes(2) = self%nblocks(2)*self%ny
      sizes(3) = self%nblocks(3)*self%nz
      subsizes(1) = self%nx
      subsizes(2) = self%ny
      subsizes(3) = self%nz
      starts(1) = 0 + self%ncoords(1)*subsizes(1)
      starts(2) = 0 + self%ncoords(2)*subsizes(2)
      starts(3) = 0 + self%ncoords(3)*subsizes(3)
      ntot = self%nx*self%ny*self%nz
      !
      call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,mpi_prec,filetype,self%mpi_err)
      call MPI_TYPE_COMMIT(filetype,self%mpi_err)
      call MPI_FILE_OPEN(self%mp_cart,'plot3dgrid.xyz',MPI_MODE_RDWR,MPI_INFO_NULL,mpi_io_file,self%mpi_err)
      size_integer = storage_size(1_ikind)/8
      offset = 3*size_integer
      do l=1,3
        call MPI_FILE_SET_VIEW(mpi_io_file,offset,mpi_prec,filetype,"native",MPI_INFO_NULL,self%mpi_err)
        call MPI_FILE_WRITE_ALL(mpi_io_file,grid3d(l,:,:,:),ntot,mpi_prec,istatus,self%mpi_err)
        call MPI_TYPE_SIZE(mpi_prec,size_real,self%mpi_err)
        do m=1,self%nblocks(1)*self%nblocks(2)*self%nblocks(3)
          offset = offset+size_real*ntot
        enddo
      enddo
      !
      call MPI_FILE_CLOSE(mpi_io_file,self%mpi_err)
      call MPI_TYPE_FREE(filetype,self%mpi_err)

    endif

    if (plot3dfield_) then
      if (present(w_io)) then
        nv_io = 1
        w_io_ = w_io
        file_prefix_ = file_prefix
      else
        nv_io = self%nv
        !w_io_ = self%w(:,1:self%nx,1:self%ny,1:self%nz)
        allocate(w_io_(1:self%nv,1:self%nx,1:self%ny,1:self%nz))
        w_io_(1,:,:,:) = w_aux_io(:,:,:,1)
        w_io_(2,:,:,:) = w_aux_io(:,:,:,2)
        w_io_(3,:,:,:) = w_aux_io(:,:,:,3)
        w_io_(4,:,:,:) = w_aux_io(:,:,:,4)
        w_io_(5,:,:,:) = w_aux_io(:,:,:,6)
        file_prefix_ = "flow_"
      endif

      size_real = storage_size(1.0_rkind)/8
      size_integer = storage_size(1_ikind)/8
      !
      if (self%masterproc) then
        open(unit=123, file=trim(file_prefix_)//chstore//'.q', access="stream", form="unformatted")
        if(present(w_io)) then
          write(123) self%grid%nxmax, self%grid%nymax, self%grid%nzmax, nv_io
        else
          write(123) self%grid%nxmax, self%grid%nymax, self%grid%nzmax
          write(123) mach_, 0._rkind, reynolds_, time_
        endif
        flush(123)
        close(123)
      endif
      call mpi_barrier(self%mp_cart,self%mpi_err)
      !
      sizes(1) = self%nblocks(1)*self%nx
      sizes(2) = self%nblocks(2)*self%ny
      sizes(3) = self%nblocks(3)*self%nz
      subsizes(1) = self%nx
      subsizes(2) = self%ny
      subsizes(3) = self%nz
      starts(1) = 0 + self%ncoords(1)*subsizes(1)
      starts(2) = 0 + self%ncoords(2)*subsizes(2)
      starts(3) = 0 + self%ncoords(3)*subsizes(3)
      ntot = self%nx*self%ny*self%nz
      !
      call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,mpi_prec,filetype,self%mpi_err)
      call MPI_TYPE_COMMIT(filetype,self%mpi_err)
      call MPI_FILE_OPEN(self%mp_cart,trim(file_prefix_)//chstore//'.q',MPI_MODE_RDWR,MPI_INFO_NULL,mpi_io_file,self%mpi_err)
      if(present(w_io)) then
        offset = 4*size_integer
      else
        offset = 3*size_integer+4*size_real
      endif
      do l=1,nv_io
        call MPI_FILE_SET_VIEW(mpi_io_file,offset,mpi_prec,filetype,"native",MPI_INFO_NULL,self%mpi_err)
        call MPI_FILE_WRITE_ALL(mpi_io_file,w_io_(l,:,:,:),ntot,mpi_prec,istatus,self%mpi_err)
        do m=1,self%nblocks(1)*self%nblocks(2)*self%nblocks(3)
          offset = offset+size_real*ntot
        enddo
      enddo
      !
      call MPI_FILE_CLOSE(mpi_io_file,self%mpi_err)
      call MPI_TYPE_FREE(filetype,self%mpi_err)
    endif
  end subroutine write_plot3d

  subroutine correct_bc(self, ibc)
    class(field_object), intent(inout) :: self !< The field.
    integer, dimension(6), intent(inout) :: ibc
    integer, dimension(6) :: itag
    integer :: iend, jend, kend
    itag = ibc
    ibc = 0
    iend = self%nblocks(1)-1
    jend = self%nblocks(2)-1
    kend = self%nblocks(3)-1

    if (self%ncoords(1)== 0) ibc(1) = itag(1)
    if (self%ncoords(1)==iend) ibc(2) = itag(2)
    if (self%ncoords(2)== 0) ibc(3) = itag(3)
    if (self%ncoords(2)==jend) ibc(4) = itag(4)
    if (self%ncoords(3)== 0) ibc(5) = itag(5)
    if (self%ncoords(3)==kend) ibc(6) = itag(6)
  endsubroutine correct_bc

  subroutine write_vtk(self, time, istore, w_aux_io)
    class(field_object), intent(inout) :: self !< The field.
    !Writing field (MPI I/O)
    real(rkind), optional :: time
    real(rkind) :: time_
    integer :: istore
    character(32) :: file_prefix_
    real(rkind), dimension(:,:,:,:), allocatable :: w_io_
    real(rkind), dimension(:,:,:,:) :: w_aux_io
    integer :: nv_io
    integer :: l
    integer :: mpi_io_file
    integer :: filetype
    integer,dimension(3) :: sizes ! Dimensions of the total grid
    integer,dimension(3) :: subsizes ! Dimensions of grid local to a procs
    integer,dimension(3) :: starts ! Starting coordinates
    integer :: size_real
    integer :: ntot
    integer (kind=mpi_offset_kind) :: offset
    integer (kind=mpi_offset_kind) :: offset_x,offset_y,offset_z,delta_offset_w
    character(4) :: chstore
    integer, dimension(mpi_status_size) :: istatus
    integer, parameter :: int64_kind = selected_int_kind(2*range(1))
    integer(int64_kind) :: gridsize_64
    !character(len=16) :: int2str
    !character(len=32) :: int2str_o
    character(len=65536) :: xml_part
    character(len=7) :: vtk_float
    character(len=12), dimension(5) :: names


    write(chstore(1:4), '(I4.4)') istore

    time_ = 1._rkind ; if(present(time)) time_ = time
    if (self%masterproc) print *, 'Storing VTK sol', istore,'at time', time

    call MPI_TYPE_SIZE(mpi_prec,size_real,self%mpi_err)
    if(size_real == 4) then
      vtk_float = "Float32"
    elseif(size_real == 8) then
      vtk_float = "Float64"
    else
      if(self%masterproc) write(*,*) "Error on VTK write! size_real must be either 4 or 8"
      call MPI_ABORT(MPI_COMM_WORLD,self%mpi_err,self%mpi_err)
    endif
    gridsize_64 = int(size_real,int64_kind)*int(self%grid%nxmax,int64_kind)*int(self%grid%nymax,&
    &int64_kind)*int(self%grid%nzmax,int64_kind)

    if(storage_size(gridsize_64) /= 64) then
      if(self%masterproc) write(*,*) "Error on VTK write! Size of int64_kind integers is not 8 bytes!"
      call MPI_ABORT(MPI_COMM_WORLD, self%mpi_err, self%mpi_err)
    endif

    nv_io = self%nv
    file_prefix_ = "flow_"
    !w_io_ = self%w(:,1:self%nx,1:self%ny,1:self%nz)
    allocate(w_io_(1:self%nv,1:self%nx,1:self%ny,1:self%nz))
    w_io_(1,:,:,:) = w_aux_io(:,:,:,1)
    w_io_(2,:,:,:) = w_aux_io(:,:,:,2)
    w_io_(3,:,:,:) = w_aux_io(:,:,:,3)
    w_io_(4,:,:,:) = w_aux_io(:,:,:,4)
    w_io_(5,:,:,:) = w_aux_io(:,:,:,6)

    names= [character(len=12) :: "density", "velocity_x", "velocity_y", "velocity_z", "temperature"]


    if (self%masterproc) then
      offset_x = 0
      offset_y = size_real*self%grid%nxmax + storage_size(gridsize_64)/8
      offset_z = offset_y + size_real*self%grid%nymax + storage_size(gridsize_64)/8
      delta_offset_w = gridsize_64 + storage_size(gridsize_64)/8 ! the second part is because of the header of bytes before data

      open(unit=123, file=trim(file_prefix_)//chstore//'.vtr', access="stream", form="unformatted", status="replace")
      xml_part = ' <?xml version="1.0"?> <VTKFile type="RectilinearGrid" version="1.0" byte_order="L&
      &ittleEndian" header_type="UInt64"> <RectilinearGrid WholeExtent="+1 +' //int2str(self%grid%nxmax)//&
      &' +1 +'//int2str(self%grid%nymax)//' +1 +'//int2str(self%grid%nzmax)//'"> <Piece Extent="+1 +'//&
      &int2str(self%grid%nxmax)// ' +1 +'//int2str(self%grid%nymax)//' +1 +'//int2str(self%grid%nzmax)//&
      &'"> <Coordinates> <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="X" format="appended&
      &" offset="'// int2str_o(offset_x)//'"/> <DataArray type="'//vtk_float//'" NumberOfComponents="1" Nam&
      &e="Y" format="appended" offset="'// int2str_o(offset_y)//'"/> <DataArray type="'//vtk_float//&
      &'" NumberOfComponents="1" Name="Z" format="appended" offset="'// int2str_o(offset_z)//'"/> </&
      &Coordinates> <PointData> '
      offset = offset_z + size_real*self%grid%nzmax + storage_size(gridsize_64)/8
      do l=1,nv_io
        xml_part = trim(adjustl(xml_part)) // ' <DataArray type="'//vtk_float//'" NumberOfComponents&
        &="1" Name="'//trim(names(l))//'" format="appended" offset="'//int2str_o(offset)//'"/>'
        offset = offset + delta_offset_w
      enddo
      xml_part = trim(adjustl(xml_part)) // ' </PointData> </Piece> </RectilinearGrid> <AppendedData encoding="raw"> '
      write(123) trim(adjustl(xml_part))

      write(123) "_"
      !write(123) storage_size(gridsize_64)/8*int(nxmax,int64_kind) , xg(1:nxmax)
      !write(123) storage_size(gridsize_64)/8*int(nymax,int64_kind) , yg(1:nymax)
      !write(123) storage_size(gridsize_64)/8*int(nzmax,int64_kind) , zg(1:nzmax)
      write(123) size_real*int(self%grid%nxmax,int64_kind) , self%grid%xg(1:self%grid%nxmax)
      write(123) size_real*int(self%grid%nymax,int64_kind) , self%grid%yg(1:self%grid%nymax)
      write(123) size_real*int(self%grid%nzmax,int64_kind) , self%grid%zg(1:self%grid%nzmax)
      flush(123)
      close(123)
    endif
    !
    sizes(1) = self%nblocks(1)*self%nx
    sizes(2) = self%nblocks(2)*self%ny
    sizes(3) = self%nblocks(3)*self%nz
    subsizes(1) = self%nx
    subsizes(2) = self%ny
    subsizes(3) = self%nz
    starts(1) = 0 + self%ncoords(1)*subsizes(1)
    starts(2) = 0 + self%ncoords(2)*subsizes(2)
    starts(3) = 0 + self%ncoords(3)*subsizes(3)
    ntot = self%nx*self%ny*self%nz
    !
    call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,mpi_prec,filetype,self%mpi_err)
    call MPI_TYPE_COMMIT(filetype,self%mpi_err)
    !
    do l=1,nv_io
      call MPI_BARRIER(self%mp_cart,self%mpi_err)
      if (self%masterproc) then
        open(unit=123, file=trim(file_prefix_)//chstore//'.vtr', access="stream", form="unformatted", position="append")
        write(123) gridsize_64
        flush(123)
        close(123)
      endif
      call MPI_BARRIER(self%mp_cart, self%mpi_err)
      call MPI_FILE_OPEN(self%mp_cart,trim(file_prefix_)//chstore//'.vtr',MPI_MODE_RDWR,MPI_INFO_NULL,mpi_io_file,self%mpi_err)
      call MPI_FILE_GET_SIZE(mpi_io_file, offset, self%mpi_err)
      call MPI_BARRIER(self%mp_cart, self%mpi_err)
      call MPI_FILE_SET_VIEW(mpi_io_file,offset,mpi_prec,filetype,"native",MPI_INFO_NULL,self%mpi_err)
      call MPI_FILE_WRITE_ALL(mpi_io_file,w_io_(l,1:self%nx,1:self%ny,1:self%nz),ntot,mpi_prec,istatus,self%mpi_err)
      call MPI_FILE_CLOSE(mpi_io_file,self%mpi_err)
    enddo

    call MPI_TYPE_FREE(filetype,self%mpi_err)
    if (self%masterproc) then
      open(unit=123, file=trim(file_prefix_)//chstore//'.vtr', access="stream", position="append", form="unformatted")
      write(123) ' </AppendedData> </VTKFile>'
      close(123)
    endif
  end subroutine write_vtk

endmodule streams_field_object

