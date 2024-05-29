!< STREAmS, general parameters.
module streams_parameters

  use, intrinsic :: iso_fortran_env
  use, intrinsic :: ieee_arithmetic
  use mpi
  use iso_c_binding
  use tcp

  implicit none
  private
  public :: ikind, ikind64, rkind, mpi_prec
  public :: tol_iter, tol_iter_nr
  public :: int2str_o, int2str
  public :: get_mpi_basic_info
  public :: mpi_initialize
  public :: rename_wrapper
  public :: pi
  public :: c_rkind
  public :: pol_int, locateval
  public :: invmat, detmat, fail_input_any
  public :: ieee_is_nan
  public :: error_unit
  public :: mpi_cart_shift_general
  !
  !INSITU
  public :: insitu_start, insitu_end
  public :: REAL64
  public :: INT64
  public :: byte_size

  integer, parameter :: ikind = INT32
  integer, parameter :: ikind64 = INT64
#ifdef SINGLE_PRECISION
  integer, parameter :: rkind = REAL32
  integer, parameter :: mpi_prec = mpi_real4
  real(rkind), parameter :: tol_iter    = 0.00001_rkind
  real(rkind), parameter :: tol_iter_nr = 0.00001_rkind
  integer, parameter :: c_rkind = C_FLOAT
#else
  integer, parameter :: rkind = REAL64
  integer, parameter :: mpi_prec = mpi_real8
  real(rkind), parameter :: tol_iter = 0.000000001_rkind
  real(rkind), parameter :: tol_iter_nr = 0.000000000001_rkind
  integer, parameter :: c_rkind = C_DOUBLE
#endif

  real(rkind) :: pi = acos(-1._rkind)

  interface
    function rename_wrapper(filein, fileout) bind(C, name="rename")
      import :: c_char, c_int
      integer(c_int) :: rename_wrapper
      character(kind=c_char) :: filein(*), fileout(*)
    endfunction rename_wrapper
  endinterface

  interface byte_size
    module procedure byte_size_int32, &
#ifdef SINGLE_PRECISION
    byte_size_real32
#else
    byte_size_real64
#endif
  endinterface

contains

  subroutine mpi_initialize()
    integer :: mpi_err
    call mpi_init(mpi_err)
  endsubroutine

  subroutine get_mpi_basic_info(nprocs, myrank, masterproc, mpi_err)
    integer :: nprocs, myrank, mpi_err
    logical :: masterproc
    call mpi_comm_size(mpi_comm_world, nprocs, mpi_err)
    call mpi_comm_rank(mpi_comm_world, myrank, mpi_err)
    masterproc = .false.
    if (myrank == 0) masterproc = .true.
  endsubroutine get_mpi_basic_info

  function int2str(int_num)
    implicit none
    integer :: int_num
    character(len=16) :: int2str, ret_value
    write(ret_value, "(I0)") int_num
    int2str = ret_value
  endfunction int2str

  function int2str_o(int_num)
    use mpi
    implicit none
    integer(KIND=MPI_OFFSET_KIND) :: int_num
    character(len=32) :: int2str_o, ret_value
    write(ret_value, "(I0)") int_num
    int2str_o = ret_value
  endfunction int2str_o

  subroutine pol_int(x,y,n,xs,ys)
    !
    !Polynomial interpolation using Neville's algorithm
    !Order of accuracy of the interpolation is n-1
    !
    integer, intent(in) :: n
    real(rkind), dimension(n), intent(in) :: x,y
    real(rkind), intent(in) :: xs
    real(rkind), intent(out) :: ys
    !
    integer :: i,m
    real(rkind), dimension(n)  :: v,vold
    !
    v = y
    !
    do m=2,n ! Tableu columns
      vold = v
      do i=1,n+1-m
        v(i) = (xs-x(m+i-1))*vold(i)+(x(i)-xs)*vold(i+1)
        v(i) = v(i)/(x(i)-x(i+m-1))
      enddo
    enddo
    ys = v(1)
  end subroutine pol_int
  !
  subroutine locateval(xx,n,x,ii)
    !
    integer, intent(in) :: n
    integer, intent(out) :: ii
    real(rkind), dimension(1:n), intent(in) :: xx
    real(rkind) :: x
    integer :: il,jm,juu
    !
    il=0
    juu=n+1
    do while (juu-il.gt.1)
      jm=(juu+il)/2
      if ((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm))) then
        il=jm
      else
        juu=jm
      endif
    end do
    ii=il
  end subroutine locateval
  !
  subroutine fail_input_any(msg)
    integer :: mpi_err
    character(len=*) :: msg
    write(error_unit,*) "Input Error! ", msg
    call mpi_abort(mpi_comm_world,15,mpi_err)
  endsubroutine fail_input_any
  !
  subroutine invmat(mat,n)
    !
    integer, intent(in) :: n
    real(rkind), dimension(n,n), intent(inout) :: mat
    integer :: i,j
    integer, dimension(n) :: indx
    real(rkind), dimension(n,n) :: y
    real(rkind) :: d
    !
    y = 0._rkind
    do i=1,n
      y(i,i) = 1._rkind
    enddo
    call ludcmp(mat,n,n,indx,d)
    do j=1,n
      call lubksb(mat,n,n,indx,y(1,j))
    enddo
    mat = y
    !
  end subroutine invmat
  !
  SUBROUTINE LUDCMP(A,N,NP,INDX,D)
  integer, parameter :: NMAX = 100
  real(rkind), parameter :: TINY = 1.0D-20
  integer, intent(in) :: N,NP
  integer, dimension(N), intent(out) :: INDX
  real(rkind), intent(out) :: D
  real(rkind), dimension(NP,NP), intent(inout) :: A
  real(rkind), dimension(NMAX) :: VV
  integer :: I,J,K,IMAX
  real(rkind) :: SUM,DUM,AAMAX
  !
  D=1._rkind
  DO 12 I=1,N
  AAMAX=0._rkind
  DO 11 J=1,N
  IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
  11      CONTINUE
  IF (AAMAX.EQ.0._rkind) write(*,*) 'Error, singular matrix.'
  VV(I)=1./AAMAX
  12    CONTINUE
  DO 19 J=1,N
  IF (J.GT.1) THEN
  DO 14 I=1,J-1
  SUM=A(I,J)
  IF (I.GT.1)THEN
  DO 13 K=1,I-1
  SUM=SUM-A(I,K)*A(K,J)
  13            CONTINUE
  A(I,J)=SUM
  ENDIF
  14        CONTINUE
  ENDIF
  AAMAX=0.
  DO 16 I=J,N
  SUM=A(I,J)
  IF (J.GT.1)THEN
  DO 15 K=1,J-1
  SUM=SUM-A(I,K)*A(K,J)
  15          CONTINUE
  A(I,J)=SUM
  ENDIF
  DUM=VV(I)*ABS(SUM)
  IF (DUM.GE.AAMAX) THEN
  IMAX=I
  AAMAX=DUM
  ENDIF
  16      CONTINUE
  IF (J.NE.IMAX)THEN
  DO 17 K=1,N
  DUM=A(IMAX,K)
  A(IMAX,K)=A(J,K)
  A(J,K)=DUM
  17        CONTINUE
  D=-D
  VV(IMAX)=VV(J)
  ENDIF
  INDX(J)=IMAX
  IF(J.NE.N)THEN
  IF(A(J,J).EQ.0._rkind)A(J,J)=TINY
  DUM=1./A(J,J)
  DO 18 I=J+1,N
  A(I,J)=A(I,J)*DUM
  18        CONTINUE
  ENDIF
  19    CONTINUE
  IF(A(N,N).EQ.0._rkind)A(N,N)=TINY
  RETURN
  END
  !
  SUBROUTINE LUBKSB(A,N,NP,INDX,B)
  integer, intent(in) :: N,NP
  integer, dimension(N), intent(in) :: INDX
  real(rkind), dimension(NP,NP), intent(in) :: A
  real(rkind), dimension(N), intent(inout) :: B
  integer :: II,I,LL,J
  real(rkind) :: SUM
  II=0
  DO 12 I=1,N
  LL=INDX(I)
  SUM=B(LL)
  B(LL)=B(I)
  IF (II.NE.0)THEN
  DO 11 J=II,I-1
  SUM=SUM-A(I,J)*B(J)
  11        CONTINUE
  ELSE IF (SUM.NE.0.) THEN
  II=I
  ENDIF
  B(I)=SUM
  12    CONTINUE
  DO 14 I=N,1,-1
  SUM=B(I)
  IF(I.LT.N)THEN
  DO 13 J=I+1,N
  SUM=SUM-A(I,J)*B(J)
  13        CONTINUE
  ENDIF
  B(I)=SUM/A(I,I)
  14    CONTINUE
  RETURN
  END
  !
  subroutine detmat(mat,n,det)
    !
    integer, intent(in) :: n
    real(rkind), dimension(n,n), intent(inout) :: mat
    real(rkind), intent(out) :: det
    integer :: j
    integer, dimension(n) :: indx
    real(rkind) :: d
    !
    call ludcmp(mat,n,n,indx,d)
    do j=1,n
      d = d*mat(j,j)
    enddo
    det = d
    !
  end subroutine detmat
  !
  subroutine insitu_start(fcoproc,vtkpipeline,masterproc)
    character(len=*), intent(in) :: vtkpipeline
    logical, intent(in) :: masterproc
    logical, intent(inout) :: fcoproc
    !
    inquire(file=vtkpipeline, exist=fcoproc)

    if (fcoproc) then
      if (masterproc) print*, 'Connecting to Catalyst...'
      if (masterproc) print*, 'Adding '//vtkpipeline
      CALL coprocessorinitializewithpython(vtkpipeline,len(vtkpipeline)) ! Initialize Catalyst
    else
      if (masterproc) print*, 'WARNING: ', vtkpipeline, ' is missing'
    endif
    !
  end subroutine insitu_start

  subroutine insitu_end(fcoproc)
    logical, intent(in) :: fcoproc
    !
    if (fcoproc) CALL coprocessorfinalize() ! Finalize Catalyst
    !
  end subroutine insitu_end

  elemental function byte_size_int32(x) result(bytes)
  integer, intent(in) :: x
  integer(ikind64)    :: bytes

  bytes = storage_size(x) / 8
endfunction byte_size_int32
#ifdef SINGLE_PRECISION
elemental function byte_size_real32(x) result(bytes)
real(rkind), intent(in) :: x
integer(ikind64)        :: bytes

bytes = storage_size(x) / 8
endfunction byte_size_real32
#else
elemental function byte_size_real64(x) result(bytes)
real(rkind), intent(in) :: x
integer(ikind64)        :: bytes

bytes = storage_size(x) / 8
endfunction byte_size_real64
#endif

subroutine mpi_cart_shift_general(mp_cart, shifts, rank_source, rank_dest)
integer, intent(in) :: mp_cart, rank_source
integer, intent(out) :: rank_dest
integer, dimension(3), intent(in) :: shifts
integer, dimension(3) :: ncoords_source, ncoords_dest, dims, coords
logical, dimension(3) :: pbc
integer :: ierr, i

!get communicator features
call mpi_cart_get(mp_cart, 3, dims, pbc, coords, ierr)
!print*,'Periodicity active on directions: ',pbc
!get coordinates of rank_source
call mpi_cart_coords(mp_cart, rank_source, 3, ncoords_source, ierr)
!print*,'coords/my_coords should be equal: ',coords, ncoords_source

ncoords_dest(1:3) = ncoords_source(1:3) + shifts(1:3)
do i=1,3
if (pbc(i)) then
  if (ncoords_dest(i)==-1) ncoords_dest(i) = dims(i)-1
  if (ncoords_dest(i)==dims(i)) ncoords_dest(i) = 0
endif
enddo
if(any(ncoords_dest(1:3)==-1) .or. any(ncoords_dest(1:3)==dims(1:3))) then
rank_dest = mpi_proc_null
else
call mpi_cart_rank(mp_cart,ncoords_dest,rank_dest,ierr)
endif
endsubroutine mpi_cart_shift_general
!
endmodule streams_parameters

