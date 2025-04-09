module tcp
  use iso_c_binding
  implicit none
#ifdef SINGLE_PRECISION
  integer, parameter :: c_real = c_float
#else
  integer, parameter :: c_real = c_double
#endif
  public


  interface
    subroutine createcpimagedata(nxs_ins,nxe_ins,nzs_ins,nze_ins&
                                ,nx,nz &
                                ,dx,dz &
                                ,ytop &
                                ,zone &
                                ,nrank) bind(C)
    import :: c_char, c_int, c_real
    real(c_real) :: dx, dz, ytop
    integer(c_int) :: nxs_ins, nxe_ins
    integer(c_int) :: nys_ins, nye_ins
    integer(c_int) :: nzs_ins, nze_ins
    integer(c_int) :: nx, nz
    integer(c_int) :: nrank
    character(kind=c_char) :: zone(*)
    endsubroutine createcpimagedata
  endinterface

  interface
    subroutine createcprectilineardata(nxs_ins,nxe_ins,nys_ins,nye_ins,nzs_ins,nze_ins&
                            ,nxstartg,nystartg,nzstartg &
                            ,nxendg,nyendg,nzendg &
                            ,x,y,z &
                            ,nrank,nproc &
                            ,zone) bind(C)
    import :: c_char, c_int, c_real
    real(c_real) :: x(*), y(*), z(*) 
    integer(c_int) :: nxs_ins, nxe_ins
    integer(c_int) :: nys_ins, nye_ins
    integer(c_int) :: nzs_ins, nze_ins
    integer(c_int) :: nxstartg, nystartg, nzstartg
    integer(c_int) :: nxendg, nyendg, nzendg
    integer(c_int) :: nrank, nproc
    character(kind=c_char) :: zone(*)
    endsubroutine createcprectilineardata
  endinterface

  interface
    subroutine createcpstructureddata(nxp2,nyp2,nzp2&
                            ,nxs_ins,nys_ins,nzs_ins&
                            ,nxe_ins,nye_ins,nze_ins&
                            ,nxmaxp2,nymaxp2,nzmaxp2 &
                            ,x,y,z &
                            ,nrank,nproc &
                            ,zone,xyzc) bind(C)
    import :: c_char, c_int, c_real
    real(c_real) :: x(*), y(*), z(*) 
    real(c_real) :: xyzc(*)
    integer(c_int) :: nxp2
    integer(c_int) :: nyp2
    integer(c_int) :: nzp2
    integer(c_int) :: nxs_ins, nxe_ins
    integer(c_int) :: nys_ins, nye_ins
    integer(c_int) :: nzs_ins, nze_ins
    integer(c_int) :: nxmaxp2, nymaxp2, nzmaxp2
    integer(c_int) :: nrank, nproc
    character(kind=c_char) :: zone(*)
    endsubroutine createcpstructureddata
  endinterface

  interface
    subroutine addfieldtoimage(zone,psixztop,name,nrank) bind(C)
    import :: c_char, c_int, c_real
    real(c_real) :: psixztop(*)
    integer(c_int) :: nrank
    character(kind=c_char) :: zone(*), name(*)
    endsubroutine addfieldtoimage
  endinterface

  interface
    subroutine addfieldtorectilinear(zone,psi,name,nrank) bind(C)
    import :: c_char, c_int, c_real
    real(c_real) :: psi(*)
    integer(c_int) :: nrank
    character(kind=c_char) :: zone(*), name(*)
    endsubroutine addfieldtorectilinear
  endinterface

  interface
    subroutine addfieldtostructured(zone,psi,name,nrank) bind(C)
    import :: c_char, c_int, c_real
    real(c_real) :: psi(*)
    integer(c_int) :: nrank
    character(kind=c_char) :: zone(*), name(*)
    endsubroutine addfieldtostructured
  endinterface

  interface
    subroutine coprocessorinitializewithpython(vtkpipeline,length) bind(C)
    import :: c_char, c_int, c_real
    character(kind=c_char) :: vtkpipeline(*)
    integer(c_int) :: length
    endsubroutine coprocessorinitializewithpython
  endinterface

  interface
    subroutine requestdatadescription(icyc,time,flag) bind(C)
    import :: c_char, c_int, c_real
    real(c_real) :: time
    integer(c_int) :: icyc, flag
    endsubroutine requestdatadescription
  endinterface

  interface
    subroutine needtocreategrid(flag) bind(C)
    import :: c_char, c_int, c_real
    integer(c_int) :: flag
    endsubroutine needtocreategrid
  endinterface

  interface
    subroutine coprocess() bind(C)
    import :: c_char, c_int, c_real
    endsubroutine coprocess
  endinterface

  interface
    subroutine coprocessorfinalize() bind(C)
    endsubroutine coprocessorfinalize
  endinterface

end module tcp
