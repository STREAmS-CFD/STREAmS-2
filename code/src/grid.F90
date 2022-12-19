!< STREAmS, grid class definition.
module streams_grid_object
!
  use streams_parameters
  use mpi
!
  implicit none
  private
  public :: grid_object
  public :: GRID_FROMFILE, GRID_UNIFORM, GRID_CHA, GRID_BL, GRID_SBLI
!
  integer(ikind), parameter :: GRID_FROMFILE        = 1_ikind
  integer(ikind), parameter :: GRID_UNIFORM         = 2_ikind
  integer(ikind), parameter :: GRID_CHA             = 3_ikind
  integer(ikind), parameter :: GRID_BL              = 4_ikind
  integer(ikind), parameter :: GRID_SBLI            = 5_ikind
!
  type :: grid_object
!   < Grid class definition.
    integer                    :: myrank, nprocs, mpi_err
    logical                    :: masterproc
    real(rkind), dimension(3)  :: domain_size
    integer                    :: nxmax, nymax, nzmax
    integer                    :: ng
    logical, dimension(3)      :: is_xyz_periodic = .false.
    logical                    :: is_y_staggered  = .false.
    integer(ikind)             :: grid_type
    real(rkind), allocatable, dimension(:) :: xg, yg, zg
    real(rkind), allocatable, dimension(:) :: dxg, dyg, dzg
    real(rkind), allocatable, dimension(:) :: d2xg, d2yg, d2zg
    real(rkind), allocatable, dimension(:) :: yn
    real(rkind), dimension(4)    :: metrics_cd1
    real(rkind), dimension(0:4)  :: metrics_cd2
    integer :: metrics_order
  contains
!   public methods
    procedure, pass(self) :: alloc
    procedure, pass(self) :: compute_metrics
    procedure, pass(self) :: generate_grid_uniform_x
    procedure, pass(self) :: generate_grid_uniform_y
    procedure, pass(self) :: generate_grid_uniform_z
    procedure, pass(self) :: generate_grid_cha
    procedure, pass(self) :: generate_grid_bl
    procedure, pass(self) :: generate_grid_geoprogression
    procedure, pass(self) :: generate_grid_uni_extension
    procedure, pass(self) :: get_deriv_coeffs
    procedure, pass(self) :: initialize
    procedure, pass(self) :: read_grid
    procedure, pass(self) :: write_grid
  endtype grid_object
! 
contains
! 
! public methods
  subroutine compute_metrics(self)
!   < Compute metrics of a block.
    class(grid_object), intent(inout) :: self  !< The grid.
    integer :: mm, i, j, k, l
!   
    call self%get_deriv_coeffs(self%metrics_cd1, self%metrics_cd2, self%metrics_order, self%metrics_order)
!   
    mm = self%metrics_order/2
!
    associate(xg         => self%xg,    yg    => self%yg,    zg    => self%zg,    &
      dxg        => self%dxg,   dyg   => self%dyg,   dzg   => self%dzg,   &
      d2xg       => self%d2xg,  d2yg  => self%d2yg,  d2zg  => self%d2zg,  &
      nxmax      => self%nxmax, nymax => self%nymax, nzmax => self%nzmax, &
      c          => self%metrics_cd1, &
      cc         => self%metrics_cd2, &
      masterproc => self%masterproc)
!     
      dxg = 0._rkind
      do i=1,nxmax
        do l=1,mm
          dxg(i) = dxg(i)+c(l)*(xg(i+l)-xg(i-l))
        enddo
        d2xg(i) = cc(0)*xg(i)
        do l=1,mm
          d2xg(i) = d2xg(i)+cc(l)*(xg(i+l)+xg(i-l))
        enddo
      enddo
!     
      dyg = 0._rkind
      do j=1,nymax
        do l=1,mm
          dyg(j) = dyg(j)+c(l)*(yg(j+l)-yg(j-l))
        enddo
        d2yg(j) = cc(0)*yg(j)
        do l=1,mm
          d2yg(j) = d2yg(j)+cc(l)*(yg(j+l)+yg(j-l))
        enddo
      enddo
!     
      dzg = 0._rkind
      do k=1,nzmax
        do l=1,mm
          dzg(k) = dzg(k)+c(l)*(zg(k+l)-zg(k-l))
        enddo
        d2zg(k) = cc(0)*zg(k)
        do l=1,mm
          d2zg(k) = d2zg(k)+cc(l)*(zg(k+l)+zg(k-l))
        enddo
      enddo
!     
      if(masterproc) then
        open(18,file='dxg.dat')
        do i=1,nxmax
          write(18,*) xg(i),dxg(i),d2xg(i)
        enddo
        close(18)
        open(18,file='dyg.dat')
        do j=1,nymax
          write(18,*) yg(j),dyg(j),d2yg(j)
        enddo
        close(18)
        open(18,file='dzg.dat')
        do k=1,nzmax
          write(18,*) zg(k),dzg(k),d2zg(k)
        enddo
        close(18)
      endif
!     
    endassociate
!   
  endsubroutine compute_metrics
! 
  subroutine get_deriv_coeffs(self, coeff_deriv1, coeff_deriv2, order1, order2)
    class(grid_object)  :: self
    real(rkind), dimension(4)   :: coeff_deriv1
    real(rkind), dimension(0:4) :: coeff_deriv2
    integer :: order1, order2
!   
    coeff_deriv1 = 0._rkind
    coeff_deriv2 = 0._rkind
!   
    select case (order1/2)
    case (1)
      coeff_deriv1(1) = 0.5_rkind
    case (2)
      coeff_deriv1(1) =  2._rkind/3._rkind
      coeff_deriv1(2) = -1._rkind/12._rkind
    case (3)
      coeff_deriv1(1) =  0.75_rkind
      coeff_deriv1(2) = -0.15_rkind
      coeff_deriv1(3) =    1._rkind/60._rkind
    case(4)
      coeff_deriv1(1) =  4._rkind/5._rkind
      coeff_deriv1(2) = -1._rkind/5._rkind
      coeff_deriv1(3) =  4._rkind/105._rkind
      coeff_deriv1(4) = -1._rkind/280._rkind
    end select
!   
    select case (order2/2)
    case (1)
      coeff_deriv2(0) = -2._rkind
      coeff_deriv2(1) =  1._rkind
    case (2)
      coeff_deriv2(0) = -2.5_rkind
      coeff_deriv2(1) =   4._rkind/3._rkind
      coeff_deriv2(2) =  -1._rkind/12._rkind
    case (3)
      coeff_deriv2(0) = -245._rkind/90._rkind
      coeff_deriv2(1) =   1.5_rkind
      coeff_deriv2(2) = -0.15_rkind
      coeff_deriv2(3) =    1._rkind/90._rkind
    case (4)
      coeff_deriv2(0) = -205._rkind/72._rkind
      coeff_deriv2(1) =    8._rkind/5._rkind
      coeff_deriv2(2) =   -1._rkind/5._rkind
      coeff_deriv2(3) =    8._rkind/315._rkind
      coeff_deriv2(4) =   -1._rkind/560._rkind
    endselect
!   
  endsubroutine get_deriv_coeffs
! 
  subroutine initialize(self, periodic, nxmax, nymax, nzmax, ng, grid_type, &
    domain_size_x, domain_size_y, domain_size_z, l0, &
    grid_vars, metrics_order, rebuild_ghost, ystaggering)
    class(grid_object)  :: self
    integer             :: nxmax, nymax, nzmax, ng, grid_type
    real(rkind)         :: domain_size_x, domain_size_y, domain_size_z, l0
    integer             :: nxmax_tot, nymax_tot, nzmax_tot
    integer             :: metrics_order
    logical             :: rebuild_ghost,ystaggering
    logical             :: write_grid
    logical, dimension(3) :: periodic
    real(rkind), dimension(:), allocatable :: grid_vars
!   
    call get_mpi_basic_info(self%nprocs, self%myrank, self%masterproc, self%mpi_err)
!   
    self%nxmax           = nxmax
    self%nymax           = nymax
    self%nzmax           = nzmax
    self%ng              = ng
    self%grid_type       = grid_type
    self%metrics_order   = metrics_order
    self%is_xyz_periodic = periodic
!   
    self%domain_size     = [domain_size_x,domain_size_y,domain_size_z]
!   
    call self%alloc()
!   
    write_grid = .true.
!   
    self%is_y_staggered = ystaggering
!   
    if (self%grid_type == GRID_FROMFILE) then
      call self%read_grid(rebuild_ghost)
      nxmax_tot = self%nxmax ; if (self%is_xyz_periodic(1)) nxmax_tot = nxmax_tot + 1
      nymax_tot = self%nymax ; if (self%is_xyz_periodic(2)) nymax_tot = nymax_tot + 1
      nzmax_tot = self%nzmax ; if (self%is_xyz_periodic(3)) nzmax_tot = nzmax_tot + 1
      domain_size_x = self%xg(nxmax_tot)-self%xg(1)
      domain_size_y = self%yg(nymax_tot)-self%yg(1)
      domain_size_z = self%zg(nzmax_tot)-self%zg(1)
      self%domain_size = [domain_size_x,domain_size_y,domain_size_z]
      write_grid = .false.
    else
      call self%generate_grid_uniform_x()
      call self%generate_grid_uniform_z()
      select case (self%grid_type)
      case(GRID_UNIFORM)
        call self%generate_grid_uniform_y()
      case(GRID_CHA)
        call self%generate_grid_cha(grid_vars(1), grid_vars(2),nint(grid_vars(3)), grid_vars(4))
      case(GRID_BL)
        call self%generate_grid_bl(grid_vars(1), grid_vars(2), nint(grid_vars(3)), grid_vars(4), grid_vars(5))
        if (nint(grid_vars(3))<self%nymax) then
          call self%generate_grid_geoprogression(nint(grid_vars(3)), grid_vars(4), nint(grid_vars(6)))
        else
          self%domain_size(2) = self%yg(nymax_tot)-self%yg(1)
        endif
!     case(GRID_SBLI)
!       call self%generate_grid_bl(grid_vars(1), grid_vars(2), nint(grid_vars(3)), grid_vars(4), grid_vars(5))
!       !call self%generate_grid_uni_extension(nint(grid_vars(3)), grid_vars(4), nint(grid_vars(6)))
!       self%domain_size(2) = self%yg(nymax_tot)-self%yg(1)
      end select
    endif
!   
    if (self%masterproc.and.write_grid) call self%write_grid()
!   
    self%xg = self%xg*l0
    self%yg = self%yg*l0
    self%zg = self%zg*l0
    self%yn = self%yn*l0
    self%domain_size = self%domain_size*l0
!   
    call self%compute_metrics()
!   
  endsubroutine initialize
! 
  subroutine generate_grid_bl(self, jbgrid, dyptarget, nymaxwr, rlywr, retau)
    class(grid_object), intent(inout) :: self
    real(rkind),        intent(in) :: jbgrid,dyptarget,rlywr,retau
    integer(ikind),     intent(in) :: nymaxwr
    real(rkind)                    :: ynf,rj,alfc
!   real(rkind)                    :: alfb,ceta
    integer                        :: j
!   
    associate(nymax=>self%nymax,  ng => self%ng, yg   => self%yg,    yn => self%yn, &
      domain_size => self%domain_size)
!     
!     alfb = 1.25_rkind
!     ceta = 0.8_rkind
!     alfc = alfb*ceta
      alfc = 4._rkind/3._rkind*(retau*rlywr)**0.75_rkind/nymaxwr
      do j=1,nymaxwr+ng
        rj    = real(j-1,rkind)
        yg(j) = 1._rkind/(1._rkind+(rj/jbgrid)**2)
        yg(j) = yg(j)*(rj*dyptarget+(0.75_rkind*alfc*rj)**(4._rkind/3._rkind)*(rj/jbgrid)**2)
      enddo
      ynf = yg(nymaxwr)
      do j=1,nymaxwr+ng
        yg(j) = yg(j)/ynf*rlywr
      enddo
!     
      do j=1,ng
        yg(1-j) = -yg(1+j)
      enddo
!     
    endassociate
!   
  endsubroutine generate_grid_bl
! 
  subroutine generate_grid_geoprogression(self, nymaxwr, rlywr, nsmooy)
    class(grid_object), intent(inout) :: self  !< The grid.
    real(rkind),        intent(in) :: rlywr
    integer,            intent(in) :: nymaxwr, nsmooy
    real(rkind) :: dy1, rlygp, rgp, rgpold
    integer :: nygp, j, l
    real(rkind), dimension(self%nymax+self%ng) :: tmpyg
!   
    associate(nymax=>self%nymax,  ng => self%ng,  &
      yg   => self%yg,    yn => self%yn, rly =>  self%domain_size(2), &
      domain_size => self%domain_size)
!     
!     dy1 = yg(nymaxwr+1)-yg(nymaxwr)
      dy1 = yg(nymaxwr)-yg(nymaxwr-1)
!     
      rlygp = rly-rlywr
      nygp  = nymax-nymaxwr
      rgp = 1.01_rkind
      do
        rgpold = rgp
        rgp = 1._rkind-(1._rkind-rgp)/dy1*rlygp
        rgp = rgp**(1._rkind/nygp)
        if (abs(rgp-rgpold) < tol_iter) exit
!       write(*,*) 'Geometric progression factor rgp =', rgp
      enddo
      if (self%masterproc) write(*,*) 'Geometric progression factor rgp =', rgp
!     
      do j=1,nygp+ng
        yg(nymaxwr+j) = yg(nymaxwr+j-1) + dy1*rgp**(j-1._rkind)
      enddo
!     
      do l=1,nsmooy
        tmpyg = yg(1:nymax+ng)
        do j=2,nymax+ng-1
          yg(j) = tmpyg(j-1)+4._rkind*tmpyg(j)+tmpyg(j+1)
          yg(j) = yg(j)/6._rkind
        enddo
      enddo
!     
      do j=1,ng
        yg(1-j) = -yg(1+j)
      enddo
!     
    endassociate
!   
  endsubroutine generate_grid_geoprogression
! 
  subroutine generate_grid_uni_extension(self, nymaxwr, rlywr, nsmooy)
    class(grid_object), intent(inout) :: self  !< The grid.
    real(rkind),        intent(in) :: rlywr
    integer,            intent(in) :: nymaxwr, nsmooy
    real(rkind) :: dy1, rlygp
    integer :: nygp, j, l
    real(rkind), dimension(self%nymax+self%ng) :: tmpyg
!   
    associate(nymax=>self%nymax,  ng => self%ng,  &
      yg   => self%yg,    yn => self%yn, rly =>  self%domain_size(2), &
      domain_size => self%domain_size)
!     
      dy1 = yg(nymaxwr+1)-yg(nymaxwr)
!     
      rlygp = rly-rlywr
      nygp  = nymax-nymaxwr
!     
      do j=1,nygp+ng
        yg(nymaxwr+j) = yg(nymaxwr+j-1) + dy1
      enddo
!     
      do l=1,nsmooy
        tmpyg = yg(1:nymax+ng)
        do j=2,nymax+ng-1
          yg(j) = tmpyg(j-1)+4._rkind*tmpyg(j)+tmpyg(j+1)
          yg(j) = yg(j)/6._rkind
        enddo
      enddo
!     
      do j=1,ng
        yg(1-j) = -yg(1+j)
      enddo
!     
    endassociate
!   
  endsubroutine generate_grid_uni_extension
! 
  subroutine generate_grid_cha(self, jbgrid, dyptarget, nsmooy, retaucha)
    class(grid_object), intent(inout) :: self  !< The grid.
    real(rkind), intent(in) :: jbgrid,dyptarget,retaucha
    integer, intent(in) :: nsmooy
    real(rkind) :: ynf,rj,fsm,fsl,alfc
    real(rkind), dimension(self%nymax/2+1) :: ynold
!   real(rkind) :: alfb,ceta
    integer :: j,l
!   
!   ny/2 = 4/3/alfc * Re_tau**(3./4.)
!   
    associate(nymax=>self%nymax,  ng => self%ng, yg => self%yg, yn => self%yn, &
      domain_size => self%domain_size)
!     
!     alfb = 1.25_rkind
!     ceta =  0.8_rkind
!     alfc = alfb*ceta
      alfc = 4._rkind/3._rkind*retaucha**0.75_rkind*2._rkind/nymax
      do j=1,nymax/2+1
        rj = real(j-1,rkind)
        yn(j) = 1._rkind/(1._rkind+(rj/jbgrid)**2)
        yn(j) = yn(j)*(rj*dyptarget+(0.75_rkind*alfc*rj)**(4._rkind/3._rkind)*(rj/jbgrid)**2)
      enddo
      ynf = yn(nymax/2+1)
      do j=1,nymax/2+1
        yn(j) = yn(j)/ynf
      enddo
!     
      do l=1,nsmooy
        ynold = yn(1:nymax/2+1)
        do j=2,nymax/2
          fsm = 6._rkind-2._rkind*(tanh((yn(j)-0.75_rkind)/0.15_rkind)+1._rkind)*0.5_rkind
          fsl = 0.5_rkind*(6._rkind-fsm)
          yn(j) = fsl*ynold(j-1)+fsm*ynold(j)+fsl*ynold(j+1)
          yn(j) = yn(j)/6._rkind
        enddo
      enddo
!     
      do j=1,nymax/2+1
        yn(j) = yn(j)-1._rkind
      enddo
      do j=1,nymax/2
        yn(nymax+2-j) = -yn(j)
      enddo
!     
      do j=1,nymax
        yg(j) = 0.5_rkind*(yn(j)+yn(j+1))
      enddo
      do j=1,ng
        yg(nymax+j) =  2._rkind-yg(nymax+1-j)
        yg(1-j)     = -2._rkind-yg(j)
      enddo
!     
    endassociate
!   
  endsubroutine generate_grid_cha
! 
  subroutine write_grid(self)
    class(grid_object), intent(in) :: self  !< The grid.
    integer :: i,j,k
    associate(nxmax=>self%nxmax, nymax=>self%nymax,  nzmax=>self%nzmax, ng => self%ng,  &
      xg => self%xg, yg => self%yg, zg => self%zg, yn => self%yn, is_y_staggered => self%is_y_staggered)
!     
      open(18,file='x.dat')
      do i=1-ng,nxmax+ng+1
        write(18,*) xg(i)
      enddo
      close(18)
!     
      open(18,file='y.dat')
      do j=1-ng,nymax+ng
        write(18,*) yg(j)
      enddo
      close(18)
!     
      open(18,file='z.dat')
      do k=1-ng,nzmax+ng
        write(18,*) zg(k)
      enddo
      close(18)
!     
      if (is_y_staggered) then
        open(18,file='yn.dat')
        do j=1,nymax+1
          write(18,*) yn(j)
        enddo
        close(18)
      endif
!     
    endassociate
!   
  endsubroutine write_grid
! 
  subroutine generate_grid_uniform_x(self)
    class(grid_object), intent(inout) :: self  !< The grid.
    real(rkind) :: dx
    integer :: i
!   
    associate(nxmax=>self%nxmax, ng => self%ng,  &
      xg => self%xg, domain_size => self%domain_size, &
      is_xyz_periodic => self%is_xyz_periodic)
!     
      dx = domain_size(1)/(nxmax-1)
      if (is_xyz_periodic(1)) dx = domain_size(1)/nxmax
!     
      do i=1-ng,nxmax+ng+1
        xg(i) = (i-1)*dx
      enddo
!     
    endassociate
!   
  endsubroutine generate_grid_uniform_x
! 
  subroutine generate_grid_uniform_y(self)
    class(grid_object), intent(inout) :: self  !< The grid.
    real(rkind) :: dy
    integer :: j
!   
    associate(nymax=>self%nymax,  ng => self%ng,  &
      yg => self%yg, domain_size => self%domain_size, &
      is_xyz_periodic => self%is_xyz_periodic)
!     
      dy = domain_size(2)/(nymax-1)
      if (is_xyz_periodic(2)) dy = domain_size(2)/nymax
!
      do j=1-ng,nymax+ng
        yg(j) = (j-1)*dy
      enddo
!     
    endassociate
!   
  endsubroutine generate_grid_uniform_y
! 
  subroutine generate_grid_uniform_z(self)
    class(grid_object), intent(inout) :: self  !< The grid.
    real(rkind) :: dz
    integer :: k
!
    associate(nzmax=>self%nzmax, ng => self%ng,  &
      zg => self%zg, domain_size => self%domain_size, &
      is_xyz_periodic => self%is_xyz_periodic)
!     
      dz = domain_size(3)/(nzmax-1)
      if (is_xyz_periodic(3)) dz = domain_size(3)/nzmax
!     
      do k=1-ng,nzmax+ng
        zg(k) = (k-1)*dz
      enddo
!     
    endassociate
!   
  endsubroutine generate_grid_uniform_z
! 
  subroutine alloc(self)
    class(grid_object), intent(inout) :: self  !< The grid.
!   
    associate(nxmax=>self%nxmax, nymax=>self%nymax,  nzmax=>self%nzmax, ng => self%ng)
!     
      allocate(self%xg(1-ng:nxmax+ng+1), self%yg(1-ng:nymax+ng), self%zg(1-ng:nzmax+ng))
      allocate(self%dxg(1:nxmax), self%dyg(1:nymax), self%dzg(1:nzmax))
      allocate(self%d2xg(1:nxmax), self%d2yg(1:nymax), self%d2zg(1:nzmax))
      allocate(self%yn(1:nymax+1))
!     
    endassociate
!   
  endsubroutine alloc
! 
  subroutine read_grid(self, rebuild_ghost)
    class(grid_object), intent(inout) :: self  !< The grid.
    logical, optional, intent(in)  ::  rebuild_ghost
    logical ::  rebuild_ghost_
    integer :: i,j,k
    integer :: i1,i2,j1,j2,k1,k2
!   
    associate(nxmax=>self%nxmax, nymax=>self%nymax,  nzmax=>self%nzmax, ng => self%ng,  &
      xg => self%xg, yg => self%yg,  zg => self%zg, is_y_staggered => self%is_y_staggered, yn=>self%yn)
!     
      rebuild_ghost_ = .false.
      if(present(rebuild_ghost)) rebuild_ghost_ = rebuild_ghost
!     
      i1 = 1-ng
      i2 = nxmax+ng+1
      j1 = 1-ng
      j2 = nymax+ng
      k1 = 1-ng
      k2 = nzmax+ng
      if (rebuild_ghost_) then
        i1 = 1
        i2 = nxmax
        j1 = 1
        j2 = nymax
        k1 = 1
        k2 = nzmax
      endif
!     
      open(10,file='x.dat')
      do i=i1,i2
        read(10,*) xg(i)
      enddo
      close(10)
      open(10,file='y.dat')
      do j=j1,j2
        read(10,*) yg(j)
      enddo
      close(10)
      open(10,file='z.dat')
      do k=k1,k2
        read(10,*) zg(k)
      enddo
      close(10)
!     
      if (is_y_staggered) then
        open(10,file='yn.dat')
        do j=1,nymax+1
          read(10,*) yn(j)
        enddo
        close(10)
      endif
!     
      if (rebuild_ghost_) then
        do i=1,ng
          xg(1-i)     = 2._rkind*xg(2-i)-xg(3-i)
          xg(nxmax+i) = 2._rkind*xg(nxmax+i-1)-xg(nxmax+i-2)
        enddo
        xg(nxmax+ng+1) = 2._rkind*xg(nxmax+ng)-xg(nxmax+ng-1)
        do j=1,ng
          yg(1-j)     = 2._rkind*yg(2-j)-yg(3-j)
          yg(nymax+j) = 2._rkind*yg(nymax+j-1)-yg(nymax+j-2)
        enddo
        do k=1,ng
          zg(1-k)     = 2._rkind*zg(2-k)-zg(3-k)
          zg(nzmax+k) = 2._rkind*zg(nzmax+k-1)-zg(nzmax+k-2)
        enddo
      endif
!     
    endassociate
!   
  end subroutine read_grid
! 
endmodule streams_grid_object
