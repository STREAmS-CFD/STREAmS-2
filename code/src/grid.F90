!< STREAmS, grid class definition.
module streams_grid_object

  use streams_parameters
  use mpi

  implicit none
  private
  public :: grid_object
  public :: GRID_FROMFILE, GRID_UNIFORM, GRID_CHA, GRID_BL, GRID_SBLI, GRID_AIRFOIL

  integer(ikind), parameter :: GRID_FROMFILE = 1_ikind
  integer(ikind), parameter :: GRID_UNIFORM = 2_ikind
  integer(ikind), parameter :: GRID_CHA = 3_ikind
  integer(ikind), parameter :: GRID_BL = 4_ikind
  integer(ikind), parameter :: GRID_SBLI = 5_ikind
  integer(ikind), parameter :: GRID_AIRFOIL = 6_ikind

  type :: grid_object
    !< Grid class definition.
    !2 or 3, for 2D or 3D grid, grid_dim instead defines curvilinear or not
    integer :: ndim

    integer :: myrank, nprocs, mpi_err
    logical :: masterproc
    real(rkind), dimension(3) :: domain_size
    integer :: nxmax, nymax, nzmax
    integer :: ng
    logical, dimension(3) :: is_xyz_periodic = .false.
    logical :: is_y_staggered = .false.
    integer(ikind) :: grid_type
    integer(ikind) :: grid_dim
    integer(ikind) :: grid2d_par
    !grid_dim = 1 (Cartesian), 2 (2D Curvilinear), 3 (3D Curvilinear)
    real(rkind), allocatable, dimension(:) :: xg, yg, zg
    real(rkind), allocatable, dimension(:) :: dxg, dyg, dzg
    real(rkind), allocatable, dimension(:) :: d2xg, d2yg, d2zg
    real(rkind), allocatable, dimension(:) :: yn
    real(rkind), allocatable, dimension(:,:) :: xc2g, yc2g
    real(rkind), dimension(4) :: metrics_cd1
    real(rkind), dimension(0:4) :: metrics_cd2
    integer :: metrics_order
    integer :: ite, itu, ile, icl, icu ! for airfoil grid
    integer, dimension(:), allocatable :: wall_tagg
    real(rkind), dimension(:), allocatable :: xwall, ywall
    real(rkind) :: angle, R_curv, dyp_target, Retaucha ! for curved channel grid
  contains
    !public methods
    procedure, pass(self) :: alloc
    procedure, pass(self) :: compute_metrics
    procedure, pass(self) :: generate_grid_uniform_x
    procedure, pass(self) :: generate_grid_uniform_y
    procedure, pass(self) :: generate_grid_uniform_z
    procedure, pass(self) :: generate_grid_cha
    procedure, pass(self) :: generate_grid_bl
    procedure, pass(self) :: generate_grid_geoprogression
    procedure, pass(self) :: generate_grid_uni_extension
    procedure, pass(self) :: generate_yg
    procedure, pass(self) :: generate_grid_cha_curv
    procedure, pass(self) :: generate_grid_airfoil
    procedure, pass(self) :: get_deriv_coeffs
    procedure, pass(self) :: initialize
    procedure, pass(self) :: read_grid
    procedure, pass(self) :: write_grid
    procedure, pass(self) :: read_gridc2
  endtype grid_object
  !
contains
  !
  !public methods
  subroutine compute_metrics(self)
    !< Compute metrics of a block.
    !< Some refactoring should be needed. Now:
    !< Grid-1D here global metrics are computed (without considering ep_ord_change stuff) and then decomposed in field
    !< Grid-2D metrics are computed in field considering ep_ord_change stuff directly decomposed
    class(grid_object), intent(inout) :: self !< The grid.
    integer :: mm, i, j, k, l
    !
    call self%get_deriv_coeffs(self%metrics_cd1, self%metrics_cd2, self%metrics_order, self%metrics_order)
    !
    mm = self%metrics_order/2

    associate(xg => self%xg, yg => self%yg, zg => self%zg, dxg => self%dxg, dyg => self%dyg,&
    & dzg => self%dzg, d2xg => self%d2xg, d2yg => self%d2yg, d2zg => self%d2zg, nxmax => self%nxmax,&
    & nymax => self%nymax, nzmax => self%nzmax, c => self%metrics_cd1, cc => self%metrics_cd2,&
    & masterproc => self%masterproc)
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
    class(grid_object) :: self
    real(rkind), dimension(4) :: coeff_deriv1
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
      coeff_deriv1(1) = 2._rkind/3._rkind
      coeff_deriv1(2) = -1._rkind/12._rkind
    case (3)
      coeff_deriv1(1) = 0.75_rkind
      coeff_deriv1(2) = -0.15_rkind
      coeff_deriv1(3) = 1._rkind/60._rkind
    case(4)
      coeff_deriv1(1) = 4._rkind/5._rkind
      coeff_deriv1(2) = -1._rkind/5._rkind
      coeff_deriv1(3) = 4._rkind/105._rkind
      coeff_deriv1(4) = -1._rkind/280._rkind
    end select
    !
    select case (order2/2)
    case (1)
      coeff_deriv2(0) = -2._rkind
      coeff_deriv2(1) = 1._rkind
    case (2)
      coeff_deriv2(0) = -2.5_rkind
      coeff_deriv2(1) = 4._rkind/3._rkind
      coeff_deriv2(2) = -1._rkind/12._rkind
    case (3)
      coeff_deriv2(0) = -245._rkind/90._rkind
      coeff_deriv2(1) = 1.5_rkind
      coeff_deriv2(2) = -0.15_rkind
      coeff_deriv2(3) = 1._rkind/90._rkind
    case (4)
      coeff_deriv2(0) = -205._rkind/72._rkind
      coeff_deriv2(1) = 8._rkind/5._rkind
      coeff_deriv2(2) = -1._rkind/5._rkind
      coeff_deriv2(3) = 8._rkind/315._rkind
      coeff_deriv2(4) = -1._rkind/560._rkind
    endselect
    !
  endsubroutine get_deriv_coeffs
  !
  subroutine initialize(self, periodic, nxmax, nymax, nzmax, ng, grid_type, domain_size_x,&
  & domain_size_y, domain_size_z, l0, grid_vars, metrics_order, rebuild_ghost, ystaggering, grid_dim,&
  & grid2d_par)
    class(grid_object) :: self
    integer :: nxmax, nymax, nzmax, ng, grid_type, grid_dim, grid2d_par
    real(rkind) :: domain_size_x, domain_size_y, domain_size_z, l0
    integer :: nxmax_tot, nymax_tot, nzmax_tot
    integer :: metrics_order
    logical :: rebuild_ghost,ystaggering
    logical :: write_grid
    logical, dimension(3) :: periodic
    real(rkind), dimension(:), allocatable :: grid_vars
    !
    call get_mpi_basic_info(self%nprocs, self%myrank, self%masterproc, self%mpi_err)
    !
    self%nxmax = nxmax
    self%nymax = nymax
    self%nzmax = nzmax
    self%ng = ng
    self%grid_type = grid_type
    self%grid_dim = grid_dim
    self%grid2d_par = grid2d_par
    self%metrics_order = metrics_order
    self%is_xyz_periodic = periodic

    if (nzmax == 1) then
      self%ndim = 2
      if (self%masterproc) print*,'Code running in 2D mode'
    else
      self%ndim = 3
      if (self%masterproc) print*,'Code running in 3D mode'
    endif
    !
    self%domain_size = [domain_size_x,domain_size_y,domain_size_z]
    !
    call self%alloc()
    !
    write_grid = .true.
    !
    self%is_y_staggered = ystaggering
    !
    if (self%grid_dim == 1) then
      if (self%masterproc) write(*,*) '-- Cartesian grid --'

      if (self%grid_type == GRID_FROMFILE) then
        call self%read_grid(rebuild_ghost)
        write_grid = .false.
      else
        call self%generate_grid_uniform_x()
        call self%generate_grid_uniform_z()
        select case (self%grid_type)
        case(GRID_UNIFORM)
          call self%generate_grid_uniform_y()
        case(GRID_CHA)
          call self%generate_grid_cha(grid_vars(1), grid_vars(2),nint(grid_vars(3)), grid_vars(4))
        case(GRID_BL,GRID_SBLI)
          call self%generate_grid_bl(grid_vars(1), grid_vars(2), nint(grid_vars(3)), grid_vars(4), grid_vars(5))
          if (nint(grid_vars(3))<self%nymax) then
            if (self%grid_type==GRID_BL) then
              call self%generate_grid_geoprogression(nint(grid_vars(3)), grid_vars(4))
            else
              call self%generate_grid_uni_extension(nint(grid_vars(3)))
            endif
          endif
          call self%generate_yg(nint(grid_vars(6)),ystaggering)
        end select
      endif
      !
      if (self%masterproc.and.write_grid) call self%write_grid()
      !
    elseif (self%grid_dim == 2) then
      if (self%masterproc) print*,'-- Curvilinear x-y grid --'

      self%domain_size(1:2) = 1._rkind ! dummy values

      call self%generate_grid_uniform_z()
      if (self%grid_type == GRID_FROMFILE) then
        if(self%grid2d_par == 0) call self%read_gridc2()
      elseif (self%grid_type == GRID_BL) then
        if(self%grid2d_par == 0) call self%read_gridc2()
      elseif (self%grid_type == GRID_CHA) then
        self%dyp_target = grid_vars(1)
        self%Retaucha = grid_vars(2)
        self%angle = grid_vars(3)
        self%R_curv = grid_vars(4)
        if(self%grid2d_par == 0) call self%generate_grid_cha_curv(grid_vars(1), grid_vars(2), grid_vars(3), grid_vars(4))
      elseif (self%grid_type == GRID_AIRFOIL) then
        if(self%grid2d_par == 0) call self%generate_grid_airfoil()
      endif

      !xg used in recyc, y used in init, at least
      if(self%grid2d_par == 0) then
        self%xg = self%xc2g(:,1) ! 1-ng:nxmax+ng+1
        self%yg = self%yc2g(1,:)
      else
        self%xg = 0._rkind
        self%yg = 0._rkind
      endif

    endif

    if (self%ndim == 2) self%zg = 0._rkind

    self%xg = self%xg*l0
    self%yg = self%yg*l0
    self%zg = self%zg*l0
    self%yn = self%yn*l0
    self%domain_size = self%domain_size*l0
    !
    call self%compute_metrics()

    if(self%ndim == 2) then
      self%dzg = 0._rkind
      self%d2zg = 0._rkind
    endif
    !
  endsubroutine initialize
  !
  subroutine generate_grid_bl(self, jbgrid, dyptarget, nymaxwr, rlywr, retau)
    class(grid_object), intent(inout) :: self
    real(rkind), intent(in) :: jbgrid,dyptarget,rlywr,retau
    integer(ikind), intent(in) :: nymaxwr
    real(rkind) :: ynf,rj,alfc
    !real(rkind) :: alfb,ceta
    integer :: j
    !
    associate(nymax=>self%nymax, ng => self%ng, yg => self%yg, yn => self%yn, domain_size => self%domain_size)
      !
      !alfb = 1.25_rkind
      !ceta = 0.8_rkind
      !alfc = alfb*ceta
      alfc = 4._rkind/3._rkind*(retau*rlywr)**0.75_rkind/(nymaxwr+1)
      do j=1,nymaxwr+1+ng
        rj = real(j-1,rkind)
        yn(j) = 1._rkind/(1._rkind+(rj/jbgrid)**2)
        yn(j) = yn(j)*(rj*dyptarget+(0.75_rkind*alfc*rj)**(4._rkind/3._rkind)*(rj/jbgrid)**2)
      enddo
      ynf = yn(nymaxwr+1)
      do j=1,nymaxwr+ng+1
        yn(j) = yn(j)/ynf*rlywr
      enddo
      !
    endassociate
    !
  endsubroutine generate_grid_bl
  !
  subroutine generate_grid_geoprogression(self, nymaxwr, rlywr)
    class(grid_object), intent(inout) :: self !< The grid.
    real(rkind), intent(in) :: rlywr
    integer, intent(in) :: nymaxwr
    real(rkind) :: dy1, rlygp, rgp, rgpold
    integer :: nygp, j, l
    !
    associate(nymax=>self%nymax, yn => self%yn, rly => self%domain_size(2), ng => self%ng)
      !
      !dy1 = yn(nymaxwr+1)-yn(nymaxwr)
      dy1 = yn(nymaxwr+2)-yn(nymaxwr+1)
      !
      rlygp = rly-rlywr
      nygp = nymax+1-(nymaxwr+1)
      rgp = 1.01_rkind
      do
        rgpold = rgp
        rgp = 1._rkind-(1._rkind-rgp)/dy1*rlygp
        rgp = rgp**(1._rkind/nygp)
        if (abs(rgp-rgpold) < tol_iter) exit
        !write(*,*) 'Geometric progression factor rgp =', rgp
      enddo
      if (self%masterproc) write(*,*) 'Geometric progression factor rgp =', rgp
      !
      do j=1,nygp+ng
        yn(nymaxwr+1+j) = yn(nymaxwr+j) + dy1*rgp**(j-1._rkind)
      enddo
      !
    endassociate
    !
  endsubroutine generate_grid_geoprogression
  !
  subroutine generate_grid_uni_extension(self, nymaxwr)
    class(grid_object), intent(inout) :: self !< The grid.
    integer, intent(in) :: nymaxwr
    real(rkind) :: dy1
    integer :: j
    !
    associate(nymax=>self%nymax, yn => self%yn, ng => self%ng)
      !
      dy1 = yn(nymaxwr+2)-yn(nymaxwr+1)
      !
      do j=nymaxwr+2,nymax+1+ng
        yn(j) = yn(j-1)+dy1
      enddo
      !
    endassociate
    !
  endsubroutine generate_grid_uni_extension
  !
  subroutine generate_yg(self, nsmooy, ystaggering)
    class(grid_object), intent(inout) :: self !< The grid.
    integer, intent(in) :: nsmooy
    logical, intent(in) :: ystaggering
    integer :: j, l
    real(rkind), allocatable, dimension(:) :: tmpyn
    !
    associate(nymax=>self%nymax, ng => self%ng, yg => self%yg, yn => self%yn, domain_size => self%domain_size)
      !
      allocate(tmpyn(nymax+1+ng))
      !
      do l=1,nsmooy
        tmpyn = yn(1:nymax+1+ng)
        do j=2,nymax+ng
          yn(j) = tmpyn(j-1)+4._rkind*tmpyn(j)+tmpyn(j+1)
          yn(j) = yn(j)/6._rkind
        enddo
      enddo
      !
      if (ystaggering) then
        do j=1,nymax+ng
          yg(j) = 0.5_rkind*(yn(j)+yn(j+1))
        enddo
        domain_size(2) = yn(nymax+1)-yn(1)
        !ghost nodes
        do j=1,ng
          yg(1-j) = -yg(j)
        enddo
        !
      else
        do j=1,nymax+ng
          yg(j) = yn(j)
        enddo
        domain_size(2) = yg(nymax)-yg(1)
        !ghost nodes
        do j=1,ng
          yg(1-j) = -yg(1+j)
        enddo
      endif
      !
    endassociate
  endsubroutine generate_yg
  !
  subroutine generate_grid_cha(self, jbgrid, dyptarget, nsmooy, retaucha)
    class(grid_object), intent(inout) :: self !< The grid.
    real(rkind), intent(in) :: jbgrid,dyptarget,retaucha
    integer, intent(in) :: nsmooy
    real(rkind) :: ynf,rj,fsm,fsl,alfc
    real(rkind), dimension(self%nymax/2+1) :: ynold
    !real(rkind) :: alfb,ceta
    integer :: j,l
    !
    !ny/2 = 4/3/alfc * Re_tau**(3./4.)
    !
    associate(nymax=>self%nymax, ng => self%ng, yg => self%yg, yn => self%yn, domain_size => self%domain_size)
      !
      !alfb = 1.25_rkind
      !ceta = 0.8_rkind
      !alfc = alfb*ceta
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
        yg(nymax+j) = 2._rkind-yg(nymax+1-j)
        yg(1-j) = -2._rkind-yg(j)
      enddo
      !
    endassociate
    !
  endsubroutine generate_grid_cha

  subroutine generate_grid_cha_curv(self, dyp_target, retaucha, angle, R_curv)
    class(grid_object), intent(inout) :: self !< The grid.
    real(rkind), intent(in) :: dyp_target, retaucha, angle, R_curv
    real(rkind), dimension(1-self%ng:self%nxmax+self%ng) :: theta
    real(rkind), dimension(1-self%ng:self%nymax+self%ng) :: ys
    real(rkind) :: dt, theta1, deltat, dcsi, csi, b, deta, bold, eta
    integer :: i, j
    !
    associate(nxmax => self%nxmax,nymax=>self%nymax, ng => self%ng, xc2g => self%xc2g,&
    & yc2g => self%yc2g, yn => self%yn, domain_size => self%domain_size)


      !Arc of circumference in csi (stream-wise)
      dt = angle/nxmax
      theta1 = 0.5_rkind*pi+0.5_rkind*angle
      deltat = angle + dt*2*ng
      dcsi = 1._rkind/(nxmax+2*ng)
      do i=1-ng,nxmax+ng
        csi = real(i-1,rkind)*dcsi
        theta(i) = theta1 - csi * deltat
      enddo
      !
      !Tanh mapping function in eta (radial direction)
      b = 1._rkind
      deta = 1._rkind/nymax ! not nymax-1, because yn goes up to nymax+1
      do
        bold = b
        b = atanh(tanh(0.5_rkind*b)*(dyp_target/retaucha-1._rkind))/(deta-0.5_rkind)
        if (abs(b-bold) < tol_iter) exit
      enddo
      if (self%masterproc) write(*,*) 'Stretching parameter b =', b
      do j=1,nymax+1 ! so that ys(1) and ys(nymax) are both inside the domain
        eta = (real(j,rkind)-1._rkind)*deta
        yn(j) = tanh(b*(eta-0.5_rkind))/tanh(b*0.5_rkind) ! (from -1 to +1)
      enddo
      do j=1,nymax
        ys(j) = R_curv + 0.5_rkind*(yn(j)+yn(j+1)) ! staggered radius (from rly-1 to rly+1)
      enddo
      do j=1,ng ! symmetry conditions to fill ghost nodes
        ys(1-j) = 2._rkind*ys(1) - ys(1+j)
        ys(nymax+j) = 2._rkind*ys(nymax) - ys(nymax-j)
      enddo
      !
      !Composition of the 2D grid
      do j=1-ng,nymax+ng
        do i=1-ng,nxmax+ng
          xc2g(i,j) = ys(j)*cos(theta(i))
          yc2g(i,j) = ys(j)*sin(theta(i))
        enddo
      enddo

      !do j=1-ng,nymax+ng
        !do i=1-ng,nxmax+ng
          !write(21,*) xc2g(i,j), yc2g(i,j)
        !enddo
      !enddo

      if (self%masterproc) then
        write(*,*) 'Stretching parameter b =', b
        write(*,*) 'Delta x+ (centerline ) =', dt*R_curv*retaucha
        write(*,*) 'L_x (centerline) =', angle*R_curv
        write(*,*) 'theta_i =', theta(1)*180/pi, 'theta_f =', theta(nxmax)*180/pi
        open(18,file='y.dat')
        do j=1-ng,nymax+ng
          write(18,*) ys(j), yc2g(nxmax/2,j) !cannot stay here, it would be out of bounds, yn(j)
        enddo
        close(18)
        open(18,file='yn.dat') ! consistent with grid_dim=1 writings
        do j=1,nymax+1
          write(18,*) yn(j)
        enddo
        close(18)
      endif

    endassociate
    !
  endsubroutine generate_grid_cha_curv

  subroutine generate_grid_airfoil(self)
    class(grid_object), intent(inout) :: self !< The grid.
    real(rkind), dimension(1:self%nxmax) :: xx1
    real(rkind) :: tol, xtrip, errmin, error, x0tr, y0tr, x0ts, y0ts
    integer :: i, j, nxg, nyg, itr, its, l
    integer :: nxmax_read,nymax_read,nzmax_one
    logical :: file_exists
    !
    associate(nxmax => self%nxmax,nymax=>self%nymax, ng => self%ng, xc2g => self%xc2g,&
    & yc2g => self%yc2g, yn => self%yn, domain_size => self%domain_size, masterproc => self%masterproc,&
    & ite => self%ite, itu => self%itu, ile => self%ile, icl => self%icl, icu => self%icu,&
    & wall_tagg => self%wall_tagg)

      nxg = nxmax+ng
      nyg = nymax+ng
      tol = 1E-10
      !
      !Read external grid (warning: the grid goes from j=1 to j=nymax+ng)
      inquire(file='grid2d.xyz',exist=file_exists)
      if (file_exists) then
        open(10,file='grid2d.xyz',form='unformatted',access='stream')
        read(10) nxmax_read,nymax_read,nzmax_one
        if(nxmax_read /= nxmax) then
          print*,'Error! Grid size along x does not match: ',nxmax_read,nxmax
          call mpi_abort(mpi_comm_world,88,self%mpi_err)
        endif
        if(nymax_read /= nymax) then
          print*,'Error! Grid size along y does not match: ',nymax_read,nymax
          call mpi_abort(mpi_comm_world,88,self%mpi_err)
        endif
        if(nzmax_one /= 1) then
          print*,'Error! Grid along z from file must be one: ',nzmax_one
          call mpi_abort(mpi_comm_world,88,self%mpi_err)
        endif
        do j=1,nymax
          do i=1,nxmax
            read(10) xc2g(i,j)
          enddo
        enddo
        do j=1,nymax
          do i=1,nxmax
            read(10) yc2g(i,j)
          enddo
        enddo
        close(10)
      else
        open(10,file='grid_xy.dat')
        do j=1,nymax ! nyg
          do i=1,nxmax
            read(10,*) xc2g(i,j),yc2g(i,j)
          enddo
        enddo
        close(10)
      endif
      !print*,'x-minmax: ',minval(xc2g(:,:)),maxval(xc2g(:,:))
      !print*,'y-minmax: ',minval(yc2g(:,:)),maxval(yc2g(:,:))
      !if(masterproc) then
        !open(11,file='grid_streams_start.xyz',form='unformatted')
        !write(11) nxmax,nymax
        !write(11) (((xc2g(i,j), i=1,nxmax), j=1,nymax), l=1,1), ! (((yc2g(i,j), i=1,nxmax), j=1,nymax), l=1,1)
        !close(11)
      !endif
      !
      !Set ghost nodes in i direction (extrapolation)
      do j=1,nymax ! nyg
        do i=1,ng

          xc2g(1-i,j) = 2._rkind*xc2g(2-i,j)-xc2g(3-i,j)
          yc2g(1-i,j) = 2._rkind*yc2g(2-i,j)-yc2g(3-i,j)
          !mirror xc2g(1-i,j) = 2._rkind*xc2g(1,j)-xc2g(1+i,j)
          !mirror yc2g(1-i,j) = 2._rkind*yc2g(1,j)-yc2g(1+i,j)

          xc2g(nxmax+i,j) = 2._rkind*xc2g(nxmax+i-1,j)-xc2g(nxmax+i-2,j)
          yc2g(nxmax+i,j) = 2._rkind*yc2g(nxmax+i-1,j)-yc2g(nxmax+i-2,j)
        enddo
      enddo
      !
      !Set ghost nodes in j direction at outer boundary
      do j=1,ng
        do i=1-ng,nxg
          xc2g(i,nymax+j) = 2._rkind*xc2g(i,nymax+j-1)-xc2g(i,nymax+j-2)
          yc2g(i,nymax+j) = 2._rkind*yc2g(i,nymax+j-1)-yc2g(i,nymax+j-2)
        enddo
      enddo

      !Find the location of the trailing edge
      i=1
      do while (i<=nxmax)
        if (abs(yc2g(i,1)-yc2g(nxmax-i+1,1))<tol) then
          i=i+1
        else
          exit
        endif
      enddo
      ite = i-1
      itu = nxmax-ite+1
      !
      !Find the location of the leading edge
      do i=1,nxmax
        xx1(i) = xc2g(i,1)
      enddo
      ile = minloc(xx1,1)
      !
      if (masterproc) then
        print*, '======================================'
        print*, ' LE, TE, Tripping, Actuation '
        print*, '======================================'
        print*, 'Global leading edge node ', ile
        print*, 'Global trailing edge nodes ', ite,itu
      endif

      !Setting ghost nodes in j direction (wake-->swap)
      do j=1,ng
        do i=1-ng,nxg
          xc2g(i,1-j) = xc2g(nxmax-i+1,1+j)
          yc2g(i,1-j) = yc2g(nxmax-i+1,1+j)
        enddo
      enddo
      !
      !Setting ghost nodes in j direction (airfoil-->extrapolation)
      do j=1,ng
        do i=ite+1,itu-1
          !extrapolation xc2g(i,1-j) = 2._rkind*xc2g(i,2-j)-xc2g(i,3-j)
          !extrapolation yc2g(i,1-j) = 2._rkind*yc2g(i,2-j)-yc2g(i,3-j)
          xc2g(i,1-j) = 2._rkind*xc2g(i,1)-xc2g(i,1+j)
          yc2g(i,1-j) = 2._rkind*yc2g(i,1)-yc2g(i,1+j)
        enddo
      enddo

      !if(masterproc) then
        !open(11,file='grid_streams.xyz',form='unformatted')
        !write(11) nxmax+2*ng,nymax+2*ng
        !write(11) (((xc2g(i,j), i=1-ng,nxmax+ng), j=1-ng,nymax+ng), l=1,1), ! (((yc2g(i,j), i=1-ng,
!        nxmax+ng), j=1-ng,nymax+ng), l=1,1)
        !close(11)
      !endif

      !
      if (masterproc) then
        open(18,file='wall_indices.dat')
        write(18,*) ite,ile,itu
        close(18)
      endif

      !wall_tag is used to distinguish airfoil wall from wake
      wall_tagg = 1
      do i=ite,ile
        wall_tagg(i) = -1 ! pressure side
      enddo
      do i=ile+1,itu
        wall_tagg(i) = 0 ! suction side
      enddo
      wall_tagg(ite) = -2
      wall_tagg(itu) = -3
      !
      if (masterproc) then
        open(18,file='tagged_wall.dat')
        do i=1-ng,nxmax+ng
          write(18,*) i,xc2g(i,1),yc2g(i,1),wall_tagg(i)
        enddo
        close(18)
      endif

      !Find the global index of c/2 on the two sides
      !test icl = minloc(abs(xc2g(ite:ile,1)-0.5_rkind))
      !test icu = minloc(abs(xc2g(ile:itu,1)-0.5_rkind))
      !test print*,'icl, icu V1 :',icl, icu

      errmin = 100._rkind
      do i=ite,ile ! pressure side of the airfoil
        error = abs(xc2g(i,1)-.5_rkind)
        if (error < errmin) then
          errmin = error
          icl = i
        endif
      enddo
      errmin = 100._rkind
      do i=ile,itu ! suction side of the airfoil
        error = abs(xc2g(i,1)-.5_rkind)
        if (error < errmin) then
          errmin = error
          icu = i
        endif
      enddo
      !print*,'icl, icu V2 :',icl, icu

    endassociate
    !
  endsubroutine generate_grid_airfoil
  !
  subroutine write_grid(self)
    class(grid_object), intent(in) :: self !< The grid.
    integer :: i,j,k
    associate(nxmax=>self%nxmax, nymax=>self%nymax, nzmax=>self%nzmax, ng => self%ng,&
    & domain_size => self%domain_size, xg => self%xg, yg => self%yg, zg => self%zg, yn => self%yn,&
    & is_y_staggered => self%is_y_staggered)
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
          write(18,*) yn(j),yn(j+1)-yn(j)
        enddo
        close(18)
      endif
      !
      open(18,file='domain_size.dat')
      write(18,*) domain_size(1:3)
      close(18)
      !
    endassociate
    !
  endsubroutine write_grid
  !
  subroutine generate_grid_uniform_x(self)
    class(grid_object), intent(inout) :: self !< The grid.
    real(rkind) :: dx
    integer :: i
    !
    associate(nxmax=>self%nxmax, ng => self%ng, xg => self%xg, domain_size => self%domain_size,&
    & is_xyz_periodic => self%is_xyz_periodic)
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
    class(grid_object), intent(inout) :: self !< The grid.
    real(rkind) :: dy
    integer :: j
    !
    associate(nymax=>self%nymax, ng => self%ng, yg => self%yg, domain_size => self%domain_size,&
    & is_xyz_periodic => self%is_xyz_periodic)
      !
      dy = domain_size(2)/(nymax-1)
      if (is_xyz_periodic(2)) dy = domain_size(2)/nymax

      do j=1-ng,nymax+ng
        yg(j) = (j-1)*dy
      enddo
      !
    endassociate
    !
  endsubroutine generate_grid_uniform_y
  !
  subroutine generate_grid_uniform_z(self)
    class(grid_object), intent(inout) :: self !< The grid.
    real(rkind) :: dz
    integer :: k

    associate(nzmax=>self%nzmax, ng => self%ng, zg => self%zg, domain_size => self%domain_size,&
    & is_xyz_periodic => self%is_xyz_periodic)
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
    class(grid_object), intent(inout) :: self !< The grid.
    !
    associate(nxmax=>self%nxmax, nymax=>self%nymax, nzmax=>self%nzmax, ng => self%ng)
      !
      allocate(self%xg(1-ng:nxmax+ng+1), self%yg(1-ng:nymax+ng), self%zg(1-ng:nzmax+ng))
      allocate(self%dxg(1:nxmax), self%dyg(1:nymax), self%dzg(1:nzmax))
      allocate(self%d2xg(1:nxmax), self%d2yg(1:nymax), self%d2zg(1:nzmax))
      allocate(self%yn(1:nymax+1+ng))
      !
      if (self%grid_dim == 2) then
        allocate(self%xwall(1-ng:nxmax+ng+1), self%ywall(1-ng:nxmax+ng+1))
        allocate(self%wall_tagg(1-ng:nxmax+ng))
        self%wall_tagg(:) = 0 ! default used for non-airfoil case
        if (self%grid2d_par == 0) then
          allocate(self%xc2g(1-ng:nxmax+ng+1, 1-ng:nymax+ng))
          allocate(self%yc2g(1-ng:nxmax+ng+1, 1-ng:nymax+ng))
        endif
      endif
      !
    endassociate
    !
  endsubroutine alloc
  !
  subroutine read_grid(self, rebuild_ghost)
    class(grid_object), intent(inout) :: self !< The grid.
    logical, optional, intent(in) :: rebuild_ghost
    logical :: rebuild_ghost_
    integer :: i,j,k
    integer :: i1,i2,j1,j2,k1,k2
    logical :: file_exists
    !
    associate(nxmax=>self%nxmax, nymax=>self%nymax, nzmax=>self%nzmax, ng => self%ng,&
    & domain_size => self%domain_size, xg => self%xg, yg => self%yg, zg => self%zg, is_y_staggered => sel&
    &f%is_y_staggered, yn=>self%yn)
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
      inquire(file='domain_size.dat',exist=file_exists)
      if (file_exists) then
        open(18,file='domain_size.dat')
        read(18,*) domain_size(1:3)
        close(18)
      else
        domain_size(1) = xg(nxmax)-xg(1)
        domain_size(2) = yg(nymax)-yg(1)
        domain_size(3) = zg(nzmax)-zg(1)
      endif
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
          xg(1-i) = 2._rkind*xg(2-i)-xg(3-i)
          xg(nxmax+i) = 2._rkind*xg(nxmax+i-1)-xg(nxmax+i-2)
        enddo
        xg(nxmax+ng+1) = 2._rkind*xg(nxmax+ng)-xg(nxmax+ng-1)
        do j=1,ng
          yg(1-j) = 2._rkind*yg(2-j)-yg(3-j)
          yg(nymax+j) = 2._rkind*yg(nymax+j-1)-yg(nymax+j-2)
        enddo
        do k=1,ng
          zg(1-k) = 2._rkind*zg(2-k)-zg(3-k)
          zg(nzmax+k) = 2._rkind*zg(nzmax+k-1)-zg(nzmax+k-2)
        enddo
      endif
      !
    endassociate
    !
  end subroutine read_grid
  !
  subroutine old_read_gridc2(self, read_z)
    class(grid_object), intent(inout) :: self !< The grid.
    integer :: i,j,k, k1, k2
    logical, intent(in), optional :: read_z
    logical :: rebuild_ghost_, read_z_
    !
    read_z_ = .false. ; if(present(read_z)) read_z_ = read_z

    associate(nxmax => self%nxmax, nymax => self%nymax, nzmax => self%nzmax, ng => self%ng,&
    & xc2g => self%xc2g, yc2g => self%yc2g, zg => self%zg)
      !
      open(10,file='grid2d.dat') ! in i grid goes up to nxmax+ng+1
      do j=1-ng,nymax+ng
        do i=1-ng,nxmax+ng+1
          read(10,*) xc2g(i,j),yc2g(i,j)
        enddo
      enddo
      close(10)
      !
      if(read_z_) then
        rebuild_ghost_ = .false.
        !
        k1 = 1-ng
        k2 = nzmax+ng
        if (rebuild_ghost_) then
          k1 = 1
          k2 = nzmax
        endif
        !
        open(10,file='z.dat')
        do k=k1,k2
          read(10,*) zg(k)
        enddo
        close(10)
        !
        if (rebuild_ghost_) then
          do k=1,ng
            zg(1-k) = 2._rkind*zg(2-k)-zg(3-k)
            zg(nzmax+k) = 2._rkind*zg(nzmax+k-1)-zg(nzmax+k-2)
          enddo
        endif
      endif
      !
    endassociate
    !
  end subroutine old_read_gridc2

  subroutine read_gridc2(self, read_z)
    class(grid_object), intent(inout) :: self !< The grid.
    integer :: i,j,k, k1, k2
    logical, intent(in), optional :: read_z
    logical :: rebuild_ghost_, read_z_
    integer :: nxmax_read,nymax_read,nzmax_one
    !
    read_z_ = .false. ; if(present(read_z)) read_z_ = read_z

    associate(nxmax => self%nxmax, nymax => self%nymax, nzmax => self%nzmax, ng => self%ng,&
    & xc2g => self%xc2g, yc2g => self%yc2g, zg => self%zg)
      !
      open(10,file='grid2d.xyz',form='unformatted',access='stream')
      read(10) nxmax_read,nymax_read,nzmax_one
      if(nxmax_read /= nxmax) then
        print*,'Error! Grid size along x does not match: ',nxmax_read,nxmax
        call mpi_abort(mpi_comm_world,88,self%mpi_err)
      endif
      if(nymax_read /= nymax) then
        print*,'Error! Grid size along y does not match: ',nymax_read,nymax
        call mpi_abort(mpi_comm_world,88,self%mpi_err)
      endif
      if(nzmax_one /= 1) then
        print*,'Error! Grid along z from file must be one: ',nzmax_one
        call mpi_abort(mpi_comm_world,88,self%mpi_err)
      endif
      do j=1,nymax
        do i=1,nxmax
          read(10) xc2g(i,j)
        enddo
      enddo
      do j=1,nymax
        do i=1,nxmax
          read(10) yc2g(i,j)
        enddo
      enddo
      close(10)
      print*,'CHECK x/y min: ',minval(xc2g),minval(yc2g)
      print*,'CHECK x/y max: ',maxval(xc2g),maxval(yc2g)

      if (self%is_xyz_periodic(1)) then
        !Set ghost nodes in i direction (periodicity)
        do j=1,nymax ! nyg
          do i=1,ng
            xc2g(1-i,j) = xc2g(nxmax-i+1,j)
            yc2g(1-i,j) = yc2g(nxmax-i+1,j)
            xc2g(nxmax+i,j) = xc2g(i,j)
            yc2g(nxmax+i,j) = yc2g(i,j)
          enddo
          xc2g(nxmax+ng+1,j) = xc2g(ng+1,j)
          yc2g(nxmax+ng+1,j) = yc2g(ng+1,j)
        enddo
      else
        !Set ghost nodes in i direction (extrapolation)
        do j=1,nymax ! nyg
          do i=1,ng
            xc2g(1-i,j) = 2._rkind*xc2g(2-i,j)-xc2g(3-i,j)
            yc2g(1-i,j) = 2._rkind*yc2g(2-i,j)-yc2g(3-i,j)
            xc2g(nxmax+i,j) = 2._rkind*xc2g(nxmax+i-1,j)-xc2g(nxmax+i-2,j)
            yc2g(nxmax+i,j) = 2._rkind*yc2g(nxmax+i-1,j)-yc2g(nxmax+i-2,j)
          enddo
          xc2g(nxmax+ng+1,j) = 2._rkind*xc2g(nxmax+ng,j)-xc2g(nxmax+ng-1,j)
          yc2g(nxmax+ng+1,j) = 2._rkind*yc2g(nxmax+ng,j)-yc2g(nxmax+ng-1,j)
        enddo
      endif
      !
      !Set ghost nodes in j direction at outer boundary
      do j=1,ng
        do i=1-ng,nxmax+ng+1
          !extrapolation xc2g(i,1-j) = 2._rkind*xc2g(i,2-j)-xc2g(i,3-j)
          !extrapolation yc2g(i,1-j) = 2._rkind*yc2g(i,2-j)-yc2g(i,3-j)
          xc2g(i,1-j) = 2._rkind*xc2g(i,1)-xc2g(i,1+j)
          yc2g(i,1-j) = 2._rkind*yc2g(i,1)-yc2g(i,1+j)
          xc2g(i,nymax+j) = 2._rkind*xc2g(i,nymax+j-1)-xc2g(i,nymax+j-2)
          yc2g(i,nymax+j) = 2._rkind*yc2g(i,nymax+j-1)-yc2g(i,nymax+j-2)
        enddo
      enddo
      !
      if(read_z_) then
        rebuild_ghost_ = .false.
        !
        k1 = 1-ng
        k2 = nzmax+ng
        if (rebuild_ghost_) then
          k1 = 1
          k2 = nzmax
        endif
        !
        open(10,file='z.dat')
        do k=k1,k2
          read(10,*) zg(k)
        enddo
        close(10)
        !
        if (rebuild_ghost_) then
          do k=1,ng
            zg(1-k) = 2._rkind*zg(2-k)-zg(3-k)
            zg(nzmax+k) = 2._rkind*zg(nzmax+k-1)-zg(nzmax+k-2)
          enddo
        endif
      endif
      !
    endassociate
    !
  end subroutine read_gridc2
  !
endmodule streams_grid_object

