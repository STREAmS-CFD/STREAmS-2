 module comp_transform
 use parameters
 use derivatives
 contains
!
  subroutine transform(ny,y,u,naux,vaux,yi,ui,tname)
    use global_variables, only: flow_init
    implicit none
    integer, parameter  :: np = 7
    integer, intent(in) :: ny, naux
    real(rkind), dimension(ny), intent(in)  :: y, u
    real(rkind), dimension(naux,ny), intent(in) :: vaux
    character(len=20), intent(in) :: tname
    real(rkind), dimension(ny), intent(out) :: yi, ui
    real(rkind), dimension(ny) :: ff, gg
    real(rkind), dimension(np) :: yvec,fvec
    real(rkind) :: ffjm1, ggjm1, yijm1, uijm1, yjm1,ujm1
    real(rkind) :: dffjm1, d2ffjm1, dggjm1, d2ggjm1
    integer :: j
!
!   vaux(1) = rho
!   vaux(2) = nu
!   vaux(3) = mu
!
    select case(trim(tname))
    case('vanDriest')

     ff = 1._rkind
     gg = sqrt(vaux(1,:))

    case('TrettelLarsson')
!
     gg = y/(vaux(1,:)**0.5_rkind*vaux(2,:))
     call ddy(ny,gg,y(1:ny),6,ff)
     gg = ff*vaux(1,:)*vaux(2,:)
!
    case('Volpiani')
!
     ff = vaux(1,:)**0.5_rkind*vaux(3,:)**(-1.5_rkind)
     gg = vaux(1,:)**0.5_rkind*vaux(3,:)**(1._rkind-1.5_rkind)
!
    end select

    if (flow_init==0) then !Staggered mesh for channel flow
     do j=1,ny
      if (j==1) then
!      yvec = y (1:np)
!      fvec = ff(1:np)
!      call interpolate(np,fvec,yvec,0._rkind,ffjm1,dffjm1,d2ffjm1)
!      fvec = gg(1:np)
!      call interpolate(np,fvec,yvec,0._rkind,ggjm1,dggjm1,d2ggjm1)
       ffjm1 = 1._rkind
       ggjm1 = 1._rkind
       ujm1  = 0._rkind 
       yjm1  = 0._rkind
       uijm1 = 0._rkind 
       yijm1 = 0._rkind
      else
       ffjm1 = ff(j-1)
       ggjm1 = gg(j-1)
       yijm1 = yi(j-1)
       uijm1 = ui(j-1)
       ujm1  = u(j-1)
       yjm1  = y(j-1)
      endif
!
      yi(j) = yijm1 + 0.5_rkind*(ff(j) + ffjm1)*(y(j) - yjm1)
      ui(j) = uijm1 + 0.5_rkind*(gg(j) + ggjm1)*(u(j) - ujm1)
!
     enddo
    else !co-located mesh, boundary layer
     yi = 0._rkind
     ui = 0._rkind
     do j=2,ny
      yi(j) = yi(j-1) + 0.5_rkind*(ff(j) + ff(j-1))*(y(j) - y(j-1))
      ui(j) = ui(j-1) + 0.5_rkind*(gg(j) + gg(j-1))*(u(j) - u(j-1))
     enddo
    endif
  end subroutine transform
 end module 
