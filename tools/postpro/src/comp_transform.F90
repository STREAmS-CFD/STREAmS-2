 module comp_transform
 use parameters
 use derivatives
 contains
!
  subroutine transform(ystag,ny,y,u,naux,vaux,yi,ui,tname)
    use global_variables, only: flow_init
    implicit none
    integer, intent(in) :: ny, naux, ystag
    real(rkind), dimension(ny), intent(in)  :: y, u
    real(rkind), dimension(naux,ny), intent(in) :: vaux
    character(len=20), intent(in) :: tname
    real(rkind), dimension(ny), intent(out) :: yi, ui
    real(rkind), dimension(ny) :: ff, gg, muirat, mucrat, Dc, Di
    real(rkind) :: ffjm1,ggjm1,yijm1,uijm1,yjm1,ujm1,Aplus,vkc
    integer :: j
!
!   vaux(1) = rho/rho_w
!   vaux(2) = nu/nu_w
!   vaux(3) = mu/mu_w
!   vaux(4) = yplus
!   vaux(5) = ystar
!   vaux(6) = Mtau  (only one number)
!
    select case(trim(tname))
    case('vanDriest')

     ff = 1._rkind
     gg = sqrt(vaux(1,:))

    case('TrettelLarsson')
!
     gg = y/(vaux(1,:)**0.5_rkind*vaux(2,:))
     call ddy(ny,gg,y(1:ny),2,ff)
     gg = ff*vaux(1,:)*vaux(2,:)
!
    case('Volpiani')
!
     ff = vaux(1,:)**0.5_rkind*(vaux(1,:)*vaux(2,:))**(-1.5_rkind)
     gg = vaux(1,:)**0.5_rkind*(vaux(1,:)*vaux(2,:))**(-0.5_rkind)
    case('Hasan')
!
     gg = y/(vaux(1,:)**0.5_rkind*vaux(2,:))
     call ddy(ny,gg,y(1:ny),6,ff)
     gg = ff*vaux(1,:)*vaux(2,:)
     vkc = 0.41_rkind
     Aplus = 17._rkind
     Di(:) = (1._rkind-exp(-vaux(5,:)/Aplus))**2
     Aplus = 17._rkind + 19.3_rkind*vaux(6,1)
     Dc(:) = (1._rkind-exp(-vaux(5,:)/Aplus))**2
     mucrat(:) = vkc * Dc(:) * vaux(5,:)
     muirat(:) = vkc * Di(:) * vaux(5,:)
     gg = gg * (1._rkind+mucrat)/(1._rkind+muirat)
!
    end select

    if (ystag>0) then

     do j=1,ny

      if (j==1) then
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
