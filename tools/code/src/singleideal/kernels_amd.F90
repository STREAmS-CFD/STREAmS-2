module streams_kernels_amd

  use streams_parameters, only : rkind, ikind, real64, c_rkind
  use hipfort
  use hipfort_check
  implicit none

  interface
    subroutine zero_flux_kernel_wrapper(stream,nx,ny,nz,nv,fl_gpu)bind(c,name="zero_flux_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,nv
      type(c_ptr), value :: fl_gpu

    endsubroutine zero_flux_kernel_wrapper
  endinterface

  interface
    subroutine init_flux_kernel_wrapper(stream,nx,ny,nz,nv,rhodt,fl_gpu,fln_gpu)bind(c,name="init_flux_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,nv
      real(c_rkind), value :: rhodt
      type(c_ptr), value :: fl_gpu,fln_gpu

    endsubroutine init_flux_kernel_wrapper
  endinterface

  interface
    subroutine count_weno_kernel1_wrapper(stream,nv,nv_aux,nx,ny,nz,ng,eul_imin,eul_imax,eul_jmin,&
    &eul_jmax,eul_kmin,eul_kmax,weno_scheme,sensor_threshold,ep_ord_change_gpu,w_aux_gpu,count_weno_x,&
    &redn_3d_gpu)bind(c,name="count_weno_kernel1_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nv,nv_aux,nx,ny,nz,ng,eul_imin,eul_imax,eul_jmin,eul_jmax,eul_kmin,eul_kmax,weno_scheme
      real(c_rkind), value :: sensor_threshold
      type(c_ptr), value :: ep_ord_change_gpu
      type(c_ptr), value :: w_aux_gpu
      real(c_rkind) :: count_weno_x
      type(c_ptr), value :: redn_3d_gpu

    endsubroutine count_weno_kernel1_wrapper
  endinterface

  interface
    subroutine count_weno_kernel2_wrapper(stream,nv,nv_aux,nx,ny,nz,ng,eul_imin,eul_imax,eul_jmin,&
    &eul_jmax,eul_kmin,eul_kmax,weno_scheme,sensor_threshold,ep_ord_change_gpu,w_aux_gpu,count_weno_y,&
    &redn_3d_gpu)bind(c,name="count_weno_kernel2_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nv,nv_aux,nx,ny,nz,ng,eul_imin,eul_imax,eul_jmin,eul_jmax,eul_kmin,eul_kmax,weno_scheme
      real(c_rkind), value :: sensor_threshold
      type(c_ptr), value :: ep_ord_change_gpu
      type(c_ptr), value :: w_aux_gpu
      real(c_rkind) :: count_weno_y
      type(c_ptr), value :: redn_3d_gpu

    endsubroutine count_weno_kernel2_wrapper
  endinterface

  interface
    subroutine count_weno_kernel3_wrapper(stream,nv,nv_aux,nx,ny,nz,ng,eul_imin,eul_imax,eul_jmin,&
    &eul_jmax,eul_kmin,eul_kmax,weno_scheme,sensor_threshold,ep_ord_change_gpu,w_aux_gpu,count_weno_z,&
    &redn_3d_gpu)bind(c,name="count_weno_kernel3_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nv,nv_aux,nx,ny,nz,ng,eul_imin,eul_imax,eul_jmin,eul_jmax,eul_kmin,eul_kmax,weno_scheme
      real(c_rkind), value :: sensor_threshold
      type(c_ptr), value :: ep_ord_change_gpu
      type(c_ptr), value :: w_aux_gpu
      real(c_rkind) :: count_weno_z
      type(c_ptr), value :: redn_3d_gpu

    endsubroutine count_weno_kernel3_wrapper
  endinterface

  interface
    subroutine count_weno_c2_kernel1_wrapper(stream,nv,nv_aux,nx,ny,nz,ng,eul_imin,eul_imax,&
    &eul_jmin,eul_jmax,eul_kmin,eul_kmax,l_base,weno_scheme,sensor_threshold,ep_ord_change_gpu,&
    &lmax_tag_gpu,wall_tag_gpu,w_aux_gpu,count_weno_x,redn_3d_gpu)bind(c,name="count_weno_c2_kernel1_wrap&
    &per")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nv,nv_aux,nx,ny,nz,ng,eul_imin,eul_imax,eul_jmin,eul_jmax,eul_kmin,eul_kmax,l_base,weno_scheme
      real(c_rkind), value :: sensor_threshold
      type(c_ptr), value :: ep_ord_change_gpu,lmax_tag_gpu,wall_tag_gpu
      type(c_ptr), value :: w_aux_gpu
      real(c_rkind) :: count_weno_x
      type(c_ptr), value :: redn_3d_gpu

    endsubroutine count_weno_c2_kernel1_wrapper
  endinterface

  interface
    subroutine count_weno_c2_kernel2_wrapper(stream,nv,nv_aux,nx,ny,nz,ng,eul_imin,eul_imax,&
    &eul_jmin,eul_jmax,eul_kmin,eul_kmax,l_base,weno_scheme,sensor_threshold,ep_ord_change_gpu,&
    &lmax_tag_gpu,wall_tag_gpu,w_aux_gpu,count_weno_y,redn_3d_gpu)bind(c,name="count_weno_c2_kernel2_wrap&
    &per")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nv,nv_aux,nx,ny,nz,ng,eul_imin,eul_imax,eul_jmin,eul_jmax,eul_kmin,eul_kmax,l_base,weno_scheme
      real(c_rkind), value :: sensor_threshold
      type(c_ptr), value :: ep_ord_change_gpu,lmax_tag_gpu,wall_tag_gpu
      type(c_ptr), value :: w_aux_gpu
      real(c_rkind) :: count_weno_y
      type(c_ptr), value :: redn_3d_gpu

    endsubroutine count_weno_c2_kernel2_wrapper
  endinterface

  interface
    subroutine count_weno_c2_kernel3_wrapper(stream,nv,nv_aux,nx,ny,nz,ng,eul_imin,eul_imax,&
    &eul_jmin,eul_jmax,eul_kmin,eul_kmax,l_base,weno_scheme,sensor_threshold,ep_ord_change_gpu,&
    &lmax_tag_gpu,wall_tag_gpu,w_aux_gpu,count_weno_z,redn_3d_gpu)bind(c,name="count_weno_c2_kernel3_wrap&
    &per")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nv,nv_aux,nx,ny,nz,ng,eul_imin,eul_imax,eul_jmin,eul_jmax,eul_kmin,eul_kmax,l_base,weno_scheme
      real(c_rkind), value :: sensor_threshold
      type(c_ptr), value :: ep_ord_change_gpu,lmax_tag_gpu,wall_tag_gpu
      type(c_ptr), value :: w_aux_gpu
      real(c_rkind) :: count_weno_z
      type(c_ptr), value :: redn_3d_gpu

    endsubroutine count_weno_c2_kernel3_wrapper
  endinterface

  interface
    subroutine euler_x_update_kernel_wrapper(stream,nx,ny,nz,ng,nv,eul_imin,eul_imax,fhat_gpu,&
    &fl_gpu,dcsidx_gpu)bind(c,name="euler_x_update_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,eul_imin,eul_imax
      type(c_ptr), value :: fhat_gpu,fl_gpu,dcsidx_gpu

    endsubroutine euler_x_update_kernel_wrapper
  endinterface

  interface
    subroutine euler_x_update_c2_kernel_wrapper(stream,nx,ny,nz,ng,nv,eul_imin,eul_imax,fhat_gpu,&
    &fl_gpu,jac_gpu)bind(c,name="euler_x_update_c2_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,eul_imin,eul_imax
      type(c_ptr), value :: fhat_gpu,fl_gpu,jac_gpu

    endsubroutine euler_x_update_c2_kernel_wrapper
  endinterface

  interface
    subroutine euler_y_update_kernel_wrapper(stream,nx,ny,nz,ng,nv,eul_jmin,eul_jmax,fhat_gpu,&
    &fl_gpu,detady_gpu)bind(c,name="euler_y_update_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,eul_jmin,eul_jmax
      type(c_ptr), value :: fhat_gpu,fl_gpu,detady_gpu

    endsubroutine euler_y_update_kernel_wrapper
  endinterface

  interface
    subroutine euler_y_update_c2_kernel_wrapper(stream,nx,ny,nz,ng,nv,eul_jmin,eul_jmax,&
    &wall_tag_gpu,fhat_gpu,fl_gpu,jac_gpu)bind(c,name="euler_y_update_c2_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,eul_jmin,eul_jmax
      type(c_ptr), value :: wall_tag_gpu
      type(c_ptr), value :: fhat_gpu,fl_gpu,jac_gpu

    endsubroutine euler_y_update_c2_kernel_wrapper
  endinterface

  interface
    subroutine euler_z_update_kernel_wrapper(stream,nx,ny,nz,ng,nv,eul_kmin,eul_kmax,fhat_gpu,&
    &fl_gpu,dzitdz_gpu)bind(c,name="euler_z_update_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,eul_kmin,eul_kmax
      type(c_ptr), value :: fhat_gpu,fl_gpu,dzitdz_gpu

    endsubroutine euler_z_update_kernel_wrapper
  endinterface

  interface
    subroutine force_rhs_2_kernel_wrapper(stream,nx,ny,nz,ng,bulk_1,bulk_2,fluid_mask_gpu,fln_gpu,&
    &w_aux_gpu)bind(c,name="force_rhs_2_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng
      real(c_rkind), value :: bulk_1,bulk_2
      type(c_ptr), value :: fluid_mask_gpu
      type(c_ptr), value :: fln_gpu,w_aux_gpu

    endsubroutine force_rhs_2_kernel_wrapper
  endinterface

  interface
    subroutine force_rhs_2_c2_kernel_wrapper(stream,nx,ny,nz,ng,r_curv,bulk_1,bulk_5,fluid_mask_gpu,&
    &fln_gpu,yn_gpu,dxdcsinc2_gpu,dydcsinc2_gpu,w_aux_gpu)bind(c,name="force_rhs_2_c2_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng
      real(c_rkind), value :: r_curv,bulk_1,bulk_5
      type(c_ptr), value :: fluid_mask_gpu
      type(c_ptr), value :: fln_gpu,yn_gpu,dxdcsinc2_gpu,dydcsinc2_gpu,w_aux_gpu

    endsubroutine force_rhs_2_c2_kernel_wrapper
  endinterface

  interface
    subroutine force_rhs_1_kernel_wrapper(stream,nx,ny,nz,ng,fluid_mask_gpu,yn_gpu,fln_gpu,w_gpu,&
    &w_aux_gpu,bulk_1,bulk_2,bulk_3,bulk_4,bulk_5,redn_3d_gpu)bind(c,name="force_rhs_1_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng
      type(c_ptr), value :: fluid_mask_gpu
      type(c_ptr), value :: yn_gpu,fln_gpu,w_gpu,w_aux_gpu
      real(c_rkind) :: bulk_1,bulk_2,bulk_3,bulk_4,bulk_5
      type(c_ptr), value :: redn_3d_gpu

    endsubroutine force_rhs_1_kernel_wrapper
  endinterface

  interface
    subroutine force_rhs_1_c2_kernel_wrapper(stream,nx,ny,nz,ng,fluid_mask_gpu,yn_gpu,fln_gpu,&
    &dcsidxnc2_gpu,dcsidync2_gpu,jac_gpu,w_gpu,w_aux_gpu,bulk_1,bulk_2,bulk_3,bulk_4,bulk_5,&
    &redn_3d_gpu)bind(c,name="force_rhs_1_c2_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng
      type(c_ptr), value :: fluid_mask_gpu
      type(c_ptr), value :: yn_gpu,fln_gpu,dcsidxnc2_gpu,dcsidync2_gpu,jac_gpu,w_gpu,w_aux_gpu
      real(c_rkind) :: bulk_1,bulk_2,bulk_3,bulk_4,bulk_5
      type(c_ptr), value :: redn_3d_gpu

    endsubroutine force_rhs_1_c2_kernel_wrapper
  endinterface

  interface
    subroutine force_var_1_kernel_wrapper(stream,nx,ny,nz,ng,indx_cp_l,indx_cp_r,&
    &calorically_perfect,t0,tol_iter_nr,bulkt,fluid_mask_gpu,yn_gpu,fln_gpu,w_gpu,w_aux_gpu,cv_coeff_gpu,&
    &bulk_5,redn_3d_gpu)bind(c,name="force_var_1_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect
      real(c_rkind), value :: t0,tol_iter_nr,bulkt
      type(c_ptr), value :: fluid_mask_gpu
      type(c_ptr), value :: yn_gpu,fln_gpu,w_gpu,w_aux_gpu,cv_coeff_gpu
      real(c_rkind) :: bulk_5
      type(c_ptr), value :: redn_3d_gpu

    endsubroutine force_var_1_kernel_wrapper
  endinterface

  interface
    subroutine force_var_2_kernel_wrapper(stream,nx,ny,nz,ng,indx_cp_l,indx_cp_r,&
    &calorically_perfect,tbdiff,t0,fluid_mask_gpu,w_gpu,w_aux_gpu,cv_coeff_gpu)bind(c,&
    &name="force_var_2_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect
      real(c_rkind), value :: tbdiff,t0
      type(c_ptr), value :: fluid_mask_gpu
      type(c_ptr), value :: w_gpu,w_aux_gpu,cv_coeff_gpu

    endsubroutine force_var_2_kernel_wrapper
  endinterface

  interface
    subroutine force_var_1_c2_kernel_wrapper(stream,nx,ny,nz,ng,indx_cp_l,indx_cp_r,&
    &calorically_perfect,t0,tol_iter_nr,bulkt,fluid_mask_gpu,yn_gpu,fln_gpu,w_gpu,w_aux_gpu,cv_coeff_gpu,&
    &jac_gpu,dcsidxnc2_gpu,dcsidync2_gpu,bulk_5,redn_3d_gpu)bind(c,name="force_var_1_c2_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect
      real(c_rkind), value :: t0,tol_iter_nr,bulkt
      type(c_ptr), value :: fluid_mask_gpu
      type(c_ptr), value :: yn_gpu,fln_gpu,w_gpu,w_aux_gpu,cv_coeff_gpu,jac_gpu,dcsidxnc2_gpu,dcsidync2_gpu
      real(c_rkind) :: bulk_5
      type(c_ptr), value :: redn_3d_gpu

    endsubroutine force_var_1_c2_kernel_wrapper
  endinterface

  interface
    subroutine force_var_2_c2_kernel_wrapper(stream,nx,ny,nz,ng,indx_cp_l,indx_cp_r,&
    &calorically_perfect,tbdiff,t0,fluid_mask_gpu,w_gpu,w_aux_gpu,cv_coeff_gpu,jac_gpu)bind(c,&
    &name="force_var_2_c2_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect
      real(c_rkind), value :: tbdiff,t0
      type(c_ptr), value :: fluid_mask_gpu
      type(c_ptr), value :: w_gpu,w_aux_gpu,cv_coeff_gpu,jac_gpu

    endsubroutine force_var_2_c2_kernel_wrapper
  endinterface

  interface
    subroutine update_flux_kernel_wrapper(stream,nx,ny,nz,nv,gamdt,fl_gpu,fln_gpu)bind(c,name="update_flux_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,nv
      real(c_rkind), value :: gamdt
      type(c_ptr), value :: fl_gpu,fln_gpu

    endsubroutine update_flux_kernel_wrapper
  endinterface

  interface
    subroutine update_field_kernel_wrapper(stream,nx,ny,nz,nv,ng,fluid_mask_gpu,w_gpu,&
    &fln_gpu)bind(c,name="update_field_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,nv,ng
      type(c_ptr), value :: fluid_mask_gpu
      type(c_ptr), value :: w_gpu,fln_gpu

    endsubroutine update_field_kernel_wrapper
  endinterface

  interface
    subroutine visflx_div_ord2_kernel_wrapper(stream,nx,ny,nz,ng,w_aux_gpu,fl_gpu,x_gpu,y_gpu,&
    &z_gpu)bind(c,name="visflx_div_ord2_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng
      type(c_ptr), value :: w_aux_gpu,fl_gpu,x_gpu,y_gpu,z_gpu

    endsubroutine visflx_div_ord2_kernel_wrapper
  endinterface

  interface
    subroutine visflx_div_kernel_wrapper(stream,nx,ny,nz,ng,visc_order,lmax,w_aux_gpu,fl_gpu,&
    &coeff_deriv1_gpu,dcsidx_gpu,detady_gpu,dzitdz_gpu)bind(c,name="visflx_div_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,visc_order,lmax
      type(c_ptr), value :: w_aux_gpu,fl_gpu,coeff_deriv1_gpu,dcsidx_gpu,detady_gpu,dzitdz_gpu

    endsubroutine visflx_div_kernel_wrapper
  endinterface

  interface
    subroutine visflx_div_c2_kernel_wrapper(stream,nx,ny,nz,ng,visc_order,lmax,vis_tag_gpu,&
    &wall_tag_gpu,w_aux_gpu,fl_gpu,coeff_deriv1_gpu,dcsidxc2_gpu,detadyc2_gpu,dcsidyc2_gpu,detadxc2_gpu,&
    &dzitdz_gpu)bind(c,name="visflx_div_c2_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,visc_order,lmax
      type(c_ptr), value :: vis_tag_gpu,wall_tag_gpu
      type(c_ptr), value :: w_aux_gpu,fl_gpu,coeff_deriv1_gpu,dcsidxc2_gpu,detadyc2_gpu,dcsidyc2_gpu,detadxc2_gpu,dzitdz_gpu

    endsubroutine visflx_div_c2_kernel_wrapper
  endinterface

  interface
    subroutine visflx_kernel_wrapper(stream,nx,ny,nz,ng,visc_order,calorically_perfect,indx_cp_l,&
    &indx_cp_r,lmax,prandtl,u0,l0,t0,w_gpu,w_aux_gpu,wallprop_gpu,fl_gpu,coeff_deriv1_gpu,&
    &coeff_deriv2_gpu,dcsidx_gpu,detady_gpu,dzitdz_gpu,dcsidxs_gpu,detadys_gpu,dzitdzs_gpu,dcsidx2_gpu,&
    &detady2_gpu,dzitdz2_gpu,cp_coeff_gpu)bind(c,name="visflx_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,visc_order,calorically_perfect,indx_cp_l,indx_cp_r,lmax
      real(c_rkind), value :: prandtl,u0,l0,t0
      type(c_ptr), value :: w_gpu,w_aux_gpu,wallprop_gpu,fl_gpu,coeff_deriv1_gpu,coeff_deriv2_gpu,&
      &dcsidx_gpu,detady_gpu,dzitdz_gpu,dcsidxs_gpu,detadys_gpu,dzitdzs_gpu,dcsidx2_gpu,detady2_gpu,&
      &dzitdz2_gpu,cp_coeff_gpu

    endsubroutine visflx_kernel_wrapper
  endinterface

  interface
    subroutine visflx_c2_kernel_wrapper(stream,nx,ny,nz,ng,visc_order,calorically_perfect,indx_cp_l,&
    &indx_cp_r,iblock,ite_rank_x,itu_rank_x,ite_l,itu_l,lmax,prandtl,u0,l0,t0,teshk,vis_tag_gpu,&
    &wall_tag_gpu,w_gpu,w_aux_gpu,wallprop_gpu,fl_gpu,coeff_deriv1_gpu,coeff_deriv2_gpu,dzitdz_gpu,&
    &dzitdzs_gpu,dzitdz2_gpu,dcsidxc2_gpu,detadyc2_gpu,detadxc2_gpu,dcsidyc2_gpu,g1_gpu,g2_gpu,g12_gpu,&
    &jac_gpu,cp_coeff_gpu)bind(c,name="visflx_c2_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,visc_order,calorically_perfect,indx_cp_l,indx_cp_r,&
      &iblock,ite_rank_x,itu_rank_x,ite_l,itu_l,lmax
      real(c_rkind), value :: prandtl,u0,l0,t0,teshk
      type(c_ptr), value :: vis_tag_gpu,wall_tag_gpu
      type(c_ptr), value :: w_gpu,w_aux_gpu,wallprop_gpu,fl_gpu,coeff_deriv1_gpu,coeff_deriv2_gpu,&
      &dzitdz_gpu,dzitdzs_gpu,dzitdz2_gpu,dcsidxc2_gpu,detadyc2_gpu,detadxc2_gpu,dcsidyc2_gpu,g1_gpu,&
      &g2_gpu,g12_gpu,jac_gpu,cp_coeff_gpu

    endsubroutine visflx_c2_kernel_wrapper
  endinterface

  interface
    subroutine visflx_nosensor_kernel_wrapper(stream,nx,ny,nz,ng,visc_order,calorically_perfect,&
    &indx_cp_l,indx_cp_r,lmax,prandtl,u0,l0,t0,w_gpu,w_aux_gpu,wallprop_gpu,fl_gpu,coeff_deriv1_gpu,&
    &coeff_deriv2_gpu,dcsidx_gpu,detady_gpu,dzitdz_gpu,dcsidxs_gpu,detadys_gpu,dzitdzs_gpu,dcsidx2_gpu,&
    &detady2_gpu,dzitdz2_gpu,cp_coeff_gpu)bind(c,name="visflx_nosensor_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,visc_order,calorically_perfect,indx_cp_l,indx_cp_r,lmax
      real(c_rkind), value :: prandtl,u0,l0,t0
      type(c_ptr), value :: w_gpu,w_aux_gpu,wallprop_gpu,fl_gpu,coeff_deriv1_gpu,coeff_deriv2_gpu,&
      &dcsidx_gpu,detady_gpu,dzitdz_gpu,dcsidxs_gpu,detadys_gpu,dzitdzs_gpu,dcsidx2_gpu,detady2_gpu,&
      &dzitdz2_gpu,cp_coeff_gpu

    endsubroutine visflx_nosensor_kernel_wrapper
  endinterface

  interface
    subroutine visflx_nosensor_c2_kernel_wrapper(stream,nx,ny,nz,ng,visc_order,calorically_perfect,&
    &indx_cp_l,indx_cp_r,ortho,iblock,ite_rank_x,itu_rank_x,ite_l,itu_l,lmax,prandtl,u0,l0,t0,&
    &vis_tag_gpu,wall_tag_gpu,w_gpu,w_aux_gpu,wallprop_gpu,fl_gpu,coeff_deriv1_gpu,coeff_deriv2_gpu,&
    &dzitdz_gpu,dzitdzs_gpu,dzitdz2_gpu,dcsidxc2_gpu,detadyc2_gpu,detadxc2_gpu,dcsidyc2_gpu,g1_gpu,&
    &g2_gpu,g12_gpu,jac_gpu,cp_coeff_gpu)bind(c,name="visflx_nosensor_c2_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,visc_order,calorically_perfect,indx_cp_l,indx_cp_r,ortho,&
      &iblock,ite_rank_x,itu_rank_x,ite_l,itu_l,lmax
      real(c_rkind), value :: prandtl,u0,l0,t0
      type(c_ptr), value :: vis_tag_gpu,wall_tag_gpu
      type(c_ptr), value :: w_gpu,w_aux_gpu,wallprop_gpu,fl_gpu,coeff_deriv1_gpu,coeff_deriv2_gpu,&
      &dzitdz_gpu,dzitdzs_gpu,dzitdz2_gpu,dcsidxc2_gpu,detadyc2_gpu,detadxc2_gpu,dcsidyc2_gpu,g1_gpu,&
      &g2_gpu,g12_gpu,jac_gpu,cp_coeff_gpu

    endsubroutine visflx_nosensor_c2_kernel_wrapper
  endinterface

  interface
    subroutine sponge_kernel_wrapper(stream,nx,ny,nz,ng,nv,j_sponge,w_gpu,fln_gpu,wfar_gpu,&
    &f_sponge_gpu)bind(c,name="sponge_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,j_sponge
      type(c_ptr), value :: w_gpu,fln_gpu,wfar_gpu,f_sponge_gpu

    endsubroutine sponge_kernel_wrapper
  endinterface

  interface
    subroutine limiter_kernel1_wrapper(stream,nx,ny,nz,ng,iblock,kblock,calorically_perfect,&
    &indx_cp_l,indx_cp_r,t0,tol_iter_nr,rho_lim,tem_lim,rho_lim_rescale,tem_lim_rescale,w_gpu,w_aux_gpu,&
    &cv_coeff_gpu,n_limited_rho,n_limited_tem,redn_3d_gpu)bind(c,name="limiter_kernel1_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,iblock,kblock,calorically_perfect,indx_cp_l,indx_cp_r
      real(c_rkind), value :: t0,tol_iter_nr,rho_lim,tem_lim,rho_lim_rescale,tem_lim_rescale
      type(c_ptr), value :: w_gpu,w_aux_gpu,cv_coeff_gpu
      real(c_rkind) :: n_limited_rho,n_limited_tem
      type(c_ptr), value :: redn_3d_gpu

    endsubroutine limiter_kernel1_wrapper
  endinterface

  interface
    subroutine limiter_kernel2_wrapper(stream,nx,ny,nz,ng,iblock,kblock,calorically_perfect,&
    &indx_cp_l,indx_cp_r,t0,tol_iter_nr,rho_lim,tem_lim,rho_lim_rescale,tem_lim_rescale,w_gpu,w_aux_gpu,&
    &cv_coeff_gpu)bind(c,name="limiter_kernel2_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,iblock,kblock,calorically_perfect,indx_cp_l,indx_cp_r
      real(c_rkind), value :: t0,tol_iter_nr,rho_lim,tem_lim,rho_lim_rescale,tem_lim_rescale
      type(c_ptr), value :: w_gpu,w_aux_gpu,cv_coeff_gpu

    endsubroutine limiter_kernel2_wrapper
  endinterface

  interface
    subroutine filter_kernel1_wrapper(stream,nx,ny,nz,ng,calorically_perfect,indx_cp_l,indx_cp_r,&
    &jfilter,t0,tol_iter_nr,wall_tag_gpu,w_gpu,w_aux_gpu,cv_coeff_gpu,coeff_filter_gpu)bind(c,&
    &name="filter_kernel1_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,calorically_perfect,indx_cp_l,indx_cp_r,jfilter
      real(c_rkind), value :: t0,tol_iter_nr
      type(c_ptr), value :: wall_tag_gpu
      type(c_ptr), value :: w_gpu,w_aux_gpu,cv_coeff_gpu,coeff_filter_gpu

    endsubroutine filter_kernel1_wrapper
  endinterface

  interface
    subroutine filter_kernel2_wrapper(stream,nx,ny,nz,ng,calorically_perfect,indx_cp_l,indx_cp_r,&
    &jfilter,t0,tol_iter_nr,wall_tag_gpu,w_gpu,w_aux_gpu,cv_coeff_gpu,coeff_filter_gpu)bind(c,&
    &name="filter_kernel2_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,calorically_perfect,indx_cp_l,indx_cp_r,jfilter
      real(c_rkind), value :: t0,tol_iter_nr
      type(c_ptr), value :: wall_tag_gpu
      type(c_ptr), value :: w_gpu,w_aux_gpu,cv_coeff_gpu,coeff_filter_gpu

    endsubroutine filter_kernel2_wrapper
  endinterface

  interface
    subroutine visflx_reduced_ord2_kernel_wrapper(stream,nx,ny,nz,ng,calorically_perfect,indx_cp_l,&
    &indx_cp_r,update_sensor,prandtl,u0,l0,t0,w_gpu,w_aux_gpu,wallprop_gpu,fl_gpu,x_gpu,y_gpu,z_gpu,&
    &cp_coeff_gpu)bind(c,name="visflx_reduced_ord2_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,calorically_perfect,indx_cp_l,indx_cp_r,update_sensor
      real(c_rkind), value :: prandtl,u0,l0,t0
      type(c_ptr), value :: w_gpu,w_aux_gpu,wallprop_gpu,fl_gpu,x_gpu,y_gpu,z_gpu,cp_coeff_gpu

    endsubroutine visflx_reduced_ord2_kernel_wrapper
  endinterface

  interface
    subroutine sensor_kernel_wrapper(stream,nx,ny,nz,ng,u0,l0,w_aux_gpu)bind(c,name="sensor_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng
      real(c_rkind), value :: u0,l0
      type(c_ptr), value :: w_aux_gpu

    endsubroutine sensor_kernel_wrapper
  endinterface

  interface
    subroutine sensor_c2_kernel_wrapper(stream,nx,ny,nz,ng,iblock,ite_rank_x,itu_rank_x,ite_l,itu_l,&
    &jweno,u0,l0_ducros,teshk,theta_threshold,w_aux_gpu,theta_ij_gpu)bind(c,name="sensor_c2_kernel_wrappe&
    &r")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,iblock,ite_rank_x,itu_rank_x,ite_l,itu_l,jweno
      real(c_rkind), value :: u0,l0_ducros,teshk,theta_threshold
      type(c_ptr), value :: w_aux_gpu,theta_ij_gpu

    endsubroutine sensor_c2_kernel_wrapper
  endinterface

  interface
    subroutine visflx_x_kernel1_wrapper(stream,nx,ny,nz,nv,ng,calorically_perfect,indx_cp_l,&
    &indx_cp_r,prandtl,t0,fl_gpu,w_aux_gpu,fhat_gpu,x_gpu,cp_coeff_gpu)bind(c,name="visflx_x_kernel1_wrap&
    &per")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,nv,ng,calorically_perfect,indx_cp_l,indx_cp_r
      real(c_rkind), value :: prandtl,t0
      type(c_ptr), value :: fl_gpu,w_aux_gpu,fhat_gpu,x_gpu,cp_coeff_gpu

    endsubroutine visflx_x_kernel1_wrapper
  endinterface

  interface
    subroutine visflx_x_kernel2_wrapper(stream,nx,ny,nz,nv,ng,calorically_perfect,indx_cp_l,&
    &indx_cp_r,prandtl,t0,fl_gpu,w_aux_gpu,fhat_gpu,x_gpu,cp_coeff_gpu)bind(c,name="visflx_x_kernel2_wrap&
    &per")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,nv,ng,calorically_perfect,indx_cp_l,indx_cp_r
      real(c_rkind), value :: prandtl,t0
      type(c_ptr), value :: fl_gpu,w_aux_gpu,fhat_gpu,x_gpu,cp_coeff_gpu

    endsubroutine visflx_x_kernel2_wrapper
  endinterface

  interface
    subroutine visflx_y_kernel_wrapper(stream,nx,ny,nz,nv,ng,calorically_perfect,indx_cp_l,&
    &indx_cp_r,prandtl,t0,fl_gpu,w_aux_gpu,y_gpu,cp_coeff_gpu)bind(c,name="visflx_y_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,nv,ng,calorically_perfect,indx_cp_l,indx_cp_r
      real(c_rkind), value :: prandtl,t0
      type(c_ptr), value :: fl_gpu,w_aux_gpu,y_gpu,cp_coeff_gpu

    endsubroutine visflx_y_kernel_wrapper
  endinterface

  interface
    subroutine visflx_z_kernel_wrapper(stream,nx,ny,nz,nv,ng,calorically_perfect,indx_cp_l,&
    &indx_cp_r,prandtl,t0,fl_gpu,w_aux_gpu,z_gpu,cp_coeff_gpu)bind(c,name="visflx_z_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,nv,ng,calorically_perfect,indx_cp_l,indx_cp_r
      real(c_rkind), value :: prandtl,t0
      type(c_ptr), value :: fl_gpu,w_aux_gpu,z_gpu,cp_coeff_gpu

    endsubroutine visflx_z_kernel_wrapper
  endinterface

  interface
    subroutine recyc_exchange_kernel_1_wrapper(stream,irecyc,nx,ny,nz,ng,nv,w_gpu,wbuf1s_gpu)bind(c,&
    &name="recyc_exchange_kernel_1_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: irecyc,nx,ny,nz,ng,nv
      type(c_ptr), value :: w_gpu,wbuf1s_gpu

    endsubroutine recyc_exchange_kernel_1_wrapper
  endinterface

  interface
    subroutine recyc_exchange_kernel_2_wrapper(stream,n1_start_recv,n1_start_send,n1_end_recv,nx,ny,&
    &nz,ng,nv,wbuf1r_gpu,wrecyc_gpu)bind(c,name="recyc_exchange_kernel_2_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: n1_start_recv,n1_start_send,n1_end_recv,nx,ny,nz,ng,nv
      type(c_ptr), value :: wbuf1r_gpu,wrecyc_gpu

    endsubroutine recyc_exchange_kernel_2_wrapper
  endinterface

  interface
    subroutine recyc_exchange_kernel_3_wrapper(stream,n2_start_recv,n2_start_send,n2_end_recv,nx,ny,&
    &nz,ng,nv,wbuf2r_gpu,wrecyc_gpu)bind(c,name="recyc_exchange_kernel_3_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: n2_start_recv,n2_start_send,n2_end_recv,nx,ny,nz,ng,nv
      type(c_ptr), value :: wbuf2r_gpu,wrecyc_gpu

    endsubroutine recyc_exchange_kernel_3_wrapper
  endinterface

  interface
    subroutine bcextr_sub_kernel1_wrapper(stream,ilat,nx,ny,nz,ng,indx_cp_l,indx_cp_r,&
    &calorically_perfect,p0,t0,rgas0,cv_coeff_gpu,w_gpu)bind(c,name="bcextr_sub_kernel1_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: ilat,nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect
      real(c_rkind), value :: p0,t0,rgas0
      type(c_ptr), value :: cv_coeff_gpu,w_gpu

    endsubroutine bcextr_sub_kernel1_wrapper
  endinterface

  interface
    subroutine bcextr_sub_kernel2_wrapper(stream,ilat,nx,ny,nz,ng,indx_cp_l,indx_cp_r,&
    &calorically_perfect,p0,t0,rgas0,cv_coeff_gpu,w_gpu)bind(c,name="bcextr_sub_kernel2_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: ilat,nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect
      real(c_rkind), value :: p0,t0,rgas0
      type(c_ptr), value :: cv_coeff_gpu,w_gpu

    endsubroutine bcextr_sub_kernel2_wrapper
  endinterface

  interface
    subroutine bcextr_sub_kernel3_wrapper(stream,ilat,nx,ny,nz,ng,indx_cp_l,indx_cp_r,&
    &calorically_perfect,p0,t0,rgas0,cv_coeff_gpu,w_gpu)bind(c,name="bcextr_sub_kernel3_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: ilat,nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect
      real(c_rkind), value :: p0,t0,rgas0
      type(c_ptr), value :: cv_coeff_gpu,w_gpu

    endsubroutine bcextr_sub_kernel3_wrapper
  endinterface

  interface
    subroutine bcextr_sub_kernel4_wrapper(stream,ilat,nx,ny,nz,ng,indx_cp_l,indx_cp_r,&
    &calorically_perfect,p0,t0,rgas0,cv_coeff_gpu,w_gpu)bind(c,name="bcextr_sub_kernel4_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: ilat,nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect
      real(c_rkind), value :: p0,t0,rgas0
      type(c_ptr), value :: cv_coeff_gpu,w_gpu

    endsubroutine bcextr_sub_kernel4_wrapper
  endinterface

  interface
    subroutine bcextr_sub_kernel5_wrapper(stream,ilat,nx,ny,nz,ng,indx_cp_l,indx_cp_r,&
    &calorically_perfect,p0,t0,rgas0,cv_coeff_gpu,w_gpu)bind(c,name="bcextr_sub_kernel5_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: ilat,nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect
      real(c_rkind), value :: p0,t0,rgas0
      type(c_ptr), value :: cv_coeff_gpu,w_gpu

    endsubroutine bcextr_sub_kernel5_wrapper
  endinterface

  interface
    subroutine bcrecyc_kernel_1_wrapper(stream,nx,ny,nz,ng,nv,wrecycav_gpu,wrecyc_gpu)bind(c,name="bcrecyc_kernel_1_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv
      type(c_ptr), value :: wrecycav_gpu,wrecyc_gpu

    endsubroutine bcrecyc_kernel_1_wrapper
  endinterface

  interface
    subroutine bcrecyc_kernel_2_wrapper(stream,nx,ny,nz,nzmax,ng,wrecycav_gpu,wrecyc_gpu)bind(c,name="bcrecyc_kernel_2_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,nzmax,ng
      type(c_ptr), value :: wrecycav_gpu,wrecyc_gpu

    endsubroutine bcrecyc_kernel_2_wrapper
  endinterface

  interface
    subroutine bcrecyc_kernel_3_wrapper(stream,nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect,&
    &rand_type,p0,rgas0,betarecyc,glund1,t0,u0,u0_02,map_j_inn_gpu,map_j_out_gpu,map_j_out_blend_gpu,&
    &cv_coeff_gpu,cp_coeff_gpu,wmean_gpu,wrecyc_gpu,inflow_random_plane_gpu,w_gpu,weta_inflow_gpu,&
    &yplus_inflow_gpu,eta_inflow_gpu,yplus_recyc_gpu,eta_recyc_gpu,eta_recyc_blend_gpu)bind(c,&
    &name="bcrecyc_kernel_3_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect,rand_type
      real(c_rkind), value :: p0,rgas0,betarecyc,glund1,t0,u0,u0_02
      type(c_ptr), value :: map_j_inn_gpu,map_j_out_gpu,map_j_out_blend_gpu
      type(c_ptr), value :: cv_coeff_gpu,cp_coeff_gpu,wmean_gpu,wrecyc_gpu,inflow_random_plane_gpu,&
      &w_gpu,weta_inflow_gpu,yplus_inflow_gpu,eta_inflow_gpu,yplus_recyc_gpu,eta_recyc_gpu,&
      &eta_recyc_blend_gpu

    endsubroutine bcrecyc_kernel_3_wrapper
  endinterface

  interface
    subroutine bclam_kernel_wrapper(stream,nx,ny,nz,ng,nv,indx_cp_l,indx_cp_r,calorically_perfect,&
    &ilat,p0,rgas0,t0,cv_coeff_gpu,w_gpu,wmean_gpu)bind(c,name="bclam_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,indx_cp_l,indx_cp_r,calorically_perfect,ilat
      real(c_rkind), value :: p0,rgas0,t0
      type(c_ptr), value :: cv_coeff_gpu,w_gpu,wmean_gpu

    endsubroutine bclam_kernel_wrapper
  endinterface

  interface
    subroutine bcfree_kernel1_wrapper(stream,ilat,nx,ny,nz,ng,nv,winf_gpu,w_gpu)bind(c,name="bcfree_kernel1_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: ilat,nx,ny,nz,ng,nv
      type(c_ptr), value :: winf_gpu,w_gpu

    endsubroutine bcfree_kernel1_wrapper
  endinterface

  interface
    subroutine bcfree_kernel2_wrapper(stream,ilat,nx,ny,nz,ng,nv,winf_gpu,w_gpu)bind(c,name="bcfree_kernel2_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: ilat,nx,ny,nz,ng,nv
      type(c_ptr), value :: winf_gpu,w_gpu

    endsubroutine bcfree_kernel2_wrapper
  endinterface

  interface
    subroutine bcfree_kernel3_wrapper(stream,ilat,nx,ny,nz,ng,nv,winf_gpu,w_gpu)bind(c,name="bcfree_kernel3_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: ilat,nx,ny,nz,ng,nv
      type(c_ptr), value :: winf_gpu,w_gpu

    endsubroutine bcfree_kernel3_wrapper
  endinterface

  interface
    subroutine bcfree_kernel4_wrapper(stream,ilat,nx,ny,nz,ng,nv,winf_gpu,w_gpu)bind(c,name="bcfree_kernel4_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: ilat,nx,ny,nz,ng,nv
      type(c_ptr), value :: winf_gpu,w_gpu

    endsubroutine bcfree_kernel4_wrapper
  endinterface

  interface
    subroutine bcfree_kernel5_wrapper(stream,ilat,nx,ny,nz,ng,nv,winf_gpu,w_gpu)bind(c,name="bcfree_kernel5_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: ilat,nx,ny,nz,ng,nv
      type(c_ptr), value :: winf_gpu,w_gpu

    endsubroutine bcfree_kernel5_wrapper
  endinterface

  interface
    subroutine bcfree_kernel6_wrapper(stream,ilat,nx,ny,nz,ng,nv,winf_gpu,w_gpu)bind(c,name="bcfree_kernel6_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: ilat,nx,ny,nz,ng,nv
      type(c_ptr), value :: winf_gpu,w_gpu

    endsubroutine bcfree_kernel6_wrapper
  endinterface

  interface
    subroutine bcfree_sub_kernel1_wrapper(stream,ilat,nx,ny,nz,ng,indx_cp_l,indx_cp_r,&
    &calorically_perfect,ptot0,ttot0,rgas0,aoa,t0,tol_iter_nr,cosangle,sinangle,w_gpu,w_aux_gpu,&
    &cv_coeff_gpu,cp_coeff_gpu)bind(c,name="bcfree_sub_kernel1_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: ilat,nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect
      real(c_rkind), value :: ptot0,ttot0,rgas0,aoa,t0,tol_iter_nr,cosangle,sinangle
      type(c_ptr), value :: w_gpu,w_aux_gpu,cv_coeff_gpu,cp_coeff_gpu

    endsubroutine bcfree_sub_kernel1_wrapper
  endinterface

  interface
    subroutine bcfree_sub_kernel2_wrapper(stream,ilat,nx,ny,nz,ng,indx_cp_l,indx_cp_r,&
    &calorically_perfect,ptot0,ttot0,rgas0,aoa,t0,tol_iter_nr,cosangle,sinangle,w_gpu,w_aux_gpu,&
    &cv_coeff_gpu,cp_coeff_gpu)bind(c,name="bcfree_sub_kernel2_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: ilat,nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect
      real(c_rkind), value :: ptot0,ttot0,rgas0,aoa,t0,tol_iter_nr,cosangle,sinangle
      type(c_ptr), value :: w_gpu,w_aux_gpu,cv_coeff_gpu,cp_coeff_gpu

    endsubroutine bcfree_sub_kernel2_wrapper
  endinterface

  interface
    subroutine bcshock_kernel_wrapper(stream,nx,ny,nz,ng,nv,ilat,xshock_imp,shock_angle,tanhfacs,&
    &tanhlen,winf_gpu,winf_past_shock_gpu,w_gpu,x_gpu,y_gpu)bind(c,name="bcshock_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,ilat
      real(c_rkind), value :: xshock_imp,shock_angle,tanhfacs,tanhlen
      type(c_ptr), value :: winf_gpu,winf_past_shock_gpu,w_gpu,x_gpu,y_gpu

    endsubroutine bcshock_kernel_wrapper
  endinterface

  interface
    subroutine bcextr_var_kernel1_wrapper(stream,nx,ny,nz,ng,w_gpu)bind(c,name="bcextr_var_kernel1_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng
      type(c_ptr), value :: w_gpu

    endsubroutine bcextr_var_kernel1_wrapper
  endinterface

  interface
    subroutine bcextr_var_kernel2_wrapper(stream,nx,ny,nz,ng,w_gpu)bind(c,name="bcextr_var_kernel2_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng
      type(c_ptr), value :: w_gpu

    endsubroutine bcextr_var_kernel2_wrapper
  endinterface

  interface
    subroutine bcextr_var_kernel3_wrapper(stream,nx,ny,nz,ng,w_gpu)bind(c,name="bcextr_var_kernel3_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng
      type(c_ptr), value :: w_gpu

    endsubroutine bcextr_var_kernel3_wrapper
  endinterface

  interface
    subroutine bcextr_var_kernel4_wrapper(stream,nx,ny,nz,ng,w_gpu)bind(c,name="bcextr_var_kernel4_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng
      type(c_ptr), value :: w_gpu

    endsubroutine bcextr_var_kernel4_wrapper
  endinterface

  interface
    subroutine bcextr_var_kernel5_wrapper(stream,nx,ny,nz,ng,w_gpu)bind(c,name="bcextr_var_kernel5_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng
      type(c_ptr), value :: w_gpu

    endsubroutine bcextr_var_kernel5_wrapper
  endinterface

  interface
    subroutine bcextr_var_kernel6_wrapper(stream,nx,ny,nz,ng,w_gpu)bind(c,name="bcextr_var_kernel6_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng
      type(c_ptr), value :: w_gpu

    endsubroutine bcextr_var_kernel6_wrapper
  endinterface

  interface
    subroutine bcextr_airfoil_var_kernel1_wrapper(stream,nx,ny,nz,ng,ndim,ileftx,irightx,ileftz,&
    &irightz,wall_tag_gpu,w_gpu)bind(c,name="bcextr_airfoil_var_kernel1_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,ndim,ileftx,irightx,ileftz,irightz
      type(c_ptr), value :: wall_tag_gpu
      type(c_ptr), value :: w_gpu

    endsubroutine bcextr_airfoil_var_kernel1_wrapper
  endinterface

  interface
    subroutine bcextr_airfoil_var_kernel2_wrapper(stream,nx,ny,nz,ng,ndim,ileftx,irightx,ileftz,&
    &irightz,wall_tag_gpu,w_gpu)bind(c,name="bcextr_airfoil_var_kernel2_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,ndim,ileftx,irightx,ileftz,irightz
      type(c_ptr), value :: wall_tag_gpu
      type(c_ptr), value :: w_gpu

    endsubroutine bcextr_airfoil_var_kernel2_wrapper
  endinterface

  interface
    subroutine bcextr_airfoil_var_kernel3_wrapper(stream,nx,ny,nz,ng,ndim,ileftx,irightx,ileftz,&
    &irightz,wall_tag_gpu,w_gpu)bind(c,name="bcextr_airfoil_var_kernel3_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,ndim,ileftx,irightx,ileftz,irightz
      type(c_ptr), value :: wall_tag_gpu
      type(c_ptr), value :: w_gpu

    endsubroutine bcextr_airfoil_var_kernel3_wrapper
  endinterface

  interface
    subroutine bcextr_airfoil_var_kernel4_wrapper(stream,nx,ny,nz,ng,ndim,ileftx,irightx,ileftz,&
    &irightz,wall_tag_gpu,w_gpu)bind(c,name="bcextr_airfoil_var_kernel4_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,ndim,ileftx,irightx,ileftz,irightz
      type(c_ptr), value :: wall_tag_gpu
      type(c_ptr), value :: w_gpu

    endsubroutine bcextr_airfoil_var_kernel4_wrapper
  endinterface

  interface
    subroutine bcextr_airfoil_var_kernel5_wrapper(stream,nx,ny,nz,ng,ndim,ileftx,irightx,ileftz,&
    &irightz,wall_tag_gpu,w_gpu)bind(c,name="bcextr_airfoil_var_kernel5_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,ndim,ileftx,irightx,ileftz,irightz
      type(c_ptr), value :: wall_tag_gpu
      type(c_ptr), value :: w_gpu

    endsubroutine bcextr_airfoil_var_kernel5_wrapper
  endinterface

  interface
    subroutine bcextr_airfoil_var_kernel6_wrapper(stream,nx,ny,nz,ng,ndim,ileftx,irightx,ileftz,&
    &irightz,wall_tag_gpu,w_gpu)bind(c,name="bcextr_airfoil_var_kernel6_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,ndim,ileftx,irightx,ileftz,irightz
      type(c_ptr), value :: wall_tag_gpu
      type(c_ptr), value :: w_gpu

    endsubroutine bcextr_airfoil_var_kernel6_wrapper
  endinterface

  interface
    subroutine bcextr_kernel1_wrapper(stream,nx,ny,nz,ng,nv,ilat,w_gpu)bind(c,name="bcextr_kernel1_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,ilat
      type(c_ptr), value :: w_gpu

    endsubroutine bcextr_kernel1_wrapper
  endinterface

  interface
    subroutine bcextr_kernel2_wrapper(stream,nx,ny,nz,ng,nv,ilat,w_gpu)bind(c,name="bcextr_kernel2_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,ilat
      type(c_ptr), value :: w_gpu

    endsubroutine bcextr_kernel2_wrapper
  endinterface

  interface
    subroutine bcextr_kernel3_wrapper(stream,nx,ny,nz,ng,nv,ilat,w_gpu)bind(c,name="bcextr_kernel3_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,ilat
      type(c_ptr), value :: w_gpu

    endsubroutine bcextr_kernel3_wrapper
  endinterface

  interface
    subroutine bcextr_kernel4_wrapper(stream,nx,ny,nz,ng,nv,ilat,w_gpu)bind(c,name="bcextr_kernel4_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,ilat
      type(c_ptr), value :: w_gpu

    endsubroutine bcextr_kernel4_wrapper
  endinterface

  interface
    subroutine bcextr_kernel5_wrapper(stream,nx,ny,nz,ng,nv,ilat,w_gpu)bind(c,name="bcextr_kernel5_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,ilat
      type(c_ptr), value :: w_gpu

    endsubroutine bcextr_kernel5_wrapper
  endinterface

  interface
    subroutine bcextr_kernel6_wrapper(stream,nx,ny,nz,ng,nv,ilat,w_gpu)bind(c,name="bcextr_kernel6_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,ilat
      type(c_ptr), value :: w_gpu

    endsubroutine bcextr_kernel6_wrapper
  endinterface

  interface
    subroutine bcsym_kernel_wrapper(stream,nx,ny,nz,ng,ilat,w_gpu)bind(c,name="bcsym_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,ilat
      type(c_ptr), value :: w_gpu

    endsubroutine bcsym_kernel_wrapper
  endinterface

  interface
    subroutine bcsym_c2_kernel1_wrapper(stream,nx,ny,nz,ng,ilat,w_gpu,dxdcsic2_gpu,&
    &dydcsic2_gpu)bind(c,name="bcsym_c2_kernel1_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,ilat
      type(c_ptr), value :: w_gpu,dxdcsic2_gpu,dydcsic2_gpu

    endsubroutine bcsym_c2_kernel1_wrapper
  endinterface

  interface
    subroutine bcsym_c2_kernel2_wrapper(stream,nx,ny,nz,ng,ilat,w_gpu,dxdcsic2_gpu,&
    &dydcsic2_gpu)bind(c,name="bcsym_c2_kernel2_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,ilat
      type(c_ptr), value :: w_gpu,dxdcsic2_gpu,dydcsic2_gpu

    endsubroutine bcsym_c2_kernel2_wrapper
  endinterface

  interface
    subroutine bcwall_kernel_wrapper(stream,nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect,&
    &ilat,twall,t0,rgas0,tol_iter_nr,w_gpu,w_aux_gpu,cv_coeff_gpu)bind(c,name="bcwall_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect,ilat
      real(c_rkind), value :: twall,t0,rgas0,tol_iter_nr
      type(c_ptr), value :: w_gpu,w_aux_gpu,cv_coeff_gpu

    endsubroutine bcwall_kernel_wrapper
  endinterface

  interface
    subroutine bcwall_airfoil_kernel_wrapper(stream,nx,ny,nz,ng,indx_cp_l,indx_cp_r,&
    &calorically_perfect,ilat,u0,twall,t0,rgas0,tol_iter_nr,a_tw,v_bs,thic,kx_tw,om_tw,time,xtw1,xtw2,&
    &wall_tag_gpu,w_gpu,w_aux_gpu,cv_coeff_gpu,xc2_gpu,yc2_gpu,dxdetanc2_gpu,dydetanc2_gpu)bind(c,&
    &name="bcwall_airfoil_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect,ilat
      real(c_rkind), value :: u0,twall,t0,rgas0,tol_iter_nr,a_tw,v_bs,thic,kx_tw,om_tw,time,xtw1,xtw2
      type(c_ptr), value :: wall_tag_gpu
      type(c_ptr), value :: w_gpu,w_aux_gpu,cv_coeff_gpu,xc2_gpu,yc2_gpu,dxdetanc2_gpu,dydetanc2_gpu

    endsubroutine bcwall_airfoil_kernel_wrapper
  endinterface

  interface
    subroutine bcwall_staggered_kernel1_wrapper(stream,nx,ny,nz,ng,indx_cp_l,indx_cp_r,&
    &calorically_perfect,ilat,twall,t0,rgas0,tol_iter_nr,w_gpu,w_aux_gpu,cv_coeff_gpu)bind(c,&
    &name="bcwall_staggered_kernel1_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect,ilat
      real(c_rkind), value :: twall,t0,rgas0,tol_iter_nr
      type(c_ptr), value :: w_gpu,w_aux_gpu,cv_coeff_gpu

    endsubroutine bcwall_staggered_kernel1_wrapper
  endinterface

  interface
    subroutine bcwall_staggered_kernel2_wrapper(stream,nx,ny,nz,ng,indx_cp_l,indx_cp_r,&
    &calorically_perfect,ilat,twall,t0,rgas0,tol_iter_nr,w_gpu,w_aux_gpu,cv_coeff_gpu)bind(c,&
    &name="bcwall_staggered_kernel2_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect,ilat
      real(c_rkind), value :: twall,t0,rgas0,tol_iter_nr
      type(c_ptr), value :: w_gpu,w_aux_gpu,cv_coeff_gpu

    endsubroutine bcwall_staggered_kernel2_wrapper
  endinterface

  interface
    subroutine compute_residual_kernel_wrapper(stream,nx,ny,nz,ng,nv,dt,fluid_mask_gpu,fln_gpu,&
    &residual_rhou,redn_3d_gpu)bind(c,name="compute_residual_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv
      real(c_rkind), value :: dt
      type(c_ptr), value :: fluid_mask_gpu
      type(c_ptr), value :: fln_gpu
      real(c_rkind) :: residual_rhou
      type(c_ptr), value :: redn_3d_gpu

    endsubroutine compute_residual_kernel_wrapper
  endinterface

  interface
    subroutine compute_airfoil_forces_runtime_kernel_wrapper(stream,nx,ny,nz,ng,nv,p0,u0,rgas0,&
    &wall_tag_gpu,w_aux_gpu,meta_gpu,csimod_gpu,dxdcsic2_gpu,dydcsic2_gpu,n,a,pn,pa,tn,ta,&
    &redn_3d_gpu)bind(c,name="compute_airfoil_forces_runtime_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv
      real(c_rkind), value :: p0,u0,rgas0
      type(c_ptr), value :: wall_tag_gpu
      type(c_ptr), value :: w_aux_gpu,meta_gpu,csimod_gpu,dxdcsic2_gpu,dydcsic2_gpu
      real(c_rkind) :: n,a,pn,pa,tn,ta
      type(c_ptr), value :: redn_3d_gpu

    endsubroutine compute_airfoil_forces_runtime_kernel_wrapper
  endinterface

  interface
    subroutine compute_rho_t_p_minmax_kernel1_wrapper(stream,nx,ny,nz,ng,rgas0,fluid_mask_gpu,&
    &w_aux_gpu,rhomin,tmin,pmin,redn_3d_gpu)bind(c,name="compute_rho_t_p_minmax_kernel1_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng
      real(c_rkind), value :: rgas0
      type(c_ptr), value :: fluid_mask_gpu
      type(c_ptr), value :: w_aux_gpu
      real(c_rkind) :: rhomin,tmin,pmin
      type(c_ptr), value :: redn_3d_gpu

    endsubroutine compute_rho_t_p_minmax_kernel1_wrapper
  endinterface

  interface
    subroutine compute_rho_t_p_minmax_kernel2_wrapper(stream,nx,ny,nz,ng,rgas0,fluid_mask_gpu,&
    &w_aux_gpu,rhomax,tmax,pmax,redn_3d_gpu)bind(c,name="compute_rho_t_p_minmax_kernel2_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng
      real(c_rkind), value :: rgas0
      type(c_ptr), value :: fluid_mask_gpu
      type(c_ptr), value :: w_aux_gpu
      real(c_rkind) :: rhomax,tmax,pmax
      type(c_ptr), value :: redn_3d_gpu

    endsubroutine compute_rho_t_p_minmax_kernel2_wrapper
  endinterface

  interface
    subroutine compute_dt_kernel_wrapper(stream,nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect,&
    &rgas0,t0,prandtl,fluid_mask_gpu,w_gpu,w_aux_gpu,cp_coeff_gpu,dcsidx_gpu,detady_gpu,dzitdz_gpu,&
    &dcsidxs_gpu,detadys_gpu,dzitdzs_gpu,dtxi_max,dtyi_max,dtzi_max,dtxv_max,dtyv_max,dtzv_max,dtxk_max,&
    &dtyk_max,dtzk_max,redn_3d_gpu)bind(c,name="compute_dt_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect
      real(c_rkind), value :: rgas0,t0,prandtl
      type(c_ptr), value :: fluid_mask_gpu
      type(c_ptr), value :: w_gpu,w_aux_gpu,cp_coeff_gpu,dcsidx_gpu,detady_gpu,dzitdz_gpu,dcsidxs_gpu,detadys_gpu,dzitdzs_gpu
      real(c_rkind) :: dtxi_max,dtyi_max,dtzi_max,dtxv_max,dtyv_max,dtzv_max,dtxk_max,dtyk_max,dtzk_max
      type(c_ptr), value :: redn_3d_gpu

    endsubroutine compute_dt_kernel_wrapper
  endinterface

  interface
    subroutine compute_dt_c2_kernel_wrapper(stream,nx,ny,nz,ng,indx_cp_l,indx_cp_r,&
    &calorically_perfect,rgas0,t0,prandtl,fluid_mask_gpu,w_gpu,w_aux_gpu,cp_coeff_gpu,dzitdz_gpu,&
    &dzitdzs_gpu,dcsidxnc2_gpu,dcsidync2_gpu,detadxnc2_gpu,detadync2_gpu,mcsi_gpu,meta_gpu,dtxi_max,&
    &dtyi_max,dtzi_max,dtxv_max,dtyv_max,dtzv_max,dtxk_max,dtyk_max,dtzk_max,redn_3d_gpu)bind(c,&
    &name="compute_dt_c2_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect
      real(c_rkind), value :: rgas0,t0,prandtl
      type(c_ptr), value :: fluid_mask_gpu
      type(c_ptr), value :: w_gpu,w_aux_gpu,cp_coeff_gpu,dzitdz_gpu,dzitdzs_gpu,dcsidxnc2_gpu,&
      &dcsidync2_gpu,detadxnc2_gpu,detadync2_gpu,mcsi_gpu,meta_gpu
      real(c_rkind) :: dtxi_max,dtyi_max,dtzi_max,dtxv_max,dtyv_max,dtzv_max,dtxk_max,dtyk_max,dtzk_max
      type(c_ptr), value :: redn_3d_gpu

    endsubroutine compute_dt_c2_kernel_wrapper
  endinterface

  interface
    subroutine eval_aux_kernel_wrapper(stream,nx,ny,nz,ng,visc_model,istart,iend,jstart,jend,kstart,&
    &kend,visc_power,visc_sutherland,visc_no,indx_cp_l,indx_cp_r,calorically_perfect,mu0,t0,sutherland_s,&
    &t_ref_dim,powerlaw_vtexp,rgas0,tol_iter_nr,cv_coeff_gpu,w_gpu,w_aux_gpu)bind(c,name="eval_aux_kernel&
    &_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,visc_model,istart,iend,jstart,jend,kstart,kend,&
      &visc_power,visc_sutherland,visc_no,indx_cp_l,indx_cp_r,calorically_perfect
      real(c_rkind), value :: mu0,t0,sutherland_s,t_ref_dim,powerlaw_vtexp,rgas0,tol_iter_nr
      type(c_ptr), value :: cv_coeff_gpu,w_gpu,w_aux_gpu

    endsubroutine eval_aux_kernel_wrapper
  endinterface

  interface
    subroutine tripping_pressure_kernel_wrapper(stream,nx,ny,nz,nv,ng,itr1,itr2,i_rank_start,pi,&
    &x0tr,y0tr,x0ts,y0ts,lamx,lamy,lamz,lamz1,phiz,phiz1,asl,bt,wall_tag_gpu,w_gpu,fl_gpu,dxdetanc2_gpu,&
    &dydetanc2_gpu,xc2_gpu,yc2_gpu,z_gpu)bind(c,name="tripping_pressure_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,nv,ng,itr1,itr2,i_rank_start
      real(c_rkind), value :: pi,x0tr,y0tr,x0ts,y0ts,lamx,lamy,lamz,lamz1,phiz,phiz1,asl,bt
      type(c_ptr), value :: wall_tag_gpu
      type(c_ptr), value :: w_gpu,fl_gpu,dxdetanc2_gpu,dydetanc2_gpu,xc2_gpu,yc2_gpu,z_gpu

    endsubroutine tripping_pressure_kernel_wrapper
  endinterface

  interface
    subroutine tripping_suction_kernel_wrapper(stream,nx,ny,nz,nv,ng,its1,its2,i_rank_start,pi,x0tr,&
    &y0tr,x0ts,y0ts,lamx,lamy,lams,lams1,phis,phis1,asl,bt,wall_tag_gpu,w_gpu,fl_gpu,dxdetanc2_gpu,&
    &dydetanc2_gpu,xc2_gpu,yc2_gpu,z_gpu)bind(c,name="tripping_suction_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,nv,ng,its1,its2,i_rank_start
      real(c_rkind), value :: pi,x0tr,y0tr,x0ts,y0ts,lamx,lamy,lams,lams1,phis,phis1,asl,bt
      type(c_ptr), value :: wall_tag_gpu
      type(c_ptr), value :: w_gpu,fl_gpu,dxdetanc2_gpu,dydetanc2_gpu,xc2_gpu,yc2_gpu,z_gpu

    endsubroutine tripping_suction_kernel_wrapper
  endinterface

  interface
    subroutine insitu_div_kernel_wrapper(stream,nx,ny,nz,ng,visc_order,npsi,mpsi,lmax,w_aux_gpu,&
    &coeff_deriv1_gpu,dcsidx_gpu,detady_gpu,dzitdz_gpu,psi_gpu)bind(c,name="insitu_div_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,visc_order,npsi,mpsi,lmax
      type(c_ptr), value :: w_aux_gpu,coeff_deriv1_gpu,dcsidx_gpu,detady_gpu,dzitdz_gpu,psi_gpu

    endsubroutine insitu_div_kernel_wrapper
  endinterface

  interface
    subroutine insitu_omega_kernel_wrapper(stream,nx,ny,nz,ng,visc_order,npsi,mpsi,lmax,w_aux_gpu,&
    &coeff_deriv1_gpu,dcsidx_gpu,detady_gpu,dzitdz_gpu,psi_gpu)bind(c,name="insitu_omega_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,visc_order,npsi,mpsi,lmax
      type(c_ptr), value :: w_aux_gpu,coeff_deriv1_gpu,dcsidx_gpu,detady_gpu,dzitdz_gpu,psi_gpu

    endsubroutine insitu_omega_kernel_wrapper
  endinterface

  interface
    subroutine insitu_ducros_kernel_wrapper(stream,nx,ny,nz,ng,visc_order,npsi,mpsi,lmax,u0,l0,&
    &w_aux_gpu,coeff_deriv1_gpu,dcsidx_gpu,detady_gpu,dzitdz_gpu,psi_gpu)bind(c,name="insitu_ducros_kerne&
    &l_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,visc_order,npsi,mpsi,lmax
      real(c_rkind), value :: u0,l0
      type(c_ptr), value :: w_aux_gpu,coeff_deriv1_gpu,dcsidx_gpu,detady_gpu,dzitdz_gpu,psi_gpu

    endsubroutine insitu_ducros_kernel_wrapper
  endinterface

  interface
    subroutine probe_interpolation_kernel_wrapper(stream,num_probe,nx,ny,nz,ng,nv_aux,ijk_probe_gpu,&
    &w_aux_gpu,w_aux_probe_gpu,probe_coeff_gpu)bind(c,name="probe_interpolation_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: num_probe,nx,ny,nz,ng,nv_aux
      type(c_ptr), value :: ijk_probe_gpu
      type(c_ptr), value :: w_aux_gpu,w_aux_probe_gpu,probe_coeff_gpu

    endsubroutine probe_interpolation_kernel_wrapper
  endinterface

  interface
    subroutine compute_tspec_kernel_wrapper(stream,nx,ny,nz,ng,ndft,j_slice,i_win,it_win,pi,&
    &w_aux_gpu,w_tspec_gpu)bind(c,name="compute_tspec_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,ndft,j_slice,i_win,it_win
      real(c_rkind), value :: pi
      type(c_ptr), value :: w_aux_gpu,w_tspec_gpu

    endsubroutine compute_tspec_kernel_wrapper
  endinterface

  interface
    subroutine compute_psd_tspec_kernel_wrapper(stream,nx,ny,nz,ndft,i_win,dt_tspec,pi,w_tspec_gpu,&
    &w_psd_tspec_gpu)bind(c,name="compute_psd_tspec_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ndft,i_win
      real(c_rkind), value :: dt_tspec,pi
      type(c_ptr), value :: w_tspec_gpu,w_psd_tspec_gpu

    endsubroutine compute_psd_tspec_kernel_wrapper
  endinterface

  interface
    subroutine compute_wallprop_c2_kernel_wrapper(stream,nx,ny,nz,ng,wall_tag_gpu,w_aux_gpu,&
    &wallprop_gpu,meta_gpu,csimod_gpu,dxdcsic2_gpu,dydcsic2_gpu)bind(c,name="compute_wallprop_c2_kernel_w&
    &rapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng
      type(c_ptr), value :: wall_tag_gpu
      type(c_ptr), value :: w_aux_gpu,wallprop_gpu,meta_gpu,csimod_gpu,dxdcsic2_gpu,dydcsic2_gpu

    endsubroutine compute_wallprop_c2_kernel_wrapper
  endinterface

  interface
    subroutine euler_x_fluxes_hybrid_c2_kernel_wrapper(stream,nv,nx,ny,nz,ng,nv_aux,eul_imin,&
    &eul_imax,lmax_base,nkeep,indx_cp_l,indx_cp_r,calorically_perfect,weno_scheme,weno_size,weno_version,&
    &force_zero_flux_min,force_zero_flux_max,sensor_threshold,rgas0,tol_iter_nr,rho0,u0,t0,&
    &ep_ord_change_gpu,lmax_tag_gpu,w_aux_gpu,fhat_gpu,cp_coeff_gpu,coeff_deriv1_gpu,dcsidxc2_gpu,&
    &dcsidxnc2_gpu,dcsidyc2_gpu,dcsidync2_gpu,jac_gpu,mcsijac1_gpu)bind(c,name="euler_x_fluxes_hybrid_c2_&
    &kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nv,nx,ny,nz,ng,nv_aux,eul_imin,eul_imax,lmax_base,nkeep,indx_cp_l,&
      &indx_cp_r,calorically_perfect,weno_scheme,weno_size,weno_version,force_zero_flux_min,&
      &force_zero_flux_max
      real(c_rkind), value :: sensor_threshold,rgas0,tol_iter_nr,rho0,u0,t0
      type(c_ptr), value :: ep_ord_change_gpu,lmax_tag_gpu
      type(c_ptr), value :: w_aux_gpu,fhat_gpu,cp_coeff_gpu,coeff_deriv1_gpu,dcsidxc2_gpu,&
      &dcsidxnc2_gpu,dcsidyc2_gpu,dcsidync2_gpu,jac_gpu,mcsijac1_gpu

    endsubroutine euler_x_fluxes_hybrid_c2_kernel_wrapper
  endinterface

  interface
    subroutine euler_x_fluxes_hybrid_kernel_wrapper(stream,nv,nx,ny,nz,ng,nv_aux,istart_face,&
    &iend_face,lmax_base,nkeep,indx_cp_l,indx_cp_r,calorically_perfect,weno_scheme,weno_size,&
    &weno_version,force_zero_flux_min,force_zero_flux_max,sensor_threshold,rgas0,tol_iter_nr,rho0,u0,t0,&
    &ep_ord_change_gpu,w_aux_gpu,fhat_gpu,cp_coeff_gpu,coeff_deriv1_gpu,dcsidx_gpu)bind(c,&
    &name="euler_x_fluxes_hybrid_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nv,nx,ny,nz,ng,nv_aux,istart_face,iend_face,lmax_base,nkeep,&
      &indx_cp_l,indx_cp_r,calorically_perfect,weno_scheme,weno_size,weno_version,force_zero_flux_min,&
      &force_zero_flux_max
      real(c_rkind), value :: sensor_threshold,rgas0,tol_iter_nr,rho0,u0,t0
      type(c_ptr), value :: ep_ord_change_gpu
      type(c_ptr), value :: w_aux_gpu,fhat_gpu,cp_coeff_gpu,coeff_deriv1_gpu,dcsidx_gpu

    endsubroutine euler_x_fluxes_hybrid_kernel_wrapper
  endinterface

  interface
    subroutine euler_x_fluxes_hybrid_rusanov_kernel_wrapper(stream,nv,nx,ny,nz,ng,nv_aux,&
    &istart_face,iend_face,lmax_base,nkeep,indx_cp_l,indx_cp_r,calorically_perfect,weno_scheme,weno_size,&
    &weno_version,force_zero_flux_min,force_zero_flux_max,sensor_threshold,rgas0,tol_iter_nr,rho0,u0,t0,&
    &ep_ord_change_gpu,w_aux_gpu,fhat_gpu,cp_coeff_gpu,coeff_deriv1_gpu,dcsidx_gpu)bind(c,&
    &name="euler_x_fluxes_hybrid_rusanov_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nv,nx,ny,nz,ng,nv_aux,istart_face,iend_face,lmax_base,nkeep,&
      &indx_cp_l,indx_cp_r,calorically_perfect,weno_scheme,weno_size,weno_version,force_zero_flux_min,&
      &force_zero_flux_max
      real(c_rkind), value :: sensor_threshold,rgas0,tol_iter_nr,rho0,u0,t0
      type(c_ptr), value :: ep_ord_change_gpu
      type(c_ptr), value :: w_aux_gpu,fhat_gpu,cp_coeff_gpu,coeff_deriv1_gpu,dcsidx_gpu

    endsubroutine euler_x_fluxes_hybrid_rusanov_kernel_wrapper
  endinterface

  interface
    subroutine euler_z_hybrid_kernel_wrapper(stream,nv,nx,ny,nz,ng,nv_aux,eul_kmin,eul_kmax,&
    &lmax_base,nkeep,indx_cp_l,indx_cp_r,calorically_perfect,force_zero_flux_min,force_zero_flux_max,&
    &weno_scheme,weno_size,weno_version,sensor_threshold,rgas0,tol_iter_nr,rho0,u0,t0,ep_ord_change_gpu,&
    &w_aux_gpu,fhat_gpu,cp_coeff_gpu,fl_gpu,coeff_deriv1_gpu,dzitdz_gpu)bind(c,name="euler_z_hybrid_kerne&
    &l_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nv,nx,ny,nz,ng,nv_aux,eul_kmin,eul_kmax,lmax_base,nkeep,indx_cp_l,&
      &indx_cp_r,calorically_perfect,force_zero_flux_min,force_zero_flux_max,weno_scheme,weno_size,&
      &weno_version
      real(c_rkind), value :: sensor_threshold,rgas0,tol_iter_nr,rho0,u0,t0
      type(c_ptr), value :: ep_ord_change_gpu
      type(c_ptr), value :: w_aux_gpu,fhat_gpu,cp_coeff_gpu,fl_gpu,coeff_deriv1_gpu,dzitdz_gpu

    endsubroutine euler_z_hybrid_kernel_wrapper
  endinterface

  interface
    subroutine euler_z_hybrid_rusanov_kernel_wrapper(stream,nv,nx,ny,nz,ng,nv_aux,eul_kmin,eul_kmax,&
    &lmax_base,nkeep,indx_cp_l,indx_cp_r,calorically_perfect,force_zero_flux_min,force_zero_flux_max,&
    &weno_scheme,weno_size,weno_version,sensor_threshold,rgas0,tol_iter_nr,rho0,u0,t0,ep_ord_change_gpu,&
    &w_aux_gpu,fhat_gpu,cp_coeff_gpu,fl_gpu,coeff_deriv1_gpu,dzitdz_gpu)bind(c,name="euler_z_hybrid_rusan&
    &ov_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nv,nx,ny,nz,ng,nv_aux,eul_kmin,eul_kmax,lmax_base,nkeep,indx_cp_l,&
      &indx_cp_r,calorically_perfect,force_zero_flux_min,force_zero_flux_max,weno_scheme,weno_size,&
      &weno_version
      real(c_rkind), value :: sensor_threshold,rgas0,tol_iter_nr,rho0,u0,t0
      type(c_ptr), value :: ep_ord_change_gpu
      type(c_ptr), value :: w_aux_gpu,fhat_gpu,cp_coeff_gpu,fl_gpu,coeff_deriv1_gpu,dzitdz_gpu

    endsubroutine euler_z_hybrid_rusanov_kernel_wrapper
  endinterface

  interface
    subroutine euler_y_hybrid_c2_kernel_wrapper(stream,nv,nx,ny,nz,ng,nv_aux,eul_jmin,eul_jmax,&
    &lmax_base,nkeep,indx_cp_l,indx_cp_r,calorically_perfect,force_zero_flux_min,force_zero_flux_max,&
    &weno_scheme,weno_size,weno_version,sensor_threshold,rgas0,tol_iter_nr,rho0,u0,t0,wall_tag_gpu,&
    &ep_ord_change_gpu,w_aux_gpu,fhat_gpu,cp_coeff_gpu,fl_gpu,coeff_deriv1_gpu,detadxc2_gpu,&
    &detadxnc2_gpu,detadyc2_gpu,detadync2_gpu,jac_gpu,metajac1_gpu)bind(c,name="euler_y_hybrid_c2_kernel_&
    &wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nv,nx,ny,nz,ng,nv_aux,eul_jmin,eul_jmax,lmax_base,nkeep,indx_cp_l,&
      &indx_cp_r,calorically_perfect,force_zero_flux_min,force_zero_flux_max,weno_scheme,weno_size,&
      &weno_version
      real(c_rkind), value :: sensor_threshold,rgas0,tol_iter_nr,rho0,u0,t0
      type(c_ptr), value :: wall_tag_gpu,ep_ord_change_gpu
      type(c_ptr), value :: w_aux_gpu,fhat_gpu,cp_coeff_gpu,fl_gpu,coeff_deriv1_gpu,detadxc2_gpu,&
      &detadxnc2_gpu,detadyc2_gpu,detadync2_gpu,jac_gpu,metajac1_gpu

    endsubroutine euler_y_hybrid_c2_kernel_wrapper
  endinterface

  interface
    subroutine euler_y_hybrid_kernel_wrapper(stream,nv,nx,ny,nz,ng,nv_aux,eul_jmin,eul_jmax,&
    &lmax_base,nkeep,indx_cp_l,indx_cp_r,calorically_perfect,force_zero_flux_min,force_zero_flux_max,&
    &weno_scheme,weno_size,weno_version,sensor_threshold,rgas0,tol_iter_nr,rho0,u0,t0,ep_ord_change_gpu,&
    &w_aux_gpu,fhat_gpu,cp_coeff_gpu,fl_gpu,coeff_deriv1_gpu,detady_gpu)bind(c,name="euler_y_hybrid_kerne&
    &l_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nv,nx,ny,nz,ng,nv_aux,eul_jmin,eul_jmax,lmax_base,nkeep,indx_cp_l,&
      &indx_cp_r,calorically_perfect,force_zero_flux_min,force_zero_flux_max,weno_scheme,weno_size,&
      &weno_version
      real(c_rkind), value :: sensor_threshold,rgas0,tol_iter_nr,rho0,u0,t0
      type(c_ptr), value :: ep_ord_change_gpu
      type(c_ptr), value :: w_aux_gpu,fhat_gpu,cp_coeff_gpu,fl_gpu,coeff_deriv1_gpu,detady_gpu

    endsubroutine euler_y_hybrid_kernel_wrapper
  endinterface

  interface
    subroutine euler_y_hybrid_rusanov_kernel_wrapper(stream,nv,nx,ny,nz,ng,nv_aux,eul_jmin,eul_jmax,&
    &lmax_base,nkeep,indx_cp_l,indx_cp_r,calorically_perfect,force_zero_flux_min,force_zero_flux_max,&
    &weno_scheme,weno_size,weno_version,sensor_threshold,rgas0,tol_iter_nr,rho0,u0,t0,ep_ord_change_gpu,&
    &w_aux_gpu,fhat_gpu,cp_coeff_gpu,fl_gpu,coeff_deriv1_gpu,detady_gpu)bind(c,name="euler_y_hybrid_rusan&
    &ov_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nv,nx,ny,nz,ng,nv_aux,eul_jmin,eul_jmax,lmax_base,nkeep,indx_cp_l,&
      &indx_cp_r,calorically_perfect,force_zero_flux_min,force_zero_flux_max,weno_scheme,weno_size,&
      &weno_version
      real(c_rkind), value :: sensor_threshold,rgas0,tol_iter_nr,rho0,u0,t0
      type(c_ptr), value :: ep_ord_change_gpu
      type(c_ptr), value :: w_aux_gpu,fhat_gpu,cp_coeff_gpu,fl_gpu,coeff_deriv1_gpu,detady_gpu

    endsubroutine euler_y_hybrid_rusanov_kernel_wrapper
  endinterface

  interface
    subroutine bc_nr_lat_x_kernel_wrapper(stream,start_or_end,nr_type,nx,ny,nz,ng,nv,indx_cp_l,&
    &indx_cp_r,calorically_perfect,rgas0,t0,cp_coeff_gpu,fl_gpu,w_aux_gpu,w_gpu,dcsidx_gpu,&
    &winf_gpu)bind(c,name="bc_nr_lat_x_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: start_or_end,nr_type,nx,ny,nz,ng,nv,indx_cp_l,indx_cp_r,calorically_perfect
      real(c_rkind), value :: rgas0,t0
      type(c_ptr), value :: cp_coeff_gpu,fl_gpu,w_aux_gpu,w_gpu,dcsidx_gpu,winf_gpu

    endsubroutine bc_nr_lat_x_kernel_wrapper
  endinterface

  interface
    subroutine bc_nr_lat_x_c2_kernel_wrapper(stream,start_or_end,nr_type,nx,ny,nz,ng,nv,indx_cp_l,&
    &indx_cp_r,calorically_perfect,rgas0,t0,cp_coeff_gpu,fl_gpu,w_aux_gpu,w_gpu,dcsidx_gpu,winf_gpu,&
    &jac_gpu,mcsi_gpu,dcsidxnc2_gpu,dcsidync2_gpu,dcsidxc2_gpu,dcsidyc2_gpu)bind(c,name="bc_nr_lat_x_c2_k&
    &ernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: start_or_end,nr_type,nx,ny,nz,ng,nv,indx_cp_l,indx_cp_r,calorically_perfect
      real(c_rkind), value :: rgas0,t0
      type(c_ptr), value :: cp_coeff_gpu,fl_gpu,w_aux_gpu,w_gpu,dcsidx_gpu,winf_gpu,jac_gpu,&
      &mcsi_gpu,dcsidxnc2_gpu,dcsidync2_gpu,dcsidxc2_gpu,dcsidyc2_gpu

    endsubroutine bc_nr_lat_x_c2_kernel_wrapper
  endinterface

  interface
    subroutine bc_nr_lat_y_kernel_wrapper(stream,start_or_end,nr_type,nx,ny,nz,ng,nv,indx_cp_l,&
    &indx_cp_r,calorically_perfect,rgas0,t0,cp_coeff_gpu,fl_gpu,w_aux_gpu,w_gpu,detady_gpu,&
    &winf_gpu)bind(c,name="bc_nr_lat_y_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: start_or_end,nr_type,nx,ny,nz,ng,nv,indx_cp_l,indx_cp_r,calorically_perfect
      real(c_rkind), value :: rgas0,t0
      type(c_ptr), value :: cp_coeff_gpu,fl_gpu,w_aux_gpu,w_gpu,detady_gpu,winf_gpu

    endsubroutine bc_nr_lat_y_kernel_wrapper
  endinterface

  interface
    subroutine bc_nr_lat_y_c2_kernel_wrapper(stream,start_or_end,nr_type,nx,ny,nz,ng,nv,indx_cp_l,&
    &indx_cp_r,calorically_perfect,rgas0,t0,wall_tag_gpu,cp_coeff_gpu,fl_gpu,w_aux_gpu,w_gpu,detady_gpu,&
    &wfar_gpu,jac_gpu,meta_gpu,detadxnc2_gpu,detadync2_gpu,detadxc2_gpu,detadyc2_gpu)bind(c,&
    &name="bc_nr_lat_y_c2_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: start_or_end,nr_type,nx,ny,nz,ng,nv,indx_cp_l,indx_cp_r,calorically_perfect
      real(c_rkind), value :: rgas0,t0
      type(c_ptr), value :: wall_tag_gpu
      type(c_ptr), value :: cp_coeff_gpu,fl_gpu,w_aux_gpu,w_gpu,detady_gpu,wfar_gpu,jac_gpu,&
      &meta_gpu,detadxnc2_gpu,detadync2_gpu,detadxc2_gpu,detadyc2_gpu

    endsubroutine bc_nr_lat_y_c2_kernel_wrapper
  endinterface

  interface
    subroutine bc_nr_lat_z_kernel_wrapper(stream,start_or_end,nr_type,nx,ny,nz,ng,nv,indx_cp_l,&
    &indx_cp_r,calorically_perfect,rgas0,t0,cp_coeff_gpu,fl_gpu,w_aux_gpu,w_gpu,dzitdz_gpu,&
    &winf_gpu)bind(c,name="bc_nr_lat_z_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: start_or_end,nr_type,nx,ny,nz,ng,nv,indx_cp_l,indx_cp_r,calorically_perfect
      real(c_rkind), value :: rgas0,t0
      type(c_ptr), value :: cp_coeff_gpu,fl_gpu,w_aux_gpu,w_gpu,dzitdz_gpu,winf_gpu

    endsubroutine bc_nr_lat_z_kernel_wrapper
  endinterface

  interface
    subroutine insitu_swirling_kernel_wrapper(stream,nv,nx,ny,nz,visc_order,ng,npsi,mpsi,u0,&
    &dcsidx_gpu,detady_gpu,dzitdz_gpu,psi_gpu,w_aux_gpu,coeff_deriv1_gpu,x_gpu)bind(c,&
    &name="insitu_swirling_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nv,nx,ny,nz,visc_order,ng,npsi,mpsi
      real(c_rkind), value :: u0
      type(c_ptr), value :: dcsidx_gpu,detady_gpu,dzitdz_gpu,psi_gpu,w_aux_gpu,coeff_deriv1_gpu,x_gpu

    endsubroutine insitu_swirling_kernel_wrapper
  endinterface

  interface
    subroutine insitu_swirling_c2_kernel_wrapper(stream,nv,nx,ny,nz,visc_order,ng,npsi,mpsi,u0,&
    &vis_tag_gpu,wall_tag_gpu,dzitdz_gpu,dcsidxc2_gpu,detadyc2_gpu,detadxc2_gpu,dcsidyc2_gpu,psi_gpu,&
    &w_aux_gpu,coeff_deriv1_gpu,x_gpu)bind(c,name="insitu_swirling_c2_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nv,nx,ny,nz,visc_order,ng,npsi,mpsi
      real(c_rkind), value :: u0
      type(c_ptr), value :: vis_tag_gpu,wall_tag_gpu
      type(c_ptr), value :: dzitdz_gpu,dcsidxc2_gpu,detadyc2_gpu,detadxc2_gpu,dcsidyc2_gpu,psi_gpu,w_aux_gpu,coeff_deriv1_gpu,x_gpu

    endsubroutine insitu_swirling_c2_kernel_wrapper
  endinterface

  interface
    subroutine insitu_schlieren_kernel_wrapper(stream,nv,nx,ny,nz,visc_order,ng,npsi,mpsi,u0,&
    &dcsidx_gpu,detady_gpu,dzitdz_gpu,psi_gpu,w_aux_gpu,coeff_deriv1_gpu,x_gpu)bind(c,&
    &name="insitu_schlieren_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nv,nx,ny,nz,visc_order,ng,npsi,mpsi
      real(c_rkind), value :: u0
      type(c_ptr), value :: dcsidx_gpu,detady_gpu,dzitdz_gpu,psi_gpu,w_aux_gpu,coeff_deriv1_gpu,x_gpu

    endsubroutine insitu_schlieren_kernel_wrapper
  endinterface

  interface
    subroutine insitu_schlieren_c2_kernel_wrapper(stream,nv,nx,ny,nz,visc_order,ng,npsi,mpsi,u0,&
    &vis_tag_gpu,wall_tag_gpu,dzitdz_gpu,dcsidxc2_gpu,detadyc2_gpu,detadxc2_gpu,dcsidyc2_gpu,psi_gpu,&
    &w_aux_gpu,coeff_deriv1_gpu,x_gpu)bind(c,name="insitu_schlieren_c2_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nv,nx,ny,nz,visc_order,ng,npsi,mpsi
      real(c_rkind), value :: u0
      type(c_ptr), value :: vis_tag_gpu,wall_tag_gpu
      type(c_ptr), value :: dzitdz_gpu,dcsidxc2_gpu,detadyc2_gpu,detadxc2_gpu,dcsidyc2_gpu,psi_gpu,w_aux_gpu,coeff_deriv1_gpu,x_gpu

    endsubroutine insitu_schlieren_c2_kernel_wrapper
  endinterface
contains

  subroutine zero_flux_kernel(nx,ny,nz,nv,fl_gpu)
    integer :: nx,ny,nz
    integer :: nv
    real(rkind), dimension(:,:,:,:), target :: fl_gpu


    call zero_flux_kernel_wrapper(c_null_ptr,nx,ny,nz,nv,c_loc(fl_gpu))
    call hipCheck(hipDeviceSynchronize())


  endsubroutine zero_flux_kernel


  subroutine init_flux_kernel(nx,ny,nz,nv,fl_gpu,fln_gpu,rhodt)
    integer :: nx,ny,nz
    integer :: nv
    real(rkind) :: rhodt
    real(rkind), dimension(:,:,:,:), target :: fl_gpu
    real(rkind), dimension(:,:,:,:), target :: fln_gpu


    call init_flux_kernel_wrapper(c_null_ptr,nx,ny,nz,nv,rhodt,c_loc(fl_gpu),c_loc(fln_gpu))
    call hipCheck(hipDeviceSynchronize())


  endsubroutine init_flux_kernel



  subroutine count_weno_kernel(nv,nv_aux,nx,ny,nz,ng,eul_imin,eul_imax,eul_jmin,eul_jmax,eul_kmin,&
  &eul_kmax,weno_scheme,sensor_threshold,w_aux_gpu,ep_ord_change_gpu,count_weno_x,count_weno_y,&
  &count_weno_z,redn_3d_gpu)
    integer :: nv,nv_aux,nx
    integer :: ny,nz,ng
    integer :: eul_imin,eul_imax,eul_jmin
    integer :: eul_jmax,eul_kmin,eul_kmax
    integer :: weno_scheme
    real(rkind) :: sensor_threshold
    real(rkind) :: count_weno_x,count_weno_y,count_weno_z
    real(rkind), dimension(:,:,:,:), target :: w_aux_gpu
    integer, dimension(:,:,:,:), target :: ep_ord_change_gpu
    real(rkind), dimension(:,:,:), target :: redn_3d_gpu

    count_weno_x = 0._rkind
    call count_weno_kernel1_wrapper(c_null_ptr,nv,nv_aux,nx,ny,nz,ng,eul_imin,eul_imax,eul_jmin,&
    &eul_jmax,eul_kmin,eul_kmax,weno_scheme,sensor_threshold,c_loc(ep_ord_change_gpu),c_loc(w_aux_gpu),&
    &count_weno_x,c_loc(redn_3d_gpu))
    call hipCheck(hipDeviceSynchronize())
    count_weno_y = 0._rkind
    call count_weno_kernel2_wrapper(c_null_ptr,nv,nv_aux,nx,ny,nz,ng,eul_imin,eul_imax,eul_jmin,&
    &eul_jmax,eul_kmin,eul_kmax,weno_scheme,sensor_threshold,c_loc(ep_ord_change_gpu),c_loc(w_aux_gpu),&
    &count_weno_y,c_loc(redn_3d_gpu))
    call hipCheck(hipDeviceSynchronize())
    count_weno_z = 0._rkind
    call count_weno_kernel3_wrapper(c_null_ptr,nv,nv_aux,nx,ny,nz,ng,eul_imin,eul_imax,eul_jmin,&
    &eul_jmax,eul_kmin,eul_kmax,weno_scheme,sensor_threshold,c_loc(ep_ord_change_gpu),c_loc(w_aux_gpu),&
    &count_weno_z,c_loc(redn_3d_gpu))
    call hipCheck(hipDeviceSynchronize())


  endsubroutine count_weno_kernel


  subroutine count_weno_c2_kernel(nv,nv_aux,nx,ny,nz,ng,eul_imin,eul_imax,eul_jmin,eul_jmax,&
  &eul_kmin,eul_kmax,l_base,weno_scheme,sensor_threshold,w_aux_gpu,ep_ord_change_gpu,lmax_tag_gpu,&
  &wall_tag_gpu,count_weno_x,count_weno_y,count_weno_z,redn_3d_gpu)
    integer :: nv,nv_aux,nx
    integer :: ny,nz,ng
    integer :: eul_imin,eul_imax,eul_jmin
    integer :: eul_jmax,eul_kmin,eul_kmax
    integer :: l_base,weno_scheme
    real(rkind) :: sensor_threshold
    real(rkind) :: count_weno_x,count_weno_y,count_weno_z
    real(rkind), dimension(:,:,:,:), target :: w_aux_gpu
    integer, dimension(:,:,:,:), target :: ep_ord_change_gpu
    integer, dimension(:), target :: lmax_tag_gpu
    integer, dimension(:), target :: wall_tag_gpu
    real(rkind), dimension(:,:,:), target :: redn_3d_gpu

    count_weno_x = 0._rkind
    call count_weno_c2_kernel1_wrapper(c_null_ptr,nv,nv_aux,nx,ny,nz,ng,eul_imin,eul_imax,eul_jmin,&
    &eul_jmax,eul_kmin,eul_kmax,l_base,weno_scheme,sensor_threshold,c_loc(ep_ord_change_gpu),&
    &c_loc(lmax_tag_gpu),c_loc(wall_tag_gpu),c_loc(w_aux_gpu),count_weno_x,c_loc(redn_3d_gpu))
    call hipCheck(hipDeviceSynchronize())
    count_weno_y = 0._rkind
    call count_weno_c2_kernel2_wrapper(c_null_ptr,nv,nv_aux,nx,ny,nz,ng,eul_imin,eul_imax,eul_jmin,&
    &eul_jmax,eul_kmin,eul_kmax,l_base,weno_scheme,sensor_threshold,c_loc(ep_ord_change_gpu),&
    &c_loc(lmax_tag_gpu),c_loc(wall_tag_gpu),c_loc(w_aux_gpu),count_weno_y,c_loc(redn_3d_gpu))
    call hipCheck(hipDeviceSynchronize())
    count_weno_z = 0._rkind
    call count_weno_c2_kernel3_wrapper(c_null_ptr,nv,nv_aux,nx,ny,nz,ng,eul_imin,eul_imax,eul_jmin,&
    &eul_jmax,eul_kmin,eul_kmax,l_base,weno_scheme,sensor_threshold,c_loc(ep_ord_change_gpu),&
    &c_loc(lmax_tag_gpu),c_loc(wall_tag_gpu),c_loc(w_aux_gpu),count_weno_z,c_loc(redn_3d_gpu))
    call hipCheck(hipDeviceSynchronize())


  endsubroutine count_weno_c2_kernel





  subroutine euler_x_update_kernel(nx,ny,nz,ng,nv,eul_imin,eul_imax,fhat_gpu,fl_gpu,dcsidx_gpu,stream_id)
    integer :: nx,ny,nz
    integer :: ng,nv,eul_imin
    integer :: eul_imax
    real(rkind), dimension(:,:,:,:), target :: fhat_gpu
    real(rkind), dimension(:,:,:,:), target :: fl_gpu
    real(rkind), dimension(:), target :: dcsidx_gpu

    type(c_ptr), value :: stream_id

    call euler_x_update_kernel_wrapper(stream_id,nx,ny,nz,ng,nv,eul_imin,eul_imax,c_loc(fhat_gpu),c_loc(fl_gpu),c_loc(dcsidx_gpu))


  endsubroutine euler_x_update_kernel


  subroutine euler_x_update_c2_kernel(nx,ny,nz,ng,nv,eul_imin,eul_imax,fhat_gpu,fl_gpu,jac_gpu,stream_id)
    integer :: nx,ny,nz
    integer :: ng,nv,eul_imin
    integer :: eul_imax
    real(rkind), dimension(:,:,:,:), target :: fhat_gpu
    real(rkind), dimension(:,:,:,:), target :: fl_gpu
    real(rkind), dimension(:,:), target :: jac_gpu

    type(c_ptr), value :: stream_id

    call euler_x_update_c2_kernel_wrapper(stream_id,nx,ny,nz,ng,nv,eul_imin,eul_imax,c_loc(fhat_gpu),c_loc(fl_gpu),c_loc(jac_gpu))


  endsubroutine euler_x_update_c2_kernel



  subroutine euler_y_update_kernel(nx,ny,nz,ng,nv,eul_jmin,eul_jmax,fhat_gpu,fl_gpu,detady_gpu,stream_id)
    integer :: nx,ny,nz
    integer :: ng,nv,eul_jmin
    integer :: eul_jmax
    real(rkind), dimension(:,:,:,:), target :: fhat_gpu
    real(rkind), dimension(:,:,:,:), target :: fl_gpu
    real(rkind), dimension(:), target :: detady_gpu

    type(c_ptr), value :: stream_id

    call euler_y_update_kernel_wrapper(stream_id,nx,ny,nz,ng,nv,eul_jmin,eul_jmax,c_loc(fhat_gpu),c_loc(fl_gpu),c_loc(detady_gpu))


  endsubroutine euler_y_update_kernel


  subroutine euler_y_update_c2_kernel(nx,ny,nz,ng,nv,eul_jmin,eul_jmax,fhat_gpu,fl_gpu,jac_gpu,wall_tag_gpu,stream_id)
    integer :: nx,ny,nz
    integer :: ng,nv,eul_jmin
    integer :: eul_jmax
    real(rkind), dimension(:,:,:,:), target :: fhat_gpu
    real(rkind), dimension(:,:,:,:), target :: fl_gpu
    real(rkind), dimension(:,:), target :: jac_gpu
    integer, dimension(:), target :: wall_tag_gpu

    type(c_ptr), value :: stream_id

    call euler_y_update_c2_kernel_wrapper(stream_id,nx,ny,nz,ng,nv,eul_jmin,eul_jmax,&
    &c_loc(wall_tag_gpu),c_loc(fhat_gpu),c_loc(fl_gpu),c_loc(jac_gpu))


  endsubroutine euler_y_update_c2_kernel

  subroutine euler_z_update_kernel(nx,ny,nz,ng,nv,eul_kmin,eul_kmax,fhat_gpu,fl_gpu,dzitdz_gpu,stream_id)
    integer :: nx,ny,nz
    integer :: ng,nv,eul_kmin
    integer :: eul_kmax
    real(rkind), dimension(:,:,:,:), target :: fhat_gpu
    real(rkind), dimension(:,:,:,:), target :: fl_gpu
    real(rkind), dimension(:), target :: dzitdz_gpu

    type(c_ptr), value :: stream_id

    call euler_z_update_kernel_wrapper(stream_id,nx,ny,nz,ng,nv,eul_kmin,eul_kmax,c_loc(fhat_gpu),c_loc(fl_gpu),c_loc(dzitdz_gpu))


  endsubroutine euler_z_update_kernel


















  subroutine force_rhs_2_kernel(nx,ny,nz,ng,fln_gpu,w_aux_gpu,bulk5g,fluid_mask_gpu)
    integer :: nx,ny,nz
    integer :: ng
    real(rkind) :: bulk_1,bulk_2
    real(rkind), dimension(:,:,:,:), target :: fln_gpu
    real(rkind), dimension(:,:,:,:), target :: w_aux_gpu
    integer, dimension(:,:,:), target :: fluid_mask_gpu
    real(rkind), dimension(5), intent(in) :: bulk5g

    bulk_1 = bulk5g(1)
    bulk_2 = bulk5g(2)
    call force_rhs_2_kernel_wrapper(c_null_ptr,nx,ny,nz,ng,bulk_1,bulk_2,c_loc(fluid_mask_gpu),c_loc(fln_gpu),c_loc(w_aux_gpu))
    call hipCheck(hipDeviceSynchronize())


  endsubroutine force_rhs_2_kernel


  subroutine force_rhs_2_c2_kernel(nx,ny,nz,ng,fln_gpu,w_aux_gpu,bulk5g,fluid_mask_gpu,r_curv,yn_gpu,dxdcsinc2_gpu,dydcsinc2_gpu)
    integer :: nx,ny,nz
    integer :: ng
    real(rkind) :: r_curv,bulk_1,bulk_5
    real(rkind), dimension(:,:,:,:), target :: fln_gpu
    real(rkind), dimension(:), target :: yn_gpu
    real(rkind), dimension(:,:), target :: dxdcsinc2_gpu
    real(rkind), dimension(:,:), target :: dydcsinc2_gpu
    real(rkind), dimension(:,:,:,:), target :: w_aux_gpu
    integer, dimension(:,:,:), target :: fluid_mask_gpu
    real(rkind), dimension(5), intent(in) :: bulk5g

    bulk_1 = bulk5g(1)
    bulk_5 = bulk5g(5)
    call force_rhs_2_c2_kernel_wrapper(c_null_ptr,nx,ny,nz,ng,r_curv,bulk_1,bulk_5,&
    &c_loc(fluid_mask_gpu),c_loc(fln_gpu),c_loc(yn_gpu),c_loc(dxdcsinc2_gpu),c_loc(dydcsinc2_gpu),&
    &c_loc(w_aux_gpu))
    call hipCheck(hipDeviceSynchronize())


  endsubroutine force_rhs_2_c2_kernel

  subroutine force_rhs_1_kernel(nx,ny,nz,ng,yn_gpu,fln_gpu,w_gpu,w_aux_gpu,bulk5,fluid_mask_gpu,redn_3d_gpu)
    integer :: nx,ny,nz
    integer :: ng
    real(rkind) :: bulk_1,bulk_2,bulk_3
    real(rkind) :: bulk_4,bulk_5
    real(rkind), dimension(:), target :: yn_gpu
    real(rkind), dimension(:,:,:,:), target :: fln_gpu
    real(rkind), dimension(:,:,:,:), target :: w_gpu
    real(rkind), dimension(:,:,:,:), target :: w_aux_gpu
    integer, dimension(:,:,:), target :: fluid_mask_gpu
    real(rkind), dimension(:,:,:), target :: redn_3d_gpu
    real(rkind), dimension(5), intent(out) :: bulk5

    bulk_1 = 0._rkind
    bulk_2 = 0._rkind
    bulk_3 = 0._rkind
    bulk_4 = 0._rkind
    bulk_5 = 0._rkind
    call force_rhs_1_kernel_wrapper(c_null_ptr,nx,ny,nz,ng,c_loc(fluid_mask_gpu),c_loc(yn_gpu),&
    &c_loc(fln_gpu),c_loc(w_gpu),c_loc(w_aux_gpu),bulk_1,bulk_2,bulk_3,bulk_4,bulk_5,c_loc(redn_3d_gpu))
    call hipCheck(hipDeviceSynchronize())
    bulk5(1) = bulk_1
    bulk5(2) = bulk_2
    bulk5(3) = bulk_3
    bulk5(4) = bulk_4
    bulk5(5) = bulk_5

  endsubroutine force_rhs_1_kernel


  subroutine force_rhs_1_c2_kernel(nx,ny,nz,ng,yn_gpu,fln_gpu,w_gpu,w_aux_gpu,bulk5,fluid_mask_gpu,&
  &dcsidxnc2_gpu,dcsidync2_gpu,jac_gpu,redn_3d_gpu)
    integer :: nx,ny,nz
    integer :: ng
    real(rkind) :: bulk_1,bulk_2,bulk_3
    real(rkind) :: bulk_4,bulk_5
    real(rkind), dimension(:), target :: yn_gpu
    real(rkind), dimension(:,:,:,:), target :: fln_gpu
    real(rkind), dimension(:,:), target :: dcsidxnc2_gpu
    real(rkind), dimension(:,:), target :: dcsidync2_gpu
    real(rkind), dimension(:,:), target :: jac_gpu
    real(rkind), dimension(:,:,:,:), target :: w_gpu
    real(rkind), dimension(:,:,:,:), target :: w_aux_gpu
    integer, dimension(:,:,:), target :: fluid_mask_gpu
    real(rkind), dimension(:,:,:), target :: redn_3d_gpu
    real(rkind), dimension(5), intent(out) :: bulk5

    bulk_1 = 0._rkind
    bulk_2 = 0._rkind
    bulk_3 = 0._rkind
    bulk_4 = 0._rkind
    bulk_5 = 0._rkind
    call force_rhs_1_c2_kernel_wrapper(c_null_ptr,nx,ny,nz,ng,c_loc(fluid_mask_gpu),c_loc(yn_gpu),&
    &c_loc(fln_gpu),c_loc(dcsidxnc2_gpu),c_loc(dcsidync2_gpu),c_loc(jac_gpu),c_loc(w_gpu),&
    &c_loc(w_aux_gpu),bulk_1,bulk_2,bulk_3,bulk_4,bulk_5,c_loc(redn_3d_gpu))
    call hipCheck(hipDeviceSynchronize())
    bulk5(1) = bulk_1
    bulk5(2) = bulk_2
    bulk5(3) = bulk_3
    bulk5(4) = bulk_4
    bulk5(5) = bulk_5

  endsubroutine force_rhs_1_c2_kernel


  subroutine force_var_1_kernel(nx,ny,nz,ng,yn_gpu,fln_gpu,w_gpu,w_aux_gpu,bulkt,fluid_mask_gpu,&
  &cv_coeff_gpu,indx_cp_l,indx_cp_r,t0,calorically_perfect,tol_iter_nr,redn_3d_gpu)
    integer :: nx,ny,nz
    integer :: ng,indx_cp_l,indx_cp_r
    integer :: calorically_perfect
    real(rkind) :: t0,tol_iter_nr,bulkt
    real(rkind) :: bulk_5
    real(rkind), dimension(:), target :: yn_gpu
    real(rkind), dimension(:,:,:,:), target :: fln_gpu
    real(rkind), dimension(:,:,:,:), target :: w_gpu
    real(rkind), dimension(:,:,:,:), target :: w_aux_gpu
    real(rkind), dimension(:), target :: cv_coeff_gpu
    integer, dimension(:,:,:), target :: fluid_mask_gpu
    real(rkind), dimension(:,:,:), target :: redn_3d_gpu

    bulk_5 = 0._rkind
    call force_var_1_kernel_wrapper(c_null_ptr,nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect,&
    &t0,tol_iter_nr,bulkt,c_loc(fluid_mask_gpu),c_loc(yn_gpu),c_loc(fln_gpu),c_loc(w_gpu),&
    &c_loc(w_aux_gpu),c_loc(cv_coeff_gpu),bulk_5,c_loc(redn_3d_gpu))
    call hipCheck(hipDeviceSynchronize())
    bulkt = bulk_5

  endsubroutine force_var_1_kernel

  subroutine force_var_2_kernel(nx,ny,nz,ng,w_gpu,w_aux_gpu,tbdiff,fluid_mask_gpu,cv_coeff_gpu,&
  &indx_cp_l,indx_cp_r,t0,calorically_perfect)
    integer :: nx,ny,nz
    integer :: ng,indx_cp_l,indx_cp_r
    integer :: calorically_perfect
    real(rkind) :: tbdiff,t0
    real(rkind), dimension(:,:,:,:), target :: w_gpu
    real(rkind), dimension(:,:,:,:), target :: w_aux_gpu
    real(rkind), dimension(:), target :: cv_coeff_gpu
    integer, dimension(:,:,:), target :: fluid_mask_gpu


    call force_var_2_kernel_wrapper(c_null_ptr,nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect,&
    &tbdiff,t0,c_loc(fluid_mask_gpu),c_loc(w_gpu),c_loc(w_aux_gpu),c_loc(cv_coeff_gpu))
    call hipCheck(hipDeviceSynchronize())


  endsubroutine force_var_2_kernel


  subroutine force_var_1_c2_kernel(nx,ny,nz,ng,yn_gpu,fln_gpu,w_gpu,w_aux_gpu,bulkt,fluid_mask_gpu,&
  &cv_coeff_gpu,indx_cp_l,indx_cp_r,t0,calorically_perfect,tol_iter_nr,jac_gpu,dcsidxnc2_gpu,&
  &dcsidync2_gpu,redn_3d_gpu)
    integer :: nx,ny,nz
    integer :: ng,indx_cp_l,indx_cp_r
    integer :: calorically_perfect
    real(rkind) :: t0,tol_iter_nr,bulkt
    real(rkind) :: bulk_5
    real(rkind), dimension(:), target :: yn_gpu
    real(rkind), dimension(:,:,:,:), target :: fln_gpu
    real(rkind), dimension(:,:,:,:), target :: w_gpu
    real(rkind), dimension(:,:,:,:), target :: w_aux_gpu
    real(rkind), dimension(:), target :: cv_coeff_gpu
    real(rkind), dimension(:,:), target :: jac_gpu
    real(rkind), dimension(:,:), target :: dcsidxnc2_gpu
    real(rkind), dimension(:,:), target :: dcsidync2_gpu
    integer, dimension(:,:,:), target :: fluid_mask_gpu
    real(rkind), dimension(:,:,:), target :: redn_3d_gpu

    bulk_5 = 0._rkind
    call force_var_1_c2_kernel_wrapper(c_null_ptr,nx,ny,nz,ng,indx_cp_l,indx_cp_r,&
    &calorically_perfect,t0,tol_iter_nr,bulkt,c_loc(fluid_mask_gpu),c_loc(yn_gpu),c_loc(fln_gpu),&
    &c_loc(w_gpu),c_loc(w_aux_gpu),c_loc(cv_coeff_gpu),c_loc(jac_gpu),c_loc(dcsidxnc2_gpu),&
    &c_loc(dcsidync2_gpu),bulk_5,c_loc(redn_3d_gpu))
    call hipCheck(hipDeviceSynchronize())
    bulkt = bulk_5

  endsubroutine force_var_1_c2_kernel

  subroutine force_var_2_c2_kernel(nx,ny,nz,ng,w_gpu,w_aux_gpu,tbdiff,fluid_mask_gpu,cv_coeff_gpu,&
  &indx_cp_l,indx_cp_r,t0,calorically_perfect,jac_gpu)
    integer :: nx,ny,nz
    integer :: ng,indx_cp_l,indx_cp_r
    integer :: calorically_perfect
    real(rkind) :: tbdiff,t0
    real(rkind), dimension(:,:,:,:), target :: w_gpu
    real(rkind), dimension(:,:,:,:), target :: w_aux_gpu
    real(rkind), dimension(:), target :: cv_coeff_gpu
    real(rkind), dimension(:,:), target :: jac_gpu
    integer, dimension(:,:,:), target :: fluid_mask_gpu


    call force_var_2_c2_kernel_wrapper(c_null_ptr,nx,ny,nz,ng,indx_cp_l,indx_cp_r,&
    &calorically_perfect,tbdiff,t0,c_loc(fluid_mask_gpu),c_loc(w_gpu),c_loc(w_aux_gpu),&
    &c_loc(cv_coeff_gpu),c_loc(jac_gpu))
    call hipCheck(hipDeviceSynchronize())


  endsubroutine force_var_2_c2_kernel


  subroutine update_flux_kernel(nx,ny,nz,nv,fl_gpu,fln_gpu,gamdt)
    integer :: nx,ny,nz
    integer :: nv
    real(rkind) :: gamdt
    real(rkind), dimension(:,:,:,:), target :: fl_gpu
    real(rkind), dimension(:,:,:,:), target :: fln_gpu


    call update_flux_kernel_wrapper(c_null_ptr,nx,ny,nz,nv,gamdt,c_loc(fl_gpu),c_loc(fln_gpu))
    call hipCheck(hipDeviceSynchronize())


  endsubroutine update_flux_kernel


  subroutine update_field_kernel(nx,ny,nz,ng,nv,w_gpu,fln_gpu,fluid_mask_gpu)
    integer :: nx,ny,nz
    integer :: nv,ng
    real(rkind), dimension(:,:,:,:), target :: w_gpu
    real(rkind), dimension(:,:,:,:), target :: fln_gpu
    integer, dimension(:,:,:), target :: fluid_mask_gpu


    call update_field_kernel_wrapper(c_null_ptr,nx,ny,nz,nv,ng,c_loc(fluid_mask_gpu),c_loc(w_gpu),c_loc(fln_gpu))
    call hipCheck(hipDeviceSynchronize())


  endsubroutine update_field_kernel


  subroutine visflx_div_ord2_kernel(nx,ny,nz,ng,w_aux_gpu,fl_gpu,x_gpu,y_gpu,z_gpu,stream_id)
    integer :: nx,ny,nz
    integer :: ng
    real(rkind), dimension(:,:,:,:), target :: w_aux_gpu
    real(rkind), dimension(:,:,:,:), target :: fl_gpu
    real(rkind), dimension(:), target :: x_gpu
    real(rkind), dimension(:), target :: y_gpu
    real(rkind), dimension(:), target :: z_gpu

    type(c_ptr), value :: stream_id

    call visflx_div_ord2_kernel_wrapper(stream_id,nx,ny,nz,ng,c_loc(w_aux_gpu),c_loc(fl_gpu),&
    &c_loc(x_gpu),c_loc(y_gpu),c_loc(z_gpu))


  endsubroutine visflx_div_ord2_kernel


  subroutine visflx_div_kernel(nx,ny,nz,ng,visc_order,w_aux_gpu,fl_gpu,coeff_deriv1_gpu,dcsidx_gpu,detady_gpu,dzitdz_gpu,stream_id)
    integer :: nx,ny,nz
    integer :: ng,visc_order,lmax
    real(rkind), dimension(:,:,:,:), target :: w_aux_gpu
    real(rkind), dimension(:,:,:,:), target :: fl_gpu
    real(rkind), dimension(:,:), target :: coeff_deriv1_gpu
    real(rkind), dimension(:), target :: dcsidx_gpu
    real(rkind), dimension(:), target :: detady_gpu
    real(rkind), dimension(:), target :: dzitdz_gpu

    type(c_ptr), value :: stream_id
    lmax = visc_order/2
    call visflx_div_kernel_wrapper(stream_id,nx,ny,nz,ng,visc_order,lmax,c_loc(w_aux_gpu),&
    &c_loc(fl_gpu),c_loc(coeff_deriv1_gpu),c_loc(dcsidx_gpu),c_loc(detady_gpu),c_loc(dzitdz_gpu))


  endsubroutine visflx_div_kernel


  subroutine visflx_div_c2_kernel(nx,ny,nz,ng,visc_order,w_aux_gpu,fl_gpu,coeff_deriv1_gpu,&
  &dcsidxc2_gpu,dcsidyc2_gpu,detadxc2_gpu,detadyc2_gpu,dzitdz_gpu,vis_tag_gpu,wall_tag_gpu,stream_id)
    integer :: nx,ny,nz
    integer :: ng,visc_order,lmax
    real(rkind), dimension(:,:,:,:), target :: w_aux_gpu
    real(rkind), dimension(:,:,:,:), target :: fl_gpu
    real(rkind), dimension(:,:), target :: coeff_deriv1_gpu
    real(rkind), dimension(:,:), target :: dcsidxc2_gpu
    real(rkind), dimension(:,:), target :: detadyc2_gpu
    real(rkind), dimension(:,:), target :: dcsidyc2_gpu
    real(rkind), dimension(:,:), target :: detadxc2_gpu
    real(rkind), dimension(:), target :: dzitdz_gpu
    integer, dimension(:), target :: vis_tag_gpu
    integer, dimension(:), target :: wall_tag_gpu

    type(c_ptr), value :: stream_id
    lmax = visc_order/2
    call visflx_div_c2_kernel_wrapper(stream_id,nx,ny,nz,ng,visc_order,lmax,c_loc(vis_tag_gpu),&
    &c_loc(wall_tag_gpu),c_loc(w_aux_gpu),c_loc(fl_gpu),c_loc(coeff_deriv1_gpu),c_loc(dcsidxc2_gpu),&
    &c_loc(detadyc2_gpu),c_loc(dcsidyc2_gpu),c_loc(detadxc2_gpu),c_loc(dzitdz_gpu))


  endsubroutine visflx_div_c2_kernel


  subroutine visflx_kernel(nx,ny,nz,ng,visc_order,prandtl,t0,indx_cp_l,indx_cp_r,cp_coeff_gpu,&
  &calorically_perfect,u0,l0,w_gpu,w_aux_gpu,fl_gpu,coeff_deriv1_gpu,coeff_deriv2_gpu,dcsidx_gpu,&
  &detady_gpu,dzitdz_gpu,dcsidxs_gpu,detadys_gpu,dzitdzs_gpu,dcsidx2_gpu,detady2_gpu,dzitdz2_gpu,&
  &wallprop_gpu)
    integer :: nx,ny,nz
    integer :: ng,visc_order,calorically_perfect
    integer :: indx_cp_l,indx_cp_r,lmax
    real(rkind) :: prandtl,u0,l0
    real(rkind) :: t0
    real(rkind), dimension(:,:,:,:), target :: w_gpu
    real(rkind), dimension(:,:,:,:), target :: w_aux_gpu
    real(rkind), dimension(:,:,:), target :: wallprop_gpu
    real(rkind), dimension(:,:,:,:), target :: fl_gpu
    real(rkind), dimension(:,:), target :: coeff_deriv1_gpu
    real(rkind), dimension(:,:), target :: coeff_deriv2_gpu
    real(rkind), dimension(:), target :: dcsidx_gpu
    real(rkind), dimension(:), target :: detady_gpu
    real(rkind), dimension(:), target :: dzitdz_gpu
    real(rkind), dimension(:), target :: dcsidxs_gpu
    real(rkind), dimension(:), target :: detadys_gpu
    real(rkind), dimension(:), target :: dzitdzs_gpu
    real(rkind), dimension(:), target :: dcsidx2_gpu
    real(rkind), dimension(:), target :: detady2_gpu
    real(rkind), dimension(:), target :: dzitdz2_gpu
    real(rkind), dimension(:), target :: cp_coeff_gpu

    lmax = visc_order/2
    call visflx_kernel_wrapper(c_null_ptr,nx,ny,nz,ng,visc_order,calorically_perfect,indx_cp_l,&
    &indx_cp_r,lmax,prandtl,u0,l0,t0,c_loc(w_gpu),c_loc(w_aux_gpu),c_loc(wallprop_gpu),c_loc(fl_gpu),&
    &c_loc(coeff_deriv1_gpu),c_loc(coeff_deriv2_gpu),c_loc(dcsidx_gpu),c_loc(detady_gpu),&
    &c_loc(dzitdz_gpu),c_loc(dcsidxs_gpu),c_loc(detadys_gpu),c_loc(dzitdzs_gpu),c_loc(dcsidx2_gpu),&
    &c_loc(detady2_gpu),c_loc(dzitdz2_gpu),c_loc(cp_coeff_gpu))
    call hipCheck(hipDeviceSynchronize())


  endsubroutine visflx_kernel


  subroutine visflx_c2_kernel(nx,ny,nz,ng,visc_order,prandtl,t0,indx_cp_l,indx_cp_r,cp_coeff_gpu,&
  &calorically_perfect,u0,l0,w_gpu,w_aux_gpu,fl_gpu,coeff_deriv1_gpu,coeff_deriv2_gpu,dcsidxc2_gpu,&
  &detadyc2_gpu,detadxc2_gpu,dcsidyc2_gpu,dzitdz_gpu,dzitdzs_gpu,dzitdz2_gpu,g1_gpu,g2_gpu,g12_gpu,&
  &jac_gpu,wallprop_gpu,vis_tag_gpu,wall_tag_gpu,iblock,ite_rank_x,itu_rank_x,ite_l,itu_l,teshk)
    integer :: nx,ny,nz
    integer :: ng,visc_order,calorically_perfect
    integer :: indx_cp_l,indx_cp_r,iblock
    integer :: ite_rank_x,itu_rank_x,ite_l
    integer :: itu_l,lmax
    real(rkind) :: prandtl,u0,l0
    real(rkind) :: t0,teshk
    real(rkind), dimension(:,:,:,:), target :: w_gpu
    real(rkind), dimension(:,:,:,:), target :: w_aux_gpu
    real(rkind), dimension(:,:,:), target :: wallprop_gpu
    real(rkind), dimension(:,:,:,:), target :: fl_gpu
    real(rkind), dimension(:,:), target :: coeff_deriv1_gpu
    real(rkind), dimension(:,:), target :: coeff_deriv2_gpu
    real(rkind), dimension(:), target :: dzitdz_gpu
    real(rkind), dimension(:), target :: dzitdzs_gpu
    real(rkind), dimension(:), target :: dzitdz2_gpu
    real(rkind), dimension(:,:), target :: dcsidxc2_gpu
    real(rkind), dimension(:,:), target :: detadyc2_gpu
    real(rkind), dimension(:,:), target :: detadxc2_gpu
    real(rkind), dimension(:,:), target :: dcsidyc2_gpu
    real(rkind), dimension(:,:), target :: g1_gpu
    real(rkind), dimension(:,:), target :: g2_gpu
    real(rkind), dimension(:,:), target :: g12_gpu
    real(rkind), dimension(:,:), target :: jac_gpu
    real(rkind), dimension(:), target :: cp_coeff_gpu
    integer, dimension(:), target :: vis_tag_gpu
    integer, dimension(:), target :: wall_tag_gpu

    lmax = visc_order/2
    call visflx_c2_kernel_wrapper(c_null_ptr,nx,ny,nz,ng,visc_order,calorically_perfect,indx_cp_l,&
    &indx_cp_r,iblock,ite_rank_x,itu_rank_x,ite_l,itu_l,lmax,prandtl,u0,l0,t0,teshk,c_loc(vis_tag_gpu),&
    &c_loc(wall_tag_gpu),c_loc(w_gpu),c_loc(w_aux_gpu),c_loc(wallprop_gpu),c_loc(fl_gpu),&
    &c_loc(coeff_deriv1_gpu),c_loc(coeff_deriv2_gpu),c_loc(dzitdz_gpu),c_loc(dzitdzs_gpu),&
    &c_loc(dzitdz2_gpu),c_loc(dcsidxc2_gpu),c_loc(detadyc2_gpu),c_loc(detadxc2_gpu),c_loc(dcsidyc2_gpu),&
    &c_loc(g1_gpu),c_loc(g2_gpu),c_loc(g12_gpu),c_loc(jac_gpu),c_loc(cp_coeff_gpu))
    call hipCheck(hipDeviceSynchronize())


  endsubroutine visflx_c2_kernel



  subroutine visflx_nosensor_kernel(nx,ny,nz,ng,visc_order,prandtl,t0,indx_cp_l,indx_cp_r,&
  &cp_coeff_gpu,calorically_perfect,u0,l0,w_gpu,w_aux_gpu,fl_gpu,coeff_deriv1_gpu,coeff_deriv2_gpu,&
  &dcsidx_gpu,detady_gpu,dzitdz_gpu,dcsidxs_gpu,detadys_gpu,dzitdzs_gpu,dcsidx2_gpu,detady2_gpu,&
  &dzitdz2_gpu,wallprop_gpu)
    integer :: nx,ny,nz
    integer :: ng,visc_order,calorically_perfect
    integer :: indx_cp_l,indx_cp_r,lmax
    real(rkind) :: prandtl,u0,l0
    real(rkind) :: t0
    real(rkind), dimension(:,:,:,:), target :: w_gpu
    real(rkind), dimension(:,:,:,:), target :: w_aux_gpu
    real(rkind), dimension(:,:,:), target :: wallprop_gpu
    real(rkind), dimension(:,:,:,:), target :: fl_gpu
    real(rkind), dimension(:,:), target :: coeff_deriv1_gpu
    real(rkind), dimension(:,:), target :: coeff_deriv2_gpu
    real(rkind), dimension(:), target :: dcsidx_gpu
    real(rkind), dimension(:), target :: detady_gpu
    real(rkind), dimension(:), target :: dzitdz_gpu
    real(rkind), dimension(:), target :: dcsidxs_gpu
    real(rkind), dimension(:), target :: detadys_gpu
    real(rkind), dimension(:), target :: dzitdzs_gpu
    real(rkind), dimension(:), target :: dcsidx2_gpu
    real(rkind), dimension(:), target :: detady2_gpu
    real(rkind), dimension(:), target :: dzitdz2_gpu
    real(rkind), dimension(:), target :: cp_coeff_gpu

    lmax = visc_order/2
    call visflx_nosensor_kernel_wrapper(c_null_ptr,nx,ny,nz,ng,visc_order,calorically_perfect,&
    &indx_cp_l,indx_cp_r,lmax,prandtl,u0,l0,t0,c_loc(w_gpu),c_loc(w_aux_gpu),c_loc(wallprop_gpu),&
    &c_loc(fl_gpu),c_loc(coeff_deriv1_gpu),c_loc(coeff_deriv2_gpu),c_loc(dcsidx_gpu),c_loc(detady_gpu),&
    &c_loc(dzitdz_gpu),c_loc(dcsidxs_gpu),c_loc(detadys_gpu),c_loc(dzitdzs_gpu),c_loc(dcsidx2_gpu),&
    &c_loc(detady2_gpu),c_loc(dzitdz2_gpu),c_loc(cp_coeff_gpu))
    call hipCheck(hipDeviceSynchronize())


  endsubroutine visflx_nosensor_kernel


  subroutine visflx_nosensor_c2_kernel(nx,ny,nz,ng,visc_order,prandtl,t0,indx_cp_l,indx_cp_r,&
  &cp_coeff_gpu,calorically_perfect,iblock,ite_rank_x,itu_rank_x,ite_l,itu_l,u0,l0,w_gpu,w_aux_gpu,&
  &fl_gpu,coeff_deriv1_gpu,coeff_deriv2_gpu,dcsidxc2_gpu,detadyc2_gpu,detadxc2_gpu,dcsidyc2_gpu,&
  &dzitdz_gpu,dzitdzs_gpu,dzitdz2_gpu,g1_gpu,g2_gpu,g12_gpu,jac_gpu,wallprop_gpu,vis_tag_gpu,&
  &wall_tag_gpu,ortho)
    integer :: nx,ny,nz
    integer :: ng,visc_order,calorically_perfect
    integer :: indx_cp_l,indx_cp_r,ortho
    integer :: iblock,ite_rank_x,itu_rank_x
    integer :: ite_l,itu_l,lmax
    real(rkind) :: prandtl,u0,l0
    real(rkind) :: t0
    real(rkind), dimension(:,:,:,:), target :: w_gpu
    real(rkind), dimension(:,:,:,:), target :: w_aux_gpu
    real(rkind), dimension(:,:,:), target :: wallprop_gpu
    real(rkind), dimension(:,:,:,:), target :: fl_gpu
    real(rkind), dimension(:,:), target :: coeff_deriv1_gpu
    real(rkind), dimension(:,:), target :: coeff_deriv2_gpu
    real(rkind), dimension(:), target :: dzitdz_gpu
    real(rkind), dimension(:), target :: dzitdzs_gpu
    real(rkind), dimension(:), target :: dzitdz2_gpu
    real(rkind), dimension(:,:), target :: dcsidxc2_gpu
    real(rkind), dimension(:,:), target :: detadyc2_gpu
    real(rkind), dimension(:,:), target :: detadxc2_gpu
    real(rkind), dimension(:,:), target :: dcsidyc2_gpu
    real(rkind), dimension(:,:), target :: g1_gpu
    real(rkind), dimension(:,:), target :: g2_gpu
    real(rkind), dimension(:,:), target :: g12_gpu
    real(rkind), dimension(:,:), target :: jac_gpu
    real(rkind), dimension(:), target :: cp_coeff_gpu
    integer, dimension(:), target :: vis_tag_gpu
    integer, dimension(:), target :: wall_tag_gpu

    lmax = visc_order/2
    call visflx_nosensor_c2_kernel_wrapper(c_null_ptr,nx,ny,nz,ng,visc_order,calorically_perfect,&
    &indx_cp_l,indx_cp_r,ortho,iblock,ite_rank_x,itu_rank_x,ite_l,itu_l,lmax,prandtl,u0,l0,t0,&
    &c_loc(vis_tag_gpu),c_loc(wall_tag_gpu),c_loc(w_gpu),c_loc(w_aux_gpu),c_loc(wallprop_gpu),&
    &c_loc(fl_gpu),c_loc(coeff_deriv1_gpu),c_loc(coeff_deriv2_gpu),c_loc(dzitdz_gpu),c_loc(dzitdzs_gpu),&
    &c_loc(dzitdz2_gpu),c_loc(dcsidxc2_gpu),c_loc(detadyc2_gpu),c_loc(detadxc2_gpu),c_loc(dcsidyc2_gpu),&
    &c_loc(g1_gpu),c_loc(g2_gpu),c_loc(g12_gpu),c_loc(jac_gpu),c_loc(cp_coeff_gpu))
    call hipCheck(hipDeviceSynchronize())


  endsubroutine visflx_nosensor_c2_kernel


  subroutine sponge_kernel(nx,ny,nz,ng,nv,w_gpu,wfar_gpu,fln_gpu,f_sponge_gpu,j_sponge)
    integer :: nx,ny,nz
    integer :: ng,nv,j_sponge
    real(rkind), dimension(:,:,:,:), target :: w_gpu
    real(rkind), dimension(:,:,:,:), target :: fln_gpu
    real(rkind), dimension(:,:), target :: wfar_gpu
    real(rkind), dimension(:), target :: f_sponge_gpu


    call sponge_kernel_wrapper(c_null_ptr,nx,ny,nz,ng,nv,j_sponge,c_loc(w_gpu),c_loc(fln_gpu),c_loc(wfar_gpu),c_loc(f_sponge_gpu))
    call hipCheck(hipDeviceSynchronize())


  endsubroutine sponge_kernel


  subroutine limiter_kernel(nx,ny,nz,ng,w_gpu,w_aux_gpu,iblock,kblock,indx_cp_l,indx_cp_r,&
  &cv_coeff_gpu,calorically_perfect,t0,tol_iter_nr,rho_lim,tem_lim,rho_lim_rescale,tem_lim_rescale,&
  &redn_3d_gpu)
    integer :: nx,ny,nz
    integer :: ng,iblock,kblock
    integer :: calorically_perfect,indx_cp_l,indx_cp_r
    real(rkind) :: t0,tol_iter_nr,rho_lim
    real(rkind) :: tem_lim,rho_lim_rescale,tem_lim_rescale
    real(rkind) :: n_limited_rho,n_limited_tem
    real(rkind), dimension(:,:,:,:), target :: w_gpu
    real(rkind), dimension(:,:,:,:), target :: w_aux_gpu
    real(rkind), dimension(:), target :: cv_coeff_gpu
    real(rkind), dimension(:,:,:), target :: redn_3d_gpu

    n_limited_rho = 0._rkind
    n_limited_tem = 0._rkind
    call limiter_kernel1_wrapper(c_null_ptr,nx,ny,nz,ng,iblock,kblock,calorically_perfect,indx_cp_l,&
    &indx_cp_r,t0,tol_iter_nr,rho_lim,tem_lim,rho_lim_rescale,tem_lim_rescale,c_loc(w_gpu),&
    &c_loc(w_aux_gpu),c_loc(cv_coeff_gpu),n_limited_rho,n_limited_tem,c_loc(redn_3d_gpu))
    call hipCheck(hipDeviceSynchronize())
    if(n_limited_rho > 0.) print*,'warning! n_limited_rho :',n_limited_rho
    if(n_limited_tem > 0.) print*,'warning! n_limited_tem :',n_limited_tem

    if(n_limited_rho > 0. .or. n_limited_tem > 0.) then
      call limiter_kernel2_wrapper(c_null_ptr,nx,ny,nz,ng,iblock,kblock,calorically_perfect,&
      &indx_cp_l,indx_cp_r,t0,tol_iter_nr,rho_lim,tem_lim,rho_lim_rescale,tem_lim_rescale,c_loc(w_gpu),&
      &c_loc(w_aux_gpu),c_loc(cv_coeff_gpu))
      call hipCheck(hipDeviceSynchronize())
    endif

  endsubroutine limiter_kernel


  subroutine filter_kernel(nx,ny,nz,ng,w_gpu,w_aux_gpu,indx_cp_l,indx_cp_r,cv_coeff_gpu,&
  &calorically_perfect,t0,tol_iter_nr,coeff_filter_gpu,jfilter,wall_tag_gpu)
    integer :: nx,ny,nz
    integer :: ng,calorically_perfect,indx_cp_l
    integer :: indx_cp_r,jfilter
    real(rkind) :: t0,tol_iter_nr
    real(rkind), dimension(:,:,:,:), target :: w_gpu
    real(rkind), dimension(:,:,:,:), target :: w_aux_gpu
    real(rkind), dimension(:), target :: cv_coeff_gpu
    real(rkind), dimension(:), target :: coeff_filter_gpu
    integer, dimension(:), target :: wall_tag_gpu


    call filter_kernel1_wrapper(c_null_ptr,nx,ny,nz,ng,calorically_perfect,indx_cp_l,indx_cp_r,&
    &jfilter,t0,tol_iter_nr,c_loc(wall_tag_gpu),c_loc(w_gpu),c_loc(w_aux_gpu),c_loc(cv_coeff_gpu),&
    &c_loc(coeff_filter_gpu))
    call hipCheck(hipDeviceSynchronize())

    call filter_kernel2_wrapper(c_null_ptr,nx,ny,nz,ng,calorically_perfect,indx_cp_l,indx_cp_r,&
    &jfilter,t0,tol_iter_nr,c_loc(wall_tag_gpu),c_loc(w_gpu),c_loc(w_aux_gpu),c_loc(cv_coeff_gpu),&
    &c_loc(coeff_filter_gpu))
    call hipCheck(hipDeviceSynchronize())


  endsubroutine filter_kernel


  subroutine visflx_reduced_ord2_kernel(nx,ny,nz,ng,prandtl,t0,indx_cp_l,indx_cp_r,cp_coeff_gpu,&
  &calorically_perfect,u0,l0,w_gpu,w_aux_gpu,fl_gpu,x_gpu,y_gpu,z_gpu,wallprop_gpu,update_sensor)
    integer :: nx,ny,nz
    integer :: ng,calorically_perfect,indx_cp_l
    integer :: indx_cp_r,update_sensor
    real(rkind) :: prandtl,u0,l0
    real(rkind) :: t0
    real(rkind), dimension(:,:,:,:), target :: w_gpu
    real(rkind), dimension(:,:,:,:), target :: w_aux_gpu
    real(rkind), dimension(:,:,:), target :: wallprop_gpu
    real(rkind), dimension(:,:,:,:), target :: fl_gpu
    real(rkind), dimension(:), target :: x_gpu
    real(rkind), dimension(:), target :: y_gpu
    real(rkind), dimension(:), target :: z_gpu
    real(rkind), dimension(:), target :: cp_coeff_gpu


    call visflx_reduced_ord2_kernel_wrapper(c_null_ptr,nx,ny,nz,ng,calorically_perfect,indx_cp_l,&
    &indx_cp_r,update_sensor,prandtl,u0,l0,t0,c_loc(w_gpu),c_loc(w_aux_gpu),c_loc(wallprop_gpu),&
    &c_loc(fl_gpu),c_loc(x_gpu),c_loc(y_gpu),c_loc(z_gpu),c_loc(cp_coeff_gpu))
    call hipCheck(hipDeviceSynchronize())


  endsubroutine visflx_reduced_ord2_kernel



  subroutine sensor_kernel(nx,ny,nz,ng,u0,l0,w_aux_gpu)
    integer :: nx,ny,nz
    integer :: ng
    real(rkind) :: u0,l0
    real(rkind), dimension(:,:,:,:), target :: w_aux_gpu


    call sensor_kernel_wrapper(c_null_ptr,nx,ny,nz,ng,u0,l0,c_loc(w_aux_gpu))
    call hipCheck(hipDeviceSynchronize())


  endsubroutine sensor_kernel


  subroutine sensor_c2_kernel(nx,ny,nz,ng,u0,l0_ducros,w_aux_gpu,iblock,ite_rank_x,itu_rank_x,ite_l,&
  &itu_l,teshk,theta_ij_gpu,theta_threshold,jweno)
    integer :: nx,ny,nz
    integer :: ng,iblock,ite_rank_x
    integer :: itu_rank_x,ite_l,itu_l
    integer :: jweno
    real(rkind) :: u0,l0_ducros,teshk
    real(rkind) :: theta_threshold
    real(rkind), dimension(:,:,:,:), target :: w_aux_gpu
    real(rkind), dimension(:,:), target :: theta_ij_gpu


    call sensor_c2_kernel_wrapper(c_null_ptr,nx,ny,nz,ng,iblock,ite_rank_x,itu_rank_x,ite_l,itu_l,&
    &jweno,u0,l0_ducros,teshk,theta_threshold,c_loc(w_aux_gpu),c_loc(theta_ij_gpu))
    call hipCheck(hipDeviceSynchronize())


  endsubroutine sensor_c2_kernel


  subroutine visflx_x_kernel(nx,ny,nz,nv,ng,prandtl,t0,indx_cp_l,indx_cp_r,cp_coeff_gpu,&
  &calorically_perfect,x_gpu,w_aux_gpu,fl_gpu,fhat_gpu)
    integer :: nx,ny,nz
    integer :: nv,ng,calorically_perfect
    integer :: indx_cp_l,indx_cp_r
    real(rkind) :: prandtl,t0
    real(rkind), dimension(:,:,:,:), target :: fl_gpu
    real(rkind), dimension(:,:,:,:), target :: w_aux_gpu
    real(rkind), dimension(:,:,:,:), target :: fhat_gpu
    real(rkind), dimension(:), target :: x_gpu
    real(rkind), dimension(:), target :: cp_coeff_gpu


    call visflx_x_kernel1_wrapper(c_null_ptr,nx,ny,nz,nv,ng,calorically_perfect,indx_cp_l,indx_cp_r,&
    &prandtl,t0,c_loc(fl_gpu),c_loc(w_aux_gpu),c_loc(fhat_gpu),c_loc(x_gpu),c_loc(cp_coeff_gpu))
    call hipCheck(hipDeviceSynchronize())

    call visflx_x_kernel2_wrapper(c_null_ptr,nx,ny,nz,nv,ng,calorically_perfect,indx_cp_l,indx_cp_r,&
    &prandtl,t0,c_loc(fl_gpu),c_loc(w_aux_gpu),c_loc(fhat_gpu),c_loc(x_gpu),c_loc(cp_coeff_gpu))
    call hipCheck(hipDeviceSynchronize())


  endsubroutine visflx_x_kernel


  subroutine visflx_y_kernel(nx,ny,nz,nv,ng,prandtl,t0,indx_cp_l,indx_cp_r,cp_coeff_gpu,calorically_perfect,y_gpu,w_aux_gpu,fl_gpu)
    integer :: nx,ny,nz
    integer :: nv,ng,calorically_perfect
    integer :: indx_cp_l,indx_cp_r
    real(rkind) :: prandtl,t0
    real(rkind), dimension(:,:,:,:), target :: fl_gpu
    real(rkind), dimension(:,:,:,:), target :: w_aux_gpu
    real(rkind), dimension(:), target :: y_gpu
    real(rkind), dimension(:), target :: cp_coeff_gpu


    call visflx_y_kernel_wrapper(c_null_ptr,nx,ny,nz,nv,ng,calorically_perfect,indx_cp_l,indx_cp_r,&
    &prandtl,t0,c_loc(fl_gpu),c_loc(w_aux_gpu),c_loc(y_gpu),c_loc(cp_coeff_gpu))
    call hipCheck(hipDeviceSynchronize())


  endsubroutine visflx_y_kernel

  subroutine visflx_z_kernel(nx,ny,nz,nv,ng,prandtl,t0,indx_cp_l,indx_cp_r,cp_coeff_gpu,calorically_perfect,z_gpu,w_aux_gpu,fl_gpu)
    integer :: nx,ny,nz
    integer :: nv,ng,calorically_perfect
    integer :: indx_cp_l,indx_cp_r
    real(rkind) :: prandtl,t0
    real(rkind), dimension(:,:,:,:), target :: fl_gpu
    real(rkind), dimension(:,:,:,:), target :: w_aux_gpu
    real(rkind), dimension(:), target :: z_gpu
    real(rkind), dimension(:), target :: cp_coeff_gpu


    call visflx_z_kernel_wrapper(c_null_ptr,nx,ny,nz,nv,ng,calorically_perfect,indx_cp_l,indx_cp_r,&
    &prandtl,t0,c_loc(fl_gpu),c_loc(w_aux_gpu),c_loc(z_gpu),c_loc(cp_coeff_gpu))
    call hipCheck(hipDeviceSynchronize())


  endsubroutine visflx_z_kernel


  subroutine recyc_exchange_kernel_1(irecyc,w_gpu,wbuf1s_gpu,nx,ny,nz,ng,nv)
    integer :: irecyc,nx,ny
    integer :: nz,ng,nv
    real(rkind), dimension(:,:,:,:), target :: w_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf1s_gpu


    call recyc_exchange_kernel_1_wrapper(c_null_ptr,irecyc,nx,ny,nz,ng,nv,c_loc(w_gpu),c_loc(wbuf1s_gpu))
    call hipCheck(hipDeviceSynchronize())


  endsubroutine recyc_exchange_kernel_1

  subroutine recyc_exchange_kernel_2(n1_start_recv,n1_start_send,n1_end_recv,wrecyc_gpu,wbuf1r_gpu,nx,ny,nz,ng,nv)
    integer :: n1_start_recv,n1_start_send,n1_end_recv
    integer :: nx,ny,nz
    integer :: ng,nv
    real(rkind), dimension(:,:,:,:), target :: wbuf1r_gpu
    real(rkind), dimension(:,:,:,:), target :: wrecyc_gpu


    call recyc_exchange_kernel_2_wrapper(c_null_ptr,n1_start_recv,n1_start_send,n1_end_recv,nx,ny,&
    &nz,ng,nv,c_loc(wbuf1r_gpu),c_loc(wrecyc_gpu))
    call hipCheck(hipDeviceSynchronize())


  endsubroutine recyc_exchange_kernel_2

  subroutine recyc_exchange_kernel_3(n2_start_recv,n2_start_send,n2_end_recv,wrecyc_gpu,wbuf2r_gpu,nx,ny,nz,ng,nv)
    integer :: n2_start_recv,n2_start_send,n2_end_recv
    integer :: nx,ny,nz
    integer :: ng,nv
    real(rkind), dimension(:,:,:,:), target :: wbuf2r_gpu
    real(rkind), dimension(:,:,:,:), target :: wrecyc_gpu


    call recyc_exchange_kernel_3_wrapper(c_null_ptr,n2_start_recv,n2_start_send,n2_end_recv,nx,ny,&
    &nz,ng,nv,c_loc(wbuf2r_gpu),c_loc(wrecyc_gpu))
    call hipCheck(hipDeviceSynchronize())


  endsubroutine recyc_exchange_kernel_3


  subroutine bcextr_sub_kernel(ilat,nx,ny,nz,ng,p0,rgas0,w_gpu,indx_cp_l,indx_cp_r,cv_coeff_gpu,t0,calorically_perfect)
    integer :: ilat,nx,ny
    integer :: nz,ng,indx_cp_l
    integer :: indx_cp_r,calorically_perfect
    real(rkind) :: p0,t0,rgas0
    real(rkind), dimension(:), target :: cv_coeff_gpu
    real(rkind), dimension(:,:,:,:), target :: w_gpu

    if (ilat==1) then
      call bcextr_sub_kernel1_wrapper(c_null_ptr,ilat,nx,ny,nz,ng,indx_cp_l,indx_cp_r,&
      &calorically_perfect,p0,t0,rgas0,c_loc(cv_coeff_gpu),c_loc(w_gpu))
      call hipCheck(hipDeviceSynchronize())
    elseif (ilat==2) then
      call bcextr_sub_kernel2_wrapper(c_null_ptr,ilat,nx,ny,nz,ng,indx_cp_l,indx_cp_r,&
      &calorically_perfect,p0,t0,rgas0,c_loc(cv_coeff_gpu),c_loc(w_gpu))
      call hipCheck(hipDeviceSynchronize())
    elseif (ilat==3) then
    elseif (ilat==4) then
      call bcextr_sub_kernel3_wrapper(c_null_ptr,ilat,nx,ny,nz,ng,indx_cp_l,indx_cp_r,&
      &calorically_perfect,p0,t0,rgas0,c_loc(cv_coeff_gpu),c_loc(w_gpu))
      call hipCheck(hipDeviceSynchronize())
    elseif (ilat==5) then
      call bcextr_sub_kernel4_wrapper(c_null_ptr,ilat,nx,ny,nz,ng,indx_cp_l,indx_cp_r,&
      &calorically_perfect,p0,t0,rgas0,c_loc(cv_coeff_gpu),c_loc(w_gpu))
      call hipCheck(hipDeviceSynchronize())
    elseif (ilat==6) then
      call bcextr_sub_kernel5_wrapper(c_null_ptr,ilat,nx,ny,nz,ng,indx_cp_l,indx_cp_r,&
      &calorically_perfect,p0,t0,rgas0,c_loc(cv_coeff_gpu),c_loc(w_gpu))
      call hipCheck(hipDeviceSynchronize())
    endif

  endsubroutine bcextr_sub_kernel








  subroutine bcrecyc_kernel_1(nx,ny,nz,ng,nv,wrecycav_gpu,wrecyc_gpu)
    integer :: nx,ny,nz
    integer :: ng,nv
    real(rkind), dimension(:,:,:), target :: wrecycav_gpu
    real(rkind), dimension(:,:,:,:), target :: wrecyc_gpu


    call bcrecyc_kernel_1_wrapper(c_null_ptr,nx,ny,nz,ng,nv,c_loc(wrecycav_gpu),c_loc(wrecyc_gpu))
    call hipCheck(hipDeviceSynchronize())


  endsubroutine bcrecyc_kernel_1


  subroutine bcrecyc_kernel_2(nx,ny,nz,nzmax,ng,wrecycav_gpu,wrecyc_gpu)
    integer :: nx,ny,nz
    integer :: nzmax,ng
    real(rkind), dimension(:,:,:), target :: wrecycav_gpu
    real(rkind), dimension(:,:,:,:), target :: wrecyc_gpu


    call bcrecyc_kernel_2_wrapper(c_null_ptr,nx,ny,nz,nzmax,ng,c_loc(wrecycav_gpu),c_loc(wrecyc_gpu))
    call hipCheck(hipDeviceSynchronize())


  endsubroutine bcrecyc_kernel_2



  subroutine bcrecyc_kernel_3(nx,ny,nz,ng,p0,u0,rgas0,w_gpu,wmean_gpu,wrecyc_gpu,weta_inflow_gpu,&
  &map_j_inn_gpu,map_j_out_gpu,map_j_out_blend_gpu,yplus_inflow_gpu,eta_inflow_gpu,yplus_recyc_gpu,&
  &eta_recyc_gpu,eta_recyc_blend_gpu,betarecyc,glund1,inflow_random_plane_gpu,indx_cp_l,indx_cp_r,&
  &cv_coeff_gpu,cp_coeff_gpu,t0,calorically_perfect,rand_type)
    integer :: nx,ny,nz
    integer :: ng,indx_cp_l,indx_cp_r
    integer :: calorically_perfect,rand_type
    real(rkind) :: p0,rgas0,betarecyc
    real(rkind) :: glund1,t0,u0
    real(rkind) :: u0_02
    real(rkind), dimension(:), target :: cv_coeff_gpu
    real(rkind), dimension(:), target :: cp_coeff_gpu
    real(rkind), dimension(:,:,:), target :: wmean_gpu
    real(rkind), dimension(:,:,:,:), target :: wrecyc_gpu
    real(rkind), dimension(:,:,:), target :: inflow_random_plane_gpu
    real(rkind), dimension(:,:,:,:), target :: w_gpu
    real(rkind), dimension(:), target :: weta_inflow_gpu
    real(rkind), dimension(:), target :: yplus_inflow_gpu
    real(rkind), dimension(:), target :: eta_inflow_gpu
    real(rkind), dimension(:), target :: yplus_recyc_gpu
    real(rkind), dimension(:), target :: eta_recyc_gpu
    real(rkind), dimension(:), target :: eta_recyc_blend_gpu
    integer, dimension(:), target :: map_j_inn_gpu
    integer, dimension(:), target :: map_j_out_gpu
    integer, dimension(:), target :: map_j_out_blend_gpu

    u0_02 = 0.02_rkind*u0
    if (rand_type==0) u0_02 = 0._rkind
    call bcrecyc_kernel_3_wrapper(c_null_ptr,nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect,&
    &rand_type,p0,rgas0,betarecyc,glund1,t0,u0,u0_02,c_loc(map_j_inn_gpu),c_loc(map_j_out_gpu),&
    &c_loc(map_j_out_blend_gpu),c_loc(cv_coeff_gpu),c_loc(cp_coeff_gpu),c_loc(wmean_gpu),&
    &c_loc(wrecyc_gpu),c_loc(inflow_random_plane_gpu),c_loc(w_gpu),c_loc(weta_inflow_gpu),&
    &c_loc(yplus_inflow_gpu),c_loc(eta_inflow_gpu),c_loc(yplus_recyc_gpu),c_loc(eta_recyc_gpu),&
    &c_loc(eta_recyc_blend_gpu))
    call hipCheck(hipDeviceSynchronize())


  endsubroutine bcrecyc_kernel_3


  subroutine bclam_kernel(ilat,nx,ny,nz,ng,nv,w_gpu,wmean_gpu,p0,rgas0,indx_cp_l,indx_cp_r,cv_coeff_gpu,t0,calorically_perfect)
    integer :: nx,ny,nz
    integer :: ng,nv,indx_cp_l
    integer :: indx_cp_r,calorically_perfect,ilat
    real(rkind) :: p0,rgas0,t0
    real(rkind), dimension(:), target :: cv_coeff_gpu
    real(rkind), dimension(:,:,:,:), target :: w_gpu
    real(rkind), dimension(:,:,:), target :: wmean_gpu

    if (ilat==1) then
      call bclam_kernel_wrapper(c_null_ptr,nx,ny,nz,ng,nv,indx_cp_l,indx_cp_r,calorically_perfect,&
      &ilat,p0,rgas0,t0,c_loc(cv_coeff_gpu),c_loc(w_gpu),c_loc(wmean_gpu))
      call hipCheck(hipDeviceSynchronize())
    elseif (ilat==2) then
    elseif (ilat==3) then
    elseif (ilat==4) then
    elseif (ilat==5) then
    elseif (ilat==6) then
    endif

  endsubroutine bclam_kernel


  subroutine bcfree_kernel(ilat,nx,ny,nz,ng,nv,winf_gpu,w_gpu)
    integer :: ilat,nx,ny
    integer :: nz,ng,nv
    real(rkind), dimension(:), target :: winf_gpu
    real(rkind), dimension(:,:,:,:), target :: w_gpu

    if (ilat==1) then
      call bcfree_kernel1_wrapper(c_null_ptr,ilat,nx,ny,nz,ng,nv,c_loc(winf_gpu),c_loc(w_gpu))
      call hipCheck(hipDeviceSynchronize())
    elseif (ilat==2) then
      call bcfree_kernel2_wrapper(c_null_ptr,ilat,nx,ny,nz,ng,nv,c_loc(winf_gpu),c_loc(w_gpu))
      call hipCheck(hipDeviceSynchronize())
    elseif (ilat==3) then
      call bcfree_kernel3_wrapper(c_null_ptr,ilat,nx,ny,nz,ng,nv,c_loc(winf_gpu),c_loc(w_gpu))
      call hipCheck(hipDeviceSynchronize())
    elseif (ilat==4) then
      call bcfree_kernel4_wrapper(c_null_ptr,ilat,nx,ny,nz,ng,nv,c_loc(winf_gpu),c_loc(w_gpu))
      call hipCheck(hipDeviceSynchronize())
    elseif (ilat==5) then
      call bcfree_kernel5_wrapper(c_null_ptr,ilat,nx,ny,nz,ng,nv,c_loc(winf_gpu),c_loc(w_gpu))
      call hipCheck(hipDeviceSynchronize())
    elseif (ilat==6) then
      call bcfree_kernel6_wrapper(c_null_ptr,ilat,nx,ny,nz,ng,nv,c_loc(winf_gpu),c_loc(w_gpu))
      call hipCheck(hipDeviceSynchronize())
    endif

  endsubroutine bcfree_kernel


  subroutine bcfree_sub_kernel(ilat,nx,ny,nz,ng,w_gpu,w_aux_gpu,aoa,t0,ptot0,ttot0,rgas0,indx_cp_l,&
  &indx_cp_r,cv_coeff_gpu,cp_coeff_gpu,calorically_perfect,tol_iter_nr)
    integer :: ilat,nx,ny
    integer :: nz,ng,indx_cp_l
    integer :: indx_cp_r,calorically_perfect
    real(rkind) :: ptot0,ttot0,rgas0
    real(rkind) :: aoa,t0,tol_iter_nr
    real(rkind) :: cosangle,sinangle
    real(rkind), dimension(:,:,:,:), target :: w_gpu
    real(rkind), dimension(:,:,:,:), target :: w_aux_gpu
    real(rkind), dimension(:), target :: cv_coeff_gpu
    real(rkind), dimension(:), target :: cp_coeff_gpu

    cosangle = cos(aoa)
    sinangle = sin(aoa)

    if (ilat==1) then
      call bcfree_sub_kernel1_wrapper(c_null_ptr,ilat,nx,ny,nz,ng,indx_cp_l,indx_cp_r,&
      &calorically_perfect,ptot0,ttot0,rgas0,aoa,t0,tol_iter_nr,cosangle,sinangle,c_loc(w_gpu),&
      &c_loc(w_aux_gpu),c_loc(cv_coeff_gpu),c_loc(cp_coeff_gpu))
      call hipCheck(hipDeviceSynchronize())
    elseif (ilat==4) then
      call bcfree_sub_kernel2_wrapper(c_null_ptr,ilat,nx,ny,nz,ng,indx_cp_l,indx_cp_r,&
      &calorically_perfect,ptot0,ttot0,rgas0,aoa,t0,tol_iter_nr,cosangle,sinangle,c_loc(w_gpu),&
      &c_loc(w_aux_gpu),c_loc(cv_coeff_gpu),c_loc(cp_coeff_gpu))
      call hipCheck(hipDeviceSynchronize())
    endif

  endsubroutine bcfree_sub_kernel


  subroutine bcshock_kernel(ilat,nx,ny,nz,ng,nv,w_gpu,winf_gpu,winf_past_shock_gpu,xshock_imp,shock_angle,x_gpu,y_gpu,tanhfacs)
    integer :: nx,ny,nz
    integer :: ng,nv,ilat
    real(rkind) :: xshock_imp,shock_angle,tanhfacs
    real(rkind) :: tanhlen
    real(rkind), dimension(:), target :: winf_gpu
    real(rkind), dimension(:), target :: winf_past_shock_gpu
    real(rkind), dimension(:,:,:,:), target :: w_gpu
    real(rkind), dimension(:), target :: x_gpu
    real(rkind), dimension(:), target :: y_gpu

    tanhlen = 8._rkind*tanhfacs

    if (ilat==1) then
    elseif (ilat==2) then
    elseif (ilat==3) then
    elseif (ilat==4) then
      call bcshock_kernel_wrapper(c_null_ptr,nx,ny,nz,ng,nv,ilat,xshock_imp,shock_angle,tanhfacs,&
      &tanhlen,c_loc(winf_gpu),c_loc(winf_past_shock_gpu),c_loc(w_gpu),c_loc(x_gpu),c_loc(y_gpu))
      call hipCheck(hipDeviceSynchronize())
    elseif (ilat==5) then
    elseif (ilat==6) then
    endif

  endsubroutine bcshock_kernel

  subroutine bcextr_var_kernel(nx,ny,nz,ng,w_gpu)
    integer :: nx,ny,nz
    integer :: ng
    real(rkind), dimension(:,:,:,:), target :: w_gpu


    call bcextr_var_kernel1_wrapper(c_null_ptr,nx,ny,nz,ng,c_loc(w_gpu))
    call hipCheck(hipDeviceSynchronize())

    call bcextr_var_kernel2_wrapper(c_null_ptr,nx,ny,nz,ng,c_loc(w_gpu))
    call hipCheck(hipDeviceSynchronize())

    call bcextr_var_kernel3_wrapper(c_null_ptr,nx,ny,nz,ng,c_loc(w_gpu))
    call hipCheck(hipDeviceSynchronize())

    call bcextr_var_kernel4_wrapper(c_null_ptr,nx,ny,nz,ng,c_loc(w_gpu))
    call hipCheck(hipDeviceSynchronize())

    call bcextr_var_kernel5_wrapper(c_null_ptr,nx,ny,nz,ng,c_loc(w_gpu))
    call hipCheck(hipDeviceSynchronize())

    call bcextr_var_kernel6_wrapper(c_null_ptr,nx,ny,nz,ng,c_loc(w_gpu))
    call hipCheck(hipDeviceSynchronize())


  endsubroutine bcextr_var_kernel


  subroutine bcextr_airfoil_var_kernel(nx,ny,nz,ng,ndim,wall_tag_gpu,w_gpu,ileftx,irightx,ileftz,irightz)
    integer :: nx,ny,nz
    integer :: ng,ndim,ileftx
    integer :: irightx,ileftz,irightz
    real(rkind), dimension(:,:,:,:), target :: w_gpu
    integer, dimension(:), target :: wall_tag_gpu

    if (ileftx < 0) then
      call bcextr_airfoil_var_kernel1_wrapper(c_null_ptr,nx,ny,nz,ng,ndim,ileftx,irightx,ileftz,&
      &irightz,c_loc(wall_tag_gpu),c_loc(w_gpu))
      call hipCheck(hipDeviceSynchronize())
    endif
    if (irightx < 0) then
      call bcextr_airfoil_var_kernel2_wrapper(c_null_ptr,nx,ny,nz,ng,ndim,ileftx,irightx,ileftz,&
      &irightz,c_loc(wall_tag_gpu),c_loc(w_gpu))
      call hipCheck(hipDeviceSynchronize())
    endif
    call bcextr_airfoil_var_kernel3_wrapper(c_null_ptr,nx,ny,nz,ng,ndim,ileftx,irightx,ileftz,&
    &irightz,c_loc(wall_tag_gpu),c_loc(w_gpu))
    call hipCheck(hipDeviceSynchronize())

    call bcextr_airfoil_var_kernel4_wrapper(c_null_ptr,nx,ny,nz,ng,ndim,ileftx,irightx,ileftz,&
    &irightz,c_loc(wall_tag_gpu),c_loc(w_gpu))
    call hipCheck(hipDeviceSynchronize())
    if(ndim == 3) then
      if (ileftz < 0) then
        call bcextr_airfoil_var_kernel5_wrapper(c_null_ptr,nx,ny,nz,ng,ndim,ileftx,irightx,ileftz,&
        &irightz,c_loc(wall_tag_gpu),c_loc(w_gpu))
        call hipCheck(hipDeviceSynchronize())
      endif
      if (irightz < 0) then
        call bcextr_airfoil_var_kernel6_wrapper(c_null_ptr,nx,ny,nz,ng,ndim,ileftx,irightx,ileftz,&
        &irightz,c_loc(wall_tag_gpu),c_loc(w_gpu))
        call hipCheck(hipDeviceSynchronize())
      endif
    endif

  endsubroutine bcextr_airfoil_var_kernel

  subroutine bcextr_kernel(ilat,nx,ny,nz,ng,nv,w_gpu)
    integer :: nx,ny,nz
    integer :: ng,nv,ilat
    real(rkind), dimension(:,:,:,:), target :: w_gpu

    if (ilat==1) then
      call bcextr_kernel1_wrapper(c_null_ptr,nx,ny,nz,ng,nv,ilat,c_loc(w_gpu))
      call hipCheck(hipDeviceSynchronize())
    elseif (ilat==2) then
      call bcextr_kernel2_wrapper(c_null_ptr,nx,ny,nz,ng,nv,ilat,c_loc(w_gpu))
      call hipCheck(hipDeviceSynchronize())
    elseif (ilat==3) then
      call bcextr_kernel3_wrapper(c_null_ptr,nx,ny,nz,ng,nv,ilat,c_loc(w_gpu))
      call hipCheck(hipDeviceSynchronize())
    elseif (ilat==4) then
      call bcextr_kernel4_wrapper(c_null_ptr,nx,ny,nz,ng,nv,ilat,c_loc(w_gpu))
      call hipCheck(hipDeviceSynchronize())
    elseif (ilat==5) then
      call bcextr_kernel5_wrapper(c_null_ptr,nx,ny,nz,ng,nv,ilat,c_loc(w_gpu))
      call hipCheck(hipDeviceSynchronize())
    elseif (ilat==6) then
      call bcextr_kernel6_wrapper(c_null_ptr,nx,ny,nz,ng,nv,ilat,c_loc(w_gpu))
      call hipCheck(hipDeviceSynchronize())
    endif

  endsubroutine bcextr_kernel


  subroutine bcsym_kernel(ilat,nx,ny,nz,ng,w_gpu)
    integer :: nx,ny,nz
    integer :: ng,ilat
    real(rkind), dimension(:,:,:,:), target :: w_gpu

    if (ilat==1) then
    elseif (ilat==2) then
    elseif (ilat==3) then
      call bcsym_kernel_wrapper(c_null_ptr,nx,ny,nz,ng,ilat,c_loc(w_gpu))
      call hipCheck(hipDeviceSynchronize())
    elseif (ilat==4) then
    elseif (ilat==5) then
    elseif (ilat==6) then
    endif

  endsubroutine bcsym_kernel


  subroutine bcsym_c2_kernel(ilat,nx,ny,nz,ng,w_gpu,dxdcsic2_gpu,dydcsic2_gpu)
    integer :: nx,ny,nz
    integer :: ng,ilat
    real(rkind), dimension(:,:,:,:), target :: w_gpu
    real(rkind), dimension(:,:), target :: dxdcsic2_gpu
    real(rkind), dimension(:,:), target :: dydcsic2_gpu

    if (ilat==1) then
    elseif (ilat==2) then
    elseif (ilat==3) then
      call bcsym_c2_kernel1_wrapper(c_null_ptr,nx,ny,nz,ng,ilat,c_loc(w_gpu),c_loc(dxdcsic2_gpu),c_loc(dydcsic2_gpu))
      call hipCheck(hipDeviceSynchronize())
    elseif (ilat==4) then
      call bcsym_c2_kernel2_wrapper(c_null_ptr,nx,ny,nz,ng,ilat,c_loc(w_gpu),c_loc(dxdcsic2_gpu),c_loc(dydcsic2_gpu))
      call hipCheck(hipDeviceSynchronize())
    elseif (ilat==5) then
    elseif (ilat==6) then
    endif

  endsubroutine bcsym_c2_kernel


  subroutine bcwall_kernel(ilat,nx,ny,nz,ng,twall,w_gpu,w_aux_gpu,indx_cp_l,indx_cp_r,cv_coeff_gpu,&
  &t0,rgas0,calorically_perfect,tol_iter_nr)
    integer :: nx,ny,nz
    integer :: ng,indx_cp_l,indx_cp_r
    integer :: calorically_perfect,ilat
    real(rkind) :: twall,t0,rgas0
    real(rkind) :: tol_iter_nr
    real(rkind), dimension(:,:,:,:), target :: w_gpu
    real(rkind), dimension(:,:,:,:), target :: w_aux_gpu
    real(rkind), dimension(:), target :: cv_coeff_gpu

    if (ilat==1) then
    elseif (ilat==2) then
    elseif (ilat==3) then
      call bcwall_kernel_wrapper(c_null_ptr,nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect,&
      &ilat,twall,t0,rgas0,tol_iter_nr,c_loc(w_gpu),c_loc(w_aux_gpu),c_loc(cv_coeff_gpu))
      call hipCheck(hipDeviceSynchronize())
    elseif (ilat==4) then
    elseif (ilat==5) then
    elseif (ilat==6) then
    endif

  endsubroutine bcwall_kernel


  subroutine bcwall_airfoil_kernel(ilat,nx,ny,nz,ng,twall,w_gpu,w_aux_gpu,wall_tag_gpu,indx_cp_l,&
  &indx_cp_r,cv_coeff_gpu,t0,rgas0,calorically_perfect,tol_iter_nr,xc2_gpu,yc2_gpu,a_tw,v_bs,thic,&
  &kx_tw,om_tw,xtw1,xtw2,time,dxdetanc2_gpu,dydetanc2_gpu,u0)
    integer :: nx,ny,nz
    integer :: ng,indx_cp_l,indx_cp_r
    integer :: calorically_perfect,ilat
    real(rkind) :: u0,twall,t0
    real(rkind) :: rgas0,tol_iter_nr,a_tw
    real(rkind) :: v_bs,thic,kx_tw
    real(rkind) :: om_tw,time,xtw1
    real(rkind) :: xtw2
    real(rkind), dimension(:,:,:,:), target :: w_gpu
    real(rkind), dimension(:,:,:,:), target :: w_aux_gpu
    real(rkind), dimension(:), target :: cv_coeff_gpu
    real(rkind), dimension(:,:), target :: xc2_gpu
    real(rkind), dimension(:,:), target :: yc2_gpu
    real(rkind), dimension(:,:), target :: dxdetanc2_gpu
    real(rkind), dimension(:,:), target :: dydetanc2_gpu
    integer, dimension(:), target :: wall_tag_gpu

    if (ilat==1) then
    elseif (ilat==2) then
    elseif (ilat==3) then
      call bcwall_airfoil_kernel_wrapper(c_null_ptr,nx,ny,nz,ng,indx_cp_l,indx_cp_r,&
      &calorically_perfect,ilat,u0,twall,t0,rgas0,tol_iter_nr,a_tw,v_bs,thic,kx_tw,om_tw,time,xtw1,xtw2,&
      &c_loc(wall_tag_gpu),c_loc(w_gpu),c_loc(w_aux_gpu),c_loc(cv_coeff_gpu),c_loc(xc2_gpu),c_loc(yc2_gpu),&
      &c_loc(dxdetanc2_gpu),c_loc(dydetanc2_gpu))
      call hipCheck(hipDeviceSynchronize())
    elseif (ilat==4) then
    elseif (ilat==5) then
    elseif (ilat==6) then
    endif

  endsubroutine bcwall_airfoil_kernel


  subroutine bcwall_staggered_kernel(ilat,nx,ny,nz,ng,twall,w_gpu,w_aux_gpu,indx_cp_l,indx_cp_r,&
  &cv_coeff_gpu,t0,rgas0,calorically_perfect,tol_iter_nr)
    integer :: nx,ny,nz
    integer :: ng,indx_cp_l,indx_cp_r
    integer :: calorically_perfect,ilat
    real(rkind) :: twall,t0,rgas0
    real(rkind) :: tol_iter_nr
    real(rkind), dimension(:,:,:,:), target :: w_gpu
    real(rkind), dimension(:,:,:,:), target :: w_aux_gpu
    real(rkind), dimension(:), target :: cv_coeff_gpu

    if (ilat==1) then
    elseif (ilat==2) then
    elseif (ilat==3) then
      call bcwall_staggered_kernel1_wrapper(c_null_ptr,nx,ny,nz,ng,indx_cp_l,indx_cp_r,&
      &calorically_perfect,ilat,twall,t0,rgas0,tol_iter_nr,c_loc(w_gpu),c_loc(w_aux_gpu),&
      &c_loc(cv_coeff_gpu))
      call hipCheck(hipDeviceSynchronize())
    elseif (ilat==4) then
      call bcwall_staggered_kernel2_wrapper(c_null_ptr,nx,ny,nz,ng,indx_cp_l,indx_cp_r,&
      &calorically_perfect,ilat,twall,t0,rgas0,tol_iter_nr,c_loc(w_gpu),c_loc(w_aux_gpu),&
      &c_loc(cv_coeff_gpu))
      call hipCheck(hipDeviceSynchronize())
    elseif (ilat==5) then
    elseif (ilat==6) then
    endif

  endsubroutine bcwall_staggered_kernel


  subroutine compute_residual_kernel(nx,ny,nz,ng,nv,fln_gpu,dt,residual_rhou,fluid_mask_gpu,redn_3d_gpu)
    integer :: nx,ny,nz
    integer :: ng,nv
    real(rkind) :: dt
    real(rkind) :: residual_rhou
    real(rkind), dimension(:,:,:,:), target :: fln_gpu
    integer, dimension(:,:,:), target :: fluid_mask_gpu
    real(rkind), dimension(:,:,:), target :: redn_3d_gpu

    residual_rhou = 0._rkind
    call compute_residual_kernel_wrapper(c_null_ptr,nx,ny,nz,ng,nv,dt,c_loc(fluid_mask_gpu),&
    &c_loc(fln_gpu),residual_rhou,c_loc(redn_3d_gpu))
    call hipCheck(hipDeviceSynchronize())


  endsubroutine compute_residual_kernel


  subroutine compute_airfoil_forces_runtime_kernel(nx,ny,nz,ng,nv,p0,u0,rgas0,w_aux_gpu,meta_gpu,&
  &csimod_gpu,wall_tag_gpu,dxdcsic2_gpu,dydcsic2_gpu,n,a,pn,pa,tn,ta,redn_3d_gpu)
    integer :: nx,ny,nz
    integer :: ng,nv
    real(rkind) :: p0,u0,rgas0
    real(rkind) :: n,a,pn
    real(rkind) :: pa,tn,ta
    real(rkind), dimension(:,:,:,:), target :: w_aux_gpu
    real(rkind), dimension(:,:), target :: meta_gpu
    real(rkind), dimension(:,:), target :: csimod_gpu
    real(rkind), dimension(:,:), target :: dxdcsic2_gpu
    real(rkind), dimension(:,:), target :: dydcsic2_gpu
    integer, dimension(:), target :: wall_tag_gpu
    real(rkind), dimension(:,:,:), target :: redn_3d_gpu

    n = 0._rkind
    a = 0._rkind
    pn = 0._rkind
    pa = 0._rkind
    tn = 0._rkind
    ta = 0._rkind
    call compute_airfoil_forces_runtime_kernel_wrapper(c_null_ptr,nx,ny,nz,ng,nv,p0,u0,rgas0,&
    &c_loc(wall_tag_gpu),c_loc(w_aux_gpu),c_loc(meta_gpu),c_loc(csimod_gpu),c_loc(dxdcsic2_gpu),&
    &c_loc(dydcsic2_gpu),n,a,pn,pa,tn,ta,c_loc(redn_3d_gpu))
    call hipCheck(hipDeviceSynchronize())


  endsubroutine compute_airfoil_forces_runtime_kernel


  subroutine compute_rho_t_p_minmax_kernel(nx,ny,nz,ng,rgas0,w_aux_gpu,rhomin,rhomax,tmin,tmax,pmin,&
  &pmax,fluid_mask_gpu,redn_3d_gpu)
    integer :: nx,ny,nz
    integer :: ng
    real(rkind) :: rgas0
    real(rkind) :: rhomin,tmin,pmin
    real(rkind) :: rhomax,tmax,pmax
    real(rkind), dimension(:,:,:,:), target :: w_aux_gpu
    integer, dimension(:,:,:), target :: fluid_mask_gpu
    real(rkind), dimension(:,:,:), target :: redn_3d_gpu

    rhomin = huge(1._rkind)
    rhomax = -100._rkind
    tmin = huge(1._rkind)
    tmax = -100._rkind
    pmin = huge(1._rkind)
    pmax = -100._rkind
    call compute_rho_t_p_minmax_kernel1_wrapper(c_null_ptr,nx,ny,nz,ng,rgas0,c_loc(fluid_mask_gpu),&
    &c_loc(w_aux_gpu),rhomin,tmin,pmin,c_loc(redn_3d_gpu))
    call hipCheck(hipDeviceSynchronize())

    call compute_rho_t_p_minmax_kernel2_wrapper(c_null_ptr,nx,ny,nz,ng,rgas0,c_loc(fluid_mask_gpu),&
    &c_loc(w_aux_gpu),rhomax,tmax,pmax,c_loc(redn_3d_gpu))
    call hipCheck(hipDeviceSynchronize())


  endsubroutine compute_rho_t_p_minmax_kernel


  subroutine compute_dt_kernel(nx,ny,nz,ng,rgas0,prandtl,dcsidx_gpu,detady_gpu,dzitdz_gpu,&
  &dcsidxs_gpu,detadys_gpu,dzitdzs_gpu,w_gpu,w_aux_gpu,dtxi_max,dtyi_max,dtzi_max,dtxv_max,dtyv_max,&
  &dtzv_max,dtxk_max,dtyk_max,dtzk_max,indx_cp_l,indx_cp_r,cp_coeff_gpu,fluid_mask_gpu,&
  &calorically_perfect,t0,redn_3d_gpu)
    integer :: nx,ny,nz
    integer :: ng,indx_cp_l,indx_cp_r
    integer :: calorically_perfect
    real(rkind) :: rgas0,t0,prandtl
    real(rkind) :: dtxi_max,dtyi_max,dtzi_max
    real(rkind) :: dtxv_max,dtyv_max,dtzv_max
    real(rkind) :: dtxk_max,dtyk_max,dtzk_max
    real(rkind), dimension(:,:,:,:), target :: w_gpu
    real(rkind), dimension(:,:,:,:), target :: w_aux_gpu
    real(rkind), dimension(:), target :: cp_coeff_gpu
    real(rkind), dimension(:), target :: dcsidx_gpu
    real(rkind), dimension(:), target :: detady_gpu
    real(rkind), dimension(:), target :: dzitdz_gpu
    real(rkind), dimension(:), target :: dcsidxs_gpu
    real(rkind), dimension(:), target :: detadys_gpu
    real(rkind), dimension(:), target :: dzitdzs_gpu
    integer, dimension(:,:,:), target :: fluid_mask_gpu
    real(rkind), dimension(:,:,:), target :: redn_3d_gpu

    dtxi_max = 0._rkind
    dtyi_max = 0._rkind
    dtzi_max = 0._rkind
    dtxv_max = 0._rkind
    dtyv_max = 0._rkind
    dtzv_max = 0._rkind
    dtxk_max = 0._rkind
    dtyk_max = 0._rkind
    dtzk_max = 0._rkind
    call compute_dt_kernel_wrapper(c_null_ptr,nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect,&
    &rgas0,t0,prandtl,c_loc(fluid_mask_gpu),c_loc(w_gpu),c_loc(w_aux_gpu),c_loc(cp_coeff_gpu),&
    &c_loc(dcsidx_gpu),c_loc(detady_gpu),c_loc(dzitdz_gpu),c_loc(dcsidxs_gpu),c_loc(detadys_gpu),&
    &c_loc(dzitdzs_gpu),dtxi_max,dtyi_max,dtzi_max,dtxv_max,dtyv_max,dtzv_max,dtxk_max,dtyk_max,dtzk_max,&
    &c_loc(redn_3d_gpu))
    call hipCheck(hipDeviceSynchronize())


  endsubroutine compute_dt_kernel


  subroutine compute_dt_c2_kernel(nx,ny,nz,ng,rgas0,prandtl,dcsidxnc2_gpu,dcsidync2_gpu,&
  &detadxnc2_gpu,detadync2_gpu,mcsi_gpu,meta_gpu,dzitdz_gpu,dzitdzs_gpu,w_gpu,w_aux_gpu,dtxi_max,&
  &dtyi_max,dtzi_max,dtxv_max,dtyv_max,dtzv_max,dtxk_max,dtyk_max,dtzk_max,indx_cp_l,indx_cp_r,&
  &cp_coeff_gpu,fluid_mask_gpu,calorically_perfect,t0,redn_3d_gpu)
    integer :: nx,ny,nz
    integer :: ng,indx_cp_l,indx_cp_r
    integer :: calorically_perfect
    real(rkind) :: rgas0,t0,prandtl
    real(rkind) :: dtxi_max,dtyi_max,dtzi_max
    real(rkind) :: dtxv_max,dtyv_max,dtzv_max
    real(rkind) :: dtxk_max,dtyk_max,dtzk_max
    real(rkind), dimension(:,:,:,:), target :: w_gpu
    real(rkind), dimension(:,:,:,:), target :: w_aux_gpu
    real(rkind), dimension(:), target :: cp_coeff_gpu
    real(rkind), dimension(:), target :: dzitdz_gpu
    real(rkind), dimension(:), target :: dzitdzs_gpu
    real(rkind), dimension(:,:), target :: dcsidxnc2_gpu
    real(rkind), dimension(:,:), target :: dcsidync2_gpu
    real(rkind), dimension(:,:), target :: detadxnc2_gpu
    real(rkind), dimension(:,:), target :: detadync2_gpu
    real(rkind), dimension(:,:), target :: mcsi_gpu
    real(rkind), dimension(:,:), target :: meta_gpu
    integer, dimension(:,:,:), target :: fluid_mask_gpu
    real(rkind), dimension(:,:,:), target :: redn_3d_gpu

    dtxi_max = 0._rkind
    dtyi_max = 0._rkind
    dtzi_max = 0._rkind
    dtxv_max = 0._rkind
    dtyv_max = 0._rkind
    dtzv_max = 0._rkind
    dtxk_max = 0._rkind
    dtyk_max = 0._rkind
    dtzk_max = 0._rkind
    call compute_dt_c2_kernel_wrapper(c_null_ptr,nx,ny,nz,ng,indx_cp_l,indx_cp_r,&
    &calorically_perfect,rgas0,t0,prandtl,c_loc(fluid_mask_gpu),c_loc(w_gpu),c_loc(w_aux_gpu),&
    &c_loc(cp_coeff_gpu),c_loc(dzitdz_gpu),c_loc(dzitdzs_gpu),c_loc(dcsidxnc2_gpu),c_loc(dcsidync2_gpu),&
    &c_loc(detadxnc2_gpu),c_loc(detadync2_gpu),c_loc(mcsi_gpu),c_loc(meta_gpu),dtxi_max,dtyi_max,&
    &dtzi_max,dtxv_max,dtyv_max,dtzv_max,dtxk_max,dtyk_max,dtzk_max,c_loc(redn_3d_gpu))
    call hipCheck(hipDeviceSynchronize())


  endsubroutine compute_dt_c2_kernel


  subroutine eval_aux_kernel(nx,ny,nz,ng,istart,iend,jstart,jend,kstart,kend,w_gpu,w_aux_gpu,&
  &visc_model,mu0,t0,sutherland_s,t_ref_dim,powerlaw_vtexp,visc_power,visc_sutherland,visc_no,&
  &cv_coeff_gpu,indx_cp_l,indx_cp_r,rgas0,calorically_perfect,tol_iter_nr,stream_id)
    integer :: nx,ny,nz
    integer :: ng,visc_model,istart
    integer :: iend,jstart,jend
    integer :: kstart,kend,visc_power
    integer :: visc_sutherland,visc_no,indx_cp_l
    integer :: indx_cp_r,calorically_perfect
    real(rkind) :: mu0,t0,sutherland_s
    real(rkind) :: t_ref_dim,powerlaw_vtexp,rgas0
    real(rkind) :: tol_iter_nr
    real(rkind), dimension(:), target :: cv_coeff_gpu
    real(rkind), dimension(:,:,:,:), target :: w_gpu
    real(rkind), dimension(:,:,:,:), target :: w_aux_gpu

    type(c_ptr), value :: stream_id

    call eval_aux_kernel_wrapper(stream_id,nx,ny,nz,ng,visc_model,istart,iend,jstart,jend,kstart,&
    &kend,visc_power,visc_sutherland,visc_no,indx_cp_l,indx_cp_r,calorically_perfect,mu0,t0,sutherland_s,&
    &t_ref_dim,powerlaw_vtexp,rgas0,tol_iter_nr,c_loc(cv_coeff_gpu),c_loc(w_gpu),c_loc(w_aux_gpu))


  endsubroutine eval_aux_kernel


  subroutine tripping_pressure_kernel(nx,ny,nz,ng,nv,i_rank_start,pi,itr1,itr2,x0tr,y0tr,x0ts,y0ts,&
  &lamx,lamy,lamz,lamz1,phiz,phiz1,asl,bt,xc2_gpu,yc2_gpu,z_gpu,w_gpu,fl_gpu,wall_tag_gpu,&
  &dxdetanc2_gpu,dydetanc2_gpu)
    integer :: nx,ny,nz
    integer :: nv,ng,itr1
    integer :: itr2,i_rank_start
    real(rkind) :: pi,x0tr,y0tr
    real(rkind) :: x0ts,y0ts,lamx
    real(rkind) :: lamy,lamz,lamz1
    real(rkind) :: phiz,phiz1,asl
    real(rkind) :: bt
    real(rkind), dimension(:,:,:,:), target :: w_gpu
    real(rkind), dimension(:,:,:,:), target :: fl_gpu
    real(rkind), dimension(:,:), target :: dxdetanc2_gpu
    real(rkind), dimension(:,:), target :: dydetanc2_gpu
    real(rkind), dimension(:,:), target :: xc2_gpu
    real(rkind), dimension(:,:), target :: yc2_gpu
    real(rkind), dimension(:), target :: z_gpu
    integer, dimension(:), target :: wall_tag_gpu


    call tripping_pressure_kernel_wrapper(c_null_ptr,nx,ny,nz,nv,ng,itr1,itr2,i_rank_start,pi,x0tr,&
    &y0tr,x0ts,y0ts,lamx,lamy,lamz,lamz1,phiz,phiz1,asl,bt,c_loc(wall_tag_gpu),c_loc(w_gpu),&
    &c_loc(fl_gpu),c_loc(dxdetanc2_gpu),c_loc(dydetanc2_gpu),c_loc(xc2_gpu),c_loc(yc2_gpu),c_loc(z_gpu))
    call hipCheck(hipDeviceSynchronize())


  endsubroutine tripping_pressure_kernel

  subroutine tripping_suction_kernel(nx,ny,nz,ng,nv,i_rank_start,pi,its1,its2,x0tr,y0tr,x0ts,y0ts,&
  &lamx,lamy,lams,lams1,phis,phis1,asl,bt,xc2_gpu,yc2_gpu,z_gpu,w_gpu,fl_gpu,wall_tag_gpu,&
  &dxdetanc2_gpu,dydetanc2_gpu)
    integer :: nx,ny,nz
    integer :: nv,ng,its1
    integer :: its2,i_rank_start
    real(rkind) :: pi,x0tr,y0tr
    real(rkind) :: x0ts,y0ts,lamx
    real(rkind) :: lamy,lams,lams1
    real(rkind) :: phis,phis1,asl
    real(rkind) :: bt
    real(rkind), dimension(:,:,:,:), target :: w_gpu
    real(rkind), dimension(:,:,:,:), target :: fl_gpu
    real(rkind), dimension(:,:), target :: dxdetanc2_gpu
    real(rkind), dimension(:,:), target :: dydetanc2_gpu
    real(rkind), dimension(:,:), target :: xc2_gpu
    real(rkind), dimension(:,:), target :: yc2_gpu
    real(rkind), dimension(:), target :: z_gpu
    integer, dimension(:), target :: wall_tag_gpu


    call tripping_suction_kernel_wrapper(c_null_ptr,nx,ny,nz,nv,ng,its1,its2,i_rank_start,pi,x0tr,&
    &y0tr,x0ts,y0ts,lamx,lamy,lams,lams1,phis,phis1,asl,bt,c_loc(wall_tag_gpu),c_loc(w_gpu),&
    &c_loc(fl_gpu),c_loc(dxdetanc2_gpu),c_loc(dydetanc2_gpu),c_loc(xc2_gpu),c_loc(yc2_gpu),c_loc(z_gpu))
    call hipCheck(hipDeviceSynchronize())


  endsubroutine tripping_suction_kernel














  subroutine insitu_div_kernel(nx,ny,nz,ng,visc_order,npsi,mpsi,w_aux_gpu,coeff_deriv1_gpu,&
  &dcsidx_gpu,detady_gpu,dzitdz_gpu,psi_gpu)
    integer :: nx,ny,nz
    integer :: ng,visc_order,npsi
    integer :: mpsi,lmax
    real(rkind), dimension(:,:,:,:), target :: w_aux_gpu
    real(rkind), dimension(:,:), target :: coeff_deriv1_gpu
    real(rkind), dimension(:), target :: dcsidx_gpu
    real(rkind), dimension(:), target :: detady_gpu
    real(rkind), dimension(:), target :: dzitdz_gpu
    real(rkind), dimension(:,:,:,:), target :: psi_gpu

    lmax = visc_order / 2
    call insitu_div_kernel_wrapper(c_null_ptr,nx,ny,nz,ng,visc_order,npsi,mpsi,lmax,&
    &c_loc(w_aux_gpu),c_loc(coeff_deriv1_gpu),c_loc(dcsidx_gpu),c_loc(detady_gpu),c_loc(dzitdz_gpu),&
    &c_loc(psi_gpu))
    call hipCheck(hipDeviceSynchronize())


  endsubroutine insitu_div_kernel


  subroutine insitu_omega_kernel(nx,ny,nz,ng,visc_order,npsi,mpsi,w_aux_gpu,coeff_deriv1_gpu,&
  &dcsidx_gpu,detady_gpu,dzitdz_gpu,psi_gpu)
    integer :: nx,ny,nz
    integer :: ng,visc_order,npsi
    integer :: mpsi,lmax
    real(rkind), dimension(:,:,:,:), target :: w_aux_gpu
    real(rkind), dimension(:,:), target :: coeff_deriv1_gpu
    real(rkind), dimension(:), target :: dcsidx_gpu
    real(rkind), dimension(:), target :: detady_gpu
    real(rkind), dimension(:), target :: dzitdz_gpu
    real(rkind), dimension(:,:,:,:), target :: psi_gpu

    lmax = visc_order/2
    call insitu_omega_kernel_wrapper(c_null_ptr,nx,ny,nz,ng,visc_order,npsi,mpsi,lmax,&
    &c_loc(w_aux_gpu),c_loc(coeff_deriv1_gpu),c_loc(dcsidx_gpu),c_loc(detady_gpu),c_loc(dzitdz_gpu),&
    &c_loc(psi_gpu))
    call hipCheck(hipDeviceSynchronize())


  endsubroutine insitu_omega_kernel


  subroutine insitu_ducros_kernel(nx,ny,nz,ng,visc_order,npsi,mpsi,u0,l0,w_aux_gpu,coeff_deriv1_gpu,&
  &dcsidx_gpu,detady_gpu,dzitdz_gpu,psi_gpu)
    integer :: nx,ny,nz
    integer :: ng,visc_order,npsi
    integer :: mpsi,lmax
    real(rkind) :: u0,l0
    real(rkind), dimension(:,:,:,:), target :: w_aux_gpu
    real(rkind), dimension(:,:), target :: coeff_deriv1_gpu
    real(rkind), dimension(:), target :: dcsidx_gpu
    real(rkind), dimension(:), target :: detady_gpu
    real(rkind), dimension(:), target :: dzitdz_gpu
    real(rkind), dimension(:,:,:,:), target :: psi_gpu

    lmax = visc_order/2
    call insitu_ducros_kernel_wrapper(c_null_ptr,nx,ny,nz,ng,visc_order,npsi,mpsi,lmax,u0,l0,&
    &c_loc(w_aux_gpu),c_loc(coeff_deriv1_gpu),c_loc(dcsidx_gpu),c_loc(detady_gpu),c_loc(dzitdz_gpu),&
    &c_loc(psi_gpu))
    call hipCheck(hipDeviceSynchronize())


  endsubroutine insitu_ducros_kernel




  subroutine probe_interpolation_kernel(num_probe,nx,ny,nz,ng,nv_aux,ijk_probe_gpu,w_aux_probe_gpu,w_aux_gpu,probe_coeff_gpu)
    integer :: num_probe,nx,ny
    integer :: nz,ng,nv_aux
    real(rkind), dimension(:,:,:,:), target :: w_aux_gpu
    real(rkind), dimension(:,:), target :: w_aux_probe_gpu
    real(rkind), dimension(:,:,:,:), target :: probe_coeff_gpu
    integer, dimension(:,:), target :: ijk_probe_gpu


    call probe_interpolation_kernel_wrapper(c_null_ptr,num_probe,nx,ny,nz,ng,nv_aux,&
    &c_loc(ijk_probe_gpu),c_loc(w_aux_gpu),c_loc(w_aux_probe_gpu),c_loc(probe_coeff_gpu))
    call hipCheck(hipDeviceSynchronize())


  endsubroutine probe_interpolation_kernel


  subroutine compute_tspec_kernel(nx,ny,nz,ng,ndft,j_slice,i_win,it_win,w_aux_gpu,w_tspec_gpu)
    integer :: nx,ny,nz
    integer :: ng,ndft,j_slice
    integer :: i_win,it_win
    real(rkind) :: pi
    real(rkind), dimension(:,:,:,:), target :: w_aux_gpu
    real(rkind), dimension(:,:,:,:,:), target :: w_tspec_gpu

    pi = acos(-1._rkind)
    call compute_tspec_kernel_wrapper(c_null_ptr,nx,ny,nz,ng,ndft,j_slice,i_win,it_win,pi,c_loc(w_aux_gpu),c_loc(w_tspec_gpu))
    call hipCheck(hipDeviceSynchronize())


  endsubroutine compute_tspec_kernel


  subroutine compute_psd_tspec_kernel(nx,ny,nz,ndft,i_win,dt_tspec,w_tspec_gpu,w_psd_tspec_gpu)
    integer :: nx,ny,nz
    integer :: ndft,i_win
    real(rkind) :: dt_tspec,pi
    real(rkind), dimension(:,:,:,:,:), target :: w_tspec_gpu
    real(rkind), dimension(:,:), target :: w_psd_tspec_gpu

    pi = acos(-1._rkind)
    call compute_psd_tspec_kernel_wrapper(c_null_ptr,nx,ny,nz,ndft,i_win,dt_tspec,pi,c_loc(w_tspec_gpu),c_loc(w_psd_tspec_gpu))
    call hipCheck(hipDeviceSynchronize())


  endsubroutine compute_psd_tspec_kernel


  subroutine compute_wallprop_c2_kernel(nx,ny,nz,ng,w_aux_gpu,wallprop_gpu,dxdcsic2_gpu,&
  &dydcsic2_gpu,csimod_gpu,meta_gpu,wall_tag_gpu)
    integer :: nx,ny,nz
    integer :: ng
    real(rkind), dimension(:,:,:,:), target :: w_aux_gpu
    real(rkind), dimension(:,:,:), target :: wallprop_gpu
    real(rkind), dimension(:,:), target :: meta_gpu
    real(rkind), dimension(:,:), target :: csimod_gpu
    real(rkind), dimension(:,:), target :: dxdcsic2_gpu
    real(rkind), dimension(:,:), target :: dydcsic2_gpu
    integer, dimension(:), target :: wall_tag_gpu


    call compute_wallprop_c2_kernel_wrapper(c_null_ptr,nx,ny,nz,ng,c_loc(wall_tag_gpu),&
    &c_loc(w_aux_gpu),c_loc(wallprop_gpu),c_loc(meta_gpu),c_loc(csimod_gpu),c_loc(dxdcsic2_gpu),&
    &c_loc(dydcsic2_gpu))
    call hipCheck(hipDeviceSynchronize())


  endsubroutine compute_wallprop_c2_kernel


endmodule streams_kernels_amd


