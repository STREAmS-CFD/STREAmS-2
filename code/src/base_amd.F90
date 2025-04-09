module streams_base_amd_object

  use streams_field_object, only : field_object
  use streams_parameters
  use mpi
  use hipfort
  use hipfort_check

  implicit none
  interface
    subroutine bcte_step_1_kernel1_wrapper(stream,nx,ny,nz,nv,ng,nrank_x,ite_rank_x,itu_rank_x,&
    &ite_l,itu_l,w_gpu,wbuftus_gpu,wbuftes_gpu)bind(c,name="bcte_step_1_kernel1_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,nv,ng,nrank_x,ite_rank_x,itu_rank_x,ite_l,itu_l
      type(c_ptr), value :: w_gpu,wbuftus_gpu,wbuftes_gpu

    endsubroutine bcte_step_1_kernel1_wrapper
  endinterface

  interface
    subroutine bcte_step_1_kernel2_wrapper(stream,nx,ny,nz,nv,ng,nrank_x,ite_rank_x,itu_rank_x,&
    &ite_l,itu_l,w_gpu,wbuftus_gpu,wbuftes_gpu)bind(c,name="bcte_step_1_kernel2_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,nv,ng,nrank_x,ite_rank_x,itu_rank_x,ite_l,itu_l
      type(c_ptr), value :: w_gpu,wbuftus_gpu,wbuftes_gpu

    endsubroutine bcte_step_1_kernel2_wrapper
  endinterface

  interface
    subroutine bcte_step_3_kernel1_wrapper(stream,nx,ny,nz,nv,ng,nrank_x,ite_rank_x,itu_rank_x,&
    &ite_l,itu_l,w_gpu,wbuftur_gpu,wbufter_gpu)bind(c,name="bcte_step_3_kernel1_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,nv,ng,nrank_x,ite_rank_x,itu_rank_x,ite_l,itu_l
      type(c_ptr), value :: w_gpu,wbuftur_gpu,wbufter_gpu

    endsubroutine bcte_step_3_kernel1_wrapper
  endinterface

  interface
    subroutine bcte_step_3_kernel2_wrapper(stream,nx,ny,nz,nv,ng,nrank_x,ite_rank_x,itu_rank_x,&
    &ite_l,itu_l,w_gpu,wbuftur_gpu,wbufter_gpu)bind(c,name="bcte_step_3_kernel2_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,nv,ng,nrank_x,ite_rank_x,itu_rank_x,ite_l,itu_l
      type(c_ptr), value :: w_gpu,wbuftur_gpu,wbufter_gpu

    endsubroutine bcte_step_3_kernel2_wrapper
  endinterface

  interface
    subroutine bcswap_wake_step_1_kernel_wrapper(stream,nx,ny,nz,ng,nv,w_gpu,wbuf4s_gpu)bind(c,&
    &name="bcswap_wake_step_1_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv
      type(c_ptr), value :: w_gpu,wbuf4s_gpu

    endsubroutine bcswap_wake_step_1_kernel_wrapper
  endinterface

  interface
    subroutine bcswap_wake_step_3_kernel_wrapper(stream,nx,ny,nz,ng,nv,wall_tag_gpu,w_gpu,&
    &wbuf3r_gpu)bind(c,name="bcswap_wake_step_3_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv
      type(c_ptr), value :: wall_tag_gpu
      type(c_ptr), value :: w_gpu,wbuf3r_gpu

    endsubroutine bcswap_wake_step_3_kernel_wrapper
  endinterface

  interface
    subroutine bcswap_step_1_kernel1_wrapper(stream,nx,ny,nz,ng,nv,ndim,w_gpu,wbuf1s_gpu,wbuf2s_gpu,&
    &wbuf3s_gpu,wbuf4s_gpu,wbuf5s_gpu,wbuf6s_gpu)bind(c,name="bcswap_step_1_kernel1_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,ndim
      type(c_ptr), value :: w_gpu,wbuf1s_gpu,wbuf2s_gpu,wbuf3s_gpu,wbuf4s_gpu,wbuf5s_gpu,wbuf6s_gpu

    endsubroutine bcswap_step_1_kernel1_wrapper
  endinterface

  interface
    subroutine bcswap_step_1_kernel2_wrapper(stream,nx,ny,nz,ng,nv,ndim,w_gpu,wbuf1s_gpu,wbuf2s_gpu,&
    &wbuf3s_gpu,wbuf4s_gpu,wbuf5s_gpu,wbuf6s_gpu)bind(c,name="bcswap_step_1_kernel2_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,ndim
      type(c_ptr), value :: w_gpu,wbuf1s_gpu,wbuf2s_gpu,wbuf3s_gpu,wbuf4s_gpu,wbuf5s_gpu,wbuf6s_gpu

    endsubroutine bcswap_step_1_kernel2_wrapper
  endinterface

  interface
    subroutine bcswap_step_1_kernel3_wrapper(stream,nx,ny,nz,ng,nv,ndim,w_gpu,wbuf1s_gpu,wbuf2s_gpu,&
    &wbuf3s_gpu,wbuf4s_gpu,wbuf5s_gpu,wbuf6s_gpu)bind(c,name="bcswap_step_1_kernel3_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,ndim
      type(c_ptr), value :: w_gpu,wbuf1s_gpu,wbuf2s_gpu,wbuf3s_gpu,wbuf4s_gpu,wbuf5s_gpu,wbuf6s_gpu

    endsubroutine bcswap_step_1_kernel3_wrapper
  endinterface

  interface
    subroutine bcswap_c2_step_1_kernel1_wrapper(stream,nx,ny,nz,ng,nv,ndim,w_gpu,wbuf1s_gpu,&
    &wbuf2s_gpu,wbuf3s_gpu,wbuf4s_gpu,wbuf5s_gpu,wbuf6s_gpu,dcsidxnc2_gpu,dcsidync2_gpu,detadxnc2_gpu,&
    &detadync2_gpu)bind(c,name="bcswap_c2_step_1_kernel1_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,ndim
      type(c_ptr), value :: w_gpu,wbuf1s_gpu,wbuf2s_gpu,wbuf3s_gpu,wbuf4s_gpu,wbuf5s_gpu,wbuf6s_gpu,&
      &dcsidxnc2_gpu,dcsidync2_gpu,detadxnc2_gpu,detadync2_gpu

    endsubroutine bcswap_c2_step_1_kernel1_wrapper
  endinterface

  interface
    subroutine bcswap_c2_step_1_kernel2_wrapper(stream,nx,ny,nz,ng,nv,ndim,w_gpu,wbuf1s_gpu,&
    &wbuf2s_gpu,wbuf3s_gpu,wbuf4s_gpu,wbuf5s_gpu,wbuf6s_gpu,dcsidxnc2_gpu,dcsidync2_gpu,detadxnc2_gpu,&
    &detadync2_gpu)bind(c,name="bcswap_c2_step_1_kernel2_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,ndim
      type(c_ptr), value :: w_gpu,wbuf1s_gpu,wbuf2s_gpu,wbuf3s_gpu,wbuf4s_gpu,wbuf5s_gpu,wbuf6s_gpu,&
      &dcsidxnc2_gpu,dcsidync2_gpu,detadxnc2_gpu,detadync2_gpu

    endsubroutine bcswap_c2_step_1_kernel2_wrapper
  endinterface

  interface
    subroutine bcswap_c2_step_1_kernel3_wrapper(stream,nx,ny,nz,ng,nv,ndim,w_gpu,wbuf1s_gpu,&
    &wbuf2s_gpu,wbuf3s_gpu,wbuf4s_gpu,wbuf5s_gpu,wbuf6s_gpu,dcsidxnc2_gpu,dcsidync2_gpu,detadxnc2_gpu,&
    &detadync2_gpu)bind(c,name="bcswap_c2_step_1_kernel3_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,ndim
      type(c_ptr), value :: w_gpu,wbuf1s_gpu,wbuf2s_gpu,wbuf3s_gpu,wbuf4s_gpu,wbuf5s_gpu,wbuf6s_gpu,&
      &dcsidxnc2_gpu,dcsidync2_gpu,detadxnc2_gpu,detadync2_gpu

    endsubroutine bcswap_c2_step_1_kernel3_wrapper
  endinterface

  interface
    subroutine bcswap_c2_step_1_kernel4_wrapper(stream,nx,ny,nz,ng,nv,ndim,w_gpu,wbuf1s_gpu,&
    &wbuf2s_gpu,wbuf3s_gpu,wbuf4s_gpu,wbuf5s_gpu,wbuf6s_gpu,dcsidxnc2_gpu,dcsidync2_gpu,detadxnc2_gpu,&
    &detadync2_gpu)bind(c,name="bcswap_c2_step_1_kernel4_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,ndim
      type(c_ptr), value :: w_gpu,wbuf1s_gpu,wbuf2s_gpu,wbuf3s_gpu,wbuf4s_gpu,wbuf5s_gpu,wbuf6s_gpu,&
      &dcsidxnc2_gpu,dcsidync2_gpu,detadxnc2_gpu,detadync2_gpu

    endsubroutine bcswap_c2_step_1_kernel4_wrapper
  endinterface

  interface
    subroutine bcswap_c2xybc_step_1_kernel1_wrapper(stream,nx,ny,nz,ng,nv,ndim,w_gpu,wbuf1xybcs_gpu,&
    &wbuf2xybcs_gpu,wbuf3s_gpu,wbuf4s_gpu,wbuf5s_gpu,wbuf6s_gpu,dcsidxnc2_gpu,dcsidync2_gpu,&
    &detadxnc2_gpu,detadync2_gpu)bind(c,name="bcswap_c2xybc_step_1_kernel1_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,ndim
      type(c_ptr), value :: w_gpu,wbuf1xybcs_gpu,wbuf2xybcs_gpu,wbuf3s_gpu,wbuf4s_gpu,wbuf5s_gpu,&
      &wbuf6s_gpu,dcsidxnc2_gpu,dcsidync2_gpu,detadxnc2_gpu,detadync2_gpu

    endsubroutine bcswap_c2xybc_step_1_kernel1_wrapper
  endinterface

  interface
    subroutine bcswap_c2xybc_step_1_kernel2_wrapper(stream,nx,ny,nz,ng,nv,ndim,w_gpu,wbuf1xybcs_gpu,&
    &wbuf2xybcs_gpu,wbuf3s_gpu,wbuf4s_gpu,wbuf5s_gpu,wbuf6s_gpu,dcsidxnc2_gpu,dcsidync2_gpu,&
    &detadxnc2_gpu,detadync2_gpu)bind(c,name="bcswap_c2xybc_step_1_kernel2_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,ndim
      type(c_ptr), value :: w_gpu,wbuf1xybcs_gpu,wbuf2xybcs_gpu,wbuf3s_gpu,wbuf4s_gpu,wbuf5s_gpu,&
      &wbuf6s_gpu,dcsidxnc2_gpu,dcsidync2_gpu,detadxnc2_gpu,detadync2_gpu

    endsubroutine bcswap_c2xybc_step_1_kernel2_wrapper
  endinterface

  interface
    subroutine bcswap_c2xybc_step_1_kernel3_wrapper(stream,nx,ny,nz,ng,nv,ndim,w_gpu,wbuf1xybcs_gpu,&
    &wbuf2xybcs_gpu,wbuf3s_gpu,wbuf4s_gpu,wbuf5s_gpu,wbuf6s_gpu,dcsidxnc2_gpu,dcsidync2_gpu,&
    &detadxnc2_gpu,detadync2_gpu)bind(c,name="bcswap_c2xybc_step_1_kernel3_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,ndim
      type(c_ptr), value :: w_gpu,wbuf1xybcs_gpu,wbuf2xybcs_gpu,wbuf3s_gpu,wbuf4s_gpu,wbuf5s_gpu,&
      &wbuf6s_gpu,dcsidxnc2_gpu,dcsidync2_gpu,detadxnc2_gpu,detadync2_gpu

    endsubroutine bcswap_c2xybc_step_1_kernel3_wrapper
  endinterface

  interface
    subroutine bcswap_c2xybc_step_1_kernel4_wrapper(stream,nx,ny,nz,ng,nv,ndim,w_gpu,wbuf1xybcs_gpu,&
    &wbuf2xybcs_gpu,wbuf3s_gpu,wbuf4s_gpu,wbuf5s_gpu,wbuf6s_gpu,dcsidxnc2_gpu,dcsidync2_gpu,&
    &detadxnc2_gpu,detadync2_gpu)bind(c,name="bcswap_c2xybc_step_1_kernel4_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,ndim
      type(c_ptr), value :: w_gpu,wbuf1xybcs_gpu,wbuf2xybcs_gpu,wbuf3s_gpu,wbuf4s_gpu,wbuf5s_gpu,&
      &wbuf6s_gpu,dcsidxnc2_gpu,dcsidync2_gpu,detadxnc2_gpu,detadync2_gpu

    endsubroutine bcswap_c2xybc_step_1_kernel4_wrapper
  endinterface

  interface
    subroutine bcswap_corner_step_1_kernel_wrapper(stream,nx,ny,nz,ng,nv,w_gpu,wbuf1s_c_gpu,&
    &wbuf2s_c_gpu,wbuf3s_c_gpu,wbuf4s_c_gpu)bind(c,name="bcswap_corner_step_1_kernel_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv
      type(c_ptr), value :: w_gpu,wbuf1s_c_gpu,wbuf2s_c_gpu,wbuf3s_c_gpu,wbuf4s_c_gpu

    endsubroutine bcswap_corner_step_1_kernel_wrapper
  endinterface

  interface
    subroutine bcswap_corner_step_1b_kernel1_wrapper(stream,nx,ny,nz,ng,nv,w_gpu,wbuf_lxly_s_gpu,&
    &wbuf_lxry_s_gpu,wbuf_rxly_s_gpu,wbuf_rxry_s_gpu,wbuf_lylz_s_gpu,wbuf_lyrz_s_gpu,wbuf_rylz_s_gpu,&
    &wbuf_ryrz_s_gpu,wbuf_lxlylz_s_gpu,wbuf_lxlyrz_s_gpu,wbuf_lxrylz_s_gpu,wbuf_lxryrz_s_gpu,&
    &wbuf_rxlylz_s_gpu,wbuf_rxlyrz_s_gpu,wbuf_rxrylz_s_gpu,wbuf_rxryrz_s_gpu)bind(c,name="bcswap_corner_s&
    &tep_1b_kernel1_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv
      type(c_ptr), value :: w_gpu,wbuf_lxly_s_gpu,wbuf_lxry_s_gpu,wbuf_rxly_s_gpu,wbuf_rxry_s_gpu,&
      &wbuf_lylz_s_gpu,wbuf_lyrz_s_gpu,wbuf_rylz_s_gpu,wbuf_ryrz_s_gpu,wbuf_lxlylz_s_gpu,wbuf_lxlyrz_s_gpu,&
      &wbuf_lxrylz_s_gpu,wbuf_lxryrz_s_gpu,wbuf_rxlylz_s_gpu,wbuf_rxlyrz_s_gpu,wbuf_rxrylz_s_gpu,&
      &wbuf_rxryrz_s_gpu

    endsubroutine bcswap_corner_step_1b_kernel1_wrapper
  endinterface

  interface
    subroutine bcswap_corner_step_1b_kernel2_wrapper(stream,nx,ny,nz,ng,nv,w_gpu,wbuf_lxly_s_gpu,&
    &wbuf_lxry_s_gpu,wbuf_rxly_s_gpu,wbuf_rxry_s_gpu,wbuf_lylz_s_gpu,wbuf_lyrz_s_gpu,wbuf_rylz_s_gpu,&
    &wbuf_ryrz_s_gpu,wbuf_lxlylz_s_gpu,wbuf_lxlyrz_s_gpu,wbuf_lxrylz_s_gpu,wbuf_lxryrz_s_gpu,&
    &wbuf_rxlylz_s_gpu,wbuf_rxlyrz_s_gpu,wbuf_rxrylz_s_gpu,wbuf_rxryrz_s_gpu)bind(c,name="bcswap_corner_s&
    &tep_1b_kernel2_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv
      type(c_ptr), value :: w_gpu,wbuf_lxly_s_gpu,wbuf_lxry_s_gpu,wbuf_rxly_s_gpu,wbuf_rxry_s_gpu,&
      &wbuf_lylz_s_gpu,wbuf_lyrz_s_gpu,wbuf_rylz_s_gpu,wbuf_ryrz_s_gpu,wbuf_lxlylz_s_gpu,wbuf_lxlyrz_s_gpu,&
      &wbuf_lxrylz_s_gpu,wbuf_lxryrz_s_gpu,wbuf_rxlylz_s_gpu,wbuf_rxlyrz_s_gpu,wbuf_rxrylz_s_gpu,&
      &wbuf_rxryrz_s_gpu

    endsubroutine bcswap_corner_step_1b_kernel2_wrapper
  endinterface

  interface
    subroutine bcswap_corner_step_1b_kernel3_wrapper(stream,nx,ny,nz,ng,nv,w_gpu,wbuf_lxly_s_gpu,&
    &wbuf_lxry_s_gpu,wbuf_rxly_s_gpu,wbuf_rxry_s_gpu,wbuf_lylz_s_gpu,wbuf_lyrz_s_gpu,wbuf_rylz_s_gpu,&
    &wbuf_ryrz_s_gpu,wbuf_lxlylz_s_gpu,wbuf_lxlyrz_s_gpu,wbuf_lxrylz_s_gpu,wbuf_lxryrz_s_gpu,&
    &wbuf_rxlylz_s_gpu,wbuf_rxlyrz_s_gpu,wbuf_rxrylz_s_gpu,wbuf_rxryrz_s_gpu)bind(c,name="bcswap_corner_s&
    &tep_1b_kernel3_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv
      type(c_ptr), value :: w_gpu,wbuf_lxly_s_gpu,wbuf_lxry_s_gpu,wbuf_rxly_s_gpu,wbuf_rxry_s_gpu,&
      &wbuf_lylz_s_gpu,wbuf_lyrz_s_gpu,wbuf_rylz_s_gpu,wbuf_ryrz_s_gpu,wbuf_lxlylz_s_gpu,wbuf_lxlyrz_s_gpu,&
      &wbuf_lxrylz_s_gpu,wbuf_lxryrz_s_gpu,wbuf_rxlylz_s_gpu,wbuf_rxlyrz_s_gpu,wbuf_rxrylz_s_gpu,&
      &wbuf_rxryrz_s_gpu

    endsubroutine bcswap_corner_step_1b_kernel3_wrapper
  endinterface

  interface
    subroutine bcswap_step_3_kernel1_wrapper(stream,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,&
    &irightx,irighty,irightz,w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,&
    &wbuf6r_gpu)bind(c,name="bcswap_step_3_kernel1_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,irightx,irighty,irightz
      type(c_ptr), value :: w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu

    endsubroutine bcswap_step_3_kernel1_wrapper
  endinterface

  interface
    subroutine bcswap_step_3_kernel2_wrapper(stream,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,&
    &irightx,irighty,irightz,w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,&
    &wbuf6r_gpu)bind(c,name="bcswap_step_3_kernel2_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,irightx,irighty,irightz
      type(c_ptr), value :: w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu

    endsubroutine bcswap_step_3_kernel2_wrapper
  endinterface

  interface
    subroutine bcswap_step_3_kernel3_wrapper(stream,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,&
    &irightx,irighty,irightz,w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,&
    &wbuf6r_gpu)bind(c,name="bcswap_step_3_kernel3_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,irightx,irighty,irightz
      type(c_ptr), value :: w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu

    endsubroutine bcswap_step_3_kernel3_wrapper
  endinterface

  interface
    subroutine bcswap_step_3_kernel4_wrapper(stream,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,&
    &irightx,irighty,irightz,w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,&
    &wbuf6r_gpu)bind(c,name="bcswap_step_3_kernel4_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,irightx,irighty,irightz
      type(c_ptr), value :: w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu

    endsubroutine bcswap_step_3_kernel4_wrapper
  endinterface

  interface
    subroutine bcswap_step_3_kernel5_wrapper(stream,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,&
    &irightx,irighty,irightz,w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,&
    &wbuf6r_gpu)bind(c,name="bcswap_step_3_kernel5_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,irightx,irighty,irightz
      type(c_ptr), value :: w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu

    endsubroutine bcswap_step_3_kernel5_wrapper
  endinterface

  interface
    subroutine bcswap_step_3_kernel6_wrapper(stream,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,&
    &irightx,irighty,irightz,w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,&
    &wbuf6r_gpu)bind(c,name="bcswap_step_3_kernel6_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,irightx,irighty,irightz
      type(c_ptr), value :: w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu

    endsubroutine bcswap_step_3_kernel6_wrapper
  endinterface

  interface
    subroutine bcswap_c2_step_3_kernel1_wrapper(stream,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,&
    &irightx,irighty,irightz,w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,&
    &dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu)bind(c,name="bcswap_c2_step_3_kernel1_wrappe&
    &r")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,irightx,irighty,irightz
      type(c_ptr), value :: w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,&
      &dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu

    endsubroutine bcswap_c2_step_3_kernel1_wrapper
  endinterface

  interface
    subroutine bcswap_c2_step_3_kernel2_wrapper(stream,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,&
    &irightx,irighty,irightz,w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,&
    &dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu)bind(c,name="bcswap_c2_step_3_kernel2_wrappe&
    &r")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,irightx,irighty,irightz
      type(c_ptr), value :: w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,&
      &dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu

    endsubroutine bcswap_c2_step_3_kernel2_wrapper
  endinterface

  interface
    subroutine bcswap_c2_step_3_kernel3_wrapper(stream,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,&
    &irightx,irighty,irightz,w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,&
    &dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu)bind(c,name="bcswap_c2_step_3_kernel3_wrappe&
    &r")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,irightx,irighty,irightz
      type(c_ptr), value :: w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,&
      &dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu

    endsubroutine bcswap_c2_step_3_kernel3_wrapper
  endinterface

  interface
    subroutine bcswap_c2_step_3_kernel4_wrapper(stream,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,&
    &irightx,irighty,irightz,w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,&
    &dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu)bind(c,name="bcswap_c2_step_3_kernel4_wrappe&
    &r")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,irightx,irighty,irightz
      type(c_ptr), value :: w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,&
      &dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu

    endsubroutine bcswap_c2_step_3_kernel4_wrapper
  endinterface

  interface
    subroutine bcswap_c2_step_3_kernel5_wrapper(stream,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,&
    &irightx,irighty,irightz,w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,&
    &dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu)bind(c,name="bcswap_c2_step_3_kernel5_wrappe&
    &r")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,irightx,irighty,irightz
      type(c_ptr), value :: w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,&
      &dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu

    endsubroutine bcswap_c2_step_3_kernel5_wrapper
  endinterface

  interface
    subroutine bcswap_c2_step_3_kernel6_wrapper(stream,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,&
    &irightx,irighty,irightz,w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,&
    &dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu)bind(c,name="bcswap_c2_step_3_kernel6_wrappe&
    &r")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,irightx,irighty,irightz
      type(c_ptr), value :: w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,&
      &dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu

    endsubroutine bcswap_c2_step_3_kernel6_wrapper
  endinterface

  interface
    subroutine bcswap_c2_step_3_kernel7_wrapper(stream,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,&
    &irightx,irighty,irightz,w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,&
    &dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu)bind(c,name="bcswap_c2_step_3_kernel7_wrappe&
    &r")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,irightx,irighty,irightz
      type(c_ptr), value :: w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,&
      &dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu

    endsubroutine bcswap_c2_step_3_kernel7_wrapper
  endinterface

  interface
    subroutine bcswap_c2_step_3_kernel8_wrapper(stream,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,&
    &irightx,irighty,irightz,w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,&
    &dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu)bind(c,name="bcswap_c2_step_3_kernel8_wrappe&
    &r")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,irightx,irighty,irightz
      type(c_ptr), value :: w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,&
      &dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu

    endsubroutine bcswap_c2_step_3_kernel8_wrapper
  endinterface

  interface
    subroutine bcswap_c2xybc_step_3_kernel1_wrapper(stream,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,&
    &irightx,irighty,irightz,w_gpu,wbuf1xybcr_gpu,wbuf2xybcr_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,&
    &wbuf6r_gpu,dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu)bind(c,name="bcswap_c2xybc_step_3&
    &_kernel1_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,irightx,irighty,irightz
      type(c_ptr), value :: w_gpu,wbuf1xybcr_gpu,wbuf2xybcr_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,&
      &wbuf6r_gpu,dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu

    endsubroutine bcswap_c2xybc_step_3_kernel1_wrapper
  endinterface

  interface
    subroutine bcswap_c2xybc_step_3_kernel2_wrapper(stream,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,&
    &irightx,irighty,irightz,w_gpu,wbuf1xybcr_gpu,wbuf2xybcr_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,&
    &wbuf6r_gpu,dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu)bind(c,name="bcswap_c2xybc_step_3&
    &_kernel2_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,irightx,irighty,irightz
      type(c_ptr), value :: w_gpu,wbuf1xybcr_gpu,wbuf2xybcr_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,&
      &wbuf6r_gpu,dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu

    endsubroutine bcswap_c2xybc_step_3_kernel2_wrapper
  endinterface

  interface
    subroutine bcswap_c2xybc_step_3_kernel3_wrapper(stream,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,&
    &irightx,irighty,irightz,w_gpu,wbuf1xybcr_gpu,wbuf2xybcr_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,&
    &wbuf6r_gpu,dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu)bind(c,name="bcswap_c2xybc_step_3&
    &_kernel3_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,irightx,irighty,irightz
      type(c_ptr), value :: w_gpu,wbuf1xybcr_gpu,wbuf2xybcr_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,&
      &wbuf6r_gpu,dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu

    endsubroutine bcswap_c2xybc_step_3_kernel3_wrapper
  endinterface

  interface
    subroutine bcswap_c2xybc_step_3_kernel4_wrapper(stream,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,&
    &irightx,irighty,irightz,w_gpu,wbuf1xybcr_gpu,wbuf2xybcr_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,&
    &wbuf6r_gpu,dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu)bind(c,name="bcswap_c2xybc_step_3&
    &_kernel4_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,irightx,irighty,irightz
      type(c_ptr), value :: w_gpu,wbuf1xybcr_gpu,wbuf2xybcr_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,&
      &wbuf6r_gpu,dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu

    endsubroutine bcswap_c2xybc_step_3_kernel4_wrapper
  endinterface

  interface
    subroutine bcswap_c2xybc_step_3_kernel5_wrapper(stream,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,&
    &irightx,irighty,irightz,w_gpu,wbuf1xybcr_gpu,wbuf2xybcr_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,&
    &wbuf6r_gpu,dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu)bind(c,name="bcswap_c2xybc_step_3&
    &_kernel5_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,irightx,irighty,irightz
      type(c_ptr), value :: w_gpu,wbuf1xybcr_gpu,wbuf2xybcr_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,&
      &wbuf6r_gpu,dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu

    endsubroutine bcswap_c2xybc_step_3_kernel5_wrapper
  endinterface

  interface
    subroutine bcswap_c2xybc_step_3_kernel6_wrapper(stream,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,&
    &irightx,irighty,irightz,w_gpu,wbuf1xybcr_gpu,wbuf2xybcr_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,&
    &wbuf6r_gpu,dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu)bind(c,name="bcswap_c2xybc_step_3&
    &_kernel6_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,irightx,irighty,irightz
      type(c_ptr), value :: w_gpu,wbuf1xybcr_gpu,wbuf2xybcr_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,&
      &wbuf6r_gpu,dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu

    endsubroutine bcswap_c2xybc_step_3_kernel6_wrapper
  endinterface

  interface
    subroutine bcswap_c2xybc_step_3_kernel7_wrapper(stream,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,&
    &irightx,irighty,irightz,w_gpu,wbuf1xybcr_gpu,wbuf2xybcr_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,&
    &wbuf6r_gpu,dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu)bind(c,name="bcswap_c2xybc_step_3&
    &_kernel7_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,irightx,irighty,irightz
      type(c_ptr), value :: w_gpu,wbuf1xybcr_gpu,wbuf2xybcr_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,&
      &wbuf6r_gpu,dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu

    endsubroutine bcswap_c2xybc_step_3_kernel7_wrapper
  endinterface

  interface
    subroutine bcswap_c2xybc_step_3_kernel8_wrapper(stream,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,&
    &irightx,irighty,irightz,w_gpu,wbuf1xybcr_gpu,wbuf2xybcr_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,&
    &wbuf6r_gpu,dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu)bind(c,name="bcswap_c2xybc_step_3&
    &_kernel8_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,irightx,irighty,irightz
      type(c_ptr), value :: w_gpu,wbuf1xybcr_gpu,wbuf2xybcr_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,&
      &wbuf6r_gpu,dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu

    endsubroutine bcswap_c2xybc_step_3_kernel8_wrapper
  endinterface

  interface
    subroutine bcswap_corner_step_3_kernel1_wrapper(stream,nx,ny,nz,ng,nv,ileftbottom,ilefttop,&
    &irightbottom,irighttop,w_gpu,wbuf1r_c_gpu,wbuf2r_c_gpu,wbuf3r_c_gpu,wbuf4r_c_gpu)bind(c,&
    &name="bcswap_corner_step_3_kernel1_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,ileftbottom,ilefttop,irightbottom,irighttop
      type(c_ptr), value :: w_gpu,wbuf1r_c_gpu,wbuf2r_c_gpu,wbuf3r_c_gpu,wbuf4r_c_gpu

    endsubroutine bcswap_corner_step_3_kernel1_wrapper
  endinterface

  interface
    subroutine bcswap_corner_step_3_kernel2_wrapper(stream,nx,ny,nz,ng,nv,ileftbottom,ilefttop,&
    &irightbottom,irighttop,w_gpu,wbuf1r_c_gpu,wbuf2r_c_gpu,wbuf3r_c_gpu,wbuf4r_c_gpu)bind(c,&
    &name="bcswap_corner_step_3_kernel2_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,ileftbottom,ilefttop,irightbottom,irighttop
      type(c_ptr), value :: w_gpu,wbuf1r_c_gpu,wbuf2r_c_gpu,wbuf3r_c_gpu,wbuf4r_c_gpu

    endsubroutine bcswap_corner_step_3_kernel2_wrapper
  endinterface

  interface
    subroutine bcswap_corner_step_3_kernel3_wrapper(stream,nx,ny,nz,ng,nv,ileftbottom,ilefttop,&
    &irightbottom,irighttop,w_gpu,wbuf1r_c_gpu,wbuf2r_c_gpu,wbuf3r_c_gpu,wbuf4r_c_gpu)bind(c,&
    &name="bcswap_corner_step_3_kernel3_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,ileftbottom,ilefttop,irightbottom,irighttop
      type(c_ptr), value :: w_gpu,wbuf1r_c_gpu,wbuf2r_c_gpu,wbuf3r_c_gpu,wbuf4r_c_gpu

    endsubroutine bcswap_corner_step_3_kernel3_wrapper
  endinterface

  interface
    subroutine bcswap_corner_step_3_kernel4_wrapper(stream,nx,ny,nz,ng,nv,ileftbottom,ilefttop,&
    &irightbottom,irighttop,w_gpu,wbuf1r_c_gpu,wbuf2r_c_gpu,wbuf3r_c_gpu,wbuf4r_c_gpu)bind(c,&
    &name="bcswap_corner_step_3_kernel4_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,ileftbottom,ilefttop,irightbottom,irighttop
      type(c_ptr), value :: w_gpu,wbuf1r_c_gpu,wbuf2r_c_gpu,wbuf3r_c_gpu,wbuf4r_c_gpu

    endsubroutine bcswap_corner_step_3_kernel4_wrapper
  endinterface

  interface
    subroutine bcswap_corner_step_3b_kernel1_wrapper(stream,nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,lylz,&
    &lyrz,rylz,ryrz,lxlylz,lxlyrz,lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz,w_gpu,wbuf_lxly_r_gpu,&
    &wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,&
    &wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,&
    &wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu)bind(c,name="bcswap_corner_s&
    &tep_3b_kernel1_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,lylz,lyrz,rylz,ryrz,lxlylz,lxlyrz,&
      &lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz
      type(c_ptr), value :: w_gpu,wbuf_lxly_r_gpu,wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,&
      &wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,&
      &wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,&
      &wbuf_rxryrz_r_gpu

    endsubroutine bcswap_corner_step_3b_kernel1_wrapper
  endinterface

  interface
    subroutine bcswap_corner_step_3b_kernel2_wrapper(stream,nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,lylz,&
    &lyrz,rylz,ryrz,lxlylz,lxlyrz,lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz,w_gpu,wbuf_lxly_r_gpu,&
    &wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,&
    &wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,&
    &wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu)bind(c,name="bcswap_corner_s&
    &tep_3b_kernel2_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,lylz,lyrz,rylz,ryrz,lxlylz,lxlyrz,&
      &lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz
      type(c_ptr), value :: w_gpu,wbuf_lxly_r_gpu,wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,&
      &wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,&
      &wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,&
      &wbuf_rxryrz_r_gpu

    endsubroutine bcswap_corner_step_3b_kernel2_wrapper
  endinterface

  interface
    subroutine bcswap_corner_step_3b_kernel3_wrapper(stream,nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,lylz,&
    &lyrz,rylz,ryrz,lxlylz,lxlyrz,lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz,w_gpu,wbuf_lxly_r_gpu,&
    &wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,&
    &wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,&
    &wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu)bind(c,name="bcswap_corner_s&
    &tep_3b_kernel3_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,lylz,lyrz,rylz,ryrz,lxlylz,lxlyrz,&
      &lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz
      type(c_ptr), value :: w_gpu,wbuf_lxly_r_gpu,wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,&
      &wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,&
      &wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,&
      &wbuf_rxryrz_r_gpu

    endsubroutine bcswap_corner_step_3b_kernel3_wrapper
  endinterface

  interface
    subroutine bcswap_corner_step_3b_kernel4_wrapper(stream,nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,lylz,&
    &lyrz,rylz,ryrz,lxlylz,lxlyrz,lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz,w_gpu,wbuf_lxly_r_gpu,&
    &wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,&
    &wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,&
    &wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu)bind(c,name="bcswap_corner_s&
    &tep_3b_kernel4_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,lylz,lyrz,rylz,ryrz,lxlylz,lxlyrz,&
      &lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz
      type(c_ptr), value :: w_gpu,wbuf_lxly_r_gpu,wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,&
      &wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,&
      &wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,&
      &wbuf_rxryrz_r_gpu

    endsubroutine bcswap_corner_step_3b_kernel4_wrapper
  endinterface

  interface
    subroutine bcswap_corner_step_3b_kernel5_wrapper(stream,nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,lylz,&
    &lyrz,rylz,ryrz,lxlylz,lxlyrz,lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz,w_gpu,wbuf_lxly_r_gpu,&
    &wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,&
    &wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,&
    &wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu)bind(c,name="bcswap_corner_s&
    &tep_3b_kernel5_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,lylz,lyrz,rylz,ryrz,lxlylz,lxlyrz,&
      &lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz
      type(c_ptr), value :: w_gpu,wbuf_lxly_r_gpu,wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,&
      &wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,&
      &wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,&
      &wbuf_rxryrz_r_gpu

    endsubroutine bcswap_corner_step_3b_kernel5_wrapper
  endinterface

  interface
    subroutine bcswap_corner_step_3b_kernel6_wrapper(stream,nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,lylz,&
    &lyrz,rylz,ryrz,lxlylz,lxlyrz,lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz,w_gpu,wbuf_lxly_r_gpu,&
    &wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,&
    &wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,&
    &wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu)bind(c,name="bcswap_corner_s&
    &tep_3b_kernel6_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,lylz,lyrz,rylz,ryrz,lxlylz,lxlyrz,&
      &lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz
      type(c_ptr), value :: w_gpu,wbuf_lxly_r_gpu,wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,&
      &wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,&
      &wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,&
      &wbuf_rxryrz_r_gpu

    endsubroutine bcswap_corner_step_3b_kernel6_wrapper
  endinterface

  interface
    subroutine bcswap_corner_step_3b_kernel7_wrapper(stream,nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,lylz,&
    &lyrz,rylz,ryrz,lxlylz,lxlyrz,lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz,w_gpu,wbuf_lxly_r_gpu,&
    &wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,&
    &wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,&
    &wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu)bind(c,name="bcswap_corner_s&
    &tep_3b_kernel7_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,lylz,lyrz,rylz,ryrz,lxlylz,lxlyrz,&
      &lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz
      type(c_ptr), value :: w_gpu,wbuf_lxly_r_gpu,wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,&
      &wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,&
      &wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,&
      &wbuf_rxryrz_r_gpu

    endsubroutine bcswap_corner_step_3b_kernel7_wrapper
  endinterface

  interface
    subroutine bcswap_corner_step_3b_kernel8_wrapper(stream,nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,lylz,&
    &lyrz,rylz,ryrz,lxlylz,lxlyrz,lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz,w_gpu,wbuf_lxly_r_gpu,&
    &wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,&
    &wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,&
    &wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu)bind(c,name="bcswap_corner_s&
    &tep_3b_kernel8_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,lylz,lyrz,rylz,ryrz,lxlylz,lxlyrz,&
      &lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz
      type(c_ptr), value :: w_gpu,wbuf_lxly_r_gpu,wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,&
      &wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,&
      &wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,&
      &wbuf_rxryrz_r_gpu

    endsubroutine bcswap_corner_step_3b_kernel8_wrapper
  endinterface

  interface
    subroutine bcswap_corner_step_3b_kernel9_wrapper(stream,nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,lylz,&
    &lyrz,rylz,ryrz,lxlylz,lxlyrz,lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz,w_gpu,wbuf_lxly_r_gpu,&
    &wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,&
    &wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,&
    &wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu)bind(c,name="bcswap_corner_s&
    &tep_3b_kernel9_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,lylz,lyrz,rylz,ryrz,lxlylz,lxlyrz,&
      &lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz
      type(c_ptr), value :: w_gpu,wbuf_lxly_r_gpu,wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,&
      &wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,&
      &wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,&
      &wbuf_rxryrz_r_gpu

    endsubroutine bcswap_corner_step_3b_kernel9_wrapper
  endinterface

  interface
    subroutine bcswap_corner_step_3b_kernel10_wrapper(stream,nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,&
    &lylz,lyrz,rylz,ryrz,lxlylz,lxlyrz,lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz,w_gpu,wbuf_lxly_r_gpu,&
    &wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,&
    &wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,&
    &wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu)bind(c,name="bcswap_corner_s&
    &tep_3b_kernel10_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,lylz,lyrz,rylz,ryrz,lxlylz,lxlyrz,&
      &lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz
      type(c_ptr), value :: w_gpu,wbuf_lxly_r_gpu,wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,&
      &wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,&
      &wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,&
      &wbuf_rxryrz_r_gpu

    endsubroutine bcswap_corner_step_3b_kernel10_wrapper
  endinterface

  interface
    subroutine bcswap_corner_step_3b_kernel11_wrapper(stream,nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,&
    &lylz,lyrz,rylz,ryrz,lxlylz,lxlyrz,lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz,w_gpu,wbuf_lxly_r_gpu,&
    &wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,&
    &wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,&
    &wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu)bind(c,name="bcswap_corner_s&
    &tep_3b_kernel11_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,lylz,lyrz,rylz,ryrz,lxlylz,lxlyrz,&
      &lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz
      type(c_ptr), value :: w_gpu,wbuf_lxly_r_gpu,wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,&
      &wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,&
      &wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,&
      &wbuf_rxryrz_r_gpu

    endsubroutine bcswap_corner_step_3b_kernel11_wrapper
  endinterface

  interface
    subroutine bcswap_corner_step_3b_kernel12_wrapper(stream,nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,&
    &lylz,lyrz,rylz,ryrz,lxlylz,lxlyrz,lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz,w_gpu,wbuf_lxly_r_gpu,&
    &wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,&
    &wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,&
    &wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu)bind(c,name="bcswap_corner_s&
    &tep_3b_kernel12_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,lylz,lyrz,rylz,ryrz,lxlylz,lxlyrz,&
      &lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz
      type(c_ptr), value :: w_gpu,wbuf_lxly_r_gpu,wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,&
      &wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,&
      &wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,&
      &wbuf_rxryrz_r_gpu

    endsubroutine bcswap_corner_step_3b_kernel12_wrapper
  endinterface

  interface
    subroutine bcswap_corner_step_3b_kernel13_wrapper(stream,nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,&
    &lylz,lyrz,rylz,ryrz,lxlylz,lxlyrz,lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz,w_gpu,wbuf_lxly_r_gpu,&
    &wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,&
    &wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,&
    &wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu)bind(c,name="bcswap_corner_s&
    &tep_3b_kernel13_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,lylz,lyrz,rylz,ryrz,lxlylz,lxlyrz,&
      &lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz
      type(c_ptr), value :: w_gpu,wbuf_lxly_r_gpu,wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,&
      &wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,&
      &wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,&
      &wbuf_rxryrz_r_gpu

    endsubroutine bcswap_corner_step_3b_kernel13_wrapper
  endinterface

  interface
    subroutine bcswap_corner_step_3b_kernel14_wrapper(stream,nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,&
    &lylz,lyrz,rylz,ryrz,lxlylz,lxlyrz,lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz,w_gpu,wbuf_lxly_r_gpu,&
    &wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,&
    &wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,&
    &wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu)bind(c,name="bcswap_corner_s&
    &tep_3b_kernel14_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,lylz,lyrz,rylz,ryrz,lxlylz,lxlyrz,&
      &lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz
      type(c_ptr), value :: w_gpu,wbuf_lxly_r_gpu,wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,&
      &wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,&
      &wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,&
      &wbuf_rxryrz_r_gpu

    endsubroutine bcswap_corner_step_3b_kernel14_wrapper
  endinterface

  interface
    subroutine bcswap_corner_step_3b_kernel15_wrapper(stream,nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,&
    &lylz,lyrz,rylz,ryrz,lxlylz,lxlyrz,lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz,w_gpu,wbuf_lxly_r_gpu,&
    &wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,&
    &wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,&
    &wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu)bind(c,name="bcswap_corner_s&
    &tep_3b_kernel15_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,lylz,lyrz,rylz,ryrz,lxlylz,lxlyrz,&
      &lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz
      type(c_ptr), value :: w_gpu,wbuf_lxly_r_gpu,wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,&
      &wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,&
      &wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,&
      &wbuf_rxryrz_r_gpu

    endsubroutine bcswap_corner_step_3b_kernel15_wrapper
  endinterface

  interface
    subroutine bcswap_corner_step_3b_kernel16_wrapper(stream,nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,&
    &lylz,lyrz,rylz,ryrz,lxlylz,lxlyrz,lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz,w_gpu,wbuf_lxly_r_gpu,&
    &wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,&
    &wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,&
    &wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu)bind(c,name="bcswap_corner_s&
    &tep_3b_kernel16_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,lylz,lyrz,rylz,ryrz,lxlylz,lxlyrz,&
      &lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz
      type(c_ptr), value :: w_gpu,wbuf_lxly_r_gpu,wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,&
      &wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,&
      &wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,&
      &wbuf_rxryrz_r_gpu

    endsubroutine bcswap_corner_step_3b_kernel16_wrapper
  endinterface

  interface
    subroutine extr_corner_ymin_kernel1_wrapper(stream,nx,ny,nz,ng,nv,w_gpu)bind(c,name="extr_corner_ymin_kernel1_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv
      type(c_ptr), value :: w_gpu

    endsubroutine extr_corner_ymin_kernel1_wrapper
  endinterface

  interface
    subroutine extr_corner_ymin_kernel2_wrapper(stream,nx,ny,nz,ng,nv,w_gpu)bind(c,name="extr_corner_ymin_kernel2_wrapper")

      import :: c_ptr, c_rkind, c_bool, c_int
      implicit none

      type(c_ptr), value :: stream
      integer(c_int), value :: nx,ny,nz,ng,nv
      type(c_ptr), value :: w_gpu

    endsubroutine extr_corner_ymin_kernel2_wrapper
  endinterface
  private
  public :: base_amd_object

  type :: base_amd_object
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
    integer(ikind) :: gpu_bind=0_ikind
    real(rkind), dimension(:,:,:,:), pointer :: w_gpu

    real(rkind), dimension(:,:,:,:), pointer :: wbuf1s_gpu, wbuf2s_gpu, wbuf3s_gpu, wbuf4s_gpu, wbuf5s_gpu, wbuf6s_gpu
    real(rkind), dimension(:,:,:,:), pointer :: wbuf1r_gpu, wbuf2r_gpu, wbuf3r_gpu, wbuf4r_gpu, wbuf5r_gpu, wbuf6r_gpu
    real(rkind), dimension(:,:,:,:), pointer :: wbuf1xybcs_gpu, wbuf2xybcs_gpu, wbuf1xybcr_gpu, wbuf2xybcr_gpu
    real(rkind), dimension(:,:,:,:), pointer :: wbuf1s_c_gpu, wbuf2s_c_gpu, wbuf3s_c_gpu, wbuf4s_c_gpu
    real(rkind), dimension(:,:,:,:), pointer :: wbuf1r_c_gpu, wbuf2r_c_gpu, wbuf3r_c_gpu, wbuf4r_c_gpu
    real(rkind), dimension(:), pointer :: dcsidx_gpu, dcsidx2_gpu, dcsidxs_gpu
    real(rkind), dimension(:), pointer :: detady_gpu, detady2_gpu, detadys_gpu
    real(rkind), dimension(:), pointer :: dzitdz_gpu, dzitdz2_gpu, dzitdzs_gpu
    real(rkind), dimension(:), pointer :: x_gpu, y_gpu, z_gpu, yn_gpu
    real(rkind), dimension(:,:), pointer :: xc2_gpu,yc2_gpu
    real(rkind), dimension(:,:), pointer :: dcsidxc2_gpu,dcsidyc2_gpu
    real(rkind), dimension(:,:), pointer :: detadxc2_gpu,detadyc2_gpu
    real(rkind), dimension(:,:), pointer :: dcsidxnc2_gpu,dcsidync2_gpu
    real(rkind), dimension(:,:), pointer :: detadxnc2_gpu,detadync2_gpu
    real(rkind), dimension(:,:), pointer :: dxdcsic2_gpu,dydcsic2_gpu
    real(rkind), dimension(:,:), pointer :: dxdetac2_gpu,dydetac2_gpu
    real(rkind), dimension(:,:), pointer :: dxdcsinc2_gpu,dydcsinc2_gpu
    real(rkind), dimension(:,:), pointer :: dxdetanc2_gpu,dydetanc2_gpu
    real(rkind), dimension(:,:), pointer :: jac_gpu,mcsijac1_gpu,metajac1_gpu
    real(rkind), dimension(:,:), pointer :: mcsi_gpu,meta_gpu,csimod_gpu,etamod_gpu
    real(rkind), dimension(:,:), pointer :: g1_gpu,g2_gpu,g12_gpu
    real(rkind), dimension(:,:), pointer :: wbuftus_gpu, wbuftes_gpu , wbuftur_gpu, wbufter_gpu
    real(rkind), dimension(:,:), pointer :: theta_ij_gpu

    integer, dimension(:), pointer :: wall_tag_gpu

    real(rkind),dimension(:,:,:,:),pointer :: wbuf_lxly_s_gpu , wbuf_lxry_s_gpu , wbuf_rxly_s_gpu , wbuf_rxry_s_gpu
    real(rkind),dimension(:,:,:,:),pointer :: wbuf_lxlz_s_gpu , wbuf_lxrz_s_gpu , wbuf_rxlz_s_gpu , wbuf_rxrz_s_gpu
    real(rkind),dimension(:,:,:,:),pointer :: wbuf_lylz_s_gpu , wbuf_lyrz_s_gpu , wbuf_rylz_s_gpu , wbuf_ryrz_s_gpu
    real(rkind),dimension(:,:,:,:),pointer :: wbuf_lxlylz_s_gpu,wbuf_lxlyrz_s_gpu,wbuf_lxrylz_s_gpu,wbuf_lxryrz_s_gpu
    real(rkind),dimension(:,:,:,:),pointer :: wbuf_rxlylz_s_gpu,wbuf_rxlyrz_s_gpu,wbuf_rxrylz_s_gpu,wbuf_rxryrz_s_gpu
    real(rkind),dimension(:,:,:,:),pointer :: wbuf_lxly_r_gpu , wbuf_lxry_r_gpu , wbuf_rxly_r_gpu , wbuf_rxry_r_gpu
    real(rkind),dimension(:,:,:,:),pointer :: wbuf_lxlz_r_gpu , wbuf_lxrz_r_gpu , wbuf_rxlz_r_gpu , wbuf_rxrz_r_gpu
    real(rkind),dimension(:,:,:,:),pointer :: wbuf_lylz_r_gpu , wbuf_lyrz_r_gpu , wbuf_rylz_r_gpu , wbuf_ryrz_r_gpu
    real(rkind),dimension(:,:,:,:),pointer :: wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu
    real(rkind),dimension(:,:,:,:),pointer :: wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu

  contains
    procedure, pass(self) :: alloc
    procedure, pass(self) :: copy_cpu_gpu
    procedure, pass(self) :: copy_gpu_cpu
    procedure, pass(self) :: initialize
    procedure, pass(self) :: bcswap
    procedure, pass(self) :: bcswap_var
    procedure, pass(self) :: bcswap_corner
    procedure, pass(self) :: bcswap_corner_var
    procedure, pass(self) :: bcswap_wake
    procedure, pass(self) :: bcswap_wake_var
    procedure, pass(self) :: bcte
    procedure, pass(self) :: bcte_var
    procedure, pass(self) :: check_gpu_mem
    procedure, pass(self) :: bcswap_edges_corners
    procedure, pass(self) :: bcswap_edges_corners_var
  endtype base_amd_object

contains
  subroutine alloc(self)
    class(base_amd_object), intent(inout) :: self

    associate(nx => self%field%nx, ny => self%field%ny, nz => self%field%nz, ng => self%field%grid%ng, nv => self%field%nv)

      call hipCheck(hipMalloc(self%w_gpu,(nx+ng)-(1-ng)+1,(ny+ng)-(1-ng)+1,(nz+ng)-(1-ng)+1,(nv)-(1)+1))

      call hipCheck(hipMalloc(self%wbuf1s_gpu,(ng)-(1)+1,(ny)-(1)+1,(nz)-(1)+1,(nv)-(1)+1))
      call hipCheck(hipMalloc(self%wbuf2s_gpu,(ng)-(1)+1,(ny)-(1)+1,(nz)-(1)+1,(nv)-(1)+1))
      call hipCheck(hipMalloc(self%wbuf3s_gpu,(nx)-(1)+1,(ng)-(1)+1,(nz)-(1)+1,(nv)-(1)+1))
      call hipCheck(hipMalloc(self%wbuf4s_gpu,(nx)-(1)+1,(ng)-(1)+1,(nz)-(1)+1,(nv)-(1)+1))
      call hipCheck(hipMalloc(self%wbuf5s_gpu,(nx)-(1)+1,(ny)-(1)+1,(ng)-(1)+1,(nv)-(1)+1))
      call hipCheck(hipMalloc(self%wbuf6s_gpu,(nx)-(1)+1,(ny)-(1)+1,(ng)-(1)+1,(nv)-(1)+1))
      call hipCheck(hipMalloc(self%wbuf1r_gpu,(ng)-(1)+1,(ny)-(1)+1,(nz)-(1)+1,(nv)-(1)+1))
      call hipCheck(hipMalloc(self%wbuf2r_gpu,(ng)-(1)+1,(ny)-(1)+1,(nz)-(1)+1,(nv)-(1)+1))
      call hipCheck(hipMalloc(self%wbuf3r_gpu,(nx)-(1)+1,(ng)-(1)+1,(nz)-(1)+1,(nv)-(1)+1))
      call hipCheck(hipMalloc(self%wbuf4r_gpu,(nx)-(1)+1,(ng)-(1)+1,(nz)-(1)+1,(nv)-(1)+1))
      call hipCheck(hipMalloc(self%wbuf5r_gpu,(nx)-(1)+1,(ny)-(1)+1,(ng)-(1)+1,(nv)-(1)+1))
      call hipCheck(hipMalloc(self%wbuf6r_gpu,(nx)-(1)+1,(ny)-(1)+1,(ng)-(1)+1,(nv)-(1)+1))

      call hipCheck(hipMalloc(self%wbuf1xybcs_gpu,(ng)-(1)+1,(ny+ng)-(1-ng)+1,(nz)-(1)+1,(nv)-(1)+1))
      call hipCheck(hipMalloc(self%wbuf2xybcs_gpu,(ng)-(1)+1,(ny+ng)-(1-ng)+1,(nz)-(1)+1,(nv)-(1)+1))
      call hipCheck(hipMalloc(self%wbuf1xybcr_gpu,(ng)-(1)+1,(ny+ng)-(1-ng)+1,(nz)-(1)+1,(nv)-(1)+1))
      call hipCheck(hipMalloc(self%wbuf2xybcr_gpu,(ng)-(1)+1,(ny+ng)-(1-ng)+1,(nz)-(1)+1,(nv)-(1)+1))
      call hipCheck(hipMalloc(self%x_gpu,(nx+ng)-(1-ng)+1))
      call hipCheck(hipMalloc(self%y_gpu,(ny+ng)-(1-ng)+1))
      call hipCheck(hipMalloc(self%z_gpu,(nz+ng)-(1-ng)+1))

      call hipCheck(hipMalloc(self%yn_gpu,(ny+1)-(1)+1))

      call hipCheck(hipMalloc(self%dcsidx_gpu,(nx)-(1)+1))
      call hipCheck(hipMalloc(self%dcsidx2_gpu,(nx)-(1)+1))
      call hipCheck(hipMalloc(self%dcsidxs_gpu,(nx)-(1)+1))

      call hipCheck(hipMalloc(self%detady_gpu,(ny)-(1)+1))
      call hipCheck(hipMalloc(self%detady2_gpu,(ny)-(1)+1))
      call hipCheck(hipMalloc(self%detadys_gpu,(ny)-(1)+1))

      call hipCheck(hipMalloc(self%dzitdz_gpu,(nz)-(1)+1))
      call hipCheck(hipMalloc(self%dzitdz2_gpu,(nz)-(1)+1))
      call hipCheck(hipMalloc(self%dzitdzs_gpu,(nz)-(1)+1))


      call hipCheck(hipMalloc(self%wbuf1s_c_gpu,(nv)-(1)+1,(ng)-(1)+1,(ny)-(1)+1,(ng)-(1)+1))
      call hipCheck(hipMalloc(self%wbuf2s_c_gpu,(nv)-(1)+1,(ng)-(1)+1,(ny)-(1)+1,(ng)-(1)+1))
      call hipCheck(hipMalloc(self%wbuf3s_c_gpu,(nv)-(1)+1,(ng)-(1)+1,(ny)-(1)+1,(ng)-(1)+1))
      call hipCheck(hipMalloc(self%wbuf4s_c_gpu,(nv)-(1)+1,(ng)-(1)+1,(ny)-(1)+1,(ng)-(1)+1))
      call hipCheck(hipMalloc(self%wbuf1r_c_gpu,(nv)-(1)+1,(ng)-(1)+1,(ny)-(1)+1,(ng)-(1)+1))
      call hipCheck(hipMalloc(self%wbuf2r_c_gpu,(nv)-(1)+1,(ng)-(1)+1,(ny)-(1)+1,(ng)-(1)+1))
      call hipCheck(hipMalloc(self%wbuf3r_c_gpu,(nv)-(1)+1,(ng)-(1)+1,(ny)-(1)+1,(ng)-(1)+1))
      call hipCheck(hipMalloc(self%wbuf4r_c_gpu,(nv)-(1)+1,(ng)-(1)+1,(ny)-(1)+1,(ng)-(1)+1))

      call hipCheck(hipMalloc(self%wbuf_lxly_s_gpu,(nv)-(1)+1,(ng)-(1)+1,(ng)-(1)+1,(nz)-(1)+1))
      call hipCheck(hipMalloc(self%wbuf_lxry_s_gpu,(nv)-(1)+1,(ng)-(1)+1,(ng)-(1)+1,(nz)-(1)+1))

      call hipCheck(hipMalloc(self%wbuf_rxly_s_gpu,(nv)-(1)+1,(ng)-(1)+1,(ng)-(1)+1,(nz)-(1)+1))
      call hipCheck(hipMalloc(self%wbuf_rxry_s_gpu,(nv)-(1)+1,(ng)-(1)+1,(ng)-(1)+1,(nz)-(1)+1))

      call hipCheck(hipMalloc(self%wbuf_lxlz_s_gpu,(nv)-(1)+1,(ng)-(1)+1,(ny)-(1)+1,(ng)-(1)+1))
      call hipCheck(hipMalloc(self%wbuf_lxrz_s_gpu,(nv)-(1)+1,(ng)-(1)+1,(ny)-(1)+1,(ng)-(1)+1))

      call hipCheck(hipMalloc(self%wbuf_rxlz_s_gpu,(nv)-(1)+1,(ng)-(1)+1,(ny)-(1)+1,(ng)-(1)+1))
      call hipCheck(hipMalloc(self%wbuf_rxrz_s_gpu,(nv)-(1)+1,(ng)-(1)+1,(ny)-(1)+1,(ng)-(1)+1))

      call hipCheck(hipMalloc(self%wbuf_lylz_s_gpu,(nv)-(1)+1,(nx)-(1)+1,(ng)-(1)+1,(ng)-(1)+1))
      call hipCheck(hipMalloc(self%wbuf_lyrz_s_gpu,(nv)-(1)+1,(nx)-(1)+1,(ng)-(1)+1,(ng)-(1)+1))

      call hipCheck(hipMalloc(self%wbuf_rylz_s_gpu,(nv)-(1)+1,(nx)-(1)+1,(ng)-(1)+1,(ng)-(1)+1))
      call hipCheck(hipMalloc(self%wbuf_ryrz_s_gpu,(nv)-(1)+1,(nx)-(1)+1,(ng)-(1)+1,(ng)-(1)+1))

      call hipCheck(hipMalloc(self%wbuf_lxlylz_s_gpu,(nv)-(1)+1,(ng)-(1)+1,(ng)-(1)+1,(ng)-(1)+1))
      call hipCheck(hipMalloc(self%wbuf_lxlyrz_s_gpu,(nv)-(1)+1,(ng)-(1)+1,(ng)-(1)+1,(ng)-(1)+1))

      call hipCheck(hipMalloc(self%wbuf_lxrylz_s_gpu,(nv)-(1)+1,(ng)-(1)+1,(ng)-(1)+1,(ng)-(1)+1))
      call hipCheck(hipMalloc(self%wbuf_lxryrz_s_gpu,(nv)-(1)+1,(ng)-(1)+1,(ng)-(1)+1,(ng)-(1)+1))

      call hipCheck(hipMalloc(self%wbuf_rxlylz_s_gpu,(nv)-(1)+1,(ng)-(1)+1,(ng)-(1)+1,(ng)-(1)+1))
      call hipCheck(hipMalloc(self%wbuf_rxlyrz_s_gpu,(nv)-(1)+1,(ng)-(1)+1,(ng)-(1)+1,(ng)-(1)+1))

      call hipCheck(hipMalloc(self%wbuf_rxrylz_s_gpu,(nv)-(1)+1,(ng)-(1)+1,(ng)-(1)+1,(ng)-(1)+1))
      call hipCheck(hipMalloc(self%wbuf_rxryrz_s_gpu,(nv)-(1)+1,(ng)-(1)+1,(ng)-(1)+1,(ng)-(1)+1))


      call hipCheck(hipMalloc(self%wbuf_lxly_r_gpu,(nv)-(1)+1,(ng)-(1)+1,(ng)-(1)+1,(nz)-(1)+1))
      call hipCheck(hipMalloc(self%wbuf_lxry_r_gpu,(nv)-(1)+1,(ng)-(1)+1,(ng)-(1)+1,(nz)-(1)+1))

      call hipCheck(hipMalloc(self%wbuf_rxly_r_gpu,(nv)-(1)+1,(ng)-(1)+1,(ng)-(1)+1,(nz)-(1)+1))
      call hipCheck(hipMalloc(self%wbuf_rxry_r_gpu,(nv)-(1)+1,(ng)-(1)+1,(ng)-(1)+1,(nz)-(1)+1))

      call hipCheck(hipMalloc(self%wbuf_lxlz_r_gpu,(nv)-(1)+1,(ng)-(1)+1,(ny)-(1)+1,(ng)-(1)+1))
      call hipCheck(hipMalloc(self%wbuf_lxrz_r_gpu,(nv)-(1)+1,(ng)-(1)+1,(ny)-(1)+1,(ng)-(1)+1))

      call hipCheck(hipMalloc(self%wbuf_rxlz_r_gpu,(nv)-(1)+1,(ng)-(1)+1,(ny)-(1)+1,(ng)-(1)+1))
      call hipCheck(hipMalloc(self%wbuf_rxrz_r_gpu,(nv)-(1)+1,(ng)-(1)+1,(ny)-(1)+1,(ng)-(1)+1))

      call hipCheck(hipMalloc(self%wbuf_lylz_r_gpu,(nv)-(1)+1,(nx)-(1)+1,(ng)-(1)+1,(ng)-(1)+1))
      call hipCheck(hipMalloc(self%wbuf_lyrz_r_gpu,(nv)-(1)+1,(nx)-(1)+1,(ng)-(1)+1,(ng)-(1)+1))

      call hipCheck(hipMalloc(self%wbuf_rylz_r_gpu,(nv)-(1)+1,(nx)-(1)+1,(ng)-(1)+1,(ng)-(1)+1))
      call hipCheck(hipMalloc(self%wbuf_ryrz_r_gpu,(nv)-(1)+1,(nx)-(1)+1,(ng)-(1)+1,(ng)-(1)+1))

      call hipCheck(hipMalloc(self%wbuf_lxlylz_r_gpu,(nv)-(1)+1,(ng)-(1)+1,(ng)-(1)+1,(ng)-(1)+1))
      call hipCheck(hipMalloc(self%wbuf_lxlyrz_r_gpu,(nv)-(1)+1,(ng)-(1)+1,(ng)-(1)+1,(ng)-(1)+1))

      call hipCheck(hipMalloc(self%wbuf_lxrylz_r_gpu,(nv)-(1)+1,(ng)-(1)+1,(ng)-(1)+1,(ng)-(1)+1))
      call hipCheck(hipMalloc(self%wbuf_lxryrz_r_gpu,(nv)-(1)+1,(ng)-(1)+1,(ng)-(1)+1,(ng)-(1)+1))

      call hipCheck(hipMalloc(self%wbuf_rxlylz_r_gpu,(nv)-(1)+1,(ng)-(1)+1,(ng)-(1)+1,(ng)-(1)+1))
      call hipCheck(hipMalloc(self%wbuf_rxlyrz_r_gpu,(nv)-(1)+1,(ng)-(1)+1,(ng)-(1)+1,(ng)-(1)+1))

      call hipCheck(hipMalloc(self%wbuf_rxrylz_r_gpu,(nv)-(1)+1,(ng)-(1)+1,(ng)-(1)+1,(ng)-(1)+1))
      call hipCheck(hipMalloc(self%wbuf_rxryrz_r_gpu,(nv)-(1)+1,(ng)-(1)+1,(ng)-(1)+1,(ng)-(1)+1))


      if (self%field%grid%grid_dim == 2) then
        call hipCheck(hipMalloc(self%xc2_gpu,(nx+ng)-(1-ng)+1,(ny+ng)-(1-ng)+1))
        call hipCheck(hipMalloc(self%yc2_gpu,(nx+ng)-(1-ng)+1,(ny+ng)-(1-ng)+1))
        call hipCheck(hipMalloc(self%dcsidxc2_gpu,(nx+ng)-(1-ng)+1,(ny+ng)-(1-ng)+1))
        call hipCheck(hipMalloc(self%dcsidyc2_gpu,(nx+ng)-(1-ng)+1,(ny+ng)-(1-ng)+1))
        call hipCheck(hipMalloc(self%detadxc2_gpu,(nx+ng)-(1-ng)+1,(ny+ng)-(1-ng)+1))
        call hipCheck(hipMalloc(self%detadyc2_gpu,(nx+ng)-(1-ng)+1,(ny+ng)-(1-ng)+1))
        call hipCheck(hipMalloc(self%dcsidxnc2_gpu,(nx+ng)-(1-ng)+1,(ny+ng)-(1-ng)+1))
        call hipCheck(hipMalloc(self%dcsidync2_gpu,(nx+ng)-(1-ng)+1,(ny+ng)-(1-ng)+1))
        call hipCheck(hipMalloc(self%detadxnc2_gpu,(nx+ng)-(1-ng)+1,(ny+ng)-(1-ng)+1))
        call hipCheck(hipMalloc(self%detadync2_gpu,(nx+ng)-(1-ng)+1,(ny+ng)-(1-ng)+1))
        call hipCheck(hipMalloc(self%dxdcsic2_gpu,(nx+ng)-(1-ng)+1,(ny+ng)-(1-ng)+1))
        call hipCheck(hipMalloc(self%dydcsic2_gpu,(nx+ng)-(1-ng)+1,(ny+ng)-(1-ng)+1))
        call hipCheck(hipMalloc(self%dxdetac2_gpu,(nx+ng)-(1-ng)+1,(ny+ng)-(1-ng)+1))
        call hipCheck(hipMalloc(self%dydetac2_gpu,(nx+ng)-(1-ng)+1,(ny+ng)-(1-ng)+1))
        call hipCheck(hipMalloc(self%dxdcsinc2_gpu,(nx+ng)-(1-ng)+1,(ny+ng)-(1-ng)+1))
        call hipCheck(hipMalloc(self%dydcsinc2_gpu,(nx+ng)-(1-ng)+1,(ny+ng)-(1-ng)+1))
        call hipCheck(hipMalloc(self%dxdetanc2_gpu,(nx+ng)-(1-ng)+1,(ny+ng)-(1-ng)+1))
        call hipCheck(hipMalloc(self%dydetanc2_gpu,(nx+ng)-(1-ng)+1,(ny+ng)-(1-ng)+1))
        call hipCheck(hipMalloc(self%mcsi_gpu,(nx+ng)-(1-ng)+1,(ny+ng)-(1-ng)+1))
        call hipCheck(hipMalloc(self%meta_gpu,(nx+ng)-(1-ng)+1,(ny+ng)-(1-ng)+1))
        call hipCheck(hipMalloc(self%mcsijac1_gpu,(nx+ng)-(1-ng)+1,(ny+ng)-(1-ng)+1))
        call hipCheck(hipMalloc(self%metajac1_gpu,(nx+ng)-(1-ng)+1,(ny+ng)-(1-ng)+1))
        call hipCheck(hipMalloc(self%jac_gpu,(nx+ng)-(1-ng)+1,(ny+ng)-(1-ng)+1))
        call hipCheck(hipMalloc(self%g1_gpu,(nx+ng)-(1-ng)+1,(ny+ng)-(1-ng)+1))
        call hipCheck(hipMalloc(self%g2_gpu,(nx+ng)-(1-ng)+1,(ny+ng)-(1-ng)+1))
        call hipCheck(hipMalloc(self%g12_gpu,(nx+ng)-(1-ng)+1,(ny+ng)-(1-ng)+1))
        call hipCheck(hipMalloc(self%csimod_gpu,(nx+ng)-(1-ng)+1,(ny+ng)-(1-ng)+1))
        call hipCheck(hipMalloc(self%etamod_gpu,(nx+ng)-(1-ng)+1,(ny+ng)-(1-ng)+1))
        call hipCheck(hipMalloc(self%theta_ij_gpu,(nx)-(1)+1,(ny)-(1)+1))

        call hipCheck(hipMalloc(self%wall_tag_gpu,(nx+ng)-(1-ng)+1))
        call hipCheck(hipMalloc(self%wbuftus_gpu,(nz)-(1)+1,(2)-(1)+1))
        call hipCheck(hipMalloc(self%wbuftes_gpu,(nz)-(1)+1,(2)-(1)+1))
        call hipCheck(hipMalloc(self%wbuftur_gpu,(nz)-(1)+1,(2)-(1)+1))
        call hipCheck(hipMalloc(self%wbufter_gpu,(nz)-(1)+1,(2)-(1)+1))
      endif

    endassociate

    call hipCheck(hipMemcpy(self%x_gpu , self%field%x,hipMemcpyHostToDevice))
    call hipCheck(hipMemcpy(self%y_gpu , self%field%y,hipMemcpyHostToDevice))
    call hipCheck(hipMemcpy(self%z_gpu , self%field%z,hipMemcpyHostToDevice))
    call hipCheck(hipMemcpy(self%yn_gpu , self%field%yn,hipMemcpyHostToDevice))

    call hipCheck(hipMemcpy(self%dcsidx_gpu , self%field%dcsidx,hipMemcpyHostToDevice))
    call hipCheck(hipMemcpy(self%dcsidxs_gpu , self%field%dcsidxs,hipMemcpyHostToDevice))
    call hipCheck(hipMemcpy(self%dcsidx2_gpu , self%field%dcsidx2,hipMemcpyHostToDevice))
    call hipCheck(hipMemcpy(self%detady_gpu , self%field%detady,hipMemcpyHostToDevice))
    call hipCheck(hipMemcpy(self%detadys_gpu , self%field%detadys,hipMemcpyHostToDevice))
    call hipCheck(hipMemcpy(self%detady2_gpu , self%field%detady2,hipMemcpyHostToDevice))
    call hipCheck(hipMemcpy(self%dzitdz_gpu , self%field%dzitdz,hipMemcpyHostToDevice))
    call hipCheck(hipMemcpy(self%dzitdzs_gpu , self%field%dzitdzs,hipMemcpyHostToDevice))
    call hipCheck(hipMemcpy(self%dzitdz2_gpu , self%field%dzitdz2,hipMemcpyHostToDevice))

    if (self%field%grid%grid_dim == 2) then
      call hipCheck(hipMemcpy(self%xc2_gpu , self%field%xc2,hipMemcpyHostToDevice))
      call hipCheck(hipMemcpy(self%yc2_gpu , self%field%yc2,hipMemcpyHostToDevice))
      call hipCheck(hipMemcpy(self%dcsidxc2_gpu , self%field%dcsidxc2,hipMemcpyHostToDevice))
      call hipCheck(hipMemcpy(self%dcsidyc2_gpu , self%field%dcsidyc2,hipMemcpyHostToDevice))
      call hipCheck(hipMemcpy(self%detadxc2_gpu , self%field%detadxc2,hipMemcpyHostToDevice))
      call hipCheck(hipMemcpy(self%detadyc2_gpu , self%field%detadyc2,hipMemcpyHostToDevice))
      call hipCheck(hipMemcpy(self%dcsidxnc2_gpu , self%field%dcsidxnc2,hipMemcpyHostToDevice))
      call hipCheck(hipMemcpy(self%dcsidync2_gpu , self%field%dcsidync2,hipMemcpyHostToDevice))
      call hipCheck(hipMemcpy(self%detadxnc2_gpu , self%field%detadxnc2,hipMemcpyHostToDevice))
      call hipCheck(hipMemcpy(self%detadync2_gpu , self%field%detadync2,hipMemcpyHostToDevice))
      call hipCheck(hipMemcpy(self%dxdcsic2_gpu , self%field%dxdcsic2,hipMemcpyHostToDevice))
      call hipCheck(hipMemcpy(self%dydcsic2_gpu , self%field%dydcsic2,hipMemcpyHostToDevice))
      call hipCheck(hipMemcpy(self%dxdetac2_gpu , self%field%dxdetac2,hipMemcpyHostToDevice))
      call hipCheck(hipMemcpy(self%dydetac2_gpu , self%field%dydetac2,hipMemcpyHostToDevice))
      call hipCheck(hipMemcpy(self%dxdcsinc2_gpu , self%field%dxdcsinc2,hipMemcpyHostToDevice))
      call hipCheck(hipMemcpy(self%dydcsinc2_gpu , self%field%dydcsinc2,hipMemcpyHostToDevice))
      call hipCheck(hipMemcpy(self%dxdetanc2_gpu , self%field%dxdetanc2,hipMemcpyHostToDevice))
      call hipCheck(hipMemcpy(self%dydetanc2_gpu , self%field%dydetanc2,hipMemcpyHostToDevice))
      call hipCheck(hipMemcpy(self%mcsi_gpu , self%field%mcsi,hipMemcpyHostToDevice))
      call hipCheck(hipMemcpy(self%meta_gpu , self%field%meta,hipMemcpyHostToDevice))
      call hipCheck(hipMemcpy(self%mcsijac1_gpu , self%field%mcsijac1,hipMemcpyHostToDevice))
      call hipCheck(hipMemcpy(self%metajac1_gpu , self%field%metajac1,hipMemcpyHostToDevice))
      call hipCheck(hipMemcpy(self%jac_gpu , self%field%jac,hipMemcpyHostToDevice))
      call hipCheck(hipMemcpy(self%g1_gpu , self%field%g1,hipMemcpyHostToDevice))
      call hipCheck(hipMemcpy(self%g2_gpu , self%field%g2,hipMemcpyHostToDevice))
      call hipCheck(hipMemcpy(self%g12_gpu , self%field%g12,hipMemcpyHostToDevice))
      call hipCheck(hipMemcpy(self%csimod_gpu , self%field%csimod,hipMemcpyHostToDevice))
      call hipCheck(hipMemcpy(self%etamod_gpu , self%field%etamod,hipMemcpyHostToDevice))
      call hipCheck(hipMemcpy(self%wall_tag_gpu , self%field%wall_tag,hipMemcpyHostToDevice))
      call hipCheck(hipMemcpy(self%theta_ij_gpu , self%field%theta_ij,hipMemcpyHostToDevice))
    endif

  endsubroutine alloc

  subroutine bcte_step_1_kernel(nx,ny,nz,nv,ng,nrank_x,ite_rank_x,itu_rank_x,ite_l,itu_l,w_gpu,wbuftus_gpu,wbuftes_gpu)
    integer :: nx,ny,nz
    integer :: nv,ng,nrank_x
    integer :: ite_rank_x,itu_rank_x,ite_l
    integer :: itu_l
    real(rkind), dimension(:,:,:,:), target :: w_gpu
    real(rkind), dimension(:,:), target :: wbuftus_gpu
    real(rkind), dimension(:,:), target :: wbuftes_gpu

    if(nrank_x == ite_rank_x) then
      call bcte_step_1_kernel1_wrapper(c_null_ptr,nx,ny,nz,nv,ng,nrank_x,ite_rank_x,itu_rank_x,&
      &ite_l,itu_l,c_loc(w_gpu),c_loc(wbuftus_gpu),c_loc(wbuftes_gpu))
      call hipCheck(hipDeviceSynchronize())
    endif
    if(nrank_x == itu_rank_x) then
      call bcte_step_1_kernel2_wrapper(c_null_ptr,nx,ny,nz,nv,ng,nrank_x,ite_rank_x,itu_rank_x,&
      &ite_l,itu_l,c_loc(w_gpu),c_loc(wbuftus_gpu),c_loc(wbuftes_gpu))
      call hipCheck(hipDeviceSynchronize())
    endif

  endsubroutine bcte_step_1_kernel


  subroutine bcte(self, steps)
    class(base_amd_object), intent(inout) :: self
    logical, dimension(3), optional :: steps
    logical, dimension(3) :: steps_
    integer, dimension(4) :: requests

    steps_ = .true. ; if(present(steps)) steps_ = steps
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv,&
    & w_gpu => self%w_gpu, wbuftes_gpu => self%wbuftes_gpu, wbuftus_gpu => self%wbuftus_gpu,&
    & wbufter_gpu => self%wbufter_gpu, wbuftur_gpu => self%wbuftur_gpu, ite_rank_x => self%field%ite_rank&
    &_x, itu_rank_x => self%field%itu_rank_x, ite_l => self%field%ite_l, itu_l => self%field%itu_l,&
    & nrank_x => self%field%nrank_x, mp_cartx => self%field%mp_cartx, iermpi => self%mpi_err)

      if(steps_(1)) then
        call bcte_step_1_kernel(nx,ny,nz,nv,ng,nrank_x,ite_rank_x,itu_rank_x,ite_l,itu_l,w_gpu,wbuftus_gpu,wbuftes_gpu)
      endif

      if(steps_(2)) then
        if (nrank_x == ite_rank_x .and. nrank_x == itu_rank_x) then
          call hipCheck(hipMemcpy(wbufter_gpu , wbuftes_gpu,hipMemcpyDeviceToDevice))
          call hipCheck(hipMemcpy(wbuftur_gpu , wbuftus_gpu,hipMemcpyDeviceToDevice))
        else
          requests(:) = mpi_request_null
          if (nrank_x == ite_rank_x) then
            call mpi_isend(wbuftes_gpu,2*nz,mpi_prec,itu_rank_x,1,mp_cartx,requests(1),iermpi)
            call mpi_irecv(wbuftur_gpu,2*nz,mpi_prec,itu_rank_x,1,mp_cartx,requests(2),iermpi)
          elseif (nrank_x == itu_rank_x) then
            call mpi_isend(wbuftus_gpu,2*nz,mpi_prec,ite_rank_x,1,mp_cartx,requests(3),iermpi)
            call mpi_irecv(wbufter_gpu,2*nz,mpi_prec,ite_rank_x,1,mp_cartx,requests(4),iermpi)
          endif
          call mpi_waitall(4,requests,mpi_statuses_ignore,iermpi)
        endif
      endif

      if(steps_(3)) then
        call bcte_step_3_kernel(nx,ny,nz,nv,ng,nrank_x,ite_rank_x,itu_rank_x,ite_l,itu_l,w_gpu,wbuftur_gpu,wbufter_gpu)
      endif
    endassociate
  endsubroutine bcte

  subroutine bcte_var(self, w_swap, steps)
    class(base_amd_object), intent(inout) :: self
    real(rkind), dimension(:,:,:,:), target :: w_swap
    logical, dimension(3), optional :: steps
    logical, dimension(3) :: steps_
    integer, dimension(4) :: requests

    steps_ = .true. ; if(present(steps)) steps_ = steps
    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv,&
    & w_gpu => self%w_gpu, wbuftes_gpu => self%wbuftes_gpu, wbuftus_gpu => self%wbuftus_gpu,&
    & wbufter_gpu => self%wbufter_gpu, wbuftur_gpu => self%wbuftur_gpu, ite_rank_x => self%field%ite_rank&
    &_x, itu_rank_x => self%field%itu_rank_x, ite_l => self%field%ite_l, itu_l => self%field%itu_l,&
    & nrank_x => self%field%nrank_x, mp_cartx => self%field%mp_cartx, iermpi => self%mpi_err)

      if(steps_(1)) then
        call bcte_step_1_kernel(nx,ny,nz,1,ng,nrank_x,ite_rank_x,itu_rank_x,ite_l,itu_l,w_swap,wbuftus_gpu,wbuftes_gpu)
      endif

      if(steps_(2)) then
        if (nrank_x == ite_rank_x .and. nrank_x == itu_rank_x) then
          call hipCheck(hipMemcpy(wbufter_gpu , wbuftes_gpu,hipMemcpyDeviceToDevice))
          call hipCheck(hipMemcpy(wbuftur_gpu , wbuftus_gpu,hipMemcpyDeviceToDevice))
        else
          requests(:) = mpi_request_null
          if (nrank_x == ite_rank_x) then
            call mpi_isend(wbuftes_gpu,1*nz,mpi_prec,itu_rank_x,3,mp_cartx,requests(1),iermpi)
            call mpi_irecv(wbuftur_gpu,1*nz,mpi_prec,itu_rank_x,3,mp_cartx,requests(2),iermpi)
          elseif (nrank_x == itu_rank_x) then
            call mpi_isend(wbuftus_gpu,1*nz,mpi_prec,ite_rank_x,3,mp_cartx,requests(3),iermpi)
            call mpi_irecv(wbufter_gpu,1*nz,mpi_prec,ite_rank_x,3,mp_cartx,requests(4),iermpi)
          endif
          call mpi_waitall(4,requests,mpi_statuses_ignore,iermpi)
        endif
      endif

      if(steps_(3)) then
        call bcte_step_3_kernel(nx,ny,nz,1,ng,nrank_x,ite_rank_x,itu_rank_x,ite_l,itu_l,w_swap,wbuftur_gpu,wbufter_gpu)
      endif
    endassociate
  endsubroutine bcte_var

  subroutine bcte_step_3_kernel(nx,ny,nz,nv,ng,nrank_x,ite_rank_x,itu_rank_x,ite_l,itu_l,w_gpu,wbuftur_gpu,wbufter_gpu)
    integer :: nx,ny,nz
    integer :: nv,ng,nrank_x
    integer :: ite_rank_x,itu_rank_x,ite_l
    integer :: itu_l
    real(rkind), dimension(:,:,:,:), target :: w_gpu
    real(rkind), dimension(:,:), target :: wbuftur_gpu
    real(rkind), dimension(:,:), target :: wbufter_gpu

    if(nrank_x == ite_rank_x) then
      call bcte_step_3_kernel1_wrapper(c_null_ptr,nx,ny,nz,nv,ng,nrank_x,ite_rank_x,itu_rank_x,&
      &ite_l,itu_l,c_loc(w_gpu),c_loc(wbuftur_gpu),c_loc(wbufter_gpu))
      call hipCheck(hipDeviceSynchronize())
    endif

    if (nrank_x==itu_rank_x) then
      call bcte_step_3_kernel2_wrapper(c_null_ptr,nx,ny,nz,nv,ng,nrank_x,ite_rank_x,itu_rank_x,&
      &ite_l,itu_l,c_loc(w_gpu),c_loc(wbuftur_gpu),c_loc(wbufter_gpu))
      call hipCheck(hipDeviceSynchronize())
    endif

  endsubroutine bcte_step_3_kernel


  subroutine bcswap_wake_step_1_kernel(nx,ny,nz,nv,ng,w_gpu,wbuf4s_gpu)
    integer :: nx,ny,nz
    integer :: ng,nv
    real(rkind), dimension(:,:,:,:), target :: w_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf4s_gpu


    call bcswap_wake_step_1_kernel_wrapper(c_null_ptr,nx,ny,nz,ng,nv,c_loc(w_gpu),c_loc(wbuf4s_gpu))
    call hipCheck(hipDeviceSynchronize())


  endsubroutine bcswap_wake_step_1_kernel


  subroutine bcswap_wake_step_3_kernel(nx,ny,nz,nv,ng,w_gpu,wbuf3r_gpu,wall_tag_gpu)
    integer :: nx,ny,nz
    integer :: ng,nv
    real(rkind), dimension(:,:,:,:), target :: w_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf3r_gpu
    integer, dimension(:), target :: wall_tag_gpu


    call bcswap_wake_step_3_kernel_wrapper(c_null_ptr,nx,ny,nz,ng,nv,c_loc(wall_tag_gpu),c_loc(w_gpu),c_loc(wbuf3r_gpu))
    call hipCheck(hipDeviceSynchronize())


  endsubroutine bcswap_wake_step_3_kernel


  subroutine bcswap_step_1_kernel(nx,ny,nz,nv,ng,ndim,w_gpu,wbuf1s_gpu,wbuf2s_gpu,wbuf3s_gpu,wbuf4s_gpu,wbuf5s_gpu,wbuf6s_gpu)
    integer :: nx,ny,nz
    integer :: ng,nv,ndim
    real(rkind), dimension(:,:,:,:), target :: w_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf1s_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf2s_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf3s_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf4s_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf5s_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf6s_gpu


    call bcswap_step_1_kernel1_wrapper(c_null_ptr,nx,ny,nz,ng,nv,ndim,c_loc(w_gpu),&
    &c_loc(wbuf1s_gpu),c_loc(wbuf2s_gpu),c_loc(wbuf3s_gpu),c_loc(wbuf4s_gpu),c_loc(wbuf5s_gpu),&
    &c_loc(wbuf6s_gpu))
    call hipCheck(hipDeviceSynchronize())

    call bcswap_step_1_kernel2_wrapper(c_null_ptr,nx,ny,nz,ng,nv,ndim,c_loc(w_gpu),&
    &c_loc(wbuf1s_gpu),c_loc(wbuf2s_gpu),c_loc(wbuf3s_gpu),c_loc(wbuf4s_gpu),c_loc(wbuf5s_gpu),&
    &c_loc(wbuf6s_gpu))
    call hipCheck(hipDeviceSynchronize())
    if (ndim==3) then
      call bcswap_step_1_kernel3_wrapper(c_null_ptr,nx,ny,nz,ng,nv,ndim,c_loc(w_gpu),&
      &c_loc(wbuf1s_gpu),c_loc(wbuf2s_gpu),c_loc(wbuf3s_gpu),c_loc(wbuf4s_gpu),c_loc(wbuf5s_gpu),&
      &c_loc(wbuf6s_gpu))
      call hipCheck(hipDeviceSynchronize())
    endif

  endsubroutine bcswap_step_1_kernel


  subroutine bcswap_c2_step_1_kernel(nx,ny,nz,nv,ng,ndim,w_gpu,wbuf1s_gpu,wbuf2s_gpu,wbuf3s_gpu,&
  &wbuf4s_gpu,wbuf5s_gpu,wbuf6s_gpu,dcsidxnc2_gpu,dcsidync2_gpu,detadxnc2_gpu,detadync2_gpu,&
  &is_periodic)
    integer :: nx,ny,nz
    integer :: ng,nv,ndim
    real(rkind), dimension(:,:,:,:), target :: w_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf1s_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf2s_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf3s_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf4s_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf5s_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf6s_gpu
    real(rkind), dimension(:,:), target :: dcsidxnc2_gpu
    real(rkind), dimension(:,:), target :: dcsidync2_gpu
    real(rkind), dimension(:,:), target :: detadxnc2_gpu
    real(rkind), dimension(:,:), target :: detadync2_gpu
    integer, dimension(3) :: is_periodic

    if (is_periodic(1) == 1) then
      call bcswap_c2_step_1_kernel1_wrapper(c_null_ptr,nx,ny,nz,ng,nv,ndim,c_loc(w_gpu),&
      &c_loc(wbuf1s_gpu),c_loc(wbuf2s_gpu),c_loc(wbuf3s_gpu),c_loc(wbuf4s_gpu),c_loc(wbuf5s_gpu),&
      &c_loc(wbuf6s_gpu),c_loc(dcsidxnc2_gpu),c_loc(dcsidync2_gpu),c_loc(detadxnc2_gpu),&
      &c_loc(detadync2_gpu))
      call hipCheck(hipDeviceSynchronize())
    else
      call bcswap_c2_step_1_kernel2_wrapper(c_null_ptr,nx,ny,nz,ng,nv,ndim,c_loc(w_gpu),&
      &c_loc(wbuf1s_gpu),c_loc(wbuf2s_gpu),c_loc(wbuf3s_gpu),c_loc(wbuf4s_gpu),c_loc(wbuf5s_gpu),&
      &c_loc(wbuf6s_gpu),c_loc(dcsidxnc2_gpu),c_loc(dcsidync2_gpu),c_loc(detadxnc2_gpu),&
      &c_loc(detadync2_gpu))
      call hipCheck(hipDeviceSynchronize())
    endif
    call bcswap_c2_step_1_kernel3_wrapper(c_null_ptr,nx,ny,nz,ng,nv,ndim,c_loc(w_gpu),&
    &c_loc(wbuf1s_gpu),c_loc(wbuf2s_gpu),c_loc(wbuf3s_gpu),c_loc(wbuf4s_gpu),c_loc(wbuf5s_gpu),&
    &c_loc(wbuf6s_gpu),c_loc(dcsidxnc2_gpu),c_loc(dcsidync2_gpu),c_loc(detadxnc2_gpu),&
    &c_loc(detadync2_gpu))
    call hipCheck(hipDeviceSynchronize())
    if (ndim==3) then
      call bcswap_c2_step_1_kernel4_wrapper(c_null_ptr,nx,ny,nz,ng,nv,ndim,c_loc(w_gpu),&
      &c_loc(wbuf1s_gpu),c_loc(wbuf2s_gpu),c_loc(wbuf3s_gpu),c_loc(wbuf4s_gpu),c_loc(wbuf5s_gpu),&
      &c_loc(wbuf6s_gpu),c_loc(dcsidxnc2_gpu),c_loc(dcsidync2_gpu),c_loc(detadxnc2_gpu),&
      &c_loc(detadync2_gpu))
      call hipCheck(hipDeviceSynchronize())
    endif

  endsubroutine bcswap_c2_step_1_kernel


  subroutine bcswap_c2xybc_step_1_kernel(nx,ny,nz,nv,ng,ndim,w_gpu,wbuf1xybcs_gpu,wbuf2xybcs_gpu,&
  &wbuf3s_gpu,wbuf4s_gpu,wbuf5s_gpu,wbuf6s_gpu,dcsidxnc2_gpu,dcsidync2_gpu,detadxnc2_gpu,detadync2_gpu,&
  &is_periodic)
    integer :: nx,ny,nz
    integer :: ng,nv,ndim
    real(rkind), dimension(:,:,:,:), target :: w_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf1xybcs_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf2xybcs_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf3s_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf4s_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf5s_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf6s_gpu
    real(rkind), dimension(:,:), target :: dcsidxnc2_gpu
    real(rkind), dimension(:,:), target :: dcsidync2_gpu
    real(rkind), dimension(:,:), target :: detadxnc2_gpu
    real(rkind), dimension(:,:), target :: detadync2_gpu
    integer, dimension(3) :: is_periodic

    if (is_periodic(1) == 1) then
      call bcswap_c2xybc_step_1_kernel1_wrapper(c_null_ptr,nx,ny,nz,ng,nv,ndim,c_loc(w_gpu),&
      &c_loc(wbuf1xybcs_gpu),c_loc(wbuf2xybcs_gpu),c_loc(wbuf3s_gpu),c_loc(wbuf4s_gpu),c_loc(wbuf5s_gpu),&
      &c_loc(wbuf6s_gpu),c_loc(dcsidxnc2_gpu),c_loc(dcsidync2_gpu),c_loc(detadxnc2_gpu),&
      &c_loc(detadync2_gpu))
      call hipCheck(hipDeviceSynchronize())
    else
      call bcswap_c2xybc_step_1_kernel2_wrapper(c_null_ptr,nx,ny,nz,ng,nv,ndim,c_loc(w_gpu),&
      &c_loc(wbuf1xybcs_gpu),c_loc(wbuf2xybcs_gpu),c_loc(wbuf3s_gpu),c_loc(wbuf4s_gpu),c_loc(wbuf5s_gpu),&
      &c_loc(wbuf6s_gpu),c_loc(dcsidxnc2_gpu),c_loc(dcsidync2_gpu),c_loc(detadxnc2_gpu),&
      &c_loc(detadync2_gpu))
      call hipCheck(hipDeviceSynchronize())
    endif
    call bcswap_c2xybc_step_1_kernel3_wrapper(c_null_ptr,nx,ny,nz,ng,nv,ndim,c_loc(w_gpu),&
    &c_loc(wbuf1xybcs_gpu),c_loc(wbuf2xybcs_gpu),c_loc(wbuf3s_gpu),c_loc(wbuf4s_gpu),c_loc(wbuf5s_gpu),&
    &c_loc(wbuf6s_gpu),c_loc(dcsidxnc2_gpu),c_loc(dcsidync2_gpu),c_loc(detadxnc2_gpu),&
    &c_loc(detadync2_gpu))
    call hipCheck(hipDeviceSynchronize())
    if (ndim==3) then
      call bcswap_c2xybc_step_1_kernel4_wrapper(c_null_ptr,nx,ny,nz,ng,nv,ndim,c_loc(w_gpu),&
      &c_loc(wbuf1xybcs_gpu),c_loc(wbuf2xybcs_gpu),c_loc(wbuf3s_gpu),c_loc(wbuf4s_gpu),c_loc(wbuf5s_gpu),&
      &c_loc(wbuf6s_gpu),c_loc(dcsidxnc2_gpu),c_loc(dcsidync2_gpu),c_loc(detadxnc2_gpu),&
      &c_loc(detadync2_gpu))
      call hipCheck(hipDeviceSynchronize())
    endif

  endsubroutine bcswap_c2xybc_step_1_kernel


  subroutine bcswap_corner_step_1_kernel(nx,ny,nz,nv,ng,w_gpu,wbuf1s_c_gpu,wbuf2s_c_gpu,wbuf3s_c_gpu,wbuf4s_c_gpu)
    integer :: nx,ny,nz
    integer :: ng,nv
    real(rkind), dimension(:,:,:,:), target :: w_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf1s_c_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf2s_c_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf3s_c_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf4s_c_gpu


    call bcswap_corner_step_1_kernel_wrapper(c_null_ptr,nx,ny,nz,ng,nv,c_loc(w_gpu),&
    &c_loc(wbuf1s_c_gpu),c_loc(wbuf2s_c_gpu),c_loc(wbuf3s_c_gpu),c_loc(wbuf4s_c_gpu))
    call hipCheck(hipDeviceSynchronize())


  endsubroutine bcswap_corner_step_1_kernel


  subroutine bcswap_corner_step_1b_kernel(nx,ny,nz,nv,ng,w_gpu,wbuf_lxly_s_gpu,wbuf_lxry_s_gpu,&
  &wbuf_rxly_s_gpu,wbuf_rxry_s_gpu,wbuf_lylz_s_gpu,wbuf_lyrz_s_gpu,wbuf_rylz_s_gpu,wbuf_ryrz_s_gpu,&
  &wbuf_lxlylz_s_gpu,wbuf_lxlyrz_s_gpu,wbuf_lxrylz_s_gpu,wbuf_lxryrz_s_gpu,wbuf_rxlylz_s_gpu,&
  &wbuf_rxlyrz_s_gpu,wbuf_rxrylz_s_gpu,wbuf_rxryrz_s_gpu)
    integer :: nx,ny,nz
    integer :: ng,nv
    real(rkind), dimension(:,:,:,:), target :: w_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf_lxly_s_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf_lxry_s_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf_rxly_s_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf_rxry_s_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf_lylz_s_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf_lyrz_s_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf_rylz_s_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf_ryrz_s_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf_lxlylz_s_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf_lxlyrz_s_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf_lxrylz_s_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf_lxryrz_s_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf_rxlylz_s_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf_rxlyrz_s_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf_rxrylz_s_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf_rxryrz_s_gpu


    call bcswap_corner_step_1b_kernel1_wrapper(c_null_ptr,nx,ny,nz,ng,nv,c_loc(w_gpu),&
    &c_loc(wbuf_lxly_s_gpu),c_loc(wbuf_lxry_s_gpu),c_loc(wbuf_rxly_s_gpu),c_loc(wbuf_rxry_s_gpu),&
    &c_loc(wbuf_lylz_s_gpu),c_loc(wbuf_lyrz_s_gpu),c_loc(wbuf_rylz_s_gpu),c_loc(wbuf_ryrz_s_gpu),&
    &c_loc(wbuf_lxlylz_s_gpu),c_loc(wbuf_lxlyrz_s_gpu),c_loc(wbuf_lxrylz_s_gpu),c_loc(wbuf_lxryrz_s_gpu),&
    &c_loc(wbuf_rxlylz_s_gpu),c_loc(wbuf_rxlyrz_s_gpu),c_loc(wbuf_rxrylz_s_gpu),c_loc(wbuf_rxryrz_s_gpu))
    call hipCheck(hipDeviceSynchronize())

    call bcswap_corner_step_1b_kernel2_wrapper(c_null_ptr,nx,ny,nz,ng,nv,c_loc(w_gpu),&
    &c_loc(wbuf_lxly_s_gpu),c_loc(wbuf_lxry_s_gpu),c_loc(wbuf_rxly_s_gpu),c_loc(wbuf_rxry_s_gpu),&
    &c_loc(wbuf_lylz_s_gpu),c_loc(wbuf_lyrz_s_gpu),c_loc(wbuf_rylz_s_gpu),c_loc(wbuf_ryrz_s_gpu),&
    &c_loc(wbuf_lxlylz_s_gpu),c_loc(wbuf_lxlyrz_s_gpu),c_loc(wbuf_lxrylz_s_gpu),c_loc(wbuf_lxryrz_s_gpu),&
    &c_loc(wbuf_rxlylz_s_gpu),c_loc(wbuf_rxlyrz_s_gpu),c_loc(wbuf_rxrylz_s_gpu),c_loc(wbuf_rxryrz_s_gpu))
    call hipCheck(hipDeviceSynchronize())

    call bcswap_corner_step_1b_kernel3_wrapper(c_null_ptr,nx,ny,nz,ng,nv,c_loc(w_gpu),&
    &c_loc(wbuf_lxly_s_gpu),c_loc(wbuf_lxry_s_gpu),c_loc(wbuf_rxly_s_gpu),c_loc(wbuf_rxry_s_gpu),&
    &c_loc(wbuf_lylz_s_gpu),c_loc(wbuf_lyrz_s_gpu),c_loc(wbuf_rylz_s_gpu),c_loc(wbuf_ryrz_s_gpu),&
    &c_loc(wbuf_lxlylz_s_gpu),c_loc(wbuf_lxlyrz_s_gpu),c_loc(wbuf_lxrylz_s_gpu),c_loc(wbuf_lxryrz_s_gpu),&
    &c_loc(wbuf_rxlylz_s_gpu),c_loc(wbuf_rxlyrz_s_gpu),c_loc(wbuf_rxrylz_s_gpu),c_loc(wbuf_rxryrz_s_gpu))
    call hipCheck(hipDeviceSynchronize())


  endsubroutine bcswap_corner_step_1b_kernel


  subroutine bcswap_step_3_kernel(nx,ny,nz,nv,ng,ndim,w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,&
  &wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,ileftx,ilefty,ileftz,irightx,irighty,irightz)
    integer :: nx,ny,nz
    integer :: ng,nv,ndim
    integer :: ileftx,ilefty,ileftz
    integer :: irightx,irighty,irightz
    real(rkind), dimension(:,:,:,:), target :: w_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf1r_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf2r_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf3r_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf4r_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf5r_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf6r_gpu

    if (ileftx/=mpi_proc_null) then
      call bcswap_step_3_kernel1_wrapper(c_null_ptr,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,&
      &irightx,irighty,irightz,c_loc(w_gpu),c_loc(wbuf1r_gpu),c_loc(wbuf2r_gpu),c_loc(wbuf3r_gpu),&
      &c_loc(wbuf4r_gpu),c_loc(wbuf5r_gpu),c_loc(wbuf6r_gpu))
      call hipCheck(hipDeviceSynchronize())
    endif
    if (irightx/=mpi_proc_null) then
      call bcswap_step_3_kernel2_wrapper(c_null_ptr,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,&
      &irightx,irighty,irightz,c_loc(w_gpu),c_loc(wbuf1r_gpu),c_loc(wbuf2r_gpu),c_loc(wbuf3r_gpu),&
      &c_loc(wbuf4r_gpu),c_loc(wbuf5r_gpu),c_loc(wbuf6r_gpu))
      call hipCheck(hipDeviceSynchronize())
    endif
    if (ilefty/=mpi_proc_null) then
      call bcswap_step_3_kernel3_wrapper(c_null_ptr,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,&
      &irightx,irighty,irightz,c_loc(w_gpu),c_loc(wbuf1r_gpu),c_loc(wbuf2r_gpu),c_loc(wbuf3r_gpu),&
      &c_loc(wbuf4r_gpu),c_loc(wbuf5r_gpu),c_loc(wbuf6r_gpu))
      call hipCheck(hipDeviceSynchronize())
    endif
    if (irighty/=mpi_proc_null) then
      call bcswap_step_3_kernel4_wrapper(c_null_ptr,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,&
      &irightx,irighty,irightz,c_loc(w_gpu),c_loc(wbuf1r_gpu),c_loc(wbuf2r_gpu),c_loc(wbuf3r_gpu),&
      &c_loc(wbuf4r_gpu),c_loc(wbuf5r_gpu),c_loc(wbuf6r_gpu))
      call hipCheck(hipDeviceSynchronize())
    endif
    if(ndim == 3) then
      if (ileftz/=mpi_proc_null) then
        call bcswap_step_3_kernel5_wrapper(c_null_ptr,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,&
        &irightx,irighty,irightz,c_loc(w_gpu),c_loc(wbuf1r_gpu),c_loc(wbuf2r_gpu),c_loc(wbuf3r_gpu),&
        &c_loc(wbuf4r_gpu),c_loc(wbuf5r_gpu),c_loc(wbuf6r_gpu))
        call hipCheck(hipDeviceSynchronize())
      endif
      if (irightz/=mpi_proc_null) then
        call bcswap_step_3_kernel6_wrapper(c_null_ptr,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,&
        &irightx,irighty,irightz,c_loc(w_gpu),c_loc(wbuf1r_gpu),c_loc(wbuf2r_gpu),c_loc(wbuf3r_gpu),&
        &c_loc(wbuf4r_gpu),c_loc(wbuf5r_gpu),c_loc(wbuf6r_gpu))
        call hipCheck(hipDeviceSynchronize())
      endif
    endif

  endsubroutine bcswap_step_3_kernel


  subroutine bcswap_c2_step_3_kernel(nx,ny,nz,nv,ng,ndim,w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,&
  &wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,ileftx,ilefty,ileftz,irightx,irighty,irightz,dxdcsinc2_gpu,&
  &dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu,is_periodic)
    integer :: nx,ny,nz
    integer :: ng,nv,ndim
    integer :: ileftx,ilefty,ileftz
    integer :: irightx,irighty,irightz
    real(rkind), dimension(:,:,:,:), target :: w_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf1r_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf2r_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf3r_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf4r_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf5r_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf6r_gpu
    real(rkind), dimension(:,:), target :: dxdcsinc2_gpu
    real(rkind), dimension(:,:), target :: dydcsinc2_gpu
    real(rkind), dimension(:,:), target :: dxdetanc2_gpu
    real(rkind), dimension(:,:), target :: dydetanc2_gpu
    integer, dimension(3) :: is_periodic

    if (is_periodic(1) == 1) then
      if (ileftx/=mpi_proc_null) then
        call bcswap_c2_step_3_kernel1_wrapper(c_null_ptr,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,&
        &irightx,irighty,irightz,c_loc(w_gpu),c_loc(wbuf1r_gpu),c_loc(wbuf2r_gpu),c_loc(wbuf3r_gpu),&
        &c_loc(wbuf4r_gpu),c_loc(wbuf5r_gpu),c_loc(wbuf6r_gpu),c_loc(dxdcsinc2_gpu),c_loc(dydcsinc2_gpu),&
        &c_loc(dxdetanc2_gpu),c_loc(dydetanc2_gpu))
        call hipCheck(hipDeviceSynchronize())
      endif
      if (irightx/=mpi_proc_null) then
        call bcswap_c2_step_3_kernel2_wrapper(c_null_ptr,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,&
        &irightx,irighty,irightz,c_loc(w_gpu),c_loc(wbuf1r_gpu),c_loc(wbuf2r_gpu),c_loc(wbuf3r_gpu),&
        &c_loc(wbuf4r_gpu),c_loc(wbuf5r_gpu),c_loc(wbuf6r_gpu),c_loc(dxdcsinc2_gpu),c_loc(dydcsinc2_gpu),&
        &c_loc(dxdetanc2_gpu),c_loc(dydetanc2_gpu))
        call hipCheck(hipDeviceSynchronize())
      endif
    else
      if (ileftx/=mpi_proc_null) then
        call bcswap_c2_step_3_kernel3_wrapper(c_null_ptr,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,&
        &irightx,irighty,irightz,c_loc(w_gpu),c_loc(wbuf1r_gpu),c_loc(wbuf2r_gpu),c_loc(wbuf3r_gpu),&
        &c_loc(wbuf4r_gpu),c_loc(wbuf5r_gpu),c_loc(wbuf6r_gpu),c_loc(dxdcsinc2_gpu),c_loc(dydcsinc2_gpu),&
        &c_loc(dxdetanc2_gpu),c_loc(dydetanc2_gpu))
        call hipCheck(hipDeviceSynchronize())
      endif
      if (irightx/=mpi_proc_null) then
        call bcswap_c2_step_3_kernel4_wrapper(c_null_ptr,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,&
        &irightx,irighty,irightz,c_loc(w_gpu),c_loc(wbuf1r_gpu),c_loc(wbuf2r_gpu),c_loc(wbuf3r_gpu),&
        &c_loc(wbuf4r_gpu),c_loc(wbuf5r_gpu),c_loc(wbuf6r_gpu),c_loc(dxdcsinc2_gpu),c_loc(dydcsinc2_gpu),&
        &c_loc(dxdetanc2_gpu),c_loc(dydetanc2_gpu))
        call hipCheck(hipDeviceSynchronize())
      endif
    endif
    if (ilefty/=mpi_proc_null) then
      call bcswap_c2_step_3_kernel5_wrapper(c_null_ptr,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,&
      &irightx,irighty,irightz,c_loc(w_gpu),c_loc(wbuf1r_gpu),c_loc(wbuf2r_gpu),c_loc(wbuf3r_gpu),&
      &c_loc(wbuf4r_gpu),c_loc(wbuf5r_gpu),c_loc(wbuf6r_gpu),c_loc(dxdcsinc2_gpu),c_loc(dydcsinc2_gpu),&
      &c_loc(dxdetanc2_gpu),c_loc(dydetanc2_gpu))
      call hipCheck(hipDeviceSynchronize())
    endif
    if (irighty/=mpi_proc_null) then
      call bcswap_c2_step_3_kernel6_wrapper(c_null_ptr,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,&
      &irightx,irighty,irightz,c_loc(w_gpu),c_loc(wbuf1r_gpu),c_loc(wbuf2r_gpu),c_loc(wbuf3r_gpu),&
      &c_loc(wbuf4r_gpu),c_loc(wbuf5r_gpu),c_loc(wbuf6r_gpu),c_loc(dxdcsinc2_gpu),c_loc(dydcsinc2_gpu),&
      &c_loc(dxdetanc2_gpu),c_loc(dydetanc2_gpu))
      call hipCheck(hipDeviceSynchronize())
    endif
    if(ndim == 3) then
      if (ileftz/=mpi_proc_null) then
        call bcswap_c2_step_3_kernel7_wrapper(c_null_ptr,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,&
        &irightx,irighty,irightz,c_loc(w_gpu),c_loc(wbuf1r_gpu),c_loc(wbuf2r_gpu),c_loc(wbuf3r_gpu),&
        &c_loc(wbuf4r_gpu),c_loc(wbuf5r_gpu),c_loc(wbuf6r_gpu),c_loc(dxdcsinc2_gpu),c_loc(dydcsinc2_gpu),&
        &c_loc(dxdetanc2_gpu),c_loc(dydetanc2_gpu))
        call hipCheck(hipDeviceSynchronize())
      endif
      if (irightz/=mpi_proc_null) then
        call bcswap_c2_step_3_kernel8_wrapper(c_null_ptr,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,&
        &irightx,irighty,irightz,c_loc(w_gpu),c_loc(wbuf1r_gpu),c_loc(wbuf2r_gpu),c_loc(wbuf3r_gpu),&
        &c_loc(wbuf4r_gpu),c_loc(wbuf5r_gpu),c_loc(wbuf6r_gpu),c_loc(dxdcsinc2_gpu),c_loc(dydcsinc2_gpu),&
        &c_loc(dxdetanc2_gpu),c_loc(dydetanc2_gpu))
        call hipCheck(hipDeviceSynchronize())
      endif
    endif

  endsubroutine bcswap_c2_step_3_kernel


  subroutine bcswap_c2xybc_step_3_kernel(nx,ny,nz,nv,ng,ndim,w_gpu,wbuf1xybcr_gpu,wbuf2xybcr_gpu,&
  &wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,ileftx,ilefty,ileftz,irightx,irighty,irightz,&
  &dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu,is_periodic)
    integer :: nx,ny,nz
    integer :: ng,nv,ndim
    integer :: ileftx,ilefty,ileftz
    integer :: irightx,irighty,irightz
    real(rkind), dimension(:,:,:,:), target :: w_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf1xybcr_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf2xybcr_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf3r_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf4r_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf5r_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf6r_gpu
    real(rkind), dimension(:,:), target :: dxdcsinc2_gpu
    real(rkind), dimension(:,:), target :: dydcsinc2_gpu
    real(rkind), dimension(:,:), target :: dxdetanc2_gpu
    real(rkind), dimension(:,:), target :: dydetanc2_gpu
    integer, dimension(3) :: is_periodic

    if (is_periodic(1) == 1) then
      if (ileftx/=mpi_proc_null) then
        call bcswap_c2xybc_step_3_kernel1_wrapper(c_null_ptr,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,&
        &ileftz,irightx,irighty,irightz,c_loc(w_gpu),c_loc(wbuf1xybcr_gpu),c_loc(wbuf2xybcr_gpu),&
        &c_loc(wbuf3r_gpu),c_loc(wbuf4r_gpu),c_loc(wbuf5r_gpu),c_loc(wbuf6r_gpu),c_loc(dxdcsinc2_gpu),&
        &c_loc(dydcsinc2_gpu),c_loc(dxdetanc2_gpu),c_loc(dydetanc2_gpu))
        call hipCheck(hipDeviceSynchronize())
      endif
      if (irightx/=mpi_proc_null) then
        call bcswap_c2xybc_step_3_kernel2_wrapper(c_null_ptr,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,&
        &ileftz,irightx,irighty,irightz,c_loc(w_gpu),c_loc(wbuf1xybcr_gpu),c_loc(wbuf2xybcr_gpu),&
        &c_loc(wbuf3r_gpu),c_loc(wbuf4r_gpu),c_loc(wbuf5r_gpu),c_loc(wbuf6r_gpu),c_loc(dxdcsinc2_gpu),&
        &c_loc(dydcsinc2_gpu),c_loc(dxdetanc2_gpu),c_loc(dydetanc2_gpu))
        call hipCheck(hipDeviceSynchronize())
      endif
    else
      if (ileftx/=mpi_proc_null) then
        call bcswap_c2xybc_step_3_kernel3_wrapper(c_null_ptr,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,&
        &ileftz,irightx,irighty,irightz,c_loc(w_gpu),c_loc(wbuf1xybcr_gpu),c_loc(wbuf2xybcr_gpu),&
        &c_loc(wbuf3r_gpu),c_loc(wbuf4r_gpu),c_loc(wbuf5r_gpu),c_loc(wbuf6r_gpu),c_loc(dxdcsinc2_gpu),&
        &c_loc(dydcsinc2_gpu),c_loc(dxdetanc2_gpu),c_loc(dydetanc2_gpu))
        call hipCheck(hipDeviceSynchronize())
      endif
      if (irightx/=mpi_proc_null) then
        call bcswap_c2xybc_step_3_kernel4_wrapper(c_null_ptr,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,&
        &ileftz,irightx,irighty,irightz,c_loc(w_gpu),c_loc(wbuf1xybcr_gpu),c_loc(wbuf2xybcr_gpu),&
        &c_loc(wbuf3r_gpu),c_loc(wbuf4r_gpu),c_loc(wbuf5r_gpu),c_loc(wbuf6r_gpu),c_loc(dxdcsinc2_gpu),&
        &c_loc(dydcsinc2_gpu),c_loc(dxdetanc2_gpu),c_loc(dydetanc2_gpu))
        call hipCheck(hipDeviceSynchronize())
      endif
    endif
    if (ilefty/=mpi_proc_null) then
      call bcswap_c2xybc_step_3_kernel5_wrapper(c_null_ptr,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,&
      &irightx,irighty,irightz,c_loc(w_gpu),c_loc(wbuf1xybcr_gpu),c_loc(wbuf2xybcr_gpu),c_loc(wbuf3r_gpu),&
      &c_loc(wbuf4r_gpu),c_loc(wbuf5r_gpu),c_loc(wbuf6r_gpu),c_loc(dxdcsinc2_gpu),c_loc(dydcsinc2_gpu),&
      &c_loc(dxdetanc2_gpu),c_loc(dydetanc2_gpu))
      call hipCheck(hipDeviceSynchronize())
    endif
    if (irighty/=mpi_proc_null) then
      call bcswap_c2xybc_step_3_kernel6_wrapper(c_null_ptr,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,&
      &irightx,irighty,irightz,c_loc(w_gpu),c_loc(wbuf1xybcr_gpu),c_loc(wbuf2xybcr_gpu),c_loc(wbuf3r_gpu),&
      &c_loc(wbuf4r_gpu),c_loc(wbuf5r_gpu),c_loc(wbuf6r_gpu),c_loc(dxdcsinc2_gpu),c_loc(dydcsinc2_gpu),&
      &c_loc(dxdetanc2_gpu),c_loc(dydetanc2_gpu))
      call hipCheck(hipDeviceSynchronize())
    endif
    if(ndim == 3) then
      if (ileftz/=mpi_proc_null) then
        call bcswap_c2xybc_step_3_kernel7_wrapper(c_null_ptr,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,&
        &ileftz,irightx,irighty,irightz,c_loc(w_gpu),c_loc(wbuf1xybcr_gpu),c_loc(wbuf2xybcr_gpu),&
        &c_loc(wbuf3r_gpu),c_loc(wbuf4r_gpu),c_loc(wbuf5r_gpu),c_loc(wbuf6r_gpu),c_loc(dxdcsinc2_gpu),&
        &c_loc(dydcsinc2_gpu),c_loc(dxdetanc2_gpu),c_loc(dydetanc2_gpu))
        call hipCheck(hipDeviceSynchronize())
      endif
      if (irightz/=mpi_proc_null) then
        call bcswap_c2xybc_step_3_kernel8_wrapper(c_null_ptr,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,&
        &ileftz,irightx,irighty,irightz,c_loc(w_gpu),c_loc(wbuf1xybcr_gpu),c_loc(wbuf2xybcr_gpu),&
        &c_loc(wbuf3r_gpu),c_loc(wbuf4r_gpu),c_loc(wbuf5r_gpu),c_loc(wbuf6r_gpu),c_loc(dxdcsinc2_gpu),&
        &c_loc(dydcsinc2_gpu),c_loc(dxdetanc2_gpu),c_loc(dydetanc2_gpu))
        call hipCheck(hipDeviceSynchronize())
      endif
    endif

  endsubroutine bcswap_c2xybc_step_3_kernel


  subroutine bcswap_corner_step_3_kernel(nx,ny,nz,nv,ng,w_gpu,wbuf1r_c_gpu,wbuf2r_c_gpu,&
  &wbuf3r_c_gpu,wbuf4r_c_gpu,ileftbottom,ilefttop,irightbottom,irighttop)
    integer :: nx,ny,nz
    integer :: ng,nv,ileftbottom
    integer :: ilefttop,irightbottom,irighttop
    real(rkind), dimension(:,:,:,:), target :: w_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf1r_c_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf2r_c_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf3r_c_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf4r_c_gpu

    if (ileftbottom/=mpi_proc_null) then
      call bcswap_corner_step_3_kernel1_wrapper(c_null_ptr,nx,ny,nz,ng,nv,ileftbottom,ilefttop,&
      &irightbottom,irighttop,c_loc(w_gpu),c_loc(wbuf1r_c_gpu),c_loc(wbuf2r_c_gpu),c_loc(wbuf3r_c_gpu),&
      &c_loc(wbuf4r_c_gpu))
      call hipCheck(hipDeviceSynchronize())
    endif
    if (irighttop/=mpi_proc_null) then
      call bcswap_corner_step_3_kernel2_wrapper(c_null_ptr,nx,ny,nz,ng,nv,ileftbottom,ilefttop,&
      &irightbottom,irighttop,c_loc(w_gpu),c_loc(wbuf1r_c_gpu),c_loc(wbuf2r_c_gpu),c_loc(wbuf3r_c_gpu),&
      &c_loc(wbuf4r_c_gpu))
      call hipCheck(hipDeviceSynchronize())
    endif
    if (ilefttop/=mpi_proc_null) then
      call bcswap_corner_step_3_kernel3_wrapper(c_null_ptr,nx,ny,nz,ng,nv,ileftbottom,ilefttop,&
      &irightbottom,irighttop,c_loc(w_gpu),c_loc(wbuf1r_c_gpu),c_loc(wbuf2r_c_gpu),c_loc(wbuf3r_c_gpu),&
      &c_loc(wbuf4r_c_gpu))
      call hipCheck(hipDeviceSynchronize())
    endif
    if (irightbottom/=mpi_proc_null) then
      call bcswap_corner_step_3_kernel4_wrapper(c_null_ptr,nx,ny,nz,ng,nv,ileftbottom,ilefttop,&
      &irightbottom,irighttop,c_loc(w_gpu),c_loc(wbuf1r_c_gpu),c_loc(wbuf2r_c_gpu),c_loc(wbuf3r_c_gpu),&
      &c_loc(wbuf4r_c_gpu))
      call hipCheck(hipDeviceSynchronize())
    endif

  endsubroutine bcswap_corner_step_3_kernel


  subroutine bcswap_corner_step_3b_kernel(nx,ny,nz,nv,ng,w_gpu,wbuf_lxly_r_gpu,wbuf_lxry_r_gpu,&
  &wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,wbuf_ryrz_r_gpu,&
  &wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,wbuf_rxlylz_r_gpu,&
  &wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu,lxly,lxry,rxly,rxry,lylz,lyrz,rylz,ryrz,&
  &lxlylz,lxlyrz,lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz)
    integer :: nx,ny,nz
    integer :: ng,nv,lxly
    integer :: lxry,rxly,rxry
    integer :: lylz,lyrz,rylz
    integer :: ryrz,lxlylz,lxlyrz
    integer :: lxrylz,lxryrz,rxlylz
    integer :: rxlyrz,rxrylz,rxryrz
    real(rkind), dimension(:,:,:,:), target :: w_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf_lxly_r_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf_lxry_r_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf_rxly_r_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf_rxry_r_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf_lylz_r_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf_lyrz_r_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf_rylz_r_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf_ryrz_r_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf_lxlylz_r_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf_lxlyrz_r_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf_lxrylz_r_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf_lxryrz_r_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf_rxlylz_r_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf_rxlyrz_r_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf_rxrylz_r_gpu
    real(rkind), dimension(:,:,:,:), target :: wbuf_rxryrz_r_gpu

    if (lxly/=mpi_proc_null) then
      call bcswap_corner_step_3b_kernel1_wrapper(c_null_ptr,nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,lylz,&
      &lyrz,rylz,ryrz,lxlylz,lxlyrz,lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz,c_loc(w_gpu),&
      &c_loc(wbuf_lxly_r_gpu),c_loc(wbuf_lxry_r_gpu),c_loc(wbuf_rxly_r_gpu),c_loc(wbuf_rxry_r_gpu),&
      &c_loc(wbuf_lylz_r_gpu),c_loc(wbuf_lyrz_r_gpu),c_loc(wbuf_rylz_r_gpu),c_loc(wbuf_ryrz_r_gpu),&
      &c_loc(wbuf_lxlylz_r_gpu),c_loc(wbuf_lxlyrz_r_gpu),c_loc(wbuf_lxrylz_r_gpu),c_loc(wbuf_lxryrz_r_gpu),&
      &c_loc(wbuf_rxlylz_r_gpu),c_loc(wbuf_rxlyrz_r_gpu),c_loc(wbuf_rxrylz_r_gpu),c_loc(wbuf_rxryrz_r_gpu))
      call hipCheck(hipDeviceSynchronize())
    endif
    if (rxry/=mpi_proc_null) then
      call bcswap_corner_step_3b_kernel2_wrapper(c_null_ptr,nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,lylz,&
      &lyrz,rylz,ryrz,lxlylz,lxlyrz,lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz,c_loc(w_gpu),&
      &c_loc(wbuf_lxly_r_gpu),c_loc(wbuf_lxry_r_gpu),c_loc(wbuf_rxly_r_gpu),c_loc(wbuf_rxry_r_gpu),&
      &c_loc(wbuf_lylz_r_gpu),c_loc(wbuf_lyrz_r_gpu),c_loc(wbuf_rylz_r_gpu),c_loc(wbuf_ryrz_r_gpu),&
      &c_loc(wbuf_lxlylz_r_gpu),c_loc(wbuf_lxlyrz_r_gpu),c_loc(wbuf_lxrylz_r_gpu),c_loc(wbuf_lxryrz_r_gpu),&
      &c_loc(wbuf_rxlylz_r_gpu),c_loc(wbuf_rxlyrz_r_gpu),c_loc(wbuf_rxrylz_r_gpu),c_loc(wbuf_rxryrz_r_gpu))
      call hipCheck(hipDeviceSynchronize())
    endif
    if (lxry/=mpi_proc_null) then
      call bcswap_corner_step_3b_kernel3_wrapper(c_null_ptr,nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,lylz,&
      &lyrz,rylz,ryrz,lxlylz,lxlyrz,lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz,c_loc(w_gpu),&
      &c_loc(wbuf_lxly_r_gpu),c_loc(wbuf_lxry_r_gpu),c_loc(wbuf_rxly_r_gpu),c_loc(wbuf_rxry_r_gpu),&
      &c_loc(wbuf_lylz_r_gpu),c_loc(wbuf_lyrz_r_gpu),c_loc(wbuf_rylz_r_gpu),c_loc(wbuf_ryrz_r_gpu),&
      &c_loc(wbuf_lxlylz_r_gpu),c_loc(wbuf_lxlyrz_r_gpu),c_loc(wbuf_lxrylz_r_gpu),c_loc(wbuf_lxryrz_r_gpu),&
      &c_loc(wbuf_rxlylz_r_gpu),c_loc(wbuf_rxlyrz_r_gpu),c_loc(wbuf_rxrylz_r_gpu),c_loc(wbuf_rxryrz_r_gpu))
      call hipCheck(hipDeviceSynchronize())
    endif
    if (rxly/=mpi_proc_null) then
      call bcswap_corner_step_3b_kernel4_wrapper(c_null_ptr,nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,lylz,&
      &lyrz,rylz,ryrz,lxlylz,lxlyrz,lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz,c_loc(w_gpu),&
      &c_loc(wbuf_lxly_r_gpu),c_loc(wbuf_lxry_r_gpu),c_loc(wbuf_rxly_r_gpu),c_loc(wbuf_rxry_r_gpu),&
      &c_loc(wbuf_lylz_r_gpu),c_loc(wbuf_lyrz_r_gpu),c_loc(wbuf_rylz_r_gpu),c_loc(wbuf_ryrz_r_gpu),&
      &c_loc(wbuf_lxlylz_r_gpu),c_loc(wbuf_lxlyrz_r_gpu),c_loc(wbuf_lxrylz_r_gpu),c_loc(wbuf_lxryrz_r_gpu),&
      &c_loc(wbuf_rxlylz_r_gpu),c_loc(wbuf_rxlyrz_r_gpu),c_loc(wbuf_rxrylz_r_gpu),c_loc(wbuf_rxryrz_r_gpu))
      call hipCheck(hipDeviceSynchronize())
    endif

    if (lylz/=mpi_proc_null) then
      call bcswap_corner_step_3b_kernel5_wrapper(c_null_ptr,nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,lylz,&
      &lyrz,rylz,ryrz,lxlylz,lxlyrz,lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz,c_loc(w_gpu),&
      &c_loc(wbuf_lxly_r_gpu),c_loc(wbuf_lxry_r_gpu),c_loc(wbuf_rxly_r_gpu),c_loc(wbuf_rxry_r_gpu),&
      &c_loc(wbuf_lylz_r_gpu),c_loc(wbuf_lyrz_r_gpu),c_loc(wbuf_rylz_r_gpu),c_loc(wbuf_ryrz_r_gpu),&
      &c_loc(wbuf_lxlylz_r_gpu),c_loc(wbuf_lxlyrz_r_gpu),c_loc(wbuf_lxrylz_r_gpu),c_loc(wbuf_lxryrz_r_gpu),&
      &c_loc(wbuf_rxlylz_r_gpu),c_loc(wbuf_rxlyrz_r_gpu),c_loc(wbuf_rxrylz_r_gpu),c_loc(wbuf_rxryrz_r_gpu))
      call hipCheck(hipDeviceSynchronize())
    endif
    if (ryrz/=mpi_proc_null) then
      call bcswap_corner_step_3b_kernel6_wrapper(c_null_ptr,nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,lylz,&
      &lyrz,rylz,ryrz,lxlylz,lxlyrz,lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz,c_loc(w_gpu),&
      &c_loc(wbuf_lxly_r_gpu),c_loc(wbuf_lxry_r_gpu),c_loc(wbuf_rxly_r_gpu),c_loc(wbuf_rxry_r_gpu),&
      &c_loc(wbuf_lylz_r_gpu),c_loc(wbuf_lyrz_r_gpu),c_loc(wbuf_rylz_r_gpu),c_loc(wbuf_ryrz_r_gpu),&
      &c_loc(wbuf_lxlylz_r_gpu),c_loc(wbuf_lxlyrz_r_gpu),c_loc(wbuf_lxrylz_r_gpu),c_loc(wbuf_lxryrz_r_gpu),&
      &c_loc(wbuf_rxlylz_r_gpu),c_loc(wbuf_rxlyrz_r_gpu),c_loc(wbuf_rxrylz_r_gpu),c_loc(wbuf_rxryrz_r_gpu))
      call hipCheck(hipDeviceSynchronize())
    endif
    if (lyrz/=mpi_proc_null) then
      call bcswap_corner_step_3b_kernel7_wrapper(c_null_ptr,nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,lylz,&
      &lyrz,rylz,ryrz,lxlylz,lxlyrz,lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz,c_loc(w_gpu),&
      &c_loc(wbuf_lxly_r_gpu),c_loc(wbuf_lxry_r_gpu),c_loc(wbuf_rxly_r_gpu),c_loc(wbuf_rxry_r_gpu),&
      &c_loc(wbuf_lylz_r_gpu),c_loc(wbuf_lyrz_r_gpu),c_loc(wbuf_rylz_r_gpu),c_loc(wbuf_ryrz_r_gpu),&
      &c_loc(wbuf_lxlylz_r_gpu),c_loc(wbuf_lxlyrz_r_gpu),c_loc(wbuf_lxrylz_r_gpu),c_loc(wbuf_lxryrz_r_gpu),&
      &c_loc(wbuf_rxlylz_r_gpu),c_loc(wbuf_rxlyrz_r_gpu),c_loc(wbuf_rxrylz_r_gpu),c_loc(wbuf_rxryrz_r_gpu))
      call hipCheck(hipDeviceSynchronize())
    endif
    if (rylz/=mpi_proc_null) then
      call bcswap_corner_step_3b_kernel8_wrapper(c_null_ptr,nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,lylz,&
      &lyrz,rylz,ryrz,lxlylz,lxlyrz,lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz,c_loc(w_gpu),&
      &c_loc(wbuf_lxly_r_gpu),c_loc(wbuf_lxry_r_gpu),c_loc(wbuf_rxly_r_gpu),c_loc(wbuf_rxry_r_gpu),&
      &c_loc(wbuf_lylz_r_gpu),c_loc(wbuf_lyrz_r_gpu),c_loc(wbuf_rylz_r_gpu),c_loc(wbuf_ryrz_r_gpu),&
      &c_loc(wbuf_lxlylz_r_gpu),c_loc(wbuf_lxlyrz_r_gpu),c_loc(wbuf_lxrylz_r_gpu),c_loc(wbuf_lxryrz_r_gpu),&
      &c_loc(wbuf_rxlylz_r_gpu),c_loc(wbuf_rxlyrz_r_gpu),c_loc(wbuf_rxrylz_r_gpu),c_loc(wbuf_rxryrz_r_gpu))
      call hipCheck(hipDeviceSynchronize())
    endif

    if (lxlylz/=mpi_proc_null) then
      call bcswap_corner_step_3b_kernel9_wrapper(c_null_ptr,nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,lylz,&
      &lyrz,rylz,ryrz,lxlylz,lxlyrz,lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz,c_loc(w_gpu),&
      &c_loc(wbuf_lxly_r_gpu),c_loc(wbuf_lxry_r_gpu),c_loc(wbuf_rxly_r_gpu),c_loc(wbuf_rxry_r_gpu),&
      &c_loc(wbuf_lylz_r_gpu),c_loc(wbuf_lyrz_r_gpu),c_loc(wbuf_rylz_r_gpu),c_loc(wbuf_ryrz_r_gpu),&
      &c_loc(wbuf_lxlylz_r_gpu),c_loc(wbuf_lxlyrz_r_gpu),c_loc(wbuf_lxrylz_r_gpu),c_loc(wbuf_lxryrz_r_gpu),&
      &c_loc(wbuf_rxlylz_r_gpu),c_loc(wbuf_rxlyrz_r_gpu),c_loc(wbuf_rxrylz_r_gpu),c_loc(wbuf_rxryrz_r_gpu))
      call hipCheck(hipDeviceSynchronize())
    endif
    if (lxlyrz/=mpi_proc_null) then
      call bcswap_corner_step_3b_kernel10_wrapper(c_null_ptr,nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,&
      &lylz,lyrz,rylz,ryrz,lxlylz,lxlyrz,lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz,c_loc(w_gpu),&
      &c_loc(wbuf_lxly_r_gpu),c_loc(wbuf_lxry_r_gpu),c_loc(wbuf_rxly_r_gpu),c_loc(wbuf_rxry_r_gpu),&
      &c_loc(wbuf_lylz_r_gpu),c_loc(wbuf_lyrz_r_gpu),c_loc(wbuf_rylz_r_gpu),c_loc(wbuf_ryrz_r_gpu),&
      &c_loc(wbuf_lxlylz_r_gpu),c_loc(wbuf_lxlyrz_r_gpu),c_loc(wbuf_lxrylz_r_gpu),c_loc(wbuf_lxryrz_r_gpu),&
      &c_loc(wbuf_rxlylz_r_gpu),c_loc(wbuf_rxlyrz_r_gpu),c_loc(wbuf_rxrylz_r_gpu),c_loc(wbuf_rxryrz_r_gpu))
      call hipCheck(hipDeviceSynchronize())
    endif
    if (lxrylz/=mpi_proc_null) then
      call bcswap_corner_step_3b_kernel11_wrapper(c_null_ptr,nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,&
      &lylz,lyrz,rylz,ryrz,lxlylz,lxlyrz,lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz,c_loc(w_gpu),&
      &c_loc(wbuf_lxly_r_gpu),c_loc(wbuf_lxry_r_gpu),c_loc(wbuf_rxly_r_gpu),c_loc(wbuf_rxry_r_gpu),&
      &c_loc(wbuf_lylz_r_gpu),c_loc(wbuf_lyrz_r_gpu),c_loc(wbuf_rylz_r_gpu),c_loc(wbuf_ryrz_r_gpu),&
      &c_loc(wbuf_lxlylz_r_gpu),c_loc(wbuf_lxlyrz_r_gpu),c_loc(wbuf_lxrylz_r_gpu),c_loc(wbuf_lxryrz_r_gpu),&
      &c_loc(wbuf_rxlylz_r_gpu),c_loc(wbuf_rxlyrz_r_gpu),c_loc(wbuf_rxrylz_r_gpu),c_loc(wbuf_rxryrz_r_gpu))
      call hipCheck(hipDeviceSynchronize())
    endif
    if (lxryrz/=mpi_proc_null) then
      call bcswap_corner_step_3b_kernel12_wrapper(c_null_ptr,nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,&
      &lylz,lyrz,rylz,ryrz,lxlylz,lxlyrz,lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz,c_loc(w_gpu),&
      &c_loc(wbuf_lxly_r_gpu),c_loc(wbuf_lxry_r_gpu),c_loc(wbuf_rxly_r_gpu),c_loc(wbuf_rxry_r_gpu),&
      &c_loc(wbuf_lylz_r_gpu),c_loc(wbuf_lyrz_r_gpu),c_loc(wbuf_rylz_r_gpu),c_loc(wbuf_ryrz_r_gpu),&
      &c_loc(wbuf_lxlylz_r_gpu),c_loc(wbuf_lxlyrz_r_gpu),c_loc(wbuf_lxrylz_r_gpu),c_loc(wbuf_lxryrz_r_gpu),&
      &c_loc(wbuf_rxlylz_r_gpu),c_loc(wbuf_rxlyrz_r_gpu),c_loc(wbuf_rxrylz_r_gpu),c_loc(wbuf_rxryrz_r_gpu))
      call hipCheck(hipDeviceSynchronize())
    endif
    if (rxlylz/=mpi_proc_null) then
      call bcswap_corner_step_3b_kernel13_wrapper(c_null_ptr,nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,&
      &lylz,lyrz,rylz,ryrz,lxlylz,lxlyrz,lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz,c_loc(w_gpu),&
      &c_loc(wbuf_lxly_r_gpu),c_loc(wbuf_lxry_r_gpu),c_loc(wbuf_rxly_r_gpu),c_loc(wbuf_rxry_r_gpu),&
      &c_loc(wbuf_lylz_r_gpu),c_loc(wbuf_lyrz_r_gpu),c_loc(wbuf_rylz_r_gpu),c_loc(wbuf_ryrz_r_gpu),&
      &c_loc(wbuf_lxlylz_r_gpu),c_loc(wbuf_lxlyrz_r_gpu),c_loc(wbuf_lxrylz_r_gpu),c_loc(wbuf_lxryrz_r_gpu),&
      &c_loc(wbuf_rxlylz_r_gpu),c_loc(wbuf_rxlyrz_r_gpu),c_loc(wbuf_rxrylz_r_gpu),c_loc(wbuf_rxryrz_r_gpu))
      call hipCheck(hipDeviceSynchronize())
    endif
    if (rxlyrz/=mpi_proc_null) then
      call bcswap_corner_step_3b_kernel14_wrapper(c_null_ptr,nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,&
      &lylz,lyrz,rylz,ryrz,lxlylz,lxlyrz,lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz,c_loc(w_gpu),&
      &c_loc(wbuf_lxly_r_gpu),c_loc(wbuf_lxry_r_gpu),c_loc(wbuf_rxly_r_gpu),c_loc(wbuf_rxry_r_gpu),&
      &c_loc(wbuf_lylz_r_gpu),c_loc(wbuf_lyrz_r_gpu),c_loc(wbuf_rylz_r_gpu),c_loc(wbuf_ryrz_r_gpu),&
      &c_loc(wbuf_lxlylz_r_gpu),c_loc(wbuf_lxlyrz_r_gpu),c_loc(wbuf_lxrylz_r_gpu),c_loc(wbuf_lxryrz_r_gpu),&
      &c_loc(wbuf_rxlylz_r_gpu),c_loc(wbuf_rxlyrz_r_gpu),c_loc(wbuf_rxrylz_r_gpu),c_loc(wbuf_rxryrz_r_gpu))
      call hipCheck(hipDeviceSynchronize())
    endif
    if (rxrylz/=mpi_proc_null) then
      call bcswap_corner_step_3b_kernel15_wrapper(c_null_ptr,nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,&
      &lylz,lyrz,rylz,ryrz,lxlylz,lxlyrz,lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz,c_loc(w_gpu),&
      &c_loc(wbuf_lxly_r_gpu),c_loc(wbuf_lxry_r_gpu),c_loc(wbuf_rxly_r_gpu),c_loc(wbuf_rxry_r_gpu),&
      &c_loc(wbuf_lylz_r_gpu),c_loc(wbuf_lyrz_r_gpu),c_loc(wbuf_rylz_r_gpu),c_loc(wbuf_ryrz_r_gpu),&
      &c_loc(wbuf_lxlylz_r_gpu),c_loc(wbuf_lxlyrz_r_gpu),c_loc(wbuf_lxrylz_r_gpu),c_loc(wbuf_lxryrz_r_gpu),&
      &c_loc(wbuf_rxlylz_r_gpu),c_loc(wbuf_rxlyrz_r_gpu),c_loc(wbuf_rxrylz_r_gpu),c_loc(wbuf_rxryrz_r_gpu))
      call hipCheck(hipDeviceSynchronize())
    endif
    if (rxryrz/=mpi_proc_null) then
      call bcswap_corner_step_3b_kernel16_wrapper(c_null_ptr,nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,&
      &lylz,lyrz,rylz,ryrz,lxlylz,lxlyrz,lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz,c_loc(w_gpu),&
      &c_loc(wbuf_lxly_r_gpu),c_loc(wbuf_lxry_r_gpu),c_loc(wbuf_rxly_r_gpu),c_loc(wbuf_rxry_r_gpu),&
      &c_loc(wbuf_lylz_r_gpu),c_loc(wbuf_lyrz_r_gpu),c_loc(wbuf_rylz_r_gpu),c_loc(wbuf_ryrz_r_gpu),&
      &c_loc(wbuf_lxlylz_r_gpu),c_loc(wbuf_lxlyrz_r_gpu),c_loc(wbuf_lxrylz_r_gpu),c_loc(wbuf_lxryrz_r_gpu),&
      &c_loc(wbuf_rxlylz_r_gpu),c_loc(wbuf_rxlyrz_r_gpu),c_loc(wbuf_rxrylz_r_gpu),c_loc(wbuf_rxryrz_r_gpu))
      call hipCheck(hipDeviceSynchronize())
    endif

  endsubroutine bcswap_corner_step_3b_kernel


  subroutine bcswap_corner(self, steps)
    class(base_amd_object), intent(inout) :: self
    logical, dimension(3), optional :: steps
    logical, dimension(3) :: steps_
    integer :: iercuda, indc
    integer, dimension(mpi_status_size) :: istatus

    steps_ = .true. ; if(present(steps)) steps_ = steps

    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv,&
    & w_gpu => self%w_gpu, wbuf1s_c_gpu => self%wbuf1s_c_gpu, wbuf2s_c_gpu => self%wbuf2s_c_gpu,&
    & wbuf3s_c_gpu => self%wbuf3s_c_gpu, wbuf4s_c_gpu => self%wbuf4s_c_gpu, wbuf1r_c_gpu => self%wbuf1r_c&
    &_gpu, wbuf2r_c_gpu => self%wbuf2r_c_gpu, wbuf3r_c_gpu => self%wbuf3r_c_gpu, wbuf4r_c_gpu => self%wbu&
    &f4r_c_gpu, ileftbottom => self%field%ileftbottom, irightbottom => self%field%irightbottom,&
    & ilefttop => self%field%ilefttop, irighttop => self%field%irighttop, mp_cart => self%field%mp_cart,&
    & iermpi => self%mpi_err, nrank => self%field%myrank)

      if(steps_(1)) then
        call bcswap_corner_step_1_kernel(nx,ny,nz,nv,ng,w_gpu,wbuf1s_c_gpu,wbuf2s_c_gpu,wbuf3s_c_gpu,wbuf4s_c_gpu)
      endif
      if(steps_(2)) then
        indc = nv*ny*ng*ng
        if (ileftbottom == nrank) then
          call hipCheck(hipMemcpy(wbuf2r_c_gpu , wbuf1s_c_gpu,hipMemcpyDeviceToDevice))
          call hipCheck(hipMemcpy(wbuf1r_c_gpu , wbuf2s_c_gpu,hipMemcpyDeviceToDevice))
        else
          call mpi_sendrecv(wbuf1s_c_gpu,indc,mpi_prec,ileftbottom ,1, wbuf2r_c_gpu,indc,mpi_prec,&
          &irighttop ,1,mp_cart,istatus,iermpi)
          call mpi_sendrecv(wbuf2s_c_gpu,indc,mpi_prec,irighttop ,2, wbuf1r_c_gpu,indc,mpi_prec,&
          &ileftbottom ,2,mp_cart,istatus,iermpi)
        endif
        if (ilefttop == nrank) then
          call hipCheck(hipMemcpy(wbuf4r_c_gpu , wbuf3s_c_gpu,hipMemcpyDeviceToDevice))
          call hipCheck(hipMemcpy(wbuf3r_c_gpu , wbuf4s_c_gpu,hipMemcpyDeviceToDevice))
        else
          call mpi_sendrecv(wbuf3s_c_gpu,indc,mpi_prec,ilefttop ,3, wbuf4r_c_gpu,indc,mpi_prec,&
          &irightbottom,3,mp_cart,istatus,iermpi)
          call mpi_sendrecv(wbuf4s_c_gpu,indc,mpi_prec,irightbottom,4, wbuf3r_c_gpu,indc,mpi_prec,&
          &ilefttop ,4,mp_cart,istatus,iermpi)
        endif
      endif
      if(steps_(3)) then
        call bcswap_corner_step_3_kernel(nx,ny,nz,nv,ng,w_gpu,wbuf1r_c_gpu,wbuf2r_c_gpu,&
        &wbuf3r_c_gpu,wbuf4r_c_gpu,ileftbottom,ilefttop,irightbottom,irighttop)
      endif
    endassociate
  endsubroutine bcswap_corner

  subroutine bcswap_corner_var(self, w_swap, steps)
    class(base_amd_object), intent(inout) :: self
    real(rkind), dimension(:,:,:,:), target :: w_swap
    logical, dimension(3), optional :: steps
    logical, dimension(3) :: steps_
    integer :: iercuda, indc
    integer, dimension(mpi_status_size) :: istatus

    steps_ = .true. ; if(present(steps)) steps_ = steps

    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv,&
    & wbuf1s_c_gpu => self%wbuf1s_c_gpu, wbuf2s_c_gpu => self%wbuf2s_c_gpu, wbuf3s_c_gpu => self%wbuf3s_c&
    &_gpu, wbuf4s_c_gpu => self%wbuf4s_c_gpu, wbuf1r_c_gpu => self%wbuf1r_c_gpu, wbuf2r_c_gpu => self%wbu&
    &f2r_c_gpu, wbuf3r_c_gpu => self%wbuf3r_c_gpu, wbuf4r_c_gpu => self%wbuf4r_c_gpu, ileftbottom => self&
    &%field%ileftbottom, irightbottom => self%field%irightbottom, ilefttop => self%field%ilefttop,&
    & irighttop => self%field%irighttop, mp_cart => self%field%mp_cart, iermpi => self%mpi_err,&
    & nrank => self%field%myrank)

      if(steps_(1)) then
        call bcswap_corner_step_1_kernel(nx,ny,nz,1,ng,w_swap,wbuf1s_c_gpu,wbuf2s_c_gpu,wbuf3s_c_gpu,wbuf4s_c_gpu)
      endif
      if(steps_(2)) then
        indc = ny*ng*ng
        if (ileftbottom == nrank) then
          call hipCheck(hipMemcpy(wbuf2r_c_gpu , wbuf1s_c_gpu,hipMemcpyDeviceToDevice))
          call hipCheck(hipMemcpy(wbuf1r_c_gpu , wbuf2s_c_gpu,hipMemcpyDeviceToDevice))
        else
          call mpi_sendrecv(wbuf1s_c_gpu,indc,mpi_prec,ileftbottom ,1, wbuf2r_c_gpu,indc,mpi_prec,&
          &irighttop ,1,mp_cart,istatus,iermpi)
          call mpi_sendrecv(wbuf2s_c_gpu,indc,mpi_prec,irighttop ,2, wbuf1r_c_gpu,indc,mpi_prec,&
          &ileftbottom ,2,mp_cart,istatus,iermpi)
        endif
        if (ilefttop == nrank) then
          call hipCheck(hipMemcpy(wbuf4r_c_gpu , wbuf3s_c_gpu,hipMemcpyDeviceToDevice))
          call hipCheck(hipMemcpy(wbuf3r_c_gpu , wbuf4s_c_gpu,hipMemcpyDeviceToDevice))
        else
          call mpi_sendrecv(wbuf3s_c_gpu,indc,mpi_prec,ilefttop ,3, wbuf4r_c_gpu,indc,mpi_prec,&
          &irightbottom,3,mp_cart,istatus,iermpi)
          call mpi_sendrecv(wbuf4s_c_gpu,indc,mpi_prec,irightbottom,4, wbuf3r_c_gpu,indc,mpi_prec,&
          &ilefttop ,4,mp_cart,istatus,iermpi)
        endif
      endif
      if(steps_(3)) then
        call bcswap_corner_step_3_kernel(nx,ny,nz,1,ng,w_swap,wbuf1r_c_gpu,wbuf2r_c_gpu,&
        &wbuf3r_c_gpu,wbuf4r_c_gpu,ileftbottom,ilefttop,irightbottom,irighttop)
      endif
    endassociate
  endsubroutine bcswap_corner_var

  subroutine bcswap(self, steps, swap_xy_corner_bc)
    class(base_amd_object), intent(inout) :: self
    logical, dimension(3), optional :: steps
    integer, optional :: swap_xy_corner_bc
    integer :: swap_xy_corner_bc_
    logical, dimension(3) :: steps_
    integer :: iercuda, indx, indy, indz
    integer, dimension(mpi_status_size) :: istatus
    integer, dimension(3) :: is_periodic

    steps_ = .true. ; if(present(steps)) steps_ = steps
    swap_xy_corner_bc_ = 1 ; if(present(swap_xy_corner_bc)) swap_xy_corner_bc_ = swap_xy_corner_bc

    is_periodic(:) = 0
    if(self%field%grid%is_xyz_periodic(1)) is_periodic(1) = 1
    if(self%field%grid%is_xyz_periodic(2)) is_periodic(2) = 1
    if(self%field%grid%is_xyz_periodic(3)) is_periodic(3) = 1

    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv,&
    & w_gpu => self%w_gpu, wbuf1s_gpu => self%wbuf1s_gpu, wbuf2s_gpu => self%wbuf2s_gpu,&
    & wbuf1xybcs_gpu => self%wbuf1xybcs_gpu, wbuf2xybcs_gpu => self%wbuf2xybcs_gpu, wbuf3s_gpu => self%wb&
    &uf3s_gpu, wbuf4s_gpu => self%wbuf4s_gpu, wbuf5s_gpu => self%wbuf5s_gpu, wbuf6s_gpu => self%wbuf6s_gp&
    &u, wbuf1r_gpu => self%wbuf1r_gpu, wbuf2r_gpu => self%wbuf2r_gpu, wbuf1xybcr_gpu => self%wbuf1xybcr_g&
    &pu, wbuf2xybcr_gpu => self%wbuf2xybcr_gpu, wbuf3r_gpu => self%wbuf3r_gpu, wbuf4r_gpu => self%wbuf4r_&
    &gpu, wbuf5r_gpu => self%wbuf5r_gpu, wbuf6r_gpu => self%wbuf6r_gpu, ileftx => self%field%ileftx,&
    &irightx => self%field%irightx, nrank_x => self%field%nrank_x, ilefty => self%field%ilefty,&
    &irighty => self%field%irighty, nrank_y => self%field%nrank_y, ileftz => self%field%ileftz,&
    &irightz => self%field%irightz, nrank_z => self%field%nrank_z, dcsidxnc2_gpu => self%dcsidxnc2_gpu,&
    & dcsidync2_gpu => self%dcsidync2_gpu, detadxnc2_gpu => self%detadxnc2_gpu, detadync2_gpu => self%det&
    &adync2_gpu, dxdcsinc2_gpu => self%dxdcsinc2_gpu, dydcsinc2_gpu => self%dydcsinc2_gpu,&
    & dxdetanc2_gpu => self%dxdetanc2_gpu, dydetanc2_gpu => self%dydetanc2_gpu, mp_cartx => self%field%mp&
    &_cartx, mp_carty => self%field%mp_carty, mp_cartz => self%field%mp_cartz, ncoords => self%field%ncoo&
    &rds, nblocks => self%field%nblocks, ndim => self%field%grid%ndim)

      if(steps_(1)) then
        if(self%field%grid%grid_dim == 1) then
          call bcswap_step_1_kernel(nx,ny,nz,nv,ng,ndim,w_gpu,wbuf1s_gpu,wbuf2s_gpu,wbuf3s_gpu,wbuf4s_gpu,wbuf5s_gpu,wbuf6s_gpu)
        elseif(self%field%grid%grid_dim == 2) then
          if(swap_xy_corner_bc_ == 0) then
            call bcswap_c2_step_1_kernel(nx,ny,nz,nv,ng,ndim,w_gpu,wbuf1s_gpu,wbuf2s_gpu,wbuf3s_gpu,&
            &wbuf4s_gpu,wbuf5s_gpu,wbuf6s_gpu,dcsidxnc2_gpu,dcsidync2_gpu,detadxnc2_gpu,detadync2_gpu,&
            &is_periodic)
          else
            call bcswap_c2xybc_step_1_kernel(nx,ny,nz,nv,ng,ndim,w_gpu,wbuf1xybcs_gpu,&
            &wbuf2xybcs_gpu,wbuf3s_gpu,wbuf4s_gpu,wbuf5s_gpu,wbuf6s_gpu,dcsidxnc2_gpu,dcsidync2_gpu,&
            &detadxnc2_gpu,detadync2_gpu,is_periodic)
          endif
        endif
      endif
      if(steps_(2)) then
        indy = nv*nx*ng*nz
        indz = nv*nx*ny*ng
        if(self%field%grid%grid_dim == 1 .or. (self%field%grid%grid_dim == 2 .and. swap_xy_corner_bc_ == 0)) then
          indx = nv*ng*ny*nz
          if(ileftx == nrank_x) then
            call hipCheck(hipMemcpy(wbuf2r_gpu , wbuf1s_gpu,hipMemcpyDeviceToDevice))
            call hipCheck(hipMemcpy(wbuf1r_gpu , wbuf2s_gpu,hipMemcpyDeviceToDevice))
          else
            call mpi_sendrecv(wbuf1s_gpu,indx,mpi_prec,ileftx ,1,wbuf2r_gpu,indx,mpi_prec,irightx,1,mp_cartx,istatus,self%mpi_err)
            call mpi_sendrecv(wbuf2s_gpu,indx,mpi_prec,irightx,2,wbuf1r_gpu,indx,mpi_prec,ileftx ,2,mp_cartx,istatus,self%mpi_err)
          endif
        else
          indx = nv*ng*(2*ng+ny)*nz
          if(ileftx == nrank_x) then
            call hipCheck(hipMemcpy(wbuf2xybcr_gpu , wbuf1xybcs_gpu,hipMemcpyDeviceToDevice))
            call hipCheck(hipMemcpy(wbuf1xybcr_gpu , wbuf2xybcs_gpu,hipMemcpyDeviceToDevice))
          else
            call mpi_sendrecv(wbuf1xybcs_gpu,indx,mpi_prec,ileftx ,1,wbuf2xybcr_gpu,indx,mpi_prec,&
            &irightx,1,mp_cartx,istatus,self%mpi_err)
            call mpi_sendrecv(wbuf2xybcs_gpu,indx,mpi_prec,irightx,2,wbuf1xybcr_gpu,indx,mpi_prec,&
            &ileftx ,2,mp_cartx,istatus,self%mpi_err)
          endif
        endif
        if(ilefty == nrank_y) then
          call hipCheck(hipMemcpy(wbuf4r_gpu , wbuf3s_gpu,hipMemcpyDeviceToDevice))
          call hipCheck(hipMemcpy(wbuf3r_gpu , wbuf4s_gpu,hipMemcpyDeviceToDevice))
        else
          call mpi_sendrecv(wbuf3s_gpu,indy,mpi_prec,ilefty ,3,wbuf4r_gpu,indy,mpi_prec,irighty,3,mp_carty,istatus,self%mpi_err)
          call mpi_sendrecv(wbuf4s_gpu,indy,mpi_prec,irighty,4,wbuf3r_gpu,indy,mpi_prec,ilefty ,4,mp_carty,istatus,self%mpi_err)
        endif
        if(ndim == 3) then
          if(ileftz == nrank_z) then
            call hipCheck(hipMemcpy(wbuf6r_gpu , wbuf5s_gpu,hipMemcpyDeviceToDevice))
            call hipCheck(hipMemcpy(wbuf5r_gpu , wbuf6s_gpu,hipMemcpyDeviceToDevice))
          else
            call mpi_sendrecv(wbuf5s_gpu,indz,mpi_prec,ileftz ,5,wbuf6r_gpu,indz,mpi_prec,irightz,5,mp_cartz,istatus,self%mpi_err)
            call mpi_sendrecv(wbuf6s_gpu,indz,mpi_prec,irightz,6,wbuf5r_gpu,indz,mpi_prec,ileftz ,6,mp_cartz,istatus,self%mpi_err)
          endif
        endif
      endif
      if(steps_(3)) then
        if(self%field%grid%grid_dim == 1) then
          call bcswap_step_3_kernel(nx,ny,nz,nv,ng,ndim,w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,&
          &wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,ileftx,ilefty,ileftz,irightx,irighty,irightz)
        elseif(self%field%grid%grid_dim == 2) then
          if(swap_xy_corner_bc_ == 0) then
            call bcswap_c2_step_3_kernel(nx,ny,nz,nv,ng,ndim,w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,&
            &wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,ileftx,ilefty,ileftz,irightx,irighty,irightz,dxdcsinc2_gpu,&
            &dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu,is_periodic)
          else
            call bcswap_c2xybc_step_3_kernel(nx,ny,nz,nv,ng,ndim,w_gpu,wbuf1xybcr_gpu,&
            &wbuf2xybcr_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,ileftx,ilefty,ileftz,irightx,irighty,&
            &irightz,dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu,is_periodic)
            call extr_corner_ymin_kernel(ncoords,nblocks,nx,ny,nz,ng,nv,w_gpu)
          endif
        endif
      endif
    endassociate
  endsubroutine bcswap

  subroutine extr_corner_ymin_kernel(ncoords,nblocks,nx,ny,nz,ng,nv,w_gpu)
    integer :: nx,ny,nz
    integer :: ng,nv
    real(rkind), dimension(:,:,:,:), target :: w_gpu
    integer, dimension(3) :: ncoords
    integer, dimension(3) :: nblocks

    if(ncoords(1) == 0) then
      call extr_corner_ymin_kernel1_wrapper(c_null_ptr,nx,ny,nz,ng,nv,c_loc(w_gpu))
      call hipCheck(hipDeviceSynchronize())
    endif
    if(ncoords(1) == nblocks(1)-1) then
      call extr_corner_ymin_kernel2_wrapper(c_null_ptr,nx,ny,nz,ng,nv,c_loc(w_gpu))
      call hipCheck(hipDeviceSynchronize())
    endif

  endsubroutine extr_corner_ymin_kernel


  subroutine bcswap_wake(self, steps)
    class(base_amd_object), intent(inout) :: self
    logical, dimension(3), optional :: steps
    logical, dimension(3) :: steps_
    integer :: iercuda, indx, indy, indz
    integer, dimension(mpi_status_size) :: istatus

    steps_ = .true. ; if(present(steps)) steps_ = steps

    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv,&
    & w_gpu => self%w_gpu, wbuf4s_gpu => self%wbuf4s_gpu, wbuf3r_gpu => self%wbuf3r_gpu,&
    & wall_tag_gpu => self%wall_tag_gpu, mp_cart => self%field%mp_cart, nrank_wake => self%field%nrank_wa&
    &ke, iermpi => self%mpi_err, nrank => self%myrank)

      if(steps_(1)) then
        call bcswap_wake_step_1_kernel(nx,ny,nz,nv,ng,w_gpu,wbuf4s_gpu)
      endif
      if(steps_(2)) then
        indy = nv*nx*ng*nz
        if (nrank_wake==nrank) then
          call hipCheck(hipMemcpy(wbuf3r_gpu , wbuf4s_gpu,hipMemcpyDeviceToDevice))
        else
          call mpi_sendrecv(wbuf4s_gpu,indy,mpi_prec,nrank_wake,4,wbuf3r_gpu,indy,mpi_prec,nrank_wake,4,mp_cart,istatus,iermpi)
        endif
      endif
      if(steps_(3)) then
        call bcswap_wake_step_3_kernel(nx,ny,nz,nv,ng,w_gpu,wbuf3r_gpu,wall_tag_gpu)
      endif
    endassociate
  endsubroutine bcswap_wake

  subroutine bcswap_wake_var(self, w_swap, steps)
    class(base_amd_object), intent(inout) :: self
    real(rkind), dimension(:,:,:,:), target :: w_swap
    logical, dimension(3), optional :: steps
    logical, dimension(3) :: steps_
    integer :: iercuda, indx, indy, indz
    integer, dimension(mpi_status_size) :: istatus

    steps_ = .true. ; if(present(steps)) steps_ = steps

    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv,&
    & w_gpu => self%w_gpu, wbuf4s_gpu => self%wbuf4s_gpu, wbuf3r_gpu => self%wbuf3r_gpu,&
    & wall_tag_gpu => self%wall_tag_gpu, mp_cart => self%field%mp_cart, nrank_wake => self%field%nrank_wa&
    &ke, iermpi => self%mpi_err, nrank => self%myrank)

      if(steps_(1)) then
        call bcswap_wake_step_1_kernel(nx,ny,nz,1,ng,w_swap,wbuf4s_gpu)
      endif
      if(steps_(2)) then
        indy = 1*nx*ng*nz
        if (nrank_wake==nrank) then
          call hipCheck(hipMemcpy(wbuf3r_gpu , wbuf4s_gpu,hipMemcpyDeviceToDevice))
        else
          call mpi_sendrecv(wbuf4s_gpu,indy,mpi_prec,nrank_wake,4,wbuf3r_gpu,indy,mpi_prec,nrank_wake,4,mp_cart,istatus,iermpi)
        endif
      endif
      if(steps_(3)) then
        call bcswap_wake_step_3_kernel(nx,ny,nz,1,ng,w_swap,wbuf3r_gpu,wall_tag_gpu)
      endif
    endassociate
  endsubroutine bcswap_wake_var

  subroutine bcswap_edges_corners(self, steps)
    class(base_amd_object), intent(inout) :: self
    logical, dimension(3), optional :: steps
    logical, dimension(3) :: steps_
    integer :: i,j,k,m, iercuda, indc
    integer, dimension(mpi_status_size) :: istatus

    steps_ = .true. ; if(present(steps)) steps_ = steps

    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv,&
    & w_gpu => self%w_gpu, wbuf1s_c_gpu => self%wbuf1s_c_gpu, wbuf2s_c_gpu => self%wbuf2s_c_gpu,&
    & wbuf3s_c_gpu => self%wbuf3s_c_gpu, wbuf4s_c_gpu => self%wbuf4s_c_gpu, wbuf1r_c_gpu => self%wbuf1r_c&
    &_gpu, wbuf2r_c_gpu => self%wbuf2r_c_gpu, wbuf3r_c_gpu => self%wbuf3r_c_gpu, wbuf4r_c_gpu => self%wbu&
    &f4r_c_gpu, ileftbottom => self%field%ileftbottom, irightbottom => self%field%irightbottom,&
    & ilefttop => self%field%ilefttop, irighttop => self%field%irighttop, pbc => self%field%grid%is_xyz_p&
    &eriodic, mp_cart => self%field%mp_cart, iermpi => self%mpi_err, nrank => self%field%myrank,&
    & lxly => self%field%lxly, lxry => self%field%lxry, rxly => self%field%rxly, rxry => self%field%rxry,&
    &lxlz => self%field%lxlz, lxrz => self%field%lxrz, rxlz => self%field%rxlz, rxrz => self%field%rxrz,&
    &lylz => self%field%lylz, lyrz => self%field%lyrz, rylz => self%field%rylz, ryrz => self%field%ryrz,&
    &lxlylz => self%field%lxlylz, lxlyrz => self%field%lxlyrz, lxrylz => self%field%lxrylz,&
    & rxlylz => self%field%rxlylz,lxryrz => self%field%lxryrz, rxlyrz => self%field%rxlyrz,&
    & rxrylz => self%field%rxrylz, rxryrz => self%field%rxryrz,wbuf_lxly_s_gpu => self%wbuf_lxly_s_gpu ,&
    & wbuf_lxry_s_gpu => self%wbuf_lxry_s_gpu, wbuf_rxly_s_gpu => self%wbuf_rxly_s_gpu ,&
    & wbuf_rxry_s_gpu => self%wbuf_rxry_s_gpu, wbuf_lxlz_s_gpu => self%wbuf_lxlz_s_gpu ,&
    & wbuf_lxrz_s_gpu => self%wbuf_lxrz_s_gpu, wbuf_rxlz_s_gpu => self%wbuf_rxlz_s_gpu ,&
    & wbuf_rxrz_s_gpu => self%wbuf_rxrz_s_gpu, wbuf_lylz_s_gpu => self%wbuf_lylz_s_gpu ,&
    & wbuf_lyrz_s_gpu => self%wbuf_lyrz_s_gpu, wbuf_rylz_s_gpu => self%wbuf_rylz_s_gpu ,&
    & wbuf_ryrz_s_gpu => self%wbuf_ryrz_s_gpu, wbuf_lxlylz_s_gpu => self%wbuf_lxlylz_s_gpu ,&
    & wbuf_lxlyrz_s_gpu => self%wbuf_lxlyrz_s_gpu, wbuf_lxrylz_s_gpu => self%wbuf_lxrylz_s_gpu ,&
    & wbuf_lxryrz_s_gpu => self%wbuf_lxryrz_s_gpu, wbuf_rxlylz_s_gpu => self%wbuf_rxlylz_s_gpu ,&
    & wbuf_rxlyrz_s_gpu => self%wbuf_rxlyrz_s_gpu, wbuf_rxrylz_s_gpu => self%wbuf_rxrylz_s_gpu ,&
    & wbuf_rxryrz_s_gpu => self%wbuf_rxryrz_s_gpu, wbuf_lxly_r_gpu => self%wbuf_lxly_r_gpu ,&
    & wbuf_lxry_r_gpu => self%wbuf_lxry_r_gpu, wbuf_rxly_r_gpu => self%wbuf_rxly_r_gpu ,&
    & wbuf_rxry_r_gpu => self%wbuf_rxry_r_gpu, wbuf_lxlz_r_gpu => self%wbuf_lxlz_r_gpu ,&
    & wbuf_lxrz_r_gpu => self%wbuf_lxrz_r_gpu, wbuf_rxlz_r_gpu => self%wbuf_rxlz_r_gpu ,&
    & wbuf_rxrz_r_gpu => self%wbuf_rxrz_r_gpu, wbuf_lylz_r_gpu => self%wbuf_lylz_r_gpu ,&
    & wbuf_lyrz_r_gpu => self%wbuf_lyrz_r_gpu, wbuf_rylz_r_gpu => self%wbuf_rylz_r_gpu ,&
    & wbuf_ryrz_r_gpu => self%wbuf_ryrz_r_gpu, wbuf_lxlylz_r_gpu => self%wbuf_lxlylz_r_gpu ,&
    & wbuf_lxlyrz_r_gpu => self%wbuf_lxlyrz_r_gpu, wbuf_lxrylz_r_gpu => self%wbuf_lxrylz_r_gpu ,&
    & wbuf_lxryrz_r_gpu => self%wbuf_lxryrz_r_gpu, wbuf_rxlylz_r_gpu => self%wbuf_rxlylz_r_gpu ,&
    & wbuf_rxlyrz_r_gpu => self%wbuf_rxlyrz_r_gpu, wbuf_rxrylz_r_gpu => self%wbuf_rxrylz_r_gpu ,&
    & wbuf_rxryrz_r_gpu => self%wbuf_rxryrz_r_gpu)

      if(steps_(1)) then
        call bcswap_corner_step_1_kernel(nx,ny,nz,nv,ng,w_gpu,wbuf1s_c_gpu,wbuf2s_c_gpu,wbuf3s_c_gpu,wbuf4s_c_gpu)
        if(pbc(2)) then
          call bcswap_corner_step_1b_kernel(nx,ny,nz,nv,ng,w_gpu,wbuf_lxly_s_gpu,wbuf_lxry_s_gpu,&
          &wbuf_rxly_s_gpu,wbuf_rxry_s_gpu,wbuf_lylz_s_gpu,wbuf_lyrz_s_gpu,wbuf_rylz_s_gpu,wbuf_ryrz_s_gpu,&
          &wbuf_lxlylz_s_gpu,wbuf_lxlyrz_s_gpu,wbuf_lxrylz_s_gpu,wbuf_lxryrz_s_gpu,wbuf_rxlylz_s_gpu,&
          &wbuf_rxlyrz_s_gpu,wbuf_rxrylz_s_gpu,wbuf_rxryrz_s_gpu)
        endif
      endif
      if(steps_(2)) then
        indc = nv*ny*ng*ng
        if (ileftbottom == nrank) then
          call hipCheck(hipMemcpy(wbuf2r_c_gpu , wbuf1s_c_gpu,hipMemcpyDeviceToDevice))
          call hipCheck(hipMemcpy(wbuf1r_c_gpu , wbuf2s_c_gpu,hipMemcpyDeviceToDevice))
        else
          call mpi_sendrecv(wbuf1s_c_gpu,indc,mpi_prec,ileftbottom ,1, wbuf2r_c_gpu,indc,mpi_prec,&
          &irighttop ,1,mp_cart,istatus,iermpi)
          call mpi_sendrecv(wbuf2s_c_gpu,indc,mpi_prec,irighttop ,2, wbuf1r_c_gpu,indc,mpi_prec,&
          &ileftbottom ,2,mp_cart,istatus,iermpi)
        endif
        if (ilefttop == nrank) then
          call hipCheck(hipMemcpy(wbuf4r_c_gpu , wbuf3s_c_gpu,hipMemcpyDeviceToDevice))
          call hipCheck(hipMemcpy(wbuf3r_c_gpu , wbuf4s_c_gpu,hipMemcpyDeviceToDevice))
        else
          call mpi_sendrecv(wbuf3s_c_gpu,indc,mpi_prec,ilefttop ,3, wbuf4r_c_gpu,indc,mpi_prec,&
          &irightbottom,3,mp_cart,istatus,iermpi)
          call mpi_sendrecv(wbuf4s_c_gpu,indc,mpi_prec,irightbottom,4, wbuf3r_c_gpu,indc,mpi_prec,&
          &ilefttop ,4,mp_cart,istatus,iermpi)
        endif
        if(pbc(2)) then
          indc = nv*ng*ng*nz
          if (lxly == nrank) then
            call hipCheck(hipMemcpy(wbuf_rxry_r_gpu , wbuf_lxly_s_gpu,hipMemcpyDeviceToDevice))
            call hipCheck(hipMemcpy(wbuf_lxly_r_gpu , wbuf_rxry_s_gpu,hipMemcpyDeviceToDevice))
          else
            call mpi_sendrecv(wbuf_lxly_s_gpu,indc,mpi_prec,lxly ,5, wbuf_rxry_r_gpu,indc,mpi_prec,rxry ,5,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rxry_s_gpu,indc,mpi_prec,rxry ,6, wbuf_lxly_r_gpu,indc,mpi_prec,lxly ,6,mp_cart,istatus,iermpi)
          endif
          if (lxry == nrank) then
            call hipCheck(hipMemcpy(wbuf_rxly_r_gpu , wbuf_lxry_s_gpu,hipMemcpyDeviceToDevice))
            call hipCheck(hipMemcpy(wbuf_lxry_r_gpu , wbuf_rxly_s_gpu,hipMemcpyDeviceToDevice))
          else
            call mpi_sendrecv(wbuf_lxry_s_gpu,indc,mpi_prec,lxry ,7, wbuf_rxly_r_gpu,indc,mpi_prec,rxly ,7,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rxly_s_gpu,indc,mpi_prec,rxly ,8, wbuf_lxry_r_gpu,indc,mpi_prec,lxry ,8,mp_cart,istatus,iermpi)
          endif
          indc = nv*nx*ng*ng
          if (lylz == nrank) then
            call hipCheck(hipMemcpy(wbuf_ryrz_r_gpu , wbuf_lylz_s_gpu,hipMemcpyDeviceToDevice))
            call hipCheck(hipMemcpy(wbuf_lylz_r_gpu , wbuf_ryrz_s_gpu,hipMemcpyDeviceToDevice))
          else
            call mpi_sendrecv(wbuf_lylz_s_gpu,indc,mpi_prec,lylz ,9, wbuf_ryrz_r_gpu,indc,mpi_prec,ryrz ,9,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_ryrz_s_gpu,indc,mpi_prec,ryrz ,10, wbuf_lylz_r_gpu,indc,mpi_prec,&
            &lylz ,10,mp_cart,istatus,iermpi)
          endif
          if (lyrz == nrank) then
            call hipCheck(hipMemcpy(wbuf_rylz_r_gpu , wbuf_lyrz_s_gpu,hipMemcpyDeviceToDevice))
            call hipCheck(hipMemcpy(wbuf_lyrz_r_gpu , wbuf_rylz_s_gpu,hipMemcpyDeviceToDevice))
          else
            call mpi_sendrecv(wbuf_lyrz_s_gpu,indc,mpi_prec,lyrz ,11, wbuf_rylz_r_gpu,indc,mpi_prec,&
            &rylz ,11,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rylz_s_gpu,indc,mpi_prec,rylz ,12, wbuf_lyrz_r_gpu,indc,mpi_prec,&
            &lyrz ,12,mp_cart,istatus,iermpi)
          endif
          indc = nv*ng*ng*ng
          if (lxlylz == nrank) then
            call hipCheck(hipMemcpy(wbuf_rxryrz_r_gpu , wbuf_lxlylz_s_gpu,hipMemcpyDeviceToDevice))
            call hipCheck(hipMemcpy(wbuf_rxrylz_r_gpu , wbuf_lxlyrz_s_gpu,hipMemcpyDeviceToDevice))
            call hipCheck(hipMemcpy(wbuf_rxlyrz_r_gpu , wbuf_lxrylz_s_gpu,hipMemcpyDeviceToDevice))
            call hipCheck(hipMemcpy(wbuf_rxlylz_r_gpu , wbuf_lxryrz_s_gpu,hipMemcpyDeviceToDevice))
            call hipCheck(hipMemcpy(wbuf_lxryrz_r_gpu , wbuf_rxlylz_s_gpu,hipMemcpyDeviceToDevice))
            call hipCheck(hipMemcpy(wbuf_lxrylz_r_gpu , wbuf_rxlyrz_s_gpu,hipMemcpyDeviceToDevice))
            call hipCheck(hipMemcpy(wbuf_lxlyrz_r_gpu , wbuf_rxrylz_s_gpu,hipMemcpyDeviceToDevice))
            call hipCheck(hipMemcpy(wbuf_lxlylz_r_gpu , wbuf_rxryrz_s_gpu,hipMemcpyDeviceToDevice))
          else
            call mpi_sendrecv(wbuf_lxlylz_s_gpu,indc,mpi_prec,lxlylz ,13, wbuf_rxryrz_r_gpu,indc,&
            &mpi_prec,rxryrz ,13,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rxryrz_s_gpu,indc,mpi_prec,rxryrz ,14, wbuf_lxlylz_r_gpu,indc,&
            &mpi_prec,lxlylz ,14,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_lxlyrz_s_gpu,indc,mpi_prec,lxlyrz ,15, wbuf_rxrylz_r_gpu,indc,&
            &mpi_prec,rxrylz ,15,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rxrylz_s_gpu,indc,mpi_prec,rxrylz ,16, wbuf_lxlyrz_r_gpu,indc,&
            &mpi_prec,lxlyrz ,16,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rxlylz_s_gpu,indc,mpi_prec,rxlylz ,17, wbuf_lxryrz_r_gpu,indc,&
            &mpi_prec,lxryrz ,17,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_lxryrz_s_gpu,indc,mpi_prec,lxryrz ,18, wbuf_rxlylz_r_gpu,indc,&
            &mpi_prec,rxlylz ,18,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rxlyrz_s_gpu,indc,mpi_prec,rxlyrz ,19, wbuf_lxrylz_r_gpu,indc,&
            &mpi_prec,lxrylz ,19,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_lxrylz_s_gpu,indc,mpi_prec,lxrylz ,20, wbuf_rxlyrz_r_gpu,indc,&
            &mpi_prec,rxlyrz ,20,mp_cart,istatus,iermpi)
          endif
        endif
      endif
      if(steps_(3)) then
        call bcswap_corner_step_3_kernel(nx,ny,nz,nv,ng,w_gpu,wbuf1r_c_gpu,wbuf2r_c_gpu,&
        &wbuf3r_c_gpu,wbuf4r_c_gpu,ileftbottom,ilefttop,irightbottom,irighttop)
        if(pbc(2)) then
          call bcswap_corner_step_3b_kernel(nx,ny,nz,nv,ng,w_gpu,wbuf_lxly_r_gpu,wbuf_lxry_r_gpu,&
          &wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,wbuf_ryrz_r_gpu,&
          &wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,wbuf_rxlylz_r_gpu,&
          &wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu,lxly,lxry,rxly,rxry,lylz,lyrz,rylz,ryrz,&
          &lxlylz,lxlyrz,lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz)
        endif
      endif
    endassociate
  endsubroutine bcswap_edges_corners

  subroutine bcswap_edges_corners_var(self, w_swap, steps)
    class(base_amd_object), intent(inout) :: self
    real(rkind), dimension(:,:,:,:), target :: w_swap
    logical, dimension(3), optional :: steps
    logical, dimension(3) :: steps_
    integer :: i,j,k,m, iercuda, indc
    integer, dimension(mpi_status_size) :: istatus

    steps_ = .true. ; if(present(steps)) steps_ = steps

    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, nv => self%nv,&
    & wbuf1s_c_gpu => self%wbuf1s_c_gpu, wbuf2s_c_gpu => self%wbuf2s_c_gpu, wbuf3s_c_gpu => self%wbuf3s_c&
    &_gpu, wbuf4s_c_gpu => self%wbuf4s_c_gpu, wbuf1r_c_gpu => self%wbuf1r_c_gpu, wbuf2r_c_gpu => self%wbu&
    &f2r_c_gpu, wbuf3r_c_gpu => self%wbuf3r_c_gpu, wbuf4r_c_gpu => self%wbuf4r_c_gpu, ileftbottom => self&
    &%field%ileftbottom, irightbottom => self%field%irightbottom, ilefttop => self%field%ilefttop,&
    & irighttop => self%field%irighttop, pbc => self%field%grid%is_xyz_periodic, mp_cart => self%field%mp&
    &_cart, iermpi => self%mpi_err, nrank => self%field%myrank, lxly => self%field%lxly,&
    & lxry => self%field%lxry, rxly => self%field%rxly, rxry => self%field%rxry,lxlz => self%field%lxlz,&
    & lxrz => self%field%lxrz, rxlz => self%field%rxlz, rxrz => self%field%rxrz,lylz => self%field%lylz,&
    & lyrz => self%field%lyrz, rylz => self%field%rylz, ryrz => self%field%ryrz,lxlylz => self%field%lxly&
    &lz, lxlyrz => self%field%lxlyrz, lxrylz => self%field%lxrylz, rxlylz => self%field%rxlylz,&
    &lxryrz => self%field%lxryrz, rxlyrz => self%field%rxlyrz, rxrylz => self%field%rxrylz,&
    & rxryrz => self%field%rxryrz,wbuf_lxly_s_gpu => self%wbuf_lxly_s_gpu , wbuf_lxry_s_gpu => self%wbuf_&
    &lxry_s_gpu, wbuf_rxly_s_gpu => self%wbuf_rxly_s_gpu , wbuf_rxry_s_gpu => self%wbuf_rxry_s_gpu,&
    & wbuf_lxlz_s_gpu => self%wbuf_lxlz_s_gpu , wbuf_lxrz_s_gpu => self%wbuf_lxrz_s_gpu,&
    & wbuf_rxlz_s_gpu => self%wbuf_rxlz_s_gpu , wbuf_rxrz_s_gpu => self%wbuf_rxrz_s_gpu,&
    & wbuf_lylz_s_gpu => self%wbuf_lylz_s_gpu , wbuf_lyrz_s_gpu => self%wbuf_lyrz_s_gpu,&
    & wbuf_rylz_s_gpu => self%wbuf_rylz_s_gpu , wbuf_ryrz_s_gpu => self%wbuf_ryrz_s_gpu,&
    & wbuf_lxlylz_s_gpu => self%wbuf_lxlylz_s_gpu , wbuf_lxlyrz_s_gpu => self%wbuf_lxlyrz_s_gpu,&
    & wbuf_lxrylz_s_gpu => self%wbuf_lxrylz_s_gpu , wbuf_lxryrz_s_gpu => self%wbuf_lxryrz_s_gpu,&
    & wbuf_rxlylz_s_gpu => self%wbuf_rxlylz_s_gpu , wbuf_rxlyrz_s_gpu => self%wbuf_rxlyrz_s_gpu,&
    & wbuf_rxrylz_s_gpu => self%wbuf_rxrylz_s_gpu , wbuf_rxryrz_s_gpu => self%wbuf_rxryrz_s_gpu,&
    & wbuf_lxly_r_gpu => self%wbuf_lxly_r_gpu , wbuf_lxry_r_gpu => self%wbuf_lxry_r_gpu,&
    & wbuf_rxly_r_gpu => self%wbuf_rxly_r_gpu , wbuf_rxry_r_gpu => self%wbuf_rxry_r_gpu,&
    & wbuf_lxlz_r_gpu => self%wbuf_lxlz_r_gpu , wbuf_lxrz_r_gpu => self%wbuf_lxrz_r_gpu,&
    & wbuf_rxlz_r_gpu => self%wbuf_rxlz_r_gpu , wbuf_rxrz_r_gpu => self%wbuf_rxrz_r_gpu,&
    & wbuf_lylz_r_gpu => self%wbuf_lylz_r_gpu , wbuf_lyrz_r_gpu => self%wbuf_lyrz_r_gpu,&
    & wbuf_rylz_r_gpu => self%wbuf_rylz_r_gpu , wbuf_ryrz_r_gpu => self%wbuf_ryrz_r_gpu,&
    & wbuf_lxlylz_r_gpu => self%wbuf_lxlylz_r_gpu , wbuf_lxlyrz_r_gpu => self%wbuf_lxlyrz_r_gpu,&
    & wbuf_lxrylz_r_gpu => self%wbuf_lxrylz_r_gpu , wbuf_lxryrz_r_gpu => self%wbuf_lxryrz_r_gpu,&
    & wbuf_rxlylz_r_gpu => self%wbuf_rxlylz_r_gpu , wbuf_rxlyrz_r_gpu => self%wbuf_rxlyrz_r_gpu,&
    & wbuf_rxrylz_r_gpu => self%wbuf_rxrylz_r_gpu , wbuf_rxryrz_r_gpu => self%wbuf_rxryrz_r_gpu)

      if(steps_(1)) then
        call bcswap_corner_step_1_kernel(nx,ny,nz,1,ng,w_swap,wbuf1s_c_gpu,wbuf2s_c_gpu,wbuf3s_c_gpu,wbuf4s_c_gpu)
        if(pbc(2)) then
          call bcswap_corner_step_1b_kernel(nx,ny,nz,1,ng,w_swap,wbuf_lxly_s_gpu,wbuf_lxry_s_gpu,&
          &wbuf_rxly_s_gpu,wbuf_rxry_s_gpu,wbuf_lylz_s_gpu,wbuf_lyrz_s_gpu,wbuf_rylz_s_gpu,wbuf_ryrz_s_gpu,&
          &wbuf_lxlylz_s_gpu,wbuf_lxlyrz_s_gpu,wbuf_lxrylz_s_gpu,wbuf_lxryrz_s_gpu,wbuf_rxlylz_s_gpu,&
          &wbuf_rxlyrz_s_gpu,wbuf_rxrylz_s_gpu,wbuf_rxryrz_s_gpu)
        endif
      endif
      if(steps_(2)) then
        indc = 1*ny*ng*ng
        if (ileftbottom == nrank) then
          call hipCheck(hipMemcpy(wbuf2r_c_gpu , wbuf1s_c_gpu,hipMemcpyDeviceToDevice))
          call hipCheck(hipMemcpy(wbuf1r_c_gpu , wbuf2s_c_gpu,hipMemcpyDeviceToDevice))
        else
          call mpi_sendrecv(wbuf1s_c_gpu,indc,mpi_prec,ileftbottom ,1, wbuf2r_c_gpu,indc,mpi_prec,&
          &irighttop ,1,mp_cart,istatus,iermpi)
          call mpi_sendrecv(wbuf2s_c_gpu,indc,mpi_prec,irighttop ,2, wbuf1r_c_gpu,indc,mpi_prec,&
          &ileftbottom ,2,mp_cart,istatus,iermpi)
        endif
        if (ilefttop == nrank) then
          call hipCheck(hipMemcpy(wbuf4r_c_gpu , wbuf3s_c_gpu,hipMemcpyDeviceToDevice))
          call hipCheck(hipMemcpy(wbuf3r_c_gpu , wbuf4s_c_gpu,hipMemcpyDeviceToDevice))
        else
          call mpi_sendrecv(wbuf3s_c_gpu,indc,mpi_prec,ilefttop ,3, wbuf4r_c_gpu,indc,mpi_prec,&
          &irightbottom,3,mp_cart,istatus,iermpi)
          call mpi_sendrecv(wbuf4s_c_gpu,indc,mpi_prec,irightbottom,4, wbuf3r_c_gpu,indc,mpi_prec,&
          &ilefttop ,4,mp_cart,istatus,iermpi)
        endif
        if(pbc(2)) then
          indc = 1*ng*ng*nz
          if (lxly == nrank) then
            call hipCheck(hipMemcpy(wbuf_rxry_r_gpu , wbuf_lxly_s_gpu,hipMemcpyDeviceToDevice))
            call hipCheck(hipMemcpy(wbuf_lxly_r_gpu , wbuf_rxry_s_gpu,hipMemcpyDeviceToDevice))
          else
            call mpi_sendrecv(wbuf_lxly_s_gpu,indc,mpi_prec,lxly ,5, wbuf_rxry_r_gpu,indc,mpi_prec,rxry ,5,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rxry_s_gpu,indc,mpi_prec,rxry ,6, wbuf_lxly_r_gpu,indc,mpi_prec,lxly ,6,mp_cart,istatus,iermpi)
          endif
          if (lxry == nrank) then
            call hipCheck(hipMemcpy(wbuf_rxly_r_gpu , wbuf_lxry_s_gpu,hipMemcpyDeviceToDevice))
            call hipCheck(hipMemcpy(wbuf_lxry_r_gpu , wbuf_rxly_s_gpu,hipMemcpyDeviceToDevice))
          else
            call mpi_sendrecv(wbuf_lxry_s_gpu,indc,mpi_prec,lxry ,7, wbuf_rxly_r_gpu,indc,mpi_prec,rxly ,7,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rxly_s_gpu,indc,mpi_prec,rxly ,8, wbuf_lxry_r_gpu,indc,mpi_prec,lxry ,8,mp_cart,istatus,iermpi)
          endif
          indc = 1*nx*ng*ng
          if (lylz == nrank) then
            call hipCheck(hipMemcpy(wbuf_ryrz_r_gpu , wbuf_lylz_s_gpu,hipMemcpyDeviceToDevice))
            call hipCheck(hipMemcpy(wbuf_lylz_r_gpu , wbuf_ryrz_s_gpu,hipMemcpyDeviceToDevice))
          else
            call mpi_sendrecv(wbuf_lylz_s_gpu,indc,mpi_prec,lylz ,9, wbuf_ryrz_r_gpu,indc,mpi_prec,ryrz ,9,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_ryrz_s_gpu,indc,mpi_prec,ryrz ,10, wbuf_lylz_r_gpu,indc,mpi_prec,&
            &lylz ,10,mp_cart,istatus,iermpi)
          endif
          if (lyrz == nrank) then
            call hipCheck(hipMemcpy(wbuf_rylz_r_gpu , wbuf_lyrz_s_gpu,hipMemcpyDeviceToDevice))
            call hipCheck(hipMemcpy(wbuf_lyrz_r_gpu , wbuf_rylz_s_gpu,hipMemcpyDeviceToDevice))
          else
            call mpi_sendrecv(wbuf_lyrz_s_gpu,indc,mpi_prec,lyrz ,11, wbuf_rylz_r_gpu,indc,mpi_prec,&
            &rylz ,11,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rylz_s_gpu,indc,mpi_prec,rylz ,12, wbuf_lyrz_r_gpu,indc,mpi_prec,&
            &lyrz ,12,mp_cart,istatus,iermpi)
          endif
          indc = 1*ng*ng*ng
          if (lxlylz == nrank) then
            call hipCheck(hipMemcpy(wbuf_rxryrz_r_gpu , wbuf_lxlylz_s_gpu,hipMemcpyDeviceToDevice))
            call hipCheck(hipMemcpy(wbuf_rxrylz_r_gpu , wbuf_lxlyrz_s_gpu,hipMemcpyDeviceToDevice))
            call hipCheck(hipMemcpy(wbuf_rxlyrz_r_gpu , wbuf_lxrylz_s_gpu,hipMemcpyDeviceToDevice))
            call hipCheck(hipMemcpy(wbuf_rxlylz_r_gpu , wbuf_lxryrz_s_gpu,hipMemcpyDeviceToDevice))
            call hipCheck(hipMemcpy(wbuf_lxryrz_r_gpu , wbuf_rxlylz_s_gpu,hipMemcpyDeviceToDevice))
            call hipCheck(hipMemcpy(wbuf_lxrylz_r_gpu , wbuf_rxlyrz_s_gpu,hipMemcpyDeviceToDevice))
            call hipCheck(hipMemcpy(wbuf_lxlyrz_r_gpu , wbuf_rxrylz_s_gpu,hipMemcpyDeviceToDevice))
            call hipCheck(hipMemcpy(wbuf_lxlylz_r_gpu , wbuf_rxryrz_s_gpu,hipMemcpyDeviceToDevice))
          else
            call mpi_sendrecv(wbuf_lxlylz_s_gpu,indc,mpi_prec,lxlylz ,13, wbuf_rxryrz_r_gpu,indc,&
            &mpi_prec,rxryrz ,13,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rxryrz_s_gpu,indc,mpi_prec,rxryrz ,14, wbuf_lxlylz_r_gpu,indc,&
            &mpi_prec,lxlylz ,14,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_lxlyrz_s_gpu,indc,mpi_prec,lxlyrz ,15, wbuf_rxrylz_r_gpu,indc,&
            &mpi_prec,rxrylz ,15,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rxrylz_s_gpu,indc,mpi_prec,rxrylz ,16, wbuf_lxlyrz_r_gpu,indc,&
            &mpi_prec,lxlyrz ,16,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rxlylz_s_gpu,indc,mpi_prec,rxlylz ,17, wbuf_lxryrz_r_gpu,indc,&
            &mpi_prec,lxryrz ,17,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_lxryrz_s_gpu,indc,mpi_prec,lxryrz ,18, wbuf_rxlylz_r_gpu,indc,&
            &mpi_prec,rxlylz ,18,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_rxlyrz_s_gpu,indc,mpi_prec,rxlyrz ,19, wbuf_lxrylz_r_gpu,indc,&
            &mpi_prec,lxrylz ,19,mp_cart,istatus,iermpi)
            call mpi_sendrecv(wbuf_lxrylz_s_gpu,indc,mpi_prec,lxrylz ,20, wbuf_rxlyrz_r_gpu,indc,&
            &mpi_prec,rxlyrz ,20,mp_cart,istatus,iermpi)
          endif
        endif
      endif
      if(steps_(3)) then
        call bcswap_corner_step_3_kernel(nx,ny,nz,1,ng,w_swap,wbuf1r_c_gpu,wbuf2r_c_gpu,&
        &wbuf3r_c_gpu,wbuf4r_c_gpu,ileftbottom,ilefttop,irightbottom,irighttop)
        if(pbc(2)) then
          call bcswap_corner_step_3b_kernel(nx,ny,nz,1,ng,w_swap,wbuf_lxly_r_gpu,wbuf_lxry_r_gpu,&
          &wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,wbuf_ryrz_r_gpu,&
          &wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,wbuf_rxlylz_r_gpu,&
          &wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu,lxly,lxry,rxly,rxry,lylz,lyrz,rylz,ryrz,&
          &lxlylz,lxlyrz,lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz)
        endif
      endif
    endassociate
  endsubroutine bcswap_edges_corners_var


  subroutine bcswap_var(self, w_swap, steps)
    class(base_amd_object), intent(inout) :: self
    real(rkind), dimension(:,:,:,:), target :: w_swap
    logical, dimension(3), optional :: steps
    logical, dimension(3) :: steps_
    integer :: iercuda, indx, indy, indz
    integer, dimension(mpi_status_size) :: istatus
    integer, dimension(3) :: is_periodic

    steps_ = .true. ; if(present(steps)) steps_ = steps

    is_periodic(:) = 0
    if(self%field%grid%is_xyz_periodic(1)) is_periodic(1) = 1
    if(self%field%grid%is_xyz_periodic(2)) is_periodic(2) = 1
    if(self%field%grid%is_xyz_periodic(3)) is_periodic(3) = 1

    associate(nx => self%nx, ny => self%ny, nz => self%nz, ng => self%ng, wbuf1s_gpu => self%wbuf1s_&
    &gpu, wbuf2s_gpu => self%wbuf2s_gpu, wbuf3s_gpu => self%wbuf3s_gpu, wbuf4s_gpu => self%wbuf4s_gpu,&
    & wbuf5s_gpu => self%wbuf5s_gpu, wbuf6s_gpu => self%wbuf6s_gpu, wbuf1r_gpu => self%wbuf1r_gpu,&
    & wbuf2r_gpu => self%wbuf2r_gpu, wbuf3r_gpu => self%wbuf3r_gpu, wbuf4r_gpu => self%wbuf4r_gpu,&
    & wbuf5r_gpu => self%wbuf5r_gpu, wbuf6r_gpu => self%wbuf6r_gpu, ileftx => self%field%ileftx,&
    &irightx => self%field%irightx, nrank_x => self%field%nrank_x, ilefty => self%field%ilefty,&
    &irighty => self%field%irighty, nrank_y => self%field%nrank_y, ileftz => self%field%ileftz,&
    &irightz => self%field%irightz, nrank_z => self%field%nrank_z, dcsidxnc2_gpu => self%dcsidxnc2_gpu,&
    & dcsidync2_gpu => self%dcsidync2_gpu, detadxnc2_gpu => self%detadxnc2_gpu, detadync2_gpu => self%det&
    &adync2_gpu, dxdcsinc2_gpu => self%dxdcsinc2_gpu, dydcsinc2_gpu => self%dydcsinc2_gpu,&
    & dxdetanc2_gpu => self%dxdetanc2_gpu, dydetanc2_gpu => self%dydetanc2_gpu, mp_cartx => self%field%mp&
    &_cartx, mp_carty => self%field%mp_carty, mp_cartz => self%field%mp_cartz, ndim => self%field%grid%nd&
    &im)

      if(steps_(1)) then
        if(self%field%grid%grid_dim == 1) then
          call bcswap_step_1_kernel(nx,ny,nz,1,ng,ndim,w_swap,wbuf1s_gpu,wbuf2s_gpu,wbuf3s_gpu,wbuf4s_gpu,wbuf5s_gpu,wbuf6s_gpu)
        elseif(self%field%grid%grid_dim == 2) then
          call bcswap_step_1_kernel(nx,ny,nz,1,ng,ndim,w_swap,wbuf1s_gpu,wbuf2s_gpu,wbuf3s_gpu,wbuf4s_gpu,wbuf5s_gpu,wbuf6s_gpu)
        endif
      endif
      if(steps_(2)) then
        indx = ng*ny*nz
        indy = nx*ng*nz
        indz = nx*ny*ng
        if(ileftx == nrank_x) then
          call hipCheck(hipMemcpy(wbuf2r_gpu , wbuf1s_gpu,hipMemcpyDeviceToDevice))
          call hipCheck(hipMemcpy(wbuf1r_gpu , wbuf2s_gpu,hipMemcpyDeviceToDevice))
        else
          call mpi_sendrecv(wbuf1s_gpu,indx,mpi_prec,ileftx ,1,wbuf2r_gpu,indx,mpi_prec,irightx,1,mp_cartx,istatus,self%mpi_err)
          call mpi_sendrecv(wbuf2s_gpu,indx,mpi_prec,irightx,2,wbuf1r_gpu,indx,mpi_prec,ileftx ,2,mp_cartx,istatus,self%mpi_err)
        endif
        if(ilefty == nrank_y) then
          call hipCheck(hipMemcpy(wbuf4r_gpu , wbuf3s_gpu,hipMemcpyDeviceToDevice))
          call hipCheck(hipMemcpy(wbuf3r_gpu , wbuf4s_gpu,hipMemcpyDeviceToDevice))
        else
          call mpi_sendrecv(wbuf3s_gpu,indy,mpi_prec,ilefty ,3,wbuf4r_gpu,indy,mpi_prec,irighty,3,mp_carty,istatus,self%mpi_err)
          call mpi_sendrecv(wbuf4s_gpu,indy,mpi_prec,irighty,4,wbuf3r_gpu,indy,mpi_prec,ilefty ,4,mp_carty,istatus,self%mpi_err)
        endif
        if(ndim == 3) then
          if(ileftz == nrank_z) then
            call hipCheck(hipMemcpy(wbuf6r_gpu , wbuf5s_gpu,hipMemcpyDeviceToDevice))
            call hipCheck(hipMemcpy(wbuf5r_gpu , wbuf6s_gpu,hipMemcpyDeviceToDevice))
          else
            call mpi_sendrecv(wbuf5s_gpu,indz,mpi_prec,ileftz ,5,wbuf6r_gpu,indz,mpi_prec,irightz,5,mp_cartz,istatus,self%mpi_err)
            call mpi_sendrecv(wbuf6s_gpu,indz,mpi_prec,irightz,6,wbuf5r_gpu,indz,mpi_prec,ileftz ,6,mp_cartz,istatus,self%mpi_err)
          endif
        endif
      endif
      if(steps_(3)) then
        if(self%field%grid%grid_dim == 1) then
          call bcswap_step_3_kernel(nx,ny,nz,1,ng,ndim,w_swap,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,&
          &wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,ileftx,ilefty,ileftz,irightx,irighty,irightz)
        elseif(self%field%grid%grid_dim == 2) then
          call bcswap_step_3_kernel(nx,ny,nz,1,ng,ndim,w_swap,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,&
          &wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,ileftx,ilefty,ileftz,irightx,irighty,irightz)
        endif
      endif
    endassociate
  endsubroutine bcswap_var

  subroutine copy_cpu_gpu(self)
    class(base_amd_object), intent(inout) :: self
    integer :: i,j,k,iv

    call hipCheck(hipMemcpy(self%w_gpu , self%field%w,hipMemcpyHostToDevice))

  endsubroutine copy_cpu_gpu


  subroutine copy_gpu_cpu(self)
    class(base_amd_object), intent(inout) :: self
    integer :: i,j,k,iv

    call hipCheck(hipMemcpy(self%field%w , self%w_gpu,hipMemcpyDeviceToHost))

  endsubroutine copy_gpu_cpu

  subroutine initialize(self, field, gpu_bind)
    !< Initialize base backend.
    class(base_amd_object), intent(inout) :: self !< The base backend.
    class(field_object), target :: field
    integer, intent(in) :: gpu_bind
    type(hipDeviceProp_t), target :: device_properties !< Device properties.
    real(rkind) :: device_mem_avail !< Device memory available (Gb).

    self%field => field
    self%nx = self%field%nx
    self%ny = self%field%ny
    self%nz = self%field%nz
    self%ng = self%field%grid%ng
    self%nv = self%field%nv

    call get_mpi_basic_info(self%nprocs, self%myrank, self%masterproc, self%mpi_err)

    call MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, self%local_comm, self%mpi_err)

    if(gpu_bind == 0) then
      call MPI_COMM_RANK(self%local_comm, self%mydev, self%mpi_err)
    else
      self%mydev = 0
    endif
    call hipCheck(hipSetDevice(self%mydev))
    call hipCheck(hipGetDeviceProperties(device_properties, self%mydev))

    device_mem_avail = real(device_properties%totalGlobalMem, rkind)/(1024_rkind**3)
    print '(A,F12.5,A)', ' available device memory ', device_mem_avail, ' Gb'

    call self%alloc()

  endsubroutine initialize

  subroutine check_gpu_mem(self, description)
    class(base_amd_object), intent(inout) :: self !< The base backend.
    character(*) :: description
    integer :: ierr
    integer(c_size_t) :: mem_free, mem_total
    character(128) :: proc_name
    integer :: resultlen
    !call mpi_barrier(mpi_comm_world, ierr)
    call mpi_get_processor_name(proc_name, resultlen, ierr)
    call hipCheck(hipMemGetInfo(mem_free, mem_total))
    !call mpi_barrier(mpi_comm_world, ierr)
    write(error_unit, "(A,2x,A,2x,A,2x,I0,2x,I0,2x,I0)") 'GPU rank,mems: ', description,&
    &proc_name(1:resultlen),self%myrank, mem_free, mem_total
  endsubroutine check_gpu_mem

endmodule streams_base_amd_object


