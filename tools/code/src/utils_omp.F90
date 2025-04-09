module utils_omp

  use, intrinsic :: iso_c_binding
  use omp_lib, only : omp_get_default_device, omp_target_alloc, omp_target_free, omp_target_memcpy
  use streams_parameters, only : ikind, ikind64, rkind, byte_size

  implicit none

  private
  public :: omp_target_alloc_f, omp_target_free_f, omp_target_memcpy_f

  interface omp_target_alloc_f
    module procedure omp_target_alloc_f_int_1, omp_target_alloc_f_int_2, omp_target_alloc_f_int_3,&
    & omp_target_alloc_f_int_4, omp_target_alloc_f_int_5, omp_target_alloc_f_int_6, omp_target_alloc_f_in&
    &t_7, omp_target_alloc_f_real_1, omp_target_alloc_f_real_2, omp_target_alloc_f_real_3,&
    & omp_target_alloc_f_real_4, omp_target_alloc_f_real_5, omp_target_alloc_f_real_6,&
    & omp_target_alloc_f_real_7
  endinterface omp_target_alloc_f

  interface omp_target_free_f
    module procedure omp_target_free_f_int_1, omp_target_free_f_int_2, omp_target_free_f_int_3,&
    & omp_target_free_f_int_4, omp_target_free_f_int_5, omp_target_free_f_int_6, omp_target_free_f_int_7,&
    & omp_target_free_f_real_1, omp_target_free_f_real_2, omp_target_free_f_real_3, omp_target_free_f_rea&
    &l_4, omp_target_free_f_real_5, omp_target_free_f_real_6, omp_target_free_f_real_7
  endinterface omp_target_free_f

  interface omp_target_memcpy_f
    module procedure omp_target_memcpy_f_int, omp_target_memcpy_f_real
  endinterface omp_target_memcpy_f

contains

  subroutine omp_target_alloc_f_int_1(fptr_dev, ubounds, omp_dev, ierr, lbounds)
    implicit none
    integer, pointer, intent(out) :: fptr_dev(:)
    integer, intent(in) :: ubounds(1)
    integer, intent(in) :: omp_dev
    integer, intent(in), optional :: lbounds(1)
    integer, intent(out) :: ierr
    integer, pointer :: fptr(:)
    type(c_ptr) :: cptr_dev
    integer(ikind64) :: sizes

    if (present(lbounds)) then
      sizes = ubounds(1) - lbounds(1) + 1
      cptr_dev = omp_target_alloc(int(sizes * byte_size(1_ikind), c_size_t), int(omp_dev, c_int))
      if (c_associated(cptr_dev)) then
        call c_f_pointer(cptr_dev, fptr, shape=[sizes])
        fptr_dev(lbounds(1):) => fptr
        ierr = 0
      else
        fptr_dev => null()
        ierr = 1000
      endif
    else
      cptr_dev = omp_target_alloc(int(int(ubounds(1),ikind64) * byte_size(1_ikind), c_size_t), int(omp_dev, c_int))
      if (c_associated(cptr_dev)) then
        call c_f_pointer(cptr_dev, fptr_dev, shape=[ubounds(1)])
        ierr = 0
      else
        fptr_dev => null()
        ierr = 1000
      endif
    endif

  endsubroutine omp_target_alloc_f_int_1

  subroutine omp_target_alloc_f_int_2(fptr_dev, ubounds, omp_dev, ierr, lbounds)
    implicit none
    integer, pointer, intent(out) :: fptr_dev(:,:)
    integer, intent(in) :: ubounds(2)
    integer, intent(in) :: omp_dev
    integer, intent(in), optional :: lbounds(2)
    integer, intent(out) :: ierr
    integer, pointer :: fptr(:,:)
    type(c_ptr) :: cptr_dev
    integer(ikind64) :: sizes(2)

    if (present(lbounds)) then
      sizes = ubounds - lbounds + 1
      cptr_dev = omp_target_alloc(int(product(sizes) * byte_size(1_ikind), c_size_t), int(omp_dev, c_int))
      if (c_associated(cptr_dev)) then
        call c_f_pointer(cptr_dev, fptr, shape=sizes)
        fptr_dev(lbounds(1):, lbounds(2):) => fptr
        ierr = 0
      else
        fptr_dev => null()
        ierr = 1000
      endif
    else
      cptr_dev = omp_target_alloc(int(product(int(ubounds,ikind64)) * byte_size(1_ikind), c_size_t), int(omp_dev, c_int))
      if (c_associated(cptr_dev)) then
        call c_f_pointer(cptr_dev, fptr_dev, shape=ubounds)
        ierr = 0
      else
        fptr_dev => null()
        ierr = 1000
      endif
    endif

  endsubroutine omp_target_alloc_f_int_2

  subroutine omp_target_alloc_f_int_3(fptr_dev, ubounds, omp_dev, ierr, lbounds)
    implicit none
    integer, pointer, intent(out) :: fptr_dev(:,:,:)
    integer, intent(in) :: ubounds(3)
    integer, intent(in) :: omp_dev
    integer, intent(in), optional :: lbounds(3)
    integer, intent(out) :: ierr
    integer, pointer :: fptr(:,:,:)
    type(c_ptr) :: cptr_dev
    integer(ikind64) :: sizes(3)

    if (present(lbounds)) then
      sizes = ubounds - lbounds + 1
      cptr_dev = omp_target_alloc(int(product(sizes) * byte_size(1_ikind), c_size_t), int(omp_dev, c_int))
      if (c_associated(cptr_dev)) then
        call c_f_pointer(cptr_dev, fptr, shape=sizes)
        fptr_dev(lbounds(1):, lbounds(2):, lbounds(3):) => fptr
        ierr = 0
      else
        fptr_dev => null()
        ierr = 1000
      endif
    else
      cptr_dev = omp_target_alloc(int(product(int(ubounds,ikind64)) * byte_size(1_ikind), c_size_t), int(omp_dev, c_int))
      if (c_associated(cptr_dev)) then
        call c_f_pointer(cptr_dev, fptr_dev, shape=ubounds)
        ierr = 0
      else
        fptr_dev => null()
        ierr = 1000
      endif
    endif

  endsubroutine omp_target_alloc_f_int_3

  subroutine omp_target_alloc_f_int_4(fptr_dev, ubounds, omp_dev, ierr, lbounds)
    implicit none
    integer, pointer, intent(out) :: fptr_dev(:,:,:,:)
    integer, intent(in) :: ubounds(4)
    integer, intent(in) :: omp_dev
    integer, intent(in), optional :: lbounds(4)
    integer, intent(out) :: ierr
    integer, pointer :: fptr(:,:,:,:)
    type(c_ptr) :: cptr_dev
    integer(ikind64) :: sizes(4)

    if (present(lbounds)) then
      sizes = ubounds - lbounds + 1
      cptr_dev = omp_target_alloc(int(product(sizes) * byte_size(1_ikind), c_size_t), int(omp_dev, c_int))
      if (c_associated(cptr_dev)) then
        call c_f_pointer(cptr_dev, fptr, shape=sizes)
        fptr_dev(lbounds(1):, lbounds(2):, lbounds(3):, lbounds(4):) => fptr
        ierr = 0
      else
        fptr_dev => null()
        ierr = 1000
      endif
    else
      cptr_dev = omp_target_alloc(int(product(int(ubounds,ikind64)) * byte_size(1_ikind), c_size_t), int(omp_dev, c_int))
      if (c_associated(cptr_dev)) then
        call c_f_pointer(cptr_dev, fptr_dev, shape=ubounds)
        ierr = 0
      else
        fptr_dev => null()
        ierr = 1000
      endif
    endif

  endsubroutine omp_target_alloc_f_int_4

  subroutine omp_target_alloc_f_int_5(fptr_dev, ubounds, omp_dev, ierr, lbounds)
    implicit none
    integer, pointer, intent(out) :: fptr_dev(:,:,:,:,:)
    integer, intent(in) :: ubounds(5)
    integer, intent(in) :: omp_dev
    integer, intent(in), optional :: lbounds(5)
    integer, intent(out) :: ierr
    integer, pointer :: fptr(:,:,:,:,:)
    type(c_ptr) :: cptr_dev
    integer(ikind64) :: sizes(5)

    if (present(lbounds)) then
      sizes = ubounds - lbounds + 1
      cptr_dev = omp_target_alloc(int(product(sizes) * byte_size(1_ikind), c_size_t), int(omp_dev, c_int))
      if (c_associated(cptr_dev)) then
        call c_f_pointer(cptr_dev, fptr, shape=sizes)
        fptr_dev(lbounds(1):, lbounds(2):, lbounds(3):, lbounds(4):, lbounds(5):) => fptr
        ierr = 0
      else
        fptr_dev => null()
        ierr = 1000
      endif
    else
      cptr_dev = omp_target_alloc(int(product(int(ubounds,ikind64)) * byte_size(1_ikind), c_size_t), int(omp_dev, c_int))
      if (c_associated(cptr_dev)) then
        call c_f_pointer(cptr_dev, fptr_dev, shape=ubounds)
        ierr = 0
      else
        fptr_dev => null()
        ierr = 1000
      endif
    endif

  endsubroutine omp_target_alloc_f_int_5

  subroutine omp_target_alloc_f_int_6(fptr_dev, ubounds, omp_dev, ierr, lbounds)
    implicit none
    integer, pointer, intent(out) :: fptr_dev(:,:,:,:,:,:)
    integer, intent(in) :: ubounds(6)
    integer, intent(in) :: omp_dev
    integer, intent(in), optional :: lbounds(6)
    integer, intent(out) :: ierr
    integer, pointer :: fptr(:,:,:,:,:,:)
    type(c_ptr) :: cptr_dev
    integer(ikind64) :: sizes(6)

    if (present(lbounds)) then
      sizes = ubounds - lbounds + 1
      cptr_dev = omp_target_alloc(int(product(sizes) * byte_size(1_ikind), c_size_t), int(omp_dev, c_int))
      if (c_associated(cptr_dev)) then
        call c_f_pointer(cptr_dev, fptr, shape=sizes)
        fptr_dev(lbounds(1):, lbounds(2):, lbounds(3):, lbounds(4):, lbounds(5):, lbounds(6):) => fptr
        ierr = 0
      else
        fptr_dev => null()
        ierr = 1000
      endif
    else
      cptr_dev = omp_target_alloc(int(product(int(ubounds,ikind64)) * byte_size(1_ikind), c_size_t), int(omp_dev, c_int))
      if (c_associated(cptr_dev)) then
        call c_f_pointer(cptr_dev, fptr_dev, shape=ubounds)
        ierr = 0
      else
        fptr_dev => null()
        ierr = 1000
      endif
    endif

  endsubroutine omp_target_alloc_f_int_6

  subroutine omp_target_alloc_f_int_7(fptr_dev, ubounds, omp_dev, ierr, lbounds)
    implicit none
    integer, pointer, intent(out) :: fptr_dev(:,:,:,:,:,:,:)
    integer, intent(in) :: ubounds(7)
    integer, intent(in) :: omp_dev
    integer, intent(in), optional :: lbounds(7)
    integer, intent(out) :: ierr
    integer, pointer :: fptr(:,:,:,:,:,:,:)
    type(c_ptr) :: cptr_dev
    integer(ikind64) :: sizes(7)

    if (present(lbounds)) then
      sizes = ubounds - lbounds + 1
      cptr_dev = omp_target_alloc(int(product(sizes) * byte_size(1_ikind), c_size_t), int(omp_dev, c_int))
      if (c_associated(cptr_dev)) then
        call c_f_pointer(cptr_dev, fptr, shape=sizes)
        fptr_dev(lbounds(1):, lbounds(2):, lbounds(3):, lbounds(4):, lbounds(5):, lbounds(6):, lbounds(7):) => fptr
        ierr = 0
      else
        fptr_dev => null()
        ierr = 1000
      endif
    else
      cptr_dev = omp_target_alloc(int(product(int(ubounds,ikind64)) * byte_size(1_ikind), c_size_t), int(omp_dev, c_int))
      if (c_associated(cptr_dev)) then
        call c_f_pointer(cptr_dev, fptr_dev, shape=ubounds)
        ierr = 0
      else
        fptr_dev => null()
        ierr = 1000
      endif
    endif

  endsubroutine omp_target_alloc_f_int_7

  subroutine omp_target_alloc_f_real_1(fptr_dev, ubounds, omp_dev, ierr, lbounds)
    implicit none
    real(rkind), pointer, intent(out) :: fptr_dev(:)
    integer, intent(in) :: ubounds(1)
    integer, intent(in), optional :: omp_dev
    integer, intent(in), optional :: lbounds(1)
    integer, intent(out) :: ierr
    real(rkind), pointer :: fptr(:)
    type(c_ptr) :: cptr_dev
    integer(ikind64) :: sizes

    if (present(lbounds)) then
      sizes = ubounds(1) - lbounds(1) + 1
      cptr_dev = omp_target_alloc(int(sizes * byte_size(1._rkind), c_size_t), int(omp_dev, c_int))
      if (c_associated(cptr_dev)) then
        call c_f_pointer(cptr_dev, fptr, shape=[sizes])
        fptr_dev(lbounds(1):) => fptr
        ierr = 0
      else
        fptr_dev => null()
        ierr = 1000
      endif
    else
      cptr_dev = omp_target_alloc(int(int(ubounds(1),ikind64) * byte_size(1._rkind), c_size_t), int(omp_dev, c_int))
      if (c_associated(cptr_dev)) then
        call c_f_pointer(cptr_dev, fptr_dev, shape=[ubounds(1)])
        ierr = 0
      else
        fptr_dev => null()
        ierr = 1000
      endif
    endif

  endsubroutine omp_target_alloc_f_real_1

  subroutine omp_target_alloc_f_real_2(fptr_dev, ubounds, omp_dev, ierr, lbounds)
    implicit none
    real(rkind), pointer, intent(out) :: fptr_dev(:,:)
    integer, intent(in) :: ubounds(2)
    integer, intent(in), optional :: omp_dev
    integer, intent(in), optional :: lbounds(2)
    integer, intent(out) :: ierr
    real(rkind), pointer :: fptr(:,:)
    type(c_ptr) :: cptr_dev
    integer(ikind64) :: sizes(2)

    if (present(lbounds)) then
      sizes = ubounds - lbounds + 1
      cptr_dev = omp_target_alloc(int(product(sizes) * byte_size(1._rkind), c_size_t), int(omp_dev, c_int))
      if (c_associated(cptr_dev)) then
        call c_f_pointer(cptr_dev, fptr, shape=sizes)
        fptr_dev(lbounds(1):, lbounds(2):) => fptr
        ierr = 0
      else
        fptr_dev => null()
        ierr = 1000
      endif
    else
      cptr_dev = omp_target_alloc(int(product(int(ubounds,ikind64)) * byte_size(1._rkind), c_size_t), int(omp_dev, c_int))
      if (c_associated(cptr_dev)) then
        call c_f_pointer(cptr_dev, fptr_dev, shape=ubounds)
        ierr = 0
      else
        fptr_dev => null()
        ierr = 1000
      endif
    endif

  endsubroutine omp_target_alloc_f_real_2

  subroutine omp_target_alloc_f_real_3(fptr_dev, ubounds, omp_dev, ierr, lbounds)
    implicit none
    real(rkind), pointer, intent(out) :: fptr_dev(:,:,:)
    integer, intent(in) :: ubounds(3)
    integer, intent(in), optional :: omp_dev
    integer, intent(in), optional :: lbounds(3)
    integer, intent(out) :: ierr
    real(rkind), pointer :: fptr(:,:,:)
    type(c_ptr) :: cptr_dev
    integer(ikind64) :: sizes(3)

    if (present(lbounds)) then
      sizes = ubounds - lbounds + 1
      cptr_dev = omp_target_alloc(int(product(sizes) * byte_size(1._rkind), c_size_t), int(omp_dev, c_int))
      if (c_associated(cptr_dev)) then
        call c_f_pointer(cptr_dev, fptr, shape=sizes)
        fptr_dev(lbounds(1):, lbounds(2):, lbounds(3):) => fptr
        ierr = 0
      else
        fptr_dev => null()
        ierr = 1000
      endif
    else
      cptr_dev = omp_target_alloc(int(product(int(ubounds,ikind64)) * byte_size(1._rkind), c_size_t), int(omp_dev, c_int))
      if (c_associated(cptr_dev)) then
        call c_f_pointer(cptr_dev, fptr_dev, shape=ubounds)
        ierr = 0
      else
        fptr_dev => null()
        ierr = 1000
      endif
    endif

  endsubroutine omp_target_alloc_f_real_3

  subroutine omp_target_alloc_f_real_4(fptr_dev, ubounds, omp_dev, ierr, lbounds)
    implicit none
    real(rkind), pointer, intent(out) :: fptr_dev(:,:,:,:)
    integer, intent(in) :: ubounds(4)
    integer, intent(in), optional :: omp_dev
    integer, intent(in), optional :: lbounds(4)
    integer, intent(out) :: ierr
    real(rkind), pointer :: fptr(:,:,:,:)
    type(c_ptr) :: cptr_dev
    integer(ikind64) :: sizes(4)

    if (present(lbounds)) then
      sizes = ubounds - lbounds + 1
      cptr_dev = omp_target_alloc(int(product(sizes) * byte_size(1._rkind), c_size_t), int(omp_dev, c_int))
      if (c_associated(cptr_dev)) then
        call c_f_pointer(cptr_dev, fptr, shape=sizes)
        fptr_dev(lbounds(1):, lbounds(2):, lbounds(3):, lbounds(4):) => fptr
        ierr = 0
      else
        fptr_dev => null()
        ierr = 1000
      endif
    else
      cptr_dev = omp_target_alloc(int(product(int(ubounds,ikind64)) * byte_size(1._rkind), c_size_t), int(omp_dev, c_int))
      if (c_associated(cptr_dev)) then
        call c_f_pointer(cptr_dev, fptr_dev, shape=ubounds)
        ierr = 0
      else
        fptr_dev => null()
        ierr = 1000
      endif
    endif

  endsubroutine omp_target_alloc_f_real_4

  subroutine omp_target_alloc_f_real_5(fptr_dev, ubounds, omp_dev, ierr, lbounds)
    implicit none
    real(rkind), pointer, intent(out) :: fptr_dev(:,:,:,:,:)
    integer, intent(in) :: ubounds(5)
    integer, intent(in), optional :: omp_dev
    integer, intent(in), optional :: lbounds(5)
    integer, intent(out) :: ierr
    real(rkind), pointer :: fptr(:,:,:,:,:)
    type(c_ptr) :: cptr_dev
    integer(ikind64) :: sizes(5)

    if (present(lbounds)) then
      sizes = ubounds - lbounds + 1
      cptr_dev = omp_target_alloc(int(product(sizes) * byte_size(1._rkind), c_size_t), int(omp_dev, c_int))
      if (c_associated(cptr_dev)) then
        call c_f_pointer(cptr_dev, fptr, shape=sizes)
        fptr_dev(lbounds(1):, lbounds(2):, lbounds(3):, lbounds(4):, lbounds(5):) => fptr
        ierr = 0
      else
        fptr_dev => null()
        ierr = 1000
      endif
    else
      cptr_dev = omp_target_alloc(int(product(int(ubounds,ikind64)) * byte_size(1._rkind), c_size_t), int(omp_dev, c_int))
      if (c_associated(cptr_dev)) then
        call c_f_pointer(cptr_dev, fptr_dev, shape=ubounds)
        ierr = 0
      else
        fptr_dev => null()
        ierr = 1000
      endif
    endif

  endsubroutine omp_target_alloc_f_real_5

  subroutine omp_target_alloc_f_real_6(fptr_dev, ubounds, omp_dev, ierr, lbounds)
    implicit none
    real(rkind), pointer, intent(out) :: fptr_dev(:,:,:,:,:,:)
    integer, intent(in) :: ubounds(6)
    integer, intent(in), optional :: omp_dev
    integer, intent(in), optional :: lbounds(6)
    integer, intent(out) :: ierr
    real(rkind), pointer :: fptr(:,:,:,:,:,:)
    type(c_ptr) :: cptr_dev
    integer(ikind64) :: sizes(6)

    if (present(lbounds)) then
      sizes = ubounds - lbounds + 1
      cptr_dev = omp_target_alloc(int(product(sizes) * byte_size(1._rkind), c_size_t), int(omp_dev, c_int))
      if (c_associated(cptr_dev)) then
        call c_f_pointer(cptr_dev, fptr, shape=sizes)
        fptr_dev(lbounds(1):, lbounds(2):, lbounds(3):, lbounds(4):, lbounds(5):, lbounds(6):) => fptr
        ierr = 0
      else
        fptr_dev => null()
        ierr = 1000
      endif
    else
      cptr_dev = omp_target_alloc(int(product(int(ubounds,ikind64)) * byte_size(1._rkind), c_size_t), int(omp_dev, c_int))
      if (c_associated(cptr_dev)) then
        call c_f_pointer(cptr_dev, fptr_dev, shape=ubounds)
        ierr = 0
      else
        fptr_dev => null()
        ierr = 1000
      endif
    endif

  endsubroutine omp_target_alloc_f_real_6

  subroutine omp_target_alloc_f_real_7(fptr_dev, ubounds, omp_dev, ierr, lbounds)
    implicit none
    real(rkind), pointer, intent(out) :: fptr_dev(:,:,:,:,:,:,:)
    integer, intent(in) :: ubounds(7)
    integer, intent(in), optional :: omp_dev
    integer, intent(in), optional :: lbounds(7)
    integer, intent(out) :: ierr
    real(rkind), pointer :: fptr(:,:,:,:,:,:,:)
    type(c_ptr) :: cptr_dev
    integer(ikind64) :: sizes(7)

    if (present(lbounds)) then
      sizes = ubounds - lbounds + 1
      cptr_dev = omp_target_alloc(int(product(sizes) * byte_size(1._rkind), c_size_t), int(omp_dev, c_int))
      if (c_associated(cptr_dev)) then
        call c_f_pointer(cptr_dev, fptr, shape=sizes)
        fptr_dev(lbounds(1):, lbounds(2):, lbounds(3):, lbounds(4):, lbounds(5):, lbounds(6):, lbounds(7):) => fptr
        ierr = 0
      else
        fptr_dev => null()
        ierr = 1000
      endif
    else
      cptr_dev = omp_target_alloc(int(product(int(ubounds,ikind64)) * byte_size(1._rkind), c_size_t), int(omp_dev, c_int))
      if (c_associated(cptr_dev)) then
        call c_f_pointer(cptr_dev, fptr_dev, shape=ubounds)
        ierr = 0
      else
        fptr_dev => null()
        ierr = 1000
      endif
    endif

  endsubroutine omp_target_alloc_f_real_7

  subroutine omp_target_free_f_int_1(fptr_dev, omp_dev)
    implicit none
    integer, pointer, intent(out) :: fptr_dev(:)
    integer, intent(in) :: omp_dev
    type(c_ptr) :: cptr_dev

    cptr_dev = c_loc(fptr_dev)

    call omp_target_free(cptr_dev, int(omp_dev, c_int))

    nullify(fptr_dev)
  endsubroutine omp_target_free_f_int_1

  subroutine omp_target_free_f_int_2(fptr_dev, omp_dev)
    implicit none
    integer, pointer, intent(out) :: fptr_dev(:,:)
    integer, intent(in) :: omp_dev
    type(c_ptr) :: cptr_dev

    cptr_dev = c_loc(fptr_dev)

    call omp_target_free(cptr_dev, int(omp_dev, c_int))

    nullify(fptr_dev)
  endsubroutine omp_target_free_f_int_2

  subroutine omp_target_free_f_int_3(fptr_dev, omp_dev)
    implicit none
    integer, pointer, intent(out) :: fptr_dev(:,:,:)
    integer, intent(in) :: omp_dev
    type(c_ptr) :: cptr_dev

    cptr_dev = c_loc(fptr_dev)

    call omp_target_free(cptr_dev, int(omp_dev, c_int))

    nullify(fptr_dev)
  endsubroutine omp_target_free_f_int_3

  subroutine omp_target_free_f_int_4(fptr_dev, omp_dev)
    implicit none
    integer, pointer, intent(out) :: fptr_dev(:,:,:,:)
    integer, intent(in) :: omp_dev
    type(c_ptr) :: cptr_dev

    cptr_dev = c_loc(fptr_dev)

    call omp_target_free(cptr_dev, int(omp_dev, c_int))

    nullify(fptr_dev)
  endsubroutine omp_target_free_f_int_4

  subroutine omp_target_free_f_int_5(fptr_dev, omp_dev)
    implicit none
    integer, pointer, intent(out) :: fptr_dev(:,:,:,:,:)
    integer, intent(in) :: omp_dev
    type(c_ptr) :: cptr_dev

    cptr_dev = c_loc(fptr_dev)

    call omp_target_free(cptr_dev, int(omp_dev, c_int))

    nullify(fptr_dev)
  endsubroutine omp_target_free_f_int_5

  subroutine omp_target_free_f_int_6(fptr_dev, omp_dev)
    implicit none
    integer, pointer, intent(out) :: fptr_dev(:,:,:,:,:,:)
    integer, intent(in) :: omp_dev
    type(c_ptr) :: cptr_dev

    cptr_dev = c_loc(fptr_dev)

    call omp_target_free(cptr_dev, int(omp_dev, c_int))

    nullify(fptr_dev)
  endsubroutine omp_target_free_f_int_6

  subroutine omp_target_free_f_int_7(fptr_dev, omp_dev)
    implicit none
    integer, pointer, intent(out) :: fptr_dev(:,:,:,:,:,:,:)
    integer, intent(in) :: omp_dev
    type(c_ptr) :: cptr_dev

    cptr_dev = c_loc(fptr_dev)

    call omp_target_free(cptr_dev, int(omp_dev, c_int))

    nullify(fptr_dev)
  endsubroutine omp_target_free_f_int_7

  subroutine omp_target_free_f_real_1(fptr_dev, omp_dev)
    implicit none
    real(rkind), pointer, intent(out) :: fptr_dev(:)
    integer, intent(in), optional :: omp_dev
    type(c_ptr) :: cptr_dev

    cptr_dev = c_loc(fptr_dev)

    call omp_target_free(cptr_dev, int(omp_dev, c_int))

    nullify(fptr_dev)
  endsubroutine omp_target_free_f_real_1

  subroutine omp_target_free_f_real_2(fptr_dev, omp_dev)
    implicit none
    real(rkind), pointer, intent(out) :: fptr_dev(:,:)
    integer, intent(in), optional :: omp_dev
    type(c_ptr) :: cptr_dev

    cptr_dev = c_loc(fptr_dev)

    call omp_target_free(cptr_dev, int(omp_dev, c_int))

    nullify(fptr_dev)
  endsubroutine omp_target_free_f_real_2

  subroutine omp_target_free_f_real_3(fptr_dev, omp_dev)
    implicit none
    real(rkind), pointer, intent(out) :: fptr_dev(:,:,:)
    integer, intent(in), optional :: omp_dev
    type(c_ptr) :: cptr_dev

    cptr_dev = c_loc(fptr_dev)

    call omp_target_free(cptr_dev, int(omp_dev, c_int))

    nullify(fptr_dev)
  endsubroutine omp_target_free_f_real_3

  subroutine omp_target_free_f_real_4(fptr_dev, omp_dev)
    implicit none
    real(rkind), pointer, intent(out) :: fptr_dev(:,:,:,:)
    integer, intent(in), optional :: omp_dev
    type(c_ptr) :: cptr_dev

    cptr_dev = c_loc(fptr_dev)

    call omp_target_free(cptr_dev, int(omp_dev, c_int))

    nullify(fptr_dev)
  endsubroutine omp_target_free_f_real_4

  subroutine omp_target_free_f_real_5(fptr_dev, omp_dev)
    implicit none
    real(rkind), pointer, intent(out) :: fptr_dev(:,:,:,:,:)
    integer, intent(in), optional :: omp_dev
    type(c_ptr) :: cptr_dev

    cptr_dev = c_loc(fptr_dev)

    call omp_target_free(cptr_dev, int(omp_dev, c_int))

    nullify(fptr_dev)
  endsubroutine omp_target_free_f_real_5

  subroutine omp_target_free_f_real_6(fptr_dev, omp_dev)
    implicit none
    real(rkind), pointer, intent(out) :: fptr_dev(:,:,:,:,:,:)
    integer, intent(in), optional :: omp_dev
    type(c_ptr) :: cptr_dev

    cptr_dev = c_loc(fptr_dev)

    call omp_target_free(cptr_dev, int(omp_dev, c_int))

    nullify(fptr_dev)
  endsubroutine omp_target_free_f_real_6

  subroutine omp_target_free_f_real_7(fptr_dev, omp_dev)
    implicit none
    real(rkind), pointer, intent(out) :: fptr_dev(:,:,:,:,:,:,:)
    integer, intent(in), optional :: omp_dev
    type(c_ptr) :: cptr_dev

    cptr_dev = c_loc(fptr_dev)

    call omp_target_free(cptr_dev, int(omp_dev, c_int))

    nullify(fptr_dev)
  endsubroutine omp_target_free_f_real_7

  function omp_target_memcpy_f_int(fptr_dst, fptr_src, dst_off, src_off, omp_dst_dev, omp_src_dev)
    implicit none
    integer(ikind) :: omp_target_memcpy_f_int
    integer, target, intent(out) :: fptr_dst(..)
    integer, target, intent(in) :: fptr_src(..)
    integer, intent(in) :: omp_dst_dev, omp_src_dev
    integer, intent(in) :: dst_off, src_off
    integer(ikind64) :: n_elements
    integer(c_size_t) :: total_dim, omp_dst_offset, omp_src_offset
    type(c_ptr) :: cptr_dst, cptr_src
    integer(c_int) :: omp_dst_device, omp_src_device

    n_elements = size(fptr_src,kind=ikind64)

    omp_dst_offset = int(dst_off, c_size_t)
    omp_src_offset = int(src_off, c_size_t)
    omp_dst_device = int(omp_dst_dev, c_int)
    omp_src_device = int(omp_src_dev, c_int)

    cptr_dst = c_loc(fptr_dst)
    cptr_src = c_loc(fptr_src)

    total_dim = int(n_elements * byte_size(1_ikind), c_size_t)

    omp_target_memcpy_f_int = int(omp_target_memcpy(cptr_dst, cptr_src, total_dim, omp_dst_offset,&
    & omp_src_offset, omp_dst_device, omp_src_device), ikind)
  endfunction omp_target_memcpy_f_int

  function omp_target_memcpy_f_real(fptr_dst, fptr_src, dst_off, src_off, omp_dst_dev, omp_src_dev)
    implicit none
    integer(ikind) :: omp_target_memcpy_f_real
    real(rkind), target, intent(out) :: fptr_dst(..)
    real(rkind), target, intent(in) :: fptr_src(..)
    integer, intent(in) :: omp_dst_dev, omp_src_dev
    integer, intent(in) :: dst_off, src_off
    integer(ikind64) :: n_elements
    integer(c_size_t) :: total_dim, omp_dst_offset, omp_src_offset
    type(c_ptr) :: cptr_dst, cptr_src
    integer(c_int) :: omp_dst_device, omp_src_device

    n_elements = size(fptr_src,kind=ikind64)

    omp_dst_offset = int(dst_off, c_size_t)
    omp_src_offset = int(src_off, c_size_t)
    omp_dst_device = int(omp_dst_dev, c_int)
    omp_src_device = int(omp_src_dev, c_int)

    cptr_dst = c_loc(fptr_dst)
    cptr_src = c_loc(fptr_src)

    total_dim = int(n_elements * byte_size(1._rkind), c_size_t)

    omp_target_memcpy_f_real = int(omp_target_memcpy(cptr_dst, cptr_src, total_dim, omp_dst_offset,&
    & omp_src_offset, omp_dst_device, omp_src_device), ikind)
  endfunction omp_target_memcpy_f_real

endmodule utils_omp

