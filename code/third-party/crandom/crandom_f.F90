module crandom_f_mod

   use iso_c_binding
   implicit none

   interface
      subroutine init_crandom(seed) bind(C, name="init_crandom")
         import :: C_INT
         integer(C_INT), intent(in), value :: seed
      endsubroutine init_crandom
   endinterface

   interface
      subroutine get_crandom(a) bind(C, name="get_crandom")
         import :: C_DOUBLE
         real(C_DOUBLE), intent(out) :: a
      endsubroutine get_crandom
   endinterface

   interface get_crandom_f
       module procedure get_crandom_f_scalar, get_crandom_f_array1d, &
                        get_crandom_f_array2d, get_crandom_f_array3d
   endinterface get_crandom_f

   contains

   subroutine init_crandom_f(num, reproducible)
       integer(C_INT) :: num
       integer(C_INT) :: start_seed
       integer :: clock
       integer, dimension(8) :: values
       logical, optional :: reproducible
       logical :: reproducible_

       reproducible_ = .true. ; if(present(reproducible)) reproducible_ = reproducible

       if(reproducible_) then
           start_seed = num
       else
           call date_and_time(VALUES=values)
           clock = (values(6)+1)*(values(7)+1)*(values(8)+1)
           start_seed = int(real(clock,C_DOUBLE)+num,C_INT) !+ (x+1) * (/ (i - 1, i = 1, n) /)
       endif
       !print*,'start_seed = ',start_seed
       call init_crandom(start_seed)
   endsubroutine init_crandom_f

   subroutine get_crandom_f_scalar(a)
       real(C_DOUBLE) :: a
       call get_crandom(a)
   endsubroutine get_crandom_f_scalar

   subroutine get_crandom_f_array1d(a)
       real(C_DOUBLE), dimension(:) :: a
       real(C_DOUBLE) :: a_sca
       integer :: l
       do l=1,size(a)
           call get_crandom(a_sca)
           a(l) = a_sca
       enddo
   endsubroutine get_crandom_f_array1d

   subroutine get_crandom_f_array2d(a)
       real(C_DOUBLE), dimension(:,:) :: a
       real(C_DOUBLE) :: a_sca
       integer :: l, m
       do l=1,size(a,1)
       do m=1,size(a,2)
           call get_crandom(a_sca)
           a(l,m) = a_sca
       enddo
       enddo
   endsubroutine get_crandom_f_array2d

   subroutine get_crandom_f_array3d(a)
       real(C_DOUBLE), dimension(:,:,:) :: a
       real(C_DOUBLE) :: a_sca
       integer :: l, m, n
       do l=1,size(a,1)
       do m=1,size(a,2)
       do n=1,size(a,3)
           call get_crandom(a_sca)
           a(l,m,n) = a_sca
       enddo
       enddo
       enddo
   endsubroutine get_crandom_f_array3d

endmodule crandom_f_mod
