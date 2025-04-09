module cpumem_mod
  use iso_c_binding
  implicit none

  interface
    subroutine getmemory(mem_total, mem_av) bind(C, name="getmemory")
      import :: C_LONG
      integer(C_LONG), intent(in) :: mem_av, mem_total
    endsubroutine getmemory
  endinterface
endmodule cpumem_mod
