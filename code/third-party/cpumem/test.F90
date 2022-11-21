program main
    use, intrinsic :: iso_fortran_env
    use cpumem_mod

    integer(C_LONG) :: mem_total, mem_av !< CPU memory.
    integer, parameter :: rkind = REAL64
    real(rkind) :: memtotal, memav
    call getmemory(mem_total,mem_av)
    memtotal = real(mem_total,rkind)/(1024_rkind**3)
    memav = real(mem_av,rkind)/(1024_rkind**3)
    write(*,'(A,F6.2,A)') ' total memory ', memtotal, ' Gb'
    write(*,'(A,F6.2,A)') ' available memory ', memav, ' Gb'
endprogram main
