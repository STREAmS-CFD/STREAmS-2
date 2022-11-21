program test_random_f
    use, intrinsic :: iso_fortran_env
    use crandom_f_mod
    real(REAL64) :: x(3)

    call init_crandom_f(1)
    call get_crandom_f(x)
    print*,'x reprod =',x

    call init_crandom_f(1)
    call get_crandom_f(x)
    print*,'x reprod =',x

    call init_crandom_f(1,reproducible=.false.)
    call get_crandom_f(x)
    print*,'x NOT reprod =',x
endprogram test_random_f
