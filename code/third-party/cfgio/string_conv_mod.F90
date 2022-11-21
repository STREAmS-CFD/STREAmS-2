! Program :
! Author  : wansooha@gmail.com
! Date    :

module string_conv_mod

    private
    integer,parameter:: MXNSTR=256
    integer, parameter:: sp=kind(0.0)
    integer, parameter:: dp=kind(0.d0)

    ! subroutines
    interface from_string
        module procedure from_string_i
        module procedure from_string_f
        module procedure from_string_d
        module procedure from_string_c
        module procedure from_string_z
        module procedure from_string_b
    end interface
    interface to_string
        module procedure to_string_i
        module procedure to_string_f
        module procedure to_string_d
        module procedure to_string_c
        module procedure to_string_z
        module procedure to_string_b
    end interface

    ! functions
    interface tostr
        module procedure i2s
        module procedure f2s
        module procedure d2s
        module procedure c2s
        module procedure z2s
        module procedure b2s
        module procedure s2s
    end interface
    interface tolist
        module procedure to_list_i
        module procedure to_list_f
        module procedure to_list_d
        module procedure to_list_c
        module procedure to_list_z
        module procedure to_list_b
        module procedure to_list_s
    end interface

    public:: from_string,to_string
    public:: tostr,tolist
    public:: list_size,list_size_cmplx,list_size_str
    public:: quote,unquote
    public:: s2i,s2f,s2d,s2c,s2z,s2b

contains

    pure function quote(str) result(v)
    character(len=*),intent(in):: str
    character(len=:),allocatable:: v
    if(len_trim(str)==0) then
        v='""'
        return
    endif
    if(str(1:1)=='"'.or.str(1:1)=="'") then
        v=str
    else
        v='"'//trim(str)//'"'
    endif
    end function

    pure function unquote(str) result(v)
    character(len=*),intent(in):: str
    character(len=:),allocatable:: v
    integer:: l
    if(str(1:1)=='"' .or. str(1:1)=="'") then
        l=len_trim(str)
        v=str(2:l-1)
    else
        v=str
    endif
    end function

! string to other types
    pure function s2i(str) result(v)
    character(len=*),intent(in):: str
    integer:: v
    read(str,*) v
    end function
    pure function s2f(str) result(v)
    character(len=*),intent(in):: str
    real(sp):: v
    read(str,*) v
    end function
    pure function s2d(str) result(v)
    character(len=*),intent(in):: str
    real(dp):: v
    read(str,*) v
    end function
    pure function s2c(str) result(v)
    character(len=*),intent(in):: str
    complex(sp):: v
    read(str,*) v
    end function
    pure function s2z(str) result(v)
    character(len=*),intent(in):: str
    complex(dp):: v
    read(str,*) v
    end function
    pure function s2b(str) result(v)
    character(len=*),intent(in):: str
    logical:: v
    select case(str)
    case("TRUE","True","true","T","t")
        v=.true.
    case("YES","Yes","yes","Y","y")
        v=.true.
    case("ON","On","on")
        v=.true.
    case(".TRUE.",".true.")
        v=.true.
    case default
        v=.false.
    end select
    end function

    subroutine from_string_i(str,v)
    integer,intent(out):: v
    character(len=*),intent(in):: str
    read(str,*) v
    end subroutine
    subroutine from_string_f(str,v)
    real(sp),intent(out):: v
    character(len=*),intent(in):: str
    read(str,*) v
    end subroutine
    subroutine from_string_d(str,v)
    real(dp),intent(out):: v
    character(len=*),intent(in):: str
    read(str,*) v
    end subroutine
    subroutine from_string_c(str,v)
    complex(sp),intent(out):: v
    character(len=*),intent(in):: str
    read(str,*) v
    end subroutine
    subroutine from_string_z(str,v)
    complex(dp),intent(out):: v
    character(len=*),intent(in):: str
    read(str,*) v
    end subroutine
    subroutine from_string_b(str,v)
    logical,intent(out):: v
    character(len=*),intent(in):: str
    v=s2b(str)
    end subroutine

! to string
    pure function i2s(v) result(str)
    integer,intent(in):: v
    character(len=MXNSTR):: str
    write(str,*) v
    end function
    pure function f2s(v) result(str)
    real(sp),intent(in):: v
    character(len=MXNSTR):: str
    write(str,*) v
    end function
    pure function d2s(v) result(str)
    real(dp),intent(in):: v
    character(len=MXNSTR):: str
    write(str,*) v
    end function
    pure function c2s(v) result(str)
    complex(sp),intent(in):: v
    character(len=MXNSTR):: str
    write(str,*) v
    end function
    pure function z2s(v) result(str)
    complex(dp),intent(in):: v
    character(len=MXNSTR):: str
    write(str,*) v
    end function
    pure function b2s(v) result(str)
    logical,intent(in):: v
    character(len=MXNSTR):: str
    write(str,*) v
    end function
    pure function s2s(v) result(str)
    character(len=*),intent(in):: v
    character(len=MXNSTR):: str
    write(str,*) v
    end function

    subroutine to_string_i(v,str)
    integer,intent(in):: v
    character(len=*),intent(out):: str
    write(str,*) v
    end subroutine
    subroutine to_string_f(v,str)
    real(sp),intent(in):: v
    character(len=*),intent(out):: str
    write(str,*) v
    end subroutine
    subroutine to_string_d(v,str)
    real(dp),intent(in):: v
    character(len=*),intent(out):: str
    write(str,*) v
    end subroutine
    subroutine to_string_c(v,str)
    complex(sp),intent(in):: v
    character(len=*),intent(out):: str
    write(str,*) v
    end subroutine
    subroutine to_string_z(v,str)
    complex(dp),intent(in):: v
    character(len=*),intent(out):: str
    write(str,*) v
    end subroutine
    subroutine to_string_b(v,str)
    logical,intent(in):: v
    character(len=*),intent(out):: str
    write(str,*) v
    end subroutine

! array to string
    function to_list_i(arr) result(str)
    integer,intent(in):: arr(:)
    character(len=MXNSTR):: str
    character:: sep
    integer i
    str=''
    do i=1,size(arr)
        if(i==1) then
            sep=''
        else
            sep=','
        endif
        str=trim(adjustl(str))//sep//trim(adjustl(tostr(arr(i))))
    enddo
    end function
    function to_list_f(arr) result(str)
    real(sp),intent(in):: arr(:)
    character(len=MXNSTR):: str
    character:: sep
    integer i
    str=''
    do i=1,size(arr)
        if(i==1) then
            sep=''
        else
            sep=','
        endif
        str=trim(adjustl(str))//sep//trim(adjustl(tostr(arr(i))))
    enddo
    end function
    function to_list_d(arr) result(str)
    real(dp),intent(in):: arr(:)
    character(len=MXNSTR):: str
    character:: sep
    integer i
    str=''
    do i=1,size(arr)
        if(i==1) then
            sep=''
        else
            sep=','
        endif
        str=trim(adjustl(str))//sep//trim(adjustl(tostr(arr(i))))
    enddo
    end function
    function to_list_c(arr) result(str)
    complex(sp),intent(in):: arr(:)
    character(len=MXNSTR):: str
    character:: sep
    integer i
    str=''
    do i=1,size(arr)
        if(i==1) then
            sep=''
        else
            sep=','
        endif
        str=trim(adjustl(str))//sep//trim(adjustl(tostr(arr(i))))
    enddo
    end function
    function to_list_z(arr) result(str)
    complex(dp),intent(in):: arr(:)
    character(len=MXNSTR):: str
    character:: sep
    integer i
    str=''
    do i=1,size(arr)
        if(i==1) then
            sep=''
        else
            sep=','
        endif
        str=trim(adjustl(str))//sep//trim(adjustl(tostr(arr(i))))
    enddo
    end function
    function to_list_b(arr) result(str)
    logical,intent(in):: arr(:)
    character(len=MXNSTR):: str
    character:: sep
    integer i
    str=''
    do i=1,size(arr)
        if(i==1) then
            sep=''
        else
            sep=','
        endif
        str=trim(adjustl(str))//sep//trim(adjustl(tostr(arr(i))))
    enddo
    end function
    function to_list_s(arr) result(str)
    character(len=*),intent(in):: arr(:)
    character(len=MXNSTR):: str
    character:: sep
    integer i
    str=''
    do i=1,size(arr)
        if(i==1) then
            sep=''
        else
            sep=','
        endif
        str=trim(adjustl(str))//sep//trim(adjustl(quote(tostr(arr(i))))) !! quote string
    enddo
    end function

    integer function list_size(text) result(npar)
    character(len=*),intent(in):: text
    character:: sep=','
    integer i
    npar=1
    do i=1,len(text)
        if(text(i:i)==sep) then
            npar=npar+1
        endif
    enddo
    end function
    integer function list_size_str(text) result(npar)
    character(len=*),intent(in):: text
    character:: sep=','
    logical:: in_sq=.false.,in_bq=.false.
    integer i
    npar=1
    do i=1,len(text)
        if(text(i:i)=="'") in_sq=.not.in_sq
        if(text(i:i)=='"') in_bq=.not.in_bq
        if(.not.in_sq .and. .not.in_bq .and. text(i:i)==sep) then
            npar=npar+1
        endif
    enddo
    end function
    integer function list_size_cmplx(text) result(npar)
    character(len=*),intent(in):: text
    character:: sep=','
    integer i
    logical:: inside=.false.
    npar=1
    do i=1,len(text)
        if(text(i:i)=='(') inside=.true.
        if(text(i:i)==')') inside=.false.
        if(.not.inside .and. text(i:i)==sep) then
            npar=npar+1
        endif
    enddo
    end function

end module
