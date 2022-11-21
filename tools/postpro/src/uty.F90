module uty
 use parameters

 contains

 subroutine invmat(mat,n)
 !
  implicit none
 !
  integer, intent(in) :: n
  real(rkind), dimension(n,n), intent(inout) :: mat
  integer :: i,j
  integer, dimension(n) :: indx
  real(rkind), dimension(n,n) :: y
  real(rkind) :: d
 !
  y = 0._rkind
  do i=1,n
   y(i,i) = 1._rkind
  enddo
  call ludcmp(mat,n,n,indx,d)
  do j=1,n
   call lubksb(mat,n,n,indx,y(1,j))
  enddo
  mat = y
 !
 end subroutine invmat
!
      SUBROUTINE LUDCMP(A,N,NP,INDX,D)
      implicit none
      integer, parameter :: NMAX = 100
      real(rkind), parameter :: TINY = 1.0D-20
      integer, intent(in) :: N,NP
      integer, dimension(N), intent(out) :: INDX
      real(rkind), intent(out) :: D
      real(rkind), dimension(NP,NP), intent(inout) :: A
      real(rkind), dimension(NMAX) :: VV
      integer :: I,J,K,IMAX
      real(rkind) :: SUM,DUM,AAMAX
!
      D=1._rkind
      DO 12 I=1,N
        AAMAX=0._rkind
        DO 11 J=1,N
          IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
11      CONTINUE
        IF (AAMAX.EQ.0._rkind) STOP 'Singular matrix.'
        VV(I)=1./AAMAX
12    CONTINUE
      DO 19 J=1,N
        IF (J.GT.1) THEN
          DO 14 I=1,J-1
            SUM=A(I,J)
            IF (I.GT.1)THEN
              DO 13 K=1,I-1
                SUM=SUM-A(I,K)*A(K,J)
13            CONTINUE
              A(I,J)=SUM
            ENDIF
14        CONTINUE
        ENDIF
       AAMAX=0.
        DO 16 I=J,N
          SUM=A(I,J)
          IF (J.GT.1)THEN
            DO 15 K=1,J-1
              SUM=SUM-A(I,K)*A(K,J)
15          CONTINUE
            A(I,J)=SUM
          ENDIF
          DUM=VV(I)*ABS(SUM)
          IF (DUM.GE.AAMAX) THEN
            IMAX=I
            AAMAX=DUM
          ENDIF
16      CONTINUE
        IF (J.NE.IMAX)THEN
          DO 17 K=1,N
            DUM=A(IMAX,K)
            A(IMAX,K)=A(J,K)
            A(J,K)=DUM
17        CONTINUE
          D=-D
          VV(IMAX)=VV(J)
        ENDIF
        INDX(J)=IMAX
        IF(J.NE.N)THEN
          IF(A(J,J).EQ.0._rkind)A(J,J)=TINY
          DUM=1./A(J,J)
          DO 18 I=J+1,N
            A(I,J)=A(I,J)*DUM
18        CONTINUE
        ENDIF
19    CONTINUE
      IF(A(N,N).EQ.0._rkind)A(N,N)=TINY
      RETURN
      END
!
      SUBROUTINE LUBKSB(A,N,NP,INDX,B)
      implicit none
      integer, intent(in) :: N,NP
      integer, dimension(N), intent(in) :: INDX
      real(rkind), dimension(NP,NP), intent(in) :: A
      real(rkind), dimension(N), intent(inout) :: B
      integer :: II,I,LL,J
      real(rkind) :: SUM
      II=0
      DO 12 I=1,N
        LL=INDX(I)
        SUM=B(LL)
        B(LL)=B(I)
        IF (II.NE.0)THEN
          DO 11 J=II,I-1
            SUM=SUM-A(I,J)*B(J)
11        CONTINUE
        ELSE IF (SUM.NE.0.) THEN
          II=I
        ENDIF
        B(I)=SUM
12    CONTINUE
      DO 14 I=N,1,-1
        SUM=B(I)
        IF(I.LT.N)THEN
          DO 13 J=I+1,N
            SUM=SUM-A(I,J)*B(J)
13        CONTINUE
        ENDIF
        B(I)=SUM/A(I,I)
14    CONTINUE
      RETURN
      END
!
end module
