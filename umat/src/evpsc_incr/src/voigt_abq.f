c     stress, stiffness
      SUBROUTINE VOIGT_ABQ(T1,T2,C2,C4,IOPT)

      implicit none
      real*8 T1(6),T2(3,3),C2(6,6),C4(3,3,3,3)
      integer IJV(6,2),i,j,iopt,i1,i2,j1,j2,n,m
      DATA ((IJV(N,M),M=1,2),N=1,6)/1,1,2,2,3,3,1,2,1,3,2,3/
c     (11),(22),(33),(12),(13),(23) --> this is for UMAT, for full space.

      IF(IOPT.EQ.1) THEN
c     Convert the given voigt-notated tensor components in 6D
c     to Cartesian components in 3x3 matrix.
      DO 30 I=1,6
      I1=IJV(I,1)
      I2=IJV(I,2)
      T2(I1,I2)=T1(I)
   30 T2(I2,I1)=T1(I)
      ENDIF
C
      IF(IOPT.EQ.2) THEN
      DO 40 I=1,6
      I1=IJV(I,1)
      I2=IJV(I,2)
   40 T1(I)=T2(I1,I2)
      ENDIF
C
      IF (IOPT.EQ.3) THEN
      DO 10 I=1,6
      I1=IJV(I,1)
      I2=IJV(I,2)
      DO 10 J=1,6
      J1=IJV(J,1)
      J2=IJV(J,2)
      C4(I1,I2,J1,J2)=C2(I,J)
      C4(I2,I1,J1,J2)=C2(I,J)
      C4(I1,I2,J2,J1)=C2(I,J)
   10 C4(I2,I1,J2,J1)=C2(I,J)
      ENDIF
C
      IF(IOPT.EQ.4) THEN
      DO 20 I=1,6
      I1=IJV(I,1)
      I2=IJV(I,2)
      DO 20 J=1,6
      J1=IJV(J,1)
      J2=IJV(J,2)
   20 C2(I,J)=C4(I1,I2,J1,J2)
      ENDIF
C
      RETURN
      END


      SUBROUTINE VOIGT_ABQ_2_basis(aux6_voigt_abq,aux6_bb)
      implicit none
      real*8 aux6_voigt_abq(6)
      real*8 aux6_bb(6)

c     locals
      real*8 aux6_voigt_evp(6) !transit
      real*8 aux33(3,3),aux66(6,6),aux3333(3,3,3,3)

      ! convert aux6_vo to bux6_voigt_evp
      aux6_voigt_evp(1)=aux6_voigt_abq(1)
      aux6_voigt_evp(2)=aux6_voigt_abq(2)
      aux6_voigt_evp(3)=aux6_voigt_abq(3)
      aux6_voigt_evp(4)=aux6_voigt_abq(6)
      aux6_voigt_evp(5)=aux6_voigt_abq(5)
      aux6_voigt_evp(6)=aux6_voigt_abq(4)

      ! convert bux6_voigt_evp to cartesian
      call voigt(aux6_voigt_evp,aux33,aux66,aux3333,1)

      ! convert aux33 to b-basis
      call chg_basis(aux6_bb,aux33,aux66,aux3333,2,6)

      return
      end SUBROUTINE VOIGT_ABQ_2_basis


C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     SUBROUTINE VOIGT_ABQ   ---->   VERSION OF
C
C     TRANSFORMS 6X1 MATRIX T1 INTO SECOND ORDER TENSOR T2 IF IOPT=1
C     AND VICEVERSA IF IOPT=2.
C     TRANSFORMS 6X6 MATRIX C2 INTO FOURTH ORDER TENSOR C4 IF IOPT=3
C     AND VICEVERSA IF IOPT=4.
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c     strain, compliance
      SUBROUTINE VOIGT_ABQ_strain(T1,T2,C2,C4,IOPT)
      implicit none
      real*8 T1(6),T2(3,3),C2(6,6),C4(3,3,3,3),fact(6,6)
      integer IJV(6,2),n,m,i,j,i1,i2,iopt,j1,j2
      DATA ((IJV(N,M),M=1,2),N=1,6)/1,1,2,2,3,3,1,2,1,3,2,3/
      DATA ((FACT(n,m),n=1,6),m=1,6)/1d0,1d0,1d0,2d0,2d0,2d0,1d0,1d0,
     $     1d0,2d0,2d0,2d0,1d0,1d0,1d0,2d0,2d0,2d0,2d0,2d0,2d0,4d0,4d0,
     $     4d0,2d0,2d0,2d0,4d0,4d0,4d0,2d0,2d0,2d0,4d0,4d0,4d0/
      IF(IOPT.EQ.1) THEN
         DO 30 I=1,3
            I1=IJV(I,1)
            I2=IJV(I,2)
            T2(I1,I2)=T1(I)
            T2(I2,I1)=T1(I)
 30      CONTINUE
         DO 31 I=4,6
            I1=IJV(I,1)
            I2=IJV(I,2)
            T2(I1,I2)=T1(I)/2d0
            T2(I2,I1)=T1(I)/2d0
 31      CONTINUE
      ENDIF
C
      IF(IOPT.EQ.2) THEN
         DO 40 I=1,3
            I1=IJV(I,1)
            I2=IJV(I,2)
            T1(I)=T2(I1,I2)
 40      CONTINUE
         DO 41 I=4,6
            I1=IJV(I,1)
            I2=IJV(I,2)
            T1(I)=T2(I1,I2)*2d0
 41      CONTINUE
      ENDIF
C
      IF (IOPT.EQ.3) THEN
         DO 10 I=1,6
            I1=IJV(I,1)
            I2=IJV(I,2)
      DO 10 J=1,6
         J1=IJV(J,1)
         J2=IJV(J,2)
            C4(I1,I2,J1,J2)=C2(I,J)*(1d0/fact(i,j))
c            write(*,*)'Cijkl=cij',i1,i2,j1,j2
c            write(*,*)'cij',i,j
c            write(*,*)'fac',(1/fact(i,j))
            C4(I2,I1,J1,J2)=C2(I,J)*(1d0/fact(i,j))
            C4(I1,I2,J2,J1)=C2(I,J)*(1d0/fact(i,j))
            C4(I2,I1,J2,J1)=C2(I,J)*(1d0/fact(i,j))
 10   CONTINUE
      ENDIF
C
      IF(IOPT.EQ.4) THEN
         DO 20 I=1,6
            I1=IJV(I,1)
            I2=IJV(I,2)
         DO 20 J=1,6
            J1=IJV(J,1)
            J2=IJV(J,2)
            C2(I,J)=fact(i,j)*C4(I1,I2,J1,J2)
c$$$            write(*,*)'cij',i,j
c$$$            write(*,*)'Cijkl=cij',i1,i2,j1,j2
c$$$            write(*,*)'fac',fact(i,j)
 20      CONTINUE
      ENDIF

      RETURN
      END
