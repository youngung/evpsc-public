C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     SUBROUTINE VOIGT   ---->   VERSION OF 09/02/98
C
C     TRANSFORMS 6X1 MATRIX T1 INTO SECOND ORDER TENSOR T2 IF IOPT=1
C     AND VICEVERSA IF IOPT=2.
C     TRANSFORMS 6X6 MATRIX C2 INTO FOURTH ORDER TENSOR C4 IF IOPT=3
C     AND VICEVERSA IF IOPT=4.
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE VOIGT(T1,T2,C2,C4,IOPT)
      implicit none
      real*8 T1(6),T2(3,3),C2(6,6),C4(3,3,3,3)
      integer iopt,IJV(6,2),i,j,i1,i2,j1,j2

      ijv(1,1)=1
      ijv(1,2)=1
      ijv(2,1)=2
      ijv(2,2)=2
      ijv(3,1)=3
      ijv(3,2)=3
      ijv(4,1)=2
      ijv(4,2)=3
      ijv(5,1)=1
      ijv(5,2)=3
      ijv(6,1)=1
      ijv(6,2)=2

      IF(IOPT.EQ.1) THEN
         t2(:,:)=0d0
         DO I=1,6
            I1=IJV(I,1)
            I2=IJV(I,2)
            T2(I1,I2)=T1(I)
            T2(I2,I1)=T1(I)
         enddo
C
      elseIF(IOPT.EQ.2) THEN
         t1(:)=0d0
         DO I=1,6
            I1=IJV(I,1)
            I2=IJV(I,2)
            T1(I)=T2(I1,I2)
         enddo
C
      elseIF (IOPT.EQ.3) THEN
         c4(:,:,:,:)=0d0
         DO I=1,6
            I1=IJV(I,1)
            I2=IJV(I,2)
         DO J=1,6
            J1=IJV(J,1)
            J2=IJV(J,2)
            C4(I1,I2,J1,J2)=C2(I,J)
            C4(I2,I1,J1,J2)=C2(I,J)
            C4(I1,I2,J2,J1)=C2(I,J)
            C4(I2,I1,J2,J1)=C2(I,J)
         enddo
         enddo
C
      elseIF(IOPT.EQ.4) THEN
         c2(:,:)=0d0
         DO I=1,6
            I1=IJV(I,1)
            I2=IJV(I,2)
         DO J=1,6
            J1=IJV(J,1)
            J2=IJV(J,2)
            C2(I,J)=C4(I1,I2,J1,J2)
         enddo
         enddo
      else
         write(*,*) 'unexpected iopt value found in voigt'
         call xit()
      ENDIF
C
      RETURN
      END subroutine voigt
C ************************************************************
c     convert between voigt and b-basis notations
C ************************************************************
      subroutine voigt_bbasis(v6,b6,v66,b66,iopt)
c     iopt=1 v6->b6
c     iopt=2 v6<-b6
c     iopt=3 v66->b66
c     iopt=4 v66<-b66
      implicit none
      integer iopt
      real*8 v6(6),b6(6),v66(6,6),b66(6,6)
c     local
      real*8 aux6(6),aux33(3,3),aux66(6,6),aux3333(3,3,3,3)


      if (iopt.eq.1) then
         call voigt(v6,aux33,aux66,aux3333,1)
         call chg_basis(b6,aux33,aux66,aux3333,2,6)
      elseif (iopt.eq.2) then
         call chg_basis(b6,aux33,aux66,aux3333,1,6)
         call voigt(v6,aux33,aux66,aux3333,2)
      elseif (iopt.eq.3) then
         call voigt(aux6,aux33,v66,aux3333,3)
         call chg_basis(aux6,aux33,b66,aux3333,4,6)
      elseif (iopt.eq.4) then
         call chg_basis(aux6,aux33,b66,aux3333,3,6)
         call voigt(aux6,aux33,v66,aux3333,4)
      else
         write(*,*)'unexpected iop given to voigt_bbasis'
         call xit()
      endif


      return
      end subroutine voigt_bbasis

C ************************************************************
C     SUBROUTINE VOIGT10     --->     VERSION 4/JAN/07
C ************************************************************

      SUBROUTINE VOIGT10(T1,T2,IOPT)
      implicit none
      real*8 T1(10),T2(4,4)
      integer IJV(10,2),n,m
      integer,intent(in):: iopt
      DATA ((IJV(N,M),M=1,2),N=1,10)
     #  /1,1,2,2,3,3,2,3,1,3,1,2,1,4,2,4,3,4,4,4/
      integer i,i1,i2

      IF(IOPT.EQ.1) THEN
      DO I=1,10
      I1=IJV(I,1)
      I2=IJV(I,2)
      T2(I1,I2)=T1(I)
      T2(I2,I1)=T1(I)
      enddo
      ENDIF
C
      IF(IOPT.EQ.2) THEN
      DO I=1,10
      I1=IJV(I,1)
      I2=IJV(I,2)
      T1(I)=T2(I1,I2)
      enddo
      ENDIF
C
      RETURN
      END
CFEE
