C
C ************************************************************************
C     SUBROUTINE CHG_BASIS    --->   VERSION 19/JUL/01
C
C     (modif. 06/FEB/98 - same convention as SELFPOLY - C.N.T.)
C     (modif. 16/JUN/99 - same convention as Maudlin  - C.N.T.)
C     (modif. 10/MAY/01 - KDIM version - R.L.)
C
C     KDIM=5 or 6, FOR DEVIATORIC or DEV+HYDROST TENSORS, RESPECTIVELY.
C     IOPT=0: DEFINES A BASIS OF 6 SECOND ORDER TENSORS B(N).
C     IOPT=1: CALCULATES SECOND ORDER TENSOR 'C2' AS AN EXPANSION IN TERMS
C             OF VECTOR COMPONENTS CE2(KDIM) AND THE BASIS TENSORS B(KDIM).
C     IOPT=2: CALCULATES COMPONENTS OF C2 AS A VECTOR CE2(KDIM).
C     IOPT=3: CALCULATES FOURTH ORDER TENSOR 'C4' AS AN EXPANSION IN TERMS
C             OF MATRIX COMPONENTS CE4(K,K) AND THE BASIS TENSORS B(KDIM).
C     IOPT=4: CALCULATES MATRIX COMPONENTS CE4(K,K) OF TENSOR 'C4'.
C **************************************************************************

      SUBROUTINE CHG_BASIS(CE2,C2,CE4,C4,IOPT,KDIM)
      implicit none
      integer, intent(in):: iopt,kdim
      real*8 ce2(kdim),c2(3,3),ce4(kdim,kdim),c4(3,3,3,3)
c      DIMENSION B(3,3,6)
      real*8 B(3,3,6)
c     locals
      integer i,j,k,l,m,n

C *** INITIALIZES b-BASIS TENSORS B(N)
      B(1:3,1:3,1:6)=0.0

      B(1,1,2)=-1./SQRT(6.)
      B(2,2,2)=-1./SQRT(6.)
      B(3,3,2)= 2./SQRT(6.)

      B(1,1,1)=-1./SQRT(2.)
      B(2,2,1)= 1./SQRT(2.)

      B(2,3,3)= 1./SQRT(2.)
      B(3,2,3)= 1./SQRT(2.)

      B(1,3,4)= 1./SQRT(2.)
      B(3,1,4)= 1./SQRT(2.)

      B(1,2,5)= 1./SQRT(2.)
      B(2,1,5)= 1./SQRT(2.)

      B(1,1,6)= 1./SQRT(3.)
      B(2,2,6)= 1./SQRT(3.)
      B(3,3,6)= 1./SQRT(3.)


C *** CALCULATES CARTESIAN SECOND ORDER TENSOR FROM b-COMPONENTS VECTOR.
      IF(IOPT.EQ.1) THEN
        DO I=1,3
        DO J=1,3
           C2(I,J)=0d0
        DO N=1,KDIM
           C2(I,J)=C2(I,J)+CE2(N)*B(I,J,N)
        enddo
        enddo
        enddo
      ENDIF

C *** CALCULATES KDIMx1 b-COMPONENTS VECTOR FROM SECOND ORDER TENSOR.
      IF(IOPT.EQ.2) THEN
         CE2(:)=0.0
         DO N=1,KDIM
         DO I=1,3
         DO J=1,3
            CE2(N)=CE2(N)+C2(I,J)*B(I,J,N)
         enddo
         enddo
         enddo
      ENDIF

C *** CALCULATES FOURTH ORDER TENSOR FROM b-COMPONENTS MATRIX.
      IF(IOPT.EQ.3) THEN
         c4(:,:,:,:)=0d0
         DO I=1,3
         DO J=1,3
         DO K=1,3
         DO L=1,3
         DO N=1,KDIM
         DO M=1,KDIM
            C4(I,J,K,L)=C4(I,J,K,L)+CE4(N,M)*B(I,J,N)*B(K,L,M)
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo

      ENDIF

C *** CALCULATES KDIMxKDIM b-COMPONENTS MATRIX FROM FOURTH ORDER TENSOR.
      IF(IOPT.EQ.4) THEN
         CE4(:,:)=0d0
        DO N=1,KDIM
        DO M=1,KDIM

        DO I=1,3
        DO J=1,3
        DO K=1,3
        DO L=1,3
           CE4(N,M)=CE4(N,M)+C4(I,J,K,L)*B(I,J,N)*B(K,L,M)
        enddo
        enddo
        enddo
        enddo
        enddo
        enddo
      ENDIF

      RETURN
      END subroutine chg_basis
