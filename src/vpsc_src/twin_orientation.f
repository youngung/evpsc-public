C
C *****************************************************************************
C     SUBROUTINE TWIN_ORIENTATION    --->    VERSION OF 08/DEC/2011
C
C     GIVEN THE SHEAR DIRECTION 'TWB' AND TWIN PLANE NORMAL 'TWN', DEFINES
C     THE ORIENTATION OF THE TWIN RELATED CRYSTAL.
C
C     FOR TYPE I twins: MAKES A ROTATION OF 180 deg AROUND NORMAL TO TWIN
C       PLANE DEFINED BY 'TWN'
C     FOR TYPE II twins: MAKES A ROTATION OF 180 deg AROUND SHEAR DIRECTION
C       DEFINED BY 'TWB'
C
C     USES THE RODRIGUES FORMULA FOR THE CASE phi=180 deg
C          --> R=I+sin(phi)*N+(1-cos(phi))*N^2=-I+2*V^2=-I+2*Vi*Vj
C
C     THE MATRIX 'ATWIN' TRANSFORMS FROM TWIN TO PARENT CRYSTAL AXES
C *****************************************************************************

      SUBROUTINE TWIN_ORIENTATION (ITWTYPE,TWB,TWN,ATWIN)
      implicit none

      real*8 TWB(3),TWN(3),ATWIN(3,3),XID33(3,3)
      integer i,j,itwtype

      DATA XID33/1.,0.,0.,0.,1.,0.,0.,0., 1./


C *** ROTATION OF 180 deg AROUND NORMAL TO TWIN PLANE
      IF(ITWTYPE.EQ.1) THEN
        DO I=1,3
        DO J=1,3
          ATWIN(I,J)=2.*TWN(I)*TWN(J)-XID33(I,J)
        ENDDO
        ENDDO
      ENDIF

C *** ROTATION OF 180 deg AROUND SHEAR VECTOR ON TWIN PLANE
      IF(ITWTYPE.EQ.2) THEN
        DO I=1,3
        DO J=1,3
          ATWIN(I,J)=2.*TWB(I)*TWB(J)-XID33(I,J)
        ENDDO
        ENDDO
      ENDIF

      RETURN
      END
