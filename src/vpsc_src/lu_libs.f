      SUBROUTINE LU_INVERSE (ierr,A,N)
C *** INVERTS A MATRIX USING LU DECOMPOSITION
      implicit none
      integer, intent(in)::n
c      real*8, intent(inout):: a(n,n)
      real*8 a(n,n)
      integer, intent(out)::ierr
c      DIMENSION A(N,N),Y(N,N),INDX(N)      ! MAY CHOKE SOME COMPILERS

      real*8 y(n,n),d,dum,amax,z(n,n)
      integer indx(n)
      integer i,j,isingular

      ierr=0

c      write(*,*) 'A(i,j) matrix inside lu_inverse'
c      write(*,'(6e12.3)') ((a(i,j),j=1,n),i=1,n)
c     read(*,*)

C **************************************************************
C *** BLOCK ADDED 03/DEC/05 TO AVOID NUMERICALLY SINGULAR MATRIX
      AMAX=0D0
      DO I=1,N
      DO J=1,N
         DUM=ABS(A(I,J))
         z(i,j)=a(i,j)
         IF(DUM .GT. AMAX) AMAX=DUM
      ENDDO
      ENDDO

      if(amax.eq.0d0)then
         write(*,*)'** amax of zero found in sub. LU_INVERSE!'
         write(*,*)'matrix a before passing to ludcmp'
         do i=1,n
            write(*,'(20e20.12)')(a(i,j),j=1,n)
         enddo
         write(*,*)'matrix a before normalization'
         do i=1,n
            write(*,'(20e20.12)')(z(i,j),j=1,n)
         enddo
         ierr=1
c         call xit()
         return
      endif

      DO I=1,N
      DO J=1,N
         A(I,J)=A(I,J)/AMAX     ! normalize the matrix
      ENDDO
      ENDDO
C **************************************************************

      DO I=1,N
         DO J=1,N
            Y(I,J)=0d0
         ENDDO
         Y(I,I)=1d0
      ENDDO

      CALL LUDCMP(A,N,N,INDX,D,ISINGULAR)
      IF(ISINGULAR.EQ.1) THEN
         ierr=2
         write(*,*)' *** SINGULAR MATRIX IN LU_INVERSE !!'
         write(*,*)'amax:',amax
         write(*,*)'matrix a'
         do i=1,n
            write(*,'(20e13.5)')(a(i,j),j=1,n)
         enddo
         write(*,*)'matrix a (before normalization) (z)'
         do i=1,n
            write(*,'(20e20.12)')(z(i,j),j=1,n)
         enddo
         write(*,*)'---'
         return
c         call xit()
c         call w_chr(7,' *** SINGULAR MATRIX IN LU_INVERSE !!')
c         call w_mdim(7,a,n,1d0)
c         call stop_debug(0,'Stopped in subr. lu_inverse')
      ENDIF

      DO J=1,N
        CALL LUBKSB(A,N,N,INDX,Y(1,J))
      ENDDO

      DO I=1,N
      DO J=1,N
        A(I,J)=Y(I,J) /AMAX      ! renormalize the inverse
      ENDDO
      ENDDO

      RETURN
      END

C *************************************************************************

      SUBROUTINE LU_EQSYSTEM(A,B,N,ISINGULAR)
      implicit none
C *** SOLVES A*X=B USING LU DECOMPOSITION
      integer,intent(in):: n
      integer,intent(out):: isingular
      real*8 A(N,N),B(N)
      real*8 d
      integer INDX(N)                  ! MAY CHOKE SOME COMPILERS

C     DIMENSION A(5,5),B(5),INDX(5)

      CALL LUDCMP(A,N,N,INDX,D,ISINGULAR)

      IF(ISINGULAR.EQ.1) RETURN

      CALL LUBKSB(A,N,N,INDX,B)

      RETURN
      END
C
C ***********************************************************
      SUBROUTINE LU_EQSYS25(A,B,N,ISINGULAR)
      implicit none

C *** SOLVES A*X=B USING LU DECOMPOSITION

c      DIMENSION A(N,N),B(N),INDX(N)      ! MAY CHOKE SOME COMPILERS
      real*8 A(25,25),B(25)
      integer INDX(25),j,isingular,n
      real*8 d

      CALL LUDCMP(A,N,N,INDX,D,ISINGULAR)
c
      do j=1,25
        d=d*a(j,j)
      enddo

      IF(ISINGULAR.EQ.1) RETURN

      CALL LUBKSB(A,N,N,INDX,B)

      RETURN
      END
