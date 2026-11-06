C **********************************************************************
C     SUBROUTINE ESHELBY_ROT      --->      VERSION 14/APR/2018

C *** ESHELBY CALCULATION FOR EVERY PHASE OR FOR EVERY GRAIN
C *** ROTATES STIFFNESS TO ELLIPSOID AXES 'EIGB'
C *** CALLS ESHELBY TO CALCULATE DISTORTION AND ROTATION ESHELBY TENSORS
C     'E4GA' & 'R4GA'IN ELLIPSOID AXES
C *** ROTATES THEM BACK TO SAMPLE AXES 'E4SA' & 'R4SA'
C **********************************************************************

      SUBROUTINE ESHELBY_ROT(AXB,EIGB,C4SA,E4SA,R4SA,IOPTION)
      implicit none
      integer, intent(in)::ioption
      real*8, intent(in)::axb(3),eigb(3,3)

      real*8 AUX33(3,3),AUX3333(3,3,3,3)
      real*8 C4SA(3,3,3,3),C4GA(3,3,3,3)
      real*8 E4SA(3,3,3,3),E4GA(3,3,3,3),R4SA(3,3,3,3),R4GA(3,3,3,3)
c      DIMENSION AXB(3),EIGB(3,3)

c      COMMON/FLUCT2/C4GA_FLU(3,3,3,3)      ! CT APR 13 2018
c      real*8 c4ga_flu
      real*8 dummy,dummy1,dummy2


c     locals
      integer i,j,m,n,i1,j1,m1,n1
      real*8 pdil
      logical iverb,iesh_fail
      parameter(iverb=.false.)

      if (iverb) then
         write(*,*)'eigb in eshelby_rot'
         do i=1,3
            write(*,'(3f6.2)')(eigb(i,j),j=1,3)
         enddo
         write(*,*)'axb'
         write(*,'(3e11.3)')(axb(i),i=1,3)
c      if (iverb) then
         write(*,*)'c4sa(1,1,1,1):',c4sa(1,1,1,1)
         write(*,*)'c4sa(3,3,3,3):',c4sa(3,3,3,3)
      endif

c      call xit()


C *** ROTATION OF STIFFNESS 'C4' FROM SAMPLE TO ELLIPSOID AXES
      DO I=1,3
      DO J=1,3
      DO M=1,3
      DO N=1,3
      DUMMY=0d0
      DO I1=1,3
      DO J1=1,3
      DO M1=1,3
      DO N1=1,3
         DUMMY=DUMMY+EIGB(I1,I)*EIGB(J1,J)*EIGB(M1,M)
     $        *EIGB(N1,N)*C4SA(I1,J1,M1,N1)
      enddo
      enddo
      enddo
      enddo
c      write(*,*)'dummy:',dummy
      C4GA(I,J,M,N)=DUMMY
c      C4GA_FLU(I,J,M,N)=DUMMY   ! to be used later by GET_THEFLU
      enddo
      enddo
      enddo
      enddo

c      call xit()

      if (iverb) then
         e4ga(:,:,:,:)=0d0
         write(*,*)'e4ga(1,1,1,1):',e4ga(1,1,1,1)
         write(*,*)'e4ga(3,3,3,3):',e4ga(3,3,3,3)
         write(*,*)'ioption:',ioption
      endif

      CALL ESHELBY(AXB,C4GA,0d0,E4GA,R4GA,AUX33,AUX33,
     $     PDIL,AUX3333,AUX3333,IOPTION,iesh_fail)

      if (iesh_fail) then
         write(*,*)'** eshelby subroutine failed in eshelby_rot'
         call xit()
      endif


      if (iverb) then
         write(*,*)'e4ga(1,1,1,1):',e4ga(1,1,1,1)
         write(*,*)'e4ga(3,3,3,3):',e4ga(3,3,3,3)
      endif
c     call xit()



C *** ROTATES THE ESHELBY TENSOR FOR PHASE OR GRAIN BACK INTO SAMPLE AXES.
      DO I=1,3
      DO J=1,3
      DO M=1,3
      DO N=1,3
      DUMMY1=0d0
      DUMMY2=0d0
      DO I1=1,3
      DO J1=1,3
      DO M1=1,3
      DO N1=1,3
         DUMMY1=DUMMY1+EIGB(I,I1)*EIGB(J,J1)*EIGB(M,M1)
     $        *EIGB(N,N1)*E4GA(I1,J1,M1,N1)
         DUMMY2=DUMMY2+EIGB(I,I1)*EIGB(J,J1)*EIGB(M,M1)
     $        *EIGB(N,N1)*R4GA(I1,J1,M1,N1)
      enddo
      enddo
      enddo
      enddo

c     write(*,*)'dummy1:',dummy1
      E4SA(I,J,M,N)=DUMMY1
      R4SA(I,J,M,N)=DUMMY2
      enddo
      enddo
      enddo
      enddo

c     call xit()
      RETURN
      END
