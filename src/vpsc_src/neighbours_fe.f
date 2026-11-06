C ***************************************************************************
C     SUBROUTINE NEIGHBOURS_fe      --->      VERSION 13/Mar/2020
C
C *** IDENTIFIES 'NNEIGH' NEIGHBORS FOR EVERY GRAIN AND ASSIGNS WEIGHTS
C     'WNEIGH' WHICH ARE GOING TO BE USED EXCLUSIVELY FOR CO-ROTATION
C     PURPOSES INSIDE SUBROUTINE UPDATE_ORIENTATION.
C *** IF IOPTION=0 DEFINES NEIGHBORS FOR EACH GRAIN PICKED AT RANDOM FROM
C     ORIENTATION FILE (REORIENTATION WORKS ALSO FOR MULTIELEMENT RUN)
C *** IF IOPTION=1 DEFINES NEIGHBORS FOR A TWO PHASE MATERIAL WHERE THE
C     SECOND PHASE IS CRYSTALLOGRAPHICALLY ASSOCIATED TO THE FIRST PHASE
C     (i.e.: MATRIX/TWIN or PARENT/CHILD IN RECRYSTALLIZATION)
C ***************************************************************************

      SUBROUTINE NEIGHBOURS_fe(
c          in
     $     IOPTION,jran,nneigh,nneimx,ngrmx,ngr,wgt,
c          out
     $     neigh,wneigh
     $     )
      implicit none
      integer,intent(in)::ioption,jran,nneigh,nneimx,ngrmx,ngr
      integer,intent(out)::neigh(0:nneimx,ngrmx)
      real*8,intent(in)::wgt(ngrmx)
      real*8,intent(out)::wneigh(0:nneimx,ngrmx)
c     ------------------------
      real*8 WAUX(0:nneigh)
      integer i,kgx
      real*8 wtot,rand,ran2
      logical verbose
      parameter(verbose=.false.)
c      parameter(verbose=.true.)

      if (verbose) then
         write(*,*)'nneigh:',nneigh
         write(*,*)'ngr:',ngr
      endif

      DO I=0,NNEIGH
         WAUX(I)=1d0/(NNEIGH+1)  ! DEFINES EQUAL CO-ROTATION WGT's
      ENDDO

      DO KGX=1,NGR              ! SCANS ALL GRAINS IN ALL PHASES
         NEIGH(0,KGX)=KGX       ! NEIGHBOR ZERO IS THE GRAIN ITSELF
         WTOT=WGT(KGX)*WAUX(0)
         IF(NNEIGH.GT.0) THEN
            IF(IOPTION.EQ.0) THEN ! PICK NEIGHBOURS AT RANDOM
               DO I=1,NNEIGH
                  RAND=RAN2(JRAN)
                  NEIGH(I,KGX)=INT(RAND*NGR+1)
                  WTOT=WTOT+WGT(NEIGH(I,KGX))*WAUX(I)
               ENDDO
            ELSE IF(IOPTION.EQ.1) THEN ! PICK NEIGHBOURS IN ORDER
               IF (KGX.LE.NGR) NEIGH(1,KGX)=KGX+NGR
               IF (KGX.GT.NGR) NEIGH(1,KGX)=KGX-NGR
               WTOT=WTOT+WGT(NEIGH(1,KGX))*WAUX(1)
            ENDIF
         ENDIF
C *** NORMALIZATION OF CO-ROTATION WGT's

         IF(WTOT.LT.1.d-10) THEN
            DO I=0,NNEIGH
               WNEIGH(I,KGX)=0d0
            ENDDO
         ELSE
            DO I=0,NNEIGH
               WNEIGH(I,KGX)=WAUX(I)*WGT(NEIGH(I,KGX))/WTOT
            ENDDO
         ENDIF
      ENDDO                     !  END LOOP OVER ORIENTATIONS


      if (verbose) then
         write(*,*)'weigh'
         do kgx=1,ngr
            do i=0,nneigh
               write(*,*)'neigh:',i
               write(*,*)wneigh(i,kgx)
            enddo
         enddo
         write(*,*)'wtot:',wtot
         call stop_debug(0,'debugging neighbours_fe')
      endif

      END subroutine neighbours_fe
