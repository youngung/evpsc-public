CFEB
C *******************************************************************
C     SUBROUTINE VAR_VEL_GRAD_evp
C
C     IMPOSES A NON-UNIFORM DEFORMATION HISTORY TO THE AGGREGATE.
C     THE VELOCITY GRADIENT IS FULLY IMPOSED AND READ FROM 'FILEHIST'.
C *******************************************************************
      SUBROUTINE VAR_VEL_GRAD_evp(IOPTION,dtime,
     $     DBAR,dsim,iDBARv,iSBARv,ILBAR,LIJBARc,
     $     ictrl,strain_control,nsteps)
      implicit none
      integer, intent(in)::ioption
      real*8, intent(out):: dtime,DBAR(6),dsim(3,3),
     $     LIJBARc(3,3)
      integer, intent(out):: iDBARv(6),nsteps,ILBAR(3,3),ictrl
c     locals
      real*8 AUX3333(3,3,3,3),aux66(6,6),scauchy(3,3)
      integer i,j,idummy,UR_HIS,iSBARv(6),strain_control
      real*8 ctrlincr,temperat
      logical iplanestress
      parameter(iplanestress=.false.)
c      parameter(iplanestress=.true.)

      UR_HIS=99
C ***********************************************************************
      IF(IOPTION.EQ.0) THEN
         write(*,*)'UR_HIS:',UR_HIS
         READ(UR_HIS,*) NSTEPS,ICTRL,CTRLINCR,TEMPERAT
         READ(UR_HIS,*)
         CTRLINCR=CTRLINCR      ! to fool compiler
         IF(ICTRL.NE.7) THEN
            WRITE(*,*)'VARIABLE VELOCITY GRADIENT OPTION ASSUMES THAT'
            WRITE(*,*)'THE VELOCITY GRADIENT AND A TIME INCREMENT ARE'
            WRITE(*,*)'IMPOSED AT EACH STEP --> SET ICTRL=7 !!'
            call xit()
         ENDIF

C *** GENERAL SETTINGS FOR FULLY IMPOSED VARIABLE VELOCITY GRADIENT
         DO I=1,3
         DO J=1,3
            ILBAR(I,J)=1
            SCAUCHY(I,J)=0d0
         ENDDO
         ENDDO
         DO I=1,6
            iDBARv(I)=1
            iSBARv(I)=0
         ENDDO


      elseIF(IOPTION.EQ.1) THEN
C *** READS COMPONENTS OF VELOCITY GRADIENT AND TIME INCRMENT FOR THE STEP
         READ(UR_HIS,*) IDUMMY,((LIJBARc(I,J),J=1,3),I=1,3),DTIME
         ILBAR(:,:)=1
         iSBARv(:)=0
         IDUMMY=IDUMMY          ! to fool the compiler
         DO I=1,3
         DO J=1,3
            DSIM(I,J)=(LIJBARc(I,J)+LIJBARc(J,I))/2d0
         ENDDO
         ENDDO
         CALL CHG_BASIS (DBAR,DSIM,AUX66,AUX3333,2,6)
      ENDIF
C ***********************************************************************
      STRAIN_CONTROL=1


      if (iplanestress) then
c        For to relaxe a normal stress component.
         ILBAR(3,3)=0
         iDBARv(3)=0
         iSBARv(3)=1
      endif

      RETURN
      END subroutine var_vel_grad_evp
C
