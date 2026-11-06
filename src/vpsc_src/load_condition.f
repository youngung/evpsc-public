C ****************************************************************************
C     SUBROUTINE LOAD_CONDITIONS      -->      VERSION 22/MAY/2026
C
C     READS BOUNDARY CONDITIONS ON STRAIN RATE AND STRESS:
C        * NUMBER OF DEF STEPS, CONTROL VARIABLE, INCREMENT, TEMPERATURE
C        * IMPOSED AND RELAXED COMPONENTS OF VELOCITY GRADIENT AND STRESS
C     CALCULATES SYMMETRIC STRAIN RATE 'DBARc'
C     CALCULATES 5(or 6)-DIM VECTORS 'DBAR' AND 'SBAR'.
C     CHECKS WHETHER THE BOUNDARY CONDITIONS ARE CONSISTENT.
C ****************************************************************************

      SUBROUTINE LOAD_CONDITIONS(
c          in
     $     iregime,UR_HIS,interaction,
c          out
     $     ictrl,nsteps,iDBARv,iSBARv,DBAR,DBARc,iLbar,LijBARc,
     $     dsbar,sbar,sbarc,temp_ini,tempdot,eijincr,DBARctrl,
     $     strain_control,dtime
     $     )
      implicit none
      INTEGER, intent(in)::iregime,UR_HIS,interaction
      real*8,intent(out)::DBAR(6),DBARc(3,3),LijBARc(3,3),dsbar(6),
     $     sbar(6),sbarc(3,3),temp_ini,tempdot,eijincr,DBARctrl,dtime
      integer, intent(out)::nsteps,iDBARv(6),iSBARv(6),iLbar(3,3),ictrl,
     $     strain_control
c     --------------
c     locals
      real*8 SBARv(6),DBARv(6),evmincr,ctrlincr,vnorm,DBARNORM,SBARNORM,
     $     temp_fin,AUX5(5),AUX55(5,5),
     $     aux66(6,6),AUX3333(3,3,3,3),PROFAC(6)
      integer IJV(6,2),i,m,n,j,iSBARsum,iDBARsum,mdim
      DATA ((IJV(N,M),M=1,2),N=1,6)/1,1,2,2,3,3,2,3,1,3,1,2/
      character*125 prosa

      mdim=6
      if (iregime.eq.1) mdim=5

      profac(1)=1.
      profac(2)=1.
      profac(3)=1.
      profac(4)=2.
      profac(5)=2.
      profac(6)=2.

C *** WRITES LOAD CONDITIONS FILE INTO 'RUN_LOG.OUT' FILE
      DO I=1,100
        READ(UNIT=UR_HIS,END=100,FMT='(A)') PROSA
      ENDDO
  100 REWIND UR_HIS

C *** READS BOUNDARY CONDITIONS ON DEFORMATION RATE AND STRESS COMPONENTS

      READ(UR_HIS,*,err=120) NSTEPS,ICTRL,CTRLINCR,Temp_Ini,temp_fin
      go to 130
 120  write(*,*)'history file format changed; you should specify'
      write(*,*)'both the initial and the final temperature'
      call stop_debug(0,'Consult with YJ: yjeong@changwon.ac.kr')
 130  IF(CTRLINCR.LE.0) THEN
         WRITE(*,'(''STRAIN OR TIME INCREMENT (CTRLINCR) HAS TO BE'')')
         WRITE(*,'(''POSITIVE  --> THE SIGN OF THE DEFORMATION IS'')')
         WRITE(*,'(''CONTROLLED BY THE SIGN OF THE COMP IMPOSED'')')
         call exit(-1)
      ENDIF

      READ(UR_HIS,*)
      DO I=1,3
         READ(UR_HIS,*) (ILBAR(I,J),J=1,3)
      ENDDO
      READ(UR_HIS,*)
      DO I=1,3
         READ(UR_HIS,*) (LIJBARC(I,J),J=1,3)
      ENDDO
      READ(UR_HIS,*)
      READ(UR_HIS,*) iSBARv(1),iSBARv(6),iSBARv(5)
      READ(UR_HIS,*)          iSBARv(2),iSBARv(4)
      READ(UR_HIS,*)                   iSBARv(3)
      READ(UR_HIS,*)
      READ(UR_HIS,*) SBARC(1,1),SBARC(1,2),SBARC(1,3)
      READ(UR_HIS,*)              SBARC(2,2),SBARC(2,3)
      READ(UR_HIS,*)                           SBARC(3,3)
      SBARC(3,2)=SBARC(2,3)
      SBARC(3,1)=SBARC(1,3)
      SBARC(2,1)=SBARC(1,2)

C **********************************************************
C *** CHECKS WHETHER THE BOUNDARY CONDITIONS ARE CONSISTENT
C **********************************************************

      iSBARsum=0
      do i=1,6
         iSBARsum=iSBARsum+iSBARv(i)
      enddo

      write(*,*)'interaction:',interaction

      IF(INTERACTION.EQ.0 .AND. ISBARSUM.NE.0) THEN
         WRITE(*,*) ' CANNOT IMPOSE STRESS COMPONENTS FOR FC CASE'
         call xit()
      ENDIF

      if (iregime.eq.1.and.ILBAR(1,1)+ILBAR(2,2)+ILBAR(3,3).eq.2) then
         write(*,*)'** CHECK DIAGONAL BOUNDARY CONDITIONS ILBAR'
         write(*,*)'For VPSC, ENFORCING TWO DIAGONAL COMPONENTS FIXES'
         write(*,*)'THE REMAINING COMPONENT BECAUSE OF THE VOLUME'
         write(*,*)'INCOMPRESSIBILITY'
         call xit()
      endif

c     ** Check if both of the off-diagonal components are relaxed.
      do i=1,2
      do j=i+1,3
         if(iLbar(i,j)+iLbar(j,i).eq.0) then
            write(*,*) 'CHECK OFF-DIAGONAL BOUNDARY CONDITIONS'
            write(*,*) 'ILBAR CANNOT RELAX BOTH OFF-DIAGONAL '
            write(*,*) 'VEL GRAD COMPS'
            call exit(-1)
         endif
      enddo
      enddo

C *** CALCULATES SYMMETRIC STRAIN-RATE COMPONENTS FROM VELOCITY GRADIENT.
C *** TAGS IMPOSED & UNKNOWN CARTESIAN COMPONENTS
      DO I=1,3
      DO J=1,3
         DBARC(I,J)=(LIJBARC(I,J)+LIJBARC(J,I))*0.5
      ENDDO
      ENDDO

      IDBARV(1)=ILBAR(1,1)
      IDBARV(2)=ILBAR(2,2)
      IDBARV(3)=ILBAR(3,3)
      IDBARV(4)=ILBAR(2,3)*ILBAR(3,2)
      IDBARV(5)=ILBAR(1,3)*ILBAR(3,1)
      IDBARV(6)=ILBAR(1,2)*ILBAR(2,1)
      iDBARsum=0
      do i=1,6
         iDBARsum=iDBARsum+iDBARv(i)
      enddo

      DO I=1,6
         IF(iSBARv(I)*IDBARV(I).NE.0.OR.(iSBARv(I)+IDBARV(I)).NE.1) THEN
            WRITE(*,*) ' BOUNDARY CONDITIONS ON STRAIN-RATE AND STRESS'
            WRITE(*,*) ' COMPONENTS MUST BE COMPLEMENTARY'
            WRITE(*,'('' IDBARV = '',6I3)') IDBARV
            WRITE(*,'('' iSBARv = '',6I3)') iSBARv
            call exit(-1)
         ENDIF
      ENDDO

C *** WRITES 6-DIM REPRESENT OF MACRO STRAIN RATE 'DBARC' AND STRESS 'SBARC'.
C     BOTH VECTORS ARE NOT REQUIRED SIMULTANEOUSLY. THEY ARE USED TO MAKE
C     AN INITIAL GUESS FOR THE FIRST STEP OF THE SC PROCEDURE.

C ***********************************************************************
C     IN CONNECTION WITH THE OVERALL STRAIN RATE AND STRESS:
C     *  IF ICTRL=0 WILL IMPOSE A VON MISES STRAIN INCREMENT. THIS OPTION
C        VALID ONLY FOR FULLY IMPOSED STRAIN TENSOR AND PLASTICITY.
C        WILL CALCULATE THE TIME INCREMENT NECESSARY TO ACHIEVE IT
C     *  IF ICTRL=1,6 WILL EITHER
C        --> IMPOSE A STRAIN RATE COMPONENT AND A STRAIN INCREMENT &
C            CALCULATE THE TIME INCREMENT REQUIRED TO ACHIEVE IT OR,
C        --> IMPOSE A STRESS COMPONENT (CREEP) AND A TIME INCREMENT.
C            THE STRAIN INCREMENT THAT RESULTS FOLLOWS FROM THE (E)VPSC
C            CALCULATED STRAIN RATE TIMES THE TIME INCREMENT.
C     *  THE OPTION ICTRL=7 MEANT TO READ AN IMPOSED TIME INCREMENT.
C ***********************************************************************

      CALL VOIGT (SBARv,SBARc,AUX66,AUX3333,2)
      CALL VOIGT (DBARv,DBARc,AUX66,AUX3333,2)

      if (iregime.eq.1) then
         CALL CHG_BASIS(AUX5,DBARC,AUX55,AUX3333,2,5)
         DBAR(1:MDIM)=AUX5(:)
         CALL CHG_BASIS(AUX5,SBARC,AUX55,AUX3333,2,5)
         SBAR(1:MDIM)=AUX5(:)
      else
         CALL CHG_BASIS(DBAR,DBARC,AUX66,AUX3333,2,6)
         CALL CHG_BASIS(SBAR,SBARC,AUX66,AUX3333,2,6)
      endif

      dsbar(:)=sbar(:)/nsteps   ! not sure if this is relevant for any problem?
      SBARNORM=VNORM(SBAR,6)
      DBARNORM=VNORM(DBAR,6)

C *** DEFINES CONTROLLING INCREMENT AND ESTIMATES TIME INCREMENT.
C *** ULTIMATE 'TINCR' WILL BE CALCULATED AT THE END OF EACH STEP.

      IF(ICTRL.EQ.0) THEN
         if (iregime.ne.1) then
            write(*,*)' ** ICTRL=0 WILL IMPOSE A VON MISES STRAIN'
            write(*,*)' INCREMENT. THIS OPTION VALID ONLY FOR FULLY'
            write(*,*)' IMPOSED STRAIN TENSOR AND VPSC (iregime=1).'
            call xit()
         endif
         STRAIN_CONTROL=1
         IF(IDBARSUM.NE.6. OR. INTERACTION.EQ.-1) THEN
            WRITE(*,*)'CAN CONTROL VM ONLY FOR PLASTIC SIMULATION'
            WRITE(*,*)'AND WHEN ALL STRAIN COMPONENTS ARE IMPOSED !'
            call xit()
         ENDIF
         EVMINCR=CTRLINCR
         DTIME=EVMINCR/VNORM(DBAR,5)
      ELSE IF(ICTRL.GE.1 .AND. ICTRL.LE.6) THEN
         IF(IDBARV(ICTRL).EQ.1) STRAIN_CONTROL=1
         IF(iSBARv(ICTRL).EQ.1) STRAIN_CONTROL=0
         IF (STRAIN_CONTROL.EQ.1)THEN
            DBARCTRL =DBARv(ictrl)
            EIJINCR=CTRLINCR
            DTIME=EIJINCR/ABS(DBARCTRL)
         endif
      ELSE IF(ICTRL.EQ.7) THEN
         IF(DBARNORM.NE.0.) STRAIN_CONTROL=1
      ENDIF

      tempdot=(temp_fin-temp_ini)/(dtime*nsteps)

      RETURN
      END subroutine load_conditions
