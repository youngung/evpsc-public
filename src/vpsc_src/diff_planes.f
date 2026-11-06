C **********************************************************************
C     SUBROUTINE DIFF_PLANES_PX    -->   VERSION 15/MAY/2025
C **********************************************************************
C *** IOPTION=0: FOR EACH PHASE OPEN DIFFRACTION FILE & READ DIFF
C                DIRECTIONS & MILLER INDICES OF DIFFRACTING PLANES
C *** IOPTION=1: FOR EACH PHASE DETERMINE GRAINS IN EACH HKL SET AND
C                CALCULATE AVERAGE LATTICE STRAIN & INTENSITY FOR THE SET
C                SAVE THE RESULTS (EELGR AND EPS_W) TO 'DIF_STR_PHX.OUT'
C                AND 'DIF_WGT_PHX.OUT' FILES, RESPECTIVELY.
C                ALSO, THE GRAINS OF ALL POLYCRYSTAL, REGARDLESS OF ITS PHASE,
C                ARE SUBJECTED TO THE CALCULATION OF PHASE 1'S DIFFRACTION
C                CONDITIONS AND STORED THEM TO THE 0-TH AXIS OF CERTAIN ARRAYS
C                AND USE THEM WITH IOPTION=2.
C *** IOPTION=2: CALCULATE AVERAGE LATTICE STRAIN & INTENSITY FOR THE SET, BUT
C                APPLY THE 0-TH AXIS RESULTS OF CERTAIN ARRAYS
C                (I.E., IGRSETPH(NDIFMX,NGRMX,0) AND WGTSETPH(NDIFMX,0))
C                THIS OPTION IS RELEVANT FOR THE CASE OF REXTAL AND VFT TWIN
C                CALCULATION, WHERE TWO PHASES ARE ACTUALLY IDENTICAL PHASES
C                WITH THE SAME CRYSTAL STRUCTURE AND DIFFRACTION CONDITIONS.
C                THE RESULTS ARE STORED TO "DIF_STR_ALL.OUT" AND "DIF_WGT_ALL.OUT".
C ***********************************************************************
C   NDIFF(IPH): # OF CRYSTALLOGRAPHIC DIFFRACTION PLANES HKL TO PROBE FOR
C               IN PHASE 'IPH'
C   ND=1,NDIFF(IPH): INDEX OF INDIVIDUAL DIFFRACTION DIRECTION AND
C                    CORRESPONDING PLANE HKL
C   NFAMILY(ND,IPH): # OF EQUIV XTALL PLANES HKL TO SCAN FOR DIRECTION ND
C   DIFF_VS(3,ND,IPH): DIRECTION IN SAMPLE AXES TO PROBE FOR PLANE HKL
C   DIFF_VC(3,NF,ND,IPH): DIRECTION IN CRYSTAL AXES OF EQUIV PLANES HKL
C   EPS_W(ND,IPH): AVERAGE NORMAL STRAIN FOR ALL PLANES HKL ALONG DIR ND
C ***********************************************************************

      SUBROUTINE DIFF_PLANES_PX(
C     INTENT(IN)::
     $     IOPTION,IDIFF,UR_DIF,UW_DIF,IFE,NPH,
     $     NGRMX,NPHMX,NDIFMX,NGR,AG,
     $     WGT,FILEDIFF,ICRYSYMPH,IREGIME,iDBARv,SCAUCHY,SAV,BC_EL,
     $     MGEL,EPSTOT,WPH,
C     INTENT(INOUT)::
     $     NDIFF,SG,ANGDETECTOR,NFAMILY,DIFF_VS,DIFF_VC,NGRSETPH,
     $     WGTSETPH,IGRSET,
C     INTENT(OUT)::
     $     EELGR)
C     -----------------
      implicit none
      INTEGER,INTENT(IN)::NGRMX,NPHMX,
     $     IOPTION,IDIFF(NPHMX),UR_DIF,IFE,NPH,
     $     NDIFMX,NGR(0:NPHMX),UW_DIF,IREGIME,iDBARv(6),ICRYSYMPH(NPHMX)
      INTEGER,INTENT(INOUT)::NDIFF(NPHMX),NFAMILY(NDIFMX,NPHMX),
     $     NGRSETPH(NDIFMX,0:NPHMX),IGRSET(NDIFMX,NGRMX,0:NPHMX)
      CHARACTER*120,INTENT(IN):: FILEDIFF(NPHMX)
      REAL*8, INTENT(IN)::AG(3,3,NGRMX),WGT(NGRMX),SCAUCHY(3,3),SAV(6),
     $     BC_EL(6,6,NGRMX),MGEL(6,6,NGRMX),EPSTOT(3,3),WPH(NPHMX)
      REAL*8,INTENT(INOUT)::ANGDETECTOR(NPHMX),DIFF_VS(3,NDIFMX,NPHMX),
     $     DIFF_VC(3,24,NDIFMX,NPHMX),WGTSETPH(NDIFMX,0:NPHMX),
     $     SG(6,NGRMX)
      REAL*8,INTENT(OUT)::EELGR(6,NGRMX)

C     --LOCALS
      CHARACTER*120 PROSA
      REAL*8 PS(3),AUX33(3,3),AUX6(6),AUX66(6,6),AUX3333(3,3,3,3),SN(3),
     $     SNEQ(3,24),SB(4),EPS,COSPSVS,COSANGDET,PI,SAVC(3,3),
     $     PARA_W(NDIFMX),EPS_W(NDIFMX,0:NPHMX),WTOT
      INTEGER IPH,ISN(4),ISB(4),NG,ND,NF,NGRX,I,J,IDIFFX,NEQUIV,IUNIT,
     $     IDUM,IPL,IGRSETX,ICASE
      CHARACTER*12 PLANE_INDICES(100)
c     FOR VPSC, ICASE=1 CALCUALTES SG6 FROM PRESSURE
c     ICASE=2 CALCULATES ALL SG COMPONENTS FROM ELASTICITY
      PARAMETER(ICASE=1)
C
C     IREGIME: 1('VPSC') OR 2('DELTA-EVPSC'), OR 3 ('SIGMA-EVPSC')
      LOGICAL VERBOSE
C      PARAMETER(VERBOSE=.TRUE.)
      PARAMETER(VERBOSE=.FALSE.)

      IF (VERBOSE) WRITE(*,*)'INSIDE SUBR. DIFF_PLANES'
      PI=4.D0*DATAN(1.D0)
C **********************************************************************
      IF (IOPTION.EQ.0) THEN
C **********************************************************************
C *** READS INDICES OF DIFFRACTING PLANE AND ANGLES OF DIFFRACTION DIRECTION.
C *** GENERATES CRYSTALLOGRAPHICALLY EQUIVALENT PLANE NORMALS 'SNEQ' IN CAXES.
C *** CALCULATES DIFFRACTION DIRECTION 'SB' IN SAXES
         DO 100 IPH=1,NPH
            IF (IDIFF(IPH).EQ.0) GOTO 100 !! SKIP
            OPEN(UNIT=UR_DIF,FILE=FILEDIFF(IPH),STATUS='OLD')
            IF (IFE.EQ.0) THEN  ! IN THE CASE OF SA CALCULATION
               IUNIT=UW_DIF+IPH
               OPEN(IUNIT,FILE='DIF_STR_PH'//CHAR(48+IPH)//'.OUT',
     $              STATUS='UNKNOWN')
               IUNIT=UW_DIF+IPH+5
               OPEN(IUNIT,FILE='DIF_WGT_PH'//CHAR(48+IPH)//'.OUT',
     $              STATUS='UNKNOWN')
C              OPEN FILES AGAIN FOR 'ALL' GRAINS (REX, VFT)
               IF(IPH.EQ.NPH.AND.NPH.GT.1) THEN
                  IUNIT=UW_DIF
                  WRITE(*,*)'IUNIT FOR DIF_STR_ALL.OUT:',IUNIT
                  OPEN(IUNIT,FILE='DIF_STR_ALL.OUT',
     $                 STATUS='UNKNOWN')
                  IUNIT=UW_DIF+5
                  WRITE(*,*)'IUNIT FOR DIF_WGT_ALL.OUT:',IUNIT
                  OPEN(IUNIT,FILE='DIF_WGT_ALL.OUT',
     $                 STATUS='UNKNOWN')
               ENDIF
            ENDIF
C -----------------------------------------------------------------
C     WRITES MILLER INDICES OF DIFFRACTING PLANES IN OUTPUT FILES
            READ(UR_DIF,'(A)') PROSA
            READ(UR_DIF,'(A)') PROSA
            READ(UR_DIF,*) NDIFF(IPH)
            READ(UR_DIF,'(A)') PROSA
            READ(UR_DIF,'(A)') PROSA
            DO IDIFFX=1,NDIFF(IPH)
               READ(UR_DIF,'(A12)') PLANE_INDICES(IDIFFX)
            ENDDO

            IF (IFE.EQ.0) THEN
               IUNIT=UW_DIF+IPH ! DIF_STR_PH
               WRITE(IUNIT,'(7A12,6X)')'EPSTOT11','EPSTOT22',
     $              'EPSTOT33','SAV11','SAV22','SAV33','WPH'
               IUNIT=UW_DIF+IPH+5 ! DIF_WGT_PH
               WRITE(IUNIT,'(7A12,6X)')'EPSTOT11','EPSTOT22',
     $              'EPSTOT33','SAV11','SAV22','SAV33','WPH'

               IF(IPH.EQ.NPH.AND.NPH.GT.1) THEN
                  IUNIT=UW_DIF ! DIF_STR_PH
                  WRITE(*,*)'IUNIT FOR DIF_STR_ALL.OUT:',IUNIT
                  WRITE(IUNIT,'(7A12,6X)')'EPSTOT11','EPSTOT22',
     $                 'EPSTOT33','SAV11','SAV22','SAV33','WTOT'
                  IUNIT=UW_DIF+5 ! DIF_WGT_PH
                  WRITE(*,*)'IUNIT FOR DIF_WGT_ALL.OUT:',IUNIT
                  WRITE(IUNIT,'(7A12,6X)')'EPSTOT11','EPSTOT22',
     $                 'EPSTOT33','SAV11','SAV22','SAV33','WTOT'
               ENDIF
            ENDIF

            REWIND UR_DIF
C -----------------------------------------------------------------
            READ(UR_DIF,'(A)') PROSA
            READ(UR_DIF,'(A)') PROSA
            READ(UR_DIF,*) NDIFF(IPH),ANGDETECTOR(IPH)
            IF (NDIFF(IPH).GT.NDIFMX) THEN
               WRITE(*,'('' NUMBER OF DIFFRACTION PLANES IS'',I4,/,
     $'' GREATER THAN CODE DIMENSION'',I4)') NDIFF(IPH),NDIFMX
               WRITE(*,*)
               WRITE(*,'(''STOP IN ROUTINE *** DIF_PLANES ***'')')
               CALL XIT()
            ENDIF
            READ(UR_DIF,'(A)') PROSA
            READ(UR_DIF,'(A)') PROSA

            DO ND=1,NDIFF(IPH)
               IDUM=ICRYSYMPH(IPH)
               CALL CRYSTAL_SYMMETRY(3,UR_DIF,IDUM,ISN,SN,
     $              SNEQ,ISB,SB,NEQUIV)
               NFAMILY(ND,IPH)=NEQUIV
               DO I=1,3
                  DIFF_VS(I,ND,IPH)=SB(I)
                  DO NF=1,NFAMILY(ND,IPH)
                     DIFF_VC(I,NF,ND,IPH)=SNEQ(I,NF)
                  ENDDO
               ENDDO
            ENDDO               ! NDIFF
            CLOSE(UNIT=UR_DIF)
 100     ENDDO                  ! IPH LOOP 1,NPH
      ENDIF   ! END OF INITIALIZATION RUN FOR IOPTION=0

C **********************************************************************
      IF (IOPTION.EQ.1) THEN
C **********************************************************************
         IF (VERBOSE) WRITE(*,*)'DIFF IOPT.EQ.1'

C *** IDENTIFIES AND TAGS THOSE GRAINS THAT CONTRIBUTE TO EACH DIFFRACTION
C     DIRECTION (ND) IN THE PHASE (IPH) --> IGRSET(ND,NGRSETPH,IPH)

         NGRSETPH(:,0)=0
         WGTSETPH(:,0)=0.

         DO 120 IPH=1,NPH
            IF (IDIFF(IPH).EQ.0) GOTO 120 !! SKIP
            COSANGDET=ABS(COS(ANGDETECTOR(IPH)*PI/180.0))
            IF (VERBOSE) WRITE(*,*)'ANGDETECTOR(IPH)',ANGDETECTOR(IPH)
C           DETERMINES NGRSETPH AND WGTSETPH
            DO 110 ND=1,NDIFF(IPH)
               NGRSETPH(ND,IPH)=0
               WGTSETPH(ND,IPH)=0.
               IF (VERBOSE) WRITE(*,'(A9,3F9.4)')
     $              'DIFF:',(DIFF_VS(I,ND,IPH),I=1,3)
               DO 105 NG=NGR(IPH-1)+1,NGR(IPH)
                  DO IPL=1,NFAMILY(ND,IPH)
                     DO I=1,3
                        PS(I)=0.0
                     DO J=1,3
                        PS(I)=PS(I)+AG(I,J,NG)*DIFF_VC(J,IPL,ND,IPH)
                     ENDDO
                     ENDDO

                     COSPSVS=0.
                     DO I=1,3
                        COSPSVS=COSPSVS+DIFF_VS(I,ND,IPH)*PS(I)
                     ENDDO
                     COSPSVS=ABS(COSPSVS)

                     IF(COSPSVS.GE.COSANGDET) THEN
                        NGRSETPH(ND,IPH)=NGRSETPH(ND,IPH)+1
                        NGRSETPH(ND,0)=NGRSETPH(ND,0)+1

C                       FOR EACH PHASE
                        NGRX=NGRSETPH(ND,IPH)
C                       IN BELOW, THE THIRD AXIS OF IGRSET EXTENDS FROM 1 TO 4.
                        IGRSET(ND,NGRX,IPH)=NG
                        WGTSETPH(ND,IPH)=WGTSETPH(ND,IPH)+WGT(NG)
C                       FOR ENTIRE PHASE
                        NGRX=NGRSETPH(ND,0)
C                       IN BELOW, THE THIRD AXIS OF IGRSET EXTENDS FROM 1 TO 4.
                        IGRSET(ND,NGRX,0)=NG
                        WGTSETPH(ND,0)=WGTSETPH(ND,0)+WGT(NG)
                     ENDIF
                  ENDDO         ! END OF IPL
 105           ENDDO            ! END OF NG

               IF (VERBOSE) WRITE(*,'(''* WGT OF EACH SET:'',20F10.5)')
     $              (WGTSETPH(I,IPH),I=1,NDIFF(IPH))
               IF (VERBOSE) WRITE(*,'(''* NGR OF EACH SET:'',20I5)')
     $              (NGRSETPH(I,IPH),I=1,NDIFF(IPH))
 110        ENDDO               ! END OF ND
 120     ENDDO


C **********************************************************************
C *** AVERAGE CRYSTALLOGRAPHIC STRAIN FOR GRAINS IN EACH SET
C *** WHAT FOLLOWS ARE DIFFERENT WAYS OF CALCULATING ELASTIC STRAIN IN EACH
C     GRAIN DEPENDING ON THE TYPE OF PROCESS BEING SIMULATED
C     ALREADY UPDATED STRESS FROM SUBROUTINE <VPSC>
         SAVC(:,:)=SCAUCHY(:,:)
C     NOW SAV IS IN 6TH DIMENSION FOR <VPSC>
         IF (IREGIME.EQ.1) CALL CHG_BASIS(SAV,SAVC,AUX66,AUX3333,2,6)

         DO 200 IPH=1,NPH
            IF (IDIFF(IPH).EQ.0) GOTO 200 !! SKIP
            DO ND=1,NDIFF(IPH)
               PARA_W(ND)=0.D0
               DO NG=1,NGRSETPH(ND,IPH)
                  IGRSETX=IGRSET(ND,NG,IPH)
C               *** ICASE 1 OR 2 ***
C               PERFORMED ONLY WHEN SG(6) IS MISSING AS IN VPSC
                  IF (IREGIME.EQ.1) THEN !! FOR VPSC
                     IF (iDBARv(1)*iDBARv(2)*iDBARv(3).EQ.1) THEN ! VPSC
                        WRITE(*,*)'---------- SHORT NOTICE ------------'
                        WRITE(*,*)'THE FEATURE OF DIFF_PLANES TO '
                        WRITE(*,*)'CALCULATE INTERNAL STRAIN EELGR'
                        WRITE(*,*)'IS AVAILABLE FOR THE CASE OF '
                        WRITE(*,*)'iDBARv(1)*iDBARv(2)*iDBARv(3).NE.1'
                        WRITE(*,*)'WITH THE GIVEN BOUNDARY CONDITION'
                        WRITE(*,*)'VPSC CANNOT OBTAIN EELGR SINCE THE'
                        WRITE(*,*)' MACROSCOPIC PRESSURE IS UNKNOWN.'
                        WRITE(*,*)'------------------------------------'
                        CALL XIT()
                     else
                        IF (ICASE.EQ.1) THEN !! WHEN ONLY CALC. SG6
                           SG(6,IGRSETX)=0.D0
                           DO J=1,6
                              SG(6,IGRSETX)=SG(6,IGRSETX)
     $                             +BC_EL(6,J,IGRSETX)*SAV(J)
                           ENDDO
                        ELSEIF (ICASE.EQ.2) THEN !! WHEN CALC. ALL SG(:)
                           SG(:,IGRSETX)=0.D0
                           DO I=1,6
                              DO J=1,6
                                 SG(I,IGRSETX)=SG(I,IGRSETX)
     $                                +BC_EL(I,J,IGRSETX)*SAV(J)
                              ENDDO
                           ENDDO
                        ENDIF
                     endif
                  ENDIF

C *** CALCULATION OF EELGR
                  EELGR(:,IGRSETX)=0.
                  DO I=1,6
                  DO J=1,6
                     EELGR(I,IGRSETX)=EELGR(I,IGRSETX)+MGEL(I,J,IGRSETX)
     $                    *SG(J,IGRSETX)
                  ENDDO
                  ENDDO

C     TO CARTESIAN
                  AUX6(:)=EELGR(:,IGRSETX)
                  CALL CHG_BASIS(AUX6,AUX33,AUX66,AUX3333,1,6)
C     PROJECTING EEL STRAIN
                  EPS=0.D0
                  DO I=1,3
                  DO J=1,3
                     EPS=EPS+DIFF_VS(I,ND,IPH)*DIFF_VS(J,ND,IPH)
     $                    *AUX33(I,J)
                  ENDDO
                  ENDDO
C     WEIGHTED USING GRAIN FRACTION
                  PARA_W(ND)=PARA_W(ND)+EPS*WGT(IGRSETX)
               ENDDO            ! END OF 'DO NG' OVER GRAINS IN SUBSET
C            IF(WGTSETPH(ND,IPH).NE.0.D0)
               IF(WGTSETPH(ND,IPH).GE.1E-15)
C                RENORMALIZE USING SET FRACTION
     $              EPS_W(ND,IPH)=PARA_W(ND)/WGTSETPH(ND,IPH)
            ENDDO               ! END OF 'DO ND' OVER DIFFRACTING PLANES

            IF (IFE.EQ.0) THEN
               IUNIT=UW_DIF+IPH
               WRITE(IUNIT,'(7E12.4,6X,100E12.4)')EPSTOT(1,1),
     $              EPSTOT(2,2),EPSTOT(3,3),SAVC(1,1),SAVC(2,2),
     $              SAVC(3,3),WPH(IPH),(EPS_W(ND,IPH),ND=1,NDIFF(IPH))
               IUNIT=UW_DIF+IPH+5
               WRITE(IUNIT,'(7E12.4,6X,100E12.4)')EPSTOT(1,1),
     $              EPSTOT(2,2),EPSTOT(3,3),SAVC(1,1),SAVC(2,2),
     $             SAVC(3,3),WPH(IPH),(WGTSETPH(ND,IPH),ND=1,NDIFF(IPH))
            ENDIF
 200     ENDDO                  !END OF IPH LOOP OVER NPH
      ENDIF                     !END OF IOPTION=1

C **********************************************************************
      IF (IOPTION.EQ.2) THEN
         SAVC(:,:)=SCAUCHY(:,:)
         DO 300 ND=1,NDIFF(1)!USING THE DIFFRACTION CONDITION OF PHASE 1
            PARA_W(ND)=0.D0
            DO NG=1,NGRSETPH(ND,0)
               IGRSETX=IGRSET(ND,NG,0)
C *** CALCULATION OF EELGR
               EELGR(:,IGRSETX)=0.
               DO I=1,6
               DO J=1,6
                  EELGR(I,IGRSETX)=EELGR(I,IGRSETX)+MGEL(I,J,IGRSETX)
     $                 *SG(J,IGRSETX)
               ENDDO
               ENDDO

C     TO CARTESIAN
               AUX6(:)=EELGR(:,IGRSETX)
               CALL CHG_BASIS(AUX6,AUX33,AUX66,AUX3333,1,6)
C     PROJECTING EEL STRAIN
               EPS=0.D0
               DO I=1,3
               DO J=1,3
                  EPS=EPS+DIFF_VS(I,ND,1)*DIFF_VS(J,ND,1)
     $                 *AUX33(I,J)
               ENDDO
               ENDDO
C     WEIGHTED USING GRAIN FRACTION
               PARA_W(ND)=PARA_W(ND)+EPS*WGT(IGRSETX)
            ENDDO               ! END OF 'DO NG' OVER GRAINS IN SUBSET
C            IF(WGTSETPH(ND,IPH).NE.0.D0)
            IF(WGTSETPH(ND,0).GE.1E-15)
C                RENORMALIZE USING SET FRACTION
     $           EPS_W(ND,0)=PARA_W(ND)/WGTSETPH(ND,0)
 300     ENDDO                  ! END OF 'DO ND' OVER DIFFRACTING PLANES

         IF (IFE.EQ.0) THEN
            IUNIT=UW_DIF
            WTOT=0.
            DO IPH=1,NPH
               WTOT=WTOT+WPH(IPH)
            ENDDO
            WRITE(IUNIT,'(7E12.4,6X,100E12.4)')EPSTOT(1,1),
     $           EPSTOT(2,2),EPSTOT(3,3),SAVC(1,1),SAVC(2,2),
     $           SAVC(3,3),WTOT,(EPS_W(ND,0),ND=1,NDIFF(1))
            IUNIT=UW_DIF+5
            WRITE(IUNIT,'(7E12.4,6X,100E12.4)')EPSTOT(1,1),
     $           EPSTOT(2,2),EPSTOT(3,3),SAVC(1,1),SAVC(2,2),
     $           SAVC(3,3),WTOT,(WGTSETPH(ND,0),ND=1,NDIFF(1))
         ENDIF
      ENDIF                     ! END OF IOPTION=2

      RETURN
      END SUBROUTINE DIFF_PLANES_PX
