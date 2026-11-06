C **********************************************************************
C     SUBROUTINE WRITE_ACTIV_FE     --->      VERSION 28/MAR/23
C
C     WRITES ACTIVITY STATISTICS OVER DEFORMATION MODES IN 'ACT_PHN.OUT'
C     IOPT0: OPEN FILES AND WRITE COLUMN HEADS.
C     IOPT1: WRITES DATA
C **********************************************************************
      SUBROUTINE WRITE_ACTIV_FE(
     $     IOPT,IUNIT_ACT,UW_GAM,IPH,IHARDLAW,NMODES,NTWMOD,WPH,
     $     EPSVM,EPSTOT,TWFRPH,EFTWFR,PRITW,SECTW,ACTIV,GAVMOD,GAVACC,
     $     GAVPH)
      implicit none
      INTEGER, INTENT(IN)::IOPT,IUNIT_ACT,UW_GAM,IPH,IHARDLAW,NMODES,
     $     NTWMOD
      REAL*8, INTENT(IN):: WPH,EPSVM,EPSTOT(3,3),TWFRPH(NTWMOD),
     $     EFTWFR(NTWMOD),PRITW,SECTW,ACTIV,GAVMOD(NMODES),
     $     GAVACC(NMODES),GAVPH
C     LOCALS
      CHARACTER*8 MODE(NMODES),TWFR(NMODES),EFFR(NMODES)
      CHARACTER*2 FLABEL
      INTEGER I,J,JMO,KMO

      FLABEL='PH'
      DO I=1,NMODES
         MODE(I)='   MODE'//CHAR(48+I)
         TWFR(I)='   TWFR'//CHAR(48+I)
         EFFR(I)='   EFFR'//CHAR(48+I)
      ENDDO

      IF (IOPT.EQ.0) THEN
         OPEN(IUNIT_ACT,FILE='ACT_'//FLABEL//CHAR(48+IPH)//'.OUT',
     $        STATUS='UNKNOWN')
         OPEN(UW_GAM,FILE='GAM_'//FLABEL//CHAR(48+IPH)//'.OUT',
     $        STATUS='UNKNOWN')

         IF (NTWMOD.EQ.0) THEN  ! ONLY DISL. SLIPS
            WRITE(IUNIT_ACT,'(12A14)')'EPSVM','E11','E22','AVACS',
     $           'GAVPH',
     $           ('   MODE'//CHAR(48+I),I=1,NMODES)
            WRITE(UW_GAM,'(12A14)')'EPSVM','E11','E22','AVACS',
     $           'GAVPH',
     $           ('   MODE'//CHAR(48+I),I=1,NMODES)
         ELSEIF (NTWMOD.NE.0) THEN ! DISL. SLIPS + TWINS.
            IF (IHARDLAW.EQ.0.OR.IHARDLAW.EQ.20) THEN
               WRITE(IUNIT_ACT,'(32A14)')'EPSVM','E11','E22','AVACS',
     $              'GAVPH','PRITW','SECTW',
     $              ('   MODE'//CHAR(48+I),I=1,NMODES),
     $              ('   TWFR'//CHAR(48+J),'   EFFR'//CHAR(48+J),
     $              J=NMODES-NTWMOD+1,NMODES)
               WRITE(UW_GAM,'(32A14)')'EPSVM','E11','E22','AVACS',
     $              'GAVPH','PRITW','SECTW',
     $              ('   MODE'//CHAR(48+I),I=1,NMODES),
     $              ('   TWFR'//CHAR(48+J),'   EFFR'//CHAR(48+J),
     $              J=NMODES-NTWMOD+1,NMODES)
            ELSEIF(IHARDLAW.EQ.2) THEN
               WRITE(IUNIT_ACT,'(32A14)')'EPSVM','E11','E22','PHFRAC',
     $              'AVACS','GAVPH',
     $              ('   MODE'//CHAR(48+I),I=1,NMODES),
     $              ('   TWFR'//CHAR(48+J),
     $              J=NMODES-NTWMOD+1,NMODES)
               WRITE(UW_GAM,'(32A14)')'EPSVM','E11','E22','PHFRAC',
     $              'AVACS','GAVPH',
     $              ('   MODE'//CHAR(48+I),I=1,NMODES),
     $              ('   TWFR'//CHAR(48+J),
     $              J=NMODES-NTWMOD+1,NMODES)
            ENDIF
         ENDIF
         RETURN                 !! EXIT

      ELSEIF(IOPT.EQ.1) THEN
C *** FOR LOADING HISTORIES
         IF (NTWMOD.EQ.0) THEN
            WRITE(IUNIT_ACT,'(14E14.5E3)') EPSVM,EPSTOT(1,1),
     $           EPSTOT(2,2),ACTIV,GAVPH,(GAVMOD(KMO),KMO=1,NMODES)
            WRITE(UW_GAM,'(14E14.5E3)') EPSVM,EPSTOT(1,1),
     $           EPSTOT(2,2),ACTIV,GAVPH,(GAVACC(KMO),KMO=1,NMODES)
         ELSEIF (NTWMOD.NE.0) THEN
            IF(IHARDLAW.EQ.0.OR.IHARDLAW.EQ.20) THEN
c              ACT_PH?.OUT
               WRITE(IUNIT_ACT,'(36E14.5E3)')EPSVM,EPSTOT(1,1),
     $              EPSTOT(2,2),ACTIV,GAVPH,PRITW,SECTW,
     $              (GAVMOD(KMO),KMO=1,NMODES),
     $              (TWFRPH(JMO),EFTWFR(JMO),JMO=1,NTWMOD)
c              GAM_PH?.OUT
               WRITE(UW_GAM,'(36E14.5E3)')EPSVM,EPSTOT(1,1),
     $              EPSTOT(2,2),ACTIV,GAVPH,PRITW,SECTW,
     $              (GAVACC(KMO),KMO=1,NMODES),
     $              (TWFRPH(JMO),EFTWFR(JMO),JMO=1,NTWMOD)
            ELSEIF(IHARDLAW.EQ.2) THEN
c              ACT_PH?.OUT
               WRITE(IUNIT_ACT,'(36E14.5E3)') EPSVM,EPSTOT(1,1),
     $              EPSTOT(2,2),WPH,ACTIV,GAVPH,
     $              (GAVMOD(KMO),KMO=1,NMODES),
     $              (TWFRPH(JMO),JMO=1,NTWMOD)
c              GAM_PH?.OUT
               WRITE(UW_GAM,'(36E14.5E3)') EPSVM,EPSTOT(1,1),
     $              EPSTOT(2,2),WPH,ACTIV,GAVPH,
     $              (GAVACC(KMO),KMO=1,NMODES),
     $              (TWFRPH(JMO),JMO=1,NTWMOD)
            ENDIF
         ENDIF
      ENDIF                     ! END OF IFS on IOPT

      RETURN
      END SUBROUTINE WRITE_ACTIV_FE
