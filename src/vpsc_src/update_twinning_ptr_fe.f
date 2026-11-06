C ************************************************************************
C    SUBROUTINE UPDATE_TWINNING_PTR_fe   --->   VERSION OF 21/APR/2020
C
C --> SEP/2005 (CNT): NOW THE STATS ON TWINNING ARE DONE INSIDE
C     1st DO LOOP. THE TAGGING FOR REORIENTATION IS DONE INSIDE 2nd LOOP.
C --> 14/NOV/2010: OPTION FOR TYPE I or TYPE II TWIN REORIENTATION ADDED
C --> 26/MAR/2017: A MONTE-CARLO SELECTION (BASED ON PROBABILITY OF
C     RANDOMLY PICKING THE PTS) WAS ADDED OLD PTS SELECTION (BASED ON
C     MAXIMUM TWIN VOL FRACTION)
C --> 09/AUG/2017: PTS OR MONTE-CARLO SELECTION CONTROLED BY IPTS=1 OR 0
C ************************************************************************
C
C    IPTSGR(KGX):    INDEX OF PREDOMINANT TWIN SYSTEM IN THE GRAIN.
C    IPTMGR(KGX):    INDEX OF PREDOMINANT TWIN MODE IN THE GRAIN.
C    TWFRGR(TWM):    ACCUMULATED TWIN FRACTION IN EACH MODE IN GRAIN KKK
C                    IN THE STEP (relative to grain).
C
C    TWFRSY(TWS,KKK):ACCUMMULATED TWIN VOLUME FRACTION IN SYSTEM TWS
C                    IN GRAIN KKK (relative to grain).
C    TWFRPH(TWM):ACCUMULATED TWIN FRACT IN EACH MODE OF EACH PHASE
C                    (relative to phase volume).
C    EFTWFR(TWM):ACCUM TWIN FRACTION GIVEN BY THE REORIENTED GRAINS
C                    IN EACH MODE OF EACH PHASE (relative to phase volume).
C    KTWSMX(KKK):INDEX OF THE PREDOMINANT TWIN SYSTEM (MAX TWIN VOLUME)
C                IN THE GRAIN. USED FOR REORIENTING GRAIN.
C                   =0 IF THE GRAIN IS NOT TAGGED FOR REORIENTATION.
C                   >0 IF THE GRAIN IS TAGGED FOR REORIENTATION.
C                   <0 IF THE GRAIN HAS BEEN REORIENTED AND IS NOT TO BE
C                      RETAGGED FOR SECONDARY TWINNING
C
C    PRITW(IPH) :ACCUMMULATED TWIN FRACT IN EACH PHASE ASSOCIATED WITH
C                GRAINS THAT WERE TWIN REORIENTED ONCE.
C    SECTW(IPH) :ACCUMMULATED TWIN FRACT IN EACH PHASE ASSOCIATED WITH
C                GRAINS THAT WERE TWIN REORIENTED MORE THAN ONCE.
C    NTWEVENTS(KKK) :NUMBER OF TIMES THAT GRAIN 'KKK' HAS BEEN REORIENTED.
C                    0: NONE, 1: PRIMARY, 2:SECONDARY, 3:TERTIARY
C ************************************************************************

      SUBROUTINE UPDATE_TWINNING_PTR_fe(
c          in
     $     tincr,iph,igb,igt,nphmx,wph,ngr,ngrmx,nslsys,ntwmod,
     $     nslmod,twsh,nsm,nmodmx,ntwmmx,nsysmx,gamdot,ntwsmx,
     $     wgt,itwinlaw,jran,twthres,dnca,dbca,ntwsys,isectw,
     $     itwtype,scauchy,
c          inout
     $     twfrsy,twfrph,ktwsmx,eftwfr,ag,dsg,ntwevents,pritw,sectw,
     $     sch_stat,var_stat,var_vf
c          out
     $     )
      implicit none
      integer,intent(in)::nphmx,ngrmx,nmodmx,ntwmmx,ntwsmx,nsysmx
      integer,intent(in)::ngr,nslsys(nphmx),
     $     ntwmod(nphmx),nslmod(nphmx),nsm(nmodmx,nphmx),
     $     itwinlaw,jran,ntwsys(nphmx),
     $     isectw(nsysmx,nphmx),
     $     itwtype(nsysmx,nphmx),iph,igb,igt
      real*8,intent(in)::tincr,wph(nphmx),
     $     twsh(ntwmmx,nphmx),
     $     gamdot(nsysmx,ngrmx),
     $     wgt(ngrmx),twthres(2,ntwmmx,nphmx),
     $     dnca(3,nsysmx,nphmx),
     $     dbca(3,nsysmx,nphmx),
     $     scauchy(3,3)
      real*8,intent(inout)::twfrsy(ntwsmx,ngrmx),twfrph(ntwmmx,nphmx),
     $     eftwfr(ntwmmx,nphmx),ag(3,3,ngrmx),dsg(6,ngrmx),pritw(nphmx),
     $     sectw(nphmx),var_vf(6,ntwmmx,nphmx)


      integer,intent(inout)::ktwsmx(ngrmx),ntwevents(ngrmx),
     $     sch_stat(-5:15,ntwmmx,nphmx),var_stat(6,ntwmmx,nphmx)

c     locals

      integer IMARK(NGRmx),IPTSGR(NGRmx),IPTMGR(NGRmx)
      real*8 TWB(3),TWN(3),
     $     ACRYS(3,3),ATWIN(3,3),rand,gmode,twfmax,
     $     gabs,twshx,twfrthres,twfradd,ran2,aux33(3,3),bux33(3,3),
     $     aux6(6),aux3333(3,3,3,3),aux66(6,6),TWFRGR(0:NTWMMX,NGRmx)
      integer nsmx,kleft,i,index,kacum,iptm,ipts,ktw,j,k,l,is,
     $     iskip,ktmode,itm,kkk,kts,kgx,its,kpoint,igr,
     $     ngphase
C      DIMENSION aux5(5),aux33(3,3),aux55(5,5),aux3333(3,3,3,3)

      IF(WPH(iph).EQ.0.) RETURN

C *** UPDATES THE VOLUME FRACTION ASSOCIATED WITH EACH TWIN SYSTEM IN EACH GRAIN
      KGX=1
      DO KKK=igb,igt            !1,NGR
        KTS=0
        ITS=NSLSYS(iph)
        TWFRGR(0,KGX)=0.

        DO ITM=1,NTWMOD(iph)
          KTMODE=ITM+NSLMOD(iph) ! loop only over twin modes.
          GMODE=0.
          TWFRGR(ITM,KGX)=0.
          TWSHX=TWSH(ITM,iph)

          DO NSMX=1,NSM(KTMODE,iph)
            KTS=KTS+1             ! shifted counter for twin systems
            ITS=ITS+1             ! absolute counter over all systems
            IF(GAMDOT(ITS,KGX).GT.0.) THEN
              GABS=ABS(GAMDOT(ITS,KGX))*TINCR
              GMODE=GMODE+GABS

c     twin fr in syst relative to grain volume
              TWFRSY(KTS,KKK)=TWFRSY(KTS,KKK)+GABS/TWSHX
            ENDIF
            TWFRGR(ITM,KGX)=TWFRGR(ITM,KGX)+TWFRSY(KTS,KKK)
          ENDDO

C * total twin fraction in grain relative to grain volume
          TWFRGR(0,KGX)=TWFRGR(0,KGX)+TWFRGR(ITM,KGX)
C * total twin fraction in aggregate for each twin mode relative to phase volume
          TWFRPH(ITM,iph)=TWFRPH(ITM,iph)+
     #                      (GMODE/TWSHX)*(WGT(KKK)/WPH(iph))


        ENDDO      ! END OF LOOP OVER TWIN MODES IN THE PHASE

C *** PICKS EITHER OLD PTR SCHEME (ITWINLAW=1) OR NEW MONTE CARLO SCHEME (ITWINLAW=2)
C ----------------------------------------------------------------------
        IF (ITWINLAW.EQ.1) THEN
C ----------------------------------------------------------------------
C *** FIND THE PREDOMINANT TWIN SYSTEM (PTS) IN EACH GRAIN (OLD SCHEME)

        IPTSGR(KGX)=0
        IPTMGR(KGX)=0
        TWFMAX=0.
        KTS=0
        ITS=NSLSYS(iph)

        DO ITM=1,NTWMOD(iph)    ! relative index
          KTMODE=ITM+NSLMOD(iph)
          DO NSMX=1,NSM(KTMODE,iph)
            KTS=KTS+1             ! shifted counter for twin systems
            ITS=ITS+1             ! absolute counter over all systems
            rand=1.0
            rand=0.5+0.5*ran2(jran)   ! adds stochasticity to selection
            IF (TWFRSY(KTS,KKK)*rand.GT.TWFMAX) THEN
              IPTSGR(KGX)=ITS     ! twin system
              IPTMGR(KGX)=ITM     ! twin mode
              TWFMAX=TWFRSY(KTS,KKK)
            ENDIF
          ENDDO
        ENDDO      ! END OF LOOP OVER TWIN MODES IN THE PHASE
C --------------------------------------------------------------------
        ELSE IF (ITWINLAW.EQ.2) THEN
C --------------------------------------------------------------------
C *** MONTE CARLO SELECTION OF PTS IN EACH GRAIN

        IPTSGR(KGX)=0
        IPTMGR(KGX)=0
        KTS=0
        ITS=NSLSYS(iph)

        TWFRADD=0.    ! normalized twin fraction in the grain
        RAND=RAN2(jran)   ! random number in [0,1] interval
        iskip=0
        DO ITM=1,NTWMOD(iph)    ! relative index
           KTMODE=ITM+NSLMOD(iph)
           DO NSMX=1,NSM(KTMODE,iph)
              KTS=KTS+1         ! shifted counter for twin systems
              ITS=ITS+1         ! absolute counter over all systems
              TWFRADD=TWFRADD+TWFRSY(KTS,KKK)/TWFRGR(0,KGX)
              IF (TWFRADD.GT.RAND) THEN
                 IPTSGR(KGX)=ITS ! twin system
                 IPTMGR(KGX)=ITM ! twin mode
                 ISKIP=1
              ENDIF
              IF(ISKIP.EQ.1) GO TO 10
           ENDDO
c       END OF LOOP OVER TWIN MODES IN THE PHASE
        ENDDO
 10     CONTINUE
C ----------------------------------------------------------------------
        ENDIF
C --------------------------------------------------------------------
        KGX=KGX+1
      ENDDO      ! END OF DO KKK OVER GRAINS IN THE PHASE

C *********************************************************************
C    THIS LOOP SCANS ALL GRAINS IN THE PHASE/ELEM IN A RANDOM WAY AND
C    TAGS THEM FOR CHECKING AGAINST THE TWIN-REORIENTATION CRITERION.
C    IMARK(1:NGPHASE) TAGS WITH '1' THE GRAINS AS THEY ARE PICKED.
C    IN EACH LOOP ITERAT A GRAIN 'KKK' IS RANDOMLY PICKED AND REORIENTED
C    ABOUT THE PTS IF THE 'REAL' TWIN FRACTION IS NOT EXCEEDED.
C *********************************************************************

c      NGPHASE=NGR(IPH)-NGR(IPH-1)
      NGPHASE=NGR
      KLEFT  =NGPHASE
      DO I=1,NGPHASE
        IMARK(I)=0
      ENDDO

      DO INDEX=1,NGPHASE    ! tags one of the yet unchecked grains
        RAND=RAN2(JRAN)
        KPOINT=INT(RAND*KLEFT)+1
        KLEFT=KLEFT-1
        KACUM=0
        IGR=0

c                                     KPOINT: index of tagged grain
        DO WHILE(IGR.LT.NGPHASE .AND. KPOINT.NE.0)
          IGR=IGR+1
          IF(IMARK(IGR).EQ.0) THEN
            KACUM=KACUM+1
            IF(KACUM.EQ.KPOINT) THEN
c              KKK=IGR+NGR(IPH-1)     ! absolute index of grain
               KKK=IGR               ! absolute index of grain
c            index of grain relative to element
              KGX=IGR
              IMARK(IGR)=1
              IPTS=IPTSGR(KGX)       ! absolute index
              IPTM=IPTMGR(KGX)       ! relative index

c     nslsyx=nslsys(iphel)

C     LABELS THE GRAIN IF THE ACCUMULATED TWIN VOLUME IN THE PREDOMINANT
C     TWINNING SYSTEM (PTS) HAS EXCEEDED A THRESHOLD VALUE.
C     THIS GRAIN WILL BE COMPLETELY REORIENTED BY TWINNING.
C     REORIENTATION IS STOPPED WHEN THE EFFECTIVE REORIENTED FRACTION
C     REACHES THE 'REAL' TWIN FRACTION CALCULATED FROM THE SHEARS.
C     THRES1: MIN ACCUM TWIN FRACTION IN GRAIN BEFORE TWIN REOR IS SWITCHED ON.
C     THRES1+THRES2: EVENTUAL ACCUM TWIN FR IN GR REQUIRED TO TWIN REORIENT IT.

              IF(KTWSMX(KKK).EQ.0 .AND. IPTS.NE.0) THEN
                IF (EFTWFR(IPTM,iph).LT.TWFRPH(IPTM,iph)) THEN
                  TWFRTHRES=TWTHRES(1,IPTM,iph)+
     $                  TWTHRES(2,IPTM,iph)*EFTWFR(IPTM,iph)/
     $                  (TWFRPH(IPTM,iph)+1.E-6)
                  IF (TWFRGR(IPTM,KGX).GT.TWFRTHRES) THEN
                    EFTWFR(IPTM,iph)=EFTWFR(IPTM,iph)+WGT(KKK)/WPH(iph)
                    KTWSMX(KKK)=IPTS


                  ENDIF
                ENDIF      ! enf of if ktwsms(kkk).eq.0

C               KTWSMX(KKK)=0   ! activate to suppress twin reorientation always

C *** REORIENTS BY TWINNING GRAINS TAGGED WITH KTWSMX>0 AND UPDATES THE
C     ORIENTATION MATRIX.
C *** THIS IS DONE AFTER HAVING CALLED UPDATE_ORIENTATION, WHICH HAS
C     REORIENTED THE GRAIN BECAUSE OF SLIP AND TWIN SHEAR

                IF (KTWSMX(KKK).GT.0) THEN

                  KTW=KTWSMX(KKK)
                  DO I=1,3
                    DO J=1,3
                      ACRYS(I,J)=AG(I,J,KKK)
                    ENDDO
                    TWN(I)=DNCA(I,KTW,iph)
                    TWB(I)=DBCA(I,KTW,iph)
                  ENDDO

                  aux33(:,:)=ag(:,:,kkk)
                  if (.false.) then
                     CALL STAT_TWINNING_fe(IPTS,IPTM,KKK,iph,nphmx,
     $                 nsysmx,ntwmmx,nslmod(iph),nmodmx,
     $                 nsm,
     $                 nslsys,
     $                 dnca,
     $                 dbca,
     $                 aux33,
     $                 scauchy,wgt(kkk),
c                      inout
     $                 sch_stat,!(-5:15,1:ntwmmx,iph),
     $                 var_stat,!(1:6,1:ntwmmx,iph),
     $                 var_vf!(1:6,1:ntwmmx,iph)
     $                 )        ! twinning statistics
                  endif

C *** CALCULATES CRYSTAL-TO-SAMPLE ORIENTATION MATRIX FOR THE TWINNED GRAIN
c                  ATWIN: CA parent <- CA child
                  CALL TWIN_ORIENTATION (ITWTYPE(KTW,iph),TWB,TWN,ATWIN)

c                 1. store delta sg converted from the previous step in crystal axes.
                  aux6(:)=dsg(:,kkk)
                  call chg_basis(aux6,bux33,aux66,aux3333,1,6)

c                 ag: sa<-ca(parent)
c                 aux33: dsg^(parent)
                  do i=1,3
                  do j=1,3
                     aux33(i,j)=0d0
                     do k=1,3
                     do l=1,3
                        aux33(i,j)=aux33(i,j)+
     $                       ag(k,i,kkk)*ag(l,j,kkk)*bux33(k,l)
                     enddo
                     enddo
                  enddo
                  enddo
c     --
c
c                 acrys: sa<-ca(parent)
c                 ag: sa<-ca(paranet) . ca(parent)<-ca(child)
                  DO I=1,3
                  DO J=1,3
                    AG(I,J,KKK)=0.
                    DO K=1,3
                      AG(I,J,KKK)=AG(I,J,KKK)+ACRYS(I,K)*ATWIN(K,J)
                    ENDDO
                  ENDDO
                  ENDDO

c                 transform sg stored in crystal axes to new sample coordinate system
                  do i=1,3
                  do j=1,3
                     bux33(i,j)=0d0
                     do k=1,3
                     do l=1,3
                        bux33(i,j)=bux33(i,j)+
     $                       ag(i,k,kkk)*ag(j,l,kkk)*aux33(k,l)
                     enddo
                     enddo
                  enddo
                  enddo

                  call chg_basis(aux6,bux33,aux66,aux3333,2,6)
                  dsg(:,kkk)=aux6(:)


C *** RESET ACCUMULATED TWINNED FRACTIONS IN THE GRAIN
                  DO IS=1,NTWSYS(iph)
                    TWFRSY(IS,KKK)=0.
                  ENDDO
CFEB
C *** ACCUMULATES TWINNED VOLUME FRACTION IN EACH PHASE.
                  IF(NTWEVENTS(KKK).EQ.0) THEN
                    PRITW(iph)=PRITW(iph)+WGT(KKK)/WPH(iph)
                  ELSE IF(NTWEVENTS(KKK).GT.0) THEN
                    SECTW(iph)=SECTW(iph)+WGT(KKK)/WPH(iph)
                  ENDIF
CFEE
                  NTWEVENTS(KKK)=NTWEVENTS(KKK)+1

c     SUPPRESES SECTW
                  IF(ISECTW(KTW,iph).EQ.0 .OR. NTWEVENTS(KKK).GT.1)
     $                 KTWSMX(KKK)=-KTWSMX(KKK)
! ALLOWS SECTW
                  IF(ISECTW(KTW,iph).EQ.1 .AND.NTWEVENTS(KKK).LE.1)
     $                 KTWSMX(KKK)= 0

C *** RESET CRSS IN SL & TW SYSTEMS TO INITIAL (OR OTHER) VALUES
C                 GTOTGR(KKK)=0
C                 DO IS=1,NSLSYS(iph)
C                   CRSS(IS,KKK)=TAU(IS,0,IPHEL)
C                 ENDDO
C                 DO IS=NSLSYS(iph)+1,NSYST
C                   CRSS(IS,KKK)=5.0*CRSS(IS,KKK)
C                 ENDDO

                ENDIF      ! END OF IF KTWSMX(KKK)>0
              ENDIF      ! END OF IF KTWSMX(KKK)=0
            ENDIF      ! END OF IF(KACUM.EQ.KPOINT)
          ENDIF      ! END OF IF(IMARK(IGR).EQ.0)
        ENDDO      ! END OF DO WHILE(IGR.LT.NGPHASE)
      ENDDO      ! END OF DO INDEX=1,NGPHASE

      RETURN
      END
