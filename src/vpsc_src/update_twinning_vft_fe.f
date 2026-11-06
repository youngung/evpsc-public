      SUBROUTINE UPDATE_TWINNING_VFT_fe(
c     in
     $     tincr,IPH,nphmx,ngrmx,ngr,nslmod,nslsys,nsysmx,ntwmod,ntwmmx,
     $     ntwsmx,twthres,nsm,gamdot,nmodmx,dnca,dbca,itwtype,nsyst,
c     inout
     $     ktwsmx,twsh,twfrsy,twfrph,sg_old,sg,dsg,wph,ag,wgt,
     $     gtotgr,crss,
c     out
     $     PRITW
     $     )
      implicit none
c
      integer,intent(in)::iph,nphmx,ngrmx,nmodmx,nsysmx,ntwsmx,
     $     ngr(0:nphmx),nslmod(nphmx),
     $     nslsys(nphmx),ntwmod(nphmx),ntwmmx,nsm(nmodmx,nphmx),
     $     itwtype(nsysmx,nphmx),NSYST(NPHmx)
      real*8,intent(in)::twsh(ntwmmx,nphmx),twthres(2,ntwmmx,nphmx),
     $     gamdot(nsysmx,ngrmx),dnca(3,nsysmx,nphmx),
     $     dbca(3,nsysmx,nphmx)
c     --
      integer,intent(inout)::ktwsmx(ngrmx)
      real*8,intent(inout)::twfrsy(ntwsmx,ngrmx),twfrph(ntwmmx,nphmx),
     $     sg_old(6,ngrmx),sg(6,ngrmx),dsg(6,ngrmx),wph(nphmx),
     $     ag(3,3,ngrmx),wgt(ngrmx),gtotgr(ngrmx),CRSS(NSYSMX,NGRMX)
      real*8,intent(out)::PRITW(NPHMX)
      real*8, intent(in)::tincr
      integer i,j,is,nsmx,ktm,itm,kts,its,kkkch,kkk,iphel,iskip,kptm
     $     ,kkkpa,ipts,itwtypex,ind0
      real*8 twshx,rand,twfmax,thres1,wgt_frac,twin_incr,
     $     bux6(6),aach(3,3),aapa(3,3),rstar(3,3),
     $     aux33(3,3),aux6(6),bux33(3,3),
     $     TWB(3),TWN(3),ACRYS(3,3),ATWIN(3,3)

c     This way of declaring IPTSGR, KPTMGR, TWFRGR was to allow save them
c     while complying with the INTEL fortran compiler
c     See (https://software.intel.com/content/www/us/en/develop/articles/fdiag6754.html)
c     for more details
      integer,allocatable,save :: IPTSGR(:),KPTMGR(:)
      real*8,allocatable,save:: TWFRGR(:,:)
c      SAVE  TWFRGR,IPTSGR,KPTMGR

      if (.not.allocated(iptsgr)) allocate(iptsgr(NGRmx))
      if (.not.allocated(kptmgr)) allocate(kptmgr(NGRmx))
      if (.not.allocated(twfrgr)) allocate(twfrgr(0:NTWMMX,NGRmx))

      ISKIP=1        ! controls temporary printout

C *** THIS LOOP UPDATES TWIN FRACTION FOR EACH SYSTEM (TWFRSY) , FOR EACH GRAIN
C      & MODE (TWFRGR) AND FOR EACH PHASE
C *** IDENTIFIES AND TAGS THE SYSTEM CONTRIBUTING MAXIMUM TWIN VOLUME FRACTION

      IPHEL=IPH
      DO 100 KKK=NGR(IPH-1)+1,NGR(IPH)

         IF(IPH.EQ.1) THEN !! parent
            KKKPA=KKK
            KKKCH=KKK+NGR(1)    ! related twin-child in phase=2
         ELSE IF(IPH.EQ.2) THEN
            KKKPA=KKK-NGR(1)
            KKKCH=KKK
         ENDIF

         IF(KTWSMX(KKK).EQ.0) THEN
c     index of predominant twin "system" in grain kkk
            IPTSGR(KKK)=0
c     index of predominant twin "mode" in grain kkk
            KPTMGR(KKK)=0
         ENDIF

CNT      TWF_THRES=1.E-10
         ITM=NSLMOD(IPHEL) !! slip system (to skip)
         ITS=NSLSYS(IPHEL) !! slip system (to skip)
         KTS=0

         DO KTM=1,NTWMOD(IPHEL)
            ITM=ITM+1
            TWFRGR(KTM,KKK)=0d0 !! pertain to each step
            TWSHX=TWSH(KTM,IPHEL)
            THRES1=TWTHRES(1,KTM,IPH) /NGR(1) !  CNT

            DO NSMX=1,NSM(ITM,IPHEL)
               KTS=KTS+1        ! shifted counter for twin systems
               ITS=ITS+1        ! absolute counter over all systems
               IF(GAMDOT(ITS,KKK) .GT. 0d0) THEN
c                           gamdot * dt / twinshear * f
                  TWIN_INCR=GAMDOT(ITS,KKK)*TINCR/TWSHX * WGT(KKK)
c     absolute fraction / cummulative
                  TWFRSY(KTS,KKK)=TWFRSY(KTS,KKK)+TWIN_INCR
c     absolute fraction / per step / per mode
                  TWFRGR(KTM,KKK)=TWFRGR(KTM,KKK)+TWIN_INCR

               ENDIF


               IF(KTWSMX(KKK).EQ.0) THEN
c     if not tagged for reorientation
                  RAND=1.0
ccc   RAND=0.5+0.5*RAN2(jran)  : adds stochasticity to selection
                  IF (TWFRSY(KTS,KKK)*rand .GT. THRES1) THEN
c     gets assigned here but reset to at next call if KWTSMX=0
                     IPTSGR(KKK)=ITS
                     KPTMGR(KKK)=KTM
                     TWFMAX=TWFRSY(KTS,KKK)

                  ENDIF
               ENDIF

c     end of DO NSMX loop over twin systems in mode
            ENDDO

C * total twin fraction in grain relative to polycrystal
c     absolute fraction / cummulative
            TWFRGR(0,KKK)=TWFRGR(0,KKK)+TWFRGR(KTM,KKK)
C * total twin fraction in aggregate for each twin mode and for each phase
            TWFRPH(KTM,IPH)=TWFRPH(KTM,IPH)+ TWFRGR(KTM,KKK)
c     end of DO KTM loop over twin modes in the phase
         ENDDO



C     TAGS THE GRAIN IF THE ACCUMULATED TWIN VOLUME OF ALL TWINS IN THE GRAIN
C     HAS EXCEEDED A THRESHOLD VALUE.
C     THRES1: MIN ACCUM TWIN FRACTION IN GRAIN BEFORE TWIN REOR IS SWITCHED ON.

         IF(KTWSMX(KKK).EQ.0 .AND. IPTSGR(KKK).NE.0) THEN
            IPTS =IPTSGR(KKK)
            KPTM=KPTMGR(KKK)
c     to compare absolute fractions
            THRES1=TWTHRES(1,KPTM,IPHEL) *WGT(KKK)
            IF (TWFRGR(0,KKK) .GT. THRES1) KTWSMX(KKK)=IPTS
         ENDIF


C *** IF KTWSMX>0 CREATES CHILD ORIENTATION & SWITCHES SIGN OF KTWSMX TO
C        FREEZE THIS AS THE VARIANT FOR FUTURE WEIGHT TRANSFER
C *** THIS IS DONE AFTER HAVING CALLED UPDATE_ORIENTATION IN MAIN,
C *** WHICH REORIENTS THE GRAIN BECAUSE OF SHEAR (INCLUDING TWIN SHEARS)

         IF (KTWSMX(KKK).GT.0 .AND. IPH.EQ.1) THEN
c     CREATES CHILD ONLY FOR PARENT PHASE
c     write(*,*)'kkk:',kkk
            IPTS=KTWSMX(KKK)
            ITWTYPEX=ITWTYPE(IPTS,IPHEL)
            DO I=1,3
               DO J=1,3
                  ACRYS(I,J)=AG(I,J,KKK)
               ENDDO
               TWN(I)=DNCA(I,IPTS,IPHEL)
               TWB(I)=DBCA(I,IPTS,IPHEL)
            ENDDO

CCC         CALL STAT_TWINNING (IPTS,KPTM,KKK,IPH)
            CALL TWIN_ORIENTATION (ITWTYPEX,TWB,TWN,ATWIN)

c           scenario 1. (inherits child's stress to child)
            aux33(:,:)=ag(:,:,kkkch)
            call transpose(aux33,3,aach)
c           scenario 2. (inherits parent's stress to child)
            aux33(:,:)=ag(:,:,kkkpa)
            call transpose(aux33,3,aapa)

c           assign new orientation to twin (child)
            call mat_dot(acrys,atwin,aux33,3)
            ag(:,:,kkkch)=aux33(:,:)

C *** ASSIGNS WEIGHT AND CRSS OF SL & TW SYSTEMS TO CREATED CHILD
C *** RESETS ACCUMULATED TWINNED FRACTIONS IN THE CHILD

c     absolute fraction accummulated so far
            WGT_FRAC  =TWFRGR(0,KKK)
            IF(WGT_FRAC .GT. WGT(KKK)) WGT_FRAC=WGT(KKK)
            WGT(KKKCH)= WGT_FRAC
            WGT(KKK)= WGT(KKK) - WGT_FRAC
            WPH(1)=WPH(1) - WGT_FRAC
            WPH(2)=WPH(2) + WGT_FRAC

CNT   TAGS THE SAME SYSTEM IN CHILD AS IN PARENT FOR ALLOWING DETWINNING
            KTWSMX(KKKCH)=KTWSMX(KKK)
            KPTMGR(KKKCH)=KPTMGR(KKK)

            GTOTGR(KKKCH)=GTOTGR(KKK)
            DO IS=1,NSLSYS(IPHEL)
c              assigns current CRSS to slip
               CRSS(IS,KKKCH)=CRSS(IS,KKK)
            ENDDO
            DO IS=NSLSYS(IPHEL)+1,NSYST(IPHEL)
c              to suppress twinning in child
               CRSS(IS,KKKCH)=5d0*CRSS(IS,KKK)
            ENDDO

c           to facilitate detwining
            CRSS(KTWSMX(KKKCH),KKKCH)= 0.5*CRSS(KTWSMX(KKK),KKK)
c           resets twin fraction in every mode of new child
            TWFRGR( : ,KKKCH)=0d0
c ***  initial stress state of created 'twins'.
            bux33(:,:)=ag(:,:,kkkch)
c           (inherits child's stress to child)
            call mat_dot(bux33,aach,rstar,3)
            ind0=kkkch
            aux6(:)=dsg(:,ind0)
            call mat6rot(aux6,rstar,bux6)
            dsg(:,kkkch)=bux6(:)
            aux6(:)=sg_old(:,ind0)
            call mat6rot(aux6,rstar,bux6)
            sg_old(:,kkkch)=bux6(:)
            aux6(:)=sg(:,ind0)
            call mat6rot(aux6,rstar,bux6)
            sg(:,kkkch)=bux6(:)
c        END OF IF KTWSMX(KKK)>0
         ENDIF


C *** TRANSFERS WEIGHT FROM PARENT TO CHILD & UPDATES TOTAL FRACTIONS.
C *** BECAUSE ONLY ONE VARIANT CONTRIBUTES TO CHILD GROWTH, AND WE ARE NOT
C     COUNTING THE OTHER TWIN SYSTEMS, THE TOTAL TWINNED VOLUME REPRESENTED
C     BY CHILDREN WILL UNDERESTIMATE THE TOTAL TWINNED VOLUME IN THE PXTAL
C *** SHOULD INTRODUCE A 'FACTOR' TO COMPENSATE FOR THIS !!!


         IF (KTWSMX(KKK) .LT. 0) THEN
c           this parent was tagged in a previous step

            KTS=ABS(KTWSMX(KKK)) - NSLSYS(IPH)
C             WGT_FRAC= TWFRSY(KTS,KKK)       ! overall accumulated fraction in system
c           incremental / all twins in mode
            WGT_FRAC= TWFRGR(KPTMGR(KKK),KKK)
            IF(WGT_FRAC .GT. WGT(KKK)) THEN
c              this grain will disappear (zero weight)
               WGT_FRAC=WGT(KKK)
               KTWSMX(KKK)=0
               IPTSGR(KKK)=0
               KPTMGR(KKK)=0
            ENDIF

            IF(IPH.EQ.1) THEN
c               write(*,*)'twinning?'
               WGT(KKK)=WGT(KKK) - WGT_FRAC
               WGT(KKKCH)=WGT(KKKCH) + WGT_FRAC
               WPH(1)=WPH(1) - WGT_FRAC
               WPH(2)=WPH(2) + WGT_FRAC
            ELSE IF(IPH.EQ.2) THEN ! CNT added
c               write(*,*)'detwinning?'
               WGT(KKK)=WGT(KKK) - WGT_FRAC
               WGT(KKKPA)=WGT(KKKPA) + WGT_FRAC
               WPH(1)=WPH(1) + WGT_FRAC
               WPH(2)=WPH(2) - WGT_FRAC
            ENDIF

            PRITW(IPH)=WPH(2)

         ENDIF

         IF(KTWSMX(KKK) .GT. 0) THEN
c           initial wgt of child was assigned at creation
c           Assign negative KTWSMX for recording thatthe grain has been reoriented.
            KTWSMX(KKK)= -KTWSMX(KKK)
         ENDIF

 100  ENDDO                     ! END OF DO KKK

      if (iphel.eq.1) then
         write( *,'('' -------> wph(1),wph(2),twfrph(1,1)'', 3f10.5)')
     $        wph(1),wph(2),twfrph(1,1)
      endif

      RETURN
      END subroutine update_twinning_vft_fe
