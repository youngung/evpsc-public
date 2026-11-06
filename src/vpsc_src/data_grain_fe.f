C *****************************************************************************
C     SUBROUTINE DATA_GRAIN        --->      VERSION 8/FEB/2020
C
C     READS GRAIN EULER ANGLES AND CALCULATES ROTATION MATRIX 'AG'
C     READS GRAIN ELLIPSOID EULER ANGLES AND SHAPE OF EACH PHASE (AND
C     OPTIONALLY OF EVERY GRAIN) AND CALCULATES DEFORMATION GRADIENT
C     'FIJPH' (AND OPTIONALLY 'FIJGR')
C *****************************************************************************

      SUBROUTINE DATA_GRAIN_fe(
c          in
     $     UR_TX, UR_SHP, ngrmx,nphmx,iph,igb,igt,
     $     wph,ishape,eulerph,axisph,
c          inout
     $     ngr_x,
c          out
     $     wgt,ag,fijph,fijgr
     $     )
      implicit none
      integer, intent(in)::UR_TX,UR_SHP,ngrmx,nphmx,iph,igb,igt,
     $     ishape(0:nphmx)
      integer, intent(inout)::ngr_x
      real*8,intent(out)::wgt(ngrmx),ag(3,3,1:ngrmx),
     $     fijph(3,3,0:nphmx),fijgr(3,3,1:ngrmx)
      real*8,intent(in)::wph(nphmx),eulerph(3,0:nphmx),
     $     axisph(0:3,3,0:nphmx)

      real*8 AA(3,3),FIJEA(3,3),FNEW(3,3),EULANG(3),AX(3)
      CHARACTER NOMEN*1
      character prosa_x*120
      integer idum,i,j,kkk,m,k!,ngrain
      real*8 da,db,dc,eul1,eul2,eul3,totwgt
      logical verbose
c      parameter(verbose=.true.)
      parameter(verbose=.false.)

      if (verbose) call w_chrc(0,'inside data_grain_fe')

      READ(UR_TX,'(a)') PROSA_X
      READ(UR_TX,'(a)') PROSA_X
      READ(UR_TX,'(a)') PROSA_X
      READ(UR_TX,  *  ) NOMEN,idum !NGRAIN

c      NGR(IPH)=NGR(IPH-1)+NGRAIN
c      NGR_x=NGRAIN
      if (verbose) call w_ival(0,'ngr_x:',ngr_x)

      if(ngr_x.gt.ngrmx) then
         write(*,'('' number of grains exceeds dimension !!'')')
         write(*,'('' --> increase parameter NGRMX to'',i7)') ngr_x
         call xit()
      endif
C ***************************************************************************
C     READS EULER ANGLES, CONVERTS TO BUNGE NOTATION, CALCULATES ROT MATRIX
      TOTWGT=0.
      DO KKK=igb,igt
         READ(UR_TX,*) (EULANG(I),I=1,3),WGT(KKK)
         if (verbose) call w_val(0,'wgt(kkk):',wgt(kkk))
         TOTWGT=TOTWGT+WGT(KKK)
         if(nomen.eq.'B' .or. nomen.eq.'b') then
            eul1=  eulang(1)
            eul2=  eulang(2)
            eul3=  eulang(3)
         else if(nomen.eq.'K' .or. nomen.eq.'k') then
            eul1= (eulang(1)-90.)
            eul2= -eulang(2)
            eul3=(-eulang(3)-90.)
         else if(nomen.eq.'R' .or. nomen.eq.'r') then
            eul1= (eulang(1)+90.)
            eul2=  eulang(2)    ! fixed 29/apr/02
            eul3= (eulang(3)-90.)
         else
            write(*,'(/,''CANNOT IDENTIFY EULER ANGLE CONVENTION !!'')')
            call xit()
         endif

C ***********************************************************************
C *** OPTIONAL: RANDOM PERTURBATION OF EULER ANGLES BY +/- 'randrange'
c       if(KKK.EQ.NGR_x(IPH-1)+1)
c     #          print *, ' RANDOM SHIFT OF EULER ANGLES ACTIVATED  !!'
c        randrange= 2.
c       eul1= eul1+ (2.*ran2(jran)-1.)*randrange
c       eul2= eul2+ (2.*ran2(jran)-1.)*randrange
c       eul3= eul3+ (2.*ran2(jran)-1.)*randrange
C ***********************************************************************
C     CALCULATES THE TRANSFORMATION MATRIX AA WHICH TRANSFORMS FROM
C     SAMPLE TO CRYSTAL. STORES AG, WHICH TRANSFORMS FROM CRYSTAL TO SAMPLE.

        CALL EULER(2,EUL1,EUL2,EUL3,AA)
        DO J=1,3
        DO K=1,3
           AG(J,K,KKK)=AA(K,J)
        ENDDO
        ENDDO
      ENDDO     ! END OF LOOP OVER GRAINS IN EACH PHASE

C     DOUBLE RENORMALIZATION OF WEIGHTS: FIRST NORMALIZE THE WEIGHTS WITHIN
C     EACH PHASE, NEXT RENORMALIZE TO THE VOLUME FRACTION OF THE PHASE.

      DO KKK=igb,igt
         WGT(KKK)=WGT(KKK)/TOTWGT*WPH(iph)
      ENDDO

      if (verbose) then
         call w_val(0,'wph:',wph(iph))
         do kkk=igb,igt
            call w_val(0,'wgt(kkk):',wgt(kkk))
         enddo
      endif

C ***************************************************************************
C       INITIAL F TENSOR AND EIGENVECTORS:

C     * IF ISHAPE=0 ASSUMES SAME INITIAL SHAPE (axisph) AND ORIENTATION
C       (eulerph) FOR ALL THE GRAINS IN THE PHASE.
C       CALCULATES ESHELBY TENSOR FOR THE AVERAGE PHASE SHAPE.

C     * IF ISHAPE=1 ASSUMES SAME INITIAL SHAPE (axisph) AND ORIENTATION
C       (eulerph) FOR ALL THE GRAINS IN THE PHASE.
C       CALCULATES ESHELBY TENSOR FOR THE AVERAGE GRAIN PHASE.
C       CALCULATES INDIVIDUAL GRAIN SHAPE EVOLUTION.

C     * IF ISHAPE=2 ASSUMES SAME INITIAL SHAPE (axisph) AND ORIENTATION
C       (eulerph) FOR ALL THE GRAINS IN THE PHASE.
C       CALCULATES ESHELBY TENSOR FOR THE INDIVIDUAL GRAIN SHAPE AND
C       KEEPS TRACK OF INDIVIDUAL GRAIN SHAPE.

C     * IF ISHAPE=3 READS INDIVIDUAL INITIAL GRAIN AXES ORIENTATION
C       AND SHAPE FROM fileaxes (da,db,dc,ax(1),ax(2),ax(3)).
C       CALCULATES ESHELBY TENSOR FOR INDIVIDUAL GRAIN SHAPE AND
C       KEEPS TRACK OF INDIVIDUAL GRAIN SHAPE.

C       'eulerph' ANGLES OF (EA) WRT (SA).
C       'aa'    TRANSFORMS FROM (SA) TO (EA)
C       'fijea'  COLUMNS ARE GRAIN AXES EXPRESSED IN Ellipsoid system
C       'fnew'  COLUMNS ARE GRAIN AXES EXPRESSED IN SAMPLE SYSTEM

      IF(ISHAPE(iph).EQ.3) THEN

CFEB -----------------------------------------------------------------
C     WRITES MORPHOLOGIC TEXTURE FILE INTO 'RUN_LOG.OUT' FILE
         DO IDUM=1,10
            READ(UNIT=UR_SHP,END=99,FMT='(A)') PROSA_X
         ENDDO
 99      REWIND UR_SHP
CFEE -----------------------------------------------------------------

         READ(UR_SHP,'(A)') PROSA_X
         READ(UR_SHP,'(A)') PROSA_X
         READ(UR_SHP,'(A)') PROSA_X
         READ(UR_SHP,'(A)') PROSA_X
C *** 'EULERPH' AND 'AXISPH(0') ARE READ in VPSC_INPUT from VPSC.IN.
C ***  THEY MAY BE RE RE-READ FROM THE MORPHOLOGIC TEXTURE FILE IF ishape=3.
C         READ(UR_SHP,*) (AXISPH(0,I,IPH),I=1,3)
C         READ(UR_SHP,*) (EULERPH(I,IPH) ,I=1,3)
      ENDIF


      da=EULERPH(1,iph)
      db=EULERPH(2,iph)
      dc=EULERPH(3,iph)

      call euler(2,da,db,dc,aa)

      do i=1,3
      do j=1,3
         fijea(i,j)=(i/j)*(j/i)*AXISPH(0,I,iph)
      enddo
      enddo

      do j=1,3
      do i=1,3
         fnew(i,j)=0.
      do m=1,3
         fnew(i,j)=fnew(i,j)+aa(m,i)*fijea(m,j)
      enddo
      enddo
      enddo

      do i=1,3
      do j=1,3
         fijph(i,j,iph)=fnew(i,j)
      enddo
      enddo
CFEB
      IF(ISHAPE(iph).EQ.1 .OR. ISHAPE(iph).EQ.2) THEN
         DO KKK=igb,igt
         DO I=1,3
         DO J=1,3
            FIJGR(I,J,KKK)=FNEW(I,J)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

C *** ISHAPE=3 -> READ EULER ANGLES & ASPECT RATIOS OF EACH ELLIPSOID.
C *** DA,DB,DC -> EULER ANGLES OF ELLIPSOID AXES WRT SAMPLE.

      IF(ISHAPE(iph).EQ.3) THEN
        DO KKK=igb,igt
           if(ishape(iph).eq.3) then
              READ(UR_SHP,*) da,db,dc,ax(1),ax(2),ax(3)
              do i=1,3
              do j=1,3
                 fijgr(i,j,KKK)=(i/j)*(j/i)*ax(i)
              enddo
              enddo
           endif

           call euler(2,da,db,dc,aa)

C     AA TRANSFORMS FROM SAMPLE TO GRAIN AXES. AAt FROM GRAIN TO SAMPLE.

          do j=1,3
          do i=1,3
             fnew(i,j)=0.
             do m=1,3
                fnew(i,j)=fnew(i,j)+aa(m,i)*fijgr(m,j,kkk)
             enddo
          enddo
          enddo

          do i=1,3
          do j=1,3
             fijgr(i,j,KKK)=fnew(i,j)
          enddo
          enddo
        ENDDO      ! END OF LOOP OVER GRAINS
      ENDIF      ! END OF CONDITION FOR ISHAPE=3
CFEE
      RETURN
      END subroutine data_grain_fe



C *****************************************************************************
C     SUBROUTINE DATA_GRAIN_FE_SHUFFLE
C
C     READS GRAIN EULER ANGLES AND CALCULATES ROTATION MATRIX 'AG'
C     READS GRAIN ELLIPSOID EULER ANGLES AND SHAPE OF EACH PHASE (AND
C     OPTIONALLY OF EVERY GRAIN) AND CALCULATES DEFORMATION GRADIENT
C     'FIJPH' (AND OPTIONALLY 'FIJGR')
C *****************************************************************************

      SUBROUTINE DATA_GRAIN_fe_shuffle(
c          in
     $     UR_TX, UR_SHP, ngrmx,nphmx,iph,igb,igt,
     $     wph,ishape,eulerph,axisph,
c
     $     no_grain_ori,
c          out
     $     wgt,ag,fijph,fijgr
     $     )
      implicit none
      integer, intent(in)::no_grain_ori,
     $     UR_TX,UR_SHP,ngrmx,nphmx,iph,igb,igt,
     $     ishape(0:nphmx)
      real*8,intent(out)::wgt(ngrmx),ag(3,3,1:ngrmx),
     $     fijph(3,3,0:nphmx),fijgr(3,3,1:ngrmx)
      real*8,intent(in)::wph(nphmx),eulerph(3,0:nphmx),
     $     axisph(0:3,3,0:nphmx)

c     locals
      real*8 AA(3,3),FIJEA(3,3),FNEW(3,3),EULANG(3),AX(3),
     $     da,db,dc,eul1,eul2,eul3,totwgt,texarray(no_grain_ori,4)
      CHARACTER NOMEN*1
      character prosa_x*120
      integer idum,i,j,kkk,m,k,iarr(ngrmx)
      logical verbose
c      parameter(verbose=.true.)
      parameter(verbose=.false.)

      if (verbose) call w_chrc(0,'inside data_grain_fe')

      READ(UR_TX,'(a)') PROSA_X
      READ(UR_TX,'(a)') PROSA_X
      READ(UR_TX,'(a)') PROSA_X
      READ(UR_TX,  *  ) NOMEN,idum

      if (idum.ne.no_grain_ori) then
         write(*,*)'** Error: no_grain_ori should equal to idum'
         call xit()
      endif

c      NGR(IPH)=NGR(IPH-1)+NGRAIN

C ***************************************************************************
C     READS EULER ANGLES, CONVERTS TO BUNGE NOTATION, CALCULATES ROT MATRIX

      do kkk=1,no_grain_ori
         READ(UR_TX,*) (texarray(kkk,I),I=1,4)
      enddo

      call knuth_shuffle(no_grain_ori,ngrmx,iarr)

      TOTWGT=0.
      do i=igb,igt
         kkk=iarr(i)
         wgt(i)=texarray(kkk,4) !WGT(KKK)
         TOTWGT=TOTWGT+wgt(i)
         eulang(:)=texarray(kkk,1:3)
         if(nomen.eq.'B' .or. nomen.eq.'b') then
            eul1=  eulang(1)
            eul2=  eulang(2)
            eul3=  eulang(3)
         else if(nomen.eq.'K' .or. nomen.eq.'k') then
            eul1= (eulang(1)-90.)
            eul2= -eulang(2)
            eul3=(-eulang(3)-90.)
         else if(nomen.eq.'R' .or. nomen.eq.'r') then
            eul1= (eulang(1)+90.)
            eul2=  eulang(2)    ! fixed 29/apr/02
            eul3= (eulang(3)-90.)
         else
            write(*,'(/,''CANNOT IDENTIFY EULER ANGLE CONVENTION !!'')')
            call xit()
         endif

C ***********************************************************************
C *** OPTIONAL: RANDOM PERTURBATION OF EULER ANGLES BY +/- 'randrange'
c       if(KKK.EQ.NGR(IPH-1)+1)
c     #          print *, ' RANDOM SHIFT OF EULER ANGLES ACTIVATED  !!'
c        randrange= 2.
c       eul1= eul1+ (2.*ran2(jran)-1.)*randrange
c       eul2= eul2+ (2.*ran2(jran)-1.)*randrange
c       eul3= eul3+ (2.*ran2(jran)-1.)*randrange
C ***********************************************************************
C     CALCULATES THE TRANSFORMATION MATRIX AA WHICH TRANSFORMS FROM
C     SAMPLE TO CRYSTAL. STORES AG, WHICH TRANSFORMS FROM CRYSTAL TO SAMPLE.

        CALL EULER(2,EUL1,EUL2,EUL3,AA)
        DO J=1,3
        DO K=1,3
           AG(J,K,i)=AA(K,J)
        ENDDO
        ENDDO
      ENDDO     ! END OF LOOP OVER GRAINS IN EACH PHASE

C     DOUBLE RENORMALIZATION OF WEIGHTS: FIRST NORMALIZE THE WEIGHTS WITHIN
C     EACH PHASE, NEXT RENORMALIZE TO THE VOLUME FRACTION OF THE PHASE.

      DO KKK=igb,igt
         WGT(KKK)=WGT(KKK)/TOTWGT*WPH(iph)
      ENDDO

      if (verbose) then
         call w_val(0,'wph:',wph(iph))
         do kkk=igb,igt
            call w_val(0,'wgt(kkk):',wgt(kkk))
         enddo
      endif

C ***************************************************************************
C       INITIAL F TENSOR AND EIGENVECTORS:

C     * IF ISHAPE=0 ASSUMES SAME INITIAL SHAPE (axisph) AND ORIENTATION
C       (eulerph) FOR ALL THE GRAINS IN THE PHASE.
C       CALCULATES ESHELBY TENSOR FOR THE AVERAGE PHASE SHAPE.

C     * IF ISHAPE=1 ASSUMES SAME INITIAL SHAPE (axisph) AND ORIENTATION
C       (eulerph) FOR ALL THE GRAINS IN THE PHASE.
C       CALCULATES ESHELBY TENSOR FOR THE AVERAGE GRAIN PHASE.
C       CALCULATES INDIVIDUAL GRAIN SHAPE EVOLUTION.

C     * IF ISHAPE=2 ASSUMES SAME INITIAL SHAPE (axisph) AND ORIENTATION
C       (eulerph) FOR ALL THE GRAINS IN THE PHASE.
C       CALCULATES ESHELBY TENSOR FOR THE INDIVIDUAL GRAIN SHAPE AND
C       KEEPS TRACK OF INDIVIDUAL GRAIN SHAPE.

C     * IF ISHAPE=3 READS INDIVIDUAL INITIAL GRAIN AXES ORIENTATION
C       AND SHAPE FROM fileaxes (da,db,dc,ax(1),ax(2),ax(3)).
C       CALCULATES ESHELBY TENSOR FOR INDIVIDUAL GRAIN SHAPE AND
C       KEEPS TRACK OF INDIVIDUAL GRAIN SHAPE.

C       'eulerph' ANGLES OF (EA) WRT (SA).
C       'aa'    TRANSFORMS FROM (SA) TO (EA)
C       'fijea'  COLUMNS ARE GRAIN AXES EXPRESSED IN Ellipsoid system
C       'fnew'  COLUMNS ARE GRAIN AXES EXPRESSED IN SAMPLE SYSTEM

      IF(ISHAPE(iph).EQ.3) THEN

CFEB -----------------------------------------------------------------
C     WRITES MORPHOLOGIC TEXTURE FILE INTO 'RUN_LOG.OUT' FILE
         DO IDUM=1,10
            READ(UNIT=UR_SHP,END=99,FMT='(A)') PROSA_X
         ENDDO
 99      REWIND UR_SHP
CFEE -----------------------------------------------------------------
         READ(UR_SHP,'(A)') PROSA_X
         READ(UR_SHP,'(A)') PROSA_X
         READ(UR_SHP,'(A)') PROSA_X
         READ(UR_SHP,'(A)') PROSA_X
      ENDIF


      da=EULERPH(1,iph)
      db=EULERPH(2,iph)
      dc=EULERPH(3,iph)

      call euler(2,da,db,dc,aa)

      do i=1,3
      do j=1,3
         fijea(i,j)=(i/j)*(j/i)*AXISPH(0,I,iph)
      enddo
      enddo

      do j=1,3
      do i=1,3
         fnew(i,j)=0.
      do m=1,3
         fnew(i,j)=fnew(i,j)+aa(m,i)*fijea(m,j)
      enddo
      enddo
      enddo

      do i=1,3
      do j=1,3
         fijph(i,j,iph)=fnew(i,j)
      enddo
      enddo
CFEB
      IF(ISHAPE(iph).EQ.1 .OR. ISHAPE(iph).EQ.2) THEN
         DO KKK=igb,igt
         DO I=1,3
         DO J=1,3
            FIJGR(I,J,KKK)=FNEW(I,J)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

C *** ISHAPE=3 -> READ EULER ANGLES & ASPECT RATIOS OF EACH ELLIPSOID.
C *** DA,DB,DC -> EULER ANGLES OF ELLIPSOID AXES WRT SAMPLE.

      IF(ISHAPE(iph).EQ.3) THEN
        DO KKK=igb,igt
           if(ishape(iph).eq.3) then
              READ(UR_SHP,*) da,db,dc,ax(1),ax(2),ax(3)
              do i=1,3
              do j=1,3
                 fijgr(i,j,KKK)=(i/j)*(j/i)*ax(i)
              enddo
              enddo
           endif

           call euler(2,da,db,dc,aa)

C     AA TRANSFORMS FROM SAMPLE TO GRAIN AXES. AAt FROM GRAIN TO SAMPLE.

          do j=1,3
          do i=1,3
             fnew(i,j)=0.
             do m=1,3
                fnew(i,j)=fnew(i,j)+aa(m,i)*fijgr(m,j,kkk)
             enddo
          enddo
          enddo

          do i=1,3
          do j=1,3
             fijgr(i,j,KKK)=fnew(i,j)
          enddo
          enddo
        ENDDO      ! END OF LOOP OVER GRAINS
      ENDIF      ! END OF CONDITION FOR ISHAPE=3
CFEE
      RETURN
      END subroutine data_grain_fe_shuffle
