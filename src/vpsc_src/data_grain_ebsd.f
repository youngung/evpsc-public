C *****************************************************************************
C     SUBROUTINE DATA_GRAIN_EBSD        --->      VERSION 5/MAR/2024

C     READS GRAIN EULER ANGLES AND CALCULATES ROTATION MATRIX 'AG'
C     READS GRAIN ELLIPSOID EULER ANGLES AND SHAPE OF EACH PHASE (AND
C     OPTIONALLY OF EVERY GRAIN) AND CALCULATES DEFORMATION GRADIENT
C     'FIJPH' (AND OPTIONALLY 'FIJGR')
C *****************************************************************************

      SUBROUTINE DATA_GRAIN_EBSD(
c          in
     $     UR_SHP, ngrmx,nphmx,iph,igb,igt,
     $     wph,ishape,eulerph,axisph,coords,ncrds,ebsd_rad,
c          inout
     $     ngr_x,
c          out
     $     wgt,ag,fijph,fijgr
     $     )
      implicit none
      integer, intent(in)::UR_SHP,ngrmx,nphmx,iph,igb,igt,
     $     ishape(0:nphmx)
      integer, intent(inout)::ngr_x
      real*8,intent(out)::wgt(ngrmx),ag(3,3,1:ngrmx),
     $     fijph(3,3,0:nphmx),fijgr(3,3,1:ngrmx)
      real*8,intent(in)::wph(nphmx),eulerph(3,0:nphmx),
     $     axisph(0:3,3,0:nphmx)

      real*8 AA(3,3),FIJEA(3,3),FNEW(3,3),AX(3)
      CHARACTER NOMEN*1
      character prosa_x*120
      integer idum,i,j,kkk,m,k!,ngrain
      real*8 da,db,dc,eul1,eul2,eul3,totwgt
      logical verbose
c      parameter(verbose=.true.)
      parameter(verbose=.false.)

c----------------------------------------------------------------------c
c     for iodf=3 case
      real*8 ebsd_select_data(4,ngrmx)
      integer ncrds,nlines_ebsd,nrowmx
      real*8, intent(in):: coords(ncrds),ebsd_rad
      parameter(nrowmx=200000)
      real*4 ebsd_raw_data(5,nrowmx)
      common/kebsd/ebsd_raw_data,nlines_ebsd
c----------------------------------------------------------------------c

c      write(*,*)'Here data_grain_ebsd.f'
c      write(*,*)'nlines_ebsd:',nlines_ebsd
c      write(*,*)'ebsd_raw_data(1:5,1):',ebsd_raw_data(1:5,1)



      CALL select_points_from_ebsd(
c           in
     $           coords,ncrds,ebsd_rad,ngrmx,
c           out
     $           ebsd_select_data
     $       )


      if (verbose) call w_chrc(0,'inside data_grain_fe')

c      NGR(IPH)=NGR(IPH-1)+NGRAIN
c     NGR_x=NGRAIN

      NOMEN='B'

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
         WGT(KKK) = ebsd_select_data(4,KKK)
         if (verbose) call w_val(0,'wgt(kkk):',wgt(kkk))
         TOTWGT=TOTWGT+WGT(KKK)
         if (nomen.eq.'B'.or. nomen.eq.'b') then
            eul1= ebsd_select_data(1,KKK)
            eul2= ebsd_select_data(2,KKK)
            eul3= ebsd_select_data(3,KKK)
         else if(nomen.eq.'K' .or. nomen.eq.'k') then
            eul1= (ebsd_select_data(1,KKK)-90.)
            eul2= -ebsd_select_data(2,KKK)
            eul3= (-ebsd_select_data(1,KKK)-90.)
         else if(nomen.eq.'R' .or. nomen.eq.'r') then
            eul1= (ebsd_select_data(1,KKK)+90.)
            eul2= ebsd_select_data(2,KKK)
            eul3= (ebsd_select_data(3,KKK)-90.)
         else
            write(*,'(/,''CANNOT IDENTIFY'''//
     $           '''EULER ANGLE CONVENTION   !!'')')
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

      DO KKK=igb,igt!1,ngr_x
         WGT(KKK)=WGT(KKK)/TOTWGT*WPH(iph)
      ENDDO

      if (verbose) then
         call w_val(0,'wph:',wph(iph))
         do kkk=igb,igt!1,ngr_x
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
         DO KKK=igb,igt!1,ngr_x
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
        DO KKK=igb,igt!1,ngr_x
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
      END subroutine data_grain_ebsd
