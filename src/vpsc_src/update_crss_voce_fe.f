C *****************************************************************************
C     SUBROUTINE UPDATE_CRSS_VOCE_FE    --->    VERSION OF 21/NOV/2010 adpeted to dEVPSC
C
C     A VOCE LAW, FUNCTION OF THE ACCUMULATED SHEAR IN EACH GRAIN
C     MODULATES THE LATENT (LINEAR) HARDENING.
C     THE UNITS (AND THE STRENGTH) OF THE HARDENING ARE CARRIED BY THE
C     MULTIPLICATIVE FACTOR 'VOCE'.
C     THE SELF & LATENT COUPLING COEFFICIENTS 'HARD' ARE DIMENSIONLESS
C     CONSTANTS RELATIVE TO THE FACTOR 'VOCE'.
C     THE INCREMENT OF CRSS IS GIVEN BY AN ANALYTIC INTEGRAL EXPRESSION OF
C     THE VOCE FUNCTION, RATHER THAN USING FORWARD EXTRAPOLATION, AS WAS
C     DONE BEFORE 28/SET/00.
C     THE PARAMETERS IN VOCE EXPRESSION MAY ADOPT 'NON-KOSHER' VALUES (DEC/00)
C
C     IOPTION=0: READS HARDENING PARAMETERS OF SLIP AND TWIN MODES AND
C                ASSIGNS THOSE TO THE SLIP AND TWIN SYSTEMS
C     IOPTION=1: INITIALIZES THE CRSS OF EACH SYSTEM
C     IOPTION=2: UPDATES THE CRSS OF EACH SYSTEM AS A FUNCTION OF THE SHEAR
C                INCREMENTS TAKING PLACE IN ALL SYSTEMS
*******************************************************************************

      SUBROUTINE UPDATE_CRSS_VOCE_fe(
c          in
     $     UR_SX,IOPTION,iphel,nphmx,ngrmx,dtime,nmodmx,ntwmmx,igb,igt,
     $     nmodes,nsysmx,nsm,nslsys,nsyst,gamdot_x,
c          inout
     $     grsze,gtotgr_x,tau_x,thet_x,hpfac_x,hard_x,crss_x,
c          out
     $     twthres_x,nrs_x,itwinlaw_x,isectw_x
     $     )
      implicit none
      integer,intent(in)::igb,igt,UR_SX,nphmx,nmodmx,ntwmmx,
     $     nmodes(nphmx),
     $     nsysmx,nsm(nmodmx,nphmx),nslsys(nphmx),nsyst(nphmx),ngrmx
      real*8,intent(in)::dtime,gamdot_x(nsysmx,ngrmx)
      real*8,intent(inout)::grsze(nphmx),TAU_x(nsysmx,0:1,nphmx),
     $     THET_x(nsysmx,0:1,nphmx),HPFAC_x(1:nsysmx,1:nphmx),
     $     hard_x(nsysmx,nsysmx,nphmx),gtotgr_x(ngrmx),
     $     crss_x(nsysmx,ngrmx)
      real*8,intent(out)::twthres_x(2,ntwmmx,nphmx)
      integer, intent(out)::nrs_x(nsysmx,nphmx),itwinlaw_x,
     $     isectw_x(nsysmx,nphmx)

c$$$      real*8 TAU(NSYSMX_x,0:1,NPHmx),THET(NSYSMX_x,0:1,NPHmx),
c$$$     $     HPFAC(NSYSMX_x,NPHmx)
c$$$      COMMON/VOCE/TAU,THET,HPFAC
c     --

      real*8 HSELFX(NMODMX),HLATEX(NMODMX,NMODMX)

c     local
      integer kkk,is,js,iprx,iprint,ioption,ntwmodx,isys,iphel,
     $     jm,jsys,isectwx,nrsx,im
      real*8 tiny,voce,fact,tau1,thet0,thet1,tau0_,dtau,expini,
     $     deltgam,gamtotx,expdel,thres1x,thres2x,hpfacx,tau1x,
     $     thet0x,thet1x,tau0x
      logical verbose
c      parameter(verbose=.true.)
      parameter(verbose=.false.)

      IPRINT=0    ! controls printing

c      write(*,*)'ioption given to <update_crss_voce_fe>:',ioption


C ***************************************************************************

      IF (IOPTION.EQ.0) THEN
        READ(UR_SX,*) GRSZE(iphel)

        ISYS=1
        NTWMODX=0
        DO IM=1,NMODES(iphel)

           READ(UR_SX,*)          ! reminder/identifier of system
           READ(UR_SX,*) NRSX     ! rate sensitivity parameter

C *** READ EXTENDED VOCE PARAMETERS AND LATENT HARDENING PARAMETERS.
C *** SELF-HARDENING IS ARBITRARILY SET AS REFERENCE
           READ(UR_SX,*) TAU0X,TAU1X,THET0X,THET1X, HPFACX
           HSELFX(IM)=1d0
           READ(UR_SX,*) (HLATEX(IM,JM),JM=1,NMODES(iphel))

C *** CHECKS WHETHER VOCE PARAMETERS ARE KOSHER:
C         TAU0>0 , TAU1 >= 0 , THET0 >= THET1 >= 0
C         TAU1=0   CORRESPONDS TO LINEAR HARDENING.
C         THETA0=0 FORCES NO-HARDENING.
C *** IF VOCE PARAMETERS ARE NON-KOSHER CHECKS FOR ILL-POSED HARDENING.

           CALL CHECK_VOCE (IM,IPHel,TAU0X,TAU1X,THET0X,THET1X)

C *** READS PTR PARAMETERS FOR TWINNING SYSTEMS
           IF(ISYS.LE.NSLSYS(iphel)) ISECTWX=0
           IF(ISYS.GT.NSLSYS(iphel)) THEN
              READ(UR_SX,*) ISECTWX,ITWINLAW_x,THRES1X,THRES2X
              NTWMODX=NTWMODX+1
              TWTHRES_x(1,NTWMODX,iphel)=THRES1X
              TWTHRES_x(2,NTWMODX,iphel)=THRES2X
           ENDIF

           DO IS=1,NSM(IM,iphel)
              NRS_x(ISYS,iphel)   =NRSX
              TAU_x(ISYS,0,iphel) =TAU0X
              TAU_x(ISYS,1,iphel) =TAU1X
              THET_x(ISYS,0,iphel)=THET0X
              THET_x(ISYS,1,iphel)=THET1X

              ISECTW_x(ISYS,iphel)=ISECTWX !
              HPFAC_x(ISYS,iphel) =HPFACX

              JSYS=1
              DO JM=1,NMODES(iphel)
              DO JS=1,NSM(JM,iphel)
                 HARD_x(ISYS,JSYS,iphel)=HLATEX(IM,JM)
                 JSYS=JSYS+1
              ENDDO
              ENDDO
              HARD_x(ISYS,ISYS,iphel)=HSELFX(IM)
              ISYS=ISYS+1

           ENDDO                ! end of DO IS
        ENDDO                   ! end of DO IM
      ENDIF                     ! END OF IOPTION=0
C ***************************************************************************

      IF (IOPTION.EQ.1) THEN
         if(verbose)then
            call w_ival(0,'iphel:',iphel)
         endif

         if (verbose) then
            call w_val(0,'tau_x(1,0,iphel):',tau_x(1,0,iphel))
         endif


         DO KKK=igb,igt!NGR(IPHel-1)+1,NGR(IPHel)
            DO IS=1,NSYST(iphel)
               CRSS_x(IS,KKK)=TAU_x(IS,0,iphel)
     $              +HPFAC_x(IS,iphel)/SQRT(GRSZE(iphel))
            ENDDO
         ENDDO
      ENDIF   ! END OF IOPTION=1
C ***************************************************************************

      IF (IOPTION.EQ.2) THEN

         if (verbose)then
            call w_chrc(0,'** Update crss voce with ioption=2')
            call w_val(0,'dtime:',dtime)
         endif

        DO KKK=igb,igt!NGR(IPHel-1)+1,NGR(IPHel)
           GAMTOTX=GTOTGR_x(KKK)
           DELTGAM=0d0
           DO IS=1,NSYST(iphel)
              DELTGAM=DELTGAM+ABS(GAMDOT_x(IS,Kkk))*dtime
           ENDDO

c          call w_val(0,'deltgam:',deltgam)

          DO IS=1,NSYST(iphel)
            DTAU=0.
            DO JS=1,NSYST(iphel)
              DTAU=DTAU+HARD_x(IS,JS,iphel)*ABS(GAMDOT_x(JS,Kkk))*dtime
            ENDDO
            TAU0_ =TAU_x(IS,0,iphel)
            TAU1 =TAU_x(IS,1,iphel)
            THET0=THET_x(IS,0,iphel)
            THET1=THET_x(IS,1,iphel)
            TINY=1.E-4*TAU0_

            VOCE=0.0
            IF(ABS(THET0).GT.TINY) THEN
              VOCE=THET1*DELTGAM
              IF(ABS(TAU1).GT.TINY) THEN
                FACT=ABS(THET0/TAU1)
                EXPINI=EXP(-GAMTOTX*FACT)
                EXPDEL=EXP(-DELTGAM*FACT)
                VOCE  =VOCE-(FACT*TAU1-THET1)/FACT*EXPINI*
     #            (EXPDEL-1.)-THET1/FACT*EXPINI*
     #            (EXPDEL*((GAMTOTX+DELTGAM)*FACT+1.)-(GAMTOTX*FACT+1.))
              ENDIF
            ENDIF
            if (deltgam.gt.0)then
               CRSS_x(IS,KKK)=CRSS_x(IS,KKK)+DTAU*VOCE/DELTGAM
            endif
            if (isnan(crss_x(is,kkk))) then
               write(*,*)'crss is NaN dtau*voce/deltgam'
               write(*,*)'dtau:',dtau
               write(*,*)'voce:',voce
               write(*,*)'deltgam:',deltgam
               call xit()
            endif
          ENDDO
          GTOTGR_x(KKK)=GAMTOTX+DELTGAM
       ENDDO

      ENDIF   ! END OF IOPTION=2
C ***************************************************************************

      RETURN
      END subroutine update_crss_voce_fe


C *********************************************************************
C     SUBROUTINE CHECK_VOCE      --->      VERSION 12/JUN/2001
C *********************************************************************
C
C     CHECKS WHETHER VOCE PARAMETERS ARE KOSHER:
C        TAU0>0 , TAU1 >= 0 , THET0 >= THET1 >= 0
C        TAU1=0   CORRESPONDS TO LINEAR HARDENING.
C        THETA0=0 FORCES NO-HARDENING.
C     IF VOCE PARAMETERS ARE NON-KOSHER CHECKS FOR ILL-POSED HARDENING.
C *********************************************************************

      SUBROUTINE CHECK_VOCE (KOUNT,IPHel,TAU0X,TAU1X,THET0X,THET1X)
      implicit none
      real*8 thet1x,thet0x,tiny,tau1x,tau0x
      integer kount,iphel

   70 FORMAT(' *** MODE',I3,'   IN PHASE',I3)
   71 FORMAT('     TAU0.LE.0 --> ILL-POSED HARDENING  !!!')
c   72 FORMAT('     THETA1<0 --> NON-KOSHER HARDENING MAY GIVE TAU<0')
c   73 FORMAT('     TAU1<0  --> NON-KOSHER HARDENING MAY GIVE TAU<0')
   75 FORMAT('     THETA0=0 --> WILL RESET TAU1=THETA1=0')
c   76 FORMAT('     |THETA1|.GE.|THETA0| --> NON-KOSHER HARDENING')
   77 FORMAT('     THETA0<0 --> HARDENING DEPENDS ON |THETA0| ANYWAY')

      TINY=1.E-4*TAU0X
      IF(TAU0X.LE.0.0) THEN
        WRITE(*,70) KOUNT,IPHel
        WRITE(*,71)
        write(*,*) 'Err: TAU0X<0'
        call exit(-1)

      ENDIF
      IF(ABS(THET1X).LE.TINY) THEN
        THET1X=0.
      ELSE IF(THET1X.LT.0.0) THEN
c         write(*,*)
c        WRITE(*,70) KOUNT,IPHel
c        WRITE(*,72)
!        print *, 'enter c to continue'
!        read  *
      ENDIF
      IF(TAU1X.LT.0) THEN
c        WRITE(*,70) KOUNT,IPHel
c        WRITE(*,73)
!        print *, 'enter c to continue'
!        read  *
      ENDIF
      IF(ABS(THET0X).LE.TINY) THEN
        IF(ABS(TAU1X).LE.TINY) THEN
          TAU1X =0.
          THET0X=THET1X
        ENDIF
        IF(ABS(TAU1X).GT.TINY) THEN
          WRITE(*,70) KOUNT,IPHel
          WRITE(*,75)
!          print *, 'enter c to continue'
!          read  *
          TAU1X =0.
          THET0X=0.
          THET1X=0.
        ENDIF
      ENDIF
      IF(ABS(TAU1X).LE.TINY) THEN
        TAU1X =0.0
        THET0X=THET1X
      ENDIF
      IF(THET0X.LT.0.0) THEN
        WRITE(*,70) KOUNT,IPHel
        WRITE(*,77)
!        print *, 'enter c to continue'
!        read  *
        THET0X=ABS(THET0X)
      ENDIF
      IF(TAU1X.NE.0.) THEN
        IF(ABS(THET1X).GE.ABS(THET0X)) THEN
c          WRITE(*,70) KOUNT,IPHel
c          WRITE(*,76)
!          print *, 'enter c to continue'
!          read  *
        ENDIF
      ENDIF

      RETURN
      END































C *****************************************************************************
C     SUBROUTINE UPDATE_CRSS_VOCE_thermo_fe
C
c     Same as update_crss_voce_fe except that temperature effect is accounted for
c     by multiplying an exponential factor.

c     GRP edit to include temperature-dependent Voce Parameters 25 Apr 22
*******************************************************************************

      SUBROUTINE UPDATE_CRSS_VOCE_thermo_fe(
c          in
     $     UR_SX,IOPTION,iphel,nphmx,ngrmx,dtime,nmodmx,ntwmmx,igb,igt,
     $     nmodes,nsysmx,nsm,nslsys,nsyst,gamdot_x,temp,tempdot,
c          inout
     $     grsze,gtotgr_x,tau_x,thet_x,hpfac_x,hard_x,crss_x,
c          out
     $     twthres_x,nrs_x,itwinlaw_x,isectw_x
     $     )
      implicit none
      integer,intent(in)::igb,igt,UR_SX,nphmx,nmodmx,ntwmmx,
     $     nmodes(nphmx),
     $     nsysmx,nsm(nmodmx,nphmx),nslsys(nphmx),nsyst(nphmx),ngrmx
      real*8,intent(in)::dtime,gamdot_x(nsysmx,ngrmx),temp,tempdot
      real*8,intent(inout)::grsze(nphmx),TAU_x(nsysmx,0:1,nphmx),
     $     THET_x(nsysmx,0:1,nphmx),HPFAC_x(1:nsysmx,1:nphmx),
     $     hard_x(nsysmx,nsysmx,nphmx),gtotgr_x(ngrmx),
     $     crss_x(nsysmx,ngrmx)
      real*8,intent(out)::twthres_x(2,ntwmmx,nphmx)
      integer, intent(out)::nrs_x(nsysmx,nphmx),itwinlaw_x,
     $     isectw_x(nsysmx,nphmx)

c$$$      real*8 TAU(NSYSMX_x,0:1,NPHmx),THET(NSYSMX_x,0:1,NPHmx),
c$$$     $     HPFAC(NSYSMX_x,NPHmx)
c$$$      COMMON/VOCE/TAU,THET,HPFAC
c     --

      real*8 HSELFX(NMODMX),HLATEX(NMODMX,NMODMX),dtaue(nsysmx,ngrmx)
c     local
      integer kkk,is,js,ks,iprx,iprint,ioption,ntwmodx,isys,iphel,
     $     jm,jsys,isectwx,nrsx,im
      real*8 tiny,voce,fact,tau1,thet0,thet1,tau0_,dtau,expini,
     $     deltgam,gamtotx,expdel,thres1x,thres2x,hpfacx,tau1x,
     $     thet0x,thet1x,tau0x,beta(nmodmx),tref
      logical verbose
c      parameter(verbose=.true.)
      parameter(verbose=.false.)

c              for uranium case.

c      beta(:)=0d0
      ! there are only two deformation modes in 'niobium' at the moment.
c      beta(1)=250.
c      beta(2)=230.
      ! beta(3)=160.
      ! beta(4)=1000.
      tref=295.

      IPRINT=0    ! controls printing
c      write(*,*)'ioption given to <update_crss_voce_fe>:',ioption

C ***************************************************************************

      IF (IOPTION.EQ.0) THEN

        READ(UR_SX,*) GRSZE(iphel)

        ISYS=1
        NTWMODX=0
        DO IM=1,NMODES(iphel)

           READ(UR_SX,*)          ! reminder/identifier of system
           READ(UR_SX,*) NRSX     ! rate sensitivity parameter

C *** READ EXTENDED VOCE PARAMETERS AND LATENT HARDENING PARAMETERS.
C *** SELF-HARDENING IS ARBITRARILY SET AS REFERENCE
           READ(UR_SX,*) TAU0X,TAU1X,THET0X,THET1X, HPFACX


           HSELFX(IM)=1d0
           READ(UR_SX,*) (HLATEX(IM,JM),JM=1,NMODES(iphel))

C *** CHECKS WHETHER VOCE PARAMETERS ARE KOSHER:
C         TAU0>0 , TAU1 >= 0 , THET0 >= THET1 >= 0
C         TAU1=0   CORRESPONDS TO LINEAR HARDENING.
C         THETA0=0 FORCES NO-HARDENING.
C *** IF VOCE PARAMETERS ARE NON-KOSHER CHECKS FOR ILL-POSED HARDENING.

           CALL CHECK_VOCE (IM,IPHel,TAU0X,TAU1X,THET0X,THET1X)

C *** READS PTR PARAMETERS FOR TWINNING SYSTEMS
           IF(ISYS.LE.NSLSYS(iphel)) ISECTWX=0
           IF(ISYS.GT.NSLSYS(iphel)) THEN
              READ(UR_SX,*) ISECTWX,ITWINLAW_x,THRES1X,THRES2X
              NTWMODX=NTWMODX+1
              TWTHRES_x(1,NTWMODX,iphel)=THRES1X
              TWTHRES_x(2,NTWMODX,iphel)=THRES2X
           ENDIF

           DO IS=1,NSM(IM,iphel)
              NRS_x(ISYS,iphel)   =NRSX
              TAU_x(ISYS,0,iphel) =TAU0X
              TAU_x(ISYS,1,iphel) =TAU1X
              THET_x(ISYS,0,iphel)=THET0X
              THET_x(ISYS,1,iphel)=THET1X

              ISECTW_x(ISYS,iphel)=ISECTWX !
              HPFAC_x(ISYS,iphel) =HPFACX

              JSYS=1
              DO JM=1,NMODES(iphel)
              DO JS=1,NSM(JM,iphel)
                 HARD_x(ISYS,JSYS,iphel)=HLATEX(IM,JM)
                 JSYS=JSYS+1
              ENDDO
              ENDDO
              HARD_x(ISYS,ISYS,iphel)=HSELFX(IM)
              ISYS=ISYS+1

           ENDDO                ! end of DO IS
        ENDDO                   ! end of DO IM
      ENDIF                     ! END OF IOPTION=0
C ***************************************************************************

      IF (IOPTION.EQ.1) THEN

c     ** Let's initialize the CRSS here, which applies to all slip/twin systems.
        ISYS=1
        DO IM=1,NMODES(iphel)
C *** Allow the code to read in the parameters to minimize editting, but unnecessary...still need HPFACX
C *** temp in K, Voce Parameters for Nb
           TAU0X =-32.91*LOG(temp) + 247.74
           TAU1X =-63.35*LOG(temp) + 475.4
c *** need to implement high temperature mechanism change for Nb
           THET0X=-202.9*LOG(temp) + 1591
           THET1X=-7.507*LOG(temp) + 54.851
           DO IS=1,NSM(IM,iphel)
              NRS_x(ISYS,iphel)   =NRSX
              TAU_x(ISYS,0,iphel) =TAU0X
              TAU_x(ISYS,1,iphel) =TAU1X
              THET_x(ISYS,0,iphel)=THET0X
              THET_x(ISYS,1,iphel)=THET1X
           ENDDO                ! end of DO IS
        ENDDO                   ! end of DO IM

         if(verbose)then
            call w_ival(0,'iphel:',iphel)
         endif
c         if (verbose) then
            call w_val(0,'tau_x(1,0,iphel):',tau_x(1,0,iphel))
c         endif
         DO KKK=igb,igt
            is=0
            do im=1,nmodes(iphel)
               do js=1,nsm(im,iphel)
                  is=is+1
                  CRSS_x(IS,KKK)=TAU_x(IS,0,iphel)
     $                 +HPFAC_x(IS,iphel)/SQRT(GRSZE(iphel))
! let's not use this for the moment.
c$$$c     thermal contribution by applying an exponential factor (Arrhenius factor)
c$$$                  CRSS_x(IS,KKK)=TAU_x(IS,0,iphel)
c$$$     $                 *exp(-(temp+tempdot*dtime-tref)/beta(im))
c$$$c               if (kkk.eq.1.and.is.eq.1) then
c$$$c                  write(*,*)'** crss_x(is,kkk):',crss_x(is,kkk)
c$$$c               endif
                  if (isnan(crss_x(is,kkk)).or.crss_x(is,kkk).le.0) then
                     write(*,*)'something went wrong in'
                     write(*,*)'update_crss_voce_thermo_fe'
                     write(*,*)'-----------'
                     write(*,*)'iphel:',iphel
                     write(*,*)'crss_x(is,kkk):',crss_x(is,kkk)
                     write(*,*)'TAU_x(IS,0,iphel):',TAU_x(IS,0,iphel)
                     write(*,*)'HPFAC_x(IS,iphel):',HPFAC_x(IS,iphel)
                     write(*,*)'SQRT(GRSZE(iphel)):',SQRT(GRSZE(iphel))
                     write(*,*)'HPFAC_x(IS,iphel)/SQRT(GRSZE(iphel))'
                     write(*,*)HPFAC_x(IS,iphel)/SQRT(GRSZE(iphel))
                     call xit()
                  endif
               ENDDO
            enddo

            iprx=iprint*(1/kkk)*(kkk/1)
            if(iprx.ne.0) then
               write(*,'('' initial CRSS for grain # '', i4)') kkk
               write(*,'(9f7.1)') (crss_x(is,kkk),is=1,nsyst(iphel))
            endif
         ENDDO
      ENDIF                     ! END OF IOPTION=1
C ***************************************************************************

      IF (IOPTION.EQ.2) THEN
         if (verbose)then
            call w_chrc(0,'** Update crss voce with ioption=2')
            call w_val(0,'dtime:',dtime)
         endif
         DO KKK=igb,igt
            GAMTOTX=GTOTGR_x(KKK)
            DELTGAM=0d0
            DO IS=1,NSYST(iphel)
               DELTGAM=DELTGAM+ABS(GAMDOT_x(IS,Kkk))*dtime
            ENDDO

            is=0
            do im=1,nmodes(iphel)
               DO jS=1,nsm(im,iphel)!NSYST(iphel)
                  is=is+1
                  DTAU=0.
                  DO kS=1,NSYST(iphel)
                     DTAU=DTAU+HARD_x(IS,kS,iphel)*
     $                    ABS(GAMDOT_x(kS,Kkk))*dtime
                  ENDDO
                  TAU0_ =TAU_x(IS,0,iphel)
                  TAU1 =TAU_x(IS,1,iphel)
                  THET0=THET_x(IS,0,iphel)
                  THET1=THET_x(IS,1,iphel)
                  TINY=1.E-4*TAU0_

                  VOCE=0.0
                  IF(ABS(THET0).GT.TINY) THEN
                     VOCE=THET1*DELTGAM
                     IF(ABS(TAU1).GT.TINY) THEN
                        FACT=ABS(THET0/TAU1)
                        EXPINI=EXP(-GAMTOTX*FACT)
                        EXPDEL=EXP(-DELTGAM*FACT)
                        VOCE  =VOCE-(FACT*TAU1-THET1)/FACT*EXPINI*
     $                       (EXPDEL-1.)-THET1/FACT*EXPINI*
     $                       (EXPDEL*((GAMTOTX+DELTGAM)*FACT+1.)
     $                       -(GAMTOTX*FACT+1.))
                     ENDIF
                  ENDIF
                  if (deltgam.gt.0)then
                     CRSS_x(IS,KKK)=CRSS_x(IS,KKK)+DTAU*VOCE/DELTGAM
                  endif

c     Let's not use this for the moment;
c$$$c     thermal contribution by applying an exponential factor
c$$$                  if (abs(tempdot).gt.0) then
c$$$                     crss_x(is,kkk)=crss_x(is,kkk)
c$$$     $                    +tau_x(is,0,iphel)*
c$$$     $                    exp(-(temp+tempdot*dtime-tref)/beta(im))
c$$$     $                    /(-beta(im))*tempdot*dtime
c$$$                  endif

                  if (isnan(crss_x(is,kkk)).or.crss_x(is,kkk).le.0)
     $                 then
                     write(*,*)'crss is nan dtau*voce/deltgam'
                     write(*,*)'dtime:',dtime
                     write(*,*)'tempdot:',tempdot
                     write(*,*)'dtau:',dtau
                     write(*,*)'voce:',voce
                     write(*,*)'deltgam:',deltgam
                     write(*,*)'dtaue(is,kkk):',dtaue(is,kkk)
                     write(*,*)'crss_x(is,kkk):',crss_x(is,kkk)
                     call xit()
                  endif
               ENDDO! js
            enddo ! im
            GTOTGR_x(KKK)=GAMTOTX+DELTGAM
         ENDDO
      ENDIF   ! END OF IOPTION=2
C ***************************************************************************

      RETURN
      END subroutine update_crss_voce_thermo_fe



c----------------------------------------------------------------------c
      subroutine get_voce_dtdg(gt,th0,th1,tau0,tau1,dtdg)
      implicit none
      real*8, intent(in)::gt,th0,th1,tau0,tau1
      real*8, intent(out)::dtdg
c     locals
      real*8 fact,expterm
c     ------------------------------------------------------------------
c     statements to fool the compiler (not used but can be used sometimes)
      fact=tau0
c     ------------------------------------------------------------------
      if (abs(tau1).gt.1e-10) then
         fact=abs(th0/tau1)
         expterm=exp(-fact*gt)
         dtdg=th1+(fact*tau1-th1+fact*th1*gt)*expterm
      else
         dtdg=0.
c$$$         fact=abs(th0/tau1)
c$$$         expterm=exp(-fact*gt)
c$$$         dtdg=th1+(fact*tau1-th1+fact*th1*gt)*expterm
      endif

      return
      end subroutine
c----------------------------------------------------------------------c
