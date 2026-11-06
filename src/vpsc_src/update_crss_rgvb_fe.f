C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     SUBROUTINE UPDATE_CRSS_rgvb_fe
C
C *** Hardening model based on forward and reversible dislocation
c        evolution on individual slip systems
C *** coded by k.kitayama (July 2011)
C     --> K. Kitayama et al,Intl J of Plasticity 46 (2013) 54-69
C *** modified by CNT for revised constitutive law (July 2012)
C *** Modified by (WW) for Mg and cycling of LC steel (2013,2014)
C     --> W. Wen et al, Intl J of Plasticty xx (2015) to appear
C     --> W. Wen et al, under preparation
C *** Modified by CNT for adapting to latest VPSC7d version (May 2015)
c *** Added a feature (iopt=3) to determine CRSS for the given strain path
c     to probe PCYS (Y. Jeong, 2015 May)
c *** Change in back stress definition - part of the revision in
c     Wei Wen's paper Acta Mat paper.
C
C if IOPTION=0: reads parameters, calculates elastic shear moduli.
C if IOPTION=1: initializes TAU & RHO using input parameters.
C if IOPTION=2: updates disloc densities RHO and CRSS in every system
c if IOPTION=3: determine the crss for the given strain path
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE UPDATE_CRSS_rgvb_fe(
c          in
     $     UR_SX,ioption,iphel,igb,igt,dtime,nsysmx_x,ngrmx,nphmx,
     $     nmodmx_x,ntwmm_x,nsyst_x,
     $     nrs_x,
     $     nslmod_x,ntwmod_x,nslsys_x,ntwsys_x,nmodes_x,nsm_x,ngr_x,
     $     dnca_x,jran_x,wgt_x,
c     inout
     $     taue_x,crss_x,gamdot_x,ishear_tag,ibs_law,xMu,xB,
     $     xf_self,xf_lat,xD,xK,alpham,xlambda,tau0,xp,xbx,facs,pow,
     $     xRho_min,xRho_max,xRho_infl,rho_f,rho_r,taud,rho_prev,
     $     rho_deb_rgvb,tauini,rhos,rho_mod,fac,for_rhos,rev_rhos,
     $     cts,for_rhosl,rev_rhosl,AMSS,deltaD,xp_g,
     $     drho_tot,debs,rho_ps,tautwint,

c          out
     $     isectw_x,twthres_x,itwinlaw_x,
     $     rho_f_ave,rho_r_ave,rho_tot_ave,p_ave)
c     $     rho_tot_pcys_wgt,rho_for_pcys_wgt,rho_r_pcys_wgt,backsurface) !!(tot, for, rev)
      implicit none
      integer, intent(in):: UR_SX,ioption,iphel,nmodmx_x,ntwmm_x,ngr_x,
     $     ngrmx,nphmx
     $     ,nsyst_x(nphmx),nsysmx_x,ntwmod_x(nphmx),
     $     nslsys_x(nphmx),
     $     ntwsys_x(nphmx),
     $     nmodes_x(nphmx),nsm_x(nmodmx_x,nphmx),
     $     jran_x,igb,igt
      integer,intent(out):: nrs_x(nsysmx_x,nphmx),itwinlaw_x,
     $     isectw_x(nsysmx_x,nphmx)
      integer,intent(inout):: nslmod_x(nphmx),
     $     ishear_tag(NSYSMX_X,NGRMX),ibs_law ! 1/29/2016
      real*8,intent(in)::dtime,dnca_x(3,nsysmx_x,nphmx),wgt_x(ngrmx)

      real*8,intent(out)::twthres_x(2,ntwmm_x,nphmx)
      real*8,intent(inout)::taue_x(nsysmx_x,ngrmx),
     $     crss_x(nsysmx_x,1:ngrmx),gamdot_x(nsysmx_x,ngrmx),
     $     xMu(nphmx,NMODMX_X),xB(nphmx,NMODMX_X),tau0(nphmx,NMODMX_X),
     $     xf_self(nphmx,NMODMX_X),xf_lat(nphmx,NMODMX_X),
     $     xD(nphmx,NMODMX_X),
     $     xK(nphmx,NMODMX_X),
     $     alpham(nphmx,NSYSMX_X,NSYSMX_X),
     $     xlambda(nmodmx_x,NGRMX),
     $     xp(nsysmx_x,NGRMX),xbx(nsysmx_x),facs(nphmx,nsysmx_x),
     $     pow(nphmx),
     $     xRho_min(nphmx,NMODMX_X),xRho_max(nphmx,NMODMX_X),
     $     xRho_infl(nphmx,NMODMX_X),
     $     rho_f(NSYSMX_X,NGRMX),rho_r(2*NSYSMX_X,NGRMX),
     $     taud(nsysmx_x),
     $     rho_prev(NSYSMX_X,NGRMX),
     $     rho_deb_rgvb(nsysmx_x,NGRMX),
     $     tauini(nphmx,NSYSMX_X),rhos(nmodmx_x),rho_mod(nmodmx_x),
     $     fac(nphmx,nmodmx_x),for_rhos(nmodmx_x),
     $     rev_rhos(nmodmx_x),
     $     cts(NSYSMX_X,NSYSMX_X),
     $     for_rhosl(nmodmx_x),rev_rhosl(nmodmx_x),
     $     AMSS(NSYSMX_X,NSYSMX_X),deltaD(NSYSMX_X,NGRMX),
     $     xp_g(NGRMX),


     $     drho_tot(NSYSMX_X,NGRMX),
     $     debs(NGRMX),rho_ps(NSYSMX_X,NGRMX),
     $     tautwint(nsysmx_x)   ! 3/20/2015



c     -------------------------
      real*8 p_ave,p_std,rho_r_ave,drho_r_stored_ave,
     $     debs_ave,drho_r_annihi_ave,rho_f_g,rho_r_g,
     $     drho_r_stored,drho_r_annihi,gamdmax,d_gamg,bl,d_gams,
     $     rho_tot_g,rho_tot_ave,drho_r_p,drho_r_m,rho_eff,ran2,
     $     rand,drho_f,rho_tot_s,rho_f_ave,xarg,factorx,annn
      integer is,i,it,kmo,isy,kkk,im,no_pos_gamdot,no_neg_gamdot,
     $     no_neu_gamdot,istag,isp,n_sys_tot,ishear_cur,ism,nog,jm,
     $     isectwx,labelwrite,lun
CNT --> ARRAYS BELOW ARE RELATED TO CALCULATION OF SHEAR MODULUS
c      real*8 C4CA(3,3,3,3),SCHCAX(5),SCH33(3,3),aux55(5,5),
c     $     aux3333(3,3,3,3),aux6(6),aux33(3,3)
c      integer ishear_prv
      real*8 rho_tg(ngr_x),rho_for_pcys_wgt,rho_r_pcys_wgt,
     $     rho_tot_pcys_wgt,backsurface
      character(len=100) :: memo
      real*8,parameter :: gamd_tiny =1.d-10 ! 1e-6 (original);
      real*8,parameter :: rho_tiny  =1.d7
      common /krgvb/annn

c      call w_chrc(0,'** Entering inside update_crss_rgvb')
c     ibs_law (0: 2015 version of back stress, 1: 2016 version of back stress)
c-----------------------------------------------------------------------
c Reads hardening law parameters
c-----------------------------------------------------------------------
      logical verbose
      parameter(verbose=.false.)
c      parameter(verbose=.true.)

      if (verbose) write(*,*) '** Entering update_crss_rgvb'

      if (ioption.eq.0) then

         twthres_x(:,:,:)=0d0
         ISECTWX=0
         NRS_x(1:nsyst_x(iphel),iphel)=20
         read(UR_SX,*) (fac(iphel,JM),JM=1,Nslmod_X(iphel))
         read(UR_SX,*) pow(iphel),factorx,annn
         read(UR_SX,*) nog
         read(UR_SX,*) ibs_law
         DO im=1,nslmod_x(iphel)
            READ(UR_SX,*) memo
            READ(UR_SX,*) xMu(iphel,im) ! mu [MPa]
            READ(UR_SX,*) xB(iphel,im) ! b [m]
            READ(UR_SX,*) tau0(iphel,im) ! t0 at 20 deg
            READ(UR_SX,*) xK(iphel,im) ! K at 20 deg
            READ(UR_SX,*) xD(iphel,im) ! D [m]
            READ(UR_SX,*) xf_self(iphel,im) ! f_self
            READ(UR_SX,*) xf_lat(iphel,im) ! f_lat
            READ(UR_SX,*) xrho_min(iphel,im)
            READ(UR_SX,*) xrho_max(iphel,im)
         ENDDO
         IF(NTWMOD_X(iphel).NE.0) then
            DO IM=1,ntwmod_x(iphel)
               READ(UR_SX,*) memo
               READ(UR_SX,*) xMu(iphel,im) ! mu [MPa]
               READ(UR_SX,*) xB(iphel,im) ! b [m]
               READ(UR_SX,*) tau0(iphel,im) ! t0 at 20 deg.
               READ(UR_SX,*) xK(iphel,im) ! K at 20 deg.
               READ(UR_SX,*) xD(iphel,im) ! D [m]
               READ(UR_SX,*) ISECTWX,ITWINLAW_x,
     $              twthres_x(1,im,iphel),twthres_x(2,im,iphel)
cnt --> these twinning-related statements need to be checked
C          istart=0 ; iend=0
C          do ii=1,IM-1
C            istart=istart+NSM_X(IM,iphel)
C            iend=iend+NSM_X(IM+1,iphel)
C          enddo
C          ISECTW(1+istart:iend,iphel)=ISECTWX
cnt --> for now no secondary twinning is allowed
               ISECTW_x(1:nsyst_x(iphel),iphel)=0
            ENDDO
         ENDIF
C *** READS LATENT HARDENING MATRIX
         read(UR_SX,'(a)') memo     !! reads the header line
         do im=1,nslsys_x(iphel)
            read(UR_SX,*) (alpham(iphel,im,i),i=1,nslsys_x(iphel))
         enddo
C *** READS TWIN-SLIP HARDENING COEFFiCIENTS
         cts(:,:)=0d0
         if(ntwsys_x(iphel).ne.0) then
            read(UR_SX,*) memo
            do it=nslsys_x(iphel)+1,nsyst_x(iphel)
               read(UR_SX,*) (CTS(it,is),is=1,nslsys_x(iphel))
            enddo
         endif
         ISM=0
         DO KMO=1,nmodes_x(iphel) !!
         DO ISY=1,NSM_X(KMO,iphel)
            ISM=ISM+1
            facs(iphel,ism)=fac(iphel,kmo)
         enddo
         enddo
C *** CALCULATES 4TH ORDER ELASTIC STIFFNESS TENSOR FROM 6X6 TENSOR here?
         do kmo=1,nmodes_x(iphel)
            xrho_infl(iphel,kmo)=(xrho_max(iphel,kmo)
     $           +xrho_min(iphel,kmo))/2.d0
         do kkk=igb,igt
            xlambda(kmo,kkk)=xD(iphel,kmo)
         enddo
         enddo
      ENDIF                     ! END OF IOPTION.EQ.0
c------------------------------------------------------------------------
c *** Initialization of CRSS tau for the first step
c------------------------------------------------------------------------
      if(ioption.eq.1) then
         rho_f(:,:)       =rho_tiny
         rho_r(:,:)       =rho_tiny
         rho_prev(:,:)    =rho_tiny
         rho_tg(:)        =rho_tiny
         rho_deb_rgvb(:,:)=rho_tiny
         ishear_tag(:,:)  =0
c         ishear_tag_pcys0(:,:) = 0
c      The loop over phase should be an external one for FE.
         DO IS=1,NSLSYS_X(iphel)
            DO ISP=1,NSLSYS_X(iphel)
c     Factor that determines the angle between glide planes. (1-|n . n|)
               AMSS(IS,ISP)=1-ABS(DNCA_x(1,IS,iphel)*
     $              DNCA_x(1,ISP,iphel)+
     $              DNCA_x(2,IS,iphel)*DNCA_x(2,ISP,iphel)+
     $              DNCA_x(3,IS,iphel)*DNCA_x(3,ISP,iphel)) !!2/12/2014
               IF (ABS(AMSS(IS,ISP)).lt.0.00001d0) AMSS(IS,ISP)=0.d0
            enddo
         enddo
         DO KKK=igb,igt!1,ngr_x
            ISM=0
         DO KMO=1,Nmodes_X(iphel)
            XARG = -3.d0*(xrho_max(iphel,kmo)+xrho_min(iphel,kmo))
     $           /(xrho_max(iphel,kmo)-xrho_min(iphel,kmo))
c            write(*,*)'xarg:',xarg
         DO IS=1,NSM_X(KMO,iphel)
            ISM=ISM+1
            XP(ISM,KKK)=0.5d0-0.5d0 * tanh(XARG) !!2/12/2014
            TAUE_x(ISM,KKK)=tau0(iphel,KMO)
            tauini(iphel,ISM)=tau0(iphel,KMO)
            xBx(ism)=xB(iphel,kmo)
            CRSS_X(ISM,KKK)=tau0(iphel,KMO)
         enddo
         enddo
         enddo
c         call xit()
      ENDIF                     ! END OF IOPTION=1




















c-------------------------------------------------------------------------
c *** Evolves the dislocation densities
c-------------------------------------------------------------------------
      IF (IOPTION.EQ.2) THEN

         if (pow(iphel).eq.0)then
            write(*,*)'pow(iphel) is zero'
            write(*,*)'pow(iphel)',pow(iphel)
            call xit()
         endif

!     Dislocation densities and reversibility parameter p
!     averaged over the grains
         rho_f_ave  =0.d0
         rho_r_ave  =0.d0
         rho_tot_ave=0.d0
         do is=1,nslmod_x(iphel)       ! set rhos for each system
            rho_mod(is) =0.0d0
            for_rhos(is)=0.0d0
            rev_rhos(is)=0.0d0
         enddo
         drho_r_stored_ave=0d0
         drho_r_annihi_ave=0d0
         p_ave=0.d0
         p_std=0.d0
         debs_ave=0.0d0
         no_pos_gamdot=0
         no_neg_gamdot=0
         no_neu_gamdot=0
c      The loop over phase should be an external one for FE.

         DO 200 KKK=igb,igt!1,ngr_x
            do is=1,nslmod_x(iphel)    ! set rhos for each mode
               rhos(is)=0.d0
               for_rhosl(is)=0.0d0
               rev_rhosl(is)=0.0d0
            enddo
            rho_f_g  =0.d0
            rho_r_g  =0.d0
            rho_tot_g=0.d0
            drho_r_stored  =0.d0
            drho_r_annihi  =0.d0
C *** DROPS SYSTEMS WITH NEGLIGIBLE SHEAR RATES
            GAMDMAX=0.D0
            DO IS=1,NSYST_X(iphel)
               IF(ABS(GAMDOT_X(IS,Kkk)).GT.GAMDMAX)
     $              GAMDMAX=GAMDOT_X(IS,Kkk)
            ENDDO
            D_GAMG=0d0
            DO IS=1,NSYST_X(iphel)
               IF(ABS(GAMDOT_X(IS,Kkk)).LT.0.001*GAMDMAX)
     $              GAMDOT_X(IS,KKK)=0d0
               D_GAMG=D_GAMG+ABS(GAMDOT_X(IS,KKK))*dtime
            ENDDO
            ISM=0
         DO KMO=1,NSLMOD_X(iphel)
         DO ISY=1,NSM_X(KMO,iphel)
            ISM=ISM+1
            bL=xB(iphel,KMO)*xlambda(kmo,KKK)
c             Positive shear sense
            IF(GAMDOT_X(ISM,KKK) > gamd_tiny) THEN
               D_GAMS=GAMDOT_X(ism,kkk)*dtime

C *** RESETS THE SHEAR TAG IF THE SENSE IS INVERTED WRT PREVIOUS STEP
c     yj: RESET of shear tag are repeated in all gamdot - this can be replaced out of the end-if
c     and commonly shared by all cases.
               ISTAG=ISHEAR_TAG(ISM,KKK) ! previous shear tag
c     if shear tag has never been clarified or -1
               IF(ISTAG.EQ.0 .OR. ISTAG.EQ.-1) THEN
                  rho_prev(ism,kkk)=rho_f(ism,kkk)
     $                 +rho_r(2*ism-1,kkk)+rho_r(2*ism,kkk)
                  ISHEAR_TAG(ISM,KKK)=1
               ENDIF
               drho_f  = (1.d0-xp(ISM,kkk))*D_GAMS/bL
     $              -xf_self(iphel,KMO)*rho_f(ISM,KKK)*D_GAMG
c     xp: reversibility fraction pertaining to each slip system
               drho_r_p= xp(ISM,kkk)*D_GAMS/bL
     $              -xf_self(iphel,KMO)*rho_r(2*ISM-1,KKK)*D_GAMG
               drho_r_m= -D_GAMS/bL*
     $              (rho_r(2*ISM,KKK)/rho_prev(ism,kkk))**annn
c     for later stats analysis
               drho_r_stored=drho_r_stored+drho_r_p
               drho_r_annihi=drho_r_annihi+drho_r_m
c     negative shear sense
            ELSE IF(GAMDOT_X(ISM,KKK) .lt. -gamd_tiny) THEN
               D_GAMS=abs(GAMDOT_X(ISM,KKK))*dtime
C     *** RESETS THE SHEAR TAG IF THE SENSE IS INVERTED WRT PREVIOUS STEP
               ISTAG=ISHEAR_TAG(ISM,KKK)
               IF(ISTAG.EQ.0 .OR. ISTAG.EQ.1) THEN
                  rho_prev(ism,kkk)=rho_f(ism,kkk)
     $                 +rho_r(2*ism-1,kkk)+rho_r(2*ism ,kkk)
                  ISHEAR_TAG(ISM,KKK)=-1
               ENDIF
               drho_f  = (1.d0-xp(ISM,kkk))*D_GAMS/bL
     $              -xf_self(iphel,KMO)*rho_f(ISM,KKK)*D_GAMG
               drho_r_p= -D_GAMS/bL
     $              *(rho_r(2*ISM-1,KKK)/rho_prev(ism,kkk))**annn
               drho_r_m= xp(ISM,kkk)*D_GAMS/bL
     $              -xf_self(iphel,KMO)*rho_r(2*ISM,KKK)*D_GAMG
               drho_r_stored=drho_r_stored+drho_r_m
               drho_r_annihi=drho_r_annihi+drho_r_p
c     -gamd_tiny <= gamdot_x(ism,kkk) <= +gamd_tiny
            ELSE
c     if gamdot was negligible
C     *** RESETS THE SHEAR TAG IF THE SENSE CHANGED WRT PREVIOUS STEP
               ISTAG=ISHEAR_TAG(ISM,KKK)
               IF(ISTAG.EQ.1.OR.ISTAG.EQ.-1) ISHEAR_TAG(ISM,KKK)=0
c     shear tag saved for PCYS calculation
c                 ishear_tag_pcys0(ism,kkk) = ishear_tag(ism,kkk)
               drho_f=   -xf_lat(iphel,KMO)*rho_f(ism,kkk)*D_GAMG
               drho_r_p= -xf_lat(iphel,kmo)*rho_r(2*ism-1,kkk)*D_GAMG
               drho_r_m= -xf_lat(iphel,kmo)*rho_r(2*ism,kkk)*D_GAMG
            ENDIF
            rho_ps(ism,kkk)=rho_f(ism,kkk)
     $           +rho_r(2*ism-1,kkk)+rho_r(2*ism ,kkk)
            rho_f(ism,kkk)     = rho_f(ism,kkk)     + drho_f
            rho_r(2*ism-1,kkk) = rho_r(2*ism-1,kkk) + drho_r_p
            rho_r(2*ism,kkk)   = rho_r(2*ism,kkk)   + drho_r_m


            if(rho_f(ism,kkk).lt.0.d0)    rho_f(ism,kkk)    =0.d0
            if(rho_r(2*ism-1,kkk).lt.0.d0)rho_r(2*ism-1,kkk)=0.d0
            if(rho_r(2*ism,kkk).lt.0.d0)  rho_r(2*ism,kkk)  =0.d0
            drho_tot(ISM,KKK)=rho_f(ism,kkk)+rho_r(2*ism-1,kkk)
     $           +rho_r(2*ism ,kkk)-rho_ps(ism,kkk) ! 2/12/2014
            if(kmo.eq.1.or.kmo.eq.2.or.kmo.eq.3)then
               rhos(kmo)=rhos(kmo)+rho_f(ism,kkk)+
     $              rho_r(2*ism-1,kkk) + rho_r(2*ism,kkk)
               for_rhosl(kmo)=for_rhosl(kmo)+rho_f(ism,kkk)
               rev_rhosl(kmo)=rev_rhosl(kmo)
     $              +rho_r(2*ism-1,kkk) + rho_r(2*ism,kkk)
            endif
            rho_f_g=rho_f_g+rho_f(ism,kkk)
            rho_r_g=rho_r_g+rho_r(2*ism-1,kkk)+rho_r(2*ism,kkk)
         enddo                  ! ISY
         enddo                     ! KMO






c     ** udpates crss of slip systems using reclaculated disl dens
         rho_tot_g=rho_f_g+rho_r_g
         DO 100 IS=1,NSLSYS_X(iphel) ! CNT
            RHO_EFF=0d0
            rho_tot_s=0d0
            do isp=1,nslsys_x(iphel)
               RHO_EFF=RHO_EFF+alpham(iphel,is,isp)*(rho_f(isp,kkk)+
     $              rho_r(2*isp-1,kkk)+rho_r(2*isp,kkk))
            ENDDO
            if(RHO_EFF.lt.0.d0) RHO_EFF=0.d0 ! avoids negative sqrt
c     Mu was forced to be that of the first slip mode
            TAUE_x(IS,KKK)=tauini(iphel,is)
     $           +xMu(iphel,1)*xBx(is)*sqrt(RHO_EFF)
            rho_tot_s = rho_f(is,kkk)+rho_r(2*is-1,kkk)+rho_r(2*is,kkk)
c     if positive gamdot
            if(gamdot_x(is,kkk).gt.gamd_tiny) then
               no_pos_gamdot=no_pos_gamdot+1
c     *** 2015 version of the back stress
               if (ibs_law.eq.0) then
c                 ww 01/21/2016
                  CRSS_X(is,kkk)=TAUE_x(is,kkk)-facs(iphel,is)
     $                 *TAUE_x(is,kkk)*(RHO_R(2*IS,kkk)/RHO_TOT_G)
     $                 **pow(iphel)
               elseif (ibs_law.eq.1) then
c     *** 2016 version of the back stress
                  CRSS_X(is,kkk)=TAUE_x(is,kkk)-
     $                 facs(iphel,is)*TAUE_x(is,kkk)* ! ww 01/21/2016
     $                 (RHO_R(2*IS,kkk)/rho_tot_s)**pow(iphel)
               endif
            elseif(gamdot_x(is,kkk).lt.-gamd_tiny)then
               no_neg_gamdot=no_neg_gamdot+1
c     *** 2015 version of the back stress
               if (ibs_law.eq.0) then
                  CRSS_X(is,kkk)=TAUE_x(is,kkk)-
     $                 facs(iphel,is)*TAUE_x(is,kkk)* ! ww 01/21/2016
     $                 (RHO_R(2*IS-1,kkk)/RHO_TOT_G)**pow(iphel)
c     *** 2016 version of the back stress
               elseif(ibs_law.eq.1) then
                  CRSS_X(is,kkk)=TAUE_x(is,kkk)-
     $                 facs(iphel,is)*TAUE_X(is,kkk)*
     $                 (RHO_R(2*IS-1,kkk)/rho_tot_s)**pow(iphel)
               endif
            else
               no_neu_gamdot=no_neu_gamdot+1
               CRSS_X(is,kkk)=TAUE_X(is,kkk)
            endif

            taud(is)=CRSS_X(is,1)-TAUE_X(is,1)
            if(abs(crss_x(is,kkk)/TAUE_X(is,kkk)).lt.0.1) then
               write(*,*) 'is=',is, '  kkk=',kkk,
     $              ' crss=',crss_x(is,kkk),' TAUE_X=',Taue_X(is,kkk)
            endif
 100     enddo                  ! is





C     *** UPDATES CRSS of TWINNING SYSTEMS USING Beyerlein et al equation
         tautwint(:)=0d0
         IF(NTWSYS_x(iphel).GT.0) THEN
C     *** THIS ALGORITHM (ESPECIALLY THE rand PART) NEEDS TO BE CHECKED
            DO 105 IT=NSLSYS_X(iphel)+1,NSYST_X(iphel)
               tautwint(it)=0.d0
               TAUE_X(IT,KKK)=tauini(iphel,IT)
            DO IS=1,NSYST_X(iphel)
               tautwint(it)=tautwint(it)+xMu(iphel,1)*xBx(it)*xBx(is)
     $              *cts(it,is)*rho_f(is,kkk)+rho_r(2*is-1,kkk)
     $              +rho_r(2*is,kkk)
            ENDDO
            rand=2*RAN2(JRAN_x)-1.0d0 ! rand =(-1.1) uapt
            CRSS_X(IT,KKK)=TAUE_X(IT,KKK)+rand*TAUE_X(IT,KKK)*factorx
 105        ENDDO
         ENDIF



C     *** UPDATES MEAN FREE PATH IN THE GRAIN USING RECALCULATED DISLOC DENSITIES
c     xK and xD were assummed to be that of the first slip mode
         do kmo=1,nmodes_x(iphel)
            xlambda(kmo,kkk)=1.d0/
     $           (sqrt(rho_tot_g)/xK(iphel,kmo)+1.d0/xD(iphel,kmo))
         enddo



C     *** UPDATES THE DEBRIS DISLOCATION DENSITY and THE REVERSIBILITY
C         FACTOR P
         DO 110 IS=1,NSLSYS_X(iphel)
            deltaD(IS,kkk)=0.d0
         DO ISP=1,NSLSYS_X(iphel)
c     2/12/2014
            deltaD(IS,kkk)=deltaD(IS,kkk)+(rho_f(isp,kkk)
     $           +rho_r(2*isp-1,kkk)+rho_r(2*isp,kkk))*AMSS(IS,ISP)
         ENDDO
         deltaD(IS,kkk)=deltaD(IS,kkk)*drho_tot(IS,KKK)
C     rho_deb is no longer debris. It is D now!
c     2/12/2014
         rho_deb_rgvb(IS,kkk)=rho_deb_rgvb(IS,kkk)+deltaD(IS,kkk)
         if (rho_deb_rgvb(IS,kkk).lt.0d0) then
            rho_deb_rgvb(IS,kkk)=0.01d0 !! 09/12/2014 ww
         endif
c     Warning Assumes xp is regarded to the first slip mode
         xp(is,kkk)=0.5d0-0.5d0*tanh(3.d0*((xrho_max(iphel,1)+
     $        xrho_min(iphel,1))/
     $        (xrho_max(iphel,1)-xrho_min(iphel,1)))
     $        *((sqrt(rho_deb_rgvb(IS,kkk))-
     $        xrho_infl(iphel,1))/xrho_infl(iphel,1)))
 110  enddo                     ! end of is




C     *** PERFORMS SOME STATISTICS FOR OUTPUT PURPOSES
         xp_g(kkk)=0d0
         debs(kkk)=0d0
         DO IS=1, NSLSYS_X(iphel)
            xp_g(kkk)=xp_g(kkk)+xp(IS,kkk)
            debs(kkk)=debs(kkk)+rho_deb_rgvb(IS,kkk)
         ENDDO
         xp_g(kkk)=xp_g(kkk)/NSLSYS_X(iphel)
         debs(kkk)=debs(kkk)/NSLSYS_X(iphel)
         p_ave=p_ave+xp_g(KKK)   *wgt_x(kkk)
         debs_ave=debs_ave+debs(kkk)*wgt_x(kkk) !ww 09/12/2014
C        p_std=p_std+xp_g(KKK)**2*wgt(kkk)
         RHO_F_AVE=RHO_F_AVE+RHO_F_G      *WGT_X(KKK)
         RHO_R_AVE=RHO_R_AVE+RHO_R_G      *WGT_X(KKK)
         RHO_TOT_AVE=RHO_TOT_AVE+RHO_TOT_G*WGT_X(KKK)
         do is=1,NSLMOD_X(iphel)
            rho_mod(is)=rho_mod(is)+rhos(is) *WGT_X(KKK)
c     uapt 3/12/2013
            for_rhos(is)=for_rhos(is)+for_rhosl(is)*wgt_x(kkk)
c     uapt 3/12/2013
            rev_rhos(is)=rev_rhos(is)+rev_rhosl(is)*wgt_x(kkk)
         enddo
         drho_r_stored_ave=drho_r_stored_ave+drho_r_stored*wgt_x(kkk)
         drho_r_annihi_ave=drho_r_annihi_ave+drho_r_annihi*wgt_x(kkk)
 200     enddo !! end of big loop over all grains.

         p_std=sqrt(abs(p_std-p_ave**2)) ! avoid negative sqrt
      ENDIF                     ! END OF IF(IOPTION.EQ.2)






































c     Updates only the back-stress and CRSS by comparing shear_tag_pcys
c     Finds slip systems and collect dislocations that experience
      if (ioption.eq.3) then
c**************************************************
c     Based on the current shear_tag determine
c     back stress and crss -
***************************************************
         backsurface = 0d0      !! YJ 04/25/2016
c**      Work-around the problem with multiphelase case. !! YJ/HK 02/22/2017
C        Need attention later if ihardlaw is phase-specific
c      The loop over phase should be an external one for FE.
         do 300 kkk=igb,igt!1,ngr_x
            ism=0
            rho_tot_g = rho_tg(kkk)
         do kmo=1,nslmod_x(iphel)
         do isy=1,nsm_x(kmo,iphel)
            ism=ism+1
            if (gamdot_x(ism,kkk).gt.gamd_tiny) then
c *** 2015 version of the back stress
               if (ibs_law.eq.0) then
                  CRSS_X(ism,kkk)=TAUE_X(ism,kkk)-
     $                 facs(iphel,ism)
     $                 *TAUE_X(ism,kkk)*(RHO_R(2*ISM,kkk)
     $                 /RHO_TOT_G)**pow(iphel)
c     *** 2016 version of the back stress
               elseif (ibs_law.eq.1) then
                  CRSS_X(ism,kkk)=TAUE_X(ism,kkk)-
     $                 facs(iphel,ism)
     $                 *TAUE_X(ism,kkk)*(RHO_R(2*ISM,kkk)/
     $                 (rho_f(ism,kkk)+rho_r(2*ism-1,kkk)
     $                 +rho_r(2*ism,kkk)))**pow(iphel)
                  backsurface=backsurface+
     $                 (facs(iphel,ism)*TAUE_X(ism,kkk)
     $                 *(RHO_R(2*ISM,kkk)/(rho_f(ism,kkk)
     $                 +rho_r(2*ism-1,kkk)+rho_r(2*ism,kkk)))
     $                 **pow(iphel))
     $                 *wgt_x(kkk)
               endif
            elseif (gamdot_x(ism,kkk).lt.-gamd_tiny) then
c     *** 2015 version of the back stress
               if (ibs_law.eq.0) then
                  CRSS_X(ism,kkk)=TAUE_X(ism,kkk)-
     $                 facs(iphel,ism)
     $                 *TAUE_X(ism,kkk)*(RHO_R(2*ISM-1,kkk)
     $                 /RHO_TOT_G)**pow(iphel)
               elseif (ibs_law.eq.1) then
c     *** 2016 version of the back stress
                  CRSS_X(ism,kkk)=TAUE_X(ism,kkk)-
     $                 facs(iphel,ism)
     $                 *TAUE_X(ism,kkk)*(RHO_R(2*ISM-1,kkk)/
     $                 (rho_f(ism,kkk)+rho_r(2*ism-1,kkk)
     $                 +rho_r(2*ism,kkk)))**pow(iphel)
                  backsurface=backsurface+
     $                 (facs(iphel,ism)*TAUE_X(ism,kkk)
     $                 *(RHO_R(2*ISM-1,kkk)/(rho_f(ism,kkk)
     $                 +rho_r(2*ism-1,kkk)+rho_r(2*ism,kkk)))
     $                 **pow(iphel))
     $                 *wgt_x(kkk)
               endif
            else
               CRSS_X(ism,kkk)=TAUE_X(ism,kkk) ! no back stress
            endif
         enddo
         enddo
 300     enddo

c         n_shear_inm=0
         n_sys_tot = 0
c       ## stats in dislcoations densities related with Bauschinger effect.
c       write rho_tot (total), rho_tot_r (reversible)
c       that relates with Bauschinger effect.
         rho_for_pcys_wgt = 0d0
         rho_tot_pcys_wgt = 0d0
         rho_r_pcys_wgt   = 0d0
c     do 200 iphel=iphelbot,ipheltop ! 1, nph
c**      Work-around the problem with multiphelase case. !! YJ/HK 02/22/2017
C        Need attention later if ihardlaw is phase-specific
         do 310 kkk=igb,igt!1,ngr_x
            ism=0
            rho_tg(kkk) = 0.d0
         do kmo=1,nslmod_x(iphel)
         do isy=1,nsm_x(kmo,iphel)
            ism=ism+1
            n_sys_tot = n_sys_tot + 1
            rho_tg(kkk) = rho_tg(kkk) + rho_f(ism,kkk) +
     $           rho_r(2*ism-1,kkk) + rho_r(2*ism,kkk)
            rho_for_pcys_wgt = rho_for_pcys_wgt +
     $           rho_f(ism,kkk)*wgt_x(kkk)
            rho_tot_pcys_wgt = rho_tot_pcys_wgt +
     $           rho_tg(kkk) * wgt_x(kkk)
c     Check current shear_tag
            if (gamdot_x(ism,kkk).gt.gamd_tiny) then
               ishear_cur = + 1
            elseif (gamdot_x(ism,kkk).lt.-gamd_tiny) then
               ishear_cur = - 1
            else
               ishear_cur =  0
            endif
c     Current shear sense and the reverse dislocations that are
c     affecting the back stress
            if (gamdot_x(ism,kkk).gt.0) then
               rho_r_pcys_wgt=rho_r_pcys_wgt+rho_r(2*ism,kkk)*wgt_x(kkk)
            elseif (gamdot_x(ism,kkk).lt.0) then
               rho_r_pcys_wgt=rho_r_pcys_wgt+rho_r(2*ism-1,kkk)
     $              *wgt_x(kkk)
            endif
         enddo
         enddo
 310     enddo
      endif                     ! end of ioption.eq.3

      if (ioption.eq.4) then
c     Remove back-stress from CRSS.
c**      Work-around the problem with multiphelase case. !! YJ/HK 02/22/2017
C        Need attention later if ihardlaw is phase-specific
c$$$         do 110 iphel=iphelbot,ipheltop
c$$$            iphel=iphel-iphelbot+1
c      The loop over phase should be an external one for FE.
         do kkk=igb,igt!1,ngr_x
            ism=0
         do kmo=1,nslmod_x(iphel)
         do isy=1,nsm_x(kmo,iphel)
            ism=ism+1
            CRSS_X(ism,kkk)=TAUE_X(ism,kkk)
         enddo
         enddo
         enddo
c     ipcys_dd_rev=1
      endif


      RETURN
      END SUBROUTINE UPDATE_CRSS_rgvb_fe






















C *** CALCULATES 4TH ORDER ELASTIC STIFFNESS TENSOR FROM 6X6 TENSOR
c$$$         CALL VOIGT(AUX6,AUX33,CELCCV,C4CA,3) ! CNT Sep 2011
c$$$         ISX=0
c$$$c         write(*,*)'HERE1  in crss_dd_rev'
c$$$         DO IM=1,nmodes_x(iphel)
c$$$c            write(UW_LOG,*)
c$$$c            write(UW_LOG,'('' MODE'',I5,''   ELASTIC SHEAR MODULUS READ''
c$$$c     $,F13.5)') IM,xMu(iphel,IM)
c$$$
c$$$            ISX=ISX+NSM_X(IM,iphel)
c$$$C *** CALCULATES 2nd ORDER SCHMID TENSOR FROM 5-DIM VECTOR
c$$$            DO I=1,5
c$$$               SCHCAX(I)=SCHCA(I,ISX,IPHEL)
c$$$            ENDDO
c$$$            CALL CHG_BASIS(SCHCAX,SCH33,AUX55,AUX3333,1,5)
c$$$C     *** CALCULATES SHEAR MODULUS FOR MODE 'IM' BY PROJECTING ELASTIC STIFFNESS
c$$$            xMu_x(:)=0.
c$$$            do i=1,3
c$$$            do j=1,3
c$$$            do k=1,3
c$$$            do l=1,3
c$$$               xMu_x=xMu_x+SCH33(i,j)*C4CA(i,j,k,l)*SCH33(k,l)
c$$$            enddo
c$$$            enddo
c$$$            enddo
c$$$            enddo
c$$$c     $           ' CALCULATED ELASTIC SHEAR MODULUS''
c$$$c     $,F13.5)') IM,xMu_x
c$$$CNT --> NEED TO CLARIFY INCONSISTENCY BETWEEN READ AND CALC SHEAR MODULUS
c$$$
c$$$C         xMu(iphel,IM)=xMu_x   ! NOT USED FOR THE TIME BEING
c$$$         enddo
