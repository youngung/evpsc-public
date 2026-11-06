C *****************************************************************************
C     SUBROUTINE UPDATE_CRSS_DD_fe      -->     VERSION 5/FEB/2022

C     Updates the CRSS taking into account dislocation density evolution.
C
C     Version of May 2011 by Capolungo,Beyerlein,Kaschner,Tome (MSEA 2009) who
C     replaced 'mode densities' by 'system densities' in the original equations
C     of Beyerlein and Tome (IJP 2008).

C     Hall-Petch-effect algorithms commented out in this version.
C     Dimensioning of DD model parameters valid only for single phase PX.

C     Average shear nodulus S_MOD replaced (CT: 2021 12 12). Now uses the
C     anisotopic shear moduli of each mode in the shear plane & along shear
C     direction MU_MODE(KMOD,IPH) calculated from the anisotropic SX elastic
C     ctes inside SUBR DATA_CRYSTAL. Another option could be to use the average
C     shear modulus C44ISO corresponding to a random PX.


c     interfaced with dEVPSC code by Y. Jeong (FEB. 2022)
c      GRSZE is phase-specific
c      ** Temperature is assumed constant dtau/dT should be obtained to work with
c         changes in temperature on the fly.

c     dislocation static recovery was modified by Y. Jeong (Oct-Nov. 2023)
c     to enforce evenly distributed recovery rate among all slip modes to mitigate
c     the preferred recovery for those modes with a high burgers vector
C *****************************************************************************

      SUBROUTINE UPDATE_CRSS_DD_fe_new(
c     intent(in)
     $     IOPTION,ife,UR_SX,UW_DIS,UW_TAU,nphmx,nsysmx,nslmod,nslsys,
     $     nmodes,
     $     nsyst,nsm,ntwmod,temp,ngr,wgt,ngrmx,DBAR,time_incr,nmodmx,
     $     gamdot,NTWMMX,mu_mode,
c     intent(inout)
c          read from (ioption.eq.0) and not evolving
     $     grsze,chifor,chisub,RHO_initW,TEMPref,burg,actener,K1gener,
     $     DRAG,EDOT0,rho_inits,Atau0,Atau1,Atau2,Fdeb0,Fdeb1,Fdeb2,
     $     HPcoef,T0nucl,T1nucl,T2nucl,T0prop,T1prop,T2prop,ctwsl,CTS1,
     $     CTS0,CTS2,K3recov,K4recov,D0vac,Efvac,Emvac,
c         subject to 'evolution'
     $     rho_deb,rho_s_for,rho_s_deb,
c         constants
     $     edotref,
c     intent(out)
     $     itwinlaw,hard,twthres,isectw,nrs,crss,
     $     RHOSAT,TAU_PROP,TAU_NUCL,
c     intent(out) for stats of DD law
     $     rho_for_ph_mod,rho_deb_ph,rho_tot_ph,rho_deb_ph_mod,
     $     rho_recov,taue_wavg,rhonucl,pnucl
     $     )

c     rho_s_for       : foreset dislocation for each slip system
c     rho_s_deb       : debris  dislocation for each slip system
c     rho_for         : forest  dislocation for each grain
c     rho_deb         : debris  dislocation for each grain

c     rho_deb_ph      : Phase-specific wegithed average of rho_deb
c     rho_deb_ph_mod  : Phase-specific wegithed average of rho_s_deb for each slip mode
c     rho_for_ph_mod  : Phase-specific wegithed average of rho_deb
      implicit none

c     intent(in)
      integer, intent(in)::ioption,ife,nphmx,nsysmx,nmodmx,NTWMMX,
     $     UR_SX,UW_DIS,UW_TAU,ngrmx,nslsys(nphmx),nslmod(nphmx),
     $     nmodes(nphmx),nsyst(nphmx),nsm(nmodmx,nphmx),ntwmod(nphmx),
     $     ngr(0:nphmx)
      real*8, intent(in)::temp,wgt(ngrmx),DBAR(6),time_incr,
     $     gamdot(nsysmx,ngrmx),mu_mode(nmodmx,nphmx)
c     intent(inout)
      real*8, intent(inout):: grsze(nphmx),chifor,chisub,RHO_initW,
     $     TEMPref,burg(nmodmx),actener(nmodmx),K1gener(nmodmx),
     $     DRAG(nmodmx),EDOT0(nmodmx),rho_inits(nmodmx),Atau0(nmodmx),
     $     Atau1(nmodmx),Atau2(nmodmx),Fdeb0(nmodmx),Fdeb1(nmodmx),
     $     Fdeb2(nmodmx),HPcoef(nmodmx),T0nucl(nmodmx),T1nucl(nmodmx),
     $     T2nucl(nmodmx),T0prop(nmodmx),T1prop(nmodmx),T2prop(nmodmx),
     $     CTWSL(NTWMMX,NMODMX),CTS1(NTWMMX,NMODMX),CTS2(NTWMMX,NMODMX),
     $     CTS0(NTWMMX,NMODMX),K3recov,K4recov,D0vac,Efvac,Emvac,
     $     rho_deb(ngrmx),rho_s_for(NSYSMX,NGRMX),
     $     rho_s_deb(nsysmx,ngrmx)
c          constants that are not subjected to evolution.
     $     ,edotref
c     intent(out)
      integer, intent(out)::nrs(nsysmx,nphmx),isectw(nsysmx,nphmx),
     $     itwinlaw
      real*8, intent(out)::hard(nsysmx,nsysmx,nphmx),
     $     twthres(2,ntwmmx,nphmx),crss(nsysmx,ngrmx),
     $     RHOSAT(nmodmx),TAU_PROP(nmodmx),TAU_NUCL(nmodmx)
c-----------------------------------------------------------------------
c     local
      integer i,iphel,iph,nslsysx,nslmodx,nmodesx,nrsx,is,js,kk,kmo,ktm,
     $     ibeg,isectwx,ksmo,imox,ktmx,kmox,iend,kkk,ism,isy,imo,jsm,
     $     jmo
      real*8 boltz_mpa,boltz_ev,twthres2x,twthres1x,edot,vnorm,
     $     delgam(nsysmx),DELRHO_for(NSYSMX),
     $     DELRHO_deb(NSYSMX),sh_mod,rho_tot_ph(nphmx),rho_deb_ph(nphmx)
     $     ,rho_deb_ph_mod(nmodmx,nphmx),rho_for_ph_mod(nmodmx,nphmx),
     $     atau(nmodmx),
     $     TAUE(nsysmx,ngrmx),RHOnucl,FDEB(nmodmx),K2recom(nmodmx),
c          recovery
     $     factor3,factor4,Dvac,rho_deb_recov,rho_for_recov,
     $     rho_s_deb_rec,rho_s_for_rec,factor,fraction,burg_avg,mu_avg
!     6 modes x 4 contributions (Peierls, for, sub, HP) for slip
c             x 3 contributions (nucl, prop, slip) for twin
     $     ,taue_wavg(nphmx,6,4),pnucl,rho_for(ngrmx),rho_totx,
     $     rho_tot_mode(nmodmx),rho_recov_rate_mode(nmodmx),
     $     rho_recov_rate,rho_recov
      character*254 prosa
c-----------------------------------------------------------------------
c     Static recovery case study
      integer irecov_case
c     0: default
c     1: hardwiring roof slip's burger vector for all
c     2: using an averaged burger's vector
c      parameter(irecov_case=0) ! use mode-specific bugers vector in the rate covery term
c      parameter(irecov_case=1) ! use the burgers vector of 4-th slip mode
      parameter(irecov_case=2) ! use the average burgers

      taue_wavg(:,:,:)=0.
      RHOnucl=0.
      Pnucl=0.

      BOLTZ_MPA=1.380622E-29    ! units of MPa*m^3/K
      BOLTZ_eV =8.617E-05       ! units of eV/K
C ***************************************************************
C *** READS PARAMETERS OF CONSTITUTIVE MODEL
C ***************************************************************

      IF (IOPTION.EQ.0) THEN
C *** THE 'CHILD' (twin) PHASE HAS INITIAL ZERO WEIGHT, AND CRSS's
C     HAVE TO BE REDEFINED WHEN THE TWINS ARE CREATED (NOT IMPLEMENTED)
         IPHEL=1
         IPH  =1                ! algorithm limited to single phase case
         NSLSYSX=NSLSYS(IPHEL)
         NSLMODX=NSLMOD(IPHEL)
         NMODESX=NMODES(IPHEL)

         NRSX=20                ! power to use in 'rate sensitive' law
         DO IS=1,NSYST(IPH)
            NRS(IS,IPH)=NRSX
            DO JS=1,NSYST(IPH)
               HARD(IS,JS,IPH)=1.0
            ENDDO
         ENDDO

         if (ife.eq.0) then
            OPEN(unit=UW_DIS,file='DIS_DEN.OUT',status='unknown',
     $           action='write')
            OPEN(unit=UW_TAU,file='TAUE.OUT',status='unknown',
     $           action='write')
            WRITE(UW_DIS,'(14a12)') ('Rho_for_mod',kk=1,NSLMODX)
     $        ,'Rho_deb','Rho_tot',
     $        ('rho_deb_mod',kk=1,nslmodx),
     $        'drho_rec'
         !! slips
         ! peierls
         do imo=1,nslmodx
            write(UW_TAU,'(a11,i1.1)',advance='no') 'Peierls',imo
         enddo
         ! Forest
         do imo=1,nslmodx
            write(UW_TAU,'(a11,i1.1)',advance='no') 'Forest',imo
         enddo
         ! debris (substructure)
         do imo=1,nslmodx
            write(UW_TAU,'(a11,i1.1)',advance='no') 'Deb',imo
         enddo
         ! Hall-Petch
         do imo=1,nslmodx
            write(UW_TAU,'(a11,i1.1)',advance='no') 'HP',imo
         enddo


         !! twins
         ! nucleation
         do imo=nslmodx+1,nmodesx
            write(UW_TAU,'(a11,i1.1)',advance='no') 'nucl',imo
         enddo
         ! propagation
         do imo=nslmodx+1,nmodesx
            write(UW_TAU,'(a11,i1.1)',advance='no') 'propg',imo
         enddo
         ! slip-twin interaction
         do imo=nslmodx+1,nmodesx
            write(UW_TAU,'(a11,i1.1)',advance='no') 'slipintx',imo
         enddo


         !! etc.
         write(UW_TAU,'(2a12)',advance='no') 'RHOnuc','Pnucl'
         do imo=1,nslmodx
            write(UW_TAU,'(a11,i1.1)',advance='no') 'rhosat',imo
         enddo
         write(UW_TAU,*)
         endif


C *** READS HARDENING PARAMETERS FOR SLIP AND TWINS (Beyerlein & Tome, 2008)
c     elastic shear modulus redundant (as it is overriden by MU_MODE)
         READ(UR_SX,*) SH_MOD,GRSZE(iphel)
         READ(UR_SX,*) CHIfor,CHIsub
         READ(UR_SX,*) RHO_initW,TEMPref,EDOTref
c     vacancy diff at 0K, formtn & migrtn energy
         READ(UR_SX,*) D0vac,Efvac,Emvac
c     static recovery law parameters (6,0.6 were used for Fe by Kohnert and Capolungo 2022)
         READ(UR_SX,*) K3recov,K4recov

         DO KMO=1,NSLMODX
            READ(UR_SX,'(A)') PROSA
            READ(UR_SX,*) BURG(KMO),ACTENER(KMO)
            READ(UR_SX,*) K1gener(KMO),DRAG(KMO)
            READ(UR_SX,*) EDOT0(KMO)
            READ(UR_SX,*) RHO_INITS(KMO)
            READ(UR_SX,*) Atau0(KMO),Atau1(KMO),Atau2(KMO)
            READ(UR_SX,*) Fdeb0(KMO),Fdeb1(KMO),Fdeb2(KMO)
            READ(UR_SX,*) HPcoef(KMO)
         ENDDO

c        This assumes that the twin modes always appear later than
c        slip modes in single crystal file
         IF(NTWMOD(IPHEL).NE.0) THEN
            IBEG=NSLSYSX+1
            DO KMO=NSLMODX+1,NMODESX
c     from 1 update to n of twin modes
               KMOX=KMO-NSLMODX
               READ(UR_SX,'(A)') PROSA
               READ(UR_SX,*) BURG(KMO)
c     nucleation threshold
               READ(UR_SX,*) T0nucl(KMO),T1nucl(KMO),T2nucl(KMO)
c     propagation threshold
               READ(UR_SX,*) T0prop(KMO),T1prop(KMO),T2prop(KMO)
c     tw-sl coupling
               READ(UR_SX,*) (CTS0(KMOX,KSMO),KSMO=1,NSLMODX)
c     tw-sl coupling
               READ(UR_SX,*) (CTS1(KMOX,KSMO),KSMO=1,NSLMODX)
c     tw-sl coupling
               READ(UR_SX,*) (CTS2(KMOX,KSMO),KSMO=1,NSLMODX)
               READ(UR_SX,*) ISECTWX,ITWINLAW,TWTHRES1X,TWTHRES2X
               TWTHRES(1,KMOX,IPHEL)=TWTHRES1X
               TWTHRES(2,KMOX,IPHEL)=TWTHRES2X
               IEND=IBEG-1+NSM(KMO,IPHEL)
               ISECTW(IBEG:IEND,IPHEL)=ISECTWX
               IBEG=IEND+1
            ENDDO
         ENDIF

      ENDIF   ! END OF ioption=0

C **************************************************************************
C *** CALCULATES INITIAL CRSS (TEMPERATURE DEPENDENT) AND INITIALIZES
C *** DISLOCATION DENSITIES IN ALL SYSTEMS OF EACH PHASE IN ALL GRAINS
C **************************************************************************

      IF(IOPTION.EQ.1 .OR. IOPTION.EQ.2) THEN

         if (temp.le.1e-10) then
            write(*,*)'ioption:',ioption
            write(*,*)'temp:',temp
            write(*,*)' if temperature is zero, error would occur '
            write(*,*)' temp should be given in Kelvin.'
            call xit()
         endif

         IPH=1
         IPHEL=1
         NSLMODX=NSLMOD(IPHEL)
         NMODESX=NMODES(IPHEL)

         EDOT=SQRT(2./3.)*VNORM(DBAR,5) !! caution when dbar = 0.
c         EDOT=SQRT(2./3.)*VNORM(DBAR,6) !! caution when DBAR = 0.

         if (ioption.eq.1) edot=1e-3
         Dvac=D0vac*exp(-(Efvac+Emvac)/(BOLTZ_eV*TEMP))

         RHO_TOT_PH(IPHEL)=0d0
         RHO_DEB_PH(IPHEL)=0d0
         DO KMO=1,NSLMODX
            RHO_FOR_PH_mod(KMO,IPHEL)=0d0
            RHO_deb_PH_mod(KMO,IPHEL)=0d0
C --> EQ 3.17 B&T --> TAU_0(T)
            ATAU(KMO)=Atau0(KMO)+Atau1(KMO)*EXP(-TEMP/Atau2(KMO))
         ENDDO

C --> EQ 3.26 B&T --> TAU_0_beta(T) FOR THE TWIN SYSTEMS
         IF(NTWMOD(IPHEL).NE.0) THEN
            DO KTM=NSLMODX+1,NMODESX
               KTMX=KTM-NSLMODX !! iterate from 1 via # of TWIN modes
               TAU_PROP(KTM)=T0prop(KTM)+
     $              T1prop(KTM)*EXP(-TEMP/T2prop(KTM))
               TAU_NUCL(KTM)=T0nucl(KTM)+
     $              T1nucl(KTM)*EXP(-TEMP/T2nucl(KTM))
c              CTWSL=CTS1*exp((EDOT-EDOTref)/CTS2) but CTWSL=CTS1 if EDOT=EDOTref
c               CTWSL(KTMX,1:NSLMODX)=CTS1(KTMX,1:NSLMODX) ! tw-sl coupling
               ctwsl(ktmx,1:nslmodx)=cts0(ktmx,1:nslmodx)
     $              +cts1(ktmx,1:nslmodx)
     $              *exp((edot-edotref)/cts2(ktmx,1:nslmodx))
            ENDDO
         ENDIF


c$$$         write(*,*)'ctwsl'
c$$$         DO KTM=NSLMODX+1,NMODESX
c$$$            KTMX=KTM-NSLMODX    !! iterate from 1 via # of TWIN modes
c$$$            write(*,*)(ctwsl(ktmx,kmo),kmo=1,nslmodx)
c$$$         enddo
c$$$         call xit()
      ENDIF
C **************************************************************************
C   ASSIGNS INITIAL RHO & CRSS VALUES FOR SLIP AND TWIN SYSTEMS
C **************************************************************************

      IF(IOPTION.EQ.1) THEN
         taue_wavg(:,:,:)=0d0
         DO KKK=NGR(IPH-1)+1,NGR(IPH)
c     determine rho_for

            ism=0
            rho_for(kkk)=0d0
            rho_deb(kkk)=0.
            do kmo=1,nslmodx
            do is=1,nsm(kmo,iph)
               ism=ism+1
               rho_s_for(ISM,KKK)=RHO_INITS(KMO)
               rho_for(kkk)=rho_for(kkk)+rho_s_for(ism,kkk)
               RHO_DEB(KKK)=rho_deb(kkk)+RHO_initW
            enddo
            enddo

            ISM=0
            DO KMO=1,NSLMODX !! number of slip modes.
               DO IS=1,NSM(KMO,IPHEL)
                  ISM=ISM+1
c
                  rho_s_deb(ISM,KKK)=rho_initW
c                  rho_s_deb(ISM,KKK)=rho_initW/nslsys(iphel) ! 2023 07 15
C --> Eq 3.19 & ~3.20 B&T --> TAU_0 + TAU_for + TAU_deb
c                 1. Peierls friction stress
                  TAUE(ISM,KKK)=ATAU(KMO)
                  taue_wavg(iph,kmo,1)=taue_wavg(iph,kmo,1)
     $                 +atau(kmo)*wgt(kkk)/nsm(kmo,iph)
c                 2. Forest dislocations (this should be fixed)
                  TAUE(ISM,KKK)=TAUE(ISM,KKK)+
     $                 CHIfor*BURG(KMO)*MU_MODE(KMO,IPH)
     $                 *sqrt(rho_for(kkk))
                  taue_wavg(iph,kmo,2)=taue_wavg(iph,kmo,2)
     $                 +CHIfor*BURG(KMO)*MU_MODE(KMO,IPH)
     $                 *sqrt(rho_for(kkk))
     $                 *wgt(kkk)/nsm(kmo,iph)
c                 3. debris (substructure)
                  IF(RHO_DEB(KKK).GT. 0.0) then
                     TAUE(ISM,KKK)=TAUE(ISM,KKK)
     $                    -CHIsub*BURG(KMO)*MU_MODE(KMO,IPH)
     $                    *SQRT(RHO_DEB(KKK))*
     $                    LOG(BURG(KMO)*SQRT(RHO_DEB(KKK)))
                     taue_wavg(iph,kmo,3)=taue_wavg(iph,kmo,3)
     $                    -CHIsub*BURG(KMO)*MU_MODE(KMO,IPH)
     $                    *SQRT(RHO_DEB(KKK))
     $                    *LOG(BURG(KMO)*SQRT(RHO_DEB(KKK)))
     $                    *wgt(kkk)/nsm(kmo,iph)
                  endif
C --> Eq 3.21 B&T; Eq 9 Capolungo et al --> TAU_0 + TAU_for + TAU_deb + TAU_hp
c                 4. HP effect
                  taue(ism,kkk)=taue(ism,kkk)+
     $                 HPcoef(KMO)*MU_MODE(KMO,IPH)
     $                 *SQRT(BURG(KMO)/GRSZE(iph))
                  taue_wavg(iph,kmo,4)=taue_wavg(iph,kmo,4)
     $                 +HPcoef(KMO)*MU_MODE(KMO,IPH)
     $                 *SQRT(BURG(KMO)/GRSZE(iph))
     $                 *wgt(kkk)/nsm(kmo,iph)

c                 ** CRSS = Peierls + Forest + Debris + HP effect
                  CRSS(ISM,KKK)=TAUE(ISM,KKK)
                  RHO_FOR_PH_mod(KMO,IPHEL)=RHO_FOR_PH_mod(KMO,IPHEL)+
     $                 WGT(KKK)*rho_s_for(ISM,KKK)/nsm(kmo,iphel)
                  rho_deb_ph_mod(kmo,iphel)=rho_deb_ph_mod(kmo,iphel)+
     $                 wgt(kkk)*rho_s_deb(ism,kkk)/nsm(kmo,iphel)
               ENDDO            ! end of IS

c-----------------------------------------------------------------------
c     YJ: the calculation of rhosat should be done out side of kkk loop.
c              rhosat = (k1/k2)^2
c              (blow expression is actually k2/k1 not k2)
               if (abs(edot).lt.1e-10) then
                  k2recom(kmo)=0.
               else
                  k2recom(KMO)=CHIfor*BURG(KMO)/ACTENER(KMO)
c                                                         DRAG [MPa]
     $                 *(1.-BOLTZ_MPA*TEMP/(BURG(KMO)**3)/DRAG(KMO)
c                      BOLTZ_MPA/BURG**3 [MPa/K]
     $                 *LOG(EDOT/EDOT0(KMO)))
               endif
c              the above is actually k2/k1
               rhosat(kmo)=k2recom(kmo)**(-2.0)
               k2recom(kmo)=k2recom(kmo)*k1gener(kmo)
            ENDDO               ! end of KMO
c-----------------------------------------------------------------------

            RHO_DEB_PH(IPHEL)=RHO_DEB_PH(IPHEL)+wgt(KKK)*RHO_DEB(KKK)
            IF(NTWMOD(IPHEL).NE.0) THEN
               DO IMO=NSLMODX+1,NMODESX   !IMO : loop from the end of slip mode  over twin modes
                  IMOX=IMO-NSLMODX        !IMOX: loop from 1                     over twin modes
                  DO IS =1,NSM(IMO,IPHEL) !IS  : loop from 1 over twin systems
                     ISM=ISM+1            !ISM : loop from the end of slip syst. over twin system
                     TAUE(ISM,KKK) =0d0
                     RHOnucl=0d0
                     JSM=0
                     DO JMO=1,NSLMODX
                        DO JS =1,NSM(JMO,IPHEL) !
                           JSM=JSM+1
c --> Eq 14 Capolungo et al
                           TAUE(ISM,KKK)=TAUE(ISM,KKK)+
     $                          MU_MODE(JMO,IPH)*
     $                          BURG(IMO)*BURG(JMO)*
     $                          CTWSL(IMOX,JMO)*rho_s_for(JSM,KKK)
                           RHOnucl=RHOnucl+
     $                          rho_s_for(JSM,KKK)/RHOSAT(JMO)
                        ENDDO
                     ENDDO
c     slip
                     taue_wavg(iph,imo,3)=taue_wavg(iph,imo,3)+
     $                    taue(ism,kkk)*wgt(kkk)/nsm(imo,iph)
C --> Eq 3-26 B&T; Eq 15 Capolungo et al
                     TAUE(ISM,KKK)=TAUE(ISM,KKK)+TAU_PROP(IMO)+
     $                    (TAU_NUCL(IMO)-TAU_PROP(IMO))*exp(-RHOnucl)
c                    nucleation
                     taue_wavg(iph,imo,1)=taue_wavg(iph,imo,1)
     $                    +tau_nucl(imo)*(1d0-exp(-RHOnucl))*wgt(kkk)
c                    propagation
                     taue_wavg(iph,imo,2)=taue_wavg(iph,imo,2)
     $                    +tau_prop(imo)*(1d0-exp(-RHOnucl))*wgt(kkk)
C --> EQ 3.27 B&T --> TAU_0_beta + TAU_hp
                     CRSS(ISM,KKK)=TAUE(ISM,KKK)
                     if (isnan(crss(ism,kkk))) then
                        write(*,*)'crss is nan'
                        write(*,*)'ioption=1 in update_crss_dd_fe_new'
                        call xit()
                     endif
                  ENDDO
               ENDDO
            ENDIF
         ENDDO                  ! end of KKK

c     needs to fix below!!
         DO KMO=1,NSLMODX
            RHO_TOT_PH(IPHEL)=RHO_TOT_PH(IPHEL)
     $           +RHO_FOR_PH_mod(KMO,IPHEL)
         ENDDO
         RHO_TOT_PH(IPHEL)=RHO_TOT_PH(IPHEL)+RHO_DEB_PH(IPHEL)
      ENDIF      ! END OF IOPTION=1







C **************************************************************************
C *** UPDATES DISLOCATION DENSITIES AND CRSS's OF SLIP & TWIN SYSTEMS
C **************************************************************************
      IF (IOPTION.EQ.2) THEN
c         write(*,*)'ioption.eq.2 begins in devpsc'
c         write(*,'(a10,e11.3)')'rho_s_for:',rho_s_for(1,1)
c         call xit()

         taue_wavg(:,:,:)=0.0
         DO KMO=1,NSLMODX
c     I think this statement is redundant (yj)
            RHO_FOR_PH_mod(KMO,IPHEL)=0.0
c     I think this statement is redundant (yj)
            rho_deb_ph_mod(kmo,iphel)=0.0
c --> Eq 11 Knezevic et al (2013)
            FDEB(KMO)=Fdeb0(KMO)
     $           +Fdeb1(KMO)*exp((TEMPref-TEMP)/Fdeb2(KMO))
c --> Eq 3-12 B&T; Eq 10 Capolungo et al (blow expression is actually k2/k1)
c           rhosat = (k1/k2)^2
c           (blow expression is actually k2/k1 not k2)
            if (abs(edot).lt.1e-10) then
               k2recom(kmo)=0.
            else
               k2recom(KMO)=CHIfor*BURG(KMO)/ACTENER(KMO)
c               DRAG [MPa]
     $              *(1.-BOLTZ_MPA*TEMP/(BURG(KMO)**3)/DRAG(KMO)
c     BOLTZ_MPA/BURG**3 [MPa/K]
     $              *LOG(EDOT/EDOT0(KMO)))
            endif
c              the above is actually k2/k1
            rhosat(kmo)=k2recom(kmo)**(-2.0)
            k2recom(kmo)=k2recom(kmo)*k1gener(kmo)
         ENDDO

c     I think this statement is redundant (yj)
         RHO_DEB_PH(IPHEL)=0.
c     I think this statement is redundant (yj)
         RHO_TOT_PH(IPHEL)=0.
         if (.false.) then
            RHO_FOR_RECOV=0.
            RHO_DEB_RECOV=0.
         else
            rho_recov=0.
         endif


c *** Updates dislocation densities (forest, debris, and static recovery)
         DO 100 KKK=NGR(IPH-1)+1,NGR(IPH)
            RHO_DEB(KKK)=0.
c           total disl. density of the current grain (kkk)
            rho_totx=0
            rho_tot_mode(:)=0.
            ism=0
            do kmo=1,nslmodx
               do isy=1,nsm(kmo,iphel)
                  ism=ism+1
                  rho_totx=rho_totx+
     $                 rho_s_for(ism,kkk)+rho_s_deb(ism,kkk)
                  rho_tot_mode(kmo)=rho_tot_mode(kmo)+
     $                 rho_s_for(ism,kkk)+rho_s_deb(ism,kkk)
               enddo
            enddo
            ISM=0
c           getting average burgers vector and mu
c            if (.false.) then
            if (.true.) then
               burg_avg=0.
               mu_avg  =0.
               DO KMO=1,NSLMODX
                  burg_avg=burg_avg+burg(kmo)
                  mu_avg =mu_avg +mu_mode(kmo,iph)
               enddo
               burg_avg=burg_avg/nslmodx
               mu_avg =mu_avg /nslmodx
            else
c           or using the Roof's burgers vector and mu
               burg_avg=burg(4)
               mu_avg  =mu_mode(4,iph)
            endif



            rho_recov_rate=0.
            rho_recov_rate_mode(:)=0.
            DO KMO=1,NSLMODX
c  static recovery rate same in each system & based on total density (20230623)
c              use mode-specific burgers vector and shear modulus
               if (irecov_case.eq.0) then
c              * factors for recovery terms
                  FACTOR3=K3recov*Dvac/BURG(KMO)
                  FACTOR4=K4recov*MU_MODE(KMO,IPH)*BURG(KMO)**4
                  if (isnan(factor3))then
                     write(*,*)'factor3:',factor3
                     write(*,*)'dvac:',dvac
                     write(*,*)'burg(kmo):',burg(kmo)
                     call xit()
                  endif
               elseif(irecov_case.eq.1) then
c                 YJ's experiments in which 4-th slip mode (i.e., Roof in <uranium_dd_test_recovery.sx>)
c                 which has the highest burgers vector, thus giving the most significant static recovery.
                  FACTOR3=K3recov*Dvac/burg(4)
                  FACTOR4=K4recov*MU_avg*burg(4)**4
               elseif(irecov_case.eq.2) then
c                 YJ'experiments to see how anisotropy in burgers vector influence by using roof slip information
c                  FACTOR3=K3recov*Dvac/burg_avg , initial trial
c                  FACTOR4=K4recov*MU_avg*burg_avg**4

c     corrected one proposed by Carlos Nov 29, 2023
                  FACTOR3=K3recov*Dvac/burg(kmo)
                  FACTOR4=K4recov*mu_mode(kmo,iph)*burg_avg**4
               endif


c               20241203 expression for recovery by total dislocation density within each slip mode.
               rho_recov_rate_mode(kmo)=factor3*rho_tot_mode(kmo)**(1.5)
     $              *(exp((factor4*sqrt(rho_tot_mode(kmo)))/
     $              (BOLTZ_MPa*TEMP))-1.)
               rho_recov_rate=rho_recov_rate+factor3*rho_totx**(1.5)*
     $              (exp((factor4*sqrt(rho_totx))/(BOLTZ_MPa*TEMP))-1.)
            enddo

            DO KMO=1,NSLMODX

               do isy=1,nsm(kmo,iphel)
                  ism=ism+1
                  delgam(ism)=abs(gamdot(ism,kkk))*time_incr
c --> Eq 9 Capolungo et al
c      this may need to be eventually changed to have 'total' dislocation density for k1 term.
c      (according to the discussion between Laurent and Carlos)
                  DELRHO_for(ISM)=K1gener(KMO)*SQRT(rho_s_for(ISM,KKK))
c     calculate with previous
     $                           -K2recom(KMO)*    (rho_s_for(ISM,KKK))


c --> Eq 3-14 B&T; Eq 11 Capolungo et al
                  DELRHO_deb(ISM)=Fdeb(KMO)*K2recom(KMO)
c     calculate with previous (yj changed)
     $                 *rho_s_for(ISM,KKK)*BURG(KMO)*SQRT(RHO_DEB(KKK))
c     mode-fraction (updates in 20241203 by yj)
                  fraction=rho_s_for(ism,kkk)/rho_tot_mode(kmo)
c     mode-fraction (updates in 20241203 by yj)
                  rho_s_for(ism,kkk)=rho_s_for(ism,kkk)
     $                 +delrho_for(ism)*delgam(ism)
     $                 -fraction*rho_recov_rate_mode(kmo)*time_incr
c     mode-fraction (updates in 20241203 by yj)
                  fraction=rho_s_deb(ism,kkk)/rho_tot_mode(kmo)
                  rho_s_deb(ism,kkk)=rho_s_deb(ism,kkk)
     $                 +delrho_deb(ism)*delgam(ism)
     $                 -fraction*rho_recov_rate_mode(kmo)*time_incr

                  if (rho_s_for(ism,kkk).lt.rho_inits(kmo)*1e-4)
     $                 rho_s_for(ism,kkk)=rho_inits(kmo)*1e-4
                  if (rho_s_deb(ism,kkk).lt.rho_initw*1e-4) then
                     rho_s_deb(ism,kkk)=rho_initw*1e-4
                  endif

c     update in grain
c                  RHO_DEB(KKK)=RHO_DEB(KKK)+DELRHO_deb(ISM)*DELGAM(ISM)
                  RHO_DEB(KKK)=RHO_DEB(KKK)+rho_s_deb(ism,kkk)

c     update in mode & phase
                  RHO_FOR_PH_mod(KMO,IPHEL)=RHO_FOR_PH_mod(KMO,IPHEL)+
     $                 WGT(KKK)*rho_s_for(ISM,KKK)/nsm(kmo,iphel)
c     update in mode & phase
                  rho_deb_ph_mod(kmo,iphel)=rho_deb_ph_mod(kmo,iphel)+
     $                 wgt(kkk)*rho_s_deb(ism,kkk)/nsm(kmo,iphel)

                  FACTOR=TIME_INCR*WGT(KKK)

                  if (.false.) then
                     RHO_FOR_RECOV=RHO_FOR_RECOV+RHO_S_FOR_REC*FACTOR
                     RHO_DEB_RECOV=RHO_DEB_RECOV+RHO_S_DEB_REC*FACTOR
                  else
                     rho_recov=rho_recov+rho_recov_rate*factor
                  endif

               ENDDO            ! end of ISY=1,NSM
            ENDDO               ! end of KMO=1,NSLMODX

c     update in phase
            RHO_DEB_PH(IPHEL)=RHO_DEB_PH(IPHEL)+wgt(KKK)*RHO_DEB(KKK)





C *** UPDATES crss FOR THE SLIP SYSTEMS USING RECALCULATED DISLOC DENSITIES
c     determine rho_for in each grain after updating the dislocation density
            rho_for(kkk)=0d0
            ism=0
            do kmo=1,nslmodx
            do is=1,nsm(kmo,iph)
               ism=ism+1
               rho_for(kkk)=rho_for(kkk)+rho_s_for(ism,kkk)
            enddo
            enddo

            ISM=0
            DO KMO=1,NSLMODX
               DO ISY=1,NSM(KMO,IPHEL)
                  ISM=ISM+1
c --> Eqs 3-19 & 3-20 B&T
c                 1. Peierls friction stress
                  taue(ism,kkk)=atau(kmo)
                  taue_wavg(iph,kmo,1)=taue_wavg(iph,kmo,1)
     $                 +atau(kmo)*wgt(kkk)/nsm(kmo,iph)

c                 2. Forest dislocations (this should be fixed)
                  TAUE(ISM,KKK)=taue(ism,kkk)
     $                 +chifor*burg(kmo)*MU_MODE(KMO,IPH)
     $                 *sqrt(rho_for(kkk))
                  taue_wavg(iph,kmo,2)=taue_wavg(iph,kmo,2)
     $                 +CHIfor*BURG(KMO)*MU_MODE(KMO,IPH)
     $                 *sqrt(rho_for(kkk))
     $                 *wgt(kkk)/nsm(kmo,iph)
c                 3. debris (substructure)

                  IF(RHO_DEB(KKK).GT. 0.0) then
                     TAUE(ISM,KKK)=taue(ism,kkk)-
     $                    chisub*burg(kmo)*MU_MODE(KMO,IPH)
     $                    *sqrt(rho_deb(kkk))*
     $                    log(burg(kmo)*sqrt(rho_deb(kkk)))
                     taue_wavg(iph,kmo,3)=taue_wavg(iph,kmo,3)
     $                    -CHIsub*BURG(KMO)*MU_MODE(KMO,IPH)
     $                    *SQRT(RHO_DEB(KKK))
     $                    *LOG(BURG(KMO)*SQRT(RHO_DEB(KKK)))
     $                    *wgt(kkk)/nsm(kmo,iph)
                  endif
C --> Eq 3.21 B&T; Eq 9 Capolungo et al --> TAU_0 + TAU_for + TAU_deb + TAU_hp
c                 4. HP effect
                  taue(ism,kkk)=taue(ism,kkk)
     $                 +HPcoef(KMO)*MU_MODE(KMO,IPH)
     $                 *SQRT(BURG(KMO)/GRSZE(iph))
                  taue_wavg(iph,kmo,4)=taue_wavg(iph,kmo,4)
     $                 +HPcoef(KMO)*MU_MODE(KMO,IPH)
     $                 *SQRT(BURG(KMO)/GRSZE(iph))
     $                 *wgt(kkk)/nsm(kmo,iph)
c                 ** CRSS = Peierls + Forest + Debris + HP effect
                  CRSS(ISM,KKK)=TAUE(ISM,KKK)
               ENDDO            ! end of ISY=1,NSM
            ENDDO               ! end of KMO=1,NSLMODX




C *** UPDATES crss FOR THE TWINNING SYSTEMS (COUPLED WITH DISL DENS)
            IF(NTWMOD(IPHEL).NE.0) THEN
               RHOnucl=0d0
               JSM=0
               do JMO=1,nslmodx ! over all slip modes
                  do JS=1,NSM(JMO,IPHEL) ! over all slip systems
                     JSM=JSM+1
                     RHOnucl=RHOnucl+rho_s_for(JSM,KKK)/RHOSAT(JMO)
                  enddo
               enddo
               Pnucl=1d0-exp(-RHOnucl)

               DO IMO=NSLMODX+1,NMODESX ! twin mode
                  IMOX=IMO-NSLMODX
                  DO IS =1,NSM(IMO,IPHEL) ! twin system
                     ISM=ISM+1
                     TAUE(ISM,KKK)=0d0
                     JSM=0
                     DO JMO=1,NSLMODX ! loop over disl. slip modes
                        DO JS =1,NSM(JMO,IPHEL)
                           JSM=JSM+1
c --> Eq 14 Capolungo et al
                           TAUE(ISM,KKK)=TAUE(ISM,KKK)+
     $                          MU_MODE(KMO,IPH)*BURG(IMO)*BURG(JMO)*
     $                          CTWSL(IMOX,JMO)*rho_s_for(JSM,KKK)
                        ENDDO
                     ENDDO
c                    slip
                     taue_wavg(iph,imo,3)=taue_wavg(iph,imo,3)
     $                    +taue(ism,kkk)*wgt(kkk)/nsm(imo,iph)
C --> Eq 3-26 B&T; Eq 15 Capolungo et al
                     TAUE(ISM,KKK)=TAUE(ISM,KKK)+TAU_PROP(IMO)+
     $                    (TAU_NUCL(IMO)-TAU_PROP(IMO))*exp(-RHOnucl)
C --> EQ 3.27 B&T --> TAU_0_beta + TAU_hp
                     CRSS(ISM,KKK)=TAUE(ISM,KKK)
                  ENDDO
c                 nucleation
                  taue_wavg(iph,imo,1)=taue_wavg(iph,imo,1)
     $                 +tau_nucl(imo)*Pnucl*wgt(kkk)
c                 propagation
                  taue_wavg(iph,imo,2)=taue_wavg(iph,imo,2)
     $                 +tau_prop(imo)*Pnucl*wgt(kkk)
               ENDDO
            ENDIF               ! end of NTWMOD.ne.0
 100     ENDDO                  ! END OF LOOP kkk OVER GRAINS





c        I need to fix below
         DO KMO=1,NSLMODX
            RHO_TOT_PH(IPHEL)=RHO_TOT_PH(IPHEL)
     $           +RHO_FOR_PH_mod(KMO,IPHEL)
         ENDDO
         RHO_TOT_PH(IPHEL)=RHO_TOT_PH(IPHEL)+RHO_DEB_PH(IPHEL)
         if (ife.eq.0) then
         WRITE(UW_DIS,'(14e12.3e3)') (RHO_FOR_PH_mod(KMO,IPHEL),
     $        KMO=1,NSLMODX)
     $        ,RHO_DEB_PH(IPHEL),RHO_TOT_PH(IPHEL)
     $        ,(rho_deb_ph_mod(kmo,iphel),kmo=1,nslmodx)
     $        ,rho_recov

         write(UW_TAU,'(28e12.3e3)')
     $        ((taue_wavg(iph,imo,i),imo=1,nslmodx)        ,i=1,4),
     $        ((taue_wavg(iph,imo,i),imo=nslmodx+1,nmodesx),i=1,3),
     $        RHOnucl,Pnucl,
     $        (rhosat(imo),imo=1,nslmodx)
         endif

      ENDIF                     ! END OF IF(IOPTION.EQ.2)
C************************************************************************

      RETURN
      END subroutine update_crss_dd_fe_new
