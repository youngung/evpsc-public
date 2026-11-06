      SUBROUTINE EVPSC_INPUT_fe(
c     in
     $     iodf,ife,nphmx,
     $     UR_INP,UR_SX,UR_TX,UR_SHP,UR_ECTE,UW_LOG,UW_DIS,UW_TAU,
     $     dtime,
c     (YJ): Is dtime really required for EVPSC_INPUT_FE?
c     I do not think dtime passed into <EVPSC_INPUT_FE>  will affect its behavior.
     $     nsysmx,ngrmx,ntwmmx,nmodmx,ntwsmx,nneimx,mxndeg,
c     out
     $     isx_ver,ihardlaw,iratesens,errs,errd,errm,errso,
     $     itmaxext,itmaxint,itmaxso,IRSVAR,JXRSINI,JXRSFIN,
     $     JXRSTEP,ibcinv,IRECOVER,ISAVE,ICUBCOM,NWRITE,
     $     interaction,IUPDORI,IUPDSHP,IUPDHAR,NNEIGH,iflu,
     $     icompgrain,
     $
     $     gamd0s,ngr,NPH,ishape,fijph,ifrag,iflat,
     $     crit_shp,axisph,axisgr,eulerph,fijgr,idiff,
     $     selcc,nmodes,nslmod,ntwmod,nsyst,nslsys,ntwsys,
     $     twsh,nsm,wgt,schca,nschca,c6,isense,itwtype,hard,
     $     isectw,crss,twthres,nrs,itwinlaw,ag,icrysymph,
     $     nunique,nrsmin,athcc,celcc,dnca,dbca,neigh,
     $     wneigh,wph,tau_voce,thet_voce,hpfac_voce,taue,
     $     ishear_tag,ibs_law,xmu,xb,xf_self,xf_lat,
     $
     $     xD,xK,alpham,xlambda,xp,xbx,facs,pow,xRho_min,xRho_max,
     $     xRho_infl,rho_f,rho_r,taud,rho_prev,rho_deb_rgvb,tauini,rhos,
     $     rho_mod,fac,for_rhos,rev_rhos,cts,for_rhosl,rev_rhosl,AMSS,
     $     deltaD,xp_g,drho_tot,debs,rho_ps,tautwint,
     $     fstar,filediff
c$$$c          update_crss_dd(old)
c$$$     $     ,burg,actener,aKGENER,DRAG,edoteso,rho_inits,rho_initw
c$$$     $     ,yone,bone,sone,TLATENT,aQ0,aQ1,aQ2,aKOVB3,YCRIT,ACRIT
c$$$     $     ,TCRIT,YPROP,APROP,TPROP,TAU_PROP,TAU_CRIT,atau
c$$$     $     ,rho_deb,rho_s_for       !,TAUE
c$$$     $     ,alpha_int,aQ6

cc     update_crss_dd(new)
     $     ,xmu_mode,
     $     grsze,chifor,chisub,rho_initw,TEMPref,burg,
     $     actener,K1gener,DRAG,EDOT0,
     $     rho_inits,Atau0,Atau1,Atau2,
     $     Fdeb0,Fdeb1,Fdeb2,HPcoef,
     $     T0nucl,T1nucl,T2nucl,T0prop,
     $     T1prop,T2prop,CTWSL,cts1,cts0,cts2,
     $     K3recov,K4recov,D0vac,Efvac,Emvac,
     $     rho_deb,rho_s_for,rho_s_deb,
c          constants from update_crss_dd(new)
     $     edotref
c     polynomial fit of temperature-dependent el modulus and ctes.
     $     ,cijv_fit,cte_fit,ndegs_fit
c     coordinate for selecting ebsd
     $     ,coords,ebsd_rad

c     recrystallization
c          intent(out)
     $     ,irx
     $     ,rexgrowpar,rexnuclpar,radius
     $     )
c     ------------------------------------------------------------


      implicit none
      integer,intent(in)::iodf
      integer, intent(in)::ife,
     $     UR_INP,UR_SX,UR_TX,UR_SHP,
     $     UR_ECTE,UW_LOG,UW_DIS,UW_TAU,
     $     nsysmx,
     $     ngrmx,
     $     nphmx,ntwmmx,nmodmx,ntwsmx,nneimx,mxndeg

      integer, intent(out)::isx_ver,ihardlaw,iratesens,itmaxext,itmaxint
     $     ,itmaxso,IRSVAR,JXRSINI,JXRSFIN,JXRSTEP,ibcinv,IRECOVER,ISAVE
     $     ,ICUBCOM,NWRITE,interaction,IUPDORI,IUPDSHP,IUPDHAR,NNEIGH,
     $     iflu,icompgrain,ishape(0:nphmx),ifrag(nphmx),iflat(0:nphmx),
     $     idiff(nphmx),ngr(0:nphmx),nmodes(nphmx),nslmod(nphmx),
     $     ntwmod(nphmx),nsyst(nphmx),nslsys(nphmx),ntwsys(nphmx),
     $     nsm(nmodmx,nphmx),isense(nsysmx,nphmx),itwtype(nsysmx,nphmx),
     $     isectw(nsysmx,nphmx),itwinlaw,icrysymph(nphmx),nunique,nrsmin
     $     ,nrs(nsysmx,nphmx),neigh(0:nneimx,ngrmx),
     $     ishear_tag(NSYSMX,ngrmx),ibs_law,ndegs_fit(nphmx,2),irx

      real*8, intent(out)::errs,errd,errm,errso,gamd0s(nsysmx,ngrmx),
     $     fijph(3,3,0:nphmx),CRIT_SHP(NPHMX),axisph(0:3,3,0:nphmx),
     $     axisgr(0:3,3,ngrmx),eulerph(3,0:nphmx),fijgr(3,3,ngrmx),
     $     selcc(6,6,nphmx),twsh(ntwmmx,nphmx),wgt(ngrmx),
     $     schca(5,nsysmx,nphmx),nschca(6,nsysmx,nphmx),c6,
     $     hard(nsysmx,nsysmx,nphmx),crss(nsysmx,ngrmx),
     $     twthres(2,ntwmmx,nphmx),ag(3,3,ngrmx),
     $     athcc(6,nphmx),celcc(6,6,nphmx),dnca(3,nsysmx,nphmx),
     $     dbca(3,nsysmx,nphmx),wneigh(0:nneimx,ngrmx),wph(nphmx),
     $     TAUE(NSYSMX,NGRMX),tau_voce(1:nsysmx,0:1,1:nphmx),
     $     thet_voce(1:nsysmx,0:1,1:nphmx),hpfac_voce(1:nsysmx,1:nphmx),
     $     xmu(nphmx,nmodmx),xb(nphmx,nmodmx),xf_self(nphmx,1:nmodmx),
     $     xf_lat(nphmx,1:nmodmx),xD(nphmx,NMODMX),xK(nphmx,NMODMX),
     $     alpham(nphmx,NSYSMX,NSYSMX),xlambda(nmodmx,NGRMX),
     $     xp(nsysmx,NGRMX),xbx(nsysmx),facs(nphmx,nsysmx),
     $     pow(nphmx),xRho_min(nphmx,NMODMX),xRho_max(nphmx,NMODMX),
     $     xRho_infl(nphmx,NMODMX),rho_f(NSYSMX,ngrmx),
     $     rho_r(1:2*NSYSMX,ngrmx),taud(1:nsysmx),
     $     rho_prev(1:NSYSMX,ngrmx),rho_deb_rgvb(1:nsysmx,1:NGRMX),
     $     tauini(1:nphmx,1:NSYSMX),rhos(1:nmodmx),rho_mod(1:nmodmx),
     $     fac(1:nphmx,nmodmx),for_rhos(1:nmodmx),rev_rhos(1:nmodmx),
     $     cts(NSYSMX,NSYSMX),for_rhosl(nmodmx),rev_rhosl(nmodmx),
     $     AMSS(NSYSMX,NSYSMX),deltaD(NSYSMX,NGRMX),xp_g(NGRMX),
     $     drho_tot(NSYSMX,NGRMX),debs(NGRMX),rho_ps(NSYSMX,NGRMX),
     $     tautwint(nsysmx),fstar,
c          polynomial fit of temperature-dependent el modulus and ctes.
     $     cijv_fit(nphmx,0:mxndeg,6,6),cte_fit(nphmx,0:mxndeg,6)
c$$$c     dislocation dd hardening model
c$$$     $     ,burg(nmodmx),actener(nmodmx),aKGENER(nmodmx),DRAG(nmodmx)
c$$$     $     ,edoteso(nmodmx),rho_inits(nmodmx),rho_initw
c$$$     $     ,yone(nmodmx),bone(nmodmx),sone(nmodmx)
c$$$     $     ,TLATENT(NTWMMX,NMODMX),aQ0(nmodmx),aQ1(nmodmx),aQ2(nmodmx)
c$$$     $     ,aKOVB3(nmodmx),YCRIT(NMODMX),ACRIT(NMODMX)
c$$$     $     ,TCRIT(NMODMX),YPROP(NMODMX),APROP(NMODMX),TPROP(NMODMX)
c$$$     $     ,TAU_PROP(nmodmx),TAU_CRIT(nmodmx),atau(nmodmx)
c$$$     $     ,rho_deb(ngrmx),rho_s_for(NSYSMX,NGRMX)!,TAUE(nsysmx,ngrmx)
c$$$     $     ,alpha_int,aQ6,xmu_mode(nmodmx,nphmx)

c     dislocation dd hardening model (new version - 2021 Jan)
     $     ,xmu_mode(nmodmx,nphmx),grsze(nphmx),chifor,chisub,rho_initw,
     $     TEMPref,burg(nmodmx),actener(nmodmx),K1gener(nmodmx),
     $     DRAG(nmodmx),EDOT0(nmodmx),rho_inits(nmodmx),Atau0(nmodmx),
     $     Atau1(nmodmx),Atau2(nmodmx),Fdeb0(nmodmx),Fdeb1(nmodmx),
     $     Fdeb2(nmodmx),HPcoef(nmodmx),T0nucl(nmodmx),T1nucl(nmodmx),
     $     T2nucl(nmodmx),T0prop(nmodmx),T1prop(nmodmx),T2prop(nmodmx),
     $     CTWSL(NTWMMX,NMODMX),CTS1(NTWMMX,NMODMX),CTS0(NTWMMX,NMODMX),
     $     CTS2(NTWMMX,NMODMX),K3recov,K4recov,D0vac,Efvac,Emvac,
     $     rho_deb(ngrmx),
     $     rho_s_for(nsysmx,ngrmx),rho_s_deb(nsysmx,ngrmx),
c            constants
     $     edotref,
c          recrystalization
     $     rexgrowpar(5),rexnuclpar(5),radius(ngrmx)
c----------------------------------------------------------------------c
      character*120,intent(out)::filediff(nphmx)
      real*8, intent(in)::dtime

! #include <SMAAspUserSubroutines.hdr>


c     ----
      real*8 gamdot(nsysmx,ngrmx),gtotgr(ngrmx),fijgr_temp(3,3,1),
     $     tau0(nphmx,nmodmx)
c      real*8 grsze(nphmx)
c     voce parameters
      character*125 filetext,filecrys,fileaxes
c     locals
      character(len=256) prosa
      character(len=256) fn_temp_emod_cte
      real*8 aux6(6),aux66(6,6),dum1,dum2,dum3,dum4,wphtot,xid3(3,3),
     $     aux3(3),aux33(3,3),aux3333(3,3,3,3),tempdum
      integer ncompa,is,i,ioption,IPH,igb,igt,jran,nelem,ierr,
     $     NPH,no_grain_per_phase,nrem,kkk,iostat,no_grain_ori
c     ----
      logical verbose
c      parameter(verbose=.true.)
      parameter(verbose=.false.)
c----------------------------------------------------------------------c
      real*8 rhosat(nmodmx),TAU_PROP(nmodmx),TAU_NUCL(nmodmx)
c     quantities for stats analysis from update_crss_dd_fe_new
      real*8 rho_for_ph_mod(nmodmx,nphmx),rho_deb_ph(nphmx),
     $     rho_tot_ph(nphmx),rho_deb_ph_mod(nmodmx,nphmx),
     $     rho_recov,taue_wavg(nphmx,6,4),RHOnucl,pnucl
c----------------------------------------------------------------------c
c     for iodf=3 case
      integer ncrds
      parameter(ncrds=3)
      real*8, intent(in):: coords(ncrds),ebsd_rad
      integer nlines_ebsd,nrowmx
      parameter(nrowmx=200000)
      real*4 ebsd_raw_data(5,nrowmx)
      common/kebsd/ebsd_raw_data,nlines_ebsd

c----------------------------------------------------------------------c

C *********   INITIALIZATION BLOCK   ***************************
      XID3(:,:)=0d0
      do i=1,3
         xid3(i,i)=1d0
      enddo

C     SEED FOR RANDOM NUMBER GENERATOR (RAN2) (USED FOR TWINNING AND RX)
c     seed for random number - required for iodf option
c      JRAN=-1
c      call random_seed(JRAN)
      call chg_basis(aux6,aux33,aux66,aux3333,0,6)

c     write 'vpsc.in' file to 'run_log.out' file
      if (ife.eq.0) then
         do i=1,100
            read(unit=UR_INP,end=10,fmt='(a)') prosa
            write(UW_LOG,'(a)') prosa
         enddo
 10      rewind UR_INP
      endif
c$$$
C     READS # OF ELEMENTS, # OF PHASES AND RELATIVE VOLUME OF PHASE IN ELEM.
C     PHASES ARE LABELED SEQUENTIALLY FOR A MULTIELEMENT RUN (i.e. IF ELEMENT
C     #1 CONSISTS OF phase1=fcc AND phase2=bcc, THEN phase3 AND phase4 ARE
C     THE fcc AND bcc PHASES IN ELEMENT #2...AND SO ON)

      read(UR_INP,*) nelem
c      call w_ival(0,'nelem:',nelem)
      READ(UR_INP,*) NPH

      if (NPH.gt.nphmx) then
         call w_ival(0,'NPH:',NPH)
         call w_ival(0,'nphmx    :',nphmx)
         call stop_debug(0,'** Error: nphmx should be increased')
      endif

      if (NPH.gt.nphmx) then
         call w_ival(0,'NPH:',NPH)
         call w_ival(0,'nphmx    :',nphmx)
         call stop_debug(0,'** Error: Set NPHMX greater than or equal to
     $ NPH')
      endif

c      call w_ival(0,'NPH:',NPH)
      read(UR_INP,*) (wph(IPH),IPH=1,NPH) ! each fe

      wphtot=0d0
      do IPH=1,NPH
         wphtot=wphtot+wph(IPH)
      enddo
      if(abs(wphtot-1.) .gt. 1.e-5) then
         do IPH=1,NPH
            write(*,*) wph(IPH)
         enddo
         write(*,*) ' --> vol fraction of phases should add to 1 !!'
         write(*,*) ' --> will normalize phase fractions to 1 !!'
         call xit()
         do IPH=1,NPH
            wph(IPH)=wph(IPH)/wphtot
         enddo
      endif

C ***************************************************************************
C     LOOP OVER PHASES.
C     FOR EACH PHASE READS CRYSTAL SYSTEMS, TEXTURE AND GRAIN SHAPES.
C ***************************************************************************

C *** INITIALIZE DEFORMATION GRADIENT (AND ELLIPSOID) OF THE ELEMENT
      ngr(0)=0
      ishape(0)=0
      iflat(:)=0
      fijph(:,:,0)=xid3(:,:)
      do i=1,3
         do iph=0,nph
            axisph(0,i,iph)=1.
         enddo
      enddo

      if (verbose) write(*,*)'Call update_shape(0)'
c      CALL UPDATE_SHAPE (0)
c      call stop_debug(0,'stopped before entering update_shape_fe')
      call update_shape_fe(
c          in
     $     ngrmx,nphmx,0,0,0,interaction,ifrag,crit_shp,ishape,
!          inout
     $     fijph,iflat,axisph,axisgr,
c          out
     $     eulerph,fijgr,ierr)

      if (verbose) then
         write(*,*)'exit update_shape(0)'
         call w_chr(0,'axisph(0,:,0)')
         aux3(:)=axisph(0,:,0)
         call w_dim(0,aux3,3,1d0,.true.)
         call w_chr(0,'axisph(0:3,:,0)')
         aux33(:,:)=axisph(1:3,:,0)
         call w_mdim(0,aux33,3,1d0)
         call w_chr(0,'eulerph(:,0)')
         aux3(:)=eulerph(:,0)
         call w_dim(0,aux3,3,1d0,.true.)
      endif

      do IPH=1,NPH
         read(UR_INP,'(a)') prosa

         ishape(IPH)=0
         ifrag(IPH) =0
         crit_shp(IPH)=25
CFEB
         read(UR_INP,*)ishape(IPH),ifrag(IPH),crit_shp(IPH)
CFEE
         read(UR_INP,*) (axisph(0,i,IPH),i=1,3)
         read(UR_INP,*) (eulerph(i,IPH) ,i=1,3)
         read(UR_INP,'(a)') prosa
         read(UR_INP,'(a)') filetext
         read(UR_INP,'(a)') prosa
         read(UR_INP,'(a)') filecrys
         read(UR_INP,'(a)') prosa
         read(UR_INP,'(a)') fileaxes
c----- Diffraction input files
         read(UR_INP,'(a)') prosa
         read(UR_INP,*) idiff(IPH)
c         read(UR_INP,'(a)') filediff
         read(UR_INP,'(a)') filediff(IPH)
C *** READS SLIP AND TWINNING MODES FOR THE PHASE
         call mutexlock(3)
         open (action='read',unit=UR_SX,file=filecrys,status='old')

c     ----------------------------------------------------------------
c     Test version of single crystal file (isx_ver)
c     see if filecrys includes 'thermal_crystal_property.dat' file
         call sx_version_tester_nsch(UR_SX,isx_ver)

         if (isx_ver.eq.1) then
            open(unit=UR_ECTE,file=fn_temp_emod_cte,
     $           status='old',action='read')
         endif
         call data_crystal_fe( ! phase specific
c             in
     $        UR_SX,UR_ECTE,isx_ver,mxndeg,IPH,ntwmmx,nmodmx,
     $        nsysmx,ntwsmx,nphmx,
c             out
     $        celcc,athcc,selcc,nmodes,nslmod,ntwmod,nsyst,nslsys,

     $        ntwsys,twsh,nsm,dnca,dbca,schca,nschca,c6,isense,itwtype,
     $        xmu_mode,icrysymph,fstar,ndegs_fit,cijv_fit,cte_fit)

         if (isx_ver.eq.1) close(unit=UR_ECTE)

C *** READS HARDENING PARAMETERS OF THE HARDENING LAW FOR THIS PHASE
C       IHARDLAW=0 --> EXTENDED VOCE & PREDOMINANT TWIN REORIENTATION
C       IHARDLAW=1 --> MTS MODEL
C       IHARDLAW=2 --> DISLOCATION DENSITY HARDENING
c       IHARDLAW=3 --> Glenn's thermo voce hardening law
c$$$
         read(UR_SX,*)
         read(UR_SX,*) ihardlaw
         read(UR_SX,*) iratesens
         if(ihardlaw.eq.1 .and. iratesens.eq.1) then
            write(*,*)'FOR MTS MODEL RESETS IRATESENS=0 IN CONST LAW'
            iratesens=0
            call xit()
         endif
         if(ihardlaw/10.eq.2 .and. iratesens.eq.1) then
            write(*,*) ' FOR DD MODEL RESETS IRATESENS=0 IN CONST LAW'
            iratesens=0
            call xit()
         endif

c         IF(IHARDLAW.EQ.0)  CALL UPDATE_CRSS_VOCE (UR_SX,0,IPH,dtime)
         if(ihardlaw.eq.0)  then
            igb=ngr(IPH-1)+1
            igt=ngr(IPH)
            call update_crss_voce_fe(
c                in
     $           UR_SX,0,IPH,nphmx,ngrmx,dtime,nmodmx,ntwmmx,igb,igt,
     $           nmodes,nsysmx,nsm,nslsys,nsyst,gamdot,
c                inout
     $           grsze,gtotgr,tau_voce,thet_voce,
     $           hpfac_voce,hard,crss,
c                out
     $           twthres,nrs,itwinlaw,isectw)
         elseif(ihardlaw.eq.3) then
            write(*,*)'** before entering update_crss_voce_thermo_fe'
            igb=ngr(IPH-1)+1
            igt=ngr(IPH)
            call update_crss_voce_thermo_fe(
c                in
     $           UR_SX,ife,IPH,nphmx,ngrmx,dtime,nmodmx,ntwmmx,
     $           igb,igt,
     $           nmodes,nsysmx,nsm,nslsys,nsyst,gamdot,0d0,0d0,
c                inout
     $           grsze,gtotgr,tau_voce,thet_voce,
     $           hpfac_voce,hard,crss,
c                out
     $           twthres,nrs,itwinlaw,isectw)
         endif
         if(ihardlaw.eq.20) then

            if (NPH.gt.1)then
               write(*,*)'update_crss_dd_fe is not ready to be applied'
               write(*,*)'for multiphase polycrystals'
               call xit()
               call exit
               call xit()
            endif
            tempdum=273.
            call update_crss_dd_fe_new(
c     intent(in)
     $           0,ife,UR_SX,UW_DIS,UW_TAU,nphmx,nsysmx,nslmod,nslsys,
     $           nmodes,
     $           nsyst,nsm,ntwmod,tempdum,ngr,wgt,ngrmx,aux6,dtime,
     $           nmodmx,
     $           gamdot,NTWMMX,xmu_mode,
c     intent(inout)
c          read from (ioption.eq.0) and not evolving
     $           grsze,chifor,chisub,RHO_initW,TEMPref,burg,actener,
     $           K1gener,
     $           DRAG,EDOT0,rho_inits,
     $           Atau0,Atau1,Atau2,
     $           Fdeb0,Fdeb1,Fdeb2,
     $           HPcoef,T0nucl,T1nucl,T2nucl,T0prop,T1prop,T2prop,
     $           ctwsl,cts1,cts0,cts2,
     $           K3recov,K4recov,D0vac,Efvac,Emvac,
c         subject to 'evolution'
     $           rho_deb,rho_s_for,rho_s_deb,
c          constants
     $           edotref,
c     intent(out)
     $           itwinlaw,hard,twthres,isectw,nrs,crss,
     $           rhosat,tau_prop,tau_nucl,
c     intent(out) for stats of DD law
     $           rho_for_ph_mod,rho_deb_ph,rho_tot_ph,rho_deb_ph_mod,
     $           rho_recov,taue_wavg,rhonucl,pnucl
     $           )

         endif
         if (ihardlaw.eq.1) then
            call stop_debug(0,'MTS and DD laws are not available')
         endif

         if(ihardlaw.eq.22) then
            write(*,*)'ihardlaw with 22 is obsolete'
            call stop_debug(0,'** Exiting as ihardlaw=22 is obsolete')
c           CALL UPDATE_CRSS_DD_CE2 (0)
         endif
c         CLOSE(UNIT=UR_SX)

C *** READS INITIAL TEXTURE (BUNGE.ROE,KOCKS CONVENTIONS) FROM 'FILETEXT'
C *** INITIALIZES DEFORMATION GRADIENT AND ELLIPSOID FOR EACH PHASE.
C *** IF ISHAPE>2 READS INITIAL GRAIN AXES AND ORIENTATION FROM FILEAXES
C     AND INITIALIZES DEF GRADIENT & ELLIPSOID FOR INDIVIDUAL GRAINS.

         if(ishape(IPH).ge.3)then
            write(*,*) 'IPH:',IPH
            write(*,*) 'ishape(IPH):',ishape(IPH)
            write(*,*) 'fileaxes:',fileaxes
            OPEN(action='read',UNIT=UR_SHP,FILE=FILEAXES,STATUS='OLD')
         endif

         if (verbose) write(*,*)'entering data_grain'


c        obtain discrete orientation on the fly from mtex-COD data
         if (iodf.eq.1.or.iodf.eq.2) then !.and.ife.eq.0) then
            if (iodf.eq.1) then
               open(action='read',unit=UR_TX,file=filetext,status='old')
            endif
c           Caution!! number of grains for this phase is
c           hardwired as ngrmx for now. This works
c           for a single phase material, but cetainly
c           requires improvement... (2021-Dec 13 YJ)
            no_grain_per_phase = ngrmx / nphmx
            nrem = ngrmx - no_grain_per_phase * nphmx
            if (nrem.ne.0) then
               write(*,*)'** Unexpected case: Non zero remainder'
               call stop_debug(0,'****')
            endif

            ngr(IPH)=ngr(IPH-1)+no_grain_per_phase
            igb=ngr(IPH-1)+1
            igt=ngr(IPH)

            call data_grain_odf(
c                in
     $           ife,iodf,igt-igb+1,UR_TX,UR_SHP,
     $           ngrmx,nphmx,IPH,igb,igt,
     $           wph,ishape,eulerph,axisph,
c                out
     $           wgt,ag,fijph,fijgr)
            if (iodf.eq.1) close(UR_TX)
         elseif(iodf.eq.0) then
            OPEN(action='read',UNIT=UR_TX,FILE=FILETEXT,STATUS='OLD')
            do i=1,3
               read(UR_TX,'(a)')prosa
            enddo
            read(UR_TX,*) prosa, no_grain_per_phase
            rewind UR_TX
            ngr(IPH)=ngr(IPH-1)+no_grain_per_phase
            igb=ngr(IPH-1)+1
            igt=ngr(IPH)
            call data_grain_fe(UR_TX,UR_SHP,
     $           ngrmx,nphmx,IPH,igb,igt,wph,
     $           ishape,eulerph,axisph,
c             inout
     $           no_grain_per_phase,
c             out
     $           wgt,ag,fijph,fijgr)
         elseif(iodf.eq.3) then
c           Caution!! number of grains for this phase is
c           hardwired as ngrmx for now. This works
c           for a single phase material, but cetainly
c           requires improvement... (2024-Mar 07 YJ)
            no_grain_per_phase = ngrmx / nphmx
            ngr(IPH)=ngr(IPH-1) + no_grain_per_phase
            igb=ngr(IPH-1)+1
            igt=ngr(IPH)

            call data_grain_ebsd(UR_SHP,ngrmx,nphmx,IPH,igb,igt,wph,
     $           ishape,eulerph,axisph,coords,ncrds,ebsd_rad,
c             inout
     $           no_grain_per_phase,
c             out
     $           wgt,ag,fijph,fijgr)

         elseif(iodf.eq.5) then
c           discrete orientation but randomly choose.
            if (nphmx.ne.1.or.NPH.ne.1) then
               write(*,*)'** Error: iodf=5 case works for only'
               write(*,*)'   single-phase polycrystal'
               call xit()
            endif

            OPEN(action='read',UNIT=UR_TX,FILE=FILETEXT,STATUS='OLD')
            do i=1,3
               read(UR_TX,'(a)')prosa
            enddo
            read(UR_TX,*) prosa, no_grain_ori
            rewind UR_TX
            ngr(IPH)=ngr(IPH-1)+ngrmx
            igb=ngr(IPH-1)+1
            igt=ngr(IPH)

            if (no_grain_ori .lt. ngrmx * 10) then
               write(*,*)'** Error: iodf=5 case works only with'
               write(*,*)'   FNTX with a large # of ori.'
               write(*,*)'   At least 10 times larger discrete'
               write(*,*)'   orientations should be given'
               call xit()
            endif
c     ---   Test if the no of data lines is greater or equal to
c           the specified <no_grain_per_phase>
            do i=1,4
               read(UR_TX,'(a)')prosa
            enddo
            do i=1,no_grain_ori
               read(UR_TX,'(a)',err=100) prosa
            enddo
            goto 150
c           ---------------------------------------------
 100        write(*,*)'** Error: insufficient data lines'
            call xit()
c           ---------------------------------------------
 150        rewind UR_TX


            call data_grain_fe_shuffle(
c             in
     $           UR_TX,UR_SHP,ngrmx,nphmx,IPH,igb,igt,wph,
     $           ishape,eulerph,axisph,
c             inout
     $           no_grain_ori,
c             out
     $           wgt,ag,fijph,fijgr)
         endif

         if (iodf.ne.3) close(unit=UR_TX)

         if(ihardlaw.eq.23) then
            if (verbose) write(*,*)'initialize crss_rgvb_fe'
            igb=ngr(IPH-1)+1
            igt=ngr(IPH)
            call update_crss_rgvb_fe(
     $           UR_SX,0,IPH,igb,igt,dtime,nsysmx,ngrmx,nphmx,nmodmx,
     $           ntwmmx,
     $           nsyst,nrs,nslmod,ntwmod,nslsys,ntwsys,nmodes,
     $           nsm,igt-igb+1,dnca,jran,wgt,
c                inout
     $           TAUE,crss,gamdot,ishear_tag,ibs_law,xMu,Xb,xf_self,
     $           xf_lat,xD,xK,alpham,xlambda,tau0,xp,xbx,facs,pow,
     $           xRho_min,xRho_max,xRho_infl,rho_f,rho_r,taud,
     $           rho_prev,rho_deb_rgvb,tauini,rhos,rho_mod,fac,for_rhos,
     $           rev_rhos,cts,for_rhosl,rev_rhosl,AMSS,deltaD,xp_g,
     $           drho_tot,debs,rho_ps,tautwint,
c                out
     $           isectw,twthres,itwinlaw,dum1,dum2,dum3,dum4)
         endif


C *** IF RECRYSTALLIZATION IS PRESENT READ & INITIALIZE PARAMETERS AND
C     ARRAYS BUT ONLY FOR IPH=1 (PARENT)

        irx=0
        read(UR_SX,'(a)',end=626) prosa
        read(UR_SX,*,end=626,iostat=iostat) irx
        if (iostat.gt.0) then
           irx=0
           goto 626
        endif


        if(irx.eq.1) then
           read(ur_sx,'(a)')prosa
           read(ur_sx,*) (rexnuclpar(i),i=1,5)
           read(ur_sx,*) (rexgrowpar(i),i=1,5)

           do kkk=ngr(iph-1)+1,ngr(iph)
              radius(kkk)=grsze(iph)*0.5
           enddo

        endif
 626    close(unit=UR_SX)

        call mutexunlock(3)

         if (verbose) write(*,*)'update_shape_fe'

         if(ishape(IPH).ge.3) close(unit=UR_SHP)
c        INITIAL SETTING TO UPDATE THE SHAPE OF THE PHASE
c         iflat(IPH)=0
         call update_shape_fe(
c             in
     $        ngrmx,nphmx,IPH,ngr(IPH-1)+1,ngr(IPH),interaction,
     $        ifrag,crit_shp,ishape,
!     inout
     $        fijph,iflat,axisph,axisgr,
c     out
     $        eulerph,fijgr_temp,ierr)

         if(ngr(IPH).gt.ngrmx) then
            write(*,*)'ngr(IPH):',ngr(IPH)
            write(*,*)'ngrmx:',ngrmx
            WRITE(*,'(''-->INCREASE PARAMETER NGRMX'')')
            call xit()
         endif


         if (verbose) then
            write(*,*)'thermo_emod_cte (0 or -1)'
            write(*,*)'ife:',ife
         endif


c     temperature-dependent elastic modulus and coeffcient of thermal-expansion
         if (isx_ver.eq.1) then
            if (ife.eq.0) then
c               call w_chr(0,'calling thermo_emod_cte, with opt=1')
c               call w_val(0,'cijv_fit(1,1,1,1)',cijv_fit(1,1,1,1))
c               call w_ival(0,'IPH',IPH)
               call thermo_emod_cte(
c     intent(in)
     $              1,UR_ECTE,
     $              IPH,nsysmx,nmodmx,nphmx,nmodes,mxndeg,nsm,
     $              298d0,schca,
c     intent(inout) !! should be determined inside <evpsc_input_fe>
     $              ndegs_fit,cijv_fit,cte_fit,
     $              celcc,selcc,xmu_mode,athcc)
            elseif(ife.eq.1) then
               if (verbose) write(*,*)'thermo_emod_cte with iopt=1'
               call thermo_emod_cte(
c     intent(in)
     $              -1,UR_ECTE,
     $              IPH,nsysmx,nmodmx,nphmx,nmodes,mxndeg,nsm,
     $              298d0,schca,
c     intent(inout) !! should be determined inside <evpsc_input_fe>
     $              ndegs_fit,cijv_fit,cte_fit,
     $              celcc,selcc,xmu_mode,athcc)
            endif

            if (verbose) then
               write(*,*)'thermo_emod_cte1'
            endif
         endif
      enddo                     ! END OF DATA INPUT LOOP OVER ALL PHASES

C ****************************************************************************
      if (verbose) write(*,*)'READS SETTINGS FOR CONVRGN PROCEDURES.'
C *** READS SETTINGS FOR CONVERGENCE PROCEDURES.
      read(UR_INP,'(a)') prosa
      read(UR_INP,*) errs,errd,errm,errso
      read(UR_INP,*) itmaxext,itmaxint,itmaxso
      irsvar=0
      read(UR_INP,*) irsvar,jxrsini,jxrsfin,jxrstep
      read(UR_INP,*) ibcinv

C ****************************************************************************
      read(UR_INP,'(a)') prosa
      read(UR_INP,*) irecover
      read(UR_INP,*) isave
      read(UR_INP,*) icubcom
      read(UR_INP,*) nwrite
C ****************************************************************************
C *** READS MODELING CONDITIONS.
      read(UR_INP,'(a)') prosa
      read(UR_INP,*) interaction
      if (interaction.eq.5) then
         call stop_debug(0,'** 2nd order interaction is not available')
      endif
      read(UR_INP,*) iupdori,iupdshp,iupdhar
      read(UR_INP,*) nneigh
      if (nneigh.ne.0) then
         call stop_debug(0,'** Neighbouring is not available')
      endif
      iflu=0
      read(UR_INP,*) iflu
      if (nneigh.ne.0) then
         call stop_debug(0,'** Fluctuation formula is not avaialbe')
      endif
c$$$C *** CHECKS PRECISION CONTROL FOR SO CALCULATIONS
c$$$      if(interaction.eq.5) then
c$$$         if(errs.gt.1.e-7.or.errd.gt.1.e-7.or.errm.gt.1.e-7)then
c$$$            write(*,
c$$$     $       '('' SO CALCULATIONS MAY NOT CONVERGE IF TOLERANCES'',
c$$$     $'' ERRS,ERRD,ERRM ARE LARGER THAN 1.E-7'')')
c$$$            read(*,*)
c$$$         endif
c$$$         if(errso.lt.1.e-2) then
c$$$            write(*,
c$$$     $       '('' SO CALCULATIONS MAY NOT CONVERGE IF TOLERANCE'',
c$$$     $'' ERRSO IS SMALLER THAN 1.E-2'')')
c$$$            read(*,*)
c$$$         endif
c$$$      endif

C *** CHECKS IF RATE SENSITIVITY IS THE SAME FOR ALL SYSTEMS.
C *** SEARCH FOR NRSMIN (NEEDED TO GET TAUMAX INSIDE NR SUBROUTINE)
      nunique=1
      ncompa=nrs(1,1)
      nrsmin=nrs(1,1)
      if (verbose) call w_ival(0,'ncompa',ncompa)

      do IPH=1,NPH!iph0,iph1
         if (verbose) call w_ival(0,'IPH:',IPH)
      do is=1,nsyst(IPH)
         if (verbose) call w_ival(0,'nrs(is,IPH):',nrs(is,IPH))
         if(nrs(is,IPH).ne.ncompa) then
            nunique=0
            if (verbose) write(*,*)'occured for',IPH,is
         endif
         if(nrs(is,IPH).lt.nrsmin) nrsmin=nrs(is,IPH)
      enddo
         if (verbose) write(*,*)'----'
      enddo

      if(nunique.eq.0) then
         write(*,'(''-->THIS RUN CONTAINS MIXED nrs EXPONENTS !!'')')
         if(interaction.ne.1) then
            write(*,'(''-->WILL FORCE AFFINE (INTERACTION=1)'')')
            interaction=1
            read(*,*)
         endif
      endif

      if(interaction.eq.5.and.iflu.eq.0) then
         write(*,*) 'FLUCTUATIONS ARE NEEDED FOR 2nd-ORDER CALCULATION'
         write(*,*) 'IFLU IS RESET TO 1'
         iflu=1
      endif
      if(interaction.eq.0.and.iflu.eq.1) then
         write(*,*)
     $   'INTRAGRANULAR FLUCTS CANNOT BE OBTAINED IN FC CALCULATION'
         write(*,*) 'IFLU IS RESET TO 0'
         iflu=0
      endif
      if(iflu.eq.1) then
         open(unit=83,action='write',file='FLUCT.OUT',status='unknown')
         if(icubcom.eq.1)
     $        open(unit=84,action='write',file='FLCUB.OUT',
     $        status='unknown')
      endif
      if(interaction.eq.5) then
         open(unit=97,action='write',file='SO.OUT',status='unknown')
      endif
      if(iupdori.eq.0 .and. iupdshp.eq.1) then
         write(*,'('' --> IUPDSHP=1 AND IUPDORI=0'')')
         write(*,'('' --> CANNOT UPDATE SHAPE WITHOUT ALSO UPDATING'',
     $            '' ORIENTATION'')')
c         call xit()
      endif

      ICOMPGRAIN=0

C ****************************************************************************
      if (verbose) write(*,*)'INITIALIZE NEIGHBOURS'
C *** INITIALIZE NEIGHBOURS TABLE IF NNEIGH IS NOT ZERO.
C     IOPTION=0 --> RANDOM PAIRING / IOPTION=1 --> SEQUENTIAL PAIRING
      ioption=0
      call neighbours_fe(
c          in
     $     ioption,jran,nneigh,nneimx,ngrmx,ngr(NPH),wgt,
c          out
     $     neigh,wneigh)

C *******************************************************************
C *** INITIALIZE ARRAYS ASSOCIATED WITH GRAINS:
C     HARDENING, ACCUMULATED SHEAR, POWER, TWINNING PARAMETERS, etc
c      call stop_debug(0,'stopped inside evpsc_input_fe')
      if (verbose) write(*,*)'INITIALIZE CRSS of Voce, DD_FE'

      gamd0s(:,:)=1d0
      gtotgr(:)=0d0
      do IPH=1,NPH
C     *** INITIALIZES TRANSFORMATION RATE DG_TRANS(I,KKK)
C     INTRODUCED TO ACCOUNT FOR GROWTH
         if(ihardlaw.eq.0) then
            igb=ngr(IPH-1)+1
            igt=ngr(IPH)
            call update_crss_voce_fe(
c                in
     $           UR_SX,1,IPH,nphmx,ngrmx,dtime,nmodmx,ntwmmx,igb,igt,
     $           nmodes,nsysmx,nsm,nslsys,nsyst,gamdot,
c                inout
     $           grsze,gtotgr,tau_voce,thet_voce,
     $           hpfac_voce,hard,crss,
c                out
     $           twthres,nrs,itwinlaw,isectw)
         elseif(ihardlaw.eq.20) then
            aux6(:)=0. !! there is no deformation yet.
            call update_crss_dd_fe_new(
c     intent(in)
     $           1,0,UR_SX,UW_DIS,UW_TAU,nphmx,nsysmx,
     $           nslmod,nslsys,nmodes,
     $           nsyst,nsm,ntwmod,tempdum,ngr,wgt,ngrmx,aux6,dtime,
     $           nmodmx,
     $           gamdot,NTWMMX,xmu_mode,
c     intent(inout)
c          read from (ioption.eq.0) and not evolving
     $           grsze,chifor,chisub,RHO_initW,TEMPref,burg,actener,
     $           K1gener,
     $           DRAG,EDOT0,rho_inits,
     $           Atau0,Atau1,Atau2,
     $           Fdeb0,Fdeb1,Fdeb2,
     $           HPcoef,T0nucl,T1nucl,T2nucl,T0prop,T1prop,T2prop,
     $           ctwsl,cts1,cts0,cts2,
     $           K3recov,K4recov,D0vac,Efvac,Emvac,
c         subject to 'evolution'
     $           rho_deb,rho_s_for,rho_s_deb,
c          constants
     $           edotref,
c     intent(out)
     $           itwinlaw,hard,twthres,isectw,nrs,crss,
     $           rhosat,tau_prop,tau_nucl,
c     intent(out) for stats of DD law

     $           rho_for_ph_mod,rho_deb_ph,rho_tot_ph,rho_deb_ph_mod,
     $           rho_recov,taue_wavg,rhonucl,pnucl
     $           )

         elseif(ihardlaw.eq.22) then
            write(*,*)'ihardlaw with 22 is obsolete'
            call exit(-1)
         elseif(ihardlaw.eq.23) then
            igb=ngr(IPH-1)+1
            igt=ngr(IPH)
            call update_crss_rgvb_fe(UR_SX,1,IPH,igb,igt,dtime,nsysmx,
     $           ngrmx,nphmx,nmodmx,
     $           ntwmmx,
     $           nsyst,nrs,nslmod,ntwmod,nslsys,
     $           ntwsys,nmodes,nsm,igt-igb+1,dnca,jran,wgt,
c                inout
     $           TAUE,crss,gamdot,ishear_tag,ibs_law,xMu,Xb,xf_self,
     $           xf_lat,xD,xK,alpham,xlambda,tau0,xp,xbx,facs,
     $           pow,xRho_min,xRho_max,xRho_infl,rho_f,rho_r,taud,
     $           rho_prev,rho_deb_rgvb,tauini,rhos,rho_mod,fac,for_rhos,
     $           rev_rhos,cts,for_rhosl,rev_rhosl,AMSS,deltaD,xp_g,
     $           drho_tot,debs,rho_ps,tautwint,
c                out
     $           isectw,twthres,itwinlaw,dum1,dum2,dum3,dum4)
         endif
      enddo

      return
      end subroutine evpsc_input_fe
