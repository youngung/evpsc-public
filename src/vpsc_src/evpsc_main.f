C ****************************************************************************
C                  EVPSC 1-site multiphase                                   "
C ****************************************************************************
C                 Original VPSC code BY R.LEBENSOHN & C.N.TOME               *
c                  has been modified by Y. JEONG & C.N.TOME                  *
c                         to account for elasticity                          *
c      The code has been extended for FE interface via UMAT of               *
c     Abaqus/standard solver by Y. JEONG, B. JEON, C.N.TOME                  *
******************************************************************************
C                   MULTIPHASE OR MULTIELEMENT VERSION                       *
C ****************************************************************************
C  RICARDO LEBENSOHN                 |  CARLOS TOME                          *
C  MST-8 - LANL - MS G755            |  MST-8 - LANL - MS G755               *
C  LOS ALAMOS - NM 87545 - USA       |  LOS ALAMOS - NM 87545 - USA          *
C  lebenso@lanl.gov                  |  tome@lanl.gov                        *
C ****************************************************************************
c  YOUNGUNG JEONG
c  MSE, Changwon National University
c  Changwon, Republic of Korea
c  yjeong@changwon.ac.kr
c****************************************************************************
C  GENERAL REFERENCE for VPSC                                                *
C  R.A.Lebensohn & C.N.Tome, Acta metall mater 41, 2611 (1993)               *
C  - GENERAL REFERENCE for dEVPSC                                            *
C  Y. Jeong & C.N.Tome, Int. J. Plast. (2020)                                *
C  - GENERAL REFERENCE for dEVPSC-FE                                         *
C  Y. Jeong, B. Jeon, C.N.Tome, Int. J. Plast. (2021)                        *
c  - Extension to Forming limit diagram calculation                          *
c  Y. Jeong et al, MSMSE (2016)                                              *
c  - sigma-EVPSC                                                             *
c  B. Jeon et al, JMRT (2024)                                                *
c----------------------------------------------------------------------------*
C  COPYRIGHT NOTICE                                                          *
C ****************************************************************************
C Copyright (c) 2020-2026, Youngung Jeong and EVPSC contributors             *
C Copyright (c) 2022, Los Alamos National Laboratory                         *
C All rights reserved.                                                       *
C                                                                            *
C Portions of this program are derived from VPSC (c) 2022 Triad National     *
C Security, LLC, released under the BSD 3-Clause License.                    *
C This program is open source under the BSD 3-Clause License.                *
C See the LICENSE and NOTICE files distributed with this program for the     *
C full license text and the U.S. Government rights notice.                   *
C ****************************************************************************
C                              # FEATURES #                                  *
C ****************************************************************************
C  - SINGLE PHASE & MULTIPHASE AGGREGATES.                                   *
C  - SC RELATION FOR MIXED RATE-SENSITIVITY EXPONENTS.                       *
C  - ESHELBY AND PRESSURE TENSORS FOR FULLY INCOMPRESSIBLE INCLUSION & HEM.  *
C  - MIXED BOUNDARY CONDITIONS ON DISPLACEMENT GRADIENT AND STRESS           *
C  - SYMMETRIC BASIS INTERNAL REPRESENTATION OF TENSORS                      *
C  - OPTION FOR ANALYZING ROLLING FCC COMPONENTS (icubcom=1)                 *
C  - SELF & LATENT HARDENING COEFFICIENTS                                    *
C  - 'VOCE' TYPE HARDENING MODULATION FOR EACH SYSTEM (ihardlaw=0) OR        *
C  - 'VOCE-thermo' TYPE HARDENING MODULATION FOR EACH SYSTEM (ihardlaw=3)    *
C    'MTS'  HARDENING LAW FOR EACH SYSTEM (ihardlaw=1)                       *
C  - PREDOMINANT TWIN REORIENTATION SCHEME (if ihardlaw=0)                   *
C  - COMPOSITE GRAIN SCHEME (if ihardlaw=2)                                  *
C  - DISLOCATION DENSITY EVOLUTION (if ihardlaw=20-29)                       *
C  - DIFFERENT INDIVIDUAL GRAIN SHAPES.                                      *
C  - GENERAL TREATMENT OF GRAIN SHAPE EVOLUTION.                             *
C  - INCREMENTAL CONTROL: VON MISES STRAIN, SINGLE STRAIN COMPONENT OR TIME  *
C  - POSTMORTEM CAPABILITY: RECOVERS ARRAYS FROM PREVIOUS RUN                *
C  - POLYCRYSTAL YIELD SURFACE AND LANKFORD COEFFICIENTS CALCULATION         *
C  - VARIABLE DEFORMATION HISTORY ROUTINE.                                   *
C  - OPTION FOR ROTATION COUPLING BETWEEN ORIENTATIONS                       *
C  - ALL 'CANNED' SUBROUTINES ARE FROM 'NUMERICAL RECIPIES' LIBRARY          *
C  - INTERNAL Bunge CONVENTION (phi1,Phi,phi2 -> phi,theta,omega) FOR        *
C    EULER ANGLES. BUT READS ALSO Kocks and Roe ANGLES AS INPUT.             *
C                                                                            *
C ****************************************************************************
C                         # INPUT FILES #                                    *
C ****************************************************************************
C ALWAYS NEEDED:                                                             *
C                                                                            *
C       VPSC7.IN      - general input and parameters                         *
C       "FILECRYS"    - deformation modes and CRSS (one per phase)           *
C       "FILETEXT"    - initial texture file (one per phase)                 *
C       "FILEPROC"    - deformation process to be simulated                  *
C                                                                            *
C SOMETIMES NEEDED:                                                          *
C                                                                            *
C  if ishape>1     "FILEAXES"  - initial orientation of ellipsoids           *
C                               (morphologic texture, one per phase)         *
C  if irecover=1   POSTMORT.IN - initial state from previous run             *
C  if icubcom=1    CUBCOMP.IN  - ideal orientations - fcc rolling            *
C  if ivgvar=1     "FILEHIST"  - strain history                              *
C  if ihardlaw=1   MTS parameters for the material at the end of "FILECRYS"  *
C  if ihardlaw=2X  DD model parameters at the end of "FILECRYS"              *
C                                                                            *
C ****************************************************************************
C                         # OUTPUT FILES #                                   *
C ****************************************************************************
C DEPENDING ON CASE AND DRIVER:                                              *
C        RERR.OUT     - convergence history                                  *
C        STATS.OUT    - statistics                                           *
C        ACT_PHn.OUT  - activity of deformation modes                        *
C        TEX_PHn.OUT  - cryst. texture for each phase/elem at dif steps      *
C        MOR_PHn.OUT  - morph. texture for each phase/elem at dif steps      *
C        STR_STR.OUT  - macroscopic stress-strain components                 *
C        RUN_LOG.OUT  - a log containing the conditions and input to the run *
C        STAT_AX.OUT  - statistics on grain axes                             *
C                                                                            *
C     if isave=1    POSTMORT.OUT - final state in grains and PX              *
C     if icauchy=1  CAUCHY.OUT - final local cauchy stresses                 *
C     if icubcom=1  CUBCOMPn.OUT - rolling components for fcc                *
C                                                                            *
C                                                                            *
C*****************************************************************************

      PROGRAM evpsc_main
      implicit none
c**   arrays from common blocks shared with incr. evpsc-umat
      include "incr_evpsc_fe.dim"
c     -----------------------------------------------------------------c
c**   controlling integer/flags
      integer ife,iverbose
      parameter(ife=0)          !0: SA 1: FE
      parameter(iverbose=1)
      integer iproc,istep
c     -----------------------------------------------------------------c
c**   Dummy constants/variables and auxilary arrays
      integer m,n,iaux33(3,3),iaux6(6),ibux6(6)
      real*8 dum,dum1,dum2,dum3,dum4,pnewdt,dtime_ref,
     $     aux2(2),aux3(3),aux33(3,3),bux33(3,3),aux3333(3,3,3,3),
     $     aux6(6),bux6(6),cux6(6),aux330(3,3,0),aux66(6,6)
c     -----------------------------------------------------------------c
      CHARACTER ITERLBL*40, fn_snap*120
      real*8 sbarc(3,3),dsbarc(3,3),dbarc(3,3),tempdot
c     -----------------------------------------------------------------c
      real*8 dtime,time
      integer i,j,iph,kkk,indx,indy,idum,nstepsx,igb,igt,it2,
     $     itextwrite
c     -----------------------------------------------------------------c
c**   diagonotics flags and values.
      real*8 rels,reld,relsgr
      integer idiag_sg,idiag_sc,ierr_shp,idiag_nr
c**   misc
      real*8 pi
      parameter(pi=3.14159265359)
      integer ijv(6,2)
      DATA ((IJV(N,M),M=1,2),N=1,6)/1,1,2,2,3,3,2,3,1,3,1,2/
c     inter istatus(1)
      real*8 time_begin,time_end ! to measure the elapsed running time
      integer UR_INP,uw_stat,uw_alls,uw_st01,uw_elsc,uw_cntr,uw_crav,
     $     uw_esgr,UW_SCCN,uw_phpr,UW_STR,URW_FLDA,UW_FLDB,
     $     uw_act,uw_gam,uw_tex,uw_dif,uw_ssch,uw_void,
     $     UW_REX,UW_REX_STAT,UW_REX_NUCO,UW_REX_GROW,
     $     UW_SRX_TEX
c     functions defined
      real*8 vnorm              ! src/vpsc_src/funcs.f
c      real*8 anwork,aneps(6),anepsvm
c     -----------------------------------------------------------------c
c     FLD simulation associated arrays.
      real*8 LIJBARc_a_lab(3,3),epstot_a_lab(3,3),epstot_b_lab(3,3),
     $     sbar_old_b_lab(6),sbar_a_lab(6),sbar_b_lab(6),rho_fld,
     $     cvm_fld,evm_limit,psi0_fld,psi_fld,ftilde0,ftilde,t0_fld(3),
     $     R_grvlab(3,3),R_labgrv(3,3)
      character*258 fn_region_a
c     -----------------------------------------------------------------c
      real*8 total_unload_time
      logical ipx_loaded_ever,fld_limit,print_loc
      parameter(print_loc=.false.)
c      parameter(print_loc=.true.)
c     command line argument
c     -----------------------------------------------------------------c
      integer iodf,ielast,ioptfs
c     -----------------------------------------------------------------c
c     0: from discrete population, 1: from ODF (mtex); 2: random isotropic texture; 3: EBSD; 4: Euichan; 5: shuffle & choose
      parameter(iodf=0)
c     -----------------------------------------------------------------c
c     0: elasto-visco-plasticity;  1: suppressing plasticity by assigning CRSS = 1e6
      parameter(ielast=0)
c     -----------------------------------------------------------------c
c     0: Using the linear form of Jeong et al. (IJP, 2021) and Jeon et al. (JMRT, 2024)
c     1: Using the new linear form of Patra et al.
      parameter(ioptfs=0)
c     --
      integer iopt
c     --- interfacing statev array used in Abaqus for EVPSC-SA
c      integer nstatv,nvar
c      parameter(nstatv=23300000)
c     parameter(nstatv=31601990)
c      parameter(nstatv=31602902)
c      real*8 statev(nstatv)
c     ----
      real*8 xnumstats(7)
c     --- sub trials
      logical crss_init,isnan,iesh_fail
c     --------------------------
c     date time related variables
      character(8)::date_
      character(10)::time_
      character(5)::zone_
      integer :: itime_val(8)
      real*8 epstot_ph(3,3,nphmx)
c     contribution
c     temperature-dependent CTE, e-mod
      integer isx_ver
      real*8 contrib(3)
c----------------------------------------------------------------------c
      real*8 rhosat(nmodmx),TAU_PROP(nmodmx),TAU_NUCL(nmodmx)
c     quantities for stats analysis from update_crss_dd_fe_new
      real*8 rho_for_ph_mod(nmodmx,nphmx),rho_deb_ph(nphmx),
     $     rho_tot_ph(nphmx),rho_deb_ph_mod(nmodmx,nphmx),
     $     rho_recov,taue_wavg(nphmx,6,4),RHOnucl,pnucl
c     quantities for fs stats
      real*8 fel(ngrmx),fvp(ngrmx)
c----------------------------------------------------------------------c
c     for the case EBSD is directly read from .ang file (iodf.eq.3)
      integer ncrds,nlines_ebsd,nrowmx,mdim
      parameter(ncrds=3)
      real*4 f_scale
      real*8 coords(ncrds),ebsd_rad
      parameter(nrowmx=200000,ebsd_rad=30,f_scale=real(0.1,4))
      real*4 ebsd_raw_data(5,nrowmx)
      common/kebsd/ebsd_raw_data,nlines_ebsd
c----------------------------------------------------------------------c
      real*8 GAVMODG(nmodmx,ngrmx)
c----------------------------------------------------------------------c
c     void-related
      real*8 wmatl,wvoid,wgtsum,dph(6,2),fdot
c----------------------------------------------------------------------c
c     recrystallization
      integer dngr_rex
      real*8 time_srx
c----------------------------------------------------------------------c
c     other misc. variables/constants/arrays
      real*8 edot_el(3,3)

      call RANDOM_INIT(repeatable=.true., image_distinct=.true.)
c      call random_number(dum) ! to shift the psuedo-random table
c      call random_number(dum) ! to shift the psuedo-random table
c      call random_number(dum) ! to shift the psuedo-random table
c      call RANDOM_INIT(repeatable=.false., image_distinct=.true.)

      is_initialized_statev=0
c----------------------------------------------------------------------c
c     This coordinates were introduced by Jaeseong to test EBSD-based
c     texture option.
      coords(:)=0.
      coords(1)=10.
      coords(2)=15.
c----------------------------------------------------------------------c
      call cpu_time(time_begin)
      call date_and_time(date_,time_,zone_,itime_val)
      write(*,*)'--------------------------'
      write(*,*)'EVPSC numerical simulation'
      write(*,*)'--------------------------'
      write(*,'(3x,a8,2x,a7)')'date','time'
      write(*,'(3x,a8,2x,2(i2.2,a1),i2.2,a1,i3.3)')
     $     date_,itime_val(5),':',itime_val(6),':',itime_val(7),'.',
     $     itime_val(8)
      write(*,'(a8,i1)')'iodf:',iodf
      write(*,*)'--------------------------'

c     initialize Eshelby related quantities
      CALL ESHELBY(AUX3,AUX3333,0.d0,AUX3333,AUX3333,AUX33,AUX33,dum,
     $     AUX3333,AUX3333,0,iesh_fail)

      mdim=6
      if     (iregime.eq.1) then
         mdim=5
      elseif (iregime.eq.2) then
         write(*,*)'** Delta-EVPSC calculation'
      elseif (iregime.eq.3) then
         write(*,*)'** Sigma-EVPSC is deprecated.'
         call xit()
      else
         write(*,*)'** Error: Unexpected iregime given:',iregime
         call xit()
      endif

C ***********************************************************************
C     ASSIGNS # TO UNITS AND OPENS I/O FILES.
      UR_SX=  1                 ! FILECRYS      (OPEN/READ/CLOSE INSIDE EVPSC_INPUT_FE)
      UR_TX=  2                 ! FILETEXT      (OPEN/READ/CLOSE INSIDE EVPSC_INPUT_FE)
      UR_SHP= 3                 ! FILEAXES      (OPEN/READ/CLOSE INSIDE EVPSC_INPUT_FE)
      UR_STV= 4                 ! statev_px   (OPEN/CLOSE IN EVPSC_MAIN)
      UR_DIF= 7                 ! FILEDIFF(nph)
      UR_INP= 9                 ! VPSC7.IN      (OPEN/CLOSE IN EVPSC_MAIN)
      UW_LOG=10                 ! RUN_LOG.OUT   (OPEN IN EVPSC_MAIN)
      UW_ERR=12                 ! RERR.OUT
      UW_STR=13                 ! STR_STR.OUT   (OPEN IN EVPSC_MAIN)

      UW_REX=14                 ! REXDD.OUT       (OPEN IN REXDD_MAIN)
      UW_REX_STAT=15            ! REXDD_STATS.OUT (OPEN IN REXDD_MAIN)
      UW_REX_NUCO=16            ! REXDD_NUCL_ORI.OUT (OPEN IN REXDD_MAIN)
      UW_REX_GROW=17            ! REXDD_GROW_xxx_yyyy.OUT (OPEN/CLOSE in REXDD_MAIN)
      UW_SRX_TEX =18            ! REXDD_SRX_TEX_PH1.OUT (OPEN/CLOSE in EVPSC_MAIN)
C     TEX_PHX.OUT (open in <EVPSC_MAIN>)
      UW_TEX=25


C     ACT_PHX.OUT (OPEN IN <WRITE_ACTIV_FE.F> AND CLOSE IN EVPSC_MAIN)

      UW_ACT=30                 ! ACT_PH1.OUT
C            31                 ! ACT_PH2.OUT
C            32                 ! ACT_PH3.OUT
C            33                 ! ACT_PH4.OUT
      UW_GAM=40                 ! GAM_PH1.OUT
C            41                 ! GAM_PH2.OUT
C            43                 ! GAM_PH3.OUT
C            44                 ! GAM_PH4.OUT
      UW_DIS=45                 ! DIS_DEN.OUT
      UW_TAU=46                 ! TAUE.OUT
      UW_DIF=50                 ! DIF_STR_PH1
C            51                 ! DIF_STR_PH2
C            52                 ! DIF_STR_PH3
C            53                 ! DIF_STR_PH4
C            54                 ! DIF_WGT_PH1
C            55                 ! DIF_WGT_PH2
C            56                 ! DIF_WGT_PH3
C            57                 ! DIF_WGT_PH4
      UR_CUB=60                 ! CUBCOMP.IN    (OPEN/CLOSE IN EVPSC_MAIN)
      UR_HIS=61                 ! FILEHIST      (OPEN/CLOSE IN EVPSC_MAIN)
      UR_EBSD=62                ! .ANG EBSD FILE
C     BELOW OPEN/CLOSE IN FLD.F
      URW_FLDA=70
      UW_FLDB =71
c
      UR_ECTE=72                ! FN_TEMP_EMOD_CTE
      UW_SCCN=3012              ! SC_CONVERGENCE.OUT
      UW_CRAV=3015              ! CRSS_WAVG.OUT
      UW_ST01=3020              ! GRAIN_STATS_NSTD.OUT
      UW_STAT=3021              ! GRAIN_STATS.OUT
      UW_ALLS=3022              ! GRAIN_ALL.OUT
      UW_PHPR=3030              ! PH_STR?.OUT
      UW_ELSC=3042              ! ELSC.OUT - FOR YOUNG MODULUS (1/S1111)
      UW_CNTR=3053              ! CONTRIB.OUT
      UW_SSCH=3060              ! SCH_STATS.OUT + IPH
      UW_ESGR=5000              !
      UW_VOID=7000              ! VOID.OUT


c-----------------------------------------------------------------------
C     CALL SUBROUTINE VPSC_INPUT FOR READING CRYSTAL, GRAIN AND TEXTURE DATA.
C     INITIALIZES ARRAYS.
      OPEN(UR_INP,FILE='vpsc7.in',STATUS='OLD')
      OPEN(UW_LOG,FILE='RUN_LOG.OUT',STATUS='UNKNOWN')
      if (print_loc) write(*,*)"** calling evpsc_input_fe"
      include "call_evpsc_input_fe.f"

      if (NPH.eq.2) then
c        initial volid volume fraction
         wmatl=wph(1)
         ! wvoid=wph(2)
         open(uw_void,file='void.out',status='unknown')
         write(uw_void,'(3(a11,x))')'wph1','wph2','fdot'
      endif

C-----------------------------------------------------------------------
c***  Statev feature in evpsc_main is suppressed due to too large
c***  arrays are demanded when running DRX simulations with grain clusters.

C     Testing statev array with iopt=0
c$$$      if (print_loc) write(*,*)"** begins statev_io with iopt=0"
c$$$      iopt=0
c$$$      include 'call_statev_io.f'
c$$$      if (print_loc) write(*,*)"** ends statev_io with iopt=0"

C-----------------------------------------------------------------------
C     EBSD data-based texture studied by Jaeseong
      if(iodf.eq.3) then
         if (print_loc) write(*,*)"** iodf==3"
         OPEN(UR_EBSD, FILE='C:/Users/Jlee/'//
     $        'repo/umat/evpsc/matData/expData/'//
     $        'Mg10Gd/Mg10Gd_cylinder_d0_crop.ang',
     $        STATUS='UNKNOWN')
         call read_ebsd(UR_EBSD,f_scale,ebsd_raw_data,nlines_ebsd)
      endif

      OPEN(UW_ERR,FILE='RERR.OUT'   ,STATUS='UNKNOWN')
      open(uw_cntr,action='write',file='contrib.out',status='UNKNOWN')
      write(uw_cntr,'(4a14)') 'elast','vplast','therm','tot'
      call write_elsc(uw_elsc,0,aux66,aux66)

      open(UW_SCCN,file='sc_convergence.out',status='unknown')
      write(UW_SCCN,'(a3,x,4a12,2(6(2x,a4,i1,5x)),a12)')'IT2','RELSGR',
     $     'RELD','RELS','DTIME',('RELD',i,i=1,6),('RELS',i,i=1,6),
     $     'alplpx'
      call write_crav(0,uw_crav,NPH,ngr,nmodes,nphmx,ngrmx,
     $     nmodmx,nsm,nsysmx,wgt,crss)

C-----------------------------------------------------------------------
C     Initialize grain_stat files "grain_stats.out",
C     "grain_stats_all.out". This feature is optional and can be
C     suppressed. However, dE^{ref}, reference strain increment
C     might be necessary depending on dT^{deriv} scheme in use.
      if (print_loc) write(*,*)"** initialize grain_stat"
      aux6(:)=0.
      bux6(:)=0.
      dtime=0.
      idum=0
      istep=0
      call grain_stat(idum,istep,idum,uw_stat,uw_st01,uw_alls,ngr,nphmx,
     $     NPH,ngrmx,dtime,dgp,tg,dsg,sg_old,wgt,aux6,bux6,
     $     gtotgr)

C-----------------------------------------------------------------------
C     Texture and grain energy outputs (TEX_PHx.OUT, esgr_x.out)
      do iph=1,NPH
         open(UW_TEX-1+iph,file='TEX_PH'//CHAR(48+iph)//'.OUT',
     $        status='unknown')
         open(uw_esgr+iph,file='esgr_'//char(48+iph)//'.out',
     $        status='unknown')
         if (NPH.gt.1) then
            igb=ngr(iph-1) +1
            igt=ngr(iph)
            call phase_property_evp_fe(0,uw_phpr+iph-1,iph,
     $           igb,igt,ngrmx,nphmx,sg,tg,wgt,dtime,epstot_ph)
         endif
      enddo
C-----------------------------------------------------------------------
C     In case VFT scheme is used for twinning, a file that contains
C     both parent and child grains 'TEX_ALL.OUT' is created.
      if (itwinlaw.eq.3.and.NPH.eq.2) then
         write(*,*)'TEX_ALL.OUT for grains'
         open(31+2,file='TEX_ALL.OUT',status='unknown')
      endif
      if(interaction.ne.1.and.iregime.ne.1) then
         call stop_debug(0,'** Error: Only affine '//
     $        'linearization is valid for VPSC')
         call xit()
      endif
      if(interaction.eq.5.and.iregime.eq.2) then
         write(*,*)'interaction=5 cannot be used for dEVPSC'
         call xit()
      endif

      IF(IHARDLAW.EQ.1) THEN
         write(*,*)'** Error: unexpected case of mts crss law;'
         call xit()
      elseIF(IHARDLAW.EQ.22) then
         write(*,*)'** Error: ihardlaw with 22 is obsolete'
         call xit()
      endif

c-----------------------------------------------------------------------
c     Hard-wiring CRSS to suppress plasticity to emulate TESC simulations.
c     This is mainly for a particular case studied for Uranium's
c     thermo-mechanical simulations where we wanted to see how plasticity
c     contributes to the stress relaxation and thermal ratchetting.
      if (ielast.eq.1) crss(:,:)=1e9

c-----------------------------------------------------------------------
c     Initialize diffraction simulation procedure. ioption=0 reads
c     diffraction file for each phase and reads diffraction directions.
      if (print_loc) write(*,*)"** calling DIFF_PLANES_PX"
      call DIFF_PLANES_PX(
c     intent(in)::
     $     0,idiff,UR_DIF,UW_DIF,0,NPH,ngrmx,nphmx,ndifmx,ngr,ag,
     $     wgt,filediff,
     $     ICRYSYMph,iregime,iDBARv,aux33,aux6,BC_EL,mgel,aux33,wph,
c     intent(inout)::
     $     ndiff,sg,angdetector,nfamily,diff_vs,diff_vc,NGRSETPH,
     $     WGTSETPH,IGRSET,
c     intent(out)::
     $     eelgr)

c-----------------------------------------------------------------------
c     other various initialization of arrays / variables.
      time=0.
      temp=298.                 ! unless later specified
      tempdot=0.
      alplpx=0.
      EPSVM =0.
      epstot(:,:) =0.
      epstot_el(:,:)=0.
      epstot_ph(:,:,:)=0.
      plwork=0.
      DBAR0(:)=0.
      DBAR(:)=0.
      sbar(:)=0.

      ngr(0)=0
      ncluster=0                ! will be initialized in <rexdd>
      eulerph(:,0)=eulerph(:,0)
      sg_old(:,:)=0.
      sg(:,:)=0.
      dsg(:,:)=0.
      sbar_old(:)=0.
      sbar_old_b_lab(:)=0.

      crss_init=.false.
      fld_limit=.false.

      fel(:)=1.
      fvp(:)=0.
      gtotgr(:)=0.
      dg_trans(:,:)=0.
      dge(:,:)=0.
      dgp(:,:)=0.
c     grsze(:)=0. ! this might be determined in <evpsc_input_fe.f> for the case of DD hardening

c     ptr twin stat
      var_vf(:,:,:)=0.
      sch_stat(:,:,:)=0
      var_stat(:,:,:)=0
c     gavacc init
      gavacc(:,:)=0.
c     eigen strain-rate due to thermal expansion/contraction


      evm_limit=0.
      PRITW(:)=0.
      SECTW(:)=0.
c      EGVM(:)=0.
      WORK(:)=0.
c     starts treating parent & child independently
      ICOMPGR(:)=0
      KTWSMX(:)=0
      NTWEVENTS(:)=0
      eftwfr(:,:)=0.
      twfrph(:,:)=0.
      twfrsy(:,:)=0.
c     rextal
      dngr_rex=0

C ***********************************************************************
c     calculate elasto-visco-plastic compliance tensors
      if (print_loc) write(*,*)"** ELSC_EVP"
      iopt=0
      include "call_elsc_evp_fe.f"
      iopt=1
      include "call_elsc_evp_fe.f"


C ***********************************************************************
C     READS SX & PX STATE FROM PREVIOUS RUN (sav,xmsec,sg,ag,crss,..)
C     WHEN IRECOVER=1.
      IF(IRECOVER.EQ.1) THEN
         call stop_debug(0,'postmortem is not working for this version')
      ENDIF
      IF(ICUBCOM.EQ.1) THEN
         write(*,*)'icubcom is not available in dEVPSC'
         call xit()
      ENDIF
C ***********************************************************************
C     OPEN OUTPUT FILES
C ***********************************************************************
      OPEN(UW_STR,file='STR_STR.OUT',status='UNKNOWN')
c     writing only the column heads
      if (print_loc) write(*,*)"** INITIALIZE WRITE_STRESS_STRAIN"
      CALL WRITE_STRESS_STRAIN_evp(UW_STR,0,dtime,epsvm,svm,epstot,
     $     temp,LIJBARc,scauchy,press,presswgt,plwork,epstot_el,edot_el)

      DO I=1,NPH
        IF(ISHAPE(I).NE.0) THEN
          OPEN(16,FILE='STAT_AX.OUT',STATUS='UNKNOWN')
          OPEN(40+I,FILE='MOR_PH'//CHAR(48+I)//'.OUT',
     $         STATUS='UNKNOWN')
        ENDIF
      ENDDO


c     OPEN ACT_PHx.OUT and GAM_PHx.OUT
      do iph=1,NPH
         call write_activ_fe(0,uw_ACT+iph,uw_GAM+iph,iph,ihardlaw,
     $        nmodes(iph),ntwmod(iph),wph(iph),EPSVM,epstot,
     $        twfrph(1:ntwmod(iph),iph),eftwfr(1:ntwmod(iph),iph),
     $        pritw(iph),sectw(iph),activ(iph),
     $        gavmod(1:nmodes(iph),iph),gavacc(1:nmodes(iph),iph),
     $        GAVPH(IPH))
      enddo

C *******************************************************************
C     DO LOOP OVER PROCESSES
C *******************************************************************
C *** READS NUMBER OF PROCESSES TO RUN SEQUENTIALLY.

C     List of processes (IVGVAR)
C     0 : Uniform loading condition with history file
C     1 : Variable loading condition where vel. grad.
C     2 :
C     3 :
C     4 :
C     5 : Store EVPSC arrays to STATEV file

C     9 : SRX (should be followed from a DRX simulation?)

C     30: FLD for region A
C     31: FLD for region B

C     -4: In-plane rotation
C     -6: Unloading (applicable only to EVPSC)
C     -7: Stress hold

      READ(UR_INP,'(A)') PROSA
      READ(UR_INP,*)     NPROC
      READ(UR_INP,'(A)') PROSA

      if (print_loc)write(*,*)"** total number of processes:", NPROC

      ipx_loaded_ever=.false.
      DO 3000 IPROC=1,NPROC
         write(*,*)
         write(*,'(a8,x,i3.3)')'** iproc:',iproc
C *** READS LOAD CONDITIONS: STRESS, VELOCITY GRADIENT, TEMP, INCR, NSTEPS.
C *** CALCULATES SYMMETRIC STRAIN RATE 'DSIM' AND 'DBAR'
         READ(UR_INP,*) IVGVAR
c     ** Begins reading process information **
         IF(IVGVAR.EQ.0) THEN
            READ(UR_INP,'(A)') FILEHIST
            OPEN(UR_HIS,file=FILEHIST,status='OLD')
            CALL LOAD_CONDITIONS(
c                intent(in)::
     $           iregime,UR_HIS,interaction,
c                intent(out)::
     $           ictrl,nsteps,iDBARv,iSBARv,DBAR,dsim,ILBAR,LIJBARc,
     $           dsbar,sbar,sbarc,temp,tempdot,eijincr,dsimctrl,
     $           strain_control,dtime)
            CLOSE(UNIT=UR_HIS)
            nstepsx=nsteps
c     for uranium, use temperature from thermo_mech_conditions
c     to initiate CRSS
            if (.not.crss_init.and.ielast.eq.1) crss_init=.True.
            if (.not.crss_init.and.IHARDLAW.EQ.0.and.ielast.eq.0) then
               crss_init=.True.
               DO IPH=1,NPH
                  igb=ngr(iph-1)+1
                  igt=ngr(iph)
                  CALL UPDATE_CRSS_VOCE_fe(
c                   intent(in)::
     $                 UR_SX,1,iph,nphmx,ngrmx,dtime,nmodmx,ntwmmx,
     $                 igb,igt,nmodes,nsysmx,nsm,nslsys,nsyst,gamdot,
c                   intent(inout)::
     $                 grsze,gtotgr,tau_voce,thet_voce,hpfac_voce,hard,
     $                 crss,
c                   intent(out)::  (but valid only for ioption.eq.0)
     $                 twthres,nrs,itwinlaw,isectw)
               enddo
            endif
c     assign dg_trans for each grains based on tempdot and athcc
            dg_trans(:,:)=0.
            if (abs(tempdot).ge.1e-20) then
               do iph=1,NPH
                  cux6(:)=tempdot*athcc(:,iph)
                  DO KKK=NGR(IPH-1)+1,NGR(IPH)
                     aux33(:,:)=ag(:,:,kkk)
                     call mat6rot(cux6,aux33,dg_trans(:,kkk))
                  enddo
               enddo
            endif
      ELSEIF(IVGVAR.EQ.1) THEN
         READ(UR_INP,'(A)') FILEHIST
         OPEN(UR_HIS,file=FILEHIST,status='old')
         CALL VAR_VEL_GRAD_evp(0,dtime,DBAR,dsim,
     $        iDBARv,iSBARv,ILBAR,LIJBARc,ictrl,strain_control,nsteps)
         NSTEPSX=NSTEPS !! as istep starts with zero
         write(*,*)'nstepsx:',nstepsx
      elseif(ivgvar.eq.9) then
         read(UR_INP,*) time_srx, dtime, tempdot
         iph=1
         OPEN(UW_SRX_TEX-1+iph,file='REXDD_SRX_TEX_PH'
     $        //char(48+iph)//'.OUT',status='unknown')
         include "call_rexdd_srx.f"
         close(UW_SRX_TEX-1+iph)
         go to 3000
      else if(ivgvar.eq.-4) then
         read(UR_INP,*) dum
         dum1=0.
         dum2=0.
         call euler(2, dum, dum1, dum2, aux33) ! transformation matrix
         write(*,'(3f8.3)') ((aux33(i,j),j=1,3),i=1,3)
         call texture_rotation_fe(aux33,interaction,ngr,NPH,
     $        nphmx,ngrmx,ishape,ifrag,crit_shp,iflat,axisph,axisgr,
     $        ag,eulerph,fijph,fijgr)
         go to 3000
      elseIF(IVGVAR.EQ.-6) THEN !! unloading condition
         READ(UR_INP,*) nsteps, total_unload_time
         iopt=0
         CALL UNLOAD_CONDITION(iopt,total_unload_time,dtime,
     $        ILBAR,nsteps,strain_control,LIJBARc,iSBARv,
     $        iDBARv,dsim,scauchy,sbar,dsbar,DBAR)
         nstepsx=nsteps
         ictrl=7
      elseif(ivgvar.eq.-7) then !! stress hold
         READ(UR_INP,*) nsteps, total_unload_time
         iopt=0
         CALL stresshold_CONDITION(iopt,total_unload_time,dtime,
     $        ILBAR,nsteps,strain_control,LIJBARc,iSBARv,
     $        iDBARv,dsim,scauchy,sbar,dsbar,DBAR)
         nstepsx=nsteps
         ictrl=7
      ELSE IF(IVGVAR.EQ.30) THEN
         write(*,*) '-----------------------------------'
         write(*,*) ' Forming Limit Diagram Path Gen #1'
         write(*,*) '-----------------------------------'
         read(UR_INP,*) rho_fld, cvm_fld, eijincr, evm_limit
         call fld_a_bc(
c             intent(in)
     $        rho_fld,cvm_fld,eijincr,
c             intent(out)
     $        LIJBARc,dsim,DBAR,scauchy,sbar,dsbar,dsimctrl,dtime,
     $        strain_control,ictrl,ILBAR,iDBARv,iSBARv)
c        open up 'region_a.bin' file, ready to store the simulated responses of region A.
         call fld_a_write(0,istep,URW_FLDA,sbar,dsbar,LIJBARc,DBAR,
     $        epstot,
     $        fijph(:,:,0),dtime)
      ELSE IF(IVGVAR.EQ.31) THEN
         write(*,*) '---------------------------------------------'
         write(*,*) 'Forming Limit Probe (Region B) Mixed Boundary'
         write(*,*) '---------------------------------------------'
         read(UR_INP,'(a)') fn_region_a
         read(UR_INP,*) ftilde0, psi0_fld
         epstot_a_lab(:,:)=0.
         psi0_fld=psi0_fld*pi/180.
         ftilde=ftilde0
         psi_fld=psi0_fld
         open(URW_FLDA,file=trim(fn_region_a),form='unformatted',
     $        access='sequential',status='old',action='read')
         call fld_b_bc(
c             intent(in)
c                **  sbar_old is referred in groove axes.
     $        0,URW_FLDA,psi_fld,ftilde,sbar_old,
c             intent(inout)
     $        t0_fld,
c             intent(out)
     $        LIJBARc,LIJBARc_a_lab,dsim,DBAR,scauchy,
     $        sbar,dsbar,sbar_a_lab,dtime,
     $        strain_control,ictrl,ILBAR,iDBARv,iSBARv)

         aux33(:,:)=0.
         aux6(:)=0.
         call fld_b_write(0,UW_FLDB,aux33,aux33,aux6,aux6,aux33,aux33,
     $        psi_fld,ftilde)
c$$$      ELSE IF(IVGVAR.EQ.2) THEN
c$$$         call stop_debug(0,'PCYS is not available')
c$$$      ELSE IF(IVGVAR.EQ.3) THEN
c$$$         call stop_debug(0,'lankford is not available')
c$$$      ELSE IF(IVGVAR.EQ.4) THEN
c$$$         write(*,*)'igvgvar.eq.4 is not working in dEVPSC'
c$$$         write(*,*)'considering using ivgvar.eq.-4'
c$$$         call xit()
c$$$      else if (ivgvar.eq.5) then !! write
c$$$c***  Statev feature in evpsc_main is suppressed due to too large
c$$$c***  arrays are demanded when running DRX simulations with grain clusters.
c$$$c
c$$$c$$$         iopt=1                 !store EVPSC arrays to file
c$$$c$$$         include 'call_statev_io.f'
c$$$c$$$         read(UR_INP,'(a)') fn_snap
c$$$c$$$         open(UR_STV,file=fn_snap,status='unknown')
c$$$c$$$         iopt=3                 !store EVPSC arrays to file
c$$$c$$$         include 'call_statev_io.f'
c$$$c$$$         close(UR_STV)
c$$$c$$$         go to 3000
c$$$      else if (ivgvar.eq.-5) then !! read
c$$$c***  Statev feature in evpsc_main is suppressed due to too large
c$$$c***  arrays are demanded when running DRX simulations with grain clusters.
c$$$c$$$         read(UR_INP,'(a)') fn_snap
c$$$c$$$         write(*,*)'fn_snap:',fn_snap
c$$$c$$$         open(UR_STV,file=fn_snap,status='old')
c$$$c$$$         iopt=4                 !read statev from a file
c$$$c$$$         include 'call_statev_io.f'
c$$$c$$$         iopt=2                 !statev -> EVPSC arrays
c$$$c$$$         include 'call_statev_io.f'
c$$$c$$$         sbar(:)=sbar_old(:)
c$$$c$$$         write(*,*)'sbar:',sbar
c$$$c$$$         write(*,*)'dsbar:',dsbar
c$$$c$$$         write(*,*)'sbar_old:',sbar_old
c$$$c$$$         ipx_loaded_ever=.true.
c$$$c$$$         close(UR_STV)
c$$$c$$$         go to 3000
      ENDIF






C *******************************************************************
C     DO LOOP OVER DEFORMATION STEPS
C *******************************************************************
      IF((IVGVAR.EQ.0.or.IVGVAR.eq.1).AND. NSTEPS.LT.NWRITE) THEN
         WRITE(*,'(/,'' Error *** NWRITE='',I3,''  > NSTEPS='',I3,
     $'' --> WILL NOT WRITE TEXTURE !!!'',/)') NWRITE,NSTEPS
         call xit()
      ENDIF

      IF(INTERACTION.EQ.0) ITERLBL=' *** FC     '
      IF(INTERACTION.EQ.1) ITERLBL=' *** AFFINE '
      IF(INTERACTION.EQ.2) ITERLBL=' *** SECANT '
      IF(INTERACTION.EQ.3) ITERLBL=' *** N-EFF  '
      IF(INTERACTION.EQ.4) ITERLBL=' *** TANGENT'

c     loop 2500: incremental steps and updates in material states
      istep=1
      do 2500 while (
     $         (ivgvar.eq.0 .or.ivgvar.eq.1 ).and.(istep.le.nstepsx)
     $     .or.(ivgvar.eq.30.or.ivgvar.eq.31).and.(.not.fld_limit)
     $     )

         WRITE(UW_ERR,'(/,''*******   STEP'',I4,5X,A)') ISTEP,ITERLBL
         WRITE(UW_LOG,'(/,''*******   STEP'',I4,5X,A)') ISTEP,ITERLBL
         WRITE(  *,'(/,''*******   STEP'',I8.6,''/'',i8.6,A)')
     $        ISTEP,nstepsx,ITERLBL


c ** Calculate self-consistent elasto-visco-plastic modulus tensor
c    for the polycrystal and single crystal elastic stiffness tensor
c    for individual grains
c    maybe it should be relocated to the end of incremental step.
         if      (
     $        IVGVAR.eq. 0
     $    .or.ivgvar.eq. 1
     $    .or.ivgvar.eq.-6
     $    .or.ivgvar.eq.-7
     $    .or.ivgvar.eq.30
     $    .or.ivgvar.eq.31
     $        ) then
            iopt=1
            include"call_elsc_evp_fe.f"
            call write_elsc(uw_elsc,1,selav,celav)
         endif

C *** IMPOSE VELOCITY GRADIENT (and/or stress) AT EACH STEP:
C       IF IMPOSING MONOTONIC thermo/mechanical LOADING (IVGVAR= -1 -> changes to 0)
C       IF IMPOSING MONOTONIC LOADING                   (IVGVAR= 0)
C       IF IMPOSING A HISTORY OF DEFORMATION            (IVGVAR= 1)
C       IF CALCULATING A 5D POLYCRYSTAL YIELD SURFACE   (IVGVAR=-2)
c       IF dsbar is imposed (mostly for Zr-ZrH case)    (IVGVAR=-10)
c       IF unloading is required                        (IVGVAR=-6)
c       IF stress hold is required                      (IVGVAR=-7)
c       IF running Marciniak-Kuczynski FLD calculations (IVGVAR=30 or ivgvar=31)

         IF(IVGVAR.EQ.1) then
            dsbar(:)=0.
            CALL VAR_VEL_GRAD_evp(1,dtime,DBAR,dsim,iDBARv,iSBARv,
     $           ILBAR,LIJBARc,idum,strain_control,idum)
         elseIF(IVGVAR.EQ.-6) then
            CALL unload_condition(1,dum1,dtime,iaux33,-1,idum,aux33,
     $           iaux6,iaux33,dsim,scauchy,sbar,dsbar,DBAR)
         elseIF(IVGVAR.EQ.-7) then
            CALL stresshold_condition(1,dum1,dtime,iaux33,-1,idum,aux33
     $           ,iaux6,iaux33,dsim,scauchy,sbar,dsbar,DBAR)
         elseif(ivgvar.eq.31) then

c          -------------------------------------------------------------------
c           Boundary condition for FLD-B simulations
            call fld_b_psi2transf(psi_fld,R_grvlab,R_labgrv)
c           Rotate such that the new reference axes align with
c     groove axes of region B.
            call texture_rotation_fe(R_labgrv,interaction,ngr,
     $           NPH,nphmx,ngrmx,ishape,ifrag,crit_shp,iflat,
     $           axisph,axisgr,
c     intent(out)
     $           ag,eulerph,fijph,fijgr)
c          eulerph, fijph, ag are all referenced to groove axes.

c          maybe this subroutine should be relocated elsewhere?
            call fld_b_bc(
c                intent(in)
     $           1,URW_FLDA,psi_fld,ftilde,
c                sbar_old (of region B) is referred in groove
     $           sbar_old_b_lab,
c                intent(inout)
     $           t0_fld,
c                intent(out)
c                 ** those quantities are w.r.t. grv axes
     $           LIJBARc,LIJBARc_a_lab,dsim,DBAR,scauchy,
     $           sbar,dsbar,sbar_a_lab,
     $
     $           dtime,
     $           strain_control,ictrl,ILBAR,iDBARv,iSBARv
     $           )
        endif

C *******************************************************************
c     Refect in grain orientation (AG) to update schmid tensor referred
c     in the sample axes (i.e., sch) from that in crystal axes (schca)
        do iph=1,NPH
           igb=ngr(iph-1)+1     ! grain bottom
           igt=ngr(iph)         ! grain top
           CALL UPDATE_SCHMID_fe(
c               in
     $          iph,igb,igt,nphmx,ngrmx,nsysmx,nsyst,schca,nschca,ag,
c               out
     $          sch,nsch)
        enddo
C *******************************************************************
C     FOR VOCE NORMALIZES ALL ARRAYS WITH UNITS OF STRESS, STRAIN RATE &
C     COMPLIANCE TO WORK WITH n POWERS OF ORDER ONE AND IMPROVE CONVERGENCE.
C     RENORMALIZES AT THE END OF EACH STEP.
C     SCALES GAMD0 IF IRATESENS=0 TO ELIMINATE RATE SENSITIVITY ASSOCIATED .
C     WITH 'n'.

        IF((IHARDLAW.EQ.0.or.IHARDLAW.EQ.20.or.ihardlaw.eq.3)
     $       .and.IRATESENS.EQ.0) THEN
           CALL SCALE_GAMD0S_evp_fe(nsysmx,nphmx,
     $          ngrmx,NPH,ngr,nsyst,DBAR,tg,gamd0s)
        ENDIF
        if (isx_ver.eq.1) then
           do iph=1,NPH
              call thermo_emod_cte(
c                  intent(in)
     $             1,UR_ECTE,iph,nsysmx,nmodmx,nphmx,nmodes,mxndeg,nsm,
     $             temp,schca,
c                  intent(inout) !! should be determined inside <evpsc_input_fe>
     $             ndegs_fit,cijv_fit,cte_fit,
c                  subjected to updates
     $             celcc,selcc,xmu_mode,athcc)
           enddo
        endif

C *******************************************************************
C     FOR DISL DENS LAW SCALE GAMD0 TO ELIMINATE n-ASSOCIATED RATE SENSITIVITY.
C     CALCULATE INITIAL CRRS's ASSOCIATED WITH PROCESS TEMP FOR FIRST STEP.
cnt --> will not work for IHARDLAW=20 and strain path changes !!
        if(ihardlaw.ge.20) then ! .and. istep.eq.1
c$$$           if (ihardlaw.ne.23) then
              call scale_gamd0s_evp_fe(nsysmx,nphmx,ngrmx,NPH,ngr
     $             ,nsyst,DBAR,tg,gamd0s)
c$$$           endif
           if(ihardlaw.eq.20.and.istep.eq.1.and.iproc.eq.1) then
              if (ielast.eq.0) then
                 call update_crss_dd_fe_new(
c     intent(in)
     $                1,0,UR_SX,UW_DIS,UW_TAU,nphmx,nsysmx,nslmod,
     $                nslsys,nmodes,nsyst,nsm,ntwmod,temp,ngr,wgt,ngrmx,
     $                DBAR,dtime,nmodmx,gamdot,ntwmmx,xmu_mode,
c     intent(inout)
c          read from (ioption.eq.0) and not evolving
     $                grsze,chifor,chisub,RHO_initW,TEMPref,burg,actener
     $                ,K1gener,DRAG,EDOT0,rho_inits,Atau0,Atau1,Atau2,
     $                Fdeb0,Fdeb1,Fdeb2,HPcoef,T0nucl,T1nucl,T2nucl,
     $                T0prop,T1prop,T2prop,ctwsl,cts1,cts0,cts2,K3recov,
     $                K4recov,D0vac,Efvac,Emvac,
c     subject to 'evolution'
     $                rho_deb,rho_s_for,rho_s_deb,
c     constants
     $                edotref,
c     intent(out)
     $                itwinlaw,hard,twthres,isectw,nrs,crss,rhosat,
     $                TAU_PROP,TAU_NUCL,
c     intent(out) for stats of DD law
     $                rho_for_ph_mod,rho_deb_ph,rho_tot_ph,
     $                rho_deb_ph_mod,rho_recov,taue_wavg,rhonucl,pnucl
     $                )
              endif
           endif
        endif


C *** UPDATES ORIENTATION OF TRANSFORMATION TENSOR DG_TRANS(I,KKK)
c          ! update dg_trans in accordance with the new orientation
c           note that athcc should be in b-basis
c           (Given as constant in Voigt notation and conversion is done inside data_crystal_fe.f)
        do iph=1,NPH
c          ** thermal eigen-strain in crystal axes
           cux6(:)=tempdot*athcc(:,iph)
           DO KKK=NGR(IPH-1)+1,NGR(IPH)
              aux33(:,:)=ag(:,:,kkk)
              call mat6rot(cux6,aux33,aux6)
              dg_trans(:,kkk)=aux6(:)
              bux6(:)=bux6(:)+aux6(:)*wgt(kkk)
           enddo
        enddo
c        CALL UPDATE_TRANS_RATE (1)      ! cnt !     Commented-out by YJ... 2017-08-10

C *******************************************************************
C     IF IMPOSING STRAIN, MAKES A ELASTIC GUESS FOR THE FIRST STEP WHEN
C     IRECOVER=0 OR STARTS FROM RECOVERED STRESS STATE 'SG' WHEN IRECOVER=1.
C     IF IMPOSING STRESS MAKES A SACHS GUESS FOR EVERY STEP.

        if(
     $       (.not.ipx_loaded_ever).and.
     $          ((istep.eq.1.and.irecover.eq.0)
     $       .or.(strain_control.eq.0.and.iproc.eq.1))
     $       ) then
           ipx_loaded_ever=.true.
           include"call_initial_guess.f"
        endif
c          initial guess
        if (iregime.eq.2) then
           dsg(:,:)=sgdot(:,:)*dtime
           sg(:,:)=sg_old(:,:)+dsg(:,:)
           DBAR0(:)=0.
        endif

c *************************************************************************
c     delta or sigma EVPSC calculation of stress and strain-rate tensors
c     for every grain in the modeled polycrsystal aggregate.
      if (print_loc) write(*,*)"** call_evpsc"
      include "call_evpsc.f"

C ***************************************************************************
C     (E)VPSC PROVIDES THE MACROSCOPIC 'TOTAL' RATE 'DSIM'.
C     CALCULATES MACROSCOPIC ROTATION RATE AND MACROSCOPIC VELOCITY GRADIENT
C     DEPENDING ON THE BOUNDARY CONDITIONS.
      if (print_loc) write(*,*)"** WBARc and LIJBARc is completed"
      WBARc(:,:)=0.
      DO I=1,3
         WBARc(I,I)=0.          ! TAKES CARE OF DIAGONAL ELEMENTS
         DO J=1,3
            IF(ILBAR(I,J).EQ.1 .AND. ILBAR(J,I).EQ.1) THEN
               WBARc(I,J)=(LIJBARc(I,J)-LIJBARc(J,I))/2.
            ENDIF
            IF(ILBAR(I,J).EQ.1 .AND. ILBAR(J,I).EQ.0) THEN
               WBARc(I,J)=LIJBARc(I,J)-DBARc(I,J)
               WBARc(J,I)=-WBARc(I,J)
            ENDIF
            LIJBARc(I,J)=DBARc(I,J)+WBARc(I,J)
         ENDDO
      ENDDO



C ***********************************************************************
C     GIVEN THE OVERALL STRAIN RATE AND STRESS THEN:
C     *  IF ICTRL=0 CALCULATES TIME INCREMENT NECESSARY TO ACHIEVE
C        THE IMPOSED VON MISES STRAIN INCREMENT.
C     *  IF ICTRL=1,6 CALCULATES TIME INCREMENT NECESSARY TO ACHIEVE
C        THE IMPOSED STRAIN COMPONENT INCREMENT.
C     *  IF ICTRL=7 AND STRAIN_CONTROL=1 A TIME INCREMENT IS IMPOSED.
C        --> CALCULATES THE VON MISES STRAIN INCR ASSOCIATED WITH THE
C            IMPOSED STRAIN-RATE.
C     *  IF STRAIN_CONTROL=0 IT CORRESPONDS TO A CREEP SIMULATION
C        --> CALCULATES THE STRAIN INCREMENT ASSOCIATED WITH THE IMPOSED
C            STRESS COMPONENT AND TIME INCREMENT
C ***********************************************************************

C     Equivalent stress and strain, and total work integration
      if (print_loc) write(*,*)"** Equivalent stress, strain "
      SVM=0.
      DVM=0.
      aux6(:)=dsbar(:)+sbar_old(:)
      DO I=1,5
         if (iregime.eq.1) then
            svm=svm+sbar(i)*sbar(i)
            DVM=DVM+sbar(i)*DBAR(i)
         else
            SVM=SVM+aux6(i)*aux6(i)
            DVM=DVM+aux6(i)*DBAR(i)
         endif
      ENDDO
      SVM=SQRT(SVM*3./2.)
      DVM=ABS(DVM)             ! TO AVOID <0 WHEN GROWTH DOMINATES CREEP
      IF(SVM.NE.0.) DVM=DVM/SVM       ! WORK CONJUGATE DEFINITION OF DVM
      IF(SVM.EQ.0.) DVM=SQRT(2./3.)*VNORM(DBAR(1:5),5)
      EPSVM=EPSVM+DVM*dtime

c     ----------------------------
      if(ictrl.eq.0) call stop_debug(0,'** Error: ictrl=0 @evpsc_main')

c     YJ is not sure if the below block is necessary.
      IF(ICTRL.EQ.0) THEN
         dtime  =EVMINCR/DVM
      ELSEIF(ICTRL.GE.1.AND.ICTRL.LE.6
     $        .AND. STRAIN_CONTROL.EQ.1) THEN
         dtime  =EIJINCR/ABS(DSIMCTRL)
      ELSE IF(ICTRL.EQ.7 .OR. STRAIN_CONTROL.EQ.0) THEN
         EVMINCR=DVM*dtime
      ENDIF


c---------------------------------------------------------------------------
      if (print_loc) write(*,*)"** Stress and strain integrations"
C *** Integration of time, temperature, stress, strain, and plastic work
!     for ivgvar.eq.31, both dsbar, sbar_old are grv axesc

c     time & temperature integration
      time=time+dtime
      temp=temp+tempdot*dtime

      if (iregime.eq.2) then
         sbar(:)=dsbar(:)+sbar_old(:)
         sbar_old(:)=sbar(:)
         call chg_basis(sbar,scauchy,aux66,aux3333,1,6)
         call chg_basis(sbar(1:5),sdeviat,aux66,aux3333,1,5)
         sg(:,:)    =sg_old(:,:)+dsg(:,:)
         sg_old(:,:)=sg(:,:)
      endif
c     Polycrystal strain integration
      EPSTOT(:,:)=EPSTOT(:,:)+DBARC(:,:)*dtime
      if (ivgvar.eq.30) then
         if (epsvm.ge.evm_limit) fld_limit=.true.
      endif
      if (ivgvar.eq.31) then ! FLD simulation
         do i=1,3
         do j=1,3
            aux33(i,j)=(LIJBARc_a_lab(i,j)+LIJBARc_a_lab(j,i))/2.
         enddo
         enddo
         epstot_a_lab(:,:)=epstot_a_lab(:,:)+aux33(:,:)*dtime
      endif

      call matvect_dot(selav,sbar,aux6,6) ! aux6_i = selav_ij sbar_j
      call chg_basis(aux6,epstot_el,aux66,aux3333,1,6)
c     estimate the elastic strain rate
      bux6(:)=dsbar(:)/dtime
      call matvect_dot(selav,bux6,aux6,6)
      call chg_basis(aux6,edot_el,aux66,aux3333,1,6)
c     the plastic work.
      do i=1,6
         plwork=plwork+0.5*(2*sbar(i)-dsbar(i))*DBAR(i)*dtime
         do j=1,6
            plwork=plwork+0.5*(2*sbar(i)-dsbar(i))*selav(i,j)*dsbar(j)
         enddo
      enddo




C *******************************************************************
C     STATISTICS & OUTPUT
C *******************************************************************
      IF(IVGVAR.LE.1
     $     .or.ivgvar.eq.30.or.ivgvar.eq.31 !! (FLD-a or FLD-b)
     $     .or.ivgvar.eq.-6.or.ivgvar.eq.-7) then

         call chg_basis(aux6,aux33,mbar,aux3333,3,6) ! compliance
         dum=1d0/aux3333(1,1,1,1) ! young's modulus
         call chg_basis(DBAR0,aux33,aux66,aux3333,1,6) ! DBAR011
         dum1=aux33(1,1)
c           hardwiring press and presswgt as zero - they are not used. (2021-July-27; YJ, MJ)
         press=0d0
         presswgt=0d0
c           ------------------------------------------------------------------
c           phase-specific stress, strain outputs
         if (NPH.gt.1) then
            do iph=1,NPH
               igb=ngr(iph-1) +1
               igt=ngr(iph)
               call phase_property_evp_fe(1,uw_phpr+iph-1,iph,
     $              igb,igt,ngrmx,nphmx,sg,tg,wgt,dtime,epstot_ph)
            enddo
         endif
      endif

c     ----------------------------------------------------------------------------------
      call stat_shear_activity_px(
c             in
     $     NPH,nphmx,ngrmx,ngr,nmodmx,nsysmx,wph,wgt,nmodes,
     $     nsyst,nsm,gamdot,
c             out
     $     activ,gavmod,gavacc,GAVMODG,GAVPH)

c     Report mean and standard deviation of granular stress/strain.
      idum=0
      iopt=1

      call grain_stat(iopt,istep,idum,uw_stat,uw_st01,uw_alls,ngr,nphmx,
     $     NPH,ngrmx,dtime,dgp,tg,dsg,sg_old,wgt,aux6,bux6,gtotgr)

C *************************************************************************
c     *** calculate internal strain on diffracting planes if idiff(iph).ne.0
      call DIFF_PLANES_PX(
c     intent(in)::
     $     1,idiff,UR_DIF,UW_DIF,0,NPH,ngrmx,nphmx,ndifmx,ngr,ag,wgt,
     $     filediff,ICRYSYMph,iregime,iDBARv,scauchy,sav,BC_EL,mgel,
     $     epstot,wph,
c     intent(inout)::
     $     ndiff,sg,angdetector,nfamily,diff_vs,diff_vc,NGRSETPH,
     $     WGTSETPH,IGRSET,
c     intent(out)::
     $     eelgr)

      if (itwinlaw.eq.3.or.irx.ge.1) then
c        Combined weight and internal lattice strain.
         call DIFF_PLANES_PX(
c     intent(in)::
     $        2,idiff,UR_DIF,UW_DIF,0,NPH,ngrmx,nphmx,ndifmx,ngr,ag,wgt,
     $        filediff,ICRYSYMph,iregime,iDBARv,scauchy,sav,BC_EL,mgel,
     $        epstot,wph,
c     intent(inout)::
     $        ndiff,sg,angdetector,nfamily,diff_vs,diff_vc,NGRSETPH,
     $        WGTSETPH,IGRSET,
c     intent(out)::
     $        eelgr)
      endif

c     skip hardening / updates.
      IF(IVGVAR.gt.2.and.(ivgvar.ne.30.and.ivgvar.ne.31)) GOTO 2500




C ***********************************************************************
C     UPDATES ORIENTATION, HARDENING AND SHAPE OF EVERY GRAIN MAKING
C     A LINEAR EXTRAPOLATION OF THE CALCULATED RATE TIMES 'dtime'.
C     (SKIPS UPDATE FOR PCYS OR LANKFORD RUN)
C ***********************************************************************

      aux330(:,:,:)=0d0
      aux33(:,:)=LIJBARc(:,:)
      CALL UPDATE_FIJ_fe(
c     in
     $     nphmx,ngrmx,0,0,0,dtime,aux33,ishape,aux330,
c     out
     $     fijph,fijgr)         ! UPDATES DEFORMATION TENSOR OF ELEMENT
      call update_shape_fe(ngrmx,nphmx,0,0,0,interaction,ifrag,crit_shp,
     $     ishape,fijph,iflat,axisph,axisgr,eulerph,fijgr,ierr_shp)

      DO IPH=1,NPH       !IPHBOT,IPHTOP
         IF(NTWMOD(IPH).NE.0) THEN
            igb=ngr(iph-1)+1
            igt=ngr(iph)
              IF(ITWINLAW.LE.2) then
                 CALL UPDATE_TWINNING_PTR_fe(
c                     in
     $                dtime,iph,igb,igt,nphmx,wph,igt-igb+1,ngrmx,
     $                nslsys,ntwmod,nslmod,twsh,nsm,nmodmx,ntwmmx,nsysmx
     $                ,gamdot,ntwsmx,wgt,itwinlaw,jran,twthres,dnca,dbca
     $                ,ntwsys,isectw,itwtype,scauchy,
c                     inout
     $                twfrsy,twfrph,ktwsmx,eftwfr,ag,dsg,ntwevents,
     $                pritw,sectw,sch_stat,var_stat,var_vf)
              elseif(itwinlaw.eq.3) then
                 call update_twinning_vft_fe(
c                     in
     $                dtime,iph,nphmx,ngrmx,ngr,nslmod,nslsys,nsysmx,
     $                ntwmod,ntwmmx,ntwsmx,twthres,nsm,gamdot,nmodmx,
     $                dnca,dbca,itwtype,nsyst,
c                     inout
     $                ktwsmx,twsh,twfrsy,twfrph,sg_old,sg,dsg,wph,
     $                ag,wgt,gtotgr,crss,
c                     out
     $                pritw)
c                 sgdot(:,igb:igt)=dsg(:,igb:igt)/time_incr_deriv
              endif
         ENDIF
      ENDDO

      if (iupdori.eq.1) then
         do iph=1,NPH
            igb=ngr(iph-1)+1
            igt=ngr(iph)
            CALL UPDATE_ORIENTATION_fe(
c     in
     $           dtime,interaction,nsysmx,ngrmx,nphmx,igb,igt,iph,
     $           ishape,wph,WBARc,asph,asgr,tg,DBAR,dnca,dbca,gamdot,
     $           wgt,nneigh,nneimx,neigh,wneigh,
c     inout
     $           ag,
c     out
     $           xlijph,xlijgr)
         enddo
      endif
      DO IPH=1,NPH
         igb=ngr(iph-1)+1
         igt=ngr(iph)
c          ** for non-zero phase, use xlijph not LIJBARc
         aux33(:,:)=xlijph(:,:,iph)
         CALL UPDATE_FIJ_fe(
c     in
     $        nphmx,ngrmx,iph,igb,igt,dtime,aux33,ishape,xlijgr,
c     out
     $        fijph,fijgr)
      enddo
      DO IPH=1,NPH
         igb=ngr(iph-1)+1
         igt=ngr(iph)
         IF(IUPDSHP.EQ.1) then
            CALL UPDATE_SHAPE_fe(
c                in
     $           ngrmx,nphmx,iph,ngr(iph-1)+1,ngr(iph),interaction,ifrag
     $           ,crit_shp,ishape,
!                inout
     $           fijph,iflat,axisph,axisgr,
c                out
     $           eulerph,fijgr,ierr_shp)
         endif
         if (ierr_shp.ne.0) then
            write(*,*)'** Error: ierr_shp is not 0'
            call xit()
         endif
      enddo

c     ----------------------------------------------------------------------------------
c     FLD print outs
      aux2(:)=0d0
      if (ivgvar.eq.30) then    ! region A
         call fld_a_write(1,istep,URW_FLDA,sbar,dsbar,LIJBARc,DBAR,
     $        epstot,fijph(:,:,0),dtime)
      elseif (ivgvar.eq.31) then
c     -------------------------------------------------------------------
c     Boundary condition for FLD-B simulations
         call fld_b_psi2transf(psi_fld,R_grvlab,R_labgrv)
         call matrot(epstot,R_labgrv,epstot_b_lab)
         call matrot(fijph(:,:,0),R_labgrv,aux33)
         call mat6rot(sbar_old,R_labgrv,sbar_old_b_lab)
         call mat6rot(sbar,R_labgrv,sbar_b_lab)
c *** update FLDb variables (psi, f)
         call fld_b_updates(aux33,t0_fld,epstot_a_lab,
     $        epstot_b_lab,ftilde0,
c             intent(out)
     $        ftilde,psi_fld)
c *** rotates back to the lab axes from grv axes.
         call fld_b_psi2transf(psi_fld,R_grvlab,R_labgrv)
c     groove axes of region B.
         call texture_rotation_fe(R_grvlab,interaction,ngr,
     $        NPH,nphmx,ngrmx,ishape,ifrag,crit_shp,iflat,
     $        axisph,axisgr,
c     intent(out)
     $        ag,eulerph,fijph,fijgr)
c     Examine if the forming limit criteria is met.
            write(*,*)'** forming limit is reached.'
            write(*,*)'LIJBARc and LIJBARc_a_lab'
            do i=1,3
               write(*,'(2(3f10.5,3x))')(LIJBARc(i,j),j=1,3),
     $              (LIJBARc_a_lab(i,j),j=1,3)
            enddo
         if (-LIJBARc(3,3).ge.LIJBARc_a_lab(3,3)*(-10.)) then
            fld_limit=.true.
         endif
         call fld_b_write(1,UW_FLDB,epstot_a_lab,epstot_b_lab,
     $        sbar_a_lab,sbar_b_lab,LIJBARc_a_lab,LIJBARc,
     $        psi_fld,ftilde)
      endif

      IF(IUPDHAR.EQ.1) THEN
C *** VOCE HARDENING PLUS PREDOMINANT TWIN REORIENTATION SCHEME
         do iph=1,NPH           !iphbot,iphtop
            igb=ngr(iph-1)+1
            igt=ngr(iph)
            IF(IHARDLAW.EQ.0.and.ielast.eq.0) then
               aux6(1)=1d0
               CALL UPDATE_CRSS_VOCE_fe(
     $              UR_SX,2,iph,nphmx,ngrmx,dtime,
c                     in
     $              nmodmx,ntwmmx,igb,igt,nmodes,nsysmx,nsm,nslsys,
     $              nsyst,gamdot,
c                     inout
     $              grsze,gtotgr,tau_voce,thet_voce,hpfac_voce,
     $              hard,crss,
c                     out
     $              twthres,nrs,itwinlaw,isectw)
c     for Glenn
            elseIF(IHARDLAW.EQ.3) then
               CALL UPDATE_CRSS_VOCE_thermo_fe(
c                     in
     $              UR_SX,2,iph,nphmx,ngrmx,dtime,nmodmx,ntwmmx,igb,
     $              igt,nmodes,nsysmx,nsm,nslsys,nsyst,gamdot,temp,
     $              tempdot,
c                     inout
     $              grsze,gtotgr,tau_voce,thet_voce,hpfac_voce,
     $              hard,crss,
c                     out
     $              twthres,nrs,itwinlaw,isectw)
            endif
         enddo                  ! end of phase

C *** MECHANICAL THRESHOLD STRESS HARDENING (NO TWINNING, ONLY ONE PHASE)
         IF(IHARDLAW.EQ.1) THEN
            write(*,*)'unexpected case of mts crss law'
            call xit()
c            TEMPX=TEMP
c            EDOT=SQRT(2./3.)*VNORM(DBAR(1:5),5)
c            CALL UPDATE_CRSS_MTS (EDOT,TEMPX,3)
c            TEMP=TEMPX
         ENDIF

C *** DISLOCATION DENSITY STRESS HARDENING (ONLY ONE PHASE)
         IF(IHARDLAW.EQ.20) then
            if (ielast.eq.0) then
               call update_crss_dd_fe_new(
c                   intent(in)
     $              2,0,UR_SX,UW_DIS,UW_TAU,nphmx,nsysmx,nslmod,nslsys,
     $              nmodes,nsyst,nsm,ntwmod,temp,ngr,wgt,ngrmx,DBAR,
     $              dtime,nmodmx,gamdot,NTWMMX,xmu_mode,
c                   intent(inout)
     $              grsze,chifor,chisub,RHO_initW,TEMPref,burg,actener,
     $              K1gener,DRAG,EDOT0,rho_inits,Atau0,Atau1,Atau2,
     $              Fdeb0,Fdeb1,Fdeb2,HPcoef,T0nucl,T1nucl,T2nucl,
     $              T0prop,T1prop,T2prop,ctwsl,cts1,cts0,cts2,K3recov,
     $              K4recov,D0vac,Efvac,Emvac,
c                   ** subject to 'evolution'
     $              rho_deb,rho_s_for,rho_s_deb,
c                   constants
     $              edotref,
c                   intent(out)
     $              itwinlaw,hard,twthres,isectw,nrs,crss,rhosat,
     $              TAU_PROP,TAU_NUCL,
c                   intent(out) for stats of DD law
     $              rho_for_ph_mod,rho_deb_ph,rho_tot_ph,rho_deb_ph_mod,
     $              rho_recov,taue_wavg,rhonucl,pnucl)
            endif
         endif

         do iph=1,NPH    !iphbot,iphtop
            igb=ngr(iph-1)+1
            igt=ngr(iph)
            tau0(:,:)=0d0 !!

            if (ihardlaw.eq.23.and.ielast.eq.0) then
               call update_crss_rgvb_fe(
c                   intent(in)
     $              UR_SX,2,iph,igb,igt,dtime,
     $              nsysmx,ngrmx,nphmx,nmodmx,ntwmmx,nsyst,nrs,nslmod,
     $              ntwmod,nslsys,ntwsys,nmodes,nsm,igt-igb+1,dnca,jran,
     $              wgt,
c                   inout
     $              TAUE,crss,gamdot,ishear_tag,ibs_law,xMu,Xb,xf_self,
     $              xf_lat,xD,xK,alpham,xlambda,tau0,xp,xbx,facs,pow,
     $              xRho_min,xRho_max,xRho_infl,rho_f,rho_r,taud,
     $              rho_prev,rho_deb_rgvb,tauini,rhos,rho_mod,fac,
     $              for_rhos,rev_rhos,cts,for_rhosl,rev_rhosl,AMSS,
     $              deltaD,xp_g,drho_tot,debs,rho_ps,tautwint,
c                   out
     $              isectw,twthres,itwinlaw,dum1,dum2,dum3,dum4)
            endif
         enddo                  ! phase loop
      ENDIF                     ! end of iupdhrd

C *******************************************************************
C *** Dynamic recrystallization
C *******************************************************************
      if (irx.eq.1)then
         iph=1
         iopt=0                 ! Dynamic recrystallization
         include "call_rexdd.f"
      endif

C *******************************************************************
C *** END OF FORWARD UPDATE OF STRAIN, SHAPE, ORIENTATION & CRSS
C *******************************************************************




C *******************************************************************
C     WRITE TEXTURE FILES FOR EACH PHASE
C *******************************************************************
c     --study to evaluate the weight-averaged stress for each mode
      call write_crav(1,uw_crav,NPH,ngr,nmodes,nphmx,ngrmx,
     $     nmodmx,nsm,nsysmx,wgt,crss)
      CALL WRITE_STRESS_STRAIN_evp(UW_STR,ISTEP+1,dtime,epsvm,svm,
     $     epstot,temp,LIJBARc,scauchy,press,presswgt,plwork,epstot_el,
     $     edot_el)

      if (print_loc) write(*,*)"** calc_contribution"

      call calc_contribution(
c          intent(in)
     $     uw_cntr,ife,ngrmx,nphmx,NPH,ngr,dgp,dge,dg_trans,wgt,
c          intent(out)
     $     contrib)

      if (print_loc) write(*,*)"** write_active_fe"

      do iph=1,NPH
         call write_activ_fe(1,uw_ACT+iph,uw_GAM+iph,iph,ihardlaw,
     $        nmodes(iph),ntwmod(iph),wph(iph),epsvm,epstot,
     $        twfrph(1:ntwmod(iph),iph),eftwfr(1:ntwmod(iph),iph),
     $        pritw(iph),sectw(iph),activ(iph),
     $        gavmod(1:nmodes(iph),iph),gavacc(1:nmodes(iph),iph),
     $        GAVPH(IPH))
      enddo

      IF(IVGVAR.LE.1.or.ivgvar.eq.31.or.ivgvar.eq.30) THEN
         ItextWRITE=0
         IF(NWRITE.EQ.0) THEN
            IF(ISTEP.EQ.NSTEPS) ItextWRITE=1
         ELSE IF(NWRITE.NE.0) THEN
            IF(MOD(ISTEP,NWRITE).EQ.0) ItextWRITE=1
         ENDIF
         if (ivgvar.eq.30) itextwrite=1

         IF(ItextWRITE.EQ.1)then
            if (print_loc) write(*,*)"** write_active_fe"
            do iph=1,NPH
               igb=ngr(iph-1)+1
               igt=ngr(iph)
               CALL WRITE_TEXTURE_fe(UW_TEX-1+iph,iph,ngrmx,nphmx,
     $              igb,igt,EPSVM,ag,wgt,axisph,eulerph)
            enddo
            if (itwinlaw.eq.3.and.NPH.eq.2) then
               CALL WRITE_TEXTURE_fe(UW_TEX+nph,1,ngrmx,nphmx,1,
     $              ngr(NPH),epsvm,ag,wgt,axisph,eulerph)
            endif
         endif
      ENDIF

      IF(ISAVE.EQ.ISTEP) THEN
         call stop_debug(0,'postmortem is not working for this version')
      ENDIF

c     This is to check if 'void' growth can be simulated by EVPSC.
c     adjust phase volume fraction by
      if (NPH.gt.1.and..false.) then
c$$$         aux33(:,:)=fijph(:,:,1)
c$$$         detm=det(aux33)
c$$$         aux33(:,:)=fijph(:,:,2)
c$$$         detv=det(aux33)
c$$$c        update in the void volume fraction
c$$$         dum=w0v/(1.-w0v)
c$$$         dum=detv*dum/(detm+detv*dum)

         dph(:,:)=0.
         do iph=1,NPH
            dph(:,iph)=0.
            do kkk=ngr(iph-1)+1,ngr(iph)
               do i=1,6
                  dph(i,iph)=dph(i,iph)+tg(i,kkk)*wgt(kkk)
               enddo
            enddo
            dph(:,iph)=dph(:,iph)/wph(iph)
         enddo

c         dph(6,1)*sqrt(3)       !! trace
         write(*,*)'trace of ph1:',dph(6,1)*sqrt(3.)
         write(*,*)'trace of ph2:',dph(6,2)*sqrt(3.)

         wmatl=wmatl+(1-wmatl)*dph(6,1)*sqrt(3.)*dtime
         fdot=(1-wvoid)*dph(6,2)*sqrt(3.)
         wvoid=wvoid+fdot*dtime

c         dum=wmatl+wvoid
c         wmatl=wmatl/dum
c         wvoid=wvoid/dum
         ! wph(2)=wvoid
         ! wph(1)=1.-wph(2)
c         wph(1)=wmatl

         do iph=1,NPH
            wgtsum=0.
            do kkk=ngr(iph-1)+1,ngr(iph)
               wgtsum=wgtsum+wgt(kkk)
            enddo
            do kkk=ngr(iph-1)+1,ngr(iph)
               wgt(kkk)=wgt(kkk)*wph(iph)/wgtsum
            enddo
         enddo

         write(uw_void,'(3(e11.3,x))')1-wvoid,wvoid,fdot
         write(*,*)'wph change!'
         write(*,'(3(e14.6,x))')1-wvoid,wvoid,fdot
      endif ! end of void growth experiments

c     END OF DO LOOP OVER DEFORMATION INCREMENTS
      write(*,'(2a10)')'EPSVM','SVM'
      write(*,'(2f10.3)') EPSVM,SVM
      istep=istep+1
C     WRITES GRAIN'S & PX STATES INTO BINARY FILE 'POSTMORT.OUT'
 2500 CONTINUE
c     END OF DO LOOP OVER PROCESSES
 3000 CONTINUE



C *******************************************************************
c     end of process and write down the texture.
c     (Particularly relevant when texture is rotated by ivgvar.eq.-4)
C *******************************************************************
      do iph=1,NPH
c         write(*,*)'writing texture at the end of process'
         igb=ngr(iph-1)+1
         igt=ngr(iph)
         CALL WRITE_TEXTURE_fe(UW_TEX-1+iph,iph,ngrmx,nphmx,igb,igt,
     $        EPSVM,ag,wgt,axisph,eulerph)
      enddo
      close(3014)

c     close opened files...
      DO IPH=1,NPH
         IF (IDIFF(IPH).EQ.1) THEN
            CLOSE(UW_DIF+IPH)   !! 'DIF_STR_PHX.OUT'
            CLOSE(UW_DIF+IPH+5) !! 'DIF_WGT_PHX.OUT'
         ENDIF
         CLOSE(UW_ESGR+IPH)
         CLOSE(UW_ACT+IPH)
         CLOSE(UW_GAM+IPH)
      ENDDO

      IF (IVGVAR.EQ.30) CLOSE(URW_FLDA)
      IF (IVGVAR.EQ.31) CLOSE(UW_FLDB)
      CLOSE(UW_ELSC)
      CLOSE(UW_STR)
      CLOSE(UW_SCCN)
      CLOSE(UW_CRAV)
      CLOSE(UW_ST01)
      CLOSE(UW_STAT)
      CLOSE(UW_ALLS)
      if (iodf.eq.3) close(UR_EBSD)
C *******************************************************************

      call cpu_time(time_end)
      call w_val(0,'Total elapsed time (s):', time_end-time_begin)
      call date_and_time(date_,time_,zone_,itime_val)
      write(*,'(a,2x,2(i2.2,a1),i2.2,a1,i3.3)')
     $     date_,itime_val(5),':',itime_val(6),':',itime_val(7),
     $     '.',itime_val(8)

C *******************************************************************
c     Ending message to tell if the simulation was successful or not.
      write(*,*)'********************'
      write(*,*)'dEVPSC run completed'
      write(*,*)'********************'

      end program evpsc_main

C*************************************************************************
C     INCLUDE THE LIBRARY OF CANNED SUBROUTINES 'LIBRARY.SUB'
C*************************************************************************
      INCLUDE 'library7.sub'
