c----------------------------------------------------------------------c
c     This calling block of subroutine <initial_state_guess>  is shared
c     by Stand-alone and FE-interface versions.
c----------------------------------------------------------------------c
      CALL EVPSC_INPUT_fe(
     $     iodf,    ! iodf
     $     ife,     ! ife  0: standalone run, 1: Abaqus FE run
     $     nphmx,UR_INP,UR_SX,UR_TX,UR_SHP,UR_ECTE,UW_LOG,
     $     UW_DIS,UW_TAU,dtime,
     $     nsysmx,ngrmx,ntwmmx,nmodmx,ntwsmx,nneimx,mxndeg,
c          out
     $     isx_ver,ihardlaw,iratesens,errs,errd,errm,errso,itmaxext,
     $     itmaxint,itmaxso,IRSVAR,JXRSINI,JXRSFIN,JXRSTEP,ibcinv,
     $     IRECOVER,ISAVE,
     $     ICUBCOM,NWRITE,interaction,IUPDORI,IUPDSHP,IUPDHAR,NNEIGH,
     $     iflu,icompgrain,gamd0s,ngr,NPH,ishape,fijph,ifrag,
     $     iflat,crit_shp,axisph,axisgr,eulerph,fijgr,idiff,selcc,
     $     nmodes,nslmod,ntwmod,nsyst,nslsys,ntwsys,twsh,nsm,wgt,schca,
     $     nschca,c6,isense,itwtype,hard,isectw,crss,twthres,
     $     nrs,itwinlaw,ag,icrysymph,nunique,nrsmin,athcc,celcc,
     $     dnca,dbca,neigh,wneigh,wph,
c          voce related
     $     tau_voce,thet_voce,hpfac_voce,
c          rgvb related
     $     taue,ishear_tag,ibs_law,xmu,xb,xf_self,xf_lat,xD,xK,alpham,
     $     xlambda,xp,xbx,facs,pow,xRho_min,xRho_max,xRho_infl,rho_f,
     $     rho_r,taud,rho_prev,rho_deb_rgvb,tauini,rhos,rho_mod,fac,
     $     for_rhos,rev_rhos,cts,for_rhosl,rev_rhosl,AMSS,deltaD,xp_g,
     $     drho_tot,debs,rho_ps,tautwint,fstar,filediff
c          update_crss_dd_fe_new
     $     ,xmu_mode,grsze,chifor,chisub,rho_initw,TEMPref,burg,actener,
     $     K1gener,DRAG,EDOT0,rho_inits,Atau0,Atau1,Atau2,Fdeb0,Fdeb1,
     $     Fdeb2,HPcoef,T0nucl,T1nucl,T2nucl,T0prop,T1prop,T2prop,CTWSL,
     $     cts1,cts0,cts2,
     $     K3recov,K4recov,D0vac,Efvac,Emvac,
     $     rho_deb,rho_s_for,rho_s_deb,
c                --constants from update_crss_dd_fe_new
     $     edotref

c          polynomial fit of temperature-dependent el modulus and ctes.
     $     ,cijv_fit,cte_fit,ndegs_fit
c          coordinate for selecting ebsd
     $     ,coords,ebsd_rad
c     recrystallization
     $     ,irx
     $     ,rexgrowpar,rexnuclpar,radius
     $     )
