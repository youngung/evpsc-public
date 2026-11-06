c----------------------------------------------------------------------c
c     This calling block of subroutine <initial_state_guess>  is shared
c     by Stand-alone and FE-interface versions.
c----------------------------------------------------------------------c
c     iopt: option to determine the behavior of the subroutine
c           ( 0): test if nstatv is properly sized.
c           ( 1): Save dEVPSC arrays to STATEV            (devpsc->statev)
c           ( 2): Load dEVPSC arrays from STATEV          (statev->devpsc)
c           ( 3): Save entire state variables to a file   (statev->file  )
c           ( 4): Load entire state variables from a file (file  ->statev)

      call statev_io(iopt,nvar,statev,nstatv,UR_STV
c     FE/IO input related condition
c     $     ,is_initialized_statev
c     AVERAGE
     $        ,DSBAR,sbardot,EPSVM,SVM
c     DATACRYST
     $        ,DNCA,DBCA,SCHCA,NSCHCA,C6,HARD,NRS,ISENSE,ITWTYPE
c     DATAGRAIN
     $        ,AG,WGT,NGR
c     DATAMODES
     $        ,TWTHRES,TWSH,NMODES,NSLMOD,NTWMOD,NSYST,NTWSYS,NSLSYS
     $        ,NSM,ISECTW
c     DATAPHASE
     $        ,WPH,EULERPH,AXISPH,FIJPH,CRIT_SHP,NPH
     $        ,NELEM,icrysymph,ISHAPE,IFLAT,IFRAG
c     grainprops
     $        ,CRSS,SG,dsg,sgdot,GTOTGR,GAMD0S,tg,dge
c     RUNCOND
     $        ,ITMAXEXT,ITMAXINT,ITMAXSO,IRECOVER,ISAVE,NPROC,IVGVAR
     $        ,NWRITE,INTERACTION,ICUBCOM,ICS,NUNIQUE,NRSMIN,IHARDLAW
     $        ,IRATESENS,ICOMPGRAIN,IUPDORI,IUPDSHP,IUPDHAR,IBCINV
     $        ,IRSVAR,IFLU,ITWINLAW,jxrsini,jxrsfin,jxrstep,ERRS,ERRD,
     $        ERRM,ERRSO,WORK,plwork
c     STATTW
     $        ,EFTWFR,TWFRPH,TWFRSY,PRITW,SECTW,KTWSMX,NTWEVENTS
c     ELASTIC
     $        ,CELCC,SELCC,celav,sbar_old
c     DIFFRACT
     $        ,idiff,ndiff,angdetector,nfamily,diff_vc,diff_vs
c     EVPMOD
     $        ,dg0
c     EVP
     $        ,sg_old,dav,DBAR0,LBAR,mbar,mbardot,mgtg,DBAR0dot
c     Voce hardening
     $        ,tau_voce,thet_voce,hpfac_voce
c     RGVB hardening
     $        ,ishear_tag,ibs_law,xMu,xB,tau0,xf_self,xf_lat,xD,xK,
     $        alpham,xlambda,xp,xbx,facs,pow,xRho_min,xRho_max,xRho_infl
     $        ,rho_f,rho_r,taud,rho_prev,rho_deb_rgvb,tauini,rhos,
     $        rho_mod,fac,for_rhos,rev_rhos,for_rhosl,rev_rhosl,AMSS,
     $        deltaD,xp_g,drho_tot,debs,rho_ps,cts,tautwint,


c     dd hardening
     $     grsze,chifor,chisub,RHO_initW
     $     ,TEMPref,burg,actener,K1gener
     $     ,DRAG,EDOT0,rho_inits,Atau0
     $     ,Atau1,Atau2,Fdeb0,Fdeb1
     $     ,Fdeb2,HPcoef,T0nucl,T1nucl
     $     ,T2nucl,T0prop,T1prop,T2prop
     $     ,CTWSL,CTS1,CTS2
     $     ,CTS0,K3recov,K4recov,D0vac,Efvac,Emvac
     $     ,rho_deb,rho_s_for
     $     ,rho_s_deb
c          constants that are not subjected to evolution.
     $     ,edotref,xmu_mode

c     EVPSC & fstar, alplpx,
     $     ,fstar,alplpx
c     Accumulated rotation matrix
     $     ,rot_acum
c     accumulated time
     $     ,acctimefe,kount_umat_call
c     accumulated activity
     $     ,aa,ac
c     polynomial fit
     $     ,cijv_fit,cte_fit,ndegs_fit
     $     )
