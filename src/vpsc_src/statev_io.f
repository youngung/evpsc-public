      subroutine statev_io(iopt,nvar,statev,nstatv,iofn
c     FE/IO input related condition
c     $     ,is_initialized_statev
c     AVERAGE
     $     ,DSBAR,sbardot,EPSVM,SVM
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
     $        ,IRSVAR,IFLU,ITWINLAW,jxrsini,jxrsfin,jxrstep
     $        ,ERRS,ERRD,ERRM,ERRSO,WORK,plwork
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
     $        ,ishear_tag,ibs_law,xMu,xB,tau0,xf_self,xf_lat,xD
     $        ,xK,alpham,xlambda,xp,xbx,facs,pow,xRho_min,xRho_max
     $        ,xRho_infl,rho_f,rho_r,taud,rho_prev,rho_deb_rgvb,tauini
     $        ,rhos,rho_mod,fac,for_rhos,rev_rhos,for_rhosl,rev_rhosl
     $        ,AMSS,deltaD,xp_g,drho_tot,debs,rho_ps,cts,tautwint,

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


c     EVPSC & fstar
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
c     arguments list
c     iopt: option to determine the behavior of the subroutine
c           ( 0): test if nstatv is properly sized.
c           ( 1): Save state variables to STATEV
c           ( 2): Load state variables from STATEV
c           ( 3): Save entire state variables to a file
c           ( 4): Load entire state variables from a file

c     statev:        state variables of Abaqus platform
c     nstatv: No. of state variables of Abaqus platform
c      implicit none
      implicit  none
      include "incr_evpsc_fe.dim"
      integer, intent(in)::iopt,nstatv,iofn
      real*8, intent(inout)::statev(nstatv)
      integer, intent(inout)::nvar
c     local
      integer i,j,k,l,iv,iphel,kkk,isys,jsys,idum

      integer iverb
      parameter(iverb=-1)
c      parameter(iverb=0)        ! default
c      parameter(iverb=1)        ! activate this when debugging nstatv and iv
c      parameter(iverb=2)        ! print out for the case of iopt=2.
c     parameter(iverb=3)        ! print out for either iopt=2 or iopt=3
      logical irgvb,idd
      parameter(irgvb=.false.)
      parameter(idd=.false.)

c     to fool compiler flags / options for debugging
c      dtref=dtref
c      idtderiv=idtderiv
c      dtref=dtref
c      dum=dtref
c      idum=idtderiv
c      idum=iregime

      if (irgvb.and.idd) then
         call w_chr(0,'Cannot have both irgvb and idd .true.')
         call xit()
      endif

      iv = 0


      if (iverb.gt.1) write(*,*)'iofn:',iofn
      if (iverb.gt.1) write(*,*)'iopt',iopt

c     one should determine the size state variable per each polycrystal
c      write(*,*)'I am inside subr. <statev_io>'
      if (iopt.eq.0) then
         nvar=0
c----------------------------------------------------------------------
c        sbar_old(6)
         nvar=nvar+6
c$$   aa(nmodmx),ac(nphmx)
         nvar=nvar+nmodmx+nphmx
c----------------------------------------------------------------------
c     FE/IO input related condition
c         is_initialized_statev,
c         nvar=nvar+1
c----------------------------------------------------------------------
c     average
         nvar=nvar+1 ! EPSVM
         nvar=nvar+1 ! svm
         nvar=nvar+6 !dsbar
         nvar=nvar+6 !sbardot
c----------------------------------------------------------------------
c     DATACRYST
         do iphel=1,nphmx
            do isys=1,nsysmx
c              DNCA(3,NSYSMX,NPHMX)
               nvar=nvar+3
c              DBCA(3,NSYSMX,NPHMX),
               nvar=nvar+3
c              SCHCA(5,NSYSMX,NPHMX)
               nvar=nvar+5
c              NSCHCA(6,NSYSMX,NPHMX)
               nvar=nvar+6
c              **integers
c              NRS(NSYSMX,NPHMX),ISENSE(NSYSMX,NPHMX),ITWTYPE(NSYSMX,NPHMX)
               nvar=nvar+3
               do jsys=1,nsysmx
c                  HARD(NSYSMX,NSYSMX,NPHMX)
                  nvar=nvar+1
               enddo
            enddo
         enddo
c        C6
         nvar=nvar+1
c----------------------------------------------------------------------
c     DATAGRAIN
c         write(*,*)'nvar up before data grain',nvar
c         write(*,*)
         iv=0
         do kkk=1,ngrmx
c           AG(3,3,NGRMX)
            nvar=nvar+9
c           WGT(NGRMX)
            nvar=nvar+1
            iv=iv+10
         enddo
c           ** integer
c           NGR(0:NPHMX)
         nvar=nvar+(1+nphmx)
c----------------------------------------------------------------------
c     DATAMODES
c         write(*,*)'nvar up to data modes',nvar
c         write(*,*)
         do iphel=1,nphmx
            if (ntwmmx.gt.0) then
               do i=1,ntwmmx
c           TWTHRES(2,NTWMMX,NPHmx),TWSH(NTWMMX,NPHmx)
                  nvar=nvar+3
               enddo
            endif
c           **integer
c           NMODES(NPHmx),NSLMOD(NPHmx),NTWMOD(NPHmx),
            nvar=nvar+3
c           NSYST(NPHmx),NTWSYS(NPHmx),NSLSYS(NPHmx),
            nvar=nvar+3
            do i=1,nmodmx
c              NSM(NMODMX,NPHmx)
               nvar=nvar+1
            enddo
            do i=1,nsysmx
c              ISECTW(NSYSMX,NPHmx)
               nvar=nvar+1
            enddo
         enddo
c----------------------------------------------------------------------
c     DATAPHASE
c         write(*,*)'nvar up to data phase',nvar
c         write(*,*)
         do iphel=1,nphmx
c           WPH(NPHMX),
            nvar=nvar+1
c           CRIT_SHP(NPHMX)
            nvar=nvar+1
         enddo
         do iphel=0,nphmx
c           AXISPH(0:3,3,0:NPHMX)
            nvar=nvar+12
c           EULERPH(3,0:NPHMX)
            nvar=nvar+3
c           FIJPH(3,3,0:NPHMX)
            nvar=nvar+9
         enddo
c     ** integer
c        NPH,NELEM,icrysymph(nphmx),ISHAPE(0:NPHMX),
         nvar=nvar+2+nphmx+(1+nphmx)
c        IFLAT(0:NPHMX),IFRAG(NPHMX)
         nvar=nvar+(1+nphmx)+nphmx
c----------------------------------------------------------------------
c         write(*,*)' nvar up before grain props',nvar
c         write(*,*)
c     grainprops
         do kkk=1,ngrmx
            do i=1,nsysmx
c              CRSS(NSYSMX,NGRMX), gamd0s(NSYSMX,NGRMX)
               nvar=nvar+2
               iv=iv+2
            enddo
C           dsg(6,ngrmx),sg(6,ngrmx),sgdot(6,ngrmx),tg(6,ngrmx),dge(6,ngrmx),gtotgr(ngrmx)
            nvar=nvar+6*5+1
            iv=iv+6*5+1
         enddo
c----------------------------------------------------------------------
c         write(*,*)'************* nvar up before runcond',nvar
c     RUNCOND
c     integer
c        ITMAXEXT,ITMAXINT,ITMAXSO,IRECOVER,ISAVE,NPROC,IVGVAR,
c        NWRITE,INTERACTION,ICUBCOM,ICS,NUNIQUE,NRSMIN,IHARDLAW,
c        IRATESENS,ICOMPGRAIN,IUPDORI,IUPDSHP,IUPDHAR,IBCINV,
c        IRSVAR,IFLU,ITWINLAW,jxrsini,jxrsfin,jxrstep
         nvar=nvar+26
c        ERRS,ERRD,ERRM,ERRSO
         nvar=nvar+4
c         write(*,*)'nvar upto errso:',nvar
         nvar=nvar+ngrmx !! work
         nvar=nvar+1 !! plwork
c----------------------------------------------------------------------
c     STATTW
c        real*8 EFTWFR(NTWMMX,NPHMX),TWFRPH(NTWMMX,NPHMX),
         nvar=nvar+ntwmmx*nphmx*2
c        TWFRSY(NTWSMX,NGRMX),PRITW(NPHMX),SECTW(NPHMX)
         nvar=nvar+NTWSMX*NGRMX+NPHMX*2
         iv=iv+NTWSMX*NGRMX
c     integer
c        KTWSMX(NGRMX),NTWEVENTS(NGRMX)
         nvar=nvar+ngrmx*2
         iv=iv+ngrmx*2
c----------------------------------------------------------------------
c         write(*,*)'nvar before twin:',nvar
c     TWIN
c     integer
c        ICOMPGR(NGRMX)
c	 nvar=nvar+NGRMX
c----------------------------------------------------------------------
         if(iverb.eq.1)write(*,*)'nvar up before elastic',nvar
         if(iverb.eq.1)write(*,*)
c     ELASTIC
c        CELCC(6,6,NPHmx),SELCC(6,6,NPHmx)
         nvar=nvar+6*6*nphmx*2

c----------------------------------------------------------------------
c     DIFFRACT
c     integer
c        idiff(nphmx)
         nvar=nvar+nphmx
c        ndiff(nphmx)
         nvar=nvar+nphmx
c         angdetector(nphmx)
         nvar=nvar+nphmx
c     NFAMILY(NDIFMX,NPHMX),
         nvar=nvar+NDIFMX*NPHMX
c     DIFF_VC(3,24,NDIFMX,NPHMX)
         nvar=nvar+3*24*ndifmx*nphmx
c     DIFF_VS(3,NDIFMX,NPHMX)
         nvar=nvar+3*ndifmx*nphmx
c     ----------------------------------------------------------------------
c     EVPMOD
c        dg0(6,ngrmx) ! required in update_macro
         nvar=nvar+6*ngrmx
         iv=iv+6*ngrmx
c     ----------------------------------------------------------------------
c     EVP
c     integer
         if(iverb.eq.1)write(*,*)'nvar up evp',nvar
         if(iverb.eq.1)write(*,*)
c        sg_old(6,ngrmx),dav(6) !! dav should be in orientation?
         nvar=nvar+6*ngrmx+6
c        DBAR0(6),DBAR0dot(6)
         nvar=nvar+12
c        LBAR(6,6),mbar(6,6),mbardot(6,6),celav(6,6),mgtg(6,6,ngrmx)
c     mgtg is required in update_macro_property
         nvar=nvar+36*4+36*ngrmx
         iv=iv+36*ngrmx
         if(iverb.eq.1) write(*,*)'stat for grain props:',iv
c     ----------------------------------------------------------------------
c     Voce hardening
c     tau_voce(1:nsysmx,0:1,1:nphmx),
         nvar=nvar+nsysmx*2*nphmx
c     thet_voce(1:nsysmx,0:1,1:nphmx),
         nvar=nvar+nsysmx*2*nphmx
c     hpfac_voce(1:nsysmx,1:nphmx)
         nvar=nvar+nsysmx*nphmx
c     ----------------------------------------------------------------------

         if (irgvb) then
c     RGVB hardening
c     integer
c        ishear_tag(NSYSMX,NGRMX),ibs_law
            nvar=nvar+NSYSMX*NGRMX+1
c        xMu(nphmx,NMODMX),xB(nphmx,NMODMX),tau0(nphmx,NMODMX)
            nvar=nvar+nphmx*NMODMX*3
c$$$     $     xf_self(nphmx,NMODMX),xf_lat(nphmx,NMODMX),xD(nphmx,NMODMX),
            nvar=nvar+nphmx*NMODMX*3
c$$$     $     xK(nphmx,NMODMX),alpham(nphmx,NSYSMX,NSYSMX),
            nvar=nvar+nphmx*NMODMX+nphmx*NSYSMX*NSYSMX
c$$$     $     xlambda(nmodmx,NGRMX),xp(nsysmx,NGRMX),xbx(nsysmx),
            nvar=nvar+NMODMX*ngrmx+nsysmx*ngrmx+nsysmx
c$$$     $     facs(nphmx,nsysmx),pow(nphmx),xRho_min(nphmx,NMODMX),
            nvar=nvar+nphmx*nsysmx+nphmx+nphmx*NMODMX
c$$$     $     xRho_max(nphmx,NMODMX),xRho_infl(nphmx,NMODMX),
            nvar=nvar+nphmx*NMODMX*2
c$$$     $     rho_f(NSYSMX,NGRMX),rho_r(2*NSYSMX,NGRMX),taud(nsysmx),
            nvar=nvar+NSYSMX*NGRMX+2*NSYSMX*NGRMX+nsysmx
c$$$     $     rho_prev(NSYSMX,NGRMX),rho_deb_rgvb(nsysmx,NGRMX),
            nvar=nvar+NSYSMX*NGRMX*2
c$$$     $     tauini(nphmx,NSYSMX),rhos(nmodmx),rho_mod(nmodmx),
            nvar=nvar+nphmx*NSYSMX+nmodmx*2
c$$$     $     fac(nphmx,nmodmx),for_rhos(nmodmx),rev_rhos(nmodmx),
            nvar=nvar+nphmx*nmodmx+nmodmx*2
c$$$     $     for_rhosl(nmodmx),rev_rhosl(nmodmx),AMSS(NSYSMX,NSYSMX),
            nvar=nvar+nmodmx*2+NSYSMX*NSYSMX
c$$$     $     deltaD(NSYSMX,NGRMX),xp_g(NGRMX),drho_tot(NSYSMX,NGRMX),
            nvar=nvar+NSYSMX*ngrmx*2+ngrmx
c$$$     $     debs(NGRMX),rho_ps(NSYSMX,NGRMX),cts(NSYSMX,NSYSMX),
            nvar=nvar+ngrmx+NSYSMX*ngrmx+NSYSMX*NSYSMX
c$$$     $     tautwint(nsysmx)
            nvar=nvar+nsysmx
         endif ! end of rgvb
         if(iverb.eq.1)write(*,*)'nvar after rgvb',nvar
         if(iverb.eq.1)write(*,*)

         if (idd) then
c$$$      real*8, intent(inout):: grsze(nphmx),chifor,chisub,RHO_initW,
            nvar=nvar+nphmx+3
c$$$     $     TEMPref,burg(nmodmx),actener(nmodmx),K1gener(nmodmx),
            nvar=nvar+1+nmodmx*3
c$$$c$$$     $     DRAG(nmodmx),EDOT0(nmodmx),rho_inits(nmodmx),Atau0(nmodmx),
            nvar=nvar+nmodmx*4
c$$$c$$$     $     Atau1(nmodmx),Atau2(nmodmx),Fdeb0(nmodmx),Fdeb1(nmodmx),
            nvar=nvar+nmodmx*4
c$$$c$$$     $     Fdeb2(nmodmx),HPcoef(nmodmx),T0nucl(nmodmx),T1nucl(nmodmx),
            nvar=nvar+nmodmx*4
c$$$c$$$     $     T2nucl(nmodmx),T0prop(nmodmx),T1prop(nmodmx),T2prop(nmodmx),
            nvar=nvar+nmodmx*4
c$$$     $     CTWSL(NTWMMX,NMODMX),CTS1(NTWMMX,NMODMX),CTS2(NTWMMX,NMODMX),
            nvar=nvar+NTWMMX*nmodmx*3
c$$$     $     CTS0(NTWMMX,NMODMX),K3recov,K4recov,D0vac,Efvac,Emvac,
            nvar=nvar+NTWMMX*nmodmx
            nvar=nvar+5
c$$$     $     rho_deb(ngrmx),rho_s_for(NSYSMX,NGRMX),
            nvar=nvar+ngrmx+ngrmx*NSYSMX
c$$$     $     rho_s_deb(nsysmx,ngrmx)
            nvar=nvar+nsysmx*ngrmx
c$$$c          constants that are not subjected to evolution.
c$$$     $     ,edotref
            nvar=nvar+1
c              xmu_mode
            nvar=nvar+nmodmx*nphmx
         endif
         if(iverb.eq.1)write(*,*)'nvar after dd',nvar
         if(iverb.eq.1)write(*,*)

c$$$     rot_acum
         nvar=nvar+9
c$$$     fstar & alplpx
         nvar=nvar+2
c$$$     acctimefe
         nvar=nvar+1
c$$$     kount_umat_call
         nvar=nvar+1


c$$$  cijv and cte polynomial fits
c        ndegs_fit
         nvar=nvar+nphmx*2
c        cijv_fit
         nvar=nvar+nphmx*(mxndeg+1)*6*6
c        cte_fit
         nvar=nvar+nphmx*(mxndeg+1)*6

         if (nvar.gt.nstatv)then
            write(*,*)'** Error: <Nvar> exceeds <nstatv>'
            write(*,'(2(a7,i8.8,x))')'nvar:',nvar,'nstatv:',nstatv
c$$$            write(*,'(x,a33,x,i8.8)')
c$$$     $           '** Consider increasing Depvar from ',
c$$$     $           nstatv
c$$$            write(*,'(a37,i8.8)')
c$$$     $           ' to a value equal to or greater than ',
c$$$     $           nvar
            call xit
c            call stop_debug(0,'nvstatv is less than nvar')
         endif
         if(iverb.eq.1) then
            write(*,*)
            write(*,*)'---------------------'
            write(*,*)'iopt.eq.0 ****nvar:',nvar
            write(*,*)'---------------------'
         endif
      elseif (iopt.eq.1) then
         if (ihardlaw.eq.23.and..not.(irgvb)) then
            write(*,*)'**Error: Conflicts between ihardlaw and irgvb'
            write(*,*)'ihardlaw:',ihardlaw
            write(*,*)'irgvb:',irgvb
            call xit()
         endif
         if (ihardlaw.ne.23.and.irgvb) then
            write(*,*)'**Error: Conflicts between ihardlaw and irgvb'
            write(*,*)'ihardlaw:',ihardlaw
            write(*,*)'irgvb:',irgvb
            call xit()
         endif
         if (ihardlaw.eq.20.and..not.idd)then
            write(*,*)'**Error: Conflicts between ihardlaw and idd'
            write(*,*)'ihardlaw:',ihardlaw
            write(*,*)'idd:',idd
            call xit()
         endif
         if (ihardlaw.ne.20.and.idd) then
            write(*,*)'**Error: Conflicts between ihardlaw and idd'
            write(*,*)'ihardlaw:',ihardlaw
            write(*,*)'idd:',idd
            call xit()
         endif


         if (iverb.gt.1) write(*,*)'iopt 1 - 1'
c----------------------------------------------------------------------
c     Linear extrapolation
         iv=0
c$$$  sbar_old
         do i=1,6
            iv=iv+1
            statev(iv)=sbar_old(i)
         enddo
         do i=1,nmodmx
c----------------------------------------------------------------------
c     Relative activity
c$$      aa(nmodmx)
            iv=iv+1
            statev(iv)=aa(i)
         enddo
         do i=1,nphmx
c$$      ac(nphmx)
            iv=iv+1
            statev(iv)=ac(i)
         enddo

c----------------------------------------------------------------------
c        state variables other than the above-listed ones
c        are not stored to statev

c----------------------------------------------------------------------
c     FE/IO input related condition
c         iv=iv+1                   ! is_initialized
c         statev(iv) = real(is_initialized_statev,8)
c     average
         iv=iv+1
         statev(iv)=epsvm
         iv=iv+1
         statev(iv)=svm
         do i=1,6
            iv=iv+1
            statev(iv)=dsbar(i)
            iv=iv+1
            statev(iv)=sbardot(i)
         enddo
c----------------------------------------------------------------------
c     DATACRYST
         if (iverb.gt.1) write(*,*)'iopt 1 - 2'
         do i=1,nphmx
            do j=1,nsysmx
               do k=1,3
                  iv=iv+1
                  statev(iv)=DNCA(k,j,i)
                  iv=iv+1
                  statev(iv)=DBCA(k,j,i)
               enddo
               do k=1,5
                  iv=iv+1
                  statev(iv)=SCHCA(k,j,i)
               enddo
               do k=1,6
                  iv=iv+1
                  statev(iv)=NSCHCA(k,j,i)
               enddo
               do k=1,nsysmx
                  iv=iv+1
                  statev(iv)=HARD(k,j,i)
               enddo
               iv=iv+1
               statev(iv)=real(nrs(j,i),8)
               iv=iv+1
               statev(iv)=real(isense(j,i),8)
               iv=iv+1
               statev(iv)=real(ITWTYPE(j,i),8)
            enddo               ! nsysmx (j)
         enddo                  ! nphmx (i)
c        C6
         iv=iv+1
         statev(iv)=C6
c----------------------------------------------------------------------
c     DATAGRAIN
c         write(*,*)'iv up before datagrain',iv
c         write(*,*)
         if (iverb.gt.1) write(*,*)'iopt 1 - 3'
         do kkk=1,ngrmx
            do i=1,3
            do j=1,3
               iv=iv+1
               statev(iv)=AG(i,j,kkk)
            enddo
            enddo
            iv=iv+1
            statev(iv)=WGT(kkk)
         enddo                  ! ngrmx
         do iphel=0,nphmx
            iv=iv+1
            statev(iv)=real(NGR(iphel),8)
c            write(*,'(a12,f5.1)')'statev(iv):',statev(iv)
c            write(*,'(a12,i5.5)')'ngr(iphel):',ngr(iphel)
         enddo                  ! nphmx
c----------------------------------------------------------------------
c     DATAMODES
c         write(*,*)'iv up before data modes',iv
c         write(*,*)
         if (iverb.gt.1) write(*,*)'iopt 1 - 4'
         do iphel=1,nphmx
            if (ntwmmx.gt.0) then
               do j=1,ntwmmx
                  do i=1,2
                     iv=iv+1
                     statev(iv)=TWTHRES(i,j,iphel)
                  enddo
                  iv=iv+1
                  statev(iv)=TWSH(j,iphel)
               enddo            ! ntwmmx
            endif
            iv=iv+1
            statev(iv)=real(NMODES(iphel),8)
            iv=iv+1
            statev(iv)=real(NSLMOD(iphel),8)
            iv=iv+1
            statev(iv)=real(NTWMOD(iphel),8)
            iv=iv+1
            statev(iv)=real(NSYST(iphel),8)
            iv=iv+1
            statev(iv)=real(NTWSYS(iphel),8)
            iv=iv+1
            statev(iv)=real(NSLSYS(iphel),8)
            do j=1,nmodmx
               iv=iv+1
               statev(iv)=real(NSM(j,iphel),8)
            enddo
            do j=1,nsysmx
               iv=iv+1
               statev(iv)=real(ISECTW(j,iphel),8)
            enddo
         enddo                  ! nphmx
c----------------------------------------------------------------------
c     DATAPHASE
c         write(*,*)'iv up before  phase',iv
c         write(*,*)
         if (iverb.gt.1) write(*,*)'iopt 1 - 5'
         do iphel=1,nphmx
            iv=iv+1
            statev(iv)=WPH(iphel)
            iv=iv+1
            statev(iv)=CRIT_SHP(iphel)
         enddo                  ! nphmx
c         write(*,*)'##############',iv
         do iphel=0,nphmx
            do i=1,3
               iv=iv+1
               statev(iv)=EULERPH(i,iphel)
            do j=0,3
               iv=iv+1
               statev(iv)=AXISPH(j,i,iphel)
            enddo
            do j=1,3
               iv=iv+1
               statev(iv)=FIJPH(j,i,iphel)
            enddo
            enddo
         enddo
c$$$      integer
         iv=iv+1
         statev(iv)=real(NPH,8)
         iv=iv+1
         statev(iv)=real(NELEM,8)
         do iphel=1,nphmx
            iv=iv+1
            statev(iv)=real(icrysymph(iphel),8)
            iv=iv+1
            statev(iv)=real(ifrag(iphel),8)
         enddo
         do iphel=0,nphmx
            iv=iv+1
            statev(iv)=real(ISHAPE(iphel),8)
            iv=iv+1
            statev(iv)=real(iflat(iphel),8)
         enddo
c----------------------------------------------------------------------
c     grainprops
c         write(*,*)'iv up before grain props',iv
c         write(*,*)
         if (iverb.gt.1) write(*,*)'iopt 1 - 6'
         do kkk=1,ngrmx
            do i=1,nsysmx
               iv=iv+1
               statev(iv)=CRSS(i,kkk)
               iv=iv+1
               statev(iv)=GAMD0S(i,kkk)
            enddo
            do i=1,6
               iv=iv+1
               statev(iv)=dSG(i,kkk)
               iv=iv+1
               statev(iv)=SG(i,kkk)
               iv=iv+1
               statev(iv)=SGdot(i,kkk)
               iv=iv+1
               statev(iv)=tg(i,kkk)
               iv=iv+1
               statev(iv)=dge(i,kkk)
            enddo
            iv=iv+1
            statev(iv)=GTOTGR(kkk)
         enddo
c----------------------------------------------------------------------
         if (iverb.gt.1) write(*,*)'iopt 1 - 7'
c         write(*,*)'************ iv up before runcond',iv
c         write(*,*)
c     RUNCOND
c     ** Maybe most of the below are not required.
         iv=iv+1
         statev(iv)=real(ITMAXEXT,8)
         iv=iv+1
         statev(iv)=real(ITMAXINT,8)
         iv=iv+1
         statev(iv)=real(ITMAXSO,8)
         iv=iv+1
         statev(iv)=real(IRECOVER,8)
         iv=iv+1
         statev(iv)=real(ISAVE,8)
         iv=iv+1
         statev(iv)=real(NPROC,8)
         iv=iv+1
         statev(iv)=real(IVGVAR,8)
         iv=iv+1
         statev(iv)=real(NWRITE,8)
         iv=iv+1
         statev(iv)=real(INTERACTION,8)
         iv=iv+1
         statev(iv)=real(ICUBCOM,8)
         iv=iv+1
         statev(iv)=real(ICS,8)
         iv=iv+1
         statev(iv)=real(NUNIQUE,8)
         iv=iv+1
         statev(iv)=real(NRSMIN,8)
         iv=iv+1
         statev(iv)=real(IHARDLAW,8)
         iv=iv+1
         statev(iv)=real(IRATESENS,8)
         iv=iv+1
         statev(iv)=real(ICOMPGRAIN,8)
         iv=iv+1
         statev(iv)=real(IUPDORI,8)
         iv=iv+1
         statev(iv)=real(IUPDSHP,8)
         iv=iv+1
         statev(iv)=real(IUPDHAR,8)
         iv=iv+1
         statev(iv)=real(IBCINV,8)
         iv=iv+1
         statev(iv)=real(IRSVAR,8)
         iv=iv+1
         statev(iv)=real(IFLU,8)
         iv=iv+1
         statev(iv)=real(ITWINLAW,8)
         iv=iv+1
         statev(iv)=real(jxrsini,8)
         iv=iv+1
         statev(iv)=real(jxrsfin,8)
         iv=iv+1
         statev(iv)=real(jxrstep,8)
         iv=iv+1
         statev(iv)=ERRS
         iv=iv+1
         statev(iv)=ERRD
         iv=iv+1
         statev(iv)=ERRM
         iv=iv+1
         statev(iv)=ERRSO
c         write(*,*)'iv upto errso:',iv

         do kkk=1,ngrmx
            iv=iv+1
            statev(iv)=WORK(kkk)
         enddo
         iv=iv+1
         statev(iv)=plwork
c----------------------------------------------------------------------
c     STATTW
         if (iverb.gt.1) write(*,*)'iopt 1 - 8'
         do iphel=1,nphmx
            iv=iv+1
            statev(iv)=PRITW(IPHEL)
            iv=iv+1
            statev(iv)=SECTW(IPHEL)
            if (ntwmmx.gt.0) then
               do i=1,ntwmmx
                  iv=iv+1
                  statev(iv)=EFTWFR(i,IPHEL)
                  iv=iv+1
                  statev(iv)=TWFRPH(i,IPHEL)
               enddo
            endif
         enddo
         do kkk=1,ngrmx
            iv=iv+1
            statev(iv)=real(KTWSMX(KKK),8)
            iv=iv+1
            statev(iv)=real(NTWEVENTS(KKK),8)
            if (ntwsmx.gt.0)then
               do i=1,ntwsmx
                  iv=iv+1
                  statev(iv)=TWFRSY(i,KKK)
               enddo
            endif
         enddo
c----------------------------------------------------------------------
c         write(*,*)'iv before twin:',iv
c     TWIN
c	 do kkk=1,ngrmx
c            iv=iv+1
c	    statev(iv)=real(ICOMPGR(kkk),8)
c	 enddo
c----------------------------------------------------------------------
c     ELASTIC
         if(iverb.eq.1) write(*,*)'iv up before elastic',iv
         if(iverb.eq.1) write(*,*)
         if (iverb.gt.1) write(*,*)'iopt 1 - 9'
         do i=1,6
            do iphel=1,nphmx
               do j=1,6
                  iv=iv+1
                  statev(iv)=CELCC(i,j,iphel)
                  iv=iv+1
                  statev(iv)=SELCC(i,j,iphel)
               enddo
            enddo
c	    do kkk=1,ngrmx
c	       iv=iv+1
c	       statev(iv)=athcs(i,kkk)
c	       do j=1,6
c		  iv=iv+1
c		  statev(iv)=CELCS(i,j,kkk)
c		  iv=iv+1
c		  statev(iv)=mgel(i,j,kkk)
c		  iv=iv+1
c		  statev(iv)=AC_EL(i,j,kkk)
c		  iv=iv+1
c		  statev(iv)=bc_el(i,j,kkk)
c	       enddo
c	    enddo
         enddo
c----------------------------------------------------------------------
c     DIFFRACT
         do iphel=1,nphmx
            iv=iv+1
            statev(iv)=real(idiff(iphel),8)
            iv=iv+1
            statev(iv)=real(ndiff(iphel),8)
            iv=iv+1
            statev(iv)=angdetector(iphel)
            do i=1,ndifmx
               iv=iv+1
               statev(iv)=real(nfamily(i,iphel),8)
               do j=1,3
                  do k=1,24
                     iv=iv+1
                     statev(iv)=diff_vc(j,k,i,iphel)
                  enddo
                  iv=iv+1
                  statev(iv)=diff_vs(j,i,iphel)
               enddo
            enddo
         enddo
c----------------------------------------------------------------------
c     EVPMOD
         do kkk=1,ngrmx
            do i=1,6
               iv=iv+1
               statev(iv)=dg0(i,kkk)
            enddo
         enddo
c----------------------------------------------------------------------
c     EVP
         if(iverb.eq.1)write(*,*)'iv up evp',iv
         if(iverb.eq.1)write(*,*)

         do i=1,6
            iv=iv+1
            statev(iv)=dav(i)
            iv=iv+1
            statev(iv)=DBAR0(i)
            iv=iv+1
            statev(iv)=DBAR0dot(i)
            do j=1,6
               iv=iv+1
               statev(iv)=LBAR(i,j)
               iv=iv+1
               statev(iv)=mbar(i,j)
               iv=iv+1
               statev(iv)=mbardot(i,j)
               iv=iv+1
               statev(iv)=celAV(i,j)
            enddo
         enddo
         do kkk=1,ngrmx
            do i=1,6
               iv=iv+1
               statev(iv)=sg_old(i,kkk)
               do j=1,6
                  iv=iv+1
                  statev(iv)=mgtg(i,j,kkk)
               enddo
            enddo
         enddo
c----------------------------------------------------------------------
c     Voce hardening
         do iphel=1,nphmx
            do j=1,nsysmx
               do i=0,1
                  iv=iv+1
                  statev(iv)=tau_voce(j,i,iphel)
                  iv=iv+1
                  statev(iv)=thet_voce(j,i,iphel)
               enddo
               iv=iv+1
               statev(iv)=hpfac_voce(j,iphel)
            enddo
         enddo
c----------------------------------------------------------------------
c     RGVB hardening
         if (irgvb) then
            iv=iv+1
            statev(iv)=real(ibs_law,8)
            do kkk=1,ngrmx
               iv=iv+1
               statev(iv)=debs(kkk)
               iv=iv+1
               statev(iv)=xp_g(kkk)
               do i=1,nsysmx
                  iv=iv+1
                  statev(iv)=real(ishear_tag(i,kkk),8)
                  iv=iv+1
                  statev(iv)=rho_ps(i,kkk)
                  iv=iv+1
               statev(iv)=xp(i,kkk)
               iv=iv+1
               statev(iv)=deltad(i,kkk)
               iv=iv+1
               statev(iv)=drho_tot(i,kkk)
               iv=iv+1
               statev(iv)=rho_f(i,kkk)
               iv=iv+1
               statev(iv)=rho_r(2*i-1,kkk)
               iv=iv+1
               statev(iv)=rho_r(2*i,kkk)
               iv=iv+1
               statev(iv)=rho_prev(i,kkk)
               iv=iv+1
               statev(iv)=rho_deb_rgvb(i,kkk)
            enddo
            do i=1,nmodmx
               iv=iv+1
               statev(iv)=xlambda(i,kkk)
            enddo
         enddo                  ! NGRMX

         do iphel=1,nphmx
            iv=iv+1
            statev(iv)=pow(iphel)
            do i=1,nmodmx
               iv=iv+1
               statev(iv)=xMu(iphel,I)
               iv=iv+1
               statev(iv)=xB(iphel,I)
               iv=iv+1
               statev(iv)=tau0(iphel,I)
               iv=iv+1
               statev(iv)=xf_self(iphel,I)
               iv=iv+1
               statev(iv)=xf_lat(iphel,I)
               iv=iv+1
               statev(iv)=xD(iphel,I)
               iv=iv+1
               statev(iv)=xK(iphel,I)
               iv=iv+1
               statev(iv)=xRho_min(iphel,i)
               iv=iv+1
               statev(iv)=xRho_max(iphel,i)
               iv=iv+1
               statev(iv)=xRho_infl(iphel,i)
               iv=iv+1
               statev(iv)=fac(iphel,i)
            enddo
            do i=1,nsysmx
               iv=iv+1
               statev(iv)=facs(iphel,i)
               iv=iv+1
               statev(iv)=tauini(iphel,i)
c               if (iphel.eq.1.and.i.lt.3) then
c                  write(*,'(a12,f7.1)')'tauini(i,kkk):',tauini(iphel,i)
c                  write(*,'(a12,f7.1)')'statev(iv):',statev(iv)
c               endif
               do j=1,nsysmx
                  iv=iv+1
                  statev(iv)=alpham(iphel,i,j)
               enddo
            enddo
         enddo                  ! nphmx

         do i=1,nmodmx
            iv=iv+1
            statev(iv)=rhos(i)
            iv=iv+1
            statev(iv)=rho_mod(i)
            iv=iv+1
            statev(iv)=for_rhos(i)
            iv=iv+1
            statev(iv)=rev_rhos(i)
            iv=iv+1
            statev(iv)=for_rhosl(i)
            iv=iv+1
            statev(iv)=rev_rhosl(i)
         enddo
         do i=1,nsysmx
            do j=1,nsysmx
               iv=iv+1
               statev(iv)=AMSS(i,j)
               iv=iv+1
               statev(iv)=cts(i,j)
            enddo
            iv=iv+1
            statev(iv)=xbx(i)
            iv=iv+1
            statev(iv)=taud(i)
            iv=iv+1
            statev(iv)=tautwint(i)
         enddo
      endif                     ! end of rgvb
      if(iverb.eq.1)write(*,*)'iv after rgvb',iv
      if(iverb.eq.1)write(*,*)
c     ---------------------------
c     DD hardening
      if (idd) then
         do iphel=1,nphmx
            iv=iv+1
            statev(iv)=grsze(iphel)
         enddo
         iv=iv+1
         statev(iv)=chifor
         iv=iv+1
         statev(iv)=chisub
         iv=iv+1
         statev(iv)=RHO_initW
         iv=iv+1
         statev(iv)=TEMPref
         do i=1,nmodmx
            iv=iv+1
            statev(iv)=burg(i)
            iv=iv+1
            statev(iv)=actener(i)
            iv=iv+1
            statev(iv)=K1gener(i)
            iv=iv+1
            statev(iv)=DRAG(i)
            iv=iv+1
            statev(iv)=EDOT0(i)
            iv=iv+1
            statev(iv)=rho_inits(i)
            iv=iv+1
            statev(iv)=Atau0(i)
            iv=iv+1
            statev(iv)=Atau1(i)
            iv=iv+1
            statev(iv)=Atau2(i)
            iv=iv+1
            statev(iv)=Fdeb0(i)
            iv=iv+1
            statev(iv)=Fdeb1(i)
            iv=iv+1
            statev(iv)=Fdeb2(i)
            iv=iv+1
            statev(iv)=HPcoef(i)
            iv=iv+1
            statev(iv)=T0nucl(i)
            iv=iv+1
            statev(iv)=T1nucl(i)
            iv=iv+1
            statev(iv)=T2nucl(i)
            iv=iv+1
            statev(iv)=T0prop(i)
            iv=iv+1
            statev(iv)=T1prop(i)
            iv=iv+1
            statev(iv)=T2prop(i)
            if (ntwmmx.gt.0) then
               do j=1,ntwmmx
                  iv=iv+1
                  statev(iv)=CTWSL(J,I)
                  iv=iv+1
                  statev(iv)=CTS1(J,I)
                  iv=iv+1
                  statev(iv)=CTS2(J,I)
                  iv=iv+1
                  statev(iv)=CTS0(J,I)
               enddo            ! end of ntwmmx (j)
            endif
         enddo ! end of nmodmx (i)
         iv=iv+1
         statev(iv)=K3recov
         iv=iv+1
         statev(iv)=K4recov
         iv=iv+1
         statev(iv)=D0vac
         iv=iv+1
         statev(iv)=Efvac
         iv=iv+1
         statev(iv)=Emvac

         do i=1,ngrmx
            iv=iv+1
            statev(iv)=rho_deb(i)
            do j=1,nsysmx
               iv=iv+1
               statev(iv)=rho_s_for(j,i)
               iv=iv+1
               statev(iv)=rho_s_deb(j,i)
            enddo
         enddo
         iv=iv+1
         statev(iv)=edotref
         do i=1,nmodmx
            do j=1,nphmx
               iv=iv+1
               statev(iv)=xmu_mode(i,j)
            enddo
         enddo
      endif                     ! end of dd


c     rot_acum
         do i=1,3
            do j=1,3
               iv=iv+1
               statev(iv)=rot_acum(i,j)
            enddo
         enddo
c     fstar & alplpx
         iv=iv+1
         statev(iv)=fstar
         iv=iv+1
         statev(iv)=alplpx

         if (iverb.gt.1) write(*,*)'iopt 1 - 11'

c     acctimefe
         iv=iv+1
         statev(iv)=acctimefe

c     kount_umat_call
         if (iverb.gt.1) write(*,*)'iopt 1 - 12'
         iv=iv+1
         statev(iv)=real(kount_umat_call,8)

ccc   polyfits of cijv, cte
c     with mxndeg, ndegs_fit
c         b) ndegs_fit (nphmx,2)
         do i=1,nphmx
         do j=1,2
               iv=iv+1
               statev(iv)=real(ndegs_fit(i,j))
         enddo
c         c) cijv_fit, cte_fit
         do j=0,mxndeg
            do k=1,6
               do l=1,6
                  iv=iv+1
                  statev(iv)=cijv_fit(i,j,k,l)
               enddo
               iv=iv+1
               statev(iv)=cte_fit(i,j,k)
            enddo
         enddo
         enddo

         if (iverb.gt.1) write(*,*)'iopt 1 - 13'
         if (iv.gt.nstatv)then
            if (iverb.gt.1) write(*,*)'iopt 1 - 14'
            if (iverb.gt.1) write(*,*)'iopt 1 - 15'
            write(*,'(2(a15,i10.10))')'total iv:',iv,'nstatv:',
     $           nstatv
            call xit()
            if (iverb.gt.1) write(*,*)'iopt 1 - 16'
c     call stop_debug(0,'nvstatv is less than nvar')
         endif
         if (iverb.gt.1) write(*,*)'iopt 1 - 17'

         if(iverb.eq.1) then
         write(*,*)
         write(*,*)'---------------------'
         write(*,*)'iopt.eq.1 ****iv:',iv
         write(*,*)'---------------------'
         endif

      elseif (iopt.eq.2) then
c----------------------------------------------------------------------
         iv=0
c$$$  sbar_old
         do i=1,6
            iv=iv+1
            sbar_old(i)=statev(iv)
         enddo
c$$      aa(nmodmx)
         do i=1,nmodmx
            iv=iv+1
            aa(i)=statev(iv)
         enddo
         do i=1,nphmx
c$$      ac(nphmx)
            iv=iv+1
            ac(i)=statev(iv)
         enddo

c----------------------------------------------------------------------
c     FE/IO input related condition
c         is_initialized_statev,
c         iv=iv+1
c         is_initialized_statev= int(statev(iv))
c----------------------------------------------------------------------
c     average
         iv=iv+1
         epsvm=statev(iv)
         iv=iv+1
         svm=statev(iv)
         do i=1,6
            iv=iv+1
            dsbar(i)=statev(iv) ! initial guess is required for evpsc
            iv=iv+1
            sbardot(i)=statev(iv)
         enddo
c----------------------------------------------------------------------
c     DATACRYST
         do i=1,nphmx
            do j=1,nsysmx
               do k=1,3
                  iv=iv+1
                  DNCA(k,j,i)=statev(iv)
                  iv=iv+1
                  DBCA(k,j,i)=statev(iv)
               enddo
               do k=1,5
                  iv=iv+1
                  SCHCA(k,j,i)=statev(iv)
               enddo
               do k=1,6
                  iv=iv+1
                  NSCHCA(k,j,i)=statev(iv)
               enddo
               do k=1,nsysmx
                  iv=iv+1
                  HARD(k,j,i)=statev(iv)
               enddo
               iv=iv+1
               nrs(j,i)=int(statev(iv))
               iv=iv+1
               isense(j,i)=int(statev(iv))
               iv=iv+1
               ITWTYPE(j,i)=int(statev(iv))
            enddo               ! nsysmx (j)
         enddo                  ! nphmx (i)
c        C6
         iv=iv+1
         C6=statev(iv)
c----------------------------------------------------------------------
c     DATAGRAIN
         do kkk=1,ngrmx
            do i=1,3
            do j=1,3
               iv=iv+1
               AG(i,j,kkk)=statev(iv)
            enddo
            enddo
            iv=iv+1
            WGT(kkk)=statev(iv)
         enddo                  ! ngrmx
         do iphel=0,nphmx
            iv=iv+1
            NGR(iphel)=int(statev(iv))
         enddo                  ! nphmx
c----------------------------------------------------------------------
c     DATAMODES
         do iphel=1,nphmx
            if (ntwmmx.gt.0) then
               do j=1,ntwmmx
                  do i=1,2
                     iv=iv+1
                     TWTHRES(i,j,iphel)=statev(iv)
                  enddo
                  iv=iv+1
                  TWSH(j,iphel)=statev(iv)
               enddo            ! ntwmmx
            endif
            iv=iv+1
            NMODES(iphel)=int(statev(iv))
            iv=iv+1
            NSLMOD(iphel)=int(statev(iv))
            iv=iv+1
            NTWMOD(iphel)=int(statev(iv))
            iv=iv+1
            NSYST(iphel)=int(statev(iv))
            iv=iv+1
            NTWSYS(iphel)=int(statev(iv))
            iv=iv+1
            NSLSYS(iphel)=int(statev(iv))
            do j=1,nmodmx
               iv=iv+1
               NSM(j,iphel)=int(statev(iv))
            enddo
            do j=1,nsysmx
               iv=iv+1
               ISECTW(j,iphel)=int(statev(iv))
            enddo
         enddo                  ! nphmx
c----------------------------------------------------------------------
c     DATAPHASE
         do iphel=1,nphmx
            iv=iv+1
            WPH(iphel)=statev(iv)
            iv=iv+1
            CRIT_SHP(iphel)=statev(iv)
         enddo                  ! nphmx

         do iphel=0,nphmx
            do i=1,3
               iv=iv+1
               EULERPH(i,iphel)=statev(iv)
            do j=0,3
               iv=iv+1
               AXISPH(j,i,iphel)=statev(iv)
            enddo
            do j=1,3
               iv=iv+1
               FIJPH(j,i,iphel)=statev(iv)
            enddo
            enddo
         enddo
c$$$      integer
         iv=iv+1
         NPH=int(statev(iv))
         iv=iv+1
         NELEM=int(statev(iv))
         do iphel=1,nphmx
            iv=iv+1
            icrysymph(iphel)=int(statev(iv))
            iv=iv+1
            ifrag(iphel)=int(statev(iv))
         enddo
         do iphel=0,nphmx
            iv=iv+1
            ISHAPE(iphel)=int(statev(iv))
            iv=iv+1
            iflat(iphel)=int(statev(iv))
         enddo
c----------------------------------------------------------------------
c     grainprops
         do kkk=1,ngrmx
            do i=1,nsysmx
               iv=iv+1
               CRSS(i,kkk)=statev(iv)
               iv=iv+1
               GAMD0S(i,kkk)=statev(iv)
            enddo
            do i=1,6
               iv=iv+1
               dSG(i,kkk)=statev(iv)
               iv=iv+1
               SG(i,kkk)=statev(iv)
               iv=iv+1
               sgdot(i,kkk)=statev(iv)
               iv=iv+1
               tg(i,kkk)=statev(iv)
               iv=iv+1
               dge(i,kkk)=statev(iv)
            enddo
            iv=iv+1
            GTOTGR(kkk)=statev(iv)
         enddo
c----------------------------------------------------------------------
c     MISCEL
c	 iv=iv+1
c	 jran=int(statev(iv))

c----------------------------------------------------------------------
c     RUNCOND
c     ** Maybe most of the below are not required.
         iv=iv+1
         ITMAXEXT=int(statev(iv))
         iv=iv+1
         ITMAXINT=int(statev(iv))
         iv=iv+1
         ITMAXSO=int(statev(iv))
         iv=iv+1
         IRECOVER=int(statev(iv))
         iv=iv+1
         ISAVE=int(statev(iv))
         iv=iv+1
         NPROC=int(statev(iv))
         iv=iv+1
         IVGVAR=int(statev(iv))
         iv=iv+1
         NWRITE=int(statev(iv))
         iv=iv+1
         INTERACTION=int(statev(iv))
         iv=iv+1
         ICUBCOM=int(statev(iv))
         iv=iv+1
         ICS=int(statev(iv))
         iv=iv+1
         NUNIQUE=int(statev(iv))
         iv=iv+1
         NRSMIN=int(statev(iv))
         iv=iv+1
         IHARDLAW=int(statev(iv))
         iv=iv+1
         IRATESENS=int(statev(iv))
         iv=iv+1
         ICOMPGRAIN=int(statev(iv))
         iv=iv+1
         IUPDORI=int(statev(iv))
         iv=iv+1
         IUPDSHP=int(statev(iv))
         iv=iv+1
         IUPDHAR=int(statev(iv))
         iv=iv+1
         IBCINV=int(statev(iv))
         iv=iv+1
         IRSVAR=int(statev(iv))
         iv=iv+1
         IFLU=int(statev(iv))
         iv=iv+1
         ITWINLAW=int(statev(iv))
         iv=iv+1
         jxrsini=int(statev(iv))
         iv=iv+1
         jxrsfin=int(statev(iv))
         iv=iv+1
         jxrstep=int(statev(iv))
         iv=iv+1
         ERRS=statev(iv)
         iv=iv+1
         ERRD=statev(iv)
         iv=iv+1
         ERRM=statev(iv)
         iv=iv+1
         ERRSO=statev(iv)

         do kkk=1,ngrmx
            iv=iv+1
            WORK(kkk)=statev(iv)
         enddo
         iv=iv+1
         plwork=statev(iv)
c----------------------------------------------------------------------
c     STATTW
         do iphel=1,nphmx
            iv=iv+1
            PRITW(IPHEL)=statev(iv)
            iv=iv+1
            SECTW(IPHEL)=statev(iv)
            if (ntwmmx.gt.0) then
               do i=1,ntwmmx
                  iv=iv+1
                  EFTWFR(i,IPHEL)=statev(iv)
                  iv=iv+1
                  TWFRPH(i,IPHEL)=statev(iv)
               enddo
            endif
         enddo
         do kkk=1,ngrmx
            iv=iv+1
            KTWSMX(KKK)=int(statev(iv))
            iv=iv+1
            NTWEVENTS(KKK)=int(statev(iv))
            if (ntwsmx.gt.0) then
               do i=1,ntwsmx
                  iv=iv+1
                  TWFRSY(i,KKK)=statev(iv)
               enddo
            endif
         enddo
c----------------------------------------------------------------------
c     TWIN
c	 do kkk=1,ngrmx
c	    iv=iv+1
c	    ICOMPGR(kkk)=int(statev(iv))
c	 enddo
c----------------------------------------------------------------------
c     ELASTIC
         do i=1,6
            do iphel=1,nphmx
               do j=1,6
                  iv=iv+1
                  CELCC(i,j,iphel)=statev(iv)
                  iv=iv+1
                  SELCC(i,j,iphel)=statev(iv)
               enddo
            enddo
c	    do kkk=1,ngrmx
c	       iv=iv+1
c	       athcs(i,kkk)=statev(iv)
c	       do j=1,6
c                  iv=iv+1
c		  CELCS(i,j,kkk)=statev(iv)
c		  iv=iv+1
c		  mgel(i,j,kkk)=statev(iv)
c		  iv=iv+1
c		  AC_EL(i,j,kkk)=statev(iv)
c		  iv=iv+1
c		  bc_el(i,j,kkk)=statev(iv)
c	       enddo
c	    enddo
         enddo

c----------------------------------------------------------------------
c     DIFFRACT
         do iphel=1,nphmx
            iv=iv+1
            idiff(iphel)=int(statev(iv))
            iv=iv+1
            ndiff(iphel)=int(statev(iv))
            iv=iv+1
            angdetector(iphel)=statev(iv)
            do i=1,ndifmx
               iv=iv+1
               nfamily(i,iphel)=int(statev(iv))
               do j=1,3
                  do k=1,24
                     iv=iv+1
                     diff_vc(j,k,i,iphel)=statev(iv)
                  enddo
                  iv=iv+1
                  diff_vs(j,i,iphel)=statev(iv)
               enddo
            enddo
         enddo
c----------------------------------------------------------------------
c     EVPMOD
         do kkk=1,ngrmx
            do i=1,6
               iv=iv+1
               dg0(i,kkk)=statev(iv)
            enddo
         enddo
c----------------------------------------------------------------------
c     EVP

         do i=1,6
            iv=iv+1
            dav(i)=statev(iv)
            iv=iv+1
            DBAR0(i)=statev(iv)
            iv=iv+1
            DBAR0dot(i)=statev(iv)
            do j=1,6
               iv=iv+1
               LBAR(i,j)=statev(iv)
               iv=iv+1
               mbar(i,j)=statev(iv)
               iv=iv+1
               mbardot(i,j)=statev(iv)
               iv=iv+1
               CelAV(i,j)=statev(iv)
            enddo
         enddo
         do kkk=1,ngrmx
            do i=1,6
               iv=iv+1
               sg_old(i,kkk)=statev(iv)
               do j=1,6
                  iv=iv+1
                  mgtg(i,j,kkk)=statev(iv)
               enddo
            enddo
         enddo
c----------------------------------------------------------------------
c     Voce hardening
         do iphel=1,nphmx
            do j=1,nsysmx
               do i=0,1
                  iv=iv+1
                  tau_voce(j,i,iphel)=statev(iv)
                  iv=iv+1
                  thet_voce(j,i,iphel)=statev(iv)
               enddo
               iv=iv+1
               hpfac_voce(j,iphel)=statev(iv)
            enddo
         enddo
c----------------------------------------------------------------------
c     RGVB hardening
         if (irgvb)then
            iv=iv+1
            ibs_law=int(statev(iv))
         do kkk=1,ngrmx
            iv=iv+1
            debs(kkk)=statev(iv)
            iv=iv+1
            xp_g(kkk)=statev(iv)
            do i=1,nsysmx
               iv=iv+1
               ishear_tag(i,kkk)=int(statev(iv))
               iv=iv+1
               rho_ps(i,kkk)=statev(iv)
               iv=iv+1
               xp(i,kkk)=statev(iv)
               iv=iv+1
               deltad(i,kkk)=statev(iv)
               iv=iv+1
               drho_tot(i,kkk)=statev(iv)
               iv=iv+1
               rho_f(i,kkk)=statev(iv)
               iv=iv+1
               rho_r(2*i-1,kkk)=statev(iv)
               iv=iv+1
               rho_r(2*i,kkk)=statev(iv)
               iv=iv+1
               rho_prev(i,kkk)=statev(iv)
               iv=iv+1
               rho_deb_rgvb(i,kkk)=statev(iv)
            enddo
            do i=1,nmodmx
               iv=iv+1
               xlambda(i,kkk)=statev(iv)
            enddo
         enddo                  ! NGRMX

         do iphel=1,nphmx
            iv=iv+1
            pow(iphel)=statev(iv)
            do i=1,nmodmx
               iv=iv+1
               xMu(iphel,I)=statev(iv)
               iv=iv+1
               xB(iphel,I)=statev(iv)
               iv=iv+1
               tau0(iphel,I)=statev(iv)
               iv=iv+1
               xf_self(iphel,I)=statev(iv)
               iv=iv+1
               xf_lat(iphel,I)=statev(iv)
               iv=iv+1
               xD(iphel,I)=statev(iv)
               iv=iv+1
               xK(iphel,I)=statev(iv)
               iv=iv+1
               xRho_min(iphel,i)=statev(iv)
               iv=iv+1
               xRho_max(iphel,i)=statev(iv)
               iv=iv+1
               xRho_infl(iphel,i)=statev(iv)
               iv=iv+1
               fac(iphel,i)=statev(iv)
            enddo
            do i=1,nsysmx
               iv=iv+1
               facs(iphel,i)=statev(iv)
               iv=iv+1
               tauini(iphel,i)=statev(iv)
c               if (iphel.eq.1.and.i.lt.3) then
c                  write(*,'(a12,f7.1)')'tauini(i,kkk):',tauini(iphel,i)
c                  write(*,'(a12,f7.1)')'statev(iv):',statev(iv)
c               endif
               do j=1,nsysmx
                  iv=iv+1
                  alpham(iphel,i,j)=statev(iv)
               enddo
            enddo
         enddo                  ! nphmx
         do i=1,nmodmx
            iv=iv+1
            rhos(i)=statev(iv)
            iv=iv+1
            rho_mod(i)=statev(iv)
            iv=iv+1
            for_rhos(i)=statev(iv)
            iv=iv+1
            rev_rhos(i)=statev(iv)
            iv=iv+1
            for_rhosl(i)=statev(iv)
            iv=iv+1
            rev_rhosl(i)=statev(iv)
         enddo
         do i=1,nsysmx
            do j=1,nsysmx
               iv=iv+1
               AMSS(i,j)=statev(iv)
               iv=iv+1
               cts(i,j)=statev(iv)
            enddo
            iv=iv+1
            xbx(i)=statev(iv)
            iv=iv+1
            taud(i)=statev(iv)
            iv=iv+1
            tautwint(i)=statev(iv)
         enddo
      endif                    !rgvb

c     ---------------------------
c     DD hardening
      if (idd) then
         do iphel=1,nphmx
            iv=iv+1
            grsze(iphel)=statev(iv)
         enddo
         iv=iv+1
         chifor=statev(iv)
         iv=iv+1
         chisub=statev(iv)
         iv=iv+1
         RHO_initW=statev(iv)
         iv=iv+1
         TEMPref=statev(iv)
         do i=1,nmodmx
            iv=iv+1
            burg(i)=statev(iv)
            iv=iv+1
            actener(i)=statev(iv)
            iv=iv+1
            K1gener(i)=statev(iv)
            iv=iv+1
            DRAG(i)=statev(iv)
            iv=iv+1
            EDOT0(i)=statev(iv)
            iv=iv+1
            rho_inits(i)=statev(iv)
            iv=iv+1
            Atau0(i)=statev(iv)
            iv=iv+1
            Atau1(i)=statev(iv)
            iv=iv+1
            Atau2(i)=statev(iv)
            iv=iv+1
            Fdeb0(i)=statev(iv)
            iv=iv+1
            Fdeb1(i)=statev(iv)
            iv=iv+1
            Fdeb2(i)=statev(iv)
            iv=iv+1
            HPcoef(i)=statev(iv)
            iv=iv+1
            T0nucl(i)=statev(iv)
            iv=iv+1
            T1nucl(i)=statev(iv)
            iv=iv+1
            T2nucl(i)=statev(iv)
            iv=iv+1
            T0prop(i)=statev(iv)
            iv=iv+1
            T1prop(i)=statev(iv)
            iv=iv+1
            T2prop(i)=statev(iv)
            if (ntwmmx.gt.0) then
            do j=1,ntwmmx
               iv=iv+1
               CTWSL(J,I)=statev(iv)
               iv=iv+1
               CTS1(J,I)=statev(iv)
               iv=iv+1
               CTS2(J,I)=statev(iv)
               iv=iv+1
               CTS0(J,I)=statev(iv)
            enddo ! end of ntwmmx (j)
            endif

         enddo ! end of nmodmx (i)
            iv=iv+1
            K3recov=statev(iv)
            iv=iv+1
            K4recov=statev(iv)
            iv=iv+1
            D0vac=statev(iv)
            iv=iv+1
            Efvac=statev(iv)
            iv=iv+1
            Emvac=statev(iv)

         do i=1,ngrmx
            iv=iv+1
            rho_deb(i)=statev(iv)
            do j=1,nsysmx
               iv=iv+1
               rho_s_for(j,i)=statev(iv)
               iv=iv+1
               rho_s_deb(j,i)=statev(iv)
            enddo
         enddo
         iv=iv+1
         edotref=statev(iv)
         do i=1,nmodmx
            do j=1,nphmx
               iv=iv+1
               xmu_mode(i,j)=statev(iv)
            enddo
         enddo
      endif                     ! end of dd


c     rot_acum
         do i=1,3
            do j=1,3
               iv=iv+1
               rot_acum(i,j)=statev(iv)
            enddo
         enddo
c     fstar & alplpx
         iv=iv+1
         fstar=statev(iv)
         iv=iv+1
         alplpx=statev(iv)

c     acctimefe
         iv=iv+1
         acctimefe=statev(iv)
c     kount_umat_call
         iv=iv+1
         kount_umat_call=int(statev(iv))

ccc   polyfits of cijv, cte
c     with ndegs_fit
c         b) ndegs_fit (nphmx,2)
         do i=1,nphmx
         do j=1,2
               iv=iv+1
               ndegs_fit(i,j)=int(statev(iv))
         enddo
c         c) cijv_fit, cte_fit
         do j=0,mxndeg
            do k=1,6
               do l=1,6
                  iv=iv+1
                  cijv_fit(i,j,k,l)=statev(iv)
               enddo
               iv=iv+1
               cte_fit(i,j,k)=statev(iv)
            enddo
         enddo
         enddo


         if(iverb.eq.1) then
         write(*,*)
         write(*,*)'---------------------'
         write(*,*)'iopt.eq.2 ****iv:',iv
         write(*,*)'---------------------'
         endif

         if (ihardlaw.eq.23.and..not.(irgvb)) then
            write(*,*)'**Error: Conflicts between ihardlaw and irgvb'
            call xit()
         endif
         if (ihardlaw.eq.20.and..not.(idd)) then
            write(*,*)'**Error: Conflicts between ihardlaw and idd'
            call xit()
         endif

      elseif(iopt.eq.3)then
c           ( 3): Save state variables to a file
c         write(*,*)'iofn:',iofn
c         write(*,*)'total nstatv:',nstatv
         write(iofn,*)(statev(i),i=1,nstatv)
      elseif(iopt.eq.4)then
c           ( 4): Load state variables from a file
c         write(*,*)'iofn:',iofn
c         write(*,*)'total nstatv:',nstatv
         read(iofn,*)(statev(i),i=1,nstatv)
      else
         write(*,*)'iopt:',iopt
         write(*,*)'nothing happened in statev_io'
         call xit()

      endif
      return
      end subroutine statev_io
