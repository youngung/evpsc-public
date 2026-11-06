c======================================================================c
c     Subroutine <rexdd_init> initialize the dislocation density of
c     given grain id 'kch' belonging to phase id (iph)
      subroutine rexdd_init_state(
c     intent(in)
c          * integer
     $     iorilaw,iregime,mdim,nmodmx,ntwmmx,nsysmx,
     $     nmodes,ntwmod,
     $     nslmod,nsm,nslsys,nsyst,isense,nrs
c          * real
     $     ,schca,nschca,gamd0s_par
     $     ,sbar,dsbar,dbar,dbar0,temp,celcc,radius_chi
     $     ,agpar,mgel,mgtgpar
     $     ,dg0par,sgpar,dsgpar,sgtrypar,dsgtrypar

     $     ,RHO_initS,RHO_initW,mu_mode
     $     ,atau0,atau1,atau2,chifor,chisub,burg,HPcoef
     $     ,ctwsl,RHOSAT,tau_prop,tau_nucl
c     intent(inout)
c          * integer
c          * real
     $     ,mgtgchi,dg0chi
     $     ,sgchi,dsgchi,sgtrychi,dsgtrychi
     $
c     intent(out)
     $     ,agchi,fijgrch
     $     ,rho_s_for_chi,rho_s_deb_chi,crss_chi
     $     )
      implicit none
c     intent(in)
      integer, intent(in):: iorilaw,iregime,mdim,
     $     nmodmx,ntwmmx,nsysmx,nmodes,ntwmod,nslmod,
     $     nsm(nmodmx),nslsys,nsyst,isense(nsysmx),nrs(nsysmx)
      real*8, intent(in)::
     $     schca(5,nsysmx)
     $     ,nschca(6,nsysmx)
     $     ,gamd0s_par(nsysmx),

     $     sbar(6),dsbar(6),dbar(6),dbar0(6),
     $     temp,celcc(6,6),radius_chi,
     $     agpar(3,3),mgel(6,6),mgtgpar(mdim,mdim),
     $     dg0par(mdim),
     $     sgpar(6),dsgpar(6),sgtrypar(6),dsgtrypar(6),
     $     RHO_initS(nmodmx),RHO_initW,mu_mode(nmodmx),
     $     atau0(nmodmx),atau1(nmodmx),atau2(nmodmx),
     $     chifor,chisub,burg(nmodmx),HPcoef(nmodmx),
     $     ctwsl(ntwmmx,nmodmx),
     $     RHOSAT(nmodmx),tau_prop(nmodmx),tau_nucl(nmodmx)
c     intent(inout)
      real*8, intent(inout)::mgtgchi(mdim,mdim),dg0chi(mdim),
     $     sgchi(6),dsgchi(6),sgtrychi(6),dsgtrychi(6)
c     intent(out)
      real*8, intent(out):: agchi(3,3),fijgrch(3,3),
     $     rho_s_for_chi(nsysmx),rho_s_deb_chi(nsysmx),
     $     crss_chi(nsysmx)

c     locals
      integer i,j,k,ism,jsm,kmo,kmox,JMO,
     $     is,js,interactionx,ishapex,idiag_sg,idiag_nr
      real*8 vnorm,dbarnorm,aa(3,3)
      real*8 atau,rho_g_for,rho_g_deb,fact,
     $     bux33(3,3),aux6(6),bux6(6),aux33(3,3),aux66(6,6),bux66(6,6),
     $     auxm(mdim),auxmm(mdim,mdim),aux3333(3,3,3,3),
     $     schsa(5,nsysmx),nschsa(6,nsysmx)

c     ---yj
      real*8 edot,
     $     DASDBARx(mdim),SASTAVx(mdim),
     $     DG_TRANSx(6),
     $     SBAR_OLD(6),
     $     SG_OLDx(6),
     $     MASTPHx(6,6),MASTGRx(6,5),
     $     c6,dtime
c     ---yj
      real*8 RHOnucl



      aux6(:)=sbar(:)
      aux6(:)=dsbar(:)
c      aux6(:)=dbar(:)
      aux6(:)=dbar0(:)
      aux6(:)=sgtrypar(:)

c     edot is used for CRSS determination of twins
      EDOT=SQRT(2./3.)*VNORM(DBAR,5) !! caution when dbar = 0.

c     Orientation (random orientation or inherits parent orientation?)
      if (iorilaw.eq.0) then
         agchi(:,:)=agpar(:,:)
      elseif(iorilaw.eq.1) then
         call rexdd_gen_rand_ori(aux33)
         agchi(:,:)=aux33(:,:)
      elseif(iorilaw.eq.2) then

c     if (.true.) then !! testing for plane-strain, rolling condition.
c     call chg_basis(sg(:,kpa),aux33,aux66,aux3333,1,5) !! represent 5D deviatoric stress using 3x3 components
c     sg(6,kpa)=-aux33(2,2)*sqrt(3.)
c     sg(6,kpa)=0.
c      write(*,*)'sg(6,kpa):',sg(6,kpa)
c     sg(6,kpa)=0.
c     endif

         aux6(:)=0.
         do i=1,6               !mdim
         do j=1,6               !mdim
            aux6(i)=aux6(i)+mgel(i,j)*sgpar(j)
         enddo
         enddo
         call rexdd_ori_energy(agpar,celcc,aux6,agchi)
      endif
c     accumulated shear strains (Gamma) in each grain
c$$$      gtotgr(kch)=0. ! maybe not necessary?

c     determine deformation gradient
      fijgrch(:,:)=0d0
      do i=1,3
         fijgrch(i,i)=1d0
      enddo

c        inherits parent's mgtg with considering different crystal orientation
      call rexdd_transf(auxm,mgtgpar,auxm,mgtgchi,agpar,agchi,1,mdim)
c         mgtgchi(:,:)=0.
c        inherits parent's dg0 with considering different crystal orientation
      call rexdd_transf(dg0par,auxmm,dg0chi,auxmm,agpar,agchi,0,mdim)
c         dg0chi(:)=0.
c         dg0chi(1:mdim)=dbar0(1:mdim)


c----------------------------------------------------------------------c
c     Determine the stress state of newly created child.
      if (iregime.eq.1) then
         sgchi(:)=0.
         sgtrychi(:)=0.
         call rexdd_transf(sgpar,aux66,sgchi,bux66,agpar,agchi,0,6)
         call rexdd_transf(sgpar,aux66,sgtrychi,bux66,agpar,agchi,0,6)
c         call rexdd_transf(sgtrypar,aux66,sgtrychi,bux66,agpar,
c     $        agchi,0,6)

c         sgchi(:)=sgpar(:)
c         sgtrychi(:)=sgtrypar(:)
c         sgchi(1:mdim)=sbar(1:mdim)
c         sgtrychi(1:mdim)=sbar(1:mdim)
c        maybe... it is better to use the macroscopic stress as the initial guess?

      elseif(iregime.eq.2) then
         call rexdd_transf(dsgpar,aux66,dsgchi,bux66,agpar,agchi,0,6)
         call rexdd_transf(dsgtrypar,aux66,dsgtrychi,bux66,agpar,
     $        agchi,0,6)

         dsgchi(1:mdim)=dsbar(1:mdim)
         dsgtrychi(1:mdim)=dsbar(1:mdim)
      endif

c     Dislocation density of individual slip systems
      ISM=0
      rho_g_for = 0.
      rho_g_deb = 0.
      DO KMO=1,NSLMOD
      DO IS=1,NSM(KMO)
         ISM=ISM+1
         rho_s_for_chi(ism)=RHO_initS(kmo)
         rho_s_deb_chi(ism)=RHO_initW/nslsys
         rho_g_for=rho_g_for+rho_s_for_chi(ism)
         rho_g_deb=rho_g_deb+rho_s_deb_chi(ism)
      ENDDO
      ENDDO

      ISM=0
      DO KMO=1,NSLMOD
c        Peierls stress
         atau=Atau0(KMO)+Atau1(KMO)*Exp(-TEMP/Atau2(KMO))
         DO IS=1,NSM(KMO)
            ISM=ISM+1
C --> Eq 3.19 & ~3.iph0 B&T --> TAU_0 + TAU_for + TAU_deb
c           1. Peierls friction stress
            crss_chi(ism)=atau
c           2. forest dislocations
            crss_chi(ism)=crss_chi(ism)+
     $           Chifor*mu_mode(kmo)*burg(kmo)
     $           *sqrt(rho_g_for)
c           3. debris (substructure)
            if (rho_g_deb.gt.0) then
               crss_chi(ism)=crss_chi(ism)+
     $              CHIsub*MU_MODE(KMO)*BURG(KMO)
     $              *SQRT(RHO_g_DEB)*LOG(BURG(KMO)
     $              *SQRT(RHO_g_DEB))
c           4. HP effect is missing.
            endif
C --> Eq 3.iph1 B&T; Eq 9 Capolungo et al --> TAU_0 + TAU_for + TAU_deb + TAU_hp
            CRSS_chi(ISM)=crss_chi(ISM)
     $           + HPcoef(KMO)*MU_MODE(KMO)
     $           *SQRT(BURG(KMO)/(2.*radius_chi))
         ENDDO                  ! end of IS
      ENDDO                     ! end of KMO
c$$$
      IF(NTWMOD.NE.0) THEN
c        KMO : loop from the end of slip mode over twin modes
c        KMOX: loop from 1                    over twin modes
         DO KMO=NSLMOD+1,NMODES
            KMOX=KMO-NSLMOD
         DO IS=1,NSM(KMO)
            ISM=ISM+1
            crss_chi(ism)=0.
            RHOnucl=0.
            JSM=0
            DO JMO=1,NSLMOD
            DO JS=1,NSM(JMO)
               JSM=JSM+1
               CRSS_chi(ISM)=MU_MODE(JMO)*BURG(KMO)*BURG(JMO)*
     $              CTWSL(KMOX,JMO)*RHO_S_FOR_CHI(JSM)
               RHOnucl=RHOnucl+
     $              RHO_S_FOR_CHI(JSM)/RHOSAT(JMO)
            ENDDO
            ENDDO
            crss_chi(ism)=crss_chi(ism)+tau_prop(kmo)+
     $           (tau_nucl(kmo)-tau_prop(kmo))*exp(-RHOnucl)
         ENDDO ! is
         ENDDO ! kmo
      ENDIF


c----------------------------------------------------------------------c
c     Obtain schmid tensor, with applying the new orientation of
c     the nucleated child.

c     ## A matrix that transforms from sample to child's crystal axes.
      do i=1,3
      do j=1,3
         aa(j,i)=agchi(i,j)
      enddo
      enddo
      call rexdd_schca2sa(5,nsysmx,nsyst,schca,aa,schsa)
      call rexdd_schca2sa(6,nsysmx,nsyst,nschca,aa,nschsa)


      if (iregime.eq.1) then
c        Let's give it a better guess on grain stress?
         interactionx=0
         ishapex=0

         DASDBARx(:)=0.
         SASTAVx(:)=0.

         MASTPHx(:,:)=0.
         MASTGRx(:,:)=0.
         DG_TRANSX(:)=0.

         SBAR_OLD(:)=0.
         SG_OLDX(:)=0.

         c6=0.
         dtime=1. ! not being used for VPSC.

         sgtrychi(:)=sgtrychi(:)/3.

c         dbarnorm=vnorm(dbar,6)
c         SGTRYchi(:)=SBAR(:)/2.

         call grain_stress_evp(
c        in
     $        interactionx,mdim,
     $        iregime,ishapex,isense,nrs,
     $        nsyst,nsysmx,

     $        DASDBARx,SASTAVx,
     $        DBAR,DBAR0,DG_TRANSx,dtime,
     $        SBAR,SBAR_OLD,SG_OLDx,

     $        MASTPHx,MASTGRx,

     $        schsa,nschsa,     !SCHx,NSCHx,
     $        c6,CRSS_chi,
     $        gamd0s_par,       !GAMD0Sx,
     $        mgel,             !MGELx,
     $        dsgtrychi,sgtrychi, !dSGTRYx,SGTRYx,
c     out
     $        dsgchi,           !dSGx,
     $        sgchi,            !SGx,
     $        k,idiag_nr,idiag_sg)

      endif

      return
      end subroutine rexdd_init_state

c======================================================================c
c     Random initial guess for gradient descending algorithm.
      subroutine rexdd_gen_rand_ori(a)
      implicit none
      real*8, intent(out):: a(3,3)
c     locals
      real*8 phi1,phi,phi2
      real*8 mxphi1,mxphi,mxphi2
      parameter(mxphi1=360.,mxphi=180.,mxphi2=360.)
c      call rexdd_gen_rand_eul(360d0,180d0,360d0,phi1,phi,phi2)
      call rexdd_gen_rand_eul(mxphi1,mxphi,mxphi2,phi1,phi,phi2)
      call euler(2,phi1,phi,phi2,a)
      return
      end subroutine
c======================================================================c
c     Generate random orientation and return the corresponding Euler angles.
      subroutine rexdd_gen_rand_eul(p1mx,pmx,p2mx,phi1,phi,phi2)
      implicit none
      real*8, intent(in):: p1mx,pmx,p2mx
      real*8, intent(out):: phi1,phi,phi2
c     locals
      real*8 random1,random2,random3,pi
      pi=4.*atan(1.)
      call random_number(random1)
      call random_number(random2)
      call random_number(random3)
      phi1= p1mx*random1
      phi =acos(random2)*pmx/pi
      phi2= p2mx*random3
      return
      end subroutine rexdd_gen_rand_eul
