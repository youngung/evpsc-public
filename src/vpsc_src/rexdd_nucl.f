c     Argument order: Numerical CONTROL & CONDITION, DIMENSION, PX, PHASE, GRAIN, SYSTEM
c     Then, if applicable, alphabetical order is used.
      subroutine rexdd_nucl(
c          intent(in):
c          * integer
     $     iregime,istep,mdim,nphmx,ngrmx,ngidmx,nsysmx,nmodmx,ntwmmx,
     $     nmodes,ntwmod,nslmod,nsm,
     $     nslsys,
     $     iph,nclustermx,nchildmx,ncluster,uw_rex_nuco,

     $     nsystph,isenseph,nrsph,
c          * real
     $     schca,               !sch(5,nsysmx)
     $     nschca,              !nsch(6,nsysmx)
     $     gamd0s,              !gamd0s(nsysmx,ngrmx)

     $     rexnuclpar,dtime,temp,vpx,
     $     sbar,dsbar,dbar,dbar0,
     $     epsvm,celcc,
     $     RHO_initS,RHO_initW,mu_mode,
     $     atau0,atau1,atau2,chifor,chisub,burg,HPcoef,

     $     ctwsl,RHOSAT,tau_prop,tau_nucl,

c          intent(inout):
c          * integer
     $     nchildset,ngr,
     $     idclusterset,
     $     idclusterset_reverse,
     $     ngid_total,
     $     gid_of_slot,slot_of_gid,
     $     gid_birth_step,gid_death_step,gid_parent,
c          * real
     $     wgt,wgt_old,ag,sg,dsg,sgtry,dsgtry,dg0,
     $     fijgr,mgel,mgtg,
     $     rho_s_for,rho_s_deb,crss,

c          intent(out):
c          * integer
     $     nkount,nuc_kount,
c          * real
     $     dwgt_nucl,pnucl,radius,
c          * logical
     $     inucleated_px,inucleated_cl
     $
     $     )
      implicit none

      integer, intent(in)::iregime,istep,mdim,
     $     nphmx,ngrmx,ngidmx,nsysmx,nmodmx,ntwmmx,
     $     nmodes,ntwmod,
     $     nslmod(nphmx),nsm(nmodmx,nphmx),nslsys(nphmx),
     $     iph,nclustermx,nchildmx,ncluster,uw_rex_nuco


     $     ,nsystph
     $     ,isenseph(nsysmx)
     $     ,nrsph(nsysmx)

      real*8, intent(in):: rexnuclpar(5),dtime,temp,vpx,
     $     sbar(6),dsbar(6),dbar(6),dbar0(6),
     $     epsvm,
     $     celcc(6,6),
     $     RHO_initS(nmodmx),RHO_initW,mu_mode(nmodmx),
     $     atau0(nmodmx),atau1(nmodmx),atau2(nmodmx),
     $     chifor,chisub,burg(nmodmx),HPcoef(nmodmx),

     $     ctwsl(ntwmmx,nmodmx),
     $     RHOSAT(nmodmx),tau_prop(nmodmx),tau_nucl(nmodmx),

     $     schca(5,nsysmx),
     $     nschca(6,nsysmx),
     $     gamd0s(nsysmx,ngrmx)

      integer, intent(inout)::nchildset(nclustermx),ngr(0:nphmx),
     $     idclusterset(nclustermx,0:nchildmx),
     $     idclusterset_reverse(ngrmx),
     $     ngid_total,
     $     gid_of_slot(ngrmx),
     $     slot_of_gid(ngidmx),
     $     gid_birth_step(ngidmx),
     $     gid_death_step(ngidmx),
     $     gid_parent(ngidmx)
      real*8, intent(inout):: wgt(ngrmx),wgt_old(ngrmx),ag(3,3,ngrmx),
     $     sg(6,ngrmx),dsg(6,ngrmx),sgtry(6,ngrmx),dsgtry(6,ngrmx),
     $     dg0(6,ngrmx),
     $     fijgr(3,3,ngrmx),mgel(6,6,ngrmx),
     $     mgtg(6,6,ngrmx),
     $     rho_s_for(nsysmx,ngrmx),rho_s_deb(nsysmx,ngrmx),
     $     crss(nsysmx,ngrmx)

      integer, intent(out)::nkount(nclustermx),nuc_kount
      real*8, intent(out):: dwgt_nucl(ngrmx),pnucl(ngrmx),
     $     radius(ngrmx)
      logical inucleated_cl(nclustermx),inucleated_px

c     locals
      real*8 pnu                ! nucleation probability
      real*8 rho_cri            ! critical disl. density
      real*8 ndot               ! nucleation frequency
      real*8 lamb               ! unitless time for poisson prob. dist.
      real*8 c,Qp,boltz,rgas,r0_nuc,rand,pi
      integer i,kkk,ism,imo,is,kpa,kpanew,icluster,ich,kch,
     $     iorilaw,ichmx,nch,imem,parent_gid
      real*8 rho_par,rho_par_mx,rho_cur,pnux,rho_g(ngrmx),aux33(3,3),
     $     phi1,phi,phi2,
     $     gamd0s_par(nsysmx)
!     Boltzmann and gas constants
      parameter(boltz=1.380649e-23) ! unit of J/K
      parameter(rgas=8.314462)  ! unit of J/(K.mol)
      integer idgr_nucs(ngrmx),nuc_kount_cl

c     locals particularly for <rexdd_init_state>
      integer nslmodph,nsmph(nmodmx),nslsysph
      real*8 mgelpar(6,6),
     $     agpar(3,3),mgtgpar(mdim,mdim),
     $     agchi(3,3),mgtgchi(mdim,mdim),
     $     dg0par(mdim),sgpar(6),dsgpar(6),sgtrypar(6),dsgtrypar(6),
     $     dg0chi(mdim),sgchi(6),dsgchi(6),sgtrychi(6),dsgtrychi(6),
     $     fijgrchi(3,3),rho_s_forchi(nsysmx),rho_s_debchi(nsysmx),
     $     crsschi(nsysmx)

c----------------------------------------------------------------------c
c     to fool the compiler
c      write(*,*)nsystph,isenseph(1),nrsph(1),
c     $     gamd0s(1,1)
c----------------------------------------------------------------------c

      nslmodph=nslmod(iph)
      nsmph(:)=nsm(:,iph)
      nslsysph=nslsys(iph)

      pi=4d0*dATAN(1d0)

c----------------------------------------------------------------------c
      do icluster=1,ncluster
         if (nchildset(icluster).ge.nchildmx) then
            write(*,*)'** Warning: nchildmx reached; thus cannot create'
            write(*,*)'   child nucleus any more for cluster',icluster
            call xit()
         endif
      enddo
c----------------------------------------------------------------------c
c     Obtain grain dislocation densities (rho_g)
      rho_g(:)=0.
      do kkk=ngr(iph-1)+1,ngr(iph)
         ism=0
         do imo=1,nslmod(iph)
         do is=1,nsm(imo,iph)
            ism=ism+1
            rho_g(kkk)=rho_g(kkk)+rho_s_for(ism,kkk)+rho_s_deb(ism,kkk)
         enddo
         enddo
      enddo

c----------------------------------------------------------------------c
c     Make sure that the first child (0) is the one with the maximum
c     disl. density. If not, select the one with the highest disl. den.
      do icluster=1,ncluster
c        Grain ID with the maximum dis. density from the previous step.
         kpa=idclusterset(icluster,0)
         kpanew=kpa

         ichmx=-1
         rho_par_mx = rho_g(kpa)

c        Only when both parent and children exist
         if (nchildset(icluster).gt.0) then
c           Find if any 'child' member has a higher disl. density
            do ich=1, nchildset(icluster)
               kch=idclusterset(icluster,ich)
c              If current child has higher disl. density,
c              update <kpanew>, <ichmx>, and <rho_par_mx>
               if (rho_g(kch).gt.rho_par_mx) then
                  kpanew=kch ! new grain that has the maximum dd
                  ichmx=ich
                  rho_par_mx=rho_g(kch)
               endif
            enddo
c           If it turns out that a child (ichmx) has the highest
c           disl. density, switch 'rank' of this child and
c           the current one.
            if (ichmx.ne.-1) then
               idclusterset(icluster,0)=kpanew
               idclusterset(icluster,ichmx)=kpa
            endif
         endif
      enddo
c----------------------------------------------------------------------c
c      call random_seed()
c     nucleation parameters
      r0_nuc=rexnuclpar(1)
      rho_cri=rexnuclpar(2)
      ndot=rexnuclpar(3)
      Qp=rexnuclpar(4)
      iorilaw=int(rexnuclpar(5))

c----------------------------------------------------------------------c
c     Main loop over cluster begins here.
c     Each cluster is tested if the nucleation criterion is met.

c     ndot: nucleation frequency per unit grain boundary area
c     nucleation frequency according to Ding and Guo
c     pnux: nucleation probability
      pnux=exp(-Qp/(rgas*temp))
      dwgt_nucl(:)=0.
      inucleated_px=.false.     ! A global boolean var.
      inucleated_cl(:)=.false.  ! A local boolean flag array
      nuc_kount=0               ! number of nucleated grains for the current time increment


      do icluster=1,ncluster
c        Initial number of members within this cluster
c        before nucleating any additional member
         nch=nchildset(icluster)
         nuc_kount_cl=0
c         nch=0 !! to force only parent to nucleate?
         do imem=0,nch ! Allows any grain to nucleate a child
            kpa=idclusterset(icluster,imem)
            rho_par=rho_g(kpa)
            pnu=pnux*exp(-rho_cri/rho_par)
c           lamb: unitless time for poisson distribution
            if (.False.) then
c           if (.True.) then
               lamb = ndot * dtime * pi * radius(kpa) **2
            else
               lamb = ndot * dtime
            endif
            pnucl(kpa)=1.-exp(-pnu*lamb) !! according to Morris and Klahn, J. App. Phy. (1973)
c----------------------------------------------------------------------c
c       Main loop to determine nucleated grains and its initial properties
c        The nucleation criterion is tested.
c        disl. density in each parent

            if (pnucl(kpa).lt.1e-30) pnucl(kpa)=0.

            call random_number(rand)
            if ((pnucl(kpa).ge.rand)
     $           .and.(nchildset(icluster).lt.nchildmx-1)
     $           ) then
            inucleated_px=.true.
            inucleated_cl(icluster)=.true.

c           Check permanent GID capacity
            if (ngid_total.ge.ngidmx) then
               write(*,*)'** Error: ngid_total reached ngidmx'
               write(*,*)'ngid_total:',ngid_total
               write(*,*)'ngidmx    :',ngidmx
               call xit()
            endif

c           The next permanent GID must still be in NOT-BORN state.
c           This check prevents accidental reuse of an issued GID.
            if (gid_death_step(ngid_total+1).ne.-2) then
               write(*,*)'** Error: next GID is not unused'
               write(*,*)'next GID   :',ngid_total+1
               write(*,*)'death state:',
     $              gid_death_step(ngid_total+1)
               call xit()
            endif

c           Check dynamic grain slot capacity
            if (ngr(iph).ge.ngrmx) then
               write(*,*)'** Error: no dynamic grain slot available'
               call xit()
            endif

c           Remember permanent GID of the parent
            parent_gid=gid_of_slot(kpa)

            if (parent_gid.le.0) then
               write(*,*)'** Error: invalid parent GID'
               write(*,*)'kpa:',kpa
               call xit()
            endif

c           nuc_kount_cl: number of nucleated grains for the current cluster
c           nuc_kount   : number of nucleated grains for the current time increment
c           nkidstot    : total accumulated number of nucleated grains
c           nkount      : total accumulated number of nucleated grains for each cluster
            nuc_kount_cl=nuc_kount_cl+1
            nuc_kount=nuc_kount+1
c            nkidstot=nkidstot+1 (deprecated)
            nkount(icluster)=nkount(icluster)+1
c           Update the number of children in the current cluster i
            nchildset(icluster)=nchildset(icluster)+1

c           Update the pointer array for members
c           Extend the bottom stack of grain id pointer
c           and the bottom-most slot is assigned to the nucleated grain

            ngr(iph)=ngr(iph)+1 !! this might need to be carefuly revised
c
            kch=ngr(iph)
            idgr_nucs(nuc_kount)=kch
            idclusterset(icluster,nchildset(icluster))=kch
            idclusterset_reverse(kch)=icluster

c----------------------------------------------------------------------c
c           Assign permanent GID to newly nucleated grain
            ngid_total=ngid_total+1

c           Connect dynamic slot <kch> and new permanent GID.
            gid_of_slot(kch)=ngid_total
            slot_of_gid(ngid_total)=kch

c           Change permanent GID state:
c             NOT BORN (-2) -> ALIVE (-1)
            gid_birth_step(ngid_total)=istep
            gid_death_step(ngid_total)=-1

c           Store permanent GID of the parent grain.
            gid_parent(ngid_total)=parent_gid

c----------------------------------------------------------------------c

c           Determining the initial state of nucleated grain
c           Inherits the properties of parent grain to kch
            wgt(kch)=0.
            wgt_old(kch)=0.
c$$$            dwgt_grow(kch)=0.
            dwgt_nucl(kch)=4*pi/3./vpx*(r0_nuc)**3
            radius(kch)=r0_nuc
            dwgt_nucl(kpa)=-dwgt_nucl(kch)

c
            mgelpar(:,:)=mgel(:,:,kpa)

            agpar=ag(:,:,kpa)
            mgtgpar(:,:)=mgtg(1:mdim,1:mdim,kpa)
            dg0par(:)=dg0(1:mdim,kpa)
            sgpar(:)=sg(:,kpa)
            dsgpar(:)=dsg(:,kpa)
            sgtrypar(:)=sgtry(:,kpa)
            dsgtrypar(:)=dsgtry(:,kpa)

c            schcaph(:,:)
c            nschcaph(:,:)
            gamd0s_par(:)=gamd0s(:,kpa)

            mgtgchi(:,:)=mgtg(1:mdim,1:mdim,kch)
            dg0chi(:)=dg0(1:mdim,kch)
            sgchi(:)=sg(:,kch)
            dsgchi(:)=dsg(:,kch)
            sgtrychi(:)=sgtry(:,kch)
            dsgtrychi(:)=dsgtry(:,kch)

            call rexdd_init_state(
c     intent(in)
c                * integer
     $           iorilaw,iregime,mdim,nmodmx,ntwmmx,
     $           nsysmx,
     $           nmodes,ntwmod,
     $           nslmodph,nsmph,nslsysph

     $           ,nsystph
     $           ,isenseph      ! (nsysmx)
     $           ,nrsph         ! (nsysmx)
c$$$c          * real
     $           ,schca       !sch(5,nsysmx)
     $           ,nschca        !nsch(6,nsysmx)
     $           ,gamd0s_par        !gamd0s(nsysmx,ngrmx)

     $           ,sbar,dsbar,dbar,dbar0,temp,celcc,radius(kch)

     $           ,agpar,mgelpar,mgtgpar
     $           ,dg0par,sgpar,dsgpar,sgtrypar,dsgtrypar

     $           ,RHO_initS,RHO_initW,mu_mode
     $           ,atau0,atau1,atau2,chifor,chisub,burg,HPcoef

     $           ,ctwsl
     $           ,RHOSAT,tau_prop,tau_nucl
c     intent(inout)
c                * real
     $           ,mgtgchi,dg0chi,sgchi,dsgchi,sgtrychi,dsgtrychi
c     intent(out)
c                * real
     $           ,agchi,fijgrchi,rho_s_forchi,rho_s_debchi,crsschi
     $           )

c           intent:: inout
            mgtg(1:mdim,1:mdim,kch)=mgtgchi(:,:)
            dg0(1:mdim,kch)=dg0chi(:)
            sg(:,kch)=sgchi(:)
            dsg(:,kch)=dsgchi(:)
            sgtry(:,kch)=sgtrychi(:)
            dsgtry(:,kch)=dsgtrychi(:)

c           intent:: out
            ag(:,:,kch)=agchi(:,:)
            fijgr(:,:,kch)=fijgrchi(:,:)
            rho_s_for(:,kch)=rho_s_forchi(:)
            rho_s_deb(:,kch)=rho_s_debchi(:)
            crss(:,kch)=crsschi(:)

c$$$            call neighbours(0)
         endif                  ! If condition for nucleation is met
         enddo ! do loop over initial members
      enddo                     ! Do loop over clusters
c----------------------------------------------------------------------c
c     Download Euler angles of nucleated grains
      do i=1,nuc_kount
         kkk=idgr_nucs(i)
         aux33(:,:)=ag(:,:,kkk)
         call euler(1,phi1,phi,phi2,aux33)
         if (i.eq.1) then
c           Download to "REXDD_NUCL_ORI.OUT"
            write(uw_rex_nuco,*)'------------'
            write(uw_rex_nuco,'(i7,f8.3)')nuc_kount,epsvm
         endif
         write(uw_rex_nuco,'(2i7,4e11.3)')
     $        i,kkk,phi1,phi,phi2,dwgt_nucl(kkk)
      enddo

c----------------------------------------------------------------------c
c$$$c     Block for checking volume transfer after nucleation
c$$$      if (inucleated_px.and..false.) then
c$$$c      if (inucleated_px.and..true.) then
c$$$         write(*,*)'nucleated'
c$$$         do icluster=1,ncluster
c$$$            if (inucleated_cl(icluster)) then
c$$$               call rexdd_display(1,icluster)
c$$$            endif
c$$$         enddo
c$$$      endif
      return
      end subroutine rexdd_nucl
