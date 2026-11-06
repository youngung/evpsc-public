c     Argument order: Numerical CONTROL & CONDITION, DIMENSION, PX, PHASE, GRAIN, SYSTEM
c     Then, if applicable, alphabetical order is used.

c======================================================================c
c     Recrystalization model based on disl. density and clusters of
c     grains. Within cluster, the grains complete against each other
c     to grow. Consequently, some grains will grow by consuming others
c     while grains with high dislocation density can shrink by being
c     consumed.
      subroutine rexdd(
c          intent(in):
c          * integer
     $     iopt,iregime,mdim,iproc,istep,nphmx,nph,iph,
     $     nclustermx,nchildmx,nchildlimit,
     $     ngrmx,ngidmx,nsysmx,nmodmx,ntwmmx,nmodes,ntwmod,
     $     nslmod,nsm,nslsys,
     $     uw_rex,UW_REX_STAT,uw_rex_nuco,uw_rex_grow,

     $     nsystph,isenseph,nrsph,
c          * real
     $     schca,nschca,gamd0s,
     $     time,dtime,temp,rexgrowpar,rexnuclpar,sbar,dsbar,dbar,dbar0,
     $     epsvm,svm,celcc,
     $     RHO_initS,RHO_initW,mu_mode,
     $     atau0,atau1,atau2,chifor,chisub,burg,HPcoef,

     $     ctwsl,RHOSAT,tau_prop,tau_nucl,
c          intent(inout):
c          * integer
     $     ncluster,nchildset,ngr,dngracc_rex,
     $     ngid_total,gid_of_slot,slot_of_gid,
     $     gid_birth_step,gid_death_step,gid_parent,
c          * real
     $     wgt,wgt_old,ag,sg,dsg,sgtry,dsgtry,sg_old,
     $     fijgr,axisgr,dg_trans,gtotgr,
     $     mgel,mgtg,dg0,radius,
     $     rho_s_for,rho_s_deb,crss,

c          intent(out):
c          * integer
     $     ngref,igref,nkount,
     $     idclusterset,idclusterset_reverse,
c          * real
     $     vpx,rexdd_tiny,dwgt_grow,dwgt_nucl,
c          * logical
     $     inucleated_px,inucleated_cl
     $     )
      implicit none
      integer,intent(in)::
     $     iopt,iregime,mdim,iproc,istep,nphmx,nph,iph,
     $     nclustermx,nchildmx,nchildlimit,ngrmx,ngidmx,
     $     nsysmx,nmodmx,ntwmmx,nmodes,ntwmod,
     $     nslmod(nphmx),nsm(nmodmx,nphmx),
     $     nslsys(nphmx),
     $     uw_rex,uw_rex_stat,uw_rex_nuco,uw_rex_grow

     $     ,nsystph
     $     ,isenseph(nsysmx)
     $     ,nrsph(nsysmx)


      real*8, intent(in)::
     $     schca(5,nsysmx),
     $     nschca(6,nsysmx),
     $     gamd0s(nsysmx,ngrmx),

     $     time,dtime,temp,
     $     rexgrowpar(5),rexnuclpar(5),
     $     sbar(6),dsbar(6),dbar(6),dbar0(6),
     $     epsvm,svm,celcc(6,6),
     $     RHO_initS(nmodmx),RHO_initW,mu_mode(nmodmx),
     $     atau0(nmodmx),atau1(nmodmx),atau2(nmodmx),
     $     chifor,chisub,burg(nmodmx),HPcoef(nmodmx),
     $
     $     ctwsl(ntwmmx,nmodmx),
     $     RHOSAT(nmodmx),tau_prop(nmodmx),tau_nucl(nmodmx)

      integer,intent(inout)::ncluster,
     $     nchildset(nclustermx),ngr(0:nphmx),
     $     dngracc_rex,
     $     ngid_total,
     $     gid_of_slot(ngrmx),
     $     slot_of_gid(ngidmx),
     $     gid_birth_step(ngidmx),
     $     gid_death_step(ngidmx),
     $     gid_parent(ngidmx)
      real*8, intent(inout)::wgt(ngrmx),wgt_old(ngrmx),ag(3,3,ngrmx),
     $     sg(6,ngrmx),dsg(6,ngrmx),sgtry(6,ngrmx),dsgtry(6,ngrmx),
     $     sg_old(6,ngrmx),
     $     fijgr(3,3,ngrmx),axisgr(0:3,3,ngrmx),dg_trans(6,ngrmx),
     $     gtotgr(ngrmx),
     $     mgel(6,6,ngrmx),mgtg(6,6,ngrmx),dg0(6,ngrmx),
     $     radius(ngrmx),
     $     rho_s_for(nsysmx,ngrmx),rho_s_deb(nsysmx,ngrmx),
     $     crss(nsysmx,ngrmx)

      integer,intent(out)::ngref,igref(ngrmx),nkount(nclustermx),
     $     idclusterset(nclustermx,0:nchildmx),
     $     idclusterset_reverse(ngrmx)
      real*8, intent(out)::dwgt_grow(ngrmx),dwgt_nucl(ngrmx),
     $     vpx,rexdd_tiny
      logical, intent(out):: inucleated_px,inucleated_cl(nclustermx)
c     locals
      integer i,j,imo,is,ism,imem,icluster,nuc_kount,
     $     nkidstot,kkk,kpa,kch,kount_van_all,
     $     nrem,idvanish(ngrmx),list(ngrmx),neach,nums(ngrmx),
     $     ngrph_k,ngrph_kp1
      logical iadd,icheck,ivanish,iwritegrow
      real*8 rho_g(ngrmx),velgrow(ngrmx),pnucl(ngrmx)
     $     ,drxfrac,pi,fact,wgtsum,crss1(ngrmx),wdot(ngrmx),
     $     avgmem
      character*5 cistep

c      parameter(iwritegrow=.true.)
c      parameter(iwritegrow=.false.)
      integer iwritecl
      parameter(iwritecl=1)


      ngrph_k=ngr(iph)-ngr(iph-1)
      pi=4d0*dATAN(1d0)
c     Initialization for REXTAL and verify the dimension limits
      if (iproc.eq.1.and.istep.eq.1) then
c----------------------------------------------------------------------c
c        Volume of entire polycrystal (VPX) to correlate
c        the volume fraction (wgt) to grain size (radius)
         vpx=ngr(nph)*4./3.*pi*(radius(1)**3)

         if (nclustermx*nchildmx.gt.ngrmx) then
            write(*,*)'nclustermx:',nclustermx
            write(*,*)'nchildmx:',nchildmx
            write(*,*)'ngrmx:',ngrmx
            write(*,*)'** Error: nclustermx*nchildmx.gt.ngrmx'
            write(*,*)'          ngrmx should be increased to ',
     $           nclustermx*nchildmx
            call xit()
         endif
c----------------------------------------------------------------------c
c     Files declared for the initial process & step
         open(UW_REX,file='REXDD.OUT',status='UNKNOWN')
         write(UW_REX,'(5a10)',advance='no')'time','temp',
     $        'Evm','Svm','DRXVF'
         write(UW_REX,'(2a7)',advance='no')'Pavg','Pstd'
         write(UW_REX,'(5a7)',advance='no')
     $        'Ngr','Ncl','knuc','nvan','avgmem'
         write(UW_REX,'(4a11)',advance='no')
     $        'RavgA','RavgW','RstdA','RstdW'
         write(UW_REX,'(4a11)',advance='no')
     $        'VavgA','VavgW','VstdA','VstdW'
         write(UW_REX,'(2a11)',advance='no')
     $        'RHOavgW','RHOstdW'
         write(UW_REX,'(6a11)',advance='no')'CRSSmn','CRSSmx',
     $        'CRSSavgW','CRSSstdW','CRSSavgA','CRSSstdA'
c        line breaker
         write(UW_REX,*)

         open(UW_REX_STAT,file='REXDD_STATS.OUT',status='unknown')
         open(UW_REX_NUCO,file='REXDD_NUCL_ORI.OUT',status='UNKNOWN')
c$$$         open(120,file='REXDD_HADAMARD.OUT',status='unknown')
c----------------------------------------------------------------------c
c     Initialize permanent grain ID tracking
         if (ngr(1).gt.ngidmx) then
            write(*,*)'** Error: initial ngr > ngidmx'
            call xit()
         endif

         ngid_total=ngr(1)

c        --------------------------------------------------------------
c        Permanent GID state convention
c
c        gid_birth_step = -1 : this GID has never been issued
c
c        gid_death_step = -2 : NOT BORN / never issued
c        gid_death_step = -1 : ALIVE
c        gid_death_step >= 0 : DEAD at the recorded step
c
c        slot_of_gid = 0      : no active dynamic slot
c        slot_of_gid > 0      : currently active
c        --------------------------------------------------------------

         gid_of_slot(:)=0
         slot_of_gid(:)=0

c        All unused GID entries are initially "not born".
         gid_birth_step(:)=-1
         gid_death_step(:)=-2

         gid_parent(:)=0

         do kkk=1,ngr(1)

c           Initial grains use identical slot and permanent GID.
            gid_of_slot(kkk)=kkk
            slot_of_gid(kkk)=kkk

c           Initial grains already exist before DRX evolution starts.
            gid_birth_step(kkk)=0

c           -1 means this issued GID is currently alive.
            gid_death_step(kkk)=-1

c           Initial grains have no parent GID.
            gid_parent(kkk)=0

         enddo
c----------------------------------------------------------------------c
c        Set up clusters
c        1. First, determine how many clusters are required.
c        2. Distribute the grains across all clusters.

         ncluster=ceiling(ngr(nph)*1./(nchildlimit))
         neach   =ceiling(ngr(nph)*1./ncluster)

         if (ncluster.eq.0) then
            write(*,*)'ncluster is zero:',ncluster
            call stop_debug(0,'** Error: unexpected case')
         endif
         if (ncluster*neach.lt.ngr(nph)) then
            write(*,*)'** Error: unexpected case'
            call stop_debug(0,'** Error: unexpected case')
         endif

c        Initial list of integers
         do kkk=ngr(iph-1)+1,ngr(iph)
            nums(kkk)=kkk
         enddo
c        Random shuffle order
         call shuffle(nums(1:ngr(1)),ngr(1))

c        The change in recrystallized volume fraction is
c        tracked down the initial pristine grains.
c        <ngref> : Number of initial grains.
c        <igref> : Grain ID array of them. This is necessary
c                 because grain stack (thus the grain ID)
c                 can be rearranged in the subroutine
c                 <rexdd_maintain>
         nrem=ngr(1)
         ngref=0
         do i=0,ncluster-1
            icluster=i+1
            nchildset(icluster)=-1
            imem=0
            nkount(icluster)=0
            do j=1,neach
               ngref=ngref+1
               kkk=nums(i*neach+j)
c               write(*,'(2i5)') kkk,nrem
               idclusterset(icluster,imem)=kkk
               igref(ngref)=kkk
               idclusterset_reverse(kkk)=icluster
               nchildset(icluster)=nchildset(icluster)+1
               nrem=nrem-1
               imem=imem+1
               if (nrem.eq.0) goto 10
            enddo
         enddo
c     ** The end of initialization / validation
 10      write(*,*)
      endif ! iproc.eq.1 and istep.eq.1


c----------------------------------------------------------------------c
c     rexdd_tiny: this tiny value is used to determine if
c     a grain has a sufficiently small weight, eligible to be
c     considered effectively as vanished. A grain is set as
c     'vanished' one if abs(wgt)<rexdd_tiny.
      rexdd_tiny=1./ngrph_k*1e-10
c----------------------------------------------------------------------c
      inucleated_px=.false.
      inucleated_cl(:)=.false.
      dwgt_grow(:)=0.
      dwgt_nucl(:)=0.

      if (iph.ne.nph) then
         write(*,*)'** Error: subroutine rexdd_nucl updates ngr pointer'
         write(*,*)'    and the algorithm does not work for multi phase'
         call xit()
      endif

c        Nucleation
      call rexdd_nucl(
c          intent(in):
c          * integer
     $     iregime,istep,mdim,nphmx,ngrmx,ngidmx,nsysmx,nmodmx,
     $     ntwmmx,nmodes,ntwmod,nslmod,nsm,
     $     nslsys,iph,nclustermx,nchildmx,ncluster,uw_rex_nuco,

     $     nsystph,isenseph,nrsph,

c          * real
     $     schca,nschca,gamd0s,
     $     rexnuclpar,dtime,temp,vpx,
     $     sbar,dsbar,dbar,dbar0,
     $     epsvm,celcc,
c              - DD law related parameters.
     $     RHO_initS,RHO_initW,mu_mode,
     $     atau0,atau1,atau2,chifor,chisub,burg,HPcoef,

     $     ctwsl,RHOSAT,tau_prop,tau_nucl,
c              - particularly for twins?
c               1. <ctwsl>
c               2. I need, RHOnucl, which is function of
c                  2-1 rho_s_for,
c                  2-2 <RHOSAT>
c               3. I need <tau_prop>
c               4. I need <tau_nucl>

c       Overall, I need, 0. ntwmmx                   - done
c                        1. ctwsl(ntwmmx,nmodmx)
c                        2. RHOSAT(nmodmx)
c                        3. tau_prop(nmodmx)
c                        4. tau_nucl(nmodmx)


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
     $     )

c     rho_s_for and rho_s_deb are updated after nucleation.
c     Accordingly, rho_g is updated in below:
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

c        Growth
c         iwritegrow=.false.
         iwritegrow=.true.
c         if (istep.eq.100) iwritegrow=.true.

         if (iwritegrow) then
            write(cistep,'(i5.5)') istep
            if (iopt.eq.0) then
               open(UW_REX_GROW,file='REXDD_DRX_GROW_'//cistep//'.OUT',
     $              status='UNKNOWN')
            elseif(iopt.eq.1) then
               open(UW_REX_GROW,file='REXDD_SRX_GROW_'//cistep//'.OUT',
     $              status='UNKNOWN')
            endif
         endif
      call rexdd_grow(
c          intent(in)::
c          * integer
     $     iwritecl,uw_rex_grow,iph,nphmx,ngrmx,ngr,nchildmx,nmodmx,
     $     nclustermx,ncluster,idclusterset,nchildset,
c          * real
     $     temp,dtime,vpx,rexgrowpar,burg,mu_mode,rho_g,
     $     dwgt_nucl,wgt,
c          intent(inout)::
c          * real
     $     radius,dwgt_grow,
c          intent(out)::
c          * real
     $     velgrow
     $     )
      if (iwritegrow) close(UW_REX_GROW)




c----------------------------------------------------------------------c
c     Some grains may vanish as a result of growth and shrink
      call rexdd_vanish(
c          intent(in)::
c          * integer
     $     istep,nphmx,ngrmx,ngidmx,ngid_total,nsysmx,
     $     nclustermx,nchildmx,ncluster
c          * real
     $     ,rexdd_tiny

c          intent(inout)::
c          * integer
     $     ,ngr,ngref,igref,nkidstot,
     $     idclusterset,idclusterset_reverse,nchildset
     $     ,gid_of_slot,slot_of_gid,gid_death_step
c          * real
     $     ,pnucl,dwgt_nucl,dwgt_grow,velgrow,radius
c          ** grain data & properties
     $     ,ag,wgt,fijgr,axisgr,crss,sg,dsg,sg_old,sgtry,dsgtry
     $     ,dg0,dg_trans,gtotgr,mgtg
c          ** dislocation densities
     $     ,rho_s_for,rho_s_deb,rho_g
c          intent(out)::
     $     ,kount_van_all
     $     )

c----------------------------------------------------------------------c
c     Display grain stack when nucleation or vanishing occurs.
c----------------------------------------------------------------------c
      if (inucleated_px.or.kount_van_all.gt.0) then

         crss1(:)=crss(1,:)

         call rexdd_disp_grs(
     $        nphmx,ngrmx,nclustermx,nchildmx,
     $        ngr,nchildset,idclusterset,idclusterset_reverse,
     $        gid_of_slot,
     $        rexdd_tiny,rho_g,wgt,dwgt_nucl,dwgt_grow,radius,
     $        crss1,pnucl)

      endif


c$$$      if (inucleated_px) then
c$$$         write(*,*)'nucleated?','ngr:',ngr(1)
c$$$         call stop_debug(0,'** stopped for nucleation')
c$$$      endif

c----------------------------------------------------------------------c
c     Statistical analysis
      call rexdd_stats(
c          intent(in):
c          * integer
     $     nphmx,iph,
     $
     $     ngr,
     $     ngrmx,nsysmx,nmodmx,nslmod,nsm,
     $
     $     ncluster,nclustermx,nchildset,
     $     nuc_kount,kount_van_all,
     $     ngref,igref,
     $
     $     uw_rex,
c          * real
     $     time,temp,epsvm,svm,
     $     wgt,
     $     radius,velgrow,
     $     rho_s_for,rho_s_deb,crss,pnucl,

c          intent(out):
c          * integer
c          * real
     $     dwgt_grow,dwgt_nucl,
c          intent(out)::
     $     drxfrac
     $     )

c     --- avg mem
      avgmem=0.
      do icluster=1,ncluster
         avgmem=avgmem+nchildset(icluster)+1
      enddo
      avgmem=avgmem/(ncluster*1.)

      ngrph_kp1=ngr(iph)-ngr(iph-1)
      dngracc_rex=dngracc_rex+ngrph_kp1-ngrph_k
      write(*,'(6a10)')'ngr_0','ngr_(k)','ngr_(k+1)','dNgr','dNgr_acc',
     $     'avgmem'
      write(*,'(5i10,f10.1)')ngrph_kp1-dngracc_rex,
     $     ngrph_k,ngrph_kp1,ngrph_kp1-ngrph_k,dngracc_rex,avgmem

c----------------------------------------------------------------------c
c     weight updates
      do kkk=ngr(iph-1)+1,ngr(iph)
         wgt(kkk)=wgt(kkk)+dwgt_nucl(kkk)+dwgt_grow(kkk)
      enddo

c----------------------------------------------------------------------c
c     weight renormalization to remove any numerical drift
      fact=3./4./pi*vpx
      wgtsum=0.
      do kkk=ngr(iph-1)+1,ngr(iph)
         wgtsum=wgtsum+wgt(kkk)
      enddo
      do kkk=ngr(iph-1)+1,ngr(iph)
         wgt(kkk)=wgt(kkk)/wgtsum
         radius(kkk)=(wgt(kkk)*fact)**(1./3.)
         if (isnan(radius(kkk))) then
            write(*,*)'in subroutine rexdd_main just after rexdd_stats'
            write(*,*)'nan was found in radius for grain',kkk
            write(*,*)'radius(kkk):',radius(kkk)
            call stop_debug(0,'** Stopped in rexdd_main (3)')
         endif
      enddo

c----------------------------------------------------------------------c
c     writes to 'REXDD_STATS.OUT'
      write(UW_REX_STAT,'(a)')'---'
      write(UW_REX_STAT,'(i6,2(a5,e11.3))')ngr(iph)-ngr(iph-1),
     $     'evm:',epsvm,'svm:',svm
      write(UW_REX_STAT,'(a6,3a11)')'KKK','DEN','WGT','R'
      do kkk=ngr(iph-1)+1,ngr(iph)
         write(UW_REX_STAT,'(i6,3e11.3)') kkk,rho_g(kkk),
     $        wgt(kkk),radius(kkk)
      enddo
c      if (inucleated_px) write(*,*)'** Nucleation event occured'
      write(*,'(a10,f9.3,a1)')'Rex. VF:',drxfrac*100,'%'

c----------------------------------------------------------------------c
c     Create new cluster to prevent the number of grains in a cluster
c     from exceeding a certain number.
      call rexdd_maintain
     $     (
c         intent(in)::
c         * integer
     $     ngrmx,nchildmx,nclustermx,nchildlimit
c         * real
     $     ,dwgt_nucl
c         intent(inout)::
c         * integer
     $     ,ncluster,nchildset,idclusterset,idclusterset_reverse
     $     )

      if (istep.eq.1.and.iopt.eq.1.and..false.) then
         write(*,*)'** Resets NGref, IGref, and DRXfrac'
         write(*,*)'   for the first step of SRX'
         ngref=0
         do icluster=1,ncluster
            do imem=0,nchildset(icluster)
               kkk=idclusterset(icluster,imem)
               ngref=ngref+1
               igref(ngref)=kkk
            enddo
         enddo
         drxfrac=0.
      endif

      return
      end subroutine rexdd
