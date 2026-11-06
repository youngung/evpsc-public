c======================================================================c
      subroutine rexdd_vanish(
c          intent(in)::
c          * integer
     $     istep,nphmx,ngrmx,ngidmx,ngid_total,nsysmx,
     $     nclustermx,nchildmx,ncluster
c          * real
     $     ,rexdd_tiny


c          intent(inout)::
c          * integer
     $     ,ngr,ngref,igref,nkidstot,
     $     idclusterset,idclusterset_reverse,nchildset,
     $     gid_of_slot,slot_of_gid,gid_death_step
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
      implicit none
c     intent(in)::
      integer, intent(in)::
     $     istep,nphmx,ngrmx,ngidmx,ngid_total,nsysmx,
     $     nclustermx,nchildmx,ncluster
      real*8, intent(in):: rexdd_tiny
c     intent(inout)::
      integer, intent(inout)::
     $     ngr(0:nphmx),ngref,igref(ngrmx),nkidstot,
     $     idclusterset(nclustermx,0:nchildmx),
     $     idclusterset_reverse(ngrmx),
     $     nchildset(nclustermx),
     $     gid_of_slot(ngrmx),
     $     slot_of_gid(ngidmx),
     $     gid_death_step(ngidmx)
      real*8, intent(inout)::
c          ** REXD
     $     pnucl(ngrmx),dwgt_nucl(ngrmx),dwgt_grow(ngrmx)
     $     ,velgrow(ngrmx),radius(ngrmx)
c          ** grain data & properties
     $     ,ag(3,3,ngrmx),wgt(ngrmx),fijgr(3,3,ngrmx)
     $     ,axisgr(0:3,3,ngrmx),crss(nsysmx,ngrmx)
     $     ,sg(6,ngrmx),dsg(6,ngrmx),sg_old(6,ngrmx),sgtry(6,ngrmx)
     $     ,dsgtry(6,ngrmx),dg0(6,ngrmx)
     $     ,dg_trans(6,ngrmx),gtotgr(ngrmx),mgtg(6,6,ngrmx)
c          ** dislocation densities
     $     ,rho_s_for(nsysmx,ngrmx),rho_s_deb(nsysmx,ngrmx)
     $     ,rho_g(ngrmx)
c     intent(out)::
      integer, intent(out):: kount_van_all
c     locals
      integer i,j,imem,jmem,idmn,imemmx,kkk,kjj,kkb,idmx,itmp,nvan,
     $     kount_van_px,kount_van_cl,icluster,iclusterb,
     $     nch,idclusterset_vanish(nclustermx),
     $     idgrain_vanish(ngrmx),idmem_vanish(nclustermx),
     $     nvan_cluster(nclustermx),nvan_cluster_id(ngrmx),
     $     neli,ieli(ngrmx),ind,itemp,jnd,
     $     dead_gid,move_gid
      real*8 wgtmp,wgtmp_mn,rho_tmp,rhomx,pi
      logical ivan

      pi=4d0*dATAN(1d0)

c======================================================================c
c     Determine the number of vanishing members in each cluster.
      kount_van_px=0
      kount_van_all=0
      idclusterset_vanish(:)=-1
      idgrain_vanish(:)=-1
      nvan_cluster(:)=0
      nvan_cluster_id(:)=-1
      do 100 icluster=1,ncluster
         nch=nchildset(icluster)
         ivan=.false.

c        Find the cluster stack ID (imemmx) and grain stack ID (idmn) of vanishing grain
         wgtmp_mn=1.0
         idmn=-1                ! grain stack ID
         imemmx=-1              ! cluster stack ID
         do imem=0,nch
            kkk=idclusterset(icluster,imem)
            wgtmp=wgt(kkk)+dwgt_nucl(kkk)+dwgt_grow(kkk)
            if (wgtmp.le.rexdd_tiny)then !! potentially vanishing member
               kount_van_all=kount_van_all+1
               nvan_cluster(icluster)=nvan_cluster(icluster)+1
c              record the stack id of grain in cluster that will vanish
               nvan_cluster_id(kount_van_all)=imem
               ivan=.true.
               if (wgtmp.lt.wgtmp_mn) then
                  wgtmp_mn=wgtmp
                  idmn=kkk
                  imemmx=imem
               endif
            endif
         enddo
c----------------------------------------------------------------------c
         if (ivan) then
            kount_van_px=kount_van_px+1
            idclusterset_vanish(kount_van_px)=icluster
            idgrain_vanish(kount_van_px)=idmn
            idmem_vanish(kount_van_px)=imemmx
         endif
 100  enddo

c======================================================================c
c     If the vanishing member is currently the parent of that cluster,
c     determine new parent (among others) that will not be vanishing
c     at the end of current time step.
      j=0
      do icluster=1,ncluster
c     Iteration using the number of vanshing member for each cluster
      do i=1,nvan_cluster(icluster) !
         j=j+1
c        cluster slot number for the vanishing member (j)
         imem=nvan_cluster_id(j)
c        If the current member is parent, that is in the '0'-th slot,
c        Find new parent by selecting one that has the next highest
c        dislocation density and that will not be vanished.
         if (imem.eq.0) then
            idmx=-1             ! the slot of next parent
            rhomx=-1.
            do jmem=1,nchildset(icluster)
               kjj=idclusterset(icluster,jmem)
               wgtmp=wgt(kjj)+dwgt_nucl(kjj)+dwgt_grow(kjj)
               rho_tmp=rho_g(kjj) !rho_for(kjj)+rho_deb(kjj)
               if (rho_tmp.gt.rhomx.and.wgtmp.gt.rexdd_tiny) then
                  idmx=jmem
                  rhomx=rho_tmp
               endif
            enddo

c           Change the parent grain.
            itmp=idclusterset(icluster,imem) !
            idclusterset(icluster,0)=idclusterset(icluster,idmx)
            idclusterset(icluster,idmx)=itmp
         endif
      enddo
      enddo

c======================================================================c
c     Re-determination of the number of vanishing members in each cluster
      kount_van_px=0            ! number of vanishing grains overall polycrystal
      kount_van_all=0
      idclusterset_vanish(:)=-1
      idgrain_vanish(:)=-1
      nvan_cluster(:)=0
      nvan_cluster_id(:)=-1
      do 200 icluster=1,ncluster
         nch=nchildset(icluster)
         ivan=.false.

c        Find the cluster stack ID (imemmx) and grain stack ID (idmn) of vanishing grain
         wgtmp_mn=1.0
         idmn=-1                ! grain stack ID
         imemmx=-1              ! cluster stack ID

c        Count the number of vanishing member in the current cluster
         kount_van_cl=0
         do imem=nch,0,-1       !! iteration counting the cluster stack backward.
            kkk=idclusterset(icluster,imem)
            wgtmp=wgt(kkk)+dwgt_nucl(kkk)+dwgt_grow(kkk)
            if (wgtmp.le.rexdd_tiny)then !! potentially vanishing member
               kount_van_all=kount_van_all+1
               kount_van_cl=kount_van_cl+1
               nvan_cluster(icluster)=nvan_cluster(icluster)+1
c              record the stack id of grain in cluster that will vanish
               nvan_cluster_id(kount_van_all)=imem
               ivan=.true.
               if (wgtmp.lt.wgtmp_mn) then
                  wgtmp_mn=wgtmp
                  idmn=kkk
                  imemmx=imem
               endif
            endif
         enddo

         if (ivan) then
            kount_van_px=kount_van_px+1
            idclusterset_vanish(kount_van_px)=icluster
            idgrain_vanish(kount_van_px)=idmn
            idmem_vanish(kount_van_px)=imemmx
         endif
 200  enddo

c======================================================================c
c     Now that none of the vanishing member is parent,
c     Let's rearrange the cluster stack.

      j=0 ! index to each vanishing grain in order
      do icluster=1,ncluster
         do i=1,nvan_cluster(icluster)
c           The number of children can be updated
c           within this loop as the slot for vanishing member
c           pushed to the bottom and is removed.
            nch=nchildset(icluster)
            j=j+1
            imem=nvan_cluster_id(j)
            kkk=idclusterset(icluster,imem)
c           If the vanishing member is not at the bottom
c           switch its slot with the one at the bottom.
            if(imem.ne.nch) then
               kkb=idclusterset(icluster,nch)
               idclusterset(icluster,imem)=
     $              idclusterset(icluster,nch)
            endif
            nchildset(icluster)=nchildset(icluster)-1
            nkidstot=nkidstot-1
         enddo                  ! loop over the members to vanish within the current cluster
c         if (nvan_cluster(icluster).gt.0) then
c            write(*,*)'cluster stack after rearrange/removing'
c            call rexdd_display(1,icluster)
c         endif
      enddo                     ! loop over cluster
c      write(*,*)'here?'
c      read(*,*)

      neli=0
      do ind=1,ngref
         kkk=igref(ind)
c        Add the list of grains to be eliminted from reference grain
         wgtmp=wgt(kkk)+dwgt_nucl(kkk)+dwgt_grow(kkk)
         if (wgtmp.le.rexdd_tiny) then
            neli=neli+1
            ieli(neli)=ind
         endif
      enddo

c======================================================================c
c     Eliminate vanishing grain from the reference grains stack that is
c     used to track the recrystallization

      if (neli.gt.0) then
c        The work-flow is below:
c         1. Loops over initial number of ngref (ngref0)
c         2. See if it is depleted.
c           2-1. If depleted, see if it is at the bottom of 'survived' ones.
c                2-1-1. If it is at the bottom, simply remove it.
c                2-1-2. If is NOT at the bottom, switch with the next one.
c           2-1. If NOT depleted, move on.

         do ind=neli,1,-1
            jnd=ieli(ind)
            kkk=igref(ind)
            wgtmp=wgt(kkk)+dwgt_nucl(kkk)+dwgt_grow(kkk)
            if (ind.ne.ngref) then
               itemp=igref(ngref)
               igref(jnd)=itemp
            endif
            ngref=ngref-1
            do i=1,ngref
               kkk=igref(i)
               wgtmp=wgt(kkk)+dwgt_nucl(kkk)+dwgt_grow(kkk)
            enddo
         enddo
      endif


c     --------------------------------------------------------
c     Let's find grain slots to vacate from bottom towards top
      nvan=0
      idgrain_vanish(:)=-1
      kkk=ngr(1)

c     Stop searching as soon as all vanishing grains are found.
      do while(nvan.lt.kount_van_all.and.kkk.gt.0)

         wgtmp=wgt(kkk)+dwgt_nucl(kkk)+dwgt_grow(kkk)

         if (wgtmp.le.rexdd_tiny) then
            nvan=nvan+1
            idgrain_vanish(nvan)=kkk
         endif

         kkk=kkk-1
      enddo

c     Check whether all expected vanishing grains were found.
      if (nvan.ne.kount_van_all) then
         write(*,*)'nvan:',nvan,
     $        'kount_van_all:',kount_van_all
         write(*,*)'Error: nvan .ne. kount_van_all'
         call xit()
      endif


c      write(*,*)
c      write(*,*)'--- Rearrange grain stack'
c      write(*,*)'grain stack before rearrangement'
c      call rexdd_display(0,-1)

c     --------------------------------------------------------
c     Rearrange the dynamic grain stack.
c
c     kkk : dynamic slot of a grain that will vanish
c     kkb : current bottom-most active dynamic grain slot
c
c     If kkk /= kkb, the grain currently stored at kkb is moved
c     into the vacant slot kkk using <rexdd_inherits>.
c
c     Permanent GID does NOT change during this relocation.
c     --------------------------------------------------------
      kkb=ngr(1)

      do i=1,kount_van_all

c        Dynamic slot to be removed
         kkk=idgrain_vanish(i)

c         write(*,*)kkk,' is to vanish'

         icluster =idclusterset_reverse(kkk)
         iclusterb=idclusterset_reverse(kkb)


c----------------------------------------------------------------------c
c        1. Record the permanent GID of the grain that disappears.
c----------------------------------------------------------------------c

         dead_gid=gid_of_slot(kkk)

c        A valid active slot must have a valid permanent GID.
         if (dead_gid.le.0.or.dead_gid.gt.ngid_total) then
            write(*,*)'** Error: invalid dead GID'
            write(*,*)'slot     :',kkk
            write(*,*)'dead_gid :',dead_gid
            write(*,*)'ngid_total:',ngid_total
            call xit()
         endif

c        Check consistency of the reverse mapping before modification.
         if (slot_of_gid(dead_gid).ne.kkk) then
            write(*,*)'** Error: inconsistent dead-grain mapping'
            write(*,*)'GID:',dead_gid
            write(*,*)'gid_of_slot(',kkk,'):',gid_of_slot(kkk)
            write(*,*)'slot_of_gid(',dead_gid,'):',
     $           slot_of_gid(dead_gid)
            call xit()
         endif

c        --------------------------------------------------------------
c        Check permanent life state before recording death.
c
c        -2 : NOT BORN
c        -1 : ALIVE
c        >=0: already DEAD
c        --------------------------------------------------------------

c        A not-born GID must never appear in an active grain slot.
         if (gid_death_step(dead_gid).eq.-2) then
            write(*,*)'** Error: NOT-BORN GID found in active slot'
            write(*,*)'GID:',dead_gid
            write(*,*)'slot:',kkk
            call xit()
         endif

c        Only an ALIVE grain is allowed to vanish.
         if (gid_death_step(dead_gid).ne.-1) then
            write(*,*)'** Error: grain already marked as dead'
            write(*,*)'GID:',dead_gid
            write(*,*)'death step:',gid_death_step(dead_gid)
            call xit()
         endif

c        Permanent history update:
c        ALIVE (-1) -> DEAD (istep)
         gid_death_step(dead_gid)=istep

c        Zero means that this permanent GID no longer has
c        an active dynamic slot.
         slot_of_gid(dead_gid)=0


c----------------------------------------------------------------------c
c        2. Update the existing reference-grain stack <igref>.
c
c        <igref> still stores DYNAMIC grain slots.
c        If the current bottom slot kkb is going to be moved to kkk,
c        its reference slot must also be changed from kkb to kkk.
c----------------------------------------------------------------------c

c$$$         write(*,*)'kkb:',kkb,'kkk:',kkk
c$$$         write(*,*)'ngref:',ngref
c$$$         write(*,'(2(a5,x))')'ind','idref'
c$$$         do ind=1,ngref
c$$$            write(*,'(2(i5,x))')ind,igref(ind)
c$$$         enddo

         do ind=1,ngref

            if (kkb.eq.igref(ind)) then

c$$$               write(*,*)'kkb is in igref'
c$$$               write(*,*)kkb,ind,igref(ind)
c$$$
c$$$               write(*,*)'before changing igref'
c$$$               do jnd=1,ngref
c$$$                  write(*,'(2(i5,x))')jnd,igref(jnd)
c$$$               enddo

               igref(ind)=kkk

c$$$               write(*,*)'after changing igref'
c$$$               do jnd=1,ngref
c$$$                  write(*,'(2(i5,x))')jnd,igref(jnd)
c$$$               enddo

c              Only one corresponding reference slot is expected.
               exit
            endif

         enddo


c----------------------------------------------------------------------c
c        3. If the vanished grain is not already at the bottom,
c           move the bottom grain into the vacant dynamic slot.
c----------------------------------------------------------------------c

         if (kkk.ne.kkb) then

c           Permanent GID of the grain currently stored at the bottom.
            move_gid=gid_of_slot(kkb)

            if (move_gid.le.0.or.move_gid.gt.ngid_total) then
               write(*,*)'** Error: invalid moving GID'
               write(*,*)'slot     :',kkb
               write(*,*)'move_gid :',move_gid
               write(*,*)'ngid_total:',ngid_total
               call xit()
            endif

c           Check reverse pointer before relocation.
            if (slot_of_gid(move_gid).ne.kkb) then
               write(*,*)'** Error: inconsistent moving-grain mapping'
               write(*,*)'GID:',move_gid
               write(*,*)'gid_of_slot(',kkb,'):',
     $              gid_of_slot(kkb)
               write(*,*)'slot_of_gid(',move_gid,'):',
     $              slot_of_gid(move_gid)
               call xit()
            endif

c           A grain that is moved within the dynamic stack
c           must still be alive.
            if (gid_death_step(move_gid).ne.-1) then
               write(*,*)'** Error: moving grain is not alive'
               write(*,*)'GID:',move_gid
               write(*,*)'death state:',
     $              gid_death_step(move_gid)
               call xit()
            endif


c           Move all physical/state variables from slot kkb to slot kkk.
            call rexdd_inherits(
c                intent(in)::
     $           kkk,kkb,ngrmx,nsysmx,nclustermx,nchildmx,nchildset

c                intent(inout)::
c                ** REXDD related quantities
     $           ,idclusterset,idclusterset_reverse
     $           ,pnucl,dwgt_nucl,dwgt_grow,velgrow,radius

c                ** grain data & properties
     $           ,ag,wgt,fijgr,axisgr,crss,sg,dsg,sg_old,sgtry,dsgtry
     $           ,dg0,dg_trans,gtotgr,mgtg

c                ** dislocation densities
     $           ,rho_s_for,rho_s_deb,rho_g
     $           )


c----------------------------------------------------------------------c
c           4. Update the permanent-ID <-> dynamic-slot mapping.
c
c           The grain itself keeps <move_gid>.
c           Only its dynamic slot changes from kkb to kkk.
c----------------------------------------------------------------------c

            gid_of_slot(kkk)=move_gid
            slot_of_gid(move_gid)=kkk

         endif


c----------------------------------------------------------------------c
c        5. Remove the old bottom dynamic slot.
c
c        After this point, kkb is outside the active grain stack.
c----------------------------------------------------------------------c

         gid_of_slot(kkb)=0

         kkb=kkb-1
         ngr(1)=ngr(1)-1

      enddo

c$$$      write(*,*)'grain stack after rearrangement'
c$$$      call rexdd_display(0,-1)
c$$$      write(*,*)'cluster stack after rearrangement'
c$$$      do icluster=1,ncluster
c$$$         call rexdd_display(1,icluster)
c$$$      enddo
c$$$      write(*,*)'reference grain list'
c$$$      write(*,*)'ngref:',ngref
c$$$      do i=1,ngref
c$$$         write(*,*)igref(i)
c$$$      enddo
c      read(*,*)
c$$$
c$$$      if (kount_van_all.gt.0) then
c$$$         write(*,*)'nvan:',nvan,'kount_van_all:',kount_van_all
c$$$         write(*,*)'* stopped for checking the rearranged cluster stack'
c$$$c         call xit()
c$$$      endif
c$$$c      call xit()
c$$$c      write(*,*)'** end of subr. rexdd_cluster_vanish'



      end subroutine rexdd_vanish
