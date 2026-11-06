c     subroutine rexdd_maintain tries to maintain a certain number
c     of cluster within each cluster to maintain consistency across
c     all clusters considered in the simulations.

c     It finds clusters that have either too many grains or too few
c     grains and reassign grains such that each cluster has nearly the
c     same number of members throughout the calculation.
      subroutine rexdd_maintain
     $     (
c         intent(in)::
c         * integer
     $     ngrmx,nchildmx,nclustermx,nchildlimit
c         * real
     $     ,dwgt_nucl
c         intent(inout)::
c         * integer
     $     ,ncluster,nchildset,idclusterset,idclusterset_reverse
c         intent(out)::
c         * integer
     $     )
      implicit none
c     intent(in)::
      integer,intent(in)::ngrmx,nchildmx,nclustermx,nchildlimit
      real*8, intent(in)::dwgt_nucl(ngrmx)

c     intent(inout)::
      integer,intent(inout)::ncluster,nchildset(nclustermx)
     $     ,idclusterset(nclustermx,0:nchildmx)
     $     ,idclusterset_reverse(ngrmx)
c     intent(out)::
c     locals
      integer icluster,ncluster0,nmem,nextra,ndextra,nextraset(ngrmx),
     $     iextra,navailset(ngrmx),navail,iavail,nfill,nncl,i,j,k,itop,
     $     ibot,ichosen,kkk,itrial,nrem,ieach,imem,neach,ngr_tot
      real*8 rand,avgmem_init,avgmem_fin
c     arrays/variables used to find full clusters and dismiss some
c     to minimize 'unfilled' clusters.
      integer nempty,icl_select,mxiter,nempty_clist(ngrmx),ndmx,
     $     ncluster_init,l,jcl_select,jmem,jkk,ikk
      parameter(mxiter=10000,ndmx=100000)

c--------------------------------------------------------------------
c     - Determine how many clusters have abundant grains; and how many
c         grains are abundant (nextra, nextraset)?
c     - Determine how many and which clusters have rooms to accommodate
c       new grains in their pool (navail, navailset)?
c     (**) Only when there is an extra member (nmem>nchildlimit)
c     - Randomly choose grains to be relocated to maintain mem=nchildlimit
c--------------------------------------------------------------------
 1    format(a15,i5)
c     Loop over clusters to find
c       1) clusters that have 'extra members'
c       2) clusters that have unoccupied slots that can accomodate more members
      nextra=0                  ! No of grains to be relocated
      navail=0                  ! No of available slots
      nextraset(:)=-1           ! List of grain IDs that are to be relocated
      navailset(:)=-1           ! List of available slots' cluster ID
      ncluster0=ncluster        ! init. number of clusters before adding any addition

      do icluster=1,ncluster0
c        nmem: the number of total members (including both parent and children)
         nmem=nchildset(icluster)+1
         if (nmem.lt.nchildlimit) then
c           Available slots
            nrem=nchildlimit-nmem
            do i=1,nrem
               navail=navail+1
c              cluster ID for each spot.
               navailset(navail)=icluster
            enddo
         elseif (nmem.gt.nchildlimit) then
c           Randomly choose nrem grains among nmem grains.
            nrem=nmem-nchildlimit
c           itop: top of the cluster stack - it is always 0 corresponding to parent.
c           ibot: bottom of the cluster stack
            itop=0
            ibot=nchildset(icluster) ! It is updated.
            do i=1,nrem
c              between 0 and nch
               itrial=0
 10            call rand_i0i1(itop,ibot,ichosen)
               kkk=idclusterset(icluster,ichosen)
c              Pass if it is newly nucleated grain, and
c              choose another one that is not newly created ones.
               if (dwgt_nucl(kkk).gt.0) then
                  itrial=itrial+1
                  if (itrial.gt.1000) then
                     write(*,*)'** Error, too many trials to avoid '
                     write(*,*)' new grains'
                     call rexdd_display(1,icluster)
                     call xit()
                  endif
                  goto 10
               endif
               nextra=nextra+1
               nextraset(nextra)=kkk
               if (ichosen.ne.ibot) then
c                 if the chosen grain is NOT at the 'current' bottom in
c                 that cluster, switch to-be-emptied slot with
c                 the bottom member not to leave any empty spot
c                 in the middle of cluster stack.
                  idclusterset(icluster,ichosen)=
     $                 idclusterset(icluster,ibot)
               endif
c              Caution requires for the next statement:
c              the extra members are not yet reallocated to
c              either to available slots or to new clusters yet.
               nchildset(icluster)=nchildset(icluster)-1
               ibot=ibot-1          ! Shifts up (decrease) the bottom ID
            enddo
         endif
      enddo                     ! loop over ncluster0

c--------------------------------------------------------------------
c     If there's no extra grains to be relocated, exit the subroutine
c--------------------------------------------------------------------
      if (nextra.eq.0) goto 200
c--------------------------------------------------------------------
c     0. If there are extra grains to be relocated, shuffle the orders
c        of (NEXTRASET) and (NAVAILSET)
c     1. Firstly, relocate them to the 'available' slots in the cluster.
c     2. Secondly, relocate the remaining grains to 'new' clusters.
c       2-1. Determine how many clusters should be created.
c       2-2. Determine how many remaining grains should be relocated
c            to each cluster.
c       2-3. Relocate the remaining ones.
c--------------------------------------------------------------------

c     0. Shuffle
c     0-1. Shuffle nextraset
      if (nextra.gt.1) call shuffle(nextraset(1:nextra),nextra)
c     0-2. Shuffle navailset
      if (navail.gt.1) call shuffle(navailset(1:navail),navail)

c     1. First relocate the extra grains to available slots
      if (nextra.gt.0.and.navail.ge.nextra) then
c        The number of navail greater than or equal to nextra
         nfill=nextra
         nrem=0
      elseif (nextra.gt.0.and.navail.lt.nextra) then
c        The number of navail is less than nextra
         nfill=navail
         nrem=nextra-navail
      endif

c     These indices (iavail and iextra) used for the loop over nfill
c     runs from bottom (navail and nextra) towards top (0)
      iavail=navail
      iextra=nextra
      if (nfill.gt.0) then
         do i=1,nfill
            icluster=navailset(iavail)
            kkk=nextraset(iextra)
            imem=nchildset(icluster)+1
            idclusterset(icluster,imem)=kkk
            idclusterset_reverse(kkk)=icluster
            nchildset(icluster)=nchildset(icluster)+1
            nextra=nextra-1
            navail=navail-1
            iavail=iavail-1
            iextra=iextra-1
         enddo
c         write(*,*)'check if filling grain out of nfill'
c         write(*,*)'is done correctly'
c         call xit()
      endif

c     2. Relocate the remaining grains with creating new clusters.
c     2-1. Determine how many new clusters are required? (nncl)
      neach=-1
      if (nrem.gt.0) then
         nncl=ceiling(nrem*1./nchildlimit*1.)
c     2-2. Determine how many remaining grains should be relocated
c          to each cluster (neach)
         neach=ceiling(nrem*1./(nncl*1.))
      endif
c     2-3. Relocate the remaining ones (nrem)
      do 100 while (nrem.gt.0)
c$$$         write(*,*)'relocating remaining grains to new clusters'
c$$$         write(*,1)'icluster:',icluster
         icluster=ncluster+1
         ncluster=ncluster+1
         nchildset(icluster)=-1
c        for the current icluster, fill neach members
c        as long as nrem stays higher than 0.
c        We fill from top (imem=0) towards bottom.
         imem=0
         do ieach=1,neach ! maximum neach, may quit if nrem becomes 0
            kkk=nextraset(iextra)
            idclusterset(icluster,imem)=kkk
            idclusterset_reverse(kkk)=icluster
            nchildset(icluster)=nchildset(icluster)+1
            nextra=nextra-1
            navail=navail-1
            iavail=iavail-1
            iextra=iextra-1
            nrem=nrem-1
            if (nrem.eq.0) then
c               write(*,*)'Completed the relocation'
               goto 100
            endif
            imem=imem+1
         enddo
 100  continue

c     sanity check
      if (nextra.gt.0)then
         write(*,1)'nextra after done:',nextra
         call xit()
      endif
c     Clean idclusterset by filling -1 for unfilled slots.
      do icluster=1,ncluster
         nmem=nchildset(icluster) ! 0,....,nmem
         idclusterset(icluster,nmem+1:)=-1
      enddo

c-----------------------------------------------------------------------
c     --- avg mem
 200  avgmem_init=0.
      do icluster=1,ncluster
         avgmem_init=avgmem_init+nchildset(icluster)+1
      enddo
      avgmem_init=avgmem_init/(ncluster*1.)

      nempty=0
      do icluster=1,ncluster
         nempty=nempty+nchildlimit-(nchildset(icluster)+1)
      enddo
      if (nempty.ge.nchildmx) then
         write(*,1)'old nempty:',nempty
         call random_number(rand)
         icl_select = 1 + floor(ncluster * rand)
         write(*,1)'icl_select:',icl_select
c        Now, renew nempty, excluding the ones in the chosen cluster.
c         nempty=nempty-(nchildlimit-(nchildset(icl_select)+1))
c         write(*,1)'new nempty:',nempty
         call rexdd_get_empty_clist(nchildmx,nclustermx,nchildlimit,
     $        ncluster,nchildset,icl_select,
     $        ndmx,
     $        nempty,nempty_clist)
      endif

c     initial number of clusters before removing any cluster.
      ncluster_init=ncluster

c     Relocate members in the chosen clsuter only if this new <nempty>
c     is greater than or equal to nchildmx
      if (nempty.ge.nchildmx) then
         k=1 ! k is the number of clusters that are going to be dismissed.
         do 300 while(nempty.ge.nchildmx.and.k.lt.mxiter)
            write(*,*)'** STEP 3; reassign grains in the chosen cluster'

c           i: refers to the cluster to be emptied.
c           j: refers to the cluster to be filled.
c           l: refers to the number of such 'relocations' from one cluster to another
c           members in icl_select
            write(*,*)'members in icl_select'
            write(*,*)'nchildset(jcl_select):',nchildset(icl_select)
            do i=0,nchildset(icl_select)
               write(*,*)i,idclusterset(icl_select,i)
            enddo
            imem=nchildset(icl_select)+1

            l=1
            write(*,*)'----'
            write(*,'(6a7)')'k','l','imem','icl','jcl','nempty'
            do 250 while(nempty.gt.0.and.imem.gt.0.and.l.lt.mxiter)
c              cluster ID to be emptied is   <ICL_SELECT>.
c              cluster ID to be filled in is <JCL_SELECT>.
               jcl_select=nempty_clist(l)
               write(*,'(6i7)')k,l,imem,icl_select,jcl_select,nempty
c              current number of members in that cluster
               jmem=nchildset(jcl_select)+1 ! shift one to the available slot

c              Reassign the grain in the old cluster <ICL> to new <JCL> cluster.
               ikk=idclusterset(icl_select,imem-1) ! imem in (0,...,nchildset)
               idclusterset(jcl_select,jmem)=ikk
               idclusterset(icl_select,imem-1)=-1  ! ensure it now is an empty slot
c              The number of members in <ICL> and <JCL> are updated.
               nchildset(jcl_select)=nchildset(jcl_select)+1
               nchildset(icl_select)=nchildset(icl_select)-1
c              The grain <IKK> is now included in the new cluster, thus
c              idclusterset_reverse is updated for this grain <IKK>.
               idclusterset_reverse(ikk)=jcl_select
               imem=imem-1
               nempty=nempty-1
               l=l+1
 250        continue
            write(*,*)'----'
            write(*,*)'** STEP 4; reorder cluster'
            if (icl_select.eq.ncluster) then
               write(*,*)'it happended to be the last cluster'
            else
               jcl_select=ncluster
               do jmem=0,nchildset(jcl_select)
                  jkk=idclusterset(jcl_select,jmem)
                  write(*,*)jmem,jkk
                  idclusterset(icl_select,jmem)=jkk
                  nchildset(icl_select)=nchildset(icl_select)+1
                  idclusterset_reverse(jkk)=icl_select
               enddo
            endif
            ncluster=ncluster-1

            write(*,*)'** STEP 5; randomly pick a cluster'
            call random_number(rand)
            icl_select = 1 + floor(ncluster * rand)
c$$$            write(*,1)'icl_select:',icl_select,
c$$$     $           'its nch:',nchildset(icl_select)

            write(*,*)'** STEP 6; renew empty slots'
            call rexdd_get_empty_clist(nchildmx,nclustermx,nchildlimit,
     $           ncluster,nchildset,icl_select,
     $           ndmx,
     $           nempty,nempty_clist)
            write(*,1)'new nempty:',nempty
            k=k+1
 300     continue

c     --- avg mem
         avgmem_fin=0.
         do icluster=1,ncluster
            avgmem_fin=avgmem_fin+nchildset(icluster)+1
         enddo
         avgmem_fin=avgmem_fin/(ncluster*1.)

         write(*,'(4a10)')'ncl_init','ncl_fin','avgmemi','avgmemf'
         write(*,'(2i10,2f10.1)')
     $        ncluster_init,ncluster,avgmem_init,avgmem_fin
         if (avgmem_init.gt.avgmem_fin) then
            write(*,*)'** Error: unexpected outcome in <rexdd_maintain>'
            write(*,*)'average number of members in a cluster reduced'
         endif
         if (k.eq.mxiter) then
            write(*,*)'** Error: reached mix iter in <rexdd_main>'
            call xit()
         endif
      endif

 1000 return
      end subroutine rexdd_maintain
c-----------------------------------------------------------------------
c     This subroutine used only in <rexdd_maintain>.
c     This subroutine perform two tasks:
c     1) It creates the new inter array that collects the each empty slot.
c        It spcifically records only the cluster id of such, and they are
c        allowed to be duplicated.
c     2) In doing 1), it excludes a single cluster (icl_exclude).
c
c     Consequently, the resulting integer array <nempty_clist> collects
c     the cluster ID of such empty slots (again, they can be duplicated),
c     and the number of empty slots is recorded to <nempty>.
c

      subroutine rexdd_get_empty_clist(
c          intent(in)
     $     nchildmx,nclustermx,nchildlimit,
     $     ncluster,nchildset,icl_exclude,ndmx,
c          intent(out)
     $     nempty,nempty_clist
     $     )
      implicit none
      integer, intent(in):: nchildmx,nclustermx,nchildlimit,ncluster,
     $     nchildset(nclustermx),icl_exclude,ndmx
      integer, intent(out):: nempty,nempty_clist(ndmx)
c     local
      integer i,icluster,nmem
      nempty=0
      nempty_clist(:)=-1
      do 10 icluster=1,ncluster
         if (icluster.eq.icl_exclude) goto 10
         nmem=nchildlimit-(nchildset(icluster)+1)
         if (nmem.eq.0) goto 10

         do i=1,nmem
            nempty=nempty+1
            if (nempty.gt.ndmx) then
               write(*,*)'** Err: increase <ndmx>'
               call xit()
            endif
            nempty_clist(nempty)=icluster
         enddo
 10   continue

c     Return shuffled ones?
      call shuffle(nempty_clist(1:nempty),nempty)

      return
      end subroutine rexdd_get_empty_clist
