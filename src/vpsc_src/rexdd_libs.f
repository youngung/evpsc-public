c     shuffle the given 'integer' array
c     according to Fisher-Yates algorithm
      subroutine shuffle(arr,n)
      implicit none
      integer, intent(in):: n
      integer, intent(inout) :: arr(n)
      integer :: i, j, temp
      real :: r
c      call random_seed()
      do i = n, 2, -1
         call random_number(r)
         j = 1 + int(r * i)     ! Random index from 1 to i
         temp = arr(i)
         arr(i) = arr(j)
         arr(j) = temp
      end do
      end subroutine


c---
      subroutine rand_i0i1(i0,i1,ichosen)
      implicit none
      integer, intent(in):: i0,i1
      integer, intent(out):: ichosen
c
      real*8 rand
c      call random_seed()
      call random_number(rand)
      ichosen=i0+int((i1-i0+1)*rand)
      return
      end subroutine rand_i0i1
c======================================================================c
c     Arguments
c     ---------
c     mdim:   dimension of Schmid tensor in b-basis (5 or 6?)
c     nsysmx: maximum dimension of Schca
c     nsyst: actual number of slip systems
c     schca: schmid tensor referred in the crystal axes (of phase)
c     schsa: schmid tensor referred in the sample axes, transformed by rotation matrix R.
c     R    : the transformation matrix that transforms from the crystal axes to
c            the sample axes.
      subroutine rexdd_schca2sa(mdim,nsysmx,nsyst,schca,r,schsa)
      implicit none
      integer, intent(in):: mdim,nsysmx,nsyst
      real*8, intent(in)::  schca(mdim,nsysmx), r(3,3)
      real*8, intent(out):: schsa(mdim,nsysmx)
c     locals
      integer is,i,j,k,l
      real*8 auxmm(mdim,mdim),aux3333(3,3,3,3),
     $     schbb(mdim),schca_cart(3,3),schsa_cart(3,3)

      do is=1,nsyst
         do i=1,mdim
            schbb(i)=schca(i,is)
         enddo
c        from b-basis to Cartesian
         call chg_basis(schbb,schca_cart,auxmm,aux3333,1,mdim)
         schsa_cart(:,:)=0.
         do i=1,3
         do j=1,3
            do k=1,3
            do l=1,3
               schsa_cart(i,j)=schsa_cart(i,j)
     $              +schca_cart(k,l)*r(i,k)*r(j,l)
            enddo
            enddo
         enddo
         enddo

         call chg_basis(schbb,schsa_cart,auxmm,aux3333,2,mdim)

         do i=1,mdim
            schsa(i,is)=schbb(i)
         enddo

      enddo

      return
      end subroutine rexdd_schca2sa


c======================================================================c
c     Given the parent's tensor with b-basis components, referred in
c     the sa, transform it to its child's crystal, then to sample axes.

c     agpar : transforms from sample axes to crystal axes of parent
c     agchi : transforms from sample axes to crystal axes of child

c     iopt0: inherits 2nd rank tensor
c     iopt1: inherits 4th rank tensor

c     sa_ch<-ca_pa<-sa_pa
c     R_ij = agchi_ki agpar_kj : sa_ch <- sa_pa
      subroutine rexdd_transf(tmpar_bb2,tmpar_bb4,tmpchi_bb2,
     $     tmpchi_bb4,agpar,agchi,iopt,mdim)
      implicit none
      integer, intent(in):: iopt, mdim
      real*8,intent(in):: tmpar_bb2(mdim),tmpar_bb4(mdim,mdim),
     $     agpar(3,3),agchi(3,3)
      real*8, intent(out):: tmpchi_bb2(mdim),tmpchi_bb4(mdim,mdim)
c     locals
      real*8 aux33(3,3),bux33(3,3),aux3333(3,3,3,3),bux3333(3,3,3,3),
     $     auxm(mdim),auxmm(mdim,mdim),R(3,3)
      integer i,j,k,l,m,n,o,p

      R(:,:)=0d0
      do i=1,3
      do j=1,3
         do k=1,3
            R(i,j)=R(i,j)+agchi(k,i)*agpar(k,j)
         enddo
      enddo
      enddo

      if (iopt.eq.0) then !! inherits (mdim) 2nd rank tensor
c        converts it from b-basis to cartesian components (aux33)
         call chg_basis(tmpar_bb2,aux33,auxmm,aux3333,1,mdim)
c        transform it (aux33) to crystal axes of parent (bux33),
c        and transform it to sample axes of child
c        b_ij = R_ik R_jl a_kl
         bux33(:,:)=0d0
         do i=1,3
         do j=1,3
c
         do k=1,3
         do l=1,3
            bux33(i,j)=bux33(i,j)+R(i,k)*R(j,l)*aux33(k,l)
         enddo
         enddo
         enddo
         enddo
c        converts it (bux33) back to b-basis components
         call chg_basis(tmpchi_bb2,bux33,auxmm,aux3333,2,mdim)
      elseif(iopt.eq.1) then !! inherits (mdim,mdim) 4th rank tensor
c        converts it to cartesian coorindates
         call chg_basis(auxm,aux33,tmpar_bb4,aux3333,3,mdim)
c        transform it (aux3333) to crystal axes of parent (bux3333).
         bux3333(:,:,:,:)=0d0
         do i=1,3
         do j=1,3
         do k=1,3
         do l=1,3
c
         do m=1,3
         do n=1,3
         do o=1,3
         do p=1,3
            bux3333(i,j,k,l)=bux3333(i,j,k,l)
     $           +R(i,m)*R(j,n)*R(k,o)*R(l,p)*aux3333(m,n,o,p)
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo
c        converts it back to b-basis components
         call chg_basis(auxm,aux33,tmpchi_bb4,bux3333,4,mdim)
      endif
      return
      end subroutine


c======================================================================c
c     displays grain stacks (iopt=0) or cluster stack (iopt=1)
c     Also, option 2 display the reference grain stack used to track down
c     the recrystallized volume fraction
      subroutine rexdd_display(iopt,icl2show)
c$$$      use vpsc8dim
c$$$      use disden
c$$$      use rextaldd_cluster
      integer,intent(in):: iopt,icl2show
c     locals
      integer kkk
      if (iopt.eq.-1) then
         write(*,'(4a5,8a11)')'kkk','iclu','imem','nch','rho',
     $        'wgtold','dw_nuc','dw_gro','wgtnew','wdot','rad[um]',
     $        'crss(1,kkk)'
         return
      endif
c$$$
c$$$      if (iopt.eq.0) then !! sorting according to grain id
c$$$         write(*,*)'-- Grain stack display'
c$$$         write(*,'(4a5,8a11)')'kkk','iclu','imem','nch','rho',
c$$$     $        'wgtold','dw_nuc','dw_gro','wgtnew','wdot','rad[um]',
c$$$     $        'crss(1,kkk)'
c$$$         do kkk=1,ngr(1)
c$$$            wgtnew=wgt(kkk)+dwgt_nucl(kkk)+dwgt_grow(kkk)
c$$$            icluster=idclusterset_reverse(kkk)
c$$$c           find the member?
c$$$            nch=nchildset(icluster)
c$$$            do imem=0,nch
c$$$               if (idclusterset(icluster,imem).eq.kkk) goto 10
c$$$            enddo
c$$$ 10         write(*,'(4i5,5e11.3,3e11.3)',advance='no')kkk,icluster,
c$$$     $           imem,nch,rho_for(kkk)+rho_deb(kkk),wgt(kkk),
c$$$     $           dwgt_nucl(kkk),dwgt_grow(kkk),wgtnew,wgtdot(kkk),
c$$$     $           radius(kkk)*1e6,crss(1,kkk)
c$$$            if (wgtnew.lt.0) write(*,'(a5)',advance='no')'neg'
c$$$            if (dwgt_nucl(kkk).gt.0)write(*,'(a5)',advance='no')'nuc+'
c$$$            if (dwgt_nucl(kkk).lt.0)write(*,'(a5)',advance='no')'nuc-'
c$$$            if (wgtnew.lt.rexdd_tiny) write(*,'(a5)',advance='no')'van*'
c$$$            write(*,*)
c$$$            if (nch.lt.0) write(*,*)'*Error, negative nch'
c$$$         enddo
c$$$         write(*,*)'(0) nkidstot:',nkidstot,'ngr(1):',ngr(1)
c$$$         write(*,*)'avg nkids per parent:',nkidstot/(ncluster*1.)
c$$$      elseif(iopt.eq.1) then !! sorting aoccrding to cluster id
c$$$         write(*,*)'-- cluster stack display'
c$$$         write(*,'(6a5,5a11,2x,a,2x,6a11)')
c$$$     $        'iclu','icrev','imem','nch','kpa','kkk',
c$$$     $        'rho','wgtold','dwgt_nuc','dwgt_grow','wgtnew','|',
c$$$     $        'wgtold','dwgt_nuc','dwgt_grow','wgtnew','rad[um]',
c$$$     $        'crss(1,kkk)'
         if (icl2show.lt.0)  then
c$$$            do icluster=1,ncluster
c$$$               call rexdd_display_cluster(icluster)
c$$$            enddo
c         else
c$$$            call rexdd_display_cluster(icl2show)
         endif
c$$$      elseif(iopt.eq.2) then !! displays reference grain stack
c$$$         write(*,*)'-- Reference grain stack display'
c$$$         write(*,*)'ngref:',ngref
c$$$         call rexdd_display_grain_colhead
c$$$         do ind=1,ngref
c$$$            kkk=igref(ind)
c$$$            wgtnew=wgt(kkk)+dwgt_nucl(kkk)+dwgt_grow(kkk)
c$$$            icluster=idclusterset_reverse(kkk)
c$$$c           find the member?
c$$$            nch=nchildset(icluster)
c$$$            do imem=0,nch
c$$$               if (idclusterset(icluster,imem).eq.kkk) goto 20
c$$$            enddo
c$$$ 20         write(*,'(4i5,5e11.3,3e11.3)',advance='no')kkk,icluster,
c$$$     $           imem,nch,rho_for(kkk)+rho_deb(kkk),wgt(kkk),
c$$$     $           dwgt_nucl(kkk),dwgt_grow(kkk),wgtnew,wgtdot(kkk),
c$$$     $           radius(kkk)*1e6,crss(1,kkk)
c$$$            if (wgtnew.lt.0) write(*,'(a5)',advance='no')'neg'
c$$$            if (dwgt_nucl(kkk).gt.0)write(*,'(a5)',advance='no')'nuc+'
c$$$            if (dwgt_nucl(kkk).lt.0)write(*,'(a5)',advance='no')'nuc-'
c$$$            if (wgtnew.lt.rexdd_tiny) write(*,'(a5)',advance='no')'van*'
c$$$            write(*,*)
c$$$         enddo
c$$$      elseif(iopt.eq.3) then   ! a single grain is shown without column heads
c$$$         kkk=icl2show
c$$$         wgtnew=wgt(kkk)+dwgt_nucl(kkk)+dwgt_grow(kkk)
c$$$         icluster=idclusterset_reverse(kkk)
c$$$c        find the member?
c$$$         do imem=0,nchildset(icluster)
c$$$            if (idclusterset(icluster,imem).eq.kkk) goto 30
c$$$         enddo
c$$$ 30      write(*,'(4i5,5e11.3,3e11.3)',advance='no')kkk,icluster,
c$$$     $        imem,nchildset(icluster),rho_for(kkk)+rho_deb(kkk),
c$$$     $        wgt(kkk),dwgt_nucl(kkk),dwgt_grow(kkk),wgtnew,wgtdot(kkk),
c$$$     $        radius(kkk)*1e6,crss(1,kkk)
c$$$         if (wgtnew.lt.0) write(*,'(a5)',advance='no')'neg'
c$$$         if (dwgt_nucl(kkk).gt.0)write(*,'(a5)',advance='no')'nuc+'
c$$$         if (dwgt_nucl(kkk).lt.0)write(*,'(a5)',advance='no')'nuc-'
c$$$         if (wgtnew.lt.rexdd_tiny) write(*,'(a5)',advance='no')'van*'
c$$$         write(*,*)
c$$$         if (nch.lt.0) write(*,*)'*Error, negative nch'
c$$$      endif
      end subroutine rexdd_display
c$$$c======================================================================c
      subroutine rexdd_display_grain_colhead
      write(*,'(4a5,5a11,2a11)')'kkk','iclu','imem','nch','rho',
     $     'wgtold','dw_nuc','dw_gro','wgtnew','wdot','rad[um]'
      return
      end subroutine rexdd_display_grain_colhead
c$$$c======================================================================c
c$$$      subroutine rexdd_display_cluster(icluster)
c$$$      use vpsc8dim
c$$$      use disden
c$$$      use rextaldd_cluster
c$$$      integer,intent(in):: icluster
c$$$      kpa=idclusterset(icluster,0)
c$$$      nch=nchildset(icluster)
c$$$      nchtot=nchtot+nch
c$$$      wgt_newc=0.
c$$$      wgt_oldc=0.
c$$$      write(*,*)'--- no of children:',nch
c$$$      wgt_clu=0.
c$$$      do imem=0,nch
c$$$         kkk=idclusterset(icluster,imem)
c$$$         wgt_clu=wgt_clu+wgt(kkk)
c$$$      enddo
c$$$
c$$$      do imem=0,nch
c$$$         kkk=idclusterset(icluster,imem)
c$$$         iclrev=idclusterset_reverse(kkk)
c$$$         wgtnew=wgt(kkk)+dwgt_nucl(kkk)+dwgt_grow(kkk)
c$$$         wgt_oldc=wgt_oldc+wgt(kkk)
c$$$         wgt_newc=wgt_newc+wgtnew
c$$$         write(*,'(6i5,5e11.3,2x,a,2x,5f11.3,e11.3)',advance='no')
c$$$     $        icluster,iclrev,imem,
c$$$     $        nch,kpa,kkk,rho_for(kkk)+rho_deb(kkk),
c$$$     $        wgt(kkk),dwgt_nucl(kkk),dwgt_grow(kkk),wgtnew,
c$$$     $        '|',
c$$$     $        wgt(kkk)/wgt_clu,dwgt_nucl(kkk)/wgt_clu,
c$$$     $        dwgt_grow(kkk)/wgt_clu,wgtnew/wgt_clu,
c$$$     $        radius(kkk)*1e6,crss(1,kkk)
c$$$         if (wgtnew.lt.0) write(*,'(a5)',advance='no')'neg'
c$$$         if (dwgt_nucl(kkk).gt.0)write(*,'(a5)',advance='no')'nuc+'
c$$$         if (dwgt_nucl(kkk).lt.0)write(*,'(a5)',advance='no')'nuc-'
c$$$         if (wgtnew.lt.rexdd_tiny) write(*,'(a5)',advance='no')'van*'
c$$$         write(*,*)
c$$$         if (kkk.gt.ngr(1)) write(*,*)'Error! kkk exceeded ngr(1)'
c$$$         if (kpa.gt.ngr(1)) write(*,*)'Error! kpa exceeded ngr(1)'
c$$$      enddo                     ! Loop over member
c$$$      if (nch.lt.0) then
c$$$         write(*,*)'*Error, negative nch'
c$$$         call xit()
c$$$      endif
c$$$      if (abs(wgt_oldc-wgt_newc).gt.1e-10) then
c$$$         write(*,*)'** warn: wgt within cluster is not preserved'
c$$$      endif
c$$$      write(*,'(2(a10,e11.3))')'wgt_oldc:',wgt_oldc,'wgt_newc:',wgt_newc
c$$$      write(*,*)'--'
c$$$      end subroutine rexdd_display_cluster
