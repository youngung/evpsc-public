      subroutine rexdd_disp_grs(
     $     nphmx,ngrmx,nclustermx,nchildmx,
     $     ngr,nchildset,idclusterset,idclusterset_reverse,
     $     gid_of_slot,
     $     rexdd_tiny,rho_g,wgt,dwgt_nucl,dwgt_grow,radius,
     $     crss1,pnucl
     $     )
      implicit none
      integer, intent(in):: nphmx,ngrmx,ngr(0:nphmx),
     $     nclustermx,nchildmx,nchildset(nclustermx),
     $     idclusterset(nclustermx,0:nchildmx),
     $     idclusterset_reverse(ngrmx),gid_of_slot(ngrmx)
      real*8, intent(in):: wgt(ngrmx),dwgt_nucl(ngrmx),dwgt_grow(ngrmx),
     $     rho_g(ngrmx),radius(ngrmx),rexdd_tiny,crss1(ngrmx),
     $     pnucl(ngrmx)
c     locals
      integer kkk,icluster,nch,imem
      real*8 wgtnew

      write(*,*)'-- Grain stack display'
      write(*,'(5a5,8a11)')'slot','GID','iclu','imem','nch',
     $     'rho','wgtold','dw_nuc','dw_gro','wgtnew','rad[um]',
     $     'crss1','Pnucl'
      do kkk=1,ngr(1)
         wgtnew=wgt(kkk)+dwgt_nucl(kkk)+dwgt_grow(kkk)
         icluster=idclusterset_reverse(kkk)
c     find the member?
         nch=nchildset(icluster)
         do imem=0,nch
            if (idclusterset(icluster,imem).eq.kkk) goto 10
         enddo
 10      write(*,'(5i5,8e11.3)',advance='no')
     $     kkk,gid_of_slot(kkk),icluster,imem,nch,
     $        rho_g(kkk),wgt(kkk),dwgt_nucl(kkk),dwgt_grow(kkk),wgtnew,
     $        radius(kkk)*1e6,crss1(kkk),pnucl(kkk)
         if (wgtnew.lt.0) write(*,'(a5)',advance='no')'neg'
         if (dwgt_nucl(kkk).gt.0)write(*,'(a5)',advance='no')'nuc+'
         if (dwgt_nucl(kkk).lt.0)write(*,'(a5)',advance='no')'nuc-'
         if (wgtnew.lt.rexdd_tiny) write(*,'(a5)',advance='no')'van*'
         write(*,*)
         if (nch.lt.0) write(*,*)'*Error, negative nch'
      enddo

      return
      end subroutine rexdd_disp_grs
