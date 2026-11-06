c     ----------------------------------------------------------------------------
c     subroutine calculates the contributions from
c     elasticity, visco-plasticity, dg_trans (e.g., thermal expansion/contraction)
c     ----------------------------------------------------------------------------
      subroutine calc_contribution(
c          intent(in)
     $     uw,ife,ngrmx,nphmx,nph,ngr,dgp,dge,dg_trans,wgt,
c          intent(out)
     $     contrib)
      implicit none
      integer, intent(in)::ngrmx,nphmx,nph,ngr(0:nphmx),uw,ife
      real*8, intent(in)::dgp(5,ngrmx),dge(6,ngrmx),dg_trans(6,ngrmx),
     $     wgt(ngrmx)
      real*8, intent(out)::contrib(3)
c     locals
      integer i,kkk,iph
      real*8 vnorm,eel(6),vp(5),eth(6),a,b,c,tot
      eel(:)=0d0
      vp(:)=0d0
      eth(:)=0d0
c     Valid up to 6th dimension for elasticity and thermal expansion/contribution
      do i=1,6
      do iph=1,nph
      do kkk=ngr(iph-1)+1,ngr(iph)
         eel(i)            =eel(i)+wgt(kkk)*dge(i,kkk)
         eth(i)            =eth(i)+wgt(kkk)*dg_trans(i,kkk)
      enddo
      enddo
      enddo


c     Valid up to 5th dimension for visco-plasticity
      do i=1,5
      do iph=1,nph
      do kkk=ngr(iph-1)+1,ngr(iph)
         vp(i)=vp(i)+wgt(kkk)*dgp(i,kkk)
      enddo
      enddo
      enddo

      a=vnorm(eel,6)
      b=vnorm(vp,5)
      c=vnorm(eth,6)

      tot=a+b+c
      if (tot.gt.0) then
         contrib(1)=a/tot
         contrib(2)=b/tot
         contrib(3)=c/tot
      else
         contrib(1)=0.
         contrib(2)=0.
         contrib(3)=0.
      endif

      if (ife.eq.0) then
         write(uw,'(4e14.5e3)',advance='no') (contrib(i),i=1,3),tot
         if (tot.gt.0) then
            write(uw,'(e14.5e3)')(a+b+c)/tot
         else
            write(uw,'(e14.5e3)')0.
         endif
      endif

      return
      end subroutine calc_contribution
