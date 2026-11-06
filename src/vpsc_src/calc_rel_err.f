c     -----------------------------------------------------------------
      subroutine calc_rel_err(
c          intent(in)
     $     iregime,nphmx,ngrmx,nph,ngr,dsbar,sbar,dbar,
     $     dg,dsg,sg,sg_old,dsgtry,sgtry,wgt,
c          intent(inout)
     $     sbarc,dbarc,
c          intent(out)
     $     relsgr,rels,reld,rels6,reld6,sav,dsav,dav)
      implicit none
      integer,intent(in)::iregime,nphmx,ngrmx,nph,ngr(0:nphmx)
      real*8,intent(in)::dsbar(6),sbar(6),dbar(6),dg(6,ngrmx),
     $     dsg(6,ngrmx),
     $     sg(6,ngrmx),sg_old(6,ngrmx),dsgtry(6,ngrmx),sgtry(6,ngrmx),
     $     wgt(ngrmx)
      real*8, intent(inout)::sbarc(3,3),dbarc(3,3)
      real*8, intent(out)::relsgr,rels,reld,rels6(6),reld6(6),dsav(6),
     $     sav(6),dav(6)
c     -- local
      integer i,j,kkk,mdim,iph
      real*8 aux6(6),sdev,vmismatch,vnorm,strain_rate_sum,stress_sum,
     $     near_zero_mx,snorm,dnorm,tnorm
      parameter(near_zero_mx=1e-6)
      mdim=6
      if (iregime.eq.1) mdim=5

      sdev=0.
      dsav(:)=0.
      sav(:)=0.
      if (iregime.eq.1) then
         do i=1,6
         do iph=1,nph
         do kkk=ngr(iph-1)+1,ngr(iph)
            sav(i) =sav(i)+sg(i,kkk)*wgt(kkk)
            SDEV=SDEV+((sg(I,KKK)-sgtry(I,kkk))**2)*WGT(KKK)
         enddo
         enddo
         enddo
      elseif (iregime.eq.2) then
         do i=1,6
         do iph=1,nph
         do kkk=ngr(iph-1)+1,ngr(iph)
            dsav(i)=dsav(i)+dsg(i,kkk)*wgt(kkk)
            sav(i) =sav(i)+(dsg(i,kkk)+sg_old(i,kkk))*wgt(kkk)
            SDEV=SDEV+((dsg(I,KKK)-dsgtry(I,kkk))**2)*WGT(KKK)
         enddo
         enddo
         enddo
      elseif(iregime.eq.3) then
         do i=1,6
         do iph=1,nph
         do kkk=ngr(iph-1)+1,ngr(iph)
            dsav(i)=dsav(i)+(sg(i,kkk)-sg_old(i,kkk))*wgt(kkk)
            sav(i) =sav(i)+sg(i,kkk)*wgt(kkk)
            SDEV=SDEV+((sg(I,KKK)-sgtry(I,kkk))**2)*WGT(KKK)
         enddo
         enddo
         enddo
      endif

      dav(:)=0.
      do i=1,6
      do iph=1,nph
      do kkk=ngr(iph-1)+1,ngr(iph)
         dav(i)=dav(i)+dg(i,kkk)*wgt(kkk) ! <edot>
      enddo
      enddo
      enddo

C *** ROUNDS-UP TENSOR COMPONENTS
      SNORM=TNORM(SBARc,3) *1.E-6
      DNORM=TNORM(DBARc,3) *1.E-6
      DO I=1,3
      DO J=1,3
        IF(ABS(SBARc(I,J)) .LT. SNORM) SBARc(I,J)=0.
        IF(ABS(DBARc(I,J)) .LT. DNORM) DBARc(I,J)=0.
      ENDDO
      ENDDO

c     ** check consistency
c     1. RELSGR
      RELSGR=SQRT(SDEV)/VNORM(SBAR,MDIM)

c     2. RELD <sum_g{edot^g}w^g>=dbar
      aux6(:)=dbar(:)+dav(:)
      strain_rate_sum=vnorm(aux6,6)
      if (strain_rate_sum.lt.near_zero_mx) then
         reld=vnorm(aux6,6)
      else
         RELD=VMISMATCH(DBAR,DAV,MDIM)
      endif

c     3. RELS
c     Overrides the RELS definition in case near zero
c     stress (or stress increment)

      if (iregime.eq.1) then
c        This is for VPSC.
         aux6(:)=sav(:)+sbar(:)
         stress_sum=vnorm(aux6,6)
         if (stress_sum.lt.near_zero_mx) then
            rels=vnorm(aux6(1:mdim),mdim)
         else
            rels=vmismatch(sav(1:mdim),sbar(1:mdim),mdim)
         endif
      elseif (iregime.eq.2) then
c        This is for delta-EVPSC
         aux6(:)=sav(:)+sbar(:)
         stress_sum=vnorm(aux6,6)
         if (stress_sum.lt.near_zero_mx) then
            rels=vnorm(aux6(1:mdim),mdim)
         else
            rels=vmismatch(dsav(1:mdim),dsbar(1:mdim),mdim)
         endif
      endif


c     4. Component-wise relative error
      do i=1,6
c        4-1. Stress (either increment or full)
         if (iregime.eq.2) then ! delta
            if (dsav(i).eq.0.) then
               rels6(i)=abs(dsav(i)-dsbar(i))
            else
               rels6(i)=abs((dsav(i)-dsbar(i))/dsav(i))
            endif
         elseif(iregime.eq.1.or.iregime.eq.3) then ! VPSC or sigma
            if (sav(i).eq.0.) then
               rels6(i)=abs(sav(i)-sbar(i))
            else
               rels6(i)=abs((sav(i)-sbar(i))/sav(i))
            endif
         endif
c        4-2. Strain rate
         if (dav(i).eq.0.) then
            reld6(i)=abs(dav(i)-dbar(i))
         else
            reld6(i)=abs((dav(i)-dbar(i))/dav(i))
         endif
      enddo

      return
      end subroutine calc_rel_err
