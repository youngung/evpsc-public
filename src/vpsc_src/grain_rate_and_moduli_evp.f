c ****************************************************************************
c     Subroutine grain_rate_and_moduli_evp for delta-EVPSC and sigma-EVPSC
c     - Youngung Jeong
c ****************************************************************************
      subroutine grain_rate_and_moduli_evp(
c          in
     $     interaction,iregime,ioptfs,nrs,gamd0s,taux,mgel,isense,
     $     sch,nsch,c6,sg_old,dsg,sg,nsyst,nsysmx,dtime,dg_trans,
c          out
     $     dge,dgp,tg,gamdot,mgtg,dg0
     $     )
c     ------------------------------------------------------------------------
c     one can distinguish the use of the proper id based on how the associated
c     array is defined in the dimension file
      implicit none
      integer,intent(in)::interaction,nsysmx,iregime,nsyst,ioptfs,
     $     isense(nsysmx),nrs(nsysmx)
      real*8,intent(in)::dtime,sg_old(6),dsg(6),sg(6),sch(5,nsysmx),
     $     nsch(6,nsysmx),c6,gamd0s(nsysmx),taux(nsysmx),mgel(6,6),
     $     dg_trans(6)
      real*8,intent(out)::gamdot(nsysmx),dge(6),dgp(5),tg(6),
     $     mgtg(6,6),dg0(6)
c     0: Using the linear form of Jeong et al. (IJP, 2021) and Jeon et al. (JMRT, 2024)
c     1: Using the new linear form of Patra et al.
c     ---- locals
      real*8 mgsec(6,6),mgtgvp(6,6),mgsecvp(6,6),mgsecevp(6,6),
     $     rss(nsyst),rssx,hyd,signhyd
      integer is,i,j,xid(6,6),mdim

      if (iregime.ne.1.and.iregime.ne.2) then
         write(*,*)'grain_rate_and_moduli_evp works for VPSC or dEVPSC'
         call xit()
      endif

      mdim=6
      if (iregime.eq.1) mdim=5 ! VPSC

      xid(:,:)=0 ! kronecker delta
      do i=1,6
         xid(i,i)=1
      enddo

c     ** gets resolved shear stresses 'rssx' and shear rates 'gamdot'.
c     sign(gamdot)=sign(rssx).
c     nrs can be even or odd.
c     rss is always > 0 and is used to calculate elasto-viscous compliance.
      if (iregime.eq.2) then
         hyd=(dsg(6)+sg_old(6))/3.
      else
         hyd=0.
      endif
      signhyd=sign(1.,hyd)
      gamdot(:)=0.


      do is=1,nsyst
         rssx=0d0
         if(iregime.eq.1) then ! VPSC
            do i=1,5
               rssx=rssx+(sch(i,is)+nsch(i,is))*sg(i)
            enddo
         elseif (iregime.eq.2) then !! dEVPSC
            do i=1,5
               rssx=rssx
     $              +(sch(i,is)+nsch(i,is))*(dsg(i)+sg_old(i))
            enddo
         endif
c        ! (ds + s(n))
         if (iregime.eq.2) then
            rssx=rssx+nsch(6,is)*(dsg(6)+sg_old(6))
         endif
         rssx=rssx+c6*abs(hyd)
         if(.not.(rssx.gt.0 .or. isense(is).eq.1)) rssx=0d0
         rssx=rssx/taux(is)
         rss(is)  =gamd0s(is)*abs(rssx**(nrs(is)-1))/taux(is)
         gamdot(is)=gamd0s(is)*abs(rssx**nrs(is))*sign(1.,rssx)
      enddo

c     * calculate pure visco-plastic strain rate
c     that is the stress dependent term
      dgp(:)=0d0                ! due only to slip (or twin)
      do i=1,5                  ! deviatoric thus 5d.
      do is=1,nsyst
         dgp(i)=dgp(i)+sch(i,is)*gamdot(is)
      enddo
      enddo
      dgp(:)=dgp(:)

c     * calculate total VP strain rate
      tg(:)=0d0
      tg(:5)=dgp(:5)
      tg(:)=tg(:)+dg_trans(:)

c     ** calculate elastic strain rate using local elastic laws
c      dge = m:(sg - sg_old)/dt
c      dge = m:dsg / dt
      dge(:)=0d0
      do i=1,6
      do j=1,6
         dge(i)=dge(i)+mgel(i,j)*dsg(j)/dtime ! ds
      enddo
      enddo

c     Add elastic strain rate, only when EVPSC is used.
      if (iregime.eq.2) tg(:)=tg(:)+dge(:)

c     ** calculate the ratio of visco-plastic strain rate to overall strain rate


c *** calculate crystal compliance m_tg, m_sc, m_evp, and eps_dot_0
      mgsec(:,:)=0.
      mgtg(:,:)=0. !! tangent (affine) EVP compliance
      mgtgvp(:,:)=0. !! it is 6x6 to consider non-Schmid effect
      mgsecvp(:,:)=0. !! it is 6x6 to consider non-Schmid effect
      mgsecevp(:,:)=0.

c     Visco-plastic compliance tensors
c     (note), EVPSC compliance of delta-EVPSC is equivalent to that of sigma-EVPSC
      do i=1,5
c        loop over 1..6. 6th is included to account for the non-schmid effect
      do j=1,6
      do is=1,nsyst
         if (j.eq.6) then
            mgsecvp(i,j)=mgsecvp(i,j)
     $           +sch(i,is)*(nsch(j,is)+c6/3.*signhyd)*rss(is)
            mgtgvp(i,j)=mgtgvp(i,j)
     $           +sch(i,is)*(nsch(j,is)+c6/3.*signhyd)*rss(is)*nrs(is)
         else
            mgsecvp(i,j)=mgsecvp(i,j)+sch(i,is)*(sch(j,is)+nsch(j,is))
     $           *rss(is)
            mgtgvp(i,j)=mgtgvp(i,j)+sch(i,is)*(sch(j,is)+nsch(j,is))
     $           *rss(is)*nrs(is)
         endif
      enddo ! loop over is=1,nsyst
      enddo ! loop over j=1..6
      enddo ! loop over i=1..5

c *** EVP tangent & EVP secant compliances not vp compliances!
c *** Also, for ioptfs.eq.6 and 7, mgsecvp and mgtgvp are mgsecvp* and mgtgvp*


c ***
c *** Complete the linearized form (MGTG and DG0)
c ***

      if (interaction.eq.1.or.interaction.eq.5) then
c        ** Determine the tangent compliance tensor
c        VPSC
         if (iregime.eq.1) then
            mgtg(:,:)=mgtgvp(:,:)
         elseif(iregime.eq.2) then
            if (ioptfs.eq.0) then
               mgtg(:,:)=mgtgvp(:,:)+mgel(:,:)/dtime
            elseif(ioptfs.eq.1) then
               mgtg(:,:)=mgsecvp(:,:)+mgel(:,:)/dtime
            endif
         endif

c        ** Determine the back-extrapolated strain rate tensor
         dg0(:)=0d0
c        VPSC
         if (iregime.eq.1) then
            do i=1,5
            do j=1,5
               dg0(i)=dg0(i)+(mgsecvp(i,j)-mgtgvp(i,j))*sg(j)
            enddo
            enddo
c        dEVPSC
         elseif(iregime.eq.2) then
            dg0(:)=0
            if(ioptfs.eq.0) then
               do i=1,6
                  do j=1,6
                     dg0(i)=dg0(i)+mgsecvp(i,j)*(sg_old(j)+dsg(j))
     $                    +mgel(i,j)/dtime*dsg(j)
     $                    -mgtg(i,j)*dsg(j)
                  enddo
               enddo
            elseif (ioptfs.eq.1) then
               do i=1,5
                  do j=1,5
                     dg0(i)=dg0(i)+mgsecvp(i,j)*sg_old(j)
                  enddo
               enddo
            endif
         endif !! end of regime
      elseif(interaction.eq.0.or.interaction.eq.2.or.
     $        interaction.eq.3.or.interaction.eq.4
     $        .and.iregime.eq.1) then
         dg0(:)=0.
         do i=1,5
         do j=1,5
            mgtg(i,j)=mgsecvp(i,j)
         enddo
         enddo
      else
         write(*,*)'Unexpected combintation of interaction and iregime'
         write(*,*)'interaction:',interaction
         write(*,*)'iregime:',iregime
         call xit()
      endif





      return
      end subroutine grain_rate_and_moduli_evp
