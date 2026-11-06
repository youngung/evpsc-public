      subroutine rexdd_stats(
c          intent(in):
c          * integer
     $     nphmx,iph,

     $     ngr,
     $     ngrmx,nsysmx,nmodmx,nslmod,nsm,

     $     ncluster,nclustermx,nchildset,

     $     nuc_kount,kount_van_all,
     $     ngref,igref,

     $     uw_rex,
c          * real
     $     time,temp,epsvm,svm,
     $     wgt,
     $     radius,velgrow,
     $     rho_s_for,rho_s_deb,crss,pnucl,

c          intent(out):
c          * integer
     $
c          * real
     $     dwgt_grow,dwgt_nucl,

c          intent(out)::
     $     drxfrac
     $     )
      implicit none
c     intent(in)
c     * integer
      integer, intent(in):: nphmx,iph,
     $     ngr(0:nphmx),
     $     ngrmx,nsysmx,nmodmx,
     $     nslmod(nphmx),nsm(nmodmx,nphmx),

     $     ncluster,
     $     nclustermx,nchildset(nclustermx),

     $     nuc_kount,kount_van_all,ngref,igref(ngrmx),

     $     uw_rex
c     * real
      real*8, intent(in)::
     $     time,temp,epsvm,svm,
     $     wgt(ngrmx),dwgt_grow(ngrmx),dwgt_nucl(ngrmx),
     $     radius(ngrmx),velgrow(ngrmx),
     $     rho_s_for(nsysmx,ngrmx),rho_s_deb(nsysmx,ngrmx),
     $     crss(nsysmx,ngrmx),pnucl(ngrmx)
c     intent(out)
c     * integer
c     * real
      real*8, intent(out):: drxfrac
c     locals
      logical isnan
      integer ind,kkk,imo,is,ism,icluster,nch,ngrph
      real*8 wgtmp,wgtph,wgt_new(ngrmx),avgmem,wgt_ref,
     $     radius_wgt_avg,radius_ari_avg,
     $     radius_wgt_std,radius_ari_std,
     $     rho_wgt_avg,rho_wgt_std,rho_g(ngrmx),
     $     pnu_wgt_avg,pnu_wgt_std,
     $     velgrow_ari_avg,velgrow_wgt_avg,
     $     velgrow_ari_std,velgrow_wgt_std
      real*8 crss_mn,crss_mx,
     $     crss_wgt_avg,crss_wgt_std,
     $     crss_ari_avg,crss_ari_std,crss_sum(ngrmx)

c     --- new volume fraction (a local array, not passed)
      ngrph=ngr(iph)-ngr(iph-1)
      wgt_new(:)=0.
      wgtph=0.
      do kkk=ngr(iph-1)+1,ngr(iph)
         wgt_new(kkk)=wgt(kkk)+dwgt_grow(kkk)+dwgt_nucl(kkk)
         wgtph=wgtph+wgt_new(kkk)
      enddo

c     -- disl. density for each grain
      do kkk=ngr(iph-1)+1,ngr(iph)
         rho_g(kkk)=0.
         ism=0
         do imo=1,nslmod(iph)
         do is=1,nsm(imo,iph)
            ism=ism+1
            rho_g(kkk)=rho_g(kkk)
     $           +rho_s_for(ism,kkk)+rho_s_deb(ism,kkk)
         enddo
         enddo
      enddo

c     -- CRSS
      crss_mn=1e10
      crss_mx=-1e10
      crss_wgt_avg=0.
      crss_ari_avg=0.
      crss_wgt_std=0.
      crss_ari_std=0.
      do kkk=ngr(iph-1)+1,ngr(iph)
         ism=0
         crss_sum(kkk)=0.
         do imo=1,nslmod(iph)
         do is=1,nsm(imo,iph)
            ism=ism+1
            if (crss(ism,kkk).lt.crss_mn) crss_mn=crss(ism,kkk)
            if (crss(ism,kkk).gt.crss_mx) crss_mx=crss(ism,kkk)
            crss_sum(kkk)=crss_sum(kkk)+crss(ism,kkk)
         enddo
         enddo
         crss_sum(kkk)=crss_sum(kkk)/ism
         crss_wgt_avg=crss_wgt_avg+crss_sum(kkk)*wgt(kkk)
         crss_ari_avg=crss_ari_avg+crss_sum(kkk)
      enddo
      crss_wgt_avg=crss_wgt_avg/wgtph
      crss_ari_avg=crss_ari_avg/ngrph

      crss_wgt_std=0.
      crss_ari_std=0.
      do kkk=ngr(iph-1)+1,ngr(iph)
         crss_wgt_std=crss_wgt_std+(crss_wgt_avg-crss_sum(kkk))**2*
     $        wgt(kkk)
         crss_ari_std=crss_ari_std+(crss_ari_avg-crss_sum(kkk))**2
      enddo
      crss_wgt_std=sqrt(crss_wgt_std/wgtph)
      crss_ari_std=sqrt(crss_ari_std/ngrph)

c     --- DRXFRAC - the volume fraction of recrystalized grains.
      drxfrac=1d0
      do ind=1,ngref
         kkk=igref(ind)
         wgtmp=wgt(kkk)+dwgt_nucl(kkk)+dwgt_grow(kkk)
         drxfrac=drxfrac-wgtmp
      enddo

      if (drxfrac.lt.0) drxfrac=0.

c     --- grain radii and radial migration speed (velgrow)
      radius_wgt_avg=0. ! weighted avg
      radius_ari_avg=0. ! arithmetic avg
      velgrow_wgt_avg=0.
      velgrow_ari_avg=0.
      do kkk=ngr(iph-1)+1,ngr(iph)
         radius_ari_avg=radius_ari_avg+radius(kkk)
         radius_wgt_avg=radius_wgt_avg+radius(kkk)*wgt_new(kkk)

         velgrow_ari_avg=velgrow_ari_avg+abs(velgrow(kkk))
         velgrow_wgt_avg=velgrow_wgt_avg+abs(velgrow(kkk))*wgt_new(kkk)
      enddo                     ! loop over grains

      radius_ari_avg=radius_ari_avg/ngrph
      radius_wgt_avg=radius_wgt_avg/wgtph
      velgrow_ari_avg=velgrow_ari_avg/ngrph
      velgrow_wgt_avg=velgrow_wgt_avg/wgtph

      radius_wgt_std=0.
      radius_ari_std=0.
      velgrow_wgt_std=0.
      velgrow_ari_std=0.
      do kkk=ngr(iph-1)+1,ngr(iph)
         radius_wgt_std=radius_wgt_std+
     $        (radius(kkk)-radius_wgt_avg)**2*wgt_new(kkk)
         radius_ari_std=radius_ari_std+
     $        (radius(kkk)-radius_ari_avg)**2

         velgrow_wgt_std=velgrow_wgt_std+
     $        (velgrow_wgt_avg-abs(velgrow(kkk)))**2*wgt_new(kkk)
         velgrow_ari_std=velgrow_ari_std+
     $        (velgrow_ari_avg-abs(velgrow(kkk)))**2
      enddo
      radius_wgt_std=sqrt(radius_wgt_std/wgtph)
      radius_ari_std=sqrt(radius_ari_std/ngrph)
      velgrow_wgt_std=sqrt(velgrow_wgt_std/wgtph)
      velgrow_ari_std=sqrt(velgrow_ari_std/ngrph)

c     --- avg mem
      avgmem=0.
      do icluster=1,ncluster
         nch=nchildset(icluster)
         avgmem=avgmem+nch+1
      enddo
      avgmem=avgmem/(ncluster*1.)

c     --- dislocation stats
c     forward
      rho_wgt_avg=0.
      do kkk=ngr(iph-1)+1,ngr(iph)
         rho_wgt_avg=rho_wgt_avg+rho_g(kkk)*wgt_new(kkk)
      enddo
      rho_wgt_avg=rho_wgt_avg/wgtph

      rho_wgt_std=0.
      do kkk=ngr(iph-1)+1,ngr(iph)
         rho_wgt_std=rho_wgt_std+(rho_wgt_avg-rho_g(kkk))**2
     $        *wgt_new(kkk)
      enddo
      rho_wgt_std=sqrt(rho_wgt_std/wgtph)

c     --- pnucl avg, std,
      pnu_wgt_avg=0.
      pnu_wgt_std=0.
      wgt_ref=0.
      do kkk=ngr(iph-1)+1,ngr(iph)
         if (dwgt_nucl(kkk).le.0) then
            wgt_ref=wgt_ref+wgt_new(kkk)
            pnu_wgt_avg=pnu_wgt_avg+pnucl(kkk)*wgt_new(kkk)
         endif
      enddo
      if (wgt_ref.gt.0) then
         pnu_wgt_avg=pnu_wgt_avg/wgt_ref
         do kkk=ngr(iph-1)+1,ngr(iph)
            if (dwgt_nucl(kkk).le.0) then
               pnu_wgt_std=pnu_wgt_std
     $              +(pnu_wgt_avg-pnucl(kkk))**2*wgt_new(kkk)
            endif
         enddo
         pnu_wgt_std=sqrt(pnu_wgt_std/wgtph)
      else
         pnu_wgt_avg=0.
         pnu_wgt_std=0.
      endif

c     --- writing stats
      write(uw_rex,'(5f10.4)',advance='no') time,temp,epsvm,svm,drxfrac
      write(uw_rex,'(2f7.3)',advance='no') pnu_wgt_avg,pnu_wgt_std
      write(uw_rex,'(i7)',advance='no')ngrph
      write(uw_rex,'(i7)',advance='no') ncluster
      write(uw_rex,'(i7)',advance='no') nuc_kount
      write(uw_rex,'(i7)',advance='no') kount_van_all
      write(uw_rex,'(f7.1)',advance='no') avgmem

      write(uw_rex,'(e11.3)',advance='no') radius_ari_avg*1e6
      write(uw_rex,'(e11.3)',advance='no') radius_wgt_avg*1e6
      write(uw_rex,'(e11.3)',advance='no') radius_ari_std*1e6
      write(uw_rex,'(e11.3)',advance='no') radius_wgt_std*1e6


      write(uw_rex,'(e11.3)',advance='no') velgrow_ari_avg
      write(uw_rex,'(e11.3)',advance='no') velgrow_wgt_avg
      write(uw_rex,'(e11.3)',advance='no') velgrow_ari_std
      write(uw_rex,'(e11.3)',advance='no') velgrow_wgt_std

      write(uw_rex,'(e11.3)',advance='no') rho_wgt_avg
      write(uw_rex,'(e11.3)',advance='no') rho_wgt_std

      write(uw_rex,'(6e11.3)',advance='no') crss_mn,crss_mx,
     $     crss_wgt_avg,crss_wgt_std,
     $     crss_ari_avg,crss_ari_std

c     line breaker
      write(uw_rex,*)

      write(*,'(2a10, 10x ,a10)')'GB vel_w','GB vel_a','Pnucl'
      write(*,'(2f10.2,a10,e10.2)')
     $     velgrow_wgt_avg*1e6,
     $     velgrow_ari_avg*1e6,'[um/s]',
     $     Pnu_wgt_avg



      if (
     $     isnan(radius_ari_avg)
     $     .or.isnan(radius_wgt_avg)
     $     .or.isnan(radius_ari_std)
     $     .or.isnan(radius_wgt_std)
     $     ) then
         write(*,*)'nan found!!'
         write(*,*)'wgtph:',wgtph
         write(*,*)'ngrph:',ngrph
         do kkk=ngr(iph-1)+1,ngr(iph)
            if (isnan(radius(kkk))) then
               write(*,*)'nan was found in radius for grain',kkk
               write(*,*)'radius(kkk):',radius(kkk)
            endif
         enddo
         call stop_debug(0,'** Stopped since nan is found')
      endif

      return
      end subroutine rexdd_stats
