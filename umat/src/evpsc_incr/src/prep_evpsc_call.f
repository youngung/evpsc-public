c     calc ELSC, update schmid tensor, assign gamd0s
c ---------------------------------------------------------------
c     ELSC calculation
c        (inout: celav)
c        (out: ac_el, bc_el, athcs,celscs,mgel,selav,celav)
      if (print_loc) write(*,*) '** entering elsc_evp_fe'
      if (is_all_zero(celav,6)) then !! required when restart the job
         write(*,*)'** Error: all components in celav are zero!'
         call xit()
      endif

      iopt=0
      include "../../../src/vpsc_src/call_elsc_evp_fe.f"
      iopt=1
      include "../../../src/vpsc_src/call_elsc_evp_fe.f"
      if (print_loc) write(*,*) '** Exiting elsc_evp_fe'
c ---------------------------------------------------------------
c     Update schmid tensor (in: schca, ag; out: sch)
      do iphel=1,NPH
         igb=ngr(iphel-1)+1
         igt=ngr(iphel)
         call update_schmid_fe(
c     in
     $        iphel,igb,igt,nphmx,ngrmx,nsysmx,nsyst,schca,
     $        nschca,ag,
c     out
     $        sch,nsch)
      enddo
c ---------------------------------------------------------------
c     gamd0g is replaced by gamd0s that may potentially differ from
c     a slip system to another
c     out: gamd0s
      if(print_loc) write(*,*) '** Begins scale_gam0ds'
c      if(ihardlaw.eq.0) then
c         if(iratesens.eq.0) then
      IF((IHARDLAW.EQ.0.or.IHARDLAW.EQ.20.or.ihardlaw.eq.3)
     $     .and.IRATESENS.EQ.0) THEN
c         write(*,*)'** Error: expected iratesens given'
c         call xit()

         if (vnorm(DBAR,6).eq.0.and.ife.eq.0) then
            write(*,*)'** Unexpected case'
            write(*,*)'DBAR is zero, and cannot scale gamd0s:'
            write(*,*)'DBAR:',DBAR
            call xit()
         endif

         call scale_gamd0s_evp_fe(nsysmx,nphmx,ngrmx,NPH,ngr,
     $        nsyst,DBAR,tg,gamd0s)
c         endif
c      elseif(ihardlaw.ge.20)then !.and.jstep(1).eq.1) then !.and.kinc.eq.1) then
c         call scale_gamd0s_evp_fe(nsysmx,nphmx,ngrmx,NPH,ngr,
c     $        nsyst,DBAR,gamd0s)
      endif

c ---------------------------------------------------------------
c     Initial guess:   (out): DBAR0, DBAR, sbar_old, dsg, sg, sg_old, dsgtry,sbar,mbar,dsbar ...)
      if (vnorm(sbar_old,6).eq.0) then
         if (print_loc) then
            write(*,*)'** Begins initial_guess as sbar_old=0'
         endif
         sbar(:)=0. ! sbar is used in initial_state_guess to estimate stress rate
         include"../../../src/vpsc_src/call_initial_guess.f"
         if (print_loc) then
            write(*,*)'** ends initial_guess'
         endif
      else
         if (print_loc) then
            write(*,*)'** Skipping initial guess, intead extrapolate'
            write(*,*)'   stress increment and stress from rates'
         endif

c        if iregime.eq.2 !! Delta-EVPSC
c        if iregime.eq.3 !! sigma-EVPSC

         if (iregime.eq.2) then
c         mbar(:,:)=mbardot(:,:)*dtime
            dsbar(:)=sbardot(:)*dtime
c         DBAR0(:)=DBAR0dot(:)*dtime
            dsg(:,:)=sgdot(:,:)*dtime
            dsgtry(:,:)=dsg(:,:)
            sg(:,:)=sg_old(:,:)+dsg(:,:)
            sgtry(:,:)=sg(:,:)
         endif
      endif
