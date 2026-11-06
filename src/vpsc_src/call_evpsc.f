c----------------------------------------------------------------------c
c     This calling block of subroutine <initial_state_guess>  is shared
c     by Stand-alone and FE-interface versions.
c----------------------------------------------------------------------c
c     Argument order: Numerical CONTROL & CONDITION, DIMENSION, PX, PHASE, GRAIN, SYSTEM
c     Then, if applicable, alphabetical order is used.
      idiag_nr=0
      xnumstats(:)=0.
      call chg_basis(dsbar,dsbarc,aux66,aux3333,1,6)


      ! deps_ref = 1e-5           !! reference strain increment.
      write(*,*)'dtime:',dtime
c      dtime_ref=0.02
      dtime_ref=dtime

      if (print_loc) write(*,*)"** call subroutine evpsc"
      CALL evpsc(
c          intent(in)
c          ** logical in
     $     print_loc,
c          ** integers in
     $     mdim,ife,iregime,interaction,iverbose,ioptfs,itmaxint,
     $     itmaxext,UW_LOG,UW_SCCN,ISTEP,nphmx,ngrmx,nsysmx,iDBARv,
     $     iSBARv,NPH,ngr,nsyst,ishape,isense,nrs,
c          ** floats in
     $     errs,errd,errm,dtime_ref,fstar,alplpx,c6,sbar_old,axisph,
     $     axisgr,wgt,sg_old,dg_trans,mgel,crss,gamd0s,sch,nsch,
c          intent(inout)
c             ** floats inout
     $     mbar,LBAR,DBAR0,
     $     sbar,dsbar,sbarc,dsbarc,DBAR,sgtry,dsgtry,sg,
     $     dsg,dg0,mgtg,
c          intent(out)
c             ** integers out
     $     it2,idiag_sc,idiag_sg,idiag_nr,
c             ** floats out
     $     xnumstats,dbarc,factor,fel,fvp,alpl,dge,dgp,dg,tg,gamdot,
     $     asph,asgr,
c     I am not sure if these are required to be passed out as an argument
     $     sav,dav)
      if (print_loc) write(*,*)"** right after evpsc 1"

      relsgr=xnumstats(1)
      rels=xnumstats(2)
      reld=xnumstats(3)
c     nritmn=xnumstats(4) !! minimum NR iterations of the last iteration loop.
c     nritmx=xnumstats(5) !! maximum NR iterations of the last iteration loop.
c     nritavg=xnumstats(6) !! avg.   NR iterations of the last iteration loop.
c     nrffrac=xnumstats(7) !! The fraction of grains for which NR failed to converge for the last iteration loop.

c      if (ife.eq.0) then
c         write(*,'(1H+,x,3e10.3,x,4f8.2)')
c     $        (xnumstats(i),i=1,7)
c      endif

      if (idiag_nr.ne.0.or.idiag_sc.ne.0.and.iregime.eq.2) then
         if (idiag_nr.ne.0) write(*,*)'*** NR failed! for fstar',fstar
         if (idiag_sc.ne.0) write(*,*)'*** SC failed! for fstar',fstar
         write(*,'(a10,6f9.3)')'dsg:',dsg(:,1)
         if (ife.eq.0) then
            call xit()
         elseif(ife.eq.1.and.idiag_nr.ne.0)then
            pnewdt=0.5
            return
         endif
      endif

      sgdot(:,:)=dsg(:,:)/dtime_ref
      sbardot(:)=dsbar(:)/dtime_ref
      mbardot(:,:)=mbar(:,:)/dtime_ref
      DBAR0dot(:)=DBAR0(:)/dtime_ref

      if (iregime.eq.2) then !dEVPSC
         alplpx=0.
         do iph=1,NPH
            do kkk=ngr(iph-1)+1,ngr(iph)
               alplpx=alplpx+alpl(kkk)*wgt(kkk)
            enddo
         enddo
         if (isnan(alplpx)) alplpx=0.
      elseif (iregime.eq.1) then ! VPSC
c        VPSC is 100% fully plastic always.
         alplpx=1.
      endif

      call chg_basis(sbar,scauchy,aux66,aux3333,1,6)
      call chg_basis(sbar,sdeviat,aux66,aux3333,1,5)

      if (iregime.eq.2) then ! delta-EVPSC
         dsbar(:)=sbardot(:)*dtime
         dsg(:,:)=sgdot(:,:)*dtime

         sbar(:)=sbar_old(:)+dsbar(:)
         sg(:,:)=sg_old(:,:)+dsg(:,:)

         call chg_basis(dsbar,dsbarc,aux66,aux3333,1,6)

         write(*,*)'dsbar:',dsbar(1)
         write(*,*)'dsg(1,1):',dsbar(1)
         write(*,*)'sbar:',sbar(1)
         write(*,*)'sg(1,1):',sbar(1)
c         stop -1

c      else(iregime.eq.1) then !VPSC
c        no need to integrate stress
      endif
