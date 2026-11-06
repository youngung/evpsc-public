c     a) Update fij, shape for phase zero
      if (print_loc) write(*,*)'** begins update_fij_fe (0)'
c     update fijph but for 0-phase (that's the average deformation of inclusions)
      CALL UPDATE_FIJ_fe(
c     in
     $     nphmx,ngrmx,0,0,0,dtime,LIJBARc,ishape,xlijgr,
c     out
     $     fijph,fijgr)         ! UPDATES DEFORMATION TENSOR OF ELEMENT
      if (print_loc) write(*,*)'** begins update_shape_fe'

c     update inclusion shapes using Fij by F^(new)_ij=(Kronecker_ik + LIJBARc_ik.dt).F^(old)_kj
c     use the updated fphij to obtain axisph, axisgr, eulerph, fijph, fijgr
      call update_shape_fe(
c          in
     $     ngrmx,nphmx,0,0,0,interaction,ifrag,crit_shp,ishape,
c          inout
     $     fijph,iflat,axisph,axisgr,
c          out
     $     eulerph,fijgr,ierr_shp)

      if (ierr_shp.ne.0)then ! failed to perform shape_fe
         call mutexlock(5)
         write(*,*)'---------------------------------------------------'
         write(*,*)'*Failed to complete UPDATE_SHAPE for PHASE ELPSD(0)'
         write(*,*)'---------------------------------------------------'
         call mutexunlock(5)
         pnewdt=0.5
         return
      endif


c     b) Update twinned crystal's orientation
      if (print_loc) write(*,*)'** begins update_twinning'
      DO IPHEL=1,NPH
         IF(NTWMOD(IPHEL).NE.0) THEN
            igb=ngr(iphel-1)+1
            igt=ngr(iphel)
            if (itwinlaw.le.2) then
               call update_twinning_ptr_fe(
c                   in
     $              dtime,iphel,igb,igt,nphmx,wph,igt-igb+1,ngrmx,
     $              nslsys,ntwmod,nslmod,twsh,nsm,nmodmx,ntwmmx,nsysmx,
     $              gamdot,ntwsmx,wgt,itwinlaw,jran,twthres,dnca,dbca,
     $              ntwsys,isectw,itwtype,scauchy,
c                   inout
     $              twfrsy,twfrph,ktwsmx,eftwfr,ag,dsg,ntwevents,pritw,
     $              sectw,sch_stat,var_stat,var_vf
c!c                   out
c!     $              ,TWFRGR
     $              )
            else
                 call update_twinning_vft_fe(
c                     in
     $                dtime,iphel,nphmx,ngrmx,ngr,nslmod,nslsys,nsysmx,
     $                ntwmod,ntwmmx,ntwsmx,twthres,nsm,gamdot,nmodmx,
     $                dnca,dbca,itwtype,nsyst,
c                     inout
     $                ktwsmx,twsh,twfrsy,twfrph,sg_old,sg,dsg,wph,
     $                ag,wgt,gtotgr,crss,
c                     out
     $                pritw)
            endif
         ENDIF
      ENDDO
c     c) Update orientation
      if (print_loc) write(*,*)'** begins update_orientation'
      do iphel=1,NPH
         igb=ngr(iphel-1)+1
         igt=ngr(iphel)
         IF(IUPDORI.EQ.1) then
c           rot = (WBARc+rotloc-rotslip)*dtime
c                  * WBARc: spin externally imposed for entire polycrystal
c                  * rotslip = 1/2{u^c -(u^c)^T}
c                  * rotloc = as * (tg - DBAR)
c           ag^(new) = ag^(old) + rot
c
c           xlijgr = WBARc + rotloc + 1/2 *{u^c + (u^c)^T}
c                  ~ spin(s) + strain rate
c           xlijph = <xlijgr>
            CALL UPDATE_ORIENTATION_fe(
c                  in
     $           dtime,interaction,nsysmx,ngrmx,nphmx,igb,igt,iphel,
     $           ishape,wph,WBARc,asph,asgr,tg,
     $           DBAR,dnca,dbca,gamdot,wgt,nneigh,nneimx,
     $           neigh,wneigh,
c                  inout
     $           ag,
c     out
     $           xlijph,xlijgr)
         else
            call stop_debug(0,'** Error: texture is not updated?')
            call xit()
         endif
      enddo

c     d) Update fij, shape for non-zero phases (or grains)
      if (print_loc) write(*,*)'** begins update_fij_fe and shape'
      DO IPHel=1,NPH       !IPHBOT,IPHTOP
         igb=ngr(iphel-1)+1
         igt=ngr(iphel)
c          ** for non-zero phase, use xlijph not LIJBARc
         aux33(:,:)=xlijph(:,:,iphel)
c        update fijph for each phase
         if (print_loc) then
            write(*,*)'xijph'
            write(*,'(3e12.4)')((aux33(i,j),j=1,3),i=1,3)
         endif
         CALL UPDATE_FIJ_fe(
c             in
     $        nphmx,ngrmx,iphel,igb,igt,dtime,aux33,ishape,xlijgr,
c               out
     $        fijph,fijgr)
         IF(IUPDSHP.EQ.1) then
            if (print_loc) write(*,*)'** begins update_shape_fe'
c           update inclusion shapes using Fij by F^(new)_ij=(Kronecker_ik + LIJBARc_ik.dt).F^(old)_kj
            CALL UPDATE_SHAPE_fe(
c             in
     $           ngrmx,nphmx,iphel,ngr(iphel-1)+1,ngr(iphel),
     $           interaction,ifrag,crit_shp,ishape,
!     inout
     $           fijph,iflat,axisph,axisgr,
c     out
     $           eulerph,fijgr,ierr_shp)



      if (ierr_shp.ne.0)then ! failed to perform shape_fe
         call mutexlock(5)
         write(*,*)'---------------------------------------------------'
         write(*,*)'*Failed to complete UPDATE_SHAPE for PHASE ELPSD(1)'
         write(*,*)'---------------------------------------------------'
         pnewdt=0.5
         return
      endif



         endif
      enddo

c     e) update temperature dependent CTEs and elastic modulus
      if (.false.) then
      do iphel=1,NPH     !iphelbot,ipheltop
         call thermo_emod_cte(
c     intent(in)
     $        -1,UR_ECTE,iph,nsysmx,nmodmx,nphmx,nmodes,mxndeg,nsm,
     $        temp,schca,
c     intent(inout) !! should be determined inside <evpsc_input_fe>
     $        ndegs_fit,cijv_fit,cte_fit,
c     !! these are updated once at <evpsc_input_fe>
c     !! and should be updated in the main subroutine <EVPSC_MAIN> or <(V)UMAT>
     $        celcc,selcc,xmu_mode,athcc)
      enddo
      endif


c     f) update CRSS
      if (print_loc) write(*,*)'** begins update_crss'
      IF(IUPDHAR.EQ.1) THEN
C     *** VOCE HARDENING PLUS PREDOMINANT TWIN REORIENTATION SCHEME
         do iphel=1,NPH    !iphelbot,ipheltop
            igb=ngr(iphel-1)+1
            igt=ngr(iphel)
            IF(IHARDLAW.EQ.0) then
               aux6(1)=1d0

c              gtotgr(new) = gtotgr(old) + sum_s|gamdot_s|dt
c              crss(new)=crss(old)+dcrss
               CALL UPDATE_CRSS_VOCE_fe(
     $              UR_SX,2,iphel,nphmx,ngrmx,dtime,
c                     in
     $              nmodmx,ntwmmx,igb,igt,nmodes,nsysmx,nsm,
     $              nslsys,nsyst,gamdot,
c                   inout
     $              grsze,gtotgr,tau_voce,thet_voce,hpfac_voce,hard,
     $              crss,
c                    out (but valid only for ioption.eq.0)
     $              twthres,nrs,itwinlaw,isectw)
            elseif (ihardlaw.eq.23) then
               dum1=0d0
               dum2=0d0
               dum3=0d0
               dum4=0d0
               call update_crss_rgvb_fe(
c                   in
     $              UR_SX,2,iphel,igb,igt,dtime,nsysmx,ngrmx,nphmx,
     $              nmodmx,ntwmmx,nsyst,nrs,nslmod,ntwmod,nslsys,ntwsys,
     $              nmodes,nsm,igt-igb+1,dnca,jran,wgt,
c                   inout
     $              TAUE,crss,gamdot,ishear_tag,ibs_law,xMu,Xb,xf_self,
     $              xf_lat,xD,xK,alpham,xlambda,tau0,xp,xbx,facs,pow,
     $              xRho_min,xRho_max,xRho_infl,rho_f,rho_r,taud,
     $              rho_prev,rho_deb_rgvb,tauini,rhos,rho_mod,fac,
     $              for_rhos,rev_rhos,cts,for_rhosl,rev_rhosl,AMSS,
     $              deltaD,xp_g,drho_tot,debs,rho_ps,tautwint,
c                   out
     $              isectw,twthres,itwinlaw,dum1,dum2,dum3,dum4)
            elseif (ihardlaw.eq.20) then
               if (NPH.gt.2) then
                  write(*,*)'update_crss_dd_fe_new has not been'
                  write(*,*)'tested for multiphase polycrystal'
                  call xit()
               endif
               call update_crss_dd_fe_new(
c     intent(in)
     $              2,1,UR_SX,UW_DIS,UW_TAU,nphmx,nsysmx,nslmod,nslsys,nmodes,
     $              nsyst,nsm,ntwmod,temp,ngr,wgt,ngrmx,DBAR,dtime,
     $              nmodmx,gamdot,NTWMMX,xmu_mode,
c     intent(inout)
c          read from (ioption.eq.0) and not evolving
     $              grsze,chifor,chisub,RHO_initW,TEMPref,burg,actener,
     $              K1gener,DRAG,EDOT0,rho_inits,Atau0,Atau1,Atau2,
     $              Fdeb0,Fdeb1,Fdeb2,HPcoef,T0nucl,T1nucl,T2nucl,T0prop
     $              ,T1prop,T2prop,ctwsl,cts1,cts0,cts2,K3recov,K4recov,
     $              D0vac,Efvac,Emvac,
c         subject to 'evolution'
     $              rho_deb,rho_s_for,rho_s_deb,
c          constants
     $              edotref,
c     intent(out)
     $              itwinlaw,hard,twthres,isectw,nrs,crss,
     $              rhosat,TAU_PROP,TAU_NUCL,
c     intent(out) for stats of DD law

     $              rho_for_ph_mod,rho_deb_ph,rho_tot_ph,rho_deb_ph_mod,
     $              rho_recov,
     $              taue_wavg,rhonucl,pnucl
     $              )


c$$$               if (noel.eq.1.and.npt.eq.1) then
c$$$
c$$$
c$$$
c$$$
c$$$                  if (jstep(1).eq.1.and.kinc.eq.1) then
c$$$c                     <crss_wavg.out>
c$$$                     open(UW_TAU+1,file=OUTDIR(1:LENOUTDIR)//
c$$$     $                    '/crss_wavg.OUT',status='unknown',
c$$$     $                    action='write')
c$$$                     write(UW_TAU+1,'(a13)',advance='no')'Time'
c$$$                     do imo=1,nmodes(iphel)
c$$$                        write(UW_TAU+1,'(x,a10,i1.1,a1)',advance='no')
c$$$     $                       '<CRSS_WAVG',imo,'>'
c$$$                     enddo
c$$$                     write(UW_TAU+1,*)
c$$$
c$$$                  else
c$$$                     open(UW_TAU+1,file=OUTDIR(1:LENOUTDIR)//
c$$$     $                    '/crss_wavg.OUT',status='unknown',
c$$$     $                    position='append',action='write')
c$$$                  endif
c$$$
c$$$
c$$$                  crss_wavg(:,:)=0d0
c$$$                  DO KKK=NGR(IPHEL-1)+1,NGR(IPHEL)
c$$$                     ISM=0
c$$$                     DO KMO=1,NMODES(IPHEL)
c$$$                        DO ISY=1,NSM(KMO,IPHEL)
c$$$                           ISM=ISM+1
c$$$                           crss_wavg(kmo,iphel)=crss_wavg(kmo,iphel)
c$$$     $                          +crss(ism,kkk)*wgt(kkk)
c$$$     $                          /nsm(kmo,iphel)
c$$$                        enddo
c$$$                     enddo
c$$$                  enddo
c$$$                  write(*,*)'UW_TAU+1:',UW_TAU+1,
c$$$     $                 'nmodes(iphel):',nmodes(iphel)
c$$$                  write(UW_TAU+1,'(24e13.4)')time,
c$$$     $                 (crss_wavg(imo,iphel),IMO=1,NMODES(IPHEL)) !,
c$$$                  write(*,'(24e13.4)')time,
c$$$     $                 (crss_wavg(imo,iphel),IMO=1,NMODES(IPHEL)) !,
c$$$               endif
c$$$
c$$$               if (noel.eq.1.and.npt.eq.1) then
c$$$                  if (jstep(1).eq.1.and.kinc.eq.1) then
c$$$                     open(UW_DIS,file=OUTDIR(1:LENOUTDIR)//'/DIS_DEN.OUT',
c$$$     $                    status='unknown',
c$$$     $                    action='write')
c$$$                     open(UW_TAU,file=OUTDIR(1:LENOUTDIR)//'/taue.OUT',
c$$$     $                    status='unknown',
c$$$     $                    action='write')
c$$$
c$$$                     write(UW_DIS,'(a12)',advance='no')'Time'
c$$$                     do imo=1,nslmod(iphel)
c$$$                        write(UW_DIS,'(a11,i1.1)',advance='no')
c$$$     $                       'RHO_for_md',imo
c$$$                     enddo
c$$$                     write(UW_DIS,'(a12)',advance='no') 'RHO_deb_ph'
c$$$                     write(UW_DIS,'(a12)',advance='no') 'RHO_tot_ph'
c$$$                     do imo=1,nslmod(iphel)
c$$$                        write(UW_DIS,'(a11,i1.1)',advance='no')
c$$$     $                       'RHO_deb_md',imo
c$$$                     enddo
c$$$                     write(UW_DIS,'(a12)',advance='no') 'rho_recov'
c$$$                     write(UW_DIS,*)
c$$$
c$$$
c$$$
c$$$c                    column heads of taue.out
c$$$                     write(UW_TAU,'(a12)',advance='no')'Time'
c$$$                     do imo=1,nslmod(iphel)
c$$$                        write(UW_TAU,'(a11,i1.1)',advance='no') 'Peierls',
c$$$     $                       imo
c$$$                     enddo
c$$$! Forest
c$$$                     do imo=1,nslmod(iphel)
c$$$                        write(UW_TAU,'(a11,i1.1)',advance='no') 'Forest',
c$$$     $                       imo
c$$$                     enddo
c$$$! debris (substructure)
c$$$                     do imo=1,nslmod(iphel)
c$$$                        write(UW_TAU,'(a11,i1.1)',advance='no') 'Deb',imo
c$$$                     enddo
c$$$! Hall-Petch
c$$$                     do imo=1,nslmod(iphel)
c$$$                        write(UW_TAU,'(a11,i1.1)',advance='no') 'HP',imo
c$$$                     enddo
c$$$!! twins
c$$$                     write(UW_TAU,'(a12)',advance='no') 'RHOnucl'
c$$$! propagation
c$$$                     write(UW_TAU,'(a12)',advance='no') 'Pnucl'
c$$$                     write(UW_TAU,*)
c$$$                  else
c$$$                     open(UW_DIS,file=OUTDIR(1:LENOUTDIR)//'/DIS_DEN.OUT',
c$$$     $                    status='unknown',
c$$$     $                    position='append',action='write')
c$$$                     open(UW_TAU,file=OUTDIR(1:LENOUTDIR)//'/taue.OUT',
c$$$     $                    status='unknown',
c$$$     $                    position='append',action='write')
c$$$                  endif
c$$$
c$$$!! etc.
c$$$
c$$$                  WRITE(UW_DIS,'(14e12.3e3)') time,
c$$$     $                 (RHO_FOR_PH_mod(KMO,IPHEL),
c$$$     $                 KMO=1,NSLMOD(IPHEL))
c$$$     $                 ,RHO_DEB_PH(IPHEL),RHO_TOT_PH(IPHEL)
c$$$     $                 ,(rho_deb_ph_mod(kmo,iphel),kmo=1,nslmod(iphel))
c$$$     $                 ,rho_recov
c$$$c     $        ,RHO_FOR_RECOV,RHO_DEB_RECOV
c$$$
c$$$                  write(UW_TAU,'(28e12.3e3)') time,
c$$$     $                 ((taue_wavg(iphel,imo,i),imo=1,nslmod(iphel))
c$$$     $                 ,i=1,4),
c$$$     $                 ((taue_wavg(iphel,imo,i),imo=nslmod(iphel)+1,
c$$$     $                 nmodes(iphel)),i=1,3),
c$$$     $                 RHOnucl,Pnucl
c$$$
c$$$                  close(UW_DIS)
c$$$                  close(UW_TAU)
c$$$               endif
            endif               ! end of ihardlaw
         enddo                  ! end of iphel
      endif                     ! end of iuphard

      acctimeFE=acctimeFE+dtime
      kount_umat_call =  kount_umat_call + 1

      if (print_loc) then
         write(*,*)'end of update_evpsc_statev.f'
      endif
