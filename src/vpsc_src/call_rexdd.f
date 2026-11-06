         if (iregime.eq.1) then
            do kkk=ngr(iph-1)+1,ngr(iph)
               sg(6,kkk)=0.
            do j=1,6
               sg(6,kkk)=sg(6,kkk)+bc_el(6,j,kkk)*sbar(j)
            enddo
            enddo
         endif
         call rexdd(
c             intent(in):
c             * integer
     $        iopt,iregime,mdim,iproc,istep,nphmx,nph,iph,nclustermx,
     $        nchildmx,nchildlimit,ngrmx,ngidmx,nsysmx,nmodmx,ntwmmx,
     $        nmodes(iph),ntwmod(iph),nslmod,nsm,nslsys,uw_rex,
     $        UW_REX_STAT,uw_rex_nuco,uw_rex_grow,nsyst(iph),
     $        isense(:,iph),nrs(:,iph),
c             * real
     $        schca(:,:,iph),nschca(:,:,iph),gamd0s,time,dtime,temp,
     $        rexgrowpar,rexnuclpar,sbar,dsbar,dbar,dbar0,epsvm,svm,
     $        celcc(:,:,iph),RHO_initS,RHO_initW,xmu_mode,
     $        atau0,atau1,atau2,chifor,chisub,burg,HPcoef,
     $        ctwsl,RHOSAT,tau_prop,tau_nucl,
c             intent(inout):
c             * integer
     $        ncluster,nchildset,ngr,dngr_rex,
     $        ngid_total,
     $        gid_of_slot,slot_of_gid,
     $        gid_birth_step,gid_death_step,gid_parent,
c             * real
     $        wgt,wgt_old,ag,sg,dsg,sgtry,dsgtry,sg_old,fijgr,axisgr,
     $        dg_trans,gtotgr,mgel,mgtg,dg0,radius,rho_s_for,rho_s_deb,
     $        crss,
c             intent(out):
c             * integer
     $        ngref,igref,nkount,idclusterset,idclusterset_reverse,
c             * real
     $        vpx,rexdd_tiny,dwgt_grow,dwgt_nucl,
c             * logical
     $        inucleated_px,inucleated_cl)
         write(*,*)'=============================='
