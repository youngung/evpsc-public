c     A simplied initial state guess subroutine for EVPSC codes
      subroutine initial_state_guess(
c          in
     $     mdim,iregime,strain_control,interaction,
     $     iDBARv,iSBARv,isense,nph,nphmx,NSYSMX,ngrmx,
     $     ngr,nsyst,
     $     nrs,dtime,gamd0s,crss,MGEL,sch,nsch,c6,sg_old,
     $     wgt,dg_trans,selav,ioptfs,
c          inout
     $     dbar,sbar,dsbar,sbar_old,
c          out
     $     sav,dav,dbar0,tg,dgel,mgtg,dg0,mbartg,LBARTG,
     $     dsg,sg,dSGTRY,SGTRY)

      implicit none
      integer, intent(in)::mdim,nphmx,NSYSMX,
     $     iregime,strain_control,interaction,
     $     iDBARv(6),iSBARv(6),isense(NSYSMX,nphmx),nph,
     $     ngrmx,ngr(0:nphmx),nsyst(nphmx),nrs(NSYSMX,nphmx),ioptfs
      real*8, intent(in)::dtime,gamd0s(NSYSMX,ngrmx),
     $     crss(NSYSMX,ngrmx),MGEL(6,6,ngrmx),sch(5,NSYSMX,ngrmx),
     $     nsch(6,NSYSMX,ngrmx),c6,sg_old(6,ngrmx),wgt(ngrmx),
     $     dg_trans(6,ngrmx),selav(6,6)
      real*8, intent(inout)::dbar(6),sbar(6),dsbar(6),sbar_old(6)
      real*8, intent(out)::dbar0(6),tg(6,ngrmx),dgel(6,ngrmx),
     $     mgtg(6,6,ngrmx),dg0(6,ngrmx),mbartg(6,6),LBARTG(6,6),
     $     dsg(6,ngrmx),sg(6,ngrmx),dSGTRY(6,ngrmx),SGTRY(6,ngrmx),
     $     sav(6),dav(6)
c     locals
c     ------------------
      real*8 DASDBARx(mdim),SASTAVx(mdim)
      real*8 MASTPHx(6,6),MASTGRx(6,6),DG_TRANSx(6),
     $     MGELx(6,6),
     $     DSGTRYx(6),SGTRYx(6),sg_oldx(6),
     $     SCHx(5,NSYSMX),NSCHx(6,NSYSMX),
     $     GAMD0Sx(NSYSMX),
     $     CRSSx(NSYSMX),GAMDOTx(NSYSMX),
     $     DSGx(6),SGx(6),MGTGx(6,6)
      integer ISENSEx(NSYSMX),NRSx(NSYSMX),ISHAPEx,NSYSTx,INTERACTIONx
c     --
      real*8 aux6(6),cux6(6),aux5(5),aux55(5,5),
     $     dux6(6),eux6(6),fux6(6),gux6(6),
     $     hux6(6),aux66(6,6),bux6(6),bux66(6,6),iux6(6),jux6(6),
     $     dgpl(5),aux33(3,3),aux3333(3,3,3,3),bux33(3,3),
     $     dbarnorm,vnorm,sdot(6)
      integer i,j,k,
     $     iph,kkk,inverr,idiag_sg,idiag_nr
      logical isnan
      DBARNORM=VNORM(DBAR(1:mdim),MDIM)

      if (iregime.eq.2.or.iregime.eq.3) then
         dbar0(:)=0.
c        solve the linear elastic constitutive equation
         sdot(:)=sbar(:)/dtime  !! sbar-dot
         call chg_basis(dbar,aux33,aux66,aux3333,1,6)
         call chg_basis(sdot,bux33,aux66,aux3333,1,6)
         call CHG_BASIS(aux6,aux33,selav,aux3333,3,6)
         call state_6x6(IDBARv,aux33,iSBARv,bux33,aux3333)
         call CHG_BASIS(sdot,bux33,aux66,aux3333,2,6)
         call CHG_BASIS(dbar,aux33,aux66,aux3333,2,6)
         sbar(:)=sdot(:)*dtime
         dsbar(:)=sbar(:)
         sbar(:) =sbar_old(:)+dsbar(:)
      endif

c     ---------------------------------------------------------------------------
c     grain stress, strain rate, and moduli
c     ## PX
      DASDBARx(:)=0.
      SASTAVx(:)=0.
      do iph=1,nph
c        ## phase-specific
         NSYSTx=NSYST(IPH)
         ISENSEx(:NSYSTx)=isense(:NSYSTx,iph)
         NRSx(:NSYSTx)=nrs(:NSYSTx,iph)
         MASTPHx(:,:)=0.
         ISHAPEx=0
      do kkk=ngr(iph-1)+1,ngr(iph)
c        ## grain stress, grain strain rate, and moduli
         if (iregime.eq.1) then !! VPSC
            IF(STRAIN_CONTROL.EQ.1) THEN !! strain control
c     grain specific
               MASTGRx(:,:)=0.
               DG_TRANSx(:)=0.
               SG_OLDx(:)=SG_OLD(:,KKK)
c              --
               SGTRY(:,kkk)=DBAR(:)/DBARNORM
c              --
               SGTRYx(:)=SGTRY(:,kkk)
               dSGTRYx(:)=0.    ! for VPSC

               MGELx(:,:)=MGEL(:,:,kkk)
               SCHx(:,:NSYSTx)=sch(:,:NSYSTx,kkk)
               NSCHx(:,:NSYSTx)=nsch(:,:NSYSTx,kkk)
               GAMD0Sx(:NSYSTx)=gamd0s(:NSYSTx,kkk)
               CRSSx(:NSYSTx)=crss(:NSYSTx,kkk)

               INTERACTIONx=0
               call grain_stress_evp(
c                  in
     $              INTERACTIONx,mdim,
     $              IREGIME,ISHAPEx,ISENSEx,NRSx,
     $              NSYSTx,NSYSMX,

     $              DASDBARx,SASTAVx,
     $              DBAR,DBAR0,DG_TRANSx,dtime,
     $              SBAR,SBAR_OLD,SG_OLDx,

     $              MASTPHx,MASTGRx,

     $              SCHx,NSCHx,c6,CRSSx,
     $              GAMD0Sx,
     $              MGELx,
     $              dSGTRYx,SGTRYx,
c                   out
     $              dSGx,SGx,k,idiag_nr,idiag_sg)
               if (idiag_nr.ne.0.or.idiag_sg.ne.0) then
                  write(*,*)'** Error, grain stress not good'
                  write(*,*)'iph:',iph,'kkk:',kkk
                  write(*,*)'idig_nr:',idiag_nr,'k:',k
                  write(*,*)'idig_sg:',idiag_sg
                  write(*,*)'iregime:',iregime
                  write(*,*)'dbar(:):',dbar(:),'dbarnorm:',dbarnorm
                  write(*,*)'SGTRY(:,kkk):',SGTRY(:,kkk)
                  write(*,*)'dsg(:,kkk):',dsg(:,kkk)
                  write(*,*)'sg(:,kkk):',sg(:,kkk)
                  call xit()
               endif
               if (iregime.eq.1) then
                  SG(:,kkk)=SGx(:)
                  dSG(:,kkk)=SGx(:)-SG_old(:,kkk)
               elseif (iregime.eq.2) then
                  dsg(:,kkk)=dSGx(:) ! out
                  sg(:,kkk)=SGx(:) ! out
               endif

c           Ends grain_stress_evp

            ELSEIF(STRAIN_CONTROL.EQ.0) THEN !! stress control
               write(*,*)"** Need to here do Sach's guess"
               call xit()
            endif

         elseif(iregime.eq.2) then !! EVPSC
c           Sachs assumption - uniform stress.
            SGTRY(:,kkk)=sbar(:)
            dSGTRY(:,kkk)=dsbar(:)
            sg(:,kkk) =SGTRY(:,kkk)
            dsg(:,kkk)=dSGTRY(:,kkk)
         endif                  ! iregime

c     *** for each grain calculates DG (using power law), MGTG and DG0
         GAMD0Sx(1:NSYSTx)=gamd0s(1:NSYSTx,kkk)
         MGELx(:,:)=MGEL(:,:,kkk)
         SG_OLDx(:)=sg_old(:,kkk)
         dSGx(:)=dsg(:,kkk)
         SGx(:)=sg(:,kkk)
         DG_TRANSx(:)=dg_trans(:,kkk)
         SCHx(:,1:NSYSTx)=sch(:,1:NSYSTx,kkk)
         NSCHx(:,1:NSYSTx)=nsch(:,1:NSYSTx,kkk)
         CRSSx(1:NSYSTx)=crss(1:NSYSTx,kkk)

c        below two lines are phase-specific
         NRSx(1:NSYSTx)=nrs(1:NSYSTx,iph)
         ISENSEx(1:NSYSTx)=isense(1:NSYSTx,iph)

         call grain_rate_and_moduli_evp(
c     intent(in)
     $        interaction,iregime,ioptfs,NRSx,GAMD0Sx,
     $        CRSSx,MGELx,
     $        ISENSEx,SCHx,NSCHx,c6,SG_OLDx,
     $        dSGx,
     $        SGx,
     $        NSYSTx,NSYSMX,dtime,DG_TRANSx,
c     intent(out)
     $        aux6,aux5,cux6,GAMDOTx,MGTGx,dux6
     $        )
c        gamdot(1:NSYSTx,kkk)=GAMDOTx(1:NSYSTx) ** gamdot is saved in <evpsc>
         mgtg(:,:,kkk)=MGTGx(:,:)
         dgel(:,kkk)=aux6(:)
         dgpl(:)=aux5(:)
         tg(:,kkk)=cux6(:)
         dg0(:,kkk)=dux6(:)

      enddo                     ! grain
      enddo                     ! phase


      DO I=1,MDIM
        DAV(I)=0.
        SAV(I)=0.
        DO IPH=1,nph
        DO KKK=NGR(IPH-1)+1,NGR(IPH)
          DAV(I)=DAV(I)+tG(I,KKK)*WGT(KKK)
          SAV(I)=SAV(I)+SG(I,KKK)*WGT(KKK)
        ENDDO
        ENDDO
      ENDDO


C     CALCULATE INITIAL GUESS FOR MACROSCOPIC MODULI 'Mtg' AS THE INVERSE
C     OF THE AVERAGE OF THE GRAIN'S STIFFNESSES --> MBARTG=<LBARTGc>^-1
C     ACTUALLY, WHEN INTERACTION=2,3,4 IT AVERAGES THE SECANT MODULI

      dbar0(:)=0.
      LBARTG(:,:)=0.
      MBARTG(:,:)=0.

      if (iregime.eq.1) then !!VPSC
         do iph=1,nph
         do kkk=ngr(iph-1)+1,ngr(iph)
            AUX55(1:mdim,1:mdim)=MGTG(1:mdim,1:mdim,KKK) ! MGTG=MGSEC if INTERACTION=0,2,3,4
            CALL LU_INVERSE(inverr,AUX55,MDIM) ! (yj): it is not a great practice to create a temporary array
            if (inverr.ne.0) then
               write(*,*)'could not properly inverse aux66 in init'
               call xit()
            endif
            DO I=1,MDIM
               DBAR0(I)=DBAR0(I)+(DG0(I,KKK)+DG_TRANS(I,KKK))*WGT(KKK)
            DO J=1,MDIM
               LBARTG(I,J)=LBARTG(I,J)+AUX55(I,J)*WGT(KKK)
            ENDDO
            ENDDO
         enddo
         enddo
         aux55(:,:)=lbartg(1:mdim,1:mdim)
         call lu_inverse(inverr,aux55,mdim)
         mbartg(1:mdim,1:mdim)=aux55(:,:)
      elseif (iregime.ge.2) then !!
c     Assume mbartg=<mgtg>, and obtain LBARTG=mbartg^-1
         do iph=1,nph
         do kkk=ngr(iph-1)+1,ngr(iph)
            mbartg(:,:)=mbartg(:,:)+mgtg(:,:,kkk)*wgt(kkk)
         enddo
         enddo
         LBARTG(:,:)=mbartg(:,:)
         call lu_inverse(inverr,LBARTG,6)
      endif
      return
      end subroutine initial_state_guess
