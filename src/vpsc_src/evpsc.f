C********************************************************************
C     SUBROUTINE EVPSC
c     - The numerical code based on the delta-EVPSC & sigma-EVPSC
c
c        ** delta-EVPSC code
c          published in Jeong & Tome (2020) IJP
c          published in Jeong, Jeon, Tome (2021) IJP
c        ** enchanced delta-EVPSC and sigma-EVPSC models
c          published in B. Jeon, Y. Jeong, C. N. Tome (2024), JMRT
c
c     Y. Jeong @ Changwon National University, email: yjeong@changwon.ac.kr
c
c      -----------------------------------------
c      |    (Outer loop, IT2)                  |
c      |   * Update Macro Quantities           |
c      |      - MASTph or MASTgr               |
c      |      - MBARTG, dbar0                  |
c      |  ------------------------------------ |
c      |  |  (inner loop, IT1)               | |
c      |  |  * SOLVE Boundary Condition      | |
c      |  |     - Obtain DBAR, dSBAR, SBAR   | |
c      |  |  * Obtain grain stress           | |
c      |  |     - SOLVE NEWTON-RAPHSON       | |
c      |  |     - Obtain dSG (or sg)         | |
c      |  |  * Obtain grain rate             | |
c      |  |     - TG                         | |
c      |  |     - mgtg, dg0                  | |
c      |  -----------------------------------| |
c      |   * Check consistency condition       |
c      -----------------------------------------
C********************************************************************

c     Argument order: Numerical CONTROL & CONDITION, DIMENSION, PX, PHASE, GRAIN, SYSTEM
c     Then, if applicable, alphabetical order is used.
      subroutine evpsc(
c          intent(in)
c             ** logical in
     $     print_loc,
c             ** integers in
     $     mdim,ife,iregime,interaction,iverbose,ioptfs,
     $     itmaxint,itmaxext,UW_LOG,UW_SCCN,istep,
     $     nphmx,ngrmx,NSYSMX,iDBARv,iSBARv,
     $     nph,ngr,nsyst,ishape,isense,nrs,

c             ** floats in
     $     errs,errd,errm,dtime,fstar,alplpx,c6,
     $     sbar_old,axisph,axisgr,wgt,sg_old,dg_trans,mgel,
     $     crss,gamd0s,sch,nsch,

c          intent(inout)
c             ** floats inout
     $     MBARTG,LBARTG,dbar0,
     $     sbar,dsbar,sbarc,dsbarc,dbar,
     $     sgtry,dsgtry,sg,dsg,dg0,mgtg,

c          intent(out)
c             ** integers out
     $     it2,idiag_sc,idiag_sg,idiag_nr,
c             ** floats out
     $     xnumstats,dbarc,
     $     factor,fel,fvp,alpl,
     $     dge,dgp,dg,tg,gamdot,asph,asgr,

c     I am not sure if these are required to be passed out as an argument
     $     sav,dav
     $     )
      implicit none
c             ** logical in
      logical,intent(in):: print_loc
c             ** integers in
      integer,intent(in):: mdim,iregime,ife,interaction,iverbose,ioptfs,
     $     itmaxint,itmaxext,UW_LOG,UW_SCCN,istep,
     $     nphmx,ngrmx,NSYSMX,
     $     iDBARv(6),iSBARv(6),
     $     nph,ngr(0:nphmx),nsyst(nphmx),ishape(0:nphmx),
     $     isense(NSYSMX,nphmx),nrs(NSYSMX,nphmx)
c             ** floats in
      real*8, intent(in)::errs,errd,errm,dtime,fstar,alplpx,c6,
     $     sbar_old(6),
     $     axisph(0:3,3,0:nphmx),axisgr(0:3,3,ngrmx),
     $     wgt(ngrmx),sg_old(6,ngrmx),dg_trans(6,ngrmx),mgel(6,6,ngrmx),
     $     crss(NSYSMX,ngrmx),gamd0s(NSYSMX,ngrmx),
     $     sch(5,NSYSMX,ngrmx),nsch(6,NSYSMX,ngrmx)
c             ** floats inout
      real*8, intent(inout)::
     $     MBARTG(6,6),LBARTG(6,6),dbar0(6),
     $     sbar(6),dsbar(6),dsbarc(3,3),sbarc(3,3),dbar(6),
     $     sgtry(6,ngrmx),dsgtry(6,ngrmx),sg(6,ngrmx),dsg(6,ngrmx),
     $     dg0(6,ngrmx),mgtg(6,6,ngrmx)
c             ** integers out
      integer, intent(out):: it2,idiag_sc,idiag_sg,idiag_nr
c             ** floats out
      real*8, intent(out)::xnumstats(7),dbarc(3,3),

     $     factor(ngrmx),fel(ngrmx),fvp(ngrmx),
     $     alpl(ngrmx),
     $     dge(6,ngrmx),dgp(5,ngrmx),
     $     dg(6,ngrmx),tg(6,ngrmx),
     $
     $     gamdot(NSYSMX,ngrmx),

     $     asph(3,3,3,3,nphmx),
     $     asgr(3,3,3,3,ngrmx)


c     ----------------------------------
c     locals
      logical iesh_fail
      integer iph,kkk,i,j,k,inverr_upmacro,iDBARb(6),iSBARb(6)
      real*8 tnorm6x6,chi,t0,t1,elc,vpc,relsgr,rels,reld,
     $     sastmix,SASTAVX,SASDBARX,
     $     rels6(6),reld6(6),dsav(6),sav(6),dav(6),
     $     aux5(5),aux6(6),cux6(6),dux6(6),aux66(6,6),aux55(5,5),
     $     aux33(3,3),aux3333(3,3,3,3),
     $     MASTph(mdim,mdim,nphmx),MASTgr(mdim,mdim,ngrmx),
     $
     $     MBARTG_old(6,6),dbar0_old(6),DBAUX(mdim),DBAUXc(3,3),
     $     sastav(mdim),SASDBAR(6),
     $     DASDBAR(mdim),sg0av(mdim),bginv(mdim,mdim),
     $     auxm(mdim),auxmm(mdim,mdim)
      integer iauxm1(mdim),iauxm2(mdim)
c     ------------------------------------------------------------------
c     Arrays to copy temporary arrays / values for grain stress & grain rate
      real*8 MASTPHx(mdim,mdim),MASTGRx(mdim,mdim),DG_TRANSx(6),
     $     MGELx(6,6),
     $     DSGTRYx(6),SGTRYx(6),SG_OLDx(6),
     $     SCHx(5,NSYSMX),NSCHx(6,NSYSMX),
     $     GAMD0Sx(NSYSMX),
     $     CRSSx(NSYSMX),GAMDOTx(NSYSMX),
     $     DSGx(6),SGx(6),MGTGx(6,6)
      integer ISENSEx(NSYSMX),NRSx(NSYSMX),ISHAPEx,NSYSTx
c     ------------------------------------------------------------------
      real*8 xnriters(ngrmx),vnorm
      logical verbose
c      parameter(verbose=.true.)
      parameter(verbose=.false.)

c     ------------------------------------------------------------------
c     statements to fool the compiler (not used but can be used sometimes)
      i=istep
      i=UW_LOG
c     ------------------------------------------------------------------
      idiag_sc=0

      if (ife.eq.0) call cpu_time(t0)
c     ==== outer-loop over stress and compliances in grains

c     Polycrystal properties pertaining to the previous step
c    * MBARTG: Overall elasto-visco-plastic compliance tensor
c    * dbar0: Overall back-extrapolated strain rate tensor
      MBARTG_old(:,:)=MBARTG(:,:)
      dbar0_old(:)=dbar0(:)

c      xnumstats(:)=0.
      relsgr=errs*2d0
      rels=errs*2d0
      reld=errd*2d0

      if (verbose) write(*,*)'** beginning of outer loop'
c     Iteration counter and Relative errors set by two times the given tolerances.
      it2=0

c     Loop over phases within this integration point
      if (iverbose.ge.1) then
         write(*,'(2x,a2,x,3a10,x,4a8,x,a3)')'it2','relsgr','rels',
     $        'reld','nritmn','nritmx','nritav','failfrc','irg'
      endif

      do 1000 while((
     $     relsgr.gt.errd
     $     .or.rels.gt.errs
     $     .or.reld.gt.errd)
     $     .and.(it2.lt.itmaxext))
         it2=it2+1
c     ==== inner-loop over overall compliance.

c-----------------------------------------------------------------------
c     *** f-transition from 1 to fstar ***
         factor(:)=1.
         if (iregime.eq.2.OR.iregime.eq.3) then ! introduced for dEVPSC and sEVPSC formulations
c            if (it2.eq.1) write(*,'(a10,f8.4)')'alplpx:',alplpx
c            if (.false.) then
            if (.true.) then
               do iph=1,nph
                  do kkk=ngr(iph-1)+1,ngr(iph)
c                     chi = vnorm(tg(:,kkk),6)/vnorm(dge(:,kkk),5)
                     MGELx(:,:)=mgel(:,:,kkk)
                     chi=tnorm6x6(mgtg(:,:,kkk))/tnorm6x6(MGELx)
                     if (chi.lt.0) chi=1.
                     factor(kkk)=fstar-(fstar-1)*exp(1-chi)
                  enddo
               enddo
            else
               if (alplpx.gt.0.1) then
                  factor(:)=fstar
               else
                  factor(:)=1.
               endif
            endif
         endif

c     ** Calculate
c        calculate Eshelby tensor S and FS=(I-S)^(-1) S
         if (print_loc) write(*,*)"** calling update_macro_property"
         call update_macro_property(
c             in
     $        mdim,iregime,errm,itmaxint,nphmx,ngrmx,
     $        nph,ngr,axisph,axisgr,
     $        ishape,dg0,mgtg,wgt,factor,dg_trans,
c             inout
     $        LBARTG,MBARTG,dbar0,
c             out
     $        sg0av,bginv,asph,asgr,MASTph,MASTgr,iesh_fail,
     $        inverr_upmacro)

         if (inverr_upmacro.ne.0)then
            idiag_sc=2
            return
         elseif(iesh_fail) then
            write(*,*)'** eshelby fail with it2:',it2
            write(*,*)'mbartg'
            do i=1,mdim
               write(*,'(6e11.3)')(mbartg(i,j),j=1,mdim)
            enddo
            write(*,*)'lbartg'
            do i=1,mdim
               write(*,'(6e11.3)')(lbartg(i,j),j=1,mdim)
            enddo
            write(*,*)'dbar0'
            write(*,'(6e11.3)')(dbar0(i),i=1,mdim)
            call xit()
            idiag_sc=3
            return
         endif

C ************************************************************************************
C     SOLVES (DBAR-DBAR0)=MBAR * SBAR ACCOUNTING FOR MIXED BOUNDARY CONDITIONS
C     IN STRESS AND STRAIN RATE COMPONENTS.

C     IF MBAR IS VISCO-PLASTIC INCOMPRESSIBILITY DICTATES THAT EITHER 3 OR 1
C     DIAGONAL COMPONENTS OF D(I,J) CAN BE ENFORCED (0 OR 2 DIAGONAL STRESS COMPS)
C     (SINCE D1~(d22-d11) & D2~(d33-d11)+(d33-d22) THEN EITHER BOTH OR NONE OF THE
C     FIRST TWO b-BASIS COMPONENTS ARE ENFORCED, RESPECTIVELY).
C     --> MBAR IS 5x5 DEVIATORIC, DBAR IS 5X1 DEVIATORIC
C     --> IF 3 DIAGONAL COMPS OF DBAR ARE ENFORCED ONLY THE DEVIATORIC PART OF SBAR
C         CAN BE OBTAINED. THE CONSTITUTIVE EQUATION CAN BE SOLVED IN THE b-BASIS
C         REPRESENTATION AS A 5x5 SYSTEM
C     --> IF ONE DIAGONAL COMPONENT OF SBAR IS ENFORCED THE OVERALL PRESSURE CAN BE
C         CALCULATED (ALTHOUGH THE BEHAVIOR IS STILL INCOMPRESSIBLE AND INDEPENDENT
C         OF THE PRESSURE).
C         THE FLAGS ON THE 'CARTESIAN COMPONENT' CONTROLING BOUNDARY CONDITIONS
C         CANNOT BE WRITTEN IN TERMS OF b-BASIS FLAGS AND THE SYSTEM HAS TO BE SOLVED
C         IN THE CARTESIAN REPRESENTATION WRITEN AS A 6x6 SYSTEM IN VOIGT NOTATION

C     IF MBAR IS ELASTO-VISCO-PLASTIC: MBAR, SBAR & DBAR ARE 6X6 & 6X1 NON-DEVIAT.
C     THE CONSTITUTIVE EQUATION IS 6x6 IN b-BASIS REPRESENTATION.
C     FBECAUSE THE CARTESIAN FLAGS ON IMPOSED BC's CANNOT BE TRANSLATED INTO b-BASIS
C     FLAGS THE PROBLEM HAS TO BE SOLVED AS A 6x6 EQUATION IN TERMS OF VOIGT COMPS.
C ************************************************************************************

C *** SUBSTRACT INDEPENDENT TERM, SOLVE LINEAR EQS, ADD IT AFTERWARDS
         DBAUX(1:mdim)= DBAR(1:mdim)-DBAR0(1:mdim)
         IF (IDBARv(1)*IDBARv(2)*IDBARv(3).EQ.1 .AND. iregime.eq.1) THEN

            IDBARb(1)=1
            IDBARb(2)=1
            IDBARb(6)=1 ! dilation is set to zero, thus is known.
            IDBARb(3)=IDBARv(4)
            IDBARb(4)=IDBARv(5)
            IDBARb(5)=IDBARv(6)
            ISBARb(1)=0
            ISBARb(2)=0
            ISBARb(6)=0 ! hydrostatic pressure is not known
            ISBARb(3)=iSBARv(4)
            ISBARb(4)=iSBARv(5)
            ISBARb(5)=iSBARv(6)

            iauxm1(:)=iDBARb(1:mdim)
            iauxm2(:)=iSBARb(1:mdim)
            auxm(:)=sbar(1:mdim)
            auxmm(:,:)=MBARTG(1:mdim,1:mdim)
            CALL STATE_NxN(iauxm1,DBAUX,iauxm2,auxm,auxmm,mdim)
            SBAR(1:mdim)=auxm(:)
            DBAR(1:mdim)=DBAUX(1:mdim)+DBAR0(1:mdim)
            CALL CHG_BASIS (DBAR,DBARc  ,AUX55,AUX3333,1,5)
            CALL CHG_BASIS (SBAR,SBARC  ,AUX55,AUX3333,1,5)
            SG(6,:)=0.
            DG(6,:)=0.

         ELSE IF(IDBARv(1)*IDBARv(2)*IDBARv(3).EQ.0.OR.iregime.ge.2)THEN
            CALL CHG_BASIS (DBAUX,DBAUXc,AUX66,AUX3333,1, MDIM)
            auxmm(:,:)=MBARTG(1:mdim,1:mdim)
            CALL CHG_BASIS (AUXM,AUX33,auxmm,AUX3333,3, MDIM)

            if (iregime.eq.1) then
               CALL STATE_6x6 (IDBARv,DBAUXc,iSBARv,SBARc,AUX3333)
               CALL CHG_BASIS (SBAR, SBARc, AUX66,AUX3333,2,6)
            elseif (iregime.eq.2) then
               CALL STATE_6x6 (IDBARv,DBAUXc,iSBARv,dSBARc,AUX3333)
               CALL CHG_BASIS (dSBAR,dSBARc, AUX66,AUX3333,2,6)
            endif

            if(iregime.eq.1) then
               dsbar(:)=sbar(:)-sbar_old(:) ! VPSC
            elseif (iregime.eq.2) then
               sbar(:)=dsbar(:)+sbar_old(:) ! dEVPSC, sbar_old is zero for VPSC.
            endif


c            CALL CHG_BASIS (DBAUX,DBAUXc,AUX66,AUX3333,2,6)
            CALL CHG_BASIS (DBAUX,DBAUXc,AUXmm,AUX3333,2,mdim)
            DBAR(1:mdim)=DBAUX(1:mdim)+DBAR0(1:mdim) ! DBAR is deviatoric when MDIM=5
            CALL CHG_BASIS (DBAR(1:mdim),DBARc,AUX66,AUX3333,1,MDIM)

c$$$            write(*,*)'dbarc and sbarc'
c$$$            do i=1,3
c$$$               write(*,'(2(3e11.3,5x))')(dbarc(i,j),j=1,3),
c$$$     $              (sbarc(i,j),j=1,3)
c$$$            enddo


C *** THE FOLLOWING IS A CONSISTENCY CHECK ON MACRO PRESSURE
c            CALL CHG_BASIS (SBAR, AUX33, AUX55,AUX3333,1,5) ! deviatoric comps of SBAR
c$$$         PBAR1=0.
c$$$         PBAR2=0.
c$$$         DO I=1,3
c$$$            PBAR1=PBAR1+iSBARv(I)*(SBARc(I,I)-AUX33(I,I))
c$$$            PBAR2=PBAR2+SBARc(I,I)
c$$$         ENDDO
c$$$         PBAR1=PBAR1 / (iSBARv(1)+iSBARv(2)+iSBARv(3))
c$$$         PBAR2=PBAR2 /3.

C      write(UW_LOG,'('' CONSISTENCY CHECK ON MACRO PRESSURE'')')
C      write(UW_LOG,'('' SBARc ='',9f12.4)') SBARc
C      write(UW_LOG,'('' PBAR1  &  PBAR2 ='',9f12.4)') PBAR1,PBAR2
C      write(UW_LOG,'('' DBARc ='',9f12.4)') DBARc

C *** CALCULATES HYDROSTATIC STRESS IN GRAINS USING THE ELASTIC LOCALIZATION
C      EQUATION AND THE 'RECONSTRUCTED' MACROSCOPIC CAUCHY STRESS.
C *** ASSUMES THAT THE DEVIATORIC STRESS COMPONENTS GIVEN BY THE VISCO-
C      PLASTIC FORMALISM ARE CLOSE TO THE ONES THAT AN ELASTO-VISCO-PLASTIC
C      FORMALISM WOULD PREDICT.

            IF(MDIM.EQ.5.and.iregime.eq.1) THEN ! applicable only to VPSC.

               DO KKK=1,NGR(NPH)
                  SG(6,KKK)=0.
                  DG(6,KKK)=0.
c                  DO J=1,6
c                     SG(6,KKK)=SG(6,KKK)+BG_EL(6,J,KKK)*SBAR(J)
c                  ENDDO
               ENDDO
            ENDIF

         endif




C *************************************************************************
C     S* & D* REQUIRED FOR DIFFERENT GRAIN SHAPES
C     SASTAV(I) IS USED INSIDE GRAIN_STRESS
C     SASDBAR(I) IS USED BELOW TO CALCULATE DASDBAR(I)

         SASTMIX=0.
         IF (INTERACTION.NE.3.AND.INTERACTION.NE.4) SASTMIX=0.75 ! why neff & tg only?

         if (iregime.eq.1) then
            DO I=1,MDIM
               SASTAV(I) =0.
               SASDBAR(I)=0.
               DO J=1,MDIM
                  SASTAV(I) =SASTAV(I) +BGINV(I,J)*(SAV(J) -SG0AV(J)) ! Eq 5-34a
                  SASDBAR(I)=SASDBAR(I)+BGINV(I,J)*(SBAR(J)-SG0AV(J))
               ENDDO
               SASTAVX  =SASTAV(I)
               SASDBARX=SASDBAR(I)
               SASTAV(I)  =(1.-SASTMIX)*SASTAVX+    SASTMIX*SASDBARX
               SASDBAR(I)=     SASTMIX*SASTAVX+(1.-SASTMIX)*SASDBARX
            ENDDO

            DO I=1,MDIM
               DASDBAR(I)=0.
               DO J=1,MDIM
                  DASDBAR(I)=DASDBAR(I)+MBARTG(I,J)*SASDBAR(J) ! DBAR0 added to DBAR or DASDBAR inside N-R
               ENDDO
            ENDDO

         elseif(iregime.ge.2) then ! (YJ) for dEVPSC, the above trick (within if mdim.eq.5 condition) is not used.
            DASDBAR(:)=0
            do i=1,mdim
               if (iregime.eq.2) then
                  SASTAV(i)=(dsbar(i)-sbar_old(i))
               elseif (iregime.eq.3) then
                  SASTAV(i)=(sbar(i)-sbar_old(i))
               endif
               do j=1,mdim
                  if (iregime.eq.2) then
                     DASDBAR(i)=DASDBAR(i)
     $                    +MBARTG(I,J)*(dSBAR(J)-SBAR_OLD(J))
                  elseif (iregime.eq.3) then
                     DASDBAR(i)=DASDBAR(i)
     $                    +MBARTG(I,J)*(SBAR(J)-SBAR_OLD(J))
                  endif
               enddo
            enddo
         endif


C *** SAVE CURRENT STRESS IN GRAIN AS 'SGTRY'.
C *** STARTING FROM CONVERGED 'MBARTG' AND 'SGTRY' RECALCULATES 'SG' USING
C     THE INTERACTION EQUATION.
C *** USES NEW 'SG' TO CALCULATE STRAIN-RATE 'DG', SHEAR RATES 'GAMDOT',
C     AND MODULI 'MGTG' FOR EVERY GRAIN.


c        Loop over phases within this integration point
         xnriters(:)=0.         !! count the iterations used in NR procedure.
         relsgr=-1.
         rels=-1.
         reld=-1.

c     ---------------------------------------------------------------------------
c     grain stress, strain rate, and moduli
c        ## PX
c        DASDBAR
c        SASTAV
         do iph=1,nph
c           ## phase-specific
            NSYSTx=NSYST(IPH)
            ISENSEx(:NSYSTx)=isense(:NSYSTx,iph)
            NRSx(:NSYSTx)=nrs(:NSYSTx,iph)
            MASTPHx(:,:)=mastph(:,:,iph)
            ISHAPEx=ISHAPE(iph)
            do kkk=ngr(iph-1)+1,ngr(iph)
c              ## grain stress
               MASTGRx(:,:)=mastgr(:,:,kkk)
               DG_TRANSx(:)=dg_trans(:,kkk)
               SG_OLDx(:)=sg_old(:,kkk)
c              --
               sgtry(:,kkk)=sg(:,kkk)
               dsgtry(:,kkk)=dsg(:,kkk)
c              --
               SGTRYx(:)=sgtry(:,kkk)
               DSGTRYx(:)=dsgtry(:,kkk)

               MGELx(:,:)=mgel(:,:,kkk) ! obtained from elsc_evp_fe
               SCHx(:,:NSYSTx)=sch(:,:NSYSTx,kkk) ! updated from update_schmid
               NSCHx(:,:NSYSTx)=nsch(:,:NSYSTx,kkk)
               GAMD0Sx(:NSYSTx)=gamd0s(:NSYSTx,kkk)
               CRSSx(:NSYSTx)=crss(:NSYSTx,kkk)

c               write(*,*)'** calling grain stress for kkk:',kkk
               call grain_stress_evp(
c                  in
     $              INTERACTION,mdim,
     $              IREGIME,ISHAPEx,ISENSEx,NRSx,
     $              NSYSTx,NSYSMX,


     $              DASDBAR,SASTAV,
     $              DBAR,DBAR0,DG_TRANSx,dtime,
     $              SBAR,SBAR_OLD,SG_OLDx,
     $
     $              MASTPHx,MASTGRx,

     $              SCHx,NSCHx,c6,CRSSx,
     $              GAMD0Sx,
     $              MGELx,
     $              DSGTRYx,SGTRYx,
c                   out
     $              dSGx,SGx,k,idiag_nr,idiag_sg)
c                   k being nriter; l being idiag_nr

               xnriters(kkk)=k*1.0
               if (idiag_sg.gt.0.or.idiag_nr.gt.0.and.ife.eq.0) then
c$$$                  write(*,'(3(a15,i3))')'idiag_sg:',idiag_sg,
c$$$     $                 'idiag_nr:',idiag_nr,'kkk:',kkk
                  xnumstats(7)=xnumstats(7)+1
c                  call xit()
                  return
               endif
               if (iregime.eq.1) then
                  SG(:,kkk)=SGx(:)
                  dSG(:,kkk)=SGx(:)-SG_old(:,kkk)
               elseif (iregime.eq.2) then
                  dsg(:,kkk)=dSGx(:) ! out
                  sg(:,kkk)=SGx(:) ! out
               endif
c           Ends grain_stress_evp

c              ## grain strain rate, and moduli
               call grain_rate_and_moduli_evp(
c                   intent(in)
     $              interaction,iregime,ioptfs,NRSx,GAMD0Sx,
     $              CRSSx,MGELx,
     $              ISENSEx,SCHx,NSCHx,c6,SG_OLDx,
     $              dSGx,SGx,
     $              NSYSTx,NSYSMX,dtime,DG_TRANSx,
c                   intent(out)
     $              aux6,aux5,cux6,GAMDOTx,MGTGx,dux6
     $              )
               gamdot(1:NSYSTx,kkk)=GAMDOTx(1:NSYSTx)
               mgtg(:,:,kkk)=MGTGx(:,:)
               dge(:,kkk)=aux6(:)
               dgp(:,kkk)=aux5(:)
               tg(:,kkk)=cux6(:)
               dg0(:,kkk)=dux6(:)
c
               vpc=vnorm(dgp,5)
               elc=vnorm(dge,6)
               fvp(kkk)=vpc/(vpc+elc)
               fel(kkk)=1.-fvp(kkk)
               alpl(kkk)=vpc/(vpc+elc)
            enddo               ! loop over grains
         enddo                  ! loop over phases


*     updates macro quantities and check the consistency in solution.
         if (print_loc) write(*,*)"** calling calc_rel_err"
         call calc_rel_err(
c             in
     $        iregime,nphmx,ngrmx,nph,ngr,dsbar,sbar,dbar,
     $        tg,dsg,sg,sg_old,dsgtry,sgtry,wgt,
c             inout
     $        sbarc,dbarc,
c             out
     $        relsgr,rels,reld,rels6,reld6,sav,dsav,dav)
         if (print_loc) write(*,*)"** updates xnumstats"

         xnumstats(1)=relsgr
         xnumstats(2)=rels
         xnumstats(3)=reld

         xnumstats(4)=1000000.  ! NR iter. min
         xnumstats(5)=-10.      ! NR iter. max
         xnumstats(6)=0.        ! NR iteration counter average per the last it2 iteration
         do iph=1,nph
         do kkk=ngr(iph-1)+1,ngr(iph)
            if (xnumstats(4).gt.xnriters(kkk))xnumstats(4)=xnriters(kkk)
            if (xnumstats(5).lt.xnriters(kkk))xnumstats(5)=xnriters(kkk)
            xnumstats(6)=xnumstats(6)+xnriters(kkk)
         enddo
         enddo
         xnumstats(6)=xnumstats(6)/ngr(nph)
         xnumstats(7)=xnumstats(7)/ngr(nph) ! total NR iteration average

 1000 continue
      if (print_loc) write(*,*)"** end of the outer loop in EVPSC"

      if (ife.eq.0) then        !! stand-alone
         if (print_loc) write(*,*)"** calling cpu_time"
         call cpu_time(t1)
         if (print_loc) write(*,*)"** write to UW_SCCN"
         write(UW_SCCN,'(i3.2,x,4e12.4,2(6e12.4),e12.4)')
     $        IT2,relsgr,reld,rels,t1-t0,
     $        (reld6(i),i=1,mdim),(rels6(i),i=1,mdim),alplpx
         write(*,'(x,i2.2,x,3e10.3,x,4f8.2,x,i3)')
     $        it2,(xnumstats(i),i=1,7),iregime
      endif

      if (print_loc) write(*,*)"** checking itmaxext"

      if (it2.ge.itmaxext) then
         if (ife.eq.0) then ! stand-alone
            idiag_sc=1
            write(*,*)'** Error: exceeded max external iteration'
c            call xit()
         elseif(ife.eq.1) then
            idiag_sc=1
            call mutexlock(2)
            write(*,*)'** Warning: exceeded max external iteration'
            write(*,'(2x,a2,x,3a10,x,4a8)')'it2','relsgr','rels',
     $           'reld','nritmn','nritmx','nritav','failfrc'
            write(*,'(x,i2.2,x,3e10.3,x,4f8.2)')
     $           it2,(xnumstats(i),i=1,7)
            write(*,*)'dbar'
            write(*,'(6e11.3)')(dbar(i),i=1,mdim)
            write(*,*)'dav'
            write(*,'(6e11.3)')(dav(i),i=1,mdim)
            write(*,*)'sbar'
            write(*,'(6e11.3)')(sbar(i),i=1,mdim)
            write(*,*)'sav'
            write(*,'(6e11.3)')(sav(i),i=1,mdim)
            call mutexunlock(2)
         endif
      endif
      if (print_loc) write(*,*)"** just before return statement"
      return
      end subroutine evpsc
