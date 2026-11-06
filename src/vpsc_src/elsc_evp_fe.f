C ************************************************************************
c     SUBROUTINE ELSC_EVP   -->      VERSION 17/JUL/2018 (YJ ver. for EVPSC)
C
C     CALCULATES ITERATIVELY SELF-CONSISTENT ELASTIC MODULI OF POLYCRYSTAL
C         CELAV=<Cg:Ag>:<Ag>^(-1)
C         with Ag =<(Cg+Cg~)^(-1)>:(CELAV+Cg~)
C         with Cg~=CELAV:(Eshg^(-1)-I)
C     CALCULATES SELF-CONSISTENT THERMAL MODULI OF POLYCRYSTAL
C         ATHAV=CELAV^(-1)::<Ag^(-1)>^(-1):=<Ag^(-1):Cg:ATHg>
C
C     ITERATIVE PROCEDURE. USES b-BASIS REPRESENTATION OF TENSORS.
C       IOPTION=0 : CALCULATES UPPER BOUND ELASTIC AND THERMAL TENSORS
C                  'CELUB' AND 'ATHUB'
C       IOPTION=1 : CALCULATES SELF CONSISTENT ELASTIC AND THERMAL TENSORS
C                  'CELAV' AND 'ATHAV'
C       IOPTION=2 : CALCULATES STRESS IN GRAINS INDUCED BY ELASTIC LOAD
C                   AND/OR TEMPERATURE CHANGE
C ************************************************************************

      SUBROUTINE ELSC_evp_fe(
c          in
     $     IOPTION,nphmx,nph,athcc,celcc,selcc,
     $     ngr_whole,ngrmx,ag,wgt,itmaxint,ishape,axisph,
     $     axisgr,
c          inout
     $     celav,ac_el,bc_el,
c          out
     $     athcs,celcs,mgel,selav
     $     )
      implicit none
      integer, intent(in)::ioption,nphmx,nph,ngrmx,itmaxint,
     $     ishape(0:nphmx),ngr_whole(0:nphmx)
      real*8,intent(in)::athcc(6,nphmx),celcc(6,6,nphmx),
     $     selcc(6,6,nphmx),ag(3,3,ngrmx),
     $     wgt(ngrmx),axisph(0:3,3,0:nphmx),axisgr(0:3,3,ngrmx)
      real*8,intent(out)::athcs(6,ngrmx),celcs(6,6,ngrmx),
     $     mgel(6,6,ngrmx),selav(6,6)
      real*8,intent(inout)::AC_EL(6,6,ngrmx),bc_el(6,6,ngrmx)
     $     ,celav(6,6)

c     ----
      integer i,j,iph,kkk,k,it,iprint,inverr
      real*8
     $     xmtilde_el(6,6,ngrmx),celavv(6,6),selavv(6,6),
     $     athavv(6),AUX6(6),AUX33(3,3),AUX66(6,6),AUX3333(3,3,3,3),
     $     bux66(6,6),cux66(6,6),C4XC(3,3,3,3),C4XS(3,3,3,3),
     $     S4XC(3,3,3,3),S4XS(3,3,3,3),A2XS(3,3),A2XC(3,3),CNEW(6,6),
     $     CNEWX(6,6),CELUB(6,6),CELUBv(6,6),C4SA(3,3,3,3),
     $     E4SA(3,3,3,3),R4SA(3,3,3,3),ATHUB(6),ATHUBv(6),ATHAUX(6),
     $     AXB(3),EIGB(3,3),ctilde(6,6),AC1(6,6),AC2(6,6),
     $     AC3(6,6),AC4(6,6),ACAV(6,6),ACAVINV(6,6),ACINV(6,6),
     $     ACINVAV(6,6),ACINVAVINV(6,6),
     $     athav(6),TC(6,6,ngrmx)
c     dimensions introduced to study sg(6) scheme:
c      real*8 sav_sg6(6)
      real*8 tmismatch6x6,rer

      logical verbose
c      SAVE      CELAVOLD
      parameter(verbose=.false.)
c      parameter(verbose=.true.)

C ************************************************************************
C       IOPTION=0 : CALCULATES UPPER BOUND 'CUB' (VOIGT) FOR ELASTIC AND
C       THERMAL MODULI OF POLYCRYSTAL
C       USES AS INITIAL GUESS FOR THE SELF-CONSISTENT CALCULATION
C ************************************************************************


      if (verbose) then
         call w_chrc(0,'** entering suborutine elsc_evp')
         write(*,*)'ioption:',ioption
         call w_ival(0,'nph:',nph)
         do iph=1,nph
            call w_ival(0,'ngr_whole(iph):',ngr_whole(iph))
            do kkk=ngr_whole(iph-1)+1,ngr_whole(iph)
               call w_val(0,'wgt(kkk)):',wgt(kkk))
            enddo
         enddo
c         call stop_debug(0,'** stopped inside elsc_evp')
      endif


      iprint=1   ! controls detailed printing

C ************************************************************************
      IF(IOPTION.EQ.0   .OR.   IOPTION.eq.1) THEN
C ************************************************************************
C *** ROTATES SINGLE CRYSTAL STIFFNESS, COMPLIANCE AND THERMAL EXPANSION
C     TENSORS (CELCC, SELCC, ATHCC) TO SAMPLE AXES.
C *** CALCULATES UPPER BOUND (VOIGT) STIFFNESS 'CELUB' FOR THE POLYCRYSTAL.
C *** CALCULATES UPPER BOUND THERMAL EXPANSION TENSOR 'ATHUB' FOR THE PX.
C *** BOUNDS USED AS INITIAL GUESS FOR SELF-CONSISTENT CALCULATION
C ************************************************************************

         CELUB(:,:)=0d0
         ATHUB(:)  =0d0

c         DO IPH=IPHBOT,IPHTOP
         do iph=1,nph

C *** TRANSFORMS (ATHCC,CELCC,SELCC) FROM b-BASIS TO CARTESIAN (A2XC,C4XC,S4XC)
C     AND ROTATES INTO SAMPLE AXES (A2XS,C4XS,S4XS) FOR EVERY GRAIN

!     phase-specific properties converted to cartesian.
            aux6(:)=athcc(:,iph)
            CALL CHG_BASIS(aux6,A2XC,AUX66,AUX3333,1,6)
            aux66(:,:)=celcc(:,:,iph)
            CALL CHG_BASIS(AUX6,AUX33,aux66,C4XC,3,6)
            aux66(:,:)=selcc(:,:,iph)
            CALL CHG_BASIS(AUX6,AUX33,aux66,S4XC,3,6)

            if (verbose) then
               call w_chrc(0,'celcc(iph)')
               aux66(:,:)=celcc(:,:,iph)
c               call w_mdim(0,aux66,6,1d0)
               write(*,'(6e20.13)')((aux66(i,j),j=1,6),i=1,6)
               call w_chrc(0,'selcc(iph)')
               aux66(:,:)=selcc(:,:,iph)
c               call w_mdim(0,aux66,6,1d0)
               write(*,'(6e20.13)')((aux66(i,j),j=1,6),i=1,6)
c               call stop_debug(0,'** celcc and selcc')
            endif


!     Convert the phase-specific properties to grain-specific ones
         DO KKK=ngr_whole(iph-1)+1,ngr_whole(iph)
            aux33(:,:)=ag(:,:,kkk)
            if (verbose) then
               call w_chrc(0,'aux33(ag)')
c               call w_mdim(0,aux33,3,1d0)
               write(*,'(3e20.13)')((aux33(i,j),j=1,3),i=1,3)
               call w_chrc(0,'a2xc')
c               call w_mdim(0,a2xc,3,1d0)
               write(*,'(3e20.13)')((a2xc(i,j),j=1,3),i=1,3)
c               call stop_debug(0,'** stopped inside elsc_evp')
            endif
            call rot_2nd_tensor(a2xc,aux33,a2xs)

c            call w_chrc(0,'a2xs')
c     call w_mdim(0,a2xc,3,1d0)
c            write(*,'(3e20.13)')((a2xs(i,j),j=1,3),i=1,3)

            CALL CHG_BASIS(aux6,A2XS,AUX66,AUX3333,2,6)
            athcs(:,kkk)=aux6(:)
            call rot_4th_tensor(c4xc,aux33,c4xs) ! ag: sa<-ca
            call rot_4th_tensor(s4xc,aux33,s4xs) ! ag: sa<-ca

C *** FROM CARTESIAN TO B-BASIS
            CALL CHG_BASIS(AUX6,AUX33,aux66,C4XS,4,6)
            if (verbose) then
               call w_chrc(0,'aux66')
               call w_mdim(0,aux66,6,1d0)
            endif
            CELCS(:,:,kkk)=aux66(:,:)
            call chg_basis(aux6,aux33,aux66,s4xs,4,6)
            mgel(:,:,kkk)=aux66(:,:)

            if (verbose) then
               call w_chrc(0,'celcs(:,:,kkk)')
               write(*,'(6e20.13)')((celcs(i,j,kkk),j=1,6),i=1,6)
            endif
            do i=1,6
            do j=1,6
               aux66(i,j)=0d0
            do k=1,6
               aux66(i,j)=aux66(i,j)+celcs(i,k,kkk)*mgel(k,j,kkk)
            enddo
            enddo
            enddo
            if (verbose) then
               call w_chrc(0,'aux66')
               call w_mdim(0,aux66,6,1d0)
            endif
            ATHUB(:)  =ATHUB(:)  +ATHCS(:,KKK)*WGT(KKK)
            CELUB(:,:)=CELUB(:,:)+CELCS(:,:,KKK)*WGT(KKK)
         enddo
         enddo
CFEB
         if (ioption.eq.0) then

C *** SAVES UPPER BOUND IN B-BASIS FOR USING AS NEXT GUESS
c            CELAV(:,:)=CELUB(:,:) ! (E)VPSC standalone
            if (verbose) then
               call w_ival(0,'ioption:',ioption)
               call w_chrc(0,'celub:')
               call w_mdim(0,celub,6,1d0)
            endif
            celav(:,:)=celub(:,:) ! FE version

C *** FROM B-BASIS TO CARTESIAN TO VOIGT

            CALL CHG_BASIS(ATHUB ,AUX33,AUX66,AUX3333,1,6)
            CALL VOIGT    (ATHUBv,AUX33,AUX66,AUX3333,2)
            CALL CHG_BASIS(AUX6,AUX33,CELUB ,AUX3333,3,6)
            CALL VOIGT    (AUX6,AUX33,CELUBv,AUX3333,4)
         endif ! END OF IOPTION=0
CFEE
      ENDIF                     ! END OF IOPTION=0 or IOPTION=1

C *************************************************************************
C    IF IOPTION=1 : USES PREVIOUS CONVERGED MODULI 'CELAVOLD' AS INITIAL
C    GUESS FOR A SELF-CONSISTENT CALCULATION OF THE ELASTIC STIFFNESS AND
C    THE THERMAL MODULI
C *************************************************************************

      IF(IOPTION.EQ.1) THEN
C *************************************************************************
         IT=0
         RER=2*1.E-05
         DO WHILE(RER.GT.1.E-05 .AND. IT.LT.ITMAXINT)

            IT=IT+1
            CNEWX(:,:)   =0d0
            ACAV(:,:)    =0d0
            ACINVAV(:,:) =0d0
            ATHAUX(:)     =0d0

            if (verbose) then
               write(*,*)'CELAV'
               do i=1,6
                  write(*,'(6e15.7)')(celav(i,j),j=1,6)
               enddo
            endif
            aux66(:,:)=celav(:,:)
            CALL CHG_BASIS(AUX6,AUX33,aux66,C4SA,3,6)
            if (verbose) then
               write(*,'(a,e11.3)')'c4sa(1,1,1,1):',c4sa(1,1,1,1)
               write(*,'(6e11.3)')((celav(i,j),j=1,6),i=1,6)
            endif

c         DO IPH=IPHBOT,IPHTOP
            do iph=1,nph

C *** SAME SHAPE FOR ALL GRAINS IN THE PHASE
               IF(ISHAPE(IPH).LE.1) THEN
                  AXB(:)     =AXISPH(0,:,IPH)
                  EIGB(1:3,:)=AXISPH(1:3,:,IPH)
               ENDIF

               DO KKK=ngr_whole(iph-1)+1,ngr_whole(iph)
C *** DIFFERENT SHAPE FOR EACH GRAIN IN THE PHASE
                  IF(ISHAPE(IPH).GE.2) THEN ! INDIVIDUAL GRAIN SHAPE
                     AXB(:)     =AXISGR(0,:,kkk)
                     EIGB(1:3,:)=AXISGR(1:3,:,kkk)
                  ENDIF

C *** ISHAPE=0,1 --> SAME SHAPE FOR ALL GRAINS IN THE PHASE
C *** ISHAPE=2,3 -->DIFFERENT SHAPE FOR EACH GRAIN IN THE PHASE

                  if (verbose) then
                     write(*,*)'iph:',iph
                     write(*,*)'axisph(0,:,iph)'
                     write(*,*)(axisph(0,i,iph),i=1,3)
                  endif


                  IF(ISHAPE(IPH).GE.2
     $                 .OR. KKK.EQ.NGR_whole(IPH-1)+1) THEN

                     if (verbose) then
                        write(*,*)'axb'
                        write(*,*)(axb(i),i=1,3)
                        call w_chrc(0,'**eigb')
                        call w_mdim(0,eigb,3,1d0)
                        write(*,'(a,e11.4)') 'c4sa(1,1,1,1):',
     $                       c4sa(1,1,1,1)
                     endif

c                    R4SA is not used
                     call eshelby_rot(axb,eigb,c4sa,e4sa,r4sa,1)

                     if (verbose) then
                        write(*,'(a,e11.4)') 'e4sa(1,1,1,1):',
     $                       e4sa(1,1,1,1)
                     endif

C *** ROTATES THE ESHELBY TENSOR FOR PHASE OR GRAIN BACK INTO SAMPLE AXES.
c               CALL CHG_BASIS(AUX6,AUX33,AUX66,E4SA,0,6) ! suppressed since this has been done
                     aux66(:,:)=0d0
                     CALL CHG_BASIS(AUX6,AUX33,AUX66,E4SA,4,6)

                     if (verbose) then
                        write(*,*)'here1-1'
                        write(*,*)'aux66'
                        write(*,'(6e15.7)')((aux66(i,j),j=1,6),i=1,6)
                     endif

                     call lu_inverse(inverr,AUX66,6) ! Eshg^(-1)
                     if (inverr.ne.0) then
                        write(*,*)'** error in elsc_evp_fe(1)'
                        write(*,*)'here1-2'
                        write(*,*)'aux66'
                        write(*,'(6e15.7)')((aux66(i,j),j=1,6),i=1,6)
                        write(*,*)'CELAV'
                        do i=1,6
                           write(*,'(6e15.7)')(celav(i,j),j=1,6)
                        enddo
                        write(*,*)'iph:',iph
                        write(*,*)'AXISPH(0,:,IPH)'
                        write(*,*)(AXISPH(0,i,IPH),i=1,3)
                        write(*,*)'axb'
                        write(*,*)(axb(i),i=1,3)
                        write(*,*)'AXISPH(1:3,:,IPH)'
                        do i=1,3
                           write(*,*)(AXISPH(i,j,IPH),j=1,3)
                        enddo
                        call w_chrc(0,'**eigb')
                        call w_mdim(0,eigb,3,1d0)
                        call xit()
                     endif

                     DO I=1,6
                        DO J=1,6
                           CTILDE(I,J)=0d0
                           DO K=1,6
                              CTILDE(I,J)=CTILDE(I,J)
     $                             +CELAV(I,K)*
     $                             (AUX66(K,J)-(K/J)*(J/K)) !XID6(K,J))
                           ENDDO
                           AC2(I,J)=CELAV(I,J)+CTILDE(I,J) ! (Cav+C~)
                        enddo
                     enddo

c                 EOI(ISHAPE(IPH).GE.2.OR.KKK.EQ.NGR_whole(IPH-1)+1)
                  ENDIF

c            xctilde_el(:,:,kkk)=ctilde(:,:)
                  aux66(:,:)=ctilde(:,:)
                  call lu_inverse(inverr,aux66,6)
                  if (inverr.ne.0) then
                     write(*,*)'** error in elsc_evp_fe(3)'
                     call xit()
                  endif
                  xmtilde_el(:,:,kkk)=aux66(:,:)

                  DO I=1,6
                     DO J=1,6
                        AC1(I,J)=CELCS(I,J,kkk)+CTILDE(I,J) ! (Cg+C~)
                     enddo
                  enddo
                  call lu_inverse(inverr,AC1,6) !(Cg+C~)^-1
                  if (inverr.ne.0) then
                     write(*,*)'** error in elsc_evp_fe(4)'
                     call xit()
                  endif
                  DO I=1,6
                     DO J=1,6
                        AC3(I,J)=0d0
                        AC4(I,J)=0d0
                        DO K=1,6
c                          Ag=(Cg+C~)^(-1):(Cav+C~) ! interaction matrix
                           AC3(I,J)=AC3(I,J)+AC1(I,K)*AC2(K,J)
c                          Tg=(Cg+C~)^(-1): C~
                           AC4(I,J)=AC4(I,J)+AC1(I,K)*CTILDE(K,J)
                        enddo
                     enddo
                  enddo

c                 ACINV used for calculating thermal moduli
                  ACINV(:,:)=AC3(:,:)
                  CALL LU_INVERSE(inverr,ACINV,6)
                  if (inverr.ne.0) then
                     write(*,*)'** error in elsc_evp_fe(5)'
                     write(*,*)'celav'
                     write(*,'(6e12.4)')((CELAV(i,j),j=1,6),i=1,6)
                     write(*,*)'celcs(:,:,kkk)'
                     write(*,'(6e14.4)')((celcs(i,j,kkk),j=1,6),i=1,6)
                     call xit()
                  endif

                  DO I=1,6
                     DO J=1,6
                        AUX66(I,J)=0d0
                        DO K=1,6
                           AUX66(I,J)=AUX66(I,J)+
     $                          CELCS(I,K,kkk)*AC3(K,J)
                        ENDDO
c                       <Cg*Ag>
                        CNEWX(I,J)=CNEWX(I,J)+AUX66(I,J)     *WGT(KKK)
c                       <Ag>
                        ACAV(I,J) =ACAV(I,J) +AC3(I,J)       *WGT(KKK)
c                       <Ag^(-1)>
                        ACINVAV(I,J) =ACINVAV(I,J)+ACINV(I,J)*WGT(KKK)
                     enddo
                  enddo

C *** CALCULATES TENSORS RELATED TO THE AVERAGE THERMAL MODULI
                  DO I=1,6
                     DO J=1,6
                        DO K=1,6
                           ATHAUX(I)=ATHAUX(I)+
     $                          ACINV(I,J)*CELCS(J,K,kkk)*
     $                          ATHCS(K,kkk)*WGT(kkk)
                        enddo
c                       strain concentration tensor 'AC_EL'
                        AC_EL(I,J,Kkk)=AC3(I,J)
c                       Thermal Concentration tensor 'TC'
                        TC(I,J,Kkk)=AC4(I,J)
                     enddo
                  enddo
               enddo
            enddo

            ACAVINV(:,:)=ACAV(:,:)
            call lu_inverse(inverr,ACAVINV,6)
            if (inverr.ne.0) then
               write(*,*)'** error in elsc_evp_fe(5)'
               write(*,*)'ctilde'
               write(*,'(6e12.4)')((ctilde(i,j),j=1,6),i=1,6)
               write(*,*)'celav'
               write(*,'(6e12.4)')((CELAV(i,j),j=1,6),i=1,6)
               write(*,*)'celcs(:,:,kkk)'
               write(*,'(6e14.4)')((celcs(i,j,kkk),j=1,6),i=1,6)
               write(*,*)'** error in elsc_evp_fe(6)'
               call xit()
            endif

            DO I=1,6
               DO J=1,6
                  CNEW(I,J)=0d0
                  DO K=1,6
                     CNEW(I,J)=CNEW(I,J)+CNEWX(I,K)*ACAVINV(K,J)
                  enddo
               enddo
            enddo

c         RER=TMISMATCH6x6(CELAV,CNEW) ! stand-alone version
            do i=1,6
               do j=1,6
                  aux66(i,j)=celav(i,j)
               enddo
            enddo
            RER=TMISMATCH6x6(aux66,CNEW) ! stand-alone version


C *** SAVES THE SYMMETRIC TENSOR AS THE NEW GUESS
C *** At least within significant 'cnew' is symmetric
            DO I=1,6
               DO J=I,6
                  celav(i,j)=0.5d0*(CNEW(I,J)+CNEW(J,I))
                  celav(j,i)=celav(i,j)
               enddo
            enddo

c         WRITE( *,'(1H+,''  EL ITER'',I4,''    REL ERR'',E11.3)')IT,RER
c         call xit()
         enddo                  ! END OF (DO..WHILE)

c         CELAVOLD(:,:)=CELAV(:,:)   ! SAVES THE SELF CONSISTENT MODULI
c         SELAV(:,:)=CELAV(:,:)
         do i=1,6
            do j=1,6
               aux66(i,j)=CELAV(i,j)
            enddo
         enddo
         call lu_inverse(inverr,aux66,6)
         if (inverr.ne.0) then
            write(*,*)'** error in elsc_evp_fe(7)'
            write(*,*)' Could not invert celav'
            call xit()
         endif
         do i=1,6
            do j=1,6
               SELAV(i,j)=aux66(i,j)
            enddo
         enddo

c *** Calculates the stress concentration tensor 'Bc_el'
         cux66(:,:)=0d0


         do iph=1,nph
         DO KKK=ngr_whole(iph-1)+1,ngr_whole(iph)
            aux66(:,:)=0d0
            bux66(:,:)=0d0


            if (verbose) then
               write(*,*)'mgel'
               do i=1,6
                  write(*,'(6e20.13)')(mgel(i,j,kkk),j=1,6)
               enddo
               write(*,*)'xmtilde_el'
               do i=1,6
                  write(*,'(6e20.13)')(xmtilde_el(i,j,kkk),j=1,6)
               enddo
            endif

            aux66(:,:)=mgel(:,:,kkk)+xmtilde_el(:,:,kkk)
            call lu_inverse(inverr,aux66,6)
               if (inverr.ne.0) then
                  write(*,*)'** error in elsc_evp_fe(8)'
                  call xit()
               endif
            bux66(:,:)=SELAV(:,:)+xmtilde_el(:,:,kkk)


            if (verbose) then
               write(*,*)'aux66'
               do i=1,6
                  write(*,'(6e20.13)')(aux66(i,j),j=1,6)
               enddo
               write(*,*)'bux66'
               do i=1,6
                  write(*,'(6e20.13)')(bux66(i,j),j=1,6)
               enddo
            endif
            do i=1,6
               do j=1,6
                  bc_el(i,j,kkk)=0d0
                  do k=1,6
                     bc_el(i,j,kkk)=bc_el(i,j,kkk)
     $                    +aux66(i,k)*bux66(k,j)
                  enddo
               enddo
            enddo

            if (verbose) then
               write(*,*)'bc_el'
               do i=1,6
                  write(*,'(6e15.8)')(bc_el(i,j,kkk),j=1,6)
               enddo
c               call stop_debug(0,'** stopped inside elsc_evp')
            endif
         enddo
         enddo

C *** CALCULATES THE SELF-CONSISTENT AVERAGE THERMAL MODULUS 'ATHAV'
         ACINVAVINV(:,:)=ACINVAV(:,:)
         CALL LU_INVERSE(inverr,ACINVAVINV,6)
               if (inverr.ne.0) then
                  write(*,*)'** error in elsc_evp_fe(9)'
                  call xit()
               endif

         DO I=1,6
            ATHAV(I)=0d0
         DO J=1,6
         DO K=1,6
            ATHAV(I)=ATHAV(I)+SELAV(I,J)*ACINVAVINV(J,K)*ATHAUX(K)
         enddo
         enddo
         enddo

      ENDIF                     ! END OF IOPTION.EQ.1


c$$$C ************************************************************************
c$$$C *** IOPTION=2: CALCULATES STRESS, ELASTIC STRAIN, THERMAL STRAIN IN
c$$$C                POLYCRYSTAL AND GRAINS BASED ON LOAD CONDITIONS
c$$$C ************************************************************************
c$$$      IF(IOPTION.EQ.2) THEN !! this should be done for 'VPSC',
c$$$c                     when the right type of loading condition is imposed.
c$$$C ************************************************************************
c$$$
c$$$      iprint=1
c$$$
c$$$      TEMP_INCR =TEMP/NSTEPS   ! for using with cooling simulations
c$$$C      TIME_INCR =TIME_INCR     ! TIME_INCR=eijincr/abs(dsimctrl) calculated in subr load_conditions
c$$$C      DSIM(:,:) =DSIM(:,:)     ! total strain rate tensor input in subr load_conditions
c$$$
c$$$      SAVc(:,:)    =SCAUCHY(:,:)   ! may adopt later SAV instead of SCAUCHY
c$$$
c$$$      SAVcINCR(:,:)= SAVc(:,:)      ! valid only if imposing zero stress BCs (subr load_conditions)
c$$$      EELcINCR(:,:)=(DSIM(:,:)-DAV_TRANSc(:,:))*TIME_INCR   ! accounts for ETWcINCR(:,:) if applicablE
c$$$
c$$$C *** IN SUBROUTINE STATE_6x6 TENSORS 'dsim & scauchy & selav' ENTER IN
c$$$C     CARTESIAN, GET INTERNALLY CONVERTED TO VOIGT, AND EXIT IN CARTESIAN.
c$$$C *** IMPOSED 'dsim & scauchy' COMPONENTS MAY BE MIXED
c$$$C *** THIS CALL MAY ONLY BE VALID FOR ZERO STRESS BOUNDARY CONDITIONS !!
c$$$
c$$$        CALL VOIGT(AUX6,AUX33,SELAVv,SELAVc,3)
c$$$        CALL STATE_6x6(iDBARv,EELcINCR,iSBARv,SAVcINCR,SELAVc)
c$$$
c$$$        DELAVc(:,:)=EELcINCR(:,:)/TIME_INCR
c$$$        SAVc(:,:) =SAVc(:,:)+SAVcINCR(:,:)
c$$$        SCAUCHY(:,:)=SAVc(:,:)      ! may adopt later SAV instead of SCAUCHY
c$$$
c$$$        CALL CHG_BASIS(DELAV,DELAVc ,AUX66,AUX3333,2,6)
c$$$        CALL CHG_BASIS(SAV  ,SAVc   ,AUX66,AUX3333,2,6)
c$$$
c$$$        DO I=1,6
c$$$          ETHAV(I)=ETHAV(I)+ATHAV(I)*TEMP_INCR
c$$$          EELAV(I)=EELAV(I)+DELAV(I)*TIME_INCR
c$$$        ENDDO
c$$$        CALL CHG_BASIS(EELAV,EELAVc,AUX66,AUX3333,1,6)
c$$$        CALL CHG_BASIS(ETHAV,ETHAVc,AUX66,AUX3333,1,6)
c$$$
c$$$C --> TEMPORARY & FOR OUTPUT PURPOSES
c$$$        SBAR(ipx,1:5)=SAV(1:5)
c$$$        CALL CHG_BASIS(DBAR,DSIM,AUX66,AUX3333,2,5)
c$$$
c$$$
c$$$        sigavxB(:)=0d0   ! for pressure of 2-phase
c$$$        eelavxB(:)=0d0
c$$$        SAVPH(:,:)=0d0
c$$$        SAVPHDEV(:,:)=0d0
c$$$
c$$$c        DO IPH=IPHBOT,IPHTOP
c$$$        do iph=iph0,iph1
c$$$        DO KKK=ngr_whole(iph-1)+1,ngr_whole(iph)
c$$$
c$$$          DELGR(:,kkk)=0d0
c$$$          DO I=1,6
c$$$          DO J=1,6
c$$$            DELGR(I,kkk)=DELGR(I,kkk)+AC_EL(I,J,kkk)*DELAV(J)-
c$$$     #      TC(I,J,kkk)*(ATHCS(J,Kkk)-ATHAV(J))*TEMP_INCR/TIME_INCR
c$$$          ENDDO
c$$$c            EELGR(I,kk)=EELGR(I,kkk)+DELGR(I,kkk)*TIME_INCR
c$$$          ENDDO
c$$$
c$$$c$$$          DO I=1,6
c$$$c$$$            SG(I,kkk)=0.
c$$$c$$$            DO J=1,6
c$$$c$$$              SG(I,kkk)=SG(I,kkk)+CELCS(I,J,kkk)*EELGR(J,kkk)
c$$$c$$$            ENDDO
c$$$c$$$          ENDDO
c$$$
c$$$C *** CHECKING CONSISTENCY OF GRAIN STRESS AND STRAIN   *********
c$$$          sigavxB(:)=sigavxB(:)+sg(:,KKK)*wgt(KKK)
c$$$          eelavxB(:)=eelavxB(:)+eelgr(:,KKK)*wgt(KKK)
c$$$
c$$$          WGTXX= WGT(KKK)/WPH(IPH)
c$$$          SAVPH(:,IPH)   =SAVPH(:,IPH)   +sg(:,KKK)*wgtxx
c$$$          SAVPHDEV(:,IPH)=SAVPHDEV(:,IPH)+sg(:,KKK)**2*wgtxx
c$$$
c$$$        ENDDO      ! END OF DO KKK
c$$$
c$$$        SAVPHDEV(:,IPH)=sqrt(ABS(SAVPHDEV(:,IPH)-SAVPH(:,IPH)**2))   ! for pressure of 2-phase
c$$$
c$$$        ENDDO      ! END OF DO IPH
c$$$
c$$$      if(iprint.eq.1) then
c$$$c        CALL CHG_BASIS(EELAVxb,EELAVxc,AUX66,AUX3333,1,6)
c$$$c        CALL CHG_BASIS(SIGAVxb,SIGAVxc,AUX66,AUX3333,1,6)
c$$$      endif
c$$$
c$$$      ENDIF      ! END OF IF(IOPTION.EQ.2)

c----------------------------------------------------------------
c     A study to compare with sg6 scheme of VPSC 2018/08/13
c     Calculates sg(6) and estimate the hydrostatic pressure.
c     Writes local cauchy stress component <sg_11> and pressure,
c     and the grain weight.
c----------------------------------------------------------------
c$$$      if (ioption.eq.3) then !! meant to be used in combination with VPSC
c$$$
c$$$!     savc(:,:)=scauchy(:,:)
c$$$         aux33(:,:)=scauchy(:,:)
c$$$         call chg_basis(sav_sg6,aux33,aux66,aux3333,2,6)
c$$$
c$$$c         do iph=iphbot,iphtop
c$$$         do iph=iph0,iph1
c$$$         do kkk=ngr_whole(iph-1)+1,ngr_whole(iph)
c$$$            sg6(:5)=sg(:5,kkk)
c$$$            sg6(6)=0d0
c$$$            do j=1,6
c$$$               sg6(6)=sg6(6)+bc_el(6,j,kkk)*sav_sg6(j)
c$$$            enddo
c$$$c           Chg_b-basis sg6 to cartesian
c$$$            call chg_basis(sg6,aux33,aux66,aux3333,1,6)
c$$$c            aux33(1,1)          !! sg_11
c$$$         enddo
c$$$         enddo
c$$$      endif

      if (verbose) call w_chrc(0,'** Exiting suborutine elsc_evp')

      RETURN

      END subroutine elsc_evp_fe
