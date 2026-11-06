c     -----------------------------------------------------------------
      subroutine update_macro_property(
c     in
     $     mdim,
     $     iregime,
     $     errm,
     $     itmaxint,
     $     nphmx,
     $     ngrmx, !! this is in fact ngrmx from <evpsc.f>
     $     nph,
     $     ngr,
     $     axisph,
     $     axisgr,
     $     ishape,
     $     dg0,
     $     mgtg,
     $     wgt,factor,
     $     dg_trans,

c     inout
     $     LBARTG,
     $     mbartg,
     $     dbar0,

c     out

     $     SG0AV,
     $     BGINV,
     $     asph,
     $     asgr,
     $     mastph,
     $     mastgr,
     $     iesh_fail,
     $     inverr
     $     )
      implicit none
      integer, intent(in)::mdim,iregime,itmaxint,nph,nphmx,ngrmx,
     $     ishape(0:nphmx),ngr(0:nphmx)
      integer, intent(out):: inverr
      real*8,intent(in)::dg0(6,ngrmx),mgtg(6,6,ngrmx),
     $     wgt(ngrmx),errm,axisph(0:3,3,0:nphmx),axisgr(0:3,3,ngrmx),
     $     factor(ngrmx),dg_trans(6,ngrmx)
      real*8,intent(inout)::LBARTG(6,6),mbartg(6,6),dbar0(6)
      real*8,intent(out)::asph(3,3,3,3,nphmx),
     $     asgr(3,3,3,3,ngrmx),SG0AV(mdim),BGINV(mdim,mdim),
     $     mastph(mdim,mdim,nphmx),mastgr(mdim,mdim,ngrmx)
      logical, intent(out)::iesh_fail
c     locals
      integer it1tg,iph,kkk,i,j,k
      real*8 rer1tg,tmismatch
      real*8 bgav(mdim,mdim),BGTAV(mdim,mdim),BGTINV(mdim,mdim),
     $     BGTD0AV(6),bc1(mdim,mdim),bc2(mdim,mdim),
     $     bgx(mdim,mdim),sg0(mdim),bgtx(6,6),
     $     mbarnew(mdim,mdim),mast(mdim,mdim),
     $     MGTGBG(6,6),MGTGBGAV(mdim,mdim),
     $     auxm(mdim),auxmm(mdim,mdim)
      inverr=0

      it1tg=0
      rer1tg=2*errm

      do 500 while(rer1tg.gt.errm.and.it1tg.lt.itmaxint)
         it1tg=it1tg+1
c        --
c        Calculates both symmetric and skew-symmetric part of Eshelby tensor
c        as well as the interaction compliance Sevp~ = (I-E)^-1:E:<Sevp>
c        * mastph (ishape.le.1)
c        * mastgr (ishape.ge.2)
c        * asph
c        * asgr

         call calc_evptilde(
     $        mdim,
     $        iregime,          ! in
     $        LBARTG,             ! in
     $        mbartg,             ! in
     $        nphmx,
     $        nph,              ! in

c             (this is actually ngrmx from <evpsc> via <update_macro_property>
     $        ngrmx,            ! in
     $        ngr,              ! in
     $        factor,           ! in
     $        axisph,           ! in
     $        axisgr,           ! in
     $        ishape,           ! in
     $        asph,             ! out
     $        asgr,             ! out
     $        mastph,         ! out
     $        mastgr,         ! out
     $        iesh_fail         ! out
     $        )
         if (iesh_fail) return


C **************************************************************************
C *** LOOP #2 OVER PHASES AND GRAINS
C **************************************************************************

C *** SOLVES ITERATIVELY SC EQUATION FOR GENERAL CASE OF DIFFERENT SHAPES.
C *** REQUIRED FOR MULTI-PHASE WHEN GRAINS OF EACH PHASE HAVE DIFFERENT SHAPE,
C     OR FOR ISHAPE>1, WHEN EVERY GRAIN HAS DIFFERENT SHAPE, OR BOTH.
C *** EQUATION NUMBERS REFER TO Section 1.5 IN VPSC MANUAL

c        ------------------------------------------------
c        Some weighted average estimates are initialized
         BGTD0AV(:)=0.          !<B*:(e0 + d_trans)>
         SG0AV(:)=0.            !<b>
         bgav(:,:)=0.           !<B>
         BGTAV(:,:)=0.          !<B*>
         MGTGBGAV(:,:)=0.

c        ------------------------------------------------

         do iph=1,nph
            if (ishape(iph).le.1) then
               mast(:,:)=mastph(:,:,iph)
               bc2(:mdim,:mdim)=mast(:mdim,:mdim)+mbartg(:mdim,:mdim)
            endif

            do kkk=ngr(iph-1)+1,ngr(iph)

               if (ishape(iph).ge.2) then ! update for each grain
                  mast(:,:)=mastgr(:,:,kkk)
                  bc2(:mdim,:mdim)=mast(:mdim,:mdim)+mbartg(:mdim,:mdim)
               endif
c              --
c              Calculate the big B tensor; First obtain bc1 and bc2, multiplicative result of which gives
c              Bgx = bc1 x bc2 - the double b tensor
c              bc1 = (Mc+M~)^-1; bc2 = (Mbartg+M~)
               do i=1,mdim
               do j=1,mdim
                  bc1(i,j)=mast(i,j)+mgtg(i,j,kkk)
               enddo
               enddo
               call lu_inverse(inverr,bc1,mdim)

               bgx(:,:)=0d0
               bgtx(:,:)=0d0 ! B*
               do i=1,mdim
               do j=1,mdim
               do k=1,mdim
c                 bgx = (Mc+M~)^-1 : (Mbartg+M~)
c                 bgx     =          bc1     :bc2
                  bgx(i,j) = bgx(i,j)+bc1(i,k)*bc2(k,j)
                  bgtx(i,j)=bgtx(i,j)+bc2(i,k)*bc1(k,j)
               enddo
c$$$          BG(I,J,KGX) =BGX(I,J)        ! Eq (5-34)  ! used in GRAIN_STRESS_ALT
               enddo
               enddo

c$$$C *** SEPARATE 5-D CALCULATION FOR SECOND ORDER APPROACH ARRAYS
c$$$        DO I=1,5
c$$$        DO J=1,5
c$$$          CHIFLU(I,J,KGX)=BC1(I,J)                        ! (M* + Mc,tg)^-1
c$$$          DUMMY2=0.
c$$$          DO K=1,5
c$$$            DUMMY2=DUMMY2 +MGTG(I,K,KGX)*BC1(K,J)
c$$$          ENDDO
c$$$          BETFLU(I,J,KGX)=DUMMY2
c$$$        ENDDO
c$$$        ENDDO

c              --
c              Calculate small b tensor
c              sg0 = (Mc+M~)^-1:(dbar0-dg0-dg_trans)
               sg0(:)=0d0
               do i=1,mdim
               do j=1,mdim
                  sg0(i)=sg0(i)+bc1(i,j)
     $                 * (dbar0(j)-dg0(j,kkk)-dg_trans(j,kkk))
                  BGTD0AV(i)=BGTD0AV(i) +bgtx(i,j)*
     $                 (dg0(j,kkk)+dg_trans(j,kkk))*wgt(kkk)
               enddo
               enddo

               DO I=1,MDIM
               DO J=1,MDIM
                  MGTGBG(I,J) =0.
                  DO K=1,MDIM
c                     MGTGBG(I,J)=MGTGBG(I,J)+MGTG(I,K,KGX)*BG(K,J,KGX) ! Eq 5-40a (yj: equation label is not consistent with VPSC8 book)
                     MGTGBG(I,J)=MGTGBG(I,J)+MGTG(I,K,KKK)*BGX(K,J)
                  ENDDO
                  ENDDO
               ENDDO

c              --
c              Below weighted averages are calculated:
c              _______________________________________
c              SG0AV     = <b>
c              BGTD0AV= <B*:(d0 + d_trans)>
c              bgav      = <Bgx>
c              BGTAV     = <Bgtx>
               do i=1,mdim
                  SG0AV(i)     = SG0AV(i)     +     sg0(i)*wgt(kkk)
               do j=1,mdim
                  BGAV(I,J)    = BGAV(I,J)    + BGX(I,J)*WGT(KKK)
                  BGTAV(i,j)   = BGTAV(i,j)   +bgtx(i,j)*wgt(kkk)
                  MGTGBGAV(I,J)=MGTGBGAV(I,J)+MGTGBG(I,J) *WGT(Kkk) ! <Mg_tg+Bg> in Eqs 5-40a/5-41a
               enddo
               enddo
c           =====
            enddo               ! loop over grains
         enddo                  ! loop over phases


C ****************************************************************
C *** END OF LOOP #2 OVER GRAINS AND PHASES
C ****************************************************************

c     BGINV   = <Bgx>^-1
c     BGTINV =  <Bgtx>^-1
         BGINV(1:MDIM,1:MDIM) =BGAV(1:MDIM,1:MDIM)
         BGTINV(1:MDIM,1:MDIM)=BGTAV(1:MDIM,1:MDIM)

         call lu_inverse(inverr,BGINV,MDIM)
         call lu_inverse(inverr,BGTINV,MDIM)

         mbarnew(:,:)=0d0
         do i=1,mdim
         do j=1,mdim
         do k=1,mdim
            mbarnew(i,j)=mbarnew(i,j)+MGTGBGAV(i,k)*BGINV(k,j)
         enddo
         enddo
         enddo

         auxmm(1:mdim,1:mdim)=mbartg(1:mdim,1:mdim)
         rer1tg = tmismatch(auxmm,mbarnew,mdim)
c     updates mbartg and dbar0

         mbartg(1:mdim,1:mdim)=mbarnew(1:mdim,1:mdim)
         auxmm(1:mdim,1:mdim)=mbarnew(1:mdim,1:mdim)
         call lu_inverse(inverr,auxmm,mdim)
         lbartg(1:mdim,1:mdim)=auxmm(1:mdim,1:mdim)

         dbar0(:)=0d0
         do i=1,mdim
         do j=1,mdim
            dbar0(i) = dbar0(i) + BGTINV(i,j)*BGTD0AV(j)
         enddo
         enddo

c     iteration on estimates on overall elasto-visco-plastic estimates
 500  continue

      return
      end subroutine update_macro_property
