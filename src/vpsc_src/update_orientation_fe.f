C **************************************************************************
C     SUBROUTINE UPDATE_ORIENTATION_fe    --->   VERSION 12/MAR/2020
C
C     USES SHEAR RATES AND GRAIN SHAPE TO CALCULATE GRAIN VELOCITY GRADIENT.
C     UPDATES GRAIN ORIENTATIONS DUE TO CRYSTALLOGRAPHIC SHEAR (SLIP &
C     TWINNING) BUT DOES NOT PERFORM TWIN REORIENTATION.
C
C     CNT: modified to use only phase-associated tensors when ISHAPE.LE.1
C     RAL: split in 3 DO loops to deal with co-rotations (17/02/00)
C     CNT: moved last loop to subroutine UPDATE_TWINNING (05/06/02)
C     CNT: fixed bug in LIJ calculation when IFLAT=1     (21/JUN/2008)
C     CNT: replaced subr REORIENT GRAIN with RODRIGUES   (OCT/2009)
C **************************************************************************

      subroutine update_orientation_fe(
c          in
     $     dtime,interaction,nsysmx,ngrmx,nphmx,igb,igt,iph,
     $     ishape,wph,WBARc,asph,asgr,tg,DBAR,
     $     dnca,dbca,gamdot,wgt,nneigh,nneimx,neigh,wneigh,
c          inout
     $     ag,
c          out
     $     xlijph,xlijgr
     $     )
      implicit none
      integer,intent(in)::nneimx,nsysmx,ngrmx,interaction,nneigh,
     $     nphmx,ishape(0:nphmx),neigh(0:nneimx,ngrmx),igb,igt,iph
      real*8, intent(in)::dtime,WBARc(3,3),asph(3,3,3,3,nphmx),
     $     asgr(3,3,3,3,ngrmx),tg(6,ngrmx),DBAR(6),wph(nphmx),
     $     dnca(3,nsysmx,nphmx),dbca(3,nsysmx,nphmx),
     $     gamdot(nsysmx,ngrmx),wgt(ngrmx),
     $     wneigh(0:nneimx,ngrmx)
      real*8,intent(inout)::ag(3,3,ngrmx)
      real*8, intent(out)::xlijph(3,3,nphmx),xlijgr(3,3,ngrmx)
      real*8 as(3,3,3,3),east(3,3),dnsa(3),dbsa(3),rotslip(3,3),
     $     rotloc(3,3),corot(3,3),aux5(5),aux33(3,3),aux55(5,5),
     $     aux6(6),aux66(6,6),aux3333(3,3,3,3),rot(3,3,NGRmx),aa(3,3),
     $     arot(3,3),XLIJGRX(3,3),XLIJGR0(3,3)
     $
      integer i,j,kkk,in,k,l,is
      logical verbose
      parameter(verbose=.false.)
c      parameter(verbose=.true.)


!     to fool the compiler
      if (.false.) write(*,*) neigh(0,1)

      if(verbose)then
         call w_chrc(0,'** Entering UPDATE_ORIENTATION')
         call w_chrc(0,'WBARc')
         call w_mdim(0,WBARc,3,1d0)
         call w_chrc(0,'WBARc')
         do i=1,3
            write(*,'(3e20.13)')(WBARc(i,j),j=1,3)
         enddo
c      call xit()
      endif

cc    ROTPLASTAV=0.
cc    ROTLOCALAV=0.
cc    ROTTOTALAV=0.


c      KGX=1
      DO I=1,3
      DO J=1,3
         XLIJPH(I,J,iph)=0d0
      ENDDO
      ENDDO

      IF(ISHAPE(iph).LE.1) THEN
         DO I=1,3
         DO J=1,3
         DO K=1,3
         DO L=1,3
            AS(I,J,K,L)=ASPH(I,J,K,L,iph)
         ENDDO
         ENDDO
         ENDDO
         ENDDO
      ENDIF

C *** THE FOLLOWING TO AVOID EMPTY 'CHILD' PHASE IN COMP GRAIN MODEL
      if(wph(iph).lt.1.e-6) then
         write(*,*)'negligible phase orientation'
         go to 2000
      endif

      do 1000 kkk=igb,igt

C *** CALCULATES LOCAL ROTATION ROTLOC=PI*S**(-1)*(DG-DAV) FOR EVERY GRAIN.
C *** ROTLOC IS ZERO FOR TAYLOR CALCULATION.

         if(interaction.le.0) then
            do i=1,3
            do j=1,3
               rotloc(i,j)=0d0
            enddo
            enddo
         else if(interaction.gt.0) then
CFEB
            if(ishape(iph).ge.2) then
               do i=1,3
               do j=1,3
               do k=1,3
               do l=1,3
                  as(i,j,k,l)=asgr(i,j,k,l,kkk)
               enddo
               enddo
               enddo
               enddo
            endif
cfee
c           Not sure if we should use full 6D
            if (.true.) then ! original for VPSC
               do i=1,5
                  aux5(i)=tg(i,kkk)-DBAR(i)
               enddo
               call chg_basis(aux5,east,aux55,aux3333,1,5)
            else ! maybe more legitimate for dEVPSC.
               do i=1,6
                  aux6(i)=tg(i,kkk)-DBAR(i)
               enddo
               call chg_basis(aux6,east,aux66,aux3333,1,6)
            endif

            do i=1,3
            do j=1,3
               rotloc(i,j)=0d0
            do k=1,3
            do l=1,3
               rotloc(i,j)=rotloc(i,j)+as(i,j,k,l)*east(k,l)
            enddo
            enddo
            enddo
            enddo
         endif

C *** CALCULATES VELOCITY GRADIENT IN EACH PHASE AND EACH GRAIN

         do i=1,3
         do j=1,3
            aa(i,j)=ag(i,j,kkk)
            xlijgr0(i,j)=0d0
         enddo
         enddo

         if (verbose)then
            call w_chrc(0,'aa')
            do i=1,3
               write(*,'(3e20.13)')(aa(i,j),j=1,3)
            enddo
            call w_chrc(0,'dnca')
            do is=1,nsysmx
               write(*,'(3e20.13)')(dnca(i,is,iph),i=1,3)
            enddo
            call w_chrc(0,'dbca')
            do is=1,nsysmx
               write(*,'(3e20.13)')(dbca(i,is,iph),i=1,3)
            enddo
            call w_chrc(0,'gamdot')
            do is=1,nsysmx
               write(*,'(e20.13)')gamdot(is,kkk)
            enddo
c            call stop_debug(0,'* stopped inside update_ori_fe')
         endif

         do is=1,nsysmx
            do i=1,3
               dnsa(i)=0d0
               dbsa(i)=0d0
            do j=1,3
               dnsa(i)=dnsa(i)+aa(i,j)*dnca(j,is,iph)
               dbsa(i)=dbsa(i)+aa(i,j)*dbca(j,is,iph)
            enddo
            enddo
            do i=1,3
            do j=1,3
c            XLIJGR0(I,J)=XLIJGR0(I,J)+DBSA(I)*DNSA(J)*GAMDOT(IS,KGX)
c              D^c
               xlijgr0(i,j)=xlijgr0(i,j)+dbsa(i)*dnsa(j)*gamdot(is,kkk)
            enddo
            enddo
         enddo

         if (verbose) then
           call fill_line(0,'-',72)
            call w_chrc(0,'rotloc')
            do i=1,3
               write(*,'(3e20.13)')(rotloc(i,j),j=1,3)
            enddo
            call w_chrc(0,'xlijgr0')
            do i=1,3
               write(*,'(3e20.13)')(xlijgr0(i,j),j=1,3)
            enddo
         endif

         do i=1,3
         do j=1,3
            xlijgrx(i,j)=WBARc(i,j)+rotloc(i,j)+
     $           (xlijgr0(i,j)+xlijgr0(j,i))/2d0
            xlijph(i,j,iph)=xlijph(i,j,iph)+
     $           xlijgrx(i,j)*wgt(kkk)/wph(iph)
         enddo
         enddo

         if (verbose) then
         call w_chrc(0,'xlijph in update_orientation')
         do i=1,3
            write(*,'(3e20.13)')(xlijph(i,j,iph),j=1,3)
         enddo
         endif
cfeb
         if(ishape(iph).gt.0) then
            do i=1,3
            do j=1,3
               xlijgr(i,j,kkk)=xlijgrx(i,j)
            enddo
            enddo
         endif
cfee
c *** crystallographic grain rotation (rigid minus plastic)
         do i=1,3
         do j=1,3
            rotslip(i,j)=(xlijgr0(i,j)-xlijgr0(j,i))/2d0
         enddo
         enddo

         if (verbose)then
            call w_chrc(0,'rotslip')
            do i=1,3
               write(*,'(3e20.13)')(rotslip(i,j),j=1,3)
            enddo
         endif

         do i=1,3
         do j=1,3
c          rot(i,j,kgx)=(WBARc(i,j)+rotloc(i,j)-rotslip(i,j))*tincr
            rot(i,j,kkk)=(WBARc(i,j)+rotloc(i,j)-rotslip(i,j))*dtime
         enddo
         enddo
         if (verbose) then
            call w_val(0,'dtime:',dtime)
            call w_chrc(0,'XLIJGR0')
            do i=1,3
               write(*,'(3e20.13)')(xlijgr0(i,j),j=1,3)
            enddo
            call w_chrc(0,'rotloc')
            do i=1,3
               write(*,'(3e20.13)')(rotloc(i,j),j=1,3)
            enddo
            call w_chrc(0,'rotslip')
            do i=1,3
               write(*,'(3e20.13)')(rotslip(i,j),j=1,3)
            enddo
            do i=1,3
               do j=1,3
                  aux33(i,j)=rot(i,j,kkk)
               enddo
            enddo
            call w_chrc(0,'rot')
            do i=1,3
               write(*,'(3e20.13)')(rot(i,j,kkk),j=1,3)
            enddo
c
         endif

C *** AVERAGE PLASTIC, LOCAL & TOTAL ROTATION (FOR STATISTICAL PURPOSES ONLY)
cc      ROTPLASTAV=ROTPLASTAV+
cc   #    SQRT(ROTSLIP(3,2)**2+ROTSLIP(1,3)**2+ROTSLIP(2,1)**2)*WGT(KKK)
cc      ROTLOCALAV=ROTLOCALAV+
cc   #    SQRT(ROTLOC(3,2)**2+ROTLOC(1,3)**2+ROTLOC(2,1)**2)*WGT(KKK)
cc      ROTTOTALAV=ROTTOTALAV+SQRT(ROT(3,2,KGX)**2+
cc   #    ROT(1,3,KGX)**2+ROT(2,1,KGX)**2)*WGT(KKK)/TINCR

c        KGX=KGX+1
 1000 continue                  ! END OF DO LOOP OVER GRAINS
 2000 continue


c$$$      call w_chrc(0,'rot')
c$$$      do kkk=1,ngr
c$$$         do i=1,3
c$$$            write(*,'(3e20.13)')(rot(i,j,kkk),j=1,3)
c$$$         enddo
c$$$      enddo

c      call stop_debug(0,'** checking rot inside update_ori_fe')

cc    WRITE( *,'('' PLA, LOC, TOT '',3F12.4)')
cc   #              ROTPLASTAV,ROTLOCALAV,ROTTOTALAV

c      KGX=1

      if (verbose) then
         call w_ival(0,'nneigh:',nneigh)
         call w_chrc(0,'wneigh')
         DO KKK=igb,igt
            do in=0,nneigh
               write(*,'(e20.13)')wneigh(in,kkk)
            enddo
         enddo
         call w_chrc(0,'wgt')
         DO KKK=igb,igt
            write(*,'(e20.13)')wgt(kkk)
         enddo
c         call stop_debug(0,'** checking wneight inside update_ori_fe')
      endif

      do kkk=igb,igt
         do i=1,3
         do j=1,3
            aa(i,j)=ag(i,j,kkk)
            corot(i,j)=0d0
c          write(*,*)'nneigh:',nneigh
         do in=0,nneigh
c             write(*,*)'kkk:',kkk,'neigh:',neigh(in,kkk)
c            corot(i,j)=corot(i,j)+rot(i,j,neigh(in,kgx))*wneigh(in,kgx)
c            corot(i,j)=corot(i,j)+rot(i,j,neigh(in,kkk))*wneigh(in,kkk) !!

c            work-around... - neigh is not properly given.
            corot(i,j)=corot(i,j)+rot(i,j,kkk)*1d0
         enddo
         enddo
         enddo

C *** CALCULATE THE NEW TRASFORMATION MATRIX AND UPDATE
         if (verbose) then
            call w_chrc(0,'corot')
            do i=1,3
               write(*,'(3e20.13)')(corot(i,j),j=1,3)
            enddo
         endif

         call rodrigues (corot,arot)

         do i=1,3
            do j=1,3
               ag(i,j,kkk)=0d0
               do k=1,3
                  ag(i,j,kkk)=ag(i,j,kkk)+arot(i,k)*aa(k,j)
               enddo
            enddo
         enddo


c$$$         call w_chrc(0,'ag')
c$$$         do i=1,3
c$$$            write(*,'(3e20.13)')(ag(i,j,kkk),j=1,3)
c$$$         enddo



         if (verbose) then

            write(*,*)'nneigh:',nneigh

            call w_chrc(0,'wgt(kkk)')
            write(*,'(e20.13)')wgt(kkk)

            call w_chrc(0,'arot')
            do i=1,3
               write(*,'(3e20.13)')(arot(i,j),j=1,3)
            enddo

            call w_chrc(0,'aa')
            do i=1,3
               write(*,'(3e20.13)')(aa(i,j),j=1,3)
            enddo

            call w_chrc(0,'ag')
            do i=1,3
               write(*,'(3e20.13)')(ag(i,j,kkk),j=1,3)
            enddo

            call w_chrc(0,'rot')
            do i=1,3
               write(*,'(3e20.13)')(rot(i,j,kkk),j=1,3)
            enddo

            call w_chrc(0,'corot')
            do i=1,3
               write(*,'(3e20.13)')(corot(i,j),j=1,3)
            enddo
            call stop_debug(0,'** checking corot inside update_ori')
         endif


c        kgx=kgx+1
      enddo                     ! end of do loop over grains
c      enddo      ! end of do loop over phases


      return
      end subroutine update_orientation_fe
