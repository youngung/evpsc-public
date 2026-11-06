C *****************************************************************************
C     SUBROUTINE GRAIN_STRESS_EVP
C
C     GIVEN A GUESS increment (or non-incremental) STRESS 'dX(I)' AND GRAIN
C     INDEX 'KKK', SOLVES INTERACTION EQUATION TO FIND THE increment
C     (or nonincremental) STRESS 'X' COMPATIBLE WITH MACROSCOPIC COMPLIANCE.
C     Solves INTERACTION EQUATION FOR A SELF-CONSISTENT CASE.
c      incremental-EVPSC
c        (dg-DBAR) = -M~:(dSG - dSBAR)
c     nonincremental-EVPSC
c        (dg-DBAR) = -M~:(SG  -  SBAR)
C *****************************************************************************

      subroutine grain_stress_evp(
c          in
     $     interaction,mdim,
     $     iregime,ISHAPEx,ISENSEx,NRSx,
     $     NSYSTx,NSYSMX,

     $     DASDBAR,SASTAV,

     $     DBAR,DBAR0,DG_TRANSx,dtime,
     $     SBAR,SBAR_OLD,
     $     SG_OLDx,

     $     MASTPHx,MASTGRx,

     $     SCHx,NSCHx,c6,
     $     CRSSx,GAMD0Sx,MGELx,
     $     dSGTRYx,SGTRYx,
c          out
     $     dSG,SG,nriter,idiag_nr,idiag_SG)

      implicit none
      integer,intent(in)::interaction,mdim,NSYSMX,
     $     iregime,ISHAPEx,ISENSEx(NSYSMX),NRSx(NSYSMX),
     $     NSYSTx
      real*8,intent(in)::SBAR(6),SBAR_OLD(6),
     $     DBAR(6),DBAR0(6),DASDBAR(mdim),SASTAV(mdim),
     $     DG_TRANSx(6),
     $     SG_OLDx(6),MASTPHx(mdim,mdim),MASTGRx(mdim,mdim),
     $     SCHx(5,NSYSMX),NSCHx(6,NSYSMX),c6,CRSSx(NSYSMX),
     $     GAMD0Sx(NSYSMX),MGELx(6,6),dSGTRYx(6),SGTRYx(6)
      integer, intent(out)::nriter,idiag_nr,idiag_SG
      real*8,  intent(out)::dSG(6),SG(6)
c     -----------------
      real*8 taumax,rssx,f0nr(6),xmastx(mdim,mdim),
     $     dtime,xori(6),x(6),dx(6),
     $     SGX(6),eps,hyd,taulim,vnorm
      integer is,i,j,ITMX
      parameter(ITMX=1000,eps=5e-04)


      idiag_SG=0

C     EMPIRIC ALGORITHM TO GET TAULIM FOR NR SUBROUTINE (RL: 01/FEB/00)
      if (gamd0sx(1).eq.0) then
         write(*,*)'gamd0sx of zero found'
         call xit()
         taulim=2.
      else
         taulim=2.*(Vnorm(DBAR,mdim)/GAMD0Sx(1))**0.05 ! yj-fix 5->mdim
      endif

      if(taulim.lt.2.) taulim=2.

c     -----------------------------------------------------------------
c     ** determine stress (or stress increment)
c     ** corrects delta stress 'x' if it exceeds the yield
c        surface (not actually performed in the case of EVPSC)
      if(iregime.eq.1) then
         SGX(:)=SGTRYx(:)
      elseif (iregime.eq.2) then
         SGX(:)=dSGTRYx(:)+SG_OLDx(:)
      endif

      hyd=SGX(6)/3.


C *** CORRECTS STRESS 'X' IF IT EXCEEDS THE YIELD SURFACE
c     - This correction is performed only for VPSC
      taumax=0.                 ! to fool
      do is=1,NSYSTx
         rssx=0d0
         do i=1,5
            rssx=rssx+(SCHx(i,is)+NSCHx(i,is))*SGX(i)
         enddo
         rssx=rssx+NSCHx(6,is)*SGX(6)
         rssx=rssx+c6*abs(hyd)
         if(.not.(rssx.gt.0.or.ISENSEx(is).eq.1)) rssx=0d0
         rssx=rssx/CRSSx(is)
         if (abs(rssx).gt.taumax) taumax=abs(rssx)
      enddo


c       tolerance might need re-adjustment.
      if (taumax .lt. 1.d-50.and.iregime.eq.1) then
         write(*,'(''taumax:'',6e12.4)')taumax
         write(*,'(''taumax<1e-20 inside subroutine grain_stress'')')
         write(*,'(''SGx  ='',6e12.4)')(SGx(i),i=1,6)
         write(*,'(''dSGTRYx  ='',6e12.4)')(dSGTRYx(i),i=1,6)
         write(*,'(''SGTRYx  ='',6e12.4)')(SGTRYx(i),i=1,6)
         write(*,'(''SG_OLDx='',6e12.4)')(SG_OLDx(i),i=1,6)
         write(*,'(''CRSSx  ='',24e12.4)')(CRSSx(is),is=1,NSYSTx)
         write(*,'(''SCHx(1)='',24e12.4)')(SCHx(1,is),is=1,NSYSTx)
         call xit()
      endif



      IF(TAUMAX.GT.TAULIM .OR. INTERACTION.EQ.0) THEN
         if (iregime.eq.1) then
            SGX(:)=SGX(:)/TAUMAX
         endif
      ENDIF

      IF(INTERACTION.EQ.0) THEN     ! interaction=-1 was used for RC in the past!
         DO I=1,mdim            ! yj-fix mdim
            F0NR(I)=DG_TRANSX(I)-DBAR(I)
            DO J=1,mdim         ! yj-fix mdim
               XMASTX(I,J)=0.   ! allows to use the same algorithm in Newt-Raph for FC & SC
            ENDDO
         ENDDO
      ELSE IF(INTERACTION.GT.0) THEN
          IF(ISHAPEX.LE.1) THEN
             DO I=1,MDIM        ! yj-fix mdim
             DO J=1,MDIM     ! yj-fix mdim
                XMASTX(I,J)=MASTPHX(I,J) ! this is M~ calc in VPSC & includes neff factor
             ENDDO
             ENDDO
        ELSE IF(ISHAPEX.GT.1) THEN
          DO I=1,MDIM ! yj-fix mdim
          DO J=1,MDIM ! yj-fix mdim
            XMASTX(I,J)=MASTGRX(I,J)      ! this is M~ calc in VPSC & includes neff factor
          ENDDO
          ENDDO
        ENDIF

        if (MDIM.eq.5) then !! VPSC
           DO I=1,MDIM
              F0NR(I)=DG_TRANSX(I)-(DASDBAR(I)+DBAR0(I)) ! DASDBAR=MBARTG*SASDBAR, calc inside VPSC
              DO J=1,MDIM
                 F0NR(I)=F0NR(I)-XMASTX(I,J)*SASTAV(J)
              ENDDO
           ENDDO
        elseif (MDIM.eq.6.and.iregime.eq.2) THEN !delta-EVPSC
           do i=1, MDIM
              F0NR(I)=DG_TRANSX(I)-DBAR(I)
              do j=1,MDIM
                 F0NR(I)=F0NR(i)-XMASTX(I,J)*(SBAR(j)-SBAR_OLD(j))
              enddo
           enddo
        elseif (MDIM.eq.6.and.iregime.eq.3) THEN !sigma-EVPSC
           do i=1, MDIM
              F0NR(I)=DG_TRANSX(I)-DBAR(I)
              do j=1,MDIM
                 F0NR(I)=F0NR(i)-XMASTX(I,J)*SBAR(j)
     $                -MGELx(i,j)*SG_OLDx(j)/dtime
              enddo
           enddo
        endif

      endif

c     -----------------------------------------------------------------
c     ** Call Newton-Raphson-EVP subroutine

c     before calling NR, store initial trial stress (or stress increment)
c     If NR fails, this tried strss (or stress increment) will be returned
c     from subroutine grain_stress
      if(iregime.eq.1) then
         xori(:)=SGx(:)
         x(:)=SGx(:)
         dx(:)=0d0 !! irrelevant
      elseif (iregime.eq.2) then
         xori(:)=dSGTRYx(:)
         dx(:)=dSGTRYx(:)
         x(:)=0d0 !! irrelevant
      endif

c     ------------------------------------------------------------------------
c     ** solve the non-linear equation to determine grain stress
c     calling newton_raphson_evp
      call newton_raphson_evp(
c     intent(in)
c          integer in
     $     iregime,ITMX,mdim,
     $     NSYSMX,NSYSTx,ISENSEx,NRSx,
c          floats in
     $     eps,dtime,SG_OLDx,
     $     xmastx,MGELx,f0nr,
     $     CRSSx,GAMD0Sx,
     $     SCHx,nSCHx,c6,taulim,
c     intent(inout)
c          floats inout
     $     dx,x,
c     intent(out)
c          integers out
     $     idiag_nr,nriter)
c     idiag_nr: diagnosis indicator
c       0: no issue
c       1: encountered nan value from objective array (F - not used)
c       2: Encountered singularity from <LU_EQSYSTEM>
c       3: Reached max iteration without satisfying the convergence criterion.

      if (idiag_SG.gt.0.or.idiag_nr.gt.0)then
c$$$         write(*,*)'idiag_SG:',idiag_SG
c$$$         write(*,*)'idiag_nr',idiag_nr
         if (iregime.eq.1)  then
            x(:)=xori(:)
         elseif (iregime.eq.2)  then
            dx(:)=xori(:)
         endif
c$$$         write(*,*)'could not converge'
c$$$         write(*,'(''SGx  ='',6e12.4)')(SGx(i),i=1,6)
c$$$         write(*,'(''dSGTRYx  ='',6e12.4)')(dSGTRYx(i),i=1,6)
c$$$         write(*,'(''SGTRYx  ='',6e12.4)')(SGTRYx(i),i=1,6)
c$$$         write(*,'(''SG_OLDx='',6e12.4)')(SG_OLDx(i),i=1,6)
c$$$         write(*,'(''CRSSx  ='',24e12.4)')(CRSSx(is),is=1,NSYSTx)
c$$$         write(*,'(''SCHx(1)='',24e12.4)')(SCHx(1,is),is=1,NSYSTx)
c         call xit()
         return
      endif

c     ------------------------------------------------------------------------
c     **  update stress by the one obtained from newton raphson iteration
      if (iregime.eq.1) then
         SG(:)=x(:)
         dSG(:)=SG(:)-SG_OLDx(:)
      elseif (iregime.eq.2) then
         dSG(:)=dx(:)
         SG(:)=SG_OLDx(:)+dSG(:)
      endif

      return
      end subroutine grain_stress_evp
