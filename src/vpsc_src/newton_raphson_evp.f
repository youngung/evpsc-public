c     newton raphson procedure to solve the grain stress (or grain stress increment) for DELTA and SIGMA EVPSC

c     Argument order: Numerical CONTROL & CONDITION, DIMENSION, PX, PHASE, GRAIN, SYSTEM
c     Then, if applicable, alphabetical order is used.
      subroutine newton_raphson_evp(
c     intent(in)
c          integer in
     $     iregime,ITMX,mdim,
     $     NSYSMX,NSYSTx,ISENSEx,NRSx,
c          floats in
     $     eps,dtime,SG_OLDx,
     $     xmastx,MGELx,f0nr,
     $     CRSSx,GAMD0X,
     $     SCHx,nSCHx,c6,taulim,
c     intent(inout)
c          floats inout
     $     dSGx,SGx,
c     intent(out)
c          integers out
     $     idiag_nr,nriter)
c     --
c     First, check if stress overshoots.
c     Second, calculate the objective function
c     Third, calculate the Jacobian (df/d(dSGx)); fgrad
c     Solve the linear equation and obtain new estimate of X
c     --

c     idiag_nr: diagnosis indicator
c       0: no issue
c       1: encountered nan value from objective array (F) (not used)
c       2: Encountered singularity from <LU_EQSYSTEM>
c       3: Reached max nriteration without satisfying the convergence crnriterion.
      implicit none
c     --
      integer, intent(in):: iregime,ITMX,mdim,NSYSTx,NSYSMX,
     $     ISENSEx(NSYSMX),NRSx(NSYSMX)
      real*8, intent(in):: eps,dtime,SG_OLDx(6),
     $     xmastx(mdim,mdim),MGELx(6,6),f0nr(6),
     $     CRSSx(NSYSMX),GAMD0X(NSYSMX),
     $     SCHx(5,NSYSMX),nSCHx(6,NSYSMX),
     $     c6,taulim
      real*8, intent(inout)::dSGx(6),SGx(6)
      integer, intent(out):: idiag_nr,nriter

c     -- local variables
      integer is,i,j,k,isingular
      real*8 fgrad(6,6),f(6),fold(6),dSGxold(6),SGxold(6),rss(NSYSTx),
     $     gd(NSYSTx),gdstar(NSYSTx),sfactor,relerror,coef,
     $     dsgx0(6),sgx0(6),
     $     hyd,signhyd,vmismatch,mgsecvp(6,6),mgstar(6,6,6),
     $     aux5(5),aux6(6),aux33(3,3),aux3333(3,3,3,3),aux66(6,6),
     $     aux55(5,5)
      logical isnan
c     ------------------------------------------------
c     arrays may not be used now but can be used later
      aux33(:,:)=0.
      aux3333(:,:,:,:)=0.
c      aux5(5)=aux5(5)*0.
c     ------------------------------------------------

      dsgx0(:)=0.
      sgx0(:)=0.
      dsgx0(1:mdim)=dsgx(1:mdim)
      sgx0(1:mdim)=sgx(1:mdim)

      coef=0.2
      SGXOLD(:)=0.
      SGXOLD(1:MDIM)=SGX(1:MDIM)

      isingular=0
      idiag_nr=0


c     ! in case stress overshoot the yield surface
c     ! the trial stress is used to reestimate the
c     ! trial 'incremental' stress

c      dtau(:)=0.
      do 1000 nriter=1, ITMX
         if (iregime.eq.2) then !! Delta-EVPSC
            aux6(:)=dSGx(:)+SG_OLDx(:)
         elseif(iregime.eq.1.or.iregime.eq.3) then !! SIGMA-EVPSC
            aux6(:)=SGx(:)
         endif

         do i=1,mdim
            if (isnan(aux6(i))) then
c$$$               call mutexlock(2)
c$$$               write(*,*)'aux6 has nan value in it'
c$$$               write(*,*)'SG_OLDx:',SG_OLDx
c$$$               write(*,*)'dSGx:',dSGx
c$$$               write(*,*)'SGx:',SGx
c$$$               call mutexunlock(2)
c$$$c               call xit()
               idiag_nr=2
               return
            endif
         enddo

         hyd=aux6(6)/3.
         signhyd=sign(1.,hyd)

C *** RESOLVED SHEARS CALCULATION & OUTSIDE Y.S. ERROR MANAGING BLOCK.
C *** NRSx MAY BE EVEN OR ODD:
C     GD ALWAYS > 0 --> GD IS USED TO GET DERIVATIVES TO BUILD THE
C     COEFFICIENT MATRIX FOR N-R METHOD IN F(I) CALCULATION (INDEPENDENT
C     TERM FOR N-R METHOD)   RSS*GD=GAMDOT -> SIGN(RSS*GD)=SIGN(RSS)

         do is=1, NSYSTx
            rss(is)=0.
            do i=1,5
               rss(is)=rss(is)+(SCHx(i,is)+nSCHx(i,is))*aux6(i)
            enddo
            rss(is)=rss(is)+nSCHx(6,is)*aux6(6)+c6*abs(hyd)


            if (.not.(rss(is).gt.0..or.ISENSEx(is).eq.1))rss(is)=0.
            rss(is)=rss(is)/CRSSx(is)

c         yj remarks: I am not sure of these corrections in NR.
c$$$          if(ihardlaw.ne.31 .and. ihardlaw.ne.32.and.mdim.eq.5) then     ! AP:
            if (iregime.eq.1) then
               IF(ABS(RSS(IS)).GT.TAULIM) THEN
c$$$                  call mutexlock(2)
c$$$                  write(*,*)'nriter:',nriter,'taulim:',taulim
c$$$                  write(*,*)'warning: rss is exceedingly large'
c$$$                  write(*,'(a6,6e10.2)')'rss(1:6)', (rss(i),i=1,6)
c$$$                  write(*,'(a6,6e10.2)')'dSGx:',dSGx
c$$$                  write(*,'(a6,6e10.2)')'SGx:',SGx
c$$$                  call mutexunlock(2)
                  DO I=1,5
                     SGX(I)=SGXOLD(I)+COEF*(SGX(I)-SGXOLD(I))
                  ENDDO
c                  call xit()
                  GO TO 1000
               ENDIF
            endif

            if(abs(rss(is)).lt.1.e-10) rss(is)=0.
            gd(is) = GAMD0X(is)
            if(NRSx(is).ne.1) gd(is)=gd(is)*abs(rss(is)**(NRSx(is)-1))
            gdstar(is)=gamd0x(is)         *abs(rss(is)**(nrsx(is)-2))*
     $           sign(1.,rss(is))
         enddo                  ! loop over system


         mgsecvp(:,:)=0d0
         do i=1,5
         do j=1,5
            do is=1,NSYSTx
               mgsecvp(i,j)=mgsecvp(i,j)+SCHx(i,is)
     $              *SCHx(j,is)/CRSSx(is)*gd(is)
            enddo
         enddo
         enddo
         mgstar(:,:,:)=0.
         do i=1,5               ! multiplied with strain rate
         do j=1,5
         do k=1,5
            do is=1,NSYSTx
               mgstar(i,j,k)=mgstar(i,j,k)+
     $              (nrsx(is)-1)*gdstar(is)/((CRSSx(is))**2)
     $              *SCHx(i,is)*SCHx(j,is)*SCHx(k,is)
            enddo
         enddo
         enddo
         enddo


         do i=1,NSYSTx
            if (isnan(rss(i))) then
c$$$               call mutexlock(2)
c$$$               write(*,*)'rss(i) has nan value in it'
c$$$               write(*,'(a6,6e10.2)')'rss:',(rss(j),j=1,mdim)
c$$$               call mutexunlock(2)
c$$$c               call xit()
               idiag_nr=2
               return
            endif
         enddo


c***  The objective function
         f(:)=0.
         do i=1,mdim
            f(i)=f0nr(i)
         enddo

c          ** Define objective function array
c$$$         if (iregime.eq.2) then ! Delta
c$$$            do i=1,6
c$$$               do j=1,6
c$$$                  f(i)=f(i)+(MGELx(i,j)/dtime+mast(i,j))*dSGx(j)
c$$$               enddo
c$$$            enddo
c$$$         elseif(iregime.eq.3) then ! Sigma
c$$$            do i=1,6
c$$$               do j=1,6
c$$$                  f(i)=f(i)+(MGELx(i,j)/dtime+mast(i,j))*SGx(j)
c$$$               enddo
c$$$            enddo
c$$$         endif

        if (mdim.eq.6) then
           if (iregime.eq.2) then
              do i=1,mdim
              do j=1,mdim
                 F(i)=F(i)+MGELx(i,j)*dsgx(j)/dtime
              enddo
              enddo
c$$$           elseif(iregime.eq.3)then
c$$$              do i=1,mdim
c$$$              do j=1,mdim
c$$$                 F(i)=F(i)+MGELx(i,j)*sgx(j)/dtime
c$$$              enddo
c$$$              enddo
           endif
        endif

        DO I=1,5
           DO IS=1,NSYSTX
              F(I)=F(I)+SCHX(I,IS)*RSS(IS)*GD(IS)
           ENDDO
        ENDDO

        if (iregime.eq.2) then
           DO I=1,MDIM
           DO J=1,MDIM
              F(I)=F(I)+XMASTX(I,J)*dSGX(J) ! M~:dSc
           ENDDO
           ENDDO
        elseif (iregime.eq.1.or.iregime.eq.3) then
           DO I=1,MDIM
           DO J=1,MDIM
              F(I)=F(I)+XMASTX(I,J)*SGX(J) ! M~:Sc
           ENDDO
           ENDDO
        endif

c        write(*,'(a10,6f10.4)')'F:',F(:)
c        call xit()


         do i=1,mdim
            if (isnan(f(i)).or.isnan(f0nr(i))) then
c$$$               call mutexlock(2)
c$$$               write(*,*)'nriter:',nriter
c$$$               write(*,*)'nan value at line 230'
c$$$               write(*,'(a6,6e10.2)')'f',(f(j),j=1,mdim)
c$$$               write(*,'(a6,6e10.2)')'f0nr',(f0nr(j),j=1,mdim)
c$$$               write(*,'(5a10)')'is','CRSSx','rss','gd','SCHx(...)'
c$$$               do is=1,NSYSTx
c$$$                  write(*,'(i10,8e10.2)')is,CRSSx(is),rss(is),gd(is),
c$$$     $                 (SCHx(j,is),j=1,5)
c$$$               enddo
c$$$               write(*,'(a6,6e10.2)')'dSGx:',dSGx
c$$$               write(*,'(a6,6e10.2)')'SGx:',SGx
c$$$               write(*,'(a6,6e10.2)')'dSGx0:',dSGx0
c$$$               write(*,'(a6,6e10.2)')'SGx0:',SGx0
c$$$               write(*,*)'dtime:',dtime
c$$$               call w_chr(0,'xmast')
c$$$               do k=1,mdim
c$$$                  write(*,'(6e10.2)')(xmastx(k,j),j=1,mdim)
c$$$               enddo
c$$$               call mutexunlock(2)
c$$$c               call xit()
               idiag_nr=2
               return
            endif
         enddo

c          ** Define gradient array
         fgrad(:,:)=0d0
         if (.true.) then
c         if (.false.) then
            do is=1,NSYSTx
               sfactor=NRSx(is)*gd(is)/CRSSx(is)
               do i=1,5
               do j=1,5
                  fgrad(i,j)=fgrad(i,j)
     $                 +sfactor*SCHx(i,is)*(SCHx(j,is)+nSCHx(j,is))
               enddo
               fgrad(i,6)=fgrad(i,6)+sfactor*SCHx(i,is)
     $              *(nSCHx(6,is)+c6/3.*signhyd)
               enddo
            enddo
         else
            do i=1,5
            do j=1,5
               do k=1,5
                  fgrad(i,j)=fgrad(i,j)+mgstar(i,k,j)*aux6(k)
               enddo
               fgrad(i,j)=fgrad(i,j)+mgsecvp(i,j)
            enddo
            enddo
         endif

         DO I=1,MDIM
         DO J=1,MDIM
            FGRAD(I,J)=fgrad(i,j)+XMASTX(I,J)
         ENDDO
         ENDDO

c       contribution from elasticity
        if (iregime.eq.2.or.iregime.eq.3) then
           do i=1,mdim
              do j=1,mdim
                 fgrad(i,j)=fgrad(i,j)+MGELx(i,j)/dtime
              enddo
           enddo
        endif






C *** SOLVES MDIM X MDIM LINEAR SYSTEM TO GET DEVIATORIC STRESS COMPONENTS
C          F(SG)=F(SG0)+FGRAD(SG0)*(SG-SG0) = 0
C          (SG-SG0) = (-FGRAD(SG0))^(-1)  F(SG0)

c      in order to observe how convergence is made
         dSGxold(:)=dSGx(:)
         fold(:)=f(:)

         SGxold(:)=SGx(:)
         do i=1,6
         do j=1,6
            FGRAD(i,j) = -FGRAD(i,j)
         enddo
         enddo
         aux66(:,:)=0.
         aux66(1:mdim,1:mdim)=fgrad(1:mdim,1:mdim)

c        solves jacob_ij X_j=f_i for X
         if (iregime.eq.1) then ! with mdim being 5
            do i=1,mdim
               aux5(i)=f(i)
            do j=1,mdim
               aux55(i,j)=fgrad(i,j)
            enddo
            enddo

            call lu_eqsystem(aux55,aux5,mdim,isingular)
c$$$            write(*,*)'iregime,mdim:',iregime,mdim
c$$$            write(*,'(a3,i5,a10,a50,a10)')'nriter:',nriter,
c$$$     $           'f','fgrad','f after'
c$$$            do i=1,mdim
c$$$               write(*,'(8x,e10.2,x,5e10.2,x,e10.2)')
c$$$     $              f(i),(fgrad(i,j),j=1,mdim),aux5(i)
c$$$            enddo
c$$$            write(*,'(a6,6e10.2)')'dSGx:',dSGx
c$$$            write(*,'(a6,6e10.2)')'SGx:',SGx
c$$$            write(*,'(a6,6e10.2)')'dSGx0:',dSGx0
c$$$            write(*,'(a6,6e10.2)')'SGx0:',SGx0
c$$$            write(*,'(a6)',advance='no')'rss:'
c$$$            do is=1,NSYSTx
c$$$               write(*,'(e10.2)',advance='no')rss(is)
c$$$            enddo
c$$$            write(*,*)
c$$$            write(*,'(a6)',advance='no')'crssx:'
c$$$            do is=1,NSYSTx
c$$$               write(*,'(e10.2)',advance='no')crssx(is)
c$$$            enddo
c$$$            write(*,*)
c$$$            write(*,'(a6)',advance='no')'gd'
c$$$            do is=1,NSYSTx
c$$$               write(*,'(e10.2)',advance='no')gd(is)
c$$$            enddo
c$$$            write(*,*)
c$$$            write(*,*)'--'
c            call xit()
            fgrad(1:mdim,1:mdim)=aux55(:,:)
            f(1:mdim)=aux5(:)
         else
            call lu_eqsystem(fgrad,f,mdim,isingular)
         endif

         if(isingular.eq.1) then
c$$$            call mutexlock(2)
c$$$            write(*,*)'------------------------------------------'
c$$$            write(*,*)'** Error in NEWTON_RAPHSON_EVP subroutine.'
c$$$            write(*,*)'   Singularity encoutered when solving '
c$$$            write(*,*)'    system of equations via LU_EQSYSTEM in '
c$$$            write(*,*)'    subroutine NEWTON_RAPHSON_EVP'
c$$$            write(*,*)'------------------------------------------'
c$$$            write(*,*)'nriter:',nriter
c$$$            write(*,'(a10,6e10.2)')'SGx:',SGx
c$$$            write(*,'(a10,6e10.2)')'dSGx:',dSGx
c$$$            write(*,'(a10,6e10.2)')'dSGxold:',dSGxold
c$$$            write(*,'(a10,6e10.2)')'SGxold:',SGxold
c$$$            write(*,'(a10,6e10.2)')'fold:',fold
c$$$            write(*,'(a10,6e10.2)')'f0nr:', f0nr
c$$$
c$$$c     $     SCHx(5,NSYSMX),nSCHx(6,NSYSMX),
c$$$            write(*,'(a10,5e10.2)')'schx(:,1):',SCHx(:,1)
c$$$            write(*,'(a10,6e10.2)')'gamd0x(1:6)', (gamd0x(is),is=1,6)
c$$$            write(*,'(a10,6e10.2)')'gd(1:6)', (gd(is),is=1,6)
c$$$            write(*,'(a10,6e10.2)')'rss(1:6)', (rss(is),is=1,6)
c$$$            write(*,*)
c$$$            write( *,'('' f    ='',6e12.3)') f
c$$$            call w_chr(0,'fgrad')
c$$$            do i=1,mdim
c$$$               write(*,'(6e10.2)')(fgrad(i,j),j=1,mdim)
c$$$            enddo
c$$$            write(*,*)'dtime:',dtime
c$$$            call w_chr(0,'xmast')
c$$$            do i=1,mdim
c$$$               write(*,'(6e10.2)')(xmastx(i,j),j=1,mdim)
c$$$            enddo
c$$$            call mutexunlock(2)
c$$$c            call xit()
            idiag_nr=2
            return
         endif

         do i=1,mdim
            if (isnan(f(i))) then
c$$$               call mutexlock(2)
c$$$               write(*,*)'f(i) has nan value in it'
c$$$               write(*,'(a6,6e10.2)')'f:',(f(j),j=1,mdim)
c$$$               write(*,'(a6,6e10.2)')'fold:',(fold(j),j=1,mdim)
c$$$               write(*,*)'fgrad'
c$$$               do j=1,mdim
c$$$                  write(*,'(6e10.2)')(fgrad(j,k),k=1,mdim)
c$$$               enddo
c$$$               write(*,*)'fgrad_old'
c$$$               do j=1,mdim
c$$$                  write(*,'(6e10.2)')(aux66(j,k),k=1,mdim)
c$$$               enddo
c$$$               write(*,'(a6,6e10.2)')'dSGx0:',dSGx0
c$$$               write(*,'(a6,6e10.2)')'SGx0:',SGx0
c$$$               write(*,*)'dtime:',dtime
c$$$               call mutexunlock(2)
               idiag_nr=2
               return
            endif
         enddo


         if(iregime.eq.1) then
            SGx(:)=SGxold(:)+f(:)
            relerror = vmismatch(SGx,SGxold,6)
         elseif (iregime.eq.2) then
            dSGx(:)=dSGxold(:)+f(:) ! xnew = xold - f*fgrad^-1
            relerror = vmismatch(dSGx,dSGxold,6)

         endif

         if (relerror.lt.eps) then
            relerror=0d0
            do i=1,mdim
               relerror=relerror+f(i)**2
            enddo
            relerror=sqrt(relerror)
c            write(*,'(a10,f10.7)')'relerror:',relerror
c            if (relerror.lt.eps) then
c               call chg_basis(dsgx,aux33,aux66,aux3333,1,6)
c               call voigt(aux6,aux33,aux66,aux3333,2)
c               write(*,'(1(a15,6f9.3))')'dsg at NR:',aux6(:)
c               return
c            endif
            return
         endif
 1000 continue                  ! end of big loop

c$$$      write(*,*)'could not converge within the big loop'
c$$$      write(*,'(a6,6e10.2)')'sgx',sgx
c$$$      write(*,'(a6,6e10.2)')'f:',(f(j),j=1,mdim)
c$$$      write(*,*)'fgrad'
c$$$      do j=1,mdim
c$$$         write(*,'(6e10.2)')(fgrad(j,k),k=1,mdim)
c$$$      enddo

      idiag_nr=3                  ! in case NR did not converge.
      return
      end subroutine newton_raphson_evp
