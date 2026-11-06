      SUBROUTINE STATE_NxN (ibc_d,BC_D,ibc_s,BC_S,XMS,ndim)
      implicit none
      integer, intent(in):: ndim,ibc_d(ndim),ibc_s(ndim)
      real*8, intent(inout):: bc_d(ndim),bc_s(ndim)
c     locals
      real*8 aux11(ndim),aux21(ndim,ndim),XMS(ndim,ndim)
      integer i,j,ier

      do i=1,ndim
        aux11(i)=-1.d0*ibc_d(i)*bc_d(i)
        do j=1,ndim
          aux11(i)=aux11(i)+xms(i,j)*ibc_s(j)*bc_s(j)
          aux21(i,j)=ibc_s(j)*(i/j)*(j/i)-ibc_d(j)*xms(i,j)
        enddo
      enddo

      CALL LU_EQSYSTEM(AUX21,AUX11,ndim,IER)

      if(ier.eq.1) then
      write(*,*) 'SINGULAR SYSTEM IN STATE_5X5'
      endif

      do i=1,ndim
        bc_d(i)=ibc_d(i)*bc_d(i)+ibc_s(i)*aux11(i)
        bc_s(i)=ibc_s(i)*bc_s(i)+ibc_d(i)*aux11(i)
      enddo

      return
      end

C
C *****************************************************************************
C     SUBROUTINE STATE_6x6      --->      VERSION 31/AUG/2017

C     SOLVES A 6x6 SYSTEMS FOR CARTESIAN STRESS-STRAIN COMPONENTS EXPRESSED
C     IN VOIGT NOTATION

C *** BASED ON SECANT EQUATION: D(i)-D0(i)=MS(i,j)*S(j), WHERE MS IS COMPLIANCE,
C     SHIFTS KNOWN COMPONENTS TO INDEPENDENT TERM, CALCULATES COEFFICIENTS
C     FOR UNKNOWN COMPONENTS AND INVERTS THE LINEAR SYSTEM.
C *** COMPONENTS OD 'D', 'S' and 'MS' ENTER IN CARTESIAN NOTATION BUT GET
C     TRANSFORMED INTO 6x1 OR 6x6 VOIGT MATRICES.

C *** CALLED AT THE END OF EACH INCREMENTAL STEP IN SUBROUTINE VPSC IF
C     IDBARv1*IDBARv2*IDBARv3=0 (NON FULLY ENFORCED DIAGONAL)
C *****************************************************************************
c
      subroutine state_6x6 (ibc_d,AUX_D,ibc_s,AUX_S,AUX_MS)
      implicit none
      integer, intent(in):: ibc_d(6),ibc_s(6)
      real*8,intent(inout):: AUX_D(3,3),AUX_S(3,3)
      real*8, intent(in):: AUX_MS(3,3,3,3)
c     locals
      real*8 aux6v(6),aux66v(6,6),bc_s(6),bc_d(6),profac(6),
     $     aux33(3,3),aux3333(3,3,3,3),aux66(6,6)
      integer i,j,ier

      do i=1,6
        profac(i)=1.d0+(i/4)
      enddo

C *** CONVERTS CARTESIAN 'D','S' and 'MS' TO VOIGT.
      CALL VOIGT (BC_D,AUX_D,AUX66,AUX3333,2)
      CALL VOIGT (BC_S,AUX_S,AUX66,AUX3333,2)
      CALL VOIGT (AUX6v,AUX33,AUX66,AUX_MS,4)

      do i=1,6
        aux6v(i)=-1.d0*ibc_d(i)*bc_d(i)
        do j=1,6
          aux6v(i)=aux6v(i)+aux66(i,j)*ibc_s(j)*bc_s(j)*profac(j)
          aux66v(i,j)=ibc_s(j)*(i/j)*(j/i)-ibc_d(j)*aux66(i,j)*profac(j)
        enddo
      enddo

      CALL LU_EQSYSTEM(AUX66v,AUX6v,6,IER)

      if(ier.eq.1) then
        write(*,*) 'SINGULAR SYSTEM IN STATE_6X6'
        call xit()
      endif

      do i=1,6
        bc_d(i)=ibc_d(i)*bc_d(i)+ibc_s(i)*aux6v(i)
        bc_s(i)=ibc_s(i)*bc_s(i)+ibc_d(i)*aux6v(i)
      enddo

C *** CONVERTS VOIGT 'D','S' TO CARTESIAN
      CALL VOIGT (BC_D,AUX_D,AUX66,AUX3333,1)
      CALL VOIGT (BC_S,AUX_S,AUX66,AUX3333,1)

      return
      end

c-----------------------------------------------------------
c     subroutine SOLVE_BC_EVP
c     This subrotine solves the macroscopic constitutive law
c     that is given as:
c     E=Mevp:dS+E0
c-----------------------------------------------------------
c     tincr   : time increment
c     iDBARv and iSBARv: flag to if the components are imposed (1) or not (0) - following voigt notation
c     DBAR    : total strain rate in B-basis (6D) - elasto-visco-plastic strain rate
c     DBAR0   : back-extrapolated total strain rate in B-basis
c     mbar  : macro elasto-visco-plastic compliance
c     dsbar   : macro delta stress dS
c     verbose : make the run verbose
c-----------------------------------------------------------
c     Version 20191108 yj

      subroutine solve_bc_evp(iDBARv,iSBARv,DBAR0,mbar,DBAR,dsbar)
      implicit none
      real*8,  intent(in):: DBAR0(6),mbar(6,6)
      real*8,intent(inout):: DBAR(6),dsbar(6)
      integer, intent(in):: iDBARv(6),iSBARv(6)
c     locals
      real*8 dbaux(6),dbauxc(3,3),aux6(6),aux33(3,3),aux66(6,6),
     $     aux3333(3,3,3,3),sbarc(3,3)!,sdev(6)
      logical verbose
c      parameter(verbose=.true.)
      parameter(verbose=.false.)

      if (verbose)then
         call w_chrc(0,'* inside subroutine solve_bc_evp')
      endif

c     need to convert dsbar(b-basis) to sbarc (cartesian)
      call chg_basis(dsbar,sbarc,aux66,aux3333,1,6)

c     solve the macroscopic constitutive equation
      dbaux(:)=DBAR(:)-DBAR0(:)
      call chg_basis(dbaux,dbauxc,aux66 ,aux3333,1,6)
      CALL CHG_BASIS(AUX6 ,AUX33 ,mbar,AUX3333,3,6)

      if (verbose) then
         call w_chrc(0,'iDBARv')
         write(*,'(6i4)')iDBARv
         call w_chrc(0,'iSBARv')
         write(*,'(6i4)') iSBARv

         call w_chrc(0,'sbarc')
         call w_mdim(0,sbarc,3,1d0)
         write(*,*)'aux3333(1,1,1,1)',aux3333(1,1,1,1)

         call w_chrc(0,'dbauxc')
         call w_mdim(0,dbauxc,3,1d0)
      endif

      call state_6x6_cnt(iDBARv,dbauxc,iSBARv,sbarc,aux3333)

      call chg_basis(dsbar,sbarc,aux66,aux3333,2,6)
c      call chg_basis(dsbar,sdev,aux66,aux3333,1,5)
      call chg_basis(dbaux,dbauxc,aux66,aux3333,2,6)
      DBAR(:)=dbaux(:)+DBAR0(:)

      end subroutine solve_bc_evp

C *****************************************************************************
C     SUBROUTINE STATE_6x6      --->      VERSION 31/AUG/2017
C
C *** BASED ON SECANT EQUATION: D(i)-D0(i)=MS(i,j)*S(j), WHERE MS IS COMPLIANCE,
C     SHIFTS KNOWN COMPONENTS TO INDEPENDENT TERM, CALCULATES COEFFICIENTS
C     FOR UNKNOWN COMPONENTS AND INVERTS THE LINEAR SYSTEM ASSUMING VOIGT
C *** 'D', 'S' and 'MS' ENTER IN CARTESIAN NOTATION.

C *** CURRENTLY CALLED AT THE END OF EACH INCREMENTAL STEP IN SUBROUTINE VPSC
C     IF IDBAR1*IDBAR2*IDBAR3=0 (NON FULLY ENFORCED DIAGONAL)
C *****************************************************************************
c
      subroutine state_6x6_cnt (ibc_d,AUX_D,ibc_s,AUX_S,AUX_MS)
      implicit none
c
      integer, intent(in):: ibc_d(6),ibc_s(6)
      real*8, intent(inout)::aux_d(3,3),aux_s(3,3)
      real*8 aux_ms(3,3,3,3)

      real*8 aux6v(6),aux66v(6,6)
c      real*8 AUX_D(3,3),AUX_S(3,3),AUX_MS(3,3,3,3)
      real*8 bc_s(6),bc_d(6)

      real*8 profac(6)
      real*8 aux33(3,3),aux3333(3,3,3,3),aux66(6,6)
      real*8 dum
      integer i,j,ier
      logical verbose,is_all_zero
c      parameter(verbose=.true.)
      parameter(verbose=.false.)
      integer kro(6,6)
      kro(:,:)=0
      do i=1,6
         kro(i,i)=1
      enddo
c
      do i=1,6
         profac(i)=1d0+(i/4)
      enddo

      if(verbose) then
         write(*,'(6f5.2)') (profac(i),i=1,6)
         write(*,*)'kronecker delta'
         do i=1,6
            write(*,'(6i2.2)') (kro(i,j),j=1,6)
         enddo
      endif
C
C *** CONVERTS CARTESIAN 'D','S' and 'MS' TO VOIGT.
      CALL VOIGT (BC_D,AUX_D,AUX66,AUX3333,2)
      CALL VOIGT (BC_S,AUX_S,AUX66,AUX3333,2)
      CALL VOIGT (AUX6v,AUX33,AUX66,AUX_MS,4)


      if(verbose) then
         write(*,*)'bc_d'
         write(*,'(6e11.4)')(bc_d(i),i=1,6)
         write(*,*)'bc_s'
         write(*,'(6e11.4)')(bc_s(i),i=1,6)
         write(*,*)'ibc_s'
         write(*,'(6i2)')(ibc_s(i),i=1,6)
         write(*,*)'ibc_d'
         write(*,'(6i2)')(ibc_d(i),i=1,6)

         write(*,*)'aux_ms in 66 (voigt'
         do i=1,6
            write(*,'(6e11.4)') (aux66(i,j),j=1,6)
         enddo
      endif

      if (is_all_zero(aux66,6)) then
         write(*,*)'aux_ms in 66 (voigt) are all zero'
         do i=1,6
            write(*,'(6e11.4)') (aux66(i,j),j=1,6)
         enddo
         call xit()
      endif

      aux66v(:,:)=0d0
      do i=1,6
         aux6v(i)=-1.d0*ibc_d(i)*bc_d(i)
         do j=1,6
            dum=0d0
            aux6v(i)=aux6v(i)+ aux66(i,j)*ibc_s(j)*bc_s(j)*profac(j)
            dum =ibc_s(j)*kro(i,j) !(i/j)*(j/i)
     $           -ibc_d(j)*aux66(i,j)*profac(j)

            if ((i.eq.1.and.j.eq.2 .or. i.eq.2.and.j.eq.1)
     $           .and.verbose) then
            write(*,*)'ibc_s(j):',ibc_s(j)
            write(*,*)'ibc_s(j):',ibc_s(j)
            write(*,*)'ibc_d(j):',ibc_d(j)
            write(*,*)'ibc_d(j):',ibc_d(j)
            write(*,*)'profac(j)',profac(j)
            write(*,*)'kro(i,j)',kro(i,j)

               write(*,'(a4,x,a1,1i2,a1,1i2,a1,x,e11.4)')'dum:',
     $              '(',i,',',j,')',dum
            endif


            aux66v(i,j)=dum
         enddo
      enddo
c      call xit()

      if (verbose) then
         call w_chrc(0,'aux66v')
         call w_mdim(0,aux66v,6,1d0)
         call w_chrc(0,'aux6v')
         call w_dim(0,aux6v,6,1d0,.true.)
      endif

      CALL LU_EQSYSTEM(AUX66v,AUX6v,6,IER)


      if(ier.eq.1) then
         write(*,*) 'SINGULAR SYSTEM IN STATE_6X6'

         write(*,*)'bc_d'
         write(*,'(6e11.4)')(bc_d(i),i=1,6)
         write(*,*)'bc_s'
         write(*,'(6e11.4)')(bc_s(i),i=1,6)
         write(*,*)'ibc_s'
         write(*,'(6i2)')(ibc_s(i),i=1,6)
         write(*,*)'ibc_d'
         write(*,'(6i2)')(ibc_d(i),i=1,6)

         write(*,*)'aux_ms in 66 (voigt'
         do i=1,6
            write(*,'(6e11.4)') (aux66(i,j),j=1,6)
         enddo
         call xit()
      endif

c      write(*,*)'reached?'
c      call xit()

      do i=1,6
         bc_d(i)=ibc_d(i)*bc_d(i)+ibc_s(i)*aux6v(i)
         bc_s(i)=ibc_s(i)*bc_s(i)+ibc_d(i)*aux6v(i)
      enddo

C *** CONVERTS VOIGT 'D','S' TO CARTESIAN
      CALL VOIGT (BC_D,AUX_D,AUX66,AUX3333,1)
      CALL VOIGT (BC_S,AUX_S,AUX66,AUX3333,1)

      return
      end subroutine state_6x6_cnt



c-----------------------------------------------------------------------
      logical function is_all_zero(array,narray)
      implicit none
      integer narray
      real*8 array(narray,narray)
      integer i,j

      do i=1,narray
         do j=1,narray
            if (abs(array(i,j)).ne.0) then
               is_all_zero=.false.
               return
            endif
         enddo
      enddo

      is_all_zero=.true.
      return
      end function
