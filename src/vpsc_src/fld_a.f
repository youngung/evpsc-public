c------------------------------------------------------------c
c     subroutine fld_a_bc
c     sets up the boundary condition
c     pertaining to region A of MK model
c------------------------------------------------------------c
      subroutine fld_a_bc(
c          intent(in)
     $     rho,cvm,deij,
c          intent(out)
     $     LIJBARc_lab,dsim_lab,DBAR_lab,scauchy_lab,
     $     sbar_lab,dsbar_lab,dsimctrl,dtime,
     $     strain_control,ictrl,ILBAR,iDBARv,iSBARv
     $     )
      implicit none
      real*8, intent(in)::rho,cvm,deij
      real*8, intent(out):: LIJBARc_lab(3,3),dsim_lab(3,3),
     $     scauchy_lab(3,3),DBAR_lab(6),sbar_lab(6),dsbar_lab(6),
     $     dtime
      integer ,intent(out)::strain_control,
     $     ictrl,ILBAR(3,3),iDBARv(6),iSBARv(6)
c     locals
      real*8 aux66(6,6),aux3333(3,3,3,3),dsimctrl
      integer IJV(6,2),i,j,m,n
      DATA ((IJV(N,M),M=1,2),N=1,6)/1,1,2,2,3,3,2,3,1,3,1,2/

c     sets up boundary condition
      call fld_a_rho2LIJBARc(rho,ictrl,LIJBARc_lab)
      ILBAR(:,:)=0
      ILBAR(1,1)=1
      ILBAR(2,1)=1
      ILBAR(2,2)=1
      ILBAR(3,1)=1
      ILBAR(3,2)=1
      LIJBARc_lab(1,1)=LIJBARc_lab(1,1)*cvm
      LIJBARc_lab(2,2)=LIJBARc_lab(2,2)*cvm
      do i=1,3
      do j=1,3
         dsim_lab(i,j)=(LIJBARc_lab(i,j)+LIJBARc_lab(j,i))/2.
      enddo
      enddo

      call chg_basis(DBAR_lab,dsim_lab,aux66,aux3333,2,6)
      scauchy_lab(:,:)=0.
      call chg_basis(sbar_lab,scauchy_lab,aux66,aux3333,2,6)
      dsbar_lab(:)=sbar_lab(:)
      iDBARv(:2)=1
      iSBARv(:2)=0
      iDBARv(3:)=0
      iSBARv(3:)=1
      strain_control=1
      dsimctrl=dsim_lab(ijv(ictrl,1),ijv(ictrl,2))
      dtime=deij/abs(dsimctrl)

      return
      end subroutine fld_a_bc
c------------------------------------------------------------c
      subroutine fld_a_write(iopt,istep,uw,sbar_lab,dsbar_lab,
     $     LIJBARc_lab,DBAR_lab,epstot_lab,f_lab,dtime)
      implicit none
      integer, intent(in)::iopt,istep,uw
      real*8, intent(in)::sbar_lab(6),dsbar_lab(6),
     $     LIJBARc_lab(3,3),DBAR_lab(6),epstot_lab(3,3),
     $     f_lab(3,3),dtime
c     locals
      real*8 sbar_old_lab(6)
c     opens up file
      if (iopt.eq.0) then
         open(uw,file='region_a.bin',form='unformatted',
     $        access='sequential',status='unknown',action='write')
      elseif(iopt.eq.1) then
         sbar_old_lab(:)=sbar_lab(:)-dsbar_lab(:)
         write(uw) istep,dtime,
     $        LIJBARc_lab,DBAR_lab,epstot_lab,
     $        dsbar_lab,sbar_old_lab,f_lab
      endif
      return
      end subroutine fld_a_write
c------------------------------------------------------------c
c     In the code, range of the given rho defines the strain
c     paths spanning the uniaxial tensile to balanced biaxial.
c      (-0.5 to 1.0) 1-axis being major tensile direction
c      (1.0  to 2.0) 2-axis being major tensile direction
c     Note that (1.0 to 2.0) is transformed to (-0.5 to 1.0)
c     in which 2-axis is being the major tensile direction

c     With transformed rho_t value to match.
c     In the range (1.0 to +),
c      1)  rho = 1   corresponds to rho_t =  1.  (Biaxial)
c      2)  rho = 2.0 corresponds to rho_t =  0.  (plane strain)
c      3)  rho = 2.5 corresponds to rho_t = -0.5 (uniaxial-like)
c      4)  rho = 2.6 corresponds to rho_t = -0.6 (uniaxial-like)

c     Thus rho_t is linearly transformed rho
      subroutine fld_a_rho2LIJBARc(rho0,ictrl,LIJBARc)
      implicit none
      real*8, intent(in):: rho0
      real*8, intent(out):: LIJBARc(3,3)
      integer,intent(out):: ictrl
c     local
      real*8 rho_t
      call fld_a_rho_t(rho0,rho_t)

      if (rho0.le.1.) then
!     conventional LIJBARc is for major axes 1
         LIJBARc(:,:) = 0d0
         LIJBARc(1,1) = 1d0
         LIJBARc(2,2) = rho0
c         LIJBARc(3,3) = -LIJBARc(1,1)-LIJBARc(2,2)
      elseif (rho0.gt.1.) then
c     The above is now valid for (- to 1.0)
c     In the range (1.0 to +) the above is transformed to
         LIJBARc(:,:) = 0.d0
         LIJBARc(1,1) = rho_t
         LIJBARc(2,2) = 1.d0
c         LIJBARc(3,3) = -LIJBARc(1,1)-LIJBARc(2,2)
      endif
      if (abs(LIJBARc(1,1)).ge.abs(LIJBARc(2,2))) then
         ictrl = 1
      elseif (abs(LIJBARc(1,1)).lt.abs(LIJBARc(2,2))) then
         ictrl = 2
      else
         write(*,*) 'Unexpected LIJBARc reported in fld_a_rho_transform'
      endif

      return
      end subroutine fld_a_rho2LIJBARc

c------------------------------------------------------------c
      subroutine fld_a_rho_t(rho0,rho_t)
      implicit none
      real*8,intent(in):: rho0
      real*8,intent(out)::rho_t
      if (rho0.le.1) then       ! Major axis: RD
         rho_t = rho0
      elseif(rho0.gt.1) then    ! Major axis: TD
         rho_t = -1 * (rho0 - 1.) + 1.
      endif
      return
      end subroutine fld_a_rho_t
