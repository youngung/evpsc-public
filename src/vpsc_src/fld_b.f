c-----------------------------------------------------------------------
c     fld_b_bc subroutine
c     (iopt=0)
c       - Initialize the t0 vector based on psi
c     (iopt=1)
c       - Read responses of region A from file and define the boundary
c         condition of region B.
c-----------------------------------------------------------------------
      subroutine fld_b_bc(
c          intent(in)
     $     iopt,ur,psi,ftilde,
     $     sbar_old_lab,
c          intent(inout)
     $     t0,
c          intent(out)
     $     LIJBARc_grv,LIJBARc_a_lab,dsim_grv,DBAR_grv,scauchy_grv,
     $     sbar_grv,dsbar_grv,sbar_a_lab,dtime,
     $     strain_control,ictrl,ILBAR,iDBARv,iSBARv
     $     )
      implicit none
      integer,intent(in)::iopt,ur
      real*8,intent(in)::psi,ftilde,sbar_old_lab(6)
      real*8,intent(inout)::t0(3)
      real*8,intent(out)::LIJBARc_grv(3,3),LIJBARc_a_lab(3,3),
     $     dsim_grv(3,3),DBAR_grv(6),scauchy_grv(3,3),sbar_grv(6),
     $     dsbar_grv(6),sbar_a_lab(6),dtime
      integer,intent(out)::strain_control,ictrl,ILBAR(3,3),iDBARv(6),
     $     iSBARv(6)
c     locals
      integer istep,i,j,k,l
      real*8 pi,DBAR_a_lab(6),
     $     dsbar_a_lab(6),epstot_a_lab(3,3),
     $     sbar_old_a_lab(6),
     $     sbar_old_grv(6),fgrad(3,3),t(3),n(3),
     $     R_labgrv(3,3),R_grvlab(3,3),
     $     dscauchy_grv(3,3),scauchy_a_lab(3,3),aux66(6,6),
     $     aux3333(3,3,3,3),scauchy_old_grv(3,3)
      parameter(pi=3.14159265359)

c     option 0 to determine the initial t vector based on the given psi
      if (iopt.eq.0) then !! initial t0 in lab axes.
         t0(1)=cos(psi+pi/2.)
         t0(2)=sin(psi+pi/2.)
         t0(3)=0.
         return
      elseif(iopt.eq.1) then
c        read region A's state
         read(ur,end=1000) istep,dtime,LIJBARc_a_lab,
     $        DBAR_a_lab,epstot_a_lab,dsbar_a_lab,sbar_old_a_lab,Fgrad
         sbar_a_lab(:)=sbar_old_a_lab(:)+dsbar_a_lab(:)
         call chg_basis(sbar_a_lab,scauchy_a_lab,aux66,aux3333,1,6)

c        update n and t vectors based on deformation gradient tensor of region A.
         do i=1,2
         do j=1,2
            t(i)=t(i)+fgrad(i,j)*t0(j)
         enddo
         enddo
         n(1)=1./sqrt(t(1)**2+t(2)**2) * t(2)
         n(2)=1./sqrt(t(1)**2+t(2)**2) * (-t(1))
c        obtain coordinate transformation matrices
         call fld_b_psi2transf(psi,R_grvlab,R_labgrv)

c        transform sbar_old_lab to sbar_old_grv (region B)
         call mat6rot(sbar_old_lab,R_grvlab,sbar_old_grv)
         call chg_basis(sbar_old_grv,scauchy_old_grv,aux66,aux3333,1,6)

c-----------------------------------------------------------------------
c        Setting up boundary condition for groove axes.
c-----------------------------------------------------------------------
c        1. Setting flags of ILBAR, iDBARv, and iSBARv
         ILBAR(:,:)=1
         ILBAR(1:2,1)=0
         ILBAR(:,3)=0
         iDBARv(:)=0
         iSBARv(:)=1
         iDBARv(2)=1
         iSBARv(2)=0

c        2. Assign known components in vel. grad., stress & strain rate tensors
         LIJBARc_grv(:,:)=0.
         scauchy_grv(:,:)=0.
         dscauchy_grv(:,:)=0.
         do k=1,3
         do l=1,3
c           vel grad.
            LIJBARc_grv(1,2)=LIJBARc_grv(1,2)
     $           +R_grvlab(1,k)*R_grvlab(2,l)*LIJBARc_a_lab(k,l)
            LIJBARc_grv(2,2)=LIJBARc_grv(2,2)
     $           +R_grvlab(2,k)*R_grvlab(2,l)*LIJBARc_a_lab(k,l)
c           stress (for dt)
            scauchy_grv(1,1)=scauchy_grv(1,1)+1./ftilde*
     $           R_grvlab(1,k)*R_grvlab(1,l)*scauchy_a_lab(k,l)
            scauchy_grv(1,2)=scauchy_grv(1,2)+1./ftilde*
     $           R_grvlab(1,k)*R_grvlab(2,l)*scauchy_a_lab(k,l)
         enddo
         enddo

         scauchy_grv(2,1)=scauchy_grv(1,2)
c        incremental cauchy stress tensor
         dscauchy_grv(1,1)=scauchy_grv(1,1)-scauchy_old_grv(1,1)
         dscauchy_grv(1,2)=scauchy_grv(1,2)-scauchy_old_grv(1,2)
         dscauchy_grv(2,1)=dscauchy_grv(1,2)
c        strain rate tensor, i.e., sym part of vel grad.
         do i=1,3
         do j=1,3
            dsim_grv(i,j)=(LIJBARc_grv(i,j)+LIJBARc_grv(j,i))*0.5
         enddo
         enddo

c        3. controling over time (ictrl=7)
         ictrl=7
         strain_control=1

c        4. Change basis from Cartesian to b-basis.
         call chg_basis(DBAR_grv,dsim_grv,aux66,aux3333,2,6)
         call chg_basis(sbar_grv,scauchy_grv,aux66,aux3333,2,6)
         call chg_basis(dsbar_grv,dscauchy_grv,aux66,aux3333,2,6)
      endif                     ! iopt (0, 1)
      return

 1000 write(*,*)'************************************************'
      write(*,*)' EVPSC-FLD for region B failed possibly because '
      write(*,*)' A region did not stretch enough.'
      write(*,*)'************************************************'
      call xit()
      end subroutine fld_b_bc
c-----------------------------------------------------------------------
c     subroutine that calculates the transformation matrix that rotates
c     from the laboratory axes to the groove axes.
c     Basis vectors of laboratory axes are e_grv, and they are aligned
c     parallel with n, t and nxt vectors, where n and t are unit vectors
c     that are perpendicular to and parallel with the groove.
      subroutine fld_b_psi2transf(psi,R_grvlab,R_labgrv)
      implicit none
      real*8,intent(in):: psi                ! radian
      real*8,intent(out):: R_grvlab(3,3),R_labgrv(3,3)
c     locals
      real*8 t(3),n(3),c(3),pi,e_lab(3,3),e_grv(3,3)
      integer i,j,k
      parameter(pi=3.14159265359)
      !! determine t0, n0
      t(1)=cos(psi+pi/2.)
      t(2)=sin(psi+pi/2.)
      t(3)=0.

      n(1)=cos(psi)
      n(2)=sin(psi)
      n(3)=0.

c     components n, t, and c are referenced to the lab axes.
      call crossproduct(n,t,c)
      e_lab(:,:)=0.
      do i=1,3
         e_lab(i,i)=1.
      enddo
      e_grv(1,:)=n(:)
      e_grv(2,:)=t(:)
      e_grv(3,:)=c(:)

      do i=1,3
      do j=1,3
         R_labgrv(i,j)=0.
         do k=1,3
            R_labgrv(i,j)=R_labgrv(i,j)+e_lab(i,k)*e_grv(j,k)
         enddo
      enddo
      enddo
      call transpose(R_labgrv,3,R_grvlab)
      return
      end subroutine
c-----------------------------------------------------------------------
c     subroutine that updates n, t vectors and psi
      subroutine fld_b_updates(fgrad_lab,t0,
     $     epstot_a_lab,epstot_b_lab,ftilde0,
c          intent(out)
     $     ftilde,psi_fld)
      implicit none
      real*8, intent(in)::fgrad_lab(3,3),t0(3),epstot_a_lab(3,3),
     $     epstot_b_lab(3,3),ftilde0
      real*8, intent(out)::psi_fld,ftilde
c     locals
      real*8 t_lab(3),n_lab(3),fact
      integer i,j
c     obtain t_lab
      t_lab(:)=0.
      do i=1,2
      do j=1,2
         t_lab(i)=t_lab(i)+fgrad_lab(i,j)*t0(j)
      enddo
      enddo
c     obtain n_lab
      fact = 1./sqrt(t_lab(1)**2+t_lab(2)**2)
      n_lab(1)=fact*t_lab(2)
      n_lab(2)=-fact*t_lab(1)
      n_lab(3)=0.
c     updates psi
      psi_fld=atan2(n_lab(2),n_lab(1))
c     updates ftilde
      ftilde=ftilde0*exp(epstot_b_lab(3,3)-epstot_a_lab(3,3))

      return
      end subroutine fld_b_updates
c-----------------------------------------------------------------------
      subroutine fld_b_write(
     $     iopt,uw,epstot_a_lab,epstot_b_lab,sbar_a_lab,sbar_b_lab,
     $     LIJBARc_a_lab,LIJBARc_b_lab,psi,ftilde)
      implicit none
      integer, intent(in):: uw, iopt
      real*8, intent(in):: epstot_a_lab(3,3),epstot_b_lab(3,3),
     $     sbar_a_lab(6),sbar_b_lab(6),LIJBARc_a_lab(3,3),
     $     LIJBARc_b_lab(3,3),psi,ftilde
c     locals
      integer k,i,j,ndim,IJV(6,2)
      DATA ((IJV(i,j),j=1,2),i=1,6)/1,1,2,2,3,3,2,3,1,3,1,2/
      real*8 pi,sal(3,3),sbl(3,3),aux66(6),aux3333(3,3,3,3)
      parameter(ndim=6)
      parameter(pi=3.14159265359)
 100     format(6x,a3,2i1)
 105     format(5x,a6)
 106     format(1x,a7)
 110     format(1x,e10.3)
 116     format(1x,f7.3)

      if (iopt.eq.0) then
c        open file
         open(uw,file='FLD.out',status='unknown',action='write')
c        write header
         do k=1,ndim
            write(uw,100,advance='no')'EAL',ijv(k,1),ijv(k,2)
         enddo
         do k=1,ndim
            write(uw,100,advance='no')'EBL',ijv(k,1),ijv(k,2)
         enddo
         do k=1,ndim
            write(uw,100,advance='no')'SAL',ijv(k,1),ijv(k,2)
         enddo
         do k=1,ndim
            write(uw,100,advance='no')'SBL',ijv(k,1),ijv(k,2)
         enddo
         write(uw,105,advance='no')'EDAL33'
         write(uw,105,advance='no')'EDBL33'
         write(uw,105,advance='no')'EB/EA'
         write(uw,105,advance='no')'psi'
         write(uw,106,advance='no')'f'
         write(uw,*)
      elseif(iopt.eq.1)then
c 120    format(a15,6f10.2)
c         write(*,120)'sbar_a_lab',(sbar_a_lab(i),i=1,ndim)
c         write(*,120)'sbar_b_lab',(sbar_b_lab(i),i=1,ndim)
         call chg_basis(sbar_a_lab,sal,aux66,aux3333,1,6)
         call chg_basis(sbar_b_lab,sbl,aux66,aux3333,1,6)
         do k=1,ndim
            write(uw,110,advance='no')epstot_a_lab(ijv(k,1),ijv(k,2))
         enddo
         do k=1,ndim
            write(uw,110,advance='no')epstot_b_lab(ijv(k,1),ijv(k,2))
         enddo
         do k=1,ndim
            write(uw,110,advance='no')sal(ijv(k,1),ijv(k,2))
         enddo
         do k=1,ndim
            write(uw,110,advance='no')sbl(ijv(k,1),ijv(k,2))
         enddo
         write(uw,110,advance='no')LIJBARc_a_lab(3,3)
         write(uw,110,advance='no')LIJBARc_b_lab(3,3)
         write(uw,110,advance='no')LIJBARc_b_lab(3,3)/LIJBARc_a_lab(3,3)
         write(uw,110,advance='no')psi*180./pi
         write(uw,116,advance='no')ftilde
         write(uw,*)
      endif

      return
      end subroutine fld_b_write
c-----------------------------------------------------------------------
