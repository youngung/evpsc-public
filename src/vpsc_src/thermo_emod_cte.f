c     --------------------------------------------------------------
c     This subroutine is used to employ the temperature-dependency
c     in elastic modulus and thermal expansion coefficient of single
c     crystal by providing the polynomial fits
c     --------------------------------------------------------------
c
c     Youngung Jeong at MSE, CWNU
c     yjeong@changwon.ac.kr
c
c     iopt=0:
c         - read the order of polynomial fit and c_ij (voigt)
c         - read the order of polynomial fit and CTE_i
c     iopt=1:
c         - estimate c_ij for the given temperature and convert it to b-basis
c         - estimate CTE_i for the given temperature
c     --------------------------------------------------------------
      subroutine thermo_emod_cte(
c     intent(in)
     $     iopt,iunit,iph,nsysmx,nmodmx,nphmx,nmodes,mxndeg,nsm,
     $     temp,schca,
c     intent(inout)
     $     ndegs,cijv_fit,cte_fit,
     $     celcc,selcc,xmu_mode,athcc)
      implicit none

      integer, intent(in):: iopt,iunit,iph,nsysmx,nmodmx,nphmx,
     $     nmodes(nphmx),nsm(nmodmx,nphmx),mxndeg
c     polynomial orders for emod & CTE, respectively.
      integer,intent(inout):: ndegs(nphmx,2)
      real*8,intent(in):: temp,schca(5,nsysmx,nphmx)
      real*8,intent(inout):: cijv_fit(nphmx,0:mxndeg,6,6),
     $     cte_fit(nphmx,0:mxndeg,6),celcc(6,6,nphmx),selcc(6,6,nphmx),
     $     xmu_mode(nmodmx,nphmx),athcc(6,nphmx)
c     locals
c     degree of polynomial function
      integer ideg,i,j
c      parameter(iunit=3014)
      character(len=256) prosa
      real*8 aux66(6,6),cijb(6,6)
c      logical inquiry_exist

      if (iopt.eq.-1) then
         !! preserved case for FE simulation
c        in case of FE simulation, the initialization
c        is preferred to be done inside sdvini for the case
c        of UMAT. For VUMAT, it is not clear where it should be
c        at the moment.
c         call w_chrc(0,'** running in thermo_emod_cte')
c         call w_chrc(0,'** calc_const_elast_cte')
         call calc_const_elast_cte(
     $        nphmx,iph,mxndeg,celcc,athcc,ndegs,cijv_fit,cte_fit)
c         call w_chrc(0,'** calc_elast_cte')
         call calc_elast_cte(
c             intent(in)
     $        temp,iph,nphmx,nsysmx,nmodmx,nsm,nmodes,
     $        ndegs,mxndeg,cijv_fit,cte_fit,schca,
c             intent(out)
     $        cijb,aux66,xmu_mode,athcc)
         celcc(:,:,iph)=cijb(:,:)
         selcc(:,:,iph)=aux66(:,:)

c         call w_chrc(0,'cijb')
c         call w_mdim(0,cijb,6,1d0)
c         call w_chrc(0,'cte')
c         write(*,*)'athcc:',athcc(:,iph)
c         call w_dim(0,athcc(:,iph),6,1d0)
c         call xit()

      elseif (iopt.eq.0) then

c         write(*,*)'iph:',iph
         cijv_fit(iph,:,:,:)=0d0
         cte_fit(iph,:,:)=0d0
c$$$         inquire(file='/tmp/thermal_crystal_property.dat',
c$$$     $        exist=inquiry_exist)
c$$$         if (.not.inquiry_exist) then
c$$$            inquire(file='thermal_crystal_property.dat',
c$$$     $           exist=inquiry_exist)
c$$$            if (.not.inquiry_exist) then !! no temperature-dependent behavior?
c$$$c              since no temperature dependency is given,
c$$$c              cij(t)=cij_constant; thus polynomial function of order zero.
c$$$               call calc_const_elast_cte(
c$$$     $              nphmx,iph,mxndeg,celcc,athcc,ndegs,cijv_fit,cte_fit)
c$$$               call calc_elast_cte(
c$$$c             intent(in)
c$$$     $              temp,iph,nphmx,nsysmx,nmodmx,nsm,nmodes,
c$$$     $              ndegs,mxndeg,cijv_fit,cte_fit,schca,
c$$$c             intent(out)
c$$$     $              cijb,aux66,xmu_mode,athcc)
c$$$               celcc(:,:,iph)=cijb(:,:)
c$$$               selcc(:,:,iph)=aux66(:,:)
c$$$               return
c$$$            else
c$$$               fnmod='thermal_crystal_property.dat'
c$$$               write(*,*)'3',inquiry_exist
c$$$            endif
c$$$         else
c$$$            fnmod='/tmp/thermal_crystal_property.dat'
c$$$         endif
c        READ from <fnmod>
c         open(action='read',unit=iunit,file=fnmod,status='old')

c     ** polynomial fit of elastic modulus and CTE
         read(iunit,'(a)')prosa
cccc     elastic  modulus
c     * Polynomial order of elastic modulus; remember # of coeffs is always order + 1
         read(iunit,'(a)')prosa
c     degree of polynomial used for emod fit
         read(iunit,*)ndegs(iph,1)
         read(iunit,'(a)')prosa
         read(iunit,'(a)')prosa
         do ideg=ndegs(iph,1),0,-1
            read(iunit,*) ((cijv_fit(iph,ideg,i,j),j=1,6),i=1,6)
         enddo
cccc     CTE
c     * Polynomial order of elastic modulus; remember # of coeffs is always order + 1
         read(iunit,'(a)')prosa
c     degree of polynomial used for CTE fit
         read(iunit,*)ndegs(iph,2)
         read(iunit,'(a)')prosa
         read(iunit,'(a)')prosa
         do ideg=ndegs(iph,2),0,-1
            read(iunit,*) (cte_fit(iph,ideg,i),i=1,6)
         enddo
         call calc_elast_cte(
c             intent(in)
     $        temp,iph,nphmx,nsysmx,nmodmx,nsm,nmodes,
     $        ndegs,mxndeg,cijv_fit,cte_fit,schca,
c             intent(out)
     $        cijb,aux66,xmu_mode,athcc)
         celcc(:,:,iph)=cijb(:,:)
         selcc(:,:,iph)=aux66(:,:)

      elseif(iopt.eq.1)then
c         call w_chr(0,'begins calc_elast_cte')
c        updates El modulus and CTE based on the given temperature and fits.
c         call w_val(0,'cijv_fit(1,1,1,1)',cijv_fit(1,1,1,1))
         call calc_elast_cte(
c             intent(in)
     $        temp,iph,nphmx,nsysmx,nmodmx,nsm,nmodes,
     $        ndegs,mxndeg,cijv_fit,cte_fit,schca,
c             intent(out)
     $        cijb,aux66,xmu_mode,athcc)
c         call w_chr(0,'exited calc_elast_cte')
c         call xit()
         celcc(:,:,iph)=cijb(:,:)
         selcc(:,:,iph)=aux66(:,:)
      endif
      end subroutine thermo_emod_cte



c--------------------------------------------------------------------------
c     Given celcc, athcc, assin constant (temperature-insensitive)
c     values (say, y=y_0); meaning that the fit with polynomial func.
c     with 0-th order.
      subroutine calc_const_elast_cte(
c     intent(in)
     $     nphmx,iph,mxndeg,celcc,athcc,
c     intent(out)
     $     ndegs,cijv_fit,cte_fit
     $     )
      implicit none
      integer, intent(in)::nphmx,iph,mxndeg
      real*8, intent(in)::celcc(6,6,nphmx),athcc(6,nphmx)
      real*8, intent(out)::cijv_fit(nphmx,0:mxndeg,6,6),
     $     cte_fit(nphmx,0:mxndeg,6)
      integer, intent(out)::ndegs(nphmx,2)
c     locals
      real*8 aux6(6),aux33(3,3),aux66(6,6),aux3333(3,3,3,3),cijv(6,6)
c     0-th order polynomials
      ndegs(iph,1)=0
      ndegs(iph,2)=0
!! single crystal elastic modulus with b-basis
      aux66(:,:)=celcc(:,:,iph)
!! change that to cartesian
      call chg_basis(aux6,aux33,aux66,aux3333,3,6)
!! change that to voigt
      call voigt(aux6,aux33,cijv,aux3333,4)
      cijv_fit(iph,0,:,:)=cijv(:,:) !! aux66 should be voigt e-emodulus
!! single crystal CTE with b-basis
      aux6(:)=athcc(:,iph)
!! change that to cartesian
      call chg_basis(aux6,aux33,aux66,aux3333,1,6)
!! change that to voigt
      call voigt(aux6,aux33,aux66,aux3333,2)
      cte_fit(iph,0,:)=aux6(:)
      return
      end subroutine calc_const_elast_cte


c--------------------------------------------------------------------------
c     subroutine that calculates below following the polynomial fits
c              1) el modulus
c              2) cte
      subroutine calc_elast_cte(
c         intent(in)
     $     temp,iph,nphmx,nsysmx,nmodmx,nsm,nmodes,
     $     ndegs,mxndeg,cijv_fit,cte_fit,schca,
c         intent(out)
     $     cijb,sijb,xmu_mode,athcc)
      implicit none
      integer, intent(in)::iph,nphmx,mxndeg,nsysmx,nmodmx
     $     ,nsm(nmodmx,nphmx),nmodes(nphmx),ndegs(nphmx,2)
      real*8, intent(in):: temp,cijv_fit(nphmx,0:mxndeg,6,6),
     $     schca(5,nsysmx,nphmx),cte_fit(nphmx,0:mxndeg,6)
      real*8, intent(out)::cijb(6,6),sijb(6,6),xmu_mode(nmodmx,nphmx),
     $     athcc(6,nphmx)
c     locals
      integer ideg,i,j,inverr,imode,isys0,k,l
      real*8 cijv(6,6),aux6(6),aux33(3,3),aux66(6,6),aux3333(3,3,3,3),
     $     aux5(5),aux55(5,5),schmid_33(3,3),cijkl(3,3,3,3),athccv(6)

c     1. elastic modulus for the current temperature
      cijv(:,:)=0d0
      do i=1,6
      do j=1,6
         do ideg=ndegs(iph,1),0,-1 ! (n), (n-1) ... 3, 2, 1, 0
            cijv(i,j)=cijv(i,j)+cijv_fit(iph,ideg,i,j)*(temp**(ideg))
         enddo
      enddo
      enddo

c        voigt to b-basis
      call voigt_bbasis(aux6,aux6,cijv,cijb,3)
      sijb(:,:)=cijb(:,:)
      call lu_inverse(inverr,sijb,6)
      if (inverr.ne.0) then
         write(*,*)'** inside calc_elast_cte within <thermo_emod_cte.f>'
         write(*,*)'temp:',temp
         write(*,*)'ndegs(iph,1):',ndegs(iph,1)
         call w_chrc(0,'cijv')
         call w_mdim(0,cijv,6,1d0)
         call w_chrc(0,'sijb')
         call w_mdim(0,sijb,6,1d0)
         write(*,*)'** Error with lu_inverse in thermo_emod_cte'
         call xit()
      endif

      call chg_basis(aux6,aux33,cijb,cijkl,3,6)

c     2. Retrieve schmid tensor of the first system in each slip/twin mode
      isys0=1
      do imode=1,nmodes(iph)
            ! schmid tensor
         aux5(:)=schca(:,isys0,iph)
         call chg_basis(aux5,schmid_33,aux55,aux3333,1,5)
c           shear modulus of each slip/twin mode
         xmu_mode(imode,iph)=0d0
         do i=1,3
         do j=1,3
         do k=1,3
         do l=1,3
            xmu_mode(imode,iph)=xmu_mode(imode,iph)+
     $           schmid_33(i,j)*cijkl(i,j,k,l)*schmid_33(k,l)
         enddo
         enddo
         enddo
         enddo
!! update isys0
         isys0=isys0+nsm(imode,iph)
      enddo

c     3. determine CTE
      athccv(:)=0d0
      do i=1,6
         do ideg=ndegs(iph,1),0,-1 ! (n), (n-1) ... 3, 2, 1, 0
            athccv(i)=athccv(i)+cte_fit(iph,ideg,i)*(temp**(ideg))
         enddo
      enddo

      call voigt_bbasis(athccv,aux6,aux66,aux66,1)
      athcc(:,iph)=aux6(:)

      return
      end subroutine
