c----------------------------------------------------------------------c
c     Youngung Jeong @ Changwon National University
c     email: yjeong@changwon.ac.kr
c----------------------------------------------------------------------c
c     User material subroutine (UMAT) of EVPSC model is built with
c     sharing some of the subroutines the stand-alone version.
c     Certain subroutines were developed explicitly for interfacing
c     with Abaqus UMAT.
c----------------------------------------------------------------------c
c     Do not try to utilize <umat.f> as it is, which may not be working
c     properly. Build & running of dEVPSC-FE are done via <run.py>
c     python script. Take a look at example input file written for
c     <run.py> script.
c     Also, refer to <docs_src/build_and_run.md> and <finite_element_simulation.md>
c     for more detailed explanations on the usage of the UMAT capacity of
c     EVPSC code for finite element simulations.
c----------------------------------------------------------------------c
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT,
     $     DRPLDE,DRPLDT,STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,
     $     DPRED,CMNAME,NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,
     $     PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
c----------------------------------------------------------------------c
      implicit none
#include <SMAAspUserSubroutines.hdr>
      CHARACTER*80 CMNAME
      dimension STRESS(NTENS),STATEV(NSTATV),DDSDDE(NTENS,NTENS),
     $     DDSDDT(NTENS),DRPLDE(NTENS),STRAN(NTENS),DSTRAN(NTENS),
     $     TIME(2),PREDEF(1),DPRED(1),PROPS(NPROPS),COORDS(3),
     $     DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),JSTEP(4)
c      integer kspt,nstatv,ntens,nshr,ndi,kinc,nprops,layer,
c     $     npt,noel
      real*8 dtemp,dtime,sse,spd,scd,rpl,drpldt,pnewdt,celent,dtime_ref
c----------------------------------------------------------------------c
c     local variables
      INTEGER lenoutdir,iphel,i,j,k,kkk,iopt,i0,igb,igt,istep,
     $     isb,ist,iverbose,nvar,iph,nbeta,ierr_shp,ijv_abq(6,2)
      DATA ((IJV_abq(i,j),j=1,2),i=1,6)/1,1,2,2,3,3,1,2,1,3,2,3/
      character*256 OUTDIR
      real*8 aux6(6),bux6(6),aux33(3,3),aux66(6,6),aux3333(3,3,3,3),
     $     bux3333(3,3,3,3),dsbar_abq_vo(6),dsbarc(3,3),sbarc(3,3),
     $     stress_old(ntens),vnorm,dum1,dum2,dum3,dum4,rot_acum_t(3,3),
     $     bux33(3,3),cux33(3,3),aux5(5),dbarc(3,3),
     $     aux55(5,5),tempdot,rels,reld,xnumstats(7),alplpx_old,
     $     alplpx_diff,pdil
      logical idevwrite,is_all_zero
c     --------------------------------------------------
c     Read inqmx from <umat-par.f>
      integer inqmx
      parameter (inqmx=496)
c      include "umat-par.f"
c     --------------------------------------------------
      integer inquiry(inqmx,2),ninquiry,ielem_type,noel0,noel1,
     $     idiag_nr,idiag,idiag_sc,inverr,ife
      parameter(ife=1) ! 0 is reserved for stand-alone; 1 is for FE.
      logical inquiry_exist,itxwrite
      character*256 path_inquiry,fnout_stv
      character*4 cnoel,cnpt,cstep,cinc
      common /kfe/ inquiry,ninquiry,noel0,noel1,path_inquiry!,ielem_type
c     element number bottom and up.
      integer noel_bot,noel_top,mdim
      common/ktest/noel_bot,noel_top
      real*8 alplstats(7,2)     ! min and max of "alplpx", "RELS", "RELD", "alplpx_diff", "EPSVM", "SVM", "|DSTRAN|"
      common/kstats/alplstats
c     -----------------------------------------------------------------c
c     arrays from common blocks shared with stand-alone evpsc
      include "../../../src/vpsc_src/incr_evpsc_fe.dim"
      integer ioptfs
c     0: Using the linear form of Jeong et al. (IJP, 2021) and Jeon et al. (JMRT, 2024)
c     1: Using the new linear form of Patra et al.
      parameter(ioptfs=1)
c     -----------------------------------------------------------------c
      real*8 ctt(20),fel(ngrmx),fvp(ngrmx)
c     -----------------------------------------------------------------c
c     date and time related
      character(8)::date_
      character(10)::time_
      character(5)::zone_
      integer :: itime_val(8)
c----------------------------------------------------------------------c
      real*8 rhosat(nmodmx),TAU_PROP(nmodmx),TAU_NUCL(nmodmx)
c     quantities for stats analysis from update_crss_dd_fe_new
      real*8 crss_wavg(nmodmx,nphmx), rho_for_ph_mod(nmodmx,nphmx),
     $     rho_deb_ph(nphmx),rho_tot_ph(nphmx),RHOnucl,pnucl,
     $     rho_deb_ph_mod(nmodmx,nphmx),rho_recov,taue_wavg(nphmx,6,4)
c----------------------------------------------------------------------c
c     The variable print_loc is used to switch on/off the printout of
c     locations in UMAT in order to find where the potential error occurs.
      logical iesh_fail,iys,print_loc,isnan
      parameter(print_loc=.false.)
c      parameter(print_loc=.true.)
      integer,parameter:: iprecison = kind(dstran)
c----------------------------------------------------------------------c
      integer aux_nrs(nsysmx),aux_isense(nsysmx),iaux6(6),ibux6(6),
     $     aux_nrs_ph(nsysmx,nphmx)
c----------------------------------------------------------------------c
      real*8 contrib(3),GAVMODG(ngrmx,nmodmx),xid(3,3)
c----------------------------------------------------------------------c
c     Begin to determine FILE I/O UNIT
      include "src/io_head.f"
      include "umat-io.f"
c----------------------------------------------------------------------c
      mdim=6
      if (iregime.eq.1) mdim=5
      if (print_loc) write(*,*)'* 1 - Beginning of dEVPSC-FE UMAT'
      write(cnoel,'(i4.4)') noel
      write(cnpt,'(i4.4)') npt
      write(cstep,'(i4.4)') jstep(1)
      write(cinc,'(i4.4)') kinc
      call cpu_time(ctt(1))
c----------------------------------------------------------------------c
      dum1=0.
      if (ntens.eq.6 .and. ndi.eq.3 .and. nshr.eq.3) then
         ielem_type=0           ! solid
         ishear_comp=4          ! e12, s12
      elseif (ntens.eq.3 .and. ndi.eq.2 .and. nshr.eq.1) then
         ielem_type=1           ! shell
         ishear_comp=3          ! e12, s12
      else
         write(*,*)'Unexpected ntens/ndi/nshr combonation'
         write(*,'(a8,i3)')'ntens:',ntens
         write(*,'(a8,i3)')'ndi:', ndi
         write(*,'(a8,i3)')'nshr:',nshr
         call xit()
      endif
c----------------------------------------------------------------------c
c     End of FILE I/O UNIT declaration
c----------------------------------------------------------------------c
      xid(:,:)=0.
      do i=1,3
         xid(i,i)=1.
      enddo
      jran=-1
c      alplpx=0.
c      alpl(:)=0.
      props(1)=1e-3
      nneigh=0     ! <nneigh> is not considered in delta-EVPSC code,
c                  ! thus is not stored. Sometimes it causes troubles
      dg_trans(:,:)=0d0
      var_vf(:,:,:)=0d0 ! real*8
      sch_stat(:,:,:)=0 ! integer
      var_stat(:,:,:)=0 ! integer

      if (temp.lt.1e-3) temp=298. !! set it as room temp.
      tempdot=dtemp/dtime
      if (jstep(1).eq.1.and.kinc.eq.0) return
      if (iverbose.ge.1.and.jstep(1).eq.1.and.
     $     ((kinc.eq.0.or.kinc.eq.1))) then
         call mutexlock(2)
         if (ielem_type.eq.0) then ! Assumed as SOLID
            noel0=10000         !! common
            noel1=0             !! common
         endif
         call mutexunlock(2)
      endif

c     For restart jobs, <jstep>, <kinc> may not be 1 or 0, and <noel_bot> is fixed to 1
      if (noel_bot.eq.0.and.(jstep(1).gt.1.or.kinc.gt.1)) then
         call mutexlock(1)
         noel_bot=1             ! common
         call mutexunlock(1)
      endif

      call cpu_time(ctt(2))
      call getoutdir(outdir, lenoutdir)

      if (jstep(1).eq.1.and.kinc.eq.1) then
         call mutexlock(3)
         include'src/init_uw_print.f'
         call mutexunlock(3)
      endif

      call cpu_time(ctt(3))

c     Load state variables to EVPSC's arrays
      iopt=2 ! Load EVPSC state variables from STATEV
      include "../../../src/vpsc_src/call_statev_io.f"
      alplpx_old=alplpx
      if (print_loc) write(*,*)'** Restore EVPSC arrays from STATEV'
      if (iregime.eq.3) then
         sbar(:)=sbar_old(:)
      endif
c      gamd0s(:,:)=1d0 !(yj) is this a bug?

      if (iverbose.eq.1) then
         call date_and_time(date_,time_,zone_,itime_val)
         write(*,'(a,2x,2(i2.2,a1),i2.2,a1,i3.3)')
     $        date_,itime_val(5),':',itime_val(6),':',itime_val(7),
     $        '.',itime_val(8)
      endif

c----------------------------------------------------------------------c
c**** Constitutive description starts here
c     Skip calculation for unloaded px is imposed with zero strain increments.
      if ((vnorm(sbar_old,6).eq.0.and.vnorm(dstran,ntens).eq.0))then
         if (print_loc) write(*,*)'** Running elsc_evp for checking'
c     --

c        ** Obtain self-consistent elastic modulus and compliance tensors
         iopt=0
         include "../../../src/vpsc_src/call_elsc_evp_fe.f"
         iopt=1
         include "../../../src/vpsc_src/call_elsc_evp_fe.f"
c        ** Calculates ddsdde by using C^{el}/dt.
c           Notice that the reference frame is not changed, thus is
c           in whatever reference frame assumed for EVPSC.
         if (print_loc) write(*,*)'** Using CELAV for jacobian'
         call convert_evpmod_jac(6,celav,dtime,ddsdde,
     $        ntens,ielem_type)
c--      Save EVPSC state variables to STATEV
         iopt=1
         if (print_loc) write(*,*)'** save dEVPSC arrays to STATEv'
         include "../../../src/vpsc_src/call_statev_io.f"
c--
         stress(:)=0
         return
      endif

      if (dtime.le.1e-30) then
         write(*,*)"** Error: Nearly zero time increment?!"
         write(*,'(a20,e12.4)')"dtime:",dtime
         write(*,*)"          probably something is wrong."
         call xit()
      endif

      if (jstep(1).eq.1.and.kinc.lt.2) then
         call mutexlock(5)
         if (noel.le.noel0) noel0=noel
         if (noel.ge.noel1) noel1=noel
         call mutexunlock(5)
      endif

      if (vnorm(dstran,ntens).gt.1e-2) then
         call mutexlock(5)
         write(*,*)'*Requesting pnewdt change due to too large deps'
         write(*,*)'noel:',noel
         write(*,*)'npt:',npt
         write(*,*)'|dstran:',vnorm(dstran,ntens)
         call mutexunlock(5)
         pnewdt=0.5
         return
      endif


c----------------------------------------------------------------------c
c**** Updates in accumulated rotation matrix
      if (print_loc) write(*,*)'** accumulated rotation matrix'
      aux33(:,:)=rot_acum(:,:)
c     eq 11 of Jeong et al. IJP (2021)
      call mat_dot(drot,aux33,rot_acum,3)
      call mat_transpose_fe(rot_acum,3,rot_acum_t)
c     --
c     obtain self-consistent elastic modulus and compliance tensors
      if (print_loc) write(*,*)'** ELSC modulus and compliance (0)'
      iopt=0
      include "../../../src/vpsc_src/call_elsc_evp_fe.f"
      if (print_loc) write(*,*)'** ELSC modulus and compliance (1)'
      iopt=1
      include "../../../src/vpsc_src/call_elsc_evp_fe.f"
c     end of elsc calculation

c----------------------------------------------------------------------c
c**** Loading condition
      if (ielem_type.eq.0) then ! solid elements (full space)
         if (print_loc) write(*,*)'** SOLID El. loading cond.'
         call update_eacc(ielem_type,dstran,ntens,rot_acum_t,aux6)
         dstran(:)=aux6(:)     ! dstran now in 6D b-basis in initial reference axes.
c--      loading condition (full space)
         if (print_loc) write(*,*)'loading condition fe for solid'
         call loading_condition_fe(
c             intent(in)
c                 dstran components in b-basis
     $        iregime,mdim,ntens,ielem_type,drot,dtime,dfgrd0,dfgrd1,
     $        dstran,
c             intent(out)
     $        strain_control,iDBARv,iSBARv,LIJBARc,DBAR,dsim,WBARc)
c--       end loading condition
      elseif(ielem_type.eq.1) then ! shell (plane-stress)
         if (print_loc) write(*,*)'** SHELL El. loading cond.'
         call update_eacc(ielem_type,dstran,ntens,rot_acum_t,aux6)
         call chg_basis(aux6,aux33,aux66,aux3333,1,6)
         call chg_basis(aux6,aux33,aux66,aux3333,2,6)
         if (print_loc) then
            write(*,*)'loading condition fe for shell'
            write(*,*)'dstran:',dstran
            write(*,*)'dstranb:',aux6
         endif
c--      loading condition (full space)
         call loading_condition_fe(
c             intent(in)
c                 dstran components in b-basis
     $        iregime,mdim,ntens,ielem_type,drot,dtime,dfgrd0,dfgrd1,
     $        aux6,
c             intent(out)
     $        strain_control,iDBARv,iSBARv,LIJBARc,DBAR,dsim,WBARc)
c--      end loading condition
      endif

c ---------------------------------------------------------------
      if (print_loc) then
         write(*,*)'** pre_devpsc'
         write(*,'(a10,6e10.2)')'dsg(:,1)',dsg(:,1)
      endif
      call cpu_time(ctt(4))
c     calc ELSC, update schmid tensor, assign gamd0s
c     and Initial guess
      include "src/prep_evpsc_call.f"
c$$$      if (iregime.eq.2) then
c$$$         dsg(:,:)=sgdot(:,:)*dtime
c$$$         sg(:,:)=sg_old(:,:)+dsg(:,:)
c$$$      endif
c ---------------------------------------------------------------
      if (iverbose.ge.1) then
         call mutexlock(5)
         write(*,*)
         do i=1,128
            write(*,'(a1)',advance='no')'-'
         enddo
         write(*,*)
         write(*,'(3a5,2a3)',advance='no')'st','inc','el','np','sp'
         write(*,'(4a9)',advance='no')'de','edot','evm','svm'
         write(*,'(4a9)',advance='no')'temp','dt','time1','time2'
         write(*,'(2a4)',advance='no')'irg','ifs'
         write(*,'(3a9)',advance='no')'alplpx','appx_old','gamd0s'
         write(*,*)
         write(*,'(3i5,2i3)',advance='no')jstep(1),kinc,noel,npt,kspt
         write(*,'(4e9.2)',advance='no')vnorm(dstran,ntens),
     $        vnorm(dstran,ntens)/dtime,epsvm,svm
         write(*,'(4e9.2)',advance='no')temp,dtime,time(1),time(2)
         write(*,'(2i4)',advance='no')iregime,ioptfs
         write(*,'(3e9.2)',advance='no')alplpx,alplpx_old,gamd0s(1,1)
         write(*,*)
         do i=1,128
            write(*,'(a1)',advance='no')'-'
         enddo
         write(*,*)
         call mutexunlock(5)
      endif
c ---------------------------------------------------------------

      if(print_loc) then
         write(*,*)'* Begins evpsc'
         write(*,'(a10,6e10.2)')'dsg(:,1)',dsg(:,1)
      endif
      istep=-1
      call cpu_time(ctt(5))



      include"../../../src/vpsc_src/call_evpsc.f"

c      if (iregime.eq.1.and.ielem_type.eq.1.and..false.) then
      if (iregime.eq.1.and.ielem_type.eq.1) then
c        enforce the plane-stress condition
         call chg_basis(sbar,sbarc,aux66,aux3333,1,6)
         if (print_loc) then
            write(*,*)'sbar before enforcing plane-stress'
            write(*,'(6e10.2)')sbar
            write(*,*)'sbarc from call_evpsc'
            do i=1,3
               write(*,'(3e10.2)')(sbarc(i,j),j=1,3)
            enddo
         endif
         s33=sbarc(3,3)
         do i=1,3
            sbarc(i,i)=sbarc(i,i)-s33
         enddo
         call chg_basis(sbar,sbarc,aux66,aux3333,2,6)
         if (print_loc) then
            write(*,*)'sbar after enforcing plane-stress'
            write(*,'(6e10.2)')sbar
            write(*,*)'sbarc from call_evpsc'
            do i=1,3
               write(*,'(3e10.2)')(sbarc(i,j),j=1,3)
            enddo
         endif
      endif

      include"include/post_evpsc_diag.f"
      call cpu_time(ctt(6))
      if (print_loc) write(*,*) '* Exiting evpsc'

c     ------------------------------------------------------
c     Update sg_old (sg_old are stored to statev)
      if (print_loc)write(*,*) '* updating sg_old, sbar_old'
      sg_old(:,:)=sg(:,:)
      aux6(:)=dsbar(:)+sbar_old(:) ! in s^(n-l) e^0 basis
      sbar_old(:)=aux6(:)       ! b-basis and e^0 basis
      if (print_loc)write(*,*)'* store sbar as cartesian to Scauchy'
      call chg_basis(aux6,scauchy,aux66,aux3333,1,6) ! e^0 basis
      bux33(:,:)=scauchy(:,:)

c     For solid elements, stress solved in e^(t=0) basis,
c     is rotated back to the corotating frame, i.e., e^(t=tk+dt)
c     For shell elements, I expect rot_acum is just unity (?)
      call tens2rank_transf(bux33,rot_acum,aux33)

      call voigt_abq(aux6,aux33,aux66,aux3333,2)
      if (print_loc) write(*,*) '* updating stress_old'
      if (ielem_type.eq.0) then ! solid
         stress(:)=aux6(:)      ! stress in voigt of abaqus
         stress_old(:)=stress(:)
      elseif (ielem_type.eq.1) then ! shell
c       store updated stress components referred in e^(n+1) basis to aux6
!           stress in voigt convention
         stress(1)=aux6(1)
         stress(2)=aux6(2)
         stress(3)=aux6(4)
c        stress_old in e^(n+1) basis
         stress_old(:)=0d0
         stress_old(1)=stress(1) ! sig_11
         stress_old(2)=stress(2) ! sig_22
         stress_old(3)=stress(3) ! sig_12
      endif

c     ------------------------------------------------------
c     update stress
      if (print_loc) write(*,*) '* update stress'
c     ------------------------------------------------------
c     Update von Mises stress and strain
      if (ielem_type.eq.0) then
         aux6(1:ntens)=stress_old(1:ntens) ! e^(t+dt) basis
      elseif (ielem_type.eq.1.and.mdim.eq.6) then
         aux6(:)=0d0
         aux6(1:2)=stress_old(1:2) ! e^(t+dt) basis
         aux6(4)=stress_old(3)
      endif
      call voigt_abq(aux6,aux33,aux66,aux3333,1) ! e^(t+dt) basis
      call chg_basis(aux6,aux33,aux66,aux3333,2,6)

      if (ielem_type.eq.1.and.mdim.eq.5) then
         aux6(:)=0.
         aux6(1:mdim)=sbar(1:mdim)
      endif

      svm=0d0
      dum2=0d0
      do i=1,5                  !! use only the deviatoric stress, thus upto 5th component
         svm=svm+aux6(i)*aux6(i)
         dum2=dum2+aux6(i)*DBAR(i)*dtime ! delta work
      enddo
      svm=sqrt(svm*3./2.)       ! von mises stress
      if (svm.ne.0) then
         dum2=dum2/svm          ! delta work / seq = delta evm
         epsvm=epsvm+dum2       ! Inrementally integrate von Mises strain
      endif

c     -----------------------------------------------------------------
c     *** calculate internal strain on diffracting planes if idiff.ne.0
      call cpu_time(ctt(7))
      if (print_loc) write(*,*) '* begins diff'
      call DIFF_PLANES_PX(
c          intent(in)::
     $     0,idiff,UR_DIF,-1,1,NPH,ngrmx,nphmx,ndifmx,ngr,ag,wgt,
     $     filediff,ICRYSYMph,iregime,iDBARv,scauchy,sav,BC_EL,mgel,
     $     epstot,wph,
c          intent(inout)::
     $     ndiff,sg,angdetector,nfamily,diff_vs,diff_vc,NGRSETPH,
     $     WGTSETPH,IGRSET,
c          intent(out)::
     $     eelgr)

      call cpu_time(ctt(8))
      if (print_loc) write(*,*) '* ends diff'
c     ------------------------------------------------------
c     Estimate activity of slip/twin modes
      if (print_loc) write(*,*) '* update activity'

      call stat_shear_activity_px(
c          in
     $     NPH,nphmx,ngrmx,ngr,nmodmx,nsysmx,wph,wgt,nmodes,
     $     nsyst,nsm,gamdot,
c          out
     $     activ,gavmod,gavacc,GAVMODG,GAVPH)

c$$$c     --- Report mean and standard deviation of granular stress/strain.
c$$$      call grain_stat(...)
c$$$      if (iverbose.ge.1) then
c$$$         write(*,'(6a11,a10)')'s(ed11)','s(ed22)','s(ed33)',
c$$$     $        's(s11)','s(s22)','s(s33)','deref'
c$$$         write(*,'(6e11.3,e10.3)')(aux6(i),i=1,3),(bux6(i),i=1,3),deref
c$$$      endif
c$$$c     ---

c     -----------------------------------------------------
c     Calculate contribution
      call calc_contribution(
c             intent(in)
     $     -1,1,ngrmx,nphmx,NPH,ngr,dgp,dge,dg_trans,wgt,
c             intent(out)
     $     contrib)

c     ------------------------------------------------------
c     Updating state variables
      call cpu_time(ctt(9))
      include 'src/update_evpsc_statev.f'
c     ------------------------------
c     store EVPSC arrays to STATEV
      if (print_loc) then
         write(*,*)'store EVPSC arrays to STATEV'
         write(*,*)'sbar:',sbar
      endif

      iopt=1
      include "../../../src/vpsc_src/call_statev_io.f"
c     ------------------------------------------------------
c     write inquiry
      include 'src/write_inquiry_umat.f'
c     ------------------------------------------------------
c     Define UVARM
      call cpu_time(ctt(10))
      call mutexlock(6)
      include 'src/uvarm_custom.f'
      call mutexunlock(6)
      call cpu_time(ctt(11))
c     ------------------------------------------------------
c     update jacobian
      if (print_loc) write(*,*) '* update jacobian'
c      if (iregime.eq.1) then
c         write(*,*)'the jacobian should be well thought'
c         call xit()
c      endif

c     I am not sure if we need to get LBAR again here.

      if (mdim.eq.6) then
c     aux66(1:mdim,1:mdim)=mbar(1:mdim,1:mdim)
c         call lu_inverse(inverr,aux66,mdim)
         aux66(1:mdim,1:mdim)=lbar(1:mdim,1:mdim)
         call chg_basis(aux6,aux33,aux66,bux3333,3,mdim)
c         lbar(1:mdim,1:mdim)=aux66(1:mdim,1:mdim)
      elseif(mdim.eq.5) then
c         aux55(1:mdim,1:mdim)=mbar(1:mdim,1:mdim)
c         call lu_inverse(inverr,aux55,mdim)
         aux55(1:mdim,1:mdim)=lbar(1:mdim,1:mdim)
         call chg_basis(aux5,aux33,aux55,bux3333,3,mdim)
c         lbar(1:mdim,1:mdim)=aux55(1:mdim,1:mdim)
      endif

      call tens4rank_transf(bux3333,rot_acum,aux3333)

      if (mdim.eq.6) then
         call chg_basis(aux6,aux33,aux66,aux3333,4,mdim)
         call convert_evpmod_jac(mdim,aux66,dtime,ddsdde,
     $        ntens,ielem_type)
      elseif(mdim.eq.5) then
         call chg_basis(aux5,aux33,aux55,aux3333,4,mdim)
         call convert_evpmod_jac(mdim,aux55,dtime,ddsdde,
     $        ntens,ielem_type)
      endif

c----------------------------------------------------------------------c
      call cpu_time(ctt(12))
      if (print_loc) then
         write(*,*)'** End of UMAT, svm:',Svm,' and Elapsed time (sec):'
         write(*,'(20f6.3)',advance='no') (ctt(i)-ctt(i-1),i=2,12)
         write(*,*)
      endif


c      if (noel.eq.1 .and. npt.eq.1.and.(layer.eq.1)) then
      if (idevwrite.and..false.) then
c        need to write DROT, ROT_ACUM
         write(*,*)'uw_defgrad:',uw_defgrad
         write(*,'(5a7,2(3a11,x))')'step','kinc','noel','npt','layer'
     $        ,'defgrad0','defgrad1'
         write(*,'(5i7)',advance='no') jstep(1),kinc,noel,npt,layer
         write(uw_defgrad,'(4a7,2a11)') 'step','kinc','noel','npt',
     $        't','t+dt'
         write(uw_defgrad,'(4i7,2e11.3)') jstep(1),kinc,noel,npt,
     $        time(2),time(2)+dtime
         write(uw_defgrad,'(4(a33,x),6a11)')
     $        'defgrad0-dij','defgrad1-dij','DROT-dij','ROT_ACUM-dij',
     $        'e11','e22','e12','s11','s22','s12'

         do i=1,3
            write(uw_defgrad,'(4(3e11.3,x),6e11.3)')
     $           (dfgrd0(i,j)-xid(i,j),j=1,3),
     $           (dfgrd1(i,j)-xid(i,j),j=1,3),
     $           (drot(i,j)-xid(i,j),j=1,3),
     $           (rot_acum(i,j)-xid(i,j),j=1,3),
     $           (stran(j)+dstran(j),j=1,2),
     $           (stran(ishear_comp)+dstran(ishear_comp))*0.5,
     $           (stress(j),j=1,2),stress(ishear_comp)
         enddo
      endif

      return
      END subroutine umat











c----------------------------------------------------------------------c
      subroutine update_eacc(ielem_type,dstran,ntens,rot_acum_t,deps_B)
c
c     Given the voigt-notated incremental strain components (dstran),
c     obtain its Cartesian tensor components (in 3x3),
c     then apply the coordinate transformation (from the current co-rotating basis to
c     the initial reference basis), and store the incremental strain components
c     in the B-basis components.


c     (Arguments)
c     ielem_type: 0: solid, 1: shell
c     dstran    : abaqus incremental strain (Voigt) referred in the co-rotating coordinates
c     ntens     : number of indepedent tensor component
c     rot_acum_t: transforms from e0 basis vectors to e(n) basis vectors
c     deps_B: strain increment referenced to the initial reference frame, and stored
c              in B-basis convention
      implicit none

      integer, intent(in)::ielem_type,ntens
      real*8, intent(in)::dstran(ntens),rot_acum_t(3,3)
      real*8, intent(out)::deps_B(6)
c     locals
      integer i,j
      real*8 aux6(6),aux33(3,3),bux33(3,3),aux66(6,6),aux3333(3,3,3,3)
      parameter(iverbose=.false.)
c      parameter(iverbose=.true.)

      if (iverbose) then
         write(*,*)'ielem_type:',ielem_type
         write(*,*)'dstran:'
         write(*,'(6f8.5)')(dstran(i),i=1,ntens)
         write(*,*)'ntens:',ntens
         write(*,*)'rot_acum_t:'
         do i=1,3
            write(*,'(3f5.2)')(rot_acum_t(i,j),j=1,3)
         enddo
      endif

      if (ielem_type.eq.0) then ! solid
c        Convert 6 dimension dstran to Cartesian
         bux33(1,1)=dstran(1)
         bux33(2,2)=dstran(2)
         bux33(3,3)=dstran(3)
         bux33(1,2)=dstran(4)*0.5
         bux33(2,1)=dstran(4)*0.5
         bux33(1,3)=dstran(5)*0.5
         bux33(3,1)=dstran(5)*0.5
         bux33(2,3)=dstran(6)*0.5
         bux33(3,2)=dstran(6)*0.5

      elseif(ielem_type.eq.1) then ! shell
         bux33(:,:)=0.
         bux33(1,1)=dstran(1)
         bux33(2,2)=dstran(2)
         bux33(1,2)=dstran(3)*0.5
         bux33(2,1)=dstran(3)*0.5
      endif

c     Apply the coordinate transformation to the 2nd order strain tensor
c     aux33[i,j] = rot_acum_t[i,k] rot_acum_t[j,l] bux33_[k,l]
      call tens2rank_transf(bux33,rot_acum_t,aux33) ! dstran in e^0 basis vectors
      call chg_basis(aux6,aux33,aux66,aux3333,2,6) ! aux6 in b-basis
      deps_B(:)=aux6(:)

      return
      end subroutine update_eacc


c----------------------------------------------------------------------c
c     subroutines developed only for UMAT interface
      include "src/calc_devm.f"
      include "src/lib_fe.f"        ! some misc. subroutines
      include "src/lib_write_fe.f"  ! output subroutines for diagnosis
      include "src/save_diffraction.f"
      include "src/statev_beta.f"

      include "src/uvarm.f"         ! User defined output accessible from Abaqus CAE (placeholder)
      include "src/voigt_abq.f"     ! to convert stress, strain, compliance (and stiffness) for Abaqus
      include "sdvini.f"
c     external
      include "src/uexternaldb.f"
c     subroutine just borrowed from vpsc
      include "src/calc_spin.f"
      include "src/loading_condition_fe.f"
c     subroutines shared with stand-alone (incr. E)VPSC
      include "../../../src/vpsc_src/calc_contribution.f"
      include "../../../src/vpsc_src/calc_evptilde.f"
      include "../../../src/vpsc_src/calc_rel_err.f"
      include "../../../src/vpsc_src/chg_basis.f"
      include "../../../src/vpsc_src/cod2grains.f"
      include "../../../src/vpsc_src/crossop.f"
      include "../../../src/vpsc_src/crystal_symmetry.f"

      include "../../../src/vpsc_src/data_crystal_fe.f"
      include "../../../src/vpsc_src/data_grain_ebsd.f"
      include "../../../src/vpsc_src/data_grain_fe.f"
      include "../../../src/vpsc_src/diff_planes.f"

      include "../../../src/vpsc_src/elsc_evp_fe.f"
      include "../../../src/vpsc_src/eshelby.f"
      include "../../../src/vpsc_src/eshelby_rot.f"
      include "../../../src/vpsc_src/euler.f"
      include "../../../src/vpsc_src/evpsc_input_fe.f"
      include "../../../src/vpsc_src/evpsc.f"

      include "../../../src/vpsc_src/funcs.f"

      include "../../../src/vpsc_src/grain_rate_and_moduli_evp.f"
      include "../../../src/vpsc_src/grain_stat.f"
      include "../../../src/vpsc_src/grain_stress_evp.f"

      include "../../../src/vpsc_src/initial_state_guess.f"

      include "../../../src/vpsc_src/knuth_shuffle.f"

      include "../../../src/vpsc_src/library7.sub"
      include "../../../src/vpsc_src/lu_libs.f"

      include "../../../src/vpsc_src/neighbours_fe.f"
      include "../../../src/vpsc_src/newton_raphson_evp.f"

      include "../../../src/vpsc_src/read_ebsd.f"
      include "../../../src/vpsc_src/rodrigues.f"

      include "../../../src/vpsc_src/scale_gamd0s_evp_fe.f"
      include "../../../src/vpsc_src/select_points_from_ebsd.f"
      include "../../../src/vpsc_src/solve_evp_bc.f"
      include "../../../src/vpsc_src/statev_io.f"     ! Subroutine to store and restore state variables
      include "../../../src/vpsc_src/stat_shear_activity.f"
      include "../../../src/vpsc_src/stat_twinning_fe.f"
      include "../../../src/vpsc_src/sx_version_tester.f"

      include "../../../src/vpsc_src/thermo_emod_cte.f"
      include "../../../src/vpsc_src/twin_orientation.f"

      include "../../../src/vpsc_src/update_crss_dd_fe_new.f"
      include "../../../src/vpsc_src/update_crss_rgvb_fe.f"
      include "../../../src/vpsc_src/update_crss_voce_fe.f"
      include "../../../src/vpsc_src/update_fij_fe.f"
      include "../../../src/vpsc_src/update_macro_property.f"
      include "../../../src/vpsc_src/update_orientation_fe.f"
      include "../../../src/vpsc_src/update_schmid_fe.f"
      include "../../../src/vpsc_src/update_shape_fe.f"

      include "../../../src/vpsc_src/update_twinning_ptr_fe.f"
      include "../../../src/vpsc_src/update_twinning_vft_fe.f"

      include "../../../src/vpsc_src/voigt.f"
      include "../../../src/vpsc_src/write_texture_fe.f"
