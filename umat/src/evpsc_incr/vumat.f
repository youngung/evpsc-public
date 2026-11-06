c     VUMAT of Delta elasto-visco-plastic self-consistent formulation

c     Youngung Jeong at Changwon National University, Republic of Korea
c     email: yjeong@changwon.ac.kr

c -------------------------------------------------------
c     NOTE: This VUMAT is incomplete and does not work yet.
c     It is included only as a placeholder for future
c     development. Use the UMAT (umat.f) with Abaqus/Standard.
c -------------------------------------------------------

c -------------------------------------------------------
c     Vectorized representation of F, stress, and strain
c     tensors are utilized in the code. use the below
c     table when referring to the specific components of
c     the tensorial quantities in use.
c -------------------------------------------------------
c  deformation gradient (defgradOld, defgradNew)
c -------------------------------------------------------
c component   plane stress   plane strain      3D
c                          (or axisymmetric)
c -------------------------------------------------------
c   1           F11             F11           F11
c   2           F22             F22           F22
c   3           F12             F33           F33
c   4           F21             F12           F12
c   5                           F21           F23
c   6                                         F31
c   7                                         F21
c   8                                         F32
c   9                                         F13
c -------------------------------------------------------
c stress/strain tensors (strainInc, stressOld, StressNew)
c -------------------------------------------------------
c component plane stress  plane strain 3D(VUMAT) 3D(UMAT)
c -------------------------------------------------------
c     1        11            11           11        11
c     2        22            22           22        22
c     3        12            33           33        33
c     4                      12           12        12
c     5                                   23        13
c     6                                   13        23
c -------------------------------------------------------
c  ielem_type  -1            1            0
c -------------------------------------------------------

      subroutine vumat(
C Read only (unmodifiable) variables -
     1  nblock, ndir, nshr, nstatv, nfieldv, nprops, jInfoArray,
     2  stepTime, totalTime, dtArray, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     5  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
C Write only (modifiable) variables -
     7  stressNew, stateNew, enerInternNew, enerInelasNew )
C
c     integer nfieldv,nprops,lanneal,nblock,ndi,nshr,nstatv ! rem
c     default---
      include 'VABA_PARAM.INC'
      parameter (i_info_AnnealFlag = 1,
     *     i_info_Intpt    = 2, ! Integration station number
     *     i_info_layer  = 3,   ! Layer number
     *     i_info_kspt   = 4,   ! Section point number in current layer
     *     i_info_effModDefn = 5, ! =1 if Bulk/ShearMod need to be defined
     *     i_info_ElemNumStartLoc   = 6) ! Start loc of user element number

#include <SMAAspUserSubroutines.hdr>
c     ---------------------------------------------------------------------------
c     Caution! Below variables can be passed in as single or double precision.
C
      dimension props(nprops), density(nblock), coordMp(nblock,*),
     1  charLength(nblock), dtArray(2*(nblock)+1),
     $  strainInc(nblock,ndir+nshr),
     2  relSpinInc(nblock,nshr), tempOld(nblock),
     3  stretchOld(nblock,ndir+nshr),
     4  defgradOld(nblock,ndir+nshr+nshr),
     5  fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6  stateOld(nblock,nstatv), enerInternOld(nblock),
     7  enerInelasOld(nblock), tempNew(nblock),
     8  stretchNew(nblock,ndir+nshr),
     8  defgradNew(nblock,ndir+nshr+nshr),
     9  fieldNew(nblock,nfieldv),
     1  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatv),
     2  enerInternNew(nblock), enerInelasNew(nblock), jInfoArray(*)
c
      character*80 cmname
c     default---
      pointer (ptrjElemNum, jElemNum)
      dimension jElemNum(nblock)

c     -- with include 'VABA_PARAM.INC'
c$$$      real*8 stepTime,dt,totalTime
c$$$      real*8 props,density,coordMp,strainInc,relSpinInc,tempOld,
c$$$     $     stretchOld,defgradOld,fieldOld,stressOld,stateOld,
c$$$     $     enerInternOld,enerInternNew,enerInelasOld,tempNew,stretchNew,
c$$$     $     defgradNew,fieldNew,stressNew,stateNew,enrInterNew,
c$$$     $     enerInelasNew
c$$$      character charLength
c$$$
c$$$      integer nelement(nblock),nMatPoint(nblock),
c$$$     $     nLayer(nblock),nSecPoint(nblock)
c     --

      integer,parameter:: iprecison =kind(dtArray(1))
C
c     ------------------------------------------------------------------
c     evpsc-specific parameters
c      dimension statev(nstatv)
      real*8 statev(nstatv),vnorm,dtime,svmfe,dwork,dstran(ndir+nshr),
     $     stran(ndir+nshr),pdil,deacc_mag
      integer iproc,istep,ntens

c     local auxilary arrays
      real*8 aux33(3,3),aux_av(ndir+nshr+nshr),bux_av(ndir+nshr+nshr),
     $     aux_sv(ndir+nshr),bux_sv(ndir+nshr),bux33(3,3),aux66(6,6),
     $     aux6(6),aux3(3),dum1,aux3333(3,3,3,3),stress_old(ndir+nshr),
     $     dum2,dum3,dum4,u(3,3),rot_acum_t(3,3),drot(3,3),
     $     stress(ndir+nshr),bux3333(3,3,3,3),tempdot
      integer i,i0,km,j,k,iphel,nph,ngr_whole,nvar,isb,ist
      include"umat-par.f"
c     ------------------------------------------------------------------
      integer ife
      parameter(ife=1) ! 0 is reserved for stand-alone; 1 is for FE.
      integer inquiry(inqmx,2),ninquiry,ielem_type,i12_comp,noel0,noel1,
     $     idiag_nr,idiag_sc
      logical inquiry_exist,itxwrite,is_all_zero
      character*256 path_inquiry,fnout_stv
      character*4 cnoel,cnpt,cstep,cinc
      common /kfe/ inquiry,ninquiry,noel0,noel1,path_inquiry,ielem_type
c     -----------------------------------------------------------------
      real*8 ctt(20)
      integer noel,npt,iverbose
      logical verbose,idevwrite
c     other misc variables
      integer lenoutdir,igb,igt,kkk,iopt
      character*256 OUTDIR

      character(8)::date_
      character(10)::time_
      character(5)::zone_
      integer :: itime_val(8)
c----------------------------------------------------------------------c
c     Flags & params that control the behavior of numerical computation.
c     1. <print_loc> is used to switch on/off the printout of locations
c        in UMAT in order to find where the potential error occurs.
c     2. <iextra> is the switch to determine whether or not the
c         extrapolation of the material properties is allowed.
c     3. <idif_calc> is the switch to turn off or on the diffraction
c        calculation. (default is off for vumat)
c     4. <de_crit>: full self-consistent calculation is allowed only when
c          <deacc_mag> exceeds <de_crit>.
c     5. <iodf> is to determine how discrete orientation is represented
c         in FE mesh/integration points.
c     6. <ibc> is to control how the boundary condition is determined.
      logical print_loc,ifullsc,iextra,idif_calc
      parameter(print_loc=.false.)
c      parameter(print_loc=.true.)
      parameter(iextra=.false.) !! use extrapolation to speed up.
c      parameter(iextra=.true.)
      parameter(idif_calc=.false.)
      parameter(ibc=3) !! boundary condition is determined by using F
      integer iodf
      parameter(iodf=0) ! 0: from discrete population, 1: from ODF (mtex); 2: random isotropic texture
      real*8 de_crit
      parameter(de_crit=1e-6)
c----------------------------------------------------------------------c
      include "../../../src/vpsc_src/incr_evpsc_fe.dim"
c----------------------------------------------------------------------c
      integer aux_nrs(nsysmx),aux_isense(nsysmx),
     $     aux_nrs_ph(nsysmx,nphmx)
c----------------------------------------------------------------------c
c     --- sub trials
      real*8 dgam_acc(nsysmx,ngrmx),dsbar_acc(6),dDBAR_acc(6),
     $     dLIJBARc_acc(3,3),dtime_trial,t_remainder
      integer k_trial_sub, k_trial_overall
      logical isuccess,iesh_fail,iys
c----------------------------------------------------------------------c
      include "src/io_head.f"
c----------------------------------------------------------------------c
      real*8 pnewdt
c----------------------------------------------------------------------c
      real*8 rels, reld
c----------------------------------------------------------------------c
      real*8 crss_wavg(nmodmx,nphmx)
c     quantities for stats analysis from update_crss_dd_fe_new
      real*8 rho_for_ph_mod(nmodmx,nphmx),rho_deb_ph(nphmx),
     $     rho_tot_ph(nphmx),rho_deb_ph_mod(nmodmx,nphmx),
     $     rho_recov,taue_wavg(nphmx,6,4),RHOnucl,pnucl
c----------------------------------------------------------------------c

      dtime=real(dtArray(1),8)
      if (dtime.eq.0) then
         write(*,*)'*** Unexpected case: zero incr. time is given'
         call xplb_exit()
      endif

      nneigh=0          ! nneigh is not stored as statev, since it
                        ! sometimes causes troubles.
c                       ! The use of neighbor-coupling is suppressed...
      dg_trans(:,:)=0d0 ! this could be extended to account for
!                         thermal expansion / contraction later.
      lAnneal = jInfoArray(i_info_AnnealFlag)
      iLayer = jInfoArray(i_info_layer)
      kspt   = jInfoArray(i_info_kspt)
      intPt  = jInfoArray(i_info_Intpt)
      iUpdateEffMod = jInfoArray(i_info_effModDefn)
      iElemNumStartLoc = jInfoArray(i_info_ElemNumStartLoc)
      ptrjElemNum = loc(jInfoArray(iElemNumStartLoc))

      iverbose=0 ! quiet
      idevwrite=.false.

      if (print_loc) then
         write(*,*)'*** Beginning'
         write(*,*)'*** mutex initialization'
      endif
c     mutex(1) given to ielem_type
      call mutexinit(1)         ! initialization on ielem_type
      call mutexinit(2)
      call mutexinit(3)
      call mutexinit(4) ! input

      var_vf(:,:,:)=0d0 ! real*8
      sch_stat(:,:,:)=0 ! integer
      var_stat(:,:,:)=0 ! integer

      iregime=2

c     Hereinafter statements are written

      ntens=ndir+nshr


      if (print_loc) then
         write(*,'(a30,4e10.3)')'totaltime,steptime,dtime:',
     $        totaltime,steptime,dtime
         write(*,'(2(a5,i5.3,x))')'ndir:',ndir,'nshr:',nshr
      endif

c     -----------------------------------------
c     ielem_type is a common variable.
      if (totaltime.eq.0) then
         call mutexlock(1)
         ielem_type=-1
         if (ndir.eq.3 .and. nshr.eq.3) then
            ielem_type=0        ! 3D
         elseif (ndir.eq.3 .and. nshr.eq.1) then
            ielem_type=1        !
         else
            write(*,*)'unexpected case'
            call xplb_exit()
         endif
         call mutexunlock(1)
         if (print_loc) then
            write(*,'(a15,i2.2)')'*** ielem_type:',ielem_type
            write(*,'(a15,i2.2)')'*** ndir:',ndir
            write(*,'(a15,i2.2)')'*** nshr:',nshr
         endif
      endif
c     -----------------------------------------
c     initialization of kount_umat_call
c      kount_umat_call = 0
      if (TotalTime.eq.0) then
         write(*,*)'******************************************'
         write(*,*)'* delta-EVPSC VUMAT for Abaqus/Explicit  *'
         write(*,*)'******************************************'
         write(*,*)"*   Youngung Jeong@Changwon Nat'l Univ.  *"
         write(*,*)'*                 yjeong@changwon.ac.kr  *'
         write(*,*)'******************************************'
         write(*,*)'** VUMAT initial calculations and check **'
         write(*,*)'******************************************'
         write(*,'(a26,x,a10)')"*** Mat'l name in input:",cmname(1:20)
         if (print_loc) write(*,*)'*** begins <call_statev_io> for test'
         iopt=0                 ! Test if nstatv is properly sized.
         include "../../../src/vpsc_src/call_statev_io.f"
         if (print_loc) write(*,*)'*** ends <call_statev_io>'
      endif

      if (iverbose.eq.1) then
         call date_and_time(date_,time_,zone_,itime_val)
         write(*,'(a,2x,2(i2.2,a1),i2.2,a1,i3.3)')
     $        date_,itime_val(5),':',itime_val(6),':',itime_val(7),
     $        '.',itime_val(8)
      endif
c     get outdirectory
      call vgetoutdir(outdir, lenoutdir)
      call mutexlock(2)
      CALL ESHELBY(AUX3,AUX3333,0d0,AUX3333,AUX3333,AUX33,AUX33,
     $     pdil,AUX3333,AUX3333,0,iesh_fail)
      call mutexunlock(2)

      do 100 km=1,nblock
         tempdot=(tempnew(km)-tempold(km))/dtime
         write(cnoel,'(i4.4)') intPt
         write(cnpt,'(i4.4)') km

         iverbose=0
         if (km.eq.1) iverbose=1

         dstran(:)=strainInc(km,:)  !strain increment
         statev(:)=stateOld(km,:)   !state variables
c     sv:  symmetric voigt (stress, strain, stretch)
c     av: asymmetric voigt (deformation gradient, spin, rotation)

         aux_sv(:)=stretchNew(km,:) !New stretch tensor (U)
         bux_sv(:)=stretchOld(km,:) !Old stretch tensor (U)
         aux_av(:)=defgradnew(km,:) !New Deformation Grad
         bux_av(:)=defgradold(km,:) !Old Deformation Grad


         if (print_loc) then
            write(*,*)'*** strainInc'
            do i=1,ndir+nshr
               write(*,'(f9.4)',advance='no')strainInc(km,i)
            enddo
            write(*,*)
            write(*,*)'*** old and new stretch tensors'
            do i=1,ndir+nshr
               write(*,'(f9.4)',advance='no')stretchOld(km,i)
            enddo
            write(*,'(x,a1,x)',advance='no')'|'
            do i=1,ndir+nshr
               write(*,'(f9.4)',advance='no')stretchNew(km,i)
            enddo
            write(*,*)
            write(*,*)'*** old and new def. grad. tensors'
            do i=1,ndir+nshr+nshr
               write(*,'(f9.4)',advance='no')defgradOld(km,i)
            enddo
            write(*,'(x,a1,x)',advance='no')'|'
            do i=1,ndir+nshr+nshr
               write(*,'(f9.4)',advance='no')defgradNew(km,i)
            enddo
         endif



         is_initialized_statev=0

         if (print_loc) then
            write(*,*)'*** begins <call_statev_io.f>, read from statev'
         endif
         iopt=2                 ! Load EVPSC state variables from STATEV
         statev(:)=stateold(km,:)
         include "../../../src/vpsc_src/call_statev_io.f"
         if (print_loc) write(*,*)'*** ends <call_statev_io.f>'


         ischme=1
c     ------------------------------
         include "vumat-io.f"   ! (file IO) is dependent on <km>


         if (totaltime.eq.0.or.is_initialized_statev.ne.1) then
            if (iverbose.ge.1) then
               write(*,*)
               write(*,*)'----------------------------------'
               write(*,*)'*** Reading polycrystal states ***'
               write(*,*)'----------------------------------'
               write(*,*)
            endif
            call cpu_time(t0)
            ipxx=0 !! hardwired to supress reading statev_px in the case of vumat
            include"src/load_statev.f"
            call cpu_time(t10)
            include 'src/init_uw_print.f' ! done after 'vumat-io.f' is performed

c           store block-specific deformation gradient & stretch tensors
c           for the first time - at the first increment when no state variables
c           are initialized.
            do i=1,ndir+nshr+nshr
               f_sc(i)=defgradOld(km,i)
            enddo
            do i=1,ndir+nshr
               u_sc(i)=stretchOld(km,i)
            enddo
         endif                  ! of if (totaltime.eq.0.or.is_initialized_statev.ne.1) then
c        end of initial increment in which state variables are initialized.





c----------------------------------------------------------------------c
c     Constitutive description starts here
         gamd0s(:,:)=1d0

c----------------------------------------------------------------------c
c**** Updates in accumulated rotation matrix
         if (print_loc) then
            write(*,*)'*** accumulate rotation matrix'
            write(*,*)'*** rot_acum so far'
            write(*,'(3e11.3)')((rot_acum(i,j),j=1,3),i=1,3)
         endif
         aux33(:,:)=rot_acum(:,:)
c     obtain drot = R^new . (R^old)^T
         call calc_drot_explicit(
c     intent(in)
     $        ielem_type,ndir,nshr,
     $        aux_sv, ! stretchnew
     $        bux_sv, ! stretchold
     $        bux_av, ! fold
     $        aux_av, ! fnew
c     intent(out)
     $        drot)
         if (print_loc) then
            write(*,*)'*** drot'
            write(*,'(3e11.3)')((drot(i,j),j=1,3),i=1,3)
         endif
         if (print_loc) write(*,*)'*** end of <calc_drot_explicit>'
         call mat_dot(drot,aux33,rot_acum,3)
         call mat_transpose_fe(rot_acum,3,rot_acum_t)
c     accumulated dtimes
         dtacc=dtacc + dtime
c     accumulated dstrains in B-basis referred to e^0
c        Calculate deacc_b, in the global coordinate system
         call update_eacc_explicit(ielem_type,dstran,ntens,rot_acum_t,
     $        aux6)
c      write(*,'(a8,6e11.2)')'aux6:',(aux6(i),i=1,6)
         deacc(:)=deacc(:)+aux6(:) ! dstran in e0 and B-basis
c     accumulated drot_acc
         aux33(:,:)=drot_acc(:,:)
         drot_acc(:,:)=0d0
         do i=1,3
            do j=1,3
               do k=1,3
                  drot_acc(i,j)=drot_acc(i,j)+drot(i,k)*aux33(k,j)
               enddo
            enddo
         enddo

         if (print_loc) write(*,*)'*** obtain deacc'
c     obtain deacc
         call calc_devm(deacc,ntens,deacc_mag) ! dum1 stores evm_acc (in b-basis)
         if (print_loc) write(*,*)'deacc_mag:',deacc_mag
c     obtain theta_acc
         dtheta_acc=0d0
         do i=1,3
            dtheta_acc=dtheta_acc+drot_acc(i,i)
         enddo
         dtheta_acc=acos((dtheta_acc-1d0)/2d0)*180d0/3.14159265359d0
c     obtain dth_tot
         dth_tot=0d0
         do i=1,3
            dth_tot=dth_tot+rot_acum_t(i,i)
         enddo
         dth_tot=acos((dth_tot-1d0)/2d0)*180d0/3.14159265359d0

         if (print_loc) write(*,*)'*** determine ifullsc'
         ifullsc=.true.
c###     determine if linear extrapolation or full SC calculation to perform
c      if (kinc.eq.4.or.kinc.eq.5) then
c         ifullsc=.false.
c         iextkount=iextkount+1
c      endif
         if (iextra.and.abs(deacc_mag).lt.de_crit
     $        .and.kount_umat_call.gt.2) then
            ifullsc=.false.
            iextkount=iextkount+1
         endif

         if (print_loc) write(*,*)'ifullsc:',ifullsc
c     --
         if (ifullsc) then
            dtime=dtacc*1d0
            dstran(:)=deacc(:)  ! in b-basis
            drot(:,:)=drot_acc(:,:)
         else                   !! extrapolation
c            dtime=dtime
            drot(:,:)=0d0
            do i=1,3
               drot(i,i)=1d0
            enddo
            call update_eacc_explicit(ielem_type,dstran,ntens,rot_acum_t,
     $           aux6)
            dstran(:)=aux6(:)   ! dstran in e0 and B-basis
         endif

c     --
         if (ifullsc) then
c     obtain self-consistent elastic modulus and compliance
            if (print_loc) write(*,*)'*** calculate celav'
            iopt=0
            include "../../../src/vpsc_src/call_elsc_evp_fe.f"
            if (print_loc) write(*,*)'end of elsc_evp_fe(0)'
            iopt=1
            include "../../../src/vpsc_src/call_elsc_evp_fe.f"
            if (print_loc) write(*,*)'end of elsc_evp_fe(1)'
c     end of elsc calculation
            dtime=dtacc*1d0
            dstran(:)=deacc(:)  ! in b-basis
            drot(:,:)=drot_acc(:,:)
c        store DBAR0_sc, mbar_dt_sc to calculate rates of changes in
c     DBAR0 and (mbar)/(Delta t) that are used for extrapolation steps.
         else                   !! extrapolation
            dtime=dtime*1d0     ! to fool compiler
            drot(:,:)=0d0
            do i=1,3
               drot(i,i)=1d0
            enddo
            call update_eacc_explicit(ielem_type,dstran,ntens,rot_acum_t
     $           ,aux6)
            dstran(:)=aux6(:)   ! dstran in e0 and B-basis
         endif
         aux_sv(:)=stretchNew(km,:) !New stretch tensor (U)
         bux_sv(:)=stretchOld(km,:) !New stretch tensor (U)
         aux_av(:)=defgradNew(km,:)
         bux_av(:)=defgradOld(km,:)
         if (print_loc) then
            write(*,*)'*** loading condition_fe_explicit'
            write(*,*)'strainInc'
            write(*,'(6e11.3)')(strainInc(km,i),i=1,ntens)
            write(*,*)'dstran'
            write(*,'(6e11.3)')(dstran(i),i=1,ntens)
            write(*,*)'dstran/dtime'
            write(*,'(6e11.3)')(dstran(i)/dtime,i=1,ntens)
         endif
c$$$         if (km.eq.1) then
c$$$            write(*,*)'strainInc'
c$$$            write(*,'(6e11.3)')(strainInc(km,i),i=1,ntens)
c$$$         endif

c**** loading condition
c     There are multiple ways for determining boundary conditions.
c     ibc=0: using strainincr/dt as the LIJBARc
c     ibc=1:
c     ibc=2: use deformation gradient
c     ibc=3: use deformation gradient (decomposed by in-house code)

         if (ifullsc) then
c           use Fnew, F_sc, Unew, U_sc, dtacc
c           overwrite stretchOld and defgradOld
c$$$            bux_av(:)=F_sc(:)
c$$$            bux_sv(:)=U_sc(:)
            do i=1,ndir
               bux_av(i)=f_sc(i)
               bux_sv(i)=u_sc(i)
            enddo
            do i=ndir+1,ndir+nshr
               bux_sv(i)=u_sc(i)
            enddo
            do i=ndir+1,ndir+nshr+nshr
               bux_av(i)=f_sc(i)
            enddo

c            write(*,*)'calling loading_condition_fe_explicit'
c            write(*,*)'ielem_type:',ielem_type
            call loading_condition_fe_explicit(
c     intent(in)
c     when ioption.eq.2, dstran is not used.
     $           ibc,ielem_type,dstran,aux_av,bux_av,aux_sv,bux_sv,
     $           dtime,ndir,nshr,
c     intent(out)
     $           DBAR,dsim,scauchy,WBARc,LIJBARc,ILBAR,iDBARv,iSBARv,
     $           rot_acum)
c            write(*,*)'exiting loading_condition_fe_explicit'
         else ! extrapolation
c           use Fnew, Fold, Unew, Uold
            call loading_condition_fe_explicit(
c     intent(in)
c     when ioption.eq.2, dstran is not used.
     $           ibc,ielem_type,dstran,aux_av,bux_av,aux_sv,bux_sv,
     $           dtime,ndir,nshr,
c     intent(out)
     $           DBAR,dsim,scauchy,WBARc,LIJBARc,ILBAR,iDBARv,iSBARv,
     $           rot_acum)
         endif



c**** end loading condition
         if (print_loc) then
            write(*,*)'** DBAR'
            write(*,'(6e11.3)')(DBAR(i),i=1,ntens)
         endif

c         write(*,*)'iverbose at 489:,',iverbose

         if (iverbose.ge.1
     $        ) then
            write(*,*)
            write(*,*)'------------------------------------------'
            if(ifullsc) then
               write(*,*)'*** Full SC calculation  ***',
     $              'intPt:',intPt
c     $              'iElemNumStartLoc:',iElemNumStartLoc
            else
               write(*,*)'*** Linear extrapolation ***',
     $              'intPt:',intPt
c     $              'iElemNumStartLoc:',iElemNumStartLoc
            endif
            write(*,'(2a12,2a6,5a11)')'tTime','sTime','noel','npt',
     $           'dt','dt_acc','|deacc|','dth_acc','dth_tot'
            write(*,'(2e12.4,2i6,3e11.3,2f11.7)')TotalTime,StepTime,
     $           noel,npt,dtime,dtacc,deacc_mag,dtheta_acc,dth_tot
            write(*,*)'------------------------------------------'
         endif

c ---------------------------------------------------------------
c     Streaming out the information of current incremental step.
         if (print_loc) write(*,*)'*** streaming out the information'
         if (iverbose.ge.1)then
            write(*,'(a5,e10.3,x,a4,e10.3,x,a6,e9.2,x,a7,e11.4,x,
     $a3,e10.3,x,a3,e10.3,x,a9,e10.3)')
     $           'dtime',dtime,'|de|',vnorm(dstran,ntens),
     $           '|edot|',vnorm(dstran,ntens)/dtime,'Ttime: ',TotalTime,
     $           'Evm',epsvm,'Svm',svm,'acctimeFE',acctimeFE
            write(*,*)'drot and WBARc'
            do i=1,3
c            write(*,'(2(3f11.3,5x))')
               write(*,'(2(3e11.3,5x))')
     $              (drot(i,j),j=1,3),
     $              (WBARc(i,j),j=1,3)
            enddo

c            r32=relSpinInc(km,1)
c            r13=relSpinInc(km,2)
c            r21=relSpinInc(km,1)
c$$$            write(*,'(a20,3e11.3)')'r32,r13,r21:',
c$$$     $           r32,r13,r21
c$$$            write(*,'(a20,3e11.3)')'r32,r13,r21/dt',
c$$$     $           r32/dtime,r13/dtime,r21/dtime

         endif

         if (ifullsc) then
            iextkount=0         !! initialize the number of 'continued' extrapolation steps
c ---------------------------------------------------------------
c     calc ELSC, update schmid tensor, assign gamd0s
            if (print_loc) write(*,*)'*** pre_devpsc begins'
            include "src/prep_evpsc_call.f"
            if (print_loc) write(*,*)'*** pre_devpsc ends'
c     ------------------------------------------------------
c     Solve the elasto-visco-plastic equation

c     (Maybe) better initial guess on mbar (20210427)
            if (( TotalTime.eq.0.and.StepTime.gt.0).or.
     $           (TotalTime.gt.0).and..false.) then
               if (dtime.ne.0) then
                  mbar(:,:)=mbar_dt_sc(:,:)/dtime
               endif

               if (iverbose.ge.1.and..false.)then
                  write(*,*)'Initial guess on DBAR0'
                  write(*,'(6e10.2)')(DBAR0(i),i=1,6)
                  write(*,*)'Initial guess on mbar'
                  do i=1,1      !6
                     write(*,'(6f10.6)')(mbar(i,j),j=1,6)
                  enddo
               endif
            endif

c ---------------------------------------------------------------
c     Initial guess:   (out): DBAR0, DBAR, sbar_old, dsg, sg, sg_old, dsgtry,sbar,mbar,dsbar ...)
            if (TotalTime.eq.0.and.StepTime.eq.0) then
               call initial_state_guess_el_evp(
c                  in
     $              iregime,iDBARv,iSBARv,dtime,time_incr_deriv,
     $              NPH,nphmx,ngrmx,
     $              ngr(NPH)-ngr(0),ngr,bc_el,selav,
c                  inout
     $              DBAR0,DBAR,sbar_old,
c                  out
     $              dsg,sg,sg_old,dsgtry,sgtry,sbar,mbar,dsbar)
            else
               istep=-1
               if (print_loc) write(*,*) '*** Begins evpsc'
               call cpu_time(ctt(5))
               idiag=0
               iproc=0
               include"../../../src/vpsc_src/call_evpsc.f"
c              Update sg_old (sg_old are stored to statev)
               DO IPHEL=1,NPH !IPHELBOT,IPHELTOP
                  DO KKK=NGR(IPHEL-1)+1,NGR(IPHEL)
                     sg_old(:,kkk)=sg(:,kkk)
                  enddo
               enddo
               call cpu_time(ctt(6))
            endif

            if (iverbose.ge.1.and.print_loc)then
               write(*,*)'*** full self-consistent calculation'
               write(*,*)'- DBAR0 for the current increment'
               write(*,'(6e10.2)')(DBAR0(i),i=1,6)
               write(*,*)'- mbar for the current increment'
               do i=1,1         !6
                  write(*,'(6f10.6)')(mbar(i,j),j=1,6)
               enddo
            endif

            if (print_loc) write(*,*) '*** Exiting evpsc'
c     ------------------------------------------------------
c     assign new sbar to scauchy and deviat
            if (iverbose.ge.1.and.print_loc) then
               write(*,*)'- DBAR'
               write(*,'(6e11.3)') (DBAR(i),i=1,6)
               write(*,*)'- dsbar'
               write(*,'(6e11.3)') (dsbar(i),i=1,6)
               write(*,*)'- sbar_old_sc'
               write(*,'(6e11.3)') (sbar_old_sc(i),i=1,6)
            endif
            aux6(:)=dsbar(:)+sbar_old_sc(:) ! in s^(n-l) e^0 basis
            sbar_old(:)=aux6(:) ! b-basis and e^0 basis
            call chg_basis(aux6,scauchy,aux66,aux3333,1,6) ! e^0 basis

            if (iverbose.ge.1.and.print_loc)then
               write(*,*)'- sbar (b-basis) after current increment'
               write(*,'(6e11.3)') (sbar_old(i),i=1,6)
               write(*,*)'- scauchy after current increment'
               write(*,'(3e11.3)') ((scauchy(i,j),j=1,3),i=1,3)
            endif

            !call chg_basis(aux6,bux33,aux66,aux3333,1,6)
            call tens2rank_transf(scauchy,rot_acum,aux33)
            call voigt_exp(aux_av,aux33,aux_sv,ndir,nshr,ielem_type,4)

c        !! sbar at the last self-consistent calculation
c        !! mbar_sc/dt  at the last self-consistent calculation
            sbar_old_sc(:)=sbar_old(:)

            stressnew(km,:)=aux_sv(:)
            stress_old(:)  =aux_sv(:)

            if (iverbose.ge.1.and.print_loc)then
               write(*,*)'- stress new:'
               write(*,'(6e11.3)')(stressnew(km,i),i=1,ndir+nshr)
            endif

            if (TotalTime.eq.0.and.StepTime.eq.0) then
               if (print_loc) then
                  write(*,*)'*** skipping the rest for the first incr.'
               endif
               goto 100
            endif
!        Case of property extrapolation
         else
            if (iverbose.ge.1.and.print_loc) then
               write(*,*)
               write(*,*)'*** solving extrapolation dEVPSC equation'
               write(*,'(a9,e11.3)')'dtacc:',dtacc
               write(*,*)'- mbar_dt_dot'
               do i=1,1         !6
                  write(*,'(6e11.3)')(mbar_dt_dot(i,j),j=1,6)
               enddo
               write(*,*)'- mbar_dt_dot*dtacc'
               do i=1,1         !6
                  write(*,'(6e11.3)')(mbar_dt_dot(i,j)*dtacc,j=1,6)
               enddo
            endif

            DBAR0(:)=DBAR0_sc(:)+DBAR0_dot(:)*dtacc
            mbar(:,:)=(mbar_dt_sc(:,:)+mbar_dt_dot(:,:)*dtacc)
     $           /dtime

            if (iverbose.ge.1.and.print_loc) then
               write(*,*)'- mbar*dtime now'
               do i=1,1         !6
                  write(*,'(6e11.3)')(mbar(i,j)*dtime,j=1,6)
               enddo
               write(*,*)'- mbar now for extrapolation'
               do i=1,1         !6
                  write(*,'(6e11.3)')(mbar(i,j),j=1,6)
               enddo
               write(*,*)'- DBAR0 now for extrapolation'
               write(*,'(6e11.3)')(DBAR0(i),i=1,6)
            endif

            call solve_bc_evp(
c     in
     $           iDBARv,iSBARv,DBAR0,mbar,
c     inout
     $           DBAR,dsbar)

            aux6(:)=dsbar(:)+sbar_old(:) ! in e^0 basis
            sbar_old(:)=aux6(:)
            call chg_basis(aux6,scauchy,aux66,aux3333,1,6) ! e^(0) basis.

            if (iverbose.ge.1.and.print_loc)then
               write(*,*)'- sbar (b-basis) after current increment'
               write(*,'(6e11.3)') (sbar_old(i),i=1,6)
               write(*,*)'- scauchy after current increment'
               write(*,'(3e11.3)') ((scauchy(i,j),j=1,3),i=1,3)
            endif

            call tens2rank_transf(scauchy,rot_acum,aux33) ! in e^(n+1) basis
            call voigt_exp(aux_av,aux33,aux_sv,ndir,nshr,ielem_type,4)

            stressNew(km,:)=aux_sv(:)
            stress_old(:)  =aux_sv(:)

            if (iverbose.ge.1.and.print_loc)then
               write(*,*)'- stress new:'
               write(*,'(6e11.3)')(stressnew(km,i),i=1,ndir+nshr)
            endif

         endif ! end of self-consistent (or extrapolation) calculation

c     ------------------------------------------------------
c     update stress
         if (print_loc) write(*,*) '*** update VM stress / VM strain'
c     ------------------------------------------------------
c     Update von Mises stress and strain
         aux6(:)=0d0
         if (ielem_type.eq.0) then
            do i=1,ntens
               aux6(i)=stress_old(i) ! e^(t+dt) basis
            enddo
         elseif (ielem_type.eq.1) then
            aux6(:)=0d0
            aux6(1:2)=stress_old(1:2) ! e^(t+dt) basis
            aux6(4)=stress_old(4)
         endif
         call voigt_abq(aux6,aux33,aux66,aux3333,1) ! e^(t+dt) basis
         call chg_basis(aux6,aux33,aux66,aux3333,2,6)
         svm=0d0
         dum2=0d0
         do i=1,5
            svm=svm+aux6(i)*aux6(i)
            dum2=dum2+aux6(i)*DBAR(i)*dtime ! delta work
         enddo
         svm=sqrt(svm*3./2.)
         if (svm.ne.0) then
            dum2=dum2/svm       ! delta work / seq = delta evm
c     von mises strain
            epsvm=epsvm+dum2
         endif

c     -----------------------------------------------------------------
c     *** calculate internal strain on diffracting planes if idiff.ne.0

         if (ifullsc) then
            if (idif_calc) then
               if (print_loc) write(*,*) '*** begins diff'
               call DIFF_PLANES_PX(
c     intent(in)::
     $              1,idiff,UR_DIF,-1,1,NPH,ngrmx,nphmx,ndifmx,ngr,
     $              ag,wgt,filediff,ICRYSYMph,iregime,iDBARv,scauchy,sav,
     $              BC_EL,mgel,epstot,wph,
c     intent(inout)::
     $              ndiff,sg,angdetector,nfamily,diff_vs,diff_vc,
     $              NGRSETPH,WGTSETPH,IGRSET,
c     intent(out)::
     $              eelgr)

               if (print_loc) write(*,*) '*** ends diff'
            endif
c     ------------------------------------------------------
c     Estimate activity of slip/twin modes
            if (print_loc) write(*,*) '*** update activity'

            call stat_shear_activity_px(
c             in
     $           NPH,nphmx,ngrmx,ngr,nmodmx,nsysmx,wph,wgt,nmodes,
     $           nsyst,nsm,gamdot,
c             out
     $           activ,gavmod,gavacc,GAVMODG)

c     ------------------------------------------------------
c     Writing activities

            if (idevwrite) then
               if (print_loc) write(*,*) '*** writing uwdeb'
               write(uwdeb,'(6(i3,x),9e11.3)')noel,npt,1,NPH,
     $              1,ngr(NPH),acctimeFE,totaltime,steptime,
     $              dtime,
     $              vnorm(stran,ntens),vnorm(dstran,ntens),
     $              vnorm(stress,ntens)
     $              ,vnorm(DBAR,6),vnorm(dsbar,6)
               if (print_loc) write(*,*) '*** writing texture'
               do iphel=1,NPH
                  igb=ngr(iphel-1)+1
                  igt=ngr(iphel)
                  call write_texture_fe(uwfe01,iphel,ngrmx,nphmx,igb,igt
     $                 ,epsvm,ag,wgt,axisph,eulerph)
               enddo
            endif

c     ------------------------------------------------------
c     Updating state variables
            if (print_loc) write(*,*) '*** begins update_evpsc_state.f'
            include 'src/update_evpsc_statev.f'
            if (print_loc) write(*,*) '*** ends update_evpsc_state.f'

c     ------------------------------
c     store EVPSC arrays to STATEV
            if (print_loc) write(*,*) '*** begins <call_statev_io.f>'
            is_initialized_statev=1
            iopt=1
            include "../../../src/vpsc_src/call_statev_io.f"
            stateNew(km,:)=statev(:)
            if (print_loc) write(*,*) '*** ends <call_statev_io.f>'

c     ------------------------------------------------------
c     write inquiry
            if (print_loc) write(*,*) '*** begins <write_inquiry>'
            include 'src/write_inquiry_vumat.f'
            if (print_loc) write(*,*) '*** ends <write_inquiry>'
         endif                  !! end of full self-consistent calculation

c     ------------------------------------------------------
c     Define UVARM
c$$$         if (print_loc) write(*,*) '* begins uvarm_custom'
c$$$         call mutexlock(5)
c$$$         include 'src/uvarm_custom.f'
c$$$         call mutexunlock(5)
c$$$         if (print_loc) write(*,*) '* end uvarm_custom'

c$$$      if (iverbose.ge.1)then
c$$$         write(*,*)'end of full self-consistent calculations?'
c$$$      endif

c     ------------------------------
c        Obtain new d(EVP compliance)/(dtacc), d(DBAR0)/(dtacc)
         if (print_loc) write(*,*) '*** obtain new dEVP'
         if (steptime+totaltime.eq.0) then
            mbar_dt_dot(:,:)=0d0
            DBAR0_dot(:)=0d0
         endif

c     ------------------------------------------------------
c     update jacobian
         if (print_loc) write(*,*) '*** update jacobian'
         aux66(:,:)=mbar(:,:)
         call lu_inverse(inverr,aux66,6)
         LBAR(:,:)=aux66(:,:)

c----------------------------------------------------------------------c
c     Updates DBAR0_dot, mbar_dt_dot, dtacc, deacc, drot_acc
         if (print_loc) then
            write(*,*)'*** update DBAR0_dot, mbar_dt_dot, dtacc'
         endif

         if (ifullsc) then
            if (steptime+totaltime.eq.0) then
c           trick by hard-wiring for the same values for the first increment.
               mbar_dt_sc(:,:)=mbar(:,:)*dtime
               DBAR0_sc(:)=DBAR0(:)
            endif

            if (iverbose.ge.1.and.print_loc)then
               write(*,*)
               write(*,*)
               write(*,'(a32,f8.5)')'*** updating dots... with dtacc:',
     $              dtacc
               write(*,*)'- mbar*dtime from the last full SC calc'
               do i=1,1         !6
                  write(*,'(6e13.5)')(mbar_dt_sc(i,j),j=1,6)
               enddo
               write(*,*)'- mbar for the current increment.'
               do i=1,1         !6
                  write(*,'(6f10.6)')(mbar(i,j),j=1,6)
               enddo
               write(*,*)'- mbar*dtime for the current increment.'
               do i=1,1         !6
                  write(*,'(6e13.5)')(mbar(i,j)*dtime,j=1,6)
               enddo
               write(*,*)'- mbar(i,:)*dtime-mbar_dt_sc(i,:)'
               do i=1,1         !6
                  write(*,'(6e13.5)')
     $                 ((mbar(i,j)*dtime-mbar_dt_sc(i,j)),j=1,6)
               enddo
            endif

c        Just as a reminder, dtime=dtacc in case of full SC calculation
            mbar_dt_dot(:,:)=(mbar(:,:)*dtime-mbar_dt_sc(:,:))
     $           /dtacc
            DBAR0_dot(:)=(DBAR0(:)-DBAR0_sc(:))/dtacc
            mbar_dt_sc(:,:)=mbar(:,:)*dtime
            DBAR0_sc(:)=DBAR0(:)

c        Store deformation gradient at the last self-consistent calculation
c            f_sc(:)=defgradnew(km,:)
c            u_sc(:)=stretchnew(km,:)

            do i=1,ndir
               f_sc(i)=defgradNew(km,i)
               u_sc(i)=stretchNew(km,i)
            enddo
            do i=ndir+1,ndir+nshr
               u_sc(i)=stretchNew(km,i)
            enddo
            do i=ndir+1,ndir+nshr+nshr
               f_sc(i)=defgradNew(km,i)
            enddo


            if (iverbose.ge.1.and.print_loc) then
               write(*,*)'- mbar_dt_dot'
               do i=1,1         !6
                  write(*,'(6e13.5)')(mbar_dt_dot(i,j),j=1,6)
               enddo
            endif

c        initialize accumulated strain increment.
c        Initialize accumulated time increment.
            dtacc=0d0
            deacc(:)=0d0
            drot_acc(:,:)=0d0
            do i=1,3
               drot_acc(i,i)=1d0
            enddo
         endif

         if (print_loc)write(*,*)'*** store updated dtacc deacc'
         is_initialized_statev=1
         iopt=1                 ! store updated dtacc and deacc
         include '../../../src/vpsc_src/call_statev_io.f'
         stateNew(km,:)=statev(:)
         if (print_loc) write(*,*)'*** end of loop over block'

c        determine internal energy
         dw=0d0
         do i=1,ndir
            dw=dw+stressNew(km,i)*strainInc(km,i)
         enddo
         do i=ndir+1,ndir+nshr
            dw=dw+stressNew(km,i)*strainInc(km,i)*2.
         enddo
         if (density(km).gt.0) then
            dw=dw/density(km)
         else
            dw=0d0
         endif

c$$$         if (.not.isnan(dw))then
c$$$            enerInternNew(km)=enerInternold(km)+dw
c$$$         else
c$$$            enerInternNew(km)=enerInternold(km)
c$$$         endif

c        determine elastic internal energy
c        call chg_basis(dsbar,aux33,aux66,aux3333,1,6)
c        call chg_basis( sbar,bux33,aux66,aux3333,1,6)

         dwel=0d0
         do i=1,6
            do j=1,6
               dwel=dwel+sbar(i)*dsbar(j)*selav(i,j)
            enddo
         enddo
         if (density(km).gt.0) then
            dwel=dwel/density(km)
         else
            dwel=0d0
         endif

c$$$     if (.not.isnan(dwel).and..not.(isnan(dw)))then
c$$$        enerInelasNew(km)=enerInelasOld(km)+dw-dwel
c$$$     else
c$$$        enerInelasNew(km)=enerInelasOld(km)
c$$$     endif

c     s : d^el
c     s : ds m^el
 100  continue                  ! end of km block

      if (print_loc) write(*,*)'*** end of vumat ***'

      return
      end subroutine vumat



      include "src/update_eacc_explicit.f"
      include "src/lib_write_fe.f" ! output subroutines for diagnosis
      include "src/voigt_exp.f" ! voigt notation and basis transformation for Abaqus/Explicit
      include "src/statev_io.f"     ! Subroutine to store and restore state variables
c      include "src/uvarm.f"         ! User defined output accessible from Abaqus CAE (placeholder)
      include "src/save_diffraction.f"
      include "src/hencky.f"
      include "src/voigt_abq.f" ! voigt notation and basis transformation for Abaqus/Explicit and Abaqus/Standard
      Include "src/loading_condition_fe_explicit.f"
      include "src/rodrigues.f" ! rodrigues
      include "src/twin_orientation_fe.f"
      include "src/lib_fe.f"
      include "src/calc_drot_explicit.f"
      include "src/calc_devm.f"
      include "src/statev_beta.f"
c      include "sdvini_exp.f"
      include "src/polar_dcmp"

c     subroutines shared with stand-alone (incr. E)VPSC
      include "../../../src/vpsc_src/evpsc_input_fe.f"
      include "../../../src/vpsc_src/eshelby_rot.f"
      include "../../../src/vpsc_src/elsc_evp_fe.f"
      include "../../../src/vpsc_src/chg_basis.f"
      include "../../../src/vpsc_src/eshelby.f"
      include "../../../src/vpsc_src/lu_libs.f"
      include "../../../src/vpsc_src/funcs.f"
      include "../../../src/vpsc_src/euler.f"
      include "../../../src/vpsc_src/crystal_symmetry.f"
      include "../../../src/vpsc_src/voigt.f"
      include "../../../src/vpsc_src/neighbours_fe.f"
      include "../../../src/vpsc_src/library7.sub"
      include "../../../src/vpsc_src/initial_state_guess_el_evp.f"
      include "../../../src/vpsc_src/solve_evp_bc.f"
      include "../../../src/vpsc_src/initial_state_guess_evp.f"
      include "../../../src/vpsc_src/grain_rate_and_moduli_evp.f"
      include "../../../src/vpsc_src/evpsc.f"
      include "../../../src/vpsc_src/evpsc_implicit.f"
      include "../../../src/vpsc_src/grain_stat.f"
      include "../../../src/vpsc_src/grain_stress_evp.f"
      include "../../../src/vpsc_src/newton_raphson_evp.f"
      include "../../../src/vpsc_src/calc_evptilde.f"
      include "../../../src/vpsc_src/update_macro_property.f"
      include "../../../src/vpsc_src/update_schmid_fe.f" ! contains in fact subr <update_schmid_fe>
      include "../../../src/vpsc_src/update_shape_fe.f"
      include "../../../src/vpsc_src/data_crystal_fe.f"
      include "../../../src/vpsc_src/update_crss_voce_fe.f"
      include "../../../src/vpsc_src/update_crss_rgvb_fe.f"
      include "../../../src/vpsc_src/data_grain_fe.f"
      include "../../../src/vpsc_src/data_grain_ebsd.f"
      include "../../../src/vpsc_src/select_points_from_ebsd.f"
      include "../../../src/vpsc_src/scale_gamd0s_evp_fe.f"
      include "../../../src/vpsc_src/update_fij_fe.f"
      include "../../../src/vpsc_src/update_orientation_fe.f"
      include "../../../src/vpsc_src/update_twinning_ptr_fe.f"
      include "../../../src/vpsc_src/update_twinning_vft_fe.f"
      include "../../../src/vpsc_src/stat_twinning_fe.f"
      include "../../../src/vpsc_src/write_texture_fe.f"
      include "../../../src/vpsc_src/stat_shear_activity.f"
      include "../../../src/vpsc_src/diff_planes.f"
      include "../../../src/vpsc_src/update_crss_dd_fe.f"
      include "../../../src/vpsc_src/update_crss_dd_fe_new.f"
      include "../../../src/vpsc_src/cod2grains.f"
      include "../../../src/vpsc_src/thermo_emod_cte.f"
      include "../../../src/vpsc_src/read_ebsd.f"
c     below is to bypass the issue of xit function that is not available in vumat
      include "../../../src/vpsc_src/fake_xit.f"
