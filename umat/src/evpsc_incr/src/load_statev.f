      if (ipxx.eq.0) then
         inquiry_exist=.false.
      elseif (ipxx.eq.1.or.ipxx.eq.2) then
         inquire(file=OUTDIR(1:LENOUTDIR)//'/statev_px',
     $        exist=inquiry_exist)
      endif

      if (.not.inquiry_exist) then
         if (print_loc) write(*,*)'** load statev'

c     ------------------------------------------
c     call mutexlock(3)

         if (print_loc) write(*,*)"** getting in evpsc_input_fe"
         open(action='read',unit=UR_INP,
     $        file=OUTDIR(1:LENOUTDIR)//'/vpsc7.in',status='old')
         include "../../../src/vpsc_src/call_evpsc_input_fe.f"
         close(UR_INP)
         if (print_loc) write(*,*)'** end evpsc_input_fe'
c         call mutexunlock(3)

c     ------------------------------------------
c        initialization of rotation matrix
         rot_acum(:,:)=0.
         do i=1,3
            rot_acum(i,i)=1d0
         enddo
c        initialization of kount_umat_call
         kount_umat_call = 0
         acctimefe=0.

c        PX level states
         svm=0.
         epsvm=0.
         sbar(:)=0.
         dsbar(:)=0.
         sbar_old(:)=0.
         DBAR0(:)=0.
         work(:)=0.
         aa(:)=0.              !! accumulated activity
         ac(:)=0.              !! accumulated activity
         alplpx=0.

c        Grain level states
         sg(:,:)=0.
         dsg(:,:)=0.
         sg_old(:,:)=0.
         sgdot(:,:)=0.
         gtotgr(:)=0.
c            grsze(:)=0.

c        Slip system level
         tau0(:,:)=0.          ! caution

         deref=-1.


c     *** IF idiff=1, Reads diffraction directions and crystal planes along
c     them for calculating average elastic strain on those planes.
         if (idif_calc) then
            call DIFF_PLANES_PX(
c     intent(in)::
     $           0,idiff,UR_DIF,-1,1,NPH,ngrmx,nphmx,ndifmx,ngr,ag,
     $           wgt,filediff,ICRYSYMph,iregime,iDBARv,aux33,aux6,
     $           BC_EL,mgel,aux33,wph,
c     intent(inout)::
     $           ndiff,sg,angdetector,nfamily,diff_vs,diff_vc,NGRSETPH,
     $           WGTSETPH,IGRSET,
c     intent(out)::
     $           eelgr)
c-----------------------------------------------------------------------
c     *** Place holder for NIST diffraction experiments
c         if (idiff(iphel).eq.2) call ....
c           call diff_planes for NIST bench mark
c-----------------------------------------------------------------------
c$$$               endif
c$$$            enddo
         endif
      else !! This feature was mostly for Mooyeong's 3-point bending FE simulation.
c              with strips that have undergone prestrains, during which texture
c              and dislocation density evolved.
         write(*,*)'------------------------------------------'
         if (ipxx.eq.1) then
            write(*,*)'** Restoring statev from statev_px file **'
         elseif(ipxx.eq.2) then
            write(*,*)'**  Restoring statev from px list file  **'
         endif
         write(*,*)'------------------------------------------'

         UR_STV=UR_STV
         write(*,*)'UR_STV:',UR_STV
         if (ipxx.eq.1) then
c     Load entire state variables from a single file
            open(unit=UR_STV,
     $           file=OUTDIR(1:LENOUTDIR)//'/statev_px')
c        load statev from a file

         elseif(ipxx.eq.2) then

            !! Find the proper file name from the list
            open(unit=UR_STV,
     $           file=OUTDIR(1:LENOUTDIR)//'/statev_px')

            !! flush away
!            npt_per_el = 8 ! for example for element type with 8 integration points
!            do i=1,(noel-1)*npt_per_el+npt
            do i=1,noel
               read(UR_STV,'(a)')
            enddo
            !! read the correct one
            read(UR_STV,'(a)')statev_fn
            close(UR_STV)

            !! open the file...
            open(unit=UR_STV,file=statev_fn)
         endif


         iopt=4
         call mutexlock(3)
         include "../../../src/vpsc_src/call_statev_io.f"
         close(unit=UR_STV)



         pdil=0d0
         CALL ESHELBY(AUX3,AUX3333,0d0,AUX3333,AUX3333,AUX33,AUX33
     $        ,PDIL,AUX3333,AUX3333,0,iesh_fail)
         call mutexunlock(3)


c        load dEVPSC array from statev
         iopt=2
         include "../../../src/vpsc_src/call_statev_io.f"

c     ------------------------------------------
c     Initialization for stress-extrapolation (assumes that the material has been unloaded)
         acctimefe=0d0
         sbar_old(:)=0d0
         sg_old(:,:)=0d0
c        initialization of rotation matrix
         rot_acum(:,:)=0d0
         do i=1,3
            rot_acum(i,i)=1d0
         enddo
c        initialization of kount_umat_call
         kount_umat_call = 0
c        ------------------------------------------

         ! epsvm=0d0 maybe, it's better to use evpsc from state_px for the sake of recording the entire history.
         dsbar(:)=0d0
         gtotgr(:)=0d0
c        grsze(:)=0d0
c        tau0(:,:)=0d0          ! caution
         aa(:)=0d0              !! accumulated activity
         ac(:)=0d0              !! accumulated activity
         iopt=1
         include "../../../src/vpsc_src/call_statev_io.f"
      endif

      do iphel=1,NPH
         if (ishape(iphel).ge.2) then
            write(*,*)'Warning!!'
            write(*,*)'Grains are assumed to have the same shape and'
            write(*,*)'orientations for reducing memory'
            write(*,*)'asgr and axisgr are not stored in statev'
            call xit()
         endif
      enddo


      if (print_loc) write(*,*)'** save dEVPSC arrays to statev'
      is_initialized_statev=1
      iopt=1
      include "../../../src/vpsc_src/call_statev_io.f"
