c     subroutine sdvini to initialize state variables of dEVPSC-FE
c     interfaced with Abaqus/standard solver via subroutine evpsc_input_fe
c
c     yjeong@changwon.ac.kr
c     One should use <*initial conditions, type=solution, user> keyword
c     at the beginning of STEP in Abaqus input file.

c     (2023-September)
c     An addditional block addded to hardwire a location-dependent rotation
c     of crystallographic texture for an extruded Mg10Gd sample for Jaeseong Lee
      SUBROUTINE SDVINI(STATEV,COORDS,NSTATV,NCRDS,NOEL,NPT,LAYER,KSPT)
C
      implicit none
c #include <SMAAspUserSubroutines.hdr>
C
      DIMENSION STATEV(NSTATV),COORDS(NCRDS)
      real*8 ebsd_rad
c     locals
      character*256 OUTDIR,statev_fn
      character*9 cnoel,ckspt
      INTEGER lenoutdir
      real*8 dtime,aux3(3),aux3333(3,3,3,3),pdil,t0,t10,aux6(6),
     $     aux33(3,3),vnorm,theta, R(3,3),agr(3,3),pi
      logical inquiry_exist,print_loc,idif_calc,iesh_fail
c      parameter(print_loc=.true.)
      parameter(print_loc=.false.)
      parameter(idif_calc=.true.)
      integer iodf
      parameter(iodf=1)
      integer ipxx
      parameter(ipxx=0) ! 0: no state px; 1: single unique statepx ; 2: list of single_px
      integer ife
      parameter(ife=1) ! 0: stand-alone; 1: FE

c----------------------------------------------------------------------c
      integer ndi,kkk,j,k
c----------------------------------------------------------------------c
c     Jaeseong's EBSD loader
      integer nlines_ebsd,nrowmx
      parameter(nrowmx=200000)
      real*4 ebsd_raw_data(5,nrowmx)
      common/kebsd/ebsd_raw_data,nlines_ebsd
      integer noel_bot,noel_top
      common/ktest/noel_bot,noel_top
c----------------------------------------------------------------------c
c     Auxilary arrays required to initialize the b-basis matrices
c     in subr. chg_basis
      real*8 aux66(6,6)         !,aux3333(3,3,3,3),aux6(6),aux33(3,3),
c----------------------------------------------------------------------c

      include "src/io_head.f" !! read nelmx; determine if nelmx is large enough.
      include "../../../src/vpsc_src/incr_evpsc_fe.dim"
      include "umat-par.f"
      ndi=-1 !! hardwiring ndi=-1 just to pass through statements in <umat-io.f> without warning
      include "umat-io.f"

c      write(*,*)'beginning of subr. sdvini'
      ebsd_rad=50e-3            ! (as for mm-SI unit) 100e-3 = 0.1 mm = 100 um
      call chg_basis(aux6,aux33,aux66,aux3333,0,6) !! do not remove this line

      write(cnoel,'(i7.7)') noel
      write(ckspt,'(i7.7)') kspt

      call mutexlock(1)
      if (noel.lt.noel_bot) noel_bot=noel
      if (noel.gt.noel_top) noel_top=noel
      call mutexunlock(1)

      if (iverbose.ge.1.and..false.) then
         write(*,*)
         write(*,*)'------------------------------------------'
         write(*,*)'** Reading polycrystal states variables **'
         write(*,*)'**   via subroutine sdvini              **'
         if (iodf.eq.3) then
            write(*,*)'** Running iodf.eq.3 option by reading'
            write(*,*)'   texture directly from EBSD ang file'
         endif
         write(*,*)'------------------------------------------'
      endif

      call getoutdir(outdir, lenoutdir)
      call cpu_time(t0)
      if (print_loc) then
         write(*,*)'** Begins load_statev'
      endif
      include"src/load_statev.f"
      if (print_loc.ge.1) then
         write(*,*)'** ended load_statev'
      endif
      call cpu_time(t10)


c     --------------------------------------------------------------------
c$$$c     begins applying the proper transformation to reorient Mg10Gd texture
c$$$c     for Jaeseong's C-ring.
c$$$      if (.false.) then
c$$$c     if (.true.) then
c$$$      pi=4.0*atan(1.0)
c$$$      theta=atan2(coords(2),coords(1))
c$$$      R(:,:)=0.
c$$$      R(1,1)=cos(theta)
c$$$      R(1,2)=cos( 90.*pi/180. - theta)
c$$$      R(2,1)=cos( 90.*pi/180. - theta)
c$$$      R(2,2)=cos(180.*pi/180. - theta)
c$$$      R(3,3)=-1.
c$$$
c$$$c$$$      if ((noel.eq.3 .or. noel .eq. 1803).and.npt.eq.1) then
c$$$c$$$         write(*,*) 'theta:', theta*180./pi
c$$$c$$$         write(*,*) 'noel:',noel
c$$$c$$$         write(*,*)'R matr'
c$$$c$$$         do i=1,3
c$$$c$$$            write(*,'(3f5.2)')(R(i,j),j=1,3)
c$$$c$$$         enddo
c$$$c$$$      endif
c$$$
c$$$      do iphel=1,NPH
c$$$         do kkk=ngr(iphel-1)+1,ngr(iphel)
c$$$            agr(:,:)=0.
c$$$            do i=1,3
c$$$               do j=1,3
c$$$                  do k=1,3
c$$$                     AGR(i,j)=AGR(i,j)+R(i,k)*AG(k,j,kkk)
c$$$                  enddo
c$$$               enddo
c$$$            enddo
c$$$c$$$            if ((noel.eq.3 .or. noel .eq. 1803).and.npt.eq.1
c$$$c$$$     $           .and.kkk.eq.1) then
c$$$c$$$               write(*,*)'ag matr'
c$$$c$$$               do i=1,3
c$$$c$$$                  write(*,'(3f5.2)')(ag(i,j,kkk),j=1,3)
c$$$c$$$               enddo
c$$$c$$$               write(*,*)'agr matr'
c$$$c$$$               do i=1,3
c$$$c$$$                  write(*,'(3f5.2)')(agr(i,j),j=1,3)
c$$$c$$$               enddo
c$$$c$$$            endif
c$$$            ag(:,:,kkk)=agr(:,:)
c$$$         enddo
c$$$      enddo
c$$$      iopt=1 !! store dEVPSC arrays to state variables
c$$$      include "../../../src/vpsc_src/call_statev_io.f"
c$$$      endif
c     --------------------------------------------------------------------


      if (iverbose.ge.1) then
         do iphel=1,NPH
c$$$            write(*,*)'uwfe01+iphel-1 in sdvini:',uwfe01+iphel-1
            open(action='write',
     $           unit=uwfe01+iphel-1,file=OUTDIR(1:LENOUTDIR)//
     $           '/TEX_PH'//CHAR(48+iphel)//
     $           '_EL_'//trim(cnoel)//
     $           '_KT_'//trim(ckspt)//'.OUT',status='unknown')
         enddo
      endif

      if(npt.eq.1.and.kspt.eq.1.and.mod(noel,100).eq.0.and..false.)then
         write(*,'(a32,f7.2,a18,i7.7,x,i2.2,x,i2.2,x,a7,f9.3)')
     $        '** SDVINI, elapsed time [ms]:',(t10-t0)*1e3,
     $        '(noel,npt,kspt):',noel,npt,kspt,'epsvm:',epsvm
      endif

      if (print_loc) then
         write(*,*)'end of sdvini'
      endif

      return
      end subroutine sdvini
