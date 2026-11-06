c     lop
c      0: At the start of analysis
c      1: At the start of current anlaysis increment
c      2: At the end   of current anlaysis increment
c      3: At the end   of analysis
c      4: At the start of restart
      subroutine uexternaldb(lop,lrestart,time,dtime,kstep,kinc)
c
      implicit none
c
c#include <SMAAspUserSubroutines.hdr>
      dimension time(2)
c     cod common block
      integer nphimx,nphi1mx,nphi2mx
c     cubic without sample symmetry using the angular resolution of 5.
c     phi1: 72+1=360/5+1
c     phi : 18+1=90/5+1
c     phi2: 18+1=90/5+1
      include "../../../src/vpsc_src/incr_evpsc_fe.dim"

      parameter(nphi1mx=72+1,nphimx=18+1,nphi2mx=18+1)
      real*8 cod(nphi1mx,nphimx,nphi2mx)
      logical is_cod_avail(5), does_ebsd_f_exist
      real*8 phi1mx,phimx,phi2mx,res
      common/kcod/cod,phi1mx,phimx,phi2mx,res,is_cod_avail
      integer t1,t0_incr,t0_anlys,t0_step
      common/ktime/t0_incr,t0_anlys,t0_step
      real*4 f_scale
      integer nlines_ebsd,nrowmx, iodf, lenoutdir
      character*256 ebsd_filename, outdir
      parameter(iodf=0)

c     block dedicated to texture input directly from EBSD map in *.ang form (iodf=3)
      parameter(ebsd_filename='c:/tmp/ebsd_raw.ang')
      parameter(nrowmx=200000)
      real*4 ebsd_raw_data(5,nrowmx)
      common/kebsd/ebsd_raw_data,nlines_ebsd

c     this block is to figure out the beginning (noel_bot)
C     and the end (noel_up) of element numbers
      integer noel_bot,noel_up
      common/ktest/noel_bot,noel_up
      real*8 alplstats(6,2)     ! min and max of "alplpx", "RELS", "RELD", "it2", "EPSVM", "SVM"
      common/kstats/alplstats
c     eshelby initiation
      real*8 aux3(3),aux6(6),aux3333(3,3,3,3),aux33(3,3),aux66(6,6),PDIL
      logical iesh_fail
c     end of cod common block

c
c
c     user coding to set up the Fortran environment, open files, close files,
c     calculate user-defined model-independent history information,
c     write history information to external files,
c     recover history information during restart analyses, etc.
c     do not include calls to utility routine XIT
c
c
      call mutexinit(1)
      call mutexinit(2)
      call mutexinit(3)
      call mutexinit(4)
      call mutexinit(5)
      call mutexinit(6)

      call getoutdir(outdir, lenoutdir)


c     At the beginning of a regular analysis or a restart job,
c     initialize arrays required in Eshelby and chg_basis subroutines.
      if (lop.eq.0 .or. lop.eq.4) then
C        INITIALIZE GAUSS-LEGENDRE OR GAUSS-LOBATO COORDINATES AND ARRAYS USED
C        IN THE DOUBLE INTEGRALS GIVING THE ESHELBY TENSORS
         pdil=0d0
         CALL ESHELBY(AUX3,AUX3333,0d0,AUX3333,AUX3333,AUX33,AUX33,PDIL,
     $        AUX3333,AUX3333,0,iesh_fail)
c        Make sure that b-basis converting matrix is initilized.
         call chg_basis(aux6,aux33,aux66,aux3333,0,6)
      endif

      if (lop.eq.0) then
         write(*,*)
         write(*,*)
         write(*,*)'**************************************************'
         write(*,*)'* delta/sigma EVPSC model interfaced with Abaqus *'
         write(*,*)'**************************************************'
         write(*,*)"*      Youngung Jeong@Changwon Nat'l Univ.       *"
         write(*,*)'*                 yjeong@changwon.ac.kr          *'
         write(*,*)'**************************************************'
         write(*,*)'*               Start of analysis                *'
         write(*,*)'**************************************************'
         write(*,*)
         write(*,*)


         noel_bot=1000000
         noel_up=0

         call system_clock(t0_anlys)
         is_cod_avail(:)=.false.


         if (iodf.eq.3) then
            write(*,*)
            write(*,*)
            write(*,*)'***********************************************'
            write(*,*)'** IODF of 3 was chosen, raw EBSD ang data file'
            write(*,*)'    is used. The ebsd_filename: '
            write(*,*) outdir(1:lenoutdir)//'/'//trim(ebsd_filename)
            f_scale=0.001       !! FE usually follows mmSI unit, while
c                                ! An "ang file" makes use of 'micro mm'.
            UR_EBSD=300             !! hard-wired file ID chosen by Jaeseong.
            inquire(file=outdir(1:lenoutdir)//'/'//trim(ebsd_filename),
     $           exist=does_ebsd_f_exist)
            if (.not.(does_ebsd_f_exist)) then
               write(*,*)'** Requested ebsd file <',
     $              outdir(1:lenoutdir)//'/'//trim(ebsd_filename),
     $              '> could not be found.'
               call xit()
            endif
            OPEN(UR_EBSD, FILE=outdir(1:lenoutdir)//'/'//
     $           trim(ebsd_filename),STATUS='UNKNOWN')

            call mutexlock(1)
            call read_ebsd(UR_EBSD,f_scale,ebsd_raw_data,nlines_ebsd)
            call mutexunlock(1)
            ! write(*,*)'nlines_ebsd:',nlines_ebsd
            CLOSE(unit=UR_EBSD)
            write(*,*)'***********************************************'
            write(*,*)
            write(*,*)
         endif
      elseif(lop.eq.1) then
         write(*,*)'** At the start of current analysis increment'
         call system_clock(t0_incr)
         alplstats(:,1)= 1e10   !min
         alplstats(:,2)=-1e10   !max
      elseif(lop.eq.2) then
         write(*,*)'** At the end   of current analysis increment'
         call system_clock(t1)
         write(*,"(a22,x,i10,a6,i2,a6,i5)")
     $        '*** elapsed time incr:',
     $        t1-t0_incr,'kstep',kstep,'kinc',kinc
         write(*,"(6(2a9))")'mn alpl','mx alpl',
     $        'mn RELS','mx RELS',
     $        'mn RELD','mx RELD',
     $        'mn it2','mx it2',
     $        'mn EVM','mx EVM',
     $        'mn SVM','mx SVM'
         write(*,"(6(2e9.2))")((alplstats(j,i),i=1,2),j=1,6)
      elseif(lop.eq.3) then
         write(*,*)'** At the end   of analysis'
         call system_clock(t1)
c$$$         write(*,"(a22,x,i10,a5,f10.2,a5)")
c$$$     $        '*** elapsed time anly:',
c$$$     $        t1-t0_anlys,'secs',(t1-t0_anlys)/3600.,'hrs'
      elseif (lop.eq.4) then
         write(*,*)
         write(*,*)
         write(*,*)'******************************************'
         write(*,*)'*delta-EVPSC model interfaced with Abaqus*'
         write(*,*)'******************************************'
         write(*,*)"*   Youngung Jeong@Changwon Nat'l Univ.  *"
         write(*,*)'*                 yjeong@changwon.ac.kr  *'
         write(*,*)'******************************************'
         write(*,*)'*           Start of restart analysis    *'
         write(*,*)'******************************************'
         write(*,*)
         write(*,*)


c$$$      elseif(lop.eq.5) then
c$$$         call system_clock(t0_step)
c$$$      elseif(lop.eq.6) then
c$$$         call system_clock(t1)
c$$$c$$$         write(*,"(a22,x,f10.2)")
c$$$c$$$     $        '--- elapsed time step:',(t1-t0_step)*1.0
c$$$c$$$         write(*,'(a38)')'--------------------------------------'
      endif

      return
      end subroutine uexternaldb
