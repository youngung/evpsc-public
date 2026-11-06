c     subroutine sdvini to initialize state variables of dEVPSC-FE
c     interfaced with Abaqus/standard solver via subroutine evpsc_input_fe
c
c     yjeong@changwon.ac.kr
c     One should use <*initial conditions, type=solution, user> keyword
c     at the beginning of STEP in Abaqus input file.
      SUBROUTINE SDVINI(STATEV,COORDS,NSTATV,NCRDS,NOEL,NPT,LAYER,KSPT)
C
      INCLUDE 'ABA_PARAM.INC'
c #include <SMAAspUserSubroutines.hdr>
C
      DIMENSION STATEV(NSTATV),COORDS(NCRDS)
c     locals
      character*256 OUTDIR
      INTEGER lenoutdir
      real*8 dtime,aux3(3),aux3333(3,3,3,3),pdil,t0,t10,aux6(6),
     $     aux33(3,3),vnorm
      logical inquiry_exist,print_loc,idif_calc,iesh_fail
c      parameter(print_loc=.true.)
      parameter(print_loc=.false.)
      parameter(idif_calc=.true.)
c$$$      include "../../../src/vpsc_src/incr_evpsc_fe.dim"
c$$$      include "umat-par.f"
c$$$      include "src/io_head.f"
c$$$      include "umat-io.f"
c$$$
c$$$      call mutexinit(1)
c$$$      call mutexinit(3)
c$$$      if (iverbose.ge.1) then
c$$$         write(*,*)
c$$$         write(*,*)'------------------------------------------'
c$$$         write(*,*)'** Reading polycrystal states variables **'
c$$$         write(*,*)'**   via subroutine sdvini              **'
c$$$         write(*,*)'------------------------------------------'
c$$$         write(*,*)
c$$$      endif
c$$$
c$$$c      write(*,*)'** begins sdvini, noel,npt:',noel,npt
c$$$
c$$$      call getoutdir(outdir, lenoutdir)
c$$$      UR_INP_vpscinp=102+((noel-1)*12+(npt-1))
c$$$      UR_SX_sx=UR_INP_vpscinp+1
c$$$      UR_TX_tx=UR_SX_sx+1
c$$$      UR_SHP_fa=UR_TX_tx+1
c$$$      UR_STV=UR_SHP_fa+1
c$$$      UR_DIF_fe=UR_STV+1
c$$$      UW_LOG_fe=UR_DIF_fe+1
c$$$
c$$$
c$$$      call cpu_time(t0)
c$$$      include"src/load_statev.f"
c$$$      call cpu_time(t10)
c$$$
c$$$      if (iverbose.ge.1) then
c$$$         do iphel=1,NPH
c$$$            open(uwfe01+iphel-1,file=OUTDIR(1:LENOUTDIR)//
c$$$     $           '/TEX_PH_EL'//
c$$$     $           CHAR(48+iphel)//'.OUT',status='unknown')
c$$$         enddo
c$$$      endif
c$$$
c$$$c      write(*,'(a32,f13.6,a12,2(i4.4,x))')
c$$$c     $     '** SDVINI, elapsed time [sec]:',t10-t0,
c$$$c     $     '(noel,npt):',noel,npt

      RETURN
      END subroutine sdvini
