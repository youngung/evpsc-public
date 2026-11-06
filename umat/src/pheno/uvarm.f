      SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     $     NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
     $     JMAC,JMATYP,MATLAYO,LACCFLA)
      include 'ABA_PARAM.INC'
      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3 FLGRAY(nuvarm)
      DIMENSION UVAR(NUVARM)
      dimension DIRECT(3,3),T(3,3),TIME(2)
      DIMENSION ARRAY(nuvarm),JARRAY(nuvarm),JMAC(*),JMATYP(*),COORD(*)
c      real*8 uvar(nuvarm)

C     The dimensions of the variables FLGRAY, ARRAY and JARRAY
C     must be set equal to or greater than 15.

c     local
c$$$      include "../../../../src/vpsc_src/incr_evpsc_fe.dim" ! for mapping statev to EVPSC arrays
c$$$c      integer nstatv
c$$$      include "../umat-nstv.f"
c$$$      real*8 statev(nstatv)
c$$$      integer iopt,iv0,iv1,igb,igt,iphel,i
c$$$      include "src/io_head.f"
c$$$      include "../umat-par.f"
c$$$      include "../umat-io.f"
c$$$
c$$$      if (nuvarm.ne.20) then
c$$$         write(*,*)'** NUVARM should be 20'
c$$$         call xit()
c$$$      endif
c$$$
c$$$
c$$$      if (nuvarm.ne.nuvarmx) then
c$$$         write(*,*)'** NUVARMMX should be 20'
c$$$         call xit()
c$$$      endif
c$$$
c$$$      do i=1,nuvarm
c$$$         uvar(i)=uvarm_custom(noel,npt,i)
c$$$      enddo

c      if (noel.eq.1..and.npt.eq.1) then
c         write(*,*)'uvar for ',noel,'npt',npt
c         write(*,'(4e20.10)')(uvar(i),i=1,4)
c      endif

c      write(*,*)'EPSTOT from SDV in subr. UVARM'
c      write(*,*) epstot(1,1),epstot(2,2),epstot(3,3)

c      call uvarm_rgvb()
c      call uvarm_activ(field,nfield,iv0,iv1,statev,nstatv)

      RETURN
      END
