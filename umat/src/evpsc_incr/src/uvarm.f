      SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     $     NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
     $     JMAC,JMATYP,MATLAYO,LACCFLA)
      implicit none
      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3 FLGRAY(nuvarm)
      DIMENSION UVAR(NUVARM)
      dimension DIRECT(3,3),T(3,3),TIME(2)
      DIMENSION ARRAY(nuvarm),JARRAY(nuvarm),JMAC(*),JMATYP(*),COORD(*)
c     element number bottom and up.
      integer noel_bot,noel_top
      common/ktest/noel_bot,noel_top

C     The dimensions of the variables FLGRAY, ARRAY and JARRAY
C     must be set equal to or greater than 15.

c     local
      include "../../../../src/vpsc_src/incr_evpsc_fe.dim" ! for mapping statev to EVPSC arrays
      integer iopt,iv0,iv1,igb,igt,iphel,i
      include "src/io_head.f"
      include "../umat-par.f"
      include "../umat-io.f"

      if (nuvarm.ne.30) then
         write(*,*)'-------------------------------------------------'
         write(*,*)'** NUVARM should be 30 (hardwired in <io_head.f>)'
         write(*,*)'** Consider changing nuvarm in your Abaqus job.'
         write(*,*)'**  Look at <USER OUTPUT VARIABLES>'
         write(*,*)'-------------------------------------------------'
         write(*,*)
         write(*,*)
         call xit()
      endif

      if (nuvarm.ne.nuvarmx) then
         write(*,*)'-------------------------------------------------'
         write(*,*)'** NUVARMMX should be equal to nuvarm            '
         write(*,*)'** Consider changing nuvarm in your Abaqus job.  '
         write(*,*)'**  Look at <USER OUTPUT VARIABLES>              '
         write(*,*)'-------------------------------------------------'
         write(*,*)
         write(*,*)

         call xit()
      endif

      do i=1,nuvarm
         uvar(i)=uvarm_custom(noel,npt,i)
      enddo

      RETURN
      END
