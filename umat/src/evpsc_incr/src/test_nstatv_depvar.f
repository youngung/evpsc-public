C     This program is to update *DEPVAR keyword in *.inp files
c     found in the working directory.
C
      program get_no_depvar
      implicit none
c     ----------------- arrays from common blocks shared with stand-alone evpsc
      integer nstatv,iofn,nf,i,nvar,iopt
c     command line argument
      integer narg,isolver_type ! (0: Abaqus/standard, 1: Abaqus/Explicit)
      parameter(narg=2)
      real*8 xmass_scale
      integer nstatevmx,istatus(narg)
      character*12 com_buf(narg)
c     --
      parameter(nstatevmx=10000000)
      real*8 statev(nstatevmx)
      character*125 nstatv_char

      include "../../../../src/vpsc_src/incr_evpsc_fe.dim"

c     1st argument (0: Abaqus/Standard input;   1: Abaqus/Explicit input))
c     2nd argument (mass scaling factor - applicable only for Abaqus/Explicit)
      istatus(:)=0
      do i=1,iargc()
         istatus(i) = 1
         call get_command_argument(i,com_buf(i))
      enddo

c     Default solver type is 0 : Abaqus/Standard
      isolver_type=0
      if (istatus(1).eq.1) then
         read(com_buf(1),'(i10)') isolver_type
      endif

      if (istatus(1).eq.1.and.isolver_type.eq.1) then
c     write(*,'(a)')'** Test runs for Abaqus/Explicit input'
c     write(*,*)'mass scale factor:',com_buf(2)
         read(com_buf(2),*) xmass_scale
         write(*,'(a21,x,e10.3)')'** Mass scale factor:',xmass_scale
      else
         write(*,*)'** Test runs for Abaqus/Standard input'
      endif
      !write(*,*)"======================================================"
c      write(*,*)"        First run call statev to determine"
c      write(*,*)"       the minimum *nstatv* (a.k.a. depvar)"
      nstatv=0
      iopt=0
      include "../../../../src/vpsc_src/call_statev_io.f" ! this run is to determine nvar
      write(nstatv_char,*)nvar
      iopt=1
      include "../../../../src/vpsc_src/call_statev_io.f" ! this run is to determine nvar
      write(nstatv_char,*)nvar
      iopt=2
      include "../../../../src/vpsc_src/call_statev_io.f" ! this run is to determine nvar
      write(nstatv_char,*)nvar
      nstatv=0
      iopt=0
      include "../../../../src/vpsc_src/call_statev_io.f" ! this run is to determine nvar
      write(nstatv_char,*)nvar
      iofn=1
c      write(*,*)'------------------------------------------------------'

      if (isolver_type.eq.0)then
         call system("python diag.py -n "//trim(nstatv_char)
     $        //" -i 0")
      elseif(isolver_type.eq.1)then
         call system("python diag.py -n "//trim(nstatv_char)
     $        //" -i 1 -msf "//com_buf(2))
      endif
      open(10,file='list_abq_inputs.txt',status='old')
      read(10,*) prosa
      read(10,*) nf
      ! write(*,*)'-------------------------------------------------'
      ! write(*,*)'** List individual input files'
      do i=1,nf
         read(10,*) nstatv,prosa
         iopt=0
         include "../../../../src/vpsc_src/call_statev_io.f"
      enddo
      end program get_no_depvar


c     ----
c     This is a fake subroutine.
c     This subroutine is required when call statev_io
      subroutine xit
      return
      end subroutine xit
