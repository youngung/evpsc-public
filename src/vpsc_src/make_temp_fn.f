c     Source code from
c     https://software.intel.com/en-us/forums/topic/289743
c     Only minor modification added.

c$$$      program main
c$$$      character(72) filename
c$$$      integer exitstat
c$$$
c$$$c$$$c     mktemp
c$$$c$$$      call execute_command_line(
c$$$c$$$     $     "mktemp /tmp/XXXXXX",
c$$$c$$$     $     wait=.true.)
c$$$
c$$$c      call init_random_seed()
c$$$      do i=1,2
c$$$         call generate_temp_name('/tmp/fld0_', filename)
c$$$      enddo
c$$$      end
c----------------------------------------------------------------------c

      subroutine generate_temp_name(prefix, outname)
      implicit none
      character(*), intent(in) :: prefix
      character(*), intent(out) :: outname
      real mrand
      integer irand
      logical ex, verbose
      parameter(verbose=.false.)

      do
! Get a random number
      call random_number(mrand)
! Get the low 24 bits
      irand = transfer(mrand, irand)
      irand = iand(irand,Z'FFFFFFFFF')
      write (outname,101) trim(prefix),irand
 101  format (A,Z8.7,'.tmp')
      inquire (file=trim(outname),exist=ex)
      if (.not.ex) then
         open(8888,file=trim(outname))
         close(8888)
         exit
      endif
      end do
      outname=trim(outname)
      if (verbose) print *, trim(outname)

! We have a name that does not exist (at least right now)
      return
      end subroutine generate_temp_name
