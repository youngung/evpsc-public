c-----------------------------------------------------------------------
c     This subroutine is to check the format of single crystal
c     file used in dEVPSC.

c     1) In certain types of single crystal files, the thermal
c     expansion coefficient and the elastic constants are given
c     as polynomial functions. The coefficients of such polynomial
c     are given in file, and the file name of that coefficients may
c     be given in the single crystal, which makes the format of
c     single crystal different than the usual format.
c
c     2) In other types of single crystal files, non-schmid effect
c     is accounted for where non-schmid factors are listed.
c
c
c     ** Table of isx_ver **
c     ------------------------------------------------------------------
c      isx_ver   CTE/EMOD file  CTE/EMOD file exists   non-schmid factor
c     ------------------------------------------------------------------
c         0           x                   -                   x
c         1           o                   o                   x
c         2           o                   x                   x
c         3           x                   -                   o
c         4           o                   o                   o
c         5           o                   x                   o
c     ------------------------------------------------------------------

      subroutine sx_version_tester_nsch(
c     intent(in)
     $     UR_SX
c     intent(out)
     $     ,isx_ver
     $     )
      implicit none
      integer, intent(in):: UR_SX
      integer, intent(out):: isx_ver
c     locals
      integer i,nline
      character*128 prosa
      real*8 aux6(6)

      call sx_version_tester_old(UR_SX,isx_ver)

      if (isx_ver.eq.0) then
c     ** without thermo-dependent CTE/EMOD file line'
         nline=13
      elseif (isx_ver.eq.1) then
c     ** with thermo-dependent CTE/EMOD file line'
         nline=15
      elseif (isx_ver.eq.2) then
c     ** with thermo-dependent CTE/EMOD file line but without actual file'
         nline=15
      endif
      rewind(UR_SX)

      do i=1,nline
         read(UR_SX,'(a)')prosa
      enddo
      read(UR_SX,*,err=20) (aux6(i),i=1,6)
      goto 30

c     block of statements for sx without nsch factors.
 20   rewind(UR_SX)
      return

 30   isx_ver=isx_ver+3
      rewind(UR_SX)
      return
      end subroutine sx_version_tester_nsch

c     ----------------------------------------------------------------
      subroutine sx_version_tester_old(
c     intent(in)
     $     UR_SX
c     intent(out)
     $     ,isx_ver
     $     )
      implicit none
      integer, intent(in):: UR_SX
      integer, intent(out):: isx_ver
c     locals
      real*8 aux6(6)
      character(len=256) fn_temp_emod_cte,prosa
      integer i
      logical inquiry_exist

c     ----------------------------------------------------------------
c     Determine the version of single crystal file (isx_ver);
c     See if filecrys includes 'thermal_crystal_property.dat' file
      fn_temp_emod_cte='None'
      isx_ver=-1
      inquiry_exist=.false.
      do i=1,4
         read(UR_SX,'(a)') prosa  !! skip next four lines
      enddo
c     See if the next line is the first 6 elastic constants.
      do i=1,6
         read(UR_SX,*,err=14) aux6(i)
      enddo

c     FNSX is in its usual file format without cte/emod file.
c     so that the block between 14-16 are skipped.
      isx_ver=0

      goto 16
c     ------------------------------------------------------
c     This block of statements are relevant only for single
c     crystal files with the <fn_temp_emod_cte> is given.
 14   isx_ver=1
      rewind UR_SX
      do i=1,4
         read(UR_SX,'(a)') prosa  !! skip next four lines
      enddo
      read(UR_SX,'(a)') fn_temp_emod_cte
      inquire(file=fn_temp_emod_cte,exist=inquiry_exist)

      if (.not.inquiry_exist) then
         write(*,*)'** Warning: ignored <fn_temp_emod_cte> file'
         write(*,*)'   The file <',trim(fn_temp_emod_cte),'>'
         write(*,*)'   could not be found.'
         isx_ver=2
      else
         write(*,*)'** The file below'
         write(*,*) trim(fn_temp_emod_cte)
         write(*,*)' exists indeed.'
      endif
c     ------------------------------------------------------
 16   rewind UR_SX

      return
      end subroutine sx_version_tester_old
