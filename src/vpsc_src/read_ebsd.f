c---------------------------------------------------------------
c     Arguments
c     ---------
c     ur: file id (300) used in <evpsc_main> and <uexternaldb.f>,
c              for SA and FE versions.
c     f_scale:   is usually set to 0.001.
c               ** An FE model in Abaqus usually follows the unit system of mmSI,
c                 while an "ang file" makes use of 'micro mm'.c
c     ebsd_map: ebsd_map of phi1, Phi, phi2, X, Y
c     nlines_ebsd: number of data lines written in *.ang file.
c                this also corresponds to the meaningful number of rows in
c                <ebsd_map> array.
      SUBROUTINE read_ebsd(
c          in
     $     ur,f_scale,
c          out
     $     ebsd_map,nlines_ebsd)

      implicit none

      real*4, intent(in):: f_scale
      integer, intent(out):: nlines_ebsd
      integer nrowmx
      parameter (nrowmx=200000)
      real*4, intent(out):: ebsd_map(5,nrowmx)
c     local variables
      integer ur, i, j, no_headlines !=kount
      character*120 line
      real*8 phi1, PHI, phi2, x, y
      logical iexit

c---------------------------------------------------------------
c     1. Count the number of headlines that are not necessary
c---------------------------------------------------------------
      iexit=.false.
      no_headlines=0
      do while(.not.iexit)
         no_headlines=no_headlines+1
         read(ur,*,iostat=i) phi1,phi,phi2,x,y
         if (i.eq.0) then
            read(ur,'(a)')line
            no_headlines=no_headlines-1
            iexit=.true.
         endif
         if (no_headlines.gt.100) then
            write(*,*)'** Error, maybe EBSD data file is not in a'
            write(*,*)'   correct form.'
            call xit()
         endif
      enddo
      rewind(ur)

c---------------------------------------------------------------
c     2. Check the maximum number of row
c---------------------------------------------------------------
      do i=1,no_headlines
         read(ur,'(a)') line
      enddo

      nlines_ebsd = 0
      do
         read(ur,*,iostat=i) phi1, PHI, phi2, x, y
         if (i.NE.0) exit
         nlines_ebsd = nlines_ebsd + 1
      end do

      write(*,'(a35,i3)')'** Info: No. of head lines:',no_headlines
      write(*,'(a35,i7)')'** Info: No. of data lines:',nlines_ebsd
      write(*,'(a35,i7)')'** Info: No. of max  lines:',nrowmx

c     Check the nrowmx is reasqonable
      IF(nrowmx.LT.nlines_ebsd) THEN
         write(*,*)'WARNING: NROWMX IS LESS THAN NUMBER OF DATA'
         write(*,*)'INCREASE NROWMX IN read_ebsd'
         call xit()
      ENDIF
      rewind (unit = ur)

c---------------------------------------------------------------
c     3. Read ebsd_map arrays from the file.
c---------------------------------------------------------------
c     Remove the unnecessary data and save the phi1, PHI, phi2,
c     x, y data in ebsd_map
      do i = 1, no_headlines
         read(ur,'(a)') line
      end do
c     read data into ebsd_map array
      do i = 1, nlines_ebsd
         read(ur,*) (ebsd_map(j,i),j=1,5)
      end do
c     scale X and Y coordinates by multiplying f_scale factor.
      do i = 1, nlines_ebsd
         ebsd_map(4,i) = ebsd_map(4,i)*f_scale
         ebsd_map(5,i) = ebsd_map(5,i)*f_scale
      end do
      return
      END SUBROUTINE read_ebsd
