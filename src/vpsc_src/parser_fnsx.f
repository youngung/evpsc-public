c-----------------------------------------------------------------------
c     program to parse single crystal file.
c     This program is primarily used to obtain 'nmodmx','ntwmmx','nsysmx',
c     and 'ntwsmx' from a given single crystal file.
c
c     The executable can be obtained by
c     make parsers
c     in the root directory of this repository.
c
c     ./sx_parser <file-name-of-sx>
c
c     Example run:
c     ./sx_parser matData/vpscData/AA_jinhwa/sx_Al6061_affine.sx
c     nmodmx:           1
c     ntwmmx:           0
c     nsysmx:          12
c     ntwsmx:           0
c     ihardlaw:         ?
c     isx_ver:          0
c-----------------------------------------------------------------------
      program parser_fnsx
      implicit none
      character*128 com_buf(1),prosa
      integer i,isx_ver,UR_SX
      parameter(UR_SX=1)
      integer nmodmx,ntwmmx,nsysmx,ntwsmx,ihardlaw

      if (iargc().eq.0) then
         write(*,*)'**Error: No argument is given to FNSX <parser>'
         stop -1
      endif

      do i=1,iargc()
         call get_command_argument(i,com_buf(i))
      enddo
      open(UR_SX,file=trim(com_buf(1)))
c      call sx_version_tester(UR_SX,isx_ver)
      call sx_version_tester_nsch(UR_SX,isx_ver)
      rewind(UR_SX)
c     Find nmodmx,ntwmmx,nsysmx,ntwsmx
      call get_dimension(UR_SX,isx_ver,nmodmx,ntwmmx,nsysmx,ntwsmx)
      do i=1,3
         read(UR_SX,*)prosa
      enddo
      read(UR_SX,*) ihardlaw
      close(UR_SX)

 10   format(a10,x,i2)
      write(*,10)'nmodmx:',nmodmx
      write(*,10)'ntwmmx:',ntwmmx
      write(*,10)'nsysmx:',nsysmx
      write(*,10)'ntwsmx:',ntwsmx
      write(*,10)'ihardlaw:',ihardlaw
      write(*,10)'isx_ver:',isx_ver

      end program parser_fnsx
c-----------------------------------------------------------------------
      subroutine get_dimension(UR_SX,isx_ver,
     $     nmodmx,ntwmmx,nsysmx,ntwsmx)
      implicit none
      integer, intent(in):: UR_SX,isx_ver
      integer, intent(out)::nmodmx,ntwmmx,nsysmx,ntwsmx
c     locals
      integer nskip_lines,i,nmodesx,nmodes,imod
      integer mode(10)
      integer modex,nsmx,isensex,itwtypex,mloop,mcount,is,isys
      real*8 twshx,c1,c2,c3,c4,c5
      character*128 prosa
      mode(:)=-1

      if (isx_ver.eq.0) then
         nskip_lines=13
      elseif(isx_ver.eq.1) then
         nskip_lines=15
      elseif(isx_ver.eq.2) then
         nskip_lines=15
      elseif(isx_ver.eq.3) then
         nskip_lines=13
      elseif(isx_ver.eq.4) then
         nskip_lines=16
      elseif(isx_ver.eq.5) then
         nskip_lines=16
      else
         write(*,*)'**Error: unexpected isx_ver given'
         stop -1
      endif
      do i=1,nskip_lines
         read(UR_SX,'(a)')prosa
      enddo

      if (isx_ver.eq.3.or.isx_ver.eq.4.or.isx_ver.eq.5) then
         read(UR_SX,*)c1,c2,c3,c4,c5
      endif

      write(*,*)'isx_ver:',isx_ver
      write(*,*)'nskip_lines:',nskip_lines
      write(*,*)'*** prosa:',prosa
c      read(UR_SX,*)prosa
c      write(*,*)'**beging ',trim(prosa)
c      stop -1
      read(UR_SX,*)nmodesx
      read(UR_SX,*)nmodes
      read(UR_SX,*)(mode(i),i=1,nmodes)

      mcount=1
      nsysmx=0
      ntwmmx=0
      ntwsmx=0
      do 100 mloop=1,nmodesx
         read(UR_SX,'(a)')prosa
         read(UR_SX,*)modex,nsmx,isensex,itwtypex
         if (itwtypex.ne.0) then
            read(UR_SX,*)twshx
         endif
         if (mloop.ne.mode(mcount)) then
            do is=1,nsmx
               read(UR_SX,*)
            enddo
            goto 100
         endif
         nsysmx=nsysmx+nsmx
         if (itwtypex.eq.0) then
         elseif(itwtypex.ne.0) then
            ntwmmx=ntwmmx+1
            ntwsmx=ntwsmx+nsmx
         endif
         mcount=mcount+1
         do isys=1,nsmx
            read(UR_SX,*)
         enddo
 100  continue
      nmodmx=nmodes
      return
      end subroutine
