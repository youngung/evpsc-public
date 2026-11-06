c---------------------------------------------------------------
c     program <diff> used to post-analysis <stv> files resulting
c     from EVPSC-FE simulations to obtain <texture>,
c     <lattice strain> and weights
c---------------------------------------------------------------
c     This program is to use 'stv' files generated from
c     evpsc-umat for calculating diffraction patterns
c     and for downloading texture file(s) from simulated polycrystals
c     pertaining to an integration point.
c
c     While the diffraction feature is still under development,
c     the feature for downloading texture file(s) is available.
c---------------------------------------------------------------
c     List of arguments
c      1. ioption
c
c     In case iopt is 0:
c      2. fnsx
c      3. diff file
c
c     In case iopt is 1 or 2:
c      2. input file
c      3. output file
c      4. iph (optional, if not give, iph=1 is assumed)
c
c     In case iopt is 3, (n+1) arguments are required to process (n) stress files.
c      2-(n+1) : n number of input stress files to be combined
c      (n+1)   : output file
c---------------------------------------------------------------
c     (Ex 00)
c     For diffraction calculation
c        c:\> diffraction 0 <fninp> <filesx> <filediff>
c
c
c     (Ex 01)
c     A typical usage of diff for extracting texture files
c        from *.stv files generated from EVPSC-FE simulations
c
c        c:\> diffraction 1 <input-stv-file-name> <output-file> <iph>
c
c     The above will generated <output> file for phase <iph>
c
c
c     (Ex 02)
c     Generate a file including grain stresses (sg), texture (ag), weights
c     for a specific phase
c
c        c:\> diffraction 2 <input-stv-file-name> <output-file> <iph>
c
c
c     (Ex 03)
c     Combine <output-files> generated from (Ex 02).
c
c        c:\> diffraction 3 <file 1> <file 2> ... <file n> <output-file>
c
c
c
c     (Ex 04) a typical combination of works:
c	c:\>diffraction 1  d:\temp\20250717\144711\98841c72\umat\src\evpsc_incr/stvdir/0001_0001_00067.stv dum
c	c:\>diffraction 2  d:\temp\20250717\144711\98841c72\umat\src\evpsc_incr/stvdir/0001_0001_00067.stv MMM1
c	c:\>diffraction 2  d:\temp\20250717\144711\98841c72\umat\src\evpsc_incr/stvdir/0001_0001_00065.stv MMM2
c	c:\>diffraction 2  d:\temp\20250717\144711\98841c72\umat\src\evpsc_incr/stvdir/0001_0001_00064.stv MMM3
c	c:\>diffraction 3  MMM1 MMM2 MMM3 MMMALL.txt
c	c:\>diffraction 0  MMMALL.txt d:\temp\20250717\144711\98841c72\umat\src\evpsc_incr\sx_0  d:\temp\20250717\144711\98841c72\umat\src\evpsc_incr\diff_0
c
c
c     (Ex 05) Or, use the python script <abq_rst_diff.py>.
c      c:\> python -u abq_rst_diff.py --noel 1 2 3 4 --stvdir <stvidr> --fnout <results> --fnsx <fnsx> --fndiff <diff condition file>
c
c      The above will combine grain stress files of all integrations points belonging
c      to elements, 1, 2, 3, and 4. Afterwards, the lattice strain & weights are obtained
c      according to the diffraction condition specified in <fndiff>. When doing so,
c      the elastic constants and crystal symmetry specified in <fnsx> are used.
c---------------------------------------------------------------
      program diff
      implicit none
      include "../../../../src/vpsc_src/incr_evpsc_fe.dim"
      integer io,id,iphel,i,j,kkk,icrysym,isn(4),isb(4),sn(3),npoles
      real*8 sneq(3,24),sb(3),wgtall
      logical verbose,exist
      parameter(verbose=.true.)
c      parameter(verbose=.false.)
      integer nargsmx,iph,igb,igt,nmxgrc
      parameter(nargsmx=100,nmxgrc=20000)
      integer istatus(nargsmx),iopt,ngivenarg,nfiles
      character(len=325) com_buf(nargsmx),fnstv_in,fnout,
     $     fnstr(nargsmx-1),dumchar,fninp,filesx,fndiff
      real*8 sgc(6,nmxgrc),wgtc(nmxgrc),agc(3,3,nmxgrc),wgtot,aux6(6),
     $     aux33(3,3),aux66(6,6),aux3333(3,3,3,3),celv(6,6),
     $     mgelb(6,6,nmxgrc)
c     ---
c     initialize b-basis converting matrices.
      call chg_basis(aux6,aux33,aux66,aux3333,0,6)
      write(*,*)'** diffraction.exe'
      write(*,*)'** Read command line arguments'
      write(*,*)'** iargc():',iargc()
c     Determine the iopt first, and if that is 3, acts differently.
c     Otherwise, the number of arguments should be 4.

      istatus(:)=0
      ngivenarg=0
      iph=1
      do i=1,iargc()
         call getarg(i,com_buf(i))
         if (com_buf(i).ne.'')then
            istatus(i)=1
            ngivenarg=ngivenarg+1
         endif
      enddo

c------------------------------------------------------------------
      if (ngivenarg.eq.0) then
         write(*,*) "This program is to use 'stv' files generated from"
         write(*,*) "EVPSC-FE for calculating diffraction patterns"
         write(*,*) "and for downloading texture file(s) from simulated"
         write(*,*) "polycrystals pertaining to an integration point."
         write(*,*) "For detailed instruction, read the head lines of"
         write(*,*) "<src/diffraction.f> source code."
         stop 0
      endif
c------------------------------------------------------------------
c     If argument is given
      read(com_buf(1),*) iopt
      write(*,*)'** iopt:',iopt
c------------------------------------------------------------------
c     Process arguments
c     iopt 1 and 2 share common procedure
c     iopt 0 and 3 have their own procedures
      if (iopt.eq.1 .or. iopt.eq.2) then
         do i=1+1,iargc()
            if (i.eq.2) then
               read(com_buf(2),'(a300)') fnstv_in
               fnstv_in=fnstv_in(:len_trim(fnstv_in))
            elseif (i.eq.3) then
               read(com_buf(3),'(a)') fnout
            elseif (i.eq.4) then
               read(com_buf(4),'(i3)') iph
            endif
         enddo
      elseif(iopt.eq.0) then
         write(*,*)'** Ioption 0 simulation diffraction'
         write(*,*)'iargc():',iargc()
         if (iargc().ne.4) then
            write(*,*)'** Err: iopt 0 needs 4 arguments including iopt'
            stop -1
         endif

         read(com_buf(2),'(a300)') dumchar
         fninp=dumchar(:len_trim(dumchar))
         read(com_buf(3),'(a300)') dumchar
         filesx=dumchar(:len_trim(dumchar))
         read(com_buf(4),'(a300)') dumchar
         fndiff=dumchar(:len_trim(dumchar)) !! phase 1 is hardwired.

         write(*,*)'fninp:',trim(fninp)
         write(*,*)'filesx:',trim(filesx)
         write(*,*)'fndiff:',trim(fndiff)

c     obtain data crystal information

         write(*,*)'** Warning: isx_ver is not considered'
         io=10
         open(io,file=filesx,status='unknown')
         do i=1,4
            read(io,'(a)')prosa
         enddo
         do i=1,6
            read(io,*)(celv(i,j),j=1,6) !! elastic constants in crystal axes
         enddo
         close(io)

c     obtain the number of grains, stress, ag, wgt
         io=1
         open(io,file=fninp,status='unknown')
         call get_info_diff(io,wgtot,j)
         write(*,*)'wgtot:',wgtot
         write(*,*)'ngrmx:',j
         igb=1
         igt=igb+j-1
         if (ngrmx.lt.j) then
            write(*,*)'** Error: Cannot combine files:'
            write(*,*)'   j is less than ngrmx'
            write(*,*)'ngrmx: ',ngrmx,'j:',j
            stop -1
         endif
         call read_stress_ori_wgt_io(io,j,
     $        sgc(:,igb:igt),agc(:,:,igb:igt),wgtc(igb:igt))
         close(io)


c        coordinate transformation applies to bux3333 to obtain mgel
         write(*,*)'celv:',celv(1,1)
         call transform_cel_mgel(igt-igb+1,celv,agc(:,:,igb:igt),
     $        mgelb(:,:,igb:igt))

         i=1 !! filesx
         io=7 !! fndiff
         open(i,      file=trim(filesx),status='unknown')
         open(io,     file=trim(fndiff),status='unknown')

         j=igt-igb+1
         call diff_stv_str(i,io,fndiff,ndifmx,j,
     $        sgc(:,igb:igt),agc(:,:,igb:igt),wgtc(igb:igt),
     $        mgelb(:,:,igb:igt))
         close(i)
         close(io)
         stop 0

      elseif(iopt.eq.3) then
         write(*,*)'** Ioption 3 combines stress / ori / weight files'
         write(*,*)'iargc():',iargc()
         do i=1,iargc()-2
            read(com_buf(i+1),'(a300)') dumchar
            fnstr(i)=dumchar(:len_trim(dumchar))
         enddo
         read(com_buf(iargc()),'(a300)')dumchar
         fnout=dumchar(:len_trim(dumchar))
         do i=1,iargc()-2
            write(*,*)'** fnstr:',trim(fnstr(i))
         enddo
         write(*,*)'** fnout:',fnout
      endif
c------------------------------------------------------------------
      if (iopt.eq.1.or.iopt.eq.2) then
         io=1
         open(io,file=fnstv_in,access='sequential')
         call save_diffraction(io,1,NPH,ngr,ag,sg,mgel,wgt,
     $        axisph,eulerph,epsvm)
         close(io)
         io=2
         igb=ngr(iph-1)+1
         igt=ngr(iph)
         open(io,file=trim(fnout))
         write(*,*)trim(fnout)
         write(*,*)'iph:',iph
         write(*,*)'igb,igt:',igb,igt
      endif
      if (iopt.eq.1) then
         write(*,*)'** Extract texture to file'

         call write_texture_fe(io,iph,ngrmx,nphmx,igb,igt,
     $        epsvm,ag,wgt,axisph,eulerph)
         close(io)
         write(*,*)'** texture were extracted to ',trim(fnout)
         stop 0
      elseif (iopt.eq.2) then
         write(*,*)'** Extract grain stress (sg), texture(ag)'
         write(*,*)'   and weights (wgt) for a specific phase'
         call write_stress_ori_wgt_io(io,ngrmx,igb,igt,sg,ag,wgt)
         close (io)
         stop 0
      elseif (iopt.eq.3) then
         !!
         nfiles = iargc()-2
         igb=1
         do i=1,nfiles
!           open that file, determine no of grains to read,
            io=i*1
            open(io,file=trim(fnstr(i)),status='unknown')
            call get_info_diff(io,wgtot,j)
            if (ngrmx.lt.j) then
               write(*,*)'** Error: Cannot combine files:'
               write(*,*)'   j is less than ngrmx'
               write(*,*)'ngrmx: ',ngrmx,'j:',j
               stop -1
            endif
            igt=igb+j-1
            call read_stress_ori_wgt_io(io,j,
     $           sgc(:,igb:igt),agc(:,:,igb:igt),wgtc(igb:igt))
            igb=igt+1
            close(io)
         enddo

         wgtall=0.
         do i=1,igt
            wgtall=wgtall+wgtc(i)
         enddo
         do i=1,igt
            wgtc(i)=wgtc(i)/wgtall
         enddo
c     now, save the files to fnout
         io=2
         open(io,file=trim(fnout),status='unknown')
         call write_stress_ori_wgt_io(io,nmxgrc,1,igt,sgc,agc,wgtc)
         close(io)
         stop 0
      endif

      write(*,*)'** Error, diffraction capabilities of diffraction.exe'
      write(*,*)'    is not completed yet'
      call xit()
      write(*,*)' below is reserved for postmortem DIFFRACTION'
      write(*,*)' capability, which is not completed yet'

      end program diff

c     -----------------------------------------------------------------
c     Write stress / orientation / weights
c     -----------------------------------------------------------------
      subroutine write_stress_ori_wgt_io(uw,ngrmx,igb,igt,sg,ag,wgt)
      implicit none
      integer, intent(in):: uw,ngrmx,igb,igt
      real*8, intent(in):: sg(6,ngrmx),ag(3,3,ngrmx),wgt(ngrmx)
c     local
      integer kkk,i,j
      real*8 aat(3,3),gda,gdb,gdc,wgtot

c     calculate total number of grains and total vol. frac
      wgtot=0d0
      do KKK=igb,igt            !1,ngrx
         wgtot=wgtot+wgt(kkk)
      enddo
      write(uw,'(a)')        '** cols 1-6  : Grain stress in b-basis'
      write(uw,'(a)')        '** cols 7-9  : Euler angles'
      write(uw,'(a)')        '** cols 10   : Weight'
      write(uw,'(i8)')  igt-igb+1 !,'** No. grains'
      write(uw,'(e12.4)') wgtot     !,'** Total vf  '
      do KKK=igb,igt
         do i=1,3
            do j=1,3
               aat(i,j)=ag(j,i,KKK)
            enddo
         enddo
         call euler(1,gda,gdb,gdc,aat)
         write(uw,'(6e12.4,5x,4e12.4)') (sg(i,kkk),i=1,6),
     $        gda,gdb,gdc,wgt(KKK)
      enddo
      return
      end subroutine
c     -----------------------------------------------------------------
c     Read information
c     -----------------------------------------------------------------
      subroutine get_info_diff(uw,wgtot,ngr)
      implicit none
      integer, intent(in) :: uw
      real*8, intent(out) :: wgtot
      integer, intent(out) :: ngr
c     locals
      character(324) prosa
c
      read(uw,'(a)')prosa
      read(uw,'(a)')prosa
      read(uw,'(a)')prosa
      read(uw,'(i8)')ngr
      read(uw,*)wgtot
      return
      end subroutine
c     -----------------------------------------------------------------
c     Write stress / orientation / weights
c     -----------------------------------------------------------------
      subroutine read_stress_ori_wgt_io(uw,ngr,sg,ag,wgt)
      implicit none
      integer, intent(in):: uw,ngr
      real*8, intent(out):: sg(6,ngr),ag(3,3,ngr),wgt(ngr)
c     locals
      integer i,j,kkk
      real*8 aux(10),eul(3),at(3,3)

      ag(:,:,:)=0.

      do kkk=1,ngr
         read(uw,*) aux(:)
         sg(:,kkk)=aux(:6)
         eul(:)=aux(7:9)
         wgt(kkk)=aux(10)
         call euler(2,eul(1),eul(2),eul(3),at)
         do i=1,3
            do j=1,3
               ag(i,j,kkk)=at(j,i)
            enddo
         enddo
      enddo

      return
      end subroutine

c     -----------------------------------------------------------------
c     Voigt / Cartesian conversion
c     -----------------------------------------------------------------

      subroutine voigt(t1,t2,c2,c4,iopt)
      implicit none
c      real*8,intent(inout):: t1(6),t2(3,3),c2(6,6),c4(3,3,3,3)
      real*8 t1(6),t2(3,3),c2(6,6),c4(3,3,3,3)
c      integer, intent(in):: iopt
      integer iopt
      integer ijv(6,2)
      integer i,j,i1,i2,j1,j2,n,m
      data ((ijv(n,m),m=1,2),n=1,6)/1,1,2,2,3,3,2,3,1,3,1,2/


      if(iopt.eq.1) then
         t2(:,:)=0d0
         do i=1,6
            i1=ijv(i,1)
            i2=ijv(i,2)
            t2(i1,i2)=t1(i)
            t2(i2,i1)=t1(i)
         enddo
c
      elseif(iopt.eq.2) then
         t1(:)=0d0
         do i=1,6
            i1=ijv(i,1)
            i2=ijv(i,2)
            t1(i)=t2(i1,i2)
         enddo
c
      elseif (iopt.eq.3) then
         c4(:,:,:,:)=0d0
         do i=1,6
            i1=ijv(i,1)
            i2=ijv(i,2)
         do j=1,6
            j1=ijv(j,1)
            j2=ijv(j,2)
            c4(i1,i2,j1,j2)=c2(i,j)
            c4(i2,i1,j1,j2)=c2(i,j)
            c4(i1,i2,j2,j1)=c2(i,j)
            c4(i2,i1,j2,j1)=c2(i,j)
         enddo
         enddo
c
      elseif(iopt.eq.4) then
         c2(:,:)=0d0
         do i=1,6
            i1=ijv(i,1)
            i2=ijv(i,2)
         do j=1,6
            j1=ijv(j,1)
            j2=ijv(j,2)
            c2(i,j)=c4(i1,i2,j1,j2)
         enddo
         enddo
      else
         write(*,*)'** Error: unexpected iopt value found in voigt'
         call xit()
      endif
c
      return
      end subroutine voigt

c     -----------------------------------------------------------------
c     Apply coordinate transformation to crystal elastic constants
c       given in voigt by using <ag>, which transforms from crystal
c       to sample coordinate system. The resulting elastic constants
c       are given in the b-basis convention.
c     -----------------------------------------------------------------
      subroutine transform_cel_mgel(ngr,celv,ag,mgelb)
      implicit none

      integer, intent(in):: ngr
      real*8, intent(in):: celv(6,6),ag(3,3,ngr)
      real*8, intent(out):: mgelb(6,6,ngr)
c     locals
      real*8 celb(6,6),melb(6,6),mgelc(3,3,3,3),aux6(6),aux33(3,3),
     $     cel3333(3,3,3,3),aux3333(3,3,3,3)
      integer inverr,i,j,k,l,i0,j0,k0,l0,kkk

c     convert from voigt to b-basis
      call voigt(aux6,aux33,celv,cel3333,3)
      call chg_basis(aux6,aux33,melb,cel3333,4,6)

c     inverse
      call lu_inverse(inverr,melb,6)
      call chg_basis(aux6,aux33,melb,mgelc,3,6)

c     transformation from crystal axes to sample axes
      mgelb(:,:,1:ngr)=0d0
      do kkk=1,ngr

         aux3333(:,:,:,:)=0.d0

      do i0=1,3
      do j0=1,3
      do k0=1,3
      do l0=1,3
         do i=1,3
         do j=1,3
         do k=1,3
         do l=1,3
            aux3333(i0,j0,k0,l0)=aux3333(i0,j0,k0,l0)+mgelc(i,j,k,l)
     $           * ag(i0,i,kkk) * ag(j0,j,kkk)
     $           * ag(k0,k,kkk) * ag(l0,l,kkk)
         enddo
         enddo
         enddo
         enddo
      enddo
      enddo
      enddo
      enddo

      call chg_basis(aux6,aux33,melb,aux3333,4,6)

      mgelb(:,:,kkk)=melb(:,:)

      enddo                     ! end of loop over grains

      end subroutine

      include "../../src/vpsc_src/library7.sub"
