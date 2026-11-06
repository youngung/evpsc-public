      subroutine stop_debug(iunit,chr)
      implicit none
      integer iunit
      character(len=*) chr
      call fill_line(iunit,'*',52)
      call w_chr(iunit,chr)
      call w_chr(iunit, '--  Reached safely at the debugging point  --')
      call fill_line(iunit,'*',52)
      call exit(-1)
      end subroutine stop_debug
c-----------------------------------------------------------------------
      subroutine w_val(iunit,str,v)
c     Arguments
c     iunit: file ID
c     str : chracter that describe the passed value
c     v   : the value
      implicit none
      integer iunit,nc
      character*80 fmt,get_fmt,ncc
      character(len=*) str
      real*8 v
      fmt = get_fmt(v)
      nc = len(str)
      write(ncc,'(i80)') nc
      ncc=adjustl(ncc)
      fmt = trim('(a')//trim(ncc)//',x,'//trim(fmt)//')'
      if (iunit.eq.0) then
         write(*,    trim(fmt)) trim(str),v
      else
         write(iunit,trim(fmt)) trim(str),v
      endif
      return
      end subroutine w_val
c-----------------------------------------------------------------------
      subroutine w_vdim(iunit,array,ndi,fact)
c     Arguments
c     ---------
c     iunit: file ID
c     array: 1-Dimensional array
c     ndi: size of the array
c     fact:multiplicative factor to scale the elements in the array
      implicit none
      integer, intent(in) :: iunit, ndi
      dimension array(ndi)
      real*8, intent(in) ::  array
      dimension aux33(3,3),auxbb(ndi,ndi),aux3333(3,3,3,3),aux6(6)
      real*8 aux33,auxbb,aux3333,aux6
      real*8, intent(in) :: fact
      call chg_basis(array,aux33,auxbb,aux3333,1,ndi)
      call voigt(aux6,aux33,auxbb,aux3333,2)
      call w_dim(iunit,aux6,6,fact,.true.)
      return
      end subroutine w_vdim
c-----------------------------------------------------------------------
      subroutine w_dim(iunit,array,ndi,fact,ibr)
c     iunit: file ID
c     array: 1-Dimensional array
c     ndi: size of the array
c     fact:multiplicative factor to scale the elements in the array
c     ibr: flag to insert line-breaker
      implicit none
      integer, intent(in) :: iunit, ndi
      dimension array(ndi)
      real*8, intent(in) ::  array
      logical, intent(in) :: ibr
c     --
      dimension brray(ndi)
      real*8 brray
      character*80 fmt,clen,get_fmt
      integer i
      real*8 fact,get_mx_bk,mxv
c     --
c     color coding for the terminal
      logical color
c      parameter(color=.false.)
      parameter(color=.true.)
      character(len=20) head,tail

      if (color) then
         head=trim(achar(27)//'[37;41;234m')
         tail=trim(achar(27)//'[0m')
      else
         head=''
         tail=''
      endif

      do i=1,ndi
         brray(i) = array(i) * fact
      enddo
      mxv = get_mx_bk(brray,ndi)
      fmt = get_fmt(mxv)

      if (iunit.ne.0) then
         fmt = '(a,'//trim(fmt)//',a)'
         write(iunit,fmt,advance='no') head,(brray(i),i=1,ndi),tail
         if (ibr) write(iunit,*)
      else
         write(clen,"(i2)") ndi
         fmt = '(a,'//trim(clen)//trim(fmt)//',a)'
         if (ibr) then
            write(*   ,fmt) head,(brray(i),i=1,ndi),tail
         else
            write(*   ,fmt,advance='no') head,(brray(i),i=1,ndi),tail
         endif
      endif

      return
      end subroutine w_dim
c-----------------------------------------------------------------------
      subroutine w_chrc(iunit,str)
c     iunit: file ID
c     chr : chracter
      implicit none
      integer, intent(in) :: iunit
      character(len=*), intent(in) ::  str
      character*80 fmt,get_fmt_str

c     color coding for the terminal
      logical color
      parameter(color=.true.)
      character(len=20) head,tail

      if (color) then
         head=trim(achar(27)//'[1;43;234m')
         tail=trim(achar(27)//'[0m')
      else
         head=''
         tail=''
      endif

      fmt = get_fmt_str(head//str//tail)
      if (iunit.eq.0) then
         write(*,    trim(fmt),advance='no') head//trim(str)//tail
         write(*,*)
      else
         write(iunit,trim(fmt),advance='no') head//trim(str)//tail
      endif
      return
      end subroutine w_chrc
c-----------------------------------------------------------------------
      subroutine w_ival(iunit,str,ivv)
c     Arguments
c     iunit: file ID
c     str : chracter that describe the passed value
c     ivv : the integer value
      implicit none
      integer, intent(in):: ivv, iunit
      character(len=*), intent(in) ::  str
c     --
      character*80 fmt,get_fmt_int,ncc,fmt_tmp
      integer nc
c     color coding for the terminal
      logical color
      parameter(color=.true.)
      character(len=20) head,tail
      character(len=80) get_fmt_str


      if (ivv.le.0) then
         write(*,*)'ivv less than or equal to zero (w_ival)'
         call xit()
      endif

      if (color) then
         head=trim(achar(27)//'[4;231m')
         tail=trim(achar(27)//'[0m')
      else
         head=''
         tail=''
      endif

      fmt = get_fmt_int(ivv)
      nc = len(str)
      write(ncc,'(i10)') nc
      ncc=adjustl(ncc)
      fmt = trim('(a')//trim(ncc)//',x,'//trim(fmt)//')'
      if (iunit.eq.0) then
         fmt_tmp=get_fmt_str(head)
         write(*,fmt_tmp,advance='no')head
         write(*,    trim(fmt),advance='no') trim(str),ivv
         fmt_tmp=get_fmt_str(tail)
         write(*,fmt_tmp,advance='no')tail
         write(*,*)
      else
         write(iunit,trim(fmt)) trim(str),ivv
      endif
      return
      end subroutine w_ival
c-----------------------------------------------------------------------
c     print out a 2nd order square matrix array
      subroutine w_mdim(iunit,array,ndi,fact)
c     iunit: file ID
c     array: 2-Dimensional array
c     ndi: size of the (ndi x ndi) array
c     fact: multiplicative factor to scale the elements in the array
      implicit none
      character*80 fmt,clen,get_fmt
      integer, intent(in):: iunit,ndi
      real*8, intent(in):: array(ndi,ndi),fact
c     locals
      real*8 brray(ndi,ndi),get_mmx,mxv
      integer i,j
      brray(:,:)=0d0
      do i=1,ndi
      do j=1,ndi
         brray(i,j) = array(i,j) * fact
      enddo
      enddo
      mxv = get_mmx(brray,ndi,ndi)
      fmt = get_fmt(mxv) ! common filter

      write(clen,"(i2)") ndi
      fmt = '('//trim(clen)//trim(fmt)//')'

      if (iunit.ne.0) then
         do i=1,ndi
            write(iunit,fmt) (brray(i,j),j=1,ndi)
         enddo
      else
         do i=1,ndi
            write(*    ,fmt) (brray(i,j),j=1,ndi)
         enddo
      endif
      return
      end subroutine w_mdim
c-----------------------------------------------------------------------
      real*8 function get_mx_bk(array,ndi)
      implicit none
      integer, intent(in):: ndi
      dimension array(ndi)
      real*8, intent(in):: array

      integer i

      get_mx_bk=array(1)
      do i=2,ndi
         if (array(i).gt.get_mx_bk) then
            get_mx_bk=array(i)
         endif
      enddo
      return
      end function get_mx_bk
c-----------------------------------------------------------------------
      subroutine w_chr(iunit,str)
c     iunit: file ID
c     chr : chracter
      implicit none
      integer,intent(in):: iunit
      character(len=*), intent(in):: str
      character*80 fmt,get_fmt_str

c     --
c     color coding for the terminal
      logical color
      parameter(color=.true.)
      character(len=20) head,tail

      if (color) then
         head=trim(achar(27)//'[0;42;234m')
         tail=trim(achar(27)//'[0m')
      else
         head=''
         tail=''
      endif

      fmt = get_fmt_str(head//str//tail)
      if (iunit.eq.0) then
         write(*,    trim(fmt)) head//trim(str)//tail
      else
         write(iunit,trim(fmt)) head//trim(str)//tail
      endif
      return
      end subroutine w_chr
c-----------------------------------------------------------------------w
c     Fill a line(s) with the given <chr> by repeating it for <n> times
      subroutine fill_line(imsg,chr,n)
c     imsg: file unit - if imsg =0 use std (*)
c     chr : symbol that will fill the lines
c     n   : the number of repetition
      implicit none
      character*1, intent(in) :: chr
      integer, intent(in):: imsg,n
      integer i
c     --
      if (imsg.eq.0) then
         do i=1,n
            write(*,   '(a1)',advance='no') chr
         enddo
         write(*,*)
      else
         do i=1,n
            write(imsg,'(a1)',advance='no') chr
         enddo
         write(imsg,*)
      endif

      return
      end subroutine fill_line
c-----------------------------------------------------------------------
c     Given the character value <str> find a suitible format to write
c     this value
      character*80 function get_fmt_str(str)
c     Argument
c     str: the subjected character with any arbitrary size
      implicit none
      character(len=*) str
      character*80 ncc
      integer nc
      nc = len(str)
      write(ncc,'(i80)') nc
      ncc = adjustl(ncc)
      get_fmt_str='(a'//trim(ncc)//',x)'
      return
      end function get_fmt_str
c$$$c-----------------------------------------------------------------------
c     Given the integer value <ival> find a suitible format to write
c     this value.
      character*80 function get_fmt_int(ival)
c     Argument
c     ival: integer value
      implicit none
      integer, intent(in) :: ival
      integer nc
      character*80 ncc
      character*20 nc_no
      write(ncc,'(i10)')ival
      ncc = adjustl(ncc)
      nc=len(trim(ncc)) ! number of character for this integer value
      write(nc_no,'(i10)') nc
      get_fmt_int='(i'//trim(adjustl(nc_no))//',x)'
      return
      end function get_fmt_int
c-----------------------------------------------------------------------
      real*8 function get_mmx(array,ncol,nrow)
      implicit none
      integer, intent(in) :: ncol,nrow
      dimension array(ncol,nrow)
      real*8, intent(in) ::  array
      integer i,j
      get_mmx=array(1,1)
      do i=1,ncol
      do j=1,nrow
         if (array(i,j).gt.get_mmx) then
            get_mmx=array(i,j)
         endif
      enddo
      enddo
      return
      end function get_mmx
c-----------------------------------------------------------------------
c     Get the suitable format to print the given value <val>
      character*80 function get_fmt(val)
c     Argument
c     val: the value
      implicit none
      real*8 val,aval
c      logical is_inf!, isnan

c      if (isnan(val)) then
c         write(*,*)'Warning: Given value is a NaN'
c         call xit()
c      endif

      aval = dabs(val)

c      if (isnan(aval)) then
c         write(*,*)'Warning: Given value is a NaN'
c         call xit()
c      endif

c      if     (is_inf(aval))                      then

      if (.false.) then
         get_fmt='f10.0,x'
      elseif (                   (aval.lt.1e-3)) then
         get_fmt="e12.4,x"
      elseif ((aval.ge.1e-3).and.(aval.lt.1e0 )) then
         get_fmt="f9.5,x"
      elseif ((aval.ge.1e0 ).and.(aval.lt.1e3 )) then
         get_fmt="f9.3,x"
      elseif ((aval.ge.1e3 ).and.(aval.lt.1e6 )) then
         get_fmt="f9.0,x"
      elseif ((aval.ge.1e6 ))                    then
         get_fmt="e12.4,x"
      else
         write(*,*) 'Err: Unexpected case of max value'
         write(*,*) 'was given to get_fmt'
         write(*,*) 'given value:', val
c         call exit(-1)
         call xit()
      endif

      get_fmt="e12.4,x" ! override

      return
      end function get_fmt
c-----------------------------------------------------------------------
