c-----------------------------------------------------------------------
c     write an empty line(s)
      subroutine w_empty_lines(imsg,n)
c     Arguments
c     imsg: file unit - if imsg =0 use std (*)
c     n   : the number of lines to be emptied
      implicit none
      integer imsg,n,i
      if (imsg.eq.0) then
         do i=1,n
            write(*,*)
         enddo
      else
         do i=1,n
            write(imsg,*)
         enddo
      endif
      return
      end subroutine w_empty_lines
c-----------------------------------------------------------------------
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
c     Get the suitable format to print the given value <val>
      character*80 function get_fmt(val)
c     Argument
c     val: the value
      implicit none
      real*8 val,aval
      logical is_inf, isnan

      if (isnan(val)) then
         write(*,*)'Warning: Given value is a NaN'
         call xit()
      endif

      aval = abs(val)

      if (isnan(aval)) then
         write(*,*)'Warning: Given value is a NaN'
         call xit()
      endif

      if     (is_inf(aval))                      then
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
         call exit(-1)
      endif
      return
      end function get_fmt
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
c-----------------------------------------------------------------------
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
      subroutine w_vals(iunit,v)
c     Arguments
c     iunit: file ID
c     v   : the value
      implicit none
      integer iunit
      character*80 fmt,get_fmt
      real*8 v
      fmt = '('//trim(get_fmt(v))//')'
      if (iunit.eq.0) then
         write(*,    trim(fmt)) v
      else
         write(iunit,trim(fmt)) v
      endif
      return
      end subroutine w_vals
c-----------------------------------------------------------------------
      subroutine w_valsc(iunit,v)
c     Arguments
c     iunit: file ID
c     v   : the value
      implicit none
      integer iunit
      character*80 fmt,get_fmt
      real*8 v
      fmt = '('//trim(get_fmt(v))//')'
      if (iunit.eq.0) then
         write(*,    trim(fmt),advance='no') v
      else
         write(iunit,trim(fmt),advance='no') v
      endif
      return
      end subroutine w_valsc
c-----------------------------------------------------------------------
      subroutine w_valc(iunit,str,v)
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
         write(*,    trim(fmt),advance='no') trim(str),v
      else
         write(iunit,trim(fmt),advance='no') trim(str),v
      endif
      return
      end subroutine w_valc
c-----------------------------------------------------------------------
      subroutine w_ival(iunit,str,iv)
c     Arguments
c     iunit: file ID
c     str : chracter that describe the passed value
c     iv  : the integer value
      implicit none
      integer, intent(in):: iv, iunit
      character(len=*), intent(in) ::  str
c     --
      character*80 fmt,get_fmt_int,ncc,fmt_tmp
      integer nc
c     color coding for the terminal
      logical color
      parameter(color=.true.)
      character(len=20) head,tail
      character(len=80) get_fmt_str

      if (color) then
         head=trim(achar(27)//'[4;231m')
         tail=trim(achar(27)//'[0m')
      else
         head=''
         tail=''
      endif

      fmt = get_fmt_int(iv)
      nc = len(str)
      write(ncc,'(i10)') nc
      ncc=adjustl(ncc)
      fmt = trim('(a')//trim(ncc)//',x,'//trim(fmt)//')'
      if (iunit.eq.0) then
         fmt_tmp=get_fmt_str(head)
         write(*,fmt_tmp,advance='no')head
         write(*,    trim(fmt),advance='no') trim(str),iv
         fmt_tmp=get_fmt_str(tail)
         write(*,fmt_tmp,advance='no')tail
         write(*,*)
      else
         write(iunit,trim(fmt)) trim(str),iv
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
      fmt = get_fmt(mxv) ! filter

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
c     print out a 2nd order non-square matrix array
      subroutine w_mndim(iunit,array,ndi,nci,fact)
c     iunit: file ID
c     array: 2-Dimensional array
c     ndi: size of the (ndi x ndi) array
c     fact: multiplicative factor to scale the elements in the array
      implicit none
      character*80 fmt,clen,get_fmt
      integer ndi,nci,iunit,i,j
      dimension array(ndi,nci),brray(ndi,nci)
      real*8 array,brray,get_mmx,mxv,fact
      do i=1,ndi
      do j=1,nci
         brray(i,j) = array(i,j) * fact
      enddo
      enddo
      mxv = get_mmx(brray,ndi,nci)
      fmt = get_fmt(mxv) ! filter

      write(clen,"(i2)") nci
      fmt = '('//trim(clen)//trim(fmt)//')'

      if (iunit.ne.0) then
         do i=1,ndi
            write(iunit,fmt) (brray(i,j),j=1,nci)
         enddo
      else
         do i=1,ndi
            write(*   ,fmt) (brray(i,j),j=1,nci)
         enddo
      endif
      return
      end subroutine w_mndim
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
      real*8 fact,get_mx,mxv
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
      mxv = get_mx(brray,ndi)
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
c     print header
      subroutine print_head(i)
      implicit none
c     i: file unit (use std if 0, use imsg elsewhere)
      integer i
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

      if (i.eq.0) then ! USE std
         call w_empty_lines(0,2)
         write(*,*) head//'*------------------*'//tail
         write(*,*) head//'|       UMAT       |'//tail
         write(*,*) head//'*------------------*'//tail
         call w_empty_lines(0,1)
      else
         call w_empty_lines(i,2)
         write(i,*) head//'*------------------*'//tail
         write(i,*) head//'|       UMAT       |'//tail
         write(i,*) head//'*------------------*'//tail
         call w_empty_lines(i,1)
      endif

      return
      end subroutine print_head
c-----------------------------------------------------------------------
c     print footer
      subroutine print_foot(i)
      implicit none
c     i: file unit (use std if 0, use imsg elsewhere)
      integer i

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

      if (i.eq.0) then ! USE std
         call w_empty_lines(0,2)
         write(*,*)head//'*------------------*'//tail
         write(*,*)head//'|       END        |'//tail
         write(*,*)head//'*------------------*'//tail
         call w_empty_lines(0,1)
      else
         call w_empty_lines(i,2)
         write(i,*)head//'*------------------*'//tail
         write(i,*)head//'|       END        |'//tail
         write(i,*)head//'*------------------*'//tail
         call w_empty_lines(i,1)
      endif

      return
      end subroutine print_foot
c-----------------------------------------------------------------------
      subroutine stop_debug(iunit,chr)
      implicit none
      integer iunit
      character(len=*) chr
      call fill_line(iunit,'*',52)
      call w_chr(iunit,chr)
      call w_chr(iunit, '--  Reached safely at the debugging point  --')
      call fill_line(iunit,'*',52)
      call exit(-1)
      end subroutine
c-----------------------------------------------------------------------
      real*8 function get_mx(array,ndi)
      implicit none
      dimension array(ndi)
      real*8, intent(in):: array
      integer, intent(in):: ndi
      integer i

      get_mx=array(1)
      do i=2,ndi
         if (array(i).gt.get_mx) then
            get_mx=array(i)
         endif
      enddo
      return
      end function get_mx
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
      subroutine pjoin(pardir,fn,fullfn)
      implicit none
      character(len=255) :: path,fn,fullfn,pardir
      character(len=1) :: sep
      call get_environment_variable('PATH',path)
      sep=path(1:1)
      fullfn = trim(pardir) // sep // trim(fn)
      end subroutine
c-----------------------------------------------------------------------
!     various testing subroutines/functions
c-----------------------------------------------------------------------
      logical function isnan_local(val)
c     Argument
c     --------
c     val
      implicit none
      real*8, intent(in)::val
      isnan_local=.false.
      if (val.ne.val) isnan_local=.true.
      return
      end function isnan_local
c-----------------------------------------------------------------------
      logical function is_inf(val)
c     check if the given value is 'infinite'
c     val: value subjected to test
      use, intrinsic :: ieee_arithmetic
      implicit none
      real*8 val
      is_inf = .not.(ieee_is_finite(val))
      return
      end function is_inf
c-----------------------------------------------------------------------
      logical function is_finite(val)
c     check if the given value is 'finite'
c     val: value subjected to test
      use, intrinsic :: ieee_arithmetic
      implicit none
      real*8, intent(in):: val
      is_finite = ieee_is_finite(val)
      return
      end function is_finite
c-----------------------------------------------------------------------
      logical function isnan_in_arr(array,ndim)
c     check if there is any element that is NaN in the given array.

c     Arguments
c     ---------
c     array
c     ndim
      implicit none
      integer i,ndim
      real*8 array(ndim)
      isnan_in_arr=.false.
      do i=1,ndim
         if (isnan(array(i))) then
            isnan_in_arr=.true.
            goto 10
         endif
      enddo
 10   return
      end function isnan_in_arr
c-----------------------------------------------------------------------
      logical function isnan_in_marr(array,ncol,nrow)
c     Check if there is any element that is NaN in the
c     given multiaxial array.

c     Arguments
c     ---------
c     array
c     ncol
c     nrow
      implicit none
      integer i,j,ncol,nrow
      real*8 array(ncol,nrow)
      isnan_in_marr=.false.
      do i=1,ncol
      do j=1,nrow
         if (isnan(array(i,j))) then
            isnan_in_marr=.true.
            goto 10
         endif
      enddo
      enddo
 10   return
      end function isnan_in_marr
c-----------------------------------------------------------------------
      logical function is_all_element_finite(array,narray,verbose)
c     Arguments                                                  intents
c     ---------                                                  -------
c     array                                                       <in>
c     narray                                                      <in>
c     verbose                                                     <in>
      implicit none
      integer, intent(in)::narray
      dimension array(narray)
      real*8, intent(in)::array
      logical, intent(in)::verbose
c     locals
      integer i,idia
      logical idiaw
      logical isnan_local,isnan
      character(len=10)int_char1,int_char2
      parameter(idia=0,idiaw=.true.)
c**   f2py
cf2py intent(in) array,narray,verbose
cf2py depend(narray) array
c**   f2py
      is_all_element_finite=.true.

      if (idiaw) then
         call fill_line(idia,'-',72)
         call w_chr(idia,'** Entering function is_all_element_finite')
         call int2chr(narray,int_char2)
      endif

      do i=1,narray
         if (idiaw) then
            call int2chr(i,int_char1)
            write(*,*) trim(int_char1)//'th element out of '//int_char2
         endif

         if (isnan_local(array(i)).or.isnan(array(i))) then
            is_all_element_finite=.false.
            if (verbose) call w_chr(0,'** Found a NaN value in'//
     $           ' the given array, the '//char(48+i)//'th element')
            goto 10
         endif
      enddo

      if (idiaw) then
         call w_chr(idia,'** Exiting function is_all_element_finite')
         call fill_line(idia,'-',72)
      endif
 10   return
      end function is_all_element_finite
c-----------------------------------------------------------------------
      logical function is_all_element_zero(array,narray,verbose)
c     Arguments                                                  intents
c     ---------                                                  -------
c     array                                                       <in>
c     narray                                                      <in>
c     verbose                                                     <in>
      implicit none
      integer, intent(in)::narray
      dimension array(narray)
      real*8, intent(in)::array
      logical, intent(in)::verbose
c     locals
      integer i,idia
      logical idiaw
      !logical is_zero,isnan_local,isnan
      character(len=10)int_char1,int_char2
      parameter(idia=0,idiaw=.true.)
c**   f2py
cf2py intent(in) array,narray,verbose
cf2py depend(narray) array
c**   f2py
      is_all_element_zero=.true.

      if (idiaw) then
         call fill_line(idia,'-',72)
         call w_chr(idia,'** Entering function is_all_element_zero')
         call int2chr(narray,int_char2)
      endif

      ! if (.false.) write(*,*) is_finite, verbose

      do i=1,narray
         if (idiaw) then
            call int2chr(i,int_char1)
            write(*,*) trim(int_char1)//'th element out of '//int_char2
         endif
         if (array(i).eq.0) then
            is_all_element_zero=.false.
            if (verbose) call w_chr(0,'** Found a NaN value in'//
     $           ' the given array, the '//char(48+i)//'th element')
            goto 10
         endif
      enddo

 10   if (idiaw) then
         call w_chr(idia,'** Exiting function is_all_element_finite')
         call fill_line(idia,'-',72)
      endif

      return
      end function is_all_element_zero

c-----------------------------------------------------------------------
c     various converting routines
c-----------------------------------------------------------------------
c     33 bases "stress" tensor to 6 bases stress tensor
c-----------------------------------------------------------------------
      subroutine tens33_trans(a33,b33)
      implicit none
      dimension a33(3,3),b33(3,3)
      real*8 a33,b33
      integer i,j
      b33(:,:) = 0d0
      do i=1,3
      do j=1,3
         b33(i,j) = a33(j,i)
      enddo
      enddo
      end subroutine
c-----------------------------------------------------------------------
c     33 bases "stress" tensor to 6 bases "stress" tensor
c-----------------------------------------------------------------------
      subroutine voigt1(a33,a6)
      implicit none
      dimension a33(3,3),a6(6)
      real*8 a33, a6
      integer i
cf2py intent(in) a33
cf2py intent(out) a6
      a6(:) = 0.d0
      do i=1,3
         a6(i) = a33(i,i)
      enddo
      a6(4) = a33(2,3)
      a6(5) = a33(1,3)
      a6(6) = a33(1,2)
      return
      end subroutine voigt1
c-----------------------------------------------------------------------
c     6 bases "stress" tensor to 33 bases stress tensor
c-----------------------------------------------------------------------
      subroutine voigt2(a6,a33)
      implicit none
      dimension a33(3,3),a6(6)
      real*8 a33, a6
      integer i,j,k,ijv(6,2),n,m
cf2py intent(in)  a6
cf2py intent(out) a33
      data ((ijv(n,m),m=1,2),n=1,6)/1,1,2,2,3,3,2,3,1,3,1,2/
      a33(:,:) = 0.d0
      do k=1,6
         i=ijv(k,1)
         j=ijv(k,2)
         a33(i,j) = a6(k)
         a33(j,i) = a6(k)
      enddo
      return
      end subroutine voigt2
c-----------------------------------------------------------------------
c     33 bases "strain" tensor to 6 bases strain tensor
c-----------------------------------------------------------------------
c     !! Caution, strain tensor is converted such that
c     the shear strains are multiplied by 2 to follow Abaqus convention
      subroutine voigt3(a33,a6)
      implicit none
      dimension a33(3,3),a6(6)
      real*8, intent(in) :: a33
      real*8, intent(out):: a6
      integer i
cf2py intent(in)  a33
cf2py intent(out) a6
      a6(:) = 0.d0
      do i=1,3
         a6(i) = a33(i,i)
      enddo
      a6(4) = a33(2,3)*2.
      a6(5) = a33(1,3)*2.
      a6(6) = a33(1,2)*2.
      return
      end subroutine voigt3
c-----------------------------------------------------------------------
c     6 bases "strain" tensor to 33 bases strain tensor
c-----------------------------------------------------------------------
c     !! Caution, strain tensor is converted such that
c     the shear strains are multiplied by 1/2. to follow Abaqus convention
      subroutine voigt4(a6,a33)
      implicit none
      dimension a33(3,3), a6(6)
      real*8 a33, a6, fact
      integer i,j,k,ijv(6,2),n,m
cf2py intent(in)  a6
cf2py intent(out) a33
      data ((ijv(n,m),m=1,2),n=1,6)/1,1,2,2,3,3,2,3,1,3,1,2/
      a33(:,:) = 0.d0
      do k=1,6
         i=ijv(k,1)
         j=ijv(k,2)
         if (i.ne.j) fact=0.5d0
         if (i.eq.j) fact=1.0d0
         a33(i,j) = a6(k)*fact
         a33(j,i) = a6(k)*fact
      enddo
      return
      end subroutine voigt4
c-----------------------------------------------------------------------
      subroutine reduce_6to3(a6,a3)
      dimension a6(6),a3(3)
      real*8 a6,a3
cf2py intent(in) a6
cf2py intent(out) a3
      a3(:) = 0d0
      a3(1) = a6(1)
      a3(2) = a6(2)
      a3(3) = a6(6)
      return
      end subroutine reduce_6to3
c-----------------------------------------------------------------------
      subroutine cnv_3to6_dev(a,b)
      dimension a(3),b(6)
      real*8, intent(in) ::  a
      real*8, intent(out)::  b
cf2py intent(in) a
cf2py intent(out) b
      b(:)=0d0
      b(1:2) = a(1:2)
      b(3)   = -b(1)-b(2)
      b(6)   = a(3)
      end subroutine cnv_3to6_dev
c-----------------------------------------------------------------------
      subroutine reduce_3to6(a3,a6)
      dimension a6(6),a3(3)
      real*8, intent(out):: a6
      real*8, intent(in) :: a3
cf2py intent(in) a3
cf2py intent(out) a6
      a6(:) = 0d0
      a6(1) = a3(1)
      a6(2) = a3(2)
      a6(6) = a3(3)
      return
      end subroutine reduce_3to6
c-----------------------------------------------------------------------
      subroutine reduce_33to66(a33,a66)
      dimension a66(6,6),a33(3,3)
      real*8 a66,a33
cf2py intent(in) a33
cf2py intent(out) a66
      a66(:  ,:  )=0d0
      a66(1:2,1:2)=a33(1:2,1:2)
      a66(1:2,6  )=a33(1:2,3  )
      a66(  6,1:2)=a33(  3,1:2)
      a66(  6,6  )=a33(  3,3  )
      return
      end subroutine reduce_33to66
c-----------------------------------------------------------------------
      subroutine reduce_66to33(a66,a33)
      dimension a66(6,6),a33(3,3)
      real*8 a66,a33
cf2py intent(in) a33
cf2py intent(out) a66
      a33(:  ,:  )=0d0
      a33(1:2,1:2)=a66(1:2,1:2)
      a33(1:2,3  )=a66(1:2,6  )
      a33(  3,1:2)=a66(  6,1:2)
      a33(  3,3  )=a66(  6,6  )
      return
      end subroutine reduce_66to33
c-----------------------------------------------------------------------
c     Convert 4 dimensional plane-stress tensor
c     - convention:
c      vec  ij component
c        1:      11
c        2:      22
c        3:      33
c        4:      12
      subroutine reduce_6to4(a6,a4)
      dimension a6(6),a4(4)
      real*8 a6,a4
      a4(1) = a6(1)
      a4(2) = a6(2)
      a4(3) = 0d0
      a4(4) = a6(6)
      return
      end subroutine reduce_6to4
c-----------------------------------------------------------------------
      subroutine reduce_4to6(a4,a6)
      dimension a6(6),a4(4)
      real*8 a6,a4
      a6(1) = a4(1)
      a6(2) = a4(2)
      a6(3) = 0d0
      a6(4) = 0d0
      a6(5) = 0d0
      a6(6) = a4(4)
      return
      end subroutine reduce_4to6
c-----------------------------------------------------------------------
      subroutine int2chr(in,chr_out)
c     Convert the given integer <in> to a character variable

c     Argument                                                   intents
c     --------                                                   -------
c     in                                                          <in>
c     chr_out                                                     <out>
      implicit none
      integer in
      character(len=*) chr_out
      character(len=10) tmp
      write(tmp,'(i10)') in
      chr_out=trim(adjustl(tmp))
      return
      end subroutine
c-----------------------------------------------------------------------
c$$$      character (len=*) function int2chr_func(in)
c$$$c
c$$$      call int2chr(in,int2chr_func)
c$$$      return
c$$$      end function
