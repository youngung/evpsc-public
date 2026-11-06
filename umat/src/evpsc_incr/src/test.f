      program main
      real*8 f(3,3),r(3,3),u(3,3)
c$$$      f(1,1)=1
c$$$      f(1,2)=0.495
c$$$      f(1,3)=0.5
c$$$      f(2,1)=-0.333
c$$$      f(2,2)=1
c$$$      f(2,3)=-0.247
c$$$      f(3,1)=0.959
c$$$      f(3,2)=0
c$$$      f(3,3)=1.5


      f(1,1)=1.
      f(1,2)=2.
      f(1,3)=2.
      f(2,1)=4.
      f(2,2)=4.
      f(2,3)=4.
      f(3,1)=4.
      f(3,2)=1.
      f(3,3)=4.
      call polar_dcmp(f,r,u)
      end program


      subroutine mat_transpose_fe(a,ndim,b)
c     Arguments
c     ---------
c     a    : original array
c     ndim : dimension
c     b    : resulting array
c      implicit none
      integer, intent(in)::ndim
      real*8, intent(in)::a(ndim,ndim)
      real*8, intent(out)::b(ndim,ndim)
c     -- local
      integer i,j
c     -----------
      do i=1,ndim
      do j=1,ndim
         b(i,j)=a(j,i)
      enddo
      enddo
      return
      end subroutine mat_transpose_fe


      include "polar_dcmp.f"
      include "../../../../src/vpsc_src/lu_libs.f"
      include "../../../../src/vpsc_src/library7.sub"
