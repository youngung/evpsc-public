      subroutine calc_devm(dstran,ntens,devm)
      implicit none
      integer, intent(in)::ntens
      real*8,intent(in)::dstran(ntens)
      real*8,intent(out)::devm
c
      integer i

      devm=0d0

      do i=1,ntens
         devm=devm+dstran(i)**2
      enddo
      devm=sqrt(2d0/3d0*devm)


      return
      end subroutine calc_devm
