c-----------------------------------------------------------------------
c     subroutine knuth shuffle_top_array
c     This subroutine applies the knuth shuffle of (1:n),
c     then choose elements upto top 'nchoice'.
c     That way, a random collection of "nchoice" integers among the given
c     "n" integers is obtained.
      subroutine knuth_shuffle(n,nchoice,iarr)
      implicit none
      integer, intent(in):: n,nchoice
      integer, intent(out):: iarr(nchoice)
c     locals
      real :: r
      integer :: i, j, narr(n), itemp
c     initialization
      do i=1, n
         narr(i)=i
      enddo
c     Knuth's shuttfle
      do i=n,2,-1
         call random_number(r)
         j=int(r*i)+1 ! scales it to [1,i]
         ! swap arr(i) <-> arr(j)
         itemp=narr(i)
         narr(i)=narr(j)
         narr(j)=itemp
      enddo
      iarr(:)=narr(:nchoice)
      return
      end subroutine knuth_shuffle
