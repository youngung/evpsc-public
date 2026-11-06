c     from stretch tensor, calculate true (logarithmic) strain  a.k.a. hencky strain
      subroutine hencky(U6,ielem_type,true)
c     U6 is given, calculate true (logarithmic) strain.

c     Arguments
c     ---------
      implicit none
      real*8, intent(in)::U6(6) ! stretch tensor U33
      real*8, intent(out)::true(3,3)
      integer, intent(in)::ielem_type
c     locals
      real*8 engi(3,3),aux33(3,3), U33(3,3),aux_av(9),aux6(6)
      integer kro(3,3),i,j
c     kronecker delta
      kro(:,:)=0d0
      do i=1,3
         kro(i,i)=1d0
      enddo
c     stretch tensor U in 3x3
      aux6(:)=U6(:)
      call voigt_exp(aux_av,U33,aux6,3,3,ielem_type,3)
c     engineering strain
      do i=1,3
      do j=1,3
         engi(i,j)=U33(i,j)-kro(i,j)
      enddo
      enddo
c     logarithmic strain
      do i=1,3
      do j=1,3
         true(i,j)=kro(i,j)*log(1d0+engi(i,j))+(1-kro(i,j))*engi(i,j)
      enddo
      enddo
      return
      end subroutine hencky
