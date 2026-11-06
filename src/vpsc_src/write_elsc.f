c     subroutine write_elsc
c     writing self-consistent elastic compliance and modulus
      subroutine write_elsc(uw,iopt,selav,celav)
      implicit none
      integer,intent(in):: uw,iopt
      real*8, intent(in)::selav(6,6),celav(6,6)
c     locals
      integer i,j,k,l
      real*8 aux6(6),aux33(3,3),aux3333(3,3,3,3),bux3333(3,3,3,3)

      if (iopt.eq.0) then !! create file
         open(uw,action='write',file='ELSC.OUT',status='UNKNOWN')
      elseif(iopt.eq.1) then !! save ELSC to files

         call chg_basis(aux6,aux33,selav,aux3333,3,6) !! compliance
         call chg_basis(aux6,aux33,celav,bux3333,3,6) !! stiffness
         write(uw,'(e18.10,8x,36e18.10,8x,4e18.10)',
     $        advance='no')
     $        1d0/aux3333(1,1,1,1),((selav(i,j),j=1,6),i=1,6),
     $        bux3333(1,2,1,2), bux3333(1,3,1,3),bux3333(2,3,2,3),
     $        bux3333(1,3,2,3)
         do i=1,3
         do j=1,3
         do k=1,3
         do l=1,3
c           ELSC compliance
            write(uw,'(e18.10)',advance='no')aux3333(i,j,k,l)
         enddo
         enddo
         enddo
         enddo
         do i=1,3
         do j=1,3
         do k=1,3
         do l=1,3
c           ELSC stiffness
            write(uw,'(e18.10)',advance='no')bux3333(i,j,k,l)
         enddo
         enddo
         enddo
         enddo
         write(uw,*)
      endif
      return
      end subroutine
