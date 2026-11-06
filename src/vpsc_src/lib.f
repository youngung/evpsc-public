c-----------------------------------------------------------------------
      subroutine w_bdim(iunit,array,ndi,fact)!,ibr)
c     Arguments
c     ---------
c     iunit: file ID
c     array: 1-Dimensional array
c     ndi: size of the array
c     fact:multiplicative factor to scale the elements in the array
c     ibr: flag to insert line-breaker
      implicit none
      integer, intent(in) :: iunit, ndi
      dimension array(ndi)
      real*8, intent(in) ::  array
c      logical, intent(in) :: ibr
      dimension aux33(3,3),auxbb(ndi,ndi),aux3333(3,3,3,3)
      real*8 aux33,auxbb,aux3333
      real*8, intent(in) :: fact

      call chg_basis(array,aux33,auxbb,aux3333,1,ndi)
      call w_mdim(iunit,aux33,3,fact)

      return
      end subroutine w_bdim
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
c$$$      subroutine w_vdim33(iunit,array,ndi,fact)
c$$$c     Arguments
c$$$c     ---------
c$$$c     iunit: file ID
c$$$c     array: 1-Dimensional array
c$$$c     ndi: size of the array
c$$$c     fact:multiplicative factor to scale the elements in the array
c$$$      implicit none
c$$$      integer, intent(in) :: iunit, ndi
c$$$      real*8, intent(in) ::  array(3,3)
c$$$      dimension aux33(3,3),auxbb(ndi,ndi),aux3333(3,3,3,3),aux6(6)
c$$$      real*8 aux33,auxbb,aux3333,aux6
c$$$      real*8, intent(in) :: fact
c$$$
c$$$      call voigt(aux6,array,auxbb,aux3333,2)
c$$$      call w_dim(iunit,aux6,6,fact,.true.)
c$$$
c$$$      return
c$$$      end subroutine w_vdim33
c-----------------------------------------------------------------------
      subroutine w_mvdim(iunit,array,ndi,fact)
c     Arguments
c     ---------
c     iunit: file ID
c     array: 1-Dimensional array
c     ndi: size of the array
c     fact: multiplicative factor to scale the elements in the array
      implicit none
      integer, intent(in) :: iunit, ndi
      dimension array(ndi,ndi)
      real*8, intent(in) :: array
      dimension aux33(3,3),aux3333(3,3,3,3),aux6(6),
     $     aux66(6,6)
      real*8 aux33,aux3333,aux6,aux66
      real*8, intent(in) :: fact
      call chg_basis(aux6(:ndi),aux33,array,aux3333,3,ndi)
      call voigt(aux6,aux33,aux66,aux3333,4)
c      call w_mdim(iunit,aux66,ndi,fact)
      call w_mdim(iunit,aux66,6,fact)
      return
      end subroutine w_mvdim
