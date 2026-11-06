c----------------------------------------------------------------------c
      subroutine update_eacc_explicit(
c     intent(in)
     $     ielem_type,dstran,ntens,rot_acum_t,
c     intent(out)
     $     deacc_B)
c     ielem_type: 0: solid, 1: shell
c     dstran: abaqus incremental strain (Voigt) referred in the co-rotating coordinates
c     ntens
c     rot_acum_t: transforms from e0 basis vectors to e(n) basis vectors
c     deacc_B: strain increment in B-basis
      implicit none

      integer, intent(in)::ielem_type,ntens
      real*8, intent(in)::dstran(ntens),rot_acum_t(3,3)
      real*8, intent(out)::deacc_B(6)
c      real*8, intent(inout)::
c     locals
      real*8 aux6(6),aux33(3,3),bux33(3,3),aux66(6,6),aux3333(3,3,3,3)
      logical verbose
      parameter(verbose=.false.)
c      parameter(verbose=.true.)
      integer i,j

      if (verbose) then
         write(*,*)'** beginning of update_eacc_explicit'
         write(*,*)'ielem_type:',ielem_type
         write(*,*)'dstran:'
         write(*,'(6f8.5)')(dstran(i),i=1,ntens)
         write(*,*)'ntens:',ntens
         write(*,*)'rot_acum_t:'
         do i=1,3
            write(*,'(3f5.2)')(rot_acum_t(i,j),j=1,3)
         enddo
      endif

      if (ielem_type.eq.0) then ! solid
         call voigt_abq(dstran,bux33,aux66,aux3333,1)
      elseif(ielem_type.eq.1) then ! shell
         aux6(:)=0d0            ! voigt
         aux6(1)=dstran(1)      ! deps in co-rotating basis
         aux6(2)=dstran(2)      ! deps in co-rotating basis
         aux6(4)=dstran(3)      ! deps in co-rotating basis
         call voigt_abq(aux6,bux33,aux66,aux3333,1)
      endif

      call tens2rank_transf(bux33,rot_acum_t,aux33) ! dstran in e^0 basis vectors
      call chg_basis(aux6,aux33,aux66,aux3333,2,6) ! aux6 in b-basis
      deacc_B(:)=aux6(:)
      return
      end subroutine update_eacc_explicit
