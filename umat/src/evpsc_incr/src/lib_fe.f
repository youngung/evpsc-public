c-----------------------------------------------------------------------
      subroutine convert_evpmod_jac(
     $     mdim,LBAR,dtime,ddsdde,ntens,ielem_type)
c     Arguments
c     ---------
c     MDIM: either 5 or 6 that determines the dimension of stress/strain in b-basis
c     LBAR: modulus tensor in b-basis in (mdim, mdim) shape
c     dtime: delta time
c     ddsdde: Jacobian to be passed following Abaqus convention
c     ntens: dimension
c     ielem_type: type of element (0: solid, 1: shell)
      implicit none
      integer,intent(in):: mdim,ntens,ielem_type
      real*8, intent(in):: LBAR(mdim,mdim), dtime
      real*8, intent(out)::ddsdde(ntens,ntens)
c     locals
      real*8 auxm(mdim),aux33(3,3),auxmm(mdim,mdim),aux3333(3,3,3,3)
      integer i,j,i0,j0,k0,l0,ijv(6,2),m,n
      DATA ((IJV(N,M),M=1,2),N=1,6)/1,1,2,2,3,3,1,2,1,3,2,3/

      call chg_basis(auxm,aux33,lbar,aux3333,3,mdim) ! b-basis to cartesian

      if (ielem_type.eq.0) then
         do i=1,ntens
c           stress components indices
            i0=ijv(i,1)
            j0=ijv(i,2)
            do j=1,ntens
c              strain components indices
               k0=ijv(j,1)
               l0=ijv(j,2)
               ddsdde(i,j)=aux3333(i0,j0,k0,l0)/dtime
            enddo
         enddo
      elseif(ielem_type.eq.1) then !! shell elements (S4R, S4 ...)
         if (mdim.eq.6) then
            ddsdde(1,1)=aux3333(1,1,1,1)/dtime
            ddsdde(2,2)=aux3333(2,2,2,2)/dtime
            ddsdde(1,2)=aux3333(1,1,2,2)/dtime
            ddsdde(2,1)=aux3333(2,2,1,1)/dtime

            ddsdde(1,3)=aux3333(1,1,1,2)/dtime
            ddsdde(3,1)=aux3333(1,2,1,1)/dtime
            ddsdde(2,3)=aux3333(2,2,1,2)/dtime
            ddsdde(3,2)=aux3333(1,2,2,2)/dtime
            ddsdde(3,3)=aux3333(1,2,1,2)/dtime
         elseif(mdim.eq.5) then
            ddsdde(1,1)=aux3333(1,1,1,1)
            ddsdde(2,2)=aux3333(2,2,2,2)
            ddsdde(1,2)=aux3333(1,1,2,2)
            ddsdde(2,1)=aux3333(2,2,1,1)

            ddsdde(1,3)=aux3333(1,1,1,2)
            ddsdde(3,1)=aux3333(1,2,1,1)
            ddsdde(2,3)=aux3333(2,2,1,2)
            ddsdde(3,2)=aux3333(1,2,2,2)
            ddsdde(3,3)=aux3333(1,2,1,2)
         endif

      endif
      return
      end subroutine
c-----------------------------------------------------------------------
      subroutine mat_transpose_fe(a,ndim,b)
c     Arguments
c     ---------
c     a    : original array
c     ndim : dimension
c     b    : resulting array
      implicit none
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
c-----------------------------------------------------------------------
      subroutine tens2rank_transf(a,b,c)
c     c_ij = b_ik b_jl a_kl
c
c     Arguments
c     ---------
c     a     : array in (3,3)
c     b     : transformation matrix in (3,3)
c     c     : resulting array in (3,3)
      implicit none
      real*8, intent(in):: a(3,3),b(3,3)
      real*8, intent(out):: c(3,3)
      integer i,j,k,l

      c(:,:)=0d0
      do i=1,3
      do j=1,3
         do k=1,3
         do l=1,3
            c(i,j)=c(i,j)+b(i,k)*b(j,l)*a(k,l)
         enddo
         enddo
      enddo
      enddo

      return
      end subroutine tens2rank_transf
c-----------------------------------------------------------------------
      subroutine tens4rank_transf(a,b,c)
c     c_ijkl = b_im b_jn b_ko b_lp a_mnop
c
c     Arguments
c     ---------
c     a     : array in (3,3,3,3)
c     b     : transformation matrix in (3,3)
c     c     : resulting array in (3,3,3,3)
      implicit none
      real*8, intent(in):: a(3,3,3,3),b(3,3)
      real*8, intent(out):: c(3,3,3,3)
      integer i,j,k,l,m,n,o,p

      c(:,:,:,:)=0d0
      do i=1,3
      do j=1,3
      do k=1,3
      do l=1,3
         do m=1,3
         do n=1,3
         do o=1,3
         do p=1,3
            c(i,j,k,l)=c(i,j,k,l)+b(i,m)*b(j,n)*b(k,o)*b(l,p)*a(m,n,o,p)
         enddo
         enddo
         enddo
         enddo
      enddo
      enddo
      enddo
      enddo

      return
      end subroutine tens4rank_transf
c-----------------------------------------------------------------------
