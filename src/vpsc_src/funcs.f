C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     FUNCTIONS OF VPSC
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      real*8 function det(a)
      implicit none
      real*8 a(3,3)

      det=a(1,1)*a(2,2)*a(3,3)
     $   +a(1,2)*a(2,3)*a(3,1)
     $   +a(1,3)*a(2,1)*a(3,2)
     $   -a(1,3)*a(2,2)*a(3,1)
     $   -a(2,3)*a(3,2)*a(1,1)
     $   -a(1,2)*a(2,1)*a(3,3)
      return
      end

C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     FUNCTION VNORM    ---->   VERSION OF 23/NOV/2011
C
C     CALCULATES THE NORM OF A VECTOR V WITH N COMPONENTS
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      real*8 FUNCTION VNORM (V,N)
      implicit none
      integer i,n
      real*8 V(N)


      VNORM=0d0
      DO I=1,N
        VNORM=VNORM+V(I)*V(I)
      ENDDO
      VNORM=SQRT(VNORM)

      RETURN
      END

C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     FUNCTION TNORM5x5    ---->   VERSION OF 23/NOV/2011
C
C     CALCULATES THE NORM OF A 5x5 MATRIX
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      real*8 FUNCTION TNORM5x5 (T)
      implicit none
      real*8 T(25),tnorm
      integer i

      TNORM=0.
      DO I=1,25
        TNORM=TNORM+T(I)*T(I)
      ENDDO
      TNORM5x5=SQRT(TNORM)

      RETURN
      END

C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     FUNCTION TNORM6x6    ---->   VERSION OF 23/NOV/2011
C
C     CALCULATES THE NORM OF A 6x6 MATRIX
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      real*8 FUNCTION TNORM6x6 (T)
      implicit none
      real*8 T(36),tnorm
      integer i

      TNORM=0.
      DO I=1,36
        TNORM=TNORM+T(I)*T(I)
      ENDDO
      TNORM6x6=SQRT(TNORM)

      RETURN
      END

C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     FUNCTION TNORM    ---->   VERSION OF 05/sep/2019
C
C     CALCULATES THE NORM OF A NxN MATRIX
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      real*8 FUNCTION TNORM (T,N)
      implicit none
      integer, intent(in)::n
      real*8, intent(in):: T(N*N)
c     locals
      integer i

      TNORM=0.
      DO I=1,N*N
        TNORM=TNORM+T(I)*T(I)
      ENDDO
      TNORM=SQRT(TNORM)

      RETURN
      END

C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     FUNCTION VMISMATCH   ---->   VERSION OF 23/NOV/2011
C
C     CALCULATES RELATIVE DIFFERENCE BETWEEN TWO VECTORS WITH N COMPONENTS.
C     THE DIFFERENCE IS RELATIVE TO THE NORM OF THE ARITHMETIC AVERAGE.
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      real*8 FUNCTION VMISMATCH (V1,V2,N)
      implicit none
      integer i,n
      real*8 V1(N),V2(N),V_dif(N),V_ave(N),vnorm


      DO I=1,N
        V_dif(I)=V1(I)-V2(I)
        V_ave(I)=0.5d0*(V1(I)+V2(i))
      ENDDO
      if (VNORM(V_ave,N).ne.0) then
         VMISMATCH =VNORM(V_dif,N)/VNORM(V_ave,N)
      else
         vmismatch=0d0
      endif
      RETURN
      END


c$$$C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c$$$C     FUNCTION WVMISMATCH   ----> AP:  VERSION OF 29March23
c$$$C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c$$$
c$$$      real*8 FUNCTION WVMISMATCH (V1,V2,N)
c$$$      implicit none
c$$$      integer, intent(in):: n
c$$$      real*8, intent(in):: V1(N),V2(N)
c$$$c     locals
c$$$      integer i
c$$$      real*8 V_dif(N),V_ave(N),error,vmax
c$$$
c$$$      vmax = 1d-30
c$$$
c$$$      DO I=1,N
c$$$         if (abs(V1(I)) .gt. vmax) vmax = abs(V1(I))
c$$$         if (abs(V2(I)) .gt. vmax) vmax = abs(V2(I))
c$$$      END DO
c$$$
c$$$      V_dif(:) = 0.e0
c$$$      error = 0.e0
c$$$      DO I=1,N
c$$$         V_dif(I)=((abs(V1(I))/vmax)*(V1(I)-V2(I))/
c$$$     $        (V1(I)+V2(I))/0.5)**2
c$$$         error = error + V_dif(I)
c$$$! V_ave(I)=0.5*(V1(I)+V2(i))
c$$$      ENDDO
c$$$      WVMISMATCH = sqrt(error)/N
c$$$
c$$$      RETURN
c$$$      END

C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     FUNCTION TMISMATCH5x5   ---->   VERSION OF 23/NOV/2011
C
C     CALCULATES RELATIVE DIFFERENCE BETWEEN TWO 5x5 MATRICES.
C     THE DIFFERENCE IS RELATIVE TO THE NORM OF THE ARITHMETIC AVERAGE.
C     REPLACES PREVIOUS 'TMISMATCH' AND 'TMISMATCH5'
C     MAKE SURE THAT DIMENSION OF INPUT ARRAY IN CALLING ROUTINE IS 5x5
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      real*8 FUNCTION TMISMATCH5x5 (T1,T2)
      implicit none
      real*8 T1(25),T2(25),T_dif(25),T_ave(25),vnorm
      integer i

      DO I=1,25
        T_dif(I)=T1(I)-T2(I)
        T_ave(I)=0.5*(T1(I)+T2(i))
      ENDDO
      TMISMATCH5x5=VNORM(T_dif,25)/VNORM(T_ave,25)

      RETURN
      END

C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     FUNCTION TMISMATCH6x6   ---->   VERSION OF 23/NOV/2011
C
C     CALCULATES RELATIVE DIFFERENCE BETWEEN TWO 6x6 MATRICES.
C     THE DIFFERENCE IS RELATIVE TO THE NORM OF THE ARITHMETIC AVERAGE.
C     REPLACES PREVIOUS 'TMISMATCH'
C     MAKE SURE THAT DIMENSION OF INPUT ARRAY IN CALLING ROUTINE IS 6x6
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      real*8 FUNCTION TMISMATCH6x6 (T1,T2)
      implicit none
      real*8 vnorm
      real*8 T1(36),T2(36)
      real*8 T_dif(36),T_ave(36)
      integer i

      DO I=1,36
         T_dif(I)=T1(I)-T2(I)
         T_ave(I)=0.5d0*(T1(I)+T2(i))
      ENDDO
      TMISMATCH6x6=VNORM(T_dif,36)/VNORM(T_ave,36)

      RETURN
      END

C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     FUNCTION TMISMATCH   ---->   VERSION OF 13/SEP/2019
C
C     CALCULATES RELATIVE DIFFERENCE BETWEEN TWO NxN MATRICES.
C     THE DIFFERENCE IS RELATIVE TO THE NORM OF THE ARITHMETIC AVERAGE.
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      real*8 FUNCTION TMISMATCH (T1,T2,N)
      implicit none
      integer, intent(in):: n
      real*8, intent(in):: T1(N*N),T2(N*N)
c     locals
      integer i
      real*8 vnorm,T_dif(N*N),T_ave(N*N)


      DO I=1,N*N
        T_dif(I) =T1(I)-T2(I)
        T_ave(I)=0.5*(T1(I)+T2(i))
      ENDDO
      TMISMATCH=VNORM(T_dif,N*N)/VNORM(T_ave,N*N)

      RETURN
      END


c----------------------------------------------------------------------c
      subroutine mat6rot(a,rot,b)
      implicit none
      real*8,intent(in)::a(6),rot(3,3)
      real*8,intent(out)::b(6)
      real*8 aux33(3,3),
     $     bux33(3,3),aux3333(3,3,3,3),aux66(6,6)
cf2py intent(in) a, rot
cf2py intent(out) b
      call chg_basis(a,aux33,aux66,aux3333,1,6)
      call matrot(aux33,rot,bux33)
      call chg_basis(b,bux33,aux66,aux3333,2,6)
      return
      end subroutine mat6rot
c----------------------------------------------------------------------c
c     Transform a matrix by rotation mat rot and returns the result as b
c     b = rot  a  rot^T
c     b_ij = rot_ik a_kl rot_jl

c     rot{b_axes<-a_axes}
      subroutine matrot(a,rot,b)
      implicit none
      real*8 a(3,3), rot(3,3), b(3,3)
      integer i,j,k,l
cf2py intent(in) a, rot
cf2py intent(out) b
      do i=1,3
      do j=1,3
         b(i,j) = 0.d0
      do k=1,3
      do l=1,3
         b(i,j) = b(i,j) + rot(i,k) * a(k,l) * rot(j,l)
      enddo
      enddo
      enddo
      enddo
      return
      end subroutine matrot

c-----------------------------------------------------------------------
c     c_ij = a_ik b_kj
c
c     Arguments
c     ---------
c     a    : array
c     b    : array
c     c    : resulting array
c     ndim : dimension

      subroutine mat_dot(a,b,c,ndim)
      implicit none
      integer, intent(in)::ndim
      real*8, intent(in)::a(ndim,ndim),b(ndim,ndim)
      real*8, intent(out)::c(ndim,ndim)
c     - locals
      integer i,j,k
c     -----
      do i=1,ndim
      do j=1,ndim
         c(i,j)=0d0
      do k=1,ndim
         c(i,j)=c(i,j)+a(i,k)*b(k,j)
      enddo
      enddo
      enddo
      return
      end subroutine mat_dot
c-----------------------------------------------------------------------
      subroutine matvect_dot(a,b,c,ndim)
c     c_i = a_ij b_j
      implicit none
      integer, intent(in):: ndim
      real*8, intent(in):: a(ndim,ndim),b(ndim)
      real*8, intent(out):: c(ndim)
c     locals
      integer i,j
      do i=1,ndim
         c(i)=0.
         do j=1,ndim
            c(i)=c(i)+a(i,j)*b(j)
         enddo
      enddo
      return
      end subroutine matvect_dot
c-----------------------------------------------------------------------
      subroutine transpose(a,ndim,b)
c     Arguments
c     ---------
c     a    : original array
c     ndim : dimension
c     b    : resulting array
      implicit none
      integer, intent(in)::ndim
      real*8, intent(in)::a(ndim,ndim)
      real*8, intent(out)::b(ndim,ndim)
c     locals
      integer i,j
      do i=1,ndim
      do j=1,ndim
         b(i,j)=a(j,i)
      enddo
      enddo
      return
      end subroutine transpose
