c     This subroutine calculates the linear operator corresponding to
c     a vector cross product between, for example, vectors a and b: i.e., a x b
c     The cross product can be represented by a linear operation:
c              | 0  -a3  a2 |  | b1 |
c     a x b  = | a3  0  -a1 |  | b2 |
c              |-a2  a1  0  |  | b3 |

c     The matrix corresponding to [a x] operation is then:
c             | 0  -a3  a2 |
c     MAT =   | a3  0  -a1 |
c             |-a2  a1  0  |
c     This matrix is returned as 'MAT' argument by this subroutine, provided that
c     the argument 'a' is the vector.

      SUBROUTINE CROSSOP(a, Mat)
      implicit none
      REAL*8, INTENT(IN) :: a(3)
      REAL*8, INTENT(OUT) :: Mat(3,3)

      Mat = 0.0D0
      Mat(1, 2) = -a(3)
      Mat(1, 3) =  a(2)
      Mat(2, 1) =  a(3)
      Mat(2, 3) = -a(1)
      Mat(3, 1) = -a(2)
      Mat(3, 2) =  a(1)

      RETURN
      END SUBROUTINE CROSSOP

c-----------------------------------------------------------------------
      subroutine crossproduct(a,b,c)
      real*8, intent(in):: a(3),b(3)
      real*8, intent(out):: c(3)
      c(1)= a(2)*b(3)-a(3)*b(2)
      c(2)=-a(1)*b(3)+a(3)*b(1)
      c(3)= a(1)*b(2)-a(2)*b(1)
      return
      end subroutine
