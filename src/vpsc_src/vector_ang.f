c     This subroutine calculates a 3x3 rotation matrix 'R'
c     that represents a rotation of an angle 'TH' about axis defined by vector 'U'.
c     The rotation is applied using the Rodrigues' rotation formula.

c     The vector 'U' is normalized to ensure it is a unit vector.
c     The angle 'TH' is used to compute the cosine(CT) and
c     sine(ST) of the rotation angle.

c     The cross product matrix 'CX' for the vector 'U' is
c     calculated using the 'CROSSOP' subroutine, by which vector cross operation
c     is represented by the linear operation.

c     The rotation matrix 'R' is returned, which can be used to
c     rotate vectors in 3-dimensional space by the angle 'TH' around
c     the axis defined by the unit vector 'U'.
      SUBROUTINE VECTOR_ANG(U, TH, R)
      REAL*8 U(3), R(3,3)
      REAL*8 TH
      REAL*8 CT, ST
      REAL*8 MAG
      INTEGER I, J
      REAL*8 CX(3,3)  ! Cross product matrix
      REAL*8 ID(3,3)  ! Identity matrix

      ID(:,:) = 0.0D0
      DO I = 1, 3
            ID(I,I) = 1.0D0
      END DO

      ! Normalize vector U
      MAG = SQRT(U(1)**2 + U(2)**2 + U(3)**2)
      U(1) = U(1) / MAG
      U(2) = U(2) / MAG
      U(3) = U(3) / MAG

      CT = COS(TH)
      ST = SIN(TH)

      CALL CROSSOP(U, CX)

      ! Calculate Rotation Matrix
      DO I = 1, 3
       DO J = 1, 3
        R(I,J) = CT*ID(I,J) + ST*CX(I,J) + (1.0D0-CT)*(U(I)*U(J))
       END DO
      END DO
      RETURN
      END
