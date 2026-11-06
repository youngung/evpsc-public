c     This subroutine generates 'N' 3x3 rotation matrix 'R'.
c     The rotation matrix is based on random orientations
c     and angle derived from a normal distribution with mean 'MU'
c     and standard deviation 'SIGMA'.

c     The steps include:
c     - Generating random spherical coordinates (Phi and Delta)
c       and converting them into a unit vector using the 'phi_delta_to_vect_subroutine'.
c     - Computing a random rotation angle 'TH' using a normal
c       distribution derived from the Box-Muller method.
c     - Generating a 3x3 rotation matrix for each iteration using
c       the 'VECTOR_ANG' subroutine. This matrix rotates vectors
c       around the random axis defined by the unit vector.
c     - Writing each rotation matrix to a file ('rot_mat.txt')

      SUBROUTINE gaus(N, MU, SIGMA, R)
      implicit none
      INTEGER, INTENT(IN) :: N
      REAL*8, INTENT(IN) :: MU, SIGMA
      REAL*8, INTENT(OUT) :: R(N, 3, 3)
      REAL*8 :: Vector(3), TH, R1, R2, PI
      INTEGER :: I, J

      PI = 3.14159265358979323846D0
c      CALL RANDOM_SEED()

      OPEN(UNIT=10, FILE='rot_mat.txt', STATUS='UNKNOWN')

      DO I = 1, N
c     Generate random spherical coordinates and convert them to a unit vector
          CALL phi_delta2vect_subroutine(Vector)

c     Generate random numbers for Box-Muller Transform
          CALL RANDOM_NUMBER(R1)
          CALL RANDOM_NUMBER(R2)

c     Box-Muller Transform: Converts two uniformly distributed
c     random numbers (R1, R2) into a normally distributed
c     random variable with mean MU and standard deviation SIGMA.
          TH=MU+SIGMA*SQRT(-2.*LOG(R1))*COS(2.*PI*R2)
          TH = MODULO(TH + PI, 2.0D0 * PI) - PI

c     Generate the rotation matrix using the unit vector and random angle
          CALL VECTOR_ANG(Vector, TH, R(I,:,:))

          WRITE(10, *) 'Rotation Matrix for Iteration:', I
          DO J = 1, 3
              WRITE(10, *) R(I, J, 1), R(I, J, 2), R(I, J, 3)
          END DO
          CONTINUE
      END DO

      CLOSE(10)

      RETURN
      END SUBROUTINE gaus
