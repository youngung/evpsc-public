c     This subroutine generates random angles Phi and Delta,
c     then converts the vector in the spherical coordinates (phi, delta, 1)
c     that of Cartesian coordinates.

c     - Phi is randomly chosen in the range [-pi, pi].
c     - Delta is randomly chosen in the range [0, pi].

c     The random number generation uses the RANDOM_NUMBER function, and
c     the seed for the random number generator is initialized by calling
c     RANDOM_SEED.

      SUBROUTINE phi_delta2vect_subroutine(Vector)
      implicit none
      REAL*8, INTENT(OUT) :: Vector(3)
      REAL*8 :: R, PI, Phi, Delta
      REAL*8 :: X, Y, Z, Norm

c      CALL RANDOM_SEED()
      PI=4d0*ATAN(1d0)

c     Generate random Phi in the range [-pi, pi]
      CALL RANDOM_NUMBER(R)
      Phi = 2.0D0 * PI * R - PI

c     Generate random Delta in the range [0, pi]
      CALL RANDOM_NUMBER(R)
      Delta = 2.0D0 * R - 1.0D0
      Delta = ACOS(Delta)

c     Convert spherical vector (Phi, Delta, 1) to Cartesian vector (X, Y, Z)
      X = COS(Phi)
      Y = SIN(Phi)
      Z = SIN(Delta)

c     Normalize the Cartesian vector
      Norm = SQRT(X**2 + Y**2 + Z**2)
      Vector = [X, Y, Z] / Norm

      RETURN
      END SUBROUTINE phi_delta2vect_subroutine
