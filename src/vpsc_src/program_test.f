c     (compile)
c     gfortran -c src/vpsc_src/program_test.f -o objdir/program_test.o
c     (link)
c     gfortran objdir/crystal_symmetry.o objdir/fake_xit.o objdir/program_test.o objdir/crossop.o objdir/euler.o objdir/gaus.o objdir/vector_ang.o objdir/phi_delta2vect.o -o euichan

c ==============================================================
C Program: euichan
C Purpose:
C   - Determine crystal symmetry from an SX file and calculate
C     related Cartesian vectors.
C   - Read Miller indices and Burgers vector, and compute the
C     corresponding normal vectors.
C Compilation Instructions:
C   - Compile: gfortran -c src/vpsc_src/program_test.f -o objdir/program_test.o
C   - Link: gfortran objdir/crystal_symmetry.o objdir/fake_xit.o objdir/program_test.o -o euichan
C Inputs:
C   - SX file containing material data (e.g., DP440_test.sx)
C Outputs:
C   - Crystal symmetry type, Cartesian vectors, Miller indices,
C     Burgers vector, and their normal vectors.
C ==============================================================

      PROGRAM euichan
      end program euichan

      SUBROUTINE gaussian_dispersion(sx_file, N, isn, isb, sigma, A_rot)

C   Variable Declarations

      implicit none
      character*100 sx_file
      integer ioption, icrysym, UR_SX, nequiv, iopt, isn(4), isb(4)
      integer i, j, k, l, N
      real*8 sn(3), sneq(3,24), sb(3), CVEC(3,3), Mat(3,3), t(3)
      real*8 A(3,3), A_rot(3,3)
      real*8 phi1, Phi, phi2
      real*8, allocatable:: R(:,:,:)

c ======================================================
c     Step 1 : Open the SX file
      sx_file = 'matData/vpscData/DP440/DP440_test.sx'

      open(unit=UR_SX, file=sx_file, status='old', action='read')

c ======================================================
c   Step 2 : calculate Cartesian Vector (CVEC)
      ioption = 1

      call crystal_symmetry(ioption, UR_SX, icrysym, isn, sn, sneq,
     &                      isb, sb, nequiv)
      print *, 'Crystal symmetry type (icrysym): ', icrysym
c ========================================================
c   Step 3 : calculate normal vector (sn, sb)
      ioption = 2

!      isn(1)=1
!      isn(2)=2
!      isn(3)=3
!      isb(1)=1
!      isb(2)=1
!      isb(3)=-1

!      print *, 'Miller indices (isn): ', isn(1), isn(2), isn(3)
!      print *, 'Burgers vector (isb): ', isb(1), isb(2), isb(3)

      call crystal_symmetry(ioption, UR_SX, icrysym, isn, sn, sneq,
     &                      isb, sb, nequiv)

      print *, 'Normal vector for Miller indices (sn): '
     &          ,sn(1),sn(2),sn(3)
      print *, 'Normal vector for Burgers vector (sb): '
     &          ,sb(1),sb(2),sb(3)

      close(UR_SX)
c ============================================================
c   Step 4 : calculate t = sn X sb
      call CROSSOP(sn, Mat)
      t(:) = 0.0D0
      do i = 1, 3
         do j = 1, 3
            t(i) = t(i) + Mat(i, j) * sb(j)
         enddo
      enddo

      print *, 't = sn X sb : ', t(1), t(2), t(3)
c =============================================================
c     Step 5 : Calculate transformation matrix [A]
      do i = 1, 3
         A(i, 1) = sb(i)
         A(i, 2) = t(i)
         A(i, 3) = sn(i)
      enddo
      print *, 'Transformation matrix [A] : '
      do i = 1, 3
         print *, (A(i,j), j=1, 3)
      enddo
c ==============================================================
c     Step 6 : Transform [A] to Euler angles
      iopt=1

      call euler(iopt, phi1, Phi, phi2, A)
      print *, 'Euler angle (phi1, Phi, phi2 in degrees ): '
      print *, 'phi1   =', phi1
      print *, 'Phi    =', Phi
      print *, 'phi2   =', phi2
c ===============================================================
c     Step 7 : Apply gaussian random rotations
      mu = 0.0D0
!      sigma = 0.3D0
!      N = 1000
      allocate(R(N,3,3))

      call gaus(N, mu, sigma, R)
c =================================================================
c     Step 8 : Apply rotaiton to A and compute A_rot
      open(unit = 10, file='D:/jupyter/tex.out', status='replace')

      do l = 1, N
         do i = 1,3
               do j = 1,3
                     A_rot(i,j)=0.0D0
                     do k = 1,3
                        A_rot(i,j)=A_rot(i,j)+A(i,k)*R(l,k,j)
                     enddo
               enddo
         enddo


!         print *, 'Rotated Transforamtion matrix [A_rot]: '
!        do i = 1,3
!               print*, (A_rot(i,j),j = 1,3)
!         enddo

c     Compute Euler angles from rotated transformation matrix [A_rot]

         iopt =1
         CALL euler(iopt, phi1, Phi, phi2, A_rot)

         write(10, '(3F10.4)') phi1, Phi, phi2

!         PRINT *, 'Euler angles (phi1, Phi, phi2) for rotation ', l, ':'
!         PRINT *, 'phi1   =', phi1
!         PRINT *, 'Phi    =', Phi
!         PRINT *, 'phi2   =', phi2
!         print *, ' '
      enddo

      close(10)
      deallocate(R)
      END SUBROUTINE gaussian_dispersion
