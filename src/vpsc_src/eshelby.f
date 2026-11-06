C ***********************************************************************
C     SUBROUTINE ESHELBY      --->      VERSION 15/NOV/07
C
C     IOPTION=0: Initialize arrays assoc. with Gauss integration points.
C     IOPTION=1: Calculate elastic Eshelby tensor for elastic inclusion.
C     IOPTION=2: Calculate incompressible Eshelby tensors ESIM (strain-
C                rate) & ESCR (spin-rate) associated with the visco-
C                plastic inclusion.
C     IOPTION=3: Calculate incompressible and hydrostatic Eshelby tensors
C                PESH (deviatoric pressure), PDIL (spherical pressure) &
C                DESH (dilatation) for visco-plastic inclusion.
C     IOPTION=4: Calculates 1st term in d(S)/d(M)=d(T)/d(M)*L+T*d(L)/d(M)
C     IOPTION=5: Calculates 2nd term in d(S)/d(M)=d(T)/d(M)*L+T*d(L)/d(M)
C
C     Options 2-3-4-5 admit a non-zero visco-plastic bulk modulus KEFF.
C
C     Algorithms are based in Lebensohn et al, MSMSE 6 (1998) p.447.
C
C     Uses explicit matrix inversion and explicit Voigt notation (when
C     possible) to optimize computing time.
C
C     Modified oct/2005 to adapt number of integration points to the shape
C     of the ellipsoid in order to keep Eshelby tensor within a certain
C     tolerance (based on analysis done by Gwenaelle Proust).
C     Aspect ratio criterion was adopted for the case AXIS(2)>AXIS(1)>AXIS(3)
C
C     Modified oct/2007 to use a Gauss-Lobatto integration with fix (ten)
C     integration points and weights, chosen to optimize the integration
C     within (eleven) different domains of the ellipsoid aspect ratios.
C     (based on analysis made by Laurent Capolungo).
C        if IGAUSSLEG=0 uses Gauss Lobatto  (cases 1 to 11) (hardwired)
C        if IGAUSSLEG=1 uses Gauss Legendre (case 12) (kept as an option)
C ***********************************************************************

      SUBROUTINE ESHELBY (axis,c4,keff,esim,escr,desh,pesh,pdil,dldm,
     $     dsddm,ioption,ifail)
c      use,intrinsic ::ieee_arithmetic
      implicit none
c     Arguments
c     ---------
c     axis:  ratio of ellipoid (intent: in)
c     c4  :                    (intent: in)
c     pesh : P_ij:L_ijkl^(tan)
c     desh : (P_ij + P_ji)/2.

      real*8 c4(3,3,3,3),esim(3,3,3,3),escr(3,3,3,3)
      real*8 p(3,3),pesh(3,3),desh(3,3)
      real*8 c2(6,6),gamma2(6,6),gamma4(3,3,3,3)
      real*8 axis(3),x1(10),a1(10),a1inv(10)
      real*8 aa1x(6),aa2x(3,3),aaww1x(6),aaww2x(3,3)

      real*8 dldm(3,3,3,3),dsddm(3,3,3,3),dldm2(6,6)
CFEB
      real*8 da1(6),da(4,4),dainv(4,4),ainv(4,4)
      real*8 aux6(6),aux33(3,3),aux66(6,6),aux3333(3,3,3,3),aux44(4,4)
CFEE
      integer ngaumx, ngaumx2
      PARAMETER (ngaumx=40,ngaumx2=1600)
      real*8 xph(ngaumx),xth(ngaumx),wph(ngaumx),wth(ngaumx)
      integer ngaussph, ngaussth
      real*8 alpha,aa1,aww,aaww1,ww
      COMMON/kESHELBY1/ngaussph(12),ngaussth(12)
      COMMON/kESHELBY2/alpha(12,3,ngaumx2),aa1(12,6,ngaumx2),
     $     aww(12,3,ngaumx2),aaww1(12,6,ngaumx2),ww(12,ngaumx2)

C *** Integration points and weights, and limits of aspect ratios
      real*8 punti, pesi, dte,puntigl,pesigl
      COMMON/kESHELBY3/punti(10,12),pesi(10,12),dte(0:10),
     $     puntigl(ngaumx),pesigl(ngaumx)

      REAL*8     keff
      INTEGER   case,i,ny,j,l,n,m,k,iph,ith,ioption,ngau_c12,npoints,
     $     igaussleg,ierr_inv
      real*8 dumscr,dumsim,dummy,pdil,ro3,ratio2,abcoro3,abc,ratio1
     $     ,costh,sinth,pi,simbtet
c      logical isnan
      logical ifail
c,
      ifail=.false.

      pi=4.*atan(1.0)

      IGAUSSLEG=0     ! hardwires the Gauss-Lobatto option
C     IGAUSSLEG=1     ! hardwires the Gauss-Legendre option
c     make sure that ngau_c12 for Gauss-Leg is .le. ngaumx
      NGAU_c12 =40

C ***********************************************************************
C     INITIALIZATION RUN
C     Calculates Gauss-Legendre integration points and weights in the
c     interval [0,pi]. Initializes Gauss-Lobatto points and weights.
C     Initializes arrays associated with each point to avoid repeating
C     its calculation at every call. All values of the points and
C     weights were calculated to optimize the error.
C ***********************************************************************

      if(ioption.eq.0) then

         do i=1,11
            ngaussph(i)=10      ! Gauss-Lobatto
            ngaussth(i)=10
         end do
         ngaussph(12)=ngau_c12  ! Gauss-Legendre
         ngaussth(12)=ngau_c12

c****************************
c  CASE 1
c****************************
         punti(1,1)=4.71236594E-02
         punti(2,1)=0.241774723d0
         punti(3,1)=0.565131843d0
         punti(4,1)=0.968887568d0
         punti(5,1)=1.37937832d0
         punti(6,1)=1.76221442d0
         punti(7,1)=2.17270517d0
         punti(8,1)=2.57646084d0
         punti(9,1)=2.89981818d0
         punti(10,1)=3.09446883d0

         pesi(1,1)=0.120191820d0
         pesi(2,1)=0.264987558d0
         pesi(3,1)=0.373805553d0
         pesi(4,1)=0.420841277d0
         pesi(5,1)=0.390970200d0
         pesi(6,1)=0.390970260d0
         pesi(7,1)=0.420841366d0
         pesi(8,1)=0.373805553d0
         pesi(9,1)=0.264987499d0
         pesi(10,1)=0.120192111d0

c****************************
c  CASE 2
c****************************
         punti(1,2)=1.57080423d-02
         punti(2,2)=0.144995824d0
         punti(3,2)=0.425559640d0
         punti(4,2)=0.829968274d0
         punti(5,2)=1.31460333d0
         punti(6,2)=1.82698941d0
         punti(7,2)=2.31162453d0
         punti(8,2)=2.71603298d0
         punti(9,2)=2.99659705d0
         punti(10,2)=3.12588477d0

         pesi(1,2)=5.41692823E-02
         pesi(2,2)=0.207461149d0
         pesi(3,2)=0.348739326d0
         pesi(4,2)=0.452716887d0
         pesi(5,2)=0.507709801d0
         pesi(6,2)=0.507709682d0
         pesi(7,2)=0.452716798d0
         pesi(8,2)=0.348738998d0
         pesi(9,2)=0.207461327d0
         pesi(10,2)=5.41692935E-02

c****************************
c  CASE 3
c****************************
         punti(1,3)=3.76990959E-02
         punti(2,3)=0.198626831d0
         punti(3,3)=0.483041346d0
         punti(4,3)=0.871647120d0
         punti(5,3)=1.32964790d0
         punti(6,3)=1.81194484d0
         punti(7,3)=2.26994562d0
         punti(8,3)=2.65855122d0
         punti(9,3)=2.94296598d0
         punti(10,3)=3.10389376d0

         pesi(1,3)=9.68142375E-02
         pesi(2,3)=0.224478707d0
         pesi(3,3)=0.341134071d0
         pesi(4,3)=0.430180043d0
         pesi(5,3)=0.478189558d0
         pesi(6,3)=0.478189170d0
         pesi(7,3)=0.430180043d0
         pesi(8,3)=0.341134191d0
         pesi(9,3)=0.224478647d0
         pesi(10,3)=9.68143344E-02

c****************************
c  CASE 4
c****************************
         punti(1,4)=3.45576368E-02
         punti(2,4)=0.187556863d0
         punti(3,4)=0.468425453d0
         punti(4,4)=0.859980166d0
         punti(5,4)=1.32527423d0
         punti(6,4)=1.81631863d0
         punti(7,4)=2.28161263d0
         punti(8,4)=2.67316723d0
         punti(9,4)=2.95403576d0
         punti(10,4)=3.10703516d0

         pesi(1,4)=8.95763785E-02
         pesi(2,4)=0.217725381d0
         pesi(3,4)=0.341026783d0
         pesi(4,4)=0.435772508d0
         pesi(5,4)=0.486694932d0
         pesi(6,4)=0.486695170d0
         pesi(7,4)=0.435772508d0
         pesi(8,4)=0.341026902d0
         pesi(9,4)=0.217725128d0
         pesi(10,4)=8.95764604E-02

c****************************
c  CASE 5
c****************************
         punti(1,5)= 3.14158052E-02
         punti(2,5)=0.177928671
         punti(3,5)= 0.457155794
         punti(4,5)= 0.851592362
         punti(5,5)= 1.32222414
         punti(6,5)= 1.81936860
         punti(7,5)=2.29000044
         punti(8,5)=2.68443704
         punti(9,5)=2.96366405
         punti(10,5)=3.11017680

         pesi(1,5)=8.26927349E-02
         pesi(2,5)=0.213228315
         pesi(3,5)=0.342008322
         pesi(4,5)=0.440196186
         pesi(5,5)=0.492670894
         pesi(6,5)=0.492670983
         pesi(7,5)=0.440195888
         pesi(8,5)=0.342008322
         pesi(9,5)=0.213227972
         pesi(10,5)=8.26930404E-02

c****************************
c  CASE 6
c****************************
         punti(1,6)= 2.98452154E-02
         punti(2,6)=0.173592165
         punti(3,6)=0.452448040
         punti(4,6)=0.848216832
         punti(5,6)=1.32101476
         punti(6,6)=1.82057810
         punti(7,6)= 2.29337597
         punti(8,6)=2.68914461
         punti(9,6)=2.96800065
         punti(10,6)=3.11174774

         pesi(1,6)=7.93928578E-02
         pesi(2,6)=0.211627841
         pesi(3,6)=0.342669785
         pesi(4,6)=0.442057431
         pesi(5,6)=0.495048553
         pesi(6,6)=0.495048642
         pesi(7,6)=0.442057490
         pesi(8,6)=0.342670023
         pesi(9,6)=0.211627468
         pesi(10,6)=7.93929026E-02

c****************************
c  CASE 7
c****************************
         punti(1,7)=2.67036632E-02
         punti(2,7)=0.165752888
         punti(3,7)=0.444431901
         punti(4,7)=0.842614472
         punti(5,7)=1.31902647
         punti(6,7)= 1.82256627
         punti(7,7)=2.29897833
         punti(8,7)=2.69716072
         punti(9,7)=2.97583985
         punti(10,7)=3.11488938

         pesi(1,7)=7.30879456E-02
         pesi(2,7)=0.209402516
         pesi(3,7)=0.344104946
         pesi(4,7)=0.445234656
         pesi(5,7)=0.498966068
         pesi(6,7)= 0.498966306
         pesi(7,7)=0.445234746
         pesi(8,7)= 0.344104946
         pesi(9,7)=0.209402665
         pesi(10,7)=7.30878562E-02

c****************************
c  CASE 8
c****************************
         punti(1,8)=2.67036632E-02
         punti(2,8)=0.165752888
         punti(3,8)=0.444431901
         punti(4,8)=0.842614472
         punti(5,8)=1.31902647
         punti(6,8)=1.82256627
         punti(7,8)=2.29897833
         punti(8,8)=2.69716072
         punti(9,8)=2.97583985
         punti(10,8)=3.11488938

         pesi(1,8)=7.30879456E-02
         pesi(2,8)=0.209402516
         pesi(3,8)=0.344104946
         pesi(4,8)=0.445234656
         pesi(5,8)=0.498966068
         pesi(6,8)=0.498966306
         pesi(7,8)=0.445234746
         pesi(8,8)=0.344104946
         pesi(9,8)=0.209402665
         pesi(10,8)= 7.30878562E-02

c****************************
c  CASE 9
c****************************
         punti(1,9)=2.43473575E-02
         punti(2,9)=0.160516247
         punti(3,9)=0.439386278
         punti(4,9)=0.839168847
         punti(5,9)=1.31781363
         punti(6,9)=1.82377899
         punti(7,9)=2.30242372
         punti(8,9)=2.70220637
         punti(9,9)=2.98107672
         punti(10,9)=3.11724544

         pesi(1,9)=6.86219111E-02
         pesi(2,9)=0.208388865
         pesi(3,9)=0.345189095
         pesi(4,9)=0.447236270
         pesi(5,9)=0.501360059
         pesi(6,9)=0.501359940
         pesi(7,9)=0.447236151
         pesi(8,9)=0.345189214
         pesi(9,9)=0.208388969
         pesi(10,9)=6.86219335E-02

c****************************
c  CASE 10
c****************************
         punti(1,10)=2.19910536E-02
         punti(2,10)=0.155757755
         punti(3,10)=0.434985727
         punti(4,10)=0.836206555
         punti(5,10)=1.31677616
         punti(6,10)= 1.82481658
         punti(7,10)=2.30538607
         punti(8,10)=2.70660710
         punti(9,10)=2.98583508
         punti(10,10)=3.11960149

         pesi(1,10)=6.43825606E-02
         pesi(2,10)=0.207786217
         pesi(3,10)=0.346235514
         pesi(4,10)=0.448981822
         pesi(5,10)=0.503410578
         pesi(6,10)= 0.503410578
         pesi(7,10)=0.448981792
         pesi(8,10)=0.346235693
         pesi(9,10)=0.207785636
         pesi(10,10)= 6.43827692E-02

c****************************
c  CASE 11
c****************************
         punti(1,11)=2.04204638E-02
         punti(2,11)=0.152822554
         punti(3,11)=0.432348520
         punti(4,11)=0.834448099
         punti(5,11)=1.31616223
         punti(6,11)=1.82543063
         punti(7,11)=2.30714464
         punti(8,11)=2.70924401
         punti(9,11)=2.98877001
         punti(10,11)=3.12117243

         pesi(1,11)=6.16818815E-02
         pesi(2,11)=0.207559645d0
         pesi(3,11)=0.346902698d0
         pesi(4,11)=0.450027168d0
         pesi(5,11)=0.504624724d0
         pesi(6,11)= 0.504624426d0
         pesi(7,11)=0.450027317d0
         pesi(8,11)=0.346902847d0
         pesi(9,11)=0.207559645d0
         pesi(10,11)=6.16819337E-02

C *****************************************************************
C  CASE 12: GAULEG generates the integration points
C ****************************************************************

         call gauleg(0.d0,pi,puntigl,pesigl,ngaussph(12))

C *****************************************************************
C *** Calculates and saves arrays that depend on integration points

         do case=1,12

            if (case.eq.12) then ! Gauss-Lobatto
               do i=1,ngau_c12
                  xph(i)=puntigl(i)
                  xth(i)=puntigl(i)
                  wph(i)= pesigl(i)
                  wth(i)= pesigl(i)
               end do
            else                 ! Gauss-Legendre
               do i=1,10
                  xph(i)=punti(i,case)
                  xth(i)=punti(i,case)
                  wph(i)= pesi(i,case)
                  wth(i)= pesi(i,case)
               end do
            end if
c--------------------------------------------------------------
c *** integration [0,pi][0,pi] adds a factor 2 in Eqs. B11 & B14.

            do ith=1,ngaussth(case)
               sinth=sin(xth(ith))
               costh=cos(xth(ith))
               simbtet=wth(ith)*sinth/(2d0*pi) ! added the factor of 2

            do iph=1,ngaussph(case)
               ny=iph+(ith-1)*ngaussph(case)
               ww(case,ny)=simbtet*wph(iph)
               alpha(case,1,ny)=sinth*cos(xph(iph))
               alpha(case,2,ny)=sinth*sin(xph(iph))
               alpha(case,3,ny)=costh

               do i=1,3
               do j=1,3
                  aa2x(i,j)  =alpha(case,i,ny)*alpha(case,j,ny)
                  aaww2x(i,j)=aa2x(i,j)*ww(case,ny)
               enddo
               enddo
c              aa2x(3,3)->aa1x(6)
               call voigt(aa1x  ,aa2x  ,aux66,aux3333,2)
c              aaww2x(3,3)->aaww1x(6)
               call voigt(aaww1x,aaww2x,aux66,aux3333,2)
               do i=1,6
                  aa1(case,i,ny)  =aa1x(i)
                  aaww1(case,i,ny)=aaww1x(i)
c                  if (aaww1(case,i,ny).eq.0)then
                  if (abs(aaww1(case,i,ny)).le.1e-25)then
                     write(*,*)'Warning: zero aaww1 may cause problems'
                  endif
               enddo

c *** Array AWW is used only if ICAUCHY=1.
c            do i=1,3
c              aww(case,i,ny)=alpha(case,i,ny)*ww(case,ny)
c            enddo

            enddo ! end of iph
            enddo ! end of ith
         enddo                  ! end of do case=1,12
      endif                     ! ENDIF FOR IOPTION=0

C************************************************************************
C     End of initialization
C************************************************************************

C ***********************************************************************
C     CALCULATION OF ESHELBY TENSORS FOR STIFFNESS 'C4' AND ELLIPSOID
C     AXES 'AXIS'
C     ASSUMES: AXIS2 > AXIS1 > AXIS3  --> RATIO1 > RATIO2 > 1
C ***********************************************************************

      if(ioption.ge.1) then

         abc=axis(1)*axis(2)*axis(3)
         ratio1=axis(2)/axis(3)
         ratio2=axis(1)/axis(3)

         if (igaussleg.eq.1) then
            case=12
         else
            dte(0) = 0d0
            dte(1) =-0.7d0*ratio1+7d0
            dte(2) =-ratio1+17d0
            dte(3) =-ratio1+23d0
            dte(4) =-ratio1+26d0
            dte(5) =-ratio1+29.3d0
            dte(6) =-ratio1+32d0
            dte(7) =-ratio1+34.85d0
            dte(8) =-ratio1+37d0
            dte(9) =-ratio1+41.9d0
            dte(10)=-ratio1+44.5d0
            case=11
            do i=1,10
               if(ratio2.ge.dte(i-1) .and. ratio2.lt.dte(i)) case=i
            enddo
         endif

c--------------------------------------------------------

         npoints=ngaussph(case)*ngaussth(case)

         pdil=0d0
         p(1:3,1:3)=0d0
         gamma2(1:6,1:6)=0d0

c       c4   (3,3,3,3) -> c2    (6,6)
         call voigt(aa1x,aa2x,c2,c4,4)
c       dldm (3,3,3,3) -> dldm2 (6,6)
         IF(IOPTION.EQ.5) CALL VOIGT(AA1X,AA2X,DLDM2,DLDM,4)

         do ny=1,npoints

c   Compute Voigt components A1(1)-A(6) of tensor A(3,3) defined by Eq.B3:
c   --->  A(i,j)=L(i,j,k,l)*a(j)*a(l)

            do i=1,6
               aa1x(i)=aa1(case,i,ny) ! (1x6) Voigt matrix of a(j)*a(l)
            enddo
c           alpha_l alpha_j L_ijkl
c     A (1x6) matrix can be obtained from a(i) a(j), saves to aa1x
c     Following the Voigt notation: aa1x(6) L(6,6) -> a1
            do i=1,6
               aux6(i)=a1(i)
            enddo
c            call esh_mult_voigt(c2,aa1x,a1(1:6)) ! c2*aa1x=a1
            call esh_mult_voigt(c2,aa1x,aux6) ! c2*aa1x=a1
            do i=1,6
               a1(i)=aux6(i)
            enddo

c     Therefore, a1(1:6) stores L(i,j,k,l) * a(j) * a(l)

            IF(IOPTION.EQ.1) THEN

c   If solving an elastic inclusion invert the system
c   --> A(3,3) x X(3,3) = C(3,3)
c   Inverts A(3,3) using explicit Voigt notation.
c   Uses explicit form of C(3,3) to calculate solution in Voigt notation.

               call esh_inv3_voigt(a1,a1inv,ierr_inv)
               if (ierr_inv.ne.0) then
                  ifail=.true.
                  return
               endif
               do i=1,6
                  x1(i)=a1inv(i)
               enddo




c               write(*,*)'------------------------'
c               call w_chrc(0,'x1')
c               write(*,'(10e12.4)')(x1(i),i=1,10)
c               call xit()

c     x1 is the inverse of A matrix in Eq. B6

            ELSE IF(IOPTION.GE.2) THEN

c   If solving a visco-plastic inclusion defines components A1(7) to A1(10).
c   Solves the system given by Eq.B4 --> A(4,4) x X(4,4) = C(4,4)
c   Inverts A(4,4) using explicit Voigt notation.
c   Uses explicit form of C(4,4) to calculate solution in Voigt notation.
c   The solution is symmetric. Numerical deviation from symmetry is averaged.

c            |   A11   A12   A13   aph1  |    | a1  a6  a5  a7  |
c            |   A21   A22   A23   aph2  |    |     a2  a4  a8  |
c            |   A31   A32   A33   aph3  |    |         a3  a9  |
c            |  aph1  aph2  aph3      0  |    |             a10 |

               a1(7) = alpha(case,1,ny)
               a1(8) = alpha(case,2,ny)
               a1(9) = alpha(case,3,ny)
               a1(10)= 0d0
               if(keff.gt.0.) a1(10)=-1d0/keff

               call esh_inv4_voigt(a1,a1inv) ! (A4)^-1

               do i=1,10
                  x1(i)=a1inv(i)
               enddo

            ENDIF               ! ioption.ge.2

CFEB
            IF(IOPTION.EQ.5) THEN

               CALL VOIGT10(A1INV,AINV,1)

               call esh_mult_voigt(dldm2,aa1x,da1) ! dldm2*aa1x=da1
               call voigt(da1,aux33,aux66,aux3333,1)

               do i=1,3
               do j=1,3
                  da(i,j)=aux33(i,j)
               enddo
               enddo

c
               do i=1,3
                  da(i,4)=0d0
                  da(4,i)=0d0
               enddo
               da(4,4)=0d0
c
c           dAinv/dp = - Ainv  : dA/dp : Ainv
c                x1  = - a1inv :  da1  : a1inv
c
               do i=1,4
               do j=1,4
                  dummy=0d0
                  do k=1,4
                     dummy=dummy+da(i,k)*ainv(k,j)
                  enddo
                  aux44(i,j)=dummy
               enddo
               enddo


               do i=1,4
               do j=1,4
                  dummy=0d0
                  do k=1,4
                     dummy=dummy+ainv(i,k)*aux44(k,j)
                  enddo
                  dainv(i,j)=-dummy
               enddo
               enddo

               CALL VOIGT10(X1,DAINV,2)

            ENDIF               ! ioption==5
c            rhobar = (axis_1 alpha_1, axis_2 alpha_2, axis_3 alpha_3)
CFEE        !rho^3?
            ro3=((alpha(case,1,ny)*axis(1))**2+
     $           (alpha(case,2,ny)*axis(2))**2+
     $           (alpha(case,3,ny)*axis(3))**2)**1.5d0
c           abcoro3 = a*b*c / (rho^3)
            abcoro3=abc/ro3

            if (abs(ro3).le.1e-20) then
               write(*,*)'cannot divide by zero!'
               write(*,*)'ro3:',ro3
               write(*,*)'(alpha(case,:,ny)',alpha(case,:,ny)
               write(*,*)'axis(:)',axis(:)
               call xit()
            endif

c   Compute the Eshelby integral Eq.B11 defining:
c         Gamma(m,j,n,i)=T(m,n,i,j)=a(m)*a(j)*G(n,i)
c   with the property:
c         Gamma(m,j,n,i)=Gamma(j,m,n,i)=Gamma(m,j,i,n)

c     in the Voigt notation, the above is expressed as
c         Gamma(i,j) = a(i) G(j)
            do i=1,6
            do j=1,6
               gamma2(i,j)=gamma2(i,j)+aaww1(case,i,ny)*x1(j)*abcoro3
            enddo
            enddo

c            call xit()

c   Compute the pressure related Eshelby integral Eq.B14
            if(ioption.eq.3) then
               do j=1,3
               do i=1,3
                  p(i,j)=p(i,j)+aww(case,j,ny)*x1(i+6)*abcoro3
               enddo
               enddo
               pdil=pdil+ww(case,ny)*x1(10)*abcoro3
            endif

         enddo                  ! end of loop over double integration

c ********************************************************************
c   Go back to the 3*3*3*3 notation
         call voigt(aa1x,aa2x,gamma2,gamma4,3)


c   Compute symmetric (distortion) Eshelby tensor from Eq.B9.
c       esim(n,m,k,l)=0.5*(gamma(m,j,n,i)+gamma(n,j,m,i))*c4(i,j,k,l)
c   Compute anti-symmetric (rotation) Eshelby tensor from Eq.B9.
c       escr(n,m,k,l)=0.5*(gamma(m,j,n,i)-gamma(n,j,m,i))*c4(i,j,k,l)

         do l=1,3
         do k=1,3
         do m=1,3
         do n=1,3
c
            dumsim=0d0
            dumscr=0d0
c
            do j=1,3
            do i=1,3
c
               IF(IOPTION.NE.4) THEN
                  dummy=(gamma4(m,j,n,i)+gamma4(n,j,m,i))
     $                 *c4(i,j,k,l)
                  dumsim=dumsim+dummy
                  dummy=(gamma4(m,j,n,i)-gamma4(n,j,m,i))
     $                 *c4(i,j,k,l)
                  dumscr=dumscr+dummy
               ELSE
                  dumsim=dumsim+(gamma4(m,j,n,i)+gamma4(n,j,m,i))
     $                 *dldm(i,j,k,l)
                  dumscr=dumscr+(gamma4(m,j,n,i)-gamma4(n,j,m,i))
     $                 *dldm(i,j,k,l)
               ENDIF

            enddo
            enddo
c
            IF(IOPTION.LT.4) THEN
               esim(n,m,k,l)=0.5d0*dumsim
               escr(n,m,k,l)=0.5d0*dumscr
            ELSE
               dsddm(n,m,k,l)=0.5d0*dumsim

            ENDIF
c
         enddo
         enddo
         enddo
         enddo

c         call xit()

c   Compute pressure & dilatation related Eshelby tensors (Eq.B13)

         if(ioption.eq.3) then
            do l=1,3
            do k=1,3
               pesh(k,l)=0d0
            do j=1,3
            do i=1,3
               pesh(k,l)=pesh(k,l)+p(i,j)*c4(i,j,k,l)
            enddo
            enddo
            enddo
            enddo
            do j=1,3
            do i=1,3
               desh(i,j)=(p(i,j)+p(j,i))/2d0
            enddo
            enddo
         endif

      endif                     !  endif for IOPTION.GE.1

      RETURN
      END subroutine eshelby
C
C *************************************************************************
C     SUBROUTINE ESH_INV3_VOIGT   --->   version 23/jul/01
C
C     Inverts the 3x3 symmetric matrix 'A' using explicit Voigt notation:
C     11->1, 22->2, 33->3, 23=32->4, 31=13->5, 12=21->6
C *************************************************************************

      SUBROUTINE ESH_INV3_VOIGT (A,AINV,ierr)
      implicit none
      integer ierr
      real*8 A(10),AINV(10), det
      DET = A(1)*A(2)*A(3) + 2d0*A(4)*A(5)*A(6) - A(1)*A(4)*A(4)
     $     - A(2)*A(5)*A(5) - A(3)*A(6)*A(6)
      ierr=0

      if (abs(det).le.1e-30) then
         ierr=1
         return
c$$$         write(*,*)'** det is too small **'
c$$$         write(*,*)'det:',det
c$$$         write(*,*)'a(1:6):',a(1:6)
c$$$         write(*,*)'-----------------------------------------------'
c$$$         write(*,*)'A(1)*A(2)*A(3):',A(1)*A(2)*A(3)
c$$$         write(*,*)'2d0*A(4)*A(5)*A(6):',2d0*A(4)*A(5)*A(6)
c$$$         write(*,*)'A(1)*A(4)*A(4):',A(1)*A(4)*A(4)
c$$$         write(*,*)'A(2)*A(5)*A(5):',A(2)*A(5)*A(5)
c$$$         write(*,*)'A(3)*A(6)*A(6):',A(3)*A(6)*A(6)
c$$$         write(*,*)'==============================================='
c         call stop_debug(0,'det is too small')
      endif

      if (abs(det).le.1d-26) det=1d-26   !!just to avoid det=0 case.

      AINV(1) = ( A(2)*A(3) - A(4)*A(4))/DET
      AINV(2) = ( A(1)*A(3) - A(5)*A(5))/DET
      AINV(3) = ( A(1)*A(2) - A(6)*A(6))/DET
      AINV(4) = (-A(1)*A(4) + A(5)*A(6))/DET
      AINV(5) = ( A(4)*A(6) - A(2)*A(5))/DET
      AINV(6) = (-A(3)*A(6) + A(4)*A(5))/DET

      RETURN
      END
C
C **********************************************************************
C     SUBROUTINE ESH_INV4_VOIGT   --->   VERSION 20/JUL/01

C     Inverts the 4*4 symmetric matrix 'A' using explicit Voigt notation:
C     11-->1, 22-->2, 33-->3, 23=32-->4, 31=13-->5, 12=21-->6
C     14-->7, 24-->8, 34-->9, 44-->10.
C **********************************************************************

      SUBROUTINE ESH_INV4_VOIGT (A,AINV)
      implicit none

      real*8 A(10),AINV(10),det
      integer i
      ainv(1) = a(2)*a(3)*a(10)+2*a(4)*a(8)*a(9) -
     $     a(2)*a(9)*a(9)-a(3)*a(8)*a(8)-a(4)*a(4)*a(10)

      ainv(2) = a(1)*a(3)*a(10)+2*a(5)*a(7)*a(9) -
     $     a(1)*a(9)*a(9)-a(3)*a(7)*a(7)-a(5)*a(5)*a(10)

      ainv(3) = a(1)*a(2)*a(10)+2*a(6)*a(7)*a(8) -
     $     a(1)*a(8)*a(8)-a(2)*a(7)*a(7)-a(6)*a(6)*a(10)

      ainv(4) = a(1)*a(4)*a(10)+a(5)*a(7)*a(8)+a(6)*a(7)*a(9) -
     $     a(1)*a(8)*a(9)-a(4)*a(7)*a(7)-a(5)*a(6)*a(10)
      ainv(4) =-ainv(4)

      ainv(5) = a(4)*a(6)*a(10)+a(2)*a(7)*a(9)+a(5)*a(8)*a(8) -
     $     a(4)*a(7)*a(8)-a(6)*a(8)*a(9)-a(2)*a(5)*a(10)

      ainv(6) = a(3)*a(6)*a(10)+a(5)*a(8)*a(9)+a(4)*a(7)*a(9) -
     $     a(3)*a(7)*a(8)-a(6)*a(9)*a(9)-a(4)*a(5)*a(10)
      ainv(6) =-ainv(6)

      ainv(7) = a(4)*a(6)*a(9)+a(4)*a(5)*a(8)+a(2)*a(3)*a(7) -
     $     a(4)*a(4)*a(7)-a(2)*a(5)*a(9)-a(3)*a(6)*a(8)
      ainv(7) =-ainv(7)

      ainv(8) = a(1)*a(4)*a(9)+a(5)*a(5)*a(8)+a(3)*a(6)*a(7) -
     $     a(4)*a(5)*a(7)-a(5)*a(6)*a(9)-a(1)*a(3)*a(8)

      ainv(9) = a(1)*a(2)*a(9)+a(5)*a(6)*a(8)+a(4)*a(6)*a(7) -
     $     a(2)*a(5)*a(7)-a(6)*a(6)*a(9)-a(1)*a(4)*a(8)
      ainv(9) =-ainv(9)

      ainv(10)=a(1)*a(2)*a(3)+2*a(4)*a(5)*a(6) -
     $     a(1)*a(4)*a(4)-a(2)*a(5)*a(5)-a(3)*a(6)*a(6)

      det=   a(1)*ainv(1)+   a(2)*ainv(2)+   a(3)*ainv(3)+
     $     2.*a(4)*ainv(4)+2.*a(5)*ainv(5)+2.*a(6)*ainv(6)+
     $     2.*a(7)*ainv(7)+2.*a(8)*ainv(8)+2.*a(9)*ainv(9)+
     $     a(10)*ainv(10)
      det=det/4.

      do i=1,10
         ainv(i)=ainv(i)/det
      enddo

      return
      end
C
C ***********************************************************
C     SUBROUTINE ESH_MULT_VOIGT   -->   version 16/MAY/2012
C ***********************************************************
      SUBROUTINE ESH_MULT_VOIGT(B,C,A)
      implicit none

C     Performs the multiplication:
C        A(i,k)=B(i,j,k,l)*C(j,l) using Voigt's notation
C        B is a 6*6 symmetric matrix
C        C is a 3*3 symmetric tensor
C        A will be a 3*3 symmetric tensor

      real*8 B(6,6),C(6),A(6)

      A(1)=B(1,1)*C(1)+B(6,6)*C(2)+B(5,5)*C(3)
     #    +2*(B(5,6)*C(4)+B(1,5)*C(5)+B(1,6)*C(6))

      A(2)=B(6,6)*C(1)+B(2,2)*C(2)+B(4,4)*C(3)
     #    +2*(B(2,4)*C(4)+B(4,6)*C(5)+B(2,6)*C(6))

      A(3)=B(5,5)*C(1)+B(4,4)*C(2)+B(3,3)*C(3)
     #    +2*(B(3,4)*C(4)+B(3,5)*C(5)+B(4,5)*C(6))

      A(4)=B(5,6)*C(1)+B(2,4)*C(2)+B(3,4)*C(3)
     #      +(B(2,3)+B(4,4))*C(4)
     #      +(B(3,6)+B(4,5))*C(5)
     #      +(B(4,6)+B(2,5))*C(6)

      A(5)=B(1,5)*C(1)+B(4,6)*C(2)+B(3,5)*C(3)
     #      +(B(3,6)+B(4,5))*C(4)
     #      +(B(1,3)+B(5,5))*C(5)
     #      +(B(1,4)+B(5,6))*C(6)

      A(6)=B(1,6)*C(1)+B(2,6)*C(2)+B(4,5)*C(3)
     #      +(B(4,6)+B(2,5))*C(4)
     #      +(B(1,4)+B(5,6))*C(5)
     #      +(B(1,2)+B(6,6))*C(6)

      RETURN
      END
