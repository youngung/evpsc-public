c ----------------------------------------------------------------------
c     Vectorized representation of F, stress, and strain
c     tensors are utilized in the code. use the below
c     table when referring to the specific components of
c     the tensorial quantities in use.
c ----------------------------------------------------------------------
c  deformation gradient (defgradOld, defgradNew)
c ----------------------------------------------------------------------
c component  plane stress   plane strain      3D
c                          (or axisymmetric)
c ----------------------------------------------------------------------
c   1           F11             F11           F11
c   2           F22             F22           F22
c   3           F12             F33           F33
c   4           F21             F12           F12
c   5                           F21           F23
c   6                                         F31
c   7                                         F21
c   8                                         F32
c   9                                         F13
c ----------------------------------------------------------------------
c stress/strain tensors (strainInc, stressOld, StressNew)
c ----------------------------------------------------------------------
c component plane stress  plane strain 3D(VUMAT) 3D(UMAT)
c ----------------------------------------------------------------------
c     1        11            11           11        11
c     2        22            22           22        22
c     3        12            33           33        33
c     4                      12           12        12
c     5                                   23        13
c     6                                   13        23
c ----------------------------------------------------------------------
c  ielem_type  -1            1            0
c ----------------------------------------------------------------------
c     Subroutine is meant to be used in VUMAT to convert the components
c     referred in symmetric or asymmetric Voigt notation and
c     Cartesian notation
      subroutine voigt_exp(aux_av,aux33,aux_sv,ndi,nshr,ielem_type,iopt)
c     ielem_type==0 (solid, 3D)
c     iopt 1:    aux_av(9) -> aux33(3,3)    Asymmetric Voigt to Cartesian
c     iopt 2:    aux33(3,3) -> aux_av(9)    Cartesian to Asymmetric Voigt
c     iopt 3:    aux_sv(6) -> aux33(3,3)    Symmetric Voigt to Cartesian
c     iopt 4:    aux33(3,3) -> aux_sv(6)    Cartesian to Symmetric Voigt
c ----------------------------------------------------------------------
c     ielem_type==1 (shell)
c     iopt 1:    aux_av(5) -> aux33(3,3)    Asymmetric Voigt to Cartesian
c     iopt 2:    aux33(3,3) -> aux_av(5)    Cartesian to Asymmetric Voigt
c     iopt 3:    aux_sv(4) -> aux33(3,3)    Symmetric Voigt to Cartesian
c     iopt 4:    aux33(3,3) -> aux_sv(4)    Cartesian to Symmetric Voigt
      implicit none
      integer, intent(in):: iopt,ndi,nshr,ielem_type
      real*8, intent(inout):: aux_av(ndi+nshr+nshr),aux33(3,3),
     $     aux_sv(ndi+nshr)
c     local
      integer i
      logical verbose
c      parameter(verbose=.true.)
      parameter(verbose=.false.)

      if (verbose) then
         write(*,*)'** entering voigt_exp'
         write(*,*)'ndi:',ndi
         write(*,*)'nshr:',nshr
         write(*,*)'iopt:',iopt
         write(*,*)'aux_av'
         write(*,'(9e15.7)')(aux_av(i),i=1,ndi+nshr+nshr)
         write(*,*)'aux_sv'
         write(*,'(9e15.7)')(aux_sv(i),i=1,ndi+nshr)
         write(*,*)'ielem_type:',ielem_type
      endif

      if ((ndi.ne.3.or.nshr.ne.3).and.ielem_type.eq.0) then
         write(*,*)'* unexpected combination of <ndi,nshr,ielem_type>
     $in <voigt_exp>'
         write(*,*)'ndi,nshr,ielem_type:',ndi,nshr,ielem_type
         call xplb_exit
      elseif ((ndi.ne.3.or.nshr.ne.1).and.ielem_type.eq.1) then
         write(*,*)'** unexpected combination of <ndi,nshr,ielem_type>
     $in <voigt_exp>'
         write(*,*)'ndi,nshr,ielem_type:',ndi,nshr,ielem_type
         call xplb_exit
      endif


      if (ielem_type.eq.0) then ! 3D
         if (iopt.eq.1)then ! e.g., F9 -> F33
            aux33(1,1)=aux_av(1)
            aux33(2,2)=aux_av(2)
            aux33(3,3)=aux_av(3)
c     shear components
            aux33(1,2)=aux_av(4)
            aux33(2,3)=aux_av(5)
            aux33(3,1)=aux_av(6)
            aux33(2,1)=aux_av(7)
            aux33(3,2)=aux_av(8)
            aux33(1,3)=aux_av(9)
         elseif(iopt.eq.2)then ! e.g., F33 -> F9
            aux_av(1)=aux33(1,1)
            aux_av(2)=aux33(2,2)
            aux_av(3)=aux33(3,3)
c     shear components
            aux_av(4)=aux33(1,2)
            aux_av(5)=aux33(2,3)
            aux_av(6)=aux33(3,1)
            aux_av(7)=aux33(2,1)
            aux_av(8)=aux33(3,2)
            aux_av(9)=aux33(1,3)
         elseif(iopt.eq.3)then ! e.g., S6->S33
            aux33(1,1)=aux_sv(1)
            aux33(2,2)=aux_sv(2)
            aux33(3,3)=aux_sv(3)
c     shear components
            aux33(1,2)=aux_sv(4)
            aux33(2,3)=aux_sv(5)
            aux33(1,3)=aux_sv(6)

            aux33(2,1)=aux_sv(4)
            aux33(3,2)=aux_sv(5)
            aux33(3,1)=aux_sv(6)
         elseif(iopt.eq.4)then ! e.g., S33->S6
            aux_sv(1)=aux33(1,1)
            aux_sv(2)=aux33(2,2)
            aux_sv(3)=aux33(3,3)
c     shear components
            aux_sv(4)=aux33(1,2)
            aux_sv(5)=aux33(2,3)
            aux_sv(6)=aux33(1,3)
         endif
      elseif(ielem_type.eq.1) then
c     ----------------------------------------------------------------
c     ielem_type==1 (ndir=3, nshr=1)
c     iopt 1:    aux_av(5) -> aux33(3,3)    Asymmetric Voigt to Cartesian
c     iopt 2:    aux33(3,3) -> aux_av(5)    Cartesian to Asymmetric Voigt
c     iopt 3:    aux_sv(4) -> aux33(3,3)    Symmetric Voigt to Cartesian
c     iopt 4:    aux33(3,3) -> aux_sv(4)    Cartesian to Symmetric Voigt
         if (iopt.eq.1) then ! e.g., F5 -> F33
            aux33(:,:)=0d0
c     Asymmetric Voigt to Cartesian
            aux33(1,1)=aux_av(1)
            aux33(2,2)=aux_av(2)
            aux33(3,3)=aux_av(3)
            aux33(1,2)=aux_av(4)
            aux33(2,1)=aux_av(5)
         elseif (iopt.eq.2) then ! e.g., F33 -> F5
            aux_av(:)=0d0
c     Cartesian to Asymmetric Voigt
            aux_av(1)=aux33(1,1)
            aux_av(2)=aux33(2,2)
            aux_av(3)=aux33(3,3)
            aux_av(4)=aux33(1,2)
            aux_av(5)=aux33(2,1)
         elseif (iopt.eq.3) then ! e.g., S4 -> S33
            aux33(:,:)=0d0
c     Symmetric Voigt to Cartesian
            aux33(1,1)=aux_sv(1)
            aux33(2,2)=aux_sv(2)
            aux33(3,3)=aux_sv(3)
            aux33(1,2)=aux_sv(4)
            aux33(2,1)=aux_sv(4)
         elseif (iopt.eq.4) then !e.g., S33 -> S4
            aux_sv(:)=0d0
c     Cartesian to Asymmetric Voigt
            aux_sv(1)=aux33(1,1)
            aux_sv(2)=aux33(2,2)
            aux_sv(3)=aux33(3,3)
            aux_sv(4)=aux33(1,2)
         endif
      endif

      if (verbose) then
         write(*,*)'** Exiting voigt_exp'
      endif

      return
      end subroutine voigt_exp
