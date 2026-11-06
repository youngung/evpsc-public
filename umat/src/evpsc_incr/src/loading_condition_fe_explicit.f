c -------------------------------------------------------
c     Vectorized representation of F, stress, and strain
c     tensors are utilized in the code. use the below
c     table when referring to the specific components of
c     the tensorial quantities in use.
c -------------------------------------------------------
c  deformation gradient (defgradOld, defgradNew)
c -------------------------------------------------------
c component  plane stress   plane strain      3D
c                          (or axisymmetric)
c -------------------------------------------------------
c   1           F11             F11           F11
c   2           F22             F22           F22
c   3           F12             F33           F33
c   4           F21             F12           F12
c   5                           F21           F23
c   6                                         F31
c   7                                         F21
c   8                                         F32
c   9                                         F13
c -------------------------------------------------------
c stress/strain tensors (strainInc, stressOld, StressNew)
c -------------------------------------------------------
c component plane stress  plane strain 3D(VUMAT) 3D(UMAT)
c -------------------------------------------------------
c     1        11            11           11        11
c     2        22            22           22        22
c     3        12            33           33        33
c     4                      12           12        12
c     5                                   23        13
c     6                                   13        23
c -------------------------------------------------------
c  ielem_type  -1            1            0
c -------------------------------------------------------
      subroutine loading_condition_fe_explicit(
c     intent(in)
     $     iopt,ielem_type,strainInc,Fnew,Fold,stretchNew,stretchOld,dt,
     $     ndi,nshr,
c     intent(out)
     $     DBAR,dsim,scauchy,WBARc,LIJBARc,ILBAR,iDBARv,iSBARv,rot_acum)
c-----------------------------------------------------------------------
c     list of arguments
c-----------------------------------------------------------------------
c      *** intent(in)
c     <iopt>: default is 3, which determines the boundary condition
c           solely based on deformation gradient tensor.
c     <ielem_type>: element type (1 or 0)
c     <strainInc>: Given strain increment from Abaqus
c     <Fnew>: Given new deformation gradient
c     <Fold>: Given old deformation gradient
c     <StretchNew>: Given stretch tensor
c     <StretchOld>: Given stretch tensor
c     <dt>: Given time increment,
c     <ndi,nshr>: number of normal and shear components
c
c      *** intent(out)
c     <DBAR>    : strain rate to be imposed to polycrystal
c     <dsim>    : strain rate derived from LIJBARc (velocity gradient tensor)
c     <scauchy> : cauchy stress tensor (not updated within this subroutine)
c     <WBARc>  : spin
c     <LIJBARc>    : velocity gradient tensor
c     <iDBARv>   : flag that control the relaxed component of strain rate
c     <iSBARv>   : flag that control the relaxed component of stress
c     <rot_acum>: accumulative rotation
c-----------------------------------------------------------------------
      implicit none
c     input
      integer, intent(in)::ndi,nshr,iopt,ielem_type ! 0: solid, 1: shell
      real*8,intent(in)::fnew(ndi+nshr+nshr),fold(ndi+nshr+nshr),
     $     strainInc(ndi+nshr),stretchNew(ndi+nshr),stretchOld(ndi+nshr)
     $     ,dt
c     output
      real*8, intent(out):: DBAR(6),dsim(3,3),scauchy(3,3),WBARc(3,3),
     $     LIJBARc(3,3),rot_acum(3,3)
      integer, intent(out)::ILBAR(3,3),iDBARv(6),iSBARv(6)
c     local
      real*8 aux_av(ndi+nshr+nshr),aux33(3,3),bux33(3,3),cux33(3,3),
     $     dux33(3,3),eux33(3,3),aux66(6,6),aux_sv(ndi+nshr),
     $     aux3333(3,3,3,3),aux22(2,2),u(3,3),drot(3,3)

      integer i,j,k,inverr
      logical verbose
      parameter(verbose=.false.)
c      parameter(verbose=.true.)

c      call w_val(0,'dt:',dt)
c      call w_chrc(0,'strainInc')
c      call w_dim(0,strainInc,ndi+nshr,1d0,.true.)
c      call w_chrc(0,'strainInc/dt')
c      call w_dim(0,strainInc/dt,ndi+nshr,1d0,.true.)

      if (verbose) then
         write(*,*)'Inside subr. loading_condition_fe_explicit'
         write(*,*)'iopt:',iopt
c        determine deformation gradient.
c	 call w_val(0,'dt:',dt)
         write(*,'(a3,x,e20.11)')'dt:',dt
         write(*,*)'ndi:',ndi
         write(*,*)'nshr:',nshr
         write(*,*)'iopt:',iopt
         write(*,*)'ielem_type:',ielem_type
         if (ielem_type.eq.0) then
            write(*,*)'solid element'
         elseif (ielem_type.eq.1) then
            write(*,*)'shell element'
         endif
c
         write(*,*)'Fnew'
         write(*,'(9e15.7)')(fnew(i),i=1,ndi+nshr+nshr)
         write(*,*)'Fold'
         write(*,'(9e15.7)')(fold(i),i=1,ndi+nshr+nshr)
         write(*,*)'Fnew-Fold'
         write(*,'(9e15.7)')(fnew(i)-fold(i),i=1,ndi+nshr+nshr)
         call w_chrc(0,'strainInc/dt')
         call w_dim(0,strainInc/dt,ndi+nshr,1d0,.true.)


c         call xplb_exit()
      endif

      if (iopt.eq.0)then        ! Boundary condition for initial check-up of VUMAT
c     LIJBARc is based on strain increment tensor and time.
         write(*,*)'this option is deprecated'
         call xplb_exit()
         aux_sv(:)=strainInc(:)/dt
         call voigt_exp(aux_av,aux33,aux_sv,ndi,nshr,ielem_type,3)
         LIJBARc(:,:)=aux33(:,:)
         if (verbose) then
            call w_chrc(0,'strainInc')
            call w_dim(0,strainInc,ndi+nshr,1d0,.true.)
            call w_chrc(0,'LIJBARc')
            call w_mdim(0,LIJBARc,3,1d0)
c            call xplb_exit()
         endif

      elseif(iopt.eq.1) then
c
c     the usual kinematics based on deformation gradient is imposed.


         do i=1,ndi+nshr+nshr
            aux_av(i)=(Fnew(i)-Fold(i))/dt
         enddo

         if (verbose) then
            write(*,*)'aux_av'
            write(*,'(9e15.7)')aux_av
         endif

         call voigt_exp(aux_av,aux33,aux_sv,ndi,nshr,ielem_type,1)

         if (verbose) then
            write(*,*)'fdot'
            do i=1,3
               write(*,'(3f15.7)')(aux33(i,j),j=1,3)
            enddo
c            call xplb_exit()
         endif

         write(*,*)'strainInc'
         write(*,'(6e15.7)')(strainInc(i),i=1,ndi+nshr)


         do i=1,ndi+nshr+nshr
            aux_av(i)=Fnew(i)
         enddo

c        iopt 1:    aux_av -> bux33 Asymmetric Voigt to Cartesian
         call voigt_exp(aux_av,bux33,aux_sv,ndi,nshr,ielem_type,1)

         if (verbose) then
            write(*,*)'aux33 (Fdot)'
            do i=1,3
               write(*,'(3e15.7)')(aux33(i,j),j=1,3)
            enddo
            write(*,*)'bux33 (Fnew)'
            do i=1,3
               write(*,'(3e15.7)')(bux33(i,j),j=1,3)
            enddo
         endif

         if (ielem_type.eq.0) then
            cux33(:,:)=bux33(:,:)
            call lu_inverse(inverr,bux33,3) ! F^-1
            if (inverr.ne.0) then
               write(*,*)'** error in loading_condition_fe_exp.f(1)'
               call xit()
            endif
c           F^-1.F
            do i=1,3
            do j=1,3
               dux33(i,j)=0d0
               do k=1,3
                  dux33(i,j)=dux33(i,j)+cux33(i,k)*bux33(k,j)
               enddo
            enddo
            enddo

            if (verbose) then
               call w_chrc(0,'dux33')
               call w_mdim(0,dux33,3,1d0)
               call w_chrc(0,'aux33')
               call w_mdim(0,aux33,3,1d0)
               call w_chrc(0,'bux33')
               call w_mdim(0,bux33,3,1d0)
            endif

            do i=1,3
               do j=1,3
                  LIJBARc(i,j)=0d0
                  do k=1,3
                     LIJBARc(i,j)=LIJBARc(i,j)+aux33(i,k)*bux33(k,j)
                  enddo
               enddo
            enddo
         elseif(ielem_type.eq.1) then
            write(*,*)'after voigt_exp'

            do i=1,2
               do j=1,2
                  aux22(i,j)=bux33(i,j)
               enddo
            enddo
            call lu_inverse(inverr,aux22,2) ! F^-1
            if (inverr.ne.0) then
               write(*,*)'** error in loading_condition_fe_exp.f(2)'
               call xit()
            endif
            write(*,*)'Finv'
            do i=1,2
               write(*,'(2e15.7)')(aux22(i,j),j=1,2)
            enddo
            call xplb_exit()
         endif

      elseif(iopt.eq.2) then

c     inputs: fnew, fold, unew, uold

c     boundary condition defined by deformation gradient and stretch tensors.

c     following step is done.
c     1. obtain rot_acum by using f^new and u^new by applying polar decomposition
c           r^new=f^(new) . (u^(new))^-1  (rot_acum array)
c           r^old=f^(old) . (u^(old))^-1
c     2. obtain incremental rotation tensor (dr)
c          dr = r^(new) . (r^(old))       (drot array)
c     3. calculate LIJBARc (velocity gradient)
c          LIJBARc = (f^(new) - f^(old)) / dt  . (f^(new))^-1
c               =  fdot . f^-1

         aux_sv(:)=stretchNew(:)
         call voigt_exp(aux_av,u,aux_sv,ndi,nshr,ielem_type,3)
         call hencky(stretchnew,ielem_type,aux33)
         aux_av(:)=fNew(:)
         call voigt_exp(aux_av,bux33,aux_sv,ndi,nshr,ielem_type,1)
         aux33(:,:)=u(:,:)
         call lu_inverse(inverr,aux33,3) ! u^(-1)
         if (inverr.ne.0) then
            write(*,*)'** error in loading_condition_fe_exp.f(3)'
            call xit()
         endif
c     ** rot_acum transforms e^0 to e^(n+1) basis.
         call mat_dot(bux33,aux33,rot_acum,3)
c     *** old deformationi gradient tensor
c         f = r.u
c         r = f.u^-1
         aux_sv(:)=stretchOld(:)
         call voigt_exp(aux_av,u,aux_sv,ndi,nshr,ielem_type,3)
         aux_av(:)=fOld(:)
         call voigt_exp(aux_av,bux33,aux_sv,ndi,nshr,ielem_type,1)
         aux33(:,:)=u(:,:)
         call lu_inverse(inverr,aux33,3) ! u^(-1)
         if (inverr.ne.0) then
            write(*,*)'** error in loading_condition_fe_exp.f(4)'
            call xit()
         endif
c     ** cux33 transforms e^0 to e^(n) basis.
         call mat_dot(bux33,aux33,cux33,3)
c     ** incremental rotation tensor
c     dr = r^(new). (r^(old))^t
         call mat_transpose_fe(cux33,3,aux33)
         call mat_dot(rot_acum,aux33,drot,3)
c        let's calculate LIJBARc
         do i=1,ndi+nshr+nshr
            aux_av(i)=(fnew(i)-fold(i))/dt
         enddo
         call voigt_exp(aux_av,aux33,aux_sv,ndi,nshr,ielem_type,1)
         do i=1,ndi+nshr+nshr
            aux_av(i)=fnew(i)
         enddo
         call voigt_exp(aux_av,bux33,aux_sv,ndi,nshr,ielem_type,1)! bux33: fnew


         if (ielem_type.eq.0) then ! solid element
            cux33(:,:)=bux33(:,:)! cux33: fnew
            call lu_inverse(inverr,bux33,3) ! bux33: f^-1
            if (inverr.ne.0) then
               write(*,*)'** error in loading_condition_fe_exp.f(5)'
               call xit()
            endif

c            call mat_dot(cux33,bux33,dux33,3) ! fnew . f^-1
c           LIJBARc = fdot . f^-1
            call mat_dot(aux33,bux33,LIJBARc,3) ! fdot. bux33 = LIJBARc
         else                   ! shell element (maybe required to cover other types too?)
            write(*,*)'unexected case in loading_condition_fe_explicit'
            call xplb_exit()
         endif
      elseif(iopt.eq.3) then

c     inputs: fnew, fold, Unew and Uold are calculated by decompsing Fnew and Fold
c     boundary condition defined by deformation gradient...

c     Following steps are performed:
c     1. obtain r^new (rot_acum) by decomposing f^new
c     2. obtain r^old            by decomposing f^old
c          dr = r^(new) . (r^(old))^T       (drot array - not used)
c     3. calculate LIJBARc (velocity gradient, from which strain rate is calculated)
c          LIJBARc = (f^(new) - f^(old)) / dt  . (f^(new))^-1
c               =  fdot . f^-1

         aux_av(:)=fNew(:)
         call voigt_exp(aux_av,aux33,aux_sv,ndi,nshr,ielem_type,1)
c$$$         write(*,*)'fnew'
c$$$         do i=1,3
c$$$            write(*,'(3f7.2)')(aux33(i,j),j=1,3)
c$$$         enddo
         call polar_dcmp(aux33,bux33,cux33) ! bux33: R, cux33: U
c$$$         write(*,*)'r new and u new'
c$$$         do i=1,3
c$$$            write(*,'(2(3f7.2,2x))')
c$$$     $           (bux33(i,j),j=1,3),(cux33(i,j),j=1,3)
c$$$         enddo
         rot_acum(:,:)=bux33(:,:)
         aux_av(:)=fOld(:)
         call voigt_exp(aux_av,bux33,aux_sv,ndi,nshr,ielem_type,1) ! fold -> bux33
         call polar_dcmp(bux33,cux33,eux33) ! cux33: r^old
         eux33(:,:)=(aux33(:,:)-bux33(:,:))/dt ! fdot: dux33
         dux33(:,:)=eux33(:,:)
c$$$         write(*,*)'fold'
c$$$         do i=1,3
c$$$            write(*,'(3f7.2)')(aux33(i,j),j=1,3)
c$$$         enddo
c$$$         write(*,*)'r old and u old'
c$$$         do i=1,3
c$$$            write(*,'(2(3f7.2,2x))')
c$$$     $           (bux33(i,j),j=1,3),(cux33(i,j),j=1,3)
c$$$         enddo
c$$$         write(*,*)'fdot'
c$$$         do i=1,3
c$$$            write(*,'(3e15.7)')(dux33(i,j),j=1,3)
c$$$         enddo
         call mat_transpose_fe(cux33,3,eux33) !(r^old)^T
         call mat_dot(rot_acum,eux33,drot,3)
         if (ielem_type.eq.0) then ! solid element
            cux33(:,:)=aux33(:,:)!cux33: fnew
            call lu_inverse(inverr,cux33,3) ! cux33: (f^new)^-1
            if (inverr.ne.0) then
               write(*,*)'** error in loading_condition_fe_exp.f(5)'
               call xit()
            endif
c           fdot . f^-1
            call mat_dot(dux33,cux33,LIJBARc,3)!! aux33: fdot . f^-1
c            call mat_dot(aux33,bux33,LIJBARc,3)
         else                   ! shell element (maybe required to cover other types too?)
            write(*,*)'unexected case in loading_condition_fe_explicit'
            call xplb_exit()
         endif


c         write(*,*)'** LIJBARc'
c         do i=1,3
c            write(*,'(3f15.7)')(LIJBARc(i,j),j=1,3)
c         enddo


      endif

      if (verbose)then
         call w_chrc(0,'LIJBARc')
         call w_mdim(0,LIJBARc,3,1d0)
      endif

      aux33(:,:)=0d0
      bux33(:,:)=0d0
      do i=1,3
         do j=1,3
            aux33(i,j)=(LIJBARc(i,j)+LIJBARc(j,i))*0.5d0 !DBAR
            bux33(i,j)=(LIJBARc(i,j)-LIJBARc(j,i))*0.5d0
         enddo
      enddo

      do i=1,3
         do j=1,3
            ILBAR(i,j)=1
         enddo
      enddo

      iDBARv(1)=ILBAR(1,1)
      iDBARv(2)=ILBAR(2,2)
      iDBARv(3)=ILBAR(3,3)
      iDBARv(4)=ILBAR(2,3)*ILBAR(3,2)
      iDBARv(5)=ILBAR(1,3)*ILBAR(3,1)
      iDBARv(6)=ILBAR(1,2)*ILBAR(2,1)

      iSBARv(:)=0
      scauchy(:,:)=0d0

      DO I=1,3
         DO J=1,3
            DSIM(I,J)=(LIJBARc(I,J)+LIJBARc(J,I))/2d0
         ENDDO
      ENDDO

      CALL CHG_BASIS(DBAR,DSIM,AUX66,AUX3333,2,6)

      if (verbose) then
         write(*,*)'dsim'
         do i=1,3
            write(*,'(3e15.7)')(dsim(i,j),j=1,3)
         enddo
      endif


      WBARc(:,:)=0d0
      DO I=1,3
         DO J=1,3
            IF(ILBAR(I,J).EQ.1.AND.ILBAR(J,I).EQ.1) THEN
               RODBAR(I,J)=(LIJBARc(I,J)-LIJBARc(J,I))/2d0
            ELSE IF(ILBAR(I,J).EQ.1) THEN
               RODBAR(I,J)=LIJBARc(I,J)-DSIM(I,J)
               RODBAR(J,I)=-RODBAR(I,J)
            ENDIF
            LIJBARc(I,J)=DSIM(I,J)+RODBAR(I,J)
         enddo
      enddo

      return
      end subroutine
