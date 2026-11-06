c     This subroutine is to obtain drot.
      subroutine calc_drot_explicit(
c     intent(in)
     $     ielem_type,ndir,nshr,stretchNew,stretchOld,fold,fnew,
c     intent(out)
     $     drot)
      implicit none
      integer,intent(in):: ielem_type,ndir,nshr
      real*8,intent(in)::stretchNew(ndir+nshr),stretchOld(ndir+nshr),
     $     fold(ndir+nshr+nshr),fnew(ndir+nshr+nshr)
      real*8,intent(out):: drot(3,3)
c     -- locals
c     - av: asymmetric voigt
c     - sv: symmetric voigt
      real*8 aux_av(ndir+nshr+nshr),aux_sv(ndir+nshr),u(3,3),aux33(3,3),
     $     bux33(3,3),cux33(3,3),rot_acum(3,3)
      integer i,j,k,inverr
      logical verbose
      parameter(verbose=.false.)
c      parameter(verbose=.true.)

      if (verbose) then
         write(*,*)'** beginning of calc_drot_explicit'
      endif

c     *** New deformation gradient tensor
c         F = R.U
c         R = F.U^-1
      aux_sv(:)=stretchNew(:)
c         * voigt_exp (iopt=3) converts aux_sv(6) -> u(3,3) if ielem_type=0 (solid)
c                              converts aux_sv(4) -> u(3,3) if ielem_type=1 (shell)
      call voigt_exp(aux_av,u,aux_sv,ndir,nshr,ielem_type,3)
      if (verbose) then
         write(*,*)'** strechnew'
         do i=1,3
            write(*,'(3e11.3)')(u(i,j),j=1,3)
         enddo
         write(*,*)'** calling <hecky> inside <calc_drot_explicit>'
      endif
      call hencky(stretchNew,ielem_type,aux33)
      if (verbose) then
         write(*,*)'** exitng hencky'
         write(*,*)'true strain in cartesian'
         write(*,*)aux33
c$$$         write(*,*)'strain increment'
c$$$         write(*,*)strainInc
         write(*,*)'fnew'
         write(*,*)fnew
         write(*,*)'Fnew'
         write(*,'(9e15.7)')(fnew(i),i=1,ndir+nshr+nshr)
      endif
      aux_av(:)=fnew(:)
      call voigt_exp(aux_av,bux33,aux_sv,ndir,nshr,ielem_type,1)
      if (verbose) then
         write(*,*)'Deformation gradient new'
         write(*,*)bux33
      endif
      aux33(:,:)=u(:,:)
      call lu_inverse(inverr,aux33,3)  ! U^(-1)
            if (inverr.ne.0) then
               write(*,*)'** error in cacl_drot_explicit.f(1)'
               call xit()
            endif
c     ** rot_acum transforms e^0 to e^(n+1) basis.
      call mat_dot(bux33,aux33,rot_acum,3)
c     *** Old deformationi gradient tensor
c     F = R.U
c     R = F.U^-1
      aux_sv(:)=stretchOld(:)
      call voigt_exp(aux_av,u,aux_sv,ndir,nshr,ielem_type,3)
      if (verbose) then
         write(*,*)'strechOld'
         do i=1,3
            write(*,'(3e11.3)')(u(i,j),j=1,3)
         enddo
      endif
      aux_av(:)=fold(:)
      call voigt_exp(aux_av,bux33,aux_sv,ndir,nshr,ielem_type,1)
      if (verbose) then
         write(*,*)'Deformation gradient old'
         do i=1,3
            write(*,'(3e11.3)')(bux33(i,j),j=1,3)
         enddo
      endif
      aux33(:,:)=u(:,:)
      call lu_inverse(inverr,aux33,3)  ! U^(-1)
            if (inverr.ne.0) then
               write(*,*)'** error in cacl_drot_explicit.f(2)'
               call xit()
            endif
c     ** cux33 transforms e^0 to e^(n) basis.
      call mat_dot(bux33,aux33,cux33,3)

c     ** incremental rotation tensor
c     dR = R^(new). (R^(old))^T
      call mat_transpose_fe(cux33,3,aux33)
      call mat_dot(rot_acum,aux33,drot,3)

      if (verbose) then
         write(*,*)'** drot'
         do i=1,3
            write(*,'(3e11.3)')(drot(i,j),j=1,3)
         enddo
      endif

      return
      end subroutine calc_drot_explicit
