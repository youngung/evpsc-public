      subroutine polar_dcmp(f,r,u)
c      implicit none
      real*8 f(3,3),r(3,3),u(3,3),aux33(3,3),rinvt(3,3),rnew(3,3)
      real*8 tol, err
      integer inverr,i,j,kount,k
      parameter(tol=1e-4)

      err=tol*2
      kount=1
      !! initial guess on rotation matrix.
      r(:,:)=f(:,:)
      do while(err.gt.tol)
         aux33(:,:)=r(:,:)
         call lu_inverse(inverr,aux33,3)
         call mat_transpose_fe(aux33,3,rinvt)
         rnew(:,:)=0.5*(r(:,:)+rinvt(:,:))
         aux33(:,:)=rnew(:,:)-r(:,:)

         err=0d0
         do i=1,3
            do j=1,3
               err=err+abs(aux33(i,j))
            enddo
         enddo
         r(:,:)=rnew(:,:)
c         write(*,*)'kount:',kount, 'new r with err:',err
c         write(*,'(3e24.16)')((r(i,j),j=1,3),i=1,3)
         kount=kount+1
      enddo

c     Calculate U by
c     F=RU - > U = R^(inv). F
      aux33(:,:)=r(:,:)
      call lu_inverse(inverr,aux33,3)
c     F=RU - > U = R^(T). F
c      call mat_transpose_fe(aux33,3,aux33)

      u(:,:)=0d0
      do i=1,3
         do j=1,3
            do k=1,3
               u(i,j)=u(i,j)+aux33(i,k)*F(k,j)
            enddo
         enddo
      enddo

c$$$      write(*,*)'kount:',kount-1, 'new r with err:',err
c$$$      do i=1,3
c$$$         write(*,'(3(3f8.3,3x))')(f(i,j),j=1,3),(r(i,j),j=1,3),
c$$$     $        (u(i,j),j=1,3)
c$$$      enddo
c$$$
c$$$
c$$$      aux33(:,:)=0d0
c$$$      do i=1,3
c$$$         do j=1,3
c$$$            do k=1,3
c$$$               aux33(i,j)=aux33(i,j)+r(i,k)*u(k,j)
c$$$            enddo
c$$$         enddo
c$$$      enddo
c$$$
c$$$      do i=1,3
c$$$         write(*,'(2(3f8.3,3x))')(f(i,j),j=1,3),(aux33(i,j),j=1,3)
c$$$      enddo


      return
      end subroutine polar_dcmp
