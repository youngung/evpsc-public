c======================================================================c
      subroutine rexdd_ori_energy(ag_pa,cel,eps_sa,ag_ch)
      implicit none
c     ag_pa: pa<-sa
      real*8, intent(in):: ag_pa(3,3),cel(6,6),eps_sa(6)
      real*8, intent(out):: ag_ch(3,3)
c     locals
      real*8 eps_ca(6),aux6(6),aux33(3,3),bux33(3,3),aux66(6,6),
     $     aux3333(3,3,3,3),x(3),x_ori(3),dx(3)
      integer i,j,k,l
c     numerical recipe parameters
      real*8 h,dx_ref,tol,df,df_frac_tol,f0,grad(3),hess(3,3),f0new,mag
      integer kount,kmax,isingular
c     controls
      logical verbose
      parameter(kmax=100)
c      parameter(verbose=.true.)
      parameter(verbose=.false.)

c     transform strain in sample axes to crystal axes of parent.
      call chg_basis(eps_sa,aux33,aux66,aux3333,1,6)
      bux33(:,:)=0.
      do i=1,3
      do j=1,3
      do k=1,3
      do l=1,3
         bux33(i,j)=bux33(i,j)+ag_pa(k,i)*ag_pa(l,j)*aux33(k,l)
      enddo
      enddo
      enddo
      enddo
      call chg_basis(eps_ca,bux33,aux66,aux3333,2,6)

c     initial guess confined to an Euler space.
c      call rexdd_gen_rand_eul(90.,90.,90.,x(1),x(2),x(3)) !! x is the Euler angle vector: (phi1,phi,phi2)
      call rexdd_gen_rand_eul(360.,360.,360.,x(1),x(2),x(3)) !! x is the Euler angle vector: (phi1,phi,phi2)
      x_ori(:)=x(:)

c     Gradient descend begins here.
c     Numerical recipe
      h=1e-2                    ! degree
      dx_ref=3.                 ! degree
      kount=0
      mag=0.
      df_frac_tol=1e-6

      call phi_el_ener(x_ori(1),x_ori(2),x_ori(3),cel,eps_ca,f0)
      df=f0
c
c     YJ @ April 12, 2026:  Maybe, hess is not necessary
      call grad_and_hess(x(1),x(2),x(3),h,cel,eps_ca,grad,hess)

      if (verbose) then
         write(*,'(a3,a11,  3a8,  x,3a8, x,3a8,  a9,  a11)')
     $        'k','mag','dx1','dx2','dx3','x1','x2','x3',
     $        'grad1','grad2','grad3','f0','df'
         write(*,'(i3,e11.3,3f8.2,x,3f8.2,x,3f8.2,f9.4,e11.3)')
     $        kount,mag,(dx(i),i=1,3),(x(i),i=1,3),(grad(i),i=1,3),f0,df
      endif

      do while(kount.le.kmax .and. abs(df).gt.abs(f0*df_frac_tol))
         kount=kount+1
c$$$         ! solve H.dx=grad
c$$$         if (.false.) then
c$$$            call lu_eqsystem(Hess,grad,3,isingular)
c$$$            if (isingular.eq.1) then
c$$$               write(*,*)'** Error: Singularity in solving LU system'
c$$$               write(*,*)'H matrix'
c$$$               do i=1,3
c$$$                  write(*,'(3e11.3,x,a1,i1,x,e11.3)')
c$$$     $                 Hess(i,1),Hess(i,2),Hess(i,3),'x',
c$$$     $                 i,grad(i)
c$$$               enddo
c$$$               call xit()
c$$$            endif
c$$$         endif
         mag=sqrt(grad(1)**2+grad(2)**2+grad(3)**2)
         do i=1,3
            dx(i)=grad(i)/mag*dx_ref
         enddo
         x(:)=x(:)-dx(:)
         call phi_el_ener(x(1),x(2),x(3),cel,eps_ca,f0new)

         df=f0new-f0
         f0=f0new
         if (verbose) then
            write(*,'(i3,e11.3,3f8.2,x,3f8.2,x,3f8.2,f9.4,e11.3)')
     $           kount,mag,(dx(i),i=1,3),(x(i),i=1,3),(grad(i),i=1,3),
     $           f0,df
c            write(*,'(i3,e11.3,3f8.2,x,3f8.2,f9.4,e11.3)')kount,mag,
c     $           (dx(i),i=1,3),(x(i),i=1,3),f0,df
         endif
         call grad_and_hess(x(1),x(2),x(3),h,cel,eps_ca,grad,hess)
      enddo
      call euler(2,x(1),x(2),x(3),aux33) !! ch<-pa
      ag_ch(:,:)=0.
      do i=1,3
      do j=1,3
         do k=1,3
            ag_ch(i,j)=ag_ch(i,j)+aux33(i,k)*ag_pa(k,j)
         enddo
      enddo
      enddo
      return
      end subroutine
c======================================================================c
      subroutine phi_el_ener(phi1,phi,phi2,cel,eps,energy)
      implicit none
      real*8,intent(in):: phi1,phi,phi2,cel(6,6),eps(6)
      real*8,intent(out)::energy
c     locals
      integer i,j,k,l
      real*8 R(3,3),aux33(3,3),bux33(3,3),aux66(6,6),
     $     aux3333(3,3,3,3),E(6)

c     cartesian components
      call chg_basis(eps,aux33,aux66,aux3333,1,6)

c     euler angle to transformation matrix.
      call euler(2,phi1,phi,phi2,R)
c     transform
      bux33(:,:)=0.
      do i=1,3
      do j=1,3
         do k=1,3
         do l=1,3
            bux33(i,j)=bux33(i,j)+R(i,k)*R(j,l)*aux33(k,l)
         enddo
         enddo
      enddo
      enddo
c     back to b-basis
      call chg_basis(E,bux33,aux66,aux3333,2,6)

      energy=0.
      do i=1,6
      do j=1,6
         energy=energy+cel(i,j)*E(i)*E(j)
      enddo
      enddo
      energy=energy*0.5
      return
      end subroutine
c======================================================================c
      subroutine grad_and_hess(x,y,z,h,cel,eps,grad,hess)
c     cel: elastic modulus in crystal axes
c     eps: strain tensor of the parent in crystal axes.
      implicit none
      real*8, intent(in) :: x,y,z,h,cel(6,6),eps(6)
      real*8, intent(out):: grad(3),hess(3,3)
c     locals
      real*8 h2inv,hsqinv
      real*8 f0,fxp,fxm,fyp,fym,fzp,fzm,fxy,fxz,fyz
c     --
      call phi_el_ener(x  ,y  ,z  ,cel,eps,f0)
      call phi_el_ener(x+h,y  ,z  ,cel,eps,fxp)
      call phi_el_ener(x-h,y  ,z  ,cel,eps,fxm)
      call phi_el_ener(x  ,y+h,z  ,cel,eps,fyp)
      call phi_el_ener(x  ,y-h,z  ,cel,eps,fym)
      call phi_el_ener(x  ,y  ,z+h,cel,eps,fzp)
      call phi_el_ener(x  ,y  ,z-h,cel,eps,fzm)
c     --
      call phi_el_ener(x+h,y+h,z  ,cel,eps,fxy)
      call phi_el_ener(x+h,y  ,z+h,cel,eps,fxz)
      call phi_el_ener(x  ,y+h,z+h,cel,eps,fyz)
c     --
c     Gradient
      h2inv=1./(2.*h)
      grad(1)=(fxp-fxm)*h2inv
      grad(2)=(fyp-fym)*h2inv
      grad(3)=(fzp-fzm)*h2inv
c     Hessian matrix
      hsqinv=1./(h**2)
c     normals
      hess(1,1)=(fxp+2*f0-fxm)*hsqinv
      hess(2,2)=(fyp+2*f0-fym)*hsqinv
      hess(3,3)=(fzp+2*f0-fzm)*hsqinv
c     shears
      hess(1,2)=(fxy-fxp-fyp+f0)*hsqinv
      hess(2,1)=hess(1,2)
      hess(1,3)=(fxz-fxp-fzp+f0)*hsqinv
      hess(3,1)=hess(1,3)
      hess(2,3)=(fyz-fyp-fzp+f0)*hsqinv
      hess(3,2)=hess(2,3)
      return
      end subroutine
