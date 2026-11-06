      subroutine loading_condition_fe(
c     intent(in)
     $     iregime,mdim,ntens,ielem_type,drot,dtime,
     $     dfgrd0,dfgrd1,dstran_b,
c     intent(out)
     $     strain_control,iDBARv,iSBARv,LIJBARc,DBAR,DBARc,WBARc)
c     dstran_b: dstrain in b-basis; also in e^0 co-rotating coordinate system

      implicit none
      integer,intent(in)::iregime,mdim,ntens,ielem_type
      integer,intent(out)::strain_control,iDBARv(6),iSBARv(6)
      real*8,intent(in)::drot(3,3),dtime,dfgrd0(3,3),dfgrd1(3,3),
     $     dstran_b(6)
      real*8,intent(out)::DBAR(6),DBARc(3,3),WBARc(3,3),
     $     LIJBARc(3,3)
c     locals
      integer i,j
      real*8 aux5(5),aux55(5,5),bux33(3,3),
     $     auxm(mdim),auxmm(mdim,mdim),
     $     aux6(6),aux66(6,6),aux33(3,3),
     $     aux3333(3,3,3,3),vnorm
      integer inverr

c     comparing with LIJBARc estimated from deformation gradients
      aux33(:,:)=(dfgrd1(:,:)-dfgrd0(:,:))/dtime ! Fdot
      bux33(:,:)=dfgrd1(:,:)
      call lu_inverse(inverr,bux33,3)
      if (inverr.ne.0) then
         write(*,*)'** Error in loading_condition_fe.f(1)'
         write(*,*)'   to inverse def. grad. F1'
         write(*,*)'write dtime:',dtime
         write(*,'(2a30)')'dfgrd0','dfgrd1'
         do i=1,3
            write(*,'(2(3e10.2))')
     $           (dfgrd0(i,j),j=1,3),(dfgrd1(i,j),j=1,3)
         enddo
         call xit()
      endif
c           velocity gradient tensor = fdot f^-1
      call mat_dot(aux33,bux33,LIJBARc,3)
      WBARc(:,:)=0.
      DBARc(:,:)=0.
      do i=1,3
         do j=1,3
            DBARc(i,j)=0.5*(LIJBARc(i,j)+LIJBARc(j,i))
            WBARc(i,j)=0.5*(LIJBARc(i,j)-LIJBARc(j,i))
         enddo
      enddo

      dbar(:)=0.
      if (iregime.eq.1) then
         dbarc(3,3)=-dbarc(1,1)-dbarc(2,2)
         call chg_basis(AUX5,DBARc,aux55,aux3333,2,5)
         dbar(1:5)=aux5(:)
      else
         call chg_basis(DBAR,DBARc,aux66,aux3333,2,6)
      endif

      strain_control=1
      if (ielem_type.eq.0) then !! solid elements
c        fully strain control
         iDBARv(:)=1
         iSBARv(:)=0
      elseif(ielem_type.eq.1) then !! shell elements (plane-stress)
c        mixed boundary condition
         iDBARv(:)=0
         iSBARv(:)=0


c        for EVPSC with mdim=6
         if (mdim.eq.6) then
c        e11 e12 e13    ev1 ev6 ev5     k k u
c            e22 e23        ev2 ev4       k u
c                e33            ev3         u

c        s11 s12 s13    sv1 sv6 sv5     u u k
c            s22 s23        sv2 sv4       u k
c                s33            sv3         k
            iDBARv(1)=1
            iDBARv(2)=1
            iDBARv(6)=1

            iSBARv(3)=1
            iSBARv(4)=1
            iSBARv(5)=1
         elseif(mdim.eq.5) then
c        for VPSC with mdim=5
c        e11 e12 e13    ev1 ev6 ev5     k k u
c            e22 e23        ev2 ev4       k u
c                e33            ev3         k

c        s11 s12 s13    sv1 sv6 sv5     u u k
c            s22 s23        sv2 sv4       u k
c                s33            sv3         u
            iDBARv(1)=1
            iDBARv(2)=1
            iDBARv(3)=1

            iDBARv(6)=1

            iSBARv(4)=1
            iSBARv(5)=1

c$$$            write(*,*)'dbarc diag. sum:',
c$$$     $           dbarc(1,1)+dbarc(2,2)+dbarc(3,3)
c$$$            write(*,'(2a15)')'LijBARc','dbarc'
c$$$            do i=1,3
c$$$               write(*,'(2(3e15.7,2x))')
c$$$     $              (LijBARc(i,j),j=1,3),(dbarc(i,j),j=1,3)
c$$$            enddo
c$$$            write(*,'(a10,6e10.2)')'dbar:',dbar
         endif


      endif
      return
      end subroutine loading_condition_fe
