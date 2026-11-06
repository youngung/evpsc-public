C ***************************************************************************
C     SUBROUTINE TEXTURE_ROTATION_FE     -->      VERSION OF 17/April/2020
C
C     USER PROVIDES A ROTATION MATRIX 'ROTMAT' THAT OPERATES ON
C     CRYSTALLOGRAPIC (crystal axes) AND MORPHOLOGIC (ellipsoid
C     axes) TEXTURE AS FOLLOWS:
C
C     1- CHANGE OF REFERENCE SYSTEM   ***
C        COLUMNS OF 'ROTMAT' ARE 'old' SYSTEM AXES EXPRESSED IN 'new' SYSTEM.
C       'ROTMAT' TRANSFORMS VECTORS AND TENSORS EXPRESED IN 'old' SET OF
C        AXES, AND EXPRESSES THEIR COMPONENTS IN 'new' SET OF AXES:
C
C          Vnew(i)  =ROTMAT(i,j)*Vold(j)
C          Tnew(i,j)=ROTMAT(i,k)*ROTMAT(j,l)*Told(k,l)
C
C     2- RIGID ROTATION OF THE SAMPLE WITH RESPECT TO FIX LABORATORY AXES ***
C        'ROTMAT' ROTATES THE SAMPLE (AND SO THE ATTACHED CRYSTAL AND
C        ELLIPSOID ORIENTATIONS), FROM 'OLD' TO 'NEW' POSITION.
C        COLUMNS OF 'ROTMAT' ARE THE SAMPLE AXES AFTER ROTATION, EXPRESSED
C        IN FIX LABORATORY AXES.
C ***************************************************************************

      SUBROUTINE TEXTURE_ROTATION_fe(
c          in
     $     ROTMAT,interaction,ngr,nph,nphmx,ngrmx,ishape,ifrag,
     $     crit_shp,iflat,axisph,axisgr,
c          inout
     $     ag,eulerph,fijph,fijgr)
      implicit none
      real*8,intent(in):: ROTMAT(3,3),crit_shp(nphmx)
      integer, intent(in)::nph,ngrmx,nphmx,ngr(0:nphmx),ishape(0:nphmx),
     $     interaction,ifrag(nphmx)
      integer iflat(0:nphmx)
      real*8, intent(inout):: ag(3,3,ngrmx),fijph(3,3,0:nphmx),
     $     fijgr(3,3,ngrmx),axisph(0:3,3,0:nphmx),axisgr(0:3,3,ngrmx),
     $     eulerph(3,0:nphmx)
c     locals
      real*8 AUX33(3,3)
      integer iphel,kgx,i,j,k,igb,igt,iverbose,ierr_shp
      parameter(iverbose=0)! silent
c      parameter(iverbose=1)! verbose


      DO IPHEL=1,NPH
         DO KGX=NGR(IPHEL-1)+1,NGR(IPHEL)
C *** 'AG' TRANSFORMS FROM CRYSTAL TO 'OLD' SAMPLE AXES.
C     THE COLUMNS OF 'AG' ARE THE CRYSTAL AXES EXPRESSED IN 'OLD' SYSTEM.
C     'ROTMAT' ROTATES THEM TO 'NEW' SYSTEM.

c     ag^{crystal     <- lab}

c     new ag^{crystal <- new lab} = ag_{ij} R{j}
            DO I=1,3
            DO J=1,3
               AUX33(I,J)=0d0
               DO K=1,3
                  AUX33(I,J)=AUX33(I,J)+ROTMAT(I,K)*AG(K,J,KGX)
               ENDDO
            ENDDO
            ENDDO
            DO I=1,3
            DO J=1,3
               AG(I,J,KGX)=AUX33(I,J)
            ENDDO
            ENDDO
         ENDDO                  ! END OF LOOP OVER GRAINS IN THE PHASE

      if (iverbose.ge.1) write(*,*)'here1'

C *** COLUMNS OF 'FIJPH' ARE GRAIN AXES EXPRESSED IN 'OLD' SAMPLE SYSTEM.
C     SUBR 'UPDATE_SHAPE' CALCULATES EIGENVECTORS AND EIGENVALUES USED
C     INSIDE 'SUBR VPSC' TO ROTATE 'C4SA' INTO 'C4GA' AND CALC ESHELBY TENSOR.

         DO I=1,3
            DO J=1,3
               AUX33(I,J)=0d0
               DO K=1,3
                  AUX33(I,J)=AUX33(I,J)+ROTMAT(I,K)*FIJPH(K,J,IPHEL)
               ENDDO
            ENDDO
         ENDDO
      if (iverbose.ge.1) then
         write(*,*)'here2'
         write(*,*)
      endif

         DO I=1,3
            DO J=1,3
               FIJPH(I,J,IPHEL)=AUX33(I,J)
            ENDDO
         ENDDO

      if (iverbose.ge.1) then
         write(*,*)'here3'
         write(*,*)
      endif

         IF(ISHAPE(IPHEL).GE.1) THEN
            if (iverbose.ge.1) then
               write(*,*)'here3.1'
               write(*,*)
            endif

            DO KGX=NGR(IPHEL-1)+1,NGR(IPHEL)

               if (iverbose.ge.1) then
                  write(*,*)'here3.2'
                  write(*,*)
               endif

               DO I=1,3
                  DO J=1,3
                     AUX33(I,J)=0d0
                     DO K=1,3
                        AUX33(I,J)=AUX33(I,J)+ROTMAT(I,K)*FIJGR(K,J,KGX)
                     ENDDO
                  ENDDO
               ENDDO

               if (iverbose.ge.1) then
                  write(*,*)'here3.3'
                  write(*,*)
               endif

               DO I=1,3
                  DO J=1,3
                     if (iverbose.ge.1) then
                        write(*,*)'here3.35'
                        write(*,*)'ngrmx:',ngrmx
                     endif
                     FIJGR(I,J,KGX)=AUX33(I,J)
                     if (iverbose.ge.1) then
                        write(*,*)'here3.36'
                        write(*,*)
                     endif
                  ENDDO
               ENDDO


               if (iverbose.ge.1) then
                  write(*,*)'here3.4'
                  write(*,*)
               endif

            ENDDO
         ENDIF

      if (iverbose.ge.1) write(*,*)'here4'

         igb=ngr(iphel-1)+1
         igt=ngr(iphel)

         if (iverbose.ge.1)then
            call w_chrc(0,'calling subr. update_shape_fe')
         endif
            CALL UPDATE_SHAPE_fe(
c             in
     $           ngrmx,
     $           nphmx,
     $           iphel,           ! iph
     $           ngr(iphel-1)+1,ngr(iphel), ! igb,igt
     $           interaction,   ! iph
c     $             NPH,  ! per an integration point
     $           ifrag,    !
     $           crit_shp,
     $           ishape,   !
c     $           ngr(iphel)-ngr(iphel-1), ! number of grains
!     inout
     $           fijph,
     $           iflat,
     $           axisph,
     $           axisgr,
c     out
     $           eulerph,
     $           fijgr,ierr_shp
     $           )

c$$$         CALL UPDATE_SHAPE_fe(
c$$$c             in
c$$$     $        iphel,            ! iphel
c$$$     $        interaction,      ! iphel
c$$$c     $             NPH,  ! per an integration point
c$$$     $        ifrag(iphel),     !
c$$$     $        crit_shp(iphel),  !
c$$$     $        ishape(iphel),    !
c$$$     $        ngr(iphel)-ngr(iphel-1), ! number of grains
c$$$!     inout
c$$$     $        fijph(:,:,iphel), ! (dummy)
c$$$     $        iflat(iphel),
c$$$     $        axisph(0:3,:,iphel),
c$$$     $        axisgr(0:3,:,igb:igt),
c$$$c     out
c$$$     $        eulerph(:,iphel),
c$$$     $        fijgr(:,:,igb:igt) !  (dummy) when initializing
c$$$     $        )
c         CALL UPDATE_SHAPE (IPHEL) ! updates also grain shape if ishape>0
      ENDDO                     ! END OF LOOP OVER PHASES



      if (iverbose.ge.1) then
         call w_chrc(0,'Leaving subr. texture_rotation_fe')
      endif

      RETURN
      END subroutine texture_rotation_fe
CFEE
