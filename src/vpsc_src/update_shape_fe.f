C **************************************************************************
C     SUBROUTINE UPDATE_SHAPE      --->      VERSION 08/Feb/2020
C
C     CALCULATES THE DIRECTION AND LENGTH (EIGENVECTORS AND EIGENVALUES OF
C     DEFORMATION TENSOR FIJ) OF THE AXES OF THE ELLIPSOID ASSOCIATED WITH
C     AVERAGE, PHASE AND/OR GRAIN ACCUMULATED DEFORMATION.
C     CALCULATES THE EULER ANGLES OF ELLIPSOID AXES WRT SAMPLE AXES.
C
C     REPLACES PREVIOUS SUBR. GRAXES_AVERAGE & SUBR. GRAXES_LOCAL  (SEP/04)
C
C     GRAIN-FRAGMENTATION SCHEME ADDED ON 12/NOV/07 (CNT)
C
C     SHAPE 'FREEZE' SCHEME MODIFIED ON 21/JUN/08 (CNT).
C     NOW IT APPLIES TO THE DEFINITION OF FIJ (INSTEAD OF LIJ). WORST GRAIN
C     OFFENDER CRITERION SUPPRESSED & REPLACED BY ASPECT RATIO OF THE PHASE.
C **************************************************************************
      SUBROUTINE UPDATE_SHAPE_fe(
c          in
     $     ngrmx,nphmx,IPH_x,igb,igt,interaction_x,ifrag_x,
     $     crit_shp_x,ishape_x,
c          inout
     $     fijph_x,iflat_x,axisph_x,axisgr_x,
c          out
     $     eulerph_x,fijgr_x,ier
     $     )
      implicit none
      integer, intent(in)::ngrmx,nphmx,iph_x,interaction_x,
     $     ifrag_x(nphmx),ishape_x(0:nphmx),igb,igt
      real*8,intent(in)::crit_shp_x(nphmx)
      integer, intent(inout):: iflat_x(0:nphmx)
      integer, intent(out):: ier
      real*8,intent(inout):: fijph_x(3,3,0:nphmx),
     $     axisph_x(0:3,3,0:nphmx),
     $     axisgr_x(0:3,3,ngrmx)
      real*8,intent(out)::eulerph_x(3,0:nphmx),fijgr_x(3,3,ngrmx)
c     -------------------------------------------
      real*8 W(3),BX(3,3),B(3,3),BT(3,3),FIJX(3,3)
      integer i,j,k,l,kkk,nrot,iphx,xid3(3,3)
      real*8 ang1,ang2,ang3,ratmax,ratmin,det,sign,exchange

c     locals and dummies
      real*8 aux33(3,3)
      logical verbose
c      parameter(verbose=.true.)
      parameter(verbose=.false.)

      ier=0

      xid3(:,:)=0
      do i=1,3
         xid3(i,i)=1
      enddo


C *** IF IPH_x=0 ELLIPSOID REPRESENTS AVERAGE DEFORMATION IN ELEMENT
C *** IF IPH>0 ELLIPSOID REPRESENTS AVERAGE DEFORMATION IN PHASE 'IPH'
C *** CALCULATES EIGENVALUES, EIGENVECTORS & EULER ANGLES OF ELEMENT GRAIN,
C     PHASE GRAIN, OR INDIVIDUAL GRAIN
C *** 'AXISPH_X' TRANSFORMS FROM ELLIPSOID TO SAMPLE AXES.

C *** uses LIJBARc to prevent numerical instability in SO procedure
      IPHX=IPH_x ! if iph_x=0, initialize deformation gradient.
      IF(INTERACTION_x.EQ.5) IPHX=0

c     obtain BX = (F.F^T)
      DO I=1,3
      DO J=1,3
         BX(I,J)=0d0
      DO K=1,3
         BX(I,J)=BX(I,J)+FIJPH_x(I,K,iph_x)*FIJPH_x(J,K,iph_x)
      ENDDO
      ENDDO
      ENDDO

      if (verbose) then
         aux33(:,:)=fijph_x(:,:,iph_x)
         call w_chrc(0,'fijph_x(:,:,iph_x)')
         call w_mdim(0,aux33,3,1d0)
      endif

c     Obtain eigenvectors and engenvalues of (F.F^T)
c     that are stored in B(3,3) and W(3) arrays.
      CALL JACOBI(BX,3,3,W,B,NROT,IER)
      CALL EIGSRT(W,B,3,3)
      if (verbose)then
         call w_chrc(0,'bx:')
         call w_mdim(0,bx,3,1d0)
         call w_chrc(0,'w:')
         call w_dim(0,w,3,1d0,.true.)
      endif

c$$$      IF (IER.EQ.1) THEN
c$$$         call w_chrc(0,'bx:')
c$$$         do i=1,3
c$$$            write(*,'(3e14.6)')(bx(i,j),j=1,3)
c$$$         enddo
c$$$
c$$$         call w_ival(0,'iph_x:',iph_x)
c$$$         call w_chrc(0,'FIJPH_x:')
c$$$         do i=1,3
c$$$            write(*,'(3e14.6)')(FIJPH_x(I,K,iph_x),k=1,3)
c$$$         enddo
c$$$
c$$$         WRITE(*,*)'** ERR IN UPDATE_SHAPE FOR PHASE ELLIPSOID',
c$$$     $        IPH_x
c$$$c         call stop_debug(0,'error in update_shape')
c$$$c         call xit()
c$$$      ENDIF

C *** EIGENVALUES (AND ASSOC EIGENVECTORS) ARE ORDERED FROM LARGER TO SMALLER.
C *** REDEFINE AXIS(2) TO BE THE LARGEST IN ORDER TO IMPROVE ACCURACY IN THE
C     CALCULATION OF THE ESHELBY TENSOR.
C *** IF DET(B)<0 MEANS THAT THE SYSTEM IS LEFT HANDED. IT IS MADE RIGHT
C     HANDED BY EXCHANGING 1 AND 2.

      SIGN=-1.
      IF(DET(B).LE.0.) SIGN=1.
      DO I=1,3
         EXCHANGE=B(I,1)
         B(I,1)=B(I,2)
         B(I,2)=EXCHANGE*SIGN
      ENDDO
      EXCHANGE=W(1)
      W(1)=W(2)
      W(2)=EXCHANGE

      RATMAX=SQRT(W(2)/W(3))
      RATMIN=SQRT(W(1)/W(3))

      DO I=1,3
         IF(IFLAT_x(iph_x).EQ.0) AXISPH_X(0,I,iph_x)=SQRT(W(I))
c         keep previous
         IF(IFLAT_x(iph_x).EQ.1) AXISPH_X(0,I,iph_x)=AXISPH_X(0,I,iph_x)
         DO J=1,3
            AXISPH_X(I,J,iph_x)=B(I,J)
            BT(I,J)        =B(J,I)
         ENDDO
      ENDDO

c *** recalculates deform grad in sample keeping ellipsoid axes invariant
      IF(IFLAT_x(iph_x).EQ.1) then
         DO I=1,3
         DO J=1,3
            FIJX(I,J)=XID3(I,J)*AXISPH_X(0,I,iph_x)
         ENDDO
         ENDDO
         DO I=1,3
         DO J=1,3
            FIJPH_x(I,J,iph_x)=0d0
         DO K=1,3
         DO L=1,3
            FIJPH_X(I,J,iph_x)=FIJPH_X(I,J,iph_x)
     $           +B(I,K)*B(J,L)*FIJX(K,L)
         ENDDO
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      CALL EULER(1,ANG1,ANG2,ANG3,BT)
      EULERPH_x(1,iph_x)=ANG1
      EULERPH_x(2,iph_x)=ANG2
      EULERPH_x(3,iph_x)=ANG3


      IF (IPH_x.NE.0) THEN
c     run except for iph=0 (the initialization)

         IF (IFRAG_x(iph_x).EQ.0) THEN
C *************************************************************************
C *** STOPS UPDATING THE ELLIPSOID SHAPE FOR A LIMIT ASPECT RATIO

         IF(RATMAX.GT.crit_shp_x(iph_x) .AND. IFLAT_x(iph_x).EQ.0) THEN
CFEB
            aux33(:,:)=fijph_x(:,:,iph_x)
            call w_chrc(0,'fijph_x(:,:,iph_x)')
            call w_mdim(0,aux33,3,1d0)

            write(*,*)'IPH_x:',IPH_x
            write(*,*)'RATMAX:',RATMAX
            write(*,*)'CRIT_SHP_X',CRIT_SHP_X(iph_x)
            write(*,*)'IFLAT_x:',IFLAT_x(iph_x)
            write(*, *) '*********************************************'
            write(*, *) '  STOPS UPDATING GRAIN SHAPE FOR PHASE ',iph_x
            write(*, *) '*********************************************'
C     FEE
            IFLAT_X(iph_x)=1
         ENDIF
C *************************************************************************

         ELSE IF (IFRAG_x(iph_x).EQ.1) THEN

C *************************************************************************
C *** SPLITS AVERAGE GRAIN WHEN IT EXCEEDS A CRITICAL ASPECT RATIO.
C *** REDEFINES LENGTH OF ELLIPSOID AXES 'AXISPH_X(0,I)' BUT NOT THE
C     ELLIPSOID ORIENTATION 'AXIS(J,I,IPH_X)'
C *** RECALCULATES DEFORMATION GRADIENT OF THE PHASE 'FIJPH_X' IN SAMPLE AXES

         IF(RATMAX.GE.CRIT_SHP_X(iph_x)) THEN
            W(2)=W(2)/4d0
            IF(RATMIN.GE.CRIT_SHP_X(iph_x)/2.) W(1)=W(1)/4.

ccc       w(1)=1.    ! resetting to equiaxed
ccc       w(2)=1.
ccc       w(3)=1.

            DO I=1,3
               AXISPH_X(0,I,iph_x)=SQRT(W(I))
            DO J=1,3
               FIJX(I,J)=XID3(I,J)*AXISPH_X(0,I,iph_x)
            ENDDO
            ENDDO
            DO I=1,3
            DO J=1,3
               FIJPH_X(I,J,iph_x)=0.
            DO K=1,3
            DO L=1,3
               FIJPH_X(I,J,iph_x)=FIJPH_X(I,J,iph_x)
     $              +B(I,K)*B(J,L)*FIJX(K,L)
            ENDDO
            ENDDO
            ENDDO
            ENDDO
c            END OF IF RATMAX>CRIT_RAT
         ENDIF
C ********************************************************************
c     END OF IF IFRAG_x
         ENDIF

c     END OF IF IPH_X.NE.0 (run excpet for initial)
      ENDIF
CFEB
      if (verbose) then
         call w_ival(0,'iph_x:',iph_x)
         call w_chrc(0,'fijph_x(1:3,1:3,iph_x)')
         call w_mdim(0,fijph_x(1:3,1:3,iph_x),3,1d0)
         call w_chrc(0,'EIGNVAL')
         call w_dim(0,AXISPH_X(0,:,iph_x),3,1d0,.true.)
         call w_chrc(0,'EIGNVEC')
         call w_mdim(0,axisph_x(1:3,1:3,iph_x),3,1d0)
      endif

C ********************************************************************
C *** UPDATES ELLIPSOIDS (EIGENVALUES & EIGENVECTORS) OF INDIVIDUAL GRAINS
C     'AXISGR_X(I,J)' TRANSFORMS FROM ELLIPSOID TO SAMPLE AXES.

      IF(ISHAPE_x(iph_x).GT.0) THEN
C ********************************************************************

      DO KKK=igb,igt
         DO I=1,3
         DO J=1,3
            BX(I,J)=0d0
         DO K=1,3
            BX(I,J)=BX(I,J)+FIJGR_X(I,K,KKK)*FIJGR_X(J,K,KKK)
         ENDDO
         ENDDO
         ENDDO

         CALL JACOBI(BX,3,3,W,B,NROT,IER)
         CALL EIGSRT(W,B,3,3)
c$$$         IF (IER.EQ.1) THEN
c$$$            WRITE(*,'(''ERROR IN UPDATE_SHAPE FOR GRAIN'',I5,
c$$$     $''  IN PHASE'',I3)') KKK,IPH_X
c$$$c            call xit()
c$$$         ENDIF

C *** EIGENVALUES (AND ASSOC EIGENVECTORS) ARE ORDERED FROM LARGER TO SMALLER.
C *** REDEFINE AXIS(2) TO BE THE LARGEST IN ORDER TO IMPROVE ACCURACY IN THE
C     CALCULATION OF THE ESHELBY TENSOR.
C *** IF DET(B)<0 MEANS THAT THE SYSTEM IS LEFT HANDED. IT IS MADE RIGHT
C     HANDED BY EXCHANGING 1 AND 2.

         SIGN=-1.
         IF(DET(B).LE.0.) SIGN=1.
         DO I=1,3
            EXCHANGE=B(I,1)
            B(I,1)=B(I,2)
            B(I,2)=EXCHANGE*SIGN
         ENDDO
         EXCHANGE=W(1)
         W(1)=W(2)
         W(2)=EXCHANGE

         RATMAX=SQRT(W(2)/W(3))
         RATMIN=SQRT(W(1)/W(3))

         DO I=1,3
            IF(IFLAT_X(iph_x).EQ.0) AXISGR_X(0,I,KKK)=SQRT(W(I))
            IF(IFLAT_X(iph_x).EQ.1) AXISGR_X(0,I,KKK)=AXISGR_X(0,I,KKK)
            DO J=1,3
               AXISGR_X(I,J,KKK)=B(I,J)
            ENDDO
         ENDDO

C *** updates def gradient keeping ellipsoid axes invariant   !cnt june/11/08
         IF(IFLAT_X(iph_x).EQ.1) THEN
            DO I=1,3
            DO J=1,3
               FIJX(I,J)=XID3(I,J)*AXISGR_X(0,I,KKK)
            ENDDO
            ENDDO
            DO I=1,3
            DO J=1,3
               FIJGR_X(I,J,KKK)=0.
            DO K=1,3
            DO L=1,3
               FIJGR_X(I,J,KKK)=FIJGR_X(I,J,KKK)+B(I,K)*B(J,L)*FIJX(K,L)
            ENDDO
            ENDDO
            ENDDO
            ENDDO
         ENDIF

         IF (IFRAG_X(iph_x).EQ.1) THEN

C ********************************************************************
C *** SPLITS IN HALF AXES OF THE GRAIN WHEN LIMIT ASPECT RATIOS ARE REACHED

            IF(RATMAX.GE.CRIT_SHP_X(iph_x)) THEN
               W(2)=W(2)/4.
               IF(RATMIN.GE.CRIT_SHP_X(iph_x)/2.) W(1)=W(1)/4.

               DO I=1,3
                  AXISGR_X(0,I,KKK)=SQRT(W(I))
               DO J=1,3
                  FIJX(I,J)=XID3(I,J)*AXISGR_X(0,I,KKK)
               ENDDO
               ENDDO
               DO I=1,3
               DO J=1,3
                  FIJGR_X(I,J,KKK)=0.
               DO K=1,3
               DO L=1,3
                  FIJGR_X(I,J,KKK)=FIJGR_X(I,J,KKK)+
     $                 B(I,K)*B(J,L)*FIJX(K,L)
               ENDDO
               ENDDO
               ENDDO
               ENDDO
            ENDIF               ! END OF IF RATMAX
C ********************************************************************
         ENDIF                  ! END OF IF IFRAG

      ENDDO                     ! END OF DO KKK

C ********************************************************************
      ENDIF                     ! END OF IF ISHAPE_x.GT.0
C ********************************************************************
CFEE
      RETURN
      END subroutine update_shape_fe
