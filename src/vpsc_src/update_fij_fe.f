C *****************************************************************************
C     SUBROUTINE UPDATE_FIJ_fe      --->      VERSION 12/Mar/2020
C
C     USES THE VELOCITY GRADIENT (AVERAGE, PHASE or GRAIN) IN THE STEP
C     TO UPDATE INCREMENTALLY THE CORRESPONDING DEFORMATION TENSOR 'FIJ'
C
C     REPLACES PREVIOUS SUBR. UPDFIJ_AVERAGE & SUBR. UPDFIJ_LOCAL  (SEP/04)
C *****************************************************************************

      SUBROUTINE UPDATE_FIJ_fe(
c          in
     $     nphmx,ngrmx,iph,igb,igt,
     $     dtime,LIJBARc_or_xlijph,ishape,xlijgr,
c          inout
     $     fijph,fijgr)
      implicit none
      integer,intent(in)::nphmx,ngrmx,iph,igb,igt,ishape(0:nphmx)
      real*8,intent(in)::dtime,LIJBARc_or_xlijph(3,3),xlijgr(3,3,ngrmx)
      real*8,intent(inout)::fijph(3,3,0:nphmx),fijgr(3,3,ngrmx)
c     locals
      real*8 FNEW(3,3),xid3(3,3)
      integer i,j,k,kkk
      xid3(:,:)=0d0
      do i=1,3
         xid3(i,i)=1d0
      enddo
C *** UPDATES THE DEFORM GRAD IN THE ELEMENT 'FIJPH(i,j,0)' USING THE
C     MACROSCOPIC VELOCITY GRADIENT 'LIJBARc'
C *** UPDATES THE DEFORM GRAD IN EACH PHASE 'FIJPH(i,j,IPH)' USING THE
C     AVERAGE VELOCITY GRADIENT FOR THE PHASE 'LIJPH' CALCULATED INSIDE
C     SUBROUTINE UPDATE_ORIENTATION.
C *** LIJPH ACCOUNTS FOR ROTATIONS BUT NOT FOR STRETCH WHEN IFLAT(iph)=1.
C *** FIJPH COINCIDES WITH FIJ OF ELEMENT IF NPH=1 AND IFLAT(1)=0.
c     LIJBARc or xlijph
      FNEW(:,:)=0d0
      DO I=1,3
      DO J=1,3
c        IF(IPH.EQ.0) THEN
          DO K=1,3
            FNEW(I,J)=FNEW(I,J)+(dtime*LIJBARc_or_xlijph(I,K)+XID3(I,K))
     $            *FIJPH(K,J,iph)
          ENDDO
c        ELSE IF(IPH.GT.0) THEN
c          DO K=1,3
c            FNEW(I,J)=FNEW(I,J)+(dtime*XLIJPH(I,K)+XID3(I,K))
c     $            *FIJPH(K,J)
c          ENDDO
c        ENDIF
      ENDDO
      ENDDO

      FIJPH(:,:,iph)=FNEW(:,:)
C *** UPDATES THE DEFORM GRAD IN EACH GRAIN 'FIJGR(i,j,KKK)' USING THE
C     VELOCITY GRADIENT FOR THE GRAIN 'LIJGR' CALCULATED INSIDE SUBROUTINE
C     UPDATE_ORIENTATION.
      IF (ISHAPE(iph).GT.0) THEN
         DO KKK=igb,igt
            FNEW(:,:)=0d0
            DO I=1,3
            DO J=1,3
            DO K=1,3
               FNEW(I,J)=FNEW(I,J)+(dtime*XLIJGR(I,K,KKK)+XID3(I,K))
     $              *FIJGR(K,J,KKK)
            ENDDO
            ENDDO
            ENDDO
            FIJGR(:,:,KKK)=FNEW(:,:)
         ENDDO                  ! END OF DO KKK
      ENDIF
      RETURN
      END subroutine update_fij_fe
