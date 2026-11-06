C **************************************************************************
C     SUBROUTINE UPDATE_SCHMID_FE
C
C     ROTATES SCHMID TENSORS OF EACH GRAIN FROM CRYSTAL TO SAMPLE AXES
C **************************************************************************

      SUBROUTINE UPDATE_SCHMID_fe(
c     in
     $     iph,igb,igt,nphmx,ngrmx,nsysmx,nsyst,schca,nschca,ag,
c     out
     $     sch,nsch)
      implicit none
      integer,intent(in)::iph,nphmx,ngrmx,nsysmx,igb,igt,
     $     nsyst(nphmx)
      real*8,intent(in)::schca(5,nsysmx,nphmx),nschca(6,nsysmx,nphmx),
     $     ag(3,3,ngrmx)
      real*8,intent(out)::sch(5,nsysmx,ngrmx),nsch(6,nsysmx,ngrmx)

c     local
      real*8 aux5(5),aux33(3,3),aux55(5,5),aux3333(3,3,3,3),aux33r(3,3),
     $       bux6(6),bux33(3,3),bux33r(3,3)
      integer i1,j1,i,j,kkk,is

      do kkk=igb,igt
         do is=1,nsyst(iph)
            DO J=1,5
               AUX5(J)=SCHCA(J,IS,iph)
            ENDDO
            DO J=1,6
               BUX6(J)=NSCHCA(J,IS,iph)
            ENDDO
            ! write(*,'(a6,5f5.2)')'schca',aux5
            CALL CHG_BASIS(AUX5,AUX33,AUX55,AUX3333,1,5)
            CALL CHG_BASIS(BUX6,BUX33,AUX55,AUX3333,1,6)
            DO I=1,3
               DO J=1,3
                  AUX33R(I,J)=0d0
                  BUX33R(I,J)=0d0
                  DO I1=1,3
                     DO J1=1,3
                        AUX33R(I,J)=AUX33R(I,J)
     $                       +AG(I,I1,KKK)*AG(J,J1,KKK)*AUX33(I1,J1)
                        BUX33R(I,J)=BUX33R(I,J)
     $                       +AG(I,I1,KKK)*AG(J,J1,KKK)*BUX33(I1,J1)
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
            CALL CHG_BASIS(AUX5,AUX33R,AUX55,AUX3333,2,5)
            CALL CHG_BASIS(BUX6,BUX33R,AUX55,AUX3333,2,6)
            DO J=1,5
               SCH(J,IS,kkk)=AUX5(J)
            ENDDO
            DO J=1,6
               NSCH(J,IS,kkk)=BUX6(J)
            ENDDO
         ENDDO
      !write(*,*)'NSCH(1,1,',kkk,'):',NSCH(1,1,kkk)
      ENDDO

      RETURN
      END subroutine update_schmid_fe

c     ---------- below subroutine is to be separated once completed.
c     iopt=0: initialize the column,
c     iopt=1: writes.
      subroutine stat_non_schmid(iopt,uw,nphmx,ngrmx,nsysmx,
     $     ngr,nsyst,dnca,dbca,ag,wgt)
      implicit none
      integer, intent(in):: iopt,uw,nphmx,ngrmx,nsysmx,ngr(0:nphmx),
     $     nsyst(nphmx)
      real*8, intent(in)::dnca(3,nsysmx,nphmx),dbca(3,nsysmx,nphmx),
     $     ag(3,3,ngrmx),wgt(ngrmx)
c     locals
      integer i,is,kkk,iph
      real*8 n(3),b(3),t(3),r(3,3),sn(3),sb(3),st(3)
      parameter(iph=1)
c      parameter(is=1)
      if (iopt.eq.0) then
         write(*,*)'entering iopt 0'
         open(uw,file='stat_nonschmid.out',status='unknown')
         write(uw,'(30a11)')
     $        'n1','n2','n3',
     $        'b1','b2','b3',
     $        't1','t2','t3',
     $        'n1n1','n2n2','n3n3',
     $        'b1b1','b2b2','b3b3',
     $        't1t1','t2t2','t3t3',


c              nb correlation
     $        'n1b1',
     $        'n2b2',
     $        'n3b3',
c              nt correlation
     $        'n1t1',
     $        'n2t2',
     $        'n3t3',
c              bt correlation
     $        'b1t1',
     $        'b2t2',
     $        'b3t3',
c
     $        'n3n3+b3b3+t3t3',
c
     $        'n3b3+n3t3+b3t3',

     $        'wgt'
         return
      elseif(iopt.eq.1) then
         write(*,*)'entering iopt 1, ngr:',ngr
         write(*,*)'iph:',iph
         do is=1,nsyst(iph)
            n(:)=dnca(:,is,iph)
            b(:)=dbca(:,is,iph)
c        ag transforms from ca to sa
            call crossproduct(n,b,t)
            do kkk=ngr(iph-1)+1,ngr(iph)
               r(:,:)=ag(:,:,kkk)
               call matvect_dot(r,n,sn,3)
               call matvect_dot(r,b,sb,3)
               call matvect_dot(r,t,st,3)
               write(uw,'(30e11.3)')
     $              (sn(i),i=1,3),
     $              (sb(i),i=1,3),
     $              (st(i),i=1,3),

     $              (sn(i)**2,i=1,3),
     $              (sb(i)**2,i=1,3),
     $              (st(i)**2,i=1,3),

c                   n x b
     $              (sn(i)*sb(i),i=1,3),
c                   n x t
     $              (sn(i)*st(i),i=1,3),
c                   b x t
     $              (sb(i)*st(i),i=1,3),

c               'n3n3+b3b3+t3t3=1'
     $              (n(3)**2+b(3)**2+t(3)**2),
c               'n3b3+n3t3+b3t3',
     $              n(3)*b(3)+n(3)*t(3)+b(3)*t(3),
     $              wgt(kkk)
            enddo
         enddo
      else
         write(*,*)'** Error: Unexpected iopt to <stat_non_schmid>'
      endif
      return
      end subroutine stat_non_schmid
