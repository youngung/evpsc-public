C **************************************************************************
C     SUBROUTINE WRITE_TEXTURE_FE      --->     VERSION 19/FEB/20
C
C     WRITE CRYSTALLOGRAPHIC TEXTURE (CRYSTAL ORIENTATION AND WEIGHT) FOR
C     EACH PHASE ASSOCIATE WITH ELEMENT
C **************************************************************************

      SUBROUTINE WRITE_TEXTURE_FE (uw,iph,ngrmx,nphmx,igb,igt,
     $     EVMX,ag,wgt,axisph0,eulerph)
      implicit none
      integer,intent(in)::uw,ngrmx,nphmx,igb,igt
      real*8,intent(in)::evmx,axisph0(0:3,3,0:nphmx),eulerph(3,0:nphmx),
     $     ag(3,3,ngrmx),wgt(ngrmx)
      real*8 AAT(3,3),gda,gdb,gdc
      integer iph,iunit,i,j,kkk

      IUNIT=uw+iph-1
      WRITE(IUNIT,'('' TEXTURE FOR PHASE'',I3,
     $''  AT VM STRAIN'',e11.3)') IPH,EVMX
c      WRITE(*,*) (AXISPH0(0,I,iph),I=1,3) ,
c     $     '  <-- axes of phase ellipsoid'
      WRITE(IUNIT,'(3e12.3,A)') (AXISPH0(0,I,iph),I=1,3) ,
     $     '  <-- axes of phase ellipsoid'
      WRITE(IUNIT,'(3e12.2,A)') (EULERPH(I,iph),I=1,3),
     $     '  <-- Euler angles of phase ellipsoid (deg)'
      WRITE(IUNIT,'(''B'',I10)') igt-igb+1
      do KKK=igb,igt!1,ngrx
         do i=1,3
            do j=1,3
               aat(i,j)=ag(j,i,KKK)
            enddo
         enddo
         call euler(1,gda,gdb,gdc,aat)
         write(IUNIT,'(3F8.2,1e12.4)') gda,gdb,gdc,wgt(KKK)
      enddo
      RETURN
      END subroutine write_texture_fe
