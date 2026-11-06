C *****************************************************************************
C     SUBROUTINE RODRIGUES (ex-REORIENT_GRAIN)   --->   VERSION 03/NOV/09

C     BUILDS INCREMENTAL ROTATION MATRIX 'AROT' BASED ON RODRIGUES FORMULA.
C             --> AROT=I+sin(phi)*C+(1-cos(phi))*C^2
C     'C'   :INCREMENTAL LATTICE SPIN (ANTI-SYMMETRIC) (normalized)
C     'AROT':INCREMENTAL TRANSFORMATION FROM INITIAL TO FINAL ORIENTATION.
C *****************************************************************************

      SUBROUTINE RODRIGUES (C,AROT)
      implicit none

      real*8 c(3,3),c2(3,3),v(3),xid33(3,3),arot(3,3)
      integer i,j,k
      real*8 vnorm,coef1,coef2
      data xid33/1.,0.,0.,0.,1.,0.,0.,0.,1./

c ** v(i) is the Rodrigues spin axis. The norm is the rotation angle phi.
      v(1)=c(3,2)
      v(2)=c(1,3)
      v(3)=c(2,1)
      vnorm=sqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))
      if(vnorm.lt.1.e-06) then
cx        print *, ' zero rotation inside Rodrigues'
cx        print *, ' press c to continue'
cx        read  *
        do i=1,3
        do j=1,3
          arot(i,j)=xid33(i,j)
        enddo
        enddo
        return
      endif

      coef1=sin(vnorm)/vnorm
      coef2=(1.-cos(vnorm))/vnorm**2

      do i=1,3
      do j=1,3
        c2(i,j)=0.
        do k=1,3
          c2(i,j)=c2(i,j)+c(i,k)*c(k,j)
        enddo
      enddo
      enddo

      do i=1,3
      do j=1,3
        arot(i,j)=xid33(i,j)+coef1*c(i,j)+coef2*c2(i,j)
      enddo
      enddo

      return 
      end
