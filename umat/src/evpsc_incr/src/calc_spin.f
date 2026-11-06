c     following Eqs. 37, 38 of Zecevic et al. IJP (2017)
      subroutine calc_spin(drot,dtime,WBARc)
      implicit none
      real*8,intent(in)::drot(3,3),dtime
      real*8, intent(out)::WBARc(3,3)

c     locals
      real*8 omega(3),da,n(3),mag,sum
      integer i,j

      sum = drot(1,1)+drot(2,2)+drot(3,3)
      sum=0.5d0*(sum-1d0)

c      write(*,*)'sum:',sum

      if (abs(sum-1d0).le.1e-10) then
         da=0d0
      else
         da=dacos(sum)
      endif

c      write(*,*)'sum:',sum
c      write(*,*)'da:',da

      if (dabs(da).lt.1e-10) THEN
         WBARc(:,:)=0d0
         return
      endif

      n(1)=-drot(2,3)+drot(3,2)
      n(2)= drot(1,3)-drot(3,1)
      n(3)=-drot(1,2)+drot(2,1)
c     normalize n
      mag=sqrt(n(1)**2+n(2)**2+n(3)**2)

c      write(*,*)'mag:',mag

      do i=1,3
         n(i)=n(i)/mag
      enddo
c      write(*,*)'n'
c      write(*,'(3e12.4)')(n(j),j=1,3)

c      write(*,*)'drot'
c      do i=1,3
c         write(*,'(3e12.4)')(drot(i,j),j=1,3)
c      enddo

c      write(*,*)'dtime:',dtime

c     Define WBARc
      do i=1,3
         omega(i)=da/dtime * n(i)
         WBARc(i,i)=0d0
      enddo

c      write(*,'(a10,3e12.4)')'omega:',(omega(i),i=1,3)

      WBARc(1,2)=-omega(3)
      WBARc(2,1)= omega(3)
      WBARc(1,3)= omega(2)
      WBARc(3,1)=-omega(2)
      WBARc(2,3)=-omega(1)
      WBARc(3,2)= omega(1)

      return
      end subroutine calc_spin
