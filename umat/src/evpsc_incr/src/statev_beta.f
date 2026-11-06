c     calculate the dimension of beta array
      subroutine calc_nbeta(nbeta,ngrmx,nsysmx)
      implicit none
      integer, intent(out)::nbeta
      integer, intent(in)::ngrmx,nsysmx
c           gtotgr,   crss,      ag, LIJBARc,    sch
      nbeta=ngrmx+nsysmx*ngrmx+9*ngrmx+9+  5*nsysmx*ngrmx
      return
      end subroutine calc_nbeta

c     iopt(1): k arrays stored to beta(1,:)
c     iopt(2): k arrays stored to beta(2,:)
      subroutine statev_beta(nbeta,beta,ngrmx,nsysmx,nphmx,
     $     nph_x,ngr,
     $     LIJBARc,
     $     crss,gtotgr,ag,sch,!fijph
     $     iopt)
      implicit none
      integer,intent(in)::nbeta,iopt,ngrmx,nsysmx,nphmx,nph_x,
     $     ngr(0:nphmx)
      real*8,intent(inout):: beta(2,nbeta),crss(nsysmx,ngrmx),
     $     gtotgr(ngrmx),ag(3,3,ngrmx),sch(5,nsysmx,ngrmx),LIJBARc(3,3)
c     $     fijph(3,3,0:nphmx)
c     locals
      integer iphel,i,j,kkk,isys,k,igb,igt

      if(iopt.eq.1.or.iopt.eq.2) then
	 k=0
	 do iphel=1,nph_x
	    igb=ngr(iphel-1)+1
	    igt=ngr(iphel)
	    do kkk=igb,igt
	       do isys=1,nsysmx
		  k=k+1
		  beta(iopt,k)=crss(isys,kkk)
	       enddo
	       k=k+1
	       beta(iopt,k)=gtotgr(kkk)
	       do i=1,3
		  do j=1,3
		     k=k+1
		     beta(iopt,k)=ag(i,j,kkk)
		  enddo
	       enddo
	       do i=1,5
		  do j=1,nsysmx
		     k=k+1
		     beta(iopt,k)=sch(i,j,kkk)
		  enddo
	       enddo
	    enddo               ! kkk
	 enddo
	 do i=1,3
	    do j=1,3
	       k=k+1
	       beta(iopt,k)=LIJBARc(i,j)
	    enddo
	 enddo
      else
	 write(*,*)'unexpected case in statev_beta'
	 call xit()
      endif

      return
      end subroutine statev_beta
