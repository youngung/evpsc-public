c-----------------------------------------------------------------------
      subroutine phase_property_evp_fe(iopt,io,iph,igb,igt,ngrmx,
     $     nphmx,sg,tg,wgt,tincr,epstot_ph)
      implicit none
      integer,intent(in):: iopt,io,iph,ngrmx,igb,igt,nphmx
      real*8,intent(in):: sg(6,ngrmx),tg(6,ngrmx),wgt(ngrmx),tincr
      real*8,intent(inout)::epstot_ph(3,3,nphmx)
c     locals
      integer kkk,i,n,m
      real*8 sgph(6),tgph(6),sph33(3,3),tph33(3,3),
     $     aux66(6,6),aux3333(3,3,3,3),wph
      integer ijv(6,2)
      DATA ((IJV(N,M),M=1,2),N=1,6)/1,1,2,2,3,3,2,3,1,3,1,2/

      if (iopt.eq.0) then
         open(io,file='PH_STR'//CHAR(48+iph)//'.OUT',status='unknown')
         write(io,'(18(7x,a3,2i1))') ('sig',ijv(i,1),ijv(i,2),i=1,6),
     $        ('eps',ijv(i,1),ijv(i,2),i=1,6),
     $        ('edt',ijv(i,1),ijv(i,2),i=1,6)
      elseif(iopt.eq.1) then
         sgph(:)=0d0
         tgph(:)=0d0
         wph=0d0
         do kkk=igb,igt
            do i=1,6
               sgph(i)=sgph(i)+sg(i,kkk)*wgt(kkk)
               tgph(i)=tgph(i)+tg(i,kkk)*wgt(kkk)
            enddo
            wph=wph+wgt(kkk)
         enddo
         sgph(:)=sgph(:)/wph
         tgph(:)=tgph(:)/wph
         call chg_basis(sgph,sph33,aux66,aux3333,1,6)
         call chg_basis(tgph,tph33,aux66,aux3333,1,6)
         epstot_ph(:,:,iph)=epstot_ph(:,:,iph)+tph33(:,:)*tincr
         write(io,'(18e12.4)')(sph33(ijv(i,1),ijv(i,2)),i=1,6),
     $        (epstot_ph(ijv(i,1),ijv(i,2),iph),i=1,6),
     $        (tph33(ijv(i,1),ijv(i,2)),i=1,6)
      endif

      return
      end subroutine phase_property_evp_fe
