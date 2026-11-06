      subroutine write_crav(iopt,uw,NPH,
     $     ngr,nmodes,nphmx,ngrmx,nmodmx,nsm,nsysmx,wgt,crss)
      implicit none
      integer,intent(in)::iopt,uw,NPH,nphmx,ngrmx,nmodmx,
     $     nsysmx,ngr(0:nphmx),nmodes(nphmx),nsm(nmodmx,nphmx)
      real*8,intent(in)::wgt(ngrmx),crss(nsysmx,ngrmx)
c     locals
      integer iph,ism,kmo,isy,imo,kkk
      real*8 crss_wavg(nmodmx,nphmx)
      if (iopt.eq.0) then
         open(uw,file='crss_wavg.out',status='unknown')
         do iph=1,NPH
            do imo=1,nmodes(iph)
               write(uw,'(2x,a2,i2.2,a1,a3,i2.2)',advance='no')
     $              'ph',iph,'-','mod',imo
            enddo
         enddo
         write(uw,*)
      elseif(iopt.eq.1) then
         crss_wavg(:,:)=0d0
         do iph=1,NPH
         DO KKK=NGR(IPH-1)+1,NGR(IPH)
            ISM=0
         DO KMO=1,NMODES(IPH)
         DO ISY=1,NSM(KMO,IPH)
            ISM=ISM+1
            crss_wavg(kmo,iph)=crss_wavg(kmo,iph)
     $           +crss(ism,kkk)*wgt(kkk)
     $           /nsm(kmo,iph)
         enddo
         enddo
         enddo
         write(uw,'(24e12.4)')
     $        (crss_wavg(imo,iph),IMO=1,NMODES(IPH))
         enddo
      endif

      return
      end subroutine write_crav
