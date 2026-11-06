c     ------------------------------------------------------------------
c     defining uvarm_custom to utilize user variable output in Abaqus CAE
c     diffraction

      if (noel.gt.nelmx) then
         write(*,*)'-----------------------------------------'
         write(*,*)'*** Error: Noel greater than nelmx found'
         write(*,*)'noel:',noel
         write(*,*)'nelmx:',nelmx
         call xit()
      endif



c$$$      if (.false.) then
c$$$c         do i=1,2               !limits to only first three planes
c$$$c            uvarm_custom(noel,npt,i)=eps_w(i,1)
c$$$c         end do
c$$$         do i=1,2               !limits to only first three planes
c$$$            uvarm_custom(noel,npt,i+2)=wgtsetph(i,1)
c$$$         end do
c$$$c     In case one wants to look at dislocation density from RGVB
c$$$c      uvarm_custom(noel,npt,7)=dum1
c$$$c      uvarm_custom(noel,npt,8)=dum2
c$$$c      uvarm_custom(noel,npt,9)=dum3 !dislocation density; valid only if rgvb model is used.
c$$$
c$$$
c$$$c     jacobian, 1/compliance1111, DBAR011
c$$$
c$$$         call chg_basis(aux6,aux33,mbar,aux3333,3,6) ! compliance
c$$$         uvarm_custom(noel,npt,6)=1d0/(aux3333(1,1,1,1)*dtime)
c$$$         uvarm_custom(noel,npt,7)=1d0/aux3333(1,1,1,1) ! young's modulus
c$$$         call chg_basis(DBAR0,aux33,aux66,aux3333,1,6) ! DBAR011
c$$$         uvarm_custom(noel,npt,8)=aux33(1,1)
c$$$         if (kinc.eq.1.and.jstep(1).eq.1) uvarm_custom(noel,npt,9)=0d0
c$$$         uvarm_custom(noel,npt,9)=dth_tot

      if(.true.) then
c        store velocity gradient to uvarm_custom
c        LIJBARc stored in uvarm is used to create history for stand-alone
c        validation via <gen_hist.py> used in <lookup.py>
         k=1
         do i=1,3
         do j=1,3
            uvarm_custom(noel,npt,k)=LIJBARc(i,j) !1-9
            k=k+1
         enddo
         enddo
         uvarm_custom(noel,npt,k)=epsvm!10
         k=k+1
         uvarm_custom(noel,npt,k)=svm !11
         k=k+1

         if (ihardlaw.eq.20) then

            if (NPH.ne.1)then
               write(*,*)'** uvarm_custom for ihardlaw of 20'
               write(*,*)'   expects only single phase cases'
               call xit()
            endif
            iphel=1

            do kmo=1,nslmod(iphel)
               uvarm_custom(noel,npt,k)=rho_for_ph_mod(kmo,iphel)
               k=k+1
            enddo
            uvarm_custom(noel,npt,k)=rho_deb_ph(iphel)
            k=k+1
            uvarm_custom(noel,npt,k)=rho_tot_ph(iphel)
            k=k+1
            do kmo=1,nslmod(iphel)
               uvarm_custom(noel,npt,k)=rho_deb_ph_mod(kmo,iphel)
               k=k+1
            enddo
            uvarm_custom(noel,npt,k)=rho_recov
            k=k+1
c$$$         elseif (ihardlaw.eq.23) then
c$$$            write(*,*)'** uvarm_custom for ihardlaw of 23'
c$$$            write(*,*)' is not ready'
c$$$            call xit()
         endif


c$$$         k=1
c$$$         uvarm_custom(noel,npt,k)=alplpx
c$$$         k=k+1
c$$$         uvarm_custom(noel,npt,k)=rels
c$$$         k=k+1
c$$$         uvarm_custom(noel,npt,k)=reld
c$$$         k=k+1
c$$$         uvarm_custom(noel,npt,k)=vnorm(DBAR0,6)
c$$$
c$$$         dum1=0.
c$$$         do i=1,6
c$$$            do j=1,6
c$$$               dum1=dum1+mbar(i,j)**2
c$$$            enddo
c$$$         enddo
c$$$         dum1=sqrt(dum1)
c$$$
c$$$         k=k+1
c$$$         uvarm_custom(noel,npt,k)=dum1
c$$$
c$$$         dum1=0.
c$$$         do i=1,6
c$$$            do j=1,6
c$$$               dum1=dum1+LBAR(i,j)**2
c$$$            enddo
c$$$         enddo
c$$$         dum1=sqrt(dum1)
c$$$
c$$$         k=k+1
c$$$         uvarm_custom(noel,npt,k)=dum1
c$$$
c$$$



c$$$c         uvarm_custom(noel,npt,k)=rels !10
c$$$c         k=k+1
c$$$c         uvarm_custom(noel,npt,k)=reld !11
c$$$c         k=k+1
c$$$
c$$$         uvarm_custom(noel,npt,k)=dum1 !10 nstd(edot11)
c$$$         k=k+1
c$$$         uvarm_custom(noel,npt,k)=dum2 !11 nstd(sig11)
c$$$         k=k+1
c$$$
c$$$c         do i=1,2               !limits to only first two planes
c$$$c            uvarm_custom(noel,npt,k)=eps_w(i,1) !12,13,
c$$$c            k=k+1
c$$$c         enddo
c$$$         do i=1,2               !limits to only first two planes
c$$$            uvarm_custom(noel,npt,k)=wgtsetph(i,1) !14,15
c$$$            k=k+1
c$$$         enddo
c$$$
c$$$         do iphel=1,nphmx !
c$$$c     ** Accumulated relative activity
c$$$            do i=1,nmodmx       ! 16~20
c$$$               if (gavmod(i,iphel).lt.0) then
c$$$                  write(*,*)'negative gavmod(i,iphel) is found:',
c$$$     $                 gavmod(i,iphel)
c$$$                  write(*,*)'i:',i
c$$$                  write(*,*)'noel:',noel
c$$$                  write(*,*)'npt:',npt
c$$$c                  call xit()
c$$$               endif
c$$$
c$$$               if (gavmod(i,iphel).gt.1) then
c$$$                  write(*,*)'gavmod(i,iphel) greater than one is found',
c$$$     $                 gavmod(i,iphel)
c$$$                  write(*,*)'i:',i
c$$$                  write(*,*)'noel:',noel
c$$$                  write(*,*)'npt:',npt
c$$$c                  call xit()
c$$$               endif
c$$$
c$$$c               aa(i)=aa(i)+gavmod(i,iphel)*dtime
c$$$               aa(i)=aa(i)+contrib(2)*gavmod(i,iphel)*dtime
c$$$
c$$$               if (aa(i).lt.0) then
c$$$                  write(*,*)'negative aa(i) is found',aa(i)
c$$$                  write(*,*)'i:',i
c$$$                  write(*,*)'noel:',noel
c$$$                  write(*,*)'npt:',npt
c$$$                  ! call xit()
c$$$               endif
c$$$
c$$$               uvarm_custom(noel,npt,k)=aa(i)/(time(2)+dtime)
c$$$               if (uvarm_custom(noel,npt,k).gt.1.and..false.) then
c$$$                  write(*,*)'aa(i)/(t+dt) higher than 1 is found'
c$$$                  write(*,*)'val:',uvarm_custom(noel,npt,k)
c$$$                  write(*,*)'i:',i
c$$$                  write(*,*)'noel:',noel
c$$$                  write(*,*)'npt:',npt
c$$$                  ! call xit()
c$$$               endif
c$$$               k=k+1
c$$$            enddo
c$$$
c$$$c     ** Twin fraction
c$$$            do i=1,ntwmmx !21
c$$$               uvarm_custom(noel,npt,k)=twfrph(i,iphel)
c$$$               k=k+1
c$$$            enddo
c$$$c     ** sum of all accumulated relative activity for this increment
c$$$            uvarm_custom(noel,npt,k)=0d0
c$$$            i0 = k - nmodmx - ntwmmx
c$$$            do i=1,nmodmx !22
c$$$               uvarm_custom(noel,npt,k)=uvarm_custom(noel,npt,k)+
c$$$     $              uvarm_custom(noel,npt,i0+i)
c$$$            enddo
c$$$            k=k+1               !sum of accumulated activity !11
c$$$c     ** sum of all instant activity (would be 1, if fully plastic)
c$$$            uvarm_custom(noel,npt,k)=0d0
c$$$            do i=1,nmodmx !23
c$$$               uvarm_custom(noel,npt,k)=uvarm_custom(noel,npt,k)+
c$$$     $              gavmod(i,iphel)
c$$$            enddo
c$$$            k=k+1
c$$$c     **
c$$$            ! 24
c$$$            ac(iphel)=ac(iphel)+activ(iphel)*dtime
c$$$            uvarm_custom(noel,npt,k)=ac(iphel)/(time(2)+dtime)
c$$$            k=k+1               !sum of activity !12
c$$$c           ! 25
c$$$            uvarm_custom(noel,npt,k)=ac(iphel)/(time(2)+dtime)
c$$$            k=k+1
c$$$         enddo
      endif
