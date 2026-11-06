         write(*,*)'=============================='
         write(*,*)' Static Recrystallization'
c        External arguments specified:
c          - iph
c          - time_srx
c          - dtime
c          - tempdot

         nsteps=int(time_srx/dtime)
         if (iregime.eq.1) then
c           in case of VPSC, there's no need to solve the equation?
         elseif (iregime.eq.2) then
c           in case of EVPSC, the stress free condition may lead to
c           thermal expansion, thus is necessary to solve the EVP
c           constitutive equation
         endif

         igb=ngr(iph-1)+1
         igt=ngr(iph)
         iopt=1                 ! Static recrystallization
         write(*,'(a6,2a11)')'istep','time','temp'
         do istep=1,nsteps
            write(*,'(i6,2f11.2)')istep,time,temp
            temp=temp+tempdot*dtime
            time=time+dtime
            include "call_rexdd.f"
            CALL WRITE_TEXTURE_fe(UW_SRX_TEX-1+iph,
     $           iph,ngrmx,nphmx,igb,igt,
     $           EPSVM,ag,wgt,axisph,eulerph)
         enddo
         write(*,*)'=============================='
