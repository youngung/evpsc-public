      alplpx_diff=alplpx-alplpx_old ! to be removed, along with alplpx_old?

c     -----
c     checking out quantities associated with numerical stability.
      call mutexlock(5)
      if (alplpx.lt.alplstats(1,1))  alplstats(1,1)=alplpx
      if (alplpx.gt.alplstats(1,2))  alplstats(1,2)=alplpx
      if (RELS.lt.alplstats(2,1))  alplstats(2,1)=RELS
      if (RELS.gt.alplstats(2,2))  alplstats(2,2)=RELS
      if (RELD.lt.alplstats(3,1))  alplstats(3,1)=RELD
      if (RELD.gt.alplstats(3,2))  alplstats(3,2)=RELD
      if (it2.lt.alplstats(4,1)) alplstats(4,1)=it2*1.
      if (it2.gt.alplstats(4,2)) alplstats(4,2)=it2*1.
      if (epsvm.lt.alplstats(5,1)) alplstats(5,1)=epsvm
      if (epsvm.gt.alplstats(5,2)) alplstats(5,2)=epsvm
      if (svm.lt.alplstats(6,1)) alplstats(6,1)=svm
      if (svm.gt.alplstats(6,2)) alplstats(6,2)=svm
      if (vnorm(dstran,ntens).lt.alplstats(7,1))
     $     alplstats(7,1)=vnorm(dstran,ntens)
      if (vnorm(dstran,ntens).gt.alplstats(7,2))
     $     alplstats(7,2)=vnorm(dstran,ntens)
      call mutexunlock(5)
c     -----

      if (idiag_nr.ne.0) then
         call mutexlock(5)
         write(*,*)'---------------------------------------------------'
         write(*,*)'** NR failed to obtain converged GRAIN STRESS state'
         write(*,*)'---------------------------------------------------'
         call mutexunlock(5)
      endif
      if (idiag_sc.eq.1) then
         call mutexlock(5)
         write(*,*)'---------------------------------------------------'
         write(*,*)'** Exceeded max external iteration'
         write(*,*)'---------------------------------------------------'
         call mutexunlock(5)
      endif
      if (idiag_sc.eq.2)then
         call mutexlock(5)
         write(*,*)'---------------------------------------------------'
         write(*,*)'** inverse mat. in update_mac..'
         write(*,*)'---------------------------------------------------'
         call mutexunlock(5)
      endif
      if (idiag_sc.eq.3) then
         call mutexlock(5)
         write(*,*)'---------------------------------------------------'
         write(*,*)'** Failed eshelby'
         write(*,*)'---------------------------------------------------'
         call mutexunlock(5)
      endif


      if (idiag_nr.ne.0.or.idiag_sc.gt.1) then
         call mutexlock(5)
         write(*,*)
         do i=1,128
            write(*,'(a1)',advance='no')'-'
         enddo
         write(*,*)
         write(*,'(3a5,2a3)',advance='no')'st','inc','el','np','sp'
         write(*,'(4a9)',advance='no')'de','edot','evm','svm'
         write(*,'(4a9)',advance='no')'temp','dt','time1','time2'
         write(*,'(2a4)',advance='no')'irg','ifs'
         write(*,'(3a9)',advance='no')'alplpx','appx_old','gamd0s'
         write(*,*)
         write(*,'(3i5,2i3)',advance='no')jstep(1),kinc,noel,npt,kspt
         write(*,'(4e9.2)',advance='no')vnorm(dstran,ntens),
     $        vnorm(dstran,ntens)/dtime,epsvm,svm
         write(*,'(4e9.2)',advance='no')temp,dtime,time(1),time(2)
         write(*,'(2i4)',advance='no')iregime,ioptfs
         write(*,'(3e9.2)',advance='no')alplpx,alplpx_old,gamd0s(1,1)
         write(*,*)
         do i=1,128
            write(*,'(a1)',advance='no')'-'
         enddo
         write(*,*)
         call mutexunlock(5)
         pnewdt=0.8
         return
      endif
