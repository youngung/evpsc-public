c     ------------------------------------------------------------------
c     subroutine grain_stat calculates grain statitics such as the mean,
c     standard deviation, and normalized standard deviation

c     iopt=0 is valid only for stand-alone simulation, creates
c     'grain_stats.out' and 'grain_stats_nstd.out' files.
c     iopt=1 calculates the grain stats and save them to these created
c     files in case the run is using the stand-alone executable.
c
c     If parameter 'iall' is set to .true., additional file 'grain_stats_all.out'
c     is initilized (iopt=0) and the stress and strain states of individual
c     grains are all written down along with its weights.
c     This file can be used to obtain the 'histogram' - distribution of
c     files.
c     ------------------------------------------------------------------
      subroutine grain_stat(iopt,istep,ife,uw,uw_nstd,uw_alls,ngr,nphmx,
     $     nph,ngrmx,dtime,dgp,tg,dsg,sg_old,wgt,nstd_edot,nstd_sig,
     $     gtotgr)
      implicit none
c     ife==0 is reserved for stand-alone; 1 is for FE.
      integer iopt,nphmx,nph,ngr(0:nphmx),uw,uw_nstd,ife,
     $     ngrmx,istep,uw_alls
      real*8 dgp(5,ngrmx),tg(6,ngrmx),dsg(6,ngrmx),sg_old(6,ngrmx),
     $     wgt(ngrmx),dtime,nstd_edot(6),nstd_sig(6),
     $     gtotgr(ngrmx)
c     locals
      integer i,iv,j,k,kkk,iph,kind,icomp,ijv(6,2),n,m
      DATA ((IJV(N,M),M=1,2),N=1,6)/1,1,2,2,3,3,2,3,1,3,1,2/
      parameter(kind=5)
      logical iall
      parameter(iall=.false.)
c      parameter(iall=.true.)
      real*8 s(kind,6),u(kind,6),v(kind,6),aux33(3,3),bux33(3,3),aux6(6)
     $     ,bux6(6),aux66(6,6),p,aux5(5),aux3333(3,3,3,3),
     $     pressure_mean,pressure_v,pressure_std,
     $     devmag_mean,devmag_v,devmag_std,vnorm,
     $     gtot_mean,gtot_v,gtot_std,
     $     tg33(3,3),vpg33(3,3),dsg33(3,3),sg33(3,3),sgd33(3,3)


      if (iopt.eq.0.and.ife.eq.0) then
c     --------------------------------------------------------------
         open(uw,file='grain_stats.out',status='unknown')

         if (iall.and.istep.eq.0) then
            open(uw_alls,file='grain_stats_all.out',status='unknown')
         endif

         write(uw,'(a14)',advance='no') 'dt'
         do k=1,6
            i=ijv(k,1)
            j=ijv(k,2)
            write(uw,'(3x,a9,2i1.1)',advance='no')'tg_mean',i,j
            write(uw,'(3x,a9,2i1.1)',advance='no')'dsg_mean',i,j
            write(uw,'(3x,a9,2i1.1)',advance='no')'sg_mean',i,j
            write(uw,'(3x,a9,2i1.1)',advance='no')'dsgr_mean',i,j
         enddo
         do k=1,6
            i=ijv(k,1)
            j=ijv(k,2)
            write(uw,'(3x,a9,2i1.1)',advance='no')'tg_std',i,j
            write(uw,'(3x,a9,2i1.1)',advance='no')'dsg_std',i,j
            write(uw,'(3x,a9,2i1.1)',advance='no')'sg_std',i,j
            write(uw,'(3x,a9,2i1.1)',advance='no')'dsgr_std',i,j
         enddo
         do k=1,6
            i=ijv(k,1)
            j=ijv(k,2)
            write(uw,'(3x,a9,2i1.1)',advance='no')'tg_nstd',i,j
            write(uw,'(3x,a9,2i1.1)',advance='no')'dsg_nstd',i,j
            write(uw,'(3x,a9,2i1.1)',advance='no')'sg_nstd',i,j
            write(uw,'(3x,a9,2i1.1)',advance='no')'dsgr_nstd',i,j
         enddo
         write(uw,'(3x,a11)',advance='no')'press_mean'
         write(uw,'(3x,a11)',advance='no')'press_std'
         write(uw,'(3x,a11)',advance='no')'press_nstd'
         write(uw,'(3x,a11)',advance='no')'devmag_mean'
         write(uw,'(3x,a11)',advance='no')'devmag_std'
         write(uw,'(3x,a11)',advance='no')'devmag_nstd'
         write(uw,*)


         open(uw_nstd,file='grain_stats_nstd.out',status='unknown')
         do k=1,3
            i=ijv(k,1)
            j=ijv(k,2)
            write(uw_nstd,'(a12,2i1.1)',advance='no')'edot',i,j
         enddo
         do k=1,3
            i=ijv(k,1)
            j=ijv(k,2)
            write(uw_nstd,'(a12,2i1.1)',advance='no')'sig',i,j
         enddo
         do k=1,3
            i=ijv(k,1)
            j=ijv(k,2)
            write(uw_nstd,'(a12,2i1.1)',advance='no')'edotvp',i,j
         enddo
         write(uw_nstd,'(2a14)')'gtotgr_mean','gtotgr_nstd'

      elseif (iopt.eq.1) then
c     --------------------------------------------------------------
c     mean, standard deviation, normalized standard deviation on b-basis components
         u(:,:)=0.
         v(:,:)=0.
         s(:,:)=0.

c        weighted average values using cartesian components (not b-basis components)
         do iph=1,nph
         do kkk=ngr(iph-1)+1,ngr(iph)
            call chg_basis(tg(:,kkk),tg33,aux66,aux3333,1,6)
            call chg_basis(dgp(:,kkk),vpg33,aux66,aux3333,1,5)
            call chg_basis(dsg(:,kkk),dsg33,aux66,aux3333,1,6)
            call chg_basis(sg_old(:,kkk),sg33,aux66,aux3333,1,6)
            sg33(:,:)=sg33(:,:)+dsg33(:,:)
            sgd33(:,:)=dsg33(:,:)/dtime

            do iv=1,6
               i=ijv(iv,1)
               j=ijv(iv,1)
c           tg
               u(1,iv)=u(1,iv)+tg33(i,j)*wgt(kkk)
c           dsg
               u(2,iv)=u(2,iv)+dsg33(i,j)*wgt(kkk)
c           dsg+sg_old
               u(3,iv)=u(3,iv)+sg33(i,j)*wgt(kkk)
c           dsg/dtime
               u(4,iv)=u(4,iv)+sgd33(i,j)*wgt(kkk)
c           vp strain rate
               u(5,iv)=u(5,iv)+vpg33(i,j)*wgt(kkk)
            enddo
         enddo
         enddo


c        standard deviations
         do iph=1,nph
         do kkk=ngr(iph-1)+1,ngr(iph)
            call chg_basis(tg(:,kkk),tg33,aux66,aux3333,1,6)
            call chg_basis(dsg(:,kkk),dsg33,aux66,aux3333,1,6)
            call chg_basis(sg_old(:,kkk),sg33,aux66,aux3333,1,6)
            sg33(:,:)=sg33(:,:)+dsg33(:,:)
            sgd33(:,:)=dsg33(:,:)/dtime
            do iv=1,6
               i=ijv(iv,1)
               j=ijv(iv,1)
c               do kkk=1,ngr
c           tg
               v(1,iv)=v(1,iv)+(tg33(i,j)-u(1,iv))**2*wgt(kkk)
c           dsg
               v(2,iv)=v(2,iv)+(dsg33(i,j)-u(2,iv))**2*wgt(kkk)
c           dsg+sg_old
               v(3,iv)=v(3,iv)+
     $              (sg33(i,j)-u(3,iv))**2*wgt(kkk)
c           dsg/dtime
               v(4,iv)=v(4,iv)+(sgd33(i,j)-u(4,iv))**2*wgt(kkk)
c           vp strain rate
               v(5,iv)=v(5,iv)+(vpg33(i,j)-u(5,iv))**2*wgt(kkk)
            enddo
         enddo
         enddo



c         normalized standard deviations
         do i=1,6
            do j=1,kind
               v(j,i)=sqrt(v(j,i))
            enddo
            do j=1,kind
               if (abs(u(j,i)).gt.1e-10) then
                  s(j,i)=v(j,i)/u(j,i)
               else
                  s(j,i)=0.0
               endif
            enddo
         enddo ! end of loop over 6 components

c         write(*,*)'avg of u:',u(1,:)
c         write(*,*)'nstd of tg:',s(1,:)


c     --------------------------------------------------------------
ccc   pressure mean, std, normalized std
c     pressure mean  (mean)
         pressure_mean=0.
         gtot_mean=0.
         do iph=1,nph
         do kkk=ngr(iph-1)+1,ngr(iph)
c         do kkk=1,ngr
            aux6(:)=dsg(:,kkk)+sg_old(:,kkk)
            call chg_basis(aux6,aux33,aux66,aux3333,1,6)
            p=(aux33(1,1)+aux33(2,2)+aux33(3,3))/3.
            pressure_mean=pressure_mean+p*wgt(kkk)
            gtot_mean=gtot_mean+gtotgr(kkk)*wgt(kkk)
         enddo
         enddo
c     pressure std   (v)
         pressure_v=0.
         gtot_v=0.
         do iph=1,nph
         do kkk=ngr(iph-1)+1,ngr(iph)
            aux6(:)=dsg(:,kkk)+sg_old(:,kkk)
            call chg_basis(aux6,aux33,aux66,aux3333,1,6)
            p=(aux33(1,1)+aux33(2,2)+aux33(3,3))/3.
            pressure_v=pressure_v+(p-pressure_mean)**2*wgt(kkk)
            gtot_v=gtot_v+(gtotgr(kkk)-gtot_mean)**2*wgt(kkk)
         enddo
         enddo
         pressure_v=sqrt(pressure_v)
         gtot_v=sqrt(gtot_v)
c        pressure normalized std (std)
         pressure_std=pressure_v/pressure_mean
         if (abs(gtot_mean).gt.1e-10) then
            gtot_std=gtot_v/gtot_mean
         else
            gtot_std=0.
         endif
c     --------------------------------------------------------------
ccc   magnitude of deviatoric stress (mean, std, and normalized std)
c     magnitude of deviatoric stress (mean; mean value)
         devmag_mean=0.
         do iph=1,nph
         do kkk=ngr(iph-1)+1,ngr(iph)
            aux5(:)=dsg(:5,kkk)+sg_old(:5,kkk)
            devmag_mean=devmag_mean+vnorm(aux5,5)*wgt(kkk)
         enddo
         enddo
c     magnitude of deviatoric stress (v; standard deviation)
         devmag_v=0.
         do iph=1,nph
         do kkk=ngr(iph-1)+1,ngr(iph)
            aux5(:)=dsg(:5,kkk)+sg_old(:5,kkk)
            devmag_v=devmag_v+(vnorm(aux5,5)-devmag_mean)**2*wgt(kkk)
         enddo
         enddo
         devmag_v=sqrt(devmag_v)
c     magnitude of deviatoric stress (std; normalized standard deviation)
         devmag_std=devmag_v/devmag_mean

c     --------------------------------------------------------------
         nstd_edot(:)=s(1,:)    !normalized std of edot11
         nstd_sig(:)=s(3,:)     !normalized std of sig11

         if (ife.eq.0) then
c           grain_stats.out
            write(uw,'(130e14.6)')dtime,((u(i,icomp),i=1,kind),
     $           icomp=1,6),((s(i,icomp),i=1,kind),icomp=1,6),
     $           ((v(i,icomp),i=1,kind),icomp=1,6),
     $           pressure_mean,pressure_std,pressure_v,
     $           devmag_mean,devmag_std,devmag_v
c           grain_stats_nstd.out
            write(uw_nstd,'(11e14.6)')(s(1,i),i=1,3),(s(3,i),i=1,3),
     $           (s(5,i),i=1,3),gtot_mean,gtot_std

            if (iall) then
               if (istep.eq.1) then
                  write(uw_alls,'(a7,a4,a6,8a12)')'istep','iph','kkk',
     $                 'tg11','tg22','tg33','sg11','sg22','sg33',
     $                 'evmd',
     $                 'wgt'
               endif
               do iph=1,nph
               do kkk=ngr(iph-1)+1,ngr(iph)
                  aux6(:)= tg(:6,kkk)
                  bux6(:)=dsg(:6,kkk)+sg_old(:6,kkk)

                  call chg_basis(aux6,aux33,aux66,aux3333,1,6)
                  call chg_basis(bux6,bux33,aux66,aux3333,1,6)
                  write(uw_alls,'(i7,i4,i6)',advance='no')
     $                 istep,iph,kkk
                  do k=1,3
                     i=ijv(k,1)
                     j=ijv(k,2)
                     write(uw_alls,'(e12.4)',advance='no')aux33(i,j)
                  enddo
                  do k=1,3
                     i=ijv(k,1)
                     j=ijv(k,2)
                     write(uw_alls,'(e12.4)',advance='no')bux33(i,j)
                  enddo
                  write(uw_alls,'(2e12.4)')vnorm(bux6,6),wgt(kkk)
               enddo
               enddo
            endif
         endif

      endif

      return
      end subroutine grain_stat
