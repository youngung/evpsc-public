C ***************************************************************************
C     SUBROUTINE GRAIN_INFO ---> VERSION 01/SEP/2023 adapted to dEVPSC code
C      ** This subroutine is not actively used in EVPSC but may be useful later.
C
C     FOR EACH GRAIN CALCULATES INCR & ACCUM VON MISES STRAIN,
C     ACCUMULATED PLASTIC WORK, AND 'TAYLOR' FACTOR DEFINED AS:
C          Mg=(Sg(i).DBAR(i))/TAU(1,g)/|DBAR|
C     CALCULATES PLASTIC WORK RATE USING MACROSCOPIC STRESS & STRAIN RATE
C     (wrate1) AND AVERAGE OF GRAIN STRESS & GRAIN RATE TENSORS (wrate2)
C ***************************************************************************

c     ndim: the actual dimension of stress (or strain) tensor components.
c           This can be either 5 (VPSC) or 6 (dEVPSC), that is given as
c           the argument to the subroutine.
c     nuse: The number of components used for the actual calculation
c           The number of arguments used in the actual calculation.
c           Since plasticity involves only deviatoric components,
c           when b-basis used, the first 5 components are important.

      subroutine grain_info(
c     intent(in)
     $     ndim,nphmx,ngrmx,nsysmx,ngr,time_incr,dav,dbar,sav,sg_old,dg,
     $     crss,wgt,iphbot,iphtop
c     intent(inout)
     $     ,workgr
c     intent(out)
     $     ,sgvm,egvm,egvmstep,taylgr
     $     )
      implicit none
      integer, intent(in)::ndim,ngrmx,nsysmx,nphmx,ngr(0:nphmx),iphbot,
     $     iphtop
      real*8, intent(in):: time_incr,dav(ndim),dbar(ndim),sav(ndim),
     $     sg_old(ndim,ngrmx),dg(ndim,ngrmx),wgt(ngrmx)
      real*8, intent(inout)::workgr(ngrmx),crss(nsysmx,ngrmx)
      real*8, intent(out):: sgvm(ngrmx),egvm(ngrmx),egvmstep(ngrmx),
     $     taylgr(ngrmx)
c     locals (yet)
      integer nuse
      parameter(nuse=5)
      real*8 vnorm,dbarnorm,dgvmx,sgvmx,taylav,taylx,workav,
     $     wrate,wrate1,wrate2
      integer i,iph,kkk,kgx

      if (nuse.gt.ndim) then
         write(*,*)'----------------------------------'
         write(*,*)'** Error: nuse is higher than ndim'
         write(*,*)'----------------------------------'
         call xit()
      endif

      dbarnorm=vnorm(dbar(:nuse),nuse)    ! only deviatoric components
      taylav=0.
      workav=0.

      wrate1=0.
      do i=1,nuse                  ! only deviatoric components
        wrate1=wrate1+sav(i)*dav(i)
      enddo

      wrate2=0.
      kgx=1
      do iph=iphbot,iphtop
      do kkk=ngr(iph-1)+1,ngr(iph)
        sgvmx=0.
        wrate=0.
        taylx=0.
        do i=1,nuse                ! only deviatoric components
          sgvmx=sgvmx+sg_old(i,kkk)**2
          wrate=wrate+sg_old(i,kkk)*dg(i,kgx)
          taylx=taylx+sg_old(i,kkk)*dbar(i)
        enddo
c     deviatoric von mises stress
        sgvm(kgx)=sqrt(sgvmx*3./2.)
c     defined as the conjugate vm
        if(sgvmx.ne.0.) dgvmx=wrate/sgvm(kgx)
        if(sgvmx.eq.0.) dgvmx=0.

c     deviatoric von mises
        egvmstep(kkk) =dgvmx*time_incr

        if (kkk.le.3.and..false.) then
           write(*,'(a17,x,e11.4)')'** egvmstep(kkk):',egvmstep(kkk)
           write(*,'(a17,x,e11.4)')'** dgvmx        :',dgvmx
        endif

c     deviatoric von mises
        egvm(kkk)  =egvm(kkk)   +egvmstep(kkk)
c     deviatoric non von mises
        workgr(kkk)=workgr(kkk) +wrate*time_incr
c     deviatoric non von mises
        taylgr(kkk)=taylx/crss(1,kkk)/dbarnorm

c *** the following are all non vm. get written in subr write_stress_strain
        taylav=taylav+taylgr(kkk)*wgt(kkk)
        workav=workav+workgr(kkk)*wgt(kkk)
        wrate2=wrate2+wrate      *wgt(kkk)

        kgx=kgx+1
      enddo
      enddo

      return
      end subroutine grain_info
