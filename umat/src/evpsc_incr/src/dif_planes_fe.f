c     this subroutine is used in combination with <diffraction>
      subroutine diff_stv_str(UR_SX,UR_DIF,filediff,ndifmx,
     $     ngr,sg,ag,wgt,mgel)
      implicit none
      integer, intent(in) :: UR_SX,UR_DIF,ndifmx,ngr
      real*8, intent(in) :: sg(6,ngr),ag(3,3,ngr),wgt(ngr),
     $     mgel(6,6,ngr)
      character(len=325), intent(in):: filediff

c     locals
      character*12 plane_indices(100)
      character(len=123) prosa
      integer idiffx,ndiff,nfamily(ndifmx),isn(4),isb(4),nd,i,j,nf,
     $     icrysym,nequiv,igrset(ndifmx,ngr),ngrsetph(ndifmx),
     $     ng,ipl,ngrx,igrsetx
      real*8 angdetector,sn(3),sb(3),sneq(3,24),eps,cospsvs,cosangdet,pi
     $     ,diff_vc(3,24,ndifmx),diff_vs(3,ndifmx),wgtsetph(ndifmx),
     $     ps(3),para_w(ndifmx),eps_w(ndifmx),eelgr(6,ngr),
     $     aux33(3,3),aux66(6,6),aux3333(3,3,3,3),sav(3,3),savg6(6),
     $     wgtall
c     average stress
      wgtall=0.
      do ng=1,ngr
         wgtall=wgtall+wgt(ng)
      enddo

      do i=1,6
         savg6(i)=0.
         do ng=1,ngr
            savg6(i)=savg6(i)+sg(i,ng)*wgt(ng)
         enddo
         savg6(i)=savg6(i)/wgtall
      enddo

      call chg_basis(savg6,sav,aux66,aux3333,1,6)

      pi=4.d0*datan(1.d0)

c     initialize crystal symmetry
      call crystal_symmetry(1,UR_SX,icrysym,isn,sn,sneq,isb,sb,nequiv)

c     read diffraction data
      read(UR_DIF,'(a)') prosa
      read(UR_DIF,'(a)') prosa
      read(UR_DIF,*) ndiff,angdetector
      if (ndiff.gt.ndifmx) then
         write(*,'('' number of diffraction planes is'',i4,/,
     $'' greater than code dimension'',i4)') ndiff,ndifmx
         write(*,*)
         write(*,'(1h ,''stop in routine *** dif_planes ***'')')
         call xit()
      endif

      read(UR_DIF,'(a)') prosa
      read(UR_DIF,'(a)') prosa
      do idiffx=1,ndiff
         read(UR_DIF,'(a12)') plane_indices(idiffx)
      enddo
      rewind(UR_DIF)
      do i=1,5
         read(UR_DIF,'(a)') prosa
         write(*,*)'prosa:',trim(prosa)
      enddo


c     read isn from UR_DIF and obtain <sn>, <nequiv>, and <sb
      do nd=1,ndiff
         call crystal_symmetry(3,UR_DIF,icrysym,isn,sn,
     $        sneq,isb,sb,nequiv)

         nfamily(nd)=nequiv
         do i=1,3
            diff_vs(i,nd)=sb(i) ! diff direction in sa; defined by chi, eta
            do nf=1,nfamily(nd)
               diff_vc(i,nf,nd)=sneq(i,nf) ! plane normal in ca, defined by (hkl)
            enddo
         enddo
      enddo


c     now that diffraction condition is defined (vs, vc),
c     began diffraction simulation and calculate lattice strains & weights

      cosangdet=abs(cos(angdetector*pi/180.0))

      do nd=1,ndiff
         ngrsetph(nd)=0
         wgtsetph(nd)=0.
         do ng=1,ngr
            do ipl=1,nfamily(nd)
               do i=1,3
                  ps(i)=0.0
                  do j=1,3
                     ps(i)=ps(i)+ag(i,j,ng)*diff_vc(j,ipl,nd)
                  enddo
               enddo

               cospsvs=0.
               do i=1,3
                  cospsvs=cospsvs+diff_vs(i,nd)*ps(i)
               enddo
               cospsvs=abs(cospsvs)

               if(cospsvs.ge.cosangdet) then
                  ngrsetph(nd)=ngrsetph(nd)+1
                  ngrx=ngrsetph(nd)
                  igrset(nd,ngrx)=ng !! the third axis of igrset extends from 1 to 4.
                  wgtsetph(nd)=wgtsetph(nd)+wgt(ng)
               endif

            enddo
         enddo
      enddo

c     now obtain elastic lattice strains & weights

      do nd=1,ndiff
         para_w(nd)=0.d0
         do ng=1,ngrsetph(nd)
            igrsetx=igrset(nd,ng)
c *** calculation of eelgr
            eelgr(:,igrsetx)=0.
            do i=1,6
               do j=1,6
                  eelgr(i,igrsetx)=eelgr(i,igrsetx)+mgel(i,j,igrsetx)*
     $                 sg(j,igrsetx)
               enddo
            enddo
            call chg_basis(eelgr(:,igrsetx),aux33,aux66,aux3333,1,6) ! to cartesian
            eps=0.d0
            do i=1,3
               do j=1,3
                  eps=eps+diff_vs(i,nd)*diff_vs(j,nd)*aux33(i,j)
               enddo
            enddo
            para_w(nd)=para_w(nd)+eps*wgt(igrsetx) ! weighted using grain fraction
         enddo                  ! end of 'do ng' over grains in subset
         if(wgtsetph(nd).ne.0.d0)
     $        eps_w(nd)=para_w(nd)/wgtsetph(nd) ! renormalize using set fraction
      enddo                     ! end of 'do nd' over diffracting planes

      write(*,*)'totngr:',ngr


      do i=1,3
         write(*,'(7x,a3,2i1)',advance='no')'sig',i,i
      enddo
      do nd=1,ndiff
         write(*,'(5x,a3,a1,a5,a1)',advance='no')
     $        'eps','(',trim(plane_indices(nd)),')'
      enddo
      do nd=1,ndiff
         write(*,'(2x,a6,a1,a5,a1)',advance='no')
     $        'ngrset','(',trim(plane_indices(nd)),')'
      enddo
      do nd=1,ndiff
         write(*,'(2x,a6,a1,a5,a1)',advance='no')
     $        'wgtset','(',trim(plane_indices(nd)),')'
      enddo
      do nd=1,ndiff
         write(*,'(1x,a7,a1,a5,a1)',advance='no')
     $        'ngr/tot','(',trim(plane_indices(nd)),')'
      enddo
      write(*,*)



      do i=1,3
         write(*,'(e12.4)',advance='no')sav(i,i)
      enddo


      do nd=1,ndiff
         write(*,'(3x,e12.4)',advance='no')eps_w(nd)
      enddo
      do nd=1,ndiff
         write(*,'(3x,i12)',advance='no')ngrsetph(nd)
      enddo
      do nd=1,ndiff
         write(*,'(3x,e12.4)',advance='no')wgtsetph(nd)
      enddo
      do nd=1,ndiff
         write(*,'(3x,e12.4)',advance='no')(ngrsetph(nd)*1.0)/(ngr*1.0)
      enddo
      write(*,*)

c$$$      do nd=1,ndiff
c$$$         write(*,'(e12.4,3x,i6,2f9.4)')
c$$$     $        eps_w(nd),ngrsetph(nd),wgtsetph(nd),
c$$$     $        (ngrsetph(nd)*1.0)/(ngr*1.0)
c$$$      enddo


      return
      end subroutine


c     -----------------------------------------------------------------




c     below is old one that is no longer used.
      subroutine diff_planes_fe(iph,filediff,
     $     diff_vs,diff_vc,ndiff,nfamily,
     $     icrysym,ioption)
      implicit none
      include "../../../../src/vpsc_src/incr_evpsc_fe.dim"
c      character*125 filediff(nphmx)

      real*8 ps(3),aux(3),aux33(3,3),aux66(6,6),aux3333(3,3,3,3)
      integer isn(4),isb(4),ng,nd,nf,ngrx,i,j,idiffx,iph,icrysym,
     $     nequiv,ioption,iunit,idum,ipl,icase,igrsetx
      real*8 sn(3),sneq(3,24),sb(4),eps,cospsvs,cosangdet,pi

      character*12 plane_indices(100)
      real*8 savc(3,3),para_w(ndifmx),eps_w(ndifmx,nphmx)
c     $     ,eelgr(6,ngrmx),wgtsetph(ndifmx,nphmx),
c     $     diff_vc(3,24,ndifmx,nphmx),diff_vs(3,ndifmx,nphmx),
c     $     angdetector(nphmx)
c      integer igrset(ndifmx,ngrmx,nphmx),ngrsetph(ndifmx,nphmx),
c     $     nfamily(ndifmx,nphmx),ndiff(nphmx)
c     iregime: 1('vpsc') or 2('delta-evpsc'), or 3 ('sigma-EVPSC')


      logical verbose
      parameter(verbose=.true.)
c      parameter(verbose=.false.)

      pi=4.d0*datan(1.d0)

      if (verbose) write(*,*)'elasticity'

c **********************************************************************
      if (ioption.eq.0) then
c **********************************************************************
c *** reads indices of diffracting plane and angles of diffraction direction.
c *** generates crystallographically equivalent plane normals 'sneq' in caxes.
c *** calculates diffraction direction 'sb' in saxes

        UR_DIF=7
c        open(unit=UR_DIF,file=filediff,status='old')
        write(*,*)'filediff(iph):',filediff
        open(unit=UR_DIF,file=filediff(iph),status='old') ! yj 2018.07.18
        iunit=19+iph
        open(iunit,file='dif_str_ph'//char(48+iph)//'.out',
     #                   status='unknown')
        iunit=22+iph
        open(iunit,file='dif_wgt_ph'//char(48+iph)//'.out',
     #                   status='unknown')
c -----------------------------------------------------------------
c     writes miller indices of diffracting planes in output files
        read(UR_DIF,'(a)') prosa
        read(UR_DIF,'(a)') prosa
        read(UR_DIF,*) ndiff(iph)
        read(UR_DIF,'(a)') prosa
        read(UR_DIF,'(a)') prosa
        do idiffx=1,ndiff(iph)
           read(UR_DIF,'(a12)') plane_indices(idiffx)
        enddo
        iunit=19+iph ! dif_str_ph
        write(iunit,'(''    epstot11    epstot22    epstot33'',
     $''       sav11       sav22       sav33'',
     $6x,100a12)') (plane_indices(i),i=1,ndiff(iph))
        iunit=22+iph ! dif_wgt_ph
        write(iunit,'(''    epstot11    epstot22    epstot33'',
     $''       sav11       sav22       sav33       wph'',
     $6x,100a12)') (plane_indices(i),i=1,ndiff(iph))

        rewind UR_DIF

c -----------------------------------------------------------------

        read(UR_DIF,'(a)') prosa
        read(UR_DIF,'(a)') prosa
        read(UR_DIF,*) ndiff(iph),angdetector(iph)
        if (ndiff(iph).gt.ndifmx) then
          write(*,'('' number of diffraction planes is'',i4,/,
     #    '' greater than code dimension'',i4)') ndiff(iph),ndifmx
          write(*,*)
          write(*,'(1h ,''stop in routine *** dif_planes ***'')')
          stop
        endif
        read(UR_DIF,'(a)') prosa
        read(UR_DIF,'(a)') prosa

        write(*,*)'ndiff(iph):',ndiff(iph)

        write(*,*)'nequiv:',nequiv

        do nd=1,ndiff(iph)
           write(*,*)'nd:',nd

           write(*,*)'nequiv (before):',nequiv

           call crystal_symmetry(3,UR_DIF,icrysym,isn,sn,
     $          sneq,isb,sb,nequiv)
           write(*,*)'nequiv (after) :',nequiv

           nfamily(nd,iph)=nequiv

        do i=1,3
           diff_vs(i,nd,iph)=sb(i) ! diff direction in sa
        do nf=1,nfamily(nd,iph)
           diff_vc(i,nf,nd,iph)=sneq(i,nf) ! plane normal in ca
        enddo
        enddo
        enddo

        close(unit=UR_DIF)

      endif   ! end of initialization run for ioption=0

c **********************************************************************
      if (ioption.eq.1) then
c **********************************************************************

         if (iregime.eq.1) then ! vpsc
            if(iDBARv(1)*iDBARv(2)*iDBARv(3).eq.1) then
               write(*,*)'---------- short notice ------------'
               write(*,*)'the feature of diff_planes to '
               write(*,*)'calculate internal strain eelgr'
               write(*,*)'is available for the case of '
               write(*,*)'iDBARv(1)*iDBARv(2)*iDBARv(3).ne.1'
               write(*,*)'with the given boundary condition'
               write(*,*)'vpsc cannot obtain eelgr since the'
               write(*,*)' macroscopic pressure is unknown.'
               write(*,*)'------------------------------------'
               call exit(-1)
            endif
         endif

c *** identifies and tags those grains that contribute to each diffraction
c     direction (nd) in the phase (iph) --> igrset(nd,ngrsetph,iph)

         cosangdet=abs(cos(angdetector(iph)*pi/180.0))

         do nd=1,ndiff(iph)
            ngrsetph(nd,iph)=0
            wgtsetph(nd,iph)=0.
         do ng=ngr(iph-1)+1,ngr(iph)
         do ipl=1,nfamily(nd,iph)
            do i=1,3
               ps(i)=0.0
            do j=1,3
               ps(i)=ps(i)+ag(i,j,ng)*diff_vc(j,ipl,nd,iph)
            enddo
            enddo

            cospsvs=0.
            do i=1,3
               cospsvs=cospsvs+diff_vs(i,nd,iph)*ps(i)
            enddo
            cospsvs=abs(cospsvs)

            if(cospsvs.ge.cosangdet) then
               ngrsetph(nd,iph)=ngrsetph(nd,iph)+1
               ngrx=ngrsetph(nd,iph)
               igrset(nd,ngrx,iph)=ng !! the third axis of igrset extends from 1 to 4.
               wgtsetph(nd,iph)=wgtsetph(nd,iph)+wgt(ng)
            endif

         enddo
         enddo
         enddo



c **********************************************************************
c *** average crystallographic strain for grains in each set
c *** what follows are different ways of calculating elastic strain in each
c     grain depending on the type of process being simulated
c
c     1- for an elastic proces the strain 'eelgr' has been calculated in the
c        subroutine 'elsc' and can be used without further calculation
c        (this is the case ivgvar.eq.-1 is the process, which is not available
c        in vpsc-fld-yld package of youngung)
c     2. except the first case, the elastic strain (eelgr) should be calculated.
c        prior to calculating the elastic strain, the stress state of each grain
c        can be categorized into three different scenarios.
c        2-1. the full stress tensor is available (that's with iregime=2 or 3),
c           no further treatment to obtain the full stress state of grains.
c        2-2. only tye deviatoric stress tensor is available without knowing the
c             hydrostatic stress. there are cases that the macroscopic hydrostatic
c             pressure is known, even for vpsc, if the boundary condition is imposed
c             such that the diagonal stress component is known. that that case,
c             the macroscopic stress <sav> is valid upto the 6th compoent.
c             under such circumstances, there are two possible options that allow
c             us to estimate the full grain stress.
c           2-2-1. calculate only the hydrostatic component from elasticity, sg(6).
c                  to do so, following is used:
c                  sg(6,ng) = bc_el(6, j, ng) * sav(j)
c           2-2-2. calculate the entire stress from elasticity. sg(1:6), that is:
c                  sg(i,ng) = bc_el(i, j, ng) * sav(j)
c     3. with sg obtained from (2-1, or 2-2-1, or 2-2-2), the internal strain can
c        be obtained from
c                  eelgr(i,igrsetx) = mgel(i,j,igrsetx) sg(j)

c                  icase=1       ! case 1
c                 icase=21 (i.e., iregime=1)
         icase=221              ! case 2-2-1
c         icase=222              ! case 2-2-2

         savc(:,:)=scauchy(:,:) ! already updated stress from subroutine <vpsc>
         if (iregime.eq.1) then ! vpsc
            call chg_basis(sav,savc, aux66,aux3333,2,6) !! now sav is in 6th dimension
         endif

         do nd=1,ndiff(iph)
            para_w(nd)=0.d0

            do ng=1,ngrsetph(nd,iph)
               igrsetx=igrset(nd,ng,iph)

c *** case 221 or 222 ***
               if (iregime.eq.1) then ! performed only when sg(6) is missing as in vpsc
                  if (icase.eq.221) then
                     sg(6,igrsetx)=0.d0
                     do j=1,6
                        sg(6,igrsetx)=sg(6,igrsetx)
     $                       +bc_el(6,j,igrsetx)*sav(j)
                     enddo
                  elseif (icase.eq.222) then
                     sg(:,igrsetx)=0.d0
                     do i=1,6
                     do j=1,6
                        sg(i,igrsetx)=sg(i,igrsetx)
     $                       +bc_el(i,j,igrsetx)*sav(j)
                     enddo
                     enddo
                  endif
               endif
c *** calculation of eelgr
               eelgr(:,igrsetx)=0.
               do i=1,6
               do j=1,6
                  eelgr(i,igrsetx)=eelgr(i,igrsetx)+mgel(i,j,igrsetx)*

c    $                 (sg(j,igrsetx)+sg_old(j,igrsetx)) ! sg now represents the 'delta' sigma in evp mod
     $                 sg(j,igrsetx)
               enddo


               enddo
               call chg_basis(eelgr(:,igrsetx),aux33,aux66,aux3333,1,6) ! to cartesian
               eps=0.d0
               do i=1,3
               do j=1,3
                  eps=eps+diff_vs(i,nd,iph)*diff_vs(j,nd,iph)*aux33(i,j)
               enddo
               enddo

               para_w(nd)=para_w(nd)+eps*wgt(igrsetx) ! weighted using grain fraction

            enddo               ! end of 'do ng' over grains in subset
            if(wgtsetph(nd,iph).ne.0.d0)
     $           eps_w(nd,iph)=para_w(nd)/wgtsetph(nd,iph) ! renormalize using set fraction
         enddo                  ! end of 'do nd' over diffracting planes

        iunit=19+iph
        write(iunit,'(7e12.4,6x,100e12.4)')
     $       epstot(1,1),epstot(2,2),epstot(3,3),
c    $    epstot(2,3),epstot(1,3),epstot(1,2),
     $       savc(1,1),savc(2,2),savc(3,3),wph(iph),
     $       (eps_w(nd,iph),nd=1,ndiff(iph))

        iunit=22+iph
        write(iunit,'(7e12.4,6x,100e12.4)')
     $       epstot(1,1),epstot(2,2),epstot(3,3),
c    $    epstot(2,3),epstot(1,3),epstot(1,2),
     $       savc(1,1),savc(2,2),savc(3,3),wph(iph),
     $       (wgtsetph(nd,iph),nd=1,ndiff(iph))

      endif            !end of ioption=1
c **********************************************************************

      return
      end subroutine diff_planes_fe
