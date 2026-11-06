c     ----------------------------------------------------------------------------
c     subroutine calc_evptilde
c     ----------------------------------------------------------------------------
c     S: compliance
c     C: stiffness
c     evp: denote elasto-visco-plastic quantities
c     av:  denote avergage quantities
c     ----------------------------------------------------------------------------
c     The tilde (interaction) compliance is used in the elasto-viscoplastic
c     interaction equation; Here, <MAST> array denotes the interaction
c     compliance.
c
c     This subroutine calculates Eshelby tensor - both its symmetric and
c     skew-symmetric parts. From the symmetric part of Eshelby tensor
c     MAST tensor is estimated, which will be used for interaction
c     and localization equations.

c     Eshelby tensor is a function of
c         1) each grain shape, i.e., axisgr(0,:,kkk)
c         or 2) phase-specific average shape, i.e., axisph(0,:,iph)
c     Which of the above is used is determined by ISHAPE flag.
c     Therefore, depending on the ISHAPE flag, the elasto-visco-plastic
c     tilde (interaction) compliance could be either
c          MASTph (phase-specific) or MASTgr (grain-specific).

c     Likewise, the skew-symmetric part of Eshelby tensor could be a function
c     of either 1) each grain shape; or 2) phase-specific average shape.
c     Depending on ISHAPE flag, either 1) ASGR or 2) ASPH stores the
c     skew-symmetric part of the Eshelby tensor.
c     Since ISHAPE determines the type of array to be calculated, its four
c     separate cases are all specificed in below:
c     ----------------------------------------------------------------------------
c     ishape=0 : average grain shape is used for each grain
c     ishape=1 : average grain shape is used for each grain, but keeps track of
c                shape evolution (not used for calculation)
c     ishape=2 : use individial grain shape to calculate eshelby tensor
c     ishape=3 : same as ishape=2 but read initial grain ellipsoid orientation
c                from file
c     ishape=4 : same as ishape=2 but read initial grain ellipsoid orientation
c                and shape (principal axial lengths) from file.
c
c     Youngung Jeong@CWNU - yjeong@changwon.ac.kr
c     ----------------------------------------------------------------------------
      subroutine calc_evptilde(mdim,iregime,LBARTG,MBARTG,nphmx,nph,
     $     ngrmx,ngr,factor,axisph,axisgr,ishape,asph,asgr,MASTph,
     $     MASTgr,iesh_fail)
      implicit none
      integer,intent(in)::mdim,iregime,nphmx,nph,ngrmx,
     $     ishape(0:nphmx),ngr(0:nphmx)
      real*8, intent(in) ::LBARTG(6,6),MBARTG(6,6),factor(ngrmx),
     $     axisph(0:3,3,0:nphmx),axisgr(0:3,3,ngrmx)
      real*8, intent(out)::asph(3,3,3,3,nphmx),
     $     asgr(3,3,3,3,ngrmx),
     $     MASTph(mdim,mdim,nphmx),MASTgr(mdim,mdim,ngrmx)
      logical, intent(out):: iesh_fail
c     -------------------------------------------------------------
c     locals
      real*8 ESA(mdim,mdim)
      real*8 axb(3),eigb(3,3),
     $     eigbt(3,3),LBARc(3,3,3,3),ImSINV(mdim,mdim),e4sa(3,3,3,3),
     $     e4ga(3,3,3,3),r4ga(3,3,3,3),r4sa(3,3,3,3),ESA_inv(mdim,mdim),
     $     c4ga(3,3,3,3),ESAINVc(3,3,3,3),fs(mdim,mdim),mast(mdim,mdim),
     $     pdil
      real*8 auxm(mdim),auxmm(mdim,mdim),aux33(3,3),bux33(3,3),
     $     aux3333(3,3,3,3),bux3333(3,3,3,3)
      integer iph,i,j,k,l,k1,l1,kkk,xidmm(mdim,mdim),inverr,iopt_Esh
c     kronecker delta(should be defined only once)
      xidmm(:,:)=0
      do i=1,mdim
         xidmm(i,i)=1
      enddo

      if (iregime.eq.1)  then !VPSC
         iopt_Esh=2
      elseif (iregime.eq.2) THEN   !EVPSC
         iopt_Esh=1
      endif

      auxmm(:,:)=LBARTG(1:mdim,1:mdim)
      call chg_basis(auxm,aux33,auxmm,LBARc,3,mdim)        !LBARTG to Cartesian form

c     loop over phases
      do iph=1,nph
         if(ishape(iph).le.1) then
            do i=1,3
               axb(i)=axisph(0,i,iph) ! initial ellipsoid ratio
            do j=1,3
               eigb(i,j)=axisph(i,j,iph) ! transf. matrix
               eigbt(j,i)=axisph(i,j,iph)
            enddo
            enddo
         endif

c     -- loop over grains
      do kkk=ngr(iph-1)+1,ngr(iph)
         if(ishape(iph).ge.2) then
            do i=1,3
               axb(i)=axisgr(0,i,kkk) ! initial ellipsoid ratio
            do j=1,3
               eigb(i,j)=axisgr(i,j,kkk) ! transf. matrix
               eigbt(j,i)=axisgr(i,j,kkk)
            enddo
            enddo
         endif

c     -- Begin Eshelby calculation for each phase (or each grain)
         if (ishape(iph).ge.2.or.kkk.eq.ngr(iph-1)+1) then
c     ---------------------------------------
            call rot_4th_tensor(LBARc,eigbt,c4ga)
            pdil=0d0
            call eshelby(axb,c4ga,0d0,e4ga,r4ga,aux33,bux33,pdil,
     $           aux3333,bux3333,iopt_esh,iesh_fail)

            if (iesh_fail) then
               write(*,*)'iesh_fail occur in calc_evptilde.f'
               write(*,*)'axb:',axb
               write(*,*)'c4ga:',c4ga
               write(*,*)'lbarc:',lbarc
               return
            endif

            call rot_4th_tensor(e4ga,eigb,e4sa)
            call rot_4th_tensor(r4ga,eigb,r4sa)
c     ---------------------------------------
c           e4sa -> ESA
            call chg_basis(auxm,aux33,ESA,e4sa,4,mdim)

c           ESA; eshelby tensor in b-basis (6x6)

c     -- Now, calculates Sevp tilde
            ImSINV(:,:)=0d0
            do i=1,mdim
            do j=1,mdim
               ESA_inv(i,j)=ESA(i,j)
               ImSINV(i,j)=xidmm(i,j)-ESA(i,j) ! I-S
            enddo
            enddo

            call lu_inverse(inverr,ImSINV,mdim) !(I-S)^-1
            call lu_inverse(inverr,ESA_inv,mdim) ! S^-1

            call chg_basis(auxm,aux33,
     $           ESA_inv,ESAINVc,3,mdim)

c     -- Calculate (I-S)^(-1)*S
            DO I=1,MDIM
            DO J=1,MDIM
               FS(I,J)=0.
               DO K=1,MDIM
                  FS(I,J)=FS(I,J)+ImSINV(I,K)*ESA(K,J) ! (I-S)^(-1)*S
               ENDDO
            ENDDO
            ENDDO

c     -- Now, calculate sevp_tilde
c     -- Sevp~ = (I-S)^-1 : S : <Sevp> (this can be either grain- or phase-specific)
            mast(:,:)=0d0
            do i=1,mdim
            do j=1,mdim
c     --
            do k=1,mdim
               mast(i,j)=mast(i,j)+FS(i,K)*MBARTG(k,j)
            enddo
c            if(interaction.eq.4.AND.irsvar.eq.1) then
c               MAST(I,J)=jxrs*MAST(I,J) ! legacy option
c            else
               !(yj) M~=n M~
c               MAST(I,J)=NEFFGR(KGX)*MAST(I,J) ! Eq 5-32b
c            endif
            enddo
            enddo

c     --- delta-EVPSC and sigma-EVPSC --- This is likely to be replaced eventually by neffgr(ngx)?
           if (iregime.eq.2) then
c              mast(:,:)=mast(:,:)*factor(kkk)
c             Affecting only 5x5 matrix to allow transition from EL to VP regime.

c             Indices i and j here were deliberately set to 5.
              do i=1,5
              do j=1,5
                 mast(i,j)=mast(i,j)*factor(kkk)
              enddo
              enddo
           endif

c     -- calculates asph and asgr used in
C      'ASPH' OR 'ASGR' ARE USED IN SUBR. UPDATE_ORIENTATION TO CALCULATE
C      SPIN-RATE DEVIATIONS.


            if (ishape(iph).le.1) then
               MASTph(:,:,iph)=mast(:,:)
c              calculates only for average grain shapec
               asph(:,:,:,:,iph)=0d0
               do i=1,3
               do j=1,3
               do k=1,3
               do l=1,3
c     --
               do k1=1,3
               do l1=1,3
                  asph(i,j,k,l,iph)=asph(i,j,k,l,iph)+
     $                 r4sa(i,j,k1,l1)*ESAINVc(k1,l1,k,l)
               enddo
               enddo
               enddo
               enddo
               enddo
               enddo
            elseif (ishape(iph).ge.2) then
               MASTgr(:,:,kkk)=mast(:,:)
c              calculate grain-specific rotations
               asgr(:,:,:,:,kkk)=0d0
               do i=1,3
               do j=1,3
               do k=1,3
               do l=1,3
c     --
               do k1=1,3
               do l1=1,3
                  asgr(i,j,k,l,kkk)=asgr(i,j,k,l,kkk)
     $                 +r4sa(i,j,k1,l1)*ESAINVc(k1,l1,k,l)
               enddo
               enddo
               enddo
               enddo
               enddo
               enddo
            endif
         endif ! if (ishape(iph).ge.2.or.kkk.eq.ngr(iph-1)+1) then
c     -- End eshelby calculation for each phase (or each grain)
      enddo                     ! Loop over grains
      enddo                     ! Loop over phases

      return
      end subroutine calc_evptilde
