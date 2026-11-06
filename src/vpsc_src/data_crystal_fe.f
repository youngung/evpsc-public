c     ------------------------------------------------------------------
c     ** Table of isx_ver **
c     ------------------------------------------------------------------
c      isx_ver   CTE/EMOD file  CTE/EMOD file exists   non-schmid factor
c     ------------------------------------------------------------------
c         0           x                   -                   x
c         1           o                   o                   x
c         2           o                   x                   x
c         3           x                   -                   o
c         4           o                   o                   o
c         5           o                   x                   o
c     ------------------------------------------------------------------
c
c     subroutine data crystal reads crystal information including
c     crystal structure, symmetry, thermal expansion coefficient,
c     elastic constants, and constitutive (hardening) behavior on the
c     granular scale.
      SUBROUTINE DATA_CRYSTAL_fe(
c          in
     $     UR_SX,UR_ECTE,isx_ver,mxndeg,IPH,ntwmmx,
     $     nmodmx,nsysmx,ntwsmx,nphmx,
c          out
     $     celcc,athcc,selcc,nmodes,nslmod,ntwmod,nsyst,NSLSYS,ntwsys,
     $     twsh,nsm,dnca,dbca,schca,nschca,c6,isense,itwtype,xmu_mode,
     $     icrysym,fstar,ndegs_fit,cijv_fit,cte_fit)
c     when called in FE, iph means iphel
      implicit none
      integer, intent(in):: UR_SX,UR_ECTE,
     $     isx_ver,mxndeg,iph,ntwmmx,nmodmx,
     $     nsysmx,ntwsmx,nphmx
      integer, intent(out):: nmodes(nphmx),nslmod(nphmx),ntwmod(nphmx),
     $     nsyst(nphmx),NSLSYS(nphmx),ntwsys(nphmx),icrysym(nphmx),
     $     nsm(nmodmx,nphmx),isense(nsysmx,nphmx),ITWTYPE(NSYSMX,nphmx),
     $     ndegs_fit(nphmx,2)
      real*8,intent(out):: celcc(6,6,nphmx),athcc(6,nphmx),
     $     selcc(6,6,nphmx),twsh(ntwmmx,nphmx),dnca(3,nsysmx,nphmx),
     $     dbca(3,nsysmx,nphmx),schca(5,nsysmx,nphmx),
     $     nschca(6,nsysmx,nphmx),C6,xmu_mode(nmodmx,nphmx),fstar,
     $     cijv_fit(nphmx,0:mxndeg,6,6),cte_fit(nphmx,0:mxndeg,6)
c     locals
      integer MODE(50)
      DATA MODE/50*0/
      integer ISN(4),ISB(4),inverr,mcount,mloop,modex,i,j,ns,icomp,k,l,
     $     is,nind,nsmx,isensex,itwtypex,nmodesx,npoles,iclosepos,
     $     icloseneg,isys,isbot,istop,ics
      real*8 SN(3),SNEQ(3,24),SB(3),celccv(6,6),athccv(6),selccv(6,6),
     $     ST(3),temp,AUX5(5),AUX33(3,3),AUX3333(3,3,3,3),aux6(6),
     $     aux66(6,6),bux3333(3,3,3,3),PS(3,3),PNS(3,3),BUX6(6),
     $     BUX33(3,3),gam,bulk,bet,c44iso,c11iso,c12iso,xnu,prod,twshx,
     $     residual,alf,snbt(nsysmx,3,3),dum
      logical verbose
      parameter(verbose=.false.)
c      parameter(verbose=.true.)
      character*120 prosa
c     ---- non-Schmid factors c1,..,c6.
      integer inonschlaw
      parameter(inonschlaw=1)
c     inonschlaw 0: Using Yalcinkaya et al.'s formulation
c     inonschlaw 1: Using Groger's formulation
c     C5 is void for Groger's formulation
c     Commonly for both Yalcinkaya and Groger,
c     C6 is for the influence of hydrostatic pressure.
      real*8 C0,C1,C2,C3,C4,C5!,C6

c     used for Groger's non-Schmid formulation
      real*8 pi,nprime(3),tprime(3),wij(3,3),rij(3,3),c60,s60
      parameter(C0=0.5)

      pi=4d0*atan(1d0)

C -----------------------------------------------------------------
C *** READS CRYSTAL SYMMETRY & UNIT CELL PARAMETERS. CALCULATES CELL VECTORS.
C *** GENERATES ALL SYMMETRY OPERATIONS ASSOCIATED WITH 'CRYSYM'.

      if (nphmx.eq.0 .or.iph.eq.0)then
         write(*,*)'Error, neither nphmx nor iph should be zero'
         call xit()
      endif

      CALL CRYSTAL_SYMMETRY(1,UR_SX,icrysym(iph),
     $     ISN,SN,SNEQ,ISB,SB,NPOLES)
      NIND=3
      IF(ICRYSYM(iph).EQ.2 .OR. ICRYSYM(iph).EQ.3) NIND=4
      celccv(:,:)=0d0
      if (isx_ver.eq.2) then
         write(*,*)'** ignoring thermal_crystal_property.dat file'
         write(*,*)'   because that file could not be found'
         read(UR_SX,'(a)') prosa
         read(UR_SX,'(a)') prosa
      endif

      if (isx_ver.eq.0.or.isx_ver.eq.2
     $     .or.isx_ver.eq.3.or.isx_ver.eq.5) then
C *** READS SINGLE CRYSTAL ELASTIC STIFFNESS
         READ(UR_SX,'(A)') PROSA
c         ** C: stiffness, S: compliance
         READ(UR_SX,*) ((celccv(I,J),J=1,6),i=1,6)
C *** READS SINGLE CRYSTAL THERMAL EXPANSION COEFFICIENTS
         READ(UR_SX,'(A)') PROSA
         READ(UR_SX,*)     (athccv(I),I=1,6)
C *** Calculates & Stores SX MODULI (b-basis & Voigt notation)
      elseif(isx_ver.eq.1.or.isx_ver.eq.4) then
c     Here, although b-basis-notated elastic modulus and CTEs are
c     calculated from subroutine thermo_emod_cte,
c     voigt-notated elastic modulus and CTE are calculated from b-basis
c     quantities. That way, the legacy lines based on the original
c     single crystal file format are reused.
         do i=1,11
            read(UR_SX,'(a)')prosa
         enddo
         temp=298.
         call thermo_emod_cte(
c     intent(in)
     $        0,UR_ECTE,iph,nsysmx,nmodmx,nphmx,nmodes,mxndeg,nsm,
     $        temp,schca,
c     intent(inout)
     $        ndegs_fit,cijv_fit,cte_fit,
c         celcc and celcc are in b-basis
     $        celcc,selcc,xmu_mode,athcc
     $        )
         aux66(:,:)=celcc(:,:,iph)
         call chg_basis(aux6,aux33,aux66,aux3333,3,6)
         call voigt(aux6,aux33,celccv,aux3333,4)
         aux66(:,:)=selcc(:,:,iph)
         call chg_basis(aux6,aux33,aux66,aux3333,3,6)
         call voigt(aux6,aux33,selccv,aux3333,4)
      else
         write(*,*)'**Error: Unexpected isx_ver <data_crystal_fe.f>'
         call xit()
      endif

      celcc(:,:,iph)=0d0
      bux3333(:,:,:,:)=0d0
      aux3333(:,:,:,:)=0d0
c     convert voigt notation to cartesian
      call voigt(aux6,aux33,celccv,bux3333,3)
c     convert cartesian to in b basis
      call chg_basis(aux6,aux33,aux66,bux3333,4,6)
      celcc(:,:,iph)=aux66(:,:)
      inverr=0
      call lu_inverse(inverr,aux66,6)
      if (inverr.ne.0) then
         write(*,*)'** error in data_crystal_fe (1)'
         write(*,*)'celccv'
         write(*,'(6e12.4)') ((celccv(i,j),j=1,6),i=1,6)
         write(*,*)'c_1111, c_1122, c_1212 after first voigt'
         write(*,'(3e12.4)') bux3333(1,1,1,1),bux3333(1,1,2,2),
     $        bux3333(1,2,1,2)
         write(*,*)'aux66 (celcc(:,:,iph))'
         do i=1,6
            write(*,'(6e12.4)') (aux66(i,j),j=1,6)
         enddo
         write(*,*)'celcc(:,:,iph)'
         do i=1,6
            write(*,'(6e12.4)') (celcc(i,j,iph),j=1,6)
         enddo
         write(*,*)'UR_SX:',UR_SX
         call xit()
      endif
      selcc(:,:,iph)=aux66(:,:)
      call chg_basis(aux6,aux33,aux66,aux3333,3,6)
      call voigt(aux6,aux33,selccv,aux3333,4) ! compliance
      call voigt(athccv,aux33,aux66,aux3333,1)
      call chg_basis(aux6,aux33,aux66,aux3333,2,6) ! thermal
      ATHCC(:,iph)=aux6(:)

C --> VERIFICATION OF CORRECT INVERSION
      RESIDUAL=0d0
      DO I=1,6
      DO J=1,6
         AUX66(I,J)=0d0
      DO K=1,6
         AUX66(I,J)=AUX66(I,J)+celcc(I,K,iph)*SELCC(K,J,iph)
      ENDDO
      RESIDUAL=RESIDUAL+ABS(AUX66(I,J)-(I/J)*(J/I))
      ENDDO
      ENDDO

      if(verbose) WRITE(*,'('' CHECKING THAT CELCC*SELCC-ID6=0 '',
     $E15.7)')  RESIDUAL

C *** UPPER BOUND ISOTROPIC ELASTIC CONSTANTS FOR A RANDOM PX
C *** UPPER BOUND BULK MODULUS (Delta_P/Delta_V/V) FOR A
C *** TEXTURED POLYCRYSTAL
      ALF=(Celccv(1,1)+Celccv(2,2)+Celccv(3,3))/3.
      BET=(Celccv(1,2)+Celccv(1,3)+Celccv(2,3))/3.
      GAM=(Celccv(4,4)+Celccv(5,5)+Celccv(6,6))/3.
      BULK=(3.*ALF+6.*BET)/9.
      C11ISO=(3.*ALF+2.*BET+4.*GAM)/5.
      C12ISO=(   ALF+4.*BET-2.*GAM)/5.
      C44ISO=(C11ISO-C12ISO)/2.
      XNU=(C11ISO-2.*C44ISO)/(2.*C11ISO-2.*C44ISO)

C *** ESTIMATES ELASTIC BULK AND COMPRESSIBILITY MODULI (PROVISORY)
C     ALF=(CELCCV(1,1)+CELCCV(2,2)+CELCCV(3,3))/3.
C     BET=(CELCCV(1,2)+CELCCV(1,3)+CELCCV(2,3))/3.
C     BULK=(3.*ALF+6.*BET)/9.


c           Convention following Yalcinkaya et al.
            ! C0=0.5
c          -----------------------------------------------------------
c                                 |  Yalcinkaya (b->m,n->n,z=mxn)
c                                 ------------------------------------
c                                 |  Symbol              value for FE
c          -----------------------------------------------------------
c$$$            C1=0 ! b x b          |  m x m               0.0544
c$$$            C2=0 ! n x n          |  n x n              -0.0293
c$$$            C3=0 ! t x t          |  z x z              -0.0267
c$$$            C4=0 ! b x t + t x b  |  m x z + z x m       0.
c$$$            C5=0 ! n x t + t x n  |  n x z + z x n       0.
c          -----------------------------------------------------------


c          Convention following Groger's formulation
c          -----------------------------------------------------------
c               c1   !  n` x  b
c               c2   !  n  x  n  -  t  x  t
c               c3   !  n` x  n` -  t` x  t`
c               c4   !  b  x  b
c               c5   !  -
c               c6   !  c6/3
c

C *** READS INFORMATION ABOUT SLIP AND TWINNING SYSTEMS
      READ(UR_SX,'(A)') PROSA
      !write(*,*)'prosa:',trim(prosa)
c      write(*,*)'isx_ver:',isx_ver
      if (isx_ver.eq.3.or.isx_ver.eq.4.or.isx_ver.eq.5) then
         READ(UR_SX,*) C1,C2,C3,C4,C5,C6
c$$$         write(*,*)' C1,C2,C3,C4,C5,C6'
c$$$         write(*,'(6f7.3)')C1,C2,C3,C4,C5,C6
      else
         c1=0.
         c2=0.
         c3=0.
         c4=0.
         c5=0.
         c6=0.
      endif
      READ(UR_SX,*)     NMODESX !! the number of 'written' modes
      READ(UR_SX,*)     NMODES(iph) !! the number of actually used modes
c     youngung 2020 03 04
      mode(:)=0 ! possibly a bug: mode should be re-initialized?
      READ(UR_SX,*)     (MODE(I),I=1,NMODES(iph)) !! pointer to NMODES

      IF(NMODES(iph).GT.NMODMX) THEN
         WRITE(*,'('' NMODES IN PHASE'',I3,'' IS'',I3)')
     $        IPH,NMODES(iph)
         WRITE(*,*)"CHANGE PARAMETER NMODMX passed to subr.",
     $        " data_crystal_fe"
         write(*,*)'iph:',iph
         write(*,*)'nmodes(iph):',nmodes(iph)
         write(*,*)'nmodmx:',nmodmx
         call xit()
      ENDIF

      ICS=1       ! ICS=1 --> centro_symmetric SCYS
      nslmod(iph)=0
      ntwmod(iph)=0
      nsyst(iph) =0
      nslsys(iph)=0
      ntwsys(iph)=0
      MCOUNT=1    ! counter for the number of modes

C *** READS DEFORMATION MODES AND ASSOCIATED PARAMETERS FROM FILECRYS

      DO 100 MLOOP=1,NMODESX !! loop over 'written' modes
         READ(UR_SX,'(a)') PROSA
         READ(UR_SX,*)     MODEX,NSMX,ISENSEX,ITWTYPEX
         IF(MODEX.NE.MLOOP) THEN
            WRITE(*,*)'WARNING !!!'
            WRITE(*,*)'MODE NUMBERS MUST BE SEQUENTIAL IN CRYSTAL FILE'
            call xit()
         ENDIF
         IF(ITWTYPEX.NE.0) THEN
            IF(ISENSEX.EQ.1) THEN
               WRITE(*,*)'WARNING: ISENSEX HAS TO BE 0 FOR TWIN SYSTEMS'
               call xit()
            ENDIF
            READ(UR_SX,*) TWSHX
         ENDIF

C *** SKIPS MODE IF IT IS NOT IN THE LIST
         IF(MLOOP.NE.MODE(MCOUNT)) THEN
            DO IS=1,NSMX
               READ(UR_SX,*)
            ENDDO
            GO TO 100
         ENDIF

         IF(ISENSEX.EQ.0) ICS=0 ! ICS=0 --> non-centro-sym SCYS

         NSM(MCOUNT,iph)=NSMX

C *** VERIFICATION OF TWINNING DATA TO BE SURE PROGRAM WILL RUN PROPERLY

         IF(ITWTYPEX.EQ.0 .AND. NTWMOD(iph).NE.0) THEN
            WRITE(*,'(''  SLIP MODE'',I4,''  IN PHASE'',I4)') MLOOP,IPH
            WRITE(*,*) '  SLIP MODES MUST PRECEDE TWIN MODES'
            WRITE(*,*) '  -->   REORDER CRYSTAL FILE'
            call xit()
         ENDIF

         NSYST(iph)=NSYST(iph)+NSMX
         IF(ITWTYPEX.EQ.0) THEN
            NSLMOD(iph)=NSLMOD(iph)+1
            NSLSYS(iph)=NSLSYS(iph)+NSMX
         ELSE IF(ITWTYPEX.NE.0) THEN
            NTWMOD(iph)=NTWMOD(iph)+1
            NTWSYS(iph)=NTWSYS(iph)+NSMX
            TWSH(NTWMOD(iph),iph) =TWSHX
         ENDIF

         IF(NSYST(iph).GT.NSYSMX) THEN
            WRITE(*,'(''NSYST IN PHASE'',I3,'' IS'',I3)') IPH,
     $           NSYST(iph)
            write(*,*)'nsysmx:',nsysmx
            WRITE(*,*)"**CHANGE PARAMETER NSYSMX passed IN subr.",
     $           " data_crystal_fe"
            call xit()
         ENDIF

         IF(NTWMOD(iph).GT.NTWMMX) THEN
            WRITE(*,'('' NTWMOD IN PHASE'',I3,'' IS'',I3)')
     $           IPH,NTWMOD(iph)
            WRITE(*,*)"CHANGE PARAMETER NTWMMX passed IN subr.",
     $           " data_crystal_fe"
            call xit()
         ENDIF

         IF(NTWSYS(iph).GT.NTWSMX) THEN
            WRITE(*,'('' NTWSYS IN PHASE'',I3,'' IS'',I3)') IPH,
     $           NTWSYS(iph)
            WRITE(*,*)"CHANGE PARAMETER NTWSMX passed IN subr.",
     $           " data_crystal_fe"
            call xit()
         ENDIF

         ISBOT=NSYST(iph)-NSMX+1
         ISTOP=NSYST(iph)

         DO ISYS=ISBOT,ISTOP

            ISENSE(ISYS,iph)=ISENSEX
            ITWTYPE(ISYS,iph)=ITWTYPEX

C *** READS MILLER INDICES AND CALCULATES CARTESIAN COMPONENTS OF NORMAL
C *** (SN) AND SHEAR (SB) VECTORS OF SLIP OR TWIN SYSTEMS.
C *** CALCULATES SCHMID TENSORS IN CRYSTAL AXES 'SCHCA' FOR EACH SYSTEM.

            READ(UR_SX,*) (ISN(I),I=1,NIND),(ISB(I),I=1,NIND)
            CALL CRYSTAL_SYMMETRY(2,UR_SX,ICRYSYM(iph),
     $           ISN,SN,SNEQ,ISB,SB,
     $           NPOLES)

            call crossproduct(sn,sb,st)

            snbt(isys,1,:)=sn(:)
            snbt(isys,2,:)=sb(:)
            snbt(isys,3,:)=st(:)

            PROD=SN(1)*SB(1)+SN(2)*SB(2)+SN(3)*SB(3)
            PROD=PROD+SN(1)*ST(1)+SN(2)*ST(2)+SN(3)*ST(3)
            PROD=PROD+SB(1)*ST(1)+SB(2)*ST(2)+SB(3)*ST(3)

            IF(ABS(PROD) .GE. 1.E-4) THEN
               WRITE(*,'('' SYSTEM IS NOT ORTHOGONAL !!'')')
               WRITE(*,'('' ISN='',4I7)') (ISN(J),J=1,NIND)
               WRITE(*,'('' ISB='',4I7)') (ISB(J),J=1,NIND)
               WRITE(*,'(''   N='',3F7.3)') (SN(J),J=1,3)
               WRITE(*,'(''   B='',3F7.3)') (SB(J),J=1,3)
               WRITE(*,'(''   t='',3F7.3)') (ST(J),J=1,3)
               call xit()
            ENDIF


            if (inonschlaw.eq.1) then
c           ** Groger's non-Schmid formulation
c              Calculate t` and n` for Groger's formulation
               wij(:,:)=0d0
               wij(1,2)=-sb(3)
               wij(1,3)= sb(2)
               wij(2,1)= sb(3)
               wij(2,3)=-sb(1)
               wij(3,1)=-sb(2)
               wij(3,2)= sb(1)
c            rij(:,:)=0d0
               rij(:,:)=0d0
               c60=cos(60.*pi/180d0)
               s60=sin(60.*pi/180d0)
               do i=1,3
               do j=1,3
                  rij(i,j)=c60*i/j*j/i+s60*wij(i,j)+(1-c60)*sb(i)*sb(j)
               enddo
               enddo
               call matvect_dot(rij,sn,nprime,3)
               call crossproduct(nprime,sb,tprime)
            endif

            DO I=1,3
               DNCA(I,ISYS,iph)=SN(I)
               DBCA(I,ISYS,iph)=SB(I)
               DO J=1,3
                  PS(I,J)=C0*(SB(I)*SN(J)+SB(J)*SN(I)) !! schmid tensor
c                 inonschlaw is Yalcinkaya et al.
                  if (inonschlaw.eq.0) then
                     PNS(i,j)=c1*sb(i)*sb(j)
     $                    +c2*sn(i)*sn(j)
     $                    +c3*st(i)*st(j)
     $                    +c4*(sb(i)*st(j)+sb(j)*st(i))
     $                    +c5*(sn(i)*st(j)+sn(j)*st(i))
                  elseif(inonschlaw.eq.1) then
                     PNS(i,j)=c1*nprime(i)*sb(j)
     $                    +c2*(sn(i)*sn(j)-st(i)*st(j))
     $                    +c3*(nprime(i)*nprime(j)-tprime(i)*tprime(j))
     $                    +c4*sb(i)*sb(j)
                  else
                     write(*,*)'** Error: unexpected inonschlaw '
                     write(*,*)'          given in data_crystal_fe.f'
                     call xit()
                  endif
                  AUX33(I,J)=PS(I,J)
                  BUX33(I,J)=PNS(I,J)
               ENDDO
            ENDDO

            aux5(:)=0d0
            CALL CHG_BASIS(AUX5,AUX33,AUX66,AUX3333,2,5)
            DO I=1,5
               schca(I,ISYS,iph)=AUX5(I)
            ENDDO
            bux6(:)=0d0
            CALL CHG_BASIS(BUX6,BUX33,AUX66,AUX3333,2,6)
            DO I=1,6
               nschca(I,ISYS,iph)=BUX6(I)
            ENDDO
         enddo                  ! END OF LOOP OVER DEFORMATION MODES

C *** USE INFO FROM LAST SYSTEM TO CALCULATE SHEAR MODULUS OF THE SLIP MODE
         aux66(:,:)=celccv(:,:) ! from voigt to cartesian
         CALL VOIGT(AUX6,AUX33,aux66,AUX3333,3)
         xMU_MODE(MCOUNT,iph)=0d0
         DO I=1,3
         DO J=1,3
         DO K=1,3
         DO L=1,3
            xMU_MODE(MCOUNT,iph)=xMU_MODE(MCOUNT,iph)+
     $           AUX33(I,J)*AUX3333(I,J,K,L)*AUX33(K,L)
         ENDDO
         ENDDO
         ENDDO
         ENDDO
         MCOUNT=MCOUNT+1
 100  CONTINUE       ! END OF LOOP OVER ALL MODES IN PHASE 'IPH'


c**   YJ experiments on non-Schmid effect on Neumann principle.
      if (.false.) then
         write(*,*)'** Check any unwanted bias in n,b,t vectors'
         write(*,'(2a8)')'isbot','istop'
         write(*,'(3i8)')isbot,istop,istop-isbot+1

         do i=1,3               ! n,b,t
!     write(*,*)'i-th vector (n,b,t),',i
            if (i.eq.1) write(*,*)'n'
            if (i.eq.2) write(*,*)'b'
            if (i.eq.3) write(*,*)'t'
            do j=1,3            ! x,y,z component
               dum=0.
               do isys=isbot,istop
                  dum=dum+snbt(isys,i,j)
               enddo
               write(*,'(a8,x)',advance='no')'sum of'
               if (j.eq.1) write(*,'(a2,x)',advance='no') 'x:'
               if (j.eq.2) write(*,'(a2,x)',advance='no') 'y:'
               if (j.eq.3) write(*,'(a2,x)',advance='no') 'z:'
               write(*,'(f9.6)')dum
            enddo
         enddo
      endif

C *** CHECKS WHETHER THE SINGLE CRYSTAL YIELD SURFACE IS OPEN

      DO ICOMP=1,5
         ICLOSEPOS=0
         ICLOSENEG=0
         DO NS=1,NSLSYS(iph)
            IF(ABS(schca(ICOMP,NS,iph)).GT.1.E-3) THEN
               ICLOSEPOS=1
               ICLOSENEG=1
            ENDIF
         ENDDO
         IF(NTWSYS(iph).NE.0) THEN
            DO NS=NSLSYS(iph)+1,NSYST(iph)
               IF(schca(ICOMP,NS,iph).GT. 1.E-3) ICLOSEPOS=1
               IF(schca(ICOMP,NS,iph).LT.-1.E-3) ICLOSENEG=1
            ENDDO
         ENDIF
         IF(ICLOSEPOS.NE.1 .OR. ICLOSENEG.NE.1) THEN
            WRITE(*,'('' WARNING ! THE SCYS IS OPEN FOR PHASE'',I5,
     $'' ALONG DIRECTION'',I5)')IPH,ICOMP
            call stop_debug(0,'stopped in <data_crystal> due to opened
     $ SCYS')
         ENDIF
      ENDDO

      read(UR_SX,'(a)')prosa
      READ(UR_SX,*) fstar

      RETURN
      END subroutine data_crystal_fe
