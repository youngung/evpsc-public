C **************************************************************************
C    SUBROUTINE STAT_TWINNING_fe   --->   VERSION OF 21/APR/2020
C
C --> TWINNING STATISTICS: SCHMID FACTOR DISTRIBUTION & VARIANT DISTRIBUTION
C **************************************************************************
C    SCH_STAT(ISLOT,ITWMODE,IPH): FREQUENCY DISTRIBUTION OF TWIN MODE
C                    'ITWMODE' OF PHASE 'IPH' AT SCHMID FACTOR BIN 'ISLOT'
C    VAR_STAT(IV,ITWMODE,IPH): VARIANT FREQUENCY OF TWIN MODE
C                    'ITWMODE' OF PHASE 'IPH' CORRESPONDS TO VARIANT 'IV'
C **************************************************************************
C     WARNING : DIMENSIONING ARRAYS TO '6' IMPLICITLY ASSUMES HCP TWINS
C **************************************************************************

      SUBROUTINE STAT_TWINNING_fe(
     $     KTW, ITWMODE, KKK,iph,nphmx,
     $     nsysmx,ntwmmx,nslmod,nmodmx,nsm,nslsys,dnca,dbca,
     $     ag,scauchy,wgt,
c         inout
     $     sch_stat,var_stat,var_vf
     $     )
      implicit none

      integer, intent(in)::nsysmx,ntwmmx,nslmod,nmodmx,nphmx,
     $     nsm(nmodmx,nphmx),nslsys(nphmx),iph
      integer, intent(inout)::sch_stat(-5:15,ntwmmx,nphmx),
     $     var_stat(6,ntwmmx,nphmx)
      real*8,intent(in)::dnca(3,nsysmx,nphmx),dbca(3,nsysmx,nphmx),
     $     ag(3,3),
     $     scauchy(3,3),wgt
      real*8,intent(inout)::var_vf(6,ntwmmx,nphmx)
      real*8 TWB(3),TWN(3),TWBX(3),TWNX(3),TRSS(6,3),AUX6(6),
     $     schx,rssmax,vnorm

      INTEGER ITWMODE,its,ktmode,i,j,kkk,islot,nsmx,kts,ktw,k

      kkk=kkk                   ! to fool the compiler
C *** CALCULATE AND ASSIGNS SCHMID FACTOR OF PTS TO 0.05 WIDE BINS FOR STATISTICS

      TWN(:) =DNCA(:,KTW,iph)
      TWB(:) =DBCA(:,KTW,iph)
      TWNX(:)=0.
      TWBX(:)=0.
      DO I=1,3
        DO J=1,3
          TWNX(I)=TWNX(I)+AG(I,J)*TWN(J)
          TWBX(I)=TWBX(I)+AG(I,J)*TWB(J)
        ENDDO
      ENDDO

C ---> TEMPORARY FIX WHEN BC'S DO NOT GIVE SCAUCHY
c      IF(iDBARv(1)*iDBARv(2)*iDBARv(3).EQ.1) SCAUCHY(:,:)=SDEVIAT(:,:)

      SCHX=0.
      DO I=1,3
        DO J=1,3
c     generalized Schmid factor definition
          SCHX=SCHX+TWNX(I)*TWBX(J)*SCAUCHY(I,J)
        ENDDO
      ENDDO
      SCHX=SCHX/VNORM(SCAUCHY,9)

      ISLOT=INT(SCHX/0.05) + 1         ! islot=1,2,...,10
      IF(ISLOT.GT.10 .OR. ISLOT.LT.1) THEN
c         WRITE(*,'('' SCAUCHY'',9F8.3)') SCAUCHY
c         WRITE(*,'('' GRAIN'',I5,'' SCHMID FACTOR'',F8.3,'' ISLOT='',I3,
c     $'' EXCEEDS DIMENSION IN SCH_STAT ARRAY'')') KKK,SCHX,ISLOT
         IF(ISLOT.GT.10) ISLOT=10 ! QUICK & DIRTY FIX FOR THE TIME BEING
         IF(ISLOT.LT. 1) ISLOT= 1 ! QUICK & DIRTY FIX FOR THE TIME BEING
      ENDIF
      SCH_STAT(ISLOT,ITWMODE,iph)=SCH_STAT(ISLOT,ITWMODE,iph) + 1

C *** TWIN VARIANT IDENTIFICATION FOR STATISTICS
        ITS=NSLSYS(iph) + (ITWMODE-1)*6
        KTMODE=ITWMODE+NSLMOD
            KTS = 0
        DO NSMX=1,NSM(KTMODE,iph)
            KTS=KTS+1
            ITS=ITS+1             ! absolute counter over all systems

            TWN(:) =DNCA(:,ITS,iph)
            TWB(:) =DBCA(:,ITS,iph)
            TWNX(:)=0.
            TWBX(:)=0.
            DO I=1,3
              DO J=1,3
                TWNX(I)=TWNX(I)+AG(I,J)*TWN(J)
                TWBX(I)=TWBX(I)+AG(I,J)*TWB(J)
              ENDDO
            ENDDO

            TRSS(KTS,:)=0        ! 'KTS' is relative
            TRSS(KTS,1)=ITS      ! 'ITS' is absolute
            DO I=1,3
              DO J=1,3
                TRSS(KTS,2)=TRSS(KTS,2)+TWNX(I)*TWBX(J)*SCAUCHY(I,J)
              ENDDO
            ENDDO
        ENDDO

            AUX6(:) = TRSS(:, 2)
            DO I = 1, 6
              RSSMAX = -1E+10
              DO J = 1, 6
              IF (RSSMAX .LT. AUX6(J)) THEN
                 K = J
                 RSSMAX = AUX6(J)
              ENDIF
           ENDDO
           AUX6(K) = -1E+5
           TRSS(K,3) = I
           IF (KTW .EQ. INT(TRSS(K,1))) THEN
              VAR_STAT(I,ITWMODE,iph) = VAR_STAT(I,ITWMODE,iph) + 1
              VAR_VF(I,ITWMODE,iph) = VAR_VF(I,ITWMODE,iph) + WGT
           ENDIF
          ENDDO

      RETURN
      END
