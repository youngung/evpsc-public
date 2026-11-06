c     subroutine var_dsbar
c     load the history of incremental stress from a file (filehist)
c     Impose the non-monotonic stress history to the aggregate.
      subroutine var_dsbar(
c          intent(in)
     $     ioption,

c          intent(inout)
     $     temperat,
c          intent(out)
     $     tincr,dsbar,dsim,
     $     iDBARv,iSBARv,ILBAR,ictrl,strain_control,nsteps,tempdot
     $     )
      implicit none
      integer, intent(in)::ioption
      real*8, intent(out):: tincr,dsbar(6),dsim(3,3),tempdot
      integer, intent(out):: iDBARv(6),nsteps,ILBAR(3,3),ictrl
      real*8, intent(inout)::temperat
      real*8 AUX3333(3,3,3,3),aux66(6,6),scauchy(3,3)
      integer i,j,UR_HIS,iSBARv(6),strain_control
c      real*8 ctrlincr

c     locals
      real*8 dsbar33(3,3)
      integer iostat

      UR_HIS=99

      IF(IOPTION.EQ.0) THEN
         write(*,*)'UR_HIS:',UR_HIS
         READ(UR_HIS,*) NSTEPS,ICTRL,temperat
         READ(UR_HIS,*)
c         CTRLINCR=CTRLINCR      ! to fool compiler
         IF(ICTRL.NE.7) THEN
            WRITE(*,*)'VARIABLE VELOCITY GRADIENT OPTION ASSUMES THAT'
            WRITE(*,*)'THE VELOCITY GRADIENT AND A TIME INCREMENT ARE'
            WRITE(*,*)'IMPOSED AT EACH STEP --> SET ICTRL=7 !!'
            call xit()
         ENDIF

c***  General settings for fully imposed variable stress increment
         DO I=1,3
            DO J=1,3
               ILBAR(I,J)=0
               SCAUCHY(I,J)=0d0
            ENDDO
         ENDDO
         ILBAR(1,2)=1
         ILBAR(1,3)=1
         ILBAR(2,3)=1

         DO I=1,6
            iDBARv(I)=0
            iSBARv(I)=1
         ENDDO
      elseif(ioption.eq.1) then
         write(*,*)'** inside var_dsbar'
         read(UR_HIS,*,iostat=iostat) tincr,dsbar33(1,1),dsbar33(2,2),
     $        dsbar33(3,3),dsbar33(2,3),dsbar33(1,3),dsbar33(1,2),
     $        tempdot

         if (iostat.lt.0) then
            call stop_debug(0,'** end of file in var_dsbar')
         endif

         dsbar33(3,2)=dsbar33(2,3)
         dsbar33(3,1)=dsbar33(1,3)
         dsbar33(2,1)=dsbar33(1,2)

c        Cauchy stress
         scauchy(:,:)=dsbar33(:,:)+scauchy(:,:)
c        incremental stress (dsbar) in b-basis
         call chg_basis(dsbar,dsbar33,aux66,aux3333,2,6)
         dsim(:,:)=0d0
      endif

      strain_control=0

      return
      end subroutine var_dsbar
