c----------------------------------------------------------------------c
c     Subroutine to unload the polycrystal: impose 'stress rate'
c     by giving an incremental stress (dsbar) and incremental time.
c     The incremental time is hardwired at the moment (as 1e-3)
c     Whatever the last stress, we divide the stress by (nstep+1), where
c     nstep should be given when the function is called.
c     Youngung Jeong
c----------------------------------------------------------------------c
      SUBROUTINE UNLOAD_CONDITION(
     $     iopt,time_unload,tincr,ILBAR,nstep,
     $     strain_control,
     $     LIJBARc,iSBARv,iDBARv,dsim,scauchy,sbar,dsbar,DBAR)
      implicit none

      INTEGER,intent(in)::iopt,nstep
      integer,intent(inout)::ILBAR(3,3),iSBARv(6),iDBARv(6)
      integer,intent(out)::strain_control
      real*8,intent(out)::tincr,dsbar(6),DBAR(6)
      real*8,intent(inout)::LIJBARc(3,3),dsim(3,3),scauchy(3,3)
      real*8,intent(in)::time_unload,sbar(6)
c     --------------
c     locals
      real*8 aux33(3,3),AUX3333(3,3,3,3),aux66(6,6),sbardot(6)
      real*8 dscauchy(3,3)
      save dscauchy,sbardot ! defined with iopt.eq.0 used when iopt.eq.1

      if (iopt.eq.0) then
         ILBAR(:,:)=0
         ILBAR(1,2)=1
         ILBAR(1,3)=1
         ILBAR(2,3)=1
         strain_control=0 ! fully stress-controlled unloading
         iSBARv(:)=1
         iDBARv(:)=0
         LIJBARc(:,:)=0.
         dsim(:,:)=0.
         tincr=time_unload/(nstep)
c        the constant stress rate for the entire unloading process is
         sbardot(:)=-sbar(:)/time_unload
      elseif (iopt.eq.1) then
         call chg_basis(sbar,scauchy,aux66,aux3333,1,6)
         call chg_basis(sbardot,aux33,aux66,aux3333,1,6)
         dscauchy(:,:)=aux33(:,:)*tincr
         scauchy(:,:)=scauchy(:,:)+dscauchy(:,:)
c     next trial differential in stress increment (dsbar)
         call chg_basis(dsbar,dscauchy,aux66,aux3333,2,6)
         call chg_basis(sbar,scauchy,aux66,aux3333,2,6)
         CALL CHG_BASIS(DBAR,DSIM,AUX66,AUX3333,2,6)
      endif

      RETURN
      END SUBROUTINE UNLOAD_CONDITION

c----------------------------------------------------------------------c
c     stress hold condition
      subroutine stresshold_condition(
     $     iopt,time_unload,tincr,ILBAR,nstep,
     $     strain_control,
     $     LIJBARc,iSBARv,iDBARv,dsim,scauchy,sbar,dsbar,DBAR)
      implicit none

      INTEGER,intent(in)::iopt,nstep
      integer,intent(inout)::ILBAR(3,3),iSBARv(6),iDBARv(6)
      integer,intent(out)::strain_control
      real*8,intent(out)::tincr,dsbar(6),DBAR(6)
      real*8,intent(inout)::LIJBARc(3,3),dsim(3,3),scauchy(3,3)
      real*8,intent(in)::time_unload,sbar(6)
c     --------------
c     locals
      real*8 aux33(3,3),AUX3333(3,3,3,3),aux66(6,6),sbardot(6)
      real*8 dscauchy(3,3)
      save dscauchy,sbardot ! defined with iopt.eq.0 used when iopt.eq.1

      if (iopt.eq.0) then
         ILBAR(:,:)=0
         ILBAR(1,2)=1
         ILBAR(1,3)=1
         ILBAR(2,3)=1
         strain_control=0 ! fully stress-controlled unloading
         iSBARv(:)=1
         iDBARv(:)=0
         LIJBARc(:,:)=0.
         dsim(:,:)=0.
         tincr=time_unload/(nstep)
c        the constant stress rate for the entire unloading process is
         sbardot(:)=0.
      elseif (iopt.eq.1) then
         call chg_basis(sbar,scauchy,aux66,aux3333,1,6)
         call chg_basis(sbardot,aux33,aux66,aux3333,1,6)
         dscauchy(:,:)=aux33(:,:)*tincr
         scauchy(:,:)=scauchy(:,:)+dscauchy(:,:)
c        next trial stress increment (dsbar)
         call chg_basis(dsbar,dscauchy,aux66,aux3333,2,6)
         call chg_basis(sbar,scauchy,aux66,aux3333,2,6)
         CALL CHG_BASIS(DBAR,DSIM,AUX66,AUX3333,2,6)
      endif

      RETURN
      END SUBROUTINE

c----------------------------------------------------------------------c
