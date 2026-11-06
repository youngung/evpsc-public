      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     $     RPL,DDSDDT,DRPLDE,DRPLDT,
     $     STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     $     NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     $     CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      implicit none
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     $     DDSDDE(NTENS,NTENS),
     $     DDSDDT(NTENS),DRPLDE(NTENS),
     $     STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     $     PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
      DIMENSION DSTRES(6),D(3,3)

      if (noel.eq.1.and.npt.eq.1)then
         call w_chrc(0,'kstep,kinc')
         call w_ival(0,'kstep',kstep)
      endif

C
C  EVALUATE NEW STRESS TENSOR
C
      EV = 0.
      DEV = 0.
      DO K1=1,NDI
         EV = EV + STRAN(K1)
         DEV = DEV + DSTRAN(K1)
      END DO
C
      TERM1 = .5*DTIME + PROPS(5)
      TERM1I = 1./TERM1
      TERM2 = (.5*DTIME*PROPS(1)+PROPS(3))*TERM1I*DEV
      TERM3 = (DTIME*PROPS(2)+2.*PROPS(4))*TERM1I
C
      DO K1=1,NDI
         DSTRES(K1) = TERM2+TERM3*DSTRAN(K1)
     $        +DTIME*TERM1I*(PROPS(1)*EV
     $        +2.*PROPS(2)*STRAN(K1)-STRESS(K1))
         STRESS(K1) = STRESS(K1) + DSTRES(K1)
      END DO
C
      TERM2 = (.5*DTIME*PROPS(2) + PROPS(4))*TERM1I
      I1 = NDI
      DO K1=1,NSHR
         I1 = I1+1
         DSTRES(I1) = TERM2*DSTRAN(I1)+
     $        DTIME*TERM1I*(PROPS(2)*STRAN(I1)-STRESS(I1))
         STRESS(I1) = STRESS(I1)+DSTRES(I1)
      END DO
C
C  CREATE NEW JACOBIAN
C
      TERM2 = (DTIME*(.5*PROPS(1)+PROPS(2))+PROPS(3)+
     $     2.*PROPS(4))*TERM1I
      TERM3 = (.5*DTIME*PROPS(1)+PROPS(3))*TERM1I
      DO K1=1,NTENS
         DO K2=1,NTENS
            DDSDDE(K2,K1) = 0.
         END DO
      END DO
C
      DO K1=1,NDI
         DDSDDE(K1,K1) = TERM2
      END DO
C
      DO K1=2,NDI
         N2 = K11
         DO K2=1,N2
            DDSDDE(K2,K1) = TERM3
            DDSDDE(K1,K2) = TERM3
         END DO
      END DO
      TERM2 = (.5*DTIME*PROPS(2)+PROPS(4))*TERM1I
      I1 = NDI
      DO K1=1,NSHR
         I1 = I1+1
         DDSDDE(I1,I1) = TERM2
      END DO
C
C  TOTAL CHANGE IN SPECIFIC ENERGY
C
      TDE = 0.
      DO K1=1,NTENS
         TDE = TDE + (STRESS(K1)-.5*DSTRES(K1))*DSTRAN(K1)
      END DO
C
C  CHANGE IN SPECIFIC ELASTIC STRAIN ENERGY
C
      TERM1 = PROPS(1) + 2.*PROPS(2)
      DO  K1=1,NDI
         D(K1,K1) = TERM1
      END DO
      DO K1=2,NDI
         N2 = K1-1
         DO K2=1,N2
            D(K1,K2) = PROPS(1)
            D(K2,K1) = PROPS(1)
         END DO
      END DO
      DEE = 0.
      DO K1=1,NDI
         TERM1 = 0.
         TERM2 = 0.
         DO K2=1,NDI
            TERM1 = TERM1 + D(K1,K2)*STRAN(K2)
            TERM2 = TERM2 + D(K1,K2)*DSTRAN(K2)
         END DO
         DEE = DEE + (TERM1+.5*TERM2)*DSTRAN(K1)
      END DO
      I1 = NDI
      DO K1=1,NSHR
         I1 = I1+1
         DEE = DEE + PROPS(2)*(STRAN(I1)+.5*DSTRAN(I1))*DSTRAN(I1)
      END DO
      SSE = SSE + DEE
      SCD = SCD + TDE  DEE
      RETURN
      END



c----------------------------------------------------------------------c
c     subroutines developed only for UMAT interface
      include "src/lib_fe.f"        ! some misc. subroutines
      include "src/lib_write_fe.f"  ! output subroutines for diagnosis
      include "src/voigt_abq.f"     ! to convert stress, strain, compliance (and stiffness) for Abaqus
      include "src/statev_io.f"     ! Subroutine to store and restore state variables
      include "src/uvarm.f"         ! User defined output accessible from Abaqus CAE (placeholder)
      include "src/save_diffraction.f"
      include "src/calc_devm.f"
      include "sdvini.f"
      include "src/statev_beta.f"
c     external
      include "src/uexternaldb.f"
c     subroutine just borrowed from vpsc
c      include "src/rodrigues.f" ! rodrigues
c      include "src/twin_orientation_fe.f"
      include "src/calc_spin.f"
      include "src/loading_condition_fe.f"
c     subroutines shared with stand-alone (incr. E)VPSC
      include "../../../src/vpsc_src/evpsc_input_fe.f"
      include "../../../src/vpsc_src/eshelby_rot.f"
      include "../../../src/vpsc_src/elsc_evp_fe.f"
      include "../../../src/vpsc_src/chg_basis.f"
      include "../../../src/vpsc_src/eshelby.f"
      include "../../../src/vpsc_src/lu_libs.f"
      include "../../../src/vpsc_src/funcs.f"
      include "../../../src/vpsc_src/euler.f"
      include "../../../src/vpsc_src/crystal_symmetry.f"
      include "../../../src/vpsc_src/voigt.f"
      include "../../../src/vpsc_src/neighbours_fe.f"
      include "../../../src/vpsc_src/library7.sub"
      include "../../../src/vpsc_src/initial_state_guess_el_evp.f"
      include "../../../src/vpsc_src/solve_evp_bc.f"
      include "../../../src/vpsc_src/initial_state_guess_evp.f"
      include "../../../src/vpsc_src/grain_rate_and_moduli_evp.f"
      include "../../../src/vpsc_src/evpsc.f"
      include "../../../src/vpsc_src/grain_stat.f"
      include "../../../src/vpsc_src/evpsc_implicit.f"
      include "../../../src/vpsc_src/grain_stress_evp.f"
      include "../../../src/vpsc_src/newton_raphson_evp.f"
      include "../../../src/vpsc_src/calc_evptilde.f"
      include "../../../src/vpsc_src/update_macro_property.f"
      include "../../../src/vpsc_src/update_schmid_fe.f"
      include "../../../src/vpsc_src/rodrigues.f"
      include "../../../src/vpsc_src/twin_orientation.f"
      include "../../../src/vpsc_src/update_shape_fe.f"
      include "../../../src/vpsc_src/data_crystal_fe.f"
      include "../../../src/vpsc_src/update_crss_voce_fe.f"
      include "../../../src/vpsc_src/update_crss_rgvb_fe.f"
      include "../../../src/vpsc_src/data_grain_fe.f"
      include "../../../src/vpsc_src/scale_gamd0s_evp_fe.f"
      include "../../../src/vpsc_src/update_fij_fe.f"
      include "../../../src/vpsc_src/update_orientation_fe.f"
      include "../../../src/vpsc_src/update_twinning_ptr_fe.f"
      include "../../../src/vpsc_src/update_twinning_vft_fe.f"
      include "../../../src/vpsc_src/stat_twinning_fe.f"
      include "../../../src/vpsc_src/write_texture_fe.f"
      include "../../../src/vpsc_src/stat_shear_activity.f"
      include "../../../src/vpsc_src/diff_planes.f"
c      include "../../../src/vpsc_src/update_crss_dd_fe.f"
      include "../../../src/vpsc_src/cod2grains.f"
      include "../../../src/vpsc_src/update_crss_dd_fe_new.f"
      include "../../../src/vpsc_src/thermo_emod_cte.f"
