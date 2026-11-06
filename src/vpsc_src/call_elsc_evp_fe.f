c----------------------------------------------------------------------c
c     This calling block of subroutine <initial_state_guess>  is shared
c     by Stand-alone and FE-interface versions.
c----------------------------------------------------------------------c

C       IOPTION=0 : CALCULATES UPPER BOUND ELASTIC AND THERMAL TENSORS
C                  'CELUB' AND 'ATHUB'
C       IOPTION=1 : CALCULATES SELF CONSISTENT ELASTIC AND THERMAL TENSORS
C                  'CELAV' AND 'ATHAV_X'
C       IOPTION=2 : CALCULATES STRESS IN GRAINS INDUCED BY ELASTIC LOAD
C                   AND/OR TEMPERATURE CHANGE
      call elsc_evp_fe(
c          in
     $     iopt,nphmx,NPH,athcc,celcc,selcc,ngr,ngrmx,ag,wgt,
     $     itmaxint,ishape,axisph,axisgr,
c          inout
     $     celav,ac_el,bc_el,
c          out
     $     athcs,celcs,mgel,selav)
