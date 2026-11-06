c----------------------------------------------------------------------c
c     This calling block of subroutine <initial_state_guess>  is shared
c     by Stand-alone and FE-interface versions.
c----------------------------------------------------------------------c
      call initial_state_guess(
c                  in
     $     mdim,iregime,strain_control,interaction,
     $     iDBARv,iSBARv,isense,NPH,nphmx,nsysmx,ngrmx,
     $     ngr,nsyst,
     $     nrs,dtime,gamd0s,crss,mgel,sch,nsch,c6,sg_old,
     $     wgt,dg_trans,selav,ioptfs,
c                  inout
     $     DBAR,sbar,dsbar,sbar_old,
c                  out
     $     sav,dav,DBAR0,tg,dge,mgtg,dg0,mbar,LBAR,
     $     dsg,sg,dsgtry,sgtry)
