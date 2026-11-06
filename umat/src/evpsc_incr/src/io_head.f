      integer NIPMX,NELMX       ! No of maximum integration points and elements
      parameter(NIPMX=8,NELMX=15750)
      integer ipx,UR_INP,
     $     UR_FEINP,uwfe01,UW_SCCN,uw_defgrad
      integer NUVARMX
      parameter(nuvarmx=30)
      real*8 uvarm_custom(nelmx,nipmx,nuvarmx)
      common /kuvarm/uvarm_custom
