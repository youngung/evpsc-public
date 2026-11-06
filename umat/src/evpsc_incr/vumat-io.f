c      write(*,*)'** inside vumat-io'
c      noel=km+iElemNumStartLoc -1
      noel=km+intPt -1
      npt=km
      ipx         = (noel-1)*50+npt ! polycrystal ID / integration point ID
      uwdev       = 105+52*(ipx-1)
      uwdeb       = uwdev+1
      UR_INP = uwdeb+1         ! 'vpsc7.in'
      UR_SX      = UR_INP+1 ! 'single crytsal file'
      UR_TX      = UR_SX+1    ! filetexture
      UR_SHP      = UR_TX+1    ! fileaxes (3 seems to be reserved?)
      UR_STV     = UR_SHP + 1  ! stv files
      UR_FEINP      = UR_STV + 1 ! 'evpsc_fe.in'
      uwfe01      = UR_FEINP+1    ! texture output
      uwjac       = uwfe01 + 1  ! File jacobian
      uwvel       = uwjac + 1   ! File velocity gradient
      UR_DIF      = uwvel + 1   ! file for single diffraction
      UW_LOG     = UR_DIF + 1  ! RUN_LOG.OUT
c     End of FILE I/O UNIT declaration
      if (noel.eq.noel_start.and.km.eq.1)then
         idevwrite=.true.
         iverbose=1
      else
         idevwrite=.false.
         iverbose=0
      endif
      iverbose=0                ! overwrite to make vumat.f completely silent
      idevwrite=.false.         ! overwrite to write to files.
