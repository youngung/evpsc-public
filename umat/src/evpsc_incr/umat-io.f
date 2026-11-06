c     Begin to determine FILE I/O UNIT
      ipx         = (noel-1)*NIPMX+npt+kspt ! polycrystal ID / integration point ID
      uw_defgrad  = 101 + 100*(ipx-1) + layer
      uwdev       = uw_defgrad + 1
      uwdeb       = uwdev + 1
      UR_INP      = uwdeb + 1   ! 'vpsc7.in'
      UR_SX       = UR_INP + 1  ! 'single crytsal file'
      UR_TX       = UR_SX + 1   ! filetexture
      UR_SHP      = UR_TX + 1   ! fileaxes (3 seems to be reserved?)
      UR_STV      = UR_SHP + 1  ! stv files
      UR_FEINP    = UR_STV + 1  ! 'evpsc_fe.in'
      uwjac       = UR_FEINP + 1 ! File jacobian
      uwvel       = uwjac + 1   ! File velocity gradient
      UR_DIF      = uwvel + 1   ! file for single diffraction
      UR_ECTE     = UR_DIF + 1  ! file for temperature-dependent Emod/CTE
      UR_EBSD     = UR_ECTE + 1
      UW_LOG      = UR_EBSD + 1 ! RUN_LOG.OUT
      uwfe01      = UW_LOG+(kspt-1)*5 ! texture output
      UW_DIS      = uwfe01 + 1
      UW_TAU      = UW_DIS + 1
      UW_SCCN     = UW_TAU + 1  ! sc_convergence.out (not used in FE)

c     End of FILE I/O UNIT declaration

      if (ipx.lt.1) then
         write(*,*)'noel:',noel,'nipmx:',nipmx,'npt:',npt
         write(*,*)'ipx is less than 1',ipx
         call xit()
      endif
      if (uwdev.lt.1) then
         write(*,*)'noel:',noel,'nipmx:',nipmx,'npt:',npt
         write(*,*)'uwdev is less than 1',uwdev
         call xit()
      endif
      if (uwdeb.lt.1) then
         write(*,*)'noel:',noel,'nipmx:',nipmx,'npt:',npt
         write(*,*)'uwdeb is less than 1',uwdeb
         call xit()
      endif
      if (UR_INP.lt.1) then
         write(*,*)'noel:',noel,'nipmx:',nipmx,'npt:',npt
         write(*,*)'UR_INP is less than 1',UR_INP
         call xit()
      endif
      if (UR_SX.lt.1) then
         write(*,*)'noel:',noel,'nipmx:',nipmx,'npt:',npt
         write(*,*)'UR_SX is less than 1',UR_SX
         call xit()
      endif
      if (UR_TX.lt.1) then
         write(*,*)'noel:',noel,'nipmx:',nipmx,'npt:',npt
         write(*,*)'UR_TX is less than 1',UR_TX
         call xit()
      endif
      if (UR_SHP.lt.1) then
         write(*,*)'noel:',noel,'nipmx:',nipmx,'npt:',npt
         write(*,*)'UR_SHP is less than 1',UR_SHP
         call xit()
      endif
      if (UR_STV.lt.1) then
         write(*,*)'noel:',noel,'nipmx:',nipmx,'npt:',npt
         write(*,*)'UR_STV is less than 1',UR_STV
         call xit()
      endif
      if (UR_FEINP.lt.1) then
         write(*,*)'noel:',noel,'nipmx:',nipmx,'npt:',npt
         write(*,*)'UR_FEINP is less than 1',UR_FEINP
         call xit()
      endif
      if (uwjac.lt.1) then
         write(*,*)'noel:',noel,'nipmx:',nipmx,'npt:',npt
         write(*,*)'uwjac is less than 1',uwjac
         call xit()
      endif
      if (uwvel.lt.1) then
         write(*,*)'noel:',noel,'nipmx:',nipmx,'npt:',npt
         write(*,*)'uwvel is less than 1',uwvel
         call xit()
      endif
      if (UR_DIF.lt.1) then
         write(*,*)'noel:',noel,'nipmx:',nipmx,'npt:',npt
         write(*,*)'UR_DIF is less than 1',UR_DIF
         call xit()
      endif
      if (UR_EBSD.lt.1) then
         write(*,*)'noel:',noel,'nipmx:',nipmx,'npt:',npt
         write(*,*)'UR_EBSD is less than 1',UR_EBSD
         call xit()
      endif
      if (UW_LOG.lt.1) then
         write(*,*)'noel:',noel,'nipmx:',nipmx,'npt:',npt
         write(*,*)'UW_LOG is less than 1',UW_LOG
         call xit()
      endif
      if (uwfe01.lt.1) then
         write(*,*)'noel:',noel,'nipmx:',nipmx,'npt:',npt
         write(*,*)'uwfe01 is less than 1',uwfe01
         call xit()
      endif



c      if ((noel.eq.noel_top).and.npt.eq.1)
      if ((noel.eq.noel_top.or.noel.eq.noel_bot).and.npt.eq.1)
c      if ((noel.eq.441.or.noel.eq.noel_bot).and.npt.eq.1)
     $     then
         idevwrite=.true.
         iverbose=1
      else
         idevwrite=.false.
         iverbose=0
      endif
c     if shell, restricted to the first cross section point.
      if (ndi.eq.2.and.nshr.eq.1.and.iverbose.eq.1) then
         if (kspt.ne.1) then
            iverbose=0
            idevwrite=.false.
         endif
      endif
