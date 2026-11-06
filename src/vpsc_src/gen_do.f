c     generate discrete orientation
      program main
      implicit none
      include 'incr_evpsc_fe.dim'
      character(len=256)filetext

      integer igb,igt

      !read(*,'(a)')filetext
      filetext='matData/vpscData/EDDQ/eddq.odf'

      OPEN(action='read',UNIT=UR_TX,FILE=FILETEXT,STATUS='OLD')
c     Caution!! number of grains for this phase is
c     hardwired as ngrmx for now. This works
c     for a single phase material, but cetainly
c     requires improvement... (2021-Dec 13 YJ)

      igb=1
      igt=ngrmx
      call data_grain_odf(
c     in
     $     0,1,ngrmx,UR_TX,UR_SHP,ngrmx,nphmx,1,igb,igt,wph,
     $     ishape,eulerph,axisph,
c     out
     $     wgt,ag,fijph,fijgr)
      close(UR_TX)

      end program
