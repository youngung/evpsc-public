      if (ninquiry.gt.0) then
         do i=1,ninquiry
            if (inquiry(i,1).eq.noel.and.inquiry(i,2).eq.npt) then
               write(prosa,'(i5.5)')kount_umat_call
               write(fnout_stv,*)trim(path_inquiry)//'/'//cnoel
     $              //'_'//cnpt//'_'//trim(adjustl(prosa))//'.stv'
               open(UR_STV,file=trim(fnout_stv),status='unknown')
               call save_diffraction(UR_STV,0,NPH,ngr,ag,sg,
     $              mgel,wgt,axisph,eulerph,epsvm)
               close(UR_STV)
            endif
         enddo
      endif
