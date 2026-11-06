      if (ninquiry.gt.0) then
         do i=1,ninquiry
            if (inquiry(i,1).eq.noel.and.inquiry(i,2).eq.npt) then
               write(prosa,'(i8.8)')kount_umat_call
               prosa=adjustl(prosa)
               write(fnout_stv,*)trim(path_inquiry)//'/'//cnoel
     $              //'_'//cnpt//'_'//trim(prosa)//'.stv'
               open(UR_STV,file=trim(fnout_stv),status='unknown')
               call save_diffraction(UR_STV,0,NPH,ngr,ag,sg,
     $              mgel,wgt,axisph,eulerph,epsvm)
               close(UR_STV)

c              Write to UWVEL ---------
               call chg_basis(DBAR,aux33,aux66,aux3333,1,6)
               dsim(:,:)=aux33(:,:)
c              -- Need to check the below definition on velocity gradient.
               do k=1,3
                  do j=1,3
                     LIJBARc(k,j) = dsim(k,j)+WBARc(k,j)
                  enddo
               enddo
               write(uwvel,'(11e20.11)')
     $              totalTime,dtime,((LIJBARc(k,j),j=1,3),k=1,3)
c              End of writing to UWVEL ---------

c        dev.log
c               if (ielem_type.eq.0) then
c                  i12_comp=4
c               elseif (ielem_type.eq.1) then
c                  i12_comp=3
c               endif
c               include'src/write_inquiry_vumat_call.f'
            endif
         enddo
      endif
