c     This block of statements is called within <umat.f>

      if (idevwrite) then
         inquire(file=OUTDIR(1:LENOUTDIR)//'/evpsc_fe.in',
     $        exist=inquiry_exist)
         ninquiry = 0
         if (inquiry_exist) then
            open(UR_FEINP,file=OUTDIR(1:LENOUTDIR)//'/evpsc_fe.in',
     $           status='old')
            read(UR_FEINP,'(a)') prosa
            read(UR_FEINP,'(a)') path_inquiry ! to which files are saved
            read(UR_FEINP,*) ninquiry

            if (ninquiry.gt.inqmx) then
               write(*,*)'** Error: ninquiry exceeds inqmx; need a fix'
               write(*,'(a25,i5.5)')'ninquiry:',ninquiry
               write(*,'(a25,i5.5)')'inqmx',inqmx
               call xit()       ! for UMAT
c	       call xplb_exit() ! for VUMAT
            endif

            do i=1,ninquiry
               read(UR_FEINP,*) inquiry(i,1),inquiry(i,2)
            enddo
            close(UR_FEINP)
         endif



         open(uw_defgrad,file=OUTDIR(1:LENOUTDIR)//'/degrad.out',
     $        status='unknown')

      endif
