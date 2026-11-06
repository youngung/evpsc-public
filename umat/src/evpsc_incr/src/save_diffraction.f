c----------------------------------------------------------------------c
c     Save those state variables required for diffraction analysis, or vice versa
      subroutine save_diffraction(io,iopt,
     $     NPH,ngr,ag,sg,mgel,wgt,axisph,eulerph,epsvm)
      implicit none
      include "../../../../src/vpsc_src/incr_evpsc_fe.dim"
      integer io,iopt,iphel,i,j,kkk
c        NPH
c        ngr
c        AG
c        sg
c        wgt
c        mgel
      if (iopt.eq.0) then
         write(io,*) NPH
         write(io,*) (ngr(iphel),iphel=0,NPH)
         do iphel=1,NPH
            do kkk=ngr(iphel-1)+1,ngr(iphel)
               write(io,*) epsvm,
     $              ((ag(i,j,kkk),j=1,3),i=1,3),
     $              (sg(i,kkk),i=1,6),
     $              ((mgel(i,j,kkk),j=1,6),i=1,6),
     $              wgt(kkk)
            enddo
         enddo
      elseif (iopt.eq.1) then
         read(io,*) NPH
         !write(*,*)'NPH inside save_diffraction:',NPH
         read(io,*) (ngr(iphel),iphel=0,NPH)
         do iphel=1,NPH
            do kkk=ngr(iphel-1)+1,ngr(iphel)
               read(io,*) epsvm,
     $              ((ag(i,j,kkk),j=1,3),i=1,3),
     $              (sg(i,kkk),i=1,6),
     $              ((mgel(i,j,kkk),j=1,6),i=1,6),
     $              wgt(kkk)
            enddo
         enddo
      endif
      return
      end subroutine save_diffraction
