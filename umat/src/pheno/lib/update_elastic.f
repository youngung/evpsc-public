c-----------------------------------------------------------------------
      subroutine update_elastic(idia,idiaw,iyld_law,ntens,nyldp,nstatv,
     $     ddsdde,cel,stran,stran_el_ns,stran_pl_ns,dstran,stress,
     $     eeq_ns,deq,yldp_ns,statev,kinc)
      include 'ABA_PARAM.INC'
      integer idia,iyld_law,nstatv,ntens,nyldp,kinc
      logical idiaw
      dimension ddsdde(ntens,ntens),cel(ntens,ntens),stran(ntens),
     $     stran_el_ns(0:1,ntens),stran_pl_ns(0:1,ntens),dstran(ntens),
     $     stress(ntens),yldp_ns(0:1,nyldp),statev(nstatv),eeq_ns(0:1)
      real*8 ddsdde,cel,stran,stran_el_ns,stran_pl_ns,dstran,
     $     stress,yldp_ns,statev,eeq_ns,deq,auxntens(ntens),
     $     buxntens(ntens),auxnyldp(nyldp)
      real*8 empa,gpa
      empa = 1d6
      gpa  = 1d9
      deq=0d0
      if (idiaw) then
         call w_chr(idia,'** inside update_elastic')
         call w_chr(idia,'* stran at step n')
         call w_dim(idia,stran,ntens,1d0/empa,.true.)
         call w_chr(idia,'* stress at step n')
         call w_dim(idia,stress,ntens,1d0/empa,.true.)
      endif

c$$$  1. Save jacobian as elastic moduli
      ddsdde(:,:) = Cel(:,:)
c$$$  2. Update strain.
      stran_el_ns(1,:) = stran_el_ns(0,:) + dstran(:)
c     call add_array(stran,dstran,ntens)
c$$$  3. Updates stress
      auxntens(:)=stran_el_ns(1,:)
      call mult_array(ddsdde,auxntens,ntens,stress)

      if (idiaw) then
         call w_chr( idia,'* stress at n+1')
         call w_dim( idia,stress,ntens,1d0/empa,.true.)
c     call w_chr( idia,'* stran n+1')
c     call w_dim( idia,stran,ntens,1d0/empa,.true.)
         call w_chr( idia,'* ddsdde')
         call w_mdim(idia,ddsdde,  ntens,1d0/gpa)
         call w_chr( idia,'* stran_el_ns(n+1)')
         auxntens(:)=stran_el_ns(1,:)
         call w_dim( idia,auxntens,ntens,1d0,.true.)
      endif

c$$$  4. Update all other state varaiables
      eeq_ns(1)        = eeq_ns(0) + deq
      stran_pl_ns(1,:) = stran_pl_ns(0,:)
      call update_yldp(iyld_law,yldp_ns,nyldp,deq)

c$$   5. Store updated state variables to statev
      auxntens(:)=stran_el_ns(1,:)
      buxntens(:)=stran_pl_ns(1,:)
      auxnyldp(:)=yldp_ns(1,:)
      call restore_statev(statev,nstatv,eeq_ns(1),auxntens,
     $     buxntens,ntens,auxnyldp,nyldp,1,.false.,idia,
     $     .false.,-1,-1,-1,-1d0,stress)
      stran_el_ns(1,:)=auxntens(:)
      stran_pl_ns(1,:)=buxntens(:)
      yldp_ns(1,:)=auxnyldp(:)
      if (idiaw) then
         call w_chr(idia,'** Leaving update_elastic')
      endif
      return
      end subroutine update_elastic
