c======================================================================c
c     Calculates grain growth speed in terms of the rate of changes in
c     grain volume fraction (wdot). This requires mobility M and driving
c     pressure P since the grain boundary migration speed (v) is:
c        v = M P              Eq.1
c     P may comprise of the difference in stored energy and the
c     grain boundary surface tension. The latter contribution is not
c     considered in the current version as the Bailey-Hirsch criterion
c     leads to an absurdly high critical radius of nucleus, requiring
c     much higher stored energy contributed from the stored energy
c     alone. As a result, nuclei with the initial radius of submicron
c     will not spontaneously grow, with its extremely high Laplace.
c     Notably, the assumption of spherical shape for grain may
c     significantly deviate from observed grain or nuclei morphology,
c     which tends to be faceted to minimize the surface energy.
c     Eq.1 leads to a system of linear equations:
c
c     A_ij wdot_j = b_i              Eq.2
c
c     This equation is solved in subroutine <rexdd_wdot>.
c
c     A sub-stepping
c     While the Eq. 2 is given in the 'rate' form, the volume fraction (w)
c     may non-linearly evolve since the right vector b_i is a function of
c     w (or radius). As a result, Eq. 2 is solved, and the results should
c     is linearly extrapolated for a sufficiently small time increment.
c     To ensure the robustness and computational efficience, a substepping
c     algorithm is implemented to prevent an exceedingly large change
c     in radius (or w). Also, in doing so, we prevent any volume fraction
c     being negative at the end of sub-step increment.
c
c     For the time increment that begins t^k and ends t^(k+1) with
c     the time increment dt, we may subdivide dt to several substeps,
c     dt_hat. For each substep dt_hat, the updates in grain volume
c     fraction (wgt) and the grain radius (r) are performed at the
c     end of each sub-step.

c     It should be noted that the grain volume fraction (array wgt)
c     is NOT updated here within this subroutine. It is 'dwgt_grow'
c     which corresponds to the change in volume fraction due to growth
c     accumulates with the time increment t^(k+1)-t^(k). WGT^(k+1) is
c     updated in <rexdd_main>.
c
      subroutine rexdd_grow(
c          intent(in)::
c          * integer
     $     iwritecl,uw_rex_grow,iph,nphmx,ngrmx,ngr,nchildmx,nmodmx,
     $     nclustermx,ncluster,idclusterset,nchildset,
c          * real
     $     temp,dtime,vpx,rexgrowpar,burg,mu_mode,rho_g,dwgt_nucl,wgt,

c         intent(inout)::
c         *real
     $     radius,dwgt_grow,
c         intent(out)::
c         * real
     $     velgrow)
      implicit none
c     Grains that are newly nucleated are not subject to growth mechanisms.
c     Those newly nucleated grains will be considered for growth from
c     the next incremental step onwards.

c     intent(in)
      integer, intent(in):: iwritecl,uw_rex_grow,iph,nphmx,nchildmx,
     $     ngrmx,ngr(0:nphmx),nmodmx,nclustermx,ncluster,
     $     idclusterset(nclustermx,0:nchildmx),nchildset(nclustermx)
      real*8, intent(in):: temp,dtime,vpx,rexgrowpar(5),burg(nmodmx),
     $     mu_mode(nmodmx),rho_g(ngrmx),dwgt_nucl(ngrmx),wgt(ngrmx)
c     intent(inout)
      real*8, intent(inout):: radius(ngrmx),dwgt_grow(ngrmx)
c     intent(out)
      real*8, intent(out):: velgrow(ngrmx)
c     locals
      integer idg(nchildmx)     !: grain pointer array that collects saved grains
      integer nmem_l            !: current number of members at l-th sub-step
c                                  that participate in growth competition
      integer ind,uwx,i,count_l,count_lmx,count_lmn,kk,kkk,icluster,
     $     imem,nmem_at_tk,istat(nchildmx),count_l_stat(nclustermx)
      real*8 count_l_avg,fact,c,dwgt_acc(ngrmx),dt_acc,dt_hat,wgtmp,
     $     wgt_l,mob,gam,pextra,wgt_ref,rvel,wgt_new,deltaQ0b,pi,Qb,
     $     eta,radius_old(ngrmx),Mij(nchildmx,nchildmx),
     $     GAMij(nchildmx,nchildmx),wdot(nchildmx),bound_mn,
     $     drad_hat(nchildmx)
      real*8 drc,drad(0:nchildmx),dc,Br(nchildmx),Bw(nchildmx)
      logical is_bw_bound
c     The maximum allowed increment is set to 1 micro meter.
c     Should it be lower than this?
      parameter(drc=1e-06)
c     --
      pi=4d0*dATAN(1d0)
      uwx=uw_rex_grow           ! file ID
c     uwx=6                     ! standard output

c     several constants
      fact=3./4./pi*vpx
      c=4.*pi/vpx ! not sure if this is ever used in this subroutine.
      dc=fact**(-1/3.)*drc
c     --

c     mobility (mob)
      Qb=rexgrowpar(1)
      deltaQ0b=rexgrowpar(2)
      gam=rexgrowpar(3)
      pextra=rexgrowpar(4)

      call rexdd_mob(
c          intent(in)
c          * real
     $     temp,Qb,deltaQ0b,
c          ** caution **
c           Burgers vector of the first-appearing slip mode
c           in the single crystal file is used
     $     burg(1),
c          intent(out)
c          * real
     $     mob)

c     Determine dislocation line energy (eta)
      eta = 0.5d0 * mu_mode(1) * burg(1)**2 ! xmus in [MPa], burgers vector in [m]
c     Attention: xmus is typically given in MPa (if not need to change/remove the below line)
      eta = eta * 1e6
c      write(*,'(a5,e11.3)')'eta:',eta      ! now, eta is in [Pa/m^2]

c     volume fraction change accumulated from sub-steps
c     Radius pertaining to the beginning of time increment, t^k,
c     is stored to <radius_old> to later estimate the average
c     radial growth rate (rvel). Rvel is not equivalent to
c     velgrow estimated in

      do kkk=ngr(iph-1)+1,ngr(iph)
         radius_old(kkk)=radius(kkk)
         dwgt_acc(kkk)=0.
         dwgt_grow(kkk)=0.
      enddo

      count_l_stat(:)=0 ! necessary since some cluster may skip growth competition

c     ==================================================================
c     loop over all clusters
      do 100 icluster=1,ncluster
c        the number of grains at the beginning of time increment, i.e., t^{(k)}
         nmem_at_tk=nchildset(icluster)+1
c        status of members.
c        0: not participating from the get-go since it is nucleated
c        1: the member is still active and competes
c       -1: the member has vanished during sub-step
         istat(:)=0
         do imem=0,nchildset(icluster)
            kkk=idclusterset(icluster,imem)
            if (dwgt_nucl(kkk).gt.0) then
               istat(imem+1)=0
            else
               istat(imem+1)=1
            endif
         enddo

c        First of all, trial sub-step time increment is using the
c        t_(k+1)-t_(k), i.e., dtime.
         dt_hat=dtime

c        Initialize the accumulated time <dt_acc> and
c        the number of sub steps (count_l)
         dt_acc=0.
         count_l=0

c        ---------------------------------------------------------------
c        Sub-step time integration
         do while(dt_acc.lt.dtime)

c           [find out subjected members]
            idg(:)=-1
            wdot(:)=0.
            drad_hat(:)=0.

            nmem_l=0
            do imem=0,nchildset(icluster)
               if (istat(imem+1).eq.1) then
                  kkk=idclusterset(icluster,imem)
                  wgtmp=wgt(kkk)+dwgt_nucl(kkk)+dwgt_acc(kkk)
                  if (wgtmp.gt.0) then
                     nmem_l=nmem_l+1
                     idg(nmem_l)=kkk
                  else
                     istat(imem+1)=-1
                  endif
               endif
            enddo

            if (nmem_l.le.1) then
               write(*,*)'There cannot be any competition'
               goto 100
            endif

c           ------------------------------------------------------------
c           [calculate wdot]
c           calculate the growth rates of competing members in the cluster
            Mij(:,:)=0.
            GAMij(:,:)=0.
            Mij(:nmem_l,:nmem_l)=mob
            GAMij(:nmem_l,:nmem_l)=gam
            if (nmem_l.eq.1)write(*,*)'unsure if we need to obtain wdot'
c$$$            write(*,*)
c$$$            write(*,*)'------------------------------------'
c$$$            write(*,*)'c:',c,'nmem_l:',nmem_l,'gam:',gam,'eta:',eta
            call rexdd_wdot(
c                intent(in)
     $           ngrmx,nchildmx,nmem_l,idg,c,eta,Mij,GAMij,rho_g,radius,
c                intent(out)
     $           wdot)
c$$$            do i=1,nmem_l
c$$$               kkk=idg(i)
c$$$               write(*,'(a6,i2.2,a2,e13.5,a10,e13.5)')
c$$$     $              'w_dot(',i,'):',wdot(i),
c$$$     $              'radius:',radius(kkk)
c$$$            enddo

c           ------------------------------------------------------------
c           [determine dt_hat that is within the bounds]
c           ** New clearer strategy implementation (may 6, 2026)

c           two bounding values
            br(:)=1e10          ! bounding value from |dr|<|dr_c|
            bw(:)=1e10          ! bounding value from w>0
            bound_mn=1e10

            ind=-1
            is_bw_bound=.false.
            do imem=1,nmem_l
               kkk=idg(imem)
c              wgt^{(l)}: vol. frac. to be accmulated if were to use the
c              given trial dt_hat.
               wgt_l=wgt(kkk)+dwgt_nucl(kkk)+dwgt_acc(kkk)
c              Bw: Bound by weight applies to members with negative growth (wdot<0)
               if (wdot(imem).lt.0) then
                  bw(imem)=-wgt_l/wdot(imem)
                  if (bound_mn.gt.bw(imem)) then
                     bound_mn=bw(imem)
                     ind=imem
                     is_bw_bound=.true.
                  endif
               endif
c              Br: Bound by the increment in grain radius change
               br(imem)=wgt_l**(1./3.)+sign(1.,wdot(imem))*dc
               br(imem)=(br(imem)**3-wgt_l)/wdot(imem)
               if (bound_mn.gt.br(imem)) then
                  bound_mn=br(imem)
                  ind=imem
                  is_bw_bound=.false.
               endif
            enddo
c           Adjust dt_hat, if necessary.
            if (bound_mn.lt.dt_hat) dt_hat=bound_mn

c           Now that the current substep time increment (dt_hat) has
c           been determined using the two bounds, Br and Bw, the grain
c           volume fraction w^l and grain radius r^l pertaining to the
c           end of time increment dt_hat are updated.
c
c           Note here that, if the dt_hat was bounded by Bw (not Br), the
c           corresponding grain should vanish at the end of the current
c           substep (l) and should be excluded from wdot calculation.
c           ------------------------------------------------------------
c           [Updates in accumulated quantities from the substeps]
c              ** dt_acc
c              ** dwgt_acc
c              ** r^{(l)}
c              ** Next guess on dt_hat = dt - dt_acc

c           update dt_acc
            dt_acc=dt_acc+dt_hat

            do imem=1,nmem_l
               kkk=idg(imem)
c              update w
               if (is_bw_bound.and.imem.eq.ind) then
                  dwgt_acc(kkk)=-wgt(kkk)-dwgt_nucl(kkk)
               else
                  dwgt_acc(kkk)=dwgt_acc(kkk)+wdot(imem)*dt_hat
               endif

c              If nmem_l is not greater than 1, there should not be
c              any intergranular competition for growth/shrinkage.
               if (nmem_l.gt.1) then
                  wgtmp=wgt(kkk)+dwgt_nucl(kkk)+dwgt_acc(kkk)
                  drad_hat(imem)=(wgtmp*fact)**(1./3.)-radius(kkk)
                  radius(kkk)=radius(kkk)+drad_hat(imem)
               else
                  wgtmp=0.
                  radius(kkk)=0.
               endif
            enddo

c           substep counts
            count_l=count_l + 1
c           [Write some quantities to the file]
            if (icluster.eq.iwritecl) then
               write(uwx,'(6a10)')'nmem_tk','count_l','nmem_l',
     $              'dt_hat','dt_acc','dtime'
               write(uwx,'(3i10)',advance='no')nmem_at_tk,count_l,nmem_l
               write(uwx,'(3f10.5)')dt_hat,dt_acc,dtime
               write(uwx,'(13a13)')'icl','imem','kkk','status','w_k',
     $              'dwgt_nuc','rho','dwgt_acc','wdot_l',
     $              'rad_k[um]','radius[um]','dR_l[um]','Rdot[mu/s]'
               i=1 ! i is the count for only those active
               do imem=1,nmem_at_tk
                  kkk=idclusterset(icluster,imem-1)
                  write(uwx,'(4i13)',advance='no')icluster,imem,
     $                 kkk,istat(imem)
                  write(uwx,'(4e13.5)',advance='no')
     $                 wgt(kkk),dwgt_nucl(kkk),rho_g(kkk),dwgt_acc(kkk)
                  if (istat(imem).eq.1) then
                     write(uwx,'(e13.5)',advance='no') wdot(i)
                     i=i+1
                  elseif(istat(imem).eq.0) then
                     write(uwx,'(a13)',advance='no')'*nucleated*'
                  elseif(istat(imem).eq.-1) then
                     write(uwx,'(a13)',advance='no')'*vanished*'
                  endif
                  write(uwx,'(2f13.5)',advance='no')
     $                 radius_old(kkk)*1e6,radius(kkk)*1e6
                  if (istat(imem).eq.1) then
                     write(uwx,'(f13.5,e13.5)') drad_hat(i-1)*1e6,
     $                    drad_hat(i-1)*1e6/dt_hat
                  else
                     write(uwx,'(2a13)') '*NA*','*NA*'
                  endif
               enddo

            endif

            if (count_l.gt.1000) then
               write(*,*)'** A bit too many sub steps required.'
               write(*,*)'count_l:',count_l
               write(*,*)'nmem_l:',nmem_l,'nmem_at_tk:',nmem_at_tk
               write(*,*)'dtime:',dtime
               write(*,*)'dt_acc:',dt_acc
               write(*,*)'dt_hat:',dt_hat
               call stop_debug(0,'** Error in <rexdd_grow>')
            endif
c           Guess dt_hat for the next sub step = dt-dt_acc
            dt_hat = dtime-dt_acc
         enddo                  !! end of sub-steps
         count_l_stat(icluster)=count_l
 100  continue  !! end of loop over clusters

      count_l_avg=0.
      count_lmn=100000
      count_lmx=-count_lmn
      do icluster=1,ncluster
         count_l_avg=count_l_avg+count_l_stat(icluster)
         if (count_lmn.gt.count_l_stat(icluster))
     $        count_lmn=count_l_stat(icluster)
         if (count_lmx.lt.count_l_stat(icluster))
     $        count_lmx=count_l_stat(icluster)
      enddo
      count_l_avg=count_l_avg/ncluster
      write(*,'(4a10)')'ncl','min l','max l','avg l'
      write(*,'(3i10,f10.1)')ncluster,count_lmn, count_lmx,count_l_avg

c     Determine grain boundary migration velocity (velgrow)
c     which can be positive (growth) or negative (shrinkage)
c     The overal average  g.b. migration velocity (rvel)
c     should be greater than or equal to zero.
      kk=0
      wgt_ref=0.
      rvel=0. !! common
      do icluster=1,ncluster
         drad(:)=0.
         do imem=0,nchildset(icluster)
            kkk=idclusterset(icluster,imem)
            dwgt_grow(kkk)=dwgt_acc(kkk)
            drad(imem)=radius(kkk)-radius_old(kkk)
            velgrow(kkk)=drad(imem)/dtime
            kk=kk+1
            rvel=rvel+drad(imem)**2*wgt(kkk)
            wgt_ref=wgt_ref+wgt(kkk)
         enddo
      enddo
      rvel=sqrt(rvel/wgt_ref)/dtime
      write(*,'(a20,f10.2,a10)')'GB vel:',rvel*1e6,'[um/s]'

      return
      end subroutine rexdd_grow

c======================================================================c
c     Solving the growth equation, A wdot = b
      subroutine rexdd_wdot(
c          intent(in)
c          * integer
     $     ngrmx,nchildmx,nmem,idg,
c          * real
     $     c,eta,Mij,GAMij,rho_g,radius,
c          intent(out)
c          * real
     $     wdot)
      implicit none
      integer, intent(in):: ngrmx,nchildmx,nmem,idg(nchildmx)
      real*8, intent(in):: c,eta,
     $     Mij(nchildmx,nchildmx),GAMij(nchildmx,nchildmx),
     $     rho_g(ngrmx),radius(ngrmx)
      real*8, intent(out):: wdot(nchildmx)
c     locals
      real*8 Aij(nmem,nmem),bi(nmem),wdotsum,lambda,rho_i,rho_j,sum,
     $     r_i,r_j
      parameter(lambda=1e15)     !! penalty multiplier
      integer i,j,ikkk,jkkk,isingular

c     Create Aij and Bi matrices
      bi(:)=0.
      Aij(:,:)=lambda           ! off diagonal component = lambda
      do i=1,nmem
         Aij(i,i)=Aij(i,i)+1.    ! making the diagonal component 1+lambda
      enddo
c     Bi matrix
      do i=1,nmem
         ikkk=idg(i)
         rho_i=rho_g(ikkk)
c        sum
         sum=0.
         r_i=radius(ikkk)
         do j=1,nmem
            if (j.ne.i) then
               jkkk=idg(j)
               rho_j=rho_g(jkkk)
               r_j=radius(jkkk)
               sum=sum+Mij(i,j)*(eta*(rho_j-rho_i)
     $              +2.*GAMij(i,j)*(1./r_j-1./r_i)
     $              )
            endif
         enddo                  ! loop over j
         bi(i)=c*(r_i**2)*sum
      enddo                     ! loop over i

      call lu_eqsystem(aij,bi,nmem,isingular)
      IF(ISINGULAR.EQ.1) THEN
         call stop_debug(0,'Error: failed to solve the growth equation')
      ENDIF

c     Correcting the wdot by distributing any residuals.
c       * getting the residual (in principle it should be zero)
      wdotsum=0.
      do i=1,nmem
         wdotsum=wdotsum+bi(i)
      enddo
      wdot(:)=0.
      do i=1,nmem
         wdot(i)=bi(i)-wdotsum/nmem
      enddo

      return
      end subroutine rexdd_wdot

c======================================================================c
c     Mobility calculation using Eqs 18 and 19 of
c     Ding and Guo, Acta Mater. 49 (2001)
      subroutine rexdd_mob(
c          intent(in)::
     $     temp,Qb,deltaQ0b,burg,
c          intent(out)::
     $     mob)
      implicit none
      real*8, intent(in)::temp,Qb,deltaQ0b,burg
      real*8, intent(out)::mob
c     locals
      real*8 boltz,rgas,rt_inv,kt_inv
      parameter(boltz=1.380649e-23) ! unit of J/K
      parameter(rgas=8.314462) ! unit of J/(K.mol)
      rt_inv=1./(rgas*temp)
      kt_inv=1./(boltz*temp)
      mob=deltaQ0b*exp(-Qb*rt_inv)*burg*kt_inv
      return
      end subroutine rexdd_mob
