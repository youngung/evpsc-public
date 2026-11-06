c======================================================================c
c     Assign the properties of kk1 to kk2
      subroutine rexdd_inherits(
c          intent(in)::
     $     kk2,kk1,ngrmx,nsysmx,nclustermx,nchildmx,nchildset

c          intent(inout)::
c          ** REXDD related quantities
     $     ,idclusterset,idclusterset_reverse
     $     ,pnucl,dwgt_nucl,dwgt_grow,velgrow,radius
c          ** grain data & properties
     $     ,ag,wgt,fijgr,axisgr,crss,sg,dsg,sg_old,sgtry,dsgtry
     $     ,dg0,dg_trans,gtotgr,mgtg
c          ** dislocation densities
     $     ,rho_s_for,rho_s_deb,rho_g
     $     )
      implicit none
c     intent(in)::
      integer, intent(in):: kk2,kk1,ngrmx,nsysmx,
     $     nclustermx,nchildmx,nchildset(nclustermx)

c     intent(inout)::
      integer, intent(inout)::idclusterset(nclustermx,0:nchildmx),
     $     idclusterset_reverse(ngrmx)
      real*8, intent(inout)::
c          ** REXD
     $     pnucl(ngrmx),dwgt_nucl(ngrmx),dwgt_grow(ngrmx)
     $     ,velgrow(ngrmx),radius(ngrmx)
c          ** grain data & properties
     $     ,ag(3,3,ngrmx),wgt(ngrmx),fijgr(3,3,ngrmx)
     $     ,axisgr(0:3,3,ngrmx),crss(nsysmx,ngrmx)
     $     ,sg(6,ngrmx),dsg(6,ngrmx),sg_old(6,ngrmx),sgtry(6,ngrmx)
     $     ,dsgtry(6,ngrmx),dg0(6,ngrmx)
     $     ,dg_trans(6,ngrmx),gtotgr(ngrmx),mgtg(6,6,ngrmx)
c          ** dislocation densities
     $     ,rho_s_for(nsysmx,ngrmx),rho_s_deb(nsysmx,ngrmx),rho_g(ngrmx)

c     local
      integer imem,icluster

c     ** rexdd
      idclusterset_reverse(kk2)=idclusterset_reverse(kk1)
      idclusterset_reverse(kk1)=-1 !! the room of kk1 is vacant now.
      icluster=idclusterset_reverse(kk2)

      do imem=0,nchildset(icluster)
         if (idclusterset(icluster,imem).eq.kk1) then
            idclusterset(icluster,imem)=kk2
         endif
      enddo

      pnucl(kk2)=pnucl(kk1)
      dwgt_nucl(kk2)=dwgt_nucl(kk1)
      dwgt_grow(kk2)=dwgt_grow(kk1)
      velgrow(kk2)=velgrow(kk1)
      radius(kk2)=radius(kk1)

c     ** grain data & properties
      ag(:,:,kk2)=ag(:,:,kk1)
      wgt(kk2)=wgt(kk1)
      wgt(kk1)=0.


      axisgr(:,:,kk2)=axisgr(:,:,kk1)
      fijgr(:,:,kk2)=fijgr(:,:,kk1)

c     --
      crss(:,kk2)=crss(:,kk1)
      sg(:,kk2)=sg(:,kk1)
      dsg(:,kk2)=dsg(:,kk1)
      sg_old(:,kk2)=sg_old(:,kk1)
      sgtry(:,kk2)=sgtry(:,kk1)
      dsgtry(:,kk2)=dsgtry(:,kk1)

      dg0(:,kk2)=dg0(:,kk1) ! used in updat_macro_property
      dg_trans(:,kk2)=dg_trans(:,kk1) ! used in update_macro_property
      gtotgr(kk2)=gtotgr(kk1)
      mgtg(:,:,kk2)=mgtg(:,:,kk1) ! used in update_macro_property

c     ** dislocation densities
      rho_s_for(:,kk2)=rho_s_for(:,kk1)
      rho_s_deb(:,kk2)=rho_s_deb(:,kk1)
      rho_g(kk2)=rho_g(kk1)

      end subroutine rexdd_inherits
