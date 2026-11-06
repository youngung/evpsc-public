c     SUBROUTINE select_points_from_ebsd
c     Arguments
c     ---------
c     <coords>  : Coordinates of FE int. point from Abaqus
c     <ncrds>   : Number of coordiantes (2 for shell, 3 for 3D solid)
c                 This determines the dimension of array <coords>
c     <ebsd_rad>: radius, within which the point of EBSD data is
c                 subjected as the pool for discrete orientations
c     <ngrmx>   : number of data point (this is fixed as 'NGRMX' for IODF=3 option)
c     <ebsd_select_data> : array of <phi1,Phi,phi2,weight>
      SUBROUTINE select_points_from_ebsd(
c          in
     $     coords,ncrds,ebsd_rad,ngrmx,
c          out
     $     ebsd_select_data)
      implicit none

      integer, intent(in):: ngrmx, ncrds
      real*8, intent(in):: ebsd_rad,coords(ncrds)
      real*8, intent(out):: ebsd_select_data(4,ngrmx)

      integer nlines_ebsd,nrowmx
c      parameter(nrowmx=200000)
      parameter(nrowmx=200)
      real*4 ebsd_raw_data(5,nrowmx)
      common/kebsd/ebsd_raw_data,nlines_ebsd

c     local
      integer i,nselect_ebsd,nrd,irecipe
      real*8 pi,dist, dist_x, dist_y,x,y,z,span,ebsd_within_r(6,nrowmx),
     $     rand_num,bm_begx,bm_endx,dist_mn,xyz_mx(3),xyz_mn(3)
c     (1)   Mg10Gd of Jaeseong for C-ring compression;
c     (2-6) POSCO-welds, correspoding respectively to FBW00, 05, 10, 15, and 20.
      parameter(irecipe=2)

      pi = ACOS(-1.)
      nselect_ebsd=0
      if (irecipe.eq.1) then
c        Mg10Gd's C-ring compression simulation of Jaeseong
         x = coords(1)
         y = coords(2)
         do i=1,nlines_ebsd
            dist_x = (ebsd_raw_data(4,i)-x)**2
            dist_y = (ebsd_raw_data(5,i)-y)**2
            dist = sqrt(dist_x + dist_y)

            IF (dist.LE.ebsd_rad) THEN
               nselect_ebsd = nselect_ebsd +1
               ebsd_within_r(1:3,nselect_ebsd) = ebsd_raw_data(1:3,i)*
     $              (180./pi)
               ebsd_within_r(4:5,nselect_ebsd) = ebsd_raw_data(4:5,i)
               ebsd_within_r(6,nselect_ebsd) = dist
            ENDIF
         enddo

c     POSCO-welds + a reference case
      elseif (irecipe.ge.2.and.irecipe.le.7) then
         if (irecipe.eq.2) then
            bm_begx=1490e-3
            bm_endx=3133e-3     !! FBW00 condition
         elseif (irecipe.eq.3) then
            bm_begx=1956e-3
            bm_endx=2689e-3     !! FBW05 condition
         elseif (irecipe.eq.4) then
            bm_begx=2472e-3
            bm_endx=3106e-3     !! FBW10 condition
         elseif (irecipe.eq.5) then
            bm_begx=4020e-3
            bm_endx=5015e-3     !! FBW15 condition
         elseif (irecipe.eq.6) then
            bm_begx=5531e-3
            bm_endx=6606e-3     !! FBW20 condition
         elseif (irecipe.eq.7) then
c           a special case to use base material EBSD scan of FBW20.
            bm_begx=5531e-3
            bm_endx=6606e-3
         endif
c         write(*,*)"EBSD recipe for POSCO's case with welds"
c         write(*,*)'ebsd_rad:',ebsd_rad
c         write(*,*)'nlines_ebsd:',nlines_ebsd
         x = coords(1)
         y = coords(2)
         z = coords(3)
c         write(*,'(a10,3f9.2)')'x,y,z:',x,y,z

         if (x.gt.bm_endx-ebsd_rad.or.irecipe.eq.7) then
c     write(*,*)'-- Base material zone'
c     The trick is to reassign a random x coordinate
c     within the base material EBSD zone.
c     The beginning x coordinate of base material region
c     is 1490~bm_endx. Since we should also consider the
c     probing radius (ebsd_rad), the possible range of
c     reassigned x is between (1490+ebsd_rad) ~ (bm_endx-ebsd_rad)
            call random_number(rand_num)
            span = bm_endx - bm_begx - ebsd_rad - ebsd_rad
            if (span.le.0) then
               write(*,*)'** Error: span cannot be negative or zero.'
               write(*,*)'span:',span
               call xit()
            elseif(bm_endx-bm_begx .lt. ebsd_rad * 4) then
               write(*,*)"** Error: ebsd_rad is too large for"
               write(*,*)"   given base-metal' EBSD span."
               write(*,*)"   <ebsd_rad> should be less than",
     $              (bm_endx-bm_begx)/4.
               write(*,*)"   But <ebsd_rad> of",ebsd_rad,"was used."
               call xit()
            endif
            x=rand_num*span + bm_begx + ebsd_rad
c$$$            write(*,*)'** span:',span
c$$$            write(*,*)'** perturbed x coordinate:',x
c$$$            write(*,*)'** rand_num:',rand_num
         endif

         if (x.gt.bm_endx-ebsd_rad) then
            write(*,*)'** Error: x cannot exceed ',bm_endx-ebsd_rad
            call xit()
         endif

         dist_mn=1e10
         xyz_mx(:)=-1e10
         xyz_mn(:)=+1e10

         do i=1,nlines_ebsd

            !! mx in x,z
            if (ebsd_raw_data(4,i).gt.xyz_mx(1)) then
               xyz_mx(1)=ebsd_raw_data(4,i)
            endif
            if (ebsd_raw_data(5,i).gt.xyz_mx(3)) then
               xyz_mx(3)=ebsd_raw_data(5,i)
            endif
            !! mn in x,z
            if (ebsd_raw_data(4,i).lt.xyz_mn(1)) then
               xyz_mn(1)=ebsd_raw_data(4,i)
            endif
            if (ebsd_raw_data(5,i).lt.xyz_mn(3)) then
               xyz_mn(3)=ebsd_raw_data(5,i)
            endif
            !! --


            dist_x = (ebsd_raw_data(4,i)-x)**2
            dist_y = (ebsd_raw_data(5,i)-z)**2
            dist = sqrt(dist_x + dist_y)
            if (dist.lt.dist_mn) dist_mn=dist
            IF (dist.LE.ebsd_rad) THEN
               nselect_ebsd = nselect_ebsd +1
c              Euler angles in degree
               ebsd_within_r(1:3,nselect_ebsd) = ebsd_raw_data(1:3,i)
c              EBSD coordinates
               ebsd_within_r(4:5,nselect_ebsd) = ebsd_raw_data(4:5,i)
c              Distance (for reference)
               ebsd_within_r(6,nselect_ebsd) = dist
            ENDIF
         enddo
      endif

c      write(*,*)'nselect_ebsd:',nselect_ebsd

      if (nselect_ebsd.eq.0) then
         write(*,*)'-----------------------------------'
         write(*,*)'** Error : unexpected nselect_ebsd'
         write(*,*)'**    nselect_ebsd is zero'
         write(*,*)'nlines_ebsd:',nlines_ebsd
         write(*,*)'nrowmx:',nrowmx
         write(*,*)'irecipe:',irecipe
         write(*,*)'bm_begx:',bm_begx
         write(*,*)'bm_endx:',bm_endx
         write(*,*)'x range in EBSD:',xyz_mn(1),'-',xyz_mx(1)
c         write(*,*)'y range',xyz_mn(2),'-',xyz_mx(2)
         write(*,*)'z range in EBSD:',xyz_mn(3),'-',xyz_mx(3)
         write(*,*)'-----------------------------------'
         write(*,'(3(a2,f10.4))')'x:',x,'z:',z, 'y:',y
         write(*,'(2(a15,f10.3))')'dist_mn:',dist_mn,
     $        'ebsd_rad:',ebsd_rad
         write(*,*)'-----------------------------------'
         call xit()
      endif

      do i=1,ngrmx
         call random_number(rand_num)
         nrd = FLOOR((nselect_ebsd-1)*rand_num)+1
         ebsd_select_data(1:3,i) = ebsd_within_r(1:3,nrd)
         ebsd_select_data(4,i) = 1.0
      enddo
      END SUBROUTINE select_points_from_ebsd
