c----------------------------------------------------------------------------
c     This source code contains three subroutines:
c        subroutine data_grain_odf
c        subroutine weighted_grains
c        subroutine trilinear
c     The main purpose of these subroutines is to obtain
c      a population of discrete orientations with uniformly weights, thus isotropic,
c      or weights from experimental COD files (cod array in common block/kcod/).
c     The switch 'iodf' given into the subroutine <data_grain_odf>, passed to into
c     subr. <weighted_grains> determines the behavior
c      (iodf=2 for the isotropic texture, iodf=1 for the case using experimental COD)
c----------------------------------------------------------------------------
      subroutine data_grain_odf(
c          in
     $     ife,iodf,ngr,UR_TX,UR_SHP,ngrmx,nphmx,iphel,igb,igt,wph,
     $     ishape,eulerph,axisph,
c          out
     $     wgt,ag,fijph,fijgr)
c     in case ife is 0 -> Call from standalone
c     in case ife is 1 -> Call from UMAT
      implicit none
      integer,intent(in):: ife,iodf,ngr,UR_TX,UR_SHP,
     $     ngrmx,nphmx,iphel,igb,
     $     igt,ishape(0:nphmx)
      real*8,intent(in):: wph(nphmx),eulerph(3,0:nphmx),
     $     axisph(0:3,1:3,0:nphmx)
      real*8,intent(out)::wgt(ngrmx),ag(3,3,1:ngrmx),fijph(3,3,0:nphmx),
     $     fijgr(3,3,1:ngrmx)

c     locals
      real*8 px(ngrmx,4),eul1,eul2,eul3,aa(3,3),totwgt,da,db,dc,
     $     fijx(3,3),fnew(3,3),ax(3),px_ph(ngr,4)
      integer kkk,i,j,k,m,idum
c      character(len=1) nomen
      character(len=256) prosa
c      parameter(nomen='B')

c     below is to download the obtained "sampled" discrete orientation
c     to files
      logical idownload
      integer iudown
c     flag to download the generated discretized orientations.
      parameter(idownload=.false.)
      parameter(iudown=9832)

      if (iodf.eq.0) then
         write(*,*)'Error: iodf==0 was given in subr. data_grain_odf'
         call xit()
      endif

c     to prevent the random_number uses the same starting position
      call weighted_grains(iodf,UR_TX,iphel,ngr,px_ph)
      px(igb:igt,:)=px_ph(1:ngr,:)
      totwgt=0.
      do kkk=igb,igt
         eul1=px(kkk,1)
         eul2=px(kkk,2)
         eul3=px(kkk,3)
         wgt(kkk)=px(kkk,4)
         totwgt=totwgt+wgt(kkk)
         call euler(2,eul1,eul2,eul3,aa)

         do j=1,3
            do k=1,3
               ag(j,k,kkk)=aa(k,j)
            enddo
         enddo
      enddo

      do kkk=igb,igt
         wgt(kkk)=wgt(kkk)/totwgt
      enddo

      if (ife.eq.0.and.idownload) then
         open(action='write',unit=iudown,
     $        file='initial_texture_from_COD.tex',status='unknown')
         write(iudown,*)
         write(iudown,*)
         write(iudown,*)
         write(iudown,*) 'B',igt-igb+1
         do kkk=igb,igt
            write(iudown,'(4e12.4)')px(kkk,1),px(kkk,2),px(kkk,3),
     $           wgt(kkk)
         enddo
         close(iudown)
      endif

      do kkk=igb,igt
         wgt(kkk)=wgt(kkk)*wph(iphel)
      enddo

      IF(ISHAPE(iphel).EQ.3) THEN
CFEB -----------------------------------------------------------------
C     WRITES MORPHOLOGIC TEXTURE FILE INTO 'RUN_LOG.OUT' FILE
         DO IDUM=1,10
            READ(UNIT=UR_SHP,END=99,FMT='(A)') PROSA
         ENDDO
 99      REWIND UR_SHP
CFEE -----------------------------------------------------------------
         READ(UR_SHP,'(A)') PROSA
         READ(UR_SHP,'(A)') PROSA
         READ(UR_SHP,'(A)') PROSA
         READ(UR_SHP,'(A)') PROSA
C *** 'EULERPH' AND 'AXISPH(0') ARE READ in VPSC_INPUT from VPSC.IN.
C ***  THEY MAY BE RE RE-READ FROM THE MORPHOLOGIC TEXTURE FILE IF ishape=3.
C         READ(UR_SHP,*) (AXISPH(0,I,IPHel),I=1,3)
C         READ(UR_SHP,*) (EULERPH(I,IPHel) ,I=1,3)
      ENDIF

      da=EULERPH(1,iphel)
      db=EULERPH(2,iphel)
      dc=EULERPH(3,iphel)

      call euler(2,da,db,dc,aa)

      do i=1,3
      do j=1,3
         fijx(i,j)=(i/j)*(j/i)*AXISPH(0,I,iphel)
      enddo
      enddo

      do j=1,3
      do i=1,3
         fnew(i,j)=0.
      do m=1,3
         fnew(i,j)=fnew(i,j)+aa(m,i)*fijx(m,j)
      enddo
      enddo
      enddo

      do i=1,3
      do j=1,3
         fijph(i,j,iphel)=fnew(i,j)
      enddo
      enddo

      IF(ISHAPE(iphel).EQ.1 .OR. ISHAPE(iphel).EQ.2) THEN
         DO KKK=igb,igt
         DO I=1,3
         DO J=1,3
            FIJGR(I,J,KKK)=FNEW(I,J)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

C *** ISHAPE=3 -> READ EULER ANGLES & ASPECT RATIOS OF EACH ELLIPSOID.
C *** DA,DB,DC -> EULER ANGLES OF ELLIPSOID AXES WRT SAMPLE.

      IF(ISHAPE(iphel).EQ.3) THEN
        DO KKK=igb,igt
           if(ishape(iphel).eq.3) then
              READ(UR_SHP,*) da,db,dc,ax(1),ax(2),ax(3)
              do i=1,3
              do j=1,3
                 fijgr(i,j,KKK)=(i/j)*(j/i)*ax(i)
              enddo
              enddo
           endif
           call euler(2,da,db,dc,aa)
C     AA TRANSFORMS FROM SAMPLE TO GRAIN AXES. AAt FROM GRAIN TO SAMPLE.
          do j=1,3
          do i=1,3
             fnew(i,j)=0.
             do m=1,3
                fnew(i,j)=fnew(i,j)+aa(m,i)*fijgr(m,j,kkk)
             enddo
          enddo
          enddo
          do i=1,3
          do j=1,3
             fijgr(i,j,KKK)=fnew(i,j)
          enddo
          enddo
        ENDDO      ! END OF LOOP OVER GRAINS
      ENDIF      ! END OF CONDITION FOR ISHAPE=3

      return
      end subroutine data_grain_odf



c     ---------------------------------------------------------------
      subroutine weighted_grains(iodf,UR_TX,iphel,ngrainran,px)
      implicit none

      integer, intent(in)::iodf,UR_TX,ngrainran,iphel
c      character(len=256),intent(in)::fin
      real*8, intent(out)::px(ngrainran,4)

      integer i,j,k,nphi,nphi1,nphi2,ix,iy,iz
      integer ngrmx,nr
      parameter(ngrmx=30000)


      real*8 phi,phi1,phi2,phi1_pre,phi1_next,val,pi,
     $     r(3,2),c(0:1,0:1,0:1),aux3(3)
      real omer(ngrmx),phir(ngrmx),ther(ngrmx)
      integer ir(3,2)
      real random1,random2,random3
      logical verbose
      parameter(verbose=.false.)
      character(len=256) prosa

c     common blocks
      integer nphimx,nphi1mx,nphi2mx
      parameter(nphi1mx=72+1,nphimx=18+1,nphi2mx=18+1)
      real*8 cod(nphi1mx,nphimx,nphi2mx)
      real*8 phi1mx,phimx,phi2mx,res
      logical is_cod_avail(5) !! polycrystal consisting of 5 phase
      common/kcod/cod,phi1mx,phimx,phi2mx,res,is_cod_avail


      if (ngrainran.gt.ngrmx)then
         call stop_debug(0,
     $        'ngrainran given in subr. weighted_grains exceeds ngrmx')
      endif


      PI=4d0*dATAN(1d0)
      val=0d0

ccc   find the resolution
      !! heading lines
      call mutexlock(4)
c      if (iodf.eq.2) write(*,*)'** random isotropic texture'

      if (.not.is_cod_avail(iphel).and.iodf.eq.1) then
         do i=1,4
            read(UR_TX,'(a)') prosa
         enddo
c     res=phi1_next - phi1_pre
         do i=1,2
            read(UR_TX,*)phi1_pre,phi,phi2,val
            read(UR_TX,*)phi1_next,phi,phi2,val
            res=phi1_next-phi1_pre
         enddo

c     Find the maximum phi1, phi, phi2 value
         k=1
         do while(.true.)
            read(UR_TX,*,end=20)phi1mx,phimx,phi2mx,val
            k=k+1
            if (k.gt.100000)then
               write(*,*)'** line exceeded 100000'
               write(*,*)'** Too many lines... something went wrong?'
               call xit()
            endif
         enddo
 20      rewind UR_TX


         if (verbose) then
            write(*,'(a9,f6.2)')'reso:',res
            write(*,'(a9,f6.2)')'phi1mx:',phi1mx
            write(*,'(a9,f6.2)')'phimx :',phimx
            write(*,'(a9,f6.2)')'phi2mx:',phi2mx
         endif

         nphi1 = int((phi1mx+res) / res) + 1
         nphi  = int(phimx / res) + 1
         nphi2 = int((phi2mx+res) / res) + 1

         if (verbose) then
            write(*,'(a9,i6)')'nphi1:',nphi1
            write(*,'(a9,i6)')'nphi :',nphi
            write(*,'(a9,i6)')'nphi2:',nphi2
         endif
         phi1mx=phi1mx+res
         phi2mx=phi2mx+res

         do i=1,4
            read(UR_TX,*)
         enddo
         do i=1,nphi2-1
            do j=1,nphi
               do k=1,nphi1-1
                  read(UR_TX,*)phi1,phi,phi2,val
                  cod(k,j,i)=val
               enddo
            enddo
         enddo
c        cod values at the edge of Euler spaces
         do j=1,nphi
            do k=1,nphi2
               cod(nphi1,j,k)=cod(1,j,k)
            enddo
         enddo
         do i=1,nphi1
            do j=1,nphi
               cod(i,j,nphi2)=cod(i,j,1)
            enddo
         enddo
         is_cod_avail(iphel)=.true.
c      else
c         write(*,*)'cod(1,1,1):',cod(1,1,1)
      endif

      if (.not.is_cod_avail(iphel).and.iodf.eq.2) then
         write(*,*)'check for iodf.eq.2'
         write(*,*)'ngrainran:',ngrainran
         phi1mx=360.
         phi2mx=360.
         write(*,*)'phi1mx:',phi1mx
         write(*,*)'phi2mx:',phi2mx
         is_cod_avail(iphel)=.true.
      endif


      call mutexunlock(4)



      DO NR=1,ngrainran
        call random_number(random1)
        call random_number(random2)
        call random_number(random3)
        PHIR(NR)=phi1mx*RANDOM1
        THER(NR)=ACOS(RANDOM2)*180./PI
        OMER(NR)= phi2mx*RANDOM3 ! cubic
c        OMER(NR)= 90.*RANDOM3 ! cubic
c$$$        IF(ICRYSYM.EQ.1) OMER(NR)= 90.*RANDOM3 ! cubic
c$$$        IF(ICRYSYM.EQ.2) OMER(NR)= 60.*RANDOM3 ! hexagonal
c$$$        IF(ICRYSYM.EQ.5) OMER(NR)=180.*RANDOM3 ! orthogonal

        if (iodf.eq.1) then
           if ((abs(phir(nr)-phi1mx).le.1e-10)
     $          .or.(abs(phir(nr)).le.1e-10))
     $          then
              write(*,*)'phir(nr) is on the edge'
              call xit()
           endif
           if (abs(ther(nr)-phimx).le.1e-10
     $          .or.(abs(ther(nr)).le.1e-10)) then
              write(*,*)'ther(nr) is on the edge'
              call xit()
           endif
           if (abs(omer(nr)-phi2mx).le.1e-10
     $          .or.(abs(omer(nr)).le.1e-10)) then
              write(*,*)'omer(nr) is on the edge'
              call xit()
           endif
        endif

      ENDDO

      do nr=1,ngrainran
         if (verbose) then
            write(*,'(a9,f6.2)')'phir(nr):',phir(nr)
            write(*,'(a9,f6.2)')'ther(nr):',ther(nr)
            write(*,'(a9,f6.2)')'omer(nr):',omer(nr)
         endif

         if (iodf.eq.1) then
            r(1,1)=phir(nr)-mod(phir(nr),res)
            r(1,2)=phir(nr)-mod(phir(nr),res)+res
            r(2,1)=ther(nr)-mod(ther(nr),res)
            r(2,2)=ther(nr)-mod(ther(nr),res)+res
            r(3,1)=omer(nr)-mod(omer(nr),res)
            r(3,2)=omer(nr)-mod(omer(nr),res)+res
            do i=1,3
               if (verbose) write(*,'(a9,2f6.1)')
     $              'r(i,1:2):',(r(i,j),j=1,2)
               do j=1,2
                  ir(i,j)=int((r(i,j)+res)/res)
               enddo
            enddo

            do i=1,2
               do j=1,2
                  do k=1,2
                     ix=ir(1,i)
                     iy=ir(2,j)
                     iz=ir(3,k)
                     phi1=(ix-1)*res
                     phi =(iy-1)*res
                     phi2=(iz-1)*res
                     if (verbose) then
                        write(*,*)'(p1,p,p2,cod)',
     $                       phi1,phi,phi2,cod(ix,iy,iz)
                     endif
                     if (iy.gt.19) then
                        write(*,*)'iy cannot be greater than 19'
                        write(*,*)'ix,iy,iz:',ix,iy,iz
                        write(*,*)'nr:',nr,'ngrainran:',ngrainran
                        write(*,'(a9,f6.2)')'phir(nr):',phir(nr)
                        write(*,'(a9,f6.2)')'ther(nr):',ther(nr)
                        write(*,'(a9,f6.2)')'omer(nr):',omer(nr)
                        call xit()
                     endif
                     c(i-1,j-1,k-1)=cod(ix,iy,iz)
                  enddo
               enddo
            enddo




            aux3(1)=phir(nr)
            aux3(2)=ther(nr)
            aux3(3)=omer(nr)
            call trilinear(c,aux3,r(1,1),r(1,2),r(2,1),r(2,2),r(3,1),
     $           r(3,2),val)
            if (verbose) write(*,'(a9,f9.5)')'val:',val
         elseif (iodf.eq.2) then
            val=1d0             ! random isotropic texture
         endif

         px(nr,1)=phir(nr)
c         px(nr,1)=phir(nr) + 45
c         px(nr,1)=phir(nr) + 90
         px(nr,2)=ther(nr)
         px(nr,3)=omer(nr)
         px(nr,4)=val
      enddo

      return
      end subroutine weighted_grains

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine trilinear(cijk,xyz,x0,x1,y0,y1,z0,z1,c)
      implicit none
      real*8,intent(in)::cijk(0:1,0:1,0:1),x0,x1,y0,y1,z0,z1,
     $     xyz(3)
      real*8,intent(out)::c
c     locals
      real*8 x,y,z,xd,yd,zd,c00,c01,c10,c11,c0,c1
      x=xyz(1)
      y=xyz(2)
      z=xyz(3)

c     distance
      xd=(x-x0)/(x1-x0)
      yd=(y-y0)/(y1-y0)
      zd=(z-z0)/(z1-z0)

c     Apply linear interpolation on the top and bottom faces
c     to obtain two separate interpolants on the respective faces.
      c00=cijk(0,0,0)*(1-xd)+cijk(1,0,0)*xd
      c01=cijk(0,0,1)*(1-xd)+cijk(1,0,1)*xd
      c10=cijk(0,1,0)*(1-xd)+cijk(1,1,0)*xd
      c11=cijk(0,1,1)*(1-xd)+cijk(1,1,1)*xd

c     linear interpolation using the two interpolants
c     on the top and bottom faces
      c0=c00*(1-yd)+c10*yd
      c1=c01*(1-yd)+c11*yd

      c=c0*(1-zd)+c1*zd
      return
      end subroutine trilinear
