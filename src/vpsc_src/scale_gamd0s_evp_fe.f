C *************************************************************************
C     SUBROUTINE SCALE_GAMD0S   ---->   VERSION 20/MAR/2020
C
C     THIS SUBROUTINE IS MEANT TO BE USED IN COMBINATION WITH MTS OR DD
C     MODELS WHERE RATE EFFECTS ARE ACCOUNTED FOR IN THE FUNCTIONAL FORM
c     OF THE CRSS.
C     GAMD0G (SAME FOR EVERY SYSTEM IN THE GRAIN) WAS REPLACED BY GAMD0S
C     (POTENTIALLY DIFFERENT FOR EVERY SYSTEM IN THE GRAIN)
C
C     BY REDEFINING THE FACTOR 'GAMD0' TO BE OF THE ORDER OF THE STRAIN
C     RATE THE RATIO (RSS/CRSS) WILL BE OF ORDER ONE AND THERE ARE NO RATE
C     SENSITIVITY EFFECTS ASSOCIATED WITH THE N'th POWER.
C     IN SUCH CASE THE N'th POWER IS ONLY A WAY TO AVOID AMBIGUITIES AND
C     GIVES UNIQUE SHEAR RATES IN EACH SYSTEM.
c
c     Extension to VPSCX and incr. EVPSC (and for its FE interface)
C *************************************************************************

      SUBROUTINE SCALE_GAMD0S_evp_fe(nsysmx,nphmx,ngrmx,nph,ngr,nsyst,
     $     DBAR,tg,gamd0s)
      implicit none
      integer,intent(in)::nsysmx,nphmx,ngrmx,nph,ngr(0:nphmx),
     $     nsyst(nphmx)
      real*8,intent(in)::DBAR(6),tg(6,ngrmx)
      real*8,intent(out)::gamd0s(nsysmx,ngrmx)
c     local
      integer iph,kkk,is,kgx,iopt
      real*8 refrate,vnorm
      parameter(iopt=0) !! the usual case
c      parameter(iopt=1) !! ad-hoc condition to mimic stress-relaxation for uranium

      REFRATE=VNORM(DBAR(:5),5)
c     to allow the stress relaxation for the case of iratesens=0 and DBAR=0
c      if (refrate.lt.1e-10) then
      if (refrate.eq.0.) then
         if (iopt.eq.0) then ! Carlos' suggestion
            do iph=1,nph
            do kkk=ngr(iph-1)+1,ngr(iph)
               gamd0s(:,kkk)=vnorm(tg(1:5,kkk),5)
            enddo
            enddo
         elseif(iopt.eq.1) then ! to visualize the stress relaxation.
            refrate=1d0
c           refrate=0.1
c           refrate=0.01
            do iph=1,nph
            do kkk=ngr(iph-1)+1,ngr(iph)
               gamd0s(:,kkk)=refrate
            enddo
            enddo
         endif
      else
         KGX=1
         DO IPH=1,nph
         DO KKK=NGR(IPH-1)+1,NGR(IPH)
            DO IS=1,NSYST(IPH)
               GAMD0S(IS,KKK)=REFRATE
            ENDDO
            KGX=KGX+1
         ENDDO
         ENDDO
      endif

      RETURN
      END subroutine scale_gamd0s_evp_fe
