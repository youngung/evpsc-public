c----------------------------------------------------------------------c
c     Youngung Jeong @ Changwon National University
c     yjeong@changwon.ac.kr
c----------------------------------------------------------------------c
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT,
     $     DRPLDE,DRPLDT,STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,
     $     DPRED,CMNAME,NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,
     $     PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
      implicit none
c      implicit none
      CHARACTER*80 CMNAME
      dimension STRESS(NTENS),STATEV(NSTATV),DDSDDE(NTENS,NTENS),
     $     DDSDDT(NTENS),DRPLDE(NTENS),STRAN(NTENS),DSTRAN(NTENS),
     $     TIME(2),PREDEF(1),DPRED(1),PROPS(NPROPS),COORDS(3),
     $     DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),JSTEP(4)
c     locals
      integer kdelta(3,3),i,j,k,l,i0,j0,i1,j1
      real lamb, g, poi,dsig(3,3),deps(3,3),bulk,xmod,stiff(3,3,3,3),
     $     aux6(6),aux33(3,3),aux66(6,6)
      integer ijv(6,2),n,m
      data ((ijv(n,m),m=1,2),n=1,6)/1,1,2,2,3,3,1,2,1,3,2,3/


c$$$      write(*,'(a2,2(a1,x))')'vo','i','j'
c$$$      do n=1,6
c$$$         write(*,'(i1,a1,i1,a1,i1)')n,':',ijv(n,1),',',ijv(n,2)
c$$$      enddo
c      call xit()

      write(*,'(2a10)')'kinc','jstep(1)'
      write(*,'(2i10)')kinc,jstep(1)
      write(*,'(6e10.2)')(dstran(i),i=1,ntens)


c     pure elastic model.

      xmod=200e3                 ! 200 GPa
      poi=0.3

c     kronecker delta
      kdelta(:,:)=0
      do i=1,3
         kdelta(i,i)=1
      enddo

      do i=1,3
         deps(i,i)=dstran(i)
      enddo
      deps(1,2)=dstran(4)
      deps(2,1)=dstran(4)
      deps(1,3)=dstran(5)
      deps(3,1)=dstran(5)
      deps(2,3)=dstran(6)
      deps(3,2)=dstran(6)

c     mod and poi to Lame first and shear modulus
      lamb=xmod*poi/(1+poi)/(1-2*poi)
      g=xmod/2./(1+poi)
      bulk=xmod/3./(1.-2.*poi)

      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  stiff(i,j,k,l)=bulk*kdelta(i,j)*kdelta(k,l)+
     $              g*(kdelta(i,k)*kdelta(j,l)+kdelta(j,l)*kdelta(j,k)
     $                 -2./3.*kdelta(i,j)*kdelta(k,l))
               enddo
            enddo
         enddo
      enddo


c     dsig_ij = lamb d_ij de_kk + 2mu de_ij
      do i=1,3
         do j=1,3
            dsig(i,j)=lamb * kdelta(i,j)*
     $           (deps(1,1)+deps(2,2)+deps(3,3))
     $           + 2 * g * deps(i,j)
         enddo
      enddo

c      write(*,*)'stress and dsig'
c      write(*,'(6e10.2,x)',advance='no')
c     $     (stress(i),i=1,ntens)
      do i=1,6
         i0=ijv(i,1)
         j0=ijv(i,2)
         write(*,'(e10.2)',advance='no')dsig(i0,j0)
      enddo
      write(*,*)
      


c     update stress
      do i=1,3
         stress(i)=dsig(i,i)+stress(i)
      enddo
      stress(4)=dsig(1,2)+stress(4)
      stress(5)=dsig(1,3)+stress(5)
      stress(6)=dsig(2,3)+stress(6)



c     obtain jacobian (~= elastic modulus)
      call voigt_abq_strain(aux6,aux33,aux66,stiff,4)
      ddsdde(:,:)=aux66(:,:)
c$$$      do i=1,6
c$$$         do j=1,6
c$$$c            write(*,*)'i,j:',i,j
c$$$            i0=ijv(i,1)
c$$$            j0=ijv(i,2)
c$$$            i1=ijv(j,1)
c$$$            j1=ijv(j,2)
c$$$c$$$            write(*,*)'i0,j0,i1,j1:',i0,j0,i1,j1
c$$$c$$$            write(*,*)'ijv(1,i):',ijv(i,1)
c$$$c$$$            write(*,*)'ijv(2,i):',ijv(i,2)
c$$$c$$$            write(*,*)'ijv(1,j):',ijv(j,1)
c$$$c$$$            write(*,*)'ijv(2,j):',ijv(j,2)
c$$$            ddsdde(i,j)=stiff(i0,j0,i1,j1)
c$$$         enddo
c$$$      enddo
      
      return
      end subroutine




      SUBROUTINE VOIGT_ABQ_strain(T1,T2,C2,C4,IOPT)
      implicit none
      real T1(6),T2(3,3),C2(6,6),C4(3,3,3,3),fact(6,6)
      integer IJV(6,2),n,m,i,j,i1,i2,iopt,j1,j2
      DATA ((IJV(N,M),M=1,2),N=1,6)/1,1,2,2,3,3,1,2,1,3,2,3/
      DATA ((FACT(n,m),n=1,6),m=1,6)/1d0,1d0,1d0,2d0,2d0,2d0,1d0,1d0,
     $     1d0,2d0,2d0,2d0,1d0,1d0,1d0,2d0,2d0,2d0,2d0,2d0,2d0,4d0,4d0,
     $     4d0,2d0,2d0,2d0,4d0,4d0,4d0,2d0,2d0,2d0,4d0,4d0,4d0/
      IF(IOPT.EQ.1) THEN
         DO 30 I=1,3
            I1=IJV(I,1)
            I2=IJV(I,2)
            T2(I1,I2)=T1(I)
            T2(I2,I1)=T1(I)
 30      CONTINUE
         DO 31 I=4,6
            I1=IJV(I,1)
            I2=IJV(I,2)
            T2(I1,I2)=T1(I)/2d0
            T2(I2,I1)=T1(I)/2d0
 31      CONTINUE
      ENDIF
C     
      IF(IOPT.EQ.2) THEN
         DO 40 I=1,3
            I1=IJV(I,1)
            I2=IJV(I,2)
            T1(I)=T2(I1,I2)
 40      CONTINUE
         DO 41 I=4,6
            I1=IJV(I,1)
            I2=IJV(I,2)
            T1(I)=T2(I1,I2)*2d0
 41      CONTINUE
      ENDIF
C
      IF (IOPT.EQ.3) THEN
         DO 10 I=1,6
            I1=IJV(I,1)
            I2=IJV(I,2)
      DO 10 J=1,6
         J1=IJV(J,1)
         J2=IJV(J,2)        
            C4(I1,I2,J1,J2)=C2(I,J)*(1d0/fact(i,j))
c            write(*,*)'Cijkl=cij',i1,i2,j1,j2
c            write(*,*)'cij',i,j
c            write(*,*)'fac',(1/fact(i,j))
            C4(I2,I1,J1,J2)=C2(I,J)*(1d0/fact(i,j))
            C4(I1,I2,J2,J1)=C2(I,J)*(1d0/fact(i,j))
            C4(I2,I1,J2,J1)=C2(I,J)*(1d0/fact(i,j))
 10   CONTINUE
      ENDIF
C
      IF(IOPT.EQ.4) THEN
         DO 20 I=1,6
            I1=IJV(I,1)
            I2=IJV(I,2)
         DO 20 J=1,6
            J1=IJV(J,1)
            J2=IJV(J,2)
            C2(I,J)=fact(i,j)*C4(I1,I2,J1,J2)
c            write(*,*)'cij',i,j
c            write(*,*)'Cijkl=cij',i1,i2,j1,j2
c            write(*,*)'fac',fact(i,j)
            
 20      CONTINUE
      ENDIF

      RETURN 
      END
