c     subroutine write_stress_strain for EVPSC models.

c     if istep==0, write only the heads.
c     else, writhe stress, strain and other quantities.

      SUBROUTINE WRITE_STRESS_STRAIN_evp (
c          in
     $     iu,ISTEP,tincr,
     $     epsvm,svm,epstot,temperat,LIJBARc,scauchy,press,presswgt,
     $     plwork,
     $     epstot_el,edot_el
     $     )
      implicit none
      integer,intent(in)::iu,istep
      real*8,intent(in)::tincr,epsvm,svm,epstot(3,3),temperat,
     $     LIJBARc(3,3),
     $     scauchy(3,3),press,presswgt,plwork,epstot_el(3,3),
     $     edot_el(3,3)
c     locals
      integer i,j,k,n,m
      integer IJV(6,2)
      DATA ((IJV(N,M),M=1,2),N=1,6)/1,1,2,2,3,3,2,3,1,3,1,2/
      if(istep.eq.0) then
         write(IU,'(14a15)',advance='no')'EVM','SVM','E11','E22','E33',
     $        'E23','E13','E12','SCAU11','SCAU22','SCU33','SCAU23',
     $        'SCAU13','SCAU12'
         write(IU,'(26a15)')'L11','L12','L13','L21','L22',
     $        'L23','L31','L32','L33','EE11','EE22','EE33','EE23',
     $        'EE13','EE12',
     $        'EDE11','EDE22','EDE33','EDE23','EDE13','EDE12',
     $        'tincr','press','<p>','TEMP','plwork'
         return
      endif

      write(IU,'(2e15.6e3)',advance='no')epsvm,svm
      do k=1,6
         i=ijv(k,1)
         j=ijv(k,2)
         write(IU,'(e15.6e3)',advance='no')epstot(i,j)
      enddo
      do k=1,6
         i=ijv(k,1)
         j=ijv(k,2)
         write(IU,'(e15.6e3)',advance='no')scauchy(i,j)
      enddo
c     upto 2+6+6=14
      do i=1,3
         do j=1,3
            write(IU,'(es15.6e3)',advance='no')LIJBARc(i,j)
         enddo
      enddo
      do k=1,6
         i=ijv(k,1)
         j=ijv(k,2)
         write(IU,'(es15.6e3)',advance='no')epstot_el(i,j)
      enddo
      do k=1,6
         i=ijv(k,1)
         j=ijv(k,2)
         write(IU,'(es15.6e3)',advance='no')edot_el(i,j)
      enddo
      write(IU,'(5es15.6e3)')tincr,press,presswgt,temperat,plwork
c     upto +9+6+6+5=26; total 14+26=40


      RETURN
      END subroutine write_stress_strain_evp
