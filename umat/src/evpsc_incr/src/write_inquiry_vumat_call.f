c$$$	       if (ndi.eq.2) then
c$$$		  write(uwdev,'(6(i4,x),3e11.3,7(3e13.4e3,1x),(5e13.4e3),
c$$$     $69e13.4e3)')
c$$$     $                 noel,npt,jstep(1),kinc,NPH,ngr(NPH)
c$$$     $                 ,time(1),
c$$$     $                 time(2),dtime,stran(1),stran(2),stran(i12_comp),
c$$$     $                 dstran(1),dstran(2),dstran(i12_comp),stress(1),
c$$$     $                 stress(2),stress(i12_comp),ddsdde(1,1),
c$$$     $                 ddsdde(2,2),ddsdde(i12_comp,i12_comp),
c$$$     $                 ((drot(k,j),j=1,3),k=1,3),
c$$$     $                 epsvm,svmfe,dstran(1)/dtime,dstran(2)/dtime,
c$$$     $                 dstran(i12_comp)/dtime,
c$$$
c$$$     $                 (stran_e0_abq_voigt(j),j=1,6),
c$$$     $                 (stress_e0_abq_voigt(j),j=1,6),
c$$$     $                 (stran(j),j=1,ntens),
c$$$     $                 (stress(j),j=1,ntens),
c$$$     $                 ((dfgrd0(k,j),j=1,3),k=1,3),
c$$$     $                 ((dfgrd1(k,j),j=1,3),k=1,3),
c$$$     $                 ((rot_acum(k,j),j=1,3),k=1,3),
c$$$     $                 ((LIJBARc(k,j),j=1,3),k=1,3)
c	       elseif(ndi.eq.3) then
		  write(uwdev,'(6(i4,x),3e11.3,7(3e13.4e3,1x),(6e13.4e3),
     $69e13.4e3)')
     $     noel,npt,-1,-1,NPH,ngr(NPH),-1,
     $     -1,dtime,0d0,0d0,0d0,strainInc(km,1),
     $     strainInc(km,2),strainInc(km,i12_comp),
     $     stressNew(km,1),stressNew(km,2),stressNew(km,i12_comp),
     $
     $     0d0,0d0,0d0,!ddsdde(1,1),ddsdde(2,2),ddsdde(i12_comp,i12_comp)
c     $     ((drot(k,j),j=1,3),k=1,3),
     $     0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,
     $     epsvm,svmfe,
     $
     $     ! dstran(1)/dtime,dstran(2)/dtime,dstran(i12_comp)/dtime,dstran(3)/dtime,
     $     0d0,0d0,0d0,0d0,

     $     0d0,0d0,0d0,0d0,0d0,0d0, !(stran_e0_abq_voigt(j),j=1,6),
     $     0d0,0d0,0d0,0d0,0d0,0d0,!(stress_e0_abq_voigt(j),j=1,6),
c     $     (stran(j),j=1,ntens),

     $     0d0,0d0,0d0,0d0,0d0,0d0,
     $     (stressNew(km,j),j=1,ntens)
c     $     ((dfgrd0(k,j),j=1,3),k=1,3),
c     $     ((dfgrd1(k,j),j=1,3),k=1,3),
c     $     ((rot_acum(k,j),j=1,3),k=1,3),
c     $     ((LIJBARc(k,j),j=1,3),k=1,3)
