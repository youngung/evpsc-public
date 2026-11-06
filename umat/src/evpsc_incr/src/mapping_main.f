      PROGRAM mapping_main

      implicit none

      integer rowmx,i,j,num_rows, num_select, selmx, UR_mapping
      parameter (rowmx=300000)
      parameter (selmx=1000)
      real*8 arr_ebsd(5,rowmx), arr_select(5,selmx),x,y,r

      UR_mapping = 999
      
      x=350.
      y=800.
      r=5.

      open(unit=UR_mapping,file='Mg10Gd_Rx250.ang',status='OLD')
      
      CALL read_ebsd(
c     in
     $     UR_mapping,
c     out
     $     arr_ebsd,num_rows)

      close(UR_mapping)


      
      CALL select_points_from_ebsd(
c     in
     $     arr_ebsd,x,y,r,num_rows,
c     out
     $     arr_select,num_select)

      END PROGRAM mapping_main
