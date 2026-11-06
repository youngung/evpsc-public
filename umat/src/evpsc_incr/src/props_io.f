      subroutine props_io(iopt,nprops,props,edotref)
c     arguments list
c     iopt: option to determine the behavior of the subroutine
c         ( 0): Test if nprops is properly sized.
c         ( 1): Save EVPSC constants to props
c         ( 2): Load EVPSC constants from props
c         ( 3): Save entire props to a file
c         ( 4): Load entire props from a file

c     props:  properties of Abaqus
c     nprops: # of properties of Abaqus

      implicit none
      integer, intent(in)::iopt,nprops
      real*8, intent(inout):: edtoref


      return
      end subroutine props_io(iopt)
