c     Notably Abaqus does not allow to use 'stop' statement, instead it uses subroutine called 'xit'.
c     Mutexlock and Mutexunlock are two subroutines used in umat to hold 'threaded' calculation.

c     In this file, fake 'xit', 'mutexlock', and 'mutexunlock' are provided so that these names
c     can be used in both Stand-alone and UMAT. Notably, in the case of former, these are not genuine routines.

      subroutine xit
      write(*,*)'** Requested to exit from fake xit'
c      stop -1
c     stop statement is not compatible with win2019 server at lab (Dec 26, 2023)
      call exit(-1)
      end subroutine xit

      subroutine mutexlock(id)
      integer id
c      id=id ! to fool
      write(*,*)id
      end subroutine mutexlock

      subroutine mutexunlock(id)
      integer id
c      id=id ! to fool
      write(*,*)id
      end subroutine mutexunlock
