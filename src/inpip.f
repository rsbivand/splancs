C Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
C Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
C R port: copyright 1998-2000 by Roger S. Bivand


      subroutine inpip(xpts,ypts,npts,xpoly,ypoly,npoly,lind)

      real*8 xpts(npts),ypts(npts)
      real*8 xpoly(npoly+1),ypoly(npoly+1)
C lind converted to integer RSB 140930
      integer lind(npts)

      do i=1,npts
        if(ipippa(xpts(i),ypts(i),xpoly,ypoly,npoly).eq.0) then
          lind(i)=0
        else
          lind(i)=1
        end if
      end do

      end
