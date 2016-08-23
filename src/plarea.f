C Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
C Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
C R port: copyright 1998-2000 by Roger S. Bivand

c-------------------------------------------------------------------------
      function plarea(xp,yp,np)
c-------------------------------------------------------------------------
c
c find the area of the polygon defined by points in xp,yp
c
      implicit double precision (a-h,o-z)

      dimension xp(np+1),yp(np+1)

      totare=0

      do is=1,np

        x1=xp(is)
        y1=yp(is)

        if(is.eq.np)then
          x2=xp(1)
          y2=yp(1)
        else
          x2=xp(is+1)
          y2=yp(is+1)
        end if

c Find the area of the trapezium
        totare = totare + (x2-x1)*(y2+y1)/2.0

      end do

c return a positive value
      plarea = abs(totare)

      end

