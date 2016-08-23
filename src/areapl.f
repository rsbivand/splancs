C Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
C Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
C R port: copyright 1998-2000 by Roger S. Bivand

      subroutine areapl(xp,yp,np,plar)
c
c return the area of the polygon as one of the arguments of the sub
c 
      implicit double precision(a-h,o-z)
      dimension xp(np+1),yp(np+1)
      plar=plarea(xp,yp,np)
      return
      end
      
