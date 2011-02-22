C Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
C Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
C R port: copyright 1998-2000 by Roger S. Bivand



      function iplace(s,ns,t)
c
c which of the variable width bins s is t in?
c
      implicit real*8 (a-h,o-z)

      dimension s(ns)

      do ib=1,ns
        if(s(ib).ge.t)then
          iplace=ib
          return
        end if
      end do
c
c if it is outside the range of s
c
      iplace=ns+1

      return
      end
