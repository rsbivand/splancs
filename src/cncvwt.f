C Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
C Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
C R port: copyright 1998-2000 by Roger S. Bivand


      function cncvwt(x,y,r,xp,yp,np)
c
c compute the weight given to a point at x,y according to how much
c of a circle of radius r is inside the bounding polygon
c 
c
      implicit double precision (a-h,o-z)
      include 'bounds.h'
      dimension xp(np+1),yp(np+1)
      parameter(pi=3.141592654d0)
c store circle/poly intersections here
      parameter(maxcrs=40)
      dimension cross(maxcrs+1)
      parameter(tiny=1.0e-7)
c set count of crossing points to zero
      ncross = 0

c first loop over the boundary and find the crossing points
      do ic=1,np
c work with the trial point at origin
        x1=xp(ic)-x
        y1=yp(ic)-y
        x2=xp(ic+1)-x
        y2=yp(ic+1)-y

        cx=x2-x1
        cy=y2-y1
 
c these are the coefficients of the quadratic giving the intercept of
c line and circle.
        a=cx*cx+cy*cy
        b=2*(x1*cx+y1*cy)
        c=x1*x1+y1*y1-r*r

c find out if real solutions exist...
        b2m4ac=b*b-4*a*c

c ... and if they do, find them.
        if (b2m4ac.ge.0) then
          t1=(-b+sqrt(b2m4ac))/(2*a)
          t2=(-b-sqrt(b2m4ac))/(2*a)

c see if the solutions lie in the line segments
          if ((t1.gt.tiny).and.(t1-1.0.le.tiny)) then
             ncross=ncross+1
c find the angle to this point on thecircle
             ctemp=atan2(y1+t1*cy,x1+t1*cx)
             if(ctemp.lt.0)ctemp=2*pi+ctemp
             cross(ncross)=ctemp
c check crossing of circle with vertex
           else if (abs(t1).le.tiny)then
c compare this polygon segment's direction with that of the previous one
             nprev = (mod((ic+ (np-2)),np)+1)
             x0 = xp(nprev) - x
             y0 = yp(nprev) - y
             idp1 = isig8((x2-x1)*x1+ (y2-y1)*y1,tiny)
             idp2 = isig8((x1-x0)*x1+ (y1-y0)*y1,tiny)
c see if the polygon passes through the circle here
             if ((idp1-idp2).ne.1 .and.
     +        abs(idp1+idp2).ne.2) then
               ncross = ncross + 1
               ctemp = atan2(y1+t1*cy,x1+t1*cx)
               if (ctemp.lt.0.0) ctemp = 2*pi + ctemp
               cross(ncross) = ctemp
             end if
          end if

          if ((t2.gt.tiny).and.(t2-1.0.lt.tiny)) then
             ncross=ncross+1
             ctemp=atan2(y1+t2*cy,x1+t2*cx)
             if(ctemp.lt.0)ctemp=2*pi+ctemp
             cross(ncross)=ctemp
c check crossing of circle with vertex
           else if (abs(t2).le.tiny)then
c compare this polygon segment's direction with that of the previous one
             nprev = (mod((ic+ (np-2)),np)+1)
             x0 = xp(nprev) - x
             y0 = yp(nprev) - y
             idp1 = isig8((x2-x1)*x1+ (y2-y1)*y1,tiny)
             idp2 = isig8((x1-x0)*x1+ (y1-y0)*y1,tiny)
c see if the polygon passes through the circle here
             if ((idp1-idp2).ne.1 .and.
     +        abs(idp1+idp2).ne.2) then
               ncross = ncross + 1
               ctemp = atan2(y1+t2*cy,x1+t2*cx)
               if (ctemp.lt.0.0) ctemp = 2*pi + ctemp
               cross(ncross) = ctemp
             end if
          end if
        end if
      end do

c now we have all the crossing point angles stored in cross(1:ncross)

c if ncross = 0 then the total angle within the poly is 2*pi unless the
c circle is large and spans the polygon. this should be checked
c beforehand so it's okay to assume 2*pi here.

      if (ncross.eq.0) then
        totang=2*pi
      else

c sort into ascending order
        call sort2(cross,ncross)

c fix the ncross+1'th element to be the first plus 2pi so that the
c   list is circular...
        cross(ncross+1)=cross(1)+2*pi

c check that the number of crossings is even - if not then error.
        if (mod(ncross,2).ne.0) then
          cncvwt=-1
          return
        end if
c now find a nice spot to do the point-in-poly search
        sepmax=0.0
        icm=0

        do ic=1,ncross
          if (cross(ic+1)-cross(ic).gt.sepmax) then
            sepmax=cross(ic+1)-cross(ic)
            icm=ic
          end if
        end do
  
c icm is now the index of the crossing with the largest gap between it
c and the next crossing point

c test for point in poly of the point on the circle between these points
        angtes=(cross(icm)+cross(icm+1))/2.

        xtest=x+r*cos(angtes)
        ytest=y+r*sin(angtes)

c find out if test point is in the polygon boundary
        linpol=ipippa(xtest,ytest,xp,yp,np)

c find the total angle between (odd-even) crossings (i.e. 1-2 + 3-4 + ...
        totang = 0.
        do ic=1,ncross-1,2
          totang = totang + (cross(ic+1)-cross(ic))
        end do

c If the point we tested for p-i-p was on an odd-even
c section and was in the poly, then totang is the amount of circle inside the
c polygon. if the point was outside the polygon, then we need to subtract
c totang from 2*pi radians to get the angle inside the polygon. conversely,
c if the point tested was between even-odd crossings and outside the polygon,
c then totang is the angle we want, and if inside the polygon then again we 
c have to do 2*pi-totang

        if ( (((mod(icm,2).eq.1).and.(linpol.eq.0))  .or.
     &        ((mod(icm,2).eq.0).and.(linpol.eq.1)) ) ) then
          totang = 2*pi-totang
        end if

      end if
c now totang is the angle contained in the polygon

c weight is proportion of total angle in the poly
      cncvwt = (2*pi)/(totang)
      return
      end

      integer function isig8(value,tiny)
c return the sign (+1,0,-1) of a value
            double precision tiny,value
          if (value.gt.tiny) then
            isig8 = 1
          else if (value.lt.-tiny) then
            isig8 = -1
          else
            isig8 = 0
          end if
        return
      end
      
