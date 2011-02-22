C Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
C Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
C R port: copyright 1998-2000 by Roger S. Bivand

      subroutine krnqne(x,y,n,xp,yp,np,h0,a1,a2,b1,b2,nx,ny
     &                 ,xgrid,ygrid,zgrid)
      implicit real*8 (a-h,o-z)
      dimension x(n),y(n),xp(np+1),yp(np+1)
      dimension xgrid(nx),ygrid(ny),zgrid(nx,ny)
c
c     version for quartic kernel, no edge correction
c
c    returns x and y vectors, and z giving smoothed values.

c
      parameter(pi=3.14159265358979323846264338)
      
      xw=a2-a1
      yw=b2-b1

      xh=xw/dfloat(nx)
      yh=yw/dfloat(ny)
      do ix=1,nx
        xgrid(ix)=.5*xh+a1+xh*(ix-1)
      end do
      do iy=1,ny
        ygrid(iy)=.5*yh+b1+yh*(iy-1)
      end do

      do ix=1,nx
        xc=xgrid(ix)
        do iy=1,ny
          yc=ygrid(iy)
          if(ipippa(xc,yc,xp,yp,np).ne.0)then
            rk=0.0
            do ip=1,n
              xpo=x(ip)
              ypo=y(ip)
              u2=((xpo-xc)**2+(ypo-yc)**2)/(h0*h0)
              if(u2.lt.2)rk=rk+(1-u2/2)**2
            end do
            zgrid(ix,iy)=(1/h0)*rk
          else
            zgrid(ix,iy)=-1
          end if
        end do
      end do
      
      return
      end

