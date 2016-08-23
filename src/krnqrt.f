C Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
C Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
C R port: copyright 1998-2000 by Roger S. Bivand

      subroutine krnqrt(x,y,n,xp,yp,np,h0,a1,a2,b1,b2,nx,ny
     &                 ,xgrid,ygrid,zgrid)
      implicit double precision (a-h,o-z)
      dimension x(n),y(n),xp(np+1),yp(np+1)
      dimension xgrid(nx),ygrid(ny),zgrid(nx,ny)
c
c     version for quartic kernel
c
c    returns x and y vectors, and z giving smoothed values.

      a=a2-a1
      b=b2-b1
c
      h=h0*dsqrt(2.0d0)
      h2=h*h
      rh2=1.0d0/h2
      const=0.954929658d0*rh2

      call xsort(x,y,n)
      xinc=a/dfloat(nx)
      x1=a1
      x2=a2
      yinc=b/dfloat(ny)
      y1=b1
      y2=b2
      xc=x1-0.5d0*xinc

      i1=1
      do 1 i=1,nx
      xc=xc+xinc
      xgrid(i)=xc
      d1=dmin1(xc-a1,a2-xc)
      xcmh=xc-h
      xcph=xc+h
    3 m=i1
      if (m.eq.n) goto 4
    5 xm=x(m)
      if (xm.ge.xcmh) goto 2
      i1=i1+1
      goto 3
    2 if (xm.ge.xcph) goto 4
      m=m+1
      if (m.eq.n) goto 4
      goto 5
    4 i2=m
      yc=b1-0.5d0*yinc
      do 6 j=1,ny
      yc=yc+yinc
      ygrid(j)=yc
      if(ipippa(xc,yc,xp,yp,np).ne.0)then
        d2=dmin1(yc-b1,b2-yc)
        ycmh=yc-h
        ycph=yc+h
        count=0.0d0
        do 7 m=i1,i2
        ym=y(m)
        if ((ym.lt.ycmh).or.(ym.ge.ycph)) goto 7
        d=(x(m)-xc)**2+(ym-yc)**2
        if (d.lt.h2) count=count+((1.0d0-d*rh2)**2)/ssarea(d1,d2,h)
    7   continue
        zgrid(i,j)=const*count
      else
        zgrid(i,j)=-1
      end if
    6 continue
    1 continue

      end
      function ssarea(d1t,d2t,h)
      implicit double precision (a-h,o-z)
c ouch. ssarea was modifying its args. bad news 
      d1=d1t
      d2=d2t
      if ((d1.ge.h).and.(d2.ge.h)) goto 1
      d1=dmin1(1.0d0,d1/h)
      d2=dmin1(1.0d0,d2/h)
      theta1=dacos(d1)
      theta2=dacos(d2)
      if ((d1*d1+d2*d2).lt.(h*h)) goto 5
      ssarea=1.0d0-0.318309886d0*(theta1+theta2)
     1      +2.0d0*(arzz(d1,theta1)+arzz(d2,theta2))
      return
    1 ssarea=1.0d0
      return
    5 theta3=datan(d1/d2)
      theta4=1.50796327d0-theta3
      ssarea=0.75d0-0.159154943d0*(theta1+theta2)
     1      +arzz(d1,theta1)+arzz(d2,theta2)
     2      +arzz(d2,theta3)+arzz(d1,theta4)
      return
      end
      function arzz(d,theta)
      implicit double precision (a-h,o-z)
      t=dtan(theta)
      d2=d*d
      d4=d2*d2
      d6=d2*d4
      arzz=t*d2*(1.0d0-d2+d4/3.0d0)+
     1  (t**3)*d4*(d2/3.0d0-1.0d0)/3.0d0
     2  +(t**5)*d6/15.0d0
      arzz=arzz*0.477464829d0
      return
      end
      subroutine xsort(x,y,n)                                                  
c                                                                           
c     shellsort algorithm                                                   
c     n     : number of elements to be sorted                               
c     x     : on input an array of dimension at least n containing          
c             real numbers                                                  
c             on output first n elements of x are sorted from smallest      
c             to largest                                                    
c     y     : on input an array of dimension at least n containing
c             real numbers
c             on output first n elements of y are shuffled to
c             correspond to sorting on x
c                                                                           
      implicit double precision (a-h,o-z)
      dimension x(n),y(n)                                                        
      i=1                                                                   
    1 i=i+1                                                                 
      if (i.le.n) goto 1                                                    
      m=i-1                                                                 
    2 m=m/2                                                                 
      if (m.eq.0) return                                                    
      k=n-m                                                                 
      do 4 j=1,k                                                            
      kk=j                                                                  
    3 if (kk.lt.1) goto 4                                                   
      if (x(kk+m).ge.x(kk)) goto 4                                          
      w=x(kk+m)                                                             
      x(kk+m)=x(kk)                                                         
      x(kk)=w                                                               
      v=y(kk+m)
      y(kk+m)=y(kk)
      y(kk)=v
      kk=kk-m                                                               
      goto 3                                                                
    4 continue                                                              
      goto 2                                                                
      end                                                              
