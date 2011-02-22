C Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
C Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
C R port: copyright 1998-2000 by Roger S. Bivand

      subroutine mse2d(x,y,n,a1,a2,b1,b2,nsmse,hsmse,amse,t)

      implicit real*8 (a-h,o-z)

      dimension x(n),y(n),t(nsmse),g(1000),hd(1000),amse(nsmse)
      common /anu/ hcell,h(1000),aval(1000),ncell
      pi=3.141592654

      b=b2-b1
      a=a2-a1
      ncell=nsmse*2
      hcell=hsmse
      do 10 i=1,n
      x(i)=x(i)-a1
      y(i)=y(i)-b1
   10 continue
      call khat(n,a,b,x,y)
      al1=(a*b)/float(n)
      nc2=ncell/2
      do 1 i=1,nc2
      g(i)=sqrt(h(i)/pi)
      t(i)=hcell*i
      hd(i)=h(i)-t(i)*t(i)*pi
      d2=g(i)-t(i)
      ic2=i*2
      sum=h(1)*form(t(i),hcell)
      do 11 jc=2,ic2
      r=hcell*float(jc)
      jc1=jc-1
      sum=sum+(h(jc)-h(jc1))*form(t(i),r)
   11 continue
      sum=sum/(al1*al1)
      amse(i)=(al1-2.0*h(i))/(pi*t(i)*t(i))+
     &            al1*al1*0.1013211*sum/(t(i)**4)
    1 continue
      
      end

      subroutine khat(n,a,b,x,y)
      implicit real*8 (a-h,o-z)                              
      dimension x(n),y(n)                                          
      common /anu/ hcell,h(1000),aval(1000),ncell
c                                                                           
c     reduced second moment measure, corresponding to empirical             
c     distribution function of weighted inter-event distances.              
c     ncell : no. of cells for tabulation of edf(enter)                     
c     hcell : width of each cell. lower limit of tabulation is zero (enter) 
c     h     : resulting estimate of k-function (output)                     
c                                                                           
      r=(float(ncell)*hcell)**2                                            
      hc1=1.0/hcell                                                       
      do 10 i=1,ncell                                                       
   10 h(i)=0.0
      do 20 i=2,n                                                           
      i1=i-1                                                                
      xi=x(i)                                                               
      yi=y(i)                                                               
      do 20 j=1,i1                                                          
      x1=xi-x(j)                                                            
      y1=yi-y(j)                                                            
      t=(x1*x1+y1*y1)                                                       
      if (t.ge.r) goto 20                                                   
      t=sqrt(t)                                                            
      it=1+int(t*hc1)                                                     
      h(it)=h(it)+fn2(xi,yi,t,a,b)+fn2(x(j),y(j),t,a,b)                     
   20 continue                                                              
      do 40 i=2,ncell                                                       
   40 h(i)=h(i)+h(i-1)                                                      
      f1=(a*b)/float(n*n)                                                  
      do 50 i=1,ncell                                                       
 50   h(i)=h(i)*f1                                                          
      return                                                                
      end                                                                   
      function fn2(x,y,t,a,b)  
      implicit real*8 (a-h,o-z)                                             
c                                                                           
c     weight function for khat                                              
c     result is reciprocal of proportion of circumference                   
c     of circle,centre (x,y) and radius t, which is                         
c     contained in (0,a)x(0,b)                                              
c     n.b. 1:t.le.0.5*min(a,b)                                              
c          2:fn2=1 if t.le.0.01                                             
c                                                                           
      w(p)=6.28318057/(6.28318057-p)                                    
      fn2=1.0
      if (t.le.0.01) return                                               
      x1=min(x,a-x)                                                       
      y1=min(y,b-y)                                                       
      if (t.le.min(x1,y1)) return                                         
      r1=sqrt(x1*x1+y1*y1)                                                 
      if (t.lt.r1) goto 1                                                   
      fn2=w(1.570796327+acos(x1/t)+acos(y1/t))                        
      return                                                                
    1 if (t.gt.y1) goto 2                                                   
      fn2=w(2.0*acos(x1/t))                                             
      return                                                                
    2 p=acos(y1/t)                                                        
      if (t.gt.x1) p=p+acos(x1/t)                                         
      fn2=w(2.0*p)                                                        
      return                                                                
      end                                                                   
      function form(t,r)
      implicit real*8 (a-h,o-z)
      if (r.lt.(2.0*t)) goto 10
      form=0.0
      return
   10 continue
      tsq=t*t
      form=2.0*tsq*acos(0.5*r/t)-r*sqrt(tsq-0.25*r*r)
      return
      end
