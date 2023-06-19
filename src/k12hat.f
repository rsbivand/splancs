C Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
C Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
C R port: copyright 1998-2000 by Roger S. Bivand


      subroutine k12hat(x1,y1,n1,x2,y2,n2,xp,yp,np,s,ns,h12,h21)

      implicit double precision(a-h,o-z)

      include 'bounds.h'

      dimension x1(n1),y1(n1),x2(n2),y2(n2)
      dimension xp(np+1),yp(np+1),s(ns),h12(ns),h21(ns)


      area = plarea(xp,yp,np)

      pi=3.141592654d0

      tmax=(s(ns))**2  

      do i=1,ns
        h12(i)=0.0d0
        h21(i)=0.0d0
      end do
c--

      do i=1,n1                                                           
        xi=x1(i)                                                              
        yi=y1(i)                                                              
        do j=1,n2                                                           
          xj=xi-x2(j)                                                           
          yj=yi-y2(j)                                                           
          t=xj*xj+yj*yj                                                         
          if (t.lt.tmax) then                                                 
            t=dsqrt(t)   
            it=iplace(s,ns,t)
            wij=weight(xi,yi,t,xp,yp,np)
            h12(it)=h12(it)+wij
          end if
        end do
      end do                                                   

      do i=1,n2                                                           
        xi=x2(i)                                                              
        yi=y2(i)                                                              
        do j=1,n1                                                           
          xj=xi-x1(j)                                                           
          yj=yi-y1(j)                                                           
          t=xj*xj+yj*yj                                                         
          if (t.lt.tmax) then                                                 
            t=dsqrt(t)             
            it=iplace(s,ns,t)              
            wij=weight(xi,yi,t,xp,yp,np)
            h21(it)=h21(it)+wij
          end if
        end do
      end do
                                                   
      do i=2,ns                                                        
        i1=i-1                                                                
        h12(i)=h12(i)+h12(i1)                                                 
        h21(i)=h21(i)+h21(i1)                                                 
      end do

      dn11=dble(n1*(n1-1))                                                    
      dn22=dble(n2*(n2-1))                                                    
      dn12=dble((n1-1)*(n2-1))
      alpha=dble(n2)/dble(n1+n2)                                        

      do i=1,ns                                                        
        h12(i)=(alpha*h12(i)+(1.0d0-alpha)*h21(i))*area/dn12                    
      end do

c--

      return                                                                
      end

