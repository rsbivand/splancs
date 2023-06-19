C Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
C Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
C R port: copyright 1998-2002 by Roger S. Bivand

      subroutine trykh(x,y,n,xp,yp,np,s,ns,hkhat,icounts,hkhats,nptns)

      implicit double precision(a-h,o-z)

      include 'bounds.h'

      dimension x(n),y(n),xp(np+1),yp(np+1),s(ns),hkhat(ns)
      dimension icounts(nptns), hkhats(nptns)

      area = plarea(xp,yp,np)

      pi=3.141592654d0

      tmax=(s(ns))**2  

      do i=1,ns
        hkhat(i)=0.0d0
      end do


      do i=2,n
        i1=i-1                                                                
        xi=x(i)                                                              
        yi=y(i)                                                              
        do  j=1,i1                                                           
          xj=xi-x(j)                                                           
          yj=yi-y(j)                                                           
          t=xj*xj+yj*yj                                                         
          if (t.lt.tmax) then                                                 

            t=dsqrt(t)
            it=iplace(s,ns,t) 
 
            if(it.le.ns) then
              wij=weight(xi,yi,t,xp,yp,np)
              wji=weight(x(j),y(j),t,xp,yp,np)

              hkhat(it)=hkhat(it)+wij+wji
              ipos=n*(it-1)
              hkhats(i+ipos)=hkhats(i+ipos) + wij
              hkhats(j+ipos)=hkhats(j+ipos) + wji
              icounts(i+ipos)=icounts(i+ipos) + 1
              icounts(j+ipos)=icounts(j+ipos) + 1
            end if
          end if
        end do
      end do

      do i=2,ns
        hkhat(i)=hkhat(i)+hkhat(i-1)
        do j=1,n
          jj=j+(n*(i-1))
          jj1=j+(n*(i-2))
          hkhats(jj)=hkhats(jj)+hkhats(jj1)
        end do
      end do

      dn=dble(n)*dble(n-1)
      adn=area/dn

      do i=1,ns                                              
        hkhat(i)=hkhat(i)*adn
        do j=1,n
          jj=j+(n*(i-1))
          hkhats(jj)=hkhats(jj)*adn
        end do
      end do

      return                                                                
      end
