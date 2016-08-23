C Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
C Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
C R port: copyright 1998-2000 by Roger S. Bivand

      subroutine dosekh(x,y,n,n1,n2,xp,yp,np,s,ns,se)

      implicit double precision (a-h,o-z)

      include 'bounds.h'

      dimension x(n),y(n),s(ns),se(ns),xp(np+1),yp(np+1)

c
c     sample standard errors of K11(s)-K22(s) under random labelling
c

      area=plarea(xp,yp,np)

      dn1=dfloat(n1)
      dn2=dfloat(n2)
      dn=dfloat(n)
      dn12=dn1*(dn1-1.0d0)
      dn22=dn2*(dn2-1.0d0)


      do ise=1,ns
c        write(*,*)'Doing secal for cell ',icell,' of ',ncell
        d2crit=s(ise)**2
        bigw=0.0d0
        s1=0.0d0
        s2=0.0d0
        do i=2,n
          iminus=i-1
          do j=1,iminus
            dij=(x(i)-x(j))**2+(y(i)-y(j))**2
            if (dij.le.d2crit) then
              dij=dsqrt(dij)
              wij=weight(x(i),y(i),dij,xp,yp,np)
              wji=weight(x(j),y(j),dij,xp,yp,np)
              vij=wij+wji
              bigw=bigw+vij
              s1=s1+vij**2
            end if
          end do
        end do
        do i=1,n
          t1=0.0d0
          do j=1,n
            if (i.ne.j) then
              dij=(x(i)-x(j))**2+(y(i)-y(j))**2
              if (dij.le.d2crit) then
                dij=dsqrt(dij)
                wij=weight(x(i),y(i),dij,xp,yp,np)
                wji=weight(x(j),y(j),dij,xp,yp,np)
                vij=wij+wji
                t1=t1+vij
              end if
            end if
          end do
          s2=s2+t1**2
        end do
        c1=s1
        c2=s2-2.0d0*s1
        c3=bigw*bigw+s1-s2
        dmu1=c1*(dn1/dn)*((dn1-1.0d0)/(dn-1.0d0))
     &    +c2*(dn1/dn)*((dn1-1.0d0)/(dn-1.0d0))*((dn1-2.0d0)/(dn-2.0d0))
     &    +c3*(dn1/dn)*((dn1-1.0d0)/(dn-1.0d0))*
     &        ((dn1-2.0d0)/(dn-2.0d0))*((dn1-3.0d0)/(dn-3.0d0))
        dmu2=c1*(dn2/dn)*((dn2-1.0d0)/(dn-1.0d0))
     &    +c2*(dn2/dn)*((dn2-1.0d0)/(dn-1.0d0))*((dn2-2.0d0)/(dn-2.0d0))
     &    +c3*(dn2/dn)*((dn2-1.0d0)/(dn-1.0d0))*
     &        ((dn2-2.0d0)/(dn-2.0d0))*((dn2-3.0d0)/(dn-3.0d0))
        dmu3=c3*(dn1/dn)*((dn1-1.0d0)/(dn-1.0d0))*(dn2/(dn-2.0d0))*
     &        ((dn2-1.0d0)/(dn-3.0d0))
        se(ise)=dmu1/(dn12*dn12)+dmu2/(dn22*dn22)-
     &            2.0d0*dmu3/(dn12*dn22)
        se(ise)=area*dsqrt(se(ise))
      end do

      return
      end
