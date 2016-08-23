C Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
C Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
C R port: copyright 1998-2000 by Roger S. Bivand

      subroutine khvmat(x,y,n,n1,n2,xp,yp,np,s,ns,Amat,Bvec,Covmat)

      implicit double precision (a-h,o-z)

      include 'bounds.h'

      dimension x(n),y(n),s(ns+1),Covmat(ns,ns),xp(np+1),yp(np+1)
      dimension Amat(n,ns),Bvec(ns)
c
c     sample variances  of K11(s)-K22(s) under random labelling
c     Algorithm A G Chetwynd  P J Diggle

      area=plarea(xp,yp,np)
      area=area**2

      dn1=dfloat(n1)
      dn2=dfloat(n2)
      dn=dfloat(n)
      dn12=dn1*(dn1-1.0d0)
      dn22=dn2*(dn2-1.0d0)
      do i=2,n
         do j=1,i-1
            dij=(x(i)-x(j))**2+(y(i)-y(j))**2
            if (dij.le.(s(ns))**2) then
               rdij=dsqrt(dij)
               wij=weight(x(i),y(i),rdij,xp,yp,np)
               wji=weight(x(j),y(j),rdij,xp,yp,np)
               vij=(wij+wji)/2
               m=iplace(s,ns,rdij)
c               call intpr('m= ',-1,m,1)
               do k=m,ns
                  Amat(i,k)=Amat(i,k)+vij
                  Bvec(k)=Bvec(k)+vij**2
                  Amat(j,k)=Amat(j,k)+vij
               end do

            end if
         end do
      end do
      do ise=1,ns
         bigv=0
         do i=1,n
            bigv=bigv+Amat(i,ise)
         end do
         do k=1,ise
            bigw=0
	    bigz=0
            do i=1,n
               bigw=bigw+Amat(i,k)
  	       comp=Amat(i,ise)*Amat(i,k)
	       bigz=bigz+comp
	    end do
            bigx=Bvec(k)*2.
            c3=(bigw*bigv)-4*(bigz)+2*(bigx)
            c2=4*(bigz-bigx)
	    c1=2*bigx

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
        Covmat(ise,k)=dmu1/(dn12*dn12)+dmu2/(dn22*dn22)-
     &            2.0d0*dmu3/(dn12*dn22)
        Covmat(ise,k)=area*(Covmat(ise,k))
c      call dblepr('Covmat = ',-1,Covmat(ise,k),1)
         end do
      end do
      return
      end

