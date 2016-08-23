C Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
C Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
C R port: copyright 1998-2000 by Roger S. Bivand

      subroutine stsecal(x,y,n,xp,yp,np,s,ns,Smat,Svec,z,tlow,tupp
     &     ,t,nt,Tmat,Tvec,long,Covmat)

      implicit double precision (a-h,o-z)

      include 'bounds.h'

      dimension x(n),y(n),s(ns),Covmat(long,long),xp(np+1),yp(np+1)
      dimension Smat(n,ns),Svec(ns),z(n),t(nt),Tmat(n,nt),Tvec(nt)
c     
c     The covariance matrix for the space time interaction
c     Algorithm A G Chetwynd  P J Diggle

      area=plarea(xp,yp,np)
      area=area**2

      capT=tupp-tlow
      capT=capT**2
      
      dn=dfloat(n)
      c1=(1/(dn*(dn-1.0d0))**2)
      c2=(1/(dn*(dn-1.0d0)*(dn-2.0d0)*(dn-3.0d0)))
      c3=4*(1/(dn*(dn-1.0d0)*(dn-2.0d0)))
      c4=2*(1/(dn*(dn-1.0d0)))
      
      
      const=area*capT*c1

      smax=s(ns)
      tmax=t(nt)

      do i=2,n
         do j=1,i-1
            dij=(x(i)-x(j))**2+(y(i)-y(j))**2
            if (dij.le.(smax)**2) then
               rdij=dsqrt(dij)
               wij=weight(x(i),y(i),rdij,xp,yp,np)
               wji=weight(x(j),y(j),rdij,xp,yp,np)
               w=(wij+wji)/2
               m=iplace(s,ns,rdij)
c               call intpr('m= ',-1,m,1)
               do k=m,ns
                  Smat(i,k)=Smat(i,k)+w
                  Svec(k)=Svec(k)+w**2
                  Smat(j,k)=Smat(j,k)+w
               end do
            end if
         end do
      end do
      do i=1,ns
c     call dblepr('Svec=',-1,2*Svec(i),1)
      end do
c      call dblepr('tmax',-1,tmax,1)
      do i=2,n
         do j=1,i-1
            uij=dabs(z(i)-z(j))
c            call dblepr('uij',-1,uij,1)
            if (uij.le.tmax) then
               vij=1.0d0
	       if ((((z(i)-tlow).le.uij).or.((tupp-z(i)).le.uij))) vij=2.0d0
c               call dblepr('vij=',-1,vij,1)
               vji=1.0d0
               if ((((z(j)-tlow).le.uij).or.((tupp-z(j)).le.uij))) vji=2.0d0
c               call dblepr('vji=',-1,vji,1)
               v=(vij+vji)/2
               m=iplace(t,nt,uij)
c               call intpr('m= ',-1,m,1)
               do k=m,nt
                  Tmat(i,k)=Tmat(i,k)+v
                  Tvec(k)=Tvec(k)+v**2
                  Tmat(j,k)=Tmat(j,k)+v
               end do
            end if
         end do
      end do
      do i=1,nt
c      call dblepr('Tvec=',-1,2*Tvec(i),1)
      end do
      irow=0

      
      do ia=1,ns
         w1=0
	 do ix=1,n
	    w1=w1+Smat(ix,ia)
         end do
c         call dblepr('w1',-1,w1,1)
	 sval=s(ia)

	 do ib=1,nt
	    irow=irow+1
	    jcol=0
   	    v1=0
	    do ix=1,n
	       v1=v1+Tmat(ix,ib)
            end do
c           call dblepr('v1',-1,v1,1)
            tval=t(ib)

            do ic =1,ns
	       w2=0
	       wdash=0
	       sdash=s(ic)
	       do ix=1,n
		  wdash=wdash+Smat(ix,ic)
		  comp=Smat(ix,ia)*Smat(ix,ic)
		  w2=w2+comp
	       end do
c                 call dblepr('wdash',-1,wdash,1)
c                 call dblepr('w2',-1,w2,1)
	       if (sval.le.sdash) then
                 w3=2*Svec(ia)
 	       else
		 w3=2*Svec(ic)
	       end if
c               call dblepr('w3',-1,w3,1)
               do id=1,nt
		  jcol=jcol+1
 	       if (jcol.le.irow) then
		     tdash=t(id)
                     vdash=0
  		     v2=0
  		     do ix=1,n
          	        vdash=vdash+Tmat(ix,id)
                        comp=Tmat(ix,ib)*Tmat(ix,id)
                        v2=v2+comp
		     end do
c                  call dblepr('vdash',-1,vdash,1) 
c                  call dblepr('v2',-1,v2,1)
                     if (tval.le.tdash) then
       		        v3=2*Tvec(ib)
		     else
    			v3=2*Tvec(id)
		     end if 
c			call dblepr('v3',-1,v3,1)
                     A=((w1*wdash) - (4*w2) +(2*w3))*
     &               ((v1*vdash) - (4*v2) + (2*v3))
                     B=(w2-w3)*(v2-v3)
		     C=w3*v3
		     D=w1*v1*wdash*vdash              
c		     write(*,*)'row ',irow,'column ',jcol
c	             call intpr('row= ',-1,irow,1)
c		     call intpr('col= ',-1,jcol,1)
                     Covmat(jcol,irow)=const*(c2*A + c3*B + c4*C - c1*D)
                  else
                     goto 190
                  end if		
               end do
            end do
  190   continue
        end do
      end do

      return
      end


