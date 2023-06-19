C Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
C Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
C R port: copyright 1998-2000 by Roger S. Bivand

      subroutine stkhat(x1,y1,t1,n1,xc,yc,nc,s,ns,t,nt,tlow,thigh,
     &                  hs,ht,hst)
c
c     space-time K-function calculation
c     assuming data on an arbitrary polygon
c  
c     See Diggle(1983, chapters 4,5) for methodological details.
c
      implicit double precision(a-h,o-z)
      dimension x1(n1),y1(n1),t1(n1),hs(ns),ht(nt),hst(ns,nt),
     +          t(nt),s(ns),
     1          xc(nc+1),yc(nc+1)

           pi=3.141592654d0
c

      area=plarea(xc,yc,nc)

      call ksthat(nc,xc,yc,x1,y1,t1,n1,
     1     s,ns,t,nt,hs,ht,hst,area,tlow,thigh)             

      end



      subroutine ksthat(nc,xc,yc,x1,y1,t1,n1,                 
     1     s,ns,t,nt,hs,ht,hst,area,tlow,thigh)                   
c                                                                           
c     space-time reduced second moment measure                               
c     data stored as n1 locations (x1,y1,t1)                          
c     in convex polygon. reduced second moment measures                           
c     returned as arrays hs,ht,hst, tabulation in                           
c     nscell by ntcell cells of size hscell by htcell,starting from                          
c     zero.                                                                 
c                                                                           
      implicit double precision (a-h,o-z)                                             
      dimension x1(n1),y1(n1),t1(n1),hs(ns),                      
     1          ht(nt),hst(ns,nt),t(nt),s(ns),
     2          xc(nc+1),yc(nc+1)

      smax=s(ns)
      tmax=t(nt)

c      call intpr('nt ',-1,nt,1)

      do 1 i=1,ns
      hs(i)=0.0d0                                                           
    1 continue
      do 101 j=1,nt
      ht(j)=0.0d0
      do 102 i=1,ns
      hst(i,j)=0.0d0
  102 continue
  101 continue
      do 2 i=2,n1                                                           
      i1=i-1                                                                
      xi=x1(i)                                                              
      yi=y1(i)                                                              
      ti=t1(i)
      do 2 j=1,i1                                                           
      xj=xi-x1(j)                                                           
      yj=yi-y1(j)                                                           
      sj=dsqrt(xj*xj+yj*yj)
      if (sj.ge.smax) goto 1002                                                 
      is=iplace(s,ns,sj)                                                    
      wij=weight(xi,yi,sj,xc,yc,nc)
      wji=weight(x1(j),y1(j),sj,xc,yc,nc)
      hs(is)=hs(is)+wij+wji
 1002 continue
      tj=dabs(ti-t1(j))
      if (tj.ge.tmax) goto 2                                                 
      it=iplace(t,nt,tj)
c      call intpr('i  ',-1,i,1)
c      call intpr('j  ',-1,j,1)
c      call intpr('it ',-1,it,1)
      vij=1.0d0
      if (((ti-tlow).le.tj).or.((thigh-ti).le.tj)) vij=2.0d0
      vji=1.0d0
      if (((t1(j)-tlow).le.tj).or.((thigh-t1(j)).le.tj)) vji=2.0d0
      ht(it)=ht(it)+vij+vji
      if (sj.le.smax) hst(is,it)=hst(is,it)+wij*vij+wji*vji
    2 continue                                                              
      do 3 is=2,ns
      is1=is-1
    3 hs(is)=hs(is)+hs(is1)
c      call dblepr('pre-acc ',-1,ht,nt)
      do 4 it=2,nt
      it1=it-1
    4 ht(it)=ht(it)+ht(it1)
c      call dblepr('post-acc ',-1,ht,nt)
      do 5 is=1,ns
      do 6 it=2,nt
      it1=it-1
      hst(is,it)=hst(is,it)+hst(is,it1)
    6 continue
    5 continue
      do 7 it=1,nt
      do 8 is=2,ns
      is1=is-1
      hst(is,it)=hst(is,it)+hst(is1,it)
    8 continue
    7 continue
      dn11=dble(n1*(n1-1))                                                    
      do 17 i=1,ns
      hs(i)=hs(i)*area/dn11                                                   
   17 continue                                                              
      do 18 j=1,nt
      ht(j)=ht(j)*(thigh-tlow)/dn11
      do 108 i=1,ns
  108 hst(i,j)=hst(i,j)*area*(thigh-tlow)/dn11
   18 continue
      return                                                                
      end                                                                   


