C Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
C Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
C R port: copyright 1998-2000 by Roger S. Bivand

      subroutine tribble(iflags,vars,npts,ndvars,ndpars,ncovars,
     &			iwhich,
     &                   pstart,pfin,
     &                   steps,
     &     reqmin,icount,kcode,dlogl)

      implicit real*8 (a-h,o-z)

      dimension iflags(npts),iwhich(ndvars)
      dimension vars(npts,ndvars+ncovars)
      dimension pstart(ndpars*2+ncovars+1)
      dimension steps(ndpars*2+ncovars+1)
      dimension pfin(ndpars*2+ncovars+1)

      nallpars = ndpars*2 + ncovars + 1

c      call intpr('ndvars',-1,ndvars,1)
c      call intpr('ndpars',-1,ndpars,1)
c      call intpr('ncovars',-1,ncovars,1)
c      call intpr('nallpars',-1,nallpars,1)

      call logem(pstart,nallpars,ndpars,1.00d0,0.01d0)

      konvge=5

      call nelmin(nallpars,pstart,pfin,ynewlo,reqmin,steps,konvge,
     & icount,iflags,vars,npts,nallpars,ndvars,iwhich,ndpars)

      kcode=konvge

      call unlog(pfin,nallpars,ndpars,1.00d0,0.01d0)

      call trblik(iflags,vars,npts,nallpars,ndvars,iwhich,ndpars
     &            ,pfin,dlogl)

      return
      end


	subroutine trblik(iflags,vars,npts,nallpars,ndvars,
     &      iwhich,ndpars,pars,dlogl)
	implicit real*8 (a-h,o-z)
	dimension iflags(npts),iwhich(ndvars)
        dimension vars(npts,nallpars-ndvars)
        dimension pars(nallpars)

c	call intpr('LOGLIK',-1,0,1)

c input
c iflags(1..npts) 1 == case, 0 = control
c vars(npts,nallpars) value of variates
c pars(nallpars)  (1..ndvars) alpha values (ndvars+1...ndvars*2) betas
c                 (ndvars*2+1 ... nallpars) thetas

c	output
c
c	   dlogl                   : log-likelihood
c

c	call intpr('nallpars',-1,nallpars,1)
c 	call intpr('NDVARS',-1,ndvars,1)

        result=0.0

c       call dblepr('EXtreme Vars',-1,vars(npts,nallpars-ndvars-1),1)

        do ic = 1,npts
          f=1
          do iv=1,ndvars
	  alpha=pars(iwhich(iv))
	  beta=pars(iwhich(iv)+ndpars)
          f=f*disfn(vars(ic,iv),alpha,beta)
          end do
c          call dblepr('f = ',-1,f,1)
c	call dblepr('DISVARS DONE, f = ',-1,f,1)
          ncovars=nallpars-(ndvars*2)-1
          sumcv=0.0
c        call intpr('Number of covariates = ',-1,ncovars,1)
          do iv = 1,ncovars
          sumcv=sumcv+(vars(ic,iv+ndvars)*pars(iv+ndpars*2))
          end do
          f=f*dexp(sumcv)
c	call intpr('YVARS DONE ',-1,0,1)
          rho=pars(nallpars)
c         call dblepr('rho = ',-1,rho,1)
          p=(rho*f)/(1+rho*f)
          if(iflags(ic).eq.0)p=1-p
          result = result + dlog(p)
        end do
	dlogl=result
c        call dblepr('LOGLIK',-1,result,1)
	return
	end

      subroutine nelmin(n,start,min,ynewlo,reqmin,step,konvge,icount,
     &              iflags,vars,npts,nallpars,ndvars,
     &              iwhich,ndpars)
      implicit real*8(a-h,o-z)
      dimension iflags(npts),iwhich(ndvars)
      dimension vars(npts,nallpars-ndvars)
       
c                                                                           
c     Simplex function minimisation procedure due to Nelder+Mead(1965),
c     as implemented by O'Neill(1971,Appl.Statist. 20, 338-45), with
c     subsequent comments by Chambers+Ertel(1974, 23, 250-1), Benyon(1976,
c     25, 97) and Hill(1978, 27, 380-2)
c
      real*8 min, rcoeff, ecoeff, ccoeff
      dimension p(20,21),pstar(20),p2star(20),pbar(20),y(20),        
     1   start(n),min(n),step(n)                               
c                                                                           
c      data rcoeff/1d0/,ecoeff/2d0/,ccoeff/5d-1/     
c       (changed to assignment RSB 16/9/2004; g77 data oddity)                  
c        reflection,extension and contraction coefficients.                 
c                                  
c	call intpr('NELMIN',-1,0,1)
                                         
c       validity checks on input.                                           
c                                                                           
      rcoeff=1.0d0
      ecoeff=2.0d0
      ccoeff=5.0d-1
      kcount=icount                                                         
      icount=0                                                              
      if(reqmin.le.0d0) icount=icount-1                                     
      if (n.gt.20) icount=icount-10                                         
      if (konvge.le.0) icount=icount-100                                    
      if (icount.lt.0) then
         konvge=1
         return  
      end if                                             
c                                                                           
      jcount=konvge                                                         
      dn=dfloat(n)                                                          
      nn=n+1                                                                
      dnn=dfloat(nn)                                                        
      del=1d0                                                               
c                                                                           
c        construction of initial simplex.                                   
c                                                                           
1001  do 1 i=1,n                                                            
1     p(i,nn)=start(i)                                                      
      call funct(n,start,z,
     &              iflags,vars,npts,nallpars,ndvars,iwhich,ndpars)

      y(nn)=z                                                               
      sum=z                                                                 
      summ=z*z                                                              
      do 2 j=1,n                                                            
      x=start(j)
      start(j)=start(j)+step(j)*del                                         
      do 3 i=1,n                                                            
3     p(i,j)=start(i)                                                       
      call funct(n,start,z,
     &              iflags,vars,npts,nallpars,ndvars,iwhich,ndpars)
      y(j)=z                                                                
      sum=sum+z                                                             
      summ=summ+z*z                                                         
2     start(j)=x
      icount=icount+nn
c                                                                           
c       simplex construction complete                                       
c                                                                           
c        find highest and lowest y values.ynewlo (=y(ihi) ) indicates       
c       the vertex of the simplex to be replaced.                           
c                                                                           
1000  ylo=y(1)                                                              
      ynewlo=ylo                                                            
      ilo=1                                                                 
      ihi=1                                                                 
      do 5 i=2,nn                                                           
      if(y(i).ge.ylo) goto 4                                                
      ylo=y(i)                                                              
      ilo=i                                                                 
4     if (y(i).le.ynewlo) goto 5                                            
      ynewlo=y(i)                                                           
      ihi=i                                                                 
5     continue                                                              
      sum=sum-ynewlo                                                        
      summ=summ-ynewlo*ynewlo                                               
c                                                                           
c      calculate pbar,the centroid of the simplex vertices                  
c          excepting that with y value ynewlo.                              
c                                                                           
      do 7 i=1,n                                                            
      z=0d0                                                                 
      do 6 j=1,nn                                                           
6     z=z+p(i,j)                                                            
      z=z-p(i,ihi)                                                          
7     pbar(i)=z/dn                                                          
c                                                                           
c      reflection through the centroid                                      
c                                                                           
      do 8 i=1,n                                                            
8     pstar(i)=(1d0+rcoeff)*pbar(i)-rcoeff*p(i,ihi)                         
      call funct(n,pstar,ystar,
     &              iflags,vars,npts,nallpars,ndvars,iwhich,ndpars)

      icount=icount+1                                                       
      if (ystar.ge.ylo) goto 12                                             
c                                                                           
c      successful reflection,so extension                                   
c                                                                           
      do 9 i=1,n                                                            
9     p2star(i)=ecoeff*pstar(i)+(1d0-ecoeff)*pbar(i)                        
      call funct(n,p2star,y2star,
     &              iflags,vars,npts,nallpars,ndvars,iwhich,ndpars)

      icount=icount+1                                                       
c                                                                           
c       retain extension or contraction.                                    
c                                                                           
      if (y2star.ge.ystar) goto 19                                            
10    do 11 i=1,n                                                           
11    p(i,ihi)=p2star(i)                                                    
      y(ihi)=y2star                                                         
      goto 900                                                              
c                                                                           
c     no extension                                                          
c                                                                           
12    l=0                                                                   
      do 13 i=1,nn                                                          
      if (y( i).gt.ystar) l=l+1                                             
13    continue                                                              
      if (l.gt.1) goto 19                                                   
      if (l.eq.0) goto 15                                                   
c                                                                           
c     contraction on the reflection side of the centroid.                   
c                                                                           
      do 14 i=1,n                                                           
14    p(i,ihi)=pstar(i)                                                     
      y(ihi)=ystar                                                          
c                                                                           
c      contraction on the  y(ihi) side of the centriod.                     
c                                                                           
15    do 16 i=1,n                                                           
16    p2star(i)=ccoeff*p(i,ihi)+(1d0-ccoeff)*pbar(i)                        
      call funct(n,p2star,y2star,
     &              iflags,vars,npts,nallpars,ndvars,iwhich,ndpars)

      icount=icount+1                                                       
      if (y2star.le.y(ihi)) goto 10                                         
c                                                                           
c       contract whole simplex.                                             
c                                                                           
      sum=0d0                                                               
      summ=0d0                                                              
      do 18 j=1,nn                                                          
      do 17 i=1,n                                                           
      p(i,j)=(p(i,j)+p(i,ilo))*5d-1                                         
17    min(i)=p(i,j)                                                         
      call funct(n,min,z,
     &              iflags,vars,npts,nallpars,ndvars,iwhich,ndpars)

      y(j)=z                                                                
      sum=sum+y(j)                                                          
18    summ=summ+y(j)*y(j)                                                   
      icount=icount+nn                                                      
      goto 901
c
c     retain reflection
c
19    do 20 i=1,n                                                           
20    p(i,ihi)=pstar(i)                                                     
      y(ihi)=ystar                                                          
900   sum=sum+y(ihi)                                                        
      summ=summ+y(ihi)*y(ihi)                                               
901   jcount=jcount-1                                                       
      if (jcount.ne.0) goto 1000                                            
c                                                                           
c     check to see if minimum reached.                                      
c                                                                           
      if (icount.gt.kcount) goto 22                                         
      jcount=konvge                                                         
      curmin=(summ-(sum*sum)/dnn)/dn                                        
c                                                                           
c     curmin is the variance of the fn values at the vertex.                
c                                                                           
      if (curmin.ge.reqmin) goto 1000                                       
c                                                                           
c       factorial tests to check that ynewlo is alocal minimum.             
c                                                                           
22    if (y(ihi).gt.y(ilo)) ihi=ilo                                         
      do 23 i=1,n                                                           
23    min(i)=p(i,ihi)                                                       
      ynewlo=y(ihi)                                                         
      if (icount.gt.kcount) then
         konvge=2
         return
      end if
      do 24 i=1,n
      del=step(i)*1.d-3
      min(i)=min(i)+del
      call funct(n,min,z,
     &              iflags,vars,npts,nallpars,ndvars,iwhich,ndpars)

      icount=icount+1
      if (z.lt.ynewlo) goto 25
      min(i)=min(i)-del-del
      call funct(n,min,z,
     &              iflags,vars,npts,nallpars,ndvars,iwhich,ndpars)

      icount=icount+1
      if (z.lt.ynewlo) goto 25
   24 min(i)=min(i)+del
      konvge=3
      return
c
c     restart procedure
c
   25 do 26 i=1,n
   26 start(i)=min(i)
      del=1.d-3
c
c     kcount=kcount+10000
c
c     icount=icount+10000
c
c     force termination after "icount" function evaluations
c
      goto 1001
       end                                                                  
	subroutine funct(ndim,start,z,
     &              iflags,vars,npts,nallpars,ndvars,iwhich,ndpars)
	implicit real*8 (a-h,o-z)
      dimension iflags(npts),iwhich(ndpars)
      dimension vars(npts,nallpars-ndvars)


	dimension start(ndim)

c        call intpr('FUNCT',-1,ndim,1)

c        call dblepr('EXtreme Vars',-1,vars(npts,nallpars-ndvars),1)


        call unlog(start,nallpars,ndpars,1.00d0,0.01d0)
c        call dblepr('FUNCT a,b,r',-1,start,3)
        call trblik(iflags,vars,npts,nallpars,ndvars,iwhich,ndpars
     &             ,start,dlogl)
	z=-dlogl
c        call dblepr('logl = ',-1,dlogl,1)


        call logem(start,nallpars,ndpars,1.00d0,0.01d0)
	return
	end


        subroutine logem(pars,nallpars,ndvars,aplus,bplus)
c
c transform the alphas and betas, and the rho
c
        implicit real*8 (a-h,o-z)
        dimension pars(nallpars)


        do i =1,ndvars
          pars(i)=dlog(pars(i)+aplus)
          pars(i+ndvars)=dlog(pars(i+ndvars)-bplus)
        end do
          pars(nallpars)=dlog(pars(nallpars)-.001)
        return

        end

        subroutine unlog(pars,nallpars,ndvars,aplus,bplus)  
c reverse transformation
        implicit real*8 (a-h,o-z)                                               
        dimension pars(nallpars)  
        do i =1,ndvars
          pars(i)=dexp(pars(i))-aplus
          pars(i+ndvars)=dexp(pars(i+ndvars))+bplus
        end do
          pars(nallpars)=dexp(pars(nallpars))+.001
        return       
        end

        function disfn(d2,alpha,beta)
        implicit real*8 (a-h,o-z)                                               
c                                                                               
c       input                                                                   
c
c          d2         : squared distance from putative source                   
c          alpha,beta : parameters to define profile of excess
c                       incidence around source
c
c       result
c
c          1+alpha*exp(-beta*d2)
c                                                                               
        disfn=1.0d0          
        arg=beta*d2
        if (arg.gt.20.0d0) return
        disfn=disfn+alpha*dexp(-arg)
        return
        end 
