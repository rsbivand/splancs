C Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
C Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
C R port: copyright 1998-2000 by Roger S. Bivand

      subroutine sort(x,n)                                                  
c                                                                           
c     shellsort algorithm                                                   
c     n     : number of elements to be sorted                               
c     x     : on enter an array of dimension at least n containing          
c             real numbers                                                  
c             on output first n elements of x are sorted from smallest      
c             to largest                                                    
c                                                                           
      implicit real*8 (a-h,o-z)                                             
      dimension x(n)                                                        
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
      kk=kk-m                                                               
      goto 3                                                                
    4 continue                                                              
      goto 2                                                                
      end                                                                   

