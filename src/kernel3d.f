C Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
C Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
C R port: copyright 1998-2000 by Roger S. Bivand


      subroutine kern3d(xpts,ypts,zpts,npts,xg,nx,yg
     &   ,ny,zg,nz,hxy,hz,rkarr)
      implicit double precision (a-h,o-z)
      dimension xpts(npts),ypts(npts),zpts(npts)
      dimension xg(nx),yg(ny),zg(nz)
      dimension rkarr(nx,ny,nz)

c now loop....
      do ix = 1,nx
        xgrid=xg(ix)
        do iy = 1,ny
          ygrid=yg(iy)
          do iz = 1,nz
            zgrid=zg(iz)
            rkarr(ix,iy,iz)=0
            do ip = 1,npts
              
              dstxys=(sqrt((xpts(ip)-xgrid)**2 +
     &             (ypts(ip)-ygrid)**2))/hxy
              dstzs=abs(zpts(ip)-zgrid)/hz

              if(dstxys.lt.1 .and. dstzs.lt.1) then

                xyk = (dstxys)**4-2*(dstxys)**2+1
                 zk = (dstzs)**4-2*(dstzs)**2+1
c                 call dblepr('xyk*zk',-1,xyk*zk,1)
                 rkarr(ix,iy,iz)=rkarr(ix,iy,iz)+
     &                (1/(hxy*hz))*(xyk*zk)
              end if
            end do
          end do
        end do
      end do


      end
