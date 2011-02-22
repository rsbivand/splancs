# Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
# Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
# R port: copyright 1998-2000 by Roger S. Bivand

sbox <- function(pts,xfrac=.1,yfrac=.1)
{
	if(!is.points(pts))stop('Points argument not valid point data')
	
	xr <- range(pts[,1],na.rm=TRUE)
	yr <- range(pts[,2],na.rm=TRUE)
	
	xw <- xr[2]-xr[1]
	xr[1] <- xr[1]-xfrac*xw
	xr[2] <- xr[2]+xfrac*xw

	yw <- yr[2]-yr[1]
	yr[1] <- yr[1]-yfrac*yw
	yr[2] <- yr[2]+yfrac*yw
	
	cbind(c(xr[1], xr[2], xr[2], xr[1]), c(yr[1], yr[1], yr[2], yr[2]))
}

