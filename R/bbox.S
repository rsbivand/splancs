# Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
# Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
# R port: copyright 1998-2000 by Roger S. Bivand

#bbox <- function(pts)
#{
#	if(!is.points(pts))stop('Points argument not valid point data')
#	
#	xr <- range(pts[,1],na.rm=TRUE)
#	yr <- range(pts[,2],na.rm=TRUE)
#	
#	cbind(c(xr[1],xr[2],xr[2],xr[1]),c(yr[1],yr[1],yr[2],yr[2]))
#}

bboxx <- function(obj) {
	cbind(c(obj[1,1],obj[1,2],obj[1,2],obj[1,1]),
	    c(obj[2,1],obj[2,1],obj[2,2],obj[2,2]))
}
