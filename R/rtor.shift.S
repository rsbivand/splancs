# Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
# Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
# R port: copyright 1998-2000 by Roger S. Bivand

rtor.shift <- function(pts,rect)
{
	if(missing(rect))rect <- bboxx(bbox(pts))
	xsc <- max(rect[,1])-min(rect[,1])
	ysc <- max(rect[,2])-min(rect[,2])
	xsh <- runif(1)*xsc
	ysh <- runif(1)*ysc
	tor.shift(pts,xsh,ysh,rect)
}
