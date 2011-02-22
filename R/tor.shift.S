# Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
# Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
# R port: copyright 1998-2003 by Roger S. Bivand

tor.shift <- function(pts,xsh=0.0,ysh=0.0,rect)
{
	if(missing(rect))rect <- bboxx(bbox(pts))
	
	xoff <- min(rect[,1])
	yoff <- min(rect[,2])
	
	xsc <- (max(rect[,1])-xoff)
	ysc <- (max(rect[,2])-yoff)
	
	pts[,1] <- pts[,1]-xoff
	pts[,2] <- pts[,2]-yoff
	pts <- Shift(pts,xsh,ysh)
	pts[,1] <- (pts[,1] %% xsc )+xoff
	pts[,2] <- (pts[,2] %% ysc )+yoff
	pts
}
