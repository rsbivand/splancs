# Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
# Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
# R port: copyright 1998-2000 by Roger S. Bivand

Shift <- function(pts,xsh=0.0,ysh=0.0)
{
	pts[,1] <- pts[,1]+xsh
	pts[,2] <- pts[,2]+ysh
	pts
}
