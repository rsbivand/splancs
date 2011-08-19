# Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
# Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
# R port: copyright 1998-2001 by Roger S. Bivand

inpip <- function (pts, poly, bound=NULL, quiet=TRUE)
{
# coercion to points 110819
        if (!is.points(pts)) pts <- as.points(pts)
        if (!is.points(poly)) poly <- as.points(poly)
	seq(1:npts(pts))[inout(pts, poly, bound, quiet)]
}
