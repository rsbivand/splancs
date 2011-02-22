# Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
# Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
# R port: copyright 1998-2000 by Roger S. Bivand

secal <- function(pts1, pts2, poly, s)
{
	sqrt(diag(khvmat(pts1, pts2, poly, s)))
}
