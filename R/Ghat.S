# Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
# Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
# R port: copyright 1998-2000 by Roger S. Bivand

Ghat <- function(pts,s)
{
	Gdists <- nndistG(pts)$dists
	plt(Gdists,s)
}

