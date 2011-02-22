# Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
# Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
# R port: copyright 1998-2000 by Roger S. Bivand

n2dist <- function(pts1,pts2)
{
#	library.dynam("splancs","nndist.o")
        storage.mode(pts1) <- "double"
        storage.mode(pts2) <- "double"
	npts1 <- npts(pts1)
	npts2 <- npts(pts2)

# reserve arrays for the result
	dists <- rep(0,npts2)
	neighs <- rep(0,npts2)
        storage.mode(dists) <- "double"

	nnlist <- .Fortran("n2dist", pts1[,1],
				 pts1[,2],
	                         as.integer(npts1),
	                         pts2[,1],
	                         pts2[,2],
	                         as.integer(npts2),
	                         dists,
				as.integer(neighs),
				PACKAGE="splancs")
	list(dists=nnlist[[7]],neighs=nnlist[[8]])
}


