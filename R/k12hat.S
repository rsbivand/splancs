# Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
# Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
# R port: copyright 1998-2000 by Roger S. Bivand

k12hat <- function(pts1,pts2,poly, s)
{
#	mathlib.dynam("splancs","k12hat.o")
        storage.mode(poly) <- "double"
        storage.mode(pts1) <- "double"
        storage.mode(pts2) <- "double"
	ptsx1 <- pts1[,1]
	ptsy1 <- pts1[,2]
	npt1 <- npts(pts1)
	ptsx2 <- pts2[,1]
	ptsy2 <- pts2[,2]
	npt2 <- npts(pts2)

	ns <- length(s)
	s <- sort(s)

	np <- length(poly[,1])
	polyx <- c(poly[,1],poly[1,1])
	polyy <- c(poly[,2],poly[1,2])

	h12 <- rep(0,times=ns)
	h21 <- h12

	klist <- .Fortran("k12hat",
		ptsx1,
		ptsy1,
		as.integer(npt1),
		ptsx2,
		ptsy2,
		as.integer(npt2),
		polyx,
		polyy,
		as.integer(np),
		as.double(s),
		as.integer(ns),
		as.double(h12),
		as.double(h21),
		PACKAGE="splancs")
	klist[[12]]
}
