# Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
# Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
# R port: copyright 1998-2000 by Roger S. Bivand

stkhat <- function(pts, times, poly, tlimits, s, tm)
{
#	mathlib.dynam('splancs','stkhat.o')
# space-time khat function
        storage.mode(poly) <- "double"
        storage.mode(pts) <- "double"
        storage.mode(s) <- "double"
        storage.mode(times) <- "double"
        storage.mode(tm) <- "double"
	tlow <- min(tlimits)
	thigh <- max(tlimits)
	s <- sort(s)
	ns <- length(s)
	tm <- sort(tm)
	nt <- length(tm)
	x1 <- pts[, 1]
	y1 <- pts[, 2]
	t1 <- times
	n1 <- length(x1)
	nc <- npts(poly)
	xc <- c(poly[, 1], poly[1, 1])
	yc <- c(poly[, 2], poly[1, 2])
	hs <- rep(0, ns)
	ht <- rep(0, nt)
	hst <- matrix(0, ns, nt)
	storage.mode(hs) <- "double"
	storage.mode(ht) <- "double"
	storage.mode(hst) <- "double"
	ktlist <- .Fortran("stkhat",
		x1,
		y1,
		t1,
		as.integer(n1),
		xc,
		yc,
		as.integer(nc),
		s,
		as.integer(ns),
		tm,
		as.integer(nt),
		as.double(tlow),
		as.double(thigh),
		hs = hs,
		ht = ht,
		hst = hst,
		PACKAGE="splancs")

	list(s=s,t=tm,ks = ktlist$hs, kt = ktlist$ht, kst = ktlist$hst)
}
