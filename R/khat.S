# Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
# Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
# R port: copyright 1998-2002 by Roger S. Bivand

khat <- function(pts, poly, s, newstyle=FALSE, checkpoly=TRUE)
{
        storage.mode(poly) <- "double"
        storage.mode(pts) <- "double"
        storage.mode(s) <- "double"
	ptsx <- pts[,1]
	ptsy <- pts[,2]
	npt <- npts(pts)

	ns <- length(s)
	s <- sort(s)

	np <- length(poly[,1])
	polyx <- c(poly[,1],poly[1,1])
	polyy <- c(poly[,2],poly[1,2])
	if (checkpoly) {
	    xar <- areapl(poly)
	    yar <- areapl(bboxx(bbox(poly)))
	    zar <- areapl(poly[chull(poly),])
	    if (xar < (yar/5) || xar < (zar/5)) 
	        stop("Polygon argument may be malformed")
	}

	hkhat <- rep(0,times=ns)
	counts <- matrix(0, nrow=npt, ncol=ns)
	mode(counts) <- "integer"
	hkhats <- matrix(0, nrow=npt, ncol=ns)
        storage.mode(hkhat) <- "double"
        storage.mode(hkhats) <- "double"

	klist <- .Fortran("trykh",
		ptsx,
		ptsy,
		as.integer(npt),
		polyx,
		polyy,
		as.integer(np),
		s,
		as.integer(ns),
		hkhat,
		as.integer(counts),
		hkhats,
		as.integer(npt*ns),
		PACKAGE="splancs")
	if (newstyle) {
		res <- list(khat=as.numeric(klist[[9]]),
			counts=matrix(klist[[10]], nrow=npt, ncol=ns),
			khats=matrix(klist[[11]], nrow=npt, ncol=ns),
			s=s)
			class(res) <- "khat"
	} else {
		res <- as.numeric(klist[[9]])
	}
	res
}

print.khat <- function(x, ...) {
	print(x$khat, ...)
	invisible(x)
}

plot.khat <- function(x, ...) {
	plot(x=x$s, y=x$khat, type="l", ylab="khat", xlab="distance", ...)
	invisible(x)
}


