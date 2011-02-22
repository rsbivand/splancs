# Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
# Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
# R port: copyright 1998-2000 by Roger S. Bivand

stvmat <- function(pts,times,poly,tlim,s,tm)
{
#	mathlib.dynam('splancs','stsecal.o')
# full spacetime variance/covariance matrix
        storage.mode(pts) <- "double"
        storage.mode(times) <- "double"
        storage.mode(poly) <- "double"
        storage.mode(s) <- "double"
        storage.mode(tm) <- "double"
	tlow <- min(tlim)
        tupp <- max(tlim)

	ns <- length(s)
        nt <- length(tm)
        tm <- sort(tm)

	np <- npts(poly)
        xp <- c(poly[, 1],poly[1,1])
        yp <- c(poly[, 2],poly[1,2])

	n <- npts(pts)
	

	x <- pts[, 1]
	y <- pts[, 2]
        z <- times

        svec <- vector(mode="numeric", length=ns)
        tvec <- vector(mode="numeric", length=nt)

	smat <- matrix(0,ncol=ns,nrow=n)
        tmat <- matrix(0,ncol=nt,nrow=n)

        long <- ns*nt
        cmat <- matrix(0,ncol=long,nrow=long)

slist <- .Fortran("stsecal",
                x,
                y,
                as.integer(n),
                xp,
                yp,
                as.integer(np),
                s,
                as.integer(ns),
		as.double(smat),
		as.double(svec),
		as.double(z),
		as.double(tlow),
		as.double(tupp),
		as.double(tm),
		as.integer(nt),
		as.double(tmat),
		as.double(tvec),
		as.integer(long),
		varmat=as.double(cmat),
		PACKAGE="splancs")

 ans <- matrix(slist$varmat,nrow=long,byrow=TRUE)
 ans <- array(ans,dim=c(nt,ns,nt,ns))
ans
}

