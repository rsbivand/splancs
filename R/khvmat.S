# Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
# Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
# R port: copyright 1998-2000 by Roger S. Bivand

khvmat <- function(pts1,pts2,poly,s)
{
#	mathlib.dynam("splancs","khvmat.o")
        storage.mode(poly) <- "double"
        storage.mode(pts1) <- "double"
        storage.mode(pts2) <- "double"
        storage.mode(s) <- "double"
	ns <- length(s)

	np <- npts(poly)
        xp <- c(poly[, 1],poly[1,1])
        yp <- c(poly[, 2],poly[1,2])

	n1 <- npts(pts1);n2 <- npts(pts2);n <- n1+n2
	

	smax <- max(s)
	x <- c(pts1[, 1], pts2[, 1])
	y <- c(pts1[, 2], pts2[, 2])
        bvec <- vector(mode="numeric", length=ns)
	table <- matrix(0,ncol=ns,nrow=n)
        cmat <- matrix(0,ncol=ns,nrow=ns)

slist <- .Fortran("khvmat",
                x,
                y,
                as.integer(n),
                as.integer(n1),
                as.integer(n2),
                xp,
                yp,
                as.integer(np),
                s,
                as.integer(ns),
		as.double(table),
		as.double(bvec),
		varmat=as.double(cmat),
		PACKAGE="splancs")

ans <- matrix(slist$varmat,nrow=ns,byrow=TRUE)
ans <- ans+t(ans)
diag(ans) <- diag(ans)/2

ans

}
