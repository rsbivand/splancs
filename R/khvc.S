# Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
# Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
# R port: copyright 1998-2000 by Roger S. Bivand

khvc <- function(pts1,pts2,poly,s)
{
#	mathlib.dynam('splancs','khvc.o')
        storage.mode(poly) <- "double"
        storage.mode(pts1) <- "double"
        storage.mode(pts2) <- "double"
        storage.mode(s) <- "double"
	anslist <- list(0,0,0,0)
        names(anslist) <- c('varmat','k11','k12','k22')
	ns <- length(s)

	np <- npts(poly)
        xp <- c(poly[, 1],poly[1,1])
        yp <- c(poly[, 2],poly[1,2])

	n1 <- npts(pts1);n2 <- npts(pts2);n <- n1+n2
	
	x <- c(pts1[, 1], pts2[, 1])
	y <- c(pts1[, 2], pts2[, 2])

        bvec <- vector(mode="numeric", length=ns)
	table <- matrix(0,ncol=ns,nrow=n)
        cmat <- matrix(0,ncol=ns,nrow=ns)
        v11 <- vector(mode='numeric',length=ns)
        v12 <- vector(mode='numeric',length=ns)
	v22 <- vector(mode='numeric',length=ns)

slist <- .Fortran("khvc",
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
		vark1=as.double(v11),
		vark12=as.double(v12),
		vark2=as.double(v22),
		covmat=as.double(cmat),
		PACKAGE="splancs")

 ans <- matrix(slist$covmat,nrow=ns,byrow=TRUE)
anslist$varmat <-  ans
anslist$k11 <- slist$vark1
anslist$k12 <- slist$vark12
anslist$k22 <- slist$vark2
anslist
}

