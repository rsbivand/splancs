# Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
# Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
# R port: copyright 1998-2000 by Roger S. Bivand

mse2d <- function(pts,poly,nsmse,range)
{
#	library.dynam("splancs","twodimmse.o")
        storage.mode(poly) <- "double"
        storage.mode(pts) <- "double"
	ptsx <- pts[,1]
	ptsy <- pts[,2]
	nptsm <- npts(pts)
	
	a1 <- min(poly[,1],na.rm=TRUE)
	a2 <- max(poly[,1],na.rm=TRUE)
	b1 <- min(poly[,2],na.rm=TRUE)
	b2 <- max(poly[,2],na.rm=TRUE)
	
	hsmse <- range/nsmse
	amse <- rep(0,nsmse)
	t <- rep(0,nsmse)
        storage.mode(amse) <- "double"
        storage.mode(t) <- "double"
	
	mselist <- .Fortran("mse2d",
		ptsx,
		ptsy,
		as.integer(nptsm),
		a1,
		a2,
		b1,
		b2,
		as.integer(nsmse),
		as.double(hsmse),
		amse,
		t,
		PACKAGE="splancs")
		
	list(mse=mselist[[10]],h=mselist[[11]])
}
