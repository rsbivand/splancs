# Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
# Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
# R port: copyright 1998-2000 by Roger S. Bivand

kernel2d <- function(pts,poly,h0,nx=20,ny=20,kernel='quartic')
{
	if(!is.points(pts))stop('Invalid points argument')
	
	if(!is.points(poly))stop('Invalid poly argument')
	
	nptsk <- npts(pts)
	npoly <- length(poly[,1])
	poly <- rbind(poly,c(poly[1,1],poly[1,2]))
        storage.mode(poly) <- "double"
        storage.mode(pts) <- "double"
	
	
	xrang <- range(poly[,1],na.rm=TRUE)
	yrang <- range(poly[,2],na.rm=TRUE)
	
	bb <- bboxx(bbox(poly))
	a1 <- xrang[1]
	a2 <- xrang[2]
	b1 <- yrang[1]
	b2 <- yrang[2]
	cat("Xrange is ",a1,a2,"\n")
	cat("Yrange is ",b1,b2,"\n")
	xgrid <- rep(0,nx)
	ygrid <- rep(0,ny)
	zgrid <- matrix(0,nx,ny)
        storage.mode(xgrid) <- "double"
        storage.mode(ygrid) <- "double"
        storage.mode(zgrid) <- "double"
	if(kernel=='quartic')
	{
		cat('Doing quartic kernel\n')
#		library.dynam('splancs','krnqrt.o')
		storage.mode(zgrid) <- "double"
		klist <- .Fortran("krnqrt",
			pts[,1],
			pts[,2],
			as.integer(nptsk),
			poly[,1],
			poly[,2],
			as.integer(npoly),
			as.double(h0),
			as.double(a1),
			as.double(a2),
			as.double(b1),
			as.double(b2),
			as.integer(nx),
			as.integer(ny),
			xgrid=xgrid,
			ygrid=ygrid,
			zgrid=zgrid,
			PACKAGE="splancs")
		klist$zgrid[klist$zgrid<0] <- NA
		res <- list(x=klist$xgrid,y=klist$ygrid,z=klist$zgrid,
			 h0=h0,kernel=kernel)
	}
	else
	{
		stop('Invalid kernel function specification')
	}
res
}

	
	
