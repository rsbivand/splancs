# Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
# Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
# R port: copyright 1998-2000 by Roger S. Bivand

"kernel3d"<-
function(pts,times, xgr, ygr, zgr, hxy, hz)
{

#        library.dynam('splancs','kernel3d.o')
        pts3 <- cbind(pts,times)
	nx <- length(xgr)
	ny <- length(ygr)
	nz <- length(zgr)
	mat3 <- array(data = 0, dim = c(nx, ny, nz))
	storage.mode(mat3) <- "double"
        storage.mode(pts3) <- "double"
        storage.mode(xgr) <- "double"
        storage.mode(ygr) <- "double"
        storage.mode(zgr) <- "double"
	ans <- .Fortran("kern3d",
		pts3[, 1],
		pts3[, 2],
		pts3[, 3],
		as.integer(length(pts3[, 3])),
		xgr,
		as.integer(length(xgr)),
		ygr,
		as.integer(length(ygr)),
		zgr,
		as.integer(length(zgr)),
		as.double(hxy),
		as.double(hz),
		kernarr = mat3,
		PACKAGE="splancs")
	list(xgr = xgr, ygr = ygr, zgr = zgr, hxy = hxy, hz = hz, v = ans$
		kernarr)
}


# Local Variables:
# mode:S
# S-temp-buffer-p:t
# End:
