# Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
# Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
# R port: copyright 1998-2001 by Roger S. Bivand

inout <- function (pts, poly, bound = NULL, quiet=TRUE) 
{
# coercion to points 110819
    if (!is.points(pts)) pts <- as.points(pts)
    if (!is.points(poly)) poly <- as.points(poly)
    storage.mode(poly) <- "double"
    storage.mode(pts) <- "double"
    xp <- c(poly[, 1], poly[1, 1])
    yp <- c(poly[, 2], poly[1, 2])
    np <- length(xp)
    nptsi <- npts(pts)
    ind <- logical(length = nptsi)
    if (!is.null(bound)) {
#
# sets in-polygon criterion for points equal to polygon
# boundaries (suggestion by Rainer Hurling <rhurlin@gwdg.de>)
#
	   if (!is.logical(bound)) 
		   stop("bound must be NULL, TRUE, or FALSE")
	   result <- integer(length=nptsi)
	   bb <- as.vector(apply(sbox(as.points(xp, yp)), 2, range))
	   za <- .C("ptinpoly1", as.integer(result), pts[, 1],
		pts[, 2], xp, yp,
		as.integer(np), as.double(bb), as.integer(nptsi),
		PACKAGE="splancs")

	   z <- ind
	   if (!quiet) {
	       bpts <- which(za[[1]] == 0)
	       if (length(bpts > 0)) {
	           cat("Points on boundary:\n")
	           print(bpts)
	       }
	       else cat("No points on boundary\n")
	   }
	   if(bound) z[which(za[[1]] <= 0)] <- TRUE
	   else z[which(za[[1]] < 0)] <- TRUE
    } else {
# ind converted to integer RSB 140930
        piplist <- .Fortran("inpip", pts[, 1], pts[,2], 
            as.integer(nptsi), xp, yp, 
            as.integer(np), as.integer(ind),
		PACKAGE="splancs")
        z <- ifelse(piplist[[7]] == 0, FALSE, TRUE)
    }
    z
}

