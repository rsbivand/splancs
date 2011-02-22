# Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
# Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
# R port: copyright 1998-2000 by Roger S. Bivand

#area <- function(poly)
areapl <- function(poly)
{
#	library.dynam("splancs","areapl.o")

	np <- length(poly[,1])
        storage.mode(poly) <- "double"
	polyx <- c(poly[,1],poly[1,1])
	polyy <- c(poly[,2],poly[1,2])
	
	answer <- 0

	alist <- .Fortran(	"areapl",
			polyx,
			polyy,
			as.integer(np),
			as.double(answer),
			PACKAGE="splancs")
	alist[[4]]	
	
}
