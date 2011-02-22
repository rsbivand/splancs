# Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
# Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
# R port: copyright 1998-2000 by Roger S. Bivand

csr <- function(poly,npoints)
{
	ranpts(NULL,poly,npoints)[1:npoints,]
}

ranpts <- function(pts,poly,nprq)
{
# recursively keep generating points until we have more than enough...
#	cat("generating ",nprq," points\n")
	genpts <- gen(poly,nprq)
	npgen <- length(genpts)/2
#	cat("generated ",npgen," points\n")
	if(npgen != 0)pts <- rbind(pts,genpts)
	if (npgen < nprq)
	{
		pts <- ranpts(pts,poly,nprq-npgen)
	}
	pts
}

gen <- function(poly,npoints)
{
	areap <- areapl(poly)
	areab <- areapl(bboxx(bbox(poly)))

	xmin <- min(poly[,1],na.rm=TRUE)
	xmax <- max(poly[,1],na.rm=TRUE)
	ymin <- min(poly[,2],na.rm=TRUE)
	ymax <- max(poly[,2],na.rm=TRUE)

	aratio <- areab/areap
	ngen <- npoints*aratio
	
	xc <- xmin+runif(ngen)*(xmax-xmin)
	yc <- ymin+runif(ngen)*(ymax-ymin)
	
	xy <- cbind(xc,yc)
	pip(xy,poly)
	
}
