# Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
# Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
# R port: copyright 1998-2000 by Roger S. Bivand

mpoint <- function(...,cpch,add=FALSE,type="p")
{
	arglist <- list(...)

	nptsets <- length(arglist)
	if(missing(cpch)){
		alphanum <- c(as.character(1:9),LETTERS)
		cpch <- alphanum[1:nptsets]
		}
	cpch <- rep(cpch,length=nptsets)
	ipch <- 1
	xmax <- NA;ymax <- NA;xmin <- NA;ymin <- NA
	if(add==FALSE){
		
		for(ia in 1:nptsets){
			if(is.points(arglist[[ia]]))
			{
				ptsc <- arglist[[ia]]
				ymax <- max(ymax,ptsc[,2],na.rm=TRUE)
				ymin <- min(ymin,ptsc[,2],na.rm=TRUE)
				xmax <- max(xmax,ptsc[,1],na.rm=TRUE)
				xmin <- min(xmin,ptsc[,1],na.rm=TRUE)
			}
		}
		pointmap(cbind(c(xmin,xmax,xmin,xmax),c(ymin,ymax,ymax,ymin)),type="n")
	}	
	for(ia in 1:nptsets){
		ptsc <- arglist[[ia]]
		pointmap(ptsc,add=TRUE,pch=cpch[ipch],type=type)
		ipch <- ipch+1
	}
	invisible(0)
}	
