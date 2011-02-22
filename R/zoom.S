# Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
# Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
# R port: copyright 1998-2000 by Roger S. Bivand

zoom <- function(quiet=FALSE,out=FALSE,...)
{
	if(!out){
		if(!quiet){
			cat("Enter 2 points as corners of zoom area\n")
		}
		pointmap(bboxx(bbox(as.points(locator(2,type="l")))),type="n",...)
	}
	if(out){
		usr <- par()$usr
		xwidth <- usr[2]-usr[1]
		ywidth <- usr[4]-usr[3]
		xmin <- usr[1]-xwidth
		xmax <- usr[2]+xwidth
		ymin <- usr[3]-ywidth
		ymax <- usr[4]+ywidth
		pointmap(cbind(c(xmin,xmax),c(ymin,ymax)),type='n')
	}
	invisible(0)
}

