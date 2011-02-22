# Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
# Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
# R port: copyright 1998-2000 by Roger S. Bivand

pointmap <- function(pts,add=FALSE,axes=TRUE,xlab="",ylab="",asp=1,...)
{
	if(add) 
		{
		points(pts[,1:2,drop=FALSE],...)
		}
	else	
		{
		plot(pts[,1:2,drop=FALSE],
		axes=axes,xlab=xlab,ylab=ylab,asp=asp,...)
	}
}
