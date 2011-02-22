# Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
# Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
# R port: copyright 1998-2000 by Roger S. Bivand
# !is.finite() bugfix by Konstantin Krivoruchko 2005

kernrat <- function(pts1,pts2,poly,h1,h2,nx=20,ny=20,kernel='quartic')
{
	kern1 <- kernel2d(pts1,poly,h1,nx,ny,kernel)
	kern2 <- kernel2d(pts2,poly,h2,nx,ny,kernel)
	z <- kern1$z/kern2$z
        is.na(z) <- !is.finite(z)
	
	list(x=kern1$x,y=kern1$y,z=z,h=c(h1,h2),kernel=kernel)
}


	
	
