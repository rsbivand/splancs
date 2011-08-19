# Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
# Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
# R port: copyright 1998-2000 by Roger S. Bivand

as.points <- function(...)
{
	nv <- nargs()
	fargs <- list(...)
	if(nv==2){
		l1 <- length(fargs[[1]])
		l2 <- length(fargs[[2]])
		if(l1==l2){
			pts <- cbind(fargs[[1]],fargs[[2]])
		}
		else{
			stop("Cannot make points from different length vectors")
		}
		
	}
	else{
		if(nv==1){
			if(is.list(fargs[[1]])){
				fargs <- fargs[[1]]
				if(any(names(fargs)=="x") & any(names(fargs)=="y")){
					arx <- fargs$x
					ary <- fargs$y
					if(length(arx)!=length(ary)){
						stop("Cannot make points from different length x and y list components!")
					}
					else{
						pts <- cbind(arx,ary)
					}
				} else if (ncol(fargs) == 2) {
# 110819
                                        pts <- as.matrix(fargs)
                                        warning("Incorrectly structured list coerces to points")
				} else {
					stop("Cannot make points from list without x and y components.")
				}
			}
			else{
				if(is.points(fargs[[1]]))pts <- fargs[[1]]
				else stop("Cannot make points from this object")
			}
		
		}
		else{
			stop("Cannot make object into points!")
		}
	}

	pts
}
