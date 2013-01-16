# Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
# Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
# R port: copyright 1998-2000 by Roger S. Bivand

#.First.lib <- function(lib,pkg)
#.onLoad <- function(lib,pkg){
.onAttach <- function(lib,pkg){
#	cat("\nSpatial Point Pattern Analysis Code in S-Plus\n")
#	cat("\n Version 2 - Spatial and Space-Time analysis\n")
        Smess <- "\nSpatial Point Pattern Analysis Code in S-Plus\n"
        Smess <- paste(Smess,
            "\n Version 2 - Spatial and Space-Time analysis\n")
        packageStartupMessage(Smess)
#	library.dynam("splancs", pkg, lib)
	invisible(0)
}

