/* Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
 Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
 R port: copyright 1998-2001 by Roger S. Bivand

This code added June 2001 to cover Rainer Hurling's <rhurlin@gwdg.de>
concern about points on polygon boundaries - the polygon is rescaled to fit
within [-1,+1] and boundary hit is set to within +/- 0.000001 at that
scaling. Scaling code taken from VR/spatial krc.c,  by W. N. Venables 
and B. D. Ripley  Copyright (C) 1994-9. Original C code from Barry
Rowlingson.

*/

#define FALSE 0
#define TRUE  !FALSE
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

void
frset_pip(double xl, double xu, double yl, double yu);
static void
dscale_pip(double xo, double yo, double *xs, double *ys);
void ptinpoly (int *presult, double xpt, double ypt, double *xbuf,
	double *ybuf, int numpts);

static double pxl1, pxu1, pyl1, pyu1, pxl2, pyl2;


void ptinpoly1 (int *presult, double *xpt, double *ypt, double *xbuf,
	double *ybuf, int *numpts, double *bb, int *npts)
{
	int i;
	double *xs, *ys;
	double xpts, ypts;
	xs = (double *) R_alloc((size_t) numpts[0], sizeof(double));
	ys = (double *) R_alloc((size_t) numpts[0], sizeof(double));

	frset_pip(bb[0], bb[1], bb[2], bb[3]);
	    
	for (i=0; i<numpts[0]; i++)
		dscale_pip(xbuf[i], ybuf[i], &xs[i], &ys[i]);
	for (i=0; i<npts[0]; i++)
	{
		dscale_pip(xpt[i], ypt[i], &xpts, &ypts);
        	ptinpoly (&presult[i], xpts, ypts, xs, ys, numpts[0]);
	}
	return;

}


void
frset_pip(double xl, double xu, double yl, double yu)
{
    pxl1 = xl;
    pyl1 = yl;
    pxu1 = xu;
    pyu1 = yu;
    pxl2 = (pxu1 + pxl1) / 2;
    pyl2 = (pyu1 + pyl1) / 2;
	
}

static void
dscale_pip(double xo, double yo, double *xs, double *ys)
{

/* Scales (xo, yo) to (xs, ys) ( -1 < xs, ys < 1 ) within bounding box */
    *xs = (xo - pxl2) / (pxu1 - pxl2);
    *ys = (yo - pyl2) / (pyu1 - pyl2);
}

/* Function to determine whether a point lies inside, outside or on the
 * boundary of a polygon.
 * Count the number of times a horizontal line from the point to minus
 * infinity crosses segments of the polygon. Special cases arise when
 * the polygon points lie on this horizontal line and the direction
 * (in Y) of the lines on either side of the offending points have to
 * be considered. Further complications arise if whole polygon segments
 * lie along the horizontal line. In this case the directions of the
 * polygon segments at either end of the one or more horizontal
 * segments have to be considered. Assume the 1st and last polygon 
 * points are the same.
 */
void ptinpoly (int *presult, double xpt, double ypt, double *xbuf,
	double *ybuf, int numpts)
{
/*    int   i;*/
    int   numcrosses;
    int   ptr;
    double ratio;
    double xcross;
    double ydif;
    int   thisyup,lastyup;
    double fmax2(double x, double y), fmin2(double x, double y);
        
    /* First decide on the direction of the segments leading from the
       last points to the first point in the polygon (taking any horizontal
       segments into account. */
    ptr = numpts - 2;
    while ((ybuf[0] == ybuf[ptr]) && (ptr != 0)) ptr--;
    lastyup = FALSE;	/* down */
    if (ybuf[0] > ybuf[ptr]) lastyup = TRUE;	/* up */
    
    /* Loop through each segment of the polygon - only stopping 
       prematurely if the point is found to lie on one of the polygon
       segments. */
    numcrosses = 0;
    *presult = 1;
    ptr = 0;
    while ((ptr != numpts-1) && (*presult != 0))
    {
    	/* Does this polygon segment go up or down? */
    	if (ybuf[ptr] < ybuf[ptr+1]) thisyup = TRUE;
    	if (ybuf[ptr] > ybuf[ptr+1]) thisyup = FALSE;
    	
        /* Does the point lie within the Y bounds of the segment? */
        if ((ypt < fmax2 (ybuf[ptr],ybuf[ptr+1]))
         && (ypt > fmin2 (ybuf[ptr],ybuf[ptr+1])))
        {
	    /* Could the horz line from the point possibly cross
	       the segment? */
	    if (xpt >= fmin2 (xbuf[ptr],xbuf[ptr+1]))
	    {
	    	/* Does the horz line from the point definitely cross
	    	   the segment? */
	    	if (xpt <= fmax2 (xbuf[ptr],xbuf[ptr+1]))
	    	{
	    	    /* Work out whether the horz line from the point does
	    	       in fact cross the polygon segment. If the segment
	    	       is horizontal then the point must lie on the
	    	       polygon segment. If it is not horizontal then the
	    	       crossing point of the two lines must be worked
	    	       out and examined to see if it is to the right or
	    	       the left of the data point. Rounding errors are
	    	       significant and have to be dealt with. */
	    	    ydif = ybuf[ptr+1] - ybuf[ptr];
	    	    if (ydif != 0.0)
	    	    {
	    	    	ratio = (ypt - ybuf[ptr]) / ydif;
	    	    	xcross = xbuf[ptr] + 
	    	    		 (ratio * (xbuf[ptr+1] - xbuf[ptr]));
	    	        if (xcross < xpt) numcrosses++;
	    	        if ((xcross-xpt <  0.000001)
	    	         && (xcross-xpt > -0.000001)) *presult = 0;
	    	    }
	    	    else
	    	    {
	    	    	numcrosses++;
	    	    	*presult = 0;
	    	    }
	    	}
	    	else numcrosses++;
	    }
        }
        else
        {
            /* The Y value of the data point does not lie inside those of
               the current polygon boundary segment. Does the current 
               polygon boundary point have the same Y-value as the data
               point? */
            if (ypt == ybuf[ptr])
            {
            	/* Yes. Is the data point the same as the current polygon
            	   boundary point? */
            	if (xpt == xbuf[ptr]) *presult = 0;
            	else
            	{
		    /* If the segment is horizontal then work out whether
		       the data point lies on the line? If not then we
		       just ignore this segment. */
		    if (ybuf[ptr] == ybuf[ptr+1])
            	    {
            	    	if ((xpt >= fmin2 (xbuf[ptr],xbuf[ptr+1]))
            	    	 && (xpt <= fmax2 (xbuf[ptr],xbuf[ptr+1])))
            	    	    *presult = 0;
            	    }
            	    else
            	    {
            	    	/* Could the horz line possibly cross the segment? */
            	    	if (xpt > xbuf[ptr])
            	    	{
            	    	    /* The segment is not horz so simply check
            	    	       whether this polygon boundary point is an
            	    	       extrema or not. */
            	    	    if (thisyup == lastyup) numcrosses++;
            	    	}
            	    }
            	}
            }
        }
	lastyup = thisyup;    
    	ptr++;
    }
    
    /* Now we can decide whether the point is inside, outside, or on the
       edge of the polygon. */
    if (*presult != 0) 
    {
    	if (numcrosses % 2 == 0) *presult = 1;
    			    else *presult = -1;
    }
}  /* end ptinpoly */


