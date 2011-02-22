c---------------------------------------------------------------
c
c this common block stores the area and convex/concave status of 
c polygon.
c it gets included by the "include 'bounds.cmn'" lines
c 
      common /bounds/area,iconvx

c area of polygon
      real*8 area

c if the polygon is convex, iconvx is set to 1, else its concave.
      integer iconvx

c---------------------------------------------------------------
