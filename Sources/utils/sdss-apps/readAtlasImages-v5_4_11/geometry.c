#include <stdio.h>
#include <stdlib.h>
#include "dervish.h"
#include "phSpanUtil.h"
#include "phGeometry.h"

/*****************************************************************************/
/*
 * Create/destroy a general polygon
 */
POLYGON *
phPolygonNew(int n)
{
   POLYGON *pol = shMalloc(sizeof(POLYGON));

   pol->n = n;
   pol->p = shMalloc(pol->n*sizeof(POINT));
   pol->is_convex = 0;

   return(pol);
}

void
phPolygonDel(POLYGON *p)
{
   if(p == NULL) return;

   shFree(p->p);
   shFree(p);
}

/*****************************************************************************/
/*
 * Create a set of polygons
 */
POLYGONS *
phPolygonsNew(int size)			/* initial number of polygons */
{
   int i;
   POLYGONS *poly = shMalloc(sizeof(POLYGONS));

   poly->size = size;
   poly->p = shMalloc(poly->size*sizeof(POLYGON *));
   poly->n = 0;

   for(i = 0; i < poly->size; i++) {
      poly->p[i] = NULL;
   }

   return(poly);
}

void
phPolygonsDel(POLYGONS *poly)
{
   int i;
   
   if(poly == NULL) return;

   for(i = 0; i < poly->n; i++) {
      phPolygonDel(poly->p[i]);
   }
   
   shFree(poly->p);
   shFree(poly);
}

/*****************************************************************************/
   
POLYGON *
phHullNewFromPoints(int n,		/* number of points */
		    const POINT **P)	/* the points */
{
   POLYGON *hull = shMalloc(sizeof(POLYGON));
   int i;
   
   hull->p = shMalloc(n*sizeof(POINT));
   hull->n = n;
   hull->is_convex = 1;

   for(i = 0; i < n; i++) {
      hull->p[i] = *P[i];
   }

   return(hull);
}

/*****************************************************************************/
/*
 * Find the two-dimensional convex hull
 *
 * the results should be "robust", and not return a wildly wrong hull,
 *	despite using floating point
 * works in O(n log n); I think a bit faster than Graham scan;
 * 	somewhat like Procedure 8.2 in Edelsbrunner's
 *  "Algorithms in Combinatorial Geometry".
 *
 * N.b. converted to not use floating point by RHL
 *
 * Ken Clarkson wrote this.  Copyright (c) 1996 by AT&T..
 * Permission to use, copy, modify, and distribute this software for any
 * purpose without fee is hereby granted, provided that this entire notice
 * is included in all copies of any software which is or includes a copy
 * or modification of this software and in all copies of the supporting
 * documentation for such software.
 * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED
 * WARRANTY.  IN PARTICULAR, NEITHER THE AUTHORS NOR AT&T MAKE ANY
 * REPRESENTATION OR WARRANTY OF ANY KIND CONCERNING THE MERCHANTABILITY
 * OF THIS SOFTWARE OR ITS FITNESS FOR ANY PARTICULAR PURPOSE.
 */
/*****************************************************************************/

static int
ccw(const POINT **P,
    int i,
    int j,
    int k)
{
   int a = P[i]->x - P[j]->x;
   int b = P[i]->y - P[j]->y;
   int c = P[k]->x - P[j]->x;
   int d = P[k]->y - P[j]->y;
   
   return((a*d - b*c <= 0) ? 1 : 0);	/* true if points i, j, k are
					   counterclockwise */
}

#define CMPM(Z,A,B) \
	v = (*A)->Z - (*B)->Z;\
	if(v > 0) return 1;\
	if(v < 0) return -1;

static int
cmpl(const POINT **a,
     const POINT **b)
{
   int v; 

   CMPM(x, a, b);
   CMPM(y, b, a);

   return 0;
}

static int
cmph(const POINT **a,
     const POINT **b)
{
   return cmpl(b,a);
}

static int
make_chain(const POINT **V,		/* the points in question */
	   int n,			/* number of points */
	   int (*cmp)(const POINT **, const POINT **))
{
   int i, j, s = 1;
   const POINT *tmp;			/* used in swapping V[] */
   
   qsort(V, n, sizeof(POINT *), (int (*)(const void *, const void *))cmp);
   
   for(i = 2; i < n; i++) {
      for(j = s; j >= 1 && ccw(V, i, j, j-1); j--) continue;
      
      s = j+1;
      tmp = V[s]; V[s] = V[i]; V[i] = tmp;
   }

   return s;
}

/*****************************************************************************/
/*
 * Find the convex hull of a set of points. The points are _not_ copied
 * into the POLYGON, and may be freed when the routine's returned
 */
POLYGON *
phHullFind(const POINT *points,		/* input points */
	   int n)			/* number of points */
{
   POLYGON *hull;				/* desired convex hull */
   int i;
   const POINT **P;			/* pointers to points[] */
   int u;

   if(n == 0) {
      return(NULL);
   }
/*
 * Setup pointers to points[]
 */
   P = shMalloc((n + 1)*sizeof(POINT));	/* the extra position is used */
   for(i = 0; i < n; i++) {
      P[i] = &points[i];
   }
/*
 * Find points in convex hull
 */
   u = make_chain(P, n, cmpl);		/* lower hull */
   P[n] = P[0];
   u += make_chain(&P[u], n-u+1, cmph);	/* upper hull */
/*
 * Create desired struct, and clean up
 */
   hull = phHullNewFromPoints(u, P);

   shFree(P);

   return(hull);
}

/*****************************************************************************/
/*
 * Return a set of convex polygons that approximate an OBJMASK
 */
POLYGONS *
phPolygonsFromObjmask(const OBJMASK *om)
{
   int i, j;
   int n;				/* number of endpoints in om */
   POINT *points;			/* endpoints in om */
   POLYGONS *polygons;			/* desired polygons */
   SPAN *s;

   if(om == NULL) {
      return(NULL);
   }

   n = 2*om->nspan;
   points = shMalloc(n*sizeof(POINT));

   s = om->s;
   for(i = j = 0; j < om->nspan; j++, s++) {
      points[i].x = s->x1;
      points[i++].y = s->y;
      
      points[i].x = s->x2;
      points[i++].y = s->y;
   }

   polygons = phPolygonsNew(1);
   polygons->p[0] = phHullFind(points, n);
   polygons->n = 1;

   shFree(points);

   return(polygons);
}

/*****************************************************************************/

void
phPolygonPrint(FILE *fd,		/* output descriptor, or NULL */
	       const char *hdr,		/* header string, or NULL */
	       const POLYGON *p)	/* polygon to print */
{
   int i;

   if(fd == NULL) {
      fd = stdout;
   }

   if(p == NULL) {
      fprintf(fd, "# No vertices\n");
      return;
   }
   
   fprintf(fd, "# %d vertices: area %g", p->n, phPolygonArea(p));
   if(hdr != NULL) {
      fprintf(fd, "  hdr \"%s\"", hdr);
   }
   fprintf(fd, "\n");
   
   for(i = 0; i < p->n; i++) {
      fprintf(fd, "%d %d\n", p->p[i].x, p->p[i].y);
   }
}

void
phPolygonsPrint(FILE *fd,		/* output descriptor, or NULL */
	       const char *hdr,		/* header string, or NULL */
	       const POLYGONS *poly)	/* polygons to print */
{
   int i;

   if(fd == NULL) {
      fd = stdout;
   }

   if(poly == NULL) {
      fprintf(fd, "# No polygons\n");
      return;
   }
   
   fprintf(fd, "# %d polygons: area %g", poly->n, phPolygonsArea(poly));
   if(hdr != NULL) {
      fprintf(fd, "  hdr \"%s\"", hdr);
   }
   fprintf(fd, "\n");
   
   for(i = 0; i < poly->n; i++) {
      phPolygonPrint(fd, NULL, poly->p[i]);
      if(i < poly->n - 1) {
	 fprintf(fd, "\n");
      }
   }
}

/*****************************************************************************/
/*
 * Return the area of a polygon or set of polygons
 */
double
phPolygonArea(const POLYGON *p)
{
   double area;				/* desired area */
   int i;

   if(p == NULL || p->n < 2) {
      return(0.0);
   }
   
   area = p->p[p->n - 1].x*p->p[0].y - p->p[p->n - 1].y*p->p[0].x;
   for(i = 1; i < p->n; i++) {
      area += p->p[i - 1].x*p->p[i].y - p->p[i - 1].y*p->p[i].x;
   }

   return(0.5*area);
}

double
phPolygonsArea(const POLYGONS *poly)
{
   double area = 0.0;			/* desired area */
   int i;

   if(poly == NULL) {
      return(area);
   }
   
   for(i = 0; i < poly->n; i++) {
      area += phPolygonArea(poly->p[i]);
   }

   return(area);
}

/*****************************************************************************/
/*
 * Create a stand-alone test executable
 */
#if 0
#include <stdio.h>
#include <assert.h>
#include <time.h>

#define SEC(T) (double)(T)/(double)CLOCKS_PER_SEC

#define N 10000
POINT points[N];

static int
read_points(void)
{
   int n = 0;
   char buf[100];

   while (fgets(buf, sizeof(buf), stdin)) {
      if(buf[0] == '#') {
	 continue;
      }
      
      if(sscanf(buf, "%d %d", &points[n].x, &points[n].y) != 2) {
	 fprintf(stderr,"Read error: %s", buf);
	 exit(1);
      }
      
      assert(++n <= N);
   }

   return n;
}

int
main(int argc,
     char **argv)
{
   POLYGON *hull;
   int n;
   clock_t t0,t1;
   
   n = read_points();
   t0 = clock();

   hull = phHullFind(points, n);
   
   if(hull != NULL) {
      phPolygonPrint(stdout, NULL, hull);
   }

   t1 = clock();
   printf("Convex hull computed in %f sec\n",SEC(t1-t0) );
   
   exit(0);
}
#endif
