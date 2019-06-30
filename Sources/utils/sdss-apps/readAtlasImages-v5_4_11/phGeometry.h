#if !defined(PHGEOMETRY_H)
#define PHGEOMETRY_H 1

typedef struct {
   int x, y;				/* coordinates of point */
} POINT;

typedef struct {
   int n;				/* number of points */
   POINT *p;				/* the points in the polygon */
   int is_convex;			/* is this polygon convex? */
} POLYGON;

typedef struct {
   int n;				/* number of polygons */
   int size;				/* dimension of p[] */
   POLYGON **p;				/* array of polygons */
} POLYGONS;

/*****************************************************************************/

POLYGON *phPolygonNew(int n);
void phPolygonDel(POLYGON *p);
double phPolygonArea(const POLYGON *p);
void phPolygonPrint(FILE *fp, const char *hdr, const POLYGON *p);

POLYGONS *phPolygonsNew(int size);
void phPolygonsDel(POLYGONS *poly);
double phPolygonsArea(const POLYGONS *poly);
void phPolygonsPrint(FILE *fp, const char *hdr, const POLYGONS *poly);

POLYGON *phHullFind(const POINT *points, int n);
POLYGON *phHullNewFromPoints(int n, const POINT **P);

#if defined(PHSPANUTIL_H)
   POLYGONS *phPolygonsFromObjmask(const OBJMASK *om);
#endif

#endif
