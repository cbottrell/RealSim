#if !defined(PHSPANUTIL_H)
#define PHSPANUTIL_H

#include "dervish.h"
#include "phConsts.h"
/*
 * typedefs and prototypes to deal with spans and masks made from
 * lists of spans.
 */

typedef struct {
   short y, x1, x2;
} SPAN;					/* pragma SCHEMA */

typedef struct {
   const int id;			/* unique ID number for OBJMASK */
   int refcntr;				/* no. of times OBJMASK is referenced*/
   int size;			/* number of spans allocated */
   int nspan;			/* actual number of spans */
   int row0, col0;		/* offset from parent region */
   int rmin, rmax;		/* bounding box */
   int cmin, cmax;
   int npix;				/* number of pixels in object */
   SPAN *s;				/* pointer to spans */
   PIX *data;				/* data of pixels in SPANs */
   float sum;				/* sum of all pixel values in OBJMASK*/
   void *user;				/* scratch pointer; available to user*/
} OBJMASK;				/* pragma SCHEMA */

   /* This enumerated type lists the values for each type in a REGION's
      SPANMASK. Note that only the first 8 will be displayed on an saoimage
      unless you fiddle with the display proc
    */

typedef enum {
   S_MASK_INTERP = 0,		/* pixel's value has been interpolated */
   S_MASK_SATUR,		/* pixel is/was saturated  */
   S_MASK_NOTCHECKED,		/* pixel was NOT examined for an object */
   S_MASK_OBJECT,		/* pixel is part of some object */
   S_MASK_BRIGHTOBJECT,		/* pixel is part of bright object */
   S_MASK_BINOBJECT,		/* pixel is part of binned object */
   S_MASK_CATOBJECT,		/* pixel is part of a catalogued object */
   S_MASK_SUBTRACTED,		/* model has been subtracted from pixel */
   S_MASK_GHOST,		/* pixel is part of a ghost */
   S_MASK_CR,			/* pixel is part of a cosmic ray */
   S_NMASK_TYPES		/* number of types; MUST BE LAST */
} S_MASKTYPE;

#define NMASK_TYPES 10			/* must equal S_NMASK_TYPES */

/* Let's make sure we have moved interfaces */
#define SPAN_COOKIE 0xc00cee

typedef struct spanmask {
/*
 * These fields are private.
 *
 * Note that the cookie must appear at offset < sizeof(MASK) if it's not to
 * be fooled by malloc reusing memory
 */
    unsigned int cookie;
    const struct spanmask *parent;
    int nchild;
/*
 * and these are public
 */
    CHAIN *masks[NMASK_TYPES];
    int nrow, ncol;			/* size of mask */
} SPANMASK;				/* pragma SCHEMA */

/*****************************************************************************/

void initSpanObjmask(void);		/* MUST be called first */

OBJMASK *
phObjmaskNew(
	    int nspans			/* number of spans */
	    );
OBJMASK *
phObjmaskRealloc(OBJMASK *sv,		/* Old OBJMASK */
		 int nspans);		/* new number of spans */
void
phObjmaskDel(OBJMASK *sv);

OBJMASK *
phObjmaskCopy(const OBJMASK *sv, float dr, float dc);

CHAIN *
phObjmaskChainCopy(const CHAIN *chain, float dr, float dc);

void
phObjmaskChainDel(
		  CHAIN *mask
		  );

OBJMASK *
phObjmaskFromRect(int x1, int y1, int x2, int y2);

OBJMASK *
phObjmaskFromCircle(int yc, int xc, int r);

SPANMASK *
phSpanmaskNew(int nrow, int ncol);

SPANMASK *
phSubSpanmaskNew(const SPANMASK *sm);

SPANMASK *
phSpanmaskCopy(const SPANMASK *ism, float dr, float dc);

void
phSpanmaskDel(
	      SPANMASK *mask
	      );
void
phObjmaskAddToSpanmask(const OBJMASK *om,
		       SPANMASK *mask,
		       S_MASKTYPE which);
CHAIN *
phSpanmaskFromRect(
		   int x1, 
		   int y1, 
		   int x2, 
		   int y2
		   );

SPANMASK *
phSpanmaskResetCorner(
		     SPANMASK *mask,
		     int row0,
		     int col0
		     );
OBJMASK *
phObjmaskResetCorner(
		     OBJMASK *mask,
		     int row0,
		     int col0
		     );

int phObjmaskIsCanonical(const OBJMASK *sv);

void
phCanonizeObjmask(
		  OBJMASK *sv,		/* OBJMASK to work on */
		  int nearly_sorted	/* is sv almost sorted already? */
		  );
void
phCanonizeObjmaskChain(CHAIN *chain,	/* CHAIN of OBJMASKs to work on */
		       int canonize_objmasks, /* canonize each OBJMASK too? */
		       int nearly_sorted); /* are OBJMASKS almost sorted? */

OBJMASK *phObjmaskTranspose(const OBJMASK *om);
void
phObjmaskBBSet(OBJMASK *om);

int
phObjmaskIntersect(const OBJMASK *sv1, 
		   const OBJMASK *sv2,
		   float dy,		/* offset to be applied to sv2 */
		   float dx);
OBJMASK *
phObjmaskMerge(OBJMASK *sv1,		/* target OBJMASK */
	       const OBJMASK *sv2,	/* input OBJMASK */
	       float dy,		/* Shift input by this amount before
					   merging */
	       float dx);
OBJMASK *
phObjmaskMergeWithCircle(OBJMASK *om,	/* mask to set (or create if NULL) */
			 float rowc, float colc, /* centre of circle */
			 float r);	/* radius */
void
phObjmaskClearWithCircle(OBJMASK *om,	/* mask to clear part of */
			 float rowc, float colc, /* centre of circle */
			 float r);	/* radius */

OBJMASK *
phMergeObjmaskChain(
		    const CHAIN *chain /* array of OBJMASK */
		    );
CHAIN *
phObjmaskChainMerge(
		   CHAIN *c1,		/* array of OBJMASK */
		   CHAIN *c2		/* array of OBJMASK */
		   );
int
phRectIntersectMask(
		    const CHAIN *chain, /* array of OBJMASK */
		    int x1,		/* boundaries of rectangle */
		    int y1,
		    int x2,
		    int y2
		    );
int
phObjmaskIntersectMask(const CHAIN *chain, /* array of OBJMASK */
		       const OBJMASK *sv); /* the OBJMASK in question */
int
phPixIntersectObjmask(
		      const OBJMASK *sv, /* input OBJMASK */
		      int x,		/* pixel position */
		      int y
		      );

int
phPixIntersectMask(
		   const CHAIN *chain,	/* array of OBJMASK */
		   int x,		/* pixel position */
		   int y
		   );
CHAIN *
phFindSVChainBB(
		const CHAIN *chain,	/* array of OBJMASK */
		int x1,			/* boundaries of rectangle */
		int y1,
		int x2,
		int y2
		);
CHAIN *
phFindSVChainBBCircle(
		      const CHAIN *chain, /* array of OBJMASK */
		      int x,		/* boundaries of rectangle */
		      int y,
		      int rad
		      );
CHAIN *
phTrimMaskToRect (
		  const CHAIN *iMask,	/* array of OBJMASK */
		  int x1,		/* boundaries of rectangle */
		  int y1,
		  int x2,
		  int y2
		  );
CHAIN *
phMaskIntersection(
		const CHAIN *mask1, 
		const CHAIN *mask2
);
CHAIN *phMaskNotIntersection(const CHAIN *mask1, const CHAIN *mask2);
OBJMASK *phObjmaskIntersection(const OBJMASK *obj, const CHAIN *mask2);
CHAIN *phObjmaskIntersectionChain(const OBJMASK *om1, const CHAIN *mask2);
OBJMASK *phObjmaskUnion(const OBJMASK *om1, const OBJMASK *om2);
OBJMASK *phObjmaskOrObjmask(OBJMASK *om1, const OBJMASK *om2);
OBJMASK *phObjmaskAndObjmask(OBJMASK *sv1, const OBJMASK *sv2);
OBJMASK *phObjmaskAndNotObjmask(OBJMASK *sv1, const OBJMASK *sv2);
OBJMASK *phObjmaskNotIntersectionObjmaskChain(const OBJMASK *obj,
					      const CHAIN *mask2);
OBJMASK *phObjmaskNotIntersectionObjmask(const OBJMASK *om1,
					  const OBJMASK *om2);
void
phMaskSetFromObjmask(const OBJMASK *om, /* an OBJMASK */
		     MASK *mask,	/* mask to set */
		     const int val);	/* value to OR into mask */
void
phMaskSetFromObjmaskChain(
			 const CHAIN *chain, /* array of OBJMASK */
			 MASK *mask,	/* mask of REGION objs were found in */
			 const int val	/* value to OR into mask */
			 );

OBJMASK **
phObjmaskSplitSimple(const OBJMASK *om,	/* initial OBJMASK */
		     int *nmask,	/* number of resulting OBJMASKs */
		     int alloc_one);	/* allocate a copy of om to return
					   if *nmask == 1? */

void
phSpanmaskSetAsPix(
		   SPANMASK *sm,	/* SPANMASK to set */
		   int row,		/* pixel location */
		   int col,
		   const S_MASKTYPE val /* type of mask */
		   );
int phObjmaskNumPix(
		    const OBJMASK *mask
		    );

int
phObjmaskFindSpan(const OBJMASK *om,	/* the mask */
		  int rc, int cc,	/* desired pixel */
		  int i0);		/* initial guess for index, or -1 */

OBJMASK *
phObjmaskGrow(const OBJMASK *om,	/* OBJMASK to grow */
	      const REGION *reg,	/* REGION that om lives in */
	      int n);			/* by how many pixels to grow */
void
phObjmaskChainGrow(CHAIN *chain,	/* the chain of OBJMASKs to grow */
		  const REGION *reg,	/* REGION that objects live in */
		  int n);		/* by how many pixels to grow */
OBJMASK *
phObjmaskGrowToSuperpixel(const OBJMASK *om, int n);

/*****************************************************************************/
/*
 * utilities to exchange pixels between REGIONs and OBJMASKs
 */
void
phObjmaskSetFromRegion(OBJMASK *om,	/* the OBJMASK to set */
		       const REGION *reg); /* the region wherein it dwells */
void
phRegionSetFromObjmask(REGION *reg,	/* the region to set */
		       const OBJMASK *om); /* and the OBJMASK */
void
phRegionSetValFromObjmask(REGION *reg,	/* the region to set */
			  const OBJMASK *om, /* and the OBJMASK */
			  int val);	/* to this value */
void
phRegionSetValFromObjmaskChain(REGION *reg, /* the region to set */
			       const CHAIN *chain, /* the CHAIN of OBJMASKs
						      specifying pixels */
			       int val); /* to set to this value */
#endif
